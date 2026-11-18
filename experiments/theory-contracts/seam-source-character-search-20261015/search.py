#!/usr/bin/env python3
"""Exact, bounded internal-character search; no state/time optimization or GPU."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import shutil
import time

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
FOLD=HERE.parent/'v16-operator-pullback-20261014'
SOURCE=HERE.parent/'seam-joint-ward-state-20261014/PROOF.txt'
PARENT_PACKET=HERE.parent/'seam-constructive-search-20261015/atlas_packet.json'
D=tuple(F(x,16) for x in (-8,0,1,-8,0,-3,-4,-4))
R=tuple(F(x) for x in (0,0,1,0,0,1,0,0))
LAMBDA=tuple(F(x) for x in (-4,0,0,-4,0,-2,-2,-2))
Q=tuple(2*x for x in R)


def encoded(value): return json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=False)
def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def strings(v): return [str(x) for x in v]
def dot(a,b): return sum((x*y for x,y in zip(a,b)),F(0))
def add(a,b): return tuple(x+y for x,y in zip(a,b))
def scale(t,a): return tuple(t*x for x in a)
def identity(): return [[F(int(i==j)) for j in range(8)] for i in range(8)]
def matvec(a,v): return tuple(dot(row,v) for row in a)
def mm(a,b): return [[sum((a[i][k]*b[k][j] for k in range(8)),F(0)) for j in range(8)] for i in range(8)]


def rank(rows):
    a=[[F(x) for x in row] for row in rows]
    pivot=0
    for col in range(8):
        hit=next((j for j in range(pivot,len(a)) if a[j][col]),None)
        if hit is None: continue
        a[pivot],a[hit]=a[hit],a[pivot]
        divisor=a[pivot][col]
        a[pivot]=[x/divisor for x in a[pivot]]
        for j in range(len(a)):
            if j!=pivot and a[j][col]:
                t=a[j][col]
                a[j]=[x-t*y for x,y in zip(a[j],a[pivot])]
        pivot+=1
        if pivot==len(a): break
    return pivot


def roots():
    """All E8 roots in doubled integer coordinates, with the even-spinor class."""
    out=[]
    for i,j in itertools.combinations(range(8),2):
        for a,b in itertools.product((-2,2),repeat=2):
            v=[0]*8;v[i]=a;v[j]=b;out.append(tuple(v))
    out += [s for s in itertools.product((-1,1),repeat=8) if sum(x<0 for x in s)%2==0]
    assert len(out)==240 and len(set(out))==240 and all(dot(v,v)==8 for v in out)
    return sorted(out)


def phase(h,root2): return (dot(h,root2)/2)%1
def phases(h,rs): return tuple(phase(h,a) for a in rs)
def character_order(h,rs):
    from math import lcm
    return lcm(*(x.denominator for x in phases(h,rs)))


def dynkin_type(simple,cartan):
    n=len(simple)
    edges={i:[j for j in range(n) if j!=i and cartan[i][j]==-1] for i in range(n)}
    degrees=[len(edges[i]) for i in range(n)]
    if n==1: return 'A1'
    if max(degrees)<=2: return f'A{n}'
    centers=[i for i in range(n) if degrees[i]==3]
    if len(centers)!=1 or sum(degrees)!=2*(n-1): return 'unclassified'
    center=centers[0];arms=[]
    for start in edges[center]:
        length=1;old=center;node=start
        while len(edges[node])==2:
            nxt=next(j for j in edges[node] if j!=old)
            old,node=node,nxt;length+=1
        arms.append(length)
    arms.sort()
    if arms==[1,1,n-3]: return f'D{n}'
    return {(1,2,2):'E6',(1,2,3):'E7',(1,2,4):'E8'}.get(tuple(arms),'unclassified')


def components(neutral):
    """Actual nonorthogonality components, exact ranks and simple-root graphs."""
    remaining=set(neutral);answer=[]
    regular=(1,2,4,8,16,32,64,128)
    while remaining:
        first=min(remaining);remaining.remove(first);component={first};todo=[first]
        while todo:
            a=todo.pop()
            linked={b for b in remaining if dot(a,b)!=0}
            remaining-=linked;component|=linked;todo+=sorted(linked)
        positive={a for a in component if dot(a,regular)>0}
        assert len(positive)*2==len(component)
        simple=sorted(a for a in positive if not any(tuple(x-y for x,y in zip(a,b)) in positive for b in positive))
        cartan=[[int(dot(a,b)/4) for b in simple] for a in simple]
        dimension=rank(component)
        assert len(simple)==dimension
        assert all(cartan[i][i]==2 for i in range(dimension))
        assert all(cartan[i][j] in (-1,0) for i in range(dimension) for j in range(dimension) if i!=j)
        answer.append({'type':dynkin_type(simple,cartan),'root_count':len(component),'rank':dimension,
                       'simple_roots_doubled':simple,'cartan_matrix':cartan,'roots':sorted(component)})
    return sorted(answer,key=lambda c:(c['rank'],c['type']))


def projector_from_moment(neutral):
    return [[(sum(a[i]*a[j] for a in neutral)-32*int(i==j))/F(32) for j in range(8)] for i in range(8)]


def verify_projector(neutral,projector):
    p=[[F(x) for x in row] for row in projector]
    expected=projector_from_moment(neutral)
    if p!=expected or mm(p,p)!=p or rank(p)!=5 or sum(p[i][i] for i in range(8))!=5:
        raise ValueError('Projector is not the rank-five neutral-root moment projector')
    for c in components(neutral):
        if c['type'] not in ('D5','A3'): raise ValueError('Wrong neutral-root components')
        for a in c['roots']:
            wanted=tuple(F(x) for x in a) if c['type']=='D5' else (F(0),)*8
            if matvec(p,a)!=wanted: raise ValueError('Projector does not resolve actual root components')
    return True


def field_phases(h,matrix):
    # Columns M e_j are actual original folded source-field charges in this dictionary.
    return [dot(h,[matrix[i][j] for i in range(8)])%1 for j in range(8)]


def bilinear_source_components(h,matrix,rs):
    source_h=tuple(dot(h,[matrix[i][j] for i in range(8)]) for j in range(8))
    integer_roots=[a for a in rs if all(x%2==0 for x in a)]
    neutral=[a for a in integer_roots if phase(source_h,a)==0]
    cs=components(neutral)
    return {'root_count':len(neutral),'semisimple_rank':sum(c['rank'] for c in cs),
            'central_Cartan_rank':8-sum(c['rank'] for c in cs),
            'actual_connected_components':[{k:v for k,v in c.items() if k!='roots'} for c in cs],
            'phase_block_grouping':'D2+D2+A3 when the four phase multiplicities are two; D2=A1+A1 is reducible',
            'scope':'neutral original D8 bilinear root algebra; distinct from the completed E8 neutral algebra'}


def describe(h,rs,with_roots=False):
    neutral=[a for a in rs if phase(h,a)==0]
    cs=components(neutral)
    return {'order':character_order(h,rs),
            'root_phase_counts':dict(sorted(Counter(str(x) for x in phases(h,rs)).items())),
            'components':[{k:v for k,v in c.items() if k!='roots' or with_roots} for c in cs],
            'neutral':neutral}


def group_data(rs):
    rows=[];signatures=set()
    for m,n in itertools.product(range(8),repeat=2):
        h=add(scale(m,D),scale(F(n,8),R))
        signature=phases(h,rs)
        if signature in signatures: raise ValueError('Source character group enumeration collided')
        signatures.add(signature)
        rows.append({'m':m,'n':n,'order':character_order(h,rs)})
    assert len(signatures)==64
    # Generator multiplication: 8d = lambda1+r/2; n therefore gains four.
    for m,n in itertools.product(range(8),repeat=2):
        dm=(m+1)%8;dn=(n+4*int(m==7))%8
        assert phases(add(scale(m+1,D),scale(F(n,8),R)),rs)==phases(add(scale(dm,D),scale(F(dn,8),R)),rs)
        assert phases(add(scale(m,D),scale(F(n+1,8),R)),rs)==phases(add(scale(m,D),scale(F((n+1)%8,8),R)),rs)
    return {'order':64,'order_counts':{str(k):v for k,v in sorted(Counter(row['order'] for row in rows).items())},
            'elements':rows,'generator_carry_relation':'8d=lambda1+r/2; m wrap adds n+=4 mod8'}


def source_controls(matrix,rs):
    cases=[]
    t=scale(F(1,4),LAMBDA)
    pstar=matvec(matrix,(F(1,2),)*4+(F(-1,2),)*4)
    for alpha in (F(0),F(1,4),F(3,4)):
        mu=tuple((F(j)+alpha)/8 for j in range(8))
        delta=tuple(x-F(1,2) for x in mu)
        nu=scale(2,mu)
        da=matvec(matrix,delta)
        nua=matvec(matrix,nu)
        assert nua==add(scale(2,da),scale(2,R))
        if alpha==F(1,4): assert da==D
        # This new internal marker is constructed, not substituted into K or the original clock.
        constructed=add(scale(2,da),scale(-alpha/2,R))
        assert constructed==t
        clock=describe(nua,rs)
        energy=dot(pstar,pstar)/2+dot(da,pstar)
        assert energy==0
        cases.append({'alpha':str(alpha),'M_delta':strings(da),'M_nu':strings(nua),
                      'physical_g2_order':clock['order'],'physical_g2_root_phase_counts':clock['root_phase_counts'],
                      'physical_g2_components':clock['components'],'constructed_internal_t':strings(constructed),
                      'original_quarter_boundary_phase':str(alpha),'K_kept':'N+|p|^2/2+(M delta_alpha).p',
                      'known_ground_pair_energy_check':['0',str(energy)]})
    return {'clock_definitions':{'mu':'(j+alpha)/8','delta':'mu-q/2','nu':'2mu; clock covector distinct from delta'},
            'canonical_internal_construction':'t=2 M delta_alpha-(alpha/2)r=lambda1/4; extra internal normalization/compensation choice',
            'cases':cases,'known_ground_pair_target_charges':[strings((F(0),)*8),strings(pstar)],
            'ground_pair_status':'existing Ls result only rechecked at its two supplied charges; no new minimization, state selection or source repair'}


def blind_source_classification():
    """Pure source-frame calculation: no M parameter, target constants or files."""
    alpha=F(1,4)
    delta=tuple((F(j)+alpha)/8-F(1,2) for j in range(8))
    source_r=(F(1,2),)*8
    source_lambda=add(scale(8,delta),scale(F(-1,2),source_r))
    assert source_lambda==tuple(F(j-4) for j in range(8))
    rs=roots();cases=[];selected=[]
    for u,v in itertools.product(range(4),repeat=2):
        h=add(scale(F(u,4),source_lambda),scale(F(v,4),source_r))
        row=describe(h,rs);neutral=row.pop('neutral')
        row.update(u=u,v=v,h_source=strings(h))
        types=sorted((c['type'],c['rank'],c['root_count']) for c in row['components'])
        row['D5_A3_candidate']=types==[('A3',3,12),('D5',5,40)]
        if row['D5_A3_candidate']:
            projector=projector_from_moment(neutral)
            verify_projector(neutral,projector)
            row['projector_source']=[[str(x) for x in r] for r in projector]
            row['projector_rank']=rank(projector)
            selected.append(row)
        cases.append(row)
    return {'alpha':str(alpha),'delta_source':strings(delta),'r_source':strings(source_r),
            'lambda1_source':strings(source_lambda),'case_count':len(cases),
            'primitive_order_four_count':sum(c['order']==4 for c in cases),
            'D5_A3_candidate_count':len(selected),'candidates':selected,'all_cases':cases,
            'reads_M_or_target_projector':False,
            'scope':'M-independent classification inside the stipulated original-frame E8 extension, not a derivation of that extension'}


def calculate(fold_path=FOLD/'fold-results.json'):
    # Finish the source-only classification before reading even the M input.
    blind=blind_source_classification()
    raw=json.loads(Path(fold_path).read_text())
    matrix=[[F(x) for x in row] for row in raw['explicit_Weyl_matrix']]
    transpose=list(map(list,zip(*matrix)))
    assert mm(matrix,transpose)==identity()
    assert tuple(F(x) for x in raw['transported_source_U1_charge'])==Q
    assert matvec(matrix,(F(1),)*8)==Q
    assert add(scale(8,D),scale(F(-1,2),R))==LAMBDA
    rs=roots()
    assert set(tuple(int(2*x) for x in matvec(matrix,scale(F(1,2),a))) for a in rs)==set(rs)
    group=group_data(rs)
    cases=[];accepted=[];seen=set()
    for u,v in itertools.product(range(4),repeat=2):
        h=add(scale(F(u,4),LAMBDA),scale(F(v,4),R))
        signature=phases(h,rs)
        assert signature not in seen;seen.add(signature)
        row=describe(h,rs)
        neutral=row.pop('neutral')
        row.update(u=u,v=v,h=strings(h),m=2*u,n=(2*v-u)%8)
        types=sorted((c['type'],c['rank'],c['root_count']) for c in row['components'])
        row['D5_A3_candidate']=types==[('A3',3,12),('D5',5,40)]
        if row['D5_A3_candidate']:
            p=projector_from_moment(neutral)
            verify_projector(neutral,p)
            fp=field_phases(h,matrix)
            row.update(projector_rank=rank(p),projector_from_neutral_moment=[[str(x) for x in r] for r in p],
                       original_complex_source_field_phases=strings(fp),
                       original_complex_source_field_phase_counts=dict(sorted(Counter(str(x) for x in fp).items())),
                       original_real_source_field_phase_counts=dict(sorted(Counter(str(x) for x in fp+[-x%1 for x in fp]).items())),
                       original_source_bilinear_neutral=bilinear_source_components(h,matrix,rs),
                       scope='conditional internal marking; not elementary 5+3 source CAR and not a physical source selection')
            accepted.append(row)
        cases.append(row)
    assert len(seen)==16
    assert sum(row['order']<=4 for row in group['elements'])==16
    primitive=sum(row['order']==4 for row in cases)
    wrong=describe(scale(F(1,4),Q),rs);wrong.pop('neutral')
    controls=source_controls(matrix,rs)
    for row in accepted:
        h=tuple(F(x) for x in row['h'])
        row['known_ground_pair_marker_phases']=[str(dot(h,tuple(F(x) for x in p))%1) for p in controls['known_ground_pair_target_charges']]
    # Only now use M and the fixed coordinate P5, strictly for an audit.
    by_uv={(c['u'],c['v']):c for c in accepted}
    comparison=[]
    coordinate_p5=[[F(int(i==j and i<5)) for j in range(8)] for i in range(8)]
    for source_case in blind['candidates']:
        target_case=by_uv[(source_case['u'],source_case['v'])]
        source_p=[[F(x) for x in row] for row in source_case['projector_source']]
        transported=mm(mm(matrix,source_p),transpose)
        target_p=[[F(x) for x in row] for row in target_case['projector_from_neutral_moment']]
        assert transported==target_p
        assert matvec(matrix,tuple(F(x) for x in source_case['h_source']))==tuple(F(x) for x in target_case['h'])
        comparison.append({'u':source_case['u'],'v':source_case['v'],'m':target_case['m'],'n':target_case['n'],
                           'transported_source_projector_equals_target_moment':True,
                           'equals_standard_coordinate_P5':transported==coordinate_p5})
    assert len(blind['candidates'])==len(accepted)
    blind['post_classification_M_comparison']=comparison
    blind['fixed_coordinate_P5_used_only_after_classification']=True
    return {'verdict':'CONDITIONAL_INTERNAL_MARKING_CANDIDATES_PHYSICAL_SOURCE_OPEN',
            'domain':{'d':strings(D),'r':strings(R),'Q':strings(Q),'lambda1':strings(LAMBDA),'lambda2':strings(R),
                      'index_convention':'r=e3+e6 in one-based mathematical coordinates',
                      'lattice_plane':'span_R{d,r} intersect E8 = Z lambda1 + Z r; saturation argument remains in parent proof',
                      'candidate_formula':'h=(u lambda1+v r)/4, u,v in0..3; no input rank-five projector',
                      'matrix_source_key':'explicit_Weyl_matrix from fold-results.json'},
            'E8_roots':{'count':len(rs),'doubled_integer_coordinate_sha256':hashlib.sha256(encoded(rs).encode()).hexdigest()},
            'source_generated_group':group,'order_dividing_four_count':len(cases),'primitive_order_four_count':primitive,
            'D5_A3_candidate_count':len(accepted),'candidates':accepted,'all_order_dividing_four_cases':cases,
            'blind_original_source_frame_control':blind,
            'wrong_character_Q_over_4':wrong,'state_and_time_scope':controls,
            'proof_obligations':['physical origin of the source plane and quotient/localization',
                                 'physical identification of the constructed internal marker without erasing the given quarter character',
                                 'compatible actual source-field/intertwiner map, not dimension or root-count matching',
                                 'selection of the existing Ls ground-doublet state and an actual event law',
                                 'raw P1/P2 and T1-T8 remain open'],
            'unchanged_source_obstruction':'Existing unchanged-K Ls ground doublet and missing physical source selection are retained, not rediscovered or repaired by this internal character search',
            'physical_source_selected':False,'new_H_or_state_chosen':False}


def provenance(paths):
    return [{'path':str(Path(p).resolve()),'sha256':digest(p)} for p in paths]


def fresh_atlas_snapshot(path):
    """Reuse the existing deterministic MCP client, never its older saved packet."""
    client_path=HERE.parent/'seam-constructive-search-20261015/atlas_inputs.py'
    spec=importlib.util.spec_from_file_location('source_character_atlas_client',client_path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    client=module.Atlas()
    try:
        status=client.call('graph_status',{})
        queries=('problem:common-source-state-event-selection','problem:seam-partner-charge-state-time')
        states=[client.call('research_state',{'query':q,'lang':'en','max_chars':30000}) for q in queries]
        refs={}
        for query,state in zip(queries,states):
            if state.get('budget',{}).get('truncated') or not state.get('problems'):
                raise ValueError('Incomplete fresh Atlas problem packet: '+query)
            for problem in state['problems']:
                for finding in problem.get('findings',[]):
                    if finding.get('source_state')!='current':
                        raise ValueError('Stale original finding: '+finding['id'])
                    for ref in finding.get('sources',[]):
                        source=(REPO/ref['path']).resolve()
                        if not source.is_relative_to(REPO) or not source.is_file() or digest(source)!=ref['sha256']:
                            raise ValueError('Changed original Atlas finding source: '+str(source))
                        refs[(ref['path'],ref['sha256'])]=ref
        revisions={s['revision'] for s in states}
        if len(revisions)!=1: raise ValueError('Atlas revision changed during snapshot; rerun search')
        packet={'schema':1,'collection':'fresh MCP graph_status and complete exact research_state calls',
                'collected_at_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
                'atlas_revision':states[0]['revision'],'graph_status':status,'research_states':states,
                'sources':list(refs.values()),'source_verification_scope':'only original sources from reviewed findings',
                'proof_boundary':'navigation and typed prior evidence; no automatic physical-source conclusion'}
        Path(path).write_text(json.dumps(packet,indent=2,ensure_ascii=False)+'\n')
        return {'path':str(Path(path).resolve()),'sha256':digest(path),'revision':states[0]['revision'],
                'collection':packet['collection'],'original_finding_source_count':len(refs)}
    finally:
        client.close()


def verify_sources(result):
    for ref in result['source_manifest']:
        path=Path(ref['path'])
        if not path.is_file() or digest(path)!=ref['sha256']:
            raise ValueError('Source changed: '+str(path))
    snapshot=result.get('atlas_snapshot')
    if snapshot:
        path=Path(snapshot['path'])
        if digest(path)!=snapshot['sha256']: raise ValueError('Atlas snapshot changed')
        packet=json.loads(path.read_text())
        for ref in packet.get('sources',[]):
            path=(REPO/ref['path']).resolve()
            if not path.is_relative_to(REPO) or not path.is_file() or digest(path)!=ref['sha256']:
                raise ValueError('Atlas original source changed: '+str(path))
    return True


def search(output,fold_path=FOLD/'fold-results.json',atlas_packet=None):
    result=calculate(fold_path)
    inputs=[Path(__file__),Path(fold_path),FOLD/'FOLD-AUDIT.md',SOURCE]
    if atlas_packet is None:
        snap=Path(output).parent/'atlas_packet.snapshot.json'
        result['atlas_snapshot']=fresh_atlas_snapshot(snap)
        inputs.append(snap)
    elif atlas_packet is not False:
        atlas_packet=Path(atlas_packet)
        snap=Path(output).parent/'atlas_packet.snapshot.json'
        if atlas_packet.resolve()!=snap.resolve(): shutil.copyfile(atlas_packet,snap)
        result['atlas_snapshot']={'path':str(snap.resolve()),'sha256':digest(snap),
                                  'scope':'navigation/source packet only, not additional proof'}
        packet=json.loads(snap.read_text())
        if 'atlas_revision' in packet:
            result['atlas_snapshot'].update(revision=packet['atlas_revision'],
                                            collection=packet.get('collection','explicit input packet'),
                                            original_finding_source_count=len(packet.get('sources',[])))
        inputs.append(atlas_packet)
    result['source_manifest']=provenance(inputs)
    verify_sources(result)
    result['computed_payload_sha256']=hashlib.sha256(encoded({k:v for k,v in result.items() if k not in ('source_manifest','atlas_snapshot')}).encode()).hexdigest()
    Path(output).write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    return result


def verify(path):
    saved=json.loads(Path(path).read_text())
    verify_sources(saved)
    fold_path=next(ref['path'] for ref in saved['source_manifest'] if Path(ref['path']).name=='fold-results.json')
    actual=calculate(fold_path)
    for key,value in actual.items():
        if encoded(saved.get(key))!=encoded(value): raise ValueError('Result replay mismatch: '+key)
    checksum=hashlib.sha256(encoded(actual).encode()).hexdigest()
    if saved['computed_payload_sha256']!=checksum: raise ValueError('Payload digest mismatch')
    return {'verified':True,'candidates':actual['D5_A3_candidate_count'],'payload_sha256':checksum}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('search','verify'))
    parser.add_argument('--output',type=Path,default=HERE/'results.json')
    parser.add_argument('--fold-results',type=Path,default=FOLD/'fold-results.json')
    parser.add_argument('--atlas-packet',type=Path,default=None,help='Explicit existing input packet; default collects a fresh MCP snapshot')
    args=parser.parse_args()
    if args.command=='verify': print(json.dumps(verify(args.output)))
    else:
        result=search(args.output,args.fold_results,args.atlas_packet)
        print(json.dumps({k:result[k] for k in ('verdict','order_dividing_four_count','primitive_order_four_count','D5_A3_candidate_count')}))


if __name__=='__main__': main()
