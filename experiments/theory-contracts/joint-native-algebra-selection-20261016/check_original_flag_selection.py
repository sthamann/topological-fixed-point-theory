"""Exact original charged seam/state stabilizers on actual native 60 flags.

No selected clock, desired projector, coupling, or optimization rule.
The existing J-compatible Cartan isometry is used literally, not fitted.
"""
from pathlib import Path
from collections import deque, Counter
import sys, json, itertools as it
import sympy as s
R=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0,str(R))
from tfpt_explorer import cartan_source as cs
from tfpt_explorer.neutral_source_response import _native_neutral_branch

def need(ok,msg):
    if not ok: raise ArithmeticError(msg)
def cl(M):return M.applyfunc(s.simplify)
def compose(p,q):return tuple(p[q[i]] for i in range(len(p)))
def orbit_partition(elements, points, action):
    remaining=set(points);out=[]
    while remaining:
        x=min(remaining);orb={action(p,x) for p in elements}
        need(orb <= remaining,'Complete group orbits, no truncated subgroup')
        remaining-=orb;out.append(sorted(orb))
    return sorted(out,key=lambda q:(len(q),q))

def main():
    ld=cs._lattice_data(); roots=ld['roots'];ri={r:i for i,r in enumerate(roots)}
    T=s.Matrix.hstack(*[s.Matrix(r) for r in cs._TARGET_SIMPLE])
    targets=[tuple(T*s.Matrix(z)) for z in ld['root_coordinates']]
    need(set(targets)==cs._standard_doubled_roots(),'Literal original isometry maps all roots')
    x=s.symbols('x0:4')
    polys=[sum(z**4 for z in x),x[0]**2*x[1]**2+x[2]**2*x[3]**2,
           x[0]**2*x[2]**2+x[1]**2*x[3]**2,x[0]**2*x[3]**2+x[1]**2*x[2]**2,s.prod(x)]
    raw=s.Matrix([[2,2,-1,-1,-1,-1],[0,0,-1,-1,1,1],
                  [0,0,-1,1,-1,1],[0,0,1,-1,-1,1],[1,-1,0,0,0,0]])
    A=[s.expand(sum(raw[i,j]*[s.Rational(1,24),s.Rational(1,4),s.Rational(1,4),s.Rational(1,4),1][i]*polys[i] for i in range(5))) for j in range(6)]
    native_rays=cs._rays(roots);genidx=(0,13,2,3,1);gens=[];logicalgens=[]
    for idx in genidx:
        ray=native_rays[idx];p=bytes(ri[cs._native_reflection(r,ray)] for r in roots)
        need(sorted(p)==list(range(240)),'Native root permutation')
        # Build physical real matrix from its images of the lattice basis.
        real=s.Matrix.hstack(*[s.Matrix(cs._native_reflection(r,ray)) for r in cs._SIMPLE])*ld['basis'].inv()
        g=s.Matrix(4,4,lambda i,j:real[2*i,2*j]+s.I*real[2*i+1,2*j])
        sub=dict(zip(x,g.T*s.Matrix(x)))
        acts=[s.expand(a.subs(sub,simultaneous=True)) for a in A]
        lp=tuple(A.index(a) for a in acts)
        need(sum(lp[i]!=i for i in range(6))==2,'Actual quartic action is a transposition')
        gens.append(p);logicalgens.append(lp)
    # Left-compose via bytes.translate (exact permutation, all 240 roots).
    translations=[bytes(list(p)+list(range(240,256))) for p in gens]
    identity=bytes(range(240));group={identity:tuple(range(6))};queue=deque([identity])
    while queue:
        p=queue.popleft();lp=group[p]
        for trans,lg in zip(translations,logicalgens):
            prod=p.translate(trans);lprod=compose(lg,lp)
            if prod not in group:group[prod]=lprod;queue.append(prod)
            else:need(group[prod]==lprod,'Faithful common root/quartic group action')
    need(len(group)==46080 and len(set(group.values()))==720,'Full original G31 and S6, no guessed quotient')
    # chi_v(alpha)=exp(2pi i v.alpha), v=(0^5,1/2^3).
    char=tuple(int(sum(t[5:]))%4 for t in targets)
    deck=tuple(c%2 for c in char)
    need(Counter(char)=={0:52,2:60,1:64,3:64},'Actual original Z4 root sectors')
    # The source state Omega is fixed by all canonical lattice lifts; lift
    # characters act trivially on the neutral oscillator fields and Omega.
    source=_native_neutral_branch(0)
    retained=tuple(int(t in source['retained']) for t in targets)
    ar=source['R_matrix']
    ardiag=tuple((s.Matrix(t).T*ar*s.Matrix(t))[0] for t in targets)
    # Complete source objects, never only their restriction to the code.
    signatures={'unmarked':tuple(0 for _ in roots),'halfway_deck':deck,
                'oriented_quarter_character':char,
                'quarter_plus_neutral_R':tuple(zip(char,ardiag)),
                'quarter_plus_retained_48':tuple(zip(char,retained))}
    for omitted in range(1,4):
        alternate=_native_neutral_branch(omitted)
        retained_alt=tuple(int(t in alternate['retained']) for t in targets)
        signatures['quarter_plus_retained_48_omitted_'+str(omitted)]=tuple(zip(char,retained_alt))
    fold=json.loads((R/'experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json').read_text())
    M=s.Matrix(fold['explicit_Weyl_matrix']);need(M.T*M==s.eye(8),'Existing fold M, not fitted map')
    fold_details={}
    for alpha in (s.Rational(1,4),s.Rational(3,4)):
        mu=s.Matrix([(j+alpha)/8 for j in range(8)])
        d=M*(mu-s.ones(8,1)/2);nu=M*(2*mu)
        phase=tuple(s.Mod((nu.T*s.Matrix(t))[0]/2,1) for t in targets)
        zeroexpect=tuple((d.T*s.Matrix(t))[0]/2 for t in targets)
        key=str(alpha)
        signatures['actual_fold_g2_'+key]=phase
        signatures['fold_state_diagnostic_d_'+key]=zeroexpect
        signatures['fold_d_plus_internal_character_'+key]=tuple(zip(zeroexpect,char))
        fold_details[key]={'d_alpha':[str(v) for v in d],'g2_root_phase_multiplicities':{str(k):v for k,v in Counter(phase).items()},
                           'same_as_internal_Z4_character':all(phase[i]==s.Rational(char[i],4) for i in range(240)),
                           'd_norm_squared':str((d.T*d)[0])}
    flags=[(q,tuple(S)) for q in range(6) for S in it.combinations([j for j in range(6) if j!=q],2)]
    def fact(p,f):return (p[f[0]],tuple(sorted(p[j] for j in f[1])))
    results={}
    for name,sig in signatures.items():
        stab=[lp for p,lp in group.items() if all(sig[p[i]]==sig[i] for i in range(240))]
        limage=sorted(set(stab));qo=orbit_partition(limage,list(range(6)),lambda p,x:p[x]);fo=orbit_partition(limage,flags,fact)
        pairsets=list(it.combinations(range(6),2))
        fixedpairs=[S for S in pairsets if all(tuple(sorted(p[j] for j in S))==S for p in limage)]
        q5stab=[p for p in limage if p[5]==5]
        q5fixedpairs=[S for S in it.combinations(range(5),2) if all(tuple(sorted(p[j] for j in S))==S for p in q5stab)]
        results[name]={'native_stabilizer_order':len(stab),'logical_image_order':len(limage),
                       'q_orbits':qo,'flag_orbit_sizes':[len(o) for o in fo],
                       'fixed_q':[o[0] for o in qo if len(o)==1],
                       'fixed_flags':[o[0] for o in fo if len(o)==1],
                       'fixed_unrooted_pairs':fixedpairs,'if_q5_supplied_fixed_pairs':q5fixedpairs,
                       'comparison_original_flag_orbit':next(o for o in fo if (5,(0,1)) in o),
                       'flag_orbits':fo}
    # An explicit same-source witness for each residual label pair swap.
    # Root action fixes the full character and retained48, and any normalized
    # real cocycle lift fixes Omega. No C5-only compression enters this test.
    sig=signatures['quarter_plus_retained_48'];witnesses=[]
    for aa,bb in ((0,1),(2,3),(4,5)):
        target=list(range(6));target[aa],target[bb]=target[bb],target[aa];target=tuple(target)
        matches=[p for p,lp in group.items() if lp==target and all(sig[p[i]]==sig[i] for i in range(240))]
        need(bool(matches),'Actual residual independent pair swap exists')
        p=min(matches)
        gt=s.Matrix.hstack(*[s.Matrix(targets[p[ri[b]]]) for b in cs._SIMPLE])*T.inv()
        need(gt.T*gt==s.eye(8),'Exact witness source isometry')
        need(gt*ar*gt.T==ar,'Witness preserves full original neutral tensor')
        witnesses.append({'logical_transposition':[aa,bb],'source_real_matrix':[[str(v) for v in row] for row in gt.tolist()],
                          'root_permutation':list(p),'preserves_240_root_character':True,'preserves_retained_48':True,
                          'fixes_lattice_vacuum':True})
    N=cl(s.diag(2,s.sqrt(12),s.sqrt(12),s.sqrt(12),s.sqrt(24))*raw/s.sqrt(48))
    def polar(q,S):
        vv=[N[:,j]-(s.sqrt(6)-1)*N[:,q]/5 for j in S]
        return cl(sum([v*v.T for v in vv],s.zeros(5)))
    P4,P5=polar(4,(0,1)),polar(5,(0,1))
    dist=s.simplify(s.trace((P4-P5)**2))
    need(dist>0,'Residual swap changes actual polar operator, not only redundant labels')
    out={'status':'PASS_EXACT_ORIGINAL_SEAM_FLAG_STABILIZERS','date':'2026-10-17',
         'root_count':240,'native_group_order':len(group),'logical_group_order':720,
         'generator_indices':list(genidx),'actual_logical_generators':logicalgens,
         'quarter_character_multiplicities':dict(Counter(char)),
         'original_omitted_family':list(source['family'][0]),
         'original_vacuum_invariant':True,'results':results,
         'residual_pair_swap_witnesses':witnesses,'actual_fold_state_diagnostic':fold_details,
         'residual_original_pair_polar_distance_squared':str(dist),
         'scope':'Given original unitary E8 lattice source, accepted J-compatible source/native isometry and actual quartic embedding. The charged characteristic is retained on all roots; not identified with native J. Retained48/R use the already stipulated omitted family and do not derive it from P1. No fitted clock, coupling, P2 or endpoint-state selection. The physical raw P1 source/state dictionary remains separate.'}
    Path(__file__).with_name('original_flag_selection_results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'groups':[len(group),720],
        'cases':{n:{'native':r['native_stabilizer_order'],'logical':r['logical_image_order'],
                    'q_orbits':r['q_orbits'],'flag_orbit_histogram':dict(Counter(r['flag_orbit_sizes'])),
                    'fixed_flag_count':len(r['fixed_flags'])} for n,r in results.items()},
        'polar_pair_distance_squared':str(dist)},ensure_ascii=False,indent=2))
if __name__=='__main__':main()
