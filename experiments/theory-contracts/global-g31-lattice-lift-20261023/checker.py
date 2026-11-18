"""Independent replay of a short nonsplitting certificate in the original cocycle."""
import ast
from collections import Counter, deque
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import pickle
import runpy
import numpy as np
import sympy as S

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
CHECKS=[]

def need(ok,name):
    if not bool(ok):raise RuntimeError(name)
    CHECKS.append(name)

def main():
    p=runpy.run_path(str(HERE/'probe.py'))
    refs,perms,pg,chosen,orbits,full,orbit_id,rays=pickle.loads((p['OLD']/'state.pkl').read_bytes())
    old=json.loads((p['OLD']/'results.json').read_text())
    raw=json.loads((HERE/'probe-results.json').read_text())
    cert=raw['prefix_results'][-1]['inconsistency_certificate']
    B=S.Matrix(old['integral_lattice_basis']);G=S.Matrix(old['integral_gram'])
    T=S.diag(*([S.Matrix([[1,-1],[1,1]])/2]*4));Bs=T*B
    need(Bs.T*Bs==G and G.det()==1,'same original even unimodular lattice')
    source=REPO/'experiments/theory-contracts/seam-joint-ward-state-20261014/operator_transport_check.py'
    tree=ast.parse(source.read_text());names={'eps_car','split_even_lattice','eps_e8'}
    body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names]
    need({n.name for n in body}==names,'extract exactly original cocycle functions')
    ns={'S':S,'q':S.ones(8,1),'s':S.ones(8,1)/2}
    exec(compile(ast.Module(body=body,type_ignores=[]),str(source),'exec'),ns)
    @lru_cache(None)
    def eps(x,y):return int(ns['eps_e8'](S.Matrix(x)/2,S.Matrix(y)/2))
    def doubled(v):return tuple(int(2*x) for x in v)
    base=[doubled(Bs[:,j]) for j in range(8)]
    def rmat(r):
        A=S.Matrix(p['real_twice'](r).tolist())/2
        return B.inv()*A*B
    mats=[rmat(r) for r in refs[:3]]
    need(all(all(x.q==1 for x in M) and M.T*G*M==G for M in mats),
         'three actual original integral lattice reflections')
    images=[[doubled(Bs*M[:,j]) for j in range(8)] for M in mats]
    qi=[int(G[j,j]/2)%2 for j in range(8)]
    for i in range(8):
        need(eps(base[i],base[i])==(-1)**qi[i],'original square refinement '+str(i))
        for j in range(8):
            need(eps(base[i],base[j])*eps(base[j],base[i])==(-1)**int(G[i,j]),
                 'original commutator '+str((i,j)))
    @lru_cache(None)
    def ordered(n,s=-1):
        image=base if s==-1 else images[s]
        point=(0,)*8;sgn=1
        for j,k in enumerate(n):
            b=image[j];unit=1
            if k<0:unit=eps(b,tuple(-v for v in b));b=tuple(-v for v in b)
            for _ in range(abs(k)):
                sgn*=unit*eps(point,b);point=tuple(a+c for a,c in zip(point,b))
        return sgn
    @lru_cache(None)
    def eta(s,n):return ordered(n,s)*ordered(n)
    def eval_word(word,j):
        n=tuple(int(i==j) for i in range(8));co=0;sgn=1
        for s in reversed(word):
            co^=sum((n[i]%2)<<(8*s+i) for i in range(8))
            sgn*=eta(s,n)
            n=tuple(int(x) for x in mats[s]*S.Matrix(n))
        return co,int(sgn==-1),n
    original_rows=[]
    for k,row in enumerate(cert):
        w,v,j=row['left_word'],row['right_word'],row['coordinate']
        A=S.eye(8);D=S.eye(8)
        for s in w:A=A*mats[s]
        for s in v:D=D*mats[s]
        need(A==D,'literal matrix relation '+str(k))
        co,a,n=eval_word(w,j);co2,b,n2=eval_word(v,j)
        need(n==n2 and co^co2==row['coefficients'],'literal coefficient row '+str(k))
        original_rows.append({**row,'rhs':a^b,'finite_coordinate_rhs':row['rhs']})
    co=rhs=0
    for row in original_rows:co^=row['coefficients'];rhs^=row['rhs']
    need(co==0 and rhs==1,'original cocycle certificate gives zero equals one')
    need(len(original_rows)==7,'seven original equations suffice')

    # Verify the whole finite quotient multiplication, not just certificate samples.
    C=np.array(raw['C'],dtype=np.int64)
    f=np.array([[p['bilinear'](C,x,y) for y in range(256)] for x in range(256)],dtype=np.int8)
    xx=np.arange(256);xor=xx[:,None]^xx[None,:]
    bits=np.array([[(x>>i)&1 for i in range(8)] for x in range(256)],dtype=np.int64)
    need(np.array_equal(f^f.T,(bits@np.array(G.tolist(),dtype=np.int64)@bits.T)%2),
         'full finite commutator equals original lattice pairing')
    need(Counter(int(f[x,x]) for x in range(256))=={0:136,1:120},'original plus quadratic space')
    phase_changes=[]
    for s,M in enumerate(mats):
        Mn=np.array(M.tolist(),dtype=np.int64)%2
        ph=np.array(p['phase_table'](Mn,C),dtype=np.int8)
        cols=[sum(int(Mn[i,j])<<i for i in range(8)) for j in range(8)]
        actions=np.array([p['act'](cols,x) for x in range(256)])
        need(len(set(map(int,actions)))==256,'finite action invertible '+str(s))
        need(np.array_equal(ph[xor]^ph[:,None]^ph[None,:],f^f[actions[:,None],actions[None,:]]),
             'every finite product preserved by canonical lift '+str(s))
        # Translate generator phases between the two cocycle coordinates.
        K=0
        for j in range(8):
            n=tuple(int(x) for x in M[:,j])
            t=int(ordered(n)==-1)^ (sum((n[i]//2)*qi[i] for i in range(8))%2)
            K|=t<<j
        phase_changes.append(K)
    correction=sum(k<<(8*s) for s,k in enumerate(phase_changes))
    need(all(r['rhs']==(r['finite_coordinate_rhs']^((r['coefficients']&correction).bit_count()%2))
             for r in original_rows),'explicit original-to-finite coordinate gauge on all certificate rows')

    # Recompute all subgroup decisions in the same 16-element original group.
    I=(2*np.eye(4,dtype=np.int64),np.zeros((4,4),dtype=np.int64))
    elems=[I];lookup={p['key'](I):0};words=[[]];queue=deque([0])
    while queue:
        j=queue.popleft()
        for si,s in enumerate(refs[:3]):
            y=p['mul'](s,elems[j]);yk=p['key'](y)
            if yk not in lookup:
                lookup[yk]=len(elems);elems.append(y);words.append([si]+words[j]);queue.append(len(elems)-1)
    need(len(elems)==16,'actual three-reflection group has order sixteen')
    table=[[lookup[p['key'](p['mul'](a,b))] for b in elems] for a in elems]
    def closure(gs):
        h={0};todo=[0]
        while todo:
            x=todo.pop()
            for s in gs:
                y=table[s][x]
                if y not in h:h.add(y);todo.append(y)
        return frozenset(h)
    sub={frozenset({0}):[]};todo=[frozenset({0})]
    while todo:
        h=todo.pop()
        for j in range(16):
            if j in h:continue
            gs=sub[h]+[j];k=closure(gs)
            if k not in sub:sub[k]=gs;todo.append(k)
    subs=[]
    for h,gs in sorted(sub.items(),key=lambda t:(len(t[0]),sorted(t[0]))):
        ans=p['solve']([elems[j] for j in gs],np.array(B.tolist(),dtype=np.int64),C)
        need(ans['split']==(len(h)<16),'all subgroup lift decisions '+str(sorted(h)))
        if ans['split']:
            lifted={0:(np.arange(256),np.zeros(256,dtype=np.int8))};todo=[0]
            gen_actions=[];fixed_rows=[]
            for s,j in enumerate(gs):
                M=rmat(elems[j]);Mn=np.array(M.tolist(),dtype=np.int64)%2
                cols=[sum(int(Mn[a,b])<<a for a in range(8)) for b in range(8)]
                ma=np.array([p['act'](cols,x) for x in range(256)])
                ph=np.array(p['phase_table'](Mn,C),dtype=np.int8)
                l=(ans['one_solution']>>(8*s))&255
                ph^=np.array([(l&x).bit_count()%2 for x in range(256)],dtype=np.int8)
                gen_actions.append((ma,ph))
                fixed_rows.extend((Mn.T-np.eye(8,dtype=np.int64)).tolist())
            while todo:
                j=todo.pop();ma,ph=lifted[j]
                for s,k in enumerate(gs):
                    mb,pb=gen_actions[s];dest=table[k][j]
                    entry=(mb[ma],ph^pb[ma])
                    if dest not in lifted:lifted[dest]=entry;todo.append(dest)
                    else:
                        need(all(np.array_equal(a,b) for a,b in zip(entry,lifted[dest])),
                             'independent complete section closure '+str((sorted(h),j,s)))
            need(set(lifted)==set(h),'positive section covers exactly subgroup '+str(sorted(h)))
            piv={}
            for row in fixed_rows:
                v=sum((int(x)%2)<<i for i,x in enumerate(row))
                while v:
                    k=v.bit_length()-1
                    if k in piv:v^=piv[k]
                    else:piv[k]=v;break
            ans['coboundary_dimension']=len(piv)
            ans['H1_dimension']=ans['dimension_of_splittings']-len(piv)
        subs.append({'elements':sorted(h),'generator_words':[words[j] for j in gs],**ans})
    need(len(subs)==23,'all twenty-three subgroups enumerated')
    out={'contract_id':'MATH.SEAM.G31.LATTICE.03','full_G31_order':46080,
         'character_kernel_order':256,'full_lift_group_order':46080*256,
         'global_character_extension_splits':False,'witness_subgroup_order':16,
         'original_basis':[[str(x) for x in Bs.row(i)] for i in range(8)],
         'original_reflection_matrices':[[[int(x) for x in M.row(i)] for i in range(8)] for M in mats],
         'original_cocycle_certificate':original_rows,'finite_to_original_generator_phase_changes':phase_changes,
         'certificate_xor':{'coefficient_mask':co,'rhs':rhs},
         'proper_subgroup_count':22,'all_proper_subgroup_restrictions_split':True,
         'subgroups':subs,'checks_passed':len(CHECKS),'checks':CHECKS,
         'sources':[{'path':str(source.relative_to(REPO)),'sha256':sha256(source.read_bytes()).hexdigest()}]}
    binary=set();kernel=[];Bi2=np.array((2*B.inv()).tolist(),dtype=np.int64)
    for g in full.values():
        v=Bi2@p['real_twice'](g)@np.array(B.tolist(),dtype=np.int64)
        if np.any(v%4):raise RuntimeError('nonintegral full group action')
        M=(v//4)%2;binary.add(tuple(map(int,M.flat)))
        if np.array_equal(M,np.eye(8,dtype=np.int64)):kernel.append(p['key'](g))
    need(len(binary)==23040 and len(kernel)==2,'full binary image and kernel plus/minus identity')
    out['binary_image_order']=len(binary);out['binary_kernel_order']=len(kernel)
    out['checks_passed']=len(CHECKS);out['checks']=CHECKS
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['global_character_extension_splits','full_lift_group_order',
          'witness_subgroup_order','certificate_xor','proper_subgroup_count','checks_passed']},indent=2))

if __name__=='__main__':main()
