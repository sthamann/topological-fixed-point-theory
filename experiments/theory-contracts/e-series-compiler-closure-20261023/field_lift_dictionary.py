"""E7-pair Chevalley operations versus original ordered E8 lattice lifts."""
import ast
from collections import Counter
from functools import lru_cache
from hashlib import sha256
import itertools
import json
from pathlib import Path
import numpy as np
import sympy as S

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
PINS={}

def functions(relative,names,ns):
    p=REPO/relative;raw=p.read_bytes();PINS[relative]=sha256(raw).hexdigest()
    nodes=[n for n in ast.parse(raw).body if isinstance(n,ast.FunctionDef) and n.name in names]
    if {n.name for n in nodes}!=set(names):raise RuntimeError('native definitions')
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(p),'exec'),ns)
    return ns

def need(ok,label):
    if not bool(ok):raise RuntimeError(label)

ns=functions('verification/v1_e8_glue.py',{'e8_roots','e8_simple_roots'}, {'np':np,'itertools':itertools})
roots=np.rint(2*np.array(ns['e8_roots']())).astype(np.int64)
B=S.Matrix(np.rint(2*ns['e8_simple_roots']()).astype(np.int64)).T
Bi=B.inv();coords=np.array((Bi*S.Matrix(roots.T)).tolist(),dtype=np.int64).T
lookup={tuple(v):i for i,v in enumerate(roots)}
basis=[lookup[tuple(int(x) for x in B[:,j])] for j in range(8)]
base=[tuple(int(x) for x in B[:,j]) for j in range(8)]
jns=functions('verification/v634_st31_structure.py',{'J_vec'},{'np':np})
J=np.column_stack([jns['J_vec'](v) for v in np.eye(8,dtype=np.int64)])
jp=np.array([lookup[tuple(J@v)] for v in roots])
dot=roots@roots.T
neg=np.array([lookup[tuple(-v)] for v in roots])
old=functions('experiments/theory-contracts/seam-joint-ward-state-20261014/operator_transport_check.py',
 {'eps_car','split_even_lattice','eps_e8'},{'S':S,'q':S.ones(8,1),'s':S.ones(8,1)/2})

@lru_cache(None)
def eps(u,v):
    c,h=u[0]%2,v[0]%2
    d=[(x-c)//2 for x in u];f=[(x-h)//2 for x in v]
    need(sum(d)%2==0 and sum(f)%2==0,'even lattice')
    z=c*(sum(f)//2)+c*h+sum(d[i]*f[j] for i in range(8) for j in range(i))
    z+=c*h*sum(i*(d[i]+f[i]) for i in range(8))
    return -1 if z%2 else 1

for i in range(8):
 for j in range(8):
    need(eps(base[i],base[j])==int(old['eps_e8'](S.Matrix(base[i])/2,S.Matrix(base[j])/2)), 'literal native cocycle')

E=np.array([[eps(tuple(u),tuple(v)) for v in roots] for u in roots],dtype=np.int8)
sums=np.full((240,240),-1,dtype=np.int64)
for i,j in zip(*np.where(dot==-4)):sums[i,j]=lookup[tuple(roots[i]+roots[j])]
nonzero=list(zip(*np.where(dot==-4)))

@lru_cache(None)
def ordered(n,images):
    point=(0,)*8;phase=1
    for b,k in zip(images,n):
        unit=eps(b,tuple(-x for x in b)) if k<0 else 1
        v=tuple(-x for x in b) if k<0 else b
        for _ in range(abs(k)):
            phase*=unit*eps(point,v);point=tuple(a+c for a,c in zip(point,v))
    return phase

base=tuple(base)
def lift_phases(perm):
    images=tuple(tuple(roots[perm[i]]) for i in basis)
    return np.array([ordered(tuple(map(int,n)),images)*ordered(tuple(map(int,n)),base) for n in coords],dtype=np.int8)

def tits(a):
    m=dot[a]//4;perm=np.array([lookup[tuple(v-int(k)*roots[a])] for v,k in zip(roots,m)])
    k=int(E[a,neg[a]]);phase=np.ones(240,dtype=np.int8)
    for i,w in enumerate(m):
        if w==1:phase[i]=-k*E[neg[a],i]
        elif w==-1:phase[i]=E[a,i]
        elif abs(w)==2:phase[i]=-k
    need(np.array_equal(phase*phase[perm],(-1)**(m%2)),'Tits square coroot character')
    return perm,phase

Jphase=lift_phases(jp)
rows=[];square_profiles=Counter();canonical_orders=Counter();correction_counts=Counter()
orbits=sorted({tuple(sorted((a,int(jp[a]),int(jp[jp[a]]),int(jp[jp[jp[a]]])))) for a in range(240)})
for oi,orb in enumerate(orbits):
    a=orb[0];b=int(jp[a]);pa,ta=tits(a);pb,tb=tits(b)
    perm=pa[pb];phase=tb*ta[pb]
    need(np.array_equal(perm, pb[pa]),'orthogonal root reflections commute')
    need(np.array_equal(perm[jp],jp[perm]),'actual J commutation')
    sq=phase*phase[perm]
    predicted=(-1)**(((roots[a]+roots[b])@roots.T//4)%2)
    need(np.array_equal(sq,predicted),'pair square exact character')
    square_profiles[int(np.sum(sq==-1))]+=1
    # Full nonzero root brackets, including opposite-root Cartan outputs.
    for u,v in nonzero:
        need(phase[sums[u,v]]*E[u,v]==phase[u]*phase[v]*E[perm[u],perm[v]],'all root bracket signs')
    for u in range(240):
        need(phase[u]*phase[neg[u]]*E[perm[u],neg[perm[u]]]==E[u,neg[u]],'opposite root bracket')
    eta=lift_phases(perm)
    char=sum(int(phase[i]==-1)<<j for j,i in enumerate(basis))
    expected=np.array([(-1)**((char&sum((int(n[j])%2)<<j for j in range(8))).bit_count()%2) for n in coords])
    need(np.array_equal(phase,eta*expected),'full240 field character dictionary')
    canonical_orders[2 if np.all(eta*eta[perm]==1) else 4]+=1
    solutions=[]
    for chi in range(256):
        signs=np.array([(-1)**((chi&sum((int(n[j])%2)<<j for j in range(8))).bit_count()%2) for n in coords],dtype=np.int8)
        ph=phase*signs
        if np.all(ph*ph[perm]==1) and np.array_equal(Jphase*ph[jp],ph*Jphase[perm]):solutions.append(chi)
    correction_counts[len(solutions)]+=1
    rows.append({'orbit':oi,'alpha2':roots[a].tolist(),'beta2':roots[b].tolist(),
                 'tits_to_generator_normalized_character_bits':char,
                 'involutive_J_commuting_input_characters':solutions})
    if oi%15==0:print('completed root-line field dictionary',oi,flush=True)

for p,h in PINS.items():need(sha256((REPO/p).read_bytes()).hexdigest()==h,'unchanged source')
result={'scope':'all60 actual Gaussian root lines; original eps_e8; full lattice equality from eight generators and common lattice action',
        'source_sha256':PINS,'root_lines':60,'nonzero_brackets_per_line':len(nonzero)+240,
        'tits_pair_square_minus_root_counts':dict(square_profiles),
        'canonical_order_counts':dict(canonical_orders),'J_commuting_involution_correction_counts':dict(correction_counts),
        'all_lattice_equality':'T_alpha T_Jalpha = S_r D_character; character fixed by8generator values',
        'rows':rows}
(HERE/'field-lift-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
