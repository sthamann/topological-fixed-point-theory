#!/usr/bin/env python3
"""Exact original unbounded-charge word action; no finite charge cutoff.
The finite 256-state construction checks only the commutator algebra.
"""
import ast, hashlib, json
from pathlib import Path
import sympy as sp
ROOT=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
original=ROOT/'experiments/theory-contracts/local-charge-transport-round20/checker.py'
# Execute only the original Gram assignments, not its recursive prerequisite suite.
tree=ast.parse(original.read_text());main=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main')
selected=[]
for n in main.body:
 if isinstance(n,ast.Assign) and len(n.targets)==1 and isinstance(n.targets[0],ast.Name) and n.targets[0].id in ('eye','B','G'):
  selected.append(n)
ns={'sp':sp};exec(compile(ast.Module(body=selected,type_ignores=[]),str(original),'exec'),ns)
G=ns['G'];g=[[int(x) for x in row] for row in G.tolist()]
def rank2(M):
 a=[[int(x)%2 for x in row] for row in M];r=0
 for c in range(len(a[0])):
  pivot=next((k for k in range(r,len(a)) if a[k][c]),None)
  if pivot is None:continue
  a[r],a[pivot]=a[pivot],a[r]
  for k in range(len(a)):
   if k!=r and a[k][c]:a[k]=[x^y for x,y in zip(a[k],a[r])]
  r+=1
 return r
zero=(0,)*8
units=[tuple(int(i==j) for i in range(8)) for j in range(8)]
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def neg(a):return tuple(-x for x in a)
def beta(a,b):
 return (sum(a[i]*b[i]*(g[i][i]//2) for i in range(8))+sum(a[i]*b[j]*g[i][j] for i in range(8) for j in range(i)))%2
def T(profile,x,p,adj=False):
 out=list(profile)
 if adj:out[x]=add(profile[x],neg(p));phase=(-1)**beta(p,out[x])
 else:phase=(-1)**beta(p,profile[x]);out[x]=add(profile[x],p)
 return tuple(out),phase
def W(profile,x,y,p,adj=False):
 if adj:return W(profile,y,x,p)
 a,s=T(profile,x,p,True);b,t=T(a,y,p);return b,s*t
def word(profile,ops):
 sign=1
 for x,y,p,adj in ops:profile,s=W(profile,x,y,p,adj);sign*=s
 return profile,sign
checks=[]
def ck(name,ok):
 if not ok:raise AssertionError(name)
 checks.append(name)
ck('original Gram even unimodular and rank8 mod2',G.det()==1 and rank2(g)==8 and all(g[i][i]%2==0 for i in range(8)))
profiles=[(zero,zero,zero)]
for p in units:profiles.extend([(p,neg(p),zero),(neg(p),zero,p)])
profiles.append(((20,-7,9,0,0,-14,3,1),(-20,7,-9,0,0,14,-3,-1),zero))
for n in profiles:
 for i,p in enumerate(units):
  for j,q in enumerate(units):
   lhs,ls=word(n,[(1,2,q,False),(0,1,p,False)])
   rhs,rs=word(n,[(0,1,p,False),(1,2,q,False)])
   ck(f'profile {profiles.index(n)} A{i}B{j} signed action',lhs==rhs and ls==(-1)**g[i][j]*rs)
   for x,y in [(0,1),(1,2)]:
    a,sa=word(n,[(x,y,q,False),(x,y,p,False)])
    b,sb=word(n,[(x,y,p,False),(x,y,q,False)])
    ck(f'same-bond commuting {profiles.index(n)} {x} {i} {j}',a==b and sa==sb)
C=[[0]*8+row for row in g]+[row+[0]*8 for row in g]
ck('two-bond form rank16',rank2(C)==16)
# Sharp finite representation of pairwise relations only: Ai=X_i, Bj=Z^(G.column(j)).
for i in range(8):
 for j in range(8):
  for t in range(256):
   b=lambda t:(-1)**(sum(((t>>k)&1)*g[k][j] for k in range(8))%2)
   ck(f'finite sign {i},{j},{t}',b(t)==(-1)**g[i][j]*b(t^(1<<i)))
ranks={}
for N in range(2,10):
 # Neutral vectors (ei at x, -ei at last site), x=0,...,N-2.
 form=[[((int(x==y)+1)*g[i][j])%2 for y in range(N-1) for j in range(8)] for x in range(N-1) for i in range(8)]
 rank=rank2(form);expected=8*(N-1 if N%2 else N-2)
 ck(f'N{N} neutral form rank',rank==expected)
 ranks[N]={'rank':rank,'minimal_dimension_for_commutator_representation':2**(rank//2)}
# Existing negative loop preserved, not called a new result.
i,j=next((i,j) for i in range(8) for j in range(8) if g[i][j]%2)
closed,phase=word(profiles[0],[(1,2,units[j],True),(0,1,units[i],True),(1,2,units[j],False),(0,1,units[i],False)])
ck('actual neutral closed word has phase minus one',closed==profiles[0] and phase==-1)
loop_i,loop_j=i,j
# Same parity is not the same future response: exact original cubic D.
import itertools
L=3;N=L**3;sites=list(itertools.product(range(L),repeat=3));idx={v:i for i,v in enumerate(sites)}
x=idx[(0,0,0)];y=idx[(1,0,0)]
edges=[]
for i,v in enumerate(sites):
 for axis in range(3):
  w=list(v);w[axis]=(w[axis]+1)%L;edges.append((i,idx[tuple(w)]))
p=units[7];r=[zero]*N;r[x]=neg(p);r[y]=p
nu=sp.Symbol('nu',nonnegative=True)
def e(v):return (sp.Matrix(v).T*G*sp.Matrix(v))[0]/2
def D(profile):return sum(e(v) for v in profile)/N+nu*sum(e(add(profile[i],neg(profile[j]))) for i,j in edges)
Dr=sp.simplify(D(r));ck('native adjacent profile exact D coefficient',Dr==sp.Rational(2,N)+14*nu)
response=[]
for k in (-3,-1,0,1,4):
 n=[tuple(-2*k*a for a in v) for v in r];nr=[add(a,b) for a,b in zip(n,r)]
 ck(f'n{k} same zero parity and neutral',all(a%2==0 for v in n for a in v) and all(sum(v[i] for v in n)==0 for i in range(8)))
 ck(f'n{k} exact charge energy response',sp.expand(D(nr)-D(n)-(1-4*k)*Dr)==0)
 ck(f'n{k} exact scalar diagonal response',all(e(nr[i])-e(n[i])==(1-4*k)*e(r[i]) for i in range(N)))
 response.append({'k':k,'response_slope_multiplier':1-4*k})
result={'status':'PASS','checks':len(checks),'original_source':str(original.relative_to(ROOT)),'sha256':hashlib.sha256(original.read_bytes()).hexdigest(),'gram':g,'mod2_rank':8,'two_bond_rank':16,'two_bond_minimal_finite_unitary_dimension':256,'profiles_in_exact_checks':len(profiles),'same_parity_response':{'D_r':str(Dr),'examples':response,'scope':'Native response derivatives, no full dynamics simulation'},'signed_neutral_loop':{'p_index':loop_i,'q_index':loop_j,'closed':True,'phase':phase},'neutral_rank_checks':ranks,'scope':['Original unbounded integer profiles, exact phases; finite samples check the symbolic identities proved in PROOF.md.','256-state matrices witness sharpness for the abstract pairwise commutator relations only; they are not a truncation of the original Hamiltonian.','No exclusion of a partial readout, nonhomomorphic CP reduction or larger full source.','No physical particle or spacetime dimension inferred.']}
Path(__file__).with_name('neutral-transport-result.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('gram','scope','neutral_rank_checks')},indent=2))
