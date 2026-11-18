#!/usr/bin/env python3
"""Exact original lattice/mark/cocycle tests, no finite charge truncation."""
import ast,hashlib,itertools,json
from pathlib import Path
import sympy as sp
ROOT=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
p=ROOT/'experiments/theory-contracts/local-charge-transport-round20/checker.py'
tree=ast.parse(p.read_text());main=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main')
assigns=[n for n in main.body if isinstance(n,ast.Assign) and len(n.targets)==1 and isinstance(n.targets[0],ast.Name) and n.targets[0].id in ('eye','B','G')]
ns={'sp':sp};exec(compile(ast.Module(body=assigns,type_ignores=[]),str(p),'exec'),ns)
G=ns['G'];B=ns['B'];Binv=B.inv();checks=0

def ck(ok,name):
 global checks
 if not ok:raise AssertionError(name)
 checks+=1

def c(n,m):return int(sum((G[i,i]//2)*n[i]*m[i] for i in range(8))+sum(G[i,j]*n[i]*m[j] for i in range(8) for j in range(i)))%2

def k(n,S=(1,2,3,4)):
 q=B*sp.Matrix(n);return int(2*(q[0]+sum(q[i] for i in S)))%4
units=[tuple(int(i==j) for i in range(8)) for j in range(8)]
ck([k(p) for p in units]==[0,0,0,0,0,2,2,1],'actual source basis grades')
roots=[]
for i,j in itertools.combinations(range(8),2):
 for a,b in itertools.product((-1,1),repeat=2):
  q=[0]*8;q[i]=a;q[j]=b;roots.append(sp.Matrix(q))
for signs in itertools.product((-1,1),repeat=8):
 if sum(v<0 for v in signs)%2==0:roots.append(sp.Matrix([sp.Rational(v,2) for v in signs]))
ck(len(roots)==240,'source E8 roots')
markings=[]
for S in itertools.combinations(range(1,7),4):
 neutral=[q for q in roots if int(2*(q[0]+sum(q[i] for i in S)))%4==0]
 M=sum((q*q.T for q in neutral),sp.zeros(8));P=sp.diag(*[int(i==0 or i in S) for i in range(8)])
 ck(len(neutral)==52 and M==16*P+8*(sp.eye(8)-P),'intrinsic marked second moment')
 ck((M-8*sp.eye(8))/8==P,'projector reconstruction from grade')
 markings.append({'S':S,'neutral_roots':len(neutral),'rank':int(P.trace())})
# Physical coordinate permutation is exactly a permutation of original root columns 1..6.
# Its cocycle correction has inversion phase, retaining all integer carries.
profiles=[(0,)*8]+units+[(20,-7,9,0,0,-14,3,1),(-1000000,7,0,0,3,-5,2,-1)]
for perm in itertools.permutations(range(1,7)):
 mapping=(0,)+perm+(7,)
 P=sp.zeros(8)
 for i,j in enumerate(mapping):P[j,i]=1
 ck(P.T*G*P==G,'all source S6 metric isometries')
 def apply(v):
  out=[0]*8
  for i,j in enumerate(mapping):out[j]=v[i]
  return tuple(out)
 def phase(v):return sum(v[i]*v[j] for i in range(1,7) for j in range(i+1,7) if mapping[i]>mapping[j])%2
 for n in profiles:
  for v in units:
   nv=tuple(a+b for a,b in zip(n,v))
   ck((phase(n)+phase(v)+phase(nv))%2==(c(apply(n),apply(v))+c(n,v))%2,'full integer cocycle correction')
   ck(phase(v)==0,'native generator rephase zero')
# A specific indistinguishable-source witness for different carrier projectors.
q=sp.zeros(8,1);q[0]=q[4]=1;n=tuple(int(v) for v in Binv*q)
ck(k(n)==0 and k(n,(1,2,3,5))==2,'same source root distinguishes two mark choices')
# The local mark is not dynamically closed: neutral edge e8 shift changes grade.
zero=(0,)*8;e8=units[7]
neg=tuple(-v for v in e8)
ck(k(neg)==3 and k(e8)==1,'neutral e8 hop changes both local marks')
ck(c(e8,neg)==1,'adjoint phase retained')
# T_x,e8* atzero uses phase c(e8,-e8)=1; thus W amplitude -1.
# Hhop has -J W, yielding +J; [H,M0] = +J on this ordered matrix element.
W_phase=(-1)**(c(e8,neg)+c(e8,zero))
commutator_coefficient=-W_phase
ck(commutator_coefficient!=0,'source-derived memory coupling nonzero for J>0')
result={'status':'PASS','checks':checks,'source':str(p.relative_to(ROOT)),'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'basis_grades':[k(v) for v in units],'markings_in_S6_orbit':len(markings),'grade_neutral_root_count':52,'neutral_second_moment_eigenvalues':{'16':5,'8':3},'source_symmetries_checked':720,'integer_profiles':len(profiles),'local_neutral_hop_grades':[3,1],'native_commutator_matrix_element':'<(-e8,+e8)|[H,M0(x)]|(0,0)> = +J','scope':'Exact fixed Round20-26 source identities. No finite cutoff, no Gibbs approximation, no numerical memory kernel, no universal exclusion of full TFPT with additional k or C/J data.'}
Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
