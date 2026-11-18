from pathlib import Path
from itertools import combinations, permutations
import sympy as s
import json, hashlib, contextlib, io
base=Path(__file__).parent
code=(base/'native-audit-source.py').read_text()
assert hashlib.sha256((code.rstrip()+'\n').encode()).hexdigest()=='4a19f8145d5ccaef08bfdc5893f64484feae630cf61e49639df52f89f56528c4'
ns={'__file__':str(base/'native-audit-source.py')}
with contextlib.redirect_stdout(io.StringIO()): exec(code.split('# Exact complete 35 x 35 two-module model.',1)[0],ns)
G,Fnum,Knative=[ns[x] for x in ['Gs','Fnum','Recovery']]
# Rational native effect frame with metric G.
Un=G.inv()*Fnum.T/4 # unit frame is sqrt(2)*Un
Gn=2*Un.T*G*Un
edges=list(combinations(range(6),2))
A=s.Matrix(15,15,lambda i,j:int(set(edges[i]).isdisjoint(edges[j])))
M=(A+3*s.eye(15))/9
B=A+s.eye(15);T=(-B**2+15*B+34*s.eye(15))/90
assert T*T==M
roots=[s.eye(6)[:,a]-s.eye(6)[:,b] for a,b in edges]
Fs=[r*r.T for r in roots];F=s.Matrix.hstack(*(x.reshape(36,1) for x in Fs))
Gram=F.T*F
assert Gram==3*s.eye(15)+s.ones(15)-A and Gram.det()!=0
# Balanced split frame, representative with sixth component +1.
signs=[s.Matrix([1 if i in (*e,5) else -1 for i in range(6)]) for e in combinations(range(5),2)]
V=s.Matrix.hstack(*signs);Gb=V.T*V/6
Pn=[v*v.T/6 for v in signs]
K=s.Rational(1,2)*sum((s.kronecker_product(p,p) for p in Pn),s.zeros(36))
assert K*F==F*M
P=s.eye(6)-s.ones(6)/6
assert sum(Pn,s.zeros(6))==2*P
# Find exact signed permutation matching original Bell Gram, no eigenvector fitting.
dn=[1]+[int(3*Gn[0,j]) for j in range(1,10)]
db=[1]+[int(3*Gb[0,j]) for j in range(1,10)]
cn=[[int(3*Gn[i,j])*dn[i]*dn[j] for j in range(10)] for i in range(10)]
cb=[[int(3*Gb[i,j])*db[i]*db[j] for j in range(10)] for i in range(10)]
def match(p):
 if len(p)==10:return p
 i=len(p)
 for j in range(10):
  if j not in p and all(cn[i][k]==cb[j][p[k]] for k in range(i)):
   r=match(p+[j])
   if r:return r
perm=match([0]); assert perm
sgn=[dn[i]*db[perm[i]] for i in range(10)]
Vm=s.Matrix.hstack(*(sgn[i]*signs[perm[i]] for i in range(10)))
assert Vm.T*Vm/6==Gn
# L maps the native, nonorthonormal G coordinates to standard R6.
# Each actual sqrt(2)*Un column maps to Vm/sqrt(6).
# Tight frame formula: L = Vm Un^T G / sqrt(12).
L=Vm*Un.T*G/s.sqrt(12)
assert (L.T*L-G).applyfunc(s.simplify)==s.zeros(5)
assert (L*Un-Vm/s.sqrt(12)).applyfunc(s.simplify)==s.zeros(6,10)
LL=s.kronecker_product(L,L)
assert (K*LL-LL*Knative).applyfunc(s.simplify)==s.zeros(36,25)
# Full Szegedy walk using integer-rational projectors on 45 incidences.
def matchings(v):
 if not v:yield ();return
 a=v[0]
 for b in v[1:]:
  for rest in matchings([x for x in v if x not in (a,b)]):yield ((a,b),)+rest
synt=list(matchings(list(range(6))))
N=s.Matrix(15,15,lambda i,j:int(edges[i] in synt[j]))
assert N*N.T==A+3*s.eye(15)
inc=[(i,j) for i in range(15) for j in range(15) if N[i,j]]
D=s.Matrix(45,15,lambda e,i:int(inc[e][0]==i))
S=s.Matrix(45,15,lambda e,j:int(inc[e][1]==j))
PD=D*D.T/3;PS=S*S.T/3
U=(2*PS-s.eye(45))*(2*PD-s.eye(45))
assert U.T*U==s.eye(45)
x=s.symbols('x')
assert s.expand(U.charpoly(x).as_expr()-(x-1)**17*(x+1)**10*(x*x+s.Rational(2,9)*x+1)**9)==0
# A unitary compression is NOT the Markov transfer: D^* U D = 2M-I.
assert D.T*U*D/3==2*M-s.eye(15)
# Character mismatch for natural unsigned vs exterior actions.
words=[t for k in [0,2,4] for t in combinations(range(5),k)]
p=[1,0,2,3,4]; unsigned=fermion=0
for w in words:
 z=[p[i] for i in w]
 if tuple(sorted(z))==w:
  unsigned+=1; fermion+=(-1)**sum(z[i]>z[j] for i in range(len(z)) for j in range(i+1,len(z)))
assert (unsigned,fermion)==(8,0)
# Petersen Ihara claims, orientations retained in trace.
pairs=list(combinations(range(5),2));AP=s.Matrix(10,10,lambda i,j:int(set(pairs[i]).isdisjoint(pairs[j])))
arcs=[(i,j) for i in range(10) for j in range(10) if AP[i,j]]
H=s.Matrix(30,30,lambda i,j:int(arcs[i][1]==arcs[j][0] and arcs[i][0]!=arcs[j][1]))
assert s.trace(H**5)==s.trace(H**6)==120
res={'native_checks_reused':len(ns['checks']),'native_signed_permutation':perm,'native_signs':sgn,'native_isometry_exact':True,'KF_equals_FM':True,'F_Gram':'3I+J-A; positive eigenvalues 12,2,6','T_squared_equals_M':True,'walk_characteristic_polynomial':str(s.factor(U.charpoly(x).as_expr())),'compression':'D_normalized^* U D_normalized = 2M-I, not M','characters':[unsigned,fermion],'Petersen_trace_H5':120,'Petersen_trace_H6':120,'physical_P1_projectors_derived':False,'F1_complete':False}
(base/'results.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
