"""Read original Q4 effects, retain both ordered legs, and compare kernels."""
from pathlib import Path
from itertools import product,combinations
import hashlib
import sympy as s
r=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
q=r/'_newest2/TFPT_Gesamtdokumentation2_20260927.md'
src=q.read_text().split('Originalmitglied: `TFPT_Rekursive_Bindung_20260926/audit.py`',1)[1].split('```python\n',1)[1].split('```',1)[0]
assert '4a19f8145d5ccaef08bfdc5893f64484feae630cf61e49639df52f89f56528c4' in [hashlib.sha256((src.rstrip()+'\n'*j).encode()).hexdigest() for j in range(4)]
state={'__file__':str(q)}
exec(src.split('# Exact complete 35 x 35 two-module model.',1)[0],state)
F,G,K=[state[x] for x in ('Fnum','Gs','Recovery')]
I=s.eye(25)
vec=lambda x:s.Matrix(25,1,list(x))
unvec=lambda x:s.Matrix(5,5,list(x))
adj=lambda x:G.inv()*x.H*G
# Normalized projectors as actual endomorphisms in the native code Gram.
f=[F[j,:].T/4 for j in range(10)]
v=[G.inv()*x for x in f]
Pr=[2*a*b.T for a,b in zip(v,f)]
assert all(x*x==x and adj(x)==x for x in Pr)
assert sum(Pr,s.zeros(5))==2*s.eye(5)
# In nonorthonormal End(C5), use K(X)=sum_j Pi_j tr(Pi_j X)/2.
KR=lambda X:sum((P*s.trace(P*X)/2 for P in Pr),s.zeros(5))
assert all(K*vec(X*G.inv())==vec(KR(X)*G.inv())
           for X in (s.eye(5)[:,a]*s.eye(5)[b,:] for a,b in product(range(5),repeat=2)))
flat=lambda xs:s.Matrix.hstack(*[vec(x) for x in xs])
pairs=list(combinations(range(10),2))
assert flat(Pr).rank()==10
assert flat([Pr[i]*Pr[j] for i,j in product(range(10),repeat=2)]).rank()==25
assert flat([Pr[i]*Pr[j]-Pr[j]*Pr[i] for i,j in pairs]).rank()==10
assert flat(Pr+[Pr[i]*Pr[j]+Pr[j]*Pr[i] for i,j in pairs]).rank()==15
print('Exact single effect rank 10; retaining all ordered two-effect words gives rank 25.')
print('Two-effect differences span the entire antisymmetric 10; symmetric sums together with effects span the full symmetric 15.')
# Build a silent symmetric five-mode with the original spectral projector.
met=s.kronecker_product(G,G)
om=vec(G.inv())
B=om*om.T*met/5
P9=s.Rational(9,4)*(K-B)
SW=s.Matrix(25,25,lambda a,b:int(a==5*(b%5)+b//5))
P5=(I+SW)/2-B-P9
for i,j in pairs:
    A=Pr[i]*Pr[j]+Pr[j]*Pr[i]
    X=unvec(P5*vec(A*G.inv()))*G
    if X!=s.zeros(5):break
assert adj(X)==X and KR(X)==s.zeros(5) and all(s.trace(P*X)==0 for P in Pr)
assert s.trace(X*Pr[i]*Pr[j])==s.trace(X*Pr[j]*Pr[i])!=0
print('Silent symmetric mode: original effect indices',i,j,'open trace',s.trace(X*Pr[i]*Pr[j]),'norm^2',s.trace(X*X))
i,j=0,1
Y=s.I*(Pr[i]*Pr[j]-Pr[j]*Pr[i])
assert adj(Y)==Y and KR(Y)==s.zeros(5) and all(s.trace(P*Y)==0 for P in Pr)
a=s.trace(Y*Pr[i]*Pr[j]);b=s.trace(Y*Pr[j]*Pr[i])
assert a==-b!=0
print('Silent orientation mode: original effect indices',i,j,'opposite open traces',a,b)
for U in (s.eye(5),-s.eye(5)):
    assert all(U*P*U.inv()==P for P in Pr)
print('Projector word algebra is complete on the code, but a central lift sign remains invisible under adjoint action.')
