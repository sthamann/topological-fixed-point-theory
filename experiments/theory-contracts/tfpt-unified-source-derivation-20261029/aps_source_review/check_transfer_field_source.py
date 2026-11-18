"""Exact implication controls; no invented TFPT tensor or physical gap fit."""
from pathlib import Path
import itertools
import json
import sympy as s

checks = []
def need(condition, name):
    if not condition:
        raise AssertionError(name)
    checks.append(name)

def equal(a, b):
    return (a-b).applyfunc(s.simplify) == s.zeros(*a.shape)

I = s.eye(2)
X = s.Matrix([[0,1],[1,0]])
Y = s.Matrix([[0,-s.I],[s.I,0]])
Z = s.diag(1,-1)
A = [I/s.sqrt(2),X/s.sqrt(6),Y/s.sqrt(6),Z/s.sqrt(6)]
B = [(A[0]+A[1])/s.sqrt(2),(A[0]-A[1])/s.sqrt(2),A[2],A[3]]
def E(tensor, m):
    return sum((a*m*a.adjoint() for a in tensor), s.zeros(2))
for j, m in enumerate([I,X,Y,Z]):
    need(equal(E(A,m),E(B,m)), f'whole CP map unchanged on basis{j}')
    need(equal(E(A,m),m if j==0 else m/3), f'exact transfer eigenvalue{j}')
need(s.simplify(s.trace(A[0]*A[0].adjoint())/2)==s.Rational(1,2),'marked response A=1/2')
need(s.simplify(s.trace(B[0]*B[0].adjoint())/2)==s.Rational(1,3),'marked response B=1/3')

e = s.symbols('e',positive=True)
P = s.Matrix([[1-e,e],[e,1-e]])
need(s.simplify((P-s.eye(2)).det())==0,'Perron eigenvalue1 all e')
need(equal(P*s.Matrix([1,-1]),(1-2*e)*s.Matrix([1,-1])),'other eigenvalue1-2e')
units=[s.diag(1,0),s.diag(0,1)]
need(s.Matrix.hstack(*[s.Matrix(u*P*v).reshape(4,1) for u in units for v in units]).det()!=0,'Burnside matrix units all e nonzero')

# A rational member; analytic proof covers the entire 0<e<1/2 family.
a,b=s.Rational(12,13),s.Rational(5,13)
eps=b*b
lam=a*a-b*b
T=[s.Matrix([[a,b],[0,0]]),s.Matrix([[0,0],[b,a]])]
need(a*a+b*b==1,'tensor normalization')
need(equal(E(T,I),I),'MPS CP unital')
l=s.Matrix([[s.Rational(1,2),a*b],[a*b,s.Rational(1,2)]])
need(equal(sum((m.adjoint()*l*m for m in T),s.zeros(2)),l),'MPS faithful stationary density')
need(l.det()>0,'stationary density faithful')
products=[T[i]*T[j] for i,j in itertools.product(range(2),repeat=2)]
need(s.Matrix.hstack(*[p.reshape(4,1) for p in products]).det()!=0,'length2 injective')
words=list(itertools.product(range(2),repeat=3))
V=s.Matrix([list((T[i]*T[j]*T[k]).reshape(4,1)) for i,j,k in words])
need(V.rank()==4,'length3 image rank4')
h=s.eye(8)-V*(V.adjoint()*V).inv()*V.adjoint()
need(equal(h,h.adjoint()) and equal(h*h,h),'local parent orthogonal projector norm1')
L=6
bits=list(itertools.product(range(2),repeat=L))
amplitudes=[]
for w in bits:
    prod=I
    for i in w: prod=prod*T[i]
    amplitudes.append(s.trace(prod))
psi=s.Matrix(amplitudes)
norm=(psi.adjoint()*psi)[0]
need(norm==1+lam**L,'full periodic source partition1+lambda^L')
mag=s.Matrix([sum(1-2*i for i in w)*v for w,v in zip(bits,amplitudes)])
need((psi.adjoint()*mag)[0]==0,'magnetization trial orthogonal')
var=(mag.adjoint()*mag)[0]/norm
need(s.simplify(var-L*(1+lam)/(1-lam)*(1-lam**L)/(1+lam**L))==0,'exact whole ring variance')
for r in [1,2,3]:
    corr=sum(v*v*(1-2*w[0])*(1-2*w[r]) for w,v in zip(bits,amplitudes))/norm
    need(s.simplify(corr-(lam**r+lam**(L-r))/(1+lam**L))==0,f'whole ring correlation r{r}')
energy=s.S.Zero
for j in range(L):
    support=[(j+k)%L for k in range(3)]
    env=[k for k in range(L) if k not in support]
    for ew in itertools.product(range(2),repeat=L-3):
        pv=[];mv=[]
        for sw in words:
            w=[0]*L
            for k,i in zip(env,ew):w[k]=i
            for k,i in zip(support,sw):w[k]=i
            ix=bits.index(tuple(w));pv.append(psi[ix]);mv.append(mag[ix])
        need(equal(h*s.Matrix(pv),s.zeros(8,1)),f'actual local parent annihilates source{j},{ew}')
        hv=h*s.Matrix(mv)
        energy+=(hv.adjoint()*hv)[0]
rayleigh=s.factor(energy/(mag.adjoint()*mag)[0])
bound=9*(1-lam)/(1+lam)*(1+lam**L)/(1-lam**L)
need(rayleigh<=bound,'direct parent variational energy bound')
out={'result':'PASS','checks':len(checks),'checks_detail':checks,'cp_marked_responses':['1/2','1/3'],'transfer_lambda':str(lam),'direct_parent_trial_rayleigh':str(rayleigh),'parent_gap_bound':str(bound),'proof_scope':'Implication counterexamples only; original physical tensor and its gap not determined','complete_TFPT_solution':False,'ledger_changed':False}
Path(__file__).with_name('TRANSFER-FIELD-SOURCE.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='checks_detail'},indent=2))
