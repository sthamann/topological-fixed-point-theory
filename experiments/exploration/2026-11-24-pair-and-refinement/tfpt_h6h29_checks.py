#!/usr/bin/env python3
"""Exploratory checks after TFPT H6h28 (5537c03fe).
Independent model checks, not a full T4 or TOE certificate.
Requires numpy, sympy. No repository files are modified.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
import numpy as np
import sympy as sp

checks = []
def check(label, condition):
    assert condition, label
    checks.append(label)
    print('PASS:', label)

X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1, -1])
I = np.eye(2)
def kron(seq):
    a = np.array([[1]], complex)
    for x in seq:
        a = np.kron(a, x)
    return a
G = [kron([Z]*j + [p] + [I]*(4-j)) for j in range(5) for p in (X,Y)]
C = np.eye(32, dtype=complex)
for a in range(0,10,2):
    C = C @ G[a]
ix = [i for i in range(32) if i.bit_count() % 2 == 0]
M = [(C @ a)[np.ix_(ix,ix)] for a in G]
check('10 symmetric Spin(10) pairing matrices', all(np.array_equal(a,a.T) for a in M))
check('all 100 pairing Clifford relations', all(np.array_equal(a.conj().T@b+b.conj().T@a,2*(i==j)*np.eye(16)) for i,a in enumerate(M) for j,b in enumerate(M)))
check('orthogonality of pairing matrices', np.array_equal([[np.trace(a.conj().T@b) for b in M] for a in M],16*np.eye(10)))
Gamma = [np.block([[np.zeros((16,16)),a],[a.conj().T,np.zeros((16,16))]]) for a in M]
J = [-0.5j*Gamma[a]@Gamma[b] for a,b in combinations(range(10),2)]
K = J+[a/2 for a in Gamma]
check('Spin(11) one-particle Casimir 55/4',np.array_equal(sum(a@a for a in K),55/4*np.eye(32)))
P = np.eye(32)/2+J[0]
check('positive half of a Cartan generator has rank 16',np.array_equal(P@P,P) and np.trace(P)==16)
def slater_sq(A):
    return np.trace(P@A)**2+np.trace(P@A@(np.eye(32)-P)@A)
check('unique saturated Slater highest weight (8,0,0,0,0)', [np.trace(P@(-.5j*Gamma[a]@Gamma[a+1])) for a in range(0,10,2)]==[8,0,0,0,0])
check('highest Slater Casimir expectation 136',sum(slater_sq(a) for a in K)==136)
# Analytic representation-theoretic argument: all 32 spinor weights occur once;
# lambda1=8 saturates all 16 positive weights. The single saturated weight
# fixes the unique maximal irrep. Others have lambda1 <= 7.5. Interlacing
# B5 -> D5 bounds C_B-C_D <= lambda1(lambda1+9).
check('other irreps bounded strictly below first vector excitation', F(15,2)*F(33,2)<136-9)
rB=[F(9,2),F(7,2),F(5,2),F(3,2),F(1,2)]; rD=[4,3,2,1,0]
def cas(w,r):return sum(x*(x+2*y) for x,y in zip(w,r))
check('ground Casimir and first vector gap',cas([8,0,0,0,0],rB)==136 and cas([1,0,0,0,0],rD)==9)
check('branching dimension check',sum(comb(l+9,9)-(comb(l+7,9) if l>=2 else 0) for l in range(9))==comb(18,8)-comb(16,6)==35750)
check('single-fermion energy 29/2',136-cas([F(15,2)]+[F(1,2)]*4,rB)+cas([F(1,2)]*5,rD)==F(29,2))
# BCS spherical overlap
rng=np.random.default_rng(54);err=0
for _ in range(64):
    n=rng.normal(size=10);n/=np.linalg.norm(n);v=rng.normal(size=10);v/=np.linalg.norm(v)
    A=sum(x*a for x,a in zip(n,M));B=sum(x*a for x,a in zip(v,M))
    err=max(err,abs(np.linalg.det((np.eye(16)+A.conj().T@B)/2)-((1+n@v)/2)**8))
check('64 BCS overlap checks (degree eight)',err<1e-13); print('BCS max residual:',err)
# Exact non-Gaussian vector response
u=sp.symbols('u');poly=sp.gegenbauer(8,sp.Rational(9,2),u)
v=sp.diff(poly,u).subs(u,1)/poly.subs(u,1)
m4=(sp.diff(poly,u).subs(u,1)+3*sp.diff(poly,u,2).subs(u,1))/poly.subs(u,1)
check('zonal Spin(11) second and fourth moments',v==sp.Rational(68,5) and m4==442)
check('connected fourth moment',m4-3*v*v==-sp.Rational(2822,25))
a,b=sp.symbols('a b');sol=sp.solve([10*a+2*b-136*v,a-b-v],[a,b])
check('ordered tensor coefficients',sol=={a:sp.Rational(782,5),b:sp.Rational(714,5)} and sp.simplify((a+2*b).subs(sol))==m4)
E4=sp.simplify(v*v/9**3-(m4-v*v)/(9**2*20))
check('static quartic energy coefficient',E4==sp.Rational(17323,182250))
# independent SO(9) zonal sector recursion (9x9)
D,L=10,8
cs=[F((L-l)*(L+l+D-1)*(l+1)*(l+D-2),(2*l+D-2)*(2*l+D)) for l in range(L)]
A=np.diag(np.sqrt(np.array(cs,float)),1);A=A+A.T
H=np.diag([l*(l+8) for l in range(9)])
check('zonal generator spectrum',np.allclose(np.linalg.eigvalsh(A),np.arange(-8,9,2),atol=1e-12))
check('independent 9x9 moment check',abs((A@A)[0,0]-float(v))<1e-12 and abs(np.linalg.matrix_power(A,4)[0,0]-442)<1e-11)
h=.025;e=np.linalg.eigvalsh(H-h*A)[0]
check('independent static energy Taylor check',abs((e+float(v)/9*h*h)/h**4-float(E4))<1e-5)
# Native W crossing: P10 (symmetric spinor) x P_A (antisymmetric family).
# Unique equivariant image under Spin(10)xSU(4), W^*W=8 P on normalized antisymmetric pair space.
R=lambda A:sum(m@A.T@m.conj().T for m in M)/16
basis={0:[np.eye(16)],2:[],4:[]}
for n in (2,4):
    for ids in combinations(range(10),n):
        A=np.eye(32,dtype=complex)
        for i in ids:A=A@G[i]
        basis[n].append(A[np.ix_(ix,ix)])
for d,lam in [(0,F(5,8)),(2,F(-3,8)),(4,F(1,8))]:
    check(f'crossed Spin(10) degree {d}: {len(basis[d])} directions, eigenvalue {lam}',all(np.array_equal(R(A),float(lam)*A) for A in basis[d]))
check('Lorentz scalar and antisymmetric internal pair counts',136*10+120*6==2080 and 136*6+120*10==2016)
# SU4 epsilon bilinear pairs complementary 2-form labels.
fp=list(combinations(range(4),2))
def eps(ids):
    if len(set(ids))<4:return 0
    return (-1)**sum(ids[i]>ids[j] for i in range(4) for j in range(i+1,4))
Hodge=np.array([[eps(a+b) for b in fp] for a in fp]);light=[i for i,p in enumerate(fp) if 3 not in p];heavy=[i for i,p in enumerate(fp) if 3 in p]
check('real six has complementary isotropic light and heavy three-planes',np.array_equal(Hodge@Hodge,np.eye(6)) and not np.any(Hodge[np.ix_(light,light)]) and not np.any(Hodge[np.ix_(heavy,heavy)]) and np.linalg.matrix_rank(Hodge[np.ix_(light,heavy)])==3)
# Wilson parent in the exact J=infinity projected limit, not finite J interacting system.
N=12
for k in product([0,np.pi],repeat=3):
    w=sum(1-np.cos(x) for x in k);h=sum(np.sin(x)*m for x,m in zip(k,[X,Y,Z]));T=np.zeros((N,N-1));T[np.arange(N-1),np.arange(N-1)]=w;T[np.arange(1,N),np.arange(N-1)]=-1
    H=np.block([[np.kron(np.eye(N),h),np.kron(T,np.eye(2))],[np.kron(T.T,np.eye(2)),-np.kron(np.eye(N-1),h)]])
    ev=np.linalg.eigvalsh(H);v0=w**np.arange(N);v0/=np.linalg.norm(v0)
    check(f'Wilson corner {tuple(int(x>0) for x in k)}: two zero eigenvalues and exact null recursion',np.sum(abs(ev)<1e-10)==2 and np.linalg.norm(T.T@v0)<1e-12)
    print('  left/right weights:',v0[0]**2,v0[-1]**2,'chirality:',round(np.prod(np.cos(k))))
err=0
for _ in range(30):
    k=rng.uniform(-np.pi,np.pi,3);w=sum(1-np.cos(k));h=sum(np.sin(x)*m for x,m in zip(k,[X,Y,Z]));T=np.zeros((N,N-1));T[np.arange(N-1),np.arange(N-1)]=w;T[np.arange(1,N),np.arange(N-1)]=-1
    H=np.block([[np.kron(np.eye(N),h),np.kron(T,np.eye(2))],[np.kron(T.T,np.eye(2)),-np.kron(np.eye(N-1),h)]])
    e=np.linalg.eigvalsh(H);q=np.linalg.norm(np.sin(k));err=max(err,min(abs(e-q)),min(abs(e+q)))
check('30 generic momenta reproduce exact Weyl branch',err<1e-12)
# Gravity t=1 from area bilinearity, no microscopic selection claimed.
E=Hodge
def biv(P):return np.array([[np.linalg.det(P[np.ix_(a,b)]) for b in fp] for a in fp])
def elementary(X,k):return sum(np.linalg.det(X[np.ix_(i,i)]) for i in combinations(range(4),k)) if k else 1
err=0
for _ in range(24):
    P=np.eye(4)+.2*rng.normal(size=(4,4));K=np.eye(4)+.2*rng.normal(size=(4,4));D=biv(P)-biv(K)
    val=-.25*np.trace(E@D@E@D.T);Xr=np.linalg.solve(P,K);target=np.linalg.det(P)*(-1.5+.5*elementary(Xr,2)-1.5*np.linalg.det(Xr));err=max(err,abs(val-target))
check('24 general tetrad tests: quadratic area difference equals V_t=1',err<1e-12);print('area residual:',err)
a,b,t=sp.symbols('a b t');V=-2*a**4+2*a**3*b+2*a*b**3-2*b**4+t*(a-b)**4/2
check('conformal area subdivision identity',sp.simplify(V.subs(t,1)+sp.Rational(3,2)*(a*a-b*b)**2)==0)
D=[2,3,4,5];mid={(i,j):sp.Rational(1+D[i]*D[j],2) for i,j in fp}
products=[mid[(0,1)]*mid[(2,3)],mid[(0,2)]*mid[(1,3)],mid[(0,3)]*mid[(1,2)]]
check('area midpoint fails general simplicity',products==[sp.Rational(147,4),36,sp.Rational(143,4)])
# Frozen TFPT is not exact TM1 or mu-tau reflection.
phi=1/(6*np.pi)+3/(256*np.pi**4);x=phi*np.exp(-5/6);s12=1/3-phi/2
print('frozen:',{'phi':phi,'s13':x,'s12':s12,'TM1_s12_at_frozen_s13':(1-3*x)/(3*(1-x)),'e1_weight':(1-x)*(1-s12)})
check('exact TM1 incompatible with both frozen angle values',abs((1-x)*(1-s12)-2/3)>0.01)
# Native clock matrices (vector), from compiler-integral-triality-20260918/certificate.json.
# The physical identification of current-index space with lepton families is not assumed proven.
C=sp.Matrix([[-3,1,1,1,1,1,-1,-1],[-1,-1,-1,-1,-1,-1,-3,1],[-1,3,-1,-1,-1,-1,1,1],[-1,-1,3,-1,-1,-1,1,1],[-1,-1,-1,3,-1,-1,1,1],[-1,-1,-1,-1,3,-1,1,1],[-1,-1,-1,-1,-1,3,1,1],[1,1,1,1,1,1,-1,3]])/4
J=sp.diag(*([sp.Matrix([[0,-1],[1,0]])]*4))
P=sp.diag(*([1]*5+[0]*3));Q=sp.Matrix.hstack(*[sp.eye(8)[:,i] for i in (5,6,7)])
FC=Q.T*C*P*C.T*Q;FJ=Q.T*J*P*J.T*Q
u=sp.Matrix([0,1,1])/sp.sqrt(2);V=sp.Matrix.hstack(sp.Matrix([1,0,0]),sp.Matrix([0,1,-1])/sp.sqrt(2))
check('native clocks are orthogonal',C.T*C==sp.eye(8) and J.T*J==sp.eye(8))
check('native clock compressions have common fixed line',FC*u==sp.zeros(3,1) and FJ*u==sp.zeros(3,1) and FC.rank()==2)
check('native clock commutator exact norm 1/64',sp.trace((FC*FJ-FJ*FC).T*(FC*FJ-FJ*FC))==sp.Rational(1,64))
check('native two-plane block',V.T*FC*V==sp.Matrix([[sp.Rational(13,16),sp.sqrt(2)/16],[sp.sqrt(2)/16,sp.Rational(5,8)]]))
check('inverse Coxeter compression breaks common line',Q.T*C.T*P*C*Q*u!=sp.zeros(3,1))
omega=complex(-.5,np.sqrt(3)/2);Fourier=np.array([[omega**(i*j)/np.sqrt(3) for j in range(3)] for i in range(3)])
check('common line gives TM1 magnitudes with conditional cyclic charged basis',np.allclose(abs(Fourier.conj().T@np.array(u,float).ravel())**2,[2/3,1/6,1/6]))
for w in [0,1/8,1]:
    B=np.array(V.T*(FC+w*FJ)*V,float);theta=.5*np.arctan2(2*B[0,1],B[0,0]-B[1,1]);print('clock weight',w,'theta(deg)',np.degrees(theta),'s13',np.sin(theta)**2/3)
print(f'\n{len(checks)}/{len(checks)} grouped checks passed. No T4, gravity UV or TOE claim.')
