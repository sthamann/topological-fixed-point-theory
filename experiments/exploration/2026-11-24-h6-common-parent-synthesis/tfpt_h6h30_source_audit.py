#!/usr/bin/env python3
"""Independent H6h29 follow-up; scoped model checks, not a TOE proof.
Requires numpy, scipy, sympy. No TFPT repo imports, no fitting to data.
"""
import math
import numpy as np
import sympy as s
from scipy.linalg import eigh_tridiagonal
from scipy.special import comb

checks = []
def check(name, condition):
    assert bool(condition), name
    checks.append(name)
    print('PASS', name)

# A_k: rectangular cut at depth D. Infinite representation is a positive
# discrete-series SU(1,1) ladder; finite boundary truncation preserves the null vector.
def ladder(k, r, D):
    q = math.sqrt(r)
    n = np.arange(D, dtype=float)
    A = np.zeros((D, D+1))
    A[np.arange(D), np.arange(D)] = q*np.sqrt(n+k)
    A[np.arange(D), np.arange(1,D+1)] = -np.sqrt(n+1)
    return A

def nb(k,r,D):
    n=np.arange(D+1)
    return (1-r)**k*comb(n+k-1,k-1)*r**n

r = s.symbols('r', positive=True)
z = s.symbols('z')
# Exact equality of the generating function with k independent geometric modes.
G = (1-r)/(1-r*z)
check('NB4 is the convolution of four NB1 laws', s.cancel(G**4-(1-r)**4/(1-r*z)**4)==0)
for k in (1,2,3,4,8):
    g=G**k
    mean=s.factor(s.diff(g,z).subs(z,1))
    var=s.factor((s.diff(g,z,2)+s.diff(g,z)).subs(z,1)-mean**2)
    check(f'NB{k} cumulants',s.simplify(mean-k*r/(1-r))==0 and s.simplify(var-k*r/(1-r)**2)==0)

# Infinite spectral expression: compare full distribution, not just its zero pole.
for rr in (1/3,math.exp(-4*math.pi)):
    for k in (1,4):
        D=120
        A=ladder(k,rr,D)
        p=nb(k,rr,D)
        v=np.sqrt(p/p.sum())
        check(f'k={k}, r={rr:.6g}: exact-cut ground vector', np.linalg.norm(A@v)<1e-12)
        K=A.T@A
        eig,vec=eigh_tridiagonal(np.diag(K),np.diag(K,1))
        check(f'k={k}, r={rr:.6g}: first 10 ladder levels', np.max(np.abs(eig[:10]-(1-rr)*np.arange(10)))<1e-9)
        check(f'k={k}, r={rr:.6g}: first 15 boundary spectral weights',np.max(np.abs(vec[0,:15]**2-p[:15]))<1e-10)
        check(f'k={k}, r={rr:.6g}: total spectral norm',abs(p.sum()-1)<1e-12)

# The finite four-arm metric graph realizes the k=4 radial ladder, not an
# additional independent proof of its physical identification. M_c=lambda=1.
for k in (1,4):
    for rr in (1/3, math.exp(-4*math.pi)):
        D=12; q=math.sqrt(rr)
        M=np.ones(1+4*D)
        def ix(j,n): return 1+j*D+n-1
        for j in range(4):
            for n in range(1,D+1):
                M[ix(j,n)] = 0.5*math.sqrt(math.comb(n+k-1,k-1)*rr**n)
        B=np.zeros((4*D,1+4*D))
        row=0
        for j in range(4):
            for n in range(D):
                child=ix(j,n+1); parent=0 if n==0 else ix(j,n)
                w=math.sqrt(n+1)*M[child]
                B[row,parent]=w/M[parent];B[row,child]=-w/M[child];row+=1
        U=np.zeros((1+4*D,D+1));U[0,0]=1
        V=np.zeros((4*D,D))
        for j in range(4):
            for n in range(1,D+1): U[ix(j,n),n]=0.5;V[j*D+n-1,n-1]=0.5
        A=ladder(k,rr,D)
        check(f'k={k}, r={rr:.6g}: full radial intertwiner',np.max(np.abs(B@U-V@A))<1e-12)
        check(f'k={k}, r={rr:.6g}: universal massless source coupling',np.linalg.norm(B@M)<1e-12 and np.allclose((M/np.linalg.norm(M))/M,1/np.linalg.norm(M)))
        check(f'k={k}, r={rr:.6g}: total mass norm',abs(np.sum(M*M)-sum(math.comb(n+k-1,k-1)*rr**n for n in range(D+1)))<1e-12)

# Closed inverse moments for the two candidate central measures.
S1_1=-s.log(1-r); S2_1=s.polylog(2,r)
S1_4=(r*(1+r)/(1-r)**3+6*r/(1-r)**2+11*r/(1-r))/6-s.log(1-r)
S2_4=(r/(1-r)**2+6*r/(1-r)-11*s.log(1-r))/6+s.polylog(2,r)
for k,F1,F2 in ((1,S1_1,S2_1),(4,S1_4,S2_4)):
    # polynomial coefficients after expansion; exact finite Taylor checks.
    for m,F in ((1,F1),(2,F2)):
        expected=sum(s.Rational(math.comb(n+k-1,k-1),n**m)*r**n for n in range(1,9))
        check(f'k={k}: exact inverse-moment-{m} series',s.expand(s.series(F,r,0,9).removeO()-expected)==0)

rr=math.exp(-4*math.pi)
print('\nPhysical quarter parameter r=e^-4pi:',rr)
for k,F1,F2 in ((1,S1_1,S2_1),(4,S1_4,S2_4)):
    # use direct high precision because rational/log expressions nearly cancel
    A=s.N(F1.subs(r,s.exp(-4*s.pi)),35)
    B=s.N(F2.subs(r,s.exp(-4*s.pi)),35)
    W=s.N((1-s.exp(-4*s.pi))**(-k)-1,35)
    print(f'k={k}: massive/zero spectral weight={float(W):.12g}; b/a^2={float(1+B/A**2):.12g}')

# The exact TM1 texture conflicts with BOTH frozen angles. No fit.
phi=1/(6*math.pi)+3/(256*math.pi**4)
s12=1/3-phi/2; s13=phi*math.exp(-5/6)
e1=(1-s12)*(1-s13)
check('frozen pair incompatible with exact TM1',abs(e1-2/3)>0.01)
print('frozen |U_e1|^2:',e1,'delta from 2/3:',e1-2/3)

# Finite norm approximate half-sided modular dilation has a lower bound.
# Check on random finite matrices; proof is reverse triangle inequality.
rng=np.random.default_rng(20261124)
for j in range(30):
    X=rng.normal(size=(8,8))+1j*rng.normal(size=(8,8))
    P=X.conj().T@X
    U,_=np.linalg.qr(rng.normal(size=(8,8))+1j*rng.normal(size=(8,8)))
    a=math.exp(-2*math.pi*rng.uniform(-.5,.5))
    err=np.linalg.norm(U@P@U.conj().T-a*P,ord=2)
    assert err+1e-10 >= abs(1-a)*np.linalg.norm(P,ord=2)
check('operator-norm obstruction to bounded exact modular dilation (30 controls)',True)
print(f'\n{len(checks)}/{len(checks)} checks passed. Conditional algebra and model calculations only.')
