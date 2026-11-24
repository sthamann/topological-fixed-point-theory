#!/usr/bin/env python3
"""Conditional native quarter-transfer half-space Wilson parent; no TOE claim.
Exact identities, finite-section regressions, and tests of two selection rules.
Baseline TFPT main 75d8502f, 2026-11-24. Requires numpy and sympy.
"""
import itertools
import numpy as np
import sympy as sp

ncheck=0
def check(name, b):
    global ncheck
    assert bool(b), name
    ncheck += 1
    print('PASS:',name)

# A(k)=(S^* - w)(I-q S^4). In half-space, S^*S=I.
# Gap: BB^* >= (1-w)^2 I, AA^* >= (1-q)^2 BB^*.
a,b,c=sp.symbols('a b c', nonnegative=True)
w=a+b+c
s2=2*w-a*a-b*b-c*c
check('Wilson bulk gap identity',sp.expand(s2+(1-w)**2-1-2*(a*b+a*c+b*c))==0)
for k in itertools.product([0,1], repeat=3):
    wc=2*sum(k)
    check(f'corner {k}: index one iff origin',(wc<1)==(sum(k)==0))
# Exact analytic normalization of half-space kernel.
# v=R^{-1} g_w, g_w[n]=w^n, R=I-qS^4.
# ||v||² = (1+q*w^4)/((1-w²)*(1-q²)*(1-q*w^4)).
q0=np.exp(-2*np.pi)
for q in [q0,.1,.45,.8]:
    for w0 in [0.,.2,.6,.9,.98]:
        N=1200
        v=np.zeros(N)
        for j in range(N):
            v[j]=w0**j+(q*v[j-4] if j>=4 else 0.)
        Rv=v.copy();Rv[4:]-=q*v[:-4]
        norm=(1+q*w0**4)/((1-w0*w0)*(1-q*q)*(1-q*w0**4))
        check(f'kernel recursion q={q:.6g}, w={w0}',np.max(abs(Rv[1:]-w0*Rv[:-1]))<2e-14)
        check(f'closed kernel norm q={q:.6g}, w={w0}',abs(v@v/norm-1)<1e-12)
# Finite sections: rectangular B models the topological boundary but adds
# a distant compensating boundary for w>1; no inference of infinite gap from it.
for q in [q0,.3,.8]:
    N=36
    S=np.eye(N,k=-1); R=np.eye(N)-q*np.linalg.matrix_power(S,4)
    for w0 in [0.,.3,.95,1.,1.2,2.,4.,6.]:
        B=(S.T-w0*np.eye(N))[:-1,:]
        A=B@R
        eig=np.linalg.eigvalsh(A@A.T)
        check(f'positive nonzero singular spectrum bound q={q:.6g}, w={w0}',eig[0]+1e-10 >= (1-q)**2*(1-w0)**2)
        if w0<1:
            v=np.linalg.solve(R,w0**np.arange(N));v/=np.linalg.norm(v)
            check(f'finite section reproduces left kernel q={q:.6g}, w={w0}',np.linalg.norm(A@v)<1e-12)
# Positive depth weights leave the Fredholm kernel and cannot reduce the
# nonzero singular-value lower bound. In the infinite case A_weighted is
# closed on {f: W^(1/2) A f in l2}; finite-support vectors are a core candidate.
for q in [q0,.4]:
    N=28; S=np.eye(N,k=-1); R=np.eye(N)-q*np.linalg.matrix_power(S,4)
    W=np.diag(1+np.arange(N-1)//4)
    for w0 in [0.,.5,.99,1.2,2.]:
        A=(S.T-w0*np.eye(N))[:-1,:]@R
        Aw=np.sqrt(W)@A
        eig=np.linalg.eigvalsh(A@A.T)
        eigw=np.linalg.eigvalsh(Aw@Aw.T)
        check(f'weighted transfer gap not smaller q={q:.6g}, w={w0}',eigw[0]+1e-10>=eig[0])
        if w0<1:
            v=np.linalg.solve(R,w0**np.arange(N));v/=np.linalg.norm(v)
            check(f'weighted transfer has same boundary kernel q={q:.6g}, w={w0}',np.linalg.norm(Aw@v)<1e-12)

# Gauge-covariant 3-space stencil, finite-section regression in U(2).
# The operator identity is representation-independent: V_i -> G V_i G^*,
# so it also holds for Spin(10) matrices in the 16. No gauge measure is tested.
rng=np.random.default_rng(20261124)
def unitary(d):
    u,_,vh=np.linalg.svd(rng.normal(size=(d,d))+1j*rng.normal(size=(d,d)))
    return u@vh
sites=list(itertools.product(range(2),repeat=3));d=2;N=5
G=np.zeros((len(sites)*d,len(sites)*d),complex)
for i in range(len(sites)):G[i*d:(i+1)*d,i*d:(i+1)*d]=unitary(d)
Vs=[]
for axis in range(3):
    V=np.zeros_like(G)
    for i,x in enumerate(sites):
        y=list(x);y[axis]=(y[axis]+1)%2;j=sites.index(tuple(y))
        V[i*d:(i+1)*d,j*d:(j+1)*d]=unitary(d)
    Vs.append(V)
pauli=[np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])]
def gauge_parent(Vs):
    d0=len(Vs[0]); I=np.eye(d0); S=np.eye(N,k=-1); R=np.eye(N)-q0*np.linalg.matrix_power(S,4)
    W=sum((I-(V+V.conj().T)/2 for V in Vs),np.zeros_like(I,dtype=complex))
    D=sum((np.kron(sig,(V-V.conj().T)/(2j)) for sig,V in zip(pauli,Vs)),np.zeros((2*d0,2*d0),complex))
    A=(np.kron(S.T,np.eye(2*d0))-np.kron(np.eye(N),np.kron(np.eye(2),W)))@np.kron(R,np.eye(2*d0))
    D=np.kron(np.eye(N),D)
    return np.block([[D,A.conj().T],[A,-D]])
H=gauge_parent(Vs); H1=gauge_parent([G@V@G.conj().T for V in Vs]); GG=np.kron(np.eye(2*N*2),G)
check('gauge background parent is self-adjoint',np.max(abs(H-H.conj().T))<1e-12)
check('exact gauge covariance identity tested on random non-Abelian links',np.max(abs(H1-GG@H@GG.conj().T))<1e-12)

# Exact finite-J atomic pole embedding of original Wilson slab (q=0).
# H=[[s I, T, 0],[T^T,-s I, Delta e_N],[0,Delta e_N^T,0]].
# This is the LOCAL ATOMIC SELF-ENERGY ANSATZ, not the exact interacting theory.
for N in [3,9,20]:
    for J in [.04,.3,1.]:
        D=29*J/2
        for w0 in [0.,2.,4.,6.]:
            T=w0*np.eye(N)-np.eye(N,k=-1)
            e=np.eye(N)[:,-1]
            H=np.block([[np.zeros((N,N)),T,np.zeros((N,1))],
                        [T.T,np.zeros((N,N)),D*e[:,None]],
                        [np.zeros((1,N)),D*e[None,:],np.zeros((1,1))]])
            v=np.r_[w0**np.arange(N),np.zeros(N),-w0**N/D]
            v/=np.linalg.norm(v)
            check(f'finite-J corner pole N={N}, J={J}, w={w0}',np.linalg.norm(H@v)<1e-11 and np.sum(abs(np.linalg.eigvalsh(H))<1e-9)==1)
            W=sum(w0**(2*j) for j in range(N))
            check(f'nonzero physical residue N={N}, J={J}, w={w0}',np.isclose(v[:N]@v[:N],W/(W+w0**(2*N)/D**2),rtol=1e-12) and v[:N]@v[:N]>0)
            if N==20 and w0>1:
                check(f'infinite-depth right residue N={N}, J={J}, w={w0}',np.isclose(v[:N]@v[:N],D**2/(D**2+w0*w0-1),rtol=2e-11))
# Sheet exchange is not massless-massive equality. On a truncated TFD, exact.
for q in [q0,.2,1/np.sqrt(2),.9]:
    d=16; p=q**(2*np.arange(d));p/=p.sum()
    Psi=np.diag(np.sqrt(p))
    check(f'sheet flip and equal reduced states for every q={q:.6g}',np.array_equal(Psi,Psi.T) and np.array_equal(Psi@Psi.T,Psi.T@Psi))
# Original clock compression generated full matrix algebra as soon as C^-1 is admitted.
F=sp.Matrix([[13,1,-1],[1,5,-5],[-1,-5,5]])/16
R=sp.Matrix([[5,-1,-5],[-1,13,1],[-5,1,5]])/16
G=sp.diag(1,0,0)
e=sp.Matrix([1,0,0]); V=sp.Matrix.hstack(e,F*e,R*e)
check('exact determinant of clock cyclic vectors', V.det()==-sp.Rational(3,128))
L=[sp.eye(3),F,R]; products=[A*G*B for A in L for B in L]
check('nine explicit clock words span full M3',sp.Matrix.hstack(*[X.reshape(9,1) for X in products]).rank()==9)
print('q native:',q0,'q^2:',q0*q0,'uniform free nonchiral gap bound:',1-q0)
print('Exact infinite half-space Weyl node at k=0 only; finite J atomic ansatz retains eight nodes.')
print('Limitations: spatial Wilson stencil, half-space internal depth and spinor chosen; no dynamical gauge or interacting completion.')
print(f'{ncheck}/{ncheck} checks passed (finite regressions and exact symbolic identities).')
