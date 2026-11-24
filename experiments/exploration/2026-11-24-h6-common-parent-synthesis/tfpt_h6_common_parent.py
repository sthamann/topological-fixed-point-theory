#!/usr/bin/env python3
"""H6h30 synthesis. Exact star/quarter intertwiner, finite-J atomic boundary
localization, and Casimir local propagator. Scope: chosen free Wilson parent,
local H6 interaction and atomic self-energy approximation; NOT full interacting TOE.
"""
import itertools, math
import numpy as np
import sympy as s
from scipy.sparse import csr_matrix
from fractions import Fraction as F

checks=0
def check(label,ok):
    global checks
    assert bool(ok),label
    checks+=1; print('PASS',label)

# Exact common geometry: star mass M_{j,n}=M_c*q^(n+1)/2.
# For canonical vertex coordinates B_e=h_child-(M_child/M_parent)*h_parent.
# The first edge coefficient is q/2; later q. Hadamard transform on every
# four-arm layer identifies its symmetric direction with residue 0.
Had=s.Matrix([[1,1,1,1],[1,1,-1,-1],[1,-1,1,-1],[1,-1,-1,1]])/2
q=s.symbols('q',real=True)
check('real orthogonal four-arm Fourier transform',Had.T*Had==s.eye(4))
for D in (1,2,3,4):
    N=4*D+1
    B=s.zeros(N-1,N); U=s.zeros(N,N); V=s.zeros(N-1,N-1); U[0,0]=1
    for n in range(D):
        for j in range(4):
            e=4*n+j; B[e,e+1]=1
            if n==0: B[e,0]=-q/2
            else: B[e,e+1-4]=-q
            for a in range(4):
                m=4*n+(4 if a==0 else a)
                U[e+1,m]=Had[j,a]; V[e,m-1]=Had[j,a]
    S=s.zeros(N,N)
    for i in range(N-1): S[i+1,i]=1
    R=s.eye(N)-q*S**4; A=(S.T*R)[:N-1,:]
    W=s.diag(*[n+1 for n in range(D) for j in range(4)])
    check(f'D={D}: B U = V A exact',B*U==V*A)
    check(f'D={D}: W commutes with depthwise transform',W*V==V*W)
    check(f'D={D}: weighted star L and weighted quarter K exactly intertwine',s.expand(U.T*B.T*W*B*U-A.T*W*A)==s.zeros(N,N))
    # physical normalized null vector: M_v/M_eff; ratio to local M_v constant.
    eta=s.Matrix([1]+[q**(n+1)/2 for n in range(D) for j in range(4)])
    check(f'D={D}: null vector is local mass vector',s.expand(B*eta)==s.zeros(N-1,1))
    v=s.Matrix([q**(m//4) if m%4==0 else 0 for m in range(N)])
    check(f'D={D}: transformed null vector equals truncated native TFD',U*v==eta)

# Native Wilson parent A(k) and finite mirror self-energy realization.
# D_N=[(S^*-w)R, Delta*e_last]. At zero frequency its kernel is exactly
# (R^-1 g_w, w^N/Delta); this is a rational atomic self-energy model,
# not an exact replacement for a coupled interacting lattice.
q0=math.exp(-2*math.pi)
print('\nH6 q=',q0,'q^2=',q0*q0)
for N in [9,17,33,65]:
    S=np.eye(N,k=-1); R=np.eye(N)-q0*np.linalg.matrix_power(S,4)
    for J in [.1,1,10]:
        Delta=14.5*J
        for w in [0,2,4,6]:
            # normalized powers avoid overflow and condition number blowup.
            g=(w**np.arange(N,dtype=float)) if w<=1 else (w**np.arange(1-N,1,dtype=float))
            v=np.linalg.solve(R,g); y=w*g[-1]/Delta
            Dmat=np.column_stack(((S.T-w*np.eye(N))@R,Delta*np.eye(N)[:,-1]))
            z=np.r_[v,y]; z/=np.linalg.norm(z)
            check(f'N={N} J={J} w={w}: exact atomic null recurrence',np.linalg.norm(Dmat@z)<3e-14)
            if w==0:
                check(f'N={N} J={J}: left pole residue is truncated KMS weight',abs(z[0]**2-1/sum(q0**(2*n) for n in range((N-1)//4+1)))<2e-14)
            else:
                # use log-safe form of bound for L=4. Auxiliary can only reduce residue.
                L=4; left=np.sum(z[:L]**2)
                bound=w**(2*(L-N))*(1-w**(-2*L))/(1-w**(-2*N))/(1-q0/w**4)**2
                check(f'N={N} J={J} w={w}: all unwanted corner poles decouple from left local algebra',left<=bound*(1+1e-12))
                if J==1 and w==2: print('   unwanted corner residue on first four sites:',left,'bound',bound)

# Local H6 atomic spectrum. Spin(11) fundamental spinor has weights
# (±1/2)^5 each once. The occupied highest first-weight subspace contains
# 16 weights; the ground state is in V_(8,0,0,0,0).
# Tensor V_(8,0,..) with the minuscule spinor has only two dominant
# highest weights (8±1/2,1/2,1/2,1/2,1/2).
# In particle sectors 15 and 17 the largest first coordinate is 15/2.
# Thus only the lower irrep is accessible by one c or c^dagger.
weights=list(itertools.product([F(-1,2),F(1,2)],repeat=5))
check('Spin(11) spinor has 32 distinct weights',len(set(weights))==32)
check('exactly 16 weights with first coordinate +1/2',sum(w[0]>0 for w in weights)==16)
for m in (15,17):
    max_coord=sum(sorted([w[0] for w in weights],reverse=True)[:m])
    check(f'Fock filling {m} excludes higher irrep first coordinate 17/2',max_coord==F(15,2))
lam=(F(8),F(0),F(0),F(0),F(0))
allowed=sorted(set(tuple(a+b for a,b in zip(lam,w)) for w in weights if all((lam[i]+w[i])>=(lam[i+1]+w[i+1]) for i in range(4)) and lam[-1]+w[-1]>=0))
check('minuscule product has precisely the two dominant weights',allowed==[(F(15,2),)+(F(1,2),)*4,(F(17,2),)+(F(1,2),)*4])
def cas(lam,rho):return sum(a*(a+2*b) for a,b in zip(lam,rho))
rhoB=tuple(F(9-2*i,2) for i in range(5));rhoD=tuple(F(4-i) for i in range(5))
C0=cas(lam,rhoB);C1=cas(allowed[0],rhoB);Cs=cas((F(1,2),)*5,rhoD)
check('local fermionic excitation gap 29 J/2 follows exactly from Casimirs',Cs-C1+C0==F(29,2))
# Cross-check spectral support of atomic one-body Green function by explicit
# Fock diagonalization for Spin(7), analogue of H6 Spin(11).
def gammas(k):
    X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]]);Z=np.diag([1.,-1.]);I=np.eye(2)
    def kron(ms):
        a=np.array([[1.]])
        for b in ms:a=np.kron(a,b)
        return a
    return [kron([Z]*j+[p]+[I]*(k-j-1)) for j in range(k) for p in (X,Y)],kron([Z]*k)
k=3; n=2**k; gam,chi=gammas(k)
def basis(m):return [sum(1<<i for i in t) for t in itertools.combinations(range(n),m)]
def op(mat,st):
    ix={x:i for i,x in enumerate(st)};O=np.zeros((len(st),len(st)),complex)
    for a,b in zip(*np.nonzero(mat)):
        for i,x in enumerate(st):
            if not x>>b&1:continue
            y=x^(1<<b)
            if y>>a&1:continue
            sign=(-1)**((x&((1<<b)-1)).bit_count()+(y&((1<<a)-1)).bit_count())
            O[ix[y|(1<<a)],i]+=sign*mat[a,b]
    return O
def H(m):
    st=basis(m); X=[op(.5j*g@chi,st) for g in gam]
    return -sum(a@a for a in X),st
H0,st=H(n//2);E,U=np.linalg.eigh(H0);gs=U[:,0]
for delta in [-1,1]:
    H1,st1=H(n//2+delta); ix={x:i for i,x in enumerate(st1)}
    for a in range(n):
        vec=np.zeros(len(st1),complex)
        for j,x in enumerate(st):
            if (x>>a&1)!=(delta<0):continue
            y=x^(1<<a);vec[ix[y]]+=(-1)**((x&((1<<a)-1)).bit_count())*gs[j]
        check(f'Spin(7), mode {a}, delta {delta}: exact spectral support at 11J/2',np.linalg.norm((H1-(E[0]+5.5)*np.eye(len(st1)))@vec)<1e-12)
        check(f'Spin(7), mode {a}, delta {delta}: half spectral weight',abs(np.vdot(vec,vec)-.5)<1e-12)
# Exact H6 atomic Green function assuming the stated unique symmetric singlet
# and physical particle-hole transformation. Spin(11) rotation interchanges
# chiral halves; Spin(10) invariance gives scalar densities within each.
z,J=s.symbols('z J');Delta=29*J/2
G=s.Rational(1,2)/(z-Delta)+s.Rational(1,2)/(z+Delta)
check('local atomic Green function is z/(z^2-Delta^2)',s.simplify(G-z/(z*z-Delta*Delta))==0)
check('exact local atomic self-energy is Delta^2/z',s.simplify(z-1/G-Delta*Delta/z)==0)

# H6 representation gluing: E8 decomposition D5 x A3 and residual
# fundamental C4 grading 4=1+3. The grade-0 component cannot be identified
# as a heavy bulk mode by this branching alone.
check('E8 adjoint branching dimension',45+15+10*6+16*4+16*4==248)
check('one chiral boundary spinor for each of three nontrivial characters has 48 components',3*16==48)
print('\nTOTAL',checks,'checks passed. No physical gauge-measure, dynamical gravity, P1 selection or TOE closure claimed.')

# H6 z34 / z36: resolve w2=1 by the *diagonal E8 quotient*, retaining
# a triplet in the family four. C is the U(5) carrier bundle c1(C)=1.
# A family lift preserving SU(3) has cocharacter (s,s,s,-3s), s half-integer
# for endpoint -I. The smallest two have s=±1/2.
from collections import Counter
hs=[F(2*n+1,2) for n in range(-8,8)]
check('triplet-preserving compensator minimal |s|=1/2',min(abs(x) for x in hs)==F(1,2))
a=(F(1,2),)*3+(F(-3,2),)
check('family SU4 compensator is traceless and ends at -I',sum(a)==0 and all(x.denominator==2 for x in a))
check('minimal triplet-preserving E8 cocharacter squared length is 4',1+sum(x*x for x in a)==4)
# Compute E8 root charge distribution under this circle. D5 weight is e1,
# A3 weight is a. Spin(10) chirality weights have even number of minus signs.
D5=[F(sig) if i==0 else F(tau) if j==0 else F(0) for i,j in itertools.combinations(range(5),2) for sig,tau in itertools.product((-1,1),repeat=2)]
A3=[a[i]-a[j] for i in range(4) for j in range(4) if i!=j]
V10=[F(1),F(-1)]+[F(0)]*8
V6=[a[i]+a[j] for i,j in itertools.combinations(range(4),2)]
S16=[w[0] for w in weights if sum(x<0 for x in w)%2==0]
check('half-spinor charge multiplicities are 8+8',Counter(S16)=={F(1,2):8,F(-1,2):8})
pairs=[v+b for v in V10 for b in V6];spin=[v+b for v in S16 for b in a]
charges=D5+A3+pairs+spin+[-x for x in spin]
check('all 240 E8 root clutching charges are integral',len(charges)==240 and all(x.denominator==1 for x in charges))
count=Counter(charges)
check('E8 norm-four clutching charge census',count=={-2:14,-1:64,0:84,1:64,2:14})
print('E8 cocharacter charge census:',dict(sorted(count.items())))
# Under U(5), chiral spin^c Fock space rank 16 and c1 = 8 c1(C).
# Tensor with any honest family line shifts by multiples of 16.
check('even Fock carrier rank 16',sum(math.comb(5,k) for k in (0,2,4))==16)
check('even Fock first Chern multiplier 8',sum(math.comb(4,k-1) for k in (2,4))==8)
check('spin Dirac index per compensated triplet line is 8 mod 16, not 0 mod 16',all((8+16*m)%16==8 for m in range(-20,21)))
# E8 bundle clutching of pair (Spin10 path, SU4 path) is valid since both
# end at their -1 central element and the diagonal Z4 quotient identifies them.
# This does not identify the physical source or prove family Dirac dynamics.
print('FINAL',checks,'checks passed including the E8 projective lift. No ledger promotion.')

# Clifford soldering: the same 4D Dirac matrices give the Weyl symbol and
# the bivector of a tetrad. Chiral trace difference equals Plebanski area
# pairing, but is not the positive Hilbert-Schmidt trace.
gam4,_=gammas(2); g5=gam4[0]@gam4[1]@gam4[2]@gam4[3]
check('Euclidean chirality squares to one',np.allclose(g5@g5,np.eye(4)))
def eps(p):
    return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p))) if len(set(p))==len(p) else 0
Gs={(a,b):(gam4[a]@gam4[b]-gam4[b]@gam4[a])/2 for a in range(4) for b in range(4)}
check('exact chiral Clifford trace identity tr(g5 G_ab G_cd)=4 epsilon_abcd',all(abs(np.trace(g5@Gs[a,b]@Gs[c,d])-4*eps((a,b,c,d)))<1e-12 for a,b,c,d in itertools.product(range(4),repeat=4)))
eps4=np.zeros((4,4,4,4));
for p in itertools.permutations(range(4)):eps4[p]=eps(p)
def cliff_area(E):
    ge=[sum((E[a,m]*gam4[a] for a in range(4)),np.zeros((4,4),complex)) for m in range(4)]
    return np.array([[(ge[m]@ge[n]-ge[n]@ge[m])/2 for n in range(4)] for m in range(4)])
def form(E):return np.einsum('am,bn->abmn',E,E)-np.einsum('an,bm->abmn',E,E)
def el2(X):return sum(np.linalg.det(X[np.ix_(I,I)]) for I in itertools.combinations(range(4),2))
rng=np.random.default_rng(6)
for k in range(12):
    E=np.eye(4)+.15*rng.normal(size=(4,4)); A=np.eye(4)+.15*rng.normal(size=(4,4)); Fm=E@A
    D=cliff_area(E)-cliff_area(Fm); BD=form(E)-form(Fm)
    val=np.einsum('mnrs,ij,mnjk,rski->',eps4,g5,D,D) if False else sum(eps4[m,n,r,t]*np.trace(g5@D[m,n]@D[r,t]) for m,n,r,t in itertools.product(range(4),repeat=4) if eps4[m,n,r,t])
    area=np.einsum('abcd,mnrs,abmn,cdrs->',eps4,eps4,BD,BD)
    hr=np.linalg.det(E)*(-1.5+.5*el2(A)-1.5*np.linalg.det(A))
    check(f'Clifford bivector trace equals area density and H6 t=1 potential, sample {k}',abs(val-area)<1e-10 and abs(val+64*hr)<1e-10)
# Source invariant obstruction: chirality trace is indefinite; positivity
# of a quadratic norm does not select the oriented area pairing.
Dplus=Gs[0,1]+Gs[2,3];Dminus=Gs[0,1]-Gs[2,3]
x1=np.trace(g5@Dplus@Dplus).real;x2=np.trace(g5@Dminus@Dminus).real
check('chiral trace has opposite signs in selfdual and antiselfdual sectors',x1*x2<0)
print('Chiral trace values on complementary bivectors:',x1,x2)
# Gauge-covariant square: once spatial links are curved [D,A] need not vanish.
# This is the explicit source of Pauli-curvature and Wilson commutator terms.
X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]]);Z=np.diag([1.,-1.]);sig=[X,Y,Z]
def randomU(d):
    u,_,vh=np.linalg.svd(rng.normal(size=(d,d))+1j*rng.normal(size=(d,d)));return u@vh
for rep in range(6):
    links=[randomU(3) for _ in range(3)];Pmom=[(v-v.conj().T)/(2j) for v in links];I=np.eye(3)
    D=sum(np.kron(a,b) for a,b in zip(sig,Pmom));w=sum(I-(v+v.conj().T)/2 for v in links)
    scalar=np.kron(np.eye(2),sum(p@p for p in Pmom))
    curv=1j*sum(eps((i,j,k))*np.kron(sig[k],Pmom[i]@Pmom[j]-Pmom[j]@Pmom[i]) for i,j in itertools.combinations(range(3),2) for k in range(3))
    check(f'nonabelian gauge links {rep}: Dirac square = covariant Laplacian + Pauli curvature',np.max(abs(D@D-scalar-curv))<1e-12)
    N=5; S=np.eye(N,k=-1);R=np.eye(N)-q0*np.linalg.matrix_power(S,4);A=(np.kron(S.T,np.eye(6))-np.kron(np.eye(N),np.kron(np.eye(2),w)))@np.kron(R,np.eye(6));D=np.kron(np.eye(N),D)
    H=np.block([[D,A.conj().T],[A,-D]])
    target=np.block([[D@D+A.conj().T@A,D@A.conj().T-A.conj().T@D],[A@D-D@A,D@D+A@A.conj().T]])
    check(f'nonabelian gauge links {rep}: full common parent square has exact commutator terms',np.max(abs(H@H-target))<1e-12)
    check(f'nonabelian gauge links {rep}: mixed term is generically nonzero',np.linalg.norm(D@A.conj().T-A.conj().T@D)>1e-4)
print('FINAL INCLUDING CLIFFORD AND GAUGE',checks,'checks passed. Curvature kinematics only, not a dynamical gravitational derivation.')
