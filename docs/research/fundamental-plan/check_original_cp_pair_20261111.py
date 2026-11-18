"""Actual v252/v88 readouts; spectral equivalence versus fixed weak current."""
from pathlib import Path
import sys
import numpy as np
import sympy as sp
root = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0,str(root/'verification'))
import v252_full_finite_triple as f
from tfpt_constants import phi0
from v568_bit_selector_ladder import Qc

s12 = np.sqrt(float(phi0)*(1-float(phi0)))
s23 = float(phi0)/(1+s12)
s13 = s12**3/3
delta = np.pi/3+3*s12*s12
c12,c23,c13=np.sqrt(1-np.array([s12,s23,s13])**2)
z=np.exp(1j*delta)
V=np.array([[c12*c13,s12*c13,s13/z],
[-s12*c23-c12*s23*s13*z,c12*c23-s12*s23*s13*z,s23*c13],
[s12*s23-c12*c23*s13*z,-c12*s23-s12*c23*s13*z,c23*c13]])
m_u=np.diag([.0022,1.27,172.76]); m_d=np.diag([.0047,.095,4.18])
m_e=np.diag([.000511,.1057,1.777])
zero=np.zeros((3,3))
Dp=f.full_D(m_u,V@m_d,m_e,zero,zero)
Dm=f.full_D(m_u,V.conj()@m_d,m_e,zero,zero)
assert abs(np.linalg.det(V)-1)<1e-14
assert abs(np.angle(np.linalg.det(m_u)*np.linalg.det(V@m_d))) < 1e-14
assert abs(np.angle(np.linalg.det(m_u)*np.linalg.det(V.conj()@m_d))) < 1e-14

def indices(name,ch,c):
    return [f.IDX[(name,ch,c,g)] for g in range(3)]
W=np.eye(96,dtype=complex); G=np.zeros((96,96),complex)
for c in range(3):
    d,u=indices('d','L',c),indices('u','L',c)
    W[np.ix_(d,d)]=V.conj()@V.conj().T
    W[np.ix_([48+i for i in d],[48+i for i in d])]=W[np.ix_(d,d)].conj()
    G[np.ix_(u,d)]=np.eye(3)  # actual particle weak raising matrix from v252 rep
p_a=f.rep(0,np.array([[0,1],[0,0]]),np.zeros((3,3)))
quarkL=[i for c in range(3) for label in ('u','d') for i in indices(label,'L',c)]
assert np.array_equal(G[np.ix_(quarkL,quarkL)],p_a[np.ix_(quarkL,quarkL)])
assert np.linalg.norm(W@Dp@W.conj().T-Dm)<1e-12
assert np.linalg.norm(W@G-G@W)>1e-3
assert np.linalg.norm(W.conj().T@W-np.eye(96))<1e-12
assert np.max(abs(np.linalg.eigvalsh(Dp)-np.linalg.eigvalsh(Dm)))<1e-12
# Exhaust the finite real *-algebra A_F on a real basis, for the unchanged v252 map.
z3=np.zeros((3,3)); z2=np.zeros((2,2))
A=[f.rep(l,z2,z3) for l in (1,1j)]
A.extend(f.rep(0,f.quaternion(a,b),z3) for a,b in ((1,0),(1j,0),(0,1),(0,1j)))
for j,k in np.ndindex(3,3):
    for factor in (1,1j):
        t=z3.astype(complex); t[j,k]=factor; A.append(f.rep(0,z2,t))
K=f.Kswap(); gamma=f.grading()
orderzero=max(np.max(abs(a@b-b@a)) for a in A for b in (K@x.conj()@K for x in A))
assert orderzero<1e-13
for D in (Dp,Dm):
    assert np.max(abs(D-D.conj().T)) < 1e-13
    assert np.max(abs(D@gamma+gamma@D)) < 1e-13
    assert np.max(abs(K@D.conj()-D@K)) < 1e-13
    assert max(np.max(abs((D@a-a@D)@b-b@(D@a-a@D)))
               for a in A for b in (K@x.conj()@K for x in A)) < 1e-12
# The neutral gauge generators (diagonal weak, colour and complex summands)
# commute with W and its transported antiparticle action.
N=[A[0],A[1],A[2],A[3]]+A[6:]
assert max(np.max(abs(W@a-a@W)) for a in N)<1e-12
assert max(np.max(abs(W@(K@a.conj()@K)-(K@a.conj()@K)@W)) for a in N)<1e-12

def rank1(v):
    return np.outer(v,v.conj())
def loop(D):
    u,d=indices('u','L',0),indices('d','L',0)
    ud=G[np.ix_(u,d)]
    D2=D@D
    _,pu=np.linalg.eigh(D2[np.ix_(u,u)])
    _,pd=np.linalg.eigh(D2[np.ix_(d,d)])
    return np.trace(rank1(pu[:,0])@ud@rank1(pd[:,0])@ud.conj().T
                    @rank1(pu[:,1])@ud@rank1(pd[:,1])@ud.conj().T)
j1,j2=loop(Dp),loop(Dm)
j_expected=s12*c12*s23*c23*s13*c13**2*np.sin(delta)
assert abs(j1.imag-j_expected)<1e-12 and abs(j2-j1.conjugate())<1e-12

# Original covariance reader v258/v724 (same dimensionless KMS time).
# Scale is only for a well-conditioned finite readout; the identity is exact
# for every common positive scale beta, not a fit or physical time conversion.
beta=1/np.linalg.norm(Dp,2)
def cov(D):
    d,u=np.linalg.eigh(D)
    return (u*(1/(1+np.exp(beta*d))))@u.conj().T
C1,C2=cov(Dp),cov(Dm)
assert np.linalg.norm(W@C1@W.conj().T-C2)<1e-12
# The original helper is real-symmetric only; its spectral formula extends
# with adjoint to these complex Hermitian full96 covariances.
def logit(C):
    d,u=np.linalg.eigh(C)
    return (u*np.log((1-d)/d))@u.conj().T
assert np.linalg.norm(logit(C1)/beta-Dp)<2e-11
assert np.linalg.norm(logit(C2)/beta-Dm)<2e-11
assert np.linalg.norm(C1-C2)>1e-5

# Same FULL reduced state on the left chirality at equal time. At half
# thermal time the compressed transfer resolves the mass-squared projectors.
left = np.flatnonzero(np.diag(gamma)>0)
assert np.linalg.norm(C1[np.ix_(left,left)]-np.eye(len(left))/2)<1e-12
assert np.linalg.norm(C2[np.ix_(left,left)]-np.eye(len(left))/2)<1e-12
def midpoint(D):
    d,u=np.linalg.eigh(D)
    return (u*(0.5/np.cosh(beta*d/2)))@u.conj().T
def heat_proj(F, idx, s):
    v,U=np.linalg.eigh(F[np.ix_(idx,idx)])
    # Midpoint eigenvalues are in reverse mass order.
    return rank1(U[:,-1-s])
def midpoint_loop(F):
    u,d=indices('u','L',0),indices('d','L',0)
    ud=G[np.ix_(u,d)]
    return np.trace(heat_proj(F,u,0)@ud@heat_proj(F,d,0)@ud.conj().T
                    @heat_proj(F,u,1)@ud@heat_proj(F,d,1)@ud.conj().T)
F1,F2=midpoint(Dp),midpoint(Dm)
assert abs(midpoint_loop(F1)-j1)<1e-7
assert abs(midpoint_loop(F2)-j2)<1e-7
print('same complete reduced left equal-time state (C_L = I/2)')
print('half-time left-state reconstruction of marked loop',
      midpoint_loop(F1).imag,midpoint_loop(F2).imag)

# Separate candidate: the already existing Q twin fails the original fixed
# P2 compiler readout, before we ask for a dynamical source.
R=sp.Matrix([[1,3,0],[1,5,2],[2,5,3]])
Sig=sp.diag(1,-1,-1)
one=sp.Matrix([1,1,1]); a=sp.Matrix([1,1,2])
out=[]
for c in (2,0):
    Q=Qc(c); M=R+Q*Sig
    rows=sp.Matrix.vstack(one.T*M,a.T*M)
    minors=[rows[:,[i,j]].det() for i,j in [(0,1),(0,2),(1,2)]]
    out.append(5*sum(map(abs,minors))/(9*rows[0,0]))
    print('Compiler twin',c,'columns',list(one.T*Q),'quark readout',out[-1])
assert out[0]==sp.Rational(55,117) and out[0]!=out[1]
print('all 24 x 24 real basis commutators pass in both original full96 triples')
print('common spectrum, neutral conjugation, phase of determinant zero')
print('marked loop',repr(j1),repr(j2))
print('fixed-frame covariance difference (HS):',np.linalg.norm(C1-C2))
print('weak-current commutator with spectral intertwiner (HS):',np.linalg.norm(W@G-G@W))
