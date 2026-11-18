#!/usr/bin/env python3
"""Reproducible finite audit: TFPT source response, Petz binding and block recursion.

New coupling is defined on the 35-dimensional logical response spaces. This does
not identify it with an already selected microscopic TFPT field Hamiltonian.
Exact rational identities and the explicitly numerical higher-order probe are
reported separately. No network, no hidden imported result matrices, no TOE claim.
"""
from __future__ import annotations
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from itertools import product
from math import prod, factorial
from pathlib import Path
import json, time, hashlib, platform
import numpy as np
from scipy import sparse as sp
import sympy as sy

ROOT=Path(__file__).resolve().parent
T0=time.time();checks=[];results={}
def check(name, condition, kind='exact', **detail):
    if not bool(condition): raise RuntimeError(f'{name} failed: {detail}')
    checks.append(dict(name=name,kind=kind,passed=True,**detail))
def zero(M): return all(sy.simplify(x)==0 for x in M)
def clean(x):
    if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [clean(v) for v in x]
    if isinstance(x,np.ndarray):return clean(x.tolist())
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,sy.MatrixBase):return clean(x.tolist())
    if isinstance(x,sy.Basic):return str(x)
    return x

# Reconstruct the symmetric fourth power in unnormalized permutation-orbit kets.
I=np.eye(2,dtype=complex);X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex);Z=np.diag([1,-1]).astype(complex)
As=[np.kron(a,b) for a in (I,X,Y,Z) for b in (I,X,Y,Z)]
labels=[a+b for a in 'IXYZ' for b in 'IXYZ']
ns=sorted([n for n in product(range(5),repeat=4) if sum(n)==4],reverse=True)
pos={n:i for i,n in enumerate(ns)};d=35
m=np.array([24//prod(factorial(v) for v in n) for n in ns],np.int64)
Ks=[];Ds=[]
for A in As:
    K=np.zeros((d,d),complex);D=np.zeros((d,d),complex)
    for c,n in enumerate(ns):
        for j in range(4):
            if not n[j]:continue
            for i in range(4):
                if A[i,j]==0:continue
                if i==j:K[c,c]+=n[j]*A[i,j]
                else:
                    nn=list(n);nn[j]-=1;nn[i]+=1
                    K[pos[tuple(nn)],c]+=(n[i]+1)*A[i,j]
        nn=[0]*4;ph=1
        for j in range(4):
            i=int(np.flatnonzero(A[:,j])[0]);nn[i]=n[j];ph*=A[i,j]**n[j]
        D[pos[tuple(nn)],c]=ph
    check('integer_Gaussian_collective_'+str(len(Ks)),
          np.array_equal(K.real,np.rint(K.real)) and np.array_equal(K.imag,np.rint(K.imag)))
    check('metric_selfadjoint_collective_'+str(len(Ks)),np.array_equal(m[:,None]*K,K.conj().T*m[None,:]))
    Ks.append(K);Ds.append(D)
Pnum=np.sum(Ds,axis=0).real.astype(np.int64) # 16 Pi5
E=np.zeros((35,5),np.int64)
supports=[[(4,0,0,0),(0,4,0,0),(0,0,4,0),(0,0,0,4)],[(2,2,0,0),(0,0,2,2)],
          [(2,0,2,0),(0,2,0,2)],[(2,0,0,2),(0,2,2,0)],[(1,1,1,1)]]
for j,S in enumerate(supports):
    for n in S:E[pos[n],j]=1
G=E.T@np.diag(m)@E
Gs=sy.diag(*[int(x) for x in np.diag(G)]);Gi=Gs.inv()
check('code_Gram',Gs==sy.diag(4,12,12,12,24))
check('code_projector',np.array_equal(Pnum@Pnum,16*Pnum) and np.trace(Pnum)==80)
check('code_range',np.array_equal(Pnum@E,16*E))
R=[]
for a,K in enumerate(Ks[1:],1):
    TK=E.T@np.diag(m)@K@K@E
    check('return_real_'+str(a),np.count_nonzero(TK.imag)==0)
    Ra=Gi*sy.Matrix(TK.real.astype(np.int64))/16
    check('plane_projector_'+str(a),Ra*Ra==Ra and sy.trace(Ra)==2)
    check('one_source_invisible_'+str(a),np.count_nonzero(E.T@np.diag(m)@K@E)==0)
    # Both endpoints of a first response lie in the matching Pauli-character sector.
    P_a=sum([(1 if np.array_equal(As[a]@b,b@As[a]) else -1)*Db for b,Db in zip(As,Ds)])
    check('response_sector_'+str(a),np.array_equal(P_a@K@E,16*K@E))
    R.append(Ra)
check('plane_sum',sum(R,sy.zeros(5))==6*sy.eye(5))

# Reconstruct the Bell effects from the explicitly supplied code-basis contractions.
Fnum=sy.Matrix([[4,4,4,4,0],[0,8,0,0,8],[4,-4,4,-4,0],[0,0,8,0,8],
               [0,0,0,8,8],[0,0,8,0,-8],[0,0,0,8,-8],[4,4,-4,-4,0],
               [0,8,0,0,-8],[4,-4,-4,4,0]])
Recovery=sy.zeros(25)
for r in range(10):
    f=Fnum[r,:].T/4;v=Gi*f
    Recovery+=2*sy.kronecker_product(v,v)*sy.kronecker_product(f,f).T
Qpair=sum([sy.kronecker_product(Ra,Ra) for Ra in R],sy.zeros(25))
check('Petz_binding_operator_identity',Qpair==sy.Rational(3,2)*sy.eye(25)+sy.Rational(9,2)*Recovery)
check('recovery_spectrum',Recovery.eigenvals()=={sy.Integer(1):1,sy.Rational(4,9):9,sy.Integer(0):15})
check('plane_pair_spectrum',Qpair.eigenvals()=={sy.Integer(6):1,sy.Rational(7,2):9,sy.Rational(3,2):15})
results['pair_identity']={'R_rank':2,'R_count':15,'sum_R':'6 I',
 'sum_R_tensor_R':'(3/2) I25 + (9/2) vectorized(Petz o E2)',
 'pair_h':'I25 - vectorized(Petz o E2)',
 'pair_h_spectrum':{'0':1,'5/9':9,'1':15},
 'second_order_H':'-192 epsilon^2/delta I25 -576 epsilon^2/delta vectorized(Petz o E2)'}
print('response and binding identity exact',flush=True)

# Exact complete 35 x 35 two-module model. L is selfadjoint in metric m tensor m.
L=sum([sp.kron(sp.csr_matrix(K),sp.csr_matrix(K),format='csr') for K in Ks[1:]])
check('bridge_real',np.count_nonzero(L.data.imag)==0)
L=L.real.astype(np.int64)
Hnum=32*sp.eye(d*d,format='csr',dtype=np.int64)-sp.kron(sp.csr_matrix(Pnum),sp.eye(d,format='csr',dtype=np.int64))-sp.kron(sp.eye(d,format='csr',dtype=np.int64),sp.csr_matrix(Pnum))
metric=np.kron(m,m);Tcode=np.kron(E,E)
B=L@Tcode
check('bridge_first_order_zero',np.count_nonzero(Tcode.T@(metric[:,None]*B))==0)
check('all_first_responses_energy_2delta',np.array_equal(Hnum@B,32*B))
Gt=sy.kronecker_product(Gs,Gs)
Bgram=Gt.inv()*sy.Matrix(B.T@(metric[:,None]*B))
check('exact_second_order_gram',Bgram==256*Qpair)
w0=sum([(24//G[j,j])*np.kron(E[:,j],E[:,j]) for j in range(5)])
w1=L@w0
ip=lambda u,v:int(np.dot(u*metric,v))
check('singlet_cross_zero',ip(w0,w1)==0)
check('singlet_bridge_norm',ip(w1,w1)==1536*ip(w0,w0))
check('closed_two_state_block',np.array_equal(L@w1,1536*w0+16*w1))
check('singlet_H0_closed',np.array_equal(Hnum@w0,np.zeros(1225,np.int64)) and np.array_equal(Hnum@w1,32*w1))
M1=w1.reshape(35,35)
# Reduced density of the normalized excited invariant state is (I-Pi5)/30.
check('excited_invariant_flat_30',np.array_equal(480*(M1@np.diag(m)@M1.T@np.diag(m)),ip(w1,w1)*(16*np.eye(35,dtype=np.int64)-Pnum)))
t,delta=sy.symbols('t delta',real=True)
Hsing=sy.Matrix([[0,16*sy.sqrt(6)*t],[16*sy.sqrt(6)*t,2*delta+16*t]])
lam=sy.symbols('lam')
check('singlet_characteristic',sy.expand(Hsing.charpoly(lam).as_expr()-(lam**2-(2*delta+16*t)*lam-1536*t**2))==0)
# Proof constants for the conservative global ground-state range |epsilon|/delta <=1/200.
# ||L||<=80 follows directly from 4 sum_16 Swap -16 I. In the neutral complement H0=2delta.
check('ground_gap_bound_constants',sy.Rational(460800,631)-560>170)
results['exact_pair']={'space_dimension':1225,'neutral_dimension':85,
 'orthonormal_singlet_block':Hsing,'ground_energy':'delta+8epsilon-sqrt((delta+8epsilon)^2+1536epsilon^2)',
 'proved_sufficient_ground_range':'0 < abs(epsilon)/delta <= 1/200',
 'proved_gap_lower_bound':'170 epsilon^2/delta',
 'excited_weight':'(1-(delta+8epsilon)/sqrt((delta+8epsilon)^2+1536epsilon^2))/2',
 'reduced_state':'(1-eta) Pi5/5 + eta (I35-Pi5)/30',
 'scope':'exact stated 35-response-space model, not a diagonalization of the full microscopic register parent'}
print('exact two-state ground sector and conservative global bound',flush=True)

# Exact three-cell chain of the derived effective pair operator.
h=sy.eye(25)-Recovery
H3=sy.kronecker_product(h,sy.eye(5))+sy.kronecker_product(sy.eye(5),h)
G3=sy.kronecker_product(Gs,Gs,Gs)
T1=sy.zeros(125,5);T2=sy.zeros(125,5);T3=sy.zeros(125,5)
for a,b in product(range(5),repeat=2):
    T1[25*a+5*a+b,b]=1/Gs[a,a]
    T2[25*a+5*b+a,b]=1/Gs[a,a]
    T3[25*b+5*a+a,b]=1/Gs[a,a]
W6=sy.Matrix([[2,2,-1,-1,-1,-1],[0,0,-1,-1,1,1],[0,0,-1,1,-1,1],[0,0,1,-1,-1,1],[1,-1,0,0,0,0]])
T4=sy.zeros(125,5)
for q in range(6):
    u=W6[:,q]
    T4+=sy.kronecker_product(u,u,u)*(u.T*Gs)/1600
Ts=[T1,T2,T3,T4]
GG=sy.Matrix(4,4,lambda i,j:sy.trace(Gi*Ts[i].T*G3*Ts[j])/5)
HH=sy.Matrix(4,4,lambda i,j:sy.trace(Gi*Ts[i].T*G3*H3*Ts[j])/5)
Hmult=GG.inv()*HH
for j in range(4):
    check('invariant_five_multiplicity_'+str(j),H3*Ts[j]==sum([Ts[i]*Hmult[i,j] for i in range(4)],sy.zeros(125,5)))
cp=H3.charpoly(lam).as_expr()
expected=(lam-2)**35*(lam-1)**5*(lam-sy.Rational(5,3))**10*(lam-sy.Rational(4,3))**16*(lam-sy.Rational(17,9))**9*(lam-sy.Rational(16,9))**16*(lam-sy.Rational(13,9))**10*(lam-sy.Rational(11,9))**9*(lam-sy.Rational(10,9))**5*(lam**2-sy.Rational(19,9)*lam+sy.Rational(8,9))**5
check('complete_three_chain_characteristic',sy.Poly(cp-expected,lam).is_zero)
s=sy.sqrt(73);e0=(19-s)/18;gap=(s-1)/18
c=sy.Matrix([-(27+3*s)/50,-sy.Rational(12,25),-(27+3*s)/50,1])
check('ground_multiplicity_vector',zero((Hmult-e0*sy.eye(4))*c))
NN=sy.simplify((c.T*GG*c)[0])
Vrg=sum([c[i]*Ts[i] for i in range(4)],sy.zeros(125,5))
check('ground_encoder_isometry_metric',zero(Vrg.T*G3*Vrg-NN*Gs))
alpha=(49+5*s)/144;beta=sy.Rational(2,5)*(1-alpha)
alpha_middle=sy.Rational(17,72)+85*s/5256
for site in (0,1,2):
    al=alpha_middle if site==1 else alpha;be=sy.Rational(2,5)*(1-al)
    for a,Ra in enumerate(R):
        Op=sy.kronecker_product(*[Ra if i==site else sy.eye(5) for i in range(3)])
        Aop=Gi*Vrg.T*G3*Op*Vrg/NN
        check(f'boundary_intertwiner_{site}_{a}',zero(Aop-al*Ra-be*sy.eye(5)))
results['recursive_projection']={'ground_energy':e0,'gap':gap,'ground_dimension':5,
 'four_intertwiner_gram':GG,'four_multiplicity_H':Hmult,'ground_coefficients':c,'ground_norm':NN,
 'endpoint_alpha':alpha,'endpoint_alpha_squared':sy.simplify(alpha**2),
 'endpoint_alpha_numeric':float(alpha),'endpoint_factor_numeric':float(alpha**2),
 'middle_alpha':alpha_middle,
 'weak_link_projection':'alpha^2 h + (4/5)(1-alpha^2) I',
 'scope':'exact compression of each weak link; full low-energy approximation requires weak/strong hierarchy and includes higher-order corrections'}
# Monogamy certificate in normalized code coordinates.
Om=sy.Matrix([1 if a==b else 0 for a,b in product(range(5),repeat=2)])
Bell=Om*Om.T/5
Q12=sy.kronecker_product(Bell,sy.eye(5));Q23=sy.kronecker_product(sy.eye(5),Bell)
check('monogamy',Q12*Q23*Q12==Q12/25)
check('three_chain_above_monogamy_bound',bool(e0>sy.Rational(4,9)))
print('three-cell ground code and endpoint recursion exact',flush=True)

# Numerical countercheck: the exact projected family is NOT closed under all virtual corrections.
# Everything below this point is labeled floating-point, not an exact certificate.
G5diag=np.diag(np.array(Gs,float));G3diag=np.diag(np.array(G3,float))
Wort=np.sqrt(G3diag)[:,None]*np.array(Vrg.evalf(),float)/(np.sqrt(float(NN))*np.sqrt(G5diag)[None,:])
Hn=np.sqrt(G3diag)[:,None]*np.array(H3,float)/np.sqrt(G3diag)[None,:]
Eval,EV=np.linalg.eigh(Hn)
Erel=Eval-float(e0);Erel[abs(Erel)<1e-11]=0
C0=EV.T@Wort
Beff=4/3*np.kron(C0,C0)
for Ra in R:
    Ro=np.sqrt(G5diag)[:,None]*np.array(Ra,float)/np.sqrt(G5diag)[None,:]
    BA=EV.T@np.kron(np.eye(25),Ro)@Wort
    BB=EV.T@np.kron(Ro,np.eye(25))@Wort
    Beff-=2/9*np.kron(BA,BB)
den=(Erel[:,None]+Erel[None,:]).ravel();inv=np.zeros_like(den)
inv[den>1e-10]=1/den[den>1e-10]
Hsecond=-Beff.T@(inv[:,None]*Beff)
G25diag=np.kron(G5diag,G5diag)
Kn=np.sqrt(G25diag)[:,None]*np.array(Recovery,float)/np.sqrt(G25diag)[None,:]
om=np.eye(5).ravel()/np.sqrt(5);P1=np.outer(om,om)
P9=9/4*(Kn-P1)
Swap5=np.eye(25).reshape(5,5,5,5).transpose(1,0,2,3).reshape(25,25)
P10=(np.eye(25)-Swap5)/2;P5=np.eye(25)-P1-P9-P10
fit=np.zeros((25,25));coeff={}
for dim,Pj in [(1,P1),(5,P5),(9,P9),(10,P10)]:
    cj=float(np.trace(Pj@Hsecond)/dim);fit+=cj*Pj;coeff[dim]=cj
res=float(np.max(abs(Hsecond-fit)))
check('second_order_native_sector_fit',res<1e-11,'numerical',residual=res)
check('second_order_5_10_splitting',abs(coeff[5]-coeff[10])>0.001,'numerical',difference=coeff[5]-coeff[10])
results['numerical_next_order']={'coefficients_in_units_weak_squared_over_strong':coeff,
 'native_sector_fit_max_residual':res,'split_5_minus_10':coeff[5]-coeff[10],
 'conclusion':'leading two-parameter form I,h is not a closed all-order renormalization family',
 'classification':'double-precision countercheck, no interval-certified spectral bounds'}

# General large-graph countercheck: global ferromagnetic source organization is not a set of local codes.
results['global_scope']={'symmetric_dimension_for_N_sources':'binomial(N+3,3)',
 'symmetric_dimension_8_sources':sy.binomial(11,3),
 'no_global_symmetric_exact_four_code_for_N_ge_5':'proof by transporting quartic Pauli checks and anticommuting weight-two checks',
 'micro_parent_not_rederived':['P1/P2 selection','bridge strength and signs','graph and scale hierarchy','field and state dictionary','3+1D limit','gravity']}
# Exact anticommuting witness used in the general proof.
XXI=sy.kronecker_product(sy.Matrix(X),sy.Matrix(X),sy.eye(2))
IZZ=sy.kronecker_product(sy.eye(2),sy.Matrix(Z),sy.Matrix(Z))
check('global_symmetry_obstruction_witness',XXI*IZZ==-IZZ*XXI)

results['summary']={'exact_checks':sum(c['kind']=='exact' for c in checks),'numerical_checks':sum(c['kind']=='numerical' for c in checks),
 'all_passed':True,'runtime_seconds':time.time()-T0,'python':platform.python_version(),'numpy':np.__version__,'sympy':sy.__version__,
 'meaning':'finite explicitly specified continuation, not a complete TFPT proof'}
results['checks']=checks
(ROOT/'results.json').write_text(json.dumps(clean(results),indent=2,ensure_ascii=False),encoding='utf-8')
np.savez_compressed(ROOT/'matrices.npz',occupancies=np.array(ns),metric35=m,code_E=E,code_G=G,
    paulis=np.array(As),collective_K=np.array(Ks),P5_numerator16=Pnum,
    recovery25=np.array(Recovery,float),planes_R=np.array([np.array(x,float) for x in R]),
    H3=np.array(H3,float),gram3=G3diag,ground_encoder_normalized=Wort,
    numerical_second_order=Hsecond,numerical_sector_projectors=np.array([P1,P5,P9,P10]))
(ROOT/'provenance.json').write_text(json.dumps({'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'basis_source':'TFPT_Quartik_Fortsetzung_Herleitungen_20260926.md',
    'dynamic_source':'TFPT_Dynamischer_Codeabschluss_Herleitung_20260926.md',
    'new_model':'two logical 35-response spaces, inherited 5+30 gap, chosen collective source exchange',
    'not_run':'full microscopic 56-register diagonalization or full TFPT repository'},indent=2),encoding='utf-8')
print(json.dumps(results['summary'],indent=2),flush=True)

