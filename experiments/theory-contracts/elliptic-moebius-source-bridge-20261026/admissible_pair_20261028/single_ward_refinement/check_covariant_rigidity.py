#!/usr/bin/env python3
"""Exact channel rigidity with fixed input/output representations, not P1 selection."""
from pathlib import Path
import importlib.util
import itertools as it
import hashlib
import json
import numpy as np
import sympy as sp

HERE=Path(__file__).resolve().parent
PREV=HERE.parent/'research-20261028-hodge-complement-channel/hodge_complement_channel.py'
if not PREV.is_file(): PREV=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/elliptic-moebius-source-bridge-20261026/hodge_complement_20261028/hodge_complement_channel.py')
spec=importlib.util.spec_from_file_location('hodge_previous',PREV)
h=importlib.util.module_from_spec(spec); spec.loader.exec_module(h)
C=h.hodge_amplitudes(); pairs=h.PAIRS
checks=[]
def require(v,name):
 if not bool(v): raise RuntimeError(name)
 checks.append(name)
def ptout(M): return np.einsum('aiaj->ij',M.reshape(10,10,10,10))
def choi_to_super(J): return J.reshape(10,10,10,10).transpose(0,2,1,3).reshape(100,100)
I=np.eye(100,dtype=np.int64); F=np.zeros((100,100),dtype=np.int64)
for a in range(10):
 for b in range(10): F[10*b+a,10*a+b]=1
vs=[x.reshape(-1) for x in C]
M=sum(np.outer(v,v) for v in vs)
Q=[M,3*(I-F),3*(I+F)-M] # six times orthogonal projectors
D=[5,45,50]
for i,q in enumerate(Q):
 require(np.array_equal(q@q,6*q),f'projector {D[i]} exact idempotence')
 require(np.trace(q)==6*D[i],f'projector {D[i]} exact rank')
 require(np.array_equal(ptout(q)*10,6*D[i]*np.eye(10,dtype=np.int64)),f'projector {D[i]} TP partial trace')
 for j,r in enumerate(Q):
  if i!=j: require(not np.any(q@r),f'orthogonality {D[i]}/{D[j]}')
require(np.array_equal(sum(Q),6*I),'complete three-projector decomposition')
# Actual complexified dual representation on bar10 tensor bar10.
G={}
for a in range(5):
 for b in range(5):
  X=np.zeros((5,5),dtype=np.int64); X[b,a]=1
  r=-h.exterior_two_differential(X)
  G[a,b]=np.kron(r,np.eye(10,dtype=np.int64))+np.kron(np.eye(10,dtype=np.int64),r)
for a in range(5):
 for b in range(5):
  for d,q in zip(D,Q): require(np.array_equal(G[a,b]@q,q@G[a,b]),f'dual gl5 E{a}{b} commutes P{d}')
cas=5*sum(G[a,b]@G[b,a] for a in range(5) for b in range(5))-16*I
for ev,q in zip([24,64,84],Q): require(np.array_equal(cas@q,ev*q),f'quadratic Casimir eigenvalue {ev}/5')
raising=sp.Matrix(np.vstack([G[i,i+1] for i in range(4)]))
N=raising.nullspace()
require(len(N)==3,'three highest-weight lines in full tensor product')
Nmat=sp.Matrix.hstack(*N)
for d,q in zip(D,Q): require((sp.Matrix(q)*Nmat).rank()==1,f'one highest weight in P{d}; irreducible')
# Choi vertex transition law. Every population entry checked as rational.
transitions=[]
for q,d in zip(Q,D):
 vals={k:set() for k in ['stay','intersect','disjoint']}
 for a,e in enumerate(pairs):
  for b,f in enumerate(pairs):
   kind='stay' if a==b else 'disjoint' if set(e).isdisjoint(f) else 'intersect'
   vals[kind].add(sp.Rational(10*int(q[10*b+a,10*b+a]),6*d))
 expected={5:{'stay':0,'intersect':0,'disjoint':sp.Rational(1,3)},45:{'stay':0,'intersect':sp.Rational(1,9),'disjoint':sp.Rational(1,9)},50:{'stay':sp.Rational(1,5),'intersect':sp.Rational(1,10),'disjoint':sp.Rational(1,15)}}[d]
 for k in vals: require(vals[k]=={expected[k]},f'P{d} all {k} population entries')
 transitions.append({k:str(v) for k,v in expected.items()})
p5,p45,p50,s,o=sp.symbols('p5 p45 p50 s o')
sol=sp.solve([p5+p45+p50-1,s-p50/5,o-sp.Rational(2,3)*p45-sp.Rational(3,5)*p50],[p5,p45,p50])
require(sol[p50]==5*s and sol[p45]==(3*o-9*s)/2,'leakage reconstructs three covariant weights')
require(sp.simplify(2*(sol[p45]+sol[p50])-(3*o+s))==0,'diamond distance from orthogonal Choi sectors equals 3o+s')
# Analytic norm proof: T-Phi=sum_{45,50}p_i(T_i-Phi) has diamond norm
# at most 2*sum(p_i), since all three vertices are CPTP. The maximally
# entangled probe gives exactly that trace norm on the orthogonal Choi
# supports. The lower and upper bounds match, without a numerical SDP.
# SU3 plus actual five-cycle Lie algebra generation, with exact matrices.
cycle=np.zeros((5,5),dtype=np.int64)
for j in range(5): cycle[(j+1)%5,j]=1
roots=set()
for k in range(5):
 for a,b in [(0,1),(1,0),(1,2),(2,1),(0,2),(2,0)]: roots.add(((a+k)%5,(b+k)%5))
change=True
while change:
 change=False
 for a,b in list(roots):
  for c,d in list(roots):
   if b==c and a!=d and (a,d) not in roots: roots.add((a,d)); change=True
require(len(roots)==20,'SU3 and fivecycle conjugates generate all twenty sl5 root directions')
# Fixed-dictionary SM phase counterexample: rank five and same populations.
U=np.diag([-1 if all(j<3 for j in e) else 1 for e in pairs]).astype(np.int64)
Mu=sum(np.outer((c@U).reshape(-1),(c@U).reshape(-1)) for c in C)
require(np.array_equal(np.diag(Mu),np.diag(M)),'SM phase variant has identical Petersen populations')
require(np.array_equal(ptout(Mu),3*np.eye(10,dtype=np.int64)),'SM phase variant trace preserving')
require(np.array_equal(Mu@Mu,6*Mu) and np.trace(Mu)==30,'SM phase variant has Choi rank five')
require(not np.array_equal(Mu,M),'SM phase variant is a different quantum channel')
sm=[]
for a,b in it.product(range(5),repeat=2):
 if (a<3)==(b<3):
  X=np.zeros((5,5),dtype=np.int64); X[a,b]=1
  r=h.exterior_two_differential(X); sm.append(r)
  require(np.array_equal(r@U,U@r),f'phase counterexample commutes block E{a}{b}')
Sj=choi_to_super(M); Su=choi_to_super(Mu)
for j,r in enumerate(sm):
 # Complexified covariance input r, target dual -r.T.
 Lin=np.kron(r,np.eye(10,dtype=np.int64))-np.kron(np.eye(10,dtype=np.int64),r.T)
 rout=-r.T
 Lout=np.kron(rout,np.eye(10,dtype=np.int64))-np.kron(np.eye(10,dtype=np.int64),rout.T)
 require(np.array_equal(Su@Lin,Lout@Su),f'SM covariance phase variant generator {j}')
R=np.zeros((10,10),dtype=np.int64)
for j,(a,b) in enumerate(pairs):
 aa,bb=(a+1)%5,(b+1)%5; f=tuple(sorted((aa,bb)))
 R[pairs.index(f),j]=1 if aa<bb else -1
RR=np.kron(R,R)
require(np.array_equal(RR@M@RR.T,M),'Hodge actual wedge fivecycle Choi covariance')
require(not np.array_equal(RR@Mu@RR.T,Mu),'SM phase variant violates actual wedge fivecycle covariance')
prob=sp.Rational(int(np.trace(M@Mu)),180)
require(prob==sp.Rational(4,15),'neutral P5 witness separates SM phase variant with probability 4/15')
require(sp.Rational(int(np.trace(M@M)),180)==1,'neutral P5 Hodge witness probability one')
# Same invariant family can transport across every color marking. Covariance
# between different marked fibers is weaker than covariance in one fiber.
new_colors={(i+1)%5 for i in range(3)}
Unew=np.diag([-1 if set(e)<=new_colors else 1 for e in pairs]).astype(np.int64)
require(np.array_equal(Unew,R@U@R.T),'phase construction transports to moved color mark')
Mu_new=sum(np.outer((c@Unew).reshape(-1),(c@Unew).reshape(-1)) for c in C)
require(np.array_equal(Mu_new,RR@Mu@RR.T),'non-Hodge channel family passes actual fivecycle transport between marks')
# A real, unital, HS-selfadjoint variant also defeats spectral selection.
Mtwo=sum(np.outer((U@c@U).reshape(-1),(U@c@U).reshape(-1)) for c in C)
Stwo=choi_to_super(Mtwo)
AU=np.kron(U,U)
require(np.array_equal(Stwo,AU@Sj@AU),'two-sided phase variant is exact superoperator similarity in fixed real chart')
require(np.array_equal(Stwo,Stwo.T),'two-sided phase variant HS selfadjoint')
require(np.array_equal(np.diag(Mtwo),np.diag(M)),'two-sided phase variant same Petersen population')
require(np.array_equal(ptout(Mtwo),3*np.eye(10,dtype=np.int64)),'two-sided phase variant TP')
require(np.array_equal(Mtwo@Mtwo,6*Mtwo) and np.trace(Mtwo)==30,'two-sided phase variant Choi rank five')
require(not np.array_equal(Mtwo,M),'two-sided phase variant different in fixed input/output dictionary')
for c in C:
 require(np.array_equal((U@c@U).T,U@c@U),'two-sided Kraus real symmetric')
 require(np.array_equal((U@c@U).T@(U@c@U),c.T@c),'two-sided phase variant same individual effects')
prob_two=sp.Rational(int(np.trace(M@Mtwo)),180)
require(prob_two==sp.Rational(7,15),'two-sided phase variant neutral P5 probability seven fifteenths')
Mtwo_new=sum(np.outer((Unew@c@Unew).reshape(-1),(Unew@c@Unew).reshape(-1)) for c in C)
require(np.array_equal(Mtwo_new,RR@Mtwo@RR.T),'two-sided variant also passes fivecycle transport between marks')
for n in range(1,7):
 require(np.array_equal(np.linalg.matrix_power(Stwo,n),AU@np.linalg.matrix_power(Sj,n)@AU),'two-sided phase variant power similarity '+str(n))
# Source U5 is the P2 degree phase exp(2pi i X/5), not the cycle above.
X=np.diag([-2,-2,-2,3,3]).astype(np.int64)
Xtwo=h.exterior_two_differential(X)
require(all(int(q)%5==3 for q in np.diag(X)),'original U5 is scalar zeta5^3 on carrier V')
require(all(int(q)%5==1 for q in np.diag(Xtwo)),'original U5 is scalar zeta5 on duads')
require(not np.array_equal(cycle@X@cycle.T,X),'literal slot cycle moves original P2 charge mark')
require(U[0,0]*U[9,9]!=U[pairs.index((0,3)),pairs.index((0,3))]*U[pairs.index((1,4)),pairs.index((1,4))],
        'CC phaseflip violates exterior square phase relation and is not a carrier torus gauge')
# Stronger source criterion: no SU5 covariance premise is required if a
# positive TP Choi saturates its rank-five process witness. Positivity gives
# support inside P5, so J=sum K_mn |C_m>><<C_n| with K positive. Its TP map
# (actual Choi output partial trace) is tr(K)I-dGamma(K.T). The adjoint
# effect T*(I) is its transpose, tr(K)I-dGamma(K). Both TP conditions are
# equivalent, but we test the actual partial trace, including complex K.
tp_columns=[]
for m,n in it.product(range(5),repeat=2):
 E=np.zeros((5,5),dtype=np.int64); E[m,n]=1
 lhs=C[m].T@C[n]
 rhs=(np.eye(10,dtype=np.int64) if m==n else np.zeros((10,10),dtype=np.int64))-h.exterior_two_differential(E.T)
 require(np.array_equal(lhs,rhs),f'TP Choi coefficient A(E{m}{n})=trE I-dGamma(E.T)')
 require(np.array_equal(ptout(np.outer(vs[m],vs[n])),lhs),f'actual Choi partial trace coefficient {m}{n}')
 tp_columns.append(lhs.reshape(-1))
tp_matrix=sp.Matrix(np.stack(tp_columns,axis=1))
require(tp_matrix.rank()==25,'all 25 source Choi coefficients fixed injectively by trace preservation')
target=sp.Matrix(np.eye(10,dtype=np.int64).reshape(-1))
tp_solution=tp_matrix.gauss_jordan_solve(target)[0]
expected=sp.Matrix((np.eye(5,dtype=np.int64).reshape(-1))).applyfunc(lambda x:x/3)
require(tp_solution==expected,'unique supported TP Choi has K=I5/3 without covariance assumption')
K2=2*np.eye(5,dtype=np.complex128);K2[0,1]=1j;K2[1,0]=-1j
V=np.stack(vs,axis=1)
partial=ptout(V@K2@V.T)
direct=sum(K2[m,n]*(C[m].T@C[n]) for m in range(5) for n in range(5))
require(np.array_equal(partial,direct),'complex Hermitian K Choi partial trace follows transposed convention')
wrong=sum(K2[m,n]*(C[n].T@C[m]) for m in range(5) for n in range(5))
require(not np.array_equal(partial,wrong),'complex K distinguishes Choi partial trace from adjoint TP effect transpose')
# Word expansion is explicit finite rank-5 neutral pair projector, not asserted P1 image.
result={
 'scope':'Fixed H=Lambda2 C5, target dual10; exact representation/channel theorem, not original P1 source selection.',
 'choi_projector_dimensions':D,'scaled_casimir_eigenvalues':[24,64,84],
 'highest_weight_lines':len(N),'covariant_simplex_vertex_populations':transitions,
 'rigidity':'SU5 conjugate covariance and Petersen support imply unique Phi, including all coherences.',
 'weights_from_leakage':{str(k):str(v) for k,v in sol.items()},
 'diamond_distance':'3 * intersecting_leakage + stay_leakage (within exact covariant simplex)',
 'SM_plus_actual_wedge_fivecycle_generates_SU5':True,
 'SM_alone_counterexample':{'U':'minus on CC, plus on CW and WW','choi_rank':5,'same_population':True,'P5_probability':str(prob),'wedge_cycle_covariance':False},
 'mark_transport_does_not_select_channel':{'one_sided_family':True,'two_sided_family':True,'scope':'Full family covariance over Gr(3,5) is weaker than full covariance at one fixed mark.'},
 'real_selfadjoint_isospectral_counterexample':{'Kraus':'U C_m U / sqrt(3), fixed real marked chart','P5_probability':str(prob_two),'same_full_spectrum_and_effects':True,'physical_dictionary_fixed':True},
 'original_U5':{'carrier':'zeta5^3 I5','duads':'zeta5 I10','mixes_slots':False,'literal_cycle_moves_charge_mark':True},
 'stronger_no_covariance_selection':{'criterion':'T CPTP and Tr(P5 J_T/10)=1','unique_channel':'Phi','TP_Choi_partial_trace_map':'K -> tr(K) I10 - dGamma(K.T)','TP_adjoint_effect_map':'K -> tr(K) I10 - dGamma(K)','coefficient_rank':25,'unique_K':'I5/3','physical_P1_process_saturation_proved':False},
 'source_sha256':hashlib.sha256(PREV.read_bytes()).hexdigest(),
 'checks':len(checks),'full_check_labels':checks}
(HERE/'covariant-rigidity-results.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='full_check_labels'},indent=2,ensure_ascii=False))
