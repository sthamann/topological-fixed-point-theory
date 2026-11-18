"""Exact mixed response of the existing v252 source, not a source selector."""
from pathlib import Path
import sys,json,itertools
import numpy as np
R=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(R/'verification'))
import v252_full_finite_triple as f
g=f.grading();sw=f.Kswap();z=np.zeros((3,3),complex)
# Six times the actual hypercharge, using the connected original real lift.
l=f.rep(3j,np.zeros((2,2)),-1j*np.eye(3)).imag
Q=l-sw@l@sw
ys=[np.array([[1+1j,2,0],[0,3-1j,1],[1j,0,4]],complex)+j*np.eye(3) for j in range(4)]
DY=f.full_D(*ys,z)
T=g@Q;H=-1j*g@DY
words=set(itertools.permutations((0,0,1,1)))
coefficient=sum(np.trace(np.linalg.multi_dot([T if k==0 else H for k in word])) for word in words)
comm=Q@DY-DY@Q
norm=lambda a:np.vdot(a,a).real
independent=2*np.trace(Q@Q@DY@DY).real+norm(comm)
# Direct finite difference isolates the even t^2 r^2 coefficient at total degree4.
trace4=lambda a:np.trace(a@a@a@a).real
fd=(trace4(T+H)+trace4(T-H)-2*trace4(T)-2*trace4(H))/2
checks=dict(original_integer_charge=np.array_equal(Q,np.rint(Q)),original_hermitian_T=np.array_equal(T,T.conj().T),original_hermitian_H=np.array_equal(H,H.conj().T),original_grading_anticommutes_DY=not np.any(g@DY+DY@g),six_word_coefficient_real=coefficient.imag==0,independent_positive_identity=coefficient.real==independent,independent_finite_difference=fd==independent,mixed_coefficient_nonzero=coefficient.real>0,actual_gauge_higgs_commutator_nonzero=norm(comm)>0)
weights=[312,168,216,72]
checks['whole_original_yukawa_weight_identity']=independent==sum(w*norm(y) for w,y in zip(weights,ys))
unit_results=[]
for sector in range(4):
 for i in range(3):
  for j in range(3):
   for phase in (1,1j):
    args=[z.copy() for _ in range(5)];args[sector][i,j]=phase
    dy=f.full_D(*args);qdy=Q@dy
    value=2*norm(qdy)+norm(Q@dy-dy@Q)
    unit_results.append(value==weights[sector])
checks['all72_original_real_yukawa_coordinates_have_their_positive_weight']=all(unit_results)
out=dict(contract='ORIGINAL-JOINT-GAUGE-HIGGS-SOURCE-20261030',checks={k:bool(v) for k,v in checks.items()},all_passed=all(checks.values()),hypercharge_scale=6,mixed_trace4_coefficient=int(coefficient.real),charge_DY_commutator_squared_norm=int(norm(comm)),trace_Q_squared_DY_squared=int(np.trace(Q@Q@DY@DY).real),identity='coefficient(t^2 r^2) Tr(t gamma Q - i r gamma DY)^4 = 2 Tr(Q^2 DY^2) + ||[Q,DY]||_HS^2',whole_original_family_formula='312||Yu||_HS^2 +168||Yd||_HS^2 +216||Ye||_HS^2 +72||Ynu||_HS^2',scope='Actual original96 coefficient and exact general trace identity. Excludes a disjoint-support decomposition of these two original source directions and universal additive separation for all spectral test functions. Does not refute an exact cancellation in the entire original master functional or its four fixed-branch scalar equations; those require the original joint Z, reference and eta. No new physical selector or solved TFPT claim.')
Path(__file__).with_name('JOINT-GAUGE-HIGGS-SOURCE.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(out,ensure_ascii=False,indent=2));assert out['all_passed']
