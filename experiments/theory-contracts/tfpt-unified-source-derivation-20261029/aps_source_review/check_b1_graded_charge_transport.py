"""Decisive invariant for the existing actual v252 real-gauge dictionary.

Compare the specified grading/charge tuples, not all conceivable physical
carrier maps. No anomaly-free condition is substituted for source derivation.
"""
from pathlib import Path
import json
import sys
import numpy as np
import sympy as s

ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'verification'))
import v252_full_finite_triple as finite
K=finite.Kswap()
# Infinitesimal original connected unimodular lift:
# lambda=exp(it/2), m=exp(-it/6)I3, q=I2.
# Evaluate its derivative multiplied by 6, using Gaussian integers.
left6=finite.rep(3j,np.zeros((2,2)), -1j*np.eye(3)).imag
Y6=left6-K@left6@K
assert np.all(Y6==np.round(Y6))
charges6=[int(x) for x in np.diag(Y6)]
g=s.diag(*[int(x) for x in np.diag(finite.grading())])
q48=[int(name in ('u','d')) for name,ch,colour,family in finite.P]
q=s.diag(*(q48+q48))
gp=g*(2*q-s.eye(96))
Y=s.diag(*[s.Rational(x,6) for x in charges6])
def trace_power(grading,k):
    return sum(grading[i,i]*Y[i,i]**k for i in range(96))
def spectrum(grading):
    counts={}
    for i in range(96):
        key=f'grade={grading[i,i]},Y={Y[i,i]}'
        counts[key]=counts.get(key,0)+1
    return counts
actual=trace_power(g,3)
alternative=trace_power(gp,3)
sign_class=[]
for sq,sl in ((1,1),(1,-1),(-1,1),(-1,-1)):
    grading=g*(sq*q+sl*(s.eye(96)-q))
    sign_class.append({'quark_sign':sq,'lepton_sign':sl,
                       'trace3':str(trace_power(grading,3)),
                       'trace5':str(trace_power(grading,5)),
                       'same_joint_multiplicities_as_original':spectrum(grading)==spectrum(g)})
expected_particle=[]
for name,ch,colour,family in finite.P:
    expected_particle.append(({'nu':-3,'e':-3,'u':1,'d':1} if ch=='L' else {'nu':0,'e':-6,'u':4,'d':-2})[name])
checks={
    'original_real_gauge_charge_diagonal':s.Matrix([[int(x) for x in row] for row in Y6])==6*Y,
    'original_real_gauge_matches_actual_SM_field_charges':charges6==expected_particle+[-x for x in expected_particle],
    'original_total_cubic_graded_trace_zero':actual==0,
    'quark_cubic_trace_minus_9_over_2':s.trace(q*g*Y**3)==-s.Rational(9,2),
    'lepton_cubic_trace_plus_9_over_2':s.trace((s.eye(96)-q)*g*Y**3)==s.Rational(9,2),
    'alternative_cubic_graded_trace_minus_9':alternative==-9,
    'no_unitary_simultaneously_intertwines_specified_grading_and_charge':actual!=alternative,
    'joint_graded_charge_multiplicities_different':spectrum(g)!=spectrum(gp),
    'original_fifth_graded_trace_10_over_3':trace_power(g,5)==s.Rational(10,3),
    'only_original_sign_choice_preserves_fixed_joint_dictionary':sum(x['same_joint_multiplicities_as_original'] for x in sign_class)==1 and sign_class[0]['same_joint_multiplicities_as_original'],
}
result={'contract':'B1-GRADED-CHARGE-TRANSPORT-20261030','checks':checks,
        'original_cubic_graded_trace':str(actual),'alternative_cubic_graded_trace':str(alternative),
        'original_joint_multiplicities':spectrum(g),'alternative_joint_multiplicities':spectrum(gp),
        'complete_four_sector_sign_class':sign_class,
        'scope':'No unitary U with U*Y U=Y and U*gamma U=gamma(2q-I). Applies to the specified graded charge tuple only. Does not exclude alternate source maps with another transported intrinsic grading, extra sectors or inflow. Does not derive physical P_F.',
        'all_passed':all(checks.values())}
Path(__file__).with_name('B1-GRADED-CHARGE-TRANSPORT.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
assert result['all_passed']
