"""Check a claimed identity in a review answer, not a TFPT comparison model."""
import json
import sympy as s

a,b,tau=s.symbols('a b tau',positive=True,real=True)
C=lambda t:s.exp(-t*t/(4*tau))
det=s.factor(C(2*a)*C(2*b)-C(a+b)**2)
factor=C(a+b)**2*(s.exp(-(a-b)**2/(2*tau))-1)
assert s.simplify(det-factor)==0
specific=s.simplify(det.subs({a:1,b:2,tau:1}))
assert specific.is_negative is True

# Any nonzero positive measure rho on finite oscillator squared masses
# has rho([0,M])>0 for some finite M. Thus
# z*int rho(dm)/(z+m) >= rho([0,M])*z/(z+M),
# which cannot converge to zero. In contrast z*exp(-tau*z) -> 0.
z,M,w=s.symbols('z M w',positive=True)
assert s.limit(z*s.exp(-tau*z),z,s.oo)==0
assert s.limit(w*z/(z+M),z,s.oo)==w

result={
 'result_id':'TORSION20261104G',
 'scope':'Exact rejection of the proposed heat-to-positive-oscillator identification in the external review; no change to the actual TFPT source.',
 'reflection_determinant_factor':str(factor),
 'two_time_determinant_exact':str(specific),
 'two_time_determinant_value':float(specific),
 'positive_measure_resolvent_lower_bound_limit':str(w),
 'heat_scaled_limit':'0',
 'claimed_heat_covariance_RP':False,
 'claimed_positive_oscillator_representation_of_single_heat_multiplier':False,
 'full_TFPT_reflection_positivity_disproved':False,
 'all_possible_exact_auxiliary_representations_excluded':False,
 'full_original_source_evaluated':False,
 'fundamental_solution_achieved':False
}
print(json.dumps(result,ensure_ascii=False,indent=2))
