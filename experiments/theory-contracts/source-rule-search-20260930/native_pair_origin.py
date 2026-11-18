"""Derive neutral pair response from actual charged E8 source roots.

No coupling, primary projection or sewing time is selected by this check.
The root-pair OPE is the lattice-vertex exponential expansion; all root
moments and oscillator states are read from the existing source.
"""
import json
import sys
from pathlib import Path
import sympy as s

REPO = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0, str(REPO))
from tfpt_explorer.neutral_source_response import _native_neutral_branch, _quadratic_state
from tfpt_explorer.source_charge_response import exact_charge_response

n = _native_neutral_branch()
roots = [s.Matrix(r)/2 for r in sorted(n['retained'])]
assert len(roots) == 48
assert all((b.T*b)[0] == 2 for b in roots)
# With the adjoint pair normalized to positive leading coefficient,
# t^2 E_b(t)E_b^dagger(0) Omega
# = Omega+t b(-1)Omega+t^2/2*(b(-1)^2+b(-2))Omega+...
# Adjoint symmetrization cancels b(-2), not the stress component.
moment = sum((b*b.T for b in roots), s.zeros(8))
linear = sum(roots, s.zeros(8,1))
derivative_grade2 = linear/2
pair = _quadratic_state(moment)
T = n['states']['T']
gram = n['gram']
inner = lambda a,b: s.factor((a.T*gram*b)[0])
stress_coefficient = inner(T,pair)/inner(T,T)
primary = pair-stress_coefficient*T
assert moment == n['A48']
assert stress_coefficient == 12
assert primary == n['states']['R']
assert inner(pair,pair) == 624
assert inner(primary,primary) == 48
assert inner(12*T,12*T) == 576
assert inner(T,primary) == 0
assert s.trace(moment)/2 == 48
assert derivative_grade2 != s.zeros(8,1)
# The 64-root chiral spinor subset is a decisive control: the same
# construction yields pure stress, hence its primary projection vanishes.
all_roots = [s.Matrix(r)/2 for r in n['spinor']]
all_moment = sum((b*b.T for b in all_roots),s.zeros(8))
all_pair = _quadratic_state(all_moment)
assert all_pair == 16*T
# Charge weighting follows from the same pair operation, no R-only guess.
charge = exact_charge_response()
y = charge['hypercharge']
weighted = sum((((y.T*b)[0])**2*b*b.T for b in roots),s.zeros(8))
weighted_stress = s.trace(weighted)/8
weighted_primary = weighted-weighted_stress*s.eye(8)
assert weighted_primary == charge['R_Y']
assert weighted_stress == s.Rational(5,2)
assert (y.T*n['R_matrix']*y)[0] == 0
assert (y.T*weighted_primary*y)[0] == s.Rational(115,36)
# A common microscopic coupling g induces g^2 in this pair and g^4 in
# its norm. No equation here sets g=c3 or alpha=radial sewing length.
g,tau = s.symbols('g tau',real=True)
result = {
 'method':'derive charged pair first; do not install a chosen neutral port',
 'pair_expansion':'t^2 E_b(t)E_b^dagger(0) = 1+t b(-1)+t^2/2*(b(-1)^2+b(-2))+...',
 'actual_charged_roots':48,
 'ordered_pair_linear_grade1_coefficient':list(map(str,linear)),
 'ordered_pair_derivative_grade2_coefficient':list(map(str,derivative_grade2)),
 'adjoint_even_pair':'B48=12*T+R',
 'pair_norm_squared':str(inner(pair,pair)),
 'stress_norm_squared':str(inner(12*T,12*T)),
 'primary_norm_squared':str(inner(primary,primary)),
 'connected_identity_subtraction_removes_stress':False,
 'primary_projection_required':'orthogonal projection away from total Sugawara stress T',
 'all_64_control':'B64=16*T, primary projection zero',
 'charge_weighted_pair':'BY=(5/2)*T+RY',
 'R_JY_JY':str(charge['R_JY_JY']),
 'RY_JY_JY':str(charge['RY_JY_JY']),
 'unprojected_pair_heat':str(624*g**4*s.exp(-2*tau)),
 'primary_pair_heat':str(48*g**4*s.exp(-2*tau)),
 'not_derived':['physical source identification','adjoint-even rather than ordered operation','physical rule excluding total stress in this response','g=c3','tau=physical alpha','mark-to-local-vertex coupling','3+1D transfer','prime law','AI superiority'],
 'verdict':'exact native pair decomposition; no unconditional fundamental solution'
}
Path(__file__).with_name('native_pair_origin_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(result,ensure_ascii=False,indent=2))
