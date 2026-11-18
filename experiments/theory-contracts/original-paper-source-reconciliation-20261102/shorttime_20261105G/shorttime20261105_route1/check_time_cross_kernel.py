"""Exact identities for the original untruncated spectral phase coefficient.

No physical profile, Higgs measure, replacement Hamiltonian or fixture is used.
"""
import json
from pathlib import Path
import sympy as s

tau, t, p, m, chi, sigma, proper = s.symbols(
    "tau t p m chi sigma proper", positive=True, real=True
)
f = s.sqrt(s.pi) / (2 * s.sqrt(t)) * s.erfc(tau / s.sqrt(t))
tail = -2 * tau / t**2 * s.exp(-tau**2 / t)
assert s.simplify(-s.diff(f, tau, 2) - tail) == 0
assert s.simplify(s.limit(s.diff(f, tau), tau, 0, dir="+") + 1/t) == 0
assert s.integrate(2 * tau / t**2 * s.exp(-tau**2 / t), (tau, 0, s.oo)) == 1/t
assert s.simplify(2*s.integrate(f, (tau, 0, s.oo))-1) == 0

# The off-cut full Hessian includes the factor two relative to the action
# coefficient, retaining the original opposite-eigenvalue charged pair.
h_tail = -m**2 / s.pi**2 * tau / t**2 * s.exp(-t*m**2-tau**2/t)
assert s.simplify(h_tail-(m*m/(2*s.pi**2))*s.exp(-t*m*m)*tail) == 0

# Actual archive31 caller prints K_sigma/Lambda and D_ref/Lambda.
# This reduction is for the already supplied inner family K_sigma=e^-sigma D_H.
# The actual B_rel sum is not generally replaced by e^-sigma D_H.
Lambda = s.Symbol("Lambda", positive=True)
t_a = proper*s.exp(-2*sigma)/Lambda**2
t_r = proper/Lambda**2
assert s.simplify((t_a-t_r).subs(sigma, 0)) == 0
chi_geo = Lambda*s.exp(sigma)
assert s.simplify(t_a-proper/chi_geo**2) == 0

report = {
    "f_t_tau": str(f),
    "minus_second_derivative_for_tau_positive": str(tail),
    "delta_coefficient_in_minus_second_derivative": str(2/t),
    "integral_of_regular_tail_over_R": str(-2/t),
    "full_pair_hessian_regular_tail": str(h_tail),
    "archive31_matched_inner_family_t_A": str(t_a),
    "archive31_unscaled_reference_t_R": str(t_r),
    "S_even_reference_response": "0",
    "T_even_identical_whole_arguments_relative_response": "0",
    "scope": "Original spectral-action coefficient; not integrated field source or covariance.",
}
Path(__file__).with_name("time-cross-kernel-results.json").write_text(
    json.dumps(report, indent=2) + "\n"
)
print(json.dumps(report, indent=2))
