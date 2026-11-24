"""A scoped all-profile spectral-action obstruction, including Einstein terms.

Let h be nonzero, continuous, h(x)/x -> 0 at infinity, and
integral h(x)*x**(-3/2) dx = 0 with ordinary endpoint convergence.
Then 1/(kappa*x+h(x)) cannot be an unsubtracted positive Stieltjes
covariance with a positive massless residue, for any real kappa.

For kappa>0, x*C(x) increases to 1/kappa under that representation.
Thus x*C<=1/kappa, H>=kappa*x and h>=0, contradicting the zero integral.
For kappa<0, the UV sign fails; for kappa=0, H=h>0 contradicts the integral.
This is an analytic proof. The tests check constants and controls, not a
finite replacement for its universal quantifiers.
"""

import sympy as sp


def spectral_no_go_audit():
    t = sp.symbols("t", positive=True)
    first = sp.simplify(3*sp.beta(sp.Rational(1, 2), sp.Rational(1, 2)).rewrite(sp.gamma))
    second = sp.simplify(12*sp.beta(sp.Rational(3, 2), sp.Rational(3, 2)).rewrite(sp.gamma))
    return {
        "status": "ANALYTIC_CONDITIONAL_NO_GO",
        "target": "flat cosmologically tuned quadratic Dirac spectral-action profile plus kappa*x",
        "profile_class": "real f with integral_0^infinity (1+u)*abs(f(u)) du finite",
        "definitions": {
            "F1(v)": "integral_v^infinity f(u) du",
            "F2(v)": "integral_v^infinity (u-v)*f(u) du",
            "h_f(x)": "3*x*integral_0^1 F1(x*t*(1-t)) dt + 12*integral_0^1 (F2(x*t*(1-t))-F2(0)) dt",
        },
        "analytic_identities": {
            "first_beta_factor": str(first), "second_beta_factor": str(second),
            "A": "integral_0^infinity v^(-1/2)*F1(v) dv",
            "B": "integral_0^infinity v^(-3/2)*(F2(v)-F2(0)) dv = -2*A",
            "sum_rule_factor": str(sp.simplify(first-2*second)),
            "absolute_integral_bound": "integral abs(h_f(x))*x^(-3/2) dx <= 12*pi*integral sqrt(u)*abs(f(u)) du",
            "uv_limit": "h_f(x) -> -6*integral u*f(u) du, hence h_f(x)/x -> 0",
        },
        "necessary_assumptions": [
            "h_f is not identically zero",
            "C(x)=Z/x+integral rho(ds)/(x+s), Z>0, rho>=0, finite for every x>0",
            "this is the full unsubtracted scalar TT covariance in a fixed source dictionary",
        ],
        "proof_by_kappa": {
            "negative": "H/x -> kappa < 0 contradicts C>0",
            "zero": "H=h>0 contradicts the convergent zero weighted integral",
            "positive": "x*C increases to 1/kappa; hence H>=kappa*x and h>=0; the sum rule forces h=0",
        },
        "einstein_addition_cannot_restore_positive_spectral_measure": True,
        "every_einstein_augmented_kernel_has_a_euclidean_zero": False,
        "euclidean_counterexample": "H_exp(x)+2*x >= 3*x*F(x) > 0 for x>0",
        "positive_control": {
            "C": "1/x+1/(x+M), M>0", "actual_uv_slope": "1/2",
            "h": "M*x/(2*(2*x+M))",
            "weighted_integral": "pi*sqrt(M)/(2*sqrt(2)) > 0",
            "why_not_excluded": "the zero sum rule does not hold",
        },
        "not_proved": ["no other fundamental gravity action exists", "full interacting quantum positivity",
                       "a conclusion for a different background or a subtracted dispersion relation"],
    }
