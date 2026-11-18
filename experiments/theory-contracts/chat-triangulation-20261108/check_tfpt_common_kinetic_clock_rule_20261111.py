"""Exact original Fourier profile test of the canonical Higgs energy sum rule."""
import json
from pathlib import Path
import sys

import sympy as sp

repo = Path(sys.argv[1]).resolve()
domain = repo / (
    "experiments/theory-contracts/original-paper-source-reconciliation-20261102"
    "/contribution_20261106Q/domain/native-diagnostic-v1.json"
)
profile = json.loads(domain.read_text())["field_proposal"]
x = sp.Symbol("x", real=True)
fields = list(map(sp.Rational, profile["base_constant"]))
for slot in profile["slots"]:
    trig = {"cos": sp.cos, "sin": sp.sin}[slot["kind"]]
    fields[slot["component"]] += sp.Rational(slot["unit"])*trig(slot["harmonic"]*x)

r = sp.Rational(13, 100)
a = sp.Symbol("a", real=True)
fields[0] += a*sp.cos(8*x)


def fourier(expr):
    terms = {}
    for term in sp.expand(expr.rewrite(sp.exp)).as_ordered_terms():
        coefficient, wave = term.as_independent(x, as_Add=False)
        k = 0 if wave == 1 else sp.simplify(wave.args[0]/(sp.I*x))
        assert k in sp.S.Integers
        terms[int(k)] = terms.get(int(k), 0)+coefficient
    return terms


def projection(expr):
    return sp.simplify(2*fourier(expr).get(7,0).as_real_imag()[0])


norm = sum(f*f for f in fields)
mark = projection(norm)
b = sp.Rational(profile["slots"][0]["unit"])
assert mark == b*a
one = b**2/4    # In units 1/(pi R Z).
two = (b**2+r**2)/4
geometry_coefficients = {}
for sign in (-1, 1):
    at = {a: sign*r}
    u = fourier(norm.subs(at))
    assert 14 not in u and u[0] == sp.Rational(5117, 5000)
    complete = sp.simplify(u[0]+sp.re(u.get(14, sp.S.Zero)))
    # Independent direct differentiation of the original mark in every
    # real coordinate of its complete gradient, with canonical metric.
    h = sp.Symbol("h", real=True)
    basis = [(sp.S.One, 2)] + [(trig(k*x),1) for k in range(1,16)
                                for trig in (sp.cos,sp.sin)]
    gradient = []
    for idx, field in enumerate(fields):
        for wave, squared_norm in basis:
            # U7 is quadratic, so the centered difference is exact.
            trial = norm-field**2+(field+h*wave)**2
            derivative = sp.diff(projection(trial), h).subs({**at, h:0})
            if derivative:
                gradient.append((idx, str(wave), derivative, squared_norm))
    direct = sum(val**2/(4*squared_norm) for _,_,val,squared_norm in gradient)
    assert direct == complete
    assert sum(val**2/4 for idx,w,val,_ in gradient
               if idx==0 and w in ("cos(x)", "cos(8*x)")) == two
    geometry_coefficients[str(sign)] = {
        "one_moving_mode": str(one),
        "two_moving_modes": str(two),
        "whole_local_gradient": str(complete),
        "whole_over_one": str(complete/one),
        "gradient_components": [[i,w,str(val),n] for i,w,val,n in gradient],
    }

# The prior continuations have the same equal-time norm and unequal
# energy-weighted sums at the very same time unit.
mark_variance = (b*r)**2
lambda_a, lambda_b = sp.Rational(1,4), sp.Rational(1,16)
slopes_times_tau = [sp.expand_log(-mark_variance*sp.log(rate),force=True)
                    for rate in (lambda_a,lambda_b)]
assert sp.simplify(slopes_times_tau[1]/slopes_times_tau[0]) == 2

# Positive spectra can share an equal-time norm and a first moment.
# The second moment detects the obstruction to a single clock rate.
spectra = [(sp.S.One,), (sp.Rational(1,2),sp.Rational(3,2))]
moments = [[sum(energy**n for energy in spec)/len(spec) for n in range(3)]
           for spec in spectra]
assert moments[0][:2] == moments[1][:2] == [1,1]
assert [m[0]*m[2]-m[1]**2 for m in moments] == [0,sp.Rational(1,4)]

print(json.dumps({
    "U7_on_original_profile": str(mark),
    "kinetic_coefficients_in_units_one_over_pi_R_Z": geometry_coefficients,
    "prior_continuation_slopes_times_same_tau": list(map(str,slopes_times_tau)),
    "spectral_moment_counterexample": [list(map(str,m)) for m in moments],
    "scope": "Exact Fourier gradients. Application as a time sum rule requires the local canonical generator and its physical source identification.",
}, indent=2))
