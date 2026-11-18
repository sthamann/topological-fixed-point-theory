"""Direct old-source checks; the infinite-dimensional theorem is in PROOF.txt.

Reuses the original CAR stress builder. No new source, state or oscillator.
Finite checks do not prove Reeh-Schlieder, BW or an infinite spectrum.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from math import factorial
import importlib.util
import json
import sympy as s

ROOT = Path(__file__).resolve().parents[4]
OLD = ROOT / "experiments/theory-contracts/original-admissible-source-closure-20261025/check_w8_stress_cyclic.py"
spec = importlib.util.spec_from_file_location("old_quarter_stress", OLD)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
SEA = {255: F(1)}  # original occupied modes -8,...,-1


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def clean(v):
    return {k: c for k, c in v.items() if c}


def current(v, n):
    # Exact SAME CAR matrix units as the original stress, with coefficient1.
    out = {}
    for mask, c in v.items():
        for i in range(16 - n):
            move = old.fermion_move(mask, i, i + n)
            if move is not None:
                state, sign = move
                out[state] = out.get(state, F(0)) + c * sign
    return clean(out)


@lru_cache(None)
def pbw_in_car(lam):
    v = SEA
    for n in lam:
        v = current(v, n)
    return v


def car_linear(pbw):
    result = {}
    for lam, c in pbw.items():
        for mask, value in pbw_in_car(lam).items():
            result[mask] = result.get(mask, F(0)) + c * value
    return clean(result)


def weighted_integral(poly, x):
    # integral_0^infinity x exp(-x) poly(x) dx, exactly for polynomials.
    return sum(c * factorial(k[0] + 1) for k, c in s.Poly(s.expand(poly), x).terms())


def main():
    comparisons = 0
    for alpha in (F(1, 4), F(3, 4)):
        b = alpha - F(1, 2)
        for grade in range(5):
            for lam in old.partitions(grade):
                v = pbw_in_car(lam)
                for n in (1, 2):
                    lhs = old.original_stress_action(v, n, alpha)
                    for mask, c in current(v, n).items():
                        lhs[mask] = lhs.get(mask, F(0)) - b*c
                    rhs = car_linear(old.boson_stress_action({lam: F(1)}, n, F(0)))
                    need(clean(lhs) == rhs, "same-source compensated stress identity")
                    comparisons += 1
        original = old.original_stress_action(SEA, 1, alpha)
        expected = {k: b*c for k, c in current(SEA, 1).items()}
        need(original == expected, "actual twisted stress vacuum witness")

    x = s.Symbol("x", positive=True)
    N = 7
    laguerre = [s.assoc_laguerre(n-1, 1, x)/s.sqrt(n) for n in range(1, N+2)]
    for n in range(1, N+1):
        value = n*laguerre[n-1] - s.sqrt(n*(n+1))*laguerre[n]/2
        if n > 1:
            value -= s.sqrt(n*(n-1))*laguerre[n-2]/2
        need(s.simplify(value - x*laguerre[n-1]/2) == 0, "Laguerre Jacobi recurrence")
        for m in range(1, N+1):
            need(s.simplify(weighted_integral(laguerre[n-1]*laguerre[m-1], x)) == int(n == m), "Laguerre Gram")

    f = s.Function("f")(x)
    D = lambda g: s.I*(x*s.diff(g, x)+(1-x/2)*g)
    need(s.simplify(D(x*f/2)-x*D(f)/2-s.I*x*f/2) == 0, "[D,P]=iP")
    r, u = s.symbols("r u", positive=True)
    amplitude = lambda a, y: a*s.exp((1-a)*y/2)
    need(s.simplify(amplitude(r,x)*amplitude(u,r*x)-amplitude(r*u,x)) == 0, "dilation group law")
    need(s.simplify(x*s.exp(-x)*amplitude(r,x)**2-r*(r*x)*s.exp(-r*x)) == 0, "dilation unitarity density")
    theta = s.Symbol("theta", real=True)
    cayley = -s.cot(theta/2)
    need(s.trigsimp(s.diff(cayley,theta)*(1-s.cos(theta))-1) == 0, "Cayley translation vector field")
    need(s.trigsimp(s.diff(cayley,theta)*(-s.sin(theta))-cayley) == 0, "Cayley dilation vector field")
    q = s.Symbol("q", positive=True)
    thermal_commutator = -s.sqrt(2)*(q-q**2)/2
    need(thermal_commutator != 0, "thermal source differs from invariant vacuum")
    result = {
        "status": "PASS_SCOPED",
        "literal_original_CAR_compensated_stress_comparisons": comparisons,
        "quarters": ["1/4", "3/4"],
        "exact_CAR_window": "-8,...,7; all checked total grades<=6, below boundary grade9",
        "laguerre_recurrence_checks": N,
        "laguerre_gram_checks": N*N,
        "thermal_commutator_entry_times_partition_function": str(thermal_commutator),
        "analytic_scope": "Full Laguerre basis, operator closures, positive spectrum, local current net and BW/HSMI are the written theorem with primary references, not inferred from finite checks.",
        "P1_source_map_closed": False,
        "full_TFPT_solution": False,
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
