#!/usr/bin/env python3
"""Exact separated two-point function of the existing shifted-current W8.

The state and its coefficients come from check_quantum_octic_completion.py.
No state, coupling, projection, or direction is fitted here.  For each of the
two original quarter folds, the checker expands

    z^8 Y_d(W8,z)|d> = sum_(ell>=0) v_ell z^ell

through level eight, computes a_ell=||v_ell||^2, and reconstructs the entire
separated two-point function from exchange symmetry and the collision bound.
"""

from collections import Counter
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, product
from math import comb, factorial, prod
from pathlib import Path
import importlib.util
import json
import sys

import sympy as s


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
COMPLETION = HERE / "quantum-octic-completion-results.json"
OLD_SOURCE = (
    ROOT
    / "experiments/theory-contracts/native-octic-seam-origin-20261025"
    / "native_folded_octic_origin.py"
)
MAX_LEVEL = 8


def need(condition, message):
    if not condition:
        raise ArithmeticError(message)


def as_fraction(value):
    value = s.Rational(value)
    return Fraction(int(value.p), int(value.q))


def rational_text(value):
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def load_old_source():
    sys.path.insert(0, str(ROOT))
    spec = importlib.util.spec_from_file_location("old_octic_source", OLD_SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def standard_roots():
    roots = []
    for i, j in combinations(range(8), 2):
        for a, b in product((-1, 1), repeat=2):
            root = [s.S.Zero] * 8
            root[i], root[j] = a, b
            roots.append(tuple(root))
    roots.extend(
        tuple(s.Rational(sign, 2) for sign in signs)
        for signs in product((-1, 1), repeat=8)
        if sum(sign < 0 for sign in signs) % 2 == 0
    )
    need(len(set(roots)) == 240, "standard E8 root census changed")
    return tuple(roots)


def build_state():
    """Rebuild W8, taking the solved coefficients from the owned JSON."""
    result = json.loads(COMPLETION.read_text())
    coefficients = [s.Rational(value) for value in result["coefficients"]]
    h = {n: s.symbols(f"h{n}_0:8") for n in range(1, 8)}
    variables = tuple(value for n in range(1, 8) for value in h[n])
    x = h[1]

    p8 = s.S.Zero
    for indices in combinations_with_replacement(range(8), 4):
        powers = Counter(indices)
        coefficient = s.Rational(
            factorial(8),
            2 * prod(factorial(2 * powers[j]) for j in range(8)),
        )
        p8 += coefficient * prod(x[j] ** (2 * powers[j]) for j in range(8))
    p8 += 20160 * prod(x)
    for i, j in combinations(range(8), 2):
        for k in (0, 2, 4, 6, 8):
            p8 += 4 * s.binomial(8, k) * x[i] ** (8 - k) * x[j] ** k
    q = sum(value * value for value in x)
    h8 = s.expand(p8 - 30 * q**4)

    contractions = {
        (a, b): sum(h[a][j] * h[b][j] for j in range(8))
        for a in range(1, 8)
        for b in range(a, 8)
        if a + b <= 8
    }
    contraction_words = [
        tuple(tuple(pair) for pair in word)
        for word in result["basis_after_H8"]
    ]
    basis = [h8] + [
        s.expand(prod(contractions[pair] for pair in word))
        for word in contraction_words
    ]
    need(len(coefficients) == len(basis) == 17, "completion basis changed")
    state = s.Poly(s.expand(sum(c * value for c, value in zip(coefficients, basis))), *variables)
    need(len(state.terms()) == result["full_PBW_terms"] == 1539, "W8 PBW support changed")
    need(all(sum((flat // 8 + 1) * exponent for flat, exponent in enumerate(monomial)) == 8
             for monomial, _coefficient in state.terms()),
         "W8 contains a PBW monomial of the wrong quotient level")
    return h, variables, state, result


def source_constant(n, component, d):
    return ((-1) ** (n - 1)) * d[component]


def level_one_test(h, state, d):
    """Small sign and primary-Ward diagnostic before the full expansion."""
    constants = {
        h[n][j]: s.Rational(source_constant(n, j, d))
        for n in range(1, 8)
        for j in range(8)
    }
    polynomial = state.as_expr()
    onepoint = s.factor(polynomial.subs(constants, simultaneous=True))
    vector = tuple(
        s.factor(s.diff(polynomial, h[1][j]).subs(constants, simultaneous=True))
        for j in range(8)
    )
    ward = s.factor(sum(d[j] * vector[j] for j in range(8)))
    norm = s.factor(sum(value * value for value in vector))
    need(ward == 8 * onepoint, "level-one primary Ward identity failed")
    need(norm > 0, "shifted W8 has no level-one response")
    need(norm >= s.factor(64 * onepoint**2 / sum(value * value for value in d)),
         "level-one Cauchy bound failed")
    return onepoint, vector, norm, ward


def merge_keys(left, right):
    """Merge two sorted occupation lists of mode-component indices."""
    return tuple(sorted(left + right))


def linear_power_expansion(n, component, exponent, d):
    """Expand one substituted oscillator power through total level eight."""
    constant = as_fraction(source_constant(n, component, d))
    options = [(0, None, constant)] + [
        (r, (r - 1) * 8 + component, Fraction(comb(r - 1, n - 1)))
        for r in range(n, MAX_LEVEL + 1)
    ]
    expansion = {(0, ()): Fraction(1)}
    for _ in range(exponent):
        updated = {}
        for (level, key), coefficient in expansion.items():
            for added_level, index, factor in options:
                new_level = level + added_level
                if new_level > MAX_LEVEL or factor == 0:
                    continue
                new_key = key if index is None else merge_keys(key, (index,))
                token = (new_level, new_key)
                updated[token] = updated.get(token, Fraction(0)) + coefficient * factor
        expansion = {key: value for key, value in updated.items() if value}
    return expansion


def shifted_levels(state, d):
    """Sparse exact coefficients v_0,...,v_8 of z^8 Y_d(W8,z)|d>."""
    power_cache = {}
    levels = [dict() for _ in range(MAX_LEVEL + 1)]
    for monomial, sympy_coefficient in state.terms():
        coefficient = as_fraction(sympy_coefficient)
        expansion = {(0, ()): coefficient}
        for flat_index, exponent in enumerate(monomial):
            if exponent == 0:
                continue
            n, component = divmod(flat_index, 8)
            n += 1
            cache_key = (n, component, exponent)
            factor_expansion = power_cache.get(cache_key)
            if factor_expansion is None:
                factor_expansion = linear_power_expansion(n, component, exponent, d)
                power_cache[cache_key] = factor_expansion
            updated = {}
            for (level_a, key_a), coefficient_a in expansion.items():
                for (level_b, key_b), coefficient_b in factor_expansion.items():
                    level = level_a + level_b
                    if level > MAX_LEVEL:
                        continue
                    key = merge_keys(key_a, key_b)
                    token = (level, key)
                    updated[token] = updated.get(token, Fraction(0)) + coefficient_a * coefficient_b
            expansion = {key: value for key, value in updated.items() if value}
        for (level, key), coefficient in expansion.items():
            levels[level][key] = levels[level].get(key, Fraction(0)) + coefficient
    return [{key: value for key, value in level.items() if value} for level in levels]


def fock_norm(polynomial):
    total = Fraction(0)
    for key, coefficient in polynomial.items():
        occupations = Counter(key)
        monomial_norm = 1
        for flat_index, exponent in occupations.items():
            mode = flat_index // 8 + 1
            monomial_norm *= mode**exponent * factorial(exponent)
        total += coefficient * coefficient * monomial_norm
    return total


def reconstruct_pole_coefficients(level_norms):
    """Solve a_l=sum_(r<=l) binom(l+r-1,2r-1) A_r."""
    coefficients = [level_norms[0]]
    for level in range(1, MAX_LEVEL + 1):
        remainder = level_norms[level]
        for r in range(1, level):
            remainder -= comb(level + r - 1, 2 * r - 1) * coefficients[r]
        coefficients.append(remainder)
    return coefficients


def standard_root_action_check(fold, d_quarters):
    reversal = s.zeros(8)
    for j in range(8):
        reversal[7 - j, j] = 1
    action = s.simplify(-fold * reversal * fold.inv())
    roots = set(standard_roots())
    need(s.simplify(action.T * action) == s.eye(8), "quarter transport is not orthogonal")
    need({tuple(action * s.Matrix(root)) for root in roots} == roots,
         "quarter transport does not preserve the actual E8 roots")
    need(s.simplify(action * s.Matrix(d_quarters["1/4"]) - s.Matrix(d_quarters["3/4"])) == s.zeros(8, 1),
         "quarter transport does not map the two actual source charges")
    return action


def main():
    h, variables, state, completion = build_state()
    old = load_old_source()
    lattice, orthogonal, _positive, fold = old.original_transport()

    # The W8 polynomial above uses canonical standard E8 coordinates.  The
    # old root chart is O*(old physical root), while the source charge is M*delta.
    roots = set(standard_roots())
    chart_roots = {
        tuple(s.simplify(value) for value in orthogonal * (s.Matrix(root) / s.sqrt(2)))
        for root in lattice["roots"]
    }
    need(chart_roots == roots, "old O chart does not reach the canonical root set")

    d_quarters = {}
    preliminary = {}
    for alpha in (s.Rational(1, 4), s.Rational(3, 4)):
        label = str(alpha)
        delta = s.Matrix([(j + alpha) / 8 - s.Rational(1, 2) for j in range(8)])
        d = tuple(s.Rational(value) for value in fold * delta)
        d_quarters[label] = d
        onepoint, vector, norm, ward = level_one_test(h, state, d)
        expected_onepoint = s.Rational(
            completion["original_quarter_source_cases"][label]["W8_onepoint_times_z8"]
        )
        need(onepoint == expected_onepoint, "old W8 one-point coefficient changed")
        preliminary[label] = {
            "d_standard": [str(value) for value in d],
            "d_squared": str(s.factor(sum(value * value for value in d))),
            "onepoint": str(onepoint),
            "level_1_vector": [str(value) for value in vector],
            "d_dot_level_1": str(ward),
            "level_1_norm_squared": str(norm),
            "cauchy_lower_bound": str(s.factor(64 * onepoint**2 / sum(value * value for value in d))),
            "level_1_excess_over_cauchy_bound": str(s.factor(
                norm - 64 * onepoint**2 / sum(value * value for value in d)
            )),
        }

    action = standard_root_action_check(fold, d_quarters)
    level_one_left = s.Matrix([s.Rational(value) for value in preliminary["1/4"]["level_1_vector"]])
    level_one_right = s.Matrix([s.Rational(value) for value in preliminary["3/4"]["level_1_vector"]])
    need(s.simplify(action * level_one_left - level_one_right) == s.zeros(8, 1),
         "level-one vectors do not obey the actual quarter Weyl transport")

    cases = {}
    for label, d in d_quarters.items():
        levels = shifted_levels(state, d)
        level_norms = [fock_norm(level) for level in levels]
        pole_coefficients = reconstruct_pole_coefficients(level_norms)
        need(level_norms[0] == as_fraction(s.Rational(preliminary[label]["onepoint"]))**2,
             "level-zero norm is not the one-point square")
        need(level_norms[1] == as_fraction(s.Rational(preliminary[label]["level_1_norm_squared"])),
             "sparse level-one result differs from the Ward precheck")
        need(pole_coefficients[8] == as_fraction(s.Rational(completion["norm_squared"])),
             "leading collision coefficient is not the W8 vacuum norm")
        # Independent stress-tensor Ward check.  For weight8 and c=8,
        # the T coefficient in W8 W8 is (2*8/8)*||W8||^2=2||W8||^2.
        # With <T(z)>=h_d/z^2 this fixes A7=2*h_d*||W8||^2.
        h_d = sum(value*value for value in d)/2
        stress_expected = as_fraction(2*h_d*s.Rational(completion["norm_squared"]))
        need(pole_coefficients[7] == stress_expected,
             "next collision coefficient fails the independent stress Ward identity")
        need(all(value >= 0 for value in level_norms), "a shifted Fock norm is negative")
        cases[label] = {
            **preliminary[label],
            "level_support_sizes": [len(level) for level in levels],
            "level_norms": [rational_text(value) for value in level_norms],
            "pole_basis_coefficients": [rational_text(value) for value in pole_coefficients],
            "independent_stress_Ward_A7": rational_text(stress_expected),
        }

    need(cases["1/4"]["level_norms"] == cases["3/4"]["level_norms"],
         "Weyl-related quarters have different level norms")
    need(cases["1/4"]["pole_basis_coefficients"] == cases["3/4"]["pole_basis_coefficients"],
         "Weyl-related quarters have different separated two-point functions")

    output = {
        "claim": "MATH.ORIGINAL-SEAM.W8-SOURCE-TWOPOINT.01",
        "status": "PASS_EXACT_IN_EXISTING_SHIFTED_CURRENT_SOURCE",
        "radial_variable": "q_rad=w/z=exp(-tau), |q_rad|<1",
        "radial_generator": "N=Lhat_0-85/256=(Lcal_0^alpha-1/32)/8=K_alpha on the folded neutral sector",
        "original_grading": "H0=Lcal_0^alpha-1/32=8*K_alpha=8*N; h_j(-r) has N-level r and original H0-grade 8r",
        "time_scope": "tau is quotient radial-cylinder separation; this algebraic conversion does not supply an absolute physical-time calibration",
        "pole_variable": "X=q_rad/(1-q_rad)^2",
        "formula": "z^8 w^8 <W8(z)W8(w)>_d=sum_(r=0)^8 A_r X^r",
        "connected_formula": "subtract A_0=<W8>_d^2",
        "reconstruction": "a_l=||v_l||^2=sum_(r=1)^l binom(l+r-1,2r-1)A_r for l>0",
        "completeness_reason": "bosonic exchange gives F(q)=F(1/q); after z^8w^8 the only finite pole is q=1 of order at most 16, so F is a degree-at-most-8 polynomial in X",
        "root_chart": "old O maps all 240 old physical roots to the canonical standard E8 roots",
        "quarter_relation": {
            "description": "d_3/4=g d_1/4 with g=-M R M^-1 preserving the actual 240-root set",
            "matrix": [[str(action[i, j]) for j in range(8)] for i in range(8)],
        },
        "vacuum_control": {
            "d": "0",
            "a_0": "0",
            "a_1": "0",
            "first_possible_level": 8,
            "level_8_norm_squared": completion["norm_squared"],
        },
        "cases": cases,
        "not_claimed": [
            "a new state, coupling, direction, or source model",
            "an identification of the two complete CAR modules",
            "a physical measurement or preparation",
            "a full P1 source derivation",
        ],
    }
    print(json.dumps(output, indent=2) + "\n", end="")


if __name__ == "__main__":
    main()
