#!/usr/bin/env python3
"""Exact mode-filter search at the named neutral TFPT source.

Only Python's standard library is used. All coefficients and exponents are
compared exactly. Run with --outputpath FILE to write JSON, or omit it for stdout.
"""
from __future__ import annotations

import argparse
import itertools
import json
from fractions import Fraction
from math import comb
from pathlib import Path
from typing import Any

HALF_EXPONENTS = tuple(Fraction(k, 2) for k in range(-6, 7))
SQ_COEFFICIENT_POWERS = tuple(range(-8, 9))  # c^2 = 2^j, hence c = 2^(j/2)


def fstr(q: Fraction) -> str:
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def two_power(j: int) -> Fraction:
    return Fraction(2**j) if j >= 0 else Fraction(1, 2**(-j))


def extended_factor_search(K: int = 8) -> dict[str, Any]:
    """Search c*n^a*(n-1)^b*(n+1)^d against K*n*(n-1)*(n+1)."""
    hits = []
    for a, b, d in itertools.product(HALF_EXPONENTS, repeat=3):
        for j in SQ_COEFFICIENT_POWERS:
            c2 = two_power(j)
            exponents = {"n": Fraction(1) + 2*a,
                         "n-1": Fraction(1) + 2*b,
                         "n+1": Fraction(1) + 2*d}
            # Squaring makes the coefficient condition rational and exact.
            if K*c2 == 1 and exponents == {
                    "n": Fraction(-1), "n-1": Fraction(0), "n+1": Fraction(0)}:
                hits.append({"a": fstr(a), "b": fstr(b), "d": fstr(d),
                             "c_squared": fstr(c2), "j": j,
                             "positive_c": f"2^({j}/2)"})
    return {
        "identity": "K*n*(n-1)*(n+1)*s_n^2 = 1/n for all integer n>=2",
        "search_box": {"a_b_d": "half-integers in [-3,3]",
                       "coefficient": "positive c with c^2=2^j, j=-8..8"},
        "factor_forcing": {"a": "-1", "b": "-1/2", "d": "-1/2",
                           "c_squared": f"1/{K}",
                           "required_j_if_power_of_two": -3 if K == 8 else None},
        "exact_hits": hits,
        "unique_rule": ("s_n=2^(-3/2)*n^-1*(n-1)^-1/2*(n+1)^-1/2"
                        if K == 8 and len(hits) == 1 else None),
    }


def pure_power_search(K: int = 8) -> dict[str, Any]:
    """Search D-only multipliers c*n^a; exact global identity is impossible."""
    hits = []
    tested = 0
    for a in HALF_EXPONENTS:
        for j in SQ_COEFFICIENT_POWERS:
            c2 = two_power(j)
            tested += 1
            # Exact factorization: K*c^2*n^(1+2a)*(n-1)*(n+1).
            # Neither a nor c changes the nonzero exponents of n-1,n+1.
            factor_exponents = {"n-1": Fraction(1), "n+1": Fraction(1)}
            factor_exponents_match = all(value == 0 for value in factor_exponents.values())
            if K*c2 > 0 and factor_exponents_match:
                hits.append((a, j))
    return {
        "family": "s_n=c*n^a (D-only power rule)",
        "search_box": {"a": "half-integers in [-3,3]",
                       "coefficient": "positive c with c^2=2^j, j=-8..8"},
        "tested_finite_box": tested,
        "exact_hits": hits,
        "class_verdict": "NONE, including beyond the finite exponent box",
        "exact_reason": (f"After squaring, K*c^2*n^(1+2a)*(n-1)*(n+1)=1/n. "
                         "The irreducible factors n-1 and n+1 have exponent 1 "
                         "on the left and exponent 0 on the right; varying a "
                         "or c cannot remove them."),
    }


def standard_primary_mode_norm(C: int, n: int) -> Fraction:
    # Coefficient of z^(-n-2) w^(n-2) in C/(z-w)^4.
    return Fraction(C * comb(n + 1, 3))


def d_minus_two(K: int = 8, modes: tuple[int, ...] = (2, 4, 8)) -> dict[str, Any]:
    rows = []
    for n in modes:
        # s_n = c*n^-2 with c^2=1/K, so K*n*(n^2-1)*s_n^2.
        actual = Fraction(n*n-1, n**3)
        target = Fraction(1, n)
        rows.append({"n": n, "filtered_mode": fstr(actual), "target": fstr(target),
                     "relative_error_signed": fstr((actual-target)/target),
                     "relative_error_abs": fstr(abs(actual-target)/target)})
    return {
        "filter": "s_n=c*n^-2, c^2=1/8 (leading amplitude normalized)",
        "identity": "1/n - 1/n^3",
        "relative_error_formula_against_1_over_n": "-1/n^2",
        "rows": rows,
    }


def composition_check() -> dict[str, Any]:
    """Check s(p+q)=s(p)s(q) for the selected extended filter at p=q=2."""
    c2 = Fraction(1, 8)
    s4_squared = c2 * Fraction(1, 4**2) * Fraction(1, 4**2 - 1)
    s2_squared = c2 * Fraction(1, 2**2) * Fraction(1, 2**2 - 1)
    product_squared = s2_squared * s2_squared
    return {
        "condition": "s(p+q)=s(p)*s(q) for an unchanged common primitive-preserving product map",
        "test_pair": {"p": 2, "q": 2},
        "positive_multiplier_squares": {
            "s_4_squared": fstr(s4_squared),
            "(s_2*s_2)_squared": fstr(product_squared),
        },
        "equal": s4_squared == product_squared,
        "verdict": "FAILS for this candidate; the exact unequal squares suffice because both multipliers are positive.",
        "scope": ("This excludes using the same scalar modewise filter unchanged as an algebra/product homomorphism "
                  "on the whole source-field mode tower. A boundary-to-field port need not be a homomorphism; "
                  "this check alone does not exclude a separately specified port map or coupling."),
    }


def local_polynomial_exclusion(K: int = 8) -> dict[str, str]:
    return {
        "claim": "No nonzero finite polynomial multiplier P(n) gives K*n*(n^2-1)*P(n)^2=1/n for all n>=2.",
        "reason": "If deg(P)=m>=0, the left side grows as a nonzero constant times n^(3+2m), while 1/n decays.",
    }


def main() -> dict[str, Any]:
    C = 48
    samples = {str(n): fstr(standard_primary_mode_norm(C, n)) for n in (2, 3, 4, 8)}
    assert standard_primary_mode_norm(C, 2) == 48
    assert standard_primary_mode_norm(C, 3) == 192
    assert standard_primary_mode_norm(C, 4) == 480
    assert standard_primary_mode_norm(C, 8) == 4032

    extended = extended_factor_search(8)
    assert extended["exact_hits"] == [{"a": "-1", "b": "-1/2", "d": "-1/2",
                                       "c_squared": "1/8", "j": -3,
                                       "positive_c": "2^(-3/2)"}]
    pure = pure_power_search(8)
    assert pure["exact_hits"] == []
    comp = composition_check()
    assert comp["equal"] is False
    d2 = d_minus_two()
    assert [r["relative_error_signed"] for r in d2["rows"]] == ["-1/4", "-1/16", "-1/64"]

    return {
        "status": "EXACT_CONDITIONAL_CASIMIR_FILTER; SOURCE_SELECTION_OPEN",
        "scope": "Bounded exact mode-filter search at the named neutral TFPT source; no physical-origin or TOE claim.",
        "source_facts": {
            "R": "4*T_A2 - 8*T_U1",
            "vacuum_two_point": "48/(z-w)^4",
            "weight": 2,
            "primary_conditions": "In the stated grade-2 Heisenberg vacuum block, Tr(A_R)=0 gives L_2|R>=0; L_1 and higher positive Virasoro modes vanish under the stated vacuum assumptions.",
            "selection_limit": "The source notes that neutrality and primarity do not select R uniquely; the 48-branch orthogonal projection selects this representative.",
        },
        "mode_norm": {
            "convention": "R(z)=sum_m R_m*z^(-m-2), radial vacuum, no half-mode normalization.",
            "expansion": "48/(z-w)^4 = 48*z^-4*sum_{k>=0} binom(k+3,3)*(w/z)^k for |z|>|w|",
            "coefficient_match": "<0|R_n R_-n|0> = 48*binom(n+1,3) = 8*n*(n^2-1), n>=2",
            "exact_samples": samples,
            "n1": "The factor n^2-1 vanishes at n=1; the active tower starts at n=2.",
        },
        "pure_power_class": pure,
        "casimir_shift_class": extended,
        "D_minus_2": d2,
        "local_polynomial_class": local_polynomial_exclusion(8),
        "source_field_composition": comp,
        "original_mark_green_scope": {
            "source": "The original even mark-space orbit Green matrix is -c3*log(2)*C on four labels; compression gives 4*P0 and determinant 1-4*epsilon*c3*log(2).",
            "missing_mode_map": "No map is specified from these four mark labels to the R Laurent modes n>=2 or to the Casimir factors n, n-1, n+1.",
            "vertex_selection": "The original source_flavor_correlator says the boundary-to-R-vertex map/reference contraction is missing. Its conditional V=c3^2|R><p0| matches a target after that port is chosen; matching does not select it from Mark-Green data.",
            "verdict": "The inverse-Casimir rule is an exact conditional mode filter, not a filter selected by the original logarithmic Mark-Green source or P1.",
        },
        "prior_graph_check": {
            "code_graph": "No nodes for the two named tfpt_explorer files; the provided original paths were read directly.",
            "theory_graph": "OK, current generated snapshot; adjacent source/boundary contracts in the primer are PARTIAL and leave physical vertex/source selection open.",
            "phrase_search": "No direct hits for the exact neutral-source/Casimir-filter phrases.",
        },
    }


def cli() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outputpath", type=Path, help="write JSON results here; stdout if omitted")
    args = parser.parse_args()
    payload = json.dumps(main(), indent=2, ensure_ascii=False, sort_keys=True) + "\n"
    if args.outputpath:
        args.outputpath.parent.mkdir(parents=True, exist_ok=True)
        args.outputpath.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")


if __name__ == "__main__":
    cli()
