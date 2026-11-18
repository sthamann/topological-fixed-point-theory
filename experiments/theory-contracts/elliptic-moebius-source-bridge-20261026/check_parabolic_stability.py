#!/usr/bin/env python3
"""Exact finite stability check for the four rotated A0* flags.

This checks the algebraic parabolic object obtained by putting the A0* flags
in the affine splitting frame of O(-2)+O(-1)^2.  It separately checks whether
the diagonal infinity shear used to motivate that frame actually makes the
A0* connection form regular.  No assertion controls program behaviour, so
python and python -OO produce the same JSON.
"""

from __future__ import annotations

import argparse
import json
from itertools import product
from pathlib import Path

import sympy as sp


I = sp.I
SQ2 = sp.sqrt(2)
SQ5 = sp.sqrt(5)
SQ10 = sp.sqrt(10)
U = sp.diag(1, I, -I)
POINTS = (sp.Integer(1), I, sp.Integer(-1), -I)
A0 = sp.Matrix([
    [sp.Rational(1, 2), SQ2 / 6, 0],
    [SQ2 / 6, sp.Rational(1, 4), SQ5 / 12],
    [0, SQ5 / 12, sp.Rational(1, 4)],
])
V0 = sp.Matrix([SQ10 / 5, -3 * SQ5 / 5, 1])
V1 = sp.Matrix([-SQ10 / 5, SQ5 / 5, 1])
V2 = sp.Matrix([SQ10, SQ5, 1])
N = sp.diag(2, 1, 1)
M0 = sp.Matrix([
    [0, -(1 + I) / 2, (1 - I) / 2],
    [-(1 + I) / 2, -I / 2, sp.Rational(-1, 2)],
    [(1 - I) / 2, sp.Rational(-1, 2), I / 2],
])


def zero(x: sp.Expr) -> bool:
    return sp.simplify(x) == 0


def matrix_zero(m: sp.Matrix) -> bool:
    return all(zero(x) for x in m)


def enumerate_group(generators):
    def freeze(matrix):
        return tuple(sp.expand(x) for x in matrix)

    elements = {freeze(sp.eye(3)): sp.eye(3)}
    frontier = [sp.eye(3)]
    while frontier:
        following = []
        for element in frontier:
            for generator in generators:
                candidate = sp.expand(generator * element)
                key = freeze(candidate)
                if key not in elements:
                    elements[key] = candidate
                    following.append(candidate)
        frontier = following
    return list(elements.values())


def character_distribution(elements):
    distribution = {}
    for element in elements:
        character = str(sp.simplify(element.trace()))
        distribution[character] = distribution.get(character, 0) + 1
    return dict(sorted(distribution.items()))


def linear_system(equations, variables):
    matrix, _ = sp.linear_eq_to_matrix(equations, variables)
    return matrix, matrix.nullspace()


def polynomial_nonzero_on_space(nullspace, polynomial, variables) -> bool:
    if not nullspace:
        return False
    parameters = sp.symbols("t0:%d" % len(nullspace))
    vector = sum(
        (parameter * basis for parameter, basis in zip(parameters, nullspace)),
        sp.zeros(len(variables), 1),
    )
    restricted = sp.expand(polynomial.subs(dict(zip(variables, vector))))
    return restricted != 0


def line_pattern_equations(pattern, variables):
    a, b, c, d, e = variables
    equations = []
    for k, (point, threshold) in enumerate(zip(POINTS, pattern)):
        value = sp.Matrix([a, b * point + c, d * point + e])
        middle = U**k * V1
        high = U**k * V2
        if threshold >= 1:
            equations.append(sp.det(sp.Matrix.hstack(middle, high, value)))
        if threshold >= 2:
            equations.extend([
                value[0] * high[1] - value[1] * high[0],
                value[0] * high[2] - value[2] * high[0],
            ])
    return equations


def saturated_degree_minus_two_exists(nullspace, variables) -> bool:
    if not nullspace:
        return False
    if any(not zero(basis[0]) for basis in nullspace):
        return True
    _, b, c, d, e = variables
    return polynomial_nonzero_on_space(nullspace, b * e - c * d, variables)


def quotient_pattern_equations(pattern, variables):
    a, b, c, d = variables
    equations = []
    for k, (point, upper_weight) in enumerate(zip(POINTS, pattern)):
        functional = sp.Matrix([[a * point + b, c, d]])
        middle = U**k * V1
        high = U**k * V2
        if upper_weight <= 1:
            equations.append((functional * high)[0])
        if upper_weight <= 0:
            equations.append((functional * middle)[0])
    return equations


def surjective_degree_minus_one_exists(nullspace) -> bool:
    return any(not zero(basis[2]) or not zero(basis[3]) for basis in nullspace)


def actual_line_weights(a, b, c, d, e):
    weights = []
    normal0 = sp.Matrix([[-SQ10, 3 * SQ5, -5]])
    for k, point in enumerate(POINTS):
        value = sp.Matrix([a, b * point + c, d * point + e])
        high = U**k * V2
        proportional = (
            zero(value[0] * high[1] - value[1] * high[0])
            and zero(value[0] * high[2] - value[2] * high[0])
        )
        in_middle = zero((normal0 * U ** (-k) * value)[0])
        weights.append(2 if proportional else (1 if in_middle else 0))
    return tuple(weights)


def actual_quotient_weights(a, b, c, d):
    weights = []
    for k, point in enumerate(POINTS):
        functional = sp.Matrix([[a * point + b, c, d]])
        middle = U**k * V1
        high = U**k * V2
        weights.append(
            2 if not zero((functional * high)[0])
            else (1 if not zero((functional * middle)[0]) else 0)
        )
    return tuple(weights)


def run_check():
    checks = []
    failures = []

    def check(check_id, condition, detail):
        passed = bool(condition)
        checks.append({"id": check_id, "passed": passed, "detail": detail})
        if not passed:
            failures.append(check_id)

    lam = sp.symbols("lambda")
    check(
        "A0_SPECTRUM",
        zero(A0.charpoly(lam).as_expr() - lam * (lam - sp.Rational(1, 3))
             * (lam - sp.Rational(2, 3))),
        "char(A0*) = lambda(lambda-1/3)(lambda-2/3) exactly",
    )
    check(
        "A0_EIGENVECTORS",
        matrix_zero(A0 * V0)
        and matrix_zero((A0 - sp.eye(3) / 3) * V1)
        and matrix_zero((A0 - 2 * sp.eye(3) / 3) * V2),
        "V0,V1,V2 are exact eigenvectors for weights 0,1/3,2/3",
    )

    residues = [sp.simplify(U**k * A0 * (U**k).inv()) for k in range(4)]
    check(
        "C4_FLAGS",
        all(matrix_zero(residues[(k + 1) % 4] - U * residues[k] * U.inv())
            for k in range(4)),
        "A_(k+1)=U A_k U^-1; this proves C4 covariance only",
    )
    check(
        "RESIDUE_SUM",
        matrix_zero(sum(residues, sp.zeros(3)) - N),
        "sum A_k = diag(2,1,1)",
    )

    peripheral_monodromies = [sp.expand(U**k * M0 * U**(-k)) for k in range(4)]
    peripheral_group = enumerate_group(peripheral_monodromies)
    covariance_group = enumerate_group([M0, U])
    peripheral_characters = character_distribution(peripheral_group)
    covariance_characters = character_distribution(covariance_group)
    check(
        "FOUR_PUNCTURE_MONODROMY_IS_A4",
        len(peripheral_group) == 12
        and peripheral_characters == {"-1": 3, "0": 8, "3": 1}
        and sum(sp.simplify(g.trace() * sp.conjugate(g.trace()))
                for g in peripheral_group) == 12,
        "the four M_k generate irreducible A4 of order 12, not S4",
    )
    check(
        "U_EXTENDS_A4_TO_S4",
        len(covariance_group) == 24
        and covariance_characters == {"-1": 9, "0": 8, "1": 6, "3": 1},
        "adjoining the C4 covariance U extends A4 to the sign-twisted standard S4 representation",
    )

    line_variables = sp.symbols("a b c d e")
    line_feasible = []
    high_line_patterns = 0
    high_line_full_rank = 0
    for pattern in product(range(3), repeat=4):
        matrix, nullspace = linear_system(
            line_pattern_equations(pattern, line_variables), line_variables
        )
        if saturated_degree_minus_two_exists(nullspace, line_variables):
            line_feasible.append(pattern)
        if sum(pattern) >= 5:
            high_line_patterns += 1
            if matrix.rank() == 5:
                high_line_full_rank += 1
    line_max = max(sum(pattern) for pattern in line_feasible)
    line_witness = actual_line_weights(3 * SQ2 / 2, 1, 0, 0, 0)
    check(
        "RANK1_DEG_MINUS2_81_PATTERNS",
        len(tuple(product(range(3), repeat=4))) == 81
        and high_line_patterns == 31
        and high_line_full_rank == 31
        and line_max == 4
        and line_witness == (1, 1, 1, 1),
        "81 threshold patterns; all 31 with sum>=5 have rank 5; max sum=4",
    )

    y, z = sp.symbols("y z")
    degree_minus_one_feasible = []
    for pattern in product(range(3), repeat=4):
        equations = []
        for k, threshold in enumerate(pattern):
            value = sp.Matrix([0, y, z])
            middle = U**k * V1
            high = U**k * V2
            if threshold >= 1:
                equations.append(sp.det(sp.Matrix.hstack(middle, high, value)))
            if threshold >= 2:
                equations.extend([
                    value[0] * high[1] - value[1] * high[0],
                    value[0] * high[2] - value[2] * high[0],
                ])
        _, nullspace = linear_system(equations, (y, z))
        if nullspace:
            degree_minus_one_feasible.append(pattern)
    degree_minus_one_max = max(sum(pattern) for pattern in degree_minus_one_feasible)
    check(
        "RANK1_DEG_MINUS1",
        degree_minus_one_max == 2,
        "a degree -1 line lies in span(e2,e3) and has total weight at most 2/3",
    )

    quotient_variables = sp.symbols("qa qb qc qd")
    quotient_feasible = []
    low_quotient_patterns = 0
    low_quotient_full_rank = 0
    for pattern in product(range(3), repeat=4):
        matrix, nullspace = linear_system(
            quotient_pattern_equations(pattern, quotient_variables),
            quotient_variables,
        )
        if surjective_degree_minus_one_exists(nullspace):
            quotient_feasible.append(pattern)
        if sum(pattern) <= 3:
            low_quotient_patterns += 1
            if matrix.rank() == 4:
                low_quotient_full_rank += 1
    quotient_min = min(sum(pattern) for pattern in quotient_feasible)
    quotient_witness = actual_quotient_weights(-SQ2 / 2, 0, 1, 0)
    check(
        "RANK2_QUOTIENT_DEG_MINUS1_81_PATTERNS",
        len(tuple(product(range(3), repeat=4))) == 81
        and low_quotient_patterns == 31
        and low_quotient_full_rank == 31
        and quotient_min == 4
        and quotient_witness == (1, 1, 1, 1),
        "81 upper-weight patterns; all 31 with sum<=3 have rank 4; min sum=4",
    )

    q0_variables = sp.symbols("A B C D E F G")
    aa, bb, cc, dd, ee, ff, gg = q0_variables
    q0_equations = []
    for k, point in enumerate(POINTS):
        functional = sp.Matrix([[
            aa * point**2 + bb * point + cc,
            dd * point + ee,
            ff * point + gg,
        ]])
        q0_equations.extend([
            (functional * (U**k * V1))[0],
            (functional * (U**k * V2))[0],
        ])
    q0_matrix, q0_nullspace = linear_system(q0_equations, q0_variables)
    check(
        "RANK2_QUOTIENT_DEG_ZERO",
        q0_matrix.rank() == 7 and not q0_nullspace,
        "the eight all-zero-weight equations have rank 7 in seven variables",
    )

    rank1_bounds = {
        "degree_-1": "pardeg <= -1 + 2/3 = -1/3",
        "degree_-2": "pardeg <= -2 + 4/3 = -2/3",
        "degree_<=-3": "pardeg <= -3 + 8/3 = -1/3",
    }
    rank2_bounds = {
        "quotient_degree_-2": "pardeg(Q) = -2 + 8/3 = 2/3",
        "quotient_degree_-1": "pardeg(Q) >= -1 + 4/3 = 1/3",
        "quotient_degree_0": "pardeg(Q) >= 1/3",
        "quotient_degree_>0": "pardeg(Q) > 0",
    }
    check(
        "ALGEBRAIC_PARABOLIC_STABILITY",
        not failures,
        "all saturated rank-1 subbundles have pardeg<0 and all rank-2 quotients have pardeg>0",
    )

    # Infinity-patch audit.  Gamma is the connection form with finite residues +A_k.
    # For O(-2)+O(-1)^2 the diagonal frame transition is g=diag(w^2,w,w).
    # A genuine infinity-patch change leaves all finite fibres untouched, but it
    # must also make Gamma regular.  The exact transformed form retains a
    # resonant E12/w pole, so the raw A0* affine flags are not yet proved to be
    # the flags of the claimed BG-split Deligne extension.
    w = sp.symbols("w")
    gamma_affine = -sum(
        (residues[k] / (w * (1 - POINTS[k] * w)) for k in range(4)),
        sp.zeros(3),
    )
    g = sp.diag(w**2, w, w)
    gamma_infinity = sp.simplify(g.inv() * gamma_affine * g + N / w)
    infinity_residue = gamma_infinity.applyfunc(
        lambda entry: sp.simplify(sp.limit(w * entry, w, 0))
    )
    expected_resonance = sp.zeros(3)
    expected_resonance[0, 1] = -2 * SQ2 / 3
    pure_shear_holomorphic = matrix_zero(infinity_residue)
    check(
        "INFINITY_DIAGONAL_SHEAR_AUDIT",
        not pure_shear_holomorphic and matrix_zero(infinity_residue - expected_resonance),
        "g=diag(w^2,w,w) leaves the exact residue -(2sqrt(2)/3)E12; pure BG shear is insufficient",
    )

    # The previous algebraic stability checks remain valid for the declared
    # split bundle plus these flags.  What fails is only the unproved claim that
    # the v115 ODE supplies this same affine BG frame without a further
    # Levelt/Birkhoff gauge.
    overall_pass = not failures
    return {
        "verdict": "PASS_SCOPED" if overall_pass else "FAIL",
        "scope": {
            "proved": "stability of the declared split parabolic object with the four raw rotated A0* flags",
            "not_proved": "identity of that affine flag object with the Deligne extension of the v115 ODE",
            "symmetry": "C4 covariance exact; no D4 reflection linearization inferred",
            "v117_relation": "the four peripheral monodromies generate irreducible A4; U extends it to S4 and supplies the covariance data; A0* identification remains numerical",
        },
        "checks": checks,
        "failures": failures,
        "rank1": {
            "degree_minus2_patterns": 81,
            "degree_minus2_feasible_patterns": len(line_feasible),
            "sum_at_least_5_patterns": high_line_patterns,
            "sum_at_least_5_full_rank": high_line_full_rank,
            "maximum_weight_units": line_max,
            "degree_minus1_maximum_weight_units": degree_minus_one_max,
            "bounds": rank1_bounds,
        },
        "rank2": {
            "degree_minus1_quotient_patterns": 81,
            "degree_minus1_feasible_patterns": len(quotient_feasible),
            "sum_at_most_3_patterns": low_quotient_patterns,
            "sum_at_most_3_full_rank": low_quotient_full_rank,
            "minimum_weight_units": quotient_min,
            "degree_zero_all_weight_zero_matrix_rank": q0_matrix.rank(),
            "bounds": rank2_bounds,
        },
        "infinity_patch": {
            "affine_residue_sum": [2, 1, 1],
            "negative_BG_transition": "g(w)=diag(w^2,w,w)",
            "finite_flags_unchanged_if_transition_valid": True,
            "pure_diagonal_shear_holomorphic": pure_shear_holomorphic,
            "remaining_residue": [[sp.sstr(infinity_residue[i, j]) for j in range(3)] for i in range(3)],
            "consequence": "a further Levelt/Birkhoff lattice choice is required and may change the affine splitting frame",
        },
        "monodromy": {
            "four_puncture_group": "A4",
            "four_puncture_group_order": len(peripheral_group),
            "four_puncture_character_distribution": peripheral_characters,
            "covariance_extension_group": "S4",
            "covariance_extension_group_order": len(covariance_group),
            "covariance_character_distribution": covariance_characters,
            "kernel_cover": "degree 12 and genus 5 for the compactified four-puncture A4 kernel cover",
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("stability-results.json"),
    )
    args = parser.parse_args()
    result = run_check()
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    args.output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    return 0 if result["verdict"] == "PASS_SCOPED" else 1


if __name__ == "__main__":
    raise SystemExit(main())
