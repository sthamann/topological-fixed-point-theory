#!/usr/bin/env python3
"""Finite algebraic test for the three-cabled seam and its jet quotient.

The Hardy/Fredholm statement is proved analytically in the accompanying text.
This script only checks its finite polynomial truncations and the induced
c, Q, and beta characters on the explicit cokernel.
"""
import json
from pathlib import Path

import sympy as sp


SHEETS = 3
TRUNCATION_DEGREE = 6
WINDINGS = range(1, 9)
I = sp.I
ZETA8 = (1 + I) / sp.sqrt(2)


def basis_index(sheet, degree, degree_bound):
    return sheet * degree_bound + degree


def multiplication_matrix(n, degree=TRUNCATION_DEGREE):
    """Matrix of z^n: (C[z]_<degree)^3 -> (C[z]_<degree+n)^3."""
    rows = SHEETS * (degree + n)
    cols = SHEETS * degree
    matrix = sp.zeros(rows, cols)
    for sheet in range(SHEETS):
        for power in range(degree):
            matrix[basis_index(sheet, power + n, degree + n),
                   basis_index(sheet, power, degree)] = 1
    return matrix


def quotient_basis(n, degree=TRUNCATION_DEGREE):
    """Columns for the expected jet quotient z^j on each sheet, 0 <= j < n."""
    rows = SHEETS * (degree + n)
    quotient = sp.zeros(rows, SHEETS * n)
    column = 0
    for sheet in range(SHEETS):
        for power in range(n):
            quotient[basis_index(sheet, power, degree + n), column] = 1
            column += 1
    return quotient


def q_action(degree_bound, spin_phase=sp.S.One):
    """Q sends sheet a to -a and z^j to i^j z^j."""
    action = sp.zeros(SHEETS * degree_bound)
    for sheet in range(SHEETS):
        target_sheet = (-sheet) % SHEETS
        for power in range(degree_bound):
            action[basis_index(target_sheet, power, degree_bound),
                   basis_index(sheet, power, degree_bound)] = spin_phase * I**power
    return action


def c_action(degree_bound):
    """c sends sheet a to a+1 and fixes z."""
    action = sp.zeros(SHEETS * degree_bound)
    for sheet in range(SHEETS):
        target_sheet = (sheet + 1) % SHEETS
        for power in range(degree_bound):
            action[basis_index(target_sheet, power, degree_bound),
                   basis_index(sheet, power, degree_bound)] = 1
    return action


def exact_text(value):
    return str(sp.simplify(value))


w, eta = sp.symbols("w eta", nonzero=True)
x = w**3
z = eta / (1 - x)
z2_reduced = sp.cancel((1 - x**2) / (1 - x)**2)
z_q = sp.simplify((-I * eta / x) / (1 - 1 / x))

checks = {
    "curve_identity_z_squared": sp.simplify(z2_reduced - (1 + x) / (1 - x)) == 0,
    "source_Q_sends_z_to_i_z": sp.simplify(z_q - I * z) == 0,
}
cases = {}

for n in WINDINGS:
    degree = TRUNCATION_DEGREE
    mult = multiplication_matrix(n, degree)
    jets = quotient_basis(n, degree)
    combined = mult.row_join(jets)

    domain_q = q_action(degree)
    codomain_q = q_action(degree + n)
    domain_c = c_action(degree)
    codomain_c = c_action(degree + n)

    quotient_q_ref = q_action(n)
    quotient_q_spin = q_action(n, ZETA8)
    quotient_c = c_action(n)

    q_ref_det = sp.simplify(quotient_q_ref.det())
    q_spin_det = sp.simplify(quotient_q_spin.det())
    relative_q = sp.simplify(q_ref_det / q_spin_det)
    relative_beta = sp.simplify(relative_q**2)
    expected_q = sp.simplify(ZETA8**(-SHEETS * n))
    expected_beta = I**n

    case_checks = {
        "rank_is_full_domain": mult.rank() == SHEETS * degree,
        "kernel_is_zero": len(mult.nullspace()) == 0,
        "cokernel_dimension_is_3n": mult.rows - mult.rank() == SHEETS * n,
        "jet_columns_complete_the_image": combined.rank() == mult.rows,
        "jet_count_is_3n": jets.cols == SHEETS * n,
        "Q_covariance_of_z_to_n": codomain_q * mult == I**n * mult * domain_q,
        "c_equivariance": codomain_c * mult == mult * domain_c,
        "Q_conjugates_c_to_inverse": (
            quotient_q_ref * quotient_c * quotient_q_ref.inv() == quotient_c**2
        ),
        "reference_Q_fourth_power_is_identity": quotient_q_ref**4 == sp.eye(SHEETS * n),
        "spin_Q_fourth_power_is_minus_identity": (
            quotient_q_spin**4 == -sp.eye(SHEETS * n)
        ),
        "relative_Q_character": sp.simplify(relative_q - expected_q) == 0,
        "relative_beta_character": sp.simplify(relative_beta - expected_beta) == 0,
        "c_character_is_n_regular_representations": (
            quotient_c.trace() == 0 and (quotient_c**2).trace() == 0
        ),
        "c_determinant_is_one": quotient_c.det() == 1,
    }
    for name, passed in case_checks.items():
        checks[f"n_{n}_{name}"] = bool(passed)

    cases[str(n)] = {
        "matrix_shape": [mult.rows, mult.cols],
        "rank": mult.rank(),
        "kernel_dimension": len(mult.nullspace()),
        "cokernel_dimension": mult.rows - mult.rank(),
        "quotient_basis": [
            f"sheet_{sheet}:z^{power}"
            for sheet in range(SHEETS)
            for power in range(n)
        ],
        "Q_reference_determinant": exact_text(q_ref_det),
        "Q_spin_determinant": exact_text(q_spin_det),
        "relative_Q_character": exact_text(relative_q),
        "expected_relative_Q_character": exact_text(expected_q),
        "relative_beta_character": exact_text(relative_beta),
        "expected_relative_beta_character": exact_text(expected_beta),
        "c_character": {
            "identity": SHEETS * n,
            "c": int(quotient_c.trace()),
            "c_squared": int((quotient_c**2).trace()),
            "determinant": int(quotient_c.det()),
        },
    }

errors = [name for name, passed in checks.items() if not passed]
payload = {
    "scope": (
        "Finite polynomial truncation only. The Hardy/Fredholm and determinant-line "
        "statements are analytic results in ORIGINAL-BOUNDARY-NORM-ORIGIN.txt."
    ),
    "curve": "eta^2=1-w^6",
    "map": "z=eta/(1-w^3), z^2=(1+w^3)/(1-w^3)",
    "truncation_degree": TRUNCATION_DEGREE,
    "sheet_action": {
        "c": "(a,j)->(a+1,j)",
        "Q_reference": "(a,j)->(-a,j) with phase i^j",
        "Q_spin": "(a,j)->(-a,j) with phase zeta8*i^j",
    },
    "cases": cases,
    "checks": checks,
    "checks_passed": not errors,
    "errors": errors,
}

output = json.dumps(payload, indent=2, sort_keys=True) + "\n"
Path(__file__).with_name("three-cabled-seam.json").write_text(output)
print(output, end="")
raise SystemExit(0 if not errors else 1)
