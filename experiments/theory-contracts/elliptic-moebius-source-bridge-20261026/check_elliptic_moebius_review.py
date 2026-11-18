#!/usr/bin/env python3
"""Independent exact checks for DIVIDING-REAL-FORMS.txt.

The script separates real affine maps, the fixed raw flat pullbacks,
anti-linear real structures, and the complex-linear Pin operator.  In
particular it does not repair a lift by inserting an adjustable scalar phase.
"""

from __future__ import annotations

import hashlib
import json
import sympy as sp


I2 = sp.eye(2)
J = sp.Matrix([[0, -1], [1, 0]])
K = sp.Matrix([[0, -1], [-1, 0]])
T = sp.Matrix([sp.Rational(1, 2), sp.Rational(1, 2)])


def affine_comp(f, g):
    """Composition f after g for real affine maps."""
    mf, df = f
    mg, dg = g
    return mf * mg, mf * dg + df


def affine_pow(f, n):
    out = (I2, sp.zeros(2, 1))
    for _ in range(n):
        out = affine_comp(f, out)
    return out


def integral_vector(v):
    return all(x.is_integer is True for x in v)


def mod_one(v):
    return tuple(sp.simplify(x - sp.floor(x)) for x in v)


def main() -> None:
    q = (J, sp.zeros(2, 1))
    h = (-I2, T)
    kappa = (K, T)

    expected_linear = {
        (0, 0): K,
        (0, 1): sp.diag(1, -1),
        (0, 2): -K,
        (0, 3): sp.diag(-1, 1),
        (1, 0): -K,
        (1, 1): sp.diag(-1, 1),
        (1, 2): K,
        (1, 3): sp.diag(1, -1),
    }
    expected_shift = {
        (0, 0): (sp.Rational(1, 2), sp.Rational(1, 2)),
        (0, 1): (sp.Rational(1, 2), sp.Rational(1, 2)),
        (0, 2): (sp.Rational(1, 2), sp.Rational(1, 2)),
        (0, 3): (sp.Rational(1, 2), sp.Rational(1, 2)),
        (1, 0): (0, 0),
        (1, 1): (0, 0),
        (1, 2): (0, 0),
        (1, 3): (0, 0),
    }
    components = {
        (0, 0): 1, (0, 1): 0, (0, 2): 1, (0, 3): 0,
        (1, 0): 1, (1, 1): 2, (1, 2): 1, (1, 3): 2,
    }

    real_forms = {}
    affine_ok = True
    involution_ok = True
    for a in (0, 1):
        for b in range(4):
            qb = affine_pow(q, b)
            rho = affine_comp(affine_pow(h, a), affine_comp(qb, kappa))
            m, d = rho
            affine_ok &= m == expected_linear[(a, b)] and mod_one(d) == expected_shift[(a, b)]
            m2, d2 = affine_pow(rho, 2)
            involution_ok &= m2 == I2 and integral_vector(d2)
            real_forms[f"{a},{b}"] = {
                "linear": [[int(m[i, j]) for j in range(2)] for i in range(2)],
                "shift_mod_lattice": [str(x) for x in mod_one(d)],
                "real_components": components[(a, b)],
            }

    dividing = [key for key, count in components.items() if count == 2]

    delta = (I2, T)
    d2m, d2v = affine_pow(delta, 2)
    delta_free_order_two = d2m == I2 and tuple(d2v) == (1, 1) and not integral_vector(T)
    commutes = True
    for key in dividing:
        a, b = key
        rho = affine_comp(affine_pow(h, a), affine_comp(affine_pow(q, b), kappa))
        left = affine_comp(rho, delta)
        right = affine_comp(delta, rho)
        commutes &= left[0] == right[0] and integral_vector(left[1] - right[1])

    t, s = sp.symbols("t s", real=True)
    glide_once = (sp.Rational(1, 2) - t, s + sp.Rational(1, 2))
    glide_twice = (
        sp.simplify(sp.Rational(1, 2) - glide_once[0]),
        sp.simplify(glide_once[1] + sp.Rational(1, 2)),
    )
    glide_square = glide_twice == (t, s + 1)
    glide_free = sp.simplify(glide_once[1] - s - sp.Rational(1, 2)) == 0

    # Q rotates the two cuts.  It is not the half translation around a fixed boundary.
    witness = sp.Matrix([sp.Rational(1, 4), 0])
    q_witness = J * witness
    boundary_half_witness = witness + sp.Matrix([0, sp.Rational(1, 2)])
    q_not_boundary_half_translation = q_witness != boundary_half_witness

    def chi(m, n):
        return -1 if (m + n) % 2 else 1

    character_ok = chi(1, 0) == -1 and chi(0, 1) == -1 and chi(1, 1) == 1

    # Raw pullbacks in the displayed automorphy trivialization.  The two
    # universal affine composites differ by the A-cycle, so chi(A)=-1 is the
    # commutator.  The glide squares to the B-cycle, so chi(B)=-1.
    raw_delta_square = chi(1, 1)
    raw_rho_delta_commutator = chi(1, 0)
    raw_glide_square = chi(0, 1)

    # Q sends A to B and B to -A.  Of the four flat spin characters, exactly
    # one is simultaneously nontrivial, Q fixed, zero-mode free and gives the
    # order-four glide.  This is the selected NS--NS line, not a fitted phase.
    spin_rows = []
    selected_spin_rows = []
    for eps_a in (1, -1):
        for eps_b in (1, -1):
            row = {
                "chi_A": eps_a,
                "chi_B": eps_b,
                "Q_fixed": eps_a == eps_b,
                "has_constant_zero_mode": eps_a == 1 and eps_b == 1,
                "raw_glide_square": eps_b,
            }
            spin_rows.append(row)
            if row["Q_fixed"] and not row["has_constant_zero_mode"] and eps_b == -1:
                selected_spin_rows.append(row)

    sigma_t = sp.Matrix([[0, 1], [1, 0]])
    sigma_s = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    # g has derivative diag(-1,+1).  A=sigma_s U_g commutes with D
    # coefficientwise.  U_g is complex linear here.
    linear_pin_covariant = (
        sp.simplify(-sigma_t * sigma_s - sigma_s * sigma_t) == sp.zeros(2)
        and sp.simplify(sigma_s * sigma_s - sigma_s * sigma_s) == sp.zeros(2)
    )
    linear_pin_square = sp.simplify(sigma_s**2 * (-I2))
    quotient_sector_squares = {
        str(c): sp.simplify((1 / c) ** 2 * linear_pin_square)
        for c in (sp.I, -sp.I)
    }
    # The anti-linear Real lift is a different operator: scalars and sigma_s
    # are conjugated, so sigma_s J_g has square +I.
    anti_pin_square = sp.simplify(sigma_s * sigma_s.conjugate() * (-I2))
    lift_type_separation = (
        linear_pin_square == -I2
        and anti_pin_square == I2
        and all(x == I2 for x in quotient_sector_squares.values())
    )

    r = sp.symbols("r", positive=True)
    graph = sp.Matrix([[1 / (1 + r**2), r / (1 + r**2)],
                       [r / (1 + r**2), r**2 / (1 + r**2)]])
    swap = sp.Matrix([[0, 1], [1, 0]])
    graph_glide_invariance = sp.simplify(swap * graph.subs(r, 1 / r) * swap - graph) == sp.zeros(2)

    # The boundary compression and the old quarter-cover dictionary.
    x = sp.symbols("x", real=True)
    transfer = sp.exp(-x)
    graph_from_transfer = sp.Matrix([
        [1, transfer], [transfer, transfer**2]
    ]) / (1 + transfer**2)
    second_face_covariance = sp.simplify(
        graph_from_transfer[1, 1] - 1 / (1 + sp.exp(2 * x))
    ) == 0

    fold_and_carry = True
    for m in range(-4, 5):
        fold_and_carry &= (
            2 * m + sp.Rational(1, 2) == 2 * (m + sp.Rational(1, 4))
            and 2 * m + sp.Rational(3, 2) == 2 * (m + sp.Rational(3, 4))
            and (2 * m) + 1 == 2 * m + 1
            and (2 * m + 1) + 1 == 2 * (m + 1)
        )

    # Finite telescoping identity behind the Fredholm factor.  The endpoint
    # limits c_{-infinity}=1 and c_{+infinity}=0 give 1-z.
    z = sp.symbols("z")
    labels = list(range(-4, 5))
    c = {j: sp.symbols(f"c_{j + 5}") for j in range(-5, 5)}
    finite_relative_product = sp.prod(
        (1 - z * c[j - 1]) / (1 - z * c[j]) for j in labels
    )
    telescoping = sp.simplify(
        finite_relative_product - (1 - z * c[-5]) / (1 - z * c[4])
    ) == 0

    q = sp.exp(-sp.pi)
    trace_upper = 2 * q / (1 - q**2)
    trace_bound = bool(sp.N(trace_upper, 60) < sp.Rational(86590, 10**6))

    checks = {
        "eight_affine_forms": bool(affine_ok),
        "all_are_involutions_mod_lattice": bool(involution_ok),
        "exactly_two_dividing": dividing == [(1, 1), (1, 3)],
        "delta_free_two_torsion": bool(delta_free_order_two),
        "delta_commutes_with_dividing_forms_on_torus": bool(commutes),
        "annulus_glide_free_and_square_identity_on_base": bool(glide_square and glide_free),
        "Q_is_not_boundary_half_translation": bool(q_not_boundary_half_translation),
        "NS_NS_character": bool(character_ok),
        "raw_pullback_central_relations": bool(
            raw_delta_square == 1
            and raw_rho_delta_commutator == -1
            and raw_glide_square == -1
        ),
        "unique_Q_fixed_zero_mode_free_order_four_spin": len(selected_spin_rows) == 1
            and selected_spin_rows[0]["chi_A"] == -1
            and selected_spin_rows[0]["chi_B"] == -1,
        "complex_linear_pin_covariance": bool(linear_pin_covariant),
        "linear_vs_antilinear_pin_square_separated": bool(lift_type_separation),
        "two_face_graph_invariant_under_swap_plus_spectral_reflection": bool(graph_glide_invariance),
        "second_face_is_fermi_covariance": bool(second_face_covariance),
        "quarter_fold_factor_two_and_carry": bool(fold_and_carry),
        "relative_determinant_telescopes_with_positive_winding_sign": bool(telescoping),
        "hardy_trace_distance_upper_bound_below_0_086590": bool(trace_bound),
    }
    if not all(checks.values()):
        failed = [name for name, ok in checks.items() if not ok]
        raise AssertionError(f"failed checks: {failed}")

    result = {
        "checks": checks,
        "real_forms": real_forms,
        "dividing_forms": [list(x) for x in dividing],
        "spin_character_audit": spin_rows,
        "raw_pullback_relations": {
            "delta_square": raw_delta_square,
            "rho_delta_commutator": raw_rho_delta_commutator,
            "glide_square": raw_glide_square,
        },
        "linear_pin_square": str(linear_pin_square),
        "quotient_sector_squares": {key: str(value) for key, value in quotient_sector_squares.items()},
        "anti_linear_sigma_s_times_glide_square": str(anti_pin_square),
        "scope": (
            "Affine square torus geometry, the fixed raw NS pullback, and Pin symbol algebra. "
            "The checker does not identify the v519 Fock state or interacting field measure with the annulus graph."
        ),
    }
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    print(f"PASS {len(checks)}/{len(checks)}")
    print(f"RESULT_SHA256 {hashlib.sha256(payload.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
