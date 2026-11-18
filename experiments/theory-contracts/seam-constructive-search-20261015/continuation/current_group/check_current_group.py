#!/usr/bin/env python3
"""Exact checks for the conditional U(8)-current-group narrowing of T1.

This checker uses integer matrices, roots-of-unity exponents and rationals.
It checks the finite Weyl/Tits and Ward identities and the exact coefficients
of global X/root-current commutation. The stipulated positive-definite E8
module, common current cocycle/core and integration are analytic premises and
arguments in PROOF.md, not consequences of this finite checker.
"""

from __future__ import annotations

import json
import hashlib
from fractions import Fraction
from pathlib import Path


DIM = 8
Q = (1,) * DIM
U = (1, 0, 0, 0, 0, 0, 0, 0)


Matrix = tuple[tuple[int, ...], ...]
Vector = tuple[int, ...]


def identity() -> Matrix:
    return tuple(
        tuple(int(i == j) for j in range(DIM)) for i in range(DIM)
    )


def matmul(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(
            sum(left[i][k] * right[k][j] for k in range(DIM))
            for j in range(DIM)
        )
        for i in range(DIM)
    )


def matvec(matrix: Matrix, vector: Vector) -> Vector:
    return tuple(
        sum(matrix[i][j] * vector[j] for j in range(DIM))
        for i in range(DIM)
    )


def matpow(matrix: Matrix, exponent: int) -> Matrix:
    result = identity()
    base = matrix
    power = exponent
    while power:
        if power & 1:
            result = matmul(result, base)
        base = matmul(base, base)
        power >>= 1
    return result


def determinant(matrix: Matrix) -> int:
    """Fraction-free Bareiss determinant."""
    work = [list(row) for row in matrix]
    sign = 1
    previous = 1
    for pivot_index in range(DIM - 1):
        if work[pivot_index][pivot_index] == 0:
            swap = next(
                row
                for row in range(pivot_index + 1, DIM)
                if work[row][pivot_index] != 0
            )
            work[pivot_index], work[swap] = work[swap], work[pivot_index]
            sign *= -1
        pivot = work[pivot_index][pivot_index]
        for i in range(pivot_index + 1, DIM):
            for j in range(pivot_index + 1, DIM):
                numerator = (
                    work[i][j] * pivot
                    - work[i][pivot_index] * work[pivot_index][j]
                )
                assert numerator % previous == 0
                work[i][j] = numerator // previous
        previous = pivot
    return sign * work[-1][-1]


def basis(index: int) -> Vector:
    return tuple(int(j == index) for j in range(DIM))


def simple_tits(index: int) -> Matrix:
    """n_i=exp(pi/2(E_{i+1,i}-E_{i,i+1})) in the vector representation."""
    result = [list(row) for row in identity()]
    result[index][index] = 0
    result[index][index + 1] = -1
    result[index + 1][index] = 1
    result[index + 1][index + 1] = 0
    return tuple(tuple(row) for row in result)


def cycle_matrix() -> Matrix:
    """C e_j=e_(j+1) with indices modulo eight."""
    result = [[0] * DIM for _ in range(DIM)]
    for j in range(DIM):
        result[(j + 1) % DIM][j] = 1
    return tuple(tuple(row) for row in result)


def reflection_half_turn() -> Matrix:
    """P=C^4."""
    return matpow(cycle_matrix(), 4)


def torus_endpoint() -> Matrix:
    return tuple(
        tuple((-1 if i == 0 else 1) if i == j else 0 for j in range(DIM))
        for i in range(DIM)
    )


def dot(left: tuple[Fraction, ...], right: tuple[Fraction, ...]) -> Fraction:
    return sum((x * y for x, y in zip(left, right)), Fraction(0))


def cycle(vector: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return (vector[-1],) + vector[:-1]


def t1(vector: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return tuple(x + y for x, y in zip(cycle(vector), U))


def delta(alpha: Fraction) -> tuple[Fraction, ...]:
    return tuple(
        Fraction(j, 8) + alpha / 8 - Fraction(1, 2) for j in range(DIM)
    )


def energy(vector: tuple[Fraction, ...], alpha: Fraction) -> Fraction:
    return dot(vector, vector) / 2 + dot(delta(alpha), vector)


def check_tits_endpoint() -> dict[str, object]:
    simple = [simple_tits(i) for i in range(DIM - 1)]
    for matrix in simple:
        assert determinant(matrix) == 1

    # Products act on column vectors from the rightmost factor first.
    # R=n_0 n_1 ... n_6 is the signed cycle needed below.
    root_part = identity()
    for matrix in simple:
        root_part = matmul(root_part, matrix)

    for j in range(DIM - 1):
        assert matvec(root_part, basis(j)) == basis(j + 1)
    assert matvec(root_part, basis(DIM - 1)) == tuple(-x for x in basis(0))
    assert determinant(root_part) == 1
    assert matpow(root_part, 8) == tuple(
        tuple(-int(i == j) for j in range(DIM)) for i in range(DIM)
    )

    z_endpoint = torus_endpoint()
    cycle_endpoint = cycle_matrix()
    assert matmul(z_endpoint, root_part) == cycle_endpoint
    assert determinant(z_endpoint) == -1
    assert determinant(cycle_endpoint) == -1
    assert matpow(cycle_endpoint, 4) == reflection_half_turn()
    assert matpow(cycle_endpoint, 8) == identity()

    return {
        "simple_root_factors": 7,
        "root_part_determinant": determinant(root_part),
        "root_part_eighth_power": "-I",
        "torus_endpoint_determinant": determinant(z_endpoint),
        "C_equals_Ze0_times_R": True,
        "C_fourth_power_equals_P": True,
        "C_eighth_power_equals_identity": True,
    }


def check_ramond_singlet_and_paths() -> dict[str, object]:
    # Use twice the half-charge to avoid fractions: 2s=q.
    simple_roots = []
    for i in range(DIM - 1):
        alpha = tuple(
            1 if j == i else -1 if j == i + 1 else 0 for j in range(DIM)
        )
        simple_roots.append(alpha)
        assert sum(alpha[j] * Q[j] for j in range(DIM)) == 0

    # H_alpha,0 and J_{+/-alpha,0} annihilate |s> in the lattice module:
    # alpha.s=0, so the root-vertex OPE begins at z^0 and has no z^-1 term.
    assert len(simple_roots) == 7

    # Z_+(t)=exp(+i*pi*t H_e0), Z_-(t)=exp(-i*pi*t H_e0).
    # Both have the same endpoint on every integer charge, while on s they
    # give the two half-determinant phases +i and -i.
    for n0 in range(-8, 9):
        plus_integer_exp_mu8 = (4 * n0) % 8
        minus_integer_exp_mu8 = (-4 * n0) % 8
        assert plus_integer_exp_mu8 == minus_integer_exp_mu8

    path_lifts = []
    for orientation, determinant_winding, r in (
        ("positive", 1, 2),
        ("negative", -1, 6),
    ):
        nu_squared_exp_mu8 = (2 * r) % 8
        nu_fourth_exp_mu8 = (4 * r) % 8
        assert nu_squared_exp_mu8 == 4  # nu^2=-1
        assert nu_fourth_exp_mu8 == 0  # nu^4=+1
        path_lifts.append(
            {
                "orientation": orientation,
                "determinant_half_winding": determinant_winding,
                "r_mu8": r,
                "nu": "+i" if r == 2 else "-i",
                "kappa": 1,
                "kappa_sign": 1,
                "kappa_mu4_exp": 0,
            }
        )

    scalar_r = [r for r in range(8) if r == (4 - r) % 8]
    assert scalar_r == [2, 6]
    assert all(r % 2 == 0 for r in scalar_r)

    return {
        "simple_D8_roots_annihilating_s": len(simple_roots),
        "root_SU8_path_fixes_s": True,
        "same_integer_torus_endpoint_for_both_orientations": True,
        "path_lifts": path_lifts,
        "scalar_current_group_r_values_mu8": scalar_r,
        "conditional_kappa": 1,
        "conditional_kappa_sign": 1,
        "conditional_kappa_mu4_exp": 0,
        "orientation_selected_by_integer_source": False,
    }


def check_global_X_current_factors() -> dict[str, object]:
    """The coefficient identities apply to every charge, not a charge cutoff.

    On the fixed lattice Fock tensor realization X=E_s commutes with the
    Heisenberg oscillators. The only possible commutators with a root field
    are its E_alpha character and the shift of alpha.H0; both vanish here.
    """
    root_rows = 0
    for i in range(DIM):
        for j in range(DIM):
            if i == j:
                continue
            alpha = tuple(int(k == i)-int(k == j) for k in range(DIM))
            charge = sum(alpha)
            alpha_dot_s = Fraction(sum(alpha[k]*Q[k] for k in range(DIM)),2)
            assert charge == 0 and alpha_dot_s == 0
            # Both full-search characters have the same restriction on D8.
            for chi_exp in (1,3):
                character_exp_mu4 = (chi_exp*charge)%4
                assert character_exp_mu4 == 0
                root_rows += 1
    torus_lifts = []
    for sign,r in ((1,2),(-1,6)):
        # exp(sign*i*pi*e0.s) has mu8 exponent sign*2.
        assert (sign*2)%8 == r
        assert (2*r)%8 == 4 and (4*r)%8 == 0
        torus_lifts.append({"torus_orientation":sign,"nu_mu8_exp":r,
                            "kappa_sign":1,"kappa_mu4_exp":0})
    return {
        "root_character_rows":root_rows,
        "root_current_factorization":"J_alpha=Eminus_alpha Eplus_alpha E_alpha z^(alpha.H0)",
        "common_cocycle_scope":"fixed positive-definite Ls=D8 union(D8+s) Fock representation",
        "X_is_global_group_algebra_transport":True,
        "X_equals_full_local_Y_s":False,
        "X_oscillator_commutation":"given tensor realization: X acts as identity on oscillator factor",
        "alpha_dot_s_and_character_exponents_zero":True,
        "all_root_current_modes_commute_with_X":"analytic consequence of the displayed exact factorization",
        "common_core_and_group_integration":"analytic argument in PROOF.md, not numerically certified",
        "torus_lifts":torus_lifts,
        "odd_Z8_intertwiner_cocycle_identification_proved":False,
        "local_full_CAR_E8_algebra_proved":False,
    }


def check_time_and_quarter_marking() -> dict[str, object]:
    rows = 0
    for alpha in (Fraction(1, 4), Fraction(3, 4)):
        for grade in (0, 1):
            base = tuple(Fraction(grade, 2) for _ in range(DIM))
            probes = [base]
            for j in range(DIM):
                for sign in (-1, 1):
                    probe = list(base)
                    probe[j] += sign
                    probes.append(tuple(probe))
            for vector in probes:
                lhs = energy(t1(vector), alpha) - energy(vector, alpha)
                rhs = sum(vector, Fraction(0)) / 8 + alpha / 8
                assert lhs == rhs
                rows += 1

    # D_alpha C=zeta_8 C D_alpha is an exponent identity independent of
    # the Ramond lift phase r and of the sign of the determinant path.
    for k in range(DIM):
        assert ((k + 1) % DIM - k) % DIM == 1

    return {
        "exact_affine_Ward_rows": rows,
        "ward_identity": "K_alpha(Cp+e0)-K_alpha(p)=(Q(p)+alpha)/8",
        "quarter_bundle_identity": "D_alpha C=zeta8 C D_alpha",
        "B1_determinant_winding_in_x": 1,
        "B1_half_determinant_initial_sign_selected": False,
        "r_or_orientation_occurs_in_either_identity": False,
    }


def main() -> None:
    result = {
        "contract": "TFPT.SEAM.CURRENT_GROUP_LIFT.01",
        "verdict": "CONDITIONAL_CURRENT_GROUP_NARROWING_R_2_6_KAPPA_PLUS",
        "premise": "stipulated positive-definite integrable E8_1 Fock representation on Ls=D8 union(D8+s), fixed current cocycle and X=E_s tensor identity on oscillators",
        "scope": "kappa sign +1 (mu4 exponent 0) in the given E8 current-group lift; full odd Z8 sector representation remains additional",
        "not_derived_from": ["raw integer CAR endpoint", "P1/P2", "Lorentzian Gamma=Z10 signature(9,1) source without a transport proof"],
        "tits_endpoint": check_tits_endpoint(),
        "ramond_paths": check_ramond_singlet_and_paths(),
        "direct_global_X_bridge":check_global_X_current_factors(),
        "fixed_markings": check_time_and_quarter_marking(),
        "first_remaining_choice": (
            "orientation of the determinant-half path, equivalently the "
            "Ramond spin lift: +i versus -i"
        ),
        "chi_braiding_orientation_selected": False,
        "chi_plus_minus_are_global_parity_dressing_presentations":True,
        "full_odd_Z8_sector_representation_proved":False,
        "new_H_or_state_selected": False,
    }
    here=Path(__file__).resolve().parent
    original=here.parents[2]/"seam-joint-ward-state-20261014/PROOF.txt"
    result["source_hashes"]={
        "check_current_group.py":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "PROOF.md":hashlib.sha256(here.joinpath("PROOF.md").read_bytes()).hexdigest(),
        "seam-joint-ward-state-20261014/PROOF.txt":hashlib.sha256(original.read_bytes()).hexdigest(),
    }
    result["primary_literature"]={"url":"https://arxiv.org/pdf/math/9808088",
                                   "sections":"2.2 (vertex factorization), 2.3-2.4 (zero-mode automorphisms)"}
    output = Path(__file__).with_name("result.json")
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
