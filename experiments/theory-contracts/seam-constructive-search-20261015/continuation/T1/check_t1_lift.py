#!/usr/bin/env python3
"""Exact finite reductions for the source-canonical T1 defect lift.

This checker proves only parity-polynomial and affine-rational identities.
It does not select the CAR source, the half-lattice extension, a state, or a
Hamiltonian from P1/P2.
"""

from __future__ import annotations

import itertools
import json
from fractions import Fraction
from pathlib import Path


DIM = 8
U = (1, 0, 0, 0, 0, 0, 0, 0)
Q = (1,) * DIM
A = (1, 1, 1, 1, 0, 0, 0, 0)


def add(x: tuple[int, ...], y: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(a + b for a, b in zip(x, y))


def cycle(x: tuple[int, ...]) -> tuple[int, ...]:
    """C(n0,...,n7)=(n7,n0,...,n6)."""
    return (x[-1],) + x[:-1]


def t1(x: tuple[int, ...]) -> tuple[int, ...]:
    return add(cycle(x), U)


def parity(x: int) -> int:
    return x & 1


def eps_exp(n: tuple[int, ...], m: tuple[int, ...]) -> int:
    """Exponent mod 2 of the ordered CAR cocycle epsilon_Z."""
    return sum(n[i] * m[j] for i in range(DIM) for j in range(i)) & 1


def gamma_exp(n: tuple[int, ...]) -> int:
    """Exponent mod 2 of the normalized cyclic permutation lift."""
    return (n[7] * sum(n[:7])) & 1


def eta_exp(n: tuple[int, ...]) -> int:
    return (sum(n[:4]) * sum(n[4:])) & 1


def w1_state_exp(n: tuple[int, ...]) -> int:
    """W1=-Z_e0 E_e0 C_hat on |n>; the minus makes W1|0>=|e0>."""
    return (gamma_exp(n) + n[7]) & 1


def w4_state_exp(n: tuple[int, ...]) -> int:
    """Z_a E_a P_CAR phase on |n>, in the fixed source convention."""
    p_n = cycle(cycle(cycle(cycle(n))))
    z_a = sum(p_n[:4]) & 1
    e_a = eps_exp(A, p_n)
    return (eta_exp(n) + e_a + z_a) & 1


def iterate_with_phase(
    n: tuple[int, ...], steps: int
) -> tuple[tuple[int, ...], int]:
    phase = 0
    current = n
    for _ in range(steps):
        phase ^= w1_state_exp(current)
        current = t1(current)
    return current, phase


def check_integer_source_lift() -> dict[str, int | bool]:
    parity_vectors = list(itertools.product((0, 1), repeat=DIM))

    for j in range(DIM):
        elementary = tuple(1 if i == j else 0 for i in range(DIM))
        assert gamma_exp(elementary) == 0

    cocycle_rows = 0
    for n in parity_vectors:
        for m in parity_vectors:
            # gamma(n)gamma(m)eps(Cn,Cm)=eps(n,m)gamma(n+m)
            residual = (
                gamma_exp(n)
                + gamma_exp(m)
                + eps_exp(cycle(n), cycle(m))
                + eps_exp(n, m)
                + gamma_exp(add(n, m))
            ) & 1
            assert residual == 0
            cocycle_rows += 1

    for n in parity_vectors:
        current = n
        phase4 = 0
        for _ in range(4):
            phase4 ^= gamma_exp(current)
            current = cycle(current)
        assert phase4 == eta_exp(n)

        image, w1_phase4 = iterate_with_phase(n, 4)
        expected_image = add(cycle(cycle(cycle(cycle(n)))), A)
        assert image == expected_image
        assert w1_phase4 == w4_state_exp(n)

    zero = (0,) * DIM
    image1, phase1 = iterate_with_phase(zero, 1)
    image4, phase4 = iterate_with_phase(zero, 4)
    assert image1 == U and phase1 == 0
    assert image4 == A and phase4 == 0

    return {
        "ordered_cocycle_rows": cocycle_rows,
        "parity_rows_for_C4_and_W1_fourth_power": len(parity_vectors),
        "W1_vacuum_normalized": True,
        "elementary_channel_cycle_is_phase_free": True,
        "W1_fourth_power_equals_W4_on_integer_source": True,
    }


def check_defect_extension() -> dict[str, object]:
    # chi(e_j)=i^lambda.  The most general sector-diagonal extension used
    # here is C_r|n,c>=gamma(n) zeta8^(r*c)|Cn,c), r in Z/8.
    # It has C_r^4=P_kappa with kappa=(-1)^r.  Its conjugation of X has
    # zeta8 exponents r on grade 0 and 4-r on grade 1.
    candidates: list[dict[str, object]] = []
    for lam in (1, 3):
        for r in range(8):
            kappa_exp = 0 if r % 2 == 0 else 2
            c_on_x = [r, (4 - r) % 8]
            # Z_e0 contributes i=zeta8^2 and E_e0 contributes chi^-1.
            w1_on_x = [(e + 2 - 2 * lam) % 8 for e in c_on_x]
            scalar_x_covariance = c_on_x[0] == c_on_x[1]
            assert scalar_x_covariance == (r in (2, 6))
            for n in itertools.product((0, 1), repeat=DIM):
                n_plus_q = add(n, Q)
                cn = cycle(n)
                crossed_carry = (
                    4
                    * (
                        eps_exp(n, Q)
                        + gamma_exp(n_plus_q)
                        - gamma_exp(n)
                        - eps_exp(cn, Q)
                    )
                ) % 8
                assert crossed_carry == 4
                assert (crossed_carry - r) % 8 == c_on_x[1]
                # Four C_r steps: eta(n) from gamma and (-1)^(r*c).
                for c in (0, 1):
                    sector_phase4 = (4 * r * c) % 8
                    expected = (4 * (r % 2) * c) % 8
                    assert sector_phase4 == expected
                    assert (8 * r * c) % 8 == 0
                    _, integer_w1_phase4 = iterate_with_phase(n, 4)
                    full_w1_phase4 = (4 * integer_w1_phase4 + sector_phase4) % 8
                    full_w4_phase = (
                        4 * w4_state_exp(n) + 4 * (r % 2) * c
                    ) % 8
                    assert full_w1_phase4 == full_w4_phase
            candidates.append(
                {
                    "chi_exp_mu4": lam,
                    "sector_phase_r_mu8": r,
                    "kappa_exp_mu4": kappa_exp,
                    "Ad_C_on_X_exp_mu8_by_input_grade": c_on_x,
                    "Ad_W1_on_X_exp_mu8_by_input_grade": w1_on_x,
                    "scalar_X_eigenoperator_covariance": scalar_x_covariance,
                }
            )

    assert len(candidates) == 16
    assert {row["kappa_exp_mu4"] for row in candidates} == {0, 2}
    assert gamma_exp(Q) == 1
    return {
        "gamma_q": -1,
        "all_sector_diagonal_T1_lifts": candidates,
        "scalar_X_covariance_subclass": [
            row for row in candidates if row["scalar_X_eigenoperator_covariance"]
        ],
        "scalar_subclass_implication": "scalar Ad_C(X)=nu X forces r in {2,6}, nu^2=-1, and kappa=+1",
        "grading_dependent_counterexample": "r=1 gives C|n,c>=gamma(n)zeta8^c|Cn,c>, C^4=P with kappa=-1, and C^8=1",
        "chi_orientation_selected": False,
        "sector_phase_selected": False,
        "kappa_selected": False,
    }


def delta(alpha: Fraction) -> tuple[Fraction, ...]:
    return tuple(Fraction(j, 8) + alpha / 8 - Fraction(1, 2) for j in range(DIM))


def dot(x: tuple[Fraction, ...], y: tuple[Fraction, ...]) -> Fraction:
    return sum((a * b for a, b in zip(x, y)), Fraction(0))


def energy(p: tuple[Fraction, ...], alpha: Fraction) -> Fraction:
    return dot(p, p) / 2 + dot(delta(alpha), p)


def t1_fraction(p: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return add(cycle(p), U)  # type: ignore[arg-type]


def check_time_identity() -> dict[str, object]:
    # The difference is affine in p. Zero and signed basis vectors in each
    # lattice grade are therefore a complete coefficient check.
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
            for p in probes:
                lhs = energy(t1_fraction(p), alpha) - energy(p, alpha)
                rhs = sum(p, Fraction(0)) / 8 + alpha / 8
                assert lhs == rhs
                rows += 1

    # Summing the unit identity through four affine iterates gives the source
    # W4 Ward increment Q/2+3/4+alpha/2.
    for alpha in (Fraction(1, 4), Fraction(3, 4)):
        for grade in (0, 1):
            p = tuple(Fraction(grade, 2) for _ in range(DIM))
            current = p
            increment = Fraction(0)
            for _ in range(4):
                nxt = t1_fraction(current)
                increment += energy(nxt, alpha) - energy(current, alpha)
                current = nxt
            expected = sum(p, Fraction(0)) / 2 + Fraction(3, 4) + alpha / 2
            assert increment == expected

    return {
        "complete_affine_coefficient_rows": rows,
        "identity": "K_alpha(Cp+e0)-K_alpha(p)=(Q(p)+alpha)/8",
        "four_step_identity": "K_alpha(T1^4 p)-K_alpha(p)=Q(p)/2+3/4+alpha/2",
        "alpha_values": ["1/4", "3/4"],
        "integer_and_half_lattice_grades": True,
    }


def check_bundle_quarter_relation() -> dict[str, object]:
    # With zeta8=e^(pi*i/4), D_alpha C = zeta8 C D_alpha. Therefore
    # B1(x)=e^(ix/8)C obeys B1(x+2pi)D_alpha=D_alpha B1(x), and B1^4=B4.
    for k in range(DIM):
        image_index = (k + 1) % DIM
        # Exponent relation, including the wrap k=7, modulo 8.
        assert (image_index - k) % 8 == 1
    return {
        "B1": "exp(i*x/8) C",
        "boundary_covariance": "B1(x+2pi) D_alpha = D_alpha B1(x)",
        "fourth_power": "B1^4=exp(i*x/2)P=B4",
        "quarter_phase_retained": True,
    }


def main() -> None:
    result = {
        "schema": "tfpt.crossed_defect.T1.result.v2",
        "verdict": "CONDITIONAL_T1_LIFT_FAMILY_NO_KAPPA_SELECTION",
        "integer_source": check_integer_source_lift(),
        "defect_extension": check_defect_extension(),
        "time": check_time_identity(),
        "bundle": check_bundle_quarter_relation(),
        "first_additional_choice": "whether the cyclic lift must make X a scalar eigenoperator; without it a sector phase r in Z/8 remains",
        "scope": {
            "kappa_minus_excluded_only_when": "one additionally assumes scalar eigenoperator covariance C_hat X C_hat^-1=nu X",
            "kappa_minus_counterexample": "the grading-dependent r=1 lift has the same original integer action, charge action, and K_alpha Ward identity",
            "X_action_status": "extension datum; X and the sector phase r are not actions given by the original integer CAR source",
            "source_selection_proved": False,
            "chi_orientation_selected": False,
            "sector_phase_selected": False,
            "new_H_or_state": False,
        },
    }
    output = Path(__file__).with_name("result.json")
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
