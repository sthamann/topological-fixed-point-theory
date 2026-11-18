#!/usr/bin/env python3
"""Exact ground and state-algebra checks for the fixed K_alpha.

There is no cutoff scan here.  The infinite-lattice conclusion is supplied
by the separable strict-convexity inequalities recorded in PROOF.md; this
file checks their exact coefficients and the finite ground-doublet algebra.
"""

from __future__ import annotations

import json
from fractions import Fraction as F
from pathlib import Path


DIM = 8
ALPHAS = (F(1, 4), F(3, 4))
S = (F(1, 2),) * DIM
R0 = (F(1, 2),) * 4 + (F(-1, 2),) * 4

Matrix = tuple[tuple[complex, ...], ...]


def delta(alpha: F) -> tuple[F, ...]:
    return tuple((F(k) + alpha) / 8 - F(1, 2) for k in range(DIM))


def dot(x: tuple[F, ...], y: tuple[F, ...]) -> F:
    return sum((a * b for a, b in zip(x, y)), F(0))


def energy(p: tuple[F, ...], alpha: F, oscillator_level: int = 0) -> F:
    return F(oscillator_level) + dot(p, p) / 2 + dot(delta(alpha), p)


def matmul(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(2)) for j in range(2))
        for i in range(2)
    )


def add(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(left[i][j] + right[i][j] for j in range(2)) for i in range(2)
    )


def scale(value: complex, matrix: Matrix) -> Matrix:
    return tuple(tuple(value * entry for entry in row) for row in matrix)


def adjoint(matrix: Matrix) -> Matrix:
    return tuple(tuple(matrix[j][i].conjugate() for j in range(2)) for i in range(2))


def trace(matrix: Matrix) -> complex:
    return matrix[0][0] + matrix[1][1]


I2: Matrix = ((1, 0), (0, 1))
ZERO: Matrix = ((0, 0), (0, 0))
DUAL: Matrix = ((1, 0), (0, -1))
RC: Matrix = ((0, 1), (1, 0))
H_R0: Matrix = ((0, 0), (0, 2))
P0: Matrix = ((1, 0), (0, 0))
PR: Matrix = ((0, 0), (0, 1))
MIXED: Matrix = ((F(1, 2), 0), (0, F(1, 2)))
SHIFT_PLUS: Matrix = ((0, 0), (1, 0))
SHIFT_MINUS: Matrix = ((0, 1), (0, 0))


def expectation(rho: Matrix, observable: Matrix) -> complex:
    return trace(matmul(rho, observable))


def check_ground_classification() -> dict[str, object]:
    rows = []
    for alpha in ALPHAS:
        d = delta(alpha)
        assert all(abs(x) < F(1, 2) for x in d)
        assert all(x < 0 for x in d[:4])
        assert all(x > 0 for x in d[4:])

        # For integer m != 0:
        # f_k(m)-f_k(0) >= |m|(|m|/2-|delta_k|) > 0.
        integer_margins = tuple(F(1, 2) - abs(x) for x in d)
        assert all(x > 0 for x in integer_margins)

        # On Z+1/2 the unique closest point to -delta_k is
        # -sign(delta_k)/2.  Its coordinate energy is 1/8-|delta_k|/2.
        half_minimizer = tuple(F(1, 2) if x < 0 else F(-1, 2) for x in d)
        assert half_minimizer == R0
        assert sum(abs(x) for x in d) == 2
        assert energy((F(0),) * DIM, alpha) == 0
        assert energy(R0, alpha) == 0

        # The next half-integer in either direction has a strictly positive
        # increment.  Sign flip costs |delta_k|; magnitude 1/2 -> 3/2 costs
        # 1-|delta_k|.  These are the only nearest competitors per coordinate.
        sign_flip_costs = tuple(abs(x) for x in d)
        magnitude_costs = tuple(1 - abs(x) for x in d)
        assert all(x > 0 for x in sign_flip_costs)
        assert all(magnitude_costs[k] > sign_flip_costs[k] for k in range(DIM))

        opposite_spinor_gap = min(sign_flip_costs)
        assert opposite_spinor_gap == F(1, 32)

        rows.append(
            {
                "alpha": str(alpha),
                "sum_abs_delta": str(sum(abs(x) for x in d)),
                "integer_unique_minimizer": [0] * DIM,
                "half_unique_minimizer": [str(x) for x in half_minimizer],
                "both_minimum_energy": "0",
                "opposite_spinor_gap": str(opposite_spinor_gap),
            }
        )

    # r0=s+(0,0,0,0,-1,-1,-1,-1), whose integer offset has even sum.
    offset = tuple(r - s for r, s in zip(R0, S))
    assert all(x.denominator == 1 for x in offset)
    assert int(sum(offset)) % 2 == 0

    return {
        "quarter_cases": rows,
        "unbounded_lattice_argument": "separable strict convexity",
        "full_Z8_union_half_ground_count_at_N0": 2,
        "Ls_D8_union_s_plus_D8_ground_count_at_N0": 2,
        "ground_labels": ["0", "r0"],
        "r0_in_s_plus_D8": True,
        "positive_oscillator_level_excluded": True,
    }


def check_ground_doublet_symmetries() -> dict[str, object]:
    # The dual sector grading is sigma_z.  Any Rc phase can be removed by a
    # basis rephasing on this doublet, so the source involution is sigma_x.
    assert matmul(DUAL, DUAL) == I2
    assert matmul(RC, RC) == I2
    assert add(matmul(DUAL, RC), matmul(RC, DUAL)) == ZERO

    # Distinct dual-invariant zero-energy densities: no uniqueness.
    assert matmul(matmul(DUAL, P0), DUAL) == P0
    assert matmul(matmul(DUAL, PR), DUAL) == PR
    assert P0 != PR

    # Simultaneous dual and Rc invariance gives the unique density I/2.
    assert matmul(matmul(RC, MIXED), RC) == MIXED
    assert matmul(matmul(DUAL, MIXED), DUAL) == MIXED

    # Rc swaps the two pure ground states.
    assert matmul(matmul(RC, P0), RC) == PR
    assert matmul(matmul(RC, PR), RC) == P0

    # Full Cartan vacuum value <H_r0,0>=0 selects P0.  I/2 instead has value1.
    assert expectation(P0, H_R0) == 0
    assert expectation(PR, H_R0) == 2
    assert expectation(MIXED, H_R0) == 1

    return {
        "dual_and_Rc_anticommute_on_ground_doublet": True,
        "zero_energy_plus_dual_invariance_density_unique": False,
        "dual_invariant_density_family": "diag(a,1-a), 0<=a<=1",
        "zero_energy_plus_dual_plus_Rc_unique_density": "I2/2",
        "Rc_invariance_is_source_selected": False,
        "full_Cartan_integer_vacuum_density": "|0><0|",
        "Cartan_vacuum_and_Rc_invariance_compatible": False,
    }


def check_algebra_resolution() -> dict[str, object]:
    # Compression of translations to span{|0>,|r0>}.  Cocycles can attach
    # unit phases to SHIFT_PLUS/MINUS; diagonal expectations remain zero.
    probes = {
        "identity": I2,
        "E_r0": SHIFT_PLUS,
        "E_minus_r0": SHIFT_MINUS,
        "other_nonzero_shift": ZERO,
    }
    bare_values = {}
    for name, rho in (("vacuum", P0), ("half_ground", PR), ("mixed", MIXED)):
        bare_values[name] = {
            probe_name: int(expectation(rho, observable).real)
            for probe_name, observable in probes.items()
        }
    assert bare_values["vacuum"] == bare_values["half_ground"]
    assert bare_values["vacuum"] == bare_values["mixed"]

    return {
        "bare_crossed_shift_algebra": {
            "dual_invariant_ground_functional_unique": True,
            "functional": "omega(E_n X^c)=1 only for n=0,c=0 and 0 otherwise",
            "density_matrix_unique_in_ground_representation": False,
            "reason": "all diagonal ground mixtures have the same translation expectations",
            "exact_probe_values": bare_values,
        },
        "Cartan_extended_algebra": {
            "ground_states_distinguished": True,
            "H_r0_zero_mode_on_doublet": [[0, 0], [0, 2]],
            "original_integer_vacuum_Ward_value": 0,
            "unique_normal_ground_density_with_that_value": "|0><0|",
        },
        "local_lattice_VOA": {
            "contains_Cartan_and_root_modes": True,
            "VOA_vacuum_vector": "|0>",
            "r0_role": "charged weight-one primary and K-ground, not a second VOA vacuum",
            "Rc_vacuum_preserving": False,
        },
    }


def main() -> None:
    result = {
        "contract": "TFPT.SEAM.GROUND_STATE.ALGEBRA.01",
        "verdict": "EXACT_GROUND_DOUBLET_DUAL_Z2_INSUFFICIENT_CARTAN_VACUUM_UNIQUE",
        "ground_classification": check_ground_classification(),
        "symmetry_selection": check_ground_doublet_symmetries(),
        "algebra_resolution": check_algebra_resolution(),
        "state_selection_summary": {
            "K_nonnegative_alone": "selects no state",
            "normal_state_with_omega_K_zero": "supported on span{|0>,|r0>}",
            "plus_dual_Z2": "one-parameter density simplex",
            "plus_dual_Z2_and_Rc": "I2/2, only if Rc invariance is additionally imposed",
            "plus_original_full_Cartan_vacuum_values": "|0><0|",
            "Rc_and_original_Cartan_vacuum": "incompatible",
        },
        "unchanged_source_compatibility": {
            "reused_source_theorem": "UR.SEAM.QUOTIENT_LOCALITY_E8.01",
            "original_even_CAR_ground_multiplicity": 1,
            "crossed_Ls_extension_ground_multiplicity": 2,
            "depends_on_chi_or_kappa": False,
            "all_current_crossed_defect_candidates_compatible_with_unchanged_source": False,
            "classification_only_status": "conditional algebraic extensions",
            "reopening_requires": (
                "a source-derived extra sector/representation or a source-derived changed "
                "physical time/domain; choosing a density does not remove the extra zero vector"
            ),
        },
        "new_H_state_or_symmetry_assumption_added": False,
    }
    output = Path(__file__).with_name("result.json")
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
