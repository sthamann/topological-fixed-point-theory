#!/usr/bin/env python3
"""Exact finite checks for three inequivalent repeatable-observation models.

This is deliberately a one-system, finite-horizon operational comparison.  It
does not assume that observation determines probabilities, a state space, a
tensor product, locality, or a preparation-noncontextual hidden model.
"""

from __future__ import annotations

import json
from fractions import Fraction
from itertools import product
from pathlib import Path
from typing import Dict, Iterable, Mapping, Optional, Sequence, Tuple


Q = Fraction
Sign = int
Point = Tuple[Sign, Sign]
Vector = Tuple[Q, Q, Q]
Matrix = Tuple[Tuple[Q, Q, Q], Tuple[Q, Q, Q], Tuple[Q, Q, Q]]
OutcomeTable = Dict[Tuple[Sign, ...], Q]

SIGNS: Tuple[Sign, Sign] = (1, -1)
POINTS: Tuple[Point, ...] = tuple(product(SIGNS, repeat=2))
ZERO_VECTOR: Vector = (Q(0), Q(0), Q(0))


def qstr(value: Q) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def sign_string(outcomes: Sequence[Sign]) -> str:
    return "".join("+" if value == 1 else "-" for value in outcomes)


def point_name(point: Point) -> str:
    return sign_string(point)


def vector_strings(vector: Vector) -> list[str]:
    return [qstr(value) for value in vector]


def table_strings(table: Mapping[Tuple[Sign, ...], Q]) -> Dict[str, str]:
    return {
        sign_string(outcome): qstr(probability)
        for outcome, probability in sorted(table.items(), key=lambda item: sign_string(item[0]))
        if probability
    }


def add_vectors(left: Vector, right: Vector) -> Vector:
    return tuple(a + b for a, b in zip(left, right))  # type: ignore[return-value]


def matrix_add(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(left[row][column] + right[row][column] for column in range(3))
        for row in range(3)
    )  # type: ignore[return-value]


def matrix_multiply(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(sum(left[row][k] * right[k][column] for k in range(3)) for column in range(3))
        for row in range(3)
    )  # type: ignore[return-value]


def matrix_apply(matrix: Matrix, vector: Vector) -> Vector:
    return tuple(
        sum(matrix[row][column] * vector[column] for column in range(3))
        for row in range(3)
    )  # type: ignore[return-value]


ZERO_MATRIX: Matrix = (
    (Q(0), Q(0), Q(0)),
    (Q(0), Q(0), Q(0)),
    (Q(0), Q(0), Q(0)),
)


def d_branch_matrix(measurement: str, outcome: Sign) -> Matrix:
    half = Q(1, 2)
    if measurement == "X":
        return (
            (half, outcome * half, Q(0)),
            (outcome * half, half, Q(0)),
            (Q(0), Q(0), Q(0)),
        )
    if measurement == "Z":
        return (
            (half, Q(0), outcome * half),
            (Q(0), Q(0), Q(0)),
            (outcome * half, Q(0), half),
        )
    raise ValueError(f"unknown D measurement: {measurement}")


def d_branch(measurement: str, outcome: Sign, vector: Vector) -> Vector:
    return matrix_apply(d_branch_matrix(measurement, outcome), vector)


def d_nonselective_matrix(measurement: str) -> Matrix:
    return matrix_add(d_branch_matrix(measurement, 1), d_branch_matrix(measurement, -1))


def d_positive(vector: Vector) -> bool:
    t, x, z = vector
    return t >= 0 and t >= abs(x) and t >= abs(z)


def d_history_distribution(vector: Vector, word: str) -> OutcomeTable:
    branches: Dict[Tuple[Sign, ...], Vector] = {(): vector}
    for measurement in word:
        next_branches: Dict[Tuple[Sign, ...], Vector] = {}
        for history, state in branches.items():
            for outcome in SIGNS:
                next_branches[history + (outcome,)] = d_branch(measurement, outcome, state)
        branches = next_branches
    return {history: state[0] for history, state in branches.items()}


def i_branch(
    measurement: str, outcome: Sign, distribution: Mapping[Point, Q]
) -> Dict[Point, Q]:
    def value(point: Point) -> Sign:
        a, b = point
        if measurement == "X":
            return a
        if measurement == "Z":
            return b
        if measurement == "Y":
            return a * b
        raise ValueError(f"unknown I measurement: {measurement}")

    return {
        point: probability if value(point) == outcome else Q(0)
        for point, probability in distribution.items()
    }


def i_reset(distribution: Mapping[Point, Q], target: Point) -> Dict[Point, Q]:
    mass = sum(distribution.values())
    return {point: mass if point == target else Q(0) for point in POINTS}


def i_history_distribution(distribution: Mapping[Point, Q], word: str) -> OutcomeTable:
    branches: Dict[Tuple[Sign, ...], Dict[Point, Q]] = {(): dict(distribution)}
    for measurement in word:
        next_branches: Dict[Tuple[Sign, ...], Dict[Point, Q]] = {}
        for history, state in branches.items():
            for outcome in SIGNS:
                next_branches[history + (outcome,)] = i_branch(measurement, outcome, state)
        branches = next_branches
    return {history: sum(state.values()) for history, state in branches.items()}


def matrix_rank(rows: Sequence[Sequence[Q]]) -> int:
    matrix = [list(row) for row in rows]
    if not matrix:
        return 0
    row_count = len(matrix)
    column_count = len(matrix[0])
    rank = 0
    for column in range(column_count):
        pivot = next((row for row in range(rank, row_count) if matrix[row][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        pivot_value = matrix[rank][column]
        matrix[rank] = [value / pivot_value for value in matrix[rank]]
        for row in range(row_count):
            if row != rank and matrix[row][column]:
                multiplier = matrix[row][column]
                matrix[row] = [
                    value - multiplier * pivot_entry
                    for value, pivot_entry in zip(matrix[row], matrix[rank])
                ]
        rank += 1
        if rank == row_count:
            break
    return rank


def all_words(max_length: int) -> list[str]:
    return [
        "".join(letters)
        for length in range(1, max_length + 1)
        for letters in product("XZ", repeat=length)
    ]


def r_selective_read(state: int, reported: int) -> Optional[int]:
    return state if state == reported else None


def build_history_tree(
    preparation: str, initial: Vector, horizon: int
) -> Tuple[Dict[Tuple[str, Tuple[Sign, ...]], Vector], Dict[Tuple[str, Tuple[Sign, ...], str, Sign], Q]]:
    """Construct a classical prefix-memory simulation for one D preparation.

    A memory state is the preparation label plus the complete settings/outcomes
    prefix.  Zero-mass prefixes are omitted.  Thus this is an unrestricted
    finite history reconstruction, not a noncontextual ontological model.
    """

    states: Dict[Tuple[str, Tuple[Sign, ...]], Vector] = {("", ()): initial}
    transitions: Dict[Tuple[str, Tuple[Sign, ...], str, Sign], Q] = {}
    frontier = [("", (), initial)]
    for _depth in range(horizon):
        next_frontier = []
        for settings, outcomes, state in frontier:
            mass = state[0]
            assert mass > 0
            for measurement in "XZ":
                for outcome in SIGNS:
                    next_state = d_branch(measurement, outcome, state)
                    conditional = next_state[0] / mass
                    transitions[(settings, outcomes, measurement, outcome)] = conditional
                    if next_state[0] > 0:
                        next_key = (settings + measurement, outcomes + (outcome,))
                        existing = states.get(next_key)
                        if existing is not None:
                            assert existing == next_state
                        else:
                            states[next_key] = next_state
                            next_frontier.append((settings + measurement, outcomes + (outcome,), next_state))
        frontier = next_frontier
    return states, transitions


def history_tree_distribution(
    word: str,
    transitions: Mapping[Tuple[str, Tuple[Sign, ...], str, Sign], Q],
) -> OutcomeTable:
    active: Dict[Tuple[Sign, ...], Q] = {(): Q(1)}
    for index, measurement in enumerate(word):
        settings = word[:index]
        next_active: Dict[Tuple[Sign, ...], Q] = {}
        for outcomes, path_probability in active.items():
            if not path_probability:
                continue
            for outcome in SIGNS:
                conditional = transitions[(settings, outcomes, measurement, outcome)]
                next_probability = path_probability * conditional
                if next_probability:
                    next_active[outcomes + (outcome,)] = next_probability
        active = next_active
    return active


def main() -> None:
    words = all_words(4)
    assert len(words) == 30

    # R: only identity and nondisturbing selective reads on two fixed labels.
    r_reports: Dict[str, Dict[str, str]] = {}
    for initial in (0, 1):
        assert r_selective_read(initial, initial) == initial
        assert r_selective_read(initial, 1 - initial) is None
        r_reports[str(initial)] = {
            str(length): str(initial) * length for length in range(1, 5)
        }

    # Shared named preparations.
    d_corners: Dict[str, Vector] = {
        point_name(point): (Q(1), Q(point[0]), Q(point[1])) for point in POINTS
    }
    d_preparations: Dict[str, Vector] = {
        **d_corners,
        "center": (Q(1), Q(0), Q(0)),
        "diag_even": (Q(1), Q(0), Q(0)),
        "diag_odd": (Q(1), Q(0), Q(0)),
    }
    i_corners: Dict[str, Dict[Point, Q]] = {
        point_name(point): {candidate: Q(candidate == point) for candidate in POINTS}
        for point in POINTS
    }
    i_preparations: Dict[str, Dict[Point, Q]] = {
        **i_corners,
        "center": {point: Q(1, 4) for point in POINTS},
        "diag_even": {point: Q(1, 2) if point in ((1, 1), (-1, -1)) else Q(0) for point in POINTS},
        "diag_odd": {point: Q(1, 2) if point in ((1, -1), (-1, 1)) else Q(0) for point in POINTS},
    }

    # D cone positivity, instrument normalization, and repeatability.
    positivity_checks: Dict[str, Dict[str, bool]] = {}
    repeatability_checks: Dict[str, Dict[str, bool]] = {}
    for measurement in "XZ":
        positivity_checks[measurement] = {}
        repeatability_checks[measurement] = {}
        nonselective = d_nonselective_matrix(measurement)
        assert nonselective[0] == (Q(1), Q(0), Q(0))
        for corner_name, ray in d_corners.items():
            for outcome in SIGNS:
                branch = d_branch(measurement, outcome, ray)
                assert d_positive(branch)
                positivity_checks[measurement][f"{corner_name}:{sign_string((outcome,))}"] = True
        for outcome in SIGNS:
            branch_matrix = d_branch_matrix(measurement, outcome)
            same_twice = matrix_multiply(branch_matrix, branch_matrix)
            opposite_after = matrix_multiply(
                d_branch_matrix(measurement, -outcome), branch_matrix
            )
            assert same_twice == branch_matrix
            assert opposite_after == ZERO_MATRIX
            repeatability_checks[measurement][sign_string((outcome,))] = True

    # All eight square symmetries are signed coordinate permutations.
    square_symmetry_images: Dict[str, list[str]] = {}
    corner_set = set(d_corners.values())
    for swap in (False, True):
        for sign_x in SIGNS:
            for sign_z in SIGNS:
                name = f"{'swap' if swap else 'keep'}:{sign_string((sign_x, sign_z))}"
                images = []
                for corner_name, (t, x, z) in d_corners.items():
                    source_x, source_z = (z, x) if swap else (x, z)
                    image = (t, sign_x * source_x, sign_z * source_z)
                    assert image in corner_set and d_positive(image)
                    images.append(f"{corner_name}->{point_name((int(image[1]), int(image[2])))}")
                square_symmetry_images[name] = images
    assert len(square_symmetry_images) == 8

    # Exact outcome tables for all requested preparations and words.
    d_tables = {
        preparation: {
            word: table_strings(d_history_distribution(vector, word)) for word in words
        }
        for preparation, vector in d_preparations.items()
    }
    i_tables = {
        preparation: {
            word: table_strings(i_history_distribution(distribution, word)) for word in words
        }
        for preparation, distribution in i_preparations.items()
    }
    for preparation in d_preparations:
        for word in words:
            assert sum(d_history_distribution(d_preparations[preparation], word).values()) == 1
            assert sum(i_history_distribution(i_preparations[preparation], word).values()) == 1

    # First and final X agree in I but only with probability 1/2 in D after XZX.
    xzx_match: Dict[str, Dict[str, str]] = {}
    for corner_name in d_corners:
        d_distribution = d_history_distribution(d_corners[corner_name], "XZX")
        i_distribution = i_history_distribution(i_corners[corner_name], "XZX")
        d_match = sum(probability for outcomes, probability in d_distribution.items() if outcomes[0] == outcomes[2])
        i_match = sum(probability for outcomes, probability in i_distribution.items() if outcomes[0] == outcomes[2])
        assert d_match == Q(1, 2)
        assert i_match == Q(1)
        xzx_match[corner_name] = {"I": qstr(i_match), "D": qstr(d_match)}

    # The two diagonal I preparations have opposite parity.  Their D images are
    # the same center and therefore have identical tables for every tested word.
    i_parity = {}
    for preparation, expected in (("diag_even", 1), ("diag_odd", -1)):
        distribution = i_history_distribution(i_preparations[preparation], "Y")
        expectation = sum(outcomes[0] * probability for outcomes, probability in distribution.items())
        assert expectation == expected
        i_parity[preparation] = {
            "distribution": table_strings(distribution),
            "expectation": qstr(expectation),
        }
    for word in words:
        assert d_history_distribution(d_preparations["diag_even"], word) == d_history_distribution(
            d_preparations["diag_odd"], word
        )

    # D affine dependence versus I simplex independence.
    d_order = ["++", "+-", "-+", "--"]
    d_matrix_rows = [[d_corners[name][row] for name in d_order] for row in range(3)]
    i_matrix_rows = [[Q(row == column) for column in range(4)] for row in range(4)]
    d_rank = matrix_rank(d_matrix_rows)
    i_rank = matrix_rank(i_matrix_rows)
    assert d_rank == 3 and i_rank == 4
    d_relation_left = add_vectors(d_corners["++"], d_corners["--"])
    d_relation_right = add_vectors(d_corners["+-"], d_corners["-+"])
    assert d_relation_left == d_relation_right == (Q(2), Q(0), Q(0))

    # Selective D branches do not commute, while their nonselective sums do.
    x_plus = d_branch_matrix("X", 1)
    z_plus = d_branch_matrix("Z", 1)
    x_after_z = matrix_multiply(x_plus, z_plus)
    z_after_x = matrix_multiply(z_plus, x_plus)
    assert x_after_z != z_after_x
    center = d_preparations["center"]
    selective_order_outputs = {
        "X_plus_after_Z_plus": vector_strings(matrix_apply(x_after_z, center)),
        "Z_plus_after_X_plus": vector_strings(matrix_apply(z_after_x, center)),
    }
    nonselective_x = d_nonselective_matrix("X")
    nonselective_z = d_nonselective_matrix("Z")
    nonselective_xz = matrix_multiply(nonselective_x, nonselective_z)
    nonselective_zx = matrix_multiply(nonselective_z, nonselective_x)
    assert nonselective_xz == nonselective_zx
    assert matrix_apply(nonselective_xz, center) == center

    # Optional I reset is a normalized constant stochastic map.
    for preparation in i_preparations.values():
        reset = i_reset(preparation, (1, 1))
        assert sum(reset.values()) == 1 and reset[(1, 1)] == 1

    # General finite classical prefix-memory reconstruction of every D table.
    reconstruction_counts: Dict[str, Dict[str, int]] = {}
    for preparation, vector in d_preparations.items():
        states, transitions = build_history_tree(preparation, vector, 4)
        for word in words:
            reconstructed = history_tree_distribution(word, transitions)
            direct = {
                outcomes: probability
                for outcomes, probability in d_history_distribution(vector, word).items()
                if probability
            }
            assert reconstructed == direct
        reconstruction_counts[preparation] = {
            "reachable_prefix_states": len(states),
            "conditional_response_entries": len(transitions),
        }

    # A compact I discriminator in addition to parity: diagonal correlations.
    i_diagonal_xz = {
        preparation: table_strings(i_history_distribution(i_preparations[preparation], "XZ"))
        for preparation in ("diag_even", "diag_odd")
    }
    assert i_diagonal_xz["diag_even"] != i_diagonal_xz["diag_odd"]

    results = {
        "checker": {
            "arithmetic": "exact fractions",
            "random_sampling": False,
            "maximum_XZ_word_length": 4,
            "XZ_word_count": len(words),
        },
        "definitions_and_extra_assumptions": {
            "R": {
                "primitive_states": [0, 1],
                "allowed_internal_generators": ["identity", "nondisturbing selective read branches"],
                "excluded_internal_generators": ["mixtures", "discard", "reprepare", "flip"],
                "probability_postulate": False,
            },
            "I": {
                "primitive_ontic_points": [point_name(point) for point in POINTS],
                "preparations": "full probability simplex on four points",
                "measurements": ["nondisturbing X=a", "nondisturbing Z=b", "parity Y=ab"],
                "optional_operation": "constant stochastic reset",
                "extra_assumption": "classical convex probability simplex",
            },
            "D": {
                "state_cone": "v=(t,x,z), t>=|x| and t>=|z|; normalized states have t=1",
                "preparations": "all four square corners and their convex mixtures",
                "operations": [
                    "the specified X and Z instrument branches",
                    "all eight square symmetries",
                    "finite serial and adaptive composition",
                    "stochastic coarse-graining",
                ],
                "extra_assumption": "ordered real vector space and convex probabilities",
                "bipartite_structure": "not defined",
            },
        },
        "R_exact_reports": r_reports,
        "D_instrument_checks": {
            "positive_on_four_extreme_rays": positivity_checks,
            "mass_preserving_nonselective_sums": {"X": True, "Z": True},
            "repeatability_and_exclusion": repeatability_checks,
            "square_symmetries": square_symmetry_images,
            "selective_X_plus_Z_plus_noncommuting": {
                "verified": True,
                "outputs_on_center": selective_order_outputs,
            },
            "nonselective_X_Z_commuting": {
                "verified": True,
                "output_on_any_normalized_input": ["1", "0", "0"],
            },
        },
        "exact_XZ_outcome_tables": {"I": i_tables, "D": d_tables},
        "empirical_discriminators": {
            "XZX_first_last_X_match_on_each_corner": xzx_match,
            "diagonal_I_parity": i_parity,
            "diagonal_I_XZ_tables": i_diagonal_xz,
            "diagonal_D_all_XZ_words_equal_through_length_4": True,
            "D_selective_order_effect": selective_order_outputs,
            "D_nonselective_order_effect": "none; both orders send every normalized state to (1,0,0)",
        },
        "affine_structure": {
            "D_corner_relation": "v(++ ) + v(--) = v(+-) + v(-+) = (2,0,0)",
            "D_homogeneous_corner_column_rank": d_rank,
            "I_simplex_corner_column_rank": i_rank,
            "interpretation_limit": "These ranks count coordinates in the declared operational representations; they are not physical space or spacetime dimensions.",
        },
        "finite_history_tree_reconstruction": {
            "horizon": 4,
            "all_requested_D_tables_matched_exactly": True,
            "counts": reconstruction_counts,
            "memory_key": "preparation label plus complete settings/outcomes prefix",
            "scope_limit": "This unrestricted, preparation- and history-dependent simulator establishes finite table reproducibility only. It assumes no locality, bounded memory, preparation noncontextuality, or common hidden state across operationally equivalent preparations, and is not a universal classical physical model.",
        },
        "verdict": {
            "observed": "R, I, and D all contain repeatable distinctions under their declared instruments, but they have inequivalent sequential and preparation statistics.",
            "necessary_correction": "Repeatable observation alone selects neither probability, the classical simplex, nor the square state cone.",
            "countermodel_scope": "D's finite observed tables admit an exact classical history-tree reconstruction when unrestricted memory and preparation dependence are allowed. This does not identify D with a simplex model under stronger operational constraints.",
            "not_claimed": [
                "that non-simplex behavior excludes every classical hidden-memory simulation",
                "that probability follows from observation",
                "a bipartite D theory, locality, quantum theory, or physical dimension",
            ],
        },
    }

    output_path = Path(__file__).with_name("observability_results.json")
    output_path.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        "PASS: exact R/I/D checks; 30 XZ words through length 4; "
        "D history-tree reconstruction matched every requested table."
    )


if __name__ == "__main__":
    main()
