#!/usr/bin/env python3
"""Exact relational channel on conjugate(H5) tensor H5.

The only additional contract is the presence of the relational tensor sector
with the diagonal G31 action.  All representation data are rebuilt through
``doily_quartic_bridge`` from the tracked compiler sources.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import sympy as sp


CONTRACT_DIR = Path(__file__).resolve().parent
if str(CONTRACT_DIR) not in sys.path:
    sys.path.insert(0, str(CONTRACT_DIR))

import doily_quartic_bridge as bridge  # noqa: E402


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def parse_scalar(value: int | str) -> sp.Expr:
    return sp.Integer(value) if isinstance(value, int) else sp.Rational(value)


def parse_matrix(values: list[list[int | str]]) -> sp.Matrix:
    return sp.Matrix([[parse_scalar(value) for value in row] for row in values])


def partial_trace_first(matrix: sp.Matrix, dimension: int) -> sp.Matrix:
    """Trace the first factor in an orthonormal product basis."""

    return sp.Matrix(
        dimension,
        dimension,
        lambda row, column: sum(
            matrix[index * dimension + row, index * dimension + column]
            for index in range(dimension)
        ),
    )


def partial_trace_second(matrix: sp.Matrix, dimension: int) -> sp.Matrix:
    """Trace the second factor in an orthonormal product basis."""

    return sp.Matrix(
        dimension,
        dimension,
        lambda row, column: sum(
            matrix[row * dimension + index, column * dimension + index]
            for index in range(dimension)
        ),
    )


def column_vectorize(matrix: sp.Matrix) -> sp.Matrix:
    """Column vectorization with conjugate-input index first."""

    rows, columns = matrix.shape
    return sp.Matrix(
        [matrix[output_index, input_index] for input_index in range(columns) for output_index in range(rows)]
    )


def matrices_equal(left: sp.Matrix, right: sp.Matrix) -> bool:
    """Exact equality after simplifying algebraic entries."""

    if left.shape != right.shape:
        return False
    difference = (left - right).applyfunc(sp.simplify)
    return difference == sp.zeros(*left.shape)


def matrix_json(matrix: sp.Matrix) -> list[list[int | str]]:
    return bridge.matrix_json(matrix)


def scalar_json(value: sp.Expr) -> int | str:
    return bridge.rational_json(value)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_report() -> dict[str, Any]:
    bridge_report = bridge.build_report()
    require(bridge_report["status"] == "PASS", "doily quartic bridge did not pass")

    dimension = 5
    identity5 = sp.eye(dimension)
    intertwiner = parse_matrix(
        bridge_report["quartic_bridge"]["primitive_intertwiner"]
    )
    metric = parse_matrix(
        bridge_report["selection_data"]["primitive_native_invariant_positive_form"]
    )
    isometric_metric = metric / 12
    centered_metric = sp.eye(5) + sp.ones(5)
    require(
        intertwiner.T * isometric_metric * intertwiner == centered_metric,
        "the existing bridge is not isometric for the declared metrics",
    )

    native_rays = [
        sp.Matrix([parse_scalar(value) for value in ray])
        for ray in bridge_report["selection_data"]["native_ovoid_rays"]
    ]
    require(len(native_rays) == 6, "native ovoid ray count changed")
    require(
        {(ray.T * metric * ray)[0] for ray in native_rays} == {10},
        "native ovoid ray norm changed",
    )

    effects = [ray * ray.T * metric / 12 for ray in native_rays]
    require(sum(effects, sp.zeros(5)) == identity5, "six native effects do not sum to identity")
    require(
        all(metric.inv() * effect.T * metric == effect for effect in effects),
        "an effect is not self-adjoint for the native metric",
    )
    require(
        all(effect.rank() == 1 and effect.trace() == sp.Rational(5, 6) for effect in effects),
        "native effect rank or trace changed",
    )

    quartic_actions = [
        parse_matrix(row["quartic_numerator"]) / 16
        for row in bridge_report["quartic_bridge"]["per_reflection"]
    ]
    ovoid_permutations = [
        row["ovoid_permutation"]
        for row in bridge_report["quartic_bridge"]["per_reflection"]
    ]
    require(len(quartic_actions) == 60, "quartic reflection count changed")
    require(
        all(action.T * metric * action == metric for action in quartic_actions),
        "a native quartic action does not preserve the positive form",
    )

    # In the word basis the canonical tensor representing the identity map is
    # vec(G^-1).  Its norm for G tensor G is five.  This is the same ray as
    # sum_i conjugate(e_i) tensor e_i in a G-orthonormal basis.
    omega_word = sp.zeros(25, 1)
    inverse_metric = metric.inv()
    for index in range(5):
        omega_word[index * 5 + index] = inverse_metric[index, index]
    tensor_metric = sp.kronecker_product(metric, metric)
    omega_word_norm = (omega_word.T * tensor_metric * omega_word)[0]
    require(omega_word_norm == 5, "singlet norm changed")
    singlet_projector_word = omega_word * omega_word.T * tensor_metric / 5
    require(singlet_projector_word.rank() == 1, "singlet projector rank changed")
    require(singlet_projector_word**2 == singlet_projector_word, "singlet projector is not idempotent")
    require(singlet_projector_word.trace() == 1, "singlet projector trace changed")
    require(
        singlet_projector_word.T * tensor_metric
        == tensor_metric * singlet_projector_word,
        "singlet projector is not self-adjoint for the tensor metric",
    )

    invariant_constraints: list[sp.Matrix] = []
    singlet_covariance: list[bool] = []
    effect_covariance: list[bool] = []
    relational_covariance: list[bool] = []

    relational_operator = sum(
        (
            sp.kronecker_product(effect.conjugate(), effect)
            for effect in effects
        ),
        sp.zeros(25),
    )
    require(relational_operator.rank() > 0, "relational operator vanished")
    require(
        relational_operator != sp.eye(25) * relational_operator.trace() / 25,
        "relational operator became scalar",
    )

    for action, permutation in zip(quartic_actions, ovoid_permutations):
        # The first factor is the conjugate representation.  The matrices happen
        # to be real in the exact word basis, but conjugation is retained here.
        tensor_action = sp.kronecker_product(action.conjugate(), action)
        invariant_constraints.append(tensor_action - sp.eye(25))
        singlet_ok = tensor_action * omega_word == omega_word
        singlet_covariance.append(singlet_ok)
        effect_ok = all(
            action * effects[index] * action.inv()
            == effects[permutation[index]]
            for index in range(6)
        )
        effect_covariance.append(effect_ok)
        relational_ok = (
            tensor_action
            * relational_operator
            * tensor_action.inv()
            == relational_operator
        )
        relational_covariance.append(relational_ok)
    require(all(singlet_covariance), "singlet ray is not invariant under all sixty actions")
    require(all(effect_covariance), "effect family is not covariant under all sixty actions")
    require(all(relational_covariance), "relational operator is not invariant under all sixty actions")

    invariant_matrix = sp.Matrix.vstack(*invariant_constraints)
    singlet_dimension = 25 - invariant_matrix.rank()
    require(singlet_dimension == 1, "relational singlet is not unique")

    # Move to a G-orthonormal basis for ordinary density matrices and partial
    # traces.  S maps word coordinates to orthonormal coordinates.
    square_root_metric = sp.diag(1, sp.sqrt(3), sp.sqrt(3), sp.sqrt(3), sp.sqrt(6))
    require(square_root_metric.T * square_root_metric == metric, "metric square root changed")
    orthogonal_actions = [
        sp.simplify(square_root_metric * action * square_root_metric.inv())
        for action in quartic_actions
    ]
    require(
        all(action.T * action == identity5 for action in orthogonal_actions),
        "orthonormalized quartic action is not orthogonal",
    )
    orthogonal_effects = [
        sp.simplify(square_root_metric * effect * square_root_metric.inv())
        for effect in effects
    ]

    # In the orthonormal chart R is Hermitian.  The word-basis matrix is
    # rational and similar to it, so its exact ranks and annihilating
    # polynomial provide a fast exact certificate for the full spectrum.
    tensor_isometry = sp.kronecker_product(square_root_metric, square_root_metric)
    relational_operator_orthogonal = sp.simplify(
        tensor_isometry * relational_operator * tensor_isometry.inv()
    )
    require(
        relational_operator_orthogonal.H == relational_operator_orthogonal,
        "relational operator is not Hermitian in the isometric chart",
    )
    require(
        matrices_equal(
            relational_operator_orthogonal,
            sum(
                (
                    sp.kronecker_product(effect.conjugate(), effect)
                    for effect in orthogonal_effects
                ),
                sp.zeros(25),
            ),
        ),
        "relational operator did not transform to the isometric chart",
    )
    identity25 = sp.eye(25)
    eigenvalue_top = sp.Rational(5, 6)
    eigenvalue_secondary = sp.Rational(2, 3)
    relational_nullities = {
        "0": 25 - relational_operator.rank(),
        "5/6": 25 - (relational_operator - eigenvalue_top * identity25).rank(),
        "2/3": 25 - (relational_operator - eigenvalue_secondary * identity25).rank(),
    }
    require(
        relational_nullities == {"0": 19, "5/6": 1, "2/3": 5},
        "relational eigenspace multiplicities changed",
    )
    relational_annihilator = (
        relational_operator
        * (relational_operator - eigenvalue_top * identity25)
        * (relational_operator - eigenvalue_secondary * identity25)
    )
    require(
        relational_annihilator == sp.zeros(25),
        "relational spectral polynomial failed",
    )

    omega_orthogonal = sp.zeros(25, 1)
    for index in range(5):
        omega_orthogonal[index * 5 + index] = 1
    singlet_density = omega_orthogonal * omega_orthogonal.T / 5
    require(singlet_density.rank() == 1 and singlet_density.trace() == 1, "orthonormal singlet density failed")
    require(singlet_density**2 == singlet_density, "orthonormal singlet density is not pure")
    marginal_first = partial_trace_second(singlet_density, 5)
    marginal_second = partial_trace_first(singlet_density, 5)
    require(marginal_first == identity5 / 5, "first singlet marginal is not maximally mixed")
    require(marginal_second == identity5 / 5, "second singlet marginal is not maximally mixed")
    require(
        relational_operator_orthogonal * omega_orthogonal
        == eigenvalue_top * omega_orthogonal,
        "Omega is not the top relational eigenvector",
    )
    singlet_relational_expectation = sp.simplify(
        sp.trace(singlet_density * relational_operator_orthogonal)
    )
    require(
        singlet_relational_expectation == eigenvalue_top,
        "singlet relational expectation changed",
    )

    # Symmetry alone does not select the pure singlet density: the maximally
    # mixed density is invariant as well.  It has the same local marginals but
    # a different joint ovoid law, checked below.
    maximally_mixed_tensor_density = identity25 / 25
    require(
        partial_trace_second(maximally_mixed_tensor_density, 5) == identity5 / 5
        and partial_trace_first(maximally_mixed_tensor_density, 5) == identity5 / 5,
        "maximally mixed tensor marginals changed",
    )
    require(
        all(
            sp.kronecker_product(action.conjugate(), action)
            * maximally_mixed_tensor_density
            * sp.kronecker_product(action.conjugate(), action).H
            == maximally_mixed_tensor_density
            for action in orthogonal_actions
        ),
        "maximally mixed tensor density is not invariant",
    )
    mixed_relational_expectation = sp.simplify(
        sp.trace(maximally_mixed_tensor_density * relational_operator_orthogonal)
    )
    require(
        mixed_relational_expectation == sp.Rational(1, 6),
        "maximally mixed relational expectation changed",
    )

    # The unnormalised Choi operator of the channel is 5 times the normalized
    # singlet state.  Unit partial trace proves trace preservation; rank-one
    # positivity proves complete positivity.  It is the identity channel under
    # the canonical conjugate-input identification.
    choi_identity_channel = 5 * singlet_density
    require(
        partial_trace_first(choi_identity_channel, 5) == identity5
        and partial_trace_second(choi_identity_channel, 5) == identity5,
        "identity-channel Choi marginals changed",
    )

    joint_probabilities = sp.zeros(6)
    direct_probabilities = sp.zeros(6)
    for left in range(6):
        left_adjoint = metric.inv() * effects[left].T.conjugate() * metric
        for right in range(6):
            trace_formula = sp.simplify(
                sp.trace(left_adjoint * effects[right]) / 5
            )
            direct = sp.simplify(
                (
                    omega_orthogonal.T
                    * sp.kronecker_product(
                        orthogonal_effects[left].conjugate(),
                        orthogonal_effects[right],
                    )
                    * omega_orthogonal
                )[0]
                / 5
            )
            require(direct == trace_formula, "conjugate tensor expectation formula failed")
            joint_probabilities[left, right] = trace_formula
            direct_probabilities[left, right] = direct

    expected_same = sp.Rational(5, 36)
    expected_distinct = sp.Rational(1, 180)
    require(
        all(joint_probabilities[index, index] == expected_same for index in range(6)),
        "same-outcome probability changed",
    )
    require(
        all(
            joint_probabilities[left, right] == expected_distinct
            for left in range(6)
            for right in range(6)
            if left != right
        ),
        "distinct-outcome probability changed",
    )
    row_marginals = [sum(joint_probabilities[row, column] for column in range(6)) for row in range(6)]
    column_marginals = [sum(joint_probabilities[row, column] for row in range(6)) for column in range(6)]
    require(set(row_marginals) == {sp.Rational(1, 6)}, "row marginal changed")
    require(set(column_marginals) == {sp.Rational(1, 6)}, "column marginal changed")
    conditional_same = expected_same / row_marginals[0]
    require(conditional_same == sp.Rational(5, 6), "conditional same probability changed")
    require(sum(joint_probabilities) == 1, "joint probabilities are not normalized")

    mixed_joint_probabilities = sp.zeros(6)
    for left in range(6):
        for right in range(6):
            mixed_joint_probabilities[left, right] = sp.simplify(
                sp.trace(
                    maximally_mixed_tensor_density
                    * sp.kronecker_product(
                        orthogonal_effects[left].conjugate(),
                        orthogonal_effects[right],
                    )
                )
            )
    require(
        set(mixed_joint_probabilities) == {sp.Rational(1, 36)},
        "maximally mixed tensor joint law changed",
    )
    require(
        mixed_joint_probabilities != joint_probabilities,
        "invariant mixed-state counterexample did not distinguish the singlet law",
    )

    # The normalized square-root (Lueders) instrument of the already given
    # effects has Kraus operators K_a=sqrt(6/5) E_a because E_a has the sole
    # nonzero eigenvalue 5/6.  This is an operational instrument prescription,
    # not a time evolution derived from symmetry.
    kraus_scale = sp.sqrt(sp.Rational(6, 5))
    lueders_kraus = [sp.simplify(kraus_scale * effect) for effect in orthogonal_effects]
    require(
        all(effect.H == effect for effect in orthogonal_effects),
        "an isometric-chart effect is not Hermitian",
    )
    require(
        all(matrices_equal(kraus.H * kraus, effect) for kraus, effect in zip(lueders_kraus, orthogonal_effects)),
        "K_a is not the positive square root of E_a",
    )
    kraus_input_sum = sum((kraus.H * kraus for kraus in lueders_kraus), sp.zeros(5))
    kraus_output_sum = sum((kraus * kraus.H for kraus in lueders_kraus), sp.zeros(5))
    require(matrices_equal(kraus_input_sum, identity5), "Lueders channel is not trace preserving")
    require(matrices_equal(kraus_output_sum, identity5), "Lueders channel is not unital")

    lueders_liouville = sum(
        (sp.kronecker_product(kraus.conjugate(), kraus) for kraus in lueders_kraus),
        sp.zeros(25),
    )
    require(
        matrices_equal(
            lueders_liouville,
            sp.Rational(6, 5) * relational_operator_orthogonal,
        ),
        "Lueders Liouville matrix is not (6/5)R",
    )
    lueders_word = sp.Rational(6, 5) * relational_operator
    lueders_nullities = {
        "0": 25 - lueders_word.rank(),
        "1": 25 - (lueders_word - identity25).rank(),
        "4/5": 25 - (lueders_word - sp.Rational(4, 5) * identity25).rank(),
    }
    require(
        lueders_nullities == {"0": 19, "1": 1, "4/5": 5},
        "Lueders eigenspace multiplicities changed",
    )
    require(
        lueders_word
        * (lueders_word - identity25)
        * (lueders_word - sp.Rational(4, 5) * identity25)
        == sp.zeros(25),
        "Lueders spectral polynomial failed",
    )
    require(
        lueders_liouville * omega_orthogonal == omega_orthogonal,
        "identity vector is not stationary for the Lueders channel",
    )

    # Liouville and Choi matrices use different index pairings.  For column
    # vectorization vec(K)[input*5+output]=K[output,input], the exact reshuffle
    # below relates them.  In this special real rank-one POVM the two arrays
    # happen to coincide because every K_a is a symmetric outer product; their
    # tensor-index roles nevertheless remain distinct.
    lueders_choi = sum(
        (
            column_vectorize(kraus) * column_vectorize(kraus).H
            for kraus in lueders_kraus
        ),
        sp.zeros(25),
    )
    require(lueders_choi.H == lueders_choi, "Lueders Choi matrix is not Hermitian")
    require(lueders_choi.rank() == 6, "Lueders Choi rank changed")
    require(sp.trace(lueders_choi) == 5, "Lueders Choi trace changed")
    require(
        matrices_equal(partial_trace_first(lueders_choi, 5), identity5)
        and matrices_equal(partial_trace_second(lueders_choi, 5), identity5),
        "Lueders Choi partial traces changed",
    )
    reshuffle_ok = all(
        sp.simplify(
            lueders_choi[input_left * 5 + output_left, input_right * 5 + output_right]
            - lueders_liouville[output_right * 5 + output_left, input_right * 5 + input_left]
        )
        == 0
        for input_left in range(5)
        for output_left in range(5)
        for input_right in range(5)
        for output_right in range(5)
    )
    require(reshuffle_ok, "Liouville-to-Choi reshuffle failed")
    liouville_choi_array_coincidence = matrices_equal(lueders_choi, lueders_liouville)
    require(
        liouville_choi_array_coincidence,
        "real rank-one POVM lost its accidental Liouville-Choi array coincidence",
    )

    # A complex one-entry control distinguishes the required adjoint from an
    # incorrect transpose.  This directly guards the convention used above.
    control_a = sp.zeros(5)
    control_b = sp.zeros(5)
    control_a[0, 0] = sp.I
    control_b[0, 0] = 1
    control_direct = sp.simplify(
        (
            omega_orthogonal.T
            * sp.kronecker_product(control_a.conjugate(), control_b)
            * omega_orthogonal
        )[0]
        / 5
    )
    control_adjoint = sp.simplify(sp.trace(control_a.H * control_b) / 5)
    control_wrong_transpose = sp.simplify(sp.trace(control_a.T * control_b) / 5)
    require(control_direct == control_adjoint, "adjoint expectation control failed")
    require(control_direct != control_wrong_transpose, "transpose control did not fire")

    per_reflection = [
        {
            "index": index,
            "ovoid_permutation": permutation,
            "singlet_covariance": singlet_covariance[index],
            "effect_covariance": effect_covariance[index],
            "relational_operator_covariance": relational_covariance[index],
        }
        for index, permutation in enumerate(ovoid_permutations)
    ]

    return {
        "claim": "conditional canonical relational identity channel on conjugate(H5) tensor H5",
        "status": "PASS",
        "arithmetic": "exact algebraic numbers and rationals",
        "additional_contract": {
            "assumed_sector": "conjugate(H5) tensor H5",
            "assumed_action": "diagonal conjugate(native quartic G31 action) tensor native quartic G31 action",
            "derived_after_contract": [
                "one-dimensional invariant vector space and its singlet ray",
                "rank-one Choi operator of the identity channel",
                "maximally mixed marginals",
                "six-outcome correlated measurement law conditional on the singlet state",
                "relational spectral maximizer conditional on adopting its variational rule",
                "normalized square-root instrument conditional on adopting the Lueders prescription",
            ],
            "not_derived": [
                "selection, dynamical creation, or population of the singlet state",
                "a rule to maximize relational expectation or minimize a Hamiltonian",
                "coupling strength or transition rate",
                "identification with a P2 or D5 carrier",
                "selection of one individual ovoid outcome",
            ],
        },
        "sources": {
            "doily_quartic_bridge.py": file_hash(CONTRACT_DIR / "doily_quartic_bridge.py"),
            **bridge_report["sources"],
        },
        "isometric_bridge": {
            "primitive_intertwiner": matrix_json(intertwiner),
            "centered_ovoid_metric": matrix_json(centered_metric),
            "native_metric_for_exact_isometry": matrix_json(isometric_metric),
            "identity": "T^T (G/12) T = I5 + J5",
        },
        "singlet": {
            "invariant_subspace_dimension": singlet_dimension,
            "word_basis_ray_vec_G_inverse": [scalar_json(value) for value in omega_word],
            "word_basis_norm_squared": scalar_json(omega_word_norm),
            "rank": int(singlet_projector_word.rank()),
            "trace": scalar_json(singlet_projector_word.trace()),
            "orthonormal_density_matrix": matrix_json(singlet_density),
            "first_marginal": matrix_json(marginal_first),
            "second_marginal": matrix_json(marginal_second),
            "all_60_covariances": all(singlet_covariance),
            "selection_by_symmetry_alone": False,
            "invariant_density_counterexample": {
                "state": "I25/25",
                "first_marginal": matrix_json(identity5 / 5),
                "second_marginal": matrix_json(identity5 / 5),
                "joint_probability_matrix": matrix_json(mixed_joint_probabilities),
                "all_entries": "1/36",
                "relational_expectation": scalar_json(mixed_relational_expectation),
                "all_60_covariances": True,
            },
        },
        "channel": {
            "type": "identity channel on native H5 under the conjugate-input Choi identification",
            "choi_rank": int(choi_identity_channel.rank()),
            "choi_first_partial_trace": matrix_json(partial_trace_second(choi_identity_channel, 5)),
            "choi_second_partial_trace": matrix_json(partial_trace_first(choi_identity_channel, 5)),
            "completely_positive_witness": "Choi operator is five times a rank-one positive projector",
            "trace_preserving_witness": "input partial trace of the unnormalised Choi operator is I5",
        },
        "six_effect_correlations": {
            "effect_formula": "E_a = q_a q_a^T G / 12",
            "joint_probability_formula": "p(a,b)=Tr(E_a^dagger E_b)/5",
            "joint_probability_matrix": matrix_json(joint_probabilities),
            "same_outcome_probability": scalar_json(expected_same),
            "distinct_outcome_probability": scalar_json(expected_distinct),
            "row_marginals": [scalar_json(value) for value in row_marginals],
            "column_marginals": [scalar_json(value) for value in column_marginals],
            "conditional_same_outcome": scalar_json(conditional_same),
            "normalization": scalar_json(sum(joint_probabilities)),
            "all_60_effect_covariances": all(effect_covariance),
        },
        "expectation_convention": {
            "identity": "<Omega|conjugate(A) tensor B|Omega> = Tr(A^dagger B)/5",
            "direct_complex_control": str(control_direct),
            "adjoint_control": str(control_adjoint),
            "wrong_transpose_control": str(control_wrong_transpose),
            "transpose_control_fires": control_direct != control_wrong_transpose,
            "all_36_effect_pairs_direct_equal_trace_formula": direct_probabilities == joint_probabilities,
        },
        "non_scalar_invariant_relational_operator": {
            "formula": "R = sum_a conjugate(E_a) tensor E_a",
            "rank": int(relational_operator.rank()),
            "trace": scalar_json(relational_operator.trace()),
            "exact_spectrum": {"5/6": 1, "2/3": 5, "0": 19},
            "eigenspace_nullities": relational_nullities,
            "spectral_polynomial": "R (R-(5/6)I) (R-(2/3)I) = 0",
            "hermitian_in_isometric_chart": True,
            "unique_largest_eigenvalue": "5/6",
            "unique_largest_eigenray": "Omega",
            "singlet_expectation": scalar_json(singlet_relational_expectation),
            "maximally_mixed_expectation": scalar_json(mixed_relational_expectation),
            "conditional_variational_consequence": "If normalized states are required to maximize <R>, the unique maximizing ray is Omega",
            "variational_or_ground_state_rule_derived": False,
            "is_scalar": False,
            "all_60_covariances": all(relational_covariance),
        },
        "lueders_channel": {
            "instrument_prescription": "K_a = positive_sqrt(E_a) = sqrt(6/5) E_a",
            "kraus_count": len(lueders_kraus),
            "trace_preserving": True,
            "unital": True,
            "liouville_identity": "L = sum_a conjugate(K_a) tensor K_a = (6/5) R",
            "exact_liouville_spectrum": {"1": 1, "4/5": 5, "0": 19},
            "eigenspace_nullities": lueders_nullities,
            "spectral_polynomial": "L (L-I) (L-(4/5)I) = 0",
            "unique_stationary_density": matrix_json(identity5 / 5),
            "unique_stationary_operator_ray": "I5",
            "all_60_covariances": all(relational_covariance),
            "physical_time_evolution_derived": False,
        },
        "liouville_choi_distinction": {
            "R_role": "operator on conjugate(H5) tensor H5; after normalization (6/5)R is the Lueders Liouville matrix",
            "R_role_is_choi": False,
            "choi_formula": "J_L = sum_a vec(K_a) vec(K_a)^dagger",
            "column_vectorization": "vec(K)[input*5+output] = K[output,input]",
            "choi_rank": int(lueders_choi.rank()),
            "choi_trace": scalar_json(sp.trace(lueders_choi)),
            "choi_first_partial_trace": matrix_json(partial_trace_first(lueders_choi, 5)),
            "choi_second_partial_trace": matrix_json(partial_trace_second(lueders_choi, 5)),
            "same_numeric_array_after_reshuffle": liouville_choi_array_coincidence,
            "array_coincidence_reason": "each K_a is a real symmetric rank-one outer product",
            "semantic_roles_remain_distinct": True,
            "explicit_reshuffle_checked": reshuffle_ok,
        },
        "per_reflection": per_reflection,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    report = build_report()
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if arguments.output is None:
        sys.stdout.write(payload)
    else:
        arguments.output.write_text(payload, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
