#!/usr/bin/env python3
"""Exact classification of covariant instruments for the native six-simplex POVM.

The input is the already verified native H5 representation, its six rank-one
effects, and the sixty marked G31 actions from ``doily_quartic_bridge``.  No
extra carrier or dynamics is introduced.
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


def column_vectorize(matrix: sp.Matrix) -> sp.Matrix:
    rows, columns = matrix.shape
    return sp.Matrix(
        [
            matrix[output_index, input_index]
            for input_index in range(columns)
            for output_index in range(rows)
        ]
    )


def matrices_equal(left: sp.Matrix, right: sp.Matrix) -> bool:
    if left.shape != right.shape:
        return False
    return (left - right).applyfunc(sp.simplify) == sp.zeros(*left.shape)


def partial_trace_first(matrix: sp.Matrix, dimension: int) -> sp.Matrix:
    return sp.Matrix(
        dimension,
        dimension,
        lambda row, column: sum(
            matrix[index * dimension + row, index * dimension + column]
            for index in range(dimension)
        ),
    )


def partial_trace_second(matrix: sp.Matrix, dimension: int) -> sp.Matrix:
    return sp.Matrix(
        dimension,
        dimension,
        lambda row, column: sum(
            matrix[row * dimension + index, column * dimension + index]
            for index in range(dimension)
        ),
    )


def matrix_payload(matrix: sp.Matrix) -> str:
    return json.dumps(bridge.matrix_json(matrix), separators=(",", ":"), sort_keys=True)


def matrix_hash(matrix: sp.Matrix) -> str:
    return hashlib.sha256(matrix_payload(matrix).encode("utf-8")).hexdigest()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def permutation_matrix(permutation: list[int]) -> sp.Matrix:
    matrix = sp.zeros(6)
    for source, target in enumerate(permutation):
        matrix[target, source] = 1
    return matrix


def build_report() -> dict[str, Any]:
    bridge_report = bridge.build_report()
    require(bridge_report["status"] == "PASS", "doily quartic bridge did not pass")

    dimension = 5
    identity5 = sp.eye(dimension)
    metric = parse_matrix(
        bridge_report["selection_data"]["primitive_native_invariant_positive_form"]
    )
    native_rays = [
        sp.Matrix([parse_scalar(value) for value in ray])
        for ray in bridge_report["selection_data"]["native_ovoid_rays"]
    ]
    actions = [
        parse_matrix(row["quartic_numerator"]) / 16
        for row in bridge_report["quartic_bridge"]["per_reflection"]
    ]
    permutations = [
        list(row["ovoid_permutation"])
        for row in bridge_report["quartic_bridge"]["per_reflection"]
    ]
    require(len(native_rays) == 6 and len(actions) == 60, "native input counts changed")
    require(sum(native_rays, sp.zeros(5, 1)) == sp.zeros(5, 1), "simplex rays no longer sum to zero")

    gram = sp.Matrix(
        6,
        6,
        lambda row, column: (native_rays[row].T * metric * native_rays[column])[0],
    )
    expected_gram = 12 * sp.eye(6) - 2 * sp.ones(6)
    require(gram == expected_gram, "native six-simplex Gram matrix changed")

    effects = [ray * ray.T * metric / 12 for ray in native_rays]
    projectors = [sp.Rational(6, 5) * effect for effect in effects]
    require(sum(effects, sp.zeros(5)) == identity5, "native effects do not sum to identity")
    require(
        all(
            effect.rank() == 1
            and effect.trace() == sp.Rational(5, 6)
            and projector.rank() == 1
            and projector.trace() == 1
            and projector**2 == projector
            and metric.inv() * projector.T * metric == projector
            for effect, projector in zip(effects, projectors)
        ),
        "rank-one effect or projector normalization changed",
    )

    effect_covariance = []
    projector_covariance = []
    for action, permutation in zip(actions, permutations):
        inverse = action.inv()
        effect_covariance.append(
            all(
                action * effects[index] * inverse == effects[permutation[index]]
                for index in range(6)
            )
        )
        projector_covariance.append(
            all(
                action * projectors[index] * inverse == projectors[permutation[index]]
                for index in range(6)
            )
        )
    require(all(effect_covariance), "effect covariance failed")
    require(all(projector_covariance), "projector covariance failed")
    require(
        {permutation[0] for permutation in permutations} == set(range(6)),
        "the marked reflection actions are not transitive on the six outcomes",
    )

    # Exact stabilizer commutant.  The forty marked reflections fixing outcome
    # zero contain the ten transpositions of the remaining five outcomes.
    stabilizer_indices = [
        index for index, permutation in enumerate(permutations) if permutation[0] == 0
    ]
    stabilizer_permutations = {
        tuple(permutations[index]) for index in stabilizer_indices
    }
    require(len(stabilizer_indices) == 40, "marked stabilizer reflection count changed")
    require(len(stabilizer_permutations) == 10, "stabilizer transposition count changed")
    require(
        all(
            permutation[0] == 0
            and sum(permutation[index] != index for index in range(6)) == 2
            for permutation in stabilizer_permutations
        ),
        "outcome stabilizer no longer contains exactly the ten complementary transpositions",
    )
    unknowns = sp.symbols("x0:25")
    generic = sp.Matrix(5, 5, unknowns)
    equations: list[sp.Expr] = []
    for index in stabilizer_indices:
        equations.extend(list(actions[index] * generic - generic * actions[index]))
    commutant_matrix, _ = sp.linear_eq_to_matrix(equations, unknowns)
    commutant_dimension = 25 - commutant_matrix.rank()
    require(commutant_dimension == 2, "outcome stabilizer commutant dimension changed")
    require(
        all(actions[index] * projectors[0] == projectors[0] * actions[index] for index in stabilizer_indices),
        "distinguished projector is not in the stabilizer commutant",
    )
    commutant_generators = sp.Matrix.hstack(
        column_vectorize(identity5), column_vectorize(projectors[0])
    )
    require(commutant_generators.rank() == 2, "I and P_0 are not independent")

    parameter = sp.Symbol("p", real=True)
    output_states = [
        sp.simplify(
            parameter * projector
            + (1 - parameter) * (identity5 - projector) / 4
        )
        for projector in projectors
    ]
    require(all(sp.trace(state) == 1 for state in output_states), "output-state trace changed")
    require(
        all(
            metric.inv() * state.T * metric == state
            for state in output_states
        ),
        "output-state self-adjointness changed",
    )
    require(
        all(
            action * output_states[index] * action.inv()
            == output_states[permutation[index]]
            for action, permutation in zip(actions, permutations)
            for index in range(6)
        ),
        "one-parameter output-state family is not fully covariant",
    )

    # The two orthogonal blocks have eigenvalues p and (1-p)/4.  Positivity is
    # therefore exactly the closed interval 0 <= p <= 1.
    require(
        all(
            matrices_equal(state * projector, parameter * projector)
            and matrices_equal(
                state * (identity5 - projector),
                ((1 - parameter) / 4) * (identity5 - projector),
            )
            for state, projector in zip(output_states, projectors)
        ),
        "output-state block eigenvalues changed",
    )

    # Work in the exact isometric chart to compare the induced channel with
    # the previously verified Lueders channel and read its full spectrum.
    square_root_metric = sp.diag(1, sp.sqrt(3), sp.sqrt(3), sp.sqrt(3), sp.sqrt(6))
    require(square_root_metric.T * square_root_metric == metric, "metric square root changed")
    orthogonal_effects = [
        sp.simplify(square_root_metric * effect * square_root_metric.inv())
        for effect in effects
    ]
    effect_span_dimension = sp.Matrix.hstack(
        *(column_vectorize(effect) for effect in orthogonal_effects)
    ).rank()
    require(effect_span_dimension == 6, "six-effect operator span dimension changed")

    # Exact tomography counterexample.  Every effect is real symmetric, hence
    # it is blind to the Hermitian imaginary antisymmetric coherence X below.
    # The two perturbed densities are distinct and strictly positive.
    blind_coherence = sp.zeros(5)
    blind_coherence[0, 1] = sp.I
    blind_coherence[1, 0] = -sp.I
    require(
        blind_coherence.H == blind_coherence and sp.trace(blind_coherence) == 0,
        "blind coherence is not traceless Hermitian",
    )
    rho_plus = identity5 / 5 + blind_coherence / 10
    rho_minus = identity5 / 5 - blind_coherence / 10
    expected_density_spectrum = {
        sp.Rational(1, 10): 1,
        sp.Rational(1, 5): 3,
        sp.Rational(3, 10): 1,
    }
    require(rho_plus != rho_minus, "tomography counterexample densities coincide")
    require(
        rho_plus.H == rho_plus
        and rho_minus.H == rho_minus
        and sp.trace(rho_plus) == 1
        and sp.trace(rho_minus) == 1,
        "tomography counterexample is not a pair of density matrices",
    )
    require(
        rho_plus.eigenvals() == expected_density_spectrum
        and rho_minus.eigenvals() == expected_density_spectrum,
        "tomography counterexample density spectrum changed",
    )
    blind_effect_pairings = [
        sp.simplify(sp.trace(effect * blind_coherence))
        for effect in orthogonal_effects
    ]
    require(
        blind_effect_pairings == [0] * 6,
        "six-effect POVM detected the blind coherence",
    )
    first_distribution_plus = [
        sp.simplify(sp.trace(effect * rho_plus)) for effect in orthogonal_effects
    ]
    first_distribution_minus = [
        sp.simplify(sp.trace(effect * rho_minus)) for effect in orthogonal_effects
    ]
    require(
        first_distribution_plus == [sp.Rational(1, 6)] * 6
        and first_distribution_minus == first_distribution_plus,
        "tomography counterexample first-outcome distributions differ",
    )
    orthogonal_projectors = [sp.Rational(6, 5) * effect for effect in orthogonal_effects]
    orthogonal_states = [
        sp.simplify(
            parameter * projector
            + (1 - parameter) * (identity5 - projector) / 4
        )
        for projector in orthogonal_projectors
    ]
    lueders_kraus = [sp.sqrt(sp.Rational(6, 5)) * effect for effect in orthogonal_effects]
    lueders_liouville = sum(
        (sp.kronecker_product(kraus.conjugate(), kraus) for kraus in lueders_kraus),
        sp.zeros(25),
    )
    instrument_channel = sum(
        (
            column_vectorize(state) * column_vectorize(effect).T
            for state, effect in zip(orthogonal_states, orthogonal_effects)
        ),
        sp.zeros(25),
    )
    identity_vector = column_vectorize(identity5)
    trace_replace_identity = identity_vector * identity_vector.T
    lueders_coefficient = (5 * parameter - 1) / 4
    replacement_coefficient = (1 - parameter) / 4
    require(
        matrices_equal(
            instrument_channel,
            lueders_coefficient * lueders_liouville
            + replacement_coefficient * trace_replace_identity,
        ),
        "classified channel formula failed",
    )

    identity_projector = identity_vector * identity_vector.T / 5
    standard_projector = sp.Rational(5, 4) * (lueders_liouville - identity_projector)
    residual_projector = sp.eye(25) - identity_projector - standard_projector
    require(
        identity_projector.rank() == 1
        and standard_projector.rank() == 5
        and residual_projector.rank() == 19,
        "operator-space projector ranks changed",
    )
    require(
        all(
            projector**2 == projector and projector.H == projector
            for projector in (identity_projector, standard_projector, residual_projector)
        ),
        "operator-space spectral projectors failed",
    )
    require(
        all(
            left * right == sp.zeros(25)
            for left, right in (
                (identity_projector, standard_projector),
                (standard_projector, identity_projector),
                (identity_projector, residual_projector),
                (residual_projector, identity_projector),
                (standard_projector, residual_projector),
                (residual_projector, standard_projector),
            )
        ),
        "operator-space spectral projectors are not orthogonal",
    )
    channel_eigenvalue = parameter - sp.Rational(1, 5)
    require(
        matrices_equal(
            instrument_channel,
            identity_projector + channel_eigenvalue * standard_projector,
        ),
        "classified channel spectral decomposition failed",
    )

    # Choi matrix of the channel after forgetting the classical outcome.
    # The replacement map rho -> Tr(rho) I has Choi matrix I25.  For this
    # real rank-one POVM the Lueders Choi array equals its Liouville array,
    # although the two index roles remain distinct.
    lueders_choi = sum(
        (
            column_vectorize(kraus) * column_vectorize(kraus).H
            for kraus in lueders_kraus
        ),
        sp.zeros(25),
    )
    require(
        matrices_equal(lueders_choi, lueders_liouville),
        "special Lueders Choi-Liouville array equality changed",
    )
    forgotten_outcome_choi = sum(
        (
            sp.kronecker_product(effect.T, state)
            for effect, state in zip(orthogonal_effects, orthogonal_states)
        ),
        sp.zeros(25),
    )
    require(
        matrices_equal(
            forgotten_outcome_choi,
            lueders_coefficient * lueders_choi
            + replacement_coefficient * sp.eye(25),
        ),
        "forgotten-outcome Choi formula failed",
    )
    choi_top_eigenvalue = parameter
    choi_standard_eigenvalue = (15 * parameter + 1) / 20
    choi_residual_eigenvalue = (1 - parameter) / 4
    require(
        matrices_equal(
            forgotten_outcome_choi,
            choi_top_eigenvalue * identity_projector
            + choi_standard_eigenvalue * standard_projector
            + choi_residual_eigenvalue * residual_projector,
        ),
        "forgotten-outcome Choi spectral decomposition failed",
    )

    # Consecutive outcomes form the exact six-state Markov chain obtained by
    # measuring E_b on the prepared state sigma_a.
    transition = sp.Matrix(
        6,
        6,
        lambda row, column: sp.simplify(
            sp.trace(orthogonal_effects[column] * orthogonal_states[row])
        ),
    )
    transition_diagonal = sp.Rational(5, 6) * parameter
    transition_off_diagonal = (6 - 5 * parameter) / 30
    require(
        all(transition[index, index] == transition_diagonal for index in range(6)),
        "repeated-outcome diagonal probability changed",
    )
    require(
        all(
            transition[row, column] == transition_off_diagonal
            for row in range(6)
            for column in range(6)
            if row != column
        ),
        "repeated-outcome off-diagonal probability changed",
    )
    require(
        all(sum(transition[row, column] for column in range(6)) == 1 for row in range(6)),
        "repeated-outcome transition rows do not sum to one",
    )
    stationary_joint = transition / 6
    require(
        stationary_joint.subs(parameter, 1)[0, 0] == sp.Rational(5, 36)
        and stationary_joint.subs(parameter, 1)[0, 1] == sp.Rational(1, 180),
        "p=1 stationary joint law changed",
    )
    # For a word a_1...a_n the measure-and-prepare form gives the first
    # probability followed by the same transition matrix at every later
    # position.  Equal first distributions therefore prove equality for all
    # finite words by induction, without a bounded word enumeration.
    all_finite_word_laws_equal = (
        first_distribution_plus == first_distribution_minus
        and all(
            sum(transition[row, column] for column in range(6)) == 1
            for row in range(6)
        )
    )
    require(
        all_finite_word_laws_equal,
        "analytic all-word tomography counterexample failed",
    )

    # Rational six-coordinate realization of the same H5.  It is the
    # zero-sum subspace in Q^6 generated by the already present simplex rays;
    # its Gram matrix agrees with the native one divided by twelve.
    identity6 = sp.eye(6)
    ones6 = sp.ones(6, 1)
    zero_sum_projector = identity6 - ones6 * ones6.T / 6
    ambient_rays = [identity6[:, index] - ones6 / 6 for index in range(6)]
    ambient_gram = sp.Matrix(
        6,
        6,
        lambda row, column: (ambient_rays[row].T * ambient_rays[column])[0],
    )
    require(ambient_gram == gram / 12, "rational simplex frame is not the native H5 frame")
    ambient_effects = [ray * ray.T for ray in ambient_rays]
    ambient_projectors = [sp.Rational(6, 5) * effect for effect in ambient_effects]
    require(
        sum(ambient_effects, sp.zeros(6)) == zero_sum_projector,
        "ambient effects do not resolve the H5 identity",
    )
    require(
        all(
            projector**2 == projector
            and projector.T == projector
            and projector.rank() == 1
            for projector in ambient_projectors
        ),
        "ambient simplex projector failed",
    )
    require(
        all(
            (zero_sum_projector - projector).T == zero_sum_projector - projector
            and (zero_sum_projector - projector) ** 2 == zero_sum_projector - projector
            and (zero_sum_projector - projector).rank() == 4
            for projector in ambient_projectors
        ),
        "ambient complementary projector failed",
    )
    require(
        all(
            permutation_matrix(permutation) * ambient_projectors[index]
            * permutation_matrix(permutation).T
            == ambient_projectors[permutation[index]]
            for permutation in permutations
            for index in range(6)
        ),
        "rational simplex frame covariance failed",
    )

    ambient_states = [
        sp.simplify(
            parameter * projector
            + (1 - parameter) * (zero_sum_projector - projector) / 4
        )
        for projector in ambient_projectors
    ]
    endpoint_data: dict[str, Any] = {}
    for endpoint in (0, 1):
        outcome_choi = [
            sp.kronecker_product(
                effect.T,
                state.subs(parameter, endpoint),
            )
            for effect, state in zip(ambient_effects, ambient_states)
        ]
        total_choi = sum(outcome_choi, sp.zeros(36))
        expected_outcome_rank = 4 if endpoint == 0 else 1
        require(
            all(matrix.T == matrix for matrix in outcome_choi),
            f"p={endpoint} outcome Choi symmetry failed",
        )
        require(
            all(matrix.rank() == expected_outcome_rank for matrix in outcome_choi),
            f"p={endpoint} minimal Kraus ranks changed",
        )
        require(
            all(all(entry.is_Rational for entry in matrix) for matrix in outcome_choi),
            f"p={endpoint} outcome Choi matrix is not rational",
        )
        require(
            partial_trace_second(total_choi, 6) == zero_sum_projector,
            f"p={endpoint} total Choi input marginal changed",
        )
        require(
            partial_trace_first(total_choi, 6) == zero_sum_projector,
            f"p={endpoint} total Choi output marginal changed",
        )
        endpoint_data[str(endpoint)] = {
            "all_six_outcome_choi_entries_rational": True,
            "outcome_choi_rank": expected_outcome_rank,
            "total_choi_rank": int(total_choi.rank()),
            "total_choi_trace": bridge.rational_json(sp.trace(total_choi)),
            "total_choi_sha256": matrix_hash(total_choi),
            "input_partial_trace_is_H5_identity": True,
            "output_partial_trace_is_H5_identity": True,
        }

    return {
        "claim": "complete finite G31-covariant instrument classification for the native rank-one six-simplex POVM",
        "status": "PASS",
        "arithmetic": "exact rationals and algebraic numbers; no floating point",
        "sources": {
            "doily_quartic_bridge.py": file_hash(CONTRACT_DIR / "doily_quartic_bridge.py"),
            **bridge_report["sources"],
        },
        "native_input": {
            "hilbert_dimension": 5,
            "outcome_count": 6,
            "effect_rank": 1,
            "effect_trace": "5/6",
            "all_60_effect_covariances": all(effect_covariance),
            "transitive_outcome_action": True,
        },
        "classification": {
            "instrument_formula": "I_a(rho) = Tr(E_a rho) sigma_a",
            "rank_one_effect_implication": "every CP instrument operation with effect E_a is measure-and-prepare",
            "stabilizer_marked_reflections": len(stabilizer_indices),
            "stabilizer_distinct_transpositions": len(stabilizer_permutations),
            "stabilizer_commutant_dimension_over_Q": commutant_dimension,
            "commutant_basis": ["I5", "P_a"],
            "state_formula": "sigma_a(p) = p P_a + (1-p)(I5-P_a)/4",
            "projector_formula": "P_a = (6/5) E_a",
            "complete_positive_parameter_interval": "0 <= p <= 1",
            "state_block_eigenvalues": {"P_a_rank_1": "p", "orthogonal_rank_4": "(1-p)/4"},
            "minimal_kraus_rank": {
                "p=1": 1,
                "p=0": 4,
                "0<p<1": 5,
            },
            "tagged_instrument_total_kraus_rank": {
                "p=1": 6,
                "p=0": 24,
                "0<p<1": 30,
            },
            "one_kraus_or_pure_output_selects": "p=1",
            "one_kraus_or_pure_output_is_additional_assumption": True,
            "minimal_tagged_environment_selects": "p=1",
            "minimal_tagged_environment_is_additional_assumption": True,
            "primitive_P1_implies_one_kraus": False,
        },
        "channel": {
            "formula": "Phi_p = ((5p-1)/4) Phi_Lueders + ((1-p)/4) Tr(rho) I5",
            "exact_spectrum_on_End_H5": {"1": 1, "p-1/5": 5, "0": 19},
            "spectral_subspace_ranks": {
                "identity": int(identity_projector.rank()),
                "standard": int(standard_projector.rank()),
                "residual": int(residual_projector.rank()),
            },
            "unique_stationary_density_for_0_le_p_le_1": bridge.matrix_json(identity5 / 5),
            "stationarity_reason": "p-1/5 lies in [-1/5,4/5], so eigenvalue 1 is simple",
            "p=1_is_Lueders": True,
            "forgotten_outcome_choi_formula": "J_p = ((5p-1)/4) J_Lueders + ((1-p)/4) I25",
            "forgotten_outcome_choi_spectrum": {
                "p": 1,
                "(15p+1)/20": 5,
                "(1-p)/4": 19,
            },
            "forgotten_outcome_choi_rank": {
                "p=1": 6,
                "p=0": 24,
                "0<p<1": 25,
            },
            "tagged_and_forgotten_kraus_ranks_are_distinct_notions": True,
        },
        "repeated_outcomes": {
            "markov_state_count": 6,
            "transition_diagonal": "5p/6",
            "transition_off_diagonal": "(6-5p)/30",
            "uniform_stationary_distribution": ["1/6"] * 6,
            "stationary_joint_diagonal": "5p/36",
            "stationary_joint_off_diagonal": "(6-5p)/180",
            "p=1_joint_diagonal": "5/36",
            "p=1_joint_off_diagonal": "1/180",
            "p=1_matches_Omega_two_site_numbers": True,
            "temporal_match_proves_entanglement": False,
        },
        "tomography_counterexample": {
            "effect_operator_span_dimension": effect_span_dimension,
            "full_operator_space_dimension": 25,
            "blind_coherence": bridge.matrix_json(blind_coherence),
            "rho_plus": bridge.matrix_json(rho_plus),
            "rho_minus": bridge.matrix_json(rho_minus),
            "both_density_spectra": {
                "1/10": 1,
                "1/5": 3,
                "3/10": 1,
            },
            "all_six_Tr_Ea_X": [bridge.rational_json(value) for value in blind_effect_pairings],
            "both_first_outcome_distributions": ["1/6"] * 6,
            "finite_word_probability_formula": "Pr(a1...an|rho) = Tr(E_a1 rho) product_{k=2}^n T[a_{k-1},a_k]",
            "all_finite_outcome_word_laws_equal_for_rho_plus_minus": all_finite_word_laws_equal,
            "complete_state_tomography": False,
            "complete_GNS_crossword_reconstruction": False,
        },
        "rational_endpoint_choi": {
            "frame": "the same native H5 represented as the zero-sum subspace of Q^6",
            "native_to_frame_gram_identity": "Gram(u_a) = Gram(q_a)/12",
            "choi_formula": "J(I_a) = E_a^T tensor sigma_a",
            "positivity_certificate": "E_a, P_a, and Q-P_a are positive orthogonal projectors up to positive rational scale",
            "kraus_formulas": {
                "p=1": "M_a = sqrt(6/5) u_a u_a^T",
                "p=0": "M_{a,r} = (1/2) v_{a,r} u_a^T for any orthonormal basis v_{a,r} of u_a^perp intersect H5",
            },
            "endpoints": endpoint_data,
        },
        "scope": {
            "complete_within": "CP instruments on the native H5 with the fixed six effects and their actual G31 outcome covariance",
            "does_not_derive": [
                "which value of p is physically realized",
                "a one-Kraus or pure-output principle from P1 primitiveness",
                "physical time, rates, spacetime, or an entanglement source",
                "complete state tomography or reconstruction of unmeasured GNS coherences",
            ],
        },
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
