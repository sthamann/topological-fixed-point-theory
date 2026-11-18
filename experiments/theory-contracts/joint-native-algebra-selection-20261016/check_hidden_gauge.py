#!/usr/bin/env python3
"""Exact/supplemental checks for the central-sector Pauli gauge freedom.

This checker does not select a physical source law.  It asks a narrower
question: once arbitrary unitary elements of the represented source C*-algebra
are admitted as Kraus events, does the central J-character permit a channel
which is invisible on H4 and H5, preserves the H5 -> H5 tensor I4/4 erasure
instrument branch by branch, and remains visible to H20 gauge correlators?
"""

from __future__ import annotations

import argparse
import hashlib
import itertools as it
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp


REPO = Path("/Users/stefanhamann/Projekte/tfpt-theoryv4")
DEFAULT_OUT = Path(__file__).with_name("hidden_gauge_results.json")
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from tfpt_explorer import process, quartic_selection  # noqa: E402


def kron_all(*matrices: Any) -> Any:
    result = matrices[0]
    for matrix in matrices[1:]:
        result = np.kron(result, matrix)
    return result


def sp_kron_all(*matrices: sp.Matrix) -> sp.Matrix:
    result = matrices[0]
    for matrix in matrices[1:]:
        result = sp.kronecker_product(result, matrix)
    return result


def block_diag(*blocks: sp.Matrix) -> sp.Matrix:
    return sp.diag(*blocks)


def matrix_unit(n: int, row: int, column: int) -> sp.Matrix:
    result = sp.zeros(n)
    result[row, column] = 1
    return result


def max_abs(matrix: np.ndarray) -> float:
    return float(np.max(np.abs(matrix))) if matrix.size else 0.0


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def integer_matrix_sha256(matrix: sp.Matrix) -> str:
    """Hash an integer matrix in a small language-independent text format."""

    payload = "\n".join(
        ",".join(str(int(matrix[row, column])) for column in range(matrix.cols))
        for row in range(matrix.rows)
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def independent_rows_mod_prime(matrix: sp.Matrix, prime: int) -> list[int]:
    """Return the first independent row indices over F_prime."""

    echelon: dict[int, list[int]] = {}
    selected: list[int] = []
    for row_index in range(matrix.rows):
        row = [int(matrix[row_index, column]) % prime for column in range(matrix.cols)]
        for pivot in sorted(echelon):
            if row[pivot]:
                factor = row[pivot]
                basis = echelon[pivot]
                row = [(entry - factor * basis_entry) % prime for entry, basis_entry in zip(row, basis)]
        if not any(row):
            continue
        pivot = next(column for column, entry in enumerate(row) if entry)
        inverse = pow(row[pivot], -1, prime)
        row = [(inverse * entry) % prime for entry in row]
        echelon[pivot] = row
        selected.append(row_index)
        if len(selected) == matrix.cols:
            break
    return selected


def determinant_mod_prime(matrix: sp.Matrix, prime: int) -> int:
    """Bare modular Gaussian determinant, used for the replay certificate."""

    if matrix.rows != matrix.cols:
        raise ValueError("determinant requires a square matrix")
    rows = [
        [int(matrix[row, column]) % prime for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]
    determinant = 1
    for column in range(matrix.cols):
        pivot = next((row for row in range(column, matrix.rows) if rows[row][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            determinant = -determinant
        pivot_value = rows[column][column]
        determinant = determinant * pivot_value % prime
        inverse = pow(pivot_value, -1, prime)
        rows[column] = [(inverse * entry) % prime for entry in rows[column]]
        for row in range(column + 1, matrix.rows):
            if rows[row][column]:
                factor = rows[row][column]
                rows[row] = [
                    (entry - factor * pivot_entry) % prime
                    for entry, pivot_entry in zip(rows[row], rows[column])
                ]
    return determinant % prime


GaussianInteger = tuple[int, int]


def gaussian_add(left: GaussianInteger, right: GaussianInteger) -> GaussianInteger:
    return left[0] + right[0], left[1] + right[1]


def gaussian_multiply(left: GaussianInteger, right: GaussianInteger) -> GaussianInteger:
    return left[0] * right[0] - left[1] * right[1], left[0] * right[1] + left[1] * right[0]


def exact_quartic_code_action(
    twice_reflection: list[list[GaussianInteger]],
    word_sets: list[list[tuple[int, ...]]],
) -> tuple[sp.Matrix, bool]:
    """Compute (2r)^tensor4 B=B A_num in exact Gaussian-integer arithmetic."""

    support_to_column = {
        word: column
        for column, words in enumerate(word_sets)
        for word in words
    }
    if len(support_to_column) != sum(len(words) for words in word_sets):
        raise ArithmeticError("quartic word-basis supports are not disjoint")
    values: dict[tuple[int, ...], list[GaussianInteger]] = {}
    for output_word in it.product(range(4), repeat=4):
        row_values: list[GaussianInteger] = []
        for words in word_sets:
            total = (0, 0)
            for input_word in words:
                product = (1, 0)
                for output_index, input_index in zip(output_word, input_word):
                    product = gaussian_multiply(product, twice_reflection[output_index][input_index])
                total = gaussian_add(total, product)
            row_values.append(total)
        values[output_word] = row_values

    numerator = sp.zeros(5)
    for row, words in enumerate(word_sets):
        representative = words[0]
        for column, (real, imaginary) in enumerate(values[representative]):
            if imaginary != 0:
                raise ArithmeticError("quartic code action has a non-real representative coefficient")
            numerator[row, column] = real

    exact_carrier = True
    for output_word, row_values in values.items():
        support_column = support_to_column.get(output_word)
        for column, value in enumerate(row_values):
            expected = (0, 0) if support_column is None else (int(numerator[support_column, column]), 0)
            exact_carrier &= value == expected
    return numerator, exact_carrier


def sympy_paulis() -> list[sp.Matrix]:
    i2 = sp.eye(2)
    x = sp.Matrix([[0, 1], [1, 0]])
    y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    z = sp.diag(1, -1)
    return [sp.kronecker_product(left, right) for left, right in it.product((i2, x, y, z), repeat=2)]


def numpy_paulis() -> list[np.ndarray]:
    identity = np.eye(4, dtype=complex)
    native, _, _ = process._pauli_data()
    return [identity, *[np.asarray(matrix, dtype=complex) for matrix in native]]


def pauli_match(matrix: np.ndarray, paulis: list[np.ndarray]) -> tuple[int, complex, float]:
    """Return the closest projective Pauli and its Hilbert-Schmidt phase."""

    best: tuple[int, complex, float] | None = None
    for index, pauli in enumerate(paulis):
        overlap = np.trace(pauli.conj().T @ matrix) / 4
        if abs(overlap) < 1e-14:
            residual = float(np.linalg.norm(matrix - pauli))
            phase = 0.0j
        else:
            phase = overlap / abs(overlap)
            residual = float(np.linalg.norm(matrix - phase * pauli))
        candidate = (index, complex(phase), residual)
        if best is None or candidate[2] < best[2]:
            best = candidate
    assert best is not None
    return best


def erasure(matrix5: np.ndarray) -> np.ndarray:
    return np.kron(matrix5, np.eye(4) / 4)


def build_report() -> dict[str, Any]:
    checks: list[dict[str, Any]] = []

    def check(name: str, ok: bool, kind: str, actual: Any, expected: Any) -> None:
        checks.append({"name": name, "kind": kind, "ok": bool(ok), "actual": actual, "expected": expected})

    # Actual native constructors used by the existing explorer.
    v_code, word_sets = process._code_basis()
    v_code = np.asarray(v_code, dtype=complex)
    paulis_np = numpy_paulis()
    rays = process._hamming_rays()
    native_reflections = []
    native_reflections_exact = []
    native_twice_reflections_gaussian: list[list[list[GaussianInteger]]] = []
    exact_ray_denominators = []
    exact_reflection_residuals = []
    exact_twice_reflection_flags = []
    native_logical = []
    native_reduce_residuals = []
    for ray in rays:
        ray = np.asarray(ray, dtype=complex)
        exact_coordinates = []
        gaussian_coordinates: list[GaussianInteger] = []
        for coordinate in ray:
            real = int(coordinate.real)
            imaginary = int(coordinate.imag)
            if coordinate.real != real or coordinate.imag != imaginary:
                raise ArithmeticError("native Hamming ray does not have Gaussian-integer coordinates")
            exact_coordinates.append(sp.Integer(real) + sp.I * sp.Integer(imaginary))
            gaussian_coordinates.append((real, imaginary))
        exact_ray = sp.Matrix(exact_coordinates)
        exact_denominator = sp.expand((exact_ray.H * exact_ray)[0])
        exact_reflection = sp.eye(4) - 2 * exact_ray * exact_ray.H / exact_denominator
        twice_reflection_gaussian: list[list[GaussianInteger]] = []
        for row in range(4):
            twice_row: list[GaussianInteger] = []
            for column in range(4):
                outer = gaussian_multiply(
                    gaussian_coordinates[row],
                    (gaussian_coordinates[column][0], -gaussian_coordinates[column][1]),
                )
                twice_row.append((2 * int(row == column) - outer[0], -outer[1]))
            twice_reflection_gaussian.append(twice_row)
        twice_reflection_symbolic = sp.Matrix(
            [[sp.Integer(real) + sp.I * sp.Integer(imaginary) for real, imaginary in row] for row in twice_reflection_gaussian]
        )
        projector = np.outer(ray, ray.conj()) / np.vdot(ray, ray)
        reflection = np.eye(4, dtype=complex) - 2 * projector
        exact_reflection_np = np.asarray(exact_reflection.tolist(), dtype=complex)
        reflection4 = kron_all(reflection, reflection, reflection, reflection)
        logical = v_code.conj().T @ reflection4 @ v_code
        native_reflections.append(reflection)
        native_reflections_exact.append(exact_reflection)
        native_twice_reflections_gaussian.append(twice_reflection_gaussian)
        exact_twice_reflection_flags.append(
            (twice_reflection_symbolic - 2 * exact_reflection).applyfunc(sp.simplify) == sp.zeros(4)
        )
        exact_ray_denominators.append(exact_denominator)
        exact_reflection_residuals.append(float(np.linalg.norm(reflection - exact_reflection_np)))
        native_logical.append(logical)
        native_reduce_residuals.append(float(np.linalg.norm(reflection4 @ v_code - v_code @ logical)))

    check(
        "all sixty source reflections are reconstructed exactly from Gaussian-integer rays",
        len(native_reflections_exact) == 60
        and set(exact_ray_denominators) == {sp.Integer(4)}
        and all(exact_twice_reflection_flags)
        and max(exact_reflection_residuals) < 1e-14,
        "exact_gaussian_integer_source_with_supplemental_matrix_identity",
        {
            "reflection_count": len(native_reflections_exact),
            "ray_norms_squared": sorted(str(value) for value in set(exact_ray_denominators)),
            "max_float_constructor_residual": max(exact_reflection_residuals),
        },
        {"reflection_count": 60, "ray_norms_squared": ["4"], "max_float_constructor_residual": 0.0},
    )

    # Pauli invariance of the actual code, and the essential Sym^4 scope.
    pauli_quartic_residuals = []
    for pauli in paulis_np:
        pauli4 = kron_all(pauli, pauli, pauli, pauli)
        pauli_quartic_residuals.append(float(np.linalg.norm(pauli4 @ v_code - v_code)))
    tensor = v_code.reshape(4, 4, 4, 4, 5)
    symmetry_residual = max(
        float(np.linalg.norm(np.transpose(tensor, (*permutation, 4)) - tensor))
        for permutation in it.permutations(range(4))
    )
    quartic_data = quartic_selection.build_quartic_selection_data()
    degree4_dimension = next(
        row["invariant_dimension"]
        for row in quartic_data["data"]["molien"]["coefficients"]
        if row["degree"] == 4
    )
    # Exact character average on the full tensor fourth power: only I has nonzero trace.
    full_tensor_fixed_dimension = (4**4 + 15 * 0**4) // 16
    check(
        "actual V5 is fixed by all 16 projective Pauli words",
        max(pauli_quartic_residuals) < 1e-12,
        "supplemental_matrix_residual",
        max(pauli_quartic_residuals),
        0.0,
    )
    check(
        "actual V5 lies in Sym^4(V4)",
        symmetry_residual < 1e-12,
        "supplemental_matrix_residual",
        symmetry_residual,
        0.0,
    )
    check(
        "Pauli fixed dimension inside Sym^4 is exactly five",
        quartic_data["status"] == "PASS" and degree4_dimension == 5 and len(word_sets) == 5,
        "exact_existing_molien_and_constructor",
        {"molien_degree4": degree4_dimension, "code_columns": len(word_sets)},
        {"molien_degree4": 5, "code_columns": 5},
    )
    check(
        "Pauli fixed dimension on the full V4 tensor fourth power is sixteen",
        full_tensor_fixed_dimension == 16,
        "exact_negative_control",
        full_tensor_fixed_dimension,
        16,
    )

    # Exact 29-dimensional represented source algebra: H4 + H5 + H20.
    paulis_sp = sympy_paulis()
    i4, i5, i20, i29 = sp.eye(4), sp.eye(5), sp.eye(20), sp.eye(29)
    j29 = block_diag(sp.I * i4, i5, -sp.I * i20)
    e20 = sum((sp.I**power * (j29**power) for power in range(4)), sp.zeros(29)) / 4
    target_e20 = block_diag(sp.zeros(4), sp.zeros(5), i20)
    check(
        "central quarter-turn polynomial is exactly the H20 character projector",
        e20 == target_e20 and e20**2 == e20 and e20.H == e20,
        "exact_symbolic",
        {"rank": int(e20.rank()), "block_ranks": [int(e20[:4, :4].rank()), int(e20[4:9, 4:9].rank()), int(e20[9:, 9:].rank())]},
        {"rank": 20, "block_ranks": [0, 0, 20]},
    )
    check(
        "J acts as i, 1, and -i on H4, H5, and H20",
        j29[:4, :4] == sp.I * i4 and j29[4:9, 4:9] == i5 and j29[9:, 9:] == -sp.I * i20,
        "exact_symbolic",
        ["i I4", "I5", "-i I20"],
        ["i I4", "I5", "-i I20"],
    )

    source_paulis = [block_diag(pauli, i5, sp.kronecker_product(i5, pauli.conjugate())) for pauli in paulis_sp]
    gated = [(i29 - e20) + e20 * source_pauli for source_pauli in source_paulis]
    unitary_flags = [unitary.H * unitary == i29 and unitary * unitary.H == i29 for unitary in gated]
    check(
        "all sixteen centrally gated Pauli elements are unitary",
        all(unitary_flags),
        "exact_symbolic",
        sum(unitary_flags),
        16,
    )
    check(
        "gated Pauli elements are identity on H4 and H5 and Pauli on H20",
        all(
            unitary[:4, :4] == i4
            and unitary[4:9, 4:9] == i5
            and unitary[9:, 9:] == sp.kronecker_product(i5, pauli.conjugate())
            for unitary, pauli in zip(gated, paulis_sp)
        ),
        "exact_symbolic",
        True,
        True,
    )

    q_exact = sp.Rational(1, 3)
    weights = [
        (1 + 15 * q_exact) / 16,
        *[(1 - q_exact) / 16] * 15,
    ]
    normalization = sum((weight * unitary.H * unitary for weight, unitary in zip(weights, gated)), sp.zeros(29))
    check(
        "common marked records are normalized in the complete source algebra",
        normalization == i29 and sum(weights) == 1,
        "exact_symbolic",
        {"weight_sum": str(sum(weights)), "w0": str(weights[0]), "w_nonzero": str(weights[1])},
        {"weight_sum": "1", "w0": "3/8", "w_nonzero": "1/24"},
    )

    # Exact branchwise erasure and exact depolarizing formula on a matrix-unit basis.
    gauge_mixed = i4 / 4
    branchwise_gauge_flags = [pauli.conjugate() * gauge_mixed * pauli.T == gauge_mixed for pauli in paulis_sp]
    check(
        "every marked Pauli record preserves the erasure gauge state branchwise",
        all(branchwise_gauge_flags),
        "exact_symbolic",
        sum(branchwise_gauge_flags),
        16,
    )
    depolarizing_flags = []
    for row, column in it.product(range(4), repeat=2):
        unit = matrix_unit(4, row, column)
        actual = sum(
            (weight * pauli * unit * pauli.H for weight, pauli in zip(weights, paulis_sp)),
            sp.zeros(4),
        )
        expected = q_exact * unit + (1 - q_exact) * sp.trace(unit) * i4 / 4
        depolarizing_flags.append(actual == expected)
    check(
        "the coarse H20 gauge channel is exactly D_q on all matrix units",
        all(depolarizing_flags),
        "exact_symbolic",
        sum(depolarizing_flags),
        16,
    )

    # A sequence containing both an actual native reflection and a marked gauge event.
    reflection = native_reflections[0]
    logical = native_logical[0]
    gauge = paulis_np[5]
    reflection20 = np.kron(logical, reflection.conj())
    gauge20 = np.kron(np.eye(5), gauge.conj())
    probe = np.asarray([1, 1j, -1, 2, -2j], dtype=complex)
    probe /= np.linalg.norm(probe)
    rho5 = np.outer(probe, probe.conj())
    weight = float(weights[5])
    start20 = erasure(rho5)
    left_native_then_gauge = weight * gauge20 @ reflection20 @ start20 @ reflection20.conj().T @ gauge20.conj().T
    right_native_then_gauge = erasure(weight * logical @ rho5 @ logical.conj().T)
    left_gauge_then_native = weight * reflection20 @ gauge20 @ start20 @ gauge20.conj().T @ reflection20.conj().T
    right_gauge_then_native = erasure(weight * logical @ rho5 @ logical.conj().T)
    mixed_sequence_residual = max(
        float(np.linalg.norm(left_native_then_gauge - right_native_then_gauge)),
        float(np.linalg.norm(left_gauge_then_native - right_gauge_then_native)),
    )
    check(
        "branchwise erasure intertwining survives a mixed native-reflection/gauge sequence",
        mixed_sequence_residual < 1e-12,
        "supplemental_actual_matrix_residual",
        mixed_sequence_residual,
        0.0,
    )

    # Native G31 covariance: every actual source reflection normalizes the projective Pauli labels.
    normalizer_residuals = []
    native_permutations = []
    for native in native_reflections:
        row = []
        for pauli in paulis_np:
            index, _phase, residual = pauli_match(native @ pauli @ native.conj().T, paulis_np)
            row.append(index)
            normalizer_residuals.append(residual)
        native_permutations.append(row)
    all_permutations = all(sorted(row) == list(range(16)) for row in native_permutations)
    orbit = {1}
    changed = True
    while changed:
        previous = len(orbit)
        orbit |= {row[index] for row in native_permutations for index in orbit}
        changed = len(orbit) != previous
    check(
        "all sixty actual native reflections normalize and permute the projective Pauli records",
        max(normalizer_residuals) < 1e-12 and all_permutations and orbit == set(range(1, 16)),
        "supplemental_actual_matrix_residual",
        {"reflection_count": len(native_reflections), "max_residual": max(normalizer_residuals), "all_rows_permutations": all_permutations, "nonidentity_orbit_size": len(orbit)},
        {"reflection_count": 60, "max_residual": 0.0, "all_rows_permutations": True, "nonidentity_orbit_size": 15},
    )
    check(
        "all native reflections preserve the actual quartic five-space",
        max(native_reduce_residuals) < 1e-12,
        "supplemental_actual_matrix_residual",
        max(native_reduce_residuals),
        0.0,
    )

    # Decisive classification inside the original primitive 60-reflection menu.
    # H4 is expressed in the exact Pauli operator basis, where every native
    # adjoint is a signed permutation.  H5 is reconstructed directly by acting
    # with the same Gaussian-integer reflections on the unnormalized five-column
    # quartic word basis, with no event-to-S6 matching step.
    signed_actions4 = []
    exact_signed_action_flags = []
    for row_index, native in enumerate(native_reflections_exact):
        action = sp.zeros(16)
        for column, pauli in enumerate(paulis_sp):
            conjugated = native * pauli * native.H
            coefficients = [sp.simplify(sp.trace(candidate.H * conjugated) / 4) for candidate in paulis_sp]
            nonzero = [(target, coefficient) for target, coefficient in enumerate(coefficients) if coefficient != 0]
            if len(nonzero) != 1 or nonzero[0][1] not in (sp.Integer(-1), sp.Integer(1)):
                raise ArithmeticError(f"native reflection {row_index} is not an exact signed Pauli permutation")
            target, sign = nonzero[0]
            action[target, column] = sign
            exact_signed_action_flags.append(
                (conjugated - sign * paulis_sp[target]).applyfunc(sp.simplify) == sp.zeros(4)
            )
        if action.T * action != sp.eye(16):
            raise ArithmeticError(f"native reflection {row_index} signed Pauli action is not orthogonal")
        signed_actions4.append(action)
    check(
        "the exact Gaussian-integer reflections act by signed permutations on the Pauli basis",
        all(exact_signed_action_flags) and len(signed_actions4) == 60,
        "exact_symbolic_source_action",
        {"reflection_count": len(signed_actions4), "verified_basis_images": sum(exact_signed_action_flags)},
        {"reflection_count": 60, "verified_basis_images": 960},
    )
    generator_columns4 = [sp.Matrix(action - sp.eye(16)).reshape(256, 1) for action in signed_actions4]
    primitive_generator4 = sp.Matrix.hstack(*generator_columns4)

    quartic_numerators = []
    exact_quartic_carrier_flags = []
    quartic_involution_flags = []
    generator_columns5 = []
    for twice_reflection in native_twice_reflections_gaussian:
        numerator, exact_carrier = exact_quartic_code_action(twice_reflection, word_sets)
        quartic_numerators.append(numerator)
        exact_quartic_carrier_flags.append(exact_carrier)
        quartic_involution_flags.append(numerator * numerator == 256 * sp.eye(5))
        # A=numerator/16 and A^{-T}=A^T because A^2=I.  Multiplying
        # this whole H5 row block by 256 keeps the rational rank and
        # leaves an integer generator matrix.
        adjoint_scaled = sp.kronecker_product(numerator, numerator.T)
        generator_columns5.append(sp.Matrix(adjoint_scaled - 256 * sp.eye(25)).reshape(625, 1))
    check(
        "all sixty actual quartic code actions are exact integer carrier involutions",
        all(exact_quartic_carrier_flags) and all(quartic_involution_flags),
        "exact_gaussian_integer_quartic_action",
        {
            "reflection_count": len(quartic_numerators),
            "exact_carrier_identities": sum(exact_quartic_carrier_flags),
            "involutions_after_dividing_by_16": sum(quartic_involution_flags),
        },
        {
            "reflection_count": 60,
            "exact_carrier_identities": 60,
            "involutions_after_dividing_by_16": 60,
        },
    )
    primitive_generator5 = sp.Matrix.hstack(*generator_columns5)
    primitive_generator45 = primitive_generator4.col_join(primitive_generator5)
    primitive_rank4 = int(primitive_generator4.rank())
    primitive_rank5 = int(primitive_generator5.rank())
    primitive_rank45 = int(primitive_generator45.rank())
    certificate_prime = 1009
    certificate_rows = independent_rows_mod_prime(primitive_generator45, certificate_prime)
    certificate_minor = primitive_generator45.extract(certificate_rows, list(range(60)))
    certificate_determinant = determinant_mod_prime(certificate_minor, certificate_prime)
    combined_matrix_hash = integer_matrix_sha256(primitive_generator45)
    check(
        "the complete H4 and H5 generators uniquely determine all sixty primitive reflection rates",
        primitive_rank4 == 55
        and primitive_rank5 == 15
        and primitive_rank45 == 60
        and len(certificate_rows) == 60
        and certificate_determinant != 0,
        "exact_integer_linear_rank_with_modular_minor_certificate",
        {
            "primitive_columns": 60,
            "H4_rank": primitive_rank4,
            "H4_nullity": 60 - primitive_rank4,
            "H5_rank": primitive_rank5,
            "combined_H4_H5_rank": primitive_rank45,
            "combined_nullity": 60 - primitive_rank45,
            "modular_certificate": {
                "prime": certificate_prime,
                "matrix_shape": [primitive_generator45.rows, primitive_generator45.cols],
                "canonical_integer_matrix_sha256": combined_matrix_hash,
                "minor_row_indices_zero_based": certificate_rows,
                "minor_determinant_mod_prime": certificate_determinant,
            },
            "H5_source_certificate": {
                "source": "direct action of the same Gaussian-integer reflections on the unnormalized five-column quartic word basis B",
                "identity": "(2r_j)^tensor4 B = B A_num,j exactly; rho5(r_j)=A_num,j/16",
                "involution": "A_num,j^2=256 I5",
                "rank_scaling": "the H5 adjoint-generator row block is multiplied uniformly by 256 to remain integer",
            },
        },
        {
            "primitive_columns": 60,
            "H4_rank": 55,
            "H4_nullity": 5,
            "H5_rank": 15,
            "combined_H4_H5_rank": 60,
            "combined_nullity": 0,
            "minor_determinant_mod_prime": "nonzero",
        },
    )

    # Exact gauge correlation and continuous equal-rate generator.
    observable = paulis_sp[1]
    depolarized_observable = sum(
        (weight * pauli * observable * pauli.H for weight, pauli in zip(weights, paulis_sp)),
        sp.zeros(4),
    )
    correlation = sp.simplify(sp.trace(observable * depolarized_observable) / 4)
    check(
        "normalized H20 traceless gauge two-time correlation equals q",
        correlation == q_exact,
        "exact_symbolic",
        str(correlation),
        str(q_exact),
    )
    generator_flags = []
    generator_eigenvalues = []
    for observable in paulis_sp:
        generated = sum(
            (pauli * observable * pauli.H - observable for pauli in paulis_sp[1:]),
            sp.zeros(4),
        ) / 16
        expected = sp.zeros(4) if observable == i4 else -observable
        generator_flags.append(generated == expected)
        generator_eigenvalues.append("0" if observable == i4 else "-gamma")
    check(
        "fifteen equal Pauli jump rates gamma/16 give gamma(D0-id)",
        all(generator_flags),
        "exact_symbolic",
        {"verified_basis_elements": sum(generator_flags), "eigenvalues": {"identity": "0", "traceless_paulis": "-gamma"}},
        {"verified_basis_elements": 16, "eigenvalues": {"identity": "0", "traceless_paulis": "-gamma"}},
    )

    # Negative control: a central scalar jump has no internal adjoint action.
    scalar_jump = i29 - 2 * e20
    representative_blocks = [
        block_diag(paulis_sp[1], sp.zeros(5), sp.zeros(20)),
        block_diag(sp.zeros(4), matrix_unit(5, 0, 1), sp.zeros(20)),
        block_diag(sp.zeros(4), sp.zeros(5), sp.kronecker_product(matrix_unit(5, 0, 1), paulis_sp[2])),
    ]
    scalar_jump_flags = [scalar_jump * observable * scalar_jump.H == observable for observable in representative_blocks]
    check(
        "negative control: the naive central e20 sign jump has no internal block action",
        scalar_jump.H * scalar_jump == i29 and all(scalar_jump_flags),
        "exact_symbolic_negative_control",
        {"unitary": scalar_jump.H * scalar_jump == i29, "fixed_representatives": sum(scalar_jump_flags)},
        {"unitary": True, "fixed_representatives": 3},
    )

    sources = [
        REPO / "tfpt_explorer/process.py",
        REPO / "tfpt_explorer/quartic_selection.py",
        REPO / "experiments/theory-contracts/compiler-quartic-carry-20260919/checker.py",
        REPO / "experiments/theory-contracts/compiler-quartic-carry-20260919/PROOF.txt",
    ]
    passed = sum(check_row["ok"] for check_row in checks)
    return {
        "status": "PASS" if passed == len(checks) else "FAIL",
        "verdict": "EXACT_CONDITIONAL_HIDDEN_GAUGE_FAMILY_WITH_BRANCHWISE_COMMON_RECORDS",
        "checks_passed": passed,
        "checks_total": len(checks),
        "dimensions": {
            "H4": 4,
            "H5": 5,
            "H20": 20,
            "full_tensor_pauli_fixed": full_tensor_fixed_dimension,
            "symmetric_quartic_pauli_fixed": degree4_dimension,
        },
        "family": {
            "central_projector": "e20=(1/4) sum_{k=0}^3 i^k J^k",
            "source_unitaries": "u_m=(1-e20)+e20 P_m, m=0,...,15",
            "weights": "w0=(1+15q)/16, w_m=(1-q)/16 for m>0, -1/15<=q<=1",
            "represented_channels": {"H4": "id", "H5": "id", "H20": "id5 tensor D_q"},
            "continuous_generator": "L=(gamma/16) sum_{m=1}^{15}(Ad(u_m)-id), q(t)=exp(-gamma t)",
            "fixed_algebra_H20": "B(H5) tensor I4; I20/20 is stationary but not unique",
        },
        "classification_theorems": {
            "complete_logical_nondisturbance": {
                "premise": "Phi20 is CPTP and Tr_4 Phi20(X)=Tr_4 X for every joint operator X",
                "conclusion": "Phi20=id5 tensor D for a CPTP gauge channel D",
                "proof": "Phi20* fixes M5 tensor I pointwise. Its unitaries saturate both Schwarz equalities and lie in the multiplicative domain, making Phi20* an M5-bimodule map. The commutant then forces Phi20*(I tensor B)=I tensor D*(B).",
                "scope": "This premise is stronger than preservation only of the erasure range rho tensor I4/4.",
            },
            "native_covariance": {
                "premise": "The residual channel above is covariant under the diagonal native G31 product representation",
                "conclusion": "D=D_q with -1/15<=q<=1",
                "proof": "Pauli conjugations give fifteen distinct nontrivial characters and the native Clifford quotient acts transitively on the fifteen traceless Pauli labels, so the covariant superoperator has one scalar q on that sector. Complete positivity gives the stated interval.",
            },
            "primitive_menu_rigidity": {
                "premise": "Both complete H4 and H5 Poisson generators are fixed and every primitive jump is one of the actual sixty native reflections",
                "conclusion": "All sixty rate differences vanish; no hidden H20 rate deformation exists inside this restricted menu.",
                "certificate": {
                    "H4_rank": primitive_rank4,
                    "H4_H5_rank": primitive_rank45,
                    "columns": 60,
                    "prime": certificate_prime,
                    "matrix_shape": [primitive_generator45.rows, primitive_generator45.cols],
                    "canonical_integer_matrix_sha256": combined_matrix_hash,
                    "minor_row_indices_zero_based": certificate_rows,
                    "minor_determinant_mod_prime": certificate_determinant,
                },
                "source_scope": "Both H4 signed actions and H5 quartic actions are reconstructed exactly from the same Gaussian-integer rays; no numerical event-to-transposition matching enters the rank certificate.",
            },
        },
        "scope": {
            "proved": [
                "exact source-algebra normalization and unitary central gating",
                "common state-independent marked record probabilities on H4, H5, and H20",
                "branchwise and sequential erasure intertwining",
                "native-G31-covariant gauge channel and its one-parameter continuous semigroup",
                "exact integer rigidity of the original sixty-reflection primitive-rate menu once the complete H4 and H5 generators are fixed",
            ],
            "extra_menu_choice": "The u_m are unitary elements of the represented source C*-algebra, not generally single original G31 words or primitive reflections. Declaring them elementary Kraus events is an additional source-event menu choice.",
            "factorizing_scope": "Factorization follows, rather than being assumed, if complete logical nondisturbance Tr_4 Phi(X)=Tr_4 X is required for every joint X. Preservation only of erasure states rho tensor I4/4 does not imply it.",
            "record_algebra_scope": "The sixteen projective Pauli record labels transform equivariantly on the represented block algebra A_joint=M4 directsum M5 directsum M20, where sector-dependent central phases cancel under adjoint action. This is not literal branch covariance on full End(H29), whose inter-sector coherences see those relative phases. Refining to sixty-four phased Pauli records would define a different record alphabet and cross-block channel.",
            "clock_scope": "Covariance is only with individual elements of the diagonal native G31 product representation. The reviewed common-source files define the H20 lifts C_j tensor conjugate(u_j) of native events, but no H20 lift of the separate marked logical R was found. No commutation with a coherent source-word clock, P1/source-sigma clock bridge, or physical time unit is asserted.",
            "carrier_scope": "H20 is the decoded erasure carrier H5 tensor conjugate(H4); identification with a charged E8 current block requires a separate bridge.",
            "symmetry_scope": "Pauli invariance alone gives a 16-dimensional fixed space in the full fourth tensor power. The five-dimensional code additionally uses Sym^4(V4).",
        },
        "checks": checks,
        "provenance": [{"path": str(path), "sha256": sha256(path)} for path in sources],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    arguments = parser.parse_args()
    report = build_report()
    arguments.out.parent.mkdir(parents=True, exist_ok=True)
    arguments.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: report[key] for key in ("status", "verdict", "checks_passed", "checks_total")}, indent=2))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
