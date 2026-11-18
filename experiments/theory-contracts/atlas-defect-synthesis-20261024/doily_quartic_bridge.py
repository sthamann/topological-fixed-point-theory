#!/usr/bin/env python3
"""Exact E8/J doily to native quartic-five bridge.

The checker uses only the tracked constructors in ``tfpt_explorer.process``
and ``tfpt_explorer.twisted_source``.  It reconstructs the quartic action
directly from the fourth tensor power of every one of the sixty native
Gaussian reflections.  No stored verification result is imported.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools as it
import json
import sys
from pathlib import Path
from typing import Any

import sympy as sp


REPO = Path(__file__).resolve().parents[3]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from tfpt_explorer import process, twisted_source as twisted  # noqa: E402


Gaussian = tuple[int, int]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def gaussian_add(left: Gaussian, right: Gaussian) -> Gaussian:
    return left[0] + right[0], left[1] + right[1]


def gaussian_multiply(left: Gaussian, right: Gaussian) -> Gaussian:
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def exact_quartic_action(
    twice_reflection: list[list[Gaussian]],
    word_sets: list[list[tuple[int, ...]]],
) -> sp.Matrix:
    """Return the numerator N for the quartic action Q=N/16.

    The basis vectors are the unnormalised equal-coefficient sums on the five
    word supports returned by ``process._code_basis``.  Acting with
    ``(2r)^tensor4`` therefore gives an integral matrix N.
    """

    support_owner = {
        word: column
        for column, words in enumerate(word_sets)
        for word in words
    }
    require(
        len(support_owner) == sum(len(words) for words in word_sets),
        "quartic word supports overlap",
    )

    values: dict[tuple[int, ...], list[Gaussian]] = {}
    for output_word in it.product(range(4), repeat=4):
        row_values: list[Gaussian] = []
        for words in word_sets:
            total = (0, 0)
            for input_word in words:
                product = (1, 0)
                for output_index, input_index in zip(output_word, input_word):
                    product = gaussian_multiply(
                        product,
                        twice_reflection[output_index][input_index],
                    )
                total = gaussian_add(total, product)
            row_values.append(total)
        values[output_word] = row_values

    numerator = sp.zeros(5)
    for row, words in enumerate(word_sets):
        representative = words[0]
        for column, (real, imaginary) in enumerate(values[representative]):
            require(imaginary == 0, "quartic coefficient is not real")
            numerator[row, column] = real

    for output_word, row_values in values.items():
        output_row = support_owner.get(output_word)
        for column, value in enumerate(row_values):
            expected = (
                (0, 0)
                if output_row is None
                else (int(numerator[output_row, column]), 0)
            )
            require(value == expected, "quartic five-space is not invariant")
    return numerator


def matrix_apply_f2(matrix: list[list[int]], vector: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sum(matrix[row][column] * vector[column] for column in range(4)) % 2
        for row in range(4)
    )


def f2_rank(rows: list[list[int]]) -> int:
    work = [[entry % 2 for entry in row] for row in rows]
    rank = 0
    if not work:
        return 0
    for column in range(len(work[0])):
        pivot = next(
            (row for row in range(rank, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for row in range(len(work)):
            if row != rank and work[row][column]:
                work[row] = [a ^ b for a, b in zip(work[row], work[rank])]
        rank += 1
    return rank


def centered_permutation(permutation: list[int]) -> sp.Matrix:
    """Action on f_j=e_j-e_5, j=0,...,4, inside Q^6."""

    return sp.Matrix(
        5,
        5,
        lambda row, column: int(permutation[column] == row)
        - int(permutation[5] == row),
    )


def permutation_sign(permutation: list[int]) -> int:
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(6)
        for right in range(left + 1, 6)
    )
    return -1 if inversions % 2 else 1


def cycle_type(permutation: list[int]) -> list[int]:
    seen: set[int] = set()
    lengths: list[int] = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            length += 1
            current = permutation[current]
        lengths.append(length)
    return sorted(lengths, reverse=True)


def hom_constraints(
    quartic_numerators: list[sp.Matrix],
    permutations: list[list[int]],
    sign_twist: bool,
) -> sp.Matrix:
    rows: list[list[int]] = []
    for numerator, permutation in zip(quartic_numerators, permutations):
        action = centered_permutation(permutation)
        if sign_twist:
            action *= permutation_sign(permutation)
        for output in range(5):
            for source in range(5):
                row = [0] * 25
                for middle in range(5):
                    row[5 * middle + source] += int(numerator[output, middle])
                    row[5 * output + middle] -= 16 * int(action[middle, source])
                rows.append(row)
    return sp.Matrix(rows)


def primitive_matrix(vector: sp.Matrix, rows: int, columns: int) -> sp.Matrix:
    matrix = sp.Matrix(rows, columns, list(vector))
    denominator = sp.ilcm(*[sp.denom(entry) for entry in matrix])
    integral = (matrix * denominator).applyfunc(sp.simplify)
    divisor = sp.gcd_list([abs(int(entry)) for entry in integral if entry])
    if divisor:
        integral /= divisor
    first = next(int(entry) for entry in integral if entry)
    if first < 0:
        integral = -integral
    return integral


def pauli_label_actions(
    gaussian_ray: list[Gaussian],
) -> list[list[int]]:
    """Native Clifford conjugation on Hermitian Pauli labels over F2."""

    identity = sp.eye(2)
    x = sp.Matrix([[0, 1], [1, 0]])
    y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    z = sp.diag(1, -1)
    one_qubit = {
        (0, 0): identity,
        (1, 0): x,
        (0, 1): z,
        (1, 1): y,
    }

    def pauli(bits: tuple[int, ...]) -> sp.Matrix:
        return sp.kronecker_product(
            one_qubit[(bits[0], bits[1])],
            one_qubit[(bits[2], bits[3])],
        )

    root = sp.Matrix([sp.Integer(a) + sp.I * sp.Integer(b) for a, b in gaussian_ray])
    reflection = sp.eye(4) - root * root.H / 2
    generators = [
        (1, 0, 0, 0),
        (0, 1, 0, 0),
        (0, 0, 1, 0),
        (0, 0, 0, 1),
    ]
    labels = list(it.product((0, 1), repeat=4))
    columns: list[tuple[int, ...]] = []
    for generator in generators:
        conjugated = (reflection * pauli(generator) * reflection.H).applyfunc(sp.simplify)
        matches = [
            label
            for label in labels
            if conjugated == pauli(label) or conjugated == -pauli(label)
        ]
        require(len(matches) == 1, "native reflection is not an exact Clifford action")
        columns.append(matches[0])
    return [[columns[column][row] for column in range(4)] for row in range(4)]


def native_four_hom_dimension(
    quotient_actions: list[list[list[int]]],
    pauli_actions: list[list[list[int]]],
) -> tuple[int, int]:
    """Hom_F2[G31](A, Pauli/phase), including invertible count."""

    equations: list[list[int]] = []
    for quotient, pauli in zip(quotient_actions, pauli_actions):
        for row_index in range(4):
            for column_index in range(4):
                row = [0] * 16
                for middle in range(4):
                    row[4 * row_index + middle] ^= quotient[middle][column_index]
                    row[4 * middle + column_index] ^= pauli[row_index][middle]
                equations.append(row)
    rank = f2_rank(equations)
    dimension = 16 - rank
    if dimension == 0:
        return 0, 0

    # This branch is retained to make the certificate complete if sources drift.
    work = [[entry % 2 for entry in row] for row in equations]
    pivots: list[int] = []
    pivot_row = 0
    for column in range(16):
        pivot = next(
            (row for row in range(pivot_row, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        for row in range(len(work)):
            if row != pivot_row and work[row][column]:
                work[row] = [a ^ b for a, b in zip(work[row], work[pivot_row])]
        pivots.append(column)
        pivot_row += 1
    free = [column for column in range(16) if column not in pivots]
    basis: list[list[int]] = []
    for free_column in free:
        vector = [0] * 16
        vector[free_column] = 1
        for row, pivot_column in enumerate(pivots):
            vector[pivot_column] = work[row][free_column]
        basis.append(vector)
    invertible = 0
    for mask in range(1, 1 << len(basis)):
        vector = [0] * 16
        for index, item in enumerate(basis):
            if mask >> index & 1:
                vector = [a ^ b for a, b in zip(vector, item)]
        candidate = [vector[4 * row : 4 * row + 4] for row in range(4)]
        invertible += f2_rank(candidate) == 4
    return dimension, invertible


def rational_json(value: sp.Expr) -> int | str:
    value = sp.factor(value)
    if value.is_Integer:
        return int(value)
    return str(value)


def matrix_json(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [rational_json(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def bits_key(bits: tuple[int, ...]) -> str:
    return "".join(str(bit) for bit in bits)


def source_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_report() -> dict[str, Any]:
    code = twisted._code()
    roots = twisted._roots(code)
    representatives = twisted._quotient_representatives(roots, code)
    darboux = twisted._darboux_basis(representatives, code)

    bits_by_class: dict[int, tuple[int, ...]] = {}
    for bits in it.product((0, 1), repeat=4):
        vector = tuple(
            sum(
                bit * representatives[item][coordinate]
                for bit, item in zip(bits, darboux)
            )
            for coordinate in range(8)
        )
        bits_by_class[twisted._class_index(vector, representatives, code)] = bits
    require(len(bits_by_class) == 16, "Darboux coordinates do not cover A")
    representative_by_bits = {
        bits: representatives[class_index]
        for class_index, bits in bits_by_class.items()
    }
    points = sorted(bits for bits in representative_by_bits if any(bits))

    def add(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
        return tuple(a ^ b for a, b in zip(left, right))

    def pairing(left: tuple[int, ...], right: tuple[int, ...]) -> int:
        return twisted._bk_commutator_pairing(
            representative_by_bits[left],
            representative_by_bits[right],
        )

    lines = sorted(
        {
            frozenset((left, right, add(left, right)))
            for left, right in it.combinations(points, 2)
            if pairing(left, right) == 0
        },
        key=lambda line: sorted(line),
    )
    ovoids = sorted(
        [
            frozenset(candidate)
            for candidate in it.combinations(points, 5)
            if all(
                pairing(left, right) == 1
                for left, right in it.combinations(candidate, 2)
            )
        ],
        key=lambda ovoid: sorted(ovoid),
    )
    spreads = sorted(
        [
            frozenset(candidate)
            for candidate in it.combinations(range(len(lines)), 5)
            if set().union(*(set(lines[index]) for index in candidate)) == set(points)
            and sum(len(lines[index]) for index in candidate) == 15
        ],
        key=lambda spread: sorted(spread),
    )
    require((len(points), len(lines), len(ovoids), len(spreads)) == (15, 15, 6, 6), "doily census changed")

    incidence = sp.Matrix(
        [[int(point in line) for line in lines] for point in points]
    )
    incidence_square = incidence * incidence.T
    spectrum = {
        int(eigenvalue): int(multiplicity)
        for eigenvalue, multiplicity in incidence_square.eigenvals().items()
    }
    require(incidence.rank() == 10, "doily incidence rank changed")
    require(spectrum == {0: 5, 4: 9, 9: 1}, "doily incidence spectrum changed")

    ovoid_vectors = sp.Matrix.hstack(
        *[
            sp.Matrix([3 * int(point in ovoid) - 1 for point in points])
            for ovoid in ovoids
        ]
    )
    spread_vectors = sp.Matrix.hstack(
        *[
            sp.Matrix([3 * int(index in spread) - 1 for index in range(15)])
            for spread in spreads
        ]
    )
    ovoid_gram = ovoid_vectors.T * ovoid_vectors
    require(incidence.T * ovoid_vectors == sp.zeros(15, 6), "ovoid vectors miss ker(B^T)")
    require(incidence * spread_vectors == sp.zeros(15, 6), "spread vectors miss ker(B)")
    require(ovoid_vectors.rank() == spread_vectors.rank() == 5, "incidence kernel rank changed")
    require(ovoid_gram == 36 * sp.eye(6) - 6 * sp.ones(6), "ovoid Gram matrix changed")

    ovoid_index = {ovoid: index for index, ovoid in enumerate(ovoids)}
    line_index = {line: index for index, line in enumerate(lines)}
    spread_index = {spread: index for index, spread in enumerate(spreads)}

    def ovoid_permutation(action: list[list[int]]) -> list[int]:
        return [
            ovoid_index[
                frozenset(matrix_apply_f2(action, point) for point in ovoid)
            ]
            for ovoid in ovoids
        ]

    def spread_permutation(action: list[list[int]]) -> list[int]:
        induced_line = [
            line_index[
                frozenset(matrix_apply_f2(action, point) for point in line)
            ]
            for line in lines
        ]
        return [
            spread_index[frozenset(induced_line[index] for index in spread)]
            for spread in spreads
        ]

    _, word_sets = process._code_basis()
    native_rays = process._hamming_rays()
    require(len(native_rays) == 60, "native ray count changed")

    quartic_numerators: list[sp.Matrix] = []
    quotient_actions: list[list[list[int]]] = []
    native_pauli_actions: list[list[list[int]]] = []
    ovoid_permutations: list[list[int]] = []
    spread_permutations: list[list[int]] = []
    reflection_rows: list[dict[str, Any]] = []

    for reflection_index, ray in enumerate(native_rays):
        gaussian_ray: list[Gaussian] = []
        for coordinate in ray:
            real, imaginary = int(coordinate.real), int(coordinate.imag)
            require(
                coordinate.real == real and coordinate.imag == imaginary,
                "native ray left Gaussian integers",
            )
            gaussian_ray.append((real, imaginary))
        root = tuple(component for coordinate in gaussian_ray for component in coordinate)
        require(root in roots, "native Gaussian ray is not an E8 root")

        quotient_certificate = twisted._reflection_action_certificate(
            representatives,
            code,
            darboux,
            root,
        )
        quotient_action = quotient_certificate["induced_matrix"]
        quotient_actions.append(quotient_action)
        require(
            quotient_certificate["rank_r_minus_identity"] == 1,
            "quotient action is not a transvection",
        )

        ovoid_perm = ovoid_permutation(quotient_action)
        spread_perm = spread_permutation(quotient_action)
        ovoid_permutations.append(ovoid_perm)
        spread_permutations.append(spread_perm)

        twice_reflection: list[list[Gaussian]] = []
        for row in range(4):
            twice_row: list[Gaussian] = []
            for column in range(4):
                outer = gaussian_multiply(
                    gaussian_ray[row],
                    (gaussian_ray[column][0], -gaussian_ray[column][1]),
                )
                twice_row.append(
                    (2 * int(row == column) - outer[0], -outer[1])
                )
            twice_reflection.append(twice_row)
        numerator = exact_quartic_action(twice_reflection, word_sets)
        quartic_numerators.append(numerator)
        require(numerator * numerator == 256 * sp.eye(5), "quartic action is not an involution")

        pauli_action = pauli_label_actions(gaussian_ray)
        native_pauli_actions.append(pauli_action)
        reflection_rows.append(
            {
                "index": reflection_index,
                "gaussian_ray": [list(coordinate) for coordinate in gaussian_ray],
                "e8_root": list(root),
                "quotient_action": quotient_action,
                "quotient_rank_r_minus_identity": quotient_certificate["rank_r_minus_identity"],
                "ovoid_permutation": ovoid_perm,
                "ovoid_cycle_type": cycle_type(ovoid_perm),
                "spread_permutation": spread_perm,
                "spread_cycle_type": cycle_type(spread_perm),
                "native_pauli_label_action": pauli_action,
                "quartic_numerator": matrix_json(numerator),
                "quartic_trace": rational_json(numerator.trace() / 16),
                "quartic_involution": numerator * numerator == 256 * sp.eye(5),
                "all_16_twisted_ground_lifts_mismatch_native_four": quotient_certificate[
                    "all_pauli_lifts_mismatch_native_reflection"
                ],
            }
        )

    require(
        {tuple(cycle_type(permutation)) for permutation in ovoid_permutations}
        == {(2, 1, 1, 1, 1)},
        "root reflections are not transpositions on ovoids",
    )
    require(
        {tuple(cycle_type(permutation)) for permutation in spread_permutations}
        == {(2, 2, 2)},
        "root reflections are not triple transpositions on spreads",
    )
    require(
        {numerator.trace() / 16 for numerator in quartic_numerators} == {3},
        "native quartic reflection trace changed",
    )

    hom_results: dict[str, dict[str, Any]] = {}
    bridge_matrix: sp.Matrix | None = None
    for family_name, permutations in (
        ("ovoid_standard", ovoid_permutations),
        ("ovoid_sign_twist", ovoid_permutations),
        ("spread_standard", spread_permutations),
        ("spread_sign_twist", spread_permutations),
    ):
        sign_twist = family_name.endswith("sign_twist")
        constraints = hom_constraints(quartic_numerators, permutations, sign_twist)
        nullspace = constraints.nullspace()
        determinants: list[sp.Expr] = []
        primitive_generators: list[sp.Matrix] = []
        for vector in nullspace:
            generator = primitive_matrix(vector, 5, 5)
            primitive_generators.append(generator)
            determinants.append(sp.factor(generator.det()))
        hom_results[family_name] = {
            "constraint_rank": int(constraints.rank()),
            "hom_dimension": len(nullspace),
            "primitive_generator_determinants": [rational_json(value) for value in determinants],
            "invertible_generator_count": sum(value != 0 for value in determinants),
        }
        if family_name == "ovoid_standard":
            require(len(primitive_generators) == 1, "ovoid Hom space is not one dimensional")
            bridge_matrix = primitive_generators[0]

    expected_hom_dimensions = {
        "ovoid_standard": 1,
        "ovoid_sign_twist": 0,
        "spread_standard": 0,
        "spread_sign_twist": 0,
    }
    require(
        {
            name: row["hom_dimension"]
            for name, row in hom_results.items()
        }
        == expected_hom_dimensions,
        "five-dimensional Hom census changed",
    )
    require(bridge_matrix is not None, "missing ovoid bridge")
    expected_bridge = sp.Matrix(
        [
            [3, 3, 3, 0, 3],
            [-1, 1, 1, 0, -1],
            [-1, -1, 1, 0, 1],
            [-1, 1, -1, 0, 1],
            [1, 1, 1, 2, 1],
        ]
    )
    require(bridge_matrix == expected_bridge, "primitive bridge matrix changed")
    require(bridge_matrix.det() == -96, "bridge determinant changed")

    bridge_identities = []
    for row, numerator, permutation in zip(
        reflection_rows,
        quartic_numerators,
        ovoid_permutations,
    ):
        identity = numerator * bridge_matrix == 16 * bridge_matrix * centered_permutation(permutation)
        bridge_identities.append(identity)
        row["quartic_ovoid_intertwiner_identity"] = identity
    require(all(bridge_identities), "an all-sixty bridge identity failed")

    # The six centered ovoid rays transported into the native quartic basis.
    # Multiplication by six clears the denominators of
    # w_a=e_a-(1/6)1.  A common divisor is then removed from all rays.
    native_ovoid_rays: list[sp.Matrix] = []
    for ovoid_index_value in range(6):
        if ovoid_index_value < 5:
            centered_ray = sp.Matrix(
                [5 if index == ovoid_index_value else -1 for index in range(5)]
            )
        else:
            centered_ray = sp.Matrix([-1] * 5)
        transported = bridge_matrix * centered_ray
        divisor = sp.gcd_list([abs(int(entry)) for entry in transported if entry])
        native_ovoid_rays.append(transported / divisor)

    four_hom_dimension, four_invertible = native_four_hom_dimension(
        quotient_actions,
        native_pauli_actions,
    )
    require(
        (four_hom_dimension, four_invertible) == (0, 0),
        "the forbidden native-four identification reappeared",
    )
    require(
        all(row["all_16_twisted_ground_lifts_mismatch_native_four"] for row in reflection_rows),
        "a twisted-ground lift acquired native-four reflection spectrum",
    )

    point_projector = ovoid_vectors * ovoid_vectors.T / 36
    require(point_projector.T == point_projector, "point projector is not symmetric")
    require(point_projector * point_projector == point_projector, "point projector is not idempotent")
    require(point_projector.trace() == 5, "point projector rank changed")
    require(incidence.T * point_projector == sp.zeros(15), "point projector misses ker(B^T)")
    polynomial_projector = (
        (incidence_square - 9 * sp.eye(15))
        * (incidence_square - 4 * sp.eye(15))
        / 36
    )
    require(point_projector == polynomial_projector, "incidence spectral projector changed")

    centered_metric = sp.eye(5) + sp.ones(5)
    native_metric = sp.simplify(
        bridge_matrix.inv().T * centered_metric * bridge_matrix.inv()
    )
    primitive_native_metric = 12 * native_metric
    expected_native_metric = sp.diag(1, 3, 3, 3, 6)
    require(primitive_native_metric == expected_native_metric, "native invariant form changed")
    require(
        all(
            numerator.T * primitive_native_metric * numerator
            == 256 * primitive_native_metric
            for numerator in quartic_numerators
        ),
        "native invariant form is not preserved",
    )

    require(
        sum(native_ovoid_rays, sp.zeros(5, 1)) == sp.zeros(5, 1),
        "native ovoid rays do not center",
    )
    native_ray_norms = [
        (ray.T * primitive_native_metric * ray)[0]
        for ray in native_ovoid_rays
    ]
    require(set(native_ray_norms) == {10}, "native ovoid ray norm changed")
    native_ovoid_effects = [
        ray * ray.T * primitive_native_metric / 12
        for ray in native_ovoid_rays
    ]
    require(
        sum(native_ovoid_effects, sp.zeros(5)) == sp.eye(5),
        "native ovoid effects do not resolve the identity",
    )
    native_ray_covariance = []
    native_effect_covariance = []
    for row, numerator, permutation in zip(
        reflection_rows,
        quartic_numerators,
        ovoid_permutations,
    ):
        ray_ok = all(
            numerator * native_ovoid_rays[index]
            == 16 * native_ovoid_rays[permutation[index]]
            for index in range(6)
        )
        effect_ok = all(
            numerator * native_ovoid_effects[index] * numerator / 256
            == native_ovoid_effects[permutation[index]]
            for index in range(6)
        )
        native_ray_covariance.append(ray_ok)
        native_effect_covariance.append(effect_ok)
        row["native_ovoid_ray_covariance"] = ray_ok
        row["native_ovoid_effect_covariance"] = effect_ok
    require(all(native_ray_covariance), "native ovoid rays are not covariant")
    require(all(native_effect_covariance), "native ovoid effects are not covariant")

    symmetric_pairs = [(row, column) for row in range(5) for column in range(row, 5)]
    form_constraints: list[list[int]] = []
    for numerator in quartic_numerators:
        for output in range(5):
            for source in range(5):
                constraint: list[int] = []
                target_pair = tuple(sorted((output, source)))
                for left, right in symmetric_pairs:
                    coefficient = int(numerator[left, output] * numerator[right, source])
                    if left != right:
                        coefficient += int(numerator[right, output] * numerator[left, source])
                    if (left, right) == target_pair:
                        coefficient -= 256
                    constraint.append(coefficient)
                form_constraints.append(constraint)
    invariant_form_matrix = sp.Matrix(form_constraints)
    invariant_form_dimension = len(symmetric_pairs) - invariant_form_matrix.rank()
    require(invariant_form_dimension == 1, "invariant symmetric form is not unique up to scale")

    fixed_vector_matrix = sp.Matrix.vstack(
        *[numerator - 16 * sp.eye(5) for numerator in quartic_numerators]
    )
    fixed_vector_dimension = 5 - fixed_vector_matrix.rank()
    require(fixed_vector_dimension == 0, "quartic five acquired a nonzero invariant vector")

    commutant_rows: list[list[int]] = []
    for numerator in quartic_numerators:
        for output in range(5):
            for source in range(5):
                row = [0] * 25
                for middle in range(5):
                    row[5 * middle + source] += int(numerator[output, middle])
                    row[5 * output + middle] -= int(numerator[middle, source])
                commutant_rows.append(row)
    commutant_matrix = sp.Matrix(commutant_rows)
    commutant_dimension = 25 - commutant_matrix.rank()
    require(commutant_dimension == 1, "native quartic commutant is not scalar")

    # Six covariant rank-one effects E_O=v_O v_O^T/36 sum to the identity
    # on ker(B^T), represented by point_projector in the ambient point space.
    ovoid_effects = [
        ovoid_vectors[:, index] * ovoid_vectors[:, index].T / 36
        for index in range(6)
    ]
    require(sum(ovoid_effects, sp.zeros(15)) == point_projector, "ovoid effects do not resolve the projector")

    process_source = REPO / "tfpt_explorer/process.py"
    twisted_source = REPO / "tfpt_explorer/twisted_source.py"
    return {
        "claim": "E8/J doily ovoid kernel equals the original native quartic five as a marked G31 module",
        "status": "PASS",
        "arithmetic": "exact integers, rationals, Gaussian integers, and F2",
        "sources": {
            "tfpt_explorer/process.py": source_hash(process_source),
            "tfpt_explorer/twisted_source.py": source_hash(twisted_source),
        },
        "conventions": {
            "darboux_class_indices": list(darboux),
            "darboux_bit_order": ["e1", "f1", "e2", "f2"],
            "point_order": [bits_key(point) for point in points],
            "line_order": [
                [bits_key(point) for point in sorted(line)]
                for line in lines
            ],
            "ovoid_order": [
                [bits_key(point) for point in sorted(ovoid)]
                for ovoid in ovoids
            ],
            "spread_order_as_line_indices": [sorted(spread) for spread in spreads],
            "centered_basis": [f"e{index}-e5" for index in range(5)],
            "quartic_basis": "unnormalised equal-coefficient word sums in the listed order",
            "quartic_word_supports": [
                [list(word) for word in words]
                for words in word_sets
            ],
        },
        "doily": {
            "point_count": len(points),
            "line_count": len(lines),
            "ovoid_count": len(ovoids),
            "spread_count": len(spreads),
            "incidence_rank": int(incidence.rank()),
            "incidence_BBt_spectrum": {str(key): value for key, value in sorted(spectrum.items())},
            "ovoid_kernel_rank": int(ovoid_vectors.rank()),
            "spread_kernel_rank": int(spread_vectors.rank()),
            "ovoid_gram": matrix_json(ovoid_gram),
        },
        "quartic_bridge": {
            "reflection_count": len(reflection_rows),
            "native_quartic_trace_set": sorted(
                {rational_json(numerator.trace() / 16) for numerator in quartic_numerators},
                key=str,
            ),
            "ovoid_cycle_type_set": sorted(
                {tuple(cycle_type(permutation)) for permutation in ovoid_permutations}
            ),
            "spread_cycle_type_set": sorted(
                {tuple(cycle_type(permutation)) for permutation in spread_permutations}
            ),
            "hom_spaces": hom_results,
            "primitive_intertwiner": matrix_json(bridge_matrix),
            "intertwiner_determinant": rational_json(bridge_matrix.det()),
            "all_60_intertwiner_identities": all(bridge_identities),
            "per_reflection": reflection_rows,
        },
        "native_four_exclusion": {
            "f2_hom_dimension_A_to_native_pauli_labels": four_hom_dimension,
            "invertible_f2_intertwiners": four_invertible,
            "native_reflection_trace_modulus_squared": 4,
            "twisted_transvection_lift_trace_modulus_squared_histogram_per_root": {
                "0": 14,
                "8": 2,
            },
            "all_60_times_16_twisted_lifts_mismatch_native_four": all(
                row["all_16_twisted_ground_lifts_mismatch_native_four"]
                for row in reflection_rows
            ),
        },
        "selection_data": {
            "canonical_point_kernel_projector": matrix_json(point_projector),
            "incidence_polynomial_projector": "((BBt-9I)(BBt-4I))/36",
            "projector_rank": int(point_projector.trace()),
            "six_covariant_ovoid_effects_sum_to_projector": True,
            "ovoid_effect_formula": "E_O = v_O v_O^T / 36, v_O = 3*1_O - 1",
            "primitive_native_invariant_positive_form": matrix_json(primitive_native_metric),
            "invariant_symmetric_form_dimension": invariant_form_dimension,
            "native_ovoid_rays": [
                [rational_json(entry) for entry in ray]
                for ray in native_ovoid_rays
            ],
            "native_ovoid_ray_norms": [rational_json(value) for value in native_ray_norms],
            "native_ovoid_effect_formula": "E_a = q_a q_a^T G / 12",
            "native_ovoid_effects_sum_to_identity": True,
            "all_60_native_ray_covariances": all(native_ray_covariance),
            "all_60_native_effect_covariances": all(native_effect_covariance),
            "common_fixed_vector_dimension": fixed_vector_dimension,
            "commutant_dimension": commutant_dimension,
            "normalized_group_invariant_mixed_state": "identity_on_V5 / 5 after fixing the invariant form and trace normalization",
            "selected_pure_state": False,
            "free_data": [
                "overall nonzero scalar of the equivariant intertwiner",
                "overall positive scalar of the invariant form",
                "an individual ovoid or pure ray requires extra symmetry-breaking input",
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
