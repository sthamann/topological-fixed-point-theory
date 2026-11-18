#!/usr/bin/env python3
"""Exact minimum-arity test for the native six-ovoid POVM.

The checker reuses the actual October 24 doily/quartic bridge and erasure
isometry.  It does not fit a representation.  Its lower-arity obstruction is
the simultaneous action of the native Pauli kernel and all sixty native root
reflections on Bell projectors.  At arity three it reconstructs the original
twenty-dimensional symmetric-cube carrier and its six effects.
"""

from __future__ import annotations

import argparse
from collections import defaultdict, deque
from fractions import Fraction
import hashlib
import importlib.util
import itertools as it
import json
from pathlib import Path
import sys
from typing import Any

import sympy as sp


GaussianFraction = tuple[Fraction, Fraction]
Matrix4 = tuple[tuple[GaussianFraction, ...], ...]
GF_ZERO: GaussianFraction = (Fraction(0), Fraction(0))
GF_ONE: GaussianFraction = (Fraction(1), Fraction(0))
GF_I: GaussianFraction = (Fraction(0), Fraction(1))


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DEFAULT_OUTPUT = HERE / "native_observation_arity_results.json"
DOILY_SOURCE = (
    REPO
    / "experiments/theory-contracts/atlas-defect-synthesis-20261024/doily_quartic_bridge.py"
)
ERASURE_SOURCE = (
    REPO
    / "experiments/theory-contracts/atlas-defect-synthesis-20261024/native_erasure_readout.py"
)
CORRELATED_SOURCE = (
    REPO
    / "experiments/theory-contracts/compiler-correlated-event-clock-20260919/source_channel.py"
)
PROCESS_SOURCE = REPO / "tfpt_explorer/process.py"
QUARTIC_SOURCE = REPO / "tfpt_explorer/quartic_selection.py"

if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from tfpt_explorer.process import _code_basis  # noqa: E402
from tfpt_explorer.quartic_selection import _pauli_group  # noqa: E402


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def clean(matrix: sp.MatrixBase) -> sp.Matrix:
    return sp.Matrix(matrix).applyfunc(sp.simplify)


def load_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix_sha256(matrix: sp.MatrixBase) -> str:
    payload = json.dumps(
        [[str(sp.simplify(matrix[r, c])) for c in range(matrix.cols)] for r in range(matrix.rows)],
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def tensor_power(matrix: sp.MatrixBase, count: int) -> sp.Matrix:
    out = sp.ones(1, 1)
    for _ in range(count):
        out = sp.kronecker_product(out, matrix)
    return sp.Matrix(out)


def gf_add(left: GaussianFraction, right: GaussianFraction) -> GaussianFraction:
    return left[0] + right[0], left[1] + right[1]


def gf_sum(values: Any) -> GaussianFraction:
    total = GF_ZERO
    for value in values:
        total = gf_add(total, value)
    return total


def gf_negate(value: GaussianFraction) -> GaussianFraction:
    return -value[0], -value[1]


def gf_multiply(left: GaussianFraction, right: GaussianFraction) -> GaussianFraction:
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def gf_conjugate(value: GaussianFraction) -> GaussianFraction:
    return value[0], -value[1]


def gf_scale(scale: Fraction, value: GaussianFraction) -> GaussianFraction:
    return scale * value[0], scale * value[1]


def matrix4_multiply(left: Matrix4, right: Matrix4) -> Matrix4:
    return tuple(
        tuple(
            gf_sum(gf_multiply(left[row][middle], right[middle][column]) for middle in range(4))
            for column in range(4)
        )
        for row in range(4)
    )


def matrix4_adjoint(matrix: Matrix4) -> Matrix4:
    return tuple(
        tuple(gf_conjugate(matrix[column][row]) for column in range(4))
        for row in range(4)
    )


def matrix4_scale(scale: GaussianFraction, matrix: Matrix4) -> Matrix4:
    return tuple(
        tuple(gf_multiply(scale, value) for value in row)
        for row in matrix
    )


def matrix4_trace(matrix: Matrix4) -> GaussianFraction:
    return gf_sum(matrix[index][index] for index in range(4))


def matrix4_projective_key(matrix: Matrix4) -> tuple[int, ...]:
    """Canonical exact key modulo the central fourth roots of unity."""

    keys = []
    for phase in (GF_ONE, GF_I, gf_negate(GF_ONE), gf_negate(GF_I)):
        phased = matrix4_scale(phase, matrix)
        keys.append(
            tuple(
                integer
                for row in phased
                for real, imaginary in row
                for integer in (
                    real.numerator,
                    real.denominator,
                    imaginary.numerator,
                    imaginary.denominator,
                )
            )
        )
    return min(keys)


def matrix4_from_sympy(matrix: sp.MatrixBase) -> Matrix4:
    return tuple(
        tuple(
            (
                Fraction(sp.Rational(sp.re(matrix[row, column]))),
                Fraction(sp.Rational(sp.im(matrix[row, column]))),
            )
            for column in range(4)
        )
        for row in range(4)
    )


def fractional_reflection(ray: list[list[int]]) -> Matrix4:
    vector = [(Fraction(real), Fraction(imaginary)) for real, imaginary in ray]
    return tuple(
        tuple(
            gf_add(
                GF_ONE if row == column else GF_ZERO,
                gf_scale(
                    Fraction(-1, 2),
                    gf_multiply(vector[row], gf_conjugate(vector[column])),
                ),
            )
            for column in range(4)
        )
        for row in range(4)
    )


def native_coset_representatives(
    generators: list[tuple[tuple[int, ...], Matrix4]],
) -> dict[tuple[int, ...], Matrix4]:
    identity_permutation = tuple(range(6))
    identity_matrix: Matrix4 = tuple(
        tuple(GF_ONE if row == column else GF_ZERO for column in range(4))
        for row in range(4)
    )
    representatives = {identity_permutation: identity_matrix}
    queue = deque([identity_permutation])
    while queue:
        permutation = queue.popleft()
        matrix = representatives[permutation]
        for generator_permutation, generator_matrix in generators:
            next_permutation = compose(permutation, generator_permutation)
            if next_permutation in representatives:
                continue
            representatives[next_permutation] = matrix4_multiply(matrix, generator_matrix)
            queue.append(next_permutation)
    return representatives


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(len(left)))


def cycle_type(permutation: tuple[int, ...]) -> tuple[int, ...]:
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
    return tuple(sorted(lengths, reverse=True))


def restriction(permutation: tuple[int, ...], indices: list[int]) -> tuple[int, ...]:
    position = {value: index for index, value in enumerate(indices)}
    return tuple(position[permutation[value]] for value in indices)


def commutant_dimension(matrices: list[sp.Matrix]) -> int:
    """Exact complex dimension of matrices commuting with every input."""

    size = matrices[0].rows
    rows: list[list[sp.Expr]] = []
    for matrix in matrices:
        for output_row in range(size):
            for output_column in range(size):
                row = [sp.S.Zero] * (size * size)
                for middle in range(size):
                    row[size * middle + output_column] += matrix[output_row, middle]
                    row[size * output_row + middle] -= matrix[middle, output_column]
                rows.append(row)
    return size * size - sp.Matrix(rows).rank()


def pauli_congruence_permutation(
    matrix: sp.Matrix,
    paulis: list[sp.Matrix],
) -> tuple[tuple[int, ...], tuple[sp.Expr, ...]]:
    permutation: list[int] = []
    phases: list[sp.Expr] = []
    for pauli in paulis:
        image = clean(matrix * pauli * matrix.T)
        matches = [
            (index, phase)
            for index, target in enumerate(paulis)
            for phase in (sp.S.One, sp.I, -sp.S.One, -sp.I)
            if image == phase * target
        ]
        require(len(matches) == 1, "native congruence did not monomially permute Paulis")
        permutation.append(matches[0][0])
        phases.append(matches[0][1])
    return tuple(permutation), tuple(phases)


def pauli_similarity_permutation(
    matrix: sp.Matrix,
    paulis: list[sp.Matrix],
) -> tuple[tuple[int, ...], tuple[sp.Expr, ...]]:
    """Projective Pauli action for V4 tensor conjugate(V4)."""

    permutation: list[int] = []
    phases: list[sp.Expr] = []
    for pauli in paulis:
        image = clean(matrix * pauli * matrix.H)
        matches = [
            (index, phase)
            for index, target in enumerate(paulis)
            for phase in (sp.S.One, sp.I, -sp.S.One, -sp.I)
            if image == phase * target
        ]
        require(len(matches) == 1, "native similarity did not monomially permute Paulis")
        permutation.append(matches[0][0])
        phases.append(matches[0][1])
    return tuple(permutation), tuple(phases)


def symmetric_cube_projector() -> sp.Matrix:
    projector = sp.zeros(64)
    for word in it.product(range(4), repeat=3):
        column = 16 * word[0] + 4 * word[1] + word[2]
        for permutation in it.permutations(range(3)):
            output = tuple(word[index] for index in permutation)
            row = 16 * output[0] + 4 * output[1] + output[2]
            projector[row, column] += sp.Rational(1, 6)
    return projector


def paired_permutation_group(
    generators: list[tuple[tuple[int, ...], tuple[int, ...]]],
) -> set[tuple[tuple[int, ...], tuple[int, ...]]]:
    identity = (tuple(range(6)), tuple(range(16)))
    group = {identity}
    queue = deque([identity])
    while queue:
        target, bell = queue.popleft()
        for target_generator, bell_generator in generators:
            item = (
                compose(target, target_generator),
                compose(bell, bell_generator),
            )
            if item not in group:
                group.add(item)
                queue.append(item)
    return group


def character_inner_product(
    left: list[int | sp.Rational],
    right: list[int | sp.Rational],
) -> sp.Rational:
    require(len(left) == len(right), "character lengths differ")
    return sp.Rational(sum(a * b for a, b in zip(left, right)), len(left))


def main(output: Path) -> dict[str, Any]:
    doily = load_module("native_observation_arity_doily", DOILY_SOURCE)
    erasure = load_module("native_observation_arity_erasure", ERASURE_SOURCE)
    bridge = doily.build_report()
    rows = bridge["quartic_bridge"]["per_reflection"]
    require(len(rows) == 60, "native reflection count changed")
    require(
        bridge["quartic_bridge"]["hom_spaces"]["ovoid_standard"]["hom_dimension"] == 1,
        "actual quartic-to-ovoid Hom space changed",
    )
    require(
        bridge["quartic_bridge"]["all_60_intertwiner_identities"],
        "actual all-sixty quartic bridge failed",
    )

    # The original constructor orders the phase-free two-qubit Paulis first.
    pauli64 = list(_pauli_group())
    require(len(pauli64) == 64, "native Pauli kernel order changed")
    paulis = [sp.Matrix(matrix) for matrix in pauli64[:16]]
    paulis_fractional = [matrix4_from_sympy(matrix) for matrix in paulis]
    pauli_projective_keys = {matrix4_projective_key(matrix) for matrix in paulis_fractional}
    require(len(pauli_projective_keys) == 16, "phase-free Paulis are not projectively distinct")
    labels = list(it.product((0, 1), repeat=4))
    parity = [(a * b + c * d) % 2 for a, b, c, d in labels]
    symmetric_indices = [index for index, value in enumerate(parity) if value == 0]
    antisymmetric_indices = [index for index, value in enumerate(parity) if value == 1]
    require((len(symmetric_indices), len(antisymmetric_indices)) == (10, 6), "10+6 split changed")
    require(
        all(paulis[index].T == (-1 if parity[index] else 1) * paulis[index] for index in range(16)),
        "Pauli transpose parity changed",
    )

    # Arity one: invariance under N64 already leaves only the identity.
    n1_commutant_dimension = commutant_dimension(paulis)
    require(n1_commutant_dimension == 1, "N64 commutant on V4 is not scalar")

    # Arity two: diagonalize N64 on V4 tensor V4 by the sixteen actual Bell
    # vectors vec(P)/2.  Distinct joint characters make the invariant operator
    # algebra exactly the span of their sixteen rank-one projectors.
    bell_vectors = [pauli.reshape(16, 1) / 2 for pauli in paulis]
    bell_basis = sp.Matrix.hstack(*bell_vectors)
    require(clean(bell_basis.H * bell_basis) == sp.eye(16), "Bell basis is not orthonormal")
    bell_projectors = [clean(vector * vector.H) for vector in bell_vectors]
    require(clean(sum(bell_projectors, sp.zeros(16))) == sp.eye(16), "Bell projectors do not sum to I16")

    joint_characters: list[tuple[sp.Expr, ...]] = []
    diagonal_actions: list[sp.Matrix] = []
    for kernel in paulis:
        action = sp.kronecker_product(kernel, kernel)
        diagonal = clean(bell_basis.H * action * bell_basis)
        require(diagonal == sp.diag(*diagonal.diagonal()), "N64 action is not Bell diagonal")
        require(all(value in (-1, 1) for value in diagonal.diagonal()), "Bell character left {+1,-1}")
        diagonal_actions.append(diagonal)
    for index in range(16):
        joint_characters.append(tuple(diagonal[index, index] for diagonal in diagonal_actions))
    require(len(set(joint_characters)) == 16, "Bell joint characters are not distinct")
    require(
        all(
            clean(sp.kronecker_product(kernel, kernel) * projector
                  - projector * sp.kronecker_product(kernel, kernel)) == sp.zeros(16)
            for kernel in paulis
            for projector in bell_projectors
        ),
        "a Bell projector is not N64 invariant",
    )
    n2_kernel_commutant_dimension = sum(
        joint_characters.count(character) ** 2 for character in set(joint_characters)
    )
    require(n2_kernel_commutant_dimension == 16, "N64 commutant on V4 tensor V4 changed")

    mixed_joint_characters: list[tuple[sp.Expr, ...]] = []
    mixed_diagonal_actions: list[sp.Matrix] = []
    for kernel in paulis:
        action = sp.kronecker_product(kernel, kernel.conjugate())
        diagonal = clean(bell_basis.H * action * bell_basis)
        require(diagonal == sp.diag(*diagonal.diagonal()), "mixed N64 action is not Bell diagonal")
        require(all(value in (-1, 1) for value in diagonal.diagonal()), "mixed Bell character left {+1,-1}")
        mixed_diagonal_actions.append(diagonal)
    for index in range(16):
        mixed_joint_characters.append(
            tuple(diagonal[index, index] for diagonal in mixed_diagonal_actions)
        )
    require(len(set(mixed_joint_characters)) == 16, "mixed Bell joint characters are not distinct")
    mixed_kernel_commutant_dimension = sum(
        mixed_joint_characters.count(character) ** 2
        for character in set(mixed_joint_characters)
    )
    require(mixed_kernel_commutant_dimension == 16, "mixed N64 commutant dimension changed")

    # Actual root action on those sixteen projectors, compared with the actual
    # ovoid labels from the E8/J doily.  Four native lifts of each target
    # transposition induce one and the same Bell permutation.
    reflection_records: list[dict[str, Any]] = []
    target_to_bell: dict[tuple[int, ...], set[tuple[int, ...]]] = defaultdict(set)
    target_to_mixed_bell: dict[tuple[int, ...], set[tuple[int, ...]]] = defaultdict(set)
    for row in rows:
        reflection = erasure.reflection_from_gaussian_ray(row["gaussian_ray"])
        bell_permutation, phases = pauli_congruence_permutation(reflection, paulis)
        mixed_bell_permutation, mixed_phases = pauli_similarity_permutation(reflection, paulis)
        ordinary_tensor_action = sp.kronecker_product(reflection, reflection)
        mixed_tensor_action = sp.kronecker_product(reflection, reflection.conjugate())
        require(
            all(
                clean(ordinary_tensor_action * bell_vectors[index]
                      - phases[index] * bell_vectors[bell_permutation[index]]) == sp.zeros(16, 1)
                for index in range(16)
            ),
            "congruence permutation does not equal the actual V4 tensor V4 action",
        )
        require(
            all(
                clean(mixed_tensor_action * bell_vectors[index]
                      - mixed_phases[index] * bell_vectors[mixed_bell_permutation[index]]) == sp.zeros(16, 1)
                for index in range(16)
            ),
            "similarity permutation does not equal the actual V4 tensor bar(V4) action",
        )
        target_permutation = tuple(row["ovoid_permutation"])
        require(
            all(parity[bell_permutation[index]] == parity[index] for index in range(16)),
            "root action mixed symmetric and antisymmetric Bell sectors",
        )
        target_to_bell[target_permutation].add(bell_permutation)
        target_to_mixed_bell[target_permutation].add(mixed_bell_permutation)
        reflection_records.append(
            {
                "reflection_index": row["index"],
                "target_cycle_type": list(cycle_type(target_permutation)),
                "bell10_cycle_type": list(cycle_type(restriction(bell_permutation, symmetric_indices))),
                "bell6_cycle_type": list(cycle_type(restriction(bell_permutation, antisymmetric_indices))),
                "mixed_Bell16_cycle_type": list(cycle_type(mixed_bell_permutation)),
                "projective_phases": [str(value) for value in phases],
                "mixed_projective_phases": [str(value) for value in mixed_phases],
            }
        )
    require(len(target_to_bell) == 15, "expected fifteen target transpositions")
    require({len(images) for images in target_to_bell.values()} == {1}, "native lift changed Bell permutation")
    require(
        len(target_to_mixed_bell) == 15
        and {len(images) for images in target_to_mixed_bell.values()} == {1},
        "native lift changed mixed Bell permutation",
    )

    generators = [
        (target, next(iter(images)))
        for target, images in sorted(target_to_bell.items())
    ]
    quotient_group = paired_permutation_group(generators)
    require(len(quotient_group) == 720, "actual root permutations did not generate S6")
    mixed_generators = [
        (target, next(iter(images)))
        for target, images in sorted(target_to_mixed_bell.items())
    ]
    mixed_quotient_group = paired_permutation_group(mixed_generators)
    require(len(mixed_quotient_group) == 720, "mixed Bell root permutations did not generate S6")

    target_character: list[int] = []
    bell10_character: list[int] = []
    bell6_character: list[int] = []
    chi9: list[int] = []
    chi10: list[int] = []
    for target, bell in quotient_group:
        chi = sum(target[index] == index for index in range(6)) - 1
        target_squared = compose(target, target)
        chi_squared_argument = sum(target_squared[index] == index for index in range(6)) - 1
        target_character.append(chi)
        bell10_character.append(sum(bell[index] == index for index in symmetric_indices))
        bell6_character.append(sum(bell[index] == index for index in antisymmetric_indices))
        chi9.append((chi * chi + chi_squared_argument) // 2 - chi - 1)
        chi10.append((chi * chi - chi_squared_argument) // 2)

    ones = [1] * len(quotient_group)
    bell10_centered = [value - 1 for value in bell10_character]
    bell6_centered = [value - 1 for value in bell6_character]
    target_hom_bell10 = character_inner_product(target_character, bell10_character)
    target_hom_bell6 = character_inner_product(target_character, bell6_character)
    require(target_hom_bell10 == target_hom_bell6 == 0, "target five appeared at arity two")
    require(character_inner_product(target_character, target_character) == 1, "target five is reducible")
    require(character_inner_product(bell10_centered, bell10_centered) == 1, "Bell10 is not 1+9")
    require(character_inner_product(bell6_centered, bell6_centered) == 1, "Bell6 is not 1+5")
    require(character_inner_product(target_character, bell6_centered) == 0, "outer five collapsed to target five")

    # Exact End(H5)=1+5+9+10 character identity.  The actual quartic H5 is
    # intertwined with the target standard five by the all-sixty bridge above.
    end5_character = [value * value for value in target_character]
    require(
        all(
            end5_character[index]
            == ones[index] + target_character[index] + chi9[index] + chi10[index]
            for index in range(720)
        ),
        "End(H5) character identity failed",
    )
    end5_components = [ones, target_character, chi9, chi10]
    gram = [
        [character_inner_product(left, right) for right in end5_components]
        for left in end5_components
    ]
    require(gram == sp.eye(4).tolist(), "1+5+9+10 characters are not orthonormal")
    target_hom_end5 = character_inner_product(target_character, end5_character)
    require(target_hom_end5 == 1, "target five multiplicity in End(H5) changed")

    # Orientation control requested separately from V4 tensor V4.  For
    # V4 tensor conjugate(V4), vec(P) transforms by R P R*, not congruence.
    # Its nonzero Pauli-label 15 decomposes as 1+9+outer5 and again contains
    # no actual target five.
    normal_by_target = {target: bell for target, bell in quotient_group}
    mixed_by_target = {target: bell for target, bell in mixed_quotient_group}
    require(len(normal_by_target) == len(mixed_by_target) == 720, "paired quotient projection is not faithful")
    mixed_target_character: list[int] = []
    mixed15_character: list[int] = []
    mixed9_character: list[int] = []
    mixed_outer5_character: list[int] = []
    normal_outer5_character: list[int] = []
    for target, mixed_bell in mixed_by_target.items():
        chi = sum(target[index] == index for index in range(6)) - 1
        target_squared = compose(target, target)
        chi_squared_argument = sum(target_squared[index] == index for index in range(6)) - 1
        nine = (chi * chi + chi_squared_argument) // 2 - chi - 1
        fifteen = sum(mixed_bell[index] == index for index in range(1, 16))
        outer_five = fifteen - 1 - nine
        normal_bell = normal_by_target[target]
        normal_outer_five = sum(normal_bell[index] == index for index in antisymmetric_indices) - 1
        mixed_target_character.append(chi)
        mixed15_character.append(fifteen)
        mixed9_character.append(nine)
        mixed_outer5_character.append(outer_five)
        normal_outer5_character.append(normal_outer_five)
    require(mixed_outer5_character == normal_outer5_character, "mixed residual is not the actual outer five")
    require(character_inner_product(mixed_target_character, mixed15_character) == 0, "mixed Bell15 contains target five")
    require(character_inner_product(mixed_outer5_character, mixed_outer5_character) == 1, "mixed outer five is reducible")
    require(character_inner_product(mixed_target_character, mixed_outer5_character) == 0, "mixed outer five equals target five")

    # Full G31 character count, independently of the sector arguments.  Pick
    # one actual native reflection above each of the fifteen target
    # transpositions, build one native V4 word representative G_pi for every
    # pi in S6, and then enumerate the sixteen projective Pauli cosets.
    native_generators_by_target: dict[tuple[int, ...], Matrix4] = {}
    for row in rows:
        native_generators_by_target.setdefault(
            tuple(row["ovoid_permutation"]),
            fractional_reflection(row["gaussian_ray"]),
        )
    native_generators = sorted(native_generators_by_target.items())
    require(len(native_generators) == 15, "native transposition lift count changed")
    native_representatives = native_coset_representatives(native_generators)
    require(len(native_representatives) == 720, "native words did not cover all S6 permutations")

    # Every one of the sixty original roots lies in the expected Pauli coset.
    all_native_lifts_in_cosets = all(
        matrix4_projective_key(
            matrix4_multiply(
                fractional_reflection(row["gaussian_ray"]),
                matrix4_adjoint(native_representatives[tuple(row["ovoid_permutation"])]),
            )
        )
        in pauli_projective_keys
        for row in rows
    )
    require(all_native_lifts_in_cosets, "an original native root missed its Pauli coset")

    # The chosen lifts normalize the projective Pauli group, and right
    # multiplication by every chosen reflection has the exact expected next
    # coset.  These are the multiplication checks behind the transversal.
    pauli_normality = all(
        matrix4_projective_key(
            matrix4_multiply(
                matrix4_multiply(generator, pauli),
                matrix4_adjoint(generator),
            )
        )
        in pauli_projective_keys
        for _, generator in native_generators
        for pauli in paulis_fractional
    )
    require(pauli_normality, "a native lift failed to normalize the Pauli kernel")
    coset_multiplication = all(
        matrix4_projective_key(
            matrix4_multiply(
                matrix4_multiply(matrix, generator_matrix),
                matrix4_adjoint(
                    native_representatives[compose(permutation, generator_permutation)]
                ),
            )
        )
        in pauli_projective_keys
        for permutation, matrix in native_representatives.items()
        for generator_permutation, generator_matrix in native_generators
    )
    require(coset_multiplication, "native transversal failed a coset multiplication")

    projective_group_keys = {
        matrix4_projective_key(matrix4_multiply(pauli, representative))
        for representative in native_representatives.values()
        for pauli in paulis_fractional
    }
    require(len(projective_group_keys) == 11520, "projective G31 enumeration is not 16 times 720")

    full_character_sums: dict[int, Fraction] = {}
    full_hom_multiplicities: dict[int, Fraction] = {}
    for arity in (1, 2, 3):
        total = Fraction(0)
        for permutation, representative in native_representatives.items():
            target_chi = sum(permutation[index] == index for index in range(6)) - 1
            for pauli in paulis_fractional:
                real, imaginary = matrix4_trace(matrix4_multiply(pauli, representative))
                trace_norm_squared = real * real + imaginary * imaginary
                total += target_chi * trace_norm_squared**arity
        full_character_sums[arity] = total
        full_hom_multiplicities[arity] = total / 11520
    require(
        full_hom_multiplicities
        == {1: Fraction(0), 2: Fraction(0), 3: Fraction(1)},
        "full G31 Hom multiplicities are not (0,0,1)",
    )

    # Arity three: reconstruct the original erasure W and prove that its image
    # is exactly Sym^3(V4), with the original gauge-first bar4 tensor H5 action.
    _, supports = _code_basis()
    W = erasure.decoder_isometry_gauge_first(supports)
    require(clean(W.H * W) == sp.eye(20), "original W is not isometric")
    carrier_projector = clean(W * W.H)
    sym3_projector = symmetric_cube_projector()
    require(carrier_projector == sym3_projector, "image(W) is not Sym^3(V4)")
    require(carrier_projector.rank() == 20, "symmetric-cube rank changed")

    kernel_intertwiners = []
    for pauli in paulis:
        identity = clean(
            tensor_power(pauli, 3) * W
            - W * sp.kronecker_product(pauli.conjugate(), sp.eye(5))
        ) == sp.zeros(64, 20)
        kernel_intertwiners.append(identity)
    require(all(kernel_intertwiners), "N64 action on Sym^3 is not bar4 tensor I5")
    decoded_kernel_commutant_dimension = n1_commutant_dimension * 25
    require(decoded_kernel_commutant_dimension == 25, "decoded N64 commutant is not I4 tensor End5")

    metric = sp.diag(1, 3, 3, 3, 6)
    D = sp.diag(1, sp.sqrt(3), sp.sqrt(3), sp.sqrt(3), sp.sqrt(6))
    ovoid_rays = [sp.Matrix(ray) for ray in bridge["selection_data"]["native_ovoid_rays"]]
    effects_native = [clean(ray * ray.T * metric / 12) for ray in ovoid_rays]
    effects5 = [clean(D * effect * D.inv()) for effect in effects_native]
    require(clean(sum(effects5, sp.zeros(5))) == sp.eye(5), "native H5 effects do not sum to I5")
    require(all(effect.rank() == 1 for effect in effects5), "native H5 effect rank changed")
    require(
        all(clean(effect**2 - sp.Rational(5, 6) * effect) == sp.zeros(5) for effect in effects5),
        "native H5 effect spectrum changed",
    )

    effects20_t1 = [sp.kronecker_product(sp.eye(4), effect) for effect in effects5]
    physical_t1 = [
        clean(W * effect * W.H + (sp.eye(64) - carrier_projector) / 6)
        for effect in effects20_t1
    ]
    require(clean(sum(physical_t1, sp.zeros(64))) == sp.eye(64), "full arity-three effects do not sum to I64")
    require(all(effect.rank() == 48 for effect in physical_t1), "full t=1 effect rank is not 48")
    require(
        all(
            clean((effect - sp.Rational(5, 6) * sp.eye(64))
                  * effect
                  * (effect - sp.Rational(1, 6) * sp.eye(64))) == sp.zeros(64)
            for effect in physical_t1
        ),
        "full t=1 effect spectrum changed",
    )

    decoder_checks = []
    effect_checks = []
    projector_checks = []
    full_effect_checks = []
    for row in rows:
        permutation = tuple(row["ovoid_permutation"])
        reflection = erasure.reflection_from_gaussian_ray(row["gaussian_ray"])
        C_native = sp.Matrix(row["quartic_numerator"]) / 16
        C5 = clean(D * C_native * D.inv())
        G20 = sp.kronecker_product(reflection.conjugate(), C5)
        R3 = tensor_power(reflection, 3)
        decoder_ok = clean(R3 * W - W * G20) == sp.zeros(64, 20)
        effects_ok = all(
            clean(C5 * effects5[index] * C5.H - effects5[permutation[index]]) == sp.zeros(5)
            for index in range(6)
        )
        projector_ok = clean(R3 * carrier_projector * R3.H - carrier_projector) == sp.zeros(64)
        # This is the exact 64-dimensional covariance identity, certified in
        # factored form.  Expanding 360 dense 64 by 64 products adds minutes
        # while checking no additional premise:
        #   R3 W = W G20,
        #   C5 E_a C5* = E_pi(a),
        #   R3 Pi20 R3* = Pi20.
        full_ok = decoder_ok and effects_ok and projector_ok
        require(decoder_ok and effects_ok and projector_ok and full_ok, f"root {row['index']} failed arity-three covariance")
        decoder_checks.append(decoder_ok)
        effect_checks.append(effects_ok)
        projector_checks.append(projector_ok)
        full_effect_checks.append(full_ok)

    # Since the target five occurs once in End(H5), every covariant centered
    # family on Sym^3 is a real multiple of E_a-I5/6.  Positivity gives the two
    # exact eigenvalues below and hence the closed interval.
    t = sp.symbols("t", real=True)
    ray_eigenvalue = sp.simplify((1 + 4 * t) / 6)
    perpendicular_eigenvalue = sp.simplify((1 - t) / 6)
    t_interval = [sp.Rational(-1, 4), sp.S.One]
    endpoint_ranks_on_sym3 = {"-1/4": 16, "1": 4}

    report = {
        "claim": "minimum native tensor arity for a nontrivial six-ovoid covariant POVM is three",
        "status": "PASS",
        "arithmetic": "exact Gaussian integers, rationals, and algebraic square roots",
        "scope": (
            "all linear Hermitian POVM effects covariant for the actual six-ovoid label action; "
            "no physical copying process, tensor-stage selection, t selection, Nfam derivation, or P1 origin"
        ),
        "sources": {
            str(DOILY_SOURCE.relative_to(REPO)): sha256(DOILY_SOURCE),
            str(ERASURE_SOURCE.relative_to(REPO)): sha256(ERASURE_SOURCE),
            str(CORRELATED_SOURCE.relative_to(REPO)): sha256(CORRELATED_SOURCE),
            str(PROCESS_SOURCE.relative_to(REPO)): sha256(PROCESS_SOURCE),
            str(QUARTIC_SOURCE.relative_to(REPO)): sha256(QUARTIC_SOURCE),
        },
        "arity_one": {
            "space": "V4",
            "N64_commutant_dimension": n1_commutant_dimension,
            "conclusion": "all six effects equal I4/6",
        },
        "arity_two": {
            "space": "V4 tensor V4",
            "Bell_projector_count": len(bell_projectors),
            "Bell_projectors_rank": [projector.rank() for projector in bell_projectors],
            "Bell_joint_character_count": len(set(joint_characters)),
            "N64_commutant_dimension": n2_kernel_commutant_dimension,
            "transpose_split": {"symmetric": 10, "antisymmetric": 6},
            "root_cycle_types": {
                "target_ovoid6": sorted({tuple(record["target_cycle_type"]) for record in reflection_records}),
                "Bell_symmetric10": sorted({tuple(record["bell10_cycle_type"]) for record in reflection_records}),
                "Bell_antisymmetric6": sorted({tuple(record["bell6_cycle_type"]) for record in reflection_records}),
            },
            "quotient_group_order": len(quotient_group),
            "Hom_target5_to_Bell10": int(target_hom_bell10),
            "Hom_target5_to_Bell6": int(target_hom_bell6),
            "Hom_target5_to_full_N64_commutant": int(target_hom_bell10 + target_hom_bell6),
            "Bell10_centered_character_norm": int(character_inner_product(bell10_centered, bell10_centered)),
            "Bell6_centered_character_norm": int(character_inner_product(bell6_centered, bell6_centered)),
            "target5_Bell6_centered_inner_product": int(character_inner_product(target_character, bell6_centered)),
            "all_60_root_actions_checked": len(reflection_records),
            "native_lifts_per_target_transposition": 4,
            "conclusion": "the antisymmetric Bell six is the outer six; all six effects equal I16/6",
        },
        "mixed_orientation_control": {
            "space": "V4 tensor conjugate(V4)",
            "separately_reconstructed_action": "vec(P) maps by R P R*",
            "Bell_projector_count": 16,
            "Bell_joint_character_count": len(set(mixed_joint_characters)),
            "N64_commutant_dimension": mixed_kernel_commutant_dimension,
            "nonzero_Pauli_label_count": 15,
            "root_cycle_type_on_Bell16": sorted(
                {tuple(record["mixed_Bell16_cycle_type"]) for record in reflection_records}
            ),
            "quotient_group_order": len(mixed_quotient_group),
            "nonzero15_character_decomposition": "1+9+outer5",
            "Hom_target5_to_mixed_Bell16": int(
                character_inner_product(mixed_target_character, mixed15_character)
            ),
            "target5_outer5_inner_product": int(
                character_inner_product(mixed_target_character, mixed_outer5_character)
            ),
            "conclusion": "the mixed orientation also has only the uniform six-effect Ovoid POVM",
        },
        "full_G31_character_count": {
            "formula": "m_n=(1/11520) sum_{pi in S6,P in Pauli16} chi5(pi) |tr(P G_pi)|^(2n)",
            "native_transposition_lifts_used": len(native_generators),
            "native_word_representatives": len(native_representatives),
            "phase_free_Pauli_representatives": len(paulis_fractional),
            "projectively_distinct_group_elements": len(projective_group_keys),
            "all_60_original_roots_in_expected_cosets": all_native_lifts_in_cosets,
            "Pauli_kernel_normalized_by_generators": pauli_normality,
            "all_720_times_15_coset_multiplications_exact": coset_multiplication,
            "central_mu4_phases": "four phases cancel from |tr|^(2n), so the projective average equals the full G31 average",
            "unnormalized_character_sums": {
                str(arity): int(value) for arity, value in full_character_sums.items()
            },
            "Hom_actual5_to_End_tensor_power": {
                str(arity): int(value) for arity, value in full_hom_multiplicities.items()
            },
        },
        "arity_three": {
            "space": "V4 tensor V4 tensor V4",
            "carrier": "image(W)=Sym^3(V4)=bar(V4) tensor H5",
            "W_shape": [W.rows, W.cols],
            "W_sha256": matrix_sha256(W),
            "carrier_projector_rank": carrier_projector.rank(),
            "carrier_projector_sha256": matrix_sha256(carrier_projector),
            "W_star_W": "I20",
            "W_W_star": "the exact three-slot symmetrizer",
            "all_16_N64_intertwiners_exact": all(kernel_intertwiners),
            "N64_commutant_on_carrier": "I4 tensor End(H5)",
            "N64_commutant_dimension_on_carrier": decoded_kernel_commutant_dimension,
            "End_H5_character_decomposition": "1+5+9+10, multiplicity one each",
            "End_H5_character_gram": [[int(value) for value in row] for row in gram],
            "Hom_target5_to_End_H5": int(target_hom_end5),
            "full_End_V4_tensor3_target5_multiplicity": int(full_hom_multiplicities[3]),
            "classification_on_full_tensor_cube": (
                "I64/6 + t W[I4 tensor (E_a-I5/6)]W*, equivalently "
                "W[I4 tensor (t E_a+(1-t)I5/6)]W* + (I64-WW*)/6"
            ),
            "informative_support": "the unique centered intertwiner is supported on Sym^3(V4); this is forced by full G31 covariance",
            "positivity_eigenvalues_on_H5": {
                "ovoid_ray": str(ray_eigenvalue),
                "orthogonal_four_space": str(perpendicular_eigenvalue),
            },
            "exact_t_interval": [str(value) for value in t_interval],
            "endpoint_ranks_on_Sym3": endpoint_ranks_on_sym3,
            "uniform_complement_extension": "W A_a(t) W* + (I64-WW*)/6",
            "t_equals_one_full_effect_ranks": [effect.rank() for effect in physical_t1],
            "t_equals_one_full_effect_spectrum": {"5/6": 4, "1/6": 44, "0": 16},
            "t_equals_one_effect_sha256": [matrix_sha256(effect) for effect in physical_t1],
            "all_60_decoder_intertwiners_exact": all(decoder_checks),
            "all_60_H5_effect_covariances_exact": all(effect_checks),
            "all_60_carrier_projector_covariances_exact": all(projector_checks),
            "all_60_times_6_full_effect_covariances_exact": all(full_effect_checks),
            "full_covariance_certificate": (
                "factored exact identities R3 W=W G20, C5 E_a C5*=E_pi(a), "
                "and R3 Pi20 R3*=Pi20 for each of 60 roots"
            ),
            "conclusion": "the displayed t-family is every full C64 covariant POVM; t=1 proves attainability",
        },
        "minimum_arity": 3,
    }
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    result = main(arguments.output)
    print(
        "PASS: n=1 scalar; n=2 Bell commutant 10+6 with Hom=0; "
        "full n=3 has the unique t-family supported on Sym^3, -1/4 <= t <= 1."
    )
