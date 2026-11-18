#!/usr/bin/env python3
"""Exact native Sp4 to Spin5 functor and basis-independent Hodge Kraus space."""

from __future__ import annotations

from itertools import combinations, combinations_with_replacement, product
import hashlib
import json
from pathlib import Path

import sympy as sp


ROOT = Path("/Users/stefanhamann/Projekte/tfpt-theoryv4")
V783_PIN = "8f4851634b83f61671b04f3a6211059c40758d21caaf1f779302c799e1e9d6c4"


def require(value: bool, label: str, checks: list[str]) -> None:
    if not bool(value):
        raise AssertionError(label)
    checks.append(label)


def is_zero(matrix: sp.Matrix) -> bool:
    return sp.simplify(matrix) == sp.zeros(*matrix.shape)


def exterior_action(matrix: sp.Matrix, pairs: list[tuple[int, int]]) -> sp.Matrix:
    return sp.Matrix.hstack(*[
        sp.Matrix([
            matrix[i, a] * matrix[j, b] - matrix[j, a] * matrix[i, b]
            for i, j in pairs
        ])
        for a, b in pairs
    ])


def exterior_lie_action(matrix: sp.Matrix, dimension: int,
                        pairs: list[tuple[int, int]]) -> sp.Matrix:
    result = sp.zeros(len(pairs))
    for column, (a, b) in enumerate(pairs):
        for i in range(dimension):
            coefficient = matrix[i, a]
            if coefficient and i != b:
                pair = tuple(sorted((i, b)))
                sign = 1 if i < b else -1
                result[pairs.index(pair), column] += sign * coefficient
            coefficient = matrix[i, b]
            if coefficient and i != a:
                pair = tuple(sorted((a, i)))
                sign = 1 if a < i else -1
                result[pairs.index(pair), column] += sign * coefficient
    return result


def symmetric_action(matrix: sp.Matrix,
                     pairs: list[tuple[int, int]]) -> sp.Matrix:
    columns = []
    for a, b in pairs:
        tensor = sp.zeros(matrix.rows)
        tensor[a, b] = 1
        tensor[b, a] = 1
        image = matrix * tensor * matrix.T
        columns.append(sp.Matrix([image[i, j] for i, j in pairs]))
    return sp.Matrix.hstack(*columns)


def symmetric_lie_action(matrix: sp.Matrix,
                         basis: list[sp.Matrix],
                         pairs: list[tuple[int, int]]) -> sp.Matrix:
    return sp.Matrix.hstack(*[
        sp.Matrix([(matrix * tensor + tensor * matrix.T)[i, j]
                   for i, j in pairs])
        for tensor in basis
    ])


def permutation(images: list[int]) -> sp.Matrix:
    result = sp.zeros(len(images))
    for source, target in enumerate(images):
        result[target, source] = 1
    return result


def encode(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(sp.simplify(value)) for value in row] for row in matrix.tolist()]


def main() -> None:
    checks: list[str] = []
    source = ROOT / "verification/v783_two_qubit_clifford.py"
    require(hashlib.sha256(source.read_bytes()).hexdigest() == V783_PIN,
            "v783 source pin", checks)

    ii = sp.I
    sqrt2 = sp.sqrt(2)
    identity2 = sp.eye(2)
    identity4 = sp.eye(4)
    H = sp.Matrix([[1, 1], [1, -1]]) / sqrt2
    S = sp.diag(1, ii)
    Y = sp.Matrix([[0, -ii], [ii, 0]])
    Z = sp.diag(1, -1)
    CNOT = sp.Matrix([
        [1, 0, 0, 0], [0, 1, 0, 0],
        [0, 0, 0, 1], [0, 0, 1, 0],
    ])
    c_native = sp.simplify(
        -ii * sp.kronecker_product(identity2, H) * CNOT *
        sp.kronecker_product(S, identity2) *
        sp.kronecker_product(H, identity2)
    )
    g_native = sp.simplify(
        sp.kronecker_product(Y, Z) * CNOT *
        sp.kronecker_product(identity2, S) *
        sp.kronecker_product(S, identity2)
    )
    omega = sp.Matrix([
        [0, 1, 0, 0], [-1, 0, 0, 0],
        [0, 0, 0, 1], [0, 0, -1, 0],
    ])
    require(c_native.T * omega * c_native == omega and
            g_native.T * omega * g_native == ii * omega,
            "native F20 is conformal symplectic with multipliers one and i", checks)

    pairs4 = list(combinations(range(4), 2))
    pairs5 = list(combinations(range(5), 2))
    symmetric_pairs = list(combinations_with_replacement(range(4), 2))

    # W is the primitive complement to the inverse symplectic bivector in
    # Lambda2(V). The columns are Hermitian orthonormal.
    primitive_basis = sp.Matrix.hstack(
        sp.Matrix([1/sqrt2, 0, 0, 0, 0, -1/sqrt2]),
        sp.eye(6)[:, 1], sp.eye(6)[:, 2],
        sp.eye(6)[:, 3], sp.eye(6)[:, 4],
    )
    primitive_left_inverse = primitive_basis.conjugate().T
    wedge_pairing6 = sp.zeros(6)
    wedge_pairing6[0, 5] = wedge_pairing6[5, 0] = 1
    wedge_pairing6[1, 4] = wedge_pairing6[4, 1] = -1
    wedge_pairing6[2, 3] = wedge_pairing6[3, 2] = 1
    q5 = sp.simplify(primitive_basis.T * wedge_pairing6 * primitive_basis)
    require(q5.T == q5 and q5.conjugate().T * q5 == sp.eye(5) and q5**2 == sp.eye(5),
            "primitive W has compatible unitary symmetric form", checks)

    # chi(C)=1, chi(G)=i. The chi^-1 twist makes the primitive five orthogonal.
    w_c = sp.simplify(primitive_left_inverse * exterior_action(c_native, pairs4) *
                      primitive_basis)
    w_g = sp.simplify(-ii * primitive_left_inverse *
                      exterior_action(g_native, pairs4) * primitive_basis)
    require(w_c.conjugate().T * w_c == sp.eye(5) and
            w_g.conjugate().T * w_g == sp.eye(5),
            "W5 F20 action is Hermitian unitary", checks)
    require(w_c.T * q5 * w_c == q5 and w_g.T * q5 * w_g == q5,
            "W5 F20 action preserves q", checks)
    require(is_zero(q5 * w_c.conjugate() - w_c * q5) and
            is_zero(q5 * w_g.conjugate() - w_g * q5),
            "W5 action preserves the canonical real form j equals q conjugation", checks)

    p_c = permutation([(j + 1) % 5 for j in range(5)])
    p_g = permutation([(2 * j) % 5 for j in range(5)])
    sign_permutation_g = -p_g
    character_rows = []
    for b, a in product(range(4), range(5)):
        native_trace = sp.simplify(sp.trace(w_c**a * w_g**b))
        slot_trace = sp.trace(p_c**a * sign_permutation_g**b)
        require(native_trace == slot_trace,
                f"W5 equals sign tensor Perm5 character a={a} b={b}", checks)
        character_rows.append({"a": a, "b": b, "trace": str(native_trace)})

    # The slot chart itself is not unique under F20 because sign*Perm5=1+4.
    variables5 = sp.symbols("l0:25")
    slot_map = sp.Matrix(5, 5, variables5)
    equations5 = list(slot_map * w_c - p_c * slot_map)
    equations5 += list(slot_map * w_g - sign_permutation_g * slot_map)
    coefficient5, _ = sp.linear_eq_to_matrix(equations5, variables5)
    require(25 - coefficient5.rank() == 2,
            "F20 slot intertwiner space has dimension two", checks)

    # Canonical Spin5 map. For a symmetric tensor S, A=S Omega lies in sp4.
    # Let L_A be its action on W. The bivector corresponding to L_A is L_A q.
    symmetric_basis = []
    for a, b in symmetric_pairs:
        tensor = sp.zeros(4)
        tensor[a, b] = 1
        tensor[b, a] = 1
        symmetric_basis.append(tensor)
    spin5_columns = []
    source_lie_actions = []
    target_lie_actions = []
    for tensor in symmetric_basis:
        lie_element = tensor * omega
        lambda2_lie = exterior_lie_action(lie_element, 4, pairs4)
        w_lie = sp.simplify(primitive_left_inverse * lambda2_lie * primitive_basis)
        require(is_zero(w_lie.T * q5 + q5 * w_lie),
                "induced W endomorphism lies in so5", checks)
        bivector = sp.simplify(w_lie * q5)
        require(bivector.T == -bivector,
                "so5 endomorphism raises to an alternating bivector", checks)
        spin5_columns.append(sp.Matrix([bivector[i, j] for i, j in pairs5]))
        source_lie_actions.append(
            symmetric_lie_action(lie_element, symmetric_basis, symmetric_pairs)
        )
        target_lie_actions.append(exterior_lie_action(w_lie, 5, pairs5))
    spin5_map = sp.Matrix.hstack(*spin5_columns)
    source_metric = sp.diag(*[1 if a == b else 2 for a, b in symmetric_pairs])
    require(spin5_map.det() == 8, "canonical Spin5 map determinant eight", checks)
    require(spin5_map.conjugate().T * spin5_map == source_metric,
            "canonical Spin5 map is exactly metric preserving", checks)
    require(all(is_zero(spin5_map * source - target * spin5_map)
                for source, target in zip(source_lie_actions, target_lie_actions)),
            "canonical map intertwines the full sp4 Lie algebra", checks)

    # This is the decisive uniqueness test: F20 alone gave dimension six on
    # the ten, while full sp4 naturality leaves a single complex line.
    variables10 = sp.symbols("u0:100")
    unknown10 = sp.Matrix(10, 10, variables10)
    equations10 = []
    for source_action, target_action in zip(source_lie_actions, target_lie_actions):
        equations10 += list(unknown10 * source_action - target_action * unknown10)
    coefficient10, _ = sp.linear_eq_to_matrix(equations10, variables10)
    require(100 - coefficient10.rank() == 1,
            "full sp4 natural intertwiner space has dimension one", checks)

    # Check the finite native F20 equivariance of the same canonical map.
    source_c10 = symmetric_action(c_native, symmetric_pairs)
    source_g10 = -ii * symmetric_action(g_native, symmetric_pairs)
    target_c10 = exterior_action(w_c, pairs5)
    target_g10 = exterior_action(w_g, pairs5)
    require(is_zero(spin5_map * source_c10 - target_c10 * spin5_map) and
            is_zero(spin5_map * source_g10 - target_g10 * spin5_map),
            "canonical Spin5 map respects the native F20 actions", checks)

    # Hodge complement Kraus maps on an oriented orthonormal real basis of W.
    # Individual matrices depend on the basis; their five-dimensional span and
    # the normalized channel do not.
    kraus = []
    for m in range(5):
        operator = sp.zeros(10)
        for column, (a, b) in enumerate(pairs5):
            if m in (a, b):
                continue
            complement = [j for j in range(5) if j not in (m, a, b)]
            sequence = [m, a, b] + complement
            inversions = sum(sequence[i] > sequence[j]
                             for i in range(5) for j in range(i + 1, 5))
            operator[pairs5.index(tuple(complement)), column] = (-1) ** inversions
        kraus.append(operator)
    completeness = sum((operator.conjugate().T * operator for operator in kraus),
                       sp.zeros(10))
    hilbert_schmidt_gram = sp.Matrix([
        [sp.trace(left.conjugate().T * right) for right in kraus]
        for left in kraus
    ])
    kraus_vectors = sp.Matrix.hstack(*[
        sp.Matrix(100, 1, list(operator)) for operator in kraus
    ])
    require(completeness == 3 * sp.eye(10),
            "Hodge Kraus completeness sum is three I", checks)
    require(hilbert_schmidt_gram == 6 * sp.eye(5) and kraus_vectors.rank() == 5,
            "Hodge Kraus space has canonical dimension five", checks)

    output = {
        "verdict": "CANONICAL_NATIVE_SPIN5_FUNCTOR_AND_KRAUS_SPACE",
        "W5": "primitive Lambda2(V) tensor chi^(-1)",
        "W5_q": encode(q5),
        "W5_real_structure": "j(w)=q conjugate(w)",
        "W5_F20": "sign tensor Perm5",
        "W5_character": character_rows,
        "F20_slot_intertwiner_dimension": 2,
        "spin5_map": encode(spin5_map),
        "spin5_map_determinant": str(spin5_map.det()),
        "spin5_metric_identity": "F^dagger F=diag(1,2,2,2,1,2,2,1,2,1)",
        "full_sp4_intertwiner_dimension": 1,
        "metric_residual": "U(1) phase; real compact form leaves a sign",
        "kraus_space_dimension": 5,
        "kraus_hilbert_schmidt_gram": "6 I5",
        "kraus_completeness": "sum C_m^dagger C_m=3 I10",
        "normalized_channel": "Phi(rho)=(1/3) sum C_m rho C_m^dagger",
        "checks": checks,
        "scope": (
            "Full Sp4 functoriality selects a unique natural ten-intertwiner "
            "line and a basis-independent Hodge Kraus space on the native "
            "compact real form. F20 alone still leaves a two-dimensional "
            "slot-chart Hom because sign tensor Perm5 decomposes as 1 plus 4. "
            "This constructs a canonical mathematical channel; it does not "
            "prove that the original physical source uses that channel."
        ),
    }
    out = Path(__file__).with_name("native-spin5-functor-results.json")
    out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "verdict": output["verdict"],
        "checks": len(checks),
        "sp4_hom_dimension": output["full_sp4_intertwiner_dimension"],
        "kraus_dimension": output["kraus_space_dimension"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
