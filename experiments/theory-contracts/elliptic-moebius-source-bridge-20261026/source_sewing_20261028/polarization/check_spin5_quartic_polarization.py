#!/usr/bin/env python3
"""Exact Spin5 Pluecker relation versus the original five quartics."""

from __future__ import annotations

import argparse
from itertools import combinations, combinations_with_replacement
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]
ROOT = Path("/Users/stefanhamann/Projekte/tfpt-theoryv4")
PROCESS = ROOT / "tfpt_explorer/process.py"
FUNCTOR = ROOT / "experiments/theory-contracts/elliptic-moebius-source-bridge-20261026/source_recursion_20261028/check_native_spin5_functor.py"
FUNCTOR_REPORT = ROOT / "experiments/theory-contracts/elliptic-moebius-source-bridge-20261026/source_recursion_20261028/NATIVE-SPIN5-FUNCTOR-AND-KRAUS-SPACE.txt"
OCTIC_CHECKER = ROOT / "experiments/theory-contracts/native-octic-seam-origin-20261025/native_folded_octic_origin.py"
CAR_REPORT = ROOT / "experiments/theory-contracts/native-octic-seam-origin-20261025/CAR-OCTIC-DICTIONARY.txt"
OPE_REPORT = ROOT / "experiments/theory-contracts/joint-native-algebra-selection-20261016/FULL-SOURCE-HULL-AND-QUARTIC-OPE.txt"
CURRENT_REPORT = ROOT / "experiments/theory-contracts/joint-native-algebra-selection-20261016/ORIGINAL-ORDERED-CURRENT-KERNEL-AND-SEAM-GRADING.txt"
PINS = {
    str(PROCESS): "e8db5fcd999c9d4c34fbee548aa2c85e1a515cf03663a62355533439ae0f1c94",
    str(FUNCTOR): "86911d6e7634592fe9e61d232535224e9b8148308433a05a8ccf983e5a110fa2",
    str(FUNCTOR_REPORT): "3ffcad772d3f0a8fef31ff000cb716b18668f83d013ebe3a55b055d498080050",
    str(OCTIC_CHECKER): "6b075e8d618d392b0d9846bef69c953a03690cd2ce7af76cac003ca49dd072c9",
    str(CAR_REPORT): "c82fce3ad57904fb7628a1abaf2c6eacd94d230d73f2ac69d8da0fe9d49f22c2",
    str(OPE_REPORT): "1d949c2a77e996e63446255d13f41856a00cde4b4d57c78b4b7766204870b2a0",
    str(CURRENT_REPORT): "62a9eeb4519b59d4514445f9d4dd5d2059e29179956a4c97ec94b7bae9ae3094",
}


checks: list[str] = []


def require(value: bool, label: str) -> None:
    if not bool(value):
        raise AssertionError(label)
    checks.append(label)


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"import specification {name}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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
                result[pairs.index(pair), column] += coefficient * (1 if i < b else -1)
            coefficient = matrix[i, b]
            if coefficient and i != a:
                pair = tuple(sorted((a, i)))
                result[pairs.index(pair), column] += coefficient * (1 if a < i else -1)
    return result


def symmetric_lie_action(matrix: sp.Matrix, basis: list[sp.Matrix],
                         pairs: list[tuple[int, int]]) -> sp.Matrix:
    return sp.Matrix.hstack(*[
        sp.Matrix([(matrix * tensor + tensor * matrix.T)[i, j] for i, j in pairs])
        for tensor in basis
    ])


def hodge_amplitudes() -> list[sp.Matrix]:
    pairs = list(combinations(range(5), 2))
    result = []
    for m in range(5):
        operator = sp.zeros(10)
        for column, (a, b) in enumerate(pairs):
            if m in (a, b):
                continue
            complement = [j for j in range(5) if j not in (m, a, b)]
            sequence = [m, a, b] + complement
            inversions = sum(sequence[i] > sequence[j]
                             for i in range(5) for j in range(i + 1, 5))
            operator[pairs.index(tuple(complement)), column] = (-1) ** inversions
        result.append(operator)
    return result


def symmetric_vector(matrix: sp.Matrix) -> sp.Matrix:
    return sp.Matrix([matrix[a, b] for a, b in
                      combinations_with_replacement(range(matrix.rows), 2)])


def quadratic_coefficients(matrix: sp.Matrix) -> sp.Matrix:
    """Coefficients of q^T matrix q in the monomial basis q_a q_b, a<=b."""
    return sp.Matrix([
        matrix[a, a] if a == b else 2 * matrix[a, b]
        for a, b in combinations_with_replacement(range(matrix.rows), 2)
    ])


def matrix_from_quadratic(polynomial: sp.Expr,
                          variables: tuple[sp.Symbol, ...]) -> sp.Matrix:
    poly = sp.Poly(sp.expand(polynomial), *variables)
    matrix = sp.zeros(len(variables))
    for a in range(len(variables)):
        matrix[a, a] = poly.coeff_monomial(variables[a] ** 2)
        for b in range(a + 1, len(variables)):
            value = poly.coeff_monomial(variables[a] * variables[b]) / 2
            matrix[a, b] = matrix[b, a] = value
    require(sp.expand((sp.Matrix(variables).T * matrix * sp.Matrix(variables))[0] - polynomial) == 0,
            "quadratic coefficient matrix reconstructs its polynomial")
    return matrix


def multiindex(pair: tuple[int, int], dimension: int = 4) -> tuple[int, ...]:
    return tuple(pair.count(j) for j in range(dimension))


def factorial_multi(exponents: tuple[int, ...]) -> int:
    return math.prod(math.factorial(value) for value in exponents)


def encode(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(sp.simplify(value)) for value in row] for row in matrix.tolist()]


def main(output: Path | None) -> dict:
    checks.clear()
    for path, digest in PINS.items():
        require(hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest,
                "source pin " + path)

    process = load("pinned_original_process_quartics", PROCESS)
    octic = load("pinned_original_car_transport", OCTIC_CHECKER)
    _, word_sets = process._code_basis()
    require([len(words) for words in word_sets] == [4, 12, 12, 12, 24],
            "original process supplies the five stated quartic word sets")

    pairs4 = list(combinations(range(4), 2))
    pairs5 = list(combinations(range(5), 2))
    sym_pairs = list(combinations_with_replacement(range(4), 2))
    sqrt2 = sp.sqrt(2)
    ii = sp.I
    omega = sp.Matrix([[0, 1, 0, 0], [-1, 0, 0, 0],
                       [0, 0, 0, 1], [0, 0, -1, 0]])
    primitive = sp.Matrix.hstack(
        sp.Matrix([1/sqrt2, 0, 0, 0, 0, -1/sqrt2]),
        sp.eye(6)[:, 1], sp.eye(6)[:, 2], sp.eye(6)[:, 3], sp.eye(6)[:, 4],
    )
    wedge_pairing = sp.zeros(6)
    wedge_pairing[0, 5] = wedge_pairing[5, 0] = 1
    wedge_pairing[1, 4] = wedge_pairing[4, 1] = -1
    wedge_pairing[2, 3] = wedge_pairing[3, 2] = 1
    q5 = sp.simplify(primitive.T * wedge_pairing * primitive)
    real_frame = sp.Matrix.hstack(
        ii * sp.eye(5)[:, 0],
        (sp.eye(5)[:, 1] - sp.eye(5)[:, 4]) / sqrt2,
        ii * (sp.eye(5)[:, 1] + sp.eye(5)[:, 4]) / sqrt2,
        (sp.eye(5)[:, 2] + sp.eye(5)[:, 3]) / sqrt2,
        ii * (sp.eye(5)[:, 2] - sp.eye(5)[:, 3]) / sqrt2,
    )
    require(real_frame.conjugate().T * real_frame == sp.eye(5) and
            q5 * real_frame.conjugate() == real_frame,
            "metric compact real frame is explicit")

    sym_basis = []
    spin5_columns = []
    for a, b in sym_pairs:
        tensor = sp.zeros(4)
        tensor[a, b] = 1
        tensor[b, a] = 1
        sym_basis.append(tensor)
        lie = tensor * omega
        w_lie = sp.simplify(primitive.conjugate().T *
                            exterior_lie_action(lie, 4, pairs4) * primitive)
        bivector = sp.simplify(w_lie * q5)
        spin5_columns.append(sp.Matrix([bivector[i, j] for i, j in pairs5]))
    # Rebuild Lie actions after the complete source basis exists.
    source_lie_actions = [
        symmetric_lie_action(tensor * omega, sym_basis, sym_pairs)
        for tensor in sym_basis
    ]
    spin5_map = sp.Matrix.hstack(*spin5_columns)
    source_metric = sp.diag(*[1 if a == b else 2 for a, b in sym_pairs])
    require(spin5_map.conjugate().T * spin5_map == source_metric,
            "Spin5 map has the actual Sym2 metric")
    spin5_real = sp.simplify(exterior_action(real_frame, pairs5).conjugate().T * spin5_map)

    # The five congruences are bilinear forms on Sym2(V), not channel
    # conjugations on End(Bell10).
    hodge = hodge_amplitudes()
    relations_raw = [sp.simplify(spin5_real.T * operator * spin5_real)
                     for operator in hodge]
    require(all(matrix.T == matrix for matrix in relations_raw),
            "all five F transpose C F congruences are symmetric bilinear forms")

    # q_A=z^A/sqrt(A!) is the normalized Sym2 Fock basis used by the original
    # 2+2 quartic flattenings.
    pair_factorials = [factorial_multi(multiindex(pair)) for pair in sym_pairs]
    raw_from_normalized = sp.diag(*[sp.sqrt(value) for value in pair_factorials])
    relations = [sp.simplify(raw_from_normalized * matrix * raw_from_normalized)
                 for matrix in relations_raw]
    relation_vectors = sp.Matrix.hstack(*[symmetric_vector(matrix) for matrix in relations])
    require(relation_vectors.rank() == 5, "Spin5 congruences span a five-dimensional relation space")

    # Multiplication Sym2(Sym2 V)->Sym4 V in normalized monomial bases.
    domain_pairs = list(combinations_with_replacement(range(10), 2))
    quartic_tuples = list(combinations_with_replacement(range(4), 4))
    quartic_exponents = [multiindex(item, 4) for item in quartic_tuples]
    multiplication = sp.zeros(35, 55)
    for column, (a, b) in enumerate(domain_pairs):
        exponent = tuple(x + y for x, y in zip(multiindex(sym_pairs[a]),
                                                multiindex(sym_pairs[b])))
        row = quartic_exponents.index(exponent)
        multiplication[row, column] = sp.sqrt(
            sp.Rational(factorial_multi(exponent), pair_factorials[a] * pair_factorials[b])
        )
    require(multiplication.rank() == 35 and len(multiplication.nullspace()) == 20,
            "Sym2 Sym2 multiplication has image 35 and Schur22 kernel 20")
    relation_coefficients = sp.Matrix.hstack(*[quadratic_coefficients(matrix)
                                               for matrix in relations])
    require(multiplication * relation_coefficients == sp.zeros(35, 5),
            "all five Spin5 congruences lie in the Schur22 multiplication kernel")

    z = sp.symbols("z0:4")
    rank_one_sym2 = sp.Matrix([
        z[a] * z[b] / sp.sqrt(pair_factorials[index])
        for index, (a, b) in enumerate(sym_pairs)
    ])
    restricted = [sp.expand((rank_one_sym2.T * matrix * rank_one_sym2)[0])
                  for matrix in relations]
    require(restricted == [0] * 5,
            "all five Pluecker quadrics vanish identically on S equals z z transpose")

    # The scalar relation tr((S Omega)^2) supplies the Sp4 singlet in the
    # same 20-dimensional kernel. The five is invariant under all ten sp4
    # generators and independent of that line, leaving dimension 14.
    qvars = sp.symbols("q0:10")
    raw_coordinates = raw_from_normalized * sp.Matrix(qvars)
    S_matrix = sp.zeros(4)
    for coordinate, (a, b) in zip(raw_coordinates, sym_pairs):
        S_matrix[a, b] = S_matrix[b, a] = coordinate
    scalar_relation = matrix_from_quadratic(sp.trace((S_matrix * omega) ** 2), qvars)
    require(multiplication * quadratic_coefficients(scalar_relation) == sp.zeros(35, 1),
            "symplectic scalar relation lies in the Schur22 kernel")
    require(relation_vectors.row_join(symmetric_vector(scalar_relation)).rank() == 6,
            "Spin5 five is independent of the symplectic scalar relation")
    for lie_raw in source_lie_actions:
        lie = sp.simplify(raw_from_normalized.inv() * lie_raw * raw_from_normalized)
        require(sp.simplify(lie.T * scalar_relation + scalar_relation * lie) == sp.zeros(10),
                "scalar relation is sp4 invariant")
        for relation in relations:
            image = sp.simplify(lie.T * relation + relation * lie)
            require(relation_vectors.row_join(symmetric_vector(image)).rank() == 5,
                    "Spin5 relation five is sp4 invariant")

    # Reconstruct the literal original quartics from process._code_basis.
    original_from_words = []
    for words in word_sets:
        original_from_words.append(sp.expand(sum(sp.prod(z[j] for j in word) for word in words)))
    P = [
        sum(value**4 for value in z),
        z[0]**2*z[1]**2 + z[2]**2*z[3]**2,
        z[0]**2*z[2]**2 + z[1]**2*z[3]**2,
        z[0]**2*z[3]**2 + z[1]**2*z[2]**2,
        sp.prod(z),
    ]
    require(original_from_words == [P[0], 6*P[1], 6*P[2], 6*P[3], 24*P[4]],
            "literal process word sets give the five original quartic polynomials")
    normalizations = [96, 8, 8, 8, 1]
    original_hessians = []
    for polynomial, normalization in zip(P, normalizations):
        normalized = sp.Poly(polynomial / sp.sqrt(normalization), *z)
        hessian = sp.zeros(10)
        for a, exponent_a in enumerate(map(multiindex, sym_pairs)):
            for b, exponent_b in enumerate(map(multiindex, sym_pairs)):
                total = tuple(x + y for x, y in zip(exponent_a, exponent_b))
                derivative = normalized.coeff_monomial(sp.prod(z[j]**total[j] for j in range(4)))
                derivative *= factorial_multi(total)
                hessian[a, b] = sp.simplify(
                    derivative / sp.sqrt(6 * pair_factorials[a] * pair_factorials[b])
                )
        original_hessians.append(hessian)
    hessian_vectors = sp.Matrix.hstack(*[symmetric_vector(matrix) for matrix in original_hessians])
    require(hessian_vectors.rank() == 5, "original quartic 2+2 Hessians span dimension five")
    require(relation_vectors.row_join(hessian_vectors).rank() == 10,
            "Spin5 relation five and original quartic Hessian five have zero intersection")
    hessian_images = sp.Matrix.hstack(*[
        multiplication * quadratic_coefficients(matrix) for matrix in original_hessians
    ])
    require(hessian_images.rank() == 5,
            "all five original Hessians have independent nonzero Sym4 images")

    # Exact original CAR transport. AA^T=0 means the positive holomorphic
    # currents have no mutual singular OPE, so their normal product preserves
    # every polynomial identity, including the five zero relations above.
    _, orthogonal, positive, fold = octic.original_transport()
    source_to_positive = sp.simplify(positive.T * orthogonal.T * fold)
    require(sp.simplify(source_to_positive * source_to_positive.conjugate().T) == sp.eye(4),
            "original CAR transport retains the positive metric")
    require(sp.simplify(source_to_positive * source_to_positive.T) == sp.zeros(4),
            "original positive current OPE is mutually regular")

    result = {
        "verdict": [
            "SPIN5_HODGE_FIVE_IS_A_PLUECKER_RELATION_FIVE",
            "ALL_FIVE_CONGRUENCES_VANISH_ON_RANK_ONE_SYM2_TENSORS",
            "ORIGINAL_NATIVE_QUARTIC_HESSIAN_FIVE_HAS_ZERO_INTERSECTION_WITH_K5",
            "ORIGINAL_HOLOMORPHIC_POINT_SPLITTING_CANNOT_LIFT_THE_ZERO_RELATION",
            "MIXED_OR_ASYMMETRIC_DERIVATIVE_OPE_REQUIRES_AN_ADDITIONAL_TYPED_CHOICE",
        ],
        "typed_map": "D_m=F_real^T C_m F_real in Sym2(Sym2 V)^*",
        "multiplication": {
            "domain_dimension": 55,
            "sym4_image_dimension": multiplication.rank(),
            "schur22_kernel_dimension": len(multiplication.nullspace()),
            "spin5_relation_dimension": relation_vectors.rank(),
            "sp4_invariant_scalar_dimension": 1,
            "sp4_invariant_spin5_relation_dimension": 5,
            "remaining_kernel_quotient_dimension": 14,
            "rank_one_restrictions": [str(value) for value in restricted],
        },
        "original_quartics": {
            "polynomials": [str(value) for value in P],
            "fock_normalizations": normalizations,
            "hessian_span_dimension": hessian_vectors.rank(),
            "hessian_sym4_image_rank": hessian_images.rank(),
            "intersection_with_spin5_relation_five": 0,
        },
        "original_point_splitting": {
            "source_to_positive_isometry": "A A^dagger=I4",
            "positive_isotropy": "A A^T=0",
            "consequence": (
                "The four positive holomorphic currents have no mutual singular contractions. "
                "Their bosonic-current normal products preserve the five polynomial zero identities. "
                "The retained fermionic point splitting and its finite counterterms realize the same "
                "zero field; they do not turn K5 into the original Sym4 quartics."
            ),
            "remaining_typed_route": (
                "A mixed adjoint OPE or a nonuniform derivative assignment can leave the holomorphic "
                "multiplication map, but selecting its pole, Taylor order and conjugate leg is extra "
                "data. No source-reviewed original identity currently identifies such a component "
                "with the canonical Spin5 K5."
            ),
        },
        "scope": (
            "This excludes ordinary quartic polarization, the original 2+2 Hessian five, and the "
            "unchanged holomorphic CAR point-splitting realization as sources of the canonical K5. "
            "It is not a no-go theorem for all mixed, derivative or separated OPE components."
        ),
        "source_sha256": PINS,
        "checks": len(checks),
    }
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    print(json.dumps(main(args.out), indent=2, sort_keys=True))
