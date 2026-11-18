#!/usr/bin/env python3
"""Exact native Spin5 Hodge space against the actual one/two source records."""

from __future__ import annotations

import argparse
from itertools import combinations, combinations_with_replacement
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = Path("/Users/stefanhamann/Projekte/tfpt-theoryv4")
SOURCE = ROOT / "experiments/theory-contracts/compiler-correlated-event-clock-20260919/source_channel.py"
FUNCTOR = HERE / "check_native_spin5_functor.py"
FUNCTOR_REPORT = HERE / "NATIVE-SPIN5-FUNCTOR-AND-KRAUS-SPACE.txt"
PINS = {
    str(SOURCE): "3852fabaf10c143fe20238b67406df6396f85c931dfd4834fe91a2d7ef900593",
    str(FUNCTOR): "86911d6e7634592fe9e61d232535224e9b8148308433a05a8ccf983e5a110fa2",
    str(FUNCTOR_REPORT): "3ffcad772d3f0a8fef31ff000cb716b18668f83d013ebe3a55b055d498080050",
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


def exact_gaussian(array: np.ndarray) -> sp.Matrix:
    if not np.allclose(array.real, np.rint(array.real), atol=0, rtol=0):
        raise AssertionError("Gaussian integer real coordinates")
    if not np.allclose(array.imag, np.rint(array.imag), atol=0, rtol=0):
        raise AssertionError("Gaussian integer imaginary coordinates")
    return sp.Matrix([
        [sp.Integer(int(round(z.real))) + sp.I * sp.Integer(int(round(z.imag))) for z in row]
        for row in array
    ])


def vec(matrix: sp.Matrix) -> sp.Matrix:
    return sp.Matrix([matrix[i, j] for j in range(matrix.cols) for i in range(matrix.rows)])


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


def pauli_basis(bits) -> list[sp.Matrix]:
    identity = sp.eye(2)
    x = sp.Matrix([[0, 1], [1, 0]])
    y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    z = sp.diag(1, -1)
    one = {(0, 0): identity, (1, 0): x, (0, 1): z, (1, 1): y}
    return [sp.kronecker_product(one[(a, b)], one[(c, d)]) for a, b, c, d in bits]


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


def encode(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(sp.simplify(value)) for value in row] for row in matrix.tolist()]


def main(output: Path | None) -> dict:
    checks.clear()
    for path, digest in PINS.items():
        require(hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest,
                "source pin " + path)

    old = load("pinned_source_channel_native_spin5", SOURCE)
    pairs4 = list(combinations(range(4), 2))
    pairs5 = list(combinations(range(5), 2))
    symmetric_pairs = list(combinations_with_replacement(range(4), 2))
    sqrt2 = sp.sqrt(2)
    ii = sp.I

    # Frozen native V4 and its conformal symplectic F20 generators.
    identity2 = sp.eye(2)
    H = sp.Matrix([[1, 1], [1, -1]]) / sqrt2
    S = sp.diag(1, ii)
    Y = sp.Matrix([[0, -ii], [ii, 0]])
    Z = sp.diag(1, -1)
    CNOT = sp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0],
                      [0, 0, 0, 1], [0, 0, 1, 0]])
    c_native = sp.simplify(-ii * sp.kronecker_product(identity2, H) * CNOT *
                           sp.kronecker_product(S, identity2) *
                           sp.kronecker_product(H, identity2))
    g_native = sp.simplify(sp.kronecker_product(Y, Z) * CNOT *
                           sp.kronecker_product(identity2, S) *
                           sp.kronecker_product(S, identity2))
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

    # Explicit orthonormal frame of Fix(j), j(w)=q conjugate(w).
    real_frame = sp.Matrix.hstack(
        ii * sp.eye(5)[:, 0],
        (sp.eye(5)[:, 1] - sp.eye(5)[:, 4]) / sqrt2,
        ii * (sp.eye(5)[:, 1] + sp.eye(5)[:, 4]) / sqrt2,
        (sp.eye(5)[:, 2] + sp.eye(5)[:, 3]) / sqrt2,
        ii * (sp.eye(5)[:, 2] - sp.eye(5)[:, 3]) / sqrt2,
    )
    require(real_frame.conjugate().T * real_frame == sp.eye(5),
            "compact real frame is Hermitian orthonormal")
    require(q5 * real_frame.conjugate() == real_frame,
            "compact real frame is fixed by j")

    symmetric_basis = []
    spin5_columns = []
    for a, b in symmetric_pairs:
        tensor = sp.zeros(4)
        tensor[a, b] = 1
        tensor[b, a] = 1
        symmetric_basis.append(tensor)
        lie = tensor * omega
        w_lie = sp.simplify(primitive.conjugate().T *
                            exterior_lie_action(lie, 4, pairs4) * primitive)
        bivector = sp.simplify(w_lie * q5)
        spin5_columns.append(sp.Matrix([bivector[i, j] for i, j in pairs5]))
    spin5_map = sp.Matrix.hstack(*spin5_columns)
    symmetric_metric = sp.diag(*[1 if a == b else 2 for a, b in symmetric_pairs])
    require(spin5_map.conjugate().T * spin5_map == symmetric_metric,
            "Spin5 map preserves the induced Sym2 metric")

    # The source Bell basis consists of symmetric Pauli matrices of HS norm 2.
    # Dividing their Sym2 coordinate columns by two is the required isometry.
    paulis = pauli_basis(old.SYMMETRIC_BITS)
    bell_to_sym = sp.Matrix.hstack(*[
        sp.Matrix([pauli[a, b] / 2 for a, b in symmetric_pairs]) for pauli in paulis
    ])
    require(bell_to_sym.conjugate().T * symmetric_metric * bell_to_sym == sp.eye(10),
            "Bell to Sym2 chart includes the exact Hilbert Schmidt factor one half")
    wedge_real = exterior_action(real_frame, pairs5)
    native_chart = sp.simplify(wedge_real.conjugate().T * spin5_map * bell_to_sym)
    require(native_chart.conjugate().T * native_chart == sp.eye(10),
            "native Bell10 to compact Spin5 wedge chart is unitary")

    # First generator test.  The chi^-1 twist on G is part of the functor.
    w_c = sp.simplify(primitive.conjugate().T * exterior_action(c_native, pairs4) * primitive)
    w_g = sp.simplify(-ii * primitive.conjugate().T *
                      exterior_action(g_native, pairs4) * primitive)
    w_c_real = sp.simplify(real_frame.conjugate().T * w_c * real_frame)
    w_g_real = sp.simplify(real_frame.conjugate().T * w_g * real_frame)
    require(w_c_real.conjugate() == w_c_real and w_g_real.conjugate() == w_g_real,
            "native F20 five acts in the compact real frame")
    source_c = sp.simplify(bell_to_sym.inv() * symmetric_action(c_native, symmetric_pairs) *
                           bell_to_sym)
    source_g = sp.simplify(-ii * bell_to_sym.inv() *
                           symmetric_action(g_native, symmetric_pairs) * bell_to_sym)
    target_c = exterior_action(w_c_real, pairs5)
    target_g = exterior_action(w_g_real, pairs5)
    require(sp.simplify(native_chart * source_c - target_c * native_chart) == sp.zeros(10) and
            sp.simplify(native_chart * source_g - target_g * native_chart) == sp.zeros(10),
            "native chart intertwines both exact F20 generators")

    standard_hodge = hodge_amplitudes()
    native_hodge = [sp.simplify(native_chart.conjugate().T * operator * native_chart)
                    for operator in standard_hodge]
    V = sp.Matrix.hstack(*[vec(operator) for operator in native_hodge])
    require(V.rank() == 5 and
            sp.simplify(V.conjugate().T * V) == 6 * sp.eye(5),
            "native Hodge Kraus space has rank five and Gram six I")
    require(sp.simplify(sum((operator.conjugate().T * operator for operator in native_hodge),
                           sp.zeros(10))) == 3 * sp.eye(10),
            "native Hodge Kraus amplitudes normalize Phi trace preservingly")

    rays = old.source_rays()
    bell_np, _, _ = old.reflection_actions(rays)
    bell = [exact_gaussian(matrix) for matrix in bell_np]
    W1 = sp.Matrix.hstack(*[vec(matrix) for matrix in bell])
    one_joint_rank = W1.row_join(V).rank()
    one_intersection = W1.rank() + V.rank() - one_joint_rank
    require(W1.rank() == 55, "actual one-record source span has rank 55")
    require(one_intersection == 5 and one_joint_rank == W1.rank(),
            "native Hodge space is fully contained in the actual one-record span")
    frame1 = W1 * W1.conjugate().T
    frame1_eigenvalues = frame1.eigenvals()
    native_frame_residual = sp.simplify(frame1 * V - 8 * V)
    require(native_frame_residual == sp.zeros(100, 5),
            "native Hodge lies in the source frame eigenvalue eight space")

    # Exact minimum-norm coefficient matrix for W1 X1=V.  Start with the
    # unique solution on pivot columns and remove its component along ker W1.
    _, pivots = W1.rref()
    require(len(pivots) == 55, "one-record frame has 55 independent event columns")
    independent = W1[:, list(pivots)]
    pivot_solution = independent.gauss_jordan_solve(V)[0]
    X1_seed = sp.zeros(60, 5)
    for row, event in enumerate(pivots):
        X1_seed[event, :] = pivot_solution[row, :]
    nullspace1 = sp.Matrix.hstack(*W1.nullspace())
    require(nullspace1.shape == (60, 5), "one-record coefficient kernel has dimension five")
    orthogonal_projector = sp.simplify(
        sp.eye(60) - nullspace1 *
        (nullspace1.conjugate().T * nullspace1).inv() * nullspace1.conjugate().T
    )
    X1 = sp.simplify(orthogonal_projector * X1_seed)
    require(W1 * X1 == V and nullspace1.conjugate().T * X1 == sp.zeros(5),
            "exact minimum-norm one-record coefficient solution")
    coefficient_gram1 = sp.simplify(X1.conjugate().T * X1)
    gram1_eigenvalues = coefficient_gram1.eigenvals()
    require(coefficient_gram1 == sp.Rational(3, 4) * sp.eye(5),
            "one-record coefficient Gram is exactly three quarters I5")
    maximum1 = sp.Rational(3, 4)
    pmax1 = sp.simplify(1 / (20 * maximum1))
    contraction_gram1 = sp.simplify(20 * pmax1 * coefficient_gram1)
    require(pmax1 == sp.Rational(1, 15) and contraction_gram1 == sp.eye(5),
            "optimal one-record probability one fifteenth and contraction Gram I5")
    coefficient_projector = sp.simplify(sp.Rational(4, 3) * X1 * X1.conjugate().T)
    physical_effect = coefficient_projector.T
    require(coefficient_projector**2 == coefficient_projector and
            coefficient_projector.conjugate().T == coefficient_projector and
            coefficient_projector.rank() == 5,
            "minimal coefficient record projector is orthogonal of rank five")
    require(sp.simplify(W1 * physical_effect.T * W1.conjugate().T / 60 -
                       V * V.conjugate().T / 45) == sp.zeros(100),
            "physical record effect transpose convention gives exactly one fifteenth Phi")

    products_np = [second @ first for first in bell_np for second in bell_np]
    products = [exact_gaussian(matrix) for matrix in products_np]
    W2 = sp.Matrix.hstack(*[vec(matrix) for matrix in products])
    require(W2.rank() == 100, "actual two-record word span is all M10")
    # The frame operator is kept at 100 by 100 before exact inversion.
    A = np.stack([matrix.reshape(-1, order="F") for matrix in products_np], axis=1)
    S2 = exact_gaussian(A @ A.conj().T)
    require(S2.rank() == 100, "actual two-record frame is invertible")
    coefficient_gram = sp.simplify(V.conjugate().T * S2.inv() * V)
    gram_eigenvalues = coefficient_gram.eigenvals()
    require(coefficient_gram == sp.Rational(1, 64) * sp.eye(5),
            "two-record coefficient Gram is exactly one sixty fourth I5")
    maximum = sp.Rational(1, 64)
    pmax = sp.simplify(1 / (1200 * maximum))
    contraction_gram = sp.simplify(1200 * pmax * coefficient_gram)
    require(pmax == sp.Rational(4, 75) and contraction_gram == sp.eye(5),
            "optimal two-record probability four seventy fifths and contraction Gram I5")

    result = {
        "verdict": [
            "NATIVE_SPIN5_FUNCTOR_FIXES_THE_BELL10_HODGE_KRAUS_SPACE",
            "ACTUAL_ONE_RECORD_ALREADY_CONTAINS_THE_FULL_NATIVE_HODGE_KRAUS_SPACE",
            "ONE_RECORD_EXACT_NATIVE_HODGE_SUBCHANNEL_HAS_OPTIMAL_SUCCESS_1_OVER_15",
            "ACTUAL_TWO_RECORD_SPAN_CONTAINS_THE_NATIVE_HODGE_SPACE",
            "TWO_RECORD_EXACT_NATIVE_HODGE_SUBCHANNEL_HAS_OPTIMAL_SUCCESS_4_OVER_75",
        ],
        "normalization": {
            "symmetric_tensor_metric": "diag(1,2,2,2,1,2,2,1,2,1)",
            "symmetric_Pauli_HS_norm_squared": 4,
            "Bell_to_Sym2_factor": "1/2",
            "native_chart_unitary": True,
            "compact_real_frame_used": True,
        },
        "generator_test": {
            "C": "U Bell(C) = Lambda2(W(C)) U",
            "G": "U (-i Bell(G)) = Lambda2(W(G)) U",
            "exact": True,
        },
        "one_record": {
            "source_count": 60,
            "source_span_rank": W1.rank(),
            "hodge_span_rank": V.rank(),
            "joint_rank": one_joint_rank,
            "intersection_dimension": one_intersection,
            "full_hodge_support_contained": one_joint_rank == W1.rank(),
            "source_frame_eigenvalues": {str(k): v for k, v in frame1_eigenvalues.items()},
            "native_hodge_frame_eigenvalue": 8,
            "minimal_record_effect_rank": 5,
            "coefficient_record_projector_real": coefficient_projector == coefficient_projector.conjugate(),
            "physical_record_effect": "E=(4/3) conjugate(X1) X1^T; J_E=W1 E^T W1^dagger/60",
            "coefficient_gram": encode(coefficient_gram1),
            "coefficient_gram_eigenvalues": {str(k): v for k, v in gram1_eigenvalues.items()},
            "maximum_exact_hodge_success_probability": str(pmax1),
            "optimal_record_contraction_gram": encode(contraction_gram1),
            "resulting_success_map": f"({pmax1}) Phi_native",
            "input_independent_success": True,
            "record_filter_source_selected": False,
        },
        "two_record": {
            "history_count": 3600,
            "word_span_rank": W2.rank(),
            "coefficient_gram": encode(coefficient_gram),
            "coefficient_gram_eigenvalues": {str(k): v for k, v in gram_eigenvalues.items()},
            "maximum_exact_hodge_success_probability": str(pmax),
            "optimal_record_contraction_gram": encode(contraction_gram),
            "resulting_success_map": f"({pmax}) Phi_native",
            "input_independent_success": True,
            "record_filter_source_selected": False,
        },
        "scope": (
            "The native metric, compact real form and full Spin5 functor fix the target Kraus space "
            "without a duad label chart. The calculation proves exact coherent one-record and "
            "two-record support and their optimal target filters. It does not show that the source "
            "selects or implements either coherent record contraction."
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
