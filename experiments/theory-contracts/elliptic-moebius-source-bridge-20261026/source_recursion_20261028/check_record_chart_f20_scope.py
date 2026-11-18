#!/usr/bin/env python3
"""Exact scope check for the duad chart used by the two-event Hodge filter."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools as it
import json
from pathlib import Path

import numpy as np
import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = Path("/Users/stefanhamann/Projekte/tfpt-theoryv4")
SOURCE = ROOT / "experiments/theory-contracts/compiler-correlated-event-clock-20260919/source_channel.py"
HODGE = ROOT / "experiments/theory-contracts/elliptic-moebius-source-bridge-20261026/hodge_complement_20261028/hodge_complement_channel.py"
NATIVE_BRIDGE = HERE / "check_native_pauli_f20_bridge.py"
NATIVE_REPORT = HERE / "NATIVE-PAULI-F20-UNITARY-BRIDGE.txt"
RECORD_CHECKER = HERE / "source_record_hodge.py"
PINS = {
    str(SOURCE): "3852fabaf10c143fe20238b67406df6396f85c931dfd4834fe91a2d7ef900593",
    str(HODGE): "72b04a5bf8b3bd805568af49dec76bd97926fcff4009a830137e2d318807522b",
    str(NATIVE_BRIDGE): "e1200c58403e5854d91a2ad6d20d4844ff077ea2b92d4248ae8b9b1606c32149",
    str(NATIVE_REPORT): "6ddc3238986183aa3ea4b30890dc5e629be25fc00a0364905ee0479834a6ed33",
    str(RECORD_CHECKER): "5fb344b8ad49a6d85bf97228ba6d5b093576f69164fa3a3139e02eec19683119",
}

checks: list[str] = []


def require(value: bool, label: str) -> None:
    if not bool(value):
        raise AssertionError(label)
    checks.append(label)


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, "import specification " + name)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def permutation_matrix(images: tuple[int, ...] | list[int]) -> sp.Matrix:
    result = sp.zeros(len(images))
    for source, target in enumerate(images):
        result[target, source] = 1
    return result


def exterior_square(matrix: sp.Matrix) -> sp.Matrix:
    pairs = list(it.combinations(range(matrix.rows), 2))
    return sp.Matrix.hstack(*[
        sp.Matrix([
            matrix[i, a] * matrix[j, b] - matrix[j, a] * matrix[i, b]
            for i, j in pairs
        ])
        for a, b in pairs
    ])


def pauli_basis(bits) -> list[sp.Matrix]:
    identity = sp.eye(2)
    x = sp.Matrix([[0, 1], [1, 0]])
    y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    z = sp.diag(1, -1)
    one = {(0, 0): identity, (1, 0): x, (0, 1): z, (1, 1): y}
    return [sp.kronecker_product(one[(a, b)], one[(c, d)]) for a, b, c, d in bits]


def bell_congruence(matrix: sp.Matrix, scale: int, basis: list[sp.Matrix]) -> sp.Matrix:
    result = sp.zeros(10)
    for source, pauli in enumerate(basis):
        image = matrix * pauli * matrix.T
        hits = []
        for target, candidate in enumerate(basis):
            for phase in (1, sp.I, -1, -sp.I):
                if image == scale * phase * candidate:
                    hits.append((target, phase))
        require(len(hits) == 1, f"unique Bell congruence image {source}")
        target, phase = hits[0]
        result[target, source] = phase
    return result


def to_sympy_integer(array: np.ndarray) -> sp.Matrix:
    require(np.array_equal(array.real, np.rint(array.real)), "Gaussian real coordinates")
    require(np.array_equal(array.imag, np.rint(array.imag)), "Gaussian imaginary coordinates")
    return sp.Matrix([
        [sp.Integer(int(round(z.real))) + sp.I * sp.Integer(int(round(z.imag))) for z in row]
        for row in array
    ])


def vector_columns(matrices: list[np.ndarray]) -> np.ndarray:
    return np.stack([matrix.reshape(-1, order="F") for matrix in matrices], axis=1)


def coefficient_gram(bell: list[np.ndarray], chart: sp.Matrix, V: sp.Matrix) -> sp.Matrix:
    chart_np = np.array(chart.tolist(), dtype=complex)
    unitaries = [chart_np.conj().T @ U @ chart_np for U in bell]
    products = [second @ first for first in unitaries for second in unitaries]
    A = vector_columns(products)
    S = to_sympy_integer(A @ A.conj().T)
    require(S.rank() == 100, "two-event frame has full rank")
    gram = sp.simplify(V.conjugate().T * S.inv() * V)
    return gram


def main(output: Path | None) -> dict:
    checks.clear()
    for path, digest in PINS.items():
        require(hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest, "source pin " + path)

    old = load("pinned_source", SOURCE)
    hodge = load("pinned_hodge", HODGE)
    basis = pauli_basis(old.SYMMETRIC_BITS)

    # Twice the normalized C matrix keeps the congruence test integral.
    c_twice = sp.Matrix([
        [-sp.I, -sp.I, -sp.I, -sp.I],
        [-sp.I, sp.I, -sp.I, sp.I],
        [1, 1, -1, -1],
        [-1, 1, 1, -1],
    ])
    g_native = sp.Matrix([
        [0, 0, 0, sp.I],
        [0, 0, -1, 0],
        [sp.I, 0, 0, 0],
        [0, 1, 0, 0],
    ])
    bell_c = bell_congruence(c_twice, 4, basis)
    bell_g = bell_congruence(g_native, 1, basis)
    require(bell_c**5 == sp.eye(10) and bell_g**4 == sp.eye(10),
            "native Bell F20 generator orders")

    p5 = permutation_matrix([(j + 1) % 5 for j in range(5)])
    pg = permutation_matrix([(2 * j) % 5 for j in range(5)])
    hodge_c = exterior_square(p5)
    hodge_g = exterior_square(pg)

    rays = old.source_rays()
    bell, raw_six, _ = old.reflection_actions(rays)
    six, _ = old.outer_mark_action(raw_six)
    mapping = old.find_bell_partition_intertwiner(bell, six)

    def infer_six(action: sp.Matrix) -> list[tuple[int, ...]]:
        bell_perm = tuple(next(row for row in range(10) if action[row, col] != 0)
                          for col in range(10))
        return [perm for perm in it.permutations(range(6))
                if all(mapping[bell_perm[b]] == old.act_partition(mapping[b], perm)
                       for b in range(10))]

    six_c, six_g = infer_six(bell_c), infer_six(bell_g)
    require(six_c == [(4, 3, 0, 2, 1, 5)], "native C outer-six action")
    require(six_g == [(2, 3, 1, 0, 4, 5)], "native G outer-six action")
    fixed = [mark for mark in range(6) if six_c[0][mark] == mark and six_g[0][mark] == mark]
    require(fixed == [5], "native F20 fixes exactly outer mark five")

    chart_np, _ = hodge.duad_chart(5, mapping)
    raw_chart = sp.Matrix(chart_np.tolist())
    raw_c = raw_chart.T * bell_c * raw_chart
    raw_g = raw_chart.T * bell_g * raw_chart
    require(raw_c != hodge_c and sp.I * raw_g != hodge_g and -sp.I * raw_g != hodge_g,
            "raw mark-five duad chart is not an F20 intertwiner for either quarter orientation")

    # Exhaust the five-slot relabelings. Exactly one aligns both underlying
    # permutations, but its phases still fail the generator equations.
    relabel_hits = []
    for perm in it.permutations(range(5)):
        q = exterior_square(permutation_matrix(perm))
        c_aligned = q.T * raw_c * q
        g_aligned = q.T * raw_g * q
        c_support = c_aligned.applyfunc(lambda z: int(z != 0))
        g_support = g_aligned.applyfunc(lambda z: int(z != 0))
        if c_support == hodge_c.applyfunc(lambda z: int(z != 0)) and \
                g_support == hodge_g.applyfunc(lambda z: int(z != 0)):
            relabel_hits.append(perm)
    require(relabel_hits == [(4, 1, 3, 2, 0)], "unique slot relabeling aligns F20 permutations")
    q = exterior_square(permutation_matrix(relabel_hits[0]))
    aligned_c = q.T * raw_c * q
    aligned_g = q.T * raw_g * q
    require(aligned_c != hodge_c and sp.I * aligned_g != hodge_g and -sp.I * aligned_g != hodge_g,
            "slot relabeling alone does not repair Bell phases")

    d_plus = sp.diag(1, -sp.I, sp.I, 1, -sp.I, -1, -sp.I, -sp.I, sp.I, 1)
    d_minus = sp.diag(1, sp.I, -sp.I, 1, -sp.I, 1, sp.I, -sp.I, -sp.I, 1)
    require(d_plus * aligned_c * d_plus.inv() == hodge_c and
            d_plus * (sp.I * aligned_g) * d_plus.inv() == hodge_g,
            "positive quarter diagonal phase repair is an exact F20 intertwiner")
    require(d_minus * aligned_c * d_minus.inv() == hodge_c and
            d_minus * (-sp.I * aligned_g) * d_minus.inv() == hodge_g,
            "negative quarter diagonal phase repair is an exact F20 intertwiner")

    bridge_plus = raw_chart * q * d_plus.inv()
    bridge_minus = raw_chart * q * d_minus.inv()
    require(bridge_plus.conjugate().T * bridge_plus == sp.eye(10) and
            bridge_minus.conjugate().T * bridge_minus == sp.eye(10),
            "repaired Bell-to-wedge charts are unitary")
    require(bridge_plus.conjugate().T * bell_c * bridge_plus == hodge_c and
            bridge_plus.conjugate().T * (sp.I * bell_g) * bridge_plus == hodge_g,
            "positive repaired chart satisfies both generator equalities")

    Cs = hodge.hodge_amplitudes()
    V = to_sympy_integer(vector_columns(Cs).astype(complex))
    raw_gram = coefficient_gram(bell, raw_chart, V)
    plus_gram = coefficient_gram(bell, bridge_plus, V)
    minus_gram = coefficient_gram(bell, bridge_minus, V)
    expected_raw = sp.diag(sp.Rational(1, 64), sp.Rational(1, 48), sp.Rational(1, 48),
                           sp.Rational(1, 48), sp.Rational(1, 64))
    expected_covariant = sp.eye(5) / 48
    require(raw_gram == expected_raw, "raw label chart has the previous three-plus-two Gram spectrum")
    require(plus_gram == expected_covariant and minus_gram == expected_covariant,
            "both repaired F20 charts have scalar Hodge coefficient Gram")
    require(raw_gram * p5 != p5 * raw_gram,
            "raw three-plus-two eigenspaces are not covariant under the native five-cycle")
    require(plus_gram * p5 == p5 * plus_gram and plus_gram * pg == pg * plus_gram,
            "repaired coefficient Gram is F20 invariant")

    result = {
        "verdict": [
            "RAW_DUAD_CHART_IS_NOT_AN_F20_INTERTWINER",
            "ONE_SLOT_RELABELING_PLUS_A_DIAGONAL_MU4_GAUGE_REPAIRS_BOTH_F20_GENERATORS",
            "PMAX_REMAINS_1_OVER_25_FOR_THE_TWO_EXPLICIT_REPAIRED_UNITARY_CHARTS",
            "RAW_THREE_PLUS_TWO_FILTER_SPLITTING_IS_NOT_F20_COVARIANT_AND_DISAPPEARS_AFTER_REPAIR",
            "FULL_U2_TIMES_U1_TIMES_U1_TEN_INTERTWINER_GAUGE_IS_NOT_CLASSIFIED_BY_THIS_TEST",
        ],
        "native_outer_fixed_mark": 5,
        "raw_chart": {
            "f20_intertwiner": False,
            "two_event_coefficient_gram_diagonal": ["1/64", "1/48", "1/48", "1/48", "1/64"],
            "pmax": "1/25",
            "three_plus_two_invariant_under_C5": False,
        },
        "repair": {
            "slot_relabeling": list(relabel_hits[0]),
            "positive_quarter_diagonal": [str(d_plus[j, j]) for j in range(10)],
            "negative_quarter_diagonal": [str(d_minus[j, j]) for j in range(10)],
            "generator_equalities_exact": True,
            "unitary": True,
            "two_event_coefficient_gram": "I5/48 for both quarter orientations",
            "pmax": "1/25",
            "minimal_record_contraction_row_gram_at_pmax": "I5",
        },
        "scope": (
            "The value 1/25 survives the two explicit monomial unitary F20 repairs, but the raw chart "
            "used previously was not itself an F20 intertwiner. The native ten has residual U(2) times "
            "U(1) times U(1) intertwiner freedom; invariance over that full family is not proved."
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
    result = main(args.out)
    print(json.dumps(result, indent=2, sort_keys=True))
