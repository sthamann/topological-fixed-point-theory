#!/usr/bin/env python3
"""Exact lift of the six native H5 effects to the original erasure H20.

This checker reconstructs the decoder isometry from the original quartic word
supports.  It then lifts the doily effects through that concrete isometry and
checks the actual 4 x 5 output action, all sixty native reflection words, and
the 120-word S5 stabilizer average for every one of the six rays.
"""

from __future__ import annotations

import argparse
from collections import deque
import hashlib
import importlib.util
import itertools as it
import json
from pathlib import Path
import sys
from typing import Any

import sympy as sp


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DEFAULT_OUTPUT = HERE / "native_erasure_readout_results.json"
DOILY_SOURCE = HERE / "doily_quartic_bridge.py"
PROCESS_SOURCE = REPO / "tfpt_explorer/process.py"
ORIGINAL_ERASURE_SOURCE = (
    REPO
    / "experiments/theory-contracts/joint-native-algebra-selection-20261016/operator_checker.py"
)

if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from tfpt_explorer.process import _code_basis  # noqa: E402


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def clean(matrix: sp.MatrixBase) -> sp.Matrix:
    return sp.Matrix(matrix).applyfunc(sp.simplify)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def expr_json(value: sp.Expr) -> int | str:
    value = sp.simplify(value)
    if value.is_Integer:
        return int(value)
    return str(value)


def matrix_json(matrix: sp.MatrixBase) -> list[list[int | str]]:
    return [
        [expr_json(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def matrix_sha256(matrix: sp.MatrixBase) -> str:
    payload = json.dumps(matrix_json(matrix), separators=(",", ":"), sort_keys=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def sparse_matrix_json(matrix: sp.MatrixBase) -> dict[str, Any]:
    entries = []
    for row in range(matrix.rows):
        for column in range(matrix.cols):
            if matrix[row, column] != 0:
                entries.append([row, column, expr_json(matrix[row, column])])
    return {"shape": [matrix.rows, matrix.cols], "entries": entries}


def load_doily_module() -> Any:
    spec = importlib.util.spec_from_file_location("doily_quartic_bridge_native_readout", DOILY_SOURCE)
    require(spec is not None and spec.loader is not None, "cannot load doily checker")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def decoder_isometry_gauge_first(
    supports: list[list[tuple[int, ...]]], erased: int = 0
) -> sp.Matrix:
    """Return W:C4_gauge tensor H5 -> the physical three-register C64.

    The original checker uses logical-first columns 4*i+a.  The requested
    decoder order is gauge first, so column 5*a+i contains the same vector.
    """

    out = sp.zeros(64, 20)
    for logical, support in enumerate(supports):
        amplitude = sp.simplify(2 / sp.sqrt(len(support)))
        for word in support:
            gauge = word[erased]
            remaining = tuple(word[k] for k in range(4) if k != erased)
            row = remaining[0] * 16 + remaining[1] * 4 + remaining[2]
            out[row, 5 * gauge + logical] += amplitude
    return out


def reflection_from_gaussian_ray(ray: list[list[int]]) -> sp.Matrix:
    root = sp.Matrix([sp.Integer(a) + sp.I * sp.Integer(b) for a, b in ray])
    return clean(sp.eye(4) - root * root.H / 2)


def tensor_power(matrix: sp.Matrix, count: int) -> sp.Matrix:
    out = sp.ones(1, 1)
    for _ in range(count):
        out = sp.kronecker_product(out, matrix)
    return sp.Matrix(out)


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    """Permutation for the matrix product left_matrix * right_matrix."""

    return tuple(left[right[index]] for index in range(len(left)))


def group_closure(
    generators: list[tuple[int, tuple[int, ...], sp.Matrix]],
) -> dict[tuple[int, ...], tuple[sp.Matrix, list[int]]]:
    identity = tuple(range(6))
    elements: dict[tuple[int, ...], tuple[sp.Matrix, list[int]]] = {
        identity: (sp.eye(5), [])
    }
    queue: deque[tuple[int, ...]] = deque([identity])
    while queue:
        permutation = queue.popleft()
        matrix, word = elements[permutation]
        for reflection_index, generator_permutation, generator_matrix in generators:
            next_permutation = compose(permutation, generator_permutation)
            next_matrix = clean(matrix * generator_matrix)
            next_word = word + [reflection_index]
            if next_permutation in elements:
                require(
                    elements[next_permutation][0] == next_matrix,
                    "the same S6 permutation acquired two quartic matrices",
                )
                continue
            elements[next_permutation] = (next_matrix, next_word)
            queue.append(next_permutation)
    return elements


def main(output: Path) -> dict[str, Any]:
    doily = load_doily_module()
    bridge = doily.build_report()
    rows = bridge["quartic_bridge"]["per_reflection"]
    require(len(rows) == 60, "native reflection count changed")

    _, supports = _code_basis()
    support_sizes = [len(support) for support in supports]
    require(support_sizes == [4, 12, 12, 12, 24], "quartic supports changed")

    # This is the actual original erasure carrier, with only the tensor-factor
    # order changed from H5 tensor C4 to C4 tensor H5.
    W = decoder_isometry_gauge_first(supports)
    U_dec = W.H
    require(clean(W.H * W) == sp.eye(20), "decoder carrier is not isometric")
    carrier_projector = clean(W * W.H)
    require(carrier_projector.rank() == 20, "physical erasure carrier rank changed")
    require(clean(carrier_projector**2) == carrier_projector, "carrier projector is not idempotent")

    metric = sp.diag(1, 3, 3, 3, 6)
    D = sp.diag(1, sp.sqrt(3), sp.sqrt(3), sp.sqrt(3), sp.sqrt(6))
    rays = [sp.Matrix(ray) for ray in bridge["selection_data"]["native_ovoid_rays"]]
    require(len(rays) == 6, "ovoid ray count changed")

    effects_native = [clean(ray * ray.T * metric / 12) for ray in rays]
    effects5 = [clean(D * effect * D.inv()) for effect in effects_native]
    projectors_native = [clean(sp.Rational(6, 5) * effect) for effect in effects_native]
    projectors5 = [clean(sp.Rational(6, 5) * effect) for effect in effects5]
    require(clean(sum(effects5, sp.zeros(5))) == sp.eye(5), "H5 effects do not sum to identity")
    require(all(effect.rank() == 1 for effect in effects5), "an H5 effect lost rank one")
    require(
        all(clean(effect**2 - sp.Rational(5, 6) * effect) == sp.zeros(5) for effect in effects5),
        "H5 effect normalization changed",
    )
    require(
        all(clean(projector**2 - projector) == sp.zeros(5) for projector in projectors5),
        "scaled H5 effects are not projectors",
    )

    effects20 = [sp.kronecker_product(sp.eye(4), effect) for effect in effects5]
    physical_effects = [clean(W * effect * W.H) for effect in effects20]
    require(clean(sum(effects20, sp.zeros(20))) == sp.eye(20), "H20 effects do not sum to identity")
    require(
        clean(sum(physical_effects, sp.zeros(64))) == carrier_projector,
        "physical effects do not sum to the identity on H20",
    )
    require(all(effect.rank() == 4 for effect in effects20), "an H20 effect lost rank four")
    require(all(effect.rank() == 4 for effect in physical_effects), "a physical effect lost rank four")
    require(
        all(clean(U_dec * physical * W - intrinsic) == sp.zeros(20) for physical, intrinsic in zip(physical_effects, effects20)),
        "physical effects do not decode to I4 tensor E5",
    )

    # Verify the probability formula for a completely symbolic 5 by 5 matrix.
    rho_symbols = sp.symbols("rho0:25")
    rho5 = sp.Matrix(5, 5, rho_symbols)
    rho20 = sp.kronecker_product(sp.eye(4) / 4, rho5)
    probability_identities = []
    for effect5, effect20 in zip(effects5, effects20):
        lhs = sp.expand(sp.trace(rho20 * effect20))
        rhs = sp.expand(sp.trace(rho5 * effect5))
        probability_identities.append(sp.simplify(lhs - rhs) == 0)
    require(all(probability_identities), "erasure-state probability identity failed")

    # Recover every actual native 4/5/20 representation matrix and test it on
    # the concrete W.  The 5 matrices in the doily report use unnormalised word
    # sums, hence the D C D^-1 conversion to the normalized original code basis.
    reflection_checks = []
    unique_transpositions: dict[tuple[int, ...], tuple[int, sp.Matrix]] = {}
    for row in rows:
        permutation = tuple(row["ovoid_permutation"])
        C_native = sp.Matrix(row["quartic_numerator"]) / 16
        C5 = clean(D * C_native * D.inv())
        R4 = reflection_from_gaussian_ray(row["gaussian_ray"])
        G20 = sp.kronecker_product(R4.conjugate(), C5)
        R3 = tensor_power(R4, 3)
        decoder_ok = clean(R3 * W - W * G20) == sp.zeros(64, 20)
        covariance_ok = all(
            clean(G20 * effects20[index] * G20.H - effects20[permutation[index]])
            == sp.zeros(20)
            for index in range(6)
        )
        require(decoder_ok and covariance_ok, f"reflection {row['index']} failed H20 covariance")
        reflection_checks.append(
            {
                "reflection_index": row["index"],
                "ovoid_permutation": list(permutation),
                "decoder_intertwiner_exact": decoder_ok,
                "all_six_effects_covariant": covariance_ok,
                "gauge_first_H20_matrix_sha256": matrix_sha256(G20),
            }
        )
        if permutation not in unique_transpositions:
            unique_transpositions[permutation] = (row["index"], C_native)
        else:
            require(
                unique_transpositions[permutation][1] == C_native,
                "duplicate transposition has inconsistent quartic action",
            )
    require(len(unique_transpositions) == 15, "expected fifteen S6 transpositions")

    generators = [
        (reflection_index, permutation, matrix)
        for permutation, (reflection_index, matrix) in sorted(unique_transpositions.items())
    ]
    group = group_closure(generators)
    require(len(group) == 720, "native transpositions did not generate S6")

    stabilizer_records = []
    for ray_index in range(6):
        stabilizer = [
            (permutation, matrix, word)
            for permutation, (matrix, word) in group.items()
            if permutation[ray_index] == ray_index
        ]
        require(len(stabilizer) == 120, "ovoid stabilizer is not S5")
        average = clean(sum((matrix for _, matrix, _ in stabilizer), sp.zeros(5)) / 120)
        require(average == projectors_native[ray_index], "S5 average is not (6/5) E_r")
        require(
            clean(D * average * D.inv()) == projectors5[ray_index],
            "orthonormal stabilizer average changed",
        )
        words = sorted((word for _, _, word in stabilizer), key=lambda word: (len(word), word))
        stabilizer_records.append(
            {
                "ray_index": ray_index,
                "order": len(stabilizer),
                "average_equals_six_fifths_effect": True,
                "max_word_length_in_bfs_certificate": max(map(len, words)),
                "words_in_native_reflection_indices": words,
            }
        )

    report = {
        "claim": "the six doily effects are exact observables on the original erasure H20",
        "status": "PASS",
        "arithmetic": "exact integers, rationals, Gaussian integers, and algebraic square roots",
        "sources": {
            str(PROCESS_SOURCE.relative_to(REPO)): sha256(PROCESS_SOURCE),
            str(ORIGINAL_ERASURE_SOURCE.relative_to(REPO)): sha256(ORIGINAL_ERASURE_SOURCE),
            str(DOILY_SOURCE.relative_to(REPO)): sha256(DOILY_SOURCE),
        },
        "basis_conventions": {
            "physical_C64_row": "16*x1 + 4*x2 + x3 after erasing tensor slot x0",
            "decoded_H20_column": "5*gauge + logical, i.e. C4_gauge tensor H5",
            "original_checker_column": "4*logical + gauge; related by an explicit permutation only",
            "native_H5_metric": matrix_json(metric),
            "orthonormal_conversion": "D=diag(1,sqrt(3),sqrt(3),sqrt(3),sqrt(6)); C_orth=D C_native D^-1",
        },
        "decoder": {
            "W_C4_tensor_H5_to_C64": sparse_matrix_json(W),
            "U_dec_C64_to_C4_tensor_H5": sparse_matrix_json(U_dec),
            "W_shape": [64, 20],
            "U_dec_shape": [20, 64],
            "W_star_W_equals_I20": True,
            "physical_carrier_projector_rank": carrier_projector.rank(),
            "physical_carrier_projector_sha256": matrix_sha256(carrier_projector),
        },
        "effects": {
            "H5_orthonormal_matrices": [matrix_json(effect) for effect in effects5],
            "H20_gauge_first_matrices": [matrix_json(effect) for effect in effects20],
            "physical_C64_sparse_matrices": [sparse_matrix_json(effect) for effect in physical_effects],
            "H5_ranks": [effect.rank() for effect in effects5],
            "H20_ranks": [effect.rank() for effect in effects20],
            "sum_H5": "I5",
            "sum_H20": "I20",
            "sum_physical_C64": "W W* = identity on the original H20 carrier",
            "probability_identity_for_symbolic_rho5": all(probability_identities),
            "probability_formula": "Tr[W(I4/4 tensor rho5)W* · W(I4 tensor E5_r)W*] = Tr[rho5 E5_r]",
        },
        "native_action": {
            "formula": "R4^tensor3 W = W(conjugate(R4) tensor C5) in gauge-first order",
            "all_60_decoder_intertwiners_exact": all(x["decoder_intertwiner_exact"] for x in reflection_checks),
            "all_60_times_6_effect_covariances_exact": all(x["all_six_effects_covariant"] for x in reflection_checks),
            "per_reflection": reflection_checks,
        },
        "stabilizer_average": {
            "generated_group_order": len(group),
            "unique_native_transposition_generators": [
                {"reflection_index": index, "ovoid_permutation": list(permutation)}
                for index, permutation, _ in generators
            ],
            "formula": "P_r=(6/5)E_r=(1/120) sum_{g in Stab_S6(r)} C5(g)",
            "all_six_stabilizers": stabilizer_records,
        },
        "interpretation": {
            "new_dynamics_added": False,
            "auxiliary_hilbert_space_added": False,
            "original_erasure_instrument_identified_with_doily_povm": False,
            "result": "an additional covariant logical readout on the already existing physical erasure carrier",
        },
    }
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = main(args.output)
    print(
        "PASS: W=64x20 exact; six H20 effects rank 4 and sum I20; "
        "60 native actions covariant; six S5 stabilizers have 120 words."
    )
