#!/usr/bin/env python3
"""Exact support test for Hodge heralding from the original event records.

Only the source Kraus matrices and their recorded coherent span are used.  A
record effect can change the positive coefficient matrix on that span, but it
cannot add a Kraus direction outside it.  The checker therefore tests the
necessary support inclusion before any optimization over record effects.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import sympy as sp


HERE = Path(__file__).resolve().parent
REPO = Path("/Users/stefanhamann/Projekte/tfpt-theoryv4")
SOURCE_DIR = REPO / "experiments/theory-contracts/compiler-correlated-event-clock-20260919"
OLD_DIR = REPO / "experiments/theory-contracts/primitive-transfer-selection-20260912"
SOURCE = SOURCE_DIR / "source_channel.py"
SOURCE_PROOF = SOURCE_DIR / "PROOF.txt"
HODGE = REPO / "experiments/theory-contracts/elliptic-moebius-source-bridge-20261026/hodge_complement_20261028/hodge_complement_channel.py"
PLUCKER = REPO / "experiments/tfpt-discovery/s6_plucker_hadamard_probe.py"
V888 = REPO / "verification/v888_finite_pfaffian_signature.py"

PINS = {
    str(SOURCE): "3852fabaf10c143fe20238b67406df6396f85c931dfd4834fe91a2d7ef900593",
    str(SOURCE_PROOF): "bf1f85c99f5f89a0a15da1af5e3ff321c2756fda132a9f8d644f9ce0e08c2f5c",
    str(HODGE): "72b04a5bf8b3bd805568af49dec76bd97926fcff4009a830137e2d318807522b",
    str(PLUCKER): "d65e7dae7cfe8dafac2153a75b350e7c23a6640d59b40f0fa6ad951d8f0caf13",
    str(V888): "f29ae2e0d7c6e679987b2b71413a4ee9c659cbcd74d9b7e3e6297a71ff0b280c",
    str(OLD_DIR / "coherent_compiler_record.py"): "db28d0111b8455439e1abb43c961ac67f5fc6572c9efe76449f17c11b71b6b55",
    str(OLD_DIR / "coherent_record_selection.py"): "6f3ee9f4b815a397d607355c0e07872c47a834b25e654db3e73bb2ac0419120c",
    str(OLD_DIR / "primitive_record_recovery.py"): "f0c2919798ef7c79f71f2f5fa454bd2226c220fd6c2c699f2e2cdb7f626f6cbb",
    str(OLD_DIR / "COHERENT_COMPILER_RECORD.md"): "d38fa6a5fa163f4deac47cfa1e377a57139b825353b84c3eb547d0b22db5bbe8",
    str(OLD_DIR / "COHERENT_RECORD_SELECTION.md"): "41fdf493a664ee3c56119a818c032edc2d24a7b833219af0f057e404cdfdb6a9",
    str(OLD_DIR / "PRIMITIVE_RECORD_RECOVERY.md"): "1af5a9c94358a6dd30c260a01dd8422a8ff47f9ce75104ecb8a10a0b04705154",
    str(OLD_DIR / "ASSEMBLY_NOT_HISTORY.md"): "c4897a92da4aa73effff021330d436324d3d7c9a3d959b81e7ee204bac164385",
}

checks: list[str] = []


def require(ok: bool, name: str) -> None:
    if not bool(ok):
        raise RuntimeError(name)
    checks.append(name)


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"import specification {name}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def exact_gaussian(a: np.ndarray) -> sp.Matrix:
    require(np.allclose(a.real, np.rint(a.real), atol=0, rtol=0), "Gaussian integer real coordinates")
    require(np.allclose(a.imag, np.rint(a.imag), atol=0, rtol=0), "Gaussian integer imaginary coordinates")
    return sp.Matrix([
        [sp.Integer(int(round(z.real))) + sp.I * sp.Integer(int(round(z.imag))) for z in row]
        for row in a
    ])


def vector_columns(matrices: list[np.ndarray]) -> np.ndarray:
    return np.stack([a.reshape(-1, order="F") for a in matrices], axis=1)


def intersection_data(W: sp.Matrix, V: sp.Matrix) -> tuple[int, list[int]]:
    """Return intersection dimension and its coordinates in the Hodge basis."""
    joint = W.row_join(V)
    dimension = W.rank() + V.rank() - joint.rank()
    relations = W.row_join(-V).nullspace()
    hodge_parts = [r[W.cols:, :] for r in relations]
    B = sp.Matrix.hstack(*hodge_parts) if hodge_parts else sp.zeros(V.cols, 0)
    require(B.rank() == dimension, "relation nullspace recovers intersection dimension")
    support = [i for i in range(V.cols) if any(B[i, j] != 0 for j in range(B.cols))]
    require(len(support) == dimension, "Hodge intersection is a coordinate plane")
    require(all(B[i, j] == 0 for i in range(V.cols) if i not in support for j in range(B.cols)),
            "no hidden Hodge direction outside reported plane")
    return dimension, support


def no_scalar_success(Cs: list[np.ndarray], support: list[int]) -> None:
    """Even arbitrary coefficients on the accessible plane cannot sum to p I except at zero."""
    columns = []
    for m in support:
        for n in support:
            columns.append((Cs[n].conj().T @ Cs[m]).reshape(-1, order="F"))
    columns.append(-np.eye(10).reshape(-1, order="F"))
    system = exact_gaussian(np.stack(columns, axis=1))
    require(system.rank() == len(columns), "accessible Hodge plane has no nonzero scalar-effect solution")


def plucker_typing(Cs: list[np.ndarray]) -> dict:
    """Separate the wedge map from multiplication in the Pluecker quotient."""
    symmetric_basis = [(a, b) for a in range(10) for b in range(a, 10)]
    require(len(symmetric_basis) == 55, "Sym2 of the ten dimensional carrier has dimension 55")
    quadrics = sp.Matrix([
        [int(C[a, a]) if a == b else int(C[a, b] + C[b, a]) for a, b in symmetric_basis]
        for C in Cs
    ])
    require(quadrics.rank() == 5, "five independent epsilon wedge components")
    require(all(sum(int(x != 0) for x in quadrics.row(m)) == 3 for m in range(5)),
            "each epsilon component is one three-term Pluecker quadric")
    require(all(abs(int(x)) in (0, 2) for x in quadrics),
            "symmetric-coordinate Pluecker coefficients have the expected common factor two")
    require(len(quadrics.nullspace()) == 50, "wedge map has kernel dimension 50")
    return {
        "domain": "Sym^2(Lambda^2 C5)",
        "domain_dimension": 55,
        "wedge_map_image_dimension": 5,
        "wedge_map_kernel_dimension": 50,
        "holomorphic_coordinate_multiplication_kernel_dimension": 5,
        "holomorphic_coordinate_multiplication_image_dimension": 50,
        "typing": "the same five tensors are wedge components and Pluecker relations for complementary maps",
    }


def marked_jump_kraus(old, mark: int, bell, six, mapping, chart: np.ndarray) -> list[np.ndarray]:
    indices, *_ = old.marked_data(mark, bell, six, mapping)
    out = []
    for k in indices:
        U = bell[k]
        permutation = old.bell_permutation(U)
        jump = np.eye(10) - np.diag([int(permutation[a] == a) for a in range(10)])
        out.append(chart.T @ (U @ jump) @ chart)
    require(len(out) == 20, f"mark {mark}: actual twenty jump Kraus matrices")
    return out


def two_event_record(Cs: list[np.ndarray], V: sp.Matrix, bell, chart: np.ndarray,
                     one_step_plane: list[int], mark: int) -> dict:
    """Exact optimal Hodge subchannel inside two actual sixty-event steps."""
    unitaries = [chart.T @ U @ chart for U in bell]
    products = [second @ first for first in unitaries for second in unitaries]
    A = vector_columns(products)
    # S=A A^dagger is the unnormalised two-step Choi frame.  Computing S
    # first keeps the exact calculation at 100 by 100 rather than 100 by 3600.
    S = exact_gaussian(A @ A.conj().T)
    require(S.rank() == 100, f"chart {mark}: two actual events span all M10")
    Y = S.inv() * V
    require(S * Y == V, f"chart {mark}: exact minimum-norm inverse equation")
    gram = sp.simplify(V.conjugate().T * Y)
    expected = sp.diag(*[sp.Rational(1, 64) if m in one_step_plane else sp.Rational(1, 48)
                         for m in range(5)])
    require(gram == expected, f"chart {mark}: exact minimum-norm record coefficient Gram matrix")

    # C60 composed twice has source amplitude 1/60 for every word.  At success
    # probability p the record contraction has row Gram 1200*p*gram.  Hence
    # p=1/25 gives 48*gram <= I, and any larger p violates the 1/48 directions.
    contraction_gram = 48 * gram
    require(all(contraction_gram[m, n] == 0 for m in range(5) for n in range(5) if m != n),
            f"chart {mark}: optimal record contraction rows are orthogonal")
    require([contraction_gram[m, m] for m in range(5)] ==
            [sp.Rational(3, 4) if m in one_step_plane else sp.Integer(1) for m in range(5)],
            f"chart {mark}: optimal record contraction is positive and contractive")
    require(sp.eye(5) - contraction_gram ==
            sp.diag(*[sp.Rational(1, 4) if m in one_step_plane else sp.Integer(0)
                      for m in range(5)]),
            f"chart {mark}: exact complement of optimal record contraction")
    return {
        "duad_chart_mark": mark,
        "actual_history_count": 3600,
        "word_kraus_span_rank": 100,
        "minimum_norm_coefficient_gram_diagonal": [str(gram[m, m]) for m in range(5)],
        "one_step_intersection_indices": one_step_plane,
        "optimal_record_contraction_row_gram_diagonal":
            [str(contraction_gram[m, m]) for m in range(5)],
        "maximum_exact_hodge_success_probability": "1/25",
        "resulting_success_map": "(1/25) Phi",
        "input_independent_success": True,
        "filter_formula": (
            "A has columns vec(U_second U_first); X=A^dagger(AA^dagger)^-1 V; "
            "R_mw=60 sqrt((1/25)/3) X_wm"
        ),
        "filter_is_source_selected": False,
    }


def main(output: Path | None) -> dict:
    checks.clear()
    for path, digest in PINS.items():
        require(hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest, "source pin " + path)

    old = load_module("pinned_source_channel", SOURCE)
    hodge = load_module("pinned_hodge_channel", HODGE)
    rays = old.source_rays()
    bell, raw_six, _ = old.reflection_actions(rays)
    six, _ = old.outer_mark_action(raw_six)
    mapping = old.find_bell_partition_intertwiner(bell, six)
    Cs = hodge.hodge_amplitudes()
    V = exact_gaussian(vector_columns(Cs).astype(complex))
    require(V.rank() == 5, "five independent Hodge Kraus directions")
    product_typing = plucker_typing(Cs)

    marked = []
    whole = []
    two_event = []
    for mark in range(6):
        chart, _ = hodge.duad_chart(mark, mapping)

        W20 = exact_gaussian(vector_columns(marked_jump_kraus(old, mark, bell, six, mapping, chart)))
        rank20 = W20.rank()
        joint20 = W20.row_join(V).rank()
        intersection20, plane20 = intersection_data(W20, V)
        require(rank20 == 15, f"mark {mark}: exact marked jump Kraus rank 15")
        require(joint20 == 18 and intersection20 == 2,
                f"mark {mark}: marked source intersects Hodge span in dimension two")
        require(not W20.row_join(V).rank() == W20.rank(),
                f"mark {mark}: full Hodge support is not contained in marked source")
        no_scalar_success(Cs, plane20)
        marked.append({
            "mark": mark,
            "source_kraus_count": 20,
            "source_span_rank": rank20,
            "hodge_span_rank": 5,
            "joint_span_rank": joint20,
            "intersection_dimension": intersection20,
            "intersection_hodge_indices": plane20,
            "full_hodge_support_contained": False,
            "largest_exact_hodge_cp_subchannel_probability": "0",
            "input_independent_success_from_intersection_plane": False,
        })

        # This is the one actual sixty-event source, represented in one of the
        # six alternative source-derived duad charts.  The six charts are not
        # stacked into a new 360-event process.
        W60 = exact_gaussian(vector_columns([chart.T @ U @ chart for U in bell]))
        rank60 = W60.rank()
        joint60 = W60.row_join(V).rank()
        intersection60, plane60 = intersection_data(W60, V)
        require(rank60 == 55, f"chart {mark}: actual sixty-event Kraus rank 55")
        require(joint60 == 58 and intersection60 == 2,
                f"chart {mark}: whole source intersects Hodge span in dimension two")
        require(plane60 == plane20, f"chart {mark}: marked and whole-source Hodge planes agree")
        require(not W60.row_join(V).rank() == W60.rank(),
                f"chart {mark}: full Hodge support is not contained in whole source")
        no_scalar_success(Cs, plane60)
        whole.append({
            "duad_chart_mark": mark,
            "same_actual_source_kraus_count": 60,
            "source_span_rank": rank60,
            "hodge_span_rank": 5,
            "joint_span_rank": joint60,
            "intersection_dimension": intersection60,
            "intersection_hodge_indices": plane60,
            "full_hodge_support_contained": False,
            "largest_exact_hodge_cp_subchannel_probability": "0",
            "input_independent_success_from_intersection_plane": False,
        })
        two_event.append(two_event_record(Cs, V, bell, chart, plane60, mark))

    require([x["intersection_hodge_indices"] for x in marked] ==
            [[3, 4], [1, 2], [1, 2], [1, 2], [0, 4], [0, 4]],
            "six exact marked intersection planes")
    require([x["intersection_hodge_indices"] for x in whole] ==
            [[3, 4], [1, 2], [1, 2], [1, 2], [0, 4], [0, 4]],
            "six exact whole-source intersection planes")

    result = {
        "verdict": [
            "NO_SINGLE_MARK_COHERENT_RECORD_EFFECT_REALIZES_A_NONZERO_EXACT_HODGE_CHANNEL",
            "NO_ONE_EVENT_WHOLE60_RECORD_EFFECT_REALIZES_A_NONZERO_EXACT_HODGE_CHANNEL",
            "EACH_SOURCE_CHART_SHARES_EXACTLY_TWO_OF_THE_FIVE_HODGE_KRAUS_DIRECTIONS",
            "THE_ACCESSIBLE_TWO_PLANE_CANNOT_HAVE_NONZERO_INPUT_INDEPENDENT_SUCCESS",
            "TWO_ACTUAL_WHOLE60_EVENTS_DO_SUPPORT_AN_EXACT_HODGE_SUBCHANNEL_WITH_OPTIMAL_SUCCESS_1_OVER_25",
            "THE_REQUIRED_COHERENT_RECORD_FILTER_IS_TARGET_CONSTRUCTED_AND_NOT_SOURCE_SELECTED",
        ],
        "one_event_mathematical_reason": (
            "For a source Stinespring isometry sum_l |l> tensor K_l, every positive record effect E "
            "has Choi range contained in span{vec K_l}. Positivity of J_source_effect-p J_Phi for p>0 "
            "would require range(J_Phi) to lie in that source span. Exact ranks disprove this."
        ),
        "normalizations": {
            "hodge": "Phi=(1/3) sum_m Ad(C_m)",
            "marked": "normalized actual jump channel=(1/12) sum_{l in moving(mark)} Ad(U_l D_jump,l)",
            "whole": "C60=(1/60) sum_l Ad(U_l)",
            "note": "positive scalar normalizations do not alter Choi support",
        },
        "marked_twenty_records": marked,
        "actual_unmarked_sixty_record_in_each_alternative_chart": whole,
        "two_event_actual_whole_sixty_record": two_event,
        "chart_boundary": (
            "The six entries are coordinate tests of the same actual sixty matrices. They were not "
            "combined into an artificial 360-event source. Without a fixed duad chart, the old Bell10 "
            "and Lambda^2 C5 Hodge channel are not linearly identified."
        ),
        "source_boundary": (
            "The old coherent-record contracts establish only conditional record constructions and "
            "explicitly leave physical record access and dynamics underived. This calculation grants "
            "coherent record access as the changed premise. The exact obstruction applies at one event; "
            "two events become sufficient only after a target-specific coherent record filter is chosen."
        ),
        "holomorphic_product_boundary": (
            "The same five epsilon tensors span the image of the wedge map Sym^2(Lambda^2 C5) to "
            "Lambda^4 C5, but they are the five Plucker relations, hence the kernel of holomorphic "
            "coordinate multiplication Sym^2(10) to the 50-dimensional degree-two coordinate ring. "
            "The record subchannel is therefore not automatically an existing holomorphic product output."
        ),
        "plucker_wedge_typing": product_typing,
        "one_event_extra_operation_needed": (
            "At one event, producing full Phi would require a system operation that adds the three "
            "missing Hodge Kraus directions or a different source. A record-space projection alone "
            "cannot add them."
        ),
        "two_event_selection_boundary": (
            "Two actual source events already span M10, so a coherent record contraction can realize "
            "the exact subchannel (1/25)Phi. The minimum-norm contraction is fixed after a duad chart "
            "and the Hodge target are supplied. Neither that target-specific record effect nor a physical "
            "heralding apparatus is selected by the original source."
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
    payload = main(args.out)
    print(json.dumps(payload, indent=2, sort_keys=True))
