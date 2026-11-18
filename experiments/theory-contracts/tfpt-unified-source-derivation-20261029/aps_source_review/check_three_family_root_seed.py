#!/usr/bin/env python3
"""Exact Lie closure of the actual C/X seed restricted to three families."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import hashlib
import itertools as it
import json
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tfpt_explorer.current_block_geometry import _actual_current_source


SOURCE = ROOT / "tfpt_explorer/current_block_geometry.py"
FAMILY_WEIGHTS = [
    (-1, -1, -1),
    (-1, 1, 1),
    (1, -1, 1),
    (1, 1, -1),
]
F_LABELS = (1, 2, 3)


def calculate() -> dict[str, object]:
    source = _actual_current_source()
    chevalley = source["chevalley"]

    c_roots = {
        a: chevalley.ridx[(-1,) * 5 + FAMILY_WEIGHTS[a]] for a in F_LABELS
    }
    x_roots = {
        (i, a): chevalley.ridx[
            tuple(-1 if j == i else 1 for j in range(5))
            + FAMILY_WEIGHTS[a]
        ]
        for i in range(5)
        for a in F_LABELS
    }

    names: dict[int, str] = {}
    for a, index in c_roots.items():
        names[index] = f"C{a}"
        names[chevalley.opp[index]] = f"C{a}^dag"
    for (i, a), index in x_roots.items():
        names[index] = f"X{i},{a}"
        names[chevalley.opp[index]] = f"X{i},{a}^dag"

    seed = set(names)
    levels = {index: 0 for index in seed}
    words: dict[int, tuple[object, ...]] = {
        index: ("seed", names[index]) for index in seed
    }
    generated = set(seed)
    frontier = set(seed)
    closure_counts = [len(generated)]

    depth = 0
    while frontier:
        depth += 1
        new: dict[int, tuple[int, int, Fraction]] = {}
        for left in sorted(frontier):
            for right in sorted(generated):
                for target, coefficient in sorted(
                    chevalley.bracket(left, right).items()
                ):
                    if (
                        target < chevalley.nR
                        and coefficient
                        and target not in generated
                        and target not in new
                    ):
                        new[target] = (left, right, coefficient)
        if not new:
            break
        for target, (left, right, coefficient) in new.items():
            levels[target] = depth
            words[target] = ("bracket", left, right, coefficient)
        generated.update(new)
        frontier = set(new)
        closure_counts.append(len(generated))

    def root(index: int) -> tuple[int, ...]:
        return tuple(int(x) for x in chevalley.roots[index])

    def certificate(index: int) -> dict[str, object]:
        word = words[index]
        if word[0] == "seed":
            return {
                "kind": "seed",
                "name": word[1],
                "root": list(root(index)),
            }
        _, left, right, coefficient = word
        return {
            "kind": "bracket",
            "root": list(root(index)),
            "depth": levels[index],
            "coefficient": str(coefficient),
            "left": certificate(int(left)),
            "right": certificate(int(right)),
        }

    odd_spinors = {
        signs
        for signs in it.product((-1, 1), repeat=5)
        if sum(value == -1 for value in signs) % 2 == 1
    }
    cx_spinors = {(-1,) * 5} | {
        tuple(-1 if j == i else 1 for j in range(5)) for i in range(5)
    }
    missing_spinors = odd_spinors - cx_spinors
    missing_10xf = {
        chevalley.ridx[signs + FAMILY_WEIGHTS[a]]
        for signs in missing_spinors
        for a in F_LABELS
    }
    missing_10xf_dag = {chevalley.opp[index] for index in missing_10xf}
    generated_10xo = {
        chevalley.ridx[signs + FAMILY_WEIGHTS[0]] for signs in missing_spinors
    }
    generated_10xo |= {chevalley.opp[index] for index in list(generated_10xo)}
    depth_two_roots = {index for index in generated if levels[index] == 2}

    target_10x6 = {
        index
        for index, coordinates in enumerate(chevalley.roots)
        if sum(value != 0 for value in coordinates[:5]) == 1
        and sum(value != 0 for value in coordinates[5:]) == 1
        and all(value in (-2, 0, 2) for value in coordinates)
    }
    target_generated = target_10x6 & generated
    target_missing = target_10x6 - generated

    root_rank = int(sp.Matrix([root(index) for index in generated]).rank())
    bracket_closed = all(
        target >= chevalley.nR or not coefficient or target in generated
        for left in generated
        for right in generated
        for target, coefficient in chevalley.bracket(left, right).items()
    )

    chamber_vector = (1, 2, 4, 8, 16, 32, 64, 128)
    generated_coordinates = {root(index) for index in generated}
    positive_roots = {
        coordinates
        for coordinates in generated_coordinates
        if sum(a * b for a, b in zip(coordinates, chamber_vector)) > 0
    }
    simple_roots = sorted(
        coordinates
        for coordinates in positive_roots
        if not any(
            tuple(a - b for a, b in zip(coordinates, candidate)) in positive_roots
            for candidate in positive_roots
        )
    )
    cartan = sp.Matrix(
        [
            [
                sum(a * b for a, b in zip(left, right)) // 4
                for right in simple_roots
            ]
            for left in simple_roots
        ]
    )
    dynkin_edges = [
        [left, right]
        for left in range(len(simple_roots))
        for right in range(left)
        if cartan[left, right] == -1
    ]
    dynkin_degrees = sorted(
        sum(vertex in edge for edge in dynkin_edges)
        for vertex in range(len(simple_roots))
    )

    def close_roots(start: set[int]) -> tuple[set[int], list[int]]:
        closed = set(start)
        active = set(start)
        counts = [len(closed)]
        while active:
            new = {
                target
                for left in active
                for right in closed
                for target, coefficient in chevalley.bracket(left, right).items()
                if target < chevalley.nR
                and coefficient
                and target not in closed
            }
            if not new:
                break
            closed.update(new)
            active = new
            counts.append(len(closed))
        return closed, counts

    omitted_family_candidates = {
        "C0": chevalley.ridx[(-1,) * 5 + FAMILY_WEIGHTS[0]],
        **{
            f"X{i},0": chevalley.ridx[
                tuple(-1 if j == i else 1 for j in range(5))
                + FAMILY_WEIGHTS[0]
            ]
            for i in range(5)
        },
    }
    one_root_extensions = {}
    compact_pair_extensions = {}
    for name, index in omitted_family_candidates.items():
        one_closed, one_counts = close_roots(generated | {index})
        pair_closed, pair_counts = close_roots(
            generated | {index, chevalley.opp[index]}
        )
        one_root_extensions[name] = {
            "root": list(root(index)),
            "closure_counts": one_counts,
            "final_root_count": len(one_closed),
        }
        compact_pair_extensions[name] = {
            "roots": [list(root(index)), list(root(chevalley.opp[index]))],
            "closure_counts": pair_counts,
            "final_root_count": len(pair_closed),
        }

    target_examples = sorted(target_generated, key=lambda index: (levels[index], root(index)))[:3]
    depth_two_examples = sorted(
        (index for index in generated if levels[index] == 2), key=root
    )[:2]

    checks = {
        "seed_is_18_positive_plus_18_adjoint_roots": len(seed) == 36,
        "closure_counts_are_36_92_112": closure_counts == [36, 92, 112],
        "closure_is_exactly_bracket_closed": bracket_closed,
        "closed_root_rank_is_eight": root_rank == 8,
        "closed_root_count_is_D8_112": len(generated) == 112,
        "D8_has_eight_simple_roots": len(simple_roots) == 8,
        "D8_Cartan_determinant_is_four": int(cartan.det()) == 4,
        "D8_Dynkin_degree_sequence": dynkin_degrees == [1, 1, 1, 2, 2, 2, 2, 3],
        "not_full_E8": len(generated) != 240,
        "E8_complement_has_128_roots": chevalley.nR - len(generated) == 128,
        "missing_10xF_positive_count_is_30": len(missing_10xf) == 30,
        "missing_10xF_positive_unreachable": not (missing_10xf & generated),
        "missing_10xF_adjoint_unreachable": not (missing_10xf_dag & generated),
        "depth_two_is_exactly_10xO_plus_adjoint": depth_two_roots == generated_10xo,
        "full_10x6_target_has_60_roots": len(target_10x6) == 60,
        "only_30_target_roots_generated": len(target_generated) == 30,
        "other_30_target_roots_unreachable": len(target_missing) == 30,
        "all_generated_target_roots_have_depth_one": {
            levels[index] for index in target_generated
        } == {1},
        "every_single_omitted_positive_root_generates_E8": all(
            row["final_root_count"] == 240 for row in one_root_extensions.values()
        ),
        "every_single_omitted_compact_pair_generates_E8": all(
            row["final_root_count"] == 240 for row in compact_pair_extensions.values()
        ),
    }
    if not all(checks.values()):
        failed = [name for name, ok in checks.items() if not ok]
        raise AssertionError(f"failed checks: {failed}")

    return {
        "scope": "Actual Chevalley root closure of C_a and X_i,a for a in F={1,2,3}, including their compact adjoint roots. This is a current-algebra closure, not a source-access, affine-grade, state, time or finite-Phi claim.",
        "source": {
            "path": str(SOURCE.relative_to(ROOT)),
            "sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        },
        "seed": {
            "positive_C": 3,
            "positive_X": 15,
            "positive_total": 18,
            "with_adjoint_roots": 36,
            "family_labels": list(F_LABELS),
        },
        "closure": {
            "root_counts_by_depth": closure_counts,
            "final_root_count": len(generated),
            "root_rank": root_rank,
            "bracket_closed": bracket_closed,
            "subsystem_type": "D8",
            "full_E8_generated": False,
            "E8_complement_root_count": chevalley.nR - len(generated),
            "D8_Weyl_frame": {
                "chamber_vector": list(chamber_vector),
                "simple_roots": [list(row) for row in simple_roots],
                "Cartan_matrix": [list(map(int, row)) for row in cartan.tolist()],
                "Cartan_determinant": int(cartan.det()),
                "Dynkin_edges": dynkin_edges,
                "Dynkin_degree_sequence": dynkin_degrees,
            },
        },
        "missing_10_tensor_F_spinors": {
            "positive_count": len(missing_10xf),
            "adjoint_count": len(missing_10xf_dag),
            "positive_depth": "unreachable",
            "adjoint_depth": "unreachable",
        },
        "generated_at_depth_two": {
            "type": "10 tensor O plus compact adjoint",
            "root_count": len(depth_two_roots),
        },
        "target_10_tensor_6": {
            "full_root_count": len(target_10x6),
            "generated_count": len(target_generated),
            "generated_depth_distribution": {
                str(level): count
                for level, count in sorted(
                    Counter(levels[index] for index in target_generated).items()
                )
            },
            "unreachable_count": len(target_missing),
        },
        "minimal_E8_escape_from_omitted_family": {
            "complex_one_positive_root": one_root_extensions,
            "compact_root_plus_adjoint": compact_pair_extensions,
            "meaning": "The D8 closure is maximal in this E8 root realization. Any one omitted family-zero C/X root leaves the D8 subsystem and generates the 128-root complement. A compact source adds the root together with its adjoint.",
        },
        "word_certificates": {
            "generated_10_tensor_6_depth_one": [
                certificate(index) for index in target_examples
            ],
            "other_D8_roots_at_depth_two": [
                certificate(index) for index in depth_two_examples
            ],
        },
        "checks": checks,
        "check_count": len(checks),
        "verdict": "The three-family C/X seed plus compact adjoints closes to an explicitly certified rank-eight D8 subsystem with 112 roots, not E8. It never generates the missing 10 tensor F spinor components. It reaches 30 of the 60 (10,6) roots at depth one; the complementary 30 are unreachable. Any one omitted family-zero C/X root, or its compact root-adjoint pair, is a minimal algebraic escape and generates E8. The full-family Hodge and full-E8 source-closure results therefore require an additional actual source-access arrow outside the restricted seed.",
    }


def main() -> None:
    result = calculate()
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if len(sys.argv) > 1:
        Path(sys.argv[1]).write_text(text)
    else:
        sys.stdout.write(text)


if __name__ == "__main__":
    main()

