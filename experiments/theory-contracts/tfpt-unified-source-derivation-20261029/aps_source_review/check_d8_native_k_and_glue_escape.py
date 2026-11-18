#!/usr/bin/env python3
"""Exact native-D5 and marked-glue escapes from the three-family D8."""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tfpt_explorer.current_block_geometry import _actual_current_source


SOURCE = ROOT / "tfpt_explorer/current_block_geometry.py"
V1 = ROOT / "verification/v1_e8_glue.py"
V983 = ROOT / "verification/v983_simple_current_generator.py"
V474 = ROOT / "verification/v474_entropic_hodge_carrier.py"
SPIN_LIFT = ROOT / "tfpt_explorer/source_spin_lift.py"
FAMILY = [(-1, -1, -1), (-1, 1, 1), (1, -1, 1), (1, 1, -1)]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def calculate() -> dict[str, object]:
    source = _actual_current_source()
    algebra = source["chevalley"]

    names: dict[int, str] = {}
    for a in (1, 2, 3):
        roots = [(-1,) * 5 + FAMILY[a]]
        roots += [
            tuple(-1 if j == i else 1 for j in range(5)) + FAMILY[a]
            for i in range(5)
        ]
        labels = [f"C{a}"] + [f"X{i},{a}" for i in range(5)]
        for root, label in zip(roots, labels):
            index = algebra.ridx[root]
            names[index] = label
            names[algebra.opp[index]] = label + "^dag"

    def close(start: set[int]) -> tuple[set[int], list[int]]:
        closed = set(start)
        counts = [len(closed)]
        while True:
            new = {
                target
                for left in closed
                for right in closed
                for target, coefficient in algebra.bracket(left, right).items()
                if target < algebra.nR and coefficient and target not in closed
            }
            if not new:
                return closed, counts
            closed.update(new)
            counts.append(len(closed))

    seed = set(names)
    d8, d8_counts = close(seed)
    d5 = {
        index
        for index, root in enumerate(algebra.roots)
        if any(root[:5]) and all(value == 0 for value in root[5:])
    }
    target = {
        index
        for index, root in enumerate(algebra.roots)
        if sum(value != 0 for value in root[:5]) == 1
        and sum(value != 0 for value in root[5:]) == 1
        and all(value in (-2, 0, 2) for value in root)
    }
    with_d5, with_d5_counts = close(seed | d5)

    # One exact word.  The first bracket is already in the restricted D8.
    x13_dag = algebra.opp[
        algebra.ridx[(1, -1, 1, 1, 1) + FAMILY[3]]
    ]
    c2_dag = algebra.opp[algebra.ridx[(-1,) * 5 + FAMILY[2]]]
    y_root = (0, 2, 0, 0, 0, -2, 0, 0)
    y = algebra.ridx[y_root]
    d_root = (-2, -2, 0, 0, 0, 0, 0, 0)
    d = algebra.ridx[d_root]
    z_root = (-2, 0, 0, 0, 0, -2, 0, 0)
    z = algebra.ridx[z_root]
    first_coefficient = algebra.bracket(x13_dag, c2_dag).get(y, 0)
    second_coefficient = algebra.bracket(d, y).get(z, 0)

    # The marked v983 generator maps from A3 sum-zero coordinates to the
    # doubled current root (1,...,1), which is exactly C0^dag here.
    a3_to_d3 = sp.Rational(1, 2) * sp.Matrix([
        [1, 1, -1, -1],
        [1, -1, 1, -1],
        [1, -1, -1, 1],
    ])
    omega_f = sp.Matrix([
        sp.Rational(3, 4), sp.Rational(-1, 4),
        sp.Rational(-1, 4), sp.Rational(-1, 4),
    ])
    mapped_omega_f = tuple(a3_to_d3 * omega_f)
    lambda_root = (1,) * 8
    lambda_index = algebra.ridx[lambda_root]
    c0 = algebra.ridx[(-1,) * 8]
    with_lambda, with_lambda_counts = close(d8 | {lambda_index, c0})

    grade = lambda root: (-sum(root[5:])) % 4
    d8_grades = Counter(grade(algebra.roots[index]) for index in d8)
    complement = set(range(algebra.nR)) - d8
    complement_grades = Counter(grade(algebra.roots[index]) for index in complement)
    lambda_pairings = [
        sum(a * b for a, b in zip(lambda_root, algebra.roots[index]))
        for index in d8
    ]

    checks = {
        "restricted_closure_is_D8_112": len(d8) == 112 and d8_counts == [36, 92, 112],
        "restricted_D8_has_only_half_of_original_D5": len(d5 & d8) == 20 and len(d5 - d8) == 20,
        "restricted_D8_reaches_only_half_of_10x6": len(target & d8) == 30 and len(target - d8) == 30,
        "first_word_is_actual_Chevalley_bracket": first_coefficient == -1 and y in d8,
        "native_D5_root_maps_to_missing_target": second_coefficient == 1 and d not in d8 and z not in d8,
        "full_original_D5_action_closes_E8": len(with_d5) == 240 and with_d5_counts == [76, 172, 232, 240],
        "marked_lambda_maps_to_all_plus_root": mapped_omega_f == (sp.Rational(1, 2),) * 3,
        "marked_lambda_is_C0_adjoint": algebra.opp[c0] == lambda_index,
        "lambda_pair_leaves_D8_and_closes_E8": lambda_index not in d8 and c0 not in d8 and len(with_lambda) == 240,
        "lambda_is_integral_on_new_D8": all(value % 4 == 0 for value in lambda_pairings),
        "native_D5_word_preserves_original_clock_grade": grade(d_root) == 0 and grade(y_root) == grade(z_root) == 2,
        "marked_lambda_retains_original_quarter_grades": grade(lambda_root) == 1 and grade(tuple(-x for x in lambda_root)) == 3,
        "original_quarter_grade_is_not_new_D8_coset_parity": set(d8_grades) == {0, 1, 2, 3} and set(complement_grades) == {0, 1, 2, 3},
    }

    return {
        "title": "Native D5 and marked-glue escapes from the three-family D8",
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "restricted_D8": {
            "closure_counts": d8_counts,
            "D5_roots_inside_outside": [len(d5 & d8), len(d5 - d8)],
            "target_10x6_inside_outside": [len(target & d8), len(target - d8)],
            "original_mu4_grade_counts": {str(k): d8_grades[k] for k in sorted(d8_grades)},
            "complement_mu4_grade_counts": {str(k): complement_grades[k] for k in sorted(complement_grades)},
        },
        "native_D5_word": {
            "word": "[E_(-e0-e1), [X_(1,3)^dag, C_2^dag]]",
            "inner_root_doubled": list(y_root),
            "inner_coefficient": str(first_coefficient),
            "D5_root_doubled": list(d_root),
            "outer_coefficient": str(second_coefficient),
            "output_root_doubled": list(z_root),
            "clock_grades": [grade(d_root), grade(y_root), grade(z_root)],
            "closure_with_full_D5_counts": with_d5_counts,
        },
        "marked_glue": {
            "omega_f_A3": [str(value) for value in omega_f],
            "omega_f_D3": [str(value) for value in mapped_omega_f],
            "lambda_root_doubled": list(lambda_root),
            "equals": "C0^dag",
            "adjoint_root_doubled": [-1] * 8,
            "closure_with_pair_counts": with_lambda_counts,
            "pairing_numerators_over_four": dict(sorted(Counter(lambda_pairings).items())),
            "clock_grades": [grade(lambda_root), grade(tuple(-x for x in lambda_root))],
        },
        "scope": (
            "Exact current-algebra and lattice-VOA typing in the original Chevalley chart. "
            "It proves the native D5 zero-mode and the marked glue vertex are algebraic escapes. "
            "It does not derive their execution, state preparation or response policy from primitive P1."
        ),
        "sources": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (SOURCE, V1, V983, V474, SPIN_LIFT)
        },
    }


def main() -> None:
    output = Path(__file__).with_name("d8-native-k-and-glue-escape-results.json")
    output.write_text(json.dumps(calculate(), indent=2, sort_keys=True) + "\n")
    print(output)


if __name__ == "__main__":
    main()
