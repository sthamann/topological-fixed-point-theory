#!/usr/bin/env python3
"""Exact comparison of the old five-slot Hodge table with the A4 carrier maps.

This is a finite convention check, not a source or dynamics model.
"""

from itertools import combinations
import json
from pathlib import Path


def perm_sign(word):
    if len(set(word)) != len(word):
        return 0
    inversions = sum(
        word[i] > word[j]
        for i in range(len(word))
        for j in range(i + 1, len(word))
    )
    return -1 if inversions % 2 else 1


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


slots = tuple(range(5))
pairs = tuple(combinations(slots, 2))
pair_index = {pair: index for index, pair in enumerate(pairs)}

comparisons = []
matrices = []
for omitted in slots:
    rest = tuple(i for i in slots if i != omitted)
    volume_character = perm_sign((omitted,) + rest)
    matrix = [[0 for _ in pairs] for _ in pairs]
    unordered_factorizations = set()
    for input_pair in combinations(rest, 2):
        output_pair = tuple(i for i in rest if i not in input_pair)
        old_wedge_sign = perm_sign(input_pair + output_pair)
        new_epsilon_sign = perm_sign((omitted,) + input_pair + output_pair)
        require(
            new_epsilon_sign == volume_character * old_wedge_sign,
            "old Hodge sign and new epsilon sign disagree",
        )
        matrix[pair_index[output_pair]][pair_index[input_pair]] = new_epsilon_sign
        unordered_factorizations.add(tuple(sorted((input_pair, output_pair))))
        comparisons.append(
            {
                "omitted": omitted,
                "input_pair": input_pair,
                "output_pair": output_pair,
                "old_wedge_sign": old_wedge_sign,
                "volume_character": volume_character,
                "new_epsilon_sign": new_epsilon_sign,
            }
        )
    require(len(unordered_factorizations) == 3, "wrong doily factorization count")
    matrices.append(matrix)

require(len(comparisons) == 30, "expected 30 directed Hodge entries")
require(
    sum(1 for row in matrices for entries in row for x in entries if x) == 30,
    "expected 30 nonzero matrix entries",
)

# Exact gain: each carrier pair excludes three of the five slots.
gain = [[0 for _ in pairs] for _ in pairs]
for matrix in matrices:
    for i in range(len(pairs)):
        for j in range(len(pairs)):
            gain[i][j] += sum(matrix[k][i] * matrix[k][j] for k in range(len(pairs)))
require(gain == [[3 if i == j else 0 for j in range(10)] for i in range(10)],
        "sum C_m^T C_m is not 3 I_10")

result = {
    "status": "PASS",
    "scope": "finite sign and support comparison only",
    "old_doily_unordered_lines": 15,
    "directed_matrix_entries": 30,
    "volume_character_by_omitted_slot": [
        perm_sign((m,) + tuple(i for i in slots if i != m)) for m in slots
    ],
    "identity": "epsilon(m,i,j,k,l)=chi(m)*sign((i,j,k,l))",
    "gain": "sum_m C_m^T C_m = 3 I_10",
    "comparisons": comparisons,
}

out = Path(__file__).with_name("old-hodge-to-a4-results.json")
out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps({k: v for k, v in result.items() if k != "comparisons"}, sort_keys=True))
