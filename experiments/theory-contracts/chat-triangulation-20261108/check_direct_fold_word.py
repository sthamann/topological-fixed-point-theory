"""Exact two vertex modes and the pulled-back Glue phase on the original Fold chart."""

import ast
import json
import sys
from pathlib import Path

import sympy as sp

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo))
from tfpt_explorer.current_block_geometry import _actual_current_source

fold = repo / "experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json"
M = sp.Matrix(json.loads(fold.read_text())["explicit_Weyl_matrix"])
tree = ast.parse((repo / "verification/v634_st31_structure.py").read_text())
idx = next(
    ast.literal_eval(node.value)
    for node in tree.body
    if isinstance(node, ast.Assign)
    and any(isinstance(t, ast.Name) and t.id == "SIGMA_IDX" for t in node.targets)
)
S = sp.eye(8)[idx, :]
assert M.T * M == sp.eye(8)
assert S**3 == sp.eye(8)

beta = sp.eye(8)[:, 2] + sp.eye(8)[:, 3]
target = [S**j * beta for j in range(3)]
p = [M.T * b for b in target]
delta = [p[j + 1] - p[j] for j in range(2)]
assert all(x.dot(x) == 4 for x in delta)
assert all(delta[j].dot(p[j]) == -2 for j in range(2))
assert all(all(x.q == 2 for x in step) and sum(step) % 2 == 0 for step in delta)

a = M.T * sp.Matrix([1, 1, 1, 1, 1, 0, 0, 0]) / 2
indices = sp.Matrix(range(8))
shift = a - indices / 4 - sp.ones(8, 1) / 2
assert all(x.q == 1 for x in shift) and sum(shift) % 2 == 0
assert [a.dot(x) for x in p] == [1, sp.Rational(1, 2), 1]
characters = [sp.simplify(sp.exp(2 * sp.pi * sp.I * x)) for x in a]
assert characters == [-1, -sp.I, 1, sp.I] * 2

ch = _actual_current_source()["chevalley"]


def cocycle(left, right):
    u, v = (ch.rs.alpha_coords(tuple(int(2 * x) for x in b)) for b in (left, right))
    assert all(x.denominator == 1 for x in (*u, *v))
    power = sum(int(u[i]) * ch.B[i][j] * int(v[j]) for i in range(8) for j in range(8))
    return (-1)**(power % 2)


for root in ch.ridx:
    b = sp.Matrix(root) / 2
    first_five = b[:5]
    if root[0] % 2 == 0:
        grade = 2 * (sum(first_five) % 2)
    else:
        grade = 1 + 2 * (sum(x - sp.Rational(1, 2) for x in first_five) % 2)
    charge = M.T * b
    assert sp.simplify(sp.exp(2 * sp.pi * sp.I * a.dot(charge))) == sp.I**grade

signs = [cocycle(target[j + 1] - target[j], target[j]) for j in range(2)]
return_sign = cocycle(target[0] - target[2], target[2])
cases = []
for alpha in (sp.Rational(1, 4), sp.Rational(3, 4)):
    background = (indices + (alpha - 4) * sp.ones(8, 1)) / 8
    k = [x.dot(x) / 2 + background.dot(x) for x in p]
    source_translation = [8 * e - alpha * sum(x) for e, x in zip(k, p)]
    assert source_translation == [4, 6, 4]
    assert all(
        sp.simplify(sp.exp(sp.pi * sp.I * j / 2)) == sp.exp(2 * sp.pi * sp.I * a.dot(x))
        for j, x in zip(source_translation, p)
    )
    cases.append({
        "quarter": str(alpha),
        "K": list(map(str, k)),
        "phase_first": -1,
        "phase_second": 1,
        "equal_strips_reference_normalized_exponent": str(k[1] + k[2] - 2 * k[0]),
        "source_translation_generator": list(map(str, source_translation)),
        "same_state_first_stress_weight": str((p[0] + background).dot(p[0] + background)),
    })

print(json.dumps({
    "source_orbit": [list(map(str, x)) for x in p],
    "vertex_charges": [list(map(str, x)) for x in delta],
    "vertex_weight": 2,
    "vertex_mode": 0,
    "ordered_vertex_signs": signs,
    "reverse_even_vertex_sign": return_sign,
    "charge_phase_character": list(map(str, a)),
    "residue_mode_action": list(map(str, characters)),
    "character_checks_on_original_roots": len(ch.ridx),
    "cases": cases,
}, indent=2))
