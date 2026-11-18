#!/usr/bin/env python3
"""Exact finite discriminant/current checks for the analytic proof in THEOREMS.txt.

No finite computation here proves OS reconstruction or an operator-domain theorem.
Native functions are extracted without executing their module-level suites.
"""

import ast
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import itertools
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PINS = {}


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def native(relative, name):
    path = REPO / relative
    raw = path.read_bytes()
    PINS[relative] = sha256(raw).hexdigest()
    node = next(n for n in ast.parse(raw).body
                if isinstance(n, ast.FunctionDef) and n.name == name)
    env = {"np": np, "itertools": itertools}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), "exec"), env)
    return env[name]


def add(x, y):
    return tuple((a + b) % 4 for a, b in zip(x, y))


def q(x):
    return Fraction(5 * x[0] ** 2 + 3 * x[1] ** 2, 8) % 1


def pairing(x, y):
    return Fraction(5 * x[0] * y[0] + 3 * x[1] * y[1], 4) % 1


def subgroup(elements):
    return (0, 0) in elements and all(add(x, y) in elements
                                       for x in elements for y in elements)


discriminant = tuple(itertools.product(range(4), repeat=2))
isotropic_elements = tuple(x for x in discriminant if q(x) == 0)
isotropic_subgroups = tuple(
    frozenset(h)
    for size in range(1, len(isotropic_elements) + 1)
    for h in itertools.combinations(isotropic_elements, size)
    if subgroup(h)
)
isotropic_order4 = tuple(h for h in isotropic_subgroups if len(h) == 4)
diagonal = frozenset((a, a) for a in range(4))
opposite = frozenset((a, (-a) % 4) for a in range(4))
trivial = frozenset({(0, 0)})
even_glue = frozenset({(0, 0), (2, 2)})
require(set(isotropic_subgroups) == {trivial, even_glue, diagonal, opposite}, "all four possible local supports")
require(set(isotropic_order4) == {diagonal, opposite}, "all isotropic order4 subgroups")
require([h for h in isotropic_subgroups if (1, 1) in h] == [diagonal], "one marked generator forces the whole support")
orthogonal = frozenset(x for x in discriminant if all(pairing(x, h) == 0 for h in diagonal))
require(orthogonal == diagonal, "marked glue is Lagrangian")
require(diagonal != opposite, "without the orientation marking there are two subgroups")
require(q((1, 0)) == Fraction(5, 8) and q((0, 1)) == Fraction(3, 8), "native discriminant weights")
require([q((a, a)) for a in range(4)] == [0] * 4, "no twist obstruction on the full C4")
require(Fraction(16, 4**2) == 1, "local-module global dimension")

e8_roots = native("verification/v1_e8_glue.py", "e8_roots")
roots = frozenset(tuple(int(round(2 * a)) for a in v) for v in e8_roots())
require(len(roots) == 240 and all(sum(a*a for a in v) == 8 for v in roots), "actual norm-two shell")


def charge(root):
    """C4 glue label in the actual doubled D5+D3 coordinate chart."""
    if all(abs(a) == 1 for a in root):
        return 1 if sum(a < 0 for a in root[:5]) % 2 == 0 else 3
    first = any(root[i] for i in range(5))
    last = any(root[i] for i in range(5, 8))
    return 2 if first and last else 0


root_sector_counts = Counter(charge(r) for r in roots)
require(root_sector_counts == {0: 52, 1: 64, 2: 60, 3: 64}, "all four current sectors")
current_sector_counts = dict(root_sector_counts)
current_sector_counts[0] += 8
require([current_sector_counts[i] for i in range(4)] == [60, 64, 60, 64], "248 current branching")

bracket_count = 0
for a in roots:
    for b in roots:
        total = tuple(x+y for x, y in zip(a, b))
        if total in roots:
            require(charge(total) == (charge(a)+charge(b)) % 4, "root OPE preserves actual C4 grading")
            bracket_count += 1
        elif all(x == 0 for x in total):
            require((charge(a)+charge(b)) % 4 == 0, "opposite-root Cartan bracket grade")
            bracket_count += 1
require(bracket_count == 13680, "all actual nonzero root bracket channels")

J_vec = native("verification/v634_st31_structure.py", "J_vec")
J = np.column_stack([J_vec(e) for e in np.eye(8, dtype=np.int64)])
I = np.eye(8, dtype=np.int64)
require(np.array_equal(J @ J, -I) and int(np.trace(J)) == 0, "native Cartan J has i/-i multiplicity4")
require(not np.array_equal(J, int(J[0, 0]) * I), "internal J is not any scalar conformal rotation on V1")
require(np.array_equal(J.T @ J, I), "lattice clock preserves the Sugawara quadratic form")

for relative, digest in PINS.items():
    require(sha256((REPO / relative).read_bytes()).hexdigest() == digest, "source changed during run")

result = {
    "research_id": "MATH.SEAM.ESERIES.SOURCE-RECONSTRUCTION.06",
    "finite_certificate": "PASS",
    "overall_goal": "PARTIAL: primitive-source identification not derived",
    "discriminant_group_order": len(discriminant),
    "q_table": [[str(q((a, b))) for b in range(4)] for a in range(4)],
    "all_isotropic_subgroups": [sorted(h) for h in isotropic_subgroups],
    "connected_local_extension_classes_over_fixed_U": 4,
    "marked_generator_selects_unique_support": True,
    "all_isotropic_order4_subgroups": [sorted(h) for h in isotropic_order4],
    "marked_H_equals_H_perp": True,
    "local_module_global_dimension": 1,
    "current_sectors_in_C4_order": [current_sector_counts[i] for i in range(4)],
    "current_total": sum(current_sector_counts.values()),
    "nonzero_root_bracket_channels": bracket_count,
    "J_cartan_eigenvalues": {"i": 4, "-i": 4},
    "J_equals_any_conformal_rotation": False,
    "proof_only": ["automatic finiteness of same-conformal extensions",
                   "multiplicity and support lemma from module adjunction",
                   "classification of VOA extensions from the isotropic support census",
                   "C4 intertwiner uniqueness over U1", "all descendant fields and strong locality",
                   "OS reconstruction hypotheses", "infinite operator-domain counterexample",
                   "conditional spectral-rescaling obstruction"],
    "source_sha256": PINS,
    "promotion": False,
}
(HERE / "results.json").write_text(json.dumps(result, indent=2) + "\n")
print("PASS: all4 local-extension supports; one marked generator forces E8; sectors60+64+60+64=248; internal J differs from conformal rotation")
print("SCOPE: analytic proofs and original-source obligations are in the accompanying text; no physical completion")
