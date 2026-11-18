"""Exact residual symmetry of the existing marked compiler dictionary."""
import ast
import contextlib
import io
import itertools as it
import json
import math
import runpy
from hashlib import sha256
from pathlib import Path

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PINS = {}


def need(condition, label):
    if not bool(condition):
        raise RuntimeError(label)


def native(relative, names, namespace):
    path = REPO / relative
    data = path.read_bytes()
    PINS[relative] = sha256(data).hexdigest()
    nodes = []
    for node in ast.parse(data).body:
        if isinstance(node, ast.FunctionDef) and node.name in names:
            nodes.append(node)
        if isinstance(node, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id in names for t in node.targets
        ):
            nodes.append(node)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), "exec"), namespace)
    need(all(name in namespace for name in names), "native source definitions")
    return namespace


ns = native("verification/v634_st31_structure.py", {"J_vec", "SIGMA_IDX"}, {"np": np})
eye = np.eye(8, dtype=np.int64)
J = np.column_stack([ns["J_vec"](x) for x in eye])
sigma = eye[np.array(ns["SIGMA_IDX"])]
c = J @ sigma
Y6 = np.array([-2, -2, -2, 3, 3, 0, 0, 0], dtype=np.int64)

# Every H element is a block signed permutation with even sign parity in
# each block. The early source-pair test is necessary for commuting with J.
centralizer = []
marked = []
for p5 in it.permutations(range(5)):
    for p3 in it.permutations(range(5, 8)):
        perm = p5 + p3
        if any(perm[k] // 2 != perm[k + 1] // 2 for k in range(0, 8, 2)):
            continue
        for signs in it.product((-1, 1), repeat=8):
            if math.prod(signs[:5]) != 1 or math.prod(signs[5:]) != 1:
                continue
            M = np.zeros((8, 8), dtype=np.int64)
            M[np.arange(8), list(perm)] = signs
            if np.array_equal(M @ J, J @ M):
                centralizer.append(M)
                if np.array_equal(M @ Y6, Y6):
                    marked.append(M)
need(len(centralizer) == 64, "complete C_H(J) order")
need(len(marked) == 2, "complete C_H(J) fixing Y")
d = np.diag([1, 1, 1, 1, 1, 1, -1, -1]).astype(np.int64)
need({tuple(M.flat) for M in marked} == {tuple(eye.flat), tuple(d.flat)}, "exact marked stabilizer")
need(np.array_equal(d @ sigma, sigma @ d) and np.array_equal(d @ c, c @ d), "all named compiler clocks")

# Four projective axes of the native family tetrahedron. Normalize first
# nonzero sign to +1; no root or vertex is selected by this normalization.
axes = tuple((1, a, b) for a, b in it.product((-1, 1), repeat=2))


def axis_image(M, axis):
    v = np.array((0,) * 5 + axis, dtype=np.int64)
    out = (M @ v)[5:]
    return tuple(int(x * out[0]) for x in out)


axis_perm = tuple(axes.index(axis_image(d, a)) for a in axes)
need(all(axis_perm[i] != i for i in range(4)), "no tetrahedral fixed vertex")
need(all(axis_perm[axis_perm[i]] == i for i in range(4)), "two transpositions")
family_actions = {
    tuple(axes.index(axis_image(M, a)) for a in axes) for M in centralizer
}
need(len(family_actions) == 4, "C_H(J) family image has order4")
need({p[0] for p in family_actions} == set(range(4)), "transitive four-axis action")
need(any(tuple(p[p[i]] for i in range(4)) != tuple(range(4))
         and tuple(p[p[p[p[i]]]] for i in range(4)) == tuple(range(4))
         for p in family_actions), "family image contains an order4 element and is cyclic")

# Use the actual native sparse quartic supports, preserving subspace and
# component transformation law rather than requiring each vector fixed.
pn = native("tfpt_explorer/process.py", {"_code_basis"}, {"np": np, "it": it, "math": math})
_, supports = pn["_code_basis"]()
z = sp.symbols("z0:4")
polys = [sum(sp.prod(z[j] for j in word) for word in support) for support in supports]
dz = (1, 1, 1, -1)
sub = {z[j]: dz[j] * z[j] for j in range(4)}
quartic_signs = []
for p in polys:
    image = sp.expand(p.subs(sub, simultaneous=True))
    candidates = [s for s in (-1, 1) if sp.expand(image - s * p) == 0]
    need(len(candidates) == 1, "actual quartic action")
    quartic_signs.append(candidates[0])
need(quartic_signs == [1, 1, 1, 1, -1], "literal P5 action of d")
for i, p in enumerate(polys):
    for j in range(4):
        derivative = sp.diff(p, z[j])
        need(sp.expand(derivative.subs(sub, simultaneous=True)
                       - quartic_signs[i] * dz[j] * derivative) == 0,
             "all20 actual erasure derivatives transform equivariantly")

# Check the first possible extra selector: the ACTUAL full-field lift phases.
# Reuse the original joint-clock execution path and its unchanged epsilon.
clock_path = HERE / "joint_clock.py"
PINS[str(clock_path.relative_to(REPO))] = sha256(clock_path.read_bytes()).hexdigest()
with contextlib.redirect_stdout(io.StringIO()):
    clock = runpy.run_path(str(clock_path))
d_lattice = clock["Bi"] * sp.Matrix(d) * clock["B"]
compose = clock["composed_basis_phases"]
power = clock["power_basis_phases"]
commute_J = []
commute_J_c = []
involutions = []
for mask in range(256):
    if compose(d_lattice, mask, clock["J"], 0) != compose(clock["J"], 0, d_lattice, mask):
        continue
    commute_J.append(mask)
    if compose(d_lattice, mask, clock["c"], 31) != compose(clock["c"], 31, d_lattice, mask):
        continue
    commute_J_c.append(mask)
    if power(d_lattice, mask, 2) == (1,) * 8:
        involutions.append(mask)
need(len(commute_J) == 16, "all d lifts commuting with fixed Jhat")
need(commute_J_c == [126, 127, 128, 129], "all d lifts commuting with fixed Jhat and chat31")
need(involutions == [126, 127], "two exact involutive full-field d automorphisms")

for relative, digest in PINS.items():
    need(sha256((REPO / relative).read_bytes()).hexdigest() == digest, "source unchanged")

result = {
    "research_id": "MATH.SEAM.JOINT.SELECTION.05",
    "verdict": "PASS",
    "centralizer_H_J_order": len(centralizer),
    "family_permutation_image_order": len(family_actions),
    "centralizer_H_J_Y_order": len(marked),
    "residual_diagonal": list(np.diag(d).astype(int).tolist()),
    "residual_axis_permutation": axis_perm,
    "quartic_component_signs": quartic_signs,
    "erasure_components_checked": 20,
    "d_lifts_commuting_with_fixed_Jhat": commute_J,
    "d_lifts_commuting_with_fixed_Jhat_and_chat31": commute_J_c,
    "involutive_full_field_d_lifts": involutions,
    "scope": "No equivariant unique E6 vertex selection from the stated compiler diagram, even with fixed Jhat and chat31 full-field phases. Character gauge conjugates the result to the other cyclic lift representatives. A separately fixed individual quartic vector is a different marking.",
    "source_sha256": PINS,
}
print(json.dumps(result, indent=2))

