#!/usr/bin/env python3
"""Exact closure of the native D5+A3 cut dictionary under J and G31 pairs.

The calculation is deliberately orbit-sized.  It never enumerates W(D8) or
W(E8).  It extracts only the native E8 shell, simple roots, and J function;
the module-level verification suites are not executed.
"""

import ast
from collections import Counter, deque
from fractions import Fraction
from hashlib import sha256
import itertools
from math import factorial
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PINS = {}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def native(relative, name):
    path = REPO / relative
    raw = path.read_bytes()
    PINS[relative] = sha256(raw).hexdigest()
    tree = ast.parse(raw)
    fn = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == name
    )
    namespace = {"np": np, "itertools": itertools}
    exec(compile(ast.Module(body=[fn], type_ignores=[]), str(path), "exec"), namespace)
    return namespace[name]


def neg(x):
    return tuple(-a for a in x)


def add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def support(x):
    return frozenset(i for i, a in enumerate(x) if a)


def independent_roots():
    out = []
    for i, j in itertools.combinations(range(8), 2):
        for si in (-2, 2):
            for sj in (-2, 2):
                v = [0] * 8
                v[i], v[j] = si, sj
                out.append(tuple(v))
    for signs in itertools.product((-1, 1), repeat=8):
        if signs.count(-1) % 2 == 0:
            out.append(tuple(signs))
    return frozenset(out)


native_roots = native("verification/v1_e8_glue.py", "e8_roots")
ROOTS = frozenset(tuple(int(round(2 * a)) for a in root) for root in native_roots())
require(ROOTS == independent_roots() and len(ROOTS) == 240, "native E8 shell")
INTEGER_ROOTS = frozenset(root for root in ROOTS if sum(a != 0 for a in root) == 2)
HALF_ROOTS = frozenset(root for root in ROOTS if all(abs(a) == 1 for a in root))
require((len(INTEGER_ROOTS), len(HALF_ROOTS)) == (112, 128), "112+128 E8 partition")

native_simple = native("verification/v1_e8_glue.py", "e8_simple_roots")
SIMPLE_ROOTS = tuple(tuple(int(round(2 * a)) for a in root) for root in native_simple())
require(len(SIMPLE_ROOTS) == 8 and all(root in ROOTS for root in SIMPLE_ROOTS), "native E8 simple roots")

native_j = native("verification/v634_st31_structure.py", "J_vec")
J_MATRIX = np.column_stack([native_j(v) for v in np.eye(8, dtype=np.int64)]).astype(np.int64)
J_EXPECTED = np.zeros((8, 8), dtype=np.int64)
for k in range(0, 8, 2):
    J_EXPECTED[k, k + 1], J_EXPECTED[k + 1, k] = -1, 1
require(np.array_equal(J_MATRIX, J_EXPECTED), "native J matrix")
require(np.array_equal(J_MATRIX @ J_MATRIX, -np.eye(8, dtype=np.int64)), "J squared")


def matrix_action(matrix, x):
    return tuple(int(a) for a in matrix @ np.array(x, dtype=np.int64))


def reflection(alpha, x):
    numerator = dot(x, alpha)
    require(numerator % 4 == 0, "integral E8 reflection coefficient")
    coefficient = numerator // 4
    return tuple(x[i] - coefficient * alpha[i] for i in range(8))


def complex_pair_reflection(alpha, x):
    jalpha = matrix_action(J_MATRIX, alpha)
    # alpha and J alpha are orthogonal norm-two roots, so the factors commute.
    return reflection(alpha, reflection(jalpha, x))


def permutation_matrix(i, j):
    out = np.eye(8, dtype=np.int64)
    out[:, [i, j]] = out[:, [j, i]]
    return out


def sign_pair_matrix(i, j):
    out = np.eye(8, dtype=np.int64)
    out[i, i] = out[j, j] = -1
    return out


H_PERMUTATIONS = [permutation_matrix(i, i + 1) for i in range(4)]
H_PERMUTATIONS += [permutation_matrix(5, 6), permutation_matrix(6, 7)]
H_SIGNS = [sign_pair_matrix(i, i + 1) for i in range(4)]
H_SIGNS += [sign_pair_matrix(5, 6), sign_pair_matrix(6, 7)]
WD8_GENERATORS = H_PERMUTATIONS + H_SIGNS + [J_MATRIX]


def monomial_permutation(matrix):
    require(np.all(np.sum(np.abs(matrix), axis=0) == 1), "monomial columns")
    require(np.all(np.sum(np.abs(matrix), axis=1) == 1), "monomial rows")
    return tuple(int(np.flatnonzero(matrix[:, j])[0]) for j in range(8))


def even_signed(matrix):
    nonzero = matrix[matrix != 0]
    return len(nonzero) == 8 and int(np.prod(nonzero)) == 1


require(all(even_signed(g) for g in WD8_GENERATORS), "H and J lie in W(D8)")

# On coordinate permutations, J contributes the missing block crossing (4 5):
# |J|=(01)(23)(45)(67), while (01)(23)(67) already lies in H.
internal = permutation_matrix(0, 1) @ permutation_matrix(2, 3) @ permutation_matrix(6, 7)
crossing = monomial_permutation(J_MATRIX @ internal)
require(crossing == (0, 1, 2, 3, 5, 4, 6, 7), "J supplies the crossing transposition")


def gf2_rank(rows):
    rows = [sum((bit & 1) << i for i, bit in enumerate(row)) for row in rows]
    rank = 0
    for column in range(8):
        pivot = next((k for k in range(rank, len(rows)) if (rows[k] >> column) & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for k in range(len(rows)):
            if k != rank and ((rows[k] >> column) & 1):
                rows[k] ^= rows[rank]
        rank += 1
    return rank


all_pair_signs = []
for i, j in itertools.combinations(range(8), 2):
    row = [0] * 8
    row[i] = row[j] = 1
    all_pair_signs.append(row)
require(gf2_rank(all_pair_signs) == 7, "pair signs span the full even sign kernel")
WD8_ORDER = (2 ** 7) * factorial(8)
require(WD8_ORDER == 5_160_960, "W(D8) order")


def orbit(seed, actions):
    seen = {seed}
    todo = deque([seed])
    while todo:
        item = todo.popleft()
        for action in actions:
            nxt = action(item)
            if nxt not in seen:
                seen.add(nxt)
                todo.append(nxt)
    return frozenset(seen)


wd8_actions = [lambda x, matrix=matrix: matrix_action(matrix, x) for matrix in WD8_GENERATORS]
int_orbit = orbit(next(iter(INTEGER_ROOTS)), wd8_actions)
half_orbit = orbit(next(iter(HALF_ROOTS)), wd8_actions)
require(int_orbit == INTEGER_ROOTS and half_orbit == HALF_ROOTS, "two W(D8) root orbits")


def a2_subsystems(root_shell):
    shell = frozenset(root_shell)
    out = set()
    roots = sorted(shell)
    for alpha, beta in itertools.combinations(roots, 2):
        if dot(alpha, beta) != -4:
            continue
        gamma = add(alpha, beta)
        if gamma in shell:
            out.add(frozenset((alpha, neg(alpha), beta, neg(beta), gamma, neg(gamma))))
    return frozenset(out)


ALL_A2 = a2_subsystems(ROOTS)
INTEGER_A2 = a2_subsystems(INTEGER_ROOTS)
require((len(ALL_A2), len(INTEGER_A2)) == (1120, 224), "all and integer-normal A2 counts")


def subsystem_support(subsystem):
    return frozenset().union(*(support(root) for root in subsystem))


FRAMES = {}
for cut in INTEGER_A2:
    FRAMES.setdefault(subsystem_support(cut), set()).add(cut)
require(len(FRAMES) == 56, "56 coordinate-three-block frames")
require(all(len(key) == 3 and len(value) == 4 for key, value in FRAMES.items()), "four cuts per frame")


def half_projector_support(a2):
    return frozenset(root for root in HALF_ROOTS if all(dot(root, normal) == 0 for normal in a2))


PROJECTORS = {a2: half_projector_support(a2) for a2 in INTEGER_A2}
require(len(set(PROJECTORS.values())) == 224, "224 distinct half-root projectors")
require(all(len(block) == 32 for block in PROJECTORS.values()), "rank32 projectors")
for frame in FRAMES.values():
    blocks = [PROJECTORS[a2] for a2 in frame]
    require(set().union(*blocks) == set(HALF_ROOTS), "frame projectors sum to I128")
    require(
        sum(len(left & right) for left, right in itertools.combinations(blocks, 2)) == 0,
        "frame projectors are pairwise orthogonal",
    )
require(Fraction(128, 128 - 32) == Fraction(4, 3), "128/96 projector trace ratio")


def map_subsystem(matrix, subsystem):
    return frozenset(matrix_action(matrix, root) for root in subsystem)


for matrix in WD8_GENERATORS:
    for a2 in INTEGER_A2:
        target = map_subsystem(matrix, a2)
        require(target in INTEGER_A2, "W(D8) preserves integer A2 cuts")
        require(
            frozenset(matrix_action(matrix, root) for root in PROJECTORS[a2]) == PROJECTORS[target],
            "projector covariance",
        )

cut_actions = [lambda a2, matrix=matrix: map_subsystem(matrix, a2) for matrix in WD8_GENERATORS]
require(orbit(next(iter(INTEGER_A2)), cut_actions) == INTEGER_A2, "one 224-cut W(D8) orbit")

ORIGINAL_FRAME_KEY = frozenset((5, 6, 7))
original_four = FRAMES[ORIGINAL_FRAME_KEY]
j_four = {map_subsystem(J_MATRIX, a2) for a2 in original_four}
require(j_four == FRAMES[frozenset((4, 6, 7))], "J transports the original tetrahedral frame exactly")
require(not (original_four & j_four), "J does not preserve the selected original four-cut family")

# Reconstruct the sixty original Gaussian pair reflections and prove literal
# equality with s_alpha s_Jalpha on the full E8 shell.
root_index = {root: i for i, root in enumerate(sorted(ROOTS))}
root_list = tuple(sorted(ROOTS))
j_lines = []
used = set()
for alpha in root_list:
    if alpha in used:
        continue
    line = []
    current = alpha
    for _ in range(4):
        line.append(current)
        used.add(current)
        current = matrix_action(J_MATRIX, current)
    j_lines.append(tuple(line))
require(len(j_lines) == 60 and all(len(set(line)) == 4 for line in j_lines), "60 Gaussian root lines")

pair_permutations = []
for line in j_lines:
    alpha = line[0]
    images = tuple(root_index[complex_pair_reflection(alpha, root)] for root in root_list)
    require(all(images[images[i]] == i for i in range(240)), "pair reflection is involutive")
    pair_permutations.append(images)
require(len(set(pair_permutations)) == 60, "60 distinct original pair reflections")

witness_alpha = (-1, -1, -1, -1, -1, -1, -1, -1)
witness_integer = (-2, -2, 0, 0, 0, 0, 0, 0)
witness_half = complex_pair_reflection(witness_alpha, witness_integer)
require(witness_alpha in ROOTS, "witness Gaussian line is native")
require(witness_integer in INTEGER_ROOTS and witness_half in HALF_ROOTS, "G31 pair mixes D8 root sectors")
require(witness_half == (-1, -1, 1, 1, 1, 1, 1, 1), "exact mixing witness")

witness_action = lambda x: complex_pair_reflection(witness_alpha, x)
require(len(orbit(next(iter(INTEGER_ROOTS)), wd8_actions + [witness_action])) == 240, "joint action is root transitive")

# The W(E8) closure of one cut is all 1120 regular A2-normal E6 cuts.
e8_simple_actions = [lambda x, alpha=alpha: reflection(alpha, x) for alpha in SIMPLE_ROOTS]
e8_cut_actions = [
    lambda a2, alpha=alpha: frozenset(reflection(alpha, root) for root in a2)
    for alpha in SIMPLE_ROOTS
]
require(orbit(next(iter(ALL_A2)), e8_cut_actions) == ALL_A2, "one 1120-cut W(E8) orbit")


def e6_root_support(a2):
    return frozenset(root for root in ROOTS if all(dot(root, normal) == 0 for normal in a2))


FULL_PROJECTORS = {a2: e6_root_support(a2) for a2 in ALL_A2}
require(all(len(block) == 72 for block in FULL_PROJECTORS.values()), "rank72 E6 root-current projectors")

# Under full W(E8), coordinate frames close to every regular A3=D3 subsystem.
seed_d3 = frozenset(root for root in INTEGER_ROOTS if support(root) <= ORIGINAL_FRAME_KEY)
require(len(seed_d3) == 12, "seed D3 frame")
e8_frame_actions = [
    lambda frame, alpha=alpha: frozenset(reflection(alpha, root) for root in frame)
    for alpha in SIMPLE_ROOTS
]
ALL_A3_FRAMES = orbit(seed_d3, e8_frame_actions)
require(len(ALL_A3_FRAMES) == 7560, "7560 W(E8)-transported A3 frames")
incidence = Counter()
for frame in ALL_A3_FRAMES:
    local_a2 = a2_subsystems(frame)
    require(len(local_a2) == 4, "four A2 faces per A3 tetrahedron")
    incidence.update(local_a2)
require(set(incidence) == set(ALL_A2), "all A2 cuts occur in transported frames")
require(set(incidence.values()) == {27}, "each E6 cut lies in 27 transported A3 frames")

WE8_ORDER = 696_729_600
require(WE8_ORDER // WD8_ORDER == 135, "W(D8) index in W(E8)")

print("MATH.SEAM.ESERIES.JOINT.05 PASS")
print("native_sha256", " ".join(f"{path}={digest}" for path, digest in sorted(PINS.items())))
print("<W(D5)xW(D3),J> = W(D8), order", WD8_ORDER)
print("intermediate cuts: 224 integer-normal E6 cuts = 56 coordinate A3 frames x 4")
print("each frame: four orthogonal rank32 projectors sum to I128; complement rank96; 128/96=4/3")
print("J transports original frame {5,6,7} to {4,6,7}; fixed four-cut selection is not preserved")
print("G31 mixing witness alpha", witness_alpha, "maps", witness_integer, "to", witness_half)
print("<W(D5)xW(D3),J,original G31 pair reflections> = W(E8), order", WE8_ORDER)
print("full closure: 1120 E6 cuts, rank72 root projectors (rank78 on the full adjoint)")
print("full frame closure: 7560 A3 tetrahedra; four A2 faces each; 27 frames through each cut")
print("scope: transition group on transported flags; no fixed compiler mark, chiral64 block, or physical cut is selected")
