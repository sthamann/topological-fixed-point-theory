#!/usr/bin/env python3
"""Exact E7/E6 cut dictionary in the fixed native D5 + D3 chart.

This is an exploratory pure-mathematics checker.  It extracts only the two
native function definitions that own the E8 shell and Gaussian quarter clock;
their module-level suites are not executed.
"""

import ast
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import itertools
from math import gcd
from pathlib import Path

import numpy as np
import sympy as sp


REPO = Path("/Users/stefanhamann/Projekte/tfpt-theoryv4")
PINS = {}


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


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def independent_roots():
    """Independent reconstruction in the same doubled coordinate chart."""
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
    return sorted(out)


roots = np.array(native("verification/v1_e8_glue.py", "e8_roots")())
R = np.rint(2 * roots).astype(np.int64)
require(np.array_equal(2 * roots, R), "native roots are not half integral")
require(sorted(map(tuple, R)) == independent_roots(), "independent E8 reconstruction")
lookup = {tuple(v): i for i, v in enumerate(R)}
require(len(lookup) == 240, "E8 shell size")
G = R @ R.T  # four times the actual Euclidean inner product

Jvec = native("verification/v634_st31_structure.py", "J_vec")
J = np.column_stack([Jvec(v) for v in np.eye(8, dtype=np.int64)])
J0 = np.zeros((8, 8), dtype=np.int64)
for k in range(0, 8, 2):
    J0[k, k + 1], J0[k + 1, k] = -1, 1
require(np.array_equal(J, J0), "independent J reconstruction")
require(np.array_equal(J @ J, -np.eye(8, dtype=np.int64)), "J squared")
jp = np.array([lookup[tuple(J @ v)] for v in R], dtype=np.int64)

D5 = frozenset(np.flatnonzero(np.all(R[:, 5:] == 0, axis=1)))
D3 = frozenset(np.flatnonzero(np.all(R[:, :5] == 0, axis=1)))
D4_FIRST4 = frozenset(np.flatnonzero(np.all(R[:, 4:] == 0, axis=1)))
require((len(D5), len(D3), len(D4_FIRST4)) == (40, 12, 24), "fixed root blocks")


def sector(i):
    if i in D5:
        return "D5"
    if i in D3:
        return "D3"
    if np.all(np.abs(R[i]) == 1):
        return "half"
    return "mixed"


def root_type(indices):
    """ADE type of a closed simply-laced root subsystem, component by component."""
    unseen = set(indices)
    if not unseen:
        return "-"
    components = []
    while unseen:
        seed = unseen.pop()
        component = {seed}
        frontier = [seed]
        while frontier:
            u = frontier.pop()
            neighbours = {v for v in unseen if G[u, v] != 0}
            unseen.difference_update(neighbours)
            component.update(neighbours)
            frontier.extend(neighbours)
        rank = sp.Matrix(R[sorted(component)].tolist()).rank()
        count = len(component)
        if (rank, count) == (1, 2):
            name = "A1"
        elif count == rank * (rank + 1):
            name = f"A{rank}"
        elif rank >= 4 and count == 2 * rank * (rank - 1):
            name = f"D{rank}"
        elif (rank, count) in ((6, 72), (7, 126), (8, 240)):
            name = f"E{rank}"
        else:
            raise RuntimeError(f"unrecognized root component rank={rank}, roots={count}")
        components.append((rank, name))
    return "+".join(name for _, name in sorted(components, reverse=True))


def root_reflection_permutation(alpha):
    """Permutation of the E8 shell by reflection in doubled root alpha."""
    alpha = np.asarray(alpha, dtype=np.int64)
    dots = R @ alpha
    require(np.all(dots % 4 == 0), "root reflection integrality")
    images = R - (dots // 4)[:, None] * alpha[None, :]
    return np.array([lookup[tuple(v)] for v in images], dtype=np.int64)


def orbit(seed, actions):
    seen = {seed}
    frontier = [seed]
    while frontier:
        x = frontier.pop()
        for action in actions:
            y = int(action[x])
            if y not in seen:
                seen.add(y)
                frontier.append(y)
    return seen


# Simple reflections generate exactly W(D5) x W(D3)=1920 x 24.
e = np.eye(8, dtype=np.int64)
simple = [
    e[0] - e[1],
    e[1] - e[2],
    e[2] - e[3],
    e[3] - e[4],
    e[3] + e[4],
    e[5] - e[6],
    e[6] - e[7],
    e[6] + e[7],
]
root_actions = [root_reflection_permutation(2 * alpha) for alpha in simple]
family_root_actions = root_actions[5:]
H_ORDER = 1920 * 24


# ---------------------------------------------------------------------------
# E7 root-normal cuts.

root_lines = sorted(
    {tuple(sorted((i, lookup[tuple(-R[i])]))) for i in range(240)}
)
line_id = {root: k for k, line in enumerate(root_lines) for root in line}
line_actions = [
    np.array([line_id[int(p[line[0]])] for line in root_lines], dtype=np.int64)
    for p in root_actions
]
family_line_actions = line_actions[5:]
jline = np.array([line_id[int(jp[line[0]])] for line in root_lines], dtype=np.int64)
E7 = [frozenset(np.flatnonzero(G[line[0]] == 0)) for line in root_lines]

e7_rows = Counter()
for i, line in enumerate(root_lines):
    normal_sector = sector(line[0])
    e7_rows[(normal_sector, root_type(E7[i] & D5), root_type(E7[i] & D3))] += 1

expected_e7 = {
    ("D5", "A3+A1", "A3"): 20,
    ("D3", "D5", "A1"): 6,
    ("mixed", "D4", "A1+A1"): 30,
    ("half", "A4", "A2"): 64,
}
require(dict(e7_rows) == expected_e7, "E7 split dictionary")
unseen_lines = set(range(120))
e7_h_orbits = []
while unseen_lines:
    current = orbit(min(unseen_lines), line_actions)
    e7_h_orbits.append(current)
    unseen_lines.difference_update(current)
require(
    sorted(len(current) for current in e7_h_orbits) == [6, 20, 30, 64],
    "four complete H orbits of E7 cuts",
)
require(
    all(root_type(E7[i] & E7[int(jline[i])]) == "D6" for i in range(120)),
    "all J-paired E7 intersections are D6",
)
j_e7_pairs = {tuple(sorted((i, int(jline[i])))) for i in range(120)}
require(len(j_e7_pairs) == 60, "60 actual J pairs of E7 cuts")
e7_containing_d5 = {i for i in range(120) if D5 <= E7[i]}
require(len(e7_containing_d5) == 6, "six E7 cuts contain fixed D5")
require(
    orbit(min(e7_containing_d5), family_line_actions) == e7_containing_d5,
    "family Weyl group transitive on six E7 cuts containing D5",
)
e7_internal_j = {
    tuple(sorted((i, int(jline[i]))))
    for i in e7_containing_d5
    if int(jline[i]) in e7_containing_d5
}
require(len(e7_internal_j) == 1, "one internal J pair among six D5-compatible E7 cuts")


# ---------------------------------------------------------------------------
# All regular E6 cuts, represented by their orthogonal A2 root systems.

A2 = set()
for a, b in zip(*np.where(G == -4)):
    s = R[a] + R[b]
    A2.add(
        tuple(
            sorted(
                (
                    int(a),
                    lookup[tuple(-R[a])],
                    int(b),
                    lookup[tuple(-R[b])],
                    lookup[tuple(s)],
                    lookup[tuple(-s)],
                )
            )
        )
    )
A2 = sorted(A2)
a2_id = {subsystem: i for i, subsystem in enumerate(A2)}
E6 = [
    frozenset(np.flatnonzero(np.all(G[list(subsystem)] == 0, axis=0)))
    for subsystem in A2
]
require(len(A2) == len(set(E6)) == 1120, "all 1120 regular E6 cuts")
require(all(len(cut) == 72 for cut in E6), "every complement is E6")

a2_actions = [
    np.array(
        [a2_id[tuple(sorted(int(x) for x in p[list(subsystem)]))] for subsystem in A2],
        dtype=np.int64,
    )
    for p in root_actions
]
family_a2_actions = a2_actions[5:]
ja2 = np.array(
    [a2_id[tuple(sorted(int(jp[x]) for x in subsystem))] for subsystem in A2],
    dtype=np.int64,
)

unseen = set(range(1120))
h_orbits = []
while unseen:
    current = orbit(min(unseen), a2_actions)
    h_orbits.append(current)
    unseen.difference_update(current)
require(len(h_orbits) == 7, "seven W(D5)xW(D3) orbits")


def normal_signature(subsystem):
    return tuple(sorted(Counter(sector(x) for x in subsystem).items()))


e6_table = {}
for current in h_orbits:
    representative = min(current)
    signature = normal_signature(A2[representative])
    require(all(normal_signature(A2[i]) == signature for i in current), "orbit signature")
    fixed_types = {
        (root_type(E6[i] & D5), root_type(E6[i] & D3)) for i in current
    }
    require(len(fixed_types) == 1, "fixed split type constant on H orbit")
    j_intersections = Counter(
        root_type(E6[i] & E6[int(ja2[i])]) for i in current
    )
    e6_table[signature] = {
        "size": len(current),
        "stabilizer": H_ORDER // len(current),
        "cut": next(iter(fixed_types)),
        "J": dict(j_intersections),
    }

expected_e6 = {
    (("D5", 6),): (40, 1152, ("A1+A1", "A3"), {"A2+A2": 16, "D4": 24}),
    (("D5", 2), ("mixed", 4)): (120, 384, ("A3", "A1+A1"), {"A2+A2": 80, "D4": 40}),
    (("D5", 2), ("half", 4)): (320, 144, ("A2+A1", "A2"), {"A2+A2": 192, "D4": 128}),
    (("D3", 2), ("mixed", 4)): (60, 768, ("D4", "-"), {"A2+A2": 32, "D4": 28}),
    (("half", 4), ("mixed", 2)): (480, 96, ("A3", "A1"), {"A2+A2": 256, "D4": 224}),
    (("D3", 6),): (4, 11520, ("D5", "-"), {"D4": 4}),
    (("D3", 2), ("half", 4)): (96, 480, ("A4", "A1"), {"A2+A2": 64, "D4": 32}),
}
require(set(e6_table) == set(expected_e6), "E6 orbit signatures")
for signature, expected in expected_e6.items():
    row = e6_table[signature]
    got = (row["size"], row["stabilizer"], row["cut"], row["J"])
    require(got == expected, f"E6 orbit row {signature}: {got}")

j_e6_pairs = {tuple(sorted((i, int(ja2[i])))) for i in range(1120)}
require(len(j_e6_pairs) == 560 and all(i != int(ja2[i]) for i in range(1120)), "560 E6 J pairs")
j_types = Counter(root_type(E6[i] & E6[int(ja2[i])]) for i in range(1120))
require(j_types == {"A2+A2": 640, "D4": 480}, "complete E6 J-intersection types")

d5_compatible = {i for i, cut in enumerate(E6) if D5 <= cut}
require(len(d5_compatible) == 4, "exactly four E6 cuts contain D5")
require(
    orbit(min(d5_compatible), family_a2_actions) == d5_compatible,
    "the four D5-compatible E6 cuts are one W(D3) orbit",
)


def primitive_axis(subsystem):
    vectors = [R[i, 5:] // 2 for i in subsystem]
    a = next(v for v in vectors if np.any(v))
    b = next(v for v in vectors if np.linalg.matrix_rank(np.array([a, v])) == 2)
    n = np.cross(a, b)
    divisor = gcd(gcd(abs(int(n[0])), abs(int(n[1]))), abs(int(n[2])))
    n //= divisor
    first = next(int(x) for x in n if x)
    if first < 0:
        n = -n
    return tuple(int(x) for x in n)


axes = {primitive_axis(A2[i]) for i in d5_compatible}
require(axes == {(1, 1, 1), (1, 1, -1), (1, -1, 1), (1, -1, -1)}, "tetrahedral axes")
require(
    all(int(ja2[i]) not in d5_compatible for i in d5_compatible),
    "J sends every D5-compatible E6 cut outside that family",
)
require(
    all(E6[i] & E6[int(ja2[i])] == D4_FIRST4 for i in d5_compatible),
    "the four compatible J pairs share the same first-four-coordinate D4",
)


# ---------------------------------------------------------------------------
# Exact representation/character checks in the fixed chart.

integer_mixed = {i for i in range(240) if sector(i) == "mixed"}
half = {i for i in range(240) if sector(i) == "half"}
require((len(D5), len(D3), len(integer_mixed), len(half)) == (40, 12, 60, 128), "root branching counts")

v10 = {
    tuple(v)
    for v in itertools.chain.from_iterable(
        ((tuple(e5), tuple(-e5)) for e5 in np.eye(5, dtype=np.int64))
    )
}
v6 = {
    tuple(v)
    for v in itertools.chain.from_iterable(
        ((tuple(e3), tuple(-e3)) for e3 in np.eye(3, dtype=np.int64))
    )
}
mixed_weights = {(tuple((R[i, :5] // 2)), tuple((R[i, 5:] // 2))) for i in integer_mixed}
require(mixed_weights == set(itertools.product(v10, v6)), "(10,6) mixed-root character")

s5_even, s5_odd, s3_even, s3_odd = set(), set(), set(), set()
for signs in itertools.product((-1, 1), repeat=5):
    (s5_even if signs.count(-1) % 2 == 0 else s5_odd).add(signs)
for signs in itertools.product((-1, 1), repeat=3):
    (s3_even if signs.count(-1) % 2 == 0 else s3_odd).add(signs)
half_weights = {(tuple(R[i, :5]), tuple(R[i, 5:])) for i in half}
require(
    half_weights == set(itertools.product(s5_even, s3_even)) | set(itertools.product(s5_odd, s3_odd)),
    "(16,4)+(16bar,4bar) half-root character",
)

# Use the +++ tetrahedral axis to verify the E6 -> D5 x U(1) character.
chosen = next(i for i in d5_compatible if primitive_axis(A2[i]) == (1, 1, 1))
n = np.array((1, 1, 1), dtype=np.int64)


def family_charge(i):
    # R is doubled, so this is 2*n.(actual family weight).
    return int(n @ R[i, 5:])


chosen_root_sectors = Counter((sector(i), family_charge(i)) for i in E6[chosen])
require(
    chosen_root_sectors == {("D5", 0): 40, ("half", 3): 16, ("half", -3): 16},
    "E6 adjoint roots restrict as D5 roots plus the two chiral spinors",
)
remaining = set(range(240)) - set(E6[chosen]) - set(A2[chosen])
remaining_sectors = Counter((sector(i), family_charge(i)) for i in remaining)
require(
    remaining_sectors
    == {
        ("D3", 4): 3,
        ("D3", -4): 3,
        ("mixed", 2): 30,
        ("mixed", -2): 30,
        ("half", 1): 48,
        ("half", -1): 48,
    },
    "E8 -> E6 x A2 character sectors",
)

# Each D5-compatible cut gives an actual rank-16 projector in one native
# (16,4) half-root block and an actual rank-48 (16,3) complement.
T4 = {
    signs
    for signs in itertools.product((-1, 1), repeat=3)
    if signs.count(-1) % 2 == 0
}
B16x4 = {
    i
    for i in half
    if tuple(R[i, :5]) in s5_even and tuple(R[i, 5:]) in T4
}
require(len(B16x4) == 64, "native (16,4) half-root block")


def even_vertex(cut):
    axis = primitive_axis(A2[cut])
    if axis.count(-1) % 2:
        axis = tuple(-x for x in axis)
    require(axis in T4, "oriented family vertex")
    return axis


def family_step_edges(weights, a2_roots):
    weights = {tuple(w) for w in weights}
    roots3 = [tuple(R[a, 5:]) for a in a2_roots]
    return {
        (w, tuple(w[k] + root[k] for k in range(3)))
        for w in weights
        for root in roots3
        if tuple(w[k] + root[k] for k in range(3)) in weights
    }


projector_subsets = {}
off48_subsets = {}
branch81_subsets = {}
for cut in d5_compatible:
    t = even_vertex(cut)
    p16 = {i for i in B16x4 if tuple(R[i, 5:]) == t}
    q48 = B16x4 - p16
    require(len(p16) == 16 and len(q48) == 48, "16+48 projector ranks")
    require(p16 == B16x4 & E6[cut], "rank-16 block is exactly the E6-adjoint spinor subset")
    require(
        {(tuple(R[i, :5]), tuple(R[i, 5:])) for i in q48}
        == set(itertools.product(s5_even, T4 - {t})),
        "rank-48 block is exactly the (16,3) weight product",
    )

    # The cut-oriented U(1) charge makes 4=3_1+1_-3 literally.
    qcut = lambda i: -int(np.array(t, dtype=np.int64) @ R[i, 5:])
    require({qcut(i) for i in p16} == {-3} and {qcut(i) for i in q48} == {1}, "4 -> 3_1 + 1_-3")
    mixed30 = {i for i in integer_mixed if qcut(i) == -2}
    scalar3 = {i for i in D3 if qcut(i) == 4}
    require((len(q48), len(mixed30), len(scalar3)) == (48, 30, 3), "81=48+30+3")

    half_family_weights = {tuple(R[i, 5:]) for i in q48}
    vector_family_weights = {tuple(R[i, 5:]) for i in mixed30}
    scalar_family_weights = {tuple(R[i, 5:]) for i in scalar3}
    require(
        all(len(weights) == 3 for weights in (half_family_weights, vector_family_weights, scalar_family_weights)),
        "three literal A2 triplet weight sets",
    )
    require(
        all(len(family_step_edges(weights, A2[cut])) == 6 for weights in (half_family_weights, vector_family_weights, scalar_family_weights)),
        "actual A2 root operators act as a triplet on all three factors",
    )

    # Actual D5 and A2 root additions act on the separate tensor factors.
    d5_edges = set()
    a2_edges = set()
    for source in q48:
        for root in D5:
            target = lookup.get(tuple(R[source] + R[root]))
            if target in q48:
                require(tuple(R[target, 5:]) == tuple(R[source, 5:]), "D5 preserves family weight")
                d5_edges.add((source, target))
        for root in A2[cut]:
            target = lookup.get(tuple(R[source] + R[root]))
            if target in q48:
                require(tuple(R[target, :5]) == tuple(R[source, :5]), "A2 preserves D5 spinor weight")
                a2_edges.add((source, target))
    require(len(d5_edges) == 480 and len(a2_edges) == 96, "literal (16,3) root-step action")
    require(
        not any(
            lookup.get(tuple(R[source] + R[root])) in B16x4
            for source in p16
            for root in A2[cut]
        ),
        "A2 annihilates the family singlet within the (16,4) block",
    )

    require(
        all(Fraction(sum(int(x) ** 2 for x in R[i, :5]), 4) == Fraction(5, 4) for i in B16x4),
        "D5 projection norm 5/4 on the native spinor block",
    )
    require(
        all(Fraction(sum(int(x) ** 2 for x in R[i, 5:]), 4) == Fraction(3, 4) for i in B16x4),
        "A3 projection norm 3/4 on the native spinor block",
    )
    require(Fraction(len(B16x4), len(q48)) == Fraction(4, 3), "cut-independent trace ratio 64/48")

    projector_subsets[cut] = p16
    off48_subsets[cut] = q48
    branch81_subsets[cut] = q48 | mixed30 | scalar3
    require(len(branch81_subsets[cut]) == 81, "actual (27,3) off-block root subset")

require(set().union(*projector_subsets.values()) == B16x4, "four rank-16 projectors resolve the identity")
require(
    sum(len(projector_subsets[a] & projector_subsets[b]) for a, b in itertools.combinations(d5_compatible, 2)) == 0,
    "four rank-16 projectors are pairwise orthogonal on the root-current basis",
)
for action, cut_action in zip(family_root_actions, family_a2_actions):
    for cut in d5_compatible:
        target_cut = int(cut_action[cut])
        require(target_cut in d5_compatible, "family action preserves the four-cut set")
        require(
            {int(action[root]) for root in projector_subsets[cut]} == projector_subsets[target_cut],
            "S4 conjugates the four projectors canonically",
        )

c3 = Fraction(1, 8)
phi_base = Fraction(1, 6)
require(phi_base / c3 == Fraction(4, 3) == Fraction(64, 48), "existing phi_base/c3 projector ratio")

Y = tuple([Fraction(-1, 3)] * 3 + [Fraction(1, 2)] * 2 + [Fraction(0)] * 3)
require(sum(y * y for y in Y) == Fraction(5, 6), "original 3+2 mark norm")
require(
    all(
        sum(Y[k] * Fraction(int(R[root, k]), 2) for k in range(8)) == 0
        for cut in d5_compatible
        for root in A2[cut]
    ),
    "Y lies in every compatible E6 cut",
)

# A D3-normal E7 cut realizes the full D5 x A1 x U(1) character.
e7_normal = lookup[(0, 0, 0, 0, 0, 2, -2, 0)]
e7_character_cut = frozenset(np.flatnonzero(G[e7_normal] == 0))
e7_character_counts = Counter((sector(i), int(R[i, 7])) for i in e7_character_cut)
require(
    e7_character_counts
    == {
        ("D5", 0): 40,
        ("D3", 0): 2,
        ("mixed", 2): 10,
        ("mixed", -2): 10,
        ("half", 1): 32,
        ("half", -1): 32,
    },
    "E7 adjoint root character under D5 x A1 x U(1)",
)
e7_half_plus = {
    (tuple(R[i, :5]), int((R[i, 5] + R[i, 6]) // 2))
    for i in e7_character_cut
    if sector(i) == "half" and R[i, 7] == 1
}
e7_half_minus = {
    (tuple(R[i, :5]), int((R[i, 5] + R[i, 6]) // 2))
    for i in e7_character_cut
    if sector(i) == "half" and R[i, 7] == -1
}
require(e7_half_plus == set(itertools.product(s5_even, {-1, 1})), "(16,2)_1 E7 weights")
require(e7_half_minus == set(itertools.product(s5_odd, {-1, 1})), "(16bar,2)_-1 E7 weights")


def y_charge_from_doubled(v):
    return sum(Y[k] * Fraction(int(v[k]), 2) for k in range(8))


# The canonical E7 pair of the preceding contract and its missing 48 roots.
alpha = lookup[(2, 2, 0, 0, 0, 0, 0, 0)]
beta = int(jp[alpha])
pair_union = set(np.flatnonzero((G[alpha] == 0) | (G[beta] == 0)))
missing48 = set(range(240)) - pair_union
require(len(pair_union) == 192 and len(missing48) == 48, "E7-pair 192+48 split")
d6_vector = {
    tuple(v)
    for v in itertools.chain.from_iterable(
        ((tuple(unit), tuple(-unit)) for unit in np.eye(6, dtype=np.int64))
    )
}
missing_character = {
    (tuple(R[i, 2:] // 2), int(G[i, alpha] // 4), int(G[i, beta] // 4))
    for i in missing48
}
require(
    missing_character == set(itertools.product(d6_vector, {-1, 1}, {-1, 1})),
    "missing48 is exactly the (12,2,2) weight product",
)
missing_y = Counter(y_charge_from_doubled(R[i]) for i in missing48)
expected_missing_y = {
    Fraction(-5, 6): 4,
    Fraction(-2, 3): 2,
    Fraction(-1, 3): 12,
    Fraction(-1, 6): 4,
    Fraction(0): 4,
    Fraction(1, 6): 4,
    Fraction(1, 3): 12,
    Fraction(2, 3): 2,
    Fraction(5, 6): 4,
}
require(missing_y == expected_missing_y, "missing48 Y character")


def spinor_y_hist(parity):
    out = Counter()
    for signs in itertools.product((-1, 1), repeat=5):
        if signs.count(-1) % 2 == parity:
            charge = sum(Fraction(signs[k], 2) * Y[k] for k in range(5))
            out[charge] += 3  # the tempting three copies
    return out


three_spinor_even = spinor_y_hist(0)
three_spinor_odd = spinor_y_hist(1)
require(missing_y != three_spinor_even and missing_y != three_spinor_odd, "48 is not 3 times a D5 spinor character")


def fmt_signature(signature):
    return "+".join(f"{count}{name}" for name, count in signature)


def fmt_counter(counter):
    return ",".join(f"{key}:{counter[key]}" for key in sorted(counter))


print("MATH.SEAM.ESERIES.CUTS.04 PASS")
print("native_sha256", " ".join(f"{path}={digest}" for path, digest in sorted(PINS.items())))
print("E7 root-normal cuts: 120 = 20[D5] + 6[D3] + 30[mixed] + 64[half]")
print("E7 J pairs: 60; every intersection D6; D5-containing cuts: 6 (one S4 orbit, one internal J pair)")
print("E6 regular cuts: 1120; H=W(D5)xW(D3) orbits:")
for signature, expected in expected_e6.items():
    size, stabilizer, cut, jcounts = expected
    print(
        f"  normal={fmt_signature(signature):12s} size={size:3d} stab={stabilizer:5d} "
        f"cut=({cut[0]},{cut[1]}) J={fmt_counter(jcounts)}"
    )
print("E6 J pairs: 560; ordered intersections A2+A2:640, D4:480")
print("D5-compatible E6 cuts: 4; one family-S4 tetrahedron; axes", sorted(axes))
print("their J partners leave D5; all four pair intersections are the fixed first-four-coordinate D4")
print("character E8|D5xA3=(45,1)+(1,15)+(10,6)+(16,4)+(16bar,4bar)")
print("character E6|D5xU1=45_0+1_0+16_-3+16bar_3 (up to simultaneous orientation reversal)")
print("each cut: (16,4)=16_singlet + (16,3), projector ranks16+48; (27,3)=48+30+3")
print("canonical S4 projectors: pairwise orthogonal, sum I64, literal D5/A2 root-step action480/96")
print("Tr(I64)/Tr(Q48)=64/48=4/3=phi_base/c3; every spinor weight has D5 norm5/4 and A3 norm3/4")
print("Y=(-1/3,-1/3,-1/3,1/2,1/2,0,0,0), Y^2=5/6; Y fixes the unordered four-cut family")
print("missing48 character differs from both 3x16 chiralities; actual type is (12,2,2) under D6xA1xA1")
