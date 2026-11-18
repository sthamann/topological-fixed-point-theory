#!/usr/bin/env python3
"""Exact joint C12-cover and full E8 lattice-field clock certificate.

The script keeps three roles separate:

* Delta_12 is deck transport on the principal marked cover;
* Q_12 covers the base quarter rotation x -> i x;
* character-corrected lifts act on the existing lattice VOA V_E8.

No section of the full G31 lattice-lift extension is constructed here.
"""

import ast
from functools import lru_cache
from hashlib import sha256
import itertools
import json
from pathlib import Path

import numpy as np
import sympy as S


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PINS = {}


def need(condition, label):
    if not bool(condition):
        raise RuntimeError(label)


def functions(relative, names, namespace):
    path = REPO / relative
    raw = path.read_bytes()
    PINS[relative] = sha256(raw).hexdigest()
    nodes = [node for node in ast.parse(raw).body
             if isinstance(node, ast.FunctionDef) and node.name in names]
    need({node.name for node in nodes} == set(names), "native definitions: " + relative)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), "exec"), namespace)
    return namespace


def pin(relative):
    raw = (REPO / relative).read_bytes()
    PINS[relative] = sha256(raw).hexdigest()


for predecessor in (
    "experiments/theory-contracts/e-series-compiler-closure-20261023/COVERS.txt",
    "experiments/theory-contracts/e-series-compiler-closure-20261023/FIELDS.txt",
    "experiments/theory-contracts/e-series-compiler-closure-20261023/cover_dictionary.py",
    "experiments/theory-contracts/e-series-compiler-closure-20261023/field_lift_dictionary.py",
    "experiments/theory-contracts/marked-seam-monodromy-20261023/THEORY.txt",
    "experiments/theory-contracts/marked-seam-monodromy-20261023/checker.py",
    "experiments/theory-contracts/global-g31-lattice-lift-20261023/THEOREMS.txt",
):
    pin(predecessor)


def parity_mask(coordinates):
    return sum((int(value) & 1) << i for i, value in enumerate(coordinates))


def char_sign(mask, coordinates):
    return -1 if (mask & parity_mask(coordinates)).bit_count() & 1 else 1


def dual_action(matrix, mask):
    """Character a -> a composed matrix in the fixed integral basis."""
    bits = S.Matrix([(mask >> i) & 1 for i in range(8)])
    out = matrix.T * bits
    return sum((int(out[i]) & 1) << i for i in range(8))


def matrix_key(matrix):
    return tuple(int(x) for x in matrix)


def gf2_rank(columns):
    pivots = {}
    for value in columns:
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    return len(pivots)


# The actual E8 lattice basis and root shell used by the predecessor.
native = functions(
    "verification/v1_e8_glue.py",
    {"e8_roots", "e8_simple_roots"},
    {"np": np, "itertools": itertools},
)
roots = np.rint(2 * np.array(native["e8_roots"]())).astype(np.int64)
B = S.Matrix(np.rint(2 * native["e8_simple_roots"]()).astype(np.int64)).T
Bi = B.inv()
coordinates = np.array((Bi * S.Matrix(roots.T)).tolist(), dtype=np.int64).T
root_lookup = {tuple(v): i for i, v in enumerate(roots)}
basis_root_indices = [root_lookup[tuple(int(x) for x in B[:, j])] for j in range(8)]
basis_vectors = tuple(tuple(int(x) for x in B[:, j]) for j in range(8))


# The unchanged ordered lattice cocycle.  Doubled coordinates are used so
# half-integral E8 vectors remain exact integers.
old = functions(
    "experiments/theory-contracts/seam-joint-ward-state-20261014/operator_transport_check.py",
    {"eps_car", "split_even_lattice", "eps_e8"},
    {"S": S, "q": S.ones(8, 1), "s": S.ones(8, 1) / 2},
)


@lru_cache(None)
def eps(u, v):
    c, h = u[0] & 1, v[0] & 1
    d = [(x - c) // 2 for x in u]
    f = [(x - h) // 2 for x in v]
    need(sum(d) % 2 == 0 and sum(f) % 2 == 0, "even E8 lattice")
    exponent = c * (sum(f) // 2) + c * h
    exponent += sum(d[i] * f[j] for i in range(8) for j in range(i))
    exponent += c * h * sum(i * (d[i] + f[i]) for i in range(8))
    return -1 if exponent & 1 else 1


for i in range(8):
    for j in range(8):
        need(
            eps(basis_vectors[i], basis_vectors[j])
            == int(old["eps_e8"](S.Matrix(basis_vectors[i]) / 2,
                                 S.Matrix(basis_vectors[j]) / 2)),
            "literal native cocycle",
        )


@lru_cache(None)
def ordered(coordinate_tuple, image_tuple):
    point = (0,) * 8
    phase = 1
    for vector, power in zip(image_tuple, coordinate_tuple):
        unit = eps(vector, tuple(-x for x in vector)) if power < 0 else 1
        step = tuple(-x for x in vector) if power < 0 else vector
        for _ in range(abs(power)):
            phase *= unit * eps(point, step)
            point = tuple(x + y for x, y in zip(point, step))
    return phase


@lru_cache(None)
def section_phase(matrix_tuple, coordinate_tuple):
    matrix = S.Matrix(8, 8, matrix_tuple)
    images = tuple(tuple(int(x) for x in B * matrix[:, j]) for j in range(8))
    return (ordered(coordinate_tuple, images)
            * ordered(coordinate_tuple, basis_vectors))


def section_basis_phases(matrix):
    key = matrix_key(matrix)
    return tuple(section_phase(key, tuple(int(x) for x in S.eye(8)[:, i]))
                 for i in range(8))


def power_basis_phases(matrix, character, exponent):
    """Phases of (S_matrix D_character)^exponent on the eight generators."""
    key = matrix_key(matrix)
    output = []
    for i in range(8):
        coordinate = S.eye(8)[:, i]
        phase = 1
        for _ in range(exponent):
            point = tuple(int(x) for x in coordinate)
            phase *= section_phase(key, point) * char_sign(character, point)
            coordinate = matrix * coordinate
        output.append(int(phase))
    return tuple(output)


def power_phase_at(matrix, character, exponent, coordinate_tuple):
    key = matrix_key(matrix)
    coordinate = S.Matrix(coordinate_tuple)
    phase = 1
    for _ in range(exponent):
        point = tuple(int(x) for x in coordinate)
        phase *= section_phase(key, point) * char_sign(character, point)
        coordinate = matrix * coordinate
    return int(phase), tuple(int(x) for x in coordinate)


def relative_character(left_phases, right_phases):
    """Character mask taking right generator phases to left phases."""
    return sum((left_phases[i] != right_phases[i]) << i for i in range(8))


# Actual commuting compiler matrices in the original real chart.
J_standard = S.zeros(8)
for i in range(0, 8, 2):
    J_standard[i, i + 1] = -1
    J_standard[i + 1, i] = 1
sigma_indices = (4, 5, 0, 1, 2, 3, 6, 7)
sigma_standard = S.zeros(8)
for output, source in enumerate(sigma_indices):
    sigma_standard[output, source] = 1

J = Bi * J_standard * B
sigma = Bi * sigma_standard * B
c = J * sigma
need(all(x.q == 1 for x in J) and all(x.q == 1 for x in sigma), "integral actions")
need(J * sigma == sigma * J, "J and sigma commute")
need(J**4 == S.eye(8) and sigma**3 == S.eye(8), "factor orders")
need(c**4 == sigma and c**9 == J and c**12 == S.eye(8), "compiler powers")

J_phase = section_basis_phases(J)
sigma_phase = section_basis_phases(sigma)
c_phase = section_basis_phases(c)
need(power_basis_phases(J, 0, 4) == (1,) * 8, "fixed Jhat has order four")


# Classify every lift chat_a=S_c D_a whose ninth power is the fixed
# generator-normalized Jhat.  Equality on the eight generators proves equality
# on every lattice vector because both sides have the same lattice action.
solutions = []
for character in range(256):
    ninth = power_basis_phases(c, character, 9)
    if ninth != J_phase:
        continue
    twelfth = power_basis_phases(c, character, 12)
    need(twelfth == (1,) * 8, "selected chat has order twelve")
    fourth = power_basis_phases(c, character, 4)
    sigma_character = relative_character(fourth, sigma_phase)
    need(power_basis_phases(sigma, sigma_character, 3) == (1,) * 8,
         "induced sigmahat has order three")
    solutions.append({
        "c_input_character": character,
        "sigma_input_character_relative_to_generator_section": sigma_character,
    })

need(solutions, "existence of compiler field lift")
solution_masks = {row["c_input_character"] for row in solutions}

# The ninth-power condition is the affine equation N_9(a)=d.  The columns
# below are the exact F2 matrix of N_9 in the ordered dual basis.
unadjusted_ninth = power_basis_phases(c, 0, 9)
ninth_defect = relative_character(unadjusted_ninth, J_phase)
ninth_norm_columns = [
    relative_character(power_basis_phases(c, 1 << i, 9), unadjusted_ninth)
    for i in range(8)
]
need(gf2_rank(ninth_norm_columns) == 6, "ninth-power affine rank six")
need(len(solutions) == 1 << (8 - gf2_rank(ninth_norm_columns)),
     "complete ninth-power solution count")


# Character gauges that keep the fixed Jhat unchanged.  Their conjugation
# action on chat_a is a -> a+b+(b composed c).
J_stabilizer = [b for b in range(256)
                if b ^ dual_action(J, b) == 0]
c_stabilizer = [b for b in J_stabilizer
                if b ^ dual_action(c, b) == 0]
need(len(J_stabilizer) == 16, "fixed-Jhat gauge group")


def gauge_transform(character, gauge):
    return character ^ gauge ^ dual_action(c, gauge)


remaining = set(solution_masks)
orbits = []
while remaining:
    representative = min(remaining)
    orbit = {gauge_transform(representative, b) for b in J_stabilizer}
    need(orbit <= solution_masks, "fixed-Jhat gauge preserves solutions")
    orbits.append(sorted(orbit))
    remaining -= orbit

need(all(len(orbit) == len(J_stabilizer) // len(c_stabilizer) for orbit in orbits),
     "uniform character-gauge stabilizer")


# If the generator-normalized sigma section is also frozen, Bezout(4,9)=1
# permits at most one chat.  The explicit census decides whether it exists.
canonical_sigma_solutions = [row for row in solutions
                             if row["sigma_input_character_relative_to_generator_section"] == 0]
need(len(canonical_sigma_solutions) <= 1, "fixed Jhat and fixed sigmahat uniqueness")


def composed_basis_phases(left_matrix, left_character, right_matrix, right_character):
    left_key = matrix_key(left_matrix)
    right_key = matrix_key(right_matrix)
    phases = []
    for i in range(8):
        point = S.eye(8)[:, i]
        source = tuple(int(x) for x in point)
        middle = right_matrix * point
        middle_tuple = tuple(int(x) for x in middle)
        phase = section_phase(right_key, source) * char_sign(right_character, source)
        phase *= section_phase(left_key, middle_tuple) * char_sign(left_character, middle_tuple)
        phases.append(int(phase))
    return tuple(phases)


for row in solutions:
    character = row["c_input_character"]
    sigma_character = row["sigma_input_character_relative_to_generator_section"]
    need(
        composed_basis_phases(J, 0, sigma, sigma_character)
        == tuple(c_phase[i] * char_sign(character, S.eye(8)[:, i]) for i in range(8)),
        "chat=Jhat sigmahat on all lattice generators",
    )
    need(
        composed_basis_phases(sigma, sigma_character, J, 0)
        == composed_basis_phases(J, 0, sigma, sigma_character),
        "Jhat and induced sigmahat commute",
    )
    for coordinate in coordinates:
        point = tuple(int(x) for x in coordinate)
        ninth_phase, ninth_point = power_phase_at(c, character, 9, point)
        need(ninth_point == tuple(int(x) for x in J * S.Matrix(point)),
             "ninth-power root projection")
        need(ninth_phase == section_phase(matrix_key(J), point),
             "ninth-power phase on all 240 roots")
        fourth_phase, fourth_point = power_phase_at(c, character, 4, point)
        need(fourth_point == tuple(int(x) for x in sigma * S.Matrix(point)),
             "fourth-power root projection")
        need(fourth_phase
             == section_phase(matrix_key(sigma), point) * char_sign(sigma_character, point),
             "fourth-power phase on all 240 roots")
        twelfth_phase, twelfth_point = power_phase_at(c, character, 12, point)
        need(twelfth_point == point and twelfth_phase == 1,
             "order twelve on all 240 roots")


# Exact associated-bundle action.  C12 is used as a principal refinement of
# the peripheral C4 cover.  Its deck generator acts on V_L through Jhat, so
# Delta12^4 is in the fibre kernel.  R12=(Q12,Delta12^4).
def cover_power(pair, exponent):
    base, deck = pair
    return exponent * base % 4, exponent * deck % 12


R12 = (1, 4)
need(cover_power(R12, 4) == (0, 4), "R12 fourth power")
need(cover_power(R12, 9) == (1, 0), "R12 ninth power is plain base lift")
need(cover_power(R12, 12) == (0, 0), "R12 order twelve")

# The diagonal quotient is generated by (Delta12,Jhat^-1).  Since chat and
# Jhat are powers of the same order-twelve lift, the joint action commutes with
# this relation.  The descended powers are therefore vertical sigmahat and
# base-quarter-plus-fibre-Jhat respectively.
need(c**9 == J and c**4 == sigma, "field projections match cover powers")


result = {
    "contract": "MATH.SEAM.ESERIES.JOINTCLOCK.01",
    "verdict": "EXACT_JOINT_ACTION_AND_COMPLETE_FIXED_JHAT_CYCLIC_LIFT_CLASSIFICATION",
    "source_sha256": PINS,
    "cover": {
        "principal_refinement": "C12: s^12=x^4-1",
        "peripheral_representation": "Delta12 -> Jhat; factors through C4",
        "diagonal_relation": "(Delta12 p,v)~(p,Jhat v)",
        "compiler_cover_clock": "R12=Q12 Delta12^4",
        "R12_fourth": cover_power(R12, 4),
        "R12_ninth": cover_power(R12, 9),
    },
    "field_sections": {
        "Jhat_generator_phases": list(J_phase),
        "sigmahat_generator_phases": list(sigma_phase),
        "chat_generator_phases": list(c_phase),
    },
    "fixed_Jhat_lifts": {
        "count": len(solutions),
        "ninth_power_affine_defect": ninth_defect,
        "ninth_power_norm_columns": ninth_norm_columns,
        "ninth_power_norm_rank": gf2_rank(ninth_norm_columns),
        "solutions": solutions,
        "Jhat_preserving_character_gauges": J_stabilizer,
        "gauge_stabilizer": c_stabilizer,
        "gauge_orbits": orbits,
        "orbit_count": len(orbits),
    },
    "fixed_Jhat_and_generator_normalized_sigmahat": {
        "count": len(canonical_sigma_solutions),
        "solutions": canonical_sigma_solutions,
    },
    "descended_joint_relations": {
        "Chat_fourth": "[p,v] -> [p,sigmahat v]",
        "Chat_ninth": "[p,v] -> [Q12 p,Jhat v]",
        "peripheral_J": "[Delta12 p,v]=[p,Jhat v], base component 0",
        "compiler_J": "[Q12 p,Jhat v], base component 1",
    },
    "scope": (
        "Existing C12 marked cover, unchanged E8 cocycle, fixed generator-normalized Jhat, "
        "and all character lifts over the single cyclic compiler subgroup <c>. "
        "No section of the full nonsplit G31 extension and no carrier identification."
    ),
}

for relative, digest in PINS.items():
    need(sha256((REPO / relative).read_bytes()).hexdigest() == digest,
         "unchanged source: " + relative)

print(json.dumps({
    "verdict": result["verdict"],
    "ninth_power_affine_defect": ninth_defect,
    "ninth_power_norm_columns": ninth_norm_columns,
    "ninth_power_norm_rank": gf2_rank(ninth_norm_columns),
    "solution_masks": sorted(solution_masks),
    "sigma_masks": [row["sigma_input_character_relative_to_generator_section"]
                    for row in solutions],
    "fixed_Jhat_lifts": len(solutions),
    "fixed_Jhat_gauge_orbits": len(orbits),
    "fixed_Jhat_and_generator_normalized_sigmahat": len(canonical_sigma_solutions),
    "orbit_size": len(orbits[0]),
    "Jhat_preserving_gauges": len(J_stabilizer),
    "gauge_stabilizer": len(c_stabilizer),
}, sort_keys=True))
