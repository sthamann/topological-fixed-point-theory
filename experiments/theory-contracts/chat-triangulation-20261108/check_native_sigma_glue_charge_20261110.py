"""Exact marked clock census on the deployed E8 root lattice and its fold."""
import ast
from collections import Counter
from itertools import combinations, product
import json
from pathlib import Path
import sys

import sympy as sp

repo = Path(sys.argv[1]).resolve()
native = repo / "verification/v634_st31_structure.py"
tree = ast.parse(native.read_text())
index = next(ast.literal_eval(node.value) for node in tree.body
             if isinstance(node, ast.Assign)
             and any(isinstance(t, ast.Name) and t.id == "SIGMA_IDX" for t in node.targets))
fold = repo / "experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json"
m = sp.Matrix(json.loads(fold.read_text())["explicit_Weyl_matrix"])
eye = sp.eye(8)
sigma = eye[index, :]
r = m.T * sigma * m
p_d = sp.diag(1, 1, 1, 1, 1, 0, 0, 0)
ell = sp.Matrix([1, 1, 1, 1, 1, 0, 0, 0])
y6 = sp.Matrix([-2, -2, -2, 3, 3, 0, 0, 0])
y = y6 / 6
assert sigma**3 == eye and m.T * m == eye
assert sigma.T * ell - ell == eye[:, 5] - eye[:, 3]

# Doubled roots, to count every observation in exact integers.
roots = set()
for a, b in combinations(range(8), 2):
    for x, z in product((-2, 2), repeat=2):
        beta = [0] * 8
        beta[a], beta[b] = x, z
        roots.add(tuple(beta))
roots.update(b for b in product((-1, 1), repeat=8) if b.count(-1) % 2 == 0)
assert len(roots) == 240


def s(b):
    return tuple(b[i] for i in index)


def grade(b):
    return sum(b[:5]) % 4


def charge(b):
    return y.dot(sp.Matrix(b)) / 2


def source_sector(b):
    q = m.T * sp.Matrix(b) / 2
    if all(x.q == 1 for x in q):
        return "0" if int(sum(q)) % 2 == 0 else "v"
    return "s" if int(sum(q - sp.ones(8, 1)/2)) % 2 == 0 else "c"


def print_root(b):
    return list(map(str, sp.Matrix(b) / 2))


class_table = [[0] * 4 for _ in range(4)]
difference_counts = {0: Counter(), 2: Counter()}
change_type = Counter()
compensation_weights = Counter()
weight_one_compensation_charges = Counter()
family_unchanged = {}
source_grades = Counter()
joint_distributions = {sp.Rational(1, 4): Counter(), sp.Rational(3, 4): Counter()}
orbit_signatures = Counter()
unseen = set(roots)
while unseen:
    b = min(unseen)
    orbit = (b, s(b), s(s(b)))
    unseen.difference_update(orbit)
    cycled = [orbit[k:] + orbit[:k] for k in range(3)]
    data = min(tuple((grade(q), charge(q)) for q in o) for o in cycled)
    if len(set(orbit)) == 1:
        orbit_signatures["fixed_root"] += 1
    elif len({a for a, _ in data}) == 1:
        orbit_signatures["class_constant_three_cycle"] += 1
    else:
        orbit_signatures["two_grade2_edges_three_cycle"] += 1
        assert sum((grade(s(q))-grade(q)) % 4 == 2 for q in orbit) == 2
    assert sum(charge(s(q))-charge(q) for q in orbit) == 0

for b in sorted(roots):
    moved = s(b)
    assert moved in roots
    g1, g2 = grade(b), grade(moved)
    dg = (g2 - g1) % 4
    assert dg == (b[5] - b[3]) % 4 and dg in (0, 2)
    dy = charge(moved) - charge(b)
    class_table[g1][g2] += 1
    difference_counts[dg][dy] += 1
    if dg == 2:
        change_type["integer" if max(map(abs, b)) == 2 else "half_integer"] += 1
        # Exact full Cartan restoration needs a charge gamma=beta-sigma beta.
        gamma = tuple(a-c for a,c in zip(b,moved))
        h = sum(a*a for a in gamma) // 8
        assert sum(a*a for a in gamma) == 8*h
        assert sum(a*c for a,c in zip(gamma,moved)) == -4*h
        compensation_weights[h] += 1
        if h == 1:
            assert gamma in roots and grade(gamma) == 2
            weight_one_compensation_charges[abs(dy)] += 1

    source_grades[(dg, source_sector(b), source_sector(moved))] += 1
    if b[5:] == moved[5:] and dy != 0:
        family_unchanged.setdefault(dg, (b, moved, dy))
    before = m.T * sp.Matrix(b) / 2
    after = m.T * sp.Matrix(moved) / 2
    assert r * before == after
    for alpha, distribution in joint_distributions.items():
        delta = sp.Matrix([(k + alpha) / 8 - sp.Rational(1, 2) for k in range(8)])
        dk = delta.dot(after - before)
        distribution[(dg, dy, dk)] += 1

assert class_table == [[28,0,24,0],[0,32,0,32],[24,0,36,0],[0,32,0,32]]
assert change_type == {"integer": 48, "half_integer": 64}
assert compensation_weights == {1: 36, 2: 72, 3: 4}
assert weight_one_compensation_charges == {sp.Rational(1, 3): 18, sp.Rational(1, 2): 18}
assert sum(difference_counts[2].values()) == 112
assert all(c[-k] == c[k] for c in difference_counts.values() for k in c)
assert difference_counts[0][0] == 42 and difference_counts[2][0] == 0
assert set(difference_counts[0]).isdisjoint(difference_counts[2])
assert orbit_signatures == {
    "fixed_root": 12, "class_constant_three_cycle": 20, "two_grade2_edges_three_cycle": 56,
}
assert sum(difference_counts[2][q] for q in (-sp.Rational(1, 2), sp.Rational(1, 2))) == 42
assert all(d[(0, 0, 0)] == 12 and d[(2, 0, 0)] == 0 for d in joint_distributions.values())
assert all(sum(d.values()) == 240 for d in joint_distributions.values())

# No arbitrary family-only scalar function on the last-three projection
# can repair fixed hypercharge. These two grade-zero roots have the same
# family projection and different Y.
b0 = 2 * (eye[:, 0] - eye[:, 4])
b1 = sigma * b0
assert b0[5:] == b1[5:] == [0, 0, 0]
assert grade(tuple(b0)) == grade(tuple(b1)) == 0
assert charge(tuple(b0)) != charge(tuple(b1))

# A linear family Cartan countercharge has no chance: a wholly D5 orbit
# of axes already carries unequal Y weights.
assert [y.dot(eye[:, k]) for k in [0, 2, 4]] == [-sp.Rational(1, 3), -sp.Rational(1, 3), sp.Rational(1, 2)]

# The literal neutral return word additionally keeps the original CAR
# sectors (0,s,0). Replay the existing Chevalley coefficients on it.
sys.path.insert(0, str(repo))
from tfpt_explorer.current_block_geometry import _actual_current_source
chevalley = _actual_current_source()["chevalley"]
assert set(chevalley.ridx) == roots
neutral_return_words = []

def add(a, b):
    return tuple(x + y for x, y in zip(a, b))

def bracket_coeff(left, right, target):
    value = chevalley.bracket(chevalley.ridx[left], chevalley.ridx[right])
    assert set(value) == {chevalley.ridx[target]}
    return int(value[chevalley.ridx[target]])

for a,b in product((-1, 1), repeat=2):
    beta = tuple(2*(a*eye[:,2]+b*eye[:,3]))
    middle = s(beta)
    third = s(middle)
    assert [grade(k) for k in (beta,middle,third)] == [0,2,0]
    assert [source_sector(k) for k in (beta,middle,third)] == ["0","s","0"]
    w = tuple(2*b*(eye[:,3]-eye[:,5]))
    d5 = tuple(2*a*(eye[:,2]-eye[:,4]))
    assert add(add(middle,w),d5) == beta
    assert not chevalley.bracket(chevalley.ridx[w], chevalley.ridx[d5])
    weak_first = bracket_coeff(w,middle,add(w,middle))*bracket_coeff(d5,add(w,middle),beta)
    d5_first = bracket_coeff(d5,middle,add(d5,middle))*bracket_coeff(w,add(d5,middle),beta)
    assert weak_first == d5_first and abs(weak_first) == 1
    drifts = {}
    for alpha in joint_distributions:
        d = m * sp.Matrix([(k + alpha) / 8 - sp.Rational(1, 2) for k in range(8)])
        clock_shift = d.dot(sp.Matrix(middle) - sp.Matrix(beta))/2
        weak_shift = d.dot(sp.Matrix(w))/2
        d5_shift = d.dot(sp.Matrix(d5))/2
        assert clock_shift + weak_shift + d5_shift == 0
        drifts[str(alpha)] = {
            "clock_delta_K": str(clock_shift),
            "weak_delta_K": str(weak_shift),
            "D5_delta_K": str(d5_shift),
            "between_insertions_K_weak_first_minus_D5_first": str(weak_shift-d5_shift),
        }
    neutral_return_words.append({
        "beta": print_root(beta), "middle": print_root(middle), "third": print_root(third),
        "middle_hypercharge_shift": str(charge(middle)-charge(beta)),
        "compensator_weight": int(sum(z*z for z in add(w,d5)) // 8),
        "weak_current": print_root(w), "D5_current": print_root(d5),
        "actual_quarter_time_shifts": drifts,
        "common_two_current_coefficient_in_original_Chevalley_basis": weak_first,
    })
assert sum(1 for b in roots
           if [grade(k) for k in (b,s(b),s(s(b)))] == [0,2,0]
           and [source_sector(k) for k in (b,s(b),s(s(b)))] == ["0","s","0"]) == 4

# Coordinate-placement control: only two of the 560 split/role placements
# can make the same sigma preserve an untransported hypercharge vector.
invariant_hypercharge_slots = []
for color in combinations(range(8), 3):
    for weak in combinations(sorted(set(range(8)) - set(color)), 2):
        candidate = sp.Matrix([-2 if i in color else 3 if i in weak else 0 for i in range(8)])
        if sigma * candidate == candidate:
            invariant_hypercharge_slots.append({
                "color": [i + 1 for i in color],
                "weak": [i + 1 for i in weak],
                "D5": [i + 1 for i in sorted(color + weak)],
            })
assert [item["D5"] for item in invariant_hypercharge_slots] == [[1, 3, 5, 7, 8], [2, 4, 6, 7, 8]]

# Include actual original drift and charge in the simultaneous groupoid.
# The previous exhaustive P,F,Y stabilizer contains exactly these two.
stabilizer = [eye, sp.diag(1, 1, 1, 1, 1, 1, -1, -1)]
target_charge = m * sp.ones(8, 1)
drift_checks = []
for g in stabilizer:
    assert g * p_d == p_d * g and g * sigma == sigma * g and g * y == y and g * ell == ell
    values = []
    for alpha in (sp.Rational(1, 4), sp.Rational(3, 4)):
        d = m * sp.Matrix([(k + alpha) / 8 - sp.Rational(1, 2) for k in range(8)])
        values.append(g * d == d)
    drift_checks.append({"fixes_charge": g * target_charge == target_charge, "fixes_quarter_drifts": values})
assert drift_checks[0]["fixes_quarter_drifts"] == [True, True]
assert drift_checks[1]["fixes_quarter_drifts"] == [False, False]

# For the algebraic charge and time covariance use every root and the
# actual simultaneously transported marks. j is held invariant by sigma.
j = sp.diag(*[sp.Matrix([[0, -1], [1, 0]])]*4)
assert sigma * j == j * sigma and (sigma * y).dot(sigma * eye[:, 3]) == y[3]
assert (eye - p_d) * j * y != sp.zeros(8, 1)
assert ((eye - p_d) * j * y).dot((eye - p_d) * j * y) == sp.Rational(1, 4)
p_next, y_next, j_next = sigma * p_d * sigma.T, sigma * y, sigma * j * sigma.T
v_next = (eye-p_next)*j_next*y_next
assert v_next.dot(v_next) == sp.Rational(1, 4)
assert all(grade(b) == int((sigma*ell).dot(sp.Matrix(s(b)))) % 4 for b in roots)
assert all(charge(b) == (sigma*y).dot(sp.Matrix(s(b)))/2 for b in roots)
for alpha in joint_distributions:
    d = m * sp.Matrix([(k + alpha) / 8 - sp.Rational(1, 2) for k in range(8)])
    assert all(d.dot(sp.Matrix(b)) == (sigma*d).dot(sp.Matrix(s(b))) for b in roots)

print(json.dumps({
    "fixed_glue_class_table": class_table,
    "glue_change_by_root_type": change_type,
    "one_step_hypercharge_difference_counts_by_glue_difference": {
        str(g): {str(d): n for d, n in sorted(c.items())} for g, c in difference_counts.items()
    },
    "fixed_glue_and_family_projection_insufficient_counterexample": {
        "root": print_root(b0), "target": print_root(b1),
        "Y": str(charge(b0)), "target_Y": str(charge(b1)), "glue_before_and_after": 0,
    },
    "orbits": orbit_signatures,
    "fixed_Y_step_count": 42,
    "fixed_Y_and_fixed_K_step_count_each_quarter": 12,
    "unweighted_glue_phase_character": "1/15",
    "glue_phase_character_conditioned_on_Y_preservation": "1",
    "glue_changing_steps_with_higgs_size_Y_difference": 42,
    "full_Cartan_compensation_lattice_weights": dict(sorted(compensation_weights.items())),
    "actual_020_glue_0s0_source_two_current_compensators": neutral_return_words,
    "weight_one_full_Cartan_compensation_by_absolute_Y": {
        str(q): n for q, n in sorted(weight_one_compensation_charges.items())
    },
    "fold_source_sector_cross_counts": {
        f"{g}:{before}->{after}": n for (g,before,after), n in sorted(source_grades.items())
    },
    "one_step_glue_Y_K_counts": {
        str(alpha): {f"{g},{q},{dk}": n for (g,q,dk), n in sorted(c.items())}
        for alpha, c in joint_distributions.items()
    },
    "P_F_Y_stabilizer_with_actual_source_marks": drift_checks,
    "fixed_sigma_hypercharge_alternate_coordinate_placements": invariant_hypercharge_slots,
    "joint_transport_preserves_exact_glue_Y_K_root_readings": True,
    "joint_transport_preserves_J_hypercharge_split_norm_squared": "1/4",
}, indent=2))
