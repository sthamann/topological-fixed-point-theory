#!/usr/bin/env python3
"""Exact existing joint-lift/native-clock comparison; no source selection.

The comparison evaluates a supplied toral clutch loop at the SAME declared
angles. A clutch transition is not automatically connection holonomy or time.
The old full-inner conjugacy theorem is reused, not recomputed here.
"""
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, product
import argparse
import json
from pathlib import Path


BASE = Path(__file__).resolve().parent
REPO = BASE.parents[2]
PINS = {
    "verification/v1_e8_glue.py":
        "1978d66a85974e5fa189cf4dec507f9be4c4dcc113e90eb1b33a9aae9259a0d7",
    "experiments/theory-contracts/e-series-source-reconstruction-20261024/"
    "NATIVE-CLOCK-CONJUGACY.txt":
        "5f208d8d32fdc77e774ec77b5a320f378cfc9e2c8793d4d3d439983aae583a45",
    "experiments/theory-contracts/e-series-source-reconstruction-20261024/"
    "native-clock-conjugacy-results.json":
        "959c9c618320125c98587c9c958d965ade1cff41ccbdf0ba5b67a68b5326d267",
}
NATIVE_RESULTS = REPO / (
    "experiments/theory-contracts/e-series-source-reconstruction-20261024/"
    "native-clock-conjugacy-results.json"
)
CHECKS = []


def require(name, condition):
    if not condition:
        raise RuntimeError(name)
    CHECKS.append(name)


def roots_doubled():
    """Exactly the standard 112 integral plus 128 half-integral E8 roots."""
    roots = []
    for i, j in combinations(range(8), 2):
        for a, b in product((-2, 2), repeat=2):
            root = [0] * 8
            root[i], root[j] = a, b
            roots.append(tuple(root))
    roots += [s for s in product((-1, 1), repeat=8) if s.count(-1) % 2 == 0]
    return roots


def cocharacter(k, m):
    # S=det(E)^(-1/2) Lambda-even(E) has ODD D5 minus parity.
    # Therefore its paired family 4 has tetrahedron -T, not T.
    return tuple(k) + tuple(-m[0] - m[j] - 1 for j in (1, 2, 3))


def root_grades(roots, k, m):
    h = cocharacter(k, m)
    values = [sum(x * y for x, y in zip(root, h)) for root in roots]
    if any(x % 2 for x in values):
        raise RuntimeError("nonintegral E8 root pairing")
    result = Counter(x // 2 for x in values)
    result[0] += 8
    return result


def branch_grades(k, m):
    """Independent weights from the actual five E8 branches, including W_m."""
    result = Counter({0: 8})
    # (45,1): the D5 adjoint, with its five Cartans already included.
    for i, j in combinations(range(5), 2):
        result.update(a * k[i] + b * k[j] for a, b in product((-1, 1), repeat=2))
    # (1,15): beta_i-beta_j=m_i-m_j, plus three Cartans above.
    result.update(m[i] - m[j] for i in range(4) for j in range(4) if i != j)
    # (16,4) is the literal Lambda-even(E) tensor sum_j det(E)^m_j.
    for count in (0, 2, 4):
        for subset in combinations(range(5), count):
            weight = sum(k[i] for i in subset)
            for mj in m:
                result.update((weight + mj, -weight - mj))
    # (10,6), with six family weights beta_i+beta_j=m_i+m_j+1.
    for ki in k:
        for i, j in combinations(range(4), 2):
            result.update((ki + m[i] + m[j] + 1, -ki + m[i] + m[j] + 1))
    return result


def phases(grades, modulus, sign=1):
    result = [0] * modulus
    for grade, multiplicity in grades.items():
        result[(sign * grade) % modulus] += multiplicity
    return tuple(result)


def lift_residues(modulus):
    for first in product(range(modulus), repeat=3):
        yield first + (-2 - sum(first),)


def quarter_order(character):
    return tuple(character[j] for j in (0, 2, 1, 3))


def toral_character(roots, h, modulus):
    values = [sum(x * y for x, y in zip(root, h)) for root in roots]
    if any(x % 2 for x in values):
        raise RuntimeError("nonintegral toral cocharacter")
    counter = Counter(x // 2 for x in values)
    counter[0] += 8
    return phases(counter, modulus)


def parallel_split_clock(roots, target):
    """Whole native compiler commuting with EVERY supplied split connection.

    This clock is a separate bundle automorphism, not a clutch evaluation.
    All243 toral thirds centralizing fixed U3xU2 are classified, not all
    possible marked field maps or physical connections.
    """
    matches = []
    for color, weak, x, y, z in product(range(3), repeat=5):
        b = [color, color, color, weak, weak, x, y, z]
        if sum(b) % 2:
            b[5] += 3  # unique E8 lattice class with these standard mod3 residues
        h12 = tuple((6 if i < 5 else 0) + 4 * value for i, value in enumerate(b))
        if toral_character(roots, h12, 12) == target:
            matches.append(tuple(b))
    expected = {
        (color, color, color, 0, 0) + tuple(phase if i == axis else 0 for i in range(3))
        for color in (1, 2) for axis in range(3) for phase in (1, 2)
    }
    require("12 of243 split-connection centralizing thirds match native compiler", (
        len(matches) == 12 and {tuple(v % 3 for v in b) for b in matches} == expected
    ))
    b = (1, 1, 1, 0, 0, 0, 0, 1)
    h12 = tuple((6 if i < 5 else 0) + 4 * value for i, value in enumerate(b))
    kappa4 = (2, 2, 2, 2, 2, 0, 0, 0)
    require("simple globally split-compatible compiler matches actual whole character", (
        b in matches and toral_character(roots, h12, 12) == target
    ))
    require("actual ninth/fourth powers on all root fields", all(
        (9 * sum(a * x for a, x in zip(r, h12)) - 3 * sum(a * x for a, x in zip(r, kappa4))) % 24 == 0
        and (4 * sum(a * x for a, x in zip(r, h12)) - 4 * sum(a * x for a, x in zip(r, b))) % 24 == 0
        for r in roots
    ))
    split_generators = []
    for block in ((0, 1, 2), (3, 4)):
        for i, j in combinations(block, 2):
            r = [0] * 8
            r[i], r[j] = 1, -1
            split_generators.append(tuple(r))
    require("sigma centralizes actual U3xU2 root generators", all(
        sum(x * y for x, y in zip(b, r)) % 3 == 0 for r in split_generators
    ))
    tetrahedron = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
    family_phases = (1, 2, 2, 1)
    correct_field_weights = True
    for count in (0, 2, 4):
        for subset in combinations(range(5), count):
            eps = tuple(1 if i in subset else -1 for i in range(5))
            nc = len(set(subset) & {0, 1, 2})
            for j, t in enumerate(tetrahedron):
                doubled_root = eps + tuple(-s for s in t)
                q2 = sum(x * y for x, y in zip(doubled_root, b))
                correct_field_weights &= doubled_root in roots and q2 % 2 == 0
                correct_field_weights &= (q2 // 2 - nc - family_phases[j]) % 3 == 0
    require("sigma literal W_m action is color center times two-plus-two family matrix", correct_field_weights)
    require("kappa is trivial on ALL K adjoint root currents", all(
        sum(x * y for x, y in zip(r, kappa4)) % 8 == 0
        for r in roots if all(x == 0 for x in r[:5]) or all(x == 0 for x in r[5:])
    ))
    return {
        "toral_thirds_checked": 243, "native_matches": 12, "matching_b": matches,
        "simple_b": b, "simple_c_h12": h12,
        "sigma_on_Wm": "Lambda-even(diag(omega I3,I2)) tensor diag(omega,omega^2,omega^2,omega)",
        "connection_scope": "Every GIVEN U3xU2 connection and its natural determinant-family connection, for every m. No independent nonabelian family connection assumed.",
        "globality": "kappa and sigma commute with all natural split transition functions; they are global parallel bundle automorphisms.",
        "not_identified": "Clutch angle, geometric base transport, physical time, Phi, source domain and source-state map remain separate.",
    }


def bounded_color_weak_clutches(roots, target):
    """Complete minimum proof in the stated integer toral class.

    An explicit match has norm squared 100. All lower norms have |k_i|<=10
    and |m_j+1/2|<=10. Sorting only uses S3 color, S2 weak and S4 family.
    No bound on an unknown minimizer is assumed.
    """
    families = [
        m for m in combinations_with_replacement(range(-10, 10), 4)
        if sum(m) == -2 and sum((2 * x + 1) ** 2 for x in m) <= 396
    ]
    colors = [
        k for k in combinations_with_replacement(range(-10, 11), 3)
        if sum(k) == 0
    ]
    carriers = [
        k + (a, 1 - a) for k in colors for a in range(-9, 1)
        if sum(x * x for x in k) + a * a + (1 - a) ** 2 <= 99
    ]
    matches = []
    checked = 0
    all_branches_agree = True
    for k in carriers:
        for m in families:
            h = cocharacter(k, m)
            norm = sum(x * x for x in h)
            if norm > 100:
                continue
            checked += 1
            grades = root_grades(roots, k, m)
            all_branches_agree &= grades == branch_grades(k, m)
            if phases(grades, 12, -1) == target:
                matches.append({"norm_squared": norm, "k": k, "m": m, "h": h})
    require("all bounded clutches agree with literal E8 branch weights", all_branches_agree)
    require("complete bounded search has 6939 ordered-weight representatives", checked == 6939)
    require("155 full-clock matches at norm 100; none below", (
        len(matches) == 155 and {x["norm_squared"] for x in matches} == {100}
    ))
    return {"checked": checked, "norm_bound": 100, "matching_representatives": matches}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=BASE / "joint-lift-native-clock-results.json")
    args = parser.parse_args()
    for path, expected in PINS.items():
        require("source pin: " + path, sha256((REPO / path).read_bytes()).hexdigest() == expected)
    native = json.loads(NATIVE_RESULTS.read_text())
    quarter_native = tuple(native["V1_character_1_i_minus1_minusi"][j] for j in (0, 2, 1, 3))
    compiler_native = tuple(native["native_compiler_C12_character"])
    require("actual native character conventions", quarter_native == (60, 60, 64, 64))
    require("actual native C12 fixed dimension", compiler_native[0] == 22)
    require("full actual native C12 character", compiler_native == (22, 20, 19, 24, 19, 20, 22, 20, 19, 24, 19, 20))
    roots = roots_doubled()
    require("240 exact roots, doubled norm squared 8", (
        len(set(roots)) == 240 and all(sum(x * x for x in r) == 8 for r in roots)
    ))
    weak = (0, 0, 0, 1, 0)
    types = defaultdict(list)
    for m in lift_residues(4):
        grades = root_grades(roots, weak, m)
        require("quarter branches " + str(m), grades == branch_grades(weak, m))
        types[quarter_order(phases(grades, 4))].append(tuple(x % 4 for x in m))
    expected_types = {(60, 60, 64, 64): 8, (70, 66, 56, 56): 18,
                      (82, 54, 56, 56): 24, (92, 28, 64, 64): 8,
                      (134, 2, 56, 56): 6}
    require("complete 64-class quarter intersection", {k: len(v) for k, v in types.items()} == expected_types)
    native_residues = set(types[quarter_native])
    expected_residues = {tuple(0 if i == j else 2 for i in range(4)) for j in range(4)}
    expected_residues |= {tuple(3 if i == j else 1 for i in range(4)) for j in range(4)}
    require("exact eight native residue classes", native_residues == expected_residues)
    minimum_lifts = []
    for m in product(range(-3, 3), repeat=4):
        if sum(m) != -2 or tuple(x % 4 for x in m) not in native_residues:
            continue
        norm = sum(x * x for x in cocharacter(weak, m))
        if norm <= 12:
            minimum_lifts.append({"m": m, "norm_squared": norm})
    require("24 norm12 quarter lifts, none smaller", (
        len(minimum_lifts) == 24 and {x["norm_squared"] for x in minimum_lifts} == {12}
    ))
    # C=exp(-2pi i ad(h)/12), C^9=Q and C^4=exp(-2pi i ad(h)/3).
    # The fixed D4 gives an unbounded analytic exclusion, not only a scan.
    full_types = Counter()
    for m in lift_residues(12):
        full_types[phases(root_grades(roots, weak, m), 12, -1)] += 1
    require("all1728 pure weak classes miss full C12", compiler_native not in full_types)
    require("pure weak fixed dimension lower bound32 is attained", min(c[0] for c in full_types) == 32)
    weak_only_matches = 0
    weak_only_minfixed = 248
    for a in range(12):
        k = (0, 0, 0, a, 1 - a)
        for m in lift_residues(12):
            char = phases(root_grades(roots, k, m), 12, -1)
            weak_only_matches += char == compiler_native
            weak_only_minfixed = min(weak_only_minfixed, char[0])
    require("complete20736 weak-toral classes miss full compiler", weak_only_matches == 0)
    require("weak-toral fixed count alone does not suffice", weak_only_minfixed == 22)
    example_k, example_m = (-6, 0, 6, 0, 1), (-5, 1, 1, 1)
    example_h = cocharacter(example_k, example_m)
    example_grades = root_grades(roots, example_k, example_m)
    require("explicit same-Chern-class full-clock realization", (
        sum(example_k[:3]) == 0 and sum(example_k[3:]) == 1 and sum(example_m) == -2
        and sum(x * x for x in example_h) == 100
        and example_grades == branch_grades(example_k, example_m)
        and phases(example_grades, 12, -1) == compiler_native
    ))
    bounded = bounded_color_weak_clutches(roots, compiler_native)
    split_clock = parallel_split_clock(roots, compiler_native)
    # Whole64 carrier has degreezero for EVERY such k,m, independent of clocks.
    whole_degree = 4 * 8 * sum(example_k) + 16 * sum(example_m)
    require("whole joint64 clutch degree remains zero", whole_degree == 0)
    result = {
        "research_id": "MATH.ORIGINAL-SEAM.JOINT-LIFT-NATIVE-CLOCK.01",
        "verdict": "EXACT_SCOPED_COMPARISON_PHYSICAL_SOURCE_NOT_DERIVED",
        "quarter_phase_order": ["1", "-1", "i", "-i"],
        "quarter_types": [{"character": k, "classes": len(v)} for k, v in sorted(types.items())],
        "native_quarter_residues": sorted(native_residues),
        "quarter_minimum_norm_squared": 12,
        "quarter_minimum_lifts": minimum_lifts,
        "pure_weak_c12_classes": 1728, "pure_weak_c12_matches": 0,
        "pure_weak_fixed_dimension_lower_bound": 32,
        "all_weak_toral_c12_classes": 20736, "all_weak_toral_c12_matches": weak_only_matches,
        "actual_native_compiler_character": compiler_native,
        "same_Chern_class_example": {"k": example_k, "m": example_m, "h": example_h,
                                    "character": phases(example_grades, 12, -1), "norm_squared": 100},
        "full_clock_norm_minimum_in_stated_toral_class": bounded,
        "global_parallel_split_clock": split_clock,
        "whole64_determinant_winding": whole_degree,
        "conjugacy_scope": "Old complete classification gives full inner E8/VOA conjugacy transporting ALL marks; no fixed Cartan or physical source map.",
        "extra_comparison_premise": "Evaluate a SUPPLIED clutch cocharacter at the same angles. This is not a derived connection holonomy or physical time.",
        "source_sha256": {p: sha256((REPO / p).read_bytes()).hexdigest() for p in list(PINS) + [
            "experiments/theory-contracts/source-joint-spin-family-lift-20260922/PROOF.md",
            "experiments/theory-contracts/source-boundary-lift-selection-20260922/PROOF.md",
            str(NATIVE_RESULTS.relative_to(REPO))]},
        "checks": CHECKS, "promotion": False, "complete_TFPT_solution": False,
    }
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"verdict": "PASS", "checks": len(CHECKS), "quarter_matches": 8,
                      "pure_weak_c12_matches": 0, "same_Chern_class_full_compiler": True,
                      "full_clock_minimum_norm_squared": 100,
                      "minimum_weight_representatives": 155,
                      "parallel_split_native_compilers": 12, "output": str(args.output)}))


if __name__ == "__main__":
    main()
