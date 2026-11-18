#!/usr/bin/env python3
"""Exact separated-point response of the stipulated E8 level-one lattice net.

This is an executable dictionary of the stipulated boundary current theory, not a
selection of its realization, physical 4D fields, state, time or raw P1 kernel.
No Gaussian current assumption is made. Cocycles and quarter characters remain
explicit. Coordinates label an orthonormal R^8, split D5 + D3 (= A3).
"""
from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from itertools import combinations, product
import json
from pathlib import Path

import sympy as sp

ZERO = (sp.S.Zero,) * 8
HALF = sp.Rational(1, 2)
V = (sp.S.Zero,) * 5 + (HALF,) * 3
# Columns are a unimodular basis of the even E8 lattice in Euclidean coordinates.
BASIS = sp.Matrix([
    [HALF, 1, -1, 0, 0, 0, 0, 0],
    [-HALF, 1, 1, -1, 0, 0, 0, 0],
    [-HALF, 0, 0, 1, -1, 0, 0, 0],
    [-HALF, 0, 0, 0, 1, -1, 0, 0],
    [-HALF, 0, 0, 0, 0, 1, -1, 0],
    [-HALF, 0, 0, 0, 0, 0, 1, -1],
    [-HALF, 0, 0, 0, 0, 0, 0, 1],
    [HALF, 0, 0, 0, 0, 0, 0, 0],
])
GRAM = BASIS.T * BASIS
BASIS_INV = BASIS.inv()


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), sp.S.Zero)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def neg(a):
    return tuple(-x for x in a)


@lru_cache(None)
def coordinates(a):
    a = tuple(sp.Rational(x) for x in a)
    if len(a) != 8:
        raise ValueError("An E8 charge has eight Euclidean coordinates")
    m = tuple(BASIS_INV * sp.Matrix(a))
    if not all(x.is_Integer for x in m):
        raise ValueError("Charge is outside the stipulated E8 lattice")
    return tuple(int(x) for x in m)


@lru_cache(None)
def epsilon(a, b):
    """Ordered lattice cocycle; diagonal zero is a declared basis convention."""
    m, n = coordinates(a), coordinates(b)
    exponent = sum(m[i] * n[j] * int(GRAM[i, j])
                   for i in range(8) for j in range(i))
    return -1 if exponent % 2 else 1


def adjoint_sign(a):
    # In this real cocycle convention V_a^dagger = epsilon(-a,a) V_-a.
    return epsilon(neg(a), a)


def grade(a):
    x = 4 * dot(V, a)
    if not x.is_Integer:
        raise ValueError("Quarter grading is not integral on this charge")
    return int(x) % 4


def roots():
    r = []
    for i, j in combinations(range(8), 2):
        for si, sj in product((-1, 1), repeat=2):
            a = list(ZERO)
            a[i], a[j] = sp.Integer(si), sp.Integer(sj)
            r.append(tuple(a))
    for signs in product((-1, 1), repeat=8):
        if sum(x == -1 for x in signs) % 2 == 0:
            r.append(tuple(HALF * x for x in signs))
    return tuple(sorted(r))


def vertex_response(charges, positions):
    """Ordered chiral vacuum response of raw lattice vertex operators.

    Pairwise distinct real rational positions are sufficient for exact probes.
    This function is a separated-point Taylor coefficient, not an integrated
    source functional with an unspecified contact prescription.
    """
    charges = tuple(tuple(sp.Rational(x) for x in a) for a in charges)
    positions = tuple(sp.Rational(z) for z in positions)
    if len(charges) != len(positions):
        raise ValueError("One position is required for each insertion")
    if len(set(positions)) != len(positions):
        raise ValueError("Coincident insertions need an additional contact prescription")
    for a in charges:
        coordinates(a)
    if any(sum(a[j] for a in charges) != 0 for j in range(8)):
        return sp.S.Zero
    value = sp.S.One
    for i, j in combinations(range(len(charges)), 2):
        exponent = dot(charges[i], charges[j])
        if not exponent.is_Integer:
            raise AssertionError("E8 lattice locality was lost")
        value *= epsilon(charges[i], charges[j]) * (positions[i] - positions[j]) ** exponent
    return sp.factor(value)


def connected_four(charges, positions):
    value = vertex_response(charges, positions)
    for pairs in (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))):
        term = sp.S.One
        for i, j in pairs:
            term *= vertex_response((charges[i], charges[j]), (positions[i], positions[j]))
        value -= term
    return sp.factor(value)


def current_response(charges=(), positions=(), cartan=()):
    """All root/Cartan current responses, using the exact Heisenberg Ward rule.

    Each Cartan insertion is (direction_in_R8, position). This recursion uses
    H-H contractions and H-V charge poles; root vertices remain exponential
    operators rather than Gaussian current variables.
    """
    charges = tuple(tuple(sp.Rational(x) for x in a) for a in charges)
    positions = tuple(sp.Rational(z) for z in positions)
    cartan = tuple((tuple(sp.Rational(x) for x in d), sp.Rational(u)) for d, u in cartan)
    all_positions = positions + tuple(u for _, u in cartan)
    if len(set(all_positions)) != len(all_positions):
        raise ValueError("Coincident current insertions require contact data")
    if any(len(d) != 8 for d, _ in cartan):
        raise ValueError("A Cartan direction has eight coordinates")
    base = vertex_response(charges, positions)

    def ward(indices):
        if not indices:
            return sp.S.One
        i, *rest = indices
        d, u = cartan[i]
        mean = sum((dot(d, a)/(u-z) for a, z in zip(charges, positions)), sp.S.Zero)
        value = mean * ward(tuple(rest))
        for j in rest:
            e, w = cartan[j]
            value += dot(d, e)/(u-w)**2 * ward(tuple(k for k in rest if k != j))
        return value

    return sp.factor(base * ward(tuple(range(len(cartan)))))


def polynomial_product(a, b, order):
    return [sum(a[j] * b[k-j] for j in range(k+1)) for k in range(order+1)]


def ns_character(order=4):
    """Full 16-Majorana NS coefficients in x=q^(1/2), vacuum factor removed."""
    p = [1] + [0] * order
    for odd in range(1, order+1, 2):
        factor = [int(sp.binomial(16, k//odd))
                  if k % odd == 0 and k//odd <= 16 else 0 for k in range(order+1)]
        p = polynomial_product(p, factor, order)
    return p


def run():
    rr = roots()
    assert len(set(rr)) == 240 and abs(BASIS.det()) == 1
    assert GRAM.det() == 1 and all(GRAM[i, i] == 2 for i in range(8))
    for a in rr:
        coordinates(a)
        assert dot(a, a) == 2
    counts = Counter(grade(a) for a in rr)
    assert [counts[i] for i in range(4)] == [52, 64, 60, 64]
    branches = Counter((dot(a[:5], a[:5])/2, dot(a[5:], a[5:])/2, grade(a)) for a in rr)
    assert branches == {(sp.S.One, sp.S.Zero, 0): 40,
                        (sp.S.Zero, sp.S.One, 0): 12,
                        (HALF, HALF, 2): 60,
                        (sp.Rational(5, 8), sp.Rational(3, 8), 1): 64,
                        (sp.Rational(5, 8), sp.Rational(3, 8), 3): 64}
    # All root pairs: actual cocycle commutator exactly cancels the integral
    # coordinate exchange sign. No identification with fermion parity is used.
    locality_pairs = 0
    for a, b in product(rr, repeat=2):
        assert epsilon(a, b) == (-1 if int(dot(a, b)) % 2 else 1) * epsilon(b, a)
        assert grade(add(a, b)) == (grade(a)+grade(b)) % 4
        locality_pairs += 1
    conjugate_pairs = 0
    for a in rr:
        assert adjoint_sign(a) * vertex_response((neg(a), a), (2, 1)) == 1
        conjugate_pairs += 1
    a = (HALF,) * 8
    b = (sp.Integer(-1), sp.Integer(-1)) + (sp.S.Zero,) * 6
    c = neg(add(a, b))
    assert dot(a, b) == -1 and c in rr
    triple = vertex_response((a, b, c), (3, 2, 1))
    quartic = vertex_response((a, neg(a), b, neg(b)), (4, 3, 2, 1))
    connected = connected_four((a, neg(a), b, neg(b)), (4, 3, 2, 1))
    assert triple != 0 and connected != 0
    assert sum(grade(x) for x in (a, b, c)) % 4 == 0
    assert vertex_response((a, b), (2, 1)) == 0
    cartan_root_response = adjoint_sign(a) * current_response(
        (a, neg(a)), (3, 1), ((a, 4),))
    assert cartan_root_response == sp.Rational(1, 3)
    e0 = (sp.S.One,) + (sp.S.Zero,) * 7
    assert current_response(cartan=((e0, 2), (e0, 1))) == 1
    assert current_response(cartan=((e0, 4), (e0, 3), (e0, 2), (e0, 1))) == sp.Rational(169, 144)
    free_ns = ns_character()
    assert free_ns == [1, 16, 120, 576, 2076]
    dictionary = []
    for root in rr:
        g = grade(root)
        dictionary.append({
            "charge": [str(x) for x in root], "quarter_grade": g,
            "quarter_character": ["1", "i", "-1", "-i"][g],
            "weight_D5": str(dot(root[:5], root[:5])/2),
            "weight_A3": str(dot(root[5:], root[5:])/2),
            "weight_total": "1", "boundary_statistics": "bosonic_current",
            "adjoint_sign": adjoint_sign(root),
        })
    return {
        "scope": "Stipulated E8 level-one vacuum lattice representation; separated-point boundary currents",
        "selection_from_raw_P1P2": False,
        "reverse_bootstrap": {
            "source": "origin_theory.tex:954-981; v6, v15, v350-v352",
            "retained": "Discrete E8 closure fixes g_car=5 and the eight normalization witnesses within its framework",
            "not_inferred": "No automatic raw-kernel field pullback, 4D statistics, field-dependent kinetic operator or source state/time",
        },
        "quarter": {
            "character_vector": [str(x) for x in V],
            "root_grades": [counts[i] for i in range(4)],
            "current_grades_including_Cartan": [60, 64, 60, 64],
            "internal_character_not_physical_time": True,
            "original_charged_quarter_control": "charged-pair-control.json, kernel-functional-control.json (unchanged)",
        },
        "dictionary": dictionary,
        "boundary_current_blocks": [
            {"Spin10_SU4_labels": "(45,1)+(1,15)", "dimension": 60, "grade": 0,
             "meaning": "40 D5 roots +12 D3 roots +8 Cartan currents"},
            {"Spin10_SU4_labels": "(10,6)", "dimension": 60, "grade": 2,
             "weights": ["1/2", "1/2"]},
            {"Spin10_SU4_labels": "(16,4) and conjugate, label orientation conventional",
             "dimension_each": 64, "grades": [1, 3], "weights": ["5/8", "3/8"]},
        ],
        "Cartan_currents": {"count": 8, "weight": "1", "grade": 0,
                           "definition": "H^a(z)=i partial X^a(z)",
                           "OPE": "H^a(z)H^b(w) ~ delta_ab/(z-w)^2"},
        "vertex_definition": "V_alpha(z)=e_alpha :exp(i alpha.X(z)):; declared ordered lattice cocycle",
        "response_formula": "delta(sum alpha_i=0) product_(i<j) epsilon(alpha_i,alpha_j) (z_i-z_j)^(alpha_i.alpha_j)",
        "Cartan_Ward_rule": "mean_i=sum_j d_i.alpha_j/(u_i-z_j); covariance_ij=d_i.d_j/(u_i-u_j)^2; recurse over means and pairings multiplying the root-vertex response",
        "formal_current_source_series": "Z[J]/Z[0]=sum_n (1/n!) integral <product_i (sum_alpha J_alpha(z_i)V_alpha(z_i)+sum_a J_a(z_i)H^a(z_i))>; current_response computes mixed coefficients only away from diagonals",
        "checks": {"locality_root_pairs": locality_pairs, "positive_adjoint_two_point_pairs": conjugate_pairs,
                   "known_vacuum_NS_character_x": free_ns,
                   "E8_character_q_first": [1, 248, 4124],
                   "E8_character_scope": "Known original v156 control; no new equality of full NS and E8 local fields",
                   "positivity_scope": "Root-state norms in the stipulated positive E8 lattice Fock representation, not a raw P1 reflection-positivity proof"},
        "nonabelian_probe": {
            "charges": [[str(y) for y in x] for x in (a, b, c)],
            "grades": [grade(x) for x in (a, b, c)],
            "three_point_at_3_2_1": str(triple),
            "four_point_a_minus_a_b_minus_b_at_4_3_2_1": str(quartic),
            "connected_four_point": str(connected),
            "adjoint_normalized_H_a_V_a_V_minus_a_at_4_3_1": str(cartan_root_response),
            "meaning": "Current responses are non-Gaussian even in the lattice/free-boson realization; no native Yukawa or independent boson-pair transfer is inferred",
        },
        "requires_for_physical_P1": [
            "An original field insertion map into this representation, with local domains, charge, statistics and physical time",
            "The actual map from physical field values to D_rel[Phi], D_ref[Phi] and the original spectral profile",
            "Original absolute amplitude, determinant-line phase and contact/geometry prescription",
        ],
        "new_physics_claim": False,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n")
    print(json.dumps({k: result[k] for k in ("scope", "quarter", "checks", "nonabelian_probe")}, ensure_ascii=False))
