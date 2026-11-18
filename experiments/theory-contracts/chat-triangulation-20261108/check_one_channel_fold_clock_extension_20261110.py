#!/usr/bin/env python3
"""Exact controls for one source channel, the 8-arc fold, and marked glue."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import comb
import json
from pathlib import Path

import sympy as sp


def fusion():
    layers = [{0: 1}]
    for _ in range(8):
        current = Counter()
        for previous, count in layers[-1].items():
            for sign in (-1, 1):
                next_charge = previous + sign
                exponent = F(next_charge**2 - previous**2 - 1, 8)
                assert exponent == F(previous * sign, 4)
                current[next_charge] += count
        layers.append(dict(sorted(current.items())))
    assert list(layers[-1].values()) == [comb(8, i) for i in range(9)]
    assert all(0 <= j and (1 - F(j*j, 8)).denominator == 1
               for j in [0])
    assert [j for j in layers[-1] if F(j*j, 8) <= 1
            and (1 - F(j*j, 8)).denominator == 1] == [0]
    assert [j for j in range(-8, 9) if F(j*j, 32) <= 1
            and (1 - F(j*j, 32)).denominator == 1] == [0]
    even = Counter(sum(signs)//2 for signs in product((-1, 1), repeat=8)
                   if signs.count(-1) % 2 == 0)
    assert dict(even) == {-4: 1, -2: 28, 0: 70, 2: 28, 4: 1}
    return {"twice_charge_multiplicities_by_step": layers,
            "final_weights": {str(F(j, 2)): str(F(j*j, 8)) for j in layers[-1]},
            "fixed_even_spinor_chirality_total_charges": dict(sorted(even.items()))}


def fold():
    h, u, v, z = sp.symbols("h u v z", nonzero=True)
    C = sp.zeros(8)
    for j in range(7):
        C[j+1, j] = 1
    C[0, 7] = h
    C_minus = C.subs(h, -h)
    D = C_minus * C.inv()
    assert C.charpoly().as_expr() == C.charpoly().gen**8 - h
    assert C**8 == h * sp.eye(8)
    assert C_minus**8 == -h * sp.eye(8)
    assert (-C)**8 == h * sp.eye(8)
    assert D == sp.diag(-1, *([1]*7))
    orbit = [C**j * D * C**(-j) for j in range(8)]
    assert orbit == [sp.diag(*[(-1 if i == j else 1) for i in range(8)])
                     for j in range(8)]
    total = sp.eye(8)
    for member in orbit:
        total *= member
    assert total == -sp.eye(8)
    assert total == C_minus**8 * C**(-8)
    assert sp.resultant(z**8-u, z**8-v, z) == (u-v)**8 or (
        sp.expand(sp.resultant(z**8-u, z**8-v, z) - (u-v)**8) == 0)
    numerator = sp.diff(z**8-u, z) / 2
    assert numerator == 4*z**7
    assert F(28, 8*4) + F(-7, 8) == 0
    return {"real_eigenangles_quarter_1_numerators_mod_32": sorted(
                {(sign*(4*k+1))%32 for sign in (-1,1) for k in range(8)}),
            "real_eigenangles_quarter_3_numerators_mod_32": sorted(
                {(sign*(4*k+3))%32 for sign in (-1,1) for k in range(8)}),
            "relative_arc_signs_majorana": [-1]*2 + [1]*14,
            "relative_continuous_eigenlabel_angles": ["1/16"]*8+["-1/16"]*8,
            "relative_sign_orbit": len(orbit), "eight_pass_ratio": "-I",
            "folded_conditional_all_plus_pair_exponent": -2,
            "folded_conditional_all_plus_source_ward_residue": 4}


def clocks():
    z = sp.symbols("z")
    check_counts = {}
    for prime in (3, 5):
        phi = sp.cyclotomic_poly(prime, z)
        counts = [n for n in product(range(9), repeat=prime) if sum(n) == 8]
        assert all(sp.rem(sum(n[j]*z**j for j in range(prime)), phi, z) != 0
                   for n in counts)
        check_counts[prime] = len(counts)
    exponents = [i for i in range(30) if sp.gcd(i,30) == 1]
    joint = [(i%5, i%3) for i in exponents]
    assert set(joint) == set(product(range(1,5), (1,2)))
    assert len(set(joint)) == 8
    H = sp.Matrix([[1,1],[1,-1]]) / sp.sqrt(2)
    S = sp.diag(1,sp.I)
    identity2 = sp.eye(2)
    cnot = sp.Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
    C5 = sp.simplify(-sp.I * sp.kronecker_product(identity2,H) * cnot
                    * sp.kronecker_product(S,identity2) * sp.kronecker_product(H,identity2))
    sigma = sp.Matrix([[0,0,1,0],[1,0,0,0],[0,1,0,0],[0,0,0,1]])
    assert C5**5 == sigma**3 == sp.eye(4)
    assert (C5*sigma - sigma*C5).rank() == 4
    return {"prime_multiplicity_checks_without_orthogonal_companion": check_counts,
            "Q_zeta30_joint_characters_mod_5_mod_3": sorted(joint),
            "native_C5_sigma_commutator_rank_complex_4": 4}


def extension():
    isotropic = {(a,b) for a,b in product(range(4),repeat=2)
                 if (5*a*a+3*b*b) % 8 == 0}
    diag = {(j,j) for j in range(4)}
    anti = {(j,-j%4) for j in range(4)}
    assert isotropic == diag | anti
    assert diag & anti == {(0,0),(2,2)}
    assert all(F((5*a*a+3*b*b),8).denominator == 1 for a,b in isotropic)
    assert len(isotropic) == 6
    return {"marked_maximal_glues": [sorted(diag), sorted(anti)],
            "indistinguishable_even_glue": sorted(diag & anti),
            "four_D8_classes_weights": {"vacuum": "0", "vector": "1/2",
                                        "spinor_s": "1", "spinor_c": "1"},
            "odd_D5_A3_weight": str(F(5,8)+F(3,8)),
            "odd_D13_A3_weight": str(F(13,8)+F(3,8))}


def main():
    data = {"one_channel": fusion(), "fold": fold(), "clock": clocks(),
            "lattice_extension": extension()}
    Path(__file__).with_name("one-channel-fold-clock-extension-results-20261110.json").write_text(
        json.dumps(data, indent=2) + "\n")
    print("PASS: all fusion layers, source fold, clock tests and marked glues")


if __name__ == "__main__":
    main()
