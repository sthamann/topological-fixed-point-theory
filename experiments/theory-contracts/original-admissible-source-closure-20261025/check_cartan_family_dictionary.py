#!/usr/bin/env python3
"""Exact finite parts of the Cartan-subspace/family-factor dictionary.

The normalizer, bundle-reduction and cover statements are proved in the
companion text. Enumeration is not substituted for those proofs.
"""
from collections import Counter
from itertools import combinations, combinations_with_replacement, product
from pathlib import Path
import hashlib
import json

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TARGET = (22, 20, 19, 24, 19, 20, 22, 20, 19, 24, 19, 20)
SIGNS = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
CHECKS = []


def require(value, name):
    if not value:
        raise RuntimeError(name)
    CHECKS.append(name)


def roots_doubled():
    roots = []
    for i, j in combinations(range(8), 2):
        for a, b in product((-2, 2), repeat=2):
            r = [0] * 8
            r[i], r[j] = a, b
            roots.append(tuple(r))
    roots.extend(r for r in product((-1, 1), repeat=8)
                 if sum(x < 0 for x in r) % 2 == 0)
    return roots


def phase_counts(residues, degree, symmetric):
    iterator = combinations_with_replacement if symmetric else combinations
    counts = Counter(sum(residues[i] for i in c) % 3
                     for c in iterator(range(len(residues)), degree))
    return tuple(counts[j] for j in range(3))


def native_cartan(roots):
    j = sp.diag(*([sp.Matrix([[0, -1], [1, 0]])] * 4))
    sigma = sp.eye(8)[[4, 5, 0, 1, 2, 3, 6, 7], :]
    c = j * sigma
    hplus = sp.zeros(8, 4)
    for k in range(4):
        hplus[2*k, k], hplus[2*k+1, k] = 1, -sp.I
    sigma4 = sigma.extract([0, 2, 4, 6], list(range(8))) * hplus
    require(j*hplus == sp.I*hplus, "original positive Cartan polarization")
    require(sigma*hplus == hplus*sigma4, "actual sigma restricts to Cartan four")
    require(c**9 == j and c**4 == sigma and c**12 == sp.eye(8),
            "actual marked C12 power identities")
    t = sp.Symbol("t")
    require(sigma4.charpoly(t).as_expr().expand() ==
            ((t-1)**2*(t**2+t+1)).expand(), "Gaussian third character 2,1,1")
    require(c.charpoly(t).as_expr().expand() ==
            ((t**2+1)**2*sp.cyclotomic_poly(12, t)).expand(),
            "original Cartan Phi4 squared times Phi12")
    witness = hplus * sp.Matrix([1, 10, 100, 1000])
    require(all((sp.Matrix([r]).dot(witness)) != 0 for r in roots),
            "explicit regular semisimple point in original Cartan four")

    # A genuine integral E8 root basis; Smith form is not computed on Z8.
    e = sp.eye(8)
    columns = [sp.Matrix([1, -1, -1, -1, -1, -1, -1, 1])/2,
               e[:, 0]+e[:, 1], -e[:, 0]+e[:, 1]]
    columns.extend(-e[:, k]+e[:, k+1] for k in range(1, 6))
    basis = sp.Matrix.hstack(*columns)
    gram = basis.T*basis
    require(gram.det() == 1 and all(x.is_Integer for x in gram),
            "even unimodular integral E8 root basis")
    lattice_j = basis.inv()*j*basis
    require(all(x.is_Integer for x in lattice_j), "native J preserves E8 lattice")
    snf = smith_normal_form(sp.eye(8)-lattice_j, domain=sp.ZZ)
    diagonal = sorted(abs(int(snf[k, k])) for k in range(8))
    require(diagonal == [1]*4+[2]*4, "torus correction kernel is F2 to fourth")
    return {"Gaussian_third_multiplicities": [2, 1, 1],
            "Gaussian_C12_exponents_J_positive": [3, 3, 7, 11],
            "regular_Cartan_point_coefficients": [1, 10, 100, 1000],
            "E8_integral_I_minus_J_Smith": diagonal,
            "torus_isogeny_degree": int(abs((sp.eye(8)-lattice_j).det()))}


def all_native_family_characters(roots):
    matches = []
    for residues in product(range(3), repeat=8):
        b = list(residues)
        if sum(b) % 2:
            b[5] += 3
        h = [(6 if i < 5 else 0)+4*b[i] for i in range(8)]
        counts = [0]*12
        counts[0] = 8
        for r in roots:
            value = sum(x*y for x, y in zip(r, h))
            if value % 2:
                raise RuntimeError("nonintegral E8 current phase")
            counts[(value//2) % 12] += 1
        if tuple(counts) != TARGET:
            continue
        parity = sum(b[5:]) % 2
        family = []
        for signs in SIGNS:
            numerator = -sum(x*y for x, y in zip(signs, b[5:]))+3*parity
            if numerator % 2:
                raise RuntimeError("nonintegral canonical order-three family lift")
            family.append((numerator//2) % 3)
        spin = Counter(x % 3 for q in b[:5] for x in (q, -q))
        matches.append({"b_mod3": list(residues),
                        "family_order3_residues": family,
                        "family_multiplicities": [family.count(k) for k in range(3)],
                        "Spin10_vector_multiplicities": [spin[k] for k in range(3)]})
    require(len(matches) == 480, "reused complete native class has 480 toral thirds")
    require(all(m["family_multiplicities"] == [0, 2, 2] for m in matches),
            "every native third has family character 0,2,2")
    require(all(m["Spin10_vector_multiplicities"] == [4, 3, 3] for m in matches),
            "every native third has Spin10 vector character 4,3,3")
    patterns = Counter(tuple(m["family_order3_residues"]) for m in matches)
    require(len(patterns) == 6 and set(patterns.values()) == {80},
            "all six family two-plus-two patterns occur eighty times")
    return {"toral_thirds_examined": 3**8, "native_matches": len(matches),
            "family_character": [0, 2, 2], "Spin10_vector_character": [4, 3, 3],
            "family_patterns": [{"residues": list(k), "count": v}
                                for k, v in sorted(patterns.items())],
            "matching_residues": matches}


def tensor_comparison():
    gaussian, family = (0, 0, 1, 2), (1, 1, 2, 2)
    symmetric, exterior = [], []
    for n in range(11):
        a, b = phase_counts(gaussian, n, True), phase_counts(family, n, True)
        require((a == b) == (n % 3 != 1), f"symmetric character degree {n}")
        symmetric.append({"degree": n, "Gaussian": a, "family": b})
    for n in range(5):
        exterior.append({"degree": n, "Gaussian": phase_counts(gaussian, n, False),
                         "family": phase_counts(family, n, False)})
    require(exterior[2]["Gaussian"] == (2, 2, 2) and
            exterior[2]["family"] == (4, 1, 1), "exterior square detects mismatch")
    t = sp.Symbol("t")
    trace_gaussian = 1/((1-t)*(1-t**3))
    trace_family = (1-t)**2/(1-t**3)**2
    require(sp.cancel(trace_gaussian-trace_family-3*t/(1-t**3)**2) == 0,
            "all-degree symmetric trace difference identity")
    return {"symmetric_powers": symmetric, "exterior_powers": exterior,
            "all_degree_trace_difference": "3t/(1-t^3)^2",
            "direct_C3_intertwiner_isomorphism": False}


def main():
    roots = roots_doubled()
    require(len(roots) == len(set(roots)) == 240, "original E8 root inventory")
    result = {"research_id": "MATH.ORIGINAL-SEAM.CARTAN-FAMILY-DICTIONARY.01",
              "finite_verdict": "PASS_EXACT_SCOPED_DICTIONARY",
              "Cartan": native_cartan(roots),
              "whole_native_class": all_native_family_characters(roots),
              "tensor_comparison": tensor_comparison(),
              "proof_only": ["Cartan subspace maximality",
                             "little Weyl equals original C_W(E8)(J)",
                             "finite normalizer and bundle reduction obstruction",
                             "even-degree pullback removes only the stated topological obstruction"],
              "source_sha256": {}, "physical_gates_closed": [],
              "complete_TFPT_solution": False, "promotion": False}
    sources = ["verification/v1_e8_glue.py", "verification/v634_st31_structure.py",
               "verification/v690_quartic_half.py",
               "experiments/theory-contracts/e-series-compiler-closure-20261023/COVERS.txt",
               "experiments/theory-contracts/e-series-source-reconstruction-20261024/NATIVE-CLOCK-CONJUGACY.txt",
               "experiments/theory-contracts/source-joint-spin-family-lift-20260922/PROOF.md",
               str((HERE/"CARTAN-FAMILY-DICTIONARY.txt").relative_to(ROOT))]
    for source in sources:
        result["source_sha256"][source] = hashlib.sha256((ROOT/source).read_bytes()).hexdigest()
    result["checks"] = CHECKS
    (HERE/"cartan-family-dictionary-results.json").write_text(
        json.dumps(result, ensure_ascii=False, sort_keys=True, indent=2)+"\n")
    print(result["finite_verdict"])


if __name__ == "__main__":
    main()
