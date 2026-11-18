#!/usr/bin/env python3
"""Exact finite check of the A7 clock direction and the E8 -> E7 x U1 root census."""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
import json
from pathlib import Path


def require(name, condition):
    if not condition:
        raise RuntimeError(f"FAIL: {name}")


def dot(left, right):
    return sum((x * y for x, y in zip(left, right)), F(0))


def e8_roots():
    roots = []
    for i, j in combinations(range(8), 2):
        for si, sj in product((-1, 1), repeat=2):
            v = [F(0)] * 8
            v[i], v[j] = F(si), F(sj)
            roots.append(tuple(v))
    for signs in product((-1, 1), repeat=8):
        if sum(s < 0 for s in signs) % 2 == 0:
            roots.append(tuple(F(s, 2) for s in signs))
    return tuple(roots)


def as_string_dict(counter):
    return {str(k): counter[k] for k in sorted(counter)}


def main():
    roots = e8_roots()
    q = (F(1),) * 8
    rho = tuple(F(7 - 2 * j, 2) for j in range(8))
    pstar = (F(1, 2),) * 4 + (F(-1, 2),) * 4
    charge = lambda v: sum(v, F(0))
    is_integer = lambda v: all(x.denominator == 1 for x in v)

    require("240 distinct norm-two E8 roots",
            len(roots) == len(set(roots)) == 240
            and all(dot(v, v) == 2 for v in roots))

    # delta_j = c_j + a/8.  Compare constant and a coefficients exactly.
    delta_constants = tuple(F(j, 8) - F(1, 2) for j in range(8))
    rhs_constants = tuple(-rho[j] / 8 - F(1, 16) for j in range(8))
    require("delta is -rho_A7/8 plus the uniform U1 term",
            delta_constants == rhs_constants)

    require("pstar is omega4 and an E8 root",
            pstar in roots and dot(pstar, pstar) == 2
            and dot(rho, pstar) == 8 and charge(pstar) == 0)

    by_q = Counter(charge(v) for v in roots)
    expected_q = {F(-4): 1, F(-2): 56, F(0): 126, F(2): 56, F(4): 1}
    require("complete E8 charge census", dict(by_q) == expected_q)

    by_q_type = Counter((charge(v), "integer" if is_integer(v) else "half")
                        for v in roots)
    expected_type = {
        (F(-4), "half"): 1,
        (F(-2), "integer"): 28, (F(-2), "half"): 28,
        (F(0), "integer"): 56, (F(0), "half"): 70,
        (F(2), "integer"): 28, (F(2), "half"): 28,
        (F(4), "half"): 1,
    }
    require("integer and half-root refinement", dict(by_q_type) == expected_type)

    neutral = {v for v in roots if charge(v) == 0}
    a7_roots = {v for v in neutral if is_integer(v)}
    lambda4 = {v for v in neutral if not is_integer(v)}
    expected_a7 = set()
    for i, j in combinations(range(8), 2):
        for sign in (-1, 1):
            v = [F(0)] * 8
            v[i], v[j] = F(sign), F(-sign)
            expected_a7.add(tuple(v))
    expected_lambda4 = set()
    for plus in combinations(range(8), 4):
        chosen = set(plus)
        expected_lambda4.add(tuple(F(1, 2) if j in chosen else F(-1, 2)
                                   for j in range(8)))
    require("neutral shell is A7 roots plus Lambda4 weights",
            a7_roots == expected_a7 and lambda4 == expected_lambda4
            and len(a7_roots) == 56 and len(lambda4) == 70)

    # The neutral roots form a rank-seven reflection-closed subsystem.
    require("neutral root span has rank seven",
            __import__("sympy").Matrix([[int(2 * x) for x in v]
                                         for v in neutral]).rank() == 7)
    reflection_closed = True
    for alpha in neutral:
        for beta in neutral:
            reflected = tuple(beta[j] - dot(alpha, beta) * alpha[j]
                              for j in range(8))
            if reflected not in neutral:
                reflection_closed = False
                break
        if not reflection_closed:
            break
    require("neutral roots are reflection closed", reflection_closed)

    # For a root state J_{alpha,-1}|0>, the H8 commutator coefficient is
    # 8-rho.alpha+(2a-1)Q(alpha)/2.  Record the a-independent base part.
    grade_base = {}
    for Q in sorted(by_q):
        grade_base[str(Q)] = as_string_dict(Counter(
            F(8) - dot(rho, v) for v in roots if charge(v) == Q
        ))
    require("unique neutral zero-grade root is pstar",
            [v for v in neutral if F(8) - dot(rho, v) == 0] == [pstar])

    result = {
        "root_count": len(roots),
        "delta_identity": "delta=-rho_A7/8+(2a-1)q/16",
        "rho_A7": [str(x) for x in rho],
        "pstar": {
            "vector": [str(x) for x in pstar],
            "norm2": str(dot(pstar, pstar)),
            "rho_dot": str(dot(rho, pstar)),
            "Q": str(charge(pstar)),
            "H8_root_commutator_grade": "0",
        },
        "charge_census": as_string_dict(by_q),
        "charge_type_census": {
            f"Q={Q},{kind}": count
            for (Q, kind), count in sorted(by_q_type.items(), key=lambda x: (x[0][0], x[0][1]))
        },
        "neutral_decomposition": {
            "E7_roots": len(neutral),
            "A7_integer_roots": len(a7_roots),
            "Lambda4_half_weights": len(lambda4),
            "rank": 7,
            "reflection_closed": reflection_closed,
        },
        "commutator": "[H8,J_(alpha,-1)]=(8-rho_A7.alpha+(2a-1)Q(alpha)/2)J_(alpha,-1)",
        "grade_base_by_charge": grade_base,
    }
    out = Path(__file__).with_name("a7-e7-clock-results.json")
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
