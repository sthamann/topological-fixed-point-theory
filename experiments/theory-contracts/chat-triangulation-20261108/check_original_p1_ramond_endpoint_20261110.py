#!/usr/bin/env python3
"""Exact charge, time, and geometric checks for the marked endpoint comparison."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
import sys

ROOT = Path(sys.argv[1]).resolve()
FOLD = ROOT / "experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json"
M = [[F(x) for x in row] for row in json.loads(FOLD.read_text())["explicit_Weyl_matrix"]]
W = (1, -1, 0, -1, -1, -1, 0, -1)


def dot(a, b):
    return sum(x*y for x, y in zip(a, b, strict=True))


def reflect(p):
    return tuple(-x for x in reversed(p))


def center(p):
    return dot(W, p) % 3


def energy(alpha, p):
    delta = tuple((j+alpha)/8-F(1,2) for j in range(8))
    return 8 * (F(dot(p,p),2) + dot(delta,p))


def main():
    results = {}
    for name, a in (("e2-e3", 1), ("e1-e3", 0)):
        q = tuple(x-y for x,y in zip(M[a], M[2], strict=True))
        assert all(abs(x) == F(1,2) for x in q)
        assert dot(q,q) == 2 and dot(W,q) == 3 and sum(q) == -2
        reflected_root = [dot(row, reflect(q)) for row in M]
        expected = [0,0,0,0,-1,1,0,0] if a == 1 else [0,0,0,1,0,1,0,0]
        assert reflected_root == expected
        assert center(reflect(q)) == 0
        branches = {}
        for alpha in (F(1,4),F(3,4)):
            branches[str(alpha)] = [
                str(energy(alpha, tuple(sign*x for x in q))) for sign in (1,-1)
            ]
            assert energy(alpha,q) == energy(1-alpha,reflect(q))
        expected_energy = {"1/4":["15/2","17/2"],"3/4":["13/2","19/2"]} if a == 1 else {
            "1/4":["7/2","25/2"],"3/4":["5/2","27/2"]
        }
        assert branches == expected_energy
        results[name] = {
            "ward_residues_at_first_endpoint": [str(x) for x in q],
            "ward_residues_at_adjoint_endpoint": [str(-x) for x in q],
            "triality": str(center(q)),
            "cover_energy_of_plus_and_minus_from_NS_vacuum": branches,
            "particle_hole_quarter_reflection_target_root": [str(x) for x in reflected_root],
        }
    p = (1,1,0,0,0,0,0,0)
    assert sum(p)%2 == sum(reflect(p))%2 == 0
    assert center(p) == 0 and center(reflect(p)) == 1
    for alpha, j in product((F(1,4),F(3,4)), range(8)):
        assert (7-j+(1-alpha))/8-F(1,2) == -((j+alpha)/8-F(1,2))
        assert energy(alpha,p) == energy(1-alpha,reflect(p))
    branch_exponents = {
        str(m): [F(2*k+1-m,2*m) for k in range(m)]
        for m in (2,3)
    }
    assert branch_exponents == {"2":[F(-1,4),F(1,4)],"3":[F(-1,3),F(0),F(1,3)]}
    for exponents in branch_exponents.values():
        assert all((x-F(1,2)).denominator != 1 for x in exponents)
    print(json.dumps({
        "marked_pairs": results,
        "literal_center_reflection_witness": {
            "NS_even_input": p, "output": reflect(p), "input_triality": "0", "output_triality": "1",
        },
        "spin_branch_exponents_in_base_half_form": {
            m:[str(x) for x in values] for m,values in branch_exponents.items()
        },
        "actual_S0_half_turn_at_regular_fixed_fibers": {"three_q_points":"+i", "three_p_points":"-i"},
        "scope": "exact original marked comparison; no primitive P1 endpoint derived",
    },indent=2))


if __name__ == "__main__":
    main()
