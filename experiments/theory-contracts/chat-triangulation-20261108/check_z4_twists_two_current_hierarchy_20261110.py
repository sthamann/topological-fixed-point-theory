#!/usr/bin/env python3
"""Exact tests of the marked Z4 twist and ordered E8 current paths."""
import ast
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import json
import sys

import sympy as sp


REPO = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(REPO))
from tfpt_explorer.current_block_geometry import _actual_current_source

FOLD = REPO / "experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json"
M = sp.Matrix(json.loads(FOLD.read_text())["explicit_Weyl_matrix"])
V634 = REPO / "verification/v634_st31_structure.py"
INDEX = next(ast.literal_eval(node.value) for node in ast.parse(V634.read_text()).body
             if isinstance(node, ast.Assign)
             and any(isinstance(t, ast.Name) and t.id == "SIGMA_IDX" for t in node.targets))
UNIT = sp.eye(8)
SIGMA = UNIT[INDEX, :]
Y = sp.Matrix([-2,-2,-2,3,3,0,0,0])/6
CHEVALLEY = _actual_current_source()["chevalley"]
ROOTS = set(CHEVALLEY.ridx)
ALPHAS = (sp.Rational(1,4),sp.Rational(3,4))


def v2(v):
    return tuple(2*x for x in v)


def root_analysis(t):
    fixed = {r for r in ROOTS if sum(x*y for x,y in zip(t,r,strict=True))/2 % 1 == 0}
    components = []
    while fixed:
        block = {fixed.pop()}
        while True:
            neighbours = {r for r in fixed if any(sum(x*y for x,y in zip(r,b)) != 0 for b in block)}
            if not neighbours:
                break
            fixed -= neighbours
            block |= neighbours
        components.append((sp.Matrix(list(block)).rank(),len(block)))
    return sorted(components)


def closest_coset_weights(t):
    results = []
    for lattice_parity in (0,1):
        choices = []
        for value in t:
            v = value + F(lattice_parity,2)
            floor = v.numerator//v.denominator
            choices.append([(j, v-j) for j in (floor-1,floor,floor+1)])
        best = min(sum(x*x for _,x in c)/2 for c in product(*choices)
                   if sum(j for j,_ in c)%2 == 0)
        results.append(best)
    return min(results)


def z4():
    actual = (F(1,2),)*5 + (F(0),)*3
    assert root_analysis(actual) == [(3,12),(5,40)]
    cases = {}
    for label, t in (("uniform_quarter", (F(1,4),)*8),
                     ("odd_signed_quarter", (F(-1,4),)+(F(1,4),)*7),
                     ("native_glue", actual)):
        weights = [closest_coset_weights([k*x for x in t]) for k in (1,2,3)]
        cases[label] = {"fixed_root_components_rank_count":root_analysis(t),
                        "E8_shifted_module_ground_weights":[str(x) for x in weights]}
    assert cases["uniform_quarter"]["fixed_root_components_rank_count"] == [(1,2),(7,126)]
    assert cases["odd_signed_quarter"]["fixed_root_components_rank_count"] == [(7,56)]
    assert cases["native_glue"]["E8_shifted_module_ground_weights"] == ["3/8","1/2","3/8"]
    native_type = 16*F(3,8)%4
    assert native_type == 2
    isotropic = {(a,b) for a,b in product(range(4),repeat=2) if (5*a*a+3*b*b)%8 == 0}
    plus = {(a,a) for a in range(4)}
    minus = {(a,-a%4) for a in range(4)}
    assert isotropic == plus | minus
    assert {a-b for a,b in plus} == {0}
    assert {(a-b)%4 for a,b in minus} == {0,2}
    for alpha in (F(1,4),F(3,4)):
        t = tuple(F(k,2)+alpha/2 for k in range(8))
        theta = tuple(x-((x+F(1,2))//1) for x in t)
        assert sum(x*x for x in theta)/2 == F(5,16)
        weights = [closest_coset_weights([k*x for x in t]) for k in (1,2,3)]
        assert weights == [F(5,16),F(1,4),F(5,16)]
        cases[f"charged_fold_C4_alpha_{alpha}"] = {
            "complex_eigenangles":sorted(str(x) for x in theta),
            "vector_action_order":8,
            "E8_root_action_order":4,
            "E8_twist_type":str(16*weights[0]%4),
            "E8_fixed_root_components_rank_count":root_analysis(t),
            "E8_shifted_module_ground_weights":[str(x) for x in weights],
        }
    assert cases["uniform_quarter"]["E8_shifted_module_ground_weights"] == ["1/4","0","1/4"]
    assert cases["odd_signed_quarter"]["E8_shifted_module_ground_weights"] == ["1/4","1/2","1/4"]
    cases["native_glue"]["E8_twist_type"] = int(native_type)
    cases["native_glue"]["holomorphic_transversal_count"] = 0
    return cases


def bracket(a,b,target):
    result = CHEVALLEY.bracket(CHEVALLEY.ridx[v2(a)],CHEVALLEY.ridx[v2(b)])
    assert set(result) == {CHEVALLEY.ridx[v2(target)]}
    return result[CHEVALLEY.ridx[v2(target)]]


def source_sector(v):
    q = M.T*v
    if all(x.is_integer for x in q):
        return "0" if sum(q)%2 == 0 else "v"
    return "s" if sum(q-sp.ones(8,1)/2)%2 == 0 else "c"


def grade(v):
    return 2*sum(v[:5])%4


def currents():
    all_rows = []
    for s,t in product((1,-1),repeat=2):
        beta = s*UNIT[:,2]+t*UNIT[:,3]
        p = SIGMA*beta
        a = t*(UNIT[:,3]-UNIT[:,5])
        b = s*(UNIT[:,2]-UNIT[:,4])
        weak_first, color_first = p+a,p+b
        assert p == s*UNIT[:,4]+t*UNIT[:,5]
        assert weak_first+b == color_first+a == beta
        assert not CHEVALLEY.bracket(CHEVALLEY.ridx[v2(a)], CHEVALLEY.ridx[v2(b)])
        assert bracket(a,p,weak_first)*bracket(b,weak_first,beta) == 1
        assert bracket(b,p,color_first)*bracket(a,color_first,beta) == 1
        assert [source_sector(x) for x in (a,b)] == ["s","0"]
        assert [str(grade(x)) for x in (a,b)] == ["2","0"]
        states = (beta,p,weak_first,color_first)
        assert [source_sector(x) for x in states] == ["0","s","0","s"]
        assert all(q.dot(q)==2 for q in states)
        assert Y.dot(weak_first)==sp.Rational(s+t,2)
        assert Y.dot(color_first)==-sp.Rational(s,3)
        assert weak_first.dot(UNIT[:,3]-UNIT[:,4]) == t-s
        assert color_first.dot(UNIT[:,3]-UNIT[:,4]) == 0
        gaps = []
        k_orders = {}
        for alpha in ALPHAS:
            d = M*sp.Matrix([(j+alpha)/8-sp.Rational(1,2) for j in range(8)])
            assert len(set(d[:3]))==3 and len(set(d[3:5]))==2
            gaps.append(str(d.dot(weak_first-color_first)))
            k_orders[str(alpha)] = [str(1+d.dot(q)) for q in states]
        all_rows.append({
            "sign_s_t":[s,t],
            "states_in_order_beta_clock_weakfirst_colorfirst":[list(map(str,x)) for x in states],
            "target_glue_grades":[str(grade(x)) for x in states],
            "source_sectors":[source_sector(x) for x in states],
            "intermediate_hypercharge":[str(Y.dot(x)) for x in (weak_first,color_first)],
            "intermediate_weak_spin":[0 if s==t else 1,0],
            "intermediate_color_representation":["1","3" if s==1 else "bar3"],
            "quarter_K_levels":k_orders,
            "between_window_K_weak_minus_color":gaps,
            "same_Chevalley_coefficient":"+1",
        })
    assert [r["between_window_K_weak_minus_color"] for r in all_rows] == [
        ["-3/8","-5/8"], ["1/4","1/4"], ["-1/4","-1/4"], ["3/8","5/8"]
    ]
    return all_rows


def main():
    out = {"twists":z4(),"current_paths":currents()}
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
