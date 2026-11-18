#!/usr/bin/env python3
"""Exact source root, common center/time and first missing Ramond operation."""
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import json
import sys

import sympy as sp

ROOT = Path(sys.argv[1]).resolve()
sys.path.insert(0,str(ROOT))
from tfpt_explorer.current_block_geometry import _actual_current_source

fold = ROOT / "experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json"
M = sp.Matrix(json.loads(fold.read_text())["explicit_Weyl_matrix"])
CHEV = _actual_current_source()["chevalley"]
E = sp.eye(8)
ONE = sp.ones(8,1)
W = sp.Matrix([1,-1,0,-1,-1,-1,0,-1])
Y = sp.Matrix([-2,-2,-2,3,3,0,0,0])/6
QUARTERS = (sp.Rational(1,4), sp.Rational(3,4))


def drift(a,p):
    return sp.Matrix([(j+a)/8-sp.Rational(1,2) for j in range(8)]).dot(p)


def target_bracket(x,y):
    x2, y2 = tuple(2*M*x),tuple(2*M*y)
    target = tuple(2*M*(x+y))
    answer = CHEV.bracket(CHEV.ridx[x2],CHEV.ridx[y2])
    assert set(answer) == {CHEV.ridx[target]}
    return int(answer[CHEV.ridx[target]])


def encode(p):
    return {"source_charge":list(map(str,p)), "target":list(map(str,M*p)),
            "Q":str(ONE.dot(p)), "Y":str(Y.dot(M*p)),
            "color_center":str(W.dot(p)%3),
            "quarter_K_shift":[str(drift(a,p)) for a in QUARTERS]}


def neutral_seed():
    q = M.T*(E[:,7]-E[:,6])
    spinors = [sp.Matrix(s)/2 for s in product((-1,1),repeat=8) if s.count(-1)%2==0]
    passing = {tuple(x) for x in spinors if ONE.dot(x)==0 and W.dot(x)%3==0
               and all(drift(a,x).is_integer for a in QUARTERS)}
    assert passing == {tuple(q),tuple(-q)}
    assert q.dot(q)==2 and all(drift(a,q)==0 for a in QUARTERS)
    p = E[:,3]-E[:,1]
    assert p.dot(p)==2 and p.dot(q)==-1 and (p+q).dot(p+q)==2
    for x in (p,q,p+q):
        assert W.dot(x)%3==0 and ONE.dot(x)==0
    assert all(1+drift(a,p)==1+drift(a,p+q)==sp.Rational(5,4) for a in QUARTERS)
    return q, {"root":encode(q),"opposite":encode(-q), "number_of_center_neutral_Q_zero_integer_drift_spinors":len(passing),
               "excited_input":encode(p),"excited_output":encode(p+q),
               "transition_K":"5/4","Chevalley_coefficient":target_bracket(q,p)}


def roots_to_A(q):
    A = M.T*(E[:,3]-E[:,5])
    diff = A-q
    nonzero = [i for i,x in enumerate(diff) if x != 0]
    assert len(nonzero)==4 and diff.dot(diff)==4
    paths = []
    for pair in combinations(nonzero,2):
        first = sp.Matrix([x if i in pair else 0 for i,x in enumerate(diff)])
        second = diff-first
        assert first.dot(first)==second.dot(second)==2
        assert first.dot(second)==0
        assert all(x.is_integer for p in (first,second) for x in p)
        assert (q+first).dot(q+first)==2
        sign1,sign2=target_bracket(first,q),target_bracket(second,q+first)
        paths.append({"first":encode(first), "second":encode(second),
                      "middle":encode(q+first), "coefficient":sign1*sign2,
                      "center_neutral_both":W.dot(first)%3==W.dot(second)%3==0})
    center_neutral = [x for x in paths if x["center_neutral_both"]]
    assert len(paths)==6 and len(center_neutral)==2
    assert len({x["coefficient"] for x in center_neutral})==1
    assert {tuple(x["first"]["source_charge"]) for x in center_neutral}=={
        tuple(map(str,-E[:,2]-E[:,6])),tuple(map(str,E[:,3]-E[:,5]))
    }
    assert W.dot(A)%3==0 and ONE.dot(A)==-2
    assert all(not drift(a,A).is_integer for a in QUARTERS)
    return {"minimal_integral_root_steps":2,"all_ordered_shortest_paths":paths,
            "center_neutral_order_count":2,"A":encode(A)}


def old_routes():
    B = M.T*(E[:,2]-E[:,4])
    weak = M.T*(E[:,3]-E[:,4])
    assert B == E[:,3]+E[:,5] and W.dot(B)%3==1
    assert weak == E[:,3]-E[:,7]
    b = E[:,2]+E[:,3]
    p = E[:,4]+E[:,5]
    assert W.dot(M.T*p)%3==0 and W.dot(M.T*b)%3==1
    # Old neutral time-preserving global lift carries center-changing roots.
    had = [(1,1,1,1),(1,1,-1,-1),(1,-1,1,-1),(1,-1,-1,1)]
    old_r = [sp.Matrix(list(h)+[-x for x in h])/2 for h in had]
    assert all(W.dot(r)%3==1 for r in old_r)
    # Universal scalar insertion cannot alter a normalized source derivative.
    x,w0,w1,u0,u1,v0,v1=sp.symbols("x w0 w1 u0 u1 v0 v1")
    scale,word,empty=w0+w1*x,u0+u1*x,v0+v1*x
    assert sp.simplify(sp.diff((scale*word)/(scale*empty),x)-sp.diff(word/empty,x))==0
    return {"B":encode(B), "actual_weak_root":encode(weak), "old_affine_lift_R_roots":[encode(r) for r in old_r],
            "raw_p_to_beta_center_change":"1", "common_scalar_normalized_source_derivative":"exact cancellation"}


def main():
    q, seed = neutral_seed()
    data={"neutral_time_preserving_excited_seed":seed,"minimal_center_neutral_generation":roots_to_A(q),
          "source_and_weight_controls":old_routes()}
    print(json.dumps(data, indent=2))


if __name__ == "__main__":
    main()
