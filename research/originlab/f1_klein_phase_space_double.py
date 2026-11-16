"""Exact Klein/dual F2^4 five-point bridge and S5-to-S4 anchor obstruction.

This mathematical construction assumes the 4-element Klein event group
G=F2^2, its dual, and an additional odd polarization. It does not
derive the physical original R113 Q0 or distinguish the physical P1.
"""
from itertools import permutations,product
import json
import sympy as sp

def require(ok,msg):
    if not ok: raise ArithmeticError(msg)

V=list(product((0,1),repeat=4))

def q(v):
    a,b,c,d=v
    return (a*c+b*d+a+c)%2

def B(x,y):
    return (x[0]*y[2]+x[2]*y[0]+x[1]*y[3]+x[3]*y[1])%2

def controls():
    zeros=[x for x in V if any(x) and q(x)==0]
    require(zeros==[(0,0,0,1),(0,1,0,0),(0,1,1,1),
                    (1,1,0,1),(1,1,1,1)],
            "five-point quadratic zero set changed")
    require(all(B(x,y)==1 for x,y in product(zeros,repeat=2) if x!=y)
            and all(sum(v[i] for v in zeros)%2==0 for i in range(4)),
            "five q-null points are not a pairwise symplectic-one simplex")
    E=sp.Matrix.hstack(*(sp.Matrix(x) for x in zeros[:4]))
    require(abs(E.det())==1,"four null points must span the source")
    coords={v:tuple(int(t)%2 for t in E.inv()*sp.Matrix(v)) for v in V}
    maps=set()
    for p in permutations(range(5)):
        def send(v):
            c=coords[v]
            return tuple(sum(c[j]*zeros[p[j]][i] for j in range(4))%2
                         for i in range(4))
        require(send(zeros[4])==zeros[p[4]]
                and all(q(v)==q(send(v)) for v in V)
                and all(B(x,y)==B(send(x),send(y))
                        for x,y in product(V,repeat=2)),
                "one full S5 permutation failed to preserve q and B")
        maps.add(tuple(send(v) for v in V))
    require(len(maps)==120
            and sum(p[4]==4 for p in permutations(range(5)))==24
            and sum(q(v)==0 for v in V)==6,
            "odd O(4,2) isometries or S4 stabilizer census changed")
    return {
        "status":"EXACT_KLEIN_DOUBLE_CONTROL_P1_SELECTION_OPEN",
        "source_double":"G x G*=F2^4 with alternating symplectic B",
        "new_odd_quadratic_refinement":"q(g,chi)=g.chi+g1+chi1",
        "extra_spin_polarization_derived":False,
        "nonzero_null_point_count":5,
        "sum_of_five_null_points":"0000",
        "q0_nonzero_points":[list(x) for x in zeros],
        "pairwise_symplectic_one":True,
        "odd_q0_count_including_zero":6,
        "full_isometry_count":len(maps),
        "full_isometry_action":"O-(4,2)=S5, transitive on all five null points",
        "one_named_point_stabilizer_order":24,
        "one_named_point_stabilizer":"S4",
        "fifth_P1_point_selected_by_full_symmetry":False,
        "original_R113_source_or_physical_symmetry_breaking_found":False,
        "F1_solved":False
    }

if __name__=="__main__":
    print(json.dumps(controls(),indent=2,sort_keys=True))
