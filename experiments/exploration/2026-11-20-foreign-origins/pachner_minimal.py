#!/usr/bin/env python3
"""Exact finite combinatorics of the minimal asymmetric, internal Pachner move.
No TFPT constants or objects are passed as input.
"""
from itertools import combinations
from collections import Counter
import json


def move_types(dmax=9):
    return [(d, p, d + 2 - p) for d in range(1, dmax + 1)
            for p in range(1, (d + 2) // 2 + 1)]


def faces(facets, size):
    return {tuple(c) for tet in facets for c in combinations(sorted(tet), size)}


def boundary_mod2(facets):
    c = Counter(f for tet in facets for f in combinations(sorted(tet), 3))
    return {f for f, v in c.items() if v % 2 == 1}


def euler(facets):
    fv = [len(faces(facets, i)) for i in range(1, 5)]
    return fv, sum((-1)**i * n for i, n in enumerate(fv))


def oriented_boundary(facets, signs):
    c = Counter()
    for tet, coeff in zip(facets, signs):
        tet = tuple(sorted(tet))
        for i in range(4):
            face = tet[:i] + tet[i + 1:]
            c[face] += coeff * (-1) ** i
    return {f: c[f] for f in c if c[f] != 0}


def test():
    all_moves = move_types()
    eligible = [(d,p,q) for d,p,q in all_moves if p>=2 and p<q]
    smallest_sum = min(p+q for d,p,q in eligible)
    selected = [m for m in eligible if m[1]+m[2] == smallest_sum]
    assert selected == [(3,2,3)]
    # The move is the two complementary sets of tetrahedral facets of a 4-simplex.
    vertices = frozenset(range(5))
    facets = [vertices - {v} for v in sorted(vertices)]
    left, right = facets[:2], facets[2:]
    B_left, B_right = boundary_mod2(left), boundary_mod2(right)
    assert B_left == B_right and len(B_left) == 6
    assert euler(left) == ([5,9,7,2],1)
    assert euler(right) == ([5,10,9,3],1)
    # Boundary^2(4-simplex)=0, using exact signed simplicial chains over Z.
    signs = [(-1)**i for i in range(5)]
    lb = oriented_boundary(left, signs[:2])
    rb = oriented_boundary(right, signs[2:])
    assert {x: lb.get(x,0)+rb.get(x,0) for x in set(lb)|set(rb)} == {
        x:0 for x in set(lb)|set(rb)}
    # Count independent ways to orient five *labelled* facets: no distinguished 3|2 species assignment.
    labelled_splits = list(combinations(range(5),2))
    assert len(labelled_splits) == 10
    # Same-complexity null: any unordered asymmetric pair (p,q) from 2..5.
    null = [(p,q) for p in range(2,6) for q in range(p+1,6)]
    assert len(null)==6 and null.count((2,3))==1
    # Ablate the load-bearing premises to reveal the dimension 3 is not parameter-free.
    no_asym = min((m for m in all_moves if m[1]>=2), key=lambda m:m[1]+m[2])
    no_internal = min((m for m in all_moves if m[1]<m[2]), key=lambda m:m[1]+m[2])
    assert no_asym==(2,2,2) and no_internal==(1,1,2)
    return {'minimal_move': {'spatial_d':3,'before':2,'after':3,'event_spacetime_d':4},
            'four_simplex_vertices':5,'tetrahedral_facets':5,
            'boundary_triangles':len(B_left),
            'f_vector_before':euler(left)[0],'f_vector_after':euler(right)[0],
            'euler_before':1,'euler_after':1,'signed_boundary_match':True,
            'unfixed_labelled_2_plus_3_splits':len(labelled_splits),
            'null_same_grammar_successes':1,'null_same_grammar_total':len(null),
            'without_asymmetry':no_asym,'without_internality':no_internal}


if __name__ == '__main__':
    print(json.dumps(test(),indent=2,ensure_ascii=False))
