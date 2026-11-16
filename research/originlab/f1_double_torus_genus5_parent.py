"""Conditional pre-E8 source spectroscopy on overlays of elliptic double covers.

The genus-five lattice theta separation is a PUBLISHED theorem, not calculated
here. This module independently checks branch geometry and 16D root structure.
No physical original TFPT Q0, OS time, or F1 origin is derived.
"""
from __future__ import annotations
from itertools import combinations, product


def overlay_genus(branch_a: tuple[int, ...], branch_b: tuple[int, ...]) -> int | None:
    """Normalized connected V4 cover of P1 from 2 elliptic double covers.

    Return None when covers coincide, yielding disconnected copies.
    Both branch sets must contain four distinct finite points.
    """
    if len(branch_a) != 4 or len(branch_b) != 4 or len(set(branch_a)) != 4 or len(set(branch_b)) != 4:
        raise ValueError("Each elliptic cover needs four distinct branch values")
    intersection = len(set(branch_a) & set(branch_b))
    return None if intersection == 4 else 5 - intersection


def projection_ramification_count(branch_a: tuple[int, ...], branch_b: tuple[int, ...]) -> int:
    """Number of branch points of C->E_a, excluding identical covers."""
    if overlay_genus(branch_a, branch_b) is None:
        raise ValueError("Coincident covers do not give connected C")
    return 2 * len(set(branch_b) - set(branch_a))


def castelnuovo_severi_trigonal_bound() -> int:
    """Genus of a curve with independent maps degree2 -> E, degree3 -> P1.

    Since gcd(2,3)=1, the induced map to E x P1 is generically
    birational; C-S gives g <= 2*1 + 3*0 + (2-1)*(3-1) = 4.
    """
    return 4


def e8_roots_doubled() -> tuple[tuple[int, ...], ...]:
    """240 norm-two E8 roots in doubled integer coordinates."""
    roots: set[tuple[int, ...]] = set()
    for i,j in combinations(range(8),2):
        for si,sj in product((-2,2), repeat=2):
            v=[0]*8
            v[i]=si
            v[j]=sj
            roots.add(tuple(v))
    for signs in product((-1,1),repeat=8):
        if sum(x<0 for x in signs)%2==0:
            roots.add(signs)
    assert len(roots)==240
    assert all(sum(x*x for x in v)==8 for v in roots)
    return tuple(sorted(roots))


def rank16_roots(parent: str) -> tuple[tuple[int, ...], ...]:
    """All norm-two vectors of E8+E8 or D16+, doubled integer coords."""
    if parent == "E8xE8":
        e=e8_roots_doubled()
        z=(0,)*8
        return tuple(sorted(tuple(r+z for r in e)+tuple(z+r for r in e)))
    if parent == "D16plus":
        roots=[]
        for i,j in combinations(range(16),2):
            for si,sj in product((-2,2),repeat=2):
                v=[0]*16
                v[i]=si
                v[j]=sj
                roots.append(tuple(v))
        assert len(roots)==480
        return tuple(roots)
    raise ValueError("Unknown parent lattice")


def root_graph_components(roots: tuple[tuple[int,...],...]) -> tuple[int,...]:
    """Connected components of graph with edges for nonorthogonal roots."""
    todo=set(range(len(roots)))
    components=[]
    while todo:
        seed=todo.pop()
        active=[seed]
        count=0
        while active:
            i=active.pop()
            count+=1
            neighbours=[j for j in todo if any(a*b for a,b in zip(roots[i], roots[j]))
                        and sum(a*b for a,b in zip(roots[i],roots[j]))!=0]
            todo.difference_update(neighbours)
            active.extend(neighbours)
        components.append(count)
    return tuple(sorted(components))


def controls() -> dict:
    a=(-4,-3,-2,-1)
    b=(1,2,3,4)
    return {
        "disjoint_branch_overlay_genus":overlay_genus(a,b),
        "projection_ramification_count":projection_ramification_count(a,b),
        "one_shared_branch_overlay_genus":overlay_genus(a,(-1,1,2,3)),
        "trigonal_bielliptic_genus_upper_bound":castelnuovo_severi_trigonal_bound(),
        "E8xE8_roots":len(rank16_roots("E8xE8")),
        "D16plus_roots":len(rank16_roots("D16plus")),
        "E8xE8_root_components":root_graph_components(rank16_roots("E8xE8")),
        "D16plus_root_components":root_graph_components(rank16_roots("D16plus")),
        "published_genus5_partition_difference_evaluated_here":False,
        "original_Q0_found":False,
        "F1_solved":False,
    }


if __name__ == "__main__":
    import json
    print(json.dumps(controls(),indent=2))
