#!/usr/bin/env python3
"""Exact affine and graph audit for the equator lift on the square elliptic curve.

The group is generated on C/(Z+iZ) by
    Q(z)=i z,
    r(z)=i z+1/2,
    kappa(z)=-i conjugate(z)+(1+i)/2.

Group closure, element orders and fixed-component counts are computed from
the affine maps.  No real-form census table is supplied to the calculation.
"""

from __future__ import annotations

from collections import Counter, deque
from fractions import Fraction as F
from functools import reduce
from math import gcd
from pathlib import Path
import hashlib
import json
import sympy as sp


Matrix = tuple[tuple[int, int], tuple[int, int]]
Vector = tuple[F, F]
Affine = tuple[Matrix, Vector]

I: Matrix = ((1, 0), (0, 1))
J: Matrix = ((0, -1), (1, 0))
K: Matrix = ((0, -1), (-1, 0))
ZERO: Vector = (F(0), F(0))
T: Vector = (F(1, 2), F(1, 2))
IDENTITY: Affine = (I, ZERO)


def mm(a: Matrix, b: Matrix) -> Matrix:
    return tuple(
        tuple(sum(a[i][k] * b[k][j] for k in range(2)) for j in range(2))
        for i in range(2)
    )  # type: ignore[return-value]


def mv(a: Matrix, v: Vector) -> Vector:
    return (
        a[0][0] * v[0] + a[0][1] * v[1],
        a[1][0] * v[0] + a[1][1] * v[1],
    )


def add(v: Vector, w: Vector) -> Vector:
    return v[0] + w[0], v[1] + w[1]


def canonical(f: Affine) -> Affine:
    matrix, shift = f
    return matrix, (shift[0] % 1, shift[1] % 1)


def compose(f: Affine, g: Affine) -> Affine:
    """Return f after g, modulo the square lattice."""
    mf, df = f
    mg, dg = g
    return canonical((mm(mf, mg), add(mv(mf, dg), df)))


def power(f: Affine, n: int) -> Affine:
    out = IDENTITY
    for _ in range(n):
        out = compose(f, out)
    return out


def closure(generators: tuple[Affine, ...]) -> set[Affine]:
    group = {IDENTITY}
    queue = deque([IDENTITY])
    while queue:
        a = queue.popleft()
        for b in generators:
            for candidate in (compose(a, b), compose(b, a)):
                if candidate not in group:
                    group.add(candidate)
                    queue.append(candidate)
    return group


def determinant(matrix: Matrix) -> int:
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def element_order(f: Affine, limit: int = 16) -> int:
    for n in range(1, limit + 1):
        if power(f, n) == IDENTITY:
            return n
    raise AssertionError(f"order exceeds {limit}: {f}")


def fixed_components_of_real_involution(f: Affine) -> int:
    """Count circles in Fix(f) by the rank-one Smith invariant of I-M.

    For (I-M)v=d mod Z^2, a primitive integer left-null vector tests
    solvability.  If solvable, gcd of the entries of I-M is the component
    group order of its torus kernel, hence the number of fixed circles.
    """
    matrix, shift = f
    if determinant(matrix) != -1 or element_order(f) != 2:
        raise ValueError("fixed-circle classification requires a real involution")
    a = (
        (1 - matrix[0][0], -matrix[0][1]),
        (-matrix[1][0], 1 - matrix[1][1]),
    )
    nonzero_column = next(
        (a[0][j], a[1][j]) for j in range(2) if (a[0][j], a[1][j]) != (0, 0)
    )
    divisor = gcd(abs(nonzero_column[0]), abs(nonzero_column[1]))
    left_null = (nonzero_column[1] // divisor, -nonzero_column[0] // divisor)
    obstruction = left_null[0] * shift[0] + left_null[1] * shift[1]
    if obstruction.denominator != 1:
        return 0
    entries = [abs(x) for row in a for x in row if x]
    return reduce(gcd, entries)


def real_type(f: Affine) -> str:
    order = element_order(f)
    if order != 2:
        return f"anti_order_{order}"
    components = fixed_components_of_real_involution(f)
    return {0: "empty", 1: "nondividing", 2: "dividing"}[components]


def affine_row(f: Affine) -> dict[str, object]:
    matrix, shift = f
    return {
        "linear": [list(matrix[0]), list(matrix[1])],
        "shift": [str(shift[0]), str(shift[1])],
        "order": element_order(f),
        "type": real_type(f) if determinant(matrix) == -1 else "holomorphic",
    }


def fixes(f: Affine, point: Vector) -> bool:
    matrix, shift = f
    return canonical((I, add(mv(matrix, point), shift))) == canonical((I, point))


def rational_rank(matrix: list[list[int]]) -> int:
    a = [[F(x) for x in row] for row in matrix]
    rows = len(a)
    cols = len(a[0]) if rows else 0
    rank = 0
    for col in range(cols):
        pivot = next((r for r in range(rank, rows) if a[r][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        p = a[rank][col]
        a[rank] = [x / p for x in a[rank]]
        for r in range(rows):
            if r != rank and a[r][col]:
                c = a[r][col]
                a[r] = [x - c * y for x, y in zip(a[r], a[rank])]
        rank += 1
    return rank


def main() -> None:
    q: Affine = (J, ZERO)
    r: Affine = (J, (F(1, 2), F(0)))
    kappa: Affine = (K, T)
    h = power(r, 2)
    seam_s = compose(r, kappa)
    seam_hs = compose(h, seam_s)

    old_holomorphic = closure((q, h))
    old_total = closure((q, h, kappa))
    full_holomorphic = closure((q, r))
    full_total = closure((q, r, kappa))

    old_anti = {f for f in old_total if determinant(f[0]) == -1}
    full_anti = {f for f in full_total if determinant(f[0]) == -1}
    new_anti = full_anti - old_anti

    old_census = Counter(real_type(f) for f in old_anti)
    full_census = Counter(real_type(f) for f in full_anti)
    new_census = Counter(real_type(f) for f in new_anti)

    # The old holomorphic order-eight group is abelian C4 x C2, not D4.
    old_commutes = all(compose(a, b) == compose(b, a) for a in old_holomorphic for b in old_holomorphic)
    old_orders = Counter(element_order(f) for f in old_holomorphic)
    h_central = all(compose(h, a) == compose(a, h) for a in full_holomorphic)
    lifted_dihedral_relation = compose(r, compose(q, power(r, 3))) == compose(h, power(q, 3))

    # S fixes y=1/4,3/4; hS fixes x=1/4,3/4.  Construct the
    # intersection graph and its incidence matrix rather than inserting V,E,b1.
    horizontal = (F(1, 4), F(3, 4))
    vertical = (F(1, 4), F(3, 4))
    vertices = tuple((x, y) for x in vertical for y in horizontal)
    vertex_index = {v: i for i, v in enumerate(vertices)}
    edges: list[tuple[int, int, str]] = []
    for y in horizontal:
        endpoints = [vertex_index[(x, y)] for x in vertical]
        edges.extend((endpoints[0], endpoints[1], f"horizontal_y={y}:{arc}") for arc in range(2))
    for x in vertical:
        endpoints = [vertex_index[(x, y)] for y in horizontal]
        edges.extend((endpoints[0], endpoints[1], f"vertical_x={x}:{arc}") for arc in range(2))

    incidence = [[0 for _ in edges] for _ in vertices]
    adjacency = {i: set() for i in range(len(vertices))}
    for column, (source, target, _) in enumerate(edges):
        incidence[source][column] = -1
        incidence[target][column] = 1
        adjacency[source].add(target)
        adjacency[target].add(source)
    reached = {0}
    queue = deque([0])
    while queue:
        v = queue.popleft()
        for w in adjacency[v] - reached:
            reached.add(w)
            queue.append(w)
    incidence_rank = rational_rank(incidence)
    graph_b1 = len(edges) - incidence_rank
    # The two horizontal and two vertical circles cut the two cyclic coordinate
    # intervals independently, giving their Cartesian-product complement cells.
    graph_faces = len(horizontal) * len(vertical)
    torus_euler_faces = len(edges) - len(vertices)

    # The four graph vertices are exactly the fixed points of h: 2z=T mod Lambda.
    h_fixed = {
        ((T[0] + m) / 2 % 1, (T[1] + n) / 2 % 1)
        for m in (0, 1) for n in (0, 1)
    }
    old_dividing = {f for f in old_anti if real_type(f) == "dividing"}
    new_dividing = {f for f in new_anti if real_type(f) == "dividing"}
    old_dividing_branch_counts = sorted(sum(fixes(f, v) for v in h_fixed) for f in old_dividing)
    new_dividing_branch_counts = sorted(sum(fixes(f, v) for v in h_fixed) for f in new_dividing)

    # Algebraic model E: y^2=x^4-1.  Clear the x denominator after applying
    # r(x,y)=(1/x,i*y/x^2); the transformed equation is minus the original.
    x, y = sp.symbols("x y", nonzero=True)
    original_curve = y**2 - x**4 + 1
    transformed_curve_cleared = sp.expand(
        x**4 * ((sp.I * y / x**2) ** 2 - (1 / x) ** 4 + 1)
    )
    algebraic_r_preserves = sp.simplify(transformed_curve_cleared + original_curve) == 0
    r2_x = sp.simplify(1 / (1 / x))
    r2_y = sp.simplify(sp.I * (sp.I * y / x**2) / (1 / x) ** 2)

    checks = {
        "algebraic_r_preserves_curve_and_squares_to_h": bool(
            algebraic_r_preserves and r2_x == x and r2_y == -y
        ),
        "r_affine_square_is_h": h == (((-1, 0), (0, -1)), T),
        "S_is_conjugation_plus_i_over_2": seam_s == (((1, 0), (0, -1)), (F(0), F(1, 2))),
        "hS_is_minus_conjugation_plus_1_over_2": seam_hs == (((-1, 0), (0, 1)), (F(1, 2), F(0))),
        "old_holomorphic_group_order_8": len(old_holomorphic) == 8,
        "old_holomorphic_group_is_C4_times_C2_not_D4": old_commutes and old_orders == Counter({4: 4, 2: 3, 1: 1}),
        "old_total_group_order_16_with_8_anti": len(old_total) == 16 and len(old_anti) == 8,
        "full_holomorphic_group_order_16": len(full_holomorphic) == 16,
        "full_holomorphic_is_central_deck_extension_of_base_D4": h_central and lifted_dihedral_relation,
        "full_total_group_order_32_with_16_anti": len(full_total) == 32 and len(full_anti) == 16,
        "old_real_census_computed": old_census == Counter({"nondividing": 4, "empty": 2, "dividing": 2}),
        "new_coset_census_computed": new_census == Counter({"anti_order_4": 4, "empty": 2, "dividing": 2}),
        "full_real_census_computed": full_census == Counter({"nondividing": 4, "empty": 4, "dividing": 4, "anti_order_4": 4}),
        "old_dividing_OS_forms_fix_no_branch_marks": old_dividing_branch_counts == [0, 0],
        "new_dividing_forms_are_exactly_S_and_hS": new_dividing == {seam_s, seam_hs},
        "S_and_hS_each_fix_all_four_branch_marks": new_dividing_branch_counts == [4, 4],
        "seam_graph_is_connected": len(reached) == len(vertices),
        "seam_graph_V4_E8": len(vertices) == 4 and len(edges) == 8,
        "seam_graph_F4_on_torus": graph_faces == 4 and graph_faces == torus_euler_faces,
        "seam_graph_b1_5_from_incidence_rank": incidence_rank == 3 and graph_b1 == 5,
        "graph_vertices_are_h_branch_points": set(vertices) == h_fixed,
    }
    if not all(checks.values()):
        raise AssertionError([name for name, ok in checks.items() if not ok])

    result = {
        "status": "PASS",
        "checks": checks,
        "group_sizes": {
            "old_holomorphic": len(old_holomorphic),
            "old_total": len(old_total),
            "full_holomorphic": len(full_holomorphic),
            "full_total": len(full_total),
        },
        "real_census": {
            "old_eight": dict(sorted(old_census.items())),
            "new_eight": dict(sorted(new_census.items())),
            "full_sixteen": dict(sorted(full_census.items())),
        },
        "branch_mark_counts_on_dividing_forms": {
            "old_two": old_dividing_branch_counts,
            "new_S_and_hS": new_dividing_branch_counts,
        },
        "old_real_maps": sorted((affine_row(f) for f in old_anti), key=lambda row: json.dumps(row, sort_keys=True)),
        "new_real_maps": sorted((affine_row(f) for f in new_anti), key=lambda row: json.dumps(row, sort_keys=True)),
        "seam_graph": {
            "vertices": [[str(x), str(y)] for x, y in vertices],
            "edges": [[source, target, label] for source, target, label in edges],
            "V": len(vertices),
            "E": len(edges),
            "F": graph_faces,
            "incidence_rank": incidence_rank,
            "b1": graph_b1,
        },
        "scope": (
            "Exact affine square-torus group and equator-preimage graph. "
            "This does not identify a chosen oval or smoothing with the original physical P1 seam state."
        ),
    }
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    output = Path(__file__).with_name("independent-seam-graph-results.json")
    output.write_text(payload)
    print(f"PASS {len(checks)}/{len(checks)}")
    print(f"RESULT_SHA256 {hashlib.sha256(payload.encode()).hexdigest()}")
    print(f"WROTE {output}")


if __name__ == "__main__":
    main()
