#!/usr/bin/env python3
"""Exact finite geometry bridge from the derived five-site perfect code.

Enumerate GL(4,2) and automorphisms of the five syndrome lines. This only
establishes a syndrome-incidence symmetry, not its unitary action on matter.
"""

from itertools import permutations
from collections import Counter
from functools import reduce
from pathlib import Path
import json
from qec_origin import word_pauli, syndrome, span, permute_pauli, gf2_rank


def compose(g, h):
    return tuple(g[h[x]] for x in range(len(g)))


def inverse(g):
    result = [0] * len(g)
    for i, value in enumerate(g):
        result[value] = i
    return tuple(result)


def order(g):
    identity = tuple(range(len(g)))
    product = identity
    for exponent in range(1, 361):
        product = compose(product, g)
        if product == identity:
            return exponent
    raise AssertionError("Order out of declared finite bound")


def generated_by(generators, degree):
    identity = tuple(range(degree))
    members, todo = {identity}, [identity]
    while todo:
        g = todo.pop()
        for h in generators:
            gh = compose(g, h)
            if gh not in members:
                members.add(gh)
                todo.append(gh)
    return members


def derived_subgroup(group):
    inverses = {g: inverse(g) for g in group}
    commutators = {compose(compose(compose(g, h), inverses[g]), inverses[h])
                   for g in group for h in group}
    return generated_by(commutators, len(next(iter(group))))


def invertible_binary_maps():
    for a in range(1, 16):
        for b in range(1, 16):
            span2 = {0, a, b, a ^ b}
            if len(span2) != 4:
                continue
            for c in range(1, 16):
                if c in span2:
                    continue
                span3 = span2 | {v ^ c for v in span2}
                for d in range(1, 16):
                    if d in span3:
                        continue
                    columns = (a, b, c, d)
                    yield tuple(reduce(int.__xor__,
                                       (columns[i] for i in range(4) if (v >> i) & 1), 0)
                                for v in range(16))


def orbit_partition(group):
    remaining = set(range(len(next(iter(group)))))
    result = []
    while remaining:
        point = min(remaining)
        orbit = {g[point] for g in group}
        result.append(sorted(orbit))
        remaining -= orbit
    return result


def deleted_even_subset_action(permutation, value):
    # Coordinates e_i+e_4, i=0..3, for the even-weight subspace of F2^5.
    five_bits = value | ((value.bit_count() & 1) << 4)
    permuted = sum(((five_bits >> i) & 1) << permutation[i] for i in range(5))
    return permuted & 15


def even_exterior_trace(permutation):
    # det(I+tP) = product_cycles(1 - (-t)^cycle_length).
    unseen = set(range(len(permutation)))
    polynomial = [1]
    while unseen:
        first = next(iter(unseen))
        point, length = first, 0
        while True:
            unseen.remove(point)
            length += 1
            point = permutation[point]
            if point == first:
                break
        result = [0] * (len(polynomial) + length)
        coefficient = -((-1) ** length)
        for i, value in enumerate(polynomial):
            result[i] += value
            result[i + length] += coefficient * value
        polynomial = result
    return sum(polynomial[::2])


def main():
    rows = [word_pauli(word) for word in ("XZZXI", "IXZZX", "XIXZZ", "ZXIXZ")]
    lines = [frozenset(syndrome(e, rows, 5)
                       for e in (1 << i, 1 << (i + 5), (1 << i) | (1 << (i + 5))))
             for i in range(5)]
    assert set.union(*(set(line) for line in lines)) == set(range(1, 16))
    assert all(len(line) == 3 and reduce(int.__xor__, line) == 0 for line in lines)
    assert len(set.union(*(set(line) for line in lines))) == sum(map(len, lines))

    automorphisms, image, kernel = set(), set(), set()
    image_of = {}
    census = 0
    for mapping in invertible_binary_maps():
        census += 1
        mapped_lines = [frozenset(mapping[v] for v in line) for line in lines]
        if set(mapped_lines) == set(lines):
            permutation = tuple(lines.index(line) for line in mapped_lines)
            automorphisms.add(mapping)
            image.add(permutation)
            image_of[mapping] = permutation
            if permutation == tuple(range(5)):
                kernel.add(mapping)
    assert census == 20160
    assert (len(automorphisms), len(image), len(kernel)) == (360, 120, 3)
    assert image == set(permutations(range(5)))
    assert sorted(order(g) for g in kernel) == [1, 3, 3]

    derived1 = derived_subgroup(automorphisms)
    derived2 = derived_subgroup(derived1)
    derived3 = derived_subgroup(derived2)
    assert [len(group) for group in (automorphisms, derived1, derived2, derived3)] == [360, 180, 60, 60]
    perfect_core_image = {image_of[g] for g in derived2}
    assert len(perfect_core_image) == 60
    is_even = lambda p: sum(p[i] > p[j] for i in range(5) for j in range(i + 1, 5)) % 2 == 0
    assert perfect_core_image == {p for p in permutations(range(5)) if is_even(p)}
    assert derived2 & kernel == {tuple(range(16))}

    code_span = span(rows)
    physical_permutations = {
        perm for perm in permutations(range(5))
        if {permute_pauli(value, perm) for value in code_span} == code_span
    }
    assert len(physical_permutations) == 10

    # No canonical marked subset of size two or three exists for S5,
    # for the canonical A5 image, or already for physical C5/D5.
    invariant_subsets = {}
    for name, group in (("syndrome_S5", image), ("syndrome_A5", perfect_core_image),
                        ("physical_D5", physical_permutations)):
        candidates = []
        for subset in range(32):
            points = {i for i in range(5) if (subset >> i) & 1}
            if all({perm[i] for i in points} == points for perm in group):
                candidates.append(len(points))
        invariant_subsets[name] = sorted(candidates)
    assert all(sizes == [0, 5] for sizes in invariant_subsets.values())

    syndrome_orbits = orbit_partition(derived2)
    full_spread_orbits = orbit_partition(automorphisms)
    deleted_matrices = {
        g: tuple(deleted_even_subset_action(image_of[g], v) for v in range(16))
        for g in derived2
    }
    even_subset_orbits = orbit_partition(set(deleted_matrices.values()))
    assert sorted(map(len, syndrome_orbits)) == [1, 15]
    assert sorted(map(len, full_spread_orbits)) == [1, 15]
    assert sorted(map(len, even_subset_orbits)) == [1, 5, 10]

    # Solve L G = D L over GF2 for all 16 coefficients of a 4x4 matrix L.
    # Each bit vector below is one homogeneous equation in those coefficients.
    equations = []
    for g in derived2:
        d = deleted_matrices[g]
        for column in range(4):
            for row in range(4):
                equation = 0
                for k in range(4):
                    if (g[1 << column] >> k) & 1:
                        equation ^= 1 << (row * 4 + k)
                    if (d[1 << k] >> row) & 1:
                        equation ^= 1 << (k * 4 + column)
                equations.append(equation)
    equation_rank = gf2_rank(equations)
    assert equation_rank == 16

    character_checks = {}
    for g in derived2:
        g_order = order(g)
        syndrome_fixed = sum(g[v] == v for v in range(16))
        fermion_trace = even_exterior_trace(image_of[g])
        pair = (syndrome_fixed, fermion_trace)
        if g_order in character_checks:
            assert character_checks[g_order] == pair
        character_checks[g_order] = pair
    assert character_checks == {1: (16, 16), 2: (4, 0), 3: (1, 4), 5: (1, 1)}

    result = {
        "syndrome_lines": [sorted(line) for line in lines],
        "GL4_2_maps_enumerated": census,
        "spread_automorphism_order": len(automorphisms),
        "line_permutation_image_order": len(image),
        "line_permutation_image": "S5, verified as equality with all 120 permutations",
        "linewise_kernel_order": len(kernel),
        "linewise_kernel_element_orders": sorted(order(g) for g in kernel),
        "exact_sequence": "1 -> C3 -> Aut(spread) -> S5 -> 1",
        "derived_series_orders": [len(group) for group in (automorphisms, derived1, derived2, derived3)],
        "canonical_perfect_core": "A5, verified by faithful equality of its image with all even permutations",
        "pure_physical_stabilizer_permutations_order": len(physical_permutations),
        "pure_physical_stabilizer_permutation_orders": dict(sorted(Counter(map(order, physical_permutations)).items())),
        "invariant_subset_sizes": invariant_subsets,
        "full_spread_orbits_on_syndromes": full_spread_orbits,
        "canonical_A5_orbits_on_syndromes": syndrome_orbits,
        "same_A5_orbits_on_even_subsets": even_subset_orbits,
        "GF2_intertwiner_equation_rank": equation_rank,
        "GF2_intertwiner_space_dimension": 16 - equation_rank,
        "character_comparison_by_element_order": {
            element_order: {"syndrome_fixed_points": pair[0], "even_exterior_trace": pair[1]}
            for element_order, pair in sorted(character_checks.items())
        },
        "monomial_lift_obstruction": "At order3 any unitary monomial lift on syndrome basis has |trace|<=1, but Lambda-even(C5) has trace4. No such intertwiner exists even with phases.",
        "scope": "Syndrome incidence symmetry; no physical continuous representation or E8 lattice has been constructed",
        "next_free_construction": "A universal projective lift of the canonical A5 would require an independently justified physical role",
    }
    output = json.dumps(result, indent=2)
    print(output)
    Path(__file__).with_name("syndrome_spread_output.json").write_text(output + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
