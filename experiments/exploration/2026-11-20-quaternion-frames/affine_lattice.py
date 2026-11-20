#!/usr/bin/env python3
"""Exact conditional E8, torsion and spinor-edge test, requiring SymPy.

The lattice origin used here is a theorem: the smallest nonzero positive
definite even unimodular lattice has rank eight and is E8. This script
starts its enumeration with a standard realization of that resulting root
system. It tests the further least-displacement proposal, including Kac
labels excluded by the semisimple-centralizer condition at order four.
"""

from collections import Counter
from itertools import combinations, product
from math import comb, gcd

from sympy import Matrix, Rational


def twice_e8_roots():
    for i, j in combinations(range(8), 2):
        for a, b in product((-2, 2), repeat=2):
            root = [0] * 8
            root[i], root[j] = a, b
            yield tuple(root)
    yield from (root for root in product((-1, 1), repeat=8) if sum(root) % 4 == 0)


SIMPLE = Matrix([
    [1, -1, -1, -1, -1, -1, -1, 1],
    [2, 2, 0, 0, 0, 0, 0, 0],
    [-2, 2, 0, 0, 0, 0, 0, 0],
    [0, -2, 2, 0, 0, 0, 0, 0],
    [0, 0, -2, 2, 0, 0, 0, 0],
    [0, 0, 0, -2, 2, 0, 0, 0],
    [0, 0, 0, 0, -2, 2, 0, 0],
    [0, 0, 0, 0, 0, -2, 2, 0],
])
COWEIGHTS = 2 * SIMPLE.inv()
VERTEX_TYPES = (
    "D8", "A8", "A7+A1", "A5+A2+A1", "A4+A4", "D5+A3", "E6+A2", "E7+A1"
)


def grades(coefficients, labels, order):
    return [
        sum(coefficient * label for coefficient, label in zip(root, labels[1:])) % order
        for root in coefficients
    ]


def normalized_variance(buckets):
    total = sum(buckets)
    order = len(buckets)
    return Rational(sum((order * count - total) ** 2 for count in buckets),
                    order * total**2)


def stats(coefficients, labels, order):
    buckets = Counter(grades(coefficients, labels, order))
    histogram = tuple(buckets[j] for j in range(order))
    position = COWEIGHTS * Matrix(labels[1:]) / order
    return histogram, position.dot(position), normalized_variance(histogram)


def kac_labels(marks, budget):
    if not marks:
        if budget == 0:
            yield ()
        return
    first, *rest = marks
    for label in range(budget // first + 1):
        for tail in kac_labels(rest, budget - label * first):
            yield (label, *tail)


def affine_audit(coefficients):
    highest = max(coefficients, key=sum)
    assert highest == (2, 3, 4, 6, 5, 4, 3, 2)
    print(f"Roots: {len(coefficients)}. Derived affine marks: {highest}.")
    print("Full-rank vertex classes: node, type, order, root counts, squared norm, variance")

    eligible = []
    for node, (order, kind) in enumerate(zip(highest, VERTEX_TYPES), start=1):
        labels = tuple(int(j == node) for j in range(9))
        counts, norm, variance = stats(coefficients, labels, int(order))
        print(f"  {node}: {kind:9s}  {order}  {counts}  {norm}  {variance}")
        has_mixed_reality = order >= 4 and order % 2 == 0
        if has_mixed_reality:
            eligible.append((norm, variance, node, labels))

    by_distance = min(eligible, key=lambda entry: entry[0])
    by_variance = min(eligible, key=lambda entry: entry[1])
    assert by_distance[2] == by_variance[2] == 6
    print("Winner under mixed reality: D5+A3. Controls: 1 of 3 eligible classes, 1 of 8 overall.")

    print("All exact order-four Kac classes: nonzero labels, semisimple rank, roots, norm")
    order_four = []
    for labels in kac_labels((1, *map(int, highest)), 4):
        if gcd(*labels) != 1:
            continue
        counts, norm, _ = stats(coefficients, labels, 4)
        support = tuple((j, value) for j, value in enumerate(labels) if value)
        semisimple_rank = 9 - len(support)
        order_four.append((norm, semisimple_rank, support))
        print(f"  {support}: {semisimple_rank}, {counts}, {norm}")

    shortest_relaxed = min(order_four)
    assert shortest_relaxed[1] < 8 and shortest_relaxed[0] < by_distance[0]
    print(f"Without the full-rank requirement: shortest order-four squared norm "
          f"{shortest_relaxed[0]}, semisimple rank {shortest_relaxed[1]}.")
    return by_distance[3]


def projections(roots, root_grades, highest):
    d5_simple_indices = (0, 2, 3, 1, 4)
    d5_in_e8 = SIMPLE[list(d5_simple_indices), :]
    d5_model = Matrix([
        [1, -1, 0, 0, 0],
        [0, 1, -1, 0, 0],
        [0, 0, 1, -1, 0],
        [0, 0, 0, 1, -1],
        [0, 0, 0, 1, 1],
    ])
    assert d5_model * d5_model.T == d5_in_e8 * d5_in_e8.T / 4
    d5_map = d5_model.inv() * d5_in_e8 / 4
    highest_row = Matrix(1, 8, highest) * SIMPLE
    a3_in_e8 = Matrix.vstack(-highest_row, SIMPLE[7, :], SIMPLE[6, :])
    assert a3_in_e8 * d5_in_e8.T == Matrix.zeros(3, 5)
    assert a3_in_e8 * a3_in_e8.T / 4 == Matrix([[2, -1, 0], [-1, 2, -1], [0, -1, 2]])
    a3_map = a3_in_e8 / 4
    projected = [set() for _ in range(4)]
    for root, grade in zip(roots, root_grades):
        weight = Matrix(root)
        projected[grade].add((tuple(d5_map * weight), tuple(a3_map * weight)))
    return projected


def whole_grade_edges(pairs):
    a3_cartan_inverse = Matrix([[2, -1, 0], [-1, 2, -1], [0, -1, 2]]).inv()
    assert len({a3 for _, a3 in pairs}) == 4
    neighbors = {pair: set() for pair in pairs}
    for first, second in combinations(sorted(pairs), 2):
        first_d5, first_a3 = first
        second_d5, second_a3 = second
        d5_delta = [a - b for a, b in zip(first_d5, second_d5)]
        a3_delta = Matrix(first_a3) - Matrix(second_a3)
        distance = sum(value**2 for value in d5_delta) + a3_delta.dot(a3_cartan_inverse * a3_delta)
        assert distance >= 2
        if distance == 2:
            neighbors[first].add(second)
            neighbors[second].add(first)

    edge_types = Counter()
    for first, adjacent in neighbors.items():
        assert len(adjacent) == 13
        for second in adjacent:
            common_neighbors = len(adjacent & neighbors[second])
            changed_d5 = first[0] != second[0]
            changed_a3 = first[1] != second[1]
            edge_types[changed_d5, changed_a3, common_neighbors] += 1
    assert edge_types == {(True, False, 6): 640, (False, True, 2): 192}
    print("Whole-grade shortest edges by changed factor: D5 640 (6 common neighbors), "
          "A3 192 (2 common neighbors).")
    print("Maximal common-neighbor multiplicity marks one unchanged A3 weight; "
          "3 others remain. Unconditioned shortest-edge baseline: 10/13.")


def spinor_audit(roots, coefficients, labels):
    projected = projections(roots, grades(coefficients, labels, 4), max(coefficients, key=sum))
    print("Non-neutral grade: D5 weights × A3 weights = product count")
    for grade, expected in ((1, (16, 4)), (2, (10, 6)), (3, (16, 4))):
        pairs = projected[grade]
        d5, a3 = map(set, zip(*pairs))
        assert (len(d5), len(a3)) == expected and len(pairs) == len(d5) * len(a3)
        print(f"  {grade}: {len(d5)} × {len(a3)} = {len(pairs)}")

    half = Rational(1, 2)
    reference = (half,) * 5
    spinor_grade = next(grade for grade in (1, 3)
                        if any(d5 == reference for d5, _ in projected[grade]))
    spinors = {d5 for d5, _ in projected[spinor_grade]}
    assert spinors == {weight for weight in product((-half, half), repeat=5)
                       if weight.count(-half) % 2 == 0}
    whole_grade_edges(projected[spinor_grade])
    d5_roots = {d5 for d5, _ in projected[0] if any(d5)}

    distances = Counter()
    for first, second in product(sorted(spinors), repeat=2):
        if first == second:
            continue
        difference = tuple(a - b for a, b in zip(first, second))
        squared_distance = sum(value**2 for value in difference)
        changed = sum(a != b for a, b in zip(first, second))
        distances[squared_distance, changed] += 1
        assert (difference in d5_roots) == (squared_distance == 2)
        # Even signed flips map first to the all-positive reference. A permutation
        # maps the second member to a chosen representative of its changed count.
        assert tuple(sorted(a * b / half for a, b in zip(first, second))) == (
            (-half,) * changed + (half,) * (5 - changed)
        )
    assert distances == {(2, 2): 160, (4, 4): 80}
    print(f"Spinor grade {spinor_grade}. Ordered edges by (squared distance, changed slots): "
          f"{dict(sorted(distances.items()))}.")

    s = Matrix(reference)
    nearest = Matrix((half, half, half, -half, -half))
    q = nearest - s.dot(nearest) / s.dot(s) * s
    q = -q  # Assign positive charge to the smaller block.
    q /= sum(value for value in q if value > 0)
    assert tuple(q) == (Rational(-1, 3),) * 3 + (half,) * 2
    neutral = [root for root in d5_roots if Matrix(root).dot(s) == Matrix(root).dot(nearest) == 0]
    assert len(neutral) == 8 and Matrix(neutral).rank() == 3
    assert sorted(sum(bool(root[i]) for root in neutral) for i in range(5)) == [2, 2, 4, 4, 4]
    print(f"Nearest edge: 3+2, Y = {tuple(q)}, neutral D5 roots = {len(neutral)} (A2+A1).")
    print("Control: unordered block-size types 1 of 2; random oriented weight pairs 160 of 240.")

    exterior = Counter()
    for a, b in product(range(4), range(3)):
        if (a + b) % 2 == 0:
            exterior[a * q[0] + b * q[3]] += comb(3, a) * comb(2, b)
    assert sum(exterior.values()) == 16
    print(f"Even exterior powers: charge multiplicities {dict(sorted(exterior.items()))}.")


def main():
    roots = tuple(twice_e8_roots())
    coefficients = [tuple(Matrix(1, 8, root) * COWEIGHTS / 2) for root in roots]
    assert len(set(roots)) == len(coefficients) == 240
    assert all(value.is_Integer for root in coefficients for value in root)
    labels = affine_audit(coefficients)
    spinor_audit(roots, coefficients, labels)


if __name__ == "__main__":
    main()
