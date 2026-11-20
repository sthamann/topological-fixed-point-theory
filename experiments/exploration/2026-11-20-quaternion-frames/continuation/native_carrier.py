"""Exact native-code and polar-projector comparisons for the quaternion flag process."""

from collections import Counter
from itertools import combinations, product
from math import factorial, prod

from sympy import Matrix, Rational, S, diag, eye, ones, simplify, sqrt, zeros


LABELS = tuple(range(6))  # label 5 represents the point at infinity
UNIT = tuple(LABELS)
EDGES = tuple(combinations(LABELS, 2))


def compose(a, b):
    return tuple(a[i] for i in b)


def inverse(a):
    return tuple(a.index(i) for i in LABELS)


def projective_action(a, b, c, d):
    xy = [(i, 1) for i in range(5)] + [(1, 0)]
    def image(x, y):
        u, v = (a*x + b*y) % 5, (c*x + d*y) % 5
        return 5 if v == 0 else (u * pow(v, -1, 5)) % 5
    return tuple(image(x, y) for x, y in xy)


def frames():
    group = {projective_action(a, b, c, d) for a,b,c,d in product(range(5), repeat=4)
             if (a*d - b*c) % 5 == 1}
    assert len(group) == 60
    assert {compose(a, b) for a, b in product(group, repeat=2)} == group
    involutions = sorted(g for g in group if g != UNIT and compose(g, g) == UNIT)
    fours = {frozenset([UNIT, a, b, compose(a,b)]) for a,b in combinations(involutions, 2)
             if compose(a,b) == compose(b,a)}
    assert len(involutions) == 15 and len(fours) == 5
    edge_to_four = {}
    normalizers = {}
    for four in fours:
        pairs = {tuple(sorted({g[i] for g in four})) for i in LABELS}
        assert len(pairs) == 3 and all(len(p) == 2 for p in pairs)
        assert not set(edge_to_four).intersection(pairs)
        edge_to_four.update((pair, four) for pair in pairs)
        normalizers[four] = {g for g in group if
                             {compose(compose(g,h),inverse(g)) for h in four} == four}
        assert len(normalizers[four]) == 12
    assert set(edge_to_four) == set(EDGES)
    return group, edge_to_four, normalizers


def matching_marks(edges, transition):
    _, edge_to_four, _ = frames()
    matching_edges = {frame: {edge for edge, group in edge_to_four.items() if frame == group}
                      for frame in edge_to_four.values()}
    marks = Counter()
    for i, first in enumerate(edges):
        for j, second in enumerate(edges):
            if i == j or not transition[i, j]:
                continue
            origin, = set(first) & set(second)
            chord = tuple(sorted(set(first) ^ set(second)))
            moved = tuple(sorted(set(first) | set(second)))
            matching = matching_edges[edge_to_four[chord]]
            paired = next(edge for edge in matching if origin in edge)
            q, = set(paired) - {origin}
            assert q not in moved
            marks[moved,q] += 1
    expected = {(triple, q) for triple in combinations(LABELS,3) for q in LABELS if q not in triple}
    assert set(marks) == expected and set(marks.values()) == {2}
    print("Missing-chord selection: 120 ordered paths give all 60 moving-triple/anchor pairs, twice each")


def xor_charge(exponents):
    charge = 0
    for label, multiplicity in enumerate(exponents):
        if multiplicity % 2:
            charge ^= label
    return charge


def translates(exponents):
    return {tuple(exponents[label ^ shift] for label in range(4)) for shift in range(4)}


def pauli_invariants(degree):
    if degree % 4:
        return set()
    return {frozenset(translates(exponents))
            for exponents in product(range(degree + 1), repeat=4)
            if sum(exponents) == degree and xor_charge(exponents) == 0}


def quartic_code_control():
    assert [len(pauli_invariants(degree)) for degree in range(1, 5)] == [0, 0, 0, 5]
    expected = [frozenset(tuple(4 if label == index else 0 for label in range(4))
                          for index in range(4))]
    for displacement in (1, 2, 3):
        expected.append(frozenset(tuple(2 if label in (index, index ^ displacement) else 0
                                         for label in range(4)) for index in range(4)))
    expected.append(frozenset({(1, 1, 1, 1)}))
    assert set(expected) == pauli_invariants(4)
    norms = tuple(sum(prod(factorial(multiplicity) for multiplicity in monomial)
                      for monomial in support) for support in expected)
    assert norms == (96, 8, 8, 8, 1)
    assert Counter(len(monomials) for monomials in expected) == {4: 1, 2: 3, 1: 1}
    print(f"First nonconstant sharp-Heisenberg invariants: degree 4, native five quartics, Fock diagonal {norms}")


def first_exit_marks(edges, transition):
    """Two observed edges select a cycle; its first exit selects the fixed label."""
    proposals = set()
    for i, first in enumerate(edges):
        for j, second in enumerate(edges):
            if i == j or transition[i, j] == 0:
                continue
            moved = set(first) | set(second)
            assert len(moved) == 3
            outside = set(range(6)) - moved
            interior = [k for k,e in enumerate(edges) if set(e) <= moved]
            for incoming in interior:
                exits = dict.fromkeys(outside, S.Zero)
                for k, edge in enumerate(edges):
                    crossed = set(edge) & outside
                    if len(crossed) == 1:
                        exits[next(iter(crossed))] += transition[incoming, k]
                assert set(exits.values()) == {Rational(1,6)}
                assert sum(transition[incoming, k] for k in interior) == Rational(1,2)

            def swap(x, edge):
                a, b = edge
                return b if x == a else a if x == b else x

            cycle = tuple(swap(swap(x, second), first) for x in range(6))
            assert {x for x in range(6) if cycle[x] != x} == moved
            assert all(cycle[cycle[cycle[x]]] == x for x in range(6))
            proposals.update((cycle,q) for q in outside)
    assert len({cycle for cycle, _ in proposals}) == 40 and len(proposals) == 120
    print("Exit selection: 40 oriented 3-cycles, 120 pointed cycles; exit rate 1/2, each new label rate 1/6")


def polar_control():
    edges = tuple(combinations(range(6), 2))
    basis = eye(6)
    vectors = [basis[:, right] - basis[:, left] for left, right in edges]
    transition = Matrix(15, 15, lambda i, j: vectors[i].dot(vectors[j])**2 / S(12))
    edge_paths = Counter()
    for i, first in enumerate(edges):
        for j, second in enumerate(edges):
            if i == j or transition[i, j] == 0:
                continue
            assert transition[i, j] == Rational(1, 12)
            shared = set(first) & set(second)
            weak_pair = tuple(sorted(set(first) ^ set(second)))
            assert len(shared) == 1 and len(weak_pair) == 2
            edge_paths[next(iter(shared)), weak_pair] += 1

    assert len(edge_paths) == 60 and set(edge_paths.values()) == {2}
    print("Distinct allowed consecutive flags: 120 ordered paths to 60 (origin, weak pair) marks")
    matching_marks(edges, transition)
    first_exit_marks(edges, transition)

    # Native six-simplex coordinate chart serves only as an external comparison.
    raw_native = Matrix([[2, 2, -1, -1, -1, -1], [0, 0, -1, -1, 1, 1],
                         [0, 0, -1, 1, -1, 1], [0, 0, 1, -1, -1, 1], [1, -1, 0, 0, 0, 0]])
    native = diag(2, sqrt(12), sqrt(12), sqrt(12), sqrt(24)) * raw_native / sqrt(48)
    assert (native.T * native).applyfunc(simplify) == eye(6) - ones(6) / 6

    coordinate_deletion = basis[:, :5] - basis[:, 5] * ones(1, 5)
    gram = coordinate_deletion.T * coordinate_deletion
    assert gram == eye(5) + ones(5)
    normalizer = eye(5) + (1 / sqrt(6) - 1) * ones(5) / 5
    model_frame = (native * coordinate_deletion * normalizer).applyfunc(simplify)
    assert (model_frame.T * model_frame).applyfunc(simplify) == eye(5)
    model_weak = (model_frame[:, :2] * model_frame[:, :2].T).applyfunc(simplify)
    star_rays = native * coordinate_deletion
    star_effect = (star_rays * star_rays.T / 6).applyfunc(simplify)
    assert {simplify(value):count for value,count in star_effect.eigenvals().items()} == {1:1, Rational(1,6):4}
    star_root = (star_effect + eye(5)/sqrt(6))/(1+1/sqrt(6))
    # The positive square root is linear on the two distinct eigenvalues.
    assert (star_root**2 - star_effect).applyfunc(simplify) == zeros(5)
    assert (star_root * model_frame - star_rays/sqrt(6)).applyfunc(simplify) == zeros(5)

    # Original polar expression from check_polar_fibre.py.
    a = (native[:, 0] - native[:, 1]) / sqrt(2)
    e0 = sqrt(Rational(6, 5)) * native[:, 5]
    colour_sum = sum((native[:, j] for j in (2, 3, 4)), zeros(5, 1))
    e1 = sqrt(Rational(5, 6)) * (colour_sum + Rational(3, 5) * native[:, 5])
    c = -sqrt(Rational(3, 5)) * e0 + sqrt(Rational(2, 5)) * e1
    native_weak = (a * a.T + e0 * e0.T + e1 * e1.T - c * c.T).applyfunc(simplify)
    assert (model_weak - native_weak).applyfunc(simplify) == zeros(5)
    distance = (((model_weak - diag(1, 0, 0, 0, 1))**2).trace()).simplify()
    assert distance != 0
    print("Unique positive whitening of five incident edge rays: rank-2 projector agrees exactly with native P2")
    print("Star effect eigenvalues: 1, 1/6 x 4; positive filter followed by projective readout")
    print(f"Unwhitened coordinate-split counterexample: squared projector distance {distance}")


if __name__ == "__main__":
    quartic_code_control()
    polar_control()
