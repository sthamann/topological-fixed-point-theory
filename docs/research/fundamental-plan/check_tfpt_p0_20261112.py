from fractions import Fraction as Q
from itertools import combinations, permutations, product


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def det(matrix):
    rows = [[Q(entry) for entry in row] for row in matrix]
    result = Q(1)
    for i in range(len(rows)):
        pivot = next((j for j in range(i, len(rows)) if rows[j][i]), None)
        if pivot is None:
            return Q(0)
        if pivot != i:
            rows[i], rows[pivot] = rows[pivot], rows[i]
            result *= -1
        result *= rows[i][i]
        for j in range(i + 1, len(rows)):
            factor = rows[j][i] / rows[i][i]
            for k in range(i, len(rows)):
                rows[j][k] -= factor * rows[i][k]
    return result


def lattice():
    basis = [tuple(int(i == j) for i in range(8)) for j in range(8)]
    d8 = {
        tuple(a * x + b * y for x, y in zip(basis[i], basis[j]))
        for i, j in combinations(range(8), 2)
        for a, b in product((-1, 1), repeat=2)
    }
    spinors = {
        tuple(Q(s, 2) for s in signs)
        for signs in product((-1, 1), repeat=8)
        if signs.count(-1) % 2 == 0
    }
    roots = d8 | spinors
    assert len(roots) == 240
    p1 = (Q(1, 2),) * 8
    p2 = tuple(Q(-1 if i in (0, 5) else 1, 2) for i in range(8))
    gamma = add(basis[0], basis[5])
    assert p1 == add(p2, gamma) and {p1, p2, gamma} <= roots
    h = tuple(-2 * x for x in basis[7])
    a = tuple(x - y for x, y in zip(basis[1], basis[7]))
    b = tuple(-x - y for x, y in zip(basis[1], basis[7]))
    assert h == add(a, b) and a in roots and b in roots
    assert dot(a, b) == dot(a, gamma) == dot(b, gamma) == 0
    assert [dot(x, p) for x in (h, a, b) for p in (p1, p2)] == [-1, -1, 0, 0, -1, -1]
    minus_a = tuple(-x for x in a)
    assert all(add(p, x) not in roots for p in (p1, p2) for x in (a, minus_a))
    assert add(p1, b) in roots and add(p2, b) in roots
    assert add(p1, b) == add(add(p2, b), gamma)
    assert [dot(add(p, b), b) for p in (p1, p2)] == [1, 1]
    e7 = {r for r in roots if dot(r, gamma) == 0}
    e6 = {r for r in e7 if dot(r, p1) == 0}
    assert len(e7) == 126 and len(e6) == 72
    h_cartan = add(p1, p2)
    assert {r for r in roots if dot(r, h_cartan) == 3} == {p1, p2}
    positive = {r for r in e6 if next(x for x in r if x) > 0}
    simple = sorted(positive - {add(r, s) for r in positive for s in positive})
    assert len(simple) == 6
    assert det([[dot(r, s) for s in simple] for r in simple]) == 3
    marked_basis = simple + [h]
    assert det([[dot(r, s) for s in marked_basis] for r in marked_basis]) == 2
    # This determinant equals det(gamma^perp in the unimodular E8 lattice).
    for k in (-3, -1, 0, 1, 2, 3):
        assert dot(tuple(k * x for x in h), p1) == -k
        assert (-1) ** (k % 2) == (-1) ** (int(dot(h, p1) * k) % 2)
    print("Windung: 126/72 Wurzeln, Gramdeterminanten 2/3, Teilwindungen 0 und -1.")
    print("Wurzelkontraktion: a fixiert die Enden; b tauscht beide mit p_j+b.")


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def inverse(p):
    return tuple(p.index(i) for i in range(len(p)))


def cycle_power(p, n):
    result = tuple(range(len(p)))
    for _ in range(n):
        result = compose(result, p)
    return result


def clocks():
    identity = tuple(range(6))
    g = (1, 2, 3, 0, 4, 5)
    cycles = [
        p for p in permutations(identity)
        if p != identity and cycle_power(p, 5) == identity
        and compose(compose(g, p), inverse(g)) == cycle_power(p, 2)
    ]
    a = next(p for p in cycles if p[4] == 4)
    b = next(p for p in cycles if p[5] == 5)
    assert len(cycles) == 8
    generated = {identity}
    queue = [identity]
    for current in queue:
        for step in (a, b):
            image = compose(current, step)
            if image not in generated:
                generated.add(image)
                queue.append(image)
    assert len(generated) == 360
    stabilizer = [p for p in permutations(identity) if p[4] == 4]
    powers = {cycle_power(a, n) for n in range(1, 5)}
    normalizer = [p for p in stabilizer if compose(compose(p, a), inverse(p)) in powers]
    involutions = [p for p in normalizer if p != identity and cycle_power(p, 2) == identity]
    assert len(normalizer) == 20 and len(involutions) == 5
    # Number of negative eigenvalues on the real augmentation space.
    assert {sum(p[i] != i for i in range(6)) // 2 for p in involutions} == {2}
    print("Geometrie: A6 mit 360 Elementen; fünf Normalisator-Involutionen, alle Signatur (2,2).")


def point_kernel(points):
    return {
        (i, j): 1 / sum((Q(x) - Q(y)) ** 2 for x, y in zip(points[i], points[j]))
        for i, j in combinations(range(4), 2)
    }


def projective_test(c):
    a2 = 1 / (c[0, 1] * c[2, 3])
    b2 = 1 / (c[0, 2] * c[1, 3])
    c2 = 1 / (c[0, 3] * c[1, 2])
    branch = b2 - a2 - c2
    return branch, branch ** 2 - 4 * a2 * c2


def geometry():
    quarter = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    moved = [quarter[0], (Q(5, 13), Q(12, 13)), *quarter[2:]]
    assert dot(moved[1], moved[1]) == 1
    first, second = point_kernel(quarter), point_kernel(moved)
    assert projective_test(first) == (8, 0)
    assert projective_test(second) == (Q(96, 13), 0)
    local_normalization = [2, 3, 5, 7]
    scaled = {(i, j): c * local_normalization[i] * local_normalization[j]
              for (i, j), c in second.items()}
    assert projective_test(scaled)[1] == 0
    mixed = {ij: (first[ij] + second[ij]) / 2 for ij in first}
    branch, residual = projective_test(mixed)
    assert branch > 0 and residual != 0
    print(f"Kreis: beide Einzelresiduen 0, Mittelresiduum {residual} ({float(residual):.12f}).")


if __name__ == "__main__":
    lattice()
    clocks()
    geometry()
