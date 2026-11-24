"""Exact event-to-lattice map for opposite flags in the Fano plane. Python standard library."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product

POINTS = range(1, 8)
LINES = tuple(sorted({tuple(sorted((p, q, p ^ q))) for p, q in combinations(POINTS, 2)}))
FLAGS = tuple((p, line) for line in LINES for p in line)


def next_flags(f, change_point):
    p, line = f
    return ([(q, line) for q in line if q != p] if change_point else
            [(p, other) for other in LINES if p in other and other != line])


def gallery_endpoints(f, sequence):
    histories = Counter({f: 1})
    for change_point in sequence:
        histories = Counter(g for h, n in histories.items() for _ in range(n)
                            for g in next_flags(h, change_point))
    return histories


def active_pairs(u):
    return [(x, x ^ u) for x in POINTS if x != u and x < x ^ u]


def twice_projector(u):
    matrix = [[0] * 8 for _ in range(8)]
    for x, y in active_pairs(u):
        matrix[x][x] += 1
        matrix[y][y] += 1
        matrix[x][y] -= 1
        matrix[y][x] -= 1
    return matrix


def fano_roots():
    # z represents the vector z/sqrt(2) in the point basis of Construction A.
    code = {tuple(((v & x).bit_count() + b) % 2 for x in range(8))
            for v in range(8) for b in range(2)}
    roots = set()
    for i, sign in product(range(8), (-1, 1)):
        z = [0] * 8
        z[i] = 2 * sign
        roots.add(tuple(z))
    for w in code:
        support = [i for i in range(8) if w[i]]
        if len(support) == 4:
            for signs in product((-1, 1), repeat=4):
                z = [0] * 8
                for i, sign in zip(support, signs):
                    z[i] = sign
                roots.add(tuple(z))
    return roots


def root_subsystems(u, roots):
    pairs = active_pairs(u)
    family = frozenset(z for z in roots if z[0] == z[u] == 0
                       and all(z[x] == -z[y] for x, y in pairs))
    carrier = frozenset(z for z in roots if all(z[x] == z[y] for x, y in pairs))
    return family, carrier


def linear_maps():
    for a, b, c in permutations(POINTS, 3):
        if a ^ b != c:
            yield (0, a, b, a ^ b, c, a ^ c, b ^ c, a ^ b ^ c)


def main():
    apartments = {frozenset(v) for v in combinations(POINTS, 3) if v[0] ^ v[1] != v[2]}
    assert len(FLAGS) == 21 and len(apartments) == 28
    counters = Counter()
    marks = Counter()
    marks_in_apartments = Counter()
    for ap in apartments:
        x, y, z = ap
        marks_in_apartments[x ^ y ^ z] += 1
    fixed_flag_counts = set()
    for f in FLAGS:
        p, line = f
        opposite = {g for g in FLAGS if p not in g[1] and g[0] not in line}
        assert gallery_endpoints(f, (True, False, True)) == dict.fromkeys(opposite, 1)
        assert gallery_endpoints(f, (False, True, False)) == dict.fromkeys(opposite, 1)
        local = Counter()
        for q, other in opposite:
            r, = set(line) & set(other)
            ap = frozenset((p, r, q))
            matches = [a for a in apartments if p in a and q in a
                       and len(a.intersection(line)) == len(a.intersection(other)) == 2]
            assert matches == [ap]
            u = p ^ r ^ q
            counters[ap] += 1
            marks[u] += 1
            local[u] += 1
        assert set(local) == set(POINTS)-set(line)
        fixed_flag_counts.add(tuple(sorted(local.values())))
    assert fixed_flag_counts == {(2, 2, 2, 2)}
    assert len(counters) == 28 and set(counters.values()) == {6}
    assert set(marks.values()) == {24} and set(marks_in_apartments.values()) == {4}
    # Independent uniform flags, and an independent reduced return gallery.
    one_event = Fraction(1, 21**2 * 8)
    longest_return = sum(marks.values()) * one_event
    assert longest_return == Fraction(1, 21)
    assert one_event / longest_return == Fraction(1, 168)
    assert set(n*one_event/longest_return for n in marks.values()) == {Fraction(1, 7)}

    roots = fano_roots()
    frames = {u: root_subsystems(u, roots) for u in POINTS}
    assert len(roots) == 240 and len(set(frames.values())) == 7
    assert len({fam for fam, car in frames.values()}) == 7
    assert {(len(fam), len(car)) for fam, car in frames.values()} == {(12, 40)}
    projectors = {u: twice_projector(u) for u in POINTS}
    for i, j in product(range(8), repeat=2):
        doubled_sum = sum(matrix[i][j] for matrix in projectors.values())
        expected = 0 if 0 in (i, j) else 7 * int(i == j) - 1
        assert doubled_sum == expected
    for perm in linear_maps():
        for u, m in projectors.items():
            dest = projectors[perm[u]]
            assert all(dest[perm[i]][perm[j]] == m[i][j] for i,j in product(range(8),repeat=2))

    print("flags, ordered opposite pairs, apartments:",len(FLAGS),sum(marks.values()),len(counters))
    print("opposite pairs per apartment / apartments per u / pairs per u:",6,4,24)
    print("longest successful return, each conditional pair, each conditional u:",
          longest_return,one_event/longest_return,24*one_event/longest_return)
    print("fixed initial flag: four u outside initial line, conditional probabilities 1/4 each")
    print("distinct embedded A3/D5 from this process:",len(set(frames.values())),"with 12/40 roots")
    print("averaged spectrum check / constant old points / six zero-sum old directions: 0, 0, 1/2")
    print("covariant under all 168 linear maps: OK")


if __name__ == "__main__":
    main()
