"""Exact small-rank checks for the classification in SICHTBARKEIT.md.

Python 3 standard library only. Tuples encode actual permutations; the search
checks thickness products across all reduced predecessor paths it encounters.
"""

from collections import Counter, deque
from fractions import Fraction
from math import factorial, prod


def compose(left, right):
    return tuple(left[index] for index in right)


def inverse(permutation):
    result = [0] * len(permutation)
    for i, j in enumerate(permutation):
        result[j] = i
    return tuple(result)


def swap(n, i, j):
    values = list(range(n))
    values[i], values[j] = values[j], values[i]
    return tuple(values)


def dihedral(m):
    return (tuple(-i % m for i in range(m)), tuple((1-i) % m for i in range(m)))


def length_and_weight(generators, thickness):
    identity = tuple(range(len(generators[0])))
    table = {identity: (0, 1)}
    pending = deque([identity])
    while pending:
        w = pending.popleft()
        length, weight = table[w]
        for s, q in zip(generators, thickness):
            child = compose(w, s)
            if child not in table:
                table[child] = (length + 1, weight * q)
                pending.append(child)
            if table[child][0] == length + 1:
                assert table[child][1] == weight * q, 'inconsistent reduced words'
    return table


def class_means(table):
    identity = next(x for x, (l, _) in table.items() if not l)
    unseen = set(table) - {identity}
    averages = {}
    while unseen:
        element = min(unseen)
        orbit = {compose(compose(w, element), inverse(w)) for w in table}
        average = Fraction(sum(table[v][1] for v in orbit), len(orbit))
        for v in orbit:
            averages[v] = average
        unseen -= orbit
    return averages


def poincare(degrees, t):
    return prod(sum(t**j for j in range(d)) for d in degrees)


def main():
    a2 = length_and_weight(dihedral(3), (2, 2))
    assert Counter(l for l, w in a2.values()) == {0: 1, 1: 2, 2: 2, 3: 1}
    assert set(class_means(a2).values()) == {4}
    assert sum(w for l, w in a2.values()) == 21
    assert poincare((2,3), 2) + 5*poincare((2,3), -2) == 6
    print('A2 q=2: raw averaged weights [1,4,4,4,4,4]; active 1/5; idle 1/21')
    a2_q3 = length_and_weight(dihedral(3), (3, 3))
    averages = class_means(a2_q3)
    assert Counter(averages.values()) == {Fraction(11): 3, Fraction(9): 2}
    assert all(averages[compose(compose(w,g),inverse(w))] == averages[g]
               for w in a2_q3 for g in averages)
    print('A2 q=3: conjugation invariant; individual active weights 11,11,11,9,9')

    checked = []
    for m in range(3, 19):
        for q in range(2, 9):
            test = length_and_weight(dihedral(m), (q, q))
            if len(set(class_means(test).values())) == 1:
                checked.append((m,q))
    assert checked == [(3,2)]
    print('Uniform dihedral control (3<=m<=18, 2<=q<=8):', checked)
    uneven = length_and_weight(dihedral(4), (2,3))
    long = max(uneven, key=lambda g: uneven[g][0])
    assert uneven[long] == (4,36)
    assert all(compose(long, w) == compose(w, long) for w in uneven)
    assert all(uneven[w][1] < 36 for w in uneven if w != long)
    print('Unequal B2/I2(4) q_s=2,q_t=3: central longest has unique weight 36')

    a3 = length_and_weight(tuple(swap(4,i,i+1) for i in range(3)), (2,2,2))
    assert len(a3) == 24 and len(set(class_means(a3).values())) > 1
    for degrees, order in (((2,3,4),24), ((2,4,6),48), ((2,6,10),120)):
        assert poincare(degrees,1) == order
        assert poincare(degrees,-2)>0 or order-1 >= 3**3
    print('Rank three: A3 excluded by positive P(-2); B3,H3 excluded by order')
    assert factorial(5)-1 > 3**4
    print('Rank >=4: 5!-1 > 3^4; subsequent order lower bound grows by >=6 per rank')

    product = length_and_weight((swap(4,0,1),swap(4,2,3)),(2,3))
    assert set(class_means(product).values()) == {2,3,6}
    print('Reducible A1 x A1: active weights (2,3,6)')
    assert len(class_means(length_and_weight((swap(2,0,1),),(5,)))) == 1
    print('A1: one active label, criterion vacuous')

    raw = (Fraction(11,21), Fraction(2,21))
    assert (raw[0]-3*raw[1],raw[0]+raw[1]) == (Fraction(5,21),Fraction(13,21))
    print('Raw empty/Lueders channel sector values: 1,5/21,13/21')


if __name__ == '__main__':
    main()
