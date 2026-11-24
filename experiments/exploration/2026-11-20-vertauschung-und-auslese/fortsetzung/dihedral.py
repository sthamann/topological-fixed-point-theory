"""Dihedral chamber visibility, with exact class and chart controls. Python 3."""
from collections import Counter, deque
from fractions import Fraction


def mul(x, y, m):
    """Elements (k,j) represent r**k s**j; s r s = r**(-1)."""
    k, j = x
    h, i = y
    return ((k + (-1)**j*h) % m, j ^ i)


def conjugate(g, w, m):
    k, j = g
    inv = ((-(-1)**j*k) % m, j)
    return mul(mul(g, w, m), inv, m)


def spheres(m, qs, qt):
    if m < 2 or min(qs, qt) < 2 or (m % 2 and qs != qt):
        raise ValueError('Need a thick regular admissible parameter system')
    identity = (0, 0)
    # Choose s=(0,1), t=(-1,1) so that r=st=(1,0).
    neighbors = [((0, 1), qs), ((m-1, 1), qt)]
    distance, mass = {identity: 0}, {identity: 1}
    queue = deque([identity])
    while queue:
        w = queue.popleft()
        for g, q in neighbors:
            v = mul(w, g, m)
            level, weight = distance[w]+1, mass[w]*q
            if v not in distance:
                distance[v], mass[v] = level, weight
                queue.append(v)
            elif level == distance[v]:
                assert weight == mass[v]
    assert len(mass) == 2*m
    return distance, mass


def conjugacy_means(m, mass):
    """All charts preserve each class mass; the uniform chart achieves its mean."""
    remaining = set(mass) - {(0, 0)}
    result = {}
    while remaining:
        w = min(remaining)
        orbit = frozenset(conjugate(g,w,m) for g in mass)
        remaining -= orbit
        result[orbit] = Fraction(sum(mass[e] for e in orbit), len(orbit))
    return result


def finite_checks():
    for m, qs, qt in ((2,2,2),(3,2,2),(3,3,3),(4,2,4),(6,2,2)):
        distance, mass = spheres(m, qs, qt)
        means = conjugacy_means(m, mass)
        print(f'I2({m}), parameters ({qs},{qt}): total {sum(mass.values())}, '
              f'distinct active class means {sorted(set(means.values()))}, '
              f'equalizable {len(set(means.values()))==1}')
        if m >= 4:
            for k in (1, 2):
                orbit = {conjugate(g,(k,0),m) for g in mass}
                assert {mass[w] for w in orbit} == {(qs*qt)**k}
                assert {distance[w] for w in orbit} == {2*k}

    diagnostics = []
    for m in range(2, 13):
        for qs in range(2, 7):
            for qt in range(2, 7):
                if m%2 and qs!=qt:
                    continue
                _, mass = spheres(m, qs, qt)
                means = conjugacy_means(m, mass)
                diagnostics.append((m, qs, qt, len(set(means.values()))==1))
    hits = [(m,qs,qt) for m,qs,qt,passed in diagnostics if passed]
    assert hits == [(3,2,2)]
    print('parameter-system control:',len(hits),'/',len(diagnostics),hits,
          '(includes formal types without asserted existence)')

    # At (3,2,2) there are six inner automorphisms but just three choices of
    # longest reflection. Pick deliberately unequal weights inside each pair.
    m = 3
    distance, mass = spheres(m,2,2)
    assert sum(mass.values())==21
    longest, = (w for w in mass if distance[w]==3)
    groups = {}
    for g in sorted(mass):
        groups.setdefault(conjugate(g,longest,m),[]).append(g)
    assert sorted(map(len,groups.values())) == [2,2,2]
    choices = (Fraction(0),Fraction(1,12),Fraction(1,6))
    chart = {}
    for pair,first in zip(groups.values(),choices):
        chart[pair[0]],chart[pair[1]] = first,Fraction(1,3)-first
    averaged = Counter()
    for g, prob in chart.items():
        for w, count in mass.items():
            averaged[conjugate(g,w,m)] += prob*count
    assert {v for w,v in averaged.items() if w!=(0,0)} == {4}
    print('class weights 1/3 each suffice even with nonuniform six chart weights')
    print('active endpoint mass each:',Fraction(4,21),
          '; conditional on being active:',Fraction(1,5))


if __name__ == '__main__':
    finite_checks()
