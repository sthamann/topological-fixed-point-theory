"""Exact root and lattice inputs for the two compact marked stabilizers."""

from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, product
from pathlib import Path
import runpy


HERE = Path(__file__).resolve().parent
shared = runpy.run_path(str(HERE / "check_tfpt_p0_20261112.py"))
dot, add, det = (shared[name] for name in ("dot", "add", "det"))


def unit(i):
    return tuple(int(j == i - 1) for j in range(8))


def neg(v):
    return tuple(-x for x in v)


def sub(v, w):
    return add(v, neg(w))


def rank(rows):
    basis = []
    for row in rows:
        vector = list(map(Q, row))
        for pivot, previous in basis:
            factor = vector[pivot]
            vector = [x - factor * y for x, y in zip(vector, previous)]
        if any(vector):
            pivot = next(i for i, x in enumerate(vector) if x)
            basis.append((pivot, [x / vector[pivot] for x in vector]))
    return len(basis)


def gram_det(vectors):
    return det([[dot(v, w) for w in vectors] for v in vectors])


def simple_roots(roots):
    positive = {r for r in roots if next(x for x in r if x) > 0}
    return sorted(positive - {add(r, s) for r in positive for s in positive})


roots = {
    add(tuple(s*x for x in unit(i)), tuple(t*x for x in unit(j)))
    for i, j in combinations(range(1, 9), 2)
    for s, t in product((-1, 1), repeat=2)
}
roots |= {
    tuple(Q(x, 2) for x in signs)
    for signs in product((-1, 1), repeat=8)
    if signs.count(-1) % 2 == 0
}
p1 = (Q(1, 2),) * 8
gamma = add(unit(1), unit(6))
p2 = sub(p1, gamma)
a = sub(unit(7), unit(8))
b = neg(add(unit(7), unit(8)))
h = add(a, b)
grading = add(p1, p2)
assert len(roots) == 240
assert {r for r in roots if dot(r, grading) == 3} == {p1, p2}
e7 = {r for r in roots if dot(r, gamma) == 0}
e6 = {r for r in e7 if dot(r, p1) == 0}
assert Counter(dot(r, p1) for r in e7) == {-1: 27, 0: 72, 1: 27}
assert len(simple_roots(e7)) == 7 and gram_det(simple_roots(e7)) == 2
assert len(simple_roots(e6)) == 6 and gram_det(simple_roots(e6)) == 3
assert gram_det(simple_roots(e6) + [h]) == 2
for p in (p1, p2):
    assert dot(p, a) == 0 and dot(p, b) == -1
    assert add(p, a) not in roots and sub(p, a) not in roots
    assert add(p, b) in roots and sub(p, b) not in roots

# The native gauge frame and eighth-clock character from the original check.
colour_1 = sub(unit(1), unit(2))
colour_2 = sub(unit(2), unit(3))
weak = sub(unit(4), unit(5))
hypercharge = (-Q(1, 3),) * 3 + (Q(1, 2),) * 2 + (0,) * 3
fixed = (colour_1, colour_2, weak, hypercharge, gamma)
clock = (0, 2, 4, -2, 0, 4, 2, 0)
assert rank(fixed) == 5
allowed = {r for r in roots if all(dot(r, v) == 0 for v in fixed)}
assert len(allowed) == 8
assert {r for r in allowed if dot(r, clock) % 2 == 0} == {a, b, neg(a), neg(b)}
assert not {r for r in allowed if dot(r, clock) % 8 == 0}
t = tuple(Q(s, 2) for s in (1, 1, 1, 1, 1, -1, 1, -1))
assert t in roots and all(dot(t, v) == 0 for v in fixed)
assert [dot(p1, r) for r in (t, a, b, h)] == [1, 0, -1, -1]
assert [dot(clock, p) % 8 for p in (p1, p2, add(p1, b), add(p2, b))] == [5, 1, 3, 7]

# Exhaust the integral and half-integral possibilities modulo the period
# of the even-sum lattice in the three free coordinates (x,u,v).
for x, u, v in product(tuple(Q(k, 2) for k in range(-4, 5)), repeat=3):
    if len({x % 1, u % 1, v % 1}) != 1 or (4*x + u + v) % 2:
        continue
    coords = (2*x, (u - v - 2*x)/2, -(u + v)/2)
    assert all(c.denominator == 1 for c in coords)
    original = (x, x, x, x, x, -x, u, v)
    assert original == tuple(sum(c*r for c, r in zip(coords, column))
                             for column in zip(t, a, b))

print("E7 roots by endpoint degree: 27 / 72 / 27. Gram determinants E7 = 2, E6 = 3, E6 plus h = 2.")
print("Fixed SM frame and middle current: 5 independent Cartan constraints, torus rank 3.")
print("Root counts with SM and middle current, then half parity, then full clock: 8 / 4 / 0.")
print("Native cocharacter basis t,a,b. Mark degrees 1,0,-1. Endpoint clock residues (5,1) -> (3,7).")
