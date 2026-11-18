"""Exact cross-check of the two Nov 10 TFPT clock proposals. Python stdlib."""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product


def root(i, j, s=1, t=1):
    return tuple(2*s if k == i else 2*t if k == j else 0 for k in range(8))


def sigma(v):
    # Row convention of verification/v634_st31_structure.py.
    return tuple(v[k] for k in (4, 5, 0, 1, 2, 3, 6, 7))


def quarter(v):
    return tuple(w for k in range(0, 8, 2) for w in (-v[k+1], v[k]))


def clock(v):
    return quarter(sigma(v))


def power(f, v, n):
    for _ in range(n):
        v = f(v)
    return v


def glue(v):
    return sum(v[:5]) % 4


def dot(v, w):
    return F(sum(a*b for a, b in zip(v, w)), 4)


roots = {root(i, j, s, t) for i, j in combinations(range(8), 2)
         for s, t in product((-1, 1), repeat=2)}
roots |= {v for v in product((-1, 1), repeat=8) if sum(v) % 4 == 0}
assert len(roots) == 240
assert {clock(v) for v in roots} == roots
assert all(power(clock, v, 4) == sigma(v) and power(clock, v, 9) == quarter(v)
           for v in roots)
counts = [sum(glue(power(clock, v, n)) != glue(v) for v in roots)
          for n in range(1, 13)]
assert counts == [112]*5 + [128] + [112]*5 + [0]

# The Fold's q,r in the same chart, as given in TFPTDESCENDANTSOURCE20261110.
q, r = root(7, 6, 1, -1), root(4, 5)
q_orbit = {power(clock, q, n) for n in range(12)}
r_orbit = {power(clock, r, n) for n in range(12)}
unoriented = {min(v, tuple(-x for x in v)) for v in q_orbit | r_orbit}
assert (len(q_orbit), len(r_orbit), len(unoriented)) == (4, 12, 8)
assert all(dot(v, v) == 2 for v in unoriented)
assert all(dot(v, w) == 0 for v, w in combinations(unoriented, 2))
assert all(dot(v, w) == 0 for v in (q, r) for w in (clock(q), clock(r)))
perp = {v for v in roots if dot(v, q) == dot(v, r) == 0}
clock_d6 = {v for v in roots if sum(x*x for x in v[:6]) == 8}
assert len(perp) == len(clock_d6) == 60 and len(perp & clock_d6) == 26
orbit_sizes = Counter(len({power(clock, v, n) for n in range(12)}) for v in roots)

# Exact folded classical ring and its mode that the Clifford channel lacks.
ring = [[{0: F(1, 2), 1: F(2, 9), 3: F(1, 18), 5: F(2, 9)}.get((j-i)%6, F(0))
         for j in range(6)] for i in range(6)]
parity = [(-1)**i for i in range(6)]
assert all(sum(row[j]*parity[j] for j in range(6)) == 0 for row in ring)
pairs = ((1, 2), (4, 5), (0, 3))
folded = [[sum(ring[i][j] for i in a for j in b)/2 for b in pairs] for a in pairs]
assert folded == [[F(13,18), F(1,18), F(2,9)], [F(1,18), F(13,18), F(2,9)],
                  [F(2,9), F(2,9), F(5,9)]]
assert all(sum(parity[i] for i in a) == 0 for a in pairs)
clifford_eigs = {d: F(7+(-1)**d*(5-2*d), 12) for d in (0, 1, 2)}
assert set(clifford_eigs.values()) == {F(1), F(1,3), F(2,3)}
assert all(clifford_eigs.values())

print('Native clock glue changes:', counts)
print('Native clock orbit census:', {k: v//k for k, v in sorted(orbit_sizes.items())})
print('Fold q and r orbit sizes:', len(q_orbit), len(r_orbit))
print('Clock orbit lattice: A1^8; index in E8 = 16; D6 common roots =', len(perp & clock_d6))
print('Classical ring quotient:', [[str(x) for x in row] for row in folded])
print('Ring kills parity:', parity, '; each single-site projector exposes it')
print('Clifford channel spectrum:', {k: str(v) for k, v in clifford_eigs.items()})
