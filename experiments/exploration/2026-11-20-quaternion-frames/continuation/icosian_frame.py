#!/usr/bin/env python3
"""Construct a native Gaussian E8 quarter-frame from one oriented icosian axis."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product


Z = (F(0), F(0))
ONE = (F(1), F(0))
H = (F(1, 2), F(0))
PHI = (F(1, 2), F(1, 2))
PHI_CONJ = (F(1, 2), -F(1, 2))


def add(*terms):
    return tuple(map(sum, zip(*terms)))


def mul(a, b):
    return a[0] * b[0] + 5 * a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def neg(a):
    return -a[0], -a[1]


def scale(scalar, vector):
    if not isinstance(scalar, tuple):
        scalar = (F(scalar), F(0))
    return tuple(mul(scalar, entry) for entry in vector)


def qadd(*vectors):
    return tuple(add(*entries) for entries in zip(*vectors))


def qmul(a, b):
    w, x, y, z = a
    u, p, q, r = b
    return (add(mul(w, u), neg(mul(x, p)), neg(mul(y, q)), neg(mul(z, r))),
            add(mul(w, p), mul(x, u), mul(y, r), neg(mul(z, q))),
            add(mul(w, q), neg(mul(x, r)), mul(y, u), mul(z, p)),
            add(mul(w, r), mul(x, q), neg(mul(y, p)), mul(z, u)))


def dot(a, b):
    p, q = add(*(mul(x, y) for x, y in zip(a, b)))
    return 2 * (p + q)


def group_units():
    identity = (ONE, Z, Z, Z)
    generators = ((H, H, H, H), (mul(H, PHI), mul(H, PHI_CONJ), H, Z))
    units = {identity}
    queue = [identity]
    for unit in queue:
        for gen in generators:
            candidate = qmul(unit, gen)
            if candidate not in units:
                units.add(candidate)
                queue.append(candidate)
    assert len(units) == 120
    return units


def reference_gaussian_roots():
    for index, sign in product(range(8), (-2, 2)):
        yield tuple(sign if j == index else 0 for j in range(8))
    supports = [tuple(j for i in pairs for j in (2 * i, 2 * i + 1))
                for pairs in combinations(range(4), 2)]
    supports += [tuple(2 * i + k for i, k in enumerate(bits))
                 for bits in product((0, 1), repeat=4) if sum(bits) % 2 == 0]
    for support, signs in product(supports, product((-1, 1), repeat=4)):
        entries = dict(zip(support, signs))
        yield tuple(entries.get(k, 0) for k in range(8))


def derive_frame():
    units = group_units()
    one = (ONE, Z, Z, Z)
    i = (Z, ONE, Z, Z)
    j = (Z, Z, ONE, Z)
    k = (Z, Z, Z, ONE)
    assert all(unit in units and qmul(unit, unit) == scale(-1, one) for unit in (i, j, k))
    roots = sorted(units | {tuple(mul(PHI_CONJ, entry) for entry in unit) for unit in units})
    assert len(roots) == 240 and all(dot(r, r) == 2 for r in roots)
    assert {qmul(i, root) for root in roots} == set(roots)
    first = scale(mul(H, PHI_CONJ), qadd(one, i, j, scale(-1, k)))
    second = scale(mul(H, PHI_CONJ), qadd(one, i, scale(-1, j), k))
    frame = [root for stem in (first, second, k, scale(-1, i)) for root in (stem, qmul(i, stem))]
    assert all(root in roots for root in frame)
    assert all(dot(frame[j], frame[k]) == 2 * (j == k) for j in range(8) for k in range(8))
    mapped = {tuple(dot(root, basis) for basis in frame) for root in roots}
    ref = set(reference_gaussian_roots())
    assert mapped == ref
    print('Oriented icosian axis J: exact native Gaussian chart on all', len(mapped), 'roots.')

    family_axes = [scale(H, qadd(one, scale(-1, i))), scale(H, qadd(j, k)), scale(H, qadd(j, scale(-1,k)))]
    assert all(dot(first, second) == (p == q) for p, first in enumerate(family_axes)
               for q, second in enumerate(family_axes))
    projections = {root: tuple(dot(root, axis) for axis in family_axes) for root in roots}
    d5_roots = [root for root in roots if not any(projections[root])]
    a3_roots = [root for root in roots if sum(x*x for x in projections[root]) == 2]
    assert len(d5_roots) == 40 and len(a3_roots) == 12
    signs = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
    family_weights = [tuple(F(x,2) for x in item) for item in signs]
    buckets = [[root for root in roots if projections[root] == w] for w in family_weights]
    assert all(len(bucket) == 16 for bucket in buckets)
    beta = scale(-1, j)
    beta_projection = projections[beta]
    overlaps = [sum(a*b for a,b in zip(beta_projection, weight)) for weight in family_weights]
    assert overlaps == [-1, 1, 0, 0]
    mark_index = overlaps.index(max(overlaps))
    assert overlaps.count(max(overlaps)) == 1
    mark = family_weights[mark_index]
    retained = [root for n, bucket in enumerate(buckets) if n != mark_index for root in bucket]
    moment_mean = sum(dot(root,root) for root in retained)/len(frame)
    # Orthogonal projector formula is calculated from the original 48 roots.
    def remainder(v):
        return qadd(*(scale(dot(root, v), root) for root in retained), scale(-moment_mean,v))
    def family_projection(v):
        return qadd(*(scale(dot(axis,v), axis) for axis in family_axes))
    marked_vector = qadd(*(scale(value,axis) for value, axis in zip(mark, family_axes)))
    for v in frame:
        assert remainder(v) == qadd(scale(4, family_projection(v)), scale(-16*dot(marked_vector, v), marked_vector))
    rbeta = remainder(beta)
    response = qadd(scale(-1, rbeta), scale(dot(beta,rbeta)/2,beta))
    assert response == scale(8, family_axes[0])
    beta_coords = tuple(dot(beta, basis) for basis in frame)
    response_coords = tuple(dot(response, basis) for basis in frame)
    assert beta_coords == (0,0,0,0,0,2,0,0)
    assert response_coords == (0,0,0,0,0,0,8,8)
    print('Family slice roots:', len(a3_roots), 'orthogonal roots:', len(d5_roots), 'spinor fiber sizes:', [len(bucket) for bucket in buckets])
    def gaussian_coordinates(coordinates):
        return tuple((coordinates[2*d]/2, -coordinates[2*d+1]/2) for d in range(4))

    b, c = gaussian_coordinates(beta_coords), gaussian_coordinates(response_coords)
    print('Native coefficient pairs (real,imaginary): b=(0,-1) on active pair 2; c=(4,-4) on active pair 3.')
    return b, c


if __name__ == '__main__':
    derive_frame()
