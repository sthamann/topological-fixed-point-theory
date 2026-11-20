"""Finite exact checks for isotropic observable-generated precession.

Requires Python >=3.9 and SymPy. No fitted inputs, files, network or random seed.
The tested state ball is an operational model; its cone is an internal order cone.
"""
from itertools import combinations
from math import comb

import sympy as sp


def hamiltonian_intertwiners(d):
    """Dimension of Hom_SO(d)(V, wedge^2 V), via sign flips and quarter turns."""
    triples = [(i, j, k) for i in range(d) for j, k in combinations(range(d), 2)]
    planes = list(combinations(range(d), 2))
    # Every even sign flip must fix a surviving component.
    survivors = [t for t in triples
                 if all(sum(index in plane for index in t) % 2 == 0 for plane in planes)]
    if not survivors:
        return 0
    positions = {t: p for p, t in enumerate(survivors)}
    rows = []
    for a, b in planes:
        image = list(range(d))
        image[a], image[b] = b, a
        for t in survivors:
            i, j, k = (image[p] for p in t)
            sign = (-1)**t.count(b)
            if j > k:
                j, k, sign = k, j, -sign
            row = [0] * len(survivors)
            row[positions[t]] += 1
            row[positions[(i, j, k)]] -= sign
            rows.append(row)
    # The finite subgroup bounds the continuous invariants from above. For d=3
    # the Levi-Civita tensor attains that bound, being fixed by all of SO(3).
    return len(survivors) - sp.Matrix(rows).rank()


def half_turn(axis):
    return 2 * axis * axis.T / axis.dot(axis) - sp.eye(3)


def flatten(matrix):
    return matrix.reshape(9, 1)


def skew(axis):
    a, b, c = axis
    return sp.Matrix([[0, -c, b], [c, 0, -a], [-b, a, 0]])


def symmetric_pair(u, v):
    return (u * v.T + v * u.T) / sp.sqrt(2)


def incidence_bijections(adjacency, target):
    """Enumerate all graph isomorphisms of two tiny unlabelled graphs."""
    size = len(adjacency)
    chosen = []

    def extend():
        i = len(chosen)
        if i == size:
            yield tuple(chosen)
            return
        for j in range(size):
            if j in chosen:
                continue
            if all(adjacency[i][k] == target[j][image] for k, image in enumerate(chosen)):
                chosen.append(j)
                yield from extend()
                chosen.pop()

    yield from extend()


def finite_pulse_calibration():
    """Identify finite quadratic-response calibration and enumerate all pulse pairs."""
    from collections import Counter
    from itertools import permutations

    # Finite SO(3) classification leaves cyclic/dihedral (fixed axial quadrupole)
    # and three polyhedral groups. The spin-two rotation character on orders
    # 1,2,3,4,5 is 5,1,-1,-1,0 for both primitive fifth-turn angle classes.
    orders = {"A4": {1:1, 2:3, 3:8}, "S4": {1:1, 2:9, 3:8, 4:6},
              "A5": {1:1, 2:15, 3:20, 5:24}}
    chi = {1:5, 2:1, 3:-1, 4:-1, 5:0}
    norms = {g: sp.Rational(sum(count * chi[q]**2 for q, count in counts.items()), sum(counts.values()))
             for g, counts in orders.items()}
    assert list(norms.values()) == [3, 2, 1]
    identity = tuple(range(5))

    def compose(a, b):
        return tuple(a[b[i]] for i in identity)

    def order(a):
        exponent, power = 1, a
        while power != identity:
            exponent += 1
            power = compose(power, a)
        return exponent

    # Coordinate realization of the previously selected rotation group. The
    # five output slots are re-enumerated as commuting contexts, independently
    # of whether its permutation-coordinate labels are later used.
    group = [p for p in permutations(identity)
             if sum(p[i]>p[j] for i,j in combinations(identity,2)) % 2 == 0]
    pulses = [a for a in group if order(a) == 2]
    contexts = {frozenset(b for b in pulses if compose(a,b)==compose(b,a)) for a in pulses}
    assert len(contexts) == 5 and all(len(c)==3 for c in contexts)
    pair_count = Counter()
    shortest = []
    for a in pulses:
        for b in pulses:
            if a == b:
                continue
            ab, ba = compose(a,b), compose(b,a)
            k = order(ab)
            pair_count[k] += 1
            if ab != ba:
                shortest.append((k, ab, ba))
    minimum = min(k for k, _, _ in shortest)
    marks = set()
    for k, ab, ba in shortest:
        if k != minimum:
            continue
        fixed = frozenset(c for c in contexts
                          if frozenset(compose(compose(ab,h),ba) for h in c) == c)
        assert len(fixed) == 2
        marks.add(fixed)
    assert pair_count == {2:30, 3:60, 5:120} and len(marks)==10
    print("Finite polyhedral quadratic character norms:", list(norms.values()), "(A4,S4,A5)")
    print("Maximal commuting pulse contexts:", len(contexts), "; pair orders:", dict(sorted(pair_count.items())))
    print("Shortest noncommuting recurrence:", minimum, "; fixed/moved:", (2,3), "; unordered marks:",len(marks))


def dodecahedral_control():
    """Additional sampling control only: isotropy does not select this finite orbit."""
    from itertools import product

    phi = (1 + sp.sqrt(5)) / 2
    axes = [sp.Matrix([1, a, b]) for a, b in product((-1, 1), repeat=2)]
    for sign in (-1, 1):
        coords = [0, phi - 1, sign * phi]
        for k in range(3):
            axes.append(sp.Matrix(coords[k:] + coords[:k]))
    # Unit traceless quadrupoles have Gram (3*(n.m)^2-1)/2 for unit n,m.
    assert all(sp.simplify(n.dot(n)) == 3 for n in axes)
    measured = sp.Matrix(10, 10, lambda i, j: sp.simplify(axes[i].dot(axes[j])**2 / 6 - sp.S.Half))
    third = sp.Rational(1, 3)
    adjacency = [[measured[i, j] == third for j in range(10)] for i in range(10)]
    pairs = [set(p) for p in combinations(range(5), 2)]
    target = [[p.isdisjoint(q) for q in pairs] for p in pairs]
    maps = list(incidence_bijections(adjacency, target))
    assert len(maps) == 120
    labels = maps[0]
    # Original Q4 Gram in the published labelled Petersen dictionary.
    # b_ij=pair indicator; c=(4-sqrt(6))/10. Normalized Gram = (2*|p∩q|-1)/3.
    native = sp.Matrix(10, 10, lambda i, j: (2 * len(pairs[labels[i]] & pairs[labels[j]]) - 1) * third)
    identity = sp.eye(10)
    assert measured + native == 2*identity
    assert measured * native == sp.zeros(10)
    assert measured**2 == 2*measured and native**2 == 2*native
    triangles = list(combinations(range(10), 3))
    def triple(G, i, j, k):
        return G[i, j] * G[j, k] * G[k, i]
    wrong = sum(triple(measured,*t) == -triple(native,*t) for t in triangles)
    assert wrong == 120
    # Forgetting the dictionary, signed relabellings exist. For a fixed first ray,
    # triple signs with that ray characterize a switching class on all other rays.
    src = list(range(1, 10))
    source_links = [[i != j and triple(measured, 0, i, j) > 0 for j in src] for i in src]
    matches = 0
    for new_root in range(10):
        dst = [i for i in range(10) if i != new_root]
        dest_links = [[i != j and triple(native, new_root, i, j) > 0 for j in dst] for i in dst]
        matches += sum(1 for _ in incidence_bijections(source_links, dest_links))
    assert matches > 0
    print("Dodecahedral test against Q4 Gram:", "sum=2I, product=0; reversed triangle signs:", wrong)
    print("Petersen-preserving labelings:", len(maps), "; matching after signed unrestricted relabeling:", matches)


def main():
    dimensions = {d: hamiltonian_intertwiners(d) for d in range(2, 13)}
    assert list(dimensions.values()) == [0, 1] + [0] * 9
    print("SO(d) invariant precessions for d=2..12:", list(dimensions.values()))

    axes = list(sp.eye(3).columnspace())
    samples = axes + [u + v for u, v in combinations(axes, 2)]
    contrasts = [half_turn(n) + sp.eye(3) / 3 for n in samples]
    symmetric_basis = sp.Matrix.hstack(*map(flatten, contrasts)).columnspace()
    assert len(symmetric_basis) == 5
    basis = [v.reshape(3, 3) for v in symmetric_basis]
    turn = half_turn(axes[2])
    plus = sp.Matrix.hstack(*(flatten(A + turn*A*turn) for A in basis)).rank()
    minus = sp.Matrix.hstack(*(flatten(A - turn*A*turn) for A in basis)).rank()
    assert (plus, minus) == (3, 2)
    print("Centered half-turn span:", len(basis), "; marked ranks:", (plus, minus))

    skew_rank = sp.Matrix.hstack(*(flatten(skew(n)) for n in axes)).rank()
    assert skew_rank == 3
    z = sp.Symbol("z")
    ranks = {}
    angle_hits, total_angles = 0, 0
    for q in range(2, 13):
        phi = sp.Poly(sp.cyclotomic_poly(q, z), z)
        # 1-cos(theta) is zero iff (z-1)^2=0, sin(theta) iff z^2-1=0.
        has_symmetric = not sp.Poly((z-1)**2, z).rem(phi).is_zero
        has_skew = not sp.Poly(z**2-1, z).rem(phi).is_zero
        ranks[q] = len(basis) * has_symmetric + skew_rank * has_skew
        total_angles += phi.degree()
        angle_hits += phi.degree() * (ranks[q] == 5)
    assert list(ranks.values()) == [5] + [8] * 10
    print("Centered primitive rotation spans for q=2..12:", list(ranks.values()))
    print("Five-dimensional span: 1/11 order symbols;", f"{angle_hits}/{total_angles} primitive angles.")

    # All diagonal sign elements with three labelled axes, including improper ones.
    outcomes = {}
    for p in range(4):
        r_even = comb(p + 1, 2) + comb(4 - p, 2) - 1
        r_odd = p*(3 - p)
        outcomes[(r_even, r_odd)] = outcomes.get((r_even, r_odd), 0) + comb(3, p)
    assert outcomes == {(5, 0): 2, (3, 2): 6}
    print("Labeled sign-involution control (3+2):", "6/8")

    theta = sp.Symbol("theta", real=True)
    s, c = sp.sin(theta), sp.cos(theta)
    dx, dy, dz = axes
    n = c*dz + s*dx
    t = c*dx - s*dz
    reference = [symmetric_pair(dz, dx), symmetric_pair(dz, dy)]
    comparison = [symmetric_pair(n, t), symmetric_pair(n, dy)]
    gram = sp.Matrix([[(a*b).trace().trigsimp() for b in comparison] for a in reference])
    assert (gram - sp.diag(sp.cos(2*theta), c)).applyfunc(sp.trigsimp).is_zero_matrix
    x = sp.Symbol("x", real=True)
    poly = (2*x-1)**2
    # Primary record: TFPT-CLOCK-WINDUNG-PICARD-20261111.txt section 14,
    # standard Hermitian pairing and the third/fourth coordinate weak plane.
    picard_basis = sp.Matrix([[1, 1], [1, -1], [1, sp.I], [1, -sp.I], [0, 0]]) / 2
    weak_basis = sp.eye(5)[:, 2:4]
    overlap = weak_basis.H * picard_basis
    picard = sp.Rational(1, 2)
    assert picard_basis.H * picard_basis == sp.eye(2)
    assert overlap * overlap.H == picard * sp.eye(2)
    assert sp.gcd(x-picard, poly-picard) == 1
    print("Squared principal cosines:", "{x, (2*x-1)^2}; x=cos(theta)^2")
    print("Picard/weak common orbit at {1/2,1/2}: impossible (exact gcd=1)")
    finite_pulse_calibration()
    dodecahedral_control()


if __name__ == "__main__":
    main()
