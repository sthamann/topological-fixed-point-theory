"""Exact structural comparison and a conditional OPH Lorentz instrument.

Run: python3 -B tfpt_oph_ursprungsregel_20261119.py   (requires SymPy)

The OPH carrier and local terms are model/polynomial_sv/carrier.py and SPEC.md.
The edge-pulse instrument uses a prepared one-load sector, kappa=8, individually
addressed seam pulses, and coherent projection with a resolved branch record.
It does not supply autonomous scheduling or the physical TFPT boundary form.
"""

from collections import Counter
from itertools import product

import sympy as sp


EDGES = {
    (0, 1), (0, 5), (0, 7), (0, 10), (0, 11), (1, 5), (1, 7), (1, 8), (1, 9), (2, 3),
    (2, 4), (2, 6), (2, 10), (2, 11), (3, 4), (3, 6), (3, 8), (3, 9), (4, 5), (4, 9),
    (4, 11), (5, 9), (5, 11), (6, 7), (6, 8), (6, 10), (7, 8), (7, 10), (8, 9), (10, 11),
}


def subgroup_checks():
    unit = (1, 0, 0, 1)
    minus_unit = (4, 0, 0, 4)

    def mul(left, right):
        a, b, c, d = left
        e, f, g, h = right
        return ((a*e + b*g) % 5, (a*f + b*h) % 5, (c*e + d*g) % 5, (c*f + d*h) % 5)

    def inv(x):
        a, b, c, d = x
        return (d, -b % 5, -c % 5, a)

    def order(x):
        value = unit
        for n in range(1, 121):
            value = mul(value, x)
            if value == unit:
                return n
        raise ValueError('Order exceeded number of group elements')

    group = {x for x in product(range(5), repeat=4) if (x[0]*x[3] - x[1]*x[2]) % 5 == 1}
    assert len(group) == 120
    order_four = {x for x in group if order(x) == 4}
    assert len(order_four) == 30
    assert all(mul(x, x) == minus_unit for x in order_four)
    q8_sets = set()
    for x, y in product(order_four, repeat=2):
        if mul(x, y) == mul(minus_unit, mul(y, x)):
            q8_sets.add(frozenset(mul(sign, a) for sign in (unit, minus_unit) for a in (unit, x, y, mul(x, y))))
    assert len(q8_sets) == 5

    def conjugate(x, sub):
        return frozenset(mul(mul(x, y), inv(x)) for y in sub)

    assert {conjugate(x, next(iter(q8_sets))) for x in group} == q8_sets
    expected_orders = {1: 1, 2: 1, 3: 8, 4: 6, 6: 8}
    for sub in q8_sets:
        norm = {x for x in group if conjugate(x, sub) == sub}
        assert Counter(map(order, norm)) == expected_orders
        c3 = next(x for x in norm if order(x) == 3 and any(mul(x, y) != mul(y, x) for y in sub))
        assert {mul(y, x) for x in (unit, c3, mul(c3, c3)) for y in sub} == norm
    return len(q8_sets), expected_orders


def picard_checks():
    eye = sp.eye(2)
    w = (-1 + sp.I*sp.sqrt(3))/2
    sigma = sp.diag(w, sp.conjugate(w))
    mirror = sp.Matrix([[-1, -2], [-1, 1]])/sp.sqrt(3)
    q = [sp.simplify(sigma**k * sp.I * mirror * sigma**(-k)) for k in range(3)]
    metric = sp.diag(1, 2)
    assert all(sp.simplify(x*x) == -eye and sp.simplify(x.H*metric*x) == metric for x in q)
    assert all(sp.simplify(x*y + y*x) == sp.zeros(2) for i, x in enumerate(q) for y in q[i+1:])
    assert any(sp.simplify(q[0]*q[1] - sign*q[2]) == sp.zeros(2) for sign in (-1, 1))
    twist = (eye + q[0])/sp.sqrt(2)
    assert sp.simplify(twist**2) == q[0] and sp.simplify(twist**4) == -eye
    return Counter(sp.simplify(sp.trace(x*sigma**j)) for x in (eye, -eye, *q, *(-a for a in q)) for j in range(3))


def directed_motifs():
    neighbors = {i: {other for edge in EDGES if i in edge for other in edge if other != i} for i in range(12)}
    faces = {tuple(sorted((a, b, c))) for a, b in EDGES for c in neighbors[a] & neighbors[b]}
    assert len(faces) == 20
    # Fix one of the two orientations. Neighboring oriented triangles use
    # opposite directions along their common seam.
    flags = {(0, 1): 7, (1, 7): 0, (7, 0): 1}
    queue = [(0, 1, 7)]
    for a, b, c in queue:
        for i, j, k in ((a, b, c), (b, c, a), (c, a, b)):
            other, = (neighbors[i] & neighbors[j]) - {k}
            if (j, i) in flags:
                assert flags[j, i] == other
                continue
            for x, y, z in ((j, i, other), (i, other, j), (other, j, i)):
                flags[x, y] = z
            queue.append((j, i, other))
    assert len(flags) == 60
    motifs = {}
    for (a, b), c in flags.items():
        d, e = flags[c, b], flags[b, a]
        assert len({a, b, c, d, e}) == 5
        assert all(y in neighbors[x] for x, y in ((a, b), (c, d), (a, c), (b, d), (a, e)))
        motifs[a, b] = (a, b, c, d, e)
    assert motifs[0, 1] == (0, 1, 7, 8, 5)
    return motifs


def instrument(a, b, c, d, e):
    eye = sp.eye(12)
    logical = sp.Matrix.hstack(eye[:, a]-eye[:, b], eye[:, c]-eye[:, d])/sp.sqrt(2)
    assert logical.H*logical == sp.eye(2)
    p = logical*logical.H

    def pulse(u, v):
        assert tuple(sorted((u, v))) in EDGES
        diff = eye[:, u] - eye[:, v]
        # On the edge's two-site part h has eigenvalues 8 and 9. The other
        # ten one-load states have eigenvalue zero for this local term.
        return eye - (1 + sp.I)*diff*diff.T/2

    assert sp.simplify(pulse(a, b)**2 * pulse(c, d)**2) == eye - 2*p
    quarter_z = sp.simplify(logical.H*pulse(a, b)*logical)
    quarter_x = sp.simplify(logical.H*pulse(a, c)*pulse(b, d)*logical)
    assert quarter_z == sp.diag(-sp.I, 1)
    assert quarter_x == sp.Matrix([[1-sp.I, 1+sp.I], [1+sp.I, 1-sp.I]])/2
    for u in (pulse(a, b), pulse(a, c)*pulse(b, d)):
        assert sp.simplify((eye - p)*u*logical) == sp.zeros(12, 2)
    assert sp.simplify(logical.H*pulse(a, b)*pulse(c, d)*logical) == -sp.I*sp.eye(2)
    k = sp.simplify(logical.H*pulse(a, e)*logical)
    lost = sp.simplify((eye-p)*pulse(a, e)*logical)
    assert k == sp.diag((3-sp.I)/4, 1)
    assert sp.simplify(lost.H*lost + k.H*k) == sp.eye(2)
    positive = sp.simplify(k.H*k)
    assert positive == sp.diag(sp.Rational(5, 8), 1)
    swap = quarter_x**2
    assert sp.simplify(swap*positive*swap*positive) == sp.Rational(5, 8)*sp.eye(2)
    return positive


def lorentz_certificate(positive):
    q = positive[0, 0]
    z = positive/sp.sqrt(q)
    # A unitary quarter turn of the finite octahedral subgroup transports z to x.
    x = sp.Matrix([[1, 1], [1, -1]])*z*sp.Matrix([[1, 1], [1, -1]])/2
    commutator = sp.simplify(z*x*z.inv()*x.inv())
    assert commutator.det() == 1
    trace = sp.simplify(sp.trace(commutator))
    assert trace == sp.Rational(12719, 6400)
    assert -2 < trace < 2 and trace.is_Integer is False
    # Similarity and projective normalization leave (tr M)^2/det M unchanged.
    # A product of Gaussian-rational pulse and projection blocks cannot have
    # projective order five on the rational two-mode code.
    order_five_invariants = tuple(sp.simplify(4*sp.cos(k*sp.pi/5)**2) for k in (1, 2))
    assert {sp.simplify(t - sp.Rational(3, 2)) for t in order_five_invariants} == {-sp.sqrt(5)/2, sp.sqrt(5)/2}
    # Rational nonintegral trace rules out a root of unity. The accompanying
    # Lie algebra argument upgrades the resulting elliptic circle to SL(2,C).
    return q, trace


def spinor_branching():
    """Branch the OPH port Spin(16) lift and recover its missing doublet."""
    # Class order and characters from the exact audited OPH McKay receipt
    # code/a5_closure/receipts/external/oph_sl2f5_mckay_e8_4ae2148a.json.
    sizes = (1, 1, 20, 30, 12, 12, 20, 12, 12)
    u = (sp.sqrt(5)-1)/2
    v = (-sp.sqrt(5)-1)/2
    two = (2, -2, -1, 0, u, v, 1, -u, -v)
    two_conjugate = (2, -2, -1, 0, v, u, 1, -v, -u)
    four = (4, -4, 1, 0, -1, -1, -1, 1, 1)
    six = (6, -6, 0, 0, 1, 1, 0, -1, -1)
    spinors = (two, two_conjugate, four, six)
    three = (3, 3, 0, -1, -u, -v, 0, -u, -v)
    three_conjugate = (3, 3, 0, -1, -v, -u, 0, -v, -u)

    def inner(a, b):
        return sp.simplify(sum(n*x*y for n, x, y in zip(sizes, a, b))/sp.Integer(120))

    assert all(inner(a, b) == (1 if i == j else 0)
               for i, a in enumerate(spinors) for j, b in enumerate(spinors))
    # The vector action is the real twelve-port permutation module plus four
    # fixed coordinates. Rotation sectors: 6 swaps, 4 triples, 2 quintuplets.
    # A fixed plane makes the two half-spin characters equal. Odd-order
    # representatives are lifted with odd order, which fixes their sign.
    def half_spin_character(angles):
        return sp.simplify(128*sp.prod(sp.cos(sp.pi*a) for a in angles))

    chi_4 = half_spin_character((sp.Rational(1, 2),)*3)
    chi_3 = half_spin_character((sp.Rational(1, 3),)*4)
    chi_5 = half_spin_character((sp.Rational(1, 5), sp.Rational(2, 5))*2)
    assert chi_4 == 0 and chi_3 == chi_5 == 8
    half_spin = (128, -128, chi_3, chi_4, chi_5, chi_5, -chi_3, -chi_5, -chi_5)

    def multiply(a, b):
        return tuple(sp.simplify(x*y) for x, y in zip(a, b))

    mult = tuple(inner(half_spin, irrep) for irrep in spinors)
    selected = tuple(inner(multiply(four, three), irrep) for irrep in spinors)
    conjugate = tuple(inner(multiply(four, three_conjugate), irrep) for irrep in spinors)
    assert mult == (0, 0, 8, 16)
    assert selected == (1, 0, 1, 1) and conjugate == (0, 1, 1, 1)
    assert sum(n*character[0] for n, character in zip(mult, spinors)) == 128

    # Averaging the projectors onto the stationary component of each
    # order-three face rotation is a positive central operator. Its zero
    # eigenspace in the coupled quartet and triplet is the single doublet.
    face_penalty = tuple(sp.Rational(rep[0] + 2*rep[2], 3*rep[0]) for rep in spinors)
    assert face_penalty == (0, 0, sp.Rational(1, 2), sp.Rational(1, 3))
    assert sum(n*rep[0] for n, rep, energy in zip(selected, spinors, face_penalty) if energy == 0) == 2
    return mult, selected, conjugate, face_penalty


if __name__ == '__main__':
    count, element_orders = subgroup_checks()
    traces = picard_checks()
    assert traces == {2: 1, -2: 1, 0: 6, 1: 8, -1: 8}
    # These are exactly the characters of either OPH doublet on the order
    # classes 1,2,3,4,6 of each of the five binary tetrahedral subgroups.
    motifs = directed_motifs()
    q, trace = lorentz_certificate(instrument(*motifs[0, 1]))
    assert all(instrument(*motif) == sp.diag(q, 1) for motif in motifs.values())
    mult, selected, conjugate, face_penalty = spinor_branching()
    print('Five binary tetrahedral subgroups in SL(2,F5); element orders:', element_orders)
    print('Picard compact spinor trace distribution:', dict(traces))
    print('Exact OPH motifs with the same instrument:', len(motifs))
    print('Positive filter:', q, '; state-independent inverse-round-trip probability:', q**2)
    print('Exact elliptic commutator trace:', trace)
    print('Multiplicities in either half-spinor of the port lift (2, 2 prime, 4, 6):', mult)
    print('Multiplicities in the McKay quartet tensored with each port triplet:', selected, conjugate)
    print('Averaged face-invariant projector on the four spinorial irreps:', face_penalty)
