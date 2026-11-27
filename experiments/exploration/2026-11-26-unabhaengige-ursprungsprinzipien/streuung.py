"""Exact local probes of the spin-2 pencil and its projective response."""

from itertools import product
import sympy as sp


def inner(left, right):
    return (left * right).trace()


def stf(matrix):
    return matrix - sp.eye(3) * matrix.trace() / 3


def mixed_components(i, j):
    matrix = sp.zeros(3)
    matrix[i, j] = matrix[j, i] = 1
    return matrix


BASIS = [sp.diag(1, -1, 0), sp.diag(1, 1, -2)] + [
    mixed_components(i, j) for i, j in [(0, 1), (0, 2), (1, 2)]
]
T = BASIS[0]


def weyl_pulses_and_response():
    # The self-dual Weyl matrix of an elementary type N pulse has the form uu^T.
    u = sp.Matrix([1, sp.I, 0])
    v = sp.Matrix([1, -sp.I, 0])
    w = sp.Matrix([1, 0, sp.I])
    a, b, e = sp.symbols("a b e")
    U, V, W = [n * n.T for n in (u, v, w)]
    assert all(n * n == sp.zeros(3) for n in (U, V, W))
    t = a * U + b * V
    assert t.det() == 0 and sp.expand(inner(t, t)) == 8 * a * b
    assert U / 2 + V / 2 == T and (U + V + W).det() != 0
    xx = sp.symbols("x:5")
    x = sum((c * m for c, m in zip(xx, BASIS)), sp.zeros(3))

    def coefficient(polynomial, order):
        return sp.expand(polynomial).coeff(e, order)

    path = T + e * x
    assert sp.expand(coefficient((path**2).trace(), 2) - inner(x, x)) == 0
    assert sp.expand(coefficient((path**3).trace(), 2) - 3 * inner(T, x * x)) == 0
    # Reading also the two first derivatives restricts x to zero diagonal.
    assert (
        sp.expand(inner(T, x)) == 2 * xx[0] and sp.expand(inner(T * T, x)) == 2 * xx[1]
    )
    aa, bb, cc = xx[2:]
    off = x.subs({xx[0]: 0, xx[1]: 0})
    assert sp.expand(inner(off, off)) == 2 * (aa**2 + bb**2 + cc**2)
    assert sp.expand(inner(T, off * off)) == bb**2 - cc**2
    for sign_a, sign_c in product([-1, 1], repeat=2):
        sample = off.subs({aa: sign_a * sp.I * sp.sqrt(2), bb: 1, cc: sign_c})
        assert (
            sp.simplify(inner(sample, sample)) == 0
            and inner(T, sample * sample) == 0
            and sample.det() != 0
        )
    # A product invariant responds in 2nd order even at a point on the curve.
    assert coefficient(((T + e * U) ** 2).trace() ** 2, 2) == 16
    print(
        "Two nonparallel elementary Weyl pulses give J=0, I nonzero; three distinct pulses can give J nonzero."
    )
    print(
        "Hessians of tr(W²), tr(W³) reproduce the pencil. Including both linear responses leaves four projective points of full matrix rank; the cubic term distinguishes them all."
    )


def time_symmetric_vacuum():
    # Time-symmetric local data g = h^4 delta, K = 0. Harmonic h solves the constraint.
    x, y, z = sp.symbols("x y z", real=True)
    h = 1 + x * x - y * y
    grad = sp.Matrix([sp.diff(h, q) for q in (x, y, z)])
    hessian = sp.hessian(h, (x, y, z))
    assert hessian.trace() == 0
    ricci = (
        -2 * hessian / h + (6 * grad * grad.T - 2 * sp.eye(3) * grad.dot(grad)) / h**2
    )
    assert sp.simplify(ricci.trace()) == 0
    assert ricci.subs({x: 0, y: 0, z: 0}) == -4 * T
    relative = ricci.subs({x: sp.S.Half, y: 0, z: 0})
    assert relative == sp.diag(-16, 48, -32) / 25
    print(
        "Exact local vacuum data h=1+x²−y²: electric Weyl has J=0 at the origin and J≠0 at x=1/2."
    )


def soft_response():
    # Equal energies in the hard center of momentum. For the two soft rays
    # perpendicular to the scattering plane all four denominators coincide.
    t = sp.symbols("t", real=True)
    a = sp.Matrix([1, 0, 0])
    b = sp.Matrix([(1 - t**2) / (1 + t**2), 2 * t / (1 + t**2), 0])
    soft = b * b.T - a * a.T
    assert sp.simplify(soft.trace()) == 0 and soft[:, 2] == sp.zeros(3, 1)
    assert sp.simplify(inner(soft, soft) - 8 * t**2 / (1 + t**2) ** 2) == 0
    assert soft.subs(t, 0) == sp.zeros(3)
    print(
        "Soft tensor: eigenvalues +s,-s,0 for every noncollinear angle; forward response vanishes."
    )


def pencil_and_curve():
    lam, z, w = sp.symbols("lambda z w")
    gram = sp.Matrix([[inner(a, b) for b in BASIS] for a in BASIS])
    assert gram.is_positive_definite
    bilinear = sp.Matrix(
        [[inner(T, (a * b + b * a) / 2) for b in BASIS] for a in BASIS]
    )
    operator = gram.inv() * bilinear
    assert (
        sp.simplify(
            operator.charpoly(lam).as_expr()
            - lam * (lam**2 - sp.Rational(1, 4)) * (lam**2 - sp.Rational(1, 3))
        )
        == 0
    )
    assert all(
        (operator - l * sp.eye(5)).rank() == 4
        for l in [0, sp.S.Half, -sp.S.Half, 1 / sp.sqrt(3), -1 / sp.sqrt(3)]
    )
    m = sp.Matrix([w * w - z * z, sp.I * (w * w + z * z), 2 * z * w])
    point = m * m.T
    assert sp.expand((m.T * m)[0]) == 0 and (point * point).applyfunc(
        sp.expand
    ) == sp.zeros(3)
    quartics = [sp.Poly(inner(b, point), z, w) for b in BASIS]
    assert (
        sp.Matrix(
            [
                [q.coeff_monomial(z**i * w ** (4 - i)) for i in range(5)]
                for q in quartics
            ]
        ).rank()
        == 5
    )
    qt = sp.expand(inner(T, point))
    assert qt == 2 * (w**4 + z**4)
    assert sp.degree(sp.gcd(qt.subs(w, 1), sp.diff(qt.subs(w, 1), z)), z) == 0
    # Four roots modulo a common rotation, with a harmonic ordering.
    a, b, c, d = 1, -1, sp.I, -sp.I
    ratio = sp.simplify((a - c) * (b - d) / ((a - d) * (b - c)))
    j = 256 * (1 - ratio + ratio**2) ** 3 / (ratio**2 * (1 - ratio) ** 2)
    assert ratio == -1 and j == 1728
    circular = T + sp.I * BASIS[2]
    assert sp.expand(inner(circular, point)) == 4 * z**4
    print(
        "Pencil roots 0, ±1/2, ±1/√3 are simple. Embedded rational normal quartic; four simple marks. Single helicity gives one root of order four."
    )
    # The first excluded scalar (second order in T) is generically independent.
    # An explicit point of S with a nonzero extra response.
    probe = BASIS[2] + sp.I * BASIS[1] / sp.sqrt(3)
    assert (
        sp.simplify(inner(probe, probe)) == 0 and sp.simplify(inner(T, probe**2)) == 0
    )
    assert sp.simplify(inner(T**2, probe**2)) == sp.Rational(4, 3)
    print(
        "Second-order scalar tr(T² X²) on a point of S: 4/3. Hermitian norm has no nonzero complex null vectors."
    )


def base_scheme():
    # Rotate T by 45 degrees so the marked roots in this chart are ±1, ±i.
    t = BASIS[2]
    coeff = sp.symbols("y:5")
    x = sum((v * b for v, b in zip(coeff, BASIS)), sp.zeros(3))
    h = inner(t, x)
    responses = sp.Matrix(list(x * x) + list(h * x))
    jacobian = responses.jacobian(coeff)
    for z in [1, -1, sp.I, -sp.I]:
        m = sp.Matrix([1 - z * z, sp.I * (1 + z * z), 2 * z])
        point = m * m.T
        at = {v: inner(b, point) / inner(b, b) for v, b in zip(coeff, BASIS)}
        assert responses.subs(at) == sp.zeros(18, 1)
        assert jacobian.subs(at).rank() == 4  # only radial tangent in affine cone
    print(
        "Base ideal of the quadratic tensor readout: all four marked points reduced (affine-cone Jacobian rank four)."
    )


def covariants_and_symmetry():
    # Independent first-order-in-T tensor covariants; CG gives multiplicity 3.
    def channels(x, t):
        c = stf(x * x)
        return [inner(x, x) * t, stf(t * c + c * t), inner(t, x) * x]

    rows = []
    for x in BASIS + [BASIS[0] + BASIS[3]]:
        for t in BASIS:
            vals = channels(x, t)
            rows.extend(
                [[y[i, j] for y in vals] for i, j in product(range(3), repeat=2)]
            )
    assert sp.Matrix(rows).rank() == 3
    # r exchanges the two principal axes; s is the half-turn through two marks.
    r = sp.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    s = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
    rotations = [sp.eye(3), r * r, r, s, r * s]

    def character(subspace, twist):
        out = []
        for g in rotations:
            line_character = inner(T, g * T * g.T) / inner(T, T)
            out.append(
                sum(inner(b, g * b * g.T) / inner(b, b) for b in subspace)
                * (line_character if twist else 1)
            )
        return out

    assert character(BASIS[2:], False) == [3, -1, -1, 1, -1]
    assert character(BASIS[2:], True) == [3, -1, 1, -1, -1]
    assert character(BASIS, True) == [5, 1, 1, -1, 1]
    print(
        "D4 classes (1,r²,r,s,rs): triplet bare (3,-1,-1,1,-1), rational sections (3,-1,1,-1,-1). Full quotient character (5,1,1,-1,1)."
    )


def picard_glue():
    # h,e1,...,e5,f1,...,f4 in a plane blowup marking. Marking is only a basis.
    eta = sp.diag(1, *([-1] * 9))
    e = [sp.eye(10)[:, j] for j in range(10)]
    k = -3 * e[0] + sum(e[1:], sp.zeros(10, 1))
    d = [e[i] - e[i + 1] for i in range(1, 5)] + [e[0] - e[1] - e[2] - e[3]]
    a = [e[i] - e[i + 1] for i in range(6, 9)]
    mix = e[5] - e[6]
    roots = sp.Matrix.hstack(*d, *a)
    gram = roots.T * eta * roots
    assert (
        gram[:5, :5].det() == -4
        and gram[5:, 5:].det() == -4
        and gram[:5, 5:] == sp.zeros(5, 3)
    )
    assert (
        (k.T * eta * k)[0] == 0
        and k.T * eta * roots == sp.zeros(1, 8)
        and (k.T * eta * mix)[0] == 0
    )
    assert (mix.T * eta * mix)[0] == -2
    coordinates = sp.Matrix.hstack(roots, k).gauss_jordan_solve(mix)[0]
    assert sp.ilcm(*[v.q for v in coordinates]) == 4
    # Use the f4-section of the null lattice as an integral lift of the quotient.
    # h-e1-e2-e3, e1-e2,...,e7-e8 give a negative E8 Cartan matrix.
    e8 = sp.Matrix.hstack(
        e[0] - e[1] - e[2] - e[3], *[e[i] - e[i + 1] for i in range(1, 8)]
    )
    assert (e8.T * eta * e8).det() == 1 and k.T * eta * e8 == sp.zeros(1, 8)
    frame = sp.Matrix.hstack(e8, k, e[9])
    assert abs(frame.det()) == 1
    print(
        "Picard calculation: det(D5⊕A3)=16, mixed root square -2 and order 4; quotient E8(-1), determinant 1."
    )


if __name__ == "__main__":
    weyl_pulses_and_response()
    time_symmetric_vacuum()
    soft_response()
    pencil_and_curve()
    base_scheme()
    covariants_and_symmetry()
    picard_glue()
