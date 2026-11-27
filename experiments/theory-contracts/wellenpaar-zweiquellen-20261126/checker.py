"""One fixed input, independent switches, different specified physical readouts.

Symbolic assertions, finite interval quadrature, and the positive two-graviton
cut in the stated flat vacuum. There is no fit to TFPT constants.
"""

from itertools import product
from math import exp, log, sinh

import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss


def classical_switches():
    u = sp.Matrix([1, sp.I, 0])
    v = sp.conjugate(u)
    U = u * u.T
    V = v * v.T
    a, b, x, y = sp.symbols("a b x y", real=True)
    field = a * U + b * V

    def read(w):
        return [sp.trace(w**2), sp.trace(w**3), sp.trace(w * w.conjugate())]

    data = tuple(sp.expand(expr) for expr in read(field))
    assert data == (8 * a * b, 0, 4 * (a * a + b * b))
    print("Physical switch words (I,J,positive tidal contraction):")
    for left, right in product((0, 1), repeat=2):
        print(
            f"  {left}{right}:",
            tuple(value.subs({a: left, b: right}) for value in data),
        )
    # Order of the classical sources is additive, not ordered multiplication of Weyl matrices.
    assert U + V == V + U and U * V != V * U
    probe = x * U + y * V
    assert (
        sp.expand(sp.trace((U + V) * probe**2 / 2)) == 0
        and sp.expand(sp.trace(probe**2)) == 8 * x * y
    )
    assert sp.hessian(data[0], (a, b)).det() == -64
    assert sp.hessian(data[2], (a, b)).is_positive_definite
    # Degree 5 discriminant pencil of earlier route restricted to this same input.
    eigen = sp.Symbol("lambda")
    assert sp.expand(field.charpoly().as_expr() - eigen * (eigen**2 - 4 * a * b)) == 0
    z = sp.symbols("z")
    # Binary form, up to scale and exchange of the chosen spinor basis.
    assert sp.factor(sp.discriminant(a * z**4 + b, z)) == 256 * a**3 * b**3
    assert sp.I**4 == 1 and sp.I**2 == -1
    print(
        "Only I has a pure bilinear response. J and its target Hessian on the two-switch plane vanish."
    )
    print(
        "Single-input phase loop: four-cycle, sign -1. Equal frequencies freeze marks while pair energy stays positive."
    )


def canonical_pair_measure():
    # Fixed continuum plane wave normalization d^3p/[(2pi)^3 2|p|], no identical-particle divisor yet.
    phase_space = (
        (2 * sp.pi) ** 4 / (2 * sp.pi) ** 6 * (4 * sp.pi) / 4 * sp.Rational(1, 2)
    )
    assert phase_space == 1 / (8 * sp.pi)
    # At zero total spatial momentum, s^4=E^8; common constants suppressed.
    E, t = sp.symbols("E t", positive=True)
    assert sp.integrate(E**8 * sp.exp(-t * E), (E, 0, sp.oo)) == sp.factorial(8) / t**9
    ts = [sp.Integer(1), sp.Integer(2)]
    gram = sp.Matrix([[(a + b) ** -9 for b in ts] for a in ts])
    assert gram.is_positive_definite
    # In a chiral channel identical particles divide full-sphere counting by 2.
    # The conjugate helicity doubles it only for equal form-factor moduli.
    assert 2 * (phase_space / 2) == phase_space
    print(
        "On-shell pair measure ordered 1/(8π), identical same helicity 1/(16π). Curvature-squared weight ∝ |g|²s⁴ and OS time kernel ∝ t^-9."
    )


def massless_pair_kinematics():
    # Distinct projective spinors give the harmonic quartet for two nonzero amplitudes.
    # A single massless Poincare particle instead requires collinear momenta.
    e, f = sp.symbols("E F", positive=True)
    theta = sp.symbols("theta", real=True)
    k = e * sp.Matrix([1, 0, 0, 1])
    q = f * sp.Matrix([1, sp.sin(theta), 0, sp.cos(theta)])
    eta = sp.diag(1, -1, -1, -1)
    mass2 = sp.trigsimp(((k + q).T * eta * (k + q))[0])
    assert sp.simplify(mass2 - 2 * e * f * (1 - sp.cos(theta))) == 0
    assert mass2.subs(theta, 0) == 0 and mass2.subs(theta, sp.pi) == 4 * e * f
    z = sp.Symbol("z")
    a, b = sp.symbols("a b")
    collinear_form = a * z**4 + b * z**4
    assert sp.discriminant(collinear_form, z) == 0
    print(
        "Nonparallel massless inputs have strictly timelike total momentum. Collinear inputs lose the four distinct marks."
    )


def oscillator_sources(area, mode):
    variance = sp.S.One / (2 * area * mode)
    create = sp.zeros(5)
    for n in range(4):
        create[n + 1, n] = sp.sqrt(n + 1)
    q = sp.sqrt(variance) * (create + create.T)
    vacuum = sp.eye(5)[:, 0]
    single = q * vacuum
    pair = (q * q - variance * sp.eye(5)) * vacuum
    assert (
        single.dot(single) == variance
        and sp.simplify(pair.dot(pair) - 2 * variance**2) == 0
    )
    assert single.dot(pair) == 0
    return variance


def one_cone():
    s = sp.symbols("s", positive=True)
    for k in (1, 2, 3):
        oscillator_sources(s, k)
    k = sp.symbols("k", integer=True, positive=True)
    norm = sp.summation(1 / (2 * s**2 * k**2), (k, 1, sp.oo))
    assert sp.simplify(norm - sp.pi**2 / (12 * s**2)) == 0
    print(
        "Paper §4.4: q|0> has norm squared 1/(2sk), centered q²|0> has 1/(2s²k²). Positive pair norm sum π²/(12s²)."
    )


def triangle_rule(order):
    """Gauss integration on each linear half of the unit hat."""
    nodes, weights = leggauss(order)
    u = (nodes + 1) / 2
    return np.concatenate([u - 1, u]), np.concatenate(
        [weights * u / 2, weights * (1 - u) / 2]
    )


def correlation(distance, boost, order=36):
    # Profile 1-|u| for |u|<1. The rapidity dilation acts on the entire right profile, including its width.
    points, weights = triangle_rule(order)
    gaps = (1 + boost) * distance + boost * points[None, :] - points[:, None]
    return -boost * float(np.sum(weights[:, None] * weights[None, :] / gaps**2))


def paper_two_sided():
    # c_perp absorbed in source strengths. The hat is an admissible H^(1/2) smearing, a limit of smooth ones.
    variance = 4 * log(2)
    unboosted = correlation(8, 1)
    boosted = correlation(8, 2)
    assert (
        abs(unboosted - correlation(8, 1, 18)) < 1e-14
        and abs(boosted - correlation(8, 2, 18)) < 1e-14
    )
    assert 0 < abs(boosted) < abs(unboosted) < variance
    assert abs(correlation(1000, 2) / correlation(1000, 1) - 8 / 9) < 1e-7
    # Independent check of the disjoint kernel using the exact antiderivative.
    # Fourier/cosine kernel integral gives -1/d^2, which can be integrated twice per hat by parts.
    import mpmath as mp

    with mp.workdps(45):

        def primitive(u):
            return u**2 * mp.log(abs(u)) / 2 - 3 * u * u / 4

        def exact(d, b):
            return (
                sum(
                    w * v * primitive((1 + b) * d + b * j - i)
                    for i, w in [(-1, 1), (0, -2), (1, 1)]
                    for j, v in [(-1, 1), (0, -2), (1, 1)]
                )
                / b
            )

        assert abs(float(exact(8, 1)) - unboosted) < 1e-14
        assert abs(float(exact(8, 2)) - boosted) < 1e-14
    for boost, cross in [(1, unboosted), (2, boosted)]:
        words = [
            exp(-0.5 * ((i * i + j * j) * variance + 2 * i * j * cross))
            for i, j in product((0, 1), repeat=2)
        ]
        connected = log(words[-1] * words[0] / (words[1] * words[2]))
        assert abs(connected + cross) < 1e-14 and abs(words[1] - 0.25) < 1e-14
        print(
            f"Paper D.17-D.20 boost={boost}, G_LL=G_RR={variance:.10f}, G_LR={cross:.10f}, Zhat[00,01,10,11]={words}, sine product={exp(-variance) * sinh(cross):.10f}"
        )
    # Exact conformal transformation of Eq D.19. t is a log-coordinate; no choice of physical seconds is made.
    A, X, Y = sp.symbols("A X Y", positive=True)
    # X=e^(a x), Y=e^(a y). Multiplying by the two absolute Jacobians gives the half-line kernels.
    cross = -(A**2) * X * Y / (X + Y) ** 2
    same = -(A**2) * X * Y / (X - Y) ** 2
    assert sp.simplify(cross.subs({X: 1, Y: 1}) + A**2 / 4) == 0
    assert (
        sp.simplify(
            cross.subs({X: 1, Y: 2}) / cross.subs({X: 1, Y: 1}) - sp.Rational(8, 9)
        )
        == 0
    )
    assert sp.simplify(same / cross - (X + Y) ** 2 / (X - Y) ** 2) == 0
    print(
        "Opposite half-lines have sech²((x-y)a/2) shape; local rapidity KMS β=2π, β=4π for x=2η."
    )


def reciprocal_control():
    w = sp.symbols("w", real=True)
    # Self-dual part of a candidate Weyl-dependent constitutive law, kappa = i(1 + gamma T).
    block = sp.diag(1 + w, 1 - w, 1)
    differences = [(block**2)[i, i] - 1 for i in (0, 1)]
    assert sp.gcd(*differences) == w
    print(
        "Within the curvature-linear self-dual constitutive ansatz closure survives only for zero Weyl correction."
    )


if __name__ == "__main__":
    classical_switches()
    canonical_pair_measure()
    massless_pair_kinematics()
    one_cone()
    paper_two_sided()
    reciprocal_control()
