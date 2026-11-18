#!/usr/bin/env python3
"""Ninth pasted note: circle motions and the marks (18 November 2026).

Sections
    T1  lift of a motion b(x) d/dx to the cover near a mark: pole unless b(m) = 0; regular if the mark moves along
    T2  the rotation family y^2 = x^4 - a: its transport over a quarter turn is hQ, not Q
    T3  local rate at a mark: half the rate of the sphere; an eighth turn per quarter turn
    T4  weight of a mark in the free source: quarter exponents give 1/32 + 1/32 = 1/16 per spinor
    T5  the four-mark function of the sixteen-site (c = 8) source: degree and Moebius weight
"""
import itertools
import sys

import mpmath as mp
import sympy as sp

FAIL = []


def check(name, ok, info=None):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info is not None else ""))
    if not ok:
        FAIL.append(name)


# ---------------------------------------------------------------- T1
print("T1  lifting a motion of the sphere")
v, m, b0, b1, b2 = sp.symbols("v m b0 b1 b2")
b = lambda x: b0+b1*x+b2*x**2                    # generator of a Moebius motion
fixed = sp.series(b(m+v**2)/(2*v), v, 0, 2).removeO()
check("fixed mark: the lifted speed b(m + v^2)/(2v) has the pole b(m)/(2v)",
      sp.simplify(fixed.coeff(v, -1)-b(m)/2) == 0)
check("a non-zero generator has at most two zeros: it cannot vanish at four marks",
      sp.degree(b(sp.Symbol("x")), sp.Symbol("x")) == 2)
moving = sp.expand(sp.simplify((b(m+v**2)-b(m))/(2*v)))
check("mark moving with the motion: lifted speed = b'(m) v/2 + b''(m) v^3/4, regular",
      sp.simplify(moving-(sp.diff(b(m), m)*v/2+sp.diff(b(m), m, 2)*v**3/4)) == 0)

# ---------------------------------------------------------------- T2
print("T2  the rotation family")
x, y, a, th = sp.symbols("x y a theta")
check("(x, y) -> (e^{i theta} x, e^{2 i theta} y) carries y^2 = x^4 - a to y^2 = x^4 - e^{4 i theta} a",
      sp.simplify((sp.exp(2*sp.I*th)*y)**2-((sp.exp(sp.I*th)*x)**4-sp.exp(4*sp.I*th)*a)
                  - sp.exp(4*sp.I*th)*(y**2-(x**4-a))) == 0)
qt = {th: sp.pi/2}
img = (sp.simplify((sp.exp(sp.I*th)*x).subs(qt)), sp.simplify((sp.exp(2*sp.I*th)*y).subs(qt)))
check("after a quarter turn the curve is the same and the transport is (x, y) -> (i x, -y) = h Q",
      img == (sp.I*x, -y), img)
print("      the recorded clock lift is Q = (i x, y), which fixes the two points over x = 0;")
print("      h Q fixes the two points over x = infinity.  The record relates them by S Q S^-1 = h Q^-1.")

# ---------------------------------------------------------------- T3
print("T3  local rate at a mark")
rate = sp.diff(sp.I*m, m)/2                      # b(x) = i x: rotation
check("rotation b(x) = i x: the local coordinate at a moving mark turns at half the rate, i/2",
      rate == sp.I/2)
check("quarter turn of the sphere: the local coordinate turns by exp(i pi/4), an eighth "
      "(the corner phase of the fifteenth round)", sp.simplify(sp.exp(sp.I/2*sp.pi/2)-sp.exp(sp.I*sp.pi/4)) == 0)

# ---------------------------------------------------------------- T4
print("T4  weight of a mark in the free source")
nu = sp.Rational(1, 4)
check("sectors with exponents -1/4 and +1/4: twist weights nu^2/2 each, 1/32 + 1/32 = 1/16 per mark and spinor",
      nu**2/2+nu**2/2 == sp.Rational(1, 16))
check("sixteen real components (eight spinors): 8/16 = 1/2 per mark, total 2 for four marks",
      8*sp.Rational(1, 16) == sp.Rational(1, 2))

# ---------------------------------------------------------------- T5
print("T5  the four-mark function for central charge eight")
ms = sp.symbols("m1:5")
X = sp.Symbol("X")
poly = sp.Poly(sp.expand(sp.prod([X-mi for mi in ms])), X)
c4, c3, c2, c1, c0 = poly.all_coeffs()
Iinv = sp.expand(12*c4*c0-3*c3*c1+c2**2)
disc = sp.prod([(ms[i]-ms[j])**2 for i, j in itertools.combinations(range(4), 2)])
lam = sp.Symbol("lambda", positive=True)
sc = {mi: lam*mi for mi in ms}
check("I has degree 4 and the discriminant degree 12 in the marks",
      sp.simplify(Iinv.subs(sc, simultaneous=True)/Iinv-lam**4) == 0
      and sp.simplify(disc.subs(sc, simultaneous=True)/disc-lam**12) == 0)
check("I / Disc^(1/2) has degree -2 = -4 h with h = 1/2: the weight of T4",
      sp.Rational(4)-sp.Rational(12, 2) == -4*sp.Rational(1, 2))
mp.mp.dps = 25
pts = [mp.mpc("0.3", "0.2"), mp.mpc("-0.7", "0.4"), mp.mpc("1.1", "-0.5"), mp.mpc("0.2", "-0.9")]
al, be, ga, de = mp.mpc("1.2", "0.3"), mp.mpc("-0.4", "0.7"), mp.mpc("0.5", "-0.2"), None
de = (1+be*ga)/al
phi = lambda t: (al*t+be)/(ga*t+de)
dphi = lambda t: 1/(ga*t+de)**2


def G(p):
    cf = mp.taylor(lambda t: mp.fprod([t-pi for pi in p]), 0, 4)
    e0, e1, e2, e3, e4 = cf
    Iv = 12*e4*e0-3*e3*e1+e2**2
    dv = mp.fprod([(p[i]-p[j])**2 for i, j in itertools.combinations(range(4), 2)])
    return Iv**2/dv                               # square of I/Disc^(1/2), free of branch choices


lhs = G([phi(t) for t in pts])
rhs = G(pts)*mp.fprod([dphi(t) for t in pts])**(-1)
check("under a Moebius map, (I/Disc^(1/2))^2 picks up the product of phi'(m_k)^(-1): weight 1/2 at each mark",
      abs(lhs-rhs)/abs(rhs) < mp.mpf(10)**(-18), mp.nstr(abs(lhs-rhs)/abs(rhs), 3))
sq = [mp.mpc(1), mp.mpc(0, 1), mp.mpc(-1), mp.mpc(0, -1)]
tet = [mp.mpc(0)]+[mp.exp(2j*mp.pi*k/3) for k in range(3)]
cf = lambda p: mp.taylor(lambda t: mp.fprod([t-pi for pi in p]), 0, 4)
Ival = lambda p: 12*cf(p)[4]*cf(p)[0]-3*cf(p)[3]*cf(p)[1]+cf(p)[2]**2
check("I vanishes when one mark sits at the centre of an equilateral triangle, and is -12 at the square",
      abs(Ival(tet)) < mp.mpf(10)**(-20) and abs(Ival(sq)+12) < mp.mpf(10)**(-20))

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
