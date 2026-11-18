#!/usr/bin/env python3
"""The diagonal of the recorded seam covariance, and the recorded mark determinant for moved marks
(18 November 2026).

Recorded, unchanged (both run without failure in this session):
    verification/v484_seam_contact_unit.py   G(d) = -(1/pi) ln|2 sin(d/2)|; at the square
        G_marks = -4 c3 ln2 * circ(0,1,2,1), spec {4, -2, 0, -2}; one insertion per mark with weight 1/4
    verification/v485_contact_diagonal_closed.py   G_reg(0; l) = (1/pi) ln(l / 2 pi), zero at l = 2 pi;
        det(I - u C) = (1 - 4u)(1 + 2u)^2, u = eps c3 ln2
    verification/v512_seam_tau_flag.py (H2.1)   the same matrix for marks {+-1, +-e^{i delta}} with free diagonal c

Sections
    D1  the recorded diagonal: zero on the unit circle
    D2  the recorded integer matrix at the square = base-two logarithms of the four responses
    D3  the recorded mark determinant for moved marks, in closed form
    D4  what it gives at the counterwitness of v512
"""
import sys

import mpmath as mp
import sympy as sp

mp.mp.dps = 30
FAIL = []


def check(name, ok, info=None):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info is not None else ""))
    if not ok:
        FAIL.append(name)


c3 = 1/(8*sp.pi)

# ---------------------------------------------------------------- D1
print("D1  the diagonal")


def g_reg(ell, d=mp.mpf("1e-14")):
    return -(1/mp.pi)*mp.log(2*mp.sin(mp.pi*d/ell))+(1/mp.pi)*mp.log(d)


check("G_reg(0; l) = (1/pi) ln(l / 2 pi) at l = pi, 2 pi, 4 pi (flat short-distance subtraction)",
      all(abs(g_reg(l)-mp.log(l/(2*mp.pi))/mp.pi) < mp.mpf(10)**(-20) for l in (mp.pi, 2*mp.pi, 4*mp.pi)))
check("it vanishes exactly on the circle of circumference 2 pi = 1/(4 c3): the free diagonal c of v512 is zero there",
      abs(g_reg(2*mp.pi)) < mp.mpf(10)**(-20) and sp.simplify(1/(4*c3)-2*sp.pi) == 0)

# ---------------------------------------------------------------- D2
print("D2  the square")
dl = sp.symbols("delta", positive=True)
G = lambda th: -sp.log(2*sp.sin(th/2))/sp.pi


def G4(delta, diag=0):
    a, b, o = G(delta), G(sp.pi-delta), G(sp.pi)
    return sp.Matrix([[diag, a, o, b], [a, diag, b, o], [o, b, diag, a], [b, o, a, diag]])


Gsq = G4(sp.pi/2).applyfunc(sp.simplify)
C = sp.Matrix([[0, 1, 2, 1], [1, 0, 1, 2], [2, 1, 0, 1], [1, 2, 1, 0]])
check("recorded: G_marks = -4 c3 ln2 * circ(0,1,2,1)", (Gsq+4*c3*sp.log(2)*C).applyfunc(sp.simplify) == sp.zeros(4, 4))
check("recorded: spectrum of the integer matrix {4, 0, -2, -2}", C.eigenvals() == {4: 1, 0: 1, -2: 2})
chars = {"symmetric": sp.Matrix([1, 1, 1, 1]), "opposite": sp.Matrix([1, -1, 1, -1]),
         "side": sp.Matrix([1, 1, -1, -1]), "side'": sp.Matrix([1, -1, -1, 1])}
ev = {k: (C*v)[0]/v[0] for k, v in chars.items()}
resp_sq = {"symmetric": 16, "opposite": 1, "side": sp.Rational(1, 4), "side'": sp.Rational(1, 4)}
check("2^(eigenvalue) = (16, 1, 1/4, 1/4): sqrt|Disc|, the opposite pairing and the two side pairings at the square",
      all(sp.Integer(2)**ev[k] == resp_sq[k] for k in chars), ev)
uu = sp.symbols("u")
check("recorded closed form det(I - u C) = (1 - 4u)(1 + 2u)^2",
      sp.expand((sp.eye(4)-uu*C).det()-(1-4*uu)*(1+2*uu)**2) == 0)

# ---------------------------------------------------------------- D3
print("D3  moved marks")
eps = sp.symbols("epsilon")
R = {"symmetric": 16*sp.sin(dl)**2, "opposite": 1/sp.sin(dl)**2,
     "side": sp.tan(dl/2)**2/4, "side'": sp.cot(dl/2)**2/4}
Gd = G4(dl)
ok = True
for dval in (sp.Rational(2, 5), sp.Rational(9, 10), sp.Rational(13, 10)):
    M = Gd.subs(dl, dval)
    for k, v in chars.items():
        lam = (M*v)[0]/v[0]
        ok &= abs(float((sp.exp(-2*sp.pi*lam)-R[k].subs(dl, dval)).evalf())) < 1e-12
check("with the recorded diagonal zero: exp(-2 pi eigenvalue) = 16 sin^2(delta), 1/sin^2(delta), "
      "tan^2(delta/2)/4, cot^2(delta/2)/4, with no free constant", ok)
check("their product is one for every delta (trace zero)", sp.simplify(sp.prod(list(R.values()))-1) == 0)
D_closed = sp.prod([1-eps*c3*sp.log(r) for r in R.values()])
okd = True
for dval in (sp.Rational(2, 5), sp.Rational(9, 10)):
    for ev_ in (sp.Rational(1, 3), sp.Rational(-7, 5)):
        lhs = (sp.eye(4)+ev_/4*Gd.subs(dl, dval)).det()
        okd &= abs(float((lhs-D_closed.subs({dl: dval, eps: ev_})).evalf())) < 1e-12
check("recorded mark determinant, one insertion of weight 1/4 per mark: "
      "det(I + (eps/4) G_marks) = prod_k (1 - eps c3 ln R_k(delta))", okd)
check("at the square it is the recorded (1 - 4u)(1 + 2u)^2 with u = eps c3 ln2",
      sp.simplify(D_closed.subs(dl, sp.pi/2)-((1-4*uu)*(1+2*uu)**2).subs(uu, eps*c3*sp.log(2))) == 0)

# ---------------------------------------------------------------- D4
print("D4  the counterwitness of v512")
dcw = sp.atan(sp.Rational(4, 3))
Rcw = {k: sp.nsimplify(sp.simplify(r.subs(dl, dcw))) for k, r in R.items()}
check("responses at {+-1, +-(3+4i)/5}: 256/25, 25/16, 1/16, 1",
      Rcw == {"symmetric": sp.Rational(256, 25), "opposite": sp.Rational(25, 16), "side": sp.Rational(1, 16), "side'": 1}, Rcw)
Dcw = sp.prod([1-eps*c3*sp.log(r) for r in Rcw.values()])
check("mark determinant there: (1 - eps c3 ln(256/25)) (1 - eps c3 ln(25/16)) (1 + 4 eps c3 ln2)",
      sp.simplify(Dcw-(1-eps*c3*sp.log(sp.Rational(256, 25)))*(1-eps*c3*sp.log(sp.Rational(25, 16)))*(1+4*eps*c3*sp.log(2))) == 0)
second = lambda rr: -sum(sp.log(r)**2 for r in rr.values())/2
check("first order in eps absent for every delta; second-order coefficient of the logarithm: "
      "-12 (ln2)^2 at the square, as recorded",
      sp.simplify(second(resp_sq)+12*sp.log(2)**2) == 0)
val = sp.N(second(Rcw)/sp.log(2)**2, 12)
print("      second-order coefficient at the counterwitness, in units of (ln2)^2:", val)
check("it differs from the square: the recorded calculus sees the moved marks at second order", abs(val+12) > 0.1)

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
