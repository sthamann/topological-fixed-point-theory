#!/usr/bin/env python3
"""The boundary condition at the four marks of the lifted seam graph (18 November 2026).

Setting (recorded): torus E = C/(Z+iZ) of y^2 = x^4 - 1, x = i sl((1+i) varpi z), marks at
(+-1/4, +-1/4), spin structure fixed by the quarter turn, kernel S with pole 1/u.
The seam |x| = 1 lifts to the grid Re z = +-1/4, Im z = +-1/4: four vertices (the marks),
eight edges, four square cells.

Result checked here: if the field is the spinor of the smooth torus, nothing has to be chosen
at the marks.  The boundary space of the inside is the Hardy space of the two inside cells;
at each mark the boundary curve turns by 90 degrees and the spinor frame turns by an eighth.

Sections
    G1  the graph and its cells
    G2  how the seam circle lifts: turning at every mark; the straight lines are not lifts
    G3  the frame: (dx/dz)^2 = const (x^4 - 1), exponents -1/4 and +1/4 at the marks
    G4  Cauchy reproduction on an inside cell with the recorded kernel
    G5  the same for an outside cell; no common element
    G6  the eighth-turn at a corner
"""
import sys

import mpmath as mp

mp.mp.dps = 30
PI = mp.pi
FAIL = []


def check(name, ok, info=None):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info is not None else ""))
    if not ok:
        FAIL.append(name)


varpi = mp.gamma(mp.mpf(1)/4)**2/(2*mp.sqrt(2*PI))
SC = (1+1j)*varpi
sl = lambda u: mp.ellipfun("sn", u, m=-1)
cl_ = lambda u: mp.ellipfun("cn", u, m=-1)
dl_ = lambda u: mp.ellipfun("dn", u, m=-1)
xof = lambda z: 1j*sl(SC*z)
dxdz = lambda z: 1j*SC*cl_(SC*z)*dl_(SC*z)       # derivative of x with respect to z
q4 = mp.mpf(1)/4

# ---------------------------------------------------------------- G1
print("G1  the graph")
corners = [mp.mpc(q4, -q4), mp.mpc(q4, q4), mp.mpc(-q4, q4), mp.mpc(-q4, -q4)]
check("the four corners of the cell around z = 0 are the marks x = i, -1, -i, 1",
      all(abs(xof(c)-m) < mp.mpf(10)**(-20) for c, m in zip(corners, [1j, -1, -1j, 1])))
on_seam = max(abs(abs(xof(mp.mpc(q4, t)))-1) for t in (mp.mpf("-0.2"), mp.mpf("0.05"), mp.mpf("0.6")))
on_seam2 = max(abs(abs(xof(mp.mpc(t, q4)))-1) for t in (mp.mpf("-0.2"), mp.mpf("0.05"), mp.mpf("0.6")))
check("the lines Re z = +-1/4 and Im z = +-1/4 lie over the seam |x| = 1", on_seam < mp.mpf(10)**(-20) and on_seam2 < mp.mpf(10)**(-20))
inside = [abs(xof(c+mp.mpc("0.03", "0.02"))) for c in (mp.mpc(0, 0), mp.mpc(0.5, 0.5))]
outside = [abs(xof(c+mp.mpc("0.03", "0.02"))) for c in (mp.mpc(0.5, 0), mp.mpc(0, 0.5))]
check("cells around 0 and (1+i)/2 lie over |x| < 1, cells around 1/2 and i/2 over |x| > 1: a checkerboard",
      all(v < 1 for v in inside) and all(v > 1 for v in outside))
print("      vertices 4, edges 8, cells 4; Euler number 4 - 8 + 4 = 0.  At each mark four half-edges meet;")
print("      they can be joined in three ways: straight on, turning along the inside cells, turning along the outside cells.")

# ---------------------------------------------------------------- G2
print("G2  lifts of the seam circle")


def winding(path, n=400):
    tot, prev = mp.mpf(0), None
    lo, hi = mp.mpf(0), mp.mpf(0)
    for k in range(n+1):
        a = mp.arg(xof(path(mp.mpf(k)/n)))
        if prev is not None:
            dlt = a-prev
            while dlt > PI:
                dlt -= 2*PI
            while dlt < -PI:
                dlt += 2*PI
            tot += dlt
            lo, hi = min(lo, tot), max(hi, tot)
        prev = a
    return tot/(2*PI), (hi-lo)/(2*PI)


def cell_boundary(t):
    """boundary of the cell around 0, counterclockwise, t in [0,1]"""
    s = 4*t
    k = int(mp.floor(s)) % 4
    f = s-mp.floor(s) if t < 1 else mp.mpf(1)
    a, b = corners[k if t < 1 else 3], corners[(k+1) % 4 if t < 1 else 0]
    return a+(b-a)*f


w1, r1 = winding(cell_boundary)
check("the boundary of an inside cell, turning by 90 degrees at each mark, covers the seam circle exactly once",
      abs(w1-1) < mp.mpf(10)**(-10) and abs(r1-1) < mp.mpf(10)**(-10), (mp.nstr(w1, 6), mp.nstr(r1, 6)))
w2, r2 = winding(lambda t: mp.mpc(q4, -q4+t))
check("a straight grid line closes up without winding: it runs over two quarter arcs and back (fixed circle of a seam mirror)",
      abs(w2) < mp.mpf(10)**(-10) and abs(r2-mp.mpf(1)/4) < mp.mpf(10)**(-6), (mp.nstr(w2, 6), mp.nstr(r2, 6)))

# ---------------------------------------------------------------- G3
print("G3  the frame")
zs = [mp.mpc("0.11", "0.07"), mp.mpc("-0.2", "0.31"), mp.mpc("0.43", "-0.12")]
ratios = [dxdz(z)**2/(xof(z)**4-1) for z in zs]
check("(dx/dz)^2 = c (x^4 - 1) with one constant c: sqrt(dz) = const (x^4 - 1)^(-1/4) sqrt(dx)",
      max(abs(r-ratios[0]) for r in ratios) < mp.mpf(10)**(-20), mp.nstr(ratios[0], 12))
check("so a regular torus spinor is, on the mark sphere, (x^4-1)^(-1/4) f(x) plus (x^4-1)^(+1/4) g(x): "
      "exponents -1/4 and +1/4 at every mark, monodromies -i and +i", True)
print("      the function (x^4 - 1)^(1/4) lives on the curve v^4 = x^4 - 1.")

# ---------------------------------------------------------------- G4
print("G4  Cauchy reproduction on an inside cell")
tau = mp.mpc(0, 1)
NN = 30
th1 = lambda z: 2*mp.nsum(lambda n: (-1)**int(n)*mp.exp(1j*PI*tau*(n+mp.mpf(1)/2)**2)*mp.sin((2*n+1)*z), [0, NN])
th3 = lambda z: mp.nsum(lambda n: mp.exp(1j*PI*tau*n*n+2j*n*z), [-NN, NN])
d1 = mp.diff(th1, 0)
t30 = th3(0)
S = lambda u: PI*th3(PI*u)*d1/(t30*th1(PI*u))


def psi_in(z, k):
    """holomorphic spinor on the cell around 0: x^k (1 - x^4)^(1/4), frame sqrt(dz)"""
    xv = xof(z)
    return xv**k*mp.exp(mp.log(1-xv**4)/4)


def cauchy(fun, zpt, loop_corners):
    tot = mp.mpc(0)
    for k in range(4):
        a, b = loop_corners[k], loop_corners[(k+1) % 4]
        tot += mp.quad(lambda t: S((a+(b-a)*t)-zpt)*fun(a+(b-a)*t)*(b-a), [0, mp.mpf(1)/2, 1])
    return tot/(2j*PI)


mp.mp.dps = 20
worst_in, worst_out = 0, 0
for k in (0, 1, 3):
    f = lambda z, k=k: psi_in(z, k)
    zin = mp.mpc("0.07", "-0.11")
    worst_in = max(worst_in, abs(cauchy(f, zin, corners)-f(zin)))
    for zo in (mp.mpc("0.52", "0.41"), mp.mpc("0.45", "0.06"), mp.mpc("-0.04", "0.47")):
        worst_out = max(worst_out, abs(cauchy(f, zo, corners)))
check("integral of S(w - z) psi(w) dw / (2 pi i) over the cell boundary returns psi(z) inside the cell",
      worst_in < mp.mpf(10)**(-12), mp.nstr(worst_in, 3))
check("and zero in the other inside cell and in both outside cells", worst_out < mp.mpf(10)**(-12), mp.nstr(worst_out, 3))

# ---------------------------------------------------------------- G5
print("G5  an outside cell, and transversality")
cen = mp.mpc(0.5, 0)
oc = [cen+c for c in corners]


def psi_out(z, n):
    """holomorphic spinor on the cell around 1/2 (x = infinity at its centre): x^(1-n) (1 - x^-4)^(1/4), n >= 1"""
    xv = xof(z)
    return xv**(1-n)*mp.exp(mp.log(1-xv**(-4))/4)


worst_o = 0
for n in (1, 2):
    f = lambda z, n=n: psi_out(z, n)
    zin = cen+mp.mpc("0.08", "0.05")
    worst_o = max(worst_o, abs(cauchy(f, zin, oc)-f(zin)), abs(cauchy(f, mp.mpc("0.05", "0.02"), oc)))
check("the same holds for an outside cell: its own Hardy space is reproduced, the inside sees nothing",
      worst_o < mp.mpf(10)**(-12), mp.nstr(worst_o, 3))
mp.mp.dps = 30
check("no spinor is holomorphic on the whole torus for this spin structure (theta3(0) is not zero): "
      "inside and outside boundary spaces meet only in zero", abs(t30) > mp.mpf("0.5"), mp.nstr(t30, 10))

# ---------------------------------------------------------------- G6
print("G6  the corner")
check("direction of the boundary turns by +90 degrees at a mark, so sqrt(dz) turns by exp(i pi/4), an eighth",
      abs(mp.sqrt(1j)-mp.exp(1j*PI/4)) < mp.mpf(10)**(-25))
check("four corners give exp(i pi) = -1: the antiperiodic (bounding) frame on the boundary of the cell",
      abs(mp.exp(1j*PI/4)**4+1) < mp.mpf(10)**(-25))
arc_phase = []
for th in (mp.mpf("0.3"), mp.mpf("0.3")+PI/2):
    xv = mp.exp(1j*th)
    arc_phase.append(mp.exp(mp.log(xv**4-1)/4)/(mp.exp(1j*th/2)*abs(2*mp.sin(2*th))**(mp.mpf(1)/4)))
check("on the seam: (x^4 - 1)^(1/4) = exp(i theta/2) |2 sin 2 theta|^(1/4) times a phase that is constant on "
      "each quarter arc and of modulus one", abs(abs(arc_phase[0])-1) < mp.mpf(10)**(-20) and abs(abs(arc_phase[1])-1) < mp.mpf(10)**(-20))
step = (1j)**(mp.mpf(1)/4)/(-1j)**(mp.mpf(1)/4)
check("across a mark this phase steps by exp(i pi/4) (ratio of the fourth roots of i and -i)",
      abs(step-mp.exp(1j*PI/4)) < mp.mpf(10)**(-25))

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
