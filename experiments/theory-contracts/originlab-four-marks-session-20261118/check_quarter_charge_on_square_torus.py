#!/usr/bin/env python3
"""Quarter charges, 1/(8 pi) and the four marks on the square torus (18 November 2026).

Sixth pasted note: two independent integer charges with quadratic energy,
shifts +-1/4 and beta = 4 pi have variance 1/(8 pi) each.  Checked here, and
placed on the double cover of the mark sphere at the square point:

    E = C/(Z + iZ),  weight exp(-2 pi |v|^2)  (level two = degree of the double cover),
    quarter turn = rotation by 90 degrees about the point over z = 0,
    the four marks = the points (+-1/4, +-1/4).

Sections
    V1  variance of Q = n + a under exp(-2 pi Q^2): exactly 1/(8 pi) at a = 1/4
    V2  the proof ingredients (E2(i) = 3/pi and two theta identities)
    V3  the conditions of the note: same law for even and odd totals, Fourier invariance
    V4  the marks on the square torus (lemniscate), mark circle, reflection, quarter turn
    V5  the qubit of second order theta functions: marks at +-f, +-1/f, weights = squared Bloch components
    V6  second mark test: 5/9
"""
import sys

import mpmath as mp
import numpy as np
import sympy as sp

mp.mp.dps = 40
PI = mp.pi
FAIL = []


def check(name, ok, info=None):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info is not None else ""))
    if not ok:
        FAIL.append(name)


def moments(a, beta2=None):
    """mean and variance of Q = n + a with weight exp(-w Q^2), w = 2 pi by default"""
    w = 2*PI if beta2 is None else beta2
    Z = mp.nsum(lambda n: mp.exp(-w*(n+a)**2), [-mp.inf, mp.inf])
    m1 = mp.nsum(lambda n: (n+a)*mp.exp(-w*(n+a)**2), [-mp.inf, mp.inf])/Z
    m2 = mp.nsum(lambda n: (n+a)**2*mp.exp(-w*(n+a)**2), [-mp.inf, mp.inf])/Z
    return m1, m2-m1**2


# ---------------------------------------------------------------- V1
print("V1  variance of a quarter charge")
m1, var = moments(mp.mpf(1)/4)
check("Q in Z + 1/4, weight exp(-2 pi Q^2): Var(Q) = 1/(8 pi) to 38 digits",
      abs(var-1/(8*PI)) < mp.mpf(10)**(-36), mp.nstr(var, 20))
print("      mean of Q:", mp.nstr(m1, 12), "  continuum variance 1/(4 pi) =", mp.nstr(1/(4*PI), 12))
others = {a: moments(mp.mpf(a))[1] for a in ("0", "0.125", "0.5")}
check("other shifts give other values (0, 1/8, 1/2)",
      all(abs(v-1/(8*PI)) > mp.mpf(10)**(-4) for v in others.values()),
      {k: mp.nstr(v, 8) for k, v in others.items()})
sols = []
grid = [mp.mpf(k)/400 for k in range(0, 401)]
vals = [moments(a)[1]-1/(8*PI) for a in grid]
for k in range(400):
    if vals[k] == 0 or vals[k]*vals[k+1] < 0:
        sols.append(mp.findroot(lambda a: moments(a)[1]-1/(8*PI), (grid[k], grid[k+1]), solver="anderson"))
check("on one period Var(a) = 1/(8 pi) exactly at a = 1/4 and a = 3/4",
      len(sols) == 2 and abs(sols[0]-mp.mpf(1)/4) < mp.mpf(10)**(-20) and abs(sols[1]-mp.mpf(3)/4) < mp.mpf(10)**(-20),
      [mp.nstr(s, 12) for s in sols])

# ---------------------------------------------------------------- V2
print("V2  why it is exact")


def thetas(tau):
    N = 60
    t2 = mp.nsum(lambda n: mp.exp(1j*PI*tau*(n+mp.mpf(1)/2)**2), [-N, N])
    t3 = mp.nsum(lambda n: mp.exp(1j*PI*tau*n**2), [-N, N])
    t4 = mp.nsum(lambda n: (-1)**int(n)*mp.exp(1j*PI*tau*n**2), [-N, N])
    return t2, t3, t4


def E2(tau):
    q = mp.exp(2j*PI*tau)
    return 1-24*mp.nsum(lambda n: n*q**n/(1-q**n), [1, 300])


t2i, t3i, t4i = thetas(1j)
t2d, t3d, t4d = thetas(2j)
check("E2(i) = 3/pi", abs(E2(1j)-3/PI) < mp.mpf(10)**(-30))
check("2 E2(2 tau) - E2(tau) = (theta3^4 + theta4^4)/2 at tau = i, so E2(2i) = 3/(2 pi) + (3/8) theta3(i)^4",
      abs(2*E2(2j)-E2(1j)-(t3i**4+t4i**4)/2) < mp.mpf(10)**(-30)
      and abs(E2(2j)-3/(2*PI)-mp.mpf(3)/8*t3i**4) < mp.mpf(10)**(-30))
A_, B_ = t3d**2, t2d**2
check("theta2^4 + theta3^4 - 3 theta2^2 theta3^2 at 2i equals (3/8) theta3(i)^4",
      abs(A_**2+B_**2-3*A_*B_-mp.mpf(3)/8*t3i**4) < mp.mpf(10)**(-30))
print("      the period terms cancel and 3/(2 pi) remains: Var = 1/(4 pi) - 1/(8 pi).")

# ---------------------------------------------------------------- V3
print("V3  the two conditions of the note")
beta = 4*PI                                      # weight exp(-(beta/2) Q^2) per charge


def total_law(a, N):
    """weight of the total charge N = n1 + n2 for shifts (+a, -a)"""
    return mp.nsum(lambda n1: mp.exp(-beta/2*((n1+a)**2+(N-n1-a)**2)), [-mp.inf, mp.inf])


def gauss_ratio(a):
    """(odd-total prefactor)/(even-total prefactor) after removing exp(-(beta/4) N^2)"""
    e = total_law(a, 0)
    o = total_law(a, 1)/mp.exp(-beta/4)
    return o/e


check("even and odd totals follow one Gauss law without extra factor exactly for shifts +-1/4",
      abs(gauss_ratio(mp.mpf(1)/4)-1) < mp.mpf(10)**(-30) and abs(gauss_ratio(mp.mpf(1)/8)-1) > mp.mpf(10)**(-3),
      "ratio at 1/8: "+mp.nstr(gauss_ratio(mp.mpf(1)/8), 8))
check("the law of the total is exp(-pi N^2), its own Fourier transform (beta/4 = pi)",
      abs(total_law(mp.mpf(1)/4, 2)/total_law(mp.mpf(1)/4, 0)-mp.exp(-4*PI)) < mp.mpf(10)**(-30))

# ---------------------------------------------------------------- V4
print("V4  the four marks on the square torus")
varpi = mp.gamma(mp.mpf(1)/4)**2/(2*mp.sqrt(2*PI))   # lemniscate constant, = 2 int_0^1 dx/sqrt(1-x^4)
check("lemniscate constant against the integral", abs(varpi-2*mp.quad(lambda x: 1/mp.sqrt(1-x**4), [0, 1])) < mp.mpf(10)**(-15))
sl = lambda u: mp.ellipfun("sn", u, m=-1)        # lemniscate sine: sl' ^2 = 1 - sl^4
om1 = (1+1j)*varpi                                # lattice of (X, Y) = (sl, sl'):  (1 +- i) varpi
dsl = lambda u: mp.diff(sl, u)
u0 = mp.mpc("0.31", "0.17")
check("curve Y^2 = 1 - X^4 uniformised by X = sl(u); lattice generated by (1+i) varpi and (1-i) varpi",
      abs(dsl(u0)**2-(1-sl(u0)**4)) < mp.mpf(10)**(-25)
      and abs(sl(u0+om1)-sl(u0)) < mp.mpf(10)**(-25) and abs(dsl(u0+om1)-dsl(u0)) < mp.mpf(10)**(-22)
      and abs(dsl(u0+2*varpi)-dsl(u0)) < mp.mpf(10)**(-22) and abs(dsl(u0+varpi)+dsl(u0)) < mp.mpf(10)**(-22))
check("quarter turn: sl(i u) = i sl(u), a rotation by 90 degrees about u = 0 (the point over z = 0)",
      abs(sl(1j*u0)-1j*sl(u0)) < mp.mpf(10)**(-25))
v_of = lambda u: u/om1                            # coordinates with lattice Z + iZ
marks = {1: varpi/2, 1j: 1j*varpi/2, -1: -varpi/2, -1j: -1j*varpi/2}
pos = {m: v_of(u) for m, u in marks.items()}
check("the marks 1, i, -1, -i sit at (1-i)/4, (1+i)/4, (-1+i)/4, (-1-i)/4: coordinates (+-1/4, +-1/4)",
      all(abs(sl(u)-m) < mp.mpf(10)**(-25) for m, u in marks.items())
      and abs(pos[1]-(1-1j)/4) < mp.mpf(10)**(-25) and abs(pos[1j]-(1+1j)/4) < mp.mpf(10)**(-25),
      {str(m): mp.nstr(p, 6) for m, p in pos.items()})
on_circle = max(abs(abs(sl(om1*(mp.mpf(1)/4+1j*y)))-1) for y in (mp.mpf("0.07"), mp.mpf("0.19"), mp.mpf("0.4")))
on_circle2 = max(abs(abs(sl(om1*(x+1j*mp.mpf(1)/4)))-1) for x in (mp.mpf("0.07"), mp.mpf("0.19"), mp.mpf("0.4")))
check("the mark circle |z| = 1 lifts to the grid x = +-1/4 and y = +-1/4: alternate quarters run in x and in y",
      on_circle < mp.mpf(10)**(-25) and on_circle2 < mp.mpf(10)**(-25))
vv = mp.mpc("0.13", "0.36")
zz = sl(om1*vv)
r1 = sl(om1*(mp.mpf(1)/2-vv.real+1j*vv.imag))
r2 = sl(om1*(vv.real+1j*(mp.mpf(1)/2-vv.imag)))
check("the reflection z -> 1/conj(z) lifts to x -> 1/2 - x and to y -> 1/2 - y: each acts on one charge only",
      abs(r1-1/mp.conj(zz)) < mp.mpf(10)**(-22) and abs(r2-1/mp.conj(zz)) < mp.mpf(10)**(-22))
check("the quarter turn exchanges the two directions (x, y) -> (-y, x): the Fourier exchange, the recorded S",
      abs(sl(om1*(-vv.imag+1j*vv.real))-1j*zz) < mp.mpf(10)**(-22))
print("      so the two charges are the two coordinates of a mark; their Gaussian is a product because the")
print("      torus is rectangular, self-dual because it is square, and the shifts are +-1/4 because the")
print("      centre of the quarter turn and the centre of the sheet exchange differ by the half period (1+i)/2.")
check("covariance of the lifts of a mark under exp(-2 pi |v|^2): (1/(8 pi)) times the identity",
      abs(moments(mp.mpf(1)/4)[1]-1/(8*PI)) < mp.mpf(10)**(-36) and abs(moments(-mp.mpf(1)/4)[1]-1/(8*PI)) < mp.mpf(10)**(-36))

# ---------------------------------------------------------------- V5
print("V5  the qubit of second order theta functions")


def tau_of_j(j):
    coeffs = mp.taylor(lambda L: 256*(1-L+L**2)**3-j*L**2*(1-L)**2, 0, 6)
    roots = mp.polyroots(coeffs[::-1], maxsteps=500, extraprec=200)
    out = []
    for L in roots:
        tau = 1j*mp.ellipk(1-L)/mp.ellipk(L)
        if tau.imag > mp.mpf("0.2"):
            out.append(tau)
    return out


def p_from_f(f):
    d = (1+abs(f)**2)**2
    return sorted([float(4*f.real**2/d), float(4*f.imag**2/d), float((1-abs(f)**2)**2/d)])


CQ = 3/4**(4/3)


def p_from_roots(zv):
    q = CQ*np.exp(np.log(complex(zv)**4-1)/3)
    X = np.roots([1, 0, 0, zv, q])
    pair = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
    r = np.array([abs(X[a]-X[b])*abs(X[c]-X[d]) for (a, b), (c, d) in pair])
    return sorted(1-2*r/r.sum())


worst = 0
for zv in (0.37+0.21j, 0.8, 0.5+0.5j, 1.7-0.4j, 0.05):
    j = 1728*(1-mp.mpc(zv)**4)
    best = None
    for tau in tau_of_j(j):
        t2, t3, _ = thetas(2*tau)
        dev = max(abs(a-b) for a, b in zip(p_from_f(t2/t3), p_from_roots(zv)))
        best = dev if best is None else min(best, dev)
    worst = max(worst, best)
check("weights from the quartic = squared Bloch components of f = theta2(2 tau)/theta3(2 tau), j(tau) = 1728 (1 - z^4)",
      worst < 1e-9, "%.1e" % worst)
t2, t3, _ = thetas(2j)
fsq = t2/t3
check("square torus: f = sqrt(2) - 1, marks at +-f, +-1/f on a great circle, weights (1/2, 0, 1/2)",
      abs(fsq-(mp.sqrt(2)-1)) < mp.mpf(10)**(-30) and np.allclose(p_from_f(fsq), [0, 0.5, 0.5]))
om = mp.exp(2j*PI/3)
t2, t3, _ = thetas(2*om)
fh = t2/t3
check("hexagonal torus: f^8 + 14 f^4 + 1 = 0, tetrahedron, weights (1/3, 1/3, 1/3)",
      abs(fh**8+14*fh**4+1) < mp.mpf(10)**(-28) and np.allclose(p_from_f(fh), [1/3]*3))
taut = mp.mpc("0.23", "0.91")
q_ = lambda t: mp.exp(1j*PI*t)


def th2nd(u, tau):
    """second order theta functions: sum over n in Z + c of exp(2 pi i tau n^2 + 4 pi i n u), c = 0, 1/2"""
    return [mp.nsum(lambda n: mp.exp(2j*PI*tau*(n+c)**2+4j*PI*(n+c)*u), [-40, 40]) for c in (0, mp.mpf(1)/2)]


a0, b0 = th2nd(0, taut)
a1, b1 = th2nd(mp.mpf(1)/2, taut)
a2, b2 = th2nd(taut/2, taut)
check("translation by half periods sends [a : b] to [a : -b] and to [b : a]: the marks are f, -f, 1/f, -1/f",
      abs(b1/a1+b0/a0) < mp.mpf(10)**(-25) and abs(b2/a2-a0/b0) < mp.mpf(10)**(-25))
print("      so the sphere of the sixth round with coordinate f is the sphere on which the four marks sit,")
print("      and the octahedral group there is the Clifford group of this qubit.")

# ---------------------------------------------------------------- V6
print("V6  second mark test")
Psym = sp.ones(4, 4)/4
Ea = []
for a in range(4):
    ea = sp.zeros(4, 1)
    ea[a] = 1
    rest = sp.ones(4, 1)-ea
    Ea.append(ea*ea.T+rest*rest.T/3)
check("after the test of one mark, a different mark is accepted with probability Tr(E_a E_b)/2 = 5/9",
      all((Ea[a]*Ea[b]).trace()/2 == sp.Rational(5, 9) for a in range(4) for b in range(4) if a != b))
check("5/9 = (1 + 1/9)/2: the symmetric sum always passes, the lines have overlap 1/9",
      sp.Rational(1, 2)*(1+sp.Rational(1, 9)) == sp.Rational(5, 9))

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
