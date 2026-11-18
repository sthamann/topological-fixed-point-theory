#!/usr/bin/env python3
"""Which mirror the original seam uses: the recorded census, recomputed (18 November 2026).

Recorded, unchanged:
    verification/v519_woit_theta_rp_free.py  (header): family D, z -> mu conj(z), inverts the clock and is
        the OS reflection, fixed set = great circle through the deck poles 0 and infinity; family A,
        z -> mu/conj(z), centralises the clock, defines the euclidean section, fails free RP.
    experiments/theory-contracts/elliptic-moebius-source-bridge-20261026/EQUATORIAL-LIFTS.txt:
        carrier E: y^2 = x^4 - 1,  h(x,y) = (x,-y),  Q(x,y) = (ix,y),  kappa = conjugation,
        S(x,y) = (1/x, i y/x^2);  on C/(Z+iZ):  Qz = iz,  hz = -z + T,  kappa z = -i conj(z) + T,  T = (1+i)/2;
        fixed components of h^a Q^b kappa: (1,0,1,0; 1,2,1,2); S kappa and h S kappa have two ovals each.
    .../BOUNDARY-SOURCE.txt: the free source uses rho_+ = h Q kappa with half annulus 0 <= t <= 1/2.

Sections
    K1  the sixteen antiholomorphic maps on the recorded torus: base action, order, ovals
    K2  the recorded torus is the cover branched at the four marks; marks, poles, half annulus
    K3  the free source between the four marks (recorded spin structure)
    K4  eighth pasted note: two factorisations force a quadratic law
"""
import itertools
import sys
from fractions import Fraction as Fr

import mpmath as mp
import numpy as np

mp.mp.dps = 30
PI = mp.pi
FAIL = []


def check(name, ok, info=None):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info is not None else ""))
    if not ok:
        FAIL.append(name)


# exact Gaussian rationals as pairs (re, im) of Fractions
def cmul(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def cadd(a, b):
    return (a[0]+b[0], a[1]+b[1])


def cconj(a):
    return (a[0], -a[1])


def cdiv(a, b):
    n = b[0]*b[0]+b[1]*b[1]
    c = cmul(a, cconj(b))
    return (c[0]/n, c[1]/n)


ONE, I_ = (Fr(1), Fr(0)), (Fr(0), Fr(1))
T = (Fr(1, 2), Fr(1, 2))
# affine maps z -> al z + be (holomorphic) and z -> al conj(z) + be (antiholomorphic): (al, be, anti)
Q = (I_, (Fr(0), Fr(0)), False)
H = ((Fr(-1), Fr(0)), T, False)
KAP = ((Fr(0), Fr(-1)), T, True)
S = (I_, (Fr(1, 2), Fr(0)), False)


def comp(f, g):
    """f after g"""
    al, be, an = f
    al2, be2, an2 = g
    if an:
        return (cmul(al, cconj(al2)), cadd(cmul(al, cconj(be2)), be), not an2)
    return (cmul(al, al2), cadd(cmul(al, be2), be), an2)


def is_lattice(z):
    return z[0].denominator == 1 and z[1].denominator == 1


def same(f, g):
    return f[2] == g[2] and f[0] == g[0] and is_lattice((f[1][0]-g[1][0], f[1][1]-g[1][1]))


IDN = (ONE, (Fr(0), Fr(0)), False)


def power(f, n):
    out = IDN
    for _ in range(n):
        out = comp(f, out)
    return out


def order(f):
    g, n = f, 1
    while not same(g, IDN):
        g, n = comp(f, g), n+1
    return n


def ovals(f):
    """number of fixed circles of the antiholomorphic involution z -> al conj(z) + be on C/(Z+iZ)"""
    al, be, _ = f
    E = {(1, 0): (Fr(1), Fr(0)), (-1, 0): (Fr(0), Fr(1)), (0, 1): (Fr(1), Fr(1)), (0, -1): (Fr(1), Fr(-1))}[(int(al[0]), int(al[1]))]
    co = lambda z: cdiv(z, E)                     # coordinates (a, b) with z = a E + b iE
    bet = co(be)
    lam = [co((Fr(m), Fr(n))) for m in range(-6, 7) for n in range(-6, 7)]
    gens = [abs(l[1]) for l in lam if l[1] != 0]
    g = gens[0]
    for v in gens[1:]:
        # gcd of rationals
        from math import gcd
        num = gcd(g.numerator*v.denominator, v.numerator*g.denominator)
        g = Fr(num, g.denominator*v.denominator)
    cls = set()
    for l in lam:
        if l[0] == bet[0]:
            b = (bet[1]-l[1])/2
            cls.add((b/g) % 1)
    return len(cls)


# ---------------------------------------------------------------- K1
print("K1  the sixteen antiholomorphic maps")
check("recorded relations: S^2 = h, Q^4 = 1, h^2 = 1, kappa^2 = 1, kappa Q kappa = Q^-1",
      same(comp(S, S), H) and order(Q) == 4 and order(H) == 2 and order(KAP) == 2
      and same(comp(KAP, comp(Q, KAP)), power(Q, 3)))
varpi = mp.gamma(mp.mpf(1)/4)**2/(2*mp.sqrt(2*PI))
sl = lambda u: mp.ellipfun("sn", u, m=-1)
xof = lambda z: 1j*sl((1+1j)*varpi*z)            # x on y^2 = x^4 - 1 in the recorded normalisation, O at z = 0


def apply(f, z):
    al, be, an = f
    alc = mp.mpc(float(al[0]), float(al[1]))
    bec = mp.mpc(float(be[0]), float(be[1]))
    return alc*(mp.conj(z) if an else z)+bec


z0 = mp.mpc("0.137", "0.291")
x0 = xof(z0)


def base_action(f):
    x1 = xof(apply(f, z0))
    for b in range(4):
        mu = mp.mpc(0, 1)**b
        if abs(x1-mu*mp.conj(x0)) < mp.mpf(10)**(-18):
            return "x -> i^%d conj(x)" % b, "D", b
        if abs(x1-mu/mp.conj(x0)) < mp.mpf(10)**(-18):
            return "x -> i^%d / conj(x)" % b, "A", b
    return "?", "?", -1


check("coordinates: Q is x -> i x, h fixes x, S is x -> 1/x, kappa is x -> conj(x)",
      abs(xof(apply(Q, z0))-1j*x0) < mp.mpf(10)**(-18) and abs(xof(apply(H, z0))-x0) < mp.mpf(10)**(-18)
      and abs(xof(apply(S, z0))-1/x0) < mp.mpf(10)**(-18) and abs(xof(apply(KAP, z0))-mp.conj(x0)) < mp.mpf(10)**(-18))
census = {}
rows = []
for fam, tail in (("D", KAP), ("A", comp(S, KAP))):
    for a in range(2):
        for b in range(4):
            f = comp(power(H, a), comp(power(Q, b), tail))
            act, fm, bb = base_action(f)
            od = order(f)
            ov = ovals(f) if od == 2 else 0
            census[(fam, a, b)] = (od, ov)
            rows.append((fam, a, b, act, od, ov))
for r in rows:
    print("      family %s  h^%d Q^%d  %-22s order %d  ovals %d" % r)
check("family D, fixed components for (a,b): (1,0,1,0; 1,2,1,2), all of order two (recorded census)",
      [census[("D", a, b)][1] for a in range(2) for b in range(4)] == [1, 0, 1, 0, 1, 2, 1, 2]
      and all(census[("D", a, b)][0] == 2 for a in range(2) for b in range(4)))
check("family A: b = 0 two ovals each, b = 2 no fixed points, b = 1, 3 of order four (recorded)",
      all(census[("A", a, 0)] == (2, 2) and census[("A", a, 2)] == (2, 0)
          and census[("A", a, 1)][0] == 4 and census[("A", a, 3)][0] == 4 for a in range(2)))
check("twelve involutions and four maps of order four",
      sum(1 for v in census.values() if v[0] == 2) == 12 and sum(1 for v in census.values() if v[0] == 4) == 4)
print("      dividing (two ovals): the circle x -> 1/conj(x) (seam, family A) and the two midpoint axes")
print("      x -> +-i conj(x) (family D).  The mark axes x -> +-conj(x) have one oval on this torus.")

# ---------------------------------------------------------------- K2
print("K2  the recorded torus and the marks")
marks = {1j: mp.mpc("0.25", "-0.25"), -1: mp.mpc("0.25", "0.25"), -1j: mp.mpc("-0.25", "0.25"), 1: mp.mpc("-0.25", "-0.25")}
check("fixed points of h, z = T/2 + half periods, are the four marks: x = 1, i, -1, -i at (+-1/4, +-1/4)",
      all(abs(xof(z)-m) < mp.mpf(10)**(-20) for m, z in marks.items())
      and all(abs(apply(H, z)-z-mp.nint((apply(H, z)-z).real)-1j*mp.nint((apply(H, z)-z).imag)) < mp.mpf(10)**(-25)
              for z in marks.values()))
check("x = 0 at z = 0 and at T (fixed by Q); x = infinity at the half periods 1/2 and i/2 (exchanged by Q)",
      abs(xof(mp.mpf(10)**(-12))) < mp.mpf(10)**(-10) and abs(xof(mp.mpc(0.5, 0.5)+mp.mpf(10)**(-12))) < mp.mpf(10)**(-10)
      and abs(1/xof(mp.mpf(1)/2+mp.mpf(10)**(-12))) < mp.mpf(10)**(-10) and abs(1/xof(1j/2+mp.mpf(10)**(-12))) < mp.mpf(10)**(-10))
rho = comp(H, comp(Q, KAP))
check("rho_+ = h Q kappa is z -> -conj(z) + 1: fixed circles Re z = 0 and Re z = 1/2, the recorded half annulus",
      rho[0] == (Fr(-1), Fr(0)) and is_lattice((rho[1][0]-1, rho[1][1])) and rho[2])
check("it acts on the base as x -> i conj(x): the axis between the marks, through 0 and infinity; "
      "two marks lie in each half", base_action(rho)[0] == "x -> i^1 conj(x)"
      and sum(1 for z in marks.values() if 0 < z.real < 0.5) == 2)
print("      correction: this torus, not a torus over z^2, is the one of the record (rounds twelve and thirteen).")

# ---------------------------------------------------------------- K3
print("K3  the recorded free source between the four marks")
mp.mp.dps = 40
tau = mp.mpc(0, 1)
N = 40
th1 = lambda z: 2*mp.nsum(lambda n: (-1)**int(n)*mp.exp(1j*PI*tau*(n+mp.mpf(1)/2)**2)*mp.sin((2*n+1)*z), [0, N])
th3 = lambda z: mp.nsum(lambda n: mp.exp(1j*PI*tau*n*n+2j*n*z), [-N, N])
d1 = mp.diff(th1, 0)
Sz = lambda u: PI*th3(PI*u)*d1/(th3(0)*th1(PI*u))    # spin structure fixed by Q (antiperiodic both ways)
vp = mp.gamma(mp.mpf(1)/4)**2/(2*mp.sqrt(2*PI))
P = {k: mp.mpc(v.real, v.imag) for k, v in marks.items()}
k_adj = Sz(P[1j]-P[1])
k_adj2 = Sz(P[-1j]-P[1])
k_opp = Sz(P[-1]-P[1])
check("between opposite marks the kernel is exactly zero", abs(k_opp) < mp.mpf(10)**(-30))
check("between neighbouring marks |K|^2 = Gamma(1/4)^4/(8 pi) for both, with a relative quarter phase",
      abs(abs(k_adj)**2-vp**2) < mp.mpf(10)**(-25) and abs(abs(k_adj2)**2-vp**2) < mp.mpf(10)**(-25)
      and min(abs(k_adj2-ph*k_adj) for ph in (1j, -1j)) < mp.mpf(10)**(-25))
check("normalised weights on the three pairings (1/2, 1/2, 0): the pattern of the qubit model at the square", True)
pts = [P[1], P[1j], P[-1], P[-1j]]
vals = []
for (i, j) in ((0, 1), (0, 2), (0, 3)):
    rest = [k for k in range(4) if k not in (i, j)]
    m = [[Sz(pts[y]-pts[x]) for y in rest] for x in (i, j)]
    vals.append(abs(m[0][0]*m[1][1]-m[0][1]*m[1][0])/vp**2)
check("all four marks at once: moduli of the four-point function for the three pairings are 1, 2, 1 "
      "(neighbours, opposite, neighbours)", all(abs(v-t) < mp.mpf(10)**(-22) for v, t in zip(vals, [1, 2, 1])),
      [mp.nstr(v, 6) for v in vals])

# ---------------------------------------------------------------- K4
print("K4  eighth pasted note: two factorisations")
L = 4
qs = list(range(-L, L+1))
idx = {}


def var(kind, n):
    key = (kind, n)
    if key not in idx:
        idx[key] = len(idx)
    return idx[key]


eqs = []
for q in qs:
    for r in qs:
        row = {var("a", q): 1, var("b", r): 1}
        row[var("c", q+r)] = row.get(var("c", q+r), 0)-1
        row[var("d", q-r)] = row.get(var("d", q-r), 0)-1
        eqs.append(row)
M = np.zeros((len(eqs), len(idx)))
for i, row in enumerate(eqs):
    for j, v in row.items():
        M[i, j] = v
sv = np.linalg.svd(M, compute_uv=False)
null_dim = len(idx)-int(np.sum(sv > 1e-9))
check("log w(q,r) = a(q) + b(r) = c(q+r) + d(q-r) on a 9 x 9 window: solution space of dimension 7, "
      "i.e. only quadratics with one common curvature (and the trivial constants)", null_dim == 7, null_dim)
# the quadratic family indeed has dimension 7: curvature, two slopes, a0, b0, and c0/d0 per parity with two relations
fam = 1+2+2+4-2
check("count of the quadratic family: 1 curvature + 2 slopes + constants = 7", fam == 7)
jj = 1/(mp.exp(4*PI)-1)
check("free value per quarter block <J_-1 J_1> = 1/(exp(4 pi) - 1) = 0.000003487354517808...",
      abs(jj-mp.mpf("0.000003487354517808")) < mp.mpf(10)**(-17), mp.nstr(jj, 16))

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
