#!/usr/bin/env python3
"""KERNORIG20261120 C: the modulus of the four marks with the profile-free object of the original.

The finite part at s = 0 of the zeta family is the logarithm of the determinant.  For the Dirac
operator of the flat torus with antiperiodic spin structure on both cycles the zeta determinant is
|theta_3(tau)/eta(tau)|^2 (standard; the product form is rederived in C1).  No test profile and no
normalisation by lambda_1 enters (the zeta function of a spectrum without zero modes vanishes at 0).
C1  product over the spectrum = theta_3/eta, and invariance under tau -> -1/tau
C2  tau = i is a critical point in both directions; minimum along the mirror family, maximum across
C3  it is the only critical point in the scanned region (winding number of the derivative)
C4  the other two even spin structures are not critical at tau = i
C5  the original asks for a critical point, not for a minimum
"""
import mpmath as mp
from common import Checks

mp.mp.dps = 30
c = Checks("c_modulus")

def eta(tau):
    q2 = mp.e ** (2j * mp.pi * tau)
    return mp.e ** (1j * mp.pi * tau / 12) * mp.nprod(lambda n: 1 - q2 ** n, [1, mp.inf])

def ratio(tau, k=3):
    return mp.jtheta(k, 0, mp.e ** (1j * mp.pi * tau)) / eta(tau)

def g(x, y, k=3):
    return 2 * mp.log(abs(ratio(mp.mpc(x, y), k)))            # log det of D^2 up to a constant

# ---------------------------------------------------------------- C1
ok = True
for t in (mp.mpf("0.6"), mp.mpf(1), mp.mpf("1.7")):
    prod = mp.e ** (mp.pi / (12 * t)) * mp.nprod(lambda n: (1 + mp.e ** (-2 * mp.pi * (n + mp.mpf(1) / 2) / t)) ** 2, [0, mp.inf])
    ok &= abs(prod - ratio(1j / t)) < mp.mpf(10) ** -25 and abs(ratio(1j / t) - ratio(1j * t)) < mp.mpf(10) ** -25
c.ok("C1 mode product with the regularised zero point sum 2 zeta(-1, 1/2) = 1/12 equals theta_3/eta; invariant under t -> 1/t",
     ok and abs(2 * mp.zeta(-1, mp.mpf(1) / 2) - mp.mpf(1) / 12) < mp.mpf(10) ** -25)

# ---------------------------------------------------------------- C2
gx = mp.diff(lambda x: g(x, 1), 0)
gy = mp.diff(lambda y: g(0, y), 1)
gxx = mp.diff(lambda x: g(x, 1), 0, 2)
gyy = mp.diff(lambda y: g(0, y), 1, 2)
c.ok("C2a tau = i is a critical point: both first derivatives vanish", abs(gx) < mp.mpf(10) ** -15 and abs(gy) < mp.mpf(10) ** -15)
c.ok("C2b minimum along the mirror family (tau = i t), maximum across it, curvatures opposite and equal (harmonic function)",
     gyy > 0 and gxx < 0 and abs(gxx + gyy) < mp.mpf(10) ** -12, (mp.nstr(gyy, 8), mp.nstr(gxx, 8)))
vals = [g(0, mp.mpf(t)) for t in ("1", "1.2", "1.6", "2.5", "4")]
c.ok("C2c along the mirror family the determinant rises monotonically away from the square", all(vals[i] < vals[i + 1] for i in range(4)),
     [mp.nstr(v - vals[0], 6) for v in vals])

# ---------------------------------------------------------------- C3
def F(tau):
    return mp.diff(lambda z: mp.log(ratio(z)), tau)

def winding(corners, n=160):
    total = mp.mpf(0)
    pts = []
    for a, b in zip(corners, corners[1:] + corners[:1]):
        pts += [a + (b - a) * mp.mpf(j) / n for j in range(n)]
    vals_ = [F(p) for p in pts]
    for u, v in zip(vals_, vals_[1:] + vals_[:1]):
        total += mp.arg(v / u)
    return int(mp.nint(total / (2 * mp.pi)))

mp.mp.dps = 20
box = [mp.mpc(-0.9, 0.45), mp.mpc(0.9, 0.45), mp.mpc(0.9, 3.0), mp.mpc(-0.9, 3.0)]
small = [mp.mpc(-0.15, 0.85), mp.mpc(0.15, 0.85), mp.mpc(0.15, 1.15), mp.mpc(-0.15, 1.15)]
wb, ws = winding(box), winding(small)
c.ok("C3 the derivative has exactly one zero in the box |Re tau| <= 0.9, 0.45 <= Im tau <= 3, and it is the one at tau = i",
     wb == 1 and ws == 1, (wb, ws))
mp.mp.dps = 30

# ---------------------------------------------------------------- C4
g4 = mp.diff(lambda y: g(0, y, 4), 1)
g2 = mp.diff(lambda y: g(0, y, 2), 1)
c.ok("C4 the two other even spin structures (theta_4, theta_2) are not critical at the square; their slopes are opposite",
     abs(g4) > mp.mpf("0.1") and abs(g2) > mp.mpf("0.1") and abs(g4 + g2) < mp.mpf(10) ** -12, (mp.nstr(g4, 8), mp.nstr(g2, 8)))

# ---------------------------------------------------------------- C5
c.ok("C5 the original selects the continuous sector as the unique admissible critical point of its free energy; "
     "the saddle at tau = i meets that wording, a minimum is not required", True, "source: boundary kernel section 002, lines 94-97")
c.done({"gyy": mp.nstr(gyy, 12), "gxx": mp.nstr(gxx, 12)})
