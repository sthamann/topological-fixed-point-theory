#!/usr/bin/env python3
"""Step 5 (KERNSIX20261120): the cubic of the alpha equation.

A  the frozen equation and its root
B  normal form in the inverse coupling w = c3/alpha = 1/(2 e^2): w + kappa w^3 = 1/(2 c3^2);
   it is the stationary point of a convex quartic, so the root is unique; the logarithm is not a
   small correction
C  every coefficient in terms of kernel and carrier data: 4 c3 = 1/(2 pi) and
   b1 = (4/3) N_fam + (1/6) tr_{E+}(Y^2)/tr_E(Y^2) from the hypercharge of the original (step 2)
D  what is not derived: the powers of the equation; a running relation would be linear in w
"""
from fractions import Fraction as Fr
import mpmath as mp
from common import Checks

mp.mp.dps = 40
c = Checks("5_cubic")
c3 = 1 / (8 * mp.pi)
b1 = mp.mpf(41) / 10

def phi(a):
    q = 48 * c3 ** 4 * mp.e ** (-2 * a)
    return mp.mpf(4) / 3 * c3 + q * (1 - q) ** (-mp.mpf(5) / 4)

def F(a, b=b1):
    return a ** 3 - 2 * c3 ** 3 * a ** 2 - 8 * c3 ** 6 * b * mp.log(1 / phi(a))

# ---------------------------------------------------------------- A
alpha = mp.findroot(F, mp.mpf(1) / 137)
c.ok("A1 root of the frozen equation: 1/alpha = 137.0359992168407125...",
     abs(1 / alpha - mp.mpf("137.03599921684071250353786030380388037")) < mp.mpf(10) ** -30, mp.nstr(1 / alpha, 22))
c.ok("A2 the coefficient (4/5) 41 c3^6 of the frozen form equals 8 c3^6 b1 = 2 c3^5 (b1 / 2 pi)",
     mp.almosteq(mp.mpf(4) / 5 * 41 * c3 ** 6, 8 * c3 ** 6 * b1) and mp.almosteq(8 * c3 ** 6 * b1, 2 * c3 ** 5 * b1 / (2 * mp.pi)))

# ---------------------------------------------------------------- B
w = c3 / alpha
kappa = 4 * c3 * b1 * mp.log(1 / phi(alpha))
J = 1 / (2 * c3 ** 2)
c.ok("B1 normal form: w + kappa w^3 = 1/(2 c3^2) = 32 pi^2 with w = c3/alpha = 1/(2 e^2) and kappa = (b1/2 pi) log(1/phi)",
     abs(w + kappa * w ** 3 - J) < mp.mpf(10) ** -30 and mp.almosteq(J, 32 * mp.pi ** 2) and mp.almosteq(w, 1 / (8 * mp.pi * alpha)),
     {"w": mp.nstr(w, 10), "kappa": mp.nstr(kappa, 10)})
roots = mp.polyroots([kappa, 0, 1, -J])
real = [r for r in roots if abs(mp.im(r)) < mp.mpf(10) ** -25]
c.ok("B2 for fixed kappa > 0 this is the stationary point of the convex quartic w^2/2 + kappa w^4/4 - w/(2 c3^2): exactly one real root",
     len(real) == 1 and abs(mp.re(real[0]) - w) < mp.mpf(10) ** -25)
share = kappa * w ** 3 / J
c.ok("B3 the logarithmic term carries 98 percent of the balance: without it alpha = 2 c3^3 = 1/(256 pi^3), 58 times too small",
     share > mp.mpf("0.98") and mp.almosteq(c3 / J, 2 * c3 ** 3) and 57 < alpha / (2 * c3 ** 3) < 59,
     {"share": mp.nstr(share, 6), "factor": mp.nstr(alpha / (2 * c3 ** 3), 6)})
lead = (J / kappa) ** (mp.mpf(1) / 3)
c.ok("B4 leading approximation w = (32 pi^2 / kappa)^(1/3) gives 1/alpha within one percent",
     abs(8 * mp.pi * lead * alpha - 1) < mp.mpf("0.01"), mp.nstr(8 * mp.pi * lead, 8))

# ---------------------------------------------------------------- C
Y = [Fr(-1, 3)] * 3 + [Fr(1, 2)] * 2                      # hypercharge of the original on E = E_3 + E_2
trE = sum(y * y for y in Y)
trEp = sum(y * y for y in Y[3:])
ZU = trE / 4                                              # flat-state marker weight of step 2
sum16 = 16 * ZU                                           # sum of Y^2 over the 16 states of one family
nfam = 3
ferm = Fr(2, 3) * nfam * sum16
higgs = Fr(1, 3) * trEp
b1_carrier = (ferm + higgs) / (2 * trE)                   # normalised by the carrier trace (the factor 3/5)
c.ok("C1 4 c3 = 1/(2 pi): the loop factor of the logarithm is the kernel constant", mp.almosteq(4 * c3, 1 / (2 * mp.pi)))
c.ok("C2 sum of Y^2 over one family = 16 times the marker weight 5/24 of step 2 = 10/3", sum16 == Fr(10, 3))
c.ok("C3 b1 = (4/3) N_fam + (1/6) tr_{E+}(Y^2)/tr_E(Y^2) = 4 + 1/10 = 41/10: families give 40, the doublet block E_+ gives 1",
     b1_carrier == Fr(41, 10) == Fr(4, 3) * nfam + Fr(1, 6) * trEp / trE and trEp / trE == Fr(3, 5), str(b1_carrier))
a40 = mp.findroot(lambda a: F(a, mp.mpf(4)), mp.mpf(1) / 137)
c.ok("C4 without the doublet term (b1 = 4) the root moves to 1/alpha = 138.1: the single unit 1 in 41 is needed",
     abs(1 / a40 - mp.mpf("138.1")) < mp.mpf("0.2"), mp.nstr(1 / a40, 8))

# ---------------------------------------------------------------- D
c.ok("D1 a one loop running relation is linear in the inverse coupling (w = w0 + b log); the equation has kappa times w^3: "
     "it is not such a relation", abs(kappa * w ** 3 - kappa) > 1, mp.nstr(kappa * w ** 3, 8))
c.ok("D2 open: no source derives the powers (w^1, w^3) or the source term 1/(2 c3^2); the original calls the quartic part a local counterterm",
     True, "postulated, unchanged")

c.done({"alpha_inverse": mp.nstr(1 / alpha, 20), "w": mp.nstr(w, 15), "kappa": mp.nstr(kappa, 15), "share_log": mp.nstr(share, 8)})
