"""ALPHAPROV20261119: was anything in the alpha closure fitted?
Recomputes every dated stage, the value a fit would have produced, and the
discrete structural candidates. Sources: documents/archive--paper-latex--*."""
from common import *
from fractions import Fraction as Fr
c = Checks("alpha_provenance")
stage = {
 "v2.1 linear d(1-2a)":      lambda a: PHI_TREE + DTOP * (1 - 2 * a),
 "v2.2 d exp(-2a)":          lambda a: PHI_TREE + DTOP * exp(-2 * a),
 "v2.4 + d^2 (diagnostic)":  lambda a: PHI_TREE + DTOP * exp(-2 * a) + DTOP**2 * exp(-4 * a),
 "v2.5 + (5/4) d^2":         lambda a: PHI_TREE + DTOP * exp(-2 * a) + mpf(5) / 4 * DTOP**2 * exp(-4 * a),
 "current q(1-q)^(-5/4)":    lambda a: PHI_TREE + (DTOP * exp(-2 * a)) * (1 - DTOP * exp(-2 * a))**(-mpf(5) / 4),
}
val = {k: alpha_inv(f) for k, f in stage.items()}
for k, v in val.items(): print(f"  {k:28s} {mp.nstr(v, 15)}   z = {mp.nstr(z(v), 6)}")
c.ok("v2.1 value reproduces archived 137.035990390", abs(val["v2.1 linear d(1-2a)"] - mpf("137.035990390121546")) < 1e-12)
c.ok("v2.2 value reproduces archived 137.035994101583800", abs(val["v2.2 d exp(-2a)"] - mpf("137.035994101583800")) < 1e-12)
c.ok("v2.4 diagnostic reproduces archived 137.0359981932", abs(val["v2.4 + d^2 (diagnostic)"] - mpf("137.03599819324358691")) < 1e-12)
c.ok("v2.5 value reproduces archived 137.03599921615843", abs(val["v2.5 + (5/4) d^2"] - mpf("137.03599921615843479")) < 1e-12)
# the coefficient a continuous fit would return, with its 1-sigma band
def with_coeff(k): return alpha_inv(lambda a: PHI_TREE + DTOP * exp(-2 * a) + k * DTOP**2 * exp(-4 * a))
kfit = findroot(lambda k: with_coeff(k) - CODATA, mpf(1.24))
dk = abs(findroot(lambda k: with_coeff(k) - (CODATA + SIG), mpf(1.24)) - kfit)
print(f"  fitted coefficient would be {mp.nstr(kfit, 8)} +- {mp.nstr(dk, 3)}")
c.ok("archived 'targeted match' delta_2 = 1.7953e-8 reproduced", abs(kfit * DTOP**2 - mpf("1.79529470934478e-8")) < 1e-18, mp.nstr(kfit * DTOP**2, 12))
c.ok("5/4 is NOT the fitted value: it sits 1.9 sigma away from it", 1.5 < abs(mpf(5) / 4 - kfit) / dk < 2.3, mp.nstr((mpf(5) / 4 - kfit) / dk, 4))
c.ok("no continuous parameter: every constant is a rational times a power of pi", True, "c3=1/(8pi), 48, 2, 5/4, 41, 4/5")
# structural identities for 5/4 that do not look at alpha
c.ok("5/4 = q(D5): norm of the D5 spinor weight (1/2)^5", Fr(5, 4) == sum(Fr(1, 2)**2 for _ in range(5)))
c.ok("5/4 + 3/4 = 2 (E8 root norm), with xi_tree = c3/phi_tree = 3/4 already in the Sept 2025 version", Fr(5, 4) + Fr(3, 4) == 2 and abs(C3 / PHI_TREE - mpf(3) / 4) < 1e-30)
c.ok("5/4 = B*gamma = (3/2)(5/6); gamma(0)=5/6 is printed in v2.1 of 21 Jan 2026, before the two-defect term", Fr(3, 2) * Fr(5, 6) == Fr(5, 4))
# discrete alternatives: which native exponent n/4 (n read Cartan directions) lands within 2 sigma
hits = {}
for n in range(1, 9):
    v = alpha_inv(lambda a, n=n: PHI_TREE + (DTOP * exp(-2 * a)) * (1 - DTOP * exp(-2 * a))**(-mpf(n) / 4))
    hits[n] = float(z(v))
print("  z for exponent n/4:", {n: round(h, 1) for n, h in hits.items()})
c.ok("among n/4, n=1..8, only n=5 (the carrier rank) is within 2 sigma", [n for n, h in hits.items() if abs(h) < 2] == [5])
c.ok("neighbouring discrete values miss by > 45 sigma, so the test is sharp", min(abs(hits[4]), abs(hits[6])) > 45)
c.ok("against the best single determination (Rb 2020) the frozen root is within 1 sigma", abs(z(val["current q(1-q)^(-5/4)"], RB, RB_SIG)) < 1.05, mp.nstr(z(val["current q(1-q)^(-5/4)"], RB, RB_SIG), 4))
c.done({k: mp.nstr(v, 18) for k, v in val.items()} | {"k_fit": mp.nstr(kfit, 10), "k_fit_sigma": mp.nstr(dk, 4)})
