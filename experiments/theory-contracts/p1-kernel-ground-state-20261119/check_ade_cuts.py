"""Native phase exp(i pi eta/2) of the primitive action on every ADE cut S^3/Gamma, Gamma in SU(2).

eta(S^3/Gamma) = -(1/|Gamma|) sum_{g != 1} 1/(2 sin^2(phi_g/2)), g with eigenvalues exp(+-i phi_g)
(verified against the direct spectral computation for Z_p in check_eta_phase_and_period.py).
"""
from fractions import Fraction
from mpmath import mpf, pi, sin, mp
from common import Checks

c = Checks("ade_cuts")
w = lambda num, den: 1 / (2 * sin(pi * mpf(num) / den / 2) ** 2)       # phi = pi*num/den

def cyc(p):   return [(1, (2 * j, p)) for j in range(1, p)]
def bdih(n):  return [(1, (j, n)) for j in range(1, 2 * n)] + [(2 * n, (1, 2))]          # order 4n
TWO_T = [(1, (1, 1)), (6, (1, 2)), (8, (1, 3)), (8, (2, 3))]                                  # order 24
TWO_O = [(1, (1, 1)), (18, (1, 2)), (6, (1, 4)), (6, (3, 4)), (8, (1, 3)), (8, (2, 3))]       # order 48
TWO_I = [(1, (1, 1)), (30, (1, 2)), (20, (1, 3)), (20, (2, 3)), (12, (1, 5)), (12, (3, 5)), (12, (2, 5)), (12, (4, 5))]  # 120

def eta(classes):
    n = 1 + sum(m for m, _ in classes)
    val = sum(m * w(*ph) for m, ph in classes) / n
    return Fraction(str(val)).limit_denominator(10 ** 6), n

def order(e):                          # order of exp(i pi e / 2) = exp(2 pi i e/4)
    return (e / 4).denominator

c.ok("class counts: |2T| = 24, |2O| = 48, |2I| = 120, |Dic_n| = 4n",
     eta(TWO_T)[1] == 24 and eta(TWO_O)[1] == 48 and eta(TWO_I)[1] == 120 and eta(bdih(7))[1] == 28)
c.ok("cyclic check: eta(Z_p) = (p^2 - 1)/(6p) for p = 2..40",
     all(eta(cyc(p))[0] == Fraction(p * p - 1, 6 * p) for p in range(2, 41)))
E = {"E6": eta(TWO_T)[0], "E7": eta(TWO_O)[0], "E8": eta(TWO_I)[0]}
c.ok("E cuts: eta = 167/144 (E6), 383/288 (E7), 1079/720 (E8, Poincare sphere)",
     E == {"E6": Fraction(167, 144), "E7": Fraction(383, 288), "E8": Fraction(1079, 720)}, E)
oE = {k: order(v) for k, v in E.items()}
c.ok("E cuts give phase orders 576, 1152, 2880: no clock", oE == {"E6": 576, "E7": 1152, "E8": 2880}, oE)
clocks = {4, 5, 8, 12, 15, 20, 24, 30}
A = {p: order(eta(cyc(p))[0]) for p in range(2, 121)}
D = {n: order(eta(bdih(n))[0]) for n in range(2, 61)}
hitsA = {p: o for p, o in A.items() if o in clocks}; hitsD = {n: o for n, o in D.items() if o in clocks}
c.ok("among A cuts (Z_p, p <= 120) only Z_5 = A4 gives a clock order, namely 5; no D cut (Dic_n, n <= 60) does",
     hitsA == {5: 5} and hitsD == {}, f"A: {hitsA}, D: {hitsD}")
# general statement for A: order of exp(2 pi i (p^2-1)/(24 p)) is 24p/gcd(p^2-1, 24p) -> equals p iff 24 | p^2-1
from math import gcd
c.ok("A cut phase order is p exactly when 24 divides p^2 - 1 (p coprime to 6); p = 5 is the smallest, next 7, 11, 13",
     all((A[p] == p) == ((p * p - 1) % 24 == 0) for p in A) and [p for p in A if A[p] == p][:4] == [5, 7, 11, 13])
# twisted sectors of the Z_5 cut: eta_s = -[(p^2-1)/6 - s(p-s)]/p
p = 5
tw = {s: Fraction(-((p * p - 1) // 6 - s * (p - s)), p) for s in range(p)}
num = {s: -sum(mp.cos(2 * pi * j * s / p) * w(2 * j, p) for j in range(1, p)) / p for s in range(p)}
c.ok("Z_5 twisted sectors: eta_s = -4/5, 0, 2/5, 2/5, 0; relative phases are tenth roots of unity",
     all(abs(num[s] - mpf(tw[s].numerator) / tw[s].denominator) < mpf(10) ** -25 for s in tw), {s: str(v) for s, v in tw.items()})
c.done({"eta_E": {k: str(v) for k, v in E.items()}, "order_E": oE, "A_orders": A, "D_orders": D})
