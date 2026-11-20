"""Native quarter phase (i pi/2) Delta eta of the primitive action on lens-space cuts, and the period reading of beta = 2 pi.

Dirac on S^3 = SU(2): eigenvalue +(k+3/2) on V_k (x) V_{k+1}, -(k+3/2) on V_{k+1} (x) V_k  (dims (k+1)(k+2)).
Z_p = <diag(zeta, zeta)> acts through U(1) in the second SU(2).  Equivariant eta by Abel summation.
"""
from mpmath import mpf, mpc, pi, sin, expj, exp, nsum, inf
from common import Checks, C3

c = Checks("eta_phase_and_period")
tol = mpf(10) ** -30

def eta_g(phi):                        # Abel-regularised sum_k [(k+1) chi_{k+2} - (k+2) chi_{k+1}](phi)
    w = expj(phi); S0 = 1 / (1 - w); S1 = w / (1 - w) ** 2
    return ((expj(2 * phi) * (S1 + S0) - expj(phi) * (S1 + 2 * S0)).imag) / sin(phi)

def eta_lens(p, eps=1):                # eps = -1: the second spin structure (p even)
    return sum(eps ** j * eta_g(2 * pi * j / p) for j in range(1, p)) / p

# 1 cross-check of the Abel closed form against a damped direct sum
phi = 2 * pi / 5; x = mpf("0.9999")
direct = nsum(lambda k: x ** k * ((k + 1) * sin((k + 2) * phi) - (k + 2) * sin((k + 1) * phi)) / sin(phi), [0, inf])
c.ok("Abel closed form of the equivariant eta sum matches the damped direct sum (x = 0.9999)",
     abs(direct - eta_g(phi)) < mpf(10) ** -3, f"{direct} vs {eta_g(phi)}")
# 2 closed form |eta| = (p^2 - 1)/(6 p)
vals = {p: eta_lens(p) for p in range(2, 13)}
c.ok("Dirac eta of S^3/Z_p (first spin structure): |eta| = (p^2 - 1)/(6p); RP^3 gives 1/4",
     all(abs(abs(v) - mpf(p * p - 1) / (6 * p)) < tol for p, v in vals.items()) and abs(abs(vals[2]) - mpf(1) / 4) < tol,
     {p: str(abs(v))[:8] for p, v in vals.items()})
# 3 Z_4 cut (boundary of the C^2/Z_4 resolution): 5/8 and 3/8 for the two spin structures
e4 = (abs(eta_lens(4, 1)), abs(eta_lens(4, -1)))
c.ok("Z_4 cut: |eta| = 5/8 or 3/8 (two spin structures)", abs(e4[0] - mpf(5) / 8) < tol and abs(e4[1] - mpf(3) / 8) < tol, e4)
# 4 native quarter phase exp(i pi eta / 2): order as root of unity
def order(eta_abs, maxn=400):
    for n in range(1, maxn):
        if abs(expj(pi * eta_abs / 2 * n) - 1) < mpf(10) ** -25: return n
orders = {"round S^2": 1, "RP^3": order(mpf(1) / 4), "Z_4 (5/8)": order(mpf(5) / 8), "Z_4 (3/8)": order(mpf(3) / 8),
          "Z_3": order(mpf(4) / 9), "Z_5": order(mpf(4) / 5)}
c.ok("quarter phase exp(i pi eta/2): trivial on round S^2, order 16 on RP^3, order 32 on the Z_4 cut, order 5 on the Z_5 cut",
     orders["RP^3"] == 16 and orders["Z_4 (5/8)"] == 32 and orders["Z_4 (3/8)"] == 32 and orders["Z_5"] == 5, orders)
# 5 no lens cut gives integer eta: the phase term never disappears off the round sphere
c.ok("(p^2 - 1)/(6p) is never an integer (p = 2..200): every Z_p cut carries a non-trivial native phase",
     all((p * p - 1) % (6 * p) != 0 for p in range(2, 201)))
# 6 which cuts give a phase whose order is one of the E8 clocks 15, 20, 24, 30 (or 4, 5, 8, 12)?
hits = {}
for p in range(2, 61):
    n = order(mpf(p * p - 1) / (6 * p), 2000)
    if n in (4, 5, 8, 12, 15, 20, 24, 30): hits[p] = n
c.ok("clock orders among Z_p cuts (p <= 60): recorded, Z_5 gives the five-clock, Z_4 gives none",
     hits.get(5) == 5 and 4 not in hits, hits)

# --- period reading of beta = 2 pi ---------------------------------------------------------------------
# 7 integer spectrum k/R on a sphere of radius R: real-time period 2 pi R; self-dual thermal circle beta = 2 pi R
for R in (mpf(1), mpf("0.37"), mpf(12)):
    beta = 2 * pi * R
    E_osc = (nsum(lambda n: n / (exp(beta * n / R) - 1), [1, inf]) - mpf(1) / 24) / R
    if abs(E_osc * R + C3) > tol or abs(R / (4 * beta) - C3) > tol: break
c.ok("the radius drops out: with beta = own period 2 pi R, oscillator energy x R = -c3 and zero-mode T R/4 = +c3 for every R",
     abs(E_osc * R + C3) < tol)
# 8 I = 4 as the four real quadratures of C^2: ground level lambda_1 = 1 at h = 1 (hbar = 1/2pi) shared by 4 quadratures
c.ok("c3 = lambda_1 / (I beta) with lambda_1 = 1 (chi_seed), I = 4 = real dimension of C^2, beta = 2 pi = own period",
     abs(mpf(1) / (4 * 2 * pi) - C3) < tol)
# 9 off the self-dual point the two pieces do not cancel: beta = pi and 4 pi
res = {b: (nsum(lambda n: n / (exp(b * n) - 1), [1, inf]) - mpf(1) / 24) + 1 / (4 * b) for b in (pi, 2 * pi, 4 * pi)}
c.ok("the balance oscillators + zero mode = 0 holds only at beta = 2 pi (tested beta = pi, 2 pi, 4 pi)",
     abs(res[2 * pi]) < tol and abs(res[pi]) > mpf("0.01") and abs(res[4 * pi]) > mpf("0.01"), {str(k)[:6]: str(v)[:12] for k, v in res.items()})
c.done({"eta": {p: str(v) for p, v in vals.items()}, "phase_orders": orders, "clock_hits": hits})
