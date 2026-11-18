"""F1 - four-mark determinant, Dirac determinants and the missing sector weight.

The chat names no formula.  This script tests the natural free-branch candidate and says so:
  marks        e1,e2,e3,infinity (branch points of the double cover), d_ij = e_i - e_j
  mark side    Vandermonde(e)^2 = prod d_ij^2  (discriminant / 16)
  Dirac side   det D_nu = theta_nu / eta for the three even spin structures of the torus
Two independent routes: marks from the Eisenstein series (roots of 4x^3 - g2 x - g3),
Dirac side from theta constants.  80+ digits, two real and two generic complex moduli.
Then: which part of the E8 character is fixed by the marks, and what fixes the rest?
"""
import json
from mpmath import mp, mpf, mpc, pi, exp, jtheta, polyroots, sqrt, nstr, log, floor

mp.dps = 95
res = {"digits_requested": 80}

def eis(tau, k, N=400):
    q = exp(2j * pi * tau)
    c = {4: 240, 6: -504}[k]
    s = mpf(1)
    for n in range(1, N):
        s += c * n ** (k - 1) * q ** n / (1 - q ** n)
    return s

def eta(tau, N=400):
    q = exp(2j * pi * tau)
    p = exp(2j * pi * tau / 24)
    for n in range(1, N):
        p *= (1 - q ** n)
    return p

def thetas(tau):
    nome = exp(1j * pi * tau)
    return [jtheta(i, 0, nome) for i in (2, 3, 4)]

def chi_dirac(tau, w=mpf(1)):
    """16 Majorana fermions: NS sectors (theta3, theta4) plus w times the spinor sector (theta2)."""
    t2, t3, t4 = thetas(tau)
    return (t3 ** 8 + t4 ** 8 + w * t2 ** 8) / (2 * eta(tau) ** 8)

def digits(a, b):
    d = abs(a - b) / max(abs(a), abs(b))
    return 95 if d == 0 else int(floor(-log(d, 10)))

points = {"real 1": mpc(0, 1.3), "real 2": mpc(0, 0.8), "complex 1": mpc(0.23, 0.91), "complex 2": mpc(-0.37, 1.17)}
rows = {}
for name, tau in points.items():
    t2, t3, t4 = thetas(tau)
    g2 = 4 * pi ** 4 / 3 * eis(tau, 4)
    g3 = 8 * pi ** 6 / 27 * eis(tau, 6)
    e = polyroots([4, 0, -g2, -g3], maxsteps=500, extraprec=400)
    d2 = [(e[i] - e[j]) ** 2 for i in range(3) for j in range(i + 1, 3)]
    vand2 = d2[0] * d2[1] * d2[2]                                   # mark side, from Eisenstein route
    dirac_prod = (t2 / eta(tau)) * (t3 / eta(tau)) * (t4 / eta(tau))  # product of the three Dirac dets
    rhs = pi ** 12 * eta(tau) ** 24 * dirac_prod ** 8
    lam = t2 ** 4 / t3 ** 4
    # cross-ratio of the marks must equal lambda(tau) on one of its six branches
    cr = [(e[a] - e[b]) / (e[c] - e[b]) for a, b, c in ((0, 1, 2), (1, 0, 2), (0, 2, 1), (2, 0, 1), (1, 2, 0), (2, 1, 0))]
    # E8 character from the marks alone: chi^3 = 32 (sum d^2)^3 / prod d^2, cube-root branch free
    chi = chi_dirac(tau)
    mark_expr = mpf(2) ** (mpf(5) / 3) * (d2[0] + d2[1] + d2[2]) / vand2 ** (mpf(1) / 3)
    ratio = chi / mark_expr                                         # must be a cube root of unity
    rows[name] = dict(
        tau=nstr(tau, 6),
        digits_vandermonde2_vs_pi12_eta24_dirac8=digits(vand2, rhs),
        digits_dirac_product_equals_2=digits(dirac_prod, mpf(2)),
        digits_crossratio_is_lambda=max(digits(c, lam) for c in cr),
        digits_sum_d2_vs_2pi4_E4=digits(d2[0] + d2[1] + d2[2], 2 * pi ** 4 * eis(tau, 4)),
        chi_over_mark_expression=nstr(ratio, 12),
        digits_ratio_cubed_equals_1=digits(ratio ** 3, mpf(1)),
        digits_T_phase=digits(chi_dirac(tau + 1) / chi, exp(-2j * pi / 3)),
        digits_S_invariance=digits(chi_dirac(-1 / tau), chi),
    )
res["points"] = rows

# what fixes the weight w of the spinor sector?  scan w: modular S defect and vacuum multiplicity
tau = points["complex 1"]
scan = {}
for w in (0, mpf(1) / 2, 1, 2):
    c = chi_dirac(tau, mpf(w))
    scan[str(w)] = dict(S_defect=nstr(abs(chi_dirac(-1 / tau, mpf(w)) - c) / abs(c), 6))
res["weight_scan_spinor_sector"] = scan

# exact integer coefficients up to q^24, with and without the spinor sector
N = 26
def mul(a, b):
    out = [0] * N
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if i + j < N and y:
                    out[i + j] += x * y
    return out
def power(a, k):
    out = [1] + [0] * (N - 1)
    for _ in range(k):
        out = mul(out, a)
    return out
# series in q:  theta3^8 + theta4^8 has only integer powers; theta2^8 = 256 q (sum q^{n(n+1)/2})^8
th3 = [0] * (2 * N); th4 = [0] * (2 * N)          # in t = q^(1/2)
n = 0
while n * n < 2 * N:
    th3[n * n] += 1 if n == 0 else 2
    th4[n * n] += (1 if n == 0 else 2) * (-1) ** n
    n += 1
def powt(a, k):
    out = [1] + [0] * (2 * N - 1)
    for _ in range(k):
        nw = [0] * (2 * N)
        for i, x in enumerate(out):
            if x:
                for j, y in enumerate(a):
                    if i + j < 2 * N and y:
                        nw[i + j] += x * y
        out = nw
    return out
ns = [(a + b) // 2 for a, b in zip(powt(th3, 8), powt(th4, 8))][::2][:N]     # (th3^8+th4^8)/2
tri = [0] * N
n = 0
while n * (n + 1) // 2 < N:
    tri[n * (n + 1) // 2] += 1; n += 1
r = [0] + [128 * x for x in power(tri, 8)][:N - 1]                           # th2^8 / 2
inv_eta8 = [1] + [0] * (N - 1)                                               # prod (1-q^n)^-8
part = [1] + [0] * (N - 1)
for k in range(1, N):
    for i in range(k, N):
        part[i] += part[i - k]
inv_eta8 = power(part, 8)
with_sp = mul([a + b for a, b in zip(ns, r)], inv_eta8)[:25]
without = mul(ns, inv_eta8)[:25]
# independent reference: E4 / eta^8
sig3 = [1] + [240 * sum(d ** 3 for d in range(1, m + 1) if m % d == 0) for m in range(1, N)]
ref = mul(sig3, inv_eta8)[:25]
res["character_coefficients_q0_to_q5"] = with_sp[:6]
res["matches_E4_over_eta8_to_q24"] = with_sp == ref
res["weight_one_currents"] = dict(with_spinor=with_sp[1], without_spinor=without[1],
                                  split="112 D8 roots + 8 Cartan = 120, plus 128 spinor weights = 248")
res["reading"] = (
    "Mark side fixes prod d_ij^2 = pi^12 eta^24 prod_nu (det D_nu)^8 and the three sector ratios "
    "theta2^8 : theta3^8 : theta4^8 = d_23^2 : d_13^2 : d_12^2.  Not fixed by the determinant relation: "
    "(i) the overall weight 1/2 of the sector sum, (ii) the cube-root branch = central charge phase c/24 = 1/3, "
    "(iii) S+ versus S- of the spinor.  (i) follows from relabelling symmetry of the three finite marks plus a "
    "unique vacuum; the S defect vanishes only at spinor weight 1.  (ii) follows from T.  (iii) is invisible in "
    "the unmarked character.")
json.dump(res, open("out/f1_marks_dirac_weight.json", "w"), indent=1)
print(json.dumps(res, indent=1))
