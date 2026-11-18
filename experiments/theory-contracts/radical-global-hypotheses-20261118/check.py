#!/usr/bin/env python3
"""RADICAL20261118: Gegenproben zu sechs Gesamthypothesen (exakt/numerisch).
Reproduktion: python3 -B check.py   (schreibt results.json, Exit 0 bei allen PASS)
Keine Ledger- oder Statusaenderung. Messwerte sind aus dem Gedaechtnis zitierte
Literaturwerte (PDG 2024, NuFIT 6.0, CODATA 2022, Planck 2018) und im README benannt."""
import json, itertools, math, sys
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 40
R = {}; FAILS = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + (" | " + str(info) if info != "" else ""))
    if not ok: FAILS.append(name)

pi = mp.pi
c3 = 1/(8*pi)
phi0 = 1/(6*pi) + 48*c3**4

# ---------------------------------------------------------------- A: alpha root
def phis(al, k=6, m=48, e=2, p=mp.mpf(5)/4):
    q = m*c3**4*mp.e**(-e*al)
    return 1/(k*pi) + q*(1-q)**(-p)
def F(al, t=2, K=mp.mpf(164)/5, **kw):
    return al**3 - t*c3**3*al**2 - K*c3**6*mp.log(1/phis(al, **kw))
al = mp.findroot(F, mp.mpf(1)/137)
ainv = 1/al
R["alpha_inv_tfpt"] = mp.nstr(ainv, 16)
check("A1 alpha root reproduces 137.0359992168407", abs(ainv - mp.mpf("137.0359992168407")) < mp.mpf("1e-12"), mp.nstr(ainv, 16))

# ---------------------------------------------------------------- B: data pulls
lamC = mp.sqrt(phi0*(1-phi0))
data = [  # name, tfpt, exp, sigma, source
 ("alpha^-1 (CODATA 2022)", ainv, 137.035999177, 2.1e-8),
 ("alpha^-1 (Rb 2020)", ainv, 137.035999206, 1.1e-8),
 ("alpha^-1 (Cs 2018)", ainv, 137.035999046, 2.7e-8),
 ("lambda_C vs |V_us| PDG", lamC, 0.2243, 0.0008),
 ("lambda_C vs Wolfenstein fit", lamC, 0.22501, 0.00068),
 ("s23 vs |V_cb| avg", mp.mpf("0.043428"), 0.0411, 0.0012),
 ("s13 vs |V_ub| avg", mp.mpf("0.003765"), 0.00382, 0.00020),
 ("delta_CKM deg", mp.mpf("68.65"), 65.7, 3.0),
 ("sin2 th12", mp.mpf("0.306747"), 0.308, 0.012),
 ("sin2 th13", mp.mpf("0.023108"), 0.02215, 0.00057),
 ("Omega_b", (1-1/(4*pi))*phi0, 0.0493, 0.0009),
 ("beta deg (ACT DR6)", phi0/(4*pi)*180/pi, 0.20, 0.08),
 ("beta deg (Planck PR4)", phi0/(4*pi)*180/pi, 0.30, 0.11),
]
pulls = []
for n, t, x, s in data:
    pl = float((t-x)/s); pulls.append((n, float(t), x, s, round(pl, 2)))
    print(f"   {n:32s} TFPT {float(t):.10g}  exp {x} +- {s}  pull {pl:+.2f}")
R["pulls"] = pulls
maxpull = max(abs(p[4]) for p in pulls if "Cs" not in p[0])
check("B1 no pull beyond 3 sigma (Cs alpha excluded, it conflicts with Rb itself)", maxpull < 3, maxpull)
me_mu = mp.mpf("0.00483633170"); mmu_tau = mp.mpf("105.6583755")/mp.mpf("1776.93")
d1 = float(Fr(12,7)*phi0**2/me_mu - 1); d2 = float(Fr(8,7)*phi0/mmu_tau - 1)
R["lepton_ratio_dev"] = {"me/mmu": d1, "mmu/mtau": d2}
print(f"   Leptonquellverhaeltnisse gegen Polmassen: me/mmu {d1:+.4%}, mmu/mtau {d2:+.4%}")
check("B2 lepton source ratios deviate from pole ratios by 0.1..3 percent (transfer needed)", 1e-3 < abs(d1) < 3e-2 and 1e-3 < abs(d2) < 3e-2)

# ---------------------------------------------------------------- C: look-elsewhere on alpha
a_exp = 1/mp.mpf("137.035999177"); sig = mp.mpf("2.1e-8")
a_lo = 1/(mp.mpf("137.035999177")+2*sig); a_hi = 1/(mp.mpf("137.035999177")-2*sig)
a_tf = al
def Kreq(a, t, k, m, e, p):
    return (a**3 - t*c3**3*a**2)/(c3**6*mp.log(1/phis(a, k=k, m=m, e=e, p=p)))
mp.mp.dps = 25
hits = []; expected = 0.0; combos = 0
BMAX = 8
dens = sum(1 for b in range(1, BMAX+1) for a_ in range(b) if math.gcd(a_, b) == 1)  # reduced fractions per unit
for t in range(0, 5):
  for k in range(1, 13):
    for m in range(1, 97):
      for e in range(0, 5):
        for p4 in range(0, 9):
            p = mp.mpf(p4)/4
            try:
                k1 = Kreq(a_lo, t, k, m, e, p); k2 = Kreq(a_hi, t, k, m, e, p)
            except Exception:
                continue
            lo, hi = (k1, k2) if k1 < k2 else (k2, k1)
            if not (0 < lo and hi < 200): continue
            combos += 1
            expected += float(hi-lo)*dens
            for b in range(1, BMAX+1):
                A = mp.ceil(lo*b)
                if A <= hi*b and math.gcd(int(A), b) == 1:
                    hits.append((t, k, m, e, float(p), f"{int(A)}/{b}"))
mp.mp.dps = 40
R["lookelsewhere"] = {"combos": combos, "expected_chance_hits": expected, "hits": hits, "window": "CODATA 2022 +-2 sigma", "denominators<=": BMAX}
print(f"   Grammatik: {combos} Kombinationen, erwartete Zufallstreffer {expected:.3f}, tatsaechliche Treffer {len(hits)}: {hits}")
tf = (2, 6, 48, 2, 1.25, "164/5")
check("C1 TFPT tuple is inside CODATA 2 sigma window", tf in hits, len(hits))
check("C2 chance expectation in this neighbourhood grammar is order one (0.05..5)", 0.05 < expected < 5, round(expected, 3))

# ---------------------------------------------------------------- D: self-dual point
def eis(tau):
    q = mp.e**(2j*pi*tau)
    s1 = mp.nsum(lambda n: n*q**n/(1-q**n), [1, mp.inf])
    s3 = mp.nsum(lambda n: n**3*q**n/(1-q**n), [1, mp.inf])
    s5 = mp.nsum(lambda n: n**5*q**n/(1-q**n), [1, mp.inf])
    return 1-24*s1, 1+240*s3, 1-504*s5
E2, E4, E6 = eis(1j)
check("D1 E2(i) = 3/pi", abs(E2 - 3/pi) < mp.mpf("1e-30"))
check("D2 E6(i) = 0", abs(E6) < mp.mpf("1e-30"))
check("D3 c3 = E2(i)/24 exactly (1/(8 pi))", abs(E2/24 - c3) < mp.mpf("1e-30"))
check("D4 phi0 tree term 1/(6 pi) = E2(i)/18 = (4/3) c3", abs(E2/18 - 1/(6*pi)) < mp.mpf("1e-30"))
eta = lambda tau: mp.e**(2j*pi*tau/24)*mp.nprod(lambda n: 1-mp.e**(2j*pi*tau*n), [1, mp.inf])
chi_i = E4/eta(1j)**8
check("D5 chi_E8(i) = E4/eta^8 = 12", abs(chi_i - 12) < mp.mpf("1e-25"), mp.nstr(chi_i.real, 12))
# energy  q d/dq log chi = -E6/(3 E4)
def energy_chi(tau):
    e2, e4, e6 = eis(tau); return -e6/(3*e4)
check("D6 vacuum-character energy vanishes at tau=i (stationary), entropy = log 12", abs(energy_chi(1j)) < mp.mpf("1e-30"), "S=" + mp.nstr(mp.log(12), 10))
# lattice charge energy per Cartan direction via theta functions, q = exp(2 pi i tau), weight q^{x^2/2}
def th(tau):
    q = mp.e**(2j*pi*tau)
    f = lambda g, s: mp.nsum(lambda n: g(n)*q**(((n+s)**2)/2), [-mp.inf, mp.inf])
    t3 = f(lambda n: 1, 0); d3 = f(lambda n: n**2/2, 0)
    t4 = f(lambda n: (-1)**int(n), 0); d4 = f(lambda n: (-1)**int(n)*n**2/2, 0)
    t2 = f(lambda n: 1, mp.mpf(1)/2); d2 = f(lambda n: (n+mp.mpf(1)/2)**2/2, mp.mpf(1)/2)
    return t2, d2, t3, d3, t4, d4
t2, d2, t3, d3, t4, d4 = th(1j)
def Dn(n):   # root lattice D_n: (t3^n + t4^n)/2 and its q d/dq
    return (t3**n+t4**n)/2, (n*t3**(n-1)*d3 + n*t4**(n-1)*d4)/2
def E8():
    v, dv = Dn(8); return v + t2**8/2, dv + 8*t2**7*d2/2
lat = {}
v, dv = E8(); lat["E8"] = (dv/v/8).real
v, dv = Dn(8); lat["D8"] = (dv/v/8).real
v5, dv5 = Dn(5); v3, dv3 = Dn(3); lat["D5+A3 (unverklebt)"] = ((dv5/v5 + dv3/v3)/8).real
lat["Z^8"] = (d3/t3).real
for k_, x in lat.items(): print(f"   <p^2/2> pro Richtung bei tau=i, {k_:20s}: {mp.nstr(x, 15)}   / c3 = {mp.nstr(x/c3, 12)}")
R["charge_energy_over_c3"] = {k_: float(x/c3) for k_, x in lat.items()}
check("D7 E8: lattice charge energy per Cartan direction = c3 exactly", abs(lat["E8"]-c3) < mp.mpf("1e-28"))
check("D8 unglued D5+A3 and D8 miss c3; self-dual Z^8 hits it (self-duality, not evenness)", abs(lat["D8"]-c3) > 1e-4 and abs(lat["D5+A3 (unverklebt)"]-c3) > 1e-4 and abs(lat["Z^8"]-c3) < mp.mpf("1e-28"))
# torus current two-point constant: C(w) = -wp(w) + (<p^2>_dir - E2/12) = -wp(w) - E6/(12 E4)
def const_direct(tau):
    a = th(tau); t2_, d2_, t3_, d3_, t4_, d4_ = a
    v = (t3_**8+t4_**8+t2_**8)/2; dv = (8*t3_**7*d3_+8*t4_**7*d4_+8*t2_**7*d2_)/2
    e2, e4, e6 = eis(tau)
    return 2*dv/v/8 - e2/12, -e6/(12*e4)
ok = True
for tau in (1j, 1.3j, 0.2+0.9j):
    a, b = const_direct(tau); ok &= abs(a-b) < mp.mpf("1e-25")
    print(f"   Stromkonstante tau={tau}: {mp.nstr(a, 12)}")
check("D9 current two-point constant = -E6/(12 E4) for three moduli; zero at tau=i", ok and abs(const_direct(1j)[0]) < mp.mpf("1e-28"))
# direct check of the Fourier form against an independent lattice sum of wp at tau=i
import numpy as np
w = 1.1-0.7j; N = 400  # radiale Ordnung: -2pi < Im w < 0
nn = np.arange(-N, N+1); om = 2*np.pi*(nn[:, None] + 1j*nn[None, :]); om = om[om != 0]
wp = 1/w**2 + np.sum(1/(w-om)**2 - 1/om**2)
q = math.exp(-2*math.pi)
C = float(lat["E8"])*2 + sum(n*(np.exp(-1j*n*w)/(1-q**n) + q**n*np.exp(1j*n*w)/(1-q**n)) for n in range(1, 60))
R["torus_two_point_vs_minus_wp"] = [complex(C).real, complex(C).imag, complex(-wp).real, complex(-wp).imag]
check("D10 mode sum of <J J> on the tau=i torus equals -wp(w) (independent lattice sum)", abs(C + wp) < 2e-3, abs(C+wp))

# antipodal position factor: plane contract gives -1/4, the self-dual torus gives -e1 = -(varpi/(2 pi))^2
qi = mp.e**(-2*pi)
wp_pi = mp.mpf(1)/4 - E2.real/12 - 2*mp.nsum(lambda n: n*qi**n/(1-qi**n)*(-1)**int(n), [1, mp.inf])
varpi = mp.gamma(mp.mpf(1)/4)**2/(2*mp.sqrt(2*pi))
R["antipodal_kernel"] = {"plane": -0.25, "torus_tau_i": float(-wp_pi), "minus_(varpi/2pi)^2": float(-(varpi/(2*pi))**2)}
print("   Antipodaler Zweipunktkern: Ebene -1/4, Torus tau=i", mp.nstr(-wp_pi, 15))
check("D11 antipodal kernel on the tau=i torus is -(varpi/(2 pi))^2 = -0.17415..., not -1/4", abs(wp_pi - (varpi/(2*pi))**2) < mp.mpf("1e-30"), mp.nstr(-wp_pi, 12))

# ---------------------------------------------------------------- E: channel hypothesis
import sympy as sp
x = sp.symbols("x")
anh = [x, 1-x, 1/x, 1/(1-x), (x-1)/x, x/(x-1)]
def stab(val): return sum(1 for f in anh if sp.simplify(f.subs(x, val)-val) == 0)
rho = sp.Rational(1, 2)+sp.sqrt(3)*sp.I/2
R["stabilizers"] = {"harmonic": [stab(v) for v in (-1, sp.Rational(1, 2), 2)], "rho": stab(rho), "generic": stab(sp.Rational(3, 7))}
check("E1 anharmonic stabiliser: harmonic points Z2, rho Z3, generic trivial", R["stabilizers"] == {"harmonic": [2, 2, 2], "rho": 3, "generic": 1})
pair = [frozenset({frozenset({0, 1}), frozenset({2, 3})}), frozenset({frozenset({0, 2}), frozenset({1, 3})}), frozenset({frozenset({0, 3}), frozenset({1, 2})})]
ker = 0; imgs = set()
for s in itertools.permutations(range(4)):
    im = tuple(pair.index(frozenset(frozenset(s[i] for i in blk) for blk in pp)) for pp in pair)
    imgs.add(im); ker += im == (0, 1, 2)
check("E2 S4=W(A3) on four marks -> S3 on the three pairings with kernel V4", len(imgs) == 6 and ker == 4)
Rm = [[1, 3, 0], [1, 5, 2], [2, 5, 3]]
sym = 0
for P in itertools.permutations(range(3)):
    for Q in itertools.permutations(range(3)):
        if [[Rm[P[i]][Q[j]] for j in range(3)] for i in range(3)] == Rm: sym += 1
R["R_perm_symmetries"] = sym
check("E3 flavour matrix R has no family exchange symmetry (only identity)", sym == 1, sym)

# ---------------------------------------------------------------- F: carrier scale = scalaron scale?
MZ = 91.1876; Mbar = 2.435e18; aem = 1/127.951; s2 = 0.23122
a1 = (5/3)*aem/(1-s2); a2 = aem/s2
b1, b2 = 41/10, -19/6
L = (1/a1 - 1/a2)*2*math.pi/(b1-b2)
mu_x = MZ*math.exp(L); Mscal = float(c3**mp.mpf(3.5))*Mbar
R["scales"] = {"mu(sin2=3/8) 1-loop SM GeV": mu_x, "M_scal GeV": Mscal, "ratio": Mscal/mu_x, "alpha_GUT_inv": 1/a2 - b2/(2*math.pi)*L}
print(f"   sin^2 thetaW = 3/8 (SM, 1 Schleife) bei {mu_x:.3e} GeV; M_scal = c3^(7/2) Mbar = {Mscal:.3e} GeV; Verhaeltnis {Mscal/mu_x:.2f}")
check("F1 scalaron scale and 3/8 crossing agree within a factor 5 (one loop, no thresholds)", 0.2 < Mscal/mu_x < 5, round(Mscal/mu_x, 2))

json.dump(R, open("results.json", "w"), indent=1, default=str)
print("\nFAILS:", FAILS)
sys.exit(1 if FAILS else 0)
