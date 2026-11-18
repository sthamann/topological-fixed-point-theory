#!/usr/bin/env python3
"""RADICALHIST20261118: Proben an den fruehen Fassungen (Theory-Root, UFE, ImprovedCubic, Eliminating K, V25-C).
Reproduktion: python3 -B check_historical.py"""
import json, sys
import mpmath as mp
mp.mp.dps = 30
pi = mp.pi; c3 = 1/(8*pi); b1 = mp.mpf(41)/10; R = {}; FAILS = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + (" | " + str(info) if info != "" else ""))
    if not ok: FAILS.append(name)
phit = 1/(6*pi); dtop = 48*c3**4; phi0 = phit+dtop
def alpha_of(phi_fn):
    F = lambda a: a**3-2*c3**3*a**2-8*b1*c3**6*mp.log(1/phi_fn(a))
    return mp.findroot(F, mp.mpf(1)/137)
A_CODATA = mp.mpf("137.035999177"); SIG = mp.mpf("2.1e-8")
# I1: the three historical closures of phi
var = {"Baum + delta (V1.0)": lambda a: phi0,
       "delta (1 - 2 alpha) (ImprovedCubic)": lambda a: phit+dtop*(1-2*a),
       "delta exp(-2 alpha)": lambda a: phit+dtop*mp.e**(-2*a),
       "aktuell q (1-q)^(-5/4)": lambda a: phit+(dtop*mp.e**(-2*a))*(1-dtop*mp.e**(-2*a))**(-mp.mpf(5)/4)}
R["alpha_inv"] = {}
for k, f in var.items():
    ai = 1/alpha_of(f); R["alpha_inv"][k] = [mp.nstr(ai, 15), float((ai-A_CODATA)/SIG)]
    print(f"   {k:38s} alpha^-1 = {mp.nstr(ai, 15)}   Zug {float((ai-A_CODATA)/SIG):+.1f} sigma")
z = [abs(v[1]) for v in R["alpha_inv"].values()]
check("I1 only the current closure is within 3 sigma of CODATA 2022; each earlier one is off by more than 100 sigma", z[3] < 3 and min(z[:3]) > 100, [round(x, 1) for x in z])
# I2: 8 b1 = (4/5) 41 and b1 = (4/3) N_fam + N_H/10, the one-loop U(1)_Y coefficient
check("I2 8 b1 = 164/5 and b1 = (4/3)*3 + 1/10 (three families, one Higgs)", 8*b1 == mp.mpf(164)/5 and b1 == mp.mpf(4)/3*3+mp.mpf(1)/10)
# I3: Coleman-Weinberg variable y = 1/(64 pi^2 alpha) = c3^2/alpha:  8 b1 L y^3 + y/(4 pi) = 1
al = alpha_of(var["aktuell q (1-q)^(-5/4)"]); L = mp.log(1/var["aktuell q (1-q)^(-5/4)"](al)); y = c3**2/al
check("I3 cubic in y = 1/(64 pi^2 alpha): 8 b1 L y^3 + y/(4 pi) = 1", abs(8*b1*L*y**3+y/(4*pi)-1) < mp.mpf("1e-25"), mp.nstr(y, 12))
# I4: one-loop running would make 1/alpha linear in the logarithm; the cubic gives alpha^3 linear in it
d = mp.mpf("1e-6")
a1 = mp.findroot(lambda a: a**3-2*c3**3*a**2-8*b1*c3**6*(L+d), al); a0 = mp.findroot(lambda a: a**3-2*c3**3*a**2-8*b1*c3**6*L, al)
slope = (1/a1-1/a0)/d; rg = b1/(2*pi)
R["d(alpha^-1)/dL"] = float(slope); R["one_loop_b1/2pi"] = float(rg)
check("I4 d(1/alpha)/dL = -15.4, a one-loop U(1) running would give +0.65: not a perturbative running relation", abs(slope+mp.mpf("15.4")) < 0.2 and abs(rg-0.6525) < 1e-3, (mp.nstr(slope, 6), mp.nstr(rg, 6)))
# I5: alpha as a field: U(a) = a^4/4 - (2/3) c3^3 a^3 - 8 b1 c3^6 L a ; stationary point is a minimum
Upp = 3*al**2-4*c3**3*al
check("I5 the root is a minimum of the historical potential U(alpha)", Upp > 0, mp.nstr(Upp, 6))
# I6: historical k cubic (V25-C): k^3 - 2 c3^2 k^2 - 8 b1 c3^6 L = 0, claim alpha <-> 2 pi k
k = mp.findroot(lambda k_: k_**3-2*c3**2*k_**2-8*b1*c3**6*L, mp.mpf("0.0085"))
R["k"] = mp.nstr(k, 12); R["2pi k / phi0"] = float(2*pi*k/phi0)
print("   k =", mp.nstr(k, 12), " 2 pi k =", mp.nstr(2*pi*k, 12), " phi0 =", mp.nstr(phi0, 12))
check("I6 2 pi k is 0.17 percent above phi0 (not alpha): the k cubic nearly reproduces the seed", abs(2*pi*k/phi0-1) < 0.003 and abs(2*pi*k/phi0-1) > 0.001, mp.nstr(2*pi*k/phi0, 10))
# I7: closing the pair self-consistently, phi := 2 pi k and L = log(1/phi), then alpha from the alpha cubic
ks = mp.findroot(lambda k_: k_**3-2*c3**2*k_**2-8*b1*c3**6*mp.log(1/(2*pi*k_)), mp.mpf("0.0085"))
ps = 2*pi*ks; ai = 1/alpha_of(lambda a: ps)
R["selfconsistent_pair"] = {"phi": mp.nstr(ps, 12), "alpha_inv": mp.nstr(ai, 12), "pull": float((ai-A_CODATA)/SIG)}
print("   geschlossenes Paar: phi =", mp.nstr(ps, 12), " alpha^-1 =", mp.nstr(ai, 12))
check("I7 the self-consistent pair (phi = 2 pi k) misses alpha by more than 1000 sigma", abs((ai-A_CODATA)/SIG) > 1000, float((ai-A_CODATA)/SIG))
# I8: Eliminating K: xi = c3/phi0, tree value 3/4; geometric action 12 pi phi_tree = 2 = chi(S^2)
check("I8 xi_tree = c3/phi_tree = 3/4 and 12 pi phi_tree = 2", abs(c3/phit-mp.mpf(3)/4) < mp.mpf("1e-28") and abs(12*pi*phit-2) < mp.mpf("1e-28"), mp.nstr(c3/phi0, 8))
json.dump(R, open("results_historical.json", "w"), indent=1)
print("\nFAILS:", FAILS); sys.exit(1 if FAILS else 0)
