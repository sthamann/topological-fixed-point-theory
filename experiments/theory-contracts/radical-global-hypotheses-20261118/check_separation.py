#!/usr/bin/env python3
"""RADICALSEPARATION20261118: Ist q = 48 c3^4 exp(-2 alpha) ein Einsetzungsabstand, und warum nur D5?
Vorab festgelegt: Zweipunktfunktionen auf dem Zylinder (Naht = Kreis mal Zeit) und eine Liste nativer Zeiten.
Reproduktion: python3 -B check_separation.py"""
import json, sys
from pathlib import Path
import mpmath as mp
import sympy as sp
mp.mp.dps = 30
R = {}; FAILS = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + (" | " + str(info) if info != "" else ""))
    if not ok: FAILS.append(name)
pi = mp.pi; c3 = 1/(8*pi); pt = 1/(6*pi); d = 48*c3**4; b1 = mp.mpf(41)/10; A = mp.mpf("137.035999177"); S = mp.mpf("2.1e-8")
def ainv(f, dd=d):
    return 1/mp.findroot(lambda a: a**3-2*c3**3*a**2-8*b1*c3**6*mp.log(1/(pt+f(dd*mp.e**(-2*a)))), mp.mpf(1)/137)
pull = lambda f, dd=d: float((ainv(f, dd)-A)/S)
F54 = lambda x: x*(1-x)**(-mp.mpf(5)/4)
# S1 a primary two-point function on the cylinder: (z1 z2)^h (z1-z2)^(-2h) = x^h (1-x)^(-2h), x = z2/z1
h = mp.mpf(5)/8; w1, w2 = mp.mpf("0.9"), mp.mpf("-0.4"); z1, z2 = mp.e**w1, mp.e**w2; x = z2/z1
check("S1 cylinder two-point function of a weight h primary is x^h (1-x)^(-2h) with x = exp(-t)", abs((z1*z2)**h*(z1-z2)**(-2*h)-x**h*(1-x)**(-2*h)) < mp.mpf("1e-25"))
# S2 the D5 spinor alone on the cylinder: q^(5/8) (1-q)^(-5/4)
p2 = pull(lambda x: x**(mp.mpf(5)/8)*(1-x)**(-mp.mpf(5)/4)); R["pull_D5_cylinder"] = p2
check("S2 D5 spinor two-point function on the cylinder, q^(5/8)(1-q)^(-5/4), misses CODATA by more than 10^5 sigma", abs(p2) > 1e5, f"{p2:.3e}")
# S3 the full E8 glue operator (weight 1) on the cylinder: q (1-q)^(-2)
p3 = pull(lambda x: x*(1-x)**(-2)); R["pull_E8_glue_cylinder"] = p3
check("S3 full glue operator on the cylinder, q(1-q)^(-2), misses by +148 sigma", abs(p3-148.1) < 0.5, round(p3, 1))
# S4 the formula is the mixed object: leading weight 1 = 5/8 + 3/8, oscillator exponent of D5 alone
q = sp.symbols("q")
check("S4 q(1-q)^(-5/4) = [q(1-q)^(-2)] (1-q)^(3/4): full glue two-point function with the A3 oscillator factor removed",
      sp.simplify(q*(1-q)**(-sp.Rational(5, 4))-q*(1-q)**(-2)*(1-q)**sp.Rational(3, 4)) == 0)
# S5 an equally short second reading: derivative of a weight 1/8 two-point function
check("S5 q(1-q)^(-5/4) = 4 q d/dq (1-q)^(-1/4): a second reading (weight 1/8 field, derivative in alpha) is equally short",
      sp.simplify(4*q*sp.diff((1-q)**(-sp.Rational(1, 4)), q)-q*(1-q)**(-sp.Rational(5, 4))) == 0)
# S6 what alpha demands of the offset t0 = log(1/delta_top)
t0 = mp.log(1/d); eps = mp.mpf("1e-6")
dz = (pull(F54, d*mp.e**(-eps))-pull(F54))/float(eps); tol = abs(1/dz); R["t0"] = float(t0); R["t0_tolerance_1sigma"] = float(tol)
check("S6 t0 = log(1/delta_top) = 9.02548; one sigma of CODATA corresponds to a shift of t0 by 6e-7", abs(t0-mp.mpf("9.02548")) < 1e-5 and 5e-7 < tol < 7e-7, (mp.nstr(t0, 8), mp.nstr(tol, 3)))
# S7 predeclared native times and separations
varpi = mp.gamma(mp.mpf(1)/4)**2/(2*mp.sqrt(2*pi)); E4i = 3*mp.gamma(mp.mpf(1)/4)**8/(2*pi)**6
cands = {"2 pi (Modul tau = i)": 2*pi, "pi": pi, "3 pi": 3*pi, "4 pi": 4*pi, "2 varpi": 2*varpi, "varpi^2": varpi**2, "4 varpi": 4*varpi,
         "2 pi E4(i)": 2*pi*E4i, "log 240": mp.log(240), "log 248": mp.log(248), "log 2^8": 8*mp.log(2), "log 12 (chi_E8(i))": mp.log(12),
         "8 log 2 + log 30": 8*mp.log(2)+mp.log(30), "log |W(D5)| = log 1920": mp.log(1920), "log(4 pi)^4": 4*mp.log(4*pi), "2 pi sqrt 2": 2*pi*mp.sqrt(2),
         "log(1/Doppelverhältnis 1/2)": mp.log(2), "pi^2": pi**2, "9": mp.mpf(9)}
dev = {k: float(abs(v-t0)) for k, v in cands.items()}; best = min(dev, key=dev.get); R["native_times"] = dev; R["nearest"] = [best, dev[best]]
check("S7 none of 19 predeclared native times equals t0 within 100 sigma of the alpha tolerance", all(v > 100*float(tol) for v in dev.values()), (best, round(dev[best], 4), f"{dev[best]/float(tol):.0f} sigma"))
# S8 exact decompositions of delta_top
volS7 = pi**4/3
check("S8 delta_top = 48 c3^4 = 3/(256 pi^4) = 1/(2^8 vol S^7): four antipodal kernel values times 48, or a density on the unit sphere of the Cartan space",
      abs(d-3/(256*pi**4)) < mp.mpf("1e-30") and abs(d-1/(256*volS7)) < mp.mpf("1e-30") and abs(volS7-2*pi**4/mp.gamma(4)) < mp.mpf("1e-28"))
# S9 q is exponentially small where a separation would be of order one: equal-time marks on the seam give x of order one
xs = {"Gegenpunkt auf dem Kreis (x = -1)": abs((1-(-1))**(-mp.mpf(5)/4)), "Viertelkreis (x = i)": abs((1-1j)**(-mp.mpf(5)/4))}
R["equal_time_values"] = {k: float(v) for k, v in xs.items()}
check("S9 at equal time the mark separations give |1-x|^(-5/4) = 0.420 (antipode) and 0.648 (quarter circle), fixed numbers with no alpha dependence",
      abs(xs["Gegenpunkt auf dem Kreis (x = -1)"]-mp.mpf("0.4204")) < 1e-3 and abs(xs["Viertelkreis (x = i)"]-mp.mpf("0.6484")) < 1e-3, {k: round(float(v), 4) for k, v in xs.items()})
R["checks_passed"] = not FAILS
Path(__file__).with_name("results_separation.json").write_text(json.dumps(R, indent=1, ensure_ascii=False) + "\n")
print("RESULT", "ALL PASS" if not FAILS else "FAIL " + ", ".join(FAILS)); sys.exit(1 if FAILS else 0)
