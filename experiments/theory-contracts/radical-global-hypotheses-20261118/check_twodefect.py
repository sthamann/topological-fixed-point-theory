#!/usr/bin/env python3
"""RADICALTWODEFECT20261118: Zweidefektgewicht ohne Blick auf Alpha.
Vorab festgelegt: (1) Standardzaehlungen eines Defektgases, (2) Zweipunktfunktionen der
Primaerfelder von so(10)_1, su(4)_1 und (E8)_1. Erst danach Vergleich mit 5/4 und CODATA.
Reproduktion: python3 -B check_twodefect.py"""
import itertools, json, sys
from fractions import Fraction as Fr
from pathlib import Path
import mpmath as mp
import sympy as sp
mp.mp.dps = 30
R = {}; FAILS = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + (" | " + str(info) if info != "" else ""))
    if not ok: FAILS.append(name)
q = sp.symbols("q"); g = 5
def c2(f):  # coefficient of q^2 after normalising the q^1 coefficient to one
    s = sp.series(f, q, 0, 4).removeO(); return sp.nsimplify(s.coeff(q, 2)/s.coeff(q, 1))
# T1 standard countings of a defect gas (fixed list)
models = {
 "log Z verdünntes Gas (q)": q,
 "Z - 1 ununterscheidbar (e^q - 1)": sp.exp(q)-1,
 "Bose Besetzung einer Mode q/(1-q) [Fassung 2.4]": q/(1-q),
 "Fermi Besetzung q/(1+q)": q/(1+q),
 "Z - 1 für g Bosemoden (1-q)^-g - 1": (1-q)**(-g)-1,
 "Z - 1 für g Fermimoden (1+q)^g - 1": (1+q)**g-1,
 "geordnete Paare mit g Zuständen, harter Kern: q + (g-1)/g q^2": q+sp.Rational(g-1, g)*q**2,
 "ungeordnete Paare / g^2: q + (g+1)/(2g) q^2": q+sp.Rational(g+1, 2*g)*q**2,
 "Z2 Doppelüberzug, Bose: q/(1-q^2)": q/(1-q**2),
}
T1 = {k: c2(f) for k, f in models.items()}; R["standard_countings"] = {k: str(v) for k, v in T1.items()}
for k, v in T1.items(): print(f"   {k:62s} {v}")
check("T1 none of nine standard defect countings gives the q^2 coefficient 5/4", all(v != sp.Rational(5, 4) for v in T1.values()), sorted({str(v) for v in T1.values()}))
# T2 Haldane exclusion statistics: n = q + (1 - 2 gamma) q^2 ; 5/4 needs gamma = -1/8 (outside 0..1)
check("T2 exclusion statistics would need the parameter gamma = -1/8, outside the physical range [0,1]", sp.solve(sp.Eq(1-2*sp.Symbol("x"), sp.Rational(5, 4)))[0] == sp.Rational(-1, 8))
# T3 conformal weights: Sugawara h = C2/(k + h_dual) at level one, and minimal coset norms by enumeration
def h_sug(C2, hd): return Fr(C2)/(1+hd)
hs = {"so(10)_1 Vektor": h_sug(Fr(9, 2), 8), "so(10)_1 Spinor": h_sug(Fr(45, 8), 8),
      "su(4)_1 Fundamental": h_sug(Fr(15, 8), 4), "su(4)_1 Antisymmetrisch": h_sug(Fr(5, 2), 4), "(E8)_1 Strom": Fr(1)}
def min_norm_D5_spinor():
    return min(sum(Fr(s, 2)**2 for s in v) for v in itertools.product((1, -1, 3, -3), repeat=5) if sum((x-1)//2 for x in v) % 2 == 0)
def min_norm_A3_fund():  # A3 weight coset (3/4,-1/4,-1/4,-1/4) + root lattice
    best = None
    for n in itertools.product(range(-2, 3), repeat=4):
        if sum(n) == 0:
            v = [Fr(3, 4)+n[0], Fr(-1, 4)+n[1], Fr(-1, 4)+n[2], Fr(-1, 4)+n[3]]; nn = sum(x*x for x in v)
            best = nn if best is None or nn < best else best
    return best
nD, nA = min_norm_D5_spinor(), min_norm_A3_fund()
R["weights"] = {k: str(v) for k, v in hs.items()}; R["glue_norms"] = [str(nD), str(nA)]
check("T3 so(10)_1 spinor h = 5/8, su(4)_1 fundamental h = 3/8, sum = 1; glue norms 5/4 and 3/4 equal 2h",
      hs["so(10)_1 Spinor"] == Fr(5, 8) and hs["su(4)_1 Fundamental"] == Fr(3, 8) and nD == Fr(5, 4) and nA == Fr(3, 4) and nD == 2*hs["so(10)_1 Spinor"] and nA == 2*hs["su(4)_1 Fundamental"], (str(nD), str(nA)))
# T4 two-point function of a lattice vertex operator from the oscillator sum: exp(l^2 sum q^n/n) = (1-q)^(-l^2)
qq = mp.mpf("0.3"); l2 = mp.mpf(5)/4
osc = mp.e**(l2*mp.nsum(lambda n: qq**n/n, [1, mp.inf]))
check("T4 <V_l(1) V_-l(q)> = exp(l^2 sum q^n/n) = (1-q)^(-l^2), checked at q = 0.3 for l^2 = 5/4", abs(osc-(1-qq)**(-l2)) < mp.mpf("1e-25"), mp.nstr(osc, 12))
# T5 the series of q (1-q)^(-2h): coefficients 1, 2h, 2h(2h+1)/2
ser = sp.series(q*(1-q)**(-sp.Rational(5, 4)), q, 0, 5).removeO()
check("T5 q (1-q)^(-5/4) = q + 5/4 q^2 + 45/32 q^3 + 195/128 q^4", [ser.coeff(q, n) for n in (1, 2, 3, 4)] == [1, sp.Rational(5, 4), sp.Rational(45, 32), sp.Rational(195, 128)])
# T6 only now: compare the native exponents 2h with CODATA
pi = mp.pi; c3 = 1/(8*pi); pt = 1/(6*pi); d = 48*c3**4; b1 = mp.mpf(41)/10; A = mp.mpf("137.035999177"); S = mp.mpf("2.1e-8")
def pull(f):
    a = mp.findroot(lambda a: a**3-2*c3**3*a**2-8*b1*c3**6*mp.log(1/(pt+f(d*mp.e**(-2*a)))), mp.mpf(1)/137); return float((1/a-A)/S)
P = {k: pull(lambda x, s=mp.mpf(2*h.numerator)/h.denominator: x*(1-x)**(-s)) for k, h in hs.items()}
R["pulls_native_exponents"] = P
for k, v in P.items(): print(f"   2h = {str(2*hs[k]):4s} {k:26s} Zug {v:+8.2f} sigma")
check("T6 among the native exponents 2h in {3/4, 1, 5/4, 2} only the so(10)_1 spinor (5/4) lies within 3 sigma of CODATA",
      [k for k, v in P.items() if abs(v) < 3] == ["so(10)_1 Spinor"], {k: round(v, 1) for k, v in P.items()})
# T7 the full E8 glue vector factorises: (1-q)^-2 = (1-q)^-5/4 (1-q)^-3/4 ; the full vector misses
check("T7 the E8 glue two-point function (1-q)^-2 factorises into the D5 and A3 factors; taken whole it misses CODATA by more than 100 sigma",
      sp.simplify((1-q)**(-2)-(1-q)**(-sp.Rational(5, 4))*(1-q)**(-sp.Rational(3, 4))) == 0 and abs(P["(E8)_1 Strom"]) > 100, round(P["(E8)_1 Strom"], 1))
# T8 alpha cannot discriminate completions that share the q^2 coefficient
alt = {"(1-q)^-5/4": lambda x: x*(1-x)**(-mp.mpf(5)/4), "abgebrochen q + 5/4 q^2": lambda x: x+mp.mpf(5)/4*x*x,
       "Produkt über n (1-q^n)^-5/4": lambda x: x*mp.nprod(lambda n: (1-x**n)**(-mp.mpf(5)/4), [1, 40])}
PA = {k: pull(f) for k, f in alt.items()}; R["completions"] = PA
check("T8 three completions with the same q^2 coefficient differ by less than 0.1 sigma: alpha fixes one number, not the function",
      max(PA.values())-min(PA.values()) < 0.1, {k: round(v, 3) for k, v in PA.items()})
# T9 is q a native cross ratio? The four marks 1, i, -1, -i have cross ratios {-1, 2, 1/2}; q is about 1.19e-4
al = mp.findroot(lambda a: a**3-2*c3**3*a**2-8*b1*c3**6*mp.log(1/(pt+alt["(1-q)^-5/4"](d*mp.e**(-2*a)))), mp.mpf(1)/137); qv = d*mp.e**(-2*al)
z = [1, 1j, -1, -1j]; cr = (z[0]-z[1])*(z[2]-z[3])/((z[0]-z[2])*(z[1]-z[3])); R["q_value"] = float(qv); R["mark_cross_ratio"] = [float(mp.re(cr)), float(mp.im(cr))]
check("T9 q = 1.185e-4 is not the cross ratio of the four marks (1/2): the insertion separation is not supplied by the mark geometry",
      abs(qv-mp.mpf("1.185e-4")) < 1e-6 and abs(cr-0.5) < 1e-12, (mp.nstr(qv, 6), cr))
R["checks_passed"] = not FAILS
Path(__file__).with_name("results_twodefect.json").write_text(json.dumps(R, indent=1, ensure_ascii=False) + "\n")
print("RESULT", "ALL PASS" if not FAILS else "FAIL " + ", ".join(FAILS)); sys.exit(1 if FAILS else 0)
