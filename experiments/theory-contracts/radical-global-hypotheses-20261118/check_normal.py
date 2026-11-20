#!/usr/bin/env python3
"""RADICALNORMAL20261118: Die Alpha Gleichung in der Variable g = alpha/(2 c3^2) = e^2/c3.
Frage: Welche Koeffizientendeutungen (8 = Rang E8, 2 = |Z2|, c3^6) haengen nur an der Wahl der Variable?
Reproduktion: python3 -B check_normal.py"""
import json, sys
from pathlib import Path
import mpmath as mp
import sympy as sp
mp.mp.dps = 30
R = {}; FAILS = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + (" | " + str(info) if info != "" else ""))
    if not ok: FAILS.append(name)
pi = mp.pi; c3 = 1/(8*pi); pt = 1/(6*pi); d = 48*c3**4; A = mp.mpf("137.035999177"); S = mp.mpf("2.1e-8")
def seam(a): q = d*mp.e**(-2*a); return pt+q*(1-q)**(-mp.mpf(5)/4)
def root(b1=mp.mpf(41)/10, phi=seam, c=None):
    c = c3 if c is None else c
    return mp.findroot(lambda a: a**3-2*c**3*a**2-8*b1*c**6*mp.log(1/phi(a)), mp.mpf(1)/137)
al = root(); L = mp.log(1/seam(al)); b1 = mp.mpf(41)/10
# N1 exact rewriting
a, c, b, l = sp.symbols("a c b l", positive=True); g = a/(2*c**2)
check("N1 alpha^3 - 2 c3^3 alpha^2 - 8 b1 c3^6 L = 8 c3^6 [ g^2 (g - c3) - b1 L ] with g = alpha/(2 c3^2)",
      sp.simplify(a**3-2*c**3*a**2-8*b*c**6*l-8*c**6*(g**2*(g-c)-b*l)) == 0)
gv = al/(2*c3**2); R["g"] = mp.nstr(gv, 15)
check("N2 g = 32 pi^2 alpha = e^2/c3 = 2.30470 solves g^2 (g - c3) = b1 log(1/phi): the factors 8 and 2 are gone, only c3 and b1 remain",
      abs(gv-32*pi**2*al) < mp.mpf("1e-25") and abs(gv-4*pi*al/c3) < mp.mpf("1e-25") and abs(gv**2*(gv-c3)-b1*L) < mp.mpf("1e-22") and abs(gv-mp.mpf("2.30470")) < 1e-5, mp.nstr(gv, 10))
# N3 the same seam in g: base (4/3) c3, weight 48 c3^4, damping exp(-4 c3^2 g)
check("N3 phi_seam = (4/3) c3 + q (1-q)^(-5/4), q = 48 c3^4 exp(-4 c3^2 g): five numbers c3, b1, 4/3, 48, 5/4",
      abs(pt-mp.mpf(4)/3*c3) < mp.mpf("1e-30") and abs(mp.e**(-2*al)-mp.e**(-4*c3**2*gv)) < mp.mpf("1e-28"))
# N4 size of each ingredient in 1/alpha
steps = {"g = (b1 L0)^(1/3), L0 = log(6 pi)": 1/(2*c3**2*(b1*mp.log(1/pt))**(mp.mpf(1)/3)),
         "+ c3/3": 1/(2*c3**2*((b1*mp.log(1/pt))**(mp.mpf(1)/3)+c3/3)),
         "exakte Kubik, phi = 1/(6 pi)": 1/root(phi=lambda a: pt),
         "+ 48 c3^4": 1/root(phi=lambda a: pt+d),
         "+ exp(-2 alpha)": 1/root(phi=lambda a: pt+d*mp.e**(-2*a)),
         "+ (1-q)^(-5/4)": 1/al}
R["ladder"] = {k: mp.nstr(v, 13) for k, v in steps.items()}
for k, v in steps.items(): print(f"   {k:34s} {mp.nstr(v, 13):16s} {float((v-A)/A*1e6):+12.4f} ppm")
ppm = [float((v-A)/A*1e6) for v in steps.values()]
check("N4 ladder in ppm from CODATA: +5567, -219, -252, +3.67, -0.037, +0.0003", abs(ppm[0]-5567) < 2 and abs(ppm[1]+218.6) < 1 and abs(ppm[2]+252.0) < 1 and abs(ppm[3]-3.67) < 0.02 and abs(ppm[4]+0.037) < 0.001 and abs(ppm[5]) < 1e-3, [round(x, 4) for x in ppm])
# N5 sensitivities: which inputs does alpha actually test
eps = mp.mpf("1e-8")
sb = float((mp.log(root(b1=b1*(1+eps)))-mp.log(al))/eps); sc = float((mp.log(root(c=c3*(1+eps), phi=lambda a: (4*c3*(1+eps)/3)+48*(c3*(1+eps))**4*mp.e**(-2*a)*(1-48*(c3*(1+eps))**4*mp.e**(-2*a))**(-mp.mpf(5)/4)))-mp.log(al))/eps)
R["dlnalpha/dlnb1"] = sb; R["dlnalpha/dlnc3"] = sc
check("N5 d ln alpha / d ln b1 = 0.331 and d ln alpha / d ln c3 = 1.89: a change of b1 by 1/41 moves alpha by 0.8 %", abs(sb-0.331) < 0.003 and abs(sc-1.89) < 0.01, (round(sb, 4), round(sc, 4)))
# N6 neighbouring field contents
nb = {"zwei Higgs (b1 = 42/10)": 1/root(b1=mp.mpf(42)/10), "kein Higgs (b1 = 4)": 1/root(b1=mp.mpf(4)), "vier Familien (b1 = 163/30)": 1/root(b1=mp.mpf(163)/30), "zwei Familien (b1 = 83/30)": 1/root(b1=mp.mpf(83)/30)}
R["neighbours"] = {k: mp.nstr(v, 8) for k, v in nb.items()}
check("N6 other field contents give 135.9 (two Higgs), 138.2 (no Higgs), 124.8 (four families), 156.1 (two families)", abs(nb["zwei Higgs (b1 = 42/10)"]-135.9) < 0.2 and abs(nb["kein Higgs (b1 = 4)"]-138.2) < 0.2 and abs(nb["vier Familien (b1 = 163/30)"]-124.5) < 0.5 and abs(nb["zwei Familien (b1 = 83/30)"]-156.5) < 0.5, {k: mp.nstr(v, 6) for k, v in nb.items()})
# N7 the first rung without any chosen correction
check("N7 the cubic with phi = 1/(6 pi) alone is 252 ppm from CODATA; the weight 48 c3^4 brings it to 3.7 ppm", abs(abs(ppm[2])-252) < 1 and abs(abs(ppm[3])-3.67) < 0.05, (round(ppm[2], 2), round(ppm[3], 2)))
# N8 the undamped weight N c3^4 that reproduces CODATA, and the frozen root
Nc = mp.findroot(lambda N: 1/root(phi=lambda a: pt+N*c3**4)-A, 47.3); Nt = mp.findroot(lambda N: 1/root(phi=lambda a: pt+N*c3**4)-1/al, 47.3)
sN = S/((1/root(phi=lambda a: pt+48*c3**4)-1/root(phi=lambda a: pt+47*c3**4)))
R["N_codata"] = [mp.nstr(Nc, 9), float(abs(sN))]; R["N_frozen"] = mp.nstr(Nt, 9)
check("N8 without damping CODATA needs the weight 47.31150(3) c3^4 instead of 48 c3^4; damping and resummation together supply 48 -> 47.31155", abs(Nc-mp.mpf("47.311498")) < 2e-6 and abs(Nt-mp.mpf("47.311553")) < 2e-6 and 2.5e-5 < abs(sN) < 3.2e-5, (mp.nstr(Nc, 8), mp.nstr(Nt, 8), mp.nstr(abs(sN), 3)))
R["checks_passed"] = not FAILS
Path(__file__).with_name("results_normal.json").write_text(json.dumps(R, indent=1, ensure_ascii=False) + "\n")
print("RESULT", "ALL PASS" if not FAILS else "FAIL " + ", ".join(FAILS)); sys.exit(1 if FAILS else 0)
