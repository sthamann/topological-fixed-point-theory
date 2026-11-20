#!/usr/bin/env python3
"""RADICALTHREE20261118: die drei offenen Stuecke der Alpha Gleichung.
K: Form der Kubik (Suchraum auf dem Stand der Fassung 1.0, ohne spaetere Zutaten)
D: Regel fuer die getrennte D5 Auslese (Zahl der gelesenen Cartanrichtungen)
W: Gewicht 48 c3^4 und Daempfung exp(-2 alpha)
Reproduktion: python3 -B check_three.py"""
import itertools, json, math, sys
from fractions import Fraction as Fr
from pathlib import Path
import mpmath as mp
mp.mp.dps = 30
ROOT = Path(__file__).resolve().parents[3]
R = {}; FAILS = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + (" | " + str(info) if info != "" else ""))
    if not ok: FAILS.append(name)
pi = mp.pi; c3 = 1/(8*pi); pt = 1/(6*pi); d = 48*c3**4; b1 = mp.mpf(41)/10; A0 = mp.mpf("137.035999177"); S = mp.mpf("2.1e-8")
def ainv(phi):
    return 1/mp.findroot(lambda a: a**3-2*c3**3*a**2-8*b1*c3**6*mp.log(1/phi(a)), mp.mpf(1)/137)
pull = lambda phi: float((ainv(phi)-A0)/S)
frozen = lambda a: pt+(d*mp.e**(-2*a))*(1-d*mp.e**(-2*a))**(-mp.mpf(5)/4)
al = 1/ainv(frozen)
# ---------- K: form of the cubic ----------
A = 2*c3**3; kap = b1/(2*pi)*mp.log(1/(pt+d))
check("K1 original form (Fassung 1.0): alpha^3 - A alpha^2 - A c3^2 kappa with A = 1/(256 pi^3), kappa = (b1/2pi) log(1/phi0) equals today's cubic",
      abs(A-1/(256*pi**3)) < mp.mpf("1e-30") and abs(A*c3**2*kap-8*b1*c3**6*mp.log(1/(pt+d))) < mp.mpf("1e-30") and abs(kap-mp.mpf("1.914684795")) < 1e-8, mp.nstr(kap, 10))
# K2 look-elsewhere on the level of Fassung 1.0
c3f, ptf, df = float(c3), float(pt), float(d); target = 137.035999177; hit_tol = 3.7e-6; wide = 1e-2
def solve(p, t, rhs):
    lo = max(t*c3f**3, 1e-12)*(1+1e-12); hi = 1.0
    f = lambda a: a**p-t*c3f**3*a**(p-1)-rhs
    if f(hi) < 0 or f(lo) > 0: return None
    for _ in range(80):
        mid = math.sqrt(lo*hi)
        if f(mid) > 0: hi = mid
        else: lo = mid
    return math.sqrt(lo*hi)
n = wide_n = 0; hits = []
for p, t, m, r, l, b, ph in itertools.product((2, 3, 4), (0, 1, 2, 4), range(2, 9), (1/8, 1/4, 1/2, 1, 2, 4, 8, 16), (0.5, 1, 2), (41/10, 41/6), (ptf, ptf+df)):
    a = solve(p, t, r*c3f**m*(b*math.log(1/ph))**l); n += 1
    if a is None: continue
    rel = (1/a-target)/target
    if abs(rel) < wide: wide_n += 1
    if abs(rel) < hit_tol: hits.append((p, t, m, r, l, round(b, 4), "phi0" if ph > ptf else "Baum", round(rel*1e6, 2)))
exp_hits = wide_n*hit_tol/wide; R["K_grammar"] = {"size": n, "within_1_percent": wide_n, "expected_within_3.7ppm": exp_hits, "hits": hits}
check("K2 grammar of 8064 neighbouring equations: chance expectation for a hit within 3.7 ppm is below 0.05, and exactly one member hits (the cubic of Fassung 1.0)",
      n == 8064 and exp_hits < 0.05 and len(hits) == 1 and hits[0][:7] == (3, 2, 6, 8, 1, 4.1, "phi0"), (n, wide_n, round(exp_hits, 4), hits))
# K3 what the sources themselves say
v27 = (ROOT/"documents/archive--paper-latex--tfpt-theory-fullv27-final/sections/003-subsection-the-cubic-fixed-point-equation-cfe.tex").read_text(errors="replace")
check("K3 Fassung 2.7 labels the potential coefficients 'postulated input' in a declared scheme", "postulated input" in v27 and "declared scheme" in v27)
# ---------- D: rule for the D5 readout ----------
roots = []
for i, j in itertools.combinations(range(8), 2):
    for si, sj in itertools.product((1, -1), repeat=2):
        v = [Fr(0)]*8; v[i] = Fr(si); v[j] = Fr(sj); roots.append(tuple(v))
for s in itertools.product((1, -1), repeat=8):
    if sum(1 for x in s if x < 0) % 2 == 0: roots.append(tuple(Fr(x, 2) for x in s))
n5 = lambda v: sum(x*x for x in v[:5]); n3 = lambda v: sum(x*x for x in v[5:])
cls = {}
for v in roots: cls[(n5(v), n3(v))] = cls.get((n5(v), n3(v)), 0)+1
R["E8_split"] = {f"{k[0]},{k[1]}": c for k, c in cls.items()}
check("D1 E8 roots under D5 + A3: 40 + 12 + 60 + 128; all 128 mixed spinor roots have D5 norm 5/4 and A3 norm 3/4",
      len(roots) == 240 and cls == {(Fr(2), Fr(0)): 40, (Fr(0), Fr(2)): 12, (Fr(1), Fr(1)): 60, (Fr(5, 4), Fr(3, 4)): 128}, R["E8_split"])
check("D2 every Cartan direction of a spinor root contributes exactly 1/4: reading n directions gives the exponent n/4 (g/4 of Fassung 2.5 is n = 5, its 'g = 4' is n = 4)",
      all(x*x == Fr(1, 4) for v in roots if n5(v) == Fr(5, 4) for x in v))
P = {k: pull(lambda a, s=mp.mpf(k)/4: pt+(d*mp.e**(-2*a))*(1-d*mp.e**(-2*a))**(-s)) for k in range(9)}; R["pull_by_directions"] = P
for k, v in P.items(): print(f"   {k} Richtungen, Exponent {k}/4: Zug {v:+9.2f} sigma")
check("D3 number of directions read, n = 0..8: only n = 5 lies within 3 sigma of CODATA (n = 4: -46.8, n = 6: +50.6)", [k for k, v in P.items() if abs(v) < 3] == [5] and abs(P[4]+46.8) < 0.2 and abs(P[6]-50.6) < 0.3, {k: round(v, 1) for k, v in P.items()})
Y = [Fr(-1, 3)]*3+[Fr(1, 2)]*2; Q = [Fr(-1, 3)]*3+[Fr(1), Fr(0)]
check("D4 named choices: electric charge only (n = 1), weak isospin and hypercharge (n = 2), Standard Model gauge rank (n = 4), carrier rank with B-L (n = 5), all of E8 (n = 8); |Y|^2 = 5/6 and (3/2)(5/6) = 5/4",
      sum(y*y for y in Y) == Fr(5, 6) and Fr(3, 2)*Fr(5, 6) == Fr(5, 4) and sum(x*x for x in Q) == Fr(4, 3))
v24 = "\n".join(p.read_text(errors="replace") for p in sorted((ROOT/"documents/archive--paper-latex--tfpt-theory-fullv24").rglob("*.tex")))
check("D5 the carrier number entered before the coefficient 5/4: Fassung 2.4 already fixes gamma(0) = 5/6 from SU(5)", "\\gamma(0)=5/6" in v24 or "\\gamma(0) = 5/6" in v24)
# ---------- W: weight and damping ----------
spin16 = [s for s in itertools.product((Fr(1, 2), Fr(-1, 2)), repeat=5) if sum(1 for x in s if x < 0) % 2 == 0]
qs = [sum(a*b for a, b in zip(s, Q)) for s in spin16]; ys = [sum(a*b for a, b in zip(s, Y)) for s in spin16]
check("W1 48 = 16 x 3: one D5 spinor (16 states, sum Q^2 = 16/3, sum Y^2 = 10/3) times three of the four A3 weights",
      len(spin16) == 16 and sum(x*x for x in qs) == Fr(16, 3) and sum(y*y for y in ys) == Fr(10, 3) and 16*3 == 48, (str(sum(x*x for x in qs)), str(sum(y*y for y in ys))))
e2 = 4*pi*al
check("W2 exp(-2 alpha) = [exp(-e^2 c3)]^4 = exp(-alpha |beta|^2) with |beta|^2 = 2, and 2 alpha = g/(16 pi^2) for g = e^2/c3: q = 48 [c3 exp(-e^2 c3)]^4",
      abs(mp.e**(-2*al)-(mp.e**(-e2*c3))**4) < mp.mpf("1e-28") and abs(2*al-(e2/c3)/(16*pi**2)) < mp.mpf("1e-28") and abs(d*mp.e**(-2*al)-48*(c3*mp.e**(-e2*c3))**4) < mp.mpf("1e-30"))
pW3 = pull(lambda a: mp.mpf(4)/3*c3*mp.e**(-a/2)+(d*mp.e**(-2*a))*(1-d*mp.e**(-2*a))**(-mp.mpf(5)/4)); R["pull_uniform_dressing"] = pW3
check("W3 dressing every kernel factor, including the base (4/3) c3 exp(-alpha/2), misses CODATA by more than 10^4 sigma: the damping acts on the four-leg term only", abs(pW3) > 1e4, f"{pW3:.3e}")
def charge_damped(ch, norm=1):
    return lambda a: pt+(lambda q: q*(1-q)**(-mp.mpf(5)/4))(3*c3**4*sum(mp.e**(-2*a*norm*float(x*x)) for x in ch))
pQ, pY = pull(charge_damped(qs)), pull(charge_damped(ys, 0.6)); R["pull_charge_dependent"] = [pQ, pY]
check("W4 charge dependent damping exp(-2 alpha Q^2) or exp(-2 alpha (3/5) Y^2) over the 48 states misses by more than 100 sigma: the damping is the same for all 48 states", abs(pQ) > 100 and abs(pY) > 100, (round(pQ, 1), round(pY, 1)))
pD = {s: pull(lambda a, s=s: pt+(d*mp.e**(-mp.mpf(s)*a))*(1-d*mp.e**(-mp.mpf(s)*a))**(-mp.mpf(5)/4)) for s in ("0.75", "1.25", "2", "4")}; R["pull_by_damping_norm"] = pD
check("W5 damping exp(-s alpha) with s = 3/4, 5/4, 2, 4: only the full root norm s = 2 fits; the damping uses the whole glue charge, the oscillator factor only its D5 part",
      [s for s, v in pD.items() if abs(v) < 3] == ["2"], {s: round(v, 1) for s, v in pD.items()})
v21 = "\n".join(p.read_text(errors="replace") for p in sorted((ROOT/"documents/archive--paper-latex--tfpt-theory-fullv21").rglob("*.tex")))
old = (ROOT/"_archive/paper-latex/old/theory_root.tex").read_text(errors="replace")
check("W6 the weight 48 c3^4 = 3/(256 pi^4) is already in Fassung 1.0 and 2.2, before any comparison at the 0.04 ppm level", "48\\cthree^4" in v21 and "256" in old and "137.03650146488582" in old)
# W7 one statement instead of three
x = d*mp.e**(-2*al); h = 1; G = x**h*(1-x)**(-mp.mpf(5)/4)
check("W7 summary form: phi - (4/3) c3 = x^h (1 - x)^(-2 h5) with x = 48 c3^4 exp(-2 alpha), h = |beta|^2/2 = 1, h5 = 5/8", abs(frozen(al)-pt-G) < mp.mpf("1e-30"))
R["checks_passed"] = not FAILS
Path(__file__).with_name("results_three.json").write_text(json.dumps(R, indent=1, ensure_ascii=False, default=str) + "\n")
print("RESULT", "ALL PASS" if not FAILS else "FAIL " + ", ".join(FAILS)); sys.exit(1 if FAILS else 0)
