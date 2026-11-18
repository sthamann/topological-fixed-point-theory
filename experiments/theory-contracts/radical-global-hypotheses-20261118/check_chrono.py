#!/usr/bin/env python3
"""RADICALCHRONO20261118: Chronologie des Koeffizienten 5/4 in der Alpha Rueckkopplung.
Quelle: datierte Fassungen unter documents/archive--paper-latex--*.
Reproduktion: python3 -B check_chrono.py"""
import json, math, re, sys
from pathlib import Path
import mpmath as mp
mp.mp.dps = 30
ROOT = Path(__file__).resolve().parents[3] / "documents"
R = {}; FAILS = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + (" | " + str(info) if info != "" else ""))
    if not ok: FAILS.append(name)
def text(d):
    return "\n".join(p.read_text(errors="replace") for p in sorted((ROOT / d).rglob("*.tex")))
V = {k: text("archive--paper-latex--" + k) for k in ("tfpt-theory-fullv21", "tfpt-theory-fullv24", "tfpt-theory-fullv25", "tfpt-31")}
five4 = re.compile(r"frac54|frac\{5\}\{4\}|g/4|\\frac\{g\}\{4\}")
# K1 dates and first appearance
check("K1 v2.2 (21 Jan 2026) and v2.4 (22 Jan 2026) contain no 5/4 or g/4; v2.5 (27 Jan 2026) does",
      "21 January 2026" in V["tfpt-theory-fullv21"] and "22 January 2026" in V["tfpt-theory-fullv24"] and "27 January 2026" in V["tfpt-theory-fullv25"]
      and not five4.search(V["tfpt-theory-fullv21"]) and not five4.search(V["tfpt-theory-fullv24"]) and bool(five4.search(V["tfpt-theory-fullv25"])))
# K2 the measured value was in the text before 5/4 appeared
check("K2 CODATA 137.035999177 is cited in v2.2 and v2.4, before the 5/4 coefficient", all("137.035999177" in V[k] for k in ("tfpt-theory-fullv21", "tfpt-theory-fullv24")))
# K3 v2.4 two-defect term with coefficient 1, labelled diagnostic
check("K3 v2.4 carries a two-defect term labelled diagnostic", "two-defect" in V["tfpt-theory-fullv24"].lower() and "diagnostic" in V["tfpt-theory-fullv24"].lower())
# K4 v2.5 computed the delta_2 needed for an exact CODATA match and relabelled coefficient 1 as legacy g=4
m = re.search(r"TFPTdeltaTwoOptimalNum\}\{([0-9.eE+-]+)\}", V["tfpt-theory-fullv25"])
check("K4 v2.5 defines the CODATA-optimal delta_2 ('targeted match') and calls delta_2 = delta_top^2 the legacy g=4 template",
      m is not None and "targeted match" in V["tfpt-theory-fullv25"] and "legacy ordered-occupancy template" in V["tfpt-theory-fullv25"])
pi = mp.pi; c3 = 1/(8*pi); pt = 1/(6*pi); d = 48*c3**4; b1 = mp.mpf(41)/10
A = mp.mpf("137.035999177"); S = mp.mpf("2.1e-8")
def ainv(phi):
    return 1/mp.findroot(lambda a: a**3-2*c3**3*a**2-8*b1*c3**6*mp.log(1/phi(a)), mp.mpf(1)/137)
def trunc(c): return ainv(lambda a: pt+d*mp.e**(-2*a)+c*d*d*mp.e**(-4*a))
z = lambda x: float((x-A)/S)
# K5 optimal coefficient in units of delta_top^2
copt_doc = mp.mpf(m.group(1))/d**2
copt = mp.findroot(lambda c: trunc(c)-A, 1.2)
sig_c = S/((trunc(1.25)-trunc(1))/mp.mpf("0.25"))
R["c_opt_document"] = float(copt_doc); R["c_opt_recomputed"] = float(copt); R["sigma_c"] = float(sig_c)
check("K5 CODATA-optimal coefficient is 1.2404 +- 0.0051 delta_top^2 (document and recomputation agree); 5/4 lies 1.9 sigma above",
      abs(copt_doc-copt) < 1e-6 and abs(copt-mp.mpf("1.2404")) < 1e-3 and 1.8 < float((mp.mpf(5)/4-copt)/sig_c) < 1.95, (mp.nstr(copt_doc, 8), mp.nstr(copt, 8), mp.nstr(sig_c, 4)))
# K6 pulls of the two versions
z1, z54 = z(trunc(1)), z(trunc(mp.mpf(5)/4)); R["pull_coeff_1"] = z1; R["pull_coeff_5/4"] = z54
check("K6 coefficient 1 (v2.4) misses CODATA by -46.9 sigma, coefficient 5/4 (v2.5) by +1.86 sigma", abs(z1+46.85) < 0.1 and abs(z54-1.865) < 0.01, (round(z1, 2), round(z54, 3)))
# K7 look-elsewhere over simple rational coefficients
def scan(maxden):
    n = 0; hits = []
    for b in range(1, maxden+1):
        for a in range(1, 3*b+1):
            if math.gcd(a, b) == 1:
                n += 1
                if abs(z(trunc(mp.mpf(a)/b))) < 2: hits.append(f"{a}/{b}")
    return n, hits
n4, h4 = scan(4); n12, h12 = scan(12); width = float(4*sig_c)
R["scan"] = {"den<=4": [n4, h4, n4*width/3], "den<=12": [n12, h12, n12*width/3]}
check("K7 rational coefficients in (0,3]: denominators <= 4 give 18 candidates, chance expectation 0.12, hit 5/4; denominators <= 12 give 138, expectation 0.94, hit only 5/4",
      n4 == 18 and h4 == ["5/4"] and n12 == 138 and h12 == ["5/4"] and abs(n4*width/3-0.123) < 0.01 and abs(n12*width/3-0.944) < 0.02, R["scan"])
# K8 resummed form appears later and changes little
res = ainv(lambda a: pt+d*mp.e**(-2*a)*(1-d*mp.e**(-2*a))**(-mp.mpf(5)/4))
R["resummed"] = [mp.nstr(res, 16), z(res)]
check("K8 the resummed form (1-q)^(-5/4) is absent in v2.5 and present in tfpt-31; it moves the pull only from 1.865 to 1.897",
      "^{-5/4}" not in V["tfpt-theory-fullv25"] and "^{-5/4}" in V["tfpt-31"] and abs(z(res)-1.897) < 0.01, R["resummed"])
# K9 natural neighbours fail
nb = {s: z(trunc(mp.mpf(a)/b)) for s, a, b in (("6/5", 6, 5), ("4/3", 4, 3), ("9/8", 9, 8), ("3/2", 3, 2))}
R["neighbours"] = nb
check("K9 neighbouring coefficients 6/5, 4/3, 9/8, 3/2 all miss by more than 7 sigma", min(abs(v) for v in nb.values()) > 7, {k: round(v, 1) for k, v in nb.items()})
R["checks_passed"] = not FAILS
Path(__file__).with_name("results_chrono.json").write_text(json.dumps(R, indent=1, ensure_ascii=False) + "\n")
print("RESULT", "ALL PASS" if not FAILS else "FAIL " + ", ".join(FAILS)); sys.exit(1 if FAILS else 0)
