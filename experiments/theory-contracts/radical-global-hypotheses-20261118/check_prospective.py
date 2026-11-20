#!/usr/bin/env python3
"""RADICALPROSPECT20261118: Prospektive Testgroesse fuer die Alpha Gleichung.
Die Formel ist eingefroren (verification/v3_em_alpha.py). Vergleich mit den drei einzelnen
Bestimmungen statt mit dem CODATA Mittel; diese stehen in keiner Quelle des Repositorys.
Messwerte (am 18.11.2026 per Websuche bestaetigt):
  Rb Rueckstoss, Morel et al. 2020:   137.035999206(11)
  Cs Rueckstoss, Parker et al. 2018:  137.035999046(27)
  Elektron g-2 + SM, Fan et al. 2023: 137.035999166(15), g/2 = 1.00115965218059(13)
Reproduktion: python3 -B check_prospective.py"""
import json, sys
from pathlib import Path
import mpmath as mp
mp.mp.dps = 30
R = {}; FAILS = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + (" | " + str(info) if info != "" else ""))
    if not ok: FAILS.append(name)
pi = mp.pi; c3 = 1/(8*pi); pt = 1/(6*pi); d = 48*c3**4; b1 = mp.mpf(41)/10
def ainv(f):
    return 1/mp.findroot(lambda a: a**3-2*c3**3*a**2-8*b1*c3**6*mp.log(1/(pt+f(d*mp.e**(-2*a)))), mp.mpf(1)/137)
T = ainv(lambda x: x*(1-x)**(-mp.mpf(5)/4)); R["tfpt_alpha_inv"] = mp.nstr(T, 16)
M = {"Rb 2020": ("137.035999206", "1.1e-8"), "Cs 2018": ("137.035999046", "2.7e-8"), "g-2 2023": ("137.035999166", "1.5e-8"), "CODATA 2022": ("137.035999177", "2.1e-8")}
Z = {k: float((T-mp.mpf(v))/mp.mpf(s)) for k, (v, s) in M.items()}; R["pulls"] = Z
for k, z in Z.items(): print(f"   {k:12s} {M[k][0]}({M[k][1]})  Zug {z:+.2f} sigma")
check("P1 frozen root 137.0359992168407", abs(T-mp.mpf("137.0359992168407")) < mp.mpf("1e-12"), mp.nstr(T, 16))
check("P2 pulls: Rb +0.98, g-2 +3.39, Cs +6.33, CODATA +1.90 sigma", abs(Z["Rb 2020"]-0.98) < 0.02 and abs(Z["g-2 2023"]-3.39) < 0.02 and abs(Z["Cs 2018"]-6.33) < 0.02 and abs(Z["CODATA 2022"]-1.90) < 0.02, {k: round(v, 2) for k, v in Z.items()})
# P3 the three inputs disagree among themselves
rb, cs, ge = (mp.mpf(M[k][0]) for k in ("Rb 2020", "Cs 2018", "g-2 2023")); s = lambda a, b: mp.sqrt(mp.mpf(M[a][1])**2+mp.mpf(M[b][1])**2)
ten = {"Rb-Cs": float((rb-cs)/s("Rb 2020", "Cs 2018")), "Rb-g2": float((rb-ge)/s("Rb 2020", "g-2 2023")), "g2-Cs": float((ge-cs)/s("g-2 2023", "Cs 2018"))}; R["tensions"] = ten
check("P3 the measurements disagree: Rb-Cs 5.5, Rb-g2 2.2, g2-Cs 3.9 sigma; the CODATA mean hides this", abs(ten["Rb-Cs"]-5.49) < 0.05 and abs(ten["Rb-g2"]-2.15) < 0.05 and abs(ten["g2-Cs"]-3.89) < 0.05, {k: round(v, 2) for k, v in ten.items()})
# P4 weighted mean of the three and its chi^2; TFPT against the unexpanded mean
w = {k: 1/mp.mpf(M[k][1])**2 for k in ("Rb 2020", "Cs 2018", "g-2 2023")}; mean = sum(w[k]*mp.mpf(M[k][0]) for k in w)/sum(w.values()); sm = 1/mp.sqrt(sum(w.values()))
chi2 = float(sum(w[k]*(mp.mpf(M[k][0])-mean)**2 for k in w)); R["weighted_mean"] = [mp.nstr(mean, 13), float(sm), chi2, float((T-mean)/sm)]
sc = mp.sqrt(chi2/2); R["weighted_mean"].append(float((T-mean)/(sm*sc)))
check("P4 weighted mean 137.0359991778(84) with chi^2 = 31 for 2 degrees of freedom (inconsistent inputs); TFPT lies +4.6 sigma above it, +1.2 sigma after scaling the error by sqrt(chi^2/2)", abs(mean-mp.mpf("137.0359991778")) < 1e-10 and 30 < chi2 < 32 and abs(float((T-mean)/sm)-4.63) < 0.05 and abs(float((T-mean)/(sm*sc))-1.18) < 0.05, (mp.nstr(mean, 13), mp.nstr(sm, 2), round(chi2, 1), round(float((T-mean)/sm), 2), round(float((T-mean)/(sm*sc)), 2)))
# P5 electron g-2 implied by the TFPT alpha (leading sensitivity da/a = d alpha/alpha)
a_exp = mp.mpf("0.00115965218059"); sa = mp.mpf("1.3e-13"); a_pred = a_exp*(ge/T); R["a_e_pred"] = [mp.nstr(a_pred, 15), float((a_exp-a_pred)/sa)]
check("P5 with standard QED the TFPT alpha implies a_e = 0.00115965218016, 4.3e-13 below the measured value (3.3 sigma of the measurement alone)", abs(a_pred-mp.mpf("0.00115965218016")) < 2e-14 and 3.0 < float((a_exp-a_pred)/sa) < 3.6, (mp.nstr(a_pred, 15), round(float((a_exp-a_pred)/sa), 2)))
# P6 decision power: a future determination with uncertainty 3e-9
fs = mp.mpf("3e-9"); dec = {k: float(abs(T-mp.mpf(M[k][0]))/fs) for k in ("Rb 2020", "Cs 2018", "g-2 2023")}; R["future_3e-9"] = dec
check("P6 a future value with uncertainty 3e-9 separates TFPT from today's central values by 3.6 (Rb), 16.9 (g-2) and 56.9 (Cs) sigma", abs(dec["Rb 2020"]-3.6) < 0.1 and abs(dec["g-2 2023"]-16.9) < 0.2 and abs(dec["Cs 2018"]-56.9) < 0.3, {k: round(v, 1) for k, v in dec.items()})
# P7 the alternatives that were on the table give distinguishable numbers
alt = {"Koeffizient 1 (Fassung 2.4)": ainv(lambda x: x+x*x), "ohne Zweidefektterm": ainv(lambda x: x), "CODATA optimal 1,2404": ainv(lambda x: x+mp.mpf("1.2404297")*x*x)}
R["alternatives"] = {k: mp.nstr(v, 15) for k, v in alt.items()}
check("P7 the CODATA-optimal coefficient 1.2404 would give 137.035999177; 5/4 gives a value 4.0e-8 higher, so the prediction is not the fitted mean", abs(alt["CODATA optimal 1,2404"]-mp.mpf("137.035999177")) < 1e-11 and abs(T-alt["CODATA optimal 1,2404"]-mp.mpf("3.98e-8")) < 2e-9, mp.nstr(T-alt["CODATA optimal 1,2404"], 4))
# P8 other readouts use the undamped seed; the relative difference to the alpha seam value
al = 1/T; q = d*mp.e**(-2*al); seam = pt+q*(1-q)**(-mp.mpf(5)/4); rel = (pt+d-seam)/(pt+d); R["seed_vs_seam_rel"] = float(rel)
check("P8 the seed phi0 used by the other readouts and the alpha seam value differ by 3.2e-5 relative; no other readout is compared at that precision", abs(rel-mp.mpf("3.24e-5")) < 2e-6, mp.nstr(rel, 4))
R["checks_passed"] = not FAILS
Path(__file__).with_name("results_prospective.json").write_text(json.dumps(R, indent=1, ensure_ascii=False) + "\n")
print("RESULT", "ALL PASS" if not FAILS else "FAIL " + ", ".join(FAILS)); sys.exit(1 if FAILS else 0)
