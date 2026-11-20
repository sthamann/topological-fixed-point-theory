"""Reihenprüfer: viele Ursprungshypothesen gleichzeitig gegen dieselben festen Zielgrößen.

Jede Hypothese ist ein Eintrag {id, fam, text, wahl, pred}.  'pred' sagt Werte für Zielgrößen voraus:
  c3     Normierung (Ziel 1/(8 pi))            k      Dämpfung exp(-k alpha) (Ziel 2)
  expo   Exponent in (1-q)^(-expo) (Ziel 5/4)  ratio  phi_tree / c3 (Ziel 4/3)
  n48    delta_top / c3^4 (Ziel 48)            clock  Ordnung der Flächenuhr (Ziel 4)
  five   Ordnung der zweiten Uhr (Ziel 5)      lam    Doppelverhältnis der vier Marken (Ziel 1/2)
'wahl' zählt die freien Entscheidungen, die die Hypothese braucht.  Bewertet wird automatisch:
exakter Treffer je Zielgröße und, wo die Alpha Gleichung berührt ist, der Abstand in sigma.
Neue Hypothesen: eine Funktion fam_*() ergänzen, die Einträge liefert.  Aufruf: python3 -B sweep.py
"""
import json, itertools, sys
from fractions import Fraction as Fr
from mpmath import mp, mpf, mpc, pi, exp, log, sqrt, sin, cos, coth, gamma, binomial, findroot, ellipk, expj, nsum, inf, fabs
mp.dps = 30
C3 = 1 / (8 * pi); CODATA, SIG = mpf("137.035999177"), mpf("0.000000021")
TARGET = {"A": mpf(2), "B": mpf(164) / 5, "c3": C3, "k": mpf(2), "expo": mpf(5) / 4, "ratio": mpf(4) / 3, "n48": mpf(48), "clock": 4, "five": 5, "lam": mpf(1) / 2}

def alpha_z(c3=C3, k=2, expo=mpf(5) / 4, ratio=mpf(4) / 3, n48=48, m=4, A=2, B=mpf(164) / 5):
    def F(a):
        q = n48 * c3 ** m * exp(-k * a); phi = ratio * c3 + q * (1 - q) ** (-expo)
        return a ** 3 - A * c3 ** 3 * a ** 2 - B * c3 ** 6 * log(1 / phi)
    try:
        a = findroot(F, 2 * float(c3) ** 2 * 0.6 * 0 + float(c3) ** 2 * 4.6)
        return float((1 / a - CODATA) / SIG)
    except Exception:
        return None

H = []
def add(id_, fam, text, wahl, **pred): H.append({"id": id_, "fam": fam, "text": text, "wahl": wahl, "pred": pred})

# ---------- Familie A: Schnittgeometrie und Nullmoden -> c3 ----------
def f_torus(y):   # log(sqrt(y)|eta|^2) auf tau = i y
    q = exp(-2 * pi * y); p = mpf(1)
    for n in range(1, 60): p *= (1 - q ** n)
    return log(sqrt(y) * (exp(-pi * y / 12) * p) ** 2)
def fam_cut():
    tau_of = lambda lam: ellipk(1 - lam) / ellipk(lam)
    shapes = {"Quadrat (lam=1/2)": mpf(1), "Rechteck lam=1/3": tau_of(mpf(1) / 3), "Rechteck lam=1/4": tau_of(mpf(1) / 4),
              "Rechteck 2:1": mpf(2), "Sechseck": sqrt(3) / 2}
    for name, y in shapes.items():
        for I, why in ((2, "volle Nullmode"), (4, "chirale Nullmode"), (8, "chiral, zwei Blätter")):
            add(f"A:{name}:I={I}", "A Schnitt", f"Torus {name}, {why}: c3 = 1/(I*2pi*Im tau)", 2, c3=1 / (I * 2 * pi * y))
    for n in range(1, 8):   # runde Sphären: Dirac Gewicht des tiefsten positiven Kanals -> k = 4 pi w
        vol = 2 * pi ** (mpf(n + 1) / 2) / gamma(mpf(n + 1) / 2)
        w = 2 ** (n // 2) / vol / (mpf(n) / 2)
        add(f"A:S^{n}:Dirac+", "A Schnitt", f"runde S^{n}, tiefster positiver Dirac Kanal: k = 4 pi w", 1, k=4 * pi * w,
            clock=(4 if n % 2 == 0 else 8))
    add("A:Weyl2D", "A Schnitt", "beliebiger 2D Schnitt, Weylgewicht 1/(2pi) je Spektraleinheit", 1, k=mpf(2))
    for n in (3, 4, 5): add(f"A:Weyl{n}D", "A Schnitt", f"{n}D Schnitt: Gewicht wächst mit K^{n-2}, keine Konstante", 1, k=None)

# ---------- Familie B: Zustand -> Dämpfung k ----------
def fam_state():
    add("B:skalar1", "B Zustand", "ein skalarer Kugelkanal", 1, k=mpf(1))
    add("B:Dirac+", "B Zustand", "ein positiver Dirac Kanal (Calderon Seite)", 1, k=mpf(2))
    add("B:Dirac+-", "B Zustand", "beide Dirac Vorzeichen", 1, k=mpf(4))
    add("B:Vakuum_h1", "B Zustand", "Fourier-selbstdualer Grundzustand, beide Quadraturen", 2, k=mpf(2))
    add("B:Vakuum_hbar1", "B Zustand", "Grundzustand mit hbar = 1 statt h = 1", 2, k=4 * pi * 1)
    for w in ("0.5", "1", "2", "3.141592653589793"):
        add(f"B:thermisch_w={w[:4]}", "B Zustand", f"Quantenoszillator omega = {w[:4]} bei beta = 2 pi", 2, k=4 * pi * coth(pi * mpf(w)) / (2 * mpf(w)))
    add("B:klassisch_T", "B Zustand", "klassische Nullmode bei T = 1/(2 pi), Steifigkeit 1", 2, k=mpf(2))
    for T, nm in ((1 / pi, "1/pi"), (1 / (4 * pi), "1/(4pi)")):
        add(f"B:klassisch_T={nm}", "B Zustand", f"klassische Nullmode bei T = {nm}", 2, k=4 * pi * T)
    for N in (1, 8, 240, 1920):
        add(f"B:E8_Keime_N={N}", "B Zustand", f"{N} unabhängige E8 Gitterkeime", 2, k="shift", shift=139.75 / N)

# ---------- Familie C: Gitterklasse -> Exponent ----------
def fam_expo():
    for n in range(2, 9): add(f"C:D{n}_Spinor", "C Exponent", f"Spinorklasse von D{n}: Norm n/4", 1, expo=mpf(n) / 4)
    for n in range(1, 8):
        for kcl in range(1, (n + 1) // 2 + 1):
            add(f"C:A{n}_Klasse{kcl}", "C Exponent", f"Klasse {kcl} von A{n}: Norm k(n+1-k)/(n+1)", 2, expo=mpf(kcl * (n + 1 - kcl)) / (n + 1))
    add("C:E6", "C Exponent", "Klebeklasse von E6: 4/3", 1, expo=mpf(4) / 3); add("C:E7", "C Exponent", "Klebeklasse von E7: 3/2", 1, expo=mpf(3) / 2)

# ---------- Familie D: ADE Schnitte -> zweite Uhr ----------
def fam_ade():
    w = lambda num, den: 1 / (2 * sin(pi * mpf(num) / den / 2) ** 2)
    def order(classes):
        n = 1 + sum(m for m, _ in classes); v = sum(m * w(*ph) for m, ph in classes) / n
        return (Fr(str(v)).limit_denominator(10 ** 6) / 4).denominator
    for p in range(2, 31): add(f"D:A{p-1}(Z{p})", "D Singularität", f"Schnitt S^3/Z_{p}", 1, five=order([(1, (2 * j, p)) for j in range(1, p)]))
    for n in range(2, 16): add(f"D:D{n+2}(Dic{n})", "D Singularität", f"Schnitt S^3/Dic_{n}", 1, five=order([(1, (j, n)) for j in range(1, 2 * n)] + [(2 * n, (1, 2))]))
    E = {"E6": [(1, (1, 1)), (6, (1, 2)), (8, (1, 3)), (8, (2, 3))], "E7": [(1, (1, 1)), (18, (1, 2)), (6, (1, 4)), (6, (3, 4)), (8, (1, 3)), (8, (2, 3))],
         "E8": [(1, (1, 1)), (30, (1, 2)), (20, (1, 3)), (20, (2, 3)), (12, (1, 5)), (12, (3, 5)), (12, (2, 5)), (12, (4, 5))]}
    for nm, cl in E.items(): add(f"D:{nm}", "D Singularität", f"Schnitt der {nm} Singularität", 1, five=order(cl))

# ---------- Familie E: Auswahl der vier Marken -> lam ----------
def fam_marks():
    def stab(lam):   # Ordnung der Möbiusgruppe, die vier Marken mit Doppelverhältnis lam erhält; enthält sie Ordnung 4?
        orb = {1 / lam, 1 - lam, lam / (lam - 1)}
        four = any(abs(x - lam) < mpf(10) ** -20 for x in orb)
        three = abs(1 / (1 - lam) - lam) < mpf(10) ** -20
        return (24 if four and three else 8 if four else 12 if three else 4), four
    ys = [mpf(k) / 200 for k in range(40, 500)]; best_rect = ys[max(range(len(ys)), key=lambda i: f_torus(ys[i]))]
    add("E:Spiegelung_allein", "E Marken", "Spiegelung erhält die Markenmenge: j reell, Torus rechteckig ODER rhombisch; lam bleibt frei", 0, lam=None)
    add("E:Marken_auf_Naht", "E Marken", "Marken = vier Ecken auf der Naht (Doppel eines flachen Rechtecks): lam reell in (0,1), Seitenverhältnis frei", 1, lam=None)
    add("E:Klein_transitiv", "E Marken", "transitive Wirkung durch die Kleinsche Vierergruppe: gilt für jedes lam, wählt nichts", 1, lam=None,
        clock=2)
    add("E:zyklisch_Z4", "E Marken", "transitive Wirkung enthält Ordnung 4 (mu_4): erzwingt lam in {-1, 1/2, 2}", 1,
        lam=[l for l in (mpf(1) / 2, mpf(1) / 3, mpf(1) / 4) if stab(l)[1]][0], clock=4)
    add("E:Extremum_Rechteck", "E Marken", "Marken auf Naht + Extremum der freien Energie", 2, lam=mpf(1) / 2 if abs(best_rect - 1) < mpf("0.006") else None)
    add("E:Extremum_frei", "E Marken", "Extremum der freien Energie ohne Spiegelung: Sechseck, lam = exp(i pi/3)", 1, lam=mpc(cos(pi / 3), sin(pi / 3)), c3=1 / (8 * pi * sqrt(3) / 2), clock=6)
    add("E:Tausch_Raum_Zeit", "E Marken", "kein ausgezeichnetes Seitenpaar (S tauscht Raum und Zeitkreis): fester Punkt a = b", 1, lam=mpf(1) / 2, clock=4)
    for n in (4, 5, 6):   # Klebegruppe von D_n: Z4 für ungerades n, Z2 x Z2 für gerades n
        cyc = n % 2 == 1
        add(f"E:Klebegruppe_D{n}", "E Marken", f"Klebegruppe von D{n} ({'Z4, zyklisch' if cyc else 'Z2 x Z2'}) wirkt transitiv auf die vier Marken", 1,
            lam=(mpf(1) / 2 if cyc else None), clock=(4 if cyc else 2))
    add("E:Seitenverh_2", "E Marken", "Rechteck 2:1 (zwei Quadrate)", 1, lam=findroot(lambda l: ellipk(1 - l) / ellipk(l) - 2, mpf("0.03")))

# ---------- Familie F: Zahlen phi_tree und delta_top aus Strukturzahlen (Suchraumgröße) ----------
def fam_numbers():
    S = [1, 2, 3, 4, 5, 6, 8]
    ratios = sorted({Fr(a, b) for a in S for b in S if Fr(1, 2) <= Fr(a, b) <= 4})
    # Fenster für phi bei z = +-3 (phi geht nur über log(1/phi) ein)
    def phi_needed(z):
        a = 1 / (CODATA + z * SIG); return exp(-(a ** 3 - 2 * C3 ** 3 * a ** 2) / (mpf(164) / 5 * C3 ** 6))
    lo, hi = sorted((phi_needed(3), phi_needed(-3))); a0 = 1 / CODATA; hits = []; tested = 0
    for r in ratios:
        for m in (2, 3, 4, 5, 6):
            for N in range(1, 1001):
                tested += 1; q = N * C3 ** m * exp(-2 * a0)
                if q >= mpf('0.5'): continue
                phi = mpf(r.numerator) / r.denominator * C3 + q * (1 - q) ** (-mpf(5) / 4)
                if lo <= phi <= hi: hits.append((str(r), N, m))
    for r, N, m in hits:
        add(f"F:phi={r}c3,delta={N}c3^{m}", "F Zahlen", f"phi_tree = {r} c3, delta_top = {N} c3^{m}", 3, ratio=mpf(Fr(r).numerator) / Fr(r).denominator, n48=mpf(N), m=m)
    return {"getestet": tested, "treffer_3sigma": hits, "fenster_rel": float((hi - lo) / lo)}

# ---------- Familie G: Koeffizienten der Kubik ----------
def fam_cubic():
    for A in (1, 2, 3, 4): add(f"G:A={A}", "G Kubik", f"quadratischer Koeffizient {A} c3^3", 1, A=mpf(A))
    for nm, b1 in (("SM 41/10", Fr(41, 10)), ("ohne Higgs 4", Fr(4)), ("zwei Higgs 21/5", Fr(21, 5)), ("MSSM 33/5", Fr(33, 5))):
        for fac in (4, 8, 16):
            add(f"G:b1={nm},x{fac}", "G Kubik", f"Logkoeffizient {fac} mal b1 ({nm})", 2, B=mpf(fac * b1.numerator) / b1.denominator)

# ---------- eigene Hypothesen ohne Code: hypotheses_extra.json = [{"id","fam","text","wahl","pred":{...}}] ----------
def fam_extra():
    import os
    if os.path.exists("hypotheses_extra.json"):
        for h in json.load(open("hypotheses_extra.json")):
            add(h["id"], h.get("fam", "X Eigene"), h.get("text", ""), h.get("wahl", 1), **{k: (mpf(str(v)) if isinstance(v, (int, float, str)) and k not in ("clock", "five") else v) for k, v in h["pred"].items()})

# ---------- Bewertung ----------
def evaluate(h):
    p = h["pred"]; res = {}; alpha_args = {}
    for key, val in p.items():
        if key in ("shift", "m"): continue
        if val is None: res[key] = "offen"; continue
        if key == "k" and val == "shift": continue
        t = TARGET[key]
        ok = (val == t) if isinstance(t, int) else abs(mpc(val) - t) < mpf(10) ** -12
        res[key] = "trifft" if ok else "verfehlt"
        if key in ("c3", "k", "expo", "ratio", "n48", "A", "B"): alpha_args[key] = val
    z = None
    if "shift" in p: z = 1.897 + p["shift"]; res["k"] = "trifft" if abs(z) < 3 else "verfehlt"
    elif alpha_args:
        if "m" in p: alpha_args["m"] = p["m"]
        z = alpha_z(**alpha_args)
        if z is None or abs(z) >= 3:
            for kk in alpha_args: res[kk] = "verfehlt" if res.get(kk) != "trifft" else res[kk]
    vals = list(res.values())
    verdict = "fällt" if "verfehlt" in vals else "trägt" if "trifft" in vals else "wählt nichts"
    return {**{k: v for k, v in h.items() if k != "pred"}, "ergebnis": res, "z_alpha": z, "urteil": verdict}

if __name__ == "__main__":
    fam_cut(); fam_state(); fam_expo(); fam_ade(); fam_marks(); fam_cubic(); fam_extra(); space = fam_numbers()
    out = [evaluate(h) for h in H]
    by = {}
    for o in out: by.setdefault(o["fam"], {"trägt": 0, "fällt": 0, "wählt nichts": 0})[o["urteil"]] += 1
    json.dump({"anzahl": len(out), "familien": by, "suchraum_F": space, "hypothesen": out}, open("results_sweep.json", "w"), indent=1, ensure_ascii=False, default=str)
    L = ["# Rangliste des Reihenprüfers", "", f"{len(out)} Hypothesen in {len(by)} Familien, automatisch bewertet.", "",
         "| Familie | trägt | fällt | wählt nichts |", "|---|---|---|---|"]
    L += [f"| {f} | {v['trägt']} | {v['fällt']} | {v['wählt nichts']} |" for f, v in by.items()]
    L += ["", "## Was trägt", "", "| Hypothese | Aussage | freie Wahlen | Treffer | z(α) |", "|---|---|---|---|---|"]
    for o in sorted([o for o in out if o["urteil"] == "trägt"], key=lambda o: (o["fam"], o["wahl"])):
        zt = "" if o["z_alpha"] is None else f"{o['z_alpha']:+.2f}"
        L.append(f"| {o['id']} | {o['text']} | {o['wahl']} | {', '.join(k for k, v in o['ergebnis'].items() if v == 'trifft')} | {zt} |")
    L += ["", "## Was nichts auswählt", ""] + [f"- {o['id']}: {o['text']}" for o in out if o["urteil"] == "wählt nichts"]
    L += ["", "## Was fällt (Auszug, nach Familie)", ""]
    for f in by:
        fl = [o for o in out if o["fam"] == f and o["urteil"] == "fällt"]
        if fl: L.append(f"- **{f}** ({len(fl)}): " + "; ".join(o["id"].split(":", 1)[1] + ("" if o["z_alpha"] is None else f" (z={o['z_alpha']:+.0f})") for o in fl[:14]) + (" …" if len(fl) > 14 else ""))
    L += ["", "## Suchraum der Zahlen (Familie F)", "", f"Getestet: {space['getestet']} Kombinationen phi_tree = (a/b)·c3, delta_top = N·c3^m. Treffer innerhalb 3 sigma: {len(space['treffer_3sigma'])}: {space['treffer_3sigma']}.",
          f"Relatives Fenster für phi: {space['fenster_rel']:.2e}."]
    L += ["", "## Hinweise zur Lesart", "",
          "- c3 allein unterscheidet nicht zwischen (Quadrat, chirale Nullmode) und (Rechteck 2:1, volle Nullmode): beide geben 1/(8 pi). Erst die Uhr (Ordnung 4) oder die Markenlage trennt sie.",
          "- thermisch omega = pi liegt numerisch bei k = 2 (z = +2), ist aber kein exakter Treffer.",
          "- 'trägt' heißt nur: kein Widerspruch zu den festen Zielgrößen. Es ist kein Beweis."]
    open("RANKING.md", "w").write("\n".join(L) + "\n")
    print(f"{len(out)} Hypothesen:", by); print("Suchraum F:", space["getestet"], "getestet,", space["treffer_3sigma"])
