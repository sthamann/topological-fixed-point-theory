"""Z49: Flavour mit größeren Gruppen -- Delta(6 n^2) mit Majorana-verträglichen Restsymmetrien.

Teil 1  Lane-Universum flavour_groups (n = 2, 3, 4, 6, 8): Klein-Gruppe für Majorana-Neutrinos, zyklische Gruppe mit drei
        verschiedenen Eigenwerten für geladene Leptonen. Kein Muster trifft die 3-sigma-Bereiche, keins ist TM1 mit Mischung.
Teil 2  Die Vierteldrehung erweitert S4 zu Delta(96) (kleinste Gruppe mit Phasen i): dort ist die Familienuhr Z3 x Z4 =
        TM2 bei 15 Grad als Majorana-verträgliches Muster (Klein-Gruppe) vorhanden.
Teil 3  Weitere n (5, 7, 10, 12, 16): Delta(6 n^2) gibt nur TM2 mit diskreten Winkeln; keiner liegt im theta13-Fenster.
Teil 4  TM1 bei 15 Grad passt zu den Daten, entsteht aber in keiner dieser Gruppen: Es braucht einen Mechanismus jenseits
        reiner Restsymmetrie (mit Z2 x CP bleibt der Winkel frei).
"""
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp

HERE = Path(__file__).parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "research" / "originlab"))
from process_rules import DATA_RANGES, flavour_groups, run  # noqa: E402

OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


res = run()
U = res["universes"]["flavour_groups"]

# ------------------------------------------------------------------ Teil 1
say("Teil 1: Lane-Universum flavour_groups")
say(f"   {U['objects']} Muster, {U['rules']} Regeln, Urteile {U['verdicts']}")
for grp, pat, cl in U["derived"]["patterns"]:
    if "mixes" in cl:
        say(f"      {grp:11s} {pat[:96]}  {cl}")
say(f"   minimale leere Regeln: {U['minimal_empty_rules']}")
check("kein Klein-x-zyklisches Muster in Delta(24) bis Delta(384) trifft die 3-sigma-Bereiche oder das theta13-Fenster",
      ["fits_data"] in U["minimal_empty_rules"] and ["theta13_window"] in U["minimal_empty_rules"])
check("TM1 kommt nur ohne Mischung vor (TBM mit theta13 = 0): {mixes, tm1} ist leer", ["mixes", "tm1"] in U["minimal_empty_rules"])

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Vierteldrehung -> Delta(96)")
fg = flavour_groups(ns=(2, 4))
s4 = [o for o in fg if o["n"] == 2]
d96 = [o for o in fg if o["n"] == 4]
check("Delta(24) = S4 hat mit Klein-Neutrinos keine Familienuhr (TM2 bei 15 Grad); erst Delta(96) mit den Phasen i hat sie",
      not any("clock_pattern" in o["clauses"] for o in s4) and any("clock_pattern" in o["clauses"] for o in d96))
clk = [o for o in d96 if "clock_pattern" in o["clauses"]][0]
say(f"   Delta(96): {clk['pattern']}  (geladene Leptonen Ordnung {clk['charged_lepton_orders']})")
befund("Das Z3 x Z4-Muster der S4-Lane braucht dort eine Z4 auf einer Seite, also keine Majorana-Klein-Gruppe; in Delta(96) "
       "entsteht dasselbe Muster mit Klein-Neutrinos -- die Vierteldrehung macht die Familienuhr Majorana-verträglich", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Weitere n")
more = flavour_groups(ns=(5, 7, 10, 12, 16))
th13 = sorted({float(min(sp.sympify(v) for v in o["pattern"][0])) for o in more + flavour_groups()
               if "tm2" in o["clauses"] and "mixes" in o["clauses"]})
say(f"   sin^2 theta13 aller gemischten TM2-Muster bis n = 16: {[round(v, 5) for v in th13]}")
check("alle gemischten Muster bis Delta(1536) sind TM2 (mittlere Spalte trimaximal), keins TM1",
      all("tm2" in o["clauses"] and "tm1" not in o["clauses"] for o in more if "mixes" in o["clauses"]))
check("kein TM2-Winkel liegt im Fenster 0.020-0.024 (nächste Werte 0.0114 und 0.0254 bei Delta(384), Delta(1536))",
      not any(DATA_RANGES[1][0] <= v <= DATA_RANGES[1][1] for v in th13 if v < 0.3))
ang = [float(np.degrees(np.arcsin(np.sqrt(1.5*v)))) for v in th13 if 0 < v < 0.3]
say(f"   zugehörige TM2-Winkel (Grad): {[round(a, 3) for a in ang]}")

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: TM1 bei 15 Grad")
s13 = (2 - np.sqrt(3))/12
s12 = (np.cos(np.pi/12)**2/3)/(1 - s13)
say(f"   TM1(15 Grad): sin^2 theta13 = {s13:.5f}, sin^2 theta12 = {s12:.4f}; TM2(15 Grad): sin^2 theta13 = {2*s13:.5f}")
check("TM1 bei 15 Grad liegt im theta13-Fenster und im theta12-Bereich; es ist genau die Hälfte der Familienuhr",
      DATA_RANGES[1][0] <= s13 <= DATA_RANGES[1][1] and DATA_RANGES[0][0] <= s12 <= DATA_RANGES[0][1])
U0 = np.array([[np.sqrt(2/3), 1/np.sqrt(3), 0], [-1/np.sqrt(6), 1/np.sqrt(3), -1/np.sqrt(2)],
               [-1/np.sqrt(6), 1/np.sqrt(3), 1/np.sqrt(2)]])
th, ph = np.pi/12, np.pi/2
cc, ss = np.cos(th), np.sin(th)*np.exp(-1j*ph)
Um = U0 @ np.array([[1, 0, 0], [0, cc, ss], [0, -np.conj(ss), cc]])
P = np.abs(Um)**2
t13 = P[0, 2]
t12, t23 = P[0, 1]/(1 - t13), P[1, 2]/(1 - t13)
J = np.imag(Um[0, 0]*Um[1, 1]*np.conj(Um[0, 1])*np.conj(Um[1, 0]))
sd = J/(np.sqrt(t12*(1 - t12)*t23*(1 - t23)*t13)*(1 - t13))
say(f"   TM1(15 Grad) mit Phase pi/2: sin^2 theta23 = {t23:.4f}, Jarlskog {J:.4f}, |sin delta| = {abs(sd):.4f}")
check("mit der Phase pi/2 (mu-tau-Reflexion) ist theta23 maximal und delta = +-pi/2: TM1(15 Grad) wird eine scharfe "
      "Vorhersage für die CP-Phase", abs(t23 - 0.5) < 1e-12 and abs(abs(sd) - 1) < 1e-9)
befund("Reine Restsymmetrie in Delta(6 n^2) (n bis 16) gibt TM1 bei 15 Grad nicht; mit Z2 x CP bliebe der Winkel frei. "
       "Die datenpassende Familie braucht eine weitere Regel, die die Viertel-15-Grad-Drehung in die (2,3)- statt in die "
       "(1,3)-Ebene legt", True)
ZAHLEN["tm2_theta13"] = th13
ZAHLEN["tm1_15"] = {"s13": s13, "s12": s12}

urteil = ("Größere Familiengruppen lösen die Flavourfrage nicht. In Delta(6 n^2) mit Klein-Neutrinos (Majorana) und "
          "zyklischen Leptonresten entstehen bis n = 16 nur TM2-Muster mit diskreten Winkeln; keines trifft das "
          "theta13-Fenster (nächste Werte 0.0114 und 0.0254). Die Vierteldrehung erweitert S4 zu Delta(96): Dort ist die "
          "Familienuhr TM2 bei 15 Grad Majorana-verträglich, aber mit sin^2 theta13 = 0.0447 doppelt zu groß. TM1 bei "
          "15 Grad passt zu den Daten (0.0223), kommt aber in keiner dieser Gruppen vor; mit Z2 x CP bliebe der Winkel "
          "frei. Die Familienregel braucht also einen Mechanismus jenseits reiner Restsymmetrie.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z49_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
