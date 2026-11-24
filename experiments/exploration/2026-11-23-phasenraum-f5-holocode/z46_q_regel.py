"""Z46: Eine Regel für q -- der Zustand der Viertelübertragung, und die Gegenprüfung der externen OriginLab-Befunde.

Teil 0  Gegenprüfung (externer Lauf vom 24.11.): Größengleichstand der 16er- und 32er-Naht, exakter Flavourwert
        (2 - sqrt 3)/6 statt 2/45, Z3 x Z4 = TM2 bei 15 Grad; TM1 bei 15 Grad passt zu den Daten, Herkunft offen.
Teil 1  Korrektur zu z45: a = 3/10 ist eine Skalenkonvention. Allgemein a = (3/10) f(0)/Int f, b = -f'(0)/(28 Int f);
        skalenfrei ist nur b/a^2, und das hängt vom Abschneider ab. Die Dirac-Zahl legt q nicht fest.
Teil 2  Zustandsregel: Der Nullvektor eta_q der Übertragung ist die Purifikation des untersten Strommodus der originalen
        Quarter-Gaußquelle (Calderón/Szegő, b = 4 pi). Ihre Norm omega(J_-1 J_1) = 1/(e^{4 pi} - 1) legt q = e^{-2 pi} fest;
        der Hardy-Zustand (J_1 Omega = 0) gibt q = 0.
Teil 3  Folgen: Einsteinpol mit Gewicht 1 - e^{-4 pi}, versteckte schmale Spin-2-Linie mit Gewicht e^{-4 pi} = 3.5e-6.
Teil 4  Lane: Zustandsuniversum q_rules und die skalenfreie Prozesssuche (ohne die Dirac-Zahl).
"""
import json
import sys
from itertools import permutations
from pathlib import Path

import numpy as np
import sympy as sp

HERE = Path(__file__).parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "research" / "originlab"))
from process_rules import (Q_DIRAC, Q_ORIGINAL, composed_process, flavour_pairs, h_exp, h_laplace, h_step,  # noqa: E402
                           ladder_occupation, original_current_norm, run, transfer_kernel)

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
U = res["universes"]
cp = res["composed_process"]

# ------------------------------------------------------------------ Teil 0
say("Teil 0: Gegenprüfung der externen OriginLab-Befunde")
best = cp["selected_by_the_longest_rules"][0]
ties = best["tie_set_without_size"]
say(f"   längste Regel ({cp['longest_minimal_rule']} Prinzipien) vor der Größenordnung: {len(ties)} Prozesse")
for t in ties:
    say(f"      {t[:118]}")
check("die längste Regel lässt die reelle 16er- und 32er-Naht zu; 16 gewinnt erst durch die Größenordnung",
      len(ties) == 2 and "seam n=16 real reversed" in ties[0] and "seam n=32 real reversed" in ties[1])
procs = composed_process()
rule_e8 = set(best["minimal_rules"][0]) - {"positivity"} | {"e8_seam"}
e8 = [o for o in procs if rule_e8 <= o["clauses"]]
check("mit dem Prinzip E8-Naht (c = 8) statt Positivität verschwindet der Gleichstand: genau ein Prozess", len(e8) == 1)
clock = [o for o in flavour_pairs() if o["pair"] == "Z3 x Z4"][0]
vals = sorted({float(sp.sympify(v)) for row in clock["pattern"] for v in row})
say(f"   Z3 x Z4 exakt: {clock['pattern'][0]}  (kleinster Eintrag {vals[0]:.10f}; 2/45 = {2/45:.10f})")
check("der kleinste Eintrag von Z3 x Z4 ist (2 - sqrt 3)/6 = 0.0446582, nicht 2/45 (bisher gerundet gedruckt)",
      abs(vals[0] - (2 - np.sqrt(3))/6) < 1e-12 and abs(vals[0] - 2/45) > 1e-4)

U0 = np.array([[np.sqrt(2/3), 1/np.sqrt(3), 0], [-1/np.sqrt(6), 1/np.sqrt(3), -1/np.sqrt(2)],
               [-1/np.sqrt(6), 1/np.sqrt(3), 1/np.sqrt(2)]])


def tm(kind, th, ph):
    c, s = np.cos(th), np.sin(th)*np.exp(-1j*ph)
    if kind == 1:                                  # erste TBM-Spalte bleibt: Drehung in der (2,3)-Ebene
        R = np.array([[1, 0, 0], [0, c, s], [0, -np.conj(s), c]])
    else:                                          # zweite Spalte bleibt: Drehung in der (1,3)-Ebene
        R = np.array([[c, 0, s], [0, 1, 0], [-np.conj(s), 0, c]])
    return U0 @ R


def canon(P):
    return min(tuple(tuple(round(P[r][c], 9) for c in cp_) for r in rp) for rp in permutations(range(3))
               for cp_ in permutations(range(3)))


def angles(Um):
    P = np.abs(Um)**2
    s13 = P[0, 2]
    return P[0, 1]/(1 - s13), s13, P[1, 2]/(1 - s13)


P_clock = [[float(sp.sympify(v)) for v in row] for row in clock["pattern"]]
th = np.pi/12
check("Z3 x Z4 ist genau TM2 bei 15 Grad (mittlere Spalte trimaximal, sin^2 theta13 = (2/3) sin^2 15 Grad)",
      canon(np.abs(tm(2, th, 0))**2) == canon(P_clock))
R12, R13, R23 = (0.27, 0.35), (0.020, 0.024), (0.43, 0.60)     # ungefähre 3-sigma-Bereiche wie z43 (NuFIT 6.0, NO)
fit = []
for ph in np.linspace(0, 2*np.pi, 721):
    s12, s13, s23 = angles(tm(1, th, ph))
    if R12[0] <= s12 <= R12[1] and R13[0] <= s13 <= R13[1] and R23[0] <= s23 <= R23[1]:
        fit.append((ph, s12, s13, s23))
s12, s13, s23 = angles(tm(1, th, np.pi/2))
say(f"   TM1 bei 15 Grad: sin^2 theta13 = (2 - sqrt 3)/12 = {s13:.5f}, sin^2 theta12 = {s12:.4f}, bei phi = pi/2 "
    f"sin^2 theta23 = {s23:.4f}; passende Phasen: {len(fit)} von 721")
check("TM1 bei 15 Grad (eine erhaltene Z2) liegt mit sin^2 theta13 = (2 - sqrt 3)/12 = 0.02233 und sin^2 theta12 = 0.318 in "
      "allen drei Bereichen, für ein Phasenintervall um phi = pi/2", len(fit) > 0 and abs(s13 - (2 - np.sqrt(3))/12) < 1e-12
      and R23[0] <= s23 <= R23[1])
s12b, s13b, _ = angles(tm(2, th, 0))
check("TM2 bei 15 Grad (= Z3 x Z4) verfehlt theta13 um den Faktor zwei", not (R13[0] <= s13b <= R13[1])
      and abs(s13b/s13 - 2) < 1e-12)
befund("Der Winkel 15 Grad steht schon in der Familienuhr (Z3 x Z4 = TM2 bei 15 Grad); die Daten verlangen ihn in TM1. "
       "Warum die erste statt der zweiten TBM-Spalte erhalten bleibt, ist nicht hergeleitet", True)
ZAHLEN["tm1_15"] = {"s13": s13, "s12": s12, "s23_bei_pi_halbe": s23}

# ------------------------------------------------------------------ Teil 1
say("\nTeil 1: Korrektur zu z45 -- a = 3/10 ist eine Skalenkonvention")
u = sp.symbols("u", nonnegative=True)
cut = {"e^-u": (sp.exp(-u), h_exp), "(1+u)^-3": ((1 + u)**-3, h_laplace(3)), "Stufe": (None, h_step)}
xs = np.array([2e-3, 4e-3, 6e-3, 8e-3])
rows = {}
for name, (f, H) in cut.items():
    if f is None:
        f0, f1, I = 1, 0, 1
    else:
        f0, f1, I = f.subs(u, 0), sp.diff(f, u).subs(u, 0), sp.integrate(f, (u, 0, sp.oo))
    a_th, b_th = sp.Rational(3, 10)*f0/I, -sp.Rational(1, 28)*f1/I
    hx = np.asarray(H(xs), float)/xs
    co = np.polyfit(xs, hx, 2)
    rows[name] = (a_th, b_th, -co[1], co[0])
    say(f"   {name:9s} a = {str(a_th):5s} (num {-co[1]:.5f}), b = {str(b_th):5s} (num {co[0]:.4f}), "
        f"b/a^2 = {str(sp.nsimplify(b_th/a_th**2)):6s}")
check("a = (3/10) f(0)/Int f und b = -f'(0)/(28 Int f) stimmen für e^-u, (1+u)^-3 und die Stufe mit dem Lane-Kern",
      all(abs(float(r[0]) - r[2]) < 2e-3 and abs(float(r[1]) - r[3]) < 3e-2 for r in rows.values()))
check("a hängt vom Abschneider ab (3/10, 3/5, 3/10) und skaliert mit Lambda^2: eine Konvention, keine Zahl der Theorie",
      rows["(1+u)^-3"][0] == sp.Rational(3, 5) and rows["e^-u"][0] == sp.Rational(3, 10))
check("das skalenfreie b/a^2 ist abschneiderabhängig (25/63, 25/42, 0) und für alle drei < 1: kein positiver Kern trifft es",
      [sp.nsimplify(r[1]/r[0]**2) for r in rows.values()] == [sp.Rational(25, 63), sp.Rational(25, 42), 0])
q, s_ = sp.symbols("q s", positive=True)
a_q = q**2/(1 - q**2)**2
check("jedes Transfergedächtnis trifft a = 3/10 mit passender Skala mu^2 = a(q) Lambda^2/(3/10); q_dirac verlangt "
      "zusätzlich mu = Lambda", all(float(a_q.subs(q, qq)) > 0 for qq in (0.1, 0.5, 0.9)))
befund("Korrektur zu z45/H6h24: 'q durch a = 3/10 festgelegt' gilt nur bei gleicher Skala mu = Lambda und e^-u-Normierung; "
       "die Dirac-Zahl ist keine Regel für q", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Zustandsregel aus der originalen Quarter-Gaußquelle")
om = {al: original_current_norm(alpha=al) for al in (0.25, 0.75)}
exact = 1/(np.exp(4*np.pi) - 1)
say(f"   omega(J_-1 J_1) aus der Fermisumme: alpha = 1/4 {om[0.25]:.12e}, alpha = 3/4 {om[0.75]:.12e}; "
    f"1/(e^(4 pi) - 1) = {exact:.12e}")
check("die Fermisumme Sum_r c_(r+1)(1 - c_r) mit c_r = (1 + e^(4 pi (r + alpha)))^-1 gibt für beide Quarter "
      "1/(e^(4 pi) - 1) (Quelle: ELLIPTIC-CURRENT-SOURCE Gl. 6, 8)", all(abs(v/exact - 1) < 1e-10 for v in om.values()))
occ = q**2/(1 - q**2)
check("<n> = Sum n (1-q^2) q^(2n) = q^2/(1-q^2) (numerisch für q = 1/2)",
      abs(sum(n*(1 - 0.25)*0.25**n for n in range(200)) - 0.25/0.75) < 1e-14)
w = sp.symbols("w", positive=True)
sol = sp.solve(sp.Eq(occ.subs(q, sp.sqrt(w)), 1/(sp.exp(4*sp.pi) - 1)), w)
say(f"   Leiterbesetzung von eta_q: <n> = q^2/(1 - q^2); gleichgesetzt mit omega: q^2 = {sol}")
check("<n> = q^2/(1-q^2) ist streng monoton: genau eine Lösung q^2 = e^(-4 pi), also q = e^(-2 pi)",
      len(sol) == 1 and sp.simplify(sol[0] - sp.exp(-4*sp.pi)) == 0)
check("die Lane misst dieselbe Besetzung am Nullvektor der Kompression (q = e^-2pi)",
      abs(ladder_occupation(Q_ORIGINAL)/exact - 1) < 1e-8 and abs(Q_ORIGINAL - np.exp(-2*np.pi)) < 1e-18)
say("   (dieselbe Dichtematrix mit Deckgenerator 2H und Normierung 2 pi: das Verhältnis hängt nicht an der Normierung)")
check("der Hardy-Zustand (J_1 Omega = 0) gibt <n> = 0, also q = 0 und kein Gedächtnis", ladder_occupation(0.0) == 0.0)
alts = {m: np.exp(-2*np.pi*m) for m in (1, 2, 3)}
check("ordnet man den Leiterschritt dem Strommodus m zu, folgt q = e^(-2 pi m); das größte q ist e^(-2 pi)",
      all(abs(ladder_occupation(v)*(np.exp(4*np.pi*m) - 1) - 1) < 1e-6 for m, v in alts.items()))
befund("q = e^-2pi folgt unter der Zuordnung 'Übertragungsschritt S^4 = ein Quant von J_-1' und der expliziten Prämisse der "
       "freien Gaußquantisierung; die Quelle selbst nennt b = 4 pi keine physikalische Temperatur", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Folgen für Spin 2")
qo = Q_ORIGINAL
Z = 1 - qo**2
band = ((1 - qo)**2, (1 + qo)**2)
C = transfer_kernel(qo)
alpha_y = 4/3*qo**2/Z
say(f"   Einsteinpol-Gewicht 1 - e^(-4 pi) = {Z:.10f}; verstecktes Band [{band[0]:.5f}, {band[1]:.5f}] mu^2 mit Gewicht "
    f"{qo**2:.4e}; b/a^2 = 1 + e^(4 pi) = {1 + np.exp(4*np.pi):.4e}")
check("das Band ist eine schmale Linie bei mu^2 (Breite 4q = 0.0075) mit Gewicht e^(-4 pi) = 3.487e-6",
      abs(band[1] - band[0] - 4*qo) < 1e-15 and abs(qo**2 - 3.487e-6) < 1e-9)
check("x C(x) -> 1 - e^(-4 pi) bei x -> 0 (Residuum des masselosen Pols)", abs(1e-8*C(1e-8) - Z) < 1e-7)
say(f"   statischer Yukawa-Zusatz eines massiven Spin-2 (Faktor 4/3): alpha = {alpha_y:.3e} bei Reichweite 1/mu")
befund("Im originalen Zustand ist die Gravitation Einstein bis auf eine massive Spin-2-Linie mit Kopplung 4.6e-6; die Skala "
       "mu legt die Zustandsregel nicht fest", True)
ZAHLEN["q"] = "e^-2pi"
ZAHLEN["versteckt"] = qo**2
ZAHLEN["yukawa_alpha"] = alpha_y

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Zustandsuniversum und skalenfreie Prozesssuche")
qr = U["q_rules"]
say(f"   q_rules: {qr['objects']} Zustände, Urteile {qr['verdicts']}")
for name, d in qr["derived"].items():
    say(f"      {name:16s} <n> = {d['occupation']:.3e}  Polgewicht = {d['massless_residue']:.7f}  {d['clauses']}")
sel = {r["selected"]: r["minimal_rules"] for r in qr["unique_selections"]}
say(f"   eindeutige Auswahlen: {sel}")
check("je eine Zustandsregel wählt je ein q: originaler Zustand -> e^-2pi, Hardy -> 0, Gleichgewicht -> 1/sqrt 2, "
      "Dirac bei gleicher Skala -> q_dirac", sel.get("e^-2pi") == [["original_state"]] and sel.get("0 (Hardy state)") ==
      [["hardy_state"]] and sel.get("1/sqrt2") == [["equal_weight"]] and sel.get("q_dirac") == [["dirac_same_scale"]])
emp = qr["minimal_empty_rules"]
check("originaler Zustand verträgt sich weder mit gleichem Blattgewicht noch mit der Dirac-Zahl bei gleicher Skala",
      ["equal_weight", "original_state"] in emp and ["dirac_same_scale", "original_state"] in emp)
free = cp["without_scale_convention"]
fb = free["selected_by_the_longest_rules"][0]
say(f"   ohne dirac_spin2: längste Regel {free['longest_minimal_rule']} Prinzipien {fb['minimal_rules']}")
say(f"   gewählt: {fb['selected'][:150]}")
for k, v in sorted(fb["ablation"].items()):
    h = "+".join(sorted(fb["minimal_rules"][0])) + " without "
    if k.startswith(h):
        say(f"      ohne {k[len(h):]:20s}: {v[0]:9s} {str(v[1][:1])[:110]}")
check("ohne die Skalenkonvention wählen acht Prinzipien genau einen Prozess: Kern nu = 1, Transfergedächtnis im "
      "originalen Zustand, reelle 16er-Naht mit umgekehrtem Blatt, Z3 x Z4",
      free["longest_minimal_rule"] == 8 and all(p in fb["selected"] for p in
                                                ("nu=1 eps=0", "transfer memory q=e^-2pi", "n=16 real reversed", "Z3 x Z4")))
tsz = sorted(len(v) for v in fb["tie_sets_without_size"].values())
check("eine der beiden minimalen Regeln (mit E8-Naht) braucht keine Größenordnung mehr", tsz == [1, 2])
h0 = "+".join(sorted(fb["minimal_rules"][0])) + " without original_state"
check("ohne die Zustandsregel fällt die Auswahl auf reines Einstein (H = x, Hardy-Zustand) zurück",
      "Einstein H = x" in fb["ablation"][h0][1][0])
check("mit der Dirac-Zahl bleibt der Lauf aus z45 (neun Prinzipien, q_dirac) erhalten",
      cp["longest_minimal_rule"] == 9 and "q=q_dirac" in best["selected"])
befund("Die Wahl zwischen q_dirac und e^-2pi ist eine Wahl zwischen einer Skalenkonvention und einem originalen Zustand; "
       "der Zustand ist die tragfähigere Regel, bleibt aber eine Prämisse (nicht aus P1 hergeleitet)", True)
ZAHLEN["prozess_skalenfrei"] = fb["selected"]
ZAHLEN["prinzipien_skalenfrei"] = fb["minimal_rules"]

urteil = ("Die q-Regel ist eine Zustandsregel. Die Dirac-Zahl a = 3/10 legt q nicht fest: a = (3/10) f(0)/Int f hängt vom "
          "Abschneider ab und skaliert mit Lambda^2; skalenfrei ist nur b/a^2 (25/63, 25/42, 0), immer < 1. Die originale "
          "Quarter-Gaußquelle (Calderón/Szegő, b = 4 pi) hat die Stromnorm omega(J_-1 J_1) = 1/(e^{4 pi} - 1); setzt man die "
          "Leiterbesetzung q^2/(1-q^2) des Nullvektors gleich, folgt eindeutig q = e^{-2 pi}, der Hardy-Zustand gibt q = 0. "
          "Damit ist Spin 2 ein Einsteinpol mit Gewicht 1 - e^{-4 pi} und eine schmale massive Linie mit Gewicht 3.5e-6 "
          "(Yukawa-Stärke 4.6e-6, Reichweite offen). Ohne die Skalenkonvention wählen acht Prinzipien genau diesen Prozess "
          "(Kern nu = 1, Transfergedächtnis im originalen Zustand, reelle 16er-Naht, Z3 x Z4); ohne Zustandsregel fällt die "
          "Wahl auf reines Einstein. Gegenprüfung: Die 32er-Naht ist nur durch die Größenordnung ausgeschlossen (mit E8-Naht "
          "nicht mehr), Z3 x Z4 gibt exakt (2 - sqrt 3)/6 und ist TM2 bei 15 Grad; die Daten passen zu TM1 bei 15 Grad, "
          "dessen Herkunft offen ist. Die Zuordnung Leiterschritt = J_-1 und die Gaußquantisierung bleiben Prämissen.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z46_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
