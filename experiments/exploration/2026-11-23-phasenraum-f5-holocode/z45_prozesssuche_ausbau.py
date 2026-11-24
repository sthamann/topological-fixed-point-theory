"""Z45: Ausbau der Prozesssuche um die gemeinsame Übergangsfamilie (OriginLab-Lane `process_rules`, zweiter Lauf).

Liest `research/originlab/process_rules.py` mit den neuen Universen (z44) und rechnet die Spin-2-Folgerung exakt nach:
Teil 1  Übergangsuniversum a(z) = z^-d (1 - q z^k): Vierteldrehung und Kern in einem Grad legen d = 1, k in {4, 8} fest, q nie.
Teil 2  Ausleseuniversum: Kein Rang-eins-Auslesen ist zugleich masselos und verbunden, keines universell gekoppelt.
Teil 3  Spin-2-Gedächtnis der Übergangsfamilie exakt: H/x = 1 - a x + b x^2 - ..., a = q^2/(1-q^2)^2, b = a (1+q^2)/(1-q^2)^2;
        a = 3/10 legt q^2 = (8 - sqrt 55)/3 fest; der nächste Koeffizient trennt Transfergedächtnis, Schurkern und Spektralwirkung.
Teil 4  Ganzer Prozess: neun lasttragende Prinzipien wählen Kern nu = 1, Transfergedächtnis, reelle 16er-Naht, Z3 x Z4;
        das Spiegelblatt wird überflüssig; der Lauf ohne Transfer (z43) bleibt exakt reproduzierbar.
"""
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp

HERE = Path(__file__).parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "research" / "originlab"))
from process_rules import Q_DIRAC, readouts, run, spin2_kernels, transfer_kernel, transfer_processes  # noqa: E402

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


res = run(state=False)          # Lane-Stand H6h24; die Zustandsregel für q steht in z46
U = res["universes"]

# ------------------------------------------------------------------ Teil 1
say("Teil 1: Übergangsuniversum a(z) = z^-d (1 - q z^k)")
tp = transfer_processes()
say(f"   {U['transfer_processes']['objects']} Kandidaten, {U['transfer_processes']['rules']} Regeln, "
    f"Urteile {U['transfer_processes']['verdicts']}")


def winners(rule):
    hits = [o for o in tp if set(rule) <= o["clauses"]]
    least = min((o["d"], o["k"]) for o in hits)
    return least, sorted(o["q"] for o in hits if (o["d"], o["k"]) == least)


qc = {(o["d"], o["k"]) for o in tp if "quarter_covariant" in o["clauses"]}
one = {(o["d"], o["k"]) for o in tp if "kernel_one_grade" in o["clauses"]}
say(f"   Vierteldrehung a(iz) ~ a(z): (d, k) in {sorted(qc)}; Kern in einem Grad: {sorted(one)}")
check("Vierteldrehungskovarianz gilt genau für k in {4, 8}; Kern in einem Grad genau für d = 1 und k in {4, 8}",
      qc == {(1, 4), (1, 8), (2, 4), (2, 8)} and one == {(1, 4), (1, 8)})
least, qs = winners(["quarter_covariant", "kernel_one_grade", "index_one"])
say(f"   Regel Vierteldrehung + Index eins + ein Grad: kleinster Prozess d, k = {least}, offen bleibt q in {qs}")
check("die vollständige Strukturregel wählt d = 1, k = 4 und lässt alle fünf q offen; keine Regel des Pools ist eindeutig",
      least == (1, 4) and len(qs) == 5 and U["transfer_processes"]["verdicts"]["UNIQUE"] == 0)
_, qs_t = winners(["quarter_covariant", "tau_i_nome"])
check("auch die tau = i-Nomen e^-pi, e^-2pi bleiben zu zweit: der Pool kennt keine Zahl, die q festlegt", len(qs_t) == 2)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Auslesen der primitiven Viertelübertragung")
ro = readouts()
for o in ro[:5]:
    say(f"   q = {o['q']:6s} {o['kind']:34s} {sorted(o['clauses'])}")
empty = U["readouts"]["minimal_empty_rules"]
say(f"   minimale leere Regeln: {empty}")
check("minimale leere Regeln sind genau {masseloser Pol, verbunden} und {universelles Residuum}",
      empty == [["connected", "massless_pole"], ["universal_residue"]])
qv = 0.5
Z = 1 - qv**2
resid = np.array([[Z*qv**(i + j) for j in range(4)] for i in range(4)])
check("Residuenmatrix der Quellen |4i> hat Rang 1 mit Diagonale (1-q^2) q^(2i): kohärent, aber nicht universell",
      np.linalg.matrix_rank(resid, tol=1e-12) == 1 and not np.allclose(np.diag(resid), resid[0, 0]))
befund("Ein masseloses, universell gekoppeltes Spin-2-Feld folgt nicht aus einem Rang-eins-Auslesen der Übertragung; "
       "es braucht eine Ward-/Eichstruktur (im Altbestand: GRAV.LINEAR.SOURCE.WARD.01 nur für vorgeschriebene Quellen)", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Spin-2-Gedächtnis der Übergangsfamilie exakt")
x, q = sp.symbols("x q", positive=True)
K = 1 - q**2
# H = K/C mit x C^2 + (x + q^2 - 1) C - 1 = 0  =>  K^2 + (x - K) K h - x h^2 = 0 für H = x h
h = sp.symbols("h0:5")
hs = sum(h[i] * x**i for i in range(5))
eq = sp.expand(K**2 + (x - K) * K * hs - x * hs**2)
sol = {}
for i in range(4):
    ci = eq.coeff(x, i).subs(sol)
    sol[h[i]] = sp.factor(sp.solve(ci, h[i])[0])
a_q, b_q = sp.factor(-sol[h[1]]), sp.factor(sol[h[2]])
say(f"   H/x = 1 - a x + b x^2 - c x^3:  a = {a_q},  b = {b_q},  c = {sp.factor(-sol[h[3]])}")
check("a(q) = q^2/(1-q^2)^2 und b(q) = q^2 (1+q^2)/(1-q^2)^4, also b/a^2 = (1+q^2)/q^2 > 1 (Källén-Lehmann b >= a^2)",
      sp.simplify(a_q - q**2 / (1 - q**2)**2) == 0 and sp.simplify(b_q - q**2 * (1 + q**2) / (1 - q**2)**4) == 0
      and sp.simplify(b_q / a_q**2 - (1 + q**2) / q**2) == 0)
qd2 = (8 - sp.sqrt(55)) / 3
a_d = sp.nsimplify(sp.simplify(a_q.subs(q, sp.sqrt(qd2))))
b_d = sp.radsimp(sp.simplify(b_q.subs(q, sp.sqrt(qd2))))
say(f"   a = 3/10 bei q^2 = (8 - sqrt 55)/3 = {float(qd2):.6f} (q = {float(sp.sqrt(qd2)):.6f}); dort b = {b_d} = {float(b_d):.6f}")
check("q^2 = (8 - sqrt 55)/3 ist die einzige Lösung in (0, 1) von q^2/(1-q^2)^2 = 3/10, und Q_DIRAC der Lane ist genau dieses q",
      a_d == sp.Rational(3, 10) and abs(float(sp.sqrt(qd2)) - Q_DIRAC) < 1e-15
      and [s for s in sp.solve(sp.Eq(q**2 / (1 - q**2)**2, sp.Rational(3, 10)), q) if s.is_positive and s < 1] == [sp.sqrt(qd2)])
C = transfer_kernel(Q_DIRAC)
xx = np.array([1e-3, 2e-3, 4e-3])
hn = (1 - Q_DIRAC**2) / C(xx) / xx
fit = np.polyfit(xx, hn, 2)
say(f"   numerisch aus der geschlossenen Randantwort: a = {-fit[1]:.6f}, b = {fit[0]:.4f}")
check("die geschlossene Randantwort der Lane gibt dieselben a und b", abs(-fit[1] - 0.3) < 1e-4 and abs(fit[0] - float(b_d)) < 2e-2)
b_dirac, b_schur = sp.Rational(1, 28), sp.Rational(1, 2) / sp.Rational(5, 3)**2
say(f"   nächster Koeffizient b: Spektralwirkung 1/28 = {float(b_dirac):.4f} (b < a^2 = 0.09, nicht positiv), "
    f"Schurkern 9/50 = {float(b_schur):.4f}, Transfergedächtnis {float(b_d):.4f}")
check("bei gleichem a = 3/10 trennen sich die drei Kerne im nächsten Koeffizienten: 1/28 < a^2 < 9/50 < b_transfer",
      b_dirac < sp.Rational(9, 100) < b_schur < b_d and b_schur == sp.Rational(9, 50))
a_pi = float(a_q.subs(q, sp.exp(-sp.pi)))
check("das tau = i-Nomen e^-pi gibt a = 0.0019, nicht 3/10: Vierteldrehungsnomen und Dirac-Korrektur passen in einer "
      "Übertragung nicht zusammen", abs(a_pi - 0.00187) < 1e-5)
sk = {o["name"]: o for o in spin2_kernels(state=False)}
check("alle fünf Transfergedächtnisse sind Källén-Lehmann-positiv; nur q_dirac trägt dirac_low_energy",
      all("kallen_lehmann" in o["clauses"] for n, o in sk.items() if n.startswith("transfer"))
      and [n for n, o in sk.items() if "dirac_low_energy" in o["clauses"] and n.startswith("transfer")] == ["transfer memory q=q_dirac"])
ZAHLEN["q_dirac_quadrat"] = str(qd2)
ZAHLEN["b_transfer"] = str(b_d)
ZAHLEN["b_werte"] = {"spektralwirkung": "1/28", "schur": "9/50", "transfer": float(b_d)}

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Ganzer Prozess mit der gemeinsamen Übertragung")
cp = res["composed_process"]
say(f"   {cp['objects']} Prozesse, {cp['rules']} Regeln, Urteile {cp['verdicts']}")
best = cp["selected_by_the_longest_rules"]
rule = best[0]["minimal_rules"][0]
say(f"   längste minimale Regel ({cp['longest_minimal_rule']} Prinzipien): {rule}")
say(f"   gewählt: {best[0]['selected']}")
head = "+".join(sorted(rule)) + " without "
for k, v in sorted(best[0]["ablation"].items()):
    if k.startswith(head):
        say(f"      ohne {k[len(head):]:20s}: {v[0]:9s} {str(v[1][:1])[:140]}")
check("neun Prinzipien wählen genau einen ganzen Prozess: Kern nu = 1, Transfergedächtnis q_dirac, 16 reelle Majorana mit "
      "umgekehrtem Blatt, Flavour Z3 x Z4", cp["longest_minimal_rule"] == 9 and len({r["selected"] for r in best}) == 1
      and all(p in best[0]["selected"] for p in ("nu=1 eps=0", "transfer memory q=q_dirac", "n=16 real reversed", "Z3 x Z4")))
check("jedes der neun Prinzipien trägt Last; das Spiegelblatt gehört nicht dazu (Orientierungsumkehr plus eine Übertragung "
      "ersetzen es)", all(v[1] != [best[0]["selected"]] or v[0] != "UNIQUE" for k, v in best[0]["ablation"].items()
                          if k.startswith(head)) and "mirror_sheet" not in rule and {"one_transfer", "orientation_reversal"} <= set(rule))
old = run(transfer=False)["composed_process"]
check("ohne Transferuniversen reproduziert die Lane den Stand H6h23 (acht Prinzipien, Schurkern)",
      old["longest_minimal_rule"] == 8 and "Schur w=1/2 m=5/3" in old["selected_by_the_longest_rules"][0]["selected"])
met = {m["universe"]: m["selected"] for m in res["selections_that_meet_the_record"]}
say(f"   Auswahlen, die einen Rekordbaustein treffen: {met}")
check("im Spin-2-Universum wählt Dirac-Korrektur plus Transfergedächtnis eindeutig q_dirac (Rekordtreffer)",
      any(m["universe"] == "spin2_kernels" and m["selected"] == "transfer memory q=q_dirac"
          for m in res["selections_that_meet_the_record"]))
befund("q wird nur durch die eine Zahl a = 3/10 festgelegt; keine Strukturklausel der Übertragung wählt q", True)
befund("Transfergedächtnis und Schurkern sind bei a = 3/10 im nächsten Koeffizienten unterscheidbar (b = 0.553 gegen 0.18); "
       "beide sind positiv, die Spektralwirkung (1/28) nicht", True)
ZAHLEN["prozess"] = best[0]["selected"]
ZAHLEN["prinzipien"] = rule
ZAHLEN["objekte"] = res["objects_enumerated"]

urteil = ("Die ausgebaute Prozesssuche bestätigt die elfte Analyse und schärft sie. Im Übergangsuniversum legen "
          "Vierteldrehung, Index eins und ein Kerngrad d = 1, k = 4 fest, q aber nie. Im Ausleseuniversum ist kein "
          "Rang-eins-Auslesen zugleich masselos und verbunden und keines universell gekoppelt: Das Graviton braucht eine "
          "Ward-/Eichstruktur, keine Nullmode der Übertragung. Das Spin-2-Gedächtnis der Übertragung ist exakt "
          "H/x = 1 - a x + b x^2 mit a = q^2/(1-q^2)^2 und b/a^2 = (1+q^2)/q^2; a = 3/10 legt q^2 = (8 - sqrt 55)/3 fest, "
          "dann ist b = 0.553, gegen 9/50 beim Schurkern und 1/28 bei der (nicht positiven) Spektralwirkung. Das "
          "tau = i-Nomen e^-pi gibt a = 0.0019; Nomen und Dirac-Korrektur passen in einer Übertragung nicht zusammen. Im "
          "zusammengesetzten Universum (67 200 Prozesse) wählen neun lasttragende Prinzipien genau einen Prozess: Kern "
          "nu = 1, Transfergedächtnis bei q_dirac, 16 reelle Majorana mit umgekehrtem Blatt, Flavour Z3 x Z4; das "
          "Spiegelblatt wird durch Orientierungsumkehr und eine gemeinsame Übertragung ersetzt. q bleibt eine Zahl aus "
          "der Dirac-Korrektur, keine Herleitung.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z45_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
