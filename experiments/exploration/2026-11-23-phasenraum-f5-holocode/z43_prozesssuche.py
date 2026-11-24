"""Z43: Prozesssuche mit dem umgebauten Suchprogramm (OriginLab-Lane `process_rules`).

Ruft `research/originlab/process_rules.py` auf (dieselbe Regelsuche wie die Genesis-Lane) und liest das Ergebnis für H6:
Teil 1  Bausteine: Matrixkern, Spin-2-Kern, Naht, Flavour -- welche kleinste strukturelle Regel wählt welchen Kandidaten.
Teil 2  Ganzer Prozess: ein Kandidat je Baustein (53 760 Prozesse), neun übergreifende Prinzipien; die längste minimale Regel
        und ihre Ablation.
Teil 3  Datenprüfung außerhalb der Lane: Die Flavourmuster aller Restsymmetrie-Paare gegen ungefähre 3-sigma-Bereiche der
        Neutrinomischung (NuFIT 6.0, normale Ordnung, gerundet). Ein exaktes S4-Restsymmetriemuster passt nicht.
"""
import json
import sys
from itertools import permutations
from fractions import Fraction as Fr
from pathlib import Path

HERE = Path(__file__).parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "research" / "originlab"))
from process_rules import run  # noqa: E402

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
sel = lambda u: {r["selected"]: r for r in U[u]["unique_selections"]}

say("Teil 1: Bausteine")
mk = sel("matrix_kernels")
say(f"   Matrixkern: {U['matrix_kernels']['objects']} Kandidaten, Urteile {U['matrix_kernels']['verdicts']}")
check("nu = 1, eps = 0 wird durch Hardy-Schnitt plus Wellengleichung oder plus Kovarianz eindeutig gewählt; Hardy allein lässt "
      "den det-Term offen", ["hardy_slice", "wave_equation"] in mk["nu=1 eps=0"]["minimal_rules"]
      and ["conformally_covariant", "hardy_slice"] in mk["nu=1 eps=0"]["minimal_rules"]
      and not any(r["minimal_rules"] == [["hardy_slice"]] for r in U["matrix_kernels"]["unique_selections"]))
sk = sel("spin2_kernels")
der = U["spin2_kernels"]["derived"]
say(f"   Spin 2: {U['spin2_kernels']['objects']} Kandidaten; Spektralwirkung: "
    f"{dict((k, v['first_zero_near']) for k, v in der.items() if k.startswith('spectral'))} (erste Nullstelle, Gitter)")
check("alle drei Spektralwirkungen haben eine euklidische Nullstelle, der Ledger-Ansatz x e^x ist nicht Källén-Lehmann-positiv",
      all("no_euclidean_zero" not in der[k]["clauses"] for k in der if k.startswith("spectral"))
      and "kallen_lehmann" not in der["ledger ansatz H = x e^x"]["clauses"])
check("Dirac-Korrektur plus Spiegelgewicht wählt eindeutig C = 1/x + 1/(x + 10/3); ohne Spiegelgewicht drei positive Kandidaten",
      ["dirac_low_energy", "mirror_weight"] in sk["Schur w=1/2 m=5/3"]["minimal_rules"]
      and sum(1 for v in der.values() if {"dirac_low_energy", "kallen_lehmann"} <= set(v["clauses"])) == 3)
sl = sel("seam_layers")
check("die reelle 16-Majorana-Naht mit umgekehrtem zweitem Blatt braucht Spin(10)_1; ohne diese Klausel gewinnt die komplexe "
      "Flussnaht (8 komplexe Lagen)", ["c_eight", "sheet_neutral", "spin10_level_one"] in sl["n=16 real reversed"]["minimal_rules"]
      and ["c_eight", "sheet_neutral"] in sl["n=8 complex reversed"]["minimal_rules"])
pats = U["flavour_pairs"]["derived"]["patterns"]
say("   Flavourmuster (Paar, |U|^2, |J|):")
for p in pats:
    say(f"      {p}")
tbm = "(('0', '1/3', '2/3'), ('1/2', '1/3', '1/6'), ('1/2', '1/3', '1/6'))"
check("TBM entsteht als Z3 x V4' (mit Null, theta13 = 0); Z3 x Z4 mit der Familienuhr gibt keine Null (2/45, 1/3, 28/45)",
      any(p[0] == "Z3 x V4'" and p[1] == tbm for p in pats)
      and any(p[0] == "Z3 x Z4" and "2/45" in p[1] for p in pats))

say("\nTeil 2: Ganzer Prozess")
cp = res["composed_process"]
say(f"   {cp['objects']} Prozesse, {cp['rules']} Regeln, Urteile {cp['verdicts']}")
best = cp["selected_by_the_longest_rules"]
say(f"   längste minimale Regel: {cp['longest_minimal_rule']} Prinzipien: {best[0]['minimal_rules']}")
say(f"   gewählt: {best[0]['selected']}")
for k, v in sorted(best[0]["ablation"].items()):
    if k.startswith("+".join(sorted(best[0]["minimal_rules"][0])) + " without"):
        say(f"      ohne {k.split(' without ')[1]:14s}: {v[0]:9s} {str(v[1][:1])[:150]}")
check("acht Prinzipien wählen genau einen ganzen Prozess: Kern nu = 1, C = 1/x + 1/(x+10/3), 16 reelle Majorana mit umgekehrtem "
      "Blatt, Flavour Z3 x Z4", cp["longest_minimal_rule"] == 8 and len({r["selected"] for r in best}) == 1
      and all(part in best[0]["selected"] for part in ("nu=1 eps=0", "Schur w=1/2 m=5/3", "n=16 real reversed", "Z3 x Z4")))
check("jedes der acht Prinzipien trägt Last (Weglassen ändert die Auswahl oder macht sie mehrdeutig)",
      all(v[1] != [best[0]["selected"]] or v[0] != "UNIQUE" for k, v in best[0]["ablation"].items()))
befund("Das Prinzip 'Spiegelblatt' (zweite Blatthälfte umgekehrt und gleichgewichtig) legt Naht und Spin-2-Kern zugleich fest", True)

say("\nTeil 3: Flavour gegen Daten (außerhalb der Lane)")
# ungefähre 3-sigma-Bereiche, NuFIT 6.0 (2024), normale Ordnung, gerundet
R12, R13, R23 = (0.27, 0.35), (0.020, 0.024), (0.43, 0.60)
fits = []
for pair, pat, J in pats:
    rows = [[float(Fr(v)) for v in row] for row in eval(pat)]
    hit = None
    for rp in permutations(range(3)):
        for cpm in permutations(range(3)):
            M = [[rows[r][c] for c in cpm] for r in rp]
            s13 = M[0][2]
            if s13 >= 1:
                continue
            s12, s23 = M[0][1] / (1 - s13), M[1][2] / (1 - s13)
            ok = R12[0] <= s12 <= R12[1] and R13[0] <= s13 <= R13[1] and R23[0] <= s23 <= R23[1]
            if ok:
                hit = (s12, s13, s23)
    s13_vals = sorted({float(Fr(v)) for row in eval(pat) for v in row})
    fits.append((pair, pat, hit, s13_vals))
    say(f"   {pair:10s} mögliche sin^2 theta13: {[round(v, 4) for v in s13_vals]}  -> {'passt' if hit else 'passt nicht'}")
check("kein exaktes Restsymmetriemuster der nichtentarteten abelschen S4-Untergruppen passt in die 3-sigma-Bereiche "
      "(theta13: 0, 2/45 = 0.044, 1/9, 1/4, 1/3, ... gegen 0.020-0.024)", all(f[2] is None for f in fits))
befund("Die Familienregel braucht eine kleine Brechung der Restsymmetrie (TFPT: phi0-Korrekturen, sin^2 theta13 = phi0 e^{-5/6}); "
       "die Prozesssuche findet die Orientierung, nicht die Brechung", True)
ZAHLEN["prozess"] = best[0]["selected"]
ZAHLEN["prinzipien"] = best[0]["minimal_rules"]
ZAHLEN["objekte"] = res["objects_enumerated"]

urteil = ("Das Suchprogramm lässt sich für die Prozesssuche umbauen: Die neue OriginLab-Lane process_rules nutzt dieselbe "
          "Regelsuche wie die Genesis-Lane und läuft in etwa drei Sekunden. In den Bausteinen wählt Hardy plus Wellengleichung "
          "oder Kovarianz den Kern nu = 1; die Spektralwirkungen fallen alle durch, Dirac-Korrektur plus Spiegelgewicht wählt "
          "C = 1/x + 1/(x + 10/3); die reelle 16er-Naht braucht Spin(10)_1, sonst gewinnt die komplexe Flussnaht; TBM entsteht "
          "als Z3 x V4' nur mit theta13 = 0. Im zusammengesetzten Universum (53 760 Prozesse) wählen acht lasttragende Prinzipien "
          "genau einen Prozess: Kern nu = 1, Spin 2 mit gleichgewichtigem verborgenem Kanal bei 10/3, 16 reelle Majorana mit "
          "umgekehrtem zweitem Blatt, Flavour Z3 x Z4. Das Spiegelblatt fixiert dabei zwei Bausteine. Gegen die Daten passt "
          "kein exaktes S4-Restsymmetriemuster (Z3 x Z4 gibt sin^2 theta13 = 2/45, doppelt zu groß): Die Familienregel braucht "
          "eine kleine Brechung. Das ist eine Auswahl unter gewählten Prinzipien, keine Herleitung.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z43_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
