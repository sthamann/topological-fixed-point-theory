"""Z173: Abstand zur Komplettlösung als Kennzahl der Programmsuche, die offene Liste und die Kernstrategie.

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke bleiben unberührt, nichts wird befördert.
Statuskürzel: E exakt, C bedingt, L Lesart, N numerisch, O offen.

Korrektur: z171 und z172 haben am Alpha-Block gearbeitet. Alpha ist ein gelöster Block und keine Lücke. Ziel ist die
Komplettlösung: alles aus der einen Annahme A*, ohne weitere Annahmen und ohne freie Zahl.

Die Kennzahl steht in research/originlab/reconstruction_search.py (complete_solution_distance).

Teil A  Abstand der besten Szenarien zur Komplettlösung
Teil B  Die offene Liste, nach fünf Fragen geordnet
Teil C  Eine Quelle für drei Fragen: die Standardform des Nahtzustands
Teil D  Bilanz und Urteil.

Hinweis (z174): Die Prüfungen dieses Skripts beschreiben den Stand z173 (108 Szenarien, Abstand 2/1/7, elf Punkte). Seit z174 gibt das
Suchprogramm 162 Szenarien, Abstand 1/1/7 und zehn Punkte aus; die Ergebnisdatei z173_ergebnis.json bleibt als Stand vor z174 stehen.
"""
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parents[2] / "research" / "originlab"))
import clock_sector_search as cs                                           # noqa: E402
import reconstruction_search as rs                                         # noqa: E402
import whole_solution_search as ws                                         # noqa: E402

OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val, text):
    BEFUNDE[name] = val
    say(f"   -> {name}: {text}")


d = rs.complete_solution_distance()
ZAHLEN["abstand"] = d["distance"]
ZAHLEN["nach_fragen"] = d["by_question"]

say("Teil A  Abstand der besten Szenarien zur Komplettlösung")
check("A1 [E] Ziel: alle zwölf Klauseln, nur die Annahme A*, keine freie Zahl, K19 und K20 geschlossen. Bewertet werden "
      "die 108 zulässigen und in sich stimmigen Szenarien", d["scenarios"] == 108 and not d["complete"])
check("A2 [E] Abstand der besten Szenarien: 2 fehlende Klauseln, 1 freie Zahl, 7 Annahmen außer A*",
      d["distance"] == {"missing_clauses": 2, "free_numbers": 1, "premises_beyond_axiom": 7})
check("A3 [E] P_IF zählt nicht mehr als Annahme: Die zwei Verklebungen D5 + A3 -> E8 gehen durch Ladungskonjugation "
      "einer Seite ineinander über und geben beide E8 (z144, z168). Das Bit sagt nur, welche von zwei "
      "spiegelgleichen Welten links heißt", d["conventions"]["P_IF"]["holds"])
befund("Abstand", "2 + 1 + 7, dazu K19 und K20",
       "Elf offene Punkte trennen den Bestand von der Komplettlösung. Jede Runde wird ab jetzt an dieser Zahl gemessen.")

say("\nTeil B  Die offene Liste, nach fünf Fragen geordnet")
check("B1 [E] fünf Fragen: Ursprung (K19, K20), Zusammensetzung (P_ord, P_uhr, P_univ), Zuordnung (I_step, P_lift), "
      "Gravitation (P_G, P_curv, Klausel Graviton), Statistik (Klausel modulare Statistik)",
      d["by_question"] == {"origin": ["K19", "K20"], "composition": ["P_ord", "P_uhr", "P_univ"],
                           "gravity": ["P_G", "P_curv", "graviton"], "identification": ["I_step", "P_lift"],
                           "statistics": ["modular_statistics"]} and d["open_total"] == 11)
for q, items in d["by_question"].items():
    say(f"   {q}: " + ", ".join(items))
check("B2 [E] kein offener Punkt liegt in einem gelösten Block (Teilchen, Familien, Massen, Mischung, alpha, E8)",
      not any(x in json.dumps(d["by_question"]) for x in ("alpha", "flavour", "mass", "E8")))
befund("Liste", "elf Punkte in fünf Fragen",
       "Ursprung zwei, Zusammensetzung drei, Zuordnung zwei, Gravitation drei, Statistik einer.")

say("\nTeil C  Eine Quelle für drei Fragen: die Standardform des Nahtzustands")
comp = cs.completion_candidates()
diag = [o for o in comp if "diagonal" in o["clauses"]]
check("C1 [E] die Standardform (Summe H_a (x) H_a-quer über die 16 Uhrsektoren, Partner durch die Modulkonjugation) ist "
      "genau eine der 27 Vervollständigungen: die diagonale. Sie ist bosonisch, also modular invariant, mit positiven "
      "Vielfachheiten und allen 16 Sektoren", len(comp) == 27 and len(diag) == 1
      and {"bosonic", "clock_sectors_present", "graph"} <= diag[0]["clauses"])
check("C2 [E] Preis: In ihr sind die 16 und die 4 Skalare (Dimension 5/4 und 3/4), keine holomorphen Fermionen; der "
      "Spin kommt dann nicht aus der Statistik der Naht, sondern aus der Helizität auf der Geraden (z168)",
      {"scalar_16", "scalar_4"} <= diag[0]["clauses"] and "holomorphic_fermions" not in diag[0]["clauses"])
net = ws.run_composition_search()
check("C3 [E/N] derselbe Partner gibt das einzige reelle lorentzsche Netz des Katalogs (z172) und macht die "
      "Paritätserhaltung nativ (z130)", net["lorentzian"] == ["1 partner | hermitian (2, 2)"])
befund("Quelle", "Kandidat, nicht eingetragen",
       "Die Standardform des Nahtzustands (GNS und Tomita-Takesaki) liefert Partner, Zeit, positive Vielfachheiten und "
       "eine modular invariante Vervollständigung in einem. Sie ist der Kandidat für Statistik und Zusammensetzung "
       "(vier der elf Punkte). Eingetragen wird er erst, wenn Chiralität, Uhrphase und Ordnung dort gerechnet sind.")

say("\nTeil D  Bilanz und Urteil")
urteil = ("Korrektur der Richtung: alpha ist gelöst und kein Ziel. Die Programmsuche führt jetzt den Abstand zur "
          "Komplettlösung als Kennzahl: 2 fehlende Klauseln, 1 freie Zahl, 7 Annahmen außer A*, dazu K19 und K20, "
          "zusammen elf Punkte in fünf Fragen. P_IF ist als Konvention erkannt und zählt nicht mehr. Kernstrategie: "
          "alles aus der Standardform des Nahtzustands rekonstruieren; sie ist der gemeinsame Kandidat für Statistik "
          "und Zusammensetzung. Keine Komplettlösung, Prüfbank unverändert, nichts befördert.")
say("   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z173_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z173.txt").write_text("\n".join(OUT) + "\n")
