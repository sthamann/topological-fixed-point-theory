"""w1: Weltgenerator und Gesamtlösungssuche (Exploration 2026-11-27).

Ruft die Lane `research/originlab/whole_solution_search.py` (Abschnitt D) auf: Keime (Blätter k, Marken n, Träger N,
Familien f, Dimension d) werden zu ganzen Welten ausgerechnet, alle Regeln aus höchstens fünf Prinzipien wählen ihre
kleinste Welt, Lesarten werden zu Verbünden zusammengesetzt und am Rollenvertrag einer vollständigen Lösung gemessen.
Ausgabe: ausgabe_w1.txt, w1_ergebnis.json.
"""
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "research" / "originlab"))
import whole_solution_search as wss                  # noqa: E402

LINES, CHECKS = [], []


def say(text=""):
    LINES.append(text)
    print(text)


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    say(f"[{'ok' if ok else 'FEHLER'}] {name}" + (f": {detail}" if detail else ""))


res = wss.run_world_search()
rules, coh, num, asm, lev, audit = (res[k] for k in ("rules", "coherent", "numbers", "assemblies", "leverage", "audit"))

say("A  Universum und Regelsuche")
check("16 464 Keime, 26 Prinzipien (19 der aufgezeichneten Welt, 7 Gegenprinzipien), 83 682 Regeln bis Länge 5",
      rules["worlds"] == 16464 and len(rules["clauses"]) == 26 and rules["rules"] == 83682, str(rules["verdicts"]))
top = rules["uniquely_selected_worlds"]
check("Die Welt k2 n4 N5 f3 d4 wird von 6692 Regeln eindeutig gewählt, die nächste (d = 2) von 3249",
      top[0] == ("k2 n4 N5 f3 d4", 6692) and top[1] == ("k2 n4 N5 f3 d2", 3249), str(top[:4]))
check("399 minimale Regeln: 60 mit drei, 287 mit vier, 52 mit fünf Prinzipien",
      rules["minimal_rules_by_length"] == {3: 60, 4: 287, 5: 52})
pure = [r["rule"] for r in rules["minimal_rules"] if len(r["rule"]) == 3 and not set(r["rule"]) & set(wss.WORLD_TARGETS)]
say(f"    Dreierregeln ohne Messwert ({len(pure)}):")
for r in pure:
    say("      " + " + ".join(r))

say("\nB  Welche Welten schließen")
check("Flach, Träger und Familien aus dem Markendivisor, Z4-Anomalie, Verklebung, eindeutige positive Teilung: "
      "nur der Keim (2, 4, 5, 3)", coh["closing_seeds_without_dimension"] == [[2, 4, 5, 3]])
check("Die Z3-Welt (3, 3, 4, 2) scheitert an Anomalie, Verklebung und Teilung",
      len(coh["flat_divisor_worlds_and_failures"]["k3 n3 N4 f2"]) == 4)
check("Der Schluss allein lässt die Dimension frei (2 bis 8); mit den Raumzeitprinzipien bleibt d = 4",
      coh["dimensions_allowed_by_closure_alone"] == list(range(2, 9))
      and coh["closing_with_spacetime_clauses"] == ["k2 n4 N5 f3 d4"])

say("\nC  Zahlen der Welten")
t = num["tfpt"]["numbers"]
say(f"    Keim (2, 5, 3): alpha^-1 = {t['alpha_inverse']:.10f}, sqrt(phi0(1-phi0)) = {t['cabibbo']:.5f}, "
    f"1/3 - phi0/2 = {t['solar']:.5f}, (1-2c3) phi0 = {t['baryon']:.5f}")
check("Von 336 Keimen (210 mit positivem Seed) trifft einer alle vier Zahlen; alpha trifft nur (2, 5, 3)",
      num["hits_histogram"] == {0: 276, 1: 25, 2: 30, 3: 4, 4: 1} and num["per_target"]["alpha_measured"] == ["k2 N5 f3"])
check("Cabibbo und Omega_b prüfen im Wesentlichen nur N = 5, der Sonnenwinkel nur f = 3; die Trennschärfe trägt alpha",
      len(num["per_target"]["cabibbo"]) == 35 and len(num["per_target"]["solar"]) == 30)

say("\nD  Lesarten und Verbünde")
check("10 Lesarten, 386 Verbünde, kein Widerspruch, keine vollständige Lösung",
      (asm["readings"], asm["assemblies"], asm["contradictions"], asm["complete"]) == (10, 386, 0, []))
for name, row in asm["single_readings"].items():
    say(f"    {name:18s} Keim erzwungen: {'ja' if row['seed_forced'] else 'nein':4s} kleinste Welt {row['seed_smallest']}"
        f"  offen {row['distance']}")
check("Allein erzwingt nur die Lesart sky_spin2 den ganzen Keim",
      [m for m, r in asm["single_readings"].items() if r["seed_forced"]] == ["sky_spin2"])
best = asm["best"][0]
check("Bester Verbund: anomaly_minimal + background_free + real_form, Abstand 7",
      best["members"] == ["anomaly_minimal", "background_free", "real_form"] and best["distance"] == 7, str(best["open"]))
check("Keine Lesart liefert aus ihrem Grundobjekt: rp, time, continuum, einstein, area, flavour_rest",
      asm["roles_no_reading_delivers"] == ["rp", "time", "continuum", "einstein", "area", "flavour_rest"])

say("\nE  Hebel")
for name, row in sorted(lev["single"].items(), key=lambda r: -r[1]["gain"]):
    say(f"    {name:16s} Abstand {row['distance']}  Gewinn {row['gain']}")
check("Drei Aussagen (ladder, equilibrium, irreducible) senken den Abstand von 7 auf 1; offen bleibt flavour_rest",
      (lev["distance_now"], lev["lowest_distance"], lev["statements_needed"]) == (7, 1, 3)
      and lev["smallest_sets"][0]["open"] == ["flavour_rest"], str(lev["smallest_sets"][0]))

say("\nF  Kann die Suche eine vollständige Lösung enthalten?")
check("Eine eingepflanzte vollständige Lesart wird als vollständig erkannt; ohne die Rolle Zeit, ohne die "
      "Markenklausel oder mit einem unverbundenen Partner nicht",
      audit["planted_complete"] and audit["planted_without_time"] == ["time"]
      and "marks4" in audit["planted_without_marks_clause"] and audit["planted_with_unbridged_partner"] == 1)
say(f"    Aus dem Keim gerechnet: {', '.join(audit['computed_from_the_seed'])}")
say(f"    Aus dem Bestand gelesen: {', '.join(audit['read_from_the_record'])}")
lad = res["ladders"]["so(4,2)"]
check("Leiter von so(4,2): Kern positiv, Zeit positiv und kontinuierlich, Weyl-Dimension des Raums 3",
      lad["kernel_positive"] and lad["time_continuous"] and lad["spacetime_dimension"] == 4, str(lad["bottom_and_gap"]))

passed = sum(ok for _, ok in CHECKS)
say(f"\n{passed}/{len(CHECKS)} Prüfungen bestanden")
res["checks"] = {"passed": passed, "total": len(CHECKS), "failed": [n for n, ok in CHECKS if not ok]}
(HERE / "ausgabe_w1.txt").write_text("\n".join(LINES) + "\n")
(HERE / "w1_ergebnis.json").write_text(json.dumps(res, indent=1, ensure_ascii=False, default=str) + "\n")
sys.exit(0 if passed == len(CHECKS) else 1)
