"""w4: Formsuche. Der Rollenvertrag wird um die zwei Rollen erweitert, die aus den eingefrorenen Vorhersagen folgen
(Chiralitaet ohne Spiegel, paritaetsungerade Antwort des Weltzustands), und die Twistorbruecke des Ledgers steht
als Lesart in der Suche.  Aufruf: ../../../.venv/bin/python w4_formsuche.py
"""
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "research" / "originlab"))
import whole_solution_search as wss                  # noqa: E402
import ledger_map                                    # noqa: E402

LINES, CHECKS = [], []


def say(text=""):
    print(text)
    LINES.append(text)


def check(name, ok, detail=""):
    CHECKS.append(bool(ok))
    say(f"  [{'ok' if ok else 'FEHLER'}] {name}" + (f": {detail}" if detail else ""))


res = wss.run_shape_search()
asm, lev, g = res["assemblies"], res["leverage"], res["gates"]
say("w4: Formsuche mit Vorhersage-Rollen")
say()
say("Ledger-Karte")
counts = ledger_map.counts()
say("  " + ", ".join(f"{v} {k}" for k, v in counts.items()))
check("Twistorbruecke und Doppelbrechung stehen im Ledger",
      len(ledger_map.find("twistor")) >= 20 and len(ledger_map.find("birefringence")) >= 5,
      f"{len(ledger_map.find('twistor'))} und {len(ledger_map.find('birefringence'))} Behauptungen")
say()
say("Gatter (exakt, auf den 240 Wurzeln)")
check("Zentrum von Spin(10) = Zentrum von SU(2,2); halbe Periode = theta mal Uhr; volle Periode = (-1)^F",
      g["centres_identified"] and g["half_period_is_theta_times_clock"] and g["full_period_is_fermion_parity"])
check("Z4-Regel: eine Helizitaet je 16, Spiegel ausgeschlossen, Orientierung vertauscht",
      g["one_helicity_per_16"] and g["mirror_excluded"] and g["orientation_swaps"], json.dumps(g["rule"][1]))
check("keine Helizitaets-U(1) in E8, Z4 trivial auf Photonen", not g["helicity_u1_in_e8"] and g["z4_trivial_on_photons"])
say()
say("Verbuende")
check("26 Rollen, 11 Lesarten, keine vollstaendige Loesung", (res["roles"], asm["readings"], asm["complete"]) == (26, 11, []))
best = asm["best"][0]
check("bester Verbund Abstand 9", best["distance"] == 9, f"{best['members']} offen {best['open']}")
check("keine Lesart liefert Chiralitaet oder die paritaetsungerade Antwort aus ihrem Grundobjekt",
      asm["roles_no_reading_delivers"][-2:] == ["chirality", "parity_response"], str(asm["roles_no_reading_delivers"]))
say()
say("Hebel")
for name, row in sorted(lev["single"].items(), key=lambda kv: -kv[1]["gain"]):
    say(f"    {name:18s} Abstand {row['distance']}  Gewinn {row['gain']}")
top = lev["smallest_sets"][0]
check("fuenf Aussagen senken den Abstand von 9 auf 1; offen bleibt flavour_rest",
      (lev["distance_now"], lev["lowest_distance"], lev["statements_needed"]) == (9, 1, 5) and top["open"] == ["flavour_rest"],
      str(top["statements"]))
say(f"    Spannung: {res['tensions'][0]['statements']}: {res['tensions'][0]['text']}")
check("Kontrolle: eingepflanzte Loesung vollstaendig, ohne paritaetsungerade Antwort nicht",
      res["planted_complete"] and res["planted_without_parity_response"] == ["parity_response"])
say()
say(f"Pruefungen: {sum(CHECKS)}/{len(CHECKS)}")
(HERE / "ausgabe_w4.txt").write_text("\n".join(LINES) + "\n")
(HERE / "w4_ergebnis.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
raise SystemExit(0 if all(CHECKS) else 1)
