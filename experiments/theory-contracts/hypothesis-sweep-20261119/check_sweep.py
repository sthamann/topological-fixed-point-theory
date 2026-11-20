"""PASS/FAIL Prüfung des Reihenprüfers: läuft sweep.py und prüft die tragenden Aussagen."""
import json, subprocess, sys
from common import Checks
subprocess.run([sys.executable, "-B", "sweep.py"], check=True, capture_output=True)
R = json.load(open("results_sweep.json")); H = {h["id"]: h for h in R["hypothesen"]}
c = Checks("sweep")
c.ok("mindestens 140 Hypothesen in sieben Familien bewertet", R["anzahl"] >= 140 and len(R["familien"]) == 7, R["anzahl"])
c.ok("Kontrolle: die festgelegte Gleichung selbst trägt (z = +1,90) in den Familien B, C, F, G",
     all(abs(H[i]["z_alpha"] - 1.897) < 0.01 for i in ("B:Dirac+", "C:D5_Spinor", "F:phi=4/3c3,delta=48c3^4", "G:A=2", "G:b1=SM 41/10,x8")))
c.ok("Exponent: von 25 Gitterklassen trägt nur die Spinorklasse von D5", R["familien"]["C Exponent"] == {"trägt": 1, "fällt": 24, "wählt nichts": 0})
c.ok("Singularität: von 46 ADE Schnitten gibt nur Z5 = A4 die Fünferuhr", R["familien"]["D Singularität"]["trägt"] == 1 and H["D:A4(Z5)"]["urteil"] == "trägt")
c.ok("Zahlen: unter 95000 Kombinationen (a/b, N, m) liegt nur (4/3, 48, 4) innerhalb 3 sigma",
     R["suchraum_F"]["getestet"] == 95000 and R["suchraum_F"]["treffer_3sigma"] == [["4/3", 48, 4]])
c.ok("Kubik: nur Koeffizient 2 und 8 mal b1 = 41/10 tragen; ohne Higgs, zwei Higgs, MSSM fallen",
     R["familien"]["G Kubik"]["trägt"] == 2 and all(H[f"G:b1={n},x8"]["urteil"] == "fällt" for n in ("ohne Higgs 4", "zwei Higgs 21/5", "MSSM 33/5")))
c.ok("Marken: Spiegelung allein und Marken auf der Naht wählen nichts aus",
     H["E:Spiegelung_allein"]["urteil"] == "wählt nichts" and H["E:Marken_auf_Naht"]["urteil"] == "wählt nichts")
c.ok("Marken: Kleinsche Vierergruppe wählt nichts (fällt an der Uhr), eine Wirkung mit Ordnung 4 erzwingt lam = 1/2",
     H["E:Klein_transitiv"]["urteil"] == "fällt" and H["E:zyklisch_Z4"]["urteil"] == "trägt")
c.ok("Marken: Klebegruppe von D5 ist zyklisch und trägt; D4 und D6 (Z2 x Z2) fallen",
     H["E:Klebegruppe_D5"]["urteil"] == "trägt" and H["E:Klebegruppe_D4"]["urteil"] == "fällt" and H["E:Klebegruppe_D6"]["urteil"] == "fällt")
c.ok("Entartung erkannt: c3 allein trennt Quadrat mit chiraler Nullmode nicht von Rechteck 2:1 mit voller Nullmode",
     H["A:Quadrat (lam=1/2):I=4"]["urteil"] == "trägt" and H["A:Rechteck 2:1:I=2"]["urteil"] == "trägt" and H["E:Seitenverh_2"]["urteil"] == "fällt")
c.ok("freies Extremum (Sechseck) fällt an c3 und an der Uhr", H["E:Extremum_frei"]["urteil"] == "fällt")
c.done({"familien": R["familien"]})
