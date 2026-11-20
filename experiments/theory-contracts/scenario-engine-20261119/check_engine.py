"""Prüft den Objektvorrat an bekannten Tatsachen und den Lauf auf Wiederholbarkeit. Aufruf: python3 -B check_engine.py"""
import json, subprocess, sys
from fractions import Fraction as Fr
import generators as g, gen_mehr as m
ok = []
def chk(label, cond): ok.append((label, bool(cond))); print("ok  " if cond else "FAIL", label)
T = {o["id"]: o for o in g.gen_teilgitter()}; has = lambda o, a, v: any(f[0] == a and f[1] == v for f in o["feats"])
full = [o for o in T.values() if o["id"].endswith("<E8") and any(f == ("zahl", 8, "Rang") for f in o["feats"])]
chk("E8 hat 15 Teilsysteme vom vollen Rang (mit E8 selbst)", len(full) == 15)
chk("D5 + A3 in E8: Index 4, Klebung Z4, Normen 5/4 und 3/4", all(has(T["gitter:D5 + A3<E8"], *x) for x in (("zahl", 4), ("typ", "Z4"), ("bruch", "5/4"), ("bruch", "3/4"))))
chk("A4 + A4 in E8: Index 5, Klebung Z5", has(T["gitter:A4 + A4<E8"], "typ", "Z5") and any(f == ("zahl", 5, "Index im Elterngitter") for f in T["gitter:A4 + A4<E8"]["feats"]))
S = {o["id"]: o for o in g.gen_spiegelungsgruppen()}
chk("G31: Ordnung 46080, 60 Spiegelungen, Zentrum 4", all(f in S["spiegel:G31"]["feats"] for f in (("zahl", 46080, "Ordnung"), ("zahl", 60, "Zahl der Spiegelungen"), ("zahl", 4, "Ordnung des Zentrums"))))
chk("W(D5) hat Ordnung 1920, W(E8) 696729600 mit 120 Spiegelungen", ("zahl", 1920, "Ordnung") in S["spiegel:W(D5)"]["feats"] and ("zahl", 696729600, "Ordnung") in S["spiegel:G37 (E8)"]["feats"] and ("zahl", 120, "Zahl der Spiegelungen") in S["spiegel:G37 (E8)"]["feats"])
M = {o["id"]: o for o in g.gen_modulkurven()}
chk("X0(11) hat Geschlecht 1, X(7) Geschlecht 3 und 168 Symmetrien, X(5) hat 12 Spitzen", has(M["modul:X0(11)"], "zahl", 1) and ("zahl", 3, "Geschlecht") in M["modul:X(7)"]["feats"] and has(M["modul:X(7)"], "zahl", 168) and ("zahl", 12, "Spitzen") in M["modul:X(5)"]["feats"])
F = {o["id"]: o for o in g.gen_flaechen()}
chk("jede elliptische Fläche der Liste hat Eulerzahl 12", all(("zahl", 12, "Eulerzahl") in o["feats"] for o in F.values()))
C = {o["id"]: o for o in m.gen_modulformen()}
chk("Parafermionen Z5 haben zentrale Ladung 8/7, SU(2) Stufe 5 hat 15/7", has(C["cft:PF5"], "bruch", "8/7") and any(has(o, "bruch", "15/7") for o in g.gen_anyonen() if o["id"] == "su2:5"))
chk("Linsenraum Z5 mit Windung 0: eta 4/5", any(has(o, "bruch", "4/5") for o in g.gen_linsen() if o["id"] == "linse:S3/Z5,k=0"))
r1 = subprocess.run([sys.executable, "-B", "suche.py", "--fakes", "150", "--seed", "7"], capture_output=True, text=True).stdout
r2 = subprocess.run([sys.executable, "-B", "suche.py", "--fakes", "150", "--seed", "7"], capture_output=True, text=True).stdout
chk("gleicher Startwert gibt denselben Lauf", r1 == r2 and len(r1) > 500)
subprocess.run([sys.executable, "-B", "suche.py", "--fakes", "1000", "--seed", "1"], capture_output=True)
R = json.load(open("results.json")); A, B = R["laeufe"]
chk("Anker steht im Lauf über alles vorn und schlägt den Zufall", A["einzeln"][0]["id"] == "gitter:D5 + A3<E8" and A["einzeln"][0]["p_gesamt"] <= 0.002)
chk("Rest ohne Anker: kein Einzelobjekt, Paar oder Tripel unter p = 0,05", min(x["p_gesamt"] for k in ("einzeln", "paare", "tripel") for x in B[k]) > 0.05)
json.dump({"name": "engine", "passed": sum(c for _, c in ok), "total": len(ok), "checks": [{"label": l, "ok": c} for l, c in ok]}, open("results_check.json", "w"), ensure_ascii=False, indent=1)
print(f"{sum(c for _, c in ok)}/{len(ok)}"); sys.exit(0 if all(c for _, c in ok) else 1)
