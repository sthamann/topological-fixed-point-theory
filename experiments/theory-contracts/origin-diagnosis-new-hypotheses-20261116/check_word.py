#!/usr/bin/env python3
"""Geordnetes E8 Wort auf dem Materieraum, 16. November 2026 (Fortsetzung ORIGINDIAG).

Vier Marken im E8 kovarianten Stromvertrag, je Mittelstrom gamma aus (10,6):
  Z0   = Tr P                                   (reflektierte Leerreferenz)
  ZU   = Tr P h h'                              (neutraler Cartanmarker h)
  ZY   = Tr P X_-gamma P' X_gamma               (geladener Mittelstrom, hin und zurück)
  ZUY  = Tr P h X_-gamma P' h' X_gamma          (Marker vor und nach dem Strom)
P = Projektor auf 128 Materiewurzeln (voll) oder auf drei Familien (96).
P' = 1 (bedingte Erwartung, Rückweg durch den ganzen Raum) oder P (Kompression).
Alle vier Größen sind quadratisch in den Strukturkonstanten, also unabhängig
von der Wurzelsektion. Exakte Bruchrechnung, nur Standardbibliothek.
"""
from fractions import Fraction as Fr
from itertools import combinations, product
import json
import sys

OUT, FAILS = {}, []


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name}" + (f" :: {detail}" if detail != "" else ""))
    OUT[name] = {"ok": bool(cond), "detail": str(detail)}
    if not cond:
        FAILS.append(name)


R = []
for i, j in combinations(range(8), 2):
    for si, sj in product((2, -2), repeat=2):
        v = [0] * 8
        v[i], v[j] = si, sj
        R.append(tuple(v))
R += [s for s in product((1, -1), repeat=8) if s.count(-1) % 2 == 0]
RS = set(R)
add = lambda a, b: tuple(x + y for x, y in zip(a, b))
dot = lambda h, r: sum(Fr(x) * y for x, y in zip(h, r)) / 2  # r verdoppelt
is_int = lambda r: all(x % 2 == 0 for x in r)
CAR, FAM = range(5), range(5, 8)
M = [r for r in R if not is_int(r)]
BIF = [r for r in R if is_int(r) and sum(1 for i in CAR if r[i]) == 1]
fam = lambda r: tuple(r[i] for i in FAM)
FAM4 = sorted({fam(r) for r in M if [r[i] for i in CAR].count(-1) % 2 == 0})
check("128 Materiewurzeln, 60 Mittelströme, 4 Familiengewichte", (len(M), len(BIF), len(FAM4)) == (128, 60, 4))

Yh = (Fr(-1, 3), Fr(-1, 3), Fr(-1, 3), Fr(1, 2), Fr(1, 2), 0, 0, 0)
T3 = (0, 0, 0, Fr(1, 2), Fr(-1, 2), 0, 0, 0)
F1 = (0, 0, 0, 0, 0, 1, -1, 0)   # Familiencartan
F2 = (0, 0, 0, 0, 0, 1, 1, -2)
MARK = {"Y": Yh, "T3": T3, "F1": F1, "F2": F2}


def marks(P, g, mode):
    """mode: 'full' (P = M), 'cond' (Rückweg frei), 'comp' (Zwischenzustand in P)."""
    Pset = set(P)
    ng = tuple(-x for x in g)
    z0 = len(P)
    steps = [a for a in P if add(a, g) in RS and (mode != "comp" or add(a, g) in Pset)]
    zy = Fr(len(steps), z0)
    zu = {k: sum(dot(h, a) ** 2 for a in P) / z0 for k, h in MARK.items()}
    zuy = {k: sum(dot(h, a) * dot(h, add(a, g)) for a in steps) / z0 for k, h in MARK.items()}
    return zy, zu, zuy


print("\nA. Voller Materieraum (128), alle 60 Mittelströme")
full = {g: marks(M, g, "full") for g in BIF}
zy_full = {v[0] for v in full.values()}
check("ZY/Z0 = 1/4 für jeden Mittelstrom auf dem vollen Raum", zy_full == {Fr(1, 4)}, zy_full)
zu_full = full[BIF[0]][1]
print("   ZU/Z0 voll:", {k: str(v) for k, v in zu_full.items()})
check("ZU/Z0 = |h|^2/4 auf dem vollen Raum (isotrop)",
      all(zu_full[k] == sum(Fr(x) ** 2 for x in MARK[k]) / 4 for k in MARK))
# geschlossene Formel: ZUY = Cov(h) * ZY - (h.gamma)^2/4 * ZY
ok_formula = True
for g in BIF:
    zy, zu, zuy = full[g]
    steps = [a for a in M if add(a, g) in RS]
    for k, h in MARK.items():
        hg = dot(h, g)
        cov = sum((dot(h, a) + hg / 2) ** 2 for a in steps) / len(steps)
        ok_formula &= zuy[k] == zy * (cov - hg ** 2 / 4)
check("ZUY = ZY * (Kovarianz quer zu gamma  -  (h.gamma)^2/4), exakt für alle 60 Ströme", ok_formula)
vals = sorted({(str(dot(Yh, g) ** 2), str(full[g][2]["Y"])) for g in BIF})
print("   Hyperladungsmarker: ((Y.gamma)^2, ZUY/Z0) =", vals)
vals_t = sorted({(str(dot(T3, g) ** 2), str(full[g][2]["T3"])) for g in BIF})
print("   T3 Marker:          ((T3.gamma)^2, ZUY/Z0) =", vals_t)
OUT["A_full_Y"] = vals
OUT["A_full_T3"] = vals_t
conn = sorted({str(full[g][2]["Y"] - full[g][1]["Y"] * full[g][0]) for g in BIF})
print("   verbunden ZUY - ZU*ZY (Y):", conn)
check("Verbundene gemischte Antwort des Hyperladungsmarkers ist für jeden Mittelstrom strikt negativ",
      all(full[g][2]["Y"] - full[g][1]["Y"] * full[g][0] < 0 for g in BIF), conn)

print("\nB. Drei Familien: bedingte Erwartung gegen Kompression")
tab = {}
agree_cond, agree_comp = True, True
for w4 in FAM4:
    w4n = tuple(-x for x in w4)
    P = [r for r in M if fam(r) not in (w4, w4n)]
    assert len(P) == 96
    for g in BIF:
        c = marks(P, g, "cond")
        k = marks(P, g, "comp")
        f = full[g]
        key = (str(c[0]), str(k[0]),
               str(c[1]["Y"]), str(c[2]["Y"]), str(k[2]["Y"]), str(f[2]["Y"]),
               str(c[1]["F1"]), str(f[1]["F1"]))
        tab[key] = tab.get(key, 0) + 1
        agree_cond &= c[0] == f[0] and all(c[1][m] == f[1][m] and c[2][m] == f[2][m] for m in ("Y", "T3"))
        agree_comp &= k[0] == f[0]
print("   Häufigkeit | ZY cond | ZY comp | ZU(Y) cond | ZUY(Y) cond | ZUY(Y) comp | ZUY(Y) voll | ZU(F1) cond | ZU(F1) voll")
for key, n in sorted(tab.items()):
    print("  ", n, "|", " | ".join(key))
OUT["B_table"] = {" | ".join(k): n for k, n in tab.items()}
check("Bedingte Erwartung auf drei Familien liefert Z0, ZU, ZY, ZUY des vollen Raums für alle Trägermarker (Y, T3), "
      "jede Wahl der vierten Familie und jeden Mittelstrom", agree_cond)
check("Kompression liefert einen anderen Wert für ZY (1/6 statt 1/4)", not agree_comp and
      {k[1] for k in tab} == {"1/6"}, {k[1] for k in tab})
# Familienmarker unterscheiden drei von vier Familien
fam_diff = {(k[6], k[7]) for k in tab}
check("Familienmarker sehen die ausgelassene Familie: ZU(F1) weicht auf drei Familien vom vollen Wert ab",
      any(a != b for a, b in fam_diff), fam_diff)

print("\nC. Geometrischer Kern an den Marken 1, i, -1 in der Winkelkarte")
z1, z2, z3 = 1, 1j, -1
two = z1 * z3 / (z1 - z3) ** 2
three = z1 * z2 * z3 / ((z1 - z2) * (z2 - z3) * (z1 - z3))
print("   Zweipunkt:", two, " Dreipunkt:", three)
check("Zweipunktkern -1/4, Dreipunktkern -i/4 (Betrag 1/4 wie im Original)", abs(two + 0.25) < 1e-15 and abs(three + 0.25j) < 1e-15)

print("\nErgebnis:", "alle Prüfungen bestanden" if not FAILS else f"FEHLER in {FAILS}")
json.dump(OUT, open(__file__.replace("check_word.py", "results_word.json"), "w"), indent=1, ensure_ascii=False)
sys.exit(1 if FAILS else 0)
