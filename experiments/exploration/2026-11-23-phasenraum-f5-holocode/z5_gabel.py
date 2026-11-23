"""Z5: Die Gabel c = 8 (Energiebedingung im Mittel, Prämisse F) gegen c = 24 (als Operatorgleichung, F').

Drei Prüfungen.
 1. Verträgt sich die Operatorlesart mit dem Zustand von H6? H6 bezieht Zeit und Kragenbreite aus dem
    KMS-Zustand bei beta = L (K6). Ist dieser Zustand ein Eigenzustand mit Energie null?
 2. Wo sitzen die drei Familien? In einem E8 sind die drei geraden Spinstrukturen die drei Terme der
    GSO-Summe chi = (theta2^8 + theta3^8 + theta4^8)/(2 eta^8). Wie wirkt SL(2,Z) auf ihnen? Und TFPT
    bettet die Familien in ein E8 ein: 248 = 78 + 8 + 2 * 27 * 3 (E6 x SU(3)).
 3. Ist die 3D Wand aus F beobachtbar? Eine chirale Naht mit c hat das Anomaliepolynom (c/24) p1.
    Verschiebt man die Wand über ein 4D Gebiet, ändert sich die Phase um exp(2 pi i (c/24) int p1)
    = exp(2 pi i c sigma/8), denn int p1 = 3 sigma. Die Wand ist unsichtbar, wenn das für jede
    Signatur sigma eins ist.
"""
import json
from pathlib import Path

import numpy as np
from mpmath import mp, mpc, jtheta, exp, pi, qp

mp.dps = 25
OK = []


def check(name, cond):
    OK.append((name, bool(cond)))
    return "ja" if cond else "NEIN"


z4 = json.loads(Path(__file__).with_name("z4_ergebnis.json").read_text())

print("1. Operatorlesart gegen den Zustand von H6\n")
var = z4["var_L0_at_i"]
print(f"   KMS-Zustand bei beta = L: <L0> - c/24 = 0, aber Var(L0) = {var:.6f} > 0.")
print("   Er ist also kein Eigenzustand von H. Unter H|psi> = 0 hätte jeder physische Zustand die Energie null;")
print("   das Boltzmanngewicht e^(-beta H) wäre auf ihnen eins, und beta (also die Zeit, die Kragenbreite")
print("   L/2 und die Zeitneutralität aus Z3/Z4) wäre nicht mehr bestimmt.")
check("KMS-Zustand kein Eigenzustand", var > 0.1)
print("   Die Operatorlesart widerspricht damit K6 von H6. Innerhalb von H6 gilt die Lesart im Mittel.")

print("\n2. Familien in einem E8\n")


def eta(t):
    return exp(1j * pi * t / 12) * qp(exp(2j * pi * t))


def sector(a, t):
    return (jtheta(a, 0, exp(1j * pi * t)) / eta(t)) ** 8 / 2


t0 = mpc(0.17, 1.13)
names = {2: "theta2 (Spinstruktur AP)", 3: "theta3 (AA)", 4: "theta4 (PA)"}
S_map, T_map, T_phase = {}, {}, set()
for a in (2, 3, 4):
    sv = sector(a, -1 / t0)
    S_map[a] = next(b for b in (2, 3, 4) if abs(sv - sector(b, t0)) < 1e-15)
    tv = sector(a, t0 + 1)
    b = next(b for b in (2, 3, 4) if abs(abs(tv) - abs(sector(b, t0))) < 1e-15)
    T_map[a] = b
    T_phase.add(round(float(mp.arg(tv / sector(b, t0)) / (2 * pi)), 9))
print(f"   S: {S_map} (AA fest, AP und PA vertauscht, ohne Phase)")
print(f"   T: {T_map}, gemeinsame Phase / 2 pi = {T_phase}")
check("S auf Sektoren", S_map == {2: 4, 3: 3, 4: 2})
check("T auf Sektoren", T_map == {2: 2, 3: 4, 4: 3} and T_phase == {-0.333333333})
print("   Bis auf die gemeinsame Rahmungsphase wirkt SL(2,Z) auf den drei Sektoren als S3 = SL(2, Z/2).")
print("   Ein E8 trägt also genau drei Spinstruktursektoren; das passt zu C' (drei Familien). In E8 hoch 3")
print("   gäbe es drei Blöcke mit je drei Sektoren, neun statt drei.")
dims = {"E6": 78, "SU3": 8, "27x3": 27 * 3}
print(f"   TFPT: 248 = {dims['E6']} + {dims['SU3']} + 2 x {dims['27x3']} = {dims['E6'] + dims['SU3'] + 2 * dims['27x3']}"
      " (Familien als SU(3)-Triplett in einem E8, SEAM.EQUIV.02)")
check("248 Zerlegung", dims["E6"] + dims["SU3"] + 2 * dims["27x3"] == 248)
print("   Mit einem E8 je Familie (E8 hoch 3) stünden die Familien nicht mehr in (27, 3) eines E8;")
print("   TFPTs Gravitation nutzt außerdem c = 8 = g_car + N_fam.")

print("\n3. Ist die 3D Wand beobachtbar?\n")
print("    c     Phase beim Verschieben der Wand, exp(2 pi i c sigma/8), sigma = 1 .. 4")
for c in (0.5, 1, 4, 8, 12, 16, 24):
    ph = [np.exp(2j * np.pi * c * s / 8) for s in range(1, 5)]
    trivial = all(abs(p - 1) < 1e-12 for p in ph)
    print(f"   {c:4}   " + "  ".join(f"{p.real:+.2f}{p.imag:+.2f}i" for p in ph) + ("   unsichtbar" if trivial else ""))
invisible = [c for c in np.arange(0.5, 25, 0.5) if all(abs(np.exp(2j * np.pi * c * s / 8) - 1) < 1e-12 for s in range(1, 9))]
print(f"   Unsichtbar für alle orientierten 4-Mannigfaltigkeiten genau bei c in 8Z: {[float(c) for c in invisible[:3]]} ...")
check("Wand unsichtbar genau bei c in 8Z", invisible[:3] == [8.0, 16.0, 24.0])
print("   Auf Spinmannigfaltigkeiten ist sigma durch 16 teilbar, dort genügte c in Z/2. Die GSO-Summe macht die")
print("   Naht aber spinstrukturfrei, also zählen alle orientierten Mannigfaltigkeiten, und dann ist c = 8 der")
print("   kleinste Wert, bei dem die Wand unsichtbar ist. Das ist dieselbe 16 wie in Z1.")
print("   Lesart: Die zwei Blätter von P1 sind entlang eines 3D Schnitts verklebt, der an der Naht endet. Trägt")
print("   der Schnitt die invertierbare E8 Phase, ist seine Lage unbeobachtbar (wie ein Dirac-String) genau")
print("   weil c = 8. F ist dann keine zusätzliche Wand, sondern der Verzweigungsschnitt der zwei Blätter.")

print("\nEntscheidung: c = 8 mit F. Die Operatorlesart (c = 24) widerspricht dem Zustand von H6 (K6) und")
print("TFPTs Einbettung der Familien in ein E8; c = 8 ist der kleinste Wert mit unsichtbarem Schnitt.")
print("\nZusammenfassung")
for n, c in OK:
    print(f"   {'bestanden' if c else 'NICHT bestanden'}: {n}")
Path(__file__).with_name("z5_ergebnis.json").write_text(json.dumps(
    {"checks": {n: c for n, c in OK}, "S_on_sectors": S_map, "T_on_sectors": T_map, "decision": "c = 8 mit F"}, indent=1))
assert all(c for _, c in OK)
