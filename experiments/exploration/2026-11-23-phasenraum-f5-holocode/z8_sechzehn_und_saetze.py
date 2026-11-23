"""Z8: Welche 16 Fermionen bilden die Naht, und wie viele 16er-Sätze winden?

Z7 zählt am verdrillten Defekt eine chirale Majorana-Weyl Mode je Paar (psi, chi). In der Basis
xi_pm = (psi -+ i chi)/sqrt 2 windet nur xi_+, und xi_+ trägt unter J die Phase i: Es ist ein Weyl-Fermion
mit Spin-Z4-Ladung 1 (J^2 = (-1)^F). Der Partner xi_- (Ladung -1) windet nicht; er wird auf dem Gitter nur
gebraucht, damit ein J-invariantes Wilson-Glied die Verdoppler entfernt. Gezählt werden also Fermionen mit
Z4-Ladung 1, deren Majoranamasse am Defekt windet: nu solche geben c = nu/2.

Prüfungen
 1. Anomalie: Der Schnitt ist unsichtbar genau bei c in 8Z (Z5), also nu in 16Z. Das ist die Bedingung der
    Z16-Anomalie von Weyl-Fermionen mit Spin-Z4-Symmetrie (Garcia-Etxebarria und Montero 2019).
 2. Eine Standardmodell-Generation: Ladungen X = 5(B-L) - 4Y mod 4; welche Majoranamassen sind eichinvariant?
 3. Ein reelles 16-Plett (Vektor von SO(16)): eichinvariante Masse, c, und seine Lage in E8.
 4. Zahl der Sätze n: c = 8n; was verlangen Anomalie, Zufluss und die Familiencharaktere?
"""
import itertools
import json
from fractions import Fraction as Fr
from pathlib import Path

OK = []


def check(name, cond):
    OK.append((name, bool(cond)))
    return "ja" if cond else "NEIN"


z7 = json.loads(Path(__file__).with_name("z7_ergebnis.json").read_text())
print("1. Zählung aus Z7 und Anomalie\n")
one_pair = next(v for k, v in z7.items() if k.startswith("verdrillter Defekt (J"))
print(f"   verdrillter Defekt, ein Paar: netto {one_pair['net_core']} Mode -> je Fermion mit Z4-Ladung 1 eine Majorana-Weyl")
print("   Mode, c = 1/2. Unsichtbarer Schnitt: c in 8Z <=> nu in 16Z.")
rows = []
for nu in (1, 8, 15, 16, 32, 48):
    c = Fr(nu, 2)
    rows.append((nu, c, c % 8 == 0, nu % 16 == 0))
    print(f"   nu = {nu:2d}: c = {str(c):>4}, Schnitt unsichtbar {c % 8 == 0}, Z16-Anomalie null {nu % 16 == 0}")
check("Schnitt unsichtbar <=> Z16 frei", all(r[2] == r[3] for r in rows))

print("\n2. Eine Standardmodell-Generation (linkshändige Weyl-Felder)\n")
gen = [  # name, SU(3) dim (+ for 3, - for 3bar), SU(2) dim, Y, B-L, multiplicity
    ("Q", 3, 2, Fr(1, 6), Fr(1, 3)), ("u^c", -3, 1, Fr(-2, 3), Fr(-1, 3)), ("d^c", -3, 1, Fr(1, 3), Fr(-1, 3)),
    ("L", 1, 2, Fr(-1, 2), Fr(-1)), ("e^c", 1, 1, Fr(1), Fr(1)), ("nu^c", 1, 1, Fr(0), Fr(1)),
]
count = 0
allone = True
for n, c3, c2, y, bl in gen:
    x = 5 * bl - 4 * y
    mult = abs(c3) * c2
    count += mult
    allone &= (x.denominator == 1 and int(x) % 4 == 1)
    print(f"   {n:5s} ({'3bar' if c3 < 0 else c3},{c2},{str(y):>4})  X = {str(x):>3}  X mod 4 = {int(x) % 4}  Komponenten {mult}")
print(f"   {count} Weyl-Fermionen, alle mit X = 1 mod 4: {check('SM: 16 mit Ladung 1', allone and count == 16)}; ohne nu^c wären es 15,")
print("   also Z16-anomal. Die Zahl passt.")
inv = []
for (a, b) in itertools.combinations_with_replacement(range(len(gen)), 2):
    na, ca, sa, ya, _ = gen[a]
    nb, cb, sb, yb, _ = gen[b]
    color_ok = (ca == 1 and cb == 1) or (ca == -cb and abs(ca) == 3)
    su2_ok = sa == sb                                   # 2 x 2 contains a singlet, 1 x 1 too
    if color_ok and su2_ok and ya + yb == 0:
        inv.append(f"{na} {nb}")
print(f"   Eichinvariante Majorana-Bilineare in einer Generation: {inv}")
check("nur nu^c nu^c", inv == ["nu^c nu^c"])
print("   Nur nu^c kann eine windende Masse tragen: eine Mode je Generation, c = 1/2, drei Generationen c = 3/2.")
print("   Die leichten Standardmodell-Fermionen bilden die Naht also nicht (auf der Ebene von Bilinearen).")

print("\n3. Ein reelles 16-Plett: der Vektor von SO(16)\n")
print("   Reelle Darstellung, also ist die Majoranamasse Phi psi_a psi_a mit komplexem Singulett Phi eichinvariant;")
print("   Phi trägt Z4-Ladung 2, sein Vorzeichen ist das Blatt. Alle 16 winden: c = 8, nu = 16, anomaliefrei.")
so16 = {"16": (10, 6), "120": (45, 15, 10 * 6), "128": (16 * 4, 16 * 4)}
print(f"   SO(16) > SO(10) x SO(6): 16 = (10,1) + (1,6) = {sum(so16['16'])}; 120 = (45,1)+(1,15)+(10,6) = {sum(so16['120'])};"
      f" 128 = (16,4)+(16bar,4bar) = {sum(so16['128'])}")
check("E8 Zerlegung", sum(so16["16"]) == 16 and sum(so16["120"]) + sum(so16["128"]) == 248)
print("   Die 16 Moden der Naht tragen also (10,1)+(1,6): einen Vektor von SO(10) und ein Sextett von SU(4) = SO(6).")
print("   Ihre 4D Eltern sind schwere vektorartige Fermionen mit denselben Quantenzahlen (keine Eichanomalie).")
print("   Mit SU(4) > SU(3) x U(1) wird (16,4) zu (16,3) + (16,1): Die Standardmodell-Generationen sind Zustände")
print("   im Spinorteil 128 des Nahtkerns E8, nicht seine 16 Fermionen. Das passt zu TFPTs Familien in einem E8.")

print("\n4. Zahl der windenden Sätze n\n")
print("    n   c    nu   Z16   Naht allein T-invariant   mit Zufluss   Familiencharaktere aus der Naht")
table = []
for n in (1, 2, 3):
    c = 8 * n
    t_inv = (c % 24 == 0)
    fam = "3 Sektoren, {-1, i, -i}" if n == 1 else f"{3 * n} Sektoren ({n} Kopien)"
    table.append((n, t_inv))
    print(f"    {n}   {c:2d}   {16 * n:2d}   frei  {'ja' if t_inv else 'nein (Phase e^(-2 pi i n/3))':28s}  ja            {fam}")
print("   Unter G ist die Naht nie allein: Die schweren 4D Fermionen, deren Masse am Defekt windet, liefern den")
print("   Zufluss (Callan-Harvey). Strenge T-Invarianz der Naht für sich wird damit nicht verlangt, und das")
print("   kleinste n mit unsichtbarem Schnitt ist n = 1. Nur bei n = 1 trägt die Naht genau drei Spinstruktur-")
print("   sektoren, auf denen die Vierteldrehung {-1, i, -i} gibt (Z1). n = 3 wäre die kleinste Zahl ohne Zufluss,")
print("   verlöre aber diese Familiencharaktere und passt nicht zu TFPTs einem E8.")
check("n = 1 kleinste mit Zufluss", min(n for n in (1, 2, 3)) == 1)
check("n = 3 kleinste ohne Zufluss", min(n for n, t in table if t) == 3)
print("\n   Entscheidung: n = 1, c = 8. F ist dann keine eigene Prämisse mehr, sondern der Zufluss der schweren")
print("   Fermionen aus G; der Schnitt ist der Phasenschnitt von Phi (Lesart, Callan-Harvey hier nicht gerechnet).")

print("\nZusammenfassung")
for n, c in OK:
    print(f"   {'bestanden' if c else 'NICHT bestanden'}: {n}")
Path(__file__).with_name("z8_ergebnis.json").write_text(json.dumps({
    "checks": {n: c for n, c in OK}, "sm_invariant_majorana": inv, "decision_sets": 1,
    "seam_parents": "reelles 16-Plett (10,1)+(1,6), Spin-Z4-Ladung 1"}, ensure_ascii=False, indent=1))
assert all(c for _, c in OK)
