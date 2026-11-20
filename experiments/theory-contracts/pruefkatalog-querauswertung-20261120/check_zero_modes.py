"""NULLMODEN20261120: Nullmoden des Dolbeault Dirac Operators auf der markierten Sphäre.

Operator: B = sqrt(2) (dbar + dbar*) auf Lambda^{0,*} (x) (O(D) + K(D)), D = die vier Marken.
Sein Kern ist H^0 + H^1 der beiden Bündel. Auf der Sphäre ist O(D) = O(4) und K(D) = O(2).
Die Kohomologie wird exakt aus dem Cech Komplex der zwei Karten gezählt. Verglichen wird mit
den Realisierungen, die das Original selbst angibt. Kein Beweis, dass dieser Operator das
Randdatum B_Sigma des Originals ist.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SEC = ROOT/"documents/archive--tfpt-45--source_extracts--02_carrier_source/sections"
FAM = " ".join(next(SEC.glob("006-*.tex")).read_text().split())
HIG = " ".join(next(SEC.glob("007-*.tex")).read_text().split())
checks = []


def check(name, ok, detail=""):
    checks.append(bool(ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))


def cech(n, N=40):
    """(h0, h1) von O(n) auf der Sphäre: f0(z) - z^n f1(1/z) auf Monomen bis Grad N."""
    hit0 = set(range(0, N + 1))                 # z^a aus der Karte um 0
    hit1 = {n - b for b in range(0, N + 1)}     # w^b aus der Karte um unendlich
    target = set(range(n - N, N + 1))
    return len(hit0 & hit1), len(target - hit0 - hit1)


# --- 1. Nullmoden -------------------------------------------------------------------------
c = {n: cech(n) for n in (-6, -4, -1, 0, 1, 2, 4)}
check("H^0(O(D)) = 5, H^1 = 0", c[4] == (5, 0))
check("H^0(K(D)) = 3, H^1 = 0", c[2] == (3, 0))
check("der Operator hat genau 8 Nullmoden, alle einer Chiralität", c[4][0] + c[2][0] == 8 and c[4][1] + c[2][1] == 0)
check("Index = 2 mal Markenzahl", all(cech(n)[0] + cech(n - 2)[0] == 2 * n for n in range(2, 9)),
      "n Marken geben 2n Nullmoden; 8 nur bei n = 4")
check("Kontrolle: Serre Dualität h1(O(n)) = h0(O(-n-2))", all(cech(n)[1] == cech(-n - 2)[0] for n in range(-8, 9)))

# --- 2. Quotient nach z -> -z: Koordinate u = z^2 ---------------------------------------------
even = [k // 2 for k in range(5) if k % 2 == 0]      # z^(2j) = u^j
odd = [(k - 1) // 2 for k in range(5) if k % 2 == 1]  # z^(2j+1) = z u^j
check("gerader Teil des Trägers = H^0(O(2)) auf der Quotientensphäre", even == [0, 1, 2] and c[2][0] == 3)
check("ungerader Teil des Trägers = z mal H^0(O(1)) auf der Quotientensphäre", odd == [0, 1] and c[1][0] == 2)
check("damit ist der gerade Teil das symmetrische Quadrat des ungeraden", sorted(a + b for i, a in enumerate(odd) for b in odd[i:]) == even)

# --- 3. Original --------------------------------------------------------------------------
check("Original: Familienraum ist H^1 der vierfach punktierten Sphäre", r"F:=H^1(X_f^\circ,\mathbb C)" in FAM
      and r"X_f^\circ\cong \mathbb P^1\setminus\mu_4" in FAM)
check("Original: Vierteldrehung r(z) = iz und s(z) = 1/z", "r(z)=iz" in FAM and "s(z)=1/z" in FAM)
check("Original: schwaches Dublett ist H^0(S^2, O(1)) = C^2", r"H^0(S^2,\mathcal O(1))\cong\mathbb C^2" in HIG
      and r"L_2\cong\mathcal O(1)" in HIG)
check("Original: Farbklasse ist trivial, L_3 = O, kein leichtes Farbtriplett", r"L_3\cong\mathcal O" in HIG
      and "no light color triplet" in HIG)
check("Original: Determinantenklassen (c1(L_2), c1(L_3)) = (1, 0)", r"(c_1(L_2),c_1(L_3))=(1,0)" in HIG)

# --- 4. Vergleich -------------------------------------------------------------------------
check("Familien: 3 logarithmische Formen = Rang von H^1, wie im Original", c[2][0] == 4 - 1)
check("schwach: Grad 1 und Dimension 2, wie im Original", (1, c[1][0]) == (1, 2))
check("Farbe: hier Grad 2 mit drei Schnitten, im Original Grad 0", (2, c[2][0]) != (0, c[0][0]),
      "keine Übereinstimmung: das Original realisiert die Farbplätze nicht als Nullmoden")
# Determinantenwindung der Drehung auf den Trägerblöcken bei jeder einheitlichen Hebung (Gewichte k + a).
lifts = {a: (sum(k + a for k in (1, 3)), sum(k + a for k in (0, 2, 4))) for a in range(-6, 7)}
half = {a / 2: (sum(k + a / 2 for k in (1, 3)), sum(k + a / 2 for k in (0, 2, 4))) for a in range(-12, 13)}
check("keine einheitliche Hebung gibt die Determinantenwindungen (1, 0)", (1, 0) not in lifts.values() and (1, 0) not in half.values(),
      f"symmetrische Hebung a = -2 gibt {lifts[-2]}")

print(f"\n{sum(checks)} von {len(checks)} Prüfungen bestanden")
raise SystemExit(0 if all(checks) else 1)
