#!/usr/bin/env python3
"""B1 -- Die Klebegleichung: warum die selbstbeschreibende Algebra M_2(C) ist.

Fuer ein n-Niveau-System gibt es drei unabhaengige Forderungen an die Zahl m der
Plaetze:
  (i)   m = n + 1        vollstaendige unverzerrte Beschreibung (n + 1 Basen)
  (ii)  n^2 = 2^(m-1)    Operatoralgebra = gerade Cliffordalgebra der m Erzeuger
  (iii) m/4 + (n-1)/n = 2  Operator (x) Zustand hat Wurzelnorm 2
Dazu die Gitterfrage: wann ist D_m + A_(n-1), verklebt durch (Spinor, Fundamental),
gerade und unimodular? Exakte Rechnung mit Bruechen.

Ausfuehren: python3 b1_klebegleichung.py
"""
import itertools
from fractions import Fraction as Fr
from math import lcm

OK = True


def check(name, cond, detail=""):
    global OK
    OK &= bool(cond)
    print(f"  [{'ok' if cond else 'FEHLER'}] {name}{(' -- ' + detail) if detail else ''}")


print("1. Drei Forderungen, n und m bis 400")
R = range(1, 401)
c1 = {(n, m) for n in R for m in R if m == n + 1}
c2 = {(n, m) for n in R for m in R if m <= 40 and n * n == 2 ** (m - 1)}
c3 = {(n, m) for n in R for m in R if Fr(m, 4) + Fr(n - 1, n) == 2}
print(f"  (i) und (ii):   {sorted(c1 & c2)}")
print(f"  (i) und (iii):  {sorted(c1 & c3)}")
print(f"  (ii) und (iii): {sorted(c2 & c3)}")
check("alle drei zusammen: nur n = 4, m = 5; (i)+(iii) und (ii)+(iii) allein schon eindeutig",
      c1 & c2 & c3 == {(4, 5)} and c1 & c3 == {(4, 5)} and c2 & c3 == {(4, 5)})
print("  (i) + (iii) ist die quadratische Gleichung n^2 - 3n - 4 = 0 mit der positiven Wurzel 4.")
check("Standardform von M_k(C) hat Dimension n = k^2: n = 4 heisst k = 2",
      [k for k in range(1, 30) if k * k == 4] == [2])

print("\n2. Gitter D_m + A_(n-1) + Z (Spinor, Fundamental): gerade und unimodular?")


def glue_lattice(m, n):
    """Diskriminantenrechnung: Ordnung der Klebeklasse, Normen ihrer Vielfachen."""
    ord_s = 4 if m % 2 else 2
    idx = lcm(ord_s, n)
    det = Fr(4 * n, idx * idx)

    def qD(k):  # Norm der Klasse k*s in D_m^*/D_m modulo 2
        k %= ord_s
        if ord_s == 4:
            return {0: Fr(0), 1: Fr(m, 4), 2: Fr(1), 3: Fr(m, 4)}[k]
        return {0: Fr(0), 1: Fr(m, 4)}[k]

    def qA(k):  # Norm der Klasse k*f in A_(n-1)^*/A_(n-1) modulo 2
        k %= n
        return Fr(k * (n - k), n)

    even = all((qD(k) + qA(k)) % 2 == 0 for k in range(idx))
    return det, even


hits = [(m, n) for m in range(2, 25) for n in range(1, 25)
        if glue_lattice(m, n) == (Fr(1), True)]
print(f"  gerade unimodular fuer m bis 24, n bis 24: {hits}")
check("mit m = n + 1 bleibt nur (m, n) = (5, 4)", [h for h in hits if h[0] == h[1] + 1] == [(5, 4)])
check("ohne (i): genau n = 1 mit m = 0 mod 8 (D8+, D16+, D24+) und n = 4 mit m = 5 mod 8",
      set(hits) == {(m, 1) for m in (8, 16, 24)} | {(m, 4) for m in (5, 13, 21)})
print("  Die Gitterbedingung allein laesst also die Reihe D5, D13, D21 mit A3 zu;")
print("  erst m = n + 1 oder n^2 = 2^(m-1) waehlt D5 + A3 aus.")
check("mit (ii) n^2 = 2^(m-1) bleibt ebenfalls nur (5, 4)",
      [h for h in hits if h[1] ** 2 == 2 ** (h[0] - 1)] == [(5, 4)])

print("\n3. Das Gitter fuer (5, 4) ausgeschrieben: 240 Wurzeln nach Klebeklasse")
roots = []
for v in itertools.product((-1, 0, 1), repeat=8):
    if sum(a * a for a in v) == 2:
        roots.append(tuple(Fr(a) for a in v))
for s in itertools.product((Fr(1, 2), Fr(-1, 2)), repeat=8):
    if sum(1 for a in s if a < 0) % 2 == 0:
        roots.append(s)
check("240 Vektoren der Norm 2", len(roots) == 240 and all(sum(a * a for a in r) == 2 for r in roots))


def glue_class(r):
    """Klasse in Z_4: Koordinaten 0..4 tragen D5, 5..7 tragen A3 = D3."""
    if r[0].denominator == 1:
        return 2 * (int(sum(r[:5])) % 2)
    return 1 if sum(1 for a in r[:5] if a < 0) % 2 == 0 else 3


cnt = {k: sum(1 for r in roots if glue_class(r) == k) for k in range(4)}
check("Klassen 0, 1, 2, 3 tragen 52, 64, 60, 64 Wurzeln (40 + 12, 16 x 4, 10 x 6, 16 x 4)",
      cnt == {0: 52, 1: 64, 2: 60, 3: 64}, f"{cnt}")
s_norm, f_norm = Fr(5, 4), Fr(3, 4)
check("q(D5) + q(A3) = 5/4 + 3/4 = 2", s_norm + f_norm == 2)
print("  Klasse 1 ist (Operatoralgebra 16) x (Zustandsraum 4): die Modulwirkung selbst ist der Kleber.")

print("\nGESAMT:", "alle Pruefungen bestanden" if OK else "FEHLER")

raise SystemExit(0 if OK else 1)
