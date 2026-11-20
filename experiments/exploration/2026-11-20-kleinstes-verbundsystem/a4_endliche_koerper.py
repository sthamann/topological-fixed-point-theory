#!/usr/bin/env python3
"""A4 -- Arithmetik ueber dem kleinsten Koerper, mit Zahlenspielkontrolle.

Prinzip: Die Welt ist ueber F_2 definiert; Plaetze sind Punkte, Zaehlungen sind
Frobenius Fixpunkte. Geprueft werden zwei Strukturaussagen (5 = 3 + 2) und die
Punktzahlen 5 und 41 gegen eine Zufallsbasislinie. Exakte Ganzzahlrechnung.

Ausfuehren: python3 a4_endliche_koerper.py
"""
import itertools
from fractions import Fraction as Fr

OK = True


def check(name, cond, detail=""):
    global OK
    OK &= bool(cond)
    print(f"  [{'ok' if cond else 'FEHLER'}] {name}{(' -- ' + detail) if detail else ''}")


print("1. Projektive Gerade: Punkte ueber F_4 unter Frobenius")
# F_4 = {0,1,w,w+1} als Paare (a,b) = a + b w, w^2 = w + 1; Frobenius x -> x^2
F4 = [(a, b) for a in (0, 1) for b in (0, 1)]


def f4mul(x, y):
    a, b = x
    c, d = y
    return ((a * c + b * d) % 2, (a * d + b * c + b * d) % 2)


frob = {x: f4mul(x, x) for x in F4}
fixed = [x for x in F4 if frob[x] == x]
check("P^1(F_4) hat 5 Punkte: 3 rationale (0, 1, unendlich) und ein konjugiertes Paar",
      len(F4) + 1 == 5 and len(fixed) + 1 == 3 and len(F4) - len(fixed) == 2)

print("\n2. Elliptische Kurven ueber F_2")
curves = {}
for a1, a2, a3, a4, a6 in itertools.product((0, 1), repeat=5):
    b2, b4, b6 = a1 * a1 + 4 * a2, 2 * a4 + a1 * a3, a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    disc = -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6
    if disc % 2 == 0:
        continue
    pts = 1 + sum(1 for x in (0, 1) for y in (0, 1)
                  if (y * y + a1 * x * y + a3 * y - x ** 3 - a2 * x * x - a4 * x - a6) % 2 == 0)
    curves.setdefault(3 - pts, []).append((a1, a2, a3, a4, a6))
traces = sorted(curves)
ss = [t for t in traces if t % 2 == 0]
check("5 Isogenieklassen (Spuren -2..2): 3 supersingulaere + 2 gewoehnliche",
      traces == [-2, -1, 0, 1, 2] and len(ss) == 3, f"Spuren {traces}, supersingulaer {ss}")

print("\n3. Punktzahlen ueber F_(2^n), n = 1..8, gegen die Zielzahlen")
TARGETS = {3, 4, 5, 8, 16, 30, 41, 68, 240, 248}
hits, expected = [], Fr(0)
for a in traces:
    s = [2, a]
    for n in range(2, 9):
        s.append(a * s[-1] - 2 * s[-2])
    counts = [2 ** n + 1 - s[n] for n in range(1, 9)]
    print(f"  Spur {a:+d}: {counts}")
    for n, c in enumerate(counts, 1):
        q = 2 ** n
        lo = q + 1 - int((4 * q) ** 0.5)
        hi = q + 1 + int((4 * q) ** 0.5)
        expected += Fr(len([t for t in TARGETS if lo <= t <= hi]), hi - lo + 1)
        if c in TARGETS:
            hits.append((a, n, c))
print(f"  Treffer: {hits}")
print(f"  Anzahl Treffer {len(hits)}, Erwartung bei gleichverteilter Zahl im Hasse Intervall {float(expected):.2f}")
check("Treffer 5 = #E(F_4) und 41 = #E(F_32) fuer die Kurve mit Frobenius 1+i vorhanden",
      (2, 2, 5) in hits and (2, 5, 41) in hits)
big = [h for h in hits if h[2] > 16]
exp_big = Fr(0)
for a in traces:
    for n in range(1, 9):
        q = 2 ** n
        lo, hi = q + 1 - int((4 * q) ** 0.5), q + 1 + int((4 * q) ** 0.5)
        exp_big += Fr(len([t for t in TARGETS if t > 16 and lo <= t <= hi]), hi - lo + 1)
print(f"  Nur Ziele > 16: Treffer {len(big)} ({big}), Erwartung {float(exp_big):.2f}")
check("Zahlentreffer liegen nicht ueber der Basislinie (Faktor unter 3)",
      len(hits) < 3 * float(expected), f"{len(hits)} gegen {float(expected):.2f}")

print("\n4. Die 68")
q = 4
check("68 = q(q^2+1) = (q+1) + (q^3-1) bei q = 4, zugleich Index von PSL(2,4) in PSL(2,16)",
      q * (q * q + 1) == 68 == (q + 1) + (q ** 3 - 1) and (16 * 255) // 60 == 68)
small = sorted({qq * (qq * qq + 1) for qq in (2, 3, 4, 5, 7, 8, 9)})
print(f"  Dieselbe Formel fuer q = 2,3,4,5,7,8,9: {small}. Ein Treffer auf 7 Werte unter 800,")
print("  dazu ist 5 + 63 eine Identitaet der Formel. Als Zahl nicht ueber der Basislinie.")

print("\nGESAMT:", "alle Pruefungen bestanden" if OK else "FEHLER")
