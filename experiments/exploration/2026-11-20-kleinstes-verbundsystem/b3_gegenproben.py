#!/usr/bin/env python3
"""B3 -- Gegenproben zu B1 und B2: falsche Varianten muessen scheitern.

Jede Gegenprobe ersetzt genau eine Zutat des Woerterbuchs oder der Klebegleichung
durch eine naheliegende falsche und prueft, dass die Aussage dann nicht mehr gilt.
Eine Gegenprobe "feuert", wenn die falsche Variante scheitert.

Ausfuehren: python3 b3_gegenproben.py
"""
import itertools
from fractions import Fraction as Fr
from math import lcm

import numpy as np

OK = True


def fires(name, cond, detail=""):
    global OK
    OK &= bool(cond)
    print(f"  [{'feuert' if cond else 'FEUERT NICHT'}] {name}{(' -- ' + detail) if detail else ''}")


_roots = []
for v in itertools.product(range(-1, 2), repeat=8):
    if sum(a * a for a in v) == 2 and sum(v) % 2 == 0:
        _roots.append(tuple(2 * a for a in v))
for y in itertools.product((0, -1), repeat=8):
    v = tuple(2 * a + 1 for a in y)
    if sum(a * a for a in v) == 8 and sum(v) % 4 == 0:
        _roots.append(v)
RD = set(_roots)
SIGMA_IDX = [4, 5, 0, 1, 2, 3, 6, 7]


def sig_vec(x):
    return tuple(x[i] for i in SIGMA_IDX)


SMALL = [complex(a, b) for a in range(-2, 3) for b in range(-2, 3) if a * a + b * b <= 4]


def in_L(x):
    par = {(int(c.real) + int(c.imag)) % 2 for c in x}
    s = sum(x)
    return len(par) == 1 and int(s.real) % 2 == 0 and int(s.imag) % 2 == 0


LROOTS = [np.array(x) for x in itertools.product(SMALL, repeat=4)
          if sum(int(c.real) ** 2 + int(c.imag) ** 2 for c in x) == 4 and in_L(x)]
BELL = [np.array(b, dtype=complex) for b in [(0, 1, 1, 0), (0, -1j, 1j, 0), (1, 0, 0, -1), (1, 0, 0, 1)]]
COMP = [np.array(b, dtype=complex) for b in np.eye(4)]


def chart(r, frame, scale, u0):
    out = []
    for k, b in enumerate(frame):
        z = (b.conj() @ r) / scale / (u0 if k == 3 else 1)
        out += [2 * z.real, 2 * z.imag]
    if any(abs(t - round(t)) > 1e-9 for t in out):
        return None
    return tuple(int(round(t)) for t in out)


def image(frame, scale, u0):
    return {chart(r, frame, scale, u0) for r in LROOTS}


print("Gegenproben zum Woerterbuch")
good = image(BELL, 2, 1j)
print(f"  richtige Karte (Paulikomponenten, Anker = i mal Eins): Bild = TFPT Wurzeln: {good == RD}")
OK &= good == RD
bad1 = image(BELL, 2, 1)
fires("G1 Ankerachse = Eins statt i mal Eins", bad1 != RD, f"gemeinsame Wurzeln {len(bad1 & RD)} von 240")
bad2 = image(COMP, 1, 1)
fires("G2 Rechenbasis statt Paulikomponenten", bad2 != RD, f"gemeinsame Wurzeln {len(bad2 & RD)} von 240")
U3 = 0.5 * np.array([[1 - 1j, -1 - 1j], [1 - 1j, 1 + 1j]])
one_sided = np.kron(U3, np.eye(2))
two_sided = np.kron(U3, U3.conj())
ok_two = all(chart(two_sided @ r, BELL, 2, 1j) == sig_vec(chart(r, BELL, 2, 1j)) for r in LROOTS)
print(f"  richtige Familiendrehung u (x) u-quer = sigma: {ok_two}")
OK &= ok_two
bad3 = sum(1 for r in LROOTS
           if chart(one_sided @ r, BELL, 2, 1j) not in (sig_vec(chart(r, BELL, 2, 1j)),
                                                         sig_vec(sig_vec(chart(r, BELL, 2, 1j)))))
fires("G3 einseitige Drehung u (x) 1 statt u (x) u-quer", bad3 > 0, f"{bad3} von 240 Wurzeln weichen ab")
om = np.array([1, 0, 0, 1], dtype=complex)
fires("G4 die einseitige Drehung laesst den Zustand Omega nicht fest",
      not np.allclose(one_sided @ om, om) and np.allclose(two_sided @ om, om))

print("\nGegenproben zur Klebegleichung")


def glue_ok(m, n):
    ord_s = 4 if m % 2 else 2
    idx = lcm(ord_s, n)

    def qD(k):
        k %= ord_s
        return ({0: Fr(0), 1: Fr(m, 4), 2: Fr(1), 3: Fr(m, 4)} if ord_s == 4 else {0: Fr(0), 1: Fr(m, 4)})[k]

    def qA(k):
        k %= n
        return Fr(k * (n - k), n)

    return Fr(4 * n, idx * idx) == 1 and all((qD(k) + qA(k)) % 2 == 0 for k in range(idx))


print(f"  richtige Wahl D5 + A3: gerade unimodular: {glue_ok(5, 4)}")
OK &= glue_ok(5, 4)
for m, n, why in ((4, 5, "Rollen vertauscht"), (6, 4, "ein Platz mehr"), (5, 3, "Zustandsraum 3"),
                  (3, 2, "ein Qubit ohne Spiegel: erfuellt (i) und (ii), aber nicht (iii)")):
    fires(f"G5 D{m} + A{n - 1} ({why})", not glue_ok(m, n))

print("\nGESAMT:", "alle Gegenproben feuern" if OK else "FEHLER")
raise SystemExit(0 if OK else 1)
