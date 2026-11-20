#!/usr/bin/env python3
"""A6 -- Fixpunkt der Dualitaet: eine Beschreibung, die ihre eigene Pruefung ist.

Prinzip: Ein binaerer Code C mit C = C^perp (selbstdual) und Gewichten in 4Z
(doppelt gerade, damit die Phasen der zugehoerigen Quantenzustaende konsistent
sind). Gesucht ist die kleinste Laenge und ob der Code dort eindeutig ist.

Ausfuehren: python3 a6_selbstdual.py
"""
import itertools
import math

OK = True


def check(name, cond, detail=""):
    global OK
    OK &= bool(cond)
    print(f"  [{'ok' if cond else 'FEHLER'}] {name}{(' -- ' + detail) if detail else ''}")


def type2_codes(n):
    """alle doppelt geraden selbstdualen Codes der Laenge n als Mengen von Bitmasken"""
    cand = [v for v in range(1, 1 << n) if bin(v).count("1") % 4 == 0]
    level = {frozenset([0])}
    for _ in range(n // 2):
        nxt = set()
        for code in level:
            for v in cand:
                if v in code or any(bin(v & c).count("1") % 2 for c in code):
                    continue
                new = code | {v ^ c for c in code}
                if all(bin(c).count("1") % 4 == 0 for c in new):
                    nxt.add(frozenset(new))
        level = nxt
        if not level:
            return []
    return list(level)


print("1. Doppelt gerade selbstduale Codes kleiner Laenge")
counts = {n: len(type2_codes(n)) for n in (2, 4, 6, 8)}
check("keine fuer n = 2, 4, 6; genau 30 fuer n = 8", counts == {2: 0, 4: 0, 6: 0, 8: 30}, f"{counts}")
codes = type2_codes(8)
c0 = min(codes, key=sorted)


def permute(code, p):
    return frozenset(sum(((c >> i) & 1) << p[i] for i in range(8)) for c in code)


aut = sum(1 for p in itertools.permutations(range(8)) if permute(c0, p) == c0)
check("alle 30 sind aequivalent: 8!/|Aut| = 30 mit |Aut| = 1344", math.factorial(8) // aut == 30 and aut == 1344,
      f"|Aut| = {aut}")
wd = sorted(bin(c).count("1") for c in c0)
check("Gewichte 0, 4 (14 mal), 8: der erweiterte Hammingcode [8,4,4]", wd == [0] + [4] * 14 + [8])

print("\n2. Konstruktion A")
kiss = 2 * 8 + sum(2 ** 4 for c in c0 if bin(c).count("1") == 4)
check("kuerzeste Vektoren von (C + 2Z^8)/sqrt(2): 16 + 14*16 = 240", kiss == 240, f"{kiss}")
print("  Selbstdualitaet macht das Gitter unimodular, doppelt gerade macht es gerade: E8.")
print("  30 = Anzahl der Codes ist zugleich die Coxeterzahl von E8; als Zahl nicht ueber der Basislinie")
print("  (30 = 2*3*5 = (1+1)(2+1)(4+1) aus der Massenformel).")

print("\nGESAMT:", "alle Pruefungen bestanden" if OK else "FEHLER")
