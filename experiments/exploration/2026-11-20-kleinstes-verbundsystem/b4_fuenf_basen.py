#!/usr/bin/env python3
"""B4 -- Die fuenf Plaetze als fuenf unverzerrte Basen; 3 + 2 als Produkt / verschraenkt.

In der Karte aus B2 (Standardform von M_2(C)) ist die Familiendrehung sigma lokal:
u (x) u-quer. Frage: Wie wirkt sigma auf die 15 Messkontexte (Geraden der Pauli
Geometrie), auf die sechs vollstaendigen Saetze unverzerrter Basen (Spreads) und
auf die fuenf 48er Pakete des Korpus (sigma Bahnen der 15 Klassen zu je 16 Wurzeln)?

Exakt bis auf die Eigenwertpruefung der Paulis (Betraege 0 oder 1, Toleranz 1e-9).

Ausfuehren: python3 b4_fuenf_basen.py
"""
import itertools
from collections import Counter

import numpy as np

OK = True


def check(name, cond, detail=""):
    global OK
    OK &= bool(cond)
    print(f"  [{'ok' if cond else 'FEHLER'}] {name}{(' -- ' + detail) if detail else ''}")


I2 = np.eye(2, dtype=complex)
X2 = np.array([[0, 1], [1, 0]], dtype=complex)
Y2 = np.array([[0, -1j], [1j, 0]])
Z2 = np.diag([1, -1]).astype(complex)
ONE = {(0, 0): I2, (1, 0): X2, (1, 1): Y2, (0, 1): Z2}
NAME = {(0, 0): "I", (1, 0): "X", (1, 1): "Y", (0, 1): "Z"}
POINTS = [p for p in itertools.product((0, 1), repeat=4) if any(p)]


def mat(p):
    return np.kron(ONE[p[:2]], ONE[p[2:]])


def nm(p):
    return NAME[p[:2]] + NAME[p[2:]]


def symp(p, q):
    return (p[0] * q[1] + p[1] * q[0] + p[2] * q[3] + p[3] * q[2]) % 2


def add(p, q):
    return tuple((a + b) % 2 for a, b in zip(p, q))


LINES = sorted({frozenset((p, q, add(p, q))) for p in POINTS for q in POINTS
                if p != q and symp(p, q) == 0})
SPREADS = [s for s in itertools.combinations(LINES, 5) if len(frozenset().union(*s)) == 15]


def is_product(line):
    return any(p[:2] == (0, 0) or p[2:] == (0, 0) for p in line)


check("15 Kontexte: 9 Produktbasen, 6 verschraenkte; 6 Spreads", (len(LINES), sum(map(is_product, LINES)), len(SPREADS)) == (15, 9, 6))

# sigma = u (x) u-quer auf den Paulis
U3 = 0.5 * np.array([[1 - 1j, -1 - 1j], [1 - 1j, 1 + 1j]])
SIG = np.kron(U3, U3.conj())


def sig_point(p):
    m = SIG @ mat(p) @ SIG.conj().T
    for q in POINTS:
        if np.allclose(m, mat(q)) or np.allclose(m, -mat(q)):
            return q
    raise AssertionError


def sig_line(line):
    return frozenset(sig_point(p) for p in line)


print("\n1. sigma auf den 15 Kontexten")
fixed = [ln for ln in LINES if sig_line(ln) == ln]
check("genau 3 feste Kontexte, alle verschraenkt; die anderen 12 in 4 freien Dreierbahnen",
      len(fixed) == 3 and not any(map(is_product, fixed)))
for ln in fixed:
    print("    fest:", sorted(nm(p) for p in ln))
bell = frozenset(p for p in POINTS if nm(p) in ("XX", "YY", "ZZ"))
check("darunter die Bellbasis {XX, YY, ZZ}", bell in fixed)
orbs, left = [], set(LINES)
while left:
    ln = left.pop()
    o = {ln, sig_line(ln), sig_line(sig_line(ln))}
    left -= o
    orbs.append(o)
free = [o for o in orbs if len(o) == 3]
kinds = Counter("Produkt" if all(map(is_product, o)) else "verschraenkt" for o in free)
check("freie Bahnen: 3 aus Produktbasen, 1 aus verschraenkten", kinds == {"Produkt": 3, "verschraenkt": 1},
      f"{dict(kinds)}")
print("  Fuenf Pakete zu je drei Kontexten: 3 Produktbahnen + (1 verschraenkte Bahn, 1 festes Tripel).")

print("\n2. sigma auf den sechs Saetzen von fuenf unverzerrten Basen")
inv = [s for s in SPREADS if {sig_line(ln) for ln in s} == set(s)]
print(f"  sigma-invariante Saetze: {len(inv)} von 6")
pattern_ok = True
for s in inv:
    cyc = [ln for ln in s if sig_line(ln) != ln]
    fix = [ln for ln in s if sig_line(ln) == ln]
    pattern_ok &= len(cyc) == 3 and len(fix) == 2 and all(map(is_product, cyc)) and not any(map(is_product, fix))
    print("    gedreht:", [sorted(nm(p) for p in ln) for ln in cyc])
    print("    fest:   ", [sorted(nm(p) for p in ln) for ln in fix])
check("in jedem sigma-invarianten Satz dreht sigma die 3 Produktbasen und laesst die 2 verschraenkten fest",
      len(inv) >= 1 and pattern_ok)

print("\n3. Abgleich mit den Wurzeln: 15 Klassen zu 16 Wurzeln, Klasse = Kontext")
SMALL = [complex(a, b) for a in range(-2, 3) for b in range(-2, 3) if a * a + b * b <= 4]


def in_L(x):
    par = {(int(c.real) + int(c.imag)) % 2 for c in x}
    s = sum(x)
    return len(par) == 1 and int(s.real) % 2 == 0 and int(s.imag) % 2 == 0


LROOTS = [np.array(x) for x in itertools.product(SMALL, repeat=4)
          if sum(int(c.real) ** 2 + int(c.imag) ** 2 for c in x) == 4 and in_L(x)]
MATS = {p: mat(p) for p in POINTS}


def context(r):
    st = frozenset(p for p in POINTS if abs(abs(r.conj() @ MATS[p] @ r) - 4) < 1e-9)
    assert st in LINES
    return st


ctx = [context(r) for r in LROOTS]
census = Counter(ctx)
check("jede Wurzel ist ein Stabilisatorzustand; 15 Kontexte zu je 16 Wurzeln", sorted(census.values()) == [16] * 15)
ok_eq = all(context(SIG @ r) == sig_line(c) for r, c in zip(LROOTS, ctx))
check("sigma auf den Wurzeln und sigma auf den Kontexten stimmen ueberein", ok_eq)
om_root = np.array([1 + 1j, 0, 0, 1 + 1j])
check("die Wurzeln der Zustandsgeraden C Omega liegen im Kontext Bellbasis (die Ankerklasse)",
      context(om_root) == bell)
packs = sorted(sum(census[ln] for ln in o) for o in orbs if len(o) == 3) + [sum(census[ln] for ln in fixed)]
check("die fuenf Pakete tragen je 48 Wurzeln", packs == [48] * 5, f"{packs}")


def gi(c):
    return int(round(c.real)), int(round(c.imag))


def same_class(x, y):
    """x - y in (1+i) L"""
    d = (x - y) / (1 + 1j)
    return all(abs(c.real - round(c.real)) < 1e-9 and abs(c.imag - round(c.imag)) < 1e-9 for c in d) and in_L(d)


def hbar(x, y):
    """abgestiegene symplektische Form: h(x, y) / 2 modulo (1 + i)"""
    a, b = gi(np.vdot(x, y) / 2)
    return (a + b) % 2


cls_ok = all(same_class(x, y) == (cx == cy)
             for (x, cx), (y, cy) in itertools.combinations(list(zip(LROOTS, ctx)), 2))
check("Klasse in L/(1+i)L = Kontext, fuer alle 28680 Wurzelpaare", cls_ok)
reps = {}
for r, c in zip(LROOTS, ctx):
    reps.setdefault(c, r)
dual_ok = all((hbar(reps[a], reps[b]) == 0) == bool(a & b) for a, b in itertools.combinations(LINES, 2))
check("hbar = 0 genau dann, wenn zwei Kontexte einen Pauli teilen (alle 105 Paare)", dual_ok)
print("  Also ist ein hbar-Ovoid von fuenf Klassen genau ein Satz von fuenf unverzerrten Basen.")
print("  Der Ledger fuehrt die fuenf Traegerplaetze als hbar-Ovoid der sigma-erhaltenen Form q*")
print("  (NORMALFORM.CFIN.01); sie sind damit einer der drei sigma-invarianten Saetze oben.")

print("\nGESAMT:", "alle Pruefungen bestanden" if OK else "FEHLER")
raise SystemExit(0 if OK else 1)
