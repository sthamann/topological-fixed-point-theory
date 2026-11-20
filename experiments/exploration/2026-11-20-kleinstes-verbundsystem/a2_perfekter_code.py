#!/usr/bin/env python3
"""A2 -- Perfektion: jede Fehlerantwort wird genau einmal gebraucht.

Prinzip: Ein Quantenfreiheitsgrad wird gegen jeden einzelnen lokalen Fehler
geschuetzt, ohne ein einziges unbenutztes Syndrom (Quanten Hamming Schranke mit
Gleichheit). Gesucht: welche Zahlen das erzwingt. Exakte Ganzzahlrechnung.

Ausfuehren: python3 a2_perfekter_code.py
"""
import itertools

OK = True


def check(name, cond, detail=""):
    global OK
    OK &= bool(cond)
    print(f"  [{'ok' if cond else 'FEHLER'}] {name}{(' -- ' + detail) if detail else ''}")


# ---------------------------------------------------------------- 1. Diophantische Bedingung
print("1. Perfekte Codes fuer einen Fehler: D^(n-k) = 1 + (D^2 - 1) n")
sols = []
for D in range(2, 65):
    for n in range(2, 2001):
        rhs, m, pw = 1 + (D * D - 1) * n, 0, 1
        while pw < rhs:
            pw *= D
            m += 1
        if pw == rhs and 1 <= n - m:
            sols.append((D, n, n - m))
k1 = [s for s in sols if s[2] == 1]
check("genau eine Loesung mit einem logischen Freiheitsgrad (D bis 64, n bis 2000)",
      k1 == [(2, 5, 1)], f"k=1: {k1}")
print(f"  alle Loesungen mit n <= 100: {[s for s in sols if s[1] <= 100]}")
print("  Daraus ohne Wahl: D = 2, n = 5 Plaetze, n - k = 4 Pruefungen, 16 = 1 + 3*5 Syndrome.")

# ---------------------------------------------------------------- 2. Der Code
print("\n2. Der [[5,1,3]] Code")
SYM = "IXZY"  # Symbol = x + 2z, Produkt bis auf Phase = XOR
code = {c: i for i, c in enumerate(SYM)}
gens = [tuple(code[c] for c in "XZZXI"[-s:] + "XZZXI"[:-s]) for s in range(4)]


def anti(a, b):
    return sum(1 for x, y in zip(a, b) if x and y and x != y) % 2


def mul(a, b):
    return tuple(x ^ y for x, y in zip(a, b))


def wt(a):
    return sum(1 for x in a if x)


S = {(0,) * 5}
for g in gens:
    S |= {mul(s, g) for s in S}
ALL = list(itertools.product(range(4), repeat=5))
N = [p for p in ALL if not any(anti(p, g) for g in gens)]
check("vier kommutierende Erzeuger, |S| = 16, |N(S)| = 64",
      not any(anti(a, b) for a in gens for b in gens) and len(S) == 16 and len(N) == 64)
A = [sum(1 for s in S if wt(s) == w) for w in range(6)]
B = [sum(1 for s in N if wt(s) == w) for w in range(6)]
check("Gewichtszaehler A = (1,0,0,0,15,0), B = (1,0,0,30,15,18): Abstand 3",
      A == [1, 0, 0, 0, 15, 0] and B == [1, 0, 0, 30, 15, 18], f"A={A}, B={B}")
print("  B_1 = B_2 = 0: je zwei Plaetze tragen keine Information, je drei genuegen (3 + 2 Schwelle).")
ZZ = (code["Z"],) * 5
check("ZZZZZ liegt in N(S) ohne S: Paritaet ist das logische Z, der gerade Sektor hat 16 = ein Zustand je Syndrom",
      ZZ in N and ZZ not in S)


def syn(e):
    return tuple(anti(e, g) for g in gens)


single = {(j, c): syn(tuple(code[c] if i == j else 0 for i in range(5))) for j in range(5) for c in "XZY"}
check("die 15 Einzelfehler haben 15 verschiedene Syndrome ungleich 0 (perfekt)",
      len(set(single.values())) == 15 and (0, 0, 0, 0) not in single.values())
lines = [frozenset(single[(j, c)] for c in "XZY") for j in range(5)]
check("je Platz bilden die drei Syndrome eine Gerade; die 5 Geraden sind ein Spread von PG(3,2)",
      all(tuple(a ^ b for a, b in zip(single[(j, 'X')], single[(j, 'Z')])) == single[(j, 'Y')] for j in range(5))
      and len(frozenset().union(*lines)) == 15)


def matvec(M, v):
    return tuple(sum(M[i][k] * v[k] for k in range(4)) % 2 for i in range(4))


mats = [tuple(tuple(r[4 * i:4 * i + 4]) for i in range(4)) for r in itertools.product((0, 1), repeat=16)]
omega = [M for M in mats if all(matvec(M, v) in l and matvec(M, v) != v for l in lines for v in l)]
I4 = tuple(tuple(int(i == j) for j in range(4)) for i in range(4))


def mm(Aa, Bb):
    return tuple(tuple(sum(Aa[i][k] * Bb[k][j] for k in range(4)) % 2 for j in range(4)) for i in range(4))


check("genau zwei lineare Abbildungen drehen alle 5 Geraden zugleich (omega, omega^2): F_4 Struktur",
      len(omega) == 2 and all(mm(M, mm(M, M)) == I4 for M in omega), f"{len(omega)} Abbildungen")
print("  Die 5 Plaetze sind also die Punkte der projektiven Geraden ueber F_4 (desarguesscher Spread).")
alt = []
for bits in itertools.product((0, 1), repeat=6):
    Bm = [[0] * 4 for _ in range(4)]
    for (i, j), b in zip(itertools.combinations(range(4), 2), bits):
        Bm[i][j] = Bm[j][i] = b
    nondeg = len({matvec(Bm, v) for v in itertools.product((0, 1), repeat=4)}) == 16
    lag = all(sum(u[i] * Bm[i][j] * v[j] for i in range(4) for j in range(4)) % 2 == 0
              for l in lines for u in l for v in l)
    if nondeg and lag:
        alt.append(Bm)
check("genau 3 symplektische Formen machen den Spread lagrangesch: er ist ein vollstaendiger Satz von 5 "
      "unverzerrten Basen eines Zwei-Qubit-Registers", len(alt) == 3, f"{len(alt)} Formen")

# ---------------------------------------------------------------- 3. Symmetrie
print("\n3. Symmetrie des Codes auf den fuenf Plaetzen")
LC = list(itertools.permutations((1, 2, 3)))
perms = set()
for pi in itertools.permutations(range(5)):
    # je Platz eine lokale Clifford Klasse (Permutation von X,Z,Y); Platz fuer Platz pruefen
    def ok(choice):
        for g in gens:
            img = [0] * 5
            for j, sig in choice.items():
                img[pi[j]] = sig[g[j] - 1] if g[j] else 0
            # partielle Pruefung nur wenn vollstaendig
            if len(choice) == 5 and tuple(img) not in S:
                return False
        return True
    found = False
    for sigs in itertools.product(LC, repeat=5):
        if ok(dict(enumerate(sigs))):
            found = True
            break
    if found:
        perms.add(pi)
pairs = list(itertools.combinations(range(5), 2))
orb = {frozenset((p[a], p[b])) for p in perms for a, b in [pairs[0]]}
stab = [p for p in perms if {p[0], p[1]} == {0, 1}]
check("Permutationsanteil der Symmetriegruppe hat Ordnung 120 oder 60 und wirkt transitiv auf den 10 Paaren",
      len(perms) in (60, 120) and len(orb) == 10, f"Ordnung {len(perms)}, Bahn eines Paares {len(orb)}")
print(f"  Stabilisator einer 3 + 2 Teilung: Ordnung {len(stab)}. Die 10 Teilungen sind gleichwertig:")
print("  der Code liefert die Schwelle 3 + 2, die Markierung einer Teilung ist eine Symmetriebrechung.")

print("\nGESAMT:", "alle Pruefungen bestanden" if OK else "FEHLER")
