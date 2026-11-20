#!/usr/bin/env python3
"""B2 -- Woerterbuchtest gegen die TFPT Karte (v634, v638, v689).

Frage 1: Ist die TFPT Karte (E8 in doppelten Koordinaten, J = Drehung in den
Paaren, sigma = Dreierzyklus der Familienpaare, Ankerpaar {6,7}) dieselbe Struktur
wie die Standardform von M_2(C) in Paulikomponenten?
Frage 2: Wo liegt D5 + A3 relativ zu J und sigma? Erzwingt die Familiensymmetrie
sigma die Lage und damit die markierte Teilung 3 + 2 der fuenf D5 Achsen?

Wurzeln, J und sigma sind woertlich aus verification/v634_st31_structure.py
uebernommen. Exakte Ganzzahlrechnung.

Ausfuehren: python3 b2_woerterbuch.py
"""
import itertools
from collections import Counter

import numpy as np

OK = True


def check(name, cond, detail=""):
    global OK
    OK &= bool(cond)
    print(f"  [{'ok' if cond else 'FEHLER'}] {name}{(' -- ' + detail) if detail else ''}")


# ---------------------------------------------------------------- TFPT Karte (v634)
_roots = []
for v in itertools.product(range(-1, 2), repeat=8):
    if sum(a * a for a in v) == 2 and sum(v) % 2 == 0:
        _roots.append(tuple(2 * a for a in v))
for y in itertools.product((0, -1), repeat=8):
    v = tuple(2 * a + 1 for a in y)
    if sum(a * a for a in v) == 8 and sum(v) % 4 == 0:
        _roots.append(v)
RD = sorted(_roots)
IDX = {r: i for i, r in enumerate(RD)}
A = np.array(RD, dtype=np.int64)
GRAM = (A @ A.T) // 4  # Wurzelnorm 2
SIGMA_IDX = [4, 5, 0, 1, 2, 3, 6, 7]


def J_vec(x):
    out = [0] * 8
    for k in range(4):
        out[2 * k], out[2 * k + 1] = -x[2 * k + 1], x[2 * k]
    return tuple(out)


def sig_vec(x):
    return tuple(x[i] for i in SIGMA_IDX)


JP = [IDX[J_vec(r)] for r in RD]
SP = [IDX[sig_vec(r)] for r in RD]

# ---------------------------------------------------------------- 1. Woerterbuch
print("1. Standardform von M_2(C) gegen die TFPT Karte")
PHI = 1 + 1j
SMALL = [complex(a, b) for a in range(-2, 3) for b in range(-2, 3) if a * a + b * b <= 4]


def in_L(x):
    par = {(int(c.real) + int(c.imag)) % 2 for c in x}
    s = sum(x)
    return len(par) == 1 and int(s.real) % 2 == 0 and int(s.imag) % 2 == 0


LROOTS = [np.array(x) for x in itertools.product(SMALL, repeat=4)
          if sum(int(c.real) ** 2 + int(c.imag) ** 2 for c in x) == 4 and in_L(x)]
# Paulikomponenten: a = c_0 1 + c_1 X + c_2 Y + c_3 Z, Vektor a Omega mit Omega = (1,0,0,1)
PAULI_OMEGA = [np.array(b, dtype=complex) for b in
               [(1, 0, 0, 1), (0, 1, 1, 0), (0, -1j, 1j, 0), (1, 0, 0, -1)]]
U3 = 0.5 * np.array([[1 - 1j, -1 - 1j], [1 - 1j, 1 + 1j]])  # Clifford der Ordnung 3
X2 = np.array([[0, 1], [1, 0]], dtype=complex)
Y2 = np.array([[0, -1j], [1j, 0]])
Z2 = np.diag([1, -1]).astype(complex)
cyc = [np.allclose(U3 @ a @ U3.conj().T, b) for a, b in ((X2, Y2), (Y2, Z2), (Z2, X2))]
check("u der Ordnung 3 dreht X nach Y nach Z nach X", all(cyc) and np.allclose(U3 @ U3 @ U3, -np.eye(2)))
SIG2Q = np.kron(U3, U3.conj())


def to_tfpt(r, u0):
    """doppelte reelle Koordinaten: Familienpaare = X, Y, Z Komponente, Ankerpaar = Einskomponente mal u0"""
    c = [(b.conj() @ r) / 2 for b in PAULI_OMEGA]
    c0 = c[0] / u0
    out = []
    for z in (c[1], c[2], c[3], c0):
        out += [2 * z.real, 2 * z.imag]
    assert all(abs(t - round(t)) < 1e-9 for t in out)
    return tuple(int(round(t)) for t in out)


found = None
for u0, name in ((1, "1"), (1j, "i")):
    img = {to_tfpt(r, u0) for r in LROOTS}
    if img == set(RD):
        found = (u0, name)
check("Wurzeln der Standardform = die 240 TFPT Wurzeln, Koordinate fuer Koordinate",
      found is not None, f"Ankerachse = {found[1] if found else '-'} mal Eins")
u0 = found[0]
jo = all(to_tfpt(1j * r, u0) == J_vec(to_tfpt(r, u0)) for r in LROOTS)
so = [all(to_tfpt(SIG2Q @ r, u0) == f(to_tfpt(r, u0)) for r in LROOTS)
      for f in (sig_vec, lambda x: sig_vec(sig_vec(x)))]
check("Multiplikation mit i ist das TFPT J", jo)
check("Konjugation mit u (X -> Y -> Z, Eins fest) ist das TFPT sigma oder sein Inverses", any(so),
      "sigma" if so[0] else "sigma^-1")
print("  Also: Familienpaare = X, Y, Z Komponenten, Ankerpaar = Einskomponente (der Zustand Omega),")
print("  J = Phase i, sigma = Clifforddrehung der drei Pauliachsen.")
print("  Reelle Form: die geraden Koordinaten sind (X, Y, Z, i mal Eins), also i mal die Quaternionen.")

# ---------------------------------------------------------------- 2. Lage von D5 + A3
print("\n2. Alle A3 Teilsysteme von E8 und ihre Lage zu J und sigma")
a3 = set()
for b in range(240):
    nb = [a for a in range(240) if GRAM[a, b] == -1]
    for al, ga in itertools.combinations(nb, 2):
        if GRAM[al, ga] == 0:
            vs = [A[al], A[b], A[ga], A[al] + A[b], A[b] + A[ga], A[al] + A[b] + A[ga]]
            a3.add(frozenset(IDX[tuple(int(t) for t in s * v)] for v in vs for s in (1, -1)))
check("7560 A3 Teilsysteme", len(a3) == 7560, f"{len(a3)}")


def axes(root_ids, ambient_ids):
    """Achsen (Norm 1) als doppelte Vektoren r1 + r2, ganzzahlig gepaart mit allen Wurzeln"""
    out = set()
    for i, j in itertools.combinations(root_ids, 2):
        if GRAM[i, j] == 0:
            w = A[i] + A[j]
            if all((int(w @ A[k]) // 4) % 2 == 0 for k in ambient_ids):
                out.add(tuple(int(t) for t in w))
    return out


def J2(w):
    return J_vec(w)


def classify(S):
    S = sorted(S)
    d5 = [k for k in range(240) if all(GRAM[k, s] == 0 for s in S)]
    ax3, ax5 = axes(S, S), axes(d5, d5)
    partner = sum(1 for g in ax3 if J2(g) in ax5) // 2
    lines5 = sum(1 for e in ax5 if J2(e) in ax5) // 4
    lines3 = sum(1 for g in ax3 if J2(g) in ax3) // 4
    return len(d5), len(ax3), len(ax5), partner, lines5, lines3


types = Counter()
sig_inv = []
for S in a3:
    t = classify(S)
    types[t[3:]] += 1
    assert t[:3] == (40, 6, 10)
    if frozenset(SP[i] for i in S) == S:
        sig_inv.append((S, t))
print("  Typ (A3 Achsen mit J Partner in D5, komplexe Geraden in D5, komplexe Geraden in A3): Anzahl")
for k, v in sorted(types.items()):
    print(f"    {k}: {v}")
check("jede D5 Lage hat 10 Achsen, jede A3 Lage 6", True)
tsig = Counter(t[3:] for _, t in sig_inv)
print(f"  sigma-invariante A3: {len(sig_inv)}, Typen {dict(tsig)}")
block = [(S, t) for S, t in sig_inv if t[3:] == (3, 1, 0)]
check("sigma-invariante Lagen vom Blocktyp (3 Partnerachsen + 1 komplexe Gerade) existieren", len(block) >= 1,
      f"{len(block)}")
std = frozenset(IDX[tuple(2 * s * (int(k == i) + t * int(k == j)) for k in range(8))]
                for i, j in ((0, 2), (0, 4), (2, 4)) for s in (1, -1) for t in (1, -1))
check("darunter die Lage A3 = gerade Familienachsen {0, 2, 4}", any(S == std for S, _ in block))
anchor_ok = True
for S, _ in block:
    d5 = [k for k in range(240) if all(GRAM[k, s] == 0 for s in S)]
    ax5 = axes(d5, d5)
    line = {e for e in ax5 if J2(e) in ax5}
    anchor_ok &= all(sig_vec(e) == e for e in line)
check("in jeder Blocklage ist die komplexe Gerade von D5 punktweise sigma-fest (liegt im Fixraum von sigma)",
      anchor_ok)
jimg = {frozenset(JP[i] for i in S) for S, _ in block}
check("J bildet die Blocklagen aufeinander ab (Orientierung)", jimg == {S for S, _ in block})
only_anchor = [S for S, _ in block
               if all(abs(v[6]) + abs(v[7]) == 0 for v in axes(sorted(S), sorted(S)))]
print(f"  Blocklagen, deren A3 das Ankerpaar nicht beruehrt: {len(only_anchor)}")


def d5_of(S):
    return [k for k in range(240) if all(GRAM[k, s] == 0 for s in S)]


def sigma_orbits_on_axes(ax):
    """Bahnlaengen von sigma auf den Achsen modulo Vorzeichen"""
    left = {min(a, tuple(-t for t in a)) for a in ax}
    out = []
    while left:
        a = left.pop()
        n, b = 1, sig_vec(a)
        while min(b, tuple(-t for t in b)) != a:
            left.discard(min(b, tuple(-t for t in b)))
            n, b = n + 1, sig_vec(b)
        out.append(n)
    return tuple(sorted(out))


orb5 = Counter()
for S, t in sig_inv:
    d5 = d5_of(sorted(S))
    orb5[(t[3:], sigma_orbits_on_axes(axes(d5, d5)), sigma_orbits_on_axes(axes(sorted(S), sorted(S))))] += 1
print("  sigma-invariante Lagen: (Typ, sigma Bahnen auf den 5 D5 Achsen, auf den 3 A3 Achsen): Anzahl")
for k, v in sorted(orb5.items()):
    print(f"    {k}: {v}")
check("in jeder Blocklage zyklt sigma drei D5 Achsen und laesst zwei fest: markierte Teilung 3 + 2",
      all(o5 == (1, 1, 3) for (ty, o5, o3), v in orb5.items() if ty == (3, 1, 0)))
lines = {}
for S, _ in block:
    d5 = d5_of(sorted(S))
    ax5 = axes(d5, d5)
    lines.setdefault(frozenset(e for e in ax5 if J2(e) in ax5), []).append(S)
check("die sechs Blocklagen sind drei sigma-feste komplexe Geraden mal zwei Orientierungen",
      len(lines) == 3 and all(len(v) == 2 for v in lines.values()))
for ln in sorted(lines, key=sorted):
    reps = sorted(e for e in ln if e > tuple(-t for t in e))
    print("    Gerade aufgespannt von", [tuple(t // 2 for t in e) for e in reps], "/ 2")

# ---------------------------------------------------------------- 3. Symmetrie
print("\n3. Symmetrie: Cliffordgruppe in der Paulikarte, Normalisator von sigma")
P = np.array([b.conj() for b in PAULI_OMEGA[1:]] + [PAULI_OMEGA[0].conj() / u0])
PINV = np.linalg.inv(P)
I2 = np.eye(2)
HT = PHI / 2 * np.array([[1, 1], [1, -1]], dtype=complex)
SG = np.diag([1, 1j]).astype(complex)
GENS = [P @ g @ PINV for g in (np.kron(HT, I2), np.kron(I2, HT), np.kron(SG, I2), np.kron(I2, SG),
                               np.diag([1, 1, 1, -1]).astype(complex))]


def key(M):
    q = M * 8
    r = np.rint(q.real).astype(int), np.rint(q.imag).astype(int)
    assert np.allclose(r[0] + 1j * r[1], q)
    return r[0].tobytes() + r[1].tobytes()


seen = {key(np.eye(4, dtype=complex)): np.eye(4, dtype=complex)}
frontier = [np.eye(4, dtype=complex)]
while frontier:
    nxt = []
    for M in frontier:
        for g in GENS:
            N = M @ g
            k = key(N)
            if k not in seen:
                seen[k] = N
                nxt.append(N)
    frontier = nxt
G = list(seen.values())
check("Cliffordgruppe ueber Q(i) in der Paulikarte: Ordnung 46080", len(G) == 46080, f"{len(G)}")
ZR = np.array([[complex(r[2 * k], r[2 * k + 1]) for k in range(4)] for r in RD])


def root_perm(g):
    W = ZR @ g.T
    out = []
    for w in W:
        t = []
        for c in w:
            t += [int(round(c.real)), int(round(c.imag))]
        out.append(IDX[tuple(t)])
    return out


check("alle Erzeuger bilden die 240 TFPT Wurzeln auf sich ab",
      all(sorted(root_perm(g)) == list(range(240)) for g in GENS))
SIGP = P @ SIG2Q @ PINV
SIGP2 = SIGP @ SIGP
cent = [g for g in G if np.allclose(g @ SIGP, SIGP @ g)]
norm = cent + [g for g in G if np.allclose(g @ SIGP, SIGP2 @ g)]
print(f"  Zentralisator von sigma: {len(cent)}, Normalisator von <sigma>: {len(norm)}")


def orbits(group, items):
    items = set(items)
    perms = [root_perm(g) for g in group]
    out = []
    while items:
        seed = items.pop()
        orb = {frozenset(p[i] for i in seed) for p in perms}
        assert all(o == seed or o in items for o in orb)
        items -= orb
        out.append(len(orb))
    return sorted(out)


blockset = [S for S, _ in block]
ob_c, ob_n = orbits(cent, blockset), orbits(norm, blockset)
print(f"  Bahnen auf den sechs Blocklagen: Zentralisator {ob_c}, Normalisator {ob_n}")
check("alle sechs Blocklagen sind unter den Symmetrien von (Gitter, J, sigma) gleichwertig", ob_c == [6])
gen12 = [S for S, t in sig_inv if t[3:] == (0, 0, 0)]
print(f"  Bahnen auf den zwoelf generischen sigma-invarianten Lagen: {orbits(cent, gen12)}")


def tomita(x):
    """a Omega -> a* Omega in TFPT Koordinaten"""
    return (x[0], -x[1], x[2], -x[3], x[4], -x[5], -x[6], x[7])


PAULIS = [np.eye(2, dtype=complex), X2, Y2, Z2]
OM = np.array([1, 0, 0, 1], dtype=complex)
star_ok = True
for r in LROOTS:
    c = [(b.conj() @ r) / 2 for b in PAULI_OMEGA]
    a_star = sum(np.conj(ck) * s for ck, s in zip(c, PAULIS))
    star_ok &= to_tfpt(np.kron(a_star, I2) @ OM, u0) == tomita(to_tfpt(r, u0))
check("die modulare Spiegelung a Omega -> a* Omega ist diag(+,-,+,-,+,-,-,+) auf den TFPT Koordinaten",
      star_ok)
check("sie erhaelt die 240 Wurzeln, antikommutiert mit J und kommutiert mit sigma",
      {tomita(r) for r in RD} == set(RD)
      and all(tomita(J_vec(r)) == tuple(-t for t in J_vec(tomita(r))) for r in RD)
      and all(tomita(sig_vec(r)) == sig_vec(tomita(r)) for r in RD))

# ---------------------------------------------------------------- 4. TFPT Lagen
print("\n4. Die Lagen des Korpus")


def coord_pos(a3c):
    i, j, k = a3c
    return frozenset(IDX[tuple(2 * s * (int(m == p) + t * int(m == q)) for m in range(8))]
                     for p, q in ((i, j), (i, k), (j, k)) for s in (1, -1) for t in (1, -1))


for name, a3c in (("v775 kanonisch, D5 = (0,2,4,6,7), A3 = (1,3,5)", (1, 3, 5)),
                  ("v775 zweite, D5 = (1,3,5,6,7), A3 = (0,2,4)", (0, 2, 4)),
                  ("v634 Arbeitskarte, D5 = (0,1,2,3,4), A3 = (5,6,7)", (5, 6, 7))):
    S = coord_pos(a3c)
    t = classify(S)[3:]
    inv = frozenset(SP[i] for i in S) == S
    print(f"    {name}: Typ {t}, sigma-invariant {inv}, modulare Spiegelung erhaelt sie "
          f"{frozenset(IDX[tomita(RD[i])] for i in S) == S}")
check("die beiden sigma-stabilen Teilungen aus v775 sind Blocklagen, die Arbeitskarte aus v634 nicht",
      classify(coord_pos((1, 3, 5)))[3:] == (3, 1, 0) and classify(coord_pos((0, 2, 4)))[3:] == (3, 1, 0)
      and classify(coord_pos((5, 6, 7)))[3:] == (1, 2, 1))

# ---------------------------------------------------------------- 5. Familien
print("\n5. Zustandsraum 4 und Operatoralgebra 16 unter sigma (Lage D5 = (0,2,4,6,7), A3 = (1,3,5))")
spin = [r for r in RD if r[0] % 2]
c1 = [r for r in spin if sum(1 for k in (1, 3, 5) if r[k] < 0) % 2 == 0]
w4 = {tuple(r[k] for k in (1, 3, 5)) for r in c1}
w16 = {tuple(r[k] for k in (0, 2, 4, 6, 7)) for r in c1}
check("Klasse 1: 64 Wurzeln = 16 Gewichte von D5 mal 4 Gewichte von A3", (len(c1), len(w16), len(w4)) == (64, 16, 4))


def cyc_orbits(ws, cyc):
    ws, out = set(ws), []
    while ws:
        w = ws.pop()
        o = {w, cyc(w), cyc(cyc(w))}
        ws -= o
        out.append(len(o))
    return sorted(out)


o4 = cyc_orbits(w4, lambda w: (w[2], w[0], w[1]))
o16 = cyc_orbits(w16, lambda w: (w[2], w[0], w[1], w[3], w[4]))
check("4 = 1 + 3 und 16 = 4 feste + 4 Dreierbahnen", o4 == [1, 3] and o16 == [1, 1, 1, 1, 3, 3, 3, 3],
      f"{o4}, {o16}")
print("  Die 4 ist der Zustandsraum C^4 selbst: ein sigma-fester Zustand und ein Tripel.")

print("\n6. Das Ankerpaar ist die Zustandsgerade C Omega")
OMEGA = PAULI_OMEGA[0]
img = {z: to_tfpt(2 * z * OMEGA, u0) for z in (1, 1j, 1 + 1j)}
check("Vielfache von Omega liegen ganz im Ankerpaar (Koordinaten 6 und 7)",
      all(v[:6] == (0,) * 6 for v in img.values()), f"{ {str(k): v[6:] for k, v in img.items()} }")
anchor_norms = [sum(t * t for t in img[z]) // 16 for z in (1, 1j, 1 + 1j)]
check("die drei Klassen der Geraden modulo 2 haben Normen (1, 1, 2), Summe 4: der Anker von AX.P1.01",
      anchor_norms == [1, 1, 2] and sum(anchor_norms) == 4, f"{anchor_norms}")
print("  Der Ledger fuehrt c3 = 1/(2 e1(a) pi) mit a = (1, 1, 2). Hier ist a die Zustandsgerade;")
print("  der Faktor 2 pi ist damit nicht hergeleitet.")

print("\nGESAMT:", "alle Pruefungen bestanden" if OK else "FEHLER")

raise SystemExit(0 if OK else 1)
