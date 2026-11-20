#!/usr/bin/env python3
"""A1 -- Das kleinste zusammengesetzte Quantensystem (zwei Qubits) als Ursprung.

Prinzip: Eine nichtkommutative Algebra, die sich selbst beschreibt, braucht ihre
Standardform. Fuer die kleinste, M_2(C), ist das C^2 (x) C^2: System und Spiegel.
Gefragt wird, was die vollstaendige, unverzerrte Beschreibung dieses Systems
erzwingt. Kein TFPT Objekt geht ein; alle Rechnungen sind exakt (F_2, Z[i],
Brueche, dyadische Matrizen).

Ausfuehren: python3 a1_zwei_qubits.py
"""
import itertools
from fractions import Fraction as Fr

import numpy as np
import sympy as sp

OK = True


def check(name, cond, detail=""):
    global OK
    OK &= bool(cond)
    print(f"  [{'ok' if cond else 'FEHLER'}] {name}{(' -- ' + detail) if detail else ''}")


# ---------------------------------------------------------------- 1. Geometrie
print("1. Pauli Geometrie von zwei Qubits (F_2^4 mit symplektischer Form)")
PTS = [v for v in itertools.product((0, 1), repeat=4) if any(v)]  # (x1,z1,x2,z2)


def symp(v, w):
    return (v[0] * w[1] + v[1] * w[0] + v[2] * w[3] + v[3] * w[2]) % 2


def add(v, w):
    return tuple((a + b) % 2 for a, b in zip(v, w))


def weight(v):
    return int(bool(v[0] or v[1])) + int(bool(v[2] or v[3]))


LINES = sorted({frozenset((a, b, add(a, b))) for a, b in itertools.combinations(PTS, 2)
                if symp(a, b) == 0}, key=sorted)
SPREADS = [c for c in itertools.combinations(LINES, 5) if len(frozenset().union(*c)) == 15]
OVOIDS = [c for c in itertools.combinations(PTS, 5)
          if all(symp(a, b) for a, b in itertools.combinations(c, 2))]


def line_type(l):
    return "P" if any(weight(v) == 1 for v in l) else "E"


check("15 Punkte, 15 kommutierende Tripel, 6 Spreads, 6 Ovoide",
      (len(PTS), len(LINES), len(SPREADS), len(OVOIDS)) == (15, 15, 6, 6))
sp_types = {tuple(sorted(line_type(l) for l in s)) for s in SPREADS}
check("jeder vollstaendige Satz von 5 unverzerrten Basen = 3 Produkt + 2 verschraenkt",
      sp_types == {("E", "E", "P", "P", "P")}, f"Typen {sp_types}")
ov_types = {tuple(sorted(weight(v) for v in o)) for o in OVOIDS}
check("jeder Satz von 5 antikommutierenden Paulis = 2 lokal + 3 nichtlokal",
      ov_types == {(1, 1, 2, 2, 2)}, f"Gewichte {ov_types}")
print("  Basislinie (Abzaehlung der Gewicht-1-Punkte, p prim, zwei p-Niveau-Systeme):")
for p in (2, 3, 5, 7):
    print(f"    p={p}: {p * p + 1} Basen = {p + 1} Produkt + {p * p - p} verschraenkt")
print("  Nur p=2 gibt 5 = 3 + 2; ein einzelnes Qubit hat 3 Basen ohne Teilung.")

SWAP4 = np.eye(4)[[0, 2, 1, 3]]
ev = sorted(int(round(x)) for x in np.linalg.eigvalsh(SWAP4))
check("Tausch von System und Spiegel: C^4 = 3 (symmetrisch) + 1 (antisymmetrisch)", ev == [-1, 1, 1, 1])

# lokale Symmetrie eines Ovoids
GL2 = [m for m in itertools.product((0, 1), repeat=4) if (m[0] * m[3] + m[1] * m[2]) % 2]


def act2(m, x, z):
    return ((m[0] * x + m[1] * z) % 2, (m[2] * x + m[3] * z) % 2)


LOCAL = []
for a in GL2:
    for b in GL2:
        for swap in (0, 1):
            def f(v, a=a, b=b, swap=swap):
                u = act2(a, v[0], v[1]) + act2(b, v[2], v[3])
                return (u[2], u[3], u[0], u[1]) if swap else u
            LOCAL.append(f)
DIRAC = ((0, 1, 0, 0), (1, 0, 0, 0), (1, 1, 1, 0), (1, 1, 1, 1), (1, 1, 0, 1))  # ZI, XI, YX, YY, YZ
stab = [f for f in LOCAL if {f(v) for v in DIRAC} == set(DIRAC)]
loc, non = DIRAC[:2], DIRAC[2:]
perm_non = {tuple(non.index(f(v)) for v in non) for f in stab}
perm_loc = {tuple(loc.index(f(v)) for v in loc) for f in stab}
check("lokale Symmetrie des Ovoids {ZI,XI,YX,YY,YZ}: Bahnen 2 + 3, S_2 x S_3",
      (len(stab), len(perm_loc), len(perm_non)) == (12, 2, 6),
      f"|Stab|={len(stab)}, wirkt als S_{len(perm_loc)} auf dem Paar, Ordnung {len(perm_non)} auf dem Tripel")

# ---------------------------------------------------------------- 2. Operatoralgebra
print("\n2. Operatoralgebra M_4(C) = Lambda^gerade C^5 und eine Standardmodellgeneration")
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
Y = 1j * X @ Z
P1 = {(0, 0): I2, (1, 0): X, (0, 1): Z, (1, 1): Y}


def pauli(v):
    return np.kron(P1[(v[0], v[1])], P1[(v[2], v[3])])


# Plaetze 0,1,2 = nichtlokale Gammas, 3,4 = lokale Gammas
PLACES = non + loc
G = [pauli(v) for v in PLACES]
check("fuenf Gammas antikommutieren und quadrieren zu 1",
      all(np.array_equal(G[a] @ G[b] + G[b] @ G[a], 2 * (a == b) * np.eye(4))
          for a in range(5) for b in range(5)))
prod5 = np.linalg.multi_dot(G)
check("Produkt aller fuenf ist ein Skalar (Lambda^4 = Lambda^1 dual)",
      np.array_equal(prod5, prod5[0, 0] * np.eye(4)), f"Skalar {prod5[0, 0]}")
EVEN = [S for r in (0, 2, 4) for S in itertools.combinations(range(5), r)]
ops = []
for S in EVEN:
    M = np.eye(4, dtype=complex)
    for a in S:
        M = M @ G[a]
    ops.append(M.reshape(-1))
check("1 + 10 + 5 = 16 gerade Monome sind eine Basis von M_4(C)",
      np.linalg.matrix_rank(np.array(ops)) == 16)

YPL = [Fr(-1, 3)] * 3 + [Fr(1, 2)] * 2  # einzige spurfreie Klassenfunktion der 3+2 Teilung, Skala s.u.
NAMES = {(0, 0): "nu^c (1,1)", (2, 0): "u^c (3b,1)", (1, 1): "Q (3,2)", (0, 2): "e^c (1,1)",
         (3, 1): "L (1,2)", (2, 2): "d^c (3b,1)"}
table = {}
for S in EVEN:
    n3, n2 = sum(a < 3 for a in S), sum(a >= 3 for a in S)
    v = (0, 0, 0, 0)
    for a in S:
        v = add(v, PLACES[a])
    yv = sum((YPL[a] for a in S), Fr(0))
    table.setdefault((n3, n2), []).append((yv, weight(v)))
print("  Multiplett      Anzahl  Y      Pauli Gewicht (0 = Eins, 1 = lokal, 2 = verschraenkend)")
allY = []
for key in [(0, 0), (2, 0), (1, 1), (0, 2), (3, 1), (2, 2)]:
    rows = table[key]
    ys, ws = {r[0] for r in rows}, sorted({r[1] for r in rows})
    allY += [r[0] for r in rows]
    print(f"  {NAMES[key]:<15} {len(rows):>4}   {str(ys.pop()):>5}  {ws}")
idx3 = {0: Fr(0), 1: Fr(1, 2), 2: Fr(1, 2), 3: Fr(0)}
idx2 = {0: Fr(0), 1: Fr(1, 2), 2: Fr(0)}
a_grav = sum(allY)
a_y3 = sum(y ** 3 for y in allY)
a_c = sum(rows[0][0] * idx3[k[0]] * sp.binomial(2, k[1]) for k, rows in table.items())
a_w = sum(rows[0][0] * idx2[k[1]] * sp.binomial(3, k[0]) for k, rows in table.items())
check("Anomalien verschwinden: Sum Y, Sum Y^3, SU(3)^2 Y, SU(2)^2 Y",
      (a_grav, a_y3, a_c, a_w) == (0, 0, 0, 0), f"{a_grav}, {a_y3}, {a_c}, {a_w}")
check("6Y ganzzahlig mit ggT 1 (Z_6 Quotient)",
      sp.igcd(*[int(6 * y) for y in set(allY) if y]) == 1)
ymin = min(abs(y) for y in allY if y)
print("  Y ist bis auf die Skala eindeutig: spurfrei und konstant auf den 3 und den 2 Plaetzen.")
print(f"  Skala: kleinste Ladung {ymin}; e^c (das Produkt der beiden lokalen Gammas) hat Y = 1.")

# ---------------------------------------------------------------- 3. Dirac Algebra
print("\n3. Dirac Algebra: dieselbe 3 + 2 Teilung als 3 Raum + (Zeit, Chiralitaet)")
G0 = pauli(loc[0])
GK = [1j * pauli(v) for v in non]
GAM = [G0] + GK
eta = np.diag([1, -1, -1, -1])
check("Cl(1,3): {G^mu, G^nu} = 2 eta^{mu nu}",
      all(np.array_equal(GAM[m] @ GAM[n] + GAM[n] @ GAM[m], 2 * eta[m, n] * np.eye(4))
          for m in range(4) for n in range(4)))
G5 = 1j * GAM[0] @ GAM[1] @ GAM[2] @ GAM[3]
check("G^5 = i G^0 G^1 G^2 G^3 ist (bis aufs Vorzeichen) das zweite lokale Gamma",
      np.array_equal(G5, pauli(loc[1])) or np.array_equal(G5, -pauli(loc[1])))
spin_w = {weight(add(a, b)) for a, b in itertools.combinations(non, 2)}
boost_w = {weight(add(loc[0], b)) for b in non}
check("Drehgeneratoren lokal (Gewicht 1, Spin Qubit), Boosts verschraenkend (Gewicht 2)",
      (spin_w, boost_w) == ({1}, {2}))
print("  Zeit unter den drei nichtlokalen Gammas wuerde die S_3 Symmetrie des Tripels brechen;")
print("  mit Zeit im lokalen Paar bleiben drei gleichberechtigte Raumrichtungen: 3 + 1.")

# ---------------------------------------------------------------- 4. Clifford Gruppe und E8
print("\n4. Clifford Gruppe ueber Q(i) und ihr invariantes Gitter")
PHI = 1 + 1j
HT = PHI / 2 * np.array([[1, 1], [1, -1]], dtype=complex)  # zeta_8 mal Hadamard, Eintraege in Q(i)
SG = np.diag([1, 1j]).astype(complex)
GENS = [np.kron(HT, I2), np.kron(I2, HT), np.kron(SG, I2), np.kron(I2, SG),
        np.diag([1, 1, 1, -1]).astype(complex)]


def key(M):
    q = M * 4
    r = np.rint(q.real).astype(int), np.rint(q.imag).astype(int)
    assert np.array_equal(r[0] + 1j * r[1], q)
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
check("Gruppenordnung 46080 = 4 * 16 * 720 (Zentrum mu_4, Paulis, Sp(4,2) = S_6)",
      len(seen) == 46080, f"|G| = {len(seen)}")


def gint(z):
    return complex(round(z.real), round(z.imag))


def n2(c):
    return round(c.real) ** 2 + round(c.imag) ** 2


def in_L(x):
    if any(gint(c) != c for c in x):
        return False
    par = {(int(c.real) + int(c.imag)) % 2 for c in x}
    s = sum(x)
    return len(par) == 1 and int(s.real) % 2 == 0 and int(s.imag) % 2 == 0


BASIS = [np.array(b, dtype=complex) for b in
         [(1, 1, 1, 1), (0, PHI, 0, PHI), (0, 0, PHI, PHI), (0, 0, 0, 2)]]
check("alle Erzeuger bilden das Gitter L = {x in Z[i]^4: x_j kongruent mod 1+i, Summe in 2Z[i]} in sich ab",
      all(in_L(g @ b) for g in GENS for b in BASIS))
SMALL = [complex(a, b) for a in range(-2, 3) for b in range(-2, 3) if a * a + b * b <= 4]
ROOTS = [np.array(x) for x in itertools.product(SMALL, repeat=4)
         if sum(n2(c) for c in x) == 4 and in_L(x)]
shapes = {}
for r in ROOTS:
    sh = tuple(sorted(n2(c) for c in r))
    shapes[sh] = shapes.get(sh, 0) + 1
check("240 kuerzeste Vektoren = 16 + 96 + 128 (E8 Wurzeln als gausssches Gitter)",
      len(ROOTS) == 240 and sorted(shapes.values()) == [16, 96, 128], f"{shapes}")


def cls(x):
    c1 = x[0]
    c2 = (x[1] - c1) / PHI
    c3 = (x[2] - c1) / PHI
    c4 = (x[3] - c1 - PHI * c2 - PHI * c3) / 2
    return tuple((round(c.real) + round(c.imag)) % 2 for c in (c1, c2, c3, c4))


def stab_line(r):
    return frozenset(v for v in PTS if any(np.array_equal(pauli(v) @ r, s * r) for s in (1, -1)))


by_cls = {}
for r in ROOTS:
    by_cls.setdefault(cls(r), []).append(r)
cls_line = {c: {stab_line(r) for r in rs} for c, rs in by_cls.items()}
check("L/(1+i)L = F_2^4: 15 Klassen zu je 16 Wurzeln",
      len(by_cls) == 15 and all(len(v) == 16 for v in by_cls.values()))
check("jede Wurzel ist ein Stabilisatorzustand; Klasse <-> kommutierendes Tripel (Basis), bijektiv",
      all(len(s) == 1 and next(iter(s)) in LINES for s in cls_line.values())
      and len({next(iter(s)) for s in cls_line.values()}) == 15)
packets = [[r for r in ROOTS if PLACES[a] in stab_line(r)] for a in range(5)]
check("fuenf Gammas zerlegen die 240 Wurzeln in 5 x 48 = (3 + 2) x 48 = 144 + 96",
      [len(p) for p in packets] == [48] * 5 and sum(len(p) for p in packets) == 240)
hopf = {}
for r in ROOTS:
    h = tuple(int(round((r.conj() @ G[a] @ r).real)) for a in range(5))
    hopf.setdefault(h, []).append(r)
verts = sorted(hopf)
check("Hopf Abbildung <r|gamma_a|r>: Bild = 10 Ecken +-4 e_a (Vektorgewichte von D5), Faser je 24",
      len(verts) == 10 and all(sorted(map(abs, v)) == [0, 0, 0, 0, 4] for v in verts)
      and all(len(f) == 24 for f in hopf.values()))
fib = hopf[verts[0]]
nb = {sum(1 for s in fib if round((r.conj() @ s).real) == 2) for r in fib}
check("jede Faser ist ein 24-Zell (D4 Wurzelsystem: 8 Nachbarn bei Skalarprodukt 2)", nb == {8})
print("  3 + 2 auf den Gammas gibt auf den Ecken 6 + 4: SO(10) > SO(6) x SO(4), also A3 = D3 und A1 + A1.")

# ---------------------------------------------------------------- 5. Standardform, Zustand, Zeit
print("\n5. Standardform von M_2(C): Zustand, Spiegelung und Zeit auf C^2 (x) C^2")
pa, pb = sp.symbols("p_0 p_1", positive=True)  # Gewichte des treuen Zustands (Normierung unerheblich)
rho = sp.diag(pa, pb)
Om = sp.Matrix([sp.sqrt(pa), 0, 0, sp.sqrt(pb)])
D12 = sp.kronecker_product(rho ** sp.Rational(1, 2), rho ** sp.Rational(-1, 2))
SW = sp.Matrix(4, 4, lambda i, j: 1 if (i // 2 == j % 2 and i % 2 == j // 2) else 0)
good = True
for i in range(2):
    for j in range(2):
        E = sp.zeros(2, 2)
        E[i, j] = 1
        lhs = SW * (D12 * sp.kronecker_product(E, sp.eye(2)) * Om).conjugate()
        rhs = sp.kronecker_product(E.T, sp.eye(2)) * Om
        good &= sp.simplify(lhs - rhs) == sp.zeros(4, 1)
check("J Delta^(1/2) a Omega = a* Omega fuer alle Matrixeinheiten (J = Tausch mit Konjugation)", good)
lam = [sp.expand_log(sp.log(rho[i, i] / rho[j, j]), force=True) for i in range(2) for j in range(2)]
check("modulares Spektrum {0, 0, +eps, -eps}: Ladungen 0, 0, +1, -1 (Einheitswindung je Periode)",
      sp.simplify(lam[1] + lam[2]) == 0 and lam[0] == 0 and lam[3] == 0)
print("  Der Spiegel (Kommutant) ist J M J; die Zeit ist der modulare Fluss des Zustands.")
print("  Zweite Stufe: Standardform von M_4(C) ist C^16; Links- und Rechtsmultiplikation mit den")
print("  fuenf Gammas geben zehn antikommutierende Operatoren, also Cl(10) auf Lambda C^5.")

print("\nGESAMT:", "alle Pruefungen bestanden" if OK else "FEHLER")
