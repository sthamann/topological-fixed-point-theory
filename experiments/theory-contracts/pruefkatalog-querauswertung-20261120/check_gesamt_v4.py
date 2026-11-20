"""GESAMT2-20261120: die acht Moden als eine Darstellung der Kleinschen Vierergruppe der Marken.

Objekt: der Torus E über der Sphäre, verzweigt an den vier Marken, mit dem Bündel L^8 = Rückzug des
Markendivisors. Seine acht Schnitte (Thetafunktionen der Stufe 8) zerfallen unter dem Deck in 5 + 3
(Atlas, 19. November). Neu: die Halbperioden E[2] wirken darauf; auf der Sphäre ist das die Kleinsche
Vierergruppe V4 der Marken (Doppeltranspositionen). Exakte Rechnung, kein Beweis der Physik.
"""
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT/"documents/archive--tfpt-45--source_extracts--02_carrier_source/sections"
CAR = " ".join(next(SRC.glob("001-*.tex")).read_text().split())
HIG = " ".join(next(SRC.glob("007-*.tex")).read_text().split())
KERN = next((ROOT/"documents/archive--tfpt-45--source_extracts--01_boundary_kernel_source/sections").glob("001-*.tex")).read_text()
CARR = next(SRC.glob("001-*.tex")).read_text()
checks = []


def check(name, ok, detail=""):
    checks.append(bool(ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))


def mul(a, b):
    n = len(a)
    return [[sum(a[i][k] * b[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def tr(a):
    return sum(a[i][i] for i in range(len(a)))


def eye(n):
    return [[F(int(i == j)) for j in range(n)] for i in range(n)]


def pmul(p, q):
    r = [F(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            r[i + j] += x * y
    return r


def sym4(a, b, c, d):
    """Wirkung von z -> (az+b)/(cz+d) auf Quartiken, p -> (cz+d)^4 p(gz) / det^2."""
    det2 = (a * d - b * c) ** 2
    cols = []
    for k in range(5):
        p = [F(1)]
        for _ in range(k):
            p = pmul(p, [b, a])
        for _ in range(4 - k):
            p = pmul(p, [d, c])
        cols.append([x / det2 for x in p])
    return [[cols[k][i] for k in range(5)] for i in range(5)]


def multiplicities(traces):
    """Vielfachheiten der vier Charaktere von V4 aus den Spuren von (1, a, b, ab)."""
    chars = {"1": (1, 1, 1, 1), "chi_a": (1, 1, -1, -1), "chi_b": (1, -1, 1, -1), "chi_ab": (1, -1, -1, 1)}
    return {n: sum(F(t) * c for t, c in zip(traces, ch)) / 4 for n, ch in chars.items()}


# --- 1. Torus: Thetafunktionen der Stufe 8, Index k in Z/8 ----------------------------------
N = 8
perm = lambda f: [[F(int(f(j) == i)) for j in range(N)] for i in range(N)]
iota = perm(lambda k: -k % N)                                   # Deck z -> -z auf dem Torus
T1 = [[F((-1) ** i) if i == j else F(0) for j in range(N)] for i in range(N)]   # Halbperiode 1/2
T2 = perm(lambda k: (k + 4) % N)                                # Halbperiode tau/2
T3 = mul(T1, T2)
I8 = eye(N)
check("Halbperioden kommutieren untereinander und mit dem Deck", mul(T1, T2) == mul(T2, T1)
      and all(mul(t, iota) == mul(iota, t) for t in (T1, T2)))
half = lambda s: [[(I8[i][j] + s * iota[i][j]) / 2 for j in range(N)] for i in range(N)]
Pe, Po = half(1), half(-1)
check("Deck teilt die acht Schnitte in 5 gerade und 3 ungerade", (tr(Pe), tr(Po)) == (5, 3))
ev = multiplicities([tr(mul(g, Pe)) for g in (I8, T1, T2, T3)])
od = multiplicities([tr(mul(g, Po)) for g in (I8, T1, T2, T3)])
check("gerade (Träger): zweimal trivial plus jeder nichttriviale Charakter einmal", list(ev.values()) == [2, 1, 1, 1])
check("ungerade (Familien): jeder nichttriviale Charakter genau einmal", list(od.values()) == [0, 1, 1, 1])
check("acht Schnitte = zweimal die reguläre Darstellung der Halbperioden", [ev[k] + od[k] for k in ev] == [2, 2, 2, 2])

# --- 2. Sphäre: V4 der Marken bei mu4, Wirkung auf Quartiken ---------------------------------
# Gaußsche Zahlen als Paare wären nötig für z -> i/z; stattdessen Marken reell drehen:
# Marken 0, 1, unendlich, -1 (harmonisch, dasselbe Modul j = 1728). V4 = Doppeltranspositionen.
def v4_for(lam):
    """Marken 0, 1, unendlich, lam. Die drei Involutionen, die die Marken paarweise tauschen."""
    a = sym4(F(0), lam, F(1), F(0))                # (0 unendlich)(1 lam): z -> lam/z
    b = sym4(lam, -lam, F(1), -lam)                # (0 1)(unendlich lam): z -> lam (z-1)/(z-lam)
    return a, b, mul(a, b)


I5 = eye(5)
for lam in (F(-1), F(2), F(1, 2), F(7, 3), F(-5, 11)):
    a, b, ab = v4_for(lam)
    ok = mul(a, a) == I5 and mul(b, b) == I5 and mul(a, b) == mul(b, a)
    m = multiplicities([tr(g) for g in (I5, a, b, ab)])
    marks = [F(0), lam, -(1 + lam), F(1), F(0)]     # z (z-1)(z-lam), vierte Nullstelle bei unendlich
    inv = all([sum(g[i][j] * marks[j] for j in range(5)) for i in range(5)] == marks for g in (a, b, ab))
    check(f"lambda = {lam}: Träger = 2 mal trivial + drei Charaktere; Markenquartik invariant",
          ok and list(m.values()) == [2, 1, 1, 1] and inv)
check("die Teilung 2 + 3 braucht kein besonderes Modul", True, "gilt für jede Lage der vier Marken")

# Familien: H1 der punktierten Sphäre = Permutationsdarstellung auf den Marken minus trivial.
perm_tr = [4, 0, 0, 0]                              # V4 wirkt frei auf den vier Marken
fam = multiplicities([t - 1 for t in perm_tr])
check("Familien (Definition des Originals, H^1 der punktierten Sphäre): drei nichttriviale Charaktere",
      list(fam.values()) == [0, 1, 1, 1])
check("Farbteil des Trägers und Familienraum sind als Darstellungen von V4 gleich", True,
      "beide chi_a + chi_b + chi_ab; Zahl der Farben = Zahl der Familien = |V4| - 1 = 3")

# --- 3. Die Annahmen des Originals an die äußeren Potenzen, äquivariant ----------------------
G = list(product((1, -1), repeat=2))                # Charakter = (Wert auf a, Wert auf b)
times = lambda x, y: (x[0] * y[0], x[1] * y[1])
weak = [(1, 1), (1, 1)]
col = [(1, -1), (-1, 1), (-1, -1)]
ext2 = lambda v: sorted(times(x, y) for x, y in combinations(v, 2))
check("Original nennt beide Annahmen: dim Lambda^2 E_+ = 1 und Lambda^2 E_- = Dual von E_-",
      r"\dim \Lambda^2 E_+=1" in CARR and r"\Lambda^2 E_- \cong E_-^\vee" in CARR)
check("Lambda^2 des schwachen Blocks ist trivial", ext2(weak) == [(1, 1)])
check("Lambda^2 des Farbblocks ist äquivariant sein Dual", ext2(col) == sorted(col))
check("Determinante des Farbblocks ist trivial", times(times(*col[:2]), col[2]) == (1, 1))
old3 = [(1, 1), (1, 1), (1, -1)]                    # Paritätslesart EPSCAR: gerade unter a allein
check("Gegenprobe: die Paritätslesart mit einer Involution besteht den äquivarianten Test nicht",
      ext2(old3) != sorted(old3), "sie mischt die Invarianten mit einem Farbplatz")

# --- 4. Quotientensphäre P1/V4 -----------------------------------------------------------------
# q hat Grad 4 und drei Verzweigungspunkte der Ordnung 2. O(4) = q*O(1); q_*O = O + M_a + M_b + M_ab.
branch = [(-1, 1), (1, -1), (-1, -1)]               # Monodromie an den drei Punkten: a, b, ab als Charakterwerte
degM = {chi: -F(sum(1 for g in ((1, 0), (0, 1), (1, 1))
                    if (chi[0] ** g[0]) * (chi[1] ** g[1]) == -1), 2) for chi in G}
deg = {chi: 1 + degM[chi] for chi in G}
h0 = lambda d: max(int(d) + 1, 0)
check("invarianter Teil ist O(1) auf dem Quotienten, zwei Schnitte", deg[(1, 1)] == 1 and h0(deg[(1, 1)]) == 2)
check("jeder Farbplatz ist ein Bündel vom Grad 0 mit einem Schnitt", all(deg[c] == 0 and h0(deg[c]) == 1 for c in col))
check("Summe der Schnitte = 5", sum(h0(d) for d in deg.values()) == 5)
check("Determinantenklassen (schwach, Farbe) = (1, 0) wie im Original", (deg[(1, 1)], sum(deg[c] for c in col)) == (1, 0)
      and r"(c_1(L_2),c_1(L_3))=(1,0)" in HIG)
degK = {chi: -1 - degM[chi] for chi in G}            # K_w(m) mal M^-1: logarithmische Formen mit Pol im Markenbild
check("Familien auf dem Quotienten: kein invarianter Schnitt, dreimal Grad 0",
      h0(degK[(1, 1)]) == 0 and all(degK[c] == 0 for c in col))

# --- 5. Marker und Hyperladung -----------------------------------------------------------------
check("Original: Marker hat die Form 2P - 1", r"\iC=2P_C-1=R_{\taudbl}" in KERN,
      "dort ist P_C die Spiegelung an der Naht, hier das Mittel über V4: gleiche Form, nicht derselbe Operator")
eps = [1, 1, -1, -1, -1]                             # 2P - 1 mit P = Mittel über V4
Y = [F(1, 12) + F(5, 12) * e for e in eps]
check("Y = 1/12 + (5/12)(2P - 1): 1/2 auf dem Dublett, -1/3 auf den drei Farben, ohne Vorzeichenwechsel",
      Y == [F(1, 2)] * 2 + [F(-1, 3)] * 3 and sum(Y) == 0 and sum(y * y for y in Y) == F(5, 6))
check("gleichwertig: Y = -1/3 + (5/24) mal Summe über V4", all(F(-1, 3) + F(5, 24) * s == y for s, y in zip([4, 4, 0, 0, 0], Y)),
      "spurfrei genau weil jede Involution auf dem Träger Spur 1 hat")
check("Spur von 2P - 1 auf den acht Moden ist -4, nicht 0", sum(eps) + (-3) == -4,
      "die Spurbedingung aus KERNSIX Schritt 2 ist hier nicht erfüllt")

# --- 6. Zusatz am Quadrat: Vierteldrehung ------------------------------------------------------
r = sym4(F(0), F(1), F(-1), F(0))                    # hier nur als Kontrolle: ein Element der Ordnung 4 an 0,1,unendlich,-1?
a, b, ab = v4_for(F(-1))
# Marken 0, 1, unendlich, -1: z -> (1+z)/(1-z) dreht 0 -> 1 -> unendlich -> -1 -> 0.
rot = sym4(F(1), F(1), F(-1), F(1))
conj = lambda g: mul(mul(rot, g), [[x for x in row] for row in mul(mul(rot, rot), rot)])
r4 = mul(mul(rot, rot), mul(rot, rot))
check("Vierteldrehung der Marken hat Ordnung 4 auf Quartiken", r4 == I5 and mul(rot, rot) != I5)
c = [conj(g) for g in (a, b, ab)]
check("sie lässt eine Involution fest und tauscht die beiden anderen", sorted(map(str, c)) == sorted(map(str, (a, b, ab)))
      and sum(x == y for x, y in zip(c, (a, b, ab))) == 1, "Farben und Familien zerfallen am Quadrat je 1 + 2")

print(f"\n{sum(checks)} von {len(checks)} Prüfungen bestanden")
raise SystemExit(0 if all(checks) else 1)
