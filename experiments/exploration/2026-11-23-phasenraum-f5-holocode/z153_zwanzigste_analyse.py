"""Z153: Zwanzigste externe Analyse ('Das Austauschgatter steckt in E8; dieselbe Struktur erzeugt positive Krümmung').

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart, N numerisch, O offen.

Die Analyse baut auf z145 auf (das z120-Gatter als relative Drehung der E8-Faktoren D5 und A3). Ihre Programme liegen
nicht vor; alle Aussagen sind hier neu gerechnet, die Konstruktion des 'Codes' ist rekonstruiert. Behauptungen:
  (1) Die gegenüber allen Standardmodell-Generatoren neutralen E8-Wurzeln sind 20 = 12 (A3) + 8 Klebewurzeln und bilden A4.
  (2) Unter einer konkreten Hebung der Familienpermutationen gibt es einen vierdimensionalen Raum, auf dem die relative
      Drehung als diag(1,1,1,-i) = (I + iP)/(1 + i) wirkt, also als das Austauschgatter mit U^2 = P.
  (3) Die Familie dieser Unterräume ist RP^4 (Überlagerung S^4), projizierte Krümmung R_mn = 3/r^2 g_mn, also Lambda > 0.
  (4) Die 240 Wurzeln, auf den A3-Raum projiziert, geben 27 Punkte (40, 128, 60, 12) und eine Fourierreihe, die bis zur
      sechsten Ordnung drehinvariant ist: 1 - k^2/8 + k^4/160 - k^6/5760.
  (5) Härtetest: Im Zweinahtmodell hat der niveauabgeglichene Sektor die Dimensionen 1, 2, 1 je Gesamtladung; ein
      vierdimensionaler Raum fester Ladung fehlt. Die relative Drehung ist keine Automorphie von E8.

Teil A  Zentralisator: 20 Wurzeln, A4, die acht Klebewurzeln und ihre B-L-Ladung.
Teil B  Die relative Drehung auf den E8-Strömen: Gruppenelement (Involution) gegen Twist (1 und -i).
Teil C  Der vierdimensionale Code im A4-Sektor (als sl(5)): Existenz, Hebung, Wahlfreiheit, welches Zweiqubit-Gatter.
Teil D  Geometrie: projizierte Verbindung auf den Hyperebenen von R^5 = Tangentialbündel von S^4; Vorzeichen = reelle Form.
Teil E  Projektion auf A3: 27 Punkte, Fourierreihe; die Isotropie gilt für jeden Unterraum (E8-Invarianten 2, 8, ...).
Teil F  Bilanz und Urteil.
"""
import itertools
import json
import os
import time
from fractions import Fraction as Fr
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
H = Fr(1, 2)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val, text):
    BEFUNDE[name] = val
    say(f"   [Befund: {val}] {text}")


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


# ---------------------------------------------------------------------------------------------------------------------
say("Z153: Zwanzigste externe Analyse (E8-Kleber als Prozess, positive Krümmung, Kontinuum)")
say("=" * 100)
say("\nTeil A  Zentralisator der Standardmodellgruppe in E8")
roots = []
for i, j in itertools.combinations(range(8), 2):
    for si, sj in itertools.product((1, -1), repeat=2):
        v = [Fr(0)] * 8
        v[i], v[j] = Fr(si), Fr(sj)
        roots.append(tuple(v))
for s in itertools.product((1, -1), repeat=8):
    if s.count(-1) % 2 == 0:
        roots.append(tuple(H * x for x in s))
roots = sorted(roots)
check("A0 E8 = D8 mit Spinorklasse: 240 Wurzeln der Norm 2", len(roots) == 240 and all(dot(r, r) == 2 for r in roots))


def typ(r):
    a, b = dot(r[:5], r[:5]), dot(r[5:], r[5:])
    return {(2, 0): "D5", (0, 2): "A3", (1, 1): "(10,6)"}.get((a, b), "Kleber")


zahl = {t: sum(1 for r in roots if typ(r) == t) for t in ("D5", "A3", "(10,6)", "Kleber")}
check("A0 Zerlegung unter D5 x A3: 40 + 12 + 60 + 128", zahl == {"D5": 40, "A3": 12, "(10,6)": 60, "Kleber": 128})
# Standardmodell in SU(5) in SO(10): Farbe auf 1..3, schwach auf 4..5, Y = (-1/3,-1/3,-1/3,1/2,1/2)
sm_cartan = [(1, -1, 0, 0, 0, 0, 0, 0), (0, 1, -1, 0, 0, 0, 0, 0), (0, 0, 0, 1, -1, 0, 0, 0),
             (Fr(-1, 3), Fr(-1, 3), Fr(-1, 3), H, H, 0, 0, 0)]
sm_roots = [r for r in roots if typ(r) == "D5" and sum(r[:5]) == 0
            and ((r[3] == 0 and r[4] == 0) or (r[0] == 0 and r[1] == 0 and r[2] == 0))]
check("A1 Standardmodell: 8 Wurzeln (6 Farbe + 2 schwach) und Hyperladung, Rang 4", len(sm_roots) == 8)
neutral = [r for r in roots if all(dot(r, h) == 0 for h in sm_cartan)]
vert = all(dot(r, b) == 0 for r in neutral for b in sm_roots)
n_a3 = [r for r in neutral if typ(r) == "A3"]
n_kl = [r for r in neutral if typ(r) == "Kleber"]
say(f"   neutrale Wurzeln: {len(neutral)} = {len(n_a3)} (A3) + {len(n_kl)} (Kleber); andere Typen: "
    f"{len(neutral) - len(n_a3) - len(n_kl)}")
check("A2 Genau 20 Wurzeln sind gegenüber allen Standardmodell-Generatoren neutral: 12 (A3) + 8 Klebewurzeln",
      len(neutral) == 20 and len(n_a3) == 12 and len(n_kl) == 8 and vert)
f = (Fr(1, 97), Fr(2, 89), Fr(3, 83), Fr(5, 79), Fr(7, 73), Fr(11, 3), Fr(13, 7), Fr(17, 29))
pos = [r for r in neutral if dot(r, f) > 0]
einf = sorted(r for r in pos if not any(tuple(x - y for x, y in zip(r, q)) in pos for q in pos))
cart = sorted(sorted(int(2 * dot(a, b) / dot(b, b)) for b in einf) for a in einf)
check("A3 Die 20 Wurzeln bilden A4 (vier einfache Wurzeln, Cartan-Matrix einer Kette)",
      len(einf) == 4 and cart == [[-1, -1, 0, 2], [-1, -1, 0, 2], [-1, 0, 0, 2], [-1, 0, 0, 2]])
check("A4 Die acht Klebewurzeln haben gleiche D5-Koordinaten +-(1/2)^5: der SU(5)-Singulett der 16 (nu_R-Richtung) mal 4",
      all(len(set(r[:5])) == 1 for r in n_kl))
xl = sorted({sum(r[:5]) for r in n_kl})
ZAHLEN["X-Ladung (Summe der D5-Koordinaten) der acht Klebewurzeln"] = [str(x) for x in xl]
check("A5 Sie sind standardmodellneutral, tragen aber die Ladung X = +-5/2 der U(1) neben SU(5) (B - L): nicht SO(10)-neutral",
      xl == [Fr(-5, 2), Fr(5, 2)])
befund("A", "E; klassisch",
       "Behauptung (1) stimmt. Es ist die bekannte Zerlegung E8 > SU(5) x SU(5): Der Zentralisator der "
       "Standardmodell-Cartanalgebra ist die zweite SU(5), sie enthält die Familien-A3. Die acht Zusatzwurzeln sind die "
       "Richtung des rechtshändigen Neutrinos mal der 4 von A3; sie tragen B - L.")

# ---------------------------------------------------------------------------------------------------------------------
say("\nTeil B  Die relative Drehung auf den E8-Strömen")
lam = (H, H, H, H, H, 0, 0, 0)                      # Zentrum von Spin(10): i auf der 16
mu = (0, 0, 0, 0, 0, H, H, H)                       # Zentrum von SU(4)


def ph4(x):                                           # e^(2 pi i x) für x in Z/4 als Potenz von i
    q = (4 * x) % 4
    assert q.denominator == 1
    return int(q)


zus = {ph4(dot(r, lam) + dot(r, mu)) for r in roots}
check("B1 Verklebung: Zentrum(D5) mal Zentrum(A3) wirkt auf allen 240 Wurzeln trivial", zus == {0})
rel = {t: sorted({ph4(dot(r, lam) - dot(r, mu)) for r in roots if typ(r) == t}) for t in zahl}
ZAHLEN["relative Drehung als Gruppenelement: Potenz von i je Typ"] = rel
check("B2 Als Gruppenelement ist die relative Drehung eine Involution: +1 auf 112 Wurzeln, -1 auf den 128 Klebewurzeln",
      rel == {"D5": [0], "A3": [0], "(10,6)": [0], "Kleber": [2]})


def twist(r):                                         # e^(-2 pi i L0^D5) e^(+2 pi i L0^A3), L0 = Norm/2
    return ph4(-(dot(r[:5], r[:5]) - dot(r[5:], r[5:])) / 2)


tw = {t: sorted({twist(r) for r in roots if typ(r) == t}) for t in zahl}
ZAHLEN["relative Drehung als Twist: Potenz von i je Typ"] = tw
check("B3 Als Twist (z145) wirkt sie mit 1 auf 112 Wurzeln und mit -i auf den 128 Klebewurzeln",
      tw == {"D5": [0], "A3": [0], "(10,6)": [0], "Kleber": [3]})
rs = set(roots)
paar = next((a, b) for a in n_kl for b in n_kl if tuple(x + y for x, y in zip(a, b)) in rs)
summe = tuple(x + y for x, y in zip(*paar))
check("B4 Der Twist ist keine Automorphie: zwei Klebewurzeln mit Summe in A3 geben (-i)(-i) = -1 statt 1",
      (twist(paar[0]) + twist(paar[1])) % 4 == 2 and twist(summe) == 0 and typ(summe) == "A3")
check("B5 Twist^2 = Gruppenelement aus B2: Der Twist ist eine nicht multiplikative Wurzel einer inneren Automorphie",
      all((2 * twist(r)) % 4 == ph4(dot(r, lam) - dot(r, mu)) for r in roots))
befund("B", "E",
       "Der Eigenwert -i des Gatters ist die Twistphase der Klebeklasse (schon z145). Als Gruppenelement von E8 hat die "
       "relative Drehung nur die Eigenwerte +-1. Die Analyse sagt selbst, dass der Twist keine Automorphie ist.")

# ---------------------------------------------------------------------------------------------------------------------
say("\nTeil C  Der vierdimensionale Code im A4-Sektor, gerechnet in sl(5): A3 = sl(4) auf e1..e4, Kleber = E_i5, E_5i")
I5 = sp.eye(5)


def E(i, j):
    m = sp.zeros(5, 5)
    m[i - 1, j - 1] = 1
    return m


def T(m):                                             # Twist: 1 auf gl(4) + Cartan, -i auf den Klebeströmen
    out = sp.zeros(5, 5)
    for a in range(5):
        for b in range(5):
            kleber = (a == 4) != (b == 4)
            out[a, b] = m[a, b] * (-sp.I if kleber else 1)
    return out


def perm(sig):                                        # Permutationsmatrix auf e1..e4, e5 fest
    m = sp.zeros(5, 5)
    for a, b in enumerate(sig):
        m[b, a] = 1
    m[4, 4] = 1
    return m


Hs = [E(1, 1) - E(2, 2), E(2, 2) - E(3, 3), E(3, 3) - E(4, 4)]
s = E(1, 5) + E(2, 5) + E(3, 5) + E(4, 5)
ok = all(T(h) == h for h in Hs) and T(s) == -sp.I * s
check("C1 Code 1 = Cartan(A3) + Summe der vier Klebeströme: Twist = diag(1, 1, 1, -i)", ok)
S4 = [perm(p) for p in itertools.permutations(range(4))]
inv = all(P * s * P.T == s for P in S4)
chars = sorted({sum(1 for a in range(4) if p[a] == a) for p in itertools.permutations(range(4))})
check("C2 Unter der Hebung durch Permutationsmatrizen ist Code 1 S4-stabil: Spiegelungsdarstellung 3 + Invariante 1",
      inv and chars == [0, 1, 2, 4])
n12 = sp.Matrix([[0, 1, 0, 0, 0], [-1, 0, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]])  # Tits-Hebung
ns_ = n12 * s * n12.inv()
check("C3 Unter der Tits-Hebung (Vorzeichen, Determinante 1) ist die Summe nicht invariant: 'konkrete Hebung' ist eine Wahl",
      ns_ != s and ns_ != -s)
check("C4 Auf Code 1 ist Twist^2 keine Familienpermutation (alle 24 lassen die Summe fest, Twist^2 kehrt sie um)",
      T(T(s)) == -s and inv)
# Code 2: verträglich mit einer Transposition P = (12) als Tausch
P12 = perm((1, 0, 2, 3))
d = E(1, 5) - E(2, 5)
fix = [E(3, 3) - E(4, 4), E(1, 1) + E(2, 2) - E(3, 3) - E(4, 4), E(1, 2) + E(2, 1)]
ok = (P12 * d * P12.T == -d and T(d) == -sp.I * d
      and all(P12 * x * P12.T == x and T(x) == x for x in fix))
check("C5 Code 2 = (E_15 - E_25) + drei (12)-feste Ströme aus gl(4): Twist = diag(1,1,1,-i) und Twist^2 = Tausch (12)", ok)
# Wahlfreiheit: Dimension des (12)-festen Teils im Twist-Eigenraum 1 und des (12)-ungeraden Teils im Eigenraum -i
basis = [E(a, b) for a in range(1, 6) for b in range(1, 6)]


def dim_eig(tw_val, p_val):
    rows = []
    for x in basis:
        if T(x) != tw_val * x:
            continue
        y = P12 * x * P12.T
        rows.append(list((x + p_val * y).reshape(1, 25)))
    return sp.Matrix(rows).rank() if rows else 0


d1 = dim_eig(1, 1) - 1                                # ohne die Spur (gl(5) -> sl(5))
d2 = dim_eig(-sp.I, -1)
ZAHLEN["Code 2: Dimension der Auswahlräume (Twist 1 und (12)-fest; Twist -i und (12)-ungerade)"] = [d1, d2]
say(f"   Auswahl für Code 2: drei Vektoren aus einem {d1}-dimensionalen Raum, einer aus einem {d2}-dimensionalen")
check("C6 Code 2 ist nicht eindeutig: eine stetige Familie von Wahlen (Gr(3, 10) mal P^1)", d1 == 10 and d2 == 2)
# Ein echtes Gruppenelement mit demselben Spektrum: die Tits-Hebung der Transposition (Ordnung 4)
A = sp.zeros(25, 25)
for c, x in enumerate(basis):
    y = n12 * x * n12.inv()
    for r_ in range(25):
        A[r_, c] = y[r_ // 5, r_ % 5]
ev = {str(k): v for k, v in A.eigenvals().items()}
ZAHLEN["Eigenwerte der Tits-Hebung auf gl(5)"] = dict(sorted(ev.items()))
check("C7 Die Tits-Hebung der Transposition ist eine Automorphie der Ordnung 4 mit Eigenwerten 1, -1, i, -i (11, 2, 6, 6): "
      "das Spektrum (1,1,1,-i) zeichnet den Twist nicht aus", ev == {"1": 11, "-1": 2, "I": 6, "-I": 6})
# Welches Zweiqubit-Gatter? Operator-Schmidt-Rang von diag(1,1,1,-i) je nach Lage des Eigenvektors zu -i
def schmidt(U):
    R = sp.zeros(4, 4)
    for a, b, c, e in itertools.product(range(2), repeat=4):
        R[2 * a + c, 2 * b + e] = U[2 * a + b, 2 * c + e]
    return R.rank()


Pswap = sp.Matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
Uswap = (sp.eye(4) + sp.I * Pswap) / (1 + sp.I)
Uctrl = sp.diag(1, 1, 1, -sp.I)
check("C8 (I + iP)/(1 + i) hat die Eigenwerte (1,1,1,-i), U^2 = P",
      sorted(map(str, Uswap.eigenvals().items())) == sorted(map(str, {sp.Integer(1): 3, -sp.I: 1}.items()))
      and sp.simplify(Uswap * Uswap - Pswap) == sp.zeros(4, 4))
ZAHLEN["Operator-Schmidt-Rang: Wurzel aus Tausch, gesteuerte Uhr"] = [schmidt(Uswap), schmidt(Uctrl)]
check("C9 Gleiches Spektrum, verschiedene Gatter: Wurzel aus Tausch hat Schmidt-Rang 4, die gesteuerte Uhr Rang 2",
      schmidt(Uswap) == 4 and schmidt(Uctrl) == 2)
befund("C", "E",
       "Behauptung (2) stimmt als Existenzaussage: Es gibt vierdimensionale Räume im A4-Sektor, auf denen der Twist "
       "diag(1,1,1,-i) ist, auch mit Twist^2 = Familientransposition. Der Raum ist aber gewählt (stetige Familie), die "
       "Hebung ist gewählt, und die Qubit-Zerlegung kommt nicht aus E8. Mit der Kleinschen Vierergruppe als "
       "Qubit-Struktur wäre dasselbe Spektrum die gesteuerte Uhr, nicht die Wurzel aus dem Tausch.")

# ---------------------------------------------------------------------------------------------------------------------
say("\nTeil D  Geometrie der Unterraumfamilie")
r = sp.symbols("r", positive=True)
y = sp.symbols("y1:5", real=True)
x5 = sp.sqrt(r**2 - sum(v**2 for v in y))
X = sp.Matrix(list(y) + [x5])
Pm = sp.eye(5) - X * X.T / r**2                      # Projektor auf die Hyperebene x^perp
dP = [Pm.diff(v) for v in y]
null = {v: 0 for v in y}
P0 = Pm.subs(null)
dP0 = [m.subs(null) for m in dP]
ok = True
anz = 0
for m_, n_ in itertools.product(range(4), repeat=2):
    F = P0 * (dP0[m_] * dP0[n_] - dP0[n_] * dP0[m_]) * P0
    for q in range(4):
        e = sp.zeros(5, 1)
        e[q] = 1
        soll = sp.zeros(5, 1)
        if q == n_:
            soll[m_] += 1 / r**2
        if q == m_:
            soll[n_] -= 1 / r**2
        ok &= sp.simplify(F * e - soll) == sp.zeros(5, 1)
        anz += 1
check(f"D1 Projizierte Krümmung F_mn = P[d_m P, d_n P]P: alle {anz} Komponentenwirkungen sind die der runden S^4", ok and anz == 64)
ric = sp.zeros(4, 4)
for m_, q in itertools.product(range(4), repeat=2):
    tot = 0
    for n_ in range(4):
        F = P0 * (dP0[n_] * dP0[m_] - dP0[m_] * dP0[n_]) * P0
        e = sp.zeros(5, 1)
        e[q] = 1
        tot += (F * e)[n_]
    ric[m_, q] = sp.simplify(tot)
check("D2 Ricci = (3/r^2) g: das ist das Tangentialbündel der S^4 mit Levi-Civita-Verbindung, gültig für jedes R^5",
      ric == 3 / r**2 * sp.eye(4))
# Vorzeichen = reelle Form: Hyperboloid x^T eta x = -r^2 mit eta = diag(1,1,1,1,-1)
eta = sp.diag(1, 1, 1, 1, -1)
Xh = sp.Matrix(list(y) + [sp.sqrt(r**2 + sum(v**2 for v in y))])
Ph = sp.eye(5) + Xh * Xh.T * eta / r**2              # eta-orthogonaler Projektor auf den Tangentialraum
dPh = [Ph.diff(v).subs(null) for v in y]
Ph0 = Ph.subs(null)
rich = sp.zeros(4, 4)
for m_, q in itertools.product(range(4), repeat=2):
    tot = 0
    for n_ in range(4):
        F = Ph0 * (dPh[n_] * dPh[m_] - dPh[m_] * dPh[n_]) * Ph0
        e = sp.zeros(5, 1)
        e[q] = 1
        tot += (F * e)[n_]
    rich[m_, q] = sp.simplify(tot)
check("D3 Mit der reellen Form der Signatur (4,1) gibt dieselbe Konstruktion Ricci = -(3/r^2) g: Das Vorzeichen ist die Wahl "
      "der reellen Struktur", rich == -3 / r**2 * sp.eye(4))
befund("D", "E; Vorzeichen und Radius C",
       "Behauptung (3) stimmt als Rechnung, trägt aber keine E8-Information: Die Hyperebenen jedes R^5 bilden RP^4, die "
       "projizierte Verbindung ist die Levi-Civita-Verbindung der runden S^4. Das positive Vorzeichen folgt aus der "
       "Wahl der kompakten reellen Form so(5); nativ ist in H6 so(3,2) (z89). Der Radius ist frei, die volle "
       "Cartan-Krümmung bleibt null (wie P_curv im Bestand).")

# ---------------------------------------------------------------------------------------------------------------------
say("\nTeil E  Projektion der 240 Wurzeln und Isotropie")
proj = {}
for rt in roots:
    proj[rt[5:]] = proj.get(rt[5:], 0) + 1
nach_norm = {}
for p, w in proj.items():
    nach_norm[dot(p, p)] = nach_norm.get(dot(p, p), 0) + w
ZAHLEN["Projektion auf A3: Normquadrat -> Zahl der Wurzeln"] = {str(k): v for k, v in sorted(nach_norm.items())}
check("E1 27 Punkte; Gewichte 40 (Norm 0), 128 (sqrt(3)/2), 60 (1), 12 (sqrt(2))",
      len(proj) == 27 and nach_norm == {Fr(0): 40, Fr(3, 4): 128, Fr(1): 60, Fr(2): 12})
k = sp.symbols("k1:9", real=True)


def moment(n, idx):
    tot = 0
    for rt in roots:
        tot += sum(sp.Rational(rt[i].numerator, rt[i].denominator) * k[i] for i in idx) ** n
    return sp.expand(tot)


def kk(idx):
    return sum(k[i] ** 2 for i in idx)


a3 = (5, 6, 7)
koef = {}
ok = True
for n, c in ((2, 60), (4, 36), (6, 30)):
    ok &= sp.expand(moment(n, a3) - c * kk(a3) ** (n // 2)) == 0
    koef[n] = sp.Rational(c, 240) / sp.factorial(n) * (-1) ** (n // 2)
ZAHLEN["Fourierkoeffizienten k^2, k^4, k^6"] = [str(koef[2]), str(koef[4]), str(koef[6])]
check("E2 Fourierreihe 1 - k^2/8 + k^4/160 - k^6/5760, drehinvariant bis zur sechsten Ordnung",
      ok and [koef[2], koef[4], koef[6]] == [-sp.Rational(1, 8), sp.Rational(1, 160), -sp.Rational(1, 5760)])
m8 = moment(8, a3)
rest8 = sp.expand(m8 - sp.Poly(m8, *[k[i] for i in a3]).coeff_monomial(k[5] ** 8) * kk(a3) ** 4)
check("E3 In achter Ordnung ist die Projektion auf A3 richtungsabhängig", rest8 != 0)
alle = tuple(range(8))
ok = True
for n, c in ((2, 60), (4, 36), (6, 30)):
    ok &= sp.expand(moment(n, alle) - c * kk(alle) ** (n // 2)) == 0
check("E4 Dieselben Identitäten gelten in allen acht Richtungen: sum (k.alpha)^n = c_n |k|^n für n = 2, 4, 6 (c = 60, 36, 30)", ok)
bessel = [sp.Rational(-1, 2) ** n / (sp.factorial(n) * sp.rf(4, n)) for n in (1, 2, 3)]
check("E5 Die Reihe ist das Mittel über die Sphäre S^7 vom Radius sqrt(2): 0F1(; 4; -k^2/2)",
      bessel == [koef[2], koef[4], koef[6]])
d5 = (0, 1, 2, 3, 4)
gem = (0, 3, 6)                                       # ein gemischter dreidimensionaler Unterraum (je eine Koordinate)
ok = all(sp.expand(moment(n, idx) - c * kk(idx) ** (n // 2)) == 0
         for idx in (d5, gem) for n, c in ((2, 60), (4, 36), (6, 30)))
check("E6 Folge: Die Projektion auf jeden Unterraum (hier D5-Raum und ein gemischter 3-Raum) ist gleich isotrop; "
      "die Eigenschaft wählt weder A3 noch drei Dimensionen", ok)
befund("E", "E",
       "Behauptung (4) stimmt. Der Grund ist allgemein: Die Weylgruppe von E8 hat außer der Norm keine Invariante "
       "unter Grad 8, die Wurzeln sind ein sphärisches 7-Design. Jede Projektion, in jede Dimension, ist bis zur "
       "siebten Ordnung isotrop. Das ist eine gute Eigenschaft für jeden E8-Transfer, aber kein Hinweis auf A3 als Raum.")

# ---------------------------------------------------------------------------------------------------------------------
say("\nTeil F  Bilanz")
bilanz = [
    ("E", "Zentralisator", "20 neutrale Wurzeln = A4 = zweite SU(5) in E8 > SU(5) x SU(5); 12 Familienwurzeln und 8 "
     "Klebewurzeln (nu_R-Richtung mal 4, mit B - L)."),
    ("E", "Gatterspektrum", "Der Twist der relativen Drehung hat auf den Strömen die Eigenwerte 1 (112 Wurzeln) und -i "
     "(128 Klebewurzeln); sein Quadrat ist eine innere Automorphie. Das ist der Inhalt von z145, jetzt auf Stromebene."),
    ("E", "Code", "Vierdimensionale Räume mit diag(1,1,1,-i) existieren, auch mit Twist^2 = Transposition. Sie sind "
     "gewählt (stetige Familie), die Hebung ist gewählt, die Qubit-Zerlegung ist nicht aus E8."),
    ("E", "Kein Alleinstellungsmerkmal", "Die Tits-Hebung einer Familientransposition ist eine echte Automorphie mit "
     "demselben Spektrum auf geeigneten Räumen. Wurzel aus Tausch und gesteuerte Uhr teilen das Spektrum."),
    ("E; C", "Krümmung", "R_mn = 3/r^2 g ist die runde S^4 jedes R^5. Vorzeichen = Wahl der kompakten reellen Form, "
     "Radius frei, Cartan-Krümmung null."),
    ("E", "Isotropie", "Gilt bis zur siebten Ordnung für jede Projektion von E8; wählt weder A3 noch drei Dimensionen."),
    ("O", "Härtetest der Analyse", "Die Analyse nennt selbst, dass der Zweinahtsektor (1, 2, 1 je Ladung) keinen "
     "vierdimensionalen Raum fester Ladung hat. Das ist der Zweiqubit-Raum, nach Ladung zerlegt; die acht neutralen "
     "Klebeströme sind standardmodellneutral und können ihn deshalb nicht tragen. Hier nicht neu gerechnet."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
urteil = ("Alle gerechneten Aussagen der zwanzigsten Analyse stimmen; keine von ihnen ersetzt eine Prämisse. (1) Der "
          "neutrale A4-Sektor ist die zweite SU(5) der klassischen Zerlegung E8 > SU(5) x SU(5). Die acht "
          "Klebewurzeln sind die Richtung des rechtshändigen Neutrinos mal der Familien-4 und tragen B - L. (2) Der "
          "Eigenwert -i des Gatters ist die Twistphase der Klebeklasse aus z145, jetzt auf der Ebene der 248 Ströme: "
          "1 auf 120, -i auf 128. Sein Quadrat ist die innere Involution, die D8 in E8 auszeichnet. Ein "
          "vierdimensionaler Code mit diag(1,1,1,-i) und U^2 = Transposition existiert, ist aber aus einer stetigen "
          "Familie gewählt; die Hebung der Permutationen und die Zerlegung in zwei Qubits sind Zutaten. Mit der "
          "natürlichen Qubit-Struktur der Vierergruppe wäre dasselbe Spektrum die gesteuerte Uhr. Die algebraische "
          "Lücke zwischen Kleber und Gatter ist damit nicht enger als nach z145: Phase ja, Prozess nein. (3) Die "
          "positive Krümmung ist die runde S^4 der Hyperebenen eines R^5 und gilt für jeden fünfdimensionalen Raum. "
          "Das Vorzeichen ist die Wahl der kompakten reellen Form; mit Signatur (4,1) gibt dieselbe Rechnung das "
          "andere Vorzeichen. Nativ trägt H6 so(3,2). Radius und dynamisches Feld fehlen, wie die Analyse selbst sagt. "
          "(4) Die Isotropie bis zur sechsten Ordnung ist eine Eigenschaft von E8 insgesamt (keine Weyl-Invariante "
          "zwischen Grad 2 und 8) und gilt für jede Projektion. Sie ist nützlich für jeden E8-Transfer, zeichnet "
          "aber weder A3 noch drei Raumdimensionen aus. (5) Neu und brauchbar: der Twist als nicht multiplikative "
          "Wurzel einer inneren Automorphie (B5), die B - L-Ladung der neutralen Klebeströme (A5) und das 7-Design "
          "(E4). Der Vorschlag 'E8 ist der Zustandsraum, die relative Bewegung der Faktoren ist der Prozess' bleibt "
          "eine Lesart; sie wird in z149 auf Operatorebene geprüft. K12, K13, K15, K16 und die neun Prämissen sind "
          "unverändert. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z153_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z153.txt").write_text("\n".join(OUT) + "\n")
