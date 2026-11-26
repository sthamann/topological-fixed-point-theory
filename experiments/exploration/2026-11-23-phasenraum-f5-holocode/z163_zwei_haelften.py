"""Z163: Erste Entscheidungsprüfung des Wegs 'zwei Hälften' (Naht-Vertex aus H6 + selbstduale Twistor-Hälfte des Ledgers).

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart, N numerisch, O offen.

Anlass: Die Durchsicht des ganzen Bestands (Ledger, Paper, Atlas, H6) zeigt zwei Programme, die sich nie getroffen haben.
  Ledger  CELEST.SEAM.01 / SEAM.EQUIV.TWISTOR.01: holomorphe (selbstduale) Theorie auf dem Twistorraum, Eichgruppe E8,
          Axion aus der Anomalie, Niveau 1; ausdrücklich nur die selbstduale Hälfte.
  H6      z110 bis z112, z123: Die Naht auf einer Twistorgeraden liefert die Vertexhälfte (Ströme, Niveau 1); es fehlt die
          kinetische Hälfte (Prämisse P_BF) und der Zähler <ij>^4 (Prämisse P_N4).
Die bekannte Struktur, die beides verbindet, ist die Twistorwirkung: holomorpher Teil auf dem Twistorraum plus
Kopplung mal Integral über die Geraden L_x eines Funktionals der Geraden (log det) [L, Mason 2005,
Boels-Mason-Skinner 2007]. Diese Rechnung prüft, ob die Teile zueinander passen, mit exakten, billigen Prüfungen.

Teil A  Anomaliebedingung der selbstdualen Hälfte: Sp_adj X^4 proportional (X.X)^2. Welche Algebren des Bestands
        erfüllen sie?
Teil B  Uhr-Graduierung von E8 und Statistik: Können die Klassen 1 und 3 (die 16) in der Twistor-Hälfte Fermionen sein?
Teil C  Vertexhälfte: zusammenhängende Stromkorrelatoren bei Niveau k (Ward-Rekursion). Einfachspur, linear in k,
        Zähler aus dem Gewicht -4.
Teil D  Wirkungsleiter des Ledgers: wie viele unabhängige Stufen, und was sagt die duale Lesart (gegen Zahlenspielerei
        markiert)?
Teil E  Bilanz und Urteil.
"""
import itertools
import json
import math
import os
import time
from fractions import Fraction as Fr
from functools import lru_cache
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # deterministische Rundung (einfädig)

import numpy as np
import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
PI = math.pi


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val, text):
    BEFUNDE[name] = val
    say(f"   -> {name}: {text}")


# ------------------------------------------------------------------------------------------------------------------
say("Z163: Entscheidungsprüfung 'zwei Hälften' (Naht-Vertex + selbstduale Twistor-Hälfte)")
say("=" * 100)
say("\nTeil A  Anomaliebedingung der selbstdualen Hälfte")
say("   Lesart [L, Costello 2021; Costello-Paquette 2022, aus der Literatur, hier nicht nachgerechnet]: Die holomorphe")
say("   Eichtheorie auf dem Twistorraum ist auf Einschleifenniveau nur konsistent, wenn Sp_adj X^4 = lambda^2 (X.X)^2")
say("   gilt (kein unabhängiger quartischer Casimir); der Rest wird von einem Axion aufgehoben. Mit Materie R zählt")
say("   Sp_adj X^4 - Sp_R X^4 (Weyl-Fermionen) + 1/2 Sp X^4 (reelle Skalare). Exakt geprüft wird die Algebra.")


def e8_roots():
    out = []
    for i, j in itertools.combinations(range(8), 2):
        for si in (1, -1):
            for sj in (1, -1):
                w = [Fr(0)] * 8
                w[i], w[j] = Fr(si), Fr(sj)
                out.append(tuple(w))
    for eps in itertools.product((1, -1), repeat=8):
        if eps.count(-1) % 2 == 0:
            out.append(tuple(Fr(e, 2) for e in eps))
    return sorted(out)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def span_basis(vecs):
    """Rationale Basis des Spanns (Gauß, exakt)."""
    basis, red = [], []
    for v in vecs:
        w = list(v)
        for r, piv in red:
            if w[piv] != 0:
                f = w[piv] / r[piv]
                w = [a - f * b for a, b in zip(w, r)]
        if any(w):
            piv = next(i for i, a in enumerate(w) if a != 0)
            red.append((w, piv))
            basis.append(v)
    return basis


def quartic(roots):
    """P = Summe (alpha.X)^4 und Q = X.X als Polynome in den Koordinaten des Cartan-Spanns."""
    B = span_basis(roots)
    cs = sp.symbols(f"c0:{len(B)}")
    P = sp.Poly(0, *cs)
    for a in roots:
        lin = sp.Poly(sum(sp.Rational(dot(a, b).numerator, dot(a, b).denominator) * c for b, c in zip(B, cs)), *cs)
        P = P + lin ** 4
    Q = sp.Poly(sum(sp.Rational(dot(bi, bj).numerator, dot(bi, bj).denominator) * ci * cj
                    for bi, ci in zip(B, cs) for bj, cj in zip(B, cs)), *cs)
    return P, Q, cs


def square_test(roots):
    P, Q, cs = quartic(roots)
    lam = sp.Rational(P.coeff_monomial(cs[0] ** 4), Q.coeff_monomial(cs[0] ** 2) ** 2)
    return (P - lam * Q ** 2).is_zero, lam


E8 = e8_roots()
e = lambda i: tuple(Fr(1) if k == i else Fr(0) for k in range(8))
sub = lambda a, b: tuple(x - y for x, y in zip(a, b))
r0, r1 = sub(e(6), e(7)), sub(e(5), e(6))
is_int = lambda w: all(x.denominator == 1 for x in w)
supp = lambda w: {i for i, x in enumerate(w) if x != 0}
ALG = {
    "A1 = su(2)": [r0, tuple(-x for x in r0)],
    "A2 = su(3)": [sub(e(i), e(j)) for i in range(3) for j in range(3) if i != j],
    "A3 = su(4)": [sub(e(i), e(j)) for i in range(4) for j in range(4) if i != j],
    "A4 = su(5)": [sub(e(i), e(j)) for i in range(5) for j in range(5) if i != j],
    "D4 = so(8)": [w for w in E8 if is_int(w) and supp(w) <= {0, 1, 2, 3}],
    "D5 = so(10)": [w for w in E8 if is_int(w) and supp(w) <= {0, 1, 2, 3, 4}],
    "E6": [w for w in E8 if dot(w, r0) == 0 and dot(w, r1) == 0],
    "E7": [w for w in E8 if dot(w, r0) == 0],
    "E8": E8,
}
DELIGNE = {"A1 = su(2)": (2, 3), "A2 = su(3)": (3, 8), "D4 = so(8)": (6, 28), "E6": (12, 78), "E7": (18, 133),
           "E8": (30, 248)}
check("A0 Wurzelzahlen: A1 2, A2 6, A3 12, A4 20, D4 24, D5 40, E6 72, E7 126, E8 240",
      [len(ALG[k]) for k in ALG] == [2, 6, 12, 20, 24, 40, 72, 126, 240])
res = {}
for name in ALG:
    ok, lam = square_test(ALG[name])
    res[name] = (ok, lam)
    extra = ""
    if name in DELIGNE:
        h, d = DELIGNE[name]
        extra = f"; 10 h^2/(dim+2) = {sp.Rational(10 * h * h, d + 2)}, h + 6 = {h + 6}"
    say(f"   {name:12s} Sp_adj X^4 = lambda^2 (X.X)^2: {'ja' if ok else 'nein'}"
        + (f", lambda^2 = {lam}" if ok else "") + extra)
check("A1 Die Bedingung gilt genau für A1, A2, D4, E6, E7, E8 (von den geprüften neun Algebren)",
      [k for k in ALG if res[k][0]] == list(DELIGNE))
check("A2 lambda^2 = 10 h^2/(dim + 2) = h + 6 auf dieser Reihe: 8, 9, 12, 18, 24, 36",
      all(res[k][1] == sp.Rational(10 * h * h, d + 2) == h + 6 for k, (h, d) in DELIGNE.items()))
check("A3 E8: lambda^2 = 36 (der Wert im Vertrag CELEST.SEAM.01, WP5e)", res["E8"][1] == 36)
check("A4 so(10), su(4), su(5) erfüllen die Bedingung nicht (echter quartischer Casimir)",
      not res["D5 = so(10)"][0] and not res["A3 = su(4)"][0] and not res["A4 = su(5)"][0])
sm = ALG["A2 = su(3)"] + [sub(e(3), e(4)), sub(e(4), e(3))]
check("A5 su(3) + su(2) (Standardmodell): kein einzelnes Quadrat (9 s^2 + 8 t^2, kein gemischter Term)",
      not square_test(sm)[0])
befund("A_anomalie", {k: [bool(v[0]), str(v[1])] for k, v in res.items()},
       "Die selbstduale Twistor-Hälfte verlangt im Bestand E8 selbst. Spin(10), SU(5), Pati-Salam und das "
       "Standardmodell erfüllen die Bedingung nicht; für E8 ist lambda^2 = 36 = h + 6 [E; Bedeutung L].")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil B  Uhr-Graduierung von E8 und Statistik")
say("   Uhr-Klasse einer Wurzel w: 4 h.w mod 4 mit h = (1/2,1/2,1/2,1/2,1/2;0,0,0) (Verschiebung des Vertrags, WP1).")
cls = lambda w: int(2 * sum(w[:5])) % 4
KL = {k: [w for w in E8 if cls(w) == k] for k in range(4)}
check("B1 Klassen (52, 64, 60, 64); mit dem Cartan-Teil (60, 64, 60, 64) wie im Vertrag (WP1)",
      [len(KL[k]) for k in range(4)] == [52, 64, 60, 64])
check("B2 Klasse 0 = so(10) + su(4) (40 + 12 Wurzeln); Klasse 3 = -Klasse 1; Klasse 2 = -Klasse 2",
      sorted(KL[0]) == sorted(w for w in E8 if is_int(w) and (supp(w) <= {0, 1, 2, 3, 4} or supp(w) <= {5, 6, 7}))
      and sorted(KL[3]) == sorted(tuple(-x for x in w) for w in KL[1])
      and sorted(KL[2]) == sorted(tuple(-x for x in w) for w in KL[2]))
xs = sp.symbols("x1:6")
ys = sp.symbols("y1:4")
VV = list(xs) + list(ys)


def tr4(ws):
    P = sp.Poly(0, *VV)
    for w in ws:
        P = P + sp.Poly(sum(sp.Rational(c.numerator, c.denominator) * v for c, v in zip(w, VV)), *VV) ** 4
    return P


p2, p4 = sp.Poly(sum(v ** 2 for v in xs), *VV), sp.Poly(sum(v ** 4 for v in xs), *VV)
q2, q4 = sp.Poly(sum(v ** 2 for v in ys), *VV), sp.Poly(sum(v ** 4 for v in ys), *VV)
NAMES = ["p4", "p2^2", "q4", "q2^2", "p2 q2"]


def dec(P):
    """Zerlegung in p4, p2^2, q4, q2^2, p2 q2 (p: so(10)-Teil, q: su(4)-Teil); Rest muss null sein."""
    co = lambda *ex: P.coeff_monomial(sp.prod(v ** k for v, k in zip(VV, ex)))
    b = sp.Rational(co(2, 2, 0, 0, 0, 0, 0, 0), 2)
    a = co(4, 0, 0, 0, 0, 0, 0, 0) - b
    d = sp.Rational(co(0, 0, 0, 0, 0, 2, 2, 0), 2)
    c = co(0, 0, 0, 0, 0, 4, 0, 0) - d
    m = co(2, 0, 0, 0, 0, 2, 0, 0)
    rest = P - (a * p4 + b * p2 ** 2 + c * q4 + d * q2 ** 2 + m * p2 * q2)
    assert rest.is_zero
    return [a, b, c, d, m]


TK = {k: tr4(KL[k]) for k in range(4)}
DK = {k: dec(TK[k]) for k in range(4)}
for k in range(4):
    say(f"   Klasse {k} ({len(KL[k]):2d} Wurzeln): " + ", ".join(f"{n}: {v}" for n, v in zip(NAMES, DK[k])))
tot = [sum(DK[k][i] for k in range(4)) for i in range(5)]
say("   Summe (alle 240 als Eichbosonen): " + ", ".join(f"{n}: {v}" for n, v in zip(NAMES, tot)))
check("B3 Summe = 36 (p2 + q2)^2: p4 und q4 heben sich nur über alle vier Klassen auf (4 - 8 + 12 - 8 = 0)",
      tot == [0, 36, 0, 36, 72] and [DK[k][0] for k in range(4)] == [4, -8, 12, -8])
var = {}
for wname, w2 in (("Klasse 2 reeller Skalar (1/2)", sp.Rational(1, 2)), ("Klasse 2 Eichboson (1)", 1),
                  ("Klasse 2 fehlt (0)", 0), ("Klasse 2 Weyl-Fermion (-1)", -1)):
    v = [DK[0][i] - DK[1][i] + w2 * DK[2][i] for i in range(5)]
    var[wname] = v
    say(f"   Uhr-Statistik, Klasse 0 Eichboson, Klasse 1 = (16,4) Weyl, {wname}: "
        + ", ".join(f"{n}: {x}" for n, x in zip(NAMES, v)))
is_sq = lambda v: v[0] == 0 and v[2] == 0 and v[4] ** 2 == 4 * v[1] * v[3]
check("B4 Keine Uhr-Statistik gibt ein Quadrat: mit der 16 als Fermion bleibt ein quartischer Rest",
      not any(is_sq(v) for v in var.values()))
check("B5 Spin(10) mit n Weyl-16 und m reellen 10: p4-Koeffizient 4 + 2 n + m > 0, nie null",
      DK[0][0] == 4 and sp.Rational(DK[1][0], 4) == -2 and sp.Rational(DK[2][0], 6) == 2)
# Eichung der Vorzeichenkonvention an zwei bekannten Lösungen [L]: SU(2) mit 8 und SU(3) mit 9 Flavours
a_, b_ = sp.symbols("a b")
su2 = {"adj": 2 * (2 * a_) ** 4, "F": 2 * a_ ** 4}
w3 = [a_, b_, -a_ - b_]
su3 = {"adj": sp.expand(sum((w3[i] - w3[j]) ** 4 for i in range(3) for j in range(3) if i != j)),
       "F": sp.expand(sum(w ** 4 for w in w3))}
check("B6 Eichung der Konvention: Sp_adj X^4 = 16 Sp_F X^4 für SU(2), = 18 Sp_F X^4 für SU(3) (N_f = 8 bzw. 9) [E]",
      sp.expand(su2["adj"] - 16 * su2["F"]) == 0 and sp.expand(su3["adj"] - 18 * su3["F"]) == 0)
befund("B_statistik", {k: [str(x) for x in v] for k, v in var.items()},
       "In der Twistor-Hälfte müssen alle 248 Komponenten Eichbosonen sein. Die 16 kann dort nicht als Fermion in "
       "den Klassen 1 und 3 sitzen, und Spin(10) mit 16ern ist für keine Zahl von Familien und Higgs-10ern "
       "anomaliefrei. Die chirale Materie (Tor T4) braucht einen eigenen Ort: die Naht selbst (16 Majoranas, "
       "nu(Uhr) = 16), nicht das Bündel [E für die Algebra, L für die Folgerung].")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil C  Vertexhälfte: zusammenhängende Stromkorrelatoren bei Niveau k")
say("   Ward-Rekursion zur OPE J^a(z) J^b(w) = k delta/(z-w)^2 + i f^abc J^c(w)/(z-w); Erzeuger mit Sp(T^a T^b) = delta.")


def gens(N):
    T = []
    for i in range(N):
        for j in range(i + 1, N):
            m = np.zeros((N, N), complex)
            m[i, j] = m[j, i] = 1 / math.sqrt(2)
            T.append(m)
            m = np.zeros((N, N), complex)
            m[i, j], m[j, i] = -1j / math.sqrt(2), 1j / math.sqrt(2)
            T.append(m)
    for d in range(1, N):
        m = np.zeros((N, N), complex)
        for i in range(d):
            m[i, i] = 1
        m[d, d] = -d
        T.append(m / math.sqrt(d * (d + 1)))
    return T


def make(N):
    T = gens(N)
    D = len(T)
    F = np.zeros((D, D, D))
    for a in range(D):
        for b in range(D):
            cm = T[a] @ T[b] - T[b] @ T[a]
            for c in range(D):
                F[a, b, c] = (-1j * np.trace(cm @ T[c])).real
    nz = {(a, b): [(c, F[a, b, c]) for c in range(D) if abs(F[a, b, c]) > 1e-12] for a in range(D) for b in range(D)}

    @lru_cache(maxsize=None)
    def full(idx, zs, k):
        n = len(idx)
        if n == 0:
            return 1.0 + 0j
        if n == 1:
            return 0j
        a1, z1 = idx[0], zs[0]
        tot = 0j
        for j in range(1, n):
            zj = zs[j]
            if idx[j] == a1:
                tot += k / (z1 - zj) ** 2 * full(idx[1:j] + idx[j + 1:], zs[1:j] + zs[j + 1:], k)
            for c, f in nz[(a1, idx[j])]:
                tot += 1j * f / (z1 - zj) * full(idx[1:j] + (c,) + idx[j + 1:], zs[1:], k)
        return tot

    @lru_cache(maxsize=None)
    def conn(idx, zs, k):
        n = len(idx)
        if n < 2:
            return 0j
        tot = full(idx, zs, k)
        rest = list(range(1, n))
        for r in range(1, n - 2):                      # Block mit dem ersten Punkt hat r + 1 Elemente, Rest >= 2
            for comb in itertools.combinations(rest, r):
                blk = (0,) + comb
                oth = tuple(i for i in rest if i not in comb)
                tot -= (conn(tuple(idx[i] for i in blk), tuple(zs[i] for i in blk), k)
                        * full(tuple(idx[i] for i in oth), tuple(zs[i] for i in oth), k))
        return tot

    def trace_formula(idx, zs):
        n = len(idx)
        tot = 0j
        for perm in itertools.permutations(range(1, n)):
            order = (0,) + perm
            m = np.eye(N, dtype=complex)
            den = 1.0 + 0j
            for s in range(n):
                m = m @ T[idx[order[s]]]
                den *= zs[order[s]] - zs[order[(s + 1) % n]]
            tot += np.trace(m) / den
        return -tot

    return full, conn, trace_formula, D


rng = np.random.default_rng(163)
ZS = tuple(complex(round(x, 6), round(y, 6)) for x, y in rng.uniform(-2, 2, size=(6, 2)))
worst_lin, worst_tr, n_cases = 0.0, 0.0, 0
ratio_multi = None
for N, nmax, ntup in ((2, 6, 6), (3, 5, 6)):
    full, conn, trf, D = make(N)
    for n in range(3, nmax + 1):
        tuples = [tuple(int(v) for v in rng.integers(0, D, size=n)) for _ in range(ntup)]
        tuples.append(tuple((i % 2) for i in range(n)) if N == 2 else tuple((i % 3) for i in range(n)))
        for idx in tuples:
            zs = ZS[:n]
            c1 = conn(idx, zs, 1.0)
            scale = max(1.0, abs(c1))
            for k in (2.0, 3.0):
                worst_lin = max(worst_lin, abs(conn(idx, zs, k) - k * c1) / scale)
            worst_tr = max(worst_tr, abs(c1 - trf(idx, zs)) / scale)
            n_cases += 1
    if N == 2:
        idx4, zs4 = (0, 0, 1, 1), ZS[:4]
        f1, c1 = full(idx4, zs4, 1.0), conn(idx4, zs4, 1.0)
        ratio_multi = abs(f1 - c1) / abs(c1)
        d3 = full(idx4, zs4, 3.0) - 3 * full(idx4, zs4, 2.0) + 3 * full(idx4, zs4, 1.0) - full(idx4, zs4, 0.0)
        d2 = full(idx4, zs4, 2.0) - 2 * full(idx4, zs4, 1.0) + full(idx4, zs4, 0.0)
ZAHLEN.update(faelle=n_cases, abw_linear=float(f"{worst_lin:.1e}"), abw_spur=float(f"{worst_tr:.1e}"),
              mehrspur_zu_einfachspur_k1=round(ratio_multi, 4))
say(f"   {n_cases} Fälle (su(2) bis 6 Punkte, su(3) bis 5 Punkte); größte Abweichung von der Linearität in k: "
    f"{worst_lin:.1e}, von der Einfachspurformel: {worst_tr:.1e}")
check("C1 Zusammenhängender Teil ist exakt linear in k: W_n(k) = k W_n(1) für n = 3..6", worst_lin < 1e-9)
check("C2 W_n(1) = -Summe über Ordnungen Sp(T...T)/(z_12 z_23 ... z_n1): reine Einfachspur", worst_tr < 1e-9)
check("C3 Der volle Korrelator hat bei 4 Punkten einen k^2-Teil (Mehrspur), der zusammenhängende nicht",
      abs(d3) < 1e-9 and abs(d2) > 1e-6)
say(f"   Bei k = 1, Indizes (1,1,2,2): |Mehrspurteil| / |Einfachspurteil| = {ratio_multi:.4f} (nicht unterdrückt)")
worst_dec = 0.0
for n in (3, 4, 5, 6):
    s = 0j
    for perm in itertools.permutations(range(1, n)):
        order = (0,) + perm
        den = 1.0 + 0j
        for t in range(n):
            den *= ZS[order[t]] - ZS[order[(t + 1) % n]]
        s += 1 / den
    worst_dec = max(worst_dec, abs(s))
check("C4 Summe aller Ordnungen des Nenners ist null (n = 3..6): vollsymmetrische Invarianten (d-Tensoren) fallen "
      "heraus, es zählt nur f und k; das gilt damit auch für (E8)_1 ohne Fermiondarstellung", worst_dec < 1e-9)
# Gewichtszählung: Zähler <ij>^p mal Parke-Taylor-Nenner; Homogenität je Bein
n5 = 5
us, vs = sp.symbols(f"u1:{n5 + 1}"), sp.symbols(f"v1:{n5 + 1}")
br = lambda i, j: us[i] * vs[j] - vs[i] * us[j]
pt = 1 / sp.prod(br(i, (i + 1) % n5) for i in range(n5))
psym = sp.Symbol("p")
euler = lambda f, i: sp.simplify((us[i] * sp.diff(f, us[i]) + vs[i] * sp.diff(f, vs[i])) / f)
deg = [psym * euler(br(0, 2), i) + euler(pt, i) for i in range(n5)]     # Grad von <13>^p mal Nenner je Bein
say(f"   Homogenität von <13>^p / (<12><23><34><45><51>) je Bein: {[str(d) for d in deg]}")
wb = sp.Symbol("w_b")
sol = sp.solve(sp.Eq(deg[0] + 2 + wb, 0), psym)
check("C5 Beine mit Gewicht 0 (Feld a): -2 + 2 (Maß) + 0 = 0; Beine mit Gewicht w_b: p - 2 + 2 + w_b = 0, also "
      "p = -w_b; für den Serre-Partner O(-4) aus z123 ist p = 4", deg[1] == -2 and deg[0] == psym - 2 and deg[2] == psym - 2
      and sol == [-wb] and sol[0].subs(wb, -4) == 4)
befund("C_vertex", {"linear_in_k": True, "einfachspur": True, "zaehler_p": 4},
       "Der Vertex der Twistorwirkung ist das zusammenhängende Funktional W = log Z der Naht, nicht Z. W ist linear "
       "im Niveau und reine Einfachspur; der Mehrspurteil aus z110 gehört zu Z und entfällt. Der Zähler <ij>^4 folgt "
       "aus dem Gewicht -4 des Partnerfelds b (Serre-Partner von O(0), z123), ohne Supersymmetrie: P_N4 entfällt. "
       "Bedingung ist, dass es b gibt; das ist das Feld der selbstdualen Hälfte des Ledgers (P_BF = Bestand des "
       "Vertrags) [E für die Algebra, C unter der Twistorwirkung als Rahmen].")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil D  Wirkungsleiter des Ledgers (Paper 1, Wirkungsgrammatik): Stufen und duale Lesart")
c3s, als = sp.Rational(1, 8) / sp.pi, sp.Symbol("alpha_inv", positive=True)
dtop = 48 * c3s ** 4
qfun = lambda x: dtop * sp.exp(-2 * x)                 # Nahtantwort q(alpha) = delta_top e^(-2 alpha) aus der alpha-Gleichung
rho = sp.Rational(3, 4) / sp.pi ** 2 * sp.exp(-2 * als)
check("D1 rho_Lambda/M^4 = (3/4 pi^2) e^(-2/alpha) = q(1/alpha)/c3^2 mit derselben Funktion q wie in der "
      "alpha-Gleichung [E, Identität der Ledgerformeln]", sp.simplify(rho - qfun(als) / c3s ** 2) == 0)
OmL = sp.Symbol("Omega_L", positive=True)
check("D2 Die Hubble-Stufe ist Friedmann: sqrt(rho_Lambda/(3 Omega_L M^4)) = e^(-1/alpha)/(2 pi sqrt(Omega_L)); "
      "die Leiter hat zwei unabhängige Stufen (elektroschwach, Lambda)",
      sp.simplify(sp.sqrt(rho / (3 * OmL)) - sp.exp(-als) / (2 * sp.pi * sp.sqrt(OmL))) == 0)
e_ = sp.Symbol("e", positive=True)
al_m = lambda n: (2 * sp.pi * n / e_) ** 2 / (4 * sp.pi)
check("D3 Duale Kopplung: mit e g = 4 pi (n = 2) ist alpha_m = 1/alpha, mit e g = 2 pi (n = 1) ist alpha_m = 1/(4 alpha)",
      sp.simplify(al_m(2) * (e_ ** 2 / (4 * sp.pi)) - 1) == 0
      and sp.simplify(al_m(1) * (e_ ** 2 / (4 * sp.pi)) - sp.Rational(1, 4)) == 0)
AINV = 137.0359992168407
A_L = 2 * AINV - math.log(3 / (4 * PI ** 2))
A_EW_ledger = 36.85                                    # -ln(v_geo/M), Paper 1, Tabelle der Wirkungsgrammatik
A_EW_nome = (2 * AINV - math.log(48 / (8 * PI) ** 4)) / 10    # falls die Stufe q(1/alpha)^(1/10) wäre
ZAHLEN.update(A_Lambda=round(A_L, 3), A_EW_ledger=A_EW_ledger, A_EW_aus_q=round(A_EW_nome, 3),
              log10_q_dual=round((math.log(48 / (8 * PI) ** 4) - 2 * AINV) / math.log(10), 3))
say(f"   -ln(rho_Lambda/M^4) = {A_L:.3f} (Ledger 276,65); log10 q(1/alpha) = {ZAHLEN['log10_q_dual']}")
say(f"   Gegenprobe elektroschwache Stufe: -ln(v/M) im Ledger {A_EW_ledger}; als q(1/alpha)^(1/10) wären es "
    f"{A_EW_nome:.3f}")
check("D4 Lambda-Stufe reproduziert (276,65)", abs(A_L - 276.65) < 0.01)
check("D5 Gegenprobe [N negativ]: Die elektroschwache Stufe ist keine Potenz derselben Antwort (Abstand > 8)",
      abs(A_EW_ledger - A_EW_nome) > 8)
befund("D_leiter", {"unabhaengige_stufen": 2, "lambda_als_duale_antwort": True, "ew_als_duale_antwort": False},
       "Die Lambda-Stufe ist die Nahtantwort q am dualen Punkt alpha -> 1/alpha (e g = 4 pi), samt Vorfaktor; in "
       "einer Twistorwirkung für Gravitation ist Lambda G genau der Koeffizient des nicht selbstdualen Vertex. Die "
       "elektroschwache Stufe folgt so nicht. Eine Identität, eine Lesart, eine Wahl (n = 2): Hinweis, kein "
       "Mechanismus [E für die Identität, L für die Lesart].")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil E  Bilanz und Urteil")
bilanz = [
    ("E", "Anomalie", "Sp_adj X^4 = lambda^2 (X.X)^2 gilt im Bestand nur für E8 (lambda^2 = 36 = h + 6), nicht für "
     "Spin(10), SU(5), Pati-Salam, Standardmodell."),
    ("E", "Statistik", "Die Uhr-Klassen heben den quartischen Teil nur gemeinsam auf (4 - 8 + 12 - 8); keine "
     "Uhr-Statistik und kein Spin(10)-Materieinhalt erfüllt die Bedingung."),
    ("E", "Vertex", "W = log Z der Naht ist linear in k und Einfachspur; der Zähler <ij>^4 folgt aus dem Gewicht -4."),
    ("C", "Zwei Hälften", "Unter der Twistorwirkung als Rahmen (P_TW) ersetzen sich P_BF und P_N4 durch den Bestand "
     "des Ledgers (selbstduale Hälfte), und der Geltungszaun 'nur selbstdual' des Vertrags wird durch den "
     "Nahtvertex geöffnet. Offen: Kopplungskoeffizient, Signatur, Gravitationsteil, Materie."),
    ("L", "Leiter", "Lambda-Stufe = duale Nahtantwort; elektroschwache Stufe nicht [N negativ]."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
urteil = ("Die beiden Programme passen an den geprüften Stellen zusammen. Das sind Verträglichkeitsprüfungen, kein "
          "Treffer aus Zahlen. (1) Die selbstduale Hälfte des Ledgers verlangt aus eigener Konsistenz E8 und schließt die "
          "kleineren Gruppen des Bestands aus; das ist ein vierdimensionaler Grund für E8, unabhängig von der "
          "Naht. (2) In dieser Hälfte müssen alle 248 Komponenten Bosonen sein. Die 16 kann dort nicht als Fermion "
          "sitzen; die chirale Materie muss von der Naht selbst kommen (nu(Uhr) = 16). Das ist eine harte Vorgabe "
          "für Tor T4 und kein Erfolg. (3) Der Nahtvertex ist das zusammenhängende Funktional log Z: Einfachspur, "
          "linear im Niveau, Zähler aus dem Gewicht. Damit entfallen P_N4 und der Mehrspureinwand aus z110, und "
          "P_BF wird zum Feld b der Ledger-Hälfte. (4) Nicht entschieden ist der Koeffizient, mit dem der Vertex "
          "in die Wirkung eingeht (die Kopplung), ebenso Signatur und Unitarität, der Gravitationsteil mit dem "
          "nativen Unendlichkeitstwistor und der Ort der Materie. (5) Die Wirkungsleiter hat nur zwei unabhängige "
          "Stufen. Die Lambda-Stufe hat die Form der Nahtantwort am dualen Punkt; die elektroschwache nicht. "
          "Das ist ein Kandidat für einen gemeinsamen vierdimensionalen Elternteil mit benannten nächsten "
          "Rechnungen, keine Herleitung. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z163_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z163.txt").write_text("\n".join(OUT) + "\n")
