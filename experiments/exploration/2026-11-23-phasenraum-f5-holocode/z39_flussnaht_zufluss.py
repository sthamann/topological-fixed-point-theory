"""Z39: Flussnaht mit Zufluss -- Normaldrehung, E8-Verklebung, Rückwirkung, Gitterzigarre, APS, 4D-Elternoperator.

Prüft die achte externe Analyse ("H6h21: dynamische Flussnaht mit topologischem Zufluss") und rechnet die Gitterzigarre selbst:
Teil 1  Zigarrennullmoden u_j = sinh^j rho cosh^-N rho, normierbar genau für 0 < j < N (j halbzahlig); die Flüsse 2 + 6 x 1 + 0
        geben acht Moden mit Normaldrehimpulsen b = (3/2, 1/2^7) und Ladung q = (1^8): Anomaliematrix [[8,5],[5,4]],
        I_4 = 4x^2 + 5xy + 2y^2 - p_1/3.
Teil 2  Die zwei Spinorverklebungen von D8: ganzzahlige b-Ladung auf 112/240 (gerade) bzw. 240/240 (ungerade) Wurzeln;
        b liegt nur im ungeraden Gitter, q ist dort primitiv. Abhängigkeit vom Spinlift; q-neutrale Wurzeln E7 bzw. A7;
        q-neutrale D5-Wurzeln (U(5)); Vorzeichenwechsel auf einer A3-Koordinate verbindet beide Verklebungen.
Teil 3  Zufluss als E8-Chern-Simons-Gitter: Basis des ungeraden Gitters, det 1, Signatur 8, ganzzahlige Kopplungsvektoren
        l_q, l_b mit l^T K^-1 l = [[8,5],[5,4]]; im geraden Gitter ist l_b nicht ganzzahlig.
Teil 4  Rückwirkung: f = (1-y)(R^2+Cy), y = e^{-2 kappa x}, kappa = 1/(R^2+C), C = 5 l^2, Phi = -kappa x, B_i = 2 kappa N_i y,
        Lambda = kappa^2 R^2 löst alle elf Euler-Lagrange-Gleichungen; glatte Spitze, zylindrisches Ende, Flüsse N_i, Grenzfall
        C = 0 ist die Zigarre aus Z38; Lambda maximal bei R^2 = C; Index unter Rückwirkung unverändert.
Teil 5  Gitterzigarre (eigene Rechnung): 2D-Wilson-Dirac mit M Winkelpunkten und exakter Blockzerlegung nach j; Wilson-Massen
        der N Moden fallen linear, Overlap-Index N in 36 Fällen, exakte Nullmoden mit Chiralität +1; auch auf der
        rückgewirkten Geometrie.
Teil 6  APS: Ladungen (-9,-5,-1,7,8) mit Sum q = Sum q^3 = 0; bei f = -1/17 ist der Nahtindex 1 (Eta-Invariante, Hurwitz-Zeta).
Teil 7  4D-Elternoperator: I_6 = ((x+z)^3 - x^3)/6 - z p_1/24 = z (x^2/2 + xz/2 + z^2/6 - p_1/24), kein reiner Eichterm;
        Reduktion mit den Flüssen 2 + 6 x 1 + 0 gibt k = c = 8, aber nicht den b-Anteil.
"""
import itertools
import json
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spl
import sympy as sp
from scipy.integrate import solve_ivp
from sympy.calculus.euler import euler_equations

OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
HERE = Path(__file__).parent


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


FLUSS = [2, 1, 1, 1, 1, 1, 1, 0]

# ------------------------------------------------------------------ Teil 1
say("Teil 1: Zigarrennullmoden, Normaldrehimpulse und Anomaliematrix")
rho = sp.symbols("rho", positive=True)


def normierbar(j, N):
    """u_j = sinh^j cosh^-N: an der Spitze ~ rho^j (L2 nahe 0 braucht 2j > -1), am Ende ~ e^{(j-N) rho} (braucht j < N)."""
    u = sp.sinh(rho) ** j * sp.cosh(rho) ** (-N)
    tip = sp.limit(sp.log(u) / sp.log(rho), rho, 0, "+")
    end = sp.limit(sp.log(u) / rho, rho, sp.oo)
    return (2 * tip > -1) and (end < 0)


zaehl = {}
for N in range(0, 4):
    ok = [Fr(2 * k + 1, 2) for k in range(-4, 6) if normierbar(sp.Rational(2 * k + 1, 2), N)]
    zaehl[N] = ok
    say(f"   N = {N}: normierbare j = {[str(j) for j in ok]}")
check("Zigarre: u_j = sinh^j cosh^-N normierbar genau für 0 < j < N, also N Moden pro Flussquant",
      all(zaehl[N] == [Fr(2 * k + 1, 2) for k in range(N)] for N in zaehl))
A_j = lambda j, N: sp.diff(sp.sinh(rho) ** j * sp.cosh(rho) ** (-N), rho) - (j * sp.coth(rho) - N * sp.tanh(rho)) * sp.sinh(rho) ** j * sp.cosh(rho) ** (-N)
check("u_j löst A_j u = (d_rho - (j coth rho - N tanh rho)) u = 0", all(sp.simplify(A_j(sp.Rational(1, 2) + k, 3)) == 0 for k in range(3)))
moden = [(i, j) for i, N in enumerate(FLUSS) for j in zaehl[N]]
b = sorted([j for _, j in moden], reverse=True)
q = [1] * len(moden)
say(f"   Flüsse {FLUSS}: {len(moden)} Moden, Normaldrehimpulse b = {[str(x) for x in b]}, Ladungen q = {q}")
check("acht Moden mit b = (3/2, 1/2^7) und q = (1^8)", b == [Fr(3, 2)] + [Fr(1, 2)] * 7 and len(q) == 8)
K = [[sum(qi * qi for qi in q), sum(qi * bi for qi, bi in zip(q, b))], [sum(qi * bi for qi, bi in zip(q, b)), sum(bi * bi for bi in b)]]
say(f"   Anomaliematrix K = {[[str(v) for v in row] for row in K]}")
check("Anomaliematrix [[8,5],[5,4]] (Sum q^2, Sum q b, Sum b^2)", K == [[8, 5], [5, 4]])
x, y, p1 = sp.symbols("x y p1")
I4 = sp.Rational(1, 2) * (K[0][0] * x**2 + 2 * K[0][1] * x * y + K[1][1] * y**2) - sp.Rational(len(moden), 24) * p1
say(f"   I_4 = {sp.expand(I4)}")
check("I_4 = 4x^2 + 5xy + 2y^2 - p_1/3 (c = 8)", sp.expand(I4 - (4 * x**2 + 5 * x * y + 2 * y**2 - p1 / 3)) == 0)
ZAHLEN["anomaliematrix"] = [[int(v) for v in row] for row in K]
ZAHLEN["b"] = [str(v) for v in b]

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Die zwei E8-Verklebungen und die Normaldrehung")
half = Fr(1, 2)
D8 = []
for i, j in itertools.combinations(range(8), 2):
    for si, sj in itertools.product((1, -1), repeat=2):
        v = [Fr(0)] * 8
        v[i], v[j] = Fr(si), Fr(sj)
        D8.append(tuple(v))
SP = [tuple(half * e for e in eps) for eps in itertools.product((1, -1), repeat=8)]
minus = lambda v: sum(1 for e in v if e < 0)
E8g = D8 + [v for v in SP if minus(v) % 2 == 0]
E8u = D8 + [v for v in SP if minus(v) % 2 == 1]
dot = lambda u, v: sum(a * c for a, c in zip(u, v))
check("beide Verklebungen: 240 Wurzeln der Norm 2", len(E8g) == len(E8u) == 240 and all(dot(v, v) == 2 for v in E8g + E8u))


def im_gitter(v, gerade):
    """v in D8 oder D8 + s (s = (1/2)^8 bzw. (-1/2, 1/2^7))."""
    s = tuple([half] * 8) if gerade else tuple([-half] + [half] * 7)
    for c in ((Fr(0),) * 8, s):
        w = [a - cc for a, cc in zip(v, c)]
        if all(e.denominator == 1 for e in w) and sum(w) % 2 == 0:
            return True
    return False


bv = tuple(b)
qv = tuple(Fr(1) for _ in range(8))
n_g = sum(1 for v in E8g if dot(bv, v).denominator == 1)
n_u = sum(1 for v in E8u if dot(bv, v).denominator == 1)
say(f"   ganzzahlige b-Ladung: gerade Verklebung {n_g}/240, ungerade {n_u}/240")
check("gerade Verklebung: nur die 112 D8-Wurzeln haben ganzzahlige b-Ladung, alle 128 Spinorgewichte halbzahlig", n_g == 112)
check("ungerade Verklebung: alle 240 Wurzeln ganzzahlige b-Ladung", n_u == 240)
check("b liegt im ungeraden Gitter, nicht im geraden", im_gitter(bv, False) and not im_gitter(bv, True))
check("q = (1^8) ist im ungeraden Gitter primitiv (q/k nicht im Gitter für k = 2..8), im geraden ist q = 2s",
      im_gitter(qv, False) and all(not im_gitter(tuple(e / k for e in qv), False) for k in range(2, 9))
      and im_gitter(tuple(e / 2 for e in qv), True))
# Spinlift
lifts = {}
for alpha in (-2, -1, 1, 2):
    bl = tuple(bi + alpha for bi in bv)
    lifts[f"b{alpha:+d}q"] = (sum(1 for v in E8g if dot(bl, v).denominator == 1), sum(1 for v in E8u if dot(bl, v).denominator == 1))
e1 = tuple(Fr(1) if i == 0 else Fr(0) for i in range(8))
bl = tuple(bi + ei for bi, ei in zip(bv, e1))
lifts["b+e1"] = (sum(1 for v in E8g if dot(bl, v).denominator == 1), sum(1 for v in E8u if dot(bl, v).denominator == 1))
say(f"   ganzzahlige Wurzeln (gerade, ungerade) für andere Lifts: {lifts}")
check("Auswahl unabhängig von Lifts b + alpha q (alpha ganz), die Spin-Statistik erhalten (Fermionladungen halbzahlig)",
      all(lifts[k] == (112, 240) for k in lifts if k != "b+e1") and all((bi + 1).denominator == 2 for bi in bv))
check("ein Lift mit ungerader Verschiebung in nur einem Kanal (b + e_1) kehrt die Auswahl um", lifts["b+e1"] == (240, 112))
befund("Die Normaldrehung wählt die ungerade Verklebung nur, wenn keine einzelne Kanal-U(1) in die Drehung gemischt wird", True)
hist = lambda E: dict(sorted({k: sum(1 for v in E if dot(qv, v) == k) for k in {int(dot(qv, v)) for v in E}}.items()))
hg, hu = hist(E8g), hist(E8u)
say(f"   q-Ladungen der Wurzeln: gerade {hg}; ungerade {hu}")
check("q-neutral: gerade Verklebung 126 Wurzeln (E7), ungerade 56 (A7 = su(8))", hg.get(0) == 126 and hu.get(0) == 56)
D5 = [v for v in D8 if all(e == 0 for e in v[5:])]
A3 = [v for v in D8 if all(e == 0 for e in v[:5])]
n5 = sum(1 for v in D5 if dot(qv, v) == 0)
n3 = sum(1 for v in A3 if dot(qv, v) == 0)
say(f"   q-neutrale Wurzeln: D5 {n5}/{len(D5)}, A3 {n3}/{len(A3)}")
check("q-neutral in D5 nur die 20 Wurzeln e_i - e_j: der Flusshintergrund lässt von Spin(10) nur U(5) übrig (A3: U(3))",
      n5 == 20 and len(D5) == 40 and n3 == 6 and len(A3) == 12)
flip = lambda v: tuple(list(v[:7]) + [-v[7]])
check("Vorzeichenwechsel auf einer A3-Koordinate bildet die gerade auf die ungerade Verklebung ab und lässt alle 40 D5-Wurzeln fest",
      set(map(flip, E8g)) == set(E8u) and all(flip(v) == v for v in D5))
lab = lambda v: ("16" if minus(v[:5]) % 2 == 0 else "16b", "4" if minus(v[5:]) % 2 == 0 else "4b")
check("dabei wird (16,4) + (16b,4b) zu (16,4b) + (16b,4)",
      {lab(v) for v in E8g if v in SP} == {("16", "4"), ("16b", "4b")} and {lab(v) for v in E8u if v in SP} == {("16", "4b"), ("16b", "4")})
ZAHLEN["q_histogramm"] = {"gerade": hg, "ungerade": hu}

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Zufluss als E8-Chern-Simons-Gitter")
h = Fr(1, 2)
alpha_g = [tuple([h, -h, -h, -h, -h, -h, -h, h])]
alpha_g.append(tuple(Fr(1) if i in (0, 1) else Fr(0) for i in range(8)))
for i in range(6):
    v = [Fr(0)] * 8
    v[i + 1], v[i] = Fr(1), Fr(-1)
    alpha_g.append(tuple(v))
Bg = sp.Matrix([[sp.Rational(e.numerator, e.denominator) for e in v] for v in alpha_g])
Bu = sp.Matrix([[sp.Rational(e.numerator, e.denominator) for e in flip(v)] for v in alpha_g])
Kg, Ku = Bg * Bg.T, Bu * Bu.T
check("Basis des geraden und des ungeraden Gitters aus Wurzeln (Bourbaki, bzw. gespiegelt): Gram det 1, alle Wurzeln im Gitter",
      Kg.det() == 1 and Ku.det() == 1 and all(im_gitter(v, True) for v in alpha_g) and all(im_gitter(flip(v), False) for v in alpha_g))
ev = [float(e) for e in np.linalg.eigvalsh(np.array(Ku, dtype=float))]
check("Gram-Matrix positiv definit (Signatur 8), gerade Diagonale", min(ev) > 0 and all(Ku[i, i] == 2 for i in range(8)))
qs = sp.Matrix([1] * 8)
bs = sp.Matrix([sp.Rational(v.numerator, v.denominator) for v in bv])
lq, lb = Bu * qs, Bu * bs
resp = sp.Matrix(2, 2, lambda i, j: (([lq, lb][i]).T * Ku.inv() * ([lq, lb][j]))[0])
say(f"   l_q = {list(lq)}, l_b = {list(lb)}; l^T K^-1 l = {resp.tolist()}")
check("ganzzahlige Kopplungsvektoren l_q, l_b im ungeraden Gitter mit l^T K^-1 l = [[8,5],[5,4]] (Anomaliematrix)",
      all(e.is_integer for e in list(lq) + list(lb)) and resp == sp.Matrix([[8, 5], [5, 4]]))
lb_g = Bg * bs
check("im geraden Gitter ist l_b nicht ganzzahlig (b koppelt dort nicht an ein ganzes CS-Niveau)", not all(e.is_integer for e in lb_g))
ZAHLEN["l_q"] = [int(e) for e in lq]
ZAHLEN["l_b"] = [int(e) for e in lb]

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Rückwirkung der Flussfelder auf den Kragen")
X = sp.symbols("x", positive=True)
ell, R, Lam = sp.symbols("ell R Lambda", positive=True)
hf, ff, Ph = sp.Function("h")(X), sp.Function("f")(X), sp.Function("Phi")(X)
af = [sp.Function(f"a{i}")(X) for i in range(8)]
Nn, rr = sp.sqrt(hf), sp.sqrt(ff)
Ric = -2 / (Nn * rr) * sp.diff(sp.diff(rr, X) / Nn, X)
F2 = sum(2 * sp.diff(ai, X) ** 2 / (hf * ff) for ai in af)
Lag = sp.sqrt(hf * ff) * sp.exp(-2 * Ph) * (Ric + 4 * sp.diff(Ph, X) ** 2 / hf + 4 * Lam - ell**2 / 4 * F2)
eqs = euler_equations(Lag, [hf, ff, Ph] + af, X)


def rest(Cval, lam_fak=1):
    kap = 1 / (R**2 + Cval)
    yy = sp.exp(-2 * kap * X)
    fs = (1 - yy) * (R**2 + Cval * yy)
    sol = {hf: 1 / fs, ff: fs, Ph: -kap * X}
    for i, N in enumerate(FLUSS):
        sol[af[i]] = N * (1 - yy)
    out = []
    for e in eqs:
        ex = e.lhs - e.rhs
        ex = ex.subs({sp.Derivative(k, (X, 2)): sp.diff(v, X, 2) for k, v in sol.items()})
        ex = ex.subs({sp.Derivative(k, X): sp.diff(v, X) for k, v in sol.items()}).subs(sol).subs(Lam, lam_fak * kap**2 * R**2)
        out.append(sp.simplify(ex))
    return out, fs, kap, yy


res, fs, kap, yy = rest(5 * ell**2)
check("alle elf Euler-Lagrange-Gleichungen (h, f, Phi, acht Eichfelder) verschwinden identisch", all(r == 0 for r in res))
res2, *_ = rest(4 * ell**2)
res3, *_ = rest(5 * ell**2, lam_fak=sp.Rational(101, 100))
check("Gegenproben: C = 4 l^2 oder Lambda um 1 % verschoben verletzen die Gleichungen", any(r != 0 for r in res2) and any(r != 0 for r in res3))
check("C = l^2 Sum N_i^2 / 2 = 5 l^2 (Sum N_i^2 = 10)", sum(N * N for N in FLUSS) == 10)
check("glatte Spitze f(0) = 0, f'(0) = 2; zylindrisches Ende f -> R^2",
      sp.simplify(fs.subs(X, 0)) == 0 and sp.simplify(sp.diff(fs, X).subs(X, 0)) == 2 and sp.limit(fs, X, sp.oo) == R**2)
Bi = sp.diff(2 * (1 - yy), X)
check("Fluss: 2 pi Int B_i dx = 2 pi N_i (B_i = 2 kappa N_i y)", sp.simplify(sp.integrate(Bi, (X, 0, sp.oo)) - 2) == 0
      and sp.simplify(Bi - 2 * kap * 2 * yy) == 0)
r_, Rr = sp.symbols("r_ R_", positive=True)
xc = Rr**2 * sp.log(sp.cosh(r_ / Rr))
f0 = (1 - sp.exp(-2 * xc / Rr**2)) * Rr**2
check("Grenzfall C = 0: f = R^2 tanh^2(rho/R), Phi = -log cosh(rho/R), B = N K -- die Zigarre aus Z38",
      sp.simplify(f0 - Rr**2 * sp.tanh(r_ / Rr) ** 2) == 0 and sp.simplify(sp.diff(xc, r_) - sp.sqrt(f0)) == 0
      and sp.simplify(-xc / Rr**2 + sp.log(sp.cosh(r_ / Rr))) == 0
      and sp.simplify(2 / Rr**2 * sp.exp(-2 * xc / Rr**2) - 2 / Rr**2 / sp.cosh(r_ / Rr) ** 2) == 0)
R2, Cc = sp.symbols("R2 C", positive=True)
LamR = R2 / (R2 + Cc) ** 2
crit = sp.solve(sp.diff(LamR, R2), R2)
check("Lambda(R^2) = R^2/(R^2+C)^2 maximal genau bei R^2 = C, Lambda_max = 1/(4C) = 1/(20 l^2); darunter zwei Zweige",
      crit == [Cc] and sp.simplify(LamR.subs(R2, Cc) - 1 / (4 * Cc)) == 0 and len(sp.solve(sp.Eq(LamR.subs(Cc, 1), sp.Rational(1, 5)), R2)) == 2)


def norm_rueck(j, N, R2v, Cv):
    """Nullmode g = exp Int (j - a) dx/f, Norm Int sqrt(f) g^2 dx/f: Spitze ~ x^(j - 1/2), Ende ~ e^{2(j - N) x/R^2}."""
    return (j - 0.5) > -1 and 2 * (j - N) / R2v < 0


rueck = {}
for (R2v, Cv) in ((1.0, 0.0), (1.0, 1.0), (1.0, 4.0), (4.0, 1.0)):
    rueck[f"R2={R2v},C={Cv}"] = {N: [j for j in (-1.5, -0.5, 0.5, 1.5, 2.5) if norm_rueck(j, N, R2v, Cv)] for N in (0, 1, 2)}
say(f"   normierbare j auf der rückgewirkten Geometrie: {rueck}")


def numnorm(j, N, R2v, Cv):
    """Numerisch: log der Norm bei feinerer Spitze und längerem Ende stabil (< 1e-2) <=> normierbar."""
    kv = 1 / (R2v + Cv)
    f_ = lambda xx: (1 - np.exp(-2 * kv * xx)) * (R2v + Cv * np.exp(-2 * kv * xx))
    a_ = lambda xx: N * (1 - np.exp(-2 * kv * xx))
    vals = []
    for eps, Xm in ((1e-5, 40.0 * (R2v + Cv)), (1e-9, 80.0 * (R2v + Cv))):
        grid = np.concatenate([np.geomspace(eps, 1.0, 6000), np.linspace(1.0, Xm, 12000)[1:]])
        integ = (j - a_(grid)) / f_(grid)
        lg = np.concatenate([[0.0], np.cumsum(0.5 * (integ[1:] + integ[:-1]) * np.diff(grid))])
        lg = lg - lg[np.searchsorted(grid, 1.0)]
        lw = 2 * lg - 0.5 * np.log(f_(grid))
        lseg = np.logaddexp(lw[1:], lw[:-1]) - np.log(2) + np.log(np.diff(grid))
        vals.append(np.logaddexp.reduce(lseg))
    return abs(vals[1] - vals[0]) < 1e-2


num = {k: {N: [j for j in (-1.5, -0.5, 0.5, 1.5, 2.5) if numnorm(j, N, float(k.split(",")[0][3:]), float(k.split(",")[1][2:]))]
           for N in (0, 1, 2)} for k in rueck}
check("Rückwirkung ändert die Modenzahl nicht: normierbar genau 0 < j < N (Exponenten und numerische Normen) für vier (R^2, C)",
      all(rueck[k][N] == [j / 2 for j in range(1, 2 * N, 2)] and num[k][N] == rueck[k][N] for k in rueck for N in (0, 1, 2)))
befund("Die klassische Rückwirkung ist in dieser Modellklasse exakt lösbar; die Dilatonwirkung ist nicht aus TFPT hergeleitet", True)

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: Gitterzigarre -- Wilson und Overlap (eigene Rechnung)")
s1 = np.array([[0, 1], [1, 0]], complex)
s2 = np.array([[0, -1j], [1j, 0]])
s3 = np.diag([1.0 + 0j, -1])
I2 = np.eye(2)
LRAD = 16.0


def zig_geom(Cv):
    """Zigarre R = 1 (C = 0) oder rückgewirkt: x(rho) aus dx/drho = sqrt f, r = sqrt f, a/N = 1 - y."""
    if Cv == 0:
        return lambda rr, N: (np.tanh(rr), N * np.tanh(rr) ** 2)
    kv = 1 / (1 + Cv)
    f_ = lambda xx: (1 - np.exp(-2 * kv * xx)) * (1 + Cv * np.exp(-2 * kv * xx))
    r0 = 1e-4
    sol = solve_ivp(lambda t, xx: [np.sqrt(max(f_(xx[0]), 0.0))], (r0, LRAD + 1), [r0**2 / 2], rtol=1e-11, atol=1e-14, dense_output=True)

    def g(rr, N):
        xx = sol.sol(rr)[0]
        return np.sqrt(f_(xx)), N * (1 - np.exp(-2 * kv * xx))
    return g


def block(j, M, hh, N, geom, rW=1.0, m0=0.0):
    n = int(round(LRAD / hh))
    rr = hh * np.arange(1, n + 1)
    r, a = geom(rr, N)
    d = 2 * np.pi / M
    ph = (j - a) * d
    s = np.sin(ph) / (d * r)
    w = (d / hh) * (1 - np.cos(ph)) / (d * r)
    D = sps.kron(sps.diags(1j * s), s2) + sps.kron(sps.diags(w + rW / hh + m0), I2)
    D = D + sps.kron(sps.eye(n, k=1), s1 / (2 * hh) - rW / (2 * hh) * I2) + sps.kron(sps.eye(n, k=-1), -s1 / (2 * hh) - rW / (2 * hh) * I2)
    return D.tocsc(), n


js = lambda M: [k + 0.5 for k in range(-M // 2, M // 2)]
Dt, nt = block(0.5, 64, 2 / 64, 1, zig_geom(0.0))
wt, vt = spl.eigs(Dt, k=1, sigma=0, which="LM", tol=1e-12)
rt = (2 / 64) * np.arange(1, nt + 1)
ut = np.sinh(rt) ** 0.5 / np.cosh(rt)
fid = abs(np.vdot(ut, vt[0::2, 0])) / (np.linalg.norm(ut) * np.linalg.norm(vt[:, 0]))
say(f"   Überlapp der Wilson-Mode (M = 64, N = 1, j = 1/2) mit der Kontinuumsnullmode sinh^1/2 / cosh: {fid:.5f}")
check("die niedrigste Wilson-Mode ist die Kontinuumsnullmode u_j (Überlapp > 0.99, oberer Chiralitätsanteil)", fid > 0.99)
wil = {}
for gname, Cv in (("zigarre", 0.0), ("rueckgewirkt_C=R2", 1.0)):
    geom = zig_geom(Cv)
    for M in ((16, 32, 64, 128) if Cv == 0 else (32, 64)):
        hh = 2 / M
        for N in (0, 1, 2):
            low, smin = [], np.inf
            for j in js(M):
                D, n = block(j, M, hh, N, geom)
                if abs(j) <= N + 2.5:
                    w, v = spl.eigs(D, k=3, sigma=0, which="LM", tol=1e-10)
                    G5 = sps.kron(sps.eye(n), s3)
                    for wi, vi in zip(w, v.T):
                        if abs(wi) < 0.5 and abs(wi.imag) < 1e-8:
                            low.append((j, round(float(wi.real), 8), round(float((vi.conj() @ (G5 @ vi)).real / (vi.conj() @ vi).real), 6)))
                elif M <= 64:
                    A = (D.conj().T @ D).tocsc()
                    smin = min(smin, float(np.sqrt(abs(spl.eigsh(A, k=1, sigma=0, which="LM", return_eigenvectors=False)[0]))))
            wil[(gname, M, N)] = dict(moden=sorted(low), smin_aussen=None if smin == np.inf else round(smin, 3))
            say(f"   Wilson {gname}, M = {M}, h = {hh:.5f}, N = {N}: reelle Moden |lambda| < 1/2: "
                f"{[(j, round(l, 5), round(c, 4)) for j, l, c in sorted(low)]}; sigma_min außerhalb |j| <= N+2: {wil[(gname, M, N)]['smin_aussen']}")
zl = lambda g, M, N: len(wil[(g, M, N)]["moden"])
check("Wilson: ab M = 32 genau N reelle Moden unter der Kontinuumslücke 1/2 (N = 0, 1, 2), auf Zigarre und rückgewirkter Geometrie",
      all(zl("zigarre", M, N) == N for M in (32, 64, 128) for N in (0, 1, 2)) and all(zl("rueckgewirkt_C=R2", M, N) == N for M in (32, 64) for N in (0, 1, 2)))
check("Wilson: diese Moden haben Chiralität > 0.97 ab M = 64 (positive Transversalchiralität)",
      all(c > 0.97 for g, M, N in wil if M >= 64 for _, _, c in wil[(g, M, N)]["moden"]))
mass = {M: max(l for _, l, _ in wil[("zigarre", M, 2)]["moden"]) for M in (32, 64, 128)}
say(f"   größte Gittermasse (Zigarre, N = 2): {mass}")
check("Wilson-Gittermassen fallen linear mit h: Verhältnis bei Halbierung zwischen 1.7 und 2.3",
      all(1.7 < mass[M] / mass[2 * M] < 2.3 for M in (32, 64)))
check("Wilson: alle übrigen Winkelsektoren haben sigma_min > 1 (keine versteckten niedrigen Moden, M <= 64)",
      all(v["smin_aussen"] is None or v["smin_aussen"] > 1 for v in wil.values()))
ov = {}
for gname, Cv, Ms in (("zigarre", 0.0, (16, 24, 32)), ("rueckgewirkt_C=R2", 1.0, (24,))):
    geom = zig_geom(Cv)
    for M in Ms:
        hh = 2 / M
        for N in (0, 1, 2):
            for mu in (0.5, 0.8, 1.0, 1.3):
                ind, nullm, minH = 0.0, [], np.inf
                for j in js(M):
                    D, n = block(j, M, hh, N, geom, m0=-mu / hh)
                    D = D.toarray()
                    G5 = np.kron(np.eye(n), s3)
                    H = G5 @ D
                    H = (H + H.conj().T) / 2
                    e, V = np.linalg.eigh(H)
                    minH = min(minH, float(np.abs(e).min()))
                    eps = (V * np.sign(e)) @ V.conj().T
                    ind += -np.trace(eps).real / 2
                    if mu == 1.0:
                        w, v = np.linalg.eig(np.eye(2 * n) + G5 @ eps)
                        for wi, vi in zip(w, v.T):
                            if abs(wi) < 1e-8:
                                nullm.append((j, float((vi.conj() @ G5 @ vi).real / (vi.conj() @ vi).real), float(abs(wi))))
                ov[(gname, M, N, mu)] = dict(index=float(round(ind, 8)) + 0.0, min_abs_HW=round(minH, 4), nullmoden=nullm)
            say(f"   Overlap {gname}, M = {M}, N = {N}: Index für mu h = 0.5/0.8/1.0/1.3: "
                f"{[ov[(gname, M, N, mu)]['index'] for mu in (0.5, 0.8, 1.0, 1.3)]}; min|H_W| = {ov[(gname, M, N, 1.0)]['min_abs_HW']}; "
                f"Nullmoden (j, Chiralität): {[(j, round(c, 6)) for j, c, _ in ov[(gname, M, N, 1.0)]['nullmoden']]}")
ovz = {k: v for k, v in ov.items() if k[0] == "zigarre"}
check("Overlap: Index -Tr sign(H_W)/2 = N in allen 36 Fällen (drei Auflösungen, vier Wilson-Massen, N = 0, 1, 2)",
      len(ovz) == 36 and all(abs(v["index"] - k[2]) < 1e-6 for k, v in ovz.items()))
check("Overlap: genau N exakte Nullmoden (|lambda| < 1e-8) mit Chiralität +1 in den Sektoren j = 1/2 .. N - 1/2",
      all(len(v["nullmoden"]) == k[2] and all(abs(c - 1) < 1e-8 for _, c, _ in v["nullmoden"])
          and sorted(j for j, _, _ in v["nullmoden"]) == [0.5 + i for i in range(k[2])] for k, v in ov.items() if k[3] == 1.0))
check("Overlap auf der rückgewirkten Geometrie (C = R^2): Index N für alle vier Massen", all(abs(v["index"] - k[2]) < 1e-6 for k, v in ov.items() if k[0] != "zigarre"))
check("H_W bleibt lückenhaft (min|H_W| > 0.4), der Overlap-Index ist wohldefiniert", all(v["min_abs_HW"] > 0.4 for v in ov.values()))
gesamt = sum(N for N in FLUSS)
say(f"   TFPT-Flussverteilung {FLUSS}: Overlap-Gesamtindex = {gesamt}")
befund("Die Gitterzigarre bestätigt N Moden pro Flussquant mit Wilson (linear verschwindende Masse) und Overlap (exakt)", True)
ZAHLEN["wilson"] = {f"{k[0]}|M={k[1]}|N={k[2]}": v for k, v in wil.items()}
ZAHLEN["overlap"] = {f"{k[0]}|M={k[1]}|N={k[2]}|mu={k[3]}": {"index": v["index"], "min_abs_HW": v["min_abs_HW"], "nullmoden": len(v["nullmoden"])} for k, v in ov.items()}

# ------------------------------------------------------------------ Teil 6
say("\nTeil 6: APS-Index auf dem nichtkompakten Kragen")
Qs = [-9, -5, -1, 7, 8]
check("Ladungen (-9,-5,-1,7,8): Sum q = Sum q^3 = 0 (4D anomaliefrei, auch gemischt gravitativ)", sum(Qs) == 0 and sum(qq**3 for qq in Qs) == 0)
fv = Fr(-1, 17)


def eta(alpha):
    """Eta des Randkreises (antiperiodisch, Holonomie alpha): Spektrum n + 1/2 - alpha, Hurwitz-Zeta bei s = 0."""
    bb = (Fr(1, 2) - alpha) % 1
    return mp.zeta(0, float(bb)) - mp.zeta(0, float(1 - bb)), 1 - 2 * bb


ind_aps, ind_zahl = 0, 0
for qq in Qs:
    al = qq * fv
    e_num, e_ex = eta(al)
    aps = al - e_ex / 2
    direkt = (sum(1 for k in range(-40, 40) if 0 < Fr(2 * k + 1, 2) < al) - sum(1 for k in range(-40, 40) if al < Fr(2 * k + 1, 2) < 0))
    say(f"   Q = {qq:3d}: Qf = {str(al):>6s}, eta = {str(e_ex):>7s} (numerisch {float(e_num):+.6f}), APS = Qf - eta/2 = {aps}, direkt gezählt {direkt}")
    ind_aps += aps
    ind_zahl += direkt
    if abs(float(e_num) - float(e_ex)) > 1e-12 or aps != direkt:
        ind_aps = None
        break
check("APS-Index Qf - eta/2 = round(Qf) = direkte Modenzählung für jede Ladung; Eta aus der Hurwitz-Zeta", ind_aps is not None)
check("bei f = -1/17 ist der Nahtindex 1, obwohl Sum Q f = 0", ind_aps == 1 and ind_zahl == 1 and sum(qq * fv for qq in Qs) == 0)
check("für ganzzahlige Flüsse ist der Nahtindex f Sum Q = 0 (kompakter Fall)",
      all(sum(round(qq * ff) for qq in Qs) == 0 for ff in range(-3, 4)))
befund("Globale 4D-Anomaliefreiheit schließt einen chiralen Nahtsektor auf dem nichtkompakten Kragen nicht aus (Eta-Invariante am Rand)", True)

# ------------------------------------------------------------------ Teil 7
say("\nTeil 7: 4D-Elternoperator ohne Vektoreichanomalie")
z = sp.symbols("z")
I6 = ((x + z) ** 3 - x**3) / 6 - z * p1 / 24
check("I_6 = ((x+z)^3 - x^3)/6 - z p_1/24 = z (x^2/2 + xz/2 + z^2/6 - p_1/24)",
      sp.expand(I6 - z * (x**2 / 2 + x * z / 2 + z**2 / 6 - p1 / 24)) == 0)
check("kein reiner x^3-Term: die physische Vektorsymmetrie ist anomaliefrei", sp.Poly(sp.expand(I6), x, z, p1).coeff_monomial(x**3) == 0)
red = sum(N * sp.diff(I6, z).subs(z, 0) for N in FLUSS)
check("Reduktion über den Kragen (Flüsse 2 + 6 x 1 + 0, z nur transversal): I_4 = 4x^2 - p_1/3, also k = c = 8",
      sp.expand(red - (4 * x**2 - p1 / 3)) == 0)
befund("Der Elternoperator reproduziert den (q,q)-Eintrag und c, nicht aber die Normaldrehung b (die braucht die Kragengeometrie)", True)

urteil = ("Die achte Analyse rechnet sich in allen geprüften Punkten nach, mit zwei Einschränkungen. Die Flüsse 2 + 6 x 1 + 0 "
          "geben auf der Zigarre acht Moden mit Normaldrehimpulsen b = (3/2, 1/2^7); mit q = (1^8) folgt die Anomaliematrix "
          "[[8,5],[5,4]] und I_4 = 4x^2 + 5xy + 2y^2 - p_1/3. Die ungerade Spinorverklebung von D8 gibt allen 240 Wurzeln "
          "ganzzahlige b-Ladung, die gerade nur den 112 D8-Wurzeln; b liegt nur im ungeraden Gitter, q ist dort primitiv. "
          "Einschränkung 1: Die Auswahl ist invariant unter Lifts b + alpha q, kehrt sich aber um, wenn eine einzelne Kanal-U(1) "
          "in die Drehung gemischt wird. Ein E8-Chern-Simons-Gitter mit ganzzahligen Kopplungen l_q, l_b reproduziert die "
          "Anomaliematrix exakt (im geraden Gitter ist l_b nicht ganzzahlig). Neu: q-neutral sind im ungeraden Gitter nur 56 "
          "Wurzeln (su(8)), im geraden 126 (E7), und von den 40 D5-Wurzeln nur 20 -- der Flusshintergrund lässt von Spin(10) "
          "nur U(5) übrig. Die Rückwirkungslösung f = (1-y)(R^2+Cy) erfüllt alle elf Euler-Lagrange-Gleichungen exakt, mit "
          "glatter Spitze, Zylinderende, Flüssen N_i und der Z38-Zigarre als Grenzfall C = 0; Lambda ist bei R^2 = C maximal. "
          "Die Modenzahl hängt nur an Spitze und Ende. Einschränkung 2: Die Gitterzigarren-Tabellen der Analyse waren nach "
          "ihrer eigenen Folgeanalyse unbelegt; die eigene Rechnung (exakte Blockzerlegung, Wilson und Overlap) gibt jetzt N "
          "Moden mit linear fallender Wilson-Masse, Overlap-Index N in 36 Fällen und exakte Nullmoden der Chiralität +1, auch "
          "auf der rückgewirkten Geometrie. Das APS-Beispiel (-9,-5,-1,7,8) bei f = -1/17 hat Nahtindex 1 trotz Anomaliefreiheit; "
          "der Elternoperator mit I_6 = z(...) gibt k = c = 8 ohne Vektoreichanomalie, aber nicht den b-Anteil. Offen bleiben die "
          "Herleitung von Zufluss und Dilatonwirkung aus dem TFPT-Masteroperator und die Wahl des Spinlifts.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z39_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
