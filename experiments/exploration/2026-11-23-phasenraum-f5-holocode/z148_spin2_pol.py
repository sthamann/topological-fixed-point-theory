"""Z148: Hat der native Gravitationssektor von H6 einen masselosen Spin-2-Pol, und welche Prämisse braucht er genau?

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart, N numerisch, O offen.

Frage: Hat der native Gravitationssektor einen masselosen Spin-2-Pol mit zwei Helizitäten und positivem Residuum, und
welche Prämisse braucht er in kleinstmöglicher Form?

Bestand: z96/z98 (MacDowell-Mansouri: b = l^2/(64 pi G), Lambda = -3/l^2, l = R_Vak/2), z118/z127 (Uhrspur gibt
MacDowell-Mansouri, R gibt 64 eps je Index, Lambda G = -3/(1024 pi kappa), G/R_Vak^2 = 1/(4096 pi kappa); Vorzeichen =
Paarungsbit), z137 (nativer Spin-2-Primär bei Delta = 4, C_T = 1/(3 pi^4), c = 1/120 je Leiter; <TT> ~ c p^4 ln p^2;
alpha_W = 3/(64 pi G Lambda)), z124 (nativer Doppelspuranteil bei n = 4 = konforme Gravitation, s^2 A = 4 T.T, Stärke
k = 1 relativ zu Yang-Mills, unter Lesart T und P_N4), z110 (Einspuranteil mit Koeffizient -k), z115/z128 (jede native
Familie ist Maurer-Cartan-flach), z75 (kein Graviton aus der Virasoro-Antwort; Helizität +-2 nur im Dubleton-Sektor),
z84/z123 (tau = Punktspiegelung, J = theta o tau, theta = radiale Inversion), z135 (drei Lesarten der Eichkopplung).

Konventionen: Koordinaten (t, x, y, z), eta = diag(-1, 1, 1, 1), AdS4 im Poincaré-Bild e^a = (l/z) dx^a, Rand bei z = 0.
alpha_W ist der euklidische Koeffizient: S_E = alpha_W Int sqrt(g) W^2, lorentzsch S_L = -alpha_W Int sqrt(-g) W^2
(wie z137). Felder hängen von (t, z) ab: zeitartiger Randimpuls im Ruhesystem, Helizität = Ladung unter der Drehung
in der (x, y)-Ebene.

Teil A  Linearisierte MacDowell-Mansouri-Wirkung um AdS4
  A1  Hintergrund F = 0; quadratische Wirkung L2 = b eps_abcd f^ab ^ f^cd; Operator erster Ordnung, 40 Felder.
  A2  Eichfreiheit: zehn Parameter (so(3,2)); f verschwindet auf reiner Eichung; Hesse-Form der nativen Familie.
  A3  Spinzusammenhang ist Hilfsfeld (24 x 24, Rang 24); nach Elimination: linearisierte Einstein-Gleichung.
  A4  Freiheitsgrade: Jetzählung (Zuwachs je Grad), Helizitätszerlegung.
  A5  Helizität +-2: Gleichung, Exponenten Delta = 0, 3, Vorzeichen des kinetischen Terms, Residuum.
  A6  Ohne P_dyn: so(3,2)-Inhalt des nativen Fockraums; Unitaritätsschranke; D(3,2) kommt nicht vor.
Teil B  Induzierter Propagator aus <TT> = c p^4 ln p^2
  B1  Koeffizient: Fourier-Transformierte von C_T I/x^8 im Tensorkanal.
  B2  (a) ohne Baumterm: kein masseloser Pol, Tachyonpol an der Logarithmusskala, kein Newton-Gesetz.
  B3  (b) mit MacDowell-Mansouri-Baumterm: Pole über Lambert-W, Residuen, Summenregel, Unitarität.
  B4  Die Skala des Logarithmus.
Teil C  Konforme Gravitation (z124) mit I_AdS und Randbedingung
  C1  Linearisierter Weyl-Operator im Helizität-2-Sektor; Faktorisierung masselos mal teilweise masselos.
  C2  Moden, Neumann-Bedingung, Regularität; Wirkungen auf der Einstein-Mode; alpha_W = -l^2/(64 pi G).
  C3  Native Kopplung: Yang-Mills-Baumamplitude gegen Weyl-Austausch, exakt rational; alpha_W = -1/(4 g^2 k).
  C4  Die Zahl G/R_Vak^2 und der Abgleich mit dem Uhrspurgewicht.
  C5  Erzwingt ein natives Datum P_Neumann? Randspiegelung, Paritäten, Unitaritätsschranke.
Teil D  Bilanz und Urteil.
"""
import itertools
import json
import os
import random
import time
from fractions import Fraction as Fr
from math import comb, gcd
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # deterministische Rundung (einfädig)

import mpmath as mp
import sympy as sp
from flint import fmpz_mat

HERE = Path(__file__).parent
REPO = HERE.parents[2]
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
mp.mp.dps = 30


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val, text):
    BEFUNDE[name] = val
    say(f"   [Befund: {val}] {text}")


def gs(x, n=6):
    """Gerundete Zahl als Text (feste Zahl gültiger Stellen)."""
    return f"{float(x):.{n}g}"


# ====================================================================================================================
# Werkzeug: lineare Jetformen. Schlüssel (J, p, q) steht für d_t^p d_z^q Phi_J; Koeffizienten sind Funktionen von z, l.
# ====================================================================================================================
z, l = sp.symbols("z l", positive=True)
ETA = [-1, 1, 1, 1]
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
PIDX = {p: i for i, p in enumerate(PAIRS)}
EINS = sp.Integer(1)


def eps4(a, b, c, d):
    p = [a, b, c, d]
    if len(set(p)) < 4:
        return 0
    s = 1
    for i in range(4):
        for j in range(i + 1, 4):
            if p[i] > p[j]:
                s = -s
    return s


def lsum(terms):
    out = {}
    for A, c in terms:
        if c == 0:
            continue
        for k, v in A.items():
            out[k] = out.get(k, 0) + c*v
    out = {k: sp.expand(v) for k, v in out.items()}
    return {k: v for k, v in out.items() if v != 0}


def lin_scale(A, c):
    return lsum([(A, c)])


def Dt(A):
    return {(J, p + 1, q): c for (J, p, q), c in A.items()}


def Dz(A):
    out = {}
    for (J, p, q), c in A.items():
        out[(J, p, q + 1)] = out.get((J, p, q + 1), 0) + c
        dc = sp.diff(c, z)
        if dc != 0:
            out[(J, p, q)] = out.get((J, p, q), 0) + dc
    out = {k: sp.expand(v) for k, v in out.items()}
    return {k: v for k, v in out.items() if v != 0}


def dcoord(A, mu):
    return Dt(A) if mu == 0 else (Dz(A) if mu == 3 else {})


def subst(A, mapping):
    """Ersetzt das Feld J durch die Jetform mapping[J], Ableitungen eingeschlossen."""
    terms = []
    for (J, p, q), c in A.items():
        if J in mapping:
            X = mapping[J]
            for _ in range(p):
                X = Dt(X)
            for _ in range(q):
                X = Dz(X)
            terms.append((X, c))
        else:
            terms.append(({(J, p, q): EINS}, c))
    return lsum(terms)


def quad_add(Q, A, B, c):
    for k1, c1 in A.items():
        for k2, c2 in B.items():
            kk = (k1, k2) if k1 <= k2 else (k2, k1)
            Q[kk] = Q.get(kk, 0) + c*c1*c2


def quad_clean(Q):
    Q = {k: sp.expand(v) for k, v in Q.items()}
    return {k: v for k, v in Q.items() if v != 0}


def euler(Q, J):
    """Euler-Lagrange-Ausdruck dS/dPhi_J einer quadratischen Jetform (beliebige Ableitungsordnung)."""
    part = {}
    for (k1, k2), c in Q.items():
        if k1[0] == J:
            part.setdefault(k1, []).append(({k2: EINS}, c))
        if k2[0] == J:
            part.setdefault(k2, []).append(({k1: EINS}, c))
    terms = []
    for key in sorted(part):
        X = lsum(part[key])
        for _ in range(key[1]):
            X = Dt(X)
        for _ in range(key[2]):
            X = Dz(X)
        terms.append((X, (-1)**(key[1] + key[2])))
    return lsum(terms)


def ordnung(Es):
    return max((p + q) for e in Es for (J, p, q) in e)


# ------------------------------------------------------------------ Hintergrund AdS4
def ebar(a, mu):
    return l/z if a == mu else 0


def wbar(a, c, mu):
    """Spinzusammenhang omega^{ac}_mu des Hintergrunds."""
    if c == 3 and a != 3 and mu == a:
        return -1/z
    if a == 3 and c != 3 and mu == c:
        return 1/z
    return 0


def wmix(a, c, mu):
    return wbar(a, c, mu)*ETA[c]                       # omega^a_c


def anti(tab, a, c):
    """Zugriff auf eine antisymmetrische Tabelle über Paare a < c."""
    if a == c:
        return None, 0
    return (tab[PIDX[(a, c)]], 1) if a < c else (tab[PIDX[(c, a)]], -1)


def build_f(h, tau):
    """f^{ab}_{mu nu} = (D tau)^{ab} + (h^a ^ e^b - h^b ^ e^a)/l^2 für a < b, mu < nu."""
    f = {}
    for (a, b) in PAIRS:
        for mu in range(4):
            for nu in range(mu + 1, 4):
                T = tau[PIDX[(a, b)]]
                terms = [(dcoord(T[nu], mu), 1), (dcoord(T[mu], nu), -1)]
                for c in range(4):
                    t1, s1 = anti(tau, c, b)
                    if s1:
                        terms += [(t1[nu], s1*wmix(a, c, mu)), (t1[mu], -s1*wmix(a, c, nu))]
                    t2, s2 = anti(tau, a, c)
                    if s2:
                        terms += [(t2[nu], s2*wmix(b, c, mu)), (t2[mu], -s2*wmix(b, c, nu))]
                terms += [(h[a][mu], ebar(b, nu)/l**2), (h[a][nu], -ebar(b, mu)/l**2),
                          (h[b][mu], -ebar(a, nu)/l**2), (h[b][nu], ebar(a, mu)/l**2)]
                f[(a, b, mu, nu)] = lsum(terms)
    return f


def hJ(a, mu):
    return a*4 + mu


def tJ(pi, mu):
    return 16 + pi*4 + mu


SPLITS = [((0, 1), (2, 3), 1), ((0, 2), (1, 3), -1), ((0, 3), (1, 2), 1), ((1, 2), (0, 3), 1), ((1, 3), (0, 2), -1),
          ((2, 3), (0, 1), 1)]

# ------------------------------------------------------------------ Jetzählung (formale Potenzreihen bei z = 1, l = 1)
_DC = {}


def dcoef(c, n):
    key = (c, n)
    if key not in _DC:
        v = sp.diff(c, z, n) if n else c
        v = sp.Rational(v.subs({z: 1, l: 1}))
        _DC[key] = Fr(int(v.p), int(v.q))
    return _DC[key]


def jet_index(fields, N):
    idx = {}
    for J in fields:
        for d in range(N + 1):
            for p in range(d + 1):
                idx[(J, p, d - p)] = len(idx)
    return idx


def jet_rows(exprs, idx, N, ordmax):
    """Taylorkoeffizienten d_t^m d_z^n (Gleichung) im Punkt, als Zeilen in den Jetkoordinaten."""
    rows = []
    for e in exprs:
        if not e:
            continue
        for d in range(N - ordmax + 1):
            for m in range(d + 1):
                n = d - m
                row = {}
                for (J, p, q), c in e.items():
                    for j in range(n + 1):
                        v = dcoef(c, n - j)
                        if v:
                            col = idx[(J, p + m, q + j)]
                            row[col] = row.get(col, 0) + comb(n, j)*v
                row = {k: v for k, v in row.items() if v}
                if row:
                    rows.append(row)
    return rows


def rank_rows(rows, ncols, cols=None):
    if cols is not None:
        cmap = {c: i for i, c in enumerate(cols)}
        rows = [{cmap[k]: v for k, v in r.items() if k in cmap} for r in rows]
        rows = [r for r in rows if r]
        ncols = len(cols)
    if not rows or ncols == 0:
        return 0
    ent = []
    for r in rows:
        den = 1
        for v in r.values():
            den = den*v.denominator//gcd(den, v.denominator)
        line = [0]*ncols
        for k, v in r.items():
            line[k] = int(v*den)
        ent.append(line)
    return fmpz_mat(ent).rank()


def zaehlung(Eqs, fields, gauge=None, gparams=None, N=6, mmax=4):
    """Für m = 0..mmax: Dimension der m-Jets von Lösungen (fortsetzbar bis Grad N), der Eichjets und die Differenz."""
    idx = jet_index(fields, N)
    rows = jet_rows(Eqs, idx, N, ordnung(Eqs))
    rA = rank_rows(rows, len(idx))
    out = []
    for m in range(mmax + 1):
        high = sorted(i for k, i in idx.items() if k[1] + k[2] > m)
        sol = (len(idx) - len(high)) - rA + rank_rows(rows, len(idx), high)
        g = 0
        if gauge is not None:
            gidx = jet_index(gparams, m + 1)
            grow = []
            for J in fields:
                e = gauge[J]
                for d in range(m + 1):
                    for p in range(d + 1):
                        q = d - p
                        row = {}
                        for (K, a, b), c in e.items():
                            for j in range(q + 1):
                                v = dcoef(c, q - j)
                                if v:
                                    col = gidx[(K, a + p, b + j)]
                                    row[col] = row.get(col, 0) + comb(q, j)*v
                        row = {k: v for k, v in row.items() if v}
                        if row:
                            grow.append(row)
            g = rank_rows(grow, len(gidx))
        out.append((m, sol, g, sol - g))
    return out


# ====================================================================================================================
say("Teil A: Linearisierte MacDowell-Mansouri-Wirkung um AdS4")
# A1 ---------------------------------------------------------------------------------------------------------------
tors_ok, flach_ok = True, True
for a in range(4):
    for mu in range(4):
        for nu in range(mu + 1, 4):
            de = (sp.diff(ebar(a, nu), z) if mu == 3 else 0) - (sp.diff(ebar(a, mu), z) if nu == 3 else 0)
            tors = de + sum(wmix(a, c, mu)*ebar(c, nu) - wmix(a, c, nu)*ebar(c, mu) for c in range(4))
            tors_ok &= sp.simplify(tors) == 0
for (a, b) in PAIRS:
    for mu in range(4):
        for nu in range(mu + 1, 4):
            dw = (sp.diff(wbar(a, b, nu), z) if mu == 3 else 0) - (sp.diff(wbar(a, b, mu), z) if nu == 3 else 0)
            ww = sum(wmix(a, c, mu)*wbar(c, b, nu) - wmix(a, c, nu)*wbar(c, b, mu) for c in range(4))
            ee = (ebar(a, mu)*ebar(b, nu) - ebar(a, nu)*ebar(b, mu))/l**2
            flach_ok &= sp.simplify(dw + ww + ee) == 0
check("[E] Hintergrund: e^a = (l/z) dx^a mit omega^(a3) = -dx^a/z ist torsionsfrei, und die so(3,2)-Krümmung "
      "F^ab = R^ab + e^a ^ e^b/l^2 verschwindet identisch (alle 36 Komponenten). AdS4 ist der flache "
      "so(3,2)-Zusammenhang; die MacDowell-Mansouri-Dichte eps F ^ F und ihre erste Variation sind dort null",
      tors_ok and flach_ok)

H = [[{(hJ(a, mu), 0, 0): EINS} for mu in range(4)] for a in range(4)]
TAU = [[{(tJ(pi, mu), 0, 0): EINS} for mu in range(4)] for pi in range(6)]
FLIN = build_f(H, TAU)
Q2 = {}
for (a, b) in PAIRS:
    for (c, d) in PAIRS:
        e4 = eps4(a, b, c, d)
        if e4:
            for m1, m2, sg in SPLITS:
                quad_add(Q2, FLIN[(a, b) + m1], FLIN[(c, d) + m2], 4*e4*sg)
Q2 = quad_clean(Q2)                                   # L2/b, Koeffizient von dt ^ dx ^ dy ^ dz
EMM = [euler(Q2, J) for J in range(40)]
ord_mm = ordnung(EMM)
tau_abl = any(J >= 16 and p + q > 0 for e in EMM[16:] for (J, p, q) in e)
h_abl = any(J < 16 and p + q > 0 for e in EMM[:16] for (J, p, q) in e)
check(f"[E] Quadratische Wirkung: Weil F im Hintergrund verschwindet, ist L2 = b eps_abcd f^ab ^ f^cd mit "
      f"f = D tau + (h ^ e - e ^ h)/l^2 exakt (keine zweite Variation von F nötig). Im (t, z)-Sektor hat L2/b "
      f"{len(Q2)} Monome in den 40 Feldern (16 h^a_mu, 24 tau^ab_mu). Die 40 Feldgleichungen sind von erster Ordnung "
      f"(max. Ordnung {ord_mm}); der Anteil eps d tau ^ d tau ist eine totale Ableitung und trägt nichts bei. In den "
      f"tau-Gleichungen stehen keine Ableitungen von tau ({not tau_abl}), in den h-Gleichungen keine von h ({not h_abl})",
      ord_mm == 1 and not tau_abl and not h_abl and len(Q2) > 0)

# A2 ---------------------------------------------------------------------------------------------------------------
XI = [{(40 + a, 0, 0): EINS} for a in range(4)]
LAM = [{(44 + pi, 0, 0): EINS} for pi in range(6)]
HG = [[None]*4 for _ in range(4)]
for a in range(4):
    for mu in range(4):
        terms = [(dcoord(XI[a], mu), 1)]
        for c in range(4):
            terms.append((XI[c], wmix(a, c, mu)))
            L_, s = anti(LAM, a, c)
            if s:
                terms.append((L_, -s*ETA[c]*ebar(c, mu)))
        HG[a][mu] = lsum(terms)
TG = [[None]*4 for _ in range(6)]
for pi, (a, b) in enumerate(PAIRS):
    for mu in range(4):
        terms = [(dcoord(LAM[pi], mu), 1)]
        for c in range(4):
            L1, s1 = anti(LAM, c, b)
            if s1:
                terms.append((L1, s1*wmix(a, c, mu)))
            L2_, s2 = anti(LAM, a, c)
            if s2:
                terms.append((L2_, s2*wmix(b, c, mu)))
        terms += [(XI[a], -ebar(b, mu)/l**2), (XI[b], ebar(a, mu)/l**2)]
        TG[pi][mu] = lsum(terms)
FG = build_f(HG, TG)
eich_null = all(len(v) == 0 for v in FG.values())
GAUGE = {}
for a in range(4):
    for mu in range(4):
        GAUGE[hJ(a, mu)] = HG[a][mu]
for pi in range(6):
    for mu in range(4):
        GAUGE[tJ(pi, mu)] = TG[pi][mu]
eich_gl = all(len(subst(e, GAUGE)) == 0 for e in EMM)
check("[E] Eichfreiheit: Unter den zehn linearisierten so(3,2)-Transformationen delta h^a = D xi^a - lambda^a_b e^b, "
      "delta tau^ab = D lambda^ab - (xi^a e^b - xi^b e^a)/l^2 mit beliebigen Funktionen xi^a(t, z) (4) und "
      f"lambda^ab(t, z) (6) ist f identisch null ({eich_null}); alle 40 Feldgleichungen sind auf reiner Eichung erfüllt "
      f"({eich_gl}). Die vier xi sind auf der Torsionsfläche die Diffeomorphismen, die sechs lambda die lokale "
      "Lorentzgruppe", eich_null and eich_gl)
check("[E/C] Was der Bestand ohne dynamisches Vierbein liefert: Nach z115 und z128 liegt jede native Familie in einer "
      "Stabilisatorbahn und ist Maurer-Cartan-flach. Ihre Tangentialrichtungen bei AdS4 sind also reine Eichung, und auf "
      "reiner Eichung ist f = 0. Die Hesse-Form der Uhrspurwirkung verschwindet damit auf allen nativen Richtungen "
      "identisch: Die Wirkung ist auf dem nativen Konfigurationsraum konstant (Euler-Zahl) und definiert dort keinen "
      "Propagator. Ein Pol kann erst entstehen, wenn Richtungen mit f != 0 zugelassen und mit der Wirkung gewichtet "
      "werden (P_dyn)", eich_null)

# A3 ---------------------------------------------------------------------------------------------------------------
MT = sp.zeros(24, 24)
for i in range(24):
    for (J, p, q), c in EMM[16 + i].items():
        if J >= 16:
            MT[i, J - 16] = c
rang_tau = MT.rank()
MTinv = MT.inv()
rhs = [{k: c for k, c in EMM[16 + i].items() if k[0] < 16} for i in range(24)]
TAU_SOL = {16 + j: lsum([(rhs[i], -MTinv[j, i]) for i in range(24) if MTinv[j, i] != 0]) for j in range(24)}
EH = [subst(EMM[J], TAU_SOL) for J in range(16)]
check(f"[E] Der Spinzusammenhang ist ein Hilfsfeld: Die 24 tau-Gleichungen sind algebraisch in tau, ihre 24 x 24-"
      f"Koeffizientenmatrix hat Rang {rang_tau} (Determinante 2^80/z^48 in Einheiten b^24). Sie sind die linearisierte "
      "Torsionsbedingung und legen tau eindeutig durch h und dessen erste Ableitungen fest. Eingesetzt bleiben 16 "
      f"Gleichungen zweiter Ordnung für h (max. Ordnung {ordnung(EH)})", rang_tau == 24 and ordnung(EH) == 2)

g0 = sp.diag(*[ETA[m]*l**2/z**2 for m in range(4)])
g0i = g0.inv()


def d0(expr, mu):
    return sp.diff(expr, z) if mu == 3 else 0


GAM0 = [[[sum(g0i[r, s]*(d0(g0[s, n], m) + d0(g0[s, m], n) - d0(g0[m, n], s)) for s in range(4))/2 for n in range(4)]
         for m in range(4)] for r in range(4)]
RIC0 = sp.zeros(4, 4)
for m in range(4):
    for n in range(4):
        RIC0[m, n] = sp.simplify(sum(d0(GAM0[r][m][n], r) - d0(GAM0[r][r][m], n)
                                     + sum(GAM0[r][r][s]*GAM0[s][m][n] - GAM0[r][n][s]*GAM0[s][r][m] for s in range(4))
                                     for r in range(4)))
LAMC = -3/l**2
hint_ok = sp.simplify(RIC0 - LAMC*g0) == sp.zeros(4, 4)


def linearisiere(g1):
    """Linearisierte Christoffel-, Riemann- und Ricci-Größen zur Störung g1 (Matrix von Jetformen)."""
    g1up = [[lsum([(g1[a][b], g0i[m, a]*g0i[n, b]) for a in range(4) for b in range(4)]) for n in range(4)]
            for m in range(4)]
    gam1 = [[[None]*4 for _ in range(4)] for _ in range(4)]
    for r in range(4):
        for m in range(4):
            for n in range(4):
                terms = []
                for s in range(4):
                    if g0i[r, s] != 0:
                        terms += [(dcoord(g1[s][n], m), g0i[r, s]/2), (dcoord(g1[s][m], n), g0i[r, s]/2),
                                  (dcoord(g1[m][n], s), -g0i[r, s]/2)]
                    c = (d0(g0[s, n], m) + d0(g0[s, m], n) - d0(g0[m, n], s))/2
                    if c != 0:
                        terms.append((g1up[r][s], -c))
                gam1[r][m][n] = lsum(terms)
    riem1 = {}
    for r in range(4):
        for s in range(4):
            for m in range(4):
                for n in range(m + 1, 4):
                    terms = [(dcoord(gam1[r][n][s], m), 1), (dcoord(gam1[r][m][s], n), -1)]
                    for k in range(4):
                        terms += [(gam1[k][n][s], GAM0[r][m][k]), (gam1[r][m][k], GAM0[k][n][s]),
                                  (gam1[k][m][s], -GAM0[r][n][k]), (gam1[r][n][k], -GAM0[k][m][s])]
                    riem1[(r, s, m, n)] = lsum(terms)
                    riem1[(r, s, n, m)] = lin_scale(riem1[(r, s, m, n)], -1)
                riem1[(r, s, m, m)] = {}
    ric1 = [[lsum([(riem1[(r, m, r, n)], 1) for r in range(4)]) for n in range(4)] for m in range(4)]
    return g1up, gam1, riem1, ric1


G1 = [[lsum([(H[m][n], ETA[m]*l/z), (H[n][m], ETA[n]*l/z)]) for n in range(4)] for m in range(4)]
G1UP, _, _, RIC1 = linearisiere(G1)
R1 = lsum([(RIC1[a][b], g0i[a, b]) for a in range(4) for b in range(4)])
TRG1 = lsum([(G1[a][b], g0i[a, b]) for a in range(4) for b in range(4)])
XUP = [[lsum([(RIC1[a][b], g0i[m, a]*g0i[n, b]) for a in range(4) for b in range(4)]
             + [(R1, -g0i[m, n]/2), (TRG1, LAMC*g0i[m, n]/2), (G1UP[m][n], -LAMC)]) for n in range(4)] for m in range(4)]
einst_ok, c0s = True, []
for a in range(4):
    for mu in range(4):
        ziel = lin_scale(XUP[mu][a], (l**4/z**4)*ETA[a]*l/z)        # sqrt(-g) delta(G + Lambda g)^{mu nu} e_{a nu}
        lhs = EH[hJ(a, mu)]
        for k in sorted(set(ziel) | set(lhs)):
            cz, cl = ziel.get(k, 0), lhs.get(k, 0)
            if cz == 0 or cl == 0:
                einst_ok = False
            else:
                c0s.append(sp.simplify(cl/cz))
c0 = c0s[0]
einst_ok &= all(sp.simplify(c - c0) == 0 for c in c0s)
ZAHLEN["A_einstein_abgleich"] = {"Proportionalität E_h/(sqrt(-g) delta(G + Lambda g) e), Einheiten b": str(c0),
                                 "verglichene Koeffizienten": len(c0s)}
check(f"[E] Nach Elimination von tau sind die 16 h-Gleichungen exakt die linearisierten Einstein-Gleichungen mit "
      f"Lambda = -3/l^2: delta S/delta h^a_mu = ({c0}) b sqrt(-g) delta(G^(mu nu) + Lambda g^(mu nu)) e_(a nu) für "
      f"delta g = e_a h^a + h^a e_a, Koeffizient für Koeffizient ({len(c0s)} Vergleiche; Hintergrund R_mn = Lambda g_mn: "
      f"{hint_ok}). Mit delta S_EH/delta h = -(1/(8 pi G)) sqrt(-g) delta(G + Lambda g) e folgt 1/(16 pi G) = 4 b/l^2, "
      "also b = l^2/(64 pi G) (z96) und mit b = 16 kappa (z127: 64 kappa/l^2 = 1/(16 pi G)) G/R_Vak^2 = "
      "1/(4096 pi kappa). Der antisymmetrische Teil von h fällt heraus (Lorentz-Eichung)",
      einst_ok and hint_ok and sp.simplify(c0 + 8/l**2) == 0)

# A4 ---------------------------------------------------------------------------------------------------------------
zahl_voll = zaehlung(EMM, list(range(40)), GAUGE, list(range(40, 50)), N=6, mmax=4)
zahl_voll5 = zaehlung(EMM, list(range(40)), GAUGE, list(range(40, 50)), N=5, mmax=3)
zuw_voll = [zahl_voll[i + 1][3] - zahl_voll[i][3] for i in range(len(zahl_voll) - 1)]


def rot_matrix():
    """Erzeuger der Drehung in der (x, y)-Ebene auf den 40 Feldkomponenten (innere und Formindizes)."""
    r = {1: {2: 1}, 2: {1: -1}}
    Rm = sp.zeros(40, 40)

    def comp_t(a, c, mu):
        if a == c:
            return None, 0
        return (tJ(PIDX[(a, c)], mu), 1) if a < c else (tJ(PIDX[(c, a)], mu), -1)
    for a in range(4):
        for mu in range(4):
            J = hJ(a, mu)
            for a2, cf in r.get(a, {}).items():
                Rm[hJ(a2, mu), J] += cf
            for m2, cf in r.get(mu, {}).items():
                Rm[hJ(a, m2), J] += cf
    for pi, (a, c) in enumerate(PAIRS):
        for mu in range(4):
            J = tJ(pi, mu)
            for a2, cf in r.get(a, {}).items():
                K, s = comp_t(a2, c, mu)
                if s:
                    Rm[K, J] += cf*s
            for c2, cf in r.get(c, {}).items():
                K, s = comp_t(a, c2, mu)
                if s:
                    Rm[K, J] += cf*s
            for m2, cf in r.get(mu, {}).items():
                Rm[tJ(pi, m2), J] += cf
    return Rm


RM = rot_matrix()
hel_dim = {q: len((RM*RM + q*q*sp.eye(40)).nullspace()) for q in (0, 1, 2, 3)}
NS2 = (RM*RM + 4*sp.eye(40)).nullspace()
h_in_2 = sum(1 for v in NS2 if any(v[i] != 0 for i in range(16)))
map2 = {J: lsum([({(100 + s, 0, 0): EINS}, v[J]) for s, v in enumerate(NS2) if v[J] != 0]) for J in range(40)}
E2 = [e for e in (subst(e, map2) for e in EMM) if e]
zahl_2 = zaehlung(E2, list(range(100, 106)), N=6, mmax=4)
zuw_2 = [zahl_2[i + 1][3] - zahl_2[i][3] for i in range(len(zahl_2) - 1)]
ZAHLEN["A_jetzaehlung"] = {
    "voll: (Grad m, Lösungsjets, Eichjets, Differenz), N = 6": [list(t) for t in zahl_voll],
    "voll: dasselbe mit N = 5": [list(t) for t in zahl_voll5],
    "Helizität +-2: (Grad m, Lösungsjets, 0, Differenz)": [list(t) for t in zahl_2],
    "Helizitätsräume (Ladung: Dimension)": {str(k): v for k, v in sorted(hel_dim.items())}}
check(f"[E] Freiheitsgrade durch Jetzählung (formale Potenzreihen im Punkt z = l = 1, exakte Ränge über Q): Die Zahl "
      f"der m-Jets von Lösungen modulo m-Jets reiner Eichung ist {[t[3] for t in zahl_voll]} für m = 0..4 (stabil "
      f"gegen die Fortsetzungsordnung: N = 5 gibt {[t[3] for t in zahl_voll5]}). Der Zuwachs je Grad ist {zuw_voll}: "
      "vier freie Funktionen einer Variablen, also zwei Wellengleichungen zweiter Ordnung, zwei propagierende "
      "Freiheitsgrade. Von 40 Feldkomponenten bleiben nach 24 Hilfsfeldern, 10 Eichfunktionen und 4 Zwangsbedingungen "
      "genau 2", zuw_voll == [4, 4, 4, 4] and [t[3] for t in zahl_voll5] == [t[3] for t in zahl_voll[:4]])
check(f"[E] Helizitäten: Unter der Drehung in der (x, y)-Ebene zerfallen die 40 Komponenten in die Ladungen 0, +-1, +-2 "
      f"mit den Dimensionen {hel_dim[0]}, {hel_dim[1]}, {hel_dim[2]} (Ladung 3: {hel_dim[3]}). Der Ladung-2-Sektor "
      f"enthält {h_in_2} Vierbein- und {hel_dim[2] - h_in_2} Zusammenhangskomponenten und keinen Eichparameter (die zehn "
      f"Parameter tragen Ladung 0 und +-1). Seine Jetzählung ist {[t[3] for t in zahl_2]} = 6 + 4 m, Zuwachs "
      f"{zuw_2}: genau zwei Wellengleichungen (je 2 m + 1 Jets des Felds und 2 weitere des Hilfsfelds). Der Zuwachs "
      "der Sektoren 0 und +-1 ist damit null; der konstante Rest 9 - 6 = 3 sind homogene Moden der (t, z)-Reduktion "
      "ohne Zuwachs. Die zwei Freiheitsgrade sind ein Dublett der Ladung 2: Helizität +2 und -2",
      hel_dim == {0: 14, 1: 20, 2: 6, 3: 0} and h_in_2 == 2 and zuw_2 == [4, 4, 4, 4]
      and [t[3] for t in zahl_2] == [6, 10, 14, 18, 22])

# A5 ---------------------------------------------------------------------------------------------------------------
PHI = 200                                             # delta g_xy = (l^2/z^2) phi
map_phi = {J: {} for J in range(16)}
map_phi[hJ(1, 2)] = {(PHI, 0, 0): l/(2*z)}
map_phi[hJ(2, 1)] = {(PHI, 0, 0): l/(2*z)}
E_phi = lsum([(subst(EH[hJ(1, 2)], map_phi), l/(2*z)), (subst(EH[hJ(2, 1)], map_phi), l/(2*z))])
ziel_phi = {(PHI, 2, 0): -4/z**2, (PHI, 0, 2): 4/z**2, (PHI, 0, 1): -8/z**3}
phi_ok = lsum([(E_phi, 1), (ziel_phi, -1)]) == {}
andere = all(len(subst(EH[J], map_phi)) == 0 for J in range(16) if J not in (hJ(1, 2), hJ(2, 1)))
Dl = sp.symbols("Delta")
ind0 = sp.factor(sp.simplify((z**2*sp.diff(z**Dl, z, 2) - 2*z*sp.diff(z**Dl, z))/z**Dl))
w_, k_ = sp.symbols("omega k", positive=True)
mode_ein = (1 - sp.I*w_*z)*sp.exp(sp.I*w_*z)
mode_ok = sp.simplify(sp.diff(mode_ein, z, 2) - 2*sp.diff(mode_ein, z)/z + w_**2*mode_ein) == 0
check(f"[E] Helizität +-2 explizit: Für delta g_xy = (l^2/z^2) phi(t, z) (ebenso für delta g_xx - delta g_yy) ist "
      f"delta S/delta phi = 4 b z^-2 (-phi_tt + phi_zz - (2/z) phi_z) ({phi_ok}; alle anderen Gleichungen sind "
      f"identisch erfüllt: {andere}). Das ist die Euler-Lagrange-Gleichung von L = 2 b z^-2 (phi_t^2 - phi_z^2), die "
      f"masselose Wellengleichung auf AdS4. Exponenten am Rand: {ind0} = 0, also Delta = 0 und Delta = 3 = s + 1, die "
      f"Unitaritätsschranke des Spin 2 (masselose Darstellung D(3,2)); Moden (1 -+ i omega z) e^(+-i omega z) "
      f"({mode_ok})", phi_ok and andere and mode_ok and sp.simplify(ind0 - Dl*(Dl - 3)) == 0)
kap = sp.symbols("kappa", positive=True)
res_mm = sp.simplify(32*sp.pi*(l**2/(1024*sp.pi*kap)))
ZAHLEN["A_residuum"] = {"kinetischer Koeffizient": "2 b z^-2 (phi_t^2 - phi_z^2), b = l^2/(64 pi G) = 16 kappa",
                        "Residuum im flachen Grenzfall omega l >> 1": "32 pi G = l^2/(32 kappa) = R_Vak^2/(128 kappa)",
                        "G/R_Vak^2 bei kappa = 1": gs(1/(4096*mp.pi)), "G/R_Vak^2 bei kappa = 2 c_3 = 1/(4 pi)": gs(1/1024)}
check(f"[C, P_dyn] Residuum: Der kinetische Koeffizient ist 2 b z^-2 mit b = l^2/(64 pi G) = 16 kappa. Für "
      f"omega l >> 1 ist das lokal der Propagator 32 pi G/p^2 mit Residuum 32 pi G = {res_mm} = R_Vak^2/(128 kappa). "
      "Positiv ist es genau für b > 0, also für kappa mal Orientierung > 0: Das Vorzeichen ist das Paarungsbit aus "
      "z118/z127 (tau-Parität des Schnitts), keine neue Größe. Mit kappa = 1 ist G/R_Vak^2 = "
      f"{gs(1/(4096*mp.pi))}, mit 2 c_3 {gs(1/1024)} (z127 reproduziert)", sp.simplify(res_mm - l**2/(32*kap)) == 0)
befund("Spin-2-Pol der Uhrspurwirkung", "bedingt (P_dyn)",
       "Teil A: Ist das Vierbein dynamisch, so hat die linearisierte MacDowell-Mansouri-Wirkung mit den nativen "
       "Spurgewichten genau einen masselosen Spin-2-Sektor D(3,2) mit zwei Helizitäten; der Spinzusammenhang ist "
       "Hilfsfeld, zehn Eichfunktionen, Residuum 32 pi G = l^2/(32 kappa), positiv bei richtigem Paarungsbit. Ohne "
       "P_dyn verschwindet die Hesse-Form auf allen nativen Richtungen.")

# A6 ---------------------------------------------------------------------------------------------------------------
# so(3,2)-Inhalt: Charaktere als {Energie: {Spin: Vielfachheit}}
EMAX = 7


def cg(s1, s2):
    return list(range(abs(s1 - s2), s1 + s2 + 1))


def ch_mul(A, B):
    out = {}
    for e1, c1 in A.items():
        for e2, c2 in B.items():
            if e1 + e2 > EMAX:
                continue
            tgt = out.setdefault(e1 + e2, {})
            for s1, m1 in c1.items():
                for s2, m2 in c2.items():
                    for s in cg(s1, s2):
                        tgt[s] = tgt.get(s, 0) + m1*m2
    return out


def ch_add(A, B, c=1):
    out = {e: dict(v) for e, v in A.items()}
    for e, v in B.items():
        tgt = out.setdefault(e, {})
        for s, m in v.items():
            tgt[s] = tgt.get(s, 0) + c*m
    return {e: {s: m for s, m in v.items() if m} for e, v in out.items() if any(v.values())}


PIMP = {n: {s: 1 for s in range(n % 2, n + 1, 2)} for n in range(EMAX + 1)}       # Impulse: Sym^n(Spin 1)


def ch_modul(E0, s, masselos=False):
    kopf = {E0: {s: 1}}
    out = ch_mul(kopf, PIMP)
    if masselos:
        out = ch_add(out, ch_mul({E0 + 1: {s - 1: 1}}, PIMP), -1)
    return out


def zerlege(chi):
    """Zerlegt einen Charakter in tiefste Gewichte (E0, s); E0 = s + 1 mit s >= 1 wird als masselos abgezogen."""
    chi = {e: dict(v) for e, v in chi.items()}
    out = []
    for e in range(0, EMAX + 1):
        while chi.get(e):
            s = max(chi[e])
            m = chi[e][s]
            masselos = (e == s + 1 and s >= 1)
            chi = ch_add(chi, ch_modul(e, s, masselos), -m)
            if any(v < 0 for lev in chi.values() for v in lev.values()):
                return None
            out.append((e, s, m, masselos))
    return out


leiter = {m: {s: 1 for s in range(m)} for m in range(1, EMAX + 1)}                # (j, j) bei E = 2j + 1 = m
zer1 = zerlege(leiter)
zwei = ch_mul(leiter, leiter)
zer2 = zerlege(zwei)
drei = ch_mul(zwei, leiter)
zer3 = zerlege(drei)
spin2_min = {n: min(e for e, s, m, _ in zz if s == 2) for n, zz in ((1, [(3, 2, 1, False)]), (2, zer2), (3, zer3))}
d32 = {n: sum(m for e, s, m, _ in zz if (e, s) == (3, 2)) for n, zz in ((1, zer1), (2, zer2), (3, zer3))}
ZAHLEN["A_so32_inhalt"] = {
    "eine Leiter": [list(t[:3]) for t in zer1],
    "zwei Leitern, tiefste Gewichte bis E = 5": [list(t[:3]) for t in zer2 if t[0] <= 5],
    "drei Leitern, tiefste Gewichte bis E = 5": [list(t[:3]) for t in zer3 if t[0] <= 5],
    "Vielfachheit von D(3,2) für 1, 2, 3 Leitern": [d32[1], d32[2], d32[3]]}
check(f"[E, neu] so(3,2)-Inhalt des nativen Fockraums (so(3,2) = Der(I_AdS) in su(2,2), E gemeinsam, SO(3) diagonal in "
      f"SO(4)): Die native Leiter (Stufen (j, j), E = 2j + 1) zerfällt in {[(e, s) for e, s, m, _ in zer1]} = "
      "D(1,0) + D(2,0), den konform gekoppelten Skalar auf AdS4 mit beiden Randverhalten (so schon z104 über den "
      "Casimir; hier über Charaktere). Neu ist der Mehrleiterinhalt: In zwei Leitern sind die "
      f"tiefsten Gewichte mit Spin 2 {[(e, s, m) for e, s, m, _ in zer2 if s == 2 and e <= 5]} (E0, s, Vielfachheit): "
      "der erste liegt bei E0 = 4, das ist z137s Spin-2-Primär (Delta = 4), auf AdS4 ein massives Kompositum "
      "(m^2 l^2 = E0 (E0 - 3) = 4). Die masselose Darstellung D(3,2) hat die Vielfachheit "
      f"{[d32[1], d32[2], d32[3]]} in 1, 2, 3 Leitern (Charaktere bis E = {EMAX}). Allgemein: n Quanten haben Spin "
      "<= E - n; Spin 2 bei E = 3 gibt es nur für n = 1, dort als Nachkomme von D(1,0). Der native Fockraum enthält "
      "kein D(3,2)", zer1 is not None and [(e, s) for e, s, m, _ in zer1] == [(1, 0), (2, 0)]
      and zer2 is not None and zer3 is not None and d32 == {1: 0, 2: 0, 3: 0}
      and min(e for e, s, m, _ in zer2 if s == 2) == 4)

# Unitaritätsschranke: Gram-Form der ersten Stufe über einem tiefsten Gewicht (E0, s = 2)
ETA5 = sp.diag(-1, -1, 1, 1, 1)


def mgen(A, B):
    m = sp.zeros(5, 5)
    for C in range(5):
        for D in range(5):
            m[C, D] = (1 if C == A else 0)*ETA5[B, D] - (1 if C == B else 0)*ETA5[A, D]
    return m


def komm(A, B):
    return A*B - B*A


E_so = sp.I*mgen(0, 1)
J_so = [-sp.I*mgen(3, 4), -sp.I*mgen(4, 2), -sp.I*mgen(2, 3)]
LP = [sp.I*mgen(2 + i, 0) - mgen(2 + i, 1) for i in range(3)]
LM = [sp.I*mgen(2 + i, 0) + mgen(2 + i, 1) for i in range(3)]
alg_ok = all(komm(E_so, LP[i]) == LP[i] and komm(E_so, LM[i]) == -LM[i] for i in range(3)) \
    and komm(J_so[0], J_so[1]) == sp.I*J_so[2] and komm(J_so[1], J_so[2]) == sp.I*J_so[0]
lc = sp.LeviCivita
struktur_ok = all(komm(LM[i], LP[j]) == 2*(1 if i == j else 0)*E_so
                  - 2*sp.I*sum((lc(i, j, k)*J_so[k] for k in range(3)), sp.zeros(5, 5))
                  for i in range(3) for j in range(3))
E0s = sp.symbols("E_0")


def spinmat(s):
    d = 2*s + 1
    Jp, Jm, Jz = sp.zeros(d, d), sp.zeros(d, d), sp.zeros(d, d)
    for i, m in enumerate(range(-s, s + 1)):
        Jz[i, i] = m
        if m < s:
            Jp[i + 1, i] = s - m
        if m > -s:
            Jm[i - 1, i] = s + m
    return [(Jp + Jm)/2, (Jp - Jm)/(2*sp.I), Jz]


def gram(s):
    Js = spinmat(s)
    d = 2*s + 1
    G = sp.zeros(3*d, 3*d)
    for i in range(3):
        for j in range(3):
            blk = 2*(1 if i == j else 0)*E0s*sp.eye(d) - 2*sp.I*sum((lc(i, j, k)*Js[k] for k in range(3)), sp.zeros(d, d))
            G[i*d:(i + 1)*d, j*d:(j + 1)*d] = blk
    return G


G2 = gram(2)
lams = [2*(E0s - 3), 2*(E0s - 1), 2*(E0s + 2)]
minpol = sp.simplify((G2 - lams[0]*sp.eye(15))*(G2 - lams[1]*sp.eye(15))*(G2 - lams[2]*sp.eye(15))) == sp.zeros(15, 15)
vielf = [15 - (G2 - lam*sp.eye(15)).subs(E0s, 7).rank() for lam in [la.subs(E0s, 7) for la in lams]]
check(f"[E] Unitaritätsschranke der nativen so(3,2) (5 x 5-Realisierung: [E, L^+-] = +-L^+-, [L^-_i, L^+_j] = "
      f"2 delta_ij E - 2 i eps_ijk J_k: {alg_ok and struktur_ok}): Über einem tiefsten Gewicht (E0, Spin 2) hat die "
      f"Gram-Form der ersten Stufe die Eigenwerte 2 (E0 - 3), 2 (E0 - 1), 2 (E0 + 2) mit Vielfachheiten {vielf} "
      f"(Minimalpolynom: {minpol}). Positiv ist sie genau für E0 >= 3 = s + 1; bei E0 = 3 entsteht der Nullvektor der "
      "Eichung (masselos), bei E0 = 2 und E0 = 1 gibt es Zustände negativer Norm. Ein unitärer Spin-2-Sektor "
      "positiver Energie beginnt also frühestens bei D(3,2)",
      alg_ok and struktur_ok and minpol and vielf == [3, 5, 7])
befund("Masseloser Spin-2-Zustand im Bestand ohne P_dyn", "nein",
       "Teil A6: Der native Fockraum (Leitern) enthält unter so(3,2) kein D(3,2); Spin 2 beginnt als Kompositum bei "
       "E0 = 4. Das gilt für die Helizität-0-Leiter; z75s Dubleton-Zustände mit Helizität +-2 liegen außerhalb "
       "(Vierquantenkomposita der Oszillatoren, nicht Dynamik des Originals).")

# ====================================================================================================================
say("\nTeil B: Induzierter Propagator aus dem nativen <TT>")
# B1 ---------------------------------------------------------------------------------------------------------------
epsr = sp.symbols("epsilon", positive=True)


def n_ft(a):
    """Fourier-Transformierte von (x^2)^-a in vier Dimensionen: n_ft(a) (p^2)^(a - 2)."""
    return sp.pi**2*2**(4 - 2*a)*sp.gamma(2 - a)/sp.gamma(a)


kft = sp.Rational(1, 2)*(n_ft(4 - epsr) + 2*2*2*(3 - epsr)*n_ft(5 - epsr) + 8*4*(4 - epsr)*(3 - epsr)*n_ft(6 - epsr))
pol_ft = sp.simplify(sp.series(kft, epsr, 0, 1).removeO().coeff(epsr, -1))
CT_nat = 1/(3*sp.pi**4)
c_nat = sp.simplify(sp.pi**4*CT_nat/40)
check(f"[E] Koeffizient des Logarithmus: <T_xy(x) T_xy(0)> = C_T I_(xy,xy)/x^8 mit I_(xy,xy) = (1/2)(1 - 2 (x_1^2 + "
      f"x_2^2)/x^2 + 8 x_1^2 x_2^2/x^4) hat für Impuls senkrecht zur (1,2)-Ebene die Fourier-Transformierte "
      f"-({pol_ft}) C_T p^4 ln p^2 (analytische Fortsetzung im Exponenten, Polterm 1/epsilon). Mit c = pi^4 C_T/40 ist das "
      f"-(c/(16 pi^2)) p^4 ln p^2, und die quadratische effektive Wirkung des Tensorkanals ist (1/2) h Pi h mit "
      f"Pi(p^2) = (c/(32 pi^2)) p^4 ln(p^2/mu^2), gleich z137s Gamma = -(L/(16 pi^2)) c Int W^2. Nativ c = {c_nat} je "
      "Leiter", sp.simplify(pol_ft - sp.pi**2/640) == 0 and c_nat == sp.Rational(1, 120))

# B2 ---------------------------------------------------------------------------------------------------------------
xs = sp.symbols("x", positive=True)                   # x = p^2/Lambda_W^2 (euklidisch)
Da = 1/(xs**2*sp.log(xs))
lim0 = sp.limit(xs*Da, xs, 0, "+")
res_a = sp.limit((xs - 1)*Da, xs, 1)
rr = sp.symbols("r", positive=True)
pot = -rr/(8*sp.pi)
lap = sp.simplify(sp.diff(rr**2*sp.diff(pot, rr), rr)/rr**2)
mono_a = sp.simplify(sp.diff(xs*Da, xs))
check(f"[E] (a) Ohne Baumterm ist der resummierte Propagator D = (32 pi^2/c)/(p^4 ln(p^2/Lambda_W^2)). Bei p^2 = 0 "
      f"gibt es keinen Pol: p^2 D geht gegen {lim0} (kein endliches Residuum), und D ist für p < Lambda_W negativ. Die "
      f"einzige Polstelle liegt bei euklidischem p^2 = Lambda_W^2 (Residuum {res_a} in Einheiten 32 pi^2/(c "
      "Lambda_W^2)), also bei raumartigem Impuls: ein Tachyon. Die Spektraldichte auf dem zeitartigen Schnitt, "
      "(32 pi^2/c)/(s^2 (ln^2(s/Lambda_W^2) + pi^2)), ist positiv, aber bei s -> 0 nicht integrierbar; p^2 D ist nicht "
      "monoton wachsend, wie es ein positives Källén-Lehmann-Maß verlangt. Statisch gibt 1/k^4 das Potential "
      f"-r/(8 pi) (Laplace davon: {lap} = -1/(4 pi r)), kein 1/r-Gesetz. Der Bestand ohne Baumterm hat also keinen "
      "masselosen Spin-2-Pol (z137 bestätigt), sondern eine Antwort ohne Teilchendeutung",
      lim0 == -sp.oo and res_a == 1 and sp.simplify(lap + 1/(4*sp.pi*rr)) == 0
      and mono_a.subs(xs, sp.Rational(1, 2)) < 0 and mono_a.subs(xs, sp.Rational(1, 10)) > 0)


# B3 ---------------------------------------------------------------------------------------------------------------
def pole(a):
    """Nullstellen von a + x ln x (x = p^2/Lambda_W^2) auf dem Hauptblatt: x = exp(W_k(-a)), k = 0, -1."""
    us = [mp.lambertw(-a, 0), mp.lambertw(-a, -1)]
    return [(mp.exp(u), 1/(mp.exp(u)*(u + 1))) for u in us]


def rho_int(a):
    f = lambda t: mp.exp(t)/((a - mp.exp(t)*t)**2 + mp.pi**2*mp.exp(2*t))        # s = e^t
    return mp.quad(f, [-mp.inf, -40, -10, -2, 0, 2, 10, 40, mp.inf])


tab_b = {}
sum_ok = True
for a_ in (mp.mpf("0.2"), mp.mpf(5), mp.mpf(1000)):
    ps = pole(a_)
    rest = max(abs(a_ + x*mp.log(x)) for x, _ in ps)
    summe = 1/a_ + sum(zr for _, zr in ps) + rho_int(a_)
    sum_ok &= rest < mp.mpf(10)**-20 and abs(summe) < mp.mpf(10)**-8
    tab_b[gs(a_, 4)] = {"Pole x = p^2/Lambda_W^2": [f"{gs(x.real)} + {gs(x.imag)} i" for x, _ in ps],
                        "Residuen": [f"{gs(zr.real)} + {gs(zr.imag)} i" for _, zr in ps],
                        "masseloses Residuum 1/a": gs(1/a_), "Integral der Spektraldichte": gs(rho_int(a_)),
                        "Summenregel (soll 0)": gs(abs(summe), 2) if abs(summe) > 1e-12 else "0"}
a_kr = sp.exp(-1)
zweig = all(abs(mp.lambertw(-mp.mpf(5), kk).imag) > mp.pi for kk in (1, -2, 2, -3))
haupt = all(abs(mp.lambertw(-mp.mpf(5), kk).imag) < mp.pi for kk in (0, -1))
p_tach = pole(mp.mpf("0.2"))
p_kompl = pole(mp.mpf(5))
ZAHLEN["B_pole"] = tab_b
check(f"[E/N] (b) Mit Baumterm (P_dyn) ist D = 1/(p^2 g), g = 1/(32 pi G) + (c/(32 pi^2)) p^2 ln(p^2/Lambda_W^2). Bei "
      f"p^2 = 0 liegt der masselose Pol mit Residuum 32 pi G > 0. Die weiteren Pole lösen x ln x = -a mit "
      f"a = pi/(c G Lambda_W^2), exakt x = exp(W_k(-a)), auf dem Hauptblatt genau k = 0 und k = -1 ({haupt}; andere "
      f"Zweige liegen daneben: {zweig}). Für a < 1/e = {gs(a_kr)} sind es zwei Tachyonen (a = 0,2: x = "
      f"{gs(p_tach[0][0].real)} mit Residuum {gs(p_tach[0][1].real)} und x = {gs(p_tach[1][0].real)} mit "
      f"{gs(p_tach[1][1].real)}, eines davon ein Geist), für a > 1/e ein konjugiert komplexes Paar (a = 5: x = "
      f"{gs(p_kompl[0][0].real)} +- {gs(abs(p_kompl[0][0].imag))} i, Residuum {gs(p_kompl[0][1].real)} -+ "
      f"{gs(abs(p_kompl[0][1].imag))} i). Summenregel: Weil p^2 D im Unendlichen verschwindet, ist masseloses Residuum + "
      "Polresiduen + Integral der (positiven) Spektraldichte = 0, numerisch bestätigt für a = 0,2, 5, 1000. Das "
      "positive Graviton und das positive Kontinuum erzwingen also Pole mit negativem Realteil des Residuums: Das "
      "resummierte Spektrum ist für kein G und kein Lambda_W unitär; als Niederenergietheorie gilt es bis "
      "p^2 von der Ordnung des Polbetrags",
      sum_ok and haupt and zweig and p_tach[0][1].real > 0 and p_tach[1][1].real < 0
      and abs(p_kompl[0][0].imag) > 0 and p_kompl[0][1].real < 0)


def a_nativ(kappa, c, y):
    return 1024*mp.pi**2*kappa/(c*y**2)                # a = pi/(c G Lambda_W^2) mit G = l^2/(1024 pi kappa), y = Lambda_W l


y_kr = mp.sqrt(1024*mp.pi**2*mp.e*120)
nat_b = {}
for y in (1, 100, 10000):
    a_ = a_nativ(1, mp.mpf(1)/120, mp.mpf(y))
    ps = pole(a_)
    nat_b[f"Lambda_W l = {y}"] = {"a": gs(a_), "Art": "komplexes Paar" if abs(ps[0][0].imag) > 0 else "zwei Tachyonen",
                                  "|p| l der Pole": [gs(mp.sqrt(abs(x))*y) for x, _ in ps]}
ZAHLEN["B_nativ"] = {"Schwelle Lambda_W l für Tachyonen (kappa = 1, c = 1/120)": gs(y_kr), "Fälle": nat_b}
check(f"[C] Native Zahlen für (b): Mit G = l^2/(1024 pi kappa) ist a = 1024 pi^2 kappa/(c (Lambda_W l)^2); für "
      f"kappa = 1, c = 1/120 liegt die Schwelle a = 1/e bei Lambda_W l = {gs(y_kr)}. Darunter (etwa Lambda_W l = 1: "
      f"|p| l = {nat_b['Lambda_W l = 1']['|p| l der Pole'][0]}) ist es das komplexe Paar, darüber (Lambda_W l = 10^4: "
      f"|p| l = {nat_b['Lambda_W l = 10000']['|p| l der Pole']}) zwei Tachyonen. Die Pole liegen bei p l >> 1, die "
      "flache Analyse ist dort selbstkonsistent. Welche Art vorliegt, entscheidet allein die freie Zahl Lambda_W l",
      nat_b["Lambda_W l = 1"]["Art"] == "komplexes Paar" and nat_b["Lambda_W l = 10000"]["Art"] == "zwei Tachyonen")

# B4 ---------------------------------------------------------------------------------------------------------------
al0, cc, mu0, sig, p2s = sp.symbols("alpha_0 c mu_0 sigma p2", positive=True)
aW = lambda a0, m: a0 + (cc/(32*sp.pi**2))*sp.log(p2s/m**2)
rg_ok = sp.simplify(sp.expand_log(aW(al0, mu0) - aW(al0 + (cc/(16*sp.pi**2))*sig, mu0*sp.exp(sig)), force=True)) == 0
g2s, ks = sp.symbols("g2 k", positive=True)
LamW = mu0*sp.exp(-16*sp.pi**2*(-1/(4*g2s*ks))/cc)
expo = {f"c = {cn}, g^2 k = {gn}": gs(4*mp.pi**2/(mp.mpf(cv)*gv)/mp.log(10), 5)
        for cn, cv in (("1/120", mp.mpf(1)/120), ("6/5", mp.mpf(6)/5), ("12/5", mp.mpf(12)/5))
        for gn, gv in (("2", 2), ("1/4", mp.mpf(1)/4))}
ZAHLEN["B_skala"] = {"Lambda_W/mu_0 = exp(4 pi^2/(c g^2 k)), log10": expo}
check(f"[E/C] Die Skala des Logarithmus: mu ist eine Kontaktterm-Freiheit, mu -> mu e^sigma ist dasselbe wie "
      f"alpha_0 -> alpha_0 + c sigma/(16 pi^2) ({rg_ok}). Invariant ist allein Lambda_W = mu_0 exp(-16 pi^2 alpha_0/c), "
      "die Skala, bei der der laufende Weyl-Koeffizient null wird; in (a) ist sie die Tachyonmasse. Nativ ist sie "
      "nicht: Die konforme Leiter trägt c, aber keinen Kontaktterm, und die einzige native Länge ist l = R_Vak/2 "
      "(z112), sodass Lambda_W l eine freie reine Zahl bleibt (P_cut). Setzt man den nativen Baumwert aus Teil C, "
      f"alpha_0 = -1/(4 g^2 k) bei mu_0, so ist Lambda_W/mu_0 = exp(4 pi^2/(c g^2 k)) (log10: {expo}): eine "
      "dimensionale Transmutation mit freiem mu_0 und freiem g, keine native Zahl",
      rg_ok and sp.simplify(sp.log(LamW/mu0) - 4*sp.pi**2/(cc*g2s*ks)) == 0)
befund("Induzierter Propagator", "kein masseloser Pol ohne Baumterm; mit Baumterm nicht unitär",
       "Teil B: (a) Ohne Baumterm kein masseloser Pol, ein Tachyon bei der Logarithmusskala, kein Newton-Gesetz. (b) Mit "
       "dem MacDowell-Mansouri-Baumterm ein masseloser Pol mit Residuum 32 pi G, dazu zwangsläufig ein komplexes "
       "Paar oder zwei Tachyonen (Summenregel). Die Skala Lambda_W ist nicht nativ (P_cut).")

# ====================================================================================================================
say("\nTeil C: Konforme Gravitation mit I_AdS und Randbedingung")
# C1 ---------------------------------------------------------------------------------------------------------------
g1w = [[{} for _ in range(4)] for _ in range(4)]
g1w[1][2] = {(PHI, 0, 0): l**2/z**2}
g1w[2][1] = {(PHI, 0, 0): l**2/z**2}
G1WUP, _, RIEM1, RIC1W = linearisiere(g1w)
R1W = lsum([(RIC1W[a][b], g0i[a, b]) for a in range(4) for b in range(4)])
TRW = lsum([(g1w[a][b], g0i[a, b]) for a in range(4) for b in range(4)])
dR = lsum([(R1W, 1), (TRW, -LAMC)])
P0 = (LAMC/6)*g0
P1 = [[lsum([(RIC1W[m][n], sp.Rational(1, 2)), (dR, -g0[m, n]/12), (g1w[m][n], -4*LAMC/12)]) for n in range(4)]
      for m in range(4)]


def riem0(lam, s, m, n):
    return -((1 if lam == m else 0)*g0[s, n] - (1 if lam == n else 0)*g0[s, m])/l**2


WEYL1 = {}
for r in range(4):
    for s in range(4):
        for m in range(4):
            for n in range(4):
                terms = [(RIEM1[(k, s, m, n)], g0[r, k]) for k in range(4) if g0[r, k] != 0]
                terms += [(g1w[r][k], riem0(k, s, m, n)) for k in range(4)]
                terms += [(g1w[r][m], -P0[s, n]), (P1[s][n], -g0[r, m]), (g1w[r][n], P0[s, m]), (P1[s][m], g0[r, n]),
                          (g1w[s][m], P0[r, n]), (P1[r][n], g0[s, m]), (g1w[s][n], -P0[r, m]), (P1[r][m], -g0[s, n])]
                WEYL1[(r, s, m, n)] = lsum(terms)
spurfrei = all(len(lsum([(WEYL1[(r, s, r, n)], g0i[r, r]) for r in range(4)])) == 0 for s in range(4) for n in range(4))
QW = {}
for key, wv in WEYL1.items():
    if wv:
        r, s, m, n = key
        quad_add(QW, wv, wv, (l**4/z**4)*g0i[r, r]*g0i[s, s]*g0i[m, m]*g0i[n, n])
QW = quad_clean(QW)                                    # sqrt(-g) W^2 in zweiter Ordnung
EW = euler(QW, PHI)                                    # delta/delta phi von Int sqrt(-g) W^2
ziel_w = {(PHI, 4, 0): sp.Integer(2), (PHI, 2, 2): sp.Integer(-4), (PHI, 0, 4): sp.Integer(2)}
weyl_ok = lsum([(EW, 1), (ziel_w, -1)]) == {}
fT = sp.Function("f")(sp.Symbol("t"), z)
ts = sp.Symbol("t")
D0op = lambda u: z**2*(sp.diff(u, z, 2) - 2*sp.diff(u, z)/z - sp.diff(u, ts, 2))
box = lambda u: sp.diff(u, z, 2) - sp.diff(u, ts, 2)
fakt_ok = sp.simplify(D0op(D0op(fT) + 2*fT) - z**4*box(box(fT))) == 0
ind2 = sp.factor(sp.simplify((z**2*sp.diff(z**Dl, z, 2) - 2*z*sp.diff(z**Dl, z) + 2*z**Dl)/z**Dl))
check(f"[E] Linearisierter Weyl-Operator auf AdS4 im Helizität-2-Sektor (delta g_xy = (l^2/z^2) phi; linearisierter "
      f"Weyl-Tensor spurfrei: {spurfrei}): sqrt(-g) W^2 = (Box phi)^2 mit Box = d_z^2 - d_t^2 bis auf totale "
      f"Ableitungen, Euler-Lagrange-Ausdruck 2 Box^2 phi ohne jede z-Abhängigkeit und ohne l ({weyl_ok}): die "
      "Weyl-Invarianz, explizit. Mit dem masselosen AdS-Operator D_0 = z^2 (d_z^2 - (2/z) d_z - d_t^2) aus A5 gilt "
      f"z^4 Box^2 = D_0 (D_0 + 2) ({fakt_ok}): masseloses Graviton mal teilweise masseloses Spin-2-Feld mit "
      f"m^2 l^2 = -2, Randexponenten {ind2} = 0, also Delta = 1, 2. Lorentzsch ist delta S_W/delta phi = "
      "-2 alpha_W z^-4 D_0 (D_0 + 2) phi", spurfrei and weyl_ok and fakt_ok and sp.simplify(ind2 - (Dl - 1)*(Dl - 2)) == 0)
Xs = sp.symbols("X")
pb = sp.apart(1/(Xs*(Xs + 2)), Xs)
check(f"[E] Aus dem Doppelpol wird ein einfacher Pol, sobald l endlich ist: 1/(D_0 (D_0 + 2)) = {pb}. Nahe der "
      "masselosen Schale ist der Weyl-Kern -4 alpha_W z^-4 D_0, der Einstein-Kern aus A5 ist 4 b z^-4 D_0. Der "
      "masselose Anteil der konformen Gravitation ist also Einstein-Gravitation mit b = -alpha_W, das heißt "
      "l^2/(64 pi G) = -alpha_W (z137: alpha_W = 3/(64 pi G Lambda)). Der teilweise masselose Anteil trägt das "
      "entgegengesetzte Vorzeichen (Geist relativ zum Graviton). Im flachen Grenzfall l -> unendlich verschmelzen "
      "beide Pole zu 1/p^4, und G alpha_W/l^2 = -1/(64 pi) lässt bei festem alpha_W kein endliches G übrig",
      sp.simplify(pb - (1/(2*Xs) - 1/(2*(Xs + 2)))) == 0)

# C2 ---------------------------------------------------------------------------------------------------------------
m_e_g = sp.cos(w_*z) + w_*z*sp.sin(w_*z)               # Einstein, gerade, Delta = 0
m_e_u = sp.sin(w_*z) - w_*z*sp.cos(w_*z)               # Einstein, ungerade, Delta = 3
m_p_u = z*sp.cos(w_*z)                                 # teilweise masselos, Delta = 1
m_p_g = z*sp.sin(w_*z)                                 # teilweise masselos, Delta = 2
d0w = lambda u: z**2*(sp.diff(u, z, 2) - 2*sp.diff(u, z)/z + w_**2*u)
moden = {"Einstein gerade": m_e_g, "Einstein ungerade": m_e_u, "PM ungerade": m_p_u, "PM gerade": m_p_g}
tab_m = {}
for name, f in moden.items():
    ser = sp.series(f, z, 0, 5).removeO()
    fuehr = min(i for i in range(5) if ser.coeff(z, i) != 0)
    tab_m[name] = {"löst D_0": sp.simplify(d0w(f)) == 0, "löst D_0 + 2": sp.simplify(d0w(f) + 2*f) == 0,
                   "führende Potenz": fuehr, "phi'(0) = 0": sp.diff(f, z).subs(z, 0) == 0,
                   "Parität unter z -> -z": "gerade" if sp.simplify(f.subs(z, -z) - f) == 0 else "ungerade"}
wron = sp.simplify(sp.Matrix([[sp.diff(f, z, i) for f in moden.values()] for i in range(4)]).det())
neumann_raum = [n for n, t in tab_m.items() if t["phi'(0) = 0"]]
ZAHLEN["C_moden"] = {k: {kk: str(vv) for kk, vv in v.items()} for k, v in tab_m.items()}
check(f"[E, neu] Moden bei zeitartigem Randimpuls (phi ~ e^(-i omega t)): Der Lösungsraum von Box^2 phi = 0 ist "
      f"vierdimensional (Wronski-Determinante {wron}); alle vier Lösungen sind bei z = 0 analytisch, der Rand ist für "
      "die konforme Gravitation durchlässig. Einstein (D_0): cos + omega z sin (Delta = 0, gerade) und sin - omega z "
      "cos (Delta = 3, ungerade). Teilweise masselos (D_0 + 2): z cos (Delta = 1, ungerade) und z sin (Delta = 2, "
      f"gerade). Die Neumann-Bedingung phi'(0) = 0 erfüllen {neumann_raum}: Sie entfernt nur die Quelle des teilweise "
      "masselosen Felds (Delta = 1), nicht dessen Zustand mit Delta = 2. Lorentzsch bleibt mit Dirichlet + Neumann "
      "neben dem Graviton ein Geistturm mit Delta = 2. Die Neumann-Bedingung allein gibt also keine Einstein-"
      "Gravitation; es braucht eine zweite Bedingung",
      wron != 0 and neumann_raum == ["Einstein gerade", "Einstein ungerade", "PM gerade"]
      and all(tab_m[n]["löst D_0"] for n in ("Einstein gerade", "Einstein ungerade"))
      and all(tab_m[n]["löst D_0 + 2"] for n in ("PM ungerade", "PM gerade"))
      and [tab_m[n]["führende Potenz"] for n in moden] == [0, 3, 1, 2])
a_, b_ = sp.symbols("a b")
reg = (a_ + b_*z)*sp.exp(-k_*z)                        # im Inneren reguläre euklidische Lösungen
reg_ok = sp.simplify(sp.diff(reg, z, 4) - 2*k_**2*sp.diff(reg, z, 2) + k_**4*reg) == 0
bsol = sp.solve(sp.diff(reg, z).subs(z, 0), b_)[0]
ein_e = reg.subs(b_, bsol).subs(a_, 1)
ein_ok = sp.simplify(sp.diff(ein_e, z, 2) - 2*sp.diff(ein_e, z)/z - k_**2*ein_e) == 0
SW = sp.simplify(sp.integrate(sp.simplify((sp.diff(ein_e, z, 2) - k_**2*ein_e)**2), (z, 0, sp.oo)))
ee = sp.symbols("varepsilon", positive=True)
rand = sp.simplify((ein_e*sp.diff(ein_e, z)/z**2).subs(z, ee))
SEser = sp.series(-rand, ee, 0, 1).removeO()
SE_endl = sp.simplify(SEser.coeff(ee, 0))
check(f"[E/C, P_Neumann + P_reg] Maldacenas Aussage in ihrem Geltungsbereich: Euklidisch und im Inneren regulär "
      f"(P_reg) bleibt (a + b z) e^(-k z) ({reg_ok}); die Neumann-Bedingung gibt b = {bsol}, also genau die Einstein-"
      f"Mode (1 + k z) e^(-k z) ({ein_ok}). Auf ihr ist die Weyl-Wirkung alpha_W Int (phi'' - k^2 phi)^2 dz = "
      f"alpha_W ({SW}), endlich ohne Gegenterm, und die Einstein-Wirkung 2 b Int z^-2 (phi'^2 + k^2 phi^2) dz hat nach "
      f"Abzug des lokalen Terms k^2/epsilon den Wert 2 b ({SE_endl}). Gleichheit verlangt alpha_W = -b = "
      "-l^2/(64 pi G), dieselbe Relation wie aus dem Kern. Für G > 0 braucht AdS also alpha_W < 0, dS alpha_W > 0",
      reg_ok and ein_ok and sp.simplify(bsol - a_*k_) == 0 and sp.simplify(SW - 2*k_**3) == 0
      and sp.simplify(SE_endl + k_**3) == 0)

# C3 ---------------------------------------------------------------------------------------------------------------
GM = [[0, Fr(1, 2), 0, 0], [Fr(1, 2), 0, 0, 0], [0, 0, 0, -Fr(1, 2)], [0, 0, -Fr(1, 2), 0]]
GMI = [[0, 2, 0, 0], [2, 0, 0, 0], [0, 0, 0, -2], [0, 0, -2, 0]]
NZ = [(m, n) for m in range(4) for n in range(4) if GM[m][n]]


def vdot(a, b):
    return sum(GM[m][n]*a[m]*b[n] for m, n in NZ)


def vec(lam, lt):
    """Bispinor lam lt^T als Vektor (11, 22, 12, 21); vdot ist die halbe Polarisierung der Determinante."""
    return [lam[0]*lt[0], lam[1]*lt[1], lam[0]*lt[1], lam[1]*lt[0]]


def ang(a, b):
    return a[0]*b[1] - a[1]*b[0]


def kreuz(u, v):
    return [u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0]]


def cdot(u, v):
    return sum(x*y for x, y in zip(u, v))


E4 = [[Fr(int(i == j)) for j in range(4)] for i in range(4)]
E3 = [[Fr(int(i == j)) for j in range(3)] for i in range(3)]
PAARUNGEN = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]


def kinematik(rng):
    while True:
        lam = [[Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(2)] for _ in range(4)]
        lt = [[Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(2)] for _ in range(4)]
        d = ang(lam[2], lam[3])
        if d == 0:
            continue
        for b in range(2):
            r = [-(lam[0][a]*lt[0][b] + lam[1][a]*lt[1][b]) for a in range(2)]
            lt[2][b] = (r[0]*lam[3][1] - r[1]*lam[3][0])/d
            lt[3][b] = (lam[2][0]*r[1] - lam[2][1]*r[0])/d
        if all(ang(lam[i], lam[j]) != 0 and ang(lt[i], lt[j]) != 0 for i in range(4) for j in range(i)):
            return lam, lt


def v3(L1, L2, L3):
    """Kubischer Vertex aus S_3 = -g f (d_mu A_nu) A^mu A^nu, ohne den Faktor -i g."""
    legs = [L1, L2, L3]
    tot = 0
    for a, b, c in itertools.permutations(range(3)):
        f = cdot(legs[a][2], kreuz(legs[b][2], legs[c][2]))
        if f:
            tot += f*vdot(legs[a][0], legs[b][1])*vdot(legs[a][1], legs[c][1])
    return tot


def strom(L1, L2):
    K = [-(x + y) for x, y in zip(L1[0], L2[0])]
    return [[v3(L1, L2, (K, E4[r], E3[e])) for r in range(4)] for e in range(3)]


def v4(legs):
    tot = 0
    for a, b, c, d in itertools.permutations(range(4)):
        f = cdot(kreuz(legs[a][2], legs[b][2]), kreuz(legs[c][2], legs[d][2]))
        if f:
            tot += f*vdot(legs[a][1], legs[c][1])*vdot(legs[b][1], legs[d][1])
    return -Fr(1, 4)*tot


def m_ym(legs):
    """Baumamplitude als effektive Wirkung (g = 1): S_4 + Summe J.J/s; die Ströme tragen (-i)^2 = -1."""
    tot = v4(legs)
    for (a, b), (c, d) in PAARUNGEN:
        s = 2*vdot(legs[a][0], legs[b][0])
        J1, J2 = strom(legs[a], legs[b]), strom(legs[c], legs[d])
        tot += -sum(J1[e][r]*GMI[r][q]*J2[e][q] for e in range(3) for r, q in NZ)/s
    return tot


def feld(k, e):
    kl = [sum(GM[m][n]*k[n] for n in range(4)) for m in range(4)]
    el = [sum(GM[m][n]*e[n] for n in range(4)) for m in range(4)]
    return [[kl[m]*el[n] - kl[n]*el[m] for n in range(4)] for m in range(4)]      # F_mn ohne den Faktor -i


def fdot(A, B):
    return sum(A[m][n]*GMI[m][r]*GMI[n][q]*B[r][q] for m, r in NZ for n, q in NZ)


def spannung(Fa, Fb):
    ff = fdot(Fa, Fb)
    return [[sum(Fa[m][a]*GMI[a][b]*Fb[n][b] + Fb[m][a]*GMI[a][b]*Fa[n][b] for a, b in NZ) - Fr(1, 2)*GM[m][n]*ff
             for n in range(4)] for m in range(4)]


SIGMA = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]


def xmat(c):
    return sum((sp.Rational(c[a].numerator, c[a].denominator)*SIGMA[a]/2 for a in range(3)), sp.zeros(2, 2))


def rat(x):
    return sp.Rational(x.numerator, x.denominator)


rng = random.Random(148)
cym_set, q_set, eich_max, null_max, spur_T, erh_T, faelle = set(), set(), 0, 0, 0, 0, 0
for probe in range(3):
    lam, lt = kinematik(rng)
    cols = [[Fr(rng.randint(-3, 3)) for _ in range(3)] for _ in range(4)]
    kv = [vec(lam[i], lt[i]) for i in range(4)]
    assert all(sum(k[m] for k in kv) == 0 for m in range(4))
    Xm = [xmat(c) for c in cols]
    for neg in itertools.combinations(range(4), 2):
        while True:
            mus = [[Fr(rng.randint(1, 7)), Fr(rng.randint(-7, 7), 3)] for _ in range(4)]
            if all(ang(lt[i], mus[i]) != 0 and ang(mus[i], lam[i]) != 0 for i in range(4)):
                break
        pol_ = [[x/ang(lt[i], mus[i]) for x in vec(lam[i], mus[i])] if i in neg
                else [x/ang(mus[i], lam[i]) for x in vec(mus[i], lt[i])] for i in range(4)]
        legs = [(kv[i], pol_[i], cols[i]) for i in range(4)]
        M = m_ym(legs)
        eich_max = max(eich_max, abs(m_ym([(kv[0], kv[0], cols[0])] + legs[1:])))
        num = ang(lam[neg[0]], lam[neg[1]])**4
        B = 0
        for perm in itertools.permutations([1, 2, 3]):
            o = (0,) + perm
            tr = (Xm[o[0]]*Xm[o[1]]*Xm[o[2]]*Xm[o[3]]).trace()
            pt = 1
            for j in range(4):
                pt *= ang(lam[o[j]], lam[o[(j + 1) % 4]])
            B += tr*rat(num)/rat(pt)
        cym_set.add(sp.nsimplify(sp.simplify(rat(M)/B)))
        Fs = [feld(kv[i], pol_[i]) for i in range(4)]
        for (a, b), (c, d) in PAARUNGEN:
            s = 2*vdot(kv[a], kv[b])
            Tab, Tcd = spannung(Fs[a], Fs[b]), spannung(Fs[c], Fs[d])
            spur_T = max(spur_T, abs(sum(GMI[m][n]*Tab[m][n] for m, n in NZ)))
            kab = [kv[a][m] + kv[b][m] for m in range(4)]
            erh_T = max(erh_T, max(abs(sum(kab[m]*Tab[m][n] for m in range(4))) for n in range(4)))
            tt = fdot(Tab, Tcd)                        # (-i)^4 = 1: Wert des wahren T_ab . T_cd
            adt = num/(ang(lam[a], lam[b])**2*ang(lam[c], lam[d])**2)
            if (neg[0] in (a, b)) != (neg[1] in (a, b)):
                q_set.add(adt*s*s/tt)
                faelle += 1
            else:
                null_max = max(null_max, abs(tt))
cym = sorted(cym_set, key=str)[0]
Qw = sorted(q_set)[0]
NX = sp.Rational(1, 2)
alphaW_g2k = -1/(4*cym*rat(Qw)*NX**2)                  # alpha_W g^2 k
ZAHLEN["C_kopplung"] = {"M_YM/(g^2 Summe Tr(X..) <ij>^4/PT), X = sigma/2": str(cym),
                        "<ij>^4/(<ab>^2 <cd>^2) / (T_ab.T_cd/s^2)": str(Qw), "Fälle im Gravitonkanal": faelle,
                        "alpha_W g^2 k": str(alphaW_g2k)}
check(f"[E] Yang-Mills als Maßstab (exakt rational, (2,2)-Kinematik, su(2) mit X = sigma/2, 3 Punkte x 6 "
      f"Helizitätslagen, zufällige Farben und Referenzspinoren): Die Baumamplitude aus L = -(1/4) F^a F^a (Kontaktterm + "
      f"drei Austauschkanäle, als effektive Wirkung M = L_eff) ist eichinvariant (e_1 -> k_1 gibt {eich_max}) und "
      f"gleich M = ({cym}) g^2 Summe_sigma Tr(X_1 X_s2 X_s3 X_s4) <ij>^4/(<1 s2><s2 s3><s3 s4><s4 1>) in allen "
      f"{3*6} Fällen. Das legt das Vorzeichen des Einspuranteils fest: Unitäres Yang-Mills hat +g^2 vor Parke-Taylor",
      eich_max == 0 and cym_set == {sp.Integer(1)})
check(f"[E] Weyl-Austausch im selben Maßstab: Mit S_L = -alpha_W Int W^2 (W^2 -> (1/2)(Box h)^2) und der Kopplung "
      f"(1/2) h T ist M_W = T_ab.T_cd/(4 alpha_W s^2) (Maxwell-Spannungstensoren spurfrei: {spur_T}, erhalten: {erh_T}). "
      f"Auf den {faelle} Fällen des Gravitonkanals ist <ij>^4/(<ab>^2 <cd>^2) = ({Qw}) T_ab.T_cd/s_ab^2 exakt (z124: 4), "
      f"im Skalarkanal ist T.T = {null_max}. Gegenprobe des Vorzeichens: Mit Einstein-Term ist der Kern "
      "s/(32 pi G) - alpha_W s^2, die massive Nullstelle liegt für alpha_W > 0 bei s = 1/(32 pi G alpha_W) > 0 "
      "(Stelles nicht tachyonischer Geist)", q_set == {Fr(4)} and null_max == 0 and spur_T == 0 and erh_T == 0
      and faelle == 36)
check(f"[E unter T, P_N4; neu] Native Stärke: Die Stromalgebra gibt je Schleife -Tr, also Einspur -k Summe Tr/PT und "
      f"Doppelspur +k^2 Summe Tr Tr/(<ab>^2 <cd>^2) (z110, z124; das relative Vorzeichen ist i^2 aus zwei Kommutatoren "
      f"und k > 0). Setzt man den Einspuranteil gleich unitärem Yang-Mills (N_0 (-k) = g^2), so ist der Doppelspuranteil "
      f"-g^2 k Tr Tr 4 T.T/s^2 mit Tr(X_a X_b) = delta_ab/2, also 1/(4 alpha_W) = -g^2 k: alpha_W = ({alphaW_g2k})/(g^2 k). "
      "Der native Weyl-Koeffizient ist negativ: Die konforme Gravitation des Doppelspursektors hat relativ zu "
      "Yang-Mills das Vorzeichen, bei dem der statische Austausch abstößt und die euklidische Weyl-Wirkung nach unten "
      "unbeschränkt ist. Genau dieses Vorzeichen verlangt AdS für G > 0",
      alphaW_g2k == -sp.Rational(1, 4))

# C4 ---------------------------------------------------------------------------------------------------------------
G_l2 = sp.simplify(-1/(64*sp.pi*(alphaW_g2k/(g2s*ks))))            # G/l^2 = -1/(64 pi alpha_W)
G_R2 = sp.simplify(G_l2/4)
les135 = {"Netzspur je Periode, 1/g^2 = 1/2": mp.mpf(2), "Niveau mit P_norm, 1/g^2 = 4": mp.mpf(1)/4,
          "W je Markenzustand, 1/g^2 = 1/4096": mp.mpf(4096)}
tab_g = {name: gs(gv/(64*mp.pi)) for name, gv in les135.items()}
kappa_bed = sp.simplify(sp.solve(sp.Eq(1/(4096*sp.pi*kap), G_R2), kap)[0])
ZAHLEN["C_G"] = {"G/l^2": str(G_l2), "G/R_Vak^2": str(G_R2), "G/R_Vak^2 für die drei Lesarten von z135 (k = 1)": tab_g,
                 "zum Vergleich Uhrspur kappa = 1": gs(1/(4096*mp.pi)), "Uhrspur 2 c_3": gs(1/1024),
                 "Abgleich beider Wege verlangt kappa =": str(kappa_bed),
                 "Ziel COSMO.LAM.01 (z98, nur Vergleich)": "10^-122,63"}
check(f"[C, P_Neumann + P_reg oder P_pos, dazu T, P_N4, P_X] Die Zahl: Aus alpha_W = -1/(4 g^2 k) und "
      f"l^2/(64 pi G) = -alpha_W folgt G/l^2 = {G_l2} und G/R_Vak^2 = {G_R2}: G > 0 auf dem nativen AdS, und die "
      "Gravitationskopplung in Einheiten des Vakuumradius ist die Eichkopplung mal Niveau, G/R_Vak^2 = alpha_g k/16 mit "
      f"alpha_g = g^2/(4 pi). Mit k = 1 und den drei Lesarten von g aus z135: {tab_g}; die Uhrspur gibt "
      f"{gs(1/(4096*mp.pi))} (kappa = 1) und {gs(1/1024)} (2 c_3). Keine zwei Werte stimmen überein. Beide Wege "
      f"fallen genau dann zusammen, wenn kappa = {kappa_bed}, also 64 kappa = 1/(g^2 k): Das Uhrspurgewicht je Index "
      "wäre die inverse Eichkopplung. Das ist eine Konsistenzbedingung zwischen zwei freien Zahlen, keine Vorhersage. "
      "Aus '1/s^2' wird '1/s' mit festem Verhältnis zur Eichkopplung, aber nicht mit fester Zahl",
      sp.simplify(G_R2 - g2s*ks/(64*sp.pi)) == 0 and sp.simplify(kappa_bed - 1/(64*g2s*ks)) == 0
      and len(set(tab_g.values()) | {gs(1/(4096*mp.pi)), gs(1/1024)}) == 5)

# C5 ---------------------------------------------------------------------------------------------------------------
tt_, x1, x2, x3 = sp.symbols("t x1 x2 x3", real=True)
co = [tt_, x1, x2, x3]
xx = -tt_**2 + x1**2 + x2**2 + x3**2
Om = (1 - xx)/(2*l)
gk = sp.diag(-1, 1, 1, 1)/Om**2
inv = [c/xx for c in co]
Jac = sp.Matrix(4, 4, lambda i, j: sp.diff(inv[i], co[j]))
xx_inv = sp.simplify(-inv[0]**2 + inv[1]**2 + inv[2]**2 + inv[3]**2)
Om_inv = sp.simplify((1 - xx_inv)/(2*l))
iso_ok = sp.simplify(Jac.T*(sp.diag(-1, 1, 1, 1)/Om_inv**2)*Jac - gk) == sp.zeros(4, 4)
om_flip = sp.simplify(Om_inv + Om/xx) == 0
# Einstein-Bedingung des konform flachen Modells: R_mn = 2 d_m d_n Om/Om + eta_mn (Box Om/Om - 3 (d Om)^2/Om^2)
eta4 = sp.diag(-1, 1, 1, 1)
grad = [sp.diff(Om, c) for c in co]
boxO = sum(eta4[i, i]*sp.diff(Om, co[i], 2) for i in range(4))
dO2 = sum(eta4[i, i]*grad[i]**2 for i in range(4))
ricK = sp.Matrix(4, 4, lambda i, j: 2*sp.diff(Om, co[i], co[j])/Om + eta4[i, j]*(boxO/Om - 3*dO2/Om**2))
ads_ok = sp.simplify(ricK + (3/l**2)*gk) == sp.zeros(4, 4)
# Abbildungen auf R^(4,2), Koordinaten (X^-1, X^0, X^1, X^2, X^3, X^4), x^mu = X^mu/(X^-1 + X^4)
tau6 = sp.diag(1, -1, -1, -1, -1, 1)
iota6 = sp.diag(-1, -1, 1, 1, 1, 1)
T6 = sp.diag(1, -1, 1, 1, 1, 1)
rho6 = sp.diag(1, 1, 1, 1, 1, -1)
prod6 = iota6*tau6*T6
proj_gleich = prod6 == -rho6


def wirkung(Mx):
    X = sp.Matrix([(1 + xx)/2, tt_, x1, x2, x3, (1 - xx)/2])
    Y = Mx*X
    return [sp.simplify(Y[i]/(Y[0] + Y[5])) for i in range(1, 5)]


rho_x = wirkung(rho6)
iota_x = wirkung(iota6)
tau_x = wirkung(tau6)
geo_ok = all(sp.simplify(rho_x[i] - co[i]/xx) == 0 for i in range(4)) \
    and sp.simplify(iota_x[0] - tt_/xx) == 0 and all(sp.simplify(iota_x[i] + co[i]/xx) == 0 for i in (1, 2, 3)) \
    and all(sp.simplify(tau_x[i] + co[i]) == 0 for i in range(4))
dets = [tau6.det(), iota6.det(), rho6.det()]
check(f"[E] Geometrie der Randspiegelung: I_AdS = I_pi - I_omega gibt den Faktor Omega = (1 - x^2)/(2 l) (z98, z112); "
      f"eta/Omega^2 ist AdS4 mit R_mn = -(3/l^2) g_mn ({ads_ok}), Rand x^2 = 1. Die Inversion x -> x/x^2 lässt den Rand "
      f"punktweise fest, ist eine Isometrie ({iso_ok}) und kehrt Omega um (Omega -> -Omega/x^2: {om_flip}): Sie ist "
      "die Spiegelung am konformen Rand, im Poincaré-Bild z -> -z, auf dem AdS-Hyperboloid die Antipodenabbildung. "
      "Auf R^(4,2) ist sie rho = diag(1,1,1,1,1,-1) mit Determinante "
      f"{dets[2]}, die Punktspiegelung tau (z84/z123) und die native Halbdrehung iota = e^(pi E) (z127) haben "
      f"Determinante {dets[0]} und {dets[1]}; es gilt rho = iota tau T projektiv ({proj_gleich}) mit der Zeitspiegelung "
      f"T (Wirkungen auf x nachgerechnet: {geo_ok}). Die Randspiegelung ist uneigentlich: Sie liegt nicht in der von "
      "tau, iota und so(3,2) erzeugten Gruppe und ist nativ nur antiunitär da, als iota J mit J = theta tau "
      "(theta = radiale Inversion, z84)", ads_ok and iso_ok and om_flip and proj_gleich and geo_ok
      and dets == [1, 1, -1])
par = {n: t["Parität unter z -> -z"] for n, t in tab_m.items()}
check(f"[E, negativ] Erzwingt ein natives Datum P_Neumann? Paritäten der vier Moden unter der Randspiegelung: {par}. "
      "Jeder Paritätssektor enthält eine Einstein- und eine teilweise masselose Mode; die Neumann-Bedingung ist "
      "keine Eigenbedingung von rho. tau und iota sind Isometrien von AdS4 und bilden jeden der beiden Lösungsräume "
      "(Kern von D_0, Kern von D_0 + 2) auf sich ab; sie stellen keine Randbedingung. J ist antiunitär und "
      "tauscht Algebra und Kommutante; aus J Omega = Omega folgt nur <d_n h> = 0 auf der Spiegelfläche, nicht "
      "d_n h = 0 als Zwangsbedingung (Schwankungen spiegelungsungerader Größen verschwinden in einem "
      "spiegelungsinvarianten Zustand nicht). Die Kappe gehört zum dS-Büschel I_pi + I_omega und wird von keiner "
      "Naht getragen (z128, z136). Der konforme Rand ist für das 4D-konforme Feld eine gewöhnliche Fläche (alle vier "
      "Moden analytisch in z). P_Neumann ist durch kein natives Datum erzwungen",
      sorted(par.values()) == ["gerade", "gerade", "ungerade", "ungerade"]
      and par["Einstein gerade"] != par["Einstein ungerade"] and par["PM gerade"] != par["PM ungerade"])
check("[C, neu] Schwächere, nativ formulierbare Fassung: P_pos = 'der ausgetauschte Spin-2-Sektor ist eine unitäre "
      "Darstellung positiver Energie der nativen so(3,2)'. Nach A6 sind die teilweise masselosen Moden (tiefste "
      "Gewichte E0 = 1, 2 bei Spin 2) dann ausgeschlossen, es bleibt D(3,2); das ersetzt P_Neumann und P_reg "
      "zusammen und braucht nur Daten, die nativ sind (E >= 0, Hilbertraum, so(3,2) = Der(I_AdS)). Die Norm der "
      "Einstein-Moden in der Weyl-Form ist die Einstein-Norm mit b = -alpha_W (C1), positiv genau für alpha_W < 0: Das "
      "native Vorzeichen alpha_W = -1/(4 g^2 k) und das native AdS-Vorzeichen (z89, z104, z118) passen zusammen; auf "
      "dS wäre dieselbe Kopplung ein Geist. Erzwungen ist P_pos nicht: Nach A6 enthält der native Fockraum D(3,2) "
      "nicht, der Sektor muss also hinzutreten (P_X aus z123), und im flachen Grenzfall wird die Norm null",
      alphaW_g2k < 0 and vielf == [3, 5, 7] and d32 == {1: 0, 2: 0, 3: 0})
befund("Einstein-Pol aus der konformen Gravitation", "bedingt (P_Neumann + P_reg oder P_pos; T, P_N4, P_X)",
       "Teil C: Auf AdS4 zerfällt der Doppelpol in Graviton und teilweise masseloses Feld. Mit Neumann-Bedingung und "
       "Regularität (oder mit P_pos) bleibt Einstein-Gravitation mit G/R_Vak^2 = g^2 k/(64 pi) > 0; die Zahl hängt an "
       "der freien Eichkopplung. Die Neumann-Bedingung allein genügt lorentzsch nicht, und kein natives Datum "
       "erzwingt sie.")

# ====================================================================================================================
say("\nTeil D: Bilanz und Urteil")
ol_pfad = REPO/"research/originlab/reconstruction_search.py"
klausel = ol_pfad.exists() and '("graviton", st["gravity"] == "supported")' in ol_pfad.read_text(encoding="utf-8")
check(f"[E] OriginLab (nur gelesen): Die Klausel 'graviton' verlangt den Status 'supported' der Gravitationsachse "
      f"({klausel}). Nach z148 bleibt der Wert bedingt: Beide Wege zum Pol brauchen eine Prämisse, die den Sektor "
      "D(3,2) erst hinzufügt (P_dyn oder P_X), und die Kopplung bleibt frei. Die Klausel bleibt leer", klausel)
bilanz = [
    ("E", "quadratischer Operator der Uhrspurwirkung",
     "L2 = b eps f ^ f, erste Ordnung, 40 Felder, 10 Eichfunktionen; tau Hilfsfeld (Rang 24); nach Elimination exakt "
     "linearisierte Einstein-Gleichung mit Lambda = -3/l^2, b = l^2/(64 pi G) = 16 kappa."),
    ("E", "Freiheitsgrade und Helizitäten",
     "Jetzählung 9 + 4 m, Zuwachs 4: zwei Freiheitsgrade, ein Dublett der Ladung 2 (Helizität +-2), D(3,2), "
     "Delta = 3."),
    ("C (P_dyn)", "Residuum",
     "32 pi G = l^2/(32 kappa), positiv bei richtigem Paarungsbit; G/R_Vak^2 = 1/(4096 pi kappa)."),
    ("E negativ", "Bestand ohne P_dyn",
     "Hesse-Form null auf allen nativen (flachen) Richtungen; der native Fockraum enthält kein D(3,2), Spin 2 beginnt "
     "als Kompositum bei E0 = 4."),
    ("E", "induzierter Propagator",
     "Pi = (c/(32 pi^2)) p^4 ln(p^2/mu^2); ohne Baumterm kein masseloser Pol, Tachyon bei Lambda_W, Potential ~ r; mit "
     "Baumterm masseloser Pol plus komplexes Paar oder zwei Tachyonen (Summenregel), nicht unitär."),
    ("O", "Skala des Logarithmus", "Lambda_W l ist eine freie reine Zahl (P_cut)."),
    ("E", "konforme Gravitation auf AdS4",
     "z^4 Box^2 = D_0 (D_0 + 2); masseloser Anteil = Einstein mit l^2/(64 pi G) = -alpha_W; teilweise masseloser "
     "Anteil (Delta = 1, 2) mit umgekehrtem Vorzeichen."),
    ("E", "Neumann-Bedingung",
     "entfernt nur Delta = 1; lorentzsch bleibt der Geistturm Delta = 2; mit Regularität (euklidisch) genau die "
     "Einstein-Mode."),
    ("E unter T, P_N4", "native Kopplung", "alpha_W = -1/(4 g^2 k): negativ, passend zu AdS."),
    ("C", "Zahl", "G/R_Vak^2 = g^2 k/(64 pi) = alpha_g k/16; Abgleich mit der Uhrspur verlangt 64 kappa = 1/(g^2 k)."),
    ("E negativ / C", "Nativität der Randbedingung",
     "kein natives Datum erzwingt P_Neumann; ersetzbar durch P_pos (unitäre so(3,2)-Darstellung positiver Energie), "
     "das seinerseits P_X voraussetzt."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
urteil = ("Der native Gravitationssektor von H6 hat keinen masselosen Spin-2-Pol; er hat zwei bedingte Wege dorthin, "
          "und für beide ist die fehlende Prämisse jetzt genau benannt. Nichts ist befördert. (1) Uhrspurwirkung: Um "
          "AdS4 ist die MacDowell-Mansouri-Wirkung exakt quadratisch in f = D tau + (h e - e h)/l^2. Der "
          "Spinzusammenhang ist ein Hilfsfeld, zehn Eichfunktionen fallen heraus, und es bleiben genau zwei "
          "Freiheitsgrade mit Helizität +-2, die masselose Darstellung D(3,2), mit Residuum 32 pi G = l^2/(32 kappa). "
          "Das Vorzeichen ist das Paarungsbit. Das alles gilt nur, wenn das Vierbein dynamisch ist (P_dyn). Ohne "
          "P_dyn verschwindet die Hesse-Form auf allen nativen Richtungen, weil jede native Familie flach ist, und "
          "der native Fockraum enthält unter so(3,2) kein D(3,2); der Spin-2-Primär aus z137 ist das massive "
          "Kompositum D(4,2). (2) Induzierter Propagator: Aus dem nativen <TT> folgt ohne Baumterm kein masseloser "
          "Pol, sondern ein Tachyon an der Logarithmusskala und ein lineares Potential. Mit dem Baumterm gibt es "
          "den masselosen Pol, aber die Summenregel erzwingt zusätzlich ein komplexes Polpaar oder zwei Tachyonen. "
          "Die Skala Lambda_W ist eine Kontaktterm-Freiheit und nicht nativ. (3) Konforme Gravitation: Auf AdS4 "
          "zerfällt der Doppelpol 1/s^2 in ein masseloses Graviton und ein teilweise masseloses Spin-2-Feld mit "
          "umgekehrtem Vorzeichen. Maldacenas Neumann-Bedingung entfernt davon nur das Randverhalten Delta = 1; "
          "lorentzsch bleibt ein Geistturm mit Delta = 2, und erst zusammen mit der Regularität im Inneren (P_reg) "
          "bleibt genau die Einstein-Mode. Der native Koeffizient ist alpha_W = -1/(4 g^2 k), negativ, also genau das "
          "Vorzeichen, das AdS für G > 0 braucht. Die Zahl ist G/R_Vak^2 = g^2 k/(64 pi): Aus 1/s^2 wird 1/s mit "
          "festem Verhältnis zur Eichkopplung, nicht mit fester Zahl; mit den drei Lesarten von g aus z135 liegt sie "
          "zwischen 0,0012 und 20, wieder ein Vakuum von Planck-Größe. Kein natives Datum erzwingt die "
          "Randbedingung: Die Randspiegelung ist uneigentlich und nativ nur antiunitär vorhanden, ihre "
          "Paritätssektoren mischen Einstein- und Geistmoden, tau und iota sind Isometrien. (4) Für K13 und die "
          "Klausel 'Graviton': Es fehlt in kleinster Form eine einzige Aussage, nämlich dass ein Sektor D(3,2) der "
          "nativen so(3,2) existiert und mit positiver Norm an die Nähte koppelt. Das sind zwei reelle Funktionen "
          "(ein Dublett der Helizität +-2), nicht ein ganzes Vierbein. P_dyn liefert sie über die Uhrspurwirkung, "
          "P_X + P_pos über den Doppelspursektor. P_curv ist damit nicht durch eine native Aussage ersetzbar, wohl "
          "aber durch eine schwächere: P_h2 = 'es gibt eine dynamische spurfreie transversale Richtung des Vierbeins "
          "mit f != 0', die linearisierte Form von P_curv, oder gleichwertig P_X + P_pos. Neu ist eine Bindung: "
          "Gelten beide Wege, so ist 64 kappa = 1/(g^2 k), und P_G fällt mit der offenen Eichkopplung von K15 "
          "zusammen. K13 bleibt L, die Klausel 'Graviton' bleibt leer.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z148_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z148.txt").write_text("\n".join(OUT) + "\n")
