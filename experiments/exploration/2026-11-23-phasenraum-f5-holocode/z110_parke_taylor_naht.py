"""Z110: Wechselwirkung als Antwort der Naht -- Stromkorrelator = Parke-Taylor (MHV)? (Tiefenwelle T1, Hypothese T).

Leitlinie 'Rekonstruktion statt Auswahl': Die Korrelatoren werden aus den nativen Nahtströmen gerechnet (Fockraum im
Quarterrahmen, Wick, Wardidentität), nicht aus einer gewünschten Amplitude. Parke-Taylor-Faktor, Nairs Stromalgebra auf CP^1
(1988), Twistorstring (Witten 2003; Berkovits-Witten 2004: Mehrspurterme = konforme Supergravitation) und die DDM-Zerlegung
sind Vergleichsziele, keine Prämissen. Firewall DIMENSION.UPLIFT.FIREWALL.01: Die Naht ist eine chirale 2D-Struktur; 4D tritt
nur als Lesart T auf (Naht = Twistorgerade L_x, Nahtlagen = Raumzeitpunkte, z61/z85), ausgewertet als Amplitude, nicht als
lokales Feld. Gelöste Blöcke (E8, 16, Flavour, alpha) werden benutzt, nicht geöffnet; nichts wird befördert. Statuskürzel:
E exakt, C bedingt, L Lesart, O offen.

Bestand: z58 (c(x) = 1/(1 + e^{4 pi x}), x in Z + alpha, alpha = 1/4, 3/4, c = 1 je Quarter, Twistgewicht (alpha - 1/2)^2/2),
z38 (Flussladung q = theta.beta, theta = (1/2)^8: E8 = 1 + 56 + 126 + 56 + 1, Zentralisator e7 + u(1)), z75 (<J J J> = A f,
verbundener Vierpunkt != 0; Quarterrahmen der 16 Majorana: alle acht komplexen Fermionen bei alpha = 1/4, U(8)), z89 (I_pi bis auf die Modulphase nativ), z93 (Bosonisierung bei beta = 4 pi: G = Twist x Theta x Oszillator), z83/z97
(Schrittart S_a: CCR-Spinor), z99 (zwei Nähte, Gitter Z^2, su(2)_1 nach Frenkel-Kac, gleiche und verschiedene Quarter).

Teil 1  Native Stromkorrelatoren auf CP^1 im Vakuum. (a) Fockraum zweier Nähte (gleiche Quarter 1/4, 1/4 und verschiedene
        1/4, 3/4): Modenkorrelatoren der su(2)-Ströme gegen die Laurentkoeffizienten der Wickformel mit den nativen
        Propagatoren (z/w)^(1/2 - alpha)/(z - w); der Zwist kürzt sich in jeder Schleife gleicher Quarter, nicht in gemischten.
        (b) n = 3..6 exakt (rational): Wardrekursion mit dem Niveau k als Symbol gegen Fermionschleifen (Dirac: su(2) aus z99,
        su(8) c u(8) des Originals; Majorana: so(16)_1); Einspuranteil = zyklische Parke-Taylor-Summe mit Koeffizient -k,
        DDM-Form, Photonentkopplung; Mehrspuranteile k^c; das ungetwistete Paarungs-su(2) des Originalpaars. (c) E8 im
        Quarterrahmen: Zwist v = (1/4)^8 = theta/2 auf den Wurzeln, also (-1)^q mit der Flussladung aus z38; Fixalgebra,
        getwistete Ströme, Grundzustände von E8 + v.
Teil 2  KMS-Zustand bei beta = 4 pi: thermische Propagatoren aus c(x) (Modensummen gegen die Thetaform aus z93), Abweichung
        des Einspuranteils von Parke-Taylor, erste Ordnung c(1/4) = 1/(1 + e^pi), Exponent aus der beta-Abhängigkeit,
        Gegenprobe mit Projektion auf Gesamtladung 0.
Teil 3  4D-Lesart auf Twistorgeraden: z_i -> lambda_i, projektive Gewichte, SL(2,C); Helizitätszähler <ij>^4: Grassmann-
        Integral über 2N Nullmoden gibt <ij>^N und Helizität 1 - N/2; native Kandidaten (Nahtfermionen, S_a, Quarterpaar,
        eta-Dublett), Lorentzinvarianz der Paarung.
Teil 4  Kopplung (Niveau, KMS, Uhr, Grad), Bilanz und Urteil zu Hypothese T; Folgen für K15 und T3/T5.
"""
import functools
import itertools
import json
import time
from collections import Counter
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.sparse as sps
import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
rng = np.random.default_rng(110)
mp.mp.dps = 30
PI = mp.pi
BETA = 4*PI
TOL = 1e-12


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def r(x, n=12):
    return float(np.round(float(x), n))


def gs(x, n=6):
    return float(f"{float(x):.{n}g}")


# ------------------------------------------------------------------ Bausteine: Permutationen, Schleifen, Parke-Taylor
def zyklen(perm):
    seen, out = set(), []
    for s in range(len(perm)):
        if s in seen:
            continue
        c, j = [], s
        while j not in seen:
            seen.add(j)
            c.append(j)
            j = perm[j]
        out.append(c)
    return out


def derangements(n):
    for perm in itertools.permutations(range(n)):
        if all(perm[i] != i for i in range(n)):
            yield perm


def glieder(c):
    return [(c[l], c[(l + 1) % len(c)]) for l in range(len(c))]


def pt(zs, ordnung):
    """Parke-Taylor-Kette 1/(z_o1o2 z_o2o3 ... z_ono1)."""
    v = Fr(1)
    for a, b in glieder(list(ordnung)):
        v /= (zs[a] - zs[b])
    return v


def produkt(Ms):
    return functools.reduce(lambda A, B: A @ B, Ms)


def spur(M):
    return int(np.trace(M))


FORMEN = {"Dirac": lambda A, B: spur(A @ B), "Majorana": lambda A, B: Fr(spur(A @ B), 2)}
SCHLEIFE = {"Dirac": lambda W: Fr(spur(W)), "Majorana": lambda W: Fr(spur(W), 2)}


def ward(ops, form):
    """<J^X1(z1) ... J^Xn(zn)> im SL(2)-invarianten Vakuum allein aus der Wardidentität (OPE J^X(z) J^Y(w) ~
    k form(X, Y)/(z - w)^2 + J^[X,Y](w)/(z - w), <J> = 0) als Polynom in k: {Potenz: Fraction}."""
    if not ops:
        return {0: Fr(1)}
    if len(ops) == 1:
        return {}
    (X0, z0), rest = ops[0], ops[1:]
    out = {}
    for i, (Xi, zi) in enumerate(rest):
        d = z0 - zi
        f = form(X0, Xi)
        if f:
            for p, c in ward(rest[:i] + rest[i + 1:], form).items():
                out[p + 1] = out.get(p + 1, 0) + c*f/(d*d)
        C = X0 @ Xi - Xi @ X0
        if C.any():
            for p, c in ward(rest[:i] + [(C, zi)] + rest[i + 1:], form).items():
                out[p] = out.get(p, 0) + c/d
    return {p: c for p, c in out.items() if c != 0}


def wick(ops, schleife):
    """Fermionschleifen: Summe über fixpunktfreie Permutationen (psi(z_i) -> psi*(z_perm(i))), je Schleife -k schleife(Wort)
    mal Parke-Taylor-Kette; Rückgabe {Zahl der Schleifen: Fraction} und die Zahl der Terme je Schleifenzahl."""
    n = len(ops)
    zs = [z for _, z in ops]
    out, anzahl = {}, Counter()
    for perm in derangements(n):
        cs = zyklen(perm)
        val = Fr(1)
        for c in cs:
            val *= -schleife(produkt([ops[i][0] for i in c]))*pt(zs, c)
        out[len(cs)] = out.get(len(cs), 0) + val
        anzahl[len(cs)] += 1
    return {p: c for p, c in out.items() if c != 0}, dict(anzahl)


def einspur_pt(ops, schleife):
    """-sum_{sigma in S_(n-1)} schleife(X_1 X_sigma2 ... X_sigman) PT(1, sigma): zyklische Parke-Taylor-Summe."""
    n = len(ops)
    zs = [z for _, z in ops]
    return sum((-schleife(produkt([ops[i][0] for i in (0,) + s]))*pt(zs, (0,) + s)
                for s in itertools.permutations(range(1, n))), Fr(0))


def ddm(ops, schleife):
    """-sum_{sigma in S_(n-2)} schleife([[..[X_1, X_s2], ..], X_s(n-1)] X_n) PT(1, sigma, n): nur Kommutatoren (f-Kette)."""
    n = len(ops)
    zs = [z for _, z in ops]
    tot = Fr(0)
    for s in itertools.permutations(range(1, n - 1)):
        N = ops[0][0]
        for i in s:
            N = N @ ops[i][0] - ops[i][0] @ N
        tot += -schleife(N @ ops[n - 1][0])*pt(zs, (0,) + s + (n - 1,))
    return tot


def punkte(n):
    zs = set()
    while len(zs) < n:
        zs.add(Fr(int(rng.integers(-40, 41)), int(rng.integers(1, 9))))
    out = list(zs)
    rng.shuffle(out)
    return out


EP = np.array([[0, 1], [0, 0]], dtype=np.int64)
EM = np.array([[0, 0], [1, 0]], dtype=np.int64)
HH = np.array([[1, 0], [0, -1]], dtype=np.int64)


def su2_zufall():
    a, b, c = (int(x) for x in rng.integers(-9, 10, size=3))
    return a*EP + b*EM + c*HH


def gl8_spurfrei():
    X = rng.integers(-2, 3, size=(8, 8)).astype(np.int64)
    X[7, 7] -= np.trace(X)
    return X


def so16_zufall():
    A = rng.integers(-1, 2, size=(16, 16))*(rng.random((16, 16)) < 0.3)
    return (np.triu(A, 1) - np.triu(A, 1).T).astype(np.int64)


# ------------------------------------------------------------------ Bausteine: Nähte im Fockraum (wie z99)
def cdag_matrix(j, nb):
    """Jordan-Wigner-Erzeuger der Mode j auf nb Moden (Vorzeichen = Parität der besetzten Moden unter j)."""
    dim = 1 << nb
    s = np.arange(dim)
    src = s[((s >> j) & 1) == 0]
    tief = src & ((1 << j) - 1)
    par = np.zeros(len(src), dtype=np.int64)
    for b in range(j):
        par ^= (tief >> b) & 1
    return sps.csr_matrix(((1 - 2*par).astype(float), (src | (1 << j), src)), shape=(dim, dim))


class Naehte:
    """Nähte zu je 2K Moden x = alpha + n, n = -K..K-1, eine gemeinsame CAR, Vakuum = Fermisee x < 0."""

    def __init__(self, alphas, K):
        self.alphas = [Fr(a) for a in alphas]
        self.N = len(alphas)
        self.modes = [(i, a + n) for i, a in enumerate(self.alphas) for n in range(-K, K)]
        self.NB = len(self.modes)
        self.dim = 1 << self.NB
        self.idx = {md: j for j, md in enumerate(self.modes)}
        self.cd = [cdag_matrix(j, self.NB) for j in range(self.NB)]
        self.c = [m.T.tocsr() for m in self.cd]
        s = np.arange(self.dim)
        occ = (s[:, None] >> np.arange(self.NB)[None, :]) & 1
        unten = np.array([x < 0 for _, x in self.modes])
        naht = np.array([i for i, _ in self.modes])
        self.Q = [occ[:, (naht == i) & ~unten].sum(1) - (1 - occ[:, (naht == i) & unten]).sum(1) for i in range(self.N)]
        vi = int(sum(1 << j for j in range(self.NB) if unten[j]))
        self.vac = np.zeros(self.dim)
        self.vac[vi] = 1.0
        self.cache = {}

    def strom(self, X, n):
        """J^X_n = sum_ab X_ab sum_x psi*_(a, x-n) psi_(b, x) (Fock-Normalordnung, für a = b und n = 0 die Ladung Q_a)."""
        key = (X.tobytes(), Fr(n))
        if key in self.cache:
            return self.cache[key]
        M = sps.csr_matrix((self.dim, self.dim))
        for a in range(self.N):
            for b in range(self.N):
                if X[a, b] == 0:
                    continue
                if a == b and n == 0:
                    M = M + float(X[a, a])*sps.diags(self.Q[a].astype(float))
                    continue
                for (i, x), j in self.idx.items():
                    t = self.idx.get((a, x - n)) if i == b else None
                    if t is not None:
                        M = M + float(X[a, b])*(self.cd[t] @ self.c[j])
        self.cache[key] = M.tocsr()
        return self.cache[key]

    def korr(self, liste):
        v = self.vac
        for X, m in reversed(liste):
            v = self.strom(X, m) @ v
        return float(self.vac @ v)


def wick_terme(Xs, alphas):
    """Wickformel der Fock-normalgeordneten Ströme mit den nativen Propagatoren G_a(z, w) = (z/w)^(1/2 - alpha_a)/(z - w):
    Liste (Koeffizient, Zwistexponenten je Punkt, Glieder)."""
    n = len(Xs)
    terme = []
    for perm in derangements(n):
        cs = zyklen(perm)
        links = tuple(l for c in cs for l in glieder(c))
        zt = []
        for c in cs:
            L = len(c)
            ct = []
            for fl in itertools.product(range(len(alphas)), repeat=L):
                coef = -1
                for l in range(L):
                    coef *= int(Xs[c[l]][fl[l - 1], fl[l]])
                if coef == 0:
                    continue
                ex = [Fr(0)]*n
                for l in range(L):
                    a = Fr(1, 2) - alphas[fl[l]]
                    ex[c[l]] += a
                    ex[c[(l + 1) % L]] -= a
                ct.append((coef, ex))
            zt.append(ct)
        for combo in itertools.product(*zt):
            coef, ex = 1, [Fr(0)]*n
            for cf, e in combo:
                coef *= cf
                ex = [p + q for p, q in zip(ex, e)]
            terme.append((coef, tuple(ex), links))
    return terme


@functools.lru_cache(maxsize=None)
def laurent(links, ziel, B=6):
    """Koeffizient von prod z_v^ziel_v in prod_(i,j) 1/(z_i - z_j), entwickelt für |z_0| > |z_1| > ... (radiale Ordnung)."""
    n = len(ziel)
    tot = 0
    for ks in itertools.product(range(B + 1), repeat=len(links)):
        e = [0]*n
        s = 1
        for (i, j), k in zip(links, ks):
            if i < j:
                e[i] -= 1 + k
                e[j] += k
            else:
                e[j] -= 1 + k
                e[i] += k
                s = -s
        if tuple(e) == ziel:
            tot += s
    return tot


def laurent_korr(liste, alphas):
    Xs = [X for X, _ in liste]
    ms = [Fr(m) for _, m in liste]
    tot = 0
    for coef, ex, links in wick_terme(Xs, alphas):
        ziel = [-m - 1 - e for m, e in zip(ms, ex)]
        if all(z.denominator == 1 for z in ziel):
            tot += coef*laurent(links, tuple(int(z) for z in ziel))
    return tot


# =================================================================== Teil 1
say("Teil 1: Native Stromkorrelatoren auf CP^1 im Vakuum")
say("  (a) Fockraum zweier Nähte gegen die Wickformel in Ortslage")
S = Naehte([Fr(1, 4), Fr(1, 4)], 4)                                                 # gleiche Quarter (z99)
M = Naehte([Fr(1, 4), Fr(3, 4)], 4)                                                 # verschiedene Quarter (Originalpaar)
say(f"   Fockraum 2 x 8 Moden, Dimension {S.dim}; Aufbau {time.time() - T0:.1f} s")
T3 = [[(EP, 1), (HH, 1), (EM, -2)], [(HH, 2), (EP, -1), (EM, -1)], [(EP, 1), (EM, 1), (HH, -2)], [(EM, 2), (HH, -1), (EP, -1)],
      [(EP, 2), (EM, -1), (HH, -1)], [(HH, 1), (HH, 1), (HH, -2)]]
T4 = [[(EP, 1), (EM, 1), (EP, -1), (EM, -1)], [(EP, 1), (EP, 1), (EM, -1), (EM, -1)], [(HH, 1), (EP, 1), (EM, -1), (HH, -1)],
      [(EP, 2), (HH, -1), (EM, 1), (HH, -2)], [(HH, 1), (HH, 1), (HH, -1), (HH, -1)], [(EP, 2), (EM, -1), (HH, 1), (HH, -2)]]
h = Fr(1, 2)
TM = [[(EP, h), (EM, -h)], [(EM, h), (EP, -h)], [(EP, 3*h), (EM, -3*h)], [(EM, 3*h), (EP, -3*h)], [(HH, 1), (HH, -1)],
      [(EM, h), (EP, h), (HH, -1)], [(HH, 1), (EM, -h), (EP, -h)], [(EM, 3*h), (HH, -1), (EP, -h)], [(EP, 3*h), (EM, -h), (HH, -1)],
      [(HH, 2), (EM, -h), (EP, -3*h)], [(EP, h), (EM, -h), (EP, h), (EM, -h)], [(EM, h), (EP, h), (EM, -h), (EP, -h)]]
NAME = {EP.tobytes(): "E+", EM.tobytes(): "E-", HH.tobytes(): "H"}
fock_tab, fock_ok = {}, True
for sysname, sys_, tupel in (("gleich", S, T3 + T4), ("verschieden", M, TM)):
    for liste in tupel:
        fv = sys_.korr(liste)
        lv = laurent_korr(liste, sys_.alphas)
        fock_ok &= abs(fv - lv) < TOL
        fock_tab[f"{sysname}: " + " ".join(f"{NAME[X.tobytes()]}_{m}" for X, m in liste)] = [r(fv), int(lv)]
zwist_gleich = {ex for liste in T3 + T4 for _, ex, _ in wick_terme([X for X, _ in liste], S.alphas)}
zwist_versch = {ex for liste in TM for _, ex, _ in wick_terme([X for X, _ in liste], M.alphas)}
ZAHLEN["fock_gegen_ortslage"] = fock_tab
zfmt = lambda e: "(" + ", ".join(str(x) for x in e) + ")"
ZAHLEN["zwistexponenten"] = {"gleiche Quarter": sorted({zfmt(e) for e in zwist_gleich}),
                             "verschiedene Quarter (Auswahl)": sorted({zfmt(e) for e in zwist_versch})[:8]}
nz_versch = [e for e in zwist_versch if any(x != 0 for x in e)]
gleich_null = all(x == 0 for e in zwist_gleich for x in e)
check("[E] (a) Native Ströme in Ortslage: Die Modenkorrelatoren <0|J_m1 ... J_mn|0> der su(2)-Ströme E^+-, H zweier Nähte "
      f"(Fockraum 2 x 8 Moden, {len(T3 + T4)} Tupel n = 3, 4 bei gleichen, {len(TM)} Tupel n = 2..4 bei verschiedenen Quartern) "
      "sind exakt die Laurentkoeffizienten (|z_1| > ... > |z_n|) der Wickformel mit dem nativen Propagator "
      "<psi(z) psi*(w)> = sum_(x > 0) z^(-x-1/2) w^(x-1/2) = (z/w)^(1/2 - alpha)/(z - w) (Fock-Normalordnung: nur fixpunktfreie "
      "Permutationen). Bei gleichen Quartern kürzt sich der Zwist in jeder Schleife (alle Zwistexponenten "
      f"{sorted({zfmt(e) for e in zwist_gleich})}): die Korrelatoren sind rationale Funktionen aus 1/(z_i - z_j), unabhängig "
      "von alpha. Bei verschiedenen Quartern bleiben Faktoren (z_i/z_j)^(+-1/2) "
      f"({len(nz_versch)} verschiedene Zwistmuster), E^+- haben halbzahlige Moden und <E^+_(1/2) E^-_(-1/2)> = 0 neben "
      "<E^-_(1/2) E^+_(-1/2)> = 1: ein Z_2-getwistetes su(2) mit Verzweigung an 0 und unendlich", fock_ok and
      gleich_null and len(nz_versch) > 0
      and abs(fock_tab["verschieden: E+_1/2 E-_-1/2"][0]) < TOL and abs(fock_tab["verschieden: E-_1/2 E+_-1/2"][0] - 1) < TOL)


def paar_plus(sys_, n):
    """P^+_n = sum_x psi*_(0,x) psi*_(1,-n-x): Paarstrom psi*_0 psi*_1 des Originalpaars (ganzzahlige Moden für 1/4 + 3/4)."""
    Mx = sps.csr_matrix((sys_.dim, sys_.dim))
    for (i, x), j in sys_.idx.items():
        t = sys_.idx.get((1, -n - x)) if i == 0 else None
        if t is not None:
            Mx = Mx + sys_.cd[j] @ sys_.cd[t]
    return Mx.tocsr()


Pp = {n: paar_plus(M, n) for n in (-2, -1, 0, 1, 2)}
Pm = {n: Pp[-n].T.tocsr() for n in Pp}
Qs = sps.diags((M.Q[0] + M.Q[1]).astype(float)).tocsr()
paar_zw = {m: (float(M.vac @ (Pp[m] @ (Pm[-m] @ M.vac))), float(M.vac @ (Pm[m] @ (Pp[-m] @ M.vac)))) for m in (1, 2)}
e_zw = {m: (M.korr([(EP, m - h), (EM, -m + h)]), M.korr([(EM, m - h), (EP, -m + h)])) for m in (1, 2)}
testv = [M.vac, Pp[-1] @ M.vac, Pm[-1] @ M.vac, M.strom(EP, -h) @ M.vac, Pp[-2] @ M.vac]
komm = max(np.abs(Pp[m_] @ (Pm[-m_] @ v) - Pm[-m_] @ (Pp[m_] @ v) - Qs @ v - m_*v).max()
           for m_ in (-2, -1, 0, 1, 2) for v in testv)
ZAHLEN["paarungs_su2_originalpaar"] = {"<P+_m P-_-m>, <P-_m P+_-m> (m = 1, 2)": paar_zw,
                                      "<E+_(m-1/2) E-_(-m+1/2)>, <E-_(m-1/2) E+_(-m+1/2)>": e_zw,
                                      "Defekt [P+_m, P-_-m] = Q_0 + Q_1 + m (m = -2..2, fünf Zustände)": r(komm, 15)}
check("[E] (a) Originalpaar (1/4, 3/4): Der Paarstrom P^+ = psi*_0 psi*_1, P^- = (P^+)^+ und Q_0 + Q_1 haben ganzzahlige Moden "
      f"und erfüllen [P^+_m, P^-_-m] = Q_0 + Q_1 + m (Defekt {komm:.0e}); ihre Zweipunktantworten sind symmetrisch "
      f"{ {m: v for m, v in paar_zw.items()} } = m, ohne Hintergrundladung ((1/4 - 1/2) + (3/4 - 1/2) = 0): das ungetwistete "
      f"su(2)_1 des Originals (nach Konjugation der 3/4-Naht das E^+- gleicher Quarter). Der Transferstrom E^+- ist "
      f"getwistet: <E^+_mu E^-_-mu>, <E^-_mu E^+_-mu> = {e_zw} = mu -+ 1/2 (mu = 1/2, 3/2)", komm < TOL and all(abs(a_ - m) < TOL and abs(b_ - m) < TOL for m, (a_, b_) in paar_zw.items())
      and all(abs(e_zw[m][0] - (m - 1)) < TOL and abs(e_zw[m][1] - m) < TOL for m in (1, 2)))
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

say("  (b) n = 3..6 exakt: Wardrekursion gegen Fermionschleifen gegen Parke-Taylor")
FAELLE = [("su(2)_1 zweier Nähte (z99, Dirac)", "Dirac", su2_zufall),
          ("su(8)_1 c u(8) des Originals (Quarterrahmen z75, Dirac)", "Dirac", gl8_spurfrei),
          ("so(16)_1 der 16 Majorana (Vakuum, z75)", "Majorana", so16_zufall)]
tab_b, ok_ward, ok_pt, ok_ddm, ok_k = {}, True, True, True, True
WARD = {}
STICH = 3
for name_, art, gen in FAELLE:
    for n in range(3, 7):
        for stich in range(STICH):
            for versuch in range(20):                                                # zufällige Entartung meiden
                zs = punkte(n)
                ops = [(gen(), z) for z in zs]
                W = ward(ops, FORMEN[art])
                if set(W) == set(range(1, n//2 + 1)):
                    break
            K, anzahl = wick(ops, SCHLEIFE[art])
            P1 = einspur_pt(ops, SCHLEIFE[art])
            D1 = ddm(ops, SCHLEIFE[art])
            ok_ward &= W == K
            ok_pt &= W.get(1, 0) == P1 and P1 != 0
            ok_ddm &= D1 == P1
            ok_k &= set(W) == set(range(1, n//2 + 1))
            if stich == 0:
                WARD[(name_, n)] = W
                tab_b[f"{name_}, n = {n}"] = {"Potenzen von k": sorted(W), "Terme je Schleifenzahl": anzahl,
                                             "|W_c|/|W_1| (k = 1)": {c: gs(abs(W[c])/abs(W[1])) for c in sorted(W)}}
ZAHLEN["ward_wick_pt"] = tab_b
check("[E] (b) Wardrekursion (nur OPE mit Niveau k als Symbol und Vakuum) = Fermionschleifen (Summe über fixpunktfreie "
      "Permutationen, je Schleife -k mal Spur mal Parke-Taylor-Kette), exakt rational für n = 3, 4, 5, 6 (je drei zufällige "
      "Stromsätze und Punkte) bei su(2)_1 (z99), "
      "su(8)_1 (Original im Quarterrahmen) und so(16)_1 (Majorana, Spur/2): die Potenz k^c zählt genau die c Fermionschleifen",
      ok_ward and ok_k)
check("[E] (b) Der verbundene Einschleifenanteil (k^1) ist exakt die zyklische Parke-Taylor-Summe "
      "-k sum_(sigma in S_(n-1)) Tr(X_1 X_sigma2 ... X_sigman)/(z_1sigma2 z_sigma2sigma3 ... z_sigman1) mit festem Koeffizienten "
      "-k und k = 1 (Spur der Fundamentaldarstellung bzw. Vektorspur/2 = Niveau 1), ohne freien Parameter, und gleich der "
      "DDM-Form -k sum_(S_(n-2)) Tr([[..[X_1, X_s2], ..], X_s(n-1)] X_n) PT(1, sigma, n): nur Strukturkonstanten, die d-Anteile "
      "der Spuren fallen heraus", ok_pt and ok_ddm)
phot_ok = True
for n in range(3, 8):
    zs = punkte(n)
    phot_ok &= sum((pt(zs, list(range(n - 1))[:p] + [n - 1] + list(range(n - 1))[p:]) for p in range(1, n)), Fr(0)) == 0
check("[E] (b) Photonentkopplung der Parke-Taylor-Ketten (sum der Einfügungen eines Punkts in alle Lagen = 0, n = 3..7, "
      "exakt): der Grund, warum die zyklische Spursumme nur f-Strukturen enthält (DDM) und das abelsche Stück verschwindet",
      phot_ok)
vz = {f"n = {n}": {c: ("+" if (-1)**c > 0 else "-") + f"k^{c}" for c in range(1, n//2 + 1)} for n in range(3, 7)}
ZAHLEN["mehrspur_koeffizienten"] = vz
rat = {k_: v["|W_c|/|W_1| (k = 1)"] for k_, v in tab_b.items() if int(k_[-1]) >= 4}
check("[E] (b) Mehrspuranteile (c >= 2 Schleifen = c Spuren) haben je Spurstruktur den Koeffizienten (-1)^c k^c, also "
      "relativ zum Einspuranteil (-1)^(c-1) k^(c-1) = (-1)^(c-1) bei Niveau 1: fest und ungedämpft. Zahl der Strukturen "
      "n = 4: 6 + 3, n = 5: 24 + 20, n = 6: 120 + 90 + 40 + 15 (= 265 fixpunktfreie Permutationen); an den Testpunkten "
      f"|W_2|/|W_1| = {[v.get(2) for v in rat.values()]} (k = 1). Im Twistorstring sind genau diese Terme die konforme "
      "Supergravitation; nativ gibt es keinen Parameter, der sie unterdrückt", ok_ward and
      all(len(v) >= 2 and v[2] > 1e-3 for v in rat.values()))
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

say("  (c) E8 im Quarterrahmen")
WURZELN = []
for i, j in itertools.combinations(range(8), 2):
    for si in (1, -1):
        for sj in (1, -1):
            w = [Fr(0)]*8
            w[i], w[j] = Fr(si), Fr(sj)
            WURZELN.append(tuple(w))
for sg in itertools.product((1, -1), repeat=8):
    if sg.count(-1) % 2 == 0:
        WURZELN.append(tuple(Fr(s, 2) for s in sg))
dot = lambda a, b: sum(x*y for x, y in zip(a, b))
VQ = tuple([Fr(1, 4)]*8)
THETA = tuple([Fr(1, 2)]*8)
ungetw = [b for b in WURZELN if dot(VQ, b).denominator == 1]
getw = [b for b in WURZELN if dot(VQ, b).denominator != 1]
e7 = [b for b in ungetw if dot(b, THETA) == 0]
su2t = [b for b in ungetw if dot(b, THETA) != 0]
ell = [Fr(97 - 11*i, 7 + i) for i in range(8)]
pos = [b for b in e7 if dot(ell, b) > 0]
posset = set(pos)
einfach = [b for b in pos if not any(tuple(x - y for x, y in zip(b, c)) in posset for c in pos)]
Cm = sp.Matrix(len(einfach), len(einfach), lambda i, j: dot(einfach[i], einfach[j]))
grad = [sum(1 for j in range(len(einfach)) if j != i and Cm[i, j] != 0) for i in range(len(einfach))]
dublett = Counter(dot(b, THETA) for b in getw)
minnorm = Counter()
for lam in [tuple([Fr(0)]*8)] + WURZELN:
    d = tuple(x - y for x, y in zip(VQ, lam))
    minnorm[dot(d, d)] += 1
nmin = min(minnorm)
q_gleich = all((dot(VQ, b).denominator == 1) == (dot(THETA, b) % 2 == 0) for b in WURZELN)
ZAHLEN["e8_quarterrahmen"] = {"ungetwistete Wurzeln": len(ungetw), "getwistete Wurzeln": len(getw),
                              "Fixalgebra dim": len(ungetw) + 8, "e7-Wurzeln (theta-orthogonal)": len(e7),
                              "su(2)-Wurzeln": [str(dot(b, THETA)) for b in su2t], "einfache Wurzeln e7": len(einfach),
                              "det Cartan e7": int(Cm.det()), "Knotengrade": sorted(grad),
                              "theta-Ladungen der getwisteten Wurzeln": {str(k): v for k, v in dublett.items()},
                              "kleinste Norm in E8 + v": str(nmin), "Zahl der Grundzustände": minnorm[nmin],
                              "Twistgewicht 8 (alpha - 1/2)^2/2": str(8*Fr(1, 4)**2/2)}
check("[E] (c) Der Quarterzwist v = (1/4)^8 = theta/2 der acht komplexen Fermionen (z75) wirkt auf die 240 E8-Ströme als "
      "Involution e^(2 pi i v.beta) = (-1)^q mit der Flussladung q = theta.beta aus z38 (auf allen Wurzeln geprüft): "
      f"{len(ungetw)} Wurzeln ungetwistet, {len(getw)} getwistet; Fixalgebra der Dimension "
      f"{len(ungetw) + 8} = e7 + su(2): {len(e7)} theta-orthogonale Wurzeln mit {len(einfach)} einfachen, det Cartan = "
      f"{int(Cm.det())}, Knotengrade {sorted(grad)} (E7), dazu su(2) mit den Wurzeln +-theta = +-(1/2)^8; die getwisteten "
      f"Ströme sind (56, 2) (theta-Ladungen {dict((str(k), v) for k, v in dublett.items())}); z38s Zentralisator e7 + u(1) "
      f"ist um die Wurzeln q = +-2 zu e7 + su(2) erweitert. Grundzustand von E8 + v: Norm {nmin} bei {minnorm[nmin]} "
      "Vektoren (v und v - theta), h = 1/4 = 8 x 1/32 = Gewicht des su(2)_1-Dubletts j(j+1)/(k+2): der Quarterrahmen setzt an 0 "
      "und unendlich das su(2)_1-Dublett mal e7-Singulett; der Zustand ist e7 + u(1)-invariant (für die Fermionen allein "
      "u(8) = su(8) + u(1), z75, mit su(8) c e7)", q_gleich and len(ungetw) == 128 and len(getw) == 112 and len(e7) == 126 and len(einfach) == 7 and int(Cm.det()) == 2
      and sorted(grad) == [1, 1, 1, 2, 2, 2, 3] and dublett == {1: 56, -1: 56} and nmin == Fr(1, 2) and minnorm[nmin] == 2
      and 8*Fr(1, 4)**2/2 == Fr(1, 2)*Fr(3, 2)/3)
befund("Teil 1: Die nativen Nahtströme geben auf CP^1 genau Nairs Struktur. Bei gleichen Quartern (z99 und das Original im "
       "Quarterrahmen) kürzt sich der Zwist in jeder Fermionschleife; der verbundene Einschleifenanteil ist die zyklische "
       "Parke-Taylor-Summe mit festem Koeffizienten -k, k = 1, ohne freien Parameter (Wardrekursion, Fockraum und DDM-Form "
       "stimmen exakt überein). Die Mehrspuranteile stehen mit (-1)^(c-1) k^(c-1) = +-1 daneben und sind nativ nicht "
       "unterdrückbar. Exakt Parke-Taylor ist es nur für ungetwistete Ströme: im Quarterrahmen von E8 die e7-Ströme; die "
       "Ströme (56, 2) und das su(2) des Dubletts sehen die Zwistfelder an 0 und unendlich", True)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

# =================================================================== Teil 2
say("\nTeil 2: KMS-Zustand bei beta = 4 pi")
say("  (a) Thermischer Propagator aus den nativen Besetzungen c(x) = 1/(1 + e^(4 pi x)), x in Z + 1/4")
A14 = mp.mpf(1)/4
QNOM = mp.exp(-BETA)


def c_occ(x, beta=BETA, th=0):
    if th == 0:
        return 1/(1 + mp.exp(beta*x))
    return 1/(1 + mp.exp(beta*x - 1j*th))


def dG(Li, Lj, alpha=A14, beta=BETA, th=0, R=16):
    """G_beta - G_0 des zeitgeordneten Propagators <T psi(z_i) psi*(z_j)>, für beide Ordnungen dieselbe meromorphe Funktion:
    (z_i z_j)^(-1/2) [-sum_(x > 0) c(x) (z_j/z_i)^x + sum_(x < 0) (1 - c(x)) (z_j/z_i)^x] (feste Logarithmen L = log z)."""
    s = mp.mpc(0)
    for n in range(R):
        x = alpha + n
        s -= c_occ(x, beta, th)*mp.exp(x*(Lj - Li))
        x = alpha - 1 - n
        s += (1 - c_occ(x, beta, th))*mp.exp(x*(Lj - Li))
    return mp.exp(-(Li + Lj)/2)*s


def G0(zi, zj, Li, Lj, alpha=A14):
    return mp.exp((mp.mpf(1)/2 - alpha)*(Li - Lj))/(zi - zj)


def dG1(Li, Lj, alpha=A14):
    """Koeffizient von c(alpha) in dG: das Teilchen der tiefsten Quartermode x = alpha."""
    return -mp.exp(-(Li + Lj)/2)*mp.exp(alpha*(Lj - Li))


def G_direkt(zi, zj, Li, Lj, alpha=A14, beta=BETA, R=400):
    """Zeitgeordneter thermischer Propagator direkt als Modensumme (radiale Ordnung nach |z|)."""
    xs = [alpha + n for n in range(-R, R)]
    if abs(zi) > abs(zj):
        return mp.exp(-(Li + Lj)/2)*mp.fsum((1 - c_occ(x, beta))*mp.exp(x*(Lj - Li)) for x in xs)
    return -mp.exp(-(Li + Lj)/2)*mp.fsum(c_occ(x, beta)*mp.exp(x*(Lj - Li)) for x in xs)


def G_theta(zi, zj, Li, Lj):
    """Geschlossene Form aus z93 (Bosonisierung): (z_i z_j)^(-1/2) e^(i a u) Theta_1/4(u/2pi)/Theta_1/4(0) (q;q)^2/
    ((e^(iu);q)(q e^(-iu);q)), e^(iu) = z_j/z_i, q = e^(-4 pi), Theta_1/4(w) = theta_3(pi w - i pi/2 | e^(-2 pi))."""
    iu = Lj - Li
    u = -1j*iu
    th = lambda w: mp.jtheta(3, PI*w - 0.5j*PI, mp.exp(-2*PI))
    t = mp.exp(iu)
    osz = mp.qp(QNOM, QNOM)**2/(mp.qp(t, QNOM)*mp.qp(QNOM/t, QNOM))
    return mp.exp(-(Li + Lj)/2)*mp.exp(A14*iu)*th(u/(2*PI))/th(0)*osz


PUNKTE = [mp.mpf(rr)*mp.expj(mp.mpf(ph)) for rr, ph in ((1.35, 0.4), (0.95, 2.2), (0.72, -1.3), (1.12, -2.5), (0.85, 1.05))]
LOGS = [mp.log(z) for z in PUNKTE]
prop_def, theta_def = mp.mpf(0), mp.mpf(0)
for i, j in itertools.permutations(range(4), 2):
    zi, zj, Li, Lj = PUNKTE[i], PUNKTE[j], LOGS[i], LOGS[j]
    gb = G0(zi, zj, Li, Lj) + dG(Li, Lj)
    prop_def = max(prop_def, abs(gb - G_direkt(zi, zj, Li, Lj))/abs(gb))
    if abs(zi) > abs(zj):
        theta_def = max(theta_def, abs(gb - G_theta(zi, zj, Li, Lj))/abs(gb))
C1 = c_occ(A14)
ZAHLEN["thermischer_propagator"] = {"Defekt Modensumme (beide Ordnungen) gegen G_0 + dG": gs(prop_def, 3),
                                    "Defekt gegen Thetaform z93": gs(theta_def, 3), "c(1/4) = 1/(1 + e^pi)": gs(C1, 12),
                                    "c(3/4)": gs(c_occ(3*A14), 6), "c(1/2) (NS-Vergleich)": gs(c_occ(2*A14), 6)}
check("[E] (a) Der zeitgeordnete thermische Propagator des Originals (Modensumme über c(x) bzw. 1 - c(x), beide radialen "
      "Ordnungen, 12 Punktpaare) ist eine meromorphe Funktion G_0 + dG mit dem Vakuumanteil (z/w)^(1/4)/(z - w) und einer "
      "schnell konvergenten Korrektur (Defekt "
      f"{float(prop_def):.0e}); er stimmt mit der Thetaform aus z93 überein (Twist x Theta_1/4 x Oszillatorfaktor, q = e^(-4 pi), "
      f"Defekt {float(theta_def):.0e}). Die tiefste Anregung ist das Quarterteilchen x = 1/4 mit c(1/4) = 1/(1 + e^pi) = "
      f"{float(C1):.10f} (die Gewichte von rho_1), das Loch x = -3/4 hat c(3/4) = {float(c_occ(3*A14)):.2e}",
      prop_def < mp.mpf(10)**-25 and theta_def < mp.mpf(10)**-25)


def korr_mp(Xs, prop, nur_einspur=False):
    """Wick mit skalarem (flavourdiagonalem) Propagator prop(i, j); spurfreie X, also keine Einerzyklen."""
    n = len(Xs)
    tot = mp.mpc(0)
    for perm in derangements(n):
        cs = zyklen(perm)
        if nur_einspur and len(cs) > 1:
            continue
        val = mp.mpc(1)
        for c in cs:
            v = -spur(produkt([Xs[i] for i in c]))
            for a, b in glieder(c):
                v *= prop(a, b)
            val *= v
        tot += val
    return tot


def einspur_erste_ordnung(Xs, P0, P1):
    """Koeffizient von c(1/4) im Einspuranteil: je Schleife ein Glied durch dG1 ersetzt."""
    n = len(Xs)
    tot = mp.mpc(0)
    for s in itertools.permutations(range(1, n)):
        c = (0,) + s
        gl = glieder(list(c))
        w = -spur(produkt([Xs[i] for i in c]))
        for l in range(len(gl)):
            v = mp.mpc(w)
            for l2, (a, b) in enumerate(gl):
                v *= P1(a, b) if l2 == l else P0(a, b)
            tot += v
    return tot


say("  (b) Abweichung von Parke-Taylor im KMS-Zustand")
P0 = lambda a, b: G0(PUNKTE[a], PUNKTE[b], LOGS[a], LOGS[b])
P1 = lambda a, b: dG1(LOGS[a], LOGS[b])


def Pb(beta, th=0):
    return lambda a, b: G0(PUNKTE[a], PUNKTE[b], LOGS[a], LOGS[b]) + dG(LOGS[a], LOGS[b], beta=beta, th=th)


zw_vak = korr_mp([EP, EM], P0)
zw_b = korr_mp([EP, EM], Pb(BETA))
zw_1 = C1/(PUNKTE[0]*PUNKTE[1])
zw_rest = abs(zw_b - zw_vak - zw_1)/abs(zw_1)
pt_vak_ok = abs(zw_vak - 1/(PUNKTE[0] - PUNKTE[1])**2) < mp.mpf(10)**-25
XS = {3: [3*EP - 2*EM + HH, -EP + 4*EM + 2*HH, 2*EP + EM - 3*HH],
      4: [3*EP - 2*EM + HH, -EP + 4*EM + 2*HH, 2*EP + EM - 3*HH, EP - EM + 2*HH]}
dev = {}
for n, Xs in XS.items():
    s0 = korr_mp(Xs, P0, True)
    sb = korr_mp(Xs, Pb(BETA), True)
    d1 = einspur_erste_ordnung(Xs, P0, P1)
    s8 = korr_mp(Xs, Pb(2*BETA), True)
    c8 = c_occ(A14, 2*BETA)
    dev[n] = {"|S_beta - S_0|/|S_0|": gs(abs(sb - s0)/abs(s0), 6),
              "(S_beta - S_0)/(c(1/4) D_1) - 1, beta = 4 pi": gs(abs((sb - s0)/(C1*d1) - 1), 4),
              "dasselbe bei beta = 8 pi": gs(abs((s8 - s0)/(c8*d1) - 1), 4)}
ZAHLEN["kms_abweichung"] = {"Zweipunkt: Rest nach c(1/4) Tr/(z1 z2), relativ": gs(zw_rest, 4), "Einspur n": dev}
check("[E] (b) Zweipunktantwort bei beta = 4 pi: <J^X(z_1) J^Y(z_2)>_beta = Tr(XY) [1/z_12^2 + c(1/4)/(z_1 z_2) + "
      f"O(c(1/4)^2)] (Rest relativ {float(zw_rest):.3f}, also zweite Ordnung ~ c(1/4)); der Vakuumanteil ist exakt 1/z_12^2. Die "
      "erste Korrektur ist die Ladungsschwankung der Nullmode, Var Q ~ c(1/4) (z93), und hängt an 0 und unendlich",
      pt_vak_ok and zw_rest < 0.1)
check("[E] (b) Einspuranteil n = 3, 4 (su(2), zwei Nähte gleicher Quarter): relative Abweichung von Parke-Taylor "
      f"{[dev[n]['|S_beta - S_0|/|S_0|'] for n in XS]} an den Testpunkten; sie ist in erster Ordnung genau c(1/4) D_1 mit "
      "D_1 = Summe über die Glieder der Schleife, eines durch -(z_i z_j)^(-1/2) (z_j/z_i)^(1/4) ersetzt (Rest "
      f"{[dev[n]['(S_beta - S_0)/(c(1/4) D_1) - 1, beta = 4 pi'] for n in XS]}, bei beta = 8 pi "
      f"{[dev[n]['dasselbe bei beta = 8 pi'] for n in XS]}): fester Koeffizient c(1/4) = 1/(1 + e^pi), kein Parameter",
      all(0.002 < dev[n]["|S_beta - S_0|/|S_0|"] < 0.5 and dev[n]["(S_beta - S_0)/(c(1/4) D_1) - 1, beta = 4 pi"] < 0.2
          and dev[n]["dasselbe bei beta = 8 pi"] < 0.1*dev[n]["(S_beta - S_0)/(c(1/4) D_1) - 1, beta = 4 pi"] for n in XS))

say("  (c) Ordnung der Abweichung: Exponent aus der beta-Abhängigkeit, Gegenprobe mit Ladung 0")


def Zfak(beta, th, R=16):
    z1 = mp.mpc(1)
    for n in range(R):
        z1 *= (1 + mp.exp(-beta*(A14 + n) + 1j*th))*(1 + mp.exp(-beta*(1 - A14 + n) - 1j*th))
    return z1**2                                                                    # zwei Nähte, gemeinsames theta


def projiziert(Xs, beta, Mth=24, ladung=False):
    num, den = mp.mpc(0), mp.mpc(0)
    for j in range(Mth):
        th = 2*PI*j/Mth
        Z = Zfak(beta, th)
        if ladung:
            val = 2*mp.fsum(c_occ(A14 + n, beta, th) - (1 - c_occ(A14 - 1 - n, beta, th)) for n in range(16))
        else:
            val = korr_mp(Xs, Pb(beta, th))
        num += Z*val
        den += Z
    return num/den


expo = {}
b1, b2 = BETA, BETA + PI/2
for n in (2, 3):
    Xs = [EP, EM] if n == 2 else XS[3]
    s0 = korr_mp(Xs, P0)
    d_nat = [korr_mp(Xs, Pb(b)) - s0 for b in (b1, b2)]
    d_q0 = [projiziert(Xs, b) - s0 for b in (b1, b2)]
    expo[n] = {"E_eff nativ": gs(mp.log(abs(d_nat[0])/abs(d_nat[1]))/(b2 - b1), 6),
               "E_eff Ladung 0": gs(mp.log(abs(d_q0[0])/abs(d_q0[1]))/(b2 - b1), 6),
               "|Abweichung|/|PT| nativ": gs(abs(d_nat[0])/abs(s0), 4),
               "|Abweichung|/|PT| Ladung 0": gs(abs(d_q0[0])/abs(s0), 4)}
q_proj = projiziert(None, BETA, ladung=True)
steig_c1 = mp.log(c_occ(A14, b1)/c_occ(A14, b2))/(b2 - b1)
ZAHLEN["exponent"] = {"n": expo, "<Q>_(Ladung 0)": gs(abs(q_proj), 3), "Steigung von ln c(1/4)": gs(steig_c1, 6),
                      "e^(-pi)": gs(mp.exp(-PI), 6),
                      "e^(-4 pi)": gs(mp.exp(-4*PI), 6)}
check("[E] (c) Die Abweichung skaliert wie e^(-beta E_eff) mit E_eff = "
      f"{[expo[n]['E_eff nativ'] for n in expo]} (n = 2, 3; beta = 4 pi gegen 4.5 pi), der Steigung von ln c(1/4) "
      f"({float(steig_c1):.4f}, etwa (1 - c)/4) bis auf die zweite Ordnung: die Quarterenergie 1/4, also Ordnung "
      f"e^(-pi) = {float(mp.exp(-PI)):.4f} bei beta = 4 pi, nicht e^(-4 pi). Projiziert man das Original auf Gesamtladung 0 "
      f"(nicht nativ; <Q> = {float(abs(q_proj)):.0e}), fallen die Einquarkzustände weg: E_eff = "
      f"{[expo[n]['E_eff Ladung 0'] for n in expo]}, Abweichung {[expo[n]['|Abweichung|/|PT| Ladung 0'] for n in expo]} "
      f"~ e^(-4 pi) = {float(mp.exp(-4*PI)):.2e}",
      all(abs(expo[n]["E_eff nativ"] - float(steig_c1)) < 0.01 and abs(expo[n]["E_eff Ladung 0"] - 1) < 0.001
          for n in expo)
      and abs(q_proj) < 1e-20)
befund("Teil 2: Im KMS-Zustand bei beta = 4 pi weicht der Stromkorrelator von Parke-Taylor um die Ordnung e^(-pi) ab, mit "
       "festem Koeffizienten c(1/4) = 1/(1 + e^pi) = 4,14 % (Quarterteilchen x = 1/4, dieselbe Zahl wie die Gewichte von "
       "rho_1). Die Korrektur ist kein Vielfaches von Parke-Taylor, sondern ein Ruhesystemterm an 0 und unendlich (Zweipunkt: "
       "c(1/4) Tr/(z_1 z_2)). e^(-4 pi) ergibt sich nur mit Projektion auf Ladung 0, die nativ nicht vorliegt (z93: der "
       "Ladungssektor ist Teil des Zustands)", True)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

# =================================================================== Teil 3
say("\nTeil 3: 4D-Lesart auf Twistorgeraden")
say("  (a) z_i -> lambda_i: Parke-Taylor-Nenner, projektive Gewichte, Lorentzgruppe")
lam_ok, gew_ok = True, True
for n in range(3, 7):
    a = sp.symbols(f"a0:{n}")
    b = sp.symbols(f"b0:{n}")
    t = sp.symbols(f"t0:{n}")
    ang = lambda i, j: a[i]*b[j] - b[i]*a[j]
    zz = [b[i]/a[i] for i in range(n)]
    PTz = 1/sp.prod([zz[i] - zz[(i + 1) % n] for i in range(n)])
    PTl = (-1)**n/sp.prod([ang(i, (i + 1) % n) for i in range(n)])
    lam_ok &= sp.cancel(PTz/sp.prod([a[i]**2 for i in range(n)]) - PTl) == 0
    sub = {**{a[i]: t[i]*a[i] for i in range(n)}, **{b[i]: t[i]*b[i] for i in range(n)}}
    gew_ok &= sp.cancel(PTl.subs(sub, simultaneous=True)/PTl*sp.prod([t[i]**2 for i in range(n)])) == 1
s_ = sp.Symbol("s")
fa, fb = sp.Function("a")(s_), sp.Function("b")(s_)
mass_ok = sp.simplify(fa*sp.diff(fb, s_) - fb*sp.diff(fa, s_) - fa**2*sp.diff(fb/fa, s_)) == 0
bra = lambda l1, l2: l1[0]*l2[1] - l1[1]*l2[0]
lams = [np.array([complex(*rng.normal(size=2)), complex(*rng.normal(size=2))]) for _ in range(5)]
Mz = rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2))
Mz = Mz/np.sqrt(np.linalg.det(Mz))
Md = np.diag([1.7 - 0.4j, 1/(1.7 - 0.4j)])
ptl = lambda L: 1/np.prod([bra(L[i], L[(i + 1) % len(L)]) for i in range(len(L))])
O_, I_ = np.array([1, 0]), np.array([0, 1])
ruhe = lambda L: bra(L[0], L[1])**2/(bra(O_, L[0])*bra(I_, L[0])*bra(O_, L[1])*bra(I_, L[1]))
lor = {"PT unter SL(2,C)": abs(ptl([Mz @ l for l in lams])/ptl(lams) - 1),
       "KMS-Term <12>^2/(<o1><i1><o2><i2>) unter SL(2,C)": abs(ruhe([Mz @ l for l in lams])/ruhe(lams) - 1),
       "KMS-Term unter diag(s, 1/s)": abs(ruhe([Md @ l for l in lams])/ruhe(lams) - 1)}
ZAHLEN["twistorgerade"] = {k: gs(v, 3) for k, v in lor.items()}
check("[E] (a) Mit z_i = lambda_i1/lambda_i0 und J(z) dz = J(lambda) <lambda d lambda> (<lambda d lambda> = lambda_0^2 dz, "
      "symbolisch) wird der Einspuranteil exakt (-1)^n/(<12><23>...<n1>) (sympy, n = 3..6): jeder Strom trägt das projektive "
      "Gewicht -2 (Skalierung lambda_i -> t_i lambda_i gibt prod t_i^(-2)), das Gewicht eines Gluons positiver Helizität; der "
      f"Nenner ist SL(2,C)-invariant (Defekt {lor['PT unter SL(2,C)']:.0e}). Der farbgeordnete MHV-Nenner ist damit nativ. "
      "Der KMS-Term aus Teil 2 wird c(1/4) Tr <12>^2/(<o1><iota1><o2><iota2>) relativ zum Nenner: er braucht die Spinoren "
      f"o, iota der Pole (Ruhesystem), ändert sich unter SL(2,C) (relativ {lor['KMS-Term <12>^2/(<o1><i1><o2><i2>) unter SL(2,C)']:.2f}) "
      f"und ist nur unter der Diagonalgruppe invariant (Defekt {lor['KMS-Term unter diag(s, 1/s)']:.0e}): Drehung und Boost "
      "um die Achse des Zustands", lam_ok and gew_ok and mass_ok and lor["PT unter SL(2,C)"] < 1e-12
      and lor["KMS-Term <12>^2/(<o1><i1><o2><i2>) unter SL(2,C)"] > 1e-3 and lor["KMS-Term unter diag(s, 1/s)"] < 1e-12)

say("  (b) Helizitätszähler: Grassmann-Nullmoden der Geraden")


def gmul(A, B):
    out = {}
    for ma, ca in A.items():
        for mb, cb in B.items():
            if ma & mb:
                continue
            s, mm = 0, mb
            while mm:
                j = (mm & -mm).bit_length() - 1
                s += bin(ma >> (j + 1)).count("1")
                mm &= mm - 1
            out[ma | mb] = out.get(ma | mb, 0) + (-1)**s*ca*cb
    return {k: v for k, v in out.items() if v != 0}


def berezin(lams_, traeger, N):
    """int d^(2N) theta prod_k prod_(A in traeger_k) (theta^A . lambda_k): Nullmoden psi^A = theta^(A alpha) lambda_alpha."""
    prod = {0: 1}
    for k, As in enumerate(traeger):
        for A in As:
            prod = gmul(prod, {1 << (2*A): lams_[k][0], 1 << (2*A + 1): lams_[k][1]})
    return prod.get((1 << (2*N)) - 1, 0)


lq = [(Fr(int(rng.integers(-9, 10))), Fr(int(rng.integers(-9, 10)))) for _ in range(5)]
br = lambda i, j: lq[i][0]*lq[j][1] - lq[i][1]*lq[j][0]
nair_ok, nair_tab = True, {}
for N in range(1, 5):
    alle = list(range(N))
    for nm in range(0, 4):
        for minus in itertools.combinations(range(5), nm):
            val = berezin(lq, [alle if k in minus else [] for k in range(5)], N)
            if nm == 2:
                i, j = minus
                nair_ok &= abs(val) == abs(br(i, j)**N) and val != 0
            else:
                nair_ok &= val == 0
    nair_tab[f"N = {N}"] = {"Zähler": f"<ij>^{N}", "Helizität der eta^N-Einsetzung": str(1 - Fr(N, 2)),
                            "Gewicht des Maßes 4 - N": 4 - N}
gem = berezin(lq, [[0, 1, 2, 3], [0, 1], [2, 3], []], 4)
nair_ok &= abs(gem) == abs(br(0, 1)**2*br(0, 2)**2)
a4 = sp.symbols("p0:4")
b4 = sp.symbols("q0:4")
tt = sp.Symbol("t")
ang4 = lambda i, j: a4[i]*b4[j] - b4[i]*a4[j]
hel_ok = True
for N in range(1, 5):
    A_ = ang4(0, 1)**N/sp.prod([ang4(i, (i + 1) % 4) for i in range(4)])
    for k, soll in ((0, N - 2), (2, -2)):
        sk = sp.cancel(A_.subs({a4[k]: tt*a4[k], b4[k]: tt*b4[k]}, simultaneous=True)/A_)
        hel_ok &= sp.simplify(sk - tt**soll) == 0
ZAHLEN["nair_nullmoden"] = nair_tab
check("[E] (b) Vergleichsziel nachgerechnet: 2N Grassmann-Nullmoden der Geraden (psi^A = theta^(A alpha) lambda_alpha, "
      "A = 1..N, Exterioralgebra) geben bei genau zwei eta^N-Einsetzungen i, j den Zähler <ij>^N und sonst 0 (alle Plus, eine, "
      "drei; N = 1..4, fünf Punkte); gemischte Verteilungen geben Produkte (<01>^2 <02>^2 für N = 4). Mit dem nativen Nenner hat "
      "die eta^N-Einsetzung das Gewicht N - 2, also Helizität 1 - N/2, die übrigen +1 (sympy). Ein Gluon negativer Helizität "
      "verlangt N = 4; genau dann hat das Maß d^4Z d^N psi das Gewicht 4 - N = 0 (Calabi-Yau, CP^(3|4))", nair_ok
      and hel_ok)

say("  (c) Native Kandidaten für die Nullmoden")
nullmoden_quarter = {str(a): sum(1 for n in range(-50, 50) if a + n == 0) for a in (Fr(1, 4), Fr(3, 4), Fr(1, 2), Fr(0))}
skal = mp.mpf(3)
gew_psi = -mp.log(abs((G0(skal*PUNKTE[0], skal*PUNKTE[1], LOGS[0] + mp.log(skal), LOGS[1] + mp.log(skal))
                       / G0(PUNKTE[0], PUNKTE[1], LOGS[0], LOGS[1]))))/mp.log(skal)/2
inv = []
Ms = []
for _ in range(3):
    Mr = rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2))
    Ms.append(Mr/np.sqrt(np.linalg.det(Mr)))
zeilen = []
for Mr in Ms:
    for p_, q_ in itertools.product(range(2), repeat=2):
        zeile = np.zeros(4, dtype=complex)
        for a_, b_ in itertools.product(range(2), repeat=2):
            zeile[2*a_ + b_] += Mr[a_, p_]*Mr[b_, q_]
        zeile[2*p_ + q_] -= 1
        zeilen.append(zeile)
_, sv, vh = np.linalg.svd(np.array(zeilen))
kern = vh[-1].reshape(2, 2) if sv[-1] < 1e-10 else None
dim_inv = int(np.sum(sv < 1e-10))
sym_teil = np.abs(kern + kern.T).max() if kern is not None else 1.0
herm = lambda l1, l2: np.vdot(l1, l2)
herm_def = abs(herm(Mz @ lams[0], Mz @ lams[1])/herm(lams[0], lams[1]) - 1)
ZAHLEN["native_kandidaten"] = {
    "Nahtfermion: Gewicht aus G_0 (Skalierung)": gs(gew_psi, 10), "Bündel auf der Geraden": "O(-1), h^0 = 0",
    "Nullmoden x = 0 je Rahmen alpha": nullmoden_quarter,
    "SL(2,C)-invariante Bilinearformen: Dimension": dim_inv, "|Symmetrieanteil| der invarianten Form": gs(sym_teil, 3),
    "CCR-Hilbertform lambda^+ lambda unter SL(2,C), relativ": gs(herm_def, 3),
    "Quarterpaar": "zwei komplexe Fermionen (c = 1 je Quarter): N <= 2, Helizität >= 0",
    "eta-Dublett": "eine Mode {Omega, psi*_1/4 Omega}: N <= 1, Helizität 1/2",
    "Schrittart S_a": "CCR (bosonisch, z83), Paarung delta_ab (Hilbert), Sym^4 S hebt die Helizität kinematisch"}
check("[E] (c) Native Kandidaten liefern die vier Grassmann-Nullmoden nicht: (i) Die Nahtfermionen haben Gewicht "
      f"{float(gew_psi):.6f} (aus der Skalierung des nativen Propagators), sind also Schnitte von O(-1) mit h^0 = 0; im "
      f"Quarterrahmen gibt es keine Mode x = 0 ({nullmoden_quarter}); nur der R-Rahmen alpha = 0 hätte eine, ohne lambda-"
      "Abhängigkeit. Sie bauen den Nenner und können nicht zugleich die O(1)-Moduln sein. (ii) Die Schrittart S_a ist CCR, "
      "also Grassmann-gerade; eine gaußsche Paarung ist symmetrisch, die einzige SL(2,C)-invariante Bilinearform auf C^2 ist "
      f"aber eps (Dimension {dim_inv}, Symmetrieanteil {sym_teil:.0e}); die native CCR-Paarung delta_ab ist die Hilbertform "
      f"lambda^+ lambda und nicht lorentzinvariant (relativ {herm_def:.2f}). (iii) Quarterpaar (zwei Spezies) und eta-Dublett "
      "(eine Mode) gäben höchstens N = 2 bzw. 1, also Einsetzungen der Helizität 0 bzw. 1/2, keine Gluonen",
      abs(gew_psi - 0.5) < 1e-20 and nullmoden_quarter["1/4"] == 0 and nullmoden_quarter["3/4"] == 0
      and nullmoden_quarter["0"] == 1 and dim_inv == 1 and sym_teil < 1e-10 and herm_def > 1e-3)
befund("Teil 3: Auf der Twistorgeraden ist der native Einspuranteil exakt der farbgeordnete MHV-Nenner 1/(<12>...<n1>) mit "
       "Gewicht -2 je Strom, lorentzinvariant. Der Helizitätszähler <ij>^4 fehlt: Er verlangt vier Grassmann-ungerade "
       "O(1)-Richtungen (N = 4, das Calabi-Yau-Maß von CP^(3|4)); dieselbe Struktur liefert die MHV-Auswahlregel (genau zwei "
       "negative Helizitäten auf der Geraden). Ohne Nullmoden (reine Stromalgebra, Wellenfunktionen vom Gewicht 0) gäbe die "
       "Gerade nur Allplus-Einsetzungen, also PT delta^4(P), was kein Yang-Mills-Baum ist (Vergleich: Allplus verschwindet "
       "dort). Nativ gibt es keine solche Richtung: Nahtfermionen sind O(-1), die Schrittart ist bosonisch mit Hilbertpaarung, "
       "Quarterpaar und eta-Dublett zählen höchstens zwei", True)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

# =================================================================== Teil 4
say("\nTeil 4: Kopplung und Urteil")
skal_tab, skal_ok = {}, True
for (name_, n), W in WARD.items():
    if not name_.startswith("su(8)"):
        continue
    exps = {c: str(Fr(c) - Fr(n, 2)) for c in W}
    skal_tab[f"n = {n}"] = {"k-Potenzen nach J -> J/sqrt(k)": exps,
                            "Einspur": str(1 - Fr(n, 2)), "Zweispur/Einspur": "k" if 2 in W else "-"}
    skal_ok &= set(W) == set(range(1, n//2 + 1))
moden_gleich = sorted({str(Fr(1, 4) - Fr(1, 4) + m) for m in range(-2, 3)})
moden_versch = sorted({str(Fr(3, 4) - Fr(1, 4) + m) for m in range(-2, 2)})
uhr = {"gleiche Quarter: Moden der Ströme": moden_gleich, "R = e^(2 pi i L_0) auf Strommoden": 1,
       "verschiedene Quarter: Moden von E^+-": moden_versch, "R auf E^+-": -1}
grad = {d: {"bosonische Moduln 4d": 4*d, "Grassmann-Nullmoden N(d+1), N = 4": 4*(d + 1),
            "negative Helizitäten (d + 1)": d + 1} for d in (1, 2, 3)}
ZAHLEN["kopplung"] = {"Niveauskalierung (su(8)_1)": skal_tab, "Uhr": uhr, "Grad der Kurve": grad}
check("[E] (a) Niveau: Nach J -> J/sqrt(k) trägt der c-Spurterm der n-Punkt-Antwort k^(c - n/2) (alle c = 1..n/2 kommen "
      "vor, n = 3..6): Einspur k^(1 - n/2) = g^(n-2) mit g_eff = k^(-1/2), Zweispur relativ k. Nativ ist k = 1 (su(2)_1, "
      "su(8)_1, E8_1), also g_eff = 1 in der Normierung der Stromzweipunktfunktion und die Mehrspur- (Gravitations-)Terme "
      "gleich stark. Das fixiert g nur unter der Prämisse P_norm (4D-Zustände = stromnormierte Einsetzungen)", skal_ok)
check("[E] (b) Uhr: R = e^(2 pi i L_0) wirkt auf die ganzzahligen Strommoden gleicher Quarter als 1 und auf die getwisteten "
      f"(E^+- bei verschiedenen Quartern, Moden {moden_versch}) als -1: eine Z_2-Auswahlregel (gerade Zahl getwisteter "
      "Ströme), kein Betrag", uhr["R auf E^+-"] == -1 and all(Fr(m).denominator == 2 for m in moden_versch)
      and all(Fr(m).denominator == 1 for m in moden_gleich))
befund("Teil 4: g ist nativ nur über das Niveau festgelegt (g_eff = k^(-1/2) = 1 unter P_norm), mit gleich starken "
       "Mehrspurtermen. Die Uhr gibt nur eine Z_2-Regel, der KMS-Zustand nur den festen additiven Ruhesystemterm c(1/4), kein "
       "Vielfaches g^(n-2). Die native Naht ist eine Gerade (Grad 1, vier bosonische Moduln = Raumzeit); mit N = 4 hätte sie "
       "8 Grassmann-Nullmoden, also genau zwei negative Helizitäten (MHV). N^kMHV bräuchte Kurven vom Grad k + 1 mit einem "
       "relativen Gewicht (Instantonfaktor), das nativ nicht vorkommt (P_deg)", True)

STATUS = {
    "Stromkorrelatoren der Naht = Nair-Struktur": ("E", "Fockraum = Wickformel in Ortslage; bei gleichen Quartern kürzt sich "
                                                        "der Zwist, Einspuranteil = zyklische Parke-Taylor-Summe, Koeffizient "
                                                        "-k, k = 1, DDM-Form, Photonentkopplung, n = 3..6 exakt"),
    "Mehrspuranteile": ("E", "(-1)^c k^c je Struktur, relativ (-1)^(c-1) bei k = 1; im Twistorstring konforme "
                             "Supergravitation; nativ nicht unterdrückbar"),
    "Quarterrahmen von E8": ("E", "Zwist (-1)^q der Flussladung aus z38: Fixalgebra e7 + su(2) (136), getwistet (56, 2) "
                                  "(112); Grundzustand = su(2)_1-Dublett (h = 1/4 = 8/32), Zustand e7 + u(1)-invariant: "
                                  "exakt Parke-Taylor nur für e7-Ströme (darin su(8))"),
    "verschiedene Quarter (Originalpaar 1/4, 3/4)": ("E", "E^+- halbzahlig, Z_2-getwistet, Verzweigung an 0 und unendlich; "
                                                          "das ungetwistete su(2)_1 ist nach Ladungskonjugation das "
                                                          "Paarungs-su(2)"),
    "KMS-Abweichung": ("E", "Ordnung e^(-pi) (E_eff = 1/4), fester Koeffizient c(1/4) = 1/(1 + e^pi) = 4,14 %, "
                            "Ruhesystemterm (Zweipunkt c(1/4) Tr/(z_1 z_2)); e^(-4 pi) nur mit Ladungsprojektion"),
    "MHV-Nenner auf der Geraden": ("E unter Lesart T", "1/(<12>...<n1>), Gewicht -2 je Strom, SL(2,C)-invariant"),
    "Helizitätszähler <ij>^4": ("O / Prämisse P_N4", "verlangt vier Grassmann-ungerade O(1)-Nullmodenrichtungen "
                                                     "(CP^(3|4)); nativ keine: Nahtfermionen O(-1), S_a bosonisch mit "
                                                     "Hilbertpaarung, Quarterpaar/eta-Dublett N <= 2 (Helizität >= 0)"),
    "Kopplung g": ("C", "Niveau k = 1 gibt g_eff = k^(-1/2) = 1 nur unter P_norm; Uhr und KMS fixieren keinen Betrag; "
                        "Gradgewichte (N^kMHV) nicht nativ"),
    "Hypothese T im Eichsektor": ("C", "Nenner und Farbstruktur nativ mit festem Koeffizienten; Zähler (P_N4), "
                                       "Normierung (P_norm) und Gradgewichte fehlen; selbst mit ihnen gäbe die Naht N = 4 SYM "
                                       "plus konforme Supergravitation gleicher Stärke, nicht reines Yang-Mills"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
urteil = (
    "Hypothese T ist für den Eichsektor bedingt, nicht nativ. Nativ und exakt: Die Nahtströme liefern auf CP^1 genau Nairs "
    "Struktur. Der Fockraum zweier Nähte reproduziert die Wickformel in Ortslage mit dem nativen Propagator "
    "(z/w)^(1/2 - alpha)/(z - w); bei gleichen Quartern (z99, und das Original im Quarterrahmen) kürzt sich der Zwist in jeder "
    "Fermionschleife, und der verbundene Einschleifenanteil ist für n = 3..6 exakt die zyklische Parke-Taylor-Summe "
    "-k sum Tr(X_1 X_s2 ... X_sn)/(z_1s2 ... z_sn1) mit k = 1, ohne freien Parameter (Wardrekursion, Fermionschleifen und "
    "DDM-Form stimmen rational überein, auch für so(16)_1 der 16 Majorana). Die Mehrspuranteile stehen mit (-1)^(c-1) k^(c-1) "
    "= +-1 daneben; im Twistorstring sind sie die konforme Supergravitation, und nativ gibt es keinen Parameter, der sie "
    "dämpft. Beim Originalpaar (1/4, 3/4) ist das Paarungs-su(2)_1 (psi*_0 psi*_1, Q_0 + Q_1) ungetwistet, der Transferstrom "
    "E^+- getwistet. Der Quarterrahmen der E8-Naht ist v = (1/4)^8 = theta/2, also der Zwist (-1)^q mit der Flussladung aus "
    "z38: Fixalgebra e7 + su(2), getwistete Ströme (56, 2), Grundzustand das su(2)_1-Dublett (h = 1/4 = 8 x 1/32); exakt "
    "Parke-Taylor sind nur die e7-Ströme (darin su(8)). Im KMS-Zustand bei beta = 4 pi weicht die Antwort um die Ordnung e^(-pi) ab (Exponent 1/4 aus "
    "der beta-Abhängigkeit), mit dem festen Koeffizienten c(1/4) = 1/(1 + e^pi) = 4,14 % des Quarterteilchens; die Korrektur "
    "ist ein Ruhesystemterm an den Polen (Zweipunkt c(1/4) Tr/(z_1 z_2), auf der Geraden <12>^2/(<o1><iota1><o2><iota2>)), "
    "nur unter Drehung und Boost um die Zustandsachse invariant. e^(-4 pi) entsteht erst mit einer Projektion auf Ladung 0, "
    "die nativ nicht vorliegt. Auf der Twistorgeraden wird der Einspuranteil exakt der farbgeordnete MHV-Nenner "
    "1/(<12>...<n1>) mit Gewicht -2 je Strom. Es fehlt der Helizitätszähler <ij>^4: Er verlangt vier Grassmann-ungerade "
    "O(1)-Nullmodenrichtungen (Calabi-Yau-Maß von CP^(3|4), Gewicht 4 - N = 0), die zugleich die MHV-Auswahlregel geben. "
    "Nativ fehlen sie: Die Nahtfermionen sind O(-1) ohne Nullmode im Quarterrahmen und bauen schon den Nenner; die Schrittart "
    "S_a ist bosonisch mit Hilbertpaarung (die einzige lorentzinvariante Paarung eps ist antisymmetrisch); Quarterpaar und "
    "eta-Dublett gäben höchstens N = 2 bzw. 1, also Einsetzungen der Helizität 0 bzw. 1/2. Die Kopplung ist nur über das "
    "Niveau festgelegt, g_eff = k^(-1/2) = 1 unter der Prämisse P_norm; Uhr und KMS fixieren keinen Betrag, Gradgewichte für "
    "N^kMHV kommen nativ nicht vor (die Naht ist eine Gerade). Benannte Prämissen: P_N4 (vier Grassmann-ungerade "
    "O(1)-Richtungen), P_norm, für Nicht-MHV P_deg. Selbst mit ihnen gäbe die Naht N = 4 SYM plus konforme Supergravitation "
    "gleicher Stärke, nicht reines Yang-Mills. Folgen: K15 bleibt L, präzisiert: Farbstruktur und MHV-Nenner der "
    "Gluonwechselwirkung sind native Antworten der Naht mit festem Koeffizienten, Zähler und Entkopplung der Gravitation "
    "nicht. Für T3/T5: Hypothese T liefert die fehlende Vorwärtsregel nur für die Nenner; die feste Dreikörperstruktur "
    "f^abc/(z_12 z_23 z_31) ist da, aber als Amplitude auf komplexer bzw. (2,2)-Kinematik, nicht als lokaler Vertex (das "
    "No-go aus z99 bleibt). Nichts davon ist befördert.")
say("\nUrteil\n   " + urteil)
say(f"\n   (Laufzeit {time.time() - T0:.1f} s)")
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z110_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
