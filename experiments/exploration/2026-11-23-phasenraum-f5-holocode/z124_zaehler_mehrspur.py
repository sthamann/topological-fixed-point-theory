"""Z124: Zähler und Mehrspur -- Gradzwist der Nahtfermionen, Pati-Salam-Quartett als eta^A, Mehrspurterme als konforme
Supergravitation (Folge zu z110, Hypothese T).

Leitlinie 'Rekonstruktion statt Auswahl': Gerechnet wird aus den nativen Nahtfermionen (Moden x in Z + alpha, Fermikante
x < 0, Quarterrahmen alpha = 1/4, z58/z110), der nativen Uhr R = e^(2 pi i L_0) (z113), der Teilchen-Loch-Konjugation (der
CAR-Anteil von Tomitas J, z116 J_F), der nativen Stromalgebra auf CP^1 (Wardrekursion, Fermionschleifen, z110) und der
Cliffordrealisierung der 16 auf fünf Moden (z116). Vergleichsziele, keine Prämissen: Nairs Grassmann-Nullmoden auf der
Twistorgeraden, der Twistorstring (Witten 2003, Berkovits 2004: CP^(3|4)), Berkovits-Witten 2004 (Mehrspurterme =
konforme Supergravitation), die topologische Verdrillung T -> T + a dJ, die Reduktion der 10D-Majorana-Weyl-16 auf
Spin(4) x Spin(6), DDM/BCJ/KLT sowie Maxwell-Spannungstensoren mit Einstein- (1/p^2) und Weyl-Propagator (1/p^4). Firewall
DIMENSION.UPLIFT.FIREWALL.01 wie z110: 4D tritt nur als Lesart T auf (Naht = Twistorgerade). H6 und die gelösten Blöcke
(Flavour, Massen, alpha, die 16, E8) werden benutzt, nicht geöffnet; nichts wird befördert. Statuskürzel: E exakt,
C bedingt, L Lesart, O offen.

Bestand: z110 (Parke-Taylor-Nenner nativ mit Koeffizient -k, k = 1; Zähler <ij>^4 braucht P_N4, vier Grassmann-ungerade
O(1)-Richtungen, 4 - N = 0; Nahtfermionen O(-1) ohne Mode x = 0 im Quarterrahmen; S_a bosonisch mit Hilbertpaarung;
Quarterpaar N <= 2, eta-Dublett N <= 1; Mehrspuranteile (-1)^c k^c bei k = 1 ungedämpft; E8 im Quarterrahmen e7 + su(2),
(56, 2), Dublett), z112 (Hodges = KLT aus nativen Klammern, Prämisse P_diag; Original x Partner = Phi), z113 (Uhr =
2D-Drehung, e^(2 pi i L_0) = R; native GSO R = R-quer, (-1)^(F+F-quer) = +1; P_GSO nicht nativ), z114 (O(-1) einziges
azyklisches Bündel auf CP^1), z116 (16 = Lambda^gerade C^5, R = i Gamma_11, Hyperladung y = (-1/3)^3 (1/2)^2), z64/z70
(O(1) + O(1) als T_C (x) pi_* theta, Bündelebene), z100 (native Cartan-Verbindungen flach), z111 (U(2)-Farbklebung).

Teil 1  Gradzwist statt Holonomiezwist. (a) Cech-Kohomologie von O(d) auf CP^1. (b) Nullmoden der Nahtfermionen ohne und
        mit Gradzwist q (Eingangs- und Ausgangsgrundzustand von L_0 -+ q Q über der nativen Fermikante) in den Rahmen
        alpha = 1/2, 1/4, 3/4; Wellenfunktionen z^(1/2 - alpha) (1, z); SL(2)-Abgeschlossenheit. (c) Berezin-Zähler mit
        diesen Nullmoden. (d) Verdrillte Virasoro-Algebra: c' = c - 12 a^2 k, Gewichte 1/2 -+ a. (e) Native Kandidaten auf dem
        Fockraum des Originalpaars: Uhr R, Teilchen-Loch-J, Spin-Z4 = Holonomie. (f) Ladungsanomalie: Mit Gradzwist
        verschwinden neutrale Stromkorrelatoren.
Teil 2  Pati-Salam-Quartett als eta^A. (a) 16 = (4,2,1) + (4-quer,1,2) aus der Cliffordrealisierung; Supermultiplett
        Lambda(C^4) nur unter SU(4) c SU(5). (b) Statistik: R^2 = -1 auf der 16; bosonische Kontrolle. (c) Homogenität:
        Berezin-Zähler für h = +1, 0, -1, dazu die E8-Ströme (16,4) (h = -2, bosonisch). (d) Die Gleichsetzung
        SU(2)_L = SL(2)_lambda ist eine topologische Verdrillung der Materie. (e) Gluonzähler: genau zwei
        eta^4-Einsetzungen, GL(4)-Kovarianz (eta als Hilfsgrößen), N = 4 nötig.
Teil 3  Mehrspurterme. (a) Niveau: Gewicht k^(c-1), M Nahtpaare diagonal = Niveau M, großes Niveau = freier Grenzfall.
        (b) Uhr/GSO: Doppelspur = Vakuumaustausch, native Projektoren lassen ihn stehen; die Uhr entfernt die getwisteten
        Transferströme selbst. (c) Farbstruktur: Rang der Spurstrukturen, U(1) der U(2) koppelt nur über Mehrspur.
        (d) P_diag/Doppelkopie: DDM-Projektion, BCJ, KLT für su(2) gegen su(3). (e) n = 4: Doppelspur = konforme
        Supergravitation (Weyl-Graviton 4 T.T/s^2, Skalar 4 (F.F)(F.F)/s^2, Doppelpol), nicht Einstein.
Teil 4  Bilanz und Urteil; Folgen für K15 (Wechselwirkung) und K13 (Gravitation).
"""
import functools
import itertools
import json
import time
from collections import Counter
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
import scipy.sparse as sps
import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
rng = np.random.default_rng(124)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def gs(x, n=6):
    return float(f"{float(x):.{n}g}")


def ang(a, b):
    return a[0]*b[1] - a[1]*b[0]


# ------------------------------------------------------------------ Grassmann-Algebra (wie z110)
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


def gprod(faktoren):
    return functools.reduce(gmul, faktoren, {0: 1})


def oben(A, ngen):
    return A.get((1 << ngen) - 1, 0)


# ------------------------------------------------------------------ Bausteine Stromalgebra (wie z110)
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


def pt_z(zs, ordnung):
    v = Fr(1)
    for a, b in glieder(list(ordnung)):
        v /= (zs[a] - zs[b])
    return v


def produkt(Ms):
    return functools.reduce(lambda A, B: A @ B, Ms)


def spur_i(M):
    return int(np.round(np.trace(M).real))


def ward(ops, form):
    """<J^X1(z1) ... J^Xn(zn)> aus der Wardidentität allein (OPE k form/(z-w)^2 + J^[X,Y]/(z-w)), Polynom in k."""
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


def wick(ops):
    """Fermionschleifen (Dirac): Summe über fixpunktfreie Permutationen, je Schleife -Tr(Wort) mal Parke-Taylor-Kette."""
    n = len(ops)
    zs = [z for _, z in ops]
    out = {}
    for perm in derangements(n):
        cs = zyklen(perm)
        val = Fr(1)
        for c in cs:
            val *= -spur_i(produkt([ops[i][0] for i in c]))*pt_z(zs, c)
        out[len(cs)] = out.get(len(cs), 0) + val
    return {p: c for p, c in out.items() if c != 0}


FORM = lambda A, B: spur_i(A @ B)
EP = np.array([[0, 1], [0, 0]], dtype=np.int64)
EM = np.array([[0, 0], [1, 0]], dtype=np.int64)
HH = np.array([[1, 0], [0, -1]], dtype=np.int64)


def su2_zufall():
    a, b, c = (int(x) for x in rng.integers(-6, 7, size=3))
    return a*EP + b*EM + c*HH


def gl8_spurfrei():
    X = rng.integers(-2, 3, size=(8, 8)).astype(np.int64)
    X[7, 7] -= np.trace(X)
    return X


def punkte(n):
    zs = set()
    while len(zs) < n:
        zs.add(Fr(int(rng.integers(-40, 41)), int(rng.integers(1, 9))))
    out = list(zs)
    rng.shuffle(out)
    return out


# ------------------------------------------------------------------ Kinematik (Vergleichsziele, wie z111/z112)
SIG = [np.eye(2, dtype=complex), np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]),
       np.array([[1, 0], [0, -1]], complex)]
ETA = np.diag([1., -1., -1., -1.])


def zk():
    return rng.normal(size=2) + 1j*rng.normal(size=2)


def kin_komplex(n):
    lam = rng.normal(size=(n, 2)) + 1j*rng.normal(size=(n, 2))
    lt = rng.normal(size=(n, 2)) + 1j*rng.normal(size=(n, 2))
    rest = sum(np.outer(lam[i], lt[i]) for i in range(2, n))
    T = -np.linalg.solve(np.array([lam[0], lam[1]]).T, rest)
    lt[0], lt[1] = T[0], T[1]
    return lam, lt


def vierer(M):
    return np.array([0.5*np.trace(SIG[m] @ M) for m in range(4)])


def feldstaerke(lam, lt, h, mu):
    """Linearisierte Feldstärke F^(mu nu) = p^mu e^nu - p^nu e^mu einer Photon-/Gluonwelle der Helizität h (z112-Konvention)."""
    p = vierer(np.outer(lam, lt))
    e = vierer(np.outer(lam, mu)/ang(lt, mu)) if h < 0 else vierer(np.outer(mu, lt)/ang(mu, lam))
    return np.outer(p, e) - np.outer(e, p)


def fdot(Fa, Fb):
    return np.einsum("mn,mn->", Fa, ETA @ Fb @ ETA)


def spannung(Fa, Fb):
    """Maxwell-Spannungstensor (Matrixelement, symmetrisiert): T = F_a F_b^T eta + (a <-> b) - eta^-1 (F_a.F_b)/2."""
    t = Fa @ ETA @ Fb.T + Fb @ ETA @ Fa.T
    return t - 0.5*np.linalg.inv(ETA)*fdot(Fa, Fb)


def tdot(Ta, Tb):
    return np.einsum("mn,mn->", Ta, ETA @ Tb @ ETA)


# =================================================================== Teil 1
say("Teil 1: Gradzwist statt Holonomiezwist")
say("  (a) Kohomologie von O(d) auf CP^1 (Cech, Laurentmonome)")


def cech(d, K):
    """h^0, h^1 von O(d): C^0 = Monome z^k (k = 0..K) über U_0 und z^(d-j) (j = 0..K) über U_inf, C^1 = Laurentmonome."""
    es = list(range(min(0, d - K), max(K, d) + 1))
    ix = {e: i for i, e in enumerate(es)}
    spalten = []
    for k in range(K + 1):
        v = np.zeros(len(es))
        v[ix[k]] = 1.0
        spalten.append(v)
    for j in range(K + 1):
        v = np.zeros(len(es))
        v[ix[d - j]] = -1.0
        spalten.append(v)
    D = np.array(spalten).T
    rk = np.linalg.matrix_rank(D)
    return D.shape[1] - rk, D.shape[0] - rk


coh = {d: tuple(int(v) for v in cech(d, 8)) for d in range(-4, 4)}
coh_stabil = all(tuple(int(v) for v in cech(d, 13)) == coh[d] for d in coh)
rr = all(h0 - h1 == d + 1 for d, (h0, h1) in coh.items())
serre = all(coh[d][1] == coh[-2 - d][0] for d in coh if -2 - d in coh)
azyklisch = [d for d, (h0, h1) in coh.items() if h0 == 0 and h1 == 0]
ZAHLEN["kohomologie_O(d)"] = {str(d): list(v) for d, v in coh.items()}
check("[E] (a) Cech-Kohomologie von O(d), d = -4..3 (Laurentmonome, Rang der Differenz; stabil bei größerem Abschneider: "
      f"{coh_stabil}): (h^0, h^1) = {ZAHLEN['kohomologie_O(d)']}; Riemann-Roch h^0 - h^1 = d + 1 ({rr}), Serre "
      f"h^1(O(d)) = h^0(O(-2-d)) ({serre}); azyklisch nur d = {azyklisch}. Die nativen Nahtfermionen (Gewicht 1/2, K^(1/2) = "
      "O(-1)) haben keine Kohomologie, daher den eindeutigen Szegő-Kern und den Parke-Taylor-Nenner; ihr Gradzwist nach "
      "O(1) = K^(1/2) (x) O(2) hätte h^0 = 2, sein Partner O(-3) h^1 = 2: genau die Nullmoden theta^alpha lambda_alpha",
      coh_stabil and rr and serre and azyklisch == [-1] and coh[1] == (2, 0) and coh[-3] == (0, 2))

say("  (b) Nullmoden: Holonomiezwist (Quarterrahmen) gegen Gradzwist q (L_0 -+ q Q über der nativen Fermikante)")


def grundladungen(alpha, q, seite, M=30, Qmax=6):
    """Ladungen der Grundzustände von L_0 - seite q Q im Rahmen alpha (native Fermikante: x < 0 besetzt). seite = +1:
    Eingang (Marke 0), seite = -1: Ausgang (Marke unendlich, BPZ-Bild)."""
    xp = sorted(alpha + n for n in range(-M, M) if alpha + n >= 0)          # x = 0 (R-Rahmen): Füllen kostet nichts
    xn = sorted((alpha + n for n in range(-M, M) if alpha + n < 0), reverse=True)

    def E(Q):
        return sum(xp[:Q]) if Q >= 0 else sum(-x for x in xn[:-Q])
    werte = {Q: E(Q) - seite*q*Q for Q in range(-Qmax, Qmax + 1)}
    m = min(werte.values())
    return [Q for Q, v in werte.items() if v == m], xp, xn


def nullmoden(alpha, q):
    """Moden, deren Besetzung zwischen Eingangs- und Ausgangsgrundzustand wechselt (Nullmoden des Pfadintegrals)."""
    qe, xp, xn = grundladungen(alpha, q, +1)
    qa, _, _ = grundladungen(alpha, q, -1)
    if len(qe) != 1 or len(qa) != 1:
        return None, qe, qa

    def besetzt(Q):
        s = set(x for x in xn)
        if Q >= 0:
            s |= set(xp[:Q])
        else:
            s -= set(xn[:-Q])
        return s
    return sorted(besetzt(qe[0]) ^ besetzt(qa[0])), qe, qa


RAHMEN = [Fr(1, 2), Fr(1, 4), Fr(3, 4)]
nm_tab, nm_ok, wellen = {}, True, {}
for q in (0, 1):
    for al in RAHMEN:
        modes, qe, qa = nullmoden(al, Fr(q))
        lam_strich = Fr(1, 2) - q                                           # Gewicht nach dem Gradzwist
        expo = sorted(-x - lam_strich for x in modes) if modes is not None else None
        nm_tab[f"q = {q}, alpha = {al}"] = {"Q_ein, Q_aus": [qe, qa], "Nullmoden x": [str(x) for x in modes or []],
                                            "Exponenten z^e": [str(e) for e in expo or []]}
        if q == 0:
            nm_ok &= modes == [] and qe == [0] and qa == [0]
        else:
            nm_ok &= modes is not None and len(modes) == 2 and qe == [1] and qa == [-1]
            nm_ok &= expo is not None and expo[1] - expo[0] == 1 and expo[0] == Fr(1, 2) - al
            wellen[al] = [float(e) for e in expo]
ram = grundladungen(Fr(0), Fr(0), +1)[0]
ZAHLEN["nullmoden"] = nm_tab


def sl2_rest(expo, lam_strich, g, n=40):
    """Abgeschlossenheit von span{z^e} unter f -> (cz + d)^(-2 lam') f((az + b)/(cz + d)) (Hauptzweig, g nahe 1)."""
    a, b, c, d = g
    zs = 0.8 + 0.6*rng.random(n) + 1j*(rng.random(n) - 0.5)
    B = np.array([zs**e for e in expo]).T
    rest = 0.0
    for e in expo:
        f = (c*zs + d)**(-2*lam_strich)*((a*zs + b)/(c*zs + d))**e
        co = np.linalg.lstsq(B, f, rcond=None)[0]
        rest = max(rest, np.linalg.norm(B @ co - f)/np.linalg.norm(f))
    return rest


def sl2_nahe_eins(eps=0.25):
    X = eps*(rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2)))
    X -= np.trace(X)/2*np.eye(2)
    w, V = np.linalg.eig(X)
    g = V @ np.diag(np.exp(w)) @ np.linalg.inv(V)
    g = g/np.sqrt(np.linalg.det(g))
    return g[0, 0], g[0, 1], g[1, 0], g[1, 1]


t_diag = 1.3 - 0.4j
sl2 = {}
for al, ex in wellen.items():
    voll = max(sl2_rest(ex, -0.5, sl2_nahe_eins()) for _ in range(3))
    diag = sl2_rest(ex, -0.5, (t_diag, 0, 0, 1/t_diag))
    sl2[str(al)] = {"SL(2,C) (drei g)": gs(voll, 3), "Diagonalgruppe": gs(diag, 3)}
ZAHLEN["sl2_abgeschlossenheit"] = sl2
check("[E] (b) Nullmoden aus Eingangs- und Ausgangsgrundzustand von L_0 -+ q Q über der nativen Fermikante (Moden x in "
      f"Z + alpha): Ohne Gradzwist (q = 0) sind beide Grundzustände das Fermimeer (Q = 0) in den Rahmen alpha = 1/2, 1/4, "
      f"3/4: keine Nullmode, Index 0 (im R-Rahmen alpha = 0 ist der Grundzustand entartet, Q in {ram}: das Ramond-Paar aus z110, "
      "Index 0). Spin-Z4-Strukturen und Quarterrahmen sind genau diese flachen Holonomien alpha in Z/4: Grad 0, keine "
      "Nullmoden. Mit Gradzwist q = 1 (Gewicht -1/2) sind die Grundzustände eindeutig Q = +1 (Eingang) und Q = -1 (Ausgang): "
      f"in jedem Rahmen genau zwei Nullmoden ({ {k: v['Nullmoden x'] for k, v in nm_tab.items() if k.startswith('q = 1')} }), "
      "Index 2 = deg O(1) + 1, topologisch. Ihre Wellenfunktionen sind z^(1/2 - alpha) (1, z) "
      f"({ {k: v['Exponenten z^e'] for k, v in nm_tab.items() if k.startswith('q = 1')} }): das lambda-Dublett nur bei "
      f"alpha = 1/2. SL(2,C)-Abgeschlossenheit des Nullmodenraums (Restfehler {sl2}): exakt nur ungetwistet, im "
      "Quarterrahmen nur unter der Diagonalgruppe (Drehung und Boost um die Markenachse)",
      nm_ok and len(ram) == 2 and sl2["1/2"]["SL(2,C) (drei g)"] < 1e-9 and sl2["1/4"]["SL(2,C) (drei g)"] > 1e-3
      and sl2["3/4"]["SL(2,C) (drei g)"] > 1e-3 and all(v["Diagonalgruppe"] < 1e-9 for v in sl2.values()))

say("  (c) Berezin-Zähler mit den Gradzwist-Nullmoden (vier Felder A = 1..4, zwei eta^4-Einsetzungen)")


def zaehler_tab(wf, lams, n_ein):
    """Für alle Lagen von eta^4 (Gluonen: eta^0 oder eta^4 je Punkt): int d^8 theta prod_(eta^4-Punkte) prod_A psi^A,
    psi^A(i) = theta^(A,0) wf0(i) + theta^(A,1) wf1(i). Rückgabe {Lage: Wert}."""
    out = {}
    for lage in itertools.product((0, 1), repeat=n_ein):
        fak = []
        for i, s in enumerate(lage):
            if s:
                w0, w1 = wf(lams[i])
                for A in range(4):
                    fak.append({1 << (2*A): w0, 1 << (2*A + 1): w1})
        out[lage] = oben(gprod(fak), 8)
    return out


lq = [np.array([1.0, z]) for z in (0.7 + 0.2j, 1.3 - 0.5j, 0.9 + 0.9j, 1.6 + 0.1j)]
wf_ns = lambda l: (l[0], l[1])                                           # O(1): theta^alpha lambda_alpha
wf_q = lambda l: (l[0]*(l[1]/l[0])**0.25, l[1]*(l[1]/l[0])**0.25)        # alpha = 1/4: z^(1/4) (1, z), homogen
wf_q3 = lambda l: (l[0]*(l[1]/l[0])**-0.25, l[1]*(l[1]/l[0])**-0.25)     # alpha = 3/4: z^(-1/4) (1, z)
tab_ns = zaehler_tab(wf_ns, lq, 4)
tab_q = zaehler_tab(wf_q, lq, 4)
tab_q3 = zaehler_tab(wf_q3, lq, 4)
nz_ns = sorted(sum(l) for l, v in tab_ns.items() if abs(v) > 1e-12)
nz_q = sorted(sum(l) for l, v in tab_q.items() if abs(v) > 1e-12)
mhv_ns = max(abs(tab_ns[l]/ang(lq[i], lq[j])**4 - 1) for l in tab_ns if sum(l) == 2
             for i, j in [tuple(k for k in range(4) if l[k])])
quot_q = [complex(tab_q[l]/ang(lq[i], lq[j])**4) for l in tab_q if sum(l) == 2 for i, j in [tuple(k for k in range(4) if l[k])]]
zz = [complex(lq[i][1]*lq[j][1]) for l in tab_q if sum(l) == 2 for i, j in [tuple(k for k in range(4) if l[k])]]
quot_ok = max(abs(a_ - b_) for a_, b_ in zip(quot_q, zz))
quot3 = max(abs(tab_q3[l]/ang(lq[i], lq[j])**4*lq[i][1]*lq[j][1] - 1) for l in tab_q3 if sum(l) == 2
            for i, j in [tuple(k for k in range(4) if l[k])])
Mg = np.array(sl2_nahe_eins(0.6)).reshape(2, 2)
lq_g = [Mg @ l for l in lq]
lor_q = abs(zaehler_tab(wf_q, lq_g, 4)[(1, 1, 0, 0)]/tab_q[(1, 1, 0, 0)] - 1)
lor_ns = abs(zaehler_tab(wf_ns, lq_g, 4)[(1, 1, 0, 0)]/tab_ns[(1, 1, 0, 0)] - 1)
Dg = np.diag([t_diag, 1/t_diag])
char_q = complex(zaehler_tab(wf_q, [Dg @ l for l in lq], 4)[(1, 1, 0, 0)]/tab_q[(1, 1, 0, 0)])
ZAHLEN["berezin_gradzwist"] = {"Lagen != 0 (Zahl der eta^4), ungetwistet": nz_ns, "Quarterrahmen": nz_q,
                               "max |Zähler/<ij>^4 - 1| ungetwistet": gs(mhv_ns, 3),
                               "Quarter: Zähler/<ij>^4 - z_i z_j": gs(quot_ok, 3),
                               "alpha = 3/4: |Zähler z_i z_j/<ij>^4 - 1|": gs(quot3, 3),
                               "SL(2,C)-Defekt ungetwistet / Quarter": [gs(lor_ns, 3), gs(lor_q, 3)],
                               "Quarter unter diag(t, 1/t): Faktor, t^-4": [str(np.round(char_q, 9)),
                                                                           str(np.round(t_diag**-4, 9))]}
check("[E] (c) Berezin-Integral über die acht Gradzwist-Nullmoden (vier Felder, je zwei): In beiden Rahmen tragen genau die "
      f"Lagen mit zwei eta^4-Einsetzungen (ungetwistet {Counter(nz_ns)}, Quarter {Counter(nz_q)}; je 6 von 16): die "
      f"MHV-Auswahlregel folgt aus dem Index allein. Ungetwistet ist der Zähler <ij>^4 (Abweichung {mhv_ns:.0e}) und "
      f"SL(2,C)-invariant ({lor_ns:.0e}). Im Quarterrahmen ist er <ij>^4 z_i z_j (Abweichung {quot_ok:.0e}; bei alpha = 3/4 "
      f"<ij>^4/(z_i z_j), {quot3:.0e}): nicht lorentzinvariant (relativ {lor_q:.2f}), unter diag(t, 1/t) nur kovariant "
      f"mit dem Charakter t^-4 ({np.round(char_q, 6)} gegen {np.round(t_diag**-4, 6)}): eine Drehimpulsladung entlang "
      "der Markenachse. Der "
      "Quarterzwist und der Gradzwist sind unverträglich: <ij>^4 verlangt den ungetwisteten Rahmen",
      nz_ns == [2]*6 and nz_q == [2]*6 and mhv_ns < 1e-9 and lor_ns < 1e-9 and quot_ok < 1e-9 and quot3 < 1e-9
      and lor_q > 1e-3 and abs(char_q - t_diag**-4) < 1e-9)

say("  (d) Gradzwist = topologische Verdrillung T -> T + a dJ (Virasoro mit U(1)-Strom, symbolisch)")
m_, n_, r_, a_, c_, k_, h_ = sp.symbols("m n r a c k h")
# [L_m, L_n] = (m - n) L_(m+n) + c/12 (m^3 - m) delta, [L_m, J_n] = -n J_(m+n), [J_m, J_n] = k m delta; L'_m = L_m - a (m+1) J_m
J_koeff = sp.expand(-a_*(n_ + 1)*(-n_) - a_*(m_ + 1)*m_)          # J_(m+n)-Koeffizient von [L'_m, L'_n]
J_soll = sp.expand(-a_*(m_ - n_)*(m_ + n_ + 1))                   # aus (m - n) L'_(m+n)
zentral = sp.expand((c_/12*(m_**3 - m_) + a_**2*(m_ + 1)*(-m_ + 1)*k_*m_))
c_strich = sp.simplify(12*zentral/(m_**3 - m_))
psi_koeff = sp.expand(((h_ - 1)*m_ - r_) - a_*(m_ + 1))             # [L'_m, psi_r] für Ladung +1
psi_soll = sp.expand(((h_ - a_) - 1)*m_ - (r_ + a_))
psib_koeff = sp.expand(((h_ - 1)*m_ - r_) + a_*(m_ + 1))            # Ladung -1
psib_soll = sp.expand(((h_ + a_) - 1)*m_ - (r_ - a_))
werte_c = {str(av): int(c_strich.subs({c_: 1, k_: 1, a_: av})) for av in (0, sp.Rational(1, 2), 1)}
ZAHLEN["verdrillung"] = {"c'": str(c_strich), "c' für c = k = 1": werte_c,
                         "Gewichte (psi, psi*) bei a = 1": ["-1/2", "3/2"], "bei a = 1/2": ["0", "1"]}
check("[E] (d) Der Gradzwist ist die topologische Verdrillung L'_m = L_m - a (m + 1) J_m (T' = T + a dJ): symbolisch ist "
      f"[L'_m, L'_n] = (m - n) L'_(m+n) + c'/12 (m^3 - m) delta mit c' = {c_strich} (J-Rest {sp.simplify(J_koeff - J_soll)}); "
      "ein Primärfeld der Ladung +-1 erhält das Gewicht h -+ a bei umbenannten Moden r -+ a. Für eine komplexe Nahtfermion "
      f"(c = k = 1): c' = {werte_c}; a = 1 gibt (psi, psi*) = (O(1), O(-3)) mit c' = -11, a = 1/2 die Twistorstring-Fermionen "
      "(Gewicht 0, 1; c = -2), deren Homogenität 1 aus dem Grad-1-Bündel der GL(1)-Eichung kommt. Allgemein hat "
      "K^(1/2 - a) (x) O(d) den Grad -1 + 2a + d; O(1) verlangt 2a + d = 2, also a != 0 oder ein Bündel vom Grad d != 0",
      sp.simplify(J_koeff - J_soll) == 0 and sp.simplify(c_strich - (c_ - 12*a_**2*k_)) == 0
      and sp.simplify(psi_koeff - psi_soll) == 0 and sp.simplify(psib_koeff - psib_soll) == 0
      and werte_c == {"0": 1, "1/2": -2, "1": -11})

say("  (e) Native Kandidaten auf dem Fockraum des Originalpaars (1/4, 3/4): Uhr R, Teilchen-Loch-J, Spin-Z4")
KF = 4
MODEN = [(i, al + n) for i, al in enumerate((Fr(1, 4), Fr(3, 4))) for n in range(-KF, KF)]
NB = len(MODEN)
ZUST = np.arange(1 << NB)
BES = ((ZUST[:, None] >> np.arange(NB)[None, :]) & 1).astype(np.int64)
XS = np.array([float(x) for _, x in MODEN])
NAHT = np.array([i for i, _ in MODEN])
ANGEREGT = np.where(XS[None, :] > 0, BES, 1 - BES)
L0 = (ANGEREGT*np.abs(XS)[None, :]).sum(1)
QN = [(np.where(XS[None, :] > 0, BES, 0) - np.where(XS[None, :] < 0, 1 - BES, 0))[:, NAHT == i].sum(1) for i in range(2)]
R_uhr = np.exp(2j*np.pi*L0)
uhr_def = np.abs(R_uhr - 1j**((QN[0] - QN[1]) % 4)).max()
QT = QN[0] + QN[1]
verdr = {str(qv): gs(np.abs(np.exp(2j*np.pi*(L0 - float(qv)*QT)) - R_uhr).max(), 2) for qv in (1, 2, -1)}
halb = float(np.abs(np.exp(2j*np.pi*(L0 - 0.5*QT)) - R_uhr.conj()).max())
bild = {j: MODEN.index((1 - i, -x)) for j, (i, x) in enumerate(MODEN)}
CZ = np.zeros(len(ZUST), dtype=np.int64)
for j in range(NB):
    CZ |= (1 - BES[:, j]) << bild[j]
ph_L0 = float(np.abs(L0[CZ] - L0).max())
ph_Q = float(np.abs(QT[CZ] + QT).max())
ph_tw = {str(qv): float(np.abs((L0 - qv*QT)[CZ] - (L0 + qv*QT)).max()) for qv in (1, 0.5)}
ZAHLEN["native_kandidaten_fock"] = {"Fockraum": f"2 x {2*KF} Moden, {len(ZUST)} Zustände",
                                    "|e^(2 pi i L_0) - i^(Q_a - Q_b)|": gs(uhr_def, 3),
                                    "|e^(2 pi i (L_0 - q Q)) - R| (q ganz)": verdr,
                                    "|e^(2 pi i (L_0 - Q/2)) - R^-1|": gs(halb, 3),
                                    "J: L_0 invariant, Q -> -Q": [ph_L0, ph_Q], "J: L_0 - q Q -> L_0 + q Q": ph_tw}
check("[E] (e) Uhr: Auf dem Fockraum des Originalpaars (alpha = 1/4, 3/4, je acht Moden, L_0 über der Fermikante) gilt "
      f"e^(2 pi i L_0) = i^(Q_a - Q_b) = R (Abweichung {uhr_def:.0e}, wie z113). Die Uhr ist blind für ganzzahlige "
      f"Gradzwiste: e^(2 pi i (L_0 - q Q)) = R für q = 1, 2, -1 ({verdr}); der halbe Zwist bildet sie auf R^-1 ab "
      f"({halb:.0e}). R sieht nur L_0 mod 1, liefert also weder den Zwist noch wählt sie ihn. J (Teilchen-Loch-Anteil, "
      f"Naht a <-> b, x -> -x): L_0 invariant ({ph_L0:.0e}), Q -> -Q ({ph_Q:.0e}), also L_0 - q Q -> L_0 + q Q ({ph_tw}): J "
      "kehrt den Zwist um; ein J-kovarianter Zwist braucht eine Orientierungswahl (ein Bit wie P_vol, z103, oder P_RJ, z114). "
      "Spin-Z4 und Quarterrahmen sind flache Holonomien in Z/4 (Teil 1b: Index 0); auf CP^1 ohne Marken ist jedes Z4-Bündel "
      "trivial (pi_1(CP^1) = 0), und nach z100 sind die nativen Verbindungen flach: kein Grad",
      uhr_def < 1e-9 and max(verdr.values()) < 1e-9 and halb < 1e-9 and ph_L0 < 1e-12 and ph_Q == 0
      and max(ph_tw.values()) < 1e-12)

say("  (f) Ladungsanomalie: Mit Gradzwist verschwinden neutrale Stromkorrelatoren")


def anomalie(nullmoden_da, einsetz, mz=3):
    """Gaußsches Grassmann-Modell: Nullmoden theta_0, theta_1 (nur psi), Nicht-Nullmoden a_r (psi), b_r (psi*) mit Kinetik
    exp(b K a). Felder an Punkten: psi(z) = [theta_0 + theta_1 z] + sum f_r(z) a_r, psi*(z) = sum g_r(z) b_r."""
    nz = 2 if nullmoden_da else 0
    ngen = nz + 2*mz
    K = rng.normal(size=(mz, mz)) + 1j*rng.normal(size=(mz, mz))
    expo = {0: 1.0}
    for r1 in range(mz):
        for s1 in range(mz):
            expo = gmul(expo, {0: 1.0, (1 << (nz + mz + r1)) | (1 << (nz + s1)): K[r1, s1]})
    fs = rng.normal(size=(8, mz))
    gs_ = rng.normal(size=(8, mz))

    def psi(i, z):
        out = {}
        if nullmoden_da:
            out[1 << 0] = 1.0
            out[1 << 1] = z
        for r1 in range(mz):
            out[1 << (nz + r1)] = fs[i, r1]
        return out

    def psib(i):
        return {1 << (nz + mz + r1): gs_[i, r1] for r1 in range(mz)}
    fak = [expo]
    zs = [0.3, -0.7, 1.1, 0.5, -1.4, 0.9]
    for i, art in enumerate(einsetz):
        if art == "J":
            fak += [psib(i), psi(i, zs[i])]
        elif art == "psi":
            fak.append(psi(i, zs[i]))
    return oben(gprod(fak), ngen)


anom = {"ungetwistet, <J J>": anomalie(False, ["J", "J"]), "ungetwistet, <J J J>": anomalie(False, ["J", "J", "J"]),
        "Gradzwist, <J J>": anomalie(True, ["J", "J"]), "Gradzwist, <J J J>": anomalie(True, ["J", "J", "J"]),
        "Gradzwist, <J psi psi>": anomalie(True, ["J", "psi", "psi"]),
        "Gradzwist, <J J psi psi>": anomalie(True, ["J", "J", "psi", "psi"])}
ZAHLEN["ladungsanomalie"] = {k: gs(abs(v), 4) for k, v in anom.items()}
check("[E] (f) Ladungsanomalie im gaußschen Grassmann-Modell (zwei psi-Nullmoden des Gradzwists, drei Paare Nicht-Nullmoden "
      f"mit zufälliger Kinetik): |Korrelatoren| {ZAHLEN['ladungsanomalie']}. Ohne Nullmoden sind die neutralen Stromkorrelatoren "
      "<J J>, <J J J> nicht null (Nenner), mit Gradzwist sind sie exakt null; nur Einsetzungen mit Nettoladung 2 (zwei psi) "
      "überleben. Dieselben Fermionen können daher nicht zugleich den Parke-Taylor-Nenner (neutrale Ströme) und den Zähler "
      "(Nullmoden) tragen: P_N4 verlangt vier zusätzliche, verdrillte Richtungen neben den Stromfermionen",
      abs(anom["ungetwistet, <J J>"]) > 1e-6 and abs(anom["ungetwistet, <J J J>"]) > 1e-6
      and anom["Gradzwist, <J J>"] == 0 and anom["Gradzwist, <J J J>"] == 0 and abs(anom["Gradzwist, <J psi psi>"]) > 1e-6
      and abs(anom["Gradzwist, <J J psi psi>"]) > 1e-6)
befund("Teil 1: Der Gradzwist ist mathematisch genau der richtige Mechanismus. Die Kopplung der U(1)-Ladung an "
       "K^-1 = O(2) (das Inverse des Quadrats der tautologischen Linie, gleichwertig die topologische Verdrillung "
       "T + dJ) macht aus den Nahtfermionen O(1)-Felder mit Index 2 und den Nullmoden theta^alpha lambda_alpha; c fällt "
       "je komplexer Fermion von 1 auf -11. Nativ liefert ihn nichts. Uhr, Quarterrahmen und Spin-Z4 sind flache "
       "Z/4-Holonomien vom Grad 0; die Uhr sieht nur L_0 mod 1. J kehrt den Zwist um, und die nativen Verbindungen sind "
       "flach (z100). Zudem ist der Quarterrahmen mit ihm unverträglich (Zähler <ij>^4 z_i z_j, nur diagonal kovariant). "
       "Die verdrillten Fermionen töten nach der Ladungsanomalie die neutralen Stromkorrelatoren, also den Nenner. Die "
       "z64/z70-Konstruktion O(1) + O(1) = T_C (x) pi_* theta ist derselbe Zwist mit T_C = O(2) auf Bündelebene, nicht als "
       "Dynamik der Naht", True)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

# =================================================================== Teil 2
say("\nTeil 2: Pati-Salam-Quartett als eta^A")
say("  (a) 16 = (4,2,1) + (4-quer,1,2) aus der Cliffordrealisierung (fünf Moden, z116) und das Supermultiplett Lambda(C^4)")
s0, sx_, sy_, sz_ = SIG[0], SIG[1], SIG[2], SIG[3]


def kr(*m):
    out = m[0]
    for x in m[1:]:
        out = np.kron(out, x)
    return out


GI = []
for j in range(5):
    GI.append(kr(*([sz_]*j + [sx_] + [s0]*(4 - j))))
    GI.append(kr(*([sz_]*j + [sy_] + [s0]*(4 - j))))
G11 = np.eye(32, dtype=complex)
for g in GI:
    G11 = G11 @ g
for ph in (1, 1j, -1, -1j):
    if np.allclose((ph*G11) @ (ph*G11), np.eye(32)):
        G11 = ph*G11
        break
A5 = [(GI[2*j] + 1j*GI[2*j + 1])/2 for j in range(5)]                    # a_j wie z116
ND = [a.conj().T @ a for a in A5]
PAR = functools.reduce(lambda A, B: A @ B, [np.eye(32) - 2*n for n in ND])
P16 = (np.eye(32) + G11)/2
if not np.allclose(P16, (np.eye(32) + PAR)/2):
    P16 = (np.eye(32) - G11)/2
Mgen = lambda a, b: 0.5j*GI[a] @ GI[b]
C6 = sum(Mgen(a, b) @ Mgen(a, b) for a in range(6) for b in range(a + 1, 6))
G6 = functools.reduce(lambda A, B: A @ B, GI[:6])
for ph in (1, 1j, -1, -1j):
    if np.allclose((ph*G6) @ (ph*G6), np.eye(32)):
        G6 = ph*G6
        break
TA = [A5[3].conj().T @ A5[4], A5[4].conj().T @ A5[3], (ND[3] - ND[4])/2]            # SU(2)_L: Einteilchensektor der Moden 3, 4
TB = [A5[3].conj().T @ A5[4].conj().T, A5[4] @ A5[3], (ND[3] + ND[4] - np.eye(32))/2]  # SU(2)_R: Sektor {0, 2}
cas = lambda T: T[0] @ T[1]/2 + T[1] @ T[0]/2 + T[2] @ T[2]
CA, CB = cas(TA), cas(TB)


def dim(*P):
    return int(round(np.trace(functools.reduce(lambda A, B: A @ B, P)).real))


PL = CA/0.75
PR = CB/0.75
P6p, P6m = (np.eye(32) + G6)/2, (np.eye(32) - G6)/2
spek_ok = (np.allclose(PL @ PL, PL) and np.allclose(PR @ PR, PR) and np.allclose(C6 @ P16, (15/4)*P16)
           and np.allclose(G6 @ C6, C6 @ G6))
zerl = {"(4,2,1)": dim(P16, P6p, PL), "(4-quer,2,1)": dim(P16, P6m, PL), "(4,1,2)": dim(P16, P6p, PR),
        "(4-quer,1,2)": dim(P16, P6m, PR), "(.,1,1) oder (.,2,2)": 16 - dim(P16, PL) - dim(P16, PR)}
if zerl["(4,2,1)"] == 0:                                                 # Benennung 4 <-> 4-quer ist Konvention
    zerl = {"(4,2,1)": zerl["(4-quer,2,1)"], "(4-quer,2,1)": 0, "(4,1,2)": zerl["(4-quer,1,2)"],
            "(4-quer,1,2)": zerl["(4,1,2)"], "(.,1,1) oder (.,2,2)": zerl["(.,1,1) oder (.,2,2)"]}
N4 = sum(ND[:4])
U4 = [A5[i].conj().T @ A5[j] for i in range(4) for j in range(4)]
CU4 = sum(U4[4*i + j] @ U4[4*j + i] for i in range(4) for j in range(4)) - N4 @ N4/4
ew_ps = Counter(float(v) for v in np.round(np.linalg.eigvalsh(P16 @ C6 @ P16 + 7*(np.eye(32) - P16)), 6))
ew_su5 = Counter(float(v) for v in np.round(np.linalg.eigvalsh(P16 @ CU4 @ P16 + 7*(np.eye(32) - P16)), 6))
ew_ps.pop(7.0, None)
ew_su5.pop(7.0, None)
schichten = Counter(int(round(v)) for v in np.diag(P16 @ N4).real[np.abs(np.diag(P16)) > 0.5])
ZAHLEN["pati_salam"] = {"Zerlegung der 16": zerl, "Casimir SU(4)_PS auf der 16": {str(k): v for k, v in ew_ps.items()},
                        "Casimir SU(4) c SU(5) auf der 16": {str(k): v for k, v in ew_su5.items()},
                        "Lambda^k(C^4)-Schichten unter SU(4) c SU(5)": dict(sorted(schichten.items()))}
check("[E] (a) Cliffordrealisierung der 16 (Lambda^gerade C^5, z116): Spin(6) aus Gamma_0..5 (Moden 0, 1, 2), Spin(4) = "
      "SU(2)_L x SU(2)_R aus den Moden 3, 4 (Einteilchen- bzw. {0, 2}-Sektor): " + f"{zerl}" + ", also 16 = (4,2,1) + "
      "(4-quer,1,2), mit der Spin(6)-Chiralität an die SU(2)-Seite gekoppelt. Dieselbe Zerlegung ist die der 16 Superladungen "
      "von N = 4 (Q^(alpha A) + Q-quer_(alpha-punkt A)) unter SU(4)_R x Lorentz, wenn Spin(4) die Lorentzgruppe ist "
      "(Reduktion der 10D-Majorana-Weyl-16): (4,2,1) hat genau die acht Komponenten theta^(A alpha). Das Supermultiplett "
      "Lambda(eta) = 1 + 4 + 6 + 4-quer + 1 steckt in der 16 nur unter SU(4) c SU(5) (Moden 0..3; Schichten "
      f"{dict(sorted(schichten.items()))}, Casimirwerte {dict(ew_su5)}), nicht unter SU(4)_PS (ein einziger "
      f"Casimirwert {dict(ew_ps)}: nur 4 und 4-quer)",
      spek_ok and zerl["(4,2,1)"] == 8 and zerl["(4-quer,1,2)"] == 8 and zerl["(4-quer,2,1)"] == 0
      and zerl["(4,1,2)"] == 0 and zerl["(.,1,1) oder (.,2,2)"] == 0 and len(ew_ps) == 1
      and sorted(ew_su5.values()) == [2, 6, 8] and dict(sorted(schichten.items())) == {0: 1, 1: 4, 2: 6, 3: 4, 4: 1})

say("  (b) Statistik: R^2 auf der 16 und bosonische Kontrolle")
R16 = 1j*G11
r2 = float(np.abs(R16 @ R16 + np.eye(32)).max())
lb = [zk() for _ in range(4)]
Gsym = rng.normal(size=(2, 2))
Gsym = Gsym @ Gsym.T + np.eye(2)
bos = lambda L, i, j: (L[i] @ Gsym @ L[j])**4
Mb = np.array(sl2_nahe_eins(0.6)).reshape(2, 2)
bos_def = abs(bos([Mb @ l for l in lb], 0, 1)/bos(lb, 0, 1) - 1)
ferm_def = abs(ang(Mb @ lb[0], Mb @ lb[1])**4/ang(lb[0], lb[1])**4 - 1)
ZAHLEN["statistik"] = {"|R^2 + 1| auf 32": gs(r2, 3), "bosonische Gaußpaarung: SL(2,C)-Defekt": gs(bos_def, 3),
                       "Grassmann <ij>^4: Defekt": gs(ferm_def, 3)}
check(f"[E] (b) Statistik: R = i Gamma_11 hat auf allen 32 Zuständen R^2 = -1 (Abweichung {r2:.0e}), alle Komponenten der 16 "
      "sind also R^2-ungerade (z106: die CAR gibt R^2-Statistik, unter P_CAR ungerade): die Grassmann-Parität wäre gestellt. "
      "Kontrolle, dass sie nötig ist: Mit kommutierenden theta und gaußscher (symmetrischer) Paarung wird der Zähler "
      f"(lambda_i G lambda_j)^4, nicht lorentzinvariant (relativ {bos_def:.2f}), während <ij>^4 invariant ist ({ferm_def:.0e})",
      r2 < 1e-12 and bos_def > 1e-3 and ferm_def < 1e-9)

say("  (c) Homogenität: Berezin-Zähler für die native Lage der Farbquartette")
u_fest = zk()
u_pkt = [zk() for _ in range(4)]
real = {"h = +1 (theta^(A alpha) lambda_alpha, verlangt)": lambda i, L: L[i],
        "h = 0, innere SU(2)_L-Polarisation u_i": lambda i, L: u_pkt[i],
        "h = 0, feste innere Polarisation u": lambda i, L: u_fest}


def zaehler_real(wf, L, lage):
    fak = []
    for i, s in enumerate(lage):
        if s:
            w = wf(i, L)
            for A in range(4):
                fak.append({1 << (2*A): w[0], 1 << (2*A + 1): w[1]})
    return oben(gprod(fak), 8)


hom = {}
for name, wf in real.items():
    v = zaehler_real(wf, lb, (1, 1, 0, 0))
    t_ = 1.7 - 0.3j
    Ls = [l.copy() for l in lb]
    Ls[0] = t_*Ls[0]
    gew = np.log(abs(zaehler_real(wf, Ls, (1, 1, 0, 0))/v))/np.log(abs(t_)) if abs(v) > 1e-14 else None
    hom[name] = {"|Zähler(1, 2)|": gs(abs(v), 4), "Gewicht in lambda_1": None if gew is None else gs(gew, 6),
                 "Helizität des eta^4-Beins (mit Nenner -2)": None if gew is None else gs(1 - gew/2, 6)}
hom["h = -1 (Nahtfermion O(-1))"] = {"Nullmoden": 0, "Zähler": "nur Allplus (PT delta^4(P)), jede eta^4-Lage 0"}
WURZELN = [tuple(Fr(s1) if k == i else Fr(s2) if k == j else Fr(0) for k in range(8))
           for i, j in itertools.combinations(range(8), 2) for s1 in (1, -1) for s2 in (1, -1)]
WURZELN += [tuple(Fr(x, 2) for x in sg) for sg in itertools.product((1, -1), repeat=8) if sg.count(-1) % 2 == 0]
e8_typ = Counter()
for b in WURZELN:
    vorn, hinten = b[:5], b[5:]
    if all(x.denominator == 2 for x in b):
        e8_typ["(16,4)" if sum(1 for x in vorn if x < 0) % 2 == 0 else "(16-quer,4-quer)"] += 1
    elif any(vorn) and any(hinten):
        e8_typ["(10,6)"] += 1
    elif any(vorn):
        e8_typ["(45,1)"] += 1
    else:
        e8_typ["(1,15)"] += 1
hom["h = -2 (E8-Ströme (16,4), Gewicht 1, bosonisch)"] = {
    "Wurzeln": dict(e8_typ), "Statistik": "bosonisch (Stromalgebra), Homogenität -2 wie jeder Strom (z110)"}
beres = {"h = +1": 4 - 1*4, "h = 0": 4 - 0*4, "h = -1": 4 + 4, "h = -2": 4 + 8}
ZAHLEN["homogenitaet"] = hom
ZAHLEN["berezin_gewicht 4 - h N (N = 4)"] = beres
check("[E] (c) Homogenität: Der Berezin-Zähler zweier eta^4-Einsetzungen hängt an der Lage der Quartett-Richtungen auf "
      f"der Geraden: {hom}. Nur h = +1 mit dem Dublettindex als lambda-Index gibt <ij>^4 (Gewicht 4, Helizität -1) und "
      f"ein Calabi-Yau-Maß (Gewicht 4 - hN: {beres}). Als innere Etiketten (h = 0) ist der Zähler lambda-unabhängig "
      "(Gewicht 0: das Bein bliebe Helizität +1) oder null, als Nahtfermionen (h = -1) gibt es keine Nullmoden. Auf der "
      "E8-Naht erscheint die 16 mit einem Vierer-Index nur als Strom (16,4) + (16-quer,4-quer) (SO(10) x SU(4)' c E8: "
      f"{dict(e8_typ)}; der Vierer ist das A3 aus z97/z102): bosonisch, Gewicht 1, Homogenität -2. Keine native Lage der "
      "Farbquartette hat Homogenität +1",
      abs(hom["h = +1 (theta^(A alpha) lambda_alpha, verlangt)"]["Gewicht in lambda_1"] - 4) < 1e-9
      and abs(hom["h = 0, innere SU(2)_L-Polarisation u_i"]["Gewicht in lambda_1"]) < 1e-9
      and hom["h = 0, feste innere Polarisation u"]["|Zähler(1, 2)|"] < 1e-12
      and dict(e8_typ) == {"(45,1)": 40, "(1,15)": 12, "(10,6)": 60, "(16,4)": 64, "(16-quer,4-quer)": 64}
      and len(WURZELN) == 240)

say("  (d) Die Gleichsetzung SU(2)_L = SL(2)_lambda ist eine topologische Verdrillung der Materie")
Jd = [np.kron(S/2, s0) + np.kron(s0, S/2) for S in SIG[1:]]
Cd = sum(J @ J for J in Jd)
spins = Counter(float(v) for v in np.round(np.linalg.eigvalsh(Cd), 9))
Yop = sum(y*n for y, n in zip([-1/3]*3 + [1/2]*2, ND))
BL = Yop - TB[2]
bl_werte = sorted(set(float(v) for v in np.round(np.diag(P16 @ BL).real[np.abs(np.diag(P16)) > 0.5], 9)))
Qem = TA[2] + Yop
komm_L = max(np.abs((P16 @ (Yop @ T - T @ Yop) @ P16)).max() for T in TA[:2])
komm_R = max(np.abs((P16 @ (Yop @ T - T @ Yop) @ P16)).max() for T in TB[:2])
q_rot = np.abs(P16 @ (Qem @ TA[0] - TA[0] @ Qem - TA[0]) @ P16).max()
ZAHLEN["verdrillung_materie"] = {"Casimir der Diagonal-SU(2) auf S_L (x) 2_L": {str(k): v for k, v in spins.items()},
                                 "(B-L)/2 = Y - T3_R auf der 16": [str(x) for x in bl_werte],
                                 "|[Y, SU(2)_L]|, |[Y, SU(2)_R]|": [gs(komm_L, 3), gs(komm_R, 3)],
                                 "|[Q, T_L^+] - T_L^+|": gs(q_rot, 3)}
check("[E] (d) Wäre der Dublettindex des Quartetts der lambda-Index, so wirkte auf die Materie die Diagonale von Lorentz-SL(2) "
      f"und SU(2)_L: S_L (x) 2_L hat die Casimirwerte {dict(spins)}, also Spin 0 + 1 (ganzzahlig, bei Grassmann-ungeraden "
      "Feldern die Verdrillung von Marcus/Vafa-Witten). In der gelösten 16 ist SU(2)_L dagegen schwacher Isospin: "
      f"Y = T3_R + (B-L)/2 mit (B-L)/2 in {[str(x) for x in bl_werte]} vertauscht mit SU(2)_L ({komm_L:.0e}), nicht mit "
      f"SU(2)_R ({komm_R:.1f}), und Q = T3_L + Y erfüllt [Q, T_L^+] = T_L^+ (Defekt {q_rot:.0e}): als Lorentzgenerator "
      "drehte T_L die elektrische Ladung. Die Gleichsetzung widerspricht dem gelösten Block (16 = eine Generation, Spin 1/2)",
      dict(spins) == {0.0: 1, 2.0: 3} and komm_L < 1e-12 and komm_R > 0.1 and q_rot < 1e-12
      and bl_werte == [-0.5, -0.166666667, 0.166666667, 0.5])

say("  (e) Gluonzähler: genau zwei eta^4-Einsetzungen; eta als Hilfsgrößen (GL(4)); N = 4 nötig")
gl_ok, gl_tab = True, {}
for n in (4, 5, 6):
    lr = [(Fr(int(rng.integers(-9, 10))), Fr(int(rng.integers(-9, 10)))) for _ in range(n)]
    while any(ang(lr[i], lr[j]) == 0 for i in range(n) for j in range(i + 1, n)):
        lr = [(Fr(int(rng.integers(-9, 10))), Fr(int(rng.integers(-9, 10)))) for _ in range(n)]
    nz = Counter()
    for lage in itertools.product((0, 1), repeat=n):
        fak = []
        for i, s in enumerate(lage):
            if s:
                for A in range(4):
                    fak.append({1 << (2*A): lr[i][0], 1 << (2*A + 1): lr[i][1]})
        v = oben(gprod(fak), 8)
        if v != 0:
            nz[sum(lage)] += 1
            i, j = [k for k in range(n) if lage[k]]
            gl_ok &= v == ang(lr[i], lr[j])**4
    gl_tab[f"n = {n}"] = {"Lagen": 2**n, "nicht null (Zahl der eta^4)": dict(nz)}
    gl_ok &= dict(nz) == {2: n*(n - 1)//2}
Mi = rng.integers(-3, 4, size=(4, 4))
while round(np.linalg.det(Mi)) == 0:
    Mi = rng.integers(-3, 4, size=(4, 4))
detM = int(round(np.linalg.det(Mi)))
lr = [(Fr(3), Fr(-2)), (Fr(1), Fr(5)), (Fr(-4), Fr(1))]
fak = []
for i in (0, 1):
    for A in range(4):
        fak.append({1 << (2*B): Fr(int(Mi[A, B]))*lr[i][0] for B in range(4)} | {1 << (2*B + 1): Fr(int(Mi[A, B]))*lr[i][1]
                                                                                  for B in range(4)})
gl4 = oben(gprod(fak), 8)
gl4_ok = gl4 == detM**2*ang(lr[0], lr[1])**4
hel_N = {}
for N in (3, 4, 5):
    v1 = ang(lb[0], lb[1])**N
    v2 = ang(1.9*lb[0], lb[1])**N
    hel_N[N] = gs(1 - (np.log(abs(v2/v1))/np.log(1.9))/2, 6)
ZAHLEN["gluonzaehler"] = {"Lagen": gl_tab, "GL(4): Zähler/<12>^4, det^2": [str(gl4/ang(lr[0], lr[1])**4), detM**2],
                          "Helizität des eta^N-Beins": hel_N}
check("[E] (e) Gluonen (je Punkt eta^0 oder eta^4) mit acht Grassmann-Nullmoden theta^(A alpha) lambda_alpha, exakt rational: "
      f"{gl_tab}; nicht null genau bei zwei eta^4-Einsetzungen, dann <ij>^4, alle anderen Helizitätslagen (keine, eine, drei "
      "und mehr negative) verschwinden. Unter eta -> M eta (beliebiges M in GL(4), det M = "
      f"{detM}) wird der Zähler det(M)^2 <ij>^4 ({gl4_ok}): Gluonbäume sehen nur das eps-Singulett eta^1 eta^2 eta^3 eta^4, die "
      "SU(4)-Etiketten sind Hilfsgrößen (N = 4 und reines Yang-Mills haben dieselben Gluonbäume). Die Zahl vier ist nicht "
      f"frei: Helizität des eta^N-Beins {hel_N} (N = 3, 4, 5); nur N = 4 gibt ein Gluon",
      gl_ok and gl4_ok and hel_N[4] == -1 and hel_N[3] == -0.5 and hel_N[5] == -1.5)
befund("Teil 2: Das Pati-Salam-Quartett kann die eta^A nativ nicht stellen. Gestellt wären Zahl und Statistik: vier "
       "Farben, mit dem SU(2)_L-Index genau acht Komponenten (4,2,1), und unter P_CAR R^2-ungerade. Die Zerlegung "
       "(4,2,1) + (4-quer,1,2) ist sogar exakt die der 16 Superladungen von N = 4, falls Spin(4)_PS die Lorentzgruppe ist "
       "(10D-Reduktion). Es fehlt die Homogenität +1: Nativ sind die Farben innere Etiketten (h = 0) oder O(-1)-Felder "
       "(h = -1). Die nötige Gleichsetzung SU(2)_L = SL(2)_lambda ist eine topologische Verdrillung, die die Materie "
       "ganzzahlig macht und die elektrische Ladung dreht; sie widerspricht dem gelösten 16-Block. Da die eta nur als "
       "eps-Singulett in Gluonbäume eingehen, dürfen sie Hilfsgrößen sein; vier ungerade O(1)-Richtungen bleiben "
       "aber nötig, und das ist P_N4 unverändert", True)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

# =================================================================== Teil 3
say("\nTeil 3: Mehrspurterme")
say("  (a) Niveau: Gewicht k^(c-1), M Nahtpaare diagonal = Niveau M, großes Niveau = freier Grenzfall")
diag_ok, diag_tab = True, {}
for n in (4, 5):
    for M in (1, 2, 3):
        for versuch in range(30):                                            # zufällige Entartung meiden
            zs = punkte(n)
            ops2 = [(su2_zufall(), z) for z in zs]
            W = ward(ops2, FORM)
            if set(W) == {1, 2}:
                break
        opsM = [(np.kron(X, np.eye(M, dtype=np.int64)), z) for X, z in ops2]
        K_ = wick(opsM)
        soll = {c: v*M**c for c, v in W.items()}
        diag_ok &= K_ == soll and len(W) >= 2
        diag_tab[f"n = {n}, M = {M}"] = {"k-Potenzen": sorted(W), "Wick(M Kopien) = Ward(k = M)": K_ == soll}
frei, anteile = {}, {}
for n in (4, 6):
    for versuch in range(30):
        zs = punkte(n)
        ops8 = [(gl8_spurfrei(), z) for z in zs]
        W = ward(ops8, FORM)
        if set(W) == set(range(1, n//2 + 1)):
            break
    anteile[n] = []
    for kk in (1, 10, 100, 10000):
        tot = sum(abs(float(v))*kk**c for c, v in W.items())
        anteile[n].append(abs(float(W[n//2]))*kk**(n//2)/tot)
        frei[f"n = {n}, k = {kk}"] = {"Anteil maximale Spurzahl c = n/2": gs(abs(float(W[n//2]))*kk**(n//2)/tot, 4),
                                     "|Zweispur|/|Einspur|": gs(abs(float(W[2]))*kk/abs(float(W[1])), 4)}
frei_ok = all(all(b_ > a_ for a_, b_ in zip(v, v[1:])) and v[-1] > 0.9 for v in anteile.values()) and all(
    abs(frei[f"n = {n}, k = 10"]["|Zweispur|/|Einspur|"]/frei[f"n = {n}, k = 1"]["|Zweispur|/|Einspur|"] - 10) < 1e-3
    for n in (4, 6))
ZAHLEN["niveau"] = {"Diagonale aus M Nahtpaaren": diag_tab, "großes Niveau": frei,
                    "native Niveaus": {"su(2) zweier Nähte (z99)": 1, "su(8) des Originals": 1, "so(16)": 1, "E8": 1}}
check("[E] (a) Der c-Spurterm trägt k^c, relativ zum Einspurterm k^(c-1), unabhängig von jeder Kopplungsnormierung (sie ist für "
      "alle Spurzahlen dieselbe). M Nahtpaare diagonal (je ein su(2)_1 wie z99; X -> X (x) 1_M, Fermionschleifen) ergeben "
      f"exakt die Wardantwort bei k = M (n = 4, 5, M = 1, 2, 3: {diag_ok}): Viele Nähte erhöhen das Niveau und damit das "
      "Mehrspurgewicht. Großes "
      f"Niveau ist der freie, abelsche Grenzfall (Anteil der maximalen Spurzahl wächst monoton gegen 1): {frei}. Nativ sind alle Niveaus 1 (su(2)_1, su(8)_1, so(16)_1, E8_1): Eine "
      "Trennung durch das Niveau gibt es nicht, und wo das Niveau wächst, wächst der Mehrspuranteil",
      diag_ok and frei_ok)

say("  (b) Uhr/GSO: Doppelspur = Vakuumaustausch; Projektoren auf dem Fockraum")


def cdag_matrix(j, nb):
    dim_ = 1 << nb
    s = np.arange(dim_)
    src = s[((s >> j) & 1) == 0]
    tief = src & ((1 << j) - 1)
    par = np.zeros(len(src), dtype=np.int64)
    for b in range(j):
        par ^= (tief >> b) & 1
    return sps.csr_matrix(((1 - 2*par).astype(float), (src | (1 << j), src)), shape=(dim_, dim_))


class Naehte:
    """Nähte zu je 2K Moden x = alpha + n, gemeinsame CAR, Vakuum = Fermisee x < 0 (wie z99/z110)."""

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
        xs = np.array([float(x) for _, x in self.modes])
        unten = xs < 0
        naht = np.array([i for i, _ in self.modes])
        self.Q = [occ[:, (naht == i) & ~unten].sum(1) - (1 - occ[:, (naht == i) & unten]).sum(1) for i in range(self.N)]
        self.L0 = (np.where(unten[None, :], 1 - occ, occ)*np.abs(xs)[None, :]).sum(1)
        vi = int(sum(1 << j for j in range(self.NB) if unten[j]))
        self.vac = np.zeros(self.dim)
        self.vac[vi] = 1.0
        self.cache = {}

    def strom(self, X, n):
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

    def paar(self, n):
        """P^+_n = sum_x psi*_(0,x) psi*_(1,-n-x) (Paarstrom des Originalpaars, z110)."""
        Mx = sps.csr_matrix((self.dim, self.dim))
        for (i, x), j in self.idx.items():
            t = self.idx.get((1, -n - x)) if i == 0 else None
            if t is not None:
                Mx = Mx + self.cd[j] @ self.cd[t]
        return Mx.tocsr()

    def korr(self, liste, P=None, nach=2):
        v = self.vac
        for pos, (X, m) in enumerate(reversed(liste)):
            if P is not None and pos == len(liste) - nach:
                v = P*v
            v = self.strom(X, m) @ v
        return float(self.vac @ v)


S2 = Naehte([Fr(1, 4), Fr(1, 4)], 3)
R_S = np.exp(2j*np.pi*S2.L0)
QS = S2.Q[0] + S2.Q[1]
uhr_S = np.abs(R_S - 1j**(QS % 4)).max()
PROJ = {"R = 1": (np.abs(R_S - 1) < 1e-9).astype(float), "R^2 = (-1)^F = 1": (QS % 2 == 0).astype(float),
        "Q = 0": (QS == 0).astype(float)}
rang_P = {k: f"{int(P.sum())}/{S2.dim}" for k, P in PROJ.items()}
STROEME = (EP, EM, HH)
modul = [S2.vac] + [S2.strom(X, -m) @ S2.vac for X in STROEME for m in (1, 2)]
modul += [S2.strom(X, -1) @ (S2.strom(Y, -1) @ S2.vac) for X in STROEME for Y in STROEME]
modul = [v for v in modul if np.linalg.norm(v) > 1e-12]
dim_modul = int(np.linalg.matrix_rank(np.array(modul)))
leck = {k: gs(max(np.linalg.norm((1 - P)*v)/np.linalg.norm(v) for v in modul), 3) for k, P in PROJ.items()}
LISTEN = [[(EP, 1), (EM, -1), (HH, 1), (HH, -1)], [(EP, 1), (HH, 1), (EM, -1), (HH, -1)],
          [(EM, 1), (EP, 1), (EP, -1), (EM, -1)], [(HH, 2), (EP, -1), (EM, -1)]]
NAME = {EP.tobytes(): "E+", EM.tobytes(): "E-", HH.tobytes(): "H"}
proj_tab, proj_ok = {}, True
for liste in LISTEN:
    voll = S2.korr(liste)
    werte = {k: S2.korr(liste, P, nach=len(liste) - 2) for k, P in PROJ.items()}
    vak = S2.korr(liste, S2.vac, nach=len(liste) - 2)
    proj_ok &= all(abs(v - voll) < 1e-12 for v in werte.values())
    proj_tab[" ".join(f"{NAME[X.tobytes()]}_{m}" for X, m in liste)] = {"voll": gs(voll), **{k: gs(v) for k, v in werte.items()},
                                                                     "nur Vakuum": gs(vak)}
O2 = Naehte([Fr(1, 4), Fr(3, 4)], 3)
R_O = sps.diags(np.exp(2j*np.pi*O2.L0))
R_Oi = sps.diags(np.exp(-2j*np.pi*O2.L0))
h = Fr(1, 2)
anti = max(abs(R_O @ O2.strom(X, m) @ R_Oi + O2.strom(X, m)).max() for X, m in ((EP, h), (EM, -h), (EP, 3*h), (EM, h)))
PR1 = (np.abs(np.exp(2j*np.pi*O2.L0) - 1) < 1e-9).astype(float)
e_zust = O2.strom(EP, -h) @ O2.vac
weg = float(np.linalg.norm(PR1*e_zust)/np.linalg.norm(e_zust))
pp = {m_v: O2.paar(m_v) for m_v in (-1, 1)}
komm_P = max(abs(R_O @ pp[m_v] @ R_Oi - pp[m_v]).max() for m_v in pp)
vak_P = float(O2.vac @ (pp[-1].T @ (pp[-1] @ O2.vac)))
ZAHLEN["projektion"] = {"gleiche Quarter: |e^(2 pi i L_0) - i^(Q_a+Q_b)|": gs(uhr_S, 3), "Rang der Projektoren": rang_P,
                        "Strommodul bis Stufe 2: Dimension": dim_modul, "Leck (1 - P) auf dem Strommodul": leck,
                        "Korrelatoren mit Projektor im Kanal (2, 2)": proj_tab,
                        "Originalpaar: |R E R^-1 + E|": gs(anti, 3), "|P_(R=1) E^+_(-1/2) Omega|/|E^+_(-1/2) Omega|": gs(weg, 3),
                        "|R P R^-1 - P|": gs(komm_P, 3), "<P^-_1 P^+_-1> (Vakuumaustausch des Paarungs-su(2))": gs(vak_P, 6)}
check("[E] (b) Zwei Nähte gleicher Quarter (Fockraum 2 x 6 Moden): Die Uhr ist dort e^(2 pi i L_0) = i^(Q_a + Q_b) "
      f"(Abweichung {uhr_S:.0e}). Die nativen Projektoren auf R = 1, R^2 = (-1)^F = 1 und Q = 0 sind auf dem Fockraum "
      f"echt ({rang_P}), wirken aber auf dem ganzen Strommodul (Vakuum und alle su(2)-Stromzustände bis Stufe 2, Dimension "
      f"{dim_modul} = 1 + 3 + 4 wie der su(2)_1-Vakuumcharakter) als Eins (Leck {leck}): Jeder Stromkorrelator bleibt "
      f"unter jeder solchen Projektion gleich, auch mit dem Projektor in einem Kanal ({proj_ok}: {proj_tab}). Die "
      "Doppelspur ist der Vakuumaustausch (Kanal mit Moden "
      "(m, -m): Korrelator = nur Vakuum; Teil 3e: k^2-Anteil = Produkt der Zweipunktantworten), also der Einsoperator im "
      "zentralen OPE-Term k/(z - w)^2, und der Einsoperator überlebt jede Projektion. Beim Originalpaar antikommutiert R "
      f"mit den getwisteten Transferströmen (|R E R^-1 + E| = {anti:.0e}), die Projektion R = 1 löscht E^+_(-1/2) Omega "
      f"ganz (Rest {weg:.0e}); die Paarströme vertauschen mit R ({komm_P:.0e}, Vakuumaustausch <P^- P^+> = {vak_P:.0f}). "
      "Die Uhr entfernt also Operatoren, nicht Spurstrukturen",
      uhr_S < 1e-9 and proj_ok and all(v < 1e-12 for v in leck.values()) and dim_modul == 8
      and all(int(v.split("/")[0]) < S2.dim for v in rang_P.values()) and anti < 1e-12 and weg < 1e-12
      and komm_P < 1e-12 and abs(vak_P - 1) < 1e-12
      and abs(proj_tab["E+_1 E-_-1 H_1 H_-1"]["nur Vakuum"] - proj_tab["E+_1 E-_-1 H_1 H_-1"]["voll"]) < 1e-12)

say("  (c) Farbstruktur: Rang der Spurstrukturen, U(1) der U(2)")


def gens_su(N):
    G = []
    for i in range(N):
        for j in range(i + 1, N):
            E = np.zeros((N, N), complex)
            E[i, j] = E[j, i] = 0.5
            G.append(E)
            E = np.zeros((N, N), complex)
            E[i, j], E[j, i] = -0.5j, 0.5j
            G.append(E)
    for k in range(1, N):
        D = np.zeros((N, N), complex)
        D[:k, :k] = np.eye(k)
        D[k, k] = -k
        G.append(D/np.sqrt(2*k*(k + 1)))
    return G


PAARUNGEN = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))


def strukturen(Xs):
    out = [np.trace(Xs[0] @ Xs[s[0]] @ Xs[s[1]] @ Xs[s[2]]) for s in itertools.permutations([1, 2, 3])]
    out += [np.trace(Xs[a] @ Xs[b])*np.trace(Xs[c] @ Xs[d]) for (a, b), (c, d) in PAARUNGEN]
    kom = lambda A, B: A @ B - B @ A
    out += [np.trace(kom(kom(Xs[0], Xs[1]), Xs[2]) @ Xs[3]), np.trace(kom(kom(Xs[0], Xs[2]), Xs[1]) @ Xs[3])]
    return out


rang = {}
for name, G in (("su(2)", gens_su(2)), ("u(2)", gens_su(2) + [np.eye(2)/np.sqrt(2)]), ("su(3)", gens_su(3)),
                ("su(4)", gens_su(4)), ("su(8)", gens_su(8))):
    zeilen = np.array([strukturen([sum(rng.normal()*g for g in G) for _ in range(4)]) for _ in range(70)])
    rk = lambda Mx: int(np.sum(np.linalg.svd(Mx, compute_uv=False) > 1e-9*np.linalg.svd(Mx, compute_uv=False)[0]))
    rang[name] = {"9 Spurstrukturen": rk(zeilen[:, :9]), "2 f-Ketten + 3 Doppelspuren": rk(zeilen[:, 6:])}
ID2 = np.eye(2, dtype=np.int64)
u1_tab, u1_ok = {}, True
for n in (4, 5):
    for versuch in range(30):
        zs = punkte(n)
        ein = ward([(ID2, zs[0])] + [(su2_zufall(), z) for z in zs[1:]], FORM)
        zwei = ward([(ID2, zs[0]), (ID2, zs[1])] + [(su2_zufall(), z) for z in zs[2:]], FORM)
        if zwei:
            break
    u1_tab[f"n = {n}"] = {"ein Photon: k-Potenzen": sorted(ein), "zwei Photonen: k-Potenzen": sorted(zwei)}
    u1_ok &= 1 not in ein and 1 not in zwei and 2 in zwei
ZAHLEN["farbstruktur"] = {"Rang (70 Zufallssätze)": rang, "U(1) der U(2)": u1_tab}
check(f"[E] (c) Rang der n = 4-Spurstrukturen über zufälligen Generatoren: {rang}. Für su(2), den spurfreien Teil der "
      "U(2)-Transferfarbe zweier Nähte (die eigentlichen Gluonen), spannen die sechs Einspur- und drei "
      "Doppelspurstrukturen nur drei Dimensionen auf, und die f-Ketten liegen im Spann der Doppelspuren: Ein- und Mehrspur "
      "sind farblich nicht trennbar. Mit dem U(1)-Generator steigt der Rang auf 9, aber nur über Beine des U(1). Erst ab "
      "su(3) sind f-Ketten und Doppelspuren unabhängig (5), ab su(4) alle neun. Das U(1) der U(2) (Gesamtladung, z111) "
      "entkoppelt exakt vom Einspuranteil, nicht "
      f"vom Mehrspuranteil (Wardrekursion, exakt: {u1_tab}): Mit U(2)-Farbe koppelt das 'Photon' nur über die Mehrspur- "
      "(Gravitations-)Terme. Die U(2)-Vollständigkeit klebt Spuren ohne neue Doppelspuren (z111), entfernt aber nicht die "
      "der Vertizes",
      rang["su(2)"] == {"9 Spurstrukturen": 3, "2 f-Ketten + 3 Doppelspuren": 3}
      and rang["su(3)"]["2 f-Ketten + 3 Doppelspuren"] == 5
      and rang["su(4)"]["9 Spurstrukturen"] == 9 and rang["su(8)"]["9 Spurstrukturen"] == 9 and u1_ok)

say("  (d) P_diag/Doppelkopie: DDM-Projektion, BCJ, KLT mit Doppelspur (su(2) gegen su(3))")
L4 = "abcd"


def spur_t(T4, o):
    return np.einsum("".join(L4[i] for i in o) + "->abcd", T4)


def dd_t(tr2, p, q, r1, t):
    return np.einsum(f"{L4[p]}{L4[q]},{L4[r1]}{L4[t]}->abcd", tr2, tr2)


def farbamp(G, lam, k=1.0, mit_dt=True, neg=(0, 1)):
    """Farbbekleidete n = 4-Antwort mit P_N4-Zähler (Lesart T): -k sum_sigma Tr(1 sigma) A_PT + k^2 sum Tr Tr
    <ij>^4/(<ab>^2 <cd>^2); Zerlegung nach f-Ketten (DDM) und Doppelspuren."""
    T = np.array(G)
    T4 = np.einsum("aij,bjk,ckl,dli->abcd", T, T, T, T)
    tr2 = np.einsum("aij,bji->ab", T, T)
    A = lambda i, j: ang(lam[i], lam[j])
    num = A(*neg)**4
    amp = np.zeros(T4.shape, complex)
    for sg in itertools.permutations([1, 2, 3]):
        o = (0,) + sg
        amp -= k*spur_t(T4, o)*num/(A(o[0], o[1])*A(o[1], o[2])*A(o[2], o[3])*A(o[3], o[0]))
    if mit_dt:
        for (p, q), (r1, t) in PAARUNGEN:
            amp += k*k*dd_t(tr2, p, q, r1, t)*num/(A(p, q)**2*A(r1, t)**2)
    C1234 = spur_t(T4, (0, 1, 2, 3)) - spur_t(T4, (1, 0, 2, 3)) - spur_t(T4, (2, 0, 1, 3)) + spur_t(T4, (2, 1, 0, 3))
    C1324 = spur_t(T4, (0, 2, 1, 3)) - spur_t(T4, (2, 0, 1, 3)) - spur_t(T4, (1, 0, 2, 3)) + spur_t(T4, (1, 2, 0, 3))
    dds = [dd_t(tr2, p, q, r1, t) for (p, q), (r1, t) in PAARUNGEN]
    basis = [C1234, C1324, sum(dds)] if len(G) == 3 else [C1234, C1324] + dds
    B = np.array([b.ravel() for b in basis]).T
    rkB = np.linalg.matrix_rank(B)
    co = np.linalg.lstsq(B, amp.ravel(), rcond=None)[0]
    rest = np.linalg.norm(B @ co - amp.ravel())/np.linalg.norm(amp)
    return co, rest, rkB == B.shape[1]


ddm_tab, ddm_ok = {}, True
for name, G in (("su(2)", gens_su(2)), ("su(3)", gens_su(3))):
    for probe in range(3):
        lam, lt = kin_komplex(4)
        A = lambda i, j: ang(lam[i], lam[j])
        s = lambda i, j: ang(lam[i], lam[j])*ang(lt[i], lt[j])
        pt1234 = A(0, 1)**4/(A(0, 1)*A(1, 2)*A(2, 3)*A(3, 0))
        pt1324 = A(0, 1)**4/(A(0, 2)*A(2, 1)*A(1, 3)*A(3, 0))
        for mit in (False, True):
            co, rest, eindeutig = farbamp(G, lam, mit_dt=mit)
            a1, a2 = co[0]/(-1.0), co[1]/(-1.0)
            bcj = abs(s(0, 1)*a1 - s(0, 2)*a2)/abs(s(0, 1)*a1)
            a3 = -a1 - a2
            klt = abs((-s(0, 1)*a1*a3)/(-s(0, 2)*a2*a3) - 1)
            dev = max(abs(a1/pt1234 - 1), abs(a2/pt1324 - 1))
            ddm_tab[f"{name}, Punkt {probe}, {'mit' if mit else 'ohne'} Doppelspur"] = {
                "Rest": gs(rest, 2), "eindeutig": bool(eindeutig), "|A/PT - 1|": gs(dev, 3), "BCJ-Rest": gs(bcj, 3),
                "|M_KLT(1234)/M_KLT(1324) - 1|": gs(klt, 3)}
            ddm_ok &= rest < 1e-10 and eindeutig
            if not mit or name == "su(3)":
                ddm_ok &= dev < 1e-10 and bcj < 1e-10 and klt < 1e-10
            else:
                ddm_ok &= dev > 1e-3 and bcj > 1e-3 and klt > 1e-3
ZAHLEN["doppelkopie"] = ddm_tab
kurz = {k.replace(", ", " "): (v["|A/PT - 1|"], v["BCJ-Rest"]) for k, v in ddm_tab.items() if "Punkt 0" in k}
check("[E, Lesart T] (d) P_diag (z112) wirkt auf Helizitäten je Bein und ist farbblind; die Doppelkopie braucht farbgeordnete "
      "Teilamplituden, also die DDM-Koeffizienten der f-Ketten. Zerlegung der farbbekleideten n = 4-Antwort (Einspur + "
      "Doppelspur, P_N4-Zähler, je drei komplexe Kinematikpunkte; Basis eindeutig): Für su(2) (Transferfarbe zweier Nähte) "
      "werden die DDM-Koeffizienten von der Doppelspur verunreinigt, BCJ (s_12 A(1234) = s_13 A(1324)) ist verletzt und KLT "
      "hängt von der Basis ab; ohne Doppelspur oder mit su(3) sind sie exakt Parke-Taylor, BCJ und KLT-Basisunabhängigkeit "
      f"gelten (|A/PT - 1|, BCJ-Rest am ersten Punkt: {kurz}). Ab su(3) trennt die Farbe die Doppelspur ab (sie hat keine "
      "kubische f f-Struktur, Rang 5 = 2 + 3, also kein Bild in der Doppelkopie), entfernt sie aber nicht aus der Eichantwort",
      ddm_ok)

say("  (e) n = 4: Doppelspur gegen Austausch mit Einstein- und Weyl-Propagator (Maxwell-Spannungstensoren)")
nat_ok = True
for versuch in range(3):
    zs = punkte(4)
    ops8 = [(gl8_spurfrei(), z) for z in zs]
    W = ward(ops8, FORM)
    dt_soll = sum(Fr(spur_i(ops8[p][0] @ ops8[q][0])*spur_i(ops8[r1][0] @ ops8[t][0]))
                  / ((zs[p] - zs[q])**2*(zs[r1] - zs[t])**2) for (p, q), (r1, t) in PAARUNGEN)
    nat_ok &= W.get(2, 0) == dt_soll and dt_soll != 0
def kin_gut(n):
    """Komplexe Kinematik ohne fast kollineare Paare (Konditionierung der Vergleichsrechnung)."""
    while True:
        lam, lt = kin_komplex(n)
        sv = [abs(ang(lam[i], lam[j])*ang(lt[i], lt[j])) for i, j in itertools.combinations(range(n), 2)]
        if min(sv) > 0.05*max(sv):
            return lam, lt


def referenz(l, h, lt):
    """Gut konditionierter Referenzspinor (fern vom eigenen Spinor), zufällig gestört."""
    v = lt if h < 0 else l
    return np.array([-np.conj(v[1]), np.conj(v[0])]) + 0.4*np.linalg.norm(v)*zk()


kanal = {"Graviton": [], "Skalar": []}
ref_def, tt_null = 0.0, 0.0
for probe in range(5):
    lam, lt = kin_gut(4)
    A = lambda i, j: ang(lam[i], lam[j])
    s = lambda i, j: ang(lam[i], lam[j])*ang(lt[i], lt[j])
    for neg in itertools.combinations(range(4), 2):
        hel = [-1 if i in neg else 1 for i in range(4)]
        Fa = [feldstaerke(lam[i], lt[i], hel[i], referenz(lam[i], hel[i], lt[i])) for i in range(4)]
        Fb = [feldstaerke(lam[i], lt[i], hel[i], referenz(lam[i], hel[i], lt[i])) for i in range(4)]
        for (p, q), (r1, t) in PAARUNGEN:
            adt = A(*neg)**4/(A(p, q)**2*A(r1, t)**2)
            if set(neg) == {p, q} or set(neg) == {r1, t}:
                Y = fdot(Fa[p], Fa[q])*fdot(Fa[r1], Fa[t])
                Yb = fdot(Fb[p], Fb[q])*fdot(Fb[r1], Fb[t])
                X = tdot(spannung(Fa[p], Fa[q]), spannung(Fa[r1], Fa[t]))
                tt_null = max(tt_null, abs(X)/abs(Y))
                ref_def = max(ref_def, abs(Y/Yb - 1))
                kanal["Skalar"].append((complex(s(p, q)**2*adt/Y), complex(s(p, q)*adt/Y)))
            else:
                X = tdot(spannung(Fa[p], Fa[q]), spannung(Fa[r1], Fa[t]))
                Xb = tdot(spannung(Fb[p], Fb[q]), spannung(Fb[r1], Fb[t]))
                ref_def = max(ref_def, abs(X/Xb - 1))
                kanal["Graviton"].append((complex(s(p, q)**2*adt/X), complex(s(p, q)*adt/X)))
weyl = {k: (np.mean([v[0] for v in L]), max(abs(v[0] - L[0][0]) for v in L)) for k, L in kanal.items()}
einst = {k: max(abs(v[1] - L[0][1]) for v in L)/abs(L[0][1]) for k, L in kanal.items()}
lam, lt = kin_komplex(4)
pol = {}
for eps_ in (1e-3, 1e-5):
    zstern = -ang(lam[0], lam[2])/ang(lam[1], lam[2])*(1 + eps_)
    l2, t2 = lam.copy(), lt.copy()
    l2[0] = lam[0] + zstern*lam[1]
    t2[1] = lt[1] - zstern*lt[0]
    sv = ang(l2[0], l2[2])*ang(t2[0], t2[2])
    adt = ang(l2[0], l2[1])**4/(ang(l2[0], l2[2])**2*ang(l2[1], l2[3])**2)
    pol[eps_] = (sv*adt, sv**2*adt)
pol_s1 = abs(pol[1e-5][0]/pol[1e-3][0])
pol_s2 = abs(pol[1e-5][1]/pol[1e-3][1] - 1)
ZAHLEN["konforme_gravitation_n4"] = {
    "Ward k^2 = sum Tr Tr/(z_ab^2 z_cd^2) (su(8), drei Sätze)": nat_ok,
    "Gravitonkanal: s^2 A_dt/(T.T), Re, Im, Streuung": [gs(weyl["Graviton"][0].real, 12), gs(weyl["Graviton"][0].imag, 3),
                                                        gs(weyl["Graviton"][1], 3)],
    "Skalarkanal: s^2 A_dt/((F.F)(F.F)), Re, Im, Streuung": [gs(weyl["Skalar"][0].real, 12), gs(weyl["Skalar"][0].imag, 3),
                                                            gs(weyl["Skalar"][1], 3)],
    "Einsteinform s A_dt/(T.T): relative Streuung": {k: gs(v, 3) for k, v in einst.items()},
    "Skalarkanal: |T.T|/|(F.F)(F.F)|": gs(tt_null, 3), "Referenzspinoren: max relative Änderung": gs(ref_def, 3),
    "Pol: |s A(eps=1e-5)/s A(eps=1e-3)|, |s^2 A Verhältnis - 1|": [gs(pol_s1, 5), gs(pol_s2, 3)],
    "Zahl der Kanäle (5 Punkte x 6 Helizitätslagen x 3 Paarungen)": {k: len(v) for k, v in kanal.items()}}
check("[E unter Lesart T und P_N4] (e) Der native Doppelspuranteil bei n = 4 ist k^2 sum Tr(X_a X_b) Tr(X_c X_d)/"
      f"(z_ab^2 z_cd^2) (Wardrekursion, su(8), exakt: {nat_ok}), mit P_N4-Zähler auf der Geraden A_dt = "
      "<ij>^4/(<ab>^2 <cd>^2). Vergleich mit Feldtheorie "
      "(Maxwell-Spannungstensoren T aus linearisierten Feldstärken, referenzunabhängig bis "
      f"{ref_def:.0e}), alle 6 Helizitätslagen x 3 Paarungen an 5 Kinematikpunkten: Liegt in jeder Spur eine negative "
      f"Helizität (Gravitonkanal, {len(kanal['Graviton'])} Fälle), so ist s_ab^2 A_dt = c T_ab.T_cd mit c = "
      f"{weyl['Graviton'][0].real:.9f} (Streuung {weyl['Graviton'][1]:.0e}): Austausch eines Spin-2-Feldes mit "
      f"Weyl-Propagator 1/p^4; die Einsteinform s A_dt/(T.T) streut relativ {einst['Graviton']:.2f}. Liegen beide negativen "
      f"in einer Spur (Skalarkanal, {len(kanal['Skalar'])} Fälle), ist T.T = 0 (relativ {tt_null:.0e}: das Graviton koppelt "
      f"nicht an (-,-)), aber s^2 A_dt = c (F_a.F_b)(F_c.F_d) mit c = {weyl['Skalar'][0].real:.9f} (Streuung "
      f"{weyl['Skalar'][1]:.0e}): ein Skalar mit 1/p^4, gekoppelt an F_-^2 und F_+^2. Pol der Ordnung 2 (s A wächst um "
      f"{pol_s1:.4g} bei eps 1e-3 -> 1e-5, s^2 A konstant bis {pol_s2:.0e}). Das ist der Austausch der N = 4 konformen "
      "Supergravitation (Weyl-Graviton und phi F_-^2 + phi-quer F_+^2, Berkovits-Witten), mit fester Stärke k = 1 relativ "
      "zu Yang-Mills, nicht Einsteingravitation (1/p^2, und im Skalarkanal null)",
      nat_ok and abs(weyl["Graviton"][0] - weyl["Skalar"][0]) < 1e-6 and weyl["Graviton"][1] < 1e-6
      and weyl["Skalar"][1] < 1e-6 and einst["Graviton"] > 0.1 and tt_null < 1e-12 and ref_def < 1e-6
      and abs(pol_s1 - 100) < 1 and pol_s2 < 1e-3)
befund("Teil 3: Keine native Struktur entfernt oder trennt die Mehrspurterme. Das Niveau ist 1, und größeres Niveau (etwa M "
       "Nähte diagonal) verstärkt sie. Die Doppelspur ist Vakuumaustausch; Uhr, (-1)^F und Ladungsprojektion halten ihn fest, "
       "die Uhr entfernt nur die getwisteten Transferströme selbst. Mit dem su(2) der U(2)-Transferfarbe sind Ein- und Mehrspur "
       "farblich untrennbar (Rang 3), das U(1) koppelt nur über Mehrspur, und P_diag bekäme verunreinigte Teilamplituden (BCJ "
       "verletzt). Ab su(3) trennt die Farbe sauber, die Doppelspur bleibt aber als eigener Sektor. Bei n = 4 ist dieser "
       "Sektor exakt der Austausch der konformen Supergravitation (Weyl-Graviton und Skalar, Doppelpol, dimensionslose "
       "Kopplung k = 1), nicht Einstein", True)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

# =================================================================== Teil 4
say("\nTeil 4: Bilanz und Urteil")
STATUS = {
    "Gradzwist O(-1) -> O(1)": ("E (Mathematik) / nein (nativ)",
                                "Kopplung der U(1)-Ladung an K^-1 = O(2) = topologische Verdrillung T + dJ: Index 2, Nullmoden "
                                "theta^alpha lambda_alpha, c = 1 -> -11 je komplexer Fermion; nativ nicht geliefert"),
    "Uhr R": ("E (negativ)", "sieht nur L_0 mod 1: blind für ganzzahlige Zwiste, halber Zwist -> R^-1; liefert und wählt keinen"),
    "Quarterrahmen, Spin-Z4": ("E (negativ)", "flache Z/4-Holonomien, Grad 0, Index 0; mit Gradzwist Zähler <ij>^4 z_i z_j, "
                                              "nur diagonal kovariant (Charakter t^-4): unverträglich"),
    "Tomita J": ("E (negativ)", "Teilchen-Loch-Anteil: L_0 - q Q -> L_0 + q Q, kehrt den Zwist um; J-kovariant nur mit "
                                "Orientierungsbit"),
    "Ladungsanomalie": ("E", "verdrillte Fermionen machen neutrale Stromkorrelatoren null: Nenner und Zähler brauchen "
                             "verschiedene Fermionen"),
    "Pati-Salam-Quartett als eta^A": ("E (negativ, gescopt)",
                                      "Zahl (4, mit SU(2)_L acht Komponenten = Form der 16 Superladungen) und R^2-Statistik "
                                      "passen; Homogenität +1 fehlt (innere Etiketten h = 0, Nahtfelder h = -1); "
                                      "SU(2)_L = SL(2)_lambda wäre eine topologische Verdrillung der Materie gegen die gelöste 16"),
    "eta als Hilfsgrößen": ("E", "Gluonzähler nur über das eps-Singulett: genau zwei eta^4 -> <ij>^4, sonst 0; GL(4) gibt "
                                 "det^2; vier ungerade O(1)-Richtungen bleiben nötig (P_N4 unverändert)"),
    "Mehrspur: Niveau": ("E (negativ)", "Gewicht k^(c-1), natives k = 1; M Nahtpaare diagonal = Niveau M verstärkt"),
    "Mehrspur: Uhr/GSO": ("E (negativ)", "Doppelspur = Vakuumaustausch, von allen nativen Projektoren gehalten; die Uhr "
                                         "entfernt die getwisteten Transferströme selbst"),
    "Mehrspur: U(2)-Farbe": ("E (negativ)", "Spurstrukturen des su(2)-Teils Rang 3: Ein- und Mehrspur untrennbar; U(1) "
                                            "koppelt nur über Mehrspur"),
    "Mehrspur: P_diag": ("E (negativ)", "farbblind; su(2): DDM-Koeffizienten verunreinigt, BCJ verletzt, KLT basisabhängig; "
                                        "ab su(3) sauber getrennt, nicht entfernt"),
    "Doppelspur n = 4": ("E unter Lesart T und P_N4",
                         "= Austausch der N = 4 konformen Supergravitation: Weyl-Graviton (4 T.T/s^2) und Skalar an F_-^2, "
                         "F_+^2 (4 (F.F)(F.F)/s^2), Doppelpol, Stärke k = 1; nicht Einstein"),
    "K15": ("L (präzisiert)", "Nenner und Farbe nativ; P_N4 weder durch Gradzwist noch durch Pati-Salam geliefert; selbst "
                              "mit P_N4 bleibt konforme Supergravitation gleicher Stärke untrennbar beigemischt"),
    "K13": ("L (präzisiert)", "der native Mehrspursektor ist unter T + P_N4 konforme (Weyl-)Gravitation ohne Länge "
                              "(dimensionslos, 1/p^4); Einstein nur über P_diag (z112) oder eine Einstein-Einschränkung"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
urteil = (
    "Alle drei Wege, den Zähler <ij>^4 nativ zu machen oder die Mehrspurterme loszuwerden, sind gescopt verneint; der "
    "Mehrspursektor ist bei n = 4 exakt konforme Supergravitation. Teil 1 (Gradzwist): Mathematisch ist der Gradzwist der "
    "richtige Mechanismus. Die Kopplung der Nahtladung an K^-1 = O(2), also das Inverse des Quadrats der tautologischen "
    "Linie, ist die topologische Verdrillung T' = T + dJ. Sie macht aus den Nahtfermionen O(1)-Felder mit Index 2 "
    "(Cech: O(1) hat h^0 = 2, O(-3) h^1 = 2, O(-1) ist azyklisch) und liefert genau die Nullmoden theta^alpha lambda_alpha. "
    "Je komplexer Fermion fällt c von 1 auf -11; a = 1/2 mit Grad-1-Bündel ist die Twistorstring-Variante (c = -2). "
    "Nativ liefert ihn nichts. Ohne Zwist ist der Index in jedem Rahmen 0. Uhr, Quarterrahmen und Spin-Z4 sind flache "
    "Z/4-Holonomien vom Grad 0, und R = e^(2 pi i L_0) sieht nur L_0 mod 1: blind für ganzzahlige Zwiste, der halbe Zwist "
    "ergibt R^-1. Der Teilchen-Loch-Anteil von J kehrt den Zwist um (L_0 - q Q -> L_0 + q Q), und die nativen "
    "Verbindungen sind flach (z100). Der Quarterrahmen ist mit dem Zwist sogar unverträglich: Die Nullmoden tragen "
    "z^(1/2 - alpha), der Zähler wird <ij>^4 z_i z_j und ist nur unter der Diagonalgruppe mit Charakter t^-4 kovariant. "
    "Nach der Ladungsanomalie machen verdrillte Fermionen die neutralen Stromkorrelatoren zu null; Nenner und Zähler "
    "brauchen also verschiedene Fermionen. Teil 2 (Pati-Salam): Die 16 zerfällt in (4,2,1) + (4-quer,1,2), das ist "
    "exakt die Form der 16 Superladungen Q^(alpha A) + Q-quer, sobald Spin(4)_PS die Lorentzgruppe ist (10D-Reduktion), "
    "und unter P_CAR ist alles R^2-ungerade. Das Supermultiplett Lambda(eta) = 1 + 4 + 6 + 4-quer + 1 steckt nur unter "
    "SU(4) c SU(5) in der 16. Es fehlt die Homogenität +1: Als innere Etiketten geben die Farben einen lambda-freien "
    "Zähler (Helizität +1) oder null, als O(-1)-Nahtfelder keine Nullmoden. Die nötige Gleichsetzung SU(2)_L = "
    "SL(2)_lambda ist eine topologische Verdrillung, die die Materie zu Spin 0 + 1 macht und die elektrische Ladung dreht. "
    "Das widerspricht der gelösten 16. Die eta dürfen Hilfsgrößen sein, denn Gluonbäume sehen nur das eps-Singulett "
    "(genau zwei eta^4 -> <ij>^4, sonst 0, GL(4) -> det^2). Vier Grassmann-ungerade O(1)-Richtungen bleiben aber nötig: "
    "P_N4 ist unverändert. Teil 3 (Mehrspur): Das Mehrspurgewicht ist k^(c-1) bei nativem k = 1. M Nahtpaare diagonal "
    "ergeben Niveau M und verstärken es; großes Niveau ist der freie Grenzfall. Die Doppelspur ist der Austausch des "
    "Einsoperators im zentralen OPE-Term, und Uhr, (-1)^F und Ladungsprojektion halten ihn fest. Die Uhr entfernt beim "
    "Originalpaar die getwisteten Transferströme selbst, nicht deren Mehrspuranteile. Mit dem su(2) der U(2)-Transferfarbe "
    "sind Ein- und Mehrspur untrennbar (Rang 3 der Spurstrukturen), und das U(1) koppelt nur über Mehrspur. P_diag ist "
    "farbblind und bekäme für su(2) verunreinigte Teilamplituden (BCJ verletzt, KLT basisabhängig). Ab su(3) trennt "
    "die Farbe sauber, ohne die Doppelspur aus der Eichantwort zu entfernen. Bei n = 4 ist die native Doppelspur mit "
    "P_N4-Zähler feldtheoretisch genau ein Austausch mit 1/p^4: Mit einer negativen Helizität je Spur ist s^2 A = 4 "
    "T.T der Maxwell-Spannungstensoren (Weyl-Graviton). Mit beiden negativen in einer Spur ist T.T = 0 und s^2 A = 4 "
    "(F.F)(F.F), ein Skalar an F_-^2 und F_+^2. Der Pol hat Ordnung 2. Das ist die N = 4 konforme Supergravitation "
    "nach Berkovits-Witten, mit dimensionsloser Kopplung fester Stärke k = 1 relativ zu Yang-Mills, nicht Einstein. "
    "Folgen: K15 bleibt L, präzisiert. Nenner und Farbstruktur sind nativ; P_N4 ist weder durch einen nativen Gradzwist "
    "noch durch das Pati-Salam-Quartett gedeckt und reduziert sich auf vier zusätzliche, ladungsverdrillte "
    "Fermionrichtungen. Selbst mit P_N4 bleibt konforme Supergravitation gleicher Stärke beigemischt, reines "
    "Yang-Mills ist nativ nicht erreichbar. K13 bleibt L, präzisiert: Der native Gravitationsanteil der Lesart T ist "
    "konforme (Weyl-)Gravitation, vierte Ordnung, ohne Länge (passend zu z112: keine native Länge für G). "
    "Einsteingravitation gibt es nur über die blattdiagonale Doppelkopie (P_diag, z112), deren Eingang für "
    "su(2)/U(2)-Farbe selbst verunreinigt ist, oder über eine Einstein-Einschränkung der konformen Gravitation "
    "(Maldacena 2011, Adamo-Mason) mit einer Unendlichkeitsstruktur. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
say(f"\n   (Laufzeit {time.time() - T0:.1f} s)")
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z124_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z124.txt").write_text("\n".join(OUT) + "\n", encoding="utf-8")
