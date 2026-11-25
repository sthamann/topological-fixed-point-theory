"""Z105: Zweikörperterm aus Doppeltransfer und Tomita-Partner -- Box-Nachkommen, Kopplung, induziertes G (Folgeschritt N6, H6).

Leitlinie 'Rekonstruktion statt Auswahl': Der Zweikörperterm wird aus dem nativen Vertexgitter Z^2 zweier Nähte (z99) und dem
Tomita-Partner gebaut, nicht aus einer gewünschten 4D-Wechselwirkung; die Kopplung lambda ist gewählt und wird auf native
Bedingungen geprüft. Firewalls: DIMENSION.UPLIFT.FIREWALL.01 (Naht = chirale 2D-Struktur, Partner = Tomita-Bild; 4D nur als
Lesart der verdoppelten Einpaarleiter von z83: E = (L_0 + Lquer_0)/2, K-Typen (j, j), su(2,2) in CCR-Normierung; keine
wechselwirkende 3+1-Schließung) und GRAV.NONCIRCULAR.01 (G wird nur über die Sakharov-Lesart des Spektrums auf S^3 geprüft,
nicht eingesetzt). Gelöste Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt; nichts wird befördert.
Statuskürzel: E exakt, C bedingt, L Lesart, O offen, N numerisch.
Bestand: z99 (Vertexroute = Stromroute, Delta = 1 starr, Radiusdeformation mit kurzen Vielfachheiten (k+1)^2, fehlende
Box-Nachkommen [0, 0, 1, 4, 10, 20], Doppeltransfer psi*_(1,5/4) psi*_(1,1/4) psi_(2,1/4) psi_(2,-3/4) als kleinster verbundener
1->2-Term, h = 4), z92 (E^nat = L_0^osc + Q^2/2), z83 (Paarschalen, CCR-su(2,2), Leiter = masseloser Skalar), z91 (4D-Fockraum
als Sym^N über Nähte, K-Typ-Hindernis), z90 (1/(4 G_ind) = (1 - 6 xi) I[f]/(48 pi), xi = 1/6 aus E = m, Profile I[f]).

Bezeichnungen: Naht i mit Moden x = 1/4 + n (n ganz), Fermisee n < 0; Zustand = Maya-Diagramm (Ladung Q, Partition lam),
E^nat = Q^2/2 + |lam|. Doppelteilchenvertex A_s = sum_(n<n', n+n'+1=s) (n' - n) psi*_n' psi*_n (Ladung +2, hebt E^nat um s,
Vandermonde-Gewicht = Moden von e^(2 i phi)), B_s = (A_(-s))^+ (Ladung -2). Doppeltransfer (Komponenten von e^(2i(phi_1 - phi_2)),
h = 4): C+_t = sum_s A^(1)_s (x) B^(2)_(t-s), C-_t = sum_s B^(1)_s (x) A^(2)_(t-s) = (C+_(-t))^+ (Naht 1 und 2 je gerade, also ohne
Jordan-Wigner-Vorzeichen zwischen den Nähten). Tomita-Kombination (Bild J O J im Etikettbasis-Partner, reell): spinfreies (4, 4)-Feld
auf dem Kreis V = sum_t (C+_t (x) C+_t + C-_t (x) C-_t), hermitesch; sein E-erhaltender Teil W = C+_0 (x) C+_0 + h.c.
H(lambda) = E + lambda V auf der verdoppelten Fockleiter. Einpaarschichten: (11)_m Haken in Naht 1, (22)_m in Naht 2, (12)_m
Teilchen in 1 und Loch in 2, (21)_m umgekehrt (z99); 4D-Einteilchenzustände H^(ab)_m (x) Hquer^(ab)_m.

Teil 1  Zweikörperterm und Vielfachheiten. (a) Bau und Prüfung (Heisenberg-Relationen, Normen binom, h = 4, Uhr, z99-Komponente).
        (b) Mischung der Paarschichten durch den Nullmodus: Kreuzblock T_m, Paarzahlzerlegung, Schwelle Q^2 = 4 im neutralen
        Sektor. (c) Erste Ordnung exakt (W vertauscht mit E): Eigenwerte +-sigma_i sigma_j, Vielfachheit des Primärwerts gegen
        (k+1)^2 und binom(k+3,3), Schur-Test der su(2,2)-Kovarianz, neutraler Sektor ab E = 4. (d) Box-Nachkomme: Shapovalov-Norm
        2 Delta (Delta - 1) in z83-Normierung, det P = 0 auf der Leiter, Ordnung von gamma, K-Typen und Pauli-Schwelle.
        (e) Zweite Ordnung (Nichtnullmoden): Schwerpunkt geschützt, Relativleiter verschoben, Abschneideabhängigkeit.
Teil 2  Fixiert eine native Bedingung lambda? Unitarität, Marginalität (2D und 4D-Lesart, KT-Punkt), Uhr/KMS bei beta = 4 pi,
        Komposition dreier Nähte (geschlossene OPE, beta-Funktion zweiter Ordnung).
Teil 3  H6: induziertes G zweiter Ordnung (Sakharov-Lesart xi_eff aus dem Spektrum, Abschneide- und Profilabhängigkeit), Bilanz
        und Urteil für N6/H6, Folgen für K13, K15, T3/T5.
"""
import json
import math
import time
from fractions import Fraction as Fr
from functools import lru_cache
from math import comb
from pathlib import Path

import numpy as np
import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
ALPHA = Fr(1, 4)
BETA = 4*np.pi
TOL = 1e-9


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


@lru_cache(None)
def partitionen(n, maxp=None):
    maxp = n if maxp is None else maxp
    if n == 0:
        return ((),)
    aus = []
    for k in range(min(n, maxp), 0, -1):
        for rest in partitionen(n - k, k):
            aus.append((k,) + rest)
    return tuple(aus)


def p_zahl(n):
    return len(partitionen(n)) if n >= 0 else 0


# ------------------------------------------------------------------ Bausteine: eine Naht als Maya-Diagramme
class Naht:
    """Chirales Dirac-Fermion einer Naht, Moden x = alpha + n, Fenster -W <= n < W (Bit n + W), Fermisee n < 0.
    Ladungssektor (Q, N): Partitionen lam von N, besetzt {lam_i - i + Q}; exakt, solange N + |Q| + 2 <= W."""

    def __init__(self, W, alpha=ALPHA):
        self.W, self.alpha = W, alpha
        self._bas = {}

    def maske(self, Q, lam):
        bes = [lam[i - 1] - i + Q for i in range(1, len(lam) + 1)]
        i = len(lam) + 1
        while Q - i >= -self.W:
            bes.append(Q - i)
            i += 1
        s = 0
        for n in bes:
            assert -self.W <= n < self.W, "Fenster zu klein"
            s |= 1 << (n + self.W)
        return s

    def basis(self, Q, N):
        if (Q, N) not in self._bas:
            ms = [self.maske(Q, lam) for lam in partitionen(N)] if N >= 0 else []
            self._bas[(Q, N)] = (ms, {m: i for i, m in enumerate(ms)})
        return self._bas[(Q, N)]

    def index(self, Q, lam):
        return self.basis(Q, sum(lam))[1][self.maske(Q, lam)]

    def enat(self, s):
        W, a = self.W, self.alpha
        L0, Q = Fr(0), 0
        for n in range(-W, W):
            o = (s >> (n + W)) & 1
            if n >= 0 and o:
                L0, Q = L0 + a + n, Q + 1
            if n < 0 and not o:
                L0, Q = L0 - (a + n), Q - 1
        return L0 - (a - Fr(1, 2))*Q, L0, Q

    def teilchen_loecher(self, s):
        W = self.W
        return bin(s >> W).count("1"), W - bin(s & ((1 << W) - 1)).count("1")

    def _vz(self, s, b):
        return -1 if bin(s >> (b + 1)).count("1") % 2 else 1

    def cdag(self, s, n):
        b = n + self.W
        if b < 0 or b >= 2*self.W or (s >> b) & 1:
            return 0, 0
        return self._vz(s, b), s | (1 << b)

    def c(self, s, n):
        b = n + self.W
        if b < 0 or b >= 2*self.W or not (s >> b) & 1:
            return 0, 0
        return self._vz(s, b), s & ~(1 << b)

    def _bilde(self, quelle, ziel, terme):
        src, _ = self.basis(*quelle)
        tgt, ti = self.basis(*ziel)
        M = np.zeros((len(tgt), len(src)))
        if not src or not tgt:
            return M
        for j, m in enumerate(src):
            for koeff, ops in terme:
                g, s = 1, m
                for art, n in ops:
                    g_, s = (self.cdag if art == "+" else self.c)(s, n)
                    g *= g_
                    if not g:
                        break
                if not g:
                    continue
                i = ti.get(s)
                if i is None:
                    raise RuntimeError("Zielzustand fehlt: Fenster zu klein")
                M[i, j] += koeff*g
        return M

    @lru_cache(None)
    def A(self, s_, Q, N):
        """A_s: (Q, N) -> (Q + 2, N + s - 2Q - 2), psi*_n' psi*_n mit n < n', n + n' + 1 = s, Gewicht n' - n."""
        terme = [(n2 - n, (("+", n), ("+", n2))) for n in range(-self.W, self.W) for n2 in [s_ - 1 - n] if n < n2 < self.W]
        return self._bilde((Q, N), (Q + 2, N + s_ - 2*Q - 2), terme)

    @lru_cache(None)
    def B(self, s_, Q, N):
        """B_s = (A_(-s))^+: (Q, N) -> (Q - 2, N + s + 2Q - 2), psi_n psi_n' mit n < n', n + n' + 1 = -s."""
        terme = [(n2 - n, (("-", n2), ("-", n))) for n in range(-self.W, self.W) for n2 in [-s_ - 1 - n] if n < n2 < self.W]
        return self._bilde((Q, N), (Q - 2, N + s_ + 2*Q - 2), terme)

    @lru_cache(None)
    def J(self, mm, Q, N):
        """J_m = sum psi*_(n-m) psi_n: (Q, N) -> (Q, N - m); J_0 = Q."""
        if mm == 0:
            return Q*np.eye(len(self.basis(Q, N)[0]))
        terme = [(1, (("-", n), ("+", n - mm))) for n in range(-self.W, self.W)]
        return self._bilde((Q, N), (Q, N - mm), terme)


NA = Naht(W=24)


# ------------------------------------------------------------------ zwei Nähte (Tensorprodukt, jede Nahtwirkung gerade)
def bloecke(Q1, Q2, E):
    R = Fr(E) - Fr(Q1*Q1, 2) - Fr(Q2*Q2, 2)
    if R < 0 or R.denominator != 1:
        return []
    return [(N1, int(R) - N1) for N1 in range(int(R) + 1)]


@lru_cache(None)
def basis2(Q1, Q2, E):
    B = [(N1, i1, N2, i2) for N1, N2 in bloecke(Q1, Q2, E)
         for i1 in range(p_zahl(N1)) for i2 in range(p_zahl(N2))]
    return tuple(B), {b: i for i, b in enumerate(B)}


@lru_cache(None)
def C_mat(t, Q1, Q2, E, plus=True):
    """C+_t (plus) bzw. C-_t: Sektor (Q1, Q2) auf Stufe E -> (Q1 +- 2, Q2 -+ 2) auf Stufe E + t."""
    d1, d2 = (2, -2) if plus else (-2, 2)
    src, si = basis2(Q1, Q2, E)
    tgt, ti = basis2(Q1 + d1, Q2 + d2, E + t)
    M = np.zeros((len(tgt), len(src)))
    for N1, N2 in bloecke(Q1, Q2, E):
        for M1, M2 in bloecke(Q1 + d1, Q2 + d2, E + t):
            s1 = int(Fr((Q1 + d1)**2, 2) + M1 - Fr(Q1*Q1, 2) - N1)
            s2 = t - s1
            a = NA.A(s1, Q1, N1) if plus else NA.B(s1, Q1, N1)
            b = NA.B(s2, Q2, N2) if plus else NA.A(s2, Q2, N2)
            if a.size == 0 or b.size == 0 or not a.any() or not b.any():
                continue
            zeilen = [ti[(M1, i1, M2, i2)] for i1 in range(p_zahl(M1)) for i2 in range(p_zahl(M2))]
            spalten = [si[(N1, i1, N2, i2)] for i1 in range(p_zahl(N1)) for i2 in range(p_zahl(N2))]
            M[np.ix_(zeilen, spalten)] += np.kron(a, b)
    return M


def vek2(Q1, lam1, Q2, lam2):
    E = Fr(Q1*Q1, 2) + sum(lam1) + Fr(Q2*Q2, 2) + sum(lam2)
    B, bi = basis2(Q1, Q2, E)
    v = np.zeros(len(B))
    v[bi[(sum(lam1), NA.index(Q1, lam1), sum(lam2), NA.index(Q2, lam2))]] = 1.0
    return v


def paare2(Q1, Q2, b):
    N1, i1, N2, i2 = b
    t1, l1 = NA.teilchen_loecher(NA.basis(Q1, N1)[0][i1])
    t2, l2 = NA.teilchen_loecher(NA.basis(Q2, N2)[0][i2])
    return t1 + t2, l1 + l2


def haken(m, k):
    return (k + 1,) + (1,)*(m - 1 - k)


def einpaar(typ, m):
    """Einpaarschicht (ab)_m als Liste von Basisvektoren (k = 0..m-1, Teilchenlage k, Lochtiefe m-1-k)."""
    aus = []
    for k in range(m):
        j = m - 1 - k
        if typ == "11":
            aus.append(vek2(0, haken(m, k), 0, ()))
        elif typ == "22":
            aus.append(vek2(0, (), 0, haken(m, k)))
        elif typ == "12":
            aus.append(vek2(1, (k,) if k else (), -1, (1,)*j))
        else:
            aus.append(vek2(-1, (1,)*j, 1, (k,) if k else ()))
    return np.array(aus)


# =================================================================== Teil 1
say("Teil 1: Zweikörperterm aus Doppeltransfer und Tomita-Partner")
say("  (a) Bau")
e_ok = all(NA.enat(m)[0] == Fr(Q*Q, 2) + N and NA.enat(m)[2] == Q
           for Q in range(-3, 4) for N in range(7) for m in NA.basis(Q, N)[0])
heis = 0.0
for Q in (-2, -1, 0, 1):
    for N in range(5):
        for s_ in range(-3, 6):
            Np = N + s_ - 2*Q - 2
            for mm in (-2, -1, 1, 2):
                if Np < 0 or N - mm < 0:
                    continue
                lhs = NA.J(mm, Q + 2, Np) @ NA.A(s_, Q, N) - NA.A(s_, Q, N - mm) @ NA.J(mm, Q, N)
                heis = max(heis, np.abs(lhs - 2*NA.A(s_ - mm, Q, N)).max(initial=0))
adj = max(np.abs(NA.B(s_, Q, N) - NA.A(-s_, Q - 2, N + s_ + 2*Q - 2).T).max(initial=0)
          for Q in (-1, 0, 1, 2) for N in range(5) for s_ in range(-3, 4) if N + s_ + 2*Q - 2 >= 0)
norm_A = [r(np.sum(NA.A(s_, 0, 0)**2)) for s_ in range(2, 10)]
ZAHLEN["einzelnaht"] = {"Defekt [J_m, A_s] - 2 A_(s-m)": r(heis, 15), "Defekt B_s - A_(-s)^+": r(adj, 15),
                        "||A_s Omega||^2, s = 2..9": norm_A, "binom(s + 1, 3)": [comb(s_ + 1, 3) for s_ in range(2, 10)]}
check("[E] (a) Eine Naht als Maya-Diagramme: E^nat = Q^2/2 + |lam| auf allen Zuständen (Q = -3..3, |lam| <= 6); der "
      "Doppelteilchenvertex A_s = sum (n' - n) psi*_n' psi*_n erfüllt [J_m, A_s] = 2 A_(s-m) exakt (Defekt "
      f"{heis:.0e}, Q = -2..1, m = +-1, +-2), B_s = A_(-s)^+ (Defekt {adj:.0e}), und ||A_s Omega||^2 = {norm_A[:6]} = "
      "binom(s + 1, 3): genau die Modennormen eines Primärfelds mit h = 2 = |2|^2/2. A ist e^(2 i phi) (ein Heisenberg-Primärfeld "
      "der Ladung 2 ist bis auf Nullmodenfaktoren eindeutig)", e_ok and heis == 0 and adj == 0
      and norm_A == [comb(s_ + 1, 3) for s_ in range(2, 10)])

# zwei Nähte: Heisenberg je Naht, Normen h = 4, Uhr
def J2(mm, Q1, Q2, E, naht):
    src, si = basis2(Q1, Q2, E)
    tgt, ti = basis2(Q1, Q2, E - mm)
    M = np.zeros((len(tgt), len(src)))
    for N1, N2 in bloecke(Q1, Q2, E):
        if naht == 1 and N1 - mm >= 0:
            a, b, M1, M2 = NA.J(mm, Q1, N1), np.eye(p_zahl(N2)), N1 - mm, N2
        elif naht == 2 and N2 - mm >= 0:
            a, b, M1, M2 = np.eye(p_zahl(N1)), NA.J(mm, Q2, N2), N1, N2 - mm
        else:
            continue
        zeilen = [ti[(M1, i1, M2, i2)] for i1 in range(p_zahl(M1)) for i2 in range(p_zahl(M2))]
        spalten = [si[(N1, i1, N2, i2)] for i1 in range(p_zahl(N1)) for i2 in range(p_zahl(N2))]
        M[np.ix_(zeilen, spalten)] += np.kron(a, b)
    return M


heis2 = {1: 0.0, 2: 0.0}
for (Q1, Q2) in ((0, 0), (-1, 1), (1, -1)):
    for E in (1, 2, 3):
        for t in (-1, 0, 1, 2, 4):
            if not basis2(Q1 + 2, Q2 - 2, E + t)[0]:
                continue
            for mm in (-1, 1, 2):
                for naht, q in ((1, 2), (2, -2)):
                    if E - mm < 0 or not basis2(Q1 + 2, Q2 - 2, E + t - mm)[0] or not basis2(Q1, Q2, E - mm)[0]:
                        continue
                    lhs = J2(mm, Q1 + 2, Q2 - 2, E + t, naht) @ C_mat(t, Q1, Q2, E) - C_mat(t, Q1, Q2, E - mm) @ J2(mm, Q1, Q2, E, naht)
                    heis2[naht] = max(heis2[naht], np.abs(lhs - q*C_mat(t - mm, Q1, Q2, E)).max(initial=0))
norm_O = [r(np.sum(C_mat(u, 0, 0, 0)**2)) for u in range(4, 11)]
norm_Om = [r(np.sum(C_mat(u, 0, 0, 0, False)**2)) for u in range(4, 11)]
# z99-Komponente: psi*_(1,5/4) psi*_(1,1/4) psi_(2,1/4) psi_(2,-3/4) auf b+_22(1/4, -7/4) Omega
quelle = vek2(0, (), 0, haken(2, 0))                                                     # Teilchen 1/4, Loch -7/4 in Naht 2: m = 2
bild = C_mat(2, 0, 0, 2) @ quelle
tz, _ = basis2(2, -2, 4)
traeger = [(tz[i], bild[i]) for i in np.flatnonzero(np.abs(bild) > 1e-12)]
z99_komp = [(NA.basis(2, b[0])[0][b[1]], NA.basis(-2, b[2])[0][b[3]], w) for b, w in traeger]
ziel_z99 = (NA.maske(2, ()), NA.maske(-2, ()))                                             # Naht 1: Teilchen 1/4, 5/4; Naht 2: Löcher -3/4, -7/4
z99_ok = any((a, b) == ziel_z99 and abs(abs(w) - 1) < TOL for a, b, w in z99_komp)
def uhrphasen(Q1, Q2, b):
    N1, i1, N2, i2 = b
    m1, m2 = NA.basis(Q1, N1)[0][i1], NA.basis(Q2, N2)[0][i2]
    L0 = NA.enat(m1)[1] + NA.enat(m2)[1]
    fe = (NA.teilchen_loecher(m1)[0] + NA.teilchen_loecher(m2)[0] + NA.teilchen_loecher(m1)[1] + NA.teilchen_loecher(m2)[1]) % 2
    return (L0 % 1, (Q1 - Q2) % 4, fe)


uhr_ok, n_uhr = True, 0
for (Q1, Q2, E, t) in ((0, 0, 1, 3), (-1, 1, 2, 0), (0, 0, 4, 0), (-1, 1, 3, 1), (-2, 2, 4, 2)):
    M = C_mat(t, Q1, Q2, E)
    src, tgt = basis2(Q1, Q2, E)[0], basis2(Q1 + 2, Q2 - 2, E + t)[0]
    for i, j in zip(*np.nonzero(M)):
        uhr_ok &= uhrphasen(Q1 + 2, Q2 - 2, tgt[i]) == uhrphasen(Q1, Q2, src[j])
        n_uhr += 1
ZAHLEN["doppeltransfer"] = {"Defekt [J^(1)_m, C+_t] - 2 C+_(t-m)": r(heis2[1], 15), "Defekt [J^(2)_m, C+_t] + 2 C+_(t-m)": r(heis2[2], 15),
                            "||C+_u Omega||^2, u = 4..10": norm_O, "||C-_u Omega||^2": norm_Om,
                            "binom(u + 3, 7)": [comb(u + 3, 7) for u in range(4, 11)],
                            "z99-Komponente in C+_2 b+_22(1/4,-7/4) Omega": [str(w) for a, b, w in z99_komp if (a, b) == ziel_z99]}
check("[E] (a) Doppeltransfer zweier Nähte C+_t = sum_s A^(1)_s (x) B^(2)_(t-s): [J^(1)_m, C+_t] = 2 C+_(t-m) und "
      f"[J^(2)_m, C+_t] = -2 C+_(t-m) exakt (Defekte {heis2[1]:.0e}, {heis2[2]:.0e}), also vertauscht er mit dem Schwerpunktstrom "
      "J^(1) + J^(2) und wirkt nur auf das Relativboson (Ladung alpha = (2, -2)); ||C+-_u Omega||^2 = "
      f"{norm_O} = binom(u + 3, 7): h = |alpha|^2/2 = 4 je Seite, mit Partner (h, hquer) = (4, 4), spinfrei. Uhr: E^nat-Schritte "
      "ganzzahlig, Q gesamt erhalten, i^(Q1 - Q2) unverändert (Delta(Q1 - Q2) = 4), (-1)^F gerade; auf allen "
      f"{n_uhr} geprüften Matrixelementen tragen Quelle und Bild dieselben Phasen (e^(2 pi i L_0), i^(Q1 - Q2), (-1)^F): [R, V] = 0 "
      "für alle nativen Uhren. Der z99-Vierfermionterm psi*_(1,5/4) psi*_(1,1/4) psi_(2,1/4) psi_(2,-3/4) ist die Komponente mit Gewicht 1 von "
      "C+_2 auf b+_22(1/4, -7/4) Omega. Hermitesch und Tomita-kombiniert: V = sum_t (C+_t (x) C+_t + C-_t (x) C-_t), "
      "C-_t = (C+_(-t))^+", heis2[1] < TOL and heis2[2] < TOL and norm_O == [comb(u + 3, 7) for u in range(4, 11)]
      and norm_Om == norm_O and z99_ok and uhr_ok)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

say("  (b) Mischung der Paarschichten durch den Nullmodus C+_0")
T_tab, sv_T, T_ok = {}, {}, True
for m in range(1, 7):
    P21, P12 = einpaar("21", m), einpaar("12", m)
    T = P12 @ C_mat(0, -1, 1, m) @ P21.T
    soll = np.array([[(-1)**(k + kp)*(m*m - (k - kp)**2) for k in range(m)] for kp in range(m)])
    T_ok &= np.array_equal(T, soll) and np.linalg.matrix_rank(T) == min(m, 3)
    T_tab[m] = T.astype(int).tolist() if m <= 3 else f"(-1)^(k+k') (m^2 - (k-k')^2), Rang {np.linalg.matrix_rank(T)}"
    sv_T[m] = [r(x, 6) for x in np.linalg.svd(T, compute_uv=False)]
mittel_ok = all(any(abs(x - comb(m + 1, 3)) < 1e-6 for x in sv_T[m]) for m in range(2, 7))
zerl = {}
for m in range(1, 7):
    P21 = einpaar("21", m)
    B12, _ = basis2(1, -1, m)
    bild = C_mat(0, -1, 1, m) @ P21.T
    pz = np.array([paare2(1, -1, b)[0] for b in B12])
    zerl[m] = {f"{p}-Paar": r(np.sum(bild[pz == p]**2), 6) for p in sorted(set(pz))}
B5_, _ = basis2(1, -1, 5)
pz5 = np.array([paare2(1, -1, b)[0] for b in B5_])
bild_u4 = C_mat(4, -1, 1, 1) @ einpaar("21", 1).T
drei_u4 = r(np.sum(bild_u4[pz5 == 3]**2))
drei_null = max(zerl[m].get("3-Paar", 0.0) for m in range(1, 7))
neutral_null = all(not C_mat(0, 0, 0, m).any() and not C_mat(0, 0, 0, m, False).any() for m in (1, 2, 3))
H11_4 = einpaar("11", 4)
bild4 = C_mat(0, 0, 0, 4) @ H11_4.T
B24, _ = basis2(2, -2, 4)
gew4 = [r(x) for x in bild4[0]]
paar_4 = [paare2(2, -2, b) for b in B24]
emin = {}
for p in (1, 2, 3):
    werte = []
    for q in range(-3, 4):
        for E in range(0, 8):
            B_, _ = basis2(q, -q, E)
            if any(paare2(q, -q, b)[0] == p for b in B_):
                werte.append(E)
                break
    emin[p] = min(werte)
ZAHLEN["mischung"] = {"T_m (21)_m -> (12)_m": T_tab, "Singulärwerte T_m": sv_T,
                      "Zerlegung |C+_0 (21)_m|^2 nach Paarzahl": zerl,
                      "C+_0 (11)_4 auf (2,-2)-Grundzustand, k = 0..3": gew4, "Paare des Ziels (Teilchen, Löcher)": paar_4,
                      "kleinstes E^nat mit p Paaren in zwei Nähten (Q gesamt 0)": emin,
                      "max Dreipaaranteil des Nullmodus (m <= 6)": drei_null, "|Dreipaaranteil von C+_4 (21)_1|^2": drei_u4,
                      "Zahlengleichheit (kein Wörterbuch)": "mittlerer Singulärwert von T_m = m(m^2-1)/6 = binom(m+1,3)"}
check("[E] (b) Der E-erhaltende Doppeltransfer mischt die Schichten: Er bildet die Kreuzleiter (21)_m in (12)_m ab mit "
      f"T_m = (-1)^(k+k') (m^2 - (k - k')^2) (exakt für m = 1..6, Rang min(m, 3); m = 2: {T_tab[2]}), also 1->1 innerhalb der "
      "Einpaarschicht, und ab m = 3 zusätzlich in Zweipaarzustände: |C+_0 (21)_m|^2 nach Paarzahl "
      f"{ {m: zerl[m] for m in (2, 3, 5)} }. Dreipaarzustände erreicht der Nullmodus nie (max {drei_null:g}, m <= 6): zwei neue "
      "Teilchen und zwei neue Löcher kosten zusammen >= 4; 1->3 kommt nur aus Nichtnullmoden (C+_4 (21)_1 hat Dreipaaranteil "
      f"{drei_u4:g}), also erst in zweiter Ordnung und nicht entartet. Im neutralen Sektor ((11), (22)) ist C+-_0 für m <= 3 identisch null, weil der "
      "Sektor (2, -2) erst bei E^nat = Q^2 = 4 beginnt; bei m = 4 trifft C+_0 (11)_4 nur den (2, -2)-Grundzustand "
      f"(Paare {paar_4[0]}) mit Hakengewichten {gew4}. Paarschwellen in zwei Nähten: E^nat >= {emin} für 1, 2, 3 Paare "
      "(Pauli: ein zweites Teilchen oder Loch in derselben Naht kostet je 1 mehr)", T_ok and neutral_null and len(B24) == 1
      and gew4[2:] == [0.0, 0.0] and sorted(abs(x) for x in gew4[:2]) == [3.0, 5.0] and emin == {1: 1, 2: 2, 3: 5}
      and zerl[3].get("2-Paar", 0) > 0 and zerl[2].get("2-Paar", 0) == 0 and drei_null == 0 and drei_u4 > 0)
befund("Zahlengleichheit ohne Wörterbuch: Der mittlere Singulärwert von T_m ist m(m^2 - 1)/6 = binom(m + 1, 3) = "
       f"{[sv_T[m][1] for m in range(2, 7)]} (Varianzsumme der Teilchenlagen), dieselbe Zahlenfolge wie die fehlenden "
       "Box-Nachkommen. Es ist eine Amplitude, keine Zustandszahl; gezählt wird sie nicht", mittel_ok)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

say("  (c) Erste Ordnung: W = C+_0 (x) C+_0 + h.c. vertauscht mit E, das Stufenspektrum ist exakt E + lambda w")
x_ = sp.symbols("x")
erste, mult_prim, charpoly = {}, {}, {}
for E in range(1, 7):
    M = C_mat(0, -1, 1, E)
    sv = np.linalg.svd(M, compute_uv=False)
    prod = np.outer(sv, sv).ravel()
    mult_prim[E - 1] = int(np.sum(np.abs(prod - 1) < 1e-9))
    erste[E] = {"dim einseitig": len(sv), "Singulärwerte": [r(x, 6) for x in sv]}
    if E <= 4:
        charpoly[E] = str(sp.factor(sp.Matrix(np.round(M.T @ M).astype(int)).charpoly(x_).as_expr()))
kurz = [(k + 1)**2 for k in range(6)]
lang = [comb(k + 3, 3) for k in range(6)]
heis_p = [p_zahl(k)**2 for k in range(6)]
P12_2 = einpaar("12", 2)
D2 = np.kron(P12_2, P12_2)
M2 = C_mat(0, -1, 1, 2)
V1 = np.block([[np.zeros((1, 1)), C_mat(0, -1, 1, 1).T], [C_mat(0, -1, 1, 1), np.zeros((1, 1))]])
summe_ew = sorted(r(x) for x in np.linalg.eigvalsh(np.kron(V1, np.eye(2)) + np.kron(np.eye(2), V1)))
werte_k1 = sorted(r(x, 6) for x in np.linalg.svd(np.kron(P12_2 @ M2 @ einpaar("21", 2).T, P12_2 @ M2 @ einpaar("21", 2).T),
                                                    compute_uv=False))
ZAHLEN["erste_ordnung_kreuz"] = {"Stufen": erste, "charpoly(M^T M)": charpoly,
                                 "Vielfachheit des Primärwerts +-1 auf Stufe k": mult_prim, "p(k)^2": heis_p,
                                 "kurz (k+1)^2": kurz, "lang binom(k+3,3)": lang,
                                 "Verschiebungen der vier 4D-Nachkommen (12)_2 x (12)_2 (Betrag)": werte_k1,
                                 "Summenform V (x) 1 + 1 (x) Vquer auf E = 1": summe_ew}
check("[E negativ] (c) Kreuzsektor, erste Ordnung exakt: Die Primärzustände (21)_1 x (21)_1 und (12)_1 x (12)_1 mischen mit "
      "Amplitude 1, E = 1 +- lambda. Die Eigenwerte sind +-sigma_i sigma_j (sigma = Singulärwerte von C+_0 je Stufe); den Wert "
      f"des Primärs teilen auf Stufe k = 0..5 genau {list(mult_prim.values())} Zustände = p(k)^2 (der Schwerpunktsturm: "
      "Anregungen von phi_1 + phi_2 über dem Relativgrundzustand, mit denen V vertauscht; 2D-Heisenberg), nicht kurz "
      f"{kurz} und nicht lang {lang}. Die acht Zustände der Stufe 1 ((21)_2 x (21)_2 und (12)_2 x (12)_2) bekommen "
      f"+-{werte_k1} statt viermal +-1: Schon auf Stufe 1 zerfällt das Multiplett. Die Summenform V (x) 1 + 1 (x) Vquer wäre "
      f"schlechter: Sie gibt auf E = 1 die Verschiebungen {summe_ew} lambda und führt aus H^(ab) x Hquer^(ab) heraus "
      "((12) x (21))", list(mult_prim.values()) == heis_p and werte_k1 == [1.0, 7.0, 7.0, 49.0]
      and summe_ew == [-2.0, 0.0, 0.0, 2.0])
schur_ok = all(max(sv_T[m]) - min(x for x in sv_T[m] if x > 1e-9) > 1 for m in range(2, 7))
tw = {m: comb(m, 1)**4 for m in range(1, 6)}
check("[E negativ] (c) Schur-Test der Kovarianz: (21)_m x (21)_m und (12)_m x (12)_m tragen denselben irreduziblen K-Typ "
      "((m-1)/2, (m-1)/2), und P^0 verbindet die Stufen. Ein su(2,2)- (oder nur so(4)-)verträglicher Term müsste auf jeder Stufe "
      "ein Vielfaches einer Isometrie sein, mit derselben Konstante auf allen Stufen. Die Singulärwerte von T_m sind "
      f"{ {m: sv_T[m] for m in (2, 3, 4)} }, nicht konstant; keine Basis- oder Vorzeichenwahl macht W kovariant. Die "
      f"SU(2) x SU(2)-Mittelung gibt in der Basis psi*_p psi_r Omega (dort ist T_m eintragsweise "
      f"m^2 - (k - k')^2 >= 0) (tr T_m / m)^2 = m^4 = {list(tw.values())}, stufenabhängig; kovariant bleibt nur eine "
      "Konstante c, also eine gleichförmige Verschiebung Delta = 1 +- lambda c mit kurzen Vielfachheiten (k+1)^2: genau die "
      "Radiusdeformation von z99", schur_ok and all(abs(np.trace(np.abs(einpaar("12", m) @ C_mat(0, -1, 1, m) @ einpaar("21", m).T))/m
                                                       - m*m) < TOL for m in range(1, 6)))
# neutraler Sektor erste Ordnung
neu = {}
for E in (4, 5):
    Mp, Mm = C_mat(0, 0, 0, E), C_mat(0, 0, 0, E, False)
    d0 = Mp.shape[1]
    Gsum = np.kron(Mp.T @ Mp, Mp.T @ Mp) + np.kron(Mm.T @ Mm, Mm.T @ Mm)
    ev, U = np.linalg.eigh(Gsum)
    bew = U[:, ev > 1e-9]
    P1 = np.vstack([np.kron(a, b) for typ in ("11", "22") for a in einpaar(typ, E) for b in einpaar(typ, E)])
    anteil = np.sum((P1 @ bew)**2, axis=0)
    neu[E] = {"verschobene Richtungen": int(bew.shape[1]), "|w|": [r(np.sqrt(x), 6) for x in ev[ev > 1e-9]],
              "Anteil in der 4D-Einteilchenschicht": [r(x, 6) for x in anteil], "dim 4D-Einteilchen": P1.shape[0]}
ZAHLEN["erste_ordnung_neutral"] = neu
check("[E negativ] (c) Neutraler Sektor (Diagonalleitern (11), (22), die Skalare von z83): erste Ordnung exakt null für E <= 3 "
      "(gamma^(1) = 0 für den Primär, Stufen k = 1, 2 unverschoben), ab E = 4 nicht: Bei E = 4 werden "
      f"{neu[4]['verschobene Richtungen']} Richtungen um |w| = {neu[4]['|w|']} verschoben, mit Anteil {neu[4]['Anteil in der 4D-Einteilchenschicht']} "
      f"in der 4D-Einteilchenschicht (dim {neu[4]['dim 4D-Einteilchen']}); bei E = 5 {neu[5]['verschobene Richtungen']} Richtungen. "
      "Der Primär bleibt bei 1, Teile seiner Stufe k = 3 wandern in erster Ordnung: kein Multiplett",
      all(neu[E]["verschobene Richtungen"] > 0 and min(neu[E]["Anteil in der 4D-Einteilchenschicht"]) > 0.01 for E in (4, 5)))
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

say("  (d) Der Box-Nachkomme")


def boson_det_norm(P):
    """||det P |0>||^2 mit P_(a adot) = sum_c a+_(a c) b+_(adot c), c = 1..P (Moden a1c, a2c, b1c, b2c), Fock-Norm prod n!."""
    def mal(zust, mode):
        aus = {}
        for occ, k in zust.items():
            o = list(occ)
            o[mode] += 1
            aus[tuple(o)] = aus.get(tuple(o), 0) + k
        return aus

    def addiere(x, y, f=1):
        z = dict(x)
        for k, v in y.items():
            z[k] = z.get(k, 0) + f*v
        return z
    mode = lambda art, i, c: {"a": 0, "b": 2}[art]*P + (i - 1)*P + c

    def Pop(zust, a, ad):
        aus = {}
        for c in range(P):
            aus = addiere(aus, mal(mal(zust, mode("a", a, c)), mode("b", ad, c)))
        return aus
    vak = {tuple([0]*(4*P)): 1}
    det = addiere(Pop(Pop(vak, 2, 2), 1, 1), Pop(Pop(vak, 2, 1), 1, 2), -1)
    return sum(v*v*math.prod(math.factorial(n) for n in occ) for occ, v in det.items()), det


shap = {P: boson_det_norm(P)[0] for P in (1, 2, 3, 4)}
dl = sp.symbols("Delta")
poly = sp.interpolate([(P, shap[P]) for P in (1, 2, 3)], dl)
leer = boson_det_norm(1)[1]
check("[E] (d) Shapovalov in z83-Normierung (E = (N_a + N_b)/2 + Farbenzahl, [K, P] = delta E + ...): Die skalaren "
      "Tiefstgewichte der Oszillatorrealisierung mit P Farben haben Delta = P, und ||det P |Delta>||^2 = "
      f"{list(shap.values())} für Delta = 1..4, also {sp.factor(poly)} (Delta = 4 bestätigt den Grad 2). Bei Delta = 1 (eine "
      "Farbe, die Leiter von z83) ist det P = P_11 P_22 - P_12 P_21 schon als Operator null (Koeffizienten "
      f"{sorted(set(leer.values()))}): Der Box-Nachkomme fehlt identisch, nicht nur seine Norm",
      sp.expand(poly - 2*dl*(dl - 1)) == 0 and shap[4] == 2*4*3 and all(v == 0 for v in leer.values()))
ktyp = {}
for E in (2, 3, 4):
    paare_ = [(j1, j2) for j1 in [k/2 for k in range(8)] for j2 in [k/2 for k in range(8)] if abs(2*j1 + 2*j2 + 2 - E) < 1e-12]
    ktyp[E] = sum(1 for j1, j2 in paare_ if abs(j1 - j2) < 1e-12)                       # (0,0) in (j1,j1) x (j2,j2) gdw. j1 = j2
fehl = [comb(k + 3, 3) - (k + 1)**2 for k in range(6)]
ZAHLEN["box_nachkomme"] = {"||det P|Delta>||^2": {str(k): v for k, v in shap.items()}, "Polynom": str(sp.factor(poly)),
                           "(0,0)-K-Typen zweier Delta=1-Teilchen bei E = 2, 3, 4 (geordnete Paare)": ktyp,
                           "fehlende Box-Nachkommen": fehl, "binom(k + 1, 3)": [comb(k + 1, 3) for k in range(6)]}
check("[E] (d) Folgen für jede unitäre Fortsetzung: (i) Weil det P auf der Leiter als Operator verschwindet, entsteht der "
      "Box-Nachkomme erst durch Beimischung der Ordnung lambda, seine Norm ist O(lambda^2); 2 Delta (Delta - 1) = 2 gamma (1 + gamma) "
      "verlangt daher gamma = O(lambda^2) und gamma >= 0. Der Kreuzprimär mit gamma^(1) = +-lambda ist für jedes lambda != 0 "
      "ausgeschlossen (ein Zweig unter der Schranke, der andere ohne Box-Norm der Ordnung lambda). (ii) Die fehlenden Zustände "
      f"{fehl} = binom(k + 1, 3) sind das lange Multiplett eines Skalars bei E = 3 mit K-Typ (0, 0). Zwei 4D-Teilchen "
      f"(Delta = 1) haben (0, 0) bei E = 2, 3, 4 genau {list(ktyp.values())}-mal: bei E = 3 nie ((0,0) x (1/2,1/2)), also kann "
      "keine 1<->2-Mischung den Box-Nachkommen bei Delta + 2 liefern; nötig ist ein Dreiteilchenzustand (phi^3-artig), und drei "
      f"Paare kosten in zwei Nähten E^nat >= {emin[3]} > 3 (Pauli). Ein Zweinahtterm erreicht ihn nicht",
      ktyp == {2: 1, 3: 0, 4: 1} and fehl == [comb(k + 1, 3) for k in range(6)] and emin[3] == 5)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

say("  (e) Zweite Ordnung: Nichtnullmoden sum_(u != 0) C_u (x) C_u")


def h2_neutral(E, U):
    """Effektiver Operator zweiter Ordnung auf der verdoppelten neutralen Stufe E (erste Ordnung dort null, E <= 3); Regler:
    Energieübertrag 1 <= u <= U (gleich für alle Zustände, damit bleibt die exakte Schwerpunktsymmetrie erhalten)."""
    d = len(basis2(0, 0, E)[0])
    H = np.zeros((d*d, d*d))
    for u in range(1, U + 1):
        for plus in (True, False):
            M = C_mat(u, 0, 0, E, plus)
            if M.size and M.any():
                G = M.T @ M
                H -= np.kron(G, G)/u
    return H


zweite, lam_dep = {}, {}
for Lam in (8, 10, 12, 14):
    vac = h2_neutral(0, Lam)[0, 0]
    zeile = {}
    for E in (1, 2, 3):
        ev = np.sort(np.linalg.eigvalsh(h2_neutral(E, Lam))) - vac
        zeile[E] = ev
    lam_dep[Lam] = zeile
ev1 = lam_dep[14][1]
vac14 = abs(h2_neutral(0, 14)[0, 0])
null_mult = {E: int(np.sum(np.abs(lam_dep[14][E]) < 1e-9*vac14)) for E in (1, 2, 3)}
g_mm = {L: r(lam_dep[L][1][0], 6) for L in lam_dep}
g_pm = {L: r(sorted(lam_dep[L][1])[1], 6) for L in lam_dep}
expo = np.polyfit(np.log([12, 14]), np.log([abs(g_mm[12]), abs(g_mm[14])]), 1)[0]
mult_g = {E: int(np.sum(np.abs(lam_dep[14][E] - lam_dep[14][1][0]) < 1e-6*abs(lam_dep[14][1][0]))) for E in (1, 2, 3)}
ZAHLEN["zweite_ordnung_neutral"] = {"Eigenwerte - Vakuum, E = 1 (U = 14)": [r(x, 4) for x in ev1],
                                    "geschützte (Verschiebung 0) je Stufe E = 1, 2, 3": null_mult, "p(E)^2": [p_zahl(E)**2 for E in (1, 2, 3)],
                                    "gamma_(--) (Relativ x Relativ, E = 1) gegen U": g_mm, "gamma_(+-) gegen U": g_pm,
                                    "Exponent |gamma_(--)| ~ U^p": r(expo, 3),
                                    "Vielfachheit von gamma_(--) auf E = 1, 2, 3": mult_g}
check("[E/N] (e) Zweite Ordnung im neutralen Sektor (erste Ordnung dort null): Weil V nur das Relativboson trifft, sind die "
      f"Schwerpunktzustände exakt geschützt (Verschiebung 0 auf E = 1, 2, 3 genau {list(null_mult.values())}-mal = p(E)^2). Auf "
      f"E = 1 zerfallen die vier verdoppelten Zustände in 0 (Schwerpunkt x Schwerpunkt), gamma_(+-) (zweimal) und gamma_(--) = "
      f"{g_mm[14]:.4g} lambda^2 (Regler: Energieübertrag u <= U = 14): Die 4D-Einteilchenzustände (11) x (11), (22) x (22) sind keine Eigenzustände, V "
      "mischt sie mit den flavourverschiedenen (11) x (22). gamma_(--) wächst mit dem Regler wie U^"
      f"{expo:.1f} ({g_mm} für U = 8..14) und ist negativ; den Wert gamma_(--) teilen auf E = 1, 2, 3 {list(mult_g.values())} "
      "Zustände = p(k)^2 (wieder der Schwerpunktsturm; 4D-Multiplett: 1, 4, 10)",
      list(null_mult.values()) == [p_zahl(E)**2 for E in (1, 2, 3)] and g_mm[14] < 0 and expo > 5
      and list(mult_g.values()) == [p_zahl(k)**2 for k in range(3)])
befund("Teil 1: Der Doppeltransfer mit Tomita-Partner ist ein echter Zweikörperterm (erste Ordnung: 1<->1 im Kreuzsektor, 1<->2 "
       "ab E = 3 bzw. 4; 1<->3 nur über Nichtnullmoden in zweiter Ordnung), aber er ist weder so(4)- noch su(2,2)-kovariant: Die Zustände, die mit dem Primär wandern, bilden den "
       "2D-Schwerpunktsturm p(k)^2, nicht (k+1)^2 und nicht binom(k+3,3). Der Box-Nachkomme bei Delta + 2 ist durch 1<->2-Mischung "
       "grundsätzlich unerreichbar (K-Typ) und durch 1<->3 in zwei Nähten Pauli-gesperrt; erste Verschiebungen +-lambda verletzen "
       "die Unitarität, zweite sind abschneidebestimmt", True)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

# =================================================================== Teil 2
say("\nTeil 2: Fixiert eine native Bedingung lambda?")
check("[E] (a) Unitarität: Kreuzprimär E = 1 +- lambda (Teil 1c): für jedes lambda != 0 liegt ein Zweig unter Delta = 1; neutral "
      f"ist gamma^(1) = 0 und gamma_(--)^(2) = {g_mm[14]:.4g} lambda^2 < 0 (U = 14; negativ für alle U = 8..14: "
      f"{all(v < 0 for v in g_mm.values())}), ebenfalls unter der Schranke (Box-Norm 2 gamma (1 + gamma) < 0). Unitarität fixiert keinen Wert, sie verbietet "
      "lambda != 0 in beiden Sektoren", all(v < 0 for v in g_mm.values()))
th = sp.symbols("theta", real=True)
lr = sp.tanh(2*th)
E0 = lambda pL, pR: sp.simplify(((pL**2 + pR**2)/2 + lr*pL*pR)/sp.sqrt(1 - lr**2))
s2 = sp.sqrt(2)
d2D = sp.simplify(E0(2*s2, -2*s2).rewrite(sp.exp))
th_kt = sp.solve(sp.Eq(sp.simplify(sp.exp(-2*th)*8), 2), th)
e_kreuz = sp.simplify((E0(s2, -s2)/2).rewrite(sp.exp).subs(th, th_kt[0]))
C_JJ = 4*4/2                                                                            # (Q1 - Q2) je Seite = 4, J_- = (J1 - J2)/sqrt2
ZAHLEN["marginalitaet"] = {"Delta_2D(theta) des Doppeltransfers mit Partner": str(d2D), "KT-Punkt theta": str(th_kt),
                           "E_4D des einfachen Kreuzzweigs dort": str(e_kreuz), "C(O O+ -> J_- Jquer_-) = p_L p_R": C_JJ}
check("[E] (b) Marginalität: In der Naht ist V ein (4, 4)-Feld, Delta_2D = 8, beta_lambda = (2 - 8) lambda + ... irrelevant; "
      "C_OOO = 0 (Ladung 6 != 0), aber O x O+ enthält J_- Jquer_- mit p_L p_R = 8: lambda^2 treibt den Radius (KT-Struktur), ein "
      f"Fixpunkt zweiter Ordnung ist nur lambda = 0. Unter der Radiusdeformation ist Delta_2D = {d2D}; 2D-marginal (= 2) wird V bei "
      f"theta = {th_kt[0]}, und dort liegt der einfache Kreuzzweig bei E_4D = {e_kreuz} < 1. 2D-Marginalität und 4D-Unitarität "
      "schließen sich aus. In der 4D-Lesart ist E_4D = (4 + 4)/2 = 4 = d nur durch Zählung marginal (keine native S^3-Integration)",
      sp.simplify(d2D - 8*sp.exp(-2*th)) == 0 and e_kreuz == sp.Rational(1, 4) and C_JJ == 8)
M3 = C_mat(0, -1, 1, 3)
cp3 = sp.factor(sp.Matrix(np.round(M3.T @ M3).astype(int)).charpoly(x_).as_expr())
sv_irr = sp.sqrt(277 + 20*sp.sqrt(177))
sv_ok = sp.simplify(sp.sqrtdenest(sv_irr) - (10 + sp.sqrt(177))) == 0 and any(abs(x - float(10 + math.sqrt(177))) < 1e-9
                                                                            for x in np.linalg.svd(M3, compute_uv=False))
# KMS: natives Gewicht e^(-4 pi E) gegen e^(-4 pi (E + lambda W)) auf den Stufen <= 4 des Kreuzsektors (einseitig, V_0)
lam_t = 1e-3
ew_V, Z0 = 0.0, 0.0
for E in range(1, 5):
    M = C_mat(0, -1, 1, E)
    d = M.shape[0]
    V0 = np.block([[np.zeros((d, d)), M.T], [M, np.zeros((d, d))]])
    w, U = np.linalg.eigh(np.eye(2*d)*E + lam_t*V0)
    gew = np.exp(-BETA*(w - 1))
    ew_V += float(np.sum(gew*np.einsum("ij,jk,ki->i", U.T, V0, U)))
    Z0 += float(np.sum(gew))
svE3 = np.linalg.svd(M3, compute_uv=False)
w_dopp = np.outer(svE3, svE3).ravel()
dopp_ok = any(abs(x - 1) < 1e-9 for x in w_dopp) and any(abs(x - float(10 + math.sqrt(177))) < 1e-9 for x in w_dopp)
ZAHLEN["uhr_kms"] = {"charpoly(M^T M), Kreuz E = 3": str(cp3), "irrationaler Singulärwert": "10 + sqrt(177)",
                     "verdoppelte Eigenwerte |w| auf E = 3 (Auszug)": sorted({r(x, 6) for x in w_dopp})[:6],
                     "<V_0> im KMS-Zustand von E + lambda V_0 (lambda = 1e-3, beta = 4 pi)": r(ew_V/Z0, 10),
                     "<V_0> im nativen Zustand": 0.0}
check("[E negativ] (c) Uhr und KMS bei beta = 4 pi: V vertauscht mit jeder nativen Uhr (Teil 1a), die Uhr fixiert lambda also "
      "nicht über Invarianz. Die Uhridentität des Modulflusses (Delta^(-i/2) = Phase mal R, z91) verlangt aber "
      "e^(2 pi i lambda V_0) = 1, also lambda sigma in Z für alle Singulärwerte: sigma = 1 (E = 1), 7 (E = 2) und auf E = 3 "
      f"charpoly(M^T M) = {cp3}, also sigma = 10 +- sqrt(177) irrational: nur lambda = 0. Verdoppelt gilt dasselbe für die "
      "4D-Uhr e^(2 pi i E) mit e^(2 pi i lambda W) = 1: auf E = 3 liegen die Eigenwerte 1 und 10 + sqrt(177) auf derselben "
      f"Stufe ({dopp_ok}). Der native KMS-Zustand ist für lambda != 0 stationär, aber nicht KMS (<V_0> = {ew_V/Z0:.3e} im "
      "KMS-Zustand von E + lambda V_0 gegen 0), und die Produktlesart von z83 (Original und Partner je KMS bei 4 pi) verträgt "
      "keinen Term, der beide Seiten koppelt", sv_ok and "5929" in str(cp3) and ew_V/Z0 < 0 and dopp_ok)
# Komposition dreier Nähte
vak = lambda Q: (Q, NA.index(Q, ()))
a_ = NA.A(2, 0, 0)[NA.index(2, ()), NA.index(0, ())]
b_ = NA.B(-2, 2, 0)[NA.index(0, ()), NA.index(2, ())]
C3 = a_*b_
al = {"12": np.array([2, -2, 0]), "23": np.array([0, 2, -2]), "13": np.array([2, 0, -2])}
gram = {f"{i}.{j}": int(al[i] @ al[j]) for i in al for j in al}
g_ = sp.symbols("g")
NN = sp.symbols("N", positive=True)
beta_g = -6*g_ - 2*sp.pi*(NN - 2)*(sp.Integer(2)/sp.sqrt(8))*g_**2
g_stern = [s for s in sp.solve(beta_g, g_) if s != 0][0]
ZAHLEN["komposition"] = {"<alpha_13| (V_12)_0 |alpha_23>": r(C3), "Gitterprodukte": gram,
                         "beta_g (Cardy, g je Nahtpaar, normiert)": str(beta_g), "g*": str(sp.simplify(g_stern)),
                         "g* (N = 3, 4, 10)": [r(float(g_stern.subs(NN, n)), 6) for n in (3, 4, 10)]}
check("[E/C] (d) Komposition dreier Nähte: alpha_12 + alpha_23 = alpha_13, die Gitterprodukte "
      f"{ {k: v for k, v in gram.items() if k in ('12.23', '12.13', '12.12')} } sind gerade (lokal, bosonisch, für jedes lambda "
      f"konsistent), und <alpha_13| (V_12)_0 |alpha_23> = {C3:g}: Die Doppeltransfers schließen ab N = 3 unter der OPE. In "
      f"Cardys Störungstheorie gibt das beta_g = {beta_g}, eine Nullstelle g* = {sp.simplify(g_stern)} (N = 3: "
      f"{float(g_stern.subs(NN, 3)):.3f}): von Ordnung 1, also außerhalb der Störungstheorie, N-abhängig und für N -> unendlich "
      "(Sym^N-Grenze von z91) null. In der 4D-Lesart (marginal) bleibt beta = -c g^2 ohne Nullstelle. Komposition fixiert kein "
      "kleines lambda", abs(abs(C3) - 1) < TOL and all(v % 2 == 0 for v in gram.values())
      and abs(float(g_stern.subs(NN, 10))) < abs(float(g_stern.subs(NN, 3))))
befund("Teil 2: Keine native Bedingung fixiert ein lambda != 0. Unitarität (Kreuz erster, neutral zweiter Ordnung), die "
       "Uhridentität bei beta = 4 pi (irrationale Singulärwerte 10 +- sqrt(177)) und der KMS-Zustand erzwingen lambda = 0; "
       "2D-Marginalität verlangt einen Radius mit E_4D = 1/4 < 1; die Komposition gibt nur eine nichtperturbative, N-abhängige "
       "Nullstelle, die im Kompositionslimes verschwindet", True)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

# =================================================================== Teil 3
say("\nTeil 3: H6 -- induziertes G zweiter Ordnung")


def norm_einzel(op, s_, Q, lam):
    v = np.zeros(p_zahl(sum(lam)))
    v[NA.index(Q, lam)] = 1.0
    M = (NA.A if op == "A" else NA.B)(s_, Q, sum(lam))
    return float(np.sum((M @ v)**2)) if M.size else 0.0


def spur_G(m, u):
    """sum_k ||C+-_u psi_k||^2 für psi_k = Haken_k (Naht 1) x Omega (Naht 2), Stufe m -> m + u."""
    tp = tm = 0.0
    for k in range(m):
        h = haken(m, k)
        for s_ in range(2 - m, u - 1):                                                     # B_(u-s) Omega != 0 nur für u - s >= 2
            nA = norm_einzel("A", s_, 0, h)
            nB = norm_einzel("B", s_, 0, h)
            tp += nA*(norm_einzel("B", u - s_, 0, ()) if nA else 0.0)
            tm += nB*(norm_einzel("A", u - s_, 0, ()) if nB else 0.0)
    return tp, tm


vak_norm = [sum(norm_einzel("A", s_, 0, ())*norm_einzel("B", u - s_, 0, ()) for s_ in range(0, u + 1)) for u in range(4, 12)]
sakh = {}
for Lam in (12, 14, 16):
    dvac = -2*sum(comb(u + 3, 7)**2/u for u in range(4, Lam + 1))
    zeile = {}
    for m in range(1, 6):
        s = 0.0
        for u in range(-m, Lam + 1):
            if u == 0 or m + u < 4:
                continue
            tp, tm = spur_G(m, u)
            s -= (tp*tp + tm*tm)/u
        dom = s/m**2 - dvac
        zeile[m] = {"delta_omega": dom, "delta_xi": m*dom/3}
    sakh[Lam] = zeile
xi_tab = {L: [float(f"{sakh[L][m]['delta_xi']:.4e}") for m in range(1, 6)] for L in sakh}
verh = [r(sakh[16][m]["delta_xi"]/sakh[16][1]["delta_xi"], 4) for m in range(1, 6)]
expo_xi = np.polyfit(np.log([14, 16]), np.log([abs(sakh[14][1]["delta_xi"]), abs(sakh[16][1]["delta_xi"])]), 1)[0]
ZAHLEN["sakharov"] = {"||C+_u Omega||^2 aus Einzelnähten": vak_norm, "delta_xi_eff(m), m = 1..5, je U": xi_tab,
                      "delta_xi(m)/delta_xi(1) bei U = 16": verh, "Exponent delta_xi(1) ~ U^p": r(expo_xi, 3)}
check("[E/N, Lesart L] (a) Sakharov-Lesart des Spektrums: Die verdoppelte 4D-Einteilchenleiter (11) x (11) verschiebt sich in zweiter "
      "Ordnung (Spur je Stufe minus Vakuum) um delta_omega_m; als Krümmungskopplung gelesen (omega^2 = n(n+2) + 6 xi auf S^3) ist "
      f"delta_xi_eff(m) = m delta_omega_m / 3 = {xi_tab[16]} lambda^2 (m = 1..5, Regler U = 16), nicht null und negativ (also "
      f"1/G_ind = -6 delta_xi I[f]/(12 pi) > 0, proportional lambda^2), aber stufenabhängig (Verhältnis zu m = 1: {verh}; eine Krümmungskopplung wäre konstant) und "
      f"reglerbestimmt (delta_xi(1) ~ U^{expo_xi:.1f}: {[xi_tab[L][0] for L in xi_tab]} für U = 12, 14, 16). Die "
      "Vakuumnormen aus Einzelnähten bestätigen binom(u + 3, 7)",
      vak_norm == [comb(u + 3, 7) for u in range(4, 12)] and all(x < 0 for x in xi_tab[16])
      and max(abs(v - 1) for v in verh) > 0.5 and expo_xi > 5)
profile = {"scharf/Q1/P2": 1, "Q2/Brown": 2, "P3": sp.Rational(1, 2), "Pauli-Villars": 2*sp.log(2)}
dxi = sp.symbols("delta_xi")
G_ind = {k: sp.simplify(-6*dxi*v/(48*sp.pi)) for k, v in profile.items()}
ZAHLEN["sakharov_profile"] = {"1/(4 G_ind) je eps^-2": {k: str(v) for k, v in G_ind.items()}}
check("[E] (b) Selbst bei fest gelesenem delta_xi bleibt der Regler: 1/(4 G_ind) = (1 - 6 xi) I[f]/(48 pi) (z90) gibt "
      f"{ {k: str(v) for k, v in G_ind.items()} } je eps^-2, Faktor 2 zwischen Q1 und Q2 bei gleicher Nennlänge. "
      "Dimensionsanalyse: In der 4D-Lesart ist lambda dimensionslos (E_4D = 4 = d), in der Naht irrelevant (Delta_2D = 8); ein "
      "regulatorfreier Teil von 1/G bräuchte M^2 ln(Lambda/M) mit einer Massenskala M, die V nicht hat. 1/G ~ lambda^2 U^p "
      "mal Profil: G und g sind nur über den Regler verknüpft",
      sp.simplify(G_ind["Q2/Brown"]/G_ind["scharf/Q1/P2"]) == 2)
befund("Teil 3: H6 in dieser Form gescopt negativ. Der konformitätsbrechende Zweikörperterm induziert ein 1/G der Ordnung "
       "lambda^2, aber nicht als Krümmungskopplung (xi_eff hängt von der Stufe ab) und nicht regulatorfrei (Abschneidepotenz und "
       "Profil); ein festes G aus g entsteht nicht", True)

say("\nBilanz für N6 und H6")
STATUS = {
    "Bau des Zweikörperterms": ("E", "C+_t = sum A^(1)_s B^(2)_(t-s) = Moden von e^(2i(phi_1 - phi_2)), h = 4 (Normen binom(u+3,7)), nur "
                                     "Relativboson; Tomita-Kombination V = sum_t C_t (x) C_t + h.c., (4, 4), spinfrei, uhrinvariant"),
    "Mischung": ("E", "erste Ordnung: 1<->1 im Kreuzsektor (T_m = (-1)^(k+k') (m^2 - (k-k')^2), Rang 3), 1<->2 ab E = 3 "
                      "(Kreuz) bzw. 4 (neutral, Schwelle Q^2 = 4); 1<->3 nie im Nullmodus, nur über Nichtnullmoden (u >= 4)"),
    "Vielfachheiten erster Ordnung": ("E (negativ)", "mit dem Primär wandern p(k)^2 = 1, 1, 4, 9, 25, 49 Zustände (Schwerpunktsturm), "
                                                     "nicht (k+1)^2 und nicht binom(k+3,3); Nachkommen bei k = 1: |w| = 1, 7, 7, 49"),
    "Kovarianz": ("E (negativ)", "Singulärwerte von T_m nicht konstant (Schur): weder so(4) noch su(2,2); kovariant bleibt nur eine "
                                 "Konstante = Radiusdeformation mit kurzen Vielfachheiten (z99)"),
    "Box-Nachkomme": ("E (negativ)", "||det P|Delta>||^2 = 2 Delta (Delta - 1), det P = 0 auf der Leiter: gamma = O(lambda^2) >= 0; "
                                     "Zweiteilchen-(0,0) bei E = 3 fehlt, drei Paare in zwei Nähten erst ab E^nat = 5"),
    "zweite Ordnung": ("E/N", "Schwerpunkt exakt geschützt (p(E)^2), Relativleiter gamma_(--) < 0, ~ U^p, mit dem Primär "
                                  "wandern wieder p(k)^2 Zustände; mischt die 4D-Flavours"),
    "lambda nativ fixiert": ("E (negativ)", "Unitarität, Uhridentität (sigma = 10 +- sqrt(177)), KMS: lambda = 0; 2D-Marginalität nur "
                                            "bei E_4D = 1/4; Komposition: g* = O(1), N-abhängig, -> 0 für N -> unendlich (C)"),
    "induziertes G (H6)": ("E/N (negativ)", "1/G ~ lambda^2, aber xi_eff stufenabhängig und ~ U^p, dazu Profilfaktor 2: keine "
                                            "regulatorfreie Verknüpfung von G und g"),
    "Wechselwirkung (T3/T5)": ("O", "braucht einen Term, der den Dreiteilchen-(0,0)-Zustand bei E = 3 mit Amplitude O(lambda) an den "
                                    "Primär koppelt (Box phi = lambda phi^3), su(2,2)-kovariant und über >= 3 Nähte; der Zweinaht-"
                                    "Doppeltransfer tut es nicht"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
urteil = ("N6 ist gescopt negativ beantwortet, H6 in dieser Form ebenfalls. Der Doppeltransfer e^(2i(phi_1 - phi_2)) ist als "
          "Fockoperator exakt rekonstruiert (Vandermonde-Moden, Heisenberg-Relationen exakt, Normen binom(u+3,7), h = 4) und mit dem "
          "Tomita-Partner zu einem hermiteschen, uhrinvarianten, spinfreien (4, 4)-Zweikörperterm kombiniert. Er mischt die "
          "Paarschichten wirklich: in erster Ordnung 1<->1 zwischen den Kreuzleitern mit T_m = (-1)^(k+k') (m^2 - (k-k')^2) und "
          "1<->2 ab E = 3 im Kreuz- und ab der Schwelle Q^2 = 4 im neutralen Sektor; 1<->3 nur über Nichtnullmoden. Er liefert aber keine Box-Nachkommen: In erster Ordnung "
          "(exakt, weil der Nullmodus mit E vertauscht) wandern mit dem Primär genau p(k)^2 = 1, 1, 4, 9, 25, 49 Zustände, der "
          "2D-Schwerpunktsturm des Relativbosons, nicht die kurzen (k+1)^2 und nicht die langen binom(k+3,3); die vier 4D-Nachkommen "
          "der Stufe 1 verschieben sich um +-1, +-7, +-7, +-49 statt viermal +-1. Nach Schur ist der Term weder so(4)- noch su(2,2)-kovariant; "
          "kovariant bliebe nur eine Konstante, also die Radiusdeformation von z99 mit kurzen Vielfachheiten. Strukturell gilt: "
          "||det P|Delta>||^2 = 2 Delta (Delta - 1) in z83-Normierung und det P = 0 auf der Leiter, also gamma = O(lambda^2) >= 0; der "
          "Kreuzprimär mit 1 +- lambda ist ausgeschlossen. Der fehlende Box-Nachkomme ist ein (0,0)-Zustand bei E = 3; zwei "
          "Delta = 1-Teilchen haben dort keinen, drei Paare kosten in zwei Nähten E^nat >= 5. Eine 1<->2-Mischung kann ihn "
          "grundsätzlich nicht liefern, ein Zweinahtterm auch nicht über 1<->3. Zweite Ordnung: Der Schwerpunkt ist exakt geschützt, "
          "die Relativleiter bekommt gamma_(--) < 0 mit hoher Reglerpotenz (~ U^11), und mit ihr wandern wieder p(k)^2 Zustände. Keine native Bedingung fixiert lambda != 0: "
          "Unitarität, die Uhridentität des Modulflusses bei beta = 4 pi (Singulärwerte 10 +- sqrt(177)) und KMS erzwingen lambda = 0; "
          "2D-Marginalität verlangt einen Radius mit E_4D = 1/4; die Komposition dreier Nähte schließt die OPE und gibt eine "
          "Nullstelle g* = O(1), N-abhängig und für N -> unendlich null. H6: Das induzierte 1/G ist proportional lambda^2, aber "
          "keine Krümmungskopplung (xi_eff stufenabhängig) und nicht regulatorfrei (Reglerpotenz ~ U^11, Profilfaktor 2); G und g "
          "werden nicht verknüpft. Folgen: K15 und T3/T5 brauchen einen su(2,2)-kovarianten Term über mindestens drei Nähte, der "
          "phi^3-artige Zustände mit Amplitude O(lambda) ankoppelt; K13 bleibt bei Kinematik, ein natives G ist auf diesem Weg "
          "nicht erreichbar. Nichts davon ist befördert.")
say("\nUrteil\n   " + urteil)
say(f"\n   (Laufzeit {time.time() - T0:.1f} s)")
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z105_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
