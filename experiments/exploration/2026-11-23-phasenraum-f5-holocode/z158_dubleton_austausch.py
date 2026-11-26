"""Z158: Wird der D(3,2)-Sektor des Dubletons zwischen zwei Nähten ausgetauscht, und wenn ja, als eigenes Teilchen
(Delta = 3) oder nur als Produkt von vier Stromquanten? Test von P_X in enger Form (z156) und Prüfung des Kriteriums.

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart bzw. Literatur, N numerisch, O offen.

Bestand: z75 (Dubleton aus zwei Strommoden a_i = J^(i)_m/sqrt m und Partnern; Linienantwort C(n + Delta - 1, 2 Delta - 1)),
z120/z133/z145 (Transfer B = e^(-i pi J_x/2), Gatter U(phi) = B^+ e^(-i phi H_Q) B = e^(-i phi (Q_+^2 + Q_-^2)/2); bei gerader
Gesamtladung Wurzel des Nahttauschs), z152 (Naht als Halbkette, Korrelationslänge 0,567), z156 (Sektor D = +-4 trägt D(3,2);
Fockcharakter (Di + Rac)^2; Kriterium 'null oder faktorisiert => P_X widerlegt'), z136/z137 (Nähte je Horizont; G als
Anzahl), z96 (l^2/G = 32 pi^3 e^(2/alpha)).

Literatur [L]: Flato-Fronsdal (1978), Howe (duale Paare), Osborn-Petkou (C_T freier Felder), Klebanov-Polyakov / Vasiliev
(freie O(N)-Randtheorie und Höhere-Spin-Gravitation), Steinberg (Stabilisatoren in Spiegelungsgruppen).

Modell (wie z133): zwei Nähte A, B, je Stromart ein chirales Dirac-Fermion mit alpha = 1/4, Niveaus x = y + 1/4, See y < 0
gefüllt. Behalten werden alle Zustände bis zu einer Gesamtenergie E_max; alle benutzten Operatorketten bleiben darunter,
es gibt also keinen Abschneidefehler. Ströme J^(ab)_(-n) = Summe_y c^+_(a,y+n) c_(b,y); J^0 = (J^A + J^B)/2,
J^z = (J^A - J^B)/2, J^x, J^y gemischt. Antwort r_k = <a_B^k Phi((a_A^+)^k)>/<a^k (a^+)^k> im Heisenberg-Bild.

Teil 1  Freier Transfer: Bild der Strommoden, exakte Zwei-Naht-Antwort (Stromalgebra und Fockraum), Quintett, thermisch,
        Kette.
Teil 2  Wechselwirkendes Gatter: Kommutatoren, Spektralzerlegung, Antwort im Vakuum und thermisch, zwei Stromarten.
Teil 3  Welche Dimension hat der Austausch: Linienantwort im Rahmen der Naht (Möbius) und im Rahmen der so(3,2).
Teil 4  Singleton-Lesart: Charaktere, Singulett, native Gruppen, zweite Quantisierung, Zählung für G.
Teil 5  Bilanz und Urteil.
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

import mpmath as mp
import numpy as np
import scipy.sparse as sps
import sympy as sp

HERE = Path(__file__).parent
REPO = HERE.parents[2]
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
TOL = 1e-9


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")
    return bool(cond)


def befund(name, val, text):
    BEFUNDE[name] = val
    say(f"   [Befund: {val}] {text}")


def gz(x, n=8):
    """Reelle Zahl mit n geltenden Ziffern, ohne negatives Nullzeichen."""
    x = float(x)
    if abs(x) < 1e-13:
        x = 0.0
    return f"{x:.{n}g}"


def gc(z, n=8):
    z = complex(z)
    re, im = (0.0 if abs(z.real) < 1e-12 else z.real), (0.0 if abs(z.imag) < 1e-12 else z.imag)
    return f"{re:.{n}g}{'+' if im >= 0 else '-'}{abs(im):.{n}g}i"


def binom(n, k):
    return math.comb(n, k) if 0 <= k <= n else 0


# ====================================================================================================================
# Fockmodell: nk Kanäle (Naht, Stromart), Niveaus y in -M..M-1, alle Zustände mit Energie <= Emax4/4
# ====================================================================================================================
def kanal(M, Emax4):
    """Zustände eines Kanals als (4 E, Maske, Ladung); Teilchen y >= 0 kosten y + 1/4, Löcher y < 0 kosten -y - 1/4."""
    exc = sorted([(4*y + 1, y, 1) for y in range(M)] + [(-4*y - 1, y, -1) for y in range(-M, 0)])
    out = []

    def rec(i, mask, E, Q):
        out.append((E, mask, Q))
        for j in range(i, len(exc)):
            e, y, q = exc[j]
            if E + e > Emax4:
                break
            rec(j + 1, mask ^ (1 << (y + M)), E + e, Q + q)

    rec(0, (1 << M) - 1, 0, 0)
    return sorted(out)


class Modell:
    def __init__(self, nk, M, Emax4):
        self.nk, self.M, self.Emax4 = nk, M, Emax4
        ein = kanal(M, Emax4)

        def spreize(mask, ch):
            s, l = 0, 0
            while mask:
                if mask & 1:
                    s |= 1 << (nk*l + ch)
                mask >>= 1
                l += 1
            return s

        zust = []

        def rec(ch, s, E, Qs):
            if ch == nk:
                zust.append((E, s, Qs))
                return
            for e, mask, q in ein:
                if E + e > Emax4:
                    break
                rec(ch + 1, s | spreize(mask, ch), E + e, Qs + (q,))

        rec(0, 0, 0, ())
        zust.sort()
        self.E4 = np.array([t[0] for t in zust])
        self.bits = [t[1] for t in zust]
        self.Q = np.array([t[2] for t in zust])
        self.Qt = self.Q.sum(axis=1)
        self.idx = {b: i for i, b in enumerate(self.bits)}
        self.dim = len(zust)
        self._cache = {}
        self.vak = np.zeros(self.dim, complex)
        self.vak[0] = 1.0

    def T(self, a, b, n):
        """Summe_y c^+_(a,y+n) c_(b,y) (ohne Diagonalterm); erhöht die Energie um n."""
        key = (a, b, n)
        if key in self._cache:
            return self._cache[key]
        if n < 0:
            out = self.T(b, a, -n).T.tocsr()
        else:
            nk, M = self.nk, self.M
            rows, cols, vals = [], [], []
            for i, s in enumerate(self.bits):
                for y in range(-M, M - n):
                    q, p = nk*(y + M) + b, nk*(y + n + M) + a
                    if p == q or not (s >> q) & 1 or (s >> p) & 1:
                        continue
                    lo, hi = (p, q) if p < q else (q, p)
                    zw = ((1 << hi) - 1) & ~((1 << (lo + 1)) - 1)
                    j = self.idx.get(s ^ (1 << q) ^ (1 << p))
                    if j is not None:
                        rows.append(j)
                        cols.append(i)
                        vals.append(-1.0 if bin(s & zw).count("1") & 1 else 1.0)
            out = sps.csr_matrix((vals, (rows, cols)), shape=(self.dim, self.dim), dtype=complex)
        self._cache[key] = out
        return out

    def fblock(self, K, f):
        """f(K) für hermitesches K, das Energie und Gesamtladung erhält (blockweise diagonalisiert)."""
        lab = {}
        for i in range(self.dim):
            lab.setdefault((int(self.E4[i]), int(self.Qt[i])), []).append(i)
        K = K.tocsr()
        rows, cols, vals = [], [], []
        for key in sorted(lab):
            ix = np.array(lab[key])
            w, v = np.linalg.eigh(K[ix][:, ix].toarray())
            F = (v*f(w)) @ v.conj().T
            r, c = np.meshgrid(ix, ix, indexing="ij")
            rows.append(r.ravel())
            cols.append(c.ravel())
            vals.append(F.ravel())
        return sps.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                              shape=(self.dim, self.dim))

    def diag(self, v):
        return sps.diags(np.asarray(v, dtype=complex)).tocsr()


def kette(ops, v):
    for O in reversed(ops):
        v = O @ v
    return v


def mnorm(A):
    return float(abs(A).max()) if A.nnz else 0.0


# ====================================================================================================================
# Exakte Stromalgebra u(2)_1 = u(1) + su(2)_1 auf dem Vakuum (Basis 0, 3, +, -), rationale Rekursion
# ====================================================================================================================
def komm(X, Y):
    """[X, Y] für X = (a, m), Y = (b, n): Liste (Koeffizient, Operator oder None für die Zahl)."""
    (a, m), (b, n) = X, Y
    out = []
    if a == "0" or b == "0":
        if a == b and m + n == 0:
            out.append((Fr(m, 2), None))
        return out
    if a == "3" and b == "3":
        if m + n == 0:
            out.append((Fr(m, 2), None))
    elif a == "3":
        out.append((Fr(1 if b == "+" else -1), (b, m + n)))
    elif b == "3":
        out.append((Fr(-1 if a == "+" else 1), (a, m + n)))
    elif a != b:
        out.append((Fr(2 if a == "+" else -2), ("3", m + n)))
        if m + n == 0:
            out.append((Fr(m), None))
    return out


@lru_cache(maxsize=None)
def vev(wort):
    if not wort:
        return Fr(1)
    if wort[-1][1] >= 0 or wort[0][1] <= 0:
        return Fr(0)
    p = max(i for i, o in enumerate(wort) if o[1] >= 0)
    tot = Fr(0)
    for j in range(p + 1, len(wort)):
        for c, op in komm(wort[p], wort[j]):
            tot += c*vev(wort[:p] + wort[p + 1:j] + ((op,) if op is not None else ()) + wort[j + 1:])
    return tot


def erw_alg(faktoren):
    tot = Fr(0)
    for combo in itertools.product(*[sorted(f.items()) for f in faktoren]):
        c = Fr(1)
        for _, x in combo:
            c *= x
        tot += c*vev(tuple(k for k, _ in combo))
    return tot


def r_frei(k, m):
    """<a_B^k (J^0 + J^x)_(-m)^k>/(m^k k!): Antwort von k Quanten unter dem freien Vierteltransfer (Vakuum)."""
    aB = {("0", m): Fr(1), ("3", m): Fr(-1)}
    aT = {("0", -m): Fr(1), ("+", -m): Fr(1, 2), ("-", -m): Fr(1, 2)}
    return erw_alg([aB]*k + [aT]*k)/(Fr(m)**k*math.factorial(k))


def G_alg(k, m):
    z = {("3", m): Fr(1)}
    x = {("+", -m): Fr(1, 2), ("-", -m): Fr(1, 2)}
    return erw_alg([z]*k + [x]*k)


# ====================================================================================================================
say("Z158: Austausch des D(3,2)-Sektors des Dubletons zwischen zwei Nähten")
say("=" * 100)
say("\nTeil 1: Der Austausch unter dem freien Transfer")

mod = Modell(2, 6, 24)                                 # zwei Nähte, eine Stromart, E <= 6
A_, B_ = 0, 1
JA = lambda n: mod.T(A_, A_, -n)                       # J^A_n senkt die Energie um n
JB = lambda n: mod.T(B_, B_, -n)
J0 = lambda n: (JA(n) + JB(n))/2
Jz = lambda n: (JA(n) - JB(n))/2
Jx = lambda n: (mod.T(A_, B_, -n) + mod.T(B_, A_, -n))/2
Jy = lambda n: (mod.T(A_, B_, -n) - mod.T(B_, A_, -n))/2j
Jz0 = mod.diag((mod.Q[:, 0] - mod.Q[:, 1])/2)
Jx0, Jy0 = Jx(0), Jy(0)
vak = mod.vak
tief = np.where(mod.E4 <= 8)[0]                         # Spalten mit E <= 2: alle Ketten bis +4 bleiben im Raum

# 1.1 Stromalgebra
res = []
JJ = {"x": Jx, "y": Jy, "z": Jz}
J00 = {"x": Jx0, "y": Jy0, "z": Jz0}
eps = {("x", "y"): ("z", 1), ("y", "z"): ("x", 1), ("z", "x"): ("y", 1),
       ("y", "x"): ("z", -1), ("z", "y"): ("x", -1), ("x", "z"): ("y", -1)}
for (a, b), (c, sg) in sorted(eps.items()):
    for m_, n_ in ((1, -1), (2, -1), (1, -2), (2, -2), (1, 1)):
        Xa, Xb = JJ[a](m_), JJ[b](n_)
        ziel = (J00[c] if m_ + n_ == 0 else JJ[c](m_ + n_))*(1j*sg)
        res.append(mnorm(((Xa @ Xb - Xb @ Xa) - ziel)[:, tief]))
for a in "xyz":
    for m_ in (1, 2):
        Xa, Xb = JJ[a](m_), JJ[a](-m_)
        res.append(mnorm(((Xa @ Xb - Xb @ Xa) - (m_/2)*sps.identity(mod.dim, format="csr"))[:, tief]))
        res.append(mnorm((J0(m_) @ Xb - Xb @ J0(m_))[:, tief]))
    res.append(mnorm(((J0(1) @ J0(-1) - J0(-1) @ J0(1)) - 0.5*sps.identity(mod.dim, format="csr"))[:, tief]))
normen = []
v = vak.copy()
for k in range(1, 5):
    v = JA(-1) @ v
    normen.append(float(np.vdot(v, v).real))
check(f"[E] Modell: zwei Nähte, {mod.dim} Zustände mit E <= 6 (abgeschlossen, kein Abschneidefehler). Die Bilineare "
      f"bilden u(2)_1 = u(1) + su(2)_1: [J^a_m, J^b_n] = i eps_abc J^c_(m+n) + (m/2) delta_ab delta_(m+n), J^0 zentral "
      f"dazu (Rest < 1e-12 auf allen Zuständen mit E <= 2); |(J^A_-1)^k|0>|^2 = {[round(x) for x in normen]} = k! "
      "(Cartan-Strom gaußsch, z58)", max(res) < 1e-12 and all(abs(x - math.factorial(k + 1)) < 1e-10
                                                           for k, x in enumerate(normen)))

# 1.2 Bild der Strommode unter dem Transfer
Bf = mod.fblock(Jx0, lambda w: np.exp(-1j*np.pi/2*w))            # B = e^(-i pi J_x/2)
Bfd = Bf.conj().T.tocsr()
d12 = max(mnorm((Bfd @ JA(-m_) @ Bf - (J0(-m_) + Jy(-m_)))[:, tief]) for m_ in (1, 2))
d12b = max(mnorm((Bf @ JA(-m_) @ Bfd - (J0(-m_) - Jy(-m_)))[:, tief]) for m_ in (1, 2))
vy = Jy(-1) @ vak
ladung_y = sorted({(int(mod.Q[i, 0]), int(mod.Q[i, 1])) for i in np.where(abs(vy) > 1e-12)[0]})
ueb = [abs(np.vdot(JA(-1) @ vak, vy)), abs(np.vdot(JB(-1) @ vak, vy))]
norm_y = float(np.vdot(vy, vy).real)
kom_zy = mnorm(((Jz(1) @ Jy(-1) - Jy(-1) @ Jz(1)) + 1j*Jx0)[:, tief])
check("[E] Bild einer Strommode: B^+ J^A_(-m) B = (J^A_(-m) + J^B_(-m))/2 + J^y_(-m) (und B J^A B^+ mit -J^y), "
      "m = 1, 2. Das ist nicht linear in den Strömen der beiden Nähte: J^y = (c_A^+ c_B - c_B^+ c_A)/2i ist ein "
      f"gemischtes Bilinear. J^y_(-1)|0> liegt ganz in den Ladungssektoren (Q_A, Q_B) = {ladung_y}, ist orthogonal zu "
      f"beiden Nahtströmen und trägt die Norm {gz(norm_y, 6)} = m/2: Von einem Stromquant bleibt 1/4 auf A, 1/4 geht "
      "nach B, 1/2 in den Transferstrom. Dieser vertauscht nicht mit den Nahtströmen: [J^z_m, J^y_(-m)] = -i J^x_0, "
      "die Transferladung selbst",
      d12 < 1e-12 and d12b < 1e-12 and ladung_y == [(-1, 1), (1, -1)] and max(ueb) < 1e-12
      and abs(norm_y - 0.5) < 1e-12 and kom_zy < 1e-12)

# 1.3 Exakte Antworten
RF = {m_: [r_frei(k, m_) for k in range(1, 5)] for m_ in range(1, 7)}
GF = {m_: [G_alg(k, m_) for k in range(1, 5)] for m_ in range(1, 7)}
geschl = all(RF[m_] == [Fr(1, 2), Fr(m_ + 1, 4*m_), Fr(m_ + 3, 8*m_), Fr(3*m_**3 + 18*m_**2 + 9*m_ - 8, 48*m_**3)]
             for m_ in RF)
geschl_G = all(GF[m_] == [0, Fr(m_, 2), 0, Fr(m_*(9*m_ - 8), 2)] for m_ in GF)
fock = {}
for m_, kmax in ((1, 4), (2, 3)):
    for k in range(1, kmax + 1):
        v = kette([J0(-m_) + Jy(-m_)]*k, vak)
        w = kette([JB(-m_)]*k, vak)
        fock[(m_, k)] = np.vdot(w, v)/(m_**k*math.factorial(k))
fock_ok = all(abs(fock[(m_, k)] - float(RF[m_][k - 1])) < 1e-12 for (m_, k) in fock)
ZAHLEN["freier Transfer r_k(m), Vakuum"] = {f"m = {m_}": [str(x) for x in RF[m_]] for m_ in RF}
check("[E] Zwei-Naht-Antwort unter dem freien Vierteltransfer im Vakuum, r_k = <a_B^k B (a_A^+)^k B^+>/k!, exakt aus der "
      f"Stromalgebra (rationale Rekursion, m = 1..6) und im Fockraum (m = 1: k <= 4; m = 2: k <= 3; Abweichung < 1e-12): "
      f"m = 1: {[str(x) for x in RF[1]]}, m = 2: {[str(x) for x in RF[2]]}, m = 3: {[str(x) for x in RF[3]]}. "
      "Geschlossen: r_1 = 1/2, r_2 = (m + 1)/(4m), r_3 = (m + 3)/(8m), r_4 = (3m^3 + 18m^2 + 9m - 8)/(48m^3). Die Antwort "
      "ist nicht null", geschl and fock_ok)
verb = {m_: [RF[m_][k - 1] - Fr(1, 2)**k for k in range(1, 5)] for m_ in RF}
check("[E] Faktorisierung: Vier unabhängige Stromantworten gäben r_1^4 = 1/16. Tatsächlich ist r_4 - r_1^4 = "
      "3/(8m) + 3/(16m^2) - 1/(6m^3) (m = 1: 19/48, m = 2: 41/192), also faktorisiert die Antwort nur für m -> unendlich. "
      "Der verbundene Anteil kommt allein aus dem gemischten Strom: <(J^z_m)^2 (J^x_-m)^2> = m/2 und "
      "<(J^z_m)^4 (J^x_-m)^4> = m(9m - 8)/2 statt 0 (ungerade Ordnungen 0). Es ist die Nichtgaußizität der "
      "nichtabelschen Niveau-1-Ströme aus z75, in den Fermionen bleibt alles Wick (der Transfer ist quasifrei)",
      geschl_G and verb[1][3] == Fr(19, 48) and verb[2][3] == Fr(41, 192)
      and all(verb[m_][3] == Fr(3, 8*m_) + Fr(3, 16*m_**2) - Fr(1, 6*m_**3) for m_ in verb))

# 1.4 Quintett unter dem freien Transfer (zwei Stromarten, Transfer je Stromart)
quint = {m_: [RF[m_][3], RF[m_][2]*RF[m_][0], RF[m_][1]**2] for m_ in (1, 2, 3, 6)}
ZAHLEN["Quintett frei (Komponenten 4+0, 3+1, 2+2)"] = {f"m = {m_}": [str(x) for x in q] for m_, q in quint.items()}
check("[E] Die fünf Spin-2-Komponenten (a_1^+)^(4-j) (a_2^+)^j (zwei Stromarten, der Transfer wirkt je Stromart): Die "
      f"Antwort ist diagonal mit r_4, r_3 r_1, r_2^2, r_1 r_3, r_4 = {[str(x) for x in quint[1]]} (m = 1), "
      f"{[str(x) for x in quint[2]]} (m = 2), {[str(x) for x in quint[3]]} (m = 3) für die ersten drei. Wäre die SU(2) "
      "des Strompaars (die Drehgruppe der so(3,2), z156) eine Symmetrie des Transfers, wären alle fünf gleich. Sie "
      "unterscheiden sich in Ordnung 1/m (Koeffizienten 3/8, 3/16, 1/8; bei m = 1 fallen die beiden letzten zufällig "
      "zusammen): Der freie Transfer überträgt das Quintett nicht als Multiplett",
      len(set(quint[1])) == 2 and all(len(set(quint[m_])) == 3 for m_ in (2, 3, 6))
      and all(quint[m_][0] - Fr(1, 16) - Fr(3, 8*m_) == Fr(3, 16*m_**2) - Fr(1, 6*m_**3) for m_ in quint)
      and all(quint[m_][1] == Fr(1, 16) + Fr(3, 16*m_) for m_ in quint)
      and all(quint[m_][2] == Fr(1, 16) + Fr(1, 8*m_) + Fr(1, 16*m_**2) for m_ in quint))

# 1.5 Thermisch (beta = 4 pi), m = 1
th = np.where(mod.E4 <= 8)[0]
gew = np.exp(-np.pi*mod.E4[th].astype(float))
Z_modell = float(gew.sum())
Z_exakt = float(mp.nprod(lambda y: (1 + mp.e**(-4*mp.pi*(y + mp.mpf(1)/4)))*(1 + mp.e**(-4*mp.pi*(y + mp.mpf(3)/4))),
                         [0, 30])**2)
gew = gew/Z_modell
VTH = np.zeros((mod.dim, len(th)), complex)
VTH[th, np.arange(len(th))] = 1.0


def erw_th(ops, g=gew, V=VTH):
    X = kette(ops, V)
    return complex(np.sum(g*np.einsum("ij,ij->j", V.conj(), X)))


nq = 1/(math.exp(4*math.pi) - 1)
n_th = [erw_th([JA(1)]*k + [JA(-1)]*k) for k in range(1, 5)]
r_th_frei = [erw_th([JB(1)]*k + [Bf] + [JA(-1)]*k + [Bfd])/n_th[k - 1] for k in range(1, 5)]
ZAHLEN["freier Transfer r_k thermisch (m = 1)"] = [gz(x.real, 10) for x in r_th_frei]
check(f"[N] Nativer Zustand (beta = 4 pi, m = 1; {len(th)} Zustände mit E <= 2, fehlendes Gewicht "
      f"{gz(1 - Z_modell/Z_exakt, 2)}): <a^k (a^+)^k> = k! (1 + n)^k mit n = 1/(e^(4 pi) - 1) (Abweichung < 1e-9), und "
      f"r_k = {[gz(x.real, 8) for x in r_th_frei]} gegen {[gz(float(x), 8) for x in RF[1]]} im Vakuum. Die Abweichung "
      "kommt vom Niveau x = 1/4 (Besetzung 1/(1 + e^pi) = 0,0414). Qualitativ dasselbe: nicht null, nicht faktorisiert. "
      "Wegen (d^+)^4 Omega = s^4 (a^+)^4 Omega (z156) ist r_4 die Antwort der tiefsten Zustände des Sektors D = +4 "
      "(der Transfer wirkt auf der Purifikation in Standardform, B mal J B J); D = -4 ist das Tomita-Bild",
      abs(1 - Z_modell/Z_exakt) < 1e-9
      and all(abs(n_th[k - 1] - math.factorial(k)*(1 + nq)**k) < 1e-9 for k in range(1, 5))
      and all(abs(x.imag) < 1e-12 for x in r_th_frei) and abs(r_th_frei[0] - 0.5) < 1e-12
      and abs(r_th_frei[3].real - r_th_frei[0].real**4) > 0.3)

# 1.6 Kette (z152): Lücke, exponentieller Abfall
mp.mp.dps = 40
kap = mp.asinh(1)
NK = 400


def c_kette(l):
    """Gleichuntergitter-Korrelator der Kette t_1 = t_2 = t, m = 2t im Abstand 2l (trapezförmig, periodisch)."""
    return -sum(mp.cos(2*l*(-mp.pi/2 + mp.pi*j/NK))/mp.sqrt(1 + mp.cos(-mp.pi/2 + mp.pi*j/NK)**2)
                for j in range(NK))/(2*NK)


cl = {l: c_kette(l) for l in (20, 30)}
steig = (mp.log(abs(cl[20])) + mp.log(20)/2 - mp.log(abs(cl[30])) - mp.log(30)/2)/10
mp.mp.dps = 15
ZAHLEN["Kette"] = {"kappa = asinh 1": gz(kap, 12), "Abfall je Platz e^-kappa = sqrt 2 - 1": gz(mp.e**(-kap), 12),
                   "gemessene Steigung je Zelle": gz(steig, 8), "1/(2 kappa)": gz(1/(2*kap), 6)}
check(f"[N] Kette aus z152 (t_1 = t_2 = t, m = 2t): Der Fermionkorrelator fällt wie e^(-kappa r) mit kappa = asinh 1 = "
      f"ln(1 + sqrt 2) = {gz(kap, 10)} (Verzweigungspunkt der Dispersion; gemessen je Zweiplatzzelle {gz(steig, 8)} "
      f"gegen 2 kappa = {gz(2*kap, 8)}); 1/(2 kappa) = {gz(1/(2*kap), 4)} ist die Korrelationslänge aus z152. Jede "
      "verbundene Antwort ladungsneutraler Operatoren zweier getrennter Stücke fällt mindestens wie e^(-2 kappa r). "
      "Ein Potenzgesetz, also eine Linienantwort mit einem Delta, gibt es in der Kette nicht",
      abs(steig - 2*kap) < 2e-3 and abs(1/(2*kap) - mp.mpf("0.567")) < 1e-3)
befund("freier Transfer", "übertragen, aber nicht als Teilchen",
       "E: Unter dem freien Transfer geht eine Strommode in (J^A + J^B)/2 plus einen gemischten Transferstrom. Die "
       "Vierstromantwort ist nicht null; sie faktorisiert in führender Ordnung in vier Antworten 1/2 und weicht in "
       "Ordnung 1/m ab, weil die gemischten Ströme su(2)_1 bilden. Das Quintett wird nicht als Multiplett übertragen. "
       "In der Kette fällt jede Antwort exponentiell.")

# ====================================================================================================================
say("\nTeil 2: Der Austausch unter dem wechselwirkenden Gatter")
PHI = np.pi/2
QA, QB = mod.diag(mod.Q[:, 0]), mod.diag(mod.Q[:, 1])
Uy = mod.fblock(Jy0, lambda w: np.exp(-1j*PHI*w*w))
U = (mod.diag(np.exp(-1j*PHI*mod.Qt**2/4)) @ Uy).tocsr()
Ud = U.conj().T.tocsr()
HQ = (QA @ QA + QB @ QB)/2
d21 = mnorm(Bfd @ HQ @ Bf - (mod.diag(mod.Qt**2/4) + Jy0 @ Jy0))
d21b = mnorm(Bfd @ QA @ Bf - (mod.diag(mod.Qt/2) + Jy0))
kom0 = max(mnorm((U @ J0(-m_) - J0(-m_) @ U)[:, tief]) for m_ in (1, 2))
komA = float(np.linalg.norm((U @ JA(-1) - JA(-1) @ U) @ vak))
komQ = mnorm(mod.diag(mod.Qt) @ JA(-1) - JA(-1) @ mod.diag(mod.Qt))
check("[E] Gatter: B^+ Q_A B = Q/2 + J_y und B^+ H_Q B = Q^2/4 + J_y^2 mit der Transferladung J_y = J^y_0, also "
      "U(phi) = e^(-i phi Q^2/4) e^(-i phi J_y^2) (z145). Die Gesamtladung vertauscht mit allen Strommoden, J_y nur mit "
      "dem Gesamtstrom J^A + J^B: [U, J^0_(-m)] = 0, aber |[U, J^A_(-1)]|0>| = "
      f"{gz(komA, 8)}. Das Gatter hängt nur von den Uhrwerten der zirkularen Moden ab (P_null) und sieht trotzdem die "
      "Strommoden einer Naht, weil die Transferladung den relativen Strom J^z dreht. Der Sektor D = +-4 einer Naht ist "
      "unter dem Gatter nicht invariant", d21 < 1e-12 and d21b < 1e-12 and kom0 < 1e-12 and abs(komA - 1) < 1e-12
      and komQ < 1e-12)

# 2.2 Spektralzerlegung
Gam = mod.fblock(Jy0, lambda w: np.exp(1j*np.pi*w))             # Vierteldrehung A -> B -> -A (z145)
Tp = mod.fblock(Jy0, lambda w: np.exp(1j*np.pi*w/2))
Tm = Tp.conj().T.tocsr()
ger = (mod.Qt % 2 == 0)
Pg, Pu = mod.diag(ger.astype(float)), mod.diag((~ger).astype(float))
c0, c1 = (1 - 1j)/2, (1 + 1j)/2
Ein = sps.identity(mod.dim, format="csr", dtype=complex)
Uzer = mod.diag(np.exp(-1j*PHI*mod.Qt**2/4)) @ (Pg @ (c0*Ein + c1*Gam) + np.exp(-1j*np.pi/8)*Pu @ (Tp + Tm)/np.sqrt(2))
d22 = mnorm(U - Uzer)
d22b = max(mnorm(Gam @ JA(-m_) @ Gam.conj().T - JB(-m_)) for m_ in (1, 2)) + float(np.linalg.norm(Gam @ vak - vak))
check(f"[E] Auf allen {mod.dim} Zuständen: Bei gerader Gesamtladung ist U(pi/2) = e^(-i pi Q^2/8) [(1 - i)/2 + (1 + i)/2 "
      "Gamma_y] mit dem Nahttausch Gamma_y = e^(i pi J_y) (Gamma_y J^A_n Gamma_y^+ = J^B_n); bei ungerader ist "
      "U = e^(-i pi Q^2/8) e^(-i pi/8) (T_+ + T_-)/sqrt 2 mit den beiden freien Vierteltransfers T_+- = e^(+-i pi J_y/2). "
      "Das Gatter ist also eine Überlagerung von höchstens zwei freien (einkörperigen) Nahtdrehungen "
      "(z133: Wurzel des Nahttauschs; hier für alle Sektoren)", d22 < 1e-12 and d22b < 1e-12)

# 2.3 Antwort im Vakuum
u0 = complex(np.vdot(vak, U @ vak))
proben = {}
for m_, kmax in ((1, 4), (2, 3)):
    for k in range(1, kmax + 1):
        w = kette([JB(-m_)]*k, vak)
        proben[f"(a^+_{m_})^{k}"] = np.vdot(w, U @ kette([JA(-m_)]*k, vak))/(np.vdot(w, w)*u0)
for name, opsA, opsB in (("J_-2 J_-1", [JA(-2), JA(-1)], [JB(-2), JB(-1)]),
                         ("J_-3 J_-1", [JA(-3), JA(-1)], [JB(-3), JB(-1)]),
                         ("J_-4", [JA(-4)], [JB(-4)]),
                         ("J_-2 (J_-1)^2", [JA(-2), JA(-1), JA(-1)], [JB(-2), JB(-1), JB(-1)])):
    w = kette(opsB, vak)
    proben[name] = np.vdot(w, U @ kette(opsA, vak))/(np.vdot(w, w)*u0)
# ein Operator, der kein Strom ist: ein einzelnes Teilchen-Loch-Paar c^+_(y=1) c_(y=-2) der Naht
nk_, M_ = mod.nk, mod.M


def paar(ch):
    """c^+_(ch, y = 1) c_(ch, y = -2)|0> mit dem Vorzeichen der Fermionordnung (wie in Modell.T)."""
    v = np.zeros(mod.dim, complex)
    q, p = nk_*(-2 + M_) + ch, nk_*(1 + M_) + ch
    zw = ((1 << p) - 1) & ~((1 << (q + 1)) - 1)
    v[mod.idx[mod.bits[0] ^ (1 << q) ^ (1 << p)]] = -1.0 if bin(mod.bits[0] & zw).count("1") & 1 else 1.0
    return v


proben["Teilchen-Loch-Paar"] = np.vdot(paar(B_), U @ paar(A_))/u0
ZAHLEN["Gatterantwort Vakuum"] = {k_: gc(v_, 6) for k_, v_ in proben.items()}
Upi = mod.diag(np.exp(-1j*np.pi*mod.Qt**2/4)) @ mod.fblock(Jy0, lambda w: np.exp(-1j*np.pi*w*w))
w4 = kette([JB(-1)]*4, vak)
r4_pi = np.vdot(w4, Upi @ kette([JA(-1)]*4, vak))/(24*np.vdot(vak, Upi @ vak))
check("[E] Zwei-Naht-Antwort unter dem Gatter im Vakuum: <X_B^+ U X_A>/(<X^+ X> <U>) = (1 + i)/2 für (a^+)^k, k = 1..4 "
      "(m = 1), k = 1..3 (m = 2), für J_-2 J_-1, J_-3 J_-1, J_-4, J_-2 (J_-1)^2 und für ein einzelnes Teilchen-Loch-Paar, "
      f"das kein Strom ist ({len(proben)} Operatoren, Abweichung < 1e-12). Für ladungsneutrale Operatoren ohne "
      "Vakuumanteil hängt die Antwort weder von der Zahl der Quanten noch von der Mode noch vom Operator ab: Es ist die "
      "Amplitude des Tauschs der ganzen Naht. Bei "
      f"phi = pi ist sie {gz(abs(r4_pi), 6)} (voller Tausch)",
      all(abs(v_ - c1) < 1e-12 for v_ in proben.values()) and abs(abs(r4_pi) - 1) < 1e-12)
check("[E] Faktorisierung unter dem Gatter: Vier unabhängige Stromantworten gäben ((1 + i)/2)^4 = -1/4, gemessen ist "
      "(1 + i)/2. Die Antwort faktorisiert nicht, aber aus einem trivialen Grund: U = c_0 + c_1 Gamma_y ist die "
      "Überlagerung zweier freier Zweige (nichts geschieht; alles wird getauscht), Operator-Schmidt-Rang 2. Heisenberg-"
      "Form auf dem geraden Sektor: U X_A U^+ = (X_A + X_B)/2 + (i/2)(X_B - X_A) Gamma_y (geprüft für X = (J^A_-1)^4). "
      "Ein gebundener Zustand oder Pol entsteht nicht",
      abs(c1**4 + 0.25) < 1e-15 and mnorm(((U @ kette([JA(-1)]*4, Ein) @ Ud)
                                           - ((kette([JA(-1)]*4, Ein) + kette([JB(-1)]*4, Ein))/2
                                              + 0.5j*(kette([JB(-1)]*4, Ein) - kette([JA(-1)]*4, Ein)) @ Gam))
                                          [:, np.where(ger & (mod.E4 <= 8))[0]]) < 1e-12)

# 2.4 Thermisch
r_th_gat = [erw_th([JB(1)]*k + [U] + [JA(-1)]*k + [Ud])/n_th[k - 1] for k in range(1, 5)]
p_ung = float(np.sum(gew[~ger[th]]))
rein = float(mp.nprod(lambda y: (1 + mp.e**(-8*mp.pi*(y + mp.mpf(1)/4)))*(1 + mp.e**(-8*mp.pi*(y + mp.mpf(3)/4)))
                      / ((1 + mp.e**(-4*mp.pi*(y + mp.mpf(1)/4)))*(1 + mp.e**(-4*mp.pi*(y + mp.mpf(3)/4))))**2, [0, 30]))
spurGam = erw_th([Gam])                                 # = Spur rho_A^2, die Reinheit einer Naht
ZAHLEN["Gatterantwort thermisch (m = 1)"] = [gc(x, 8) for x in r_th_gat]
ZAHLEN["thermisch"] = {"P(Q ungerade)": gz(p_ung, 8), "<Gamma_y> = Reinheit einer Naht": gz(spurGam.real, 8)}
check(f"[N] Nativer Zustand (beta = 4 pi, m = 1): r_k = {[gc(x, 6) for x in r_th_gat]}. Der Imaginärteil ist bis auf "
      f"1e-5 gleich <Gamma_y>/2 = {gz(spurGam.real/2, 6)}; <Gamma_y> = Spur rho_A^2 = {gz(rein, 6)} ist die Reinheit "
      "einer Naht (der Tausch trägt nur bei gleichen Nahtzuständen bei). Der Realteil enthält die ungeraden Sektoren "
      f"(Gewicht {gz(p_ung, 5)}), in denen das Gatter aus den freien Vierteltransfers von Teil 1 besteht; die "
      "Abhängigkeit von k kommt von dort", all(abs(x.imag - spurGam.real/2) < 2e-5 for x in r_th_gat)
      and 0.05 < p_ung < 0.1 and abs(spurGam.real - rein) < 1e-9 and abs(spurGam.imag) < 1e-12
      and abs(r_th_gat[0].real - r_th_gat[3].real) > 0.01)

# 2.5 Zwei Stromarten: das Gatter tauscht gemeinsam
mod2 = Modell(4, 2, 8)                                 # Kanäle (A1, B1, A2, B2), E <= 2
jx2 = sum((mod2.T(2*f, 2*f + 1, 0) + mod2.T(2*f + 1, 2*f, 0))/2 for f in range(2))
jy2 = sum((mod2.T(2*f, 2*f + 1, 0) - mod2.T(2*f + 1, 2*f, 0))/2j for f in range(2))
U2 = mod2.diag(np.exp(-1j*PHI*mod2.Qt**2/4)) @ mod2.fblock(jy2, lambda w: np.exp(-1j*PHI*w*w))
B2 = mod2.fblock(jx2, lambda w: np.exp(-1j*np.pi/2*w))
vA = kette([mod2.T(0, 0, 1), mod2.T(2, 2, 1)], mod2.vak)
vB = kette([mod2.T(1, 1, 1), mod2.T(3, 3, 1)], mod2.vak)
g11 = np.vdot(vB, U2 @ vA)/np.vdot(mod2.vak, U2 @ mod2.vak)
f11 = np.vdot(vB, B2 @ vA)/np.vdot(mod2.vak, B2 @ mod2.vak)
check(f"[E] Zwei Stromarten ({mod2.dim} Zustände, E <= 2): Für die Quintett-Komponente a_1^+ a_2^+ gibt das Gatter "
      f"{gc(g11, 6)} = (1 + i)/2, nicht ((1 + i)/2)^2 = i/2, der freie Transfer {gz(f11.real, 6)} = r_1^2. Das Gatter "
      "hängt an der gesamten Transferladung und tauscht alle Stromarten gemeinsam. Unter dem Gatter ist die Antwort "
      "auf dem Quintett also (1 + i)/2 mal Gram-Matrix (SU(2)-invariant, aber ohne jede Auswahl), unter dem freien "
      "Transfer nach 1.4 nicht", abs(g11 - c1) < 1e-12 and abs(f11 - 0.25) < 1e-12)
befund("Gatter", "Tausch der ganzen Naht, kein Teilchen",
       "E: Das Gatter vertauscht nicht mit den Strommoden einer Naht. Seine Zwei-Naht-Antwort ist im Vakuum für jeden "
       "ladungsneutralen Operator dieselbe Zahl (1 + i)/2: Es überlagert 'nichts' und 'ganze Naht getauscht'. Das ist weder null noch ein "
       "Produkt von vier Stromantworten, aber auch keine Antwort eines Teilchens: Sie unterscheidet weder D noch Spin "
       "noch Mode, und zwischen den Nähten gibt es keine Strecke, entlang derer ein Delta messbar wäre.")

# ====================================================================================================================
say("\nTeil 3: Welche Dimension hat der Austausch?")


# Bargmann-Polynome in (d_1, d_2, e_1, e_2): Erzeuger = Multiplikation, Vernichter = Ableitung
def p_mul(p, i):
    out = {}
    for e, c in p.items():
        f = list(e)
        f[i] += 1
        out[tuple(f)] = out.get(tuple(f), 0) + c
    return out


def p_abl(p, i):
    out = {}
    for e, c in p.items():
        if e[i]:
            f = list(e)
            f[i] -= 1
            out[tuple(f)] = out.get(tuple(f), 0) + c*e[i]
    return out


def p_add(p, q, s=1):
    out = dict(p)
    for e, c in q.items():
        out[e] = out.get(e, 0) + s*c
    return {e: c for e, c in out.items() if c != 0}


def p_norm2(p):
    return sum(Fr(c)**2*math.prod(math.factorial(n) for n in e) for e, c in p.items())


HEB = {"L_1": [((0, 3), 1), ((1, 2), 1)],              # d_1^+ e_2^+ + d_2^+ e_1^+
       "L_2": [((0, 3), 1), ((1, 2), -1)],             # (bis auf i) d_1^+ e_2^+ - d_2^+ e_1^+
       "L_3": [((0, 2), 1), ((1, 3), -1)]}             # d_1^+ e_1^+ - d_2^+ e_2^+


def heben(p, name):
    out = {}
    for (i, j), s in HEB[name]:
        out = p_add(out, p_mul(p_mul(p, i), j), s)
    return out


def senken(p, name):
    out = {}
    for (i, j), s in HEB[name]:
        out = p_add(out, p_abl(p_abl(p, i), j), s)
    return out


def energie(p):
    return {e: c*(Fr(sum(e), 2) + 1) for e, c in p.items()}


# 3.1 Rahmen der Naht (Möbius-sl(2) von L_0, L_+-1): :J^4: hat Gewicht 4
def part4(n, k=4, kleinst=1):
    if k == 0:
        return [()] if n == 0 else []
    return [(a,) + r for a in range(kleinst, n // k + 1) for r in part4(n - a, k - 1, a)]


lin2d = {}
for n in range(4, 17):
    tot = Fr(0)
    for lam in part4(n):
        mult = math.prod(math.factorial(lam.count(x)) for x in set(lam))
        tot += Fr(24, mult)**2*math.prod(lam)*mult
    lin2d[n] = tot
faltung = {n: sum(a*b*c*(n - a - b - c) for a in range(1, n) for b in range(1, n - a) for c in range(1, n - a - b))
           for n in range(4, 17)}
komp_ok = all(sum(binom(n1 + (4 - j) - 1, 2*(4 - j) - 1) * binom(n - n1 + j - 1, 2*j - 1)
                  for n1 in range(4 - j, n - j + 1)) == binom(n + 3, 7) for j in range(1, 4) for n in range(4, 17))
check("[E] Rahmen der Naht: Der Vierstromoperator :J^4: ist quasiprimär vom Gewicht 4 unter der Möbius-sl(2) der Naht "
      "(L_1 (J_-1)^4|0> = 0). Die Norm seiner Mode n ist 24 C(n + 3, 7) (n = 4..16, über Partitionen), also die "
      "Linienantwort von Delta = 4, und sie ist exakt die vierfache Faltung der Stromantwort n (Delta = 1): 1/x^8 = "
      "(1/x^2)^4. (a_m^+)^4 ist darin der eine Term der Mode n = 4m mit Gewicht 24 m^4. Dasselbe gilt für alle fünf "
      "Komponenten :J_1^(4-j) J_2^j:. 'Delta = 4' und 'vier Stromquanten' sind im Rahmen der Naht dieselbe Aussage",
      all(lin2d[n] == 24*binom(n + 3, 7) for n in lin2d) and all(faltung[n] == binom(n + 3, 7) for n in faltung)
      and komp_ok)

# 3.2 Rahmen der so(3,2): tiefste Gewichte des Sektors D = +4 und die sl(2) aus E und L^+-_k
alg_ok, lin32 = True, {}
for j in range(5):
    v0 = {(4 - j, j, 0, 0): 1}
    for name in sorted(HEB):
        p = dict(v0)
        alg_ok &= senken(p, name) == {}
        for n in range(0, 7):
            if n:
                p = heben(p, name)
            lin32[(j, name, n)] = p_norm2(p)/(Fr(math.factorial(n))**2*p_norm2(v0))
        q = heben(v0, name)                             # [L^-, L^+] = 2 E auf einem Zustand prüfen
        alg_ok &= p_add(senken(heben(q, name), name), heben(senken(q, name), name), -1) == \
            {e: 2*c for e, c in energie(q).items()}
check("[E] Rahmen der so(3,2): Für jede der drei nativen Hebungen L^+_k = d^+ sigma_k e^+ gilt [L^-_k, L^+_k] = 2E; mit E "
      "bilden sie eine sl(2). Die fünf tiefsten Zustände (d_1^+)^(4-j) (d_2^+)^j Omega = s^4 (a_1^+)^(4-j) (a_2^+)^j Omega "
      "werden von L^- vernichtet, und |(L^+_k)^n v|^2/(n!^2 |v|^2) = C(n + 5, 5) = C(N + 2, 5) mit N = 3 + n "
      "(15 Fälle, n = 0..6, exakt): die Linienantwort von Delta = 3. Das folgt allein aus E_0 = 3 und ist in D(3,2) "
      "von selbst erfüllt", alg_ok and all(v_ == binom(n + 5, 5) for (j, name, n), v_ in lin32.items()))

# 3.3 Vergleich
tab = {N: (binom(N + 2, 5), binom(N + 3, 7), binom(N + 5, 11)) for N in range(3, 10)}
ZAHLEN["Linienantworten Delta = 3 | 4 | 6 für n = 3..9"] = {str(N): list(v_) for N, v_ in tab.items()}
x_, y_ = sp.symbols("x y")
P_ = 1/(x_ - y_)
cO = [(sp.Rational(1, 2), 0), (-2, 1), (sp.Rational(1, 2), 2)]          # Koeffizient, Ableitungen an phi_alpha
zwei = sum(c*d*sp.diff(P_, x_, k, y_, l)*sp.diff(P_, x_, 2 - k, y_, 2 - l) for c, k in cO for d, l in cO)
zwei_konst = sp.simplify(zwei*(x_ - y_)**6)
null0 = sp.simplify(sum(c*sp.diff(P_, x_, k)*sp.diff(P_, x_, 2 - k) for c, k in cO))
check("[E] Die beiden Antworten widersprechen sich nicht, sie gehören zu zwei verschiedenen sl(2). L_(-1) der Naht "
      "verschiebt die Mode m, L^+ der so(3,2) erzeugt bei fester Mode ein Paar aus Original und Partner. Ein Stromquant "
      "trägt L_0 = m, aber so(3,2)-Energie 1/2: E_0 = 3 = 1 (Nullpunkt = zwei Singletons zu je 1/2) + 4 mal 1/2. Für "
      f"n = 3..6: Delta = 3 gibt {[tab[N][0] for N in range(3, 7)]}, Delta = 4 gibt {[tab[N][1] for N in range(3, 7)]}. "
      "Vier Quanten als vier so(3,2)-Teilchen D(3/2,1/2) gäben Delta = 6, was auch nicht vorkommt",
      [tab[N][0] for N in range(3, 7)] == [1, 6, 21, 56] and [tab[N][1] for N in range(3, 7)] == [0, 1, 8, 36]
      and energie({(4, 0, 0, 0): 1}) == {(4, 0, 0, 0): 3} and energie({(1, 0, 0, 0): 1}) == {(1, 0, 0, 0): Fr(3, 2)}
      and energie({(0, 0, 0, 0): 1}) == {(0, 0, 0, 0): 1} and [tab[N][2] for N in range(3, 8)] == [0, 0, 0, 1, 12])
check("[E] Flato-Fronsdal auf der Linie nachgerechnet: Zwei Felder vom Gewicht 1/2 (die sl(2)-Gewichte der beiden "
      "Singleton-Vakua; Propagator 1/(x - y)) haben auf Stufe 2 genau ein Quasiprimärfeld "
      "O = phi_a phi_b''/2 - 2 phi_a' phi_b' + phi_a'' phi_b/2. Seine Zweipunktfunktion ist nach Wick ein Produkt zweier "
      f"Singleton-Propagatoren und zugleich exakt {zwei_konst}/(x - y)^6, also die Antwort von Delta = 3; die Mischung "
      "mit phi_a phi_b verschwindet. Faktorisierung nach Wick und 'Antwort eines Teilchens mit Delta = 3' schließen "
      "sich nicht aus. So ist es in der Dualität zwischen freien Randfeldern und Höherer-Spin-Gravitation [L]",
      zwei_konst.is_number and zwei_konst != 0 and null0 == 0)
befund("Kriterium aus z156", "falsch gestellt",
       "E: 'Faktorisiert in vier Stromquanten (Delta = 4)' und 'D(3,2) mit Delta = 3' gelten zugleich, das eine im "
       "Rahmen der Naht (Möbius), das andere im Rahmen der so(3,2) (Original + Partner). Aus Sicht der so(3,2) sind die "
       "vier Quanten kein Vierteilchenzustand, sondern die Anregung eines Zwei-Singleton-Zustands (Flato-Fronsdal). "
       "Faktorisierung widerlegt P_X nicht, und die Linienantwort C(n + 2, 5) bestätigt P_X nicht (sie ist "
       "Darstellungstheorie). Was P_X testet, ist, ob der Singleton ein Feld mit invariantem Vakuum und Arten ist.")

# ====================================================================================================================
say("\nTeil 4: Singleton-Lesart")


# Charaktere als Gewichte {(2E, 2 J_3): Vielfachheit} und als Spins {2E: {2s: Vielfachheit}}
def cg2(a, b):
    return range(abs(a - b), a + b + 1, 2)


def sym1(n):
    return [2*l for l in range(n % 2, n + 1, 2)]


def ch_add(A, B, c=1):
    out = {e: dict(v) for e, v in A.items()}
    for e, v in B.items():
        tgt = out.setdefault(e, {})
        for s, m in v.items():
            tgt[s] = tgt.get(s, 0) + c*m
    return {e: {s: m for s, m in sorted(v.items()) if m} for e, v in sorted(out.items()) if any(v.values())}


def ch_mul(A, B, emax2):
    out = {}
    for e1, c1 in A.items():
        for e2, c2 in B.items():
            if e1 + e2 > emax2:
                continue
            tgt = out.setdefault(e1 + e2, {})
            for s1, m1 in c1.items():
                for s2, m2 in c2.items():
                    for s in cg2(s1, s2):
                        tgt[s] = tgt.get(s, 0) + m1*m2
    return ch_add(out, {})


def verma(E2, s2, emax2):
    out, n = {}, 0
    while E2 + 2*n <= emax2:
        tgt = out.setdefault(E2 + 2*n, {})
        for l2 in sym1(n):
            for j2 in cg2(s2, l2):
                tgt[j2] = tgt.get(j2, 0) + 1
        n += 1
    return out


def irrep(E2, s2, emax2):
    V = verma(E2, s2, emax2)
    if s2 >= 2 and E2 == s2 + 2:
        return ch_add(V, verma(E2 + 2, s2 - 2, emax2), -1)
    if (E2, s2) == (1, 0):
        return ch_add(V, verma(5, 0, emax2), -1)
    if (E2, s2) == (2, 1):
        return ch_add(V, verma(4, 1, emax2), -1)
    return V


def zerlege(chi, emax2):
    chi = ch_add(chi, {})
    out = []
    for e in range(0, emax2 + 1):
        while chi.get(e):
            s = max(chi[e])
            m = chi[e][s]
            chi = ch_add(chi, irrep(e, s, emax2), -m)
            if any(v < 0 for lev in chi.values() for v in lev.values()):
                return None
            out.append((e, s, m))
    return out


def gew_singleton(art, emax2):
    """Rac (gerade Quantenzahl n) bzw. Di (ungerade): E = 1/2 + n/2, Spin n/2."""
    return {(1 + n, w): 1 for n in range(0 if art == "Rac" else 1, emax2, 2) for w in range(-n, n + 1, 2)
            if 1 + n <= emax2}


def gw_mul(A, B, emax2):
    out = {}
    for (e1, w1), m1 in A.items():
        for (e2, w2), m2 in B.items():
            if e1 + e2 <= emax2:
                out[(e1 + e2, w1 + w2)] = out.get((e1 + e2, w1 + w2), 0) + m1*m2
    return out


def gw_adams(A, emax2):
    return {(2*e, 2*w): m for (e, w), m in A.items() if 2*e <= emax2}


def gw_komb(A, B, s):
    out = dict(A)
    for k_, m in B.items():
        out[k_] = out.get(k_, 0) + s*m
    return {k_: m // 2 for k_, m in out.items() if m}


def gw_spins(A):
    out = {}
    for (e, w), m in sorted(A.items()):
        if w >= 0:
            d = m - A.get((e, w + 2), 0)
            if d:
                out.setdefault(e, {})[w] = d
    return out


def hb(x2):
    return str(x2 // 2) if x2 % 2 == 0 else f"{x2}/2"


def text(zer, bis2):
    return [f"D({hb(e)},{hb(s)})" + (f" x{m}" if m != 1 else "") for e, s, m in zer if e <= bis2]


EM = 22
rac, di = gew_singleton("Rac", EM), gew_singleton("Di", EM)
zer = {}
for name, A in (("Rac", rac), ("Di", di)):
    quad = gw_mul(A, A, EM)
    zer[f"{name} x {name}"] = zerlege(gw_spins(quad), EM)
    zer[f"Sym^2 {name}"] = zerlege(gw_spins(gw_komb(quad, gw_adams(A, EM), 1)), EM)
    zer[f"Lambda^2 {name}"] = zerlege(gw_spins(gw_komb(quad, gw_adams(A, EM), -1)), EM)
zer["Rac x Di"] = zerlege(gw_spins(gw_mul(rac, di, EM)), EM)
BIS = 18                                               # vollständig bis E = 9
ZAHLEN["Zerlegungen bis E = 9"] = {k_: text(v_, BIS) for k_, v_ in sorted(zer.items())}
erw_rr = [(2*s + 2, 2*s, 1) for s in range(0, 9)]
erw_dd = sorted([(4, 0, 1)] + [(2*s + 2, 2*s, 1) for s in range(1, 9)])
kurz = lambda z: sorted(t for t in z if t[0] <= BIS)
anz32 = {k_: sum(m for e, s, m in v_ if (e, s) == (6, 4)) for k_, v_ in zer.items()}
check("[E] Flato-Fronsdal über Charaktere (bis E = 9): Rac x Rac = Summe_(s >= 0) D(s+1, s); Di x Di = D(2,0) + "
      "Summe_(s >= 1) D(s+1, s); Rac x Di = Summe halbzahliger s. Symmetrisch und antisymmetrisch getrennt: "
      f"Sym^2 Rac = {text(zer['Sym^2 Rac'], 12)}..., Lambda^2 Rac = {text(zer['Lambda^2 Rac'], 12)}..., "
      f"Lambda^2 Di = {text(zer['Lambda^2 Di'], 12)}..., Sym^2 Di = {text(zer['Sym^2 Di'], 12)}.... D(3,2) liegt je "
      "einmal in Sym^2 Rac (Bosonen) und in Lambda^2 Di (Fermionen), nie im Teil mit der falschen Statistik",
      kurz(zer["Rac x Rac"]) == erw_rr and kurz(zer["Di x Di"]) == erw_dd
      and kurz(zer["Sym^2 Rac"]) == [t for t in erw_rr if (t[1]//2) % 2 == 0]
      and kurz(zer["Lambda^2 Rac"]) == [t for t in erw_rr if (t[1]//2) % 2 == 1]
      and kurz(zer["Lambda^2 Di"]) == sorted([(4, 0, 1)] + [t for t in erw_rr if t[1] >= 4 and (t[1]//2) % 2 == 0])
      and kurz(zer["Sym^2 Di"]) == [t for t in erw_rr if (t[1]//2) % 2 == 1]
      and all(t[1] % 2 == 1 for t in kurz(zer["Rac x Di"]))
      and (anz32["Sym^2 Rac"], anz32["Lambda^2 Rac"], anz32["Lambda^2 Di"], anz32["Sym^2 Di"]) == (1, 0, 1, 0))


# 4a: Singuletts der Artengruppe O(N) in C^N x C^N
def on_invarianten(N):
    I = np.eye(N)
    zeilen = []
    for a in range(N):
        for b in range(a + 1, N):
            G = np.zeros((N, N))
            G[a, b], G[b, a] = 1, -1
            zeilen.append(np.kron(G, I) + np.kron(I, G))
    R = np.eye(N)
    R[0, 0] = -1
    zeilen.append(np.kron(R, R) - np.eye(N*N))
    S = np.zeros((N*N, N*N))
    for a in range(N):
        for b in range(N):
            S[a*N + b, b*N + a] = 1
    out = []
    for s in (1, -1):
        Mx = np.vstack(zeilen + [S - s*np.eye(N*N)])
        out.append(N*N - int(np.linalg.matrix_rank(Mx, tol=1e-9)))
    return tuple(out)


inv = {N: on_invarianten(N) for N in range(1, 8)}
ZAHLEN["O(N)-Invarianten in (Sym^2, Lambda^2) von C^N"] = {str(N): list(v_) for N, v_ in inv.items()}
check("[E] Zweite Quantisierung mit N Arten (Einteilchenraum C^N x Rac, Artengruppe O(N)): Der Zweiteilchenraum ist "
      "Sym^2 C^N x Sym^2 Rac + Lambda^2 C^N x Lambda^2 Rac. O(N)-Invarianten in (Sym^2 C^N, Lambda^2 C^N) für N = 1..7: "
      f"{[inv[N] for N in sorted(inv)]}. Der Singulettsektor ist also für jedes N genau Sym^2 Rac = Summe gerader s "
      "(für Majorana-Arten Lambda^2 Di) und enthält genau ein D(3,2), unabhängig von N. Für U(N) ist es Rac x Rac, "
      "ebenfalls mit genau einem D(3,2). In dieser Lesart ersetzte die Singulettbedingung P_mode, P_paar und P_chir",
      all(v_ == (1, 0) for v_ in inv.values()) and anz32["Rac x Rac"] == 1 and anz32["Di x Di"] == 1)

# 4a': Wo liegen die Sektoren D = +-4 des Dubletons?
al1, al2, be1, be2 = sp.symbols("alpha1 alpha2 beta1 beta2")
lage_ok = True
for j in range(5):
    pp = sp.expand(((al1 - sp.I*be1)/sp.sqrt(2))**(4 - j)*((al2 - sp.I*be2)/sp.sqrt(2))**j)      # D = +4
    pm = sp.expand(((al1 + sp.I*be1)/sp.sqrt(2))**(4 - j)*((al2 + sp.I*be2)/sp.sqrt(2))**j)      # D = -4
    tau = lambda p: sp.expand(p.subs({al1: be1, al2: be2, be1: al1, be2: al2}, simultaneous=True))
    sig = lambda p: sp.expand(p.subs({be1: -be1, be2: -be2}))
    lage_ok &= sp.expand(tau(pp) - pm) == 0 and sp.expand(sig(pp) - pm) == 0
    for p, par in ((pp + pm, 0), (pp - pm, 1)):
        po = sp.Poly(p, al1, al2, be1, be2)
        lage_ok &= all((e[2] + e[3]) % 2 == par and (e[0] + e[1]) % 2 == par for e in po.monoms())
check("[E] Lage im Dubleton (alpha, beta die beiden Singleton-Oszillatorpaare, d = (alpha + i beta)/sqrt 2): Der Tausch "
      "tau: alpha <-> beta und die Spiegelung beta -> -beta bilden D = +4 auf D = -4 ab. Die tau-gerade Kombination "
      "liegt ganz in Rac x Rac (symmetrisch, also Sym^2 Rac), die tau-ungerade ganz in Di x Di (antisymmetrisch, also "
      "Lambda^2 Di). Die Statistik wählt also in jedem Fall eine Diagonale, nie eine chirale Kopie; nach z156 dreht "
      "der Modulfluss e^(itK) (Phasen e^(+-16 pi i m t)) die Diagonalen ineinander. Die Singulettbedingung ist mit der "
      "Modulinvarianz der einzelnen Kopie nicht verträglich, P_chir würde durch eine Parität ersetzt", lage_ok)

# Tomita-Konjugation: J d_i J = e_i (antilinear); alpha_i = (d_i + e'_i)/sqrt 2, beta_i = -i (d_i - e'_i)/sqrt 2,
# e'_1 = e_2, e'_2 = -e_1 (z156). Operatoren als Koeffizientenvektoren in (d_1, d_2, e_1, e_2).
w2 = sp.sqrt(2)
alpha_v = [sp.Matrix([1, 0, 0, 1])/w2, sp.Matrix([0, 1, -1, 0])/w2]
beta_v = [-sp.I*sp.Matrix([1, 0, 0, -1])/w2, -sp.I*sp.Matrix([0, 1, 1, 0])/w2]
Jkonj = lambda v: sp.Matrix([sp.conjugate(v[2]), sp.conjugate(v[3]), sp.conjugate(v[0]), sp.conjugate(v[1])])
check("[E] Die Tomita-Konjugation tauscht die beiden Singletons: J alpha_1 J = i beta_2, J alpha_2 J = -i beta_1 (und "
      "zurück). Ein einzelner Singleton ist also weder Original noch Partner noch J-invariant, sondern eine Hälfte, die "
      "beide mischt; J-invariant ist erst das Paar. Ein Feld, dessen Einteilchenraum ein Singleton ist, liegt damit "
      "nicht in der Algebra einer Naht",
      Jkonj(alpha_v[0]) == sp.I*beta_v[1] and Jkonj(alpha_v[1]) == -sp.I*beta_v[0]
      and all(sp.simplify(Jkonj(Jkonj(v)) - v) == sp.zeros(4, 1) for v in alpha_v + beta_v))

# 4b: native Gruppen
# (i) Drehung eines Strompaars = Erzeuger der so(3,2)
drehung_ok = True
for e in itertools.product(range(3), repeat=4):
    p = {e: 1}
    lhs = p_add(senken(heben(p, "L_3"), "L_1"), heben(senken(p, "L_1"), "L_3"), -1)       # [L^-_1, L^+_3] p
    G = {}
    for (i, j), s in (((0, 1), 1), ((1, 0), -1), ((2, 3), 1), ((3, 2), -1)):             # x_i^+ x_j
        G = p_add(G, p_mul(p_abl(p, j), i), s)
    drehung_ok &= lhs == G
check("[E] Eine reelle Drehung zweier Cartan-Ströme (a_1, a_2) -> R (a_1, a_2), mit ihr (b_1, b_2), ist auf dem Dubleton "
      "der Erzeuger (d_1^+ d_2 - d_2^+ d_1) + (e_1^+ e_2 - e_2^+ e_1) = [L^-_1, L^+_3] (81 Monome), also eine Drehung der "
      "nativen so(3,2). Jede Gruppe, die Cartan-Ströme mischt (O(8), Weyl-Gruppen), wirkt auf ein Strompaar durch die "
      "Drehgruppe und nicht als Artengruppe; ihr 'Singulett' im Quintett wäre eine einzelne Polarisation, kein "
      "so(3,2)-Untermodul. Eine Artengruppe muss mit der so(3,2) vertauschen", drehung_ok)


# (ii) Howe: Kopien unter gemeinsamer so(3,2); O(N_c)-Singulett
def monome(nvar, grad):
    return [e for e in itertools.product(range(grad + 1), repeat=nvar) if sum(e) == grad]


def invarianten_dim(nc, grad):
    """Dimension der so(nc)-Invarianten in den Polynomen vom Grad grad in x^c_i (c < nc, i = 1, 2)."""
    mon = monome(2*nc, grad)
    ix = {e: k for k, e in enumerate(mon)}
    S = np.zeros((len(mon), len(mon)))
    for c in range(nc):
        for d in range(c + 1, nc):
            rows, cols, vals = [], [], []
            for k, e in enumerate(mon):
                for i in range(2):
                    for (von, nach, s) in ((2*d + i, 2*c + i, 1), (2*c + i, 2*d + i, -1)):   # x^nach d/dx^von
                        if e[von]:
                            f = list(e)
                            f[von] -= 1
                            f[nach] += 1
                            rows.append(ix[tuple(f)])
                            cols.append(k)
                            vals.append(s*e[von])
            Rm = sps.csr_matrix((vals, (rows, cols)), shape=(len(mon), len(mon)))
            S += (Rm.T @ Rm).toarray()
    w = np.linalg.eigvalsh(S)
    return int(np.sum(w < 1e-8*max(1.0, w[-1])))


howe2 = [invarianten_dim(2, 2*k) for k in range(5)]
# O(2) = SO(2) mit Spiegelung tau: gerade Invarianten zählen (D = 0 und tau = +1)


def tau_gerade_dim(grad):
    mon = monome(4, grad)                               # (alpha1, alpha2, beta1, beta2)
    ix = {e: k for k, e in enumerate(mon)}
    Dm = np.zeros((len(mon), len(mon)))
    for k, e in enumerate(mon):
        for i in range(2):
            for (von, nach, s) in ((2 + i, i, 1), (i, 2 + i, -1)):
                if e[von]:
                    f = list(e)
                    f[von] -= 1
                    f[nach] += 1
                    Dm[ix[tuple(f)], k] += s*e[von]
    Tm_ = np.zeros((len(mon), len(mon)))
    for k, e in enumerate(mon):
        Tm_[ix[(e[2], e[3], e[0], e[1])], k] = 1
    Mx = np.vstack([Dm, Tm_ - np.eye(len(mon))])
    return len(mon) - int(np.linalg.matrix_rank(Mx, tol=1e-9))


# Achtung: in tau_gerade_dim sind die Variablen (alpha_i, beta_i), in invarianten_dim(2, .) (x^1_i, x^2_i) = dasselbe.
howe2g = [tau_gerade_dim(2*k) for k in range(5)]
howe4 = [invarianten_dim(4, 2*k) for k in range(4)]
EC = 12
einz = {1 + n: {n: 1} for n in range(0, EC)}            # Di + Rac
mehr = {}
chiN = {0: {0: 1}}
for Nd in (1, 2, 3):
    chiN = ch_mul(ch_mul(chiN, einz, EC), einz, EC)
    z_ = zerlege(chiN, EC)
    s2 = sorted(t for t in z_ if t[1] == 4)
    mehr[Nd] = (min(chiN), sum(m for e, s, m in z_ if (e, s) == (6, 4)), s2[0])
ZAHLEN["Kopien unter gemeinsamer so(3,2)"] = {
    "SO(2)-Invarianten je Stufe (2 Kopien)": howe2, "O(2)-Invarianten": howe2g, "O(4)-Invarianten (4 Kopien)": howe4,
    "Dubletons 1, 2, 3: (2 E_min, Vielfachheit D(3,2), tiefstes Spin-2-Gewicht)": {str(k_): str(v_) for k_, v_ in mehr.items()}}
check("[E] Kopien als Vielfachheit einer gemeinsamen so(3,2) (nur so kann eine Gruppe auf Kopien wirken und mit ihr "
      "vertauschen): Die Kommutante ist die Howe-duale O(N_c) der N_c Singleton-Kopien. Ihr Singulettsektor hat je "
      f"Stufe die Dimension {howe2g} (N_c = 2) und {howe4} (N_c = 4) = Zahl der Monome in den drei Hebungen, ist also "
      "der Skalar D(N_c/2, 0): D(1,0) für einen Dubleton, D(2,0) für zwei. D(3,2) ist nie ein Singulett, sondern das "
      f"O(2)-Dublett mit Ladung +-4. Für 1, 2, 3 Dubletons ist die Vielfachheit von D(3,2) {[mehr[k_][1] for k_ in (1, 2, 3)]} "
      f"und die tiefste Energie {[hb(mehr[k_][0]) for k_ in (1, 2, 3)]} (z156 reproduziert). Das Singulett einer Gruppe "
      "auf den Dubleton-Kopien kann P_mode und P_paar also nicht ersetzen",
      howe2 == [1, 4, 9, 16, 25] and howe2g == [1, 3, 6, 10, 15] and howe4 == [1, 3, 6, 10]
      and [mehr[k_][1] for k_ in (1, 2, 3)] == [2, 0, 0] and [mehr[k_][0] for k_ in (1, 2, 3)] == [2, 4, 6]
      and [mehr[k_][2] for k_ in (1, 2, 3)] == [(6, 4, 2), (8, 4, 25), (10, 4, 105)])

# (iii) Welche Gruppen auf den acht Cartan-Strömen sind nativ?
wurzeln = [tuple(v) for v in itertools.product((-1, 0, 1), repeat=8) if sum(abs(x) for x in v) == 2]
wurzeln = [tuple(2*x for x in v) for v in wurzeln]      # in halben Einheiten: Einträge +-2
wurzeln += [v for v in itertools.product((-1, 1), repeat=8) if sum(1 for x in v if x < 0) % 2 == 0]
wurzeln = sorted(wurzeln)
verdr = (1,)*8                                          # Viertelverdrillung (1/4)^8, bis auf den Faktor
neutral = [w for w in wurzeln if sum(a*b for a, b in zip(w, verdr)) == 0]
perm = [w for w in neutral if sorted(abs(x) for x in w) == [0]*6 + [2, 2]]
W_E8 = 696729600
stab = W_E8 // len(wurzeln)
rng = np.random.default_rng(158)
PKT = rng.integers(-4, 5, size=(420, 8)).astype(float)/4
MON4 = np.array(monome(8, 4), dtype=float)


def mono_werte(X):
    return np.prod(X[:, None, :]**MON4[None, :, :], axis=2)


def quartische_invarianten(spiegel):
    S = np.zeros((len(MON4), len(MON4)))
    M0 = mono_werte(PKT)
    for w in spiegel:
        a = np.array(w, float)/2
        Dm = mono_werte(PKT - np.outer(PKT @ a, a)*(2/(a @ a))) - M0
        S += Dm.T @ Dm
    ev = np.linalg.eigvalsh(S)
    return int(np.sum(ev < 1e-9*ev[-1]))


halb = lambda ws: [w for w in ws if next(x for x in w if x) > 0]
qinv = {"W(E8)": quartische_invarianten(halb(wurzeln)), "W(E7)": quartische_invarianten(halb(neutral)),
        "S_8": quartische_invarianten(halb(perm))}
ZAHLEN["Gruppen auf den acht Cartan-Strömen"] = {
    "E8-Wurzeln": len(wurzeln), "davon neutral zur Viertelverdrillung": len(neutral),
    "Ordnung Stabilisator W(E7)": stab, "Ordnung S_8": math.factorial(8), "Index": stab // math.factorial(8),
    "quartische Invarianten": qinv}
check(f"[E/L] Gruppen im Bestand: Von den {len(wurzeln)} E8-Wurzeln sind {len(neutral)} neutral zur Viertelverdrillung "
      f"(E7, Bestand H6h20 und H6h21); der Stabilisator in W(E8) hat die Ordnung {W_E8}/{len(wurzeln)} = {stab} = |W(E7)| und "
      f"enthält S_8 (Vertauschung der acht Ströme, Index {stab // math.factorial(8)}). W(E8) erhält die Algebra, aber "
      "nicht die Verdrillung; O(8) ist nur eine Symmetrie des Oszillatorzustands einer Mode (wie SU(2)_L in z75). "
      f"Quartische Invarianten in den acht Strömen: {qinv['W(E8)']} (W(E8): nur (H.H)^2, das Quadrat der "
      f"Sugawara-Form), {qinv['W(E7)']} (W(E7)), {qinv['S_8']} (S_8). Keine dieser Gruppen vertauscht mit der so(3,2) "
      "eines Strompaars (4b (i)); eine Paarung der acht Ströme in vier Paare zeichnet der Bestand nicht aus (genannt "
      "sind die Teilungen 5 + 3 und 3 + 2), und zwischen den Moden m vermittelt keine zustandserhaltende Gruppe (z156: "
      "die so(3,2) der Moden unterscheiden sich um coth(pi m)). Jede Artengruppe wäre neu",
      len(wurzeln) == 240 and len(neutral) == 126 and len(perm) == 56 and stab == 2903040
      and qinv == {"W(E8)": 1, "W(E7)": 3, "S_8": 5})

# 4c: Welche Lesart verträgt die native Wirkung?
chi2 = ch_mul(einz, einz, EC)
z2 = zerlege(chi2, EC)
check("[E] Zweite Quantisierung: (1) Der Dubleton-Fockraum hat die tiefste Energie "
      f"{hb(min(chi2))} und keinen so(3,2)-invarianten Vektor (z156). Ein Vielteilchenraum über einem Einteilchenraum V "
      "hat ein invariantes Vakuum bei E = 0 und additive Energien; der Dubleton-Fockraum ist das für kein V. (2) Auf "
      "jedem Sektor ist die Zahl der Stromquanten N_a = E_ab - 1 + D/2 (z156), also die Energie innerhalb einer "
      "irreduziblen Darstellung und keine Teilchenzahl; die Hebungen L^+ = d^+ sigma e^+ ändern sie. (3) Mehrere "
      "Dubletons unter gemeinsamer so(3,2) haben die Nullpunktsenergie N statt 0 (4b). Verträglich mit der "
      "metaplektischen Wirkung ist nur: Der Fockraum ist der Raum von genau zwei Singletons (alpha und beta, je genau "
      "eines) und damit der Einteilchenraum der masselosen Teilchen aller Helizitäten (Flato-Fronsdal). Die a_i sind "
      "Feldmoden, aber eines Feldes der Naht: zweitquantisiert bezüglich der Möbius-sl(2) (Einteilchencharakter "
      "q/(1 - q), Vakuum invariant), nicht bezüglich der so(3,2). P_1T ist also nicht die falsche Lesart, sondern die "
      "einzige; was fehlt, ist die zweite Quantisierung des Singletons (neu: P_2Q)",
      min(chi2) == 2 and 0 not in chi2 and sum(m for e, s, m in z2 if (e, s) == (6, 4)) == 2
      and (2, 0, 1) in z2 and (4, 0, 1) in z2)

# 4d: Zählung für G [L]
mp.mp.dps = 30
CT1 = mp.mpf(3)/(32*mp.pi**2)                           # je reellem Skalar oder Majorana-Fermion in d = 3 [L]
S3 = 4*mp.pi                                           # Fläche der Einheitssphäre in d = 3
ct_skalar = 3/(2*S3**2)                                # Osborn-Petkou: d/((d - 1) S_d^2)
ct_dirac = 3*2/(2*S3**2)                               # d/(2 S_d^2) mal Spinordimension 2
d_ = 3
vorf = mp.mpf(d_ + 1)/(d_ - 1)*mp.gamma(d_ + 1)/(mp.pi**(mp.mpf(d_)/2)*mp.gamma(mp.mpf(d_)/2))/(8*mp.pi)   # C_T G/l^2
l2G_je_art = CT1/vorf                                  # = pi/32
ainv = mp.mpf("137.0359992168")                        # EM.FP.01 (gelöster Block, nur benutzt)
l2G_ledger = 32*mp.pi**3*mp.e**(2*ainv)
N_noetig = l2G_ledger/l2G_je_art
s0 = mp.log(2)/2
naehte = 32*mp.pi**4*mp.e**(2*ainv)/s0
nativ = {"2 (ein Dubleton: alpha, beta)": 2, "8 (Cartan-Ströme)": 8, "16 (Majorana der Naht)": 16,
         "56 (28 Strompaare)": 56, "248 (E8-Ströme)": 248}
ZAHLEN["G als Anzahl [L]"] = {
    "C_T je Art": mp.nstr(CT1, 12), "l^2/G je Art = pi/32": mp.nstr(l2G_je_art, 12),
    "nativ l^2/G": {k_: mp.nstr(v_*l2G_je_art, 6) for k_, v_ in nativ.items()},
    "Ledger l^2/G": mp.nstr(l2G_ledger, 6), "nötige Artenzahl": mp.nstr(N_noetig, 6),
    "Nähte je Horizont (z136)": mp.nstr(naehte, 6), "Verhältnis = 16 ln 2/pi^2": mp.nstr(N_noetig/naehte, 10)}
check("[L] G als Anzahl: Mit C_T = 3/(32 pi^2) je reellem Skalar oder Majorana-Fermion (Osborn-Petkou: Skalar "
      "d/((d-1) S_d^2), Dirac d/(2 S_d^2) mal 2; S_3 = 4 pi) und C_T = 3 l^2/(pi^3 G) (allgemeine Formel bei d = 3) "
      f"ist l^2/G = pi^3 C_T/3 = pi N/32 = {mp.nstr(l2G_je_art, 6)} N. Native Zählungen geben "
      f"{[mp.nstr(v_*l2G_je_art, 4) for v_ in nativ.values()]} für N = {list(nativ.values())}: planckgroß wie alle "
      f"früheren Wege. Der Ledgerwert 32 pi^3 e^(2/alpha) = {mp.nstr(l2G_ledger, 5)} verlangt N = 1024 pi^2 e^(2/alpha) "
      f"= {mp.nstr(N_noetig, 5)} Arten. Das ist das {mp.nstr(N_noetig/naehte, 6)}-Fache (= 16 ln 2/pi^2) der "
      f"{mp.nstr(naehte, 3)} Nähte je Horizont aus z136. Beide Zahlen sind dasselbe l^2/G in zwei Einheiten "
      "(pi/32 je Art gegen ln 2/2 je Naht); das Verhältnis ist keine ganze Zahl und kein Treffer",
      abs(ct_skalar - CT1) < 1e-25 and abs(ct_dirac/2 - CT1) < 1e-25 and abs(vorf - 3/mp.pi**3) < 1e-25
      and abs(l2G_je_art - mp.pi/32) < 1e-25 and abs(N_noetig/naehte - 16*mp.log(2)/mp.pi**2) < 1e-20
      and abs(naehte/mp.mpf("9.6e122") - 1) < 0.01 and abs(N_noetig/(1024*mp.pi**2*mp.e**(2*ainv)) - 1) < 1e-20)
mp.mp.dps = 15
befund("Singleton-Lesart", "bedingt (P_2Q, P_sing, P_N)",
       "E/L: In einer zweitquantisierten Singleton-Theorie mit N Arten enthält der Singulett-Zweiteilchensektor genau "
       "ein D(3,2), unabhängig von N; das ersetzte P_mode, P_paar und P_chir. Nativ ist davon nichts: Der Dubleton ist "
       "genau der Zwei-Singleton-Raum ohne Vakuum, D(3,2) ist dort kein Singulett der Kopiengruppe, jede Gruppe auf den "
       "Cartan-Strömen wirkt als Drehung, und l^2/G wäre pi N/32 mit N als neuer freier Zahl.")

# ====================================================================================================================
say("\nTeil 5: Bilanz und Urteil")
ol = (REPO / "research" / "originlab" / "reconstruction_search.py").read_text()
check("[E] OriginLab (nur gelesen): Die Klausel 'graviton' verlangt den Status 'supported' der Gravitationsachse. Nach "
      "z158 bleibt er bedingt: Kein nativer Zwei-Naht-Schritt tauscht den Sektor als Teilchen aus. Die Klausel bleibt "
      "leer, nichts wird geändert", '("graviton", st["gravity"] == "supported")' in ol)
bilanz = [
    ("E", "Freier Transfer",
     "Strommode -> (J^A + J^B)/2 + gemischter Transferstrom (Gewichte 1/4, 1/4, 1/2). r_1 = 1/2, "
     "r_4 = 1/16 + 3/(8m) + 3/(16m^2) - 1/(6m^3): nicht null, nur für m -> unendlich faktorisiert; der verbundene Teil "
     "ist die su(2)_1-Struktur der Transferströme. Das Quintett wird nicht als Multiplett übertragen."),
    ("E", "Gatter",
     "Vertauscht nicht mit den Strommoden einer Naht. Bei gerader Ladung (1 - i)/2 + (1 + i)/2 mal Nahttausch: Antwort "
     "(1 + i)/2 für jeden neutralen Operator, unabhängig von Quantenzahl, Mode, Spin. Bei ungerader Ladung zwei freie "
     "Vierteltransfers. Kein Pol, kein gebundener Zustand."),
    ("N", "Nativer Zustand und Kette",
     "Bei beta = 4 pi dasselbe Bild mit Korrekturen der Größe 0,04 bis 0,08. In der Kette aus z152 fällt der "
     "Fermionkorrelator wie (sqrt 2 - 1)^r, jede neutrale Antwort mindestens wie (sqrt 2 - 1)^(2r); ein Delta ist dort "
     "nicht definiert."),
    ("E", "Dimension",
     "Delta = 4 = vier Stromquanten im Rahmen der Naht (Möbius); Delta = 3 im Rahmen der so(3,2) (Original + Partner), "
     "dort von selbst. Kein Widerspruch: E_0 = 3 = 1 + 4/2, die vier Quanten sind die Anregung von zwei Singletons."),
    ("E", "Kriterium aus z156",
     "Falsch gestellt. Faktorisierung widerlegt P_X nicht (Flato-Fronsdal; nachgerechnet auf der Linie), und "
     "C(n + 2, 5) bestätigt P_X nicht."),
    ("E/L", "Singleton-Lesart",
     "Singulett-Zweiteilchensektor von N Arten: genau ein D(3,2) für jedes N. Nativ fehlt die zweite Quantisierung des "
     "Singletons (P_2Q); die Kopiengruppe des Dubletons hat D(3,2) als Dublett, nicht als Singulett."),
    ("E/L", "Native Gruppen",
     "S_8 in W(E7) (zustandserhaltend), W(E8), O(8): alle mischen Cartan-Ströme und wirken auf ein Paar als Drehung der "
     "so(3,2). Keine Artengruppe im Bestand."),
    ("L", "G",
     "l^2/G = pi N/32; nativ planckgroß; der Ledgerwert verlangt N = 1024 pi^2 e^(2/alpha), dasselbe wie die "
     "Nahtzahl je Horizont in anderer Einheit. P_G bleibt eine freie Zahl (als P_N)."),
]
for st, was, txt in bilanz:
    say(f"   [{st}] {was}: {txt}")
urteil = (
    "Nein als eigenes Teilchen, gescopt; P_X ist für den Dubleton weder bestätigt noch widerlegt, sondern "
    "umformuliert. (1) Übertragen wird der Sektor von beiden nativen Zwei-Naht-Schritten, aber von keinem als Teilchen. "
    "Unter dem freien Vierteltransfer geht eine Strommode in (J^A + J^B)/2 plus einen gemischten Transferstrom, der "
    "keiner Naht gehört. Die Vierstromantwort ist r_4 = 1/16 + 3/(8m) + 3/(16m^2) - 1/(6m^3) (m = 1: 11/24): führend "
    "vier unabhängige Antworten 1/2, die Abweichung ist die su(2)_1-Struktur der Transferströme, und die fünf "
    "Komponenten des Quintetts werden verschieden übertragen (die Drehgruppe der so(3,2) ist keine Symmetrie des "
    "Transfers). Unter dem z120-Gatter ist die Antwort im Vakuum für jeden neutralen Operator (1 + i)/2: Das Gatter "
    "überlagert bei gerader Ladung 'nichts' und 'ganze Naht getauscht'. Das faktorisiert nicht, unterscheidet aber "
    "weder D noch Spin noch Mode, und zwischen den Nähten liegt keine Strecke. In der Kette aus z152 fällt jede Antwort "
    "exponentiell. Gilt für: eine und zwei Stromarten, m = 1 bis 6 "
    "(Vakuum exakt), m = 1 im nativen Zustand, P_sym. (2) Das Kriterium aus z156 ist falsch gestellt. Delta = 4 (vier "
    "Stromquanten) und Delta = 3 (D(3,2)) gelten zugleich, für zwei verschiedene sl(2): die Möbius-Gruppe der Naht "
    "verschiebt die Mode, die so(3,2) erzeugt bei fester Mode Paare aus Original und Partner. Aus Sicht der so(3,2) "
    "sind die vier Quanten die Anregung eines Zwei-Singleton-Zustands (E_0 = 1 + 4/2), und ein nach Wick faktorisierter "
    "Korrelator kann exakt die Antwort von Delta = 3 sein. Faktorisierung widerlegt also nichts, und die Linienantwort "
    "C(n + 2, 5) ist reine Darstellungstheorie. (3) Was P_X tatsächlich verlangt, ist die zweite Quantisierung des "
    "Singletons: ein Randfeld mit invariantem Vakuum und N Arten (P_2Q), dessen Singuletts physikalisch sind (P_sing). "
    "Dann enthält der Zweiteilchensektor für jedes N genau ein D(3,2), und P_mode, P_paar, P_chir und P_1T entfielen "
    "(P_chir um den Preis, dass die ausgetauschte Mode eine Diagonale und nicht modulinvariant ist); l^2/G wäre "
    "pi N/32 [L], also P_G eine Anzahl (P_N), nativ planckgroß. Nativ ist das nicht: Der Dubleton ist genau "
    "der Raum zweier Singletons ohne Vakuum, D(3,2) ist dort das O(2)-Dublett mit Ladung +-4 und nie ein Singulett, "
    "mehrere Kopien unter gemeinsamer so(3,2) geben den Skalar D(N_c/2, 0), jede Gruppe im Bestand, die "
    "Cartan-Ströme mischt, wirkt als Drehung, und die Tomita-Konjugation tauscht die beiden Singletons, sodass ein "
    "einzelner in keiner Nahtalgebra liegt. (4) Prämissen: P_X (eng) wird P_2Q + P_sing; dafür entfielen P_1T, P_mode, "
    "P_paar, P_chir. Das sind zwei statt fünf Namen, aber P_2Q ist die stärkere Annahme (sie setzt die Randtheorie der "
    "Dualität voraus). P_TP oder P_CTM (physikalischer Partner) bleibt nötig, P_G bleibt frei (als P_N). P_dyn, P_h2 "
    "und P_curv sind als anderer Weg unverändert. Ein freier Austausch wäre auch unter P_2Q kein Graviton mit "
    "eigener Kopplung: Die Kopplung ist dort 1/N. K13 bleibt L, die Klausel 'graviton' bleibt leer, K15 ist nicht "
    "berührt, die neun Prämissen des Kandidaten sind unverändert; nichts ist befördert. Nächster Test: ob es im "
    "Bestand ein Objekt mit so(3,2)-invariantem Vakuum gibt, dessen Einteilchenraum ein einzelner Singleton ist (ein "
    "halber Dubleton; in der Gitterlesart eine Mode, die beide Halbketten mischt); und, falls ja, die "
    "Dreipunktkopplung zweier Nahtoperatoren an das Singulett-Bilinear, deren Quadrat 1/C_T sein müsste. Ohne ein "
    "solches Objekt ist der Weg über den Dubleton für P_X geschlossen und es bleibt P_dyn oder P_h2.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z158_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z158.txt").write_text("\n".join(OUT) + "\n")
