"""Z106: Folgeschritt N7 'Gradierte Komposition aus einer gemeinsamen CAR' (Hypothese H7).

Leitlinie 'Rekonstruktion statt Auswahl': Die Komposition der Nähte wird aus der CAR gebaut, die das Original schon für seine
beiden Quarter benutzt, nicht aus der gewünschten 4D-Statistik. Firewall DIMENSION.UPLIFT.FIREWALL.01: 4D nur als Lesart über
die radiale Standardform (z83), der Spin j über deren räumliche Drehungen; keine wechselwirkende 3+1-Schließung. Gelöste
Bausteine (16, E8, Flavour) bleiben unberührt, nichts wird befördert. Bestand: z91 (Sym^N der Einteilchenleiter = freier
4D-Fockraum bis E = N; Paaranregungen verschiedener Nähte vertauschen in einer gemeinsamen CAR; Prämisse P_ident;
W_12 = (-1)^(n_1 n_2) als Unterschied von fermionischem und gewöhnlichem Nahttausch), z97 (gekleidete Vertizes
Psi_s = V^+ (x) S_s erhalten G = U_2pi R^2; die native bosonische Komposition gibt zwei Weyl-Quanten an der Schwelle j = 1 statt
j = 0), z99 (zwei Nähte in einer gemeinsamen CAR: Gitter Z^2, Transfer psi*_1 psi_2 = su(2)_1-Strom).

Hypothese H7: In einer gemeinsamen CAR über die Nähte antikommutieren die gekleideten Quanten (ein Ladungsschritt, ungerade
Fermionzahl) von selbst; die gradierte Komposition wäre nativ.

Aufbau: Naht i mit den Quartern a (alpha = 1/4) und b (alpha = 3/4), je Quarter die Moden x = alpha - 1 (See) und x = alpha.
Ladungsschritt V^+ = U_q = psi*_alpha n_(alpha-1) + psi*_(alpha-1) (1 - n_alpha), die Maya-Verschiebung |Q> -> |Q+1> im Fenster
(ungerade, aus Moden des Quarters gebaut). Je Naht die radiale Standardform (z83/z97, CCR-Schritte, M = 2) als gerader Faktor.
Alle Moden aller Nähte in einer Jordan-Wigner-CAR; drei Nähte mit Quarter a allein (je zwei Moden).

Teil 1  Gekleidete Vertizes in der gemeinsamen CAR (zwei und drei Nähte): Moden-Einbettung der Nahtalgebren
        (Majorana-Entwicklung) = Jordan-Wigner-Kette, Austausch der gekleideten Vertizes Psi^(i)_s, Vertauschen der
        Paaranregungen und radialen Heber, Super-Tensorprodukt (gradierte Vertauschung, Koszul-Regel, Basis aus Monomen,
        drei Nähte) gegen die ungradierte Einbettung.
Teil 2  Zwei (drei) Weyl-Quanten an der Schwelle E = 3 (9/2): kanonischer Nahttausch Gamma(pi) gegen den gewöhnlichen Tausch,
        Casimir; die Familie der gaußschen Tauschoperatoren Gamma(u); Koszul-Vorzeichen gegen Spin-Statistik und G = U_2pi R^2;
        Zählung der gradierten symmetrischen Potenz gegen den freien 4D-Fockraum (z91) und gegen Sym^N.
Teil 3  Ist die gemeinsame CAR nativ? Hardy-Komposition (eine Polarisation über alle Nähte; Gamma(pi) aus den Nullmoden der
        Transferströme von z99; Dublettwirkung; Hüpfterm), Superlokalität des Gitters Z^(2N) (Kozyklus = Jordan-Wigner-Kette),
        Beobachtbarkeit der Kette (Ordnungen, Korrelationen, Kommutant und Superkommutant), gerade Daten gegen die
        Klein-Transformation; Bilanz, Urteil H7 und Folge für die Statistikprämisse der Szenario-Suche.
"""
import itertools
import json
import time
from fractions import Fraction as Fr
from math import comb
from pathlib import Path

import numpy as np
import scipy.sparse as sps
from scipy.linalg import expm

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
rng = np.random.default_rng(106)
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


def nrm(X):
    """Größter Betrag eines dünnen oder dichten Operators (0 für den Nulloperator)."""
    if sps.issparse(X):
        X = X.tocsr()
        X.eliminate_zeros()
        return float(abs(X).max()) if X.nnz else 0.0
    return float(np.abs(X).max()) if np.size(X) else 0.0


def kr(*ops):
    out = sps.csr_matrix(ops[0])
    for o in ops[1:]:
        out = sps.kron(out, sps.csr_matrix(o), format="csr")
    return out


# ------------------------------------------------------------------ radiale Standardform (wie z83/z97, CCR-Schritte)
AMP_CCR = (lambda m, k: np.sqrt(k + 1), lambda m, k: np.sqrt(m - k))
SCHRITT = (lambda k: k + 1, lambda k: k)                          # S_0 Teilchen hoch, S_1 Loch runter (z83)


def standardform(M):
    """(+)_(m, mquer <= M) H_m (x) H-quer_mquer wie z97: P[a,b] = S_a (x) S-quer_b, A = Bilineare links, B rechts, E, und die
    einseitigen Schritte SL[a] = S_a (x) 1, SR[b] = 1 (x) S-quer_b."""
    lab = [(m, mb, k, kb) for m in range(1, M + 1) for mb in range(1, M + 1) for k in range(m) for kb in range(mb)]
    ix = {s: i for i, s in enumerate(lab)}
    n = len(lab)

    def dicht(r, c, w):
        return sps.csr_matrix((w, (r, c)), shape=(n, n)).toarray()

    P = {}
    for a in range(2):
        for b in range(2):
            r, c, w = [], [], []
            for (m, mb, k, kb), i in ix.items():
                if m < M and mb < M:
                    r.append(ix[(m + 1, mb + 1, SCHRITT[a](k), SCHRITT[b](kb))])
                    c.append(i)
                    w.append(AMP_CCR[a](m, k)*AMP_CCR[b](mb, kb))
            P[(a, b)] = dicht(r, c, w)

    def bilinear(c_, a_, seite):
        r, cc, w = [], [], []
        for (m, mb, k, kb), i in ix.items():
            mm, kk = (m, k) if seite == 0 else (mb, kb)
            k1 = kk - (1 if a_ == 0 else 0)
            if mm == 1 or not 0 <= k1 <= mm - 2:
                continue
            k2 = SCHRITT[c_](k1)
            r.append(ix[(m, mb, k2, kb) if seite == 0 else (m, mb, k, k2)])
            cc.append(i)
            w.append(AMP_CCR[a_](mm - 1, k1)*AMP_CCR[c_](mm - 1, k1))
        return dicht(r, cc, w)

    A = {(c_, a_): bilinear(c_, a_, 0) for c_ in range(2) for a_ in range(2)}
    B = {(c_, a_): bilinear(c_, a_, 1) for c_ in range(2) for a_ in range(2)}
    SL, SR = {}, {}
    for a in range(2):
        r, c, w, r2, c2, w2 = [], [], [], [], [], []
        for (m, mb, k, kb), i in ix.items():
            if m < M:
                r.append(ix[(m + 1, mb, SCHRITT[a](k), kb)])
                c.append(i)
                w.append(AMP_CCR[a](m, k))
            if mb < M:
                r2.append(ix[(m, mb + 1, k, SCHRITT[a](kb))])
                c2.append(i)
                w2.append(AMP_CCR[a](mb, kb))
        SL[a], SR[a] = dicht(r, c, w), dicht(r2, c2, w2)
    return lab, ix, P, A, B, np.diag([0.5*(s[0] + s[1]) for s in lab]), SL, SR


M_RAD = 2
LAB, IXS, PS, AS, BS, ES, SL, SR = standardform(M_RAD)
NR = len(LAB)
MV, MBV = np.array([s[0] for s in LAB]), np.array([s[1] for s in LAB])
SG = np.where((MV - MBV) % 2 == 0, 1.0, -1.0)                      # U_2pi = (-1)^(2h) auf der Standardform (z97)
DG = {(c_, a_): AS[(c_, a_)] - BS[(a_, c_)] for c_ in range(2) for a_ in range(2)}   # räumliche Drehungen (z83/z97)
JR = [(DG[(0, 1)] + DG[(1, 0)])/2, (DG[(0, 1)] - DG[(1, 0)])/2j, (DG[(0, 0)] - DG[(1, 1)])/2]
E11 = np.zeros(NR)
E11[IXS[(1, 1, 0, 0)]] = 1.0                                       # skalarer Grundzustand (E, h) = (1, 0)
I9 = np.eye(NR)
SWAP_P = np.zeros((NR*NR, NR*NR))
for p0 in range(NR):
    for p1 in range(NR):
        SWAP_P[p1*NR + p0, p0*NR + p1] = 1.0
# kleiner Faktor P' = span{(1,1), (2,1)_k} für drei Nähte (unter S_s und J abgeschlossen)
SUB = [IXS[(1, 1, 0, 0)], IXS[(2, 1, 0, 0)], IXS[(2, 1, 1, 0)]]
SLs = {a: SL[a][np.ix_(SUB, SUB)] for a in range(2)}
JRs = [J[np.ix_(SUB, SUB)] for J in JR]


# ------------------------------------------------------------------ Nähte in einer gemeinsamen CAR
class Naehte:
    """Nähte 0..N-1 mit den Quartern 'quarter' (a: alpha = 1/4, b: alpha = 3/4), je Quarter die Moden n = -1 (See,
    x = alpha - 1) und n = 0 (x = alpha), alle in einer Jordan-Wigner-CAR. Die Bits sind kanonisch (Naht, Quarter, n)
    vergeben; 'reihenfolge' legt nur die Jordan-Wigner-Ordnung der Nähte fest. ungradiert=True: jede Naht mit eigener Kette,
    Moden verschiedener Nähte vertauschen (gewöhnliche Tensorkomposition)."""

    def __init__(self, N, quarter=("a", "b"), reihenfolge=None, ungradiert=False):
        self.N, self.quarter = N, tuple(quarter)
        self.labels = [(i, q, n) for i in range(N) for q in self.quarter for n in (-1, 0)]
        self.bit = {l: j for j, l in enumerate(self.labels)}
        self.NB = len(self.labels)
        self.dim = 1 << self.NB
        reihe = list(range(N)) if reihenfolge is None else list(reihenfolge)
        self.jw = [l for i in reihe for l in self.labels if l[0] == i]
        self.pos = {l: p for p, l in enumerate(self.jw)}
        s = np.arange(self.dim)
        occ = (s[:, None] >> np.arange(self.NB)[None, :]) & 1
        self.occ = occ
        self.c = {}
        for p, l in enumerate(self.jw):
            j = self.bit[l]
            vor = [self.bit[m] for m in self.jw[:p] if not ungradiert or m[0] == l[0]]
            par = occ[:, vor].sum(1) % 2 if vor else np.zeros(self.dim, dtype=np.int64)
            src = s[occ[:, j] == 1]
            self.c[l] = sps.csr_matrix(((1 - 2*par[src]).astype(float), (src ^ (1 << j), src)),
                                       shape=(self.dim, self.dim))
        self.cd = {l: X.T.tocsr() for l, X in self.c.items()}
        self.nop = {l: sps.diags(occ[:, self.bit[l]].astype(float)).tocsr() for l in self.labels}
        self.I = sps.identity(self.dim, format="csr")
        self.f = np.stack([occ[:, [self.bit[l] for l in self.labels if l[0] == i]].sum(1) for i in range(N)], 1)
        self.Qq = {(i, q): occ[:, self.bit[(i, q, 0)]] - (1 - occ[:, self.bit[(i, q, -1)]])
                   for i in range(N) for q in self.quarter}
        self.Q = np.stack([sum(self.Qq[(i, q)] for q in self.quarter) for i in range(N)], 1)
        self.vi = sum(1 << self.bit[l] for l in self.labels if l[2] == -1)
        self.vac = np.zeros(self.dim)
        self.vac[self.vi] = 1.0

    def U(self, i, q="a"):
        """Ladungsschritt V^+ auf dem Quarter q der Naht i: Maya-Verschiebung |Q> -> |Q+1> im Fenster (ungerade)."""
        hi, lo = (i, q, 0), (i, q, -1)
        return (self.cd[hi] @ self.nop[lo] + self.cd[lo] @ (self.I - self.nop[hi])).tocsr()

    def R2(self, i):
        """R_i^2 = (-1)^(Q_a + Q_b) = (-1)^F der Naht i (z82)."""
        return sps.diags((-1.0)**self.Q[:, i]).tocsr()

    def tausch(self, perm, kanonisch=True):
        """Nahtpermutation (Naht i -> perm[i]): kanonisch = zweite Quantisierung Gamma(perm) (Moden auf Moden, fermionisches
        Umordnungsvorzeichen in der Jordan-Wigner-Ordnung, auf den Fermisee normiert); sonst gewöhnliche Bitvertauschung."""
        rows, cols, vals = [], [], []
        for s in range(self.dim):
            occl = sorted((l for l in self.labels if (s >> self.bit[l]) & 1), key=lambda l: self.pos[l])
            bild = [(perm[l[0]], l[1], l[2]) for l in occl]
            rows.append(sum(1 << self.bit[b] for b in bild))
            cols.append(s)
            inv = sum(1 for x_, y_ in itertools.combinations([self.pos[b] for b in bild], 2) if x_ > y_)
            vals.append(-1.0 if (kanonisch and inv % 2) else 1.0)
        T = sps.csr_matrix((vals, (rows, cols)), shape=(self.dim, self.dim))
        return (T*T[self.vi, self.vi]).tocsr() if kanonisch else T


def majorana(sys_, labs):
    out = []
    for l in labs:
        out += [(sys_.c[l] + sys_.cd[l]).tocsr(), (1j*(sys_.cd[l] - sys_.c[l])).tocsr()]
    return out


def monome(gs):
    """Alle geordneten Majorana-Monome gamma_S (S aufsteigend) als Liste von Matrizen."""
    out = [sps.identity(gs[0].shape[0], format="csr", dtype=complex)]
    for g in gs:
        out += [(X @ g).tocsr() for X in out]
    return out


def einbettung(X, mon_lok, mon_glob):
    """Moden-Einbettung iota(X): Entwicklung nach lokalen Majorana-Monomen, Ersetzen durch die globalen Moden der Naht."""
    Xd = X.toarray() if sps.issparse(X) else np.asarray(X)
    d = Xd.shape[0]
    out = sps.csr_matrix(mon_glob[0].shape, dtype=complex)
    for G, Gg in zip(mon_lok, mon_glob):
        a = np.vdot(G.toarray().ravel(), Xd.ravel())/d
        if abs(a) > 1e-14:
            out = out + a*Gg
    return out.tocsr()


def rest_linear(X, basis):
    """Relativer Abstand eines Operators vom linearen Spann der Moden (0 = gaußsch-linear)."""
    A = np.array([B.toarray().ravel() for B in basis]).T
    x = X.toarray().ravel()
    k, *_ = np.linalg.lstsq(A, x, rcond=None)
    return float(np.linalg.norm(A @ k - x)/np.linalg.norm(x))


# =================================================================== Teil 1
say("Teil 1: Gekleidete Vertizes in einer gemeinsamen CAR (zwei und drei Nähte)")

# --- 1.1 Moden-Einbettung der Nahtalgebren und Super-Tensorprodukt
Z1, Z2, Z2u = Naehte(1), Naehte(2), Naehte(2, ungradiert=True)
mon_lok = monome(majorana(Z1, Z1.labels))
mon_glob = {i: monome(majorana(Z2, [(i, q, n) for (_, q, n) in Z1.labels])) for i in range(2)}
PI1 = sps.diags((-1.0)**Z1.f[:, 0]).tocsr()                     # Parität einer Naht (16 Zustände)
I16 = sps.identity(16, format="csr")


def zufall_hom(par):
    X = sps.csr_matrix(rng.normal(size=(16, 16)) + 1j*rng.normal(size=(16, 16)))
    return ((X + (1 if par == 0 else -1)*(PI1 @ X @ PI1))/2).tocsr()


iota = lambda X, i: einbettung(X, mon_lok, mon_glob[i])
ein_fehler, hom_fehler = 0.0, 0.0
for par in (0, 1):
    for _ in range(2):
        X = zufall_hom(par)
        ein_fehler = max(ein_fehler, nrm(iota(X, 0) - sps.kron(I16, X)), nrm(iota(X, 1) - sps.kron(X, PI1 if par else I16)))
for _ in range(2):
    X, Y = (sps.csr_matrix(rng.normal(size=(16, 16)) + 1j*rng.normal(size=(16, 16))) for _ in range(2))
    for i in range(2):
        hom_fehler = max(hom_fehler, nrm(iota(X @ Y, i) - iota(X, i) @ iota(Y, i)))
moden_ok = max(nrm(iota(Z1.c[(0, q, n)], i) - Z2.c[(i, q, n)]) for i in range(2) for (_, q, n) in Z1.labels) < TOL
gr_fehler, kos_fehler = 0.0, 0.0
for pa, pb, pc in itertools.product((0, 1), repeat=3):
    A_, B_, A2 = zufall_hom(pa), zufall_hom(pb), zufall_hom(pc)
    B2 = zufall_hom(pa)
    a0, b1, a20, b21 = iota(A_, 0), iota(B_, 1), iota(A2, 0), iota(B2, 1)
    gr_fehler = max(gr_fehler, nrm(a0 @ b1 - (-1)**(pa*pb)*(b1 @ a0)))
    kos_fehler = max(kos_fehler, nrm(a0 @ b1 @ a20 @ b21 - (-1)**(pb*pc)*(iota(A_ @ A2, 0) @ iota(B_ @ B2, 1))))
cu0, cu1 = Z2u.c[(0, "a", 0)], Z2u.c[(1, "a", 0)]
ungr_anti = nrm(cu0 @ cu1 + cu1 @ cu0)
ungr_komm = nrm(cu0 @ cu1 - cu1 @ cu0)
klein = max(nrm(Z2u.c[(1, q, n)] - sps.diags((-1.0)**Z2.f[:, 0]) @ Z2.c[(1, q, n)]) for (_, q, n) in Z1.labels)
# drei Nähte (Quarter a): Monome der Nähte multiplizieren sich zu den globalen Monomen, die eine Basis bilden
Z3 = Naehte(3, ("a",))
Z1a = Naehte(1, ("a",))
g3 = majorana(Z3, Z3.labels)
mon3 = {i: monome(g3[4*i:4*i + 4]) for i in range(3)}
mon3_glob = monome(g3)
prod_ok, spur_max = True, 0.0
for S0, X0 in enumerate(mon3[0]):
    for S1, X1 in enumerate(mon3[1]):
        for S2, X2 in enumerate(mon3[2]):
            G = mon3_glob[S0 + 16*S1 + 256*S2]
            prod_ok &= nrm(X0 @ X1 @ X2 - G) < TOL
            if S0 + S1 + S2:
                spur_max = max(spur_max, abs(G.diagonal().sum()))
par3 = [np.array([bin(S).count("1") % 2 for S in range(16)]) for _ in range(3)]
gr3, kos3 = 0.0, 0.0
for _ in range(6):
    el = []
    for i in range(3):
        p = int(rng.integers(2))
        sel = [S for S in range(16) if par3[i][S] == p]
        co = rng.normal(size=len(sel)) + 1j*rng.normal(size=len(sel))
        el.append((p, sum(c_*mon3[i][S] for c_, S in zip(co, sel)).tocsr()))
    for (i, (pi_, Xi)), (j, (pj_, Xj)) in itertools.combinations(enumerate(el), 2):
        gr3 = max(gr3, nrm(Xi @ Xj - (-1)**(pi_*pj_)*(Xj @ Xi)))
    (p0, x0), (p1, x1), (p2, x2) = el
    el2 = [(int(rng.integers(2)), None) for _ in range(3)]
    y = []
    for i, (p, _) in enumerate(el2):
        sel = [S for S in range(16) if par3[i][S] == p]
        co = rng.normal(size=len(sel)) + 1j*rng.normal(size=len(sel))
        y.append((p, sum(c_*mon3[i][S] for c_, S in zip(co, sel)).tocsr()))
    (q0, y0), (q1, y1), (q2, y2) = y
    links = (x0 @ x1 @ x2) @ (y0 @ y1 @ y2)
    rechts = (-1)**(p1*q0 + p2*q0 + p2*q1)*((x0 @ y0) @ (x1 @ y1) @ (x2 @ y2))
    kos3 = max(kos3, nrm(links - rechts))
ZAHLEN["einbettung"] = {"iota_i(X) - S^|X| (x) X (Majorana-Entwicklung, zwei Nähte)": ein_fehler,
                        "Homomorphiefehler": hom_fehler, "gradierte Vertauschung (zwei Nähte)": gr_fehler,
                        "Koszul-Regel (zwei Nähte)": kos_fehler, "ungradiert: |{c_1, c_2}|, |[c_1, c_2]|": [ungr_anti, ungr_komm],
                        "Klein: c^u_2 - (-1)^(F_1) c_2": klein, "drei Nähte: max |Tr gamma_U|, U != leer": spur_max,
                        "drei Nähte: gradierte Vertauschung, Koszul": [gr3, kos3]}
check("[E] Moden-Einbettung = Jordan-Wigner-Kette: Die Nahtalgebra B(F_i) (vier Moden, Quarter a und b) wird durch ihre "
      "Majorana-Entwicklung in die gemeinsame CAR eingebettet, die lokalen Moden durch die globalen ersetzt. Das Ergebnis ist "
      f"ein Homomorphismus (Fehler {hom_fehler:.0e}) mit iota_i(c_loc) = c_glob, und für homogene X ist iota_1(X) = "
      f"(-1)^(F_0 |X|) (x) X exakt (Fehler {ein_fehler:.0e}): die Kette ist die Parität der früheren Naht, nur an ungeraden "
      "Elementen. Die gewöhnliche Tensoreinbettung ist dasselbe mit c^u_(1,x) = (-1)^(F_0) c_(1,x) (Klein-Transformation, Fehler "
      f"{klein:.0e}); dort vertauschen Moden verschiedener Nähte ({ungr_komm:.0e}) statt zu antikommutieren ({ungr_anti:g})",
      moden_ok and ein_fehler < TOL and hom_fehler < 1e-10 and klein < TOL and ungr_komm < TOL and ungr_anti > 1)
check("[E] Super-Tensorprodukt: Für zufällige homogene Elemente zweier Nahtalgebren gilt iota_0(x) iota_1(y) = "
      f"(-1)^(|x||y|) iota_1(y) iota_0(x) (Fehler {gr_fehler:.0e}) und die Koszul-Regel (x (x) y)(x' (x) y') = "
      f"(-1)^(|y||x'|) xx' (x) yy' (Fehler {kos_fehler:.0e}); die Bilder der lokalen Moden sind genau die globalen Moden, die "
      "die ganze CAR erzeugen, und 4^4 * 4^4 = 4^8: die gemeinsame CAR ist B(F_0) (x)^ B(F_1). Drei Nähte: die 4096 Produkte "
      "der Nahtmonome sind genau die globalen Majorana-Monome (spurfrei außer 1, max |Tr| = "
      f"{spur_max:.0e}, also eine Orthogonalbasis), gradierte Vertauschung aller Paare und Koszul-Regel mit drei Faktoren "
      f"(Fehler {gr3:.0e}, {kos3:.0e}): (A_0 (x)^ A_1) (x)^ A_2 = A_0 (x)^ A_1 (x)^ A_2",
      gr_fehler < 1e-10 and kos_fehler < 1e-9 and prod_ok and spur_max < 1e-9 and gr3 < 1e-10 and kos3 < 1e-9)

# --- 1.2 gekleidete Vertizes und Paaranregungen
DIMF = Z2.dim
U0, U1 = Z2.U(0, "a"), Z2.U(1, "a")
PSI = {(0, s): kr(U0, SL[s], I9) for s in range(2)}
PSI.update({(1, s): kr(U1, I9, SL[s]) for s in range(2)})
PSIQ = {(0, t): kr(U0.T, SR[t], I9) for t in range(2)}
PSIQ.update({(1, t): kr(U1.T, I9, SR[t]) for t in range(2)})
anti_max = max(nrm(PSI[(0, s)] @ X + X @ PSI[(0, s)]) for s in range(2) for X in [PSI[(1, t)] for t in range(2)]
               + [PSIQ[(1, t)] for t in range(2)])
nicht_null = min(nrm(PSI[(0, s)] @ PSI[(1, t)]) for s in range(2) for t in range(2))
paare = [kr(Z2.cd[(i, q, 0)] @ Z2.c[(i, qq, -1)], I9, I9) for i in range(2) for q in "ab" for qq in "ab"]
heber = {(i, s, t): (PSI[(i, s)] @ PSIQ[(i, t)]).tocsr() for i in range(2) for s in range(2) for t in range(2)}
heber_form = max(nrm(heber[(0, s, t)] - kr(U0 @ U0.T, PS[(s, t)], I9)) for s in range(2) for t in range(2))
komm_paar = 0.0
for X in paare[:4] + [heber[(0, s, t)] for s in range(2) for t in range(2)]:
    for Y in paare[4:] + [heber[(1, s, t)] for s in range(2) for t in range(2)] + [PSI[(1, s)] for s in range(2)]:
        komm_paar = max(komm_paar, nrm(X @ Y - Y @ X))
# drei Nähte: Vertizes auf F_3 (x) P'^3
I3 = np.eye(3)
P3OPS = lambda i, X: [X if j == i else I3 for j in range(3)]
PSI3 = {(i, s): kr(Z3.U(i, "a"), *P3OPS(i, SLs[s])) for i in range(3) for s in range(2)}
anti3 = max(nrm(PSI3[(i, s)] @ PSI3[(j, t)] + PSI3[(j, t)] @ PSI3[(i, s)])
            for i in range(3) for j in range(3) if i != j for s in range(2) for t in range(2))
paar3 = [kr(Z3.cd[(i, "a", 0)] @ Z3.c[(i, "a", -1)], I3, I3, I3) for i in range(3)]
komm3 = max(nrm(paar3[i] @ PSI3[(j, s)] - PSI3[(j, s)] @ paar3[i]) for i in range(3) for j in range(3) if i != j
            for s in range(2))
ZAHLEN["vertizes"] = {"max |{Psi^(0), Psi^(1) bzw. Psi-quer^(1)}|": anti_max, "min |Psi^(0) Psi^(1)|": nicht_null,
                      "Heber Psi_s Psi-quer_t - (U U^+) (x) P_st": heber_form,
                      "max |[Paar/Heber^(0), Paar/Heber/Psi^(1)]|": komm_paar,
                      "drei Nähte: max |{Psi^(i), Psi^(j)}|, max |[Paar^(i), Psi^(j)]|": [anti3, komm3]}
check("[E] Gekleidete Vertizes Psi^(i)_s = V^+_a (x) S_s (Ladungsschritt mal Schrittart, z97) und Psi-quer^(i)_t = V^-_a (x) "
      "S-quer_t verschiedener Nähte antikommutieren in der gemeinsamen CAR exakt (max |{Psi^(0), Psi^(1)}| = "
      f"{anti_max:.0e} bei |Psi^(0) Psi^(1)| = {nicht_null:g}; drei Nähte, alle Paare: {anti3:.0e}). Die Paaranregungen "
      "psi*_p psi_r einer Naht und die radialen Heber Psi_s Psi-quer_t = (U U^+) (x) P_st (Fehler "
      f"{heber_form:.0e}), Produkte zweier ungerader Vertizes derselben Naht, vertauschen mit allem auf der anderen Naht "
      f"(max {komm_paar:.0e}; drei Nähte {komm3:.0e}): Die Statistik trägt allein der Ladungsschritt, die Schrittart ist gerade",
      anti_max < TOL and nicht_null > 0.5 and heber_form < TOL and komm_paar < TOL and anti3 < TOL and komm3 < TOL)
befund("Teil 1: H7 trifft in der gemeinsamen CAR zu. Die Moden-Einbettung jeder Nahtalgebra ist die Jordan-Wigner-Kette, die "
       "Gesamtalgebra ist das Super-Tensorprodukt der Nahtalgebren (zwei und drei Nähte). Gekleidete Vertizes verschiedener "
       "Nähte antikommutieren, Paaranregungen und radiale Heber sind gerade und vertauschen; die gewöhnliche Komposition ist "
       "die Klein-Transformierte mit der Parität der früheren Naht an den ungeraden Elementen", True)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

# =================================================================== Teil 2
say("\nTeil 2: Weyl-Quanten an der Schwelle, Tauschoperatoren, Spin-Statistik und Zählung")

# --- 2.1 kanonischer Nahttausch gegen gewöhnlichen
GAM_F = Z2.tausch([1, 0])
GEW_F = Z2.tausch([1, 0], kanonisch=False)
W12 = sps.diags((-1.0)**(Z2.Q[:, 0]*Z2.Q[:, 1])).tocsr()
gam_moden = max(nrm(GAM_F @ Z2.c[(i, q, n)] @ GAM_F.T - Z2.c[(1 - i, q, n)]) for (i, q, n) in Z2.labels)
gam_inv = nrm(GAM_F @ GAM_F - Z2.I)
gew_rel = nrm(GEW_F - GAM_F @ W12)
alle_moden = [Z2.c[l] for l in Z2.labels] + [Z2.cd[l] for l in Z2.labels]
gew_bild = GEW_F @ Z2.c[(0, "a", 0)] @ GEW_F.T
gew_rest = rest_linear(gew_bild, alle_moden)
gew_form = nrm(gew_bild - Z2.c[(1, "a", 0)] @ sps.diags((-1.0)**Z2.Q[:, 0]))
T_FERM = kr(GAM_F, SWAP_P)
T_GEW = kr(GEW_F, SWAP_P)
BASE = np.kron(np.kron(Z2.vac, E11), E11)
ZAHLEN["tausch"] = {"Gamma(pi) c_(i,x) Gamma^+ - c_(1-i,x)": gam_moden, "Gamma^2 - 1": gam_inv,
                    "T_gew - Gamma W_12, W_12 = (-1)^(Q_0 Q_1)": gew_rel, "T_gew c_0 T_gew: Rest zum Modenspann": gew_rest,
                    "T_gew c_0 T_gew - c_1 (-1)^(Q_0)": gew_form}
check("[E] Kanonischer Nahttausch: Gamma(pi), die zweite Quantisierung der Einteilchenvertauschung (Moden auf Moden mit "
      f"Umordnungsvorzeichen, Fermisee fest), bildet c_(0,x) <-> c_(1,x) exakt ab (Fehler {gam_moden:.0e}) und ist eine "
      f"Involution ({gam_inv:.0e}). Der gewöhnliche Tausch der Bits ist T_gew = Gamma(pi) W_12 mit W_12 = (-1)^(Q_0 Q_1) "
      f"(Fehler {gew_rel:.0e}, z91) und nicht gaußsch: T_gew c_0 T_gew = c_1 (-1)^(Q_0) (Fehler {gew_form:.0e}), relativer "
      f"Abstand zum linearen Modenspann {gew_rest:.2f}",
      gam_moden < TOL and gam_inv < TOL and gew_rel < TOL and gew_form < TOL and gew_rest > 0.5)

# --- 2.2 zwei Weyl-Quanten an der Schwelle E = 3
V2 = np.array([PSI[(0, s)] @ (PSI[(1, t)] @ BASE) for s in range(2) for t in range(2)]).T
gram = V2.conj().T @ V2
J2 = [kr(np.eye(DIMF), J, I9) + kr(np.eye(DIMF), I9, J) for J in JR]
CAS2 = sum((J @ J) for J in J2)
U2PI = kr(np.eye(DIMF), np.diag(SG), np.diag(SG))
RR = kr(Z2.R2(0) @ Z2.R2(1), I9, I9)
G_TOT = (RR @ U2PI).tocsr()


def schwelle(T, V):
    """Matrix von T auf span V, Rest außerhalb, invariante Unterräume und deren Casimir."""
    TV = np.array([T @ V[:, k] for k in range(V.shape[1])]).T
    M_ = V.conj().T @ TV
    rest = float(np.abs(V @ M_ - TV).max())
    w, v = np.linalg.eig(M_)
    inv = v[:, np.isclose(w, 1)]
    if inv.shape[1]:
        q_, _ = np.linalg.qr(inv)
        X = V @ q_
        cas = np.linalg.eigvalsh(X.conj().T @ np.array([CAS2 @ X[:, k] for k in range(X.shape[1])]).T)
    else:
        cas = np.array([])
    return M_, rest, inv.shape[1], [float(np.round(x, 9)) for x in sorted(cas.real)]


M_f, rest_f, dim_f, cas_f = schwelle(T_FERM, V2)
M_g, rest_g, dim_g, cas_g = schwelle(T_GEW, V2)
ew_U = sorted({float(np.round((V2[:, k].conj() @ (U2PI @ V2[:, k])).real, 12)) for k in range(4)})
ew_G = sorted({float(np.round((V2[:, k].conj() @ (G_TOT @ V2[:, k])).real, 12)) for k in range(4)})
E_V2 = sorted({float(np.round(V2[:, k].conj() @ (kr(np.eye(DIMF), ES, I9) + kr(np.eye(DIMF), I9, ES)) @ V2[:, k], 12).real)
               for k in range(4)})
komm_G = max(nrm(G_TOT @ T_FERM - T_FERM @ G_TOT), nrm(G_TOT @ T_GEW - T_GEW @ G_TOT))
ZAHLEN["schwelle_zwei"] = {"Gram der vier Zustände": np.round(gram.real, 12).tolist(), "E": E_V2,
                           "U_2pi, G = U_2pi R^2": [ew_U, ew_G],
                           "Gamma(pi): dim invariant, Casimir": [dim_f, cas_f],
                           "T_gew: dim invariant, Casimir": [dim_g, cas_g], "[G_tot, T] max": komm_G}
check("[E] Zwei Weyl-Quanten an der Schwelle: Psi^(0)_s Psi^(1)_t auf zwei skalaren Grundzuständen geben vier orthonormale "
      f"Zustände mit E = {E_V2}, U_2pi = {ew_U} und G = U_2pi R^2 = {ew_G} (behalten, [G_tot, Tausch] = {komm_G:.0e}). Der "
      f"kanonische Tausch Gamma(pi) (x) Tausch der Standardformen wirkt als -Flip; P_ident lässt genau {dim_f} Zustand mit "
      f"Casimir {cas_f}: j = 0, antisymmetrisch in der Schrittart, Spin-Statistik richtig. Der gewöhnliche Tausch lässt "
      f"{dim_g} Zustände mit Casimir {cas_g} (j = 1, z97s Sym^2). Die gemeinsame CAR wählt die richtige Kombination",
      np.allclose(gram, np.eye(4)) and rest_f < TOL and rest_g < TOL and dim_f == 1 and cas_f == [0.0] and dim_g == 3
      and cas_g == [2.0, 2.0, 2.0] and ew_U == [1.0] and ew_G == [1.0] and E_V2 == [3.0] and komm_G < TOL)

# --- 2.3 drei Weyl-Quanten an der Schwelle E = 9/2
PERMS3 = list(itertools.permutations(range(3)))


def p_perm3(perm):
    Pm = np.zeros((27, 27))
    for p in itertools.product(range(3), repeat=3):
        q = [0, 0, 0]
        for i in range(3):
            q[perm[i]] = p[i]
        Pm[q[0]*9 + q[1]*3 + q[2], p[0]*9 + p[1]*3 + p[2]] = 1.0
    return Pm


T3 = {pm: kr(Z3.tausch(list(pm)), p_perm3(pm)) for pm in PERMS3}
T3g = {pm: kr(Z3.tausch(list(pm), kanonisch=False), p_perm3(pm)) for pm in PERMS3}
komp = lambda s_, t_: tuple(s_[t_[i]] for i in range(3))
darst = max(nrm(T3[s_] @ T3[t_] - T3[komp(s_, t_)]) for s_ in PERMS3 for t_ in PERMS3)
darst_g = max(nrm(T3g[s_] @ T3g[t_] - T3g[komp(s_, t_)]) for s_ in PERMS3 for t_ in PERMS3)
BASE3 = kr(Z3.vac.reshape(-1, 1), *[np.eye(3)[:, [0]]]*3).toarray().ravel()
V3 = np.array([PSI3[(0, s)] @ (PSI3[(1, t)] @ (PSI3[(2, u)] @ BASE3)) for s in range(2) for t in range(2)
               for u in range(2)]).T
J3 = [sum(kr(np.eye(Z3.dim), *P3OPS(i, J)) for i in range(3)) for J in JRs]
CAS3 = sum(J @ J for J in J3)


def invariant3(T):
    Pi = sum(T[pm] for pm in PERMS3)/6
    M_ = V3.conj().T @ np.array([Pi @ V3[:, k] for k in range(8)]).T
    w, v = np.linalg.eigh((M_ + M_.conj().T)/2)
    X = V3 @ v[:, np.isclose(w, 1)]
    cas = np.linalg.eigvalsh(X.conj().T @ np.array([CAS3 @ X[:, k] for k in range(X.shape[1])]).T) if X.shape[1] else []
    return int(round(np.trace(M_).real)), sorted(float(np.round(x, 9)) for x in np.real(cas))


dim3_f, cas3_f = invariant3(T3)
dim3_g, cas3_g = invariant3(T3g)
ZAHLEN["schwelle_drei"] = {"Darstellungsfehler S_3 (kanonisch, gewöhnlich)": [darst, darst_g],
                           "kanonisch: dim invariant": dim3_f, "gewöhnlich: dim invariant, Casimir": [dim3_g, cas3_g]}
check("[E] Drei Nähte, drei Weyl-Quanten an der Schwelle E = 9/2 (acht Zustände Psi^(0) Psi^(1) Psi^(2) auf drei skalaren "
      f"Grundzuständen): Gamma(sigma) ist eine Darstellung von S_3 (Fehler {darst:.0e}); P_ident lässt {dim3_f} Zustände "
      f"(Pauli: drei gleiche Fermionen passen nicht in ein Dublett), der gewöhnliche Tausch {dim3_g} mit Casimir {cas3_g} "
      "(j = 3/2, bosonisch)",
      darst < TOL and darst_g < TOL and dim3_f == 0 and dim3_g == 4 and cas3_g == [3.75]*4)

# --- 2.4 Familie der gaußschen Tauschoperatoren Gamma(u), u = [[0, e^(i phi)], [e^(i theta), 0]]
V1 = np.array([PSI[(i, s)] @ BASE for i in range(2) for s in range(2)]).T
famil = []
for jt in range(4):
    for jp in range(4):
        th, ph = jt*np.pi/2, jp*np.pi/2
        Gu = (GAM_F @ sps.diags(np.exp(1j*(th*Z2.Q[:, 0] + ph*Z2.Q[:, 1])))).tocsr()
        Tu = kr(Gu, SWAP_P)
        inv_ = nrm(Tu @ Tu - sps.identity(Tu.shape[0])) < TOL
        lam = complex(V2[:, 2].conj() @ (Tu @ V2[:, 1]))              # Tu v_01 = lam v_10
        lam_ok = np.allclose(Tu @ V2[:, 1], lam*V2[:, 2]) and np.isclose(lam, -np.exp(1j*(th + ph)))
        M1 = V1.conj().T @ np.array([Tu @ V1[:, k] for k in range(4)]).T
        einzel = int(np.sum(np.isclose(np.linalg.eigvals(M1), 1)))
        famil.append((jt, jp, inv_, complex(np.round(lam, 12)), lam_ok, einzel))
invol = [(a_, b_) for a_, b_, i_, _, _, _ in famil if i_]
ZAHLEN["gauss_tausch"] = {"(theta, phi) in Einheiten pi/2 mit Gamma(u)^2 = 1": invol,
                          "Austauschfaktor je (theta, phi)": {f"{a_},{b_}": str(l_) for a_, b_, _, l_, _, _ in famil},
                          "invariante Einzel-Weyl-Zustände je (theta, phi)": {f"{a_},{b_}": e_ for a_, b_, _, _, _, e_ in famil}}
fam_ok = all(l_ok for *_, l_ok, _ in famil) and all(i_ == ((a_ + b_) % 4 == 0) for a_, b_, i_, *_ in famil) \
    and all((e_ == 2) == i_ and (e_ == 0) == (not i_) for _, _, i_, _, _, e_ in famil) \
    and all(np.isclose(l_, -1) for _, _, i_, l_, _, _ in famil if i_)
check("[E] Alle gaußschen Nahttausche Gamma(u), u = Vertauschung mit Phasen (16 Paare theta, phi in Vielfachen von pi/2): Der "
      "Austauschfaktor auf zwei Weyl-Quanten ist exakt det u = -e^(i(theta + phi)); Gamma(u)^2 = e^(i(theta + phi) Q). "
      f"Involutiv sind genau die {len(invol)} Fälle theta + phi = 0 mod 2 pi (u = D sigma_x D^+, Eichung der Nahtphasen), und "
      "genau dann überleben einzelne Weyl-Quanten unter P_ident (2 invariante Zustände, sonst 0). Jeder dieser gaußschen "
      "Tausche, der einzelne Fermionen zulässt, gibt -1: bosonische Statistik für Weyl-Quanten verlangt einen nicht gaußschen "
      "Tausch (T_gew) oder verwirft alle Zustände ungerader Ladung", fam_ok and len(invol) == 4)

# --- 2.5 Koszul-Vorzeichen gegen Spin-Statistik und GSO
ART = {"Skalar (1,1)": (lambda Z, i: Z.I, lambda: I9, 0, 0),
       "Weyl V+ S_0": (lambda Z, i: Z.U(i, "a"), lambda: SL[0], 1, 1),
       "Anti-Weyl V- S-quer_0": (lambda Z, i: Z.U(i, "a").T.tocsr(), lambda: SR[0], -1, 1),
       "nackte Ladung V+": (lambda Z, i: Z.U(i, "a"), lambda: I9, 0, 1),
       "nackter Schritt S_0": (lambda Z, i: Z.I, lambda: SL[0], 1, 0)}
iot = lambda name, i: kr(ART[name][0](Z2, i), *(ART[name][1](), I9)[::1 if i == 0 else -1])
tab, gso_ok = {}, True
for x_, y_ in itertools.product(ART, repeat=2):
    v = iot(x_, 0) @ (iot(y_, 1) @ BASE)
    w = iot(y_, 0) @ (iot(x_, 1) @ BASE)
    if np.linalg.norm(v) < 0.5:
        continue
    kap = complex(w.conj() @ (T_FERM @ v))
    ok_ = np.allclose(T_FERM @ v, kap*w)
    _, _, hx, px = ART[x_]
    _, _, hy, py = ART[y_]
    Gx, Gy = (-1)**(abs(hx) + px), (-1)**(abs(hy) + py)
    koszul = (-1)**(px*py)
    spin = (-1)**(abs(hx)*abs(hy))
    tab[f"{x_} | {y_}"] = {"Austausch": int(round(kap.real)), "Koszul": koszul, "Spin-Statistik": spin, "G": [Gx, Gy]}
    gso_ok &= ok_ and np.isclose(kap, koszul) and ((koszul == spin) or (Gx == -1 or Gy == -1))
diag_tab = {x_: (tab[f"{x_} | {x_}"]["Austausch"], tab[f"{x_} | {x_}"]["Spin-Statistik"], tab[f"{x_} | {x_}"]["G"][0])
            for x_ in ART}
behalten_ok = all(d[0] == d[1] for x_, d in diag_tab.items() if d[2] == 1)
verworfen_falsch = all(d[0] != d[1] for x_, d in diag_tab.items() if d[2] == -1)
ZAHLEN["koszul_spin"] = {"gleiche Quanten: (Austausch, Spin-Statistik, G)": diag_tab, "Tabelle": tab}
check("[E] Koszul gegen Spin-Statistik: Gamma(pi) tauscht zwei beliebige Quanten x, y (Skalar, Weyl, Anti-Weyl, nackte "
      "Ladung, nackter Schritt) mit dem Vorzeichen (-1)^(|x||y|), |x| = R^2-Parität (25 Paare). Für gleiche Quanten "
      f"(Austausch, Spin-Statistik, G) = {diag_tab}: auf G = +1 ist das Koszul-Vorzeichen genau (-1)^(2h 2h'), auf den beiden "
      "G = -1-Quanten (skalares Fermion, bosonischer Spinor) genau falsch. Die gemeinsame CAR gibt Ladungsstatistik "
      "(R^2), erst P_GSO macht daraus Spin-Statistik; G_tot vertauscht mit dem Tausch",
      gso_ok and behalten_ok and verworfen_falsch and len(tab) == 25)

# --- 2.6 Zählung: gradierte symmetrische Potenz gegen den freien 4D-Fockraum (z91)
EH = 60                                                            # halbe Energieeinheiten, E <= 30
KM = EH//2
BOS = {2*m: m*m for m in range(1, KM + 1)}                         # skalare Leiter: m^2 bei E = m
FER = {2*n + 1: 2*n*(n + 1) for n in range(1, KM)}                 # Weyl h = +-1/2: je n(n+1) bei E = n + 1/2


def erzeuger(bos, fer, gradiert=True):
    """dp[k][e]: Zahl der Zustände mit k Quanten und Energie e/2 (gradiert: Fermionen antisymmetrisch)."""
    dp = [[0]*(EH + 1) for _ in range(KM + 1)]
    dp[0][0] = 1
    stufen = [(e, d, 0) for e, d in bos.items()] + [(e, d, 1) for e, d in fer.items()]
    for e, d, ungerade in stufen:
        if e > EH:
            continue
        neu = [[0]*(EH + 1) for _ in range(KM + 1)]
        for k in range(KM + 1):
            for ee in range(EH + 1):
                v = dp[k][ee]
                if not v:
                    continue
                j = 0
                while k + j <= KM and ee + j*e <= EH:
                    fak = comb(d, j) if (ungerade and gradiert) else comb(d + j - 1, j)
                    if fak == 0:
                        break
                    neu[k + j][ee + j*e] += v*fak
                    j += 1
        dp = neu
    return dp


dp_g = erzeuger(BOS, FER)
dp_b = erzeuger(BOS, {})
dp_p = erzeuger(BOS, FER, gradiert=False)
fock = [sum(dp_g[k][e] for k in range(KM + 1)) for e in range(EH + 1)]
SNg = lambda N, e: sum(dp_g[k][e] for k in range(min(N, KM) + 1))
SNp = lambda N, e: sum(dp_p[k][e] for k in range(min(N, KM) + 1))
# freier 4D-Fockcharakter: prod (1 - q^m)^(-m^2) prod (1 + q^(n+1/2))^(2 n (n+1)), unabhängig über Euler-Produkte
frei = [1] + [0]*EH
for e, d in list(BOS.items()) + list(FER.items()):
    neu = [0]*(EH + 1)
    for ee, v in enumerate(frei):
        if v:
            j = 0
            while ee + j*e <= EH:
                neu[ee + j*e] += v*(comb(d, j) if e % 2 else comb(d + j - 1, j))
                j += 1
    frei = neu
Z91_FOCK = [1, 1, 5, 14, 40, 101, 266, 649, 1593, 3765, 8813, 20168, 45649, 101591, 223654, 486046]
bos_sektor = [sum(dp_b[k][2*E] for k in range(KM + 1)) for E in range(16)]
NS = (2, 3, 5, 10, 20, 30)
gf_ok = all(SNg(N, e) == frei[e] for N in NS for e in range(2*N + 1)) and fock == frei
plain_ueber = {N: SNp(N, 6) - SNg(N, 6) for N in NS}
# explizit: gradierte S_N-Invarianten auf V^(x)N (V bis E = 3), N = 2, 3
VB = [(0, 0)] + [(2*m, 0) for m in (1, 2, 3) for _ in range(m*m)] + [(2*n + 1, 1) for n in (1, 2) for _ in range(2*n*(n + 1))]


def invarianten(N, gradiert):
    cnt = {}
    perms = list(itertools.permutations(range(N)))
    for tup in itertools.product(range(len(VB)), repeat=N):
        e = sum(VB[t][0] for t in tup)
        if e > 6:
            continue
        for p in perms:
            if any(tup[p[i]] != tup[i] for i in range(N)):
                continue
            sg = 1
            if gradiert:
                img = [p[i] for i in range(N) if VB[tup[i]][1]]
                sg = -1 if sum(1 for a_, b_ in itertools.combinations(img, 2) if a_ > b_) % 2 else 1
            cnt[e] = cnt.get(e, 0) + sg
    return {e: Fr(v, len(perms)) for e, v in sorted(cnt.items())}


expl = {N: (invarianten(N, True), invarianten(N, False)) for N in (2, 3)}
expl_ok = all(expl[N][0].get(e, 0) == SNg(N, e) and expl[N][1].get(e, 0) == SNp(N, e) for N in (2, 3) for e in range(7))
ZAHLEN["zaehlung"] = {"frei 4D (Skalar + Weyl), E = 0, 1/2, ..., 6": frei[:13], "bosonischer Sektor E = 0..15": bos_sektor,
                      "gewöhnliches Sym^N minus gradiert bei E = 3": plain_ueber,
                      "explizit N = 2, 3 (gradiert, gewöhnlich), E = 0..3 in halben Einheiten":
                          {N: [[str(expl[N][0].get(e, 0)) for e in range(7)], [str(expl[N][1].get(e, 0)) for e in range(7)]]
                           for N in (2, 3)}}
check("[E] Zählung (ganzzahlig, E <= 30): Die gradierte symmetrische Potenz von (Vakuum + skalare Leiter m^2 + Weyl-Leitern "
      "h = +-1/2 je n(n+1)) über N Nähte, also P_ident mit Gamma(S_N), ist für E <= N exakt der freie 4D-Fockcharakter "
      "prod (1 - q^m)^(-m^2) prod (1 + q^(n+1/2))^(2n(n+1)) (N = 2, 3, 5, 10, 20, 30); explizit mit Koszul-Vorzeichen auf "
      f"V^(x)N für N = 2, 3 bestätigt. Der bosonische Sektor ist unverändert z91 ({bos_sektor[:8]} ...). Der gewöhnliche Tausch "
      f"zählt bei E = 3 genau {plain_ueber[3]} Zustände zu viel (Sym^2 statt Lambda^2 der vier Weyl-Zustände bei E = 3/2)",
      gf_ok and bos_sektor == Z91_FOCK and expl_ok and all(v == 4 for v in plain_ueber.values()))
befund("Teil 2: Die gemeinsame CAR wählt an der Schwelle j = 0 (Lambda^2), mit drei Quanten Pauli. Unter den gaußschen "
       "Tauschoperatoren (Vertauschung mit Nahtphasen) lassen genau die involutiven einzelne Fermionen zu, und sie geben alle "
       "-1. Das Vorzeichen folgt der "
       "R^2-Parität; mit der GSO-Regel G = 1 ist das genau Spin-Statistik, auf G = -1 genau falsch. Die Zählung ist der freie "
       "4D-Fockraum von Skalar und Weyl-Fermion, der bosonische Sektor bleibt z91", True)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

# =================================================================== Teil 3
say("\nTeil 3: Ist die gemeinsame CAR nativ? Hardy-Komposition, Superlokalität, Jordan-Wigner-Kette")

# --- 3.1 Hardy-Komposition: ein polarisierter Einteilchenraum über beide Nähte, Gamma(pi) aus den Transferströmen
EP = sum(Z2.cd[(0, q, n)] @ Z2.c[(1, q, n)] for q in "ab" for n in (-1, 0)).tocsr()   # E^+_0 = sum psi*_1 psi_2 (z99)
EM = EP.T.tocsr()
QT = Z2.Q.sum(1)
Xg = expm(1j*np.pi/2*(EP + EM).toarray()) @ np.diag(np.exp(-1j*np.pi/2*QT))
phase = Xg[Z2.vi, Z2.vi]
strom_tausch = float(np.abs(Xg - phase*GAM_F.toarray()).max())
H_ = sps.diags((Z2.f[:, 0] - Z2.f[:, 1]).astype(float))
su2 = max(nrm(EP @ EM - EM @ EP - H_), nrm(H_ @ EP - EP @ H_ - 2*EP))
dublett = max(max(nrm(EP @ Z2.cd[(1, q, n)] - Z2.cd[(1, q, n)] @ EP - Z2.cd[(0, q, n)]),
                  nrm(EP @ Z2.cd[(0, q, n)] - Z2.cd[(0, q, n)] @ EP)) for q in "ab" for n in (-1, 0))
EPu = sum(Z2u.cd[(0, q, n)] @ Z2u.c[(1, q, n)] for q in "ab" for n in (-1, 0)).tocsr()
EMu = EPu.T.tocsr()
su2_u = max(nrm(EPu @ EMu - EMu @ EPu - H_), nrm(H_ @ EPu - EPu @ H_ - 2*EPu))
moden_u = [Z2u.c[l] for l in Z2u.labels] + [Z2u.cd[l] for l in Z2u.labels]
bild_u = (EPu @ Z2u.cd[(1, "a", 0)] - Z2u.cd[(1, "a", 0)] @ EPu).tocsr()
dublett_u_rest = rest_linear(bild_u, moden_u)
dublett_u_form = nrm(bild_u - (Z2u.cd[(0, "a", 0)] - 2*Z2u.cd[(1, "a", 0)] @ EPu))
Xth = expm(1j*np.pi/4*(EP + EM).toarray())
Xth_u = expm(1j*np.pi/4*(EPu + EMu).toarray())
konj = sps.csr_matrix(Xth @ Z2.cd[(1, "a", 0)].toarray() @ Xth.conj().T)
konj_u = sps.csr_matrix(Xth_u @ Z2u.cd[(1, "a", 0)].toarray() @ Xth_u.conj().T)
rest_konj, rest_konj_u = rest_linear(konj, alle_moden), rest_linear(konj_u, moden_u)
konj_form = nrm(konj - (np.cos(np.pi/4)*Z2.cd[(1, "a", 0)] + 1j*np.sin(np.pi/4)*Z2.cd[(0, "a", 0)]))
HOP = (EP + EM).tocsr()
hop_f = nrm(GAM_F @ HOP - HOP @ GAM_F)
hop_g = nrm(GEW_F @ HOP - HOP @ GEW_F)
hop_g_form = nrm(GEW_F @ HOP @ GEW_F.T + sps.diags((-1.0)**QT) @ HOP)
ZAHLEN["hardy_komposition"] = {"Gamma(pi) - Phase exp(i pi/2 (E+ + E-)) exp(-i pi Q/2)": strom_tausch, "Phase": str(np.round(phase, 12)),
                               "su(2)-Nullmoden gemeinsam, ungradiert": [su2, su2_u],
                               "[E+, psi*_(1,x)] - psi*_(0,x) (gemeinsam)": dublett,
                               "ungradiert: Rest zum Modenspann, Form psi*_0 - 2 psi*_1 E+": [dublett_u_rest, dublett_u_form],
                               "Drehung exp(i pi/4 (E+ + E-)) psi*_1: Rest zum Modenspann gemeinsam, ungradiert":
                                   [rest_konj, rest_konj_u],
                               "Hüpfterm E+ + E-: [Gamma, H], [T_gew, H], T_gew H T_gew + (-1)^Q H": [hop_f, hop_g, hop_g_form]}
check("[E] Hardy-Komposition (ein polarisierter Einteilchenraum h_0 (+) h_1, eine Fermisee-Projektion; ihr Fockraum ist die "
      "gemeinsame CAR): Die Nahtvertauschung liegt in der U(2) der Nullmoden E^+ = sum_x psi*_(0,x) psi_(1,x), E^- = (E^+)^+ der "
      f"Transferströme von z99: Gamma(pi) = Phase * exp(i pi/2 (E^+ + E^-)) exp(-i pi Q/2) exakt (Fehler {strom_tausch:.0e}). "
      f"Nur in der gemeinsamen CAR wirkt der Transfer linear auf die Nahtfermionen, [E^+, psi*_(1,x)] = psi*_(0,x) (Fehler "
      f"{dublett:.0e}), und die Drehung exp(i pi/4 (E^+ + E^-)) bildet psi*_(1,x) auf cos psi*_(1,x) + i sin psi*_(0,x) ab "
      f"(Fehler {konj_form:.0e}): (psi_0, psi_1) ist ein U(2)-Dublett. Ungradiert ist [E^+, psi*_(1,x)] = psi*_(0,x) - 2 psi*_(1,x) "
      f"E^+ (Fehler {dublett_u_form:.0e}), nicht linear, und die Drehung führt aus dem Modenspann heraus (Rest {rest_konj_u:.3f} "
      f"gegen {rest_konj:.0e}). Die su(2)-Relationen der Nullmoden gelten in beiden ({su2:.0e}, {su2_u:.0e}). Ein Hüpfterm "
      f"E^+ + E^- ist unter Gamma(pi) invariant ({hop_f:.0e}), unter T_gew nicht ({hop_g:.2f}; T_gew H T_gew = -(-1)^Q H, Fehler "
      f"{hop_g_form:.0e})",
      strom_tausch < 1e-12 and abs(abs(phase) - 1) < 1e-12 and dublett < TOL and su2 < TOL and su2_u < TOL
      and konj_form < 1e-12 and rest_konj < 1e-12 and rest_konj_u > 0.1 and dublett_u_form < TOL and dublett_u_rest > 0.3
      and hop_f < TOL and hop_g > 1 and hop_g_form < TOL)

# --- 3.2 Superlokalität des Gitters Z^(2N): Kozyklus = Jordan-Wigner-Kette
KEYS = [(i, q) for i in range(2) for q in "ab"]


def nullmode(Z, alpha):
    X = Z.I
    for (i, q), a in zip(KEYS, alpha):
        if a == 1:
            X = X @ Z.U(i, q)
        elif a == -1:
            X = X @ Z.U(i, q).T
    return X.tocsr()


lok = {"gemeinsam": [0, 0, 0], "ungradiert": [0, 0, 0]}
fehl_innen = {"gemeinsam": 0, "ungradiert": 0}
for name_, Z in (("gemeinsam", Z2), ("ungradiert", Z2u)):
    nm = {al: nullmode(Z, al) for al in itertools.product((-1, 0, 1), repeat=4)}
    for al, be in itertools.product(nm, repeat=2):
        if any(abs(x_ + y_) > 1 for x_, y_ in zip(al, be)):
            continue
        v1, v2 = nm[al] @ (nm[be] @ Z.vac), nm[be] @ (nm[al] @ Z.vac)
        if np.linalg.norm(v1) < 0.5 or np.linalg.norm(v2) < 0.5:
            continue
        lam = float(v2 @ v1)
        soll = (-1)**(sum(x_*y_ for x_, y_ in zip(al, be)) + sum(x_*x_ for x_ in al)*sum(y_*y_ for y_ in be))
        lok[name_][0] += 1
        lok[name_][1] += int(np.isclose(lam, soll))
        jw = (-1)**sum(al[k]*be[l] for k in range(4) for l in range(4) if k != l)
        lok[name_][2] += int(np.isclose(lam, jw))
        if not np.isclose(lam, soll) and all((x_ == 0 and y_ == 0) or KEYS[k][0] == 0 for k, (x_, y_) in enumerate(zip(al, be))):
            fehl_innen[name_] += 1
ZAHLEN["superlokalitaet"] = {"(Paare, Superlokalität erfüllt, = JW-Kozyklus)": lok, "Fehler innerhalb einer Naht": fehl_innen}
check("[E] Superlokalität: Für das Gitter Z^4 der Quarterladungen zweier Nähte verlangt eine Gittervertex-Superalgebra "
      "e^alpha e^beta = (-1)^(alpha.beta + |alpha|^2 |beta|^2) e^beta e^alpha für die Nullmodenverschiebungen. In der gemeinsamen "
      f"CAR gilt das für alle {lok['gemeinsam'][0]} Paare im Fenster ({lok['gemeinsam'][1]} erfüllt), und der Faktor ist der "
      "Jordan-Wigner-Kozyklus (-1)^(sum_(k != l) alpha_k beta_l) (Klein-Faktoren = Kette). Ungradiert sind nur "
      f"{lok['ungradiert'][1]} von {lok['ungradiert'][0]} erfüllt; die Fehler liegen alle über Nähte hinweg (innerhalb einer "
      f"Naht {fehl_innen['ungradiert']}): Die Quarter a, b des Originals sind schon gradiert komponiert, die gemeinsame CAR "
      "setzt dieselbe Regel über die Nähte fort",
      lok["gemeinsam"][1] == lok["gemeinsam"][0] == lok["gemeinsam"][2] and lok["gemeinsam"][0] > 100
      and lok["ungradiert"][1] < lok["ungradiert"][0] and fehl_innen["ungradiert"] == 0 and fehl_innen["gemeinsam"] == 0)

# --- 3.3 Ist die Jordan-Wigner-Kette beobachtbar?
ZB = Naehte(2, reihenfolge=(1, 0))
D12 = sps.diags((-1.0)**(Z2.f[:, 0]*Z2.f[:, 1])).tocsr()
ord_fehler = max(nrm(D12 @ Z2.c[l] @ D12 - ZB.c[l]) for l in Z2.labels)
d_gleich_w = nrm(D12 - W12)
xs ={(i, q, n): float({"a": Fr(1, 4), "b": Fr(3, 4)}[q] + n) for (i, q, n) in Z2.labels}
fk = {l: 1/(1 + np.exp(4*np.pi*x)) for l, x in xs.items()}                 # KMS-Besetzung bei beta = 4 pi (z58)
rho = np.ones(Z2.dim)
for l in Z2.labels:
    o = Z2.occ[:, Z2.bit[l]]
    rho *= np.where(o == 1, fk[l], 1 - fk[l])
korr_diff, n_korr = 0.0, 0
for _ in range(300):
    ls = [Z2.labels[k] for k in rng.choice(len(Z2.labels), size=int(rng.integers(2, 4)), replace=False)]
    wort = [(l, d_) for l in ls for d_ in (0, 1)]                  # jede Mode einmal erzeugt, einmal vernichtet
    wort = [wort[k] for k in rng.permutation(len(wort))]
    wA, wB = Z2.I, ZB.I
    for l, d_ in wort:
        wA = wA @ (Z2.cd[l] if d_ else Z2.c[l])
        wB = wB @ (ZB.cd[l] if d_ else ZB.c[l])
    a_, b_ = wA.diagonal() @ rho, wB.diagonal() @ rho
    korr_diff = max(korr_diff, abs(a_ - b_))
    n_korr += int(abs(a_) > 1e-6)
anti_B = nrm(ZB.U(0) @ ZB.U(1) + ZB.U(1) @ ZB.U(0))
Q_tot = Z2.Q.sum(1)
sektor = bool(np.all((-1.0)**(Z2.Q[:, 0]*Z2.Q[:, 1]) == (-1.0)**((Q_tot + 1)*Z2.Q[:, 0])))
ord3 = 0.0
for pm in PERMS3:
    Zs = Naehte(3, ("a",), reihenfolge=pm)
    dz = np.zeros(Z3.dim)
    for i, j in itertools.combinations(range(3), 2):
        if pm.index(i) > pm.index(j):
            dz += Z3.f[:, i]*Z3.f[:, j]
    Dz = sps.diags((-1.0)**dz).tocsr()
    ord3 = max(ord3, max(nrm(Dz @ Z3.c[l] @ Dz - Zs.c[l]) for l in Z3.labels))
# Kommutant und Superkommutant der Naht 0 (zwei Nähte, Quarter a, 16 Zustände)
Z2a = Naehte(2, ("a",))
d16 = Z2a.dim
g0 = majorana(Z2a, [(0, "a", -1), (0, "a", 0)])
g1 = majorana(Z2a, [(1, "a", -1), (1, "a", 0)])
A1 = np.array([X.toarray().ravel() for X in monome(g1)]).T
Pi_ = np.diag((-1.0)**Z2a.f.sum(1))
Id = np.eye(d16)
Me = (np.eye(d16*d16) + np.kron(Pi_, Pi_))/2
Mo = (np.eye(d16*d16) - np.kron(Pi_, Pi_))/2
L_komm = np.vstack([np.kron(Id, G.toarray().T) - np.kron(G.toarray(), Id) for G in g0])
L_super = np.vstack([(np.kron(Id, G.toarray().T) - np.kron(G.toarray(), Id)) @ Me
                     + (np.kron(Id, G.toarray().T) + np.kron(G.toarray(), Id)) @ Mo for G in g0])


def kern(L):
    _, s, vh = np.linalg.svd(L)
    return vh[np.sum(s > 1e-9):].conj().T


K_komm, K_super = kern(L_komm), kern(L_super)
rang = lambda X: int(np.linalg.matrix_rank(X, tol=1e-8))
dim_komm, dim_super = K_komm.shape[1], K_super.shape[1]
super_gleich = rang(np.hstack([K_super, A1])) == rang(A1) == dim_super
schnitt = rang(K_komm) + rang(A1) - rang(np.hstack([K_komm, A1]))
F0 = np.diag((-1.0)**Z2a.f[:, 0])
A1_tw = np.array([(X.toarray() if bin(S).count("1") % 2 == 0 else F0 @ X.toarray()).ravel()
                  for S, X in enumerate(monome(g1))]).T
komm_gleich = rang(np.hstack([K_komm, A1_tw])) == rang(A1_tw) == dim_komm
# gerade Daten: gemeinsame und ungradierte Komposition haben dieselbe gerade Algebra
Z2au = Naehte(2, ("a",), ungradiert=True)
ger = lambda Z: np.array([X.toarray().ravel() for S, X in enumerate(monome(majorana(Z, Z.labels)))
                          if bin(S).count("1") % 2 == 0]).T
Gg, Gu = ger(Z2a), ger(Z2au)
gerade_gleich = rang(Gg) == rang(Gu) == rang(np.hstack([Gg, Gu])) == d16*d16//2
ZAHLEN["kette"] = {"Ordnung (1,0) gegen (0,1): D c D - c', D = (-1)^(F_0 F_1)": ord_fehler, "D - W_12": d_gleich_w,
                   "KMS-Korrelationen (300 Wörter, davon ungleich null)": [korr_diff, n_korr],
                   "{U_0, U_1} in Ordnung (1,0)": anti_B, "drei Nähte, alle Ordnungen: Fehler": ord3,
                   "W_12 = R_0^(2(Q+1)) auf Sektoren": sektor,
                   "dim Kommutant, dim Superkommutant, dim A_1, dim (Kommutant n A_1)": [dim_komm, dim_super, rang(A1), schnitt],
                   "gerade Algebren gleich (Rang)": [rang(Gg), rang(Gu), rang(np.hstack([Gg, Gu]))]}
check("[E] Die Jordan-Wigner-Kette ist nicht beobachtbar: Die Ordnungen (0,1) und (1,0) hängen über D = (-1)^(F_0 F_1) = W_12 "
      f"zusammen (D c D = c' für alle Moden, Fehler {ord_fehler:.0e}; D - W_12 = {d_gleich_w:.0e}), drei Nähte über "
      f"Produkte (-1)^(F_i F_j) für alle sechs Ordnungen ({ord3:.0e}). D ist diagonal, fixiert den KMS-Zustand, und auf jedem "
      f"Sektor fester Gesamtladung ist W_12 = R_0^(2(Q+1)) ({sektor}, z91): 300 KMS-Korrelationen von Modenwörtern "
      f"({n_korr} ungleich null) sind in beiden Ordnungen gleich (max {korr_diff:.0e}), die gekleideten Vertizes "
      f"antikommutieren auch in (1,0) ({anti_B:.0e}). Lokal ist die Kette als verdrehte Dualität: Der Superkommutant der "
      f"Naht 0 ist genau die Modenalgebra der Naht 1 (dim {dim_super}); der gewöhnliche Kommutant (dim {dim_komm}) ist "
      f"A_1^gerade + (-1)^(F_0) A_1^ungerade und teilt mit A_1 nur den geraden Teil (dim {schnitt})",
      ord_fehler < TOL and d_gleich_w < TOL and ord3 < TOL and sektor and korr_diff < 1e-14 and n_korr > 20
      and anti_B < TOL and dim_super == 16 and super_gleich and dim_komm == 16 and komm_gleich and schnitt == 8)
check("[E negativ] Gerade Daten wählen die Komposition nicht: Gemeinsame und ungradierte Komposition haben dieselbe gerade "
      f"Algebra (Ränge {rang(Gg)}, {rang(Gu)}, zusammen {rang(np.hstack([Gg, Gu]))} = 256/2), denselben KMS-Zustand und "
      "dieselben su(2)-Relationen der Transferströme; sie unterscheiden sich nur darin, welche ungeraden Operatoren die "
      "Fermionen der zweiten Naht sind (Klein-Faktor (-1)^(F_0) = R_0^2). Aus einer einzelnen Naht oder aus geraden Daten "
      "folgt die gemeinsame CAR also nicht; sie folgt aus der Einteilchenkomposition (Dublett, Gamma(pi) in der U(2) der "
      "Transferströme, Hüpfterm) oder aus der Superlokalität des Gitters",
      gerade_gleich and su2_u < TOL and su2 < TOL)
befund("Teil 3: Die gemeinsame CAR ist die Hardy-Komposition der Nähte, ein polarisierter Einteilchenraum über alle Nähte. "
       "Dort liegt der kanonische Tausch in der U(2) der Transferströme von z99, die Nahtfermionen bilden ein Dublett, ein "
       "Hüpfterm ist nur unter Gamma(pi) symmetrisch, und das Gitter Z^(2N) ist superlokal mit der Kette als Kozyklus. Die "
       "Kette selbst ist eine Ordnungseichung ohne beobachtbare Folge. Erzwungen ist die gemeinsame CAR aber nicht: Die "
       "ungradierte Komposition hat dieselbe gerade Algebra und denselben Zustand", True)

# =================================================================== Bilanz
say("\nBilanz und Urteil H7")
STATUS = {
    "Austausch gekleideter Vertizes in der gemeinsamen CAR": ("E", "Psi^(i) Psi^(j) = -Psi^(j) Psi^(i), zwei und drei Nähte; "
                                                                 "Paare und radiale Heber gerade, vertauschen"),
    "Struktur": ("E", "Super-Tensorprodukt der Nahtalgebren; Moden-Einbettung = Jordan-Wigner-Kette; ungradiert = Klein-"
                      "Transformierte"),
    "Schwelle E = 3": ("E", "Gamma(pi): ein Zustand, Casimir 0 (j = 0); gewöhnlicher Tausch: drei, Casimir 2 (j = 1, z97)"),
    "drei Weyl-Quanten an der Schwelle": ("E", "gradiert 0 Zustände (Pauli), gewöhnlich 4 mit j = 3/2"),
    "gaußsche Tausche": ("E", "Austauschfaktor det u; involutiv <=> einzelne Fermionen überleben <=> -1"),
    "Spin-Statistik": ("E/C", "Koszul = R^2-Statistik; = Spin-Statistik genau auf G = +1 (mit P_GSO)"),
    "Zählung": ("E", "gradiertes Sym^N = freier 4D-Fockraum (Skalar + Weyl) bis E = N; bosonischer Sektor = z91"),
    "Hardy-Komposition": ("C", "eine Polarisation über alle Nähte = gemeinsame CAR; Gamma(pi) in der U(2) der z99-Ströme"),
    "Superlokalität": ("E", "Z^(2N)-Kozyklus = Kette; ungradiert verletzt über Nähte hinweg"),
    "Jordan-Wigner-Kette": ("E", "nicht beobachtbar (Ordnungseichung W_12, Uhrpotenz auf Sektoren); verdrehte Dualität"),
    "gerade Daten": ("E negativ", "wählen die Komposition nicht (gleiche gerade Algebra und gleicher Zustand)"),
    "H7 gesamt": ("bedingt", "gradierte Komposition folgt aus der gemeinsamen CAR (Prämisse P_CAR = Hardy-Komposition der "
                             "Nähte, dieselbe Regel wie für die Quarter des Originals, Rahmen von z91 und z99)"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
ZAHLEN["H7_urteil_kurz"] = {"gradierte Komposition": "in der gemeinsamen CAR automatisch",
                            "gemeinsame CAR": "bedingt (P_CAR), nicht aus geraden Daten",
                            "Statistikprämissen danach": "P_GSO (Spin = Uhr), P_ident, P_CAR; keine eigene gradierte Komposition"}
urteil = ("H7 ist bedingt bestätigt. In einer gemeinsamen CAR ist die gradierte Komposition keine eigene Prämisse, sie folgt "
          "exakt. Die Moden-Einbettung jeder Nahtalgebra ist die Jordan-Wigner-Kette, und die Gesamtalgebra ist das "
          "Super-Tensorprodukt der Nahtalgebren (zwei und drei Nähte). Gekleidete Vertizes V^+ (x) S_s verschiedener Nähte "
          "antikommutieren, weil der Ladungsschritt ungerade ist; Paaranregungen und radiale Heber sind gerade und vertauschen. "
          "An der Schwelle E = 3 lässt der kanonische Tausch Gamma(pi) mit P_ident genau einen Zustand, j = 0 mit Casimir 0. "
          "z97s j = 1 kommt nur vom gewöhnlichen Tausch T_gew = Gamma(pi) W_12, und der ist nicht gaußsch. Drei Weyl-Quanten "
          "an der Schwelle sind verboten (Pauli). Unter den gaußschen Nahttauschen (Vertauschung mit Nahtphasen) sind genau "
          "die involutiven die, unter denen einzelne Fermionen überleben, und alle geben -1. Das Austauschvorzeichen folgt der "
          "R^2-Parität; mit der "
          "GSO-Regel G = 1 ist es genau Spin-Statistik, auf den verworfenen Quanten genau falsch. Die gradierte symmetrische "
          "Potenz ist bis E = N der freie 4D-Fockraum von Skalar und Weyl-Fermion; der bosonische Sektor ist unverändert z91. "
          "Die Kette ist nicht beobachtbar: Ordnungen unterscheiden sich um W_12, eine Uhrpotenz auf Ladungssektoren, und alle "
          "Korrelationen sind gleich. Nativ ist die gemeinsame CAR nur bedingt. Sie ist die Hardy-Komposition, ein "
          "polarisierter Einteilchenraum über alle Nähte. Dort liegt Gamma(pi) in der U(2) der Transferströme von z99, die "
          "Nahtfermionen bilden ein Dublett, und das Gitter Z^(2N) ist superlokal mit der Kette als Kozyklus; die Quarter des "
          "Originals sind schon so komponiert. Gerade Daten erzwingen sie nicht, denn die ungradierte Komposition hat dieselbe "
          "gerade Algebra und denselben KMS-Zustand. Folge für die Szenario-Suche: Der Zusatz 'über Nähte braucht es eine "
          "gradierte Komposition' bei P_GSO entfällt zugunsten der benannten Prämisse P_CAR, die z91 und z99 schon benutzen. "
          "Spin-Statistik über Nähte ist dann P_CAR mit P_GSO und P_ident; Chiralität bleibt unberührt (Index 0, z97). Nichts "
          "ist befördert; 3+1-Aussagen bleiben Kinematik der Einpaarschicht.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z106_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
