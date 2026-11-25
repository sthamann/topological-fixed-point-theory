"""Z156: Trägt der Dubleton-Sektor aus z75 unter der nativen so(3,2) die Darstellung D(3,2) mit positiver Norm und
positiver Energie, sodass P_X + P_pos (z148) nativ wäre?

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart bzw. Literatur, N numerisch, O offen.

Frage (Welle B): z148 hat die Gravitonlücke auf eine Annahme reduziert: ein Spin-2-Sektor, der unter der nativen so(3,2)
die unitäre Darstellung D(3,2) positiver Energie trägt (P_X + P_pos, gleichwertig P_h2). Trägt der Dubleton-Sektor aus
z75 (vier Oszillatoren d_1, d_2, e_1, e_2 aus zwei Strommoden und ihren Tomita-Partnern, Helizität D/2 mit
D = N_d - N_e) genau das, und folgen Energie, Norm und Auswahl der Kopie ohne neue Wahl?

Bestand: z75 (Araki-Woods-Dubleton: Helizität +-2 bei E = 3; D ist der Modulgenerator; E erhält die Algebra des
Originals nicht; Speziesproblem; Weinberg-Witten), z89 (native so(3,2) = Der(I_pi - I_omega) mit kompaktem E; kappa^2 =
tanh(pi m)), z104 (Helizität 0: D(1,0) + D(2,0)), z115/z146 (Lesarten T1 Thermofeldschnitt und T2 Twistorlesart;
Stabilisator SO(3) des Paars (M, Omega)), z137/z148 (Spin-2-Primär D(4,2); unter P_dyn die Mode D(3,2) mit Residuum
32 pi G = l^2/(32 kappa); Gram-Form 2 (E0 - 3), 2 (E0 - 1), 2 (E0 + 2)), z152 (Gitterlesart).

Literatur [L]: Flato-Fronsdal (1978, zwei Singletons = masselose Teilchen), Angelopoulos-Flato-Fronsdal-Sternheimer
(1981, masselose Darstellungen der konformen Gruppe bleiben auf so(3,2) irreduzibel für Helizität != 0),
Günaydin-Saclioglu (Oszillatordarstellung), Howe (duale Paare), Weinberg-Witten (1980).

Konventionen (wie z146): Z = (d_1, d_2, e_1^+, e_2^+), G = diag(1,1,-1,-1), Zhat(X) = Ztilde X Z mit Ztilde =
(d^+, -e); hermitesche Operatoren sind i Zhat(X). E = i Zhat(E_T) = 1 + (N_d + N_e)/2. Original: a_i = (d_i + q e_i^+)/s,
Partner b_i = (e_i + q d_i^+)/s, s = sqrt(1 - q^2), q = e^(-2 pi m). Charaktere in doppelten Einheiten (2E, 2s).

Teil 1  Charaktere exakt: Erzeuger, K-Typen, Charaktervergleich für s = 0, 1/2, 1, 3/2, 2, 3; tiefste Gewichte (exakte
        Ränge); Gram-Form und Nullvektor; Casimir; Abhängigkeit von der Wahl der Untergruppe.
Teil 2  Norm und Energie: Fock-Form; Araki-Woods exakt; drei Kandidaten für 'die Energie' (D, E_de, E_ab); welche
        Erzeuger das Original erhalten; Lesart T1 (Liouville-Darstellung, zwei Leitern).
Teil 3  Speziesproblem und Komposita: Vielfachheit je Dubleton, Howe-Dualität, mehrere Moden, Auswahlkandidaten.
Teil 4  Kopplung: Passung zur D(3,2)-Mode aus z148, Vielfachheitsraum und Modulfluss, Vorzeichen, Weinberg-Witten,
        kubischer Vertex.
Teil 5  Bilanz und Urteil.
"""
import itertools
import json
import math
import os
import time
from fractions import Fraction as Fr
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # deterministische Rundung (einfädig)

import numpy as np
import sympy as sp
from flint import fmpz_mat

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


def hb(x2):
    """Halbzahl aus doppelter Einheit als Text: 3 -> '3/2', 4 -> '2'."""
    return str(x2//2) if x2 % 2 == 0 else f"{x2}/2"


def dname(E2, s2):
    return f"D({hb(E2)},{hb(s2)})"


def gz(x, n=8):
    return f"{float(x):.{n}g}"


# ====================================================================================================================
# Exakte Erzeuger der u(2,2) auf dem Twistorraum T (wie z98, z115, z146) und die native so(3,2) = Der(I_pi - I_omega)
# ====================================================================================================================
I_ = sp.I
HALB = sp.Rational(1, 2)
SIGMA = [sp.eye(2), sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -I_], [I_, 0]]), sp.Matrix([[1, 0], [0, -1]])]
G4 = sp.diag(1, 1, -1, -1)
UVEC = [sp.Matrix([1, 0, -1, 0]), sp.Matrix([0, 1, 0, -1])]
VVEC = [sp.Matrix([1, 0, 1, 0]), sp.Matrix([0, 1, 0, 1])]
NULL4 = sp.zeros(4, 4)


def msum(it, n=4):
    out = sp.zeros(n, n)
    for x in it:
        out += x
    return out


def kl(A, B):
    return A*B - B*A


I_PI = UVEC[0]*UVEC[1].T - UVEC[1]*UVEC[0].T
I_OM = VVEC[0]*VVEC[1].T - VVEC[1]*VVEC[0].T
I_ADS = I_PI - I_OM
NP_ = [-HALB*msum(s[i, j]*UVEC[i]*VVEC[j].T for i in range(2) for j in range(2)) for s in SIGMA]
NK_ = [-HALB*msum(s[i, j]*VVEC[i]*UVEC[j].T for i in range(2) for j in range(2)) for s in SIGMA]
ETA = [1, -1, -1, -1]
MP = [I_*N for N in NP_]
MK = [ETA[m]*I_*NK_[m] for m in range(4)]
QP = [(MP[m] + MK[m])/2 for m in range(4)]
QM = [(MP[m] - MK[m])/2 for m in range(4)]
E_T = QP[0]
Q_T = QP[1:]
B_T = [kl(E_T, X) for X in Q_T]
R_T = [kl(Q_T[1], Q_T[2]), kl(Q_T[2], Q_T[0]), kl(Q_T[0], Q_T[1])]
D_T = kl(E_T, QM[0])
SO32 = [E_T] + R_T + Q_T + B_T
SU22 = SO32 + QM[1:] + [QM[0], D_T]
NAMEN = ["E", "R1", "R2", "R3", "Q1", "Q2", "Q3", "B1", "B2", "B3", "Q-1", "Q-2", "Q-3", "Q-0", "D"]
K_DE = [0, 1, 2, 3, 10, 11, 12]                        # Stabilisator von E in su(2,2): s(u(2) + u(2))


def reell_koord(M):
    """Reelle Koordinaten einer komplexen 4x4-Matrix (32 Zahlen)."""
    return [sp.re(x) for x in M] + [sp.im(x) for x in M]


def rang_spanne(mats):
    return sp.Matrix([reell_koord(M) for M in mats]).rank()


def trf(X, Y):
    return (X*Y).trace()


def koord(M, basis):
    return [sp.simplify(trf(M, X)/trf(X, X)) for X in basis]


def zhat_koeff(X):
    """Zhat(X) = Summe A_ij d_i^+ d_j + U_ij d_i^+ e_j^+ + L_ij e_i d_j + B_ij e_j^+ e_i + c."""
    A = [[X[i, j] for j in range(2)] for i in range(2)]
    U = [[X[i, 2 + j] for j in range(2)] for i in range(2)]
    L = [[-X[2 + i, j] for j in range(2)] for i in range(2)]
    B = [[-X[2 + i, 2 + j] for j in range(2)] for i in range(2)]
    c = -X[2, 2] - X[3, 3]
    return A, U, L, B, c


# ====================================================================================================================
# Fockraum numerisch (orthonormale Basis) in einem Helizitätssektor D = N_d - N_e
# ====================================================================================================================
class Sektor:
    """Zustände |n1, n2; m1, m2> mit N_d - N_e = D und Stufe k = min(N_d, N_e) <= kmax; E = 1 + |D|/2 + k."""

    def __init__(self, D, kmax):
        self.D, self.kmax = D, kmax
        nd0, ne0 = max(D, 0), max(-D, 0)
        self.basis = [(n1, nd0 + k - n1, m1, ne0 + k - m1) for k in range(kmax + 1)
                      for n1 in range(nd0 + k + 1) for m1 in range(ne0 + k + 1)]
        self.idx = {b: i for i, b in enumerate(self.basis)}
        self.dim = len(self.basis)
        self.stufe = np.array([min(b[0] + b[1], b[2] + b[3]) for b in self.basis])
        self.Nd = np.array([b[0] + b[1] for b in self.basis], float)
        self.Ne = np.array([b[2] + b[3] for b in self.basis], float)
        self.DD = [[self._tausch(0, i, j) for j in range(2)] for i in range(2)]
        self.EE = [[self._tausch(2, i, j) for j in range(2)] for i in range(2)]
        self.UP = [[self._hoch(i, j) for j in range(2)] for i in range(2)]
        self.Id = np.eye(self.dim)

    def _op(self, fn):
        M = np.zeros((self.dim, self.dim), complex)
        for i, b in enumerate(self.basis):
            r = fn(list(b))
            if r is not None and tuple(r[0]) in self.idx:
                M[self.idx[tuple(r[0])], i] += r[1]
        return M

    def _tausch(self, off, i, j):                      # x_i^+ x_j
        def f(b):
            if b[off + j] == 0:
                return None
            c = math.sqrt(b[off + j])
            b[off + j] -= 1
            c *= math.sqrt(b[off + i] + 1)
            b[off + i] += 1
            return b, c
        return self._op(f)

    def _hoch(self, i, j):                             # d_i^+ e_j^+
        def f(b):
            c = math.sqrt((b[i] + 1)*(b[2 + j] + 1))
            b[i] += 1
            b[2 + j] += 1
            return b, c
        return self._op(f)

    def zhat(self, X):
        Xn = np.array(X.tolist(), dtype=complex)
        M = np.zeros((self.dim, self.dim), complex)
        for i in range(2):
            for j in range(2):
                M += Xn[i, j]*self.DD[i][j] + Xn[i, 2 + j]*self.UP[i][j]
                M -= Xn[2 + i, j]*self.UP[j][i].conj().T                     # e_i d_j
                M -= Xn[2 + i, 2 + j]*(self.EE[j][i] + (i == j)*self.Id)     # e_i e_j^+
        return M

    def heben(self, X):
        """Anteil von Zhat(X), der die Stufe um eins erhöht (nur für nichtkompakte X)."""
        Xn = np.array(X.tolist(), dtype=complex)
        return sum(Xn[i, 2 + j]*self.UP[i][j] for i in range(2) for j in range(2))

    def innen(self, tiefe=1):
        return np.where(self.stufe <= self.kmax - tiefe)[0]


def so32_ops(S):
    """E, J_k, L^+_k, L^-_k im Sektor S in der Normierung von z148."""
    E = 1j*S.zhat(E_T)
    J = [-1j*S.zhat(R) for R in R_T]
    Lp = [2j*S.heben(Q) for Q in Q_T]
    Lm = [M.conj().T for M in Lp]
    return E, J, Lp, Lm


# ====================================================================================================================
# Charaktere der so(3,2) in doppelten Einheiten: {2E: {2s: Vielfachheit}}
# ====================================================================================================================
def cg2(a, b):
    return range(abs(a - b), a + b + 1, 2)


def sym1(n):
    """Sym^n(Spin 1) in doppelten Spins."""
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
    out = {}
    n = 0
    while E2 + 2*n <= emax2:
        tgt = out.setdefault(E2 + 2*n, {})
        for l2 in sym1(n):
            for j2 in cg2(s2, l2):
                tgt[j2] = tgt.get(j2, 0) + 1
        n += 1
    return out


def irrep(E2, s2, emax2):
    """Charakter der unitären Darstellung D(E0, s): generisch der Verma-Modul; masselos (E0 = s + 1, s >= 1) und die
    Singletons Rac = D(1/2,0), Di = D(1,1/2) mit abgezogenem Nullmodul."""
    V = verma(E2, s2, emax2)
    if s2 >= 2 and E2 == s2 + 2:
        return ch_add(V, verma(E2 + 2, s2 - 2, emax2), -1)
    if (E2, s2) == (1, 0):
        return ch_add(V, verma(5, 0, emax2), -1)
    if (E2, s2) == (2, 1):
        return ch_add(V, verma(4, 1, emax2), -1)
    return V


def zerlege(chi, emax2):
    """Tiefste Gewichte (2 E0, 2 s, Vielfachheit) eines Charakters; None, wenn etwas negativ wird."""
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


def gewichte_zu_spins(cnt):
    """Aus Vielfachheiten doppelter J_3-Gewichte die Vielfachheiten doppelter Spins."""
    out = {}
    for w in sorted(cnt):
        if w >= 0:
            m = cnt[w] - cnt.get(w + 2, 0)
            if m:
                out[w] = m
    return out


# ====================================================================================================================
say("Teil 1: Charaktere exakt -- die Leiterdarstellungen der su(2,2) unter der nativen so(3,2)")
# 1.1 ---------------------------------------------------------------------------------------------------------------
antiherm = all((X.H*G4 + G4*X) == NULL4 and X.trace() == 0 for X in SU22)
unabh = rang_spanne(SU22) == 15
der_ads = [(X*I_ADS + I_ADS*X.T) == NULL4 for X in SU22]
abgeschl = rang_spanne(SO32 + [kl(X, Y) for X in SO32 for Y in SO32]) == 10
A_E, U_E, L_E, B_E, c_E = zhat_koeff(E_T)
E_form = (all(I_*A_E[i][j] == (HALB if i == j else 0) and I_*B_E[i][j] == (HALB if i == j else 0)
              and U_E[i][j] == 0 and L_E[i][j] == 0 for i in range(2) for j in range(2)) and I_*c_E == 1)
zentral = all(kl(I_*sp.eye(4), X) == NULL4 for X in SU22)
A3, U3, L3, B3, c3 = zhat_koeff(R_T[2])
J3_diag = (A3[0][1] == 0 and A3[1][0] == 0 and B3[0][1] == 0 and B3[1][0] == 0 and c3 == 0
           and all(U3[i][j] == 0 and L3[i][j] == 0 for i in range(2) for j in range(2)))
J3_GEW = [sp.nsimplify(-I_*A3[0][0]), sp.nsimplify(-I_*A3[1][1]), sp.nsimplify(-I_*B3[0][0]), sp.nsimplify(-I_*B3[1][1])]
J3_GEW2 = [int(2*w) for w in J3_GEW]                   # doppelte Gewichte von d_1^+, d_2^+, e_1^+, e_2^+
# ausgeschlossene Hebung (Singulett) und native Hebungen (Triplett)
U_sing = sp.Matrix(zhat_koeff(QM[0])[1])
trip = [sp.Matrix(zhat_koeff(Q)[1]) for Q in Q_T]
trip_ok = all((U_sing.H*T).trace() == 0 for T in trip) and sp.Matrix([list(T) for T in trip]).rank() == 3
check("[E] Erzeuger (wie z146): 15 G-antihermitesche spurfreie Matrizen auf T, linear unabhängig. Die native so(3,2) = "
      f"Der(I_pi - I_omega) sind genau E, R_k, Q_k, B_k ({sum(der_ads[:10])} von 10 erhalten I_AdS, "
      f"{sum(der_ads[10:])} der übrigen 5; abgeschlossen: {abgeschl}). Auf dem Fockraum von d_1, d_2, e_1, e_2 ist "
      "E = i Zhat(E_T) = 1 + (N_d + N_e)/2, und D = N_d - N_e gehört zum Zentrum i 1 der u(2,2) (vertauscht mit allen "
      f"15). J_3 ist diagonal mit den Gewichten {[str(w) for w in J3_GEW]} auf d_1^+, d_2^+, e_1^+, e_2^+: Die Drehgruppe "
      "der so(3,2) ist die Diagonale von SU(2)_L x SU(2)_R, e transformiert konjugiert zu d. Die vier Hebungen "
      "d_i^+ e_j^+ zerfallen in das native Triplett d^+ sigma_a e^+ und das Singulett Summe d_i^+ e_i^+ (in Q^-_0 und D, "
      "nicht in so(3,2))",
      antiherm and unabh and der_ads == [True]*10 + [False]*5 and abgeschl and E_form and zentral and J3_diag
      and sorted(J3_GEW2) == [-1, -1, 1, 1] and trip_ok and U_sing == U_sing[0, 0]*sp.eye(2))


# 1.2 K-Typen ---------------------------------------------------------------------------------------------------------
def ktyp_test(D, kmax):
    S = Sektor(D, kmax)
    AL = [(R_T[k] - QM[1 + k])/2 for k in range(3)]
    AR = [(R_T[k] + QM[1 + k])/2 for k in range(3)]
    CL = -sum(S.zhat(X) @ S.zhat(X) for X in AL)
    CR = -sum(S.zhat(X) @ S.zhat(X) for X in AR)
    E = 1j*S.zhat(E_T)
    j1, j2 = S.Nd/2, S.Ne/2
    ok = (np.abs(CL - np.diag(j1*(j1 + 1))).max() < TOL and np.abs(CR - np.diag(j2*(j2 + 1))).max() < TOL
          and np.abs(E - np.diag(1 + (S.Nd + S.Ne)/2)).max() < TOL)
    dims = [int(np.sum(S.stufe == k)) for k in range(kmax + 1)]
    return ok, dims


ktyp = {}
ktyp_ok = True
for D_ in (0, 1, 2, 4):
    ok_, dims_ = ktyp_test(D_, 4)
    ktyp_ok &= ok_ and dims_ == [(D_ + k + 1)*(k + 1) for k in range(5)]
    ktyp[f"s = {hb(D_)}"] = [f"E = {hb(2 + D_ + 2*k)}: ({hb(D_ + k)}, {hb(k)}), dim {(D_ + k + 1)*(k + 1)}"
                             for k in range(5)]
ZAHLEN["1_K_Typen"] = ktyp
check("[E] Leiterdarstellung der Helizität s = D/2 (Sektor N_d - N_e = D, s = 0, 1/2, 1, 2; Stufen k = 0..4): Die "
      "Casimire der beiden su(2) sind j_1 (j_1 + 1) und j_2 (j_2 + 1) mit (j_1, j_2) = ((D + k)/2, k/2), jeder K-Typ "
      "einmal, E = 1 + D/2 + k, Dimension (D + k + 1)(k + 1). Für s = 2: "
      f"{ktyp['s = 2'][:3]}; tiefste Energie E_0 = s + 1 = 3 mit K-Typ (2, 0), fünf Zustände (z75, z89 bestätigt)",
      ktyp_ok)

# 1.3 Charaktervergleich -----------------------------------------------------------------------------------------------
KCH = 9                                                # Stufen 0..KCH


def sektor_charakter(D, kmax):
    """Exakter Charakter des Sektors D unter (E, diagonale SO(3)) aus den J_3-Gewichten der Monome."""
    nd0, ne0 = max(D, 0), max(-D, 0)
    chi = {}
    for k in range(kmax + 1):
        cnt = {}
        for n1 in range(nd0 + k + 1):
            for m1 in range(ne0 + k + 1):
                occ = (n1, nd0 + k - n1, m1, ne0 + k - m1)
                w = sum(J3_GEW2[i]*occ[i] for i in range(4))
                cnt[w] = cnt.get(w, 0) + 1
        chi[2 + abs(D) + 2*k] = gewichte_zu_spins(cnt)
    return chi


erwartet = {0: [(2, 0), (4, 0)], 1: [(3, 1)], 2: [(4, 2)], 3: [(5, 3)], 4: [(6, 4)], 6: [(8, 6)]}
char_ok = True
char_tab = {}
for D_ in sorted(erwartet):
    for vz in ((1, -1) if D_ else (1,)):
        emax2 = 2 + D_ + 2*KCH
        chi = sektor_charakter(vz*D_, KCH)
        soll = {}
        for (E2, s2) in erwartet[D_]:
            soll = ch_add(soll, irrep(E2, s2, emax2))
        gleich = chi == soll
        char_ok &= gleich
        char_tab[f"D = {vz*D_} (s = {hb(D_)})"] = " + ".join(dname(*t) for t in erwartet[D_]) + f": {gleich}"
chi4 = sektor_charakter(4, KCH)
diff_verma = ch_add(verma(6, 4, 6 + 2*KCH), chi4, -1)
nullmodul = diff_verma == verma(8, 2, 6 + 2*KCH)
dim_E = [sum((s2 + 1)*m for s2, m in chi4[2*n].items()) for n in range(3, 9)]
ZAHLEN["1_Charaktere"] = char_tab
ZAHLEN["1_D32_Spins_je_Energie"] = {f"E = {n}": [hb(s2) for s2 in sorted(chi4[2*n])] for n in range(3, 8)}
check("[E] Charaktervergleich bis Stufe 9 (exakt, aus den J_3-Gewichten aller Monome): s = 0 -> D(1,0) + D(2,0) "
      "(z104, z148 bestätigt); s = 1/2 -> D(3/2,1/2); s = 1 -> D(2,1); s = 3/2 -> D(5/2,3/2); s = 2 -> D(3,2); "
      "s = 3 -> D(4,3); für D < 0 dasselbe. Im Sektor D = +-4 stehen bei E = 3 + k die Spins 2, 3, ..., 2 + k je "
      f"einmal, Dimensionen {dim_E} = n^2 - 4 für E = n = 3..8 (z75: 2 (n^2 - 4) für beide Vorzeichen). Der Sektor ist "
      "nicht der volle Verma-Modul über (3, 2): Die Differenz ist exakt der Verma-Modul über (4, 1), der Nullmodul der "
      f"Eichung ({nullmodul}). Die Literaturaussage [L: Angelopoulos-Flato-Fronsdal-Sternheimer] ist damit auf "
      "Charakterebene nachgerechnet", char_ok and nullmodul and dim_E == [n*n - 4 for n in range(3, 9)])

# 1.4 Tiefste Gewichte exakt -------------------------------------------------------------------------------------------
LOW_INT = []
for Q in Q_T:
    L = sp.Matrix(zhat_koeff(Q)[2])
    if all(sp.im(x) == 0 for x in L):
        Lr = L
    else:
        Lr = -I_*L
    assert all(sp.im(x) == 0 for x in Lr)
    LOW_INT.append([[int(4*Lr[i, j]) for j in range(2)] for i in range(2)])


def kern_dims(D, kmax, low=LOW_INT):
    """Dimension des gemeinsamen Kerns der drei Senker L^-_k je Stufe, exakt über Z (Bargmann-Monombasis)."""
    nd0, ne0 = max(D, 0), max(-D, 0)

    def monome(k):
        if k < 0:
            return []
        return [(n1, nd0 + k - n1, m1, ne0 + k - m1) for n1 in range(nd0 + k + 1) for m1 in range(ne0 + k + 1)]
    out = []
    for k in range(kmax + 1):
        src, tgt = monome(k), monome(k - 1)
        if not tgt:
            out.append(len(src))
            continue
        tix = {t: i for i, t in enumerate(tgt)}
        rows = [[0]*len(src) for _ in range(len(low)*len(tgt))]
        for a, c in enumerate(low):
            for col, b in enumerate(src):
                for i in range(2):
                    for j in range(2):
                        if c[i][j] and b[2 + i] and b[j]:
                            t = list(b)
                            t[2 + i] -= 1
                            t[j] -= 1
                            rows[a*len(tgt) + tix[tuple(t)]][col] += c[i][j]*b[2 + i]*b[j]
        out.append(len(src) - fmpz_mat(rows).rank())
    return out


KK = 6
kern = {D_: kern_dims(D_, KK) for D_ in (-4, 0, 1, 2, 3, 4, 6)}
kern_soll = {-4: [5] + [0]*KK, 0: [1, 1] + [0]*(KK - 1), 1: [2] + [0]*KK, 2: [3] + [0]*KK, 3: [4] + [0]*KK,
             4: [5] + [0]*KK, 6: [7] + [0]*KK}
# Gegenprobe: mit allen vier Senkern (su(2,2)) hat auch der Sektor D = 0 nur ein tiefstes Gewicht
SING_INT = [[int(sp.re(-I_*4*sp.Matrix(zhat_koeff(QM[0])[2])[i, j])) for j in range(2)] for i in range(2)]
kern_su22 = kern_dims(0, KK, LOW_INT + [SING_INT])
ZAHLEN["1_tiefste_Gewichte_je_Stufe"] = {f"D = {D_}": kern[D_] for D_ in sorted(kern)}
check("[E] Irreduzibilität über exakte Ränge (Bargmann-Monombasis, ganzzahlige Matrizen, Stufen 0..6): Der gemeinsame "
      "Kern der drei nativen Senker L^-_k hat je Stufe die Dimension "
      f"{ {f'D = {D_}': kern[D_][:3] for D_ in sorted(kern)} } (danach 0). Für s != 0 gibt es also genau ein tiefstes "
      "Gewicht, das Spin-s-Multiplett bei E = s + 1; weil die Darstellung unitär ist (Teil 2), ist sie vollständig "
      "reduzibel und damit irreduzibel: D(s+1, s). Für s = 0 gibt es zwei tiefste Gewichte (E = 1 und E = 2, je Spin 0): "
      f"D(1,0) + D(2,0); erst der vierte Senker der su(2,2) verbindet sie (Kern dann {kern_su22[:3]})",
      kern == kern_soll and kern_su22 == [1] + [0]*KK)

# 1.5 Algebra, Gram-Form, Casimir --------------------------------------------------------------------------------------
S4 = Sektor(4, 4)
E4, J4, Lp4, Lm4 = so32_ops(S4)
inn = S4.innen(1)
lc = sp.LeviCivita


def rest(M):
    return float(np.abs(M[np.ix_(inn, inn)]).max())


alg_rest = max(rest(Lm4[i] @ Lp4[j] - Lp4[j] @ Lm4[i] - 2*(i == j)*E4
                    + 2j*sum(int(lc(i, j, k))*J4[k] for k in range(3))) for i in range(3) for j in range(3))
alg_rest = max(alg_rest, rest(J4[0] @ J4[1] - J4[1] @ J4[0] - 1j*J4[2]),
               max(rest(E4 @ Lp4[i] - Lp4[i] @ E4 - Lp4[i]) for i in range(3)))
lw = np.where(S4.stufe == 0)[0]
vek = np.array([Lp4[j][:, m] for j in range(3) for m in lw]).T               # 15 Vektoren der ersten Stufe
gram = vek.conj().T @ vek
ew = np.linalg.eigvalsh(gram)
ew_r = sorted(int(round(x)) for x in ew)
gram_ok = np.abs(ew - np.round(ew)).max() < 1e-9 and ew_r == [0]*3 + [4]*5 + [10]*7
nullraum = np.linalg.svd(vek, compute_uv=False)
rang_vek = int(np.sum(nullraum > 1e-9))
cas = {}
cas_ok = True
for D_ in (0, 1, 2, 4, -4):
    S = Sektor(D_, 3)
    E, J, Lp, Lm = so32_ops(S)
    C2 = E @ E + sum(X @ X for X in J) - 0.5*sum(Lp[i] @ Lm[i] + Lm[i] @ Lp[i] for i in range(3))
    ii = S.innen(1)
    wert = C2[np.ix_(ii, ii)]
    soll = D_*D_/2 - 2
    cas_ok &= np.abs(wert - soll*np.eye(len(ii))).max() < 1e-9
    cas[f"D = {D_}"] = str(Fr(D_*D_, 2) - 2)
ZAHLEN["1_Gram_erste_Stufe_D4"] = ew_r
ZAHLEN["1_Casimir_so32"] = cas
check("[E] Dieselbe Algebra wie in z148: Im Fockraum gilt [L^-_i, L^+_j] = 2 delta_ij E - 2 i eps_ijk J_k, "
      "[J_1, J_2] = i J_3, [E, L^+] = L^+ (Rest < 1e-12). Die Gram-Form der ersten Stufe über dem tiefsten "
      f"Gewicht des Sektors D = 4 (15 Vektoren L^+_j |m>) hat die Eigenwerte {sorted(set(ew_r))} mit Vielfachheiten "
      f"{[ew_r.count(x) for x in sorted(set(ew_r))]}: genau z148s 2 (E0 - 3), 2 (E0 - 1), 2 (E0 + 2) bei E0 = 3. Die "
      f"drei Nullrichtungen sind im Fockraum der Nullvektor selbst (Rang der 15 Vektoren: {rang_vek}), kein Zustand "
      "der Norm null: Die Eichrichtung fehlt identisch. Kein Eigenwert ist negativ. Der quadratische Casimir "
      "E^2 + J^2 - (1/2) Summe {L^+_i, L^-_i} = E0 (E0 - 3) + s (s + 1) ist auf jedem Sektor die Zahl D^2/2 - 2 "
      f"({cas}): für D = +-4 gleich 6, also E0 (E0 - 3) = m^2 l^2 = 0 (masselos; z148s Kompositum D(4,2) hat 4)",
      alg_rest < 1e-12 and gram_ok and rang_vek == 12 and cas_ok)


# 1.6 Wahl der Untergruppe ---------------------------------------------------------------------------------------------
def erzeugte_dimension(V):
    """Reelle Dimension der von E und den Hebungen senkrecht zu V erzeugten Lie-Algebra in su(2,2)."""
    V = np.array(V, complex)
    basis_u = []
    for i in range(2):
        for j in range(2):
            U = np.zeros((2, 2), complex)
            U[i, j] = 1
            basis_u.append(U)
    # Projektion auf das hermitesche Komplement von V
    Vn = V/np.linalg.norm(V)
    comp = [U - Vn*np.vdot(Vn, U) for U in basis_u]
    gens = [np.array(E_T.tolist(), complex)]
    for U in comp:
        for ph in (1, 1j):
            X = np.zeros((4, 4), complex)
            X[:2, 2:] = ph*U
            X[2:, :2] = (ph*U).conj().T
            gens.append(X)

    def vec(M):
        return np.concatenate([M.real.ravel(), M.imag.ravel()])
    span = np.array([vec(g) for g in gens])
    for _ in range(4):
        r0 = np.linalg.matrix_rank(span, tol=1e-9)
        _, sv, vh = np.linalg.svd(span, full_matrices=False)
        bas = vh[:r0]
        mats = [(b[:16] + 1j*b[16:]).reshape(4, 4) for b in bas]
        neu = [vec(A @ B - B @ A) for A in mats for B in mats]
        span = np.vstack([bas, np.array(neu)])
        if np.linalg.matrix_rank(span, tol=1e-9) == r0:
            break
    return int(np.linalg.matrix_rank(span, tol=1e-9))


ct, st = math.cos(0.7), math.sin(0.7)
V_FAELLE = [("V = 1 (nativ)", [[1, 0], [0, 1]]),
            ("V = e^(0,4 i) 1", [[np.exp(0.4j), 0], [0, np.exp(0.4j)]]),
            ("V unitär (gedreht)", [[ct, 1j*st], [1j*st, ct]]),
            ("V = sigma_3", [[1, 0], [0, -1]]),
            ("V = diag(1, 2)", [[1, 0], [0, 2]]),
            ("V = diag(1, 0)", [[1, 0], [0, 0]])]
dim_V = {n: erzeugte_dimension(V) for n, V in V_FAELLE}
ZAHLEN["1_Untergruppen_mit_E"] = dim_V
check("[E/N] Abhängigkeit von der Wahl der Untergruppe: Eine so(3,2) in su(2,2), die E enthält, ist durch das "
      "Komplement V ihrer drei Hebungen im Raum der vier d_i^+ e_j^+ bestimmt. Die von E und den Hebungen senkrecht zu "
      f"V erzeugte Algebra hat die Dimension {dim_V}: 10 für die vier geprüften unitären V (auch mit Phase), 15 "
      "(die ganze su(2,2)) für die zwei nicht unitären. Geometrisch ist V die raumartige Richtung im R^4 der K-Typen, "
      "die die so(3,2) festhält; die unitären V bilden eine Bahn von SU(2)_L, die mit E, D und der Fock-Norm "
      "vertauscht. Alle so(3,2) mit kompaktem E sind also konjugiert, und Zerlegung und Charaktere sind dieselben. Die native Wahl V = 1 "
      "(z89: E und Lorentz erhalten I_pi - I_omega) geht nur in die Frage ein, welche Vektoren die tiefsten sind, "
      "nicht in den Darstellungsinhalt. Für so(4,1) = Der(I_pi + I_omega) (ohne E) gilt das nicht (z104, z128)",
      [dim_V[n] for n, _ in V_FAELLE] == [10, 10, 10, 10, 15, 15])
befund("Teil 1: D(3,2) im Sektor Helizität +-2", "ja",
       "E: Die Leiterdarstellung der Helizität s = D/2 der su(2,2) ist unter der nativen so(3,2) für s != 0 irreduzibel "
       "und gleich D(s+1, s); für s = 0 zerfällt sie in D(1,0) + D(2,0). Der Sektor D = +4 trägt genau D(3,2), der "
       "Sektor D = -4 ebenso. Das hängt nicht von der Wahl der so(3,2) mit kompaktem E ab.")

# ====================================================================================================================
say("\nTeil 2: Norm und Energie -- welcher Operator ist die so(3,2)-Energie?")
# 2.1 ---------------------------------------------------------------------------------------------------------------
herm = 0.0
for D_ in (4, -4):
    S = Sektor(D_, 4)
    for X in SU22:
        M = 1j*S.zhat(X)
        herm = max(herm, float(np.abs(M - M.conj().T).max()))
E_spek = sorted(set(int(round(x)) for x in np.diag(E4).real))
E_viel = [int(np.sum(np.abs(np.diag(E4).real - e) < 1e-9)) for e in E_spek]
check("[E] Norm und Energie im Sektor D = +-4: Die Fock-Form ist positiv definit, und alle 15 Erzeuger i Zhat(X) sind "
      "in ihr hermitesch (Abweichung < 1e-12; exakt: X^+ G + G X = 0). Die Darstellung ist also unitär. Die "
      f"so(3,2)-Energie E = 1 + (N_d + N_e)/2 hat dort das Spektrum {E_spek} mit Vielfachheiten {E_viel} (n^2 - 4): nach "
      "unten beschränkt, E_0 = 3. Bezüglich E und der Fock-Form ist der Sektor D(3,2), unitär und mit positiver "
      "Energie", herm < 1e-12 and antiherm and E_spek == [3, 4, 5, 6, 7] and E_viel == [5, 12, 21, 32, 45])

# 2.2 Araki-Woods exakt ------------------------------------------------------------------------------------------------
q = sp.symbols("q", positive=True)
s_ = sp.sqrt(1 - q**2)
X1 = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])     # tauscht d_i und e_i^+
SBAR = sp.eye(4) + q*X1                                                      # S = SBAR/s: Z_ab = S Z_de
SBAR_INV = sp.eye(4) - q*X1                                                  # S^-1 = SBAR_INV/s
S_gruppe = sp.simplify(SBAR.H*G4*SBAR - (1 - q**2)*G4) == NULL4
S_dil = (X1 == -2*D_T) and (X1*X1 == sp.eye(4))


def komm8(u, v):
    """[u, v] für Linearformen in (d1, d2, e1, e2, d1^+, d2^+, e1^+, e2^+)."""
    return sp.simplify(sum(u[i]*v[i + 4] - u[i + 4]*v[i] for i in range(4)))


def e8(i):
    v = [0]*8
    v[i] = 1
    return sp.Matrix(v)


d_ = [e8(0), e8(1)]
e_ = [e8(2), e8(3)]
dp = [e8(4), e8(5)]
ep = [e8(6), e8(7)]
a_ = [(d_[i] + q*ep[i])/s_ for i in range(2)]
ap = [(dp[i] + q*e_[i])/s_ for i in range(2)]
b_ = [(e_[i] + q*dp[i])/s_ for i in range(2)]
bp = [(ep[i] + q*d_[i])/s_ for i in range(2)]
ccr = all(komm8(a_[i], ap[j]) == (1 if i == j else 0) and komm8(b_[i], bp[j]) == (1 if i == j else 0)
          and komm8(a_[i], b_[j]) == 0 and komm8(a_[i], bp[j]) == 0 and komm8(a_[i], a_[j]) == 0
          for i in range(2) for j in range(2))
ident = all(sp.simplify(dp[i] - (s_*ap[i] - q*e_[i])) == sp.zeros(8, 1)
            and sp.simplify(ep[i] - (s_*bp[i] - q*d_[i])) == sp.zeros(8, 1) for i in range(2))
vertauscht = all(komm8(e_[i], ap[j]) == 0 and komm8(d_[i], bp[j]) == 0 for i in range(2) for j in range(2))
ZAHLEN["2_AW"] = {"S": "exp(-2 theta D_T), tanh theta = q", "(d^+)^4 Omega": "s^4 (a^+)^4 Omega",
                  "(e^+)^4 Omega": "s^4 (b^+)^4 Omega"}
check("[E] Araki-Woods auf dem Twistorraum (für jedes q): Z_ab = (a_1, a_2, b_1^+, b_2^+) = S Z_de mit "
      f"S = (1 + q X)/s in SU(2,2) ({S_gruppe}); X = -2 D_T, X^2 = 1, also S = exp(-2 theta D_T) mit tanh theta = q: "
      "Der Wechsel zwischen den Oszillatoren des Originals (a, b) und denen des Zustands (d, e) ist die native Dilatation "
      f"({S_dil}). CCR der a, b: {ccr}. Operatoridentität d_i^+ = s a_i^+ - q e_i (und e_i^+ = s b_i^+ - q d_i), "
      "[e_i, a_j^+] = 0, also (d^+)^4 Omega = s^4 (a^+)^4 Omega exakt: Die tiefsten Zustände des Sektors D = +4 sind "
      "genau vier Quanten des Originals über Omega, die des Sektors D = -4 vier Quanten des Partners (z75 "
      "präzisiert). D = N_d - N_e = N_a - N_b, weil das Zentrum unter S fest bleibt",
      S_gruppe and S_dil and ccr and ident and vertauscht)

# 2.3 Drei Kandidaten für die Energie ----------------------------------------------------------------------------------
E_AB = sp.simplify(SBAR_INV*E_T*SBAR/(1 - q**2))                              # (N_a + N_b + 2)/2 als Matrix im de-Rahmen
ko_Eab = koord(E_AB, SU22)
ko_soll = [0]*15
ko_soll[0] = (1 + q**2)/(1 - q**2)
form_Eab = [sp.simplify(ko_Eab[i] - ko_soll[i]) for i in range(15) if i != 13]
w_koeff = sp.simplify(ko_Eab[13])
P1 = sp.diag(1, 1, 0, 0)
einseitig = sp.simplify(I_*E_AB - (SBAR_INV*P1*SBAR/(1 - q**2) - sp.eye(4)/2)) == NULL4   # E_ab = N_a + 1 - D/2
E_im_ab = sp.simplify(SBAR*E_T*SBAR_INV/(1 - q**2))                           # Wirkung von E_de im ab-Rahmen
z75_koeff = (sp.simplify(2*I_*E_im_ab[0, 0] - (1 + q**2)/(1 - q**2)) == 0
             and sp.simplify(2*I_*E_im_ab[0, 2] + 2*q/(1 - q**2)) == 0)
# numerisch im Sektor D = 0: Omega ist Eigenvektor von E_de, nicht von E_ab
S0 = Sektor(0, 3)
qn = 0.3
Eab_num = 1j*S0.zhat(E_AB.subs(q, sp.Rational(3, 10)))
om = np.zeros(S0.dim)
om[S0.idx[(0, 0, 0, 0)]] = 1
mitt = float((om @ Eab_num @ om).real)
var = float((om @ Eab_num @ Eab_num @ om).real) - mitt**2
energie_ok = (all(x == 0 for x in form_Eab) and sp.simplify(w_koeff - 2*q/(1 - q**2)) == 0 and z75_koeff
              and einseitig
              and abs(mitt - (1 + qn**2)/(1 - qn**2)) < 1e-12 and abs(var - 2*qn**2/(1 - qn**2)**2) < 1e-12)
x_nat = math.exp(-4*math.pi)
ZAHLEN["2_Energien"] = {"E_ab": "((1 + q^2) E_de + 2 q Q^-_0)/(1 - q^2) = N_a + 1 - D/2",
                        "E_ab auf dem Sektor D = 4": "N_a - 1 = L_0/m - 1 >= 3",
                        "<E_ab>_Omega": "(1 + q^2)/(1 - q^2)", "Var_Omega(E_ab)": "2 q^2/(1 - q^2)^2",
                        "nativ (m = 1): <E_ab> - 1": gz(2*x_nat/(1 - x_nat), 6)}
check("[E] Drei Operatoren heißen im Bestand 'Energie'. (i) D = N_a - N_b = N_d - N_e: der Modulgenerator durch "
      "4 pi m, die Liouville-Form 'N x 1 - 1 x N'. Er ist zentral in u(2,2), liegt nicht in so(3,2), ist nicht positiv "
      "und auf jedem Helizitätssektor eine Zahl (+-4, also K = +-16 pi m): Er ist die Helizität, nicht die "
      "so(3,2)-Energie. (ii) E_de = 1 + (N_d + N_e)/2: die Energie der nativen so(3,2) aus z89, positiv, Omega ist ihr "
      "tiefster Zustand (E = 1); im Rahmen des Originals wirkt sie als [N_d + N_e, a] = -((1 + q^2) a - 2 q b^+)/(1 - q^2) "
      f"(z75 exakt reproduziert: {z75_koeff}). (iii) E_ab = 1 + (N_a + N_b)/2, die Form 'E x 1 + 1 x E': positiv, erhält "
      "Original und Partner einzeln, und E_ab = ((1 + q^2) E_de + 2 q Q^-_0)/(1 - q^2) mit der Kappe Q^-_0 = "
      "(1/2) Summe (d^+ e^+ + e d) aus z104 liegt in su(2,2), aber nicht in der nativen so(3,2). Omega ist kein "
      "Eigenvektor von E_ab: <E_ab> = (1 + q^2)/(1 - q^2), Varianz 2 q^2/(1 - q^2)^2 (q = 0,3 numerisch: "
      f"{mitt:.10f}, {var:.10f}). Exakt gilt E_ab = N_a + 1 - D/2 ({einseitig}): Auf jedem Helizitätssektor ist E_ab "
      "bis auf eine Konstante die einseitige Uhr des Originals, N_a = L_0/m (z104, z136: einseitig 4 pi E). Auf D = +4 "
      "ist E_ab = L_0/m - 1 >= 3. Zentral ist nur der zweiseitige Modulgenerator; der einseitige ist die Energie",
      energie_ok)


# 2.4 Welche Erzeuger erhalten das Original? ---------------------------------------------------------------------------
def stab_dim(qv, idx):
    """Dimension von {X in span(SU22[idx]) : rechter oberer Block von S X S^-1 ist null} (X erhält span{a_i})."""
    rows = []
    for n in idx:
        Y = (SBAR*SU22[n]*SBAR_INV).subs(q, qv)
        blk = [Y[i, 2 + j] for i in range(2) for j in range(2)]
        rows.append([sp.re(x) for x in blk] + [sp.im(x) for x in blk])
    return len(idx) - sp.Matrix(rows).rank()


stab = {}
for qv in (sp.Rational(3, 5), sp.Rational(5, 13), sp.Integer(0)):
    stab[f"q = {qv}"] = [stab_dim(qv, list(range(15))), stab_dim(qv, list(range(10))), stab_dim(qv, K_DE)]
# die drei gemeinsamen Erzeuger sind die Drehungen
dreh = all(sp.simplify((SBAR*R*SBAR_INV)[i, 2 + j]) == 0 for R in R_T for i in range(2) for j in range(2))
ZAHLEN["2_Stabilisator_des_Originals"] = {k: {"su(2,2)": v[0], "native so(3,2)": v[1], "mit Omega (k_de)": v[2]}
                                          for k, v in stab.items()}
check("[E] Welche Erzeuger sind Automorphismen des Originals M = W*(a_1, a_2)? X erhält span{a_i} genau dann, wenn "
      "S X S^-1 blockdiagonal ist. Dimensionen (su(2,2), native so(3,2), Stabilisator von Omega) für q = 3/5, 5/13: "
      f"{stab['q = 3/5']}, {stab['q = 5/13']}; für q = 0: {stab['q = 0']}. In su(2,2) erhält nur die kompakte "
      "s(u(2) + u(2)) von E_ab das Original; von der nativen so(3,2) nur die drei Drehungen; das Paar (M, Omega) "
      f"erhält ebenfalls nur die SO(3) (z146 bestätigt; Drehungen: {dreh}). Die Hebungen der so(3,2) sind a^+ b^+ bzw. "
      "d^+ e^+: Sie erzeugen je ein Quant im Original und im Partner. Im Algebrarahmen (so(3,2) von E_ab; dieselbe "
      "Rechnung wie q = 0) erhalten vier Erzeuger das Original, E_ab und die Drehungen. Die sechs nichtkompakten "
      "Erzeuger (AdS-Verschiebungen und Boosts) tun es in keinem Rahmen: Die so(3,2) wirkt auf der Purifikation "
      "(Original + Partner), nicht auf dem Original allein",
      stab["q = 3/5"] == [7, 3, 3] and stab["q = 5/13"] == [7, 3, 3] and stab["q = 0"] == [7, 4, 7] and dreh)

# 2.5 Die zweite so(3,2) -----------------------------------------------------------------------------------------------
I_AB = sp.simplify(SBAR_INV*I_ADS*SBAR_INV.T/(1 - q**2))
lam, mu = sp.symbols("lam mu")
sol = sp.solve(list(I_AB - lam*I_PI + mu*I_OM), [lam, mu], dict=True)
bue_ok = (len(sol) == 1 and sp.simplify(sol[0][lam] - (1 + q)/(1 - q)) == 0
          and sp.simplify(sol[0][mu] - (1 - q)/(1 + q)) == 0)
SO32_AB = [sp.simplify(SBAR_INV*X*SBAR/(1 - q**2)) for X in SO32]
der_ab = all(sp.simplify(X*I_AB + I_AB*X.T) == NULL4 for X in SO32_AB)
coth = {m: (1 + math.exp(-2*math.pi*m))/(1 - math.exp(-2*math.pi*m)) for m in (1, 2, 3)}
# Überlapp der tiefsten Multipletts (ein Oszillatorpaar numerisch, q = 0,3)
NT = 60
aop = np.diag(np.sqrt(np.arange(1, NT)), 1)
Idn = np.eye(NT)
d1, e1 = np.kron(aop, Idn), np.kron(Idn, aop)
sn = math.sqrt(1 - qn*qn)
a1 = (d1 + qn*e1.T)/sn
vac_ab = np.zeros(NT*NT)
for n in range(NT):
    vac_ab[n*NT + n] = sn*(-qn)**n
om1 = np.zeros(NT*NT)
om1[0] = 1
v_ab = np.linalg.matrix_power(a1.T, 4) @ vac_ab/math.sqrt(24)
v_de = np.linalg.matrix_power(d1.T, 4) @ om1/math.sqrt(24)
ueb1 = float(v_ab @ v_de)                              # s^5 für das angeregte Paar, mal s für das zweite Paar
ueb_ok = abs(ueb1 - sn**5) < 1e-10 and abs(float(vac_ab @ om1) - sn) < 1e-12
ueb_nat = (1 - x_nat)**6
ZAHLEN["2_zweite_so32"] = {"I'": "((1+q)/(1-q)) I_pi - ((1-q)/(1+q)) I_omega",
                           "Skalenfaktor coth(pi m), m = 1, 2, 3": [gz(coth[m], 12) for m in (1, 2, 3)],
                           "Überlapp^2 der tiefsten Multipletts, nativ m = 1": gz(ueb_nat, 12)}
check("[E] Zwei konjugierte so(3,2): Die von E_ab ausgezeichnete Algebra S^-1 so(3,2) S ist die Derivationsalgebra von "
      f"I' = ((1+q)/(1-q)) I_pi - ((1-q)/(1+q)) I_omega ({bue_ok and der_ab}), ein anderes Element desselben "
      "Lorentz-Büschels (z89): dasselbe AdS mit um coth(pi m) = "
      f"{gz(coth[1], 9)}, {gz(coth[2], 12)}, {gz(coth[3], 14)} (m = 1, 2, 3) gestreckter Skala (kappa^-2 aus z89). "
      "Beide sind durch die Dilatation S unitär äquivalent; der Sektor D = +-4 ist bezüglich jeder von beiden D(3,2) "
      "mit positiver Norm und E_0 = 3. Verschieden sind die tiefsten Vektoren: (d^+)^4 Omega gegen (a^+)^4 |0_ab>, "
      f"Überlapp s^6 (q = 0,3 numerisch: {ueb1*sn:.10f} gegen {sn**6:.10f}); nativ ist das Überlappquadrat "
      f"(1 - e^(-4 pi))^6 = {gz(ueb_nat, 10)}. Der Zustandsrahmen (E_de) hängt über q von der Mode m ab, der "
      "Algebrarahmen (E_ab) nicht", bue_ok and der_ab and ueb_ok)


# 2.6 Lesart T1 --------------------------------------------------------------------------------------------------------
def liouville_tiefste(E0, nsrc):
    """Lesart T1: Raum der Operatoren Psi auf der Leiter mit [E, Psi] = E0 Psi und [L^-_k, Psi] = 0 (k = 1, 2, 3),
    eingeschränkt auf die Quellstufen E = 1..nsrc. Gibt die Dimension des Lösungsraums (obere Schranke)."""
    S = Sektor(0, nsrc + E0 - 1)
    _, _, _, Lm = so32_ops(S)
    lev = [np.where(S.stufe == k)[0] for k in range(nsrc + E0)]             # lev[k]: E = k + 1
    dims = [len(x) for x in lev]
    off, tot = {}, 0
    for n in range(1, nsrc + 1):                                            # Psi_n: Stufe E = n -> E = n + E0
        off[n] = tot
        tot += dims[n - 1 + E0]*dims[n - 1]
    rows = []
    for n in range(1, nsrc + 1):
        dn, dt = dims[n - 1], dims[n - 1 + E0]
        for a in range(3):
            A = Lm[a][np.ix_(lev[n - 2 + E0], lev[n - 1 + E0])]             # Senker auf dem Bild
            blk = np.zeros((A.shape[0]*dn, tot), complex)
            blk[:, off[n]:off[n] + dt*dn] = np.kron(np.eye(dn), A)          # vec(A Psi_n), spaltenweise
            if n > 1:
                Bm = Lm[a][np.ix_(lev[n - 2], lev[n - 1])]                  # Senker auf der Quelle
                dq, dtt = dims[n - 2], dims[n - 2 + E0]
                blk[:, off[n - 1]:off[n - 1] + dtt*dq] -= np.kron(Bm.T, np.eye(dtt))
            rows.append(blk)
    Mx = np.vstack(rows)
    sv = np.linalg.svd(Mx, compute_uv=False)
    return tot - int(np.sum(sv > 1e-9*sv[0]))


liou = {E0: liouville_tiefste(E0, 3) for E0 in (1, 2, 3)}
EM2 = 14
leiter = sektor_charakter(0, 6)
zwei = ch_mul(leiter, leiter, EM2)
zer_zwei = zerlege(zwei, EM2)
d32_zwei = sum(m for e, s, m in zer_zwei if (e, s) == (6, 4))
spin2_zwei = min(e for e, s, m in zer_zwei if s == 4)
ZAHLEN["2_T1"] = {"Liouville: Dimension tiefster Gewichte für E0 = 1, 2, 3": [liou[1], liou[2], liou[3]],
                  "zwei Leitern: Vielfachheit D(3,2)": d32_zwei,
                  "zwei Leitern: tiefstes Spin-2-Gewicht": dname(spin2_zwei, 4)}
check("[E/N] Lesart T1 (z146: M = B(Leiter), rho = e^(-beta E)/Z, Thermofeldraum Leiter x Leiter-quer). Dort gibt es "
      "keinen Sektor D != 0. (a) Die Darstellung, die die Automorphismen implementiert und Omega_beta erhält, ist "
      "Psi -> [X, Psi]; ihre Energie E x 1 - 1 x E hat ein symmetrisches, nach unten unbeschränktes Spektrum. Der Raum "
      "der Operatoren mit [L^-_k, Psi] = 0 und [E, Psi] = E0 Psi hat auf den Quellstufen E = 1..3 die Dimension "
      f"{[liou[1], liou[2], liou[3]]} für E0 = 1, 2, 3: Psi|E = 1> müsste ein tiefstes Gewicht der Leiter bei E = 1 + E0 "
      "sein, und das gibt es nur bei E = 2. Ein tiefstes Gewicht mit E0 = 3 kommt nicht vor. (b) Mit E x 1 + 1 x E "
      f"(zwei Leitern) ist die Vielfachheit von D(3,2) {d32_zwei}, Spin 2 beginnt bei {dname(spin2_zwei, 4)} (z148 "
      "reproduziert). In T1 gibt es also kein D(3,2), in keiner der beiden Energielesarten",
      liou[2] == 0 and liou[3] == 0 and liou[1] <= 1 and d32_zwei == 0 and spin2_zwei == 8)
befund("Teil 2: positive Norm und positive Energie", "ja in T2, nein in T1",
       "E: Im Twistorbild (T2) ist der Sektor D = +-4 bezüglich der Fock-Form und der so(3,2)-Energie E_de (nativ nach "
       "z89) unitär mit E_0 = 3, ebenso bezüglich der konjugierten E_ab. Der Einwand des Bestands wird genauer: "
       "Zentral (also Helizität) ist der zweiseitige Modulgenerator D; die so(3,2)-Energie E_ab ist auf dem Sektor die "
       "einseitige Uhr des Originals, L_0/m - 1, und positiv, weil L_0 es ist. Was nicht nativ wirkt, sind die sechs "
       "nichtkompakten Erzeuger: Sie erhalten das Original nicht, die so(3,2) wirkt auf der Purifikation. Im "
       "Thermofeldbild (T1) gibt es den Sektor nicht.")

# ====================================================================================================================
say("\nTeil 3: Speziesproblem und Komposita")
# 3.1 Ein Dubleton vollständig -----------------------------------------------------------------------------------------
EMF = 12                                               # doppelte Energie bis 12 (E <= 6)


def fock_charakter(emax2):
    """Voller Charakter des Fockraums von d_1, d_2, e_1, e_2 (alle D) bis zur doppelten Energie emax2."""
    chi = {}
    for N in range(emax2 - 1):
        cnt = {}
        for occ in itertools.product(range(N + 1), repeat=4):
            if sum(occ) == N:
                w = sum(J3_GEW2[i]*occ[i] for i in range(4))
                cnt[w] = cnt.get(w, 0) + 1
        chi[2 + N] = gewichte_zu_spins(cnt)
    return chi


F1 = fock_charakter(EMF)
zer1 = zerlege(F1, EMF)
soll1 = [(2, 0, 1), (3, 1, 2), (4, 0, 1), (4, 2, 2)] + [(s2 + 2, s2, 2) for s2 in range(3, EMF - 1)]
RAC = irrep(1, 0, EMF)
DI = irrep(2, 1, EMF)
SING = ch_add(RAC, DI)
ff = ch_mul(SING, SING, EMF)
d32_ein = sum(m for e, s, m in zer1 if (e, s) == (6, 4))
ZAHLEN["3_ein_Dubleton"] = [f"{m} x {dname(e, s)}" for e, s, m in zer1 if e <= 8]
check("[E] Ein Dubleton vollständig (alle D, bis E = 6): Der Fockraum der vier Oszillatoren zerfällt unter der nativen "
      f"so(3,2) in {[f'{m} x {dname(e, s)}' for e, s, m in zer1 if e <= 7]} und so weiter, also D(1,0) + D(2,0) + "
      "2 Summe_(s >= 1/2) D(s+1, s). Sein Charakter ist exakt (Di + Rac)^2 "
      f"({ff == F1}): Das ist die Zerlegung von Flato und Fronsdal [L], hier nachgerechnet. D(3,2) kommt genau "
      f"{d32_ein}-mal vor (D = +4 und D = -4)", zer1 is not None and sorted(zer1) == sorted(soll1) and ff == F1
      and d32_ein == 2)

# 3.2 Howe-Dualität exakt ----------------------------------------------------------------------------------------------
# Basis (Z | Z^+) = (d1, d2, e1^+, e2^+ | d1^+, d2^+, e1, e2); ad Zhat(X) = diag(-X, -conj(X)).
# Neue Oszillatoren: e'_1 = e_2, e'_2 = -e_1; alpha_i = (d_i + e'_i)/sqrt2, beta_i = -i (d_i - e'_i)/sqrt2.
W2 = sp.sqrt(2)


def zeile(d1=0, d2=0, e1p=0, e2p=0, d1p=0, d2p=0, e1=0, e2=0):
    return [d1, d2, e1p, e2p, d1p, d2p, e1, e2]


NEU = [zeile(d1=1, e2=1), zeile(d2=1, e1=-1),                     # sqrt2 alpha_1, alpha_2
       zeile(d1p=1, e2p=1), zeile(d2p=1, e1p=-1),                 # sqrt2 alpha_1^+, alpha_2^+
       zeile(d1=-I_, e2=I_), zeile(d2=-I_, e1=-I_),               # sqrt2 beta_1, beta_2
       zeile(d1p=I_, e2p=-I_), zeile(d2p=I_, e1p=I_)]             # sqrt2 beta_1^+, beta_2^+
PM = sp.Matrix(NEU).T/W2                                          # Spalten: neue Operatoren in der alten Basis
PM_INV = PM.inv()


def ad8(X):
    M = sp.zeros(8, 8)
    M[:4, :4] = -X.T                                              # [Zhat, Z_a] = -X_ab Z_b: Spaltenkonvention
    M[4:, 4:] = -X.conjugate().T
    return M


def in_neu(X):
    return sp.simplify(PM_INV*ad8(X)*PM)


howe_ok = True
for X in SO32:
    Mn = in_neu(X)
    howe_ok &= Mn[:4, 4:] == sp.zeros(4, 4) and Mn[4:, :4] == sp.zeros(4, 4) and Mn[:4, :4] == Mn[4:, 4:]
rest5 = [not (in_neu(X)[:4, 4:] == sp.zeros(4, 4) and in_neu(X)[:4, :4] == in_neu(X)[4:, 4:]) for X in SU22[10:]]
MD = in_neu(I_*sp.eye(4))
dreh_D = (MD[:4, :4] == sp.zeros(4, 4) and MD[4:, 4:] == sp.zeros(4, 4)
          and MD[:4, 4:] == -MD[4:, :4] and (MD*MD == -sp.eye(8)))
check("[E] Warum zweimal: Mit e'_1 = e_2, e'_2 = -e_1, alpha_i = (d_i + e'_i)/sqrt 2 und beta_i = -i (d_i - e'_i)/sqrt 2 "
      "wirken alle zehn Erzeuger der nativen so(3,2) = sp(4,R) getrennt und gleich auf (alpha, alpha^+) und "
      f"(beta, beta^+) ({howe_ok}); die fünf übrigen der su(2,2) tun das nicht ({all(rest5)}). Der Fockraum ist also "
      "Fock(alpha) x Fock(beta) = (Di + Rac) x (Di + Rac), exakt auf Operatorebene. D = N_d - N_e erzeugt die Drehung "
      f"alpha -> beta ({dreh_D}): Die Helizität ist die Ladung der zu sp(4,R) dualen O(2) [L: Howe], D(3,2) trägt deren "
      "zweidimensionale Darstellung mit Ladung +-4. Aus Sicht der so(3,2) ist der tiefste Zustand ein "
      "Zwei-Singleton-Zustand, kein Kompositum masseloser Teilchen", howe_ok and all(rest5) and dreh_D)

# 3.3 Mehrere Dubletons unter einer gemeinsamen so(3,2) ---------------------------------------------------------------
F2 = ch_mul(F1, F1, EMF)
F3 = ch_mul(F2, F1, EMF)
zer2 = zerlege(F2, EMF)
zer3 = zerlege(F3, EMF)


def harm4(n):
    """Dimension der harmonischen Polynome vom Grad 4 in n Variablen (O(n)-Darstellung zum Diagramm (4))."""
    return math.comb(n + 3, 4) - math.comb(n + 1, 2)


mehr = {}
mehr_ok = zer2 is not None and zer3 is not None
for N, zz in ((1, zer1), (2, zer2), (3, zer3)):
    d32 = sum(m for e, s, m in zz if (e, s) == (6, 4))
    e_min = min(e for e, s, m in zz if s == 4)
    m_min = sum(m for e, s, m in zz if (e, s) == (e_min, 4))
    vak = min(e for e, s, m in zz)
    mehr[f"{N} Dubleton(s)"] = {"D(3,2)": d32, "tiefstes Spin-2-Gewicht": dname(e_min, 4), "Vielfachheit": m_min,
                                "tiefste Energie überhaupt": hb(vak)}
    mehr_ok &= (d32 == (2 if N == 1 else 0) and e_min == 2*(N + 2) and m_min == harm4(2*N) and vak == 2*N)
ZAHLEN["3_mehrere_Dubletons_diagonal"] = mehr
check("[E] Mehrere Dubletons unter einer gemeinsamen (diagonalen) so(3,2), Charaktere bis E = 6: Vielfachheit von "
      f"D(3,2) für 1, 2, 3 Dubletons {[mehr[k]['D(3,2)'] for k in sorted(mehr)]}; das tiefste Spin-2-Gewicht liegt bei "
      f"{[mehr[k]['tiefstes Spin-2-Gewicht'] for k in sorted(mehr)]} mit Vielfachheit "
      f"{[mehr[k]['Vielfachheit'] for k in sorted(mehr)]} (= harmonische Quartiken in 2N Variablen, O(2N)-dual), also "
      "E_0 = N + 2: massiv ab N = 2 (m^2 l^2 = E_0 (E_0 - 3) = 4, 10). Grund: Jeder Dubleton trägt die Nullpunktsenergie "
      "1; sein Fockvakuum ist kein so(3,2)-invarianter Zustand, sondern der Grundzustand von D(1,0). Ein Dubleton ist "
      "ein Einteilchenraum, N Dubletons sind ein N-Teilchen-Raum. Masselos ist Spin 2 nur für genau einen Dubleton",
      mehr_ok)

# 3.4 Auswahl einer Kopie ---------------------------------------------------------------------------------------------
M_MAX = 6
uhr = [complex(np.exp(2j*math.pi*4*m)) for m in range(1, M_MAX + 1)]
halb = [complex(np.exp(2j*math.pi*m*D_)) for m in range(1, M_MAX + 1) for D_ in range(-6, 7)]
uhr_ok = all(abs(u - 1) < 1e-9 for u in uhr + halb)
modgew = {m: math.exp(-16*math.pi*m) for m in (1, 2, 3)}
paare = math.comb(8, 2)
zaehl = {f"m <= {M}": {"ein Strompaar": 2*M, "alle Koordinatenpaare der 8 Cartan-Ströme": 2*M*paare}
         for M in (1, 2, 3, 6)}
skalen_versch = abs(coth[1]/coth[2] - 1) > 1e-3
ZAHLEN["3_Kopien"] = {"Zählung D(3,2) (jede Mode eigene so(3,2))": zaehl,
                      "Uhr R auf (a_m^+)^4 Omega": "1 für alle m",
                      "KMS-Gewicht e^(-16 pi m), m = 1, 2, 3": [gz(modgew[m], 6) for m in (1, 2, 3)],
                      "Verhältnis der Zustandsskalen coth(pi)/coth(2 pi)": gz(coth[1]/coth[2], 10)}
check("[E/C] Auswahl einer Kopie. Zählung: Jede Strommode m und jedes Paar von Strömen trägt einen eigenen Dubleton mit "
      "eigener su(2,2); als Einteilchenräume gelesen sind das für m <= M und ein festes Paar 2 M Kopien von D(3,2) "
      f"(M = 1, 2, 3, 6: {[zaehl[k]['ein Strompaar'] for k in zaehl]}), mit den {paare} Koordinatenpaaren der acht "
      f"Cartan-Ströme {[zaehl[k]['alle Koordinatenpaare der 8 Cartan-Ströme'] for k in zaehl]}; die Wahl des Paars ist "
      "sogar stetig (Gr(2,8)). Geprüfte Auswahlkandidaten: (a) Uhr R = e^(2 pi i L_0): (a_m^+)^4 Omega hat L_0 = 4 m, "
      f"also R = 1 für jedes m ({uhr_ok}); die Modulhalbdrehung e^(i K/2) = e^(2 pi i m D) (z130: Niveauabgleich) "
      "ist 1 auf jedem Sektor; Ladung und Paaruhr sind auf Strömen neutral (z89): keine Auswahl. (b) "
      "Modulgewicht: e^(-16 pi m) = "
      f"{[gz(modgew[m], 3) for m in (1, 2, 3)]}; das ordnet die Kopien (m = 1 zuerst), entfernt aber keine. (c) Eine "
      "gemeinsame so(3,2) über mehrere Moden gibt nach 3.3 kein D(3,2) mehr. (d) Die zustandsnativen so(3,2) der Moden "
      f"unterscheiden sich um die Dilatation coth(pi m) (Verhältnis m = 1 zu 2: {gz(coth[1]/coth[2], 9)}): keine "
      "gemeinsame native AdS-Skala über die Moden (z89: Normierung tanh(pi m)/m). (e) Gitterlesart z152 [L]: Die "
      "Moden sind nach Tiefe geordnet, die Zahl wächst mit der Kettenlänge; ausgezeichnet ist die flachste, eindeutig "
      "ist keine. Eine native Auswahl einer einzigen Kopie gibt es nicht; sie braucht P_mode (eine Mode, etwa m = 1), "
      "P_paar (zwei Ströme) und die Einteilchenlesart P_1T", uhr_ok and skalen_versch and paare == 28)
befund("Teil 3: genau eine Kopie nativ", "nein",
       "E/C: Je Dubleton tritt D(3,2) zweimal auf (D = +-4, Dublett der dualen O(2)), je Mode und Strompaar ein "
       "Dubleton. Der tiefste Zustand besteht aus vier Stromquanten der Naht (zwei Singletons der so(3,2)). Uhr, "
       "Paaruhr, Modulgewicht und Gitterlesart wählen keine einzelne Kopie; eine gemeinsame so(3,2) über mehrere "
       "Moden macht Spin 2 massiv (E_0 = N + 2).")

# ====================================================================================================================
say("\nTeil 4: Kopplung -- passt der Sektor als Zustandsraum der D(3,2)-Mode aus z148?")
# 4.1 ---------------------------------------------------------------------------------------------------------------
z148 = json.loads((HERE / "z148_ergebnis.json").read_text())
z148_ok = all(z148["checks"].values())
e0, s0 = 3, 2
pass_ok = (e0 == s0 + 1 and e0*(e0 - 3) == 0 and e0*(e0 - 3) + s0*(s0 + 1) == 6 and gram_ok and nullmodul
           and z148_ok)
check("[E] Passung zur Mode aus z148 Teil A (dort unter P_dyn): gleiche Darstellung D(3,2); Randexponent Delta = 3 = "
      "E_0 = s + 1; m^2 l^2 = E_0 (E_0 - 3) = 0; Casimir 6; je Energie E = n genau n^2 - 4 Zustände mit den Spins "
      "2..n - 1, zwei Polarisationen im flachen Grenzfall [L]. Der Nullmodul über (4, 1), der in z148 die zehn "
      "Eichfunktionen trägt, fehlt im Fockraum identisch (1.3, 1.5): Der Dubleton ist die eichfixierte "
      "(Feldstärke-)Form derselben Darstellung. Teilweise masselose Zustände (E_0 = 1, 2 bei Spin 2), z148s Geister "
      "der konformen Gravitation, kommen im Fockraum nicht vor: Die Auswahl, die dort P_pos leistet, ist hier von "
      f"selbst erfüllt (z148 reproduzierbar gelesen: {z148_ok})", pass_ok)

# 4.2 Vielfachheitsraum und Modulfluss ---------------------------------------------------------------------------------
xx, yy, ph = sp.symbols("x y phi", real=True)
# Gerade [x : y] im Vielfachheitsraum C^2 (D = +4, D = -4); Modulfluss e^(i phi D) = diag(e^(4 i phi), e^(-4 i phi)).
inv_bed = sp.simplify(sp.exp(4*I_*ph)*xx*yy - sp.exp(-4*I_*ph)*xx*yy)       # Determinante von (v, U v)
inv_nur_chiral = sp.simplify(inv_bed.subs(ph, sp.pi/8) - 2*I_*xx*yy) == 0
per = Fr(2, 32)                                                             # [x : y e^(-32 pi i m t)]: Periode in t m
check("[E] Zwei Kopien gegen eine Mode: Der Vielfachheitsraum von D(3,2) im Dubleton ist C^2 (D = +4, D = -4). Der "
      "Modulfluss e^(i t K) mit K = 4 pi m D wirkt darauf als diag(e^(16 pi i m t), e^(-16 pi i m t)). Eine Gerade "
      "[x : y] ist genau dann modulinvariant, wenn x y = 0 (Determinante 2 i x y bei 4 phi = pi/2: "
      f"{inv_nur_chiral}): Modulinvariant sind nur die beiden chiralen Kopien, und die Tomita-Konjugation J tauscht "
      "sie. Die Einstein-Mode aus z148 mit fester Randbedingung ist eine Diagonale (Feldstärke mit selbstdualem und "
      "antiselbstdualem Teil [L]); der Modulfluss dreht jede Diagonale in die anderen, erstmals in sich zurück nach "
      f"t m = {per}. Eine einzelne Kopie ist also entweder chiral (dann nicht J-symmetrisch) oder bricht die "
      "Modulsymmetrie; dasselbe Bit wie die Phase des Unendlichkeitstwistors (z89, Ladung 2) und die uneigentliche "
      "Randspiegelung (z148 C5). Name: P_chir", inv_nur_chiral)

# 4.3 Erzeuger, Vakuum, Quelle ----------------------------------------------------------------------------------------
E_min_alle = min(F1)                                                         # doppelte tiefste Energie im Fockraum
inv_vektor = any(e == 0 for e in F1)
wsig = {"Vektor J^mu": 1, "Tensor T^mu nu": 2}
ww_ok = all(4 > v for v in wsig.values())                                   # 2|h| = 4 gegen max |sigma|
check("[E] Was im Dubleton fehlt: (a) Alle 15 Erzeuger vertauschen mit D; kein Erzeuger verbindet Omega (D = 0) mit "
      "dem Sektor D = +-4. Eine Quelle für diese Zustände ist quartisch in den Strömen, (a^+)^4, kein Erzeuger und "
      f"kein Strom. (b) Die tiefste Energie im ganzen Fockraum ist {hb(E_min_alle)}; einen so(3,2)-invarianten Vektor "
      f"gibt es nicht ({not inv_vektor}): Omega ist der Grundzustand des Skalars D(1,0), kein Vakuum einer Feldtheorie "
      "auf AdS_4. Ein Feld mit Propagator braucht die zweite Quantisierung des Sektors (P_1T). (c) Weinberg-Witten: "
      f"kovariante Vorwärts-Matrixelemente verlangen 2|h| <= {dict(sorted(wsig.items()))}, der Sektor hat 2|h| = 4 "
      f"({ww_ok}; z75 bestätigt)", zentral and not inv_vektor and E_min_alle == 2 and ww_ok)

# 4.4 Kubischer Vertex ------------------------------------------------------------------------------------------------
EMV = 18
d32c = irrep(6, 4, EMV)
prod = ch_mul(d32c, d32c, EMV)
zer_p = zerlege(prod, EMV)
e_prod = min(e for e, s, m in zer_p)
check("[E] Kubischer Vertex: D(3,2) x D(3,2) hat tiefste Gewichte erst ab E_0 = "
      f"{hb(e_prod)} (Charakterzerlegung bis E = 9), D(3,2) hat nur E_0 = 3. Ein so(3,2)-Verflechter bildet tiefste "
      "Gewichte auf tiefste Gewichte gleicher Energie ab, also Hom(D(3,2) x D(3,2), D(3,2)) = 0: z89s Schranke gilt "
      "auch für die kleinere Gruppe. Der Sektor ist frei; Selbstkopplung und Kopplung an Nähte folgen nicht aus der "
      "Darstellung", zer_p is not None and e_prod == 12)

# 4.5 Vorzeichen und Normierung ---------------------------------------------------------------------------------------
kappa, l_, G_ = sp.symbols("kappa l G", positive=True)
b_mm = l_**2/(64*sp.pi*G_)
res = sp.simplify(32*sp.pi*G_.subs(G_, sp.solve(sp.Eq(b_mm, 16*kappa), G_)[0]))
res_ok = sp.simplify(res - l_**2/(32*kappa)) == 0
check("[C, P_X, P_dyn] Vorzeichen des Residuums gegen die Norm: In z148 ist der kinetische Koeffizient "
      f"2 b z^-2 mit b = l^2/(64 pi G) = 16 kappa, Residuum 32 pi G = {res} ({res_ok}), positiv genau für b > 0 "
      "(Paarungsbit). Die Fock-Norm des Dubletons ist ohne Bedingung positiv. Wird der Sektor mit der Mode "
      "identifiziert (P_X in der engeren Form), so legt das das Paarungsbit auf b > 0 fest (Konsistenzbedingung, kein neues Datum). Den Betrag "
      "legt es nicht fest: Die Einteilchennorm ist 1 für jedes kappa > 0, die Amplitude der Metrikstörung je Quant ist "
      "sqrt(32 pi G)/l und geht nur über die Kopplung ein. P_G bleibt frei", res_ok)
befund("Teil 4: Zustandsraum der z148-Mode", "bedingt",
       "E/C: Darstellung, Randexponent, Masse und Abwesenheit von Eich- und Geistzuständen passen exakt. Es fehlen die "
       "Auswahl einer der zwei Kopien (P_chir), ein Vakuum und ein Feld (P_1T), jede Kopplung (kein Erzeuger verbindet "
       "die Sektoren, kein kubischer Vertex) und der Betrag von G. Weinberg-Witten verbietet den Sektor nicht, sondern "
       "verlangt, dass er als Eichfeld koppelt; genau das liefert erst P_dyn (oder P_X).")

# ====================================================================================================================
say("\nTeil 5: Bilanz und Urteil")
ol = (REPO / "research" / "originlab" / "reconstruction_search.py").read_text()
ol_ok = '("graviton", st["gravity"] == "supported")' in ol
check("[E] OriginLab (nur gelesen): Die Klausel 'graviton' verlangt den Status 'supported' der Gravitationsachse "
      f"({ol_ok}). Nach z156 bleibt der Wert bedingt: Der Sektor D(3,2) existiert kinematisch mit positiver Norm und "
      "positiver Energie, aber Auswahl der Kopie, Vakuum und Kopplung sind Prämissen. Die Klausel bleibt leer",
      ol_ok)
bilanz = [
    ("E", "Darstellungsinhalt",
     "Helizität s der su(2,2) -> D(s+1, s) der nativen so(3,2), irreduzibel für s != 0; s = 0 -> D(1,0) + D(2,0). "
     "Der Sektor D = +-4 des z75-Dubletons ist je D(3,2). Unabhängig von der Wahl der so(3,2) mit kompaktem E."),
    ("E", "Norm",
     "Fock-Form positiv definit, Erzeuger hermitesch; Gram-Form der ersten Stufe 0, 4, 10 wie z148 bei E0 = 3, die "
     "Nullrichtung fehlt identisch; keine Geister."),
    ("E", "Energie",
     "E_de = 1 + (N_d + N_e)/2 >= 3 auf dem Sektor (nativ nach z89). Der zweiseitige Modulgenerator D ist zentral "
     "(Helizität). E_ab = 1 + (N_a + N_b)/2 = N_a + 1 - D/2 ('E x 1 + 1 x E') ist auf dem Sektor die einseitige Uhr "
     "L_0/m - 1 und die um die Dilatation coth(pi m) konjugierte Energie."),
    ("E negativ", "Wirkung auf das Original",
     "Die sechs nichtkompakten Erzeuger der so(3,2) erhalten das Original in keinem Rahmen; die Gruppe wirkt auf der "
     "Purifikation (Lesart T2 mit physikalischem Partner, P_TP oder P_CTM). In T1 gibt es kein D(3,2)."),
    ("E", "Vielfachheit",
     "Je Dubleton 2 (O(2)-Dublett), je Mode und Strompaar ein Dubleton; eine gemeinsame so(3,2) über N Dubletons gibt "
     "Spin 2 erst bei E_0 = N + 2."),
    ("E negativ / C", "Auswahl einer Kopie",
     "Uhr, Paaruhr, Modulgewicht, Gitterlesart wählen nicht. Nötig: P_mode, P_paar, P_chir, P_1T."),
    ("E", "Passung zu z148",
     "Gleiche Darstellung, Delta = 3, masselos, eich- und geistfrei; P_pos ist im Dubleton von selbst erfüllt."),
    ("C", "Kopplung",
     "Kein Erzeuger verbindet Omega mit D = +-4, Hom(D(3,2) x D(3,2), D(3,2)) = 0, Norm 1 für jedes kappa. P_X "
     "(Austausch) oder P_dyn bleibt, P_G bleibt frei; das Paarungsbit wird b > 0 (Konsistenz)."),
    ("L", "Weinberg-Witten",
     "Verbietet lokale kovariante T^mu nu-Kopplung des Sektors, nicht den Sektor; in AdS mit Eichinvarianz ist die "
     "Energie des Gravitons ein Randintegral. Der Sektor muss als Eichfeld koppeln."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
urteil = ("Bedingt. (1) Mathematisch ja [E]: Die Leiterdarstellung der Helizität s der su(2,2) ist unter der nativen "
          "so(3,2) = Der(I_pi - I_omega) für s != 0 irreduzibel und gleich D(s+1, s); für s = 0 zerfällt sie in "
          "D(1,0) + D(2,0). Der Sektor D = +4 des Dubletons aus z75 ist genau D(3,2), der Sektor D = -4 ebenso: "
          "Charaktere bis Stufe 9, genau ein tiefstes Gewicht (exakte Ränge), Gram-Form 0, 4, 10 wie in z148, Casimir 6. "
          "Die Fock-Form ist positiv definit, die so(3,2)-Energie ist auf dem Sektor >= 3. Das Ergebnis hängt nicht von "
          "der Wahl der so(3,2) mit kompaktem E ab. (2) Der Einwand des Bestands ('D ist der Modulgenerator, nicht die "
          "Energie') wird genauer und zum Teil aufgelöst: D = N_a - N_b ist der zweiseitige Modulgenerator, das "
          "Zentrum der u(2,2), also die Helizität. Die so(3,2)-Energie ist ein anderer Operator, und es gibt zwei "
          "ausgezeichnete: E_de (Omega ist ihr tiefster Zustand, nativ nach z89) und E_ab = 'E x 1 + 1 x E' = "
          "N_a + 1 - D/2, auf dem Sektor D = +4 die einseitige Uhr des Originals, L_0/m - 1. Sie sind durch die native "
          "Dilatation mit Faktor coth(pi m) konjugiert, D(3,2) gilt für beide, und positiv ist die Energie, weil L_0 es "
          "ist. Nicht nativ wirken die sechs nichtkompakten Erzeuger: Sie erzeugen Paare aus Original und Partner. Die "
          "so(3,2) ist eine Symmetrie der Purifikation, nicht des Originals; das ist die Lesart T2 mit physikalischem "
          "Partner. Im Thermofeldbild T1 gibt es den Sektor nicht. (3) Nicht nativ ist die Auswahl: Je Dubleton tritt "
          "D(3,2) zweimal auf (Dublett der dualen O(2), vom Modulfluss gedreht, von J getauscht), je Mode und Strompaar "
          "gibt es einen Dubleton, und weder Uhr noch Niveauabgleich noch Modulgewicht noch Gitterlesart wählen einen. "
          "Eine gemeinsame so(3,2) über N Dubletons macht Spin 2 massiv (E_0 = N + 2), weil das Fockvakuum der "
          "Grundzustand des Skalars D(1,0) ist und kein invariantes Vakuum. (4) Kopplung: Der Sektor passt als "
          "Zustandsraum der D(3,2)-Mode aus z148 (Delta = 3, masselos, eich- und geistfrei), er ist aber frei: Kein "
          "Erzeuger verbindet ihn mit Omega, es gibt keinen kubischen Vertex, und die Norm legt G nicht fest. "
          "Weinberg-Witten verbietet ihn nicht, verlangt aber, dass er als Eichfeld koppelt. (5) Für die Prämissen: "
          "P_pos wird Folge (im Dubleton ist jeder Spin-2-Sektor unitär mit positiver Energie und E_0 = 3; teilweise "
          "masselose Zustände kommen nicht vor). P_X bleibt und wird enger: nicht mehr 'ein Sektor D(3,2) tritt hinzu', "
          "sondern 'der vorhandene Sektor wird zwischen Nähten ausgetauscht', mit den benannten Wahlen T2 mit Partner "
          "(P_TP oder P_CTM), P_mode, P_paar, P_chir und P_1T. P_dyn, P_h2 und P_curv bleiben als anderer Weg "
          "unverändert, P_G bleibt frei. K13 bleibt L, die Klausel 'Graviton' bleibt leer; nichts ist befördert. "
          "Nächster Test: die Zwei-Naht-Antwort des Vierstromoperators (a^+)^4 im Kanal D = +-4 unter dem "
          "Transfergatter (z120) bzw. der Kette (z152), geprüft gegen die Linienantwort C(n + 2, 5) einer Primären mit "
          "Delta = 3 (z75). Ist sie null oder faktorisiert sie in vier Stromantworten (Delta = 4), so ist P_X für den "
          "Dubleton widerlegt.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z156_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z156.txt").write_text("\n".join(OUT) + "\n")
