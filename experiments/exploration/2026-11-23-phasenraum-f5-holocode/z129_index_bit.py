"""Z129: Folgeschritt zu V2 -- fixiert der Nahtindex +1 von A_q das Verklebungsbit Ad Z_SU = Ad R (P_IF = P_glue)?

Leitlinie 'Rekonstruktion statt Auswahl': Ein Bit gilt nur als nativ fixiert, wenn vorhandene TFPT-Struktur den anderen Wert
ausschließt. Firewalls: DIMENSION.UPLIFT.FIREWALL.01 (Clifford-, Twistor- und Toeplitzkinematik einer Naht, keine
wechselwirkende 3+1-Schließung), GRAV.NONCIRCULAR.01 (die Uhr nur als Gewicht bzw. Symmetrie). H6 und die gelösten Bausteine
(Flavour, Massen, alpha, die 16, E8) bleiben fest; nichts wird befördert.
Bestand: z116 (A_q = S^+(1 - q S^4), q = e^-2pi, Index +1, Nullmode sqrt(1 - q^2) sum q^n |4n>, D_perp = gamma_- A +
gamma_+ A^+), z117 (Phasenraum W = T + T* der Quarteroszillatoren aus Fockoperatoren, Gamma(Q) = [[0, -*Q], [Q, 0]], tau_W,
J_H6 = J_W (x) J_F K, Schnitt P_0 = gamma_+ gamma_- = (1 + tau_W)/2, Verklebung I_T = I_F <=> Gamma_tot = +1, Spiegel,
Uhrspur), z118 (R A_q R^-1 = -i A_q, D^2 bis auf das Zuordnungsbit, Wick-Kompensator e^(i pi Gamma_2/4), Vorzeichenbits der
MM-Wirkung), z122 (P_IF: Ad Z_SU = Ad R, X = R i^(-2h) = 1; Blatttausch Theta m <-> mquer hält Moden, G, R, E, P_0 und
tauscht die Paarungen; Ladungsenergie seitenlokal nativ 1/4, falsch 3/4), z48 (Marken k = 4(r + alpha) sind die Moden der
Quarter, R|k> = i^k |k>).
Bezeichnungen: T = span{d_i, e_i^+} senkt 2h = N_d - N_e = m - mquer (z85/z89/z122), T* = span{d_i^+, e_i} hebt es;
Gamma_7 = +1 auf T; 16 = Lambda^gerade C^5 mit R_F = i Gamma_11 = i; Gamma_tot = Gamma_7 (x) Gamma_11 (S^+-: Gamma_tot = +-1);
gamma_5 = Gamma_7 tau_W (4D-Chiralität, z117); Gamma_311 = gamma_5 (x) Gamma_11 (Cl(3,11)-Weyl); Schnitt = (1 + tau_W)/2
(pi-Ebene), Komplement = omega-Ebene; Zuordnung a = +1: D_perp = gamma_- A + gamma_+ A^+ (z116), a = -1: gamma_+ A + gamma_- A^+;
Marken: Hardy k >= 0 (L_0 = k/4 >= 0), Anti-Hardy k <= -1. P, C, CP: Blatttausch Theta, innere Spiegelung 1 (x) gamma_1
(16 <-> 16-quer), ihr Produkt.

Teil 1  Nullmode und Schnitt: Nahtindex auf den Marken (Hardy gegen Anti-Hardy, Uhr der Nullmode, Hardy-Umkehr), Lichtkegelpaar
        und Schnitt auf W, D_perp im graduierten Modell (Kern, Schnitt = Hardy x Zuordnung, Wick-Richtung), Sektorinhalt der
        Nullmoden auf W (x) Marken (x) F (tau-Kopie, Cl(3,11), Z4).
Teil 2  Das Verklebungsbit: Blatttausch Theta aus der Fockrealisierung (d <-> e), seine Pin-Form (Spiegelung), die
        entscheidende Kontrolle auf der vollen Indexstruktur, innere Ladungskonjugation, Uhrladung von A_q.
Teil 3  Folgen: Wörterbuch P_IF <-> Gamma_tot (16 gegen 16-quer), MM-Vorzeichen und P_vol, Symmetrietyp des Bits (C, P, CP)
        und was es fixieren kann, Bilanz und Urteil (K16, chiral_matter, graviton).
"""
import itertools
import json
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sps

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
QQ = float(np.exp(-2*np.pi))                                       # q = e^-2pi


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")
    return bool(cond)


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def r0(x, tol=1e-12):
    """Residuum für Ausgabe und JSON: 0 unterhalb der Toleranz (Rundungsrauschen der dichten Algebra), sonst drei Stellen."""
    x = float(x)
    return 0.0 if x < tol else float(f"{x:.3g}")


def zr(x, nd=6):
    """Reelle Zahl gerundet, ohne negative Null (reproduzierbare Ausgabe)."""
    return float(np.round(float(np.real(x)), nd)) + 0.0


def ze(x):
    """Komplexe Zahl gerundet als Text, ohne Vorzeichen der Null."""
    re, im = zr(np.real(x), 6), zr(np.imag(x), 6)
    if im == 0.0:
        return f"{re:g}"
    if re == 0.0:
        return {1.0: "i", -1.0: "-i"}.get(im, f"{im:g}i")
    return f"{re:g}{'+' if im > 0 else '-'}{abs(im):g}i"


def nrm(X):
    """Größter Betrag eines dichten oder dünnen Operators (0 für den Nulloperator)."""
    if sps.issparse(X):
        X = sps.csr_matrix(X)
        return float(abs(X).max()) if X.nnz else 0.0
    return float(np.abs(X).max()) if np.size(X) else 0.0


s0 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]], complex)


def kr(*m):
    out = m[0]
    for x in m[1:]:
        out = np.kron(out, x)
    return out


def euklid(n):
    G = []
    for j in range(n):
        G.append(kr(*([sz]*j + [sx] + [s0]*(n - j - 1))))
        G.append(kr(*([sz]*j + [sy] + [s0]*(n - j - 1))))
    return G


def chiral(G):
    M = np.eye(G[0].shape[0], dtype=complex)
    for g in G:
        M = M @ g
    for ph in (1, 1j, -1, -1j):
        if np.allclose((ph*M) @ (ph*M), np.eye(M.shape[0])):
            return ph*M
    raise ValueError


def nullraum(M, tol=1e-9):
    _, s, vh = np.linalg.svd(np.atleast_2d(M))
    return vh[int(np.sum(s > tol*max(1.0, s[0]))):]


def index_von(A, tol=1e-10):
    r = int(np.sum(np.linalg.svd(A, compute_uv=False) > tol))
    return A.shape[1] - r, A.shape[0] - r


def basis_von(P, G=None):
    """Orthonormalbasis des Bildes eines hermiteschen Projektors P; mit G darin die G-Eigenbasis (aufsteigend)."""
    w, V = np.linalg.eigh((P + P.conj().T)/2)
    B = V[:, w > 0.5]
    if G is not None:
        g = B.conj().T @ G @ B
        _, U = np.linalg.eigh((g + g.conj().T)/2)
        B = B @ U
    return B


# =================================================================== Teil 1
say("Teil 1: Nullmode und Schnitt -- in welchem Sektor liegt die Indexnullmode, und legt das Indexvorzeichen ihn fest?")


# --- 1.1 Nahtindex auf den Marken
def A_hardy(N, q=QQ):
    """A_q = S^+(1 - q S^4) = S^+ - q S^3 auf den Hardy-Marken k = 0..N-1 -> 0..N-2 (rechteckig wie z116/z118)."""
    A = np.zeros((N - 1, N))
    for k in range(N):
        if k >= 1:
            A[k - 1, k] += 1.0
        if k + 3 <= N - 2:
            A[k + 3, k] -= q
    return A


def A_anti(K, q=QQ):
    """Dasselbe Symbol S^+ - q S^3 auf der Anti-Hardy-Hälfte: Marken k = -1..-K -> -1..-K-1 (rechteckig)."""
    src, trg = list(range(-1, -K - 1, -1)), list(range(-1, -K - 2, -1))
    ti = {k: i for i, k in enumerate(trg)}
    A = np.zeros((len(trg), len(src)))
    for j, k in enumerate(src):
        A[ti[k - 1], j] += 1.0
        if k + 3 <= -1:
            A[ti[k + 3], j] -= q
    return A


HARDY, ANTI = {}, {}
for N in (13, 25, 41, 65):
    A = A_hardy(N)
    kk, cc = index_von(A)
    _, _, vh = np.linalg.svd(A)
    z = vh[-1].real*np.sign(vh[-1].real[0])
    soll = np.array([np.sqrt(1 - QQ**2)*QQ**(k//4) if k % 4 == 0 else 0.0 for k in range(N)])
    Rs, Rt = np.diag([1j**k for k in range(N)]), np.diag([1j**k for k in range(N - 1)])
    HARDY[N] = {"Index": kk - cc, "Kern": kk, "Nullmode - Formel": r0(np.abs(z - soll).max(), 1e-10),
                "R A R^-1 + i A": r0(nrm(Rt @ A @ np.linalg.inv(Rs) + 1j*A)), "R z - z": r0(np.abs(Rs @ z - z).max()),
                "Träger mod 4": sorted({k % 4 for k in range(N) if abs(z[k]) > 1e-13})}
for K in (12, 24, 40, 64):
    kk, cc = index_von(A_anti(K))
    ANTI[K] = {"Index": kk - cc, "Kokern": cc}
KL = 20                                                           # zweiseitige Laurent-Abschneidung |k| <= 20
KV = np.arange(-KL, KL + 1)
R2, L02 = np.diag(1j**(KV % 4)), np.diag(KV/4)
PI = np.zeros((len(KV), len(KV)))
for i, k in enumerate(KV):
    PI[int(np.nonzero(KV == -k)[0][0]), i] = 1.0
UMKEHR = {"Pi R Pi^-1 - R^-1": r0(nrm(PI @ R2 @ PI.T - np.linalg.inv(R2))), "Pi R Pi^-1 - R": r0(nrm(PI @ R2 @ PI.T - R2)),
          "Pi L0 Pi^-1 + L0": r0(nrm(PI @ L02 @ PI.T + L02)),
          "Pi P_Hardy Pi^-1 = P_(k <= 0)": bool(np.allclose(PI @ np.diag((KV >= 0)*1.0) @ PI.T, np.diag((KV <= 0)*1.0)))}
ZAHLEN["Marken"] = {"Hardy": {str(k): v for k, v in HARDY.items()}, "Anti-Hardy": {str(k): v for k, v in ANTI.items()},
                    "Hardy-Umkehr": UMKEHR}
check("[E] Nahtindex auf den Marken (z116, z118, z48): A_q = S^+ - q S^3 hat auf den Hardy-Marken k >= 0 (L_0 = k/4 >= 0, "
      f"Abschneidungen N = 13 .. 65) Index {sorted({v['Index'] for v in HARDY.values()})} (Kern 1, Kokern 0); die Nullmode ist "
      f"sqrt(1 - q^2) sum q^n |4n> (Abweichung max {max(v['Nullmode - Formel'] for v in HARDY.values()):.0e}), trägt nur Marken "
      f"k mod 4 in {sorted({k for v in HARDY.values() for k in v['Träger mod 4']})} und hat damit die Uhr R z = z (Residuum "
      f"{max(v['R z - z'] for v in HARDY.values()):.0e}): uhrneutral, obwohl R A_q R^-1 = -i A_q exakt gilt (Residuum "
      f"{max(v['R A R^-1 + i A'] for v in HARDY.values()):.0e}). Dasselbe Symbol auf der Anti-Hardy-Hälfte k <= -1 hat Index "
      f"{sorted({v['Index'] for v in ANTI.values()})}. Kontrolle: Die Hardy-Umkehr Pi: k -> -k (zweiseitig, |k| <= {KL}) "
      f"bildet die Hardy-Hälfte auf k <= 0 ab ({UMKEHR['Pi P_Hardy Pi^-1 = P_(k <= 0)']}) und kehrt Uhr und Energie zugleich um "
      f"(Pi R Pi^-1 = R^-1, Residuum {UMKEHR['Pi R Pi^-1 - R^-1']:.0e}, Abstand zu R {UMKEHR['Pi R Pi^-1 - R']:.0f}; "
      f"Pi L_0 Pi^-1 = -L_0): Hardy-Orientierung und Uhrrichtung sind ein Datum (R = e^(2 pi i L_0) mit L_0 >= 0)",
      all(v["Index"] == 1 and v["Kern"] == 1 and v["Nullmode - Formel"] == 0.0 and v["R A R^-1 + i A"] == 0.0
          and v["R z - z"] == 0.0 and v["Träger mod 4"] == [0] for v in HARDY.values())
      and all(v["Index"] == -1 and v["Kokern"] == 1 for v in ANTI.values())
      and UMKEHR["Pi R Pi^-1 - R^-1"] == 0.0 and UMKEHR["Pi R Pi^-1 - R"] > 1 and UMKEHR["Pi L0 Pi^-1 + L0"] == 0.0
      and UMKEHR["Pi P_Hardy Pi^-1 = P_(k <= 0)"])


# --- 1.2 Phasenraum W = T + T* aus den Fockoperatoren (wie z117), Lichtkegelpaar und Schnitt
def fock(nv, k):
    a = sps.diags(np.sqrt(np.arange(1, nv)), 1, format="csr")
    Id = sps.identity(nv, format="csr")
    ops = []
    for j in range(k):
        op = None
        for m in range(k):
            Mx = a if m == j else Id
            op = Mx if op is None else sps.kron(op, Mx, format="csr")
        ops.append(op.tocsr())
    return ops


def dag(A):
    return sps.csr_matrix(A).T.conj().tocsr()


NV = 5
d1, d2, e1, e2 = fock(NV, 4)
OCC = np.array(list(itertools.product(range(NV), repeat=4)))
DIMF = NV**4
SW = sps.csr_matrix((np.ones(DIMF), (np.array([np.ravel_multi_index((o[2], o[3], o[0], o[1]), (NV,)*4) for o in OCC]),
                                     np.arange(DIMF))), shape=(DIMF, DIMF))              # linearer Tausch d <-> e
DOP = sps.diags((OCC[:, 0] + OCC[:, 1] - OCC[:, 2] - OCC[:, 3]).astype(complex)).tocsr()  # 2h = N_d - N_e
EOP = sps.diags((OCC.sum(1) + 2).astype(complex)/2).tocsr()                               # E = (N_d + N_e + 2)/2
LOW = np.where(OCC.sum(1) <= NV - 3)[0]
BAS = [d1, d2, dag(e1), dag(e2), dag(d1), dag(d2), (-e1).tocsr(), (-e2).tocsr()]         # (Z^alpha, Z-quer_alpha), z117
NAMEN = ["d_1", "d_2", "e_1^+", "e_2^+", "d_1^+", "d_2^+", "-e_1", "-e_2"]
A8 = np.array([sps.csc_matrix(X)[:, LOW].toarray().ravel() for X in BAS]).T
OMEGA_F = np.zeros(DIMF)
OMEGA_F[0] = 1.0


def koeff(X):
    y = sps.csc_matrix(X)[:, LOW].toarray().ravel()
    x = np.linalg.lstsq(A8, y, rcond=None)[0]
    return x, float(np.abs(A8 @ x - y).max())


def matrix_von(F):
    M = np.zeros((8, 8), complex)
    res = 0.0
    for k in range(8):
        x, r = koeff(F(BAS[k]))
        M[:, k] = x
        res = max(res, r)
    return M, res


TH, r_th = matrix_von(lambda X: (SW @ X @ SW).tocsr())           # Blatttausch Theta auf W (linear, z122: m <-> mquer)
MJ, r_j = matrix_von(lambda X: (SW @ X.conj() @ SW).tocsr())     # Tomita J_W (antilinear, z117)
MD, r_d = matrix_von(dag)                                         # Adjungierte
AD, r_a = matrix_von(lambda X: (DOP @ X - X @ DOP).tocsr())       # [2h, .]
AE, r_e = matrix_von(lambda X: (EOP @ X - X @ EOP).tocsr())       # [E, .]
I8, I32 = np.eye(8), np.eye(32)
S4 = np.block([[np.zeros((2, 2)), np.eye(2)], [np.eye(2), np.zeros((2, 2))]]).astype(complex)
G4 = np.diag([1., 1., -1., -1.]).astype(complex)
Z44 = np.zeros((4, 4), complex)
TAUW = MD @ MJ.conj()
GAM7 = np.diag([1.]*4 + [-1.]*4).astype(complex)
LC = np.zeros((4,)*4)
for _p in itertools.permutations(range(4)):
    LC[_p] = round(np.linalg.det(np.eye(4)[list(_p)]))


def stern(Q):
    return 0.5*np.einsum("abcd,cd->ab", LC, Q)


def Gam(Q):
    """Gamma(Q) auf W = T + T*: T -> T* über Q, T* -> T über -*Q (z117); Gamma(Q)^2 = Pf(Q)."""
    return np.block([[Z44, -stern(Q)], [Q, Z44]])


def qf(X, Y):
    return complex(((Gam(X) @ Gam(Y) + Gam(Y) @ Gam(X))/2)[0, 0])


u_om = [np.array([1, 0, -1, 0.]), np.array([0, 1, 0, -1.])]
v_pi = [np.array([1, 0, 1, 0.]), np.array([0, 1, 0, 1.])]
I_PI = (np.outer(u_om[0], u_om[1]) - np.outer(u_om[1], u_om[0])).astype(complex)       # z89/z117
I_OM = (np.outer(v_pi[0], v_pi[1]) - np.outer(v_pi[1], v_pi[0])).astype(complex)
GP, GM = Gam(I_PI)/2, Gam(I_OM)/2                                                      # gamma_+, gamma_-
P_CUT = GP @ GM                                                                         # Schnitt (z117: P_0)
G5W = GAM7 @ TAUW                                                                       # gamma_5 (z117)
Gi = euklid(5)
G11 = chiral(Gi)
a5 = [(Gi[2*j] + 1j*Gi[2*j + 1])/2 for j in range(5)]
PAR = np.eye(32, dtype=complex)
for _j in range(5):
    PAR = PAR @ (np.eye(32) - 2*a5[_j].conj().T @ a5[_j])
G11 = G11*np.sign((G11 @ PAR)[0, 0].real)                                              # +1 auf Lambda^gerade = 16
R_F = 1j*G11
JF = np.eye(32, dtype=complex)
for _j in range(5):
    JF = JF @ Gi[2*_j]
M_H6, M_MW = np.kron(MJ, JF), np.kron(MD, JF)                                           # J_H6, J_MW (antilinear, z117)
null_ok = nrm(GP @ GP) < 1e-12 and nrm(GM @ GM) < 1e-12 and nrm(GP @ GM + GM @ GP - I8) < 1e-12
cut_T = np.linalg.matrix_rank(P_CUT[:4, :4]) == 2 and np.allclose(P_CUT[:4, :4] @ np.array(v_pi).T, np.array(v_pi).T)
jmw = {s_: r0(nrm(np.kron((I8 + s_*TAUW)/2, I32) @ (M_H6 - s_*M_MW) @ np.kron((I8 + s_*TAUW)/2, I32).conj()))
       for s_ in (1, -1)}
vert = {"[gamma_5, gamma_+-]": r0(max(nrm(G5W @ GP - GP @ G5W), nrm(G5W @ GM - GM @ G5W))),
        "{Gamma_7, gamma_+-}": r0(max(nrm(GAM7 @ GP + GP @ GAM7), nrm(GAM7 @ GM + GM @ GAM7))),
        "{tau_W, gamma_+-}": r0(max(nrm(TAUW @ GP + GP @ TAUW), nrm(TAUW @ GM + GM @ TAUW)))}
ZAHLEN["Schnitt"] = {"Residuen Fock": r0(max(r_th, r_j, r_d, r_a, r_e)), "J_H6 -+ J_MW auf der tau-Kopie +-": jmw,
                     "Vertauschungen": vert, "[2h, .] auf W": [zr(x) for x in np.diag(AD)]}
check("[E] Lichtkegelpaar und Schnitt (z117, Fock 5^4, Residuum < 1e-12): gamma_+ = Gamma(I_pi)/2, gamma_- = Gamma(I_omega)/2 "
      f"sind null mit {{gamma_+, gamma_-}} = 1 ({null_ok}); der Schnitt P_0 = gamma_+ gamma_- ist (1 + tau_W)/2 "
      f"({np.allclose(P_CUT, (I8 + TAUW)/2)}), auf T genau die pi-Ebene span{{d_i + e_i^+}} (Rang 2: {cut_T}); auf der "
      f"tau-Kopie P_0 gilt J_H6 = +J_MW, auf dem Komplement J_H6 = -J_MW (Residuen {jmw}). gamma_+- antivertauschen mit "
      "Gamma_7 (T <-> T*) und tau_W (Schnitt <-> Komplement) und vertauschen daher mit gamma_5 = Gamma_7 tau_W "
      f"({vert}); [2h, .] ist -1 auf T und +1 auf T*",
      max(r_th, r_j, r_d, r_a, r_e) < 1e-12 and null_ok and np.allclose(P_CUT, (I8 + TAUW)/2) and cut_T
      and jmw[1] == 0.0 and jmw[-1] == 0.0 and all(v == 0.0 for v in vert.values())
      and np.allclose(np.diag(AD), [-1]*4 + [1]*4))

# --- 1.3 D_perp im graduierten Modell: Kern, Schnitt = Hardy x Zuordnung, Wick-Richtung
BC, BO = basis_von(P_CUT, GAM7), basis_von(I8 - P_CUT, GAM7)                           # Schnitt, Komplement (Gamma_7-Basen)


def D_grad(a, A):
    """D_perp auf (W_q (x) H_s) + (W_z (x) H_t) für A: H_s -> H_t. a = +1 (z116): q = Schnitt, D = gamma_- A + gamma_+ A^+;
    a = -1: q = Komplement, D = gamma_+ A + gamma_- A^+. Rückgabe: D, q-Basis, z-Basis, Dimension des q-Blocks."""
    Bq, Bz = (BC, BO) if a == 1 else (BO, BC)
    gh, gr = (GM, GP) if a == 1 else (GP, GM)
    Gh, Gr = Bz.conj().T @ gh @ Bq, Bq.conj().T @ gr @ Bz
    n1, n2 = 4*A.shape[1], 4*A.shape[0]
    D = np.zeros((n1 + n2, n1 + n2), complex)
    D[n1:, :n1] = np.kron(Gh, A)
    D[:n1, n1:] = np.kron(Gr, A.conj().T)
    return D, Bq, Bz, n1


def nullsektor(a, A):
    """(dim Kern im Schnitt, dim Kern im Komplement) von D_perp auf W (x) Marken."""
    D, _, _, n1 = D_grad(a, A)
    _, s, vh = np.linalg.svd(D)
    K = vh[int(np.sum(s > 1e-10)):].conj().T
    dq = int(round(float(np.real(np.trace(K[:n1].conj().T @ K[:n1])))))
    dz = K.shape[1] - dq
    return (dq, dz) if a == 1 else (dz, dq)


N1 = 13
AH, AA = A_hardy(N1), A_anti(N1 - 1)
SEKT = {(a, lab): nullsektor(a, A) for a in (1, -1) for lab, A in (("Hardy", AH), ("Anti-Hardy", AA))}
D1, _, _, _ = D_grad(1, AH)
D2_soll = sla.block_diag(np.kron(np.eye(4), AH.T @ AH), np.kron(np.eye(4), AH @ AH.T))
herm, d2res = r0(nrm(D1 - D1.conj().T)), r0(nrm(D1 @ D1 - D2_soll))
RS, RT = np.diag([1j**k for k in range(N1)]), np.diag([1j**k for k in range(N1 - 1)])
WICK = {}
for a in (1, -1):
    D, _, _, _ = D_grad(a, AH)
    tq, tz = (1, -1) if a == 1 else (-1, 1)                    # tau_W auf q- und z-Block
    for sg in (1, -1, 0):
        ph_q, ph_z = (np.exp(1j*sg*np.pi/4*tq), np.exp(1j*sg*np.pi/4*tz)) if sg else (1.0, 1.0)
        X = sla.block_diag(ph_q*np.kron(np.eye(4), RS), ph_z*np.kron(np.eye(4), RT))
        WICK[f"a = {a:+d}, " + {1: "e^(+i pi tau/4)", -1: "e^(-i pi tau/4)", 0: "ohne"}[sg]] = r0(nrm(X @ D @ np.linalg.inv(X) - D))
SEKT_TXT = {f"a = {k[0]:+d}, {k[1]}": f"Schnitt {v[0]}, Komplement {v[1]}" for k, v in SEKT.items()}
ZAHLEN["D_perp graduiert"] = {"Nullmoden (Schnitt, Komplement) auf W (x) Marken": SEKT_TXT, "D - D^+": herm,
                              "D^2 - (P_0 A^+A + (1 - P_0) AA^+)": d2res, "Wick-Kovarianz |X D X^-1 - D|": WICK}
check("[E] D_perp im graduierten Modell (Schnitt (x) H_s + Komplement (x) H_t, N = 13; D hermitesch, Residuum "
      f"{herm:.0e}; D^2 = P_0 A^+A + (1 - P_0) AA^+ exakt, Residuum {d2res:.0e}). Kern auf W (x) Marken je (Zuordnung, Marken): "
      f"{SEKT_TXT}. Die Nullmode liegt im Schnitt (pi-Ebene) genau für (a = +1, Hardy) und (a = -1, Anti-Hardy), sonst im "
      "Komplement: Schnitt = Hardy-Orientierung x Zuordnung. Das Indexvorzeichen legt den Lichtkegelsektor (z118s Zuordnungsbit) "
      "also nur zusammen mit der Zuordnung fest. Die Zuordnung ist genau die Wick-Richtung: |X D X^-1 - D| = "
      f"{WICK}; a = +1 (z116) ist nur unter e^(-i pi tau_W/4) (x) R_Marken invariant (= z118s e^(i pi Gamma_2/4) mit "
      "Gamma_2 = -tau_W), a = -1 nur unter e^(+i pi tau_W/4) (x) R_Marken, die Uhr allein bricht beide",
      SEKT[(1, "Hardy")] == (4, 0) and SEKT[(1, "Anti-Hardy")] == (0, 4) and SEKT[(-1, "Hardy")] == (0, 4)
      and SEKT[(-1, "Anti-Hardy")] == (4, 0) and herm == 0.0 and d2res == 0.0
      and WICK["a = +1, e^(-i pi tau/4)"] == 0.0 and WICK["a = -1, e^(+i pi tau/4)"] == 0.0
      and WICK["a = +1, e^(+i pi tau/4)"] > 1 and WICK["a = -1, e^(-i pi tau/4)"] > 1
      and WICK["a = +1, ohne"] > 1 and WICK["a = -1, ohne"] > 1)

# --- 1.4 Sektorinhalt der Nullmoden auf W (x) Marken (x) F
NF = 9
AF = A_hardy(NF)
DW, _, _, _ = D_grad(1, AF)
NW = DW.shape[0]


def gW(X):
    """W-Operator, der Schnitt und Komplement erhält, auf (Schnitt (x) H_s) + (Komplement (x) H_t); dazu das Leck."""
    leck = max(nrm(BO.conj().T @ X @ BC), nrm(BC.conj().T @ X @ BO))
    return sla.block_diag(np.kron(BC.conj().T @ X @ BC, np.eye(NF)), np.kron(BO.conj().T @ X @ BO, np.eye(NF - 1))), leck


G7g, _ = gW(GAM7)
TAUg, _ = gW(TAUW)
DF = np.kron(DW, I32)
DFs = sps.csr_matrix(DF)


def diag_von(X):
    """Diagonale eines (in den gewählten Basen diagonalen) Operators; dazu das Außerdiagonal-Residuum."""
    return np.diag(X).copy(), nrm(X - np.diag(np.diag(X)))


g7v, o1 = diag_von(G7g)
tauv, o2 = diag_von(TAUg)
g11v, o3 = diag_von(G11)
ONE32, ONEW = np.ones(32), np.ones(NW)
G7F, TAUF = np.kron(g7v, ONE32), np.kron(tauv, ONE32)             # Gradierungen als Diagonalvektoren auf W (x) Marken (x) F
GTOT, G311 = np.kron(g7v, g11v), np.kron(g7v*tauv, g11v)
G5F = np.kron(g7v*tauv, ONE32)
RFF = np.kron(ONEW, 1j*g11v)
P16, PT, PCF = np.kron(ONEW, (1 + g11v)/2), (1 + G7F)/2, (1 + TAUF)/2
diag_ok = max(o1, o2, o3) < 1e-12
herm_F = r0(nrm(DF - DF.conj().T))
EW, EV = np.linalg.eigh((DF + DF.conj().T)/2)
KER = EV[:, np.abs(EW) < 1e-9]
KN = np.sum(np.abs(KER)**2, axis=1)                               # Gewicht des Kernprojektors je Basisvektor
spalt = zr(np.min(np.abs(EW[np.abs(EW) >= 1e-9])), 4)


def dk(p):
    """dim(Kern ∩ Bild(p)) für einen mit dem Kernprojektor vertauschenden Diagonalprojektor p."""
    return int(round(float(np.real(p @ KN))))


def ew(p, x):
    """Erwartungswert des Diagonaloperators x auf Kern ∩ Bild(p), normiert (0, falls leer)."""
    d = float(np.real(p @ KN))
    return zr(np.real((p*x) @ KN)/d, 6) if d > 0.5 else 0.0


def spur_k(p, x):
    """Tr_(Kern ∩ Bild(p)) x für Diagonalvektoren p, x."""
    return complex((p*x) @ KN)


RMv = np.concatenate([np.tile([1j**k for k in range(NF)], 4), np.tile([1j**k for k in range(NF - 1)], 4)])
WIv = np.concatenate([np.exp(-1j*np.pi/4)*np.ones(4*NF), np.exp(1j*np.pi/4)*np.ones(4*(NF - 1))])
RMF, WIF = np.kron(RMv, ONE32), np.kron(WIv, ONE32)              # Marken-Uhr, Wick-Element e^(-i pi tau_W/4) (diagonal)
uhr_kern = ze(spur_k(np.ones_like(KN), RMF)/KER.shape[1])        # Mittel unimodularer Zahlen = 1 <=> alle = 1
ZEILEN = {"(T, 16)": PT*P16, "(T, 16q)": PT*(1 - P16), "(T*, 16)": (1 - PT)*P16, "(T*, 16q)": (1 - PT)*(1 - P16)}
TAB = {k: {"dim": dk(p), "Gamma_tot": ew(p, GTOT), "gamma_5": ew(p, G5F), "R_F": ze(spur_k(p, RFF)/max(dk(p), 1)),
           "tau_W": ew(p, TAUF)} for k, p in ZEILEN.items()}
PSP, PSM = (1 + GTOT)/2, (1 - GTOT)/2
dims = {"Kern": KER.shape[1], "Kern im Schnitt": dk(PCF), "S^+": dk(PSP), "S^-": dk(PSM), "16": dk(P16), "16q": dk(1 - P16)}
anti_tot = r0(nrm(DFs @ sps.diags(GTOT) + sps.diags(GTOT) @ DFs))
kom_311 = r0(nrm(DFs @ sps.diags(G311) - sps.diags(G311) @ DFs))
kom_R = r0(nrm(DFs @ sps.diags(RFF) - sps.diags(RFF) @ DFs))
idx_tot = int(round(float(np.real(GTOT @ KN))))
idx_311 = {s_: (dk((1 + s_*G311)/2*PCF), dk((1 + s_*G311)/2*(1 - PCF))) for s_ in (1, -1)}
LEF = {j: (ze(spur_k(PSP, RFF**j)), ze(spur_k(PSM, RFF**j))) for j in range(4)}
heli16 = {"S^+ + S^-": zr(np.real(spur_k(P16, G5F))), "S^+": zr(np.real(spur_k(PSP*P16, G5F))),
          "S^-": zr(np.real(spur_k(PSM*P16, G5F)))}
ZAHLEN["Nullmoden auf W (x) Marken (x) F"] = {"Dimension": DF.shape[0], "D - D^+": herm_F, "Spektrallücke": spalt,
                                             "Dimensionen": dims, "Marken-Uhr auf dem Kern (Mittel)": uhr_kern,
                                             "Zeilen": TAB, "{D, Gamma_tot}": anti_tot,
                                             "[D, Gamma_311]": kom_311, "[D, R_F]": kom_R, "Tr_Kern Gamma_tot": idx_tot,
                                             "Gamma_311 = +-1: (Schnitt, Komplement)": {str(k): v for k, v in idx_311.items()},
                                             "Lefschetz Tr R_F^j (S^+, S^-)": {str(j): v for j, v in LEF.items()},
                                             "Netto-gamma_5 der 16-Nullmoden": heli16}
TAB_TXT = {k: f"dim {v['dim']}, Gamma_tot {v['Gamma_tot']:+g}, gamma_5 {v['gamma_5']:+g}, R_F {v['R_F']}" for k, v in TAB.items()}
check(f"[E, negativ] Sektorinhalt der Nullmoden auf W (x) Marken (x) F (N = 9, Dimension {DF.shape[0]}, D hermitesch, "
      f"Spektrallücke {spalt}): Kern {dims['Kern']} = Schnitt (x) Nullmode (x) F, davon im Schnitt {dims['Kern im Schnitt']} "
      f"(tau-Kopie mit J_H6 = +J_MW), Marken-Uhr auf dem ganzen Kern {uhr_kern}. Die vier Zeilen: {TAB_TXT}. D_perp "
      "antivertauscht mit Gamma_tot "
      f"({anti_tot:.0e}) und vertauscht mit Gamma_311 = gamma_5 Gamma_11 ({kom_311:.0e}) und R_F ({kom_R:.0e}); der Index nach "
      f"Gamma_tot ist Tr_Kern Gamma_tot = {idx_tot}, der tau-Index je Cl(3,11)-Weyl-Sektor (Schnitt, Komplement) "
      f"{idx_311[1]} und {idx_311[-1]}; Lefschetzzahlen Tr R_F^j auf Kern ∩ S^+ und Kern ∩ S^-: {LEF}; Netto-gamma_5 der "
      f"16-Nullmoden {heli16}. Die Indexnullmode wählt keinen Cl(3,11)-Sektor: R = i paart sich mit T (S^+) und mit T* (S^-) "
      "gleich oft, und alle uhrartigen Spuren sind auf S^+ und S^- gleich",
      diag_ok and herm_F == 0.0 and uhr_kern == "1"
      and dims == {"Kern": 128, "Kern im Schnitt": 128, "S^+": 64, "S^-": 64, "16": 64, "16q": 64}
      and all(v["dim"] == 32 for v in TAB.values()) and TAB["(T, 16)"]["Gamma_tot"] == 1 and TAB["(T*, 16q)"]["Gamma_tot"] == 1
      and TAB["(T, 16q)"]["Gamma_tot"] == -1 and TAB["(T*, 16)"]["Gamma_tot"] == -1
      and TAB["(T, 16)"]["gamma_5"] == 1 and TAB["(T*, 16)"]["gamma_5"] == -1 and all(v["tau_W"] == 1 for v in TAB.values())
      and anti_tot == 0.0 and kom_311 == 0.0 and kom_R == 0.0 and idx_tot == 0 and idx_311 == {1: (64, 0), -1: (64, 0)}
      and all(v[0] == v[1] for v in LEF.values()) and LEF[0] == ("64", "64") and heli16 == {"S^+ + S^-": 0.0, "S^+": 32.0,
                                                                                            "S^-": -32.0})
befund("Teil 1: Die Indexnullmode liegt in der tau-Kopie (1 + tau_W)/2 (J_H6 = +J_MW), auf der pi-Ebene, mit Marken-Uhr 1 "
       "(uhrneutral, k = 0 mod 4). Welchen Lichtkegelsektor sie trägt, ist Hardy-Orientierung x Zuordnung; das Indexvorzeichen "
       "fixiert z118s Zuordnungsbit nur zusammen mit der Zuordnung, und die Zuordnung ist die Wick-Richtung e^(-+i pi tau_W/4). "
       "Den Cl(3,11)-Sektor wählt sie nicht: Ihr Multiplett enthält (T, 16), (T*, 16), (T, 16q), (T*, 16q) zu je 32, ist nach "
       "Gamma_tot indexfrei und auf beiden Cl(3,11)-Weyl-Sektoren gleich (64, 64)", True)

# =================================================================== Teil 2
say("\nTeil 2: Das Verklebungsbit -- der Blatttausch Theta als entscheidende Kontrolle")

# --- 2.1 Theta aus der Fockrealisierung
C_ = TH[4:, :4]


def AM(M):
    return np.block([[-M.T, Z44], [Z44, M]])


u22 = True
for a_, b_ in itertools.product(range(4), repeat=2):
    M = np.zeros((4, 4), complex)
    M[a_, b_] = 1
    u22 &= bool(np.allclose(TH @ AM(M) @ TH, AM(-C_ @ M.T @ np.linalg.inv(C_))))
fock_th = {"SW Omega - Omega": r0(np.abs(SW @ OMEGA_F - OMEGA_F).max()), "[SW, E]": r0(abs(SW @ EOP - EOP @ SW).max()),
           "SW 2h SW + 2h": r0(abs(SW @ DOP @ SW + DOP).max())}
wth = {"Theta^2 - 1": r0(nrm(TH @ TH - I8)), "Theta Gamma_7 Theta + Gamma_7": r0(nrm(TH @ GAM7 @ TH + GAM7)),
       "Theta [2h,.] Theta + [2h,.]": r0(nrm(TH @ AD @ TH + AD)), "Theta [E,.] Theta - [E,.]": r0(nrm(TH @ AE @ TH - AE)),
       "Theta - M_J": r0(nrm(TH - MJ)), "ganzzahlig": bool(np.allclose(TH, TH.real.round()))}
JK = r0(nrm(M_H6 @ np.kron(I8, JF).conj() - np.kron(TH, I32)))             # J_H6 o (K_W (x) J_F K) = Theta (x) 1
z122 = json.loads((HERE / "z122_ergebnis.json").read_text())["zahlen"]["CAR zwei Nähte"]["Theta: Moden, Vertizes, G, R, E, h"]
z122_th = {"Moden": z122[0], "G": z122[2], "R": z122[3], "E": z122[4]}
ZAHLEN["Blatttausch Theta"] = {"Fock": fock_th, "W": wth, "u(2,2)-Normalisator (16 Erzeuger)": u22,
                               "J_H6 (K_W (x) J_F K) - Theta (x) 1": JK, "z122 (CAR): Theta auf Moden, G, R, E": z122_th}
check("[E] Blatttausch Theta (z122: m <-> mquer) in der Fockrealisierung: Der lineare Tausch d <-> e hält das Fockvakuum, "
      f"vertauscht mit E = (N_d + N_e + 2)/2 und kehrt 2h = N_d - N_e um ({fock_th}). Auf W ist Theta ganzzahlig, Theta^2 = 1, "
      "bildet T auf T* ab, kehrt den Helizitätserzeuger um und hält die Energie "
      f"({wth}); es normalisiert die native u(2,2) (Theta A_M Theta = A_(-C M^T C^-1) für alle 16 Erzeuger: {u22}). Als Matrix "
      f"ist Theta = M_J, also J_W = Theta K, und Theta (x) 1 = J_H6 (K_W (x) J_F K) exakt (Residuum {JK:.0e}): Der Blatttausch "
      "ist das Produkt aus Tomita-J und dem dritten z117-Spiegel. Auf F und auf den Marken (den Quartermoden, z48) wirkt er als "
      f"Eins; z122 (CAR-Modell): Theta hält Moden, G, R und E exakt ({z122_th})",
      all(v == 0.0 for v in fock_th.values()) and all(v == 0.0 for k, v in wth.items() if k != "ganzzahlig") and wth["ganzzahlig"]
      and u22 and JK == 0.0 and all(v == 0.0 for v in z122_th.values()))

# --- 2.2 Pin-Form: Theta ist eine Spiegelung
MTH7 = TH @ GAM7
Q_C = MTH7[4:, :4]
ist_gam = r0(nrm(MTH7 - Gam(Q_C)))
pf_C = ze((Gam(Q_C) @ Gam(Q_C))[0, 0])
EBAS = []
for _i, _j in itertools.combinations(range(4), 2):
    _X = np.zeros((4, 4), complex)
    _X[_i, _j], _X[_j, _i] = 1, -1
    EBAS += [_X, 1j*_X]


def realform(psi):
    rows = [np.concatenate([(np.exp(-1j*psi)*stern(G4 @ X.conj() @ G4) - X).real.ravel(),
                            (np.exp(-1j*psi)*stern(G4 @ X.conj() @ G4) - X).imag.ravel()]) for X in EBAS]
    return [sum(c*X for c, X in zip(v.real, EBAS)) for v in nullraum(np.array(rows).T)]


def tau_teil(V, sign):
    rows = [np.concatenate([(S4 @ X @ S4 - sign*X).real.ravel(), (S4 @ X @ S4 - sign*X).imag.ravel()]) for X in V]
    return [sum(c*X for c, X in zip(v.real, V)) for v in nullraum(np.array(rows).T)]


def in_span(X, basis):
    Am = np.array([np.concatenate([B.real.ravel(), B.imag.ravel()]) for B in basis]).T
    b = np.concatenate([X.real.ravel(), X.imag.ravel()])
    c = np.linalg.lstsq(Am, b, rcond=None)[0]
    return c, float(np.linalg.norm(Am @ c - b))


def orthonorm(V):
    Gm = np.array([[qf(X, Y).real for Y in V] for X in V])
    w, U = np.linalg.eigh(Gm)
    return [(int(np.sign(w[k])), sum(U[j, k]*V[j] for j in range(len(V)))/np.sqrt(abs(w[k]))) for k in range(len(w))]


Qbild = lambda Q: C_ @ stern(Q) @ C_.T                                                  # Theta Gamma(Q) Theta^-1 = Gamma(Qbild)
RFM = {}
for psi in (0.0, np.pi):
    V = realform(psi)
    RFM[psi] = {"V": V, "R31": tau_teil(V, -1), "R11": tau_teil(V, +1)}
kov = {("0" if p == 0 else "pi"): (r0(max(nrm(TH @ Gam(X) @ TH - Gam(Qbild(X))) for X in RFM[p]["V"])),
                                   r0(max(in_span(Qbild(X), RFM[p]["V"])[1] for X in RFM[p]["V"]), 1e-10))
       for p in RFM}
C_in = {"V^0": r0(in_span(Q_C, RFM[0.0]["V"])[1], 1e-10), "R^(3,1) (tau-ungerade)": r0(nrm(S4 @ Q_C @ S4 + Q_C))}
qCC = qf(Q_C, Q_C)
spieg = max(nrm(Qbild(X - qf(X, Q_C)/qCC*Q_C) - (X - qf(X, Q_C)/qCC*Q_C)) for X in RFM[0.0]["V"])
spieg_C = r0(nrm(Qbild(Q_C) + Q_C))
L31 = np.zeros((4, 4))
for b, Xb in enumerate(RFM[0.0]["R31"]):
    L31[:, b] = in_span(Qbild(Xb), RFM[0.0]["R31"])[0]
ev31 = sorted(zr(x) for x in np.linalg.eigvals(L31).real)
det31 = zr(np.linalg.det(L31))
fest11 = r0(max(nrm(Qbild(X) - X) for X in RFM[0.0]["R11"] + [I_PI, I_OM]))
t_vek = [X for s_, X in orthonorm(RFM[0.0]["R31"]) if s_ == 1]                 # zeitartige Einheit (Gamma^2 = +1 bei psi = 0)
orthochron = bool(len(t_vek) == 1 and qf(t_vek[0], Qbild(t_vek[0])).real >= 1 - 1e-9)   # Q(t, Theta t) >= 1: gleicher Zeitkegel
fest = {"Theta gamma_+- Theta - gamma_+-": r0(max(nrm(TH @ GP @ TH - GP), nrm(TH @ GM @ TH - GM))),
        "Theta tau_W Theta - tau_W": r0(nrm(TH @ TAUW @ TH - TAUW)), "Theta P_0 Theta - P_0": r0(nrm(TH @ P_CUT @ TH - P_CUT)),
        "Theta gamma_5 Theta + gamma_5": r0(nrm(TH @ G5W @ TH + G5W))}
ZAHLEN["Pin-Form"] = {"Theta Gamma_7 = Gamma(Q_C)": ist_gam, "Q_C": [[zr(x) for x in row] for row in Q_C.real],
                      "Gamma(Q_C)^2 = Pf": pf_C, "Q_C in": C_in, "Clifford-Kovarianz (psi: Residuum, Bild in V^psi)": kov,
                      "Spiegelung an Q_C (orthogonal fest, Q_C -> -Q_C)": [r0(spieg), spieg_C],
                      "Eigenwerte auf R^(3,1)": ev31, "det auf R^(3,1)": det31, "R^(1,1), I_pi, I_omega fest": fest11,
                      "orthochron (Q(t, Theta t) >= 1)": orthochron, "fest": fest}
check("[E] Theta ist eine Spiegelung: Theta Gamma_7 = Gamma(Q_C) exakt (Residuum "
      f"{ist_gam:.0e}) mit Q_C = d_1 ^ e_1^+ + d_2 ^ e_2^+ (C_ab = Theta-Block T -> T*), also Theta = Gamma(Q_C) Gamma_7. Q_C liegt "
      f"in der Realform V^0 und ist tau-ungerade, also in R^(3,1) ({C_in}), mit Gamma(Q_C)^2 = Pf(Q_C) = {pf_C}: raumartig (bei "
      "psi = 0 hat die Zeit Gamma^2 = +1, z117). Theta ist Clifford-kovariant, Theta Gamma(Q) Theta^-1 = Gamma(C *Q C^T), und "
      f"hält beide Realformen (psi = 0, pi: {kov}); auf den sechs konformen Richtungen ist es die Spiegelung an Q_C (Q_C -> -Q_C "
      f"{spieg_C:.0e}, Orthogonalraum fest {r0(spieg):.0e}). Auf R^(3,1): Eigenwerte {ev31}, det {det31:g}; als Spiegelung an "
      f"einem raumartigen Vektor orthochron (Q(t, Theta t) >= 1 für die zeitartige Einheit t: {orthochron}), also die "
      f"Paritätsklasse; R^(1,1) mit I_pi, I_omega fest ({fest11:.0e}). Daher exakt {fest}: Der Blatttausch ist eine räumliche "
      "Spiegelung, die Lichtkegelpaar, Dilatation und Schnitt festhält, T mit T* vertauscht und gamma_5 umkehrt",
      ist_gam == 0.0 and pf_C == "-1" and C_in["V^0"] == 0.0 and C_in["R^(3,1) (tau-ungerade)"] == 0.0
      and all(v[0] == 0.0 and v[1] == 0.0 for v in kov.values()) and spieg < 1e-12 and spieg_C == 0.0
      and ev31 == [-1.0, 1.0, 1.0, 1.0] and det31 == -1.0 and fest11 == 0.0
      and orthochron and all(v == 0.0 for v in fest.values()))

# --- 2.3 Entscheidende Kontrolle: Theta, innere Spiegelung C und CP auf der vollen Indexstruktur
THg, leck_th = gW(TH)
THF = sps.kron(sps.csr_matrix(THg), sps.identity(32), format="csr")
CF = sps.kron(sps.identity(NW), sps.csr_matrix(Gi[0]), format="csr")   # 1 (x) gamma_1: 16 <-> 16-quer (innere Ladungskonjugation)
ADg, o4 = gW(AD)
adv, o5 = diag_von(ADg)
ADF = np.kron(adv, ONE32)                                        # [2h, .] als Diagonalvektor


def wirk(U, Ui, Y, unitaer):
    """+1: U Y U^-1 = Y, 'inv': = Y^-1 (nur für unitäres Y, vor -1 geprüft: für R_F ist R_F^-1 = -R_F), -1: = -Y, 0: sonst."""
    Z = U @ Y @ Ui
    if nrm(Z - Y) < 1e-10:
        return "+1"
    if unitaer and nrm(Z - Y.conj().T) < 1e-10:
        return "inv"
    if nrm(Z + Y) < 1e-10:
        return "-1"
    return "0"


SPIEGEL = {"P = Theta": (THF, THF), "C = 1 (x) gamma_1": (CF, CF), "CP = Theta C": ((THF @ CF).tocsr(), (CF @ THF).tocsr())}
inv_ok = max(nrm(U @ Ui - sps.identity(U.shape[0])) for U, Ui in SPIEGEL.values()) < 1e-12
dg = lambda v: sps.diags(v, format="csr")
DATEN = {"D_perp": (DFs, False), "D_perp^2": ((DFs @ DFs).tocsr(), False), "Schnitt": (dg(PCF), False),
         "R_Marken": (dg(RMF), True), "Wick e^(-i pi tau/4)": (dg(WIF), True), "R_F": (dg(RFF), True),
         "[2h, .] (Modulrichtung auf T)": (dg(ADF), True), "Gamma_7": (dg(G7F), True), "Gamma_tot": (dg(GTOT), True),
         "Gamma_311": (dg(G311), True)}
TAFEL = {s_: {d_: wirk(U, Ui, Y, un) for d_, (Y, un) in DATEN.items()} for s_, (U, Ui) in SPIEGEL.items()}


def bildzeile(U, quelle):
    """In welche Zeile bildet U die Kernzeile 'quelle' ab (Anteil > 0,999)?"""
    Kq = KER @ basis_von(KER.conj().T @ (ZEILEN[quelle][:, None]*KER))
    Y = U @ Kq
    gew = np.sum(np.abs(Y)**2, axis=1)
    for z_, p in ZEILEN.items():
        if float(np.real(p @ gew)) > 0.999*Kq.shape[1]:
            return z_
    return "gemischt"


ZB = {s_: {q_: bildzeile(U, q_) for q_ in ZEILEN} for s_, (U, _) in SPIEGEL.items()}
ker_inv = {}
for s_, (U, _) in SPIEGEL.items():
    UK = U @ KER
    ker_inv[s_] = r0(nrm(UK - KER @ (KER.conj().T @ UK)))
GT256 = np.kron(GAM7, G11)
J_kontrolle = {"J_H6 Gamma_tot J_H6^-1 - Gamma_tot": r0(nrm(M_H6 @ GT256.conj() @ np.linalg.inv(M_H6) - GT256)),
               "J_W gamma_+- J_W^-1 - gamma_+-": r0(max(nrm(MJ @ GP.conj() @ np.linalg.inv(MJ) - GP),
                                                        nrm(MJ @ GM.conj() @ np.linalg.inv(MJ) - GM)))}
nicht_trivial = r0(nrm(THF @ dg(GTOT) @ THF - dg(GTOT)))
ZAHLEN["Spiegeltafel der Indexstruktur"] = {"Wirkung": TAFEL, "Zeilenabbildung des Kerns": ZB, "Kern invariant": ker_inv,
                                            "Leck von Theta zwischen Schnitt und Komplement": r0(leck_th),
                                            "Kontrolle J_H6": J_kontrolle, "|Theta Gamma_tot Theta - Gamma_tot|": nicht_trivial}
inv_index = all(TAFEL[s_][d_] == "+1" for s_ in SPIEGEL for d_ in ("D_perp", "D_perp^2", "Schnitt", "R_Marken",
                                                                    "Wick e^(-i pi tau/4)"))
check("[E, negativ] Entscheidende Kontrolle auf der vollen Indexstruktur (W (x) Marken (x) F, N = 9): Theta (Blatttausch, P), "
      "C = 1 (x) gamma_1 (16 <-> 16-quer) und CP = Theta C wirken so (+1 fest, -1 umgekehrt, inv invertiert): "
      f"{TAFEL}. Alle drei halten D_perp, D_perp^2, Schnitt, Marken-Uhr und Wick-Kompensator exakt fest und bilden den Kern auf "
      f"sich ab (Residuen {ker_inv}); Theta und C kehren Gamma_tot und Gamma_311 um, CP hält sie. Theta bildet die Kernzeilen so "
      f"ab: {ZB['P = Theta']} (16 mit T <-> 16 mit T*, gleiche Uhr R_F), C: {ZB['C = 1 (x) gamma_1']}, CP: "
      f"{ZB['CP = Theta C']}. Also ist jede aus Indexdaten gebaute Auswahl Theta-symmetrisch und kann S^+ und S^- nicht "
      "trennen; mit C fällt auch der Ausweg 'Index plus Modulrichtung' (C hält [2h, .] und tauscht S^+-), mit Theta der Ausweg "
      f"'Index plus Uhr' (Theta hält R_F). Kontrollen: Theta ist nichttrivial (|Theta Gamma_tot Theta - Gamma_tot| = "
      f"{nicht_trivial:g}); das Tomita-J_H6 erhält D_perp und Gamma_tot ({J_kontrolle}), ist also kein Tausch",
      inv_index and inv_ok and max(o4, o5) < 1e-12 and all(v == 0.0 for v in ker_inv.values()) and leck_th < 1e-12
      and TAFEL["P = Theta"]["Gamma_tot"] == "-1" and TAFEL["C = 1 (x) gamma_1"]["Gamma_tot"] == "-1"
      and TAFEL["CP = Theta C"]["Gamma_tot"] == "+1" and TAFEL["P = Theta"]["Gamma_311"] == "-1"
      and TAFEL["P = Theta"]["R_F"] == "+1" and TAFEL["C = 1 (x) gamma_1"]["R_F"] == "inv"
      and TAFEL["P = Theta"]["[2h, .] (Modulrichtung auf T)"] == "-1" and TAFEL["C = 1 (x) gamma_1"]["[2h, .] (Modulrichtung auf T)"] == "+1"
      and ZB["P = Theta"] == {"(T, 16)": "(T*, 16)", "(T, 16q)": "(T*, 16q)", "(T*, 16)": "(T, 16)", "(T*, 16q)": "(T, 16q)"}
      and ZB["CP = Theta C"]["(T, 16)"] == "(T*, 16q)" and all(v == 0.0 for v in J_kontrolle.values()) and nicht_trivial == 2.0)

# --- 2.4 Uhrladung von A_q
vR = np.diag(R_F)
adR = sorted({ze(x) for x in np.round(np.outer(vR, 1/vR).ravel(), 12)})
vZ = np.diag(-1j*GAM7)                                           # Z_SU = U_(pi/2) = -i Gamma_7 auf W (z117)
adZ = sorted({ze(x) for x in np.round(np.outer(vZ, 1/vZ).ravel(), 12)})
vW = np.exp(-1j*np.pi/4*np.linalg.eigvalsh(TAUW))
adW = sorted({ze(x) for x in np.round(np.outer(vW, 1/vW).ravel(), 12)})
ad_A = ze(np.vdot(AH.ravel(), (RT @ AH @ np.linalg.inv(RS)).ravel())/np.vdot(AH.ravel(), AH.ravel()))
_, _, vhH = np.linalg.svd(AH)
zH = vhH[-1]
aequi = [ze(np.vdot(zH, np.linalg.matrix_power(RS, j) @ zH)) for j in range(4)]   # Kokern 0
XI = WIF*RMF*RFF                                                                    # kompensierte Gesamtuhr (diagonal)
LEF8 = {j: (ze(spur_k(PSP, XI**j)), ze(spur_k(PSM, XI**j))) for j in range(8)}
ZAHLEN["Uhrladung"] = {"Ad R_F auf End(F)": adR, "Ad Z_SU auf End(W)": adZ, "Ad e^(-i pi tau/4) auf End(W)": adW,
                       "R A_q R^-1 = lam A_q": ad_A, "äquivariante Indexzahlen Tr_Kern R^j - Tr_Kokern R^j": aequi,
                       "Lefschetz Xi^j (S^+, S^-), Xi = Wick (x) R_Marken (x) R_F": {str(j): v for j, v in LEF8.items()}}
check(f"[E, negativ] Uhrladung von A_q (R A_q R^-1 = {ad_A} A_q): Ad R_F hat auf End(F) nur die Eigenwerte {adR}, Ad Z_SU auf "
      f"End(W) nur {adZ}; kein Operator auf der Spinorseite trägt die Ladung +-i (z118: uhrneutral). Der einzige Kompensator ist "
      f"das Wick-Element e^(-i pi tau_W/4) mit Ad-Eigenwerten {adW}, eine achte Einheitswurzel außerhalb der Z4-Uhr; Theta, C "
      "und CP halten ihn fest (Spiegeltafel). Die äquivarianten Indexzahlen von A_q sind Tr_Kern R^j - Tr_Kokern R^j = "
      f"{aequi} für j = 0..3: der triviale Charakter, die Nullmode ist uhrneutral. Mit der kompensierten Gesamtuhr Xi = "
      "e^(-i pi tau/4) (x) R_Marken (x) R_F sind die Lefschetzzahlen Tr Xi^j auf Kern ∩ S^+ und Kern ∩ S^- für j = 0..7 gleich "
      f"({LEF8}): Die Uhrladung bringt keine Information über die Verklebung",
      ad_A == "-i" and adR == ["-1", "1"] and adZ == ["-1", "1"] and "i" in adW and "-i" in adW
      and aequi == ["1", "1", "1", "1"] and all(v[0] == v[1] for v in LEF8.values()))
befund("Teil 2: Der Nahtindex fixiert das Verklebungsbit nicht. Der Blatttausch Theta ist eine räumliche Spiegelung "
       "(Theta = Gamma(Q_C) Gamma_7, Q_C = d_i ^ e_i^+ raumartig): Er hält D_perp, D^2, Nullmode, Schnitt, Marken, beide Uhren "
       "und den Wick-Kompensator exakt fest und vertauscht S^+ mit S^-, also die beiden Paarungen. Er kehrt den Index nicht "
       "um: Die Hardy-Orientierung ist die Energie- und Uhrrichtung der Marken (Hardy-Umkehr hieße R -> R^-1), und Theta hält "
       "beide. Die innere Spiegelung C schließt zusätzlich 'Index plus Modulrichtung' aus. Die Uhrladung -i von A_q ändert "
       "nichts: Ihr Kompensator ist eine Theta-feste achte Wurzel, die Nullmode ist uhrneutral", True)

# =================================================================== Teil 3
say("\nTeil 3: Folgen -- Vorzeichen, MM-Wirkung, P_vol und was das Bit fixieren kann")

# --- 3.1 Wörterbuch P_IF <-> Gamma_tot
erz = {NAMEN[k]: bool(np.abs(BAS[k] @ OMEGA_F).max() > 0.5) for k in range(8)}
erz_T = sorted(n for k, n in enumerate(NAMEN) if erz[n] and k < 4)
erz_Ts = sorted(n for k, n in enumerate(NAMEN) if erz[n] and k >= 4)
XW = np.diag(1j**(-np.diag(AD).real.round().astype(int) % 4))              # i^(-[2h,.]) auf W: Twistorteil von X
U2 = np.diag(np.exp(1j*np.pi/2*np.diag(AD)))                              # U_(pi/2) = Ad Z_SU auf W (z117, z122)
wb = {"i^(-2h) = i Gamma_7 (= z117s I_T)": r0(nrm(XW - 1j*GAM7)),
      "U_(pi/2) = -i Gamma_7 (= z122s I_T = Ad Z_SU)": r0(nrm(U2 + 1j*GAM7)),
      "z117s I_T = (Ad Z_SU)^-1": r0(nrm(1j*GAM7 - np.linalg.inv(U2))),
      "X = i^(-2h) (x) R_F = -Gamma_tot": r0(nrm(np.kron(XW, R_F) + GT256))}
KL_MOD = {}
for z_, (t_, f_) in {"(T, 16)": (1, 1), "(T, 16q)": (1, -1), "(T*, 16)": (-1, 1), "(T*, 16q)": (-1, -1)}.items():
    dR = "i" if f_ == 1 else "-i"
    dh = Fr(-1, 2) if t_ == 1 else Fr(1, 2)
    Xw = (1j if t_ == 1 else -1j)*(1j if f_ == 1 else -1j)
    KL_MOD[z_] = {"Klasse (Delta R, Delta h)": f"({dR}, {'+' if dh > 0 else ''}{dh})", "X": ze(Xw), "Gamma_tot": t_*f_,
                  "z122": ("nativ" if abs(Xw - 1) < 1e-12 else "falsch")}
z122_IT = json.loads((HERE / "z122_ergebnis.json").read_text())["zahlen"]["native Paarung"]["(Twistorphase, Uhrphase) je Typ"]
ad_Tstern = ze(U2[4, 4])                                                   # Ad Z_SU auf T* (hebt h)
pif = (1 - GTOT)/2                                                         # P_IF = (X = 1) = (Gamma_tot = -1)
pif16 = {"dim": dk(pif*P16), "Gamma_7": ew(pif*P16, G7F), "gamma_5": ew(pif*P16, G5F)}
z117_hel = json.loads((HERE / "z117_ergebnis.json").read_text())["zahlen"]["helizitaet"]
ZAHLEN["Wörterbuch"] = {"Erzeuger (B Omega != 0)": {"in T": erz_T, "in T*": erz_Ts}, "Identitäten": wb,
                        "Modulzeilen": KL_MOD, "z122 (I_T, I_F) auf Psi": z122_IT["Psi"], "Ad Z_SU auf T*": ad_Tstern,
                        "R = i-Nullmoden unter P_IF (Schnitt a = +1)": pif16,
                        "z117: gamma_5 der R = +i-Komponente, S^- pi": z117_hel["S- | Schnitt pi (Index +1)"]["gamma_5 der R = +i-Komponente"]}
check("[E] Wörterbuch P_IF <-> Gamma_tot (Vorzeichen 16 gegen 16-quer): Die Erzeuger in T sind "
      f"{erz_T} (Partnerblatt), in T* {erz_Ts} (Originalblatt, hebt 2h). Die native Klasse (R = i, Delta h = +1/2), z122s "
      "Psi = V^+ (x) S mit dem Originalschritt, ist im Modul T* (x) 16 und liegt in Gamma_tot = -1, also in z117s S^-: "
      f"{KL_MOD}. Exakt gilt {wb}: z122s Twistorphase I_T = Ad Z_SU ist auf T* gleich {ad_Tstern} (z122 auf Psi: "
      f"{z122_IT['Psi']}), z117s I_T = i Gamma_7 ist ihr Inverses, und X = R i^(-2h) = -Gamma_tot. Also P_IF <=> X = 1 <=> "
      "Gamma_tot = -1 <=> (z122) I_T = I_F <=> (z117) I_T = -I_F: dasselbe Bit wie P_glue, aber z117s Etikett S^+ ist im "
      "Operatorbild die Spiegelpaarung (in S^+ antwortet die 16 über e^+, z117). Unter P_IF und dem Schnitt von z116 sind die "
      f"R = i-Nullmoden (T* ∩ Schnitt) (x) 16: {pif16} (Twistorhelizität +1/2, z122: N16 = +3; gamma_5 = -1 in z117s Konvention, "
      f"z117: {z117_hel['S- | Schnitt pi (Index +1)']['gamma_5 der R = +i-Komponente']:g})",
      erz_T == ["e_1^+", "e_2^+"] and erz_Ts == ["d_1^+", "d_2^+"] and all(v == 0.0 for v in wb.values())
      and KL_MOD["(T*, 16)"]["z122"] == "nativ" and KL_MOD["(T*, 16)"]["Gamma_tot"] == -1
      and KL_MOD["(T, 16q)"]["z122"] == "nativ" and KL_MOD["(T, 16)"]["z122"] == "falsch" and KL_MOD["(T*, 16q)"]["z122"] == "falsch"
      and z122_IT["Psi"] == [["i", "i"]] and ad_Tstern == "i" and pif16 == {"dim": 32, "Gamma_7": -1.0, "gamma_5": -1.0}
      and z117_hel["S- | Schnitt pi (Index +1)"]["gamma_5 der R = +i-Komponente"] == -1.0)

# --- 3.2 MM-Vorzeichen (Uhrspur) und P_vol
R31 = [Gam(X) for _, X in sorted(orthonorm(RFM[0.0]["R31"]), key=lambda t: -t[0])]    # Zeit (+1) zuerst
if not np.allclose(1j*R31[0] @ R31[1] @ R31[2] @ R31[3], G5W):
    R31[3] = -R31[3]
orient = r0(nrm(1j*R31[0] @ R31[1] @ R31[2] @ R31[3] - G5W))
EPS = np.zeros((4, 4, 4, 4))
for _p in itertools.permutations(range(4)):
    EPS[_p] = round(np.linalg.det(np.eye(4)[list(_p)]))
PAARE = [(a_, b_) for a_ in range(4) for b_ in range(a_ + 1, 4)]
R256 = np.kron(I8, R_F)
PSc = {("S+" if s_ == 1 else "S-"): (np.eye(256) + s_*GT256)/2 for s_ in (1, -1)}
PCc = {"pi": np.kron(P_CUT, I32), "omega": np.kron(I8 - P_CUT, I32)}
TH256 = np.kron(TH, I32)
SPUR, transport = {}, 0.0
for (a_, b_) in PAARE:
    for (c_, d_) in PAARE:
        X = np.kron(R31[a_] @ R31[b_] @ R31[c_] @ R31[d_], I32)
        for sS, PS_ in PSc.items():
            for sC, PC_ in PCc.items():
                SPUR.setdefault((sS, sC), []).append((EPS[a_, b_, c_, d_], complex(np.trace(PS_ @ PC_ @ R256 @ X))))
        for sC, PC_ in PCc.items():
            transport = max(transport, abs(np.trace(PSc["S+"] @ PC_ @ R256 @ X)
                                           - np.trace(PSc["S-"] @ PC_ @ R256 @ TH256 @ X @ TH256)))
cMM = {}
for k_, werte in SPUR.items():
    kon = {zr((w/e).real, 6) for e, w in werte if e != 0}
    rest = max(abs(w) for e, w in werte if e == 0)
    cMM[f"{k_[0]}, {k_[1]}"] = (sorted(kon)[0] if len(kon) == 1 and rest < 1e-9 else None)
regel = all(cMM[f"{sS}, {sC}"] == 64*(1 if sS == "S+" else -1)*(1 if sC == "pi" else -1) for sS in PSc for sC in PCc)
U4 = np.diag(np.exp(1j*np.pi/4*np.diag(AD)))
U4g, _ = gW(U4)
u4v, o6 = diag_von(U4g)
UK4 = np.kron(u4v, ONE32)[:, None]*KER
pvol = {"[U, P_0], [U, tau_W], [U, Gamma_7], [U, gamma_5]": r0(max(nrm(U4 @ Y - Y @ U4) for Y in (P_CUT, TAUW, GAM7, G5W))),
        "U hält den Nullmodenraum": r0(max(nrm(UK4 - KER @ (KER.conj().T @ UK4)), o6)),
        "|U gamma_+ U^-1 - gamma_+| (Kontrolle)": r0(nrm(U4 @ GP @ np.linalg.inv(U4) - GP)),
        "Theta U Theta - U^-1": r0(nrm(TH @ U4 @ TH - np.linalg.inv(U4)))}
ZAHLEN["MM-Vorzeichen und P_vol"] = {"Orientierung i Gamma_0123 = gamma_5": orient, "c (Uhrspur / eps)": cMM,
                                     "c = 64 x Gamma_tot x Schnitt": regel, "Theta-Transport (max)": r0(transport, 1e-9),
                                     "det Theta auf R^(3,1)": det31, "P_vol (U_(pi/4))": pvol}
check("[E] MM-Vorzeichen und P_vol: Mit dem nativ orientierten R^(3,1)-Rahmen (i Gamma_0123 = gamma_5, Residuum "
      f"{orient:.0e}) ist die Uhrspur Tr_(S ∩ Schnitt)(R_F Gamma_a Gamma_b Gamma_c Gamma_d) = c eps_abcd für alle 36 Paare mit "
      f"c = {cMM} = 64 x Gamma_tot x Schnitt ({regel}, Uhr R). Theta trägt exakt Tr_(S^+ ∩ Schnitt)(R X) = "
      f"Tr_(S^- ∩ Schnitt)(R Theta X Theta^-1) (36 Paare, beide Schnitte, Abweichung {r0(transport, 1e-9):.0e}) und kehrt die "
      f"Rahmenorientierung um (det {det31:g}): Das MM-Vorzeichen ist bei fester Orientierung Theta-ungerade, also mit "
      "Gamma_tot das Verklebungsbit; bei gegebener Verklebung ist es Schnitt x Uhr, der Schnitt Hardy x Zuordnung (Teil 1). "
      f"P_vol: U_(pi/4) vertauscht mit Schnitt, tau_W, Gamma_7 und gamma_5 und hält den Nullmodenraum ({pvol}), dreht aber "
      "gamma_+ (Modulphase, z117): Der Index sieht P_vol nicht, und Theta kehrt U_(pi/4) um",
      orient == 0.0 and regel and r0(transport, 1e-9) == 0.0 and det31 == -1.0
      and pvol["[U, P_0], [U, tau_W], [U, Gamma_7], [U, gamma_5]"] == 0.0 and pvol["U hält den Nullmodenraum"] == 0.0
      and pvol["|U gamma_+ U^-1 - gamma_+| (Kontrolle)"] > 0.1 and pvol["Theta U Theta - U^-1"] == 0.0)

# --- 3.3 Symmetrietyp des Bits und was es fixieren kann
SEITEN = ("o", "p")
VERT = [(sc, c, sl) for sc in SEITEN for c in (1, -1) for sl in SEITEN]     # (Ladungsseite, Schritt, Blattseite), z122 Modell B
E_LAD = {c: Fr(c*c, 2) + (Fr(1, 4) - Fr(1, 2))*c for c in (1, -1)}          # Quarter a, alpha = 1/4 (z82, z122)


def klasse_v(v):
    sc, c, sl = v
    dT = c if sc == "o" else -c                                              # R_phys = i^(Q_o - Q_p)
    return f"({'i' if dT == 1 else '-i'}, {'+1/2' if sl == 'o' else '-1/2'})"


NATIV_K, FALSCH_K = ["(-i, -1/2)", "(i, +1/2)"], ["(-i, +1/2)", "(i, -1/2)"]
praedikat = lambda v: v[0] == v[2] and E_LAD[v[1]] == Fr(1, 4)               # seitenlokal und tiefste Ladungsenergie
P_v = lambda v: (v[0], v[1], "p" if v[2] == "o" else "o")                    # Theta: Blatt tauschen
C_v = lambda v: (v[0], -v[1], v[2])                                           # C: Ladungsschritt umkehren
wahl = lambda pr: sorted({klasse_v(v) for v in VERT if pr(v)})
energie = {"seitenlokal nativ": str(E_LAD[1]), "seitenlokal falsch": str(E_LAD[-1])}
WAHL = {"Prädikat": wahl(praedikat), "P-Bild": wahl(lambda v: praedikat(P_v(v))), "C-Bild": wahl(lambda v: praedikat(C_v(v))),
        "CP-Bild": wahl(lambda v: praedikat(C_v(P_v(v))))}
# Theta-ungerades Datum P_IF auf den Nullmoden: chiral; sein Theta-Bild: die Spiegelpaarung
chir = {"P_IF (Gamma_tot = -1)": {"dim": dk(PSM), "16 mit T*": dk(PSM*ZEILEN["(T*, 16)"]), "16 mit T": dk(PSM*ZEILEN["(T, 16)"]),
                                  "Netto-Gamma_7 der 16": zr(np.real(spur_k(PSM*P16, G7F)))},
        "Theta-Bild (Gamma_tot = +1)": {"dim": dk(PSP), "16 mit T*": dk(PSP*ZEILEN["(T*, 16)"]), "16 mit T": dk(PSP*ZEILEN["(T, 16)"]),
                                        "Netto-Gamma_7 der 16": zr(np.real(spur_k(PSP*P16, G7F)))},
        "ohne Datum": {"Netto-Gamma_7 der 16": zr(np.real(spur_k(P16, G7F)))}}
typ_bit = {s_: TAFEL[s_]["Gamma_tot"] for s_ in SPIEGEL}
typ_index = {s_: TAFEL[s_]["D_perp"] for s_ in SPIEGEL}
typ_zuordnung = {s_: TAFEL[s_]["Wick e^(-i pi tau/4)"] for s_ in SPIEGEL}
ZAHLEN["Symmetrietyp"] = {"Verklebungsbit Gamma_tot unter (P, C, CP)": typ_bit, "D_perp unter (P, C, CP)": typ_index,
                          "Zuordnung/Wick unter (P, C, CP)": typ_zuordnung, "seitenlokale Energieordnung (z122)": WAHL,
                          "Ladungsenergien": energie, "Nullmoden mit Theta-ungeradem Datum": chir}
check("[C, Kontrolle] Symmetrietyp und was das Bit fixieren kann: Das Verklebungsbit ist P-ungerade, C-ungerade und "
      f"CP-gerade ({typ_bit}), wie die Chiralität des Standardmodells; D_perp, Nullmode, Schnitt und Zuordnung (Wick) sind unter "
      f"allen drei fest ({typ_index}, {typ_zuordnung}). Ein Datum vom Typ des Bits wählt: z122s seitenlokale Energieordnung "
      f"(seitenlokal und tiefste Ladungsenergie {energie['seitenlokal nativ']} gegen {energie['seitenlokal falsch']}) wählt genau "
      f"{WAHL['Prädikat']} (nativ), ihr P- und C-Bild genau {WAHL['P-Bild']} bzw. {WAHL['C-Bild']} (falsch), ihr CP-Bild "
      f"{WAHL['CP-Bild']}: derselbe Typ, aber nur eine Ordnung (Abstand 1/2 in L_0, e^(-2 pi)), kein Ausschluss. Auf den "
      f"Nullmoden macht das Datum P_IF das Multiplett chiral, sein Theta-Bild die Spiegelpaarung, ohne Datum ist es "
      f"vektorartig: {chir}",
      typ_bit == {"P = Theta": "-1", "C = 1 (x) gamma_1": "-1", "CP = Theta C": "+1"}
      and all(v == "+1" for v in typ_index.values()) and all(v == "+1" for v in typ_zuordnung.values())
      and WAHL["Prädikat"] == NATIV_K and WAHL["P-Bild"] == FALSCH_K and WAHL["C-Bild"] == FALSCH_K and WAHL["CP-Bild"] == NATIV_K
      and chir["P_IF (Gamma_tot = -1)"] == {"dim": 64, "16 mit T*": 32, "16 mit T": 0, "Netto-Gamma_7 der 16": -32.0}
      and chir["Theta-Bild (Gamma_tot = +1)"] == {"dim": 64, "16 mit T*": 0, "16 mit T": 32, "Netto-Gamma_7 der 16": 32.0}
      and chir["ohne Datum"]["Netto-Gamma_7 der 16"] == 0.0)
befund("Teil 3: Das MM-Vorzeichen (bei fester Orientierung) ist mit Gamma_tot das Verklebungsbit und wird vom Index so wenig "
       "fixiert wie dieses; z118s Zuordnungsbit fixiert der Index nur zusammen mit der Wick-Richtung, P_vol gar nicht. Das "
       "fehlende Datum ist ein P-ungerades, C-ungerades, CP-gerades Bit. Minimal ist P_glue selbst (gleichwertig P_IF bzw. "
       "eta_+-Kleidung); in Indexform wäre es die Gleichsetzung der Hardy-Orientierung der Naht mit der Modul- bzw. "
       "Holomorphierichtung von T (dann kehrte Theta den Index um). Der einzige vorhandene Kandidat vom richtigen Typ ist "
       "Seitenlokalität plus Energieordnung alpha = 1/4, die ordnet, aber nicht ausschließt. Im Operatorbild ist P_IF "
       "Gamma_tot = -1 (z117s S^-)", True)

# =================================================================== Bilanz
say("\nBilanz und Urteil")
STATUS = {
    "Nullmode und Schnitt": ("E", "Kern = Schnitt (x) Nullmode (x) F, tau-Kopie mit J_H6 = +J_MW, pi-Ebene, Marken-Uhr 1; "
                                  "Schnitt = Hardy x Zuordnung, Zuordnung = Wick-Richtung"),
    "z118s Zuordnungsbit": ("E (negativ)", "Indexvorzeichen allein legt es nicht fest; mit der z116-Zuordnung ist es pi"),
    "Cl(3,11)-Sektor": ("E (negativ)", "je 32 in (T, 16), (T*, 16), (T, 16q), (T*, 16q); D vertauscht mit Gamma_311, Index 64/64, "
                                       "Gamma_tot-Index 0, Lefschetz gleich"),
    "Blatttausch Theta": ("E", "Theta = Gamma(Q_C) Gamma_7, räumliche Spiegelung (Paritätsklasse); exakte Symmetrie von D_perp, "
                               "Kern, Schnitt, Marken, Uhren, Wick; vertauscht S^+ und S^-"),
    "Uhrladung von A_q": ("E (negativ)", "Ad R_F, Ad Z_SU nur +-1; Kompensator achte Wurzel, Theta-fest; Nullmode uhrneutral"),
    "Verklebungsbit aus dem Index": ("E (negativ)", "nicht fixiert: Index P-, C- und CP-gerade, Bit P- und C-ungerade"),
    "Wörterbuch": ("E", "P_IF = X = 1 = Gamma_tot = -1 (z117s S^-); dasselbe Bit wie P_glue, entgegengesetztes Etikett"),
    "MM-Vorzeichen": ("E", "= 64 x Gamma_tot x Schnitt x Uhr bei nativer Orientierung; Theta-ungerade, vom Index nicht fixiert"),
    "P_vol": ("E", "Nullmodenraum U_(pi/4)-invariant; Index blind für P_vol (Moduleichrichtung, z117)"),
    "Symmetrietyp": ("E/C", "Bit P-ungerade, C-ungerade, CP-gerade; seitenlokale Energieordnung (z122) vom selben Typ, nur Ordnung"),
    "minimales Zusatzdatum": ("C", "ein P- und C-ungerades Bit: P_glue = P_IF; in Indexform Hardy-Orientierung der Naht = "
                                   "Modul-/Holomorphierichtung von T"),
    "K16": ("L", "unverändert: Linienindex 3 je 16 unter T, P_GSO, K2, P_IF; der Nahtindex liefert P_IF nicht"),
    "chiral_matter": ("bedingt", "unverändert auf P_IF (= P_glue als ein Bit) + T + K2 + P_GSO; Nahtindex-Route gescopt negativ"),
    "graviton": ("bedingt", "MM-Vorzeichen trägt Gamma_tot; P_R, 2 c_3 und P_cap unverändert offen"),
    "gesamt": ("gescoptes Negativ", "Nahtindex +1 fixiert weder Verklebung noch (allein) Schnitt; minimal fehlt ein P-ungerades Bit"),
}
for k_, (s2, t2) in STATUS.items():
    say(f"   [{s2}] {k_}: {t2}")
ZAHLEN["status"] = STATUS
ZAHLEN["z129_urteil_kurz"] = {
    "Urteil": "gescoptes Negativ: der Nahtindex +1 fixiert das Verklebungsbit nicht (Theta-Kontrolle)",
    "K16": "bleibt L",
    "chiral_matter": "bedingt auf P_IF + T + K2 (+ P_GSO), unverändert; Nahtindex-Route geschlossen",
    "graviton": "MM-Vorzeichen = Gamma_tot x Schnitt x Uhr; vom Index nicht fixiert; übrige Prämissen unverändert"}
urteil = ("Der Nahtindex +1 von A_q fixiert das Verklebungsbit Ad Z_SU = Ad R (P_IF = P_glue) nicht; das ist ein gescoptes "
          "Negativ für die Indexroute. (1) Nullmode: A_q hat auf den Hardy-Marken Index +1 mit der Nullmode sqrt(1 - q^2) sum "
          "q^n |4n>, die nur Marken k = 0 mod 4 trägt und daher uhrneutral ist (R z = z), obwohl R A_q R^-1 = -i A_q. Im "
          "graduierten Modell ist der Kern von D_perp = gamma_- A + gamma_+ A^+ genau Schnitt (x) Nullmode: die tau-Kopie "
          "(1 + tau_W)/2, in der J_H6 = +J_MW gilt, auf T die pi-Ebene span{d_i + e_i^+}. Welchen Lichtkegelsektor die Nullmode "
          "trägt, ist das Produkt aus Hardy-Orientierung und Zuordnung (welches gamma mit A multipliziert); das Indexvorzeichen "
          "legt z118s Zuordnungsbit also nur zusammen mit der Zuordnung fest, und diese ist genau die Wick-Richtung (D_perp ist "
          "nur unter e^(-i pi tau_W/4) (x) R_Marken invariant, die Gegenzuordnung nur unter e^(+i pi tau_W/4) (x) R_Marken). "
          "Hardy-Orientierung und Uhrrichtung sind ein Datum (die Hardy-Umkehr kehrt R und L_0 zugleich um). (2) Sektor: Auf "
          "W (x) Marken (x) F (2176 Dimensionen) ist der Kern 128-dimensional und zerfällt in (T, 16), (T, 16q), (T*, 16), "
          "(T*, 16q) zu je 32; D_perp antivertauscht mit Gamma_tot, vertauscht mit Gamma_311 und R_F, hat Gamma_tot-Index 0 und "
          "auf beiden Cl(3,11)-Weyl-Sektoren den tau-Index 64; alle Lefschetzzahlen der Uhr sind auf S^+ und S^- gleich. Die "
          "Nullmode wählt keinen Cl(3,11)-Sektor. (3) Entscheidende Kontrolle: Der Blatttausch Theta (d <-> e, z122) ist auf W "
          "die Spiegelung Theta = Gamma(Q_C) Gamma_7 an der raumartigen Richtung Q_C = d_i ^ e_i^+ von R^(3,1), also eine "
          "Paritätstransformation; er hält Lichtkegelpaar, Dilatation, Schnitt, Energie, die Marken (Quartermoden) und beide "
          "Uhren fest, vertauscht T mit T* und damit S^+ mit S^-. Er ist eine exakte Symmetrie von D_perp, D^2, Nullmode, "
          "Schnitt und Wick-Kompensator und bildet (T, 16) auf (T*, 16) ab. Er kehrt den Index nicht um, denn eine Hardy-Umkehr "
          "müsste die Uhr invertieren. Die innere Spiegelung C = 1 (x) gamma_1 ist ebenfalls Symmetrie der Indexstruktur, hält "
          "die Modulrichtung und tauscht S^+-: Weder 'Index plus Uhr' noch 'Index plus Modulrichtung' entscheidet. Die "
          "Uhrladung -i von A_q ändert nichts: Ad R_F und Ad Z_SU haben nur die Eigenwerte +-1, der Kompensator ist eine "
          "Theta-feste achte Wurzel, die äquivarianten Indexzahlen sind trivial. (4) Folgen: Das Bit ist P-ungerade, C-ungerade "
          "und CP-gerade, wie die Chiralität des Standardmodells; die Indexstruktur ist unter P, C und CP fest. Das "
          "MM-Vorzeichen ist bei nativer Orientierung 64 x Gamma_tot x Schnitt x Uhr, Theta-ungerade und damit vom Index nicht "
          "fixiert; P_vol sieht der Index nicht (U_(pi/4) hält den Nullmodenraum). Wörterbuch: Die native Klasse (R = i, Delta "
          "h = +1/2) ist T* (x) 16 und liegt in Gamma_tot = -1; X = R i^(-2h) = -Gamma_tot; z117s I_T = i Gamma_7 ist das "
          "Inverse von z122s I_T = Ad Z_SU. P_IF und P_glue sind dasselbe Bit, aber im Operatorbild ist P_IF z117s S^-; z117s "
          "Etikett S^+ (die 16 antwortet über e^+) ist die Spiegelpaarung. Unter P_IF und dem z116-Schnitt sind die R = i-"
          "Nullmoden (T* ∩ Schnitt) (x) 16 mit Twistorhelizität +1/2 (z122: N16 = +3). Minimal fehlt ein P- und C-ungerades "
          "Bit: P_glue selbst; in Indexform die Gleichsetzung der Hardy-Orientierung der Naht mit der Modul- bzw. "
          "Holomorphierichtung von T, unter der Theta den Index umkehren würde. Der einzige vorhandene Kandidat vom richtigen "
          "Typ ist z122s seitenlokale Energieordnung (alpha = 1/4: nativ 1/4, falsch 3/4), die ordnet, aber nicht ausschließt. "
          "Folge für K16: bleibt L. Folge für chiral_matter: unverändert bedingt auf P_IF + T + K2 (+ P_GSO), die Nahtindex-"
          "Route ist geschlossen. Folge für graviton: Das MM-Vorzeichen trägt dasselbe Bit; P_R, 2 c_3 und P_cap bleiben offen. "
          "Nichts ist befördert; 3+1-Aussagen bleiben Kinematik der Einpaarschicht.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z129_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
