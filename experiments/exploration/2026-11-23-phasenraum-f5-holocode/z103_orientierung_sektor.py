"""Z103: Wählt eine native Orientierung den Gravitationssektor? (Folgeauftrag N4, Hypothese H4)

Leitlinie 'Rekonstruktion statt Auswahl'. Firewalls: GRAV.NONCIRCULAR.01 (der Modulfluss tritt nur als Symmetrie der nativen
Daten auf, nicht als Zeit oder Raumzeit), DIMENSION.UPLIFT.FIREWALL.01 (alles ist Twistor- und Fockkinematik einer bzw. zweier
Nähte, Kurventransport auf der Spinbrücke, Quarterfermionen und Formenalgebra; keine wechselwirkende 3+1-Schließung). Gelöste
Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt; nichts wird befördert. Bestand: z98 (Plebanski aus zwei Blättern:
B_phi = e^(i phi) Sigma_+ + e^(-i phi) Sigma_- = cos phi e^e - sin phi *(e^e) ist genau für phi in (pi/2) Z einfach, vier
Sektoren +-e^e, +-*(e^e); U_(pi/4) verbindet sie; Phase von I = Phase einer Volumenform; Modul-U(1) erhält alle nativen
Daten, I_pi trägt die Ladung 2), z97 (Z_SU = i^(2h - 2)), z95 (J R J = R^-1 auf dem Quadrattorus), z84 (Theta_n als
U(2) x U(2)-Bahn von tau, Familie D invertiert R, Familie A und Tomita-J zentralisieren R, H0(N) = V^* (x) eta), z83
(Ad_R = 1 auf den Schritten), z82 (R = i^(Q_a - Q_b), eta-Dublett, Torus erzeugt von Q_a - Q_b), z79 (R = e^(2 pi i L_0),
Quarteroszillatoren = Twistorraum, auf ihnen R = 1).
Statuskürzel: E exakt, C bedingt, L Lesart, O offen, N numerisch.

Teil 1  Zensus der Orientierungsdaten. Für jedes Datum (Marken mit Uhr R, Hardy-Orientierung, tau und Theta_n, Uhrrichtung
        R gegen R^-1, KMS-Richtung, Signatur (2,2), eta-Dublett, Zentrum Z_SU) die exakte Wirkung auf die Volumenform eps
        des Twistorraums, auf die Phase von I_pi und auf die vier Plebanski-Sektoren (Kern V = eps c^+ c^+ auf zwei Nähten,
        Fock 3^8; Sektor über die Simplizität aus z98) sowie das Verhalten unter U_(pi/4) (vertauscht, invertiert, bricht).
        Uhr auf Quarterfermionen (R = i^(Q_a - Q_b), trivial auf den Strömen) und auf der Spinbrücke H0(N) (Kurventransport).
Teil 2  Auswahl: Gruppe der nativen Sektorwirkungen und ihre invarianten Teilmengen (alle 16), Stabilisator aller nativen
        Strukturen in der Modulgruppe, die rahmenabhängigen Spiegelungen K_theta und C_theta als Modulbahn, Suche nach einem
        U_(pi/4)-ungeraden Orientierungsdatum, die Volumenform eps_psi als Selektor (Paar je Orientierung), minimale Prämisse.
Teil 3  Bilanz, Urteil zu H4 und Folge für K13.
"""
import itertools
import json
import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
import scipy.sparse as sps
import sympy as sp

HERE = Path(__file__).parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "research" / "originlab"))
from reconstruction_search import _transport  # noqa: E402

OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
rng = np.random.default_rng(103)
I2 = np.eye(2, dtype=complex)
PAULI = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.diag([1., -1.]).astype(complex)]
SIGMA = [I2] + PAULI
EPS2 = np.array([[0., 1.], [-1., 0.]], dtype=complex)
ETA = np.diag([1., -1., -1., -1.])
G4 = np.diag([1., 1., -1., -1.]).astype(complex)                        # CCR-Form auf T = span{d_1, d_2, e_1^+, e_2^+}
S4 = np.block([[np.zeros((2, 2)), I2], [I2, np.zeros((2, 2))]]).astype(complex)   # tau = Theta_1
P4 = np.array([1, 1j, -1, -1j])                                         # i^k exakt
LC = np.zeros((4,)*4)
for _p in itertools.permutations(range(4)):
    LC[_p] = round(np.linalg.det(np.eye(4)[list(_p)]))
uvec = [np.array([1, 0, -1, 0.]), np.array([0, 1, 0, -1.])]
I_PI = (np.outer(uvec[0], uvec[1]) - np.outer(uvec[1], uvec[0])).astype(complex)  # I_pi (z89, z98)


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


def laufzeit():
    say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")


def nullraum(M, tol=1e-9):
    M = np.atleast_2d(M)
    _, s, vh = np.linalg.svd(M)
    rank = int(np.sum(s > tol*max(1.0, s[0]))) if s.size else 0
    return vh[rank:]


# ------------------------------------------------------------------ Fockräume (wie z98)
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


def tausch(nv, naehte):
    """Teilchen-Loch-Tausch d <-> e in jeder Naht (reelle Permutation der Fockbasis)."""
    k = 4*naehte
    occ = np.array(list(itertools.product(range(nv), repeat=k)))
    perm = occ.copy()
    for s in range(naehte):
        b = 4*s
        perm[:, b:b + 2], perm[:, b + 2:b + 4] = occ[:, b + 2:b + 4], occ[:, b:b + 2]
    idx = np.ravel_multi_index(perm.T, (nv,)*k)
    dim = nv**k
    return sps.csr_matrix((np.ones(dim), (idx, np.arange(dim))), shape=(dim, dim)), occ


def diag_op(werte):
    return sps.diags(np.asarray(werte, dtype=complex)).tocsr()


def norm_max(A):
    A = sps.csr_matrix(A)
    return float(abs(A).max()) if A.nnz else 0.0


# eine Naht, Fock 5^4: Twistor-Matrizen der Daten
NV1 = 5
d1, d2, e1, e2 = fock(NV1, 4)
SW1, OCC1 = tausch(NV1, 1)
DV1 = OCC1[:, 0] + OCC1[:, 1] - OCC1[:, 2] - OCC1[:, 3]
LOW1 = np.where(OCC1.sum(1) <= NV1 - 3)[0]
ZT = [d1, d2, dag(e1), dag(e2)]                                         # Z^alpha
ZTd = [dag(z) for z in ZT]


def U1(th):
    return diag_op(np.exp(1j*th*DV1))


# zwei Nähte, Fock 3^8: der eps-Kern V = eps^jk c_j^(1)+ c_k^(2)+ (Modulladung 2)
NV2 = 3
q8 = fock(NV2, 8)
SW2, OCC2 = tausch(NV2, 2)
DV2 = (OCC2[:, 0] + OCC2[:, 1] - OCC2[:, 2] - OCC2[:, 3]) + (OCC2[:, 4] + OCC2[:, 5] - OCC2[:, 6] - OCC2[:, 7])
DSEAM = [OCC2[:, 0] + OCC2[:, 1] - OCC2[:, 2] - OCC2[:, 3], OCC2[:, 4] + OCC2[:, 5] - OCC2[:, 6] - OCC2[:, 7]]
CD = [[dag((q8[4*s + j] + dag(q8[4*s + 2 + j])).tocsr()) for j in range(2)] for s in range(2)]   # c^+ = d^+ + e
V0 = (CD[0][0] @ CD[1][1] - CD[0][1] @ CD[1][0]).tocsr()
V0d = dag(V0)


def U2(th):
    return diag_op(np.exp(1j*th*DV2))


def ip(X, Y):
    return complex(sps.csr_matrix(X).conj().multiply(sps.csr_matrix(Y)).sum())


NV0 = ip(V0, V0).real


def B_kern(phi):
    return (np.exp(1j*phi)*V0 + np.exp(-1j*phi)*V0d).tocsr()


def zerlege(Y):
    """Y = a V + b V^+ + Rest (V und V^+ liegen in verschiedenen Ladungsblöcken, also orthogonal)."""
    a, b = ip(V0, Y)/NV0, ip(V0d, Y)/NV0
    return a, b, norm_max(Y - a*V0 - b*V0d)


# ------------------------------------------------------------------ Plebanski-Sektoren in der Raumzeit (wie z98, Teil 2.4)
EPS_UNTEN = np.einsum("abcd,ai,bj,ck,dl->ijkl", LC, ETA, ETA, ETA, ETA)


def stern(F):
    return 0.5*np.einsum("mnrs,rs->mn", EPS_UNTEN, ETA @ F @ ETA)


def keil(F, G_):
    return 0.25*np.einsum("mnrs,mn,rs->", LC, F, G_)


SIGP00 = np.array([[(SIGMA[m] @ EPS2 @ SIGMA[n].T - SIGMA[n] @ EPS2 @ SIGMA[m].T)[0, 0] for n in range(4)] for m in range(4)])
LAM_D = complex(np.vdot(SIGP00.ravel(), stern(SIGP00).ravel())/np.vdot(SIGP00.ravel(), SIGP00.ravel()))
EV = np.zeros((4, 4, 4, 4))
for _a, _b in itertools.product(range(4), repeat=2):
    EV[_a, _b, _a, _b] += 1
    EV[_a, _b, _b, _a] -= 1
SPV = np.array([[0.5*(EV[a, b] + stern(EV[a, b])/LAM_D) for b in range(4)] for a in range(4)])
STE = np.array([[stern(EV[a, b]) for b in range(4)] for a in range(4)])


def B_raum(phi):
    return np.exp(1j*phi)*SPV + np.exp(-1j*phi)*SPV.conj()


def simpl_rest(phi):
    B = B_raum(phi)
    K = np.array([[[[keil(B[a, b], B[c, d]) for d in range(4)] for c in range(4)] for b in range(4)] for a in range(4)])
    cK = np.sum(LC*K)/np.sum(LC*LC)
    return float(np.linalg.norm(K - cK*LC))


SEKTOR_NAME = {}
for _k in range(4):
    _B = B_raum(_k*np.pi/2)
    for _name, _ref in (("e^e", EV), ("-e^e", -EV), ("*(e^e)", STE), ("-*(e^e)", -STE)):
        if np.allclose(_B, _ref, atol=1e-12):
            SEKTOR_NAME[_k] = _name
GRAV, TOPO = {0, 2}, {1, 3}


def sektor_von_phase(phi, tol=1e-9):
    k = int(np.round(phi/(np.pi/2))) % 4
    return k if abs(np.angle(np.exp(1j*(phi - k*np.pi/2)))) < tol else None


# ------------------------------------------------------------------ lorentzsche Realität des Twistorbivektors (wie z98, 1.2)
def stern_unten(I, psi):
    return 0.5*np.exp(1j*psi)*np.einsum("abcd,cd->ab", LC, I)


def quer_unten(I):
    return G4 @ I.conj() @ G4


def typ(I, psi=0.0):
    """+1: I reell (I-quer = *I, Gravitationspaar), -1: imaginär (I-quer = -*I, topologisches Paar), 0: keins."""
    q, s = quer_unten(I), stern_unten(I, psi)
    sk = max(float(np.abs(q).max()), 1e-300)
    if np.abs(q - s).max() < 1e-10*sk:
        return 1
    if np.abs(q + s).max() < 1e-10*sk:
        return -1
    return 0


def I_k(k):
    return np.exp(-1j*k*np.pi/2)*I_PI          # Twistorbivektor zum Kern e^(i k pi/2) V (I trägt die Ladung -2, V die +2)


# =================================================================== Teil 1
say("Teil 1: Zensus der Orientierungsdaten")

# --- 1.1 Twistor-Matrizen der Daten (eine Naht, Fock 5^4): F(Z^alpha) = M Z + N Z^+
BASIS8 = ZT + ZTd
A8 = np.array([sps.csc_matrix(B)[:, LOW1].toarray().ravel() for B in BASIS8]).T


def twistor_matrix(F):
    MN = np.zeros((4, 8), complex)
    res = 0.0
    for al in range(4):
        y = sps.csc_matrix(F(ZT[al]))[:, LOW1].toarray().ravel()
        x = np.linalg.lstsq(A8, y, rcond=None)[0]
        res = max(res, float(np.abs(A8 @ x - y).max()))
        MN[al] = x
    return MN[:, :4], MN[:, 4:], res


ZSU1 = diag_op(P4[(DV1 - 2) % 4])                                       # Z_SU = i^(2h - 2), 2h = D (z97)


def Jad1(X):
    return (SW1 @ sps.csr_matrix(X).conj() @ SW1).tocsr()


DATEN1 = {
    "U_(pi/4) (Modul)": (lambda X: (U1(np.pi/4) @ X @ dag(U1(np.pi/4))).tocsr(), False),
    "Z_SU = i^(2h-2)": (lambda X: (ZSU1 @ X @ dag(ZSU1)).tocsr(), False),
    "Z_SU^-1": (lambda X: (dag(ZSU1) @ X @ ZSU1).tocsr(), False),
    "tau = (Ad J)^+": (lambda X: dag(Jad1(X)), False),
    "Tomita J": (Jad1, True),
    "K (Familie D, modenerhaltend)": (lambda X: sps.csr_matrix(X).conj().tocsr(), True),
    "C (Teilchen-Loch, unitär)": (lambda X: (SW1 @ X @ SW1).tocsr(), False),
}
TMAT = {}
ok_fit = True
for name, (F, anti) in DATEN1.items():
    M, N, res = twistor_matrix(F)
    TMAT[name] = (M, N, anti)
    ok_fit &= res < 1e-10
# auf T wirkende Daten: Koeffizientenmatrix A = M^T (Z = sum z_alpha Z^alpha -> sum phi(z_alpha) F(Z^alpha))
A_U = TMAT["U_(pi/4) (Modul)"][0].T
A_Z = TMAT["Z_SU = i^(2h-2)"][0].T
A_tau = TMAT["tau = (Ad J)^+"][0].T
A_K = TMAT["K (Familie D, modenerhaltend)"][0].T
ok_T = (np.allclose(A_U, np.exp(-1j*np.pi/4)*np.eye(4)) and np.allclose(A_Z, -1j*np.eye(4)) and np.allclose(A_tau, S4)
        and np.allclose(A_K, np.eye(4))
        and all(np.allclose(TMAT[k][1], 0) for k in ("U_(pi/4) (Modul)", "Z_SU = i^(2h-2)", "tau = (Ad J)^+",
                                                     "K (Familie D, modenerhaltend)"))
        and all(np.allclose(TMAT[k][0], 0) and np.allclose(TMAT[k][1], S4) for k in ("Tomita J", "C (Teilchen-Loch, unitär)")))


def wirkung_T(A, anti):
    """Faktor auf eps (Rückzug eps o A^-1) und neues I_pi (A phi(I) A^T)."""
    I_neu = A @ (I_PI.conj() if anti else I_PI) @ A.T
    lam = complex(np.vdot(I_PI.ravel(), I_neu.ravel())/np.vdot(I_PI.ravel(), I_PI.ravel()))
    auf_linie = np.allclose(I_neu, lam*I_PI)
    return complex(1/np.linalg.det(A)), (lam if auf_linie else None)


T_TAB = {}
for name, A, anti in (("U_(pi/4) (Modul)", A_U, False), ("Z_SU = i^(2h-2)", A_Z, False), ("tau = (Ad J)^+", A_tau, False),
                      ("K (Familie D, modenerhaltend)", A_K, True)):
    f_eps, lam = wirkung_T(A, anti)
    T_TAB[name] = {"linear": not anti, "eps-Faktor (Rückzug)": str(np.round(f_eps, 12)) + (" und konjugiert" if anti else ""),
                   "I_pi -> lambda I_pi (phi -> " + ("-phi" if anti else "phi") + " + arg lambda)":
                       str(np.round(lam, 12)) if lam is not None else "verlässt die Linie"}
ZAHLEN["twistor_matrizen"] = T_TAB
check(f"Twistor-Matrizen der Daten (Fock 5^4, Fit auf T + T^+, Residuum < 1e-10): U_(pi/4) wirkt auf T als e^(-i pi/4) 1 "
      f"(det = -1, also eps -> -eps, I_pi -> -i I_pi), Z_SU = i^(D-2) als -i 1 (det 1, eps fest, I_pi -> -I_pi), tau als "
      f"S4 (det 1, I_pi fest), die modenerhaltende Konjugation K als Konjugation (eps_psi -> eps_-psi, I_pi reell fest, "
      f"Phase phi -> -phi); Tomita J und der unitäre Teilchen-Loch-Tausch C bilden T auf T^+ ab (Z^alpha -> (S4 Z)^+), sie "
      f"wirken auf T-Strukturen über tau bzw. tau K", ok_fit and ok_T)

# --- 1.2 Sektorwirkung auf dem Kern zweier Nähte (Fock 3^8) und Verhalten unter U_(pi/4)
ZSU2 = diag_op(P4[(DSEAM[0] - 2 + DSEAM[1] - 2) % 4])


def J2(X):
    return (SW2 @ sps.csr_matrix(X).conj() @ SW2).tocsr()


def K2(X):
    return sps.csr_matrix(X).conj().tocsr()


def C2(X):
    return (SW2 @ X @ SW2).tocsr()


def AdU(th):
    U = U2(th)
    return lambda X: (U @ X @ dag(U)).tocsr()


DATEN2 = {
    "U_(pi/4) (Modul)": AdU(np.pi/4),
    "Z_SU = i^(2h-2)": lambda X: (ZSU2 @ X @ dag(ZSU2)).tocsr(),
    "Z_SU^-1": lambda X: (dag(ZSU2) @ X @ ZSU2).tocsr(),
    "tau = (Ad J)^+": lambda X: dag(J2(X)),
    "Tomita J": J2,
    "Uhr R (Ad_R = 1 auf den Oszillatoren)": lambda X: sps.csr_matrix(X),
    "K (Familie D, modenerhaltend)": K2,
    "C (Teilchen-Loch, unitär)": C2,
}


def sektor_wirkung(F):
    """Permutation der vier Sektoren und affines Phasengesetz phi -> s phi + delta (aus phi = 0 und phi = 0,37)."""
    perm, ok = [], True
    for k in range(4):
        a, b, rest = zerlege(F(B_kern(k*np.pi/2)))
        ok &= rest < 1e-12 and abs(abs(a) - 1) < 1e-12 and abs(b - np.conj(a)) < 1e-12
        kk = sektor_von_phase(np.angle(a))
        ok &= kk is not None and simpl_rest(np.angle(a)) < 1e-10
        perm.append(kk)
    a0 = zerlege(F(B_kern(0.0)))[0]
    a1 = zerlege(F(B_kern(0.37)))[0]
    s = int(np.round(np.angle(a1/a0)/0.37))
    delta = float(np.mod(np.angle(a0), 2*np.pi))
    return tuple(perm), s, delta, ok


def vertauscht_mit_U(F, anti=False, th=np.pi/4):
    """F o Ad U = Ad U o F (vertauscht) oder F o Ad U = Ad U^-1 o F (invertiert), geprüft auf V, V^+ und B(0,37)."""
    AU, AUi = AdU(th), AdU(-th)
    tests = [V0, V0d, B_kern(0.37), (V0 @ V0d).tocsr()]
    r_v = max(norm_max(F(AU(X)) - AU(F(X))) for X in tests)
    r_i = max(norm_max(F(AU(X)) - AUi(F(X))) for X in tests)
    return ("vertauscht" if r_v < 1e-12 else ("invertiert" if r_i < 1e-12 else "bricht")), r_v, r_i


S_TAB = {}
ok_sek = True
for name, F in DATEN2.items():
    perm, s, delta, ok = sektor_wirkung(F)
    verh, r_v, r_i = vertauscht_mit_U(F)
    S_TAB[name] = {"Sektorpermutation (0: e^e, 1, 2, 3)": perm, "phi -> s phi + delta": (s, round(delta, 12)),
                   "U_(pi/4)": verh}
    ok_sek &= ok
ROT = (1, 2, 3, 0)
ok_sek &= (S_TAB["U_(pi/4) (Modul)"]["Sektorpermutation (0: e^e, 1, 2, 3)"] == ROT
           and S_TAB["Z_SU = i^(2h-2)"]["Sektorpermutation (0: e^e, 1, 2, 3)"] == (2, 3, 0, 1)
           and S_TAB["Z_SU^-1"]["Sektorpermutation (0: e^e, 1, 2, 3)"] == (2, 3, 0, 1)
           and all(S_TAB[k]["Sektorpermutation (0: e^e, 1, 2, 3)"] == (0, 1, 2, 3)
                   for k in ("tau = (Ad J)^+", "Tomita J", "Uhr R (Ad_R = 1 auf den Oszillatoren)"))
           and all(S_TAB[k]["Sektorpermutation (0: e^e, 1, 2, 3)"] == (0, 3, 2, 1)
                   for k in ("K (Familie D, modenerhaltend)", "C (Teilchen-Loch, unitär)"))
           and all(S_TAB[k]["U_(pi/4)"] == "vertauscht" for k in ("Z_SU = i^(2h-2)", "Z_SU^-1", "tau = (Ad J)^+", "Tomita J",
                                                                  "Uhr R (Ad_R = 1 auf den Oszillatoren)", "U_(pi/4) (Modul)"))
           and all(S_TAB[k]["U_(pi/4)"] == "invertiert" for k in ("K (Familie D, modenerhaltend)", "C (Teilchen-Loch, unitär)")))
ok_sek &= SEKTOR_NAME[0] == "e^e" and SEKTOR_NAME[2] == "-e^e" and set(SEKTOR_NAME.values()) == {"e^e", "-e^e", "*(e^e)",
                                                                                                   "-*(e^e)"}
ZAHLEN["sektoren"] = {"Namen": SEKTOR_NAME, "Daten": S_TAB}
check(f"Sektorwirkung auf dem Kern zweier Nähte (Fock 3^8, B = e^(i phi) V + e^(-i phi) V^+, Bild exakt wieder von dieser "
      f"Form, Bildphase einfach nach z98): Sektoren {SEKTOR_NAME}. U_(pi/4) dreht zyklisch {ROT} (phi -> phi + pi/2); "
      f"Z_SU und Z_SU^-1 wirken beide als Halbdrehung (2, 3, 0, 1) (B -> -B); tau, Tomita J und die Uhr R wirken "
      f"trivial; nur K (Familie D) und der unitäre Teilchen-Loch-Tausch C wirken als Spiegelung (0, 3, 2, 1) (phi -> -phi, "
      f"fest: e^e und -e^e). Verhalten unter U_(pi/4) (auf V, V^+, B, V V^+): Z_SU, tau, J und R vertauschen, K und C "
      f"invertieren den Modulfluss (K U K = U^-1, C U C = U^-1)", ok_sek)

# --- 1.3 die Daten einzeln: Strukturen, Umkehrabbildungen, Modulinvarianz
OMV1 = np.zeros(NV1**4, complex)
OMV1[0] = 1                                                             # KMS-Vektor (Araki-Woods-Vakuum, z98)
U1q = U1(np.pi/4)
ok_kms = np.allclose(U1q @ OMV1, OMV1)
Dop1 = diag_op(DV1.astype(float))
JDJ = Jad1(Dop1)
ok_kms &= norm_max(JDJ + Dop1) < 1e-14 and norm_max(Jad1(U1q) - U1q) < 1e-13           # J D J = -D, J U J = U
ok_kms &= norm_max((SW1 @ Dop1 @ SW1) + Dop1) < 1e-14 and norm_max(SW1 @ U1q @ SW1 - dag(U1q)) < 1e-13   # C D C = -D
ok_G = np.allclose(A_U.conj().T @ G4 @ A_U, G4) and np.allclose(S4.conj().T @ G4 @ S4, -G4)
P_plus = np.diag([1, 1, 0, 0]).astype(complex)                          # Hardy: Teilchen oben (G-positiv)
ok_hardy = np.allclose(A_U @ P_plus @ np.linalg.inv(A_U), P_plus) and np.allclose(S4 @ P_plus @ S4, np.eye(4) - P_plus)


def Qn(n):
    return n[0]*I2 + 1j*(n[1]*PAULI[0] + n[2]*PAULI[1] + n[3]*PAULI[2])


def Theta(n):
    return np.block([[np.zeros((2, 2)), Qn(n)], [Qn(n).conj().T, np.zeros((2, 2))]])


ok_theta = np.allclose(Theta([1, 0, 0, 0]), S4)
verlaesst = 0
NS = [v/np.linalg.norm(v) for v in rng.normal(size=(8, 4))]
for n in NS:
    Th = Theta(n)
    ok_theta &= (np.allclose(Th @ Th, np.eye(4)) and np.allclose(Th.conj().T @ G4 @ Th, -G4)
                 and abs(np.linalg.det(Th) - 1) < 1e-12 and np.allclose(Th @ A_U, A_U @ Th))
    for k in range(4):
        I_neu = Th @ I_k(k) @ Th.T
        ok_theta &= typ(I_neu) == typ(I_k(k))
        lam = np.vdot(I_PI.ravel(), I_neu.ravel())/np.vdot(I_PI.ravel(), I_PI.ravel())
        verlaesst += int(not np.allclose(I_neu, lam*I_PI))
ok_typ = [typ(I_k(k)) for k in range(4)] == [1, -1, 1, -1]
ZAHLEN["strukturen"] = {"KMS-Vektor U-invariant": bool(np.allclose(U1q @ OMV1, OMV1)), "J D J": "-D", "J U J": "U",
                        "C D C": "-D", "C U C": "U^-1", "U^+ G U": "G", "tau^+ G tau": "-G",
                        "Hardy-Projektor P_+ U-invariant, tau: P_+ -> P_-": bool(ok_hardy),
                        "Theta_n: Paartyp erhalten (8 n x 4 Sektoren)": bool(ok_theta),
                        "Theta_n: I_pi verlässt die Linie (von 32)": int(verlaesst),
                        "Paartyp der Sektoren (typ I_k, psi = 0)": [typ(I_k(k)) for k in range(4)]}
check(f"Orientierungsstrukturen und ihre Umkehr: (KMS) der KMS-Vektor ist U_(pi/4)-fest, die Richtung beta > 0 ist das "
      f"Vorzeichen des Modulgenerators D; Tomita J kehrt sie um (J D J = -D) und vertauscht doch mit U_theta (J U J = U), "
      f"der unitäre Tausch C kehrt sie um und invertiert U (C U C = U^-1). (Signatur und Hardy) U^+ G U = G und der "
      f"Teilchenprojektor P_+ ist U-fest; die Umkehr G -> -G, P_+ -> P_- leistet tau (linear, det 1, Sektoren fest). "
      f"(Theta_n) alle acht zufälligen Theta_n sind Involutionen, G-antiisometrisch, det 1, vertauschen mit U und erhalten "
      f"den Paartyp (reell/imaginär bezüglich eps_0) jedes Sektors, obwohl I_pi in {verlaesst}/32 Fällen die Linie verlässt; "
      f"Paartyp der Sektoren {[typ(I_k(k)) for k in range(4)]} (+1 = e^e-Paar)",
      ok_kms and ok_G and ok_hardy and ok_theta and ok_typ and verlaesst == 32)
laufzeit()

# --- 1.4 Uhr und Marken: Quarterfermionen (R = i^(Q_a - Q_b), trivial auf den Strömen)
MODI = [Fr(-7, 4), Fr(-3, 4), Fr(1, 4), Fr(5, 4)] + [Fr(-5, 4), Fr(-1, 4), Fr(3, 4), Fr(7, 4)]   # Quarter a, Quarter b
NM = len(MODI)
sz, sm = np.diag([1., -1.]), np.array([[0., 1.], [0., 0.]])


def jw(j):
    op = np.eye(1)
    for m in range(NM):
        op = np.kron(op, sz if m < j else (sm if m == j else np.eye(2)))
    return op


CF = [jw(j) for j in range(NM)]
NOCC = [np.real(np.diag(c.T @ c)) for c in CF]
NORM = [NOCC[j] - (1 if MODI[j] < 0 else 0) for j in range(NM)]            # :n_x: (Diracsee gefüllt)
L0 = sum(float(MODI[j])*NORM[j] for j in range(NM))
QA = sum(NORM[j] for j in range(4))
QB = sum(NORM[j] for j in range(4, 8))
R_f = np.exp(2j*np.pi*L0)
ok_R = np.allclose(R_f, P4[np.round(QA - QB).astype(int) % 4])
Rm = np.diag(R_f)
ok_strom = True
n_strom = 0
for blk in (range(4), range(4, 8)):
    for n in (-3, -2, -1, 1, 2, 3):
        Jn = sum(CF[i].T @ CF[j] for i in blk for j in blk if MODI[i] - MODI[j] == n)
        if isinstance(Jn, np.ndarray):
            ok_strom &= np.allclose(Rm @ Jn @ Rm.conj().T, Jn)
            n_strom += 1
ph_a = [complex((Rm @ CF[j] @ Rm.conj().T)[np.nonzero(CF[j])][0]/CF[j][np.nonzero(CF[j])][0]) for j in (2, 6)]
ok_R &= ok_strom and np.allclose(ph_a, [-1j, 1j]) and np.allclose(np.diag(Rm.conj()), 1/R_f)
ZAHLEN["uhr_fermionen"] = {"Moden": [str(x) for x in MODI], "Fockdimension": 2**NM, "Ströme geprüft": n_strom,
                           "R psi_x R^-1 / psi_x (x = 1/4, 3/4)": [str(np.round(p, 12)) for p in ph_a]}
check(f"Uhr auf Quarterfermionen (acht Moden x in Z + 1/4 bzw. Z + 3/4, Fock 2^8): R = e^(2 pi i L_0) = i^(Q_a - Q_b) exakt "
      f"auf allen 256 Zuständen; alle {n_strom} Strommoden J_n = sum psi*_(x+n) psi_x (n = +-1, +-2, +-3, beide Quarter) "
      f"vertauschen mit R, geladene Moden nicht (R psi R^-1 = {[str(np.round(p, 6)) for p in ph_a]} psi); die Konjugation "
      f"der Besetzungsbasis gibt R^-1 (Familie D). Die Twistoroszillatoren sind Stromoszillatoren (z79), also wirkt R auf T, "
      f"I_pi, eps und die Sektoren als Eins, und R vertauscht mit U_(pi/4)", ok_R)

# --- 1.5 Uhr und Marken auf der Spinbrücke H0(N) = V^* (x) eta (Kurventransport, OriginLab)
R_S = _transport(1, -1)
DECK_S = _transport(1, 0)
R_np = np.array(R_S.evalf(30).tolist(), dtype=complex)
DECK_np = np.array(DECK_S.evalf(30).tolist(), dtype=complex)


def faktoren(M):
    """M[i + 2j, k + 2l] = a[i, k] b[j, l] (i: V^*, j: eta, wie z84), also M = kron(b, a)."""
    Rm_ = np.zeros((4, 4), complex)
    for i, j, k, l in itertools.product(range(2), repeat=4):
        Rm_[2*i + k, 2*j + l] = M[i + 2*j, k + 2*l]
    u_, s_, vh_ = np.linalg.svd(Rm_)
    return (u_[:, 0]*np.sqrt(s_[0])).reshape(2, 2), (vh_[0]*np.sqrt(s_[0])).reshape(2, 2), int(np.sum(s_ > 1e-9*s_[0]))


gV, gE, rang_R = faktoren(R_np)
lamV, lamE = np.sqrt(complex(np.linalg.det(gV))), np.sqrt(complex(np.linalg.det(gE)))
lam_c = lamV*lamE                                                       # R = lam (gV/lamV) (x) (gE/lamE), Faktoren in SL(2)
Q_prod = np.kron(EPS2, EPS2)                                            # eps_eta (x) eps_V, Index i + 2j
c_konf = complex(np.vdot(Q_prod.ravel(), (R_np.T @ Q_prod @ R_np).ravel())/np.vdot(Q_prod.ravel(), Q_prod.ravel()))
ok_brücke = (R_S.det() == -1 and R_S**4 == sp.eye(4) and set(R_S.eigenvals().keys()) == {1, sp.I, -1, -sp.I}
             and DECK_S.det() == 1 and rang_R == 1 and np.allclose(np.kron(gE, gV), R_np)
             and abs(lam_c**4 + 1) < 1e-12 and abs(lam_c**8 - 1) < 1e-12
             and np.allclose(R_np.T @ Q_prod @ R_np, c_konf*Q_prod) and min(abs(c_konf - 1j), abs(c_konf + 1j)) < 1e-12)
# R-invariante symmetrische Bilinearformen: R^T q R = q
bil = []
for M in [np.outer(np.eye(4)[a], np.eye(4)[b]) + np.outer(np.eye(4)[b], np.eye(4)[a]) for a, b in
          itertools.combinations_with_replacement(range(4), 2)]:
    bil.append(M)
inv_q = nullraum(np.array([(R_np.T @ M @ R_np - M).ravel() for M in bil]).T)
q_inv = [sum(c*M for c, M in zip(v, bil)) for v in inv_q]
q_test = sum(rng.normal()*Q for Q in q_inv)
A_q = np.array([Q.ravel() for Q in q_inv]).T
rest_prod = float(np.linalg.norm(A_q @ np.linalg.lstsq(A_q, Q_prod.ravel(), rcond=None)[0] - Q_prod.ravel()))
R_deck = R_np @ DECK_np
ok_brücke &= ((R_S**3).det() == -1 and len(q_inv) == 3 and abs(np.linalg.det(q_test)) > 1e-6 and rest_prod > 0.5
              and abs(np.linalg.det(R_deck) + 1) < 1e-12 and np.allclose(np.linalg.matrix_power(R_deck, 4), np.eye(4)))
ZAHLEN["spinbruecke"] = {"R auf H0(N) (Basis 1, z, F, zF)": str(R_S.tolist()), "det R": str(R_S.det()),
                         "Eigenwerte": [str(v) for v in R_S.eigenvals()], "det Deck": str(DECK_S.det()),
                         "Zentralphase lambda (R = lambda gV^ (x) gE^, SL(2)-Faktoren)": str(np.round(lam_c, 12)),
                         "lambda^4": str(np.round(lam_c**4, 12)),
                         "konformer Faktor auf eps_V (x) eps_eta": str(np.round(c_konf, 12)),
                         "Dimension R-invarianter symmetrischer Formen": len(q_inv),
                         "Abstand der Produktmetrik von ihrer Spanne": rest_prod,
                         "det (R Deck) (andere Hebung)": str(np.round(np.linalg.det(R_deck), 12))}
check(f"Uhr auf der Spinbrücke (OriginLab-Transport, exakt): R hat auf H0(N) = V^* (x) eta die Eigenwerte {{1, i, -1, -i}}, "
      f"R^4 = 1 und det R = {R_S.det()} (das Deck hat det {DECK_S.det()}, die andere Hebung R Deck ebenfalls det -1); R ist ein Produkt (Realignment-Rang {rang_R}) "
      f"= lambda (gV^ (x) gE^) mit SL(2)-Faktoren und der Zentralphase lambda = {np.round(lam_c, 6)} (lambda^4 = -1, "
      f"primitive achte Einheitswurzel); die Produktmetrik eps_V (x) eps_eta geht in {np.round(c_konf, 6)} mal sich selbst "
      f"über (komplex konform). R-invariant sind nur {len(q_inv)} symmetrische Formen (Paarung der Eigenräume 1-1, -1-1, "
      f"i-(-i)), keine davon die Produktmetrik; für sie ist R uneigentlich (det -1). Wegen R^4 = 1 ist det R hebungsfrei, "
      f"und die Gegenrichtung R^-1 = R^3 hat ebenfalls det -1: keine Volumenform auf H0(N) ist R-invariant, und die "
      f"zyklische Ordnung der Marken ist an ihr nicht ablesbar. Die Faktoren V^* und eta (selbstdual gegen antiselbstdual) "
      f"erhält R (z84, Satz A); der Vorzeichenwechsel kommt aus der komplexen Konformphase, nicht aus einem Faktortausch",
      ok_brücke)

# --- 1.6 eta-Dublett: Hebungen, Deck, Volumenform eps_eta; keine äquivariante Übertragung auf die pi-Ebene
u = sp.symbols("u")


def eta_matrix(vorz):
    """Transport auf span{1, f}, f = w/(u^2 - 1), Kurve w^2 = 1 - u^4 (wie z82): (T b)(u, w) = f b(-i u, vorz w)."""
    Fq = 1/(u**2 - 1)
    basis = [(sp.Integer(1), sp.Integer(0)), (sp.Integer(0), Fq)]
    M = sp.zeros(2, 2)
    for j, (A, Bw) in enumerate(basis):
        A2, B2 = A.subs(u, -sp.I*u), vorz*Bw.subs(u, -sp.I*u)
        prod = (sp.cancel(Fq*B2*(1 - u**4)), sp.cancel(Fq*A2))
        a, c = sp.cancel(prod[0]), sp.cancel(prod[1]*(u**2 - 1))
        M[0, j], M[1, j] = a, c
    return M


G_SPIN, G_ANDERE = eta_matrix(-1), eta_matrix(1)
DECK_ETA = sp.diag(1, -1)                                               # w -> -w: f -> -f
ok_eta = (G_SPIN == sp.Matrix([[0, -1], [1, 0]]) and G_SPIN.det() == 1 and G_ANDERE.det() == -1 and DECK_ETA.det() == -1
          and np.allclose(np.kron(np.array(DECK_ETA.tolist(), dtype=float), np.eye(2)), DECK_np.real))
Rs = np.array(G_SPIN.tolist(), dtype=complex)
w_, Vw = np.linalg.eig(Rs)
Jz = Vw @ np.diag([0.5 if abs(x - 1j) < 1e-9 else -0.5 for x in w_]) @ np.linalg.inv(Vw)   # R_eta = e^(i pi J_z)
ok_eta &= np.allclose(Vw @ np.diag(np.exp(1j*np.pi*np.diag(np.linalg.inv(Vw) @ Jz @ Vw))) @ np.linalg.inv(Vw), Rs)


def intertwiner_dim(A_quelle, A_ziel):
    """Komplexe Dimension von {X: X A_quelle = A_ziel X} (lineare Abbildungen eta -> pi-Ebene)."""
    cols = []
    for k in range(4):
        X = np.zeros(4, complex)
        X[k] = 1
        X = X.reshape(2, 2)
        cols.append((X @ A_quelle - A_ziel @ X).ravel())
    return 4 - int(np.linalg.matrix_rank(np.array(cols).T, tol=1e-9))


def expJz(th):
    return Vw @ np.diag(np.exp(1j*th*np.diag(np.linalg.inv(Vw) @ Jz @ Vw))) @ np.linalg.inv(Vw)


dim_R = intertwiner_dim(Rs, I2)                                         # R auf der pi-Ebene = 1 (Ad_R = 1)
dim_mod = {th: intertwiner_dim(expJz(th), np.exp(-1j*th)*I2) for th in (0.7, np.pi/4)}
ok_eta &= dim_R == 0 and all(v == 0 for v in dim_mod.values()) and abs(np.linalg.det(expJz(0.7)) - 1) < 1e-12
ZAHLEN["eta_dublett"] = {"Spinorhebung": str(G_SPIN.tolist()), "det": str(G_SPIN.det()), "andere Hebung": str(G_ANDERE.tolist()),
                         "det andere": str(G_ANDERE.det()), "Deck auf eta": "diag(1, -1), det -1",
                         "R-äquivariante eta -> pi": dim_R, "modular äquivariante eta -> pi (theta = 0,7; pi/4)":
                             list(dim_mod.values())}
check(f"eta-Dublett (exakt): die Spinorhebung ist {G_SPIN.tolist()} mit det 1, die andere Hebung hat det -1, das Deck "
      f"(w -> -w) wirkt als diag(1, -1) mit det -1: die Volumenform eps_eta ist nativ nur bis aufs Vorzeichen bestimmt "
      f"(Blattwahl). Der Modultorus wirkt auf eta in SL(2) (e^(i theta J_z), det 1, R_eta = e^(i pi J_z)), auf der "
      f"pi-Ebene aber als e^(-i theta) 1 (det e^(-2 i theta)), und R wirkt dort als Eins: Es gibt keine R-äquivariante "
      f"(Dimension {dim_R}) und keine modular äquivariante (Dimension {list(dim_mod.values())}) lineare Abbildung "
      f"eta -> pi-Ebene, die eps_eta auf I_pi übertragen könnte", ok_eta)
eps_J, lam_J = wirkung_T(S4, False)                                     # J wirkt auf T-Strukturen über tau
eps_C, lam_C = wirkung_T(S4, True)                                      # C über tau K


def perm_von(k):
    return S_TAB[k]["Sektorpermutation (0: e^e, 1, 2, 3)"]


def kurz_c(z):
    return str(np.round(complex(z), 6)).replace("(", "").replace(")", "")


ZENSUS = {
    "1 Marken, Uhr R = i^(Q_a - Q_b)": ("eps fest (R auf T = 1)", "I_pi fest", perm_von("Uhr R (Ad_R = 1 auf den Oszillatoren)"),
                                        S_TAB["Uhr R (Ad_R = 1 auf den Oszillatoren)"]["U_(pi/4)"],
                                        f"auf H0(N) det R = det R^-1 = -1, lambda^2 = {kurz_c(c_konf)}"),
    "2 Hardy-Orientierung (P_+, Teilchen oben)": ("eps fest", "I_pi fest", "Umkehr tau: " + str(perm_von("tau = (Ad J)^+"))
                                                  + ", Umkehr C: " + str(perm_von("C (Teilchen-Loch, unitär)")),
                                                  "P_+ fest; C invertiert", "C rahmenabhängig (2.3)"),
    "3 tau = Theta_1, Theta_n": ("det 1", "tau: I_pi fest; Theta_n: Paartyp fest", perm_von("tau = (Ad J)^+"),
                                 S_TAB["tau = (Ad J)^+"]["U_(pi/4)"], f"I_pi verlässt die Linie in {verlaesst}/32"),
    "4 Uhrrichtung R gegen R^-1": (f"J: über tau, Faktor {kurz_c(eps_J)}; K: Konjugation", f"J: lambda = {kurz_c(lam_J)}; "
                                   f"K: phi -> -phi", "J: " + str(perm_von("Tomita J")) + ", K: "
                                   + str(perm_von("K (Familie D, modenerhaltend)")),
                                   "J " + S_TAB["Tomita J"]["U_(pi/4)"] + ", K "
                                   + S_TAB["K (Familie D, modenerhaltend)"]["U_(pi/4)"], "K rahmenabhängig (2.3)"),
    "5 KMS-Richtung beta > 0 (Vorzeichen von D)": ("eps fest", "I_pi fest", "Umkehr J: " + str(perm_von("Tomita J"))
                                                   + ", Umkehr C: " + str(perm_von("C (Teilchen-Loch, unitär)")),
                                                   "D fest; J vertauscht, C invertiert", f"C über tau K: eps-Faktor "
                                                   f"{kurz_c(eps_C)} und konjugiert, lambda = {kurz_c(lam_C)}"),
    "6 Signatur (2,2) der CCR-Form": ("eps fest", "I_pi fest", "Umkehr tau: " + str(perm_von("tau = (Ad J)^+")),
                                      "G fest", "tau^+ G tau = -G"),
    "7 eta-Dublett": ("anderer Faktor; eps_eta nur bis aufs Vorzeichen (Deck det -1)", "keine äquivariante "
                      "Übertragung (Dimension 0)", "keine Wirkung auf T", "Modultorus erhält die Hebung", "Spinorhebung det 1"),
    "8 Zentrum Z_SU = i^(2h-2)": ("det 1", "I_pi -> -I_pi", str(perm_von("Z_SU = i^(2h-2)")) + " (auch Z_SU^-1)",
                                  S_TAB["Z_SU = i^(2h-2)"]["U_(pi/4)"], "Ad Z_SU = Ad U_(pi/2)"),
    "Referenz U_(pi/4)": ("det -1: eps -> -eps", "I_pi -> -i I_pi", perm_von("U_(pi/4) (Modul)"), "-", "Modulfluss"),
}
ZAHLEN["zensus"] = {k: dict(zip(("Volumenform", "Phase von I_pi", "Sektoren", "U_(pi/4)", "Bemerkung"), v))
                    for k, v in ZENSUS.items()}
say("   Zensus (Volumenform | Phase von I_pi | Sektoren 0..3 -> | U_(pi/4) | Bemerkung):")
for k_, v_ in ZENSUS.items():
    say(f"     {k_}: " + " | ".join(str(x) for x in v_))
ok_zensus = (abs(eps_J - 1) < 1e-12 and abs(lam_J - 1) < 1e-12 and abs(eps_C - 1) < 1e-12 and abs(lam_C - 1) < 1e-12
             and all(ZENSUS[k][3] == "vertauscht" for k in ("1 Marken, Uhr R = i^(Q_a - Q_b)", "3 tau = Theta_1, Theta_n",
                                                             "8 Zentrum Z_SU = i^(2h-2)")))
check("Zensustabelle aus den gerechneten Größen: J und C wirken auf T-Strukturen über tau bzw. tau K mit Volumenfaktor 1 "
      "und erhalten I_pi; die Daten 1, 3 und 8 vertauschen mit U_(pi/4), die Umkehrungen von 2, 4 und 5 laufen entweder über "
      "tau oder J (Sektoren fest) oder über C bzw. K (Spiegelung, invertieren U)", ok_zensus)
befund("E: Zensus. Keines der acht Orientierungsdaten wirkt auf die Twistor-Volumenform oder wählt einen Sektor: Marken und "
       "Uhr R, Uhrrichtung (Tomita J), KMS-Richtung, Signatur (2,2), Hardy-Orientierung (tau), Theta_n, eta-Dublett und Z_SU "
       "vertauschen mit U_(pi/4); auf den Sektoren wirken sie als Eins (R, tau, J), als Halbdrehung B -> -B (Z_SU, beide "
       "Richtungen) oder paartyperhaltend (Theta_n). Spiegelnd wirken nur die modenerhaltende Konjugation K (Familie D, "
       "R -> R^-1) und der unitäre Teilchen-Loch-Tausch C (Hardy-Umkehr), und beide invertieren den Modulfluss. Neu: Auf "
       "der Spinbrücke hat die Uhr selbst det R = -1 (ebenso R^-1), sie ist eine SL(2) x SL(2)-Drehung mal der "
       "Zentralphase lambda mit lambda^4 = -1, also komplex konform mit Faktor lambda^2 = -i; dort lassen die Marken keine "
       "invariante Volumenform zu, wohl aber die Faktorzerlegung V^* (x) eta", True)
laufzeit()

# =================================================================== Teil 2
say("\nTeil 2: Auswahl eines Sektors oder Paares")

# --- 2.1 Gruppe der nativen Sektorwirkungen und invariante Teilmengen
NATIV = {k: S_TAB[k]["Sektorpermutation (0: e^e, 1, 2, 3)"] for k in
         ("U_(pi/4) (Modul)", "Z_SU = i^(2h-2)", "Z_SU^-1", "tau = (Ad J)^+", "Tomita J", "Uhr R (Ad_R = 1 auf den Oszillatoren)")}


def erzeugt(perms):
    G = {(0, 1, 2, 3)}
    rand = set(G)
    while rand:
        neu = set()
        for g in rand:
            for p in perms:
                h = tuple(p[g[i]] for i in range(4))
                if h not in G:
                    neu.add(h)
        G |= neu
        rand = neu
    return G


def invariante_teilmengen(G):
    return [set(S) for r in range(5) for S in itertools.combinations(range(4), r)
            if all({g[i] for i in S} == set(S) for g in G)]


G_nat = erzeugt(list(NATIV.values()))
inv_nat = invariante_teilmengen(G_nat)
G_ohne = erzeugt([p for k, p in NATIV.items() if k != "U_(pi/4) (Modul)"])
inv_ohne = invariante_teilmengen(G_ohne)
fix_nat = {k: [i for i in range(4) if p[i] == i] for k, p in NATIV.items()}
ZAHLEN["auswahl"] = {"Gruppe nativ (Ordnung)": len(G_nat), "invariante Teilmengen nativ": [sorted(s) for s in inv_nat],
                     "ohne U_(pi/4): Gruppe": sorted(G_ohne), "ohne U_(pi/4): invariante Teilmengen": [sorted(s) for s in inv_ohne],
                     "Fixsektoren je Datum": fix_nat}
check(f"Curie-Prüfung, explizit: Die nativen Sektorwirkungen aus 1.2 erzeugen die zyklische Gruppe der Ordnung "
      f"{len(G_nat)} (U_(pi/4) dreht, Z_SU = Drehung um zwei, R, tau, J trivial); von allen 16 Teilmengen der vier Sektoren "
      f"sind nur {[sorted(s) for s in inv_nat]} invariant. Ohne den Modulfluss blieben {[sorted(s) for s in inv_ohne]}: "
      f"Z_SU (nativ) verbietet einen einzelnen Sektor schon allein, die Paarwahl hängt allein an U_(pi/4)",
      inv_nat == [set(), {0, 1, 2, 3}] and len(G_nat) == 4 and inv_ohne == [set(), {0, 2}, {1, 3}, {0, 1, 2, 3}])

# --- 2.2 alle Kombinationen nativer Strukturen: gemeinsamer Stabilisator in {U_(j pi/4)} und erlaubte Sektormengen
ROT_POT = {j: tuple((i + j) % 4 for i in range(4)) for j in range(8)}   # U_(j pi/4) wirkt als rho^j (1.2)
ok_pot = all(erzeugt([ROT]) >= {ROT_POT[j]} for j in range(8)) and ROT_POT[1] == ROT


def A_th(th):
    return np.exp(-1j*th)*np.eye(4)


TESTS = {
    "KMS-Vektor": lambda th: np.allclose(U1(th) @ OMV1, OMV1),
    "Modulgenerator D (KMS-Richtung)": lambda th: norm_max(U1(th) @ Dop1 @ dag(U1(th)) - Dop1) < 1e-12,
    "CCR-Form G (Signatur)": lambda th: np.allclose(A_th(th).conj().T @ G4 @ A_th(th), G4),
    "Hardy-Projektor P_+": lambda th: np.allclose(A_th(th) @ P_plus, P_plus @ A_th(th)),
    "Uhr R (auf T die Eins, 1.4)": lambda th: np.allclose(A_th(th) @ np.eye(4), np.eye(4) @ A_th(th)),
    "tau": lambda th: np.allclose(A_th(th) @ S4, S4 @ A_th(th)),
    "Theta_n (acht n)": lambda th: all(np.allclose(A_th(th) @ Theta(n), Theta(n) @ A_th(th)) for n in NS),
    "Tomita J": lambda th: norm_max(Jad1(U1(th)) - U1(th)) < 1e-12,
    "Z_SU": lambda th: norm_max(U1(th) @ ZSU1 - ZSU1 @ U1(th)) < 1e-14,
}
STAB = {k: [j for j in range(8) if t(j*np.pi/4)] for k, t in TESTS.items()}
STAB_EPS = [j for j in range(8) if abs(1/np.linalg.det(A_th(j*np.pi/4)) - 1) < 1e-12]
ergebnis_ohne, ergebnis_mit = set(), set()
n_kombi = 0
for r in range(1, len(STAB) + 1):
    for kom in itertools.combinations(STAB, r):
        n_kombi += 1
        js = set.intersection(*[set(STAB[k]) for k in kom])
        ergebnis_ohne.add(tuple(tuple(sorted(x)) for x in invariante_teilmengen(erzeugt([ROT_POT[j] for j in js]))))
        js_e = js & set(STAB_EPS)
        ergebnis_mit.add(tuple(tuple(sorted(x)) for x in invariante_teilmengen(erzeugt([ROT_POT[j] for j in js_e]))))
ZAHLEN["kombinationen"] = {"Stabilisator je Struktur (j: U_(j pi/4))": STAB, "Stabilisator eps_0": STAB_EPS,
                           "Kombinationen": n_kombi, "erlaubte Sektormengen ohne eps": sorted(ergebnis_ohne),
                           "erlaubte Sektormengen mit eps": sorted(ergebnis_mit)}
check(f"Alle {n_kombi} Kombinationen der neun nativen Orientierungsstrukturen (KMS-Vektor, Modulgenerator, G, Hardy-Projektor, "
      f"Uhr, tau, Theta_n, J, Z_SU): jede einzelne ist unter allen U_(j pi/4) fest (Stabilisatoren {set(map(len, STAB.values()))} "
      f"von 8), also auch jede Kombination; die invarianten Sektormengen sind in allen Fällen {sorted(ergebnis_ohne)}. Erst "
      f"mit der Volumenform eps_0 (Stabilisator j in {STAB_EPS}, Ladung 4) werden es {sorted(ergebnis_mit)}: genau die "
      f"zwei Paare", ok_pot and all(v == list(range(8)) for v in STAB.values()) and STAB_EPS == [0, 2, 4, 6]
      and ergebnis_ohne == {((), (0, 1, 2, 3))} and ergebnis_mit == {((), (0, 2), (1, 3), (0, 1, 2, 3))})

# --- 2.3 rahmenabhängige Spiegelungen als Modulbahn
BAHN = {}
ok_bahn = True
for j, th in enumerate((0.0, np.pi/8, np.pi/4, 3*np.pi/8)):
    U_, Ui_ = U2(th), U2(-th)
    Kth = (lambda X, U_=U_, Ui_=Ui_: (U_ @ K2(Ui_ @ X @ dag(Ui_)) @ dag(U_)).tocsr())
    Cth = (lambda X, U_=U_, Ui_=Ui_: (U_ @ C2(Ui_ @ X @ dag(Ui_)) @ dag(U_)).tocsr())
    pk, sk, dk, okk = sektor_wirkung(Kth)
    pc, sc, dc, okc = sektor_wirkung(Cth)
    BAHN[f"theta = {j} pi/8"] = {"K_theta": pk, "C_theta": pc, "Fixsektoren": [i for i in range(4) if pk[i] == i],
                                 "phi -> s phi + delta": (sk, round(dk, 12))}
    ok_bahn &= okk and okc and pk == pc and sk == -1 and abs(np.angle(np.exp(1j*(dk - 4*th)))) < 1e-12
fixe = [v["Fixsektoren"] for v in BAHN.values()]
U_q = U2(np.pi/4)
r_K = norm_max((U_q @ K2(dag(U_q) @ B_kern(0.37) @ U_q) @ dag(U_q)) - K2(B_kern(0.37)))
ok_bahn &= fixe == [[0, 2], [], [1, 3], []] and r_K > 0.1
ZAHLEN["spiegelbahn"] = BAHN
check(f"Rahmenabhängige Spiegelungen: K_theta = U_theta K U_theta^-1 und C_theta wirken gleich, als phi -> -phi + 4 theta; "
      f"für theta = 0, pi/8, pi/4, 3 pi/8 sind die Fixsektoren {fixe}: Die Bahn enthält die Spiegelung, die das "
      f"e^e-Paar festhält, und die, die das *(e^e)-Paar festhält, verbunden durch U_(pi/4) (U K U^-1 != K, Abstand "
      f"{r_K:.3f}). Die Konjugation der Moden (Familie D, Zeitumkehr) und der unitäre Teilchen-Loch-Tausch sind nativ nur "
      f"als diese Bahn bestimmt; der Rahmen ist die Phase der Moden, also die Modulphase", ok_bahn)

# --- 2.4 U_(pi/4)-ungerade Daten: welche Orientierung kippt U_(pi/4)?
UNGERADE = {}
Ut_inner = np.diag(np.exp(1j*(np.pi/4)*(QA - QB)/2))                    # Modultorus im Ladungsteil (z82: J_z = (Q_a - Q_b)/2)
paare = {
    "Uhr R gegen R^-1 (Ladungsteil, Modultorus)": (np.diag(R_f), np.diag(R_f).conj(), "inner_R"),
    "Uhr R gegen R^-1 (auf T)": (np.eye(4, dtype=complex), np.eye(4, dtype=complex), "T"),
    "eta-Spinorhebung gegen andere Hebung": (Rs, np.array(G_ANDERE.tolist(), dtype=complex), "inner_eta"),
    "Z_SU gegen Z_SU^-1": (ZSU1, dag(ZSU1), "fock"),
    "D gegen -D (KMS-Richtung)": (Dop1, -Dop1, "fock"),
    "G gegen -G (Signatur)": (G4, -G4, "T"),
    "P_+ gegen P_- (Hardy)": (P_plus, np.eye(4) - P_plus, "T"),
    "Theta_n gegen Theta_-n": (Theta(NS[0]), Theta(-NS[0]), "T"),
}
for name, (X, Xf, art) in paare.items():
    if art == "T":
        Y = A_U @ X @ np.linalg.inv(A_U)
    elif art == "fock":
        Y = (U1q @ X @ dag(U1q)).toarray() if sps.issparse(X) else U1q @ X @ dag(U1q)
        X, Xf = (X.toarray() if sps.issparse(X) else X), (Xf.toarray() if sps.issparse(Xf) else Xf)
    elif art == "inner_R":
        Y = Ut_inner @ X @ Ut_inner.conj().T
    else:
        Y = expJz(np.pi/4) @ X @ expJz(-np.pi/4)
    gleich, kipp = np.allclose(Y, X), np.allclose(Y, Xf)
    UNGERADE[name] = "fest" if gleich else ("gekippt" if kipp else "verschoben")
eps0 = LC.astype(complex)
eps_bild = eps0/np.linalg.det(A_U)
UNGERADE["Volumenform eps gegen -eps"] = "gekippt" if np.allclose(eps_bild, -eps0) else "fest"
ZAHLEN["U_(pi/4)-Parität der Orientierungen"] = UNGERADE
check(f"Welche Orientierung kippt U_(pi/4) (auf T, im Fockraum, im Ladungsteil über den Modultorus, auf eta über "
      f"e^(i pi J_z/4))? {UNGERADE}: alle nativen Orientierungen sind U_(pi/4)-gerade, ungerade ist allein die "
      f"Volumenform. Kein natives Orientierungsbit ist mit dem Paar (e^e gegen *(e^e)) korreliert",
      all(v == "fest" for k, v in UNGERADE.items() if k != "Volumenform eps gegen -eps")
      and UNGERADE["Volumenform eps gegen -eps"] == "gekippt")

# --- 2.5a bedingte Verklebungen aus z97 (P_A3: R = Z_SU^-1 auf der Leiter; P_GSO: R^2 = U_2pi = U_theta=pi)
perm_U2pi, s_U2pi, d_U2pi, ok_U2pi = sektor_wirkung(AdU(np.pi))
perm_A3 = perm_von("Z_SU^-1")
ok_verkl = ok_U2pi and perm_U2pi == (0, 1, 2, 3) and perm_A3 == (2, 3, 0, 1) and perm_U2pi == ROT_POT[4]
ZAHLEN["bedingte_verklebungen"] = {"P_A3: R = Z_SU^-1 wirkt als": perm_A3, "P_GSO: R^2 = U_2pi wirkt als": perm_U2pi}
check(f"Bedingte Verklebungen (z97) ändern nichts: Unter P_A3 (Kleber R = Z_SU^-1 auf der radialen Leiter) wirkt die Uhr "
      f"auf den Sektoren als {perm_A3} (Halbdrehung B -> -B), unter P_GSO (R^2 = U_2pi, auf dem Fockraum U_theta mit "
      f"theta = pi) als {perm_U2pi} (Eins); beide liegen in <U_(pi/2)>, keine wählt ein Paar", ok_verkl)

# --- 2.5 die Volumenform als Selektor und die minimale Prämisse
SEL = {}
for j in range(8):
    psi = j*np.pi/4
    SEL[f"psi = {j} pi/4"] = [k for k in range(4) if typ(I_k(k), psi) == 1]
SEL["psi = 0,4"] = [k for k in range(4) if typ(I_k(k), 0.4) == 1]
ok_sel = (SEL["psi = 0 pi/4"] == [0, 2] and SEL["psi = 4 pi/4"] == [1, 3]
          and all(SEL[f"psi = {j} pi/4"] == [] for j in (1, 2, 3, 5, 6, 7)) and SEL["psi = 0,4"] == [])
# U_(pi/4) bildet die Auswahl von eps_0 auf die von eps_pi ab; der Stabilisator <U_(pi/2)> wirkt als {1, Z_SU}
bild = sorted(S_TAB["U_(pi/4) (Modul)"]["Sektorpermutation (0: e^e, 1, 2, 3)"][k] for k in SEL["psi = 0 pi/4"])
G_stab = erzeugt([S_TAB["Z_SU = i^(2h-2)"]["Sektorpermutation (0: e^e, 1, 2, 3)"]])
inv_stab = invariante_teilmengen(G_stab)
ok_sel &= bild == SEL["psi = 4 pi/4"] and {0, 2} in inv_stab and {1, 3} in inv_stab and {0} not in inv_stab
STAB16 = [j for j in range(16) if all(t(j*np.pi/8) for t in TESTS.values())]
STAB16_UNOR = [j for j in range(16) if min(abs(1/np.linalg.det(A_th(j*np.pi/8)) - 1),
                                            abs(1/np.linalg.det(A_th(j*np.pi/8)) + 1)) < 1e-12]
unor = sorted(set(SEL["psi = 0 pi/4"]) | set(SEL["psi = 4 pi/4"]))
ok_sel &= STAB16 == list(range(16)) and STAB16_UNOR == list(range(0, 16, 2)) and unor == [0, 1, 2, 3]
ZAHLEN["volumenform_selektor"] = {"reelle Sektoren je psi": SEL, "U_(pi/4)-Bild der Auswahl von eps_0": bild,
                                  "Stabilisator der nativen Strukturen in U_(j pi/8)": STAB16,
                                  "Stabilisator der unorientierten Volumenform {+-eps_0}": STAB16_UNOR,
                                  "invariant unter dem Stabilisator von eps": [sorted(s) for s in inv_stab]}
check(f"Die Volumenform als Selektor: Bezüglich eps_psi (Phase psi) sind die lorentzsch reellen Sektoren {SEL}. Eine "
      f"unorientierte Volumenform (psi mod pi) lässt genau die vier einfachen Sektoren zu, ihre Orientierung wählt das Paar: "
      f"psi = 0 gibt e^e und -e^e, psi = pi gibt +-*(e^e); U_(pi/4) bildet die eine Auswahl auf die andere ab ({bild}). Der "
      f"Stabilisator von eps (<U_(pi/2)>) wirkt als {{1, Z_SU}} und lässt nur die Paare invariant: Ein einzelner Sektor "
      f"(e^e gegen -e^e) bliebe auch dann offen; er ist das Vorzeichen von B, also der Wirkung. Auch die unorientierte "
      f"Volumenform (Ebene der Simplizität, Sektoren {unor}) ist nicht nativ: In U_(j pi/8) ist ihr Stabilisator "
      f"j in {STAB16_UNOR}, der aller nativen Strukturen dagegen alle {len(STAB16)} j", ok_sel)
befund("E: Kein natives Datum und keine native Kombination wählt einen Sektor oder ein Paar. Die native Symmetrie wirkt auf "
       "den vier Sektoren als zyklische Gruppe Z_4 = <U_(pi/4)> (Z_SU = U_(pi/4)^2, alles andere trivial); invariant sind "
       "nur die leere Menge und alle vier. Spiegelungen (Zeitumkehr K, Teilchen-Loch C) gibt es nur als Modulbahn, die "
       "beide Paarfixierungen enthält. Die minimale Prämisse ist genau ein Bit: die Orientierung der Twistor-Volumenform "
       "relativ zur Realität von I (I-quer = +*I), gleichwertig ein Repräsentant der Nebenklasse von <U_(pi/2)> in "
       "<U_(pi/4)>, also der Schritt U(2,2) -> SU(2,2) samt Vorzeichen", True)
laufzeit()

# =================================================================== Teil 3
say("\nTeil 3: Bilanz und Urteil")
STATUS = {
    "Marken und Uhr R": ("E", "R = i^(Q_a - Q_b), auf Strömen und T die Eins; vertauscht mit U_(pi/4); Sektoren fest"),
    "Uhr auf der Spinbrücke": ("E, neu", "det R = det R^-1 = -1 auf H0(N), R = lambda SL(2) x SL(2), lambda^4 = -1; "
                                          "Produktmetrik komplex konform (Faktor -i); keine R-invariante Volumenform, "
                                          "Faktorzerlegung V^* (x) eta erhalten"),
    "Hardy-Orientierung und Signatur (2,2)": ("E", "P_+ und G U-fest; Umkehr durch tau (Sektoren fest) oder C (Spiegelung, "
                                                   "invertiert U, rahmenabhängig)"),
    "tau, Theta_n": ("E", "linear, det 1, vertauschen mit U; tau Sektoren fest, Theta_n paartyperhaltend"),
    "Uhrrichtung R gegen R^-1": ("E", "Tomita J: vertauscht mit U, Sektoren fest; Familie D (K): invertiert U, Spiegelung, "
                                      "rahmenabhängig"),
    "KMS-Richtung beta > 0": ("E", "Vorzeichen von D; Umkehr durch J (trivial) oder C (Spiegelung, rahmenabhängig)"),
    "eta-Dublett": ("E", "eps_eta nur bis aufs Vorzeichen (Deck); keine R- oder modular äquivariante Übertragung auf I_pi"),
    "Zentrum Z_SU": ("E", "Ad Z_SU = Ad U_(pi/2): B -> -B, beide Richtungen gleich; verbietet einen Einzelsektor"),
    "Auswahl nativ": ("E, negativ", "Sektorgruppe Z_4, invariant nur leer/alles; alle Kombinationen U_(pi/4)-fest"),
    "minimale Prämisse": ("O", "P_vol: Orientierung der Twistor-Volumenform relativ zu I (ein Bit), Einzelsektor zusätzlich "
                               "das Vorzeichen der Wirkung"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
urteil = ("H4 ist als gescoptes Negativ beantwortet: Keine native Orientierung wählt den Gravitationssektor. (1) Zensus: Die "
          "acht Orientierungsdaten (Marken mit Uhr R, Hardy-Orientierung, tau und Theta_n, Uhrrichtung, KMS-Richtung, "
          "Signatur (2,2), eta-Dublett, Zentrum Z_SU) vertauschen alle mit U_(pi/4) und lassen die Twistor-Volumenform fest "
          "(det 1 bzw. Wirkung auf einem anderen Faktor). Auf den vier Plebanski-Sektoren wirken sie als Eins (R, tau, Tomita "
          "J), als Halbdrehung B -> -B (Z_SU in beiden Richtungen, gleich Ad U_(pi/2)) oder paartyperhaltend (Theta_n). Die "
          "Uhr R = i^(Q_a - Q_b) lebt im Ladungsteil und ist auf allen Strömen, also auf T, die Eins. Spiegelnd wirken nur "
          "die modenerhaltende Konjugation K (Familie D, R -> R^-1) und der unitäre Teilchen-Loch-Tausch C (Hardy- bzw. "
          "KMS-Umkehr); beide invertieren den Modulfluss, und ihre Rahmenbahn K_theta enthält die e^e-fixierende und die "
          "*(e^e)-fixierende Spiegelung (theta = 0 und pi/4). (2) Auswahl: Die native Sektorgruppe ist Z_4 = <U_(pi/4)>; "
          "von allen 16 Teilmengen sind nur die leere und die volle invariant, und jede der 511 Kombinationen der neun "
          "nativen Orientierungsstrukturen ist U_(pi/4)-fest (Stabilisator 8/8 gegen 4/8 für eps). Von allen Orientierungen kippt U_(pi/4) allein die "
          "Volumenform. Die Volumenform eps_psi ist der Selektor: unorientiert lässt sie genau die vier einfachen Sektoren "
          "zu, orientiert (psi = 0 bzw. pi) genau ein Paar, e^e/-e^e bzw. +-*(e^e). Ein Einzelsektor bleibt selbst dann "
          "offen, weil das native Z_SU = U_(pi/2) im Stabilisator liegt und B -> -B abbildet. Minimale Prämisse: P_vol, die "
          "Orientierung der Twistor-Volumenform relativ zur Realität von I (I-quer = +*I), ein Bit; das ist z98s Reduktion "
          "U(2,2) -> SU(2,2) mit Vorzeichen, keine neue Größe. Auch die unorientierte Volumenform (Ebene der "
          "Simplizität) ist nicht nativ (Stabilisator 8 von 16 in U_(j pi/8)). Neu ist zudem: Auf der Spinbrücke hat die Uhr "
          "selbst det R = det R^-1 = -1 (R = lambda mal SL(2) x SL(2), lambda^4 = -1, Produktmetrik mal -i); dort lassen "
          "die Marken keine invariante Volumenform zu, wohl aber die Faktorzerlegung V^* (x) eta. (3) Folge für K13: Die Plebanski-Sektorwahl reiht sich als P_vol in die Gravitationsprämissen "
          "ein (Dynamik, b bzw. G, dS-Vorzeichen); sie ist kein eigener offener Mechanismus mehr, sondern dasselbe Bit wie die "
          "Volumenform. K13 bleibt L. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)}); Laufzeit "
    f"{time.time() - T0:.0f} s")
(HERE / "z103_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
