"""Z127: Folgeauftrag zu U2 (z118) -- Prämissen des Gravitationswerts 'MacDowell-Mansouri als Spur mit Uhr'.
Folgt P_R aus P_clock und dem KMS-Funktional bei beta = 4 pi? Fixiert die KMS-Richtung das Paarungsbit? Ist das Gewicht
2 c_3 nativ? Welche Prämissen bleiben für die Szenario-Klausel 'graviton'?

Leitlinie 'Rekonstruktion statt Auswahl'. Firewalls: GRAV.NONCIRCULAR.01 (die Uhr R tritt nur als Gewicht bzw. Modulelement
auf, nie als Zeit), GRAV.PARENT.GAP.NO_GO.01 (kein Elternoperator wird vorausgesetzt; D_perp = gamma_- A + gamma_+ A^+ ist
das Testobjekt aus z116/z118), DIMENSION.UPLIFT.FIREWALL.01 (Clifford-, Twistor- und Toeplitzkinematik einer Naht plus
Formenalgebra; keine wechselwirkende 3+1-Schließung). Gelöste Blöcke (Flavour, Massen, alpha, 16, E8) und H6 bleiben
unberührt; GRAV.NONLINEAR.01 (G = v_geo^-2), v60 (rho_Lambda aus alpha) und c_3 = 1/(8 pi) (K5, EM.WARD.02, SEAM.KHALF.01)
werden nur gelesen; nichts wird befördert. Bestand: z118 (Spur mit Uhr, zulässige Gewichte span{Gamma_2, R}, MM mit
c0 : c1 : c2 = 1/4 : 1 : 6, Uhr auf beiden Faktoren s-abhängig, Paarungsbit offen, G nur mit 2 c_3), z117 (J_H6, Schnitt,
Helizität = Gesamtchiralität x Schnitt, U_(pi/2) (x) R = Gamma_ges), z107 (P_tn und P_clock aus A*, beta = 4 pi), z94
(P_clock: Delta^(-i/2) = R, Modulfluss auf T: e^(4 pi i m t)), z104 (Kappe O_1, Spitzen, Wick Ad(e^(+-i pi D/2)) E = +-i Q^-,
Vorzeichen von Lambda = Klasse der Zeit), z98 (MM-Algebra nativ, l = R_Vak/2, |b| = (pi^2/2) e^(2/alpha) als Ziel), z47/z48
(Polgewicht 1 - q^2, Thermofeldleiter q = e^-2pi).
Statuskürzel: E exakt, C bedingt, L Lesart, O offen.

Teil 1  Ist P_R eine Folge? Modulare Gewichte auf S^+ (Original: Uhrgruppe <Delta^(-i/2)>; Vakuumlichtkegel: Dilatation) gegen
        die endlichen Gewichte span{Gamma_2, R}; die Nahtnullmode als KMS-Leiter bei beta = 4 pi (uhrneutral); äquivariante
        Indexspur (Kern, s -> unendlich) gegen Wärmespur bei endlichem s; Wick-Vervollständigung der gemeinsamen Uhr.
Teil 2  Fixiert die KMS-Richtung das Paarungsbit? Paarung x Uhrrichtung -> Kompensator; Kompensator = Wick-Abbildung
        (Punkte x -> -+ i x, Spitzen der Kappe); Vorzeichenstruktur der radialen sl(2,R); die native Halbdrehung
        iota = e^(pi E) als Symmetrie aller KMS-Daten; native Form des Bits und die Konventionen von z117 und z118.
Teil 3  Ist das Gewicht 2 c_3 nativ? Einheitsbetrag modularer Gewichte, KMS-, Ward- und Randgewichte; native Namen von
        2 c_3 (Gauß-Bonnet-Dichte, KMS-Temperatur, Wärmekern); MacDowell-Mansouri bindet Lambda G an das Spurgewicht.
Teil 4  Bilanz: verbleibende Prämissen, OriginLab-Klausel 'graviton', K13.
"""
import itertools
import json
import re
import time
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.linalg as sla
import sympy as sp

HERE = Path(__file__).parent
REPO = HERE.parents[2]
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
rng = np.random.default_rng(127)
mp.mp.dps = 30
QV = mp.exp(-2*mp.pi)
QQ = float(QV)                                                     # q = e^-2pi
ALPHA_INV = sp.Rational("137.035999177")                           # CODATA 2022, nur für die Größenordnung
C3 = 1/(8*sp.pi)


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
    _, s, vh = np.linalg.svd(M, full_matrices=M.shape[0] < M.shape[1])
    rank = int(np.sum(s > tol*max(1.0, s[0]))) if s.size else 0
    return vh[rank:].conj()


def prop(A, B):
    """A = lam B? (lam, Residuum)"""
    a, b = np.ravel(A), np.ravel(B)
    lam = np.vdot(b, a)/np.vdot(b, b)
    return complex(lam), float(np.abs(a - lam*b).max())


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


def konj(G, s):
    """C mit C g* C^-1 = s g für alle g (kleiner Faktor, SVD), wie z116/z118; None, wenn es keine gibt."""
    d = G[0].shape[0]
    rows = [np.kron(g.conj().T, np.eye(d)) - s*np.kron(np.eye(d), g) for g in G]
    _, sv, vh = np.linalg.svd(np.vstack(rows))
    if sv[-1] > 1e-9:
        return None
    C = vh[-1].conj().reshape(d, d, order="F")
    return C/abs(np.linalg.det(C))**(1/d)


# ------------------------------------------------------------------ Twistordaten (wie z98/z104/z118)
I2 = np.eye(2)
SIGMA = [np.eye(2, dtype=complex), sx, sy, sz]
G4 = np.diag([1., 1., -1., -1.]).astype(complex)                     # CCR-Form auf T
S4 = np.block([[np.zeros((2, 2)), I2], [I2, np.zeros((2, 2))]]).astype(complex)   # tau = Theta_1 auf T
uvec = [np.array([1, 0, -1, 0.]), np.array([0, 1, 0, -1.])]
vvec = [np.array([1, 0, 1, 0.]), np.array([0, 1, 0, 1.])]
UU = np.array(uvec).T.astype(complex)
VV = np.array(vvec).T.astype(complex)
Inp = (np.outer(uvec[0], uvec[1]) - np.outer(uvec[1], uvec[0])).astype(complex)       # I_pi (z89)
Ion = (np.outer(vvec[0], vvec[1]) - np.outer(vvec[1], vvec[0])).astype(complex)       # I_omega (z89)
NP = [-0.5*sum(s[i, j]*np.outer(uvec[i], vvec[j]) for i in range(2) for j in range(2)) for s in SIGMA]
NK = [-0.5*sum(s[i, j]*np.outer(vvec[i], uvec[j]) for i in range(2) for j in range(2)) for s in SIGMA]
MP = [1j*N for N in NP]
MK = [1j*N for N in NK]
ETA = np.diag([1., -1., -1., -1.])
MKL = [ETA[m, m]*MK[m] for m in range(4)]
E_T = (MP[0] + MKL[0])/2                                              # kompakte Energie (Uhr des Originals, z104)
QM_T = (MP[0] - MKL[0])/2                                             # Kappengenerator Q^- (z104)
D_T = sla.expm(np.pi/2*E_T) @ QM_T @ sla.expm(-np.pi/2*E_T)           # Dilatation, z104: Ad(e^(pi E/2)) Q^- = D
TMS = np.array([[1, 0, -1, 0], [0, 1, 0, -1], [1, 0, 1, 0], [0, 1, 0, 1]], dtype=complex)   # Inzidenz wie z104
LC = np.zeros((4,)*4)
for _p in itertools.permutations(range(4)):
    LC[_p] = round(np.linalg.det(np.eye(4)[list(_p)]))


def gam(X):
    """Klein-Clifford wie z118: T -> T* über (1/2) eps(X, ., .), T* -> T über X."""
    M = np.zeros((8, 8), complex)
    M[4:, :4] = 0.5*np.einsum("abcd,ab->dc", LC, X)
    M[:4, 4:] = X
    return M


def stern(Q):
    return 0.5*np.einsum("abcd,cd->ab", LC, Q)


def Gam117(Q):
    """Cliffordwirkung in der Konvention von z117: T -> T* über Q, T* -> T über -*Q."""
    Z = np.zeros((4, 4), complex)
    return np.block([[Z, -stern(Q)], [Q, Z]])


def punkt(B):
    """Punkt x^mu einer Ebene (Spalten von B) über die Inzidenz von z104; None = unendlich."""
    wc = TMS @ B
    if abs(np.linalg.det(wc[2:])) < 1e-9:
        return None
    xx = -1j*wc[:2] @ np.linalg.inv(wc[2:])
    return np.array([np.trace(xx @ s)/2 for s in SIGMA])


def ebene_von_punkt(x):
    xx = sum(x[m]*SIGMA[m] for m in range(4))
    return np.linalg.solve(TMS, np.vstack([1j*xx, I2]))


def ebenen(M):
    """Eigenebenen einer 4x4-Matrix: (Eigenwert, Gram-Eigenwerte (2,2)-Form, Punkt, Basis)."""
    w, V = np.linalg.eig(M)
    out = []
    for lam in np.unique(np.round(w, 8)):
        B = np.linalg.qr(V[:, np.abs(w - lam) < 1e-6])[0]
        gram = B.conj().T @ G4 @ B
        out.append((complex(lam), np.round(np.linalg.eigvalsh((gram + gram.conj().T)/2), 8), punkt(B), B))
    return out


# ------------------------------------------------------------------ 256 Dimensionen (z116/z118): konform (8) x innen (32)
Gi = euklid(5)
G11 = chiral(Gi)
I32 = np.eye(32)
E6 = euklid(3)
hh = [E6[0], E6[1], E6[2], E6[3], 1j*E6[4], 1j*E6[5]]                 # Quadrate (+, +, +, +, -, -): 4 = Lorentzzeit, 5 = Zusatzzeit
h7 = chiral(hh)
G16 = [kr(x, I32) for x in hh] + [kr(h7, gi)*(1j if np.allclose(h7 @ h7, np.eye(8)) else 1.0) for gi in Gi]
Gt16 = chiral(G16)
Sp16 = (np.eye(256) + Gt16)/2
gp16, gm16 = (G16[3] + G16[5])/2, (G16[3] - G16[5])/2
Pplus, Pminus = gm16 @ gp16, gp16 @ gm16
Gam2 = Pplus - Pminus
LI = (0, 1, 2, 4)
g5_16 = chiral([G16[a] for a in LI])
Ruhr = kr(np.eye(8), 1j*G11)                                          # Uhr R = i Gamma_11 auf dem inneren Faktor
H7 = kr(h7, I32)
lorpaare = list(itertools.combinations(LI, 2))
GABCD = {(a, b, c, d): G16[a] @ G16[b] @ G16[c] @ G16[d] for (a, b) in lorpaare for (c, d) in lorpaare}
EPS4 = np.zeros((4, 4, 4, 4))
for _p in itertools.permutations(range(4)):
    EPS4[_p] = round(np.linalg.det(np.eye(4)[list(_p)]))
pos = {a: i for i, a in enumerate(LI)}
ETA4 = np.diag([1., 1., 1., -1.])
EPSV = np.array([EPS4[pos[k[0]], pos[k[1]], pos[k[2]], pos[k[3]]] for k in GABCD])
METV = np.array([ETA4[pos[k[0]], pos[k[3]]]*ETA4[pos[k[1]], pos[k[2]]]*(k[0] == k[3])*(k[1] == k[2])
                 - ETA4[pos[k[0]], pos[k[2]]]*ETA4[pos[k[1]], pos[k[3]]]*(k[0] == k[2])*(k[1] == k[3]) for k in GABCD])
SPp, SPm = Sp16 @ Pplus, Sp16 @ Pminus
SXP = {+1: [X @ SPp for X in GABCD.values()], -1: [X @ SPm for X in GABCD.values()]}   # Lorentz-Vierer vertauschen mit P_+-


def spuren(W):
    """Tr_(S^+)(W Y P_+) und Tr_(S^+)(W Y P_-) für die 36 Paare von Lorentz-Bivektoren Y."""
    return (np.array([np.sum(W.T*M) for M in SXP[+1]]), np.array([np.sum(W.T*M) for M in SXP[-1]]))


def endlich(tp, tm, tol=1e-9):
    return float(np.abs(tp + tm).max()) < tol


# Realstruktur J (Ladungskonjugation C16 wie z118)
C10 = {+1: np.eye(32, dtype=complex), -1: np.eye(32, dtype=complex)}
for _j in range(5):
    C10[+1] = C10[+1] @ Gi[2*_j]
    C10[-1] = C10[-1] @ Gi[2*_j + 1]
C16 = None
for _s6 in (1, -1):
    for _s10 in (1, -1):
        _Ca, _Cb = konj(hh, _s6), C10[_s10]
        if _Ca is None:
            continue
        _C = np.kron(_Ca, _Cb)
        for _s in (1, -1):
            if all(np.allclose(_C @ g.conj() @ np.linalg.inv(_C), _s*g) for g in G16):
                if np.allclose(_C @ _C.conj(), np.eye(256)) and np.allclose(_C @ Gt16.conj(), Gt16 @ _C):
                    C16 = _C


# ------------------------------------------------------------------ Hardy-Marken (z118)
def A_rechteck(N, q, p=None):
    """A = S^+ p(S^4) auf C^N -> C^(N-1); p = Koeffizienten von p(z), Standard 1 - q z."""
    p = [1.0, -q] if p is None else p
    S = np.zeros((N + 4*len(p) + 1, N + 4*len(p) + 1))
    for n in range(S.shape[0] - 1):
        S[n + 1, n] = 1
    P = sum(c*np.linalg.matrix_power(S, 4*j) for j, c in enumerate(p))
    return (S.T @ P)[:N - 1, :N]


def index_von(A, tol=1e-12):
    sv = np.linalg.svd(A, compute_uv=False)
    r = int(np.sum(sv > tol))
    return (A.shape[1] - r) - (A.shape[0] - r)


def nullmode(A):
    _, _, vh = np.linalg.svd(A)
    z = vh[-1].conj()
    return z*np.exp(-1j*np.angle(z[np.argmax(np.abs(z))]))


def heat(A, Ws, Wt, s):
    """(Tr(W_s e^(-s A^+A)), Tr(W_t e^(-s AA^+)))"""
    return (complex(np.trace(Ws @ sla.expm(-s*A.conj().T @ A))), complex(np.trace(Wt @ sla.expm(-s*A @ A.conj().T))))


def uhr(n, sgn=1):
    return np.diag([(1j**sgn)**k for k in range(n)])


# =================================================================== Teil 1
say("Teil 1: Ist P_R eine Folge von P_clock und dem KMS-Funktional bei beta = 4 pi?")

# --- 1.1 Zulässige endliche Gewichte (z118 reproduziert)
kand = {"1": np.eye(256), "gamma_5": g5_16, "Gamma_2": Gam2, "R": Ruhr, "Gamma_2 + 0,37 R (Holst)": Gam2 + 0.37*Ruhr,
        "i Gamma_2": 1j*Gam2}
SP = {k: spuren(W) for k, W in kand.items()}
fin = {k: endlich(*v) for k, v in SP.items()}
lamR_g, resR_g = prop(Ruhr @ Sp16, Gam2 @ g5_16 @ Sp16)
lamR_h, resR_h = prop(Ruhr @ Sp16, H7 @ Sp16)
cR, rR = prop(SP["R"][1], EPSV)
cG, rG = prop(SP["Gamma_2"][1], METV)
ZAHLEN["gewichte"] = {"endlich": fin, "R|S+ = lam Gamma_2 gamma_5": str(np.round(lamR_g, 9)), "R|S+ = lam (h7 x 1)": str(np.round(lamR_h, 9)),
                      "Tr(R Y P_-) = c eps": str(np.round(cR, 9)), "Tr(Gamma_2 Y P_-) = c Metrik": str(np.round(cG, 9))}
check(f"[E] Zulässige endliche Gewichte (z118 reproduziert): Die relative Spur Tr_(S^+)(W Y e^(-s D^2)) über die 36 Paare "
      f"von Lorentz-Bivektoren ist endlich (Aufhebung P_+ gegen P_-) für {[k for k, v in fin.items() if v]}, nicht für "
      f"{[k for k, v in fin.items() if not v]}; auf S^+ ist R = {np.round(lamR_g, 6)} Gamma_2 gamma_5 = {np.round(lamR_h, 6)} "
      f"(h7 x 1); Gewicht R gibt {np.round(cR.real, 6)} eps_abcd, Gewicht Gamma_2 {np.round(cG.real, 6)} x Metrikkontraktion "
      f"je Index (Residuen {rR:.0e}, {rG:.0e})",
      fin["Gamma_2"] and fin["R"] and fin["Gamma_2 + 0,37 R (Holst)"] and fin["i Gamma_2"] and not fin["1"]
      and not fin["gamma_5"] and resR_g < 1e-12 and resR_h < 1e-12 and abs(lamR_g + 1j) < 1e-12 and abs(lamR_h - 1j) < 1e-12
      and rR < 1e-9 and rG < 1e-9 and abs(cR - 64) < 1e-9 and abs(cG + 64) < 1e-9)

# --- 1.2 Modulare Gewichte: Uhrgruppe des Originals, Dilatation des Vakuumlichtkegels, beide zusammen
U_konf = {}
for m_ in (1, 2):
    for t_ in (0.5, -0.5, 1.0, 0.25):
        U_konf[(m_, t_)] = sla.expm(4j*np.pi*m_*t_*h7)                 # Modulfluss auf T + T* (z117: e^(4 pi i m t) auf T)
triv_halb = all(np.allclose(U_konf[(m_, t_)], np.eye(8)) for m_ in (1, 2) for t_ in (0.5, -0.5, 1.0))
viertel = prop(U_konf[(1, 0.25)], np.eye(8))
basisW = {}
for k_ in range(4):
    Rk = np.linalg.matrix_power(Ruhr, k_)
    basisW[("R", k_)] = spuren(Rk)
    basisW[("G2R", k_)] = spuren(Gam2 @ Rk)
uhrgruppe = {k_: endlich(*basisW[("R", k_)]) for k_ in range(4)}
treffer = []
for j_ in range(32):
    al = j_*np.pi/16
    for k_ in range(4):
        tp = np.cos(al)*basisW[("R", k_)][0] + 1j*np.sin(al)*basisW[("G2R", k_)][0]
        tm = np.cos(al)*basisW[("R", k_)][1] + 1j*np.sin(al)*basisW[("G2R", k_)][1]
        if endlich(tp, tm):
            treffer.append((j_/16, k_))
erwartet = sorted([(0.0, 1), (0.0, 3), (1.0, 1), (1.0, 3), (0.5, 0), (0.5, 2), (1.5, 0), (1.5, 2)])
vak = {}
for lab_, z_ in (("0,3", 0.3), ("i pi/4", 1j*np.pi/4), ("i pi/2", 1j*np.pi/2), ("-i pi/2", -1j*np.pi/2),
                 ("0,2 + i pi/2", 0.2 + 1j*np.pi/2), ("i pi", 1j*np.pi)):
    vak[lab_] = endlich(*spuren(sla.expm(z_*Gam2)))
werte_mod = {"R": np.round(basisW[("R", 1)][1][np.abs(EPSV) > 0][:3], 6).tolist(),
             "i Gamma_2": np.round(spuren(1j*Gam2)[1][np.abs(METV) > 0][:3], 6).tolist()}
ZAHLEN["modular"] = {"Modulfluss auf T bei t in Z/2 trivial (m = 1, 2)": triv_halb, "t = 1/4, m = 1: lam 1": str(np.round(viertel[0], 6)),
                     "Uhrgruppe R^k endlich": uhrgruppe, "e^(z Gamma_2) endlich": vak,
                     "gemeinsame Gruppe e^(i alpha Gamma_2) R^k endlich bei (alpha/pi, k)": treffer,
                     "Werte": {k: [str(x) for x in v] for k, v in werte_mod.items()}}
check(f"[E, neu] Modulare Auswahl des Gewichts: Der Modulfluss des Originals ist auf T + T* e^(4 pi i m t h7) und bei den "
      f"Uhrzeiten t in Z/2 trivial ({triv_halb}; Kontrolle t = 1/4: {np.round(viertel[0], 3)} x 1), auf dem inneren Faktor "
      f"(Z_4-Ladung) Delta^(it) = R^(-2t). Die Uhrgruppe <Delta^(-i/2)> = {{R^k}} liefert endliche Spuren genau für k ungerade "
      f"({uhrgruppe}), also +-R. Die Dilatationsfamilie des Vakuumlichtkegels e^(z Gamma_2) ist endlich nur bei z = +-i pi/2 "
      f"({vak}), also +-i Gamma_2 (topologischer Sektor). Die gemeinsame Gruppe e^(i alpha Gamma_2) R^k (alpha in pi/16 Z) "
      f"ist endlich genau bei {treffer}: nur die reinen Gewichte +-R und +-i Gamma_2, nie eine Holst-Mischung a Gamma_2 + b R "
      f"mit a b != 0. Die Gewichtswahl (z118: = P_vol-Paarwahl) folgt also aus 'Gewicht = Modulelement des Originals' plus "
      f"Endlichkeit; das Modulelement des Vakuumlichtkegels gäbe den topologischen Sektor",
      triv_halb and abs(viertel[0] + 1) < 1e-12 and viertel[1] < 1e-12 and uhrgruppe == {0: False, 1: True, 2: False, 3: True}
      and sorted(treffer) == erwartet and vak == {"0,3": False, "i pi/4": False, "i pi/2": True, "-i pi/2": True,
                                                "0,2 + i pi/2": False, "i pi": False})

# --- 1.3 Die Nahtnullmode ist die KMS-Leiter bei beta = 4 pi und uhrneutral
NM = 81
Aq81 = A_rechteck(NM, QQ)
zm = nullmode(Aq81)
k_arr = np.arange(NM)
neben = float(np.abs(zm[k_arr % 4 != 0]).max())
gibbs = np.where(k_arr % 4 == 0, (1 - np.exp(-4*np.pi))*np.exp(-np.pi*k_arr), 0.0)
gibbs_dev = float(np.abs(np.abs(zm)**2 - gibbs).max())
R81 = uhr(NM)
R_zm = complex(np.vdot(zm, R81 @ zm))


def delta_halb(beta, n):
    """Delta_beta^(-i/2) = e^(i K/2), K = beta L_0, L_0 = k/4 auf den Hardy-Marken."""
    return np.diag(np.exp(1j*beta*np.arange(n)/8))


ist_R = np.allclose(delta_halb(4*np.pi, NM), R81)
leiter = np.arange(0, NM, 4)
neutral = {f"{b_/np.pi:g} pi": bool(np.allclose(np.diag(delta_halb(b_, NM))[leiter], 1))
           for b_ in (4*np.pi, 8*np.pi, 12*np.pi, 2*np.pi, 6*np.pi, 4.04*np.pi)}
kms_wert = {}
for b_ in (4*np.pi, 2*np.pi, 4.04*np.pi):
    zb = nullmode(A_rechteck(NM, np.exp(-b_/2)))
    kms_wert[f"{b_/np.pi:g} pi"] = complex(np.round(np.vdot(zb, delta_halb(b_, NM) @ zb), 10))
ZAHLEN["nullmode"] = {"max |zm_k| (k nicht = 0 mod 4)": neben, "max | |zm|^2 - Gibbs(beta = 4 pi) |": gibbs_dev,
                      "<zm|R|zm>": str(np.round(R_zm, 12)), "Delta^(-i/2) = R (beta = 4 pi)": ist_R,
                      "Delta_beta^(-i/2) = 1 auf der Leiter |4j>": neutral, "<zm_beta|Delta_beta^(-i/2)|zm_beta>": {k: str(v) for k, v in kms_wert.items()}}
check(f"[E] Die Nahtnullmode ist die KMS-Leiter und uhrneutral: ker A_q (N = {NM}) liegt auf den Marken k = 4j "
      f"(Rest {neben:.0e}), |zm_(4j)|^2 = (1 - e^-4pi) e^(-4 pi j) sind genau die Gibbsgewichte von L_0 = j bei beta = 4 pi "
      f"(Abweichung {gibbs_dev:.0e}): <zm|X|zm> ist das KMS-Funktional auf der uhrneutralen Leiter (z48: Thermofeld). Bei "
      f"beta = 4 pi ist Delta^(-i/2) = R auf den Marken ({ist_R}), und auf der Leiter |4j> ist Delta_beta^(-i/2) genau für "
      f"beta in 4 pi Z die Eins ({neutral}); <zm|R|zm> = {np.round(R_zm.real, 12)}, Kontrollen <zm_b|Delta_b^(-i/2)|zm_b> = "
      f"{ {k: str(v) for k, v in kms_wert.items()} }. A* (z107) wählt beta = 4 pi: Die Uhr fällt auf dem Kern von selbst weg",
      neben < 1e-14 and gibbs_dev < 1e-14 and ist_R and abs(R_zm - 1) < 1e-12
      and neutral == {"4 pi": True, "8 pi": True, "12 pi": True, "2 pi": False, "6 pi": False, "4.04 pi": False}
      and abs(kms_wert["4 pi"] - 1) < 1e-10 and abs(kms_wert["2 pi"] - 1) > 1e-3 and abs(kms_wert["4.04 pi"] - 1) > 1e-8)
laufzeit()

# --- 1.4 Äquivariante Indexspur (Kern) gegen Wärmespur bei endlichem s
NH = 41
Aq41 = A_rechteck(NH, QQ)
Rs41, Rt41 = uhr(NH), uhr(NH - 1)
cRp, cRm = SP["R"]
f_s = {}
for s_ in (0.1, 1.0, 4.0, 30.0):
    ts, tt = heat(Aq41, Rs41, Rt41, s_)
    f_s[s_] = ts - tt                                                  # gemeinsame Uhr: Tr = Tr(R Y P_-) f(s)
spur_gemeinsam = {s_: cRp*heat(Aq41, Rs41, Rt41, s_)[1] + cRm*heat(Aq41, Rs41, Rt41, s_)[0] for s_ in (0.1, 30.0)}
kern_wert = cRm*R_zm                                                   # s -> unendlich: nur der Kern
pr_wert = cRm*1.0                                                      # P_R: R nur innen, Index 1
s_abh = max(abs(f_s[s_] - 1) for s_ in (0.1, 1.0, 4.0)) > 1e-7
grenz = abs(f_s[30.0] - 1) < 1e-12
form = all(abs((f_s[s_] - 1).real - (f_s[s_] - 1).imag) < 1e-12 for s_ in f_s)
kern_gleich = float(np.abs(kern_wert - pr_wert).max()) < 1e-9 and float(np.abs(spur_gemeinsam[30.0] - pr_wert).max()) < 1e-9
z4k = nullmode(A_rechteck(NM, QQ, p=[1.0, -QQ, 0.1]))
lef_kov = complex(np.vdot(z4k, R81 @ z4k))
stoer = {}
for eps_ in (0.3, 0.6):
    A3 = Aq81.copy()
    K = rng.normal(size=(6, 6))
    A3[:6, :6] += eps_*K/np.linalg.norm(K, 2)
    z3 = nullmode(A3)
    stoer[eps_] = {"Index": index_von(A3, 1e-10), "Lefschetz <z|R|z>": complex(np.round(np.vdot(z3, R81 @ z3), 6)),
                   "globale Spur (s = 4)": round(heat(A3, np.eye(NM), np.eye(NM - 1), 4.0)[0].real
                                                 - heat(A3, np.eye(NM), np.eye(NM - 1), 4.0)[1].real, 9)}
ZAHLEN["index_gegen_waerme"] = {"f(s) - 1 (gemeinsame Uhr)": {str(k): str(np.round(v - 1, 12)) for k, v in f_s.items()},
                                "Kernwert = P_R-Wert (36 Paare)": kern_gleich, "Lefschetz, Z_4-kovariante Störung": str(np.round(lef_kov, 12)),
                                "Z_4-brechende Randstörung": {str(k): {kk: str(vv) for kk, vv in v.items()} for k, v in stoer.items()}}
check(f"[E, neu] Indexspur gegen Wärmespur: Mit der gemeinsamen Uhr Delta^(-i/2) = R_innen (x) R_Marken ist die Spur bei "
      f"endlichem s s-abhängig (f(s) - 1 = {[str(np.round(f_s[s_] - 1, 9)) for s_ in (0.1, 1.0, 4.0)]}, Form (1 + i): "
      f"{form}; z118), für s -> unendlich (|f - 1| = {abs(f_s[30.0] - 1):.0e} bei s = 30) aber gleich der Kernspur, und die ist wegen der "
      f"uhrneutralen Nullmode für alle 36 Paare genau der P_R-Wert 64 eps ({kern_gleich}): Die äquivariante Indexspur "
      f"(Lefschetz-Zahl) der gemeinsamen Moduluhr braucht die Einschränkung 'nur innen' nicht. Getragen wird das von der "
      f"Z_4-Homogenität des Nahtschritts: Unter einer Z_4-kovarianten Störung (Symbol 1 - q z + 0,1 z^2 in S^4) bleibt "
      f"<z|R|z> = {np.round(lef_kov.real, 9)}, unter einer Z_4-brechenden Randstörung bleiben Index und globale Spur 1, die "
      f"Lefschetz-Zahl wandert aber ({ {k: str(v['Lefschetz <z|R|z>']) for k, v in stoer.items()} })",
      s_abh and grenz and form and kern_gleich and abs(lef_kov - 1) < 1e-12
      and all(v["Index"] == 1 and abs(v["globale Spur (s = 4)"] - 1) < 1e-8 and abs(v["Lefschetz <z|R|z>"] - 1) > 1e-3
              for v in stoer.values()))

# --- 1.5 Wick-Vervollständigung der gemeinsamen Uhr (z118: Kompensator)
gp8, gm8 = gam(Inp)/2, gam(Ion)/2
Pp8, Pm8 = gm8 @ gp8, gp8 @ gm8
G2_8 = Pp8 - Pm8
Nq = 13
Sq = np.zeros((Nq, Nq))
for _n in range(Nq - 1):
    Sq[_n + 1, _n] = 1
Aquad = Sq.T - QQ*np.linalg.matrix_power(Sq, 3)
Rq = uhr(Nq)
D_minus = np.kron(gm8, Aquad) + np.kron(gp8, Aquad.conj().T)          # z118: Nullmode in P_-
D_plus = np.kron(gp8, Aquad) + np.kron(gm8, Aquad.conj().T)           # umgekehrte Paarung: Nullmode in P_+


def symmetrisch(X8, Rm, D):
    g = np.kron(X8, Rm)
    return float(np.abs(g @ D @ np.linalg.inv(g) - D).max())


X8 = sla.expm(1j*np.pi/4*G2_8)
sym_X = symmetrisch(X8, Rq, D_minus)
sym_1 = symmetrisch(np.eye(8), Rq, D_minus)
sym_Xph = symmetrisch(np.exp(1j*np.pi/4)*X8, Rq, D_minus)
komp = {}
for s_ in (0.1, 1.0, 4.0):
    ts, tt = heat(Aq41, Rs41, Rt41, s_)
    komp[s_] = (np.exp(1j*np.pi/4)*cRp[EPSV != 0][0]*tt + np.exp(-1j*np.pi/4)*cRm[EPSV != 0][0]*ts,
                1j*cRp[EPSV != 0][0]*tt + 1.0*cRm[EPSV != 0][0]*ts)
eps0 = EPSV[EPSV != 0][0]
komp_const = max(abs(v[0] - 64*eps0*np.exp(-1j*np.pi/4)) for v in komp.values())
kompn_const = max(abs(v[1] - 64*eps0) for v in komp.values())
X_T = X8[:4, :4]
nicht_modular = prop(X_T, np.eye(4))[1] > 0.1
ZAHLEN["wick_vervollstaendigung"] = {"|g D g^-1 - D|, X = e^(i pi Gamma_2/4)": sym_X, "X = 1": sym_1, "X e^(i pi/4)": sym_Xph,
                                     "Spur (s = 0,1/1/4), X": [str(np.round(v[0], 9)) for v in komp.values()],
                                     "Spur, X' = e^(i pi P_+/2)": [str(np.round(v[1], 9)) for v in komp.values()],
                                     "X auf T Eigenwerte": [str(np.round(x, 6)) for x in np.linalg.eigvals(X_T)]}
check(f"[E] Wick-Vervollständigung: g = X (x) R_Marken mit X = e^(i pi Gamma_2/4) ist eine Symmetrie von D_perp "
      f"(|g D g^-1 - D| = {sym_X:.0e}; ohne X {sym_1:.2f}; ein globaler Phasenfaktor ändert nichts: {sym_Xph:.0e}). Mit "
      f"g = X R_innen R_Marken ist die Spur für jedes s gleich 64 e^(-i pi/4) eps (Abweichung {komp_const:.0e}), also "
      f"topologisch, aber komplex; erst die Phasennormierung X' = e^(i pi/4) X = e^(i pi P_+/2) (die Eins auf dem "
      f"Nullmodensektor) gibt 64 eps (Abweichung {kompn_const:.0e}), den P_R-Wert. X ist kein Element des Modulflusses des "
      f"Originals (auf T skalar; X hat dort die Eigenwerte e^(-+i pi/4): {nicht_modular}), sondern das imaginäre "
      f"Dilatationselement, [L] die Viertelstelle des Lichtkegelmodulflusses des Vakuums",
      sym_X < 1e-12 and sym_1 > 0.1 and sym_Xph < 1e-12 and komp_const < 1e-9 and kompn_const < 1e-9 and nicht_modular)
befund("E/C: P_R ist zum größten Teil eine Folge. (a) Die Wahl R statt Gamma_2 (und keine Holst-Mischung) folgt aus "
       "'Gewicht = Modulelement des Originals' plus Endlichkeit: Die endlichen modularen Gewichte auf S^+ sind genau +-R "
       "(Uhrgruppe <Delta^(-i/2)>, P_clock); der Vakuumlichtkegel gäbe +-i Gamma_2 (topologisch). (b) 'Nur innen' folgt in "
       "der Indexlesart: Die Nahtnullmode ist bei beta = 4 pi (A*, z107) die KMS-Leiter auf den Marken 4j, auf denen "
       "Delta^(-i/2) = R die Eins ist; die äquivariante Indexspur der gemeinsamen Moduluhr ist daher genau der P_R-Wert "
       "64 eps. (c) In der Wärmespur bei endlichem s bleibt die gemeinsame Uhr s-abhängig (Ordnung e^-4pi); topologisch wird "
       "sie dort nur mit dem Wick-Element X (Vakuumsektor) und einer Phasennormierung. P_R schrumpft so auf zwei Lesarten "
       "des Spurfunktionals: P_mod 'das Gewicht stammt aus dem Modulfluss des Originals' (das ist die Prüfhypothese P_clock + "
       "KMS selbst) und P_idx 'Wirkung = äquivariante Indexspur (s -> unendlich) der Moduluhr'", True)
laufzeit()

# =================================================================== Teil 2
say("\nTeil 2: Fixiert die KMS-Richtung das Paarungsbit?")

# --- 2.1 Paarung x Uhrrichtung -> Kompensator
tab = {}
for dname, D in (("Nullmode in P_-", D_minus), ("Nullmode in P_+", D_plus)):
    for rname, Rm in (("R", Rq), ("R^-1", np.linalg.inv(Rq))):
        tab[f"{dname}, {rname}"] = [round(a_, 4) for a_ in np.linspace(-1, 1, 17)[:-1]
                                    if symmetrisch(sla.expm(1j*a_*np.pi*G2_8), Rm, D) < 1e-12]
alpha_tab = {k: sorted({float(round(((a_ + 0.5) % 1.0) - 0.5, 4)) for a_ in v}) for k, v in tab.items()}
tab = {k: [float(a_) for a_ in v] for k, v in tab.items()}
ZAHLEN["kompensator"] = {"alpha/pi (e^(i alpha Gamma_2) (x) R^rho symmetrisch)": tab, "alpha/pi mod 1": alpha_tab}
check(f"[E] Paarung, Uhrrichtung und Kompensator: e^(i alpha Gamma_2) (x) R^rho ist Symmetrie von D_sigma genau bei "
      f"alpha/pi mod 1 = {alpha_tab}: alpha = rho sigma pi/4 (sigma = +1 für die Nullmode in P_-). Die Richtung des "
      f"imaginären Dilatationselements ist das Produkt aus Paarungsbit und Uhrrichtung; keine der beiden allein",
      alpha_tab == {"Nullmode in P_-, R": [0.25], "Nullmode in P_-, R^-1": [-0.25], "Nullmode in P_+, R": [-0.25],
                    "Nullmode in P_+, R^-1": [0.25]})

# --- 2.2 Kompensator = Wick-Abbildung; Punkte und Spitzen
X_wick = sla.expm(-1j*np.pi/2*D_T)
X_wick_i = sla.expm(1j*np.pi/2*D_T)
xt_gleich = np.allclose(X8[:4, :4], X_wick)
wick_p = prop(X_wick @ E_T @ np.linalg.inv(X_wick), QM_T)
wick_m = prop(X_wick_i @ E_T @ np.linalg.inv(X_wick_i), QM_T)
EB = {n: ebenen(M) for n, M in (("E", E_T), ("Q^-", QM_T), ("D", D_T))}
pE = {str(np.round(l_, 3)): (g_.tolist(), None if p_ is None else complex(np.round(p_[0], 9))) for l_, g_, p_, _ in EB["E"]}
pQ = {str(np.round(l_, 3)): None if p_ is None else complex(np.round(p_[0], 9)) for l_, g_, p_, _ in EB["Q^-"]}
pD = {str(np.round(l_, 3)): None if p_ is None else complex(np.round(p_[0], 9)) for l_, g_, p_, _ in EB["D"]}
c_tr = np.array([0.3, 0.1, -0.2, 0.05])
tr_pkt = punkt(sla.expm(sum(c_tr[m]*MP[m] for m in range(4))) @ VV)
tr_ok = np.allclose(tr_pkt, -c_tr)
x0 = np.array([0.3, 0.1, -0.2, 0.05])
x_wick = punkt(X_wick @ ebene_von_punkt(x0))
punkt_ok = np.allclose(x_wick, -1j*x0)
orig = [B for l_, g_, p_, B in EB["E"] if g_[0] > 0][0]
spitze_X = punkt(X_wick @ orig)
spitze_Xi = punkt(X_wick_i @ orig)
senke = punkt(sla.expm(25.0*QM_T) @ VV)
th = 0.7
bahn = punkt(X_wick @ sla.expm(th*E_T) @ VV)
bahn_kap = punkt(sla.expm(-1j*th*QM_T) @ VV)
ZAHLEN["wick"] = {"X|T = e^(-i pi D/2)": xt_gleich, "Ad(e^(-i pi D/2)) E = lam Q^-": str(np.round(wick_p[0], 9)),
                  "Ad(e^(+i pi D/2)) E = lam Q^-": str(np.round(wick_m[0], 9)), "E-Ebenen (Gram, x^0)": {k: [v[0], str(v[1])] for k, v in pE.items()},
                  "Q^--Ebenen (x^0)": {k: str(v) for k, v in pQ.items()}, "D-Ebenen (x^0)": {k: str(v) for k, v in pD.items()},
                  "X auf x0": [str(np.round(x, 6)) for x in x_wick], "Originalebene -> X / X^-1": [str(np.round(spitze_X[0], 9)), str(np.round(spitze_Xi[0], 9))],
                  "Senke des Q^--Flusses": str(np.round(senke[0], 9))}
check(f"[E] Kompensator = Wick-Abbildung: Auf T ist e^(i pi Gamma_2/4) = e^(-i pi D/2) ({xt_gleich}), Ad(e^(-+i pi D/2)) E = "
      f"{np.round(wick_p[0], 6)} Q^- bzw. {np.round(wick_m[0], 6)} Q^-. Inzidenz (z104; Translation e^(c P) schiebt den "
      f"Ursprung nach -c: {tr_ok}): E-Ebenen {pE} (definit, Vakuumpunkte), Q^--Ebenen {pQ} (Spitzen), D-Ebenen {pD} "
      f"(0 und unendlich). X wirkt auf Punkte als x -> -i x ({punkt_ok}: Wick-Drehung t -> -i t in der Karte mit P_mu als "
      f"Translationen; in der Karte mit K_mu, also um unendlich, ist es t -> +i t), die Uhrbahn geht in die Kappenbahn bei "
      f"u = -i theta über ({np.allclose(bahn, bahn_kap)}), und X "
      f"bildet die Originalebene (Gram +1) auf die Spitze t = {np.round(spitze_X[0].real, 6)} ab, die Senke des positiven "
      f"Q^--Flusses ({np.round(senke[0].real, 6)}); X^-1 auf t = {np.round(spitze_Xi[0].real, 6)} (Quelle). Paarungsbit = "
      f"Wick-Richtung = Spitze der Originalebene (bei fester Uhrrichtung)",
      xt_gleich and wick_p[1] < 1e-12 and wick_m[1] < 1e-12 and abs(wick_p[0] + 1j) < 1e-12 and abs(wick_m[0] - 1j) < 1e-12
      and tr_ok and punkt_ok and np.allclose(bahn, bahn_kap)
      and sorted((v[1].imag for v in pE.values())) == [-1, 1] and sorted(v.real for v in pQ.values()) == [-1, 1]
      and sum(v is None for v in pD.values()) == 1 and all(abs(v) < 1e-9 for v in pD.values() if v is not None)
      and abs(spitze_X[0] + 1) < 1e-9 and abs(spitze_Xi[0] - 1) < 1e-9
      and abs(senke[0] + 1) < 1e-6)

# --- 2.3 Was die KMS-Richtung festlegt
teil = {b_: [float(np.sum(np.exp(-b_*np.arange(n_)/4))) for n_ in (40, 80)] for b_ in (4*np.pi, -4*np.pi)}
konv = abs(teil[4*np.pi][1] - 1/(1 - np.exp(-np.pi))) < 1e-12 and teil[-4*np.pi][1] > 1e100
R_pos = np.allclose(delta_halb(4*np.pi, Nq), Rq)
R_neg = np.allclose(delta_halb(-4*np.pi, Nq), np.linalg.inv(Rq))
ker_A = index_von(A_rechteck(25, QQ))
beide = symmetrisch(sla.expm(1j*np.pi/4*G2_8), Rq, D_minus) < 1e-12 and symmetrisch(sla.expm(-1j*np.pi/4*G2_8), Rq, D_plus) < 1e-12
ZAHLEN["kms_richtung"] = {"Hardy-Zustandssumme (N = 40, 400), beta = +-4 pi": {str(round(k/np.pi)): [str(x) for x in v] for k, v in teil.items()},
                          "Delta^(-i/2) = R (beta = 4 pi)": R_pos, "= R^-1 (beta = -4 pi)": R_neg, "Index des Absenkschritts": ker_A,
                          "beide Paarungen mit eigenem Kompensator bei beta > 0": beide}
check(f"[E] Was die KMS-Richtung festlegt: Auf den Hardy-Marken k >= 0 konvergiert Tr e^(-beta L_0) genau für beta > 0 "
      f"({konv}); dann ist Delta^(-i/2) = R ({R_pos}; bei beta = -4 pi R^-1: {R_neg}), und der Absenkschritt A trägt den "
      f"Index +{ker_A} (Kern = KMS-Leiter). Das fixiert die Uhrrichtung rho, nicht sigma: Bei beta > 0 sind beide Paarungen "
      f"mit ihrem eigenen Kompensator e^(+-i pi Gamma_2/4) symmetrisch ({beide})",
      konv and R_pos and R_neg and ker_A == 1 and beide)

# --- 2.4 Vorzeichenstruktur der radialen sl(2,R) und native Involutionen
def br(A, B):
    return A @ B - B @ A


k1, r1 = prop(br(E_T, QM_T), D_T)
k2, r2 = prop(br(D_T, E_T), QM_T)
k3, r3 = prop(br(QM_T, D_T), E_T)
muster = [sg for sg in itertools.product((1, -1), repeat=3)
          if np.allclose(br(sg[0]*E_T, sg[1]*QM_T), sg[2]*k1*D_T) and np.allclose(br(sg[2]*D_T, sg[0]*E_T), sg[1]*k2*QM_T)
          and np.allclose(br(sg[1]*QM_T, sg[2]*D_T), sg[0]*k3*E_T)]
iota_T = sla.expm(np.pi*E_T)


def zeichen(g):
    return tuple(int(np.round(prop(g @ M @ np.linalg.inv(g), M)[0].real)) for M in (E_T, QM_T, D_T))


inv_tab = {"iota = e^(pi E)": zeichen(iota_T), "tau = S4": zeichen(S4), "iota tau": zeichen(iota_T @ S4)}
iota_G = float(np.abs(iota_T.conj().T @ G4 @ iota_T - G4).max())
iota_det = complex(np.linalg.det(iota_T))
iota_tau = prop(iota_T @ S4 @ np.linalg.inv(iota_T), S4)
iota_uv = np.linalg.matrix_rank(np.hstack([iota_T @ UU, VV]), 1e-9) == 2 and np.linalg.matrix_rank(np.hstack([iota_T @ VV, UU]), 1e-9) == 2
iota_X = np.allclose(iota_T @ X_wick @ np.linalg.inv(iota_T), X_wick_i)
ZAHLEN["vorzeichen_sl2"] = {"[E,Q^-] = k D, [D,E] = k Q^-, [Q^-,D] = k E": [str(np.round(k1, 9)), str(np.round(k2, 9)), str(np.round(k3, 9))],
                            "erlaubte Muster (s_E, s_Q, s_D)": muster, "Involutionen": inv_tab,
                            "iota: G-unitär, det": [iota_G, str(np.round(iota_det, 9))], "iota tau iota^-1 = lam tau": str(np.round(iota_tau[0], 9)),
                            "iota tauscht span u, span v": iota_uv, "iota X iota^-1 = X^-1": iota_X}
check(f"[E] Vorzeichenstruktur: [E, Q^-] = {np.round(k1.real, 6)} D, [D, E] = {np.round(k2.real, 6)} Q^-, [Q^-, D] = "
      f"{np.round(k3.real, 6)} E; Vorzeichenwechsel (s_E, s_Q, s_D) sind genau dann Automorphismen, wenn s_E s_Q s_D = +1 "
      f"({muster}). Realisiert durch native Involutionen {inv_tab}. Mit s_E = KMS-Richtung (Uhr), s_Q = Orientierung der "
      f"Kappe (Spitze), s_D = Paarung gilt also Paarung = KMS-Richtung x Kappenorientierung. iota = e^(pi E) ist G-unitär "
      f"(Rest {iota_G:.0e}, det {np.round(iota_det.real, 6)}), vertauscht mit E, kehrt tau um (iota tau iota^-1 = "
      f"{np.round(iota_tau[0].real, 6)} tau), tauscht die D-Ebenen 0 und unendlich ({iota_uv}) und X mit X^-1 ({iota_X})",
      r1 < 1e-12 and r2 < 1e-12 and r3 < 1e-12 and sorted(muster) == sorted([(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)])
      and inv_tab == {"iota = e^(pi E)": (1, -1, -1), "tau = S4": (-1, -1, 1), "iota tau": (-1, 1, -1)}
      and iota_G < 1e-12 and abs(iota_det - 1) < 1e-12 and abs(iota_tau[0] + 1) < 1e-12 and iota_uv and iota_X)
laufzeit()

# --- 2.5 iota ist Symmetrie aller KMS-Daten und kehrt die Paarung um
iota = G16[4] @ G16[5]                                                # Drehung um pi in der Ebene der beiden Zeitrichtungen
E_cl = G16[4] @ G16[5]/2                                              # Uhrgenerator (Spin-Erzeuger in der (4,5)-Ebene)


def vt(A, B, s=1):
    return float(np.abs(A @ B - s*B @ A).max())


iota_rel = {"[iota, R]": vt(iota, Ruhr), "[iota, Gamma_ges]": vt(iota, Gt16), "[iota, h7 x 1] (Modul-U(1))": vt(iota, H7),
            "[iota, E_Uhr]": vt(iota, E_cl), "{iota, Gamma_2}": vt(iota, Gam2, -1), "{iota, gamma_5}": vt(iota, g5_16, -1),
            "J: C iota* C^-1 - iota": float(np.abs(C16 @ iota.conj() @ np.linalg.inv(C16) - iota).max()),
            "Raumdrehungen": max(vt(iota, G16[a] @ G16[b]) for a, b in ((0, 1), (0, 2), (1, 2))),
            "Boosts (antivertauschen)": max(vt(iota, G16[a] @ G16[4], -1) for a in (0, 1, 2))}
gtausch = float(np.abs(iota @ gp16 @ np.linalg.inv(iota) - gm16).max())
io8 = sla.expm(np.pi*sla.block_diag(E_T, -E_T.T))
D_flip = prop(np.kron(io8, np.eye(Nq)) @ D_minus @ np.linalg.inv(np.kron(io8, np.eye(Nq))), D_plus)
YE = GABCD[(0, 1, 2, 4)]
eps_m = complex(np.sum((Ruhr @ Sp16).T*(YE @ Pminus)))
eps_p = complex(np.sum((Ruhr @ Sp16).T*(YE @ Pplus)))
Rplus = (np.eye(256) - 1j*Ruhr)/2                                     # Projektor auf R = +i (Gamma_11 = +1)
heli = {}
for lab, P_ in (("P_- (z118)", Pminus), ("P_+", Pplus)):
    Q_ = Sp16 @ P_ @ Rplus
    heli[lab] = complex(np.round(prop(g5_16 @ Q_, Q_)[0], 9))
ZAHLEN["iota_256"] = {"Relationen": iota_rel, "iota gamma_+ iota^-1 - gamma_-": gtausch, "iota_8 D_- iota_8^-1 = lam D_+": [str(np.round(D_flip[0], 9)), D_flip[1]],
                      "Tr(R gamma_0124 P_-), (P_+)": [str(eps_m), str(eps_p)], "gamma_5 der R = +i-Komponente im Nullmodensektor": {k: str(v) for k, v in heli.items()}}
check(f"[E, negativ] Die Halbdrehung iota = e^(pi E) (in 256 Dimensionen Gamma_4 Gamma_5, Drehung um pi in der Ebene der "
      f"Lorentzzeit und der Zusatzzeit) erhält alle KMS-Daten: Uhr R, S^+, Modul-U(1), Uhrgenerator, Realstruktur J und die "
      f"Raumdrehungen (Residuen {[f'{v:.0e}' for k, v in iota_rel.items() if not k.startswith('{') and 'Boosts' not in k]}); "
      f"sie kehrt Gamma_2, gamma_5 und die Boosts um, tauscht gamma_+ und gamma_- ({gtausch:.0e}) und bildet D_perp mit "
      f"Nullmode in P_- auf {np.round(D_flip[0].real, 6)} x D_perp mit Nullmode in P_+ ab (Rest {D_flip[1]:.0e}). Im festen "
      f"Rahmen kippt der eps-Term ({np.round(eps_m.real, 3)} -> {np.round(eps_p.real, 3)}) und mit ihm die Helizität der "
      f"R = +i-Materie im Nullmodensektor ({heli}). Da die Marken, A_q, beta > 0 und R unberührt bleiben, kann die "
      f"KMS-Richtung das Paarungsbit nicht festlegen",
      all(v < 1e-12 for k, v in iota_rel.items()) and gtausch < 1e-12 and abs(D_flip[0] + 1) < 1e-12 and D_flip[1] < 1e-12
      and abs(eps_m - 64) < 1e-9 and abs(eps_p + 64) < 1e-9 and heli == {"P_- (z118)": 1 + 0j, "P_+": -1 + 0j})

# --- 2.6 Native Form des Bits: tau-Parität des Schnitts; Konventionen von z117 und z118
schnitt = {}
for name, G in (("z118 (T -> T* über *X)", gam), ("z117 (T -> T* über Q)", Gam117)):
    gp_, gm_ = G(Inp)/2, G(Ion)/2
    nullpaar = np.allclose(gp_ @ gp_, 0) and np.allclose(gm_ @ gm_, 0) and np.allclose(gp_ @ gm_ + gm_ @ gp_, np.eye(8))
    PT = (gp_ @ gm_)[:4, :4]                                           # Nullmodensektor von gamma_- A + gamma_+ A^+
    auf_u = np.allclose(PT @ UU, UU) and np.allclose(PT @ VV, 0)
    auf_v = np.allclose(PT @ VV, VV) and np.allclose(PT @ UU, 0)
    schnitt[name] = {"Nullpaar": nullpaar, "Schnitt": "span u" if auf_u else ("span v" if auf_v else "?"),
                     "tau-Parität": int(np.round(prop(S4 @ (UU if auf_u else VV), UU if auf_u else VV)[0].real))}
rel_kl = [prop(gam(X), Gam117(stern(X))) for X in [np.outer(np.eye(4)[i], np.eye(4)[j]) - np.outer(np.eye(4)[j], np.eye(4)[i])
                                                    for i, j in itertools.combinations(range(4), 2)]]
stern_I = prop(stern(Inp), Ion)
M6 = np.zeros((6, 6), complex)
_paare6 = list(itertools.combinations(range(4), 2))
for _k, (_i, _j) in enumerate(_paare6):
    _X = np.zeros((4, 4), complex)
    _X[_i, _j], _X[_j, _i] = 1, -1
    _Y = stern(_X)
    M6[:, _k] = [_Y[a, b] for a, b in _paare6]
det_stern = complex(np.round(np.linalg.det(M6), 9))
ZAHLEN["schnitt_konventionen"] = {"Schnitt je Konvention": schnitt, "gam_z118(X) = lam Gam_z117(*X)": sorted({str(np.round(r[0], 9)) for r in rel_kl}),
                                  "*I_pi = lam I_omega": str(np.round(stern_I[0], 9)), "det(*) auf Lambda^2 T": str(det_stern)}
check(f"[E] Native Form des Bits und die Konventionen im Bestand: tau = Theta_1 unterscheidet die beiden D-Ebenen (span u: "
      f"tau = -1, Punkt unendlich; span v: tau = +1, Ursprung), und iota kehrt tau um (2.4). Dieselbe Formel D_perp = "
      f"gamma_- A + gamma_+ A^+ mit gamma_+ = Gamma(I_pi)/2 legt den Schnitt in den beiden Cliffordkonventionen des Bestands "
      f"verschieden: {schnitt}. Die Konventionen hängen über den Hodge-Stern zusammen (gam_z118(X) = "
      f"{sorted({str(np.round(r[0].real, 6)) for r in rel_kl})} Gam_z117(*X), Rest {max(r[1] for r in rel_kl):.0e}; *I_pi = "
      f"{np.round(stern_I[0].real, 6)} I_omega), und * ist auf Lambda^2 T = R^(4,2) uneigentlich (det = "
      f"{np.round(det_stern.real, 6)}). Das Paarungsbit ist also die tau-Parität des Schnitts bzw. die Orientierung der "
      f"Dilatationsebene R^(1,1); z117 und z118 haben es stillschweigend entgegengesetzt gesetzt",
      all(v["Nullpaar"] for v in schnitt.values()) and schnitt["z118 (T -> T* über *X)"]["tau-Parität"] == -1
      and schnitt["z117 (T -> T* über Q)"]["tau-Parität"] == 1 and all(abs(r[0] + 1) < 1e-12 and r[1] < 1e-12 for r in rel_kl)
      and abs(stern_I[0] - 1) < 1e-12 and stern_I[1] < 1e-12 and abs(det_stern + 1) < 1e-9)
befund("E/negativ: Die KMS-Richtung fixiert das Paarungsbit nicht. Sie fixiert die Uhrrichtung (Delta^(-i/2) = R, nicht R^-1) "
       "und den indextragenden Absenkschritt; die Richtung des Kompensators ist Paarung x Uhrrichtung, und der Kompensator ist "
       "die Wick-Abbildung, die die Originalebene auf eine Spitze der Kappe legt. Die Vorzeichenstruktur der radialen sl(2,R) "
       "gibt Paarung = KMS-Richtung x Kappenorientierung; die native Halbdrehung iota = e^(pi E) erhält beta > 0, R, S^+, J, "
       "Modul-U(1), Marken und A_q und kehrt Paarung, Wick-Richtung, Spitze, Lorentzzeit, gamma_5 und das Vorzeichen des "
       "eps-Terms gemeinsam um. Das Bit ist nativ ausdrückbar als tau-Parität des Schnitts (gleichwertig: Spitze der "
       "Originalebene, Orientierung der Dilatationsebene, Helizität der R = +i-Materie im Nullmodensektor bei festem S^+); "
       "z117 und z118 haben es über ihre Cliffordkonvention entgegengesetzt gesetzt", True)
laufzeit()

# =================================================================== Teil 3
say("\nTeil 3: Ist das Gewicht 2 c_3 nativ?")

# --- 3.1 Einheitsbetrag modularer Gewichte, KMS-, Ward- und Randgewichte
werte = np.concatenate([basisW[("R", 1)][1], basisW[("R", 3)][1], spuren(1j*Gam2)[1], spuren(-1j*Gam2)[1]])
nz = werte[np.abs(werte) > 1e-9]
einheit = bool(np.allclose(np.abs(nz), 64)) and all(min(abs(x/64 - u) for u in (1, 1j, -1, -1j)) < 1e-9 for x in nz)
w_grid = np.linspace(0, 10, 201)
fam = np.array([np.vdot(zm, np.diag(np.exp(-w_*np.pi*k_arr)) @ zm).real for w_ in w_grid])
fam_lo, fam_hi = float(fam.min()), float(fam.max())
ward = float(abs(zm[0])**2)
zwei_c3 = float(1/(4*np.pi))


def w_rand(s):
    return float(sla.expm(-s*Aq41.T @ Aq41)[0, 0] - sla.expm(-s*Aq41 @ Aq41.T)[0, 0])


lo, hi = 1e-4, 1.0
for _ in range(60):
    mid = (lo + hi)/2
    lo, hi = (mid, hi) if w_rand(mid) < zwei_c3 else (lo, mid)
s_stern = (lo + hi)/2
ZAHLEN["gewichte_nativ"] = {"modulare Werte |.|/64 und Phasen": einheit, "KMS-Familie <zm|Delta^w|zm>, w in [0,10]": [fam_lo, fam_hi],
                            "Ward 1 - q^2 = |zm_0|^2": ward, "2 c_3": zwei_c3, "s* (Randplatz = 2 c_3)": s_stern, "-ln(1 - 2 c_3)": float(-np.log(1 - zwei_c3))}
check(f"[E] Kein natives Spurgewicht hat den Wert 2 c_3: Alle endlichen modularen Gewichte (+-R, +-i Gamma_2) geben Werte "
      f"vom Betrag 64 mit Phasen in {{1, i, -1, -i}} je Index ({einheit}); die KMS-Erwartung der Uhr auf dem Kern ist 1 (1.3); "
      f"die KMS-Familie <zm|Delta^w|zm> (w >= 0) liegt in [{fam_lo:.8f}, {fam_hi:.1f}] = [1 - q^2, 1], das Ward-Gewicht "
      f"1 - q^2 = {ward:.11f} ist die KMS-Wahrscheinlichkeit des Grundzustands. 2 c_3 = {zwei_c3:.6f} liegt um den Faktor "
      f"{fam_lo/zwei_c3:.3f} = 4 pi (1 - q^2) darunter. Eine lokale Randspur erreicht 2 c_3 nur bei der Abschneidung "
      f"s* = {s_stern:.6f} = -ln(1 - 2 c_3) + O(q^2) ({-np.log(1 - zwei_c3):.6f}): Der Zielwert wird über s eingesetzt",
      einheit and abs(fam_hi - 1) < 1e-12 and abs(fam_lo - (1 - QQ**2)) < 1e-9 and abs(ward - (1 - QQ**2)) < 1e-12
      and fam_lo/zwei_c3 > 12 and abs(s_stern + np.log(1 - zwei_c3)) < 1e-4)

# --- 3.2 Native Namen von 2 c_3
Z2 = 2
chi = 2
kegel = 4*sp.pi                                                       # vier Kegelpunkte mit Winkel pi, Defizit je pi
GB = sp.simplify(kegel - 2*sp.pi*chi)
c3_GB = sp.simplify(1/(Z2*kegel))
dirac = sp.simplify(2*C3*kegel)
mw = sp.simplify(C3*kegel)
beta, beta_p = 4*sp.pi, 8*sp.pi
t_ = sp.Symbol("t", positive=True)
waermekern = sp.simplify((4*sp.pi*t_)**-1*t_)
F16 = kr(np.eye(8), (np.eye(32) + G11)/2)
Pi_a = Sp16 @ Pminus @ F16
Pi_b = Sp16 @ Pminus @ (np.eye(256) - F16)
rg = (int(np.linalg.matrix_rank(Pi_a)), int(np.linalg.matrix_rank(Pi_b)))
J_tauscht = np.allclose(C16 @ Pi_a.conj() @ np.linalg.inv(C16), Pi_b)
halb = (complex(np.trace(Ruhr @ YE @ Pi_a)), complex(np.trace(Ruhr @ YE @ Pi_b)))
ZAHLEN["c3_namen"] = {"Gauss-Bonnet: Kegeldefizite - 2 pi chi": str(GB), "c_3 = 1/(|Z_2| oint K)": str(c3_GB), "2 c_3 oint K (Dirac)": str(dirac),
                      "c_3 oint K (Majorana-Weyl)": str(mw), "2 c_3 = 1/beta": str(sp.simplify(2*C3 - 1/beta)), "c_3 = 1/beta'": str(sp.simplify(C3 - 1/beta_p)),
                      "Wärmekern (4 pi t)^-1 t": str(waermekern), "J-Summanden Ränge": rg, "J tauscht": J_tauscht, "eps je Summand": [str(h) for h in halb]}
check(f"[E/L] Native Namen von 2 c_3 (sympy): Auf der Naht S^2(2,2,2,2) sind die vier Kegeldefizite zusammen 4 pi = "
      f"2 pi chi(S^2) (Rest {GB}), c_3 = 1/(|Z_2| oint K) = {c3_GB}, also 2 c_3 oint K = {dirac} (Einheitsindex, beide "
      f"Blätter) und c_3 oint K = {mw} (ein Blatt, Majorana-Weyl); die Zahl 2 = |Z_2| ist die Zahl der J-verbundenen "
      f"Summanden von S^+ cap P_- (Ränge {rg}, J tauscht: {J_tauscht}, je {np.round(halb[0].real, 3)} des eps-Terms). "
      f"Zugleich ist 2 c_3 = 1/beta (KMS-Temperatur des Originals, Rest {sp.simplify(2*C3 - 1/beta)}), c_3 = 1/beta' "
      f"(Produktlesart beta' = 8 pi, z104) und 2 c_3 = (4 pi t)^-1 t ({waermekern}, 2D-Wärmekern). Lesart [L]: 2 c_3 ist "
      f"eine Dichte je Krümmungseinheit bzw. eine Temperatur; die globale Spur sieht davon nur das Integral 1",
      GB == 0 and sp.simplify(c3_GB - C3) == 0 and dirac == 1 and mw == sp.Rational(1, 2) and sp.simplify(2*C3 - 1/beta) == 0
      and sp.simplify(C3 - 1/beta_p) == 0 and sp.simplify(waermekern - 2*C3) == 0 and rg == (32, 32) and J_tauscht
      and abs(halb[0] - halb[1]) < 1e-9 and abs(halb[0] + halb[1] - 64) < 1e-9)
laufzeit()

# --- 3.3 MacDowell-Mansouri bindet Lambda G an das Spurgewicht
TfullR = np.zeros((4, 4, 4, 4))
for (a, b, c, d), v in zip(GABCD.keys(), SP["R"][1]):
    ia, ib, ic, id_ = pos[a], pos[b], pos[c], pos[d]
    for sab, (x1, x2) in ((1, (ia, ib)), (-1, (ib, ia))):
        for scd, (x3, x4) in ((1, (ic, id_)), (-1, (id_, ic))):
            TfullR[x1, x2, x3, x4] = sab*scd*v.real


def kn(S1, S2):
    return (np.einsum("ac,bd->abcd", S1, S2) + np.einsum("bd,ac->abcd", S1, S2) - np.einsum("ad,bc->abcd", S1, S2)
            - np.einsum("bc,ad->abcd", S1, S2))


def dens(Af, Bf):
    return np.einsum("mnrs,mn,rs->", LC, Af, Bf)/4


def lagr(T, Rf, EE, lam):
    F = Rf + lam*EE
    return sum(T[a, b, c, d]*dens(F[a, b], F[c, d]) for a, b, c, d in itertools.product(range(4), repeat=4)
               if T[a, b, c, d] != 0)/4


cs = []
for _ in range(4):
    Rabcd = sum(kn(S, S)/2 for S in [(lambda M: M + M.T)(rng.normal(size=(4, 4))) for _ in range(3)])
    Rup = np.einsum("ai,bj,ijcd->abcd", ETA4, ETA4, Rabcd)
    Ric = np.einsum("ac,abcd->bd", ETA4, Rabcd)
    Rsk = float(np.einsum("bd,bd->", ETA4, Ric))
    Eul = float(np.einsum("abcd,ijkl,ai,bj,ck,dl->", Rabcd, Rabcd, ETA4, ETA4, ETA4, ETA4)
                - 4*np.einsum("ab,ab->", Ric, ETA4 @ Ric @ ETA4) + Rsk**2)
    e = rng.normal(size=(4, 4)) + 2*np.eye(4)
    dete = float(np.linalg.det(e))
    Rf = np.einsum("abcd,cm,dn->abmn", Rup, e, e)
    EE = np.einsum("am,bn->abmn", e, e) - np.einsum("an,bm->abmn", e, e)
    L0, Lp, Lm = lagr(TfullR, Rf, EE, 0.0), lagr(TfullR, Rf, EE, 1.0), lagr(TfullR, Rf, EE, -1.0)
    cs.append((L0/(Eul*dete), (Lp - Lm)/2/(Rsk*dete), ((Lp + Lm)/2 - L0)/dete))
cs = np.array(cs)
c0n, c1n, c2n = (int(round(float(np.mean(cs[:, j])))) for j in range(3))
streu = float(np.abs(cs - cs.mean(axis=0)).max())
kap, lam_s, G_s, Rv, v_s = sp.symbols("kappa lambda G R_Vak v", positive=True)
G_von = 1/(16*sp.pi*kap*c1n*lam_s)                                    # kappa c1 lambda = 1/(16 pi G)
Lam_von = -sp.Integer(c2n)*lam_s/(2*c1n)
LG = sp.simplify(Lam_von*G_von)
GR = sp.simplify(G_von.subs(lam_s, 4/Rv**2)/Rv**2)                    # natives lambda = sigma/l^2, l = R_Vak/2 (z98)
nativ = {"Index (1)": sp.Integer(1), "Ward (1 - q^2)": 1 - sp.exp(-4*sp.pi), "2 c_3": 2*C3}
tab_LG = {k: [float(sp.N(LG.subs(kap, w_), 12)), float(sp.N(GR.subs(kap, w_), 12))] for k, w_ in nativ.items()}
ziel = sp.exp(-2*ALPHA_INV)/(128*sp.pi**3)                            # z98/v60: G/R_Vak^2 (nur gelesen)
k_noetig = sp.simplify(sp.solve(sp.Eq(GR, ziel), kap)[0])
log_k = float(sp.N(sp.log(k_noetig, 10), 12))
R_z118 = sp.solve(sp.Eq((2*C3/4)*v_s**2, 1/(16*sp.pi*G_von.subs({kap: 2*C3, lam_s: 4/Rv**2}))), Rv)
R_z118_alt = sp.solve(sp.Eq((2*C3/4)*v_s**2, 1/(16*sp.pi*G_von.subs({kap: 2*C3/c1n, lam_s: 4/Rv**2}))), Rv)
LG_z118 = sp.simplify(LG.subs(kap, 2*C3))
LG_ledger = float(sp.N(8*sp.pi*(sp.Rational(3, 4)/sp.pi**2)*sp.exp(-2*ALPHA_INV)/(64*sp.pi**2), 12))
ZAHLEN["mm_normierung"] = {"c0, c1, c2 (Streuung)": [c0n, c1n, c2n, streu], "Lambda G": str(LG), "G/R_Vak^2": str(GR),
                           "natives kappa: [Lambda G, G/R_Vak^2]": tab_LG, "kappa nötig für das Ziel (z98, v60)": str(k_noetig),
                           "log10 kappa nötig": log_k, "z118-Normierung mit kappa = 2 c_3: R_Vak": [str(x) for x in R_z118],
                           "dieselbe, 2 c_3 als Gesamtkoeffizient: R_Vak": [str(x) for x in R_z118_alt],
                           "Lambda G (z118: Lambda = -12/R_Vak^2, G = v^-2)": str(LG_z118), "Lambda G (Ledger, dS)": LG_ledger}
check(f"[E, neu] MacDowell-Mansouri bindet Lambda G an das Spurgewicht: Aus den Spurwerten folgt c0 : c1 : c2 = "
      f"{c0n} : {c1n} : {c2n} (Streuung {streu:.0e}); mit kappa c1 lambda = 1/(16 pi G) und Lambda = -3 lambda ist "
      f"Lambda G = {LG} und beim nativen lambda = 4/R_Vak^2 G/R_Vak^2 = {GR}, unabhängig von jeder Länge. Native Gewichte "
      f"geben [Lambda G, G/R_Vak^2] = { {k: [f'{x:.4e}' for x in v] for k, v in tab_LG.items()} }: ein planckscher AdS-Raum. "
      f"Das Ziel G/R_Vak^2 = e^(-2/alpha)/(128 pi^3) (z98, v60 gelesen) verlangt kappa = {k_noetig} = 10^{log_k:.2f}. Die "
      f"z118-Normierung (2 c_3/4) v^2 = 1/(16 pi G) erzwingt mit kappa = 2 c_3 R_Vak = {R_z118[0]} (als Gesamtkoeffizient "
      f"gelesen {R_z118_alt[0]}) und Lambda G = {LG_z118}, gegen {LG_ledger:.3e} (dS) im Ledger",
      (c0n, c1n, c2n) == (16, 64, 384) and streu < 1e-8 and sp.simplify(LG + 3/(1024*sp.pi*kap)) == 0
      and sp.simplify(GR - 1/(4096*sp.pi*kap)) == 0 and all(abs(v[0]) > 1e-4 for v in tab_LG.values())
      and 118 < log_k < 119 and len(R_z118) == 1 and sp.simplify(R_z118[0] - 32/v_s) == 0 and sp.simplify(LG_z118 + sp.Rational(3, 256)) == 0
      and 1e-123 < LG_ledger < 1e-121)
befund("E/negativ: Das Gewicht 2 c_3 ist nicht nativ, und G = v_geo^-2 folgt nicht ohne gesetztes Gewicht. Die relative "
       "Spur und jedes endliche modulare Gewicht geben 64 x Einheit x Index; KMS-, Ward- und Indexgewichte liegen in "
       "[1 - q^2, 1]. Nativ ist 2 c_3 nur als Dichte: 2 c_3 oint K = 1 (Gauß-Bonnet der zweiblättrigen Naht, 2 = |Z_2| = "
       "Zahl der J-Summanden), als KMS-Temperatur 1/beta und als 2D-Wärmekernkoeffizient; die globale Spur integriert diese "
       "Dichte zur 1. Schärfer: In MacDowell-Mansouri bestimmt das Spurgewicht kappa allein Lambda G = -3/(1024 pi kappa); "
       "z118s 'G = v^-2 mit 2 c_3' und z118s 'Lambda = -12/R_Vak^2' zusammen ergeben R_Vak = 32 v_geo^-1 und "
       "Lambda G = -3/256, einen planckschen AdS-Raum. Für das Ledger-Lambda bräuchte die Spur das Gewicht "
       "(pi^2/32) e^(2/alpha) = 10^118,5", True)
laufzeit()

# =================================================================== Teil 4
say("\nTeil 4: Bilanz")
ol = (REPO/"research/originlab/reconstruction_search.py").read_text(encoding="utf-8")
klausel = '("graviton", st["gravity"] == "supported")' in ol
grav = re.findall(r'\("gravity", "([^"]+)", "(\w+)", "(z\d+_ergebnis\.json)"', ol)
mm = [(st_, f_) for v_, st_, f_ in grav if v_.startswith("MacDowell-Mansouri")]
stati = sorted({st_ for _, st_, _ in grav})
ZAHLEN["originlab"] = {"Klausel graviton verlangt 'supported'": klausel, "Gravitationswerte (Wert, Status, Quelle)": grav, "MM-Einträge": mm}
check(f"[E] OriginLab (nur gelesen): Die Klausel 'graviton' verlangt den Status 'supported' ({klausel}); die Gravitationsachse "
      f"hat {len(grav)} Belege mit den Status {stati}, MacDowell-Mansouri ist {mm}. Nach z127 bleibt der Wert bedingt: "
      f"Abgeleitet sind die Gewichtswahl (modular) und 'nur innen' (Indexlesart), offen bleiben P_mod/P_idx, das Paarungsbit "
      f"und die Normierung kappa; die letzte steht mit dem Ledger-Lambda im Widerspruch von rund 119 Größenordnungen",
      klausel and len(grav) >= 8 and "supported" not in stati and sorted(mm) == [("conditional", "z118_ergebnis.json"),
                                                                               ("conditional", "z98_ergebnis.json")])
bilanz = [
    "[E, neu] Gewichtswahl modular: endliche Modulelemente auf S^+ sind genau +-R (Original, Uhrgruppe) und +-i Gamma_2 "
    "(Vakuumlichtkegel); keine Holst-Mischung; P_vol in Spurform ist damit 'Gewicht = Modulelement des Originals'",
    "[E, neu] Nahtnullmode = KMS-Leiter bei beta = 4 pi auf den Marken 4j, Delta^(-i/2) dort die Eins (genau für beta in 4 pi Z)",
    "[E/C, neu] Indexlesart: äquivariante Indexspur der gemeinsamen Moduluhr = P_R-Wert 64 eps; Wärmespur bei endlichem s "
    "s-abhängig (z118); P_R schrumpft auf P_idx (Wahl des Spurfunktionals); getragen von der Z_4-Homogenität von A_q",
    "[E] Wick-Vervollständigung: X = e^(i pi Gamma_2/4) macht die gemeinsame Uhr topologisch, Wert 64 e^(-i pi/4) eps; "
    "Phasennormierung gibt 64 eps; X ist Vakuumsektor, nicht Original",
    "[E, negativ] KMS-Richtung: fixiert Uhrrichtung und Absenkschritt; Kompensator = Paarung x Uhrrichtung; "
    "Paarung = KMS-Richtung x Kappenorientierung; iota = e^(pi E) erhält alle KMS-Daten und kehrt die Paarung um",
    "[E] Paarungsbit nativ = tau-Parität des Schnitts (= Spitze der Originalebene = Orientierung der Dilatationsebene = "
    "Helizität der R = +i-Nullmodenmaterie bei festem S^+); z117 und z118 setzen es über die Cliffordkonvention "
    "entgegengesetzt",
    "[E, negativ] 2 c_3: kein natives Spurgewicht (modular |.| = 1, KMS/Ward in [1 - q^2, 1]); nativ nur als Dichte "
    "(2 c_3 oint K = 1, |Z_2| = J-Summanden), Temperatur 1/beta und Wärmekern",
    "[E, neu] MM-Normierung: Lambda G = -3/(1024 pi kappa), G/R_Vak^2 = 1/(4096 pi kappa); mit 2 c_3 und G = v^-2 folgt "
    "R_Vak = 32 v_geo^-1, Lambda G = -3/256; das Ledger-Lambda braucht kappa = 10^118,5",
    "[O] verbleibend für den MM-Wert: P_mod/P_idx (Lesart des Spurfunktionals), Paarungsbit, Normierung kappa (Hierarchie), "
    "dS-Vorzeichen (P_cap), dynamisches Vierbein; 'supported' ist ausgeschlossen, K13 bleibt L",
]
ZAHLEN["bilanz"] = bilanz
for z in bilanz:
    say("   " + z)
urteil = ("Z127 beantwortet die drei Prämissenfragen zu z118; nichts ist befördert. (1) P_R ist zum größten Teil eine Folge: "
          "Die endlichen modularen Gewichte auf S^+ sind genau +-R (Uhrgruppe <Delta^(-i/2)> des Originals, P_clock) und "
          "+-i Gamma_2 (Dilatation des Vakuumlichtkegels); eine Holst-Mischung ist in keiner modularen Familie. Die Wahl R statt "
          "Gamma_2 folgt also aus 'Gewicht = Modulelement des Originals' plus Endlichkeit, die P_vol-Paarwahl von z118 ist "
          "damit modular entschieden. 'Nur innen' folgt in der Indexlesart: Die Nahtnullmode ist bei beta = 4 pi (A*) die "
          "KMS-Leiter auf den Marken 4j, wo Delta^(-i/2) = R die Eins ist; die äquivariante Indexspur der gemeinsamen Moduluhr "
          "ist exakt der P_R-Wert 64 eps, getragen von der Z_4-Homogenität des Nahtschritts. In der Wärmespur bei endlichem "
          "s bleibt z118s Befund (s-abhängig, Ordnung e^-4pi); topologisch wird sie dort nur mit dem Wick-Element "
          "X = e^(i pi Gamma_2/4) aus dem Vakuumsektor und einer Phasennormierung. P_R schrumpft auf die Lesart des "
          "Spurfunktionals: P_mod (Gewicht aus dem Modulfluss des Originals, die Prüfhypothese selbst) und P_idx (Wirkung = "
          "äquivariante Indexspur). (2) Die KMS-Richtung fixiert das Paarungsbit nicht. Sie legt die Uhrrichtung (R, nicht R^-1) und "
          "den indextragenden Absenkschritt fest; die Kompensatorrichtung ist Paarung x Uhrrichtung, der Kompensator die "
          "Wick-Abbildung x -> -i x, die die Originalebene auf eine Kappenspitze legt, und die Vorzeichenstruktur der radialen "
          "sl(2,R) gibt Paarung = KMS-Richtung x Kappenorientierung. Die native Halbdrehung iota = e^(pi E) erhält beta > 0, R, "
          "S^+, J, Modul-U(1), Marken und A_q und kehrt Paarung, Wick-Richtung, Spitze, Lorentzzeit, gamma_5 und das Vorzeichen "
          "des eps-Terms gemeinsam um. Nativ ausdrückbar ist das Bit als tau-Parität des Schnitts; z117 (tau = +1) und z118 "
          "(tau = -1) haben es über ihre Cliffordkonventionen stillschweigend entgegengesetzt gesetzt. (3) 2 c_3 ist kein natives "
          "Spurgewicht: Modulare Gewichte haben Betrag 1, KMS-, Ward- und Indexgewichte liegen in [1 - q^2, 1]. Nativ ist 2 c_3 "
          "als Dichte (2 c_3 oint K = 1 auf der zweiblättrigen Naht, 2 = |Z_2| = Zahl der J-Summanden), als KMS-Temperatur "
          "1/beta und als 2D-Wärmekernkoeffizient; die globale Spur sieht davon nur das Integral 1. G = v_geo^-2 folgt daher "
          "nicht ohne gesetztes Gewicht. Neu und schärfer: In MacDowell-Mansouri bestimmt das Spurgewicht allein "
          "Lambda G = -3/(1024 pi kappa); z118s G = v^-2 (mit 2 c_3) und Lambda = -12/R_Vak^2 ergeben zusammen R_Vak = 32 "
          "v_geo^-1 und Lambda G = -3/256, einen planckschen AdS-Raum, während das Ledger-Lambda kappa = (pi^2/32) e^(2/alpha) "
          "= 10^118,5 verlangt. (4) Bilanz: Für den Gravitationswert MacDowell-Mansouri bleiben P_mod/P_idx, das Paarungsbit, die "
          "Normierung kappa (mit der Hierarchie e^(2/alpha)), das dS-Vorzeichen (P_cap) und ein dynamisches Vierbein. Er ruht "
          "nicht nur auf A* und P_clock und kann in der Szenario-Suche nicht 'supported' werden; die Klausel 'graviton' bleibt "
          "leer, der Wert bleibt bedingt, mit einem gescopten Widerspruch zum Ledger-Lambda. K13 bleibt L.")
say("\nUrteil\n   " + urteil)
say(f"\nLaufzeit {time.time() - T0:.0f} s")
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z127_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
