"""Z150: Gravitation und Eichung aus derselben Schleifenantwort? Ein Link in Spin(3,2) x G_SM im Transferstrom.

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart, N numerisch, O offen.

Frage (neue Strategie): z135 gibt mit einem Eichlink im Transferstrom die Plakettenstruktur und das native Verhältnis
1/g_i^2 = lambda k_i (sin^2 theta_W = 3/8), aber kein natives lambda. z118/z127 geben mit der Uhr R als Spurgewicht die
MacDowell-Mansouri-Wirkung mit l^2/G = 1024 pi kappa, aber kein natives kappa. Setzt man beide in denselben Mechanismus
(ein Link in Spin(3,2) x G_SM, dieselbe Spur über denselben Träger), ist dann (l^2/G) : (1/g_U^2) eine native reine Zahl?

Bestand: z118 (256 = Cl(4,2) x Cl(10), S^+, Nullpaar P_+-, Gamma_2, R = -i Gamma_2 gamma_5 auf S^+, 64 eps), z127 (endliche
modulare Gewichte genau +-R und +-i Gamma_2; l^2/G = 1024 pi kappa), z135 (Rückkehrblock cos^4 + sin^4 U_p, drei Lesarten
des Absolutwerts, k_i = 4), z136 (Überbestimmungstabelle für G/R_Vak^2), z137, z112, z96/z98 (Ledger-Ziel in MM-Variablen).

Teil 1  Modell und native Daten
  1.1  Cl(4,2) x Cl(10) wie z118; zwei native Träger der Dimension 64 in S^+: T16 = Twistor x 16 (Dirac-Spinor von
       so(3,2), die Uhr ist dort skalar) und der Indexkern K = S^+ cap P_- (z118; dort ist die Uhr die Chiralität).
  1.2  so(3,2) in S^+ (fünfte Richtung zeitartig), Stabilisator von K, Kontrolle so(4,1).
  1.3  Eichgeneratoren von SU(3) x SU(2) x U(1)_Y auf der 16; k_i und 3/8 als Kontrolle.
  1.4  Link g = g_grav x g_eich im 2D-Netz der Vierteldrehungen: Rückkehrblock wörtlich; das Netz ist pseudo-unitär.
  1.5  4D-Anordnung mit einem dualen Ebenenpaar: Produkt zweier Rückkehrblöcke, Zwei-Plaketten-Term.
Teil 2  Quadratische Invarianten je Träger und Einschub (Gram-Tensoren), Mischterme, höhere Ordnungen.
Teil 3  Kontinuumsform: Hodge-Struktur (dieselbe Plakette) gegen Keil-Struktur (duale Plaketten); Normierung.
Teil 4  Das Verhältnis (l^2/G) : (1/g_U^2) je Einschub und je Lesart aus z135; Zahlenwert mit alpha;
        Überbestimmungstabelle.
Teil 5  Signatur und Vorzeichen.
Teil 6  Bilanz und Urteil.
"""
import itertools
import json
import os
import time
from fractions import Fraction as Fr
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # deterministische Rundung (einfädig)

import numpy as np  # noqa: E402
import sympy as sp  # noqa: E402
from scipy.linalg import expm  # noqa: E402

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
rng = np.random.default_rng(150)
ALPHA_INV = 137.03599921684                            # gelöster Block, nur als Zahl eingesetzt
TOL = 1e-9


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val, text):
    BEFUNDE[name] = val
    say(f"   [Befund: {val}] {text}")


def e0(x):
    return f"{float(abs(x)):.0e}"


def rz(z, n=6):
    """Komplexe Zahl als gerundeter Text ohne -0."""
    z = complex(z)
    re, im = round(z.real, n) + 0.0, round(z.imag, n) + 0.0
    if im == 0:
        return f"{re:g}"
    if re == 0:
        return f"{im:g}i"
    return f"{re:g}{im:+g}i"


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


def prop(A, B):
    """A = lam B? (lam, Residuum)"""
    a, b = np.ravel(A), np.ravel(B)
    lam = np.vdot(b, a)/np.vdot(b, b)
    return complex(lam), float(np.abs(a - lam*b).max())


def tr2(A, X, Y):
    """Tr(A X Y)."""
    return complex(np.sum((A @ X)*Y.T))


# ================================================================================================== Teil 1
say("Teil 1: Modell und native Daten")
say("  1.1 Cl(4,2) x Cl(10) wie z118; die zwei nativen Träger der Dimension 64")
Gi = euklid(5)
G11 = chiral(Gi)
I8, I32, I256 = np.eye(8), np.eye(32), np.eye(256)
E6 = euklid(3)
hh = [E6[0], E6[1], E6[2], E6[3], 1j*E6[4], 1j*E6[5]]                 # Quadrate (+, +, +, +, -, -)
h7 = chiral(hh)
G16 = [kr(x, I32) for x in hh] + [1j*kr(h7, gi) for gi in Gi]
Gt16 = chiral(G16)
Sp16 = (I256 + Gt16)/2
gp16, gm16 = (G16[3] + G16[5])/2, (G16[3] - G16[5])/2
Pplus, Pminus = gm16 @ gp16, gp16 @ gm16
Gam2 = Pplus - Pminus
LI = (0, 1, 2, 4)                                                    # Lorentzrichtungen quer zum Nullpaar
ETA4 = np.diag([1., 1., 1., -1.])
g5 = chiral([G16[a] for a in LI])
Ruhr = kr(I8, 1j*G11)
lam_tot, res_tot = prop(Gt16, kr(h7, G11))
lamR_g, resR_g = prop(Ruhr @ Sp16, Gam2 @ g5 @ Sp16)
lamR_h, resR_h = prop(Ruhr @ Sp16, kr(h7, I32) @ Sp16)
diag_ok = all(np.allclose(M, np.diag(np.diag(M))) for M in (h7, G11, Gt16))
idxT = np.where(np.diag(h7).real > 0)[0]
idx16 = np.where(np.diag(G11).real > 0)[0]
P_T16 = kr(np.diag((np.diag(h7).real > 0).astype(float)), np.diag((np.diag(G11).real > 0).astype(float))).astype(complex)
P_K = Sp16 @ Pminus
TRAEGER = {"T16": P_T16, "K": P_K}
rang = {k: int(round(np.trace(P).real)) for k, P in TRAEGER.items()}
rang["S+"] = int(round(np.trace(Sp16).real))
rang["T16 cap K"] = int(round(np.trace(P_T16 @ P_K).real))
proj_ok = all(np.allclose(P @ P, P) and np.allclose(P, P.conj().T) and np.allclose(Sp16 @ P, P) for P in TRAEGER.values())
# Uhr und Gamma_2 auf den Trägern
lam_RK, res_RK = prop(Ruhr @ P_K, g5 @ P_K)
lam_GK, res_GK = prop(Gam2 @ P_K, P_K)
lam_RT, res_RT = prop(Ruhr @ P_T16, P_T16)
lam_GT, res_GT = prop(Gam2 @ P_T16, g5 @ P_T16)
# z118-Kontrolle: 64 eps und -64 (eta eta - eta eta)
EPS4 = np.zeros((4, 4, 4, 4))
for p_ in itertools.permutations(range(4)):
    EPS4[p_] = round(np.linalg.det(np.eye(4)[list(p_)]))
z118_R, z118_G = [], []
for a, b, c, d in itertools.product(range(4), repeat=4):
    if a >= b or c >= d:
        continue
    X = G16[LI[a]] @ G16[LI[b]] @ G16[LI[c]] @ G16[LI[d]]
    tR, tG = complex(np.sum((Ruhr @ P_K).T*X)), complex(np.sum((Gam2 @ P_K).T*X))
    z118_R.append(abs(tR - 64*EPS4[a, b, c, d]))
    met = ETA4[a, d]*ETA4[b, c] - ETA4[a, c]*ETA4[b, d]
    z118_G.append(abs(tG + 64*met))
ZAHLEN["traeger"] = {"Rang": rang, "Gamma_ges = lam h7 x Gamma_11": rz(lam_tot), "R|S+ = lam Gamma_2 gamma_5": rz(lamR_g),
                     "R|S+ = lam (h7 x 1)": rz(lamR_h), "K: R = lam gamma_5": rz(lam_RK), "K: Gamma_2 = lam 1": rz(lam_GK),
                     "T16: R = lam 1": rz(lam_RT), "T16: Gamma_2 = lam gamma_5": rz(lam_GT),
                     "z118-Kontrolle: max |Tr_K(R gggg) - 64 eps|": float(max(z118_R)),
                     "z118-Kontrolle: max |Tr_K(Gamma_2 gggg) + 64 (eta eta - eta eta)|": float(max(z118_G))}
check(f"[E] Native Clifford-Daten wie z118 (256 = Cl(4,2) x Cl(10), S^+ Rang {rang['S+']}): Gamma_ges = {rz(lam_tot)} h7 x "
      f"Gamma_11, auf S^+ R = {rz(lamR_g)} Gamma_2 gamma_5 = {rz(lamR_h)} (h7 x 1). In S^+ liegen zwei native Träger der "
      f"Dimension 64: T16 = (h7 = +1) x (Gamma_11 = +1), der Twistorraum mal die 16, also der Dirac-Spinor von so(3,2) mal "
      f"die 16; und der Indexkern K = S^+ cap P_- aus z118 (Rang {rang['K']}, Schnitt mit T16 Rang {rang['T16 cap K']}). "
      f"Auf T16 ist die Uhr skalar, R = {rz(lam_RT)}, und Gamma_2 = {rz(lam_GT)} gamma_5 ist die Lorentz-Chiralität; auf K "
      f"ist es umgekehrt, Gamma_2 = {rz(lam_GK)} und R = {rz(lam_RK)} gamma_5. Kontrolle z118: Tr_K(R g_a g_b g_c g_d) = 64 "
      f"eps_abcd ({e0(max(z118_R))}), Tr_K(Gamma_2 ...) = -64 (eta_ad eta_bc - eta_ac eta_bd) ({e0(max(z118_G))})",
      diag_ok and proj_ok and rang == {"T16": 64, "K": 64, "S+": 128, "T16 cap K": 32} and res_tot < TOL
      and abs(lam_tot - 1) < TOL and resR_g < TOL and abs(lamR_g + 1j) < TOL and resR_h < TOL and abs(lamR_h - 1j) < TOL
      and res_RK < TOL and abs(lam_RK - 1j) < TOL and res_GK < TOL and abs(lam_GK + 1) < TOL and res_RT < TOL
      and abs(lam_RT - 1j) < TOL and res_GT < TOL and abs(lam_GT + 1) < TOL and max(z118_R) < TOL and max(z118_G) < TOL)

say("  1.2 so(3,2) in S^+")


def so_alg(idx5):
    eta = np.array([np.round((G16[a] @ G16[a])[0, 0].real) for a in idx5])
    S = [[None]*5 for _ in range(5)]
    for i in range(5):
        for j in range(5):
            S[i][j] = (G16[idx5[i]] @ G16[idx5[j]] - G16[idx5[j]] @ G16[idx5[i]])/4
    return eta, S


def strukturrest(eta, S):
    r = 0.0
    for i, j in itertools.combinations(range(5), 2):
        for k, l in itertools.combinations(range(5), 2):
            soll = ((i == l)*0 + eta[j]*(j == k)*S[i][l] - eta[i]*(i == k)*S[j][l] - eta[j]*(j == l)*S[i][k]
                    + eta[i]*(i == l)*S[j][k])
            r = max(r, float(np.abs(S[i][j] @ S[k][l] - S[k][l] @ S[i][j] - soll).max()))
    return r


PAIRS5 = list(itertools.combinations(range(5), 2))                    # Positionen in idx5
LOR = [p for p in PAIRS5 if p[1] < 4]
TRA = [p for p in PAIRS5 if p[1] == 4]
IDX_ADS = (0, 1, 2, 4, 5)                                             # fünfte Richtung G16[5], Quadrat -1
IDX_DS = (0, 1, 2, 4, 3)                                              # Kontrolle: fünfte Richtung G16[3], Quadrat +1
ETA5, S5 = so_alg(IDX_ADS)
ETA5d, S5d = so_alg(IDX_DS)
SIG = [S5[i][j] for i, j in PAIRS5]
rest_alg, rest_alg_d = strukturrest(ETA5, S5), strukturrest(ETA5d, S5d)
pp, pp_d = prop(S5[0][4] @ S5[1][4] - S5[1][4] @ S5[0][4], S5[0][1]), prop(S5d[0][4] @ S5d[1][4] - S5d[1][4] @ S5d[0][4], S5d[0][1])
kom_T16 = max(float(np.abs(X @ P_T16 - P_T16 @ X).max()) for X in SIG)
kom_R = max(float(np.abs(X @ Ruhr - Ruhr @ X).max()) for X in SIG)
kom_G2 = {p: float(np.abs(S5[p[0]][p[1]] @ Gam2 - Gam2 @ S5[p[0]][p[1]]).max()) for p in PAIRS5}
stab_K = [p for p in PAIRS5 if kom_G2[p] < TOL]
# Stabilisator von K in so(3,2): Rang der Abbildung X -> [X, Gamma_2]
kom_mat = np.array([(X @ Gam2 - Gam2 @ X).ravel() for X in SIG])
rang_kom = int(np.linalg.matrix_rank(kom_mat, tol=1e-9))
ZAHLEN["so32"] = {"eta (0,1,2,4,5)": [float(x) for x in ETA5], "Strukturrest so(3,2)": rest_alg,
                  "[Pi_a, Pi_b] = lam Sigma_ab": rz(pp[0]), "Kontrolle so(4,1): eta": [float(x) for x in ETA5d],
                  "Kontrolle so(4,1): [Pi_a, Pi_b] = lam Sigma_ab": rz(pp_d[0]), "max ||[Sigma, P_T16]||": kom_T16,
                  "max ||[Sigma, R]||": kom_R, "Erzeuger mit [Sigma, Gamma_2] = 0": [str(p) for p in stab_K],
                  "Rang von X -> [X, Gamma_2] auf so(3,2)": rang_kom}
check(f"[E] so(3,2) in S^+: Sigma_AB = [G_A, G_B]/4 über die Richtungen (0, 1, 2, 4, 5) mit eta = "
      f"{[int(x) for x in ETA5]} erfüllt die Relationen von so(3,2) (Rest {e0(rest_alg)}); die Translationen Pi_a = "
      f"Sigma_a5 geben [Pi_a, Pi_b] = {rz(pp[0])} Sigma_ab, also F^ab = R^ab + e^a e^b/l^2 (sigma = +1, AdS, wie z98/z118). "
      f"Alle zehn Erzeuger vertauschen mit der Uhr ({e0(kom_R)}) und erhalten T16 ({e0(kom_T16)}); mit Gamma_2 vertauschen "
      f"genau die sechs Lorentz-Erzeuger (Rang der Abbildung X -> [X, Gamma_2]: {rang_kom}). Auf dem Indexkern K wirkt "
      f"also nur die Lorentzgruppe, die volle so(3,2) wirkt nur auf T16. Kontrolle: Mit der fünften Richtung G_3 "
      f"(Quadrat +1) entsteht so(4,1) mit [Pi_a, Pi_b] = {rz(pp_d[0])} Sigma_ab (dS)",
      rest_alg < TOL and rest_alg_d < TOL and pp[1] < TOL and abs(pp[0] - 1) < TOL and pp_d[1] < TOL
      and abs(pp_d[0] + 1) < TOL and kom_T16 < TOL and kom_R < TOL and stab_K == LOR and rang_kom == 4)

say("  1.3 Eichgeneratoren auf der 16")
A5 = [(Gi[2*j] + 1j*Gi[2*j + 1])/2 for j in range(5)]
YQ = sorted([Fr(0)] + [Fr(1)] + [Fr(-2, 3)]*3 + [Fr(1, 6)]*6 + [Fr(1, 3)]*3 + [Fr(-1, 2)]*2)
Y5 = np.diag([-1/3]*3 + [1/2]*2).astype(complex)


def rep32(X, ops):
    return sum(ops[i].conj().T @ (X[i, j]*ops[j]) for i in range(5) for j in range(5) if X[i, j] != 0) - np.trace(X)/2*I32


OPS = None
for kand in (A5, [x.conj().T for x in A5]):
    spek = sorted(Fr(float(v)).limit_denominator(12) for v in np.diag(rep32(Y5, kand)).real[idx16])
    if spek == YQ:
        OPS = kand
        break
car = max(float(np.abs(OPS[i] @ OPS[j].conj().T + OPS[j].conj().T @ OPS[i] - (i == j)*I32).max())
          for i in range(5) for j in range(5))
GM = np.zeros((8, 3, 3), complex)
GM[0][0, 1] = GM[0][1, 0] = 1
GM[1][0, 1], GM[1][1, 0] = -1j, 1j
GM[2][0, 0], GM[2][1, 1] = 1, -1
GM[3][0, 2] = GM[3][2, 0] = 1
GM[4][0, 2], GM[4][2, 0] = -1j, 1j
GM[5][1, 2] = GM[5][2, 1] = 1
GM[6][1, 2], GM[6][2, 1] = -1j, 1j
GM[7] = np.diag([1, 1, -2])/np.sqrt(3)
PAULI = [sx, sy, sz]
T32, TNAME = [], []
for a in range(8):
    X_ = np.zeros((5, 5), complex)
    X_[:3, :3] = GM[a]/2
    T32.append(rep32(X_, OPS))
    TNAME.append(f"SU(3) {a + 1}")
for a in range(3):
    X_ = np.zeros((5, 5), complex)
    X_[3:, 3:] = PAULI[a]/2
    T32.append(rep32(X_, OPS))
    TNAME.append(f"SU(2) {a + 1}")
Y32 = rep32(Y5, OPS)
T32.append(np.sqrt(3/5)*Y32)
TNAME.append("U(1)_Y (GUT)")
NG = len(T32)
TM = [X[np.ix_(idx16, idx16)] for X in T32]                           # auf der 16
TQ = [X[np.ix_(np.setdiff1d(np.arange(32), idx16), np.setdiff1d(np.arange(32), idx16))] for X in T32]   # auf der 16-quer
T256 = [kr(I8, X) for X in T32]
gram16 = np.array([[np.trace(A_ @ B_) for B_ in TM] for A_ in TM])
trY2 = float(np.trace(Y32[np.ix_(idx16, idx16)] @ Y32[np.ix_(idx16, idx16)]).real)
trY3 = float(abs(np.trace(np.linalg.matrix_power(Y32[np.ix_(idx16, idx16)], 3))))
k_i = [Fr(float(2*gram16[i, i].real)).limit_denominator(100) for i in (0, 8, 11)]
sin2 = Fr(2)/(Fr(2) + Fr(trY2).limit_denominator(100))
kom_ge = max(float(np.abs(X @ T_ - T_ @ X).max()) for X in SIG for T_ in T256)
kom_g11 = max(float(np.abs(X @ G11 - G11 @ X).max()) for X in T32)
ZAHLEN["eichung"] = {"CAR": car, "Gram Tr_16(T_i T_j): Abstand zu 2 delta": float(np.abs(gram16 - 2*np.eye(NG)).max()),
                     "Tr_16 Y^2": trY2, "|Tr_16 Y^3|": trY3, "k_i (SU(3), SU(2), U(1)_Y)": [str(k) for k in k_i],
                     "sin^2 theta_W": str(sin2), "max ||[Sigma_AB, T_i]||": kom_ge}
check(f"[E] Eichgeneratoren von SU(3) x SU(2) x U(1)_Y in so(10) auf der 16 = (Gamma_11 = +1) (Hyperladungsspektrum wie "
      f"z82/z135): Tr_16(T_i T_j) = 2 delta_ij für alle zwölf (Abweichung {e0(np.abs(gram16 - 2*np.eye(NG)).max())}), "
      f"Tr_16 Y^2 = {Fr(trY2).limit_denominator(100)}, Tr_16 Y^3 = 0 ({e0(trY3)}); k_i = 2 Tr_16 T_i^2 = "
      f"{', '.join(str(k) for k in k_i)} und sin^2 theta_W = {sin2} (Kontrolle z135). Die Eichgeneratoren vertauschen mit "
      f"allen Erzeugern von so(3,2) ({e0(kom_ge)}): Der Link ist ein direktes Produkt g_grav x g_eich",
      OPS is not None and car < TOL and np.abs(gram16 - 2*np.eye(NG)).max() < TOL and abs(trY2 - 10/3) < TOL
      and trY3 < TOL and all(k == 4 for k in k_i) and sin2 == Fr(3, 8) and kom_ge < TOL and kom_g11 < TOL)

say("  1.4 Link in Spin(3,2) x G_SM im 2D-Netz der Vierteldrehungen (z135 wörtlich)")
S8 = [[(hh[IDX_ADS[i]] @ hh[IDX_ADS[j]] - hh[IDX_ADS[j]] @ hh[IDX_ADS[i]])/4 for j in range(5)] for i in range(5)]
ST = [S8[i][j][np.ix_(idxT, idxT)] for i, j in PAIRS5]               # so(3,2) auf dem Twistorraum (4 x 4)
G8 = 1j*hh[4] @ hh[5]                                                # invariante hermitesche Form
GT = G8[np.ix_(idxT, idxT)]
form_rest = max(float(np.abs(X.conj().T @ GT + GT @ X).max()) for X in ST)
form_sig = sorted(np.round(np.linalg.eigvalsh(GT), 6).tolist())
KOMPAKT = [n for n, (i, j) in enumerate(PAIRS5) if ETA5[i]*ETA5[j] > 0]
herm_rest = {n: float(np.abs(ST[n] + ST[n].conj().T).max()) for n in range(10)}
TH = np.pi/4
CC, SS = np.cos(TH), np.sin(TH)


def zufall_grav(skala, nur_kompakt=False):
    om = rng.normal(size=10)*skala
    if nur_kompakt:
        om = np.array([om[n] if n in KOMPAKT else 0.0 for n in range(10)])
    return sum(o*X for o, X in zip(om, ST))


def zufall_eich(skala):
    th = rng.normal(size=NG)*skala
    return 1j*sum(t*X for t, X in zip(th, TM))


def floquet(nsite, d, layers, links, linv, th):
    """Einteilchen-Floquetoperator des Netzes; Hinweg U (von B nach A), Rückweg U^-1 (zeilenweise je Bindung)."""
    c, s = np.cos(th), np.sin(th)
    F = np.eye(nsite*d, dtype=complex)
    for lay in layers:
        for ia, ib, key in lay:
            sa, sb = slice(ia*d, ia*d + d), slice(ib*d, ib*d + d)
            Fa, Fb = F[sa].copy(), F[sb].copy()
            F[sa] = c*Fa - 1j*s*(links[key] @ Fb)
            F[sb] = c*Fb - 1j*s*(linv[key] @ Fa)
    return F


def pfad(x, layers, links, linv, d):
    pos, M = x, np.eye(d, dtype=complex)
    for lay in layers:
        for ia, ib, key in lay:
            if pos == ia:
                M, pos = linv[key] @ M, ib
                break
            if pos == ib:
                M, pos = links[key] @ M, ia
                break
    return pos, M


BS2 = [(1, 0), (0, 1), (-1, 0), (0, -1)]
L2 = 4
SITES2 = [(x, y) for x in range(L2) for y in range(L2)]
IDX2 = {s_: i for i, s_ in enumerate(SITES2)}
A2 = [s_ for s_ in SITES2 if (s_[0] + s_[1]) % 2 == 0]
LAY2 = [[(IDX2[a], IDX2[((a[0] + b[0]) % L2, (a[1] + b[1]) % L2)], (j, a)) for a in A2] for j, b in enumerate(BS2)]
KEYS2 = [k for lay in LAY2 for _, _, k in lay]
NS2 = len(SITES2)
netz2 = {}
for name, kompakt in (("Spin(3,2) x G_SM", False), ("U(2) x G_SM (maximal kompakt)", True)):
    links, linv = {}, {}
    for k in KEYS2:
        gS, gM = expm(zufall_grav(0.35, kompakt)), expm(zufall_eich(0.5))
        links[k], linv[k] = np.kron(gS, gM), np.kron(np.linalg.inv(gS), gM.conj().T)
    F = floquet(NS2, 64, LAY2, links, linv, TH)
    abw, spur_abw, spur_hol = 0.0, 0.0, 0j
    for x in range(NS2):
        end, M = pfad(x, LAY2, links, linv, 64)
        blk = F[x*64:(x + 1)*64, x*64:(x + 1)*64]
        abw = max(abw, float(np.abs(blk - (CC**4*np.eye(64) + SS**4*M)).max()), 0.0 if end == x else 1.0)
        spur_hol += np.trace(M)
    Gnetz = np.kron(np.eye(NS2), np.kron(GT, np.eye(16)))
    netz2[name] = {"Rückkehrblock": abw, "Spur": float(abs(np.trace(F) - (NS2*64*CC**4 + SS**4*spur_hol))),
                   "||F^+ F - 1||": float(np.abs(F.conj().T @ F - np.eye(NS2*64)).max()),
                   "||F^+ G F - G||": float(np.abs(F.conj().T @ Gnetz @ F - Gnetz).max()),
                   "Im Tr F": float(abs(np.trace(F).imag))}
gS, gM = expm(zufall_grav(0.35)), expm(zufall_eich(0.5))
prod_abw = float(abs(np.trace(np.kron(gS, gM)) - np.trace(gS)*np.trace(gM)))
reell_S = float(abs(np.trace(gS).imag))
ZAHLEN["netz_2d"] = {"Fälle": netz2, "Form G auf T: Rest": form_rest, "Signatur von G": form_sig,
                     "kompakte Erzeuger (Positionen)": [str(PAIRS5[n]) for n in KOMPAKT],
                     "Tr(g_S x g_M) - Tr g_S Tr g_M": prod_abw, "Im Tr_T g_S": reell_S}
n_a, n_b = netz2["Spin(3,2) x G_SM"], netz2["U(2) x G_SM (maximal kompakt)"]
check(f"[E] Link g = g_grav x g_eich in Spin(3,2) x G_SM auf T16 (d = 64) im 2D-Netz der Vierteldrehungen (Torus 4 x 4, "
      f"Teilschichten +x, +y, -x, -y, Rückweg mit g^-1): Der Rückkehrblock ist wörtlich cos^4 theta + sin^4 theta U_p "
      f"({e0(n_a['Rückkehrblock'])}), Tr F = N d cos^4 + sin^4 sum_p Tr U_p ({e0(n_a['Spur'])}), und Tr U_p = Tr_T(g_grav) "
      f"Tr_16(g_eich) ist ein Produkt ({e0(prod_abw)}), Tr_T(g_grav) reell ({e0(reell_S)}). Das Netz ist aber nicht "
      f"unitär (||F^+ F - 1|| = {n_a['||F^+ F - 1||']:.3f}), nur pseudo-unitär bezüglich der su(2,2)-Form G = i G_4 G_5 "
      f"mit Signatur {form_sig} (||F^+ G F - G|| = {e0(n_a['||F^+ G F - G||'])}). Unitär ist es genau für Links in der "
      f"maximal kompakten U(2) = Spin(3) x Spin(2) ({len(KOMPAKT)} Erzeuger: ||F^+ F - 1|| = {e0(n_b['||F^+ F - 1||'])})",
      n_a["Rückkehrblock"] < 1e-9 and n_a["Spur"] < 1e-8 and n_a["||F^+ F - 1||"] > 0.05 and n_a["||F^+ G F - G||"] < 1e-8
      and n_b["Rückkehrblock"] < 1e-9 and n_b["||F^+ F - 1||"] < 1e-9 and form_rest < TOL and form_sig == [-1.0, -1.0, 1.0, 1.0]
      and len(KOMPAKT) == 4 and prod_abw < 1e-9 and reell_S < 1e-9
      and all((herm_rest[n] < TOL) == (n in KOMPAKT) for n in range(10)))

say("  1.5 4D-Anordnung mit einem dualen Ebenenpaar (xy, dann zt)")
DIR8 = [(1, 0, 0, 0), (0, 1, 0, 0), (-1, 0, 0, 0), (0, -1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (0, 0, -1, 0), (0, 0, 0, -1)]
L4 = 4
SITES4 = list(itertools.product(range(L4), repeat=4))
IDX4 = {s_: i for i, s_ in enumerate(SITES4)}
A4 = [s_ for s_ in SITES4 if sum(s_) % 2 == 0]
LAY4 = [[(IDX4[a], IDX4[tuple((a[i] + b[i]) % L4 for i in range(4))], (j, a)) for a in A4] for j, b in enumerate(DIR8)]
KEYS4 = [k for lay in LAY4 for _, _, k in lay]
NS4 = len(SITES4)
links4 = {k: expm(zufall_grav(0.3)) for k in KEYS4}
linv4 = {k: np.linalg.inv(v) for k, v in links4.items()}
netz4 = {}
for th in (TH, 0.6):
    c_, s_ = np.cos(th), np.sin(th)
    F4 = floquet(NS4, 4, LAY4, links4, linv4, th)
    abw4, einzel4 = 0.0, 0.0
    for x in range(NS4):
        e1, Mxy = pfad(x, LAY4[:4], links4, linv4, 4)
        e2, Mzt = pfad(x, LAY4[4:], links4, linv4, 4)
        blk = F4[x*4:(x + 1)*4, x*4:(x + 1)*4]
        soll = (c_**4*np.eye(4) + s_**4*Mzt) @ (c_**4*np.eye(4) + s_**4*Mxy)
        abw4 = max(abw4, float(np.abs(blk - soll).max()), 0.0 if (e1 == x and e2 == x) else 1.0)
        einzel4 = max(einzel4, float(np.abs(blk - (c_**8*np.eye(4) + s_**4*c_**4*(Mxy + Mzt))).max()))
    netz4[round(th/np.pi, 4)] = {"Abweichung vom Produkt": abw4, "Abweichung ohne Zwei-Plaketten-Term": einzel4,
                                 "a_B/a_A = 2 s^4/(c^4 + s^4)": float(2*s_**4/(c_**4 + s_**4))}
# zweite Ordnung des Blocks: a_A (Tr W A^2 + Tr W B^2) + a_B Tr(W B A)
Aq, Bq = zufall_grav(1.0), zufall_grav(1.0)
G2_8 = (hh[3] - hh[5]) @ (hh[3] + hh[5])/4 - (hh[3] + hh[5]) @ (hh[3] - hh[5])/4     # Gamma_2 auf Cl(4,2)
W5S = (1j*G2_8)[np.ix_(idxT, idxT)]                                   # i Gamma_2 auf dem Twistorraum
g2_ok = float(np.abs(kr(G2_8, I32) - Gam2).max())
zweite = {}
for nameW, Wq in (("1", np.eye(4, dtype=complex)), ("i Gamma_2", W5S)):
    def f_(eps, Wq=Wq):
        return np.trace(Wq @ (CC**4*np.eye(4) + SS**4*expm(eps*Bq)) @ (CC**4*np.eye(4) + SS**4*expm(eps*Aq)))
    h_ = 1e-3
    num = (f_(h_) + f_(-h_) - 2*f_(0.0))/(2*h_**2)
    a_A, a_B = SS**4*(CC**4 + SS**4)/2, SS**8
    soll = a_A*(np.trace(Wq @ Aq @ Aq) + np.trace(Wq @ Bq @ Bq)) + a_B*np.trace(Wq @ Bq @ Aq)
    zweite[nameW] = float(abs(num - soll)/max(abs(soll), 1e-12))
ZAHLEN["netz_4d"] = {"Fälle je theta/pi": {str(k): v for k, v in netz4.items()}, "zweite Ordnung: rel. Abweichung": zweite,
                     "a_A = s^4 (c^4 + s^4)/2 bei pi/4": float(SS**4*(CC**4 + SS**4)/2), "a_B = s^8 bei pi/4": float(SS**8)}
check(f"[E] 4D-Anordnung mit einem dualen Ebenenpaar (Torus 4^4, acht Teilschichten +x, +y, -x, -y, +z, +t, -z, -t, "
      f"Spin(3,2)-Links): Der Rückkehrblock ist das Produkt (cos^4 + sin^4 U_zt)(cos^4 + sin^4 U_xy) (Abweichung "
      f"{e0(max(v['Abweichung vom Produkt'] for v in netz4.values()))} bei theta = pi/4 und 0,6; ohne den Term U_zt U_xy "
      f"{min(v['Abweichung ohne Zwei-Plaketten-Term'] for v in netz4.values()):.3f}). In zweiter Ordnung trägt er a_A (Tr W "
      f"F_xy^2 + Tr W F_zt^2) + a_B Tr(W F_zt F_xy) mit a_A = sin^4 (cos^4 + sin^4)/2 und a_B = sin^8 (rel. Abweichung "
      f"{e0(max(zweite.values()))}); bei der nativen Vierteldrehung ist a_A = a_B = 1/16. Das Netz liefert also zwei "
      f"Strukturen: dieselbe Plakette zweimal (Hodge-Typ, F_mn F^mn) und zwei duale Plaketten (Keil-Typ, F ^ F). Der "
      f"wörtliche Block aus z135 enthält nur die erste. Die 4D-Anordnung ist eine Konstruktion (Firewall), kein Datum",
      all(v["Abweichung vom Produkt"] < 1e-9 and v["Abweichung ohne Zwei-Plaketten-Term"] > 1e-3 for v in netz4.values())
      and max(zweite.values()) < 1e-5 and abs(netz4[0.25]["a_B/a_A = 2 s^4/(c^4 + s^4)"] - 1) < 1e-12 and g2_ok < TOL)
print(f"   (Laufzeit bisher {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)

# ================================================================================================== Teil 2
say("\nTeil 2: Quadratische Invarianten je Träger und Einschub")
say("  2.1 Gram-Tensoren Tr_C(W X Y)")
ETA5m = np.diag(ETA5)
K10 = np.array([[-(ETA5[i]*(i == k)*ETA5[j]*(j == l) - ETA5[i]*(i == l)*ETA5[j]*(j == k)) for (k, l) in PAIRS5]
                for (i, j) in PAIRS5], dtype=float)                  # Killing-Typ: -(eta_AC eta_BD - eta_AD eta_BC)
E10 = np.zeros((10, 10))
for m_, (i, j) in enumerate(PAIRS5):
    for n_, (k, l) in enumerate(PAIRS5):
        if j < 4 and l < 4:
            E10[m_, n_] = EPS4[i, j, k, l]
GEW = {"1": I256, "Gamma_2": Gam2, "+i Gamma_2": 1j*Gam2, "-i Gamma_2": -1j*Gam2, "+R": Ruhr, "-R": -Ruhr, "gamma_5": g5}
NATIV = ("+i Gamma_2", "-i Gamma_2", "+R", "-R")
gram = {}
for cname, P in TRAEGER.items():
    for wname, Wm in GEW.items():
        A_ = P @ Wm
        Qgg = np.array([[tr2(A_, X, Y) for Y in SIG] for X in SIG])
        Qee = np.array([[tr2(A_, X, Y) for Y in T256] for X in T256])
        Qge = np.array([[tr2(A_, X, Y) for Y in T256] for X in SIG])
        B_ = np.array([K10.ravel(), E10.ravel()]).T
        koef, *_ = np.linalg.lstsq(B_.astype(complex), Qgg.ravel(), rcond=None)
        rest = float(np.abs(B_ @ koef - Qgg.ravel()).max())
        gam_ = np.trace(Qee)/NG
        gram[(cname, wname)] = {"alpha": complex(koef[0]), "beta": complex(koef[1]), "Rest": rest, "gamma": complex(gam_),
                                "Rest Eich": float(np.abs(Qee - gam_*np.eye(NG)).max()), "misch": float(np.abs(Qge).max()),
                                "Tr W": complex(np.trace(A_))}
soll = {("T16", "1"): (16, 0, 8), ("T16", "Gamma_2"): (0, 16j, 0), ("T16", "+i Gamma_2"): (0, -16, 0),
        ("T16", "-i Gamma_2"): (0, 16, 0), ("T16", "+R"): (16j, 0, 8j), ("T16", "-R"): (-16j, 0, -8j),
        ("T16", "gamma_5"): (0, -16j, 0), ("K", "1"): (16, 0, 8), ("K", "Gamma_2"): (-16, 0, -8),
        ("K", "+i Gamma_2"): (-16j, 0, -8j), ("K", "-i Gamma_2"): (16j, 0, 8j), ("K", "+R"): (0, 16, 0),
        ("K", "-R"): (0, -16, 0), ("K", "gamma_5"): (0, -16j, 0)}
gram_ok = all(abs(v["alpha"] - soll[k][0]) < 1e-8 and abs(v["beta"] - soll[k][1]) < 1e-8 and abs(v["gamma"] - soll[k][2]) < 1e-8
              and v["Rest"] < 1e-8 and v["Rest Eich"] < 1e-8 and v["misch"] < 1e-8 for k, v in gram.items())
ZAHLEN["gram"] = {f"{c} | {w}": {"alpha (Killing)": rz(v["alpha"]), "beta (eps)": rz(v["beta"]), "gamma (Eichung)": rz(v["gamma"]),
                                  "Rest": v["Rest"], "Mischterm": v["misch"]} for (c, w), v in gram.items()}
say("     Träger | Einschub    | alpha (Killing-Form so(3,2)) | beta (eps_abcd) | gamma (Tr T_i T_j = gamma delta) | Mischterm")
for (c, w), v in gram.items():
    say(f"     {c:<6} | {w:<11} | {rz(v['alpha']):<28} | {rz(v['beta']):<15} | {rz(v['gamma']):<32} | {e0(v['misch'])}")
check("[E] Gram-Tensoren der zweiten Ordnung, Tr_C(W F^2) mit F = (1/2) F^AB Sigma_AB + i F^i T_i, für beide Träger und "
      "sieben Einschübe: Der Gravitationsteil ist alpha mal die Killing-Form von so(3,2) (-(eta_AC eta_BD - eta_AD "
      "eta_BC), auf allen zehn Erzeugern, Translationen eingeschlossen) plus beta mal eps_abcd (nur Lorentz), der Eichteil "
      "ist gamma delta_ij, der Mischteil Tr_C(W Sigma_AB T_i) verschwindet immer. Ohne Einschub: alpha = 16, gamma = 8 "
      "auf beiden Trägern (k_i = 2 gamma = 16 = 4 mal 4, Verhältnis 3/8 unverändert), beta = 0. Die vier endlichen "
      "modularen Gewichte aus z127 sind auf jedem Träger entweder ein Skalar +-i (alpha = +-16 i, gamma = +-8 i, beta = 0) "
      "oder die Chiralität +-i gamma_5 (alpha = 0, gamma = 0, beta = +-16 reell). Auf T16 ist die Uhr der Skalar und "
      "i Gamma_2 die Chiralität, auf K umgekehrt", gram_ok)

say("  2.2 Endlichkeit der relativen Spur (Kriterium aus z118) und der Eichterm im Indexfunktional")
endlich = {}
for wname in ("1", "gamma_5", "Gamma_2", "+i Gamma_2", "+R"):
    A_ = Sp16 @ GEW[wname]
    sL = max(abs(tr2(A_, SIG[m_], SIG[n_])) for m_ in range(10) for n_ in range(10) if PAIRS5[m_] in LOR and PAIRS5[n_] in LOR)
    sE = max(abs(tr2(A_, X, Y)) for X in T256 for Y in T256)
    endlich[wname] = {"Lorentz": float(sL), "Eichung": float(sE), "endlich": bool(sL < 1e-8 and sE < 1e-8)}
holst = {}
for a_, b_ in ((1.0, 0.0), (0.0, 1.0), (1.0, 0.37), (0.5, -2.0)):
    A_ = P_K @ (a_*Gam2 + b_*Ruhr)
    holst[f"{a_:g} Gamma_2 + {b_:g} R"] = (complex(tr2(A_, T256[0], T256[0])), complex(tr2(A_, SIG[0], SIG[7])))
holst_ok = all(abs(v[0] + 8*float(k.split()[0])) < 1e-8 for k, v in holst.items())
eps_ref = E10[0, 7]                                                # Positionen (0,1) und (2,3): Richtungen 0, 1 und 2, 4
ZAHLEN["endlichkeit"] = {"Tr_(S+)(W X) je Einschub": endlich,
                         "Indexkern, W = a Gamma_2 + b R: (Tr W T_1 T_1, Tr W Sigma_01 Sigma_23)":
                             {k: [rz(v[0]), rz(v[1])] for k, v in holst.items()}}
check(f"[E] Endlichkeit wie in z118 (Aufhebung zwischen P_+ und P_-, also Tr_(S^+)(W X) = 0) für Lorentz- und "
      f"Eichprodukte: endlich sind Gamma_2, i Gamma_2 und R ({[k for k, v in endlich.items() if v['endlich']]}), nicht 1 "
      f"und gamma_5 ({[k for k, v in endlich.items() if not v['endlich']]}; für 1 ist Tr_(S^+) T_i T_i = "
      f"{endlich['1']['Eichung']:.0f}). Im Indexfunktional auf K hat W = a Gamma_2 + b R den Eichterm -8 a delta_ij und "
      f"den eps-Term 16 b eps_abcd (geprüft an Sigma_01 Sigma_23, eps_0123 = {eps_ref:g}): {', '.join(f'{k}: ({rz(v[0])}, {rz(v[1])})' for k, v in holst.items())}. "
      f"Beide zugleich gibt es nur für a b != 0, die Holst-Mischung, die nach z127 kein modulares Gewicht ist. Der "
      f"einschubfreie Träger (W = 1) ist als relative Spur nicht endlich",
      [k for k, v in endlich.items() if v["endlich"]] == ["Gamma_2", "+i Gamma_2", "+R"] and holst_ok and abs(eps_ref) == 1
      and abs(holst["0 Gamma_2 + 1 R"][1] - 16*eps_ref) < 1e-8 and abs(holst["1 Gamma_2 + 0 R"][1]) < 1e-8
      and abs(holst["1 Gamma_2 + 0.37 R"][1] - 0.37*16*eps_ref) < 1e-8)
befund("kein gemeinsamer Einschub", "ja",
       "E: Kein natives endliches Gewicht gibt zugleich den Eichterm und den eps-Term. Die Gewichte vom Skalartyp (auf K "
       "+-i Gamma_2, auf T16 +-R) geben Killing-Form plus Eichterm, rein imaginär; die vom Chiralitätstyp (auf K +-R, auf "
       "T16 +-i Gamma_2) geben den eps-Term reell und den Eichterm exakt null. MacDowell-Mansouri und Yang-Mills liegen "
       "also nie im selben Spurfunktional; ihr relatives Gewicht ist genau das Holst-Verhältnis a : b.")

say("  2.3 Mischterme und höhere Ordnungen auf T16")
ordn = {}
for wname, WS in (("1", np.eye(4, dtype=complex)), ("i Gamma_2", W5S)):
    FS, FM = zufall_grav(1.0), zufall_eich(1.0)
    Fv = np.kron(FS, np.eye(16)) + np.kron(np.eye(4), FM)
    Wv = np.kron(WS, np.eye(16))
    zeile = {}
    for n in (1, 2, 3, 4):
        tot = np.trace(Wv @ np.linalg.matrix_power(Fv, n))
        rein_g = 16*np.trace(WS @ np.linalg.matrix_power(FS, n))
        rein_e = np.trace(WS)*np.trace(np.linalg.matrix_power(FM, n))
        zeile[n] = {"gesamt": complex(tot), "rein grav": complex(rein_g), "rein eich": complex(rein_e),
                    "gemischt": complex(tot - rein_g - rein_e)}
    zeile["4: 6 Tr_S(W F_g^2) Tr_16(F_e^2)"] = complex(6*np.trace(WS @ FS @ FS)*np.trace(FM @ FM))
    ordn[wname] = zeile
ZAHLEN["ordnungen"] = {w: {str(k): ({kk: rz(vv) for kk, vv in v.items()} if isinstance(v, dict) else rz(v))
                           for k, v in z.items()} for w, z in ordn.items()}
misch_ok = all(abs(z[n]["gemischt"]) < 1e-8 for z in ordn.values() for n in (1, 2, 3)) and all(
    abs(z[4]["gemischt"] - z["4: 6 Tr_S(W F_g^2) Tr_16(F_e^2)"]) < 1e-7 and abs(z[4]["gemischt"]) > 1e-3 for z in ordn.values())
kub_ok = all(abs(z[3]["gesamt"]) < 1e-8 and abs(z[1]["gesamt"]) < 1e-8 for z in ordn.values())
check(f"[E] Mischterme zwischen Gravitation und Eichung (T16, zufällige Elemente von so(3,2) + su(3) + su(2) + u(1)): in "
      f"erster bis dritter Ordnung null (Tr_16 T_i = 0 samt Tr Y = 0; Tr_T(W Sigma) = 0; so(3,2) = sp(4) und die 16 haben "
      f"keine kubische Spur: Tr_C(W F^3) = 0), für beide Einschubtypen. Der erste Mischterm steht in vierter Ordnung und "
      f"ist das Produkt 6 Tr_T(W F_g^2) Tr_16(F_e^2) (ohne Einschub {rz(ordn['1'][4]['gemischt'], 4)}, mit Chiralität "
      f"{rz(ordn['i Gamma_2'][4]['gemischt'], 4)}); er folgt aus Tr U_p = Tr_T(g_grav) Tr_16(g_eich) und ist im Kontinuum "
      f"ein Operator der Dimension acht. Am MM-Vakuum F_grav = 0 verschwindet er", misch_ok and kub_ok)

say("  2.4 Uhr auf dem Indexkern und die Eichantwort")
om10 = rng.normal(size=(10, 10))
om10 = om10 - om10.T
X32 = (1j/4)*sum(om10[a, b]*Gi[a] @ Gi[b] for a in range(10) for b in range(10))
XSM = 1j*sum(t*X for t, X in zip(rng.normal(size=NG), T32))
uhr_e = {}
for name, Xg in (("so(10) allgemein", X32), ("SM-Unteralgebra", XSM)):
    X256 = kr(I8, Xg)
    X16 = Xg[np.ix_(idx16, idx16)]
    werte = [complex(np.trace(P_K @ Ruhr @ np.linalg.matrix_power(X256, n))) for n in range(1, 6)]
    uhr_e[name] = {"Tr_K(R X^n), n = 1..5": [rz(w, 5) for w in werte], "4 i Tr_16 X^5": rz(4j*np.trace(np.linalg.matrix_power(X16, 5)), 5),
                   "ok": bool(max(abs(w) for w in werte[:4]) < 1e-7 and abs(werte[4]) > 1e-3
                              and abs(werte[4] - 4j*np.trace(np.linalg.matrix_power(X16, 5))) < 1e-6*abs(werte[4]))}
ZAHLEN["uhr_eichung"] = uhr_e
check(f"[E] Mit der Uhr auf dem Indexkern K = (2_L x 16) + (2_R x 16-quer) ist die reine Eichantwort Tr_K(R g_eich) = "
      f"2 i (Tr_16 - Tr_16-quer)(g_eich): Sie verschwindet bis zur vierten Ordnung und beginnt mit 4 i Tr_16 X^5, dem "
      f"chiralen Pfaffschen Term aus z135 (so(10) allgemein: {uhr_e['so(10) allgemein']['Tr_K(R X^n), n = 1..5'][4]}; "
      f"SM-Unteralgebra: {uhr_e['SM-Unteralgebra']['Tr_K(R X^n), n = 1..5'][4]}). Ein Eichterm F^2 kommt mit der Uhr nicht vor",
      all(v["ok"] for v in uhr_e.values()))
print(f"   (Laufzeit bisher {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)

# ================================================================================================== Teil 3
say("\nTeil 3: Kontinuumsform: Hodge (dieselbe Plakette) gegen Keil (duale Plaketten)")
say("  3.1 Vier Zellen: Struktur mal Einschubtyp")


def tens4(P, Wm):
    T = np.zeros((4, 4, 4, 4), complex)
    A_ = P @ Wm
    for a, b, c, d in itertools.product(range(4), repeat=4):
        if a == b or c == d:
            continue
        T[a, b, c, d] = tr2(A_, S5[a][b], S5[c][d])
    return T


T_PLAIN = tens4(P_T16, I256)                                         # 16 mal Killing
T_CHIR = tens4(P_K, Ruhr)                                            # 16 eps (z118: 64 eps für gamma_a gamma_b)
t_real = max(float(np.abs(T_PLAIN.imag).max()), float(np.abs(T_CHIR.imag).max()))
T_PLAIN, T_CHIR = T_PLAIN.real, T_CHIR.real


def kn(S1, S2):
    return (np.einsum("ac,bd->abcd", S1, S2) + np.einsum("bd,ac->abcd", S1, S2) - np.einsum("ad,bc->abcd", S1, S2)
            - np.einsum("bc,ad->abcd", S1, S2))


def hodge(T, F, ginv, dete):
    return abs(dete)*0.25*np.einsum("abcd,abmn,cdrs,mr,ns->", T, F, F, ginv, ginv, optimize=True)


def keil(T, F):
    return np.einsum("abcd,mnrs,abmn,cdrs->", T, EPS4, F, F, optimize=True)/16


zellen = {"A, ohne": [], "A, chiral": [], "B, ohne": [], "B, chiral": []}
zykl_max = 0.0
for _ in range(5):
    Rabcd = sum(kn(S_, S_)/2 for S_ in [(lambda M: M + M.T)(rng.normal(size=(4, 4))) for _ in range(3)])
    zykl_max = max(zykl_max, float(np.abs(Rabcd + np.einsum("acdb->abcd", Rabcd) + np.einsum("adbc->abcd", Rabcd)).max()))
    Rup = np.einsum("ai,bj,ijcd->abcd", ETA4, ETA4, Rabcd)                # R^ab_cd
    Ric = np.einsum("ac,abcd->bd", ETA4, Rabcd)
    Rsk = float(np.einsum("bd,bd->", ETA4, Ric))
    Riem2 = float(np.einsum("abcd,ijkl,ai,bj,ck,dl->", Rabcd, Rabcd, ETA4, ETA4, ETA4, ETA4))
    Eul = Riem2 - 4*float(np.einsum("ab,ab->", Ric, ETA4 @ Ric @ ETA4)) + Rsk**2
    Pont = float(np.einsum("cdef,abcd,abef->", EPS4, Rup, Rabcd))
    e = rng.normal(size=(4, 4)) + 2*np.eye(4)                              # e^a_m
    if np.linalg.det(e) < 0:
        e[:, 0] *= -1
    dete = float(np.linalg.det(e))
    ginv = np.linalg.inv(e.T @ ETA4 @ e)
    Rf = np.einsum("abcd,cm,dn->abmn", Rup, e, e)
    EE = np.einsum("am,bn->abmn", e, e) - np.einsum("an,bm->abmn", e, e)
    for zname, fun, T in (("A, ohne", "h", T_PLAIN), ("A, chiral", "h", T_CHIR), ("B, ohne", "k", T_PLAIN),
                          ("B, chiral", "k", T_CHIR)):
        if fun == "h":
            L0, Lp, Lm = (hodge(T, Rf + l_*EE, ginv, dete) for l_ in (0.0, 1.0, -1.0))
        else:
            L0, Lp, Lm = (keil(T, Rf + l_*EE) for l_ in (0.0, 1.0, -1.0))
        zellen[zname].append({"L0": L0/dete, "L1": (Lp - Lm)/2/dete, "L2": ((Lp + Lm)/2 - L0)/dete, "R": Rsk, "Riem2": Riem2,
                              "Euler": Eul, "Pont": Pont})


def konst(xs):
    return float(max(abs(x - xs[0]) for x in xs)/max(abs(xs[0]), 1e-12))


def skala(z):
    return max(abs(m_["L0"]) for m_ in z)


ao = zellen["A, ohne"]
c_ao = ([m_["L0"]/m_["Riem2"] for m_ in ao], [m_["L1"]/m_["R"] for m_ in ao], [m_["L2"] for m_ in ao])
bc = zellen["B, chiral"]
c_bc = ([m_["L0"]/m_["Euler"] for m_ in bc], [m_["L1"]/m_["R"] for m_ in bc], [m_["L2"] for m_ in bc])
ac, bo = zellen["A, chiral"], zellen["B, ohne"]
c_ac = [m_["L0"]/m_["Pont"] for m_ in ac]
c_bo = [m_["L0"]/m_["Pont"] for m_ in bo]
null_ac = max(max(abs(m_["L1"]), abs(m_["L2"])) for m_ in ac)/skala(ac)
null_bo = max(max(abs(m_["L1"]), abs(m_["L2"])) for m_ in bo)/skala(bo)
ZAHLEN["zellen"] = {"A, ohne Einschub: (Riem^2, lambda R, lambda^2)": [round(c_ao[0][0], 9), round(c_ao[1][0], 9), round(c_ao[2][0], 9)],
                    "B, Chiralität: (Euler, lambda R, lambda^2)": [round(c_bc[0][0], 9), round(c_bc[1][0], 9), round(c_bc[2][0], 9)],
                    "A, Chiralität: L0/Pontrjagin, (L1, L2)/Skala": [round(c_ac[0], 9), float(null_ac)],
                    "B, ohne Einschub: L0/Pontrjagin, (L1, L2)/Skala": [round(c_bo[0], 9), float(null_bo)],
                    "Streuungen": [konst(c_ao[0]), konst(c_ao[1]), konst(c_ao[2]), konst(c_bc[0]), konst(c_bc[1]),
                                   konst(c_bc[2]), konst(c_ac), konst(c_bo)], "zyklische Identität": zykl_max}
streu = max(ZAHLEN["zellen"]["Streuungen"])
check(f"[E] Vier Zellen für F^ab = R^ab + lambda e^a ^ e^b (fünf zufällige algebraische Krümmungstensoren mit zyklischer "
      f"Identität, zufällige Vierbeine, g = e^T eta e; Dichten je det e; Streuung {e0(streu)}). (A, ohne Einschub) "
      f"Tr(F_mn F^mn) = {c_ao[0][0]:.0f} Riem^2 + ({c_ao[1][0]:.0f}) lambda R + ({c_ao[2][0]:.0f}) lambda^2, also "
      f"-8 (Riem^2 + 4 lambda R + 24 lambda^2): Einstein-Hilbert mit Lambda = -3 lambda, dazu Riemann-Quadrat. (B, "
      f"Chiralität) Tr(W F ^ F) = {c_bc[0][0]:.0f} Euler + {c_bc[1][0]:.0f} lambda R + {c_bc[2][0]:.0f} lambda^2: "
      f"MacDowell-Mansouri, ein Viertel der z118-Werte 16 : 64 : 384, weil F hier die Spinorkrümmung (1/2) F^ab Sigma_ab "
      f"ist, wieder Lambda = -3 lambda. (A, Chiralität) eps_abcd F^ab_mn F^cd,mn = {c_ac[0]:.0f} mal Pontrjagin, der "
      f"Einstein-Term ist der Holst-Term und verschwindet ({e0(null_ac)}). (B, ohne Einschub) F^ab ^ F_ab = "
      f"{c_bo[0]:.0f} mal Pontrjagin, kein Einstein-Term ({e0(null_bo)}). Einstein-Dynamik gibt es also genau in den "
      f"Zellen (Hodge, Skalartyp) und (Keil, Chiralitätstyp); der Eichterm F_mn F^mn steht nur in der ersten",
      streu < 1e-8 and t_real < 1e-9 and zykl_max < 1e-10 and abs(c_ao[0][0] + 8) < 1e-8 and abs(c_ao[1][0] + 32) < 1e-8
      and abs(c_ao[2][0] + 192) < 1e-7 and abs(abs(c_bc[0][0]) - 4) < 1e-8 and abs(abs(c_bc[1][0]) - 16) < 1e-8
      and abs(abs(c_bc[2][0]) - 96) < 1e-7 and abs(c_bc[2][0]/c_bc[1][0] - 6) < 1e-8 and null_ac < 1e-9 and null_bo < 1e-9
      and abs(c_ac[0]) > 0.1 and abs(c_bo[0]) > 0.1)

say("  3.2 Normierung beider Einstein-Zellen")
kap0, kaph, lam_s, G_s, g_s, sig_s, l_s = sp.symbols("kappa_0 kappahat lambda G g sigma l", real=True)
ALPHA_K, GAMMA_E, BETA_E = 16, 8, 16
# (A): S_E = -(kappa_0/4) Int Tr_V(F_mn F^mn); Eichung: Tr_V = -gamma F^i F^i => 1/g^2 = kappa_0 gamma
inv_g2 = kap0*GAMMA_E
# Gravitation: -(kappa_0/4) (-alpha/2) 4 lambda R = (alpha kappa_0/2) lambda R  =  -(1/(16 pi G)) R
G_A = sp.solve(sp.Eq(sp.Rational(ALPHA_K, 2)*kap0*lam_s, -1/(16*sp.pi*G_s)), G_s)[0]
l2G_A = sp.simplify((l_s**2/G_A).subs(lam_s, sig_s/l_s**2))
verh_A = sp.simplify(l2G_A/inv_g2)
GLam_A = sp.simplify((G_A*(-3*lam_s)).subs(kap0, 1/(GAMMA_E*g_s**2)))
# (B): S = kappahat Int Tr_K(W F ^ F): kappahat beta lambda = 1/(16 pi G)
G_B = sp.solve(sp.Eq(kaph*BETA_E*lam_s, 1/(16*sp.pi*G_s)), G_s)[0]
l2G_B = sp.simplify((l_s**2/G_B).subs(lam_s, sig_s/l_s**2))
kap_z127 = sp.Symbol("kappa", positive=True)
z127_ok = sp.simplify(l2G_B.subs({kaph: 4*kap_z127, sig_s: 1}) - 1024*sp.pi*kap_z127) == 0
verh_B = sp.simplify(l2G_B/inv_g2)
ZAHLEN["normierung"] = {"(A) 1/g^2": str(inv_g2), "(A) l^2/G": str(l2G_A), "(A) (l^2/G)/(1/g^2)": str(verh_A),
                        "(A) G Lambda": str(GLam_A), "(B) l^2/G": str(l2G_B), "(B) (l^2/G)/(1/g^2)": str(verh_B),
                        "(B) mit kappahat = 4 kappa: 1024 pi kappa": bool(z127_ok)}
check(f"[E/C] Normierung. (A) Aus der Gitterspur S = kappa_0 sum_p Re Tr_V(1 - U_p) = -(kappa_0/4) Int Tr_V(F_mn F^mn) "
      f"(z135, isotrop) folgt 1/g_i^2 = {inv_g2} und, mit dem euklidischen Einstein-Hilbert-Term -(1/16 pi G) R, "
      f"l^2/G = {l2G_A}, also (l^2/G) : (1/g_U^2) = {verh_A} und G Lambda = {GLam_A}. kappa_0 fällt heraus. Die Zahl ist "
      f"8 pi alpha/gamma mit alpha = dim(16) Tr_T-Normierung = 16 und gamma = dim(T) Tr_16 T^2 = 8: gleichwertig "
      f"k_grav : k_i = 2 alpha : 2 gamma = 32 : 16 = 2 : 1, das Gegenstück zu k_3 = k_2 = k_1. (B) Aus S = kappahat Int "
      f"Tr_K(R F ^ F) folgt l^2/G = {l2G_B}; mit kappahat = 4 kappa ist das z127s 1024 pi kappa ({bool(z127_ok)}). Das "
      f"Verhältnis ist dort {verh_B}: es enthält das relative Gewicht rho = kappahat/kappa_0 der beiden Funktionale. "
      f"Bedingungen für (A): P_exp (Antwort = Wirkung), P_hodge (der Stern des Netzes ist der Stern der Vierbeinmetrik)",
      sp.simplify(verh_A + 16*sp.pi*sig_s) == 0 and sp.simplify(GLam_A - 3*g_s**2/(16*sp.pi)) == 0 and z127_ok
      and sp.simplify(verh_B - 32*sp.pi*sig_s*kaph/kap0) == 0)
print(f"   (Laufzeit bisher {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)

# ================================================================================================== Teil 4
say("\nTeil 4: Das Verhältnis (l^2/G) : (1/g_U^2) je Einschub und Lesart")
say("  4.1 Die drei Lesarten aus z135 auf dem gemeinsamen Träger")
DV = 64
FS, FM = zufall_grav(1.0), zufall_eich(1.0)
trV2 = 16*np.trace(FS @ FS) + 4*np.trace(FM @ FM)


def antworten(eps):
    gS, gM = expm(eps*FS), expm(eps*FM)
    gSi, gMi = np.linalg.inv(gS), gM.conj().T
    tg, tgi = np.trace(gS)*np.trace(gM), np.trace(gSi)*np.trace(gMi)
    r1 = DV - tg                                                     # Tr_V(1 - U_p)
    r2 = (DV**2 - tg*tgi)/(4*DV**4)                                  # W(g) mit Tr g Tr g^-1
    t5, t5i = np.trace(W5S @ gS)*np.trace(gM), np.trace(np.linalg.inv(W5S) @ gSi)*np.trace(gMi)
    r1c = -t5                                                        # Tr_V(W_5 (1 - U_p)), Tr W_5 = 0
    r2c = -(t5*t5i)/(4*DV**4)
    return r1, r2, r1c, r2c


EPSQ = 1e-3
a1, a2, a1c, a2c = antworten(EPSQ)
b1, b2, b1c, b2c = antworten(EPSQ/2)
les_num = {"R1: Tr_V(1 - U_p)/(-eps^2 Tr_V F^2/2)": complex(a1/(-EPSQ**2*trV2/2)),
           "R2: W(g)/(-eps^2 Tr_V F^2/(4 d^3))": complex(a2/(-EPSQ**2*trV2/(4*DV**3))),
           "R1 mit Chiralität: Antwort/(-eps^2 Tr_T(W F_g^2) 16/2)": complex(a1c/(-EPSQ**2*16*np.trace(W5S @ FS @ FS)/2)),
           "R2 mit Chiralität: Ordnung in eps": float(np.log2(abs(a2c)/abs(b2c)))}
gam_c, alp_c = 0.0, 0.0                                              # Chiralitätstyp: gamma = 0, alpha = 0
LESARTEN = [
    ("R1 Netzspur je Periode (kappa_0 = sin^4 theta = 1/4)", Fr(1, 4)),
    ("R2 W je Etikettenzustand (kappa_0 = 1/(2 d^3), d = 64)", Fr(1, 2*DV**3)),
    ("R3a Niveau mit P_norm, Ströme auf dem gemeinsamen Träger (1/g^2 = 2 gamma)", Fr(2)),
]
tab_les = []
for name, k0 in LESARTEN:
    ig2 = k0*GAMMA_E
    l2G = 8*sp.pi*ALPHA_K*k0                                         # Betrag; Vorzeichen -sigma
    tab_les.append({"Lesart": name, "kappa_0": str(k0), "1/g_U^2": str(ig2), "alpha_U^-1": float(4*np.pi*ig2),
                    "|l^2/G|": float(l2G), "(|l^2/G|)/(1/g_U^2) / pi": str(sp.nsimplify(l2G/ig2/sp.pi))})
# R3b: getrennte Ströme: k_i = 2 Tr_16 T^2 = 4 auf der 16, k_grav = 2 Tr_T (i Sigma_12)^2 = 2 auf T
k_grav_T = float(2*np.trace((1j*ST[0]) @ (1j*ST[0])).real)
k_i_16 = 4
verh_R3b = sp.nsimplify(8*sp.pi*(k_grav_T/2)/(k_i_16/2))             # 8 pi alpha/gamma mit alpha = 1, gamma = 2
tab_les.append({"Lesart": "R3b Niveau mit P_norm, getrennte Ströme (k_i = 4 auf der 16, k_grav = 2 auf T)", "kappa_0": "1/2",
                "1/g_U^2": "4", "alpha_U^-1": float(16*np.pi), "|l^2/G|": float(4*verh_R3b),
                "(|l^2/G|)/(1/g_U^2) / pi": str(sp.nsimplify(verh_R3b/sp.pi))})
ZAHLEN["lesarten"] = {"numerisch": {k: (rz(v) if isinstance(v, complex) else v) for k, v in les_num.items()}, "Tabelle": tab_les,
                      "k_grav auf T": k_grav_T}
for z in tab_les:
    say(f"     {z['Lesart']:<82} | 1/g_U^2 = {z['1/g_U^2']:<9} | alpha_U^-1 = {z['alpha_U^-1']:<10.6g} | |l^2/G| = "
        f"{z['|l^2/G|']:<11.6g} | Verhältnis = {z['(|l^2/G|)/(1/g_U^2) / pi']} pi")
verh_werte = [z["(|l^2/G|)/(1/g_U^2) / pi"] for z in tab_les]
check(f"[E/C] Lesarten des Absolutwerts auf dem gemeinsamen Träger V = T16 (d = 64): Die Netzspur (R1) und die "
      f"Schleifenantwort W(g) = (d^2 - Tr g Tr g^-1)/(4 d^4) (R2) haben in zweiter Ordnung dieselbe Invariante Tr_V F^2 "
      f"(numerisch {rz(les_num['R1: Tr_V(1 - U_p)/(-eps^2 Tr_V F^2/2)'], 5)} und "
      f"{rz(les_num['R2: W(g)/(-eps^2 Tr_V F^2/(4 d^3))'], 5)}), ebenso die Niveaulesart über den gemeinsamen Träger (R3a). "
      f"In allen dreien fällt lambda heraus: (l^2/G) : (1/g_U^2) = {verh_werte[0]} pi, während 1/g_U^2 selbst um den "
      f"Faktor {float(Fr(tab_les[2]['1/g_U^2'])/Fr(tab_les[1]['1/g_U^2'])):.3g} streut. Abhängig ist das Verhältnis von einer "
      f"Wahl: Zählt man die Niveaus getrennter Ströme (R3b: k_i = 4 auf der 16, k_grav = {k_grav_T:g} auf T), ist es "
      f"{verh_werte[3]} pi. Mit dem Chiralitätseinschub gibt R1 nur den eps-Term (Verhältnis zur Formel "
      f"{rz(les_num['R1 mit Chiralität: Antwort/(-eps^2 Tr_T(W F_g^2) 16/2)'], 5)}) und keinen Eichterm; in R2 beginnt die "
      f"Antwort erst in Ordnung {les_num['R2 mit Chiralität: Ordnung in eps']:.2f} (Quadrat des eps-Terms), in R3 gibt es "
      f"kein Niveau, weil eps_abcd keine Killing-Form ist",
      abs(les_num["R1: Tr_V(1 - U_p)/(-eps^2 Tr_V F^2/2)"] - 1) < 1e-5 and abs(les_num["R2: W(g)/(-eps^2 Tr_V F^2/(4 d^3))"] - 1) < 1e-5
      and abs(les_num["R1 mit Chiralität: Antwort/(-eps^2 Tr_T(W F_g^2) 16/2)"] - 1) < 1e-4
      and abs(les_num["R2 mit Chiralität: Ordnung in eps"] - 4) < 0.01 and verh_werte[:3] == ["16", "16", "16"]
      and verh_werte[3] == "4" and abs(k_grav_T - 2) < TOL)

say("  4.2 Tabelle: Einschub mal Struktur mal Lesart")
MATRIX = [
    ("ohne Einschub / Skalartyp", "A Hodge", "Eichterm ja; Gravitation: Killing-Form, EH + Lambda + Riem^2 + Torsion^2", "16 pi (R1, R2, R3a); 4 pi (R3b)"),
    ("ohne Einschub / Skalartyp", "B Keil", "Eichung: theta-Term; Gravitation: Pontrjagin; keine Dynamik", "keines (0 : 0)"),
    ("Chiralitätstyp (MM)", "A Hodge", "Eichterm 0; Gravitation: Pontrjagin, Holst = 0", "keines"),
    ("Chiralitätstyp (MM)", "B Keil", "Eichterm 0; Gravitation: MacDowell-Mansouri (Euler + EH + Lambda)", "32 pi rho, rho = kappahat/kappa_0 frei"),
]
for z in MATRIX:
    say(f"     {z[0]:<26} | {z[1]:<8} | {z[2]:<76} | {z[3]}")
ZAHLEN["matrix"] = [{"Einschub": z[0], "Struktur": z[1], "Inhalt": z[2], "Verhältnis": z[3]} for z in MATRIX]
rho_min = float(SS**8/(4*SS**4*(CC**4 + SS**4)/2))
ZAHLEN["rho_illustration"] = {"rho = a_B/(4 a_A) im dualen Ebenenpaar bei pi/4": rho_min, "Verhältnis 32 pi rho / pi": 32*rho_min}
befund("Verhältnis nativ", "bedingt",
       f"E/C: Eine reine Zahl ohne lambda gibt es nur in der Zelle (Hodge, ohne Einschub): (l^2/G) : (1/g_U^2) = 16 pi, "
       f"gleichwertig g_grav^2 = g_U^2/2 für die so(3,2)-Yang-Mills-Kopplung. Sie ist in den Lesarten R1, R2, R3a "
       f"dieselbe und hängt nur an der Wahl 'gemeinsamer Träger' (sonst 4 pi). Für MacDowell-Mansouri (Keil, "
       f"Chiralitätstyp) gibt es kein Verhältnis: Der Eichterm fehlt dort, und das relative Gewicht rho der zwei "
       f"Funktionale ist frei. Zur Einordnung, nicht als Ergebnis: Im dualen Ebenenpaar wäre unter drei weiteren Wahlen "
       f"(gleicher Betrag beider Gewichte, isotrope Ergänzung auf drei duale Paare mit eps-Vorzeichen, lineare Lesart) "
       f"rho = a_B/(4 a_A) = {rho_min:g}, also {32*rho_min:g} pi.")

say("  4.3 Zahlenwert mit der gelösten Feinstrukturkonstante; Überbestimmungstabelle")
aU_ohne = 3*ALPHA_INV/8                                              # alpha_U^-1 = (3/8) alpha^-1 bei sin^2 = 3/8, ohne Lauf


def g_r2(aU, verh_pi=16):
    """G/R_Vak^2 = G/(4 l^2) mit l^2/G = verh * alpha_U^-1/(4 pi)."""
    return 1.0/(4*verh_pi*np.pi*aU/(4*np.pi))


z136 = json.loads((HERE / "z136_ergebnis.json").read_text(encoding="utf-8"))
ALT = [z for z in z136["zahlen"]["tabelle"] if z["G/R_Vak^2"] is not None]
TARGET = float(next(z["G/R_Vak^2"] for z in ALT if z["Quelle"].startswith("(b')")))
NEU = [
    ("gemeinsame Spur, alpha als Zahl: alpha_U^-1 = (3/8) 137,036 (ohne Lauf)", g_r2(aU_ohne), "C (P_exp, P_hodge, P_plain)"),
    ("gemeinsame Spur, alpha_U^-1 = 60 (z25-Fenster, unten)", g_r2(60.0), "C, Lauf grob"),
    ("gemeinsame Spur, alpha_U^-1 = 74 (z25-Fenster, oben)", g_r2(74.0), "C, Lauf grob"),
    ("gemeinsame Spur, Lesart R1 selbst (1/g_U^2 = 2, alpha_U^-1 = 8 pi)", g_r2(8*np.pi), "C/L"),
    ("gemeinsame Spur, Lesart R2 selbst (1/g_U^2 = 1/65536)", g_r2(4*np.pi/65536), "C/L"),
    ("gemeinsame Spur, Lesart R3a selbst (1/g_U^2 = 16)", g_r2(64*np.pi), "C/L"),
    ("getrennte Ströme R3b (Verhältnis 4 pi, 1/g_U^2 = 4)", g_r2(16*np.pi, 4), "C/L"),
    ("getrennte Ströme, alpha als Zahl (Verhältnis 4 pi)", g_r2(aU_ohne, 4), "C/L"),
]
say("   Überbestimmungstabelle (reine Zahl |G|/R_Vak^2; vorhandene Zeilen aus z136_ergebnis.json, nur gelesen):")
AUSWAHL = ("(a) z118", "(b')", "(d) z112")
zeilen_alt = [z for z in ALT if z["Quelle"].startswith(AUSWAHL)] + [
    z for z in ALT if z["Quelle"].startswith("(e)") and "flach" in z["Lesart"] and "Leiterfeld" not in z["Lesart"]]
for z in zeilen_alt:
    w = float(z["G/R_Vak^2"])
    say(f"     {z['Quelle']:<44} | {z['Lesart']:<86} | {w:<12.6g} (10^{np.log10(w):+.2f}) | {z['Status']}")
Z137 = [("(f) z137 induziert", "eine Leiter, c = 1/120: pi/(16 c L)", np.pi/(16*(1/120))), ("(f) z137 induziert", "3 x 16 Weyl, c = 6/5", np.pi/(16*1.2)),
        ("(f) z137 induziert", "3 x 16 Weyl + 12 Eichbosonen, c = 12/5", np.pi/(16*2.4))]
for q, l_, w in Z137:
    say(f"     {q:<44} | {l_ + ' (mal 1/L, L frei)':<86} | {w:<12.6g} (10^{np.log10(w):+.2f}) | C (P_ind, P_cut)")
for l_, w, st in NEU:
    say(f"     {'(g) z150 Schleifenantwort, Hodge':<44} | {l_:<86} | {w:<12.6g} (10^{np.log10(w):+.2f}) | {st}")
say(f"     {'(g) z150 Schleifenantwort, Keil (MM)':<44} | {'1/(4096 pi kappa) mit kappa = rho kappa_0/4, rho frei':<86} | "
    f"{'--':<12}            | O")
alle_alt = [(z["Quelle"] + " | " + z["Lesart"], float(z["G/R_Vak^2"])) for z in ALT] + [(q + " | " + l_, w) for q, l_, w in Z137]
gleich = [(l_, a_) for l_, w, _ in NEU for a_, wa in alle_alt if abs(w/wa - 1) < 1e-9]
naechst = min(alle_alt, key=lambda t: abs(np.log(t[1]/NEU[0][1])))
kap_aequiv = aU_ohne/(256*np.pi)                                     # 1/(4096 pi kappa) = 1/(16 alpha_U^-1)
aU_2c3 = 256*np.pi*(1/(4*np.pi))
l2G_ohne = 4*aU_ohne
l2G_ziel_log = 2*ALPHA_INV/np.log(10) + np.log10(32*np.pi**3)
ZAHLEN["zahlenwert"] = {"alpha_U^-1 = (3/8) alpha^-1": aU_ohne, "l^2/G = 4 alpha_U^-1 = (3/2) alpha^-1": l2G_ohne,
                        "G/R_Vak^2 = alpha_U/16": NEU[0][1], "log10": float(np.log10(NEU[0][1])),
                        "R_Vak/l_P = 4 sqrt(alpha_U^-1)": float(4*np.sqrt(aU_ohne)), "G Lambda = 3 alpha_U/4": 0.75/aU_ohne,
                        "Ziel G/R_Vak^2 (z136, COSMO.LAM.01)": TARGET, "Abstand in Größenordnungen": float(np.log10(NEU[0][1]/TARGET)),
                        "log10 l^2/G: nativ, Ziel 32 pi^3 e^(2/alpha)": [float(np.log10(l2G_ohne)), float(l2G_ziel_log)],
                        "neue Zeilen": [{"Lesart": l_, "G/R_Vak^2": w, "Status": st} for l_, w, st in NEU],
                        "z137-Zeilen (mal 1/L)": [{"Lesart": l_, "Wert": w} for _, l_, w in Z137],
                        "exakt gleiche vorhandene Zeilen": [[a, b] for a, b in gleich],
                        "nächste vorhandene Zeile": [naechst[0], naechst[1]],
                        "z127-Gewicht kappa, das denselben Wert gäbe": kap_aequiv,
                        "alpha_U^-1, das dem z127-Gewicht 2 c_3 entspricht": aU_2c3}
check(f"[N/C] Zahlenwert in der Zelle (Hodge, ohne Einschub) mit alpha nur als Zahl: alpha_U^-1 = (3/8) alpha^-1 = "
      f"{aU_ohne:.4f} an der Anpassungsskala (ohne Lauf) gibt l^2/G = 4 alpha_U^-1 = (3/2) alpha^-1 = {l2G_ohne:.3f}, "
      f"G/R_Vak^2 = alpha_U/16 = {NEU[0][1]:.6g} = 10^{np.log10(NEU[0][1]):.2f}, R_Vak = {4*np.sqrt(aU_ohne):.1f} "
      f"Plancklängen und G Lambda = 3 alpha_U/4 = {0.75/aU_ohne:.5f}. Der Lauf ändert wenig: Das z25-Fenster alpha_U^-1 = "
      f"60 ... 74 gibt {NEU[1][1]:.3g} ... {NEU[2][1]:.3g}. Das ist wieder ein Vakuum von Planckgröße, in derselben Gruppe "
      f"wie z127 (7,8e-5 bis 2,0e-3) und z112 (ab 0,00995); nächste vorhandene Zeile: {naechst[0]} = {naechst[1]:.6g}. "
      f"Zum Ledger-Ziel {TARGET:.4g} fehlen {np.log10(NEU[0][1]/TARGET):.1f} Größenordnungen. Der Grund ist strukturell: "
      f"Das Verhältnis ist linear in 1/alpha (log10 l^2/G = {np.log10(l2G_ohne):.2f}), das Ziel exponentiell (32 pi^3 "
      f"e^(2/alpha), log10 = {l2G_ziel_log:.2f}). Ein Spurverhältnis auf Baumniveau kann e^(2/alpha) nicht liefern. "
      f"Exakt gleiche vorhandene Zeilen: {len(gleich)}"
      + (f" ({'; '.join(a + ' = ' + b for a, b in gleich)}; Zahlen der Form 1/(2^k pi) unter verschiedenen "
         f"Voraussetzungen, kein Indiz)" if gleich else "")
      + f". Einordnung ohne Anspruch: Derselbe Wert entstünde in z127 mit dem Gewicht kappa = alpha_U^-1/(256 pi) = "
        f"{kap_aequiv:.4f}; das dortige 2 c_3 entspräche alpha_U^-1 = {aU_2c3:.0f}. Das sind zwei verschiedene Funktionale "
        f"(Hodge gegen Keil), und 2 c_3 ist kein natives Spurgewicht (z127)",
      abs(NEU[0][1] - 1/(16*aU_ohne)) < 1e-15 and 1e-4 < NEU[0][1] < 1e-2 and np.log10(NEU[0][1]/TARGET) > 119
      and abs(aU_2c3 - 64) < 1e-9 and l2G_ziel_log > 120)
print(f"   (Laufzeit bisher {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)

# ================================================================================================== Teil 5
say("\nTeil 5: Signatur und Vorzeichen")
hh_e = [E6[0], E6[1], E6[2], E6[3], 1j*E6[4], 1j*E6[5], E6[4]]        # Index 6: E6[4] mit Quadrat +1 (für so(5))
FAELLE = [("so(3,2): AdS, Lorentz-Signatur (nativ)", (0, 1, 2, 4, 5)), ("so(4,1): dS, Lorentz-Signatur", (0, 1, 2, 4, 3)),
          ("so(4,1): AdS, euklidischer Schnitt", (0, 1, 2, 3, 5)), ("so(5): dS, euklidischer Schnitt", (0, 1, 2, 3, 6))]
sig_tab = {}
for name, idx in FAELLE:
    eta = np.array([np.round((hh_e[a] @ hh_e[a])[0, 0].real) for a in idx])
    gen = {(i, j): ((hh_e[idx[i]] @ hh_e[idx[j]] - hh_e[idx[j]] @ hh_e[idx[i]])/4)[np.ix_(idxT, idxT)] for i, j in PAIRS5}
    Q = np.array([[16*np.trace(gen[p] @ gen[q]) for q in PAIRS5] for p in PAIRS5])
    ew = np.linalg.eigvalsh((Q + Q.conj().T).real/2)
    lor = sorted({int(np.sign(round(Q[m_, m_].real))) for m_, p in enumerate(PAIRS5) if p in LOR})
    tra = sorted({int(np.sign(round(Q[m_, m_].real))) for m_, p in enumerate(PAIRS5) if p in TRA})
    kom = prop(gen[(0, 4)] @ gen[(1, 4)] - gen[(1, 4)] @ gen[(0, 4)], gen[(0, 1)])
    sigma = int(round(kom[0].real))                                   # lambda = sigma/l^2
    sig_tab[name] = {"eta": [int(x) for x in eta], "Gram (+, -)": [int(np.sum(ew > 1e-9)), int(np.sum(ew < -1e-9))],
                     "Lorentzblock": lor, "Translationsblock": tra, "sigma": sigma, "Lambda": -sigma,
                     "sgn(G g^2) = -sigma": -sigma, "Rest": float(np.abs(Q.imag).max() + kom[1])}
ZAHLEN["signatur"] = sig_tab
for name, v in sig_tab.items():
    say(f"     {name:<40} | eta = {str(v['eta']):<18} | Gram (+, -) = {str(v['Gram (+, -)']):<8} | Lorentz {str(v['Lorentzblock']):<8} "
        f"| Translation {str(v['Translationsblock']):<5} | sigma = {v['sigma']:+d} | sgn Lambda = {v['Lambda']:+d} | sgn(G g^2) = "
        f"{v['sgn(G g^2) = -sigma']:+d}")
s_nat, s_e_ads, s_e_ds = (sig_tab[FAELLE[0][0]], sig_tab[FAELLE[2][0]], sig_tab[FAELLE[3][0]])
check(f"[E] Signatur der einschubfreien Invariante Tr_V F_grav^2 = 16 mal Killing-Form: Für die native so(3,2) hat sie "
      f"{s_nat['Gram (+, -)'][0]} positive (Boosts und räumliche Translationen, nichtkompakt) und {s_nat['Gram (+, -)'][1]} "
      f"negative Richtungen (Drehungen und die Zeittranslation, kompakt); der Eichterm Tr_V(i T)^2 = -8 ist negativ wie die "
      f"kompakten. Auf dem euklidischen Schnitt der AdS-Gruppe ist der Lorentzblock {s_e_ads['Lorentzblock']} (wie die "
      f"Eichung) und der Translationsblock {s_e_ads['Translationsblock']} (entgegengesetzt): Riem^2 hat das Vorzeichen des "
      f"Eichterms, Torsion^2 das falsche. Nur für die kompakte so(5) (dS) sind alle zehn Richtungen gleichsinnig "
      f"({s_e_ds['Gram (+, -)']}). Der Einstein-Term trägt sgn(G g^2) = -sigma: Bei positivem g^2 ist G für die native "
      f"AdS-Gruppe (sigma = +1, Lambda < 0) negativ und nur für dS positiv. In jedem Fall gilt G Lambda = 3 g^2/(16 pi) > 0",
      s_nat["Gram (+, -)"] == [6, 4] and s_nat["sigma"] == 1 and s_e_ads["Lorentzblock"] == [-1]
      and s_e_ads["Translationsblock"] == [1] and s_e_ads["sigma"] == 1 and s_e_ds["Gram (+, -)"] == [0, 10]
      and s_e_ds["sigma"] == -1 and sig_tab[FAELLE[1][0]]["Gram (+, -)"] == [4, 6] and sig_tab[FAELLE[1][0]]["sigma"] == -1
      and all(v["Rest"] < 1e-9 for v in sig_tab.values()))
befund("Vorzeichen", "negativ",
       "E: In der einzigen Zelle mit nativem Verhältnis ist das Vorzeichen falsch. Die gemeinsame Spur koppelt die drei "
       "Vorzeichen: sgn(G) sgn(g^2) = sgn(Lambda). Mit der nativen so(3,2) (Lambda < 0, z98/z104/z118) ist entweder G "
       "negativ oder die Eichfelder sind Geister. Dazu kommen Riem^2 (nicht topologisch, anders als der Euler-Term in MM) "
       "und ein Torsionsquadrat mit falschem Vorzeichen; auf Netzebene zeigt sich dasselbe als Nichtunitarität (Teil 1.4). "
       "In MM (Keil, Chiralität) ist das Gesamtvorzeichen frei (Paarungsbit, z127), Lambda = -3 lambda < 0 liegt fest, und "
       "nur der topologische Euler-Term ist quadratisch in der Krümmung; dort fehlt aber der Eichterm.")

# ================================================================================================== Teil 6
say("\nTeil 6: Bilanz und Urteil")
bilanz = [
    ("E", "Modell",
     "Zwei native Träger der Dimension 64 in S^+: T16 (Twistor mal 16, volle so(3,2), Uhr skalar) und der Indexkern K "
     "(nur Lorentz, Uhr = Chiralität). Der Link g_grav x g_eich gibt wörtlich den Rückkehrblock cos^4 + sin^4 U_p; "
     "Tr U_p ist ein Produkt. Das Netz ist für Spin(3,2) nur pseudo-unitär."),
    ("E", "Zweite Ordnung",
     "Ohne Einschub: 16 mal Killing-Form von so(3,2) plus 8 delta_ij für die Eichung; k_i = 16, 3/8 unverändert; kein "
     "Mischterm bis zur dritten Ordnung, in vierter das Produkt 6 Tr(W F_g^2) Tr(F_e^2)."),
    ("E", "Native Einschübe",
     "Die vier endlichen modularen Gewichte sind je Träger ein Skalar +-i (Killing plus Eichung, imaginär) oder die "
     "Chiralität +-i gamma_5 (eps-Term reell, Eichterm exakt null, erst Tr X^5). Beides zugleich nur mit der "
     "Holst-Mischung, die nicht modular ist (z127)."),
    ("E", "Struktur",
     "Einstein-Dynamik nur in (Hodge, Skalartyp): -8 (Riem^2 + 4 lambda R + 24 lambda^2), und in (Keil, "
     "Chiralitätstyp): MacDowell-Mansouri 4 : 16 : 96. Die anderen zwei Zellen sind Pontrjagin. Der wörtliche Block aus "
     "z135 ist Hodge; den Keil gibt erst ein duales Ebenenpaar (a_B = a_A = 1/16 bei pi/4)."),
    ("E/C", "Verhältnis",
     "(l^2/G) : (1/g_U^2) = 16 pi in (Hodge, ohne Einschub), gleich in R1, R2, R3a; 4 pi bei getrennten Strömen. Für "
     "MM 32 pi rho mit freiem rho. Gezählte Wahlen: Einschubtyp (2), Struktur (2), Träger gemeinsam oder getrennt (2), "
     "für MM zusätzlich rho (kontinuierlich) und die 4D-Anordnung."),
    ("N/C", "Wert",
     f"Mit alpha als Zahl: l^2/G = (3/2) alpha^-1 = {l2G_ohne:.2f}, G/R_Vak^2 = alpha_U/16 = {NEU[0][1]:.4g}; "
     f"{np.log10(NEU[0][1]/TARGET):.1f} Größenordnungen über dem Ledger-Ziel; linear statt exponentiell in 1/alpha."),
    ("E negativ", "Vorzeichen",
     "sgn(G) sgn(g^2) = sgn(Lambda): Mit der nativen AdS-Gruppe ist G negativ. Riem^2 ist nicht topologisch, "
     "Torsion^2 hat das falsche Vorzeichen."),
    ("O", "Offen",
     "Ein natives Prinzip für rho (Summe über Modulgewichte mit festem Betrag?), eine native 4D-Anordnung, und ob "
     "ein nichtperturbatives Gewicht e^(-S) mit S = 2/alpha = 8 pi/e^2 aus dem Netz folgt."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
ZAHLEN["bilanz"] = [{"status": st, "thema": was, "aussage": text} for st, was, text in bilanz]
urteil = ("Die Strategie reduziert P_G nicht auf die Eichkopplung; sie scheitert für MacDowell-Mansouri und ist für die "
          "Yang-Mills-artige Gravitation nur bedingt und mit falschem Vorzeichen erfolgreich. (1) Das Modell ist "
          "wörtlich durchführbar: Ein Link in Spin(3,2) x G_SM auf dem Träger Twistor mal 16 gibt denselben "
          "Rückkehrblock cos^4 + sin^4 U_p wie in z135, und Tr U_p ist das Produkt der beiden Spuren. (2) In zweiter "
          "Ordnung gibt die Spur ohne Einschub 16 mal die Killing-Form von so(3,2) plus 8 delta_ij für die Eichung, "
          "ohne Mischterm. Daraus folgt eine native reine Zahl, in der lambda herausfällt: g_grav^2 = g_U^2/2, "
          "gleichwertig (l^2/G) : (1/g_U^2) = 16 pi. Sie ist in allen drei Lesarten des Absolutwerts aus z135 dieselbe, "
          "solange Gravitation und Eichung über denselben Träger gespurt werden; bei getrennten Strömen ist sie 4 pi. "
          "(3) Das ist aber nicht die MacDowell-Mansouri-Wirkung aus z118/z127. Der eps-Term braucht den "
          "Chiralitätseinschub (auf dem Indexkern die Uhr), und mit ihm ist der Eichterm exakt null; die Eichantwort "
          "beginnt dann erst mit dem Pfaffschen Term fünfter Ordnung. Umgekehrt gibt der Skalartyp den Eichterm, aber "
          "keinen eps-Term. Kein natives endliches Gewicht enthält beides; das wäre die Holst-Mischung, die z127 als "
          "nicht modular ausgeschlossen hat. Zusätzlich braucht MM die Keil-Struktur zweier dualer Plaketten, der "
          "Eichterm die Hodge-Struktur derselben Plakette. Für MM bleibt deshalb ein freies relatives Gewicht rho, und "
          "kappa wird nicht festgelegt. (4) In der Zelle mit nativem Verhältnis ist die Gravitation Yang-Mills-artig: "
          "Einstein-Hilbert mit Lambda = -3 lambda, dazu ein nicht topologisches Riemann-Quadrat und ein "
          "Torsionsquadrat. Die Vorzeichen sind gekoppelt, sgn(G) sgn(g^2) = sgn(Lambda): Mit der nativen AdS-Gruppe "
          "so(3,2) ist G negativ, oder die Eichfelder sind Geister; nur die dS-Gruppe wäre konsistent. Auf Netzebene "
          "ist das die Nichtunitarität des Netzes für nichtkompakte Links. (5) Der Zahlenwert mit alpha als Zahl ist "
          f"l^2/G = (3/2) alpha^-1 = {l2G_ohne:.2f}, also G/R_Vak^2 = alpha_U/16 = {NEU[0][1]:.4g}: wieder ein Vakuum "
          f"von Planckgröße, {np.log10(NEU[0][1]/TARGET):.1f} Größenordnungen über dem Ledger-Ziel. Das Verhältnis ist "
          "linear in 1/alpha, das Ziel exponentiell; ein Spurverhältnis kann diese Lücke nicht schließen. Prämissen "
          "der bedingten Aussage: P_link (jetzt mit nichtkompakter Gruppe), P_exp, P_hodge (neu: Stern des Netzes = "
          "Stern der Vierbeinmetrik), P_plain (neu: das Gewicht vom Skalartyp reell gelesen; nativ ist es auf T16 die "
          "Uhr mit der Phase i, die einschubfreie Spur selbst ist kein endliches Gewicht) und P_gem (neu: gemeinsamer "
          "Träger statt getrennter Ströme). Für K13 heißt das: bleibt L; ein fünfter Weg mit Planckwert und nativ falschem Vorzeichen. Für K15: "
          "unverändert L; k_i und 3/8 überleben die Erweiterung des Trägers. P_G bleibt eine freie Zahl. Nichts ist "
          "befördert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z150_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n", encoding="utf-8")
(HERE / "ausgabe_z150.txt").write_text("\n".join(OUT) + "\n", encoding="utf-8")
