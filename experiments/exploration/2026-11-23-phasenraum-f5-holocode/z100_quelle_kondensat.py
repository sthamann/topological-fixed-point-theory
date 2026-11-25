"""Z100: Eine Quelle für Krümmung und Wechselwirkung? Kondensat und Kompositzustände (Folgeauftrag N1, Hypothese H1).

Leitlinie 'Rekonstruktion statt Auswahl'. Firewalls: GRAV.NONCIRCULAR.01 (der Modulfluss bleibt Uhr; die Tomita-Familie ist
Kinematik, keine Zeit), GRAV.PARENT.GAP.NO_GO.01, DIMENSION.UPLIFT.FIREWALL.01 (alles ist Fock- und Twistorkinematik einer
bzw. weniger Nähte auf Abschneidungen; keine wechselwirkende 3+1-Schließung). Gelöste Blöcke (alpha, Flavour, Massen, 16, E8,
P2) bleiben unberührt; nichts wird befördert. Bestand: z98 (Tomita-Familie J_g = U(g) J U(g)^-1 unter der nativen AdS-Gruppe
= AdS4, Levi-Civita als Nullmoden-Verbindung, MM-Krümmung F = 0, Modulgenerator auf T zentral), z99 (alle nativen Operatoren
einkörperig, verbundene 1->2-Amplitude erst mit Zweikörperterm), z89 (Impulsphasen-Kondensat exp(beta.c^+ - h.c.) zweier Nähte
wählt det(beta1, beta2) I_pi; Gravitonsektor D = -4, tiefste Zustände |0, 0; m, 4 - m>), z91 (hartkernige Komposition von N
Nähten, kappa_2 = 1/2, kappa_3 = 0, kappa_4 = -1/(2N)). Hypothese H1: Ein nicht gaußscher Zustand gibt zugleich F != 0 mit
festem b (Gravitation) und eine verbundene 1->2-Amplitude mit festem g (Wechselwirkung).
Statuskürzel: E exakt, C bedingt, L Lesart, O offen, N numerisch.

Kandidaten: (a1) Impulsphasen-Kondensat psi_beta = exp(beta.c^+ - h.c.) Omega (z89; zwei Nähte mit det(beta1, beta2) != 0),
(a2) wörtliches D = -4-Kondensat psi_lambda = Omega + lambda h mit h = |0, 0; 2, 2> (tiefstes Graviton, nicht gaußsch),
(b) hartkerniger Kompositzustand von N Nähten (z91), jede Naht im KMS-Zustand.

Teil 1  Gravitationsseite: modulare Daten (Delta, J) auf Abschneidungen (Original a, Partner b, Polarzerlegung); Nullmoden
        der Tomita-Familie J_g = U(g) J_psi U(g)^-1 unter der nativen so(3,2) (Fock-Gleichungssystem, Twistor-Kleingruppen);
        Reduktivität und Symmetrie des Paares (so(3,2), h_psi); Cartan-Verbindung (Maurer-Cartan), MM-Krümmung, Torsion;
        modulare Berry-Krümmung mit den Nullmoden von K_psi; quadratische Antwort (Fisher-Information) und b.
Teil 2  Wechselwirkungsseite: verbundene 1->2-Amplituden nach z99 für alle quadratischen Monome, die nativen Erzeuger und
        K_psi; Paritätsauswahlregel (-1)^D; der quartische Vertex des D = -4-Kondensats; Kompositzustand (einkörperiges K,
        Kumulanten, hartkerniges Leck der nativen AdS-Wirkung).
Teil 3  Beziehung zwischen b (bzw. G) und g, Bilanz (E/C/L/O), Urteil zu H1; Folgen für K13, K15, K20.
"""
import itertools
import json
import math
import time
from pathlib import Path

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sps
import scipy.sparse.linalg as spla

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
rng = np.random.default_rng(100)
I2 = np.eye(2)
SIGMA = [np.eye(2, dtype=complex), np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]),
         np.diag([1., -1.]).astype(complex)]
Q_TEST = 0.25                                                       # Prüfwert; nativ q = e^(-2 pi m) (z89)
Q_NATIV = math.exp(-2*math.pi)


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


def r(x, n=12):
    return float(np.round(float(x), n))


def nullraum(M, tol=1e-9):
    M = np.atleast_2d(M)
    _, s, vh = np.linalg.svd(M)
    rank = int(np.sum(s > tol*max(1.0, s[0]))) if s.size else 0
    return vh[rank:]


# ------------------------------------------------------------------ u(2,2) auf dem Twistorraum (wie z98)
G4 = np.diag([1., 1., -1., -1.]).astype(complex)
S4 = np.block([[np.zeros((2, 2)), I2], [I2, np.zeros((2, 2))]]).astype(complex)       # tau = (Ad J)^+ auf Z
_rows = []
for _k in range(32):
    _v = np.zeros(32)
    _v[_k] = 1
    _M = (_v[:16] + 1j*_v[16:]).reshape(4, 4)
    _C = _M.conj().T @ G4 + G4 @ _M
    _rows.append(np.concatenate([_C.real.ravel(), _C.imag.ravel()]))
U22 = [(v[:16] + 1j*v[16:]).reshape(4, 4) for v in nullraum(np.array(_rows).T)]


def teilalgebra(bed, basis=None):
    basis = U22 if basis is None else basis
    rows = []
    for M in basis:
        parts = []
        for f in bed:
            C = f(M)
            parts += [C.real.ravel(), C.imag.ravel()]
        rows.append(np.concatenate(parts))
    ns = nullraum(np.array(rows).T)
    return [sum(c*M for c, M in zip(v, basis)) for v in ns]


uvec = [np.array([1, 0, -1, 0.]), np.array([0, 1, 0, -1.])]
vvec = [np.array([1, 0, 1, 0.]), np.array([0, 1, 0, 1.])]
Inp = (np.outer(uvec[0], uvec[1]) - np.outer(uvec[1], uvec[0])).astype(complex)       # I_pi
NP = [-0.5*sum(s[i, j]*np.outer(uvec[i], vvec[j]) for i in range(2) for j in range(2)) for s in SIGMA]
NK = [-0.5*sum(s[i, j]*np.outer(vvec[i], uvec[j]) for i in range(2) for j in range(2)) for s in SIGMA]
ETA = np.diag([1., -1., -1., -1.])
MP = [1j*N for N in NP]
MKL = [ETA[m, m]*1j*NK[m] for m in range(4)]
LOR = teilalgebra([lambda M: M @ S4 - S4 @ M, lambda M: M @ Inp + Inp @ M.T])
QP = [(MP[m] + MKL[m])/2 for m in range(4)]                         # Q^+_mu, Q^+_0 = E (z98)
SO32 = LOR + QP
NAMEN = [f"L{k}" for k in range(6)] + ["Q0=E", "Q1", "Q2", "Q3"]
BASIS_R = np.array([np.concatenate([M.real.ravel(), M.imag.ravel()]) for M in SO32]).T


def koord(M):
    """Reelle Koordinaten von M in der Basis SO32 (Lorentz, Q^+)."""
    x, *_ = np.linalg.lstsq(BASIS_R, np.concatenate([M.real.ravel(), M.imag.ravel()]), rcond=None)
    return x


FSTR = np.array([[koord(A @ B - B @ A) for B in SO32] for A in SO32])               # [X_a, X_b] = f_ab^c X_c
res_fstr = max(np.abs(sum(FSTR[a, b, c]*SO32[c] for c in range(10)) - (SO32[a] @ SO32[b] - SO32[b] @ SO32[a])).max()
               for a in range(10) for b in range(10))


def klammer_r(x, y):
    return np.einsum("a,b,abc->c", x, y, FSTR)


# ------------------------------------------------------------------ Fockraum einer Naht in der Faktorisierung Original (x) Partner
NV = 11
q = Q_TEST
s_q = math.sqrt(1 - q*q)
kap = math.sqrt((1 - q)/(1 + q))


def moden(n, k):
    a = sps.diags(np.sqrt(np.arange(1, n)), 1, format="csr")
    Id = sps.identity(n, format="csr")
    out = []
    for j in range(k):
        mats = [Id]*k
        mats[j] = a
        op = mats[0]
        for Mx in mats[1:]:
            op = sps.kron(op, Mx, format="csr")
        out.append(op)
    return out


def dag(A):
    return sps.csr_matrix(A).T.conj().tocsr()


a1, a2, b1, b2 = moden(NV, 4)                                       # Reihenfolge (a1, a2 | b1, b2): H = H_a (x) H_b
AO, BO = [a1, a2], [b1, b2]
DIM, DA = NV**4, NV**2
OCC = np.array(list(itertools.product(range(NV), repeat=4)))
TOT = OCC.sum(1)
NULL = sps.csr_matrix((DIM, DIM), dtype=complex)
d_ = [((AO[i] - q*dag(BO[i]))/s_q).tocsr() for i in range(2)]       # d = (a - q b^+)/s, e = (b - q a^+)/s (z89 invertiert)
e_ = [((BO[i] - q*dag(AO[i]))/s_q).tocsr() for i in range(2)]
ZT = [d_[0], d_[1], dag(e_[0]), dag(e_[1])]                          # Z = (d_1, d_2, e_1^+, e_2^+)
ZTd = [dag(x) for x in ZT]
CC = [(d_[i] + dag(e_[i])).tocsr() for i in range(2)]               # c = d + e^+ = kappa (a + b^+)
ERZ = [dag(d_[0]), dag(d_[1]), dag(e_[0]), dag(e_[1])]              # Einquanten-Erzeuger relativ zu Omega


def fock_gen(M):
    """Hermitesches H mit [H, Z] = -i M Z, also e^(itH) Z e^(-itH) = e^(tM) Z: H = Z^+ (i G4 M) Z."""
    Am = 1j*G4 @ M
    H = NULL.copy()
    for al in range(4):
        for be in range(4):
            if abs(Am[al, be]) > 1e-14:
                H = H + Am[al, be]*(ZTd[al] @ ZT[be])
    return H.tocsr()


HS = [fock_gen(M) for M in SO32]
LOWc = np.where(TOT <= NV - 4)[0]
om1 = np.diag([s_q*q**n for n in range(NV)])
C_OM = np.kron(om1, om1)                                            # C[(a1 a2), (b1 b2)]
C_OM /= np.linalg.norm(C_OM)
OMEGA = C_OM.reshape(-1).astype(complex)
DOP = sps.diags((OCC[:, 0] + OCC[:, 1] - OCC[:, 2] - OCC[:, 3]).astype(float)).tocsr()   # D = N_a - N_b = N_d - N_e
PARITAET = sps.diags((-1.0)**TOT).tocsr()


def polar(C):
    U_, sv, Wh = np.linalg.svd(C)
    return U_ @ Wh, sv, U_


def modular_K(C):
    """K^a = -log rho_a, K^b = -log rho_b aus der Singulärwertzerlegung von C (rho_a = C C^+, rho_b = C^T conj(C))."""
    U_, sv, Wh = np.linalg.svd(C)
    lg = -2*np.log(np.maximum(sv, 1e-300))
    Ka = (U_*lg) @ U_.conj().T
    Kb = (Wh.T*lg) @ Wh.conj()
    return Ka, Kb


def K_auf(Ka, Kb, v):
    """Zweiseitiger Modulhamiltonoperator K = K^a (x) 1 - 1 (x) K^b auf einen Vektor (C -> K^a C - C K^b^T)."""
    C = v.reshape(DA, DA)
    return (Ka @ C - C @ Kb.T).reshape(-1)


def J_spalten(V, cols):
    """Spalten U_psi e_j der Tomita-Konjugation J_psi = (V (x) V^T) SWAP K (C = |C^+| V Polarzerlegung)."""
    out = np.zeros((DIM, len(cols)), complex)
    for k, j in enumerate(cols):
        ia, ib = divmod(j, DA)
        out[:, k] = np.kron(V[:, ib], V[ia, :])
    return out


def stabilisator(V, L0):
    """Reelle Lösungen von (X + i phi) U_psi = U_psi conj(X + i phi), X in so(3,2) (Fock), auf Spalten mit N_tot <= L0."""
    cols = np.where(TOT <= L0)[0]
    cols2 = np.where(TOT <= L0 + 2)[0]
    pos = {j: k for k, j in enumerate(cols2)}
    U2 = J_spalten(V, cols2)
    U1 = U2[:, [pos[j] for j in cols]]
    bl = []
    for H in HS:
        bl.append((H @ U1 + U2 @ H.conj()[cols2][:, cols].toarray()).ravel())
    bl.append((2*U1).ravel())
    Am = np.array(bl).T
    Ar = np.concatenate([Am.real, Am.imag])
    _, sv, vh = np.linalg.svd(Ar, full_matrices=False)
    return sv, vh


def koh(beta):
    Lg = sum((beta[j]*dag(CC[j]) - np.conj(beta[j])*CC[j] for j in range(2)), NULL).tocsc()
    return Lg


def erw(v, X):
    return complex(np.vdot(v, X @ v))


laufzeit()

# =================================================================== Teil 1
say("Teil 1: Gravitationsseite -- Tomita-Familie, Cartan-Verbindung, MM-Krümmung, quadratische Antwort")

# --- 1.1 Modulare Daten auf Abschneidungen
say("  1.1 Modulare Daten (Original a, Partner b, NV = 11 je Mode, q = 0,25 Prüfwert)")
gen_err = 0.0
for H, M in zip(HS, SO32):
    for al in range(4):
        lhs = (H @ ZT[al] - ZT[al] @ H)[:, LOWc]
        rhs = sum(((-1j*M[al, be])*ZT[be] for be in range(4)), NULL)[:, LOWc]
        gen_err = max(gen_err, abs(lhs - rhs).max())
herm_err = max(abs(H - dag(H)).max() for H in HS)
vak_err = max(np.abs(x @ OMEGA).max() for x in d_ + e_)
V0 = polar(C_OM)[0]
Ka0, Kb0 = modular_K(C_OM)
lowa = np.where(np.array(list(itertools.product(range(NV), repeat=2))).sum(1) <= 6)[0]
Na_diag = np.array(list(itertools.product(range(NV), repeat=2))).sum(1).astype(float)
k0_soll = 2*math.log(1/q)*Na_diag - 2*math.log(s_q**2)
ok_K0 = np.abs((Ka0 - np.diag(k0_soll))[np.ix_(lowa, lowa)]).max() < 1e-8 and \
    np.abs((Kb0 - np.diag(k0_soll))[np.ix_(lowa, lowa)]).max() < 1e-8
D_kom = max(abs((DOP @ H - H @ DOP)[:, LOWc]).max() for H in HS)
# (a1) Impulsphasen-Kondensat
BETA = np.array([0.03, 0.02j])
aa2 = moden(NV, 2)


def w_a(beta, fak=1.0):
    g = sum((fak*kap*beta[j]*dag(aa2[j]) - fak*kap*np.conj(beta[j])*aa2[j] for j in range(2)), sps.csr_matrix((DA, DA)))
    return sla.expm(g.toarray())


L_B = koh(BETA)
PSI_B = spla.expm_multiply(L_B, OMEGA)
C_B = PSI_B.reshape(DA, DA)
V_B = polar(C_B)[0]
wa1, wa2 = w_a(BETA), w_a(BETA, 2.0)
ok_polar = np.abs((V_B - wa2)[np.ix_(lowa, lowa)]).max() < 1e-10
rho_B = C_B @ C_B.conj().T
ok_rho = np.abs((rho_B - wa1 @ (C_OM @ C_OM.T) @ wa1.conj().T)[np.ix_(lowa, lowa)]).max() < 1e-12
Ka_B, Kb_B = modular_K(C_B)
ok_KB = np.abs((Ka_B - wa1 @ Ka0 @ wa1.conj().T)[np.ix_(lowa[:15], lowa[:15])]).max() < 1e-6 and \
    np.abs((Kb_B - wa1.T @ Kb0 @ wa1.conj())[np.ix_(lowa[:15], lowa[:15])]).max() < 1e-6          # w_b = w_a^T
dK = Ka_B - Ka0                                                      # affin: [dK, a_j] ist ein Vielfaches der Eins
aa2d = [x.toarray() for x in aa2]
lo2 = np.where(np.array(list(itertools.product(range(NV), repeat=2))).sum(1) <= 4)[0]
aff_res = 0.0
for A in aa2d + [x.T for x in aa2d]:
    Kc = (dK @ A - A @ dK)[np.ix_(lo2, lo2)]
    aff_res = max(aff_res, np.abs(Kc - Kc[0, 0]*np.eye(len(lo2))).max())
ZAHLEN["modulare_daten"] = {"NV je Mode": NV, "q (Prüfwert)": q, "q nativ": Q_NATIV, "Fehler [H, Z] = -i M Z": gen_err,
                            "|d Omega|, |e Omega|": vak_err, "max |[D, H_so(3,2)]|": D_kom,
                            "Kondensat beta": [str(x) for x in BETA],
                            "Polarteil V gegen w_a(2 beta) (tiefer Block)": float(np.abs((V_B - wa2)[np.ix_(lowa, lowa)]).max()),
                            "K^a_psi - K^a: Residuum der Affinität": aff_res}
check(f"[E] Modulare Daten auf Abschneidungen: Die nativen Erzeuger (Lorentz, Q^+_mu, Q_0 = E) sind als Fockoperatoren "
      f"H = Z^+ (i G4 M) Z mit d = (a - q b^+)/s, e = (b - q a^+)/s exakt ([H, Z] = -i M Z, Fehler {gen_err:.1e}, hermitesch "
      f"{herm_err:.0e}), Omega (Purifikation, KMS) wird von d, e vernichtet ({vak_err:.0e}), J_Omega = SWAP K (Polarteil 1), "
      f"K^a = -log rho_a = 2 log(1/q) N_a + const, und D = N_a - N_b vertauscht mit so(3,2) ({D_kom:.0e}): der Modulgenerator "
      f"ist auf den nativen Erzeugern zentral (z98). Kondensat (a1) psi_beta = exp(beta.c^+ - h.c.) Omega: Polarteil V = "
      f"w_a(2 beta) = exp(2 kappa (beta.a^+ - h.c.)), also J_psi = (w_a^2 (x) 1) J (w_a^2 (x) 1)^+; rho_a(psi) = w_a rho_a w_a^+, "
      f"K_psi = W K W^+ (W = w_a (x) w_b die volle Verschiebung); K^a_psi - K^a ist affin (Kommutatoren mit a, a^+ sind "
      f"Vielfache der Eins, Residuum {aff_res:.1e}): Das Kondensat ist quasifrei",
      gen_err < 1e-12 and herm_err < 1e-12 and vak_err < 1e-14 and np.abs(V0 - np.eye(DA)).max() < 1e-14 and ok_K0
      and D_kom < 1e-12 and ok_polar and ok_rho and ok_KB and aff_res < 1e-8)

# (a2) wörtliches D = -4-Kondensat
H_GRAV = (dag(e_[0]) @ dag(e_[0]) @ dag(e_[1]) @ dag(e_[1]) @ OMEGA)/2                 # |0, 0; 2, 2>, D = -4, E = 3
H_GRAV40 = (dag(e_[0]) @ dag(e_[0]) @ dag(e_[0]) @ dag(e_[0]) @ OMEGA)/math.sqrt(24)    # |0, 0; 4, 0>


def kond4(lam, h=H_GRAV):
    v = OMEGA + lam*h
    return v/np.linalg.norm(v)


LAM = 0.01                                                          # Hauptwert; die Amplitude ist frei
PSI_L = kond4(LAM)
C_L = PSI_L.reshape(DA, DA)
B_part = (dag(BO[0]) @ dag(BO[0]) @ dag(BO[1]) @ dag(BO[1]))*s_q**4/2               # h = B Omega, B im Partner
ok_partner = np.linalg.norm(H_GRAV - B_part @ OMEGA) < 1e-5 and abs(np.linalg.norm(H_GRAV) - 1) < 1e-5
ladung_h = float(np.vdot(H_GRAV, DOP @ H_GRAV).real)
Ka_L, Kb_L = modular_K(C_L)
Pa = np.diag((-1.0)**np.array(list(itertools.product(range(NV), repeat=2))).sum(1))
par_K = max(np.abs((Pa @ Ka_L - Ka_L @ Pa)[np.ix_(lowa, lowa)]).max(), np.abs((Pa @ Kb_L - Kb_L @ Pa)[np.ix_(lowa, lowa)]).max())
# Ladung-4-Anteil von K^a_psi, K^b_psi (bezüglich N_a bzw. N_b): quadratische Operatoren tragen nur Ladung 0, +-2
Nadiag = Na_diag
msk4 = np.abs(Nadiag[lowa][:, None] - Nadiag[lowa][None, :]) == 4
lad4 = {}
for lam in (2e-4, 1e-4, 0.1):
    Ka_, Kb_ = modular_K(kond4(lam).reshape(DA, DA))
    lad4[lam] = (float(np.abs(Ka_[np.ix_(lowa, lowa)][msk4]).max()), float(np.abs(Kb_[np.ix_(lowa, lowa)][msk4]).max()))
lin_a = lad4[2e-4][0]/lad4[1e-4][0]
lin_b = lad4[2e-4][1]/lad4[1e-4][1]
ok_ng = lad4[1e-4][0] > 1e-7 and abs(lin_a - 2) < 1e-3 and abs(lin_b - 2) < 1e-2 and lad4[0.1][1] > 1
ZAHLEN["d4_kondensat"] = {"lambda (Hauptwert)": LAM, "<h|D|h>": ladung_h, "h = B Omega (Partner)": bool(ok_partner),
                          "[K^a, P_a], [K^b, P_b] (tiefer Block)": float(par_K),
                          "max |Ladung-4-Anteil| (K^a, K^b) je lambda": {str(k): v for k, v in lad4.items()},
                          "Verhältnis lambda 2e-4 / 1e-4 (K^a, K^b)": [lin_a, lin_b],
                          "Störungsbereich": "lambda << q^4 (Prüfwert 0,0039; nativ q^4 = e^(-8 pi) = 1,2e-11)"}
check(f"[E] (a2) Das wörtliche D = -4-Kondensat psi = Omega + lambda |0, 0; 2, 2> (<h|D|h> = {ladung_h:.6f}) ist "
      f"psi = (1 (x) (1 + lambda B)) Omega mit B = s^4 (b_1^+)^2 (b_2^+)^2/2 im Partner (das tiefste Graviton ist eine "
      f"Partneranregung). Sein K_psi ist nicht gaußsch: K^a_psi und K^b_psi haben Matrixelemente der Ladung 4 (K^a "
      f"{lad4[1e-4][0]:.3e}, K^b {lad4[1e-4][1]:.3e} bei lambda = 1e-4, Verhältnis zu 2e-4: {lin_a:.4f}, {lin_b:.4f}), "
      f"die kein quadratischer Operator trägt (Ladungen 0, +-2); linear nur für lambda << q^4, darüber nicht störungsartig "
      f"(K^b-Anteil {lad4[0.1][1]:.2f} bei lambda = 0,1; nativ ist q^4 = 1,2e-11). K_psi vertauscht mit (-1)^N_a bzw. "
      f"(-1)^N_b (Residuum {par_K:.0e}), denn psi hat nur gerade Modulladungen (0, -4)",
      ok_partner and abs(ladung_h + 4) < 1e-6 and par_K < 1e-10 and ok_ng)
laufzeit()

# --- 1.2 Nullmoden der Tomita-Familie J_g = U(g) J_psi U(g)^-1
say("  1.2 Nullmoden (Stabilisator von J_psi in so(3,2)) und Geometrie der Bahn")
L0 = NV - 7
sv_vak, _ = stabilisator(V0, L0)
sv_koh, vh_koh = stabilisator(V_B, L0)
sv_koh2, _ = stabilisator(polar(spla.expm_multiply(koh(2*BETA), OMEGA).reshape(DA, DA))[0], L0)
lin_brech = float(sv_koh2[-3]/sv_koh[-3])
V_L = polar(C_L)[0]
sv_d4, vh_d4 = stabilisator(V_L, L0)


def nulldim(sv, rel=1e-6):
    return int(np.sum(sv < rel*sv[0]))


dims = {"Vakuum": nulldim(sv_vak), "(a1) eine Naht": nulldim(sv_koh), "(a2) D = -4": nulldim(sv_d4)}
# Twistor-Kleingruppen: [H_M, beta.c^+ - h.c.] = 0  <=>  M^T v = 0 mit v = sum conj(beta_j) v_j


def kleingruppe(betas):
    rows = []
    for M in LOR:
        parts = []
        for bt in betas:
            w = M.T @ (np.conj(bt[0])*vvec[0] + np.conj(bt[1])*vvec[1])
            parts += [w.real, w.imag]
        rows.append(np.concatenate(parts))
    return nullraum(np.array(rows).T)


kg1 = kleingruppe([BETA])
Qa = np.linalg.qr(vh_koh[-2:, :6].T)[0]
Qb = np.linalg.qr(kg1.T)[0]
winkel = np.linalg.svd(Qa.T @ Qb, compute_uv=False)
beta2_unab = np.array([0.01 - 0.02j, 0.025])
beta2_koll = (0.7 - 0.4j)*BETA
det_unab = BETA[0]*beta2_unab[1] - BETA[1]*beta2_unab[0]
zwei = {"det != 0": len(kleingruppe([BETA, beta2_unab])), "det = 0 (kollinear)": len(kleingruppe([BETA, beta2_koll]))}
# Stabilisator der D = -4-Kondensats in Lorentz-Koordinaten
st_d4 = vh_d4[-dims["(a2) D = -4"]:] if dims["(a2) D = -4"] else np.zeros((0, 11))
# Kompositzustand: N Nähte, alle im KMS-Zustand, diagonale Wirkung auf T^N mit tau^(+N)


def stab_produkt(N):
    rows = []
    for M in SO32:
        Mn = sla.block_diag(*([M]*N))
        Sn = sla.block_diag(*([S4]*N))
        C_ = Mn @ Sn - Sn @ Mn
        rows.append(np.concatenate([C_.real.ravel(), C_.imag.ravel()]))
    return len(nullraum(np.array(rows).T))


komp = {N: stab_produkt(N) for N in (1, 2, 3)}
ZAHLEN["nullmoden"] = {"kleinste Singulärwerte (Vakuum)": [float(x) for x in sv_vak[-7:]],
                       "kleinste Singulärwerte (a1)": [float(x) for x in sv_koh[-4:]],
                       "(a1) Bruch der übrigen Richtungen, Faktor bei 2 beta": lin_brech,
                       "kleinste Singulärwerte (a2)": [float(x) for x in sv_d4[-4:]], "Dimensionen": dims,
                       "cos Hauptwinkel Fock-Nullraum (a1) gegen Twistor-Kleingruppe": [float(x) for x in winkel],
                       "zwei Nähte (a1)": zwei, "det(beta1, beta2)": str(np.round(det_unab, 6)),
                       "Kompositzustand N = 1, 2, 3": komp,
                       "Stabilisator (a2) in Koordinaten": np.round(st_d4, 6).tolist()}
check(f"[E] Nullmoden der Tomita-Familie (Fock, NV = 11, Spalten N_tot <= {L0}): Vakuum {dims['Vakuum']} = Lorentz "
      f"(z98), Lücke {sv_vak[-7]:.3f} gegen {sv_vak[-6]:.0e}; Kondensat (a1) auf einer Naht {dims['(a1) eine Naht']}, und "
      f"das sind genau die Nullrotationen, die den Spinor beta in der c-Ebene festhalten (cos der Hauptwinkel "
      f"{np.round(winkel, 9).tolist()}; Twistor: M^T (conj(beta).v) = 0); die übrigen vier Lorentzrichtungen brechen linear "
      f"in beta (Singulärwerte {sv_koh[-3]:.3f}, {sv_koh[-4]:.3f}, bei 2 beta Faktor {lin_brech:.4f}). Zwei Nähte "
      f"(diagonale Wirkung) mit det(beta1, beta2) != 0: Schnitt der "
      f"Kleingruppen {zwei['det != 0']}, kollinear (det = 0): {zwei['det = 0 (kollinear)']}. Das Kondensat, das I_pi "
      f"wählt, löscht also alle Nullmoden. D = -4-Kondensat (a2): {dims['(a2) D = -4']} Nullmode(n). Kompositzustand aus "
      f"N Nähten im KMS-Zustand (J = J^(x N), diagonale Wirkung): {list(komp.values())} = Lorentz für N = 1, 2, 3",
      dims["Vakuum"] == 6 and dims["(a1) eine Naht"] == 2 and len(kg1) == 2 and winkel.min() > 1 - 1e-6
      and zwei["det != 0"] == 0 and zwei["det = 0 (kollinear)"] == 2 and all(v == 6 for v in komp.values())
      and sv_vak[-7] > 1e-2 and sv_koh[-3] > 1e-2 and abs(lin_brech - 2) < 1e-2)


# Reduktivität und Symmetrie von (so(3,2), h)
def reduktiv(hb):
    """Gibt es ein ad(h)-invariantes Komplement m? Löst P_h[h_a, Y] + [h_a, phi(Y)] = phi(P_m0[h_a, Y]) linear in phi."""
    hb = np.atleast_2d(hb)
    Hq = np.linalg.qr(hb.T)[0]
    m0 = nullraum(Hq.T)
    Mq = np.linalg.qr(m0.T)[0]
    k, mdim = Hq.shape[1], Mq.shape[1]
    Ph = lambda x: Hq.T @ x                                            # Koordinaten in h
    Pm = lambda x: Mq.T @ x
    rows, rhs = [], []
    for a in range(k):
        ha = Hq[:, a]
        for i in range(mdim):
            Y = Mq[:, i]
            br = klammer_r(ha, Y)
            const = Ph(br)
            pm = Pm(br)
            # Unbekannte phi[p, i'] (p in h, i' in m0): [h_a, phi(Y_i)] - phi(P_m0[h_a, Y_i])
            coeff = np.zeros((k, k, mdim))
            for p in range(k):
                hp = Hq[:, p]
                coeff[:, p, i] += Ph(klammer_r(ha, hp))
                for i2 in range(mdim):
                    coeff[p, p, i2] -= pm[i2]
            for c in range(k):
                rows.append(coeff[c].ravel())
                rhs.append(-const[c])
    A_, b_ = np.array(rows), np.array(rhs)
    sol, *_ = np.linalg.lstsq(A_, b_, rcond=None)
    res = float(np.linalg.norm(A_ @ sol - b_))
    phi = sol.reshape(k, mdim)
    m_basis = Mq + Hq @ phi
    return res, Hq, m_basis


def split(hb, mb, x):
    B = np.concatenate([hb, mb], axis=1)
    c = np.linalg.solve(B, x)
    return hb @ c[:hb.shape[1]], mb @ c[hb.shape[1]:]


res_vak, H_vak, m_vak = reduktiv(np.array([koord(M) for M in LOR]))
res_koh, _, _ = reduktiv(np.array([koord(sum(c*M for c, M in zip(v, LOR))) for v in kg1]))
sym_vak = max(np.linalg.norm(split(H_vak, m_vak, klammer_r(m_vak[:, i], m_vak[:, j]))[1])
              for i in range(4) for j in range(4))
st_d4_lor = st_d4[:, :10] if len(st_d4) else np.zeros((0, 10))
red_d4 = None
if len(st_d4_lor):
    res_d4, H_d4, m_d4 = reduktiv(st_d4_lor)
    tors_d4 = max(np.linalg.norm(split(H_d4, m_d4, klammer_r(m_d4[:, i], m_d4[:, j]))[1])
                  for i in range(m_d4.shape[1]) for j in range(m_d4.shape[1]))
    M_st = sum(c*M for c, M in zip(st_d4_lor[0], SO32))
    dreh = bool(np.abs(M_st[:2, 2:]).max() < 1e-6 and np.abs(M_st[2:, :2]).max() < 1e-6
                and np.abs(np.linalg.eigvals(M_st).real).max() < 1e-6)
    red_d4 = {"Residuum": res_d4, "max |P_m [m, m]| (Torsion)": tors_d4, "dim h": int(H_d4.shape[1]),
              "Stabilisator ist eine Drehung (blockdiagonal, kompakt)": dreh}
# Kleingruppe nilpotent? (ad-Spektrum)
kg_r = np.array([koord(sum(c*M for c, M in zip(v, LOR))) for v in kg1])
adk = [np.array([klammer_r(x, np.eye(10)[b]) for b in range(10)]).T for x in kg_r]
nilp = max(np.abs(np.linalg.matrix_power(A, 5)).max() for A in adk)
ZAHLEN["reduktivitaet"] = {"Vakuum (Lorentz) Residuum": res_vak, "Vakuum max |P_m [m, m]|": sym_vak,
                           "(a1) Kleingruppe Residuum": res_koh, "(a1) ad^5 der Kleingruppe": float(nilp),
                           "(a2)": red_d4, "Strukturkonstanten-Residuum": float(res_fstr)}
text_d4 = (f"Stabilisator {red_d4['dim h']}-dimensional, eine Drehung (Achse des Gravitons |0, 0; 2, 2>), reduktiv "
           f"(Residuum {red_d4['Residuum']:.0e}), nicht symmetrisch (Torsion {red_d4['max |P_m [m, m]| (Torsion)']:.3f})"
           if red_d4 else "keine Nullmoden")
check(f"[E] Geometrie der Bahnen: Vakuum: (so(3,2), Lorentz) reduktiv (Residuum {res_vak:.0e}) und symmetrisch "
      f"(|P_m[m, m]| = {sym_vak:.0e}): AdS4, torsionsfrei, z98. Kondensat (a1): die Kleingruppe ist nilpotent (ad^5 = "
      f"{nilp:.0e}) und es gibt KEIN ad(h)-invariantes Komplement (kleinstes Residuum {res_koh:.3f}): die achtdimensionale "
      f"Bahn ist nicht reduktiv, eine Zerlegung in Verbindung und Vierbein existiert nicht. Zwei Nähte mit det != 0: keine "
      f"Nullmoden, die Familie ist die Gruppe SO(3,2) selbst. D = -4-Kondensat (a2): {text_d4}",
      res_fstr < 1e-10 and res_vak < 1e-10 and sym_vak < 1e-10 and res_koh > 1e-3 and nilp < 1e-10
      and red_d4 is not None and red_d4["Residuum"] < 1e-8 and red_d4["max |P_m [m, m]| (Torsion)"] > 1e-3
      and red_d4["Stabilisator ist eine Drehung (blockdiagonal, kompakt)"])
laufzeit()

# --- 1.3 Cartan-Verbindung und MM-Krümmung
say("  1.3 Cartan-Verbindung, MM-Krümmung, modulare Berry-Krümmung")
XB, YB = QP[1] + 0.3*QP[0] + 0.2*LOR[2], QP[2] - 0.5*QP[3] + 0.4*LOR[4]


def g_st(s_, t_):
    return sla.expm(s_*XB) @ sla.expm(t_*YB)


s0, t0, hh = 0.31, -0.27, 1e-5
A_s = lambda s_, t_: np.linalg.inv(g_st(s_, t_)) @ (g_st(s_ + hh, t_) - g_st(s_ - hh, t_))/(2*hh)
A_t = lambda s_, t_: np.linalg.inv(g_st(s_, t_)) @ (g_st(s_, t_ + hh) - g_st(s_, t_ - hh))/(2*hh)
dsAt = (A_t(s0 + hh, t0) - A_t(s0 - hh, t0))/(2*hh)
dtAs = (A_s(s0, t0 + hh) - A_s(s0, t0 - hh))/(2*hh)
As0, At0 = A_s(s0, t0), A_t(s0, t0)
F_MC = dsAt - dtAs + (As0 @ At0 - At0 @ As0)
# Torsion/Krümmung der kanonischen Verbindung je Zustand (reduktive Fälle)
xs, ys = koord(As0), koord(At0)


br_vak = klammer_r(split(H_vak, m_vak, xs)[1], split(H_vak, m_vak, ys)[1])
R_vak, T_vak = split(H_vak, m_vak, -br_vak)
tors_split = None
if red_d4:
    br_d4 = klammer_r(split(H_d4, m_d4, xs)[1], split(H_d4, m_d4, ys)[1])
    R_d4, T_d4 = split(H_d4, m_d4, -br_d4)
    tors_split = {"|R|": float(np.linalg.norm(R_d4)), "|T|": float(np.linalg.norm(T_d4))}
# modulare Berry-Krümmung mit den Nullmoden von K_psi (a1): F(X, Y) = -[l_X, l_Y] = 2i Im(r_X^T G4 conj(r_Y)), r_X = M^T v


def r_vek(M, beta):
    return M.T @ (np.conj(beta[0])*vvec[0] + np.conj(beta[1])*vvec[1])


def berry_zentral(beta):
    R_ = [r_vek(M, beta) for M in SO32]
    return np.array([[2*np.imag(R_[i] @ G4 @ np.conj(R_[j])) for j in range(10)] for i in range(10)])


FB = berry_zentral(BETA)
FB2 = berry_zentral(2*BETA)
# Fock-Kontrolle: <Omega| [l_X, l_Y] |Omega> mit l_X = [iH_X, L_beta] für zwei Lorentz- und zwei Q-Richtungen
LB = L_B
idx_k = (1, 4, 6, 7)
lX = [(1j*(HS[k] @ LB - LB @ HS[k])).tocsr() for k in idx_k]
fock_FB = [[complex(np.vdot(OMEGA, (lX[i] @ lX[j] - lX[j] @ lX[i]) @ OMEGA)) for j in range(4)] for i in range(4)]
ok_fockFB = all(abs(-fock_FB[i][j] - 1j*FB[idx_k[i], idx_k[j]]) < 1e-9 for i in range(4) for j in range(4)) \
    and max(abs(FB[idx_k[i], idx_k[j]]) for i in range(2) for j in range(2, 4)) > 1e-5
zentral_res = float(np.abs((lX[0] @ lX[2] - lX[2] @ lX[0]) @ OMEGA - fock_FB[0][2]*OMEGA).max())
rang_FB = int(np.linalg.matrix_rank(FB, tol=1e-12))
blk_LL, blk_QQ, blk_LQ = float(np.abs(FB[:6, :6]).max()), float(np.abs(FB[6:, 6:]).max()), float(np.abs(FB[:6, 6:]).max())
kg_kern = max(float(np.abs(FB @ koord(sum(c*M for c, M in zip(v, LOR)))).max()) for v in kg1)
ZAHLEN["kruemmung"] = {"|F_MM| (Maurer-Cartan, Schleife)": float(np.abs(F_MC).max()),
                       "Vakuum: |R|, |T| der kanonischen Verbindung": [float(np.linalg.norm(R_vak)), float(np.linalg.norm(T_vak))],
                       "(a2): kanonische Verbindung": tors_split, "zentrale Berry-Krümmung (a1), Rang": rang_FB,
                       "zentrale Berry-Krümmung max |F| Lorentz-Lorentz, Q-Q, Lorentz-Q": [blk_LL, blk_QQ, blk_LQ],
                       "zentrale Berry-Krümmung auf der Kleingruppe": kg_kern,
                       "Skalierung beta -> 2 beta": float(np.abs(FB2).max()/np.abs(FB).max()),
                       "Fock-Kontrolle (zentral, Residuum)": zentral_res}
text_t = f"Torsion |T| = {tors_split['|T|']:.3f} != 0" if tors_split else "ohne Zerlegung"
check(f"[E] Cartan-Verbindung und Krümmung: Für jede Familie J_g = U(g) J_psi U(g)^-1 ist die Cartan-Verbindung mit "
      f"Rahmen in der nativen Gruppe die Maurer-Cartan-Form A = g^-1 dg und hängt nicht von psi ab: F_MM = dA + A ^ A = "
      f"{np.abs(F_MC).max():.1e} (Schleife e^(sX) e^(tY)), für Gauß-, Kondensat- und Kompositzustände gleich. psi wählt nur "
      f"die Zerlegung A = Gamma + e: Vakuum R = -[e, e]_h (|R| = {np.linalg.norm(R_vak):.3f}), Torsion "
      f"{np.linalg.norm(T_vak):.0e}; (a2) {text_t}; "
      f"(a1) keine invariante Zerlegung. Modulare Berry-Krümmung mit den Nullmoden von K_psi (Czech-de-Boer-Vorschrift): "
      f"Vakuum null (D zentral); Kondensat (a1) zentral, F(X, Y) = -[l_X, l_Y] = 2i Im(r_X^T G4 conj(r_Y)) mit r_X = "
      f"M_X^T (conj(beta).v) (Fock-Kontrolle {ok_fockFB}, Residuum der Zentralität {zentral_res:.0e}, Abschneidefehler), "
      f"Rang {rang_FB}; sie verschwindet auf Lorentz x Lorentz ({blk_LL:.0e}) und auf den AdS-Translationen Q x Q "
      f"({blk_QQ:.0e}), weil beide Bilder r_X in total G4-nullen Ebenen liegen, koppelt nur Lorentz an Q (max "
      f"{blk_LQ:.2e}) und verschwindet auf der Kleingruppe ({kg_kern:.0e}); "
      f"quadratisch in beta (beta -> 2 beta: Faktor {np.abs(FB2).max()/np.abs(FB).max():.6f}): eine U(1)-Krümmung (die "
      f"symplektische Fläche der Verschiebung), keine Lorentzkrümmung",
      np.abs(F_MC).max() < 1e-6 and np.linalg.norm(T_vak) < 1e-8 and np.linalg.norm(R_vak) > 1e-3 and ok_fockFB
      and zentral_res < 1e-6 and abs(np.abs(FB2).max()/np.abs(FB).max() - 4) < 1e-9 and rang_FB == 6
      and blk_LL < 1e-14 and blk_QQ < 1e-14 and blk_LQ > 1e-5 and kg_kern < 1e-14
      and (tors_split is None or tors_split["|T|"] > 1e-4))

# --- 1.4 Quadratische Antwort: Fisher-Information der Familie U(g) psi
say("  1.4 Quadratische Antwort (Fisher-Information 4 Cov_psi(H_i, H_j))")


def fisher(v):
    Hv = [H @ v for H in HS]
    mu = np.array([np.vdot(v, x) for x in Hv])
    Cm = np.array([[np.vdot(Hv[i], Hv[j]) for j in range(10)] for i in range(10)]) - np.outer(mu.conj(), mu)
    return 4*Cm.real


F_vak = fisher(OMEGA)
F_k1 = fisher(PSI_B)
F_k2 = fisher(spla.expm_multiply(koh(2*BETA), OMEGA))
ex1, ex2 = F_k1 - F_vak, F_k2 - F_vak
skal_k = float(np.abs(ex2).max()/np.abs(ex1).max())
F_l1 = fisher(kond4(0.02))
F_l2 = fisher(kond4(0.01))
skal_l = float(np.abs(F_l1 - F_vak).max()/np.abs(F_l2 - F_vak).max())
ew_vak = np.linalg.eigvalsh(F_vak)
qblock = np.round(F_vak[6:, 6:], 10)
ZAHLEN["fisher"] = {"Vakuum: Eigenwerte": np.round(ew_vak, 10).tolist(), "Vakuum: Q-Block": qblock.tolist(),
                    "Vakuum: Rang": int(np.sum(ew_vak > 1e-9)), "(a1) Überschuss, Skalierung beta -> 2 beta": skal_k,
                    "(a1) max Überschuss": float(np.abs(ex1).max()),
                    "(a2) Überschuss, Skalierung lambda 0,02 / 0,01": skal_l,
                    "Kompositzustand": "additiv: F_N = N F_1 (Produktzustand)"}
check(f"[E] Quadratische Antwort: Die Fisher-Information 4 Cov_psi(H_i, H_j) der Familie U(g) psi ist im Vakuum fest und "
      f"q-unabhängig (Rang {int(np.sum(ew_vak > 1e-9))}: E und Drehungen lassen Omega fest; Q-Block "
      f"diag{np.round(np.diag(qblock), 6).tolist()}, also verschwindet die Antwort gerade in der AdS-Zeitrichtung Q_0 = E); "
      f"beim Kondensat (a1) kommt ein Überschuss 4 ||l_X Omega||^2 hinzu, der quadratisch in beta skaliert (beta -> 2 beta: "
      f"{skal_k:.6f}), beim D = -4-Kondensat (a2) einer quadratisch in lambda ({skal_l:.4f} für 0,02 gegen 0,01); für N "
      f"Nähte im Produktzustand ist sie additiv (N F_1). Fest ist die quadratische Antwort damit nur im Vakuum, wo F = 0 "
      f"nichts zu normieren hat",
      int(np.sum(ew_vak > 1e-9)) == 6 and abs(qblock[0, 0]) < 1e-9 and abs(skal_k - 4) < 1e-3 and abs(skal_l - 4) < 0.01)
befund("E (negativ für H1, Gravitationsseite): Die Tomita-Familie eines beliebigen Zustands unter der nativen Gruppe ist eine "
       "Bahn; ihre Cartan-Verbindung ist die Maurer-Cartan-Form und F_MM = 0, gaußsch oder nicht. Der Zustand ändert nur die "
       "Nullmoden: Das Kondensat (a1) bricht Lorentz auf die Nullrotationen von beta (nicht reduktiv, keine "
       "Verbindung-Vierbein-Zerlegung), zwei Nähte mit det(beta1, beta2) != 0 löschen alle Nullmoden; das D = -4-Kondensat "
       "gibt eine reduktive, nicht symmetrische Bahn mit Torsion; der Kompositzustand gibt AdS4 unverändert. Eine "
       "gekrümmte vierdimensionale Cartan-Geometrie entsteht nie. Die einzige neue Krümmung ist die zentrale "
       "U(1)-Berry-Krümmung des Kondensats, quadratisch in der freien Amplitude beta. b wird nicht festgelegt: Es "
       "multipliziert F ^ F = 0, und die quadratische Antwort ist nur im Vakuum fest", True)
laufzeit()

# =================================================================== Teil 2
say("\nTeil 2: Wechselwirkungsseite -- verbundene 1->2-Amplituden (Kriterium z99)")
# Zweiquanten- und Dreiquantenräume relativ zu Omega
mon2 = [ERZ[i] @ (ERZ[j] @ OMEGA) for i in range(4) for j in range(i, 4)]
mon3 = [ERZ[i] @ (ERZ[j] @ (ERZ[k] @ OMEGA)) for i in range(4) for j in range(i, 4) for k in range(j, 4)]
Q2 = np.linalg.qr(np.array(mon2).T)[0]
Q3 = np.linalg.qr(np.array(mon3).T)[0]
pi2 = lambda v: float(np.linalg.norm(Q2.conj().T @ v))
pi3 = lambda v: float(np.linalg.norm(Q3.conj().T @ v))
LIN = [d_[0], d_[1], e_[0], e_[1]] + ERZ
QUAD = [(X @ Y).tocsr() for X in LIN for Y in LIN]                  # alle 64 quadratischen Monome

# --- 2.1 Kondensat (a1) im eigenen Quasiteilchenrahmen: W^+ [X, W b^+ W^+] psi = [W^+ X W, b^+] Omega
say("  2.1 Impulsphasen-Kondensat (a1)")
Wv = lambda v: spla.expm_multiply(L_B, v)
Wdv = lambda v: spla.expm_multiply(-L_B, v)
Ka_Bs, Kb_Bs = modular_K(C_B)
KB_op = lambda v: K_auf(Ka_Bs, Kb_Bs, v)
Wb = [Wv(b @ OMEGA) for b in ERZ]                                    # W b^+ Omega = b_psi^+ psi


def verbunden_koh(Xapp):
    """max ||Pi_2 W^+ [X, b_psi^+] psi|| über b^+ in ERZ, X als Abbildung v -> X v."""
    xpsi = Xapp(PSI_B)
    m = 0.0
    for b, wb in zip(ERZ, Wb):
        v = Wdv(Xapp(wb)) - b @ Wdv(xpsi)
        m = max(m, pi2(v))
    return m


v_quad = max(verbunden_koh(lambda v, X=X: X @ v) for X in QUAD)
v_nat = max(verbunden_koh(lambda v, H=H: H @ v) for H in HS)
v_K = verbunden_koh(KB_op)
X3 = (ERZ[0] @ ERZ[2] @ d_[1]).tocsr()                              # kubische Kontrolle d_1^+ e_1^+ d_2
v_3 = verbunden_koh(lambda v: X3 @ v)
ZAHLEN["verbunden_a1"] = {"64 quadratische Monome": v_quad, "10 native Erzeuger": v_nat, "K_psi (aus dem Zustand)": v_K,
                          "Kontrolle d1^+ e1^+ d2 (gewählte Kopplung 1)": v_3}
check(f"[E] (a1) Verbundene 1->2-Amplituden im eigenen Rahmen des Kondensats (b_psi^+ = W b^+ W^+, Zweiquantenanteil von "
      f"W^+ [X, b_psi^+] psi): alle 64 quadratischen Monome in d, d^+, e, e^+ max {v_quad:.1e}, die zehn nativen Erzeuger "
      f"{v_nat:.1e}, der aus dem Zustand berechnete Modulhamiltonoperator K_psi {v_K:.1e}; die kubische Kontrolle "
      f"d_1^+ e_1^+ d_2 gibt {v_3:.6f} (die gewählte Kopplung). W^+ X W = X + linear + const: Das Kondensat verschiebt "
      f"nur, es koppelt nicht", v_quad < 1e-8 and v_nat < 1e-8 and v_K < 1e-6 and abs(v_3 - 1) < 1e-3)
laufzeit()

# --- 2.2 D = -4-Kondensat (a2): Paritätsauswahlregel und quartischer Vertex
say("  2.2 D = -4-Kondensat (a2)")
par_nat = max(abs((PARITAET @ H - H @ PARITAET)[:, LOWc]).max() for H in HS)
par_psi = float(np.linalg.norm(PARITAET @ PSI_L - PSI_L))
v2_L, v3_L = {}, {}
for lam in (1e-4, 2e-4, 0.01, 0.1):
    C_ = kond4(lam).reshape(DA, DA)
    Ka_, Kb_ = modular_K(C_)
    Kv = lambda v: K_auf(Ka_, Kb_, v)
    KO = Kv(OMEGA)
    v2_L[lam] = max(pi2(Kv(b @ OMEGA) - b @ KO) for b in ERZ)
    v3_L[lam] = max(pi3(Kv(b @ OMEGA) - b @ KO) for b in ERZ)
# im Zustand: <psi| b_i b_j X b_k^+ |psi> für X = K_psi und die nativen Erzeuger (lambda = 0,01)
amp_state = 0.0
Kv = lambda v: K_auf(Ka_L, Kb_L, v)
for X in [Kv] + [lambda v, H=H: H @ v for H in HS]:
    for bk in ERZ:
        u = X(bk @ PSI_L)
        for bi in ERZ:
            for bj in ERZ:
                amp_state = max(amp_state, abs(np.vdot(bi @ (bj @ PSI_L), u)))
lin3 = v3_L[2e-4]/v3_L[1e-4]
ZAHLEN["verbunden_a2"] = {"[P, H_so(3,2)]": float(par_nat), "|P psi - psi|": par_psi,
                          "Pi_2 [K_psi, b^+] Omega je lambda": {str(k): v for k, v in v2_L.items()},
                          "Pi_3 [K_psi, b^+] Omega je lambda": {str(k): v for k, v in v3_L.items()},
                          "Verhältnis Pi_3 (2e-4 / 1e-4)": lin3,
                          "max |<psi| b b X b^+ |psi>| (X = K_psi, so(3,2))": amp_state}
check(f"[E] (a2) Paritätsauswahlregel: P = (-1)^(N_a + N_b) = (-1)^D vertauscht mit allen nativen Erzeugern ({par_nat:.0e}) "
      f"und hält psi fest ({par_psi:.0e}), also auch K_psi (1.1). Darum ist jede 1->2-Amplitude (ungerade Zahl von "
      f"Quanten) null: Pi_2 [K_psi, b^+] Omega = {max(v2_L[x] for x in (1e-4, 2e-4, 0.01)):.0e} (lambda bis 0,01; "
      f"{v2_L[0.1]:.0e} bei 0,1, Rundung im nicht störungsartigen Bereich), <psi| b b X b^+ |psi> = {amp_state:.0e} für K_psi "
      f"und so(3,2). Erlaubt ist der quartische Vertex 1->3: Pi_3 [K_psi, b^+] Omega = {v3_L[1e-4]:.3e} bei lambda = 1e-4, "
      f"Verhältnis zu 2e-4 {lin3:.4f} (linear), {v3_L[0.1]:.3f} bei lambda = 0,1: festgelegt allein durch die freie Amplitude",
      par_nat < 1e-12 and par_psi < 1e-12 and max(v2_L[x] for x in (1e-4, 2e-4, 0.01)) < 1e-9 and v2_L[0.1] < 1e-4
      and amp_state < 1e-9 and v3_L[1e-4] > 1e-6 and abs(lin3 - 2) < 1e-3 and v3_L[0.1] > 1e-2)
laufzeit()

# --- 2.3 Kompositzustand (b): N hartkernige Nähte, Original und Partner je ein Qubit, KMS-Purifikation
say("  2.3 Hartkerniger Kompositzustand (b)")
sp_ = np.array([[0, 0], [1, 0.]])                                    # sigma^+ = |1><0|
sx_ = np.array([[0, 1], [1, 0.]])
PAULI = [np.eye(2), sx_, np.array([[0, -1j], [1j, 0]]), np.diag([1., -1.])]


def naht_op(Oo, Op, i, N):
    """Operator Oo (x) Op auf Naht i (Original, Partner) im Produkt von N Nähten (Reihenfolge o_1 p_1 o_2 p_2 ...)."""
    mats = [sps.identity(4, format="csr")]*N
    mats[i] = sps.csr_matrix(np.kron(Oo, Op))
    out = mats[0]
    for Mx in mats[1:]:
        out = sps.kron(out, Mx, format="csr")
    return out.tocsr()


komp_tab, kum_tab = {}, {}
for N in (2, 3, 4, 5):
    tq = np.array([1, 0, 0, q])/math.sqrt(1 + q*q)                   # |00> + q |11>
    Psi = tq
    for _ in range(N - 1):
        Psi = np.kron(Psi, tq)
    # Anregungen relativ zu Psi: Einzelnahtbasis aus t und drei orthogonalen Anregungen; P2 = genau zwei Nähte angeregt
    Bn = np.linalg.qr(np.column_stack([tq, np.eye(4)[:, 1], np.eye(4)[:, 2], np.eye(4)[:, 3]]))[0]
    Bn[:, 0] *= np.sign(Bn[0, 0])
    BN = Bn
    for _ in range(N - 1):
        BN = np.kron(BN, Bn)
    zahl = (np.array(list(itertools.product(range(4), repeat=N))) > 0).sum(1)
    P2 = BN[:, zahl == 2]
    p2 = lambda v, P2=P2: float(np.linalg.norm(P2.conj().T @ v))
    # Einteilchenerzeuger auf Naht 0: |10><t| (Original angeregt)
    x1 = np.kron([0, 1.], [1, 0.])
    mats = [sps.identity(4, format="csr")]*N
    mats[0] = sps.csr_matrix(np.outer(x1, tq))
    b0 = mats[0]
    for Mx in mats[1:]:
        b0 = sps.kron(b0, Mx, format="csr")
    b0 = b0.tocsr()
    # einkörperige Operatoren: alle 16 Paulis je Naht, einzeln und kollektiv
    eink = []
    for i in range(N):
        for P_o in PAULI:
            for P_p in PAULI:
                eink.append(naht_op(P_o, P_p, i, N))
    koll = [sum((naht_op(P_o, P_p, i, N) for i in range(N)), sps.csr_matrix((4**N, 4**N), dtype=complex))
            for P_o in PAULI for P_p in PAULI]
    # Modulhamiltonoperator aus dem Zustand: rho_O (Original) über Partialspur
    Cn = Psi.reshape([2]*(2*N))
    perm = list(range(0, 2*N, 2)) + list(range(1, 2*N, 2))
    Cmat = np.transpose(Cn, perm).reshape(2**N, 2**N)
    U_, svn, Wh = np.linalg.svd(Cmat)
    lg = -2*np.log(svn)
    Ka_n = (U_*lg) @ U_.conj().T
    Kb_n = (Wh.T*lg) @ Wh.conj()
    k1 = np.diag(-np.log(np.array([1, q*q])/(1 + q*q)))
    soll = sum((np.kron(np.kron(np.eye(2**i), k1), np.eye(2**(N - i - 1))) for i in range(N)), np.zeros((2**N, 2**N)))
    einkoerper_K = float(np.abs(Ka_n - soll).max())
    # K auf den vollen Raum (Reihenfolge o_1 p_1 ...): über Umordnung
    def K_voll(v, Ka_n=Ka_n, Kb_n=Kb_n, N=N, perm=perm):
        T_ = np.transpose(v.reshape([2]*(2*N)), perm).reshape(2**N, 2**N)
        T2 = Ka_n @ T_ - T_ @ Kb_n.T
        inv = np.argsort(perm)
        return np.transpose(T2.reshape([2]*(2*N)), inv).reshape(-1)
    b0Psi = b0 @ Psi
    m_e = max(p2(X @ b0Psi - b0 @ (X @ Psi)) for X in eink + koll)
    m_K = p2(K_voll(b0Psi) - b0 @ K_voll(Psi))
    zwei_k = (naht_op(sx_, np.eye(2), 0, N) @ naht_op(sp_, np.eye(2), 1, N)).tocsr()
    m_2 = p2(zwei_k @ b0Psi - b0 @ (zwei_k @ Psi))
    Gp = sum((naht_op(sx_, sx_, i, N) for i in range(N)), sps.csr_matrix((4**N, 4**N))).tocsr()
    fN = 4*float(np.vdot(Psi, Gp @ (Gp @ Psi)).real - np.vdot(Psi, Gp @ Psi).real**2)
    komp_tab[N] = {"einkörperig max": m_e, "K_komp": m_K, "Zweikörper-Kontrolle": m_2, "K - Summe k_i": einkoerper_K,
                   "Fisher F_N der kollektiven Paarerzeugung Summe sigma^x_o sigma^x_p": fN}
    # Kumulanten des Kollektivfelds x = sum sigma^x_o / sqrt(2N) in Vakuum und KMS-Zustand
    for name_, rho1 in (("Vakuum", np.diag([1., 0.])), ("KMS", np.diag([1, q*q])/(1 + q*q))):
        xo = sum((np.kron(np.kron(np.eye(2**i), sx_), np.eye(2**(N - i - 1))) for i in range(N)),
                 np.zeros((2**N, 2**N)))/math.sqrt(2*N)
        rho = rho1
        for _ in range(N - 1):
            rho = np.kron(rho, rho1)
        mu = [np.trace(rho @ np.linalg.matrix_power(xo, k_)).real for k_ in range(5)]
        kum_tab[(N, name_)] = [r(mu[2]), r(mu[3]), r(mu[4] - 3*mu[2]**2)]
ZAHLEN["komposit"] = {str(k): v for k, v in komp_tab.items()}
ZAHLEN["komposit_kumulanten"] = {f"N = {k[0]}, {k[1]}": v for k, v in kum_tab.items()}
f_pro_naht = [v["Fisher F_N der kollektiven Paarerzeugung Summe sigma^x_o sigma^x_p"]/N for N, v in komp_tab.items()]
ok_komp = all(v["einkörperig max"] < 1e-12 and v["K_komp"] < 1e-10 and v["Zweikörper-Kontrolle"] > 1e-2
              and v["K - Summe k_i"] < 1e-10 for v in komp_tab.values()) and max(f_pro_naht) - min(f_pro_naht) < 1e-12
ok_kum = all(abs(v[0] - 0.5) < 1e-12 and abs(v[1]) < 1e-12 and abs(v[2] + 1/(2*k[0])) < 1e-12 for k, v in kum_tab.items())
# hartkernig auf Oszillatorebene: Leck der nativen so(3,2) aus dem Raum mit Besetzung <= 1 je Mode (d, e-Basis)
nv3 = 3
dd1, dd2, ee1, ee2 = moden(nv3, 4)
ZT3 = [dd1, dd2, dag(ee1), dag(ee2)]
ZT3d = [dag(x) for x in ZT3]
occ3 = np.array(list(itertools.product(range(nv3), repeat=4)))
hc = np.where(occ3.max(1) <= 1)[0]
aus = np.where(occ3.max(1) > 1)[0]
leck = []
for M in SO32:
    Am = 1j*G4 @ M
    H = sum((Am[al, be]*(ZT3d[al] @ ZT3[be]) for al in range(4) for be in range(4) if abs(Am[al, be]) > 1e-14),
            sps.csr_matrix((nv3**4, nv3**4), dtype=complex)).tocsr()
    leck.append(float(abs(H[aus][:, hc]).max()))
ZAHLEN["hartkern_leck"] = dict(zip(NAMEN, leck))
check(f"[E] (b) Kompositzustand aus N = 2..5 hartkernigen Nähten (Original und Partner je ein Qubit, KMS-Purifikation "
      f"|00> + q|11> je Naht): Der aus dem Zustand berechnete K ist exakt einkörperig (|K - Summe k_i| = "
      f"{max(v['K - Summe k_i'] for v in komp_tab.values()):.0e}); der Zweinahtanteil von [X, beta^+] Psi verschwindet für "
      f"alle einkörperigen X (16 je Naht, einzeln und kollektiv: max {max(v['einkörperig max'] for v in komp_tab.values()):.0e}) "
      f"und für K ({max(v['K_komp'] for v in komp_tab.values()):.0e}), die Zweikörper-Kontrolle sigma^x_1 sigma^+_2 gibt "
      f"{komp_tab[3]['Zweikörper-Kontrolle']:.6f}. Die quadratische Antwort ist extensiv (Fisher der kollektiven "
      f"Paarerzeugung F_N/N = {f_pro_naht[0]:.6f} für N = 2..5). Das Kollektivfeld x = Summe sigma^x/sqrt(2N) hat im "
      f"Vakuum und im "
      f"KMS-Zustand kappa_2 = 1/2, kappa_3 = 0, kappa_4 = -1/(2N) (z91; temperaturunabhängig, die Kurtosis von N "
      f"unabhängigen +-1-Variablen). Auf Oszillatorebene ist der hartkernige Raum (Besetzung <= 1) nicht invariant unter "
      f"der nativen so(3,2): Q_1..Q_3 lecken ({min(leck[7:]):.3f} bis {max(leck[7:]):.3f}, Paarerzeugung d^+ e^+), "
      f"Lorentz bis {max(leck[:6]):.3f}, nur E = Q_0 nicht ({leck[6]:.0e}); in der z91-Lesart (Naht = ein Leiterteilchen) "
      f"wirkt die Gruppe als 0 + Singleton, und die Komposition ist der freie Raum (z91: Sym^N(1 + V) = Summe_(k <= N) "
      f"Sym^k V, gezählt bis E = 60)",
      ok_komp and ok_kum and min(leck[7:]) > 0.1 and max(leck[:6]) > 0.1 and leck[6] < 1e-12)
befund("E (negativ für H1, Wechselwirkungsseite): Kein Kandidat gibt eine verbundene 1->2-Amplitude, weder mit den nativen "
       "Erzeugern noch mit dem eigenen Modulhamiltonoperator. (a1) ist quasifrei (W^+ X W = X + linear + const); (a2) ist "
       "nicht gaußsch, aber paritätsgerade (Modulladungen 0, -4), und P = (-1)^D verbietet jede ungerade Amplitude: es gibt "
       "nur den quartischen Vertex 1->3, linear in der freien Amplitude lambda; (b) ist ein Produkt, K und alle Erzeuger "
       "sind einkörperig, die Nichtgaußizität ist statistisch (kappa_4 = -1/(2N), kappa_3 = 0). Eine kubische Kopplung "
       "verlangt zugleich Nichtgaußizität und ungerade Modulladung; (a1) hat nur die Ladung, (a2) nur die "
       "Nichtgaußizität", True)
laufzeit()

# =================================================================== Teil 3
say("\nTeil 3: Beziehung zwischen b (bzw. G) und g, Bilanz und Urteil")
# (a1): Spur der Fisher-Überschüsse hängt nur von |beta| ab (SU(2) transitiv auf Spinorrichtungen), det(beta1, beta2) nicht


def fisher_ueberschuss_twistor(beta):
    """4 ||l_X Omega||^2 über die Basis, l_X Omega = (conj(r_1) d_1^+ + conj(r_2) d_2^+ - r_3 e_1^+ - r_4 e_2^+) Omega."""
    R_ = [r_vek(M, beta) for M in SO32]
    vec = [np.array([np.conj(x[0]), np.conj(x[1]), -x[2], -x[3]]) for x in R_]
    return 4*np.array([[np.vdot(vec[i], vec[j]).real for j in range(10)] for i in range(10)])


Gfrob = np.array([[np.trace(A.conj().T @ B).real for B in SO32] for A in SO32])
tr_inv = lambda F: float(np.trace(np.linalg.solve(Gfrob, F)))
ok_fiss = np.abs(fisher_ueberschuss_twistor(BETA) - ex1).max() < 1e-6
spuren = []
for _ in range(6):
    z = rng.normal(size=2) + 1j*rng.normal(size=2)
    z *= 0.03/np.linalg.norm(z)
    spuren.append(tr_inv(fisher_ueberschuss_twistor(z)))
norm_b = 0.03
b1v = np.array([norm_b, 0], complex)
dets = {}
for th_ in (0.0, 0.4, 1.0, np.pi/2):
    b2v = norm_b*np.array([math.cos(th_), math.sin(th_)], complex)
    dets[round(th_, 4)] = (abs(b1v[0]*b2v[1] - b1v[1]*b2v[0]),
                          tr_inv(fisher_ueberschuss_twistor(b1v)) + tr_inv(fisher_ueberschuss_twistor(b2v)))
# (a2): Verhältnis quartischer Vertex / Fisher-Überschuss hängt von der gewählten Gravitonkomponente ab
ratio_a2 = {}
for name_, h_ in (("|0,0;2,2>", H_GRAV), ("|0,0;4,0>", H_GRAV40)):
    vals = []
    for lam in (1e-4,):
        C_ = kond4(lam, h_).reshape(DA, DA)
        Ka_, Kb_ = modular_K(C_)
        Kv2 = lambda v: K_auf(Ka_, Kb_, v)
        KO = Kv2(OMEGA)
        g4 = max(pi3(Kv2(b @ OMEGA) - b @ KO) for b in ERZ)
        fex = float(np.abs(fisher(kond4(lam, h_)) - F_vak).max())
        vals = (g4, fex, g4**2/fex)
    ratio_a2[name_] = [r(x, 9) for x in vals]
# (b): Fisher additiv in N (aus dem Kompositmodell), kappa_4 = -1/(2N): N-unabhängiges Produkt
prod_b = {N: r(v["Fisher F_N der kollektiven Paarerzeugung Summe sigma^x_o sigma^x_p"]*kum_tab[(N, "KMS")][2], 12)
          for N, v in komp_tab.items()}
ZAHLEN["beziehung"] = {"(a1) Spur G^-1 F_Überschuss für zufällige Richtungen bei |beta| = 0,03": [r(x, 12) for x in spuren],
                       "(a1) |det(beta1, beta2)| und Spur bei festen Normen": {str(k): [r(a, 12), r(b, 12)]
                                                                               for k, (a, b) in dets.items()},
                       "(a2) (g4, Fisher-Überschuss, g4^2/Überschuss) je Komponente": ratio_a2,
                       "(b) F_N kappa_4(N) im Kompositmodell": prod_b}
dvals = [v[0] for v in dets.values()]
svals = [v[1] for v in dets.values()]
check(f"[E] Beziehung b <-> g: (a1) Die quadratischen Gravitationsdaten sind Summen über die Nähte und hängen nur von "
      f"|beta_s| ab (Spur G^-1 F_Überschuss für sechs zufällige Richtungen bei |beta| = 0,03: Streuung "
      f"{max(spuren) - min(spuren):.0e}, Twistorformel = Fock {ok_fiss}); die Klammerkopplung |det(beta1, beta2)| = "
      f"|beta1||beta2| sin theta durchläuft bei festen Normen {np.round(dvals, 6).tolist()} (Spur konstant "
      f"{svals[0]:.6f}): keine Funktion g(b). (a2) g4^2/Fisher-Überschuss = {ratio_a2['|0,0;2,2>'][2]:.6f} für |0,0;2,2> "
      f"gegen {ratio_a2['|0,0;4,0>'][2]:.6f} für |0,0;4,0>: lambda kürzt sich, die gewählte Gravitonkomponente nicht. (b) "
      f"Fisher extensiv (F_N = N F_1), kappa_4 = -1/(2N): F_N kappa_4 = {prod_b[2]} für N = 2..5 (die Groß-N-Regel 1/G ~ N, "
      f"g_stat ~ 1/N), aber ohne Krümmung und ohne kubischen Vertex",
      max(spuren) - min(spuren) < 1e-9 and ok_fiss and max(svals) - min(svals) < 1e-9 and max(dvals) - min(dvals) > 1e-4
      and abs(ratio_a2["|0,0;2,2>"][2] - ratio_a2["|0,0;4,0>"][2]) > 1e-3*abs(ratio_a2["|0,0;2,2>"][2])
      and max(prod_b.values()) - min(prod_b.values()) < 1e-10)
befund("L/C: Die einzige feste Beziehung ist die Groß-N-Skalierung des Kompositzustands (quadratische Antwort ~ N, "
       "statistische Kopplung kappa_4 = -1/(2N)); sie verbindet aber keine Krümmung mit keinem kubischen Vertex. Beim "
       "Kondensat sind Krümmung (zentral) und Klammerkopplung det(beta1, beta2) unabhängig; beim D = -4-Kondensat kürzt "
       "sich die Amplitude, die Richtung des Kondensats bleibt eine Wahl", True)

STATUS = {
    "Tomita-Familie und F": ("E, negativ", "jede Bahn J_g = U(g) J_psi U(g)^-1 hat die Maurer-Cartan-Form als "
                                           "Cartan-Verbindung, F_MM = 0 für jeden Zustand"),
    "Nullmoden": ("E", "Vakuum und Komposit Lorentz (AdS4); Kondensat (a1) Nullrotationen von beta, nicht reduktiv; zwei "
                       "Nähte mit det != 0 keine; D = -4-Kondensat reduktiv mit Torsion"),
    "modulare Berry-Krümmung": ("E", "Vakuum null; Kondensat zentral U(1), 2 Im(r_X^T G4 conj r_Y), quadratisch in beta, "
                                     "null auf Q x Q und Lorentz x Lorentz"),
    "quadratische Antwort und b": ("E, negativ", "Fisher fest nur im Vakuum (F = 0); beta^2 bzw. lambda-abhängig, Komposit "
                                                 "additiv; b multipliziert F ^ F = 0"),
    "1->2 (a1)": ("E, negativ", "quasifrei: W^+ X W = X + linear + const, alle Quadratiken und K_psi geben 0"),
    "1->2 (a2)": ("E, negativ", "Parität (-1)^D: ungerade Amplituden exakt 0; nur 1->3, linear in lambda"),
    "1->2 (b)": ("E, negativ", "K und Erzeuger einkörperig; kappa_3 = 0, kappa_4 = -1/(2N) statistisch"),
    "Beziehung b <-> g": ("L", "nur Groß-N-Skalierung des Komposits; Kondensat: unabhängig; D = -4: richtungsabhängig"),
    "H1": ("E negativ für (a1), (a2), (b)", "kein Kandidat gibt F != 0 oder eine kubische Kopplung, schon gar nicht "
                                           "beides mit festen Koeffizienten"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
urteil = ("H1 ist für die drei Kandidaten verneint, als gescoptes Negativ. (1) Gravitationsseite: Die Tomita-Familie "
          "J_g = U(g) J_psi U(g)^-1 eines beliebigen Zustands unter der nativen so(3,2) ist eine Bahn; ihre "
          "Cartan-Verbindung mit Rahmen in der nativen Gruppe ist die Maurer-Cartan-Form, also F_MM = 0, gaußsch oder "
          "nicht. Die Empfehlung aus z98 (deformierte Tomita-Familie) kann daher in dieser Form keine Krümmung geben. Der "
          "Zustand wählt nur die Nullmoden: Vakuum und Kompositzustand Lorentz (AdS4, torsionsfrei); das "
          "Impulsphasen-Kondensat aus z89 bricht Lorentz auf die zwei Nullrotationen, die beta in der c-Ebene festhalten "
          "(nilpotent, kein invariantes Komplement, also keine Zerlegung in Verbindung und Vierbein), und zwei Nähte mit "
          "det(beta1, beta2) != 0, also genau das Kondensat, das I_pi wählt, löschen alle Nullmoden; das wörtliche "
          "D = -4-Kondensat gibt eine reduktive, nicht symmetrische Bahn mit Torsion. Neu ist eine zentrale "
          "U(1)-Berry-Krümmung des Kondensats bezüglich der Nullmoden von K_psi, F(X, Y) = 2i Im(r_X^T G4 conj(r_Y)), "
          "quadratisch in der freien Amplitude beta; sie ist keine Lorentzkrümmung, verschwindet auf den AdS-Translationen "
          "Q x Q und auf Lorentz x Lorentz und koppelt nur Lorentz an Q. Die quadratische Antwort "
          "(Fisher-Information) ist nur im Vakuum fest (und verschwindet dort in der AdS-Zeitrichtung E), beim Kondensat "
          "beta- bzw. lambda-abhängig, beim Komposit additiv in N; b multipliziert F ^ F = 0 und bleibt unbestimmt. "
          "(2) Wechselwirkungsseite: Das Kondensat ist quasifrei (W^+ X W = X + linear + const), alle 64 Quadratiken, die "
          "nativen Erzeuger und das aus dem Zustand berechnete K_psi geben null verbundene 1->2-Anteile. Das D = -4-Kondensat "
          "ist nicht gaußsch, aber paritätsgerade; P = (-1)^D vertauscht mit allen nativen Erzeugern und mit K_psi, also "
          "verschwindet jede ungerade Amplitude exakt, und es bleibt ein quartischer 1->3-Vertex linear in lambda. Der "
          "Kompositzustand ist ein Produkt: K ist exakt einkörperig, die Nichtgaußizität ist die statistische "
          "kappa_4 = -1/(2N) bei kappa_3 = 0. Eine kubische Kopplung braucht zugleich ungerade Modulladung und "
          "Nichtgaußizität; das Kondensat hat nur die Ladung, das D = -4-Kondensat nur die Nichtgaußizität. (3) Beziehung: "
          "Beim Kondensat sind die quadratischen Gravitationsdaten Summen über die Nähte, abhängig nur von |beta_s|, während "
          "die Klammerkopplung det(beta1, beta2) bei festen Normen frei variiert; beim D = -4-Kondensat kürzt sich lambda im "
          "Verhältnis g4^2/Fisher, aber die gewählte Gravitonkomponente nicht; beim Komposit gilt nur die Groß-N-Regel "
          "(Antwort ~ N, kappa_4 ~ 1/N) als Lesart, ohne Krümmung und ohne kubischen Vertex. Folgen: K13 bleibt offen, und der "
          "Weg 'F aus einer deformierten Tomita-Bahn' ist geschlossen, für jeden Zustand; Krümmung verlangt eine Familie, "
          "die keine Bahn ist: entweder J(M_g, psi) mit festem Zustand und bewegter Algebra (nativ definiert, hier nicht "
          "gerechnet, der nächste Schritt) oder ein Zustandsprofil über AdS4 (eine Wahl). K15: Die kubische Kopplung ist auf 'ungerade Ladung und nicht gaußsch' eingegrenzt; z99s Zweikörperterm "
          "bleibt nötig. K20: Kondensat und Komposit sind keine gemeinsame Quelle; der Engpass zerfällt in zwei verschiedene "
          "Forderungen, eine nicht-Bahn-Familie (Gravitation) und einen paritätsungeraden Zweikörperterm (Wechselwirkung). "
          "Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
say(f"\n   (Laufzeit {time.time() - T0:.1f} s)")
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z100_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
