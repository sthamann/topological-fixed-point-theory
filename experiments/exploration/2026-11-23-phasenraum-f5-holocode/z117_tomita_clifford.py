"""Z117: Folgeauftrag U1 -- erzeugt der H6-Tomita-Operator die gemeinsame reelle Clifford-Struktur, und welche Chiralität
wählt die native Verklebung?

Leitlinie 'Rekonstruktion statt Auswahl'. H6 bleibt wie hergeleitet; gelöste Blöcke (E8, 16 = Lambda^gerade C^5, drei
Familien, Flavour, alpha) werden benutzt, nicht neu geöffnet. DIMENSION.UPLIFT.FIREWALL.01: Die zusätzlichen Clifford-Richtungen
sind algebraisch (Bivektoren des Twistorraums und die fünf Majoranas der 16), keine neuen physikalischen Dimensionen. Bestand:
z85/z89 (Quarteroszillatoren d, e; Z = (d_1, d_2, e_1^+, e_2^+), CCR-Form der Signatur (2,2), Tomita J = Tausch d <-> e,
tau = (Ad J)^+, pi- und omega-Ebene, U_theta = e^(i theta D)), z97 (Z_SU = i^(2h-2), Z4-Regel R = Z_SU^-1 unter P_A3),
z98 (AdS-Büschel I_pi - sigma I_omega, sigma = +1 nativ), z102 (Kleberroute zur Chiralität geschlossen), z103 (P_vol, Phase
psi der Twistor-Volumenform), z116 (Cl(3,11), Cl(4,12), J_F, R = i Gamma_11, Spur mit Uhr = +-64 eps, Verklebung I_T = I_F).
Statuskürzel: E exakt, C bedingt, L Lesart, O offen, N numerisch.

Teil 1  Natives J_H6 auf W (x) F: W = Einteilchen-Phasenraum der Quarteroszillatoren (T + T*), aus Fockoperatoren exakt
        gelesen: Tomita J_W, Adjungierte Dag, tau_W, Modulladung. Konforme Cliffordwirkung Gamma(Q) aus den Bivektoren des
        Twistorraums (Gamma^2 = Pf, kovariant unter der nativen su(2,2) der Fock-Bilinearen), Realformen V^psi (psi = Phase der
        Volumenform). Vorzeichen von J_W und Dag auf allen Gammas; J_H6 = J_W (x) J_F K: Quadrat, Gesamtchiralität,
        Verhältnis zur Majorana-Weyl-Struktur (Dag (x) J_F K, z116-C16 über einen Verflechter), Kommutanten, Typ der
        irreduziblen Kopie (reell/quaternionisch); Durchmusterung der gemeinsamen Cliffordwirkungen (Zwirn, Phase, psi).
Teil 2  Auswahl: Verklebung I_T = I_F <=> Gamma_ges = +1 <=> Z4-Regel R = Z_SU^-1 auf dem Spinormodul; Kandidaten (Naht/Hardy,
        KMS-Richtung, Uhr R gegen R^-1, Toeplitz-Index bzw. lichtartiger Schnitt, Modulfluss) unter den beiden Spiegeln, die
        S^+ und S^- tauschen; gemeinsame Moduldrehung; gemischte Antworten im KMS-Zustand.
Teil 3  Spur mit Uhr auf dem nativen Raum (Vorzeichentafel), P_vol und die Korrelation psi <-> Phase der Verklebung, Urteil.
"""
import itertools
import json
from pathlib import Path

import numpy as np
import scipy.sparse as sps

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
rng = np.random.default_rng(117)


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
    """C mit C g* C^-1 = s g für alle g (kleiner Faktor, SVD); None, wenn es keine gibt (wie z116)."""
    d = G[0].shape[0]
    rows = [np.kron(g.conj().T, np.eye(d)) - s*np.kron(np.eye(d), g) for g in G]
    _, sv, vh = np.linalg.svd(np.vstack(rows))
    if sv[-1] > 1e-9:
        return None
    C = vh[-1].conj().reshape(d, d, order="F")
    return C/abs(np.linalg.det(C))**(1/d)


def nullraum(M, tol=1e-9):
    _, s, vh = np.linalg.svd(np.atleast_2d(M))
    return vh[int(np.sum(s > tol*max(1.0, s[0]))):]


def clifford_ok(G):
    d = G[0].shape[0]
    eta = [np.round((g @ g)[0, 0], 9) for g in G]
    ok = all(np.allclose(G[m] @ G[n] + G[n] @ G[m], 2*(m == n)*eta[m]*np.eye(d)) for m in range(len(G)) for n in range(len(G)))
    ok &= all(abs(e.imag) < 1e-9 and abs(abs(e) - 1) < 1e-9 for e in eta)
    return ok, int(sum(e.real > 0 for e in eta)), int(sum(e.real < 0 for e in eta))


def zeichen(Mc, Gk):
    """antilinear C(w) = Mc conj(w): C Gk C^-1 = s Gk -> s (oder None)."""
    X = Mc @ Gk.conj() @ np.linalg.inv(Mc)
    for s in (1, -1, 1j, -1j):
        if np.allclose(X, s*Gk, atol=1e-10):
            return s
    return None


# =================================================================== Teil 1
say("Teil 1: Natives J_H6 -- Phasenraum der Quarteroszillatoren, konforme Cliffordwirkung, gemeinsame Realstruktur")


# --- 1.1 Einteilchen-Phasenraum W = T + T* aus den Fockoperatoren (z85: d_1, d_2, e_1, e_2; Fock 5^4)
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
                                     np.arange(DIMF))), shape=(DIMF, DIMF))
DV = OCC[:, 0] + OCC[:, 1] - OCC[:, 2] - OCC[:, 3]
DOP = sps.diags(DV.astype(complex)).tocsr()
LOW = np.where(OCC.sum(1) <= NV - 3)[0]
BAS = [d1, d2, dag(e1), dag(e2), dag(d1), dag(d2), (-e1).tocsr(), (-e2).tocsr()]       # (Z^alpha, Z-quer_alpha)
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


MJ, r1 = matrix_von(lambda X: (SW @ X.conj() @ SW).tocsr())      # J (w.B) J = sum conj(w_k) J B_k J
MD, r2 = matrix_von(dag)                                          # (w.B)^+ = sum conj(w_k) B_k^+
AD, r3 = matrix_von(lambda X: (DOP @ X - X @ DOP).tocsr())        # [D, .], D = N_d - N_e (Modulladung)
OMG = np.array([[complex(np.round((BAS[k] @ BAS[l] - BAS[l] @ BAS[k]).tocsc()[:, LOW][LOW, :].diagonal().mean(), 12))
                 for l in range(8)] for k in range(8)])
S4 = np.block([[np.zeros((2, 2)), np.eye(2)], [np.eye(2), np.zeros((2, 2))]]).astype(complex)
G4 = np.diag([1., 1., -1., -1.]).astype(complex)
TAUW = MD @ MJ.conj()
Z44 = np.zeros((4, 4), complex)
GAM7_W = np.diag([1.]*4 + [-1.]*4).astype(complex)                # +1 auf T, -1 auf T*
ok11 = (max(r1, r2, r3) < 1e-12 and np.allclose(MJ, MJ.real.round()) and np.allclose(MD, MD.real.round())
        and np.allclose(MJ @ MJ.conj(), np.eye(8)) and np.allclose(MD @ MD.conj(), np.eye(8))
        and np.allclose(TAUW, np.block([[S4, Z44], [Z44, -S4]])) and np.allclose(AD, np.diag([-1.]*4 + [1.]*4))
        and np.allclose(OMG, np.block([[Z44, np.eye(4)], [-np.eye(4), Z44]]))
        and np.allclose(MJ[:4, :4], 0) and np.allclose(MJ[4:, 4:], 0) and np.allclose(SW @ OMEGA_F, OMEGA_F))
ZAHLEN["phasenraum_W"] = {"M_J (Tomita, antilinear)": MJ.real.astype(int).tolist(), "M_Dag": MD.real.astype(int).tolist(),
                          "tau_W = M_Dag conj(M_J)": TAUW.real.astype(int).tolist()}
check("Phasenraum W = span{d_i, e_i^+} + span{d_i^+, e_i} = T + T* der Quarteroszillatoren, Basis (Z^alpha, Z-quer_alpha), aus "
      "den Fockoperatoren exakt gelesen (Fock 5^4, Residuum < 1e-12): Tomita J_W = Ad J (antilinear) und die Adjungierte Dag "
      "sind ganzzahlig, J_W^2 = 1 und Dag^2 = 1 exakt, beide bilden T auf T* ab; tau_W = Dag J_W = diag(S4, -S4) ist linear "
      "(auf T die z89-Matrix S4); die Modulladung [D, .] ist -1 auf T und +1 auf T*; die CCR-Form ist [[0, 1], [-1, 0]]; "
      "J Omega = Omega", ok11)

# --- 1.2 Konforme Cliffordwirkung aus den Bivektoren des Twistorraums
LC = np.zeros((4,)*4)
for _p in itertools.permutations(range(4)):
    LC[_p] = round(np.linalg.det(np.eye(4)[list(_p)]))


def stern(Q):
    return 0.5*np.einsum("abcd,cd->ab", LC, Q)


def pf(Q):
    return Q[0, 1]*Q[2, 3] - Q[0, 2]*Q[1, 3] + Q[0, 3]*Q[1, 2]


def Gam(Q):
    """Gamma(Q) auf W = T + T*: T -> T* über Q (obere Indizes), T* -> T über -*Q; Gamma(Q)^2 = Pf(Q)."""
    return np.block([[Z44, -stern(Q)], [Q, Z44]])


def AM(M):
    return np.block([[-M.T, Z44], [Z44, M]])


ZT, ZB = BAS[:4], BAS[4:]
okA = True
for _a, _b in itertools.product(range(4), repeat=2):
    _M = np.zeros((4, 4), complex)
    _M[_a, _b] = 1
    _QM = (ZB[_a] @ ZT[_b]).tocsr()
    _A, _r = matrix_von(lambda X: (_QM @ X - X @ _QM).tocsr())
    okA &= np.allclose(_A, AM(_M)) and _r < 1e-10
Qr = rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4))
Qr = Qr - Qr.T
Mr = rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4))
Mr -= np.trace(Mr)/4*np.eye(4)
ok_quad = np.allclose(Gam(Qr) @ Gam(Qr), pf(Qr)*np.eye(8))
ok_kov = np.allclose(AM(Mr) @ Gam(Qr) - Gam(Qr) @ AM(Mr), Gam(Mr @ Qr + Qr @ Mr.T))
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


def orthonorm(V):
    Gm = np.array([[(Gam(X) @ Gam(Y) + Gam(Y) @ Gam(X))[0, 0].real/2 for Y in V] for X in V])
    w, U = np.linalg.eigh(Gm)
    return [(int(np.sign(w[k])), sum(U[j, k]*V[j] for j in range(len(V)))/np.sqrt(abs(w[k]))) for k in range(len(w))]


u_om = [np.array([1, 0, -1, 0.]), np.array([0, 1, 0, -1.])]
v_pi = [np.array([1, 0, 1, 0.]), np.array([0, 1, 0, 1.])]
I_PI = (np.outer(u_om[0], u_om[1]) - np.outer(u_om[1], u_om[0])).astype(complex)       # z89/z103
I_OM = (np.outer(v_pi[0], v_pi[1]) - np.outer(v_pi[1], v_pi[0])).astype(complex)
RF = {}
for psi in (0.0, np.pi):
    V = realform(psi)
    R11 = orthonorm(tau_teil(V, +1))
    R31 = orthonorm(tau_teil(V, -1))
    R11.sort(key=lambda t: -t[0])                       # e_+ (Gamma^2 = +1) zuerst
    R31.sort(key=lambda t: (-t[0] if psi == 0 else t[0]))  # psi = 0: Zeit (+1) zuerst; psi = pi: Zeit (-1) zuerst
    RF[psi] = {"V": V, "R11": R11, "R31": R31}
ok_dim = all(len(RF[p]["V"]) == 6 and len(RF[p]["R11"]) == 2 and len(RF[p]["R31"]) == 4 for p in RF)
sig = {p: (sorted([s for s, _ in RF[p]["R31"]]), sorted([s for s, _ in RF[p]["R11"]])) for p in RF}
pi_reell = (np.allclose(stern(G4 @ I_PI.conj() @ G4), I_PI) and np.allclose(stern(G4 @ I_OM.conj() @ G4), I_OM))
spin = []
for p in RF:
    gs = [Gam(X) for _, X in RF[p]["R31"] + RF[p]["R11"]]
    gens = [(gs[a] @ gs[b] - gs[b] @ gs[a])/4 for a, b in itertools.combinations(range(6), 2)]
    # su(2,2) auf W: A_M mit M^+ G + G M = 0 und Spur 0 (hier: Bild der Spinerzeuger liegt in {A_M}, Rang 15)
    in_AM = all(np.allclose(g[:4, 4:], 0) and np.allclose(g[4:, :4], 0) and np.allclose(-g[:4, :4].T, g[4:, 4:]) and
                abs(np.trace(g[4:, 4:])) < 1e-12 and np.allclose(g[4:, 4:].conj().T @ G4 + G4 @ g[4:, 4:], 0) for g in gens)
    rang = np.linalg.matrix_rank(np.array([np.concatenate([g.real.ravel(), g.imag.ravel()]) for g in gens]))
    spin.append(in_AM and rang == 15)
e_p, e_m = RF[0.0]["R11"][0][1], RF[0.0]["R11"][1][1]
tau_dil = [c for c in (1, -1) if np.allclose(TAUW, c*Gam(e_p) @ Gam(e_m))]
ZAHLEN["konforme_clifford"] = {"Signatur R^(3,1), R^(1,1) bei psi = 0 (Gamma^2-Zählung)": [sig[0.0][0], sig[0.0][1]],
                               "dto. psi = pi": [sig[np.pi][0], sig[np.pi][1]],
                               "Pf(I_pi - I_om) (AdS, z98)": pf(I_PI - I_OM).real, "Pf(I_pi + I_om) (dS)": pf(I_PI + I_OM).real}
check(f"Cliffordwirkung aus Bivektoren: Gamma(Q) = [[0, -*Q], [Q, 0]] erfüllt Gamma^2 = Pf(Q) ({ok_quad}); die native u(2,2) "
      f"der Fock-Bilinearen Z-quer_a Z^b wirkt auf W exakt als diag(-M^T, M) ({okA}, alle 16 Erzeuger) und Gamma ist unter "
      f"ihrem spurfreien Teil kovariant ({ok_kov}); die Realformen V^psi = {{G Q-quer G = e^(i psi) *Q}} sind 6-dimensional "
      f"({ok_dim}); I_pi und I_omega sind reell bei psi = 0 ({pi_reell}); R^(1,1) = tau-fest = span{{I_pi, I_omega}}, R^(3,1) = "
      f"tau-ungerade; Gamma^2-Signaturen psi = 0: R^(3,1) {sig[0.0][0]}, R^(1,1) {sig[0.0][1]} (meist minus), psi = pi: "
      f"{sig[np.pi][0]}, {sig[np.pi][1]}; die Spinerzeuger [Gamma, Gamma]/4 liegen in der nativen su(2,2) und spannen sie "
      f"(Rang 15) für beide psi ({spin}); tau_W = {tau_dil} Gamma(e_+) Gamma(e_-), also das R^(1,1)-Volumenelement (Dilatation), "
      f"nicht die 4D-Chiralität", okA and ok_quad and ok_kov and ok_dim and pi_reell and all(spin) and len(tau_dil) == 1
      and sig[0.0] == ([-1, -1, -1, 1], [-1, 1]) and sig[np.pi] == ([-1, 1, 1, 1], [-1, 1]))

# --- 1.3 Vorzeichen der beiden antilinearen Involutionen auf den sechs Gammas
VZ = {}
for p in RF:
    VZ[p] = {"Dag R31": [zeichen(MD, Gam(X)) for _, X in RF[p]["R31"]], "Dag R11": [zeichen(MD, Gam(X)) for _, X in RF[p]["R11"]],
             "J_W R31": [zeichen(MJ, Gam(X)) for _, X in RF[p]["R31"]], "J_W R11": [zeichen(MJ, Gam(X)) for _, X in RF[p]["R11"]]}
ZAHLEN["vorzeichen"] = {("psi=0" if p == 0 else "psi=pi"): {k: [str(x) for x in v] for k, v in d.items()} for p, d in VZ.items()}
ok13 = (VZ[0.0]["Dag R31"] == [-1]*4 and VZ[0.0]["Dag R11"] == [-1]*2 and VZ[0.0]["J_W R31"] == [-1]*4
        and VZ[0.0]["J_W R11"] == [1]*2 and VZ[np.pi]["Dag R31"] == [1]*4 and VZ[np.pi]["J_W R11"] == [-1]*2)
check(f"Vorzeichen (C Gamma C^-1 = s Gamma): Dag ist auf allen sechs Gammas einheitlich (psi = 0: s = -1, psi = pi: s = +1), "
      f"also die Majorana-Struktur der konformen Cl(2,4) bzw. Cl(4,2); das Tomita-J_W hat auf R^(3,1) das Vorzeichen von Dag "
      f"(psi = 0: -1) und auf R^(1,1) das entgegengesetzte (+1): J_W = Dag tau_W ist Majorana für Lorentz und spiegelt die "
      f"beiden konformen Zusatzrichtungen ({ZAHLEN['vorzeichen']})", ok13)

# --- 1.4 Gemeinsamer Raum W (x) F (256): J_H6, Quadrat, Gesamtchiralität, Verhältnis zu J_MW
Gi = euklid(5)
G11 = chiral(Gi)
JF = np.eye(32, dtype=complex)
for _j in range(5):
    JF = JF @ Gi[2*_j]
R_F = 1j*G11
a5 = [(Gi[2*j] + 1j*Gi[2*j + 1])/2 for j in range(5)]
PAR = np.eye(32, dtype=complex)
for _j in range(5):
    PAR = PAR @ (np.eye(32) - 2*a5[_j].conj().T @ a5[_j])
G11 = G11*np.sign((G11 @ PAR)[0, 0].real)          # +1 auf Lambda^gerade (16, enthält das Fockvakuum)
R_F = 1j*G11
I32 = np.eye(32)
I8 = np.eye(8)
M_H6 = np.kron(MJ, JF)                               # J_H6 = J_W (x) J_F K
M_MW = np.kron(MD, JF)                               # J_MW = Dag (x) J_F K
GTOT = np.kron(GAM7_W, G11)
JH6_q = M_H6 @ M_H6.conj()
ok_q = np.allclose(JH6_q, np.eye(256))
ok_chir = np.allclose(M_H6 @ GTOT.conj(), GTOT @ M_H6)
U_rel = M_MW @ M_H6.conj()
ok_rel = np.allclose(U_rel, np.kron(TAUW, I32))
lor = {p: [Gam(RF[p]["R31"][a][1]) @ Gam(RF[p]["R31"][b][1]) for a, b in itertools.combinations(range(4), 2)] for p in RF}
mix = {p: [Gam(RF[p]["R31"][a][1]) @ Gam(RF[p]["R11"][b][1]) for a in range(4) for b in range(2)] for p in RF}
ok_lor = all(np.allclose(MJ @ L.conj(), L @ MJ) for p in RF for L in lor[p])
ok_su_nicht = all(np.allclose(MJ @ L.conj(), -L @ MJ) for p in RF for L in mix[p])
spin10 = [Gi[a] @ Gi[b] for a, b in itertools.combinations(range(10), 2)]
ok_int = all(np.allclose(JF @ X.conj(), X @ JF) for X in spin10) and np.allclose(JF @ R_F.conj(), R_F @ JF)
ok_mod = np.allclose(MJ @ (1j*AD).conj(), (1j*AD) @ MJ) and np.allclose(MD @ (1j*AD).conj(), (1j*AD) @ MD)
check(f"J_H6 = J_W (x) J_F K auf W (x) F (256): J_H6^2 = +1 exakt ({ok_q}); vertauscht mit der Gesamtchiralität Gamma_7 (x) "
      f"Gamma_11 ({ok_chir}), mit den sechs Lorentzerzeugern ({ok_lor}), mit Spin(10) und der Uhr R ({ok_int}) und mit dem "
      f"Modulfluss U_theta ({ok_mod}); die acht gemischten konformen Erzeuger Gamma_a Gamma_(+-) dreht es um ({ok_su_nicht}), also "
      f"nur der Kommutant von tau in su(2,2). Exakt gilt J_MW^-1 J_H6 = tau_W (x) 1 ({ok_rel}) mit J_MW = Dag (x) J_F K: das "
      f"native Tomita-J ist die Majorana-Weyl-Konjugation mal der Dilatationsspiegelung tau_W", ok_q and ok_chir and ok_rel
      and ok_lor and ok_su_nicht and ok_int and ok_mod)


# --- 1.5 Durchmusterung der gemeinsamen Cliffordwirkungen; Kommutant und Typ der irreduziblen Kopie
def kommutant(Xs, Ys):
    """Kommutant der von X (x) 1 und Y (x) gamma_i erzeugten Algebra in End(W) (x) span{1, Gamma_11}: Paare (A, B)."""
    def loese(kom_X, kom_Y):
        rows = []
        for X in Xs:
            rows.append(np.kron(np.eye(8), X) - np.kron(X.T, np.eye(8)) if kom_X else None)
        for Y in Ys:
            rows.append(np.kron(np.eye(8), Y) - np.kron(Y.T, np.eye(8)) if kom_Y else np.kron(np.eye(8), Y) + np.kron(Y.T, np.eye(8)))
        M = np.vstack([r for r in rows if r is not None])
        return [v.reshape(8, 8, order="F") for v in nullraum(M)]
    As = loese(True, True)
    Bs = loese(True, False)
    return [(A, np.zeros((8, 8), complex)) for A in As] + [(np.zeros((8, 8), complex), Bm) for Bm in Bs]


def kopientyp(Xs, Ys, MW_):
    """Reelle Form des Kommutanten unter c -> J c J^-1 (J = MW_ (x) J_F K); +: gespalten (J hält eine irreduzible Kopie, MW auf
    der Kopie), -: nicht gespalten (quaternionisch bzw. komplex, J vertauscht die Kopien)."""
    K = kommutant(Xs, Ys)
    sigma = lambda A, Bm: (MW_ @ A.conj() @ np.linalg.inv(MW_), -(MW_ @ Bm.conj() @ np.linalg.inv(MW_)))
    basis = []
    for A, Bm in K:
        basis += [(A, Bm), (1j*A, 1j*Bm)]
    rows = []
    for A, Bm in basis:
        sA, sB = sigma(A, Bm)
        rows.append(np.concatenate([(sA - A).real.ravel(), (sA - A).imag.ravel(), (sB - Bm).real.ravel(), (sB - Bm).imag.ravel()]))
    reell = []
    for v in nullraum(np.array(rows).T):
        reell.append((sum(c*b[0] for c, b in zip(v.real, basis)), sum(c*b[1] for c, b in zip(v.real, basis))))
    q = []
    for A, Bm in reell:
        A0 = A - np.trace(A)/8*np.eye(8)
        if np.abs(A0).max() + np.abs(Bm).max() < 1e-9:
            continue
        S_A, S_B = A0 @ A0 + Bm @ Bm, A0 @ Bm + Bm @ A0
        if np.allclose(S_B, 0, atol=1e-8) and np.allclose(S_A, S_A[0, 0]*np.eye(8), atol=1e-8):
            q.append(S_A[0, 0].real)
    gespalten = any(x > 1e-9 for x in q) or len(K) == 1          # irreduzibel: J wirkt auf dem Modul selbst
    return len(K), len(reell), gespalten


def gemeinsam(psi, zwirn, phi, konform):
    R31 = [Gam(X) for _, X in RF[psi]["R31"]]
    R11 = [Gam(X) for _, X in RF[psi]["R11"]]
    G7 = chiral(R31 + R11)
    G7 = G7*np.sign(G7[0, 0].real)
    C = {"Gamma_7": G7, "gamma_5 = i Gamma_0123": chiral(R31), "Gamma(e_+)": R11[0], "Gamma(e_-)": R11[1]}[zwirn]
    Xs = R31 + (R11 if konform else [])
    Ys = [phi*C]
    G = [np.kron(X, I32) for X in Xs] + [np.kron(phi*C, g) for g in Gi]
    ok, npos, nneg = clifford_ok(G)
    s_H6 = [zeichen(MJ, X) for X in Xs] + [zeichen(MJ, phi*C)]
    s_MW = [zeichen(MD, X) for X in Xs] + [zeichen(MD, phi*C)]
    einh = lambda s: None if (None in s or len(set(s)) > 1) else s[0]
    Gt = chiral(G)
    typ_H6 = kopientyp(Xs, Ys, MJ) if einh(s_H6) is not None else None
    typ_MW = kopientyp(Xs, Ys, MD) if einh(s_MW) is not None else None
    schnitt = all(np.allclose(TAUW @ X, X @ TAUW) for X in Xs + Ys)
    return {"Clifford": ok, "Signatur": (npos, nneg), "s(J_H6)": einh(s_H6), "s(J_MW)": einh(s_MW),
            "J_H6 vertauscht mit Chiralität": bool(np.allclose(M_H6 @ Gt.conj(), Gt @ M_H6)),
            "J_MW vertauscht mit Chiralität": bool(np.allclose(M_MW @ Gt.conj(), Gt @ M_MW)),
            "Typ J_H6 (dim Kommutant, dim reell, gespalten)": typ_H6, "Typ J_MW": typ_MW, "hält den Schnitt (tau_W)": schnitt}


SCAN = {}
for psi in (0.0, np.pi):
    for zw in ("Gamma_7", "gamma_5 = i Gamma_0123", "Gamma(e_+)", "Gamma(e_-)"):
        for phi in (1, 1j):
            for konform in ((True, False) if zw == "Gamma_7" else (False,)):
                SCAN[(("0" if psi == 0 else "pi"), zw, ("1" if phi == 1 else "i"), "konform" if konform else "Lorentz")] = \
                    gemeinsam(psi, zw, phi, konform)
mw_H6 = sorted([k for k, v in SCAN.items() if v["s(J_H6)"] is not None and v["Typ J_H6 (dim Kommutant, dim reell, gespalten)"]
                and v["Typ J_H6 (dim Kommutant, dim reell, gespalten)"][2] and v["J_H6 vertauscht mit Chiralität"]])
mw_MW = sorted([k for k, v in SCAN.items() if v["s(J_MW)"] is not None and v["Typ J_MW"] and v["Typ J_MW"][2]
                and v["J_MW vertauscht mit Chiralität"]])
ZAHLEN["scan_gemeinsame_clifford"] = {" | ".join(k): {kk: (str(vv) if not isinstance(vv, (bool, int, float, type(None))) else vv)
                                                     for kk, vv in v.items()} for k, v in SCAN.items()}
ZAHLEN["MW_auf_Kopie_J_H6"] = [" | ".join(k) for k in mw_H6]
ZAHLEN["MW_J_MW"] = [" | ".join(k) for k in mw_MW]
konf_H6 = [k for k in mw_H6 if k[3] == "konform"]
mw_H6_s = [k for k in mw_H6 if SCAN[k]["hält den Schnitt (tau_W)"]]
mw_H6_r11 = [k for k in mw_H6 if not SCAN[k]["hält den Schnitt (tau_W)"]]
sig_H6 = sorted(set(SCAN[k]["Signatur"] for k in mw_H6))
psiphi_H6 = sorted(set((k[0], k[2]) for k in mw_H6_s))
mw_MW_konf = [k for k in mw_MW if k[3] == "konform"]
psiphi_MW = sorted(set((k[0], k[2]) for k in mw_MW_konf))
ok15 = (all(v["Clifford"] for v in SCAN.values()) and not konf_H6 and len(mw_H6_s) == 4
        and set(k[1] for k in mw_H6_s) == {"Gamma_7", "gamma_5 = i Gamma_0123"}
        and all(s in ((11, 3), (3, 11)) for s in sig_H6) and psiphi_H6 == [("0", "1"), ("pi", "i")]
        and len(mw_MW_konf) == 2 and psiphi_MW == [("0", "1"), ("pi", "i")])
check(f"Durchmusterung ({len(SCAN)} gemeinsame Cliffordwirkungen: psi in {{0, pi}}, Zwirn in {{Gamma_7, gamma_5, Gamma(e_+), "
      f"Gamma(e_-)}}, Phase phi in {{1, i}}, konform oder nur Lorentz; alle Clifford): J_H6 ist Majorana-Weyl auf einer "
      f"irreduziblen Kopie (einheitliches Vorzeichen, vertauscht mit der Chiralität, reelle Form des Kommutanten gespalten) nur auf "
      f"Lorentz x innen, Signatur {sig_H6} = Cl(3,11) bis aufs Gesamtvorzeichen, in keinem konformen Fall ({konf_H6}). Mit "
      f"Erhalt des nativen Schnitts bleiben genau die Zwirne Gamma_7 und gamma_5 mit (psi, phi) in {psiphi_H6} ({len(mw_H6_s)} "
      f"Fälle); die R^(1,1)-Zwirne {mw_H6_r11} brechen den Schnitt. J_MW = Dag (x) J_F K ist MW auf der vollen konformen Algebra "
      f"(Cl(4,12) bzw. Cl(12,4)) genau für (psi, phi) in {psiphi_MW}", ok15)

# Verhältnis auf den Kopien: tau_W (x) 1 im Kommutanten von Cl(3,11), nicht von Cl(4,12)
R31_0 = [Gam(X) for _, X in RF[0.0]["R31"]]
R11_0 = [Gam(X) for _, X in RF[0.0]["R11"]]
G7_0 = chiral(R31_0 + R11_0)
G7_0 = G7_0*np.sign(G7_0[0, 0].real)
CL311 = [np.kron(X, I32) for X in R31_0] + [np.kron(G7_0, g) for g in Gi]
ok_kom311 = all(np.allclose(np.kron(TAUW, I32) @ X, X @ np.kron(TAUW, I32)) for X in CL311)
ok_anti = all(np.allclose(TAUW @ X, -X @ TAUW) for X in R11_0)
P_kop = {s_: np.kron((I8 + s_*TAUW)/2, I32) for s_ in (1, -1)}
gleich_auf_kopie = all(np.allclose(P_kop[s_] @ (M_H6 - s_*M_MW) @ P_kop[s_].conj(), 0) for s_ in (1, -1))
gp_, gm_ = Gam(I_PI)/2, Gam(I_OM)/2                  # lichtartige gamma_+ ~ I_pi (Unendlichkeitstwistor), gamma_- ~ I_omega
null_ok = np.allclose(gp_ @ gp_, 0) and np.allclose(gm_ @ gm_, 0) and np.allclose(gp_ @ gm_ + gm_ @ gp_, I8)
P_cut = gp_ @ gm_                                    # z116: P_- = gamma_+ gamma_-, Sektor von A^+A, trägt die Indexnullmode
pi_ebene = np.array([v_pi[0], v_pi[1]]).T
om_ebene = np.array([u_om[0], u_om[1]]).T
schnitt_pi = np.linalg.matrix_rank(P_cut[:4, :4]) == 2 and np.allclose(P_cut[:4, :4] @ om_ebene, 0) \
    and np.allclose(P_cut[:4, :4] @ pi_ebene, pi_ebene)
schnitt_tau = np.allclose(P_cut, (I8 + TAUW)/2)
check(f"Verhältnis: tau_W (x) 1 liegt im Kommutanten der Lorentz-innen-Algebra Cl(3,11) ({ok_kom311}) und antivertauscht mit "
      f"den beiden R^(1,1)-Gammas ({ok_anti}); auf den beiden tau-Kopien gilt J_H6 = +-J_MW ({gleich_auf_kopie}). Die lichtartigen "
      f"gamma_+ = Gamma(I_pi)/2, gamma_- = Gamma(I_omega)/2 ({null_ok}); mit D_perp = gamma_- A + gamma_+ A^+ (z116) liegt die "
      f"Nullmode des Index +1 im Sektor P_0 = gamma_+ gamma_- = Kern von gamma_+; P_0 ist auf T genau die pi-Ebene ({schnitt_pi}) "
      f"und auf ganz W gleich (1 + tau_W)/2 ({schnitt_tau}): Die tau-Kopie mit J_H6 = +J_MW ist der Schnitt des Nahtindex. Also: J_H6 ist auf dem Schnitt genau die MW-Struktur, auf dem Komplement ihr "
      f"Negatives; verschieden sind beide nur auf der konformen Cl(4,12)", ok_kom311 and ok_anti and gleich_auf_kopie
      and null_ok and schnitt_pi and schnitt_tau)

# --- 1.6 Vergleich mit z116: Cl(4,12)-Konstruktion aus euklidischen Faktoren, Verflechter V8 (psi = pi)
E6 = euklid(3)
hz = [E6[0], E6[1], E6[2], E6[3], 1j*E6[4], 1j*E6[5]]            # z116: drei Raum, e_4, Zeit, e_-1
R31_p = [Gam(X) for _, X in RF[np.pi]["R31"]]                    # Zeit (-1) zuerst, dann drei Raum (+1)
R11_p = [Gam(X) for _, X in RF[np.pi]["R11"]]                    # e_+ (+1), e_- (-1)
nat = [R31_p[1], R31_p[2], R31_p[3], R11_p[0], R31_p[0], R11_p[1]]
rows = [np.kron(np.eye(8), N) - np.kron(H.T, np.eye(8)) for H, N in zip(hz, nat)]     # N V - V H = 0
V8 = nullraum(np.vstack(rows))
V8 = V8[0].reshape(8, 8, order="F") if len(V8) == 1 else None
h7, n7 = chiral(hz), chiral(nat)
G16 = [kr(x, I32) for x in hz] + [kr(h7, gi)*(1j if np.allclose(h7 @ h7, np.eye(8)) else 1.0) for gi in Gi]
Gt16 = chiral(G16)
C16 = None
for s6 in (1, -1):
    for s10 in (1, -1):
        C8, C10 = konj(hz, s6), konj(Gi, s10)
        if C8 is None or C10 is None:
            continue
        C = np.kron(C8, C10)
        for s in (1, -1):
            if all(np.allclose(C @ g.conj() @ np.linalg.inv(C), s*g) for g in G16):
                if np.allclose(C @ C.conj(), np.eye(256)) and np.allclose(C @ Gt16.conj(), Gt16 @ C):
                    C16, s16 = C, s
ok16 = V8 is not None and C16 is not None
if ok16:
    V = np.kron(V8, I32)
    C16t = V @ C16 @ np.linalg.inv(V.conj())
    GT16 = V @ Gt16 @ np.linalg.inv(V)                                    # Gesamtchiralität transportiert
    ref = np.unravel_index(np.argmax(np.abs(M_MW)), M_MW.shape)
    kand_ = {"J_MW": M_MW, "J_MW Gamma_ges": M_MW @ GT16.conj()}
    rel16 = {}
    for nm_, Mk in kand_.items():
        ix = np.unravel_index(np.argmax(np.abs(Mk)), Mk.shape)
        lam_k = C16t[ix]/Mk[ix]
        rel16[nm_] = (bool(np.allclose(C16t, lam_k*Mk, atol=1e-8)), complex(np.round(lam_k, 10)))
    PS = (np.eye(256) + GT16)/2
    lamS = C16t[ref]/M_MW[ref] if abs(M_MW[ref]) > 0 else 0
    auf_Splus = np.allclose(PS @ C16t @ PS.conj(), rel16["J_MW"][1]*PS @ M_MW @ PS.conj(), atol=1e-8) if rel16["J_MW"][0] else \
        np.allclose(PS @ C16t @ PS.conj(), rel16["J_MW Gamma_ges"][1]*PS @ M_MW @ PS.conj(), atol=1e-8)
    prop_H6 = any(np.allclose(C16t, c_*M_H6, atol=1e-8) for c_ in (C16t[ix_]/M_H6[ix_] for ix_ in
                                                                     [np.unravel_index(np.argmax(np.abs(M_H6)), M_H6.shape)]))
    n7_ok = np.allclose(V8 @ h7 @ np.linalg.inv(V8), n7)
    gt_ok = np.allclose(GT16, GTOT) or np.allclose(GT16, -GTOT)
    ZAHLEN["z116_vergleich"] = {"s16": s16, "C16 transportiert": {k: [v[0], str(v[1])] for k, v in rel16.items()},
                                "auf S^+ gleich lambda J_MW": bool(auf_Splus)}
    prop_MW = any(v[0] and abs(abs(v[1]) - 1) < 1e-9 for v in rel16.values())
check(f"z116-Konstruktion (Cl(4,12) aus euklid(3) x Cl(10), C16 aus den Faktoren mit C C-quer = 1 und Chiralität, s16 = "
      f"{s16 if ok16 else None}): Der Verflechter V8 auf die nativen Gammas (psi = pi, Ordnung drei Raum, e_+, Zeit, e_-) existiert "
      f"und ist eindeutig ({V8 is not None}), bildet h7 auf Gamma_7 ab ({ok16 and n7_ok}) und die Gesamtchiralitäten aufeinander "
      f"({ok16 and gt_ok}); transportiert ist C16 = lambda J_MW bzw. lambda J_MW Gamma_ges mit |lambda| = 1 "
      f"({ZAHLEN.get('z116_vergleich', {}).get('C16 transportiert')}), auf S^+ also genau lambda J_MW ({ok16 and auf_Splus}); "
      f"proportional zu J_H6 ist es nicht ({ok16 and not prop_H6})", ok16 and prop_MW and auf_Splus and not prop_H6 and n7_ok
      and gt_ok)

# --- 1.7 Restfreiheit der Cliffordwirkung: Gamma_lambda = Konjugation mit diag(lambda, 1)
rest = {}
for lab, lam_ in (("1", 1.0), ("e^0,7i", np.exp(0.7j)), ("1,6", 1.6), ("1,6 e^0,3i", 1.6*np.exp(0.3j))):
    Dl = np.diag([lam_]*4 + [1]*4)
    gl = [Dl @ Gam(X) @ np.linalg.inv(Dl) for _, X in RF[0.0]["R31"] + RF[0.0]["R11"]]
    rest[lab] = ([zeichen(MD, g) for g in gl], [zeichen(MJ, g) for g in gl])
U_mod = np.diag(np.exp(1j*0.35*np.diag(AD)))        # U_theta auf W, theta = 0,35: diag(e^-i theta, e^i theta) ~ diag(e^0,7i, 1)
gl_mod = [U_mod @ Gam(X) @ np.linalg.inv(U_mod) for _, X in RF[0.0]["R31"] + RF[0.0]["R11"]]
Dl7 = np.diag([np.exp(-0.7j)]*4 + [1]*4)
ist_modul = all(np.allclose(a_, Dl7 @ Gam(X) @ np.linalg.inv(Dl7)) for a_, (_, X) in zip(gl_mod, RF[0.0]["R31"] + RF[0.0]["R11"]))
ok17 = (rest["e^0,7i"] == rest["1"] and all(None in rest[k][0] and None in rest[k][1] for k in ("1,6", "1,6 e^0,3i"))
        and ist_modul)
check(f"Restfreiheit: Die su(2,2)-kovarianten Cliffordwirkungen bilden die Familie Gamma_lambda (Konjugation mit diag(lambda, 1) auf "
      f"T + T*); Dag und J_W erzwingen |lambda| = 1 (für |lambda| = 1,6 kein einheitliches Vorzeichen: {rest['1,6'][0]}), die "
      f"Phase bleibt frei und ist genau der Modulfluss U_theta ({ist_modul}; z98: Phase von I), der mit allen nativen Daten "
      f"vertauscht", ok17)
befund("E: J_H6 erzeugt die gemeinsame reelle Struktur auf Lorentz x innen. Das native Tomita-J_H6 = J_W (x) J_F K hat Quadrat "
       "+1, vertauscht mit der Gesamtchiralität, Lorentz, Spin(10), Uhr und Modulfluss, und ist auf jeder tau-Kopie genau die "
       "Majorana-Weyl-Konjugation von Cl(3,11) (bis auf das Vorzeichen der Kopie; auf dem Schnitt mit +). Auf der konformen "
       "Cl(4,12) ist es nicht die "
       "MW-Konjugation: J_H6 = J_MW (tau_W (x) 1), tau_W = Gamma(e_+) Gamma(e_-) spiegelt die beiden Zusatzrichtungen (CPT-artige "
       "Modulkonjugation). Die MW-Struktur der Cl(4,12) ist Dag (x) J_F K; beide sind nativ und erzeugen einander über tau_W = "
       "Dag J_W", True)

# =================================================================== Teil 2
say("\nTeil 2: Auswahl -- Gesamtchiralität, Verklebung, Schnitt")
PP = {s_: (np.eye(256) + s_*GTOT)/2 for s_ in (1, -1)}
I_T = 1j*GAM7_W                                     # komplexe Struktur von T (i auf T, -i auf T*)
Z_SU_W = -1j*GAM7_W                                 # z103: Z_SU = i^(D-2) wirkt auf T als -i
ok_zsu = np.allclose(Z_SU_W, np.diag(np.exp(1j*np.pi/2*np.diag(AD))))   # = Moduldrehung U_(pi/2) auf W
glue_p = np.kron(I_T, I32) - np.kron(I8, R_F)
glue_m = np.kron(I_T, I32) + np.kron(I8, R_F)
dim_glue = (256 - np.linalg.matrix_rank(glue_p), 256 - np.linalg.matrix_rank(glue_m))
ok_glue = (np.allclose(glue_p @ PP[1], 0) and np.allclose(glue_m @ PP[-1], 0) and dim_glue == (128, 128)
           and np.allclose(np.kron(Z_SU_W, R_F), GTOT))
check(f"Verklebung = Gesamtchiralität = Z4-Regel: Mit I_T = i auf T (= Z_SU^-1; Z_SU wirkt auf W als die Moduldrehung "
      f"U_(pi/2) ({ok_zsu})) und I_F = R = i Gamma_11 ist der Unterraum I_T (x) 1 = 1 (x) I_F genau S^+ = (T x 16) + (T* x 16-quer) "
      f"(Dimension {dim_glue[0]}), I_T = -I_F genau S^-; zudem Z_SU (x) R = Gamma_ges exakt: Die holomorphe Verklebung ist die "
      f"z97-Regel R = Z_SU^-1, eingeschränkt auf den Spinormodul ({ok_glue})", ok_glue and ok_zsu)

# Spiegel, die S^+ und S^- tauschen: welche Daten hält welcher Spiegel fest?
PHI_F = np.kron(I8, Gi[0])                          # innere Spiegelung (Pin(10)), linear
T0g = Gam(RF[0.0]["R31"][0][1])                     # Zeitgamma (psi = 0, Gamma^2 = +1)
PHI_T = np.kron(T0g, I32)                           # Spiegelung im R^(3,1) (Pin(3,1)), linear, bildet T auf T* ab
M_PHI_K = np.kron(I8, JF)                           # K_W (x) J_F K: antilinear, Fockbasis reell
tauscht = {nm: bool(np.allclose((P_ @ PP[1].conj() if anti else P_ @ PP[1]) @ np.linalg.inv(P_), PP[-1]))
           for nm, P_, anti in (("1 (x) gamma_1", PHI_F, False), ("Gamma_t (x) 1", PHI_T, False), ("K_W (x) J_F K", M_PHI_K, True))}


def verhalten(X, P_, anti, gen=False):
    """+1: festgehalten, -1: umgekehrt (X -> -X bzw. bei Erzeugern Fluss umgekehrt), 0: sonst. gen: X ist ein Erzeuger i H."""
    Y = (P_ @ X.conj() @ np.linalg.inv(P_)) if anti else (P_ @ X @ np.linalg.inv(P_))
    if np.allclose(Y, X):
        return 1
    if np.allclose(Y, -X):
        return -1
    return 0


U_mod_W = lambda th: np.diag(np.exp(1j*th*np.diag(AD)))
DATEN = {
    "Naht/Hardy: Gamma_7 (T gegen T*)": (np.kron(GAM7_W, I32), False),
    "KMS-Richtung: Modulerzeuger i D": (np.kron(1j*AD, I32), True),
    "Schnitt des Index: P_0": (np.kron(P_cut, I32), False),
    "tau_W": (np.kron(TAUW, I32), False),
    "Uhr R = i Gamma_11": (np.kron(I8, R_F), False),
    "Tomita J_H6": (M_H6, None),
}
tab = {}
for nm, (X, gen) in DATEN.items():
    zeile = {}
    for snm, P_, anti in (("1 (x) gamma_1", PHI_F, False), ("Gamma_t (x) 1", PHI_T, False), ("K_W (x) J_F K", M_PHI_K, True)):
        if gen is None:                                # antilinearer Operator: P J P^-1 = +-J?
            Pb = P_.conj() if not anti else P_.conj()
            Y = (P_ @ X @ np.linalg.inv(Pb)) if not anti else (P_ @ X.conj() @ np.linalg.inv(P_.conj()))
            zeile[snm] = 1 if np.allclose(Y, X) else (-1 if np.allclose(Y, -X) else 0)
        else:
            zeile[snm] = verhalten(X, P_, anti)
    tab[nm] = zeile
ZAHLEN["spiegel_tafel"] = tab
# jedes Einzeldatum wird von einem S^+-S^--Tausch festgehalten, der J_H6 bis aufs Vorzeichen erhält
einzel_fest = {nm: [snm for snm, v in z.items() if v == 1] for nm, z in tab.items() if nm != "Tomita J_H6"}
jh6_ok = {snm: tab["Tomita J_H6"][snm] != 0 for snm in tauscht}
alle_einzeln = all(len(v) > 0 for v in einzel_fest.values())
# Paar (i D, R): kein Tausch hält beide; ihr Produkt ist Gamma_ges
paar = [snm for snm in tauscht if tab["KMS-Richtung: Modulerzeuger i D"][snm] == 1 and tab["Uhr R = i Gamma_11"][snm] == 1]
prod_ok = np.allclose(np.kron(U_mod_W(np.pi/2), R_F), GTOT) and np.allclose(np.kron(U_mod_W(np.pi/2), np.linalg.inv(R_F)), -GTOT)
check(f"Spiegel: 1 (x) gamma_1 (innen), Gamma_t (x) 1 (Raumzeit) und K_W (x) J_F K (antilinear) vertauschen S^+ und S^- "
      f"({tauscht}) und erhalten J_H6 bis aufs Vorzeichen ({jh6_ok}). Wirkung (+1 fest, -1 umgekehrt) auf die Kandidaten: {tab}. "
      f"Jedes Einzeldatum wird von mindestens einem Tausch festgehalten ({alle_einzeln}: {einzel_fest}), wählt also nicht. Kein "
      f"Tausch hält Modulrichtung i D und Uhr R zugleich ({paar}); ihr Produkt ist exakt die Gesamtchiralität: U_(pi/2) (x) R = "
      f"Gamma_ges, U_(pi/2) (x) R^-1 = -Gamma_ges ({prod_ok})", all(tauscht.values()) and all(jh6_ok.values()) and alle_einzeln
      and paar == [] and prod_ok)

# gemeinsamer Modulfluss: kein Element ist auf beiden Faktoren die Vierteldrehung
paritaet = all((1 + 4*k) != 4*m*(1 + 4*j) and (1 + 4*k) != 4*m*(3 + 4*j)
               for m in range(1, 9) for k in range(-50, 50) for j in range(-50, 50))
tt = np.linspace(0, 4, 40001)
treffer = {m: int(np.sum((np.abs(np.exp(4j*np.pi*m*tt) - 1j) < 1e-3) & ((np.abs(np.exp(-1j*np.pi*tt) - 1j) < 1e-3)
                                                                         | (np.abs(np.exp(-1j*np.pi*tt) + 1j) < 1e-3))))
           for m in (1, 2, 3)}
t_uhr = [np.exp(4j*np.pi*m*0.5) for m in (1, 2, 3)]
# Richtung: U_(pi/2) = sigma_t auf T bei t_T = -1/(8m), R = e^(iK/2) = sigma_t auf den Ladungen bei t_R = -1/2 (beide negativ)
gleichsinnig = all(np.isclose(np.exp(4j*np.pi*m*(-1/(8*m))), np.exp(-1j*np.pi/2)) for m in (1, 2, 3)) and \
    np.isclose(np.exp(-1j*np.pi*(-0.5)*1), 1j)
check(f"Gemeinsamer Modulfluss: Auf T wirkt sigma_t als e^(4 pi i m t) (z85: K = 4 pi m D), auf den Ladungen (K/pi = Q_a - Q_b "
      f"mod 4, z76) als e^(-i pi t (Q_a - Q_b)). Vierteldrehung auf T: t = (1 + 4k)/(8m); auf der Uhr: t = 1/2 + 2j bzw. 3/2 + 2j; "
      f"beides zugleich hieße 1 + 4k (ungerade) = 4m(1 + 4j) bzw. 4m(3 + 4j) (gerade), unmöglich ({paritaet}; Gitter t in [0, 4], "
      f"Treffer {treffer}); zur Uhrzeit t = 1/2 ist der Fluss auf T die Eins ({np.allclose(t_uhr, 1)}). U_(pi/2) (x) R ist also "
      f"kein Element des gemeinsamen Modulflusses. Beide Faktoren sind aber gleichsinnig modular: U_(pi/2) = sigma_(-1/(8m)) auf T, "
      f"R = sigma_(-1/2) auf den Ladungen ({gleichsinnig})", paritaet and all(v == 0 for v in treffer.values())
      and np.allclose(t_uhr, 1) and gleichsinnig)

# gemischte Antworten im KMS-Zustand (Produktzustand Omega (x) rho_F)
NVk = 4
k4 = fock(NVk, 4)
ZTk = [k4[0], k4[1], dag(k4[2]), dag(k4[3]), dag(k4[0]), dag(k4[1]), (-k4[2]).tocsr(), (-k4[3]).tocsr()]
om = np.zeros(NVk**4)
om[0] = 1.0
GT = np.array([[np.vdot(X @ om, Y @ om) for Y in ZTk] for X in ZTk])        # <B_k^+ B_l>_Omega
erw_B = max(abs(np.vdot(om, X @ om)) for X in ZTk)
NF = [a5[j].conj().T @ a5[j] for j in range(5)]
NTOT = np.real(np.diag(sum(NF)))
vak = np.zeros((32, 32))
vak[np.argmin(NTOT), np.argmin(NTOT)] = 1.0
vak_in_16 = abs(G11[np.argmin(NTOT), np.argmin(NTOT)] - 1) < 1e-12
HF = sum(k*NF[j] for j, k in enumerate((0, 1, 2, 3, 0)))                 # Trägercharaktere i^k (z82)
rhoG = np.diag(np.exp(-np.real(np.diag(HF))))
zustaende = {"Fockvakuum (in der 16)": vak, "Gibbs, H = sum k n_k": rhoG/np.trace(rhoG), "Spurzustand": np.eye(32)/32}
gem = {}
for nm, rho in zustaende.items():
    Cm = np.kron(GT, rho)
    ev_p = np.sort(np.linalg.eigvalsh(PP[1] @ Cm @ PP[1]))
    ev_m = np.sort(np.linalg.eigvalsh(PP[-1] @ Cm @ PP[-1]))
    verb = max(abs(np.vdot(om, X @ om)*np.trace(rho @ Gi[i]) - np.vdot(om, X @ om)*np.trace(rho @ Gi[i]))
               + abs(np.vdot(om, X @ om)) for X in ZTk for i in range(10))
    gem[nm] = {"isospektral S^+ / S^-": bool(np.allclose(ev_p, ev_m)), "max |<B (x) gamma>|": float(verb),
               "Spur S^+": float(np.trace(PP[1] @ Cm).real), "Spur S^-": float(np.trace(PP[-1] @ Cm).real)}
ZAHLEN["gemischte_antworten"] = gem
diag_T = np.round(np.diag(GT).real, 9).tolist()
check(f"Gemischte Antworten im KMS-Zustand (Twistorfaktor im Vakuum Omega von d, e, der Purifikation von a, b; innen Fockvakuum "
      f"(in der 16: {vak_in_16}), Gibbs mit den Trägercharakteren, Spurzustand): <B> = 0 ({erw_B:.0e}), also verschwinden alle "
      f"Oszillator-Majorana-Zweipunktfunktionen, und die Antwortmatrix <Psi^+ Psi'> ist auf S^+ und S^- isospektral ({gem}). "
      f"<B^+ B> auf (Z, Z-quer) = {diag_T}: In S^+ antwortet die 16 über die Partner-Erzeuger e^+ in T, in S^- über d^+ in T*; "
      f"J_T (x) K_F bildet das eine auf das andere ab", vak_in_16 and erw_B < 1e-14
      and all(v["isospektral S^+ / S^-"] and v["max |<B (x) gamma>|"] < 1e-12 for v in gem.values()))

# Schnitt x Gesamtchiralität -> Helizität der 16 (R = +i); 4D-Chiralität nativ orientiert: gamma_5 = Gamma_7 tau_W
g5 = GAM7_W @ TAUW
g5_ist = [c for c in (1, -1) if np.allclose(g5, c*chiral(R31_0))]
heli = {}
for sG in (1, -1):
    for cut_nm, Pc in (("Schnitt pi (Index +1)", P_cut), ("Schnitt omega (Index -1)", gm_ @ gp_)):
        Pphys = PP[sG] @ np.kron(Pc, I32)
        X = Pphys @ np.kron(I8, (np.eye(32) + G11)/2)
        hel = np.trace(X @ np.kron(g5, I32)).real/max(np.trace(X).real, 1e-12)
        heli[(("S+" if sG == 1 else "S-"), cut_nm)] = {"dim (komplex)": int(round(np.trace(Pphys).real)),
                                                        "gamma_5 der R = +i-Komponente": float(np.round(hel, 9)),
                                                        "dim R = +i": int(round(np.trace(X).real))}
ZAHLEN["helizitaet"] = {" | ".join(k): v for k, v in heli.items()}
hv = {k: v["gamma_5 der R = +i-Komponente"] for k, v in heli.items()}
ok_hel = (all(abs(abs(x) - 1) < 1e-9 for x in hv.values())
          and hv[("S+", "Schnitt pi (Index +1)")] == -hv[("S-", "Schnitt pi (Index +1)")]
          and hv[("S+", "Schnitt pi (Index +1)")] == -hv[("S+", "Schnitt omega (Index -1)")]
          and all(v["dim (komplex)"] == 64 and v["dim R = +i"] == 32 for v in heli.values())
          and np.allclose(M_H6 @ np.kron(P_cut, I32).conj(), np.kron(P_cut, I32) @ M_H6)
          and np.allclose(np.kron(P_cut, I32) @ GTOT, GTOT @ np.kron(P_cut, I32)) and len(g5_ist) == 1)
check(f"Schnitt und Helizität: gamma_5 := Gamma_7 tau_W ist die 4D-Chiralität ({g5_ist} x chiral(R^(3,1))); der Schnitt P_0 ist "
      f"J_H6-verträglich und vertauscht mit Gamma_ges. S^+- geschnitten hat 64 komplexe Dimensionen, mit J_H6 also 64 reelle = "
      f"eine Weyl-16 (32 komplexe Komponenten mit R = +i, eine Helizität): {ZAHLEN['helizitaet']}. Die Helizität der 16 ist das "
      f"Produkt aus Gesamtchiralität und Schnitt, zwei getrennte Bits", ok_hel)
befund("E/C: Keine native Größe wählt die Gesamtchiralität. Verklebung I_T = I_F, Gesamtchiralität Gamma_ges = +1 und die auf den "
       "Spinormodul eingeschränkte z97-Regel R = Z_SU^-1 sind exakt dieselbe Bedingung. Jedes Einzeldatum (Naht/Hardy, "
       "KMS-Richtung, Schnitt, tau, Uhr) wird von einem S^+-S^--Tausch festgehalten, der J_H6 erhält; nur das Paar (Modulrichtung "
       "auf T, Uhrrichtung innen) bestimmt sie, und zwar als Produkt U_(pi/2) (x) R = Gamma_ges. Dieses Produkt ist kein Element des "
       "gemeinsamen Modulflusses (Paritätsargument), und die gemischten Antworten des Produktzustands sind isospektral. Minimale "
       "Prämisse: ein Bit, P_glue = 'der Spinormodul ist neutral unter dem Produkt der gleichsinnigen modularen Vierteldrehungen "
       "U_(pi/2) (x) R' (gruppentheoretisch: die Symmetrie ist (SU(2,2) x Spin(10))/Z4_diag; gleichwertig: die 16 mit R = +i ist an "
       "T = span{a_i, b_i^+} verklebt). Das Vorzeichen hängt dann nur noch am Hebungszweig der Uhr auf der 16 (z82: R = z_X^(+-1)). "
       "Der Index legt unabhängig davon den Schnitt fest; die Helizität der 16 ist das Produkt beider Bits", True)

# =================================================================== Teil 3
say("\nTeil 3: Spur mit Uhr, P_vol und Urteil")
eps = np.zeros((4, 4, 4, 4))
for p in itertools.permutations(range(4)):
    eps[p] = round(np.linalg.det(np.eye(4)[list(p)]))
paare = [(a, b) for a in range(4) for b in range(a + 1, 4)]


def spurtafel(psi, sG, kopie, uhr):
    R31g = [Gam(X) for _, X in RF[psi]["R31"]]
    # native Orientierung: Gamma_0123 so, dass i Gamma_0123 = gamma_5 = Gamma_7 tau_W
    if not np.allclose(1j*R31g[0] @ R31g[1] @ R31g[2] @ R31g[3], g5):
        R31g[3] = -R31g[3]
    Gs = [np.kron(g, I32) for g in R31g]
    P_ = PP[sG] @ P_kop[kopie]
    Ru = np.kron(I8, R_F if uhr == 1 else np.linalg.inv(R_F))
    werte, metr = {}, {}
    for (a, b) in paare:
        for (c, d) in paare:
            X = Gs[a] @ Gs[b] @ Gs[c] @ Gs[d]
            werte[(a, b, c, d)] = np.trace(P_ @ Ru @ X)
            metr[(a, b, c, d)] = np.trace(P_ @ X)
    k_ = {k: v/eps[k] for k, v in werte.items() if eps[k] != 0}
    konst = set(np.round(list(k_.values()), 8))
    null_sonst = all(abs(v) < 1e-9 for k, v in werte.items() if eps[k] == 0)
    metr_nur = all(abs(v) < 1e-9 for k, v in metr.items() if {k[0], k[1]} != {k[2], k[3]})
    return (list(konst)[0] if len(konst) == 1 and null_sonst else None), metr_nur


SPUR = {}
for psi in (0.0, np.pi):
    for sG in (1, -1):
        for kop in (1, -1):
            for uhr in (1, -1):
                c_, m_ = spurtafel(psi, sG, kop, uhr)
                SPUR[(("0" if psi == 0 else "pi"), ("S+" if sG == 1 else "S-"), ("pi-Kopie" if kop == 1 else "omega-Kopie"),
                      ("R" if uhr == 1 else "R^-1"))] = (c_, m_)
ZAHLEN["spur_mit_uhr"] = {" | ".join(k): str(v[0]) for k, v in SPUR.items()}
c0 = SPUR[("0", "S+", "pi-Kopie", "R")][0]
regel = all(v[0] is not None and v[1] and abs(v[0] - c0*(1 if k[1] == "S+" else -1)*(1 if k[2] == "pi-Kopie" else -1)
                                               *(1 if k[3] == "R" else -1)) < 1e-8 for k, v in SPUR.items())
psi_blind = all(abs(SPUR[("0",) + k[1:]][0] - SPUR[("pi",) + k[1:]][0]) < 1e-8 for k in SPUR if k[0] == "0")
check(f"Spur mit Uhr auf dem nativen Raum (R^(3,1) nativ orientiert durch i Gamma_0123 = Gamma_7 tau_W): Auf S^+-+ und der pi- bzw. "
      f"omega-Kopie gilt Tr(R Gamma_a Gamma_b Gamma_c Gamma_d) = c eps_abcd für alle 36 Paare, ohne Uhr nur die Metrik; c = {c0} "
      f"für (S^+, pi, R), und c = 64 x (Gesamtchiralität) x (Schnitt) x (Uhrrichtung) in allen 16 Fällen ({regel}); unabhängig "
      f"von psi ({psi_blind}). Die Uhrspur misst die Helizität der R = +i-Materie (c/64 = gamma_5 der 16)", regel and psi_blind
      and abs(c0 - 64) < 1e-8)

# P_vol: U_(pi/4) ist eine Symmetrie der ganzen gemeinsamen Struktur und kippt psi
U4 = np.diag(np.exp(1j*np.pi/4*np.diag(AD)))


def Gam_plus(Q):
    return np.block([[Z44, stern(Q)], [Q, Z44]])        # andere Blocknormierung (kappa = +1): Gamma~^2 = -Pf


basis0 = [X for _, X in RF[0.0]["R31"] + RF[0.0]["R11"]]
V0_bild = [U4 @ Gam(X) @ np.linalg.inv(U4) for X in basis0]
kippt = all(np.allclose(Y, Gam_plus(1j*X)) for Y, X in zip(V0_bild, basis0)) and \
    all(np.allclose(stern(G4 @ (1j*X).conj() @ G4), -(1j*X)) for X in basis0)                 # i X liegt in V^pi
sig_erh = all(np.allclose(Y @ Y, Gam(X) @ Gam(X)) for Y, X in zip(V0_bild, basis0))
U4g = np.kron(U4, I32)
sym = (np.allclose(U4g @ M_H6 @ np.linalg.inv(U4g.conj()), M_H6) and np.allclose(U4g @ M_MW @ np.linalg.inv(U4g.conj()), M_MW)
       and np.allclose(U4 @ GAM7_W, GAM7_W @ U4) and np.allclose(U4 @ TAUW, TAUW @ U4) and np.allclose(U4 @ P_cut, P_cut @ U4))
# die transportierte Lorentz-innen-Algebra (phi = 1) bleibt MW für J_H6: Vorzeichen und Typ
X_bild = V0_bild[:4]
s_bild = [zeichen(MJ, X) for X in X_bild] + [zeichen(MJ, GAM7_W)]
typ_bild = kopientyp(X_bild, [GAM7_W], MJ)
mw_bild = len(set(s_bild)) == 1 and None not in s_bild and typ_bild[2]
spur_U4 = np.trace(PP[1] @ P_kop[1] @ np.kron(I8, R_F) @ np.kron(V0_bild[0] @ V0_bild[1] @ V0_bild[2] @ V0_bild[3], I32))
spur_0 = np.trace(PP[1] @ P_kop[1] @ np.kron(I8, R_F) @ np.kron(Gam(basis0[0]) @ Gam(basis0[1]) @ Gam(basis0[2]) @ Gam(basis0[3]),
                                                                 I32))
zeit_q = {("0" if p == 0 else "pi"): (Gam(RF[p]["R31"][0][1]) @ Gam(RF[p]["R31"][0][1]))[0, 0].real for p in RF}
inv_zeit = all(abs((1 if k[2] == "1" else -1) - zeit_q[k[0]]) < 1e-12 for k in mw_H6_s) and len(mw_H6_s) == 4
ok_pvol = kippt and sig_erh and sym and mw_bild and abs(spur_U4 - spur_0) < 1e-8 and inv_zeit
check(f"P_vol (z103): Die Moduldrehung U_(pi/4) bildet die Cliffordwirkung auf V^0 auf die auf V^pi ab (Gamma(X) -> Gamma~(i X) "
      f"mit i X in V^pi, andere Blocknormierung; {kippt}), erhält die Signatur ({sig_erh}) und vertauscht mit J_H6, J_MW, Gamma_7, "
      f"tau_W und dem Schnitt ({sym}); die transportierte Lorentz-innen-Algebra mit phi = 1 bleibt MW für J_H6 (Vorzeichen "
      f"{s_bild}, Typ {typ_bild}: {mw_bild}); die Uhrspur bleibt ({np.round(spur_0, 6)} -> {np.round(spur_U4, 6)}). Die Kopplung "
      f"(psi, phi) der Durchmusterung gilt bei fester Blocknormierung; invariant ist nur: Die inneren Gammas haben das "
      f"Clifford-Vorzeichen der Zeit ({inv_zeit}). P_vol ist eine Modul-Eichrichtung der gemeinsamen Struktur", ok_pvol)
befund("C: Dieselbe Wahl fixiert nicht P_vol, aber sie ersetzt es. Die Uhrspur Tr(R Gamma Gamma Gamma Gamma) = 64 x "
       "(Chiralität) x (Schnitt) x (Uhr) eps ist U_(pi/4)-invariant und damit blind für P_vol; mit Uhr ergibt sie die "
       "eps-Kontraktion (Gravitationspaar), ohne Uhr nur die Metrik (topologisches Paar), und ihr Vorzeichen ist die Helizität "
       "der R = +i-Materie. In der Form 'Wirkung = Uhrspur über den gemeinsamen MW-Spinor' wählt also die Uhr das Paar und "
       "P_glue (mit Schnitt und Uhrrichtung) den Einzelsektor; P_vol selbst ist eine Modul-Eichrichtung (U_(pi/4) ist Symmetrie "
       "von J_H6, Uhr, Schnitt und MW-Struktur). Prämisse dafür: die Uhrspur ist die Wirkung (z116 offen, hier Lesart)", True)

urteil = ("U1 beantwortet. (1) Erzeugung [E]: Das native Tomita-J_H6 = J_W (x) J_F K auf W (x) F (W = Einteilchen-Phasenraum der "
          "Quarteroszillatoren, T + T*, aus den Fockoperatoren exakt gelesen) hat J_H6^2 = +1 exakt und vertauscht mit der "
          "Gesamtchiralität, Lorentz, Spin(10), der Uhr und dem Modulfluss. Die konforme Cliffordwirkung ist nativ (Bivektoren des "
          "Twistorraums, Gamma^2 = Pf, Spinerzeuger = die su(2,2) der Fock-Bilinearen) bis auf die Modulphase, eine Symmetrie. Auf "
          "der Lorentz-innen-Algebra ist J_H6 Majorana-Weyl auf jeder irreduziblen Kopie, mit erzwungener Signatur Cl(3,11) (die "
          "inneren Gammas haben das Clifford-Vorzeichen der Zeit) und erzwungenem Zwirn (Gamma_7 oder gamma_5, schnitterhaltend). "
          "Auf der konformen Cl(4,12) ist es nicht die MW-Konjugation: J_H6 = J_MW (tau_W (x) 1), J_MW = Dag (x) J_F K ist die "
          "MW-Struktur der Cl(4,12) (z116s C16 bis auf Phase bzw. Gamma_ges, auf S^+ gleich), tau_W = Dilatationsspiegelung der "
          "beiden Zusatzrichtungen. Der Schnitt des Nahtindex ist genau die tau-Kopie (1 + tau_W)/2, und dort gilt J_H6 = J_MW: "
          "Auf dem physikalischen Schnitt erzeugt das Tomita-J die gemeinsame MW-Struktur. (2) Auswahl [E negativ]: Verklebung "
          "I_T = I_F, Gesamtchiralität und die auf den Spinormodul eingeschränkte z97-Regel R = Z_SU^-1 sind dieselbe Bedingung; "
          "jedes native Einzeldatum wird von einem J_H6-verträglichen S^+-S^--Tausch festgehalten; nur das Paar (Modulrichtung auf T, "
          "Uhr innen) bestimmt sie, als Produkt U_(pi/2) (x) R = Gamma_ges, und dieses Produkt ist kein Element des gemeinsamen "
          "Modulflusses; die gemischten Antworten sind isospektral. Der Index wählt unabhängig den Schnitt (pi-Ebene); die "
          "Helizität der 16 ist das Produkt beider Bits. (3) Urteil: Chiralität über die gemeinsame reelle Clifford-Struktur ist "
          "bedingt, mit der minimalen Prämisse P_glue (ein Bit: der Spinormodul ist neutral unter dem Produkt der gleichsinnigen "
          "modularen Vierteldrehungen U_(pi/2) (x) R, gruppentheoretisch (SU(2,2) x Spin(10))/Z4_diag, d. h. die 16 ist an T = "
          "span{a_i, b_i^+} verklebt); das Vorzeichen hängt dann nur am Hebungszweig der Uhr auf der 16 (z82). Mit P_glue und "
          "dem nativen Schnitt bleibt "
          "genau eine Weyl-16 einer Helizität; ohne P_glue ist der Inhalt 16 + 16-quer, vektorartig. Die Uhrspur ist 64 x "
          "Chiralität x Schnitt x Uhr mal eps, U_(pi/4)-invariant und blind für P_vol; als Wirkung gelesen wählt die Uhr das "
          "Gravitationspaar (eps statt Metrik) und P_glue das Vorzeichen, P_vol wird zur Modul-Eichrichtung. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z117_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
