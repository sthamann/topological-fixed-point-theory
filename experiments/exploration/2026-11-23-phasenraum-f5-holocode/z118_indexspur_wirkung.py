"""Z118: Folgeauftrag U2 -- Liefert die native Quelle den Elternoperator D_perp und die Spur mit Uhr als
Gravitationswirkung, und welches Gewicht normiert G?

Leitlinie 'Rekonstruktion statt Auswahl'. Firewalls: GRAV.NONCIRCULAR.01 (die Uhr R tritt nur als Gewicht bzw. Symmetrie
auf, nie als Zeit), GRAV.PARENT.GAP.NO_GO.01 (kein Elternoperator wird vorausgesetzt: D_perp wird aus nativen Stücken
rekonstruiert oder als nicht nativ ausgewiesen), DIMENSION.UPLIFT.FIREWALL.01 (alles ist Twistor-, Clifford- und
Toeplitzkinematik einer Naht plus Formenalgebra; keine wechselwirkende 3+1-Schließung). Gelöste Blöcke (Flavour, Massen,
alpha, 16, E8) bleiben unberührt; der Ledgerwert G = v_geo^-2 aus GRAV.NONLINEAR.01 (v359) wird nur gelesen; nichts wird
befördert. Bestand: z116 (Algebra der siebzehnten Analyse), z98 (MM-Algebra nativ, sigma = +1, l = R_Vak/2), z103 (P_vol,
Sektorgruppe Z_4, Z_SU, P_A3), z104 (Kappe, Wick-Abbildung, positive Zeit => Lambda <= 0), z100 (Maurer-Cartan), z47/z48
(Index +1, Polgewicht 1 - q^2, Randkopplung), z83 (L^nat, E = L_0, Uhr auf T trivial), z57 (Flächenregel).
Statuskürzel: E exakt, C bedingt, L Lesart, O offen.

Teil 1  D_perp nativ? Klein-Clifford-Algebra der Twistorbivektoren auf T + T* (Volumenform eps): die drei nativen
        Nullpaare der radialen sl(2,R) (Dilatation: Ebenen von I_pi und I_omega; Kappe Q^-: Spitzen; E: Original und
        Partner, komplex), lichtartiger Schnitt P_+- auf T, Chiralitätsidentität Gamma_7 = c Gamma_2 gamma_5, die Uhr auf
        S^+ (256) = Z_SU^-1 (P_A3); Uhrladungen und Indizes der nativen Schritte (Viertelschritt, A_q, Transfer
        J_-1 = S^4, Möbiussenker L_1) und die Uhrkovarianz von D_perp.
Teil 2  Wirkung: zulässige Gewichte der relativen Spur (Kommutant und Endlichkeit), die native Uhr auf beiden Faktoren,
        MacDowell-Mansouri mit nativen Daten (Euler, Einstein, Lambda, Holst), Vorzeichenbits.
Teil 3  Normierung von G: globale gegen lokale Kopplung (Index gegen Randgewicht 1 - q^2, Störungen), Ganzzahligkeit der
        Spur, J-Summanden (Dirac- gegen Majorana-Weyl-Zählung), Vergleich mit GRAV.NONLINEAR.01.
Teil 4  Bilanz und Urteil (K13, Gravitationsachse).
"""
import itertools
import json
import time
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.linalg as sla
import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
rng = np.random.default_rng(118)
mp.mp.dps = 30
QV = mp.exp(-2*mp.pi)
QQ = float(QV)                                                     # q = e^-2pi


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
    _, s, vh = np.linalg.svd(M, full_matrices=M.shape[0] < M.shape[1])  # dünn, wenn möglich (große Zeilenzahl)
    rank = int(np.sum(s > tol*max(1.0, s[0]))) if s.size else 0
    return vh[rank:].conj()                                           # M @ v = 0 für jede Zeile v


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
    """C mit C g* C^-1 = s g für alle g (kleiner Faktor, SVD), wie z116; None, wenn es keine gibt."""
    d = G[0].shape[0]
    rows = [np.kron(g.conj().T, np.eye(d)) - s*np.kron(np.eye(d), g) for g in G]
    _, sv, vh = np.linalg.svd(np.vstack(rows))
    if sv[-1] > 1e-9:
        return None
    C = vh[-1].conj().reshape(d, d, order="F")
    return C/abs(np.linalg.det(C))**(1/d)


def kommutant(gens, d):
    """Basis des Kommutanten {W : g W = W g} (vec spaltenweise, kleiner Faktor)."""
    rows = np.vstack([np.kron(np.eye(d), g) - np.kron(g.T, np.eye(d)) for g in gens])
    return [v.reshape(d, d, order="F") for v in nullraum(rows)]


# ------------------------------------------------------------------ Twistordaten (wie z98)
I2 = np.eye(2)
SIGMA = [np.eye(2, dtype=complex), sx, sy, sz]
G4 = np.diag([1., 1., -1., -1.]).astype(complex)                     # CCR-Form auf T
S4 = np.block([[np.zeros((2, 2)), I2], [I2, np.zeros((2, 2))]]).astype(complex)   # tau auf T
uvec = [np.array([1, 0, -1, 0.]), np.array([0, 1, 0, -1.])]
vvec = [np.array([1, 0, 1, 0.]), np.array([0, 1, 0, 1.])]
Inp = (np.outer(uvec[0], uvec[1]) - np.outer(uvec[1], uvec[0])).astype(complex)       # I_pi (z89)
Ion = (np.outer(vvec[0], vvec[1]) - np.outer(vvec[1], vvec[0])).astype(complex)       # I_omega (z89)
NP = [-0.5*sum(s[i, j]*np.outer(uvec[i], vvec[j]) for i in range(2) for j in range(2)) for s in SIGMA]
NK = [-0.5*sum(s[i, j]*np.outer(vvec[i], uvec[j]) for i in range(2) for j in range(2)) for s in SIGMA]
NE = np.diag([-0.5, -0.5, 0.5, 0.5]).astype(complex)
MP = [1j*N for N in NP]
MK = [1j*N for N in NK]
ETA = np.diag([1., -1., -1., -1.])
MKL = [ETA[m, m]*MK[m] for m in range(4)]
E_T = (MP[0] + MKL[0])/2                                              # kompakte Energie (= i N_E, z98)
QM_T = (MP[0] - MKL[0])/2                                             # Kappengenerator Q^- (z104)
D_T = sla.expm(np.pi/2*E_T) @ QM_T @ sla.expm(-np.pi/2*E_T)           # Dilatation, z104: Ad(e^(pi E/2)) Q^- = D
LC = np.zeros((4,)*4)
for _p in itertools.permutations(range(4)):
    LC[_p] = round(np.linalg.det(np.eye(4)[list(_p)]))
BIV = []
for _i in range(4):
    for _j in range(_i + 1, 4):
        _X = np.zeros((4, 4), complex)
        _X[_i, _j], _X[_j, _i] = 1, -1
        BIV.append(_X)
GAM7 = np.diag([1.]*4 + [-1.]*4).astype(complex)                      # T gegen T*: Chiralität der Twistorspinoren


def gam(X, psi=0.0):
    """Klein-Clifford: X in Lambda^2 T wirkt auf T + T*: Z -> (1/2) eps(X, Z, .) in T*, W -> X W in T."""
    M = np.zeros((8, 8), complex)
    M[4:, :4] = 0.5*np.exp(1j*psi)*np.einsum("abcd,ab->dc", LC, X)
    M[:4, 4:] = X
    return M


def qf(X, Y, psi=0.0):
    G1, G2 = gam(X, psi), gam(Y, psi)
    return complex(((G1 @ G2 + G2 @ G1)/2)[0, 0])


def stern(X, psi=0.0):
    return 0.5*np.exp(1j*psi)*np.einsum("abcd,cd->ab", LC, X)


def quer(X):
    return G4 @ X.conj() @ G4


def reelle_basis(psi=0.0):
    """Lorentzsche Realität X-quer = *X (z98): reeller Unterraum R^(4,2) von Lambda^2 T."""
    rows = []
    for k in range(12):
        c = np.zeros(6, complex)
        c[k % 6] = 1 if k < 6 else 1j
        X = sum(ci*B for ci, B in zip(c, BIV))
        Dx = quer(X) - stern(X, psi)
        rows.append(np.concatenate([Dx.real.ravel(), Dx.imag.ravel()]))
    out = []
    for v in nullraum(np.array(rows).T):
        c = v[:6] + 1j*v[6:]
        out.append(sum(ci*B for ci, B in zip(c, BIV)))
    return out


GRID = np.linspace(0, 2*np.pi, 2881)[:-1]


def reell_phase(X, psi=0.0):
    r = np.array([np.abs(quer(np.exp(1j*p)*X) - stern(np.exp(1j*p)*X, psi)).max() for p in GRID])/np.abs(X).max()
    k = int(np.argmin(r))
    return float(GRID[k]), float(r[k])


def ebenen(M):
    """Die beiden zweidimensionalen Eigenebenen (aufsteigend nach Re, dann Im des Eigenwerts)."""
    w, V = np.linalg.eig(M)
    order = np.lexsort((np.round(w.imag, 8), np.round(w.real, 8)))
    w, V = w[order], V[:, order]
    return [(complex(w[0]), V[:, :2]), (complex(w[2]), V[:, 2:])]


def biv_von(P):
    a, b = P[:, 0], P[:, 1]
    return np.outer(a, b) - np.outer(b, a)


RB = reelle_basis()


def reell_koord(X):
    """X = sum c_k RB_k mit reellen c? (c, Residuum)"""
    A = np.array([np.concatenate([B.real.ravel(), B.imag.ravel()]) for B in RB]).T
    b = np.concatenate([X.real.ravel(), X.imag.ravel()])
    c = np.linalg.lstsq(A, b, rcond=None)[0]
    return c, float(np.linalg.norm(A @ c - b))


def quer_raum(vecs):
    """reelle Vektoren in R^(4,2), Q-orthogonal zu allen vecs; Rückgabe: Basis und Gram-Signatur."""
    rows = []
    for v in vecs:
        z = np.array([qf(B, v) for B in RB])
        rows += [z.real, z.imag]
    ns = nullraum(np.array(rows))
    basis = [sum(c*B for c, B in zip(v.real, RB)) for v in ns]
    Gm = np.array([[qf(X, Y).real for Y in basis] for X in basis]) if basis else np.zeros((0, 0))
    ew = np.linalg.eigvalsh((Gm + Gm.T)/2) if basis else np.array([])
    return basis, (int(np.sum(ew > 1e-9)), int(np.sum(ew < -1e-9)))


def orthonormal(basis):
    """reelle Orthonormalbasis (Q = +-1) der Spanne, positive zuerst."""
    Gm = np.array([[qf(X, Y).real for Y in basis] for X in basis])
    w, V = np.linalg.eigh((Gm + Gm.T)/2)
    out = [sum(V[k, j]*basis[k] for k in range(len(basis)))/np.sqrt(abs(w[j])) for j in range(len(w))]
    sig = np.sign(w)
    order = np.argsort(-sig, kind="stable")
    return [out[j] for j in order], [int(sig[j]) for j in order]


def chir4(frame):
    g = np.eye(8, dtype=complex)
    for f in frame:
        g = g @ gam(f)
    for ph in (1, 1j, -1, -1j):
        if np.allclose((ph*g) @ (ph*g), np.eye(8)):
            return ph*g
    raise ValueError


# =================================================================== Teil 1
say("Teil 1: D_perp nativ? Nullpaare, lichtartiger Schnitt, Uhr")

# --- 1.1 Klein-Clifford-Algebra aus Twistorbivektoren und das Nullpaar der Dilatation
ok_cl, Qf = True, np.zeros((6, 6), complex)
for a in range(6):
    for b in range(6):
        Ga, Gb = gam(BIV[a]), gam(BIV[b])
        ac = Ga @ Gb + Gb @ Ga
        ok_cl &= np.allclose(ac, ac[0, 0]*np.eye(8))
        Qf[a, b] = ac[0, 0]/2
psi_t = 0.7
ok_psi = all(abs(qf(BIV[a], BIV[b], psi_t) - np.exp(1j*psi_t)*Qf[a, b]) < 1e-12 for a in range(6) for b in range(6))
Gm_RB = np.array([[qf(X, Y) for Y in RB] for X in RB])
ew_RB = np.linalg.eigvalsh((Gm_RB.real + Gm_RB.real.T)/2)
sig_RB = (int(np.sum(ew_RB > 1e-9)), int(np.sum(ew_RB < -1e-9)))
fr6, sg6 = orthonormal(RB)
prod6 = np.eye(8, dtype=complex)
for f in fr6:
    prod6 = prod6 @ gam(f)
lam7, res7 = prop(prod6, GAM7)
ok_dil = all(prop(D_T @ MP[m] - MP[m] @ D_T, MP[m])[1] < 1e-12 and abs(prop(D_T @ MP[m] - MP[m] @ D_T, MP[m])[0] - 1) < 1e-12
             for m in range(4))
null_pi, null_om, q_po = qf(Inp, Inp), qf(Ion, Ion), qf(Inp, Ion)
ph_pi, r_pi = reell_phase(Inp)
ph_om, r_om = reell_phase(Ion)
nrm = np.sqrt(2*q_po)
gp8, gm8 = gam(Inp)/nrm, gam(Ion)/nrm                                 # gamma_+ <-> I_pi, gamma_- <-> I_omega
ok_null = np.allclose(gp8 @ gp8, 0) and np.allclose(gm8 @ gm8, 0) and np.allclose(gp8 @ gm8 + gm8 @ gp8, np.eye(8))
Pp8, Pm8 = gm8 @ gp8, gp8 @ gm8                                       # P_+ = gamma_- gamma_+ (wie z116)
G2_8 = Pp8 - Pm8
UU = np.array(uvec).T.astype(complex)
VV = np.array(vvec).T.astype(complex)
ok_proj_v = (np.allclose(Pp8[:4, :4] @ VV, VV) and np.allclose(Pp8[:4, :4] @ UU, 0) and np.allclose(Pm8[:4, :4] @ UU, UU)
             and np.allclose(Pm8[:4, :4] @ VV, 0))
ok_proj_u = (np.allclose(Pp8[:4, :4] @ UU, UU) and np.allclose(Pp8[:4, :4] @ VV, 0) and np.allclose(Pm8[:4, :4] @ VV, VV)
             and np.allclose(Pm8[:4, :4] @ UU, 0))
ok_proj = ok_proj_v or ok_proj_u
ebene_p, ebene_m = (("span v = Ebene von I_omega (tau = +1)", "span u = Ebene von I_pi (tau = -1)") if ok_proj_v else
                    ("span u = Ebene von I_pi (tau = -1)", "span v = Ebene von I_omega (tau = +1)"))
Pi_u, Pi_v = (Pm8[:4, :4], Pp8[:4, :4]) if ok_proj_v else (Pp8[:4, :4], Pm8[:4, :4])
ok_PK = all(np.allclose(Pi_u @ NP[m], NP[m]) and np.allclose(NP[m] @ Pi_u, 0) and np.allclose(Pi_v @ NK[m], NK[m])
            and np.allclose(NK[m] @ Pi_v, 0) for m in range(4))
lam_dil, res_dil = prop(G2_8[:4, :4], D_T)
tauU, tauV = S4 @ UU, S4 @ VV
tau_eig = (np.allclose(tauU, -UU), np.allclose(tauV, VV))
ZAHLEN["klein_clifford"] = {"Signatur reeller Unterraum (psi = 0)": sig_RB, "Gamma_7 = lam * Produkt der sechs": str(np.round(lam7, 9)),
                            "Q(I_pi, I_omega)": str(q_po), "reelle Phase I_pi, I_omega": [ph_pi, ph_om],
                            "Gamma_2 auf T = lam * D": str(np.round(lam_dil, 9)), "tau auf (span u, span v)": [-1, 1]}
check(f"[E] Klein-Clifford: Lambda^2 T mit der Form X ^ Y = Q(X, Y) eps wirkt auf T + T* (Z -> (1/2) eps(X, Z, .), "
      f"W -> X W) mit Gamma(X) Gamma(Y) + Gamma(Y) Gamma(X) = 2 Q(X, Y) für alle 36 Basispaare ({ok_cl}); die "
      f"Volumenformphase skaliert nur Q (Q_psi = e^(i psi) Q, {ok_psi}); die lorentzsche Realität X-quer = *X (z98) gibt "
      f"R^(4,2) mit Signatur {sig_RB}, und das Produkt der sechs Orthonormalvektoren ist {np.round(lam7, 6)} Gamma_7 "
      f"(Gamma_7 = +1 auf T, -1 auf T*; Residuum {res7:.0e}). Die nativen Bivektoren I_pi und I_omega sind null "
      f"(Q = {abs(null_pi):.0e}, {abs(null_om):.0e}), reell (Phase {ph_pi:.3f}, {ph_om:.3f}) und gepaart (Q(I_pi, I_omega) = "
      f"{q_po.real:.0f}): gamma_+ = Gamma(I_pi)/2, gamma_- = Gamma(I_omega)/2 erfüllen gamma^2 = 0, {{gamma_+, gamma_-}} = 1 "
      f"({ok_null}); P_+ = gamma_- gamma_+ ist auf T der Projektor auf {ebene_p}, P_- der auf {ebene_m} "
      f"({ok_proj}); Gamma_2 = P_+ - P_- ist auf T genau {np.round(lam_dil.real, 6)} D (Residuum {res_dil:.0e}) mit der "
      f"nativen Dilatation D = Ad(e^(pi E/2)) Q^- ([D, P_mu] = P_mu: {ok_dil}); P_mu bildet die Ebene von I_omega in die "
      f"von I_pi ab und annulliert diese, K_mu umgekehrt ({ok_PK}): der Nullschnitt der radialen su(2,2) durch P und K ist "
      f"dasselbe Paar. Das lichtartige Paar der siebzehnten Analyse ist also nativ das Nullpaar der Dilatation, die Punkte "
      f"0 und unendlich",
      ok_cl and ok_psi and sig_RB in ((2, 4), (4, 2)) and res7 < 1e-10 and abs(abs(lam7) - 1) < 1e-10 and ok_dil
      and abs(null_pi) < 1e-12 and abs(null_om) < 1e-12 and abs(q_po) > 1 and r_pi < 1e-12 and r_om < 1e-12 and ok_null
      and ok_proj and res_dil < 1e-12 and abs(abs(lam_dil) - 2) < 1e-12 and all(tau_eig) and ok_PK)

# --- 1.2 Die drei nativen Nullpaare der radialen sl(2,R)
paare = {}
for name, M in (("D (Dilatation)", D_T), ("Q^- (Kappe)", QM_T), ("E (Uhr)", E_T)):
    (w1, P1), (w2, P2) = ebenen(M)
    gram = [np.round(np.array([[np.vdot(x, G4 @ y) for y in P.T] for x in P.T]), 9) for P in (P1, P2)]
    isotrop = [np.allclose(g, 0) for g in gram]
    n1, n2 = biv_von(P1), biv_von(P2)
    (f1, r1), (f2, r2) = reell_phase(n1), reell_phase(n2)
    if r1 < 1e-9:
        n1 = np.exp(1j*f1)*n1
    if r2 < 1e-9:
        n2 = np.exp(1j*f2)*n2
    qb = qf(n1, n2)
    trans, sig_t = quer_raum([n1, n2])
    ebene, sig_e = quer_raum(trans) if len(trans) == 4 else ([], None)
    t1, t2 = S4 @ n1 @ S4.T, S4 @ n2 @ S4.T
    tau_fest = abs(prop(t1, n1)[1]) < 1e-9 and abs(prop(t2, n2)[1]) < 1e-9
    tau_tausch = abs(prop(t1, n2)[1]) < 1e-9 and abs(prop(t2, n1)[1]) < 1e-9
    paare[name] = {"Eigenwerte": [str(np.round(w1, 6)), str(np.round(w2, 6))], "isotrop": isotrop,
                   "Q(n,n)": [abs(qf(n1, n1)), abs(qf(n2, n2))], "Q(n+,n-)": str(np.round(qb, 6)),
                   "reell (Residuum)": [r1, r2], "Querraum-Signatur": sig_t, "Paarebene-Signatur": sig_e,
                   "tau": "fest" if tau_fest else ("tauscht" if tau_tausch else "mischt"), "n": (n1, n2), "trans": trans}
WP = sla.expm(1j*np.pi/2*D_T)
wick_p = prop(WP @ E_T @ np.linalg.inv(WP), QM_T)
WM = sla.expm(-1j*np.pi/2*D_T)
wick_m = prop(WM @ E_T @ np.linalg.inv(WM), QM_T)
ZAHLEN["nullpaare"] = {k: {kk: vv for kk, vv in v.items() if kk not in ("n", "trans")} for k, v in paare.items()}
ZAHLEN["wick"] = {"Ad(e^(+i pi D/2)) E = lam Q^-": str(np.round(wick_p[0], 9)), "Ad(e^(-i pi D/2)) E": str(np.round(wick_m[0], 9))}
pD, pQ, pE = paare["D (Dilatation)"], paare["Q^- (Kappe)"], paare["E (Uhr)"]
check(f"[E] Native Nullpaare der radialen sl(2,R) {{E, Q^-, D}} (Eigenebenen auf T, ihre Bivektoren in R^(4,2)): "
      f"D: Ebenen isotrop {pD['isotrop']}, Bivektoren null und reell (Residuum {max(pD['reell (Residuum)']):.0e}), Querraum "
      f"{pD['Querraum-Signatur']} (lorentzsch), tau {pD['tau']}; Q^-: isotrop {pQ['isotrop']}, reell bis auf die Phase "
      f"(Residuum {max(pQ['reell (Residuum)']):.0e}), Querraum {pQ['Querraum-Signatur']}, tau {pQ['tau']} (vertauscht die "
      f"Spitzen); E: Ebenen definit (Gram +-1, nicht isotrop {pE['isotrop']}), Bivektoren null, aber für keine Phase reell "
      f"(Residuum {min(pE['reell (Residuum)']):.2f}), ihre reelle Paarebene {pE['Paarebene-Signatur']}, Querraum "
      f"{pE['Querraum-Signatur']} (euklidisch), tau {pE['tau']}. Wick: Ad(e^(+-i pi D/2)) E = "
      f"{np.round(wick_p[0], 6)} Q^- bzw. {np.round(wick_m[0], 6)} Q^- (z104): Die Richtung der Wick-Abbildung legt fest, "
      f"auf welche Spitze die Originalebene fällt. Lichtartig (reell) sind nur die D- und die Q^--Paare (und die E-Bahn "
      f"cos theta Q^- + sin theta D dazwischen); das uhrgraduierte Paar ist komplex",
      all(pD["isotrop"]) and all(pQ["isotrop"]) and not any(pE["isotrop"]) and max(pD["reell (Residuum)"]) < 1e-9
      and max(pQ["reell (Residuum)"]) < 1e-9 and min(pE["reell (Residuum)"]) > 0.3
      and sorted(pD["Querraum-Signatur"]) == [1, 3] and sorted(pQ["Querraum-Signatur"]) == [1, 3]
      and 0 in pE["Querraum-Signatur"] and pD["tau"] == "fest" and pQ["tau"] == "tauscht" and pE["tau"] == "tauscht"
      and wick_p[1] < 1e-10 and wick_m[1] < 1e-10 and abs(wick_p[0] - 1j) < 1e-10 and abs(wick_m[0] + 1j) < 1e-10)

# --- 1.3 Chiralitätsidentität, Modul-U(1) und die Uhr auf S^+
fr4, sg4 = orthonormal(pD["trans"])
g5_tw = chir4(fr4)
c7, r7 = prop(GAM7, G2_8 @ g5_tw)
nQ1, nQ2 = pQ["n"]
cq = np.sqrt(2*qf(nQ1, nQ2))
gpQ, gmQ = gam(nQ1)/cq, gam(nQ2)/cq
G2Q = gmQ @ gpQ - gpQ @ gmQ
g5_Q = chir4(orthonormal(pQ["trans"])[0])
c7Q, r7Q = prop(GAM7, G2Q @ g5_Q)
th = 0.37
Uth = sla.expm(-1j*th*GAM7)                                          # Modul-U(1): e^(-i theta) auf T (z98, z103)
ok_U_kov = all(np.allclose(Uth @ gam(X) @ np.linalg.inv(Uth), gam(np.exp(-2j*th)*X, 4*th)) for X in BIV)
gerade = [gam(fr4[i]) @ gam(fr4[j]) for i, j in itertools.combinations(range(4), 2)] + [Pp8, Pm8, G2_8, g5_tw]
ok_U_ger = all(np.allclose(Uth @ X, X @ Uth) for X in gerade)
ZSU = sla.expm(-1j*np.pi/2*GAM7)                                     # Z_SU = U_(pi/2): -i auf T (z103)
# 256 Dimensionen wie z116
Gi = euklid(5)
G11 = chiral(Gi)
I32 = np.eye(32)
E6 = euklid(3)
hh = [E6[0], E6[1], E6[2], E6[3], 1j*E6[4], 1j*E6[5]]
h7 = chiral(hh)
G16 = [kr(x, I32) for x in hh] + [kr(h7, gi)*(1j if np.allclose(h7 @ h7, np.eye(8)) else 1.0) for gi in Gi]
Gt16 = chiral(G16)
Sp16 = (np.eye(256) + Gt16)/2
Sm16 = (np.eye(256) - Gt16)/2
gp16, gm16 = (G16[3] + G16[5])/2, (G16[3] - G16[5])/2
Pplus, Pminus = gm16 @ gp16, gp16 @ gm16
Gam2 = Pplus - Pminus
LI = (0, 1, 2, 4)
g5_16 = chiral([G16[a] for a in LI])
Ruhr = kr(np.eye(8), 1j*G11)
H7 = kr(h7, I32)
lamR_g, resR_g = prop(Ruhr @ Sp16, Gam2 @ g5_16 @ Sp16)
lamR_h, resR_h = prop(Ruhr @ Sp16, H7 @ Sp16)
lam_tot, res_tot = prop(Gt16, kr(h7, G11))
gp8z, gm8z = (hh[3] + hh[5])/2, (hh[3] - hh[5])/2
lam_h7, res_h7 = prop(h7, (gm8z @ gp8z - gp8z @ gm8z) @ chiral([hh[a] for a in LI]))
# Z_SU in der Basis von z116: Z_SU = -i auf T; T = Eigenraum h7 = +1 (Orientierung) bzw. -1
ZSU_a = -1j*h7
ZSU_b = 1j*h7
R_ist_ZSUinv = np.allclose((Ruhr @ Sp16), kr(np.linalg.inv(ZSU_a), I32) @ Sp16)
ZAHLEN["chiralitaet"] = {"Gamma_7 = c Gamma_2 gamma_5 (D-Paar)": str(np.round(c7, 9)), "(Q^--Paar)": str(np.round(c7Q, 9)),
                         "R|S+ = lam Gamma_2 gamma_5": str(np.round(lamR_g, 9)), "R|S+ = lam (h7 x 1)": str(np.round(lamR_h, 9)),
                         "Gamma_ges = lam h7 x Gamma_11": str(np.round(lam_tot, 9)), "h7 = lam Gamma_2 gamma_5 (8)": str(np.round(lam_h7, 9))}
check(f"[E] Chiralität: Auf T + T* ist Gamma_7 = {np.round(c7, 6)} Gamma_2 gamma_5 (D-Paar, Residuum {r7:.0e}) und "
      f"= {np.round(c7Q, 6)} Gamma_2 gamma_5 (Q^--Paar, {r7Q:.0e}) mit der Chiralität gamma_5 des jeweiligen Querraums: "
      f"Twistorchiralität = Lichtkegelgradierung mal 4D-Chiralität. Das Modul-U(1) U_theta = e^(-i theta Gamma_7) wirkt "
      f"kovariant (U Gamma_psi(X) U^-1 = Gamma_(psi + 4 theta)(e^(-2 i theta) X): Ladung 2 für I, 4 für eps, {ok_U_kov}) "
      f"und vertauscht mit allen geraden Daten (Lorentz-Bivektoren, P_+-, Gamma_2, gamma_5: {ok_U_ger}). In 256 Dimensionen "
      f"(z116) gilt Gamma_ges = {np.round(lam_tot, 3)} h7 x Gamma_11, h7 = {np.round(lam_h7, 3)} Gamma_2 gamma_5, und auf S^+ "
      f"R = i Gamma_11 = {np.round(lamR_g, 3)} Gamma_2 gamma_5 = {np.round(lamR_h, 3)} (h7 x 1) (Residuen {resR_g:.0e}, "
      f"{resR_h:.0e}): Die Uhr ist auf S^+ die Twistorchiralität, also Z_SU^-1 = U_(-pi/2) (Orientierung T = {{h7 = +1}}: "
      f"{R_ist_ZSUinv}; sonst Z_SU). Die 14D-Weyl-Bedingung leistet damit dieselbe Verklebung R = Z_SU^-1 wie der "
      f"E8-Kleber unter P_A3 (z97); als Symmetrie dreht die Uhr dann die Plebanski-Sektoren nur halb (z103)",
      r7 < 1e-10 and r7Q < 1e-10 and abs(abs(c7) - 1) < 1e-10 and ok_U_kov and ok_U_ger and resR_g < 1e-12
      and resR_h < 1e-12 and res_tot < 1e-12 and res_h7 < 1e-12 and abs(lamR_h - 1j) < 1e-12 and R_ist_ZSUinv)
laufzeit()

# --- 1.4 Uhrladungen und Indizes der nativen Schritte; Uhrkovarianz von D_perp


def A_rechteck(N, q, p=None):
    """A = S^+ p(S^4) auf C^N -> C^(N-1) (Hardy-Marken k >= 0); p = Koeffizienten von p(z), Standard 1 - q z."""
    p = [1.0, -q] if p is None else p
    S = np.zeros((N + 4*len(p) + 1, N + 4*len(p) + 1))
    for n in range(S.shape[0] - 1):
        S[n + 1, n] = 1
    P = sum(c*np.linalg.matrix_power(S, 4*j) for j, c in enumerate(p))
    A = S.T @ P
    return A[:N - 1, :N]


def index_von(A, tol=1e-12):
    sv = np.linalg.svd(A, compute_uv=False)
    r = int(np.sum(sv > tol))
    return A.shape[1] - r, A.shape[0] - r


def rel_spur(A, W1=None, W2=None, s=1.0):
    e1 = sla.expm(-s*A.conj().T @ A)
    e2 = sla.expm(-s*A @ A.conj().T)
    W1 = np.eye(A.shape[1]) if W1 is None else W1
    W2 = np.eye(A.shape[0]) if W2 is None else W2
    return complex(np.trace(W1 @ e1) - np.trace(W2 @ e2))


N = 25
Aq = A_rechteck(N, QQ)
Rs = np.diag([1j**k for k in range(N)])
Rt = np.diag([1j**k for k in range(N - 1)])
ladung = prop(Rt @ Aq @ np.linalg.inv(Rs), Aq)
diag_off = sorted({j - i for i, j in zip(*np.nonzero(np.abs(Aq) > 1e-15))})
Sdown = np.zeros((N - 1, N))                                          # S^+: |k> -> |k-1>
for k in range(1, N):
    Sdown[k - 1, k] = 1
S4down = np.zeros((N - 4, N))
for k in range(4, N):
    S4down[k - 4, k] = 1
S4up = np.zeros((N + 4, N))
for k in range(N):
    S4up[k + 4, k] = 1
schritte = {"Viertelschritt S^+": Sdown, "A_q = S^+(1 - q S^4)": Aq, "Transfer abwärts S^+4": S4down, "Transfer J_-1 = S^4": S4up}
ind_tab, spur_tab, uhr_tab = {}, {}, {}
for name, A in schritte.items():
    k_, c_ = index_von(A)
    ind_tab[name] = k_ - c_
    spur_tab[name] = [round(rel_spur(A, s=s_).real, 10) for s_ in (0.3, 1.0, 4.0)]
    Rsrc = np.diag([1j**k for k in range(A.shape[1])])
    Rtrg = np.diag([1j**k for k in range(A.shape[0])])
    uhr_tab[name] = complex(np.round(prop(Rtrg @ A @ np.linalg.inv(Rsrc), A)[0], 9))
# Möbius-Senker L_1 auf D_h (Tiefstgewichtsmodul, z83/z91): <n-1|L_1|n> = sqrt(n (n + 2h - 1))
moeb = {}
for h_ in (0.5, 1.0, 2.0):
    L1 = np.zeros((N - 1, N))
    for n in range(1, N):
        L1[n - 1, n] = np.sqrt(n*(n + 2*h_ - 1))
    k_, c_ = index_von(L1, 1e-9)
    moeb[h_] = (k_ - c_, [round(rel_spur(L1/np.sqrt(N), s=s_).real, 10) for s_ in (0.3, 1.0, 4.0)])
# Uhrkovarianz von D_perp (Twistorbasis, quadratische Abschneidung N = 13)
Nq = 13
Sq = np.zeros((Nq, Nq))
for n in range(Nq - 1):
    Sq[n + 1, n] = 1
Aquad = Sq.T - QQ*np.linalg.matrix_power(Sq, 3)
Rq = np.diag([1j**k for k in range(Nq)])
Dperp = np.kron(gm8, Aquad) + np.kron(gp8, Aquad.conj().T)
D2 = Dperp @ Dperp
ok_D2 = np.allclose(D2, np.kron(Pp8, Aquad @ Aquad.conj().T) + np.kron(Pm8, Aquad.conj().T @ Aquad))
tau8 = sla.block_diag(S4, S4)
lor8 = sla.expm(0.6*gam(fr4[1]) @ gam(fr4[2]))
wick8 = sla.expm(1j*np.pi/4*G2_8)
kandidaten = {"native Uhr auf T (= 1, z103)": np.eye(8), "Z_SU = U_(pi/2)": ZSU, "U_(pi/4)": sla.expm(-1j*np.pi/4*GAM7),
              "tau": tau8, "Lorentzdrehung": lor8, "Wick e^(i pi Gamma_2/4)": wick8}
kov = {}
for name, X in kandidaten.items():
    Rt_ = np.kron(X, Rq)
    Ri_ = np.linalg.inv(Rt_)
    kov[name] = (float(np.abs(Rt_ @ Dperp @ Ri_ - Dperp).max()), float(np.abs(Rt_ @ D2 @ Ri_ - D2).max()))
ok_tau_cl = np.allclose(tau8 @ gam(Inp) @ np.linalg.inv(tau8), gam(S4 @ Inp @ S4.T))
L1q = np.zeros((Nq, Nq))
for n in range(1, Nq):
    L1q[n - 1, n] = np.sqrt(n*n)
D_moeb = np.kron(gm8, L1q) + np.kron(gp8, L1q.T)
RL = np.diag([np.exp(2j*np.pi*(1 + n)) for n in range(Nq)])                 # R = e^(2 pi i L_0) auf D_1
moeb_kov = float(np.abs(np.kron(np.eye(8), RL) @ D_moeb @ np.linalg.inv(np.kron(np.eye(8), RL)) - D_moeb).max())
# L_0-Kovarianz: H = H_s (x) 1 + 1 (x) L_0, H_s = -i E auf T + T* (auf T*: -E^T)
E8 = sla.block_diag(E_T, -E_T.T)
Hs = -1j*E8
ok_E8 = np.allclose(E8 @ gam(Inp) - gam(Inp) @ E8, gam(E_T @ Inp + Inp @ E_T.T))
nE1, nE2 = paare["E (Uhr)"]["n"]
cE = np.sqrt(2*qf(nE1, nE2))
gE1, gE2 = gam(nE1)/cE, gam(nE2)/cE
ok_nullE = np.allclose(gE1 @ gE1, 0) and np.allclose(gE2 @ gE2, 0) and np.allclose(gE1 @ gE2 + gE2 @ gE1, np.eye(8))
ladE = [prop(Hs @ g - g @ Hs, g) for g in (gE1, gE2)]
ladD = [prop(Hs @ g - g @ Hs, g) for g in (gp8, gm8)]
gEp, gEm = (gE1, gE2) if ladE[0][0].real > 0 else (gE2, gE1)            # E-Ladung +1 bzw. -1
L0D = np.diag([1.0 + n for n in range(Nq)])                              # L_0 auf D_1
Htot = np.kron(Hs, np.eye(Nq)) + np.kron(np.eye(8), L0D)
D_moebE = np.kron(gEp, L1q) + np.kron(gEm, L1q.T)
kovE = {"Möbius, Uhrpaar E": float(np.abs(Htot @ D_moebE - D_moebE @ Htot).max()),
        "Möbius, Dilatationspaar": float(np.abs(Htot @ D_moeb - D_moeb @ Htot).max())}
HtotA = np.kron(Hs, np.eye(Nq)) + np.kron(np.eye(8), np.diag([k/4 for k in range(Nq)]))   # L_0 = k/4 auf den Marken
for lab, (ga, gb) in (("A_q, Uhrpaar E", (gEp, gEm)), ("A_q, Uhrpaar E vertauscht", (gEm, gEp)),
                      ("A_q, Dilatationspaar", (gm8, gp8))):
    Dx = np.kron(ga, Aquad) + np.kron(gb, Aquad.conj().T)
    kovE[lab] = float(np.abs(HtotA @ Dx - Dx @ HtotA).max())
ZAHLEN["schritte"] = {"Index": ind_tab, "relative Spur s = 0,3/1/4": spur_tab, "Uhrladung R A R^-1 = lam A": uhr_tab,
                      "A_q: Nebendiagonalen (Markengrade)": diag_off, "Möbius L_1 auf D_h: (Index, Spuren)": {str(k): v for k, v in moeb.items()},
                      "Uhrkovarianz D_perp (|R D R^-1 - D|, |R D^2 R^-1 - D^2|)": kov, "Möbius-D_perp uhrinvariant": moeb_kov,
                      "E-Ladungen Uhrpaar": [str(np.round(x[0], 9)) for x in ladE], "Dilatationspaar E-Eigenoperator (Residuen)": [x[1] for x in ladD],
                      "L_0-Kovarianz |[H, D]|": kovE}
nur_wick = all((v[0] > 1e-3) == (k != "Wick e^(i pi Gamma_2/4)") for k, v in kov.items())
check(f"[E] Uhrladungen und Indizes auf den Hardy-Marken (N = {N}): R A_q R^-1 = {np.round(ladung[0], 6)} A_q (Residuum "
      f"{ladung[1]:.0e}), A_q hat die Markengrade {sorted(int(-x) for x in diag_off)} (S^+ und q S^3, beide mit Uhrphase -i): "
      f"Z_4-homogen, nicht L_0-homogen. Indizes {ind_tab}, relative Spuren s-unabhängig {spur_tab}; Uhrladungen {uhr_tab}: "
      f"Index +1 trägt allein der Viertelschritt (fermionisch, Ladung -i), der uhrneutrale Transfer J_-1 = S^4 hat Index "
      f"-+4. Der Möbiussenker L_1 auf D_h (h = 1/2, 1, 2) ist uhrneutral, L_0-homogen und hat Index +1 je Tiefstgewicht "
      f"({ {str(k): v[0] for k, v in moeb.items()} }). D_perp = gamma_- A_q + gamma_+ A_q^+ (Twistorbasis) gibt D^2 = "
      f"P_+ A A^+ + P_- A^+ A ({ok_D2}); unter R_ges = X (x) R_Marken ist D^2 für jede native Wahl von X invariant, D_perp "
      f"aber nur mit dem imaginären Dilatationselement X = e^(i pi Gamma_2/4) (auf T = e^(-i pi D/2), eine Wick-Abbildung, z104; "
      f"|R D R^-1 - D| = "
      f"{ {k: f'{v[0]:.2f}' for k, v in kov.items()} }); tau wirkt Clifford-kovariant ({ok_tau_cl}). Mit L_1 statt A_q ist "
      f"D_perp uhrinvariant (Residuum {moeb_kov:.0e}). Energiekovarianz [H_s + L_0, D] = 0 (H_s = -i E auf T + T*, "
      f"Clifford-kovariant: {ok_E8}): das Uhrpaar trägt die E-Ladungen {[str(np.round(x[0].real, 6)) for x in ladE]} "
      f"(gamma^2 = 0, {{gamma, gamma'}} = 1: {ok_nullE}), das Dilatationspaar ist kein E-Eigenoperator (Residuen "
      f"{[f'{x[1]:.2f}' for x in ladD]}); |[H, D]| = { {k: f'{v:.2e}' for k, v in kovE.items()} }: kovariant ist allein der "
      f"Möbiussenker am komplexen Uhrpaar (euklidischer Querraum), A_q mit keinem Paar",
      ladung[1] < 1e-14 and abs(ladung[0] + 1j) < 1e-14 and diag_off == [-3, 1] and ind_tab["Viertelschritt S^+"] == 1
      and ind_tab["A_q = S^+(1 - q S^4)"] == 1 and ind_tab["Transfer abwärts S^+4"] == 4 and ind_tab["Transfer J_-1 = S^4"] == -4
      and abs(uhr_tab["Viertelschritt S^+"] + 1j) < 1e-9 and abs(uhr_tab["A_q = S^+(1 - q S^4)"] + 1j) < 1e-9
      and abs(uhr_tab["Transfer abwärts S^+4"] - 1) < 1e-9 and abs(uhr_tab["Transfer J_-1 = S^4"] - 1) < 1e-9
      and all(max(abs(x - v[0]) for x in v) < 1e-9 for v in spur_tab.values()) and all(v[0] == 1 for v in moeb.values())
      and ok_D2 and nur_wick and all(v[1] < 1e-12 for v in kov.values()) and ok_tau_cl and moeb_kov < 1e-12
      and ok_E8 and ok_nullE and all(x[1] < 1e-12 for x in ladE) and sorted(round(x[0].real, 9) for x in ladE) == [-1, 1]
      and all(x[1] > 1e-3 for x in ladD) and kovE["Möbius, Uhrpaar E"] < 1e-12 and kovE["Möbius, Dilatationspaar"] > 1e-3
      and min(v for k, v in kovE.items() if k.startswith("A_q")) > 1e-3)
befund("E/C: D_perp ist kein natives Operatordatum. Nativ sind das lichtartige Paar (Nullpaar der Dilatation: Ebenen von "
       "I_pi und I_omega, Punkte 0 und unendlich; daneben das Kappenpaar, beide reell, aus einer E-Bahn, durch tau "
       "unterschieden) und der Toeplitzindex +1. Der Index sitzt aber im Viertelschritt (Uhrladung -i), während jede native "
       "Struktur auf der Spinorseite uhrneutral ist; D_perp hat daher keine Uhrladung, nur seine Uhrbahn und D^2 sind nativ "
       "invariant. Der Kompensator wäre das imaginäre Dilatationselement e^(i pi Gamma_2/4) (Wick-Abbildung), keine "
       "Symmetrie. Uhrneutrale native Schritte haben "
       "Index -+4 (Transfer) oder koppeln als Möbiussenker an das uhrgraduierte, komplexe (euklidische) Paar. Konstruierbar "
       "ist D^2 = P_- A^+A + P_+ AA^+ aus nativen Stücken bis auf die Zuordnung des Nullmodensektors", True)
laufzeit()

# =================================================================== Teil 2
say("\nTeil 2: Die Wirkung -- zulässige Gewichte der relativen Spur")

# --- 2.1 Kommutant und Endlichkeit
LOR8 = [hh[a] @ hh[b] for a, b in itertools.combinations(LI, 2)]
K8 = kommutant(LOR8, 8)
gens32 = []
for _ in range(4):
    cmat = rng.normal(size=(10, 10))
    gens32.append(sum(cmat[i, j]*(Gi[i] @ Gi[j]) for i in range(10) for j in range(i + 1, 10)))
K32 = kommutant(gens32, 32)
Wb = [kr(A, B) for A in K8 for B in K32]
spalten = [np.concatenate([(W @ Pplus - Pplus @ W).ravel(), (W @ Gt16 - Gt16 @ W).ravel()]) for W in Wb]
ns_w = nullraum(np.array(spalten).T)
Wz = [sum(c*W for c, W in zip(v, Wb)) for v in ns_w]
Wsp = [W @ Sp16 for W in Wz]
Mw = np.array([W.ravel() for W in Wsp])
_, sv_w, vh_w = np.linalg.svd(Mw @ Mw.conj().T)
rang_w = int(np.sum(sv_w > 1e-9*sv_w[0]))
# Basis der Einschränkungen auf S^+
_, _, vh_r = np.linalg.svd(Mw, full_matrices=False)
Wr = [vh_r[k].reshape(256, 256) for k in range(rang_w)]
lorpaare = [(a, b) for a, b in itertools.combinations(LI, 2)]
GABCD = {}
for (a, b) in lorpaare:
    for (c, d) in lorpaare:
        GABCD[(a, b, c, d)] = G16[a] @ G16[b] @ G16[c] @ G16[d]
PpS, PmS = Pplus @ Sp16, Pminus @ Sp16


def spuren(W, P):
    return np.array([np.sum(W.T*(X @ P)) for X in GABCD.values()])


Tpl = np.array([spuren(W, PpS) for W in Wr])
Tmi = np.array([spuren(W, PmS) for W in Wr])
ns_end = nullraum((Tpl + Tmi).T)
Wend = [sum(c*W for c, W in zip(v, Wr)) for v in ns_end]
kandW = {"Gamma_2": Gam2 @ Sp16, "R": Ruhr @ Sp16, "1": Sp16, "gamma_5": g5_16 @ Sp16}
im_end = {}
for k, W in kandW.items():
    A_ = np.array([X.ravel() for X in Wend]).T
    x_ = np.linalg.lstsq(A_, W.ravel(), rcond=None)[0]
    im_end[k] = float(np.linalg.norm(A_ @ x_ - W.ravel())/np.linalg.norm(W.ravel())) < 1e-9
A4 = np.array([X.ravel() for X in Wr]).T
im_vier = all(float(np.linalg.norm(A4 @ np.linalg.lstsq(A4, W.ravel(), rcond=None)[0] - W.ravel())) < 1e-8*np.linalg.norm(W.ravel())
              for W in kandW.values())
EPS4 = np.zeros((4, 4, 4, 4))
for p_ in itertools.permutations(range(4)):
    EPS4[p_] = round(np.linalg.det(np.eye(4)[list(p_)]))
pos = {a: i for i, a in enumerate(LI)}
ETA4 = np.diag([1., 1., 1., -1.])                                    # gamma_a^2 auf (0, 1, 2, 4), z116
TR = {k: np.sum((Ruhr @ Sp16).T*(X @ Pminus)) for k, X in GABCD.items()}
TG = {k: np.sum((Gam2 @ Sp16).T*(X @ Pminus)) for k, X in GABCD.items()}
cR = {k: v/EPS4[pos[k[0]], pos[k[1]], pos[k[2]], pos[k[3]]] for k, v in TR.items() if EPS4[pos[k[0]], pos[k[1]], pos[k[2]], pos[k[3]]] != 0}
met = {k: ETA4[pos[k[0]], pos[k[3]]]*ETA4[pos[k[1]], pos[k[2]]]*(k[0] == k[3])*(k[1] == k[2])
       - ETA4[pos[k[0]], pos[k[2]]]*ETA4[pos[k[1]], pos[k[3]]]*(k[0] == k[2])*(k[1] == k[3]) for k in GABCD}
cG = {k: TG[k]/met[k] for k in GABCD if met[k] != 0}
R_eps = len(set(np.round(list(cR.values()), 9))) == 1 and all(abs(TR[k]) < 1e-9 for k in GABCD
                                                               if EPS4[pos[k[0]], pos[k[1]], pos[k[2]], pos[k[3]]] == 0)
G_met = len(set(np.round(list(cG.values()), 9))) == 1 and all(abs(TG[k]) < 1e-9 for k in GABCD if met[k] == 0)
wR, wG = list(cR.values())[0], list(cG.values())[0]
# Dualität: Tr(R B Y) gegen Tr(Gamma_2 (*B) Y) mit (*gamma)_ab = sum_(e<f) eps_abef eta^ee eta^ff gamma_ef
dual_vals = {(a, b, c, d): sum(EPS4[pos[a], pos[b], pos[e_], pos[f_]]*ETA4[pos[e_], pos[e_]]*ETA4[pos[f_], pos[f_]]*TG[(e_, f_, c, d)]
                               for (e_, f_) in lorpaare) for (a, b) in lorpaare for (c, d) in lorpaare}
kap_d = {k: TR[k]/v for k, v in dual_vals.items() if abs(v) > 1e-9}
dual_ok = (len(set(np.round(list(kap_d.values()), 9))) == 1
           and all(abs(TR[k]) < 1e-9 for k, v in dual_vals.items() if abs(v) <= 1e-9))
kap_dual = list(kap_d.values())[0]
ZAHLEN["gewichte"] = {"dim Kommutant (Lorentz, 8)": len(K8), "dim Kommutant (Spin(10), 32)": len(K32),
                      "dim mit [W, P_+] = [W, Gamma_ges] = 0": len(Wz), "Rang auf S^+": rang_w, "dim endlich": len(Wend),
                      "in der endlichen Spanne": im_end, "R: c eps (P_-)": str(np.round(wR, 6)), "Gamma_2: c Metrik (P_-)": str(np.round(wG, 6)),
                      "Tr(R B Y) = kappa Tr(Gamma_2 *B Y)": str(np.round(kap_dual, 9))}
check(f"[E] Zulässige Gewichte: Kommutant der Lorentz-Bivektoren auf dem konformen Faktor dim {len(K8)}, der Spin(10) "
      f"auf 32 dim {len(K32)} (1, Gamma_11); mit [W, P_+] = [W, Gamma_ges] = 0 bleiben {len(Wz)}, auf S^+ Rang {rang_w}: "
      f"genau span{{1, Gamma_2, gamma_5, Gamma_2 gamma_5}} ({im_vier}). Die relative Spur Tr(W gamma_ab gamma_cd e^(-s D^2)) "
      f"ist genau dann endlich (Aufhebung zwischen P_+ und P_-), wenn W in einem {len(Wend)}-dimensionalen Raum liegt: "
      f"{ {k: v for k, v in im_end.items()} } -- also span{{Gamma_2, R}}, R = -i Gamma_2 gamma_5. Werte bei Index 1: "
      f"Gewicht R gibt {np.round(wR.real, 6)} eps_abcd ({R_eps}), Gewicht Gamma_2 (die kanonische McKean-Singer-Gradierung "
      f"von D_perp) {np.round(wG.real, 6)} (eta_ad eta_bc - eta_ac eta_bd) ({G_met}): die Indexgradierung allein gibt den "
      f"topologischen Sektor, die Uhr den Einstein-Sektor, W = a Gamma_2 + b R beide (Holst-Mischung a/b). Dualität: "
      f"Tr(R B Y) = {np.round(kap_dual.real, 6)} Tr(Gamma_2 (*B) Y) für alle 36 Paare von Lorentz-Bivektoren ({dual_ok}): "
      f"Gewicht R mit Sektor B ist Gewicht Gamma_2 mit Sektor *B. Die Gewichtswahl ist also genau die Paarwahl "
      f"{{+-e^e}} gegen {{+-*(e^e)}} aus z103 (P_vol), in Spurform",
      len(K32) == 2 and rang_w == 4 and im_vier and len(Wend) == 2 and im_end["Gamma_2"] and im_end["R"]
      and not im_end["1"] and not im_end["gamma_5"] and R_eps and G_met and abs(abs(wR) - 64) < 1e-9 and abs(abs(wG) - 64) < 1e-9
      and dual_ok and abs(abs(kap_dual) - 1) < 1e-9)

# --- 2.2 Die native Uhr wirkt auch auf den Marken
uhr_beide = {}
for qq_, lab in ((QQ, "q = e^-2pi"), (0.5, "q = 1/2")):
    for Nn in (41, 81):
        A = A_rechteck(Nn, qq_)
        Rs_ = np.diag([1j**k for k in range(Nn)])
        Rt_ = np.diag([1j**k for k in range(Nn - 1)])
        uhr_beide[f"{lab}, N = {Nn}"] = [rel_spur(A, Rs_, Rt_, s_) for s_ in (0.1, 0.3, 1.0, 4.0, 30.0)]
dev = {k: [v_ - 1 for v_ in v] for k, v in uhr_beide.items()}
s_abh = all(max(abs(x) for x in v[:4]) > 1e-8 for v in dev.values())
n_unabh = all(max(abs(a_ - b_) for a_, b_ in zip(dev[f"{lab}, N = 41"], dev[f"{lab}, N = 81"])) < 1e-6
              for lab in ("q = e^-2pi", "q = 1/2"))
form_1pi = all(abs(x.real - x.imag) < 1e-12 for v in dev.values() for x in v)
grenze = all(abs(v[-1]) < 1e-9 for k, v in dev.items() if "e^-2pi" in k)
dmax = {lab: max(abs(x) for x in dev[f"{lab}, N = 81"]) for lab in ("q = e^-2pi", "q = 1/2")}
ZAHLEN["uhr_auf_marken"] = {k: [str(np.round(x, 10)) for x in v] for k, v in dev.items()}
check(f"[E, negativ] Die native Uhr R = e^(2 pi i L_0) wirkt auch auf den Marken (R|k> = i^k |k>, z48, z82). Mit R auf "
      f"beiden Faktoren ist die relative Spur Tr(R_ges gamma_abcd e^(-s D^2)) = 64 eps f(s) mit f(s) = Tr(R e^(-s A^+A)) - "
      f"Tr(R e^(-s AA^+)) = 1 - (1 + i) delta(s) (Form (1+i): {form_1pi}): s-abhängig ({s_abh}), unabhängig von der "
      f"Abschneidung ({n_unabh}), komplex; |f - 1| bis {dmax['q = e^-2pi']:.2e} bei q = e^-2pi und {dmax['q = 1/2']:.3f} bei "
      f"q = 1/2, f -> 1 erst für s -> unendlich ({grenze}). Grund: A_q verschiebt die Uhrklasse, die Nichtnullspektren der "
      f"vier Klassen heben sich nur bis auf Randterme auf. Topologisch (Index) ist die Spur mit Uhr nur, wenn R allein auf "
      f"dem inneren Faktor wirkt",
      s_abh and n_unabh and form_1pi and grenze and dmax["q = e^-2pi"] > 1e-7)

# --- 2.3 MacDowell-Mansouri mit nativen Daten
TfullR = np.zeros((4, 4, 4, 4), complex)
TfullG = np.zeros((4, 4, 4, 4), complex)
for a, b, c, d in itertools.product(range(4), repeat=4):
    if a == b or c == d:
        continue
    X = G16[LI[a]] @ G16[LI[b]] @ G16[LI[c]] @ G16[LI[d]] @ Pminus
    TfullR[a, b, c, d] = np.sum((Ruhr @ Sp16).T*X)
    TfullG[a, b, c, d] = np.sum((Gam2 @ Sp16).T*X)


def kn(S1, S2):
    return (np.einsum("ac,bd->abcd", S1, S2) + np.einsum("bd,ac->abcd", S1, S2) - np.einsum("ad,bc->abcd", S1, S2)
            - np.einsum("bc,ad->abcd", S1, S2))


def dens(Af, Bf):
    return np.einsum("mnrs,mn,rs->", LC, Af, Bf)/4


def lagr(T, Rf, EE, lam):
    F = Rf + lam*EE
    return sum(T[a, b, c, d]*dens(F[a, b], F[c, d]) for a, b, c, d in itertools.product(range(4), repeat=4)
               if T[a, b, c, d] != 0)/4


mm = []
for _ in range(6):
    Rabcd = sum(kn(S, S)/2 for S in [(lambda M: M + M.T)(rng.normal(size=(4, 4))) for _ in range(3)])
    Rup = np.einsum("ai,bj,ijcd->abcd", ETA4, ETA4, Rabcd)                 # R^ab_cd
    zykl = float(np.abs(Rabcd + np.einsum("acdb->abcd", Rabcd) + np.einsum("adbc->abcd", Rabcd)).max())
    Ric = np.einsum("ac,abcd->bd", ETA4, Rabcd)
    Rs_sk = float(np.einsum("bd,bd->", ETA4, Ric))
    Ric_up = ETA4 @ Ric @ ETA4
    Eul = float(np.einsum("abcd,ijkl,ai,bj,ck,dl->", Rabcd, Rabcd, ETA4, ETA4, ETA4, ETA4)
                - 4*np.einsum("ab,ab->", Ric, Ric_up) + Rs_sk**2)            # Gauss-Bonnet-Kombination
    e = rng.normal(size=(4, 4)) + 2*np.eye(4)                               # e^a_mu
    dete = float(np.linalg.det(e))
    Rf = np.einsum("abcd,cm,dn->abmn", Rup, e, e)
    EE = np.einsum("am,bn->abmn", e, e) - np.einsum("an,bm->abmn", e, e)
    out = {}
    for lab, T in (("R", TfullR), ("G", TfullG)):
        L0, Lp, Lm = lagr(T, Rf, EE, 0.0), lagr(T, Rf, EE, 1.0), lagr(T, Rf, EE, -1.0)
        out[lab] = (L0, (Lp - Lm)/2, (Lp + Lm)/2 - L0)
    mm.append({"Rsk": Rs_sk, "E": Eul, "det": dete, "zykl": zykl, **out})
c1 = [m_["R"][1]/(m_["Rsk"]*m_["det"]) for m_ in mm]
c2 = [m_["R"][2]/m_["det"] for m_ in mm]
c0 = [m_["R"][0]/(m_["E"]*m_["det"]) for m_ in mm]
skala = max(abs(m_["R"][1]) for m_ in mm)
holst = max(abs(m_["G"][1]) for m_ in mm)/skala
volG = max(abs(m_["G"][2]) for m_ in mm)/skala
konst = lambda xs: max(abs(x - xs[0]) for x in xs)/abs(xs[0])
ratio = sp.nsimplify(round((c2[0]/c1[0]).real, 9))
lam_s, sig_s, l_s, G_s, Rv = sp.symbols("lambda sigma l G R_Vak", real=True)
kap = sp.Symbol("kappa")
# S = kappa (c0 Euler + c1 lam R + c2 lam^2) vol  <->  (1/16 pi G)(R - 2 Lambda) vol
Lam_s = sp.simplify(-ratio*lam_s/2)
Lam_nativ = sp.simplify(Lam_s.subs(lam_s, sig_s/l_s**2).subs({sig_s: 1, l_s: Rv/2}))
b_eul = sp.simplify(sp.nsimplify(round((c0[0]/c1[0]).real, 9))*l_s**2/(16*sp.pi*G_s))
ZAHLEN["mm"] = {"c1 = L1/(R det e)": [str(np.round(x, 9)) for x in c1[:2]], "c2/c1": str(ratio),
                "c0/c1 (Euler)": str(np.round(c0[0]/c1[0], 9)), "Lambda": str(Lam_s), "Lambda nativ (sigma = 1, l = R_Vak/2)": str(Lam_nativ),
                "Euler-Gewicht": str(b_eul), "Holst/Skala": holst, "Volumen Gamma_2/Skala": volG,
                "zyklische Identität": max(m_["zykl"] for m_ in mm)}
check(f"[E] MacDowell-Mansouri in der Spur (sechs zufällige algebraische Krümmungstensoren mit zyklischer Identität, "
      f"Rest {max(m_['zykl'] for m_ in mm):.0e}, und zufällige Vierbeine): Gewicht R: Tr_rel(R F ^ F) mit F = R + lambda "
      f"e ^ e hat lambda^0 = c0 Euler det e (c0 konstant, Streuung {konst(c0):.0e}), lambda^1 = c1 R det e ({konst(c1):.0e}), "
      f"lambda^2 = c2 det e ({konst(c2):.0e}), c2/c1 = {ratio}: die gemischte Antwort ist rein Einstein, aber beim nativen "
      f"lambda = sigma/l^2 trägt die Wirkung Euler-Gewicht {b_eul} und Lambda = {Lam_s} = {Lam_nativ} (sigma = +1, "
      f"l = R_Vak/2, z98): AdS, wie z104s Regel (positive Zeit => Lambda <= 0); dS bräuchte P_cap. Gewicht Gamma_2: "
      f"lambda^1 = Holst = {holst:.0e}, lambda^2 = {volG:.0e} (relativ), nur Pontrjagin: topologisch, ohne Einstein und "
      f"Lambda",
      konst(c0) < 1e-8 and konst(c1) < 1e-8 and konst(c2) < 1e-8 and ratio == 6 and holst < 1e-10 and volG < 1e-10
      and Lam_nativ == -12/Rv**2 and max(m_["zykl"] for m_ in mm) < 1e-12)
laufzeit()

# --- 2.4 Vorzeichenbits
basis_w = complex(np.sum((Ruhr @ Sp16).T*(GABCD[(0, 1, 2, 4)] @ Pminus)))
var = {"R -> R^-1": complex(np.sum((np.linalg.inv(Ruhr) @ Sp16).T*(GABCD[(0, 1, 2, 4)] @ Pminus))),
       "S^+ -> S^-": complex(np.sum((Ruhr @ Sm16).T*(GABCD[(0, 1, 2, 4)] @ Pminus))),
       "Nullmode in P_+ (gamma_+ A + gamma_- A^+, Index -1)": complex(np.sum((Ruhr @ Sp16).T*(GABCD[(0, 1, 2, 4)] @ Pplus))),
       "4D-Orientierung (gamma_0 <-> gamma_1)": complex(np.sum((Ruhr @ Sp16).T*(G16[1] @ G16[0] @ G16[2] @ G16[4] @ Pminus)))}
alle_kippen = all(abs(v + basis_w) < 1e-9 for v in var.values())
g5_aus = prop(g5_16 @ Sp16, 1j*Gam2 @ Ruhr @ Sp16)
C10 = {+1: np.eye(32, dtype=complex), -1: np.eye(32, dtype=complex)}   # explizit: Produkt der reellen bzw. imaginären
for j in range(5):                                                      # Majoranas von Cl(10) (statt der SVD aus z116)
    C10[+1] = C10[+1] @ Gi[2*j]
    C10[-1] = C10[-1] @ Gi[2*j + 1]
ok_C10 = all(np.allclose(C10[s] @ g.conj() @ np.linalg.inv(C10[s]), s*g) for s in (1, -1) for g in Gi)
C16 = None
for s6 in (1, -1):
    for s10 in (1, -1):
        Ca, Cb = konj(hh, s6), C10[s10]
        if Ca is None or Cb is None:
            continue
        C = np.kron(Ca, Cb)
        for s_ in (1, -1):
            if all(np.allclose(C @ g.conj() @ np.linalg.inv(C), s_*g) for g in G16):
                if np.allclose(C @ C.conj(), np.eye(256)) and np.allclose(C @ Gt16.conj(), Gt16 @ C):
                    C16 = C
J_G2 = C16 is not None and np.allclose(C16 @ Gam2.conj(), Gam2 @ C16)
J_R = C16 is not None and np.allclose(C16 @ Ruhr.conj(), Ruhr @ C16)
nat_G2 = {"tau (D-Paar)": np.allclose(tau8 @ G2_8, G2_8 @ tau8), "U_theta": np.allclose(Uth @ G2_8, G2_8 @ Uth),
          "Lorentz": all(np.allclose(X @ G2_8, G2_8 @ X) for X in gerade[:6]), "J": J_G2}
ZAHLEN["vorzeichen"] = {"Basis Tr(R gamma_0124 P_- S^+)": str(basis_w), **{k: str(v) for k, v in var.items()},
                        "gamma_5 = lam i Gamma_2 R auf S^+": str(np.round(g5_aus[0], 9)), "native Symmetrien vertauschen mit Gamma_2": nat_G2,
                        "J R J^-1 = R": J_R}
check(f"[E] Vorzeichen des eps-Terms (Basis {np.round(basis_w.real, 3)}): er kippt unter jeder der vier Umkehrungen "
      f"{ {k: round(v.real, 3) for k, v in var.items()} } ({alle_kippen}). Drei sind festgelegt: die Uhrrichtung durch die "
      f"Marken (J R J^-1 = R: {J_R}), S^+ durch die linkshändige 16 (gelöster Block), und auf S^+ ist die 4D-Orientierung "
      f"keine eigene Größe: gamma_5 = {np.round(g5_aus[0], 6)} i Gamma_2 R (Residuum {g5_aus[1]:.0e}). Offen bleibt die "
      f"Zuordnung des Nullmodensektors: alle nativen Symmetrien vertauschen mit Gamma_2 ({nat_G2}), keine tauscht P_+ und "
      f"P_-; das ist das Vorzeichen der Wirkung (z103: Einzelsektor)",
      alle_kippen and J_R and g5_aus[1] < 1e-12 and abs(g5_aus[0] - 1) < 1e-12 and all(nat_G2.values()) and ok_C10)
befund("E/C: Die Spur mit Uhr ist nativ verfügbar, aber nicht ausgezeichnet. Auf S^+ ist R = -i Gamma_2 gamma_5 = Z_SU^-1, "
       "die relative Spur mit R also der gamma_5-verdrehte McKean-Singer-Index (Lefschetz-Zahl). Die zulässigen Gewichte "
       "bilden span{Gamma_2, R}; die kanonische Indexgradierung Gamma_2 allein gibt den topologischen Sektor (Pontrjagin + "
       "Holst), R gibt MacDowell-Mansouri (Euler + Einstein + Lambda = -3 sigma/l^2), jede Mischung mit b != 0 Einstein-"
       "Dynamik. Die Wahl W = R ist eine Prämisse P_R (dieselbe Information wie P_vol, jetzt als Gewicht statt als "
       "Plebanski-Sektor: das Modul-U(1) wirkt auf die Gewichte trivial, also weder Curie-Verbot noch Auswahl), dazu muss R "
       "allein auf dem inneren Faktor wirken (sonst s-abhängig), und das Vorzeichen braucht das Zuordnungsbit. Euler bleibt "
       "(Gewicht l^2/(64 pi G) relativ), Lambda ist nativ negativ", True)
laufzeit()

# =================================================================== Teil 3
say("\nTeil 3: Normierung von G")

# --- 3.1 Globale gegen lokale Kopplung
ss = (0.01, 0.3, 1.0, 4.0, 20.0, 60.0)
Ag = A_rechteck(81, QQ)
_, _, vh_ = np.linalg.svd(Ag)
zm = vh_[-1]*np.sign(vh_[-1][0])
P0s, P0t = np.zeros((81, 81)), np.zeros((80, 80))
P0s[0, 0] = P0t[0, 0] = 1
Pzs = np.outer(zm, zm)
Pzt = np.outer(zm[:80], zm[:80])
w_glob = [rel_spur(Ag, s=s_).real for s_ in ss]
w_rand = [rel_spur(Ag, P0s, P0t, s_).real for s_ in ss]
w_zm = [rel_spur(Ag, Pzs, Pzt, s_).real for s_ in ss]
Zrand = zm[0]**2
stoer = {}
for lab, kw in (("Symbol 1 - q z + 0,1 z^2", dict(p=[1.0, -QQ, 0.1])), ("Symbol 1 - q z - 0,2 z^2", dict(p=[1.0, -QQ, -0.2])),
                ("Randblock 6 x 6 (0,3)", dict(block=0.3)), ("Randblock 6 x 6 (0,6)", dict(block=0.6))):
    A = A_rechteck(81, QQ, kw.get("p"))
    if "block" in kw:
        K = rng.normal(size=(6, 6))
        A = A.copy()
        A[:6, :6] += kw["block"]*K/np.linalg.norm(K, 2)
    k_, c_ = index_von(A, 1e-10)
    _, _, vh2 = np.linalg.svd(A)
    z2 = vh2[-1]
    stoer[lab] = {"Index": k_ - c_, "globale Spur (s = 0,3/4)": [round(rel_spur(A, s=s_).real, 9) for s_ in (0.3, 4.0)],
                  "Randgewicht Z": float(abs(z2[0])**2), "Gewicht am Abschneiderand": float(abs(z2[-1])**2)}
glob_ok = all(abs(x - 1) < 1e-9 for x in w_glob) and all(
    v["Index"] == 1 and all(abs(x - 1) < 1e-8 for x in v["globale Spur (s = 0,3/4)"]) and v["Gewicht am Abschneiderand"] < 1e-10
    for v in stoer.values())
Z_streu = [v["Randgewicht Z"] for v in stoer.values()]
ZAHLEN["kopplung"] = {"s": ss, "global": w_glob, "erster Randplatz": w_rand, "ganze Nullmode": w_zm,
                      "|zm_0|^2": Zrand, "1 - q^2": 1 - QQ**2, "Störungen": stoer}
check(f"[E] Kopplungsgewichte an den Nahtzustand (N = 81): global W = 1 gibt die relative Spur {np.round(w_glob, 9).tolist()} "
      f"(Index, s-unabhängig); erster Randplatz W = |0><0|: {np.round(w_rand, 8).tolist()} (s-abhängig, von 0 bis |zm_0|^2 "
      f"= {Zrand:.11f} = 1 - q^2 = {1 - QQ**2:.11f}, Polgewicht Z aus z47/z48); ganze Nullmode W = |zm><zm|: "
      f"{np.round(w_zm, 8).tolist()} (-> 1). Unter Symbol- und Randstörungen bleiben Index und globale Spur 1 ({glob_ok}), "
      f"das Randgewicht wandert: {[round(z, 6) for z in Z_streu]}. Nur die globale Kopplung ist indexgeschützt; das "
      f"Ward-Gewicht 1 - q^2 ist Grenzwert einer lokalen, nicht geschützten Spur. Zwei native Prinzipien wählen also "
      f"verschieden: Ward-Universalität (z48) den Randzustand mit 1 - q^2, Indexschutz die globale Kopplung mit 1 "
      f"(Unterschied e^-4pi)",
      glob_ok and abs(Zrand - (1 - QQ**2)) < 1e-12 and abs(w_rand[0]) < 0.02 and abs(w_rand[-1] - (1 - QQ**2)) < 1e-12
      and abs(w_zm[-1] - 1) < 1e-9 and max(abs(z - (1 - QQ**2)) for z in Z_streu) > 1e-3)

# --- 3.2 Ganzzahligkeit und J-Summanden
ganz = all(abs(v - round(v.real)) < 1e-9 for v in list(TR.values()) + list(TG.values()))
F16 = kr(np.eye(8), (np.eye(32) + G11)/2)
Pi_a = Sp16 @ Pminus @ F16
Pi_b = Sp16 @ Pminus @ (np.eye(256) - F16)
rg = (int(np.linalg.matrix_rank(Pi_a)), int(np.linalg.matrix_rank(Pi_b)))
J_tauscht = C16 is not None and np.allclose(C16 @ Pi_a.conj() @ np.linalg.inv(C16), Pi_b)
X0124 = GABCD[(0, 1, 2, 4)]
t_a = complex(np.trace(Ruhr @ X0124 @ Pi_a))
t_b = complex(np.trace(Ruhr @ X0124 @ Pi_b))
ZAHLEN["summanden"] = {"Ränge": rg, "J tauscht": J_tauscht, "Tr(R gamma_0124) je Summand": [str(t_a), str(t_b)], "ganzzahlig": ganz}
check(f"[E] Die relative Spur liefert nur ganze Zahlen (alle 36 Werte beider Gewichte ganzzahlig: {ganz}); ein "
      f"transzendentes Gewicht wie c_3 = 1/(8 pi) kann sie nicht erzeugen. S^+ cap P_- zerfällt in zwei Lorentz x "
      f"Spin(10)-Summanden (Ränge {rg}), die die Realstruktur J (z116, J^2 = 1) aufeinander abbildet ({J_tauscht}); jeder "
      f"trägt die Hälfte des eps-Terms ({np.round(t_a.real, 3)} + {np.round(t_b.real, 3)} = {np.round((t_a + t_b).real, 3)}). "
      f"Der Faktor 2 in '2 c_3' zählt also die beiden J-Summanden: Dirac-Zählung; die unabhängigen Majorana-Weyl-Daten "
      f"(Pfaffsche, ein Summand) geben die Hälfte",
      ganz and rg == (32, 32) and J_tauscht and abs(t_a - t_b) < 1e-9 and abs(abs(t_a + t_b) - 64) < 1e-9)

# --- 3.3 Vergleich mit GRAV.NONLINEAR.01
c3 = 1/(8*sp.pi)
v_s = sp.Symbol("v", positive=True)
Gs = sp.Symbol("G", positive=True)
G_dirac = sp.solve(sp.Eq(2*c3/4*v_s**2, 1/(16*sp.pi*Gs)), Gs)[0]
G_mw = sp.solve(sp.Eq(1*c3/4*v_s**2, 1/(16*sp.pi*Gs)), Gs)[0]
G_ledger = sp.solve(sp.Eq(8*sp.pi*Gs*v_s**2, 1/c3), Gs)[0]          # 8 pi G = 1/c_3 in v_geo-Einheiten
q2 = mp.exp(-4*mp.pi)
gew = {"ganze Nullmode (Index)": mp.mpf(1), "erster Randplatz 1 - e^-4pi": 1 - q2, "quartischer Überlapp (1 - e^-4pi)^2": (1 - q2)**2}
G_rel = {k: mp.nstr(1/w - 1, 6) for k, w in gew.items()}
ZAHLEN["G"] = {"Dirac-Zählung (2 c_3)": str(G_dirac), "Majorana-Weyl (c_3)": str(G_mw), "Ledger 8 pi G = 1/c_3 (v = 1)": str(G_ledger),
               "relative Verschiebung von G je Kopplung": G_rel}
check(f"[E/C] Normierung (sympy): (2 c_3/4) v^2 = 1/(16 pi G) gibt G = {G_dirac} (Dirac-Zählung), mit c_3 je unabhängigem "
      f"Summanden (Majorana-Weyl) G = {G_mw}; der Ledgerwert (GRAV.NONLINEAR.01, Jacobson-Gleichgewicht 8 pi G = 1/c_3 "
      f"in v_geo-Einheiten, nur gelesen) ist G = {G_ledger}. Die Kopplungsgewichte verschieben G um "
      f"{G_rel} (relativ), weit unter jeder Messbarkeit, aber verschieden: Index, Randplatz und quartischer Überlapp "
      f"normieren verschieden bei gleichem Index",
      G_dirac == v_s**-2 and G_mw == 2*v_s**-2 and G_ledger == v_s**-2 and abs(float(1/gew["erster Randplatz 1 - e^-4pi"] - 1)
                                                                        - float(q2/(1 - q2))) < 1e-15)
befund("E/C/O: G ist nicht nativ normiert. Der Index legt nur 'genau ein eps-Term' mit ganzzahligem Gewicht fest (64 je "
       "Index, 32 je J-Summand); global (Gewicht 1) ist die einzige indexgeschützte Kopplung, das Ward-Randgewicht "
       "1 - q^2 (z47/z48) ist Grenzwert einer lokalen Spur und wandert unter Störungen. c_3 = 1/(8 pi) kommt nicht aus der "
       "Spur, und der Faktor 2 ist die Dirac-Zählung der beiden J-Summanden; mit Majorana-Weyl-Zählung wäre G = 2 v_geo^-2. "
       "G = v_geo^-2 (GRAV.NONLINEAR.01) wird nur mit dem gesetzten Gewicht 2 c_3 reproduziert: Eingabe, nicht Herleitung",
       True)
laufzeit()

# =================================================================== Teil 4
say("\nTeil 4: Bilanz und Urteil")
bilanz = [
    "[E] lichtartiges Paar nativ: Nullpaar der Dilatation (Ebenen von I_pi, I_omega; P_+- auf T = deren Projektoren, "
    "Gamma_2 = -2 D auf T); daneben das Kappenpaar (tau tauscht), das Uhrpaar E ist komplex (euklidischer Querraum)",
    "[E] Gamma_7 = c Gamma_2 gamma_5; auf S^+ (256): R = -i Gamma_2 gamma_5 = i (h7 x 1) = Z_SU^-1: 14D-Weyl leistet "
    "die Verklebung R = Z_SU^-1 von P_A3 (z97)",
    "[E, negativ] D_perp ohne Uhrladung (Index im Viertelschritt, Ladung -i; Spinorseite uhrneutral); nur D^2 und die "
    "Spuren invariant; Kompensator = imaginäres Dilatationselement (Wick); uhrneutrale Schritte: Transfer Index -+4, "
    "Möbius L_1 Index 1 je D_h, energiekovariant nur am komplexen Uhrpaar",
    "[E] zulässige Gewichte span{Gamma_2, R}: Gamma_2 (Indexgradierung) topologisch, R Einstein (64 eps); Holst-Mischung; "
    "Dualität Tr(R B Y) = Tr(Gamma_2 *B Y): Gewichtswahl = Paarwahl von P_vol",
    "[E, negativ] native Uhr auf beiden Faktoren: f(s) = 1 - (1+i) delta(s), s-abhängig und komplex: P_R (R nur innen)",
    "[E] MM mit nativen Daten: Euler + Einstein + Lambda = -3 sigma/l^2 = -12/R_Vak^2 (AdS); Gamma_2: Pontrjagin, Holst = 0",
    "[E/O] Vorzeichen: vier Bits, drei fest (Marken, 16, gamma_5 = i Gamma_2 R), Zuordnung des Nullmodensektors offen",
    "[E] Kopplung: nur global indexgeschützt (1); 1 - q^2 lokal, s-abhängig, störungsabhängig",
    "[E/C] G: Spur ganzzahlig, c_3 importiert; 2 = Dirac-Zählung der J-Summanden; G = v^-2 nur mit Gewicht 2 c_3 (sonst 2 v^-2)",
]
for z in bilanz:
    say("   " + z)
urteil = ("U2 ist beantwortet, als bedingtes Ergebnis mit präzisierten Prämissen; nichts ist befördert. (1) Nativität von "
          "D_perp: Das lichtartige Paar der siebzehnten Analyse ist nativ das Nullpaar der Dilatation. In der aus "
          "Twistorbivektoren gebauten Klein-Clifford-Algebra auf T + T* sind I_pi und I_omega reelle Nullvektoren mit "
          "Q(I_pi, I_omega) = 2, gamma_+- = Gamma(I_pi)/2, Gamma(I_omega)/2, P_+- auf T die Projektoren auf die beiden "
          "Ebenen (tau = +1 bzw. -1) und Gamma_2 = -2 D auf T. Das Kappenpaar (Q^-, Spitzen) ist ebenfalls reell (tau vertauscht "
          "es), das uhrgraduierte Paar (E, Original/Partner) ist komplex mit euklidischem Querraum. Der Nahtindex +1 sitzt "
          "aber im Viertelschritt S^+, der die Uhrladung -i trägt (A_q ist Z_4-, nicht L_0-homogen), während jede native "
          "Spinorstruktur uhrneutral ist: D_perp = gamma_- A + gamma_+ A^+ hat keine Uhrladung, nur seine Uhrbahn, D^2 und "
          "die Spuren sind invariant; kompensieren könnte allein das imaginäre Dilatationselement e^(i pi Gamma_2/4). Uhrneutrale native Schritte "
          "haben Index -+4 (Transfer J_-1) oder sind Möbiussenker (Index 1 je D_h), energiekovariant nur am komplexen "
          "Uhrpaar (euklidisch). D^2 = P_- A^+A + "
          "P_+ AA^+ ist aus nativen Stücken baubar bis auf die Zuordnung des Nullmodensektors. (2) Wirkung: Auf S^+ ist "
          "die Uhr R = i Gamma_11 = -i Gamma_2 gamma_5 = Z_SU^-1 (die 14D-Weyl-Bedingung leistet dieselbe Verklebung wie "
          "P_A3, z97); die Spur "
          "mit Uhr ist der gamma_5-verdrehte McKean-Singer-Index. Die zulässigen Gewichte sind genau span{Gamma_2, R}: "
          "Die kanonische Indexgradierung Gamma_2 gibt den topologischen Sektor (Pontrjagin, Holst = 0), R den "
          "MacDowell-Mansouri-Sektor (64 eps je Index). Die Wahl R ist nicht erzwungen, sie ist die Prämisse P_R; wegen "
          "Tr(R B Y) = Tr(Gamma_2 (*B) Y) ist sie genau die Paarwahl von P_vol (z103) in Spurform, und das Modul-U(1) "
          "wirkt auf die Gewichte trivial (weder Curie-Verbot noch Auswahl). "
          "Wirkt die native Uhr wie R = e^(2 pi i L_0) auch auf den Marken, wird die Spur s-abhängig und komplex "
          "(f = 1 - (1+i) delta(s), |f - 1| bis 1,8e-6 bei q = e^-2pi): P_R verlangt die Uhr allein auf dem inneren Faktor. Die "
          "gemischte Antwort ist 2 eps R e e ohne Euler- und Volumenterm, beim nativen lambda = sigma/l^2 trägt die Wirkung "
          "aber das Euler-Gewicht und Lambda = -3 sigma/l^2 = -12/R_Vak^2 (AdS), in Übereinstimmung mit z98 und z104 "
          "(positive Zeit => Lambda <= 0); dS braucht weiter P_cap. Das Vorzeichen hängt an vier Bits, von denen drei "
          "festliegen (Marken, linkshändige 16, gamma_5 = i Gamma_2 R auf S^+); offen bleibt die Zuordnung des "
          "Nullmodensektors, keine native Symmetrie tauscht P_+ und P_-. (3) G: Die relative Spur ist ganzzahlig (64 je "
          "Index, 32 je J-Summand); global ist die einzige indexgeschützte Kopplung, das Ward-Randgewicht 1 - q^2 ist der "
          "Grenzwert einer lokalen, s- und störungsabhängigen Spur. c_3 = 1/(8 pi) kommt nicht aus der Spur, und der Faktor "
          "2 ist die Dirac-Zählung der beiden J-Summanden; mit Majorana-Weyl-Zählung wäre G = 2 v_geo^-2. G = v_geo^-2 "
          "(GRAV.NONLINEAR.01) wird nur mit dem gesetzten Gewicht 2 c_3 reproduziert. Urteil für K13 und die "
          "Gravitationsachse: bedingt, nicht nativ. Minimalprämissen: P_R (Uhr als Gewicht, nur innen; ersetzt P_vol in "
          "MM-Form), das Zuordnungsbit des Nullmodensektors (Vorzeichen) und das Gewicht 2 c_3 (G). K13 bleibt L.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)}); Laufzeit "
    f"{time.time() - T0:.0f} s")
(HERE / "z118_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil, "bilanz": bilanz}, indent=1,
    ensure_ascii=False, default=str) + "\n")
