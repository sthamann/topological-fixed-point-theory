"""Z112: Gravitation über den nativen Infinity-Twistor: Hodges-MHV, Doppelkopie und die Normierung von G (Tiefenwelle T3).

Leitlinie 'Rekonstruktion statt Auswahl'. Firewalls: GRAV.NONCIRCULAR.01 (Impulse sind die nativen P_mu = (1/2) c^+ sigma_mu c
mit der parabolischen Minkowskizeit, nicht der Modulfluss), GRAV.PARENT.GAP.NO_GO.01 (gerechnet werden nur Baumamplituden auf
masselosen Impulseigenzuständen; kein Elternoperator, kein Pol wird behauptet), DIMENSION.UPLIFT.FIREWALL.01 (alles ist
Spinor-, Twistor- und Fockkinematik einer Naht und rationale Baumamplituden als Vergleichsziel; keine wechselwirkende
3+1-Schließung). Externe Formeln (Parke-Taylor, Hodges' MHV-Formel, KLT mit dem Impulskern von Bjerrum-Bohr, Damgaard,
Søndergaard und Vanhove, Skinners Twistorstring mit Infinity-Twistor) sind Vergleichsziele, keine Prämissen. Gelöste Blöcke
(Flavour, Massen, alpha, 16, E8) bleiben unberührt; nichts wird befördert. Bestand (nur gelesen): z78/z89 (c = d + e^+
kommutativ, P_mu, I_pi nativ bis auf die Modulphase, |I| aus dem KMS-Zustand, Lorentzbüschel lambda I_pi + mu I_omega, g frei),
z98 (Phase = Eichrichtung, Sigma_+ mit eps auf dem Originalindex, Sigma_- = J Sigma_+ J, P_mu = (1/2)(J c J) sigma c, so(3,2)
mit l = R_Vak/2, Omega = 1 + det x, Lambda-G-Ziel G/R_Vak^2 = e^(-2/alpha)/(128 pi^3)), z103 (P_vol), z100 (Maurer-Cartan),
z109 Teil 7 (erste Variationen fixieren G nicht), z99 (su(2)_1-Ströme, Niveau 1), z85 (Koordinaten (w, c) = (z_j - z_j+2,
z_j + z_j+2), Inzidenz w = i x c). Statuskürzel: E exakt, C bedingt, L Lesart, O offen, N numerisch.

Teil 1  Native Klammern: c-Algebra und Impulseigenzustände (Fock 8^4), welcher native Infinity-Twistor welche Klammer gibt
        (Twistoren auf Geraden L_x, ebene Wellen in omega bzw. c), AdS-Büschel und Inzidenz, Modulphase als globale
        Kleingruppenphase, KMS-Normierung s_ij = <ij>[ij] = 2 p_i.p_j, Impulserhaltung, Dreipunktkinematik.
Teil 2  Gravitation MHV: Hodges' reduzierte Determinante mit nativen Klammern für n = 3..6 (Unabhängigkeit von Streichungen
        und Referenzspinoren, Negativkontrolle ohne Impulserhaltung), Vergleich mit KLT (Impulskern) aus Parke-Taylor,
        Helizitäts- und Modulladungen, die native Frage Original x Partner (J-Paarung je Bein gegen blattdiagonale Paarung),
        linearisierte Krümmung der Produkte (Weyl gegen spurfreien Ricci, Sigma-Spannen aus z98).
Teil 3  Normierung von G: Dimension der Doppelkopie-Kopplung L = kappa/(2 g^2), native Länge (Dilatation, Vakuum,
        |I|-Normierung), Prämissenabtastung gegen das Lambda-G-Ziel.
Teil 4  Bilanz und Urteil.
"""
import itertools
import json
import math
import time
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
rng = np.random.default_rng(112)
mp.mp.dps = 30
I2 = np.eye(2)
SIGMA = [np.eye(2, dtype=complex), np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]),
         np.diag([1., -1.]).astype(complex)]
EPS2 = np.array([[0., 1.], [-1., 0.]], dtype=complex)                  # eps^(12) = 1
ETA = np.diag([1., -1., -1., -1.])
KAP2_NATIV = {m: float(mp.tanh(mp.pi*m)) for m in (1, 2, 3)}           # kappa^2 = tanh(pi m) (z85, z89)
ALPHA_INV = mp.mpf("137.035999177")                                     # CODATA 2022, nur zur Größenordnung
LC = np.zeros((4,)*4)
for _p in itertools.permutations(range(4)):
    LC[_p] = round(np.linalg.det(np.eye(4)[list(_p)]))


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


def r(x, n=6):
    return float(round(float(x), n))


# ------------------------------------------------------------------ Fockraum einer Naht (wie z85, z89, z98)
def fock4(n):
    a = sps.diags(np.sqrt(np.arange(1, n)), 1, format="csr")
    Id = sps.identity(n, format="csr")
    ops = []
    for k in range(4):
        mats = [Id]*4
        mats[k] = a
        op = mats[0]
        for Mx in mats[1:]:
            op = sps.kron(op, Mx, format="csr")
        ops.append(op)
    return ops


def dag(A):
    return sps.csr_matrix(A).T.conj().tocsr()


def komm(A, B):
    return (A @ B - B @ A).tocsr()


def tausch(nv):
    """Vertauschung d <-> e auf dem abgeschnittenen Fockraum (reell); J = Tausch mal komplexe Konjugation (z89, z98)."""
    occ = np.array(list(itertools.product(range(nv), repeat=4)))
    dim = nv**4
    idx = [np.ravel_multi_index((o[2], o[3], o[0], o[1]), (nv,)*4) for o in occ]
    return sps.csr_matrix((np.ones(dim), (idx, np.arange(dim))), shape=(dim, dim)), occ


def J_ad(X, SWm):
    return (SWm @ sps.csr_matrix(X).conj() @ SWm).tocsr()


NV = 8
DIMV = NV**4
d1, d2, e1, e2 = fock4(NV)
DD, EE = [d1, d2], [e1, e2]
SW, OCC = tausch(NV)
DVAL = OCC[:, 0] + OCC[:, 1] - OCC[:, 2] - OCC[:, 3]
LOW = np.where(OCC.sum(1) <= NV - 3)[0]
OMV = np.zeros(DIMV, complex)
OMV[0] = 1


def klein(A, tol=1e-10):
    A = sps.csc_matrix(A)[:, LOW]
    return (abs(A).max() if A.nnz else 0.0) < tol


ZT = [d1, d2, dag(e1), dag(e2)]                                   # Z^alpha = (d_1, d_2, e_1^+, e_2^+)
CC = [(DD[i] + dag(EE[i])).tocsr() for i in range(2)]             # c = d + e^+ (pi-Ebene, tau = +1, Originalindex)
WO = [(DD[i] - dag(EE[i])).tocsr() for i in range(2)]             # omega'' = d - e^+ (tau = -1)
Pop = [0.5*sum(s[i, j]*(dag(CC[i]) @ CC[j]) for i in range(2) for j in range(2)) for s in SIGMA]
Eop = (0.5*sum((dag(X) @ X for X in DD + EE), sps.csr_matrix((DIMV, DIMV))) + sps.identity(DIMV)).tocsr()

# ------------------------------------------------------------------ Twistorraum (z-Koordinaten, wie z89/z98)
G4 = np.diag([1., 1., -1., -1.]).astype(complex)
S4 = np.block([[np.zeros((2, 2)), I2], [I2, np.zeros((2, 2))]]).astype(complex)
TMS = np.array([[1, 0, -1, 0], [0, 1, 0, -1], [1, 0, 1, 0], [0, 1, 0, 1]], float)   # (w, c) = TMS z
TMSi = np.linalg.inv(TMS)
uvec = [np.array([1, 0, -1, 0.]), np.array([0, 1, 0, -1.])]       # [c_j^+, f] = -(u_j . grad) f (z89)
vvec = [np.array([1, 0, 1, 0.]), np.array([0, 1, 0, 1.])]         # [omega''_j^+, f] = -(v_j . grad) f
Inp = (np.outer(uvec[0], uvec[1]) - np.outer(uvec[1], uvec[0])).astype(complex)      # I_pi (oberer Bivektor)
Ion = (np.outer(vvec[0], vvec[1]) - np.outer(vvec[1], vvec[0])).astype(complex)      # I_omega
NP = [-0.5*sum(s[i, j]*np.outer(uvec[i], vvec[j]) for i in range(2) for j in range(2)) for s in SIGMA]
NK = [-0.5*sum(s[i, j]*np.outer(vvec[i], uvec[j]) for i in range(2) for j in range(2)) for s in SIGMA]
NE = np.diag([-0.5, -0.5, 0.5, 0.5]).astype(complex)
MP = [1j*N for N in NP]
MKL = [ETA[m, m]*1j*NK[m] for m in range(4)]                       # K_a = eta_ab K^b (z98)
MEt = 1j*NE


def unten(I):
    """Unterer Infinity-Twistor (*I)_ab = (1/2) eps_abcd I^cd (Volumenform eps_1234 = 1 in z-Koordinaten, vgl. P_vol z103)."""
    return 0.5*np.einsum("abcd,cd->ab", LC, I)


def ang(a, b):
    return a[0]*b[1] - a[1]*b[0]


def auf_gerade(x, lam):
    """Twistor auf der Geraden L_x: (w, c) = (i x lam, lam), Inzidenz w = i x c (z98), in z-Koordinaten."""
    return TMSi @ np.concatenate([1j*(x @ lam), lam])


def zpunkt(xm):
    return sum(xm[k]*SIGMA[k] for k in range(4))


def zk():
    return rng.normal(size=2) + 1j*rng.normal(size=2)


# ------------------------------------------------------------------ Kinematik, Klammern, Amplituden (Vergleichsziele)
def psd_wurzel(P):
    w, V = np.linalg.eigh((P + P.conj().T)/2)
    return V @ np.diag(np.sqrt(np.clip(w, 0, None))) @ V.conj().T


def kin_komplex(n):
    """Komplexe Kinematik: lambda und lambda~ unabhängig (analytische Fortsetzung), Impulserhaltung linear gelöst."""
    lam = rng.normal(size=(n, 2)) + 1j*rng.normal(size=(n, 2))
    lt = rng.normal(size=(n, 2)) + 1j*rng.normal(size=(n, 2))
    rest = sum(np.outer(lam[i], lt[i]) for i in range(2, n))
    T = -np.linalg.solve(np.array([lam[0], lam[1]]).T, rest)
    lt[0], lt[1] = T[0], T[1]
    return lam, lt


def kin_reell(n):
    """Native reelle Kinematik: lambda = c-Eigenwert, lambda~ = eta conj(lambda) (J-Bild, eta = -1 für einlaufende Beine)."""
    lam = rng.normal(size=(n, 2)) + 1j*rng.normal(size=(n, 2))
    lt = lam.conj().copy()
    P = sum(np.outer(lam[i], lam[i].conj()) for i in range(2, n))
    Uq = np.linalg.qr(rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2)))[0]
    QU = psd_wurzel(P) @ Uq
    lam[0], lam[1] = QU[:, 0], QU[:, 1]
    lt[0], lt[1] = -QU[:, 0].conj(), -QU[:, 1].conj()
    return lam, lt


def kin_drei(n=3):
    """MHV-Dreipunktkinematik: alle lambda~ parallel (komplex), sum lambda_i lambda~_i = 0."""
    lam = rng.normal(size=(3, 2)) + 1j*rng.normal(size=(3, 2))
    l0 = zk()
    t = np.linalg.solve(np.array([lam[0], lam[1]]).T, -lam[2])
    lt = np.array([t[0]*l0, t[1]*l0, l0])
    return lam, lt


def klammern(lam, lt):
    n = len(lam)
    A = np.array([[ang(lam[i], lam[j]) for j in range(n)] for i in range(n)])
    S = np.array([[ang(lt[i], lt[j]) for j in range(n)] for i in range(n)])
    return A, S


def impulsrest(lam, lt):
    return float(np.abs(sum(np.outer(lam[i], lt[i]) for i in range(len(lam)))).max())


def pt(order, A, neg=(0, 1)):
    """Parke-Taylor (Vergleichsziel): <ab>^4/(<12><23>...<n1>), a, b negative Helizität."""
    den = 1
    for k in range(len(order)):
        den = den*A[order[k], order[(k + 1) % len(order)]]
    return A[neg[0], neg[1]]**4/den


def ptq(order, S, pos=(0, 1)):
    """Anti-MHV (konjugierte Helizitäten): [ab]^4/([12][23]...[n1])."""
    den = 1
    for k in range(len(order)):
        den = den*S[order[k], order[(k + 1) % len(order)]]
    return S[pos[0], pos[1]]**4/den


def hodges(lam, A, S, zeilen, spalten, x, y, neg=(0, 1)):
    """Hodges (Vergleichsziel): M_n = <ab>^8 det'(Phi), det' = (-1)^(n+1) sigma sigma' |Phi|^(ijk)_(rst)/(c_ijk c_rst)."""
    n = len(lam)
    Phi = np.zeros((n, n), complex)
    for i in range(n):
        for j in range(n):
            if i != j:
                Phi[i, j] = S[i, j]/A[i, j]
        Phi[i, i] = -sum(S[i, j]*ang(lam[j], x)*ang(lam[j], y)/(A[i, j]*ang(lam[i], x)*ang(lam[i], y))
                         for j in range(n) if j != i)
    kr = [k for k in range(n) if k not in zeilen]
    kc = [k for k in range(n) if k not in spalten]
    minor = np.linalg.det(Phi[np.ix_(kr, kc)]) if kr else 1.0
    a, b, c = zeilen
    p_, q_, s_ = spalten
    cz = A[a, b]*A[b, c]*A[c, a]
    cs = A[p_, q_]*A[q_, s_]*A[s_, p_]
    vz = (-1)**(sum(zeilen) + sum(spalten))
    return A[neg[0], neg[1]]**8*(-1)**(n + 1)*vz*minor/(cz*cs)


def kern(I, J, s):
    """Impulskern S[I|J]_1 = prod_t (s_(i_t 1) + sum_(q>t) theta(i_t, i_q) s_(i_t i_q)), theta = 1 bei umgekehrter Ordnung in J."""
    out = 1
    for t, it in enumerate(I):
        term = s[it, 0]
        for q in range(t + 1, len(I)):
            if J.index(it) > J.index(I[q]):
                term = term + s[it, I[q]]
        out = out*term
    return out


def klt(A, S, links, rechts):
    """KLT (Vergleichsziel): M_n = (-1)^(n+1) sum A~(n-1, n, sigma~, 1) S[sigma~|sigma] A(1, sigma, n-1, n); dazu sum |Terme|."""
    n = len(A)
    s = A*S
    mitte = list(range(1, n - 2))
    tot, betrag = 0, 0.0
    for sg in itertools.permutations(mitte):
        for sgt in itertools.permutations(mitte):
            term = rechts([n - 2, n - 1, *sgt, 0])*kern(list(sgt), list(sg), s)*links([0, *sg, n - 2, n - 1])
            tot = tot + term
            betrag += abs(term)
    return (-1)**(n + 1)*tot, betrag


def drehe(lam, lt, beine, al):
    """Kleingruppe je Bein: lambda -> e^(i al) lambda, lambda~ -> e^(-i al) lambda~ (Impulse fest)."""
    l2, t2 = lam.copy(), lt.copy()
    for i in beine:
        l2[i] = np.exp(1j*al)*lam[i]
        t2[i] = np.exp(-1j*al)*lt[i]
    return l2, t2


STD = ((0, 1, 2), (0, 1, 2))
XR, YR = zk(), zk()

# =================================================================== Teil 1
say("Teil 1: Native Klammern")

# --- 1.1 c-Algebra, Impulseigenzustände, KMS-Normierung, J und Modul-U(1) (Fock 8^4)
ok_komm = all(klein(komm(CC[i], dag(CC[j]))) and klein(komm(CC[i], CC[j])) for i in range(2) for j in range(2))
ok_Pc = all(klein(komm(Pop[m], CC[j])) and klein(komm(Pop[m], dag(CC[j]))) for m in range(4) for j in range(2))
ok_J = all(klein(J_ad(CC[j], SW) - dag(CC[j])) for j in range(2))
ok_JP = all(klein(J_ad(Pop[m], SW) - Pop[m]) for m in range(4))
Dop = sps.diags(DVAL.astype(float)).tocsr()
ok_JD = klein(J_ad(Dop, SW) + Dop)                                  # J D J = -D: J spiegelt Modulladung und Helizität
th = 0.7
Uth = sps.diags(np.exp(1j*th*DVAL)).tocsr()
ok_U = all(klein(Uth @ CC[j] @ dag(Uth) - np.exp(-1j*th)*CC[j], 1e-12) for j in range(2)) and np.allclose(Uth @ OMV, OMV)
G_kms = np.array([[np.vdot(CC[i] @ OMV, CC[j] @ OMV) for j in range(2)] for i in range(2)])
ok_kms = np.allclose(G_kms, I2, atol=1e-12)
spek = []
ok_spek = True
for _ in range(6):
    lam = zk()
    p = np.array([0.5*np.vdot(lam, SIGMA[m] @ lam) for m in range(4)])
    Pm = sum(p[m]*SIGMA[m] for m in range(4))
    ok_spek &= (np.allclose(Pm, np.outer(lam, lam.conj())) and abs(p[0]**2 - np.sum(p[1:]**2)) < 1e-12
                and p[0].real > 0 and np.abs(p.imag).max() < 1e-14)
    spek.append(float(p[0].real))
ZAHLEN["c_algebra"] = {"Fock": f"{NV}^4", "<Omega|c_i^+ c_j|Omega>": np.round(G_kms.real, 12).tolist(),
                       "U_theta c U_theta^+ = e^(-i theta) c": True, "J c J = c^+, J P J = P": bool(ok_J and ok_JP),
                       "p_0 = |lambda|^2/2 (Stichproben)": [r(v) for v in spek]}
check(f"c-Algebra und Impulseigenzustände (Fock {NV}^4): c_i und c_j^+ vertauschen paarweise, P_mu = (1/2) c^+ sigma_mu c "
      f"vertauscht mit c und c^+, also sind die Impulseigenzustände die gemeinsamen Eigenzustände der kommutativen "
      f"c-Algebra mit Spektrum lambda in C^2: p_mu = (1/2) lambda^+ sigma_mu lambda, sum p_mu sigma_mu = lambda lambda^+, "
      f"p^2 = 0, p_0 = |lambda|^2/2 > 0 (sechs Stichproben). J c J = c^+ (der Partnerspinor ist das J-Bild, lambda~ = "
      f"conj(lambda)), J P_mu J = P_mu und J D J = -D (J spiegelt die Modulladung); U_theta = e^(i theta D) gibt U c U^+ = "
      f"e^(-i theta) c und hält das KMS-Vakuum fest; der KMS-Zustand normiert die c-Ebene, <c_i^+ c_j> = delta_ij",
      ok_komm and ok_Pc and ok_J and ok_JP and ok_JD and ok_U and ok_kms and ok_spek)

# --- 1.2 welcher native Infinity-Twistor welche Klammer gibt (Twistoren auf Geraden, ebene Wellen)
q_pi, q_om, q_ads, q_tr = [], [], [], []
for _ in range(5):
    x = zpunkt(rng.normal(size=4))
    lams = [zk() for _ in range(4)]
    Zs = [auf_gerade(x, l) for l in lams]
    dx = np.linalg.det(x)
    a_ = zpunkt(rng.normal(size=4))
    Zt = [auf_gerade(x + a_, l) for l in lams]
    for i, j in itertools.combinations(range(4), 2):
        aij = ang(lams[i], lams[j])
        bp = Zs[i] @ unten(Inp) @ Zs[j]
        q_pi.append(bp/aij)
        q_om.append((Zs[i] @ unten(Ion) @ Zs[j])/(aij*dx))
        q_ads.append((Zs[i] @ unten(Inp - Ion) @ Zs[j])/(aij*(1 + dx)))
        q_tr.append(((Zt[i] @ unten(Inp) @ Zt[j]) - bp, (Zt[i] @ unten(Ion) @ Zt[j])/(Zs[i] @ unten(Ion) @ Zs[j])
                     - np.linalg.det(x + a_)/dx))
ok_lower = (np.allclose(q_pi, 1, atol=1e-12) and np.allclose(q_om, -1, atol=1e-12) and np.allclose(q_ads, 1, atol=1e-12)
            and max(abs(t[0]) for t in q_tr) < 1e-12 and max(abs(t[1]) for t in q_tr) < 1e-9)
xv = 1j*np.eye(2)
lv = [zk(), zk()]
om_vak = (auf_gerade(xv, lv[0]) @ unten(Inp - Ion) @ auf_gerade(xv, lv[1]))
ok_vak = abs(om_vak) < 1e-12 and abs(1 + np.linalg.det(xv)) < 1e-12
ob = []
for _ in range(5):
    lts = [zk(), zk()]
    g_ = [TMS[:2].T @ l for l in lts]                              # Gradient von exp(lambda~ . w)
    mu_ = [np.array([-(uvec[k] @ g) for k in range(2)]) for g in g_]   # c^+-Eigenwert: [c^+, f] = -(u . grad) f
    h_ = [TMS[2:].T @ l for l in lts]                              # Gradient von exp(nu . c)
    nu_ = [np.array([-(vvec[k] @ h) for k in range(2)]) for h in h_]   # omega''^+-Eigenwert
    ob.append((g_[0] @ Inp @ g_[1] - ang(mu_[0], mu_[1]), g_[0] @ Ion @ g_[1], h_[0] @ Ion @ h_[1] - ang(nu_[0], nu_[1]),
               h_[0] @ Inp @ h_[1], ang(mu_[0], mu_[1])/ang(lts[0], lts[1])))
ok_upper = all(abs(t[0]) < 1e-12 and abs(t[1]) < 1e-12 and abs(t[2]) < 1e-12 and abs(t[3]) < 1e-12
               and abs(t[4] - 4) < 1e-12 for t in ob)
ZAHLEN["klammern"] = {"*I_pi(Z_i, Z_j)/<ij>": r(np.mean(np.real(q_pi)), 12),
                      "*I_omega(Z_i, Z_j)/(<ij> det x)": r(np.mean(np.real(q_om)), 12),
                      "*(I_pi - I_omega)(Z_i, Z_j)/(<ij> (1 + det x))": r(np.mean(np.real(q_ads)), 12),
                      "Vakuumpunkt x = i 1: *(I_pi - I_omega)": float(abs(om_vak)),
                      "B_(I_pi)(e^(l~ w), e^(l~' w)) = eps(mu~, mu~'), mu~ = c^+-Eigenwert = -2 l~": True,
                      "B_(I_omega) auf c-Wellen = eps der omega''^+-Eigenwerte (K-Spinoren)": True}
check(f"Native Infinity-Twistoren und ihre Klammern (Twistoren auf Geraden L_x mit Inzidenz w = i x c, 5 Punkte x, je 6 "
      f"Paare): Der untere Twistor *I_pi gibt genau die Winkelklammer, *I_pi(Z_i, Z_j) = <ij> = eps(c-Eigenwerte) mit "
      f"Koeffizient 1 und translationsinvariant; *I_omega gibt -det(x) <ij> (ortsabhängig, Inversionsbild); das AdS-Büschel "
      f"*(I_pi - I_omega) gibt Omega(x) <ij> mit Omega = 1 + det x, genau z98s konformer Faktor, der am Vakuumpunkt x = i 1 "
      f"verschwindet ({abs(om_vak):.1e}). Der obere Twistor I_pi (Poissonklammer) gibt auf ebenen Wellen exp(l~ . w) "
      f"genau die eckige Klammer der c^+-Eigenwerte, B_(I_pi)(f_i, f_j) = eps(mu~_i, mu~_j) f_i f_j mit Koeffizient 1 "
      f"(mu~ = -2 l~, also 4 [l~ l~']), und I_omega gibt dort 0; auf c-Wellen exp(nu . c) ist es umgekehrt: I_omega gibt "
      f"die Klammer der omega''^+-Eigenwerte (Eigenzustände von K_mu), I_pi gibt 0", ok_lower and ok_vak and ok_upper)
befund("E: Nativ ist <ij> = eps(c_i, c_j), die eps-Form auf dem Originalindex (z98: I_pi = Sigma_+-Kontraktion), und "
       "[ij] = eps(c_i^+, c_j^+) = eps(J c_i J, J c_j J), die Poissonklammer B_(I_pi) = m_Z o (eps c^+ (x) c^+) aus z89, "
       "also die eps-Form des Tomita-Partners. Ein und derselbe native Bivektor I_pi liefert beide Klammern: unten (Hodge) "
       "die Winkel-, oben (Poisson) die eckige Klammer, beide mit Koeffizient 1. I_omega gehört zu den K-Eigenzuständen "
       "(Inversion), und das native AdS-Büschel I_pi - I_omega reskaliert die Winkelklammer auf L_x mit Omega(x) = 1 + det x", True)

# --- 1.3 Modulphase = globale Kleingruppenphase; KMS-Normierung s_ij = <ij>[ij] = 2 p_i.p_j; J-Paarung
lam, lt = kin_reell(5)
A, S = klammern(lam, lt)
eta = np.array([-1, -1, 1, 1, 1])
P4 = [np.array([0.5*eta[i]*np.vdot(lam[i], SIGMA[m] @ lam[i]) for m in range(4)]).real for i in range(5)]
ok_s = all(abs(A[i, j]*S[i, j] - 2*(P4[i] @ ETA @ P4[j])) < 1e-10*max(1, abs(A[i, j]*S[i, j]))
           for i, j in itertools.combinations(range(5), 2))
ok_Jp = all(abs(S[i, j] - eta[i]*eta[j]*np.conj(A[i, j])) < 1e-12 for i, j in itertools.combinations(range(5), 2))
ok_Pist = np.abs(sum(P4)).max() < 1e-12
l2, t2 = drehe(lam, lt, range(5), th)
A2, S2 = klammern(l2, t2)
ok_mod = (np.allclose(A2, np.exp(2j*th)*A) and np.allclose(S2, np.exp(-2j*th)*S) and np.allclose(A2*S2, A*S)
          and all(np.allclose(np.outer(l2[i], t2[i]), np.outer(lam[i], lt[i])) for i in range(5)))
ZAHLEN["normierung_kms"] = {"s_ij = <ij>[ij] = 2 p_i.p_j (n = 5, reell)": bool(ok_s), "[ij] = eta_i eta_j conj<ij>": bool(ok_Jp),
                            "Modulphase: <ij> Ladung, [ij] Ladung": [2, -2]}
check(f"KMS-Normierung und Modulphase: Mit den nativen p_mu = (1/2) eta lambda^+ sigma_mu lambda (eta = -1 für die zwei "
      f"einlaufenden Beine, Summe der p_mu = 0 auf {np.abs(sum(P4)).max():.1e}) gilt <ij>[ij] = 2 p_i.p_j (Minkowski) für alle "
      f"10 Paare ohne freien Faktor, und [ij] = eta_i eta_j conj<ij> (J-Paarung). Die Modulphase U_theta wirkt auf allen "
      f"Beinen zugleich als lambda -> e^(i theta) lambda, lambda~ -> e^(-i theta) lambda~: <ij> trägt die Ladung 2, [ij] "
      f"die Ladung -2, alle p_i und s_ij sind invariant. Die Modulphase von I ist also die gemeinsame Kleingruppenphase "
      f"aller Beine", ok_s and ok_Jp and ok_Pist and ok_mod)

# --- 1.4 Impulserhaltung, Inzidenz, Dreipunktkinematik
res_ip, res_sch, res_s = [], [], []
for n in (4, 5, 6):
    for gen_ in (kin_komplex, kin_reell):
        lam, lt = gen_(n)
        A, S = klammern(lam, lt)
        res_ip.append(impulsrest(lam, lt))
        M = A @ S                                                  # (A S)_ik = sum_j <ij>[jk]
        res_sch.append(float(np.abs(M).max()/np.abs(A).max()/np.abs(S).max()))
        sij = A*S
        res_s.append(float(np.abs(sij.sum(1)).max()/np.abs(sij).max()))
lam3 = zk()
P3 = np.outer(lam3, lam3.conj())
Q3 = psd_wurzel(P3) @ np.linalg.qr(rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2)))[0]
l3 = np.array([Q3[:, 0], Q3[:, 1], lam3])
A3, _ = klammern(l3, l3.conj())
drei_reell = float(np.abs(A3).max()/np.linalg.norm(lam3)**2)
lamk, ltk = kin_drei()
Ak, Sk = klammern(lamk, ltk)
drei_k = (float(np.abs(Sk).max()), float(np.abs(Ak[0, 1]*Ak[1, 2]*Ak[2, 0])))
ZAHLEN["impulserhaltung"] = {"Rest sum lambda lambda~": max(res_ip), "sum_j <ij>[jk] (rel.)": max(res_sch),
                             "sum_j s_ij (rel.)": max(res_s), "n = 3 reell: max |<ij>|/|lambda|^2": drei_reell,
                             "n = 3 komplex: max |[ij]|, |<12><23><31>|": drei_k}
check(f"Impulserhaltung und Dreipunktkinematik: Auf komplexer und nativ reeller Kinematik (n = 4, 5, 6) ist sum lambda_i "
      f"lambda~_i = 0 (Rest {max(res_ip):.1e}), sum_j <ij>[jk] = 0 für alle i, k (rel. {max(res_sch):.1e}) und sum_j s_ij = 0 "
      f"(rel. {max(res_s):.1e}). Für n = 3 erzwingt die J-Paarung (lambda~ = +-conj(lambda)) Kollinearität: alle <ij> "
      f"verschwinden (rel. {drei_reell:.1e}); die MHV-Dreipunktkinematik braucht lambda~ unabhängig von conj(lambda) "
      f"(alle [ij] = 0, <12><23><31> = {drei_k[1]:.3f} != 0), also eine komplexifizierte, nicht J-gepaarte Kinematik",
      max(res_ip) < 1e-12 and max(res_sch) < 1e-12 and max(res_s) < 1e-12 and drei_reell < 1e-12 and drei_k[0] < 1e-12
      and drei_k[1] > 1e-3)
befund("E: Die Bausteine der Gravitationsamplituden sind nativ und ohne freien Faktor normiert: <ij> aus dem Original, "
       "[ij] aus dem Partner (J-Bild), s_ij = <ij>[ij] = 2 p_i.p_j mit den nativen P_mu, Impulserhaltung als Rang-2-Bedingung "
       "sum lambda lambda~ = 0. Die Modulphase ist die globale Kleingruppenphase (keine Observable). Für n = 3 braucht es "
       "die komplexifizierte Kinematik, in der die J-Paarung lambda~ = +-conj(lambda) aufgehoben ist", True)
laufzeit()

# =================================================================== Teil 2
say("\nTeil 2: Gravitation MHV")

# --- 2.1 Hodges: Unabhängigkeit von Streichungen und Referenzspinoren, Negativkontrolle
streu = {}
anzahl = 0
for n in (4, 5, 6):
    smax = 0.0
    for gen_ in [kin_komplex]*4 + [kin_reell]*3:
        lam, lt = gen_(n)
        A, S = klammern(lam, lt)
        werte = []
        for _ in range(3):
            x, y = zk(), zk()
            for zl in itertools.combinations(range(n), 3):
                for sp_ in itertools.combinations(range(n), 3):
                    werte.append(hodges(lam, A, S, zl, sp_, x, y))
        werte = np.array(werte)
        anzahl += len(werte)
        smax = max(smax, float(np.abs(werte - werte[0]).max()/abs(werte[0])))
    streu[n] = smax
neg = {}
for n in (4, 5, 6):
    lam, lt = kin_komplex(n)
    lt[n - 1] = lt[n - 1] + 1e-3*zk()
    A, S = klammern(lam, lt)
    werte = np.array([hodges(lam, A, S, zl, sp_, x, y) for x, y in ((XR, YR), (zk(), zk()))
                      for zl in itertools.combinations(range(n), 3) for sp_ in itertools.combinations(range(n), 3)])
    neg[n] = float(np.abs(werte - werte[0]).max()/abs(werte[0]))
ZAHLEN["hodges_unabhaengigkeit"] = {"max. rel. Streuung (n: Wert)": streu, "Auswertungen": anzahl,
                                    "Negativkontrolle (Impulsrest 1e-3)": neg}
check(f"Hodges' reduzierte Determinante mit nativen Klammern: Für n = 4, 5, 6 auf je 4 komplexen und 3 nativ reellen "
      f"Kinematikpunkten, allen C(n,3)^2 Streichungen von Zeilen und Spalten und je 3 Referenzpaaren (x, y) ({anzahl} "
      f"Auswertungen) ist M_n = <12>^8 det'(Phi) eindeutig, max. rel. Streuung {max(streu.values()):.1e}; ohne "
      f"Impulserhaltung (Rest 1e-3 an einem Bein) hängt sie von der Wahl ab (Streuung "
      f"{ {k: float('%.1e' % v) for k, v in neg.items()} }): Die native Rang-2-Bedingung sum lambda lambda~ = 0 trägt die "
      f"Wohldefiniertheit (n = 3: Phi ist leer, det' = 1/(<12><23><31>)^2)",
      max(streu.values()) < 1e-9 and min(neg.values()) > 1e-6)

# --- 2.2 Hodges gegen KLT (Impulskern) aus Parke-Taylor
quot = {}
for n in (3, 4, 5, 6):
    qs = []
    gens = [kin_drei]*5 if n == 3 else [kin_komplex]*5 + [kin_reell]*3
    for gen_ in gens:
        lam, lt = gen_(n)
        A, S = klammern(lam, lt)
        H = hodges(lam, A, S, *STD, XR, YR)
        K, _ = klt(A, S, lambda o: pt(o, A), lambda o: pt(o, A))
        qs.append(H/K)
    quot[n] = float(np.abs(np.array(qs) - 1).max())
M3 = hodges(lamk, Ak, Sk, *STD, XR, YR)
M3_soll = Ak[0, 1]**6/(Ak[1, 2]**2*Ak[2, 0]**2)
ZAHLEN["hodges_klt"] = {"max |Hodges/KLT - 1| (n: Wert)": quot, "M_3 = <12>^6/(<23>^2 <31>^2)": float(abs(M3/M3_soll - 1)),
                        "KLT-Terme (n - 3)!^2": {n: math.factorial(n - 3)**2 for n in (3, 4, 5, 6)}}
check(f"Hodges = KLT: Mit s_ij = <ij>[ij] und dem Impulskern S[sigma~|sigma]_1 (Vergleichsziel) ist M_n^Hodges = "
      f"(-1)^(n+1) sum A(n-1, n, sigma~, 1) S A(1, sigma, n-1, n) aus Parke-Taylor auf beiden Seiten für n = 3, 4, 5, 6 "
      f"(5 komplexe und, ab n = 4, 3 nativ reelle Punkte; (n-3)!^2 = 1, 1, 4, 36 Terme), max |Hodges/KLT - 1| = "
      f"{max(quot.values()):.1e}; M_3 = <12>^6/(<23>^2 <31>^2) = A_3^2 (rel. {abs(M3/M3_soll - 1):.1e})",
      max(quot.values()) < 1e-8 and abs(M3/M3_soll - 1) < 1e-12)

# --- 2.3 Helizitäts- und Modulladungen
al = 0.3
lad = {}
ok_lad = True
for n in (4, 5, 6):
    lam, lt = kin_reell(n)
    A, S = klammern(lam, lt)
    H = hodges(lam, A, S, *STD, XR, YR)
    Ao = pt(list(range(n)), A)
    qH, qA = [], []
    for i in range(n):
        l2, t2 = drehe(lam, lt, [i], al)
        A2, S2 = klammern(l2, t2)
        qH.append(np.angle(hodges(l2, A2, S2, *STD, XR, YR)/H)/al)
        qA.append(np.angle(pt(list(range(n)), A2)/Ao)/al)
    l2, t2 = drehe(lam, lt, range(n), al)
    A2, S2 = klammern(l2, t2)
    QH = np.angle(hodges(l2, A2, S2, *STD, XR, YR)/H)/al if n != 4 else abs(hodges(l2, A2, S2, *STD, XR, YR)/H - 1)
    lad[n] = {"q_i Gravitation": [round(v, 9) for v in qH], "q_i Eich": [round(v, 9) for v in qA], "Q global": round(float(QH), 9)}
    soll_H = [4, 4] + [-4]*(n - 2)
    ok_lad &= np.allclose(qH, soll_H, atol=1e-9) and np.allclose(np.array(qH), 2*np.array(qA), atol=1e-9)
    ok_lad &= (abs(QH) < 1e-9) if n == 4 else abs(QH - (16 - 4*n)) < 1e-9
l2, t2 = drehe(lamk, ltk, range(3), al)
A2, S2 = klammern(l2, t2)
Q3 = np.angle(hodges(l2, A2, S2, *STD, XR, YR)/M3)/al
ok_lad &= abs(Q3 - 4) < 1e-9
lad[3] = {"Q global": round(float(Q3), 9)}
ZAHLEN["ladungen"] = lad
check(f"Helizitäts- und Modulladungen: Unter der Kleingruppe eines Beins (lambda -> e^(i a) lambda, lambda~ -> e^(-i a) "
      f"lambda~) trägt M_n die Ladung q = +4 an den Beinen negativer und -4 an denen positiver Helizität, genau das "
      f"Doppelte der Parke-Taylor-Ladung +-2 (n = 4, 5, 6, nativ reell). Die globale Modulladung ist Q(M_n) = 16 - 4n: "
      f"Q(M_3) = {Q3:.6f} (z98s Faktor e^(4 i theta) des kubischen Vertex), Q(M_4) = 0 (die Phase fällt heraus), "
      f"Q(M_5) = {lad[5]['Q global']:.6f}, Q(M_6) = {lad[6]['Q global']:.6f}", ok_lad)

# --- 2.4 Original x Partner je Bein (J-Paarung) gegen blattdiagonale Paarung
paar = {}
ok_paar = True
for n in (4, 5, 6):
    rows = []
    for _ in range(3):
        lam, lt = kin_reell(n)
        A, S = klammern(lam, lt)
        H = hodges(lam, A, S, *STD, XR, YR)
        Kmix, Bmix = klt(A, S, lambda o: pt(o, A), lambda o: np.conj(pt(o, A)))
        Kdia, _ = klt(A, S, lambda o: pt(o, A), lambda o: np.conj(ptq(o, S)))
        konjq = max(abs(np.conj(ptq(o, S)) - pt(o, A))/abs(pt(o, A)) for o in itertools.permutations(range(n)))
        qmix = []
        if n == 4:
            for i in range(n):
                l2, t2 = drehe(lam, lt, [i], al)
                A2, S2 = klammern(l2, t2)
                qmix.append(float(np.angle(klt(A2, S2, lambda o: pt(o, A2), lambda o: np.conj(pt(o, A2)))[0]/Kmix)/al))
        rows.append({"Hodges/blattdiagonal": complex(H/Kdia), "|Mix|/sum|Terme|": float(abs(Kmix)/Bmix),
                     "|Mix/Hodges|": float(abs(Kmix/H)), "arg(Mix/Hodges)": float(np.angle(Kmix/H)),
                     "arg Hodges": float(np.angle(H)),
                     "Im Mix/|Mix|": float(abs(Kmix.imag)/abs(Kmix)) if abs(Kmix) > 0 else 0.0,
                     "conj(Anti-MHV) = MHV (rel.)": float(konjq), "q_i Mix": [round(v, 9) for v in qmix]})
    dia = max(abs(rw["Hodges/blattdiagonal"] - 1) for rw in rows)
    konj = max(rw["conj(Anti-MHV) = MHV (rel.)"] for rw in rows)
    ok_paar &= dia < 1e-8 and konj < 1e-12
    if n == 4:
        ok_paar &= (all(rw["|Mix|/sum|Terme|"] > 1e-3 for rw in rows) and all(rw["Im Mix/|Mix|"] < 1e-10 for rw in rows)
                    and all(np.allclose(rw["q_i Mix"], 0, atol=1e-9) for rw in rows)
                    and all(abs(rw["|Mix/Hodges|"] - 1) < 1e-9 for rw in rows)
                    and max(abs(math.sin(rw["arg Hodges"])) for rw in rows) > 1e-2)
    else:
        ok_paar &= all(rw["|Mix|/sum|Terme|"] < 1e-10 for rw in rows)
    paar[n] = {"max |Hodges/blattdiagonal - 1|": dia, "|Mix|/sum|Terme|": [float("%.2e" % rw["|Mix|/sum|Terme|"]) for rw in rows],
               "q_i Mix (n = 4)": rows[0]["q_i Mix"] if n == 4 else None,
               "|Mix/Hodges|, arg (n = 4)": [[r(rw["|Mix/Hodges|"], 9), r(rw["arg(Mix/Hodges)"], 4)] for rw in rows]
               if n == 4 else None}
ZAHLEN["original_partner"] = paar
check(f"Original x Partner je Bein (nativ reelle Kinematik, je 3 Punkte): Die J-Paarung desselben Zustands, A_rechts = "
      f"conj(A_links) (J antiunitär, Partner mit gespiegelter Helizität), gibt im KLT-Produkt je Bein die Ladung 0: für "
      f"n = 4 (ein einziger KLT-Term) eine reelle, helizitätsneutrale Amplitude (q_i = {paar[4]['q_i Mix (n = 4)']}) vom "
      f"selben Betrag wie Hodges, aber mit anderer Phase und Helizität (|Mix/Hodges|, arg = "
      f"{paar[4]['|Mix/Hodges|, arg (n = 4)']}); für n = 5 und 6 ist sie exakt null (|Mix|/sum|Terme| = "
      f"{paar[5]['|Mix|/sum|Terme|']}, {paar[6]['|Mix|/sum|Terme|']}). Die blattdiagonale Paarung, A_rechts = "
      f"conj(Anti-MHV) = Partnerbild des helizitätsgespiegelten Zustands, ist auf J-gepaarter Kinematik exakt Parke-Taylor "
      f"(rel. {max(p_['max |Hodges/blattdiagonal - 1|'] for p_ in paar.values()):.1e} gegen Hodges): Nur sie gibt die "
      f"Gravitation", ok_paar)

# --- 2.5 linearisierte Krümmung der Produkte: Weyl gegen spurfreien Ricci, Sigma-Spannen aus z98
EI = np.linalg.inv(ETA)
EPSL = np.einsum("abcd,ai,bj,ck,dl->ijkl", LC, ETA, ETA, ETA, ETA)


def vierer(M):
    return np.array([0.5*np.trace(SIGMA[m] @ M) for m in range(4)])


def riemann(pl, h):
    """Linearisierte Riemannkrümmung einer ebenen Welle h e^(i p x) (untere Indizes)."""
    return 0.5*(np.einsum("m,s,nr->mnrs", pl, pl, h) + np.einsum("n,r,ms->mnrs", pl, pl, h)
                - np.einsum("m,r,ns->mnrs", pl, pl, h) - np.einsum("n,s,mr->mnrs", pl, pl, h))


def zerlege(R):
    Ric = np.einsum("mr,mnrs->ns", EI, R)
    Rs = np.einsum("ns,ns->", EI, Ric)
    g = ETA
    C = (R - 0.5*(np.einsum("ac,db->abcd", g, Ric) - np.einsum("ad,cb->abcd", g, Ric) - np.einsum("bc,da->abcd", g, Ric)
                  + np.einsum("bd,ca->abcd", g, Ric)) + (Rs/6)*(np.einsum("ac,db->abcd", g, g) - np.einsum("ad,cb->abcd", g, g)))
    sC = 0.5*np.einsum("abef,ei,fj,ijcd->abcd", EPSL, EI, EI, C)
    return C, Ric, Rs, sC


def im_spann_res(M, basis):
    Am = np.array([b.ravel() for b in basis]).T
    xs = np.linalg.lstsq(Am, M.ravel(), rcond=None)[0]
    return float(np.linalg.norm(Am @ xs - M.ravel())/np.linalg.norm(M))


SIGP = [np.array([[(SIGMA[m] @ EPS2 @ SIGMA[n].T - SIGMA[n] @ EPS2 @ SIGMA[m].T)[i, k] for n in range(4)] for m in range(4)])
        for i in range(2) for k in range(2)]                       # z98: eps auf dem Originalindex
SIGM = [np.array([[(SIGMA[m].T @ EPS2 @ SIGMA[n] - SIGMA[n].T @ EPS2 @ SIGMA[m])[j, l] for n in range(4)] for m in range(4)])
        for j in range(2) for l in range(2)]                       # z98: eps auf dem Partnerindex
kru = {}
ok_kru = True
for probe in range(3):
    lam, mu = zk(), zk()
    lt, mut = lam.conj(), mu.conj()
    res_ = {}
    for ph in (0.0, 0.4):
        lr, ltr = np.exp(1j*ph)*lam, np.exp(-1j*ph)*lt             # Modulphase (Referenzspinoren fest)
        p4 = vierer(np.outer(lr, ltr))
        em = vierer(np.outer(lr, mut)/ang(ltr, mut))               # Originalblatt: lambda lambda (Helizität -1)
        ep = vierer(np.outer(mu, ltr)/ang(mu, lr))                 # Partnerblatt: lambda~ lambda~ (Helizität +1)
        pl, eml, epl = ETA @ p4, ETA @ em, ETA @ ep
        for name, h in (("OO", np.outer(eml, eml)), ("PP", np.outer(epl, epl)),
                        ("OP", np.outer(eml, epl) + np.outer(epl, eml))):
            res_[(name, ph)] = riemann(pl, h)
        if ph == 0.0:
            Fo = np.outer(p4, em) - np.outer(em, p4)
            Fp = np.outer(p4, ep) - np.outer(ep, p4)
            spann = (im_spann_res(Fo, SIGM), im_spann_res(Fo, SIGP), im_spann_res(Fp, SIGP), im_spann_res(Fp, SIGM))
            pp = np.outer(pl, pl)
    zeile = {}
    for name in ("OO", "PP", "OP"):
        R0 = res_[(name, 0.0)]
        C, Ric, Rs, sC = zerlege(R0)
        nC = np.abs(C).max()
        dual = complex(np.vdot(C.ravel(), sC.ravel())/np.vdot(C.ravel(), C.ravel())) if nC > 1e-12 else 0j
        cR = complex(np.vdot(pp.ravel(), Ric.ravel())/np.vdot(pp.ravel(), pp.ravel()))
        Rph = res_[(name, 0.4)]
        ladung = float(np.angle(np.vdot(R0.ravel(), Rph.ravel())/np.vdot(R0.ravel(), R0.ravel()))/0.4)
        zeile[name] = {"|R|": float(np.abs(R0).max()), "|Ricci|": float(np.abs(Ric).max()), "R": float(abs(Rs)),
                       "|Weyl|": float(nC), "*C/C": str(np.round(dual, 9)), "Ricci - c p p": float(np.abs(Ric - cR*pp).max()),
                       "Modulladung": round(ladung, 9)}
    kru[probe] = {**zeile, "Spann F_O in Sigma_- / Sigma_+, F_P in Sigma_+ / Sigma_-": [float("%.1e" % v) for v in spann]}
    oo, pp_, op = zeile["OO"], zeile["PP"], zeile["OP"]
    ok_kru &= (oo["|Ricci|"] < 1e-12 and pp_["|Ricci|"] < 1e-12 and oo["|Weyl|"] > 1e-6 and pp_["|Weyl|"] > 1e-6
               and oo["*C/C"] != pp_["*C/C"] and abs(complex(oo["*C/C"]) + complex(pp_["*C/C"])) < 1e-9
               and abs(abs(complex(oo["*C/C"])) - 1) < 1e-9 and op["|Weyl|"] < 1e-12 and op["|Ricci|"] > 1e-6
               and op["R"] < 1e-12 and op["Ricci - c p p"] < 1e-12 and abs(oo["Modulladung"] - 4) < 1e-9
               and abs(pp_["Modulladung"] + 4) < 1e-9 and abs(op["Modulladung"]) < 1e-9
               and spann[0] < 1e-12 and spann[2] < 1e-12 and spann[1] > 0.5 and spann[3] > 0.5)
ZAHLEN["kruemmung"] = {str(k): v for k, v in kru.items()}
k0 = kru[0]
check(f"Linearisierte Krümmung der Doppelkopie-Produkte (drei reelle Nullimpulse, Polarisationen mit Referenzspinoren): "
      f"Originalblatt x Originalblatt (Helizität -1 x -1) ist reiner Weyl, Ricci = 0, *C = {k0['OO']['*C/C']} C, "
      f"Modulladung {k0['OO']['Modulladung']:.0f}; Partner x Partner ist reiner Weyl der anderen Dualität (*C = "
      f"{k0['PP']['*C/C']} C, Ladung {k0['PP']['Modulladung']:.0f}); Original x Partner hat Weyl = 0 "
      f"({k0['OP']['|Weyl|']:.1e}), R = 0 und Ricci = c p_mu p_nu != 0: reiner spurfreier Ricci Phi, modulneutral. Die "
      f"Feldstärke des Originalblatts liegt im Spann von z98s Sigma_- (eps auf dem Partnerindex), die des Partnerblatts im "
      f"Spann von Sigma_+ (Reste {k0['Spann F_O in Sigma_- / Sigma_+, F_P in Sigma_+ / Sigma_-']}): Gravitonen sind "
      f"Sigma_- x Sigma_- und Sigma_+ x Sigma_+, das gemischte Sigma_+ x Sigma_- ist Phi", ok_kru)
befund("E (gescoptes Negativ für die wörtliche Hypothese T): Original x Tomita-Partner je Bein ist nicht die "
       "Gravitations-Doppelkopie. Weil J die Helizität spiegelt (D -> -D), ist das je Bein J-gepaarte Produkt modulneutral: "
       "Es gibt den Dilaton-Axion-Sektor (Plebanski: spurfreier Ricci Phi = Sigma_+ x Sigma_-), für n = 4 eine "
       "helizitätsneutrale Amplitude und für n >= 5 exakt null, nie ein Graviton. Gerade Phi setzt die Vakuum-Einsteingleichung "
       "null", True)
befund("C (Hypothese T in blattdiagonaler Fassung T'): Die Doppelkopie ist blattdiagonal. Gravitonen sind Original x "
       "Original (Weyl Sigma_- x Sigma_-) und Partner x Partner (Sigma_+ x Sigma_+); die rechte Kopie eines Beins ist das "
       "Partnerbild des helizitätsgespiegelten Zustands, auf J-gepaarter Kinematik numerisch gleich Parke-Taylor. Die "
       "Tomita-Paarung Original x Partner tritt genau an drei Stellen auf: im Impuls P_mu = (1/2)(J c J) sigma c, im "
       "KLT-Kern s_ij = <ij>[ij] und in Hodges' Phi_ij = [ij]/<ij> (Poissonklammer des Partners über dem Propagator des "
       "Originals). Prämisse P_diag: die blattdiagonale Paarung (gleiche Helizität links und rechts je Bein), nicht aus J, "
       "das nur die gespiegelte liefert", True)
laufzeit()

# =================================================================== Teil 3
say("\nTeil 3: Normierung von G")

# --- 3.1 Dimension der Doppelkopie-Kopplung
rho = 1.7
expo = {}
ok_dim = True
for n in (3, 4, 5, 6):
    lam, lt = kin_drei() if n == 3 else kin_komplex(n)
    A, S = klammern(lam, lt)
    A2, S2 = klammern(rho*lam, rho*lt)                             # p -> rho^2 p (Massendimension 1 <-> rho^2)
    ePT = math.log(abs(pt(list(range(n)), A2)/pt(list(range(n)), A)))/math.log(rho**2)
    eH = math.log(abs(hodges(rho*lam, A2, S2, *STD, XR, YR)/hodges(lam, A, S, *STD, XR, YR)))/math.log(rho**2)
    mitte = list(range(1, n - 2))
    sg = tuple(mitte)
    eS = (math.log(abs(kern(list(sg), list(sg), A2*S2)/kern(list(sg), list(sg), A*S)))/math.log(rho**2)) if mitte else 0.0
    expo[n] = {"[A_n]": round(ePT, 9), "[M_n]": round(eH, 9), "[S]": round(eS, 9)}
    ok_dim &= abs(ePT - (4 - n)) < 1e-9 and abs(eH - 2) < 1e-9 and abs(eS - 2*(n - 3)) < 1e-9
ZAHLEN["dimension"] = {"Exponenten (Massendimension)": expo, "volle Amplituden": "A = g^(n-2) A, M = (kappa/2)^(n-2) M, "
                       "[M_voll] = 4 - n", "L = kappa/(2 g^2)": "Massendimension -1"}
check(f"Dimension der Doppelkopie-Kopplung: Unter p -> rho^2 p skalieren Parke-Taylor mit Massendimension 4 - n, Hodges' M_n "
      f"mit 2 und der Impulskern mit 2(n - 3) (n = 3..6 exakt: {expo}), also KLT(A, A) mit 2 wie M_n. Mit den vollen "
      f"Kopplungen (A g^(n-2), M (kappa/2)^(n-2), beide mit Massendimension 4 - n) gilt M_voll = (kappa/(2 g^2))^(n-2) "
      f"KLT(A_voll, A_voll): Die Doppelkopie fixiert nur das Verhältnis L = kappa/(2 g^2), und L ist eine Länge. Ohne "
      f"eine Längenskala (in Strings alpha') folgt kappa nicht aus g", ok_dim)

# --- 3.2 native Länge: Dilatation, Vakuum, |I|-Normierung
def gewicht(M, X):
    C = M @ X - X @ M
    return complex(np.vdot(X.ravel(), C.ravel())/np.vdot(X.ravel(), X.ravel())), float(np.abs(C - (np.vdot(X.ravel(), C.ravel())/np.vdot(X.ravel(), X.ravel()))*X).max())


wP = [gewicht(S4, X) for X in MP]
wK = [gewicht(S4, X) for X in MKL]
wE = gewicht(S4, MEt)
ok_w = (all(abs(w[0] - wP[0][0]) < 1e-12 and w[1] < 1e-12 for w in wP) and all(abs(w[0] + wP[0][0]) < 1e-12 and w[1] < 1e-12
                                                                               for w in wK)
        and abs(wP[0][0]) > 0.5 and wE[1] > 1e-3)
DilF = sum((dag(DD[j]) @ dag(EE[j]) - EE[j] @ DD[j] for j in range(2)), sps.csr_matrix((DIMV, DIMV))).tocsr()
ok_adc = all(klein(komm(DilF, CC[j]) + CC[j]) and klein(komm(DilF, WO[j]) - WO[j]) for j in range(2))
ok_E = np.allclose(Eop @ OMV, OMV)
dil_vak = float(np.linalg.norm(DilF @ OMV))
sq = 0.1
Om_s = spla.expm_multiply((sq*DilF).tocsc(), OMV)
cc_s = sum(np.vdot(CC[j] @ Om_s, CC[j] @ Om_s) for j in range(2)).real/2
ok_sq = abs(cc_s - math.exp(2*sq)) < 1e-6 or abs(cc_s - math.exp(-2*sq)) < 1e-6
ZAHLEN["native_laenge"] = {"Dilatationsgewicht von P, K (S4)": [r(wP[0][0].real, 9), r(wK[0][0].real, 9)],
                           "[Dil, E] nicht proportional zu E (Rest)": wE[1], "||Dil Omega||": dil_vak,
                           "<c^+ c> im gestauchten Vakuum (s = 0,1)": r(cc_s, 9), "e^(2 s)": r(math.exp(2*sq), 9),
                           "l (z98)": "R_Vak/2", "R_Vak in CCR-Einheiten (m = 1)": KAP2_NATIV[1]}
check(f"Native Länge: Die Dilatation (S4 auf T) gibt P_mu das Gewicht {wP[0][0].real:+.0f} und K_mu {wK[0][0].real:+.0f}, "
      f"die kompakte Energie E = (P_0 + K_0)/2 hat kein Gewicht (Rest {wE[1]:.2f}). Auf Fock {NV}^4 ist der native Erzeuger "
      f"Dil = sum (d^+ e^+ - e d) mit [Dil, c] = -c und [Dil, omega''] = omega'' (Zweimodenstauchung); das KMS-Vakuum ist "
      f"Grundzustand von E (E Omega = Omega), aber nicht dilatationsinvariant (||Dil Omega|| = {dil_vak:.4f}), und im "
      f"gestauchten Vakuum ist <c^+ c> = {cc_s:.6f} = e^(2 s) statt 1. Die KMS-Normierung |I| = 1 und die Vakuumgröße R_Vak "
      f"sind also ein und dasselbe Skalendatum; eine zweite native Länge gibt es nicht (z98: l = R_Vak/2)",
      ok_w and ok_adc and ok_E and dil_vak > 1 and ok_sq)

# --- 3.3 Prämissenabtastung gegen das Lambda-G-Ziel
lg = lambda v: float(mp.log10(v))
ziel = mp.e**(-2/(1/ALPHA_INV))/(128*mp.pi**3)
g2_kand = {"1/k": mp.mpf(1), "2/k": mp.mpf(2), "2 pi/k": 2*mp.pi, "4 pi/k": 4*mp.pi}         # Niveau k = 1 (z99)
L_kand = {"l = R_Vak/2": mp.mpf(1)/2, "R_Vak": mp.mpf(1), "CCR-Einheit R_Vak/tanh(pi)": 1/mp.tanh(mp.pi)}
tab = {}
for (gn, g2), (ln, Lr) in itertools.product(g2_kand.items(), L_kand.items()):
    GR = g2**2*Lr**2/(8*mp.pi)                                    # G = kappa^2/(32 pi), kappa = 2 g^2 L
    tab[f"g^2 = {gn}, L = {ln}"] = {"log10 G/R_Vak^2": round(lg(GR), 4), "R_Vak/l_Pl": round(float(1/mp.sqrt(GR)), 4)}
lgs = [v["log10 G/R_Vak^2"] for v in tab.values()]
rpl = [v["R_Vak/l_Pl"] for v in tab.values()]
L_noetig = {gn: round(lg(mp.sqrt(8*mp.pi*ziel)/g2), 4) for gn, g2 in g2_kand.items()}
R_ziel = float(mp.log10(1/mp.sqrt(ziel)))
S_dS = {"P_L (min, max)": [round(float(mp.pi/4/(10**mp.mpf(max(lgs)))), 4), round(float(mp.pi/4/(10**mp.mpf(min(lgs)))), 4)],
        "Ziel log10 pi l^2/G": round(lg(mp.pi/4/ziel), 4)}
ZAHLEN["prämissen_G"] = {"Ziel log10 G/R_Vak^2": round(lg(ziel), 4), "Kandidaten": tab,
                         "Spanne log10 G/R_Vak^2": [min(lgs), max(lgs)], "Spanne R_Vak/l_Pl": [min(rpl), max(rpl)],
                         "Ziel log10 R_Vak/l_Pl": round(R_ziel, 4), "nötiges log10 L/R_Vak je g^2": L_noetig,
                         "de-Sitter-Lesart pi l^2/G": S_dS}
check(f"Prämissenabtastung (Prämisse P_L: L = kappa/(2 g^2) ist eine native Länge mal O(1), g^2 aus dem Niveau k = 1 in "
      f"vier Normierungen, L in drei nativen Fassungen): G/R_Vak^2 = g^4 (L/R_Vak)^2/(8 pi) liegt zwischen 10^{min(lgs):.2f} "
      f"und 10^{max(lgs):.2f}, also R_Vak = {min(rpl):.2f} bis {max(rpl):.2f} Planck-Längen. Das Lambda-G-Ziel verlangt "
      f"G/R_Vak^2 = e^(-2/alpha)/(128 pi^3) = 10^{lg(ziel):.2f} (R_Vak = 10^{R_ziel:.2f} l_Pl), also L/R_Vak = 10^"
      f"{min(L_noetig.values()):.2f} bis 10^{max(L_noetig.values()):.2f}: eine zweite, um rund 60 Größenordnungen kleinere "
      f"Länge, die nativ nicht vorkommt. Abstand der nativen Fassungen zum Ziel mindestens {min(lgs) - lg(ziel):.1f} "
      f"Größenordnungen", max(lgs) < 1.5 and min(lgs) > -3 and min(lgs) - lg(ziel) > 100 and max(L_noetig.values()) < -55)
befund("E/O: Die Doppelkopie fixiert kappa nicht. Sie bestimmt nur L = kappa/(2 g^2), eine Länge; nativ gibt es genau eine "
       "Länge, das Vakuum (die KMS-Normierung von |I| ist dasselbe Datum), und der Koeffizient des Poissonvertex relativ zur "
       "Fock-Normierung bleibt frei (z89, z98). Mit jeder nativen O(1)-Fassung (Prämisse P_L) ist das Vakuum planckgroß; das "
       "Lambda-G-Ziel G/R_Vak^2 = 10^(-122,6) ist damit unverträglich und verlangt eine zweite Länge L ~ 10^(-60,6) R_Vak, "
       "also e^(-1/alpha) als nichtperturbativen Faktor, den keine Doppelkopie-Normierung enthält. Dazu bleibt das native "
       "Vorzeichen AdS (z98)", True)
befund(f"L: Holografisch gelesen ist pi l^2/G die Zahl der Freiheitsgrade (de-Sitter-Entropie, z96, z109). Unter P_L liegt "
       f"sie zwischen {S_dS['P_L (min, max)'][0]:.2f} und {S_dS['P_L (min, max)'][1]:.1f} für eine Naht; das Ziel verlangt "
       f"10^{S_dS['Ziel log10 pi l^2/G']:.2f}, also ebenso viele native Freiheitsgrade (vgl. z88: viele Teilchen brauchen "
       f"viele Nähte). Das ist eine Zählforderung, keine Herleitung", True)
laufzeit()

# =================================================================== Teil 4
say("\nTeil 4: Bilanz und Urteil")
STATUS = {
    "Winkelklammer": ("E", "<ij> = *I_pi(Z_i, Z_j) = eps(c-Eigenwerte), Originalindex, Koeffizient 1, translationsinvariant"),
    "eckige Klammer": ("E", "[ij] = B_(I_pi) = eps(c^+-Eigenwerte) = eps des Partners (J-Bild), Koeffizient 1"),
    "I_omega, AdS-Büschel": ("E", "*I_omega = -det(x) <ij> auf L_x (K-Eigenzustände); *(I_pi - I_omega) = Omega(x) <ij>, "
                                  "Omega = 1 + det x (z98), null am Vakuumpunkt"),
    "Normierung": ("E", "s_ij = <ij>[ij] = 2 p_i.p_j mit nativem P_mu, kein freier Faktor; Modulphase = globale Kleingruppe"),
    "Hodges MHV": ("E (Vergleichsziel)", "n = 3..6 eindeutig (Streichungen, Referenzspinoren), Negativkontrolle ohne "
                                          "Impulserhaltung"),
    "Hodges = KLT": ("E (Vergleichsziel)", "Impulskern mit Parke-Taylor, Verhältnis 1 für n = 3..6, komplex und nativ reell"),
    "Ladungen": ("E", "je Bein +-4 = 2 x Eich, global Q(M_n) = 16 - 4n (n = 3: e^(4 i theta) wie z98, n = 4: 0)"),
    "Original x Partner je Bein": ("E, gescoptes Negativ", "modulneutral: Phi (Dilaton-Axion), n = 4 helizitätsneutral, "
                                                            "n >= 5 null; nie Graviton"),
    "blattdiagonale Doppelkopie": ("C", "Original x Original, Partner x Partner; Tomita-Paarung im Impuls, im KLT-Kern und in "
                                        "Phi_ij; Prämisse P_diag"),
    "kappa, G": ("O", "Doppelkopie fixiert nur L = kappa/(2 g^2) (Länge); native Länge nur R_Vak; P_L gibt planckgroßes "
                      "Vakuum, unverträglich mit dem Lambda-G-Ziel 10^(-122,6)"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
urteil = (
    "Hypothese T ist in ihrer wörtlichen Form gescopt verneint, in blattdiagonaler Fassung bedingt tragfähig; G bleibt "
    "offen. (1) Native Klammern [E]: Die Impulseigenzustände sind die gemeinsamen Eigenzustände der kommutativen c-Algebra "
    "(lambda = c-Eigenwert, p_mu = (1/2) lambda^+ sigma_mu lambda), der Partnerspinor ist das J-Bild (J c J = c^+). Ein "
    "einziger nativer Bivektor I_pi liefert beide Klammern mit Koeffizient 1: unten (Hodge) <ij> = *I_pi(Z_i, Z_j) = "
    "eps(c_i, c_j) auf dem Originalindex, oben (Poisson) [ij] = eps(c_i^+, c_j^+), die eps-Form des Tomita-Partners. "
    "I_omega gibt -det(x) <ij> auf der Geraden L_x und gehört zu den K-Eigenzuständen; das native AdS-Büschel I_pi - I_omega "
    "gibt Omega(x) <ij> mit Omega = 1 + det x, z98s konformen Faktor. Mit dem KMS-normierten c gilt s_ij = <ij>[ij] = "
    "2 p_i.p_j ohne freien Faktor; die Modulphase ist die globale Kleingruppenphase (<ij> Ladung 2, [ij] Ladung -2). Für "
    "n = 3 braucht es komplexe, nicht J-gepaarte Kinematik. (2) Gravitation MHV [E als Vergleichsziel]: Mit diesen "
    "Klammern ist Hodges' reduzierte Determinante für n = 3..6 eindeutig (alle Streichungen, Referenzspinoren; ohne "
    "Impulserhaltung nicht) und gleich der KLT-Doppelkopie aus Parke-Taylor (Verhältnis 1, komplex und nativ reell). Die "
    "Ladungen verdoppeln sich je Bein (+-4 = 2 x +-2), global Q(M_n) = 16 - 4n, also e^(4 i theta) für den kubischen "
    "Vertex wie in z98 und Phasenfreiheit bei n = 4. Native Frage: Die J-Paarung desselben Zustands (rechte Kopie = conj der "
    "linken) ist je Bein modulneutral und gibt den Dilaton-Axion-Sektor, für n = 4 eine reelle, helizitätsneutrale "
    "Amplitude vom Betrag der Gravitonamplitude, für n >= 5 exakt null; linearisiert ist Original x Partner reiner spurfreier Ricci Phi = Sigma_+ x Sigma_- (Weyl = 0), "
    "genau die Komponente, die die Vakuum-Einsteingleichung null setzt. Gravitonen sind blattdiagonal: Original x Original "
    "(Weyl in Sigma_- x Sigma_-) und Partner x Partner (Sigma_+ x Sigma_+); die rechte Kopie eines Beins ist das "
    "Partnerbild des helizitätsgespiegelten Zustands, auf J-gepaarter Kinematik gleich Parke-Taylor. Die Tomita-Paarung "
    "tritt im Impuls, im KLT-Kern s_ij und in Phi_ij = [ij]/<ij> auf, nicht als zweite Eichkopie. Die Original-Partner-"
    "Struktur realisiert also nicht genau diese Doppelkopie, sondern ihren Kern; die Doppelkopie selbst braucht P_diag "
    "(gleiche Helizität links und rechts je Bein). (3) Normierung von G [E/O]: Die Doppelkopie fixiert nur L = kappa/(2 g^2), "
    "eine Länge (Skalierung exakt). Nativ gibt es eine Länge, das Vakuum; die KMS-Normierung von |I| ist dasselbe Datum (im "
    "gestauchten Vakuum ist <c^+ c> = e^(2 s)), und der Vertexkoeffizient bleibt frei. Unter jeder nativen O(1)-Prämisse P_L "
    "(g^2 aus dem Niveau 1, L = R_Vak/2, R_Vak oder die CCR-Einheit) ist G/R_Vak^2 = 10^(-2,0) bis 10^(0,8), das Vakuum also "
    "planckgroß; das Lambda-G-Ziel 10^(-122,6) verlangt eine zweite Länge L ~ 10^(-60,6) R_Vak, die nicht nativ ist, und "
    "dazu dS statt des nativen AdS. Verträglich wäre P_L nur, wenn R_Vak die Planck- und nicht die Hubbleskala ist. "
    "(4) Folgen: K13 bleibt L; die absolute Normierung der Metrik ist auf ein einziges fehlendes Datum reduziert, die Länge "
    "L = kappa/(2 g^2) (gleichwertig b bzw. G), die nicht aus der Doppelkopie folgt. Szenarioachse Gravitation: neuer Wert "
    "'blattdiagonale Doppelkopie mit Tomita-Kern' (bedingt, P_diag und G frei bzw. P_L); 'Original x Partner' ist als "
    "Gravitonquelle ausgeschlossen und gehört zum Dilaton-Axion-Wert. T7: MHV-Baumgravitation ist aus nativen Klammern "
    "reproduziert, als Kinematik ohne native Dynamik; T6: kappa bzw. G ist keine Folge von g und Niveau. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z112_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
