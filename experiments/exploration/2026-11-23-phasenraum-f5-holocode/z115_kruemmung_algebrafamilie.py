"""Z115: Krümmung aus bewegter Algebra bei festem Zustand: die Familie J(M_g, psi) (tiefe Welle T6).

Leitlinie 'Rekonstruktion statt Auswahl'. Firewalls: GRAV.NONCIRCULAR.01 (der Modulfluss bleibt Uhr; die Tomita-Familie ist
Kinematik, keine Zeit), GRAV.PARENT.GAP.NO_GO.01, DIMENSION.UPLIFT.FIREWALL.01 (Fock- und Twistorkinematik einer Naht bzw.
der Saaten von z101 auf Abschneidungen; keine wechselwirkende 3+1-Schließung). Gelöste Blöcke (alpha, Flavour, Massen, 16,
E8) bleiben unberührt; nichts wird befördert. Bestand: z98 (Bahnfamilie J_g = U(g) J U(g)^-1 = AdS4, Levi-Civita, F = 0),
z100 (Maurer-Cartan-Lemma: jede Bahn U(g) J_psi U(g)^-1 hat F_MM = 0; Empfehlung: J(M_g, psi) bei festem Zustand), z101
(kollektiver Zustand N = 6, Standardform diagonal, J zustandsunabhängig). Statuskürzel: E exakt, C bedingt, L Lesart, O offen.

Familie: M_g = U(g) M U(g)^*, psi fest, J(M_g, psi) = U(g) J(M, U(g)^* psi) U(g)^* (Kovarianz). Für gaußsche,
ladungserhaltende Daten ist alles Einteilchenalgebra auf T = span{Z}, Z = (d_1, d_2, e_1^+, e_2^+): die Algebra ist die
positive 2-Ebene ihrer Vernichter a, die Kommutante die negative Ebene b^+, der Zustand die positive Ebene seiner Vernichter;
tau = (Ad J)^+ ist eine lineare Involution auf T (z98).

Teil 1  Vakuum: Einteilchenformel für tau(M, psi) aus Algebra- und Zustandsebene, Fockkontrolle (Polarteil, J Omega = Omega);
        Lage der Familie (Menge Conj_Omega der Konjugationen, die Omega festhalten), Vergleich mit der Bahnfamilie, Kern und
        Rang, globale Karte (tan theta = coth(2r) sinh rho), Äquivarianz, Verschränkungshorizont; nativer Rahmen (Invariante
        iota der nativen Gruppe), eigenes Modell (Maurer-Cartan, Berry-Krümmung, Torsion, Signatur, Schnittkrümmungen);
        quadratische Antwort (BKM-Fisher der reduzierten Zustände).
Teil 2  Nicht invariante Zustände: Kondensat (a1) (Verschiebungslemma), D = -4-Kondensat (a2) (Fock: Fixpunkt, Rang,
        lambda-Skalierung, Fisher), KMS-Zustand des Originals bei beta = 4 pi (natives q = e^(-2 pi)), kollektiver Zustand
        von z101 (Saatraum, Qubit-Lesart, Paargruppe je Saat).
Teil 3  Grenzfall q -> 1 (Konformität der Karte AdS4 -> U(2)), Bilanz, Urteil; Folgen für K13 und die Gravitationsachse
        der Szenario-Suche.
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
rng = np.random.default_rng(115)
I2 = np.eye(2)
SIGMA = [np.eye(2, dtype=complex), np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]),
         np.diag([1., -1.]).astype(complex)]
Q_TEST = 0.25                                                       # Prüfwert wie z100
Q_NATIV = math.exp(-2*math.pi)                                      # KMS des Originals bei beta = 4 pi (q^2 = e^(-4 pi))


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


def kl(A, B):
    return A @ B - B @ A


# ------------------------------------------------------------------ u(2,2) und die native so(3,2) auf T (wie z98, z100)
G4 = np.diag([1., 1., -1., -1.]).astype(complex)
S4 = np.block([[np.zeros((2, 2)), I2], [I2, np.zeros((2, 2))]]).astype(complex)       # tau_0 = (Ad J_Omega)^+ auf Z
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
Ion = (np.outer(vvec[0], vvec[1]) - np.outer(vvec[1], vvec[0])).astype(complex)       # I_omega
NP = [-0.5*sum(s[i, j]*np.outer(uvec[i], vvec[j]) for i in range(2) for j in range(2)) for s in SIGMA]
NK = [-0.5*sum(s[i, j]*np.outer(vvec[i], uvec[j]) for i in range(2) for j in range(2)) for s in SIGMA]
ETA = np.diag([1., -1., -1., -1.])
MP = [1j*N for N in NP]
MKL = [ETA[m, m]*1j*NK[m] for m in range(4)]
LOR = teilalgebra([lambda M: M @ S4 - S4 @ M, lambda M: M @ Inp + Inp @ M.T])
QP = [(MP[m] + MKL[m])/2 for m in range(4)]                         # Q^+_mu, Q_0 = E (z98)
SO32 = LOR + QP
ROT = teilalgebra([lambda M: kl(M, S4), lambda M: M @ Inp + Inp @ M.T, lambda M: kl(M, G4)])
BOOST = teilalgebra([lambda M: kl(M, S4), lambda M: M @ Inp + Inp @ M.T, lambda M: M @ G4 + G4 @ M])
KBIG = teilalgebra([lambda M: kl(M, G4)])                           # u(2) + u(2): Stabilisator der Vakuumebene
ISYM = Inp - Ion                                                    # so(3,2) = Der(I_pi - I_omega) (z98)
ISYM_INV = np.linalg.inv(ISYM)


def iota(tau):
    """Invariante der nativen Gruppe auf Involutionen: tr(tau I tau^T I^-1), I = I_pi - I_omega."""
    return complex(np.trace(tau @ ISYM @ tau.T @ ISYM_INV))


def Rw(wort):
    """U(g) = U(e^(t1 X1)) U(e^(t2 X2)) ...; U(g) Z U(g)^* = R Z mit R = e^(tn Xn) ... e^(t1 X1)."""
    R = np.eye(4, dtype=complex)
    for X, t in wort:
        R = sla.expm(t*X) @ R
    return R


def ebenen(q):
    s = math.sqrt(1 - q*q)
    A0 = np.hstack([I2, q*I2]).astype(complex)/s                     # a = (d + q e^+)/s
    B0 = np.hstack([q*I2, I2]).astype(complex)/s                     # b^+ = (e^+ + q d)/s
    return A0, B0


P_OM = np.hstack([I2, 0*I2]).astype(complex)                        # Vernichter von Omega in T: d_1, d_2


def polar_u(L):
    U_, sv, Wh = np.linalg.svd(L)
    return U_ @ Wh, sv


def tomita1(A, B, P):
    """Einteilchen-Tomita: Algebra mit Vernichtern a = A Z, Kommutante b^+ = B Z, Zustand mit Vernichtern P Z.
    P = X A + Y B, a psi = Lam b^+ psi mit Lam = -X^-1 Y, V = Polarteil von Lam, tau(a) = V b^+, tau(b^+) = V^+ a."""
    T = np.vstack([A, B])
    XY = P @ np.linalg.inv(T)
    Lam = -np.linalg.solve(XY[:, :2], XY[:, 2:])
    V, sv = polar_u(Lam)
    tau = np.linalg.inv(T) @ np.block([[np.zeros((2, 2)), V], [V.conj().T, np.zeros((2, 2))]]) @ T
    return tau, V, sv, Lam


def Phi(R, q=Q_TEST, P=P_OM):
    """tau(M_g, psi) bei festem Zustand (Algebra bewegt: Vernichter A0 R, B0 R)."""
    A0, B0 = ebenen(q)
    return tomita1(A0 @ R, B0 @ R, P)[0]


def Phi_bahn(R):
    """Bahnfamilie z98: tau(U(g) J U(g)^-1) = R^-1 S4 R."""
    return np.linalg.inv(R) @ S4 @ R


# ------------------------------------------------------------------ Fockraum einer Naht (z100-Maschinerie)
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


a1, a2, b1, b2 = moden(NV, 4)
AO, BO = [a1, a2], [b1, b2]
DIM, DA = NV**4, NV**2
OCC = np.array(list(itertools.product(range(NV), repeat=4)))
TOT = OCC.sum(1)
OCCA = np.array(list(itertools.product(range(NV), repeat=2))).sum(1)
NULL = sps.csr_matrix((DIM, DIM), dtype=complex)


def fock_daten(qq):
    sq = math.sqrt(1 - qq*qq)
    d_ = [((AO[i] - qq*dag(BO[i]))/sq).tocsr() for i in range(2)]
    e_ = [((BO[i] - qq*dag(AO[i]))/sq).tocsr() for i in range(2)]
    ZT = [d_[0], d_[1], dag(e_[0]), dag(e_[1])]
    ZTd = [dag(x) for x in ZT]

    def fock_gen(M):
        Am = 1j*G4 @ M
        H = NULL.copy()
        for al in range(4):
            for be in range(4):
                if abs(Am[al, be]) > 1e-14:
                    H = H + Am[al, be]*(ZTd[al] @ ZT[be])
        return H.tocsr()
    om1 = np.diag([sq*qq**n for n in range(NV)])
    C_OM = np.kron(om1, om1)
    C_OM /= np.linalg.norm(C_OM)
    return d_, e_, fock_gen, C_OM.reshape(-1).astype(complex)


d_, e_, fock_gen, OMEGA = fock_daten(q)
HS = [fock_gen(M) for M in SO32]
ERZ = [dag(d_[0]), dag(d_[1]), dag(e_[0]), dag(e_[1])]
aa2 = moden(NV, 2)
lowa = np.where(OCCA <= 4)[0]


def polar(C):
    U_, sv, Wh = np.linalg.svd(C)
    return U_ @ Wh, sv, U_


def gamma1(V):
    """Fockoperator Gamma(V) auf dem Originalraum (Gamma(V) a_k^+ Gamma(V)^* = sum_i V_ik a_i^+)."""
    L = sla.logm(V)
    G = sum((L[i, j]*(dag(aa2[i]) @ aa2[j]) for i in range(2) for j in range(2)), sps.csr_matrix((DA, DA)))
    return sla.expm(G.toarray())


def J_auf(H, V, v, t=1.0):
    """J(M_g, psi) v = U (V (x) V^T) SWAP conj(U^* v), U = e^(i t H), V Polarteil von U^* psi (antilinear in v)."""
    w = spla.expm_multiply(-1j*t*H, v) if H is not None else v
    W = V @ np.conj(w).reshape(DA, DA).T @ V
    w2 = W.reshape(-1)
    return spla.expm_multiply(1j*t*H, w2) if H is not None else w2


def V_von(psi, H=None, t=1.0):
    w = spla.expm_multiply(-1j*t*H, psi) if H is not None else psi
    return polar(w.reshape(DA, DA))[0]


laufzeit()

# =================================================================== Teil 1
say("Teil 1: Vakuum -- die Familie J(M_g, Omega) bei festem Zustand")

# --- 1.1 Einteilchenformel und Fockkontrolle
say("  1.1 Einteilchenformel tau(M, psi) aus Algebra- und Zustandsebene, Fockkontrolle (NV = 11 je Mode, q = 0,25)")
A0, B0 = ebenen(q)
tau0, V00, sv00, _ = tomita1(A0, B0, P_OM)
ok_tau0 = np.abs(tau0 - S4).max() < 1e-14 and np.allclose(sv00, [q, q])
WORTE = {"Q1 (t = 0,4)": [(QP[1], 0.4)], "Boost (t = 0,3)": [(BOOST[0], 0.3)], "E (t = 0,5)": [(QP[0], 0.5)],
         "Zufall": [(sum(c*X for c, X in zip(rng.normal(size=10), SO32))*0.25, 1.0)]}
kov, pol, fix = {}, {}, {}
for name_, wort in WORTE.items():
    X, t = wort[0]
    R = Rw(wort)
    # Kovarianz: tau(M_g, Omega) = R^-1 tau(M, U^* Omega) R, U^* Omega hat die Vernichter P_Om R^-1
    t_dir = Phi(R)
    t_kov, V1, _, _ = tomita1(A0, B0, P_OM @ np.linalg.inv(R))
    kov[name_] = float(np.abs(t_dir - np.linalg.inv(R) @ t_kov @ R).max())
    H = fock_gen(X)
    VF = V_von(OMEGA, H, t)
    pol[name_] = float(np.abs((VF - VF[0, 0]*gamma1(V1))[np.ix_(lowa, lowa)]).max())     # globale Phase von U^* Omega
    fix[name_] = float(np.linalg.norm(J_auf(H, VF, OMEGA, t) - OMEGA))
ZAHLEN["einteilchen_fock"] = {"|tau(M, Omega) - S4|": float(np.abs(tau0 - S4).max()), "Kovarianz (Einteilchen)": kov,
                              "|V_Fock - Gamma(V_1)| (tiefer Block N_a <= 4)": pol, "|J(M_g, Omega) Omega - Omega|": fix}
check(f"[E] Einteilchenformel: Für gaußsche, ladungserhaltende Daten ist tau = (Ad J)^+ auf T durch zwei Ebenen bestimmt: "
      f"Algebra (Vernichter a = A Z, Kommutante b^+ = B Z) und Zustand (Vernichter P Z). Aus P = X A + Y B folgt a psi = "
      f"Lam b^+ psi, Lam = -X^-1 Y; mit dem Polarteil V von Lam ist tau(a) = V b^+, tau(b^+) = V^+ a. Für (M, Omega) gibt das "
      f"tau = S4 und Lam = q 1 (Abweichung {np.abs(tau0 - S4).max():.0e}); die Kovarianz tau(M_g, Omega) = R^-1 tau(M, U^* "
      f"Omega) R gilt exakt (max {max(kov.values()):.0e}). Fockkontrolle (vier Gruppenelemente, auch ein Zufallselement): Der "
      f"Polarteil von C(U(g)^* Omega) ist bis auf die globale Phase von U(g)^* Omega gleich Gamma(V_1) (max "
      f"{max(pol.values()):.1e} im tiefen Block, Abschneidefehler ~ q^NV; die Phase ändert J nur um einen Faktor, tau = "
      f"(Ad J)^+ nicht), und "
      f"J(M_g, Omega) hält Omega fest (max {max(fix.values()):.1e})",
      ok_tau0 and max(kov.values()) < 1e-12 and max(pol.values()) < 1e-6 and max(fix.values()) < 1e-7)
laufzeit()

# --- 1.2 Lage der Familie, Vergleich mit der Bahnfamilie, Kern, Karte, Äquivarianz, Horizont
say("  1.2 Lage der Familie: Conj_Omega = U(2), Vergleich mit der Bahnfamilie, globale Karte, Verschränkungshorizont")
form = []
for _ in range(20):
    R = Rw([(sum(c*X for c, X in zip(rng.normal(size=10), SO32))*0.5, 1.0)])
    t_ = Phi(R)
    be = t_[:2, 2:]
    form.append(max(np.abs(t_[:2, :2]).max(), np.abs(t_[2:, 2:]).max(), np.abs(be @ be.conj().T - I2).max(),
                    np.abs(t_[2:, :2] - be.conj().T).max(), np.abs(t_ @ t_ - np.eye(4)).max()))
vgl = {}
for name_, X in [("E", QP[0]), ("Drehung", ROT[0]), ("Boost", BOOST[0]), ("Q1", QP[1]), ("Q2", QP[2]), ("Q3", QP[3])]:
    vgl[name_] = [float(np.abs(Phi(Rw([(X, t)])) - Phi_bahn(Rw([(X, t)]))).max()) for t in (0.3, 0.45)]


def dPhi(R0, X, h=1e-6, q_=Q_TEST):
    """Ableitung entlang g -> g e^(eps X) (rechts; in der Zeilendarstellung R -> e^(eps X) R)."""
    return (Phi(sla.expm(h*X) @ R0, q_) - Phi(sla.expm(-h*X) @ R0, q_))/(2*h)


def rang_sv(R0, q_=Q_TEST):
    D = np.array([np.concatenate([dPhi(R0, X, q_=q_).real.ravel(), dPhi(R0, X, q_=q_).imag.ravel()]) for X in SO32])
    return np.linalg.svd(D, compute_uv=False)


sv_eins = rang_sv(np.eye(4))
kern_lor = max(np.abs(dPhi(np.eye(4), X)).max() for X in LOR)
raenge = [int(np.sum(sv_eins > 1e-6*sv_eins[0]))]
for _ in range(3):
    sv_ = rang_sv(Rw([(sum(c*X for c, X in zip(rng.normal(size=10), SO32))*0.5, 1.0)]))
    raenge.append(int(np.sum(sv_ > 1e-6*sv_[0])))
# globale Karte: g = e^(phi E) e^(rho n.Q)
c_q = (1 + q*q)/(2*q)
karte = []
for _ in range(3):
    n = rng.normal(size=3)
    n /= np.linalg.norm(n)
    Xn = sum(n[i]*QP[i + 1] for i in range(3))
    for rho in (0.2, 0.9, 2.5):
        be = Phi(Rw([(Xn, rho)]))[:2, 2:]
        th = float(np.abs(np.angle(np.linalg.eigvals(be))).max())
        karte.append((abs(math.tan(th) - c_q*math.sinh(rho))/(c_q*math.sinh(rho)), abs(np.linalg.det(be) - 1)))
phase = max(float(np.abs(Phi(Rw([(QP[0], ph), (QP[2], 0.8)]))[:2, 2:] - np.exp(1j*ph)*Phi(Rw([(QP[2], 0.8)]))[:2, 2:]).max())
            for ph in (0.3, 1.1, 2.9))
rand_th = [float(np.abs(np.angle(np.linalg.eigvals(Phi(Rw([(QP[3], r_)]))[:2, 2:]))).max()/np.pi) for r_ in (4.0, 8.0, 16.0)]
rand_sv = [float(tomita1(A0, B0, P_OM @ np.linalg.inv(Rw([(QP[3], r_)])))[2].min()) for r_ in (4.0, 8.0, 16.0)]
# Äquivarianz: links unter K = Stab(Omega) (E, Drehungen), rechts nicht unter Boosts und E
R_z = Rw([(sum(c*X for c, X in zip(rng.normal(size=10), SO32))*0.4, 1.0)])
links = max(float(np.abs(Phi(R_z @ Rw([(X, 0.3)])) - np.linalg.inv(Rw([(X, 0.3)])) @ Phi(R_z) @ Rw([(X, 0.3)])).max())
            for X in (QP[0], ROT[0], ROT[1]))
rechts_rot = float(np.abs(Phi(Rw([(ROT[1], 0.3)]) @ R_z) - Phi(R_z)).max())
rechts_boost = float(np.abs(Phi(Rw([(BOOST[1], 0.3)]) @ R_z) - Phi(R_z)).max())
rechts_E = float(np.abs(Phi(Rw([(QP[0], 0.3)]) @ R_z) - Phi(R_z)).max())
# Verschränkungshorizont: Boost mit Schnelligkeit 2 artanh(q) macht eine Mode von M_g unverschränkt
t_stern = 2*math.atanh(q)
sv_hor = lambda t: tomita1(A0, B0, P_OM @ np.linalg.inv(Rw([(BOOST[0], t)])))[2]
hor = {"t* = 2 artanh q": t_stern, "kleinster Schmidtwert bei t*": float(sv_hor(t_stern).min()),
       "bei 0,9 t*": float(sv_hor(0.9*t_stern).min()), "|beta - 1| bei 0,9 t*": float(np.abs(Phi(Rw([(BOOST[0], 0.9*t_stern)]))[:2, 2:] - I2).max()),
       "|beta - 1| bei 1,1 t*": float(np.abs(Phi(Rw([(BOOST[0], 1.1*t_stern)]))[:2, 2:] - I2).max()),
       "Eigenwerte des Polarteils bei 1,1 t*": np.round(np.linalg.eigvals(tomita1(A0, B0, P_OM @ np.linalg.inv(
           Rw([(BOOST[0], 1.1*t_stern)])))[1]), 9).tolist()}
ZAHLEN["lage"] = {"Conj_Omega-Form (20 Zufallselemente), max Abweichung": float(max(form)),
                  "|Phi - Phi_Bahn| bei t = 0,3 / 0,45": vgl, "Singulärwerte dPhi bei g = 1": np.round(sv_eins, 9).tolist(),
                  "max |dPhi| auf Lorentz": float(kern_lor), "Rang (1 und drei Zufallspunkte)": raenge,
                  "c = coth(2r) = (1 + q^2)/(2q)": c_q, "Karte: rel. Fehler tan theta / c sinh rho, |det beta - 1|":
                      [max(k[0] for k in karte), max(k[1] for k in karte)], "E-Phase e^(i phi)": phase,
                  "theta/pi bei rho = 4, 8, 16": rand_th, "kleinster Schmidtwert bei rho = 4, 8, 16": rand_sv,
                  "linke K-Äquivarianz": links, "rechts Drehung / Boost / E": [rechts_rot, rechts_boost, rechts_E],
                  "Horizont": hor}
check(f"[E] Lage: Jedes J(M_g, Omega) hält Omega fest, also liegt die Familie in Conj_Omega; gaußsch ist das genau die Menge "
      f"tau = [[0, beta], [beta^+, 0]] mit beta in U(2) (20 Zufallselemente, Abweichung {max(form):.0e}). Mit der Bahnfamilie "
      f"U(g) J U(g)^-1 stimmt sie nur auf E, den Drehungen und den Boosts überein (|Phi - Phi_Bahn| {vgl['E'][0]:.0e}, "
      f"{vgl['Drehung'][0]:.0e}, {vgl['Boost'][0]:.0e}), auf Q_1..Q_3 nicht ({vgl['Q1'][0]:.3f}, {vgl['Q2'][0]:.3f}, "
      f"{vgl['Q3'][0]:.3f}): Das Vakuum ist nur unter K = SO(3) x SO(2)_E invariant, nicht unter der AdS-Gruppe. Kern von "
      f"dPhi bei g = 1 ist Lorentz (max {kern_lor:.0e}), Rang {raenge[0]}, an drei Zufallspunkten {raenge[1:]}: eine "
      f"vierdimensionale Familie, offen in U(2). Globale Karte g = e^(phi E) e^(rho n.Q): beta = e^(i phi) exp(-i theta "
      f"n'.sigma) mit tan theta = c sinh rho, c = coth(2r) = (1 + q^2)/(2q) = {c_q:.4f} (rel. Fehler "
      f"{max(k[0] for k in karte):.0e}, det = 1 auf den Q-Strahlen, E-Phase {phase:.0e}); der AdS-Rand rho -> oo geht auf "
      f"den Äquator theta = pi/2 (theta/pi = {np.round(rand_th, 6).tolist()}, Schmidtwerte -> 1: {np.round(rand_sv, 6).tolist()}). "
      f"Links K-äquivariant ({links:.0e}), rechts nur unter Drehungen ({rechts_rot:.0e}), nicht unter Boosts und E "
      f"({rechts_boost:.3f}, {rechts_E:.3f}): Die Familie ist keine Funktion auf G/Lorentz. Verschränkungshorizont: Ein Boost "
      f"mit t* = 2 artanh q = {t_stern:.5f} macht eine Mode von M_g unverschränkt (Schmidtwert {hor['kleinster Schmidtwert bei t*']:.0e}), "
      f"Omega ist dort nicht trennend, dahinter springt J (Polarteil-Eigenwerte {hor['Eigenwerte des Polarteils bei 1,1 t*']})",
      max(form) < 1e-12 and max(vgl["E"] + vgl["Drehung"] + vgl["Boost"]) < 1e-12 and min(vgl["Q1"] + vgl["Q2"] + vgl["Q3"]) > 0.1
      and kern_lor < 1e-8 and raenge == [4, 4, 4, 4] and max(k[0] for k in karte) < 1e-10 and max(k[1] for k in karte) < 1e-12
      and phase < 1e-12 and abs(rand_th[-1] - 0.5) < 1e-6 and abs(rand_sv[-1] - 1) < 1e-6 and links < 1e-12
      and rechts_rot < 1e-12 and rechts_boost > 1e-2 and rechts_E > 1e-2 and hor["kleinster Schmidtwert bei t*"] < 1e-12
      and hor["|beta - 1| bei 0,9 t*"] < 1e-12 and hor["|beta - 1| bei 1,1 t*"] > 1)
laufzeit()

# --- 1.3 Nativer Rahmen, eigenes Modell, Berry-Krümmung, Torsion, Signatur, Schnittkrümmungen
say("  1.3 Cartan-Verbindung: nativer Rahmen (Invariante iota), eigenes Modell K_big = U(2) x U(2), Krümmung, Torsion")
Rg = Rw([(sum(c*X for c, X in zip(rng.normal(size=10), SO32))*0.7, 1.0)])
tz = Phi(Rw([(sum(c*X for c, X in zip(rng.normal(size=10), SO32))*0.5, 1.0)]))
iota_inv = abs(iota(np.linalg.inv(Rg) @ tz @ Rg) - iota(tz))
iota_bahn = max(abs(iota(Phi_bahn(Rw([(sum(c*X for c, X in zip(rng.normal(size=10), SO32))*0.6, 1.0)]))) - 4) for _ in range(5))
iota_fest = []
for rho in (0.0, 0.3, 0.8, 2.0):
    be = Phi(Rw([(QP[1], rho)]))[:2, 2:]
    th = float(np.abs(np.angle(np.linalg.eigvals(be))).max())
    iota_fest.append((rho, iota(Phi(Rw([(QP[1], rho)]))).real, 4*math.cos(2*th)))


def stab_dim(tau, basis=SO32):
    rows = [np.concatenate([kl(X, tau).real.ravel(), kl(X, tau).imag.ravel()]) for X in basis]
    sv = np.linalg.svd(np.array(rows).T, compute_uv=False)
    return int(np.sum(sv < 1e-9*sv[0]))


stab_S4 = stab_dim(S4)
stab_gen = [stab_dim(Phi(Rw([(sum(c*X for c, X in zip(rng.normal(size=10), SO32))*0.5, 1.0)]))) for _ in range(3)]
stab_kbig = stab_dim(S4, KBIG)
# Tangentialräume bei tau_0: Hauptwinkel zwischen fester Familie und Bahn (reelles Frobenius-Produkt)
Tf = np.array([np.concatenate([dPhi(np.eye(4), X).real.ravel(), dPhi(np.eye(4), X).imag.ravel()]) for X in QP]).T
Tb = np.array([np.concatenate([kl(S4, X).real.ravel(), kl(S4, X).imag.ravel()]) for X in QP]).T
cosw = np.sort(np.linalg.svd(np.linalg.qr(Tf)[0].T @ np.linalg.qr(Tb)[0], compute_uv=False))[::-1]
# Signatur der Spurform B(X, Y) = tr(X Y) auf den Erzeugern m (tau-ungerade) der beiden Familien
m_fest = [0.5*S4 @ dPhi(np.eye(4), X) for X in QP]
m_bahn = [X for X in QP]
B_fest = np.array([[np.trace(x @ y).real for y in m_fest] for x in m_fest])
B_bahn = np.array([[np.trace(x @ y).real for y in m_bahn] for x in m_bahn])
anti_ok = max(np.abs(S4 @ x + x @ S4).max() for x in m_fest)


def schnitt(x, y):
    br = kl(x, y)
    D = np.trace(x @ x).real*np.trace(y @ y).real - np.trace(x @ y).real**2
    return float(np.trace(br @ br).real/D)


K_bahn = {"Q1Q2": schnitt(m_bahn[1], m_bahn[2]), "Q1Q3": schnitt(m_bahn[1], m_bahn[3]), "EQ1": schnitt(m_bahn[0], m_bahn[1])}
K_fest = {"Q1Q2": schnitt(m_fest[1], m_fest[2]), "Q1Q3": schnitt(m_fest[1], m_fest[3]), "EQ1": schnitt(m_fest[0], m_fest[1])}
# Schleife g = e^(sX) e^(tY) (z100): Rahmen W mit tau = W^-1 S4 W; theta = dW W^-1, Gamma = P_h theta, e = P_m theta
XL, YL = QP[1] + 0.3*QP[0] + 0.2*LOR[2], QP[2] - 0.5*QP[3] + 0.4*LOR[4]


def rahmen_fest(s_, t_):
    be = Phi(Rw([(XL, s_), (YL, t_)]))[:2, 2:]
    w = sla.sqrtm(be)
    return sla.block_diag(np.linalg.inv(w), w)


def rahmen_bahn(s_, t_):
    return Rw([(XL, s_), (YL, t_)])


def d4(f, x, h):
    return (-f(x + 2*h) + 8*f(x + h) - 8*f(x - h) + f(x - 2*h))/(12*h)


def schleife(rahmen, s0=0.31, t0=-0.27, h=1e-3, H=1e-2):
    def th(s_, t_, mu):
        W = rahmen(s_, t_)
        if mu == 0:
            dW = d4(lambda x: rahmen(x, t_), s_, h)
        else:
            dW = d4(lambda x: rahmen(s_, x), t_, h)
        return dW @ np.linalg.inv(W)
    Ph = lambda X: (X + S4 @ X @ S4)/2
    Pm = lambda X: (X - S4 @ X @ S4)/2
    ths, tht = th(s0, t0, 0), th(s0, t0, 1)
    d_s_tht = d4(lambda x: th(x, t0, 1), s0, H)
    d_t_ths = d4(lambda x: th(s0, x, 0), t0, H)
    F = d_s_tht - d_t_ths - kl(ths, tht)
    Rk = Ph(d_s_tht) - Ph(d_t_ths) - kl(Ph(ths), Ph(tht))
    Tk = Pm(d_s_tht) - Pm(d_t_ths) - kl(Ph(ths), Pm(tht)) - kl(Pm(ths), Ph(tht))
    ee = kl(Pm(ths), Pm(tht))
    im_kbig = max(np.abs(ths[:2, 2:]).max(), np.abs(tht[:2, 2:]).max())
    return {"|F_MC|": float(np.abs(F).max()), "|R|": float(np.abs(Rk).max()), "|R - [e, e]|": float(np.abs(Rk - ee).max()),
            "|T|": float(np.abs(Tk).max()), "theta in k_big (Nebendiagonale)": float(im_kbig)}


sch_fest = schleife(rahmen_fest)
sch_bahn = schleife(rahmen_bahn)
ZAHLEN["cartan"] = {"iota: Invarianz unter G": float(iota_inv), "iota - 4 auf der Bahnfamilie (5 Punkte)": float(iota_bahn),
                    "iota entlang Q1 (rho, iota, 4 cos 2theta)": [[r_, round(a_, 9), round(b_, 9)] for r_, a_, b_ in iota_fest],
                    "native Stabilisatordimension bei tau_0 / generisch": [stab_S4, stab_gen],
                    "Stabilisator von tau_0 in K_big": stab_kbig,
                    "cos Hauptwinkel T_fest gegen T_Bahn bei tau_0": np.round(cosw, 12).tolist(),
                    "Spurform B fest (E, Q1..Q3)": np.round(B_fest, 9).tolist(), "Spurform B Bahn": np.round(B_bahn, 9).tolist(),
                    "Schnittkrümmung B([x,y],[x,y])/D Bahn": K_bahn, "Schnittkrümmung fest": K_fest,
                    "Schleife feste Familie (Rahmen in K_big)": sch_fest, "Schleife Bahnfamilie (nativer Rahmen)": sch_bahn}
sig_f = np.linalg.eigvalsh(B_fest)
sig_b = np.linalg.eigvalsh(B_bahn)
check(f"[E] Kein nativer Rahmen: iota(tau) = tr(tau I tau^T I^-1) mit I = I_pi - I_omega ist unter der nativen Gruppe "
      f"invariant ({iota_inv:.0e}; so(3,2) = Der(I), z98) und auf jeder Bahnfamilie konstant 4 ({iota_bahn:.0e}); entlang "
      f"der festen Familie ist iota = 4 cos 2 theta ({', '.join(f'{a_:.4f}' for _, a_, _ in iota_fest)} für rho = 0; 0,3; 0,8; "
      f"2). Die Familie schneidet also die nativen Bahnen (bei tau_0 die vierdimensionale AdS4, Stabilisator {stab_S4}; "
      f"generisch siebendimensional, Stabilisator {stab_gen[0]}) und liegt in keiner: Ein Rahmen in der nativen Gruppe "
      f"existiert nicht, F_MM im AdS-Modell ist undefiniert. Bei tau_0 teilen beide Familien nur die E-Richtung (cos der "
      f"Hauptwinkel {np.round(cosw, 9).tolist()}). Eigenes Modell: Conj_Omega ist eine Bahn von K_big = U(2) x U(2) "
      f"(Stabilisator der Vakuumebene in der ladungserhaltenden Gaußgruppe) mit Stabilisator U(2)_diag (Dimension "
      f"{stab_kbig}), das symmetrische Paar (u(2) + u(2), u(2)). Mit Rahmen in K_big (Schleife e^(sX) e^(tY)): Maurer-Cartan "
      f"{sch_fest['|F_MC|']:.0e}, Berry-Krümmung |R| = {sch_fest['|R|']:.4f} = [e, e] (Rest {sch_fest['|R - [e, e]|']:.0e}), "
      f"Torsion {sch_fest['|T|']:.0e}; die Bahnfamilie auf derselben Schleife |R| = {sch_bahn['|R|']:.4f}, Torsion "
      f"{sch_bahn['|T|']:.0e}. Spurform (z98-Normierung): Bahn Signatur {np.round(sig_b, 6).tolist()} (lorentzsch), fest "
      f"{np.round(sig_f, 6).tolist()} = -(1, c^2, c^2, c^2), definit. Schnittkrümmungen B([x,y],[x,y])/D: Bahn "
      f"{K_bahn['Q1Q2']:.6f} (Q-Ebenen), {K_bahn['EQ1']:.6f} (E-Ebenen), fest {K_fest['Q1Q2']:.6f} und {K_fest['EQ1']:.0e}: "
      f"bei definiter Metrik -B ist das +1, also S^3 statt H^3, die E-Richtung flach. Die feste Familie ist die kompakte "
      f"Duale U(2) = (S^1 x S^3)/Z_2, keine gekrümmte AdS4",
      iota_inv < 1e-12 and iota_bahn < 1e-10 and all(abs(a_ - b_) < 1e-9 for _, a_, b_ in iota_fest)
      and abs(iota_fest[-1][1] - 4) > 1 and stab_S4 == 6 and stab_gen == [3, 3, 3] and stab_kbig == 4
      and abs(cosw[0] - 1) < 1e-9 and cosw[1] < 1e-6 and anti_ok < 1e-9 and sch_fest["|F_MC|"] < 1e-6
      and sch_fest["|R - [e, e]|"] < 1e-6 and sch_fest["|T|"] < 1e-6 and sch_fest["|R|"] > 1e-2
      and sch_fest["theta in k_big (Nebendiagonale)"] < 1e-9 and sch_bahn["|F_MC|"] < 1e-6 and sch_bahn["|T|"] < 1e-6
      and np.all(sig_f < 0) and int(np.sum(sig_b > 0)) == 3 and abs(K_bahn["Q1Q2"] + 1) < 1e-9
      and abs(K_fest["Q1Q2"] + 1) < 1e-9 and abs(K_fest["EQ1"]) < 1e-12 and abs(sig_f.max() + 1) < 1e-9
      and abs(sig_f.min() + c_q**2) < 1e-6)
laufzeit()

# --- 1.4 Quadratische Antwort: BKM-Fisher der reduzierten Zustände psi|M_g (über U(g): (U(g)^* psi)|M)
say("  1.4 Quadratische Antwort (BKM-Fisher, zweite Ableitung der relativen Entropie der reduzierten Zustände)")


def bkm(psi, gens, cut=1e-13):
    C = psi.reshape(DA, DA)
    rho = C @ C.conj().T
    p, U_ = np.linalg.eigh(rho)
    keep = p > cut*p.max()
    p, U_ = p[keep], U_[:, keep]
    lp = np.log(p)
    dp = p[:, None] - p[None, :]
    gross = np.abs(dp) > 1e-6*np.maximum(p[:, None], p[None, :])
    W = np.where(gross, (lp[:, None] - lp[None, :])/np.where(gross, dp, 1), 2/(p[:, None] + p[None, :]))
    Ds = []
    for H in gens:
        dC = (-1j*(H @ psi)).reshape(DA, DA)
        drho = dC @ C.conj().T + C @ dC.conj().T
        Ds.append(U_.conj().T @ drho @ U_)
    return np.array([[float(np.sum((Di.conj()*Dj).real*W)) for Dj in Ds] for Di in Ds])


F_vak = bkm(OMEGA, HS)
ew_vak = np.linalg.eigvalsh(F_vak)
NAMEN = [f"L{k}" for k in range(6)] + ["E", "Q1", "Q2", "Q3"]
kern_K = max(np.abs(F_vak[6, :]).max(), max(np.abs(bkm(OMEGA, [fock_gen(R_)])[0, 0]) for R_ in ROT))
F_boost = [float(bkm(OMEGA, [fock_gen(B_)])[0, 0]) for B_ in BOOST]
F_Q = [float(F_vak[7 + i, 7 + i]) for i in range(3)]
ZAHLEN["bkm_vakuum"] = {"Eigenwerte (10 native Richtungen)": np.round(ew_vak, 9).tolist(), "Kern auf E und Drehungen": kern_K,
                        "Boosts": F_boost, "Q1..Q3": F_Q}
rang_bkm = int(np.sum(ew_vak > 1e-9*ew_vak.max()))
check(f"[E] Quadratische Antwort im Vakuum: Die BKM-Fisher-Metrik der reduzierten Zustände (zweite Ableitung von "
      f"S(psi|M_g || psi|M), über U(g) nach M gebracht) ist auf E und den Drehungen null ({kern_K:.0e}) und auch auf "
      f"Q_1..Q_3 ({max(F_Q):.0e}); positiv ist sie nur auf den Boosts ({np.round(F_boost, 9).tolist()}), Rang {rang_bkm}. "
      f"Die nativen Richtungen zerfallen also exakt: Drehungen halten M fest (Stab(M) = SO(3)); Boosts bewegen nur das "
      f"Modulspektrum (J fest, BKM = 2 je Einheit); E und Q bewegen nur J (die U(2)-Familie, reduzierter Zustand in erster "
      f"Ordnung fest); 3 + 4 = 7 = dim G/SO(3). Die quadratische Antwort verschwindet damit genau auf den geometrischen "
      f"Richtungen der Familie: Sie kann keinen Koeffizienten b einer Krümmung dieser Familie festlegen",
      kern_K < 1e-9 and max(F_Q) < 1e-12 and all(abs(x - 2) < 1e-6 for x in F_boost) and rang_bkm == 3)
befund("E (Vakuum): J(M_g, Omega) reproduziert AdS4 nicht und krümmt es auch nicht. Die Familie liegt in Conj_Omega, der Menge "
       "der Konjugationen, die Omega festhalten; gaußsch ist das U(2) = (S^1 x S^3)/Z_2, eine Bahn der Stabilisatorgruppe "
       "K_big = U(2) x U(2) des Zustands, nicht der nativen Gruppe. Sie stimmt mit der Bahnfamilie nur auf E, Drehungen und "
       "Boosts überein; die Q-Bewegung der Algebra dreht beta in SU(2) (tan theta = coth(2r) sinh rho, der AdS-Rand geht auf "
       "den Äquator). Die native Invariante iota = 4 cos 2 theta variiert: kein nativer Rahmen. Im eigenen Modell ist die "
       "Cartan-Verbindung Maurer-Cartan (F = 0), die Berry-Krümmung die der runden S^3, die Torsion null, die Spurform "
       "definit: die kompakte Duale der AdS-Kinematik (H^3 -> S^3), keine Dynamik. Die quadratische Antwort lebt nur auf den "
       "Boosts (J fest) und verschwindet auf der Familie", True)
laufzeit()

# =================================================================== Teil 2
say("\nTeil 2: Nicht invariante Zustände -- Kondensate, KMS-Zustand des Originals, kollektiver Zustand")

# --- 2.1 Kondensat (a1): Verschiebungslemma J(N, W psi) = W J(N, psi) W^* für jede von linearen Feldern erzeugte Algebra
say("  2.1 Impulsphasen-Kondensat (a1): Verschiebungslemma")
CC = [(d_[i] + dag(e_[i])).tocsr() for i in range(2)]               # c = d + e^+ = kappa (a + b^+)
BETA = np.array([0.03, 0.02j])


def koh(beta):
    return sum((beta[j]*dag(CC[j]) - np.conj(beta[j])*CC[j] for j in range(2)), NULL).tocsc()


L_B = koh(BETA)
PSI_B = spla.expm_multiply(L_B, OMEGA)
proben = [OMEGA, ERZ[0] @ OMEGA, ERZ[2] @ OMEGA, PSI_B]
low4 = np.where(TOT <= 4)[0]
lemma, fixB = {}, {}
for name_, (X, t) in {"Q1 (t = 0,3)": (QP[1], 0.3), "Boost (t = 0,2)": (BOOST[1], 0.2), "Zufall": (
        sum(c*Y for c, Y in zip(rng.normal(size=10), SO32))*0.2, 1.0)}.items():
    H = fock_gen(X)
    VB = V_von(PSI_B, H, t)
    V0 = V_von(OMEGA, H, t)
    dev = 0.0
    for v in proben:
        lhs = J_auf(H, VB, v, t)
        rhs = spla.expm_multiply(L_B, J_auf(H, V0, spla.expm_multiply(-L_B, v), t))
        dev = max(dev, float(np.abs((lhs - rhs)[low4]).max()))
    lemma[name_] = dev
    fixB[name_] = float(np.linalg.norm(J_auf(H, VB, PSI_B, t) - PSI_B))
# Fisher-Überschuss des Kondensats: quadratisch in beta (frei)
F_B1 = bkm(PSI_B, HS)
F_B2 = bkm(spla.expm_multiply(koh(2*BETA), OMEGA), HS)
skal_B = float(np.abs(F_B2 - F_vak).max()/np.abs(F_B1 - F_vak).max())
ZAHLEN["kondensat_a1"] = {"max |J(M_g, W Omega) - W J(M_g, Omega) W^*| (tiefer Block)": lemma, "|J psi_beta - psi_beta|": fixB,
                          "BKM-Überschuss max": float(np.abs(F_B1 - F_vak).max()), "Skalierung beta -> 2 beta": skal_B}
check(f"[E] (a1) Verschiebungslemma: Für jede Algebra N, die von linearen Feldern erzeugt wird, und jede Verschiebung W gilt "
      f"W^* N W = N, also J(N, W Omega) = W J(N, Omega) W^* und Delta ebenso. Fock (NV = 11, beta = (0,03; 0,02i), drei "
      f"bewegte Algebren M_g, vier Probevektoren): Abweichung max {max(lemma.values()):.1e} im tiefen Block, J hält psi_beta "
      f"fest ({max(fixB.values()):.0e}). Die Kondensatfamilie ist die Vakuumfamilie, mit der festen Unitären W_beta "
      f"konjugiert: dieselbe Lage (U(2)), dieselbe Krümmung, Torsion und Signatur; beta fällt aus der Geometrie heraus. Nur "
      f"die quadratische Antwort trägt beta: BKM-Überschuss quadratisch (beta -> 2 beta: Faktor {skal_B:.4f}), also frei",
      max(lemma.values()) < 1e-6 and max(fixB.values()) < 1e-6 and abs(skal_B - 4) < 2e-2)
laufzeit()

# --- 2.2 D = -4-Kondensat (a2): Fock, Tangentialabbildung bei g = 1
say("  2.2 D = -4-Kondensat (a2): Fixpunkt, Rang der Familie, lambda-Skalierung, Fisher")
H_GRAV = (dag(e_[0]) @ dag(e_[0]) @ dag(e_[1]) @ dag(e_[1]) @ OMEGA)/2                 # |0, 0; 2, 2>, D = -4
H_GRAV40 = (dag(e_[0]) @ dag(e_[0]) @ dag(e_[0]) @ dag(e_[0]) @ OMEGA)/math.sqrt(24)    # |0, 0; 4, 0>


def kond4(lam, h=H_GRAV):
    v = OMEGA + lam*h
    return v/np.linalg.norm(v)


L0c = 3
cols = np.where(TOT <= L0c)[0]
cols2 = np.where(TOT <= L0c + 2)[0]
pos2 = {j: k for k, j in enumerate(cols2)}


def J_spalten(V, cs):
    out = np.zeros((DIM, len(cs)), complex)
    for k, j in enumerate(cs):
        ia, ib = divmod(j, DA)
        out[:, k] = np.kron(V[:, ib], V[ia, :])
    return out


def dV_lyap(C, dC):
    """Ableitung des Polarteils V von C = P V: Y = dC V^+, P K + K P = -i (Y - Y^+), dV = i K V."""
    U_, sv, Wh = np.linalg.svd(C)
    V = U_ @ Wh
    Y = dC @ V.conj().T
    Rhs = U_.conj().T @ (-1j*(Y - Y.conj().T)) @ U_
    K = U_ @ (Rhs/(sv[:, None] + sv[None, :])) @ U_.conj().T
    return 1j*K @ V, V


def tangenten(psi):
    """Tangenten dL(X) der Familie L_g = U (V_g (x) V_g^T) SWAP U^T (J = L K) bei g = 1, Spalten N_tot <= 3."""
    C = psi.reshape(DA, DA)
    out = []
    for H in HS:
        dC = (-1j*(H @ psi)).reshape(DA, DA)
        dV, V = dV_lyap(C, dC)
        U2 = J_spalten(V, cols2)
        U1 = U2[:, [pos2[j] for j in cols]]
        bahn = 1j*(H @ U1) + 1j*(U2 @ H.conj()[cols2][:, cols].toarray())
        zust = np.zeros((DIM, len(cols)), complex)
        for k, j in enumerate(cols):
            ia, ib = divmod(j, DA)
            zust[:, k] = np.kron(dV[:, ib], V[ia, :]) + np.kron(V[:, ib], dV[ia, :])
        out.append(bahn + zust)
    return out


def rang_info(T):
    M = np.array([np.concatenate([t.real.ravel(), t.imag.ravel()]) for t in T]).T
    U_, sv, vh = np.linalg.svd(M, full_matrices=False)
    return sv, vh


T_vak = tangenten(OMEGA)
sv_vak, _ = rang_info(T_vak)
BR = np.array([np.concatenate([M.real.ravel(), M.imag.ravel()]) for M in SO32]).T
koord = lambda M: np.linalg.lstsq(BR, np.concatenate([M.real.ravel(), M.imag.ravel()]), rcond=None)[0]
ROT_K = np.array([koord(R_) for R_ in ROT])
BOOST_K = np.array([koord(B_) for B_ in BOOST])


def boost_daten(psi, T):
    tan_ = [float(np.linalg.norm(sum(c*T[i] for i, c in enumerate(bk)))) for bk in BOOST_K]
    fis_ = [float(bkm(psi, [fock_gen(B_)])[0, 0]) for B_ in BOOST]
    return np.array(tan_), np.array(fis_)


F_boost0 = np.array(F_boost)
fixL, sv_l, kernw, bt, bf = {}, {}, {}, {}, {}
for name_, h_ in (("|0,0;2,2>", H_GRAV), ("|0,0;4,0>", H_GRAV40)):
    for lam in (1e-4, 2e-4):
        psi = kond4(lam, h_)
        T_ = tangenten(psi)
        sv_, vh_ = rang_info(T_)
        sv_l[(name_, lam)] = sv_
        ker = vh_[-3:]
        kernw[(name_, lam)] = float(np.linalg.svd(np.linalg.qr(ker.T)[0].T @ np.linalg.qr(ROT_K.T)[0], compute_uv=False).min())
        V_ = polar(psi.reshape(DA, DA))[0]
        fixL[(name_, lam)] = float(np.linalg.norm(J_auf(None, V_, psi) - psi))
        bt[(name_, lam)], fis_ = boost_daten(psi, T_)
        bf[(name_, lam)] = fis_ - F_boost0
lin_t = [float(x) for x in bt[("|0,0;2,2>", 2e-4)]/bt[("|0,0;2,2>", 1e-4)]]
lin_f = [float(x) for x in bf[("|0,0;2,2>", 2e-4)]/bf[("|0,0;2,2>", 1e-4)]]
r_l = {k: float(bf[k].sum()/np.sum(bt[k]**2)) for k in bt}
r_einzel = {k: (bf[k]/bt[k]**2).tolist() for k in bt}
rang_l = {k: int(np.sum(v > 1e-9*v[0])) for k, v in sv_l.items()}
rang_vak = int(np.sum(sv_vak > 1e-6*sv_vak[0]))
spur = [float(np.sum(bt[(n_, 1e-4)]**2)) for n_ in ("|0,0;2,2>", "|0,0;4,0>")]
einzel_streu = max(abs(a_ - b_) for a_, b_ in zip(r_einzel[("|0,0;2,2>", 1e-4)], r_einzel[("|0,0;4,0>", 1e-4)]))
ZAHLEN["kondensat_a2"] = {"|J psi - psi|": {str(k): v for k, v in fixL.items()},
                          "Singulärwerte Vakuum": np.round(sv_vak, 12).tolist(),
                          "Singulärwerte je (Komponente, lambda)": {str(k): np.round(v, 12).tolist() for k, v in sv_l.items()},
                          "Rang je (Komponente, lambda)": {str(k): v for k, v in rang_l.items()},
                          "cos kleinster Hauptwinkel Kern gegen Drehungen": {str(k): v for k, v in kernw.items()},
                          "Boost-Tangenten (Frobenius, Spalten N_tot <= 3)": {str(k): v.tolist() for k, v in bt.items()},
                          "Boost-BKM-Überschuss": {str(k): v.tolist() for k, v in bf.items()},
                          "Verhältnis 2e-4/1e-4 Tangente, BKM": [lin_t, lin_f],
                          "Spur BKM-Überschuss / Spur Tangente^2": {str(k): v for k, v in r_l.items()},
                          "Einzelverhältnisse je Boost": {str(k): v for k, v in r_einzel.items()}}
check(f"[E] (a2) Das D = -4-Kondensat psi = Omega + lambda h (h = |0, 0; 2, 2> oder |0, 0; 4, 0>, lambda = 1e-4, 2e-4 "
      f"<< q^4) ist Fixpunkt seiner Familie (|J psi - psi| {max(fixL.values()):.0e}). Tangentialabbildung bei g = 1 (Fock, "
      f"Spalten N_tot <= 3): Vakuum Rang {rang_vak} (Lorentz null), Kondensat Rang {max(rang_l.values())} (beide Komponenten, beide lambda); der Kern "
      f"ist genau die Drehalgebra Stab(M) = so(3) (cos der Hauptwinkel >= {min(kernw.values()):.9f}), die jede Familie "
      f"J(M_g, psi) trivial festhält, weil sie M nicht bewegt. Die drei Boostrichtungen bekommen Tangenten linear in lambda "
      f"(Verhältnis {np.round(lin_t, 6).tolist()}), der BKM-Überschuss auf den Boosts ist quadratisch "
      f"({np.round(lin_f, 4).tolist()}): Die Familie entfaltet sich von 4 auf 7 Dimensionen (= dim G/SO(3), die ganze "
      f"Algebrabahn) und hat kein vierdimensionales Modell; ihre Größe ist die freie Amplitude. lambda kürzt sich im "
      f"Verhältnis Spur BKM-Überschuss / Spur Tangente^2 = {r_l[('|0,0;2,2>', 1e-4)]:.6f} (gleich für |0,0;4,0>: "
      f"{r_l[('|0,0;4,0>', 1e-4)]:.6f}, Schur auf dem Spin-2-Multiplett), die Einzelrichtungen hängen von der Komponente ab "
      f"(Streuung {einzel_streu:.4f}); die Tangentennorm ist aber eine Abschneidewahl, kein natives Maß",
      max(fixL.values()) < 1e-8 and rang_vak == 4 and set(rang_l.values()) == {7} and min(kernw.values()) > 1 - 1e-6
      and all(abs(x - 2) < 1e-3 for x in lin_t) and all(abs(x - 4) < 2e-2 for x in lin_f)
      and abs(r_l[("|0,0;2,2>", 1e-4)] - r_l[("|0,0;4,0>", 1e-4)]) < 1e-3*r_l[("|0,0;2,2>", 1e-4)] and einzel_streu > 1e-4)
laufzeit()

# --- 2.3 KMS-Zustand des Originals bei beta = 4 pi: natives q = e^(-2 pi)
say("  2.3 KMS-Zustand des Originals bei beta = 4 pi (q^2 = e^(-4 pi))")
qn = Q_NATIV
c_n = (1 + qn*qn)/(2*qn)
A0n, B0n = ebenen(qn)
karte_n = []
for rho in (1e-3, 0.5, 3.0):
    be = Phi(Rw([(QP[3], rho)]), qn)[:2, 2:]
    th = float(np.abs(np.angle(np.linalg.eigvals(be))).max())
    karte_n.append(abs(math.tan(th) - c_n*math.sinh(rho))/(c_n*math.sinh(rho)))
form_n = []
for _ in range(10):
    t_ = Phi(Rw([(sum(c*X for c, X in zip(rng.normal(size=10), SO32))*0.5, 1.0)]), qn)
    form_n.append(max(np.abs(t_[:2, :2]).max(), np.abs(t_[2:, 2:]).max(), np.abs(t_ @ t_ - np.eye(4)).max()))
sv_n = rang_sv(np.eye(4), qn)
rang_n = int(np.sum(sv_n > 1e-6*sv_n[0]))
t_stern_n = 2*math.atanh(qn)
hor_n = float(tomita1(A0n, B0n, P_OM @ np.linalg.inv(Rw([(BOOST[0], t_stern_n)])))[2].min())
# Fockkontrolle bei nativem q (Abschneidung unkritisch)
d_n, e_n, fock_gen_n, OMEGA_n = fock_daten(qn)
Hn = fock_gen_n(QP[1])
VFn = V_von(OMEGA_n, Hn, 0.4)
V1n = tomita1(A0n, B0n, P_OM @ np.linalg.inv(Rw([(QP[1], 0.4)])))[1]
pol_n = float(np.abs((VFn - VFn[0, 0]*gamma1(V1n))[np.ix_(lowa, lowa)]).max())
HS_n = [fock_gen_n(M) for M in SO32]
F_n = bkm(OMEGA_n, HS_n, cut=1e-30)
ew_n = np.linalg.eigvalsh(F_n)
ZAHLEN["kms_4pi"] = {"q = e^(-2 pi)": qn, "c = (1 + q^2)/(2q) = cosh(2 pi)": [c_n, math.cosh(2*math.pi)],
                     "Karte rel. Fehler": max(karte_n), "Conj-Form": max(form_n), "Rang dPhi": rang_n,
                     "Horizont t* = 2 artanh q": t_stern_n, "Schmidtwert bei t*": hor_n,
                     "Fock |V - Gamma(V_1)|": pol_n, "BKM Eigenwerte": np.round(ew_n, 9).tolist(),
                     "BKM Boost, Q": [float(bkm(OMEGA_n, [fock_gen_n(BOOST[0])], cut=1e-30)[0, 0]), float(F_n[7, 7])]}
check(f"[E] KMS-Zustand des Originals bei beta = 4 pi (q = e^(-2 pi) = {qn:.6e}): Er ist unter K invariant, also gilt die "
      f"Vakuumstruktur unverändert (Conj-Form {max(form_n):.0e}, Rang {rang_n}, Fockkontrolle {pol_n:.0e}). Die Karte ist "
      f"tan theta = c sinh rho mit c = (1 + q^2)/(2q) = cosh(2 pi) = {c_n:.4f} (rel. Fehler {max(karte_n):.0e}): Ein Q-Schritt "
      f"dreht beta um das {c_n:.1f}-fache des E-Schritts, der Äquator (AdS-Rand) ist schon bei rho ~ 1/c erreicht; der "
      f"Verschränkungshorizont liegt bei t* = 2 artanh q = {t_stern_n:.3e} (Schmidtwert {hor_n:.0e}). BKM-Fisher wie bei "
      f"q = 0,25 nur auf den Boosts, q-unabhängig 2 (Eigenwerte {np.round(ew_n[ew_n > 1e-9], 6).tolist()}): fest, ohne "
      f"freien Parameter, aber null auf der Familie, und im eigenen Modell bleibt F = 0",
      max(form_n) < 1e-10 and rang_n == 4 and max(karte_n) < 1e-8 and abs(c_n - math.cosh(2*math.pi)) < 1e-9
      and hor_n < 1e-10 and pol_n < 1e-8 and int(np.sum(ew_n > 1e-9*ew_n.max())) == 3
      and np.allclose(ew_n[ew_n > 1e-9*ew_n.max()], 2, atol=1e-6))
laufzeit()

# --- 2.4 Kollektiver Zustand von z101: Saatraum, Qubit-Lesart (Original und Partner je ein Qubit je Saat), Paargruppe
say("  2.4 Kollektiver Zustand N = 6 (z101, z94): Saatraum C^64 (x) C^64, Paargruppe je Saat")
NS = 6
NC = 1 << NS
POP = np.array([bin(m).count("1") for m in range(NC)])
YN = math.exp(-math.pi)
z94 = json.loads((HERE / "z94_ergebnis.json").read_text())                # nur gelesen
f94 = z94["zahlen"]["aussenbehauptung_gleiche_saaten"]["f(n) für N = 6"]
F94 = np.zeros(NC, dtype=np.int64)
for s_, v_ in f94.items():
    F94[sum(int(ch) << k for k, ch in enumerate(s_))] = v_
E94 = POP + 4*F94
p94 = YN**E94.astype(float)
p94 /= p94.sum()
pP = YN**POP.astype(float)
pP /= pP.sum()
rand94 = max(abs(p94[[m for m in range(NC) if (m >> k) & 1]].sum() - YN/(1 + YN)) for k in range(NS))
sp1 = np.zeros((NC, NC))
spk = []
for k in range(NS):
    S_ = np.zeros((NC, NC))
    for m in range(NC):
        if not (m >> k) & 1:
            S_[m | (1 << k), m] = 1
    spk.append(sps.csr_matrix(S_))
Id64 = sps.identity(NC, format="csr")
Tx = [(sps.kron(spk[k], spk[k]) + sps.kron(spk[k], spk[k]).T).tocsr() for k in range(NS)]          # Q-artig
Ty = [(1j*(sps.kron(spk[k], spk[k]) - sps.kron(spk[k], spk[k]).T)).tocsr() for k in range(NS)]     # Lorentz-artig
nk = [(spk[k] @ spk[k].T).tocsr() for k in range(NS)]
Tz = [(0.5*(sps.kron(nk[k], Id64) + sps.kron(Id64, nk[k]) - sps.identity(NC*NC))).tocsr() for k in range(NS)]
_perm = np.array([(j % NC)*NC + j//NC for j in range(NC*NC)])
SWAP = sps.csr_matrix((np.ones(NC*NC), (_perm, np.arange(NC*NC))), shape=(NC*NC, NC*NC))


def lokal_anteil(O, k):
    """Bester Frobenius-Näherer der Form o_k (x) 1 (Saat k): Teilspur über die übrigen Saaten."""
    T_ = O.reshape([2]*NS*2)
    ax_o = NS - 1 - k
    ax_i = 2*NS - 1 - k
    rest_o = [a for a in range(NS) if a != ax_o]
    rest_i = [NS + a for a in range(NS) if a != ax_o]
    Tp = np.transpose(T_, [ax_o, ax_i] + rest_o + rest_i).reshape(2, 2, NC//2, NC//2)
    o = np.einsum("abii->ab", Tp)/(NC//2)
    L_ = np.einsum("ab,ij->abij", o, np.eye(NC//2))
    perm = [ax_o, ax_i] + rest_o + rest_i
    inv = np.argsort(perm)
    return np.transpose(L_.reshape([2, 2] + [2]*(NS - 1) + [2]*(NS - 1)), inv).reshape(NC, NC)


def saat_daten(p):
    xi = np.diag(np.sqrt(p)).reshape(-1).astype(complex)
    C = xi.reshape(NC, NC)
    dVs, nl = {}, {}
    for name_, gens in (("Tx", Tx), ("Ty", Ty), ("Tz", Tz)):
        for k in range(NS):
            dC = (-1j*(gens[k] @ xi)).reshape(NC, NC)
            dV, V = dV_lyap(C, dC)
            dVs[(name_, k)] = dV
            if name_ != "Ty":
                n_ = np.linalg.norm(dV)
                nl[(name_, k)] = float(np.linalg.norm(dV - lokal_anteil(dV, k))/n_) if n_ > 0 else 0.0
    # Tangenten der Familie: L = (2 i H + dV (x) 1 + 1 (x) dV^T) SWAP (V = 1); Gram-Matrix
    tang = []
    for name_, gens in (("Tx", Tx), ("Ty", Ty), ("Tz", Tz)):
        for k in range(NS):
            H = gens[k]
            orbit = 1j*(H + SWAP @ H.T @ SWAP)                                  # dL SWAP = i (H + SWAP H^T SWAP) + ...
            dV = dVs[(name_, k)]
            zust = sps.kron(sps.csr_matrix(dV), Id64) + sps.kron(Id64, sps.csr_matrix(dV.T))
            tang.append((orbit + zust).tocsr())
    G = np.array([[np.sum((A_.conj().multiply(B_)).real) for B_ in tang] for A_ in tang])
    ew = np.linalg.eigvalsh(G)
    # BKM der reduzierten Zustände für die Q-artigen Einzelsaat-Erzeuger
    rho = C @ C.conj().T
    pv = np.real(np.diag(rho))
    lp = np.log(pv)
    dp = pv[:, None] - pv[None, :]
    gross = np.abs(dp) > 1e-6*np.maximum(pv[:, None], pv[None, :])
    W = np.where(gross, (lp[:, None] - lp[None, :])/np.where(gross, dp, 1), 2/(pv[:, None] + pv[None, :]))
    Ds, Dx = [], []
    for k in range(NS):
        dC = (-1j*(Ty[k] @ xi)).reshape(NC, NC)
        Ds.append(dC @ C.conj().T + C @ dC.conj().T)
        dCx = (-1j*(Tx[k] @ xi)).reshape(NC, NC)
        Dx.append(dCx @ C.conj().T + C @ dCx.conj().T)
    Fx = np.array([[float(np.sum((Di.conj()*Dj).real*W)) for Dj in Ds] for Di in Ds])
    fx_null = max(float(np.abs(D_).max()) for D_ in Dx)
    ty_null = max(np.abs(dVs[("Ty", k)]).max() for k in range(NS))
    return {"nichtlokal": nl, "Gram-Eigenwerte": ew, "Fx": Fx, "Ty null": float(ty_null), "Tx: |d rho|": fx_null}


sd_P = saat_daten(pP)
sd_94 = saat_daten(p94)
nlP = max(sd_P["nichtlokal"].values())
nl94 = [sd_94["nichtlokal"][("Tx", k)] for k in range(NS)]
rangP = int(np.sum(sd_P["Gram-Eigenwerte"] > 1e-9*sd_P["Gram-Eigenwerte"].max()))
rang94 = int(np.sum(sd_94["Gram-Eigenwerte"] > 1e-9*sd_94["Gram-Eigenwerte"].max()))
offP = float(np.abs(sd_P["Fx"] - np.diag(np.diag(sd_P["Fx"]))).max())
off94 = float(np.abs(sd_94["Fx"] - np.diag(np.diag(sd_94["Fx"]))).max())
diagP, diag94 = np.diag(sd_P["Fx"]), np.diag(sd_94["Fx"])
ZAHLEN["kollektiv_z101"] = {"Randfehler z94 (Saatbesetzung y/(1+y))": float(rand94),
                            "Lorentz-artige Paarrichtungen: max |dV|": [sd_P["Ty null"], sd_94["Ty null"]],
                            "Nichtlokalität dV (Tx, Tz) Produkt max": nlP, "Nichtlokalität dV (Tx) je Saat, z94": nl94,
                            "Rang der Familie (18 Paarerzeuger) Produkt / z94": [rangP, rang94],
                            "BKM Ty diagonal Produkt": np.round(diagP, 9).tolist(), "BKM Ty diagonal z94": np.round(diag94, 9).tolist(),
                            "BKM Ty Nebendiagonale max Produkt / z94": [offP, off94],
                            "Tx (Q-artig): max |d rho| Produkt / z94": [sd_P["Tx: |d rho|"], sd_94["Tx: |d rho|"]]}
check(f"[E] Kollektiver Zustand (z101) in der Qubit-Lesart (Original und Partner je ein Qubit je Saat, Purifikation "
      f"sum sqrt(p_b) |b, b>, z94-Verteilung mit Randfehler {rand94:.0e}): Die Lorentz-artigen Paarrichtungen i(s^+ s^+ - h.c.) "
      f"vertauschen mit J_0 und halten C hermitesch positiv, also bleibt J fest (max |dV| {sd_P['Ty null']:.0e} bzw. "
      f"{sd_94['Ty null']:.0e}), für jeden Zustand in Standardform. Die Q-artigen Paarrichtungen bewegen J; beim "
      f"Produktzustand ist die Antwort saatlokal (Nichtlokalität {nlP:.0e}), beim kollektiven Zustand nicht (Anteil außerhalb "
      f"der Saat {np.round(nl94, 4).tolist()}, beherrscht von den Gewichtsverhältnissen y^(E(b+k) - E(b))): eine kollektive "
      f"Relativantwort mit festem, quantisiertem Koeffizienten. Rang der J-Familie unter den 18 Paarerzeugern in beiden "
      f"Fällen {rangP} bzw. {rang94}: eine Phase je Saat, Ty null. Wie im Vakuum zerfallen die Richtungen: Q-artige bewegen "
      f"J und lassen den reduzierten Zustand in erster Ordnung fest (|d rho| {max(sd_P['Tx: |d rho|'], sd_94['Tx: |d rho|']):.0e}), "
      f"Lorentz-artige bewegen nur das Spektrum. BKM auf den Lorentz-artigen Richtungen: Produkt saatdiagonal "
      f"(Nebendiagonale {offP:.0e}, je Saat {diagP[0]:.6f}), kollektiv mit Saatkopplung (Nebendiagonale bis {off94:.6f}). "
      f"Ein natives AdS-Modell gibt es auf dem Saatraum nicht (der hartkernige Raum ist nicht so(3,2)-invariant, z100)",
      rand94 < 1e-12 and sd_P["Ty null"] < 1e-12 and sd_94["Ty null"] < 1e-12 and nlP < 1e-10 and min(nl94) > 1e-3
      and rangP == rang94 == 6 and offP < 1e-10 and off94 > 1e-6
      and max(sd_P["Tx: |d rho|"], sd_94["Tx: |d rho|"]) < 1e-12)
befund("E (nicht invariante Zustände): Keiner gibt eine gekrümmte Familie mit festem Koeffizienten. (a1) ist die "
       "Vakuumfamilie, mit W_beta konjugiert (beta nur in der quadratischen Antwort, frei). (a2) bricht die "
       "Boost-Nullrichtungen linear in der freien Amplitude lambda; die Familie entfaltet sich auf die ganze Algebrabahn "
       "G/SO(3) (Rang 7, null bleiben nur die Drehungen, die M festhalten) und hat kein vierdimensionales Modell. Der "
       "KMS-Zustand bei beta = 4 pi ist K-invariant und gibt die Vakuumstruktur mit c = cosh(2 pi), fest, aber flach im "
       "eigenen Modell. Der kollektive Zustand gibt eine feste, nicht saatlokale Relativantwort und Saatkopplungen in der "
       "quadratischen Antwort, aber kein natives AdS-Modell und eine Auswahlprämisse (z101). In allen Fällen verschwindet die "
       "quadratische Antwort auf den Richtungen, die J bewegen", True)
laufzeit()

# =================================================================== Teil 3
say("\nTeil 3: Grenzfall, Bilanz und Urteil")
# Konformitätsdefekt der Karte AdS4 -> U(2): kappa(rho) = (d theta/d rho) sinh rho / sin theta (= 1 genau für c = 1)


def theta_von(rho, q_):
    be = Phi(Rw([(QP[3], rho)]), q_)[:2, 2:]
    return float(np.abs(np.angle(np.linalg.eigvals(be))).max())


konf = {}
for q_ in (Q_NATIV, 0.25, 0.9, 0.999):
    ks = []
    for rho in (0.3, 1.0, 2.0):
        dth = (theta_von(rho + 1e-5, q_) - theta_von(rho - 1e-5, q_))/2e-5
        ks.append(dth*math.sinh(rho)/math.sin(theta_von(rho, q_)))
    c_ = (1 + q_*q_)/(2*q_)
    soll = [math.cosh(r_)/math.sqrt(1 + c_*c_*math.sinh(r_)**2) for r_ in (0.3, 1.0, 2.0)]
    konf[q_] = {"c": c_, "kappa(rho = 0,3; 1; 2)": ks, "Formel cosh rho/sqrt(1 + c^2 sinh^2 rho)": soll,
                "Defekt max |kappa - 1|": max(abs(k_ - 1) for k_ in ks)}
ZAHLEN["grenzfall_konform"] = {str(k): v for k, v in konf.items()}
check(f"[E/L] Grenzfall maximaler Verschränkung: Die Karte tan theta = c sinh rho hat das Skalenverhältnis kappa = "
      f"(d theta/d rho) sinh rho/sin theta = cosh rho/sqrt(1 + c^2 sinh^2 rho) (aus der Familie numerisch, Formel bis "
      f"{max(max(abs(a_ - b_) for a_, b_ in zip(v['kappa(rho = 0,3; 1; 2)'], v['Formel cosh rho/sqrt(1 + c^2 sinh^2 rho)'])) for v in konf.values()):.0e}); "
      f"konform (kappa = 1) ist sie genau für c = 1, also q -> 1: Defekt {konf[0.999]['Defekt max |kappa - 1|']:.1e} bei "
      f"q = 0,999, {konf[0.9]['Defekt max |kappa - 1|']:.3f} bei 0,9, {konf[0.25]['Defekt max |kappa - 1|']:.3f} bei 0,25, "
      f"{konf[Q_NATIV]['Defekt max |kappa - 1|']:.4f} bei nativem q (kappa ~ 1/(c tanh rho)). Im Grenzfall ist die feste Familie die konforme "
      f"Kompaktifizierung von AdS4 in den halben Einsteinzylinder (Halbkugelmodell von H^3, E-Zeit als Phase von U(2)); bei "
      f"nativem q ist sie stark gestaucht. Lesart: Der Abstand vom Grenzfall misst die endliche Verschränkung des Zustands, "
      f"keine Krümmung",
      konf[0.999]["Defekt max |kappa - 1|"] < 1e-5 and konf[0.25]["Defekt max |kappa - 1|"] > 0.1
      and konf[Q_NATIV]["Defekt max |kappa - 1|"] > 0.9
      and all(max(abs(a_ - b_) for a_, b_ in zip(v["kappa(rho = 0,3; 1; 2)"], v["Formel cosh rho/sqrt(1 + c^2 sinh^2 rho)"])) < 1e-6
              for v in konf.values()))
STATUS = {
    "Lage der Familie": ("E", "J(M_g, psi) hält psi fest: Familie in Conj_psi; Vakuum gaußsch U(2) = (S^1 x S^3)/Z_2"),
    "Bahnfamilie": ("E, negativ", "stimmt nur auf E, Drehungen und Boosts überein; Q dreht beta in SU(2), tan theta = "
                                  "coth(2r) sinh rho"),
    "nativer Rahmen": ("E, negativ", "iota = 4 cos 2 theta variiert: Familie schneidet die nativen Bahnen, F_MM im AdS-Modell "
                                     "undefiniert"),
    "eigenes Modell": ("E", "K_big = U(2) x U(2)-Bahn, symmetrisch: F = 0, Berry = [e, e] (runde S^3), Torsion 0, definit"),
    "Verschränkungshorizont": ("E", "Boost 2 artanh q: eine Mode von M_g unverschränkt, Omega nicht trennend, J springt"),
    "Grenzfall q -> 1": ("E/L", "Karte konform: konforme Kompaktifizierung von AdS4 im halben Einsteinzylinder U(2)"),
    "quadratische Antwort": ("E", "BKM nur auf den Boosts (J fest), dort 2 und q-unabhängig; null auf E und Q, also auf der "
                                  "Familie; Überschuss beta^2 (a1), lambda^2 (a2): normiert keine Krümmung"),
    "(a1)": ("E", "Verschiebungslemma: Vakuumfamilie konjugiert"),
    "(a2)": ("E, negativ", "Boostrichtungen linear in lambda, Rang 7 = dim G/SO(3), kein 4D-Modell, Größe frei"),
    "KMS beta = 4 pi": ("E", "Vakuumstruktur mit c = cosh(2 pi), flach"),
    "kollektiv (z101)": ("E/L", "feste nicht saatlokale Relativantwort; kein natives AdS-Modell, Auswahl offen"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
urteil = ("T6 ist beantwortet, als gescoptes Negativ mit einer neuen exakten Kinematik. (1) Vakuum: J(M_g, Omega) = U(g) "
          "J(M, U(g)^* Omega) U(g)^* hält Omega fest; die Familie liegt deshalb in der Menge Conj_Omega der Konjugationen mit "
          "J Omega = Omega, gaußsch tau = [[0, beta], [beta^+, 0]], beta in U(2). Sie reproduziert AdS4 nicht: Mit der "
          "Bahnfamilie U(g) J U(g)^-1 stimmt sie nur auf E, den Drehungen und den Boosts überein, weil das Vakuum nur unter "
          "K = SO(3) x SO(2)_E invariant ist. Die Q-Bewegung der Algebra dreht beta in SU(2): In globalen Koordinaten "
          "g = e^(phi E) e^(rho n.Q) ist beta = e^(i phi) exp(-i theta n'.sigma) mit tan theta = coth(2r) sinh rho "
          "(q = tanh r); der AdS-Rand geht auf den Äquator theta = pi/2, wo die Schmidtwerte 1 werden. Ein Boost mit "
          "Schnelligkeit 2 artanh q macht eine Mode der Algebra unverschränkt (Verschränkungshorizont, J springt). Die native "
          "Invariante iota = tr(tau I tau^T I^-1), I = I_pi - I_omega, ist auf jeder Bahnfamilie 4 und auf der festen "
          "Familie 4 cos 2 theta: Die Familie schneidet die nativen Bahnen, ein Rahmen in der nativen Gruppe existiert nicht, "
          "F_MM im AdS-Modell ist nicht definiert. Ihr eigenes Modell ist das symmetrische Paar (u(2) + u(2), u(2)): "
          "Conj_Omega ist eine Bahn der Stabilisatorgruppe K_big = U(2) x U(2) des Zustands, die Cartan-Verbindung mit Rahmen "
          "dort ist Maurer-Cartan (F = 0), die Berry-Krümmung ist [e, e] mit Schnittkrümmung +1 auf den Q-Ebenen (runde S^3) "
          "und 0 mit E, die Torsion verschwindet, die Spurform ist definit. Das ist die kompakte Duale der AdS-Kinematik "
          "(H^3 -> S^3 bei derselben Uhr E), keine gekrümmte AdS4. Das z100-Lemma verallgemeinert sich: Auch die Familie bei "
          "festem Zustand ist eine Bahn, nur der Stabilisatorgruppe des Zustands statt der Bewegungsgruppe. (2) Nicht "
          "invariante Zustände: Das Kondensat (a1) gibt exakt die Vakuumfamilie, mit der festen Verschiebung konjugiert "
          "(Verschiebungslemma); beta steckt nur in der quadratischen Antwort. Das D = -4-Kondensat (a2) bricht die "
          "Boost-Nullrichtungen linear in der freien Amplitude: Die Familie entfaltet sich auf die ganze Algebrabahn G/SO(3) "
          "(Rang 7; null bleiben nur die Drehungen, die M festhalten) und hat kein vierdimensionales Modell; der "
          "BKM-Überschuss auf den Boosts ist quadratisch in lambda, lambda kürzt sich im Spurverhältnis, aber die "
          "Tangentennorm ist kein natives Maß. Der KMS-Zustand des Originals bei beta = 4 pi ist K-invariant und gibt die "
          "Vakuumstruktur mit c = cosh(2 pi): fest, aber flach. Der kollektive Zustand (z101, Qubit-Lesart) hat "
          "Lorentz-artige Nullrichtungen wie jeder Zustand in Standardform, eine nicht saatlokale, feste Relativantwort und "
          "Saatkopplungen in der quadratischen Antwort, aber kein natives AdS-Modell und eine Auswahlprämisse. b bleibt "
          "unbestimmt: Die quadratische Antwort (BKM der reduzierten Zustände) lebt nur auf den Boosts, wo J fest ist "
          "(Wert 2, q-unabhängig), und verschwindet exakt auf E und Q, also auf der Familie selbst; wo ein Koeffizient fest "
          "ist (Vakuum, KMS, kollektiv), ist die Krümmung im eigenen Modell null, und wo die Familie aus dem Modell fällt "
          "(a2), trägt sie die freie Amplitude. (3) Urteil: Bewegte Algebra bei festem Zustand gibt keine native "
          "Quelle gravitativer Dynamik, weder unbedingt noch bedingt; gescoptes Negativ für die gaußsche Klasse, das "
          "D = -4-Kondensat und den kollektiven Zustand in der Qubit-Lesart. Folgen: K13 bleibt offen; der Weg 'Krümmung aus "
          "relativer Lage' ist für homogene Zustände geschlossen, weil die Familie dann eine Bahn der Zustandssymmetrie ist. "
          "Für die Gravitationsachse der Szenario-Suche heißt das: Der Wert 'Krümmung aus bewegter Algebra' ist kein "
          "gestützter Gravitationswert; Krümmung verlangt eine Familie, die weder eine Bahn der Bewegungsgruppe noch eine "
          "der Zustandssymmetrie ist, also ein ortsabhängiges Zustandsprofil oder eine Dynamik (z101). Neu und nativ ist die "
          "Kinematik: Die feste Familie realisiert die kompakte Duale U(2) mit der Karte tan theta = coth(2r) sinh rho und "
          "dem Verschränkungshorizont; im Grenzfall q -> 1 ist die Karte konform, die konforme Kompaktifizierung von AdS4 im "
          "halben Einsteinzylinder, bei nativem q (c = cosh 2 pi) ist sie stark gestaucht. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
say(f"\n   (Laufzeit {time.time() - T0:.1f} s)")
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z115_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
