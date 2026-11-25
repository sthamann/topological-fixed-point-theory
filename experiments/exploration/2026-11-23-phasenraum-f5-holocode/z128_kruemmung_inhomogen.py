"""Z128: Krümmung aus inhomogenem Zustand? Native ortsabhängige Profile, die Familie J(M_g, psi) und die Kappe (Folge T6).

Leitlinie 'Rekonstruktion statt Auswahl'. Firewalls: GRAV.NONCIRCULAR.01 (der Modulfluss bleibt Uhr; Tomita-Familien sind
Kinematik, keine Zeit; die Kappenzeit wird nur unter der benannten Prämisse P_cap benutzt), GRAV.PARENT.GAP.NO_GO.01,
DIMENSION.UPLIFT.FIREWALL.01 (Einteilchen- und Fockkinematik einer bzw. mehrerer Nähte auf Abschneidungen; keine
wechselwirkende 3+1-Schließung). Gelöste Blöcke (alpha, Flavour, Massen, 16, E8) bleiben unberührt; nichts wird befördert.
Bestand: z115 (J(M_g, Omega) liegt in Conj_Omega = U(2), Bahn der Zustandsstabilisatorgruppe, flach; Karte
tan theta = coth(2r) sinh rho; iota = tr(tau I tau^T I^-1); BKM nur auf den Boosts), z100 (Maurer-Cartan-Lemma), z98
(MacDowell-Mansouri-Algebra, AdS-Bahn flach, so(4,1) = Der(I_pi + I_omega) für sigma = -1), z104 (Kappe O_1, K_cap = 2 pi Q^-,
Lambda-Vorzeichen = Klasse der Zeit, Marken z^4 = 1 = Spitzen der tau-Eigenkappen, Original unter dem Kappenfluss nicht
invariant ~ e^(-8 pi)), z72 (Wiesbrock: K_cap + K_Uhr/2 = 2 pi P_0), z83 (Markenidentifikation). Statuskürzel: E exakt,
C bedingt, L Lesart, O offen.

Frage (T6, z115): Gibt ein natives ortsabhängiges Zustandsprofil (ein 'q-Gradient') eine Familie J(M_g, psi), die weder
Bahn der Bewegungsgruppe noch der Stabilisatorgruppe des Zustands ist und Krümmung mit festem Koeffizienten b trägt? Und
relativ zur Kappe mit Modulzeit P_cap (hyperbolische Klasse, dS)?

Teil 1  Native inhomogene Zustände: Zustandsscheibe der SO(3)-invarianten Gaußzustände (radiale sl(2,R) als
        Möbiusgruppe, Marken als Randfixpunkte, Schmidtwert = pseudohyperbolischer Abstand); native Profile: lokales
        Hislop-Longo-Temperaturprofil der Kappe (K1a), Wiesbrock-Profil q(a) (K1b), Markenvakua (K2), kappenverschobener
        Zustand (frei); was nativ ist und was Wahl.
Teil 2  Geometrie der Familien: q-Blindheit von J (lokales Gleichgewicht = Bahn), ESU-Struktur von Conj_Omega und
        reduzierte Cartan-Verbindungen im AdS- und dS-Modell, Gaußsches Lagelemma (jedes reine Gaußprofil liegt offen in
        einer Stabilisatorbahn), q-Gradient und quadratische Antwort (BKM), Komposita an verschiedenen Orten bzw. mit
        q-Gefälle (K3), Z4-Überlagerung der Markenvakua (K4, native Amplitude).
Teil 3  Die Kappe: dS-Bahn von tau (so(4,1), Vorzeichen der Krümmung), kein Nahtmodulfluss ist hyperbolisch (Kappe braucht
        das Typ-III-Netz), feste Familie in der dS-Klasse (SO(4) hält Omega fest, Kappenzeit -> Uhrbogen, tan phi = c sinh s,
        Grenzfall q -> 1 = konforme Kompaktifizierung von dS4 im Einsteinstreifen); Bilanz, Urteil, Folgen für K13 und die
        Gravitationsachse der Szenario-Suche.
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
import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
rng = np.random.default_rng(128)
I2 = np.eye(2)
SIGMA = [np.eye(2, dtype=complex), np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]),
         np.diag([1., -1.]).astype(complex)]
Q_TEST = 0.25                                                       # Prüfwert wie z100, z115
Q_NATIV = math.exp(-2*math.pi)                                      # KMS des Originals bei beta = 4 pi


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


def vec(M):
    return np.concatenate([np.asarray(M).real.ravel(), np.asarray(M).imag.ravel()])


# ------------------------------------------------------------------ u(2,2), so(3,2), so(4,1) auf T (wie z98, z104, z115)
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
            C = np.atleast_2d(f(M))
            parts += [C.real.ravel(), C.imag.ravel()]
        rows.append(np.concatenate(parts))
    ns = nullraum(np.array(rows).T)
    return [sum(c*M for c, M in zip(v, basis)) for v in ns]


def rest_im_spann(M, basis):
    A = np.array([vec(B) for B in basis]).T
    b = vec(M)
    x = np.linalg.lstsq(A, b, rcond=None)[0]
    return float(np.linalg.norm(A @ x - b))


uvec = [np.array([1, 0, -1, 0.]), np.array([0, 1, 0, -1.])]
vvec = [np.array([1, 0, 1, 0.]), np.array([0, 1, 0, 1.])]
Inp = (np.outer(uvec[0], uvec[1]) - np.outer(uvec[1], uvec[0])).astype(complex)       # I_pi
Ion = (np.outer(vvec[0], vvec[1]) - np.outer(vvec[1], vvec[0])).astype(complex)       # I_omega
NP = [-0.5*sum(s[i, j]*np.outer(uvec[i], vvec[j]) for i in range(2) for j in range(2)) for s in SIGMA]
NK = [-0.5*sum(s[i, j]*np.outer(vvec[i], uvec[j]) for i in range(2) for j in range(2)) for s in SIGMA]
ETA = np.diag([1., -1., -1., -1.])
MP = [1j*N for N in NP]                                             # P_mu
MKL = [ETA[m, m]*1j*NK[m] for m in range(4)]                        # K_mu (Index gesenkt)
LOR = teilalgebra([lambda M: M @ S4 - S4 @ M, lambda M: M @ Inp + Inp @ M.T])
QP = [(MP[m] + MKL[m])/2 for m in range(4)]                         # Q^+_mu (AdS), Q^+_0 = E
QM = [(MP[m] - MKL[m])/2 for m in range(4)]                         # Q^-_mu (dS), Q^-_0 = Kappe (z104)
SO32 = LOR + QP
SO41 = LOR + QM
ROT = teilalgebra([lambda M: kl(M, S4), lambda M: M @ Inp + Inp @ M.T, lambda M: kl(M, G4)])
BOOST = teilalgebra([lambda M: kl(M, S4), lambda M: M @ Inp + Inp @ M.T, lambda M: M @ G4 + G4 @ M])
ISYM = Inp - Ion                                                    # so(3,2) = Der(I_pi - I_omega)
ISYM_INV = np.linalg.inv(ISYM)
ISYMD = Inp + Ion                                                   # so(4,1) = Der(I_pi + I_omega)
ISYMD_INV = np.linalg.inv(ISYMD)
ME = QP[0]                                                          # Uhr E (elliptisch)
MQ0 = QM[0]                                                         # Kappe Q^- = (P_0 - K_0)/2 (hyperbolisch)
MD = kl(ME, MQ0)                                                    # Dilatation ([E, Q^-] = i D)
MP0 = MP[0]                                                         # P_0 = E + Q^- (parabolisch, Wiesbrock)


def iota(tau, I=ISYM, Iinv=ISYM_INV):
    """Invariante tr(tau I tau^T I^-1): mit I = I_pi - I_omega unter SO(3,2), mit I_pi + I_omega unter SO(4,1) invariant."""
    return complex(np.trace(tau @ I @ tau.T @ Iinv))


def Rw(wort):
    """U(g) = U(e^(t1 X1)) U(e^(t2 X2)) ...; U(g) Z U(g)^* = R Z mit R = e^(tn Xn) ... e^(t1 X1)."""
    R = np.eye(4, dtype=complex)
    for X, t in wort:
        R = sla.expm(t*X) @ R
    return R


def zufall(scale=0.5, basis=None):
    basis = SO32 if basis is None else basis
    return Rw([(sum(c*X for c, X in zip(rng.normal(size=len(basis)), basis))*scale, 1.0)])


def ebenen(q):
    s = math.sqrt(1 - q*q)
    return np.hstack([I2, q*I2]).astype(complex)/s, np.hstack([q*I2, I2]).astype(complex)/s


P_OM = np.hstack([I2, 0*I2]).astype(complex)                        # Vernichter von Omega in T: d_1, d_2


def ebene_z(z):
    """SO(3)-invariante positive Ebene: Vernichter d_i + z e_i^+ (Zustand im Punkt z der Scheibe)."""
    return np.hstack([I2, z*I2]).astype(complex)/math.sqrt(1 - abs(z)**2)


def alg_z(z):
    """Algebra mit Vernichterebene im Punkt z, Kommutante = G4-orthogonale negative Ebene."""
    n = math.sqrt(1 - abs(z)**2)
    return np.hstack([I2, z*I2]).astype(complex)/n, np.hstack([np.conj(z)*I2, I2]).astype(complex)/n


def z_von(P):
    return np.linalg.solve(P[:, :2], P[:, 2:])


def polar_u(L):
    U_, sv, Wh = np.linalg.svd(L)
    return U_ @ Wh, sv


def tomita1(A, B, P):
    """Einteilchen-Tomita (z115): Algebra a = A Z, Kommutante b^+ = B Z, Zustand mit Vernichtern P Z."""
    T = np.vstack([A, B])
    XY = P @ np.linalg.inv(T)
    Lam = -np.linalg.solve(XY[:, :2], XY[:, 2:])
    V, sv = polar_u(Lam)
    tau = np.linalg.inv(T) @ np.block([[np.zeros((2, 2)), V], [V.conj().T, np.zeros((2, 2))]]) @ T
    return tau, V, sv, Lam


def PhiP(R, P=P_OM, q=Q_TEST):
    """tau(M_g, psi): Algebra bewegt (A R, B R), Zustand mit Vernichterebene P fest."""
    A, B = ebenen(q)
    return tomita1(A @ R, B @ R, P)


def tau_von(R, P=P_OM, q=Q_TEST):
    return PhiP(R, P, q)[0]


def dtau(R0, X, P=P_OM, q=Q_TEST, h=1e-6):
    return (tau_von(sla.expm(h*X) @ R0, P, q) - tau_von(sla.expm(-h*X) @ R0, P, q))/(2*h)


def rang(R0, P=P_OM, q=Q_TEST, basis=None):
    basis = SO32 if basis is None else basis
    D = np.array([vec(dtau(R0, X, P, q)) for X in basis])
    sv = np.linalg.svd(D, compute_uv=False)
    return int(np.sum(sv > 1e-6*sv[0])), sv


def adapt(P):
    """G4-orthonormale Basis [P; P^perp] (P positiv, P^perp negativ): Zustandsangepasster Rahmen."""
    L = np.linalg.cholesky(P @ G4 @ P.conj().T)
    Pn = np.linalg.solve(L, P)
    Q = nullraum(Pn @ G4)                                           # Zeilen v mit Pn G4 v^+ = 0
    L2 = np.linalg.cholesky(-(Q @ G4 @ Q.conj().T))
    return np.vstack([Pn, np.linalg.solve(L2, Q)])


def conj_form(tau, P):
    T = adapt(P)
    tp = T @ tau @ np.linalg.inv(T)
    be = tp[:2, 2:]
    return float(max(np.abs(tp[:2, :2]).max(), np.abs(tp[2:, 2:]).max(), np.abs(be @ be.conj().T - I2).max(),
                     np.abs(tp[2:, :2] - be.conj().T).max(), np.abs(tp @ tp - np.eye(4)).max()))


def moebius(X, t):
    """Skalare 2x2-Blöcke von e^(t X): z -> (b + z d)/(a + z c) auf der Zustandsscheibe."""
    R = sla.expm(t*X)
    bl = [R[:2, :2], R[:2, 2:], R[2:, :2], R[2:, 2:]]
    dev = max(float(np.abs(B - B[0, 0]*I2).max()) for B in bl)
    a, b, c, d = (B[0, 0] for B in bl)
    return (a, b, c, d), dev


def d4(f, x, h):
    return (-f(x + 2*h) + 8*f(x + h) - 8*f(x - h) + f(x - 2*h))/(12*h)


def schleife(rahmen, s0=0.31, t0=-0.27, h=1e-3, H=1e-2, tau0=S4):
    """Schleife wie z115: Rahmen W mit tau = W^-1 tau0 W; theta = dW W^-1, Gamma = P_h theta, e = P_m theta."""
    def th(s_, t_, mu):
        W = rahmen(s_, t_)
        dW = d4(lambda x: rahmen(x, t_), s_, h) if mu == 0 else d4(lambda x: rahmen(s_, x), t_, h)
        return dW @ np.linalg.inv(W)
    Ph = lambda X: (X + tau0 @ X @ tau0)/2
    Pm = lambda X: (X - tau0 @ X @ tau0)/2
    ths, tht = th(s0, t0, 0), th(s0, t0, 1)
    d_s_tht = d4(lambda x: th(x, t0, 1), s0, H)
    d_t_ths = d4(lambda x: th(s0, x, 0), t0, H)
    F = d_s_tht - d_t_ths - kl(ths, tht)
    Rk = Ph(d_s_tht) - Ph(d_t_ths) - kl(Ph(ths), Ph(tht))
    Tk = Pm(d_s_tht) - Pm(d_t_ths) - kl(Ph(ths), Pm(tht)) - kl(Pm(ths), Ph(tht))
    ee = kl(Pm(ths), Pm(tht))
    return {"|F_MC|": float(np.abs(F).max()), "|R|": float(np.abs(Rk).max()), "|R - [e, e]|": float(np.abs(Rk - ee).max()),
            "|T|": float(np.abs(Tk).max())}


# ------------------------------------------------------------------ Fockraum einer Naht (z100/z115-Maschinerie)
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


class Naht:
    """Fockraum Original (a_1, a_2) x Partner (b_1, b_2), d = (a - q b^+)/s, e = (b - q a^+)/s, Omega = KMS-Purifikation."""

    def __init__(self, q, NV=11):
        self.q, self.NV = q, NV
        a1, a2, b1, b2 = moden(NV, 4)
        self.AO, self.BO = [a1, a2], [b1, b2]
        self.DIM, self.DA = NV**4, NV**2
        self.OCC = np.array(list(itertools.product(range(NV), repeat=4)))
        self.TOT = self.OCC.sum(1)
        self.NULL = sps.csr_matrix((self.DIM, self.DIM), dtype=complex)
        sq = math.sqrt(1 - q*q)
        self.d_ = [((self.AO[i] - q*dag(self.BO[i]))/sq).tocsr() for i in range(2)]
        self.e_ = [((self.BO[i] - q*dag(self.AO[i]))/sq).tocsr() for i in range(2)]
        self.ZT = [self.d_[0], self.d_[1], dag(self.e_[0]), dag(self.e_[1])]
        self.ZTd = [dag(x) for x in self.ZT]
        om1 = np.diag([sq*q**n for n in range(NV)])
        C_OM = np.kron(om1, om1)
        C_OM /= np.linalg.norm(C_OM)
        self.OMEGA = C_OM.reshape(-1).astype(complex)
        self.E0 = np.zeros(self.DIM, complex)                       # Fockvakuum der Originalmoden (|0_a 0_b>)
        self.E0[0] = 1
        self.HS = [self.gen(M) for M in SO32]
        self.HE = self.gen(ME)
        self.cE = complex(np.vdot(self.OMEGA, self.HE @ self.OMEGA))
        self.Mab = np.block([[I2, -q*I2], [-q*I2, I2]]).astype(complex)/sq          # [d; e^+] = Mab [a; b^+]

    def gauss(self, P):
        """Exakter Gaußvektor zur Vernichterebene P: exp(-a^+ Gamma b^+) |0_a 0_b>, Gamma = X^-1 Y aus P Mab = [X, Y].
        Auf dem tiefen Block exakt (der Erzeuger hebt nur an); Phase: <Omega|psi> > 0."""
        PX = P @ self.Mab
        Gam = np.linalg.solve(PX[:, :2], PX[:, 2:])
        G = sum((Gam[i, j]*(dag(self.AO[i]) @ dag(self.BO[j])) for i in range(2) for j in range(2)), self.NULL)
        v = spla.expm_multiply(-G.tocsc(), self.E0)
        ov = np.vdot(self.OMEGA, v)
        return v*abs(ov)/ov/np.linalg.norm(v)

    def gen(self, M):
        Am = 1j*G4 @ M
        H = self.NULL.copy()
        for al in range(4):
            for be in range(4):
                if abs(Am[al, be]) > 1e-14:
                    H = H + Am[al, be]*(self.ZTd[al] @ self.ZT[be])
        return H.tocsr()

    def uhr(self, theta, v):
        """U(e^(theta E)) mit Omega-Phase 1 (U Omega = Omega)."""
        return spla.expm_multiply(1j*theta*(self.HE - self.cE*sps.identity(self.DIM, format="csr")), v)

    def bkm(self, psi, gens, cut=1e-13):
        DA = self.DA
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


def dV_lyap(C, dC):
    """Ableitung des Polarteils V von C = P V (z115)."""
    U_, sv, Wh = np.linalg.svd(C)
    V = U_ @ Wh
    Y = dC @ V.conj().T
    Rhs = U_.conj().T @ (-1j*(Y - Y.conj().T)) @ U_
    K = U_ @ (Rhs/(sv[:, None] + sv[None, :])) @ U_.conj().T
    return 1j*K @ V, V


def tangenten(N, psi, L0c=3):
    """Tangenten der Familie L_g = U (V_g (x) V_g^T) SWAP U^T (J = L K) bei g = 1, Spalten N_tot <= L0c (z115)."""
    cols = np.where(N.TOT <= L0c)[0]
    cols2 = np.where(N.TOT <= L0c + 2)[0]
    pos2 = {j: k for k, j in enumerate(cols2)}
    C = psi.reshape(N.DA, N.DA)
    out = []
    for H in N.HS:
        dC = (-1j*(H @ psi)).reshape(N.DA, N.DA)
        dV, V = dV_lyap(C, dC)
        U2 = np.zeros((N.DIM, len(cols2)), complex)
        for k, j in enumerate(cols2):
            ia, ib = divmod(j, N.DA)
            U2[:, k] = np.kron(V[:, ib], V[ia, :])
        U1 = U2[:, [pos2[j] for j in cols]]
        bahn = 1j*(H @ U1) + 1j*(U2 @ H.conj()[cols2][:, cols].toarray())
        zust = np.zeros((N.DIM, len(cols)), complex)
        for k, j in enumerate(cols):
            ia, ib = divmod(j, N.DA)
            zust[:, k] = np.kron(dV[:, ib], V[ia, :]) + np.kron(V[:, ib], dV[ia, :])
        out.append(bahn + zust)
    M = np.array([vec(t) for t in out]).T
    _, sv, vh = np.linalg.svd(M, full_matrices=False)
    return sv, vh


BR = np.array([vec(M) for M in SO32]).T
koord = lambda M: np.linalg.lstsq(BR, vec(M), rcond=None)[0]
ROT_K = np.array([koord(R_) for R_ in ROT])


def cos_kern_rot(vh_kern):
    return float(np.linalg.svd(np.linalg.qr(vh_kern.T)[0].T @ np.linalg.qr(ROT_K.T)[0], compute_uv=False).min())


NAHT = Naht(Q_TEST)
laufzeit()

# =================================================================== Teil 1
say("Teil 1: Native inhomogene Zustände")

# --- 1.1 Zustandsscheibe
say("  1.1 Zustandsscheibe: SO(3)-invariante Gaußzustände, radiale sl(2,R) als Möbiusgruppe, Schmidtwert")
inv_rot = 0.0
for _ in range(5):
    z = rng.uniform(0, 0.8)*np.exp(1j*rng.uniform(0, 2*np.pi))
    for X in ROT:
        inv_rot = max(inv_rot, float(np.abs(z_von(ebene_z(z) @ sla.expm(0.7*X)) - z*I2).max()))
Pnicht = np.hstack([I2, np.diag([0.3, -0.2])]).astype(complex)
nicht_rot = float(np.abs(z_von(Pnicht @ sla.expm(0.7*ROT[0])) - z_von(Pnicht)).max())
mob = {}
for name_, X in (("E", ME), ("Q^-", MQ0), ("D", MD), ("P_0", MP0)):
    (a_, b_, c_, d_), dev = moebius(X, 0.7)
    zt = [rng.uniform(0, 0.7)*np.exp(1j*rng.uniform(0, 2*np.pi)) for _ in range(4)]
    tr_err = max(abs(z_von(ebene_z(z) @ sla.expm(0.7*X))[0, 0] - (b_ + z*d_)/(a_ + z*c_)) for z in zt)
    fp = np.roots([c_, a_ - d_, -b_]) if abs(c_) > 1e-14 else np.array([])
    mob[name_] = {"Blockabweichung": dev, "Möbius-Fehler": float(tr_err), "Fixpunkte": [complex(np.round(x, 9)) for x in fp],
                  "Bild von 0": complex(np.round(b_/a_, 9))}
rot_E = max(abs((b_ + z*d_)/(a_ + z*c_) - np.exp(0.7j)*z) for (a_, b_, c_, d_) in [moebius(ME, 0.7)[0]]
            for z in (0.3, 0.2j, 0.1 + 0.4j))
fixQ = sorted(mob["Q^-"]["Fixpunkte"], key=lambda x: x.imag)
fixD = sorted(mob["D"]["Fixpunkte"], key=lambda x: x.real)
fixP = mob["P_0"]["Fixpunkte"]
z_alg = complex(z_von(ebenen(Q_TEST)[0])[0, 0])
lam_err, sv_err = 0.0, 0.0
for _ in range(20):
    zA = rng.uniform(0, 0.85)*np.exp(1j*rng.uniform(0, 2*np.pi))
    zP = rng.uniform(0, 0.85)*np.exp(1j*rng.uniform(0, 2*np.pi))
    A_, B_ = alg_z(zA)
    _, _, sv_, Lam_ = tomita1(A_, B_, ebene_z(zP))
    lam = (zA - zP)/(1 - np.conj(zA)*zP)
    lam_err = max(lam_err, float(np.abs(Lam_ - lam*I2).max()))
    sv_err = max(sv_err, float(np.abs(sv_ - abs(lam)).max()))
sv_split = float(np.ptp(tomita1(*alg_z(0.3), Pnicht/np.sqrt(1 - 0.09))[2]))
ZAHLEN["scheibe"] = {"SO(3)-Invarianz von P_z": inv_rot, "Kontrolle nicht skalare Ebene": nicht_rot, "Möbius": mob,
                     "E dreht z -> e^(it) z": float(rot_E), "Algebra M im Punkt": z_alg,
                     "|Lambda - lambda 1|, lambda = (z_A - z_P)/(1 - conj(z_A) z_P)": lam_err, "Schmidtwerte": sv_err,
                     "Aufspaltung nicht invariantes Paar": sv_split}
check(f"[E] Zustandsscheibe: Die SO(3)-invarianten reinen Gaußzustände sind die Ebenen P_z = [1, z 1] (Vernichter "
      f"d_i + z e_i^+), |z| < 1 (Invarianz {inv_rot:.0e}; nicht skalare Kontrollebene bewegt sich um {nicht_rot:.3f}). Die "
      f"radiale sl(2,R) wirkt mit skalaren Blöcken (Abweichung {max(v['Blockabweichung'] for v in mob.values()):.0e}) als "
      f"Möbiusgruppe der Scheibe (Formelfehler {max(v['Möbius-Fehler'] for v in mob.values()):.0e}): E dreht um den Ursprung "
      f"(z -> e^(it) z, {rot_E:.0e}), die Kappe Q^- verschiebt den Ursprung nach -i tanh(t/2) (Bild {mob['Q^-']['Bild von 0']}) "
      f"mit den Randfixpunkten {fixQ}, D nach -tanh(t/2) mit {fixD}, P_0 ist parabolisch mit dem doppelten Fixpunkt "
      f"{fixP[0] if fixP else None}. Die Fixpunkte von Q^- und D sind zusammen genau z^4 = 1, die Marken bzw. Spitzen der "
      f"tau-Eigenkappen (z104). Das Vakuum Omega sitzt im Ursprung, die Algebra M im Punkt z = {z_alg.real:.2f} = q auf der "
      f"D-Achse. Relative Lage: Lambda = lambda 1 mit lambda = (z_A - z_P)/(1 - conj(z_A) z_P) (20 Zufallspaare, "
      f"{lam_err:.0e}), der Schmidtwert ist der pseudohyperbolische Abstand |lambda| ({sv_err:.0e}); ein nicht SO(3)-invariantes "
      f"Paar spaltet die Schmidtwerte auf ({sv_split:.3f})",
      inv_rot < 1e-12 and nicht_rot > 1e-2 and max(v["Blockabweichung"] for v in mob.values()) < 1e-12
      and max(v["Möbius-Fehler"] for v in mob.values()) < 1e-12 and rot_E < 1e-12
      and np.allclose(fixQ, [-1j, 1j], atol=1e-8) and np.allclose(fixD, [-1, 1], atol=1e-8)
      and len(fixP) == 2 and np.allclose(fixP, [1, 1], atol=1e-6)
      and abs(mob["Q^-"]["Bild von 0"] + 1j*math.tanh(0.35)) < 1e-8 and abs(mob["D"]["Bild von 0"] + math.tanh(0.35)) < 1e-8
      and abs(z_alg - Q_TEST) < 1e-14 and lam_err < 1e-12 and sv_err < 1e-12 and sv_split > 0.05)
laufzeit()

# --- 1.2 Native Profile
say("  1.2 Native Profile: Hislop-Longo-Temperatur (K1a), Wiesbrock-Profil (K1b), Markenvakua (K2), Kappenzeit (frei)")
ts, rs = sp.symbols("t r", real=True)
xi0, xir = (1 - ts**2 - rs**2)/2, -ts*rs                            # Kappenfeld xi_Q (z104)
norm2 = sp.factor(sp.expand(xi0**2 - xir**2))
T_r = sp.simplify(1/(2*sp.pi*sp.sqrt(norm2.subs(ts, 0))))
T_soll = 1/(sp.pi*sp.Abs(1 - rs**2))
hl_ok = sp.simplify(T_r - T_soll) == 0 and T_r.subs(rs, 0) == 1/sp.pi
wies = {}
for qq in (Q_TEST, Q_NATIV):
    A_, B_ = ebenen(qq)
    err, ph_, io_ = 0.0, [], []
    for a in (0.3, 1.0, 3.0, 10.0):
        tau_, _, sv_, _ = tomita1(A_ @ sla.expm(a*MP0), B_ @ sla.expm(a*MP0), P_OM)
        x = a/2
        err = max(err, float(np.abs(sv_ - math.sqrt((qq*qq + x*x*(1 - qq)**2)/(1 + x*x*(1 - qq)**2))).max()))
        be = tau_[:2, 2:]
        ph_.append(float(np.abs(be - be[0, 0]*I2).max()))
        io_.append(abs(iota(tau_) - 4))
    fern = float(tomita1(A_ @ sla.expm(1e4*MP0), B_ @ sla.expm(1e4*MP0), P_OM)[2].min())
    wies[qq] = {"Formelfehler q(a)": err, "beta skalar": max(ph_), "|iota - 4|": max(io_), "q(a = 10^4)": fern}
psi_uhr = NAHT.uhr(np.pi/2, NAHT.E0)                                 # U(e^(pi E/2)) |0_a 0_b> (Fock, Abschneidung)
low6 = NAHT.TOT <= 6
annih = [float(np.linalg.norm(((NAHT.ZT[i] + 1j*Q_TEST*NAHT.ZT[2 + i]) @ psi_uhr)[low6])) for i in range(2)]
annih_ko = float(np.linalg.norm(((NAHT.ZT[0] - 1j*Q_TEST*NAHT.ZT[2]) @ psi_uhr)[low6]))
P_MARK = ebene_z(1j*Q_TEST)
psi_mark = NAHT.gauss(P_MARK)                                        # exakter Gaußvektor des Markenvakuums
mark_gauss = float(np.abs((psi_uhr*np.vdot(psi_uhr, psi_mark)/abs(np.vdot(psi_uhr, psi_mark)) - psi_mark)[low6]).max())
mark_eben = float(np.abs(z_von(ebenen(Q_TEST)[0] @ sla.expm(np.pi/2*ME)) - 1j*Q_TEST*I2).max())
prof_err, prof = 0.0, []
for phi in (0.0, 0.6, 1.2, np.pi/2 - 1e-3, 2.5, -np.pi/2):
    sv_ = PhiP(sla.expm(phi*ME), P_MARK)[2]
    soll = Q_TEST*abs(np.exp(1j*phi) - 1j)/abs(1 - 1j*Q_TEST**2*np.exp(-1j*phi))
    prof_err = max(prof_err, float(np.abs(sv_ - soll).max()))
    prof.append(float(sv_[0]))
kap_z = {s: complex(z_von(P_OM @ sla.expm(s*MQ0))[0, 0]) for s in (0.3, 1.0)}
kap_err = max(abs(v + 1j*math.tanh(s/2)) for s, v in kap_z.items())
ZAHLEN["profile"] = {"|xi_Q|^2": str(norm2), "T(r) bei t = 0": str(T_r), "Wiesbrock q(a)^2 = (q^2 + x^2 (1-q)^2)/(1 + x^2 (1-q)^2)":
                     {str(k): v for k, v in wies.items()}, "Markenvakuum: |(d + i q e^+) psi|": annih,
                     "Kontrolle |(d - i q e^+) psi|": annih_ko, "Markenebene z = i q": mark_eben,
                     "exakter Gaußvektor gegen U(e^(pi E/2)) |0> (tiefer Block)": mark_gauss,
                     "Profil q(phi) (phi = 0; 0,6; 1,2; pi/2 - 10^-3; 2,5; -pi/2)": prof, "Profilformelfehler": prof_err,
                     "Kappenzustand z_s = -i tanh(s/2)": float(kap_err)}
check(f"[E/C] Native Profile: (K1a) Das Kappenfeld hat |xi_Q|^2 = {norm2}; im Vakuum (KMS bei beta = 2 pi bezüglich "
      f"K_cap, z104) ist die lokale Temperatur auf dem Schnitt t = 0 T(r) = 1/(2 pi |xi_Q|) = {T_r} (sympy, "
      f"T(0) = 1/pi): ein natives radiales Profil, bedingt auf das freie Netz; als lokales q braucht es eine gewählte "
      f"Nahtfrequenz. (K1b) Die mit der Wiesbrock-Zeit P_0 = E + Q^- verschobenen Algebren N_a sehen den festen "
      f"Nahtzustand mit q(a)^2 = (q^2 + x^2 (1 - q)^2)/(1 + x^2 (1 - q)^2), x = a/2 (Fehler "
      f"{max(v['Formelfehler q(a)'] for v in wies.values()):.0e} bei q = 0,25 und nativem q), q(a) -> 1 (bei a = 10^4: "
      f"{wies[Q_TEST]['q(a = 10^4)']:.9f}): ein natives Profil ohne freien Parameter, das zum parabolischen Fixpunkt z = 1 "
      f"hin maximal verschränkt (lokale Temperatur divergiert, Tolman-artig); J bleibt dabei auf dem AdS-Uhrkreis "
      f"(beta skalar {max(v['beta skalar'] for v in wies.values()):.0e}, "
      f"|iota - 4| {max(v['|iota - 4|'] for v in wies.values()):.0e}). (K2) Das Fockvakuum der Markenalgebra M_(pi/2) "
      f"(Viertelmarke, z83) ist U(e^(pi E/2)) |0_a 0_b> und hat die Vernichter d + i q e^+ (Fock {max(annih):.0e}; "
      f"Kontrolle d - i q e^+: {annih_ko:.3f}), also den Scheibenpunkt i q ({mark_eben:.0e}; exakter Gaußvektor "
      f"exp(-a^+ Gamma b^+) |0> bis auf die Phase gleich, {mark_gauss:.0e} im tiefen Block); die uhrgedrehten Algebren "
      f"M_phi sehen es mit q(phi) = q |e^(i phi) - i|/|1 - i q^2 e^(-i phi)| (Fehler {prof_err:.0e}; Werte "
      f"{np.round(prof, 5).tolist()}), mit einem Verschränkungshorizont an der eigenen Marke phi = pi/2. Der "
      f"kappenverschobene Zustand U(e^(s Q^-)) Omega sitzt bei -i tanh(s/2) ({kap_err:.0e}); s ist frei",
      hl_ok and max(v["Formelfehler q(a)"] for v in wies.values()) < 1e-12 and 1 - wies[Q_TEST]["q(a = 10^4)"] < 1e-6
      and max(v["beta skalar"] for v in wies.values()) < 1e-12 and max(v["|iota - 4|"] for v in wies.values()) < 1e-9
      and max(annih) < 1e-7 and annih_ko > 0.1 and mark_eben < 1e-12 and mark_gauss < 1e-7 and prof_err < 1e-12
      and prof[3] < 1e-2
      and kap_err < 1e-12)
befund("E/C: Native inhomogene Zustände existieren. Nativ ohne freien Parameter sind das Wiesbrock-Profil q(a) (fester "
       "Nahtzustand, mit P_0 verschobene Algebren) und die Markenvakua (Fockvakua der Markenalgebren, unter z83s "
       "Markenidentifikation; eine Zustandswahl: unverschränkt statt KMS) mit dem Uhrprofil q(phi) und Horizont an der "
       "eigenen Marke; nativ bedingt auf das Netz ist das Hislop-Longo-Profil T(r) = 1/(pi (1 - r^2)), dessen Übersetzung in "
       "ein lokales q eine Nahtfrequenz wählt; frei ist die Kappenzeit s des kappenverschobenen Zustands. Komposita (K3) und "
       "die Z4-Überlagerung der Markenvakua (K4) folgen in Teil 2", True)
laufzeit()

# =================================================================== Teil 2
say("\nTeil 2: Geometrie der Familien J(M_g, psi)")

# --- 2.1 q-Blindheit
say("  2.1 q-Blindheit: lokales Gleichgewicht (K1) ist die native Bahn")
blind = 0.0
for qq in list(rng.uniform(0.01, 0.95, size=10)) + [Q_NATIV]:
    A_, B_ = ebenen(qq)
    blind = max(blind, float(np.abs(tomita1(A_, B_, P_OM)[0] - S4).max()))
lgg = 0.0
for _ in range(8):
    R_ = zufall(0.6)
    qx = rng.uniform(0.02, 0.9)
    lgg = max(lgg, float(np.abs(PhiP(R_, P_OM @ R_, qx)[0] - np.linalg.inv(R_) @ S4 @ R_).max()))
# Modulgenerator der Naht: K_a = -log rho_a = omega N_a + const, also [K_a, a_i] = -omega a_i (zentral auf T)
C_ = NAHT.OMEGA.reshape(NAHT.DA, NAHT.DA)
p_, U_ = np.linalg.eigh(C_ @ C_.conj().T)
keep = p_ > 1e-12
Ka = -(U_[:, keep]*np.log(p_[keep])) @ U_[:, keep].conj().T
aa = moden(NAHT.NV, 2)
lowa = np.where(np.array(list(itertools.product(range(NAHT.NV), repeat=2))).sum(1) <= 4)[0]


def ein_teilchen_L(Ka_):
    L = np.zeros((2, 2), complex)
    res, kreuz = 0.0, 0.0
    for i in range(2):
        comm = Ka_ @ aa[i].toarray() - aa[i].toarray() @ Ka_
        A = np.array([(aa[j].toarray()[:, lowa]).ravel() for j in range(2)]).T
        Ad_ = np.array([(aa[j].T.conj().toarray()[:, lowa]).ravel() for j in range(2)]).T
        b = (-comm[:, lowa]).ravel()
        x = np.linalg.lstsq(np.hstack([A, Ad_]), b, rcond=None)[0]
        L[i] = x[:2]
        kreuz = max(kreuz, float(np.abs(x[2:]).max()))
        res = max(res, float(np.linalg.norm(np.hstack([A, Ad_]) @ x - b)/np.linalg.norm(b)))
    return L, res, kreuz


L_om, res_om, kreuz_om = ein_teilchen_L(Ka)
omega_soll = 2*math.log(1/Q_TEST)
kipp = float(np.abs(tau_von(np.eye(4), P_OM @ sla.expm(0.05*MQ0)) - S4).max())
kipp_boost = float(np.abs(tau_von(np.eye(4), P_OM @ sla.expm(0.05*BOOST[0])) - S4).max())
ZAHLEN["q_blind"] = {"|tau(M_q, Omega_q) - S4| (11 Werte von q)": blind, "lokales Gleichgewicht gegen Bahn": lgg,
                     "Einteilchen-Modulgenerator L (Fock)": np.round(L_om, 9).tolist(), "omega = 2 log(1/q)": omega_soll,
                     "Fitrest, a^+-Anteil": [res_om, kreuz_om], "Kontrolle: Zustandskippung längs der Kappe": kipp,
                     "Zustandskippung längs eines Boosts (vertauscht mit tau_0)": kipp_boost}
check(f"[E] q-Blindheit (K1): tau(M_q, Omega_q) = S4 für jedes q (11 Werte, auch nativ: {blind:.0e}); eine Familie lokaler "
      f"Gleichgewichte x -> (U(g_x) M_q(x) U^*, U(g_x) Omega_q(x)) mit beliebigem Profil q(x) hat genau die Bahnfamilie "
      f"R_x^-1 S4 R_x ({lgg:.0e}): die native AdS-Bahn, nach z100 Maurer-Cartan-flach. Das Profil sitzt allein in Delta: "
      f"Der Fock-Modulgenerator K_a = -log rho_a erfüllt [K_a, a_i] = -sum L_ij a_j mit L = omega 1 (L = "
      f"{np.round(np.diag(L_om).real, 6).tolist()}, omega = 2 log(1/q) = {omega_soll:.6f}, Fitrest {res_om:.0e}, "
      f"a^+-Anteil {kreuz_om:.0e}): ein Skalar je Ort, zentral auf T, ohne Verbindung. Kontrolle: Eine Kippung des "
      f"Zustands längs der Kappe (statt einer q-Änderung) bewegt J ({kipp:.3f}); eine Kippung längs eines Boosts nicht "
      f"({kipp_boost:.0e}), weil Lorentz mit tau_0 vertauscht",
      blind < 1e-12 and lgg < 1e-12 and np.abs(L_om - omega_soll*I2).max() < 1e-5 and res_om < 1e-6 and kreuz_om < 1e-9
      and kipp > 1e-3 and kipp_boost < 1e-12)
laufzeit()

# --- 2.2 ESU-Struktur von Conj_Omega, reduzierte Cartan-Verbindungen
say("  2.2 Struktur von Conj_Omega (Einsteinzylinder S^1_E x S^3_(Q^-)) und reduzierte Cartan-Verbindungen")
KB32 = teilalgebra([lambda M: kl(M, G4)], basis=SO32)
KB41 = teilalgebra([lambda M: kl(M, G4)], basis=SO41)
blk = lambda M: float(max(np.abs(M[:2, 2:]).max(), np.abs(M[2:, :2]).max()))
diag_ = max(blk(M) for M in [ME] + QM[1:] + ROT)
quetsch = min(blk(M) for M in QP[1:] + [MQ0] + BOOST)
c_q = (1 + Q_TEST**2)/(2*Q_TEST)
T_fam = [dtau(np.eye(4), X) for X in QP]
T_esu = [kl(X, S4) for X in [ME] + QM[1:]]
T_ads = [kl(X, S4) for X in QP]
koef = np.linalg.lstsq(np.array([vec(M) for M in T_esu]).T, np.array([vec(M) for M in T_fam]).T, rcond=None)[0]
koef_soll = np.diag([-1, c_q, c_q, c_q])
cos_esu = np.linalg.svd(np.linalg.qr(np.array([vec(M) for M in T_fam]).T)[0].T
                        @ np.linalg.qr(np.array([vec(M) for M in T_esu]).T)[0], compute_uv=False)
cos_ads = np.sort(np.linalg.svd(np.linalg.qr(np.array([vec(M) for M in T_fam]).T)[0].T
                                @ np.linalg.qr(np.array([vec(M) for M in T_ads]).T)[0], compute_uv=False))[::-1]


def komplement(gbasis):
    M = np.array([[np.trace(Y @ g).real for g in gbasis] for Y in U22])
    return [sum(c*Y for c, Y in zip(v, U22)) for v in nullraum(M.T)]


def projektor(gbasis):
    pbasis = komplement(gbasis)
    Bm = np.array([vec(M) for M in gbasis + pbasis]).T

    def P(X):
        c = np.linalg.lstsq(Bm, vec(X), rcond=None)[0]
        return (sum(ci*M for ci, M in zip(c[:len(gbasis)], gbasis)), sum(ci*M for ci, M in zip(c[len(gbasis):], pbasis)))
    return P, len(pbasis)


PR32, dimp32 = projektor(SO32)
PR41, dimp41 = projektor(SO41)
XL, YL = QP[1] + 0.3*QP[0] + 0.2*LOR[2], QP[2] - 0.5*QP[3] + 0.4*LOR[4]


def rahmen_fest(s_, t_, P=P_OM):
    tp = tau_von(Rw([(XL, s_), (YL, t_)]), P)
    if P is not P_OM:
        T = adapt(P)
        tp = T @ tp @ np.linalg.inv(T)
    w = sla.sqrtm(tp[:2, 2:])
    return sla.block_diag(np.linalg.inv(w), w)


def reduziert(PR, transl, s0=0.31, t0=-0.27, h=1e-3, H=1e-2):
    def th(s_, t_, mu):
        W = rahmen_fest(s_, t_)
        dW = d4(lambda x: rahmen_fest(x, t_), s_, h) if mu == 0 else d4(lambda x: rahmen_fest(s_, x), t_, h)
        return dW @ np.linalg.inv(W)
    A = lambda s_, t_, mu: PR(th(s_, t_, mu))[0]
    As, At = A(s0, t0, 0), A(s0, t0, 1)
    F = d4(lambda x: A(x, t0, 1), s0, H) - d4(lambda x: A(s0, x, 0), t0, H) - kl(As, At)
    Bs, Bt = PR(th(s0, t0, 0))[1], PR(th(s0, t0, 1))[1]
    gauss = float(np.abs(F - PR(kl(Bs, Bt))[0]).max())
    Bm = np.array([vec(M) for M in LOR + transl]).T
    cF = np.linalg.lstsq(Bm, vec(F), rcond=None)[0]
    R0 = Rw([(XL, s0), (YL, t0)])

    def rahmen_R(R):
        w = sla.sqrtm(tau_von(R)[:2, 2:])
        return sla.block_diag(np.linalg.inv(w), w)
    W0i = np.linalg.inv(rahmen_R(R0))
    loet = []
    for X in SO32:
        th_ = (rahmen_R(sla.expm(1e-6*X) @ R0) - rahmen_R(sla.expm(-1e-6*X) @ R0))/2e-6 @ W0i
        loet.append(np.linalg.lstsq(Bm, vec(PR(th_)[0]), rcond=None)[0][6:])
    sv_l = np.linalg.svd(np.array(loet), compute_uv=False)
    return {"|F|": float(np.abs(F).max()), "|F - P[B, B]|": gauss, "Lorentzteil": float(np.abs(cF[:6]).max()),
            "Torsion (Translationsteil)": float(np.abs(cF[6:]).max()),
            "Lötform-Singulärwerte": np.round(sv_l, 9).tolist(), "Lötrang": int(np.sum(sv_l > 1e-6*max(1.0, sv_l[0])))}


red_ads = reduziert(PR32, QP)
red_ds = reduziert(PR41, QM)
sch_fest = schleife(rahmen_fest)
ZAHLEN["esu"] = {"dim k_big n so(3,2), k_big n so(4,1)": [len(KB32), len(KB41)],
                 "E in k_big n so(3,2)": rest_im_spann(ME, KB32), "Q^-_k in k_big n so(4,1)": max(rest_im_spann(X, KB41) for X in QM[1:]),
                 "blockdiagonal (E, Q^-_k, Drehungen): max Nebenblock": diag_, "Quetscher (Q^+_k, Q^-_0, Boosts): min Nebenblock": quetsch,
                 "Tangente: Koeffizienten in (E, Q^-_1..3)": np.round(koef, 9).tolist(), "cos Hauptwinkel gegen ESU": np.round(cos_esu, 12).tolist(),
                 "cos Hauptwinkel gegen AdS-Bahn": np.round(cos_ads, 12).tolist(), "dim Komplemente": [dimp32, dimp41],
                 "reduziert AdS": red_ads, "reduziert dS": red_ds, "eigenes Modell (Rahmen in K_big)": sch_fest}
check(f"[E] Struktur von Conj_Omega: k_big n so(3,2) = so(3) + <E> (Dimension {len(KB32)}), k_big n so(4,1) = so(3) + "
      f"<Q^-_1, Q^-_2, Q^-_3> = so(4) (Dimension {len(KB41)}); E, Q^-_k und die Drehungen sind blockdiagonal (Nebenblock "
      f"{diag_:.0e}), Q^+_k, die Kappe Q^-_0 und die Boosts sind Quetscher (Nebenblock >= {quetsch:.3f}). Die Tangente der "
      f"Vakuumfamilie bei tau_0 ist dtau(E) = -[E, tau_0], dtau(Q^+_k) = c [Q^-_k, tau_0] mit c = coth(2r) = {c_q:.4f} "
      f"(Rest {np.abs(koef - koef_soll).max():.0e}; cos der Hauptwinkel gegen span{{[E, tau_0], [Q^-_k, tau_0]}} "
      f"{np.round(cos_esu, 9).tolist()}, gegen die AdS-Bahn {np.round(cos_ads, 9).tolist()}): Die AdS-Translationen der "
      f"Algebra werden auf die dS-Raumtranslationen der kompakten Dualen abgebildet. Conj_Omega = (S^1_E x S^3_(Q^-))/Z_2 "
      f"ist der kompaktifizierte Einsteinzylinder mit Zeit aus AdS und Raum aus dS. Reduzierte Cartan-Verbindung "
      f"A = P_g(dW W^-1) mit Rahmen in "
      f"K_big: im AdS-Modell |F| = {red_ads['|F|']:.4f} = P[B, B] (Rest {red_ads['|F - P[B, B]|']:.0e}), ganz im "
      f"Lorentzteil (die S^3-Berry-Krümmung von z115, |R| = {sch_fest['|R|']:.4f}), Torsion "
      f"{red_ads['Torsion (Translationsteil)']:.0e}, Lötform vom Rang {red_ads['Lötrang']} (nur E); im dS-Modell "
      f"|F| = {red_ds['|F|']:.0e}, Lötform vom Rang {red_ds['Lötrang']} (nur Q^-_k). In keinem der beiden Modelle ist die "
      f"Familie eine vierdimensionale Cartan-Geometrie; im eigenen Modell ist sie Maurer-Cartan-flach ({sch_fest['|F_MC|']:.0e})",
      len(KB32) == 4 and len(KB41) == 6 and rest_im_spann(ME, KB32) < 1e-9 and max(rest_im_spann(X, KB41) for X in QM[1:]) < 1e-9
      and diag_ < 1e-12 and quetsch > 0.3 and np.abs(koef - koef_soll).max() < 1e-6 and np.all(np.abs(cos_esu - 1) < 1e-9)
      and abs(cos_ads[0] - 1) < 1e-9 and cos_ads[1] < 1e-6 and dimp32 == 6 and dimp41 == 6
      and red_ads["|F - P[B, B]|"] < 1e-5 and abs(red_ads["|F|"] - sch_fest["|R|"]) < 1e-5 and red_ads["|F|"] > 0.1
      and red_ads["Torsion (Translationsteil)"] < 1e-6 and red_ads["Lötrang"] == 1 and red_ds["|F|"] < 1e-5
      and red_ds["Lötrang"] == 3 and sch_fest["|F_MC|"] < 1e-5)
laufzeit()

# --- 2.3 Gaußsches Lagelemma
say("  2.3 Gaußsches Lagelemma: jedes reine Gaußprofil liegt offen in einer Stabilisatorbahn Conj_psi = U(2)")
Rh = sla.expm(sum(c*X for c, X in zip(rng.normal(size=16), U22))*0.3)
ZUST = {"Markenvakuum (z = i q)": P_MARK, "Kappenzustand (s = 0,3)": P_OM @ sla.expm(0.3*MQ0), "Zufallstranslat U(h) Omega": P_OM @ Rh}
lage = {}
for name_, P in ZUST.items():
    Rx = None
    if name_.startswith("Zufall"):
        Rx = Rh
    elif name_.startswith("Kappe"):
        Rx = sla.expm(0.3*MQ0)
    kov = 0.0
    if Rx is not None:
        A_, B_ = ebenen(Q_TEST)
        for _ in range(5):
            Rg = zufall(0.5)
            S_ = Rg @ np.linalg.inv(Rx)
            kov = max(kov, float(np.abs(tomita1(A_ @ Rg, B_ @ Rg, P_OM @ Rx)[0]
                                        - np.linalg.inv(Rx) @ tomita1(A_ @ S_, B_ @ S_, P_OM)[0] @ Rx).max()))
    form = max(conj_form(tau_von(zufall(0.5), P), P) for _ in range(10))
    raenge = [rang(np.eye(4), P)[0]] + [rang(zufall(0.5), P)[0] for _ in range(3)]
    io_ = [iota(tau_von(zufall(0.5), P)).real for _ in range(8)]
    lage[name_] = {"Kovarianz": kov, "Conj-Form": form, "Rang (g = 1, drei Zufallspunkte)": raenge,
                   "iota min/max": [float(min(io_)), float(max(io_))]}
sch_mark = schleife(lambda s_, t_: rahmen_fest(s_, t_, P_MARK))
ZAHLEN["lagelemma"] = {"Zustände": lage, "Schleife Markenvakuum (Rahmen in K_big^psi)": sch_mark}
check(f"[E] Gaußsches Lagelemma: Für psi = U(h) Omega gilt tau(M_g, psi) = R_h^-1 tau(M_(h^-1 g), Omega) R_h "
      f"(Kappenzustand und Zufallstranslat, {max(v['Kovarianz'] for v in lage.values()):.0e}): Die Familie ist die "
      f"Vakuumfamilie über der verschobenen Nebenklasse h^-1 G, konjugiert. In der zustandsangepassten Basis [P; P^perp] "
      f"hat jedes Glied die Form [[0, beta], [beta^+, 0]], beta in U(2) ({max(v['Conj-Form'] for v in lage.values()):.0e}); "
      f"der Rang ist überall 4 ({ {k: v['Rang (g = 1, drei Zufallspunkte)'] for k, v in lage.items()} }): Die Familie "
      f"ist offen in Conj_psi, einer einzigen Bahn der Stabilisatorgruppe K_big^psi = U(2) x U(2). iota variiert "
      f"({ {k: np.round(v['iota min/max'], 3).tolist() for k, v in lage.items()} }): kein nativer Rahmen. Eigenes Modell "
      f"(Markenvakuum, Rahmen in K_big^psi): Maurer-Cartan {sch_mark['|F_MC|']:.0e}, Berry-Krümmung |R| = "
      f"{sch_mark['|R|']:.4f} = [e, e] (Rest {sch_mark['|R - [e, e]|']:.0e}), Torsion {sch_mark['|T|']:.0e}. Ein "
      f"Ortsprofil innerhalb der reinen Gaußzustände ändert also nur die Karte g -> beta, nie die Geometrie",
      max(v["Kovarianz"] for v in lage.values()) < 1e-12 and max(v["Conj-Form"] for v in lage.values()) < 1e-12
      and all(v["Rang (g = 1, drei Zufallspunkte)"] == [4, 4, 4, 4] for v in lage.values())
      and all(v["iota min/max"][1] - v["iota min/max"][0] > 0.1 for v in lage.values())
      and sch_mark["|F_MC|"] < 1e-5 and sch_mark["|R - [e, e]|"] < 1e-5 and sch_mark["|T|"] < 1e-5 and sch_mark["|R|"] > 1e-2)
laufzeit()

# --- 2.4 q-Gradient und quadratische Antwort
say("  2.4 q-Gradient: Kopplung von J-Bewegung und Spektrum, quadratische Antwort (BKM) auf der Familie")


def bewegungen(P, q=Q_TEST, h=1e-6):
    dT, dN = [], []
    for X in SO32:
        tp, _, _, Lp = PhiP(sla.expm(h*X), P, q)
        tm, _, _, Lm = PhiP(sla.expm(-h*X), P, q)
        dT.append(vec((tp - tm)/(2*h)))
        dN.append(vec((Lp.conj().T @ Lp - Lm.conj().T @ Lm)/(2*h)))
    dT, dN = np.array(dT), np.array(dN)
    rk = lambda M: int(np.sum(np.linalg.svd(M, compute_uv=False) > 1e-6*np.linalg.svd(M, compute_uv=False)[0]))
    rJ, rN = rk(dT), rk(dN)
    rG = rk(np.hstack([dT/np.abs(dT).max(), dN/np.abs(dN).max()]))
    return rJ, rN, rG


def mu_induziert(N, psi, P, cut=1e-13):
    """Induzierte BKM auf der Familie: min BKM(X) über alle X mit gegebener J-Tangente, relativ zur Spurmetrik von U(2)."""
    Fm = N.bkm(psi, N.HS, cut=cut)
    D = [dtau(np.eye(4), X, P, N.q) for X in SO32]
    Jm = np.array([vec(M) for M in D])
    _, sv, Vh = np.linalg.svd(Jm.T, full_matrices=False)
    r = int(np.sum(sv > 1e-7*sv[0]))
    K, Xp = Vh[r:].T, Vh[:r].T
    S = Xp.T @ Fm @ Xp - Xp.T @ Fm @ K @ np.linalg.pinv(K.T @ Fm @ K, rcond=1e-10) @ K.T @ Fm @ Xp
    Vimg = [sum(Xp[i, j]*D[i] for i in range(10)) for j in range(r)]
    G = np.array([[0.25*np.trace(a @ b).real for b in Vimg] for a in Vimg])
    return np.sort(np.real(sla.eigvals((S + S.T)/2, G))), Fm


bew = {"Vakuum": bewegungen(P_OM), "Markenvakuum": bewegungen(P_MARK), "Kappenzustand s = 0,3": bewegungen(P_OM @ sla.expm(0.3*MQ0))}
mu_vak, F_vak = mu_induziert(NAHT, NAHT.OMEGA, P_OM)
mu_mark, F_mark = mu_induziert(NAHT, psi_mark, P_MARK)
mu_kap = {s: mu_induziert(NAHT, NAHT.gauss(P_OM @ sla.expm(s*MQ0)), P_OM @ sla.expm(s*MQ0))[0] for s in (0.05, 0.1, 0.2)}
NAHT_K = Naht(0.05)
NAHT_N = Naht(Q_NATIV)
mu_k = mu_induziert(NAHT_K, NAHT_K.gauss(ebene_z(0.05j)), ebene_z(0.05j))[0]
mu_n = mu_induziert(NAHT_N, NAHT_N.gauss(ebene_z(1j*Q_NATIV)), ebene_z(1j*Q_NATIV), cut=1e-30)[0]
vak_EQ = float(np.abs(F_vak[6:, :]).max())
ZAHLEN["q_gradient"] = {"Rang J / Spektrum / gemeinsam": {k: list(v) for k, v in bew.items()},
                        "BKM Vakuum auf E, Q": vak_EQ, "BKM Markenvakuum diag (L0..L5, E, Q1..Q3)": np.round(np.diag(F_mark), 6).tolist(),
                        "mu Vakuum": mu_vak.tolist(), "mu Markenvakuum (q = 0,25)": mu_mark.tolist(),
                        "mu Kappenzustand (s = 0,05; 0,1; 0,2)": {str(k): v.tolist() for k, v in mu_kap.items()},
                        "mu/q^2 Markenvakuum (q = 0,05; nativ)": [float(mu_k[-1]/0.05**2), float(mu_n[-1]/Q_NATIV**2)],
                        "mu Markenvakuum nativ": float(mu_n[-1])}
ov = {k: v[0] + v[1] - v[2] for k, v in bew.items()}
r_kap = float(mu_kap[0.1][-1]/mu_kap[0.05][-1])
check(f"[E] q-Gradient: Bei g = 1 bewegen im Vakuum 4 Richtungen J und 3 (die Boosts) das Spektrum, gemeinsam 7, "
      f"Überlappung {ov['Vakuum']}; beim Markenvakuum und beim Kappenzustand sind es 4 und 4, gemeinsam 7, Überlappung "
      f"{ov['Markenvakuum']} bzw. {ov['Kappenzustand s = 0,3']}: genau eine Richtung (die Uhr, SO(3)-invariant) bewegt J "
      f"und q zugleich. Quadratische Antwort (Fock, BKM): im Vakuum null auf E und Q ({vak_EQ:.0e}), beim Markenvakuum "
      f"nicht (E {F_mark[6, 6]:.4f}, Q {F_mark[7, 7]:.4f}). Induziert auf der Familie (Minimum über den Kern der "
      f"J-Bewegung, relativ zur Spurmetrik von U(2)): Vakuum Rang 0 (max |mu| {np.abs(mu_vak).max():.0e}), Markenvakuum "
      f"Rang 1 mit mu = {mu_mark[-1]:.6f} bei q = 0,25 und mu/q^2 = {mu_k[-1]/0.05**2:.6f} (q = 0,05), "
      f"{mu_n[-1]/Q_NATIV**2:.9f} (nativ, mu = {mu_n[-1]:.4e} = 16 e^(-4 pi) bis O(q^4)): fest, ohne freien Parameter; "
      f"Kappenzustand Rang 1 mit mu = {mu_kap[0.05][-1]:.5f}, {mu_kap[0.1][-1]:.5f}, {mu_kap[0.2][-1]:.5f} für "
      f"s = 0,05; 0,1; 0,2 (Verhältnis {r_kap:.3f} ~ 4): quadratisch in der freien Kappenzeit. Der q-Gradient ist eine "
      f"Form vom Rang 1 längs der Uhr auf einer Familie, deren eigene Krümmung null ist (2.3): Er kann keinen Koeffizienten "
      f"b einer vierdimensionalen Krümmung festlegen",
      bew["Vakuum"] == (4, 3, 7) and bew["Markenvakuum"] == (4, 4, 7) and bew["Kappenzustand s = 0,3"] == (4, 4, 7)
      and vak_EQ < 1e-10 and F_mark[6, 6] > 0.1 and F_mark[7, 7] > 0.1 and np.abs(mu_vak).max() < 1e-8
      and np.abs(mu_mark[:-1]).max() < 1e-5 and mu_mark[-1] > 0.5 and abs(mu_k[-1]/0.05**2 - 16) < 1e-2
      and abs(mu_n[-1]/Q_NATIV**2 - 16) < 1e-6 and all(np.abs(v[:-1]).max() < 1e-5 for v in mu_kap.values())
      and 3.9 < r_kap < 4.3)
laufzeit()

# --- 2.5 Komposita (K3)
say("  2.5 Komposita: Nähte an verschiedenen Orten bzw. mit q-Gefälle (diagonale Wirkung der nativen Gruppe)")


def tau_naht(R, Rk, P, q=Q_TEST):
    A_, B_ = ebenen(q)
    return tomita1(A_ @ Rk @ R, B_ @ Rk @ R, P)[0]


def komposit(R, naehte):
    return [tau_naht(R, Rk, P, q) for Rk, P, q in naehte]


def rang_komp(R0, naehte, h=1e-6):
    D = [np.concatenate([vec((a - b)/(2*h)) for a, b in zip(komposit(sla.expm(h*X) @ R0, naehte),
                                                             komposit(sla.expm(-h*X) @ R0, naehte))]) for X in SO32]
    _, sv, vh = np.linalg.svd(np.array(D).T, full_matrices=False)
    r = int(np.sum(sv > 1e-6*sv[0]))
    return r, vh[r:]


RE = lambda th: sla.expm(th*ME)
KOMP = {"antipodale Marken (0, pi)": [(np.eye(4), P_OM, Q_TEST), (RE(np.pi), P_OM, Q_TEST)],
        "Viertelmarken (0, pi/2)": [(np.eye(4), P_OM, Q_TEST), (RE(np.pi/2), P_OM, Q_TEST)],
        "vier Marken": [(RE(k*np.pi/2), P_OM, Q_TEST) for k in range(4)],
        "Vakuum + Kappenpartner (s = 0,3)": [(np.eye(4), P_OM, Q_TEST), (np.eye(4), P_OM @ sla.expm(0.3*MQ0), Q_TEST)],
        "q-Gefälle (q = 0,25 und 0,5, ein Ort)": [(np.eye(4), P_OM, 0.25), (np.eye(4), P_OM, 0.5)]}
komp = {}
for name_, nae in KOMP.items():
    r1, kern1 = rang_komp(np.eye(4), nae)
    rz = [rang_komp(zufall(0.5), nae)[0] for _ in range(2)]
    io_ = [sum(iota(t).real for t in komposit(zufall(0.5), nae)) for _ in range(5)]
    komp[name_] = {"Rang g = 1": r1, "Rang Zufall": rz, "iota_ges Spannweite": float(np.ptp(io_)),
                   "cos Kern gegen Drehungen (g = 1)": cos_kern_rot(kern1) if len(kern1) == 3 else None}
iso = 0.0
for _ in range(5):
    b1, b2 = (t[:2, 2:] for t in komposit(zufall(0.5), KOMP["antipodale Marken (0, pi)"]))
    iso = max(iso, float(np.abs(b2 + np.linalg.det(b1)*np.linalg.inv(b1)).max()))
# q-Gefälle auf dem globalen AdS-Schnitt g = e^(phi E) e^(eps R) e^(rho n.Q): Verzerrtes Produkt
n_ = np.array([0.3, -0.5, 0.81])
n_ /= np.linalg.norm(n_)
Xn = sum(n_[i]*QP[i + 1] for i in range(3))


def metrik_naht(rho, q, h=1e-5):
    fs = {"rho": lambda x: Rw([(Xn, rho + x)]), "eps": lambda x: Rw([(ROT[0], x), (Xn, rho)]),
          "phi": lambda x: Rw([(ME, x), (Xn, rho)])}
    vs = {k: (tau_von(f(h), q=q) - tau_von(f(-h), q=q))/(2*h) for k, f in fs.items()}
    return {(a, b): 0.25*np.trace(vs[a] @ vs[b]).real for a in vs for b in vs}


war_err, quer = 0.0, []
for qq in (0.25, 0.5):
    cq = (1 + qq*qq)/(2*qq)
    for rho in (0.3, 1.0, 2.0):
        Gm = metrik_naht(rho, qq)
        th_ = math.atan(cq*math.sinh(rho))
        war_err = max(war_err, abs(Gm[("rho", "rho")] - (cq*math.cosh(rho)/(1 + cq*cq*math.sinh(rho)**2))**2),
                      abs(Gm[("phi", "phi")] - 1), max(abs(Gm[(a, b)]) for a, b in (("rho", "eps"), ("rho", "phi"), ("eps", "phi"))))
        quer.append(Gm[("eps", "eps")]/math.sin(th_)**2)
_Bq = np.array([vec(M) for M in QP[1:]]).T                          # Drehvektor von ROT[0] aus Ad auf Q^+_1..3
_om = np.array([np.linalg.lstsq(_Bq, vec(kl(ROT[0], QP[k + 1])), rcond=None)[0] for k in range(3)])
rxn2 = float(np.linalg.norm(np.cross(np.array([_om[1, 2], _om[2, 0], _om[0, 1]]), n_))**2)
rho_s = sp.symbols("rho", positive=True)


def R3(cs):
    th = [sp.atan(c*sp.sinh(rho_s)) for c in cs]
    F = sum(sp.diff(t, rho_s)**2 for t in th)
    h = sp.sqrt(sum(sp.sin(t)**2 for t in th))
    hs = sp.diff(h, rho_s)/sp.sqrt(F)
    hss = sp.diff(hs, rho_s)/sp.sqrt(F)
    return -4*hss/h + 2*(1 - hs**2)/h**2


RHO_P = (0.05, 0.3, 1.0, 2.0, 4.0)
R_gleich = [float(R3((sp.Rational(17, 8), sp.Rational(17, 8))).subs(rho_s, r_)) for r_ in RHO_P]
R_gefaelle = [float(R3((sp.Rational(17, 8), sp.Rational(5, 4))).subs(rho_s, r_)) for r_ in RHO_P]
ZAHLEN["komposita"] = {"Familien": komp, "antipodal: |beta_2 + det(beta_1) beta_1^-1|": iso,
                       "q-Gefälle: Fehler des verzerrten Produkts": war_err, "Querfaktor (konstant)": [min(quer), max(quer)],
                       "|r x n|^2 des Drehgenerators (unabhängig)": rxn2,
                       "R_3(rho) gleiche q (rho = 0,05; 0,3; 1; 2; 4)": R_gleich, "R_3(rho) q-Gefälle": R_gefaelle}
check(f"[E] Komposita (K3, diagonale Wirkung): Zwei Nähte an antipodalen Marken (0, pi) geben Rang 4 "
      f"({komp['antipodale Marken (0, pi)']['Rang g = 1']}, {komp['antipodale Marken (0, pi)']['Rang Zufall']}) mit "
      f"beta_2 = -det(beta_1) beta_1^-1 ({iso:.0e}), dem Graphen einer Isometrie von U(2): homogen, flach. Nähte an "
      f"verschiedenen Orten entfalten sich: Viertelmarken, alle vier Marken und Vakuum + Kappenpartner haben Rang "
      f"{[komp[k]['Rang g = 1'] for k in ('Viertelmarken (0, pi/2)', 'vier Marken', 'Vakuum + Kappenpartner (s = 0,3)')]} "
      f"bei g = 1 und an Zufallspunkten = dim G/SO(3), Kern genau die Drehungen (cos >= "
      f"{min(komp[k]['cos Kern gegen Drehungen (g = 1)'] for k in ('Viertelmarken (0, pi/2)', 'vier Marken')):.9f}), "
      f"iota_ges variiert (Spannweite >= {min(v['iota_ges Spannweite'] for v in komp.values()):.3f}): weder Bahn der "
      f"Bewegungs- noch der Stabilisatorgruppe (Rang 7 < 4n), aber ohne vierdimensionales Modell; die Lage der Marken ist "
      f"nativ fest. Ein q-Gefälle zwischen zwei Nähten am selben Ort hat bei g = 1 Rang "
      f"{komp['q-Gefälle (q = 0,25 und 0,5, ein Ort)']['Rang g = 1']}, an Zufallspunkten "
      f"{komp['q-Gefälle (q = 0,25 und 0,5, ein Ort)']['Rang Zufall']}: es entfaltet sich ebenfalls. Auf dem globalen "
      f"AdS-Schnitt e^(phi E) e^(rho n.Q) ist jede Naht exakt das verzerrte Produkt dphi^2 + theta'^2 drho^2 + sin^2 theta "
      f"dOmega^2 mit tan theta = c sinh rho (Fehler {war_err:.0e}; Querfaktor konstant {np.mean(quer):.9f} = |r x n|^2 "
      f"des Drehgenerators {rxn2:.9f}, also Einheitssphäre); das "
      f"Komposit hat den Skalar R_3(rho) = {np.round(R_gleich, 6).tolist()} bei gleichem q (S^1 x S^3, konform flach) und "
      f"{np.round(R_gefaelle, 6).tolist()} bei q = 0,25 und 0,5: ortsabhängige Krümmung (S^1 x M^3 mit nicht konstanter "
      f"Krümmung, Weyl != 0), deren Größe das q-Gefälle ist. Nativ tragen alle Nähte dasselbe q = e^(-2 pi)",
      komp["antipodale Marken (0, pi)"]["Rang g = 1"] == 4 and komp["antipodale Marken (0, pi)"]["Rang Zufall"] == [4, 4]
      and iso < 1e-12 and all(komp[k]["Rang g = 1"] == 7 and komp[k]["Rang Zufall"] == [7, 7]
                              for k in ("Viertelmarken (0, pi/2)", "vier Marken", "Vakuum + Kappenpartner (s = 0,3)"))
      and min(komp[k]["cos Kern gegen Drehungen (g = 1)"] for k in ("Viertelmarken (0, pi/2)", "vier Marken")) > 1 - 1e-6
      and min(v["iota_ges Spannweite"] for v in komp.values()) > 0.1
      and komp["q-Gefälle (q = 0,25 und 0,5, ein Ort)"]["Rang g = 1"] == 4
      and komp["q-Gefälle (q = 0,25 und 0,5, ein Ort)"]["Rang Zufall"] == [7, 7] and war_err < 1e-8
      and np.ptp(quer) < 1e-8 and abs(np.mean(quer) - rxn2) < 1e-8 and np.ptp(R_gleich) < 1e-9 and abs(R_gleich[0] - 3) < 1e-9
      and np.ptp(R_gefaelle) > 0.5)
laufzeit()

# --- 2.6 Native Überlagerung (K4): Z4-Überlagerung der Markenvakua
say("  2.6 Z4-Überlagerung der Markenvakua (K4): native Amplitude, Standardform, Rang der Familie")
EBENEN_MARK = [ebene_z(Q_TEST*1j**k) for k in range(4)]
z4_eben = max(float(np.abs(z_von(EBENEN_MARK[k] @ sla.expm(np.pi/2*ME)) - z_von(EBENEN_MARK[(k + 1) % 4])).max())
              for k in range(4))
PSI_K = [NAHT.gauss(P) for P in EBENEN_MARK]                         # Phase <Omega|psi_k> > 0 = Uhr-Konvention (U Omega = Omega)
A_K = sum(np.vdot(NAHT.OMEGA, p) for p in PSI_K)


def katze_bei(R0, eps=1.0, N=NAHT, psis=PSI_K, ebenen_=EBENEN_MARK):
    """U(g0)^* psi_eps, psi_eps ~ (1 - eps) A Omega + eps sum_k psi_k (eps = 1: Z4-Überlagerung); jede Komponente ist
    gaußsch, U(g0)^* psi_j = e^(i chi_j) psi(P_j R0^-1), die relativen Phasen folgen aus der Unitarität:
    e^(i (chi_j - chi_Omega)) = <Omega|psi_j>/<psi(P_Omega R0^-1)|psi(P_j R0^-1)>."""
    Ri = np.linalg.inv(R0)
    prim_om = N.gauss(P_OM @ Ri)
    v = (1 - eps)*A_K*prim_om
    kons = 0.0
    for p, P in zip(psis, ebenen_):
        prim = N.gauss(P @ Ri)
        r = np.vdot(N.OMEGA, p)/np.vdot(prim_om, prim)
        kons = max(kons, abs(abs(r) - 1))
        v = v + eps*(r/abs(r))*prim
    return v/np.linalg.norm(v), kons


PSI_Z4 = katze_bei(np.eye(4))[0]
ueb = complex(np.vdot(NAHT.OMEGA, PSI_Z4))
dev = PSI_Z4 - ueb*NAHT.OMEGA
lam_nat = float(np.linalg.norm(dev)/abs(ueb))
lam_soll = math.sqrt(5)*Q_TEST**4
low6 = NAHT.TOT <= 6
z4_inv = float(np.abs((NAHT.uhr(np.pi/2, PSI_Z4) - PSI_Z4)[low6]).max())
z4_ko = float(np.abs((NAHT.uhr(np.pi/4, PSI_Z4) - PSI_Z4)[low6]).max())
C_z4 = PSI_Z4.reshape(NAHT.DA, NAHT.DA)
std_off = float(np.abs(C_z4 - np.diag(np.diag(C_z4))).max())
std_min = float(np.diag(C_z4).real.min())
std_im = float(np.abs(np.diag(C_z4).imag).max())
sv_1, vh_1 = tangenten(NAHT, PSI_Z4)
LOR_K = np.array([koord(M) for M in LOR])
cos_lor = float(np.linalg.svd(np.linalg.qr(vh_1[-6:].T)[0].T @ np.linalg.qr(LOR_K.T)[0], compute_uv=False).min())
R_TEST = [zufall(0.15) for _ in range(2)]
sv_z, kons_z = [], []
for R0 in R_TEST:
    v, k_ = katze_bei(R0)
    sv_z.append(tangenten(NAHT, v)[0])
    kons_z.append(k_)
sv_vak_g = tangenten(NAHT, NAHT.gauss(P_OM @ np.linalg.inv(R_TEST[0])))[0]
sv_e = {e: tangenten(NAHT, katze_bei(R_TEST[0], e)[0])[0] for e in (0.05, 0.1)}
lin = (sv_e[0.1][4:7]/sv_e[0.05][4:7]).tolist()
NAHT13 = Naht(Q_TEST, 13)
PSI_K13 = [NAHT13.gauss(P) for P in EBENEN_MARK]
_Ri = np.linalg.inv(R_TEST[0])
_po = NAHT13.gauss(P_OM @ _Ri)
_v13 = sum((lambda prim, p: (np.vdot(NAHT13.OMEGA, p)/np.vdot(_po, prim))/abs(np.vdot(NAHT13.OMEGA, p)/np.vdot(_po, prim))*prim)(
    NAHT13.gauss(P @ _Ri), p) for p, P in zip(PSI_K13, EBENEN_MARK))
sv_13 = tangenten(NAHT13, _v13/np.linalg.norm(_v13))[0]
konv = float(np.abs(sv_13[4:7] - sv_z[0][4:7]).max())
rk = lambda sv: int(np.sum(sv > 1e-6*sv[0]))
ZAHLEN["z4_katze"] = {"Vierteluhr permutiert die Markenebenen": z4_eben,
                      "lambda = |psi - <Omega|psi> Omega|/|<Omega|psi>|": lam_nat, "sqrt(5) q^4": lam_soll,
                      "1 - |<Omega|psi>|^2": float(1 - abs(ueb)**2), "5 q^8": 5*Q_TEST**8,
                      "Z4-Invarianz (Fock, tiefer Block)": z4_inv, "Kontrolle Achteldrehung": z4_ko,
                      "Standardform: Nebendiagonale, min Diagonale, Imaginärteil": [std_off, std_min, std_im],
                      "Singulärwerte g = 1": np.round(sv_1, 9).tolist(), "cos Kern gegen Lorentz (g = 1)": cos_lor,
                      "Singulärwerte an zwei Zufallspunkten": [np.round(s_, 6).tolist() for s_ in sv_z],
                      "Phasenkonsistenz | |r| - 1 |": max(kons_z), "Vakuum am selben Punkt": np.round(sv_vak_g, 9).tolist(),
                      "neue Richtungen eps = 0,05 / 0,1": [sv_e[0.05][4:7].tolist(), sv_e[0.1][4:7].tolist()],
                      "Linearität (0,1 / 0,05)": lin, "Konvergenz NV = 11 gegen 13": konv,
                      "native Amplitude sqrt(5) e^(-8 pi)": math.sqrt(5)*math.exp(-8*math.pi)}
check(f"[E] Z4-Überlagerung (K4): Die Vierteluhr permutiert die vier Markenebenen zyklisch ({z4_eben:.0e}); mit der "
      f"Uhrphase <Omega|psi_k> > 0 ist psi_Z4 ~ Summe_k psi_k (exakte Gaußvektoren) Z4-invariant (Fock, tiefer Block "
      f"{z4_inv:.0e}; Achteldrehung als Kontrolle {z4_ko:.0e}) und weicht vom Vakuum genau um den Vierpaarterm ab: "
      f"lambda = {lam_nat:.7f} gegen sqrt(5) q^4 = {lam_soll:.7f}, 1 - |<Omega|psi>|^2 = {1 - abs(ueb)**2:.5e} gegen "
      f"5 q^8 = {5*Q_TEST**8:.5e}. Bei g = 1 ist sie in Standardform (Koeffizientenmatrix diagonal {std_off:.0e}, reell "
      f"{std_im:.0e}, positiv >= {std_min:.1e}), also hält Lorentz J fest wie bei jedem Standardzustand (Rang "
      f"{rk(sv_1)}, Kern = Lorentz, cos >= {cos_lor:.9f}). An generischen Punkten (exakte Gaußkomponenten mit "
      f"Phasenverfolgung, Konsistenz {max(kons_z):.0e}) entfaltet sich die Familie: Rang {[rk(s_) for s_ in sv_z]} "
      f"(Vakuum am selben Punkt {rk(sv_vak_g)}), neue Singulärwerte {np.round(sv_z[0][4:7], 4).tolist()} (NV = 13: "
      f"Abweichung {konv:.0e}), linear in der Beimischung (Verhältnis {np.round(lin, 4).tolist()} bei Verdopplung von "
      f"eps): der (a2)-Mechanismus von z115, aber mit nativ fester Amplitude sqrt(5) q^4, bei nativem q sqrt(5) e^(-8 pi) "
      f"= {math.sqrt(5)*math.exp(-8*math.pi):.2e}, von der Ordnung der Kappen-Nichtinvarianz des Originals (z104). Fest, "
      f"aber siebendimensional: kein vierdimensionales Modell",
      z4_eben < 1e-12 and z4_inv < 1e-5 and z4_ko > 1e-3 and abs(lam_nat/lam_soll - 1) < 1e-3
      and abs((1 - abs(ueb)**2)/(5*Q_TEST**8) - 1) < 1e-3 and std_off < 1e-15 and std_im < 1e-15 and std_min > 0
      and rk(sv_1) == 4 and cos_lor > 1 - 1e-9 and [rk(s_) for s_ in sv_z] == [7, 7] and max(kons_z) < 1e-5
      and rk(sv_vak_g) == 4 and min(sv_z[0][4:7]) > 0.1 and konv < 1e-4 and all(abs(x - 2) < 2e-2 for x in lin))
befund("E (Teil 2): Kein natives Profil gibt eine gekrümmte vierdimensionale Familie mit festem Koeffizienten. Lokale "
       "Gleichgewichte (K1) sind q-blind und geben die flache Bahn. Jedes reine Gaußprofil (Markenvakuum, Kappenzustand, "
       "jedes Translat) liegt offen in Conj_psi = U(2), einer Stabilisatorbahn: flach, ohne nativen Rahmen; der Ort ändert "
       "nur die Karte und koppelt J an das Spektrum (q-Gradient längs der Uhr, BKM vom Rang 1: fest 16 e^(-4 pi) beim "
       "Markenvakuum, frei ~ s^2 beim Kappenzustand). Die Struktur dahinter: Conj_Omega = S^1_E x S^3_(Q^-), Zeit aus AdS, "
       "Raum aus dS; reduziert ist die Lötform im AdS-Modell vom Rang 1, im dS-Modell vom Rang 3. Familien, die weder Bahn "
       "der Bewegungs- noch der Stabilisatorgruppe sind, entstehen erst mit mehreren Orten oder einem q-Gefälle (Komposita) "
       "bzw. mit Nichtgaußizität (Z4-Überlagerung mit nativer Amplitude sqrt(5) e^(-8 pi)); sie entfalten sich auf Rang 7 "
       "= dim G/SO(3) und haben kein vierdimensionales Modell. Ortsabhängige Krümmung zeigt sich nur auf dem globalen "
       "AdS-Schnitt eines Komposits mit q-Gefälle, mit dem Gefälle als Koeffizient; nativ gibt es kein zweites q", True)
laufzeit()

# =================================================================== Teil 3
say("\nTeil 3: Die Kappe (Modulzeit P_cap, hyperbolische Klasse, dS)")

# --- 3.1 dS-Bahn von tau
say("  3.1 dS-Bahn von tau_0 unter SO(4,1) = Der(I_pi + I_omega)")
DER41 = teilalgebra([lambda M: M @ ISYMD + ISYMD @ M.T])
in41 = max(rest_im_spann(M, DER41) for M in SO41)
qp_raus = min(rest_im_spann(M, DER41) for M in QP)
_rows41 = [vec(kl(X, S4)) for X in SO41]
_sv41 = np.linalg.svd(np.array(_rows41).T, compute_uv=False)
stab41 = int(np.sum(_sv41 < 1e-9*_sv41[0]))
ad_tau = max(float(np.abs(S4 @ X @ S4 + X).max()) for X in QM)
lor_gerade = max(float(np.abs(S4 @ X @ S4 - X).max()) for X in LOR)


def krummung(Q):
    B = np.array([[np.trace(x @ y).real for y in Q] for x in Q])
    ks, res = [], 0.0
    for a, b, c in itertools.product(range(4), repeat=3):
        lhs = kl(kl(Q[a], Q[b]), Q[c])
        rhs = B[b, c]*Q[a] - B[a, c]*Q[b]
        if np.abs(rhs).max() > 1e-9:
            k = np.vdot(rhs.ravel(), lhs.ravel()).real/np.vdot(rhs.ravel(), rhs.ravel()).real
            ks.append(k)
            res = max(res, float(np.abs(lhs - k*rhs).max()))
    s_ = -1.0 if B[0, 0] > 0 else 1.0                                 # Signatur normiert: Zeit Q_0 negativ
    br = kl(Q[1], Q[2])
    D = np.trace(Q[1] @ Q[1]).real*np.trace(Q[2] @ Q[2]).real - np.trace(Q[1] @ Q[2]).real**2
    K_B = float(np.trace(br @ br).real/D)
    return np.round(np.diag(B), 9).tolist(), float(np.mean(ks)), res, K_B/s_


B_ads, k_ads, res_ads, K_ads = krummung(QP)
B_ds, k_ds, res_ds, K_ds = krummung(QM)
io41 = []
for _ in range(5):
    R_ = zufall(0.6, SO41)
    tz = tau_von(zufall(0.4))
    io41.append(abs(iota(np.linalg.inv(R_) @ tz @ R_, ISYMD, ISYMD_INV) - iota(tz, ISYMD, ISYMD_INV)))
io_bahn41 = max(abs(iota(np.linalg.inv(R_) @ S4 @ R_, ISYMD, ISYMD_INV) - 4) for R_ in [zufall(0.6, SO41) for _ in range(5)])
io_ads_auf_ds = float(np.ptp([iota(np.linalg.inv(R_) @ S4 @ R_).real for R_ in [zufall(0.6, SO41) for _ in range(5)]]))
XLd, YLd = QM[1] + 0.3*QM[0] + 0.2*LOR[2], QM[2] - 0.5*QM[3] + 0.4*LOR[4]
sch_ds = schleife(lambda s_, t_: Rw([(XLd, s_), (YLd, t_)]))
ZAHLEN["ds_bahn"] = {"SO41 in Der(I_pi + I_omega)": in41, "Kontrolle Q^+ nicht darin": qp_raus, "Stabilisator von tau_0": stab41,
                     "Ad tau = -1 auf Q^-": ad_tau, "Lorentz gerade": lor_gerade, "B|Q^+, k, Rest, K": [B_ads, k_ads, res_ads, K_ads],
                     "B|Q^-, k, Rest, K": [B_ds, k_ds, res_ds, K_ds], "iota_dS invariant": float(max(io41)),
                     "iota_dS - 4 auf der dS-Bahn": float(io_bahn41), "Kontrolle: iota_AdS auf der dS-Bahn (Spannweite)": io_ads_auf_ds,
                     "Schleife dS-Bahn (Rahmen in SO(4,1))": sch_ds}
check(f"[E/C] dS-Bahn: span{{Lorentz, Q^-_mu}} liegt in Der(I_pi + I_omega) (Rest {in41:.0e}; Kontrolle Q^+: "
      f"{qp_raus:.3f}), der Stabilisator von tau_0 in so(4,1) hat Dimension {stab41} (Lorentz), Ad tau_0 = -1 auf Q^- "
      f"({ad_tau:.0e}) und +1 auf Lorentz ({lor_gerade:.0e}): das symmetrische Paar von dS4. [[Q_a, Q_b], Q_c] = "
      f"k (B_bc Q_a - B_ac Q_b) mit k = {k_ds:.6f} (Rest {res_ds:.0e}) wie bei AdS (k = {k_ads:.6f}), aber mit "
      f"umgekehrter Spurform B|Q^- = {B_ds} gegen B|Q^+ = {B_ads}: bei gleich normierter Signatur (Zeit negativ) ist die "
      f"Schnittkrümmung {K_ds:+.6f} gegen {K_ads:+.6f}, also K_phys = {4*K_ds:+.1f} gegen {4*K_ads:+.1f} mit l = 1/2: "
      f"Lambda = +12/R_Vak^2 (dS) gegen -12/R_Vak^2 (AdS), das Vorzeichen ist die Klasse der Zeit (z104) schon auf der "
      f"Tomita-Bahn. iota_dS = tr(tau I' tau^T I'^-1), I' = I_pi + I_omega, ist SO(4,1)-invariant ({max(io41):.0e}) und auf "
      f"der dS-Bahn 4 ({io_bahn41:.0e}); iota_AdS variiert dort (Spannweite {io_ads_auf_ds:.3f}). Die Bahnfamilie mit "
      f"Rahmen in SO(4,1) ist Maurer-Cartan ({sch_ds['|F_MC|']:.0e}), R = [e, e] (Rest {sch_ds['|R - [e, e]|']:.0e}), "
      f"Torsion {sch_ds['|T|']:.0e}: das flache MacDowell-Mansouri-Vakuum der dS-Klasse. Lambda ist unter P_cap fest "
      f"(l = R_Vak/2), b = -R_Vak^2/(256 pi G) nicht (G offen, z104)",
      in41 < 1e-9 and qp_raus > 0.5 and stab41 == 6 and ad_tau < 1e-12 and lor_gerade < 1e-12 and abs(k_ds - 1) < 1e-9
      and abs(k_ads - 1) < 1e-9 and res_ds < 1e-12 and np.allclose(B_ds, -np.array(B_ads)) and K_ds > 0 > K_ads
      and abs(K_ds + K_ads) < 1e-9 and max(io41) < 1e-12 and io_bahn41 < 1e-10 and io_ads_auf_ds > 0.1
      and sch_ds["|F_MC|"] < 1e-5 and sch_ds["|R - [e, e]|"] < 1e-5 and sch_ds["|T|"] < 1e-5)
laufzeit()

# --- 3.2 Naht gegen Kappenfluss
say("  3.2 Kein Nahtmodulfluss ist hyperbolisch: die Kappe braucht das Typ-III-Netz")
Xr = sum(c*X for c, X in zip(rng.normal(size=16), U22))*0.15
psi_r = spla.expm_multiply(1j*NAHT.gen(Xr), NAHT.OMEGA)
Cr = psi_r.reshape(NAHT.DA, NAHT.DA)
pr_, Ur_ = np.linalg.eigh(Cr @ Cr.conj().T)
kr = pr_ > 1e-12
Ka_r = -(Ur_[:, kr]*np.log(pr_[kr])) @ Ur_[:, kr].conj().T
L_r, res_r, kreuz_r = ein_teilchen_L(Ka_r)
ew_L = np.linalg.eigvals(L_r)
herm_r = float(np.abs(L_r - L_r.conj().T).max())
ew_kap = np.linalg.eigvals(sla.expm(1.0*MQ0))
ew_uhr = np.linalg.eigvals(sla.expm(1.0*ME))
ZAHLEN["naht_gegen_kappe"] = {"zufällige Naht: L": np.round(L_r, 6).tolist(), "|L - L^+|": herm_r, "Eigenwerte L": np.round(ew_L, 6).tolist(),
                              "Fitrest, a^+-Anteil": [res_r, kreuz_r], "|Eigenwerte| von e^(Q^-)": np.round(np.abs(ew_kap), 9).tolist(),
                              "|Eigenwerte| von e^(E)": np.round(np.abs(ew_uhr), 9).tolist()}
check(f"[E] Naht gegen Kappenfluss: Für jede Naht (Typ-I-Faktor, gaußscher treuer Zustand; hier das Vakuum und ein "
      f"Zufallszustand U(h) Omega) ist der Modulhamiltonoperator K_a = -log rho_a eine ladungserhaltende Quadratik: "
      f"[K_a, a_i] = -sum L_ij a_j mit L hermitesch ({herm_r:.0e}; Fitrest {res_r:.0e}, a^+-Anteil {kreuz_r:.0e}; "
      f"Eigenwerte {np.round(ew_L.real, 6).tolist()}), also wirkt Delta^(it) auf T unimodular (elliptisch, kompakt). Der "
      f"Kappenfluss e^(s Q^-) hat auf T die Eigenwerte e^(+-s/2) (|.| = {np.round(np.abs(ew_kap), 6).tolist()} bei s = 1; "
      f"die Uhr: {np.round(np.abs(ew_uhr), 6).tolist()}): hyperbolisch. Keine Naht hat die Modulzeit P_cap; die Kappe "
      f"ist ein Typ-III-Faktor des Netzes (z104, bedingt auf das freie Netz). Dort ist das Vakuum unter der konformen "
      f"Gruppe invariant, die Kappenfamilie J(A(g O_1), Omega_4D) also eine Bahn der Zustandsstabilisatorgruppe und nach "
      f"z100 flach [C/L]",
      herm_r < 1e-10 and res_r < 1e-6 and kreuz_r < 1e-9 and np.abs(ew_L.imag).max() < 1e-9
      and np.abs(np.abs(ew_uhr) - 1).max() < 1e-12 and np.abs(np.abs(ew_kap) - 1).min() > 0.3
      and np.allclose(sorted(np.abs(ew_kap)), [math.exp(-0.5)]*2 + [math.exp(0.5)]*2))
laufzeit()

# --- 3.3 Feste Familie in der dS-Klasse
say("  3.3 Feste Familie in der dS-Klasse: SO(4) hält Omega fest, die Kappenzeit wird zum Uhrbogen")
so4_fix = max(float(np.abs((P_OM @ sla.expm(0.7*X))[:, 2:]).max()) for X in QM[1:] + ROT)
so4_bahn = max(float(np.abs(tau_von(sla.expm(t*X)) - np.linalg.inv(sla.expm(t*X)) @ S4 @ sla.expm(t*X)).max())
               for X in QM[1:] for t in (0.3, 1.1))
kart = {}
for qq in (Q_TEST, Q_NATIV):
    cq = (1 + qq*qq)/(2*qq)
    err, skal, ia, ids = 0.0, 0.0, [], []
    for s in (0.1, 0.5, 1.0, 2.0, 4.0):
        t_ = tau_von(sla.expm(s*MQ0), q=qq)
        be = t_[:2, 2:]
        skal = max(skal, float(np.abs(be - be[0, 0]*I2).max()))
        err = max(err, abs(abs(np.angle(be[0, 0])) - math.atan(cq*math.sinh(s))))
        ia.append(abs(iota(t_) - 4))
        ids.append(iota(t_, ISYMD, ISYMD_INV).real)
    kart[qq] = {"c": cq, "|phi| - atan(c sinh s)": float(err), "beta skalar": skal, "|iota_AdS - 4|": float(max(ia)),
                "iota_dS (s = 0,1 .. 4)": np.round(ids, 6).tolist()}
raenge41 = [rang(np.eye(4), basis=SO41)[0]] + [rang(zufall(0.5, SO41), basis=SO41)[0] for _ in range(3)]
kern41 = max(float(np.abs(dtau(np.eye(4), X)).max()) for X in LOR)
HS41 = [NAHT.gen(M) for M in SO41]
F41 = NAHT.bkm(NAHT.OMEGA, HS41)
ew41 = np.linalg.eigvalsh(F41)
F41_Q = float(np.abs(F41[6:, :]).max())


def phi_von(s, qq):
    return abs(np.angle(tau_von(sla.expm(s*MQ0), q=qq)[0, 2]))


konf = {}
for qq in (0.999, 0.9, Q_TEST, Q_NATIV):
    cq = (1 + qq*qq)/(2*qq)
    ks, soll = [], []
    for s in (0.3, 1.0, 2.0):
        dphi = (phi_von(s + 1e-5, qq) - phi_von(s - 1e-5, qq))/2e-5
        ks.append(dphi*math.cosh(s))
        soll.append(cq*math.cosh(s)**2/(1 + cq*cq*math.sinh(s)**2))
    konf[qq] = {"kappa_dS (s = 0,3; 1; 2)": ks, "Formel c cosh^2 s/(1 + c^2 sinh^2 s)": soll,
                "Defekt max |kappa - 1|": max(abs(k_ - 1) for k_ in ks), "Formelfehler": max(abs(a_ - b_) for a_, b_ in zip(ks, soll))}
ZAHLEN["ds_feste_familie"] = {"SO(4) hält Omega fest": so4_fix, "längs SO(4): fest = Bahn": so4_bahn,
                              "Kappenkarte": {str(k): v for k, v in kart.items()}, "Rang über so(4,1)": raenge41,
                              "max |dtau| auf Lorentz": kern41, "BKM-Eigenwerte (so(4,1))": np.round(ew41, 9).tolist(),
                              "BKM auf Q^-_mu": F41_Q, "Konformität": {str(k): v for k, v in konf.items()}}
check(f"[E/L] Feste Familie in der dS-Klasse: SO(4) = so(3) + <Q^-_k> hält Omega fest ({so4_fix:.0e}), längs SO(4) ist "
      f"J(M_g, Omega) die Bahnfamilie ({so4_bahn:.0e}): die Raumkugel S^3 wird exakt abgebildet. Längs der Kappe ist beta = "
      f"e^(i phi) 1 (skalar, {max(v['beta skalar'] for v in kart.values()):.0e}) mit tan|phi| = c sinh s, c = coth(2r) "
      f"(Fehler {max(v['|phi| - atan(c sinh s)'] for v in kart.values()):.0e} bei q = 0,25 und nativ, c = "
      f"{kart[Q_NATIV]['c']:.4f} = cosh 2 pi): dieselbe Kartenfunktion wie der AdS-Raum in z115, jetzt auf dem Zeitkreis; "
      f"|phi| < pi/2 für alle s. Dabei läuft die Familie auf dem AdS-Uhrkreis (|iota_AdS - 4| "
      f"{max(v['|iota_AdS - 4|'] for v in kart.values()):.0e}), während iota_dS variiert ({kart[Q_TEST]['iota_dS (s = 0,1 .. 4)']}): "
      f"kein nativer dS-Rahmen. Rang über so(4,1) {raenge41}, Kern bei g = 1 Lorentz ({kern41:.0e}); die quadratische "
      f"Antwort lebt nur auf den Boosts (Eigenwerte {np.round(ew41[ew41 > 1e-9], 6).tolist()}) und ist auf Q^-_mu null "
      f"({F41_Q:.0e}). Skalenverhältnis kappa_dS = (dphi/ds) cosh s = c cosh^2 s/(1 + c^2 sinh^2 s) (Formelfehler "
      f"{max(v['Formelfehler'] for v in konf.values()):.0e}), konform genau für c = 1: Defekt "
      f"{konf[0.999]['Defekt max |kappa - 1|']:.1e} bei q = 0,999, {konf[0.9]['Defekt max |kappa - 1|']:.3f} bei 0,9, "
      f"{konf[Q_TEST]['Defekt max |kappa - 1|']:.3f} bei 0,25, {konf[Q_NATIV]['Defekt max |kappa - 1|']:.4f} nativ "
      f"(kappa ~ coth^2 s/c). Im "
      f"Grenzfall q -> 1 ist die Familie die konforme Kompaktifizierung von dS4 im Einsteinstreifen |phi| < pi/2 "
      f"(Penrose-Quadrat), Spiegel von z115s AdS-Halbzylinder: AdS halbiert den Raum, dS die Zeit [L]",
      so4_fix < 1e-12 and so4_bahn < 1e-12 and max(v["beta skalar"] for v in kart.values()) < 1e-12
      and max(v["|phi| - atan(c sinh s)"] for v in kart.values()) < 1e-9 and max(v["|iota_AdS - 4|"] for v in kart.values()) < 1e-9
      and np.ptp(kart[Q_TEST]["iota_dS (s = 0,1 .. 4)"]) > 1 and raenge41 == [4, 4, 4, 4] and kern41 < 1e-8
      and int(np.sum(ew41 > 1e-9*ew41.max())) == 3 and np.allclose(ew41[ew41 > 1e-9*ew41.max()], 2, atol=1e-6)
      and F41_Q < 1e-10 and max(v["Formelfehler"] for v in konf.values()) < 1e-5
      and konf[0.999]["Defekt max |kappa - 1|"] < 1e-3 and konf[Q_TEST]["Defekt max |kappa - 1|"] > 0.1
      and konf[Q_NATIV]["Defekt max |kappa - 1|"] > 0.9)
befund("E/C (Teil 3): Relativ zur Kappe entsteht keine dS-Krümmung. Die Tomita-Bahn unter SO(4,1) ist dS4 mit Lambda = "
       "+12/R_Vak^2 unter P_cap, flach als Cartan-Geometrie (F = 0); b bleibt frei (G offen). Keine Naht trägt den "
       "hyperbolischen Kappenfluss; im Typ-III-Netz ist das Vakuum konform invariant und die Kappenfamilie eine Bahn. Die "
       "feste Familie der Naht in der dS-Klasse ist dieselbe kompakte Duale U(2): SO(4) hält Omega fest (S^3 exakt), die "
       "Kappenzeit wird mit tan phi = c sinh s auf den halben Uhrkreis gestaucht; für q -> 1 ist das die konforme "
       "Kompaktifizierung von dS4. Das Vorzeichen von Lambda entscheidet, welcher Faktor des Einsteinzylinders halbiert "
       "wird (Raum bei AdS, Zeit bei dS), nicht eine Krümmung", True)
laufzeit()

# =================================================================== Bilanz und Urteil
say("\nBilanz")
STATUS = {
    "Zustandsscheibe": ("E", "SO(3)-invariante Gaußzustände = Scheibe; radiale sl(2,R) = Möbiusgruppe; Fixpunkte von Q^- "
                             "und D = Marken z^4 = 1; Schmidtwert = pseudohyperbolischer Abstand"),
    "native Profile": ("E/C", "Wiesbrock-Profil q(a) und Markenvakua nativ (q(phi) mit Horizont); Hislop-Longo T(r) "
                              "bedingt auf das Netz; Kappenzeit s frei"),
    "q-Blindheit": ("E", "J(M_q, Omega_q) = J für jedes q: lokales Gleichgewicht = native Bahn (flach); Profil nur im "
                         "zentralen Modulgenerator"),
    "ESU-Struktur": ("E", "Conj_Omega = S^1_E x S^3_(Q^-); dtau(Q^+_k) = c [Q^-_k, tau_0]; Lötform AdS Rang 1, dS Rang 3"),
    "Gaußsches Lagelemma": ("E, negativ", "jedes reine Gaußprofil offen in Conj_psi = U(2), Stabilisatorbahn: flach, "
                                         "kein nativer Rahmen; Ort ändert nur die Karte"),
    "q-Gradient": ("E", "eine Richtung (Uhr) bewegt J und q; induzierte BKM Rang 1: 16 e^(-4 pi) fest (Markenvakuum), "
                        "~ s^2 frei (Kappenzustand); keine Krümmung zu normieren"),
    "Komposita": ("E, negativ", "verschiedene Orte bzw. q-Gefälle: Rang 7, kein 4D-Modell; antipodal: Isometriegraph; "
                                "q-Gefälle auf dem AdS-Schnitt: R_3(rho) ortsabhängig, Größe = Gefälle (nicht nativ)"),
    "Z4-Überlagerung": ("E", "native Amplitude sqrt(5) q^4 = sqrt(5) e^(-8 pi); bei g = 1 Standardform (Rang 4), "
                             "generisch Rang 7, kein 4D-Modell"),
    "dS-Bahn": ("E/C", "so(4,1) = Der(I_pi + I_omega), dS4, Lambda = +12/R_Vak^2 unter P_cap, F = 0, b frei"),
    "Naht gegen Kappe": ("E", "Nahtmodulfluss elliptisch, Kappe hyperbolisch: P_cap braucht das Typ-III-Netz"),
    "feste Familie dS": ("E/L", "SO(4) fixiert Omega; Kappenzeit tan phi = c sinh s auf dem Uhrbogen; q -> 1: dS4 im "
                                "Einsteinstreifen; Lambda-Vorzeichen halbiert Zeit statt Raum"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
urteil = ("Der Folgeauftrag zu T6 ist beantwortet: gescoptes Negativ für native Quellen, mit neuer exakter Struktur und "
          "einem benannten bedingten Mechanismus. (1) Native inhomogene Zustände gibt es. Die SO(3)-invarianten reinen "
          "Gaußzustände bilden eine Scheibe, auf der die radiale sl(2,R) als Möbiusgruppe wirkt: Die Uhr E dreht um das "
          "Vakuum, die Kappe Q^- und die Dilatation D verschieben längs senkrechter Durchmesser, deren Randfixpunkte "
          "zusammen genau die Marken z^4 = 1 sind; die Algebra M sitzt bei z = q, und der Schmidtwert eines Paares ist der "
          "pseudohyperbolische Abstand. Nativ ohne freien Parameter sind das Wiesbrock-Profil q(a)^2 = (q^2 + x^2 (1 - q)^2)/"
          "(1 + x^2 (1 - q)^2) der mit P_0 verschobenen Algebren und die Markenvakua (Fockvakua der Markenalgebren) mit dem "
          "Uhrprofil q(phi) und Horizont an der eigenen Marke; das Hislop-Longo-Profil T(r) = 1/(pi (1 - r^2)) ist nativ "
          "bedingt auf das Netz, seine Übersetzung in ein lokales q eine Wahl; die Kappenzeit eines verschobenen Zustands ist "
          "frei. (2) Keines dieser Profile krümmt die Familie J(M_g, psi). J ist q-blind: Ein lokales Gleichgewicht mit "
          "beliebigem Temperaturprofil hat die native Bahnfamilie, flach nach z100; das Profil sitzt nur im zentralen "
          "Modulgenerator. Jedes reine Gaußprofil gibt eine Familie, die offen in Conj_psi = U(2) liegt (Rang 4), einer "
          "einzigen Bahn der Stabilisatorgruppe: Das Lagelemma tau(M_g, U(h) Omega) = R_h^-1 tau(M_(h^-1 g), Omega) R_h "
          "zeigt, dass der Ort nur die Karte verschiebt. Struktur: Conj_Omega ist der kompaktifizierte Einsteinzylinder "
          "(S^1_E x S^3_(Q^-))/Z_2, "
          "Zeit aus der AdS-Uhr, Raum aus den dS-Raumtranslationen; die Algebra-Translationen Q^+_k werden mit dem Faktor "
          "coth(2r) auf Q^-_k abgebildet. Reduziert auf das AdS-Modell hat die Familie eine Lötform vom Rang 1 und als "
          "Krümmung nur die S^3-Berry-Krümmung, auf das dS-Modell Rang 3 und F = 0: in keinem Modell eine Cartan-Geometrie; "
          "iota variiert, ein nativer Rahmen fehlt. Der q-Gradient ist real, aber schwach: Bei inhomogenen Gaußzuständen "
          "bewegt genau eine Richtung (die Uhr) J und das Spektrum zugleich, die quadratische Antwort induziert auf der "
          "Familie eine Form vom Rang 1, fest 16 e^(-4 pi) beim Markenvakuum und ~ s^2 frei beim Kappenzustand; sie hat "
          "keine Krümmung zu normieren. (3) Familien, die weder Bahn der Bewegungs- noch der Stabilisatorgruppe sind, "
          "entstehen erst mit mehreren Orten (Komposita an Viertelmarken, mit Kappenpartner), mit einem q-Gefälle zwischen "
          "Nähten oder mit Nichtgaußizität (Z4-Überlagerung der Markenvakua, native Amplitude sqrt(5) e^(-8 pi), von der "
          "Ordnung der Kappen-Nichtinvarianz in z104). Alle entfalten sich auf Rang 7 = dim G/SO(3) und haben kein "
          "vierdimensionales Modell, F_MM im AdS-Modell bleibt undefiniert. Ortsabhängige Krümmung zeigt sich nur auf dem "
          "globalen AdS-Schnitt eines Komposits mit q-Gefälle (S^1 x M^3, R_3(rho) variiert, Weyl != 0), mit dem Gefälle "
          "als Koeffizient; nativ tragen alle Nähte dasselbe q, und gleiche Nähte geben S^1 x S^3, konform flach. (4) Kappe: "
          "Unter P_cap ist die Tomita-Bahn unter SO(4,1) = Der(I_pi + I_omega) dS4 mit Lambda = +12/R_Vak^2, flach (F = 0); "
          "die Schnittkrümmung kehrt gegenüber AdS das Vorzeichen um, b bleibt frei. Keine Naht trägt den hyperbolischen "
          "Kappenfluss (Nahtmodulfluss elliptisch), die Kappe ist ein Typ-III-Faktor des Netzes, wo das Vakuum konform "
          "invariant und die Kappenfamilie eine Bahn ist. Die feste Familie der Naht in der dS-Klasse ist dieselbe U(2): "
          "SO(4) hält Omega fest, die Raumkugel wird exakt abgebildet, die Kappenzeit mit tan phi = coth(2r) sinh s auf den "
          "halben Uhrkreis; für q -> 1 ist das die konforme Kompaktifizierung von dS4 im Einsteinstreifen, Spiegel der "
          "AdS-Halbzylinder-Karte von z115. Das Lambda-Vorzeichen entscheidet, welcher Faktor halbiert wird, nicht eine "
          "Krümmung. Urteil: gescoptes Negativ. Ein natives Zustandsprofil liefert keine Quelle gravitativer Dynamik mit "
          "festem b; bedingt (Prämisse P_dq: Nähte mit verschiedenem q am selben Ort) entsteht ortsabhängige euklidische "
          "Krümmung, aber nur auf einem gewählten Schnitt, ohne vierdimensionales Modell der vollen Familie und mit dem "
          "q-Gefälle als freiem Koeffizienten. Folgen: K13 bleibt L, präzisiert: Die Familie J(M_g, psi) bei festem Zustand "
          "liegt für jeden Zustand in Conj_psi, einer Bahn der Stabilisatorgruppe (gaußsch geprüft; in jeder endlichen "
          "Abschneidung allgemein, weil je zwei antiunitäre Involutionen mit festem psi durch eine psi festhaltende Unitäre "
          "konjugiert sind), mit Rahmen dort ist ihre Verbindung stets Maurer-Cartan. Krümmung kann nur als induzierte "
          "Geometrie einer echten Teilfamilie entstehen, und jede solche Familie, die hier vorkommt (Komposita, "
          "Nichtgaußizität), entfaltet sich auf Rang 7. Ein Ortsprofil allein krümmt also keine vierdimensionale Familie; "
          "dafür braucht es eine Dynamik, die Algebra und Zustand nicht durch Gruppenelemente verbindet, ein zweites q gibt "
          "Krümmung nur auf einem gewählten Schnitt. Für die Gravitationsachse der Szenario-Suche: Der Wert 'Krümmung aus "
          "inhomogenem Profil' "
          "ist kein gestützter Gravitationswert; 'q-Gefälle zwischen Nähten' ist höchstens ein bedingter Wert mit der "
          "Prämisse P_dq und ohne 4D-Modell. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
say(f"\n   (Laufzeit {time.time() - T0:.1f} s)")
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z128_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n", encoding="utf-8")
(HERE / "ausgabe_z128.txt").write_text("\n".join(OUT) + "\n", encoding="utf-8")
