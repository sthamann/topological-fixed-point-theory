"""Z104: Lambda > 0 aus der modularen Zeit einer Kappe? (Folgeauftrag N5, Hypothese H5).

Leitlinie 'Rekonstruktion statt Auswahl'. Firewalls: GRAV.NONCIRCULAR.01 (Richtung Zustand -> Zeit: zuerst das freie konforme
Netz des Leiterfelds aus z83/z90 und sein Vakuum, dann der Satz von Hislop-Longo; der Kappenfluss wird nirgends als Zeit
vorausgesetzt, sondern als Modulfluss eines Paares (Algebra, Zustand) geprüft), DIMENSION.UPLIFT.FIREWALL.01 (Kinematik einer
Darstellung und eines freien Felds, keine Dynamik). Gelöste Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt;
COSMO.LAM.01 bzw. v60 werden nur gelesen; nichts wird befördert. Bestand: z96 (Q_0^dS = Modulboost des Doppelkegels mit den
Spitzen t = +-1, Hislop-Longo), z98 (natives sigma = +1, l = R_Vak/2, Lambda = -12/R_Vak^2, Q_0^- spektralsymmetrisch), z90
(Modulgenerator der Kappe = hyperbolisches Element der radialen sl(2,R), s-Welle D_1), z83 (E = L_0, radiale sl(2,R),
Produktlesart Original x Partner bei beta' = 8 pi), z72 (parabolische Zeit P_0 = L_0 + H_1 im Vakuumsektor), z70/z84
(Modulfluss des Originals periodisch, die OS-Rekonstruktion gibt ihn zurück). Hypothese H5: Die physikalische (Hubble-)Zeit
ist die Modulzeit einer Kappe, nicht die radiale Energie; dann ist die MacDowell-Mansouri-Algebra so(4,1) und Lambda > 0.
Statuskürzel: E exakt, C bedingt, L Lesart, O offen.

Teil 1  Kappe und Modulzeit in der nativen radialen su(2,2) (Fock 8^4 und Twistorraum wie z98): K_cap = 2 pi Q_0^-, die
        radiale sl(2,R) {E, Q^-, D} als Original-Partner-Quetschung, Eigenebenen (definit = komplexe Punkte, null = Spitzen),
        Killinghorizont; native Kappe aus E, tau und den Marken; Modulfluss des Vakuums (KMS auf der Kappenachse) gegen den
        des Originals (zentral, periodisch, Zweipunktfunktion nicht stationär); exakte Beziehungen zu L_0, beta = 4 pi, P_0.
Teil 2  dS-Kinematik auf der radialen Leiter: Casimir von so(4,1) und so(3,2), K-Typen, Serie; Spektrum von Q^- je
        Drehimpuls (D_(l+1), sech-Gesetz, symmetrisch, stetig); Lambda als Funktion der Zeitwahl Q_0 = a E + b Q^- und
        Positivität; Gibbons-Hawking-Temperatur in Nahteinheiten; Beziehung zu beta = 4 pi (Erzeugerklasse, Perioden,
        euklidischer Schnitt).
Teil 3  Urteil H5: Auswahltabelle der Zeitkandidaten, Prämisse, Folgen für K13 und das Lambda-G-Ziel.
"""
import itertools
import json
import time
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.linalg as sla
import scipy.sparse as sps
import scipy.sparse.linalg as spla
import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
rng = np.random.default_rng(104)
mp.mp.dps = 40
I2 = np.eye(2)
SIGMA = [np.eye(2, dtype=complex), np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]),
         np.diag([1., -1.]).astype(complex)]
ALPHA_INV = mp.mpf("137.035999177")                                     # CODATA 2022, nur zur Größenordnung
BETA_P = 8*mp.pi                                                        # Produktlesart Original x Partner (z83)


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


# ------------------------------------------------------------------ Fockraum einer Naht (wie z98)
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


NV = 8
DIMV = NV**4
d1, d2, e1, e2 = fock4(NV)
DD, EE = [d1, d2], [e1, e2]
OCC = np.array(list(itertools.product(range(NV), repeat=4)))
DVAL = OCC[:, 0] + OCC[:, 1] - OCC[:, 2] - OCC[:, 3]
LOW = np.where(OCC.sum(1) <= NV - 3)[0]
LOW4 = np.where(OCC.max(1) <= 4)[0]
SC = np.where((DVAL == 0) & (OCC.sum(1) <= NV - 3))[0]                  # Skalarsektor N_d = N_e, E = 1, 2, 3
OMV = np.zeros(DIMV, complex)
OMV[0] = 1


def klein(A, tol=1e-10, spalten=None):
    A = sps.csc_matrix(A)[:, LOW if spalten is None else spalten]
    return (abs(A).max() if A.nnz else 0.0) < tol


ETA = np.diag([1., -1., -1., -1.])
kombis = list(itertools.combinations(range(4), 2))
ZT = [d1, d2, dag(e1), dag(e2)]                                         # Z^alpha = (d_1, d_2, e_1^+, e_2^+)
CC = [(DD[i] + dag(EE[i])).tocsr() for i in range(2)]                   # c = d + e^+ (pi-Ebene)
WO = [(DD[i] - dag(EE[i])).tocsr() for i in range(2)]                   # omega'' = d - e^+


def bil(X, A, Y):
    """Summe_ij A_ij X_i^+ Y_j (Bilinearform, exakt bis auf den Abschneiderand)."""
    return sum((A[i, j]*(dag(X[i]) @ Y[j]) for i in range(2) for j in range(2) if abs(A[i, j]) > 0),
               sps.csr_matrix((DIMV, DIMV))).tocsr()


Pop = [0.5*bil(CC, s, CC) for s in SIGMA]
Kop = [0.5*bil(WO, s, WO) for s in SIGMA]
KopL = [ETA[m, m]*Kop[m] for m in range(4)]                             # K_a = eta_ab K^b
Eop = (0.5*sum((dag(X) @ X for X in DD + EE), sps.csr_matrix((DIMV, DIMV))) + sps.identity(DIMV)).tocsr()
Dop = sps.diags(DVAL.astype(float)).tocsr()
# [P_m, K^n] = (1/2)(c^+ s_m s_n omega - omega^+ s_n s_m c), exakt als Bilinearform; M_mn hermitesch
PK = {(m, n): 0.5*(bil(CC, SIGMA[m] @ SIGMA[n], WO) - bil(WO, SIGMA[n] @ SIGMA[m], CC)) for m in range(4) for n in range(4)}
MopL = {(m, n): (-0.25j*(ETA[n, n]*PK[(m, n)] - ETA[m, m]*PK[(n, m)])).tocsr() for m, n in kombis}

# ------------------------------------------------------------------ u(2,2) auf dem Twistorraum (wie z98)
G4 = np.diag([1., 1., -1., -1.]).astype(complex)
S4 = np.block([[np.zeros((2, 2)), I2], [I2, np.zeros((2, 2))]]).astype(complex)   # tau = Theta_1 auf Z
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


def im_spann(M, basis):
    if not basis:
        return np.abs(M).max() < 1e-12
    A = np.array([np.concatenate([B.real.ravel(), B.imag.ravel()]) for B in basis]).T
    b = np.concatenate([M.real.ravel(), M.imag.ravel()])
    x = np.linalg.lstsq(A, b, rcond=None)[0]
    return np.linalg.norm(A @ x - b) < 1e-9*max(1.0, np.linalg.norm(b))


def koeff(M, basis):
    A = np.array([np.concatenate([B.real.ravel(), B.imag.ravel()]) for B in basis]).T
    return np.linalg.lstsq(A, np.concatenate([M.real.ravel(), M.imag.ravel()]), rcond=None)[0]


def signatur(basis):
    Gm = np.array([[np.trace(X @ Y).real for Y in basis] for X in basis])
    ew = np.linalg.eigvalsh((Gm + Gm.T)/2)
    sk = max(1.0, np.abs(ew).max())
    return int(np.sum(ew > 1e-9*sk)), int(np.sum(ew < -1e-9*sk)), int(np.sum(np.abs(ew) <= 1e-9*sk))


uvec = [np.array([1, 0, -1, 0.]), np.array([0, 1, 0, -1.])]
vvec = [np.array([1, 0, 1, 0.]), np.array([0, 1, 0, 1.])]
Inp = (np.outer(uvec[0], uvec[1]) - np.outer(uvec[1], uvec[0])).astype(complex)       # I_pi (z89)
NP = [-0.5*sum(s[i, j]*np.outer(uvec[i], vvec[j]) for i in range(2) for j in range(2)) for s in SIGMA]
NK = [-0.5*sum(s[i, j]*np.outer(vvec[i], uvec[j]) for i in range(2) for j in range(2)) for s in SIGMA]
NE = np.diag([-0.5, -0.5, 0.5, 0.5]).astype(complex)
MP = [1j*N for N in NP]
MKL = [ETA[m, m]*1j*NK[m] for m in range(4)]
MEt = 1j*NE
MQ = (MP[0] - MKL[0])/2                                                  # Q_0^- = (P_0 - K_0)/2 (z96, z98)
MD = MEt @ MQ - MQ @ MEt                                                 # [E, Q^-] = i D  <=>  [M_E, M_Q] = M_D
LOR = teilalgebra([lambda M: M @ S4 - S4 @ M, lambda M: M @ Inp + Inp @ M.T])
Mmn = {(m, n): 0.25*((MP[m] @ MKL[n] - MKL[n] @ MP[m]) - (MP[n] @ MKL[m] - MKL[m] @ MP[n])) for m, n in kombis}
TMS_np = np.array([[1, 0, -1, 0], [0, 1, 0, -1], [1, 0, 1, 0], [0, 1, 0, 1]], dtype=complex)

# ------------------------------------------------------------------ sympy: konforme Vektorfelder (wie z98)
XS = sp.symbols("x0:4", real=True)
SIGS = [sp.eye(2), sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
XM = sum((XS[k]*SIGS[k] for k in range(4)), sp.zeros(2))
TMS = sp.Matrix([[1, 0, -1, 0], [0, 1, 0, -1], [1, 0, 1, 0], [0, 1, 0, 1]])
TMSi = TMS.inv()
DETX = sp.expand(XM.det())
ETAs = sp.diag(1, -1, -1, -1)


def zu_sympy(M):
    return sp.Matrix(4, 4, [sp.nsimplify(complex(x).real) + sp.I*sp.nsimplify(complex(x).imag) for x in np.ravel(M)])


def vektorfeld(Mz):
    """Physikalische Orientierung (P_0 = +d_0): minus der Bewegung der Geraden L_x unter Z -> Z + eps M Z (z98)."""
    Mw = TMS*zu_sympy(Mz)*TMSi
    A, B, C, D = Mw[:2, :2], Mw[:2, 2:], Mw[2:, :2], Mw[2:, 2:]
    dx = A*XM - XM*D - sp.I*B - sp.I*XM*C*XM
    return [sp.expand(-(dx*SIGS[k]).trace()/2) for k in range(4)]


def killing_rest(xi, Om):
    div = sum(sp.diff(xi[m], XS[m]) for m in range(4))
    return sp.expand(sum(xi[m]*sp.diff(Om, XS[m]) for m in range(4)) - div*Om/4)


def ricci_skalar(gd, X, punkt):
    """Ricci-Skalar einer Diagonalmetrik (wie z98), am Punkt ausgewertet."""
    n = len(X)
    gi = [1/x for x in gd]
    Gam = [[[0]*n for _ in range(n)] for _ in range(n)]
    for l in range(n):
        for m in range(n):
            for k in range(n):
                val = 0
                if l == m:
                    val += sp.diff(gd[l], X[k])
                if l == k:
                    val += sp.diff(gd[l], X[m])
                if m == k:
                    val -= sp.diff(gd[m], X[l])
                Gam[l][m][k] = gi[l]*val/2
    R = 0
    for s_ in range(n):
        Ric = 0
        for r in range(n):
            Ric += sp.diff(Gam[r][s_][s_], X[r]) - sp.diff(Gam[r][r][s_], X[s_])
            for l in range(n):
                Ric += Gam[r][r][l]*Gam[l][s_][s_] - Gam[r][s_][l]*Gam[l][r][s_]
        R += gi[s_]*Ric
    return sp.nsimplify(sp.simplify(R.subs(punkt)))


def euler_dichte(gd, X):
    """Euler-Dichte E_4 = Riem^2 - 4 Ric^2 + R^2 einer Diagonalmetrik (R^a_bcd = d_c Gam^a_db - ..., MTW)."""
    n = len(X)
    Gam = [[[0]*n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                val = 0
                if a == b:
                    val += sp.diff(gd[a], X[c])
                if a == c:
                    val += sp.diff(gd[a], X[b])
                if b == c:
                    val -= sp.diff(gd[b], X[a])
                Gam[a][b][c] = sp.simplify(val/(2*gd[a]))
    Rm = {}
    for a, b in itertools.product(range(n), repeat=2):
        for c in range(n):
            Rm[a, b, c, c] = 0
            for d in range(c + 1, n):
                v = sp.simplify(sp.diff(Gam[a][d][b], X[c]) - sp.diff(Gam[a][c][b], X[d])
                                + sum(Gam[a][c][e]*Gam[e][d][b] - Gam[a][d][e]*Gam[e][c][b] for e in range(n)))
                Rm[a, b, c, d], Rm[a, b, d, c] = v, -v
    riem2 = sum(Rm[a, b, c, d]**2*gd[a]/(gd[b]*gd[c]*gd[d]) for a, b, c, d in itertools.product(range(n), repeat=4))
    ric = [[sum(Rm[a, b, a, d] for a in range(n)) for d in range(n)] for b in range(n)]
    ric2 = sum(ric[b][d]**2/(gd[b]*gd[d]) for b in range(n) for d in range(n))
    rs = sum(ric[b][b]/gd[b] for b in range(n))
    return sp.simplify(riem2 - 4*ric2 + rs**2), sp.simplify(rs)


def eigenebenen(M):
    """Eigenebenen einer 4x4-Matrix: (Eigenwert, Gram-Eigenwerte der (2,2)-Form, Punkt x^mu oder 'unendlich')."""
    w, V = np.linalg.eig(M)
    lams = []
    for lam in w:
        if not any(abs(lam - l2) < 1e-8 for l2 in lams):
            lams.append(lam)
    aus = []
    for lam in lams:
        B = V[:, np.where(np.abs(w - lam) < 1e-8)[0]]
        B = np.linalg.qr(B)[0]
        gram = B.conj().T @ G4 @ B
        wc = TMS_np @ B
        if abs(np.linalg.det(wc[2:])) < 1e-9:
            pkt = "unendlich"
        else:
            xx = -1j*wc[:2] @ np.linalg.inv(wc[2:])
            pkt = [complex(np.round(np.trace(xx @ s)/2, 10)) for s in SIGMA]
        aus.append((complex(np.round(lam, 10)), np.round(np.linalg.eigvalsh((gram + gram.conj().T)/2), 10).tolist(), pkt))
    return aus


def punkt_von_ebene(B):
    wc = TMS_np @ B
    if abs(np.linalg.det(wc[2:])) < 1e-9:
        return None
    return -1j*wc[:2] @ np.linalg.inv(wc[2:])


def ebene_von_punkt(x):
    return np.linalg.solve(TMS_np, np.vstack([1j*x, I2]))


def vektor6(M):
    """so(4,2)-Vektordarstellung auf Bivektoren: I -> M I + I M^T (z98, Der(I))."""
    paare = list(itertools.combinations(range(4), 2))
    Bs = []
    for i, j in paare:
        B = np.zeros((4, 4), complex)
        B[i, j], B[j, i] = 1, -1
        Bs.append(B)
    V = np.zeros((6, 6), complex)
    for k, B in enumerate(Bs):
        Y = M @ B + B @ M.T
        V[:, k] = [Y[i, j] for i, j in paare]
    return V


# =================================================================== Teil 1
say("Teil 1: Kappe und Modulzeit in der nativen radialen su(2,2)")
# --- 1.1 Kappengenerator, radiale sl(2,R) und Original-Partner-Quetschung
Kp = (dag(d1) @ dag(e1) + dag(d2) @ dag(e2)).tocsr()                    # K_+ = Summe d^+ e^+ (Paar Original/Partner)
Km = dag(Kp)
Qm = ((Pop[0] - Kop[0])/2).tocsr()                                       # Q_0^- = (P_0 - K_0)/2
Dil = (0.5j*(Km - Kp)).tocsr()                                           # Dilatation
Kcap = (2*np.pi*Qm).tocsr()                                              # Hislop-Longo-Normierung der Kappe
ok_q = klein(Qm - 0.5*(Kp + Km)) and klein(Eop - (Pop[0] + Kop[0])/2) and klein(Pop[0] - Eop - Qm)
ok_sl = klein(komm(Eop, Qm) - 1j*Dil) and klein(komm(Dil, Eop) - 1j*Qm) and klein(komm(Qm, Dil) + 1j*Eop)
ok_dil = all(klein(komm(Dil, Pop[m]) - 1j*Pop[m]) and klein(komm(Dil, KopL[m]) + 1j*KopL[m]) for m in range(4)) \
    and all(klein(komm(Dil, MopL[k])) for k in kombis)
ok_tw = all(klein(komm(X, ZT[a]) - sum((-1j*M)[a, b]*ZT[b] for b in range(4)))
            for X, M in ((Qm, MQ), (Dil, MD), (Eop, MEt)) for a in range(4))
ok_tw &= all(klein(komm(MopL[k], ZT[a]) - sum((-1j*Mmn[k])[a, b]*ZT[b] for b in range(4))) for k in kombis for a in range(4))
ew = {n: np.sort_complex(np.round(np.linalg.eigvals(M), 10)).tolist() for n, M in (("E", MEt), ("Q^-", MQ), ("D", MD))}
sig_rad = signatur([MEt, MQ, MD])
spur = {n: float(np.trace(M @ M).real) for n, M in (("E", MEt), ("Q^-", MQ), ("D", MD), ("P_0", MP[0]))}
ok_ew = (np.allclose(sorted(np.linalg.eigvals(MEt).imag), [-0.5, -0.5, 0.5, 0.5])
         and np.allclose(np.linalg.eigvals(MEt).real, 0)
         and np.allclose(sorted(np.linalg.eigvals(MQ).real), [-0.5, -0.5, 0.5, 0.5])
         and np.allclose(sorted(np.linalg.eigvals(MD).real), [-0.5, -0.5, 0.5, 0.5])
         and np.abs(MP[0] @ MP[0]).max() < 1e-14 and np.abs(MP[0]).max() > 0.1)
ew6 = {n: np.round(np.sort(np.linalg.eigvals(vektor6(M)).real if n != "E" else np.linalg.eigvals(vektor6(M)).imag), 10).tolist()
       for n, M in (("E", MEt), ("Q^-", MQ))}
ok_ew &= ew6["Q^-"] == [-1, 0, 0, 0, 0, 1] and ew6["E"] == [-1, 0, 0, 0, 0, 1] and \
    np.allclose(np.linalg.eigvals(vektor6(MQ)).imag, 0) and np.allclose(np.linalg.eigvals(vektor6(MEt)).real, 0)
ZAHLEN["radiale_sl2"] = {"Twistor-Eigenwerte": {k: [str(x) for x in v] for k, v in ew.items()},
                         "Vektordarstellung (R^(4,2)): Q^- reell, E imaginär": ew6,
                         "Spurform tr(M^2)": spur, "Signatur span{E, Q^-, D}": str(sig_rad),
                         "Q^- als Quetschung": "(1/2) Summe (d_i^+ e_i^+ + e_i d_i)", "K_cap": "2 pi Q^-"}
check(f"[E] Radiale sl(2,R) und Kappengenerator (Fock 8^4, Twistorraum): Q_0^- = (P_0 - K_0)/2 = (1/2) Summe_i (d_i^+ "
      f"e_i^+ + e_i d_i) ist die Zweimoden-Quetschung von Original (d) und Partner (e), E = (P_0 + K_0)/2, P_0 = E + Q^-; mit "
      f"D = (i/2)(K_- - K_+) gilt [E, Q^-] = i D, [D, E] = i Q^-, [Q^-, D] = -i E (su(1,1)), D ist die Dilatation ([D, P_mu] = "
      f"i P_mu, [D, K_mu] = -i K_mu, vertauscht mit M_mn); Twistorbilder exakt ([X, Z] = -i M_X Z). Klassen: E elliptisch "
      f"(Eigenwerte +-i/2, tr = {spur['E']:+.0f}), Q^- und D hyperbolisch (+-1/2, tr = {spur['Q^-']:+.0f}), P_0 nilpotent "
      f"(parabolisch, tr = {spur['P_0']:.0f}); Signatur {sig_rad}. In der Vektordarstellung auf R^(4,2) hat Q^- die Eigenwerte "
      f"{ew6['Q^-']} (reell: Boost mit Eigenwert +-1 auf zwei Nullstrahlen, den Spitzen), E die Eigenwerte i x {ew6['E']}: "
      f"Die Hislop-Longo-Normierung ist damit K_cap = 2 pi Q^-",
      ok_q and ok_sl and ok_dil and ok_tw and ok_ew and sig_rad == (2, 1, 0))

# --- 1.2 Eigenebenen (Punkte) und Killinghorizont
EB = {n: eigenebenen(M) for n, M in (("E", MEt), ("Q^-", MQ), ("D", MD))}
kernP = sla.null_space(MP[0])
ok_eb = (len(EB["E"]) == 2 and all(min(abs(np.array(g))) > 0.5 and np.sign(g[0]) == np.sign(g[1]) for _, g, _ in EB["E"])
         and sorted(round(p[0].imag) for _, _, p in EB["E"]) == [-1, 1] and all(abs(p[0].real) < 1e-9 and
                                                                             max(abs(np.array(p[1:]))) < 1e-9 for _, _, p in EB["E"]))
ok_eb &= (len(EB["Q^-"]) == 2 and all(max(abs(np.array(g))) < 1e-9 for _, g, _ in EB["Q^-"])
          and sorted(round(p[0].real) for _, _, p in EB["Q^-"]) == [-1, 1]
          and all(abs(p[0].imag) < 1e-9 and max(abs(np.array(p[1:]))) < 1e-9 for _, _, p in EB["Q^-"]))
ok_eb &= (len(EB["D"]) == 2 and all(max(abs(np.array(g))) < 1e-9 for _, g, _ in EB["D"])
          and sorted(str(p) if p == "unendlich" else "null" for _, _, p in EB["D"]) == ["null", "unendlich"]
          and all(p == "unendlich" or max(abs(np.array(p))) < 1e-9 for _, _, p in EB["D"]))
ok_eb &= kernP.shape[1] == 2 and punkt_von_ebene(kernP) is None
XI = {n: vektorfeld(M) for n, M in (("E", MEt), ("Q^-", MQ), ("D", MD), ("P_0", MP[0]), ("K_0", MKL[0]))}
q2 = sum(x**2 for x in XS)
ok_vf = (all(sp.expand(a - b) == 0 for a, b in zip(XI["P_0"], [1, 0, 0, 0]))
         and all(sp.expand(a - b) == 0 for a, b in zip(XI["E"], [(1 + q2)/2] + [XS[0]*XS[k] for k in (1, 2, 3)]))
         and all(sp.expand(a - b) == 0 for a, b in zip(XI["Q^-"], [(1 - q2)/2] + [-XS[0]*XS[k] for k in (1, 2, 3)]))
         and all(sp.expand(a + b) == 0 for a, b in zip(XI["D"], XS))
         and all(sp.expand(a - b) == 0 for a, b in zip(XI["K_0"], [q2] + [2*XS[0]*XS[k] for k in (1, 2, 3)])))
t_, r_ = sp.symbols("t r", real=True)
ACHSE = {XS[0]: t_, XS[1]: r_, XS[2]: 0, XS[3]: 0}
norm = {n: sp.factor(sp.expand(sum(ETAs[m, m]*XI[n][m]**2 for m in range(4)).subs(ACHSE))) for n in ("E", "Q^-")}
ok_hor = (sp.expand(norm["Q^-"] - sp.Rational(1, 4)*(1 - (t_ + r_)**2)*(1 - (t_ - r_)**2)) == 0
          and sp.expand(norm["E"] - sp.Rational(1, 4)*(1 + (t_ + r_)**2)*(1 + (t_ - r_)**2)) == 0)
ZAHLEN["eigenebenen"] = {n: [{"Eigenwert": str(l_), "Gram (2,2)-Form": g, "Punkt x^mu": str(p)} for l_, g, p in v]
                         for n, v in EB.items()}
ZAHLEN["eigenebenen"]["P_0"] = "Kern 2-dim, eine Ebene = Punkt unendlich (parabolisch, ein doppelter Fixpunkt)"
ZAHLEN["vektorfelder"] = {n: [str(c) for c in v] for n, v in XI.items()}
ZAHLEN["horizont"] = {"|xi_Q|^2 (t, r)": str(norm["Q^-"]), "|xi_E|^2 (t, r)": str(norm["E"])}
check(f"[E] Eigenebenen und Horizont (Inzidenz omega'' = i x c, sympy exakt): Die Eigenebenen von E sind (2,2)-definit "
      f"(Gram {EB['E'][0][1]}, {EB['E'][1][1]}) und gehören zu den komplexen Punkten x^0 = +-i (Vakuumpunkte): keine reelle "
      f"Kappe. Die Eigenebenen von Q^- sind null (totalisotrop) und gehören zu den reellen Punkten x = (+-1, 0, 0, 0): die "
      f"Spitzen des Doppelkegels O_1 = {{|t| + r < 1}}; die von D sind null mit den Punkten 0 und unendlich (Lichtkegel V_+); "
      f"P_0 ist nilpotent mit dem einzigen (doppelten) Fixpunkt unendlich. Vektorfelder (P_0 = d_t): xi_Q = ((1 - t^2 - |x|^2)/2, "
      f"-t x), xi_E = ((1 + t^2 + |x|^2)/2, t x), xi_D = -x, xi_K0 = (t^2 + |x|^2, 2 t x); |xi_Q|^2 = {norm['Q^-']} verschwindet "
      f"genau auf dem Rand von O_1 (Killinghorizont, zukunftsgerichtet innen), |xi_E|^2 = {norm['E']} > 0 überall (kein "
      f"Horizont, keine Kappe)", ok_eb and ok_vf and ok_hor)

# --- 1.3 Native Kappe aus E, tau und den Marken
so3 = teilalgebra([lambda M: M @ MEt - MEt @ M], basis=LOR)
zentr = teilalgebra([lambda M, R=R: M @ R - R @ M for R in so3])
zentr_su = teilalgebra([lambda M: np.array([[np.trace(M)]])], basis=zentr)
ungerade = teilalgebra([lambda M: S4 @ M @ S4 + M], basis=zentr_su)
gerade = teilalgebra([lambda M: S4 @ M @ S4 - M], basis=zentr_su)
orth = teilalgebra([lambda M: np.array([[np.trace(M @ MEt)]])], basis=ungerade)
ok_kap = (len(so3) == 3 and len(zentr) == 4 and len(zentr_su) == 3 and im_spann(1j*np.eye(4), zentr)
          and all(im_spann(X, zentr_su) for X in (MEt, MQ, MD)) and len(ungerade) == 2 and len(gerade) == 1
          and all(im_spann(X, ungerade) for X in (MEt, MQ)) and im_spann(MD, gerade) and len(orth) == 1)
Xn = orth[0]*(0.5/np.abs(np.linalg.eigvals(orth[0])).max())
ok_kap &= min(np.abs(Xn - MQ).max(), np.abs(Xn + MQ).max()) < 1e-10
orth_voll = teilalgebra([lambda M: np.array([[np.trace(M @ MEt)]])], basis=zentr_su)
ok_kap &= len(orth_voll) == 2 and all(im_spann(X, orth_voll) for X in (MQ, MD))
ok_tau = all(np.allclose(S4 @ X @ S4, s*X) for X, s in ((MEt, -1), (MQ, -1), (MD, 1))) and \
    all(np.allclose(S4 @ X @ S4, X) for X in LOR)
xr = sum(c*s for c, s in zip(rng.normal(size=4), SIGMA))
ok_tau &= np.allclose(punkt_von_ebene(S4 @ ebene_von_punkt(xr)), -xr)                 # tau: x -> -x, Fixpunkte 0, unendlich
bahn = {th: sla.expm(th*MEt) @ MQ @ sla.expm(-th*MEt) for th in (0.4, np.pi/2)}
ok_bahn = all(np.allclose(bahn[th], np.cos(th)*MQ + np.sin(th)*MD) for th in bahn) and np.allclose(bahn[np.pi/2], MD)
ok_bahn &= all(abs(np.trace(bahn[th] @ MEt)) < 1e-12 and im_spann(bahn[th], zentr_su) for th in bahn) and \
    not np.allclose(S4 @ bahn[0.4] @ S4, bahn[0.4]) and not np.allclose(S4 @ bahn[0.4] @ S4, -bahn[0.4])
ZAHLEN["native_kappe"] = {"dim so(3) = Lorentz cap C(E)": len(so3), "dim Zentralisator von so(3) in u(2,2)": len(zentr),
                          "dim in su(2,2)": len(zentr_su), "tau-ungerade": "span{E, Q^-}", "tau-gerade": "span{D}",
                          "Killing-orthogonal zu E und tau-ungerade": "+-Q^- (eindeutig)", "Killing-orthogonal zu E": "span{Q^-, D}",
                          "E-Bahn": "Ad(e^(theta E)) Q^- = cos theta Q^- + sin theta D (Halbsphären-Rauten am Zentrum)"}
check(f"[E] Native Kappe aus E und tau: so(3) = Lorentz cap C(E) (dim {len(so3)}); der Zentralisator von so(3) in u(2,2) "
      f"ist die radiale sl(2,R) plus Zentrum (dim {len(zentr)}, spurfrei {len(zentr_su)}); tau = Theta_1 wirkt darauf als "
      f"E -> -E, Q^- -> -Q^-, D -> D (und +1 auf Lorentz), auf Punkten als x -> -x (Fixpunkte 0 und unendlich). Der tau-ungerade "
      f"Teil ist span{{E, Q^-}}; Killing-orthogonal zu E bleibt genau eine Richtung, normiert +-Q^-: die Kappe O_1 (Halbsphäre "
      f"der S^3 zur Uhrzeit 0 um den tau-Fixpunkt, Spitzen t = +-Im x_Vak) ist durch (E, tau) bis auf die Orientierung (O_1 "
      f"oder O_1') bestimmt. Der tau-gerade Teil gibt +-D: der Lichtkegel V_+- mit den tau-Fixpunkten als Spitzen. Ohne tau "
      f"bilden die zulässigen Kappen die E-Bahn cos theta Q^- + sin theta D (alle Halbsphären-Rauten am Zentrum); E allein "
      f"legt keine Kappe fest (elliptisch, komplexe Fixpunkte), Ad(e^(pi E/2)) Q^- = D", ok_kap and ok_tau and ok_bahn)

# Marken: Uhrzeit auf der Zentrumsweltlinie (t = tan(theta/2)), Nullstellen der Felder
th_ = sp.Symbol("theta", real=True)
uhr = {}
for n_ in ("E", "Q^-", "D", "P_0", "K_0"):
    f0 = XI[n_][0].subs({XS[0]: t_, XS[1]: 0, XS[2]: 0, XS[3]: 0})
    uhr[n_] = sp.lambdify(th_, (2*f0/(1 + t_**2)).subs(t_, sp.tan(th_/2)), "mpmath")
ERW = {"E": lambda x: 1, "Q^-": lambda x: mp.cos(x), "D": lambda x: -mp.sin(x), "P_0": lambda x: 1 + mp.cos(x),
       "K_0": lambda x: 1 - mp.cos(x)}
pkt_th = [mp.mpf("0.3") + mp.mpf("0.7")*k for k in range(9)]
ok_uhr = all(abs(uhr[n_](x) - ERW[n_](x)) < mp.mpf(10)**-30 for n_ in uhr for x in pkt_th)
null_sets = {"E": set(), "Q^-": {sp.pi/2, 3*sp.pi/2}, "D": {sp.Integer(0), sp.pi}, "P_0": {sp.pi}, "K_0": {sp.Integer(0)}}
FORM = {"E": sp.Integer(1), "Q^-": sp.cos(th_), "D": -sp.sin(th_), "P_0": 1 + sp.cos(th_), "K_0": 1 - sp.cos(th_)}
ok_null = all(set(sp.solveset(FORM[n_], th_, sp.Interval.Ropen(0, 2*sp.pi))) == null_sets[n_] for n_ in FORM)
MARKEN = {sp.Integer(0), sp.pi/2, sp.pi, 3*sp.pi/2}                     # z^4 = 1, z = e^(i theta)
ok_mark = (null_sets["Q^-"] | null_sets["D"]) == MARKEN and not (null_sets["Q^-"] & null_sets["D"])
ok_mark &= {sp.Mod(-x, 2*sp.pi) for x in MARKEN} == MARKEN and {sp.Mod(-x, 2*sp.pi) for x in null_sets["D"]} == null_sets["D"]
ok_mark &= {sp.Mod(x + sp.pi/2, 2*sp.pi) for x in null_sets["D"]} == null_sets["Q^-"]
ZAHLEN["marken"] = {"Felder in Uhrzeit": {"E": "1", "Q^-": "cos theta", "D": "-sin theta", "P_0": "1 + cos theta",
                                          "K_0": "1 - cos theta"},
                    "Nullstellen": {k: sorted(str(x) for x in v) for k, v in null_sets.items()},
                    "Marken z^4 = 1": "theta in {0, pi/2, pi, 3 pi/2}"}
check("[E/C] Marken: Auf der Zentrumsweltlinie in Uhrzeit theta (t = tan(theta/2), E = Drehung) sind die Felder exakt "
      "die der Naht in z72-Konvention: E = 1, Q^- = cos theta (H_1), D = -sin theta (H_2 bis aufs Vorzeichen), P_0 = 1 + cos "
      "theta, K_0 = 1 - cos theta (30 Stellen). Nullstellen auf [0, 2 pi): E keine, Q^- {pi/2, 3 pi/2} (Spitzen t = +-1), D "
      "{0, pi} (Punkte 0, unendlich = tau-Fixpunkte), P_0 {pi}, K_0 {0}. Unter z83s Identifikation (Nahtdrehung = Uhr, Marken "
      "z^4 = 1 bei theta = k pi/2, tau wirkt als theta -> -theta) [C] sind die vier Marken genau die Spitzen der beiden "
      "tau-Eigenkappen: Antipodenpaar {1, -1} = V_+-, Paar {i, -i} = O_1, O_1'. Die Vierteldrehung vertauscht die Paare "
      "(D <-> Q^-): Die Marken liefern die Spitzenmenge, unterscheiden die Kappen aber nicht; das tut tau",
      ok_uhr and ok_null and ok_mark)
befund("E/C: Eine native Kappe existiert: O_1, die Halbsphäre der S^3 zur Uhrzeit 0 um den tau-Fixpunkt (Doppelkegel mit den "
       "Spitzen t = +-Im x_Vak, den reellen 'Wick-Partnern' der Vakuumpunkte x^0 = +-i), eindeutig aus (E, tau) bis auf die "
       "Orientierung; ihr Modulgenerator in Hislop-Longo-Normierung ist K_cap = 2 pi Q^-. E allein gibt keine Kappe, tau allein "
       "nur das Zentrum (und mit D den Lichtkegel V_+), die Marken liefern unter z83s Identifikation die vier Spitzen beider "
       "tau-Eigenkappen", True)
laufzeit()

# --- 1.4 Modulfluss: Vakuum (Hislop-Longo nach dem Netz) gegen das Original
u1, u2, v_ = sp.symbols("u1 u2 v", real=True)
a_, b_ = sp.symbols("a b", positive=True)
t_u = sp.tanh(u1/2)
ok_fluss = sp.simplify(sp.diff(t_u, u1) - XI["Q^-"][0].subs({XS[0]: t_u, XS[1]: 0, XS[2]: 0, XS[3]: 0})) == 0
ta, tb = (a_ - 1)/(a_ + 1), (b_ - 1)/(b_ + 1)                           # t = tanh(u/2) mit a = e^(u1), b = e^(u2)
dta, dtb = 2*a_/(a_ + 1)**2, 2*b_/(b_ + 1)**2                           # dt/du
W_M = 1/(4*sp.pi**2*(ta - tb)**2)                                       # Wightman, masseloser Skalar Delta = 1 (z83)
F_ab = sp.simplify(W_M*dta*dtb)
G_ab = a_*b_/(4*sp.pi**2*(a_ - b_)**2)                                  # = 1/(16 pi^2 sinh^2((u1 - u2)/2))
ok_stat = sp.simplify(F_ab - G_ab) == 0 and sp.simplify(
    (sp.sinh((u1 - u2)/2)**2).rewrite(sp.exp) - (sp.exp(u1) - sp.exp(u2))**2/(4*sp.exp(u1 + u2))) == 0
Gv = 1/(16*sp.pi**2*sp.sinh(v_/2)**2)
ok_kms = sp.simplify(Gv.subs(v_, v_ - 2*sp.pi*sp.I) - Gv) == 0
G_num = sp.lambdify(v_, Gv, "mpmath")
streifen = max(abs(G_num(mp.mpf(x) - 1j*mp.pi)) for x in np.linspace(-6, 6, 61))
pol = abs(G_num(mp.mpf("1e-8") - 2j*mp.pi))
ok_kms &= streifen < 1 and pol > 1e10
# dS-Rahmen: Omega_dS = 1 - det x (z98, sigma = -1), auf der Achse 1 - t^2; Eigenzeit tau_p = artanh t = u/2
tau_p = sp.integrate(1/(1 - t_**2), t_)
ok_ds = sp.simplify(sp.diff(tau_p, t_) - 1/(1 - t_**2)) == 0 and \
    sp.simplify(sp.diff(sp.atanh(t_u), u1) - sp.Rational(1, 2)) == 0 and sp.atanh(t_u).subs(u1, 0) == 0
W_dS = sp.simplify((1 - ta**2)*(1 - tb**2)*W_M)                         # konform gekoppelt, Weylgewicht 1
H_ = sp.Integer(2)                                                       # H = 1/l, l = 1/2 (z98)
ok_ds &= sp.simplify(W_dS - 4*G_ab) == 0 and sp.simplify(
    (H_**2/(16*sp.pi**2*sp.sinh(H_*(u1/2 - u2/2)/2)**2)).rewrite(sp.exp).subs({sp.exp(u1): a_, sp.exp(u2): b_})
    - 4*G_ab) == 0
l_dS = sp.Rational(1, 2)
T_GH = 1/(2*sp.pi*l_dS)
ZAHLEN["vakuum_kms"] = {"Fluss auf der Achse": "t = tanh(u/2)", "F(u1, u2)": "1/(16 pi^2 sinh^2((u1 - u2)/2))",
                        "KMS-Streifenbreite (u)": "2 pi", "max|G| auf Im v = -pi": float(streifen),
                        "|G| bei v = -2 pi i": float(pol), "Eigenzeit dS": "tau_p = artanh t = u/2",
                        "W_dS": "H^2/(16 pi^2 sinh^2(H dtau/2)), H = 2", "T_GH (Vakuumeinheiten)": str(T_GH)}
check(f"[E, C nach dem Netz] Vakuum auf der Kappenachse (Hislop-Longo in der erlaubten Reihenfolge: freies konformes Netz "
      f"des Leiterfelds aus z83/z90, dann der Satz): Der Kappenfluss ist auf der Achse t = tanh(u/2); die Vakuum-"
      f"Wightmanfunktion des Leiterfelds (Delta = 1) im Rahmen der Modulzeit u ist exakt 1/(16 pi^2 sinh^2((u1 - u2)/2)): "
      f"stationär und KMS mit G(v - 2 pi i) = G(v), analytisch im Streifen -2 pi < Im v < 0 (|G| <= {float(streifen):.3f} auf der "
      f"Mittellinie, Pol bei v = -2 pi i): beta = 2 pi in u, also KMS bei beta = 1 bezüglich K_cap = 2 pi Q^-. Im dS-Rahmen "
      f"g = eta/(1 - det x)^2 (z98, sigma = -1) ist die Eigenzeit des Zentrumsbeobachters tau_p = artanh t = u/2, die "
      f"Zweipunktfunktion die Bunch-Davies-Form H^2/(16 pi^2 sinh^2(H dtau/2)) mit H = 1/l = 2: T_GH = 1/(2 pi l) = {T_GH}",
      ok_fluss and ok_stat and ok_kms and ok_ds)

# Original: zentraler, periodischer Modulfluss; thermische Zweipunktfunktion nicht stationär unter dem Kappenfluss
ok_zentral = all(klein(komm(Dop, X)) for X in [Eop, Qm, Dil] + Pop + KopL + list(MopL.values()))
q_sc = float(np.linalg.norm(Qm[:, SC].toarray()))
d_sc = float(np.abs(Dop.diagonal()[SC]).max())
ok_zentral &= d_sc == 0 and q_sc > 1


def W_zyl(dt):
    return 1/(16*mp.pi**2*mp.sin(dt/2)**2)


def F_beta(x1, x2, beta, K=3):
    tc = lambda u: 2*mp.atan(mp.tanh(u/2))
    dtc = tc(x1) - tc(x2)
    W = W_zyl(dtc) if beta is None else sum(W_zyl(dtc + 1j*k*beta) for k in range(-K, K + 1))
    return mp.sech(x1)*mp.sech(x2)*W                                    # d tau_c/du = sech u


U1, U2, S_ = mp.mpf("0.4"), mp.mpf("-0.3"), mp.mpf("1.2")
G_ref = 1/(16*mp.pi**2*mp.sinh((U1 - U2)/2)**2)
abw_vak = abs(F_beta(U1, U2, None) - G_ref)/G_ref
stat_vak = abs(F_beta(U1 + S_, U2 + S_, None) - F_beta(U1, U2, None))/G_ref
F0 = F_beta(U1, U2, BETA_P)
stat_orig = abs(F_beta(U1 + S_, U2 + S_, BETA_P) - F0)/abs(F0)
stat_orig2 = abs(F_beta(U1 - S_, U2 - S_, BETA_P) - F0)/abs(F0)
ok_orig = abw_vak < mp.mpf(10)**-35 and stat_vak < mp.mpf(10)**-35 and stat_orig > mp.mpf(10)**-13 and \
    stat_orig2 > mp.mpf(10)**-13 and stat_orig < mp.mpf(10)**-9
ZAHLEN["original"] = {"[N_d - N_e, alle Erzeuger]": 0.0, "|Q^- auf Skalarsektor|": q_sc,
                      "max|N_d - N_e| auf Skalarsektor": d_sc,
                      "Vakuum: Abweichung von 1/(16 pi^2 sinh^2)": mp.nstr(abw_vak, 3),
                      "Vakuum: Nichtstationarität": mp.nstr(stat_vak, 3),
                      "Original (beta' = 8 pi): Nichtstationarität (s = +1.2, -1.2)": [mp.nstr(stat_orig, 5),
                                                                                         mp.nstr(stat_orig2, 5)],
                      "e^(-8 pi)": mp.nstr(mp.e**(-8*mp.pi), 5)}
check(f"[E] Original gegen Kappenfluss: Der Modulgenerator des Originals auf Fock 8^4 ist die Ladung N_d - N_e (z98: D, "
      f"zentral), sie vertauscht mit E, Q^-, der Dilatation, P, K, M und ist auf dem Skalarsektor null, wo Q^- nicht verschwindet "
      f"(|Q^-| = {q_sc:.3f}): Q^- ist keine Funktion des Originalgenerators, dessen Fluss periodisch ist (ganzzahliges "
      f"Spektrum; z70: spec K in pi Z). Als 4D-Zustand (Produktlesart z83: KMS bezüglich E bei beta' = 8 pi, Bildsumme in "
      f"Uhrzeit) ist die Achsen-Zweipunktfunktion im Rahmen u nicht stationär: relative Änderung unter u -> u +- 1,2 = "
      f"{mp.nstr(stat_orig, 4)}, {mp.nstr(stat_orig2, 4)} (Ordnung e^(-8 pi) = {mp.nstr(mp.e**(-8*mp.pi), 3)}), das Vakuum "
      f"dagegen {mp.nstr(stat_vak, 2)}: Das Original ist unter dem Kappenfluss nicht invariant, also für ihn bei keinem beta "
      f"KMS", ok_zentral and ok_orig)

# exakte Beziehungen: K_cap = 2 pi (P_0 - E) = 2 pi P_0 - K_Uhr/2, großer Kappenlimes -> P_0
th_g = 1.7
Rg = np.exp(th_g)
MQR = sla.expm(th_g*MD) @ MQ @ sla.expm(-th_g*MD)
ok_rel = klein(Kcap - 2*np.pi*(Pop[0] - Eop)) and klein(Kcap + (4*np.pi*Eop)/2 - 2*np.pi*Pop[0])
ok_rel &= np.allclose(MQR, (Rg*MP[0] - MKL[0]/Rg)/2) and \
    np.allclose(2/Rg*MQR - MP[0], -MKL[0]/Rg**2)
spitzen_R = sorted(round(p[0].real, 8) for _, _, p in eigenebenen(MQR) if p != "unendlich")
ok_rel &= np.allclose(spitzen_R, [-Rg, Rg])
ZAHLEN["beziehungen"] = {"K_cap": "2 pi (P_0 - E) = 2 pi P_0 - K_Uhr/2, K_Uhr = 4 pi E (z72: Uhr 4 pi L_0)",
                         "Dilatierte Kappe": "Q^-_R = (R P_0 - K_0/R)/2, Spitzen t = +-R, l = R/2, Lambda = 12/R^2",
                         "Spitzen bei R = e^1.7": spitzen_R, "Grenzwert": "P_0 = lim (2/R) Q^-_R = lim K_(O_R)/(pi R)"}
check(f"[E] Exakte Beziehungen: K_cap = 2 pi Q^- = 2 pi (P_0 - E); mit der Uhr des Originals K_Uhr = beta L_0 = 4 pi E "
      f"(z72) gilt K_cap + K_Uhr/2 = 2 pi P_0, der Faktor 1/2 ist 2 pi/beta. Die Dilatation verschiebt die Kappe: "
      f"Ad(e^(theta D)) Q^- = (R P_0 - K_0/R)/2 mit Spitzen t = +-R (R = e^theta, numerisch {spitzen_R}), und (2/R) Q^-_R = "
      f"P_0 - K_0/R^2 -> P_0: z72s parabolische Zeit ist der flache Grenzfall der Kappenzeit, P_0 = lim K_(O_R)/(pi R) "
      f"(l = R/2, Lambda = 12/R^2 -> 0)", ok_rel)
befund("E/C: Der Kappenfluss ist der Modulfluss des Vakuumsektors (Hislop-Longo, bedingt auf das freie Netz des "
       "Leiterfelds), nicht der des Originals: Dessen Modulgenerator ist zentral (Helizität/Ladung, periodisch) bzw. "
       "einseitig die Uhr 4 pi E, und sein Zustand ist unter dem Kappenfluss nicht invariant (Abweichung der Ordnung e^(-8 pi)). "
       "Alle drei Zeiten liegen in einer radialen sl(2,R): Uhr E elliptisch, Minkowskizeit P_0 = E + Q^- parabolisch "
       "(flacher Grenzfall der Kappen), Kappenzeit Q^- hyperbolisch", True)
laufzeit()

# =================================================================== Teil 2
say("\nTeil 2: dS-Kinematik auf der radialen Leiter")
# --- 2.1 Casimir von so(4,1) und so(3,2), K-Typen, Serie


def casimir(Mtw, Xs):
    """Kanonischer (Killing-)Casimir: [M_i, M_j] = f_ij^k M_k, B = tr(ad ad), C = -Summe B^ij X_i X_j (rho(M) = -i X)."""
    n = len(Mtw)
    A = np.array([np.concatenate([M.real.ravel(), M.imag.ravel()]) for M in Mtw]).T
    f = np.zeros((n, n, n))
    res = 0.0
    for i in range(n):
        for j in range(n):
            Cm = Mtw[i] @ Mtw[j] - Mtw[j] @ Mtw[i]
            b = np.concatenate([Cm.real.ravel(), Cm.imag.ravel()])
            x = np.linalg.lstsq(A, b, rcond=None)[0]
            res = max(res, float(np.linalg.norm(A @ x - b)))
            f[i, j] = x
    B = np.einsum("ilk,jkl->ij", f, f)
    Bi = np.linalg.inv(B)
    C = (-sum(Bi[i, j]*(Xs[i] @ Xs[j]) for i in range(n) for j in range(n) if abs(Bi[i, j]) > 1e-14)).tocsr()
    C_tw = sum(Bi[i, j]*(Mtw[i] @ Mtw[j]) for i in range(n) for j in range(n))      # Twistor = Spinordarstellung
    return C, B, res, C_tw


CAS = {}
for sg, nm in ((1, "so(3,2)"), (-1, "so(4,1)")):
    Mtw = [(MP[m] + sg*MKL[m])/2 for m in range(4)] + [Mmn[k] for k in kombis]
    Xs = [((Pop[m] + sg*KopL[m])/2).tocsr() for m in range(4)] + [MopL[k] for k in kombis]
    C, B, res, C_tw = casimir(Mtw, Xs)
    Csc = C[:, SC].toarray()
    c0 = complex(Csc[0, 0])
    konst = float(np.abs(Csc[SC, :] - c0*np.eye(len(SC))).max())
    raus = float(np.abs(np.delete(Csc, SC, axis=0)).max())
    kommutiert = all(klein(komm(C, X), 1e-9, LOW4) for X in Xs)
    neg = int(np.sum(np.linalg.eigvalsh(B) < -1e-9))
    CAS[nm] = {"Mtw": Mtw, "Xs": Xs, "C": C, "c0": c0, "konstant": konst, "raus": raus, "kommutiert": kommutiert,
               "kompakte Richtungen": neg, "res": res, "spinor": complex(C_tw[0, 0]),
               "spinor konstant": float(np.abs(C_tw - C_tw[0, 0]*np.eye(4)).max())}
c_ads, c_ds = CAS["so(3,2)"]["c0"].real, CAS["so(4,1)"]["c0"].real
Cphys = 6*c_ads                                                          # so(5): C_phys = 6 C_Killing (Spinor 5/2, Vektor 4)
Delta = sp.symbols("Delta")
D_loes = sorted(sp.solve(sp.Eq(Delta*(Delta - 3), sp.nsimplify(round(6*c_ds, 12))), Delta))
# so(4) = maximal kompakt in so(4,1) = K-Typ-so(4) der su(2,2) (Kommutant von E)
so4 = teilalgebra([lambda M: M @ MEt - MEt @ M], basis=CAS["so(4,1)"]["Mtw"])
X4 = [sum(c*X for c, X in zip(koeff(M, CAS["so(4,1)"]["Mtw"]), CAS["so(4,1)"]["Xs"])).tocsr() for M in so4]
C4, _, _, _ = casimir(so4, X4)
c4 = {}
for Ev in (1, 2, 3):
    idx = [i for i in SC if abs(Eop.diagonal()[i].real - Ev) < 1e-9]
    blk = C4[idx, :][:, idx].toarray()
    c4[Ev] = (float(blk[0, 0].real), float(np.abs(blk - blk[0, 0]*np.eye(len(idx))).max()), len(idx))
ok_k = len(so4) == 6 and all(v[1] < 1e-10 for v in c4.values()) and abs(c4[1][0]) < 1e-10 and \
    abs(c4[3][0]/c4[2][0] - 8/3) < 1e-10 and [v[2] for v in c4.values()] == [1, 4, 9]
ok_cas = (all(v["konstant"] < 1e-10 and v["raus"] < 1e-10 and v["kommutiert"] and v["res"] < 1e-12 for v in CAS.values())
          and abs(c_ads - c_ds) < 1e-12 and abs(Cphys + 2) < 1e-10 and D_loes == [1, 2]
          and all(abs(v["spinor"] - 5/12) < 1e-12 and v["spinor konstant"] < 1e-12 for v in CAS.values())
          and CAS["so(3,2)"]["kompakte Richtungen"] == 4 and CAS["so(4,1)"]["kompakte Richtungen"] == 6)
ZAHLEN["casimir"] = {"C_Killing so(3,2) auf dem Skalar": c_ads, "C_Killing so(4,1) auf dem Skalar": c_ds,
                     "C_Killing auf Twistoren (Spinor, beide Formen)": [CAS[k]["spinor"].real for k in CAS],
                     "Eichung": "Spinor C_phys = n(n-1)/8 = 5/2 = 6 x 5/12",
                     "Abweichung von konstant (AdS, dS)": [CAS["so(3,2)"]["konstant"], CAS["so(4,1)"]["konstant"]],
                     "C_phys = 6 C_Killing": Cphys, "AdS: E_0(E_0 - 3) = -2": "E_0 = 1, 2 (D(1,0) + D(2,0))",
                     "dS: Delta(Delta - 3) = -2": [str(x) for x in D_loes], "m^2 l^2 = Delta(3 - Delta)": 2,
                     "so(4)-Casimir je Stufe E = 1, 2, 3 (Wert, Abweichung, dim)": c4,
                     "Hauptserie bräuchte Delta(3 - Delta) >= 9/4": True}
check(f"[E] dS-Darstellung auf der radialen Leiter (Skalarsektor N_d = N_e, E = 1, 2, 3; Fock 8^4): span{{Q^-_a, M_ab}} ist "
      f"so(4,1) (Killingform: {CAS['so(4,1)']['kompakte Richtungen']} kompakte Richtungen; so(3,2): "
      f"{CAS['so(3,2)']['kompakte Richtungen']}), die Darstellung ist unitär (hermitesche Erzeuger). Die kanonischen "
      f"Casimire beider Unteralgebren sind auf dem ganzen Skalarsektor konstant und gleich: C_K = {c_ads:.12f} (AdS), "
      f"{c_ds:.12f} (dS), vertauschen mit allen zehn Erzeugern; geeicht an der Twistor- (Spinor-)Darstellung, wo C_K = "
      f"{CAS['so(4,1)']['spinor'].real:.12f} = 5/12 und C_phys = 5/2 ist, gilt C_phys = 6 C_K = {Cphys:.6f} = E_0(E_0 - 3) "
      f"für E_0 = 1 und 2 "
      f"(AdS: D(1,0) + D(2,0)) und = Delta(Delta - 3) mit Delta = {D_loes} für dS. Die maximal kompakte so(4) der so(4,1) "
      f"ist die K-Typ-so(4) der su(2,2) (Kommutant von E, dim {len(so4)}): ihr Casimir ist auf den Stufen E = 1, 2, 3 "
      f"konstant ({c4[1][0]:.3f}, {c4[2][0]:.3f}, {c4[3][0]:.3f}; Verhältnis 8/3 wie (j, j) mit j = 0, 1/2, 1), jeder "
      f"sphärische K-Typ (j, j) kommt genau einmal vor. Also: die irreduzible sphärische Komplementärserie mit Delta = 1 "
      f"(äquivalent 2), m^2 l^2 = Delta(3 - Delta) = 2 = 12 xi mit xi = 1/6 (konform gekoppelter Skalar, Bunch-Davies); "
      f"die Hauptserie (Delta(3 - Delta) >= 9/4) ist ausgeschlossen", ok_cas and ok_k)

# --- 2.2 Spektrum von Q^- je Drehimpuls: D_(l+1), sech-Gesetz, symmetrisch, stetig
rot = [MopL[(1, 2)], MopL[(1, 3)], MopL[(2, 3)]]
rot_tw = [Mmn[(1, 2)], Mmn[(1, 3)], Mmn[(2, 3)]]
kap = koeff(rot_tw[0] @ rot_tw[1] - rot_tw[1] @ rot_tw[0], rot_tw)
kappa = float(np.abs(kap).max())
J2 = (sum((X @ X for X in rot), sps.csr_matrix((DIMV, DIMV)))/kappa**2).tocsr()
ok_rot = all(klein(komm(X, Eop)) for X in rot)
leiter = {}
for Ev in (1, 2, 3):
    idx = [i for i in SC if abs(Eop.diagonal()[i].real - Ev) < 1e-9]
    Kmat = Km[:, idx].toarray()
    kern = sla.null_space(Kmat, rcond=1e-10)
    vecs = np.zeros((DIMV, kern.shape[1]), complex)
    vecs[idx, :] = kern
    j2 = np.linalg.eigvalsh(vecs.conj().T @ (J2 @ vecs))
    leiter[Ev] = (kern.shape[1], np.round(j2, 10).tolist())
ok_leit = leiter[1][0] == 1 and leiter[2][0] == 3 and leiter[3][0] == 5 and \
    all(np.allclose(leiter[Ev][1], (Ev - 1)*Ev) for Ev in (1, 2, 3))
# Kettenelemente <n+1|Q^-|n> = sqrt((n+1)(n+2h))/2 auf Fock (h = 1: n = 0, 1; h = 2: n = 0)
kette_ok = True
kette_tab = {}
for h, start in ((1, OMV), (2, None)):
    if start is None:
        idx = [i for i in SC if abs(Eop.diagonal()[i].real - 2) < 1e-9]
        kern = sla.null_space(Km[:, idx].toarray(), rcond=1e-10)
        start = np.zeros(DIMV, complex)
        start[idx] = kern[:, 0]
    vs = [start/np.linalg.norm(start)]
    for n in range(2 if h == 1 else 1):
        w = Kp @ vs[-1]
        vs.append(w/np.linalg.norm(w))
    for n in range(len(vs) - 1):
        amp = np.vdot(vs[n + 1], Qm @ vs[n])
        soll = np.sqrt((n + 1)*(n + 2*h))/2
        rest = np.linalg.norm(Qm @ vs[n] - amp*vs[n + 1] - (np.vdot(vs[n - 1], Qm @ vs[n])*vs[n - 1] if n > 0 else 0))
        kette_tab[f"h={h}, n={n}"] = (float(abs(amp)), float(soll), float(rest))
        kette_ok &= abs(abs(amp) - soll) < 1e-12 and rest < 1e-12 and abs(np.vdot(vs[n], Qm @ vs[n])) < 1e-12


def dh_matrizen(h, N):
    n = np.arange(N)
    off = np.sqrt((n[:-1] + 1)*(n[:-1] + 2*h))
    L0 = sps.diags(h + n.astype(float)).tocsr()
    H1 = sps.diags([off/2, off/2], [1, -1]).tocsr()
    H2 = sps.diags([1j*off/2, -1j*off/2], [1, -1]).tocsr()
    return L0, H1, H2


sech_tab = {}
sech_ok = True
for h in (1, 2, 3):
    L0h, H1h, _ = dh_matrizen(h, 4000)
    e0 = np.zeros(4000, complex)
    e0[0] = 1
    for s in (0.5, 1.5, 3.0):
        amp = complex(spla.expm_multiply(-1j*s*H1h, e0)[0])
        soll = np.cosh(s/2)**(-2*h)
        sech_tab[f"h={h}, s={s}"] = (amp.real, amp.imag, soll)
        sech_ok &= abs(amp - soll) < 1e-10


def rho_h(lam, h):
    return 2**(2*h)*abs(mp.gamma(h + 1j*lam))**2/(2*mp.pi*mp.gamma(2*h))


ft_ok = True
ft_tab = {}
for h in (1, 2):
    norm_h = mp.quad(lambda x: rho_h(x, h), [-mp.inf, 0, mp.inf])
    for s in (0.7, 2.0):
        ft = mp.quad(lambda x: rho_h(x, h)*mp.cos(s*x), [-mp.inf, 0, mp.inf])
        ft_tab[f"h={h}, s={s}"] = (mp.nstr(ft, 15), mp.nstr(mp.cosh(s/2)**(-2*h), 15))
        ft_ok &= abs(ft - mp.cosh(s/2)**(-2*h)) < mp.mpf(10)**-20
    ft_ok &= abs(norm_h - 1) < mp.mpf(10)**-20 and abs(rho_h(mp.mpf("0.8"), h) - rho_h(mp.mpf("-0.8"), h)) < mp.mpf(10)**-30
Par = sps.diags(((-1.0)**(OCC[:, 0] + OCC[:, 1]))).tocsr()
ok_par = klein(Par @ Qm @ Par + Qm) and klein(Par @ Dil @ Par + Dil) and klein(Par @ Eop @ Par - Eop)
abst = {}
for N in (500, 2000, 8000):
    _, H1h, _ = dh_matrizen(1, N)
    ev = np.linalg.eigvalsh(H1h.toarray()) if N <= 2000 else sla.eigh_tridiagonal(
        np.zeros(N), np.sqrt((np.arange(N - 1) + 1)*(np.arange(N - 1) + 2))/2, eigvals_only=True)
    pos = np.sort(ev[ev > 0])
    abst[N] = (float(pos[1] - pos[0]), float(ev.min()), float(ev.max()))
ok_stetig = abst[500][0] > abst[2000][0] > abst[8000][0] and all(abs(v[1] + v[2]) < 1e-6*abs(v[2]) for v in abst.values())
ZAHLEN["spektrum_Q"] = {"so(3)-Struktur kappa": kappa, "Kern K_- je E (dim, J^2)": leiter, "Kette (|amp|, soll, Rest)": kette_tab,
                        "sech-Gesetz <0|e^(-i s H_1)|0> (Re, Im, cosh^(-2h)(s/2))": sech_tab,
                        "Spektraldichte rho_h = 2^(2h)|Gamma(h + i lam)|^2/(2 pi Gamma(2h)), Fourier": ft_tab,
                        "D_1-Abschneidung: (Abstand bei 0, min, max)": abst}
check(f"[E] Spektrum der Kappenenergie Q^- = K_cap/(2 pi): Auf dem Skalarsektor sind die Tiefstgewichte der radialen "
      f"sl(2,R) (Kern von K_-) bei E = 1, 2, 3 genau die Drehimpulse l = 0, 1, 2 (dim {leiter[1][0]}, {leiter[2][0]}, "
      f"{leiter[3][0]}, J^2 = l(l+1)): die Leiter ist Summe_l D_(l+1) x (2l+1), die Kettenelemente <n+1|Q^-|n> = sqrt((n+1)"
      f"(n+2h))/2 exakt. Auf D_h gilt <0|e^(-i s Q^-)|0> = cosh(s/2)^(-2h) (h = 1, 2, 3; N = 4000, Fehler < 1e-10), "
      f"Fourier-Partner der normierten, geraden Dichte rho_h(lambda) = 2^(2h)|Gamma(h + i lambda)|^2/(2 pi Gamma(2h)) "
      f"(mpmath, 1e-20): absolut stetig auf ganz R, symmetrisch; (-1)^(N_d) kehrt Q^- und D um, E nicht. Abschneidungen von "
      f"D_1: Abstand bei 0 = {abst[500][0]:.3f}, {abst[2000][0]:.3f}, {abst[8000][0]:.3f} (N = 500, 2000, 8000), Spektrum "
      f"symmetrisch [{abst[8000][1]:.0f}, {abst[8000][2]:.0f}]. Vergleich: E diskret {{h, h+1, ...}}, P_0 = [0, inf) "
      f"(Gamma-Gesetz, z72/z83), Q^- = R ohne Grundzustand",
      ok_rot and ok_leit and kette_ok and sech_ok and ft_ok and ok_par and ok_stetig)
laufzeit()

# --- 2.3 Lambda als Funktion der Zeitwahl Q_0 = a E + b Q^-, und Positivität
s_sym = sp.Symbol("s")
GEN_P = [vektorfeld(M) for M in MP]
GEN_K = [vektorfeld(M) for M in MKL]
pkt = {XS[0]: sp.Rational(1, 5), XS[1]: sp.Rational(1, 3), XS[2]: sp.Rational(-1, 4), XS[3]: sp.Rational(1, 7)}
lam_tab = {}
ok_lam = True
for a_w, b_w in ((1, 0), (sp.Rational(3, 4), sp.Rational(1, 4)), (sp.Rational(1, 2), sp.Rational(1, 2)),
                 (sp.Rational(1, 4), sp.Rational(3, 4)), (0, 1)):
    lam_w, mu_w = a_w + b_w, a_w - b_w                                   # Q_mu = (lam P_mu + mu K_mu)/2
    eqs = []
    for m in range(4):
        xi = [sp.expand(lam_w*p + mu_w*k) for p, k in zip(GEN_P[m], GEN_K[m])]
        eqs += sp.Poly(killing_rest(xi, 1 + s_sym*DETX), *XS).coeffs()
    los = sp.solve(eqs, s_sym, dict=True)
    s_w = los[0][s_sym] if los else None
    Om = 1 + s_w*DETX
    R_w = ricci_skalar([-1/Om**2] + [1/Om**2]*3, XS, pkt)
    Lam_w = R_w/4
    MQ0 = complex(a_w)*MEt + complex(b_w)*MQ
    B_w = float(np.trace(MQ0 @ MQ0).real)
    lam_tab[f"a={a_w}, b={b_w}"] = {"s (Omega = 1 + s det x)": str(s_w), "R": str(R_w), "Lambda": str(Lam_w),
                                    "tr Q_0^2": round(B_w, 12)}
    ok_lam &= (len(los) == 1 and sp.simplify(s_w - mu_w/lam_w) == 0 and sp.simplify(Lam_w + 12*mu_w/lam_w) == 0
               and int(sp.sign(Lam_w)) == int(np.sign(round(B_w, 12))))
pos_tab = {}
for a_w, b_w in ((1.0, 0.5), (1.0, 1.0), (1.0, 1.5)):
    mins = []
    for N in (200, 800, 3200):
        n = np.arange(N)
        mins.append(float(sla.eigh_tridiagonal(a_w*(1 + n.astype(float)), b_w*np.sqrt((n[:-1] + 1)*(n[:-1] + 2))/2,
                                                eigvals_only=True, select="i", select_range=(0, 0))[0]))
    pos_tab[f"a={a_w}, b={b_w}"] = [round(x, 6) for x in mins]
ok_pos = (all(abs(m - np.sqrt(0.75)) < 1e-6 for m in pos_tab["a=1.0, b=0.5"])
          and all(0 < m for m in pos_tab["a=1.0, b=1.0"]) and pos_tab["a=1.0, b=1.0"][2] < pos_tab["a=1.0, b=1.0"][0]
          and pos_tab["a=1.0, b=1.5"][0] < 0 and pos_tab["a=1.0, b=1.5"][2] < 3*pos_tab["a=1.0, b=1.5"][1] < 0)
ZAHLEN["lambda_zeitwahl"] = lam_tab
ZAHLEN["positivitaet_min_eig_D1 (N = 200, 800, 3200)"] = pos_tab
check(f"[E] Lambda als Funktion der Zeitwahl (sympy): Für Q_0 = a E + b Q^-, also Q_mu = ((a+b) P_mu + (a-b) K_mu)/2, "
      f"lassen die vier Q_mu genau Omega = 1 + ((a-b)/(a+b)) det x invariant, und der Ricci-Skalar gibt Lambda = "
      f"-12 (a-b)/(a+b) in Einheiten Im x_Vak = 1: " + ", ".join(f"({k}) -> {v['Lambda']}" for k, v in lam_tab.items())
      + f"; sign Lambda = sign tr(Q_0^2): elliptisch -> AdS, parabolisch -> flach, hyperbolisch -> dS. Positivität auf "
      f"D_1: inf spec(a L_0 + b H_1) = h sqrt(a^2 - b^2) für a > |b| ({pos_tab['a=1.0, b=0.5'][2]:.6f}), -> 0+ für a = b "
      f"({pos_tab['a=1.0, b=1.0']}), unbeschränkt nach unten für |b| > a ({pos_tab['a=1.0, b=1.5']}): positive Zeit <=> "
      f"Lambda <= 0", ok_lam and ok_pos)
laufzeit()

# --- 2.4 Gibbons-Hawking-Temperatur in Nahteinheiten und die Beziehung zu beta = 4 pi
uu = sp.Symbol("u", real=True)
tc_u = 2*sp.atan(sp.tanh(uu/2))                                         # Uhrzeit auf der Achse
dtc = sp.simplify(sp.diff(tc_u, uu))
ok_t = sp.simplify(dtc - 1/sp.cosh(uu)) == 0 and dtc.subs(uu, 0) == 1
beta_uhr = 2*sp.pi*dtc.subs(uu, 0)
beta_eigen = sp.pi                                                       # tau_p = u/2
verh = {"beta_Original/beta_dS (Uhr)": sp.nsimplify(4*sp.pi/beta_uhr), "beta'/beta_dS (Uhr)": sp.nsimplify(8*sp.pi/beta_uhr)}
VE, VQ = vektor6(MEt), vektor6(MQ)
per = {"exp(2 pi V_E) = 1 (reelle Periode der Uhr)": float(np.abs(sla.expm(2*np.pi*VE) - np.eye(6)).max()),
       "exp(2 pi i V_Q) = 1 (imaginäre Periode der Kappe)": float(np.abs(sla.expm(2j*np.pi*VQ) - np.eye(6)).max()),
       "exp(4 pi i V_E) - 1 (beta_Original imaginär)": float(np.abs(sla.expm(4j*np.pi*VE) - np.eye(6)).max()),
       "exp(4 pi V_Q) - 1 (reell)": float(np.abs(sla.expm(4*np.pi*VQ) - np.eye(6)).max()),
       "Twistor exp(2 pi M_E) = -1": float(np.abs(sla.expm(2*np.pi*MEt) + np.eye(4)).max())}
wick = sla.expm(0.5j*np.pi*MD) @ MEt @ sla.expm(-0.5j*np.pi*MD)
ok_per = (per["exp(2 pi V_E) = 1 (reelle Periode der Uhr)"] < 1e-10 and per["exp(2 pi i V_Q) = 1 (imaginäre Periode der Kappe)"] < 1e-10
          and per["exp(4 pi i V_E) - 1 (beta_Original imaginär)"] > 1e4 and per["exp(4 pi V_Q) - 1 (reell)"] > 1e4
          and per["Twistor exp(2 pi M_E) = -1"] < 1e-10 and np.allclose(wick, 1j*MQ)
          and np.allclose(sorted(np.linalg.eigvals(VQ).real), [-1, 0, 0, 0, 0, 1]))
th4 = sp.symbols("th1:4", positive=True)
ph_, l_s = sp.symbols("phi l", positive=True)
gS4 = [l_s**2, l_s**2*sp.sin(th4[0])**2, l_s**2*sp.sin(th4[0])**2*sp.sin(th4[1])**2,
       l_s**2*sp.sin(th4[0])**2*sp.sin(th4[1])**2*sp.sin(th4[2])**2]
E4_S4, R_S4 = euler_dichte(gS4, [th4[0], th4[1], th4[2], ph_])
vol_S4 = sp.integrate(l_s**4*sp.sin(th4[0])**3*sp.sin(th4[1])**2*sp.sin(th4[2]), (th4[0], 0, sp.pi), (th4[1], 0, sp.pi),
                      (th4[2], 0, sp.pi), (ph_, 0, 2*sp.pi))
chi_S4 = sp.simplify(E4_S4*vol_S4/(32*sp.pi**2))
tau_e, psi_e, the_e = sp.symbols("tau_E psi theta_E", real=True)
E4_H, _ = euler_dichte([sp.Integer(1), sp.Integer(1), sp.sin(psi_e)**2, sp.sin(psi_e)**2*sp.sin(the_e)**2],
                         [tau_e, psi_e, the_e, ph_])
ok_chi = sp.simplify(E4_S4 - 24/l_s**4) == 0 and chi_S4 == 2 and sp.simplify(E4_H) == 0 and sp.simplify(R_S4 - 12/l_s**2) == 0
ZAHLEN["temperatur"] = {"T_GH Vakuumeinheiten (Eigenzeit, l = 1/2)": "1/pi", "beta_dS Eigenzeit": str(beta_eigen),
                        "beta_dS Uhrzeit am Kappenzentrum": str(beta_uhr), "d tau_Uhr/du": str(dtc),
                        **{k: str(v) for k, v in verh.items()}, "Perioden": per,
                        "Wick: Ad(e^(i pi D/2)) E = i Q^-": True,
                        "E_4 auf S^4 (Radius l)": str(E4_S4), "chi(S^4)": str(chi_S4),
                        "E_4 auf S^1 x S^3": str(E4_H), "chi(S^1 x S^3)": 0}
check(f"[E] Temperatur in Nahteinheiten: Am Kappenzentrum ist d tau_Uhr/du = sech u = 1, also beta_dS = 2 pi in Uhrzeit "
      f"(Einheiten von E = L_0, in denen das Original beta = 4 pi hat) und pi in dS-Eigenzeit (T_GH = 1/pi bei l = Im x_Vak/2); "
      f"Verhältnisse beta_Original/beta_dS = {verh['beta_Original/beta_dS (Uhr)']}, beta'/beta_dS = "
      f"{verh['beta' + chr(39) + '/beta_dS (Uhr)']}. Mechanismustest: In der Vektordarstellung ist exp(2 pi V_E) = 1 (Uhr "
      f"periodisch in reeller Zeit) und exp(2 pi i V_Q) = 1 (Kappe periodisch in imaginärer Zeit: beta = 2 pi ist geometrisch), "
      f"aber exp(4 pi i V_E) != 1 (|.| = {per['exp(4 pi i V_E) - 1 (beta_Original imaginär)']:.2e}): beta = 4 pi ist "
      f"keine Periode der Uhr, sondern eine Zustandswahl. Die Wickabbildung Ad(e^(i pi D/2)) E = i Q^- verbindet genau die "
      f"reelle 2 pi der Uhr mit der imaginären 2 pi der Kappe. Euklidische Schnitte (Gauss-Bonnet, sympy): dS/Kappe S^4 mit "
      f"E_4 = {E4_S4}, chi = {chi_S4}; Original S^1_beta x S^3 mit E_4 = {E4_H}, chi = 0",
      ok_t and ok_per and ok_chi)
befund("E: Die Gibbons-Hawking-Temperatur der nativen Kappe ist T = 1/pi in Vakuumeinheiten (Eigenzeit, l = 1/2) bzw. "
       "1/(2 pi) in Uhrzeit am Kappenzentrum. Zu beta = 4 pi gibt es nur eine Zahlbeziehung (Faktor 2, Produktlesart "
       "Faktor 4), keinen Mechanismus: Die 2 pi der Kappe ist die imaginäre Periode eines hyperbolischen Erzeugers (glatter "
       "euklidischer Schnitt S^4, chi = 2, universell nach Hislop-Longo), die 4 pi des Originals ist ein KMS-Datum des "
       "elliptischen Erzeugers ohne imaginäre Periode (euklidischer Schnitt S^1 x S^3, chi = 0, MM-Eulerterm null). "
       "Native Wick-Dualität gibt es nur zwischen der reellen Uhrperiode 2 pi und der Kappenperiode 2 pi", True)
laufzeit()

# =================================================================== Teil 3
say("\nTeil 3: Urteil H5")
# --- 3.1 Auswahltabelle der Zeitkandidaten
KAND = {"E": (MEt, Eop), "P_0": (MP[0], Pop[0]), "Q^-": (MQ, Qm), "D": (MD, Dil)}
L0s, H1s, H2s = dh_matrizen(1, 1600)
D1_ops = {"E": L0s, "P_0": L0s + H1s, "Q^-": H1s, "D": H2s}               # s-Welle D_1 der radialen sl(2,R)
tab = {}
for nm, (M, X) in KAND.items():
    mn = float(np.linalg.eigvalsh(D1_ops[nm].toarray()).min())
    positiv = mn > -1e-9
    boost = np.allclose(M @ M, np.eye(4)/4)                              # diagonalisierbar, Eigenwerte +-1/2
    ebs = eigenebenen(M) if boost else []
    modular = (boost and len(ebs) == 2 and all(max(abs(np.array(g))) < 1e-9 for _, g, _ in ebs)
               and klein(Par @ X @ Par + X))
    mm_zeit = im_spann(M, ungerade)
    uhr_orig = im_spann(M, [MEt])
    if mm_zeit:
        c_ = koeff(M, [MEt, MQ])
        lam_v = -12*(c_[0] - c_[1])/(c_[0] + c_[1]) if abs(c_[0] + c_[1]) > 1e-12 else None
    else:
        lam_v = None
    tab[nm] = {"min spec auf D_1 (N = 1600)": round(mn, 6), "positiv": positiv, "Kappen-Modulgenerator": modular,
               "MM-Zeit (tau-ungerade, so(3)-invariant)": mm_zeit, "Uhr des Originals": uhr_orig,
               "Lambda (Im x_Vak = 1)": None if lam_v is None else round(float(lam_v), 10) + 0.0}
auswahl = {
    "MM-Zeit und Kappenmodul": [n for n, v in tab.items() if v["MM-Zeit (tau-ungerade, so(3)-invariant)"] and v["Kappen-Modulgenerator"]],
    "MM-Zeit und positiv": [n for n, v in tab.items() if v["MM-Zeit (tau-ungerade, so(3)-invariant)"] and v["positiv"]],
    "MM-Zeit und Uhr des Originals": [n for n, v in tab.items() if v["MM-Zeit (tau-ungerade, so(3)-invariant)"] and v["Uhr des Originals"]],
    "positiv und Kappenmodul": [n for n, v in tab.items() if v["positiv"] and v["Kappen-Modulgenerator"]]}
ok_tab = (auswahl["MM-Zeit und Kappenmodul"] == ["Q^-"] and auswahl["MM-Zeit und positiv"] == ["E", "P_0"]
          and auswahl["MM-Zeit und Uhr des Originals"] == ["E"] and auswahl["positiv und Kappenmodul"] == []
          and tab["Q^-"]["Lambda (Im x_Vak = 1)"] == 12 and tab["E"]["Lambda (Im x_Vak = 1)"] == -12
          and tab["P_0"]["Lambda (Im x_Vak = 1)"] == 0)
ZAHLEN["auswahltabelle"] = tab
ZAHLEN["auswahl"] = auswahl
for nm, v in tab.items():
    say(f"   {nm:4s}: positiv {str(v['positiv']):5s} | Kappenmodul {str(v['Kappen-Modulgenerator']):5s} | MM-Zeit "
        f"{str(v['MM-Zeit (tau-ungerade, so(3)-invariant)']):5s} | Uhr des Originals {str(v['Uhr des Originals']):5s} | "
        f"Lambda {v['Lambda (Im x_Vak = 1)']}")
check(f"[E] Auswahltabelle der Zeitkandidaten (positiv = Spektrum auf D_1 nach unten beschränkt; Kappenmodul = normierter "
      f"Boost M^2 = 1/4 mit Null-Eigenebenen, also reellen Spitzen, und spektral umkehrbar; MM-Zeit = tau-ungerade und so(3)-invariant, also Q_0 einer "
      f"Lorentz-kovarianten Familie; Uhr des Originals = proportional zu E): MM-Zeit und Kappenmodul -> "
      f"{auswahl['MM-Zeit und Kappenmodul']} (Lambda = +12/R_Vak^2), MM-Zeit und positiv -> {auswahl['MM-Zeit und positiv']} "
      f"(Lambda = -12/R_Vak^2 bzw. 0), MM-Zeit und Uhr des Originals -> {auswahl['MM-Zeit und Uhr des Originals']} "
      f"(AdS), positiv und Kappenmodul -> {auswahl['positiv und Kappenmodul']} (leer: Modulgeneratoren sind spektral "
      f"symmetrisch, J log Delta J = -log Delta, also nie positiv; für die Kappe, Typ III_1, ist das Spektrum ganz R). H5 ist damit genau die Prämisse P_cap: 'Q_0 ist die Modulzeit der nativen "
      f"Kappe im Vakuumsektor'", ok_tab)

# --- 3.2 H5-Kinematik und das Lambda-G-Ziel (bedingte Identitäten, v60 nur gelesen)
Rv, Gs, al = sp.symbols("R_Vak G alpha", positive=True)
l_h5 = Rv/2
Lam_h5 = 3/l_h5**2
sig_h5 = -1
b_h5 = l_h5**2/(64*sp.pi*Gs*sig_h5)
S_euler = sp.simplify(abs(b_h5)*32*sp.pi**2*chi_S4)
G_ziel = sp.exp(-2/al)/(128*sp.pi**3)*Rv**2
S_ziel = sp.simplify((sp.pi*l_h5**2/Gs).subs(Gs, G_ziel))
l2G_ziel = sp.simplify((l_h5**2/Gs).subs(Gs, G_ziel))
lg_e = 2*ALPHA_INV/mp.log(10)
zahl = {"log10 G/R_Vak^2": float(-lg_e - mp.log10(128*mp.pi**3)), "log10 S_dS": float(mp.log10(32*mp.pi**4) + lg_e),
        "log10 l^2/G": float(mp.log10(32*mp.pi**3) + lg_e)}
ok_h5 = (sp.simplify(Lam_h5 - 12/Rv**2) == 0 and sp.simplify(S_euler - sp.pi*l_h5**2/Gs) == 0
         and sp.simplify(S_ziel - 32*sp.pi**4*sp.exp(2/al)) == 0 and sp.simplify(l2G_ziel - 32*sp.pi**3*sp.exp(2/al)) == 0
         and b_h5.subs({Rv: 1, Gs: 1}) < 0)
ZAHLEN["h5_kinematik"] = {"Lambda": str(Lam_h5), "l": str(l_h5), "b (sigma = -1)": str(b_h5),
                          "|b| 32 pi^2 chi(S^4)": str(S_euler), "Ziel G/R_Vak^2": "e^(-2/alpha)/(128 pi^3)",
                          "Ziel l^2/G": str(l2G_ziel), "Ziel S_dS = pi l^2/G": str(S_ziel), **zahl}
check(f"[E bedingt] H5-Kinematik: Mit Q_0 = Q^- (P_cap) ist die MM-Algebra so(4,1), Lambda = 3/l^2 = {Lam_h5} > 0 mit "
      f"l = R_Vak/2 (Spiegelbild von z98s -12/R_Vak^2), b = l^2/(64 pi G sigma) = {b_h5} < 0, und das Eulergewicht "
      f"|b| 32 pi^2 chi(S^4) = {S_euler} ist die dS-Entropie (glatter Schnitt S^4 der Kappe, Teil 2). Das Lambda-G-Ziel "
      f"(v60 nur gelesen) bleibt eine bedingte Identität: G/R_Vak^2 = e^(-2/alpha)/(128 pi^3) (10^{zahl['log10 G/R_Vak^2']:.2f}) "
      f"<=> l^2/G = {l2G_ziel} <=> S_dS = {S_ziel} (10^{zahl['log10 S_dS']:.2f}); H5 behebt nur das Vorzeichen, G bleibt "
      f"nativ unbestimmt (z90)", ok_h5)
laufzeit()

# --- 3.3 Bilanz und Urteil
STATUS = {
    "Kappengenerator": ("E", "K_cap = 2 pi Q^-, Q^- = (1/2) Summe (d^+ e^+ + e d) = Original-Partner-Quetschung, "
                             "hyperbolisch in der radialen sl(2,R) {E, Q^-, D}; Spitzen t = +-1 (Null-Eigenebenen)"),
    "native Kappe": ("E/C", "O_1 aus (E, tau) eindeutig bis auf Orientierung (Killing-orthogonal zu E, tau-ungerade); "
                            "Marken = Spitzen von O_1 und V_+- unter z83s Identifikation [C]; E allein: keine Kappe"),
    "Modulfluss nativ?": ("C/E negativ", "Vakuumsektor ja (Hislop-Longo nach dem Netz, KMS 2 pi exakt auf der Achse); "
                                        "Original nein (zentral/periodisch; nicht invariant, Abweichung ~e^(-8 pi))"),
    "Beziehungen": ("E", "K_cap = 2 pi (P_0 - E) = 2 pi P_0 - K_Uhr/2; P_0 = lim K_(O_R)/(pi R) (flacher Grenzfall)"),
    "dS-Darstellung": ("E", "so(4,1) auf dem Skalarsektor: Komplementärserie Delta = 1 (m^2 l^2 = 2, xi = 1/6), "
                            "irreduzibel (K-Typen (j, j) je einmal); Casimir = AdS-Casimir"),
    "Energiespektrum": ("E", "Q^- je l: D_(l+1), sech-Gesetz, Dichte 2^(2h)|Gamma(h + i lam)|^2/(2 pi Gamma(2h)), R symmetrisch"),
    "Lambda und Zeitklasse": ("E", "Lambda = -12 (a-b)/(a+b) für Q_0 = a E + b Q^-; positive Zeit <=> Lambda <= 0"),
    "Temperatur": ("E/L", "T_GH = 1/pi (Vakuumeinheiten), beta = 2 pi in Uhrzeit am Zentrum; zu beta = 4 pi nur Faktor 2, "
                          "kein Mechanismus (S^4 gegen S^1 x S^3)"),
    "H5": ("C", "bedingt mit P_cap (Zeit = Vakuum-Modulfluss der nativen Kappe), nicht nativ; Positivität und die Uhr "
                "des Originals wählen Lambda <= 0"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
urteil = ("H5 hält nur bedingt, nicht nativ. (1) Die Kappe ist nativ konstruierbar: Q_0^- = (P_0 - K_0)/2 ist auf dem "
          "Fockraum die Zweimoden-Quetschung (1/2) Summe (d^+ e^+ + e d) von Original und Partner, das hyperbolische Element "
          "der radialen sl(2,R) {E, Q^-, D}; seine Eigenebenen auf dem Twistorraum sind null und gehören zu den reellen Spitzen "
          "t = +-1, sein Feld ((1 - t^2 - |x|^2)/2, -t x) hat den Rand des Doppelkegels O_1 als Killinghorizont, und K_cap = "
          "2 pi Q^- ist die Hislop-Longo-Normierung. O_1 ist durch die nativen Daten E und tau eindeutig bis auf die "
          "Orientierung bestimmt (tau-ungerade und Killing-orthogonal zu E in der radialen sl(2,R)): die Halbsphäre der S^3 "
          "zur Uhrzeit 0 um den tau-Fixpunkt, mit den Spitzen t = +-Im x_Vak als reellen Wick-Partnern der Vakuumpunkte. "
          "E allein gibt keine Kappe (definite Eigenebenen, kein Horizont); unter z83s Identifikation sind die vier Marken "
          "z^4 = 1 genau die Spitzen der beiden tau-Eigenkappen O_1 und V_+- [C]. (2) Ihr Modulfluss ist nativ nur im "
          "Vakuumsektor: Mit dem freien Netz des Leiterfelds (z83/z90, Hislop-Longo in der erlaubten Reihenfolge) ist die "
          "Achsen-Zweipunktfunktion im Kappenrahmen exakt 1/(16 pi^2 sinh^2(du/2)), stationär und KMS bei beta = 2 pi; das "
          "Original dagegen hat einen zentralen, periodischen Modulgenerator (Ladung N_d - N_e bzw. einseitig die Uhr 4 pi E) "
          "und ist unter dem Kappenfluss nicht invariant (Abweichung der Ordnung e^(-8 pi) in der Produktlesart), also für ihn "
          "bei keinem beta KMS. Exakt gilt K_cap = 2 pi (P_0 - E) = 2 pi P_0 - K_Uhr/2, und z72s parabolische Zeit ist der "
          "flache Grenzfall P_0 = lim K_(O_R)/(pi R). (3) dS-Kinematik: span{Q^-_a, M_ab} wirkt auf dem Skalarsektor der "
          "Leiter als irreduzible unitäre Komplementärserie von so(4,1) mit Delta = 1 (kanonischer Casimir gleich dem von "
          "so(3,2), C_phys = -2, m^2 l^2 = 2 = 12 xi, jeder K-Typ (j, j) einmal), also der konform gekoppelte Skalar im "
          "Bunch-Davies-Zustand. Die Kappenenergie hat je Drehimpuls l das Spektrum von D_(l+1): ganz R, absolut stetig, "
          "symmetrisch, mit der Dichte 2^(2h)|Gamma(h + i lambda)|^2/(2 pi Gamma(2h)) und dem sech-Gesetz cosh(s/2)^(-2h); "
          "keine Positivität. Die Gibbons-Hawking-Temperatur ist T = 1/pi in Vakuumeinheiten (l = 1/2) bzw. beta = 2 pi in "
          "Uhrzeit am Kappenzentrum; zu beta = 4 pi besteht nur der Zahlfaktor 2 (Produktlesart 4) ohne Mechanismus: Die 2 pi "
          "ist die imaginäre Periode des hyperbolischen Erzeugers (glatter Schnitt S^4, chi = 2), die 4 pi ist ein KMS-Datum "
          "des elliptischen ohne imaginäre Periode (Schnitt S^1 x S^3, chi = 0). (4) Exakt gilt Lambda = -12 (a-b)/(a+b)/"
          "R_Vak^2 für Q_0 = a E + b Q^-: Das Vorzeichen von Lambda ist die Konjugationsklasse der Zeit (elliptisch AdS, "
          "parabolisch flach, hyperbolisch dS), und positive Zeit ist genau Lambda <= 0. Die Auswahltabelle zeigt: "
          "MM-Zeit und Kappenmodul wählen allein Q^- (Lambda = +12/R_Vak^2), MM-Zeit und Positivität wählen E oder P_0, "
          "MM-Zeit und die Uhr des Originals wählen E (AdS); Positivität und Modularität schließen einander aus. H5 ist damit "
          "gleichwertig zur Prämisse P_cap: 'die physikalische Zeit ist der Vakuum-Modulfluss der nativen Kappe O_1', und sie "
          "verlangt den Verzicht auf eine positive Zeit. Folgen: K13 bleibt L, präzisiert um die exakte Zuordnung "
          "Vorzeichen von Lambda = Klasse der Zeit; das Lambda-G-Ziel G/R_Vak^2 = e^(-2/alpha)/(128 pi^3) bleibt unerreicht "
          "(H5 behebt nur das Vorzeichen; mit ihm wäre S_dS = pi l^2/G = 32 pi^4 e^(2/alpha) das Eulergewicht auf S^4, ohne "
          "natives G). Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z104_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
