"""Z98: Phase des Infinity-Twistors, Kopplung g und der MacDowell-Mansouri-Weg (Folgeauftrag G2).

Leitlinie 'Rekonstruktion statt Auswahl'. Firewalls: GRAV.NONCIRCULAR.01 (die Zeit ist der parabolische Erzeuger, nicht der
Modulfluss; wo eine Konstruktion den Modulfluss als Zeit bräuchte, zählt das als Befund), GRAV.PARENT.GAP.NO_GO.01,
DIMENSION.UPLIFT.FIREWALL.01 (alles ist Twistor- und Fockkinematik einer bzw. zweier Nähte, Liealgebra und Formenalgebra;
keine wechselwirkende 3+1-Schließung). Gelöste Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt; nichts wird
befördert. Bestand: z89 (I_pi nativ bis auf die Modulphase, B = m_Z o (eps c^+ (x) c^+), Lorentzbüschel, radialer Kanal
AdS, g frei), z84 (Theta_n als U(2) x U(2)-Bahn von tau, Theta_1 = tau, Reflexionspositivität), z83 (radiale su(2,2)), z85,
z29 (Plebanski-Simplizität und Realität), z65 (Tomita-J ist die Realstruktur PT -> PT*), z90 (kein G aus dem Schnitt),
Ledger GRAV.NONLINEAR.01 und v60 (rho_Lambda/Mbar^4 = (3/(4 pi^2)) e^(-2/alpha), nur gelesen). Geprüft wird die zwölfte
externe Analyse (Q_a = (P_a + sigma K_a)/2, MacDowell-Mansouri, Lambda-G-Relation, zwei Blätter, modulare Berry-Verbindung).
Statuskürzel: E exakt, C bedingt, L Lesart, O offen, N numerisch.

Teil 1  Phase von I und Kopplung g: Modul-U(1) gegen die nativen Daten (Fock 8^4, zwei Nähte 4^4 x 4^4), Realitätsbedingungen
        für e^(i phi) I_pi (tau linear, Tomita-Theta auf zwei Nähten, lorentzsche Realität mit Volumenform, euklidische
        quaternionische Struktur), Entfernbarkeit der Phase auf Gravitonzuständen, Reflexionspositivität der
        Austauschamplitude zweiter Ordnung und g.
Teil 2  MacDowell-Mansouri nativ: Q^sigma = (P + sigma K)/2 aus der nativen radialen su(2,2), Vergleich mit
        Der(I_pi - sigma I_omega), natives sigma, Radius aus dem invarianten konformen Faktor und der Krümmung (sympy), tau als
        Punktspiegelung von AdS4 (symmetrisches Paar), modulare Berry-Verbindung (Kozyklus flach, Modulgenerator zentral,
        Tomita-Familie = Levi-Civita von AdS4 mit MM-Krümmung null), Plebanski-Zweiformen aus zwei Blättern und Phasentest.
Teil 3  Skalen: MacDowell-Mansouri-Zerlegung (Koeffizienten auf Formen), bedingte Identität l^2/G = 32 pi^3 e^(2/alpha) aus v60,
        Vorzeichen und Größenordnung gegen den nativen Radius.
Teil 4  Bilanz und Urteil für G2.
"""
import itertools
import json
import math
import time
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.linalg as sla
import scipy.sparse as sps
import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
rng = np.random.default_rng(98)
mp.mp.dps = 30
I2 = np.eye(2)
SIGMA = [np.eye(2, dtype=complex), np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]),
         np.diag([1., -1.]).astype(complex)]
EPS2 = np.array([[0., 1.], [-1., 0.]], dtype=complex)                  # eps^(jk) mit eps^(12) = 1
KAP2_NATIV = {m: float(mp.tanh(mp.pi*m)) for m in (1, 2, 3)}           # kappa^2 = tanh(pi m) (z85, z89)
ALPHA_INV = mp.mpf("137.035999177")                                     # CODATA 2022, nur zur Größenordnung


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


# ------------------------------------------------------------------ Fockraum einer Naht (wie z85, z89)
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
    """Vertauschung d <-> e auf dem abgeschnittenen Fockraum (reell); J = Tausch mal komplexe Konjugation (z89)."""
    occ = np.array(list(itertools.product(range(nv), repeat=4)))
    dim = nv**4
    idx = [np.ravel_multi_index((o[2], o[3], o[0], o[1]), (nv,)*4) for o in occ]
    return sps.csr_matrix((np.ones(dim), (idx, np.arange(dim))), shape=(dim, dim)), occ


def J_ad(X, SWm):
    """J X J^-1 = SW conj(X) SW (J antiunitär, die Fockbasis ist J-reell bis auf den Tausch)."""
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
CC = [(DD[i] + dag(EE[i])).tocsr() for i in range(2)]             # c = d + e^+ (pi-Ebene, tau = +1)
WO = [(DD[i] - dag(EE[i])).tocsr() for i in range(2)]             # omega'' = d - e^+ (tau = -1)
Pop = [0.5*sum(s[i, j]*(dag(CC[i]) @ CC[j]) for i in range(2) for j in range(2)) for s in SIGMA]
Kop = [0.5*sum(s[i, j]*(dag(WO[i]) @ WO[j]) for i in range(2) for j in range(2)) for s in SIGMA]
Eop = 0.5*sum((dag(X) @ X for X in DD + EE), sps.csr_matrix((DIMV, DIMV))) + sps.identity(DIMV)
Dop = sps.diags(DVAL.astype(float)).tocsr()


def tau_op(X):
    return dag(J_ad(X, SW))


# ------------------------------------------------------------------ u(2,2) auf dem Twistorraum (z-Koordinaten, wie z89)
G4 = np.diag([1., 1., -1., -1.]).astype(complex)
S4 = np.block([[np.zeros((2, 2)), I2], [I2, np.zeros((2, 2))]]).astype(complex)   # tau = (Ad J)^+ auf Z
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


def derivationen(I):
    return teilalgebra([lambda M: M @ I + I @ M.T])


def im_spann(M, basis):
    if not basis:
        return np.abs(M).max() < 1e-12
    A = np.array([np.concatenate([B.real.ravel(), B.imag.ravel()]) for B in basis]).T
    b = np.concatenate([M.real.ravel(), M.imag.ravel()])
    x = np.linalg.lstsq(A, b, rcond=None)[0]
    return np.linalg.norm(A @ x - b) < 1e-9*max(1.0, np.linalg.norm(b))


def signatur(basis):
    Gm = np.array([[np.trace(X @ Y).real for Y in basis] for X in basis])
    ew = np.linalg.eigvalsh((Gm + Gm.T)/2)
    sk = max(1.0, np.abs(ew).max())
    return int(np.sum(ew > 1e-9*sk)), int(np.sum(ew < -1e-9*sk)), int(np.sum(np.abs(ew) <= 1e-9*sk))


uvec = [np.array([1, 0, -1, 0.]), np.array([0, 1, 0, -1.])]
vvec = [np.array([1, 0, 1, 0.]), np.array([0, 1, 0, 1.])]
Inp = (np.outer(uvec[0], uvec[1]) - np.outer(uvec[1], uvec[0])).astype(complex)       # I_pi (z89)
Ion = (np.outer(vvec[0], vvec[1]) - np.outer(vvec[1], vvec[0])).astype(complex)       # I_omega (z89)
NP = [-0.5*sum(s[i, j]*np.outer(uvec[i], vvec[j]) for i in range(2) for j in range(2)) for s in SIGMA]
NK = [-0.5*sum(s[i, j]*np.outer(vvec[i], uvec[j]) for i in range(2) for j in range(2)) for s in SIGMA]
NE = np.diag([-0.5, -0.5, 0.5, 0.5]).astype(complex)
MP = [1j*N for N in NP]
MK = [1j*N for N in NK]
MEt = 1j*NE
LOR = teilalgebra([lambda M: M @ S4 - S4 @ M, lambda M: M @ Inp + Inp @ M.T])
LC = np.zeros((4,)*4)
for _p in itertools.permutations(range(4)):
    LC[_p] = round(np.linalg.det(np.eye(4)[list(_p)]))

# =================================================================== Teil 1
say("Teil 1: Phase von I und Kopplung g")

# --- 1.1 Modul-U(1) gegen die nativen Daten (eine Naht, Fock 8^4) und Ladung des eps-Kerns (zwei Nähte)
th = 0.7
U = sps.diags(np.exp(1j*th*DVAL)).tocsr()
ok_JU = klein(J_ad(U, SW) - U, 1e-13)
ok_vak = np.allclose(U @ OMV, OMV)
ok_D = all(klein(komm(Dop, X)) for X in Pop + Kop + [Eop])
ok_cJ = all(klein(J_ad(CC[j], SW) - dag(CC[j])) and klein(J_ad(WO[j], SW) + dag(WO[j])) for j in range(2))
perm_tau = [2, 3, 0, 1]
ok_tau = all(klein(tau_op(ZT[k]) - ZT[perm_tau[k]]) for k in range(4))
ok_tau_lin = all(klein(tau_op(np.exp(0.4j)*ZT[k]) - np.exp(0.4j)*ZT[perm_tau[k]]) for k in range(4))
ok_tauU = all(klein(tau_op(U @ ZT[k] @ dag(U)) - U @ tau_op(ZT[k]) @ dag(U)) for k in range(4))
NV2 = 4
q1, q2, f1, f2 = fock4(NV2)
SW4, OCC4 = tausch(NV2)
D4 = OCC4[:, 0] + OCC4[:, 1] - OCC4[:, 2] - OCC4[:, 3]
C4 = [(q1 + dag(f1)).tocsr(), (q2 + dag(f2)).tocsr()]
V0 = (sps.kron(dag(C4[0]), dag(C4[1])) - sps.kron(dag(C4[1]), dag(C4[0]))).tocsr()   # eps^jk c_j^(1)+ c_k^(2)+
SW2 = sps.kron(SW4, SW4).tocsr()
Dtot = (D4[:, None] + D4[None, :]).ravel()
U2 = sps.diags(np.exp(1j*th*Dtot)).tocsr()
res_ladung = float(abs(U2 @ V0 @ dag(U2) - np.exp(2j*th)*V0).max())
PHIS = [0.0, 0.3, 1.1, np.pi/2, 2.5]
res_theta = max(float(abs(J_ad(np.exp(1j*ph)*V0, SW2) - dag(np.exp(1j*ph)*V0)).max()) for ph in PHIS)
ZAHLEN["modul_u1"] = {"theta": th, "[D, P_mu], [D, K_mu], [D, E]": 0.0, "Ladung eps-Kern (U V U^+ = e^(2 i theta) V), Residuum":
                      res_ladung, "Theta-Realität J V_phi J = V_phi^+ (max über phi)": res_theta,
                      "nnz eps-Kern": int(V0.nnz)}
check(f"Modul-U(1) gegen die nativen Daten (Fock 8^4): U_theta = e^(i theta D) vertauscht mit der Tomita-Konjugation "
      f"(J U J^-1 = U, J antiunitär), lässt das KMS-Vakuum fest und vertauscht mit P_mu, K_mu und E; tau = (Ad J)^+ ist "
      f"auf T = span{{Z}} linear (tau(e^(i phi) Z) = e^(i phi) tau(Z)) und vertauscht mit U_theta; J c J = c^+, "
      f"J omega'' J = -omega''^+. Zwei Nähte (4^4 x 4^4): der eps-Kern V = eps^jk c_j^(1)+ c_k^(2)+ hat die Modulladung 2 "
      f"(Residuum {res_ladung:.1e}), und die Theta-Realität J V_phi J^-1 = V_phi^+ gilt für jede Phase phi "
      f"(fünf Phasen, Residuum {res_theta:.1e}): Die native Reflexion Theta_1 = tau wählt keine Phase",
      ok_JU and ok_vak and ok_D and ok_cJ and ok_tau and ok_tau_lin and ok_tauU and res_ladung < 1e-12
      and res_theta < 1e-12 and V0.nnz > 0)

# --- 1.2 Realitätsbedingungen für I(phi) = e^(i phi) I_pi auf dem Twistorraum
ok_tau_biv = all(np.allclose(S4 @ (np.exp(1j*ph)*Inp) @ S4.T, np.exp(1j*ph)*Inp) for ph in np.linspace(0, 2*np.pi, 9))


def stern_unten(I, psi):
    """(*I)_ab = (1/2) eps_abcd I^cd mit der Volumenform eps_1234 = e^(i psi)."""
    return 0.5*np.exp(1j*psi)*np.einsum("abcd,cd->ab", LC, I)


def quer_unten(I):
    """I-quer_ab = G_aa' conj(I^a'b') G_b'b (hermitesche Twistorform, Konjugation T -> T*)."""
    return G4 @ I.conj() @ G4


lam0 = np.vdot(quer_unten(Inp).ravel(), stern_unten(Inp, 0).ravel())/np.vdot(quer_unten(Inp).ravel(),
                                                                             quer_unten(Inp).ravel())
ok_prop = np.allclose(stern_unten(Inp, 0), lam0*quer_unten(Inp)) and abs(abs(lam0) - 1) < 1e-12
GRID = np.linspace(0, 2*np.pi, 1441)[:-1]
lor_tab = {}
ok_lor = ok_prop
for psi in (0.0, 0.4, 1.3):
    res = np.array([np.abs(quer_unten(np.exp(1j*ph)*Inp) - stern_unten(np.exp(1j*ph)*Inp, psi)).max() for ph in GRID])
    ph_s = float(np.mod(-(psi + np.angle(lam0))/2, np.pi))
    r_s = [np.abs(quer_unten(np.exp(1j*p)*Inp) - stern_unten(np.exp(1j*p)*Inp, psi)).max() for p in (ph_s, ph_s + np.pi)]
    fern = np.abs(np.angle(np.exp(2j*(GRID - ph_s)))) > 0.1
    ok_lor &= max(r_s) < 1e-12 and res[fern].min() > 0.05
    lor_tab[f"psi={psi}"] = {"reelle Phase phi* (mod pi)": ph_s, "phi* + psi/2 (mod pi)": float(np.mod(ph_s + psi/2, np.pi))}
konst = {round(v["phi* + psi/2 (mod pi)"], 12) for v in lor_tab.values()}
ok_lor &= len(konst) == 1
det_zentrum = complex(np.linalg.det(np.exp(-1j*th)*np.eye(4)))
ok_lor &= abs(det_zentrum - np.exp(-4j*th)) < 1e-12
Mq = sla.block_diag(np.array([[0, -1], [1, 0]]), np.array([[0, -1], [1, 0]])).astype(complex)
quat_tab = {}
ok_q = (np.allclose(Mq @ Mq.conj(), -np.eye(4)) and np.allclose(Mq @ S4, S4 @ Mq) and np.allclose(Mq.conj().T @ G4 @ Mq, G4)
        and np.allclose(Mq @ Inp @ Mq.T, Inp))
for thq in (0.0, 0.2, 0.55):
    sq = lambda I: np.exp(-4j*thq)*(Mq @ I.conj() @ Mq.T)                # sigma_theta = e^(-2 i theta) Mq conj auf T
    res = np.array([np.abs(sq(np.exp(1j*ph)*Inp) - np.exp(1j*ph)*Inp).max() for ph in GRID])
    ph_s = float(np.mod(-2*thq, np.pi))
    r_s = max(np.abs(sq(np.exp(1j*p)*Inp) - np.exp(1j*p)*Inp).max() for p in (ph_s, ph_s + np.pi))
    fern = np.abs(np.angle(np.exp(2j*(GRID - ph_s)))) > 0.1
    ok_q &= r_s < 1e-12 and res[fern].min() > 0.05
    quat_tab[f"theta={thq}"] = ph_s
# antilineares tau: tau K_theta mit der Fock-Konjugation K im U(1)_D-Rahmen theta (K_theta = e^(-2 i theta) conj auf T)
tauK_tab = {}
ok_tK = True
for thq in (0.0, 0.3):
    tk = lambda I: np.exp(-4j*thq)*(S4 @ I.conj() @ S4.T)
    res = np.array([np.abs(tk(np.exp(1j*ph)*Inp) - np.exp(1j*ph)*Inp).max() for ph in GRID])
    ph_s = float(np.mod(-2*thq, np.pi))
    fern = np.abs(np.angle(np.exp(2j*(GRID - ph_s)))) > 0.1
    ok_tK &= max(np.abs(tk(np.exp(1j*p)*Inp) - np.exp(1j*p)*Inp).max() for p in (ph_s, ph_s + np.pi)) < 1e-12 \
        and res[fern].min() > 0.05
    tauK_tab[f"theta={thq}"] = ph_s
ZAHLEN["realitaet"] = {"tau (linear)": "jede Phase", "lorentzsch: *I_pi = lambda0 I_pi-quer, lambda0": str(np.round(lam0, 12)),
                       "lorentzsch: reelle Phase je Volumenformphase psi": lor_tab,
                       "det(e^(-i theta) 1) = e^(-4 i theta)": True,
                       "quaternionisch sigma_theta: reelle Phase phi* = -2 theta (mod pi)": quat_tab,
                       "antilineares tau K_theta: reelle Phase": tauK_tab}
check(f"Realitätsbedingungen für I(phi) = e^(i phi) I_pi (Twistorraum, Gitter 0,25 Grad): (a) tau = Theta_1 hält I(phi) "
      f"für alle phi fest, und jede Theta_n ist C-linear (z84), ihre Fixbivektoren bilden komplexe Unterräume (kein "
      f"Phasendatum); die antilineare Fassung tau K (K = Konjugation der Fockbasis) fixiert phi mod pi, aber K hängt vom "
      f"U(1)_D-Rahmen ab (K_theta = U_(2 theta) K) und fixiert dann phi* = -2 theta ({tauK_tab}); (b) die lorentzsche "
      f"Realität I-quer = *I "
      f"(Konjugation über die CCR-Form G, Hodge mit Volumenform e^(i psi)) gilt genau für phi = phi*(psi) mod pi, mit "
      f"phi* + psi/2 = {sorted(konst)[0]:.6f} mod pi für psi = 0, 0,4, 1,3: Sie fixiert die Phase relativ zur Phase der "
      f"Volumenform, und die trägt unter dem Zentrum die Ladung 4 (det(e^(-i theta) 1) = e^(-4 i theta)); (c) die "
      f"euklidische quaternionische Struktur sigma (sigma^2 = -1, vertauscht mit tau, G-verträglich) fixiert phi mod pi, "
      f"aber ihre Modulbahn sigma_theta = e^(-2 i theta) sigma hat dieselben Eigenschaften und fixiert phi* = -2 theta "
      f"({ {k: round(v, 4) for k, v in quat_tab.items()} })",
      ok_tau_biv and ok_lor and ok_q and ok_tK)
befund("E: Theta_1 = tau fixiert die Modulphase von I nicht. Jede native Struktur (KMS-Zustand, J, Delta^(it), tau, "
       "Theta_n, P_mu, K_mu, E, CCR-Form) ist Modul-U(1)-invariant, I_pi trägt die Ladung 2; nach dem Curie-Prinzip "
       "legt kein natives Datum die Phase fest. Realstrukturen, die sie fixieren würden (lorentzsch mit Volumenform, "
       "euklidisch quaternionisch), kommen in Modul-U(1)-Bahnen; die Phase von I ist gleichwertig zur Phase der "
       "Volumenform eps_abcd, also zur Reduktion U(2,2) -> SU(2,2), deren Rest Z_4 nur I -> +-I lässt", True)

# --- 1.3 Entfernbarkeit der Phase auf Gravitonzuständen (Fockwörterbuch z85/z89, exakt)
Z = sp.symbols("z1:5")


def bivektor(r1, r2):
    return r1*r2.T - r2*r1.T


U_PI = [sp.Matrix([1, 0, -1, 0]), sp.Matrix([0, 1, 0, -1])]
I_PI = bivektor(*U_PI)


def grad(f):
    return sp.Matrix([sp.diff(f, z) for z in Z])


def klammer(I, f, g):
    return sp.expand((grad(f).T*I*grad(g))[0, 0])


def norm_fak(N):
    n1, n2, m1, m2 = N
    return sp.sqrt(sp.Rational(math.factorial(n1)*math.factorial(n2), math.factorial(m1)*math.factorial(m2)))


def elementar(N):
    return norm_fak(N)*Z[2]**N[2]*Z[3]**N[3]/(Z[0]**(N[0] + 1)*Z[1]**(N[1] + 1))


def fock_von(expr):
    out = {}
    for t in sp.Add.make_args(sp.expand(expr)):
        if t == 0:
            continue
        c, ex = t, []
        for z in Z:
            c, e = c.as_coeff_exponent(z)
            ex.append(int(e))
        if ex[0] >= 0 or ex[1] >= 0:
            continue
        N = (-ex[0] - 1, -ex[1] - 1, ex[2], ex[3])
        out[N] = sp.expand(out.get(N, 0) + c/norm_fak(N))
    return {k: v for k, v in out.items() if v != 0}


def D_von(N):
    return N[0] + N[1] - N[2] - N[3]


def E_von(N):
    return sp.Rational(sum(N), 2) + 1


N5 = [(0, 0, m, 4 - m) for m in range(5)]                          # Gravitonen tiefster Energie E = 3, D = -4
PAARE = [(a, b) for a in range(5) for b in range(a + 1, 5)]
AUS = {p: fock_von(klammer(I_PI, elementar(N5[p[0]]), elementar(N5[p[1]]))) for p in PAARE}
KEYS = sorted({k for v in AUS.values() for k in v})
A_amp = np.array([[complex(AUS[p].get(k, 0)) for p in PAARE] for k in KEYS])
dsh = {D_von(k) - D_von(N5[a]) - D_von(N5[b]) for (a, b), v in AUS.items() for k in v}
Evals = sorted({E_von(k) for k in KEYS})
thr = 0.37
faktor = np.array([[np.exp(1j*thr*(D_von(k) - D_von(N5[a]) - D_von(N5[b]))) for (a, b) in PAARE] for k in KEYS])
ok_uni = np.allclose(faktor*A_amp, np.exp(4j*thr)*A_amp) and dsh == {4}
GramA = A_amp.conj().T @ A_amp
ok_inv = all(np.allclose((np.exp(1j*ph)*A_amp).conj().T @ (np.exp(1j*ph)*A_amp), GramA) for ph in PHIS)
ZAHLEN["entfernbarkeit"] = {"Paare (E = 3, D = -4)": len(PAARE), "Ausgangskomponenten": len(KEYS),
                            "D-Verschiebung": sorted(dsh), "E der Ausgänge": [str(x) for x in Evals],
                            "Rang der Amplitudenmatrix": int(np.linalg.matrix_rank(A_amp))}
check(f"Entfernbarkeit (exakt, {len(PAARE)} Gravitonpaare E = 3, {len(KEYS)} Ausgangskomponenten bei E = "
      f"{[str(x) for x in Evals]}): Die Konjugation mit U_theta auf allen drei Beinen multipliziert jede Amplitude "
      f"<k|B_I(h_a, h_b)> mit demselben Faktor e^(4 i theta) (D-Verschiebung {sorted(dsh)}), also ist der Vertex mit "
      f"Phase phi unitär äquivalent zu dem mit phi - 4 theta, über eine Symmetrie des Zustands, von J, P_mu, K_mu und E "
      f"(1.1). Alle U(1)-invarianten Größen (Gram-Matrix der Ausgänge, |Amplitude|^2) sind phasenunabhängig",
      ok_uni and ok_inv and Evals == [5, 6])

# --- 1.4 Reflexionspositivität der Austauschamplitude zweiter Ordnung (radiale OS-Form, E-Grad 5 und 6) und g
E_k = np.array([float(E_von(k)) for k in KEYS])


def M_t(t):
    return A_amp.conj().T @ np.diag(np.exp(-t*E_k)) @ A_amp


G0 = np.eye(len(PAARE))                                            # <ab|cd> = delta für a < b, c < d
rp_ok = True
eig_tab = {}
for g in (0.1, 1.0, 10.0, 1000.0):
    for t in (0.0, 0.5, 2.0):
        ew = np.linalg.eigvalsh(G0 + g**2*M_t(t))
        rp_ok &= ew.min() > 1 - 1e-9
        eig_tab[f"g={g}, t={t}"] = float(ew.min())
lam_max = float(np.linalg.eigvalsh(M_t(0.0)).max())
g_krit = 1/math.sqrt(lam_max)
ghost = [float(np.linalg.eigvalsh(G0 - (fac*g_krit)**2*M_t(0.0)).min()) for fac in (0.99, 1.01)]
Gchi = G0 + np.exp(0.5j*np.pi)*M_t(0.0)
nonherm = float(np.abs(Gchi - Gchi.conj().T).max())
rangM = int(np.linalg.matrix_rank(M_t(0.0), tol=1e-10))
ZAHLEN["reflexionspositivitaet"] = {"min. Eigenwert G0 + g^2 M(t)": eig_tab, "Rang M": rangM, "lambda_max M(0)": lam_max,
                                    "Geistvorzeichen: g_krit": g_krit, "min. EW bei 0,99 / 1,01 g_krit": ghost,
                                    "relative Phase pi/2: |G - G^+|": nonherm}
check(f"Reflexionspositivität zweiter Ordnung: Die radiale OS-Form der Austauschamplitude ist G(g, t) = G_0 + "
      f"|g|^2 V^+ e^(-t E) V (Zwischenzustände E = 5, 6; Rang {rangM}), positiv für jedes g (0,1 bis 1000) und t (0 bis 2) "
      f"und für jede Phase von I, weil sich die Phase in V^+ V heraushebt. Sie scheitert nur, wenn der gespiegelte Vertex "
      f"nicht der adjungierte ist: relative Phase pi/2 gibt eine nicht hermitesche Form (|G - G^+| = {nonherm:.3f}), das "
      f"Geistvorzeichen G_0 - g^2 M wird oberhalb g_krit = {g_krit:.4f} indefinit (min. EW {ghost[0]:.4f} bzw. "
      f"{ghost[1]:.4f}). Positivität fixiert also die Realität (Vertex und Spiegelbild konjugiert), nicht |g|",
      rp_ok and ghost[0] > 0 and ghost[1] < 0 and nonherm > 1e-3 and rangM > 0)
befund("E: Der P_0-graduierte kubische Vertex ist für jede Modulphase und jedes reelle g reflexionspositiv (auf "
       "Abschneidungen; Stabilität der vollen Theorie ist eine andere Frage). Positivität fixiert die relative Phase "
       "zwischen Vertex und Spiegelvertex, die J ohnehin liefert, aber weder die Modulphase (Eichrichtung der freien "
       "Theorie) noch |g|. Der KMS-Zustand normiert |I|, nicht den Koeffizienten der Klammer: g bleibt frei", True)
laufzeit()

# =================================================================== Teil 2
say("\nTeil 2: MacDowell-Mansouri nativ")
# --- 2.1 Q^sigma = (P + sigma K)/2 aus der nativen radialen su(2,2)
ok_nativ = all(klein(komm(Pop[mu], ZT[a_]) - sum(NP[mu][a_, b_]*ZT[b_] for b_ in range(4))) and
               klein(komm(Kop[mu], ZT[a_]) - sum(NK[mu][a_, b_]*ZT[b_] for b_ in range(4))) for mu in range(4) for a_ in range(4))
ok_nativ &= all(klein(komm(Eop, ZT[a_]) - NE[a_, a_]*ZT[a_]) for a_ in range(4))
kombis = list(itertools.combinations(range(4), 2))
ETA = np.diag([1., -1., -1., -1.])
MKL = [ETA[m, m]*MK[m] for m in range(4)]                         # K_a = eta_ab K^b (die native K^mu trägt den oberen Index)
Mmn = {(m, n): 0.25*((MP[m] @ MKL[n] - MKL[n] @ MP[m]) - (MP[n] @ MKL[m] - MKL[m] @ MP[n])) for m, n in kombis}
Q_oben = [(MP[m] + MK[m])/2 for m in range(4)]                     # Kontrolle: ohne Indexsenkung
geschl_oben = all(im_spann(A @ B - B @ A, LOR + Q_oben) for A, B in itertools.combinations(LOR + Q_oben, 2))
alg = {}
ok_alg = ok_nativ and not geschl_oben
for s in (1, -1):
    Q = [(MP[m] + s*MKL[m])/2 for m in range(4)]
    basis = LOR + Q
    geschl = all(im_spann(A @ B - B @ A, basis) for A, B in itertools.combinations(basis, 2))
    der = derivationen(Inp - s*Ion)
    gleich = len(der) == 10 and all(im_spann(X, der) for X in basis)
    qq = all(np.allclose(Q[m] @ Q[n] - Q[n] @ Q[m], s*Mmn[(m, n)]) for m, n in kombis)
    m_lor = all(im_spann(Mmn[k], LOR) and np.abs(Mmn[k]).max() > 1e-3 for k in kombis)
    pp = all(np.allclose(MP[m] @ MP[n] - MP[n] @ MP[m], 0) and np.allclose(MKL[m] @ MKL[n] - MKL[n] @ MKL[m], 0)
             for m, n in kombis)
    cs, ts = [], []
    for (m, n), Mx in Mmn.items():
        for r in range(4):
            cs.append((Mx @ Q[r] - Q[r] @ Mx).ravel())
            ts.append((ETA[n, r]*Q[m] - ETA[m, r]*Q[n]).ravel())
    cv, tv = np.concatenate(cs), np.concatenate(ts)
    f_l = complex(np.vdot(tv, cv)/np.vdot(tv, tv))
    res_l = float(np.linalg.norm(cv - f_l*tv))
    sig = signatur(basis)
    alg[s] = {"geschlossen": geschl, "= Der(I_pi - sigma I_omega)": gleich, "[Q, Q] = sigma M": qq, "M in Lorentz": m_lor,
              "[P,P] = [K,K] = 0": pp, "Lorentzwirkung f": round(f_l.real, 12), "Residuum": res_l, "Signatur": sig,
              "E enthalten": im_spann(MEt, basis)}
    ok_alg &= geschl and gleich and qq and m_lor and pp and res_l < 1e-10 and abs(f_l.imag) < 1e-12
ok_alg &= alg[1]["Signatur"] == (6, 4, 0) and alg[-1]["Signatur"] == (4, 6, 0) and alg[1]["E enthalten"] \
    and not alg[-1]["E enthalten"] and np.allclose((MP[0] + MKL[0])/2, MEt)
f_lor = alg[1]["Lorentzwirkung f"]
ZAHLEN["mm_algebra"] = {str(k): {kk: (str(vv) if isinstance(vv, tuple) else vv) for kk, vv in v.items()} for k, v in alg.items()}
ZAHLEN["mm_algebra"]["Kontrolle ohne Indexsenkung geschlossen"] = bool(geschl_oben)
check(f"MacDowell-Mansouri-Algebra nativ: Mit den nativen P_mu = (1/2) c^+ sigma_mu c und K^mu = (1/2) omega''^+ "
      f"sigma_mu omega'' (Fock 8^4: [P, Z] = N_P Z, [K, Z] = N_K Z, [E, Z] = N_E Z), K_a = eta_ab K^b und M_mn = ([P_m, "
      f"K_n] - [P_n, K_m])/4 gilt [P, P] = [K, K] = 0, M_mn liegt in der Lorentzalgebra, [Q_m, Q_n] = sigma M_mn exakt und "
      f"[M_mn, Q_r] = f (eta_nr Q_m - eta_mr Q_n) mit f = {f_lor:.6f}, eta = diag(1, -1, -1, -1) (Residuum "
      f"{alg[1]['Residuum']:.1e}). span{{Q^sigma, Lorentz}} ist geschlossen und gleich Der(I_pi - sigma I_omega): "
      f"sigma = +1 gibt so(3,2) mit Spurform {alg[1]['Signatur'][:2]} und Q_0^+ = E exakt (z89s E-invarianter "
      f"AdS-Bivektor), sigma = -1 gibt so(4,1) {alg[-1]['Signatur'][:2]} ohne E. Kontrolle: Ohne Indexsenkung (P_mu + "
      f"K^mu)/2 schließt die Spanne nicht (Inversion vertauscht die Indexlage). Die Algebra der zwölften Analyse ist "
      f"genau z89s Büschel", ok_alg)

# --- 2.2 natives sigma und Radius: invarianter konformer Faktor und Krümmung (sympy, exakt)
XS = sp.symbols("x0:4", real=True)
SIGS = [sp.eye(2), sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
XM = sum((XS[k]*SIGS[k] for k in range(4)), sp.zeros(2))
TMS = sp.Matrix([[1, 0, -1, 0], [0, 1, 0, -1], [1, 0, 1, 0], [0, 1, 0, 1]])    # (w, c) = (z_j - z_j+2, z_j + z_j+2)
TMSi = TMS.inv()


def zu_sympy(M):
    return sp.Matrix(4, 4, [sp.nsimplify(complex(x).real) + sp.I*sp.nsimplify(complex(x).imag) for x in np.ravel(M)])


def vektorfeld(Mz):
    """Inzidenz w = i x c: der Twistorfluss Z -> Z + eps M Z bewegt die Gerade L_x um dx = A x - x D - i B - i x C x."""
    Mw = TMS*zu_sympy(Mz)*TMSi
    A, B, C, D = Mw[:2, :2], Mw[:2, 2:], Mw[2:, :2], Mw[2:, 2:]
    dx = A*XM - XM*D - sp.I*B - sp.I*XM*C*XM
    return [sp.expand((dx*SIGS[k]).trace()/2) for k in range(4)]


ETAs = sp.diag(1, -1, -1, -1)
GEN = {"P": [vektorfeld(M) for M in MP], "K": [vektorfeld(M) for M in MKL], "M": [vektorfeld(M) for M in Mmn.values()],
       "E": [vektorfeld(MEt)]}
ok_reell = all(sp.expand(c - sp.conjugate(c)) == 0 for L in GEN.values() for xi in L for c in xi)


def ist_ckv(xi):
    div = sum(sp.diff(xi[m], XS[m]) for m in range(4))
    xl = [sum(ETAs[m, n]*xi[n] for n in range(4)) for m in range(4)]
    return all(sp.expand(sp.diff(xl[n], XS[m]) + sp.diff(xl[m], XS[n]) - div*ETAs[m, n]/2) == 0
               for m in range(4) for n in range(4))


ok_ckv = all(ist_ckv(xi) for L in GEN.values() for xi in L)
DETX = sp.expand(XM.det())
s_sym = sp.Symbol("s")
OMEGA = 1 + s_sym*DETX


def killing_rest(xi, Om):
    div = sum(sp.diff(xi[m], XS[m]) for m in range(4))
    return sp.expand(sum(xi[m]*sp.diff(Om, XS[m]) for m in range(4)) - div*Om/4)


s_loes = {}
for sg in (1, -1):
    eqs = []
    for m in range(4):
        xi = [sp.expand(a + sg*b) for a, b in zip(GEN["P"][m], GEN["K"][m])]
        eqs += sp.Poly(killing_rest(xi, OMEGA), *XS).coeffs()
    s_loes[sg] = sp.solve(eqs, s_sym, dict=True)
ok_s = s_loes[1] == [{s_sym: 1}] and s_loes[-1] == [{s_sym: -1}]
ok_lor_om = all(killing_rest(xi, 1 + DETX) == 0 for xi in GEN["M"])


def ricci_skalar(gd, X, punkt):
    """Ricci-Skalar einer Diagonalmetrik (R^r_snm = d_n Gam^r_ms - ..., R_sm = R^r_srm), am Punkt ausgewertet."""
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


Lt = sp.Rational(3)
kugel = ricci_skalar([1/(1 + sum(x**2 for x in XS)/(4*Lt**2))**2]*4, XS, {XS[0]: sp.Rational(1, 3), XS[1]: sp.Rational(-1, 5),
                                                                       XS[2]: sp.Rational(2, 7), XS[3]: sp.Rational(1, 2)})
pkt = {XS[0]: sp.Rational(1, 5), XS[1]: sp.Rational(1, 3), XS[2]: sp.Rational(-1, 4), XS[3]: sp.Rational(1, 7)}
R_ads = ricci_skalar([-1/(1 + DETX)**2] + [1/(1 + DETX)**2]*3, XS, pkt)        # mostly plus: g = -eta/Omega^2
R_ds = ricci_skalar([-1/(1 - DETX)**2] + [1/(1 - DETX)**2]*3, XS, pkt)
Lam_ads, Lam_ds = R_ads/4, R_ds/4
l_ads, l_ds = sp.sqrt(3/abs(Lam_ads)), sp.sqrt(3/abs(Lam_ds))
xiE = GEN["E"][0]
vak = [{XS[0]: sp.I, XS[1]: 0, XS[2]: 0, XS[3]: 0}, {XS[0]: -sp.I, XS[1]: 0, XS[2]: 0, XS[3]: 0}]
ok_vakpkt = all(all(sp.simplify(c.subs(v)) == 0 for c in xiE) for v in vak) and \
    any(sp.simplify(c.subs({XS[0]: 2*sp.I, XS[1]: 0, XS[2]: 0, XS[3]: 0})) != 0 for c in xiE)
om_vak = sp.simplify((1 + DETX).subs(vak[0]))
ZAHLEN["radius"] = {"Omega (sigma = +1)": "1 + det x = 1 + x0^2 - |x|^2", "Omega (sigma = -1)": "1 - det x",
                    "Test Kugel L = 3: R": str(kugel), "R (AdS, mostly plus)": str(R_ads), "R (dS)": str(R_ds),
                    "Lambda AdS": str(Lam_ads), "Lambda dS": str(Lam_ds), "l (Zustandseinheiten, Im x_Vak = 1)": str(l_ads),
                    "l in CCR-Einheiten (m = 1, 2, 3)": {m: KAP2_NATIV[m]/2 for m in (1, 2, 3)},
                    "Fixpunkte von E": "x = +-i (Im x^0 = +-1)", "Omega(x_Vak)": str(om_vak),
                    "z89-Konvention": "Lambda-hat = Lambda/12, |Lambda-hat| (Im x)^2 = 1"}
check(f"Natives sigma und Radius (sympy, exakt): Die Twistorflüsse wirken auf x (Inzidenz omega'' = i x c, reelle Punkte "
      f"x = x^mu sigma_mu hermitesch) als reelle konforme Killingfelder von eta = diag(1, -1, -1, -1); P + sigma K "
      f"und die Lorentzgeneratoren lassen genau g = eta/Omega^2 mit Omega = 1 + sigma det x invariant (s = "
      f"{[v[0][s_sym] if v else None for v in s_loes.values()]} für sigma = +1, -1). Krümmung (Ricci-Skalar direkt, Test Kugel L = 3: R = {kugel} = 12/9): "
      f"sigma = +1 gibt R = {R_ads}, Lambda = {Lam_ads} (AdS), sigma = -1 gibt R = {R_ds} (dS), beide mit l = {l_ads}. "
      f"Die Fixpunkte der kompakten Energie E = (P_0 + K_0)/2 sind x = +-i (Im x^0 = 1, dort Omega = {om_vak}): Der "
      f"AdS-Radius ist die Hälfte der Vakuumgröße, l = Im x_Vak/2, also Lambda = -12/(Im x_Vak)^2 und z89s Lambda-hat = "
      f"Lambda/12; in CCR-Einheiten l = tanh(pi m)/2 ({KAP2_NATIV[1]/2:.6f} für m = 1)",
      ok_reell and ok_ckv and ok_s and ok_lor_om and kugel == sp.Rational(4, 3) and R_ads == -48 and R_ds == 48
      and l_ads == sp.Rational(1, 2) and ok_vakpkt and om_vak == 0)
befund("E/C: Nativ ist sigma = +1 (AdS): Die radiale Energie E (Moden = Stufen, z83) ist kompakt und liegt nur in "
       "span{(P + K)/2, Lorentz} = so(3,2). Der Radius ist l = Im x_Vak/2 in Einheiten der Vakuumgröße; eine absolute "
       "Länge gibt es nativ nicht (Dilatationsskala P2 aus z84). Die dS-Wahl sigma = -1 braucht die hyperbolische "
       "Q_0 = (P_0 - K_0)/2 (Teil 3)", True)
laufzeit()

# --- 2.3 tau als Punktspiegelung von AdS4, modulare Berry-Verbindung
Qp = [(MP[m] + MKL[m])/2 for m in range(4)]
APLUS = LOR + Qp
H_tau = teilalgebra([lambda M: M @ S4 - S4 @ M], basis=APLUS)
ok_h = len(H_tau) == 6 and all(im_spann(X, LOR) for X in H_tau)
ok_ad = all(np.allclose(S4 @ X @ S4, X) for X in LOR) and all(np.allclose(S4 @ X @ S4, -X) for X in Qp)
ok_mm = all(im_spann(Qp[a] @ Qp[b] - Qp[b] @ Qp[a], LOR) for a, b in kombis)
Bm = np.array([[np.trace(Qp[a] @ Qp[b]).real for b in range(4)] for a in range(4)])
cs, ts = [], []
for a, b, c in itertools.product(range(4), repeat=3):
    X = Qp[a] @ Qp[b] - Qp[b] @ Qp[a]
    cs.append((X @ Qp[c] - Qp[c] @ X).ravel())
    ts.append((Bm[b, c]*Qp[a] - Bm[a, c]*Qp[b]).ravel())
cv, tv = np.concatenate(cs), np.concatenate(ts)
k_alg = complex(np.vdot(tv, cv)/np.vdot(tv, tv))
res_k = float(np.linalg.norm(cv - k_alg*tv))
# physikalische Normierung: Q_mu bewegt den Ursprung mit Geschwindigkeit 1/2 (P_mu ist die Einheitstranslation, Teil 2.2),
# also g_phys = B/4 auf m; Schnittkrümmung K_B = -k_alg, K_phys = 4 K_B = -1/l^2
xi_Q = [vektorfeld(M) for M in Qp]
v_urspr = [[sp.nsimplify(c.subs({x: 0 for x in XS})) for c in xi] for xi in xi_Q]
ok_geschw = all(abs(v_urspr[m][n]) == (sp.Rational(1, 2) if m == n else 0) for m in range(4) for n in range(4))
K_phys = -4*k_alg.real
l_alg = 1/math.sqrt(-K_phys) if K_phys < 0 else float("nan")
Ph = lambda X: (X + S4 @ X @ S4)/2
Pm = lambda X: (X - S4 @ X @ S4)/2
XB, YB = Qp[1] + 0.3*Qp[0], Qp[2] - 0.5*Qp[3]


def A_s(t):
    return sla.expm(-t*YB) @ XB @ sla.expm(t*YB)                 # g^-1 d_s g für g = e^(s X) e^(t Y); A_t = Y


t0, hh = 0.37, 1e-5
dAs = (A_s(t0 + hh) - A_s(t0 - hh))/(2*hh)
As0 = A_s(t0)
F_A = -dAs + (As0 @ YB - YB @ As0)                                # d_s A_t - d_t A_s + [A_s, A_t]
Gam_s, e_s = Ph(As0), Pm(As0)
F_Gam = -Ph(dAs)                                                  # Gamma_t = P_h(Y) = 0
torsion = -Pm(dAs) + (Gam_s @ YB - YB @ Gam_s)
riem = -(e_s @ YB - YB @ e_s)
ok_berry = (np.abs(F_A).max() < 1e-8 and np.allclose(F_Gam, riem, atol=1e-8) and np.abs(torsion).max() < 1e-8
            and np.abs(F_Gam).max() > 1e-2 and np.allclose(Ph(YB), 0))
ok_bahn = len(APLUS) - len(H_tau) == 4
# Negativkontrolle: relative Modulkette (Connes-Kozyklus) dreier Zustände ist flach


def potenz(rho, z):
    w, V = np.linalg.eigh(rho)
    return V @ np.diag(w.astype(complex)**z) @ V.conj().T


nm = 10
r1 = np.diag(0.3**np.arange(nm))
r1 /= np.trace(r1)
Wr = np.linalg.qr(rng.normal(size=(nm, nm)) + 1j*rng.normal(size=(nm, nm)))[0]
r2 = Wr @ np.diag(0.5**np.arange(nm)) @ Wr.conj().T
r2 /= np.trace(r2).real
Xr = rng.normal(size=(nm, nm)) + 1j*rng.normal(size=(nm, nm))
r3 = Xr @ Xr.conj().T + 0.1*np.eye(nm)
r3 /= np.trace(r3).real
uk = lambda ra, rb, t: potenz(ra, 1j*t) @ potenz(rb, -1j*t)
hol = max(np.abs(uk(r1, r2, t) @ uk(r2, r3, t) @ uk(r3, r1, t) - np.eye(nm)).max() for t in (0.3, 1.7, -2.2))
s_, t_ = 0.6, -1.1
koz = np.abs(uk(r1, r2, s_ + t_) - uk(r1, r2, s_) @ potenz(r2, 1j*s_) @ uk(r1, r2, t_) @ potenz(r2, -1j*s_)).max()
ok_koz = hol < 1e-10 and koz < 1e-10 and np.abs(uk(r1, r2, 0.9) - np.eye(nm)).max() > 1e-2
ZAHLEN["berry"] = {"Stabilisator von tau in so(3,2)": len(H_tau), "Bahn von tau (Dimension)": len(APLUS) - len(H_tau),
                   "Spurform B(Q_a, Q_b)": np.round(Bm, 12).tolist(), "konstante Krümmung k_alg": round(k_alg.real, 12),
                   "Residuum konstante Krümmung": res_k, "Geschwindigkeit von Q_mu im Ursprung": "1/2",
                   "Schnittkrümmung physikalisch": K_phys, "l algebraisch": l_alg,
                   "|F_A| (MM-Krümmung, Maurer-Cartan)": float(np.abs(F_A).max()),
                   "|F_Gamma| (Berry-Krümmung)": float(np.abs(F_Gam).max()), "|Torsion|": float(np.abs(torsion).max()),
                   "Kozyklus-Holonomie (Dreieck)": float(hol), "Kozyklus-Identität": float(koz)}
check(f"tau als Punktspiegelung von AdS4: Der Stabilisator von tau = Theta_1 in so(3,2) ist genau die Lorentzalgebra "
      f"(Dimension {len(H_tau)}), Ad tau = +1 auf Lorentz und -1 auf Q^+_mu, [Q, Q] liegt in Lorentz: (so(3,2), so(3,1)) "
      f"ist das symmetrische Paar von AdS4, die tau-Bahn ist vierdimensional, und [[Q_a, Q_b], Q_c] = k (B_bc Q_a - B_ac "
      f"Q_b) mit k = {k_alg.real:.6f} und B = Spurform diag{np.round(np.diag(Bm), 6).tolist()} (Residuum {res_k:.1e}): "
      f"konstante Krümmung; Q_mu bewegt den Ursprung mit Geschwindigkeit 1/2, also K_phys = -4 k = {K_phys:.6f} und "
      f"l = {l_alg:.6f}, unabhängig von Teil 2.2 bestätigt. Berry-Verbindung der Tomita-Familie tau_g = g tau g^-1 (Nullmoden = Kommutant von tau = "
      f"Lorentz) auf g = e^(sX) e^(tY): MM-Krümmung F_A = {np.abs(F_A).max():.1e} (Maurer-Cartan), Berry-Krümmung "
      f"F_Gamma = -[e, e] != 0 (|F| = {np.abs(F_Gam).max():.4f}), Torsion {np.abs(torsion).max():.1e}: die "
      f"Levi-Civita-Verbindung von AdS4 mit Vierbein e = P_m(g^-1 dg). Negativkontrolle: Der Connes-Kozyklus "
      f"(Dphi : Dpsi)_t dreier Zustände ist flach (Dreiecksholonomie {hol:.1e}, Kozyklus {koz:.1e}); der native "
      f"Modulgenerator D ist auf T zentral (1.1), seine Bahn ein Punkt, seine Berry-Verbindung null",
      ok_h and ok_ad and ok_mm and res_k < 1e-10 and abs(k_alg.imag) < 1e-12 and ok_berry and ok_bahn and ok_koz
      and ok_geschw and abs(l_alg - 0.5) < 1e-12)
befund("E/C: Die modulare Berry-Verbindung liefert nativ genau das MacDowell-Mansouri-Vakuum. Der Modulfluss selbst ist "
       "auf T zentral (keine Verbindung), die relative Modulkette ist flach; nicht trivial ist allein die Familie der "
       "Tomita-Konjugationen J_g = U(g) J U(g)^-1 (Original|Partner an jedem Punkt), auf dem Twistorraum tau_g. Unter der "
       "AdS-Gruppe ist sie AdS4 selbst, ihre Nullmoden-Verbindung ist die Levi-Civita-Verbindung (Riemann = -[e, e]) und "
       "die volle Cartan-Verbindung ist flach: F^ab = R^ab + (sigma/l^2) e^a e^b = 0. Das ist Kinematik eines starren "
       "Vakuums, keine Dynamik", True)

# --- 2.4 Plebanski-Zweiformen aus zwei Blättern und Phasentest
EPS_UNTEN = np.einsum("abcd,ai,bj,ck,dl->ijkl", LC, ETA, ETA, ETA, ETA)


def stern(F):
    return 0.5*np.einsum("mnrs,rs->mn", EPS_UNTEN, ETA @ F @ ETA)


def keil(F, G_):
    return 0.25*np.einsum("mnrs,mn,rs->", LC, F, G_)


Fr = rng.normal(size=(4, 4))
Fr = Fr - Fr.T
ok_stern2 = np.allclose(stern(stern(Fr)), -Fr)
SIGP = [[np.array([[(SIGMA[m] @ EPS2 @ SIGMA[n].T - SIGMA[n] @ EPS2 @ SIGMA[m].T)[i, k] for n in range(4)]
                   for m in range(4)]) for k in range(2)] for i in range(2)]                 # eps auf dem Originalindex (I_pi)
SIGM = [[np.array([[(SIGMA[m].T @ EPS2 @ SIGMA[n] - SIGMA[n].T @ EPS2 @ SIGMA[m])[j, l] for n in range(4)]
                   for m in range(4)]) for l in range(2)] for j in range(2)]                 # eps auf dem Partnerindex
ok_sym = all(np.allclose(SIGP[i][k], SIGP[k][i]) for i in range(2) for k in range(2))
ok_konj = all(np.allclose(SIGM[i][k], SIGP[i][k].conj()) for i in range(2) for k in range(2))
lam_d = complex(np.vdot(SIGP[0][0].ravel(), stern(SIGP[0][0]).ravel())/np.vdot(SIGP[0][0].ravel(), SIGP[0][0].ravel()))
ok_dual = all(np.allclose(stern(SIGP[i][k]), lam_d*SIGP[i][k]) and np.allclose(stern(SIGM[i][k]), np.conj(lam_d)*SIGM[i][k])
              for i in range(2) for k in range(2)) and abs(abs(lam_d) - 1) < 1e-12 and abs(lam_d.real) < 1e-12
ok_kreuz = all(abs(keil(SIGP[i][k], SIGM[j][l])) < 1e-12 for i, k, j, l in itertools.product(range(2), repeat=4))
Wt = np.array([[keil(SIGP[i][k], SIGP[j][l]) for j in range(2) for l in range(2)] for i in range(2) for k in range(2)])
Tt = np.array([[EPS2[i, j]*EPS2[k, l] + EPS2[i, l]*EPS2[k, j] for j in range(2) for l in range(2)]
               for i in range(2) for k in range(2)])
c_vol = complex(np.vdot(Tt.ravel(), Wt.ravel())/np.vdot(Tt.ravel(), Tt.ravel()))
ok_simpl = np.allclose(Wt, c_vol*Tt) and abs(c_vol.real) < 1e-12 and abs(c_vol) > 1e-3
ok_twoJ = all(klein(J_ad(CC[j], SW) - dag(CC[j])) for j in range(2)) and all(
    klein(Pop[mu] - 0.5*sum(SIGMA[mu][i, j]*(J_ad(CC[i], SW) @ CC[j]) for i in range(2) for j in range(2)))
    for mu in range(4))
# reelle Lorentzform B_phi = e^(i phi) Sigma_+ + e^(-i phi) Sigma_- und Plebanski-Simplizität
EV = np.zeros((4, 4, 4, 4))
for I_, J_ in itertools.product(range(4), repeat=2):
    EV[I_, J_, I_, J_] += 1
    EV[I_, J_, J_, I_] -= 1
Pplus = lambda F: 0.5*(F + stern(F)/lam_d)
SPV = np.array([[Pplus(EV[I_, J_]) for J_ in range(4)] for I_ in range(4)])
rang_sd = np.linalg.matrix_rank(np.array([SPV[I_, J_].ravel() for I_, J_ in kombis]
                                         + [SIGP[i][k].ravel() for i, k in ((0, 0), (0, 1), (1, 1))]), tol=1e-10)


def simpl_rest(phi):
    B = np.exp(1j*phi)*SPV + np.exp(-1j*phi)*SPV.conj()
    K = np.array([[[[keil(B[a, b], B[c, d]) for d in range(4)] for c in range(4)] for b in range(4)] for a in range(4)])
    cK = np.sum(LC*K)/np.sum(LC*LC)
    return float(np.linalg.norm(K - cK*LC)), B, cK


PHG = [j*np.pi/24 for j in range(48)]
rest = [simpl_rest(p)[0] for p in PHG]
nullen = [j for j, r in enumerate(rest) if r < 1e-10]
quot = [r/abs(math.sin(2*p)) for r, p in zip(rest, PHG) if abs(math.sin(2*p)) > 1e-6]
r0, B0, c0 = simpl_rest(0.0)
r90, B90, c90 = simpl_rest(np.pi/2)
stE = np.array([[stern(EV[I_, J_]) for J_ in range(4)] for I_ in range(4)])
ok_sek = np.allclose(B0, EV) and np.allclose(B90, stE*(1j/lam_d).real) and np.allclose(simpl_rest(np.pi)[1], -EV)
vol_phase = [abs((np.exp(2j*p)*c_vol).real) < 1e-12 for p in PHG]
ok_volph = [j for j, v in enumerate(vol_phase) if v] == [0, 12, 24, 36]
ZAHLEN["plebanski"] = {"Dualität Sigma_+ (eps auf Originalindex)": str(np.round(lam_d, 12)),
                       "Sigma_+ ^ Sigma_+ = c (eps eps + eps eps), c": str(np.round(c_vol, 12)),
                       "Rang selbstdual (vektoriell + spinoriell)": int(rang_sd), "Nullstellen der Simplizität (j pi/24)": nullen,
                       "Rest/|sin 2 phi| (min, max)": [min(quot), max(quot)],
                       "Sektoren": f"phi = 0: B = e^e; pi/2: {(1j/lam_d).real:+.0f} *(e^e); pi: -e^e; 3 pi/2: "
                                   f"{-(1j/lam_d).real:+.0f} *(e^e)",
                       "Volumenrealität Re(e^(2 i phi) c) = 0 bei j": [j for j, v in enumerate(vol_phase) if v]}
check(f"Plebanski aus zwei Blättern: J c J = c^+ und P_mu = (1/2) (J c J) sigma_mu c (Fock 8^4): das Vierbein "
      f"sigma_mu^(AA') paart den Partnerindex (A, J-Bild) mit dem Originalindex (A'). Sigma_+ = sigma_m eps sigma_n^T - "
      f"(m <-> n) (eps auf dem Originalindex, also I_pi) ist symmetrisch und selbstdual (*Sigma_+ = {lam_d:.0f} "
      f"Sigma_+), Sigma_- = conj(Sigma_+) ist genau die eps-Kontraktion des Partnerindex und antiselbstdual, Sigma_+ ^ "
      f"Sigma_- = 0, Sigma_+ ^ Sigma_+ = c (eps eps + eps eps) mit c = {c_vol:.0f} (rein imaginär; selbstduale Formen "
      f"Rang {rang_sd}). Phasentest: Die reelle Lorentzform B_phi = e^(i phi) Sigma_+ + e^(-i phi) Sigma_- = cos phi "
      f"e^e {'+' if (1j/lam_d).real > 0 else '-'} sin phi *(e^e) erfüllt die Plebanski-Simplizität B^IJ ^ B^KL ~ eps^IJKL genau bei phi = j pi/24 mit j in "
      f"{nullen} (Rest/|sin 2 phi| = {min(quot):.6f} konstant), ebenso die Volumenrealität Re(Sigma ^ Sigma) = 0: vier "
      f"Sektoren +-e^e, +-*(e^e)",
      ok_stern2 and ok_sym and ok_konj and ok_dual and ok_kreuz and ok_simpl and ok_twoJ and rang_sd == 3
      and nullen == [0, 12, 24, 36] and max(quot) - min(quot) < 1e-9 and ok_sek and ok_volph and r0 < 1e-12)
befund("C: Die zwei Blätter liefern nativ Sigma_+ (Original, eps = I_pi) und Sigma_- = J Sigma_+ J = conj(Sigma_+) "
       "(Partner) samt dem flachen Vierbein sigma_mu aus P_mu, für jede Modulphase. Die Phase von I fixiert erst die "
       "Plebanski-Simplizität der reellen Form, und zwar nur bis Z_4: phi in (pi/2) Z, vier Sektoren. Die native "
       "Modulrotation U_(pi/4) = e^(i pi D/4) bildet den Gravitationssektor e^e auf den topologischen *(e^e) ab und "
       "vertauscht mit allen nativen Daten; die Wahl des Sektors (und der Orientierung) ist eine Prämisse, ebenso die "
       "Simplizität selbst (Wirkung) und ein dynamisches Vierbein", True)
laufzeit()

# =================================================================== Teil 3
say("\nTeil 3: Skalen")
# --- 3.1 MacDowell-Mansouri-Zerlegung auf Formen
EE2 = lambda a, b: np.outer(np.eye(4)[a], np.eye(4)[b]) - np.outer(np.eye(4)[b], np.eye(4)[a])
s_eeee = sum(LC[a, b, c, d]*keil(EE2(a, b), EE2(c, d)) for a, b, c, d in itertools.permutations(range(4)))
Rr = rng.normal(size=(4, 4, 4, 4))
Rr = Rr - Rr.transpose(1, 0, 2, 3)
Rr = Rr - Rr.transpose(0, 1, 3, 2)
Rform = [[0.5*sum(Rr[a, b, e, f]*EE2(e, f) for e in range(4) for f in range(4)) for b in range(4)] for a in range(4)]
s_Ree = sum(LC[a, b, c, d]*keil(Rform[a][b], EE2(c, d)) for a, b, c, d in itertools.permutations(range(4)))
R_skalar = sum(Rr[a, b, a, b] for a in range(4) for b in range(4))
bS, sgS, lS, GS, LamS, kS = sp.symbols("b sigma l G Lambda k")
loes = sp.solve([sp.Eq(4*bS*sgS/lS**2, 1/(16*sp.pi*GS)), sp.Eq(24*bS*sgS**2/lS**4, -2*LamS/(16*sp.pi*GS))], [bS, LamS],
                dict=True)[0]
b_mm = sp.simplify(loes[bS])
Lam_mm = sp.simplify(loes[LamS].subs(sgS**2, 1).subs(sgS**3, sgS))
k_vak = sp.solve(sp.Eq(kS + sgS/lS**2, 0), kS)[0]
ok_mmz = (abs(s_eeee - 24) < 1e-12 and abs(s_Ree - 2*R_skalar) < 1e-9 and sp.simplify(b_mm - lS**2/(64*sp.pi*GS*sgS)) == 0
          and sp.simplify(Lam_mm + 3*sgS/lS**2) == 0 and sp.simplify(3*k_vak - Lam_mm) == 0)
ZAHLEN["mm_zerlegung"] = {"eps_abcd e^a e^b e^c e^d / vol": float(s_eeee), "eps R e e / (R vol)": float(s_Ree/R_skalar),
                          "b": str(b_mm), "Lambda": str(Lam_mm), "Vakuum F = 0: k": str(k_vak)}
check(f"MacDowell-Mansouri-Zerlegung (Formenalgebra, zufälliger Krümmungstensor): eps_abcd e^a e^b e^c e^d = "
      f"{s_eeee:.0f} vol, eps_abcd R^ab e^c e^d = {s_Ree/R_skalar:.6f} R vol; mit F = R + (sigma/l^2) e e ist "
      f"b eps F F = Euler + (4 b sigma/l^2) R vol + (24 b sigma^2/l^4) vol, und der Vergleich mit (1/16 pi G)(R - 2 "
      f"Lambda) vol gibt b = {b_mm}, Lambda = {Lam_mm} (wie die zwölfte Analyse); das Vakuum F = 0 hat R^ab = k e^a e^b "
      f"mit k = {k_vak}, also R_ab = 3 k g_ab = Lambda g_ab (konsistent)", ok_mmz)

# --- 3.2 bedingte Identität aus v60 (rho_Lambda) und MM
c3 = 1/(8*sp.pi)
pref = sp.simplify((8*sp.pi)**2*48*c3**4)
aS = sp.Symbol("alpha", positive=True)
rho = pref*sp.exp(-2/aS)                                           # rho_Lambda / Mbar^4 (v60, Ledger GRAV.NONLINEAR.01)
l2G = sp.simplify(24*sp.pi/rho)                                    # rho/Mbar^4 = 8 pi G |Lambda| = 24 pi G/l^2
b_betrag = sp.simplify(l2G/(64*sp.pi))
ok_v60 = (sp.simplify(pref - 3/(4*sp.pi**2)) == 0 and sp.simplify(l2G - 32*sp.pi**3*sp.exp(2/aS)) == 0
          and sp.simplify(b_betrag - sp.pi**2/2*sp.exp(2/aS)) == 0)
lg_e = 2*ALPHA_INV/mp.log(10)
zahl = {"log10 e^(2/alpha)": float(lg_e), "log10 l^2/G": float(mp.log10(32*mp.pi**3) + lg_e),
        "log10 |b|": float(mp.log10(mp.pi**2/2) + lg_e), "log10 l/sqrt(G)": float((mp.log10(32*mp.pi**3) + lg_e)/2),
        "log10 G/R_Vak^2 (bei l = R_Vak/2)": float(-lg_e - mp.log10(128*mp.pi**3))}
ZAHLEN["skalen"] = {"(8 pi)^2 48 c3^4": str(pref), "l^2/G": str(l2G), "|b|": str(b_betrag),
                    "alpha^-1 (CODATA, nur Größenordnung)": str(ALPHA_INV), **zahl}
check(f"Bedingte Identität (sympy): Aus v60 (rho_Lambda/Mbar^4 = (8 pi)^2 48 c3^4 e^(-2/alpha) = {pref} e^(-2/alpha), "
      f"c3 = 1/(8 pi), nur gelesen) und den MM-Relationen |Lambda| = 3/l^2, rho_Lambda = Lambda Mbar^2, b = l^2/(64 pi G) "
      f"folgen l^2/G = {l2G} und |b| = {b_betrag}, wie die zwölfte Analyse. Größenordnungen (alpha^-1 = 137,036): "
      f"log10(l^2/G) = {zahl['log10 l^2/G']:.3f}, log10|b| = {zahl['log10 |b|']:.3f}, l/sqrt(G) = 10^"
      f"{zahl['log10 l/sqrt(G)']:.3f}", ok_v60)

# --- 3.3 Vorzeichen und Größenordnung gegen den nativen Radius
Par = sps.diags(((-1.0)**(OCC[:, 0] + OCC[:, 1]))).tocsr()
Hhyp = (Pop[0] - Kop[0]).tocsr()                                   # 2 Q_0^- = P_0 - K_0 = Summe (d^+ e^+ + e d)
ok_hyp = klein(Par @ Hhyp @ Par + Hhyp) and not klein(Hhyp) and klein(Hhyp - sum(
    (dag(DD[i]) @ dag(EE[i]) + EE[i] @ DD[i] for i in range(2)), sps.csr_matrix((DIMV, DIMV))))
sign_nativ = int(sp.sign(Lam_ads))
sign_v60 = int(sp.sign(pref))
ok_vz = sign_nativ == -1 and sign_v60 == 1
ZAHLEN["vorzeichen"] = {"sign Lambda nativ (sigma = +1)": sign_nativ, "sign rho_Lambda v60": sign_v60,
                        "Q_0^- = (P_0 - K_0)/2 spektralsymmetrisch": True,
                        "nötiges G/R_Vak^2": "e^(-2/alpha)/(128 pi^3)"}
check(f"Vorzeichen und Skala: Das native sigma = +1 gibt Lambda = -3/l^2 < 0 und damit rho_Lambda = Lambda/(8 pi G) < 0, "
      f"v60 hat rho_Lambda > 0 (Vorzeichen {sign_nativ} gegen {sign_v60}). Die dS-Wahl sigma = -1 braucht Q_0^- = (P_0 - "
      f"K_0)/2 = (1/2) Summe (d^+ e^+ + e d): (-1)^(N_d) Q_0^- (-1)^(N_d) = -Q_0^- (Fock 8^4), das Spektrum ist "
      f"symmetrisch, nicht nach unten beschränkt, keine Moden-gleich-Stufen-Energie; sie ist der hyperbolische "
      f"(Doppelkegel-)Erzeuger. Mit l = R_Vak/2 verlangt die Relation G/R_Vak^2 = e^(-2/alpha)/(128 pi^3) = 10^"
      f"{zahl['log10 G/R_Vak^2 (bei l = R_Vak/2)']:.2f}; nativ ist G nicht bestimmt (z89: g frei, z90: induziertes "
      f"1/G = 0)", ok_vz and ok_hyp)
befund("E/O: Die Lambda-G-Relation der zwölften Analyse ist algebraisch richtig, aber nur eine bedingte Identität: Sie "
       "setzt die historische rho_Lambda-Formel (v60) und die MM-Wirkung voraus und bestimmt dann l^2/G bzw. b. Eine "
       "native Aussage ist sie nicht: Das native sigma = +1 gibt das falsche Vorzeichen (AdS statt dS), der native "
       "Radius ist die halbe Vakuumgröße ohne absolute Skala, und G (bzw. b) folgt aus keiner nativen Struktur. Sie ist "
       "ein Konsistenzziel: G/R_Vak^2 = e^(-2/alpha)/(128 pi^3) bei dS-Vorzeichen, das nativ nicht erreichbar ist. "
       "Lesart [L]: Das native so(3,2) ist Kinematik der Naht (Radius = halbe Vakuumgröße), nicht die Eichgruppe der "
       "physikalischen Lambda", True)
laufzeit()

# =================================================================== Teil 4
say("\nTeil 4: Bilanz und Urteil")
STATUS = {
    "Modulphase von I und Theta_1 = tau": ("E, negativ", "alle nativen Daten Modul-U(1)-invariant, I Ladung 2; "
                                          "Theta-Realität J V J = V^+ für jede Phase; Phase entfernbar (Eichrichtung)"),
    "Realstrukturen": ("E", "lorentzsch: Phase von I = Phase der Volumenform (SU(2,2)-Reduktion, Rest +-I); "
                            "quaternionisch: Modulbahn sigma_theta, jede fixiert eine andere Phase"),
    "Reflexionspositivität und g": ("E, negativ", "OS-Form zweiter Ordnung positiv für jedes g und jede Phase; "
                                    "fixiert nur Realität (Vertex = Spiegelvertex^+), nicht |g|"),
    "MM-Algebra": ("E", "[Q_m, Q_n] = sigma M_mn, span = Der(I_pi - sigma I_omega); sigma = +1: so(3,2) mit Q_0 = E"),
    "natives sigma, Radius": ("E/C", "sigma = +1 (E kompakt), l = Im x_Vak/2, Lambda = -12/(Im x_Vak)^2; "
                                     "z89s Lambda-hat = Lambda/12; keine absolute Länge"),
    "tau als AdS-Punktspiegelung": ("E", "Stabilisator Lorentz, Ad tau = -1 auf Q, Bahn AdS4, konstante Krümmung"),
    "modulare Berry-Verbindung": ("E/C", "Modulgenerator zentral (trivial), Connes-Kozyklus flach; Tomita-Familie gibt "
                                         "Levi-Civita von AdS4, MM-Krümmung null: Vakuum, keine Dynamik"),
    "Plebanski aus zwei Blättern": ("C", "Sigma_+ (I_pi), Sigma_- = J-Bild = conj, flaches Vierbein aus P_mu; "
                                         "Simplizität fixiert die Phase bis Z_4 (vier Sektoren), Sektorwahl Prämisse"),
    "Lambda-G-Relation": ("E bedingt, O nativ", "l^2/G = 32 pi^3 e^(2/alpha), |b| = (pi^2/2) e^(2/alpha) aus v60 + MM; "
                                                "natives Vorzeichen falsch (AdS), G nicht nativ: Konsistenzziel"),
    "Prämissen": ("O", "dynamisches Vierbein bzw. Cartan-Verbindung, Wirkung (MM/Plebanski) mit Koeffizient b, "
                       "Sektorwahl, dS-Vorzeichen, absolute Skala"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
urteil = ("G2 ist beantwortet, überwiegend als gescoptes Negativ mit einer präzisierten Kinematik. (1) Die native "
          "OS-Spiegelung Theta_1 = tau fixiert die Modulphase von I nicht: tau ist auf T linear, J vertauscht mit dem "
          "Modulfluss, und jede native Struktur (KMS-Zustand, J, Delta^(it), tau, Theta_n, P_mu, K_mu, E, CCR-Form) ist "
          "Modul-U(1)-invariant, während I_pi die Ladung 2 trägt; die Theta-Realität J V J^-1 = V^+ des eps-Kerns gilt auf "
          "zwei Nähten für jede Phase. Die Phase ist eine Eichrichtung der freien Theorie: Konjugation mit U_theta "
          "multipliziert jede Gravitonamplitude gleichförmig mit e^(4 i theta). Realstrukturen, die sie fixieren würden, "
          "kommen in Modulbahnen: lorentzsch ist die Phase von I gleich der Phase einer Volumenform (Reduktion U(2,2) -> "
          "SU(2,2), Rest +-I), quaternionisch fixiert jedes sigma_theta eine andere Phase. Der P_0-graduierte kubische "
          "Vertex ist reflexionspositiv für jede Phase und jedes reelle g (OS-Form zweiter Ordnung G_0 + |g|^2 V^+ e^(-tE) V); "
          "Positivität verlangt nur, dass der Spiegelvertex der adjungierte ist (sonst nicht hermitesch, mit Geistvorzeichen "
          "indefinit oberhalb g_krit), und fixiert g nicht. (2) Die MacDowell-Mansouri-Algebra der zwölften Analyse ist "
          "exakt nativ: Mit P_mu, K_mu aus der radialen su(2,2) gilt [Q_m, Q_n] = sigma M_mn, und span{Q^sigma, Lorentz} "
          "= Der(I_pi - sigma I_omega), also genau z89s Büschel; nativ ist sigma = +1 (Q_0 = E, Moden = Stufen, so(3,2)). "
          "Der invariante konforme Faktor ist Omega = 1 + det x, der Ricci-Skalar -48, also AdS mit l = Im x_Vak/2 "
          "(Präzisierung von z89: Radius gleich der halben Vakuumgröße, Lambda-hat = Lambda/12). tau = Theta_1 ist die "
          "Punktspiegelung dieses AdS4 (Stabilisator Lorentz, Ad tau = -1 auf Q, konstante Krümmung). Die modulare "
          "Berry-Verbindung gibt nativ nur das MM-Vakuum: Der Modulgenerator ist auf T zentral, die relative Modulkette "
          "(Connes-Kozyklus) flach, und die Familie der Tomita-Konjugationen J_g ist AdS4 mit der Levi-Civita-Verbindung als "
          "Nullmoden-Verbindung und verschwindender MM-Krümmung F = R + (sigma/l^2) e e = 0. Die zwei Blätter liefern "
          "Sigma_+ (eps auf dem Originalindex, I_pi), Sigma_- = J Sigma_+ J = conj(Sigma_+) und das flache Vierbein aus "
          "P_mu = (1/2)(J c J) sigma c; die Plebanski-Simplizität der reellen Form fixiert die Phase von I bis Z_4 (vier "
          "Sektoren +-e^e, +-*(e^e)), und die native Modulrotation U_(pi/4) verbindet sie. Prämissen bleiben: ein "
          "dynamisches Vierbein bzw. eine nicht flache Cartan-Verbindung, die Wirkung mit ihrem Koeffizienten b, die "
          "Sektorwahl. (3) Die Lambda-G-Relation l^2/G = 32 pi^3 e^(2/alpha), |b| = (pi^2/2) e^(2/alpha) ist als bedingte "
          "Identität aus v60 und den MM-Relationen richtig (b = l^2/(64 pi G sigma), Lambda = -3 sigma/l^2 auf Formen "
          "geprüft), aber keine native Aussage: Das native sigma = +1 gibt Lambda < 0 gegen rho_Lambda > 0, die dS-Wahl "
          "bräuchte den hyperbolischen, spektralsymmetrischen Erzeuger (P_0 - K_0)/2, und G bleibt frei (G/R_Vak^2 = "
          "e^(-2/alpha)/(128 pi^3) wäre nötig). (4) Nativ für die Gravitation sind damit: I_pi bis auf die Modulphase, "
          "die AdS-Kinematik so(3,2) mit l = R_Vak/2, tau als Punktspiegelung, die Plebanski-Formen Sigma_+-, das flache "
          "Vierbein und das MM-Vakuum als modulare Berry-Verbindung. Nicht nativ: Phase (bis auf eine Volumenform bzw. "
          "einen Plebanski-Sektor), g bzw. G bzw. b, das dS-Vorzeichen und jede Dynamik. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z98_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
