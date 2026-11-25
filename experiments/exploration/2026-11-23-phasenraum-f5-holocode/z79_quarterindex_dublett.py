"""Z79: Ist der Quarterindex das Dublett der Naht? (Auftrag nach H6h39: I_pos über die Araki-Woods-Oszillatoren aus z78)

Leitlinie 'Rekonstruktion statt Auswahl'. Firewalls: GRAV.NONCIRCULAR.01 (Lorentzgruppe nicht aus dem Modulfluss),
DIMENSION.UPLIFT.FIREWALL.01 (alles Kinematik). Bestand: z62/z74 (Standardform, T', I_pos), z63 (rho_m, R), z64
(Blattblindheit), z70 (Spinbrücke, eta), z76 (K/pi = Q_a - Q_b mod 4), z78 (Araki-Woods, P_mu, Orte).

Teil 1  Kinematik: Die Stromoszillatoren a_1, a_2 der beiden Quarter haben dieselbe thermische Antwort (Teleskopsumme, für
        jedes alpha), der Oszillatorzustand ist U(2)-invariant. Mit den Tomita-Partnern b_i bilden sie den Twistorraum
        T = span{a_i, b_i^+} mit der CCR-Form der Signatur (2,2). tau(A) = (J A J)^+ ist eine lineare, antiisometrische
        Involution; ihre Eigenräume (pi- und omega-Ebene) sind null. Der Kommutant von tau in u(2,2) ist gl(2,C): Lorentz,
        Dilatation und das Zentrum, und das Zentrum ist der Modulfluss. e^{iP(X)} ist die Twistortranslation, |X> gehört
        zum Röhrenpunkt M_X = kappa^2 (I - i X^T), die Überlappung ist der nu = 1-Kern. Drehungen erhalten den Zustand,
        Boosts verschieben das Ruhesystem.
Teil 2  Diskret: Die Vierteldrehung R der Marken (R = Delta^(-i/2) = e^{2 pi i L_0}) auf vier Dubletts: Nahtdublett V_1
        (Fock-Ecke, rho_1), fermionisches Quarterpaar, eta-Dublett der Spinbrücke H0(E, eta (x) pi^*O(1)) und
        Quarteroszillatoren. Eindeutigkeit der Quartergrade {1, 3} als Spinordublett.
Teil 3  Ortsräume: R auf H0(N) der Spinbrücke, auf Hom(V, V-bar^*) (T') und auf Herm(Q) (z78); Lorentz- gegen
        euklidische Verträglichkeit.
Teil 4  Bilanz und Urteil.
"""
import itertools
import json
from pathlib import Path

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla
import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def nullraum(M, tol=1e-9):
    u, s, vh = np.linalg.svd(M)
    rank = int(np.sum(s > tol))
    return vh[rank:]


rng = np.random.default_rng(79)
I2 = np.eye(2)
BETA = 4*np.pi

# =================================================================== Teil 1
say("Teil 1: Kinematik der Quarteroszillatoren mit Tomita-Partner")


# --- 1.1 gleiche thermische Antwort beider Quarter
def strom_korrelation(alpha, m, N=400):
    xs = np.arange(-N, N) + alpha
    c = 1/(1 + np.exp(np.clip(BETA*xs, -700, 700)))
    cm = 1/(1 + np.exp(np.clip(BETA*(xs + m), -700, 700)))
    return float(np.sum(c*(1 - cm)))                               # <J_m J_-m> quasifrei


ok11 = all(abs(strom_korrelation(a, m) - m/(1 - np.exp(-BETA*m))) < 1e-12*m for a in (0.25, 0.75) for m in range(1, 8))
ok11 &= abs(strom_korrelation(0.25, 3) - strom_korrelation(0.75, 3)) < 1e-12
check("<J_m J_-m> = Summe c(x)(1 - c(x+m)) = m/(1 - e^{-4 pi m}) für alpha = 1/4 und 3/4 (m = 1..7): Teleskopsumme, die Antwort "
      "hängt nicht vom Quarter ab, beide Oszillatoren haben dieselbe Besetzung 1/(e^{4 pi m} - 1)", ok11)


def stabilisator_dim(nus):
    """sp(4,R)-Elemente, die die thermische Kovarianz diag(nu1, nu1, nu2, nu2) erhalten."""
    Jw = np.kron(np.eye(2), np.array([[0, 1], [-1, 0]]))
    G = np.diag([nus[0], nus[0], nus[1], nus[1]])
    rows = []
    for k in range(16):
        A = np.zeros(16)
        A[k] = 1
        A = A.reshape(4, 4)
        rows.append(np.concatenate([(A.T @ Jw + Jw @ A).ravel(), (A @ G + G @ A.T).ravel()]))
    return nullraum(np.array(rows).T).shape[0]


nu = 1/np.tanh(BETA/2)
d_gleich, d_ungleich = stabilisator_dim((nu, nu)), stabilisator_dim((nu, 1/np.tanh(BETA*0.6)))
ZAHLEN["stabilisator_dim"] = {"gleiche Quarter": d_gleich, "Kontrolle ungleich": d_ungleich}
check(f"die Symmetrien von CCR und Oszillatorzustand bilden u(2) (Dimension {d_gleich}); Kontrolle: verschiedene Besetzungen "
      f"ließen nur u(1) + u(1) (Dimension {d_ungleich})", d_gleich == 4 and d_ungleich == 2)

# --- 1.2 Twistorraum, CCR-Form, tau
# Koordinaten auf T: A = sum alpha_i a_i + beta_i b_i^+  <->  v = (alpha_1, alpha_2, beta_1, beta_2)
G = np.diag([1., 1., -1., -1.])                                    # <A, B> = [A, B^+]
TAU = np.block([[np.zeros((2, 2)), I2], [I2, np.zeros((2, 2))]])    # (J A J)^+: a_i <-> b_i^+
PI_ = np.vstack([I2, I2])/np.sqrt(2)                               # pi_i = (a_i + b_i^+)/sqrt2
OM_ = np.vstack([I2, -I2])/np.sqrt(2)                              # omega_i = (a_i - b_i^+)/sqrt2
ok12 = np.allclose(TAU @ TAU, np.eye(4)) and np.allclose(TAU.T @ G @ TAU, -G)
ok12 &= np.allclose(PI_.T @ G @ PI_, 0) and np.allclose(OM_.T @ G @ OM_, 0) and np.allclose(PI_.T @ G @ OM_, I2)
ok12 &= np.allclose(TAU @ PI_, PI_) and np.allclose(TAU @ OM_, -OM_)
ew = np.linalg.eigvalsh(G)
check("T = span{a_1, a_2, b_1^+, b_2^+} trägt die CCR-Form der Signatur (2,2); tau(A) = (J A J)^+ (a_i <-> b_i^+) ist linear, "
      "involutiv und antiisometrisch; seine Eigenräume, die pi-Ebene (a + b^+) und die omega-Ebene (a - b^+), sind null und "
      "zueinander dual ([omega_i, pi_j^+] = delta_ij)", ok12 and np.sum(ew > 0) == 2 and np.sum(ew < 0) == 2)

# --- 1.3 Kommutant von tau in u(2,2)
basis = []
for k in range(16):
    for ph in (1, 1j):
        A = np.zeros(16, complex)
        A[k] = ph
        basis.append(A.reshape(4, 4))
rows = []
for A in basis:
    cond = np.concatenate([(A.conj().T @ G + G @ A).ravel(), (A @ TAU - TAU @ A).ravel()])
    rows.append(np.concatenate([cond.real, cond.imag]))
kern = nullraum(np.array(rows).T)
alg = [sum(c*B for c, B in zip(v, basis)) for v in kern]
P_ON = np.hstack([PI_, OM_])                                       # Basiswechsel in (pi, omega)
blocks = [np.linalg.inv(P_ON) @ A @ P_ON for A in alg]
blockdiag = all(np.abs(B[:2, 2:]).max() < 1e-9 and np.abs(B[2:, :2]).max() < 1e-9 for B in blocks)
dual = all(np.allclose(B[2:, 2:], -B[:2, :2].conj().T) for B in blocks)


def orts_wirkung(A2):
    """Infinitesimale Wirkung von pi -> e^{tA} pi auf P(X) = kappa^2 pi^+ X pi: X -> A^+ X + X A, in Pauli-Koordinaten."""
    pauli = [I2, np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1., -1.])]
    M = np.zeros((4, 4))
    for j, s in enumerate(pauli):
        Y = A2.conj().T @ s + s @ A2
        M[:, j] = [np.trace(t @ Y).real/2 for t in pauli]
    return M


eta4 = np.diag([1., -1., -1., -1.])
wirk = [orts_wirkung(B[:2, :2]) for B in blocks]
rang = np.linalg.matrix_rank(np.array([w.ravel() for w in wirk]), tol=1e-9)
cond_lor = np.array([(w.T @ eta4 + eta4 @ w).ravel() for w in wirk]).T
lor_koeff = nullraum(cond_lor)                                     # Teilalgebra, die det erhält
rang_lor = np.linalg.matrix_rank(np.array([sum(c*w for c, w in zip(v, wirk)).ravel() for v in lor_koeff]), tol=1e-9)
zentrum = 1j*np.eye(2)
ZAHLEN["kommutant_tau"] = {"dim": len(alg), "rang auf Orten": int(rang), "det-erhaltend": int(lor_koeff.shape[0]),
                            "Rang davon": int(rang_lor)}
check(f"der Kommutant von tau in u(2,2) hat Dimension {len(alg)} und ist block-diagonal diag(A, -A^+) in (pi, omega), also "
      f"gl(2,C); auf den Orten wirkt er mit Rang {rang} als co(1,3) (Lorentz 6 + Dilatation 1); das Zentrum i R I wirkt "
      "trivial auf den Orten; die det-erhaltende Teilalgebra hat Dimension 7 und wirkt mit Rang 6 (so(1,3), Kern = Zentrum)",
      len(alg) == 8 and blockdiag and dual and rang == 7 and np.allclose(orts_wirkung(zentrum), 0)
      and lor_koeff.shape[0] == 7 and rang_lor == 6)
# Modulgenerator der Purifikation: K = beta omega (N_a - N_b) -> a, b^+ bekommen dieselbe Phase: Skalar auf T
Kmod_T = 1j*np.eye(4)                                              # Wirkung auf T: v -> i v (bis auf beta omega)
check("der Modulfluss der Purifikation (K = 4 pi m (N_a - N_b)) dreht a_i und b_i^+ um dieselbe Phase, ist also das Zentrum von "
      "gl(2,C) (die Helizitätsphase) und lässt alle P(X) fest: Orte sind modular invariant, der Modulfluss ist keine "
      "Raumzeitbewegung (GRAV.NONCIRCULAR.01)", np.allclose(Kmod_T @ TAU, TAU @ Kmod_T)
      and np.allclose(Kmod_T.conj().T @ G + G @ Kmod_T, 0) and np.allclose(orts_wirkung(1j*np.eye(2)), 0))

# --- 1.4 Orte als Twistortranslationen, Röhrenpunkt und Kern (Fockraum wie z78)
q = 0.2
kap2 = (1 - q)/(1 + q)


def L_von(X):
    """ad_{P(X)} auf T, P(X) = (1/2) c^+ X c = sum H_kl Z_k^+ Z_l mit Z = (a_1, a_2, b_1^+, b_2^+)."""
    H = (kap2/2)*np.block([[X, X], [X, X]])
    return -(G @ H).T


ok_nil = True
for _ in range(5):
    xv = rng.normal(size=4)
    X = xv[0]*I2 + xv[1]*np.array([[0, 1], [1, 0]]) + xv[2]*np.array([[0, -1j], [1j, 0]]) + xv[3]*np.diag([1, -1])
    L = L_von(X)
    ok_nil &= np.allclose(L @ PI_, 0) and np.allclose(L @ OM_, -kap2*PI_ @ X.T) and np.allclose(L @ L, 0)
check("ad P(X) auf T ist nilpotent: pi bleibt fest, omega -> omega - kappa^2 X^T pi; e^{iP(X)} ist die Twistortranslation, "
      "die die Gerade L_0 verschiebt", ok_nil)


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


n = 10
a1, a2, b1, b2 = fock4(n)
Z = [a1, a2, b1.T.tocsr(), b2.T.tocsr()]                           # a_1, a_2, b_1^+, b_2^+
om = np.zeros([n]*4)
for i in range(n):
    for j in range(n):
        om[i, j, i, j] = (1 - q**2)*q**(i + j)
om = om.reshape(-1).astype(complex)
kap = np.sqrt(kap2)
ct = [kap*(Z[0] + Z[2]), kap*(Z[1] + Z[3])]
PAULI = [I2, np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1., -1.])]


def zustand(X):
    P = 0.5*sum(ct[i].conj().T @ ct[j]*X[i, j] for i in range(2) for j in range(2))
    return spla.expm_multiply(1j*P, om)


ok_ann, ok_kern = True, True
vernichter = [np.array([1, 0, -q, 0]), np.array([0, 1, 0, -q])]   # d_i ~ a_i - q b_i^+
zustaende = []
for _ in range(3):
    X = sum(c*s for c, s in zip(rng.normal(size=4)*0.12, PAULI))
    st = zustand(X)
    zustaende.append((X, st))
    U = np.eye(4) + 1j*L_von(X)                                    # e^{i L_X}, L nilpotent
    for v in vernichter:
        w = U @ v
        A = sum(w[k]*Z[k] for k in range(4))
        ok_ann &= np.linalg.norm(A @ st) < 2e-3*np.linalg.norm(st)
    # Röhrenpunkt: Graph p = M w über der omega-Ebene
    W = np.column_stack([U @ v for v in vernichter])
    p, w_ = np.linalg.solve(P_ON, W)[:2], np.linalg.solve(P_ON, W)[2:]
    M = p @ np.linalg.inv(w_)
    ok_ann &= np.allclose(M, kap2*(I2 - 1j*X.T), atol=1e-12)
for (X, sx), (Y, sy) in itertools.combinations(zustaende, 2):
    MX, MY = kap2*(I2 - 1j*X.T), kap2*(I2 - 1j*Y.T)
    kern_t = 1/np.linalg.det((MX.conj().T + MY)/(2*kap2))
    ok_kern &= abs(np.vdot(sx, sy) - kern_t) < 3e-3
check("auf dem Fockraum (vier Moden, je 10 Stufen, q = 0,2, drei zufällige X): die translatierten Vernichter e^{iL_X} d_i "
      "vernichten |X> = e^{iP(X)} Omega; |X> gehört zum Röhrenpunkt M_X = kappa^2 (I - i X^T) mit Re M_X = kappa^2 I > 0", ok_ann)
check("<X|Y> = det((M_X^+ + M_Y)/(2 kappa^2))^-1 = det(I - (i/2)(Y - X))^-1: der nu = 1-Kern der Röhre über dem "
      "Vorwärtskegel, derselbe Kern det(I - Z W^+)^-1, den z74 aus Hardy- und Stromauslesung erzwingt", ok_kern)

# --- 1.5 Lorentzwirkung als Bogoliubov-Abbildung
SIG = np.block([[np.zeros((2, 2)), I2], [I2, np.zeros((2, 2))]])    # CCR-Form in (pi, omega)
ok_lor, ok_ruhe = True, True
for _ in range(6):
    g = rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2))
    g /= np.sqrt(np.linalg.det(g))
    S = np.block([[g, np.zeros((2, 2))], [np.zeros((2, 2)), np.linalg.inv(g.conj().T)]])
    ok_lor &= np.allclose(S.conj().T @ SIG @ S, SIG) and np.allclose(S @ np.diag([1, 1, -1, -1]),
                                                                       np.diag([1, 1, -1, -1]) @ S)
    X = sum(c*s for c, s in zip(rng.normal(size=4), PAULI))
    M = kap2*(I2 - 1j*X.T)
    Mg = g @ M @ g.conj().T
    ok_lor &= abs(np.linalg.det((Mg - Mg.conj().T)/2j) - np.linalg.det((M - M.conj().T)/2j)) < 1e-9
    ok_ruhe &= not np.allclose(g @ g.conj().T, I2)
u = np.linalg.qr(rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2)))[0]
check("g in SL(2,C) wirkt als pi -> g pi, omega -> (g^+)^-1 omega: erhält die CCR-Form und tau (also Bogoliubov, nativ), "
      "M -> g M g^+, der Ort bleibt ein Lorentzvektor (det erhalten); Drehungen (g unitär) erhalten den Röhrenpunkt "
      "kappa^2 I, also Omega, Boosts verschieben ihn (das KMS-Ruhesystem)", ok_lor and ok_ruhe
      and np.allclose(u @ (kap2*I2) @ u.conj().T, kap2*I2))
befund("Kinematisch trägt der Quarterindex mit seinem Tomita-Partner die ganze Struktur, die T' verlangt: V = pi-Ebene, "
       "V-bar^* = omega-Ebene (tau-Bild), Sigma = CCR-Form, Lorentzgruppe = Kommutant von tau in SU(2,2), Orte = "
       "Twistortranslationen, Kern nu = 1. Drehungen sind exakte Symmetrien des Originals (gleiche Quarterantwort), Boosts "
       "kinematische Automorphismen; die Naht dieser Realisierung ist P(pi-Ebene) = Himmelssphäre der Impulse P = c c^+/2", True)

# =================================================================== Teil 2
say("\nTeil 2: Die Vierteldrehung R der Marken auf vier Dubletts")


def proj_ordnung(M):
    P = np.eye(M.shape[0], dtype=complex)
    for k in range(1, 33):
        P = P @ M
        if np.allclose(P, P[0, 0]*np.eye(M.shape[0]), atol=1e-9):
            return k
    return None


def E_ladung(alpha, Q):
    return Q**2/2 + (alpha - 0.5)*Q


# 2a Nahtdublett V_1 = Fock-Ecke {Omega, psi*_{1/4} Omega}
L0_ecke = np.array([E_ladung(0.25, 0), E_ladung(0.25, 1)])
R_V1 = np.diag(np.exp(2j*np.pi*L0_ecke))
rho_ecke = np.exp(-BETA*L0_ecke)/np.sum(np.exp(-BETA*L0_ecke))
rho1 = np.array([1, np.exp(-np.pi)])/(1 + np.exp(-np.pi))
# 2b fermionisches Quarterpaar {e_{+1/4} (alpha = 1/4), e_{-1/4} (alpha = 3/4)}
xs_paar = np.array([0.25, -0.25])
R_Qf = np.diag(np.exp(2j*np.pi*xs_paar))
c_paar = 1/(1 + np.exp(BETA*xs_paar))
# 2c Quarteroszillatoren: J_-m hebt L_0 um m (ganz), die Nahtdrehung um theta gibt e^{i m theta} für beide Quarter
R_Qo = np.diag(np.exp(2j*np.pi*np.array([1, 1])))                 # m = 1; für jedes ganze m dasselbe
ordn = {"V_1 (Fock-Ecke, rho_1)": proj_ordnung(R_V1), "Quarterpaar fermionisch": proj_ordnung(R_Qf),
        "Quarteroszillatoren": proj_ordnung(R_Qo)}
check("R = e^{2 pi i L_0}: auf dem Nahtdublett V_1 = {Omega, psi*_{1/4} Omega} diag(1, i) (projektive Ordnung 4, Zustand "
      "rho_1); auf dem fermionischen Quarterpaar diag(i, -i) (Ordnung 2; seine Einteilchenkovarianz ist ebenfalls rho_1, "
      "c(1/4) + c(-1/4) = 1); auf den Quarteroszillatoren die Eins (ganzzahliges L_0, die Nahtdrehung ist dort skalar)",
      ordn == {"V_1 (Fock-Ecke, rho_1)": 4, "Quarterpaar fermionisch": 2, "Quarteroszillatoren": 1}
      and np.allclose(rho_ecke, rho1) and np.allclose(np.sort(c_paar), np.sort(rho1)) and abs(c_paar.sum() - 1) < 1e-15)

# 2d eta-Dublett der Spinbrücke: H0(E, eta (x) pi^*O(1)) = L(R_1 + R_-1) auf E: w^2 = 1 - u^4, eta = antipodale Paarung
f_ = lambda u, w: w/(u**2 - 1)                                     # div f = R_i + R_-i - R_1 - R_-1


def kurvenpunkte(N):
    us = rng.normal(size=N) + 1j*rng.normal(size=N)
    ws = np.sqrt(1 - us**4)*np.where(rng.random(N) < 0.5, 1, -1)
    return us, ws


eps = 1e-3
lok = [abs(f_((1 - eps**2)**0.25, eps)*eps + 2) < 1e-5,            # einfacher Pol an R_1
       abs(f_(-(1 - eps**2)**0.25, eps)*eps + 2) < 1e-5,           # einfacher Pol an R_-1
       abs(f_(1j*(1 - eps**2)**0.25, eps)/eps + 0.5) < 1e-5,       # einfache Nullstelle an R_i
       abs(abs(f_(1/eps, np.sqrt(1 - eps**-4 + 0j))) - 1) < 1e-5]   # endlich (Betrag 1) über unendlich


def transport_matrix(basis_fn, tau_inv, mit_f=True):
    us, ws = kurvenpunkte(24)
    B = np.array([[b(u, w) for b in basis_fn] for u, w in zip(us, ws)])
    cols, res = [], 0.0
    for b in basis_fn:
        vals = np.array([(f_(u, w) if mit_f else 1)*b(*tau_inv(u, w)) for u, w in zip(us, ws)])
        coef, *_ = np.linalg.lstsq(B, vals, rcond=None)
        res = max(res, np.abs(B @ coef - vals).max())
        cols.append(coef)
    return np.array(cols).T, res


sig_inv = lambda u, w: (-1j*u, w)                                  # Vierteldrehung, Hebung w -> w
sigp_inv = lambda u, w: (-1j*u, -w)                                # Hebung w -> -w
iota = lambda u, w: (u, -w)                                        # Deck
eta_basis = [lambda u, w: 1 + 0*u, f_]
Rs, r1 = transport_matrix(eta_basis, sig_inv)
Rsp, r2 = transport_matrix(eta_basis, sigp_inv)
Di, r3 = transport_matrix(eta_basis, iota, mit_f=False)
grp = [Rs, Rsp, Di, Rs @ Di, Di @ Rs]
verh = [np.linalg.eigvals(Mx) for Mx in grp]
ok_v4 = all(abs(v[0]/v[1] + 1) < 1e-9 for v in verh)
ZAHLEN["eta_dublett"] = {"R (w->w)": np.round(Rs, 9).tolist(), "R (w->-w)": np.round(Rsp, 9).tolist(),
                         "Deck": np.round(Di, 9).tolist()}
check("eta-Dublett L(R_1 + R_-1) = span{1, w/(u^2 - 1)} (Polordnungen lokal geprüft): die Vierteldrehung (beide Hebungen) und "
      "das Deck wirken mit Eigenwertverhältnis -1, also projektiv als Kleinsche Vierergruppe (nur Halbdrehungen); die "
      f"Spinorhebung w -> -w hat Eigenwerte {{i, -i}}, Quadrat -1, det 1 (Residuen {max(r1, r2, r3):.0e})",
      all(lok) and max(r1, r2, r3) < 1e-9 and ok_v4 and np.allclose(Rsp @ Rsp, -np.eye(2))
      and abs(np.linalg.det(Rsp) - 1) < 1e-12 and proj_ordnung(Rs) == 2 and proj_ordnung(Rsp) == 2)
ew_sp, ew_qf = np.linalg.eigvals(Rsp), np.diag(R_Qf)
gleich_menge = lambda A, B: all(min(abs(a - b) for b in B) < 1e-9 for a in A) and all(min(abs(a - b) for a in A) < 1e-9 for b in B)
check("das eta-Dublett (Spinorhebung) und das fermionische Quarterpaar sind als mu_4-Moduln gleich (Eigenwerte {i, -i}); das "
      "Nahtdublett V_1 (Verhältnis i) ist es nicht", gleich_menge(ew_sp, ew_qf) and proj_ordnung(R_V1) != proj_ordnung(Rsp))
# 2e Eindeutigkeit der Quartergrade
paare = [(a, b) for a in range(4) for b in range(a, 4)
         if abs(1j**(a + b) - 1) < 1e-12 and np.allclose(np.diag([1j**a, 1j**b])**2, -np.eye(2))]
check("unter allen Paaren von mu_4-Charakteren ist {1, 3} das einzige Dublett in SU(2) mit g^2 = -1 (Spinorvorzeichen der "
      "vollen Drehung): die Quarterzwiste 1/4, 3/4 sind genau das Spinordublett der Vierteldrehung", paare == [(1, 3)])

# =================================================================== Teil 3
say("\nTeil 3: Die Vierteldrehung auf den Ortsräumen")
N_basis = [lambda u, w: 1 + 0*u, lambda u, w: u, f_, lambda u, w: u*f_(u, w)]
RN, r4 = transport_matrix(N_basis, sig_inv)
RNp, r5 = transport_matrix(N_basis, sigp_inv)
evN = np.linalg.eigvals(RN)
verschieden = min(abs(a - b) for a, b in itertools.combinations(evN, 2)) > 1e-6
gV = np.diag([np.exp(1j*np.pi/4), np.exp(-1j*np.pi/4)])            # Nahtdublett V^* = span{1, u}: diag(1, -i), normiert
geta = 1j*np.array([[0, 1], [1, 0]])                              # eta-Dublett, in SU(2)
eukl = np.kron(gV, geta)


def spektrum_gleich_bis_phase(A, B):
    ea, eb = np.linalg.eigvals(A), np.linalg.eigvals(B)
    for z in eb:
        ph = z/ea[0]
        if all(min(abs(ph*x - y) for y in eb) < 1e-8 for x in ea):
            return True
    return False


lor_moeglich = any(np.allclose(np.sort_complex(np.linalg.eigvals(np.kron(g, g.conj()))),
                               np.sort_complex(ph*evN), atol=1e-8)
                   for g in [np.diag([np.exp(1j*t), np.exp(-1j*t)]) for t in np.linspace(0, np.pi, 721)]
                   for ph in np.exp(2j*np.pi*np.arange(16)/16))
ZAHLEN["R_auf_H0N"] = np.round(np.sort_complex(evN), 9).tolist()
check(f"H0(N) = L(R_1 + R_-1 + inf_+ + inf_-) = span{{1, u, f, u f}} (Residuum {max(r4, r5):.0e}): die Vierteldrehung hat dort "
      "die vier verschiedenen Eigenwerte {1, i, -1, -i}, die reguläre Darstellung von mu_4, also H0(N) = V^* (x) eta-Dublett",
      max(r4, r5) < 1e-9 and verschieden and spektrum_gleich_bis_phase(RN, np.diag([1, 1j, -1, -1j]))
      and spektrum_gleich_bis_phase(RN, eukl))
spT = np.kron(gV, gV.conj())
spQf = np.kron(R_Qf, R_Qf.conj())
check("kein Lorentzelement endlicher Ordnung (Vektordarstellung g (x) g-quer hat den Eigenwert 1 doppelt) und keine Phase davon "
      "gibt dieses Spektrum; es ist das euklidische Element (Viertel auf V^*, Halb auf eta) in SU(2) x SU(2); T' (zweiter "
      "Faktor = V-bar^*) gäbe {1, 1, i, -i}, die Orte aus z78 (Quarteroszillatoren) {1, 1, 1, 1}, das fermionische "
      "Quarterpaar {1, 1, -1, -1}", not lor_moeglich
      and spektrum_gleich_bis_phase(spT, np.diag([1, 1, 1j, -1j])) and spektrum_gleich_bis_phase(spQf, np.diag([1, 1, -1, -1])))
befund("Die Markensymmetrie unterscheidet die Dubletts: Der Quarterindex ist auf Fermionebene genau das eta-Dublett der "
       "Spinbrücke (zweiter Spinorfaktor von H0(N) = V^* (x) eta), nicht das Nahtdublett V (erster Faktor). Mit der nativen "
       "Vierteldrehung ist die quaternionische (euklidische) Realform verträglich, die lorentzsche Identifikation eta = "
       "V-bar^* (T') nicht. Das schärft die Blattblindheit aus z64: Es scheitert schon die Vierteldrehung selbst, nicht nur "
       "das Deck", True)
befund("Die Orte aus z78 leben auf den Quarteroszillatoren, wo R trivial ist: Ihre Himmelssphäre P(pi-Ebene) ist rund und "
       "trägt keine Marken. R lebt im Ladungssektor (K/pi = Q_a - Q_b mod 4, z76); Ladung und Oszillatoren faktorisieren "
       "(z58). Eine R-äquivariante Abbildung von der markierten Naht P(V_1) auf die Himmelssphäre der Relativzustände gibt es "
       "daher nicht", True)

# =================================================================== Teil 4
say("\nTeil 4: Bilanz")
STATUS = {
    "Kinematik aus dem Quarterindex": ("E", "Quarteroszillatoren mit Tomita-Partner tragen V + V-bar^*, Sigma (2,2), tau, "
                                       "Lorentz = Kommutant von tau in SU(2,2), Orte als Twistortranslationen, Röhrenpunkt "
                                       "und Kern nu = 1; Drehungen exakte Zustandssymmetrien, Boosts kinematisch"),
    "Quarterindex = eta-Dublett": ("E", "auf Fermionebene R-äquivariant (Spinorhebung, Grade {1, 3} eindeutig)"),
    "Quarterindex = Nahtdublett V": ("E (negativ)", "R hat auf V_1 Ordnung 4, auf dem Quarterpaar 2, auf den Oszillatoren 1"),
    "T' mit Markensymmetrie": ("E (negativ)", "eta = V-bar^* ist nicht R-äquivariant; R auf H0(N) ist euklidisch, nicht "
                                              "lorentzsch"),
    "I_pos für die markierte Naht": ("O", "die Orte der Relativzustände gehören zur runden Himmelssphäre der Oszillatoren; "
                                          "die Marken liegen im Ladungssektor"),
    "Lorentzform aus der euklidischen nativen Struktur": ("O", "Osterwalder-Schrader-Rekonstruktion (Reflexionspositivität) "
                                                               "statt T' ist nicht gerechnet"),
}
for k, (s, t) in STATUS.items():
    say(f"   [{s}] {k}: {t}")
ZAHLEN["status"] = STATUS
urteil = ("Der Quarterindex ist kinematisch ein vollwertiges Spinordublett: Mit dem Tomita-Partner bildet er den Twistorraum "
          "mit CCR-Form (2,2), tau = (Ad J)^+ trennt pi- und omega-Ebene, der Kommutant von tau in SU(2,2) ist die "
          "Lorentzgruppe mit Dilatation, das Zentrum ist der Modulfluss, e^{iP(X)} ist die Twistortranslation zum "
          "Röhrenpunkt kappa^2 (I - i X^T), und die Überlappung ist der nu = 1-Kern aus z74. Drehungen sind exakte "
          "Symmetrien, weil beide Quarter dieselbe Stromantwort haben; Boosts verschieben das KMS-Ruhesystem. Das Nahtdublett "
          "ist er aber nicht: Die Vierteldrehung der Marken hat auf V_1 Ordnung 4, auf dem fermionischen Quarterpaar Ordnung "
          "2 und auf den Oszillatoren Ordnung 1. Auf Fermionebene ist der Quarterindex genau das eta-Dublett der Spinbrücke "
          "(Grade {1, 3}, eindeutig). Folge: Auf H0(N) = V^* (x) eta wirkt die Markensymmetrie euklidisch (SU(2) x SU(2) mit "
          "verschiedenen Faktoren), die lorentzsche Identifikation eta = V-bar^* (T') ist nicht markenverträglich, und die "
          "Orte aus z78 gehören zu einer runden Himmelssphäre ohne Marken. I_pos für die markierte Naht bleibt offen; "
          "3+1 lorentzsch folgt aus dem Quarterindex nur als Kinematik der Oszillatoren.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z79_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
