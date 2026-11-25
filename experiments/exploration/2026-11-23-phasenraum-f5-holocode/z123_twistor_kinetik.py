"""Z123: Twistorkinetik aus der Tomita-Paarung? Das Austauschfeld (a, b) zwischen Nahtgeraden (Folgeschritt zu z111/T2).

Leitlinie 'Rekonstruktion statt Auswahl': Die holomorphe BF-Kinetik tr(b dbar a), der CSW-Propagator dbar^-1 und die
CSW-Regeln sind VERGLEICHSZIELE, keine Prämissen. Gefragt ist, ob native Strukturen die Kinetik des in z111 benannten
Austauschfeldes liefern: a (Homogenität 0, Fock-Sektor D = N_d - N_e = -2) und b (Homogenität -4, D = +2), verbunden durch
Tomitas n -> -n - 4 (z85). Kandidat aus T2: die bilineare Tomita-Paarung omega(psi, phi) = <J psi|phi>.
Firewall DIMENSION.UPLIFT.FIREWALL.01: 4D tritt nur als Lesart auf (Naht = Twistorgerade L_x, Hypothese T aus z110) und als
Kinematik freier Einteilchenzustände; es wird keine wechselwirkende 3+1-Schließung behauptet. H6 und die gelösten Blöcke
(Flavour, Massen, alpha, 16, E8) bleiben fest, nichts wird befördert. Statuskürzel: E exakt, C bedingt, L Lesart, O offen.

Bestand: z85 (vier Oszillatoren d_1, d_2, e_1, e_2; Z = (d_1, d_2, e_1^+, e_2^+), Z-quer = (d_1^+, d_2^+, -e_1, -e_2);
Homogenität n = -D - 2; Čech-Bild |n1 n2; m1 m2> -> sqrt(n1! n2!/(m1! m2!)) z3^m1 z4^m2/(z1^(n1+1) z2^(n2+1)); Tomita J =
Tausch d <-> e mit Konjugation, J D J = -D; Delta^(1/2) = q^D, nativ q = e^(-2 pi)), z114 (dieselbe Leiter als Standardform
z83/z97, Konvention n = D - 2 = J-Bild; Z_SU = exp(i pi (E + 2 J3_R)) = i^n), z84 (Theta_n = [[0, Q(n)], [Q(n)^+, 0]],
Theta_1 = tau aus z79, antiisometrisch für diag(1, 1, -1, -1); Theta_n = U(2) x U(2)-Bahn von tau; radiale Inversion und
thermischer 4D-Kern bei beta = 4 pi reflexionspositiv), z111 (BCFW/CSW-
Vergleichsmaschine; Referenz eta aus I_pi; Cauchykern -1/(<ab> det(x - y)) = Penrose-Transformation von 1/((A.Z)(B.Z));
euklidisch 4 pi^2/P^2; Normierung dort nicht nativ), z110 (Nahtströme: Einspuranteil -k mal Parke-Taylor, auf der Geraden
1/(z12...zn1) = (-1)^n/(<12>...<n1>); g_eff = k^(-1/2) unter P_norm; Zähler <ij>^4 nur unter P_N4).

Teil 1  Tomita-Paarung zwischen a (D = -2) und b (D = +2): Modell (gl(4)-Relationen), Nichtentartung je Stufe, exakte
        Operatoridentität J E J = -(tau E tau)^+ für alle 16 Weyl-geordneten Erzeuger, naive Invarianz J X J = -X^+ und das
        No-go, Realform (tau = Theta_1 = Punktspiegelung x -> -x, Fixalgebra Lorentz + Dilatation + u(1)), Klassifikation
        sigma-verdreht invarianter Bilinearformen (sigma = id, Ad eta, Ad tau, Ad Theta_n), Čech-Bild omega = Residuen-
        paarung o tau^* und ihre sl(4)-Invarianz.
Teil 2  Off-shell-Fortsetzung über die OS-Struktur: kohärente Zustände (Fockreihe gegen geschlossene Formen), ihr Čech-Bild
        = Elementarzustand 1/((A.Z)(B.Z)) von z111; ein Kern G, zwei Spiegelungen: die OS-Spiegelung theta (Röhrenkonjugation,
        euklidisch die radiale Inversion aus z84) gibt die Hilbertform, tau = theta J die Tomita-Paarung; Penrose-
        Transformation = Tomita-Paarung mit dem tau-gespiegelten Punkt,
        Röhrenform und Translationsinvarianz, CSW-Eichung dbar^-1 zwischen zwei Geraden (delta-quer-Lokalisierung) gegen den
        a-b-Kern, Impulsraum, Normierung durch Faktorisierung und die nativen Kandidaten (Niveau, KMS, Uhr, Tomita-Phase).
Teil 3  Zusammenbau S = int tr(b dbar a) + Nahtvertizes (log det auf den Geraden): CSW gegen BCFW für alle 4- und 5-Punkt-
        MHV- und NMHV-Lagen, Negativkontrollen (Hilbert-Paarung, ungespiegelte Paarung, Referenz je Propagator); Bilanz,
        Urteil, Folgen für K15 und die Szenarioachse Wechselwirkung.
"""
import itertools
import json
import math
import time
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.sparse as sps
import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
rng = np.random.default_rng(123)
Q_NAT = math.exp(-2*math.pi)                                                           # q = e^(-2 pi m), m = 1 (z85)
C14 = 1/(1 + math.exp(math.pi))                                                         # c(1/4) bei beta = 4 pi (z110)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def r(x, n=12):
    return float(np.round(float(x), n))


def laufzeit():
    say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")


def cz(*shape):
    return rng.normal(size=shape) + 1j*rng.normal(size=shape)


# ------------------------------------------------------------------ Fockmodell (z85)
def sektor(D, K):
    """Fockzustände (n1, n2, m1, m2) mit N_d - N_e = D und Stufe k = min(N_d, N_e) <= K, nach Stufe sortiert."""
    out = []
    for k in range(K + 1):
        Nd, Ne = (k + D, k) if D >= 0 else (k, k - D)
        for n1 in range(Nd + 1):
            for m1 in range(Ne + 1):
                out.append((n1, Nd - n1, m1, Ne - m1))
    return out


stufe = lambda st: min(st[0] + st[1], st[2] + st[3])
tausch = lambda st: (st[2], st[3], st[0], st[1])                                       # J: d <-> e (dazu Konjugation)
D_von = lambda st: st[0] + st[1] - st[2] - st[3]
TAU = [2, 3, 0, 1]                                                                      # tau: Index alpha <-> alpha + 2
ETA = np.diag([1., 1., -1., -1.])
TAU_M = np.zeros((4, 4))
for _a in range(4):
    TAU_M[_a, TAU[_a]] = 1.0


def E_wirk(al, be, st):
    """Z^al Z-quer_be auf |n1 n2; m1 m2> (erst Z-quer_be, dann Z^al)."""
    s = list(st)
    c = 1.0
    if be < 2:
        c *= math.sqrt(s[be] + 1)
        s[be] += 1
    else:
        if s[be] == 0:
            return None
        c *= -math.sqrt(s[be])
        s[be] -= 1
    if al < 2:
        if s[al] == 0:
            return None
        c *= math.sqrt(s[al])
        s[al] -= 1
    else:
        c *= math.sqrt(s[al] + 1)
        s[al] += 1
    return tuple(s), c


def E_matrix(al, be, basis, idx):
    """Weyl-geordneter Erzeuger E^al_be = (Z^al Z-quer_be + Z-quer_be Z^al)/2 auf einer Sektorbasis (Ausgänge außerhalb der
    Basis fallen weg; geprüft wird nur, wo das nicht vorkommt)."""
    rows, cols, vals = [], [], []
    for j, st in enumerate(basis):
        res = E_wirk(al, be, st)
        if res is not None and res[0] in idx:
            rows.append(idx[res[0]])
            cols.append(j)
            vals.append(res[1])
        if al == be:
            rows.append(j)
            cols.append(j)
            vals.append(-0.5)
    return sps.csr_matrix((vals, (rows, cols)), shape=(len(basis), len(basis)))


def hub(al, be):
    """Stufenänderung von E^al_be: +1 (e^+ d^+, Heber), -1 (d e, Senker), 0 sonst."""
    return 1 if (al >= 2 and be < 2) else (-1 if (al < 2 and be >= 2) else 0)


def sigma_koeff(g):
    """sigma_g(E^al_be) = sum_(ga, de) g[al, ga] g^-1[de, be] E^ga_de (Twistortransformation Z -> g Z)."""
    gi = np.linalg.inv(g)
    return {(al, be): {(ga, de): g[al, ga]*gi[de, be] for ga in range(4) for de in range(4)
                       if abs(g[al, ga]*gi[de, be]) > 1e-14} for al in range(4) for be in range(4)}


PAULI = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.diag([1., -1.]).astype(complex)]
Qn = lambda v: v[0]*np.eye(2) + 1j*sum(v[k + 1]*PAULI[k] for k in range(3))
Theta = lambda v: np.block([[np.zeros((2, 2)), Qn(v)], [Qn(v).conj().T, np.zeros((2, 2))]])


# =================================================================== Teil 1
say("Teil 1: Tomita-Paarung omega(psi, phi) = <J psi|phi> zwischen a (D = -2, Homogenität 0) und b (D = +2, Homogenität -4)")
say("  (a) Modell, J und die Sektoren")
K1 = 6
S1 = {D: sektor(D, K1) for D in (-4, -2, 0, 2, 4)}
I1 = {D: {st: i for i, st in enumerate(B)} for D, B in S1.items()}
ok_J, ok_hom = True, True
for D, B in S1.items():
    for st in B:
        t = tausch(st)
        ok_J &= tausch(t) == st and D_von(t) == -D and stufe(t) == stufe(st) and t in I1[-D]
        n85, n85J = -D - 2, -D_von(t) - 2
        n114 = D - 2                                                                    # z114: n = m - mquer - 2 = D - 2
        ok_hom &= n85J == -n85 - 4 and n114 == -n85 - 4
dims = {D: [sum(1 for st in S1[D] if stufe(st) == k) for k in range(K1 + 1)] for D in (-2, 2)}
ok_dims = all(dims[D][k] == (k + 1)*(k + 3) for D in (-2, 2) for k in range(K1 + 1))
A_, B_ = S1[-2], S1[2]
iA, iB = I1[-2], I1[2]
nA, nB = len(A_), len(B_)
EA = {(al, be): E_matrix(al, be, A_, iA) for al in range(4) for be in range(4)}
EB = {(al, be): E_matrix(al, be, B_, iB) for al in range(4) for be in range(4)}
innen = np.array([stufe(st) <= K1 - 2 for st in A_])
komm_def = 0.0
for (al, be), (ga, de) in itertools.product(EA, EA):
    lhs = (EA[(al, be)] @ EA[(ga, de)] - EA[(ga, de)] @ EA[(al, be)]).toarray()
    rhs = ((al == de)*EA[(ga, be)] - (ga == be)*EA[(al, de)]).toarray()
    komm_def = max(komm_def, float(np.abs((lhs - rhs)[:, innen]).max()))
ZAHLEN["modell"] = {"Stufendimensionen D = -2, +2 (k = 0..6)": dims, "gl(4)-Relationen (256 Paare, innen)": komm_def,
                    "Zustände je Sektor bis Stufe 6": nA}
check(f"[E] (a) Modell z85/z114: Die 16 Weyl-geordneten Bilineare E^al_be = (Z^al Z-quer_be + Z-quer_be Z^al)/2 erfüllen auf "
      f"dem Sektor D = -2 die gl(4)-Relationen [E^a_b, E^c_d] = delta^a_d E^c_b - delta^c_b E^a_d (256 Paare, Defekt "
      f"{komm_def:.0e}); J (Tausch d <-> e) ist eine Involution, erhält die Stufe, kehrt D um und bildet die Homogenität "
      "n = -D - 2 auf -n - 4 ab (alle Zustände der Sektoren D = -4..4 bis Stufe 6). a (n = 0) und b (n = -4) sind genau ein "
      f"J-Paar; je Stufe k haben beide (k + 1)(k + 3) Zustände ({dims[2][:4]}...). In der Konvention von z114 (n = D - 2) "
      "ist die Rolle von a und b vertauscht, die Paarung unten ist symmetrisch und davon unabhängig",
      ok_J and ok_hom and ok_dims and komm_def < 1e-12)

say("  (b) Nichtentartung: omega paart a mit b, Stufe für Stufe")
W_AB = sps.csr_matrix((np.ones(nA), (np.arange(nA), [iB[tausch(st)] for st in A_])), shape=(nA, nB))
lvA = np.array([stufe(st) for st in A_])
lvB = np.array([stufe(st) for st in B_])
sing = []
for k in range(K1 + 1):
    blk = W_AB.toarray()[np.ix_(lvA == k, lvB == k)]
    sing.append(np.linalg.svd(blk, compute_uv=False))
aus_stufe = float(np.abs(W_AB.toarray()*(lvA[:, None] != lvB[None, :])).max())
WAA = np.array([[1.0 if tausch(s1) == s2 else 0.0 for s2 in A_] for s1 in A_])
hilbert_quer = sum(1 for st in A_ if st in iB)                                          # gemeinsame Basisvektoren
konj_quer = float(np.abs(np.array([[1.0 if s1 == s2 else 0.0 for s2 in B_] for s1 in A_])).max())
sym_def = float(np.abs(W_AB.toarray() - np.array([[1.0 if tausch(s2) == s1 else 0.0 for s2 in B_] for s1 in A_])).max())
sv_def = max(float(np.abs(s - 1).max()) for s in sing)
ZAHLEN["nichtentartung"] = {"Singulärwerte je Stufe - 1 (max)": sv_def, "Rang je Stufe": [int(np.sum(s > 0.5)) for s in sing],
                            "omega(a, a)": float(np.abs(WAA).max()), "Hilbert <a|b>: gemeinsame Basisvektoren": hilbert_quer,
                            "reine Konjugation psi^T phi zwischen a und b": konj_quer}
check(f"[E] (b) omega(psi, phi) = <J psi|phi> = psi^T W phi ist bilinear und symmetrisch (Defekt {sym_def:.0e}), paart a nur "
      f"mit b (omega(a, a) = {np.abs(WAA).max():.0f}) und nur gleiche Stufen (außerhalb {aus_stufe:.0f}); je Stufe ist der "
      f"Block eine Permutation, alle Singulärwerte 1 (Abweichung {sv_def:.0e}, Ränge {[int(np.sum(s > 0.5)) for s in sing]} = "
      f"volle Dimension): nichtentartet, mit Betrag 1 auf Fock-normierten Zuständen. Negativkontrollen: Die Hilbertform "
      f"<psi|phi> und die reine Konjugation psi^T phi verschwinden zwischen a und b (gemeinsame Basisvektoren {hilbert_quer}, "
      f"Matrix {konj_quer:.0f}); erst J verbindet die beiden Sektoren",
      sym_def == 0 and np.abs(WAA).max() == 0 and aus_stufe == 0 and sv_def < 1e-12
      and all(int(np.sum(s > 0.5)) == len(s) for s in sing) and hilbert_quer == 0 and konj_quer == 0)
laufzeit()

say("  (c) su(2,2)-Wirkung: die richtige Realform-Aussage")
gueltigB = np.array([stufe(st) <= K1 - 1 for st in B_])
gueltigA = np.array([stufe(st) <= K1 - 1 for st in A_])
Wd = W_AB.toarray()
dreh_def = 0.0
for (al, be) in EA:
    lhsB = Wd.T @ EA[(al, be)].toarray() @ Wd                                           # J E J auf b (reelle Matrizen)
    rhsB = -EB[(TAU[al], TAU[be])].toarray().T
    lhsA = Wd @ EB[(al, be)].toarray() @ Wd.T
    rhsA = -EA[(TAU[al], TAU[be])].toarray().T
    dreh_def = max(dreh_def, float(np.abs((lhsB - rhsB)[:, gueltigB]).max()), float(np.abs((lhsA - rhsA)[:, gueltigA]).max()))
naiv_vek = []
for (al, be) in EA:
    Nm = Wd.T @ EA[(al, be)].toarray() @ Wd + EB[(al, be)].toarray().T                  # J X J + X^+ (naiv)
    naiv_vek.append(Nm[:, gueltigB].ravel())
naiv_M = np.array(naiv_vek).T
sv_n = np.linalg.svd(naiv_M, compute_uv=False)
rang_n = int(np.sum(sv_n > 1e-9*sv_n[0]))
_, _, Vh = np.linalg.svd(naiv_M)
kern = Vh[rang_n:].conj()
tau_sym = max(float(np.abs(x.reshape(4, 4) - TAU_M @ x.reshape(4, 4) @ TAU_M).max()) for x in kern)
einzel_naiv = [float(np.abs(v).max()) for v in naiv_vek]
xH = np.diag([0.5, 0.5, -0.5, -0.5])                                                   # H = (N_d + N_e + 2)/2
H_A = sum(xH[a, a]*EA[(a, a)] for a in range(4)).toarray()
H_B = sum(xH[a, a]*EB[(a, a)] for a in range(4)).toarray()
E_diag = np.diag(H_A)
h_def = (H_A.T @ Wd + Wd @ H_B)                                                         # omega(H psi, phi) + omega(psi, H phi)
h_stufe0 = float(h_def[np.ix_(lvA == 0, lvB == 0)].max())
ZAHLEN["realform"] = {"J E J + (tau E tau)^+ (16 Erzeuger, max)": dreh_def, "naiv J X J = -X^+: Rang, Kern": [rang_n, 16 - rang_n],
                      "Kern tau-symmetrisch (Defekt)": tau_sym, "naiver Defekt je Basis-Erzeuger (min)": min(einzel_naiv),
                      "Energie: omega(H psi, phi) + omega(psi, H phi) auf Stufe 0": h_stufe0,
                      "Energie E auf Stufe 0": float(E_diag[lvA == 0][0])}
check(f"[E] (c) Exakte Operatoridentität für alle 16 Weyl-geordneten Erzeuger: J E^al_be J = -(E^(tau al)_(tau be))^+ mit "
      f"tau = (1 3)(2 4) (Defekt {dreh_def:.0e}, beide Richtungen a <-> b, alle Eingänge bis Stufe 5). Gleichwertig: omega ist "
      "tau-verdreht invariant, omega(E psi, phi) + omega(psi, tau(E) phi) = 0, ohne Zentralterm (die Weyl-Ordnung nimmt "
      "die Verschiebung delta_(al be) auf)", dreh_def < 1e-12)
check(f"[E negativ] (c) Die naive Aussage J X J = -X^+ (also omega(X psi, phi) + omega(psi, X phi) = 0) gilt für keinen der 16 "
      f"Basis-Erzeuger einzeln (min. Defekt {min(einzel_naiv):.2f}) und genau auf einem Unterraum der Dimension {16 - rang_n} "
      f"(Rang {rang_n}), und dieser ist exakt der tau-symmetrische Teil x = tau x tau (Defekt {tau_sym:.0e}). Die Energie "
      f"H = (N_d + N_e + 2)/2 ist tau-ungerade: omega(H psi, phi) + omega(psi, H phi) = 2E omega = {h_stufe0:.0f} auf Stufe 0. "
      "No-go: Zwischen zwei Sektoren positiver Energie gibt es keine su(2,2)-invariante Bilinearform ((E + E') beta = 0 "
      "erzwingt beta = 0); T2s Formulierung J X J = -X^+ ist für su(2,2) unmöglich",
      rang_n == 8 and tau_sym < 1e-10 and min(einzel_naiv) > 0.4 and abs(h_stufe0 - 4) < 1e-12)

# --- Realform: tau = Theta_1, Punktspiegelung, Fixalgebra
KAPPA = math.sqrt(math.tanh(math.pi))
U_PI = np.array([[1, 0, -1, 0], [0, 1, 0, -1]], dtype=complex)                        # z111: span = {pi = 0} = I_pi
V_OM = np.array([[1, 0, 1, 0], [0, 1, 0, 1]], dtype=complex)                          # z111: span = {omega = 0} = I_omega
CPO = np.block([[np.eye(2), np.eye(2)], [KAPPA**2*np.eye(2), -KAPPA**2*np.eye(2)]])/(math.sqrt(2)*KAPPA)   # (pi, omega) = CPO Z
tau_po = CPO @ TAU_M @ np.linalg.inv(CPO)
theta1 = Theta([1, 0, 0, 0])
eig_pi = float(np.abs(TAU_M @ U_PI.T + U_PI.T).max())
eig_om = float(np.abs(TAU_M @ V_OM.T - V_OM.T).max())
x4 = cz(2, 2)
Zx = np.linalg.inv(CPO) @ np.concatenate([cz(2), np.zeros(2)])                         # beliebiges pi, dann omega = i x pi
pi0 = (CPO @ Zx)[:2]
Zx = np.linalg.inv(CPO) @ np.concatenate([pi0, 1j*x4 @ pi0])
Zt = TAU_M @ Zx
punkt = float(np.abs((CPO @ Zt)[2:] - 1j*(-x4) @ (CPO @ Zt)[:2]).max())
# reelle Dimension von u(2,2) und ihrer tau-Teile
basis_r = []
for k in range(32):
    M = np.zeros(16, complex)
    M[k % 16] = 1.0 if k < 16 else 1j
    basis_r.append(M.reshape(4, 4))
real_vec = lambda M: np.concatenate([M.real.ravel(), M.imag.ravel()])
u22 = np.array([real_vec(M.conj().T @ ETA + ETA @ M) for M in basis_r]).T
fix = np.array([real_vec(TAU_M @ M @ TAU_M - M) for M in basis_r]).T
odd = np.array([real_vec(TAU_M @ M @ TAU_M + M) for M in basis_r]).T
dim_u22 = 32 - np.linalg.matrix_rank(u22)
dim_fix = 32 - np.linalg.matrix_rank(np.vstack([u22, fix]))
dim_odd = 32 - np.linalg.matrix_rank(np.vstack([u22, odd]))
# Stabilisator: ein tau-gerades Element erhält I_pi und I_omega, Translation und spezielle konforme bewegen sie
Gpo = {"Translation omega -> omega + a pi": np.block([[np.zeros((2, 2)), np.zeros((2, 2))], [cz(2, 2), np.zeros((2, 2))]]),
       "spezielle konforme pi -> pi + c omega": np.block([[np.zeros((2, 2)), cz(2, 2)], [np.zeros((2, 2)), np.zeros((2, 2))]]),
       "Lorentz/Dilatation (blockdiagonal)": np.block([[cz(2, 2), np.zeros((2, 2))], [np.zeros((2, 2)), cz(2, 2)]])}
paritaet = {}
for k_, G in Gpo.items():
    Gz = np.linalg.inv(CPO) @ G @ CPO
    paritaet[k_] = ("gerade" if np.abs(TAU_M @ Gz @ TAU_M - Gz).max() < 1e-12 else
                    "ungerade" if np.abs(TAU_M @ Gz @ TAU_M + Gz).max() < 1e-12 else "gemischt")
erhalt = lambda Gz, U: np.linalg.matrix_rank(np.hstack([U.T, Gz @ U.T]), tol=1e-9) == 2
GL = np.linalg.inv(CPO) @ Gpo["Lorentz/Dilatation (blockdiagonal)"] @ CPO
GT = np.linalg.inv(CPO) @ Gpo["Translation omega -> omega + a pi"] @ CPO
stab = (erhalt(GL, U_PI) and erhalt(GL, V_OM), erhalt(GT, U_PI), erhalt(GT, V_OM))
ZAHLEN["tau"] = {"tau^2 - 1": float(np.abs(TAU_M @ TAU_M - np.eye(4)).max()), "tau eta tau + eta": float(np.abs(TAU_M @ ETA @ TAU_M + ETA).max()),
                 "tau - Theta_(1,0,0,0)": float(np.abs(TAU_M - theta1).max()), "tau in (pi, omega)": np.round(tau_po.real, 12).tolist(),
                 "I_pi Eigenwert -1 (Defekt)": eig_pi, "I_omega Eigenwert +1 (Defekt)": eig_om, "x -> -x (Defekt)": punkt,
                 "reelle Dimensionen u(2,2), tau-gerade, tau-ungerade": [dim_u22, dim_fix, dim_odd], "Parität der Erzeuger": paritaet,
                 "Lorentz erhält I_pi, I_omega; Translation erhält I_pi; Translation erhält I_omega": stab}
check(f"[E] (c) Realform: tau ist genau Theta_(1,0,0,0) aus z84 (Defekt {np.abs(TAU_M - theta1).max():.0e}), tau^2 = 1 und "
      f"tau eta tau = -eta (antiisometrisch, nicht in U(2,2)). In den z85/z111-Koordinaten ist tau = diag(1, -1) auf (pi, "
      f"omega), also omega -> -omega: die Punktspiegelung x -> -x (Defekt {punkt:.0e}); I_pi ist Eigenraum zu -1, I_omega zu +1 "
      f"(Defekte {eig_pi:.0e}, {eig_om:.0e}). Reelle Dimensionen: u(2,2) {dim_u22} = {dim_fix} (tau-gerade) + {dim_odd} "
      f"(tau-ungerade); tau-gerade sind Lorentz, Dilatation und u(1) (der Stabilisator von I_pi und I_omega: {stab[0]}), "
      f"tau-ungerade Translationen und spezielle konforme Abbildungen ({paritaet}; die Translation erhält I_pi, nicht "
      f"I_omega: {stab[1]}, {stab[2]}). Also ist J die CPT-Abbildung der Leiter: J U(g) J = U(tau g tau), mit U(g) für Lorentz "
      "und Dilatation vertauschend, Translationen umkehrend",
      np.abs(TAU_M - theta1).max() < 1e-15 and np.abs(TAU_M @ ETA @ TAU_M + ETA).max() < 1e-15 and punkt < 1e-12
      and np.abs(tau_po - np.diag([1, 1, -1, -1])).max() < 1e-12 and eig_pi < 1e-15 and eig_om < 1e-15
      and (dim_u22, dim_fix, dim_odd) == (16, 8, 8) and stab == (True, True, False)
      and paritaet == {"Translation omega -> omega + a pi": "ungerade", "spezielle konforme pi -> pi + c omega": "ungerade",
                       "Lorentz/Dilatation (blockdiagonal)": "gerade"})
laufzeit()

say("  (d) Klassifikation: sigma-verdreht invariante Bilinearformen beta zwischen a und b (Stufen <= 3)")
L6 = 3
A6 = [st for st in A_ if stufe(st) <= L6]
B6 = [st for st in B_ if stufe(st) <= L6]
jA6, jB6 = {st: i for i, st in enumerate(A6)}, {st: i for i, st in enumerate(B6)}
n6a, n6b = len(A6), len(B6)
EA6 = {k: E_matrix(*k, A6, jA6) for k in EA}
EB6 = {k: E_matrix(*k, B6, jB6) for k in EA}
lA6 = np.array([stufe(st) for st in A6])
lB6 = np.array([stufe(st) for st in B6])
W6 = np.zeros((n6a, n6b))
for i, st in enumerate(A6):
    W6[i, jB6[tausch(st)]] = 1.0
nv = rng.normal(size=4)
nv /= np.linalg.norm(nv)
TWISTE = {"sigma = id": np.eye(4), "sigma = Ad eta (inner, Isometrie)": ETA,
          "sigma = Ad tau (Theta_1)": TAU_M.astype(complex), "sigma = Ad Theta_n (n zufällig auf S^3)": Theta(nv)}


def verdreht_kern(g):
    """Dimension des Lösungsraums von X_a^T B + B sigma_g(X)_b = 0 (alle 16 Erzeuger, nur Gleichungen ohne Abschneidefehler)
    und der Überlapp |<B, W>|/(|B| |W|) mit der Tomita-Paarung (falls eindimensional)."""
    sk = sigma_koeff(g)
    aus_block = np.abs(g[:2, :2]).max() > 1e-12                                          # blockdiagonal (id, eta) oder nicht
    zeilen = []
    for (al, be), terme in sk.items():
        X = EA6[(al, be)]
        Y = sum((c*EB6[gd] for gd, c in terme.items()), sps.csr_matrix((n6b, n6b), dtype=complex))
        hx = hub(al, be)
        hy = hx if aus_block else -hx
        ok_i = lA6 + hx <= L6
        ok_j = lB6 + hy <= L6
        M = sps.kron(sps.identity(n6b), X.T) + sps.kron(Y.T, sps.identity(n6a))
        maske = (ok_i[:, None] & ok_j[None, :]).ravel(order="F")
        zeilen.append(M.tocsr()[np.where(maske)[0]])
    Ag = sps.vstack(zeilen).tocsr()
    G = (Ag.conj().T @ Ag).toarray()
    ew, ev = np.linalg.eigh(G)
    null = ew < 1e-9*ew[-1]
    dim = int(null.sum())
    ueber = None
    if dim == 1:
        b = ev[:, 0]
        w = W6.ravel(order="F")
        ueber = float(abs(np.vdot(w, b))/(np.linalg.norm(w)*np.linalg.norm(b)))
    luecke = float(ew[dim]/ew[-1]) if dim < len(ew) else 0.0
    return dim, ueber, luecke, Ag.shape


klass = {k: verdreht_kern(g) for k, g in TWISTE.items()}
ZAHLEN["klassifikation"] = {k: {"Dimension": v[0], "Überlapp mit omega": v[1], "Spektrallücke": v[2], "Gleichungen x Unbekannte": list(v[3])}
                            for k, v in klass.items()}
check(f"[E] (d) Klassifikation (Stufen <= {L6}, {n6a} x {n6b} = {n6a*n6b} Unbekannte, alle 16 Erzeuger): sigma-verdreht "
      f"invariante Bilinearformen beta(X psi, phi) + beta(psi, sigma(X) phi) = 0 zwischen a und b gibt es für sigma = id "
      f"{klass['sigma = id'][0]}, für das innere, isometrische Ad eta {klass['sigma = Ad eta (inner, Isometrie)'][0]}, für das "
      f"äußere, eta-umkehrende Ad tau genau {klass['sigma = Ad tau (Theta_1)'][0]}, und die ist die Tomita-Paarung (Überlapp "
      f"{klass['sigma = Ad tau (Theta_1)'][1]:.12f}). Ein anderes Bahnelement Theta_n (z84) gibt ebenfalls eine Form "
      f"(Dimension {klass['sigma = Ad Theta_n (n zufällig auf S^3)'][0]}), aber nicht omega (Überlapp "
      f"{klass['sigma = Ad Theta_n (n zufällig auf S^3)'][1]:.3f}): Die b-a-Kinetik ist bis auf den Faktor eindeutig, sobald "
      "die Spiegelung gewählt ist, und J gehört genau zu Theta_1 = tau",
      klass["sigma = id"][0] == 0 and klass["sigma = Ad eta (inner, Isometrie)"][0] == 0
      and klass["sigma = Ad tau (Theta_1)"][0] == 1 and klass["sigma = Ad tau (Theta_1)"][1] > 1 - 1e-10
      and klass["sigma = Ad Theta_n (n zufällig auf S^3)"][0] == 1 and klass["sigma = Ad Theta_n (n zufällig auf S^3)"][1] < 0.99
      and min(v[2] for v in klass.values()) > 1e-6)
laufzeit()

say("  (e) Čech-Bild: omega = Residuenpaarung nach tau^*")


def cech(st):
    """z85: |n1 n2; m1 m2> -> sqrt(n1! n2!/(m1! m2!)) z3^m1 z4^m2/(z1^(n1+1) z2^(n2+1)); (Exponenten, Koeffizient^2)."""
    n1, n2, m1, m2 = st
    return (-n1 - 1, -n2 - 1, m1, m2), Fr(math.factorial(n1)*math.factorial(n2), math.factorial(m1)*math.factorial(m2))


tau_stern = lambda e: (e[2], e[3], e[0], e[1])                                         # (tau^* f)(z) = f(tau z)
K7 = 4
A7 = [st for st in A_ if stufe(st) <= K7]
B7 = [st for st in B_ if stufe(st) <= K7]
res_def, res_ohne, n_paar = 0.0, 0, 0
for s1 in A7:
    e1, c1 = cech(s1)
    for s2 in B7:
        e2, c2 = cech(s2)
        e2t = tau_stern(e2)
        prod = tuple(a + b for a, b in zip(e1, e2t))
        wert = math.sqrt(c1*c2) if prod == (-1, -1, -1, -1) else 0.0
        res_def = max(res_def, abs(wert - (1.0 if tausch(s1) == s2 else 0.0)))
        res_ohne += int(tuple(a + b for a, b in zip(e1, e2)) == (-1, -1, -1, -1))
        n_paar += 1


def lp_mul(F, G):
    out = {}
    for e1, c1 in F.items():
        for e2, c2 in G.items():
            e = tuple(a + b for a, b in zip(e1, e2))
            out[e] = out.get(e, 0) + c1*c2
    return out


def lp_E(al, be, F):
    """E^al_be im Čech-Bild: -z^al d/dz^be - delta_(al be)/2 (z85: Z = Multiplikation, Z-quer = -d/dz)."""
    out = {}
    for e, c in F.items():
        if e[be] != 0:
            e2 = list(e)
            e2[be] -= 1
            e2[al] += 1
            out[tuple(e2)] = out.get(tuple(e2), 0) - c*e[be]
        if al == be:
            out[e] = out.get(e, 0) - c/2
    return out


res_ = lambda F: F.get((-1, -1, -1, -1), Fr(0))
inv_def, res_nonzero = Fr(0), []
for _ in range(6):
    F, G = {}, {}
    for _t in range(6):
        a1, a2 = (int(x) for x in rng.integers(0, 3, size=2))
        b1 = int(rng.integers(0, a1 + a2 + 3))
        b2 = a1 + a2 + 2 - b1
        F[(-1 - a1, -1 - a2, b1, b2)] = Fr(int(rng.integers(-5, 6)))
        G[(a1, a2, -1 - b1, -1 - b2)] = Fr(int(rng.integers(-5, 6)))                  # passendes Gegenstück (Homogenität -4)
        d1, d2 = (int(x) for x in rng.integers(0, 3, size=2))
        c1_ = int(rng.integers(0, max(1, d1 + d2 - 1)))
        if d1 + d2 - 2 - c1_ >= 0:
            G[(c1_, d1 + d2 - 2 - c1_, -1 - d1, -1 - d2)] = Fr(int(rng.integers(-5, 6)))
    res_nonzero.append(res_(lp_mul(F, G)))
    for al in range(4):
        for be in range(4):
            inv_def = max(inv_def, abs(res_(lp_mul(lp_E(al, be, F), G)) + res_(lp_mul(F, lp_E(al, be, G)))))
ZAHLEN["cech"] = {"omega - Res(Č psi . tau^* Č phi) (max, Paare)": [res_def, n_paar], "Res ohne tau^*: Paare mit Residuum": res_ohne,
                  "Invarianzdefekt Res(E F G) + Res(F E G) (exakt)": str(inv_def),
                  "Res(F G) der Stichproben": [str(x) for x in res_nonzero]}
check(f"[E/L] (e) Čech-Bild (exakt, Brüche): Für alle {n_paar} Basispaare bis Stufe {K7} ist omega(psi, phi) = Res(Č psi . "
      f"tau^* Č phi), das Residuum bei (z1 z2 z3 z4)^-1 des Produkts der a-Klasse (Pole in z1, z2) mit der tau-gespiegelten "
      f"b-Klasse (Pole in z3, z4) (Defekt {res_def:.0e}); ohne tau^* gibt es kein einziges Residuum ({res_ohne}). Die "
      "Residuenpaarung selbst ist exakt invariant: Res((E F) G) + Res(F (E G)) = 0 für alle 16 Weyl-geordneten Erzeuger "
      f"-z^al d_be - delta/2 und sechs zufällige Paare F in H^1(O(0)), G in H^1(O(-4)) (Defekt {inv_def}, Res(F G) = "
      f"{[str(x) for x in res_nonzero[:3]]}...). Das ist die holomorphe sl(4,C)-Invarianz der Serre-Paarung O(0) x O(-4) -> "
      "O(-4) = K_PT, also (Lesart, eindeutig bis auf den Faktor) der On-shell-Gestalt von int D^3Z tr(b dbar a): Die "
      "Tomita-Paarung ist diese Paarung nach der Spiegelung tau = Theta_1",
      res_def < 1e-12 and res_ohne == 0 and inv_def == 0 and all(x != 0 for x in res_nonzero))
befund("Teil 1: Die Tomita-Paarung ist die native a-b-Kinetik bis auf die Spiegelung tau. omega = <J psi|phi> ist "
       "nichtentartet mit Betrag 1 je Stufe und verbindet genau Homogenität 0 mit -4. Die richtige Realform-Aussage ist "
       "J E J = -(tau E tau)^+ für alle 16 Erzeuger mit tau = Theta_1 (z84) = Punktspiegelung x -> -x; J ist die CPT-Abbildung "
       "der Leiter. Naiv invariant ist omega nur unter Lorentz + Dilatation + u(1) (Dimension 8), und eine unverdrehte "
       "su(2,2)-invariante Paarung zwischen Sektoren positiver Energie gibt es nicht (No-go); die verdrehte gibt es genau "
       "einmal, und sie ist omega. Im Čech-Bild ist omega die Residuen-(Serre-)Paarung O(0) x O(-4) nach tau^*, die exakt "
       "sl(4)-invariant und bis auf den Faktor eindeutig ist, also (Lesart) die On-shell-Gestalt von int tr(b dbar a)", True)
laufzeit()

# =================================================================== Teil 2
say("\nTeil 2: Off-shell-Fortsetzung über die OS-Struktur (theta und tau = theta J)")
say("  (a) Kohärente Zustände: Fockreihe gegen geschlossene Formen")
KC = 26
SC = {D: sektor(D, KC) for D in (-2, 0, 2)}
IC = {D: {st: i for i, st in enumerate(B)} for D, B in SC.items()}


def paar_op(D, a, b):
    """d_a^+ e_b^+ auf dem Sektor D (die Paarerzeugung hebt die Stufe um 1)."""
    rows, cols, vals = [], [], []
    for j, st in enumerate(SC[D]):
        s = list(st)
        c = math.sqrt((s[a] + 1)*(s[2 + b] + 1))
        s[a] += 1
        s[2 + b] += 1
        if tuple(s) in IC[D]:
            rows.append(IC[D][tuple(s)])
            cols.append(j)
            vals.append(c)
    return sps.csr_matrix((vals, (rows, cols)), shape=(len(SC[D]), len(SC[D])))


PAAR = {D: [[paar_op(D, a, b) for b in range(2)] for a in range(2)] for D in SC}
SWAP = {D: np.array([IC[-D][tausch(st)] for st in SC[D]]) for D in SC}


def start(D, s=None):
    v = np.zeros(len(SC[D]), complex)
    if D == 0:
        v[IC[0][(0, 0, 0, 0)]] = 1.0
    else:
        o = (lambda a, b: (a, b, 0, 0)) if D == 2 else (lambda a, b: (0, 0, a, b))     # (s.d^+)^2 bzw. (s.e^+)^2 Omega
        v[IC[D][o(2, 0)]] = s[0]**2*math.sqrt(2)
        v[IC[D][o(1, 1)]] = 2*s[0]*s[1]
        v[IC[D][o(0, 2)]] = s[1]**2*math.sqrt(2)
    return v


def koh(D, W, s=None):
    """exp(sum W_ab d_a^+ e_b^+) angewandt auf den tiefsten Zustand (Reihe bis Stufe KC)."""
    P = sum(W[a, b]*PAAR[D][a][b] for a in range(2) for b in range(2))
    v = start(D, s)
    out = v.copy()
    for k in range(1, KC + 1):
        v = P @ v/k
        out = out + v
    return out


omega_f = lambda D, psi, phi: complex(np.sum(psi*phi[SWAP[D]]))                       # omega(psi in D, phi in -D)
klein = lambda s: s*cz(2, 2)/2
herm0 = lambda W1, W2: 1/np.linalg.det(np.eye(2) - W1.conj().T @ W2)
om0 = lambda W1, W2: 1/np.linalg.det(np.eye(2) - W1 @ W2)
herm2 = lambda W1, a, W2, b: 2*(b @ np.linalg.inv(np.eye(2) - W1.conj() @ W2.T) @ a.conj())**2/np.linalg.det(np.eye(2) - W1.conj() @ W2.T)
om_ab = lambda W1, u, W2, v: 2*(u @ np.linalg.inv(np.eye(2) - W2 @ W1) @ v)**2/np.linalg.det(np.eye(2) - W1 @ W2)
koh_def = {"Hilbert D = 0": 0.0, "Tomita D = 0": 0.0, "Hilbert D = +2": 0.0, "Tomita a-b": 0.0}
for _ in range(6):
    W1, W2, u, v = klein(0.4), klein(0.4), cz(2), cz(2)
    p0, q0 = koh(0, W1), koh(0, W2)
    koh_def["Hilbert D = 0"] = max(koh_def["Hilbert D = 0"], abs(np.vdot(p0, q0)/herm0(W1, W2) - 1))
    koh_def["Tomita D = 0"] = max(koh_def["Tomita D = 0"], abs(omega_f(0, p0, q0)/om0(W1, W2) - 1))
    p2, q2 = koh(2, W1, u), koh(2, W2, v)
    koh_def["Hilbert D = +2"] = max(koh_def["Hilbert D = +2"], abs(np.vdot(p2, q2)/herm2(W1, u, W2, v) - 1))
    pa = koh(-2, W1, u)
    koh_def["Tomita a-b"] = max(koh_def["Tomita a-b"], abs(omega_f(-2, pa, q2)/om_ab(W1, u, W2, v) - 1))
ZAHLEN["kohaerent"] = {k: float(v) for k, v in koh_def.items()} | {"Stufen der Reihe": KC}
check("[E] (a) Kohärente Zustände Psi_W = exp(W_ab d_a^+ e_b^+) Psi_0 (Fockreihe bis Stufe 26, sechs zufällige komplexe W1, W2, "
      "u, v): Hilbertkerne <Psi_W1|Psi_W2> = det(1 - W1^+ W2)^-1 (Skalar) und 2 (v^T (1 - W1-quer W2^T)^-1 u-quer)^2/det(...) "
      "(D = +2), Tomitakerne omega(Psi_W1, Psi_W2) = det(1 - W1 W2)^-1 und omega(Psi^a_(W1,u), Psi^b_(W2,v)) = "
      "2 (u^T (1 - W2 W1)^-1 v)^2/det(1 - W1 W2) (max rel. Abweichungen "
      f"{ {k: f'{v:.0e}' for k, v in koh_def.items()} }): Der Tomitakern ist holomorph in beiden Argumenten, der Hilbertkern "
      "antiholomorph im ersten", max(koh_def.values()) < 1e-12)

# --- Čech-Bild des kohärenten Zustands = Elementarzustand 1/((A.Z)(B.Z)) (z111)
W9 = klein(0.6)
p9 = koh(0, W9)
el_def = 0.0
for i, st in enumerate(SC[0]):
    if stufe(st) > 6:
        continue
    n1, n2, m1, m2 = st
    cech_koeff = p9[i]*math.sqrt(math.factorial(n1)*math.factorial(n2)/(math.factorial(m1)*math.factorial(m2)))
    taylor = sum(math.comb(n1, j)*W9[0, 0]**j*W9[0, 1]**(n1 - j)*math.comb(n2, m1 - j)*W9[1, 0]**(m1 - j)*W9[1, 1]**(n2 - m1 + j)
                 for j in range(max(0, m1 - n2), min(n1, m1) + 1))
    el_def = max(el_def, abs(cech_koeff - taylor)/max(1.0, abs(taylor)))
ZAHLEN["elementarzustand"] = el_def
check("[E] (a) Das Čech-Bild des kohärenten Zustands ist der Elementarzustand von z111: Č(Psi_W) = 1/((z1 - (W z')_1)(z2 - "
      "(W z')_2)), z' = (z3, z4), also 1/((A.Z)(B.Z)) mit zwei Ebenen durch die Gerade {z = W z'} (Koeffizienten bis Stufe 6 "
      f"gegen die Taylorreihe, Defekt {el_def:.0e}). Die Fock-Normierung (Omega <-> 1/(z1 z2)) legt die Büschelwahl fest, "
      "die in z111 offen war: <ab> = 1", el_def < 1e-12)
laufzeit()

say("  (b) OS-Fortsetzung: ein Kern, zwei Spiegelungen; Penrose-Transformation = Tomita-Paarung mit dem gespiegelten Punkt")
cay = lambda Y: (Y - 1j*np.eye(2)) @ np.linalg.inv(Y + 1j*np.eye(2))                   # Röhre Im Y > 0 -> Kugel
J_def, herm_G, quat = 0.0, 0.0, 0.0
for _ in range(6):
    W = klein(0.4)
    pv = koh(0, W)
    J_def = max(J_def, float(np.abs(np.conj(pv)[SWAP[0]] - koh(0, W.conj().T)).max()))
    Y1 = 1j*np.eye(2) + 0.35*cz(2, 2)
    Y2 = 1j*np.eye(2) + 0.35*cz(2, 2)
    hk = np.vdot(koh(0, cay(Y1)), koh(0, cay(Y2)))
    g_h = -4*hk/(np.linalg.det(Y1.conj().T - 1j*np.eye(2))*np.linalg.det(Y2 + 1j*np.eye(2)))
    herm_G = max(herm_G, abs(g_h*np.linalg.det(Y1.conj().T - Y2) - 1))
    uq, _r = np.linalg.qr(cz(2, 2))
    uq = uq/np.sqrt(np.linalg.det(uq))                                                  # u in SU(2): euklidischer Schnitt W = r u
    rq = float(rng.uniform(0.3, 0.9))
    Wq = rq*uq
    quat = max(quat, float(np.abs(np.linalg.inv(Wq).conj().T - uq/rq).max()), float(np.abs(np.linalg.inv(Wq) - uq.conj().T/rq).max()))
ZAHLEN["os_struktur"] = {"J Psi_W - Psi_(W^+) (Fock)": J_def, "Hilbertkern gegen 1/det(Y1^+ - Y2)": herm_G,
                         "theta(r u) = u/r, tau(r u) = u^+/r (Defekt)": quat}
check(f"[E] (b) OS-Struktur, ein Kern G und zwei Spiegelungen: J Psi_W = Psi_(W^+) (Fockreihe, Defekt {J_def:.0e}), also J = "
      "theta o tau mit theta: W -> W^-+ (in der Röhre Y -> Y^+, die Konjugation zwischen zukünftiger und vergangener Röhre; auf "
      "dem euklidischen Schnitt W = r u, u in SU(2), die radiale Inversion r u -> u/r, z84s reflexionspositive Spiegelung) und "
      f"tau: W -> W^-1 (Y -> -Y; euklidisch r u -> u^+/r) (Defekt {quat:.0e}). Die Hilbertform ist G an der theta-Spiegelung: "
      f"-4 <Psi_Y1|Psi_Y2>/(det(Y1^+ - i) det(Y2 + i)) = 1/det(Y1^+ - Y2) (Defekt {herm_G:.0e}); die Tomita-Paarung ist "
      "derselbe Kern an der tau-Spiegelung (unten). theta macht die positive Hilbertform (OS), tau = theta J die bilineare "
      "Kinetik; zurückgespiegelt bleibt in beiden Fällen der holomorphe Propagator G",
      J_def < 1e-13 and herm_G < 1e-10 and quat < 1e-12)


def penrose(W, M, Mpkt=64):
    """(1/2 pi i) Kontur um die erste Nullstelle von f_W auf L_M = {z = M z'}, z' = pi = (1, s) (Trapezregel, wie z111)."""
    R = M - W
    s1, s2 = -R[0, 0]/R[0, 1], -R[1, 0]/R[1, 1]
    rad = 0.4*abs(s1 - s2)
    th = 2*np.pi*np.arange(Mpkt)/Mpkt
    s = s1 + rad*np.exp(1j*th)
    f = 1/((R[0, 0] + R[0, 1]*s)*(R[1, 0] + R[1, 1]*s))
    return np.mean(f*rad*np.exp(1j*th))


ratio_os, ratio_ohne = [], []
for _ in range(12):
    W = klein(0.4)
    Minv = klein(0.4)                                                                   # gespiegelter Punkt tau.M = M^-1 in der Kugel
    M = np.linalg.inv(Minv)
    phi = penrose(W, M)
    om = omega_f(0, koh(0, Minv), koh(0, W))
    ratio_os.append(phi*np.linalg.det(M)/om)
    phiY = penrose(W, Minv)                                                             # ungespiegelte Lesart: Punkt = M^-1
    ratio_ohne.append(phiY*np.linalg.det(Minv)/om)
spread = lambda v: float(np.abs(np.array(v) - v[0]).max()/abs(v[0]))
ZAHLEN["penrose_gegen_tomita"] = {"Verhältnis (gespiegelt)": str(np.round(ratio_os[0], 12)), "Streuung gespiegelt": spread(ratio_os),
                                  "Streuung ungespiegelt": spread(ratio_ohne)}
check("[E] (b) Penrose-Transformation des Elementarzustands Č(Psi_W) auf der Geraden L_M = {z = M z'} (Konturintegral wie "
      "z111, 64 Punkte) mal det M ist für zwölf zufällige Paare (W, M) gleich der Tomita-Paarung omega(Psi_(tau.M), Psi_W) "
      f"mit dem gespiegelten Punktzustand tau.M = M^-1 aus der Fockreihe (Verhältnis {complex(ratio_os[0]):.6f}, Streuung "
      f"{spread(ratio_os):.0e}): Die Auswertung am Punkt M ist die Tomita-Paarung mit dem tau-gespiegelten kohärenten "
      f"Zustand. Ungespiegelte Lesart (Punkt = M^-1) streut um {spread(ratio_ohne):.2f}. Mit z111: Das ist der Cauchykern "
      "-1/(<ab> det(M - W)) mit <ab> = 1", spread(ratio_os) < 1e-10 and abs(ratio_os[0] + 1) < 1e-10 and spread(ratio_ohne) > 0.1)

def os_kern(X1, X2, D=0, u=None, v=None):
    """Tomita-Paarung der kohärenten Zustände an Y1 = -X1 (= tau.X1) und Y2 = X2, Gewichte det(i - X1), det(X2 + i) abgezogen."""
    W1, W2 = cay(-X1), cay(X2)
    w = np.linalg.det(1j*np.eye(2) - X1)*np.linalg.det(X2 + 1j*np.eye(2))
    if D == 0:
        return -4*omega_f(0, koh(0, W1), koh(0, W2))/w
    return 8*omega_f(-2, koh(-2, W1, u), koh(2, W2, v))/w


tube_def, ab_def, trans_tw, trans_un = 0.0, 0.0, 0.0, 0.0
for _ in range(8):
    X1 = -1j*np.eye(2) + 0.35*cz(2, 2)                                                  # untere Röhre (gespiegelt)
    X2 = 1j*np.eye(2) + 0.35*cz(2, 2)                                                   # obere Röhre
    g = os_kern(X1, X2)
    tube_def = max(tube_def, abs(g*np.linalg.det(X1 - X2) - 1))
    u, v = cz(2), cz(2)
    ut, vt = (1j*np.eye(2) - X1).T @ u, (X2 + 1j*np.eye(2)) @ v
    gab = os_kern(X1, X2, 2, u, v)
    soll = (ut @ np.linalg.inv(X2 - X1) @ vt)**2/np.linalg.det(X1 - X2)
    ab_def = max(ab_def, abs(gab/soll - 1))
    a_ = 0.15*cz(2, 2)
    trans_tw = max(trans_tw, abs(os_kern(X1 + a_, X2 + a_)/g - 1))
    Y1, Y2 = -X1, X2                                                                    # ungespiegelte Lesart: Orte Y1, Y2
    g_un = -4*omega_f(0, koh(0, cay(Y1)), koh(0, cay(Y2)))/(np.linalg.det(Y1 + 1j*np.eye(2))*np.linalg.det(Y2 + 1j*np.eye(2)))
    g_un_a = (-4*omega_f(0, koh(0, cay(Y1 + a_)), koh(0, cay(Y2 + a_)))
              / (np.linalg.det(Y1 + a_ + 1j*np.eye(2))*np.linalg.det(Y2 + a_ + 1j*np.eye(2))))
    trans_un = max(trans_un, abs(g_un_a/g_un - 1))
ZAHLEN["roehre"] = {"g det(X1 - X2) - 1": tube_def, "a-b gegen (u~ (X2 - X1)^-1 v~)^2/det": ab_def,
                    "Translation gespiegelt": trans_tw, "Translation ungespiegelt": trans_un}
check("[E] (b) Röhrenform (Fockreihe, acht Paare X1 in der unteren, X2 in der oberen Röhre): Nach Abzug der konformen Gewichte "
      f"ist der OS-Kern G(X1, X2) = omega(Psi_(-X1), Psi_X2) genau 1/det(X1 - X2) (Defekt {tube_def:.0e}) und im a-b-Sektor "
      f"(u~ (X2 - X1)^-1 v~)^2/det(X1 - X2) mit u~ = (i - X1)^T u, v~ = (X2 + i) v (Defekt {ab_def:.0e}); er hängt nur von "
      f"X1 - X2 ab (gemeinsame komplexe Translation: {trans_tw:.0e}). Ungespiegelt (Orte Y1 = -X1, Y2) ist der Kern "
      f"1/det(Y1 + Y2) und nicht translationsinvariant ({trans_un:.2f}): Erst die Rückspiegelung mit tau (= theta J) macht aus "
      "der Tomita-Paarung den Propagator, holomorph in beiden Punkten, zwischen vergangener und zukünftiger Röhre",
      tube_def < 1e-10 and ab_def < 1e-10 and trans_tw < 1e-10 and trans_un > 0.05)
laufzeit()

say("  (c) CSW-Eichung: dbar^-1 zwischen zwei Geraden (delta-quer-Lokalisierung) gegen den a-b-Kern")


def lokalisiere(x, y, Zst):
    """Z* + s Z(lambda_1) + t Z(lambda_2) = 0 mit Z(lambda) = (lambda, x lambda) auf L_x und (lambda, y lambda) auf L_y:
    lineares 4x4-System für (s lambda_1, t lambda_2); dazu die holomorphe Jacobi-Determinante nach (sigma1, sigma2, s, t)."""
    Mat = np.block([[np.eye(2), np.eye(2)], [x, y]])
    sol = np.linalg.solve(Mat, -Zst)
    sl1, tl2 = sol[:2], sol[2:]
    s, t = sl1[0], tl2[0]                                                               # lambda = (1, sigma)
    l1, l2 = sl1/s, tl2/t
    e2 = np.array([0, 1])
    J = np.column_stack([s*np.concatenate([e2, x @ e2]), t*np.concatenate([e2, y @ e2]),
                         np.concatenate([l1, x @ l1]), np.concatenate([l2, y @ l2])])
    return s, t, l1, l2, np.linalg.det(J), -sl1


lok = {"parallel (I_pi)": 0.0, "lambda* - (x - y)^-1 eta (I_pi)": 0.0, "Jacobi/(s t det(x - y)) Streuung": [],
       "parallel (Z* außerhalb I_pi)": 1e9, "Translation lambda* (außerhalb)": 0.0, "a-b-Kern gegen (u~.lambda*)^2/det": [],
       "a-b-Kern mit eta != v~": []}
ang = lambda a, b: a[0]*b[1] - a[1]*b[0]
for _ in range(10):
    X1 = -1j*np.eye(2) + 0.35*cz(2, 2)
    X2 = 1j*np.eye(2) + 0.35*cz(2, 2)
    x, y = X2, X1
    u, v = cz(2), cz(2)
    ut, vt = (1j*np.eye(2) - X1).T @ u, (X2 + 1j*np.eye(2)) @ v
    s, t, l1, l2, Jd, lst = lokalisiere(x, y, np.concatenate([np.zeros(2), vt]))       # Z* = (0, eta) auf I_pi, eta = v~
    lok["parallel (I_pi)"] = max(lok["parallel (I_pi)"], abs(ang(l1, l2))/(np.linalg.norm(l1)*np.linalg.norm(l2)))
    lok["lambda* - (x - y)^-1 eta (I_pi)"] = max(lok["lambda* - (x - y)^-1 eta (I_pi)"],
                                                 float(np.abs(lst - np.linalg.inv(x - y) @ vt).max()))
    lok["Jacobi/(s t det(x - y)) Streuung"].append(Jd/(s*t*np.linalg.det(x - y)))
    gab = os_kern(X1, X2, 2, u, v)
    lok["a-b-Kern gegen (u~.lambda*)^2/det"].append(gab*np.linalg.det(X1 - X2)/(ut @ lst)**2)
    _, _, _, _, _, lst2 = lokalisiere(x, y, np.concatenate([np.zeros(2), cz(2)]))
    lok["a-b-Kern mit eta != v~"].append(gab*np.linalg.det(X1 - X2)/(ut @ lst2)**2)
    Zst_aus = np.concatenate([cz(2), vt])                                               # pi* != 0: nicht auf I_pi
    _, _, m1, m2, _, la = lokalisiere(x, y, Zst_aus)
    lok["parallel (Z* außerhalb I_pi)"] = min(lok["parallel (Z* außerhalb I_pi)"], abs(ang(m1, m2))/(np.linalg.norm(m1)*np.linalg.norm(m2)))
    a_ = 0.3*cz(2, 2)
    _, _, _, _, _, la_t = lokalisiere(x + a_, y + a_, Zst_aus)
    lok["Translation lambda* (außerhalb)"] = max(lok["Translation lambda* (außerhalb)"], float(np.abs(la_t - la).max()/np.abs(la).max()))
jac_sp = spread(lok["Jacobi/(s t det(x - y)) Streuung"])
ab_sp, ab_sp2 = spread(lok["a-b-Kern gegen (u~.lambda*)^2/det"]), spread(lok["a-b-Kern mit eta != v~"])
ZAHLEN["csw_lokalisierung"] = {"parallel (I_pi)": lok["parallel (I_pi)"], "lambda* - (x - y)^-1 eta": lok["lambda* - (x - y)^-1 eta (I_pi)"],
                               "Jacobi/(s t det(x - y))": str(np.round(lok["Jacobi/(s t det(x - y)) Streuung"][0], 12)),
                               "Jacobi-Streuung": jac_sp, "a-b-Kern/((u~.lambda*)^2/det)": str(np.round(lok["a-b-Kern gegen (u~.lambda*)^2/det"][0], 12)),
                               "Streuung (eta = v~)": ab_sp, "Streuung (eta zufällig)": ab_sp2,
                               "min |<l1 l2>| außerhalb I_pi": lok["parallel (Z* außerhalb I_pi)"],
                               "Translationsdefekt außerhalb I_pi": lok["Translation lambda* (außerhalb)"]}
check("[E] (c) CSW-Eichung mit Z* = (0, eta) auf I_pi (z111): Der Träger von dbar^-1 = int ds dt/(s t) t^4 delta-quer^4(Z* + "
      "s Z1 + t Z2) zwischen L_x und L_y ist das eine Punktpaar mit lambda_1 || lambda_2 (Defekt "
      f"{lok['parallel (I_pi)']:.0e}) und s lambda_1 = -(x - y)^-1 eta (Defekt {lok['lambda* - (x - y)^-1 eta (I_pi)']:.0e}); "
      f"die holomorphe Jacobi-Determinante ist s t det(x - y) (Verhältnis {complex(lok['Jacobi/(s t det(x - y)) Streuung'][0]):.3f}, "
      f"Streuung {jac_sp:.0e}), mit t = -s also das Gewicht 1/det(x - y), der Cauchy-Penrose-Kern von z111. Der OS-Tomita-"
      f"Kern des a-b-Sektors ist exakt (u~ . lambda*)^2/det(x - y) mit diesem lambda* bei eta = v~ (Streuung {ab_sp:.0e}); mit "
      f"einer anderen Referenz streut das Verhältnis um {ab_sp2:.2f}. Kontrolle Z* außerhalb I_pi: lambda_1, lambda_2 nicht "
      f"parallel (min |<l1 l2>| {lok['parallel (Z* außerhalb I_pi)']:.2f}), lambda* nicht translationsinvariant "
      f"({lok['Translation lambda* (außerhalb)']:.2f}). Der a-b-Kern hat also genau die Gestalt von dbar^-1 in CSW-Eichung "
      "(Gewicht 1/det(x - y), Auswertung bei lambda* = (x - y)^-1 eta), mit der b-Polarisation v~ als Referenz",
      lok["parallel (I_pi)"] < 1e-12 and lok["lambda* - (x - y)^-1 eta (I_pi)"] < 1e-12 and jac_sp < 1e-10
      and abs(abs(lok["Jacobi/(s t det(x - y)) Streuung"][0]) - 1) < 1e-10 and ab_sp < 1e-9 and ab_sp2 > 0.1
      and lok["parallel (Z* außerhalb I_pi)"] > 1e-3 and lok["Translation lambda* (außerhalb)"] > 0.05)
laufzeit()

say("  (d) Impulsraum: 4 pi^2/P^2 und lambda_P = P|v]")
xs = sp.Matrix(2, 2, sp.symbols("x00 x01 x10 x11"))
us, vs_ = sp.symbols("u0 u1"), sp.symbols("v0 v1")
Dop = lambda f: sum(us[A]*vs_[B]*sp.diff(f, xs[A, B]) for A in range(2) for B in range(2))
detx = xs.det()
lhs = Dop(Dop(1/detx))
rhs = 2*(sp.Matrix([vs_]) @ xs.inv() @ sp.Matrix(us))[0, 0]**2/detx
abl_ok = sp.simplify(sp.together(lhs - rhs)) == 0
mp.mp.dps = 20
hank = {}
for Pm in (0.5, 1.0, 2.0):
    integral = mp.quadosc(lambda w: mp.besselj(1, w), [0, mp.inf], zeros=lambda k: mp.besseljzero(1, k))
    hank[Pm] = float(4*mp.pi**2/Pm*integral/Pm*Pm**2)
res_fak = 0.0
for _ in range(10):
    lam_, lt_ = cz(2), cz(2)
    Pon = np.outer(lam_, lt_)
    uu, vv = cz(2), cz(2)
    res_fak = max(res_fak, abs((uu @ Pon @ vv)**2 - (uu @ lam_)**2*(lt_ @ vv)**2))
ZAHLEN["impulsraum"] = {"D^2(1/det x) = 2 (v^T x^-1 u)^2/det x (sympy)": bool(abl_ok),
                        "P^2 FT[1/|x|^2] für |P| = 0.5, 1, 2": [r(h_, 10) for h_ in hank.values()], "4 pi^2": r(4*np.pi**2, 10),
                        "On-shell-Faktorisierung (u P v)^2 = <u l>^2 [l~ v]^2": res_fak}
check("[E] (d) Impulsraum: Mit D = u_A v_B d/dx_(AB) gilt exakt (u x^-1 v)^2/det x = D^2(1/det x)/2 (sympy), und die "
      "euklidische Fouriertransformation von 1/|x|^2 ist 4 pi^2/|P|^2 (Hankel, z111: P^2 FT = "
      f"{[f'{h_:.10f}' for h_ in hank.values()]}, 4 pi^2 = {4*np.pi**2:.10f}). Also ist der a-b-Kern im Impulsraum "
      "-2 pi^2 (u.P.v)^2/P^2 = -2 pi^2 <u lambda_P>^2/P^2 mit lambda_P = P|v]: genau der CSW-Propagator 1/P^2 mit der "
      f"Off-shell-Fortsetzung P|eta] und eta = v; auf der Schale P = lambda lambda~ zerfällt der Zähler in <u lambda>^2 "
      f"[lambda~ v]^2, das Produkt der beiden Helizitäts-Wellenfunktionen (Defekt {res_fak:.0e})",
      abl_ok and all(abs(h_ - 4*np.pi**2) < 1e-8 for h_ in hank.values()) and res_fak < 1e-10)
laufzeit()

# ------------------------------------------------------------------ Vergleichsmaschine (wie z111)
EPS2 = np.array([[0, 1], [-1, 0]])


def kinematik(n):
    lam, lt = cz(n, 2), cz(n, 2)
    S = sum(np.outer(lam[i], lt[i]) for i in range(n - 2))
    sol = np.linalg.solve(np.array([lam[n - 2], lam[n - 1]]).T, -S)
    lt[n - 2], lt[n - 1] = sol[0], sol[1]
    return lam, lt


def impuls(lam, lt, idx):
    return sum(np.outer(lam[i], lt[i]) for i in idx)


def pt(lam, h):
    n = len(h)
    m = [i for i in range(n) if h[i] < 0]
    return ang(lam[m[0]], lam[m[1]])**4/np.prod([ang(lam[i], lam[(i + 1) % n]) for i in range(n)])


def ptq(lt, h):
    n = len(h)
    p = [i for i in range(n) if h[i] > 0]
    return ang(lt[p[0]], lt[p[1]])**4/np.prod([ang(lt[i], lt[(i + 1) % n]) for i in range(n)])


def a3(lam, lt, h):
    nm = sum(1 for x in h if x < 0)
    sa = max(abs(ang(lam[i], lam[(i + 1) % 3])) for i in range(3))
    ss = max(abs(ang(lt[i], lt[(i + 1) % 3])) for i in range(3))
    if nm == 2:
        return 0.0 if sa < 1e-9*ss else pt(lam, h)
    if nm == 1:
        return 0.0 if ss < 1e-9*sa else ptq(lt, h)
    return 0.0


def nullzerlegung(Ph):
    k = int(np.argmax(np.abs(Ph).max(0)))
    zeile = int(np.argmax(np.abs(Ph[:, k])))
    lP = Ph[:, k].copy()
    return lP, Ph[zeile]/lP[zeile]


def amp(lam, lt, h):
    """Baumamplitude aus BCFW (Rekursion bis zu den Dreipunktamplituden, Verschiebung [-, +>), Kopplung 1."""
    n = len(h)
    if n == 3:
        return a3(lam, lt, h)
    nm = sum(1 for x in h if x < 0)
    if nm < 2 or n - nm < 2:
        return 0.0
    for i in range(n):
        j = (i + 1) % n
        if (h[i], h[j]) == (-1, 1):
            return bcfw(lam, lt, h, i, j)
    return None


def bcfw(lam, lt, h, i, j):
    n = len(h)
    folge = [(j + k) % n for k in range(n)]
    summe = 0.0
    for m in range(2, n - 1):
        L, R = folge[:m], folge[m:]
        P = impuls(lam, lt, L)
        adj = np.array([[P[1, 1], -P[0, 1]], [-P[1, 0], P[0, 0]]])
        z = np.linalg.det(P)/(lt[j] @ adj @ lam[i])
        lam2, lt2 = lam.copy(), lt.copy()
        lt2[i] = lt[i] + z*lt[j]
        lam2[j] = lam[j] - z*lam[i]
        lP, ltP = nullzerlegung(impuls(lam2, lt2, L))
        for hp in (1, -1):
            aL = amp(np.array([lam2[k] for k in L] + [lP]), np.array([lt2[k] for k in L] + [-ltP]), [h[k] for k in L] + [hp])
            if aL == 0:
                continue
            aR = amp(np.array([lP] + [lam2[k] for k in R]), np.array([ltP] + [lt2[k] for k in R]), [-hp] + [h[k] for k in R])
            summe += aL/np.linalg.det(P)*aR
    return summe


def bogen(n):
    for s in range(n - 1):
        for l in range(2, n - 1):
            if s + l <= n - 1:
                L = list(range(s, s + l))
                yield L, list(range(s + l, n)) + list(range(0, s))


def vertex(spin, hel, cV, g):
    """Nahtvertex c_V g^(n-2) <ab>^4/(<12>...<n1>): Nenner nativ (z110, Hypothese T), Zähler unter P_N4; 0 außer bei genau
    zwei negativen Helizitäten."""
    m = [i for i, x in enumerate(hel) if x < 0]
    if len(m) != 2:
        return 0.0
    n = len(hel)
    return cV*g**(n - 2)*ang(spin[m[0]], spin[m[1]])**4/np.prod([ang(spin[i], spin[(i + 1) % n]) for i in range(n)])


def zusammenbau(lam, lt, h, eta, N=1.0, cV=1.0, g=1.0, leitung="tomita"):
    """S = int tr(b dbar a) + Nahtvertizes in CSW-Eichung: Einzelvertex plus Summe über die Bogenkanäle L|R mit einem
    Austausch N/P^2, lambda_P = P eta. leitung: 'tomita' (omega: h an einem Ende, -h am anderen), 'hilbert' (<psi|phi>:
    gleiche Helizität an beiden Enden), 'ungespiegelt' (omega ohne tau: Kern G(x + y), Impulsträger P_L = P_R)."""
    n = len(h)
    summe = vertex(lam, h, cV, g)
    groesste, ndiag = abs(summe), 0
    for L, R in bogen(n):
        P = impuls(lam, lt, L)
        if leitung == "ungespiegelt" and not np.allclose(P, impuls(lam, lt, R)):
            continue
        e = eta(L, R) if callable(eta) else eta
        lP = P @ e
        for hL in (1, -1):
            hR = hL if leitung == "hilbert" else -hL
            VL = vertex([lam[k] for k in L] + [lP], [h[k] for k in L] + [hL], cV, g)
            VR = vertex([lP] + [lam[k] for k in R], [hR] + [h[k] for k in R], cV, g)
            beitrag = VL*N/np.linalg.det(P)*VR
            if beitrag != 0:
                ndiag += 1
                groesste = max(groesste, abs(beitrag))
            summe += beitrag
    return summe, ndiag, groesste


I_PI = np.outer(U_PI[0], U_PI[1]) - np.outer(U_PI[1], U_PI[0])
om_von = lambda Z: KAPPA*np.array([Z[0] - Z[2], Z[1] - Z[3]])/np.sqrt(2)


def eta_ipi():
    """Referenz aus I_pi wie z111: Z* = I_pi W, eta = eps omega(Z*)."""
    return EPS2 @ om_von(I_PI @ cz(4))


LAGEN = {n: [h for h in itertools.product((-1, 1), repeat=n) if sum(x < 0 for x in h) in (2, 3)] for n in (4, 5)}
KIN = {n: [kinematik(n) for _ in range(3)] for n in (4, 5)}
REF = {(n, h, k): amp(*KIN[n][k], list(h)) for n in (4, 5) for h in LAGEN[n] for k in range(3)}
NMHV5 = [h for h in LAGEN[5] if sum(x < 0 for x in h) == 3]

say("  (e) Normierung durch Faktorisierung")
k_sym = sp.Symbol("k", positive=True)
zerl_ok = all(sp.simplify(-k_sym*(-1)**n*k_sym**(-sp.Rational(n, 2)) - (-1)*(-1/sp.sqrt(k_sym))**(n - 2)) == 0 for n in range(3, 9))
fak = {}
for k_ in (1, 2, 3):
    g_ = -1/math.sqrt(k_)
    for cV in (1.0, -1.0):
        for N in (1.0, -1.0, 2.0):
            rat = []
            for h in NMHV5:
                for kk in range(3):
                    lam, lt = KIN[5][kk]
                    A = zusammenbau(lam, lt, list(h), eta_ipi(), N, cV, g_)[0]
                    rat.append(A/(cV*g_**3*REF[(5, h, kk)]))
            fak[(k_, cV, N)] = (complex(np.mean(rat)), float(np.abs(np.array(rat) - rat[0]).max()))
fak_def = max(abs(v[0] - k[1]*k[2]) + v[1] for k, v in fak.items())
reph = 0.0
for h in LAGEN[5]:
    for kk in range(3):
        lam, lt = KIN[5][kk]
        nm = sum(x < 0 for x in h)
        A = zusammenbau(lam, lt, list(h), eta_ipi(), -1.0, 1.0, 1.0)[0]
        reph = max(reph, abs(A - (-1)**(nm - 2)*REF[(5, h, kk)])/abs(REF[(5, h, kk)]))
ZAHLEN["faktorisierung"] = {"z110: -k (-1)^n k^(-n/2) = -(-k^(-1/2))^(n-2), n = 3..8": bool(zerl_ok),
                            "A_CSW/(c_V g^(n-2) A_BCFW) = c_V N (max Defekt, k = 1, 2, 3; c_V, N = +-1, 2)": fak_def,
                            "Umphasung: N = -1 <=> (-1)^(n_- - 2) (max rel.)": reph,
                            "Beispiele (k, c_V, N) -> Verhältnis": {str(k): str(np.round(v[0], 10)) for k, v in list(fak.items())[:6]}}
check(f"[E] (e) Faktorisierung legt die Normierung fest und nur sie: Für alle zehn 5-Punkt-NMHV-Lagen (je drei Kinematikpunkte, "
      f"eta auf I_pi) ist A_CSW(N, c_V, g)/(c_V g^(n-2) A_BCFW) = c_V N exakt (Defekt {fak_def:.0e} über k = 1, 2, 3, c_V = +-1, "
      "N = +-1, 2), unabhängig von g: Das Niveau steckt ganz in g = -k^(-1/2) und kürzt sich aus N heraus. Der Nahtkoeffizient "
      "von z110 ist von dieser Form (sympy, n = 3..8): -k (-1)^n k^(-n/2) = c_V g^(n-2) mit c_V = -1, g = -k^(-1/2). "
      "Faktorisierung verlangt c_V N = 1. Das Vorzeichen ist Konvention: N = -1 ist dasselbe wie die Umphasung der negativen "
      f"Helizität, A -> (-1)^(n_- - 2) A (alle 20 Lagen, Defekt {reph:.0e}); physikalisch ist |c_V N| = 1",
      zerl_ok and fak_def < 1e-8 and reph < 1e-9)

# --- native Kandidaten
JOmega = {th: abs(np.exp(1j*th) - 1) for th in (0.0, 0.5, np.pi)}
uhr = {}
for D in (-2, 0, 2):
    vals = set()
    for st in S1[D]:
        E_ = (st[0] + st[1] + st[2] + st[3] + 2)/2
        J3R = (st[2] - st[3])/2
        vals.add(complex(np.round(np.exp(1j*np.pi*(E_ + 2*J3R)), 12)) + 0.0)
    uhr[D] = sorted({round(x.real, 12) + 0.0 for x in vals})
BETA = 4*np.pi
pkte = [(0.3, 0.5), (1.0, 2.0), (2.5, 1.0), (4.0, 0.2)]
kms_bild = [sum(1/(2*(math.cosh(t + k*BETA) - math.cos(gm))) for k in range(-8, 9))*2*(math.cosh(t) - math.cos(gm)) - 1
            for t, gm in pkte]
kms_sp = max(kms_bild)/min(kms_bild)
KANDIDATEN = {"Tomita (J antiunitär, J Omega = Omega)": 1.0, "KMS modular Delta^(1/2) = q^D auf b": Q_NAT**2,
              "KMS modular Delta^(-1/2)": Q_NAT**-2, "KMS der Naht 1 - c(1/4)": 1 - C14,
              "Niveau k = 1 im Propagator": 1.0, "Niveau k = 2 im Propagator": 2.0, "z111 4 pi^2": 4*np.pi**2}
erfuellt = {k: abs(abs(v) - 1) < 1e-12 for k, v in KANDIDATEN.items()}
ZAHLEN["native_kandidaten"] = {"|N|": {k: r(v, 12) for k, v in KANDIDATEN.items()}, "|c_V N| = 1 erfüllt": erfuellt,
                               "Uhr Z_SU = exp(i pi (E + 2 J3_R)) je Sektor D": uhr, "|e^(i theta) - 1| (J_theta Omega = Omega?)": JOmega,
                               "Ortsraum: kohärenter Kern -1/4, kanonisch 1/(4 pi^2), Verhältnis": r(np.pi**2, 10),
                               "thermische Bildsumme bei beta = 4 pi / Vakuum - 1 (z84, vier Punkte)": [float(f"{x:.4e}") for x in kms_bild]}
check("[E/C] (e) Native Kandidaten für |N| bei |c_V| = 1 (z110 unter P_norm): Die Tomita-Paarung gibt |N| = 1, weil J "
      "antiunitär ist (alle Singulärwerte 1, Teil 1b); ihre Phase legt J Omega = Omega fest (e^(i theta) J verletzt das für "
      f"theta != 0: {[r(v, 3) for v in JOmega.values()]}), betrifft aber nur das konventionelle Vorzeichen. Erfüllt: "
      f"{[k for k, v in erfuellt.items() if v]}. Die KMS-Gewichte scheitern: modular q^(+-2) = e^(-+4 pi) = "
      f"{Q_NAT**2:.2e} bzw. {Q_NAT**-2:.2e}, an der Naht 1 - c(1/4) = {1 - C14:.4f}, und die thermische Bildsumme auf "
      f"S^1_(4 pi) x S^3 (z84) weicht ortsabhängig vom Vakuumkern ab ({[f'{x:.1e}' for x in kms_bild]}), ist also gar kein "
      "konstanter Faktor. Die Uhr ist auf a und b gleich 1 "
      f"(D = -2: {uhr[-2]}, D = +2: {uhr[2]}, Skalar D = 0: {uhr[0]}) und gibt weder Betrag noch Vorzeichen. Das Niveau "
      "kürzt sich (k = 1 erfüllt nur zufällig, k = 2 nicht). 4 pi^2 aus z111 ist das Fouriermaß, keine Zustandsnormierung: "
      "Im Ortsraum hat der kohärente Kern den Faktor -1/4, das kanonisch normierte Feld 1/(4 pi^2), ihr Verhältnis pi^2 ist "
      "die Feldnormierung des Punktzustands",
      erfuellt["Tomita (J antiunitär, J Omega = Omega)"] and not erfuellt["KMS modular Delta^(1/2) = q^D auf b"]
      and not erfuellt["KMS modular Delta^(-1/2)"] and not erfuellt["KMS der Naht 1 - c(1/4)"]
      and not erfuellt["Niveau k = 2 im Propagator"] and not erfuellt["z111 4 pi^2"]
      and uhr[-2] == uhr[2] == [1.0] and uhr[0] == [-1.0] and JOmega[0.0] == 0 and kms_sp > 1e-3 and min(JOmega[0.5], JOmega[np.pi]) > 0.1)
befund("Teil 2: Die OS-Struktur macht aus der Tomita-Paarung den Propagator: J = theta o tau, die Hilbertform ist der Kern G "
       "an der OS-Spiegelung theta, die Tomita-Paarung derselbe Kern an tau. Die Penrose-Transformation ist die Tomita-"
       "Paarung mit dem tau-gespiegelten Punktzustand; in der Röhre ist der Kern 1/det(X1 - X2) (Skalar) bzw. (u~ (X2 - X1)^-1 "
       "v~)^2/det(X1 - X2) (a-b), translationsinvariant; das ist z111s Cauchykern mit <ab> = 1 und, für a-b, genau die Gestalt "
       "von dbar^-1 in CSW-Eichung mit eta = v~ (delta-quer-Lokalisierung, Gewicht 1/det(x - y)); im Impulsraum <u lambda_P>^2/P^2 mit "
       "lambda_P = P|v]. Die Normierung ist kein freier Parameter: Faktorisierung verlangt |c_V N| = 1, und |N| = 1 kommt aus "
       "der Antiunitarität von J (unter P_norm, |c_V| = 1). Niveau (kürzt sich), KMS (e^(-+4 pi), 0,959), Uhr (trivial auf "
       "a und b) und Tomita-Phase (nur Vorzeichenkonvention) fixieren den Betrag nicht", True)
laufzeit()

# =================================================================== Teil 3
say("\nTeil 3: Zusammenbau S = int tr(b dbar a) + Nahtvertizes und Test gegen BCFW")
say("  (a) native Paarung")
G_NAT = -1.0                                                                            # g = -k^(-1/2), k = 1 (z110)
nat = {}
for n in (4, 5):
    for h in LAGEN[n]:
        nm = sum(x < 0 for x in h)
        for kk in range(3):
            lam, lt = KIN[n][kk]
            werte = []
            for _ in range(3):
                A, nd, gr = zusammenbau(lam, lt, list(h), eta_ipi(), 1.0, 1.0, G_NAT)
                werte.append(A)
            ref = G_NAT**(n - 2)*REF[(n, h, kk)]
            key = (n, nm)
            if abs(ref) > 0:
                dev = max(abs(w - ref)/abs(ref) for w in werte)
            else:
                dev = max(abs(w) for w in werte)/gr
            eta_sp = float(np.abs(np.array(werte) - werte[0]).max()/max(abs(werte[0]), gr))
            alt = nat.get(key, (0.0, 0.0, set(), 0))
            nat[key] = (max(alt[0], dev), max(alt[1], eta_sp), alt[2] | {nd}, alt[3] + 1)
anti5 = max(abs(G_NAT**3*ptq(KIN[5][kk][1], list(h)) - G_NAT**3*REF[(5, h, kk)])/abs(REF[(5, h, kk)]) for h in NMHV5 for kk in range(3))
ZAHLEN["zusammenbau_nativ"] = {f"n = {k[0]}, negative = {k[1]}": {"max rel. Abweichung": r(v[0], 15), "eta-Streuung": r(v[1], 15),
                                                                   "Austauschdiagramme": sorted(v[2]), "Fälle": v[3]}
                               for k, v in nat.items()} | {"5-Punkt-NMHV gegen Anti-MHV-Formel": r(anti5, 15)}
check("[E/C] (a) Zusammenbau mit der nativen Paarung (Vertizes c_V g^(n-2) PT unter T und P_N4 mit g = -1 aus k = 1, "
      "Austausch = OS-Tomita-Kern: a <-> b, N = 1, lambda_P = P eta mit eta aus I_pi): Alle 6 MHV- und 4 Dreiminus-Lagen bei "
      "vier Punkten und alle 10 MHV- und 10 NMHV-Lagen bei fünf Punkten (je drei Kinematikpunkte, je drei eta) stimmen mit BCFW "
      f"(Kopplung g^(n-2)) überein: max rel. {max(v[0] for k, v in nat.items() if k != (4, 3)):.0e}; die Dreiminus-Amplituden "
      f"bei vier Punkten verschwinden durch Aufhebung zweier Austauschdiagramme (Rest {nat[(4, 3)][0]:.0e} der Diagrammgröße), "
      f"MHV hat keine Austauschdiagramme ({sorted(nat[(5, 2)][2])}), 5-Punkt-NMHV {sorted(nat[(5, 3)][2])} je Lage und ist "
      f"die Anti-MHV-Formel [ab]^4/([12]...[51]) ({anti5:.0e}); die Summe hängt nicht von eta ab "
      f"({max(v[1] for v in nat.values()):.0e}). Vorzeichenkonvention: Die Orientierung des P_N4-Maßes nimmt das "
      "Schleifenvorzeichen c_V = -1 aus z110 auf (Teil 2e)",
      max(v[0] for v in nat.values()) < 1e-9 and max(v[1] for v in nat.values()) < 1e-9 and anti5 < 1e-10
      and nat[(5, 2)][2] == {0} and nat[(4, 2)][2] == {0} and min(nat[(5, 3)][2]) >= 2 and min(nat[(4, 3)][2]) >= 2)

say("  (b) Negativkontrollen: falsche Paarungen")
kontr = {}
for leit in ("hilbert", "ungespiegelt"):
    for n in (4, 5):
        for nm_ in (2, 3):
            devs, nds = [], set()
            for h in LAGEN[n]:
                if sum(x < 0 for x in h) != nm_:
                    continue
                for kk in range(3):
                    lam, lt = KIN[n][kk]
                    A, nd, gr = zusammenbau(lam, lt, list(h), eta_ipi(), 1.0, 1.0, G_NAT, leitung=leit)
                    ref = G_NAT**(n - 2)*REF[(n, h, kk)]
                    devs.append(abs(A - ref)/abs(ref) if abs(ref) > 0 else abs(A)/max(gr, 1e-300))
                    nds.add(nd)
            kontr[(leit, n, nm_)] = (max(devs), min(devs), sorted(nds))
ZAHLEN["kontrollen_paarung"] = {f"{k[0]}, n = {k[1]}, negative = {k[2]}": {"max": r(v[0], 6), "min": r(v[1], 6), "Diagramme": v[2]}
                                for k, v in kontr.items()}
h_mhv = kontr[("hilbert", 5, 2)]
h_nmhv = kontr[("hilbert", 5, 3)]
u_nmhv = kontr[("ungespiegelt", 5, 3)]
check("[E negativ] (b) Falsche Paarungen scheitern. Paarung ohne J (die Hilbertform bzw. psi^T phi lebt nur innerhalb eines "
      "Sektors, a mit a und b mit b, also gleiche Helizitätsmarke an beiden Enden): MHV bekommt unechte Austauschdiagramme "
      f"({h_mhv[2]} je Lage; 4 Punkte {kontr[('hilbert', 4, 2)][2]}) und verfehlt Parke-Taylor um {h_mhv[1]:.3g} bis "
      f"{h_mhv[0]:.3g}, NMHV verliert alle "
      f"Diagramme ({h_nmhv[2]}) und gibt 0 (Abweichung {h_nmhv[0]:.2f}). Ungespiegelte Tomita-Paarung (ohne tau: Kern "
      f"G(x + y), Impulsträger P_L = P_R gegen P_L = -P_R): kein Austauschdiagramm überlebt ({u_nmhv[2]}), NMHV = 0 "
      f"(Abweichung {u_nmhv[0]:.2f}); MHV bleibt richtig, weil es keinen Austausch braucht "
      f"({kontr[('ungespiegelt', 5, 2)][0]:.0e})",
      h_mhv[1] > 1e-2 and min(h_mhv[2]) >= 1 and h_nmhv[2] == [0] and abs(h_nmhv[0] - 1) < 1e-12 and u_nmhv[2] == [0]
      and abs(u_nmhv[0] - 1) < 1e-12 and kontr[("ungespiegelt", 5, 2)][0] < 1e-10)

say("  (c) Negativkontrollen: Normierung und Referenz")
norm_kontr = {}
for name, N in (("KMS modular e^(-4 pi)", Q_NAT**2), ("KMS der Naht 1 - c(1/4)", 1 - C14), ("Fouriermaß 4 pi^2", 4*np.pi**2),
                ("Tomita-Phase pi (N = -1) ohne Umphasung", -1.0)):
    devs = []
    for h in NMHV5:
        for kk in range(3):
            lam, lt = KIN[5][kk]
            A = zusammenbau(lam, lt, list(h), eta_ipi(), N, 1.0, G_NAT)[0]
            devs.append(abs(A - G_NAT**3*REF[(5, h, kk)])/abs(REF[(5, h, kk)]))
    norm_kontr[name] = (min(devs), max(devs))
ref_kontr = []
for h in NMHV5:
    for kk in range(3):
        lam, lt = KIN[5][kk]
        tab = {}

        def eta_je(L, R, tab=tab):
            key = tuple(L)
            if key not in tab:
                tab[key] = eta_ipi()                                                    # eigene b-Polarisation je Austausch
            return tab[key]
        A = zusammenbau(lam, lt, list(h), eta_je, 1.0, 1.0, G_NAT)[0]
        ref_kontr.append(abs(A - G_NAT**3*REF[(5, h, kk)])/abs(REF[(5, h, kk)]))
ZAHLEN["kontrollen_normierung_referenz"] = {k: [r(v[0], 6), r(v[1], 6)] for k, v in norm_kontr.items()} | {
    "Referenz je Austausch (min, max)": [r(min(ref_kontr), 6), r(max(ref_kontr), 6)]}
check("[E negativ] (c) Falsche Normierungen verfehlen alle NMHV-Lagen um genau |N - 1|: KMS modular "
      f"{norm_kontr['KMS modular e^(-4 pi)'][0]:.4f}, KMS der Naht {norm_kontr['KMS der Naht 1 - c(1/4)'][0]:.4f}, "
      f"Fouriermaß 4 pi^2 {norm_kontr['Fouriermaß 4 pi^2'][0]:.2f}, N = -1 ohne Umphasung "
      f"{norm_kontr['Tomita-Phase pi (N = -1) ohne Umphasung'][0]:.2f}. Eine eigene Referenz je Austausch (die b-Polarisation "
      f"v~ des Tomita-Kerns lokal statt eines globalen Punktes auf I_pi) verfehlt um {min(ref_kontr):.2f} bis "
      f"{max(ref_kontr):.2f}: Die Tomita-Kinetik braucht zusätzlich die globale CSW-Eichung aus I_pi (z111)",
      all(abs(v[0] - abs(N_ - 1)) < 1e-9*max(1, abs(N_ - 1)) and abs(v[1] - abs(N_ - 1)) < 1e-9*max(1, abs(N_ - 1)) for v, N_ in
          zip(norm_kontr.values(), (Q_NAT**2, 1 - C14, 4*np.pi**2, -1.0))) and min(ref_kontr) > 1e-3)
laufzeit()

say("  (d) Bilanz")
STATUS = {
    "Paarung": ("E", "omega = <J psi|phi> nichtentartet zwischen a (n = 0, D = -2) und b (n = -4, D = +2), Betrag 1 je Stufe, "
                     "symmetrisch; Hilbertform und reine Konjugation verbinden a und b nicht"),
    "Realform": ("E", "J E J = -(tau E tau)^+ für alle 16 Weyl-geordneten Erzeuger; tau = Theta_1 (z84) = Punktspiegelung x -> -x, "
                      "J = CPT der Leiter; naiv invariant nur Lorentz + Dilatation + u(1) (8); No-go für unverdrehte Invarianz"),
    "Eindeutigkeit": ("E", "verdreht invariante a-b-Formen: id 0, Ad eta 0, Ad tau 1 (= omega), Ad Theta_n 1 (nicht omega)"),
    "Čech-Bild": ("E/L", "omega = Res(Č psi . tau^* Č phi) exakt; die Residuen-(Serre-)Paarung O(0) x O(-4) ist exakt "
                         "sl(4)-invariant und bis auf den Faktor eindeutig, also die On-shell-Gestalt von int D^3Z tr(b dbar a) (Lesart)"),
    "OS-Fortsetzung": ("E", "J = theta o tau (theta: OS-Spiegelung, euklidisch radiale Inversion); Hilbertform = G o theta, "
                            "omega = G o tau; Penrose-Transformation = omega mit gespiegeltem Punktzustand; Röhrenkern 1/det(X1 - X2) bzw. "
                            "(u~ (X2 - X1)^-1 v~)^2/det(X1 - X2), translationsinvariant; ungespiegelt 1/det(Y1 + Y2) nicht"),
    "CSW-Propagator": ("E", "dbar^-1 in CSW-Eichung zwischen L_x, L_y: Träger lambda* = (x - y)^-1 eta, Gewicht 1/det(x - y); "
                            "a-b-Kern hat genau diese Gestalt mit eta = v~; Impulsraum <u lambda_P>^2/P^2, lambda_P = P|v]; "
                            "eine globale Referenz (I_pi) ist zusätzlich nötig"),
    "Normierung": ("E/C", "Faktorisierung: c_V N = 1, Vorzeichen = Umphasung; |N| = 1 aus der Antiunitarität von J, |c_V| = 1 "
                          "unter P_norm (z110); Niveau kürzt sich, KMS e^(-+4 pi) bzw. 0,959 falsch, Uhr trivial, Tomita-Phase nur "
                          "Konvention, 4 pi^2 Fouriermaß"),
    "Zusammenbau": ("E/C", "CSW = BCFW für alle 4- und 5-Punkt-MHV- und NMHV-Lagen mit nativer Paarung, unter T (Vertizes), "
                           "P_N4 (Zähler, Orientierung nimmt c_V = -1 auf) und P_norm; Hilbert- und ungespiegelte Paarung, falsche "
                           "Normierung und lokale Referenz scheitern"),
    "Austausch selbst": ("O", "dass zwischen zwei Nähten tatsächlich ein gaußsches Feld mit dieser Kovarianz ausgetauscht wird, "
                              "liefert keine native Struktur (z111: gemeinsame CAR ohne Kreuzkern); Prämisse P_X statt P_BF"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
urteil = ("Die Kinetik des Austauschfeldes ist nativ bis auf die Spiegelung, und die Spiegelung ist selbst nativ: Die Tomita-"
          "Paarung omega(psi, phi) = <J psi|phi> ist die b-a-Kinetik, aber nicht in der Form J X J = -X^+, die für su(2,2) "
          "unmöglich ist (zwischen zwei Sektoren positiver Energie gibt es keine invariante Bilinearform). Exakt gilt "
          "J E J = -(tau E tau)^+ für alle 16 Weyl-geordneten Erzeuger, mit tau = Theta_1 aus z84, der Punktspiegelung "
          "x -> -x; J ist die CPT-Abbildung der Leiter, naiv invariant nur unter Lorentz, Dilatation und u(1). Die tau-"
          "verdreht invariante a-b-Form ist eindeutig und gleich omega (innere Verdrehung oder keine: 0 Lösungen, anderes "
          "Theta_n: eine andere Form), sie ist nichtentartet mit Betrag 1 je Stufe, und im Čech-Bild ist sie genau die "
          "Residuen-(Serre-)Paarung O(0) x O(-4) nach tau^*, die exakt sl(4)-invariant ist, also (Lesart, eindeutig bis auf "
          "den Faktor) die On-shell-Gestalt von int D^3Z tr(b dbar a). Die OS-Fortsetzung liefert den Propagator: J = theta o tau "
          "mit der OS-Spiegelung theta (Röhrenkonjugation, euklidisch die radiale Inversion aus z84); die Hilbertform ist der "
          "holomorphe Kern G an der theta-Spiegelung, die Tomita-Paarung derselbe Kern an der tau-Spiegelung. Die Penrose-"
          "Transformation ist die Tomita-Paarung mit dem tau-gespiegelten Punktzustand, in der Röhre ist der Kern 1/det(X1 - X2) (z111s Cauchykern, jetzt mit der "
          "Büschelnormierung <ab> = 1 aus der Fock-Normierung) und im a-b-Sektor (u~ (X2 - X1)^-1 v~)^2/det(X1 - X2), genau "
          "die Gestalt von dbar^-1 in CSW-Eichung zwischen zwei Geraden (Träger lambda* = (x - y)^-1 eta, Gewicht 1/det(x - y)) mit der "
          "Referenz eta = b-Polarisation; im Impulsraum <u lambda_P>^2/P^2 mit lambda_P = P|v] (4 pi^2/P^2 aus z111 ist nur "
          "das Fouriermaß). Die Normierung ist kein freier Parameter mehr: Faktorisierung verlangt c_V N = 1, das Vorzeichen "
          "ist eine Umphasung der negativen Helizität, und |N| = 1 folgt aus der Antiunitarität von J, sobald die "
          "Vertexzustände Fock-normiert sind (P_norm aus z110, |c_V| = 1). Von den vorgeschlagenen nativen Größen fixiert keine "
          "den Betrag: Das Niveau kürzt sich (es bestimmt nur g = -k^(-1/2)), die KMS-Gewichte geben e^(-+4 pi) bzw. "
          "1 - c(1/4) = 0,959 und scheitern, die Uhr ist auf a und b gleich 1, die Tomita-Phase legt J Omega = Omega nur das "
          "konventionelle Vorzeichen fest. Zusammengebaut (Nahtvertizes unter T und P_N4, Tomita-Kinetik, I_pi-Eichung) "
          "entsteht CSW exakt: alle 4- und 5-Punkt-MHV- und NMHV-Lagen stimmen mit BCFW überein, die Dreiminus-Amplituden "
          "bei vier Punkten heben sich weg, die Summe ist eta-unabhängig. Die Paarung ohne J (Hilbertform innerhalb eines "
          "Sektors, gleiche Helizität an beiden Enden) erzeugt unechte MHV-Diagramme und tötet NMHV, die ungespiegelte Tomita-Paarung hat den falschen Impuls"
          "träger, eine lokale Referenz je Austausch scheitert. Urteil: bedingt nativ mit benannten Prämissen. Nativ (E) sind "
          "Form, Eindeutigkeit, Spiegelung und Propagator, der Betrag der Normierung unter P_norm (E/C); Prämissen bleiben T (Naht = Gerade, "
          "Vertizes), P_N4 (Zähler, dessen Orientierung das Vorzeichen aufnimmt), P_norm und neu P_X: dass zwischen zwei "
          "Nähten tatsächlich ein freies Feld mit dieser Kovarianz ausgetauscht wird; die gemeinsame CAR hat weiterhin keinen "
          "Kreuzkern (z111). P_BF aus z111 schrumpft damit auf P_X: Kinetik und Normierung sind nicht mehr frei. Folgen: K15 "
          "bleibt L, präzisiert: Yang-Mills-Bäume = Nahtvertizes (T, P_N4) + Tomita-Kinetik (nativ) + I_pi (nativ) + P_X. "
          "Die Szenarioachse Wechselwirkung behält den bedingten Wert 'CSW-Komposition', jetzt mit den Prämissen T, P_N4, "
          "P_norm, P_X und ohne freie Normierung. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
say(f"\n   (Laufzeit {time.time() - T0:.1f} s)")
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z123_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
