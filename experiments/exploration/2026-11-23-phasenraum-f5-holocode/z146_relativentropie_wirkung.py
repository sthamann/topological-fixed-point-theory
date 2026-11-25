"""Z146: Relativentropie-Wirkung mit den nativen Tomita-Daten: Legt die Fisher-Form die Gravitationskopplung fest?

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Firewalls: GRAV.NONCIRCULAR.01 (kein G = 1, keine Flächenregel als Eingabe; Importe werden als
solche genannt), DIMENSION.UPLIFT.FIREWALL.01 (1D-Aussagen sind als 1D typisiert). Statuskürzel: E exakt, C bedingt
(Prämisse genannt), L Lesart, N numerisch, O offen.

Frage: Gibt die zweite Variation der relativen Entropie (Bogoliubov-Kubo-Mori-/Fisher-Form) des nativen KMS-Zustands
(beta = 4 pi, q = e^(-2 pi)) auf der nativen Deformationsfamilie eine quadratische Form, die (a) die linearisierte
MacDowell-Mansouri-/Einstein-Form um das native AdS-Vakuum ist und (b) eine eindeutige Zahl G/R_Vak^2 liefert?
(Idee: Lashkari-Van Raamsdonk, Fisher-Information = kanonische Energie; Hollands-Wald.)

Bestand: z115 (Familie J(M_g, psi) bei festem Zustand; BKM der reduzierten Zustände nur auf den Boosts, Wert 2), z128
(BKM vom Rang 1 längs der Uhr, 16 q^2), z136 (erster Hauptsatz für K = 4 pi L_0, Möbius-Fisher 0,2618; Modulgenerator
elliptisch; Überbestimmungstabelle für G/R_Vak^2), z127 (Uhrspur: G/R_Vak^2 = 1/(4096 pi kappa); 2 c_3 = 1/beta), z137
(C_T = 1/(3 pi^4)), z59 (Virasoro-Antworten: c = 1, h = 1/32 je Quarter), z68 (Summe k/(e^(k pi) + 1) = 1/24), z107
(beta = 4 pi spannungsfrei), z98/z118 (MM-Normierung b = l^2/(64 pi G), l = R_Vak/2).

Zwei Lesarten der 'nativen Tomita-Daten' (beide gerechnet, die Wahl ist benannt):
  T1  Thermofeldschnitt Original|Partner (z83, z86, z136): M = B(H), die native Gruppe wirkt durch Automorphismen von M,
      rho = e^(-beta E)/Z, rho_g = U(g) rho U(g)^*. Zustandsraum: Einteilchen-Leiter (P_1T), Leiterfeld (P_Feld) oder der
      Fockraum eines Quarters unter der radialen sl(2,R) (1D).
  T2  Twistorlesart (z115, z128): M = Algebra der Moden a = (d + q e^+)/s, die Gruppe mischt M und M'; verglichen werden
      die reduzierten Zustände von U(g)^* Omega und Omega auf M.

Teil A  T1: der native Zustand auf der Leiter der su(2,2) und die volle Hesse-Matrix von S(rho_g || rho) auf allen 15
        (bzw. 10) Erzeugern; geschlossene Form S = beta <E> (tr(g^+ g)/4 - 1); Leiterfeld; Kegel der Modulgeneratoren.
Teil B  Zerlegung unter der Stabilisatorgruppe: Rang, Signatur, Eigenwerte, Killing-Proportionalität, Kähler-Struktur;
        Einschränkung auf die AdS4-Tangente, auf Lorentz und auf die Minkowski-Translationen; T2 zum Vergleich.
Teil C  Vergleich mit der Gravitation: invariante Formen (MM, Holst), die native Familie als reine Eichrichtung
        (quadratische MM-Wirkung null), Hollands-Wald auf der Bahn (Masse statt Kopplung), Horizontterm (T < T_min),
        Zentralterme (Virasoro in 1D, C_T in 4D), Randinhalt unter so(3,2); neue Zeilen der Überbestimmungstabelle.
Teil D  Nichtzirkularität (GRAV.NONCIRCULAR.01).
Teil E  Bilanz und Urteil.
"""
import json
import math
import os
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # deterministische Rundung (einfädig)

import mpmath as mp
import numpy as np
import scipy.linalg as sla
import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
mp.mp.dps = 40
rng = np.random.default_rng(146)
BETA = 4*math.pi
X_NATIV = math.exp(-BETA)                                # x = q^2 = e^(-4 pi)
Q_NATIV = math.exp(-2*math.pi)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")
    return bool(cond)


def befund(name, val, text):
    BEFUNDE[name] = val
    say(f"   [Befund: {val}] {text}")


def kl(A, B):
    return A @ B - B @ A


def gr(x):
    """Größenordnung eines Residuums als Text."""
    return f"{float(abs(x)):.0e}"


# ------------------------------------------------------------------ u(2,2) und die native so(3,2) auf T (wie z98, z115)
I2 = np.eye(2)
SIGMA = [np.eye(2, dtype=complex), np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]),
         np.diag([1., -1.]).astype(complex)]
G4 = np.diag([1., 1., -1., -1.]).astype(complex)
uvec = [np.array([1, 0, -1, 0.]), np.array([0, 1, 0, -1.])]
vvec = [np.array([1, 0, 1, 0.]), np.array([0, 1, 0, 1.])]
Inp = (np.outer(uvec[0], uvec[1]) - np.outer(uvec[1], uvec[0])).astype(complex)       # I_pi
Ion = (np.outer(vvec[0], vvec[1]) - np.outer(vvec[1], vvec[0])).astype(complex)       # I_omega
ISYM = Inp - Ion                                                                      # so(3,2) = Der(I_pi - I_omega)
NP_ = [-0.5*sum(s[i, j]*np.outer(uvec[i], vvec[j]) for i in range(2) for j in range(2)) for s in SIGMA]
NK_ = [-0.5*sum(s[i, j]*np.outer(vvec[i], uvec[j]) for i in range(2) for j in range(2)) for s in SIGMA]
ETA = np.diag([1., -1., -1., -1.])
MP = [1j*N for N in NP_]                                                              # P_mu
MK = [ETA[m, m]*1j*NK_[m] for m in range(4)]                                          # K_mu
QP = [(MP[m] + MK[m])/2 for m in range(4)]                                            # Q^+_mu, Q^+_0 = E (z98)
QM = [(MP[m] - MK[m])/2 for m in range(4)]                                            # Q^-_mu, Q^-_0 = Kappe (z104)
E_T = QP[0]
Q_T = QP[1:]
B_T = [kl(E_T, X) for X in Q_T]                                                       # Boosts
R_T = [kl(Q_T[1], Q_T[2]), kl(Q_T[2], Q_T[0]), kl(Q_T[0], Q_T[1])]                    # Drehungen
D_T = kl(E_T, QM[0])                                                                  # Dilatation
SO32 = [E_T] + R_T + Q_T + B_T
SU22 = SO32 + QM[1:] + [QM[0], D_T]
NAMEN = ["E", "R1", "R2", "R3", "Q1", "Q2", "Q3", "B1", "B2", "B3", "Q-1", "Q-2", "Q-3", "Q-0", "D"]
IDX = {n: i for i, n in enumerate(NAMEN)}
K_IDX = [0, 1, 2, 3, 10, 11, 12]                                                      # Stabilisator von E in su(2,2)
P_IDX = [4, 5, 6, 7, 8, 9, 13, 14]                                                    # Nebenraum p in su(2,2)


def trf(X, Y):
    return float(np.trace(X @ Y).real)


def koord(M, basis):
    """Koordinaten von M in einer spurorthogonalen Basis."""
    return np.array([trf(M, X)/trf(X, X) for X in basis])


def ad_matrizen(basis):
    return [np.array([koord(kl(X, Y), basis) for Y in basis]).T for X in basis]


def invariante_formen(ads, n):
    """Raum der symmetrischen Bilinearformen S auf R^n mit ad^T S + S ad = 0 für alle ad."""
    paare = [(i, j) for i in range(n) for j in range(i, n)]
    basis_s = []
    for (i, j) in paare:
        S = np.zeros((n, n))
        S[i, j] = S[j, i] = 1.0
        basis_s.append(S)
    rows = np.array([np.concatenate([(a.T @ S + S @ a).ravel() for a in ads]) for S in basis_s]).T
    _, sv, vh = np.linalg.svd(rows)
    rang = int(np.sum(sv > 1e-9*max(1.0, sv[0])))
    return [sum(c*S for c, S in zip(v, basis_s)) for v in vh[rang:]]


def projektion(F, formen):
    """Frobenius-Projektion von F auf den Raum der Formen; gibt Projektion und relativen Rest."""
    if not formen:
        return np.zeros_like(F), 1.0
    A = np.array([S.ravel() for S in formen]).T
    c, *_ = np.linalg.lstsq(A, F.ravel(), rcond=None)
    P = (A @ c).reshape(F.shape)
    return P, float(np.linalg.norm(F - P)/np.linalg.norm(F))


# ------------------------------------------------------------------ Leiter der su(2,2): Helizität 0 im Fockraum von d, e
def leiter(nmax):
    """Leiterdarstellung (Delta = 1): Zustände |n1, n2; m1, m2> mit n1 + n2 = m1 + m2 = N <= nmax, E = N + 1.
    Z = (d_1, d_2, e_1^+, e_2^+), Zhat(X) = Ztilde X Z mit Ztilde = (d^+, -e); [Zhat(X), Zhat(Y)] = Zhat([X, Y])."""
    basis = [(n1, N - n1, m1, N - m1) for N in range(nmax + 1) for n1 in range(N + 1) for m1 in range(N + 1)]
    idx = {b: i for i, b in enumerate(basis)}
    dim = len(basis)

    def op(fn):
        M = np.zeros((dim, dim), complex)
        for i, b in enumerate(basis):
            r = fn(list(b))
            if r is not None and tuple(r[0]) in idx:
                M[idx[tuple(r[0])], i] += r[1]
        return M

    def tausch(off, i, j):
        def f(b):
            if b[off + j] == 0:
                return None
            c = math.sqrt(b[off + j])
            b[off + j] -= 1
            c *= math.sqrt(b[off + i] + 1)
            b[off + i] += 1
            return b, c
        return op(f)

    def hoch(i, j):                                                  # d_i^+ e_j^+
        def f(b):
            c = math.sqrt((b[i] + 1)*(b[2 + j] + 1))
            b[i] += 1
            b[2 + j] += 1
            return b, c
        return op(f)

    DD = [[tausch(0, i, j) for j in range(2)] for i in range(2)]
    EE = [[tausch(2, i, j) for j in range(2)] for i in range(2)]
    UP = [[hoch(i, j) for j in range(2)] for i in range(2)]
    Id = np.eye(dim)

    def zhat(X):
        M = np.zeros((dim, dim), complex)
        for i in range(2):
            for j in range(2):
                M += X[i, j]*DD[i][j] + X[i, 2 + j]*UP[i][j]
                M -= X[2 + i, j]*UP[j][i].conj().T                   # e_i d_j
                M -= X[2 + i, 2 + j]*(EE[j][i] + (i == j)*Id)        # e_i e_j^+
        return M
    return np.array([b[0] + b[1] for b in basis]), zhat


def kms(lev, beta):
    p = np.exp(-beta*(lev + 1.0))
    return p/p.sum()


def hesse(p, Hd, gens, beta):
    """Hess(a, b) = beta/2 Tr rho ([[H, X_a], X_b] + [[H, X_b], X_a]) für rho = diag(p), H = diag(Hd)."""
    n = len(gens)
    C = [Hd[:, None]*X - X*Hd[None, :] for X in gens]                 # [H, X_a]
    M = np.zeros((n, n))
    for a in range(n):
        pC, Cp = p[:, None]*C[a], C[a]*p[None, :]
        for b in range(n):
            M[a, b] = beta*float(((pC - Cp)*gens[b].T).sum().real)
    return (M + M.T)/2


def bkm_spektral(p, gens):
    """BKM-Form der Tangenten [X, rho]: Summe Re(conj(X_a,ij) X_b,ij) (p_i - p_j)(ln p_i - ln p_j)."""
    lp = np.log(p)
    W = (p[:, None] - p[None, :])*(lp[:, None] - lp[None, :])
    n = len(gens)
    return np.array([[float((np.conj(gens[a])*gens[b]*W).sum().real) for b in range(n)] for a in range(n)])


# ====================================================================================================================
say("Teil A: T1 -- der native KMS-Zustand auf der Leiter und die volle Hesse-Matrix der relativen Entropie")
# --- A1 Erzeuger und Darstellung
herm = max(np.abs(X.conj().T @ G4 + G4 @ X).max() for X in SU22)
spur = max(abs(np.trace(X)) for X in SU22)
rang15 = int(np.linalg.matrix_rank(np.array([np.concatenate([X.real.ravel(), X.imag.ravel()]) for X in SU22])))
BT = np.array([[trf(X, Y) for Y in SU22] for X in SU22])
bt_diag = np.round(np.diag(BT), 12).tolist()
der_I = [float(np.abs(X @ ISYM + ISYM @ X.T).max()) for X in SU22]
abschluss = max(np.abs(kl(X, Y) - sum(c*Z for c, Z in zip(koord(kl(X, Y), SO32), SO32))).max()
                for X in SO32 for Y in SO32)
J_p = max(np.abs(kl(E_T, kl(E_T, SU22[a])) + SU22[a]).max() for a in P_IDX)
k_komm = max(np.abs(kl(E_T, SU22[a])).max() for a in K_IDX)
ZAHLEN["erzeuger"] = {"Namen": NAMEN, "Spurform tr_T(X_a X_b) diagonal": bt_diag,
                      "Der(I_pi - I_omega) je Erzeuger": [round(x, 12) for x in der_I]}
check(f"[E] Erzeuger auf dem Twistorraum T (z98, z115): 15 G-antihermitesche, spurfreie Matrizen (Reste {gr(herm)}, "
      f"{gr(spur)}), linear unabhängig (Rang {rang15}), spurorthogonal mit tr_T(X^2) = -1 auf E, den Drehungen R_k und "
      f"Q^-_k (kompakt, 7) und +1 auf Q^+_k = Q_k, den Boosts B_k = [E, Q_k], der Kappe Q^-_0 und D = [E, Q^-_0] "
      f"(nichtkompakt, 8). Die ersten zehn (E, R_k, Q_k, B_k) sind die native so(3,2) = Der(I_pi - I_omega) (Rest "
      f"{gr(max(der_I[:10]))}, abgeschlossen {gr(abschluss)}); die übrigen fünf liegen nicht darin (>= "
      f"{min(der_I[10:]):.1f}). ad_E hat auf dem Nebenraum p das Quadrat -1 (Rest {gr(J_p)}): E ist auf Vektoren "
      f"einheitsnormiert (Periode 2 pi, z136), eine Einheit in p ist eine Einheit Rapidität bzw. geodätische Länge in "
      f"Einheiten von l", herm < 1e-14 and spur < 1e-14 and rang15 == 15
      and bt_diag == [-1.0]*4 + [1.0]*6 + [-1.0]*3 + [1.0]*2 and np.abs(BT - np.diag(np.diag(BT))).max() < 1e-14
      and max(der_I[:10]) < 1e-14 and min(der_I[10:]) > 1 and abschluss < 1e-14 and J_p < 1e-14 and k_komm < 1e-14)

NMAX = 10
LEV, ZHAT = leiter(NMAX)
DIM = len(LEV)
GEN = [ZHAT(X) for X in SU22]
H_OP = 1j*GEN[0]
HD = LEV + 1.0
tief = LEV < NMAX - 1
hom = max(np.abs((kl(GEN[a], GEN[b]) - ZHAT(kl(SU22[a], SU22[b])))[np.ix_(tief, tief)]).max()
          for a in range(15) for b in range(a + 1, 15))
anti = max(np.abs(X + X.conj().T).max() for X in GEN)
viel = [int(np.sum(LEV == N)) for N in range(5)]
ZAHLEN["leiter"] = {"Stufen bis": NMAX, "Dimension": DIM, "Vielfachheiten E = 1..5": viel}
check(f"[E] Leiter der su(2,2) (z83; Delta = 1): Helizität-0-Sektor des Fockraums von Z = (d_1, d_2, e_1^+, e_2^+), "
      f"Stufen N <= {NMAX} (Dimension {DIM}). Zhat(X) = Ztilde X Z ist antihermitesch ({gr(anti)}) und unterhalb der "
      f"Abschneidestufe ein Lie-Homomorphismus ohne Zentralterm (Rest {gr(hom)}); H = i Zhat(E) = N + 1 ist diagonal "
      f"({gr(np.abs(H_OP - np.diag(HD)).max())}) mit den Vielfachheiten {viel} = m^2 für E = m = 1..5. Der Nullpunkt "
      f"E = 1 ist durch die Kommutatoren festgelegt (halbeinfache Algebra, keine freie Konstante)",
      anti < 1e-13 and hom < 1e-12 and np.abs(H_OP - np.diag(HD)).max() < 1e-12 and viel == [1, 4, 9, 16, 25])

# --- A2 KMS-Zustand, erste Ordnung, geschlossene Form
P_N = kms(LEV, BETA)
RHO = np.diag(P_N)
E_MW = float((P_N*HD).sum())
x_mp = mp.exp(-4*mp.pi)
E_mp = 1 + x_mp/(1 + x_mp) + 3*x_mp/(1 - x_mp)                      # aus Z = x (1 + x)/(1 - x)^3
VAR_mp = x_mp/(1 + x_mp)**2 + 3*x_mp/(1 - x_mp)**2
NF_mp = 4*mp.pi*E_mp


def e_mittel(U):
    return float(np.trace(U @ RHO @ U.conj().T @ H_OP).real)


erste = max(abs(e_mittel(sla.expm(1e-3*X)) - e_mittel(sla.expm(-1e-3*X)))/2e-3 for X in GEN)


def s_rel_direkt(U):
    """S(rho_g || rho) spektral: Tr rho_g ln rho_g - Tr rho_g ln rho."""
    rg = U @ RHO @ U.conj().T
    w = np.linalg.eigvalsh((rg + rg.conj().T)/2)
    w = w[w > 1e-300]
    return float((w*np.log(w)).sum() - (np.diag(rg).real*np.log(P_N)).sum()), float(np.trace(rg @ H_OP).real)


geschl = {}
for nm, s_ in (("Q1", 0.3), ("B1", 0.3), ("D", 0.3), ("Q-0", 0.2)):
    Sd, Eg = s_rel_direkt(sla.expm(s_*GEN[IDX[nm]]))
    geschl[nm] = (Sd, BETA*(Eg - E_MW), BETA*E_MW*(math.cosh(s_) - 1))
koef = rng.normal(size=15)*0.07
g_zuf = sla.expm(sum(c*X for c, X in zip(koef, SU22)))
Sd_z, Eg_z = s_rel_direkt(sla.expm(sum(c*X for c, X in zip(koef, GEN))))
S_form_z = BETA*E_MW*(float(np.trace(g_zuf.conj().T @ g_zuf).real)/4 - 1)
X_gem = (SU22[IDX["Q1"]] + SU22[IDX["B2"]])/math.sqrt(2)
S_gem = BETA*(s_rel_direkt(sla.expm(0.3*ZHAT(X_gem)))[1] - E_MW)
gem_cosh, gem_rang2 = BETA*E_MW*(math.cosh(0.3) - 1), BETA*E_MW*(math.cosh(0.3*math.sqrt(2)) - 1)/2
geschl_ok = all(abs(a - c) < 1e-11 and abs(b - c) < 1e-11 for a, b, c in geschl.values())
ZAHLEN["zustand"] = {"x = e^(-4 pi)": mp.nstr(x_mp, 20), "<E> (Leiter)": mp.nstr(E_mp, 20),
                     "<E> numerisch": round(E_MW, 14), "Var(E)": mp.nstr(VAR_mp, 15),
                     "N_F = beta <E>": mp.nstr(NF_mp, 20),
                     "S(rho_g || rho) (spektral, beta d<E>, Formel)": {k: [round(v, 12) for v in t] for k, t in geschl.items()},
                     "Zufallselement: S spektral, Formel": [round(Sd_z, 11), round(S_form_z, 11)],
                     "gemischte Richtung (Q1 + B2)/sqrt 2, s = 0,3: S, cosh-Form, Rang-2-Form":
                         [round(S_gem, 10), round(gem_cosh, 10), round(gem_rang2, 10)]}
check(f"[E] Nativer Zustand (T1): rho = e^(-4 pi E)/Z auf der Leiter, Z = x (1 + x)/(1 - x)^3 mit x = e^(-4 pi) = "
      f"{mp.nstr(x_mp, 8)}, <E> = 1 + x/(1 + x) + 3x/(1 - x) = {mp.nstr(E_mp, 15)} (numerisch {E_MW:.13f}). Erste "
      f"Ordnung: d<E>/ds = 0 für alle 15 Erzeuger (zentrale Differenz {gr(erste)}). Weil rho_g = U rho U^* dieselbe Entropie hat, "
      f"ist S(rho_g || rho) = beta (<E>_g - <E>) exakt; spektral gerechnet stimmt das (Q1, B1, D bei s = 0,3; Q^-_0 bei "
      f"0,2) mit beta <E> (cosh s - 1) überein: {[round(geschl[k][0], 10) for k in ('Q1', 'B1', 'D', 'Q-0')]}",
      erste < 1e-9 and abs(E_MW - float(E_mp)) < 1e-13 and geschl_ok)
check(f"[E, neu] Geschlossene Form der Relativentropie auf der ganzen Gruppe: S(rho_g || rho) = beta <E> (tr(g^+ g)/4 - "
      f"1) mit g in SU(2,2) auf T (Zufallselement: {Sd_z:.10f} gegen {S_form_z:.10f}). Für Einsrichtungen ist das beta "
      f"<E> (cosh s - 1): die Bewegungsenergie beta M (gamma - 1) eines Boosts bzw. die Lageenergie beta M (cosh rho - 1) "
      f"einer Verschiebung im AdS-Potential, mit M l = <E>. Die Bahn hat Rang 2: für (Q1 + B2)/sqrt 2 gilt statt cosh s - 1 "
      f"(cosh(sqrt 2 s) - 1)/2 (s = 0,3: {S_gem:.9f} gegen {gem_rang2:.9f}; cosh-Form {gem_cosh:.9f})",
      abs(Sd_z - S_form_z) < 1e-10 and abs(S_gem - gem_rang2) < 1e-10 and abs(S_gem - gem_cosh) > 1e-3)

# --- A3 volle Hesse-Matrix
HESS = hesse(P_N, HD, GEN, BETA)
BKM = bkm_spektral(P_N, GEN)
PRED = BETA*E_MW*np.array([[trf(kl(E_T, X), kl(E_T, Y)) for Y in SU22] for X in SU22])
fd = {}
for nm in ("Q1", "B2", "D", "Q-0", "E", "R1", "Q-2"):
    eps = 1e-3
    fd[nm] = 2*BETA*(float(np.trace(sla.expm(eps*GEN[IDX[nm]]) @ RHO @ sla.expm(-eps*GEN[IDX[nm]]) @ H_OP).real)
                     - E_MW)/eps**2
ew15 = np.linalg.eigvalsh(HESS)
ew10 = np.linalg.eigvalsh(HESS[:10, :10])
ZAHLEN["hesse_T1"] = {"Diagonale (15 Erzeuger)": np.round(np.diag(HESS), 9).tolist(),
                      "größtes Nebendiagonalelement": float(np.abs(HESS - np.diag(np.diag(HESS))).max()),
                      "Eigenwerte su(2,2)": np.round(ew15, 9).tolist(), "Eigenwerte so(3,2)": np.round(ew10, 9).tolist(),
                      "|Hesse - BKM spektral|": float(np.abs(HESS - BKM).max()),
                      "|Hesse - beta <E> tr([E,X][E,Y])|": float(np.abs(HESS - PRED).max()),
                      "zweite Differenz (eps = 1e-3)": {k: round(v, 6) for k, v in fd.items()}}
hess_ok = (np.abs(HESS - BKM).max() < 1e-10 and np.abs(HESS - PRED).max() < 1e-10
           and all(abs(fd[k] - HESS[IDX[k], IDX[k]]) < 1e-4 for k in fd))
check(f"[E] Volle Hesse-Matrix von S(rho_g || rho) bei g = 1 auf allen 15 Erzeugern der su(2,2): Hess(X, Y) = beta/2 "
      f"<[[E, X], Y] + [[E, Y], X]> ist gleich der spektral gerechneten BKM-Form der Tangenten [X, rho] (Abweichung "
      f"{gr(np.abs(HESS - BKM).max())}) und gleich beta <E> tr_T([E, X][E, Y]) ({gr(np.abs(HESS - PRED).max())}); zweite "
      f"Differenzen bestätigen sie ({[round(fd[k], 4) for k in ('Q1', 'B2', 'D', 'Q-0', 'E', 'R1', 'Q-2')]}). Sie ist "
      f"diagonal in der Erzeugerbasis: 0 auf E, R_k, Q^-_k und N_F = beta <E> = {mp.nstr(NF_mp, 12)} auf Q_k, B_k, "
      f"Q^-_0 und D. Eigenwerte su(2,2): 7 x 0 und 8 x {ew15[-1]:.8f}; so(3,2): 4 x 0 und 6 x {ew10[-1]:.8f}",
      hess_ok and np.allclose(ew15[:7], 0, atol=1e-9) and np.allclose(ew15[7:], float(NF_mp), atol=1e-9)
      and np.allclose(ew10[:4], 0, atol=1e-9) and np.allclose(ew10[4:], float(NF_mp), atol=1e-9))

# --- A4 Kontrolle bei großem x
LEV_K, ZHAT_K = leiter(12)
GEN_K = [ZHAT_K(X) for X in SU22]
BETA_K = math.log(20.0)
P_K = kms(LEV_K, BETA_K)
E_K = float((P_K*(LEV_K + 1.0)).sum())
x_k = math.exp(-BETA_K)
HESS_K = hesse(P_K, LEV_K + 1.0, GEN_K, BETA_K)
PRED_K = BETA_K*E_K*np.array([[trf(kl(E_T, X), kl(E_T, Y)) for Y in SU22] for X in SU22])
ZAHLEN["kontrolle_x"] = {"x": 0.05, "<E>": round(E_K, 10), "Formel": round(1 + x_k/(1 + x_k) + 3*x_k/(1 - x_k), 10),
                         "beta <E>": round(BETA_K*E_K, 9), "|Hesse - Vorhersage|": float(np.abs(HESS_K - PRED_K).max())}
check(f"[E] Kontrolle bei starker thermischer Besetzung (x = 0,05, Stufen bis 12): <E> = {E_K:.9f} (Formel "
      f"{1 + x_k/(1 + x_k) + 3*x_k/(1 - x_k):.9f}), Hesse-Matrix = beta <E> tr_T([E, X][E, Y]) (Abweichung "
      f"{gr(np.abs(HESS_K - PRED_K).max())}): Die Form ist für jedes beta dieselbe invariante Form, nur die Zahl beta <E> "
      f"ändert sich ({BETA_K*E_K:.6f} statt {float(NF_mp):.6f})",
      abs(E_K - (1 + x_k/(1 + x_k) + 3*x_k/(1 - x_k))) < 1e-9 and np.abs(HESS_K - PRED_K).max() < 1e-8)

# --- A5 Leiterfeld (bosonische zweite Quantisierung) und Kontrolle der quasifreien Formel
nB = 1.0/np.expm1(BETA*HD)
HESS_F = hesse(nB, HD, GEN, BETA)
E_FOCK_mp = mp.fsum(mp.mpf(m)**3*x_mp**m/(1 - x_mp**m) for m in range(1, 60))
qth = mp.exp(-2*mp.pi)                                               # Nome der Thetafunktionen bei tau = 2 i
E4_2i = (mp.jtheta(2, 0, qth)**8 + mp.jtheta(3, 0, qth)**8 + mp.jtheta(4, 0, qth)**8)/2
NF_FOCK_mp = 4*mp.pi*E_FOCK_mp
PRED_F = float(NF_FOCK_mp)*np.array([[trf(kl(E_T, X), kl(E_T, Y)) for Y in SU22] for X in SU22])
# Kontrolle: zwei Bosonmoden, Gamma(u) mit u in U(2), exakt im Raum Gesamtzahl <= 16
NT = 16
bas2 = [(a, n - a) for n in range(NT + 1) for a in range(n + 1)]
ix2 = {b: i for i, b in enumerate(bas2)}
h1 = np.array([1.5, 2.5])
beta_t = 1.2
A1 = 0.4*np.array([[0.3j, 1.0 + 0.5j], [-1.0 + 0.5j, -0.7j]])       # antihermitesch
dG = np.zeros((len(bas2), len(bas2)), complex)
for (a, b), i in ix2.items():
    dG[i, i] += A1[0, 0]*a + A1[1, 1]*b
    if b > 0:
        dG[ix2[(a + 1, b - 1)], i] += A1[0, 1]*math.sqrt((a + 1)*b)
    if a > 0:
        dG[ix2[(a - 1, b + 1)], i] += A1[1, 0]*math.sqrt(a*(b + 1))
En2 = np.array([h1[0]*a + h1[1]*b for a, b in bas2])
p2 = np.exp(-beta_t*En2)
p2 /= p2.sum()
U2 = sla.expm(dG)
r2 = U2 @ np.diag(p2) @ U2.conj().T
w2 = np.linalg.eigvalsh((r2 + r2.conj().T)/2)
w2 = w2[w2 > 1e-300]
S2_dir = float((w2*np.log(w2)).sum() - (np.diag(r2).real*np.log(p2)).sum())
u1 = sla.expm(A1)
n1 = np.diag(1.0/np.expm1(beta_t*h1))
S2_form = beta_t*float(np.trace(n1 @ (u1.conj().T @ np.diag(h1) @ u1 - np.diag(h1))).real)
ZAHLEN["leiterfeld"] = {"<E>_Fock = Summe m^3 x^m/(1 - x^m)": mp.nstr(E_FOCK_mp, 20),
                        "(E_4(2i) - 1)/240 (Thetafunktionen)": mp.nstr((E4_2i - 1)/240, 20),
                        "N_F Fock = beta <E>_Fock": mp.nstr(NF_FOCK_mp, 15),
                        "|Hesse_Fock - Vorhersage|": float(np.abs(HESS_F - PRED_F).max()),
                        "Kontrolle zwei Moden: S direkt, Formel": [round(S2_dir, 10), round(S2_form, 10)]}
check(f"[E] Leiterfeld (P_Feld: bosonische zweite Quantisierung der Leiter, U = Gamma(u)): S(rho_g || rho) = beta "
      f"Tr_1 n_B (u^* E u - E), also dieselbe Form mit <E>_Fock = Summe m^3 x^m/(1 - x^m) = (E_4(2i) - 1)/240 = "
      f"{mp.nstr(E_FOCK_mp, 12)} (Thetafunktionen: {mp.nstr((E4_2i - 1)/240, 12)}): N_F = {mp.nstr(NF_FOCK_mp, 10)} "
      f"(Abweichung der Hesse-Matrix {gr(np.abs(HESS_F - PRED_F).max())}). Kontrolle der quasifreien Formel an zwei "
      f"Moden im vollen Fockraum: {S2_dir:.9f} gegen {S2_form:.9f}",
      abs(E_FOCK_mp - (E4_2i - 1)/240) < mp.mpf(10)**-30 and np.abs(HESS_F - PRED_F).max() < 1e-14
      and abs(S2_dir - S2_form) < 1e-9)

# --- A6 Kegel der Modulgeneratoren: Erzeugerparameter und Temperaturrichtung
lnz = {}
for th in (0.5, 2.0):
    Kop = BETA*H_OP + th*1j*GEN[IDX["Q1"]]
    wk = np.linalg.eigvalsh((Kop + Kop.conj().T)/2)
    b2 = math.sqrt(BETA**2 - th**2)
    x2 = math.exp(-b2)
    lnz[th] = (math.log(np.exp(-wk).sum()) - math.log(X_NATIV*(1 + X_NATIV)/(1 - X_NATIV)**3),
               math.log(x2*(1 + x2)/(1 - x2)**3) - math.log(X_NATIV*(1 + X_NATIV)/(1 - X_NATIV)**3))
SUSZ_mp = E_mp/(4*mp.pi)
VARK_mp = (4*mp.pi)**2*VAR_mp
ZAHLEN["kegel"] = {"ln Z(beta E + theta Q1) - ln Z (numerisch, ln Z(sqrt(beta^2 - theta^2)))":
                   {str(k): [round(a, 12), round(b, 12)] for k, (a, b) in lnz.items()},
                   "Suszeptibilität <E>/beta": mp.nstr(SUSZ_mp, 15), "2 c_3 = 1/(4 pi)": mp.nstr(1/(4*mp.pi), 15),
                   "Verhältnis": mp.nstr(SUSZ_mp*4*mp.pi, 15), "Var(K) = beta^2 Var(E)": mp.nstr(VARK_mp, 12)}
check(f"[E] Kegel der Modulgeneratoren K = beta E + theta X (X in p, Einheit): ln Z hängt nur von der Killing-Norm ab, "
      f"ln Z(K) = ln Z(sqrt(beta^2 - theta^2)) (theta = 0,5 und 2: {[round(lnz[t][0], 10) for t in (0.5, 2.0)]} gegen "
      f"{[round(lnz[t][1], 10) for t in (0.5, 2.0)]}). Dieselbe Fisher-Form im Erzeugerparameter theta = beta s ist die "
      f"Kubo-Mori-Suszeptibilität <E>/beta = {mp.nstr(SUSZ_mp, 12)} = 2 c_3 x <E> (2 c_3 = 1/(4 pi) = "
      f"{mp.nstr(1/(4*mp.pi), 12)}, Faktor {mp.nstr(E_mp, 12)}); in der Temperaturrichtung kommt Var(K) = beta^2 Var(E) = "
      f"{mp.nstr(VARK_mp, 8)} hinzu. Die Familie aller nativen Modulgeneratoren ist damit 9- bzw. 7-dimensional mit der "
      f"Fisher-Form Var(K) (d beta/beta)^2 + N_F ds^2",
      all(abs(a - b) < 1e-11 for a, b in lnz.values()) and abs(SUSZ_mp*4*mp.pi - E_mp) < mp.mpf(10)**-30)
befund("T1: Hesse-Matrix der relativen Entropie", "ja",
       "E: Auf der Bahn des nativen KMS-Zustands unter der nativen Gruppe ist die relative Entropie exakt beta mal die "
       "Änderung der Uhrenergie, S(rho_g || rho) = beta <E> (tr(g^+ g)/4 - 1); ihre Hesse-Matrix ist beta <E> "
       "tr_T([E, X][E, Y]), null auf dem Stabilisator von E und N_F = beta <E> auf allen nichtkompakten Richtungen. Die "
       "einzigen Zahlen darin sind beta und der Erwartungswert der Uhr.")

# ====================================================================================================================
say("\nTeil B: Zerlegung unter der Stabilisatorgruppe; Einschränkungen; Twistorlesart T2")
# --- B1 Stabilisator, Schur, Signatur, Killing-Proportionalität
AD15 = ad_matrizen(SU22)
AD10 = ad_matrizen(SO32)
kern15 = [a for a in range(15) if np.abs(kl(E_T, SU22[a])).max() < 1e-12]
kern10 = [a for a in range(10) if np.abs(kl(E_T, SO32[a])).max() < 1e-12]
p10 = [4, 5, 6, 7, 8, 9]
inv_p15 = invariante_formen([AD15[a][np.ix_(P_IDX, P_IDX)] for a in K_IDX], 8)
inv_p10 = invariante_formen([AD10[a][np.ix_(p10, p10)] for a in kern10], 6)
inv_g10 = invariante_formen(AD10, 10)
inv_g15 = invariante_formen(AD15, 15)
F10 = HESS[:10, :10]
_, rest_kill10 = projektion(F10, inv_g10)
_, rest_kill15 = projektion(HESS, inv_g15)
killing_p = float(np.abs(HESS[np.ix_(P_IDX, P_IDX)] - float(NF_mp)*BT[np.ix_(P_IDX, P_IDX)]).max())
ZAHLEN["zerlegung_T1"] = {"Kern su(2,2)": [NAMEN[a] for a in kern15], "Kern so(3,2)": [NAMEN[a] for a in kern10],
                          "Rang, Signatur su(2,2)": [8, "(8, 0)"], "Rang, Signatur so(3,2)": [6, "(6, 0)"],
                          "K-invariante symmetrische Formen auf p (su(2,2), so(3,2))": [len(inv_p15), len(inv_p10)],
                          "G-invariante Formen auf g (su(2,2), so(3,2))": [len(inv_g15), len(inv_g10)],
                          "relativer Rest gegen die Killing-Form auf ganz g": [round(rest_kill15, 9), round(rest_kill10, 9)]}
check(f"[E] Zerlegung (T1): Der Kern der Fisher-Form ist genau der Stabilisator von rho, k = Zentralisator von E: "
      f"{[NAMEN[a] for a in kern15]} in su(2,2) (s(u(2) + u(2)), Dimension 7), {[NAMEN[a] for a in kern10]} in so(3,2) "
      f"(so(3) + so(2), Dimension 4). Auf dem Nebenraum p (Dimension 8 bzw. 6) ist sie positiv definit, Signatur (8, 0) "
      f"bzw. (6, 0), mit einem einzigen Eigenwert N_F, und gleich N_F mal der Spurform ({gr(killing_p)}). Das erzwingt "
      f"schon Schur: p ist unter K irreduzibel, der Raum der K-invarianten symmetrischen Formen auf p hat Dimension "
      f"{len(inv_p15)} bzw. {len(inv_p10)}. Frei ist nur die Zahl N_F. Auf ganz g ist die Form nicht proportional zur "
      f"Killing-Form (relativer Rest {rest_kill15:.4f} bzw. {rest_kill10:.4f}): sie ist K-invariant, nicht G-invariant",
      kern15 == K_IDX and kern10 == [0, 1, 2, 3] and killing_p < 1e-9 and len(inv_p15) == 1 and len(inv_p10) == 1
      and len(inv_g15) == 1 and len(inv_g10) == 1 and rest_kill10 > 0.1 and rest_kill15 > 0.1)

# --- B2 Kähler-Struktur, Krümmung
OMEGA_F = np.array([[BETA*float((1j*np.trace(RHO @ kl(GEN[a], GEN[b]))).real) for b in range(15)] for a in range(15)])
J15 = AD15[0]
kaehler = float(np.abs(HESS + OMEGA_F @ J15).max())
antis = float(np.abs(OMEGA_F + OMEGA_F.T).max())


def schnitt(X, Y):
    c = kl(X, Y)
    return trf(c, c)/(trf(X, X)*trf(Y, Y) - trf(X, Y)**2)


kr = {"(Q1, Q2)": schnitt(Q_T[0], Q_T[1]), "(Q1, B1)": schnitt(Q_T[0], B_T[0]), "(B1, B2)": schnitt(B_T[0], B_T[1]),
      "(Q1, B2)": schnitt(Q_T[0], B_T[1]), "(Q-0, D)": schnitt(QM[0], D_T)}
ZAHLEN["kaehler"] = {"|Hess + omega J|": kaehler, "omega antisymmetrisch": antis,
                     "Schnittkrümmungen der Spurform": {k: round(v, 12) for k, v in kr.items()},
                     "Schnittkrümmung der Fisher-Metrik": "Spurform / N_F"}
check(f"[E] Kähler-Struktur: omega(X, Y) = i beta <[X, Y]> ist die Kirillov-Kostant-Form der Bahn durch beta <E> E "
      f"(antisymmetrisch {gr(antis)}), J = ad_E die komplexe Struktur, und Hess(X, Y) = -omega(X, J Y) (Rest "
      f"{gr(kaehler)}). Die Fisher-Geometrie ist der hermitesch-symmetrische Raum G/K mit der Bergman-Metrik in der "
      f"Normierung N_F: für so(3,2) der sechsdimensionale Raum SO(3,2)/(SO(3) x SO(2)), d. h. der Raum der negativen "
      f"2-Ebenen in R^(3,2) bzw. der zeitartigen Geodäten (Uhrweltlinien) von AdS4 [L für die Benennung]; für su(2,2) die "
      f"achtdimensionale Zukunftsröhre. Schnittkrümmungen der Spurform: {dict((k, round(v, 6)) for k, v in sorted(kr.items()))}"
      f", in der Fisher-Metrik durch N_F geteilt; Rang 2 (flache Ebene (Q1, B2)). Als Cartan-Geometrie ist eine Bahn "
      f"Maurer-Cartan-flach (z100): Das ist Krümmung des Zustandsraums, keine der Raumzeit",
      kaehler < 1e-9 and antis < 1e-12 and abs(kr["(Q1, Q2)"] + 1) < 1e-12 and abs(kr["(Q1, B1)"] + 1) < 1e-12
      and abs(kr["(Q1, B2)"]) < 1e-12)

# --- B3 Einschränkungen: AdS4-Tangente, Lorentz, Minkowski-Translationen
m_idx = [0, 4, 5, 6]
h_idx = [1, 2, 3, 7, 8, 9]
F_m = np.round(np.diag(F10[np.ix_(m_idx, m_idx)])/float(NF_mp), 9).tolist()
B_m = np.round(np.diag(BT[np.ix_(m_idx, m_idx)]), 9).tolist()
F_h = np.round(np.diag(F10[np.ix_(h_idx, h_idx)])/float(NF_mp), 9).tolist()
B_h = np.round(np.diag(BT[np.ix_(h_idx, h_idx)]), 9).tolist()
GEN_P = [ZHAT(X) for X in MP]
F_PP = hesse(P_N, HD, GEN_P, BETA)/float(NF_mp)
mm_h = max(np.abs(kl(SO32[a], SO32[b]) - sum(c*SO32[k] for c, k in zip(koord(kl(SO32[a], SO32[b]), SO32)[h_idx], h_idx))).max()
           for a in m_idx for b in m_idx)
ZAHLEN["einschraenkung"] = {"Fisher/N_F auf m = (E, Q1, Q2, Q3)": F_m, "Spurform auf m": B_m,
                            "Fisher/N_F auf h = (R_k, B_k)": F_h, "Spurform auf h": B_h,
                            "Fisher/N_F auf den Translationen P_mu": np.round(F_PP, 9).tolist()}
check(f"[E, negativ] Einschränkungen: Auf der AdS4-Tangente m = (E, Q_1, Q_2, Q_3) ([m, m] in Lorentz: {gr(mm_h)}) ist "
      f"die Fisher-Form N_F x {F_m} gegen die Spurform {B_m}: entartet längs der Zeit, die räumliche Projektion "
      f"eta + u u mit u = E. Auf der Lorentz-Algebra h = (R_k, B_k) ist sie N_F x {F_h} gegen {B_h}: Sie lebt auf "
      f"Richtungen, die den AdS-Punkt festhalten (Boosts), und verschwindet auf der Zeit. Auf den Minkowski-Translationen "
      f"P_mu ist sie N_F x delta_(mu nu) (Diagonale {np.round(np.diag(F_PP), 9).tolist()}, Nebendiagonale "
      f"{gr(np.abs(F_PP - np.diag(np.diag(F_PP))).max())}): euklidisch, nicht eta. Der Zustand zeichnet sein Ruhesystem "
      f"aus", F_m == [0.0, 1.0, 1.0, 1.0] and B_m == [-1.0, 1.0, 1.0, 1.0] and F_h == [0.0]*3 + [1.0]*3
      and B_h == [-1.0]*3 + [1.0]*3 and np.abs(F_PP - np.eye(4)).max() < 1e-9 and mm_h < 1e-13)


# --- B4 Twistorlesart T2 (z115): reduzierte Zustände auf M, Gruppe mischt M und M'
def t2_daten(q):
    s = math.sqrt(1 - q*q)
    A0 = np.hstack([I2, q*I2]).astype(complex)/s                     # a = (d + q e^+)/s
    C0 = np.diag([0, 0, 1., 1.]).astype(complex)                     # <Z_a^+ Z_b> im Vakuum der d, e
    nb = q*q/(1 - q*q)
    dn = [(A0 @ X).conj() @ C0 @ A0.T + A0.conj() @ C0 @ (A0 @ X).T for X in SU22]
    Hs = np.array([[float(np.trace(a @ b).real) for b in dn] for a in dn])/(nb*(1 + nb))

    def srel(Rm):
        A = A0 @ Rm
        n1 = A.conj() @ C0 @ A.T
        w = np.linalg.eigvalsh((n1 + n1.conj().T)/2)
        return float(sum(v*math.log(v/nb) - (1 + v)*math.log1p((v - nb)/(1 + nb)) for v in w))
    return Hs, srel


H2_test, srel_test = t2_daten(0.25)
H2_nat, srel_nat = t2_daten(Q_NATIV)
fd2 = {nm: 2*srel_test(sla.expm(1e-4*SU22[IDX[nm]]))/1e-8 for nm in ("B1", "D", "Q1", "Q-0", "E")}
fd2n = {nm: 2*srel_nat(sla.expm(2e-5*SU22[IDX[nm]]))/4e-10 for nm in ("B1", "D")}
gross = {nm: 2*srel_nat(sla.expm(1e-3*SU22[IDX[nm]]))/1e-6 for nm in ("B1", "D", "Q1")}
ew2 = np.linalg.eigvalsh(H2_nat)
ZAHLEN["hesse_T2"] = {"Diagonale (q nativ)": np.round(np.diag(H2_nat), 9).tolist(),
                      "Diagonale (q = 0,25)": np.round(np.diag(H2_test), 9).tolist(),
                      "Eigenwerte": np.round(ew2, 9).tolist(), "Rang so(3,2), su(2,2)": [3, 4],
                      "zweite Differenz q = 0,25 (eps = 1e-4)": {k: round(v, 5) for k, v in fd2.items()},
                      "zweite Differenz q nativ (eps = 2e-5)": {k: round(v, 3) for k, v in fd2n.items()},
                      "2 S/s^2 bei s = 1e-3 (q nativ)": {k: round(v, 4) for k, v in gross.items()}}
check(f"[E] Twistorlesart T2 (z115 reproduziert und auf su(2,2) erweitert): Die Hesse-Matrix von S(U(g)^* Omega|M || "
      f"Omega|M) ist diagonal mit dem Wert 2 auf den drei Boosts und auf D und 0 sonst (Diagonale "
      f"{np.round(np.diag(H2_nat), 6).tolist()}), unabhängig von q (Abweichung zwischen q = 0,25 und nativ "
      f"{gr(np.abs(H2_test - H2_nat).max())}): Rang 3 auf so(3,2), Rang 4 auf su(2,2), Signatur (3, 0) bzw. (4, 0). "
      f"Zweite Differenzen: q = 0,25 {dict((k, round(v, 4)) for k, v in sorted(fd2.items()))}, nativ "
      f"{dict((k, round(v, 2)) for k, v in sorted(fd2n.items()))}. Unter dem Stabilisator SO(3) des Paars (M, Omega) "
      f"zerfällt so(3,2) in 1 + 3 + 3 + 3 (E, R, Q, B) mit den Werten 0, 0, 0, 2: nicht proportional zur Killing-Form "
      f"auf dem Nebenraum (die Q und B gleich gewichtet), die Familie ist kein homogener Raum der nativen Gruppe. Der "
      f"quadratische Bereich endet bei |s| ~ q = e^(-2 pi) (Verschränkungshorizont, z115): bei s = 1e-3 ist 2 S/s^2 schon "
      f"{dict((k, round(v, 3)) for k, v in sorted(gross.items()))}",
      np.abs(H2_test - H2_nat).max() < 1e-9 and np.allclose(np.diag(H2_nat), [0]*7 + [2]*3 + [0]*4 + [2], atol=1e-9)
      and np.abs(H2_nat - np.diag(np.diag(H2_nat))).max() < 1e-9 and abs(fd2["B1"] - 2) < 1e-3 and abs(fd2["D"] - 2) < 1e-3
      and fd2["Q1"] < 1e-3 and abs(fd2n["B1"] - 2) < 0.02 and abs(fd2n["D"] - 2) < 0.02)
befund("Zerlegung und Killing-Proportionalität", "ja (T1) / nein (T2)",
       "E: In T1 ist die Fisher-Form die Bergman-Metrik des hermitesch-symmetrischen Raums G/K (K = Stabilisator des "
       "Zustands, maximal kompakt), Rang 6 bzw. 8, positiv definit, N_F = beta <E> mal Spurform; in T2 hat sie Rang 3 "
       "bzw. 4 mit dem festen Wert 2 und ist nicht Killing-proportional. Keine der beiden ist eine Form auf AdS4 = "
       "G/Lorentz: T1 verschwindet auf der Zeit und lebt auf den Boosts, T2 lebt nur auf den Boosts (und D).")

# ====================================================================================================================
say("\nTeil C: Vergleich mit der Gravitation")
# --- C1 invariante Formen unter Lorentz: MM (eps), Holst (Spurform)
inv_h = invariante_formen([AD10[a] for a in h_idx], 10)
inv_hh = invariante_formen([AD10[a][np.ix_(h_idx, h_idx)] for a in h_idx], 6)
EPS_h = np.zeros((6, 6))
for k in range(3):
    EPS_h[k, 3 + k] = EPS_h[3 + k, k] = 1.0                           # eps(R_k, B_k): die MM-Paarung eps_abcd auf so(3,1)
eps_inv = max(float(np.abs(AD10[a][np.ix_(h_idx, h_idx)].T @ EPS_h + EPS_h @ AD10[a][np.ix_(h_idx, h_idx)]).max())
              for a in h_idx)
Bh_mat = np.diag([-1.0]*3 + [1.0]*3)
Fh_mat = F10[np.ix_(h_idx, h_idx)]/float(NF_mp)
_, rest_h = projektion(Fh_mat, [Bh_mat, EPS_h])
c_eps = float((Fh_mat*EPS_h).sum()/(EPS_h*EPS_h).sum())
c_holst = float((Fh_mat*Bh_mat).sum()/(Bh_mat*Bh_mat).sum())
_, rest_g = projektion(F10, inv_h)
Ad_b = sla.expm(0.5*AD10[IDX["B1"]])
Ad_r = sla.expm(0.5*AD10[IDX["R1"]])
Ad_e = sla.expm(0.5*AD10[IDX["E"]])
lor_def = float(np.linalg.norm(Ad_b.T @ F10 @ Ad_b - F10)/np.linalg.norm(F10))
k_def = max(float(np.linalg.norm(A_.T @ F10 @ A_ - F10)/np.linalg.norm(F10)) for A_ in (Ad_r, Ad_e))
ZAHLEN["invariante_formen"] = {"Lorentz-invariante symmetrische Formen auf so(3,2)": len(inv_h),
                               "auf so(3,1)": len(inv_hh), "eps-Form invariant (Rest)": eps_inv,
                               "Fisher|h: Anteil eps (MM), Anteil Spurform (Holst)": [round(c_eps, 12), round(c_holst, 12)],
                               "relativer Rest auf h": round(rest_h, 9), "relativer Rest auf g": round(rest_g, 9),
                               "Defekt unter Boost e^(0,5 B1)": round(lor_def, 9), "Defekt unter K": k_def}
check(f"[E, negativ] (a) Die Fisher-Form ist nicht die MacDowell-Mansouri-Form: Die quadratische MM-Wirkung ist "
      f"b eps_abcd F^ab F^cd, ihre Fasermetrik die Paarung eps(R_k, B_k) auf Lorentz (invariant: {gr(eps_inv)}). Die "
      f"Lorentz-invarianten symmetrischen Formen auf so(3,1) sind zweidimensional (eps und die Spurform, Dimension "
      f"{len(inv_hh)}), auf so(3,2) dreidimensional ({len(inv_h)}: dazu eta auf m). Die Fisher-Form (T1) auf h hat den "
      f"eps-Anteil {c_eps:.1f} und den Spurform-Anteil {c_holst:.1f} (Holst-artig), mit relativem Rest {rest_h:.4f} "
      f"(auf ganz so(3,2): {rest_g:.4f}); sie ist unter K invariant ({gr(k_def)}), unter einem Boost e^(0,5 B) aber nicht "
      f"(relativer Defekt {lor_def:.4f}). Der Träger der MM-Form (Paare Drehung x Boost) und der Träger der Fisher-Form "
      f"(Boost x Boost, Q x Q) schneiden sich nicht. Für T2 gilt dasselbe (Wert 2 nur auf Boost x Boost)",
      len(inv_hh) == 2 and len(inv_h) == 3 and eps_inv < 1e-12 and abs(c_eps) < 1e-12 and abs(c_holst - 0.5) < 1e-9
      and rest_h > 0.5 and lor_def > 0.1 and k_def < 1e-9)


# --- C2 die native Familie ist auf der Gravitationsseite die triviale Störung
def g_schnitt(xv):
    return [sla.expm(xv[0]*E_T), sla.expm(xv[1]*Q_T[0]), sla.expm(xv[2]*Q_T[1]), sla.expm(xv[3]*Q_T[2])]


def a_form(xv):
    """A_mu = g^-1 d_mu g für g = e^(x0 E) e^(x1 Q1) e^(x2 Q2) e^(x3 Q3), exakt."""
    gs = g_schnitt(xv)
    gens_ = [E_T, Q_T[0], Q_T[1], Q_T[2]]
    out = []
    for mu in range(4):
        rest = np.eye(4, dtype=complex)
        for g_ in gs[mu + 1:]:
            rest = rest @ g_
        out.append(np.linalg.inv(rest) @ gens_[mu] @ rest)
    return out


def g_voll(xv):
    g = np.eye(4, dtype=complex)
    for g_ in g_schnitt(xv):
        g = g @ g_
    return g


X0 = np.array([0.3, 0.2, -0.4, 0.5])
hh = 1e-5
A_0 = a_form(X0)
dA = []
for mu in range(4):
    e_ = np.zeros(4)
    e_[mu] = hh
    Ap, Am = a_form(X0 + e_), a_form(X0 - e_)
    dA.append([(Ap[nu] - Am[nu])/(2*hh) for nu in range(4)])
F_mc = max(np.abs(dA[mu][nu] - dA[nu][mu] + kl(A_0[mu], A_0[nu])).max() for mu in range(4) for nu in range(4))
dA_sym = 0.0
for X in SO32:
    lam = lambda xv, X=X: np.linalg.inv(g_voll(xv)) @ X @ g_voll(xv)
    for mu in range(4):
        e_ = np.zeros(4)
        e_[mu] = hh
        dlam = (lam(X0 + e_) - lam(X0 - e_))/(2*hh)
        dA_sym = max(dA_sym, float(np.abs(dlam + kl(A_0[mu], lam(X0))).max()))
ZAHLEN["mm_vakuum"] = {"|F| der Maurer-Cartan-Form": float(F_mc), "max |delta A| über 10 Erzeuger": dA_sym,
                       "Rang der quadratischen MM-Wirkung auf der nativen Familie": 0, "Rang der Fisher-Form (T1, T2)": [6, 3]}
check(f"[E, negativ] Die native Familie ist auf der Gravitationsseite die triviale Störung: Für das native Vakuum A = "
      f"g^-1 dg (F = 0, Rest {gr(F_mc)}) lässt die globale Wirkung g -> g_0 g die Verbindung fest; die zugehörige "
      f"Eichstörung delta A = d lambda + [A, lambda] mit lambda = g^-1 X g verschwindet für alle zehn Erzeuger (max "
      f"{gr(dA_sym)}). Auch für ortsabhängiges lambda ist delta F = [F, lambda] = 0. Die quadratische MM-Wirkung und die "
      f"kanonische Energie der linearisierten Gravitation haben auf der nativen Deformationsfamilie den Rang 0, die "
      f"Fisher-Form den Rang 6 (T1) bzw. 3 (T2). 'Fisher = kanonische Energie der Metrikstörung' ist dort keine "
      f"Gleichung für G, sondern verletzt (0 gegen N_F)", F_mc < 1e-6 and dA_sym < 1e-6)

# --- C3 Hollands-Wald auf der Bahn: Masse statt Kopplung; kein Horizontterm
l_, rp, Gs = sp.symbols("l r_p G", positive=True)
T_bh = (l_**2 + 3*rp**2)/(4*sp.pi*l_**2*rp)                          # Schwarzschild-AdS4
r_min = sp.solve(sp.diff(T_bh, rp), rp)[0]
T_min = sp.simplify(T_bh.subs(rp, r_min))
T_nat = 1/(4*sp.pi*l_)                                               # beta = 4 pi in Einheiten l (E auf Vektoren normiert)
T_hp = 1/(sp.pi*l_)
verh = sp.simplify(T_min/T_nat)
S_bh = sp.pi*rp**2/Gs
C_bh = sp.simplify(T_bh*sp.diff(S_bh, rp)/sp.diff(T_bh, rp))
ZAHLEN["horizont"] = {"T_min l": str(sp.simplify(T_min*l_)), "T_nativ l": str(sp.simplify(T_nat*l_)),
                      "T_HP l": str(sp.simplify(T_hp*l_)), "T_min/T_nativ": str(verh),
                      "Wärmekapazität des Lochs": str(sp.factor(C_bh)), "M l (Leiter)": mp.nstr(E_mp, 15)}
check(f"[E/C] Hollands-Wald auf der Bahn (Lesart P_HW: der Zustand als Körper der Masse M in AdS4): Die kanonische "
      f"Energie einer Symmetrieverschiebung ist delta^2 M, also M l (s^2/2) je Einheitsrichtung; Fisher = beta x "
      f"kanonische Energie gilt dann identisch mit "
      f"M l = <E> = {mp.nstr(E_mp, 10)} (Leiter). Das ist der erste Hauptsatz von z136 in zweiter Ordnung und enthält G "
      f"nicht: G erschiene erst im Schwarzschild-Parameter 2 G M/l = 2 <E> G/l^2, den nichts festlegt. Der einzige "
      f"G-tragende Term der kanonischen Energie, delta^2 A/(4G) eines Horizonts, fehlt nativ: Ein AdS4-Loch hat T >= "
      f"T_min = sqrt(3)/(2 pi l) (bei r_+ = l/sqrt 3), der native Zustand T = 1/(4 pi l), um den Faktor 2 sqrt(3) = "
      f"{float(verh):.4f} darunter (Hawking-Page 1/(pi l)), unabhängig von G. Ein Loch hätte die Fisher-Zahl Var(K) = C = "
      f"2 pi r_+^2 (3 r_+^2 + l^2)/(G (3 r_+^2 - l^2)), proportional zu l^2/G; bei beta = 4 pi l existiert es für kein G "
      f"[C: Einstein-Gleichung als Vergleich, nicht nativ]",
      sp.simplify(T_min - sp.sqrt(3)/(2*sp.pi*l_)) == 0 and sp.simplify(verh - 2*sp.sqrt(3)) == 0
      and sp.simplify(C_bh - 2*sp.pi*rp**2*(l_**2 + 3*rp**2)/(Gs*(3*rp**2 - l_**2))) == 0
      and sp.simplify(r_min - l_/sp.sqrt(3)) == 0 and sp.simplify(T_hp/T_nat - 4) == 0)


# --- C4 Zentralterme: 1D (Virasoro) und 4D (C_T)
def n_q(x):
    return 1/(1 + mp.exp(4*mp.pi*x))


D_m = {}
for m in range(1, 6):
    D_m[m] = mp.fsum((mp.mpf(k) + mp.mpf(1)/4 + mp.mpf(m)/2)**2*(n_q(mp.mpf(k) + mp.mpf(1)/4)
                                                                 - n_q(mp.mpf(k) + mp.mpf(1)/4 + m))
                     for k in range(-80, 80))
L0_q = D_m[1]/2
F_vir = {m: 4*mp.pi*mp.mpf(m)/2*D_m[m] for m in D_m}
z136 = json.loads((HERE / "z136_ergebnis.json").read_text())
moeb_z136 = float(z136["zahlen"]["zweite_ordnung"]["Möbius F/2"])
vir_ok = all(abs(D_m[m] - mp.mpf(m)**3/12) < mp.mpf(10)**-30 for m in D_m)
ZAHLEN["virasoro_1D"] = {"D_m = <[L_m, L_-m]> (Quarter, beta = 4 pi)": {str(m): mp.nstr(D_m[m], 20) for m in D_m},
                         "<L_0> mit Twistgewicht": mp.nstr(L0_q, 25),
                         "F_m = beta (m/2) D_m": {str(m): mp.nstr(F_vir[m], 15) for m in F_vir},
                         "pi m^4/6": {str(m): mp.nstr(mp.pi*m**4/6, 15) for m in F_vir},
                         "z136 Möbius F/2": moeb_z136, "pi/12": mp.nstr(mp.pi/12, 12)}
check(f"[E, 1D-typisiert] Der Quarter-Fockraum der Naht unter der radialen sl(2,R) (z136) ist dieselbe Rechnung: "
      f"F = beta <L_0> mit dem kommutatorfesten L_0 = (1/2)[L_1, L_-1] (Twistgewicht 1/32 eingeschlossen). Bei beta = "
      f"4 pi ist <L_0> = 1/32 + 1/96 = 1/24 = c/24 exakt ({mp.nstr(L0_q - mp.mpf(1)/24, 3)}; Bestand: h = 1/32 aus z59, "
      f"1/96 aus z68, Spannungsfreiheit aus z107), also N_F = beta c/24 = pi/6 je Quarter; z136 fand F/2 = "
      f"{moeb_z136:.8f} = pi/12 = {mp.nstr(mp.pi/12, 10)} ohne die Zahl zu benennen. Über die native Algebra hinaus (Virasoro, "
      f"m >= 2) trägt die Form den Zentralterm: D_m = <[L_m, L_-m]> = 2 m <L_0> + (c/12)(m^3 - m) = (c/12) m^3 für m = "
      f"1..5 (Abweichung < 1e-30), F_m = pi c m^4/6 mit c = 1. Die Zentralladung erscheint nur außerhalb der nativen "
      f"sl(2,R)", vir_ok and abs(L0_q - mp.mpf(1)/24) < mp.mpf(10)**-30 and abs(moeb_z136 - float(mp.pi/12)) < 1e-7)

LEV6, ZHAT6 = leiter(6)
H6 = 1j*ZHAT6(E_T)
senk = [ZHAT6(Q_T[k] - 1j*B_T[k]) for k in range(3)]
senk_ok = max(np.abs(kl(H6, L_) + L_).max() for L_ in senk)
_, s3, vh3 = np.linalg.svd(np.vstack(senk))
kern3 = vh3[int(np.sum(s3 > 1e-9)):]
stufen3 = sorted(round(float(v), 9) for v in np.linalg.eigvalsh(kern3.conj() @ np.diag(LEV6 + 1.0) @ kern3.T))
rot_kern = max(float(np.abs(ZHAT6(R_) @ kern3.T).max()) for R_ in R_T)        # Spin 0: von den Drehungen annulliert
senk4 = senk + [ZHAT6(QM[0] - 1j*D_T)]
senk4_ok = np.abs(kl(H6, senk4[3]) + senk4[3]).max()
n_lw4 = int(np.sum(np.linalg.svd(np.vstack(senk4), compute_uv=False) < 1e-9))
d_s = sp.symbols("d", positive=True)
S_d = 2*sp.pi**(d_s/2)/sp.gamma(d_s/2)
CT_frei = d_s/((d_s - 1)*S_d**2)
CT4, CT3 = sp.simplify(CT_frei.subs(d_s, 4)), sp.simplify(CT_frei.subs(d_s, 3))
CT_holo = (d_s + 1)/(d_s - 1)*sp.gamma(d_s + 1)/(sp.pi**(d_s/2)*sp.gamma(d_s/2))/(8*sp.pi)     # mal L^(d-1)/G
holo3, holo4 = sp.simplify(CT_holo.subs(d_s, 3)), sp.simplify(CT_holo.subs(d_s, 4))
ZAHLEN["randinhalt"] = {"so(3,2)-Tiefstgewichte in der Leiter: Zahl, Stufen": [len(kern3), stufen3],
                        "su(2,2)-Tiefstgewichte": n_lw4, "C_T freier Skalar d = 4, d = 3": [str(CT4), str(CT3)],
                        "C_T G/L^(d-1) holographisch d = 3, d = 4 (Literatur)": [str(holo3), str(holo4)],
                        "L^3/G_5 je Leiter (P_holo5)": str(sp.simplify(CT4/holo4)),
                        "l^2/G_4 je 3D-Skalareinheit (P_holo3)": str(sp.simplify(CT3/holo3))}
check(f"[E/L] 4D: Die native Algebra ist halbeinfach, ihre Darstellung auf der Leiter hat keinen Zentralterm (A1); "
      f"deshalb enthält die Fisher-Form auf der nativen Bahn nur beta und <E>. Das 4D-Gegenstück der Zentralladung ist "
      f"C_T = {CT4} je Leiter (z137; Formel d/((d - 1) S_d^2), in d = 3: {CT3}); es erscheint erst in der BKM-Form von "
      f"Metrikquellen, außerhalb der nativen Familie, und dort mit dem Abschneidelogarithmus von z137. Ein Wörterbuch "
      f"Lashkari-Van Raamsdonk verlangte zusätzlich Holographie [L, Literaturkoeffizienten, hier nicht nachgerechnet]: "
      f"C_T = {holo4} L^3/G_5 (AdS5, falsche Dimension für MM) bzw. C_T = {holo3} l^2/G_4 für eine 3D-Randtheorie von "
      f"AdS4. Deren Inhalt ist nativ aber kein lokales Feld mit Spannungstensor: Unter so(3,2) zerfällt die Leiter in "
      f"genau {len(kern3)} Tiefstgewichtsdarstellungen mit E_0 = {stufen3} und Spin 0 (Drehungen: {gr(rot_kern)}; D(1,0) + "
      f"D(2,0), verallgemeinerte "
      f"freie Felder; unter su(2,2): {n_lw4}); ein erhaltener Spin-2-Strom mit Delta = 3 fehlt auf der Einteilchenleiter",
      hom < 1e-12 and CT4 == 1/(3*sp.pi**4) and CT3 == sp.Rational(3, 32)/sp.pi**2 and len(kern3) == 2
      and stufen3 == [1.0, 2.0] and rot_kern < 1e-9 and n_lw4 == 1 and senk_ok < 1e-9 and senk4_ok < 1e-9
      and sp.simplify(holo3 - 3/sp.pi**3) == 0 and sp.simplify(holo4 - 5/sp.pi**3) == 0)

# --- C5 neue Zeilen der Überbestimmungstabelle (reine Zahl G/R_Vak^2), Vergleich mit z136


def g_aus_w(w):
    """z118/z127/z136: w c_1 lambda = 1/(16 pi G), c_1 = 64, lambda = 4/R_Vak^2, also G/R_Vak^2 = 1/(4096 pi w)."""
    return 1/(4096*mp.pi*w)


NEU = [
    ("(f1) Hollands-Wald auf der Bahn: Fisher = beta delta^2 M, M l = <E>", None, "E",
     "G tritt nicht auf; nur 2 G M/l = 2 <E> G/l^2"),
    ("(f2) P_FK, T1, P_1T: Spurgewicht w = N_F = beta <E> (Gruppenparameter)", g_aus_w(NF_mp), "L",
     "w = 4 pi (1 + 4 x + ...)"),
    ("(f3) P_FK, T1, P_1T: w = <E>/beta (Erzeugerparameter, Kubo-Mori-Suszeptibilität)", g_aus_w(SUSZ_mp), "L",
     "= (1/1024)/<E>; z118 (a2) mit 2 c_3 mal 1/<E>"),
    ("(f4) P_FK, T1, P_Feld: w = beta <E>_Fock", g_aus_w(NF_FOCK_mp), "L", "<E>_Fock = (E_4(2i) - 1)/240"),
    ("(f5) P_FK, T1, P_Feld: w = <E>_Fock/beta", g_aus_w(E_FOCK_mp/(4*mp.pi)), "L", "Erzeugerparameter"),
    ("(f6) P_FK, T2: w = 2 (Boosts und D, q-unabhängig)", g_aus_w(mp.mpf(2)), "L", "z115"),
    ("(f7) P_FK, T1, Quarter (1D): w = beta c/24 = pi/6", g_aus_w(mp.pi/6), "L", "radiale sl(2,R); je Naht pi/3"),
    ("(f8) P_etaF: G/l^2 = 1/(c_eta N_F), c_eta frei", None, "C/O",
     "Fisher-Metrik ist dimensionslos; braucht 'Fisher-Information je Planck-Fläche' (wie P_eta aus z136)"),
    ("(f9) P_holo3 (Lashkari-Van Raamsdonk mit Flächenregel): G/R_Vak^2 = 8/(pi n)", None, "L (Import)",
     "n = C_T der 3D-Randtheorie in Skalareinheiten 3/(32 pi^2); nativ kein Randspannungstensor"),
    ("(f10) P_BH3 (Brown-Henneaux, 1D-typisiert): l_3/G_3 = 2 c/3", None, "L (Import)",
     "3D-Kopplung; G_4 = G_3 L_quer braucht eine freie Querlänge"),
]
ALT = [(z["Quelle"], z["Lesart"], mp.mpf(z["G/R_Vak^2"])) for z in z136["zahlen"]["tabelle"] if z["G/R_Vak^2"] is not None]
say("   Neue Zeilen der Überbestimmungstabelle (reine Zahl G/R_Vak^2; '--' = keine reine Zahl):")
tab_json = []
naechste = []
for lesart, wert, status, bem in NEU:
    ws = "--" if wert is None else f"{mp.nstr(wert, 9)} (10^{float(mp.log10(wert)):+.2f})"
    say(f"     (f) z146 | {lesart:<88} | {ws:<28} | {status:<10} | {bem}")
    tab_json.append({"Quelle": "(f) z146", "Lesart": lesart, "G/R_Vak^2": None if wert is None else mp.nstr(wert, 15),
                     "Status": status, "Bemerkung": bem})
    if wert is not None:
        fak, q_, l_alt = min(((max(wert/a, a/wert), q_, l_alt) for q_, l_alt, a in ALT), key=lambda t: t[0])
        naechste.append((lesart, float(fak), f"{q_}: {l_alt}"))
exakt = [t for t in naechste if abs(t[1] - 1) < 1e-9]
beste = min(naechste, key=lambda t: t[1])
werte_neu = sorted(float(w) for _, w, _, _ in NEU if w is not None)
intern = min(max(a/b, b/a) for i, a in enumerate(werte_neu) for b in werte_neu[i + 1:])
ZAHLEN["tabelle"] = tab_json
ZAHLEN["vergleich"] = {"alte numerische Zeilen (z136)": len(ALT), "exakte Übereinstimmungen": len(exakt),
                       "nächste alte Zeile je neuer Zeile (Faktor)": [[a, round(f, 9), b] for a, f, b in naechste],
                       "kleinster Faktor zwischen zwei neuen Zeilen": round(intern, 6),
                       "Spanne der neuen Zeilen (log10)": [round(math.log10(werte_neu[0]), 3), round(math.log10(werte_neu[-1]), 3)]}
check(f"[E] Vergleich mit der Tabelle von z136 ({len(ALT)} numerische Zeilen): keine der sechs neuen Zahlen stimmt mit "
      f"einer vorhandenen exakt überein ({len(exakt)} Treffer). Am nächsten liegt '{beste[0]}' an '{beste[2]}' mit dem "
      f"Faktor {beste[1]:.10f} = <E>_beta = 1 + x/(1 + x) + 3x/(1 - x): Die Kubo-Mori-Suszeptibilität der "
      f"Einteilchenleiter ist 2 c_3 mal <E>, das von Hand gesetzte Gewicht 2 c_3 = 1/beta aus z118 also bis auf die "
      f"thermische Besetzung 4 e^(-4 pi) = {mp.nstr(4*x_mp, 4)} der Wert Delta/beta des Leitergrundzustands. Das ist "
      f"kein Treffer: Es braucht drei Wahlen (P_FK, Erzeuger- statt Gruppenparameter mit Faktor beta^2 = 16 pi^2, "
      f"Einteilchenleiter statt Feld mit Faktor {mp.nstr(E_mp/E_FOCK_mp, 4)}), ist nicht exakt, und die Fisher-Form ist "
      f"nach C1 nicht die MM-Form. Die sechs Fisher-Zeilen selbst streuen über 10^{math.log10(werte_neu[0]):.2f} bis "
      f"10^{math.log10(werte_neu[-1]):.2f}, keine zwei stimmen überein (kleinster Faktor {intern:.4f})",
      len(exakt) == 0 and abs(beste[1] - float(E_mp)) < 1e-9 and "(f3)" in beste[0] and "2 c_3" in beste[2]
      and intern > 1.5)
befund("Zuordnung Fisher = kanonische Energie ohne freie Wahl", "nein",
       "E/C: Auf der nativen Familie ist die Gravitationsseite null (reine Eichrichtung, kein Horizont bei beta = 4 pi); "
       "die Fisher-Form ist die Bewegungs- und Lageenergie einer Masse M l = <E> und enthält G nicht. Eine Zahl für "
       "l^2/G entsteht nur mit einer zusätzlichen Gleichsetzung (P_FK, P_etaF) oder einem holographischen Import "
       "(P_holo3, P_BH3); jede gibt eine weitere Zeile, keine stimmt mit einer vorhandenen exakt überein.")

# ====================================================================================================================
say("\nTeil D: Nichtzirkularität (GRAV.NONCIRCULAR.01)")
EINGABEN = [
    ("beta = 4 pi, Modulgenerator K = 4 pi E (elliptisch)", "nativ (A*, z107; z136)", "kein BW/CHM, kein geometrischer Modulfluss behauptet"),
    ("Normierung von E (Periode 2 pi auf Vektoren), Einheit in p = Rapidität", "nativ (z136, hier A1)", "legt N_F als reine Zahl fest"),
    ("Zustandsraum: Einteilchenleiter, Leiterfeld, Quarter-Fockraum; T1 gegen T2", "Wahl (P_1T, P_Feld, T1/T2)", "vier verschiedene Zahlen N_F"),
    ("l = R_Vak/2, Lambda = -12/R_Vak^2", "nativ (z98, z118)", "nur für die Umrechnung in R_Vak"),
    ("b = l^2/(64 pi G) bzw. w c_1 lambda = 1/(16 pi G)", "Definition von G als MM-Koeffizient (z98, z118)", "kein G = 1; G bleibt Symbol"),
    ("P_FK: Fisher-Koeffizient = Spurgewicht w der Uhrspur", "Prämisse, nicht hergeleitet", "Zeilen (f2) bis (f7)"),
    ("Schwarzschild-AdS-Thermodynamik (T_min, C)", "Vergleich, Einstein-Gleichung vorausgesetzt", "nur zum Ausschluss des Horizontterms"),
    ("Holographische Koeffizienten (C_T ~ l^2/G, c = 3 l/2G)", "IMPORT: Flächenregel/RT bzw. Brown-Henneaux", "Zeilen (f9), (f10); nach GRAV.NONCIRCULAR.01 ASSUMPTION"),
]
for was, herkunft, rolle in EINGABEN:
    say(f"   {was}\n      Herkunft: {herkunft}\n      Rolle: {rolle}")
ZAHLEN["eingaben"] = [{"Eingabe": a, "Herkunft": b, "Rolle": c} for a, b, c in EINGABEN]
g_frei = all(z["G/R_Vak^2"] is None for z in tab_json if z["Status"] in ("E", "C/O", "L (Import)"))
check("[E] Nichtzirkularität: In die Fisher-Formen (Teile A, B) gehen weder G = 1 noch eine Flächenregel noch "
      "Bisognano-Wichmann ein; sie hängen nur an beta = 4 pi, der Normierung der Uhr und dem gewählten Zustandsraum. "
      "Eben deshalb enthalten sie G nicht. G erscheint erst (i) als Definition des MM-Koeffizienten (kein G = 1) "
      "zusammen mit der Prämisse P_FK oder (ii) über einen holographischen Import, der die Flächenregel S = A/(4G) "
      "voraussetzt (Lashkari-Van Raamsdonk und Hollands-Wald benutzen den Horizontterm delta^2 A/(4G)); (ii) wäre nach "
      "GRAV.NONCIRCULAR.01 eine ASSUMPTION und zirkulär als Herleitung von G. Die Zeilen ohne Prämisse und die "
      "Importzeilen tragen keine Zahl", g_frei and len(EINGABEN) == 8)

# ====================================================================================================================
say("\nTeil E: Bilanz und Urteil")
bilanz = [
    ("E", "T1, Relativentropie auf der Gruppe",
     "S(rho_g || rho) = beta <E> (tr(g^+ g)/4 - 1) für g in SU(2,2); Hesse-Matrix beta <E> tr_T([E, X][E, Y])."),
    ("E", "T1, Zerlegung",
     f"Kern = Stabilisator von E; auf p positiv definit, Rang 8 (su(2,2)) bzw. 6 (so(3,2)), ein Eigenwert N_F = beta <E> "
     f"= {mp.nstr(NF_mp, 12)} (Einteilchenleiter), Killing-proportional auf p (Schur), Bergman-Metrik von G/K."),
    ("E", "weitere native Normierungen",
     f"Leiterfeld N_F = beta (E_4(2i) - 1)/240 = {mp.nstr(NF_FOCK_mp, 8)}; Quarter (1D) N_F = beta c/24 = pi/6; "
     f"Erzeugerparameter <E>/beta = 2 c_3 <E>; Temperaturrichtung Var(K) = {mp.nstr(VARK_mp, 6)}."),
    ("E", "T2 (Twistorlesart)", "Rang 3 (so(3,2)) bzw. 4 (su(2,2)), Wert 2 auf Boosts und D, q-unabhängig; nicht "
     "Killing-proportional; quadratischer Bereich nur |s| << e^(-2 pi)."),
    ("E negativ", "(a) MM-/Einstein-Form",
     "Nein: falscher Raum (G/K statt G/Lorentz), auf m entartet (0, +, +, +), auf den Translationen euklidisch, nicht "
     "Lorentz-invariant, eps-Anteil 0; die Gravitationsseite ist auf der nativen Familie null."),
    ("E/C", "Hollands-Wald auf der Bahn",
     "Fisher = beta x Bewegungs- und Lageenergie einer Masse M l = <E>; kein G. Horizontterm fehlt: T_nativ = "
     "T_min/(2 sqrt 3) für jedes G."),
    ("E, 1D", "Zentralterm", "Quarter bei beta = 4 pi: <L_0> = c/24 exakt, <[L_m, L_-m]> = (c/12) m^3, F_m = pi c m^4/6; "
     "c erscheint nur außerhalb der nativen sl(2,R). In 4D entspricht dem C_T (z137), ebenfalls außerhalb."),
    ("E", "(b) Zahl für G/R_Vak^2",
     f"Keine ohne Zusatzprämisse. Unter P_FK sechs weitere Zeilen von 10^{math.log10(werte_neu[0]):.2f} bis "
     f"10^{math.log10(werte_neu[-1]):.2f}; keine exakte Übereinstimmung; nächste: (f3) = z118 (2 c_3) geteilt durch <E>."),
    ("E", "Nichtzirkularität", "Kein G = 1, keine Flächenregel in den Fisher-Formen; Importe benannt."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
urteil = (
    "Gescopt negativ: Die Relativentropie mit den nativen Tomita-Daten legt P_G nicht fest; sie gibt nur weitere Zeilen "
    "der Überbestimmungstabelle. Nichts ist befördert. (1) Rechnung: In der Thermofeldlesart T1 (die native Gruppe wirkt "
    "durch Automorphismen der Algebra) ist die relative Entropie auf der ganzen Bahn exakt S(rho_g || rho) = beta <E> "
    "(tr(g^+ g)/4 - 1). Ihre Hesse-Matrix auf allen 15 Erzeugern der su(2,2) ist beta <E> tr_T([E, X][E, Y]): null auf "
    "dem Stabilisator des Zustands (E, Drehungen, Q^-_k), positiv definit mit einem einzigen Eigenwert N_F = beta <E> "
    "auf den acht nichtkompakten Richtungen (so(3,2): sechs, Q_k und Boosts). Auf dem Nebenraum ist sie proportional "
    "zur Killing-Form, erzwungen durch Schur; sie ist die Bergman-Metrik des hermitesch-symmetrischen Raums G/K, für "
    "so(3,2) des Raums der Uhrweltlinien von AdS4. In der Twistorlesart T2 (z115) hat die Form den Rang 3 bzw. 4 mit "
    "dem festen Wert 2 auf Boosts und D. (2) Zu (a): Die Form ist nicht die linearisierte MacDowell-Mansouri- oder "
    "Einstein-Form. Sie lebt auf G/K statt auf AdS4 = G/Lorentz, ist auf der AdS-Tangente entartet (0, +, +, +), auf "
    "den Minkowski-Translationen euklidisch, unter Boosts nicht invariant und hat den MM-Anteil (eps) null. Umgekehrt "
    "ist die native Familie auf der Gravitationsseite die triviale Störung: Die Maurer-Cartan-Form ist unter der "
    "globalen Gruppe invariant, die quadratische MM-Wirkung und die kanonische Energie der Metrikstörung haben dort den "
    "Rang 0. (3) Zu (b): Eine eindeutige Zahl gibt es je Lesart, aber sie ist G-frei: N_F = beta <E> = 4 pi (1 + 4 "
    "e^(-4 pi) + ...) auf der Einteilchenleiter, beta (E_4(2i) - 1)/240 im Leiterfeld, beta c/24 = pi/6 je Quarter in "
    "1D, 2 in T2. Physikalisch ist sie beta mal die Bewegungs- und Lageenergie einer Masse M l = <E> im AdS-Potential, "
    "also Hollands-Wald für Materie, der erste Hauptsatz von z136 in zweiter Ordnung. G käme nur über einen "
    "Horizontterm (fehlt: der native Zustand liegt um 2 sqrt 3 unter der Mindesttemperatur eines AdS4-Lochs, für jedes "
    "G) oder über einen Zentralterm. Den trägt die halbeinfache native Algebra nicht; er erscheint erst außerhalb, in "
    "1D als Virasoro-Ladung (bei beta = 4 pi exakt <[L_m, L_-m]> = (c/12) m^3), in 4D als C_T von z137, und braucht "
    "dann ein holographisches Wörterbuch, das die Flächenregel importiert. (4) Tabelle: Unter der Prämisse P_FK "
    f"(Fisher-Koeffizient = Spurgewicht der Uhrspur) entstehen sechs Zeilen zwischen 10^{math.log10(werte_neu[0]):.2f} und "
    f"10^{math.log10(werte_neu[-1]):.2f}; keine stimmt "
    "mit einer vorhandenen exakt überein. Die nächste ist (1/1024)/<E>: Die Kubo-Mori-Suszeptibilität der "
    "Einteilchenleiter ist 2 c_3 <E>, das Gewicht 2 c_3 = 1/beta aus z118 also bis auf den Faktor 1 + 4 e^(-4 pi) der "
    "Wert Delta/beta des Grundzustands. Das gibt dem gesetzten Gewicht einen nativen Namen, ist aber kein Treffer "
    "(drei Wahlen, nicht exakt, falsche Invarianz). (5) Nichtzirkularität: In die Fisher-Formen geht weder G = 1 noch "
    "eine Flächenregel ein; eben deshalb enthalten sie G nicht. Folgen: K13 bleibt L, P_G bleibt eine freie Zahl, "
    "P_curv ist nicht berührt (die Krümmung von G/K ist die des Zustandsraums, die Bahn ist Maurer-Cartan-flach). K15 "
    "und die übrigen Prämissen (A*, I_step, P_lift, P_IF, P_univ, P_uhr, P_ord) sind nicht berührt. Geltungsbereich des "
    "Negativs: die Bahn der nativen Gruppe und der Kegel ihrer Modulgeneratoren, für quasifreie Zustände in T1 "
    "(Einteilchenleiter, Leiterfeld, Quarter) und T2. Nicht erfasst sind Metrikquellen in 4D (z137) und ein "
    "wechselwirkendes Netz.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z146_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "bilanz": [list(b) for b in bilanz], "urteil": urteil},
    indent=1, ensure_ascii=False, default=str) + "\n")
(HERE / "ausgabe_z146.txt").write_text("\n".join(OUT) + "\n")
