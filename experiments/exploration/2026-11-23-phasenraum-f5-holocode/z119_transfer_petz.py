"""Z119: Folgeschritt U3 der sechzehnten Analyse (25.11.): Austauschwechselwirkung aus virtuellen Ladungstransfers und
Petz-Rekonstruierbarkeit statt P_koh.

Leitlinie 'Rekonstruktion statt Auswahl': Eine Mittelung oder ein Extremalprinzip zählt nur, wenn das Original es liefert.
Firewall DIMENSION.UPLIFT.FIREWALL.01: alle Aussagen typisiert (Nahtmodell mit zwei Moden je Naht, Kanal zwischen
Paarschalen H_m -> H_(m+1)); keine wechselwirkende 3+1-Schließung. Gelöste Bausteine bleiben unberührt. Bestand: z116 Teil 8
(Transfer-SU(2), Pauli-Kopplung Spin-Singulett = Isospin 1, w_Singulett = (1 + 2 e^-4pi)/3) und Teil 9 (F_koh = (m+1)/(2m),
F_deko = (2m+1)/(3 m^2)), z99 (su(2)_1-Transferströme), z106 (gemeinsame CAR, Gamma(pi) in der U(2) der Ströme, Kette nicht
beobachtbar, P_ident), z91 (Phi_1/2, P_koh, Torus, Spiegel, L^nat-Casimir), z80 (thermische Gewichte W_m(r), Raten u_m(r)),
z83 (L^nat = Summe (k+1) psi*_(k+1) psi_k mit den Amplituden n_a + 1), z105 (phi^3 über drei Nähte), z108 (Kompositionssatz,
I_AB|C).

Teil 1  Austausch. (a) KMS gegen Twirl: Der KMS-Zustand jeder additiven Dynamik ist auf festen Ladungen ein Produkt (bei
        reiner Ladungsenergie maximal gemischt); diskrete Eichungen (Kette W_12, Tausch Gamma(pi)) lassen ihn fest; die
        kontinuierliche Transfer-SU(2) bricht die Ladungsenergie, der Twirl ändert beobachtbare Nahtladungen; auf
        transferinvarianten Observablen sind KMS und Twirl gleich, und log 3 ist der Vielfachheitsterm log dim(Isospin 1).
        Vergleich: Hüpfdynamik (nicht nativ) gibt das entgegengesetzte Vorzeichen. (b) N = 3, 4 Nähte: Schur-Gewichte,
        K_H(S) = log(dim W/f) + O(e^-4pi), Paulizerlegung, Paarkopplung J_N (N = 2..12), I_AB|C, Vergleich z105/z108.
        (c) Statistik des Dubletts bestimmt das Vorzeichen; Gamma(pi) und P_ident; keine 1->2-Amplitude.
Teil 2  Petz. (a) Satz im torusverträglichen Kanalraum mit nativen Raten: F = Tr M^2 <= F_koh, Gleichheit genau für den
        kohärenten Kanal (bis auf diagonale Eichung), für jede treue diagonale Referenz (gleichförmig, KMS-Gewichte z80 bei
        beta = 4 pi und 2, zufällig). (b) Ohne Torus: exakte Blockkanäle mit Raten, Spurerhaltung, Zustandserhaltung,
        J-Realität, Petz-Symmetrie und Spiegel übertreffen F_koh; nur mit Spurerhaltung bis 0,88 gegen 0,75. (c) Native
        Dilatation: L^nat als Operator (Rang 1, zwei Torusgewichte, nicht spurerhaltend), gewichtsaufgelöste Polarfamilie
        (native Raten genau bei t = 1/2, dann Phi_1/2), diagonaler Petz-Komposit = nativer Casimir für jeden Toruskanal.
Teil 3  Bilanz (E/C/L/O) und Urteil.
"""
import itertools
import json
import time
from fractions import Fraction as Fr
from math import comb
from pathlib import Path

import numpy as np
import scipy.linalg as sla
from scipy.optimize import minimize

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
rng = np.random.default_rng(119)
BETA = 4*np.pi
QQ = np.exp(-BETA)                                                   # e^-4pi = q^2


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def jw(n):
    I_, Zp, smm = np.eye(2), np.diag([1.0, -1.0]), np.array([[0.0, 1.0], [0.0, 0.0]])
    out = []
    for j in range(n):
        m_ = [Zp]*j + [smm] + [I_]*(n - j - 1)
        op = m_[0]
        for Mx in m_[1:]:
            op = np.kron(op, Mx)
        out.append(op)
    return out


def kr(*m):
    out = m[0]
    for x in m[1:]:
        out = np.kron(out, x)
    return out


def msqrt(Mx, inv=False):
    w, U = np.linalg.eigh((Mx + Mx.conj().T)/2)
    w = np.clip(w, 0, None)
    d = np.array([1/np.sqrt(x) if x > 1e-13*max(w.max(), 1e-300) else 0.0 for x in w]) if inv else np.sqrt(w)
    return (U*d) @ U.conj().T


def mlog(Mx):
    w, U = np.linalg.eigh((Mx + Mx.conj().T)/2)
    return (U*np.log(w)) @ U.conj().T


def entropie(r):
    w = np.linalg.eigvalsh((r + r.conj().T)/2)
    w = w[w > 1e-15]
    return float(-(w*np.log(w)).sum())


def ptrace(r, dims, keep):
    """Partielle Spur eines Zustands auf C^d1 x ... x C^dn, behalten werden die Faktoren in keep."""
    n = len(dims)
    t = r.reshape(dims + dims)
    weg = [i for i in range(n) if i not in keep]
    for i in sorted(weg, reverse=True):
        t = np.trace(t, axis1=i, axis2=i + t.ndim//2)
    d = int(np.prod([dims[i] for i in keep]))
    return t.reshape(d, d)


s0 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]], complex)
PAULI = [s0, sx, sy, sz]


# =================================================================== Teil 1
say("Teil 1: Austauschwechselwirkung aus virtuellen Ladungstransfers")


def nahtmodell(N):
    """N Nähte mit je zwei Moden (s = 0, 1) in einer gemeinsamen CAR (z106), Mode (i, s) -> 2 i + s."""
    c = jw(2*N)
    cd = [x.T for x in c]
    num = [x.T @ x for x in c]
    Ni = [num[2*i] + num[2*i + 1] for i in range(N)]
    dim = 2**(2*N)
    I_ = np.eye(dim)
    S = [sum((cd[2*i] @ c[2*i + 1] + cd[2*i + 1] @ c[2*i])/2 for i in range(N)),
         sum(-1j*(cd[2*i] @ c[2*i + 1] - cd[2*i + 1] @ c[2*i])/2 for i in range(N)),
         sum((num[2*i] - num[2*i + 1])/2 for i in range(N))]
    E = {(i, j): sum(cd[2*i + s] @ c[2*j + s] for s in range(2)) for i in range(N) for j in range(N)}
    Ntot = sum(Ni)
    idx = np.where(np.isclose(np.diag(Ntot), N))[0]                  # Sektor mit N Anregungen
    vak = np.zeros(dim)
    vak[0] = 1.0
    B = []                                                           # Basis |s_1..s_N> = c+_(1 s1) ... c+_(N sN) |0>
    for ss in itertools.product(range(2), repeat=N):
        v = vak.copy()
        for i in reversed(range(N)):
            v = cd[2*i + ss[i]] @ v
        B.append(v)
    B = np.array(B).T[idx, :]
    sub = lambda X: X[np.ix_(idx, idx)]
    return dict(N=N, c=c, cd=cd, Ni=Ni, Ntot=Ntot, S=S, E=E, idx=idx, B=B, sub=sub,
                S2=sub(sum(x @ x for x in S)).real, H=sub(sum((n_ - I_) @ (n_ - I_) for n_ in Ni)/2).real)


def spinbloecke(S2):
    ev, V = np.linalg.eigh(S2)
    out = {}
    for sv in sorted(set(np.round(ev, 8))):
        Vs = V[:, np.isclose(ev, sv)]
        out[float(sv)] = Vs @ Vs.conj().T
    return out


def twirl_schur(X, bl):
    """Mittelung über die Transfer-U(N): spininvariantes X -> Summe_S Tr(X P_S)/Tr(P_S) P_S (Howe-Dualität, vielfachfrei)."""
    return sum(np.trace(X @ P)/np.trace(P)*P for P in bl.values())


def haar_un(n):
    Zm = (rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n)))/np.sqrt(2)
    Q, Rr = np.linalg.qr(Zm)
    return Q @ np.diag(np.diag(Rr)/abs(np.diag(Rr)))


def fock_u(M, U):
    """Zweite Quantisierung Gamma(U) = exp(Summe (log U)_ij E_ij) auf dem Sektor."""
    Lg = sla.logm(U)
    G = sum(Lg[i, j]*M["sub"](M["E"][(i, j)]) for i in range(M["N"]) for j in range(M["N"]))
    return sla.expm(G)


# ---------------------------------------------------------------- (a)
say("  (a) Ist die Mittelung über die Transfer-SU(2) nativ?")
M2 = nahtmodell(2)
sub2 = M2["sub"]
bl2 = spinbloecke(M2["S2"])
D2 = np.diag(np.exp(-BETA*np.diag(M2["H"])))
rhoG = D2/np.trace(D2)
rhoT = twirl_schur(rhoG, bl2)
haar = np.zeros_like(rhoG, dtype=complex)
NH = 600
for _ in range(NH):
    G_ = fock_u(M2, haar_un(2))
    haar += G_ @ rhoG @ G_.conj().T/NH
B2 = M2["B"]
SZZ = kr(sz, sz)
S2phys = B2.T @ M2["S2"] @ B2
P_sing = spinbloecke(S2phys)[0.0]
P_trip = spinbloecke(S2phys)[2.0]


def auf_11(rho):
    r = B2.T @ rho @ B2
    return r/np.trace(r)


def kopplung(r):
    K = -mlog(r)
    return float(np.trace(K @ P_sing).real - np.trace(K @ P_trip).real/3)


rG11, rT11 = auf_11(rhoG), auf_11(rhoT)
kG, kT = kopplung(rG11), kopplung(rT11)
corrG, corrT = float(np.trace(rG11 @ SZZ).real), float(np.trace(rT11 @ SZZ).real)
td = 0.5*np.abs(np.linalg.eigvalsh(rG11 - rT11)).sum()
mi = lambda r: entropie(ptrace(r, [2, 2], [0])) + entropie(ptrace(r, [2, 2], [1])) - entropie(r)
# additive Dynamik: Ladungsenergie plus zufällige spinabhängige Nahtterme
hA, hB = [(lambda X: (X + X.conj().T)/8)(rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2))) for _ in range(2)]
cd2, c2 = M2["cd"], M2["c"]
Hadd = sum(hA[s, t]*cd2[s] @ c2[t] + hB[s, t]*cd2[2 + s] @ c2[2 + t] for s in range(2) for t in range(2)) \
    + 0.7*cd2[0] @ c2[0] @ cd2[1] @ c2[1] - 0.4*cd2[2] @ c2[2] @ cd2[3] @ c2[3]
Dadd = sla.expm(-BETA*(M2["H"] + sub2(Hadd)))
rA11 = auf_11(Dadd/np.trace(Dadd))
mi_add, mi_tw = mi(rA11), mi(rT11)
ZAHLEN["1a_kms_gegen_twirl"] = {"Kopplung KMS (1,1)": kG, "Kopplung Twirl (1,1)": kT, "<s^z s^z>*4 KMS": corrG,
                                "<s^z s^z>*4 Twirl": corrT, "Spurabstand": td, "I(A:B) additiv": mi_add,
                                "I(A:B) Twirl": mi_tw, "Haar-Abweichung (600)": float(np.abs(haar - rhoT).max())}
check(f"[E negativ] (a-i) KMS-Route: Der KMS-Zustand der Ladungsenergie Q^2/2 bei beta = 4 pi, eingeschränkt auf die festen "
      f"Ladungen (1,1), ist maximal gemischt (Abweichung von I/4 {np.abs(rG11 - np.eye(4)/4).max():.0e}): Austauschkopplung "
      f"{kG:.1e}, Korrelation <sigma^z_A sigma^z_B> = {corrG:.1e}. Der getwirlte Zustand (Schur, Haar-Stichprobe der echten "
      f"Transfer-SU(2), Abweichung {np.abs(haar - rhoT).max():.3f}) gibt auf (1,1) die Kopplung {kT:.7f} = log(3/(1 + 2 e^-4pi)) "
      f"und die Korrelation {corrT:.6f}; Spurabstand beider (1,1)-Zustände {td:.4f}. Mit zufälligen spinabhängigen "
      f"Nahttermen (additive Dynamik) bleibt der eingeschränkte KMS-Zustand ein Produkt: I(A:B) = {mi_add:.1e} gegen "
      f"{mi_tw:.4f} beim Twirl (z108: additive Komposition = Produkt). Der eingeschränkte KMS-Zustand ist nicht der getwirlte",
      np.abs(rG11 - np.eye(4)/4).max() < 1e-12 and abs(kG) < 1e-10 and abs(kT - np.log(3/(1 + 2*QQ))) < 1e-10
      and abs(corrT - 0.2) < 1e-5 and mi_add < 1e-10 and mi_tw > 0.01 and np.abs(haar - rhoT).max() < 0.03)

# diskrete Eichungen und kontinuierliche Transfer-SU(2)
NA, NB = sub2(M2["Ni"][0]), sub2(M2["Ni"][1])
W12 = np.diag((-1.0)**(np.diag(NA)*np.diag(NB)))
X01 = M2["E"][(0, 1)] + M2["E"][(1, 0)]
Gam = sla.expm(-1j*np.pi/2*M2["Ntot"]) @ sla.expm(1j*np.pi/2*X01)
fehl_gam = max(np.abs(Gam @ c2[s] @ Gam.conj().T - c2[2 + s]).max() for s in range(2)) \
    + max(np.abs(Gam @ c2[2 + s] @ Gam.conj().T - c2[s]).max() for s in range(2))
GamS = sub2(Gam)
kom_W, kom_G = np.abs(W12 @ rhoG - rhoG @ W12).max(), np.abs(GamS @ rhoG - rhoG @ GamS).max()
kG_W = kopplung(auf_11((rhoG + W12 @ rhoG @ W12)/2))
kG_G = kopplung(auf_11((rhoG + GamS @ rhoG @ GamS.conj().T)/2))
kom_HE = np.abs(M2["H"] @ sub2(M2["E"][(0, 1)]) - sub2(M2["E"][(0, 1)]) @ M2["H"]).max()
P1A = np.diag(np.isclose(np.diag(NA), 1).astype(float))
ion_G, ion_T = 1 - float(np.trace(rhoG @ P1A).real), 1 - float(np.trace(rhoT @ P1A).real)
ZAHLEN["1a_eichungen"] = {"Gamma(pi) tauscht Moden (Fehler)": float(fehl_gam), "[W_12, rho_KMS]": float(kom_W),
                          "[Gamma, rho_KMS]": float(kom_G), "||[H, E^+]||": float(kom_HE),
                          "P(N_A != 1) KMS": ion_G, "P(N_A != 1) Twirl": ion_T}
check(f"[E negativ] (a-ii) Eichlesart wie z106: Die Kette W_12 = (-1)^(N_A N_B) und der kanonische Tausch Gamma(pi) = "
      f"e^(-i pi N/2) e^(i pi/2 (E^+ + E^-)) (bildet c_A <-> c_B exakt ab, Fehler {fehl_gam:.0e}, liegt in der U(2) der Ströme) "
      f"vertauschen mit dem KMS-Zustand ({kom_W:.0e}, {kom_G:.0e}); ihre Mittelungen ändern nichts, Kopplung {kG_W:.0e} und "
      f"{kG_G:.0e}. Die kontinuierliche Transfer-SU(2) dagegen bricht die native Ladungsenergie (||[Q^2/2, E^+]|| = "
      f"{kom_HE:.1f}), und der Twirl verändert eine beobachtbare Nahtgröße: P(N_A != 1) = {ion_G:.1e} im KMS-Zustand gegen "
      f"{ion_T:.4f} = (2/3)(1 + 2 e^-4pi)/(4 + 2 e^-4pi) nach dem Twirl. Der relative Transferrahmen ist über die Nahtladung "
      "(Uhr i^(Q_A)) beobachtbar, anders als die Kette: Die Mittelung ist keine Eichredundanz der nativen Daten",
      fehl_gam < 1e-12 and kom_W < 1e-14 and kom_G < 1e-12 and abs(kG_W) < 1e-10 and abs(kG_G) < 1e-10 and kom_HE > 0.5
      and ion_G < 1e-5 and abs(ion_T - (2/3)*(1 + 2*QQ)/(4 + 2*QQ)) < 1e-10)

# eichinvariante Lesart: Vielfachheitsterm
PS = {S: np.trace(rhoG @ P).real for S, P in bl2.items()}
PS_T = {S: np.trace(rhoT @ P).real for S, P in bl2.items()}
CU = sub2(sum(M2["E"][(i, j)] @ M2["E"][(j, i)] for i in range(2) for j in range(2)))
inv_gleich = max(abs(PS[S] - PS_T[S]) for S in PS) + abs(np.trace(rhoG @ CU) - np.trace(rhoT @ CU))
k_alg = -np.log(PS[0.0]/1) + np.log(PS[2.0]/3)                       # je Spinzustand der invarianten Algebra
diff_rand = kT - k_alg
ZAHLEN["1a_eichinvariant"] = {"P(S=0)": PS[0.0], "P(S=1)": PS[2.0], "Kopplung algebraisch": k_alg,
                              "Kopplung Hilbertraum (z116)": kT, "Differenz": diff_rand}
check(f"[E] (a-iii) Eichinvariante Lesart: Auf transferinvarianten Observablen (Spinprojektoren, U(2)-Casimir) sind KMS- "
      f"und getwirlter Zustand gleich (Abweichung {inv_gleich:.0e}); P(S=0) = {PS[0.0]:.7f} = (1 + 2 e^-4pi)/(4 + 2 e^-4pi). "
      f"Der modulare Hamiltonoperator der invarianten Algebra (Spur je Spinzustand) hat die Kopplung {k_alg:.3e} = "
      f"-log(1 + 2 e^-4pi), also praktisch null und mit antiferromagnetischem Vorzeichen; die z116-Kopplung "
      f"log(3/(1 + 2 e^-4pi)) ist dieselbe Größe mit der Hilbertraumspur, die Differenz ist {diff_rand:.12f} = log 3 = "
      "log dim(Isospin 1): der Vielfachheitsterm der nichtabelschen Darstellung (wie log dim R in der Verschränkung von "
      "Eichtheorien), keine Energie. Wer den Transferrahmen als Eichung liest, darf die feste Ladung (1,1) nicht projizieren "
      "und erhält keinen Austausch; wer projiziert, braucht die Mittelung als eigene Prämisse",
      inv_gleich < 1e-12 and abs(k_alg + np.log(1 + 2*QQ)) < 1e-12 and abs(diff_rand - np.log(3)) < 1e-12)

# Vergleich: Hüpfdynamik (nicht nativ) gibt das entgegengesetzte Vorzeichen
hub = {}
for t_ in (0.05, 0.1):
    Dh = sla.expm(-BETA*(M2["H"] - t_*sub2(X01)))
    hub[t_] = kopplung(auf_11(Dh/np.trace(Dh)))
ZAHLEN["1a_huepfen_vergleich"] = {str(t_): v for t_, v in hub.items()}
check(f"[E, Vergleich] (a-iv) Mit einem gewählten Hüpfterm -t (E^+ + E^-) (z99: nicht nativ, ein chemisches Potential der "
      f"Nullmoden) gibt der eingeschränkte KMS-Zustand Superaustausch mit dem entgegengesetzten Vorzeichen: Kopplung "
      f"{ {t_: round(v, 5) for t_, v in hub.items()} } (t = 0.05, 0.1; Singulett bevorzugt, ~ -beta 4 t^2 = "
      f"{-BETA*4*0.05**2:.4f}, {-BETA*4*0.1**2:.4f}). Der Twirl ist also auch kein Grenzfall einer Transferdynamik",
      hub[0.05] < 0 and hub[0.1] < hub[0.05] and abs(hub[0.1]/(-BETA*4*0.01) - 1) < 0.25)
befund("C: Die Mittelung über die Transfer-SU(2) ist nicht nativ. Der KMS-Zustand (auch jeder additiven Nahtdynamik) ist auf "
       "festen Ladungen ein Produkt; diskrete Eichungen (Kette, Tausch) ändern nichts; die kontinuierliche Transfer-SU(2) "
       "ist durch die native Ladungsenergie gebrochen und ihr Twirl ändert beobachtbare Nahtladungen. Eichinvariant gelesen "
       "verschwindet die Kopplung (-log(1 + 2 e^-4pi)); log 3 ist der Vielfachheitsterm log dim(Isospin 1) der "
       "Hilbertraumspur. Der Austausch der Analyse 16 hängt an der Prämisse P_twirl ('Gewicht eines Zustands fester Ladung = "
       "Haar-Mittel über seine Transferbahn')", True)

# ---------------------------------------------------------------- (b)
say("  (b) N Nähte: Paarkopplung, Mehrkörperterme, z105 und z108")


def dim_un_zwei_spalten(N, a, b):
    """Dimension der U(N)-Darstellung mit zwei Spalten der Längen a >= b (Haken-Inhalt)."""
    zl = [2 if i < b else 1 for i in range(a)]
    sp_ = [a, b]
    num = den = 1
    for i, rl in enumerate(zl):
        for j in range(rl):
            num *= N + j - i
            den *= (rl - j - 1) + (sp_[j] - i - 1) + 1
    return Fr(num, den)


def f_zwei_zeilen(N, b):
    return comb(N, b) - (comb(N, b - 1) if b >= 1 else 0)


def J_formel(N):
    """Paarkopplung (Koeffizient von Summe s_i.s_j) von K(S) = log(dim W/f) im Grenzfall beta -> unendlich."""
    tot = 0.0
    tab = {}
    for b in range(N//2 + 1):
        a = N - b
        S = Fr(a - b, 2)
        dW, f = dim_un_zwei_spalten(N, a, b), f_zwei_zeilen(N, b)
        k = float(np.log(float(dW/f)))
        ss = float((S*(S + 1) - Fr(3*N, 4))/2/comb(N, 2))
        tot += k*float((2*S + 1)*f)*ss
        tab[float(S)] = (int(dW), f, k)
    return 16/(3*2**N)*tot, tab


def pauli_zerlegung(K, N):
    koeff = {}
    for ix in itertools.product(range(4), repeat=N):
        cf = np.trace(K @ kr(*[PAULI[i] for i in ix])).real/2**N
        if abs(cf) > 1e-12:
            koeff[ix] = cf
    koerper = {}
    for ix, cf in koeff.items():
        b = sum(1 for i in ix if i)
        koerper[b] = koerper.get(b, 0.0) + cf**2
    J = 4*koeff.get(tuple([1, 1] + [0]*(N - 2)), 0.0)
    return {b: float(np.sqrt(v)) for b, v in sorted(koerper.items())}, J


def vielfachfrei_ok(bl, CUN):
    """Je Spinblock eine einzige U(N)-Darstellung: der U(N)-Casimir ist auf dem Block konstant."""
    for P in bl.values():
        w_, V_ = np.linalg.eigh(P)
        Vs = V_[:, w_ > 0.5]
        Cb = Vs.conj().T @ CUN @ Vs
        if np.abs(P @ CUN - CUN @ P).max() > 1e-9 or np.ptp(np.linalg.eigvalsh((Cb + Cb.conj().T)/2)) > 1e-9:
            return False
    return True


erg_b = {}
for N in (2, 3, 4):
    M = nahtmodell(N) if N > 2 else M2
    bl = spinbloecke(M["S2"])
    D = np.diag(np.exp(-BETA*np.diag(M["H"])))
    CUN = M["sub"](sum(M["E"][(i, j)] @ M["E"][(j, i)] for i in range(N) for j in range(N)))
    vielfachfrei = vielfachfrei_ok(bl, CUN)
    TD = twirl_schur(D, bl)
    w = {S: float((np.trace(D @ P)/np.trace(P)).real) for S, P in bl.items()}
    Jf, tab = J_formel(N)
    KH = {}
    for S2v, P in bl.items():
        Sv = (-1 + np.sqrt(1 + 4*S2v))/2
        dW, f, klim = tab[round(Sv*2)/2]
        KH[float(Sv)] = (float(-np.log(w[S2v])), float(klim))
    Bp = M["B"]
    Kphys = Bp.T @ (-mlog(TD)) @ Bp
    koerper, J = pauli_zerlegung(Kphys, N)
    kom_E = max(np.abs(TD @ M["sub"](M["E"][(i, j)]) - M["sub"](M["E"][(i, j)]) @ TD).max()
                for i in range(N) for j in range(N))
    eintrag = {"Gewichte w_S": {str(round((-1 + np.sqrt(1 + 4*k))/2, 2)): v for k, v in w.items()},
               "K_H(S) gegen log(dimW/f)": {str(k): v for k, v in KH.items()}, "Körpernormen": koerper,
               "J (Pauli)": J, "J (Formel, beta->inf)": Jf, "vielfachfrei": bool(vielfachfrei), "[T(D), E_ij]": float(kom_E)}
    if N == 3:
        hs = np.zeros_like(TD, dtype=complex)
        for _ in range(300):
            G_ = fock_u(M, haar_un(3))
            hs += G_ @ D @ G_.conj().T/300
        eintrag["Haar-Abweichung (300)"] = float(np.abs(hs - TD).max())
        rt = Bp.T @ TD @ Bp
        rt = rt/np.trace(rt)
        dims = [2, 2, 2]
        I3 = (entropie(ptrace(rt, dims, [0, 2])) + entropie(ptrace(rt, dims, [1, 2])) - entropie(ptrace(rt, dims, [2]))
              - entropie(rt))
        rk = Bp.T @ D @ Bp
        rk = rk/np.trace(rk)
        I3k = (entropie(ptrace(rk, dims, [0, 2])) + entropie(ptrace(rk, dims, [1, 2])) - entropie(ptrace(rk, dims, [2]))
               - entropie(rk))
        eintrag["I(A:B|C) Twirl"], eintrag["I(A:B|C) KMS"] = I3, I3k
    erg_b[N] = eintrag
J_tab = {N: J_formel(N)[0] for N in range(2, 13)}
ZAHLEN["1b_N_naehte"] = {str(N): v for N, v in erg_b.items()}
ZAHLEN["1b_J_N"] = {str(N): v for N, v in J_tab.items()}
ok_b = (all(e["vielfachfrei"] for e in erg_b.values())
        and all(abs(a_ - b_) < 1e-4 for e in erg_b.values() for a_, b_ in e["K_H(S) gegen log(dimW/f)"].values())
        and all(abs(e["J (Pauli)"] - e["J (Formel, beta->inf)"]) < 1e-4 for e in erg_b.values())
        and set(erg_b[3]["Körpernormen"]) == {0, 2} and set(erg_b[4]["Körpernormen"]) == {0, 2, 4}
        and erg_b[4]["Körpernormen"][4] > 1e-3 and erg_b[3]["Haar-Abweichung (300)"] < 0.05
        and all(J_tab[N + 1] > J_tab[N] for N in range(2, 12)) and abs(J_tab[2] + np.log(3)) < 1e-12)
check(f"[E] (b-i) N Nähte mit je einer Anregung: Spin S und Transfer-U(N) sind Howe-dual und vielfachfrei (U(N)-Casimir je "
      f"Spinblock konstant), der Twirl (Schur; Haar-Stichprobe der U(3), Abweichung "
      f"{erg_b[3]['Haar-Abweichung (300)']:.3f}) gibt K_H(S) = log(dim W_lambda'/f_lambda) + O(e^-4pi): N = 3: "
      f"{ {k: round(v[0], 6) for k, v in erg_b[3]['K_H(S) gegen log(dimW/f)'].items()} } (log 4 = log(8/2)), N = 4: "
      f"{ {k: round(v[0], 6) for k, v in erg_b[4]['K_H(S) gegen log(dimW/f)'].items()} } (log 10, log 5, 0). Paulizerlegung "
      f"auf (C^2)^N: N = 3 nur Zweikörper (Normen {erg_b[3]['Körpernormen']}), N = 4 dazu ein Vierkörperterm "
      f"({erg_b[4]['Körpernormen']}), ungerade Körperzahlen nie. Paarkopplung J_N (Koeffizient von Summe s_i.s_j): "
      f"{ {N: round(v, 5) for N, v in J_tab.items()} }, N J_N = "
      f"{ {N: round(N*v, 3) for N, v in J_tab.items() if N in (2, 3, 4, 6, 8, 12)} }: alle-mit-allen, ferromagnetisch, "
      "mit 1/N fallend (Mittelfeld), keine feste Kopplung", ok_b)
check(f"[E negativ] (b-ii) Vergleich z105/z108: Der getwirlte Boltzmann-Operator vertauscht mit allen Transferströmen E_ij "
      f"(max {max(e['[T(D), E_ij]'] for e in erg_b.values()):.0e}) und mit jeder Besetzungszahl; er hat keine Matrixelemente "
      "zwischen Sektoren verschiedener Anregungszahl, also keine 1->3-Beimischung, gamma = 0 und kein Box-Nachkomme "
      "(||det P||^2 = 2 gamma (1 + gamma) = 0): die phi^3-Bedingung von z105 ist nicht erfüllt. Im Sinn von z108 liegt er im "
      f"nichtadditiven Teil: I(A:B|C) = {erg_b[3]['I(A:B|C) Twirl']:.4f} > 0 für N = 3 (eingeschränkter KMS-Zustand "
      f"{erg_b[3]['I(A:B|C) KMS']:.0e}), aber diagonal in den Besetzungen wie der z94-Zustand: korreliert, ohne Übergang",
      all(e["[T(D), E_ij]"] < 1e-10 for e in erg_b.values()) and erg_b[3]["I(A:B|C) Twirl"] > 1e-3
      and abs(erg_b[3]["I(A:B|C) KMS"]) < 1e-10)

# ---------------------------------------------------------------- (c)
say("  (c) Statistik des Dubletts, Gamma(pi)/P_ident, Amplituden")


def bosonmodell():
    """Zwei Bosonen in vier Moden (A hoch, A runter, B hoch, B runter): Besetzungsbasis, Einkörperoperatoren."""
    basis = [n for n in itertools.product(range(3), repeat=4) if sum(n) == 2]
    pos = {n: i for i, n in enumerate(basis)}

    def eins(i, j):
        Mx = np.zeros((len(basis), len(basis)))
        for n in basis:
            if n[j] == 0:
                continue
            amp = np.sqrt(n[j])
            m_ = list(n)
            m_[j] -= 1
            amp *= np.sqrt(m_[i] + 1)
            m_[i] += 1
            Mx[pos[tuple(m_)], pos[n]] += amp
        return Mx
    num = [eins(i, i) for i in range(4)]
    Sx = (eins(0, 1) + eins(1, 0) + eins(2, 3) + eins(3, 2))/2
    Sy = (-1j*(eins(0, 1) - eins(1, 0)) - 1j*(eins(2, 3) - eins(3, 2)))/2
    Sz = (num[0] - num[1] + num[2] - num[3])/2
    NA_, NB_ = num[0] + num[1], num[2] + num[3]
    H = ((NA_ - np.eye(len(basis))) @ (NA_ - np.eye(len(basis))) + (NB_ - np.eye(len(basis))) @ (NB_ - np.eye(len(basis))))/2
    return basis, (Sx @ Sx + Sy @ Sy + Sz @ Sz).real, H, NA_, NB_


bb, S2b, Hb, NAb, NBb = bosonmodell()
blb = spinbloecke(S2b)
Db = np.diag(np.exp(-BETA*np.diag(Hb)))
wb = {S: float(np.trace(Db @ P)/np.trace(P)) for S, P in blb.items()}
w_f = {S: float((np.trace(D2 @ P)/np.trace(P)).real) for S, P in bl2.items()}
k_bos = np.log(wb[0.0]) - np.log(wb[2.0])                            # K(Triplett) - K(Singulett)
GamP = B2.T @ GamS @ B2
eig_s = float(np.trace(GamP @ P_sing).real)
eig_t = float(np.trace(GamP @ P_trip).real/3)
Pid = (np.eye(4) + GamP)/2
rang_id = int(round(np.trace(Pid).real))
ZAHLEN["1c"] = {"Fermionen w_S": {str(k): v for k, v in w_f.items()}, "Bosonen w_S": {str(k): v for k, v in wb.items()},
                "Bosonen K(Triplett)-K(Singulett)": k_bos, "Gamma auf Singulett": eig_s, "Gamma auf Triplett": eig_t,
                "Rang P_ident auf (1,1)": rang_id}
check(f"[E] (c-i) Das Vorzeichen kommt aus der Statistik des Dubletts, nicht aus der Dynamik: fermionisch (zwei CAR-Moden je "
      f"Naht) w_S = { {f'S={round((-1 + np.sqrt(1 + 4*k))/2)}': round(v, 7) for k, v in w_f.items()} } (Spin-Singulett unterdrückt, ferromagnetisch, z116); bosonisch "
      f"(gerades Dublett) w_S = { {f'S={round((-1 + np.sqrt(1 + 4*k))/2)}': round(v, 7) for k, v in wb.items()} }: dort ist das Spin-Triplett an Isospin 1 gebunden "
      f"und unterdrückt, K = {k_bos:.7f} P_Triplett (antiferromagnetisch). Gamma(pi) hat auf (1,1) die Eigenwerte "
      f"{eig_s:+.0f} (Singulett) und {eig_t:+.0f} (Triplett), P_ident (z106) lässt genau {rang_id} Zustand, das Singulett "
      "j = 0, genau den Zustand, den der Twirl um 1/3 unterdrückt; mit P_ident trägt der Austauschterm nur eine Konstante",
      abs(w_f[0.0] - (1 + 2*QQ)/3) < 1e-12 and abs(w_f[2.0] - 1) < 1e-12 and abs(wb[2.0] - (1 + 2*QQ)/3) < 1e-12
      and abs(wb[0.0] - 1) < 1e-12 and abs(eig_s - 1) < 1e-12 and abs(eig_t + 1) < 1e-12 and rang_id == 1)
K11 = -mlog(rT11)
S11 = [B2.T @ sub2(x) @ B2 for x in M2["S"]]
kom_spin = max(np.abs(K11 @ x - x @ K11).max() for x in S11)
I16 = np.eye(16)
Dfull = np.diag(np.exp(-BETA*np.diag(sum((n_ - I16) @ (n_ - I16) for n_ in M2["Ni"])/2)))
TDf = np.zeros((16, 16), complex)
for _ in range(300):
    Lg = sla.logm(haar_un(2))
    G_ = sla.expm(sum(Lg[i, j]*M2["E"][(i, j)] for i in range(2) for j in range(2)))
    TDf += G_ @ Dfull @ G_.conj().T/300
nd = np.diag(M2["Ntot"])
quer = float(np.abs(TDf[nd[:, None] != nd[None, :]]).max())
eins = TDf[np.ix_(nd == 1, nd == 1)]
eins_streu = float(np.abs(eins - np.trace(eins)/4*np.eye(4)).max())
ZAHLEN["1c_amplituden"] = {"[K_eff, S_a]": float(kom_spin), "Twirl voller Fockraum: Elemente zwischen Anregungszahlen": quer,
                           "Ein-Anregungs-Sektor: Abweichung von konstant": eins_streu}
check(f"[E] (c-ii) Keine verbundene 1->2-Amplitude: K_eff = c (1/4 - s_A.s_B) lebt auf festen Besetzungen und vertauscht mit "
      f"dem Gesamtspin (max {kom_spin:.0e}); der über die echte Transfer-U(2) gemittelte Boltzmann-Operator auf dem ganzen "
      f"Fockraum (alle Anregungszahlen 0..4, 300 Haar-Elemente) hat zwischen verschiedenen Anregungszahlen die Elemente "
      f"max {quer:.0e} und ist auf dem Ein-Anregungs-Sektor konstant (Abweichung {eins_streu:.0e}); mit Teil b vertauscht er "
      "mit allen E_ij. Der Term ist ein Modulgewicht (Zweikörper-Energieverschiebung im modularen Hamiltonoperator), sein "
      "Fluss e^(-i t K) dreht nur Spins bei festen Besetzungen; Transfers kommen nur in zweiter Ordnung ((1,1) -> (2,0) -> "
      "(1,1)) und kehren in den Ausgangssektor zurück",
      kom_spin < 1e-12 and quer < 1e-12 and eins_streu < 1e-12)
befund("C/L: Das Seam-Dublett s_A muss ein Paar unabhängiger CAR-Moden auf derselben Ladungsstufe sein, damit die "
       "Pauli-Kopplung Spin-Singulett = Isospin 1 gilt; dann ist der Term ferromagnetisch. Das schwache Dublett C^2 c C^5 der "
       "16 erfüllt das formal (der Term wäre ein Isospin-Isospin-Modulgewicht ohne Eichboson); die Schrittart aus z83/z106 "
       "ist ein gerades Etikett am ungeraden Ladungsschritt, ihre Doppelbesetzung auf einer Naht ist nicht durch zwei Moden "
       "gleicher Stufe gegeben, die Pauli-Kopplung also nicht gezeigt; und für zwei Weyl-Quanten an der Schwelle lässt "
       "P_ident ohnehin nur j = 0. In keiner Lesart entsteht eine 1->2-Amplitude", True)

# =================================================================== Teil 2
say("\nTeil 2: Petz-Rekonstruierbarkeit statt P_koh")


def kraus_aus(Uv, Dv):
    m, r = Uv.shape
    Ks = np.zeros((r, m + 1, m), complex)
    Ks[:, np.arange(1, m + 1), np.arange(m)] = Uv.T
    Ks[:, np.arange(m), np.arange(m)] = Dv.T
    return Ks


def petz_F(Ks, sig):
    """Verschränkungsfidelität von Petz o Kanal, Referenz diag(sig), Purifikation sum sqrt(sig_k) |kk>."""
    S = np.diag(sig)
    tau = np.einsum("iab,bc,idc->ad", Ks, S, Ks.conj())
    ti = msqrt(tau, inv=True)
    Z = np.einsum("ab,ibc,cd->iad", ti, Ks, np.diag(sig**1.5))
    Mm = Ks.conj().reshape(len(Ks), -1) @ Z.reshape(len(Ks), -1).T
    return float(np.sum(np.abs(Mm)**2))


def koh_kraus(p, q):
    m = len(p)
    Uv, Dv = np.zeros((m, 2)), np.zeros((m, 2))
    Uv[:, 0], Dv[:, 1] = np.sqrt(p), np.sqrt(q)
    return kraus_aus(Uv, Dv)


def F_koh_formel(p, q, sig):
    m = len(p)
    tau = np.array([(sig[j - 1]*p[j - 1] if j >= 1 else 0) + (sig[j]*q[j] if j < m else 0) for j in range(m + 1)])
    a = sig**1.5/np.sqrt(tau[1:])
    b = sig**1.5/np.sqrt(tau[:-1])
    return float((a @ p)**2 + (b @ q)**2), a, b


def vakuum(m):
    p = np.array([(k + 1)/(m + 1) for k in range(m)])
    return p, 1 - p


def thermisch(m, beta, alpha=0.25):
    """z80: c(x) = 1/(1 + e^(beta (r + alpha))), W_m(r) = c(r) - c(r+m), u_m(r) = Summe (m-j) w(r+j)/((m+1) W_m(r))."""
    c = lambda r: 1/(1 + np.exp(np.clip(beta*(r + alpha), -700, 700)))
    rows = list(range(-m - 1, 1))
    W = np.array([c(r) - c(r + m) for r in rows])
    u = np.array([sum((m - j)*(c(r + j) - c(r + j + 1)) for j in range(m))/((m + 1)*(c(r) - c(r + m))) for r in rows])
    return u, 1 - u, W/W.sum()


def zufalls_torus(p, q):
    m = len(p)
    Xu = rng.normal(size=(m, m)) + 1j*rng.normal(size=(m, m))
    Xd = rng.normal(size=(m, m)) + 1j*rng.normal(size=(m, m))
    Xu /= np.linalg.norm(Xu, axis=1, keepdims=True)
    Xd /= np.linalg.norm(Xd, axis=1, keepdims=True)
    Uv = np.hstack([np.sqrt(p)[:, None]*Xu, np.zeros((m, m))])
    Dv = np.hstack([np.zeros((m, m)), np.sqrt(q)[:, None]*Xd])
    return kraus_aus(Uv, Dv)


# ---------------------------------------------------------------- (a) Satz im Torusraum
say("  (a) Torusverträgliche Kanäle mit nativen Raten: Satz und Prüfung")
satz = {}
ok_satz = True
for m in range(2, 7):
    refs = {"gleichförmig": (*vakuum(m), np.ones(m)/m), "zufällig": (*vakuum(m), (lambda s: s/s.sum())(rng.uniform(0.1, 1, m))),
            "KMS z80 beta=4pi": thermisch(m, BETA), "KMS z80 beta=2": thermisch(m, 2.0)}
    for name, (p, q, sig) in refs.items():
        Fk = petz_F(koh_kraus(p, q), sig)
        Ff, a, b = F_koh_formel(p, q, sig)
        zuf = [petz_F(zufalls_torus(p, q), sig) for _ in range(150)]
        ph = np.exp(2j*np.pi*rng.uniform(size=len(p)))
        ph2 = np.exp(2j*np.pi*rng.uniform(size=len(p)))
        Fg = petz_F(kraus_aus(np.stack([np.sqrt(p)*ph, 0*p], 1), np.stack([0*q, np.sqrt(q)*ph2], 1)), sig)
        n_ = len(p)
        eps = 1e-3
        Xu = np.zeros((n_, 2), complex)
        Xu[:, 0] = 1
        Xu[0, 1] = eps
        Xu /= np.linalg.norm(Xu, axis=1, keepdims=True)
        Fs = petz_F(kraus_aus(np.hstack([np.sqrt(p)[:, None]*Xu, np.zeros((n_, 1))]),
                              np.hstack([np.zeros((n_, 2)), np.sqrt(q)[:, None]])), sig)
        satz[(m, name)] = {"F_koh": Fk, "Formel": Ff, "max zufällig": max(zuf), "geeicht": Fg, "gestört (eps=1e-3)": Fs}
        ok_satz &= (abs(Fk - Ff) < 1e-10*max(1, Fk) and max(zuf) < Fk - 1e-9 and abs(Fg - Fk) < 1e-10 and Fs < Fk)
        if name == "gleichförmig":
            ok_satz &= abs(Fk - (m + 1)/(2*m)) < 1e-12
ZAHLEN["2a_satz"] = {f"{m}|{n}": v for (m, n), v in satz.items()}
check("[E] (a) Satz: Im torusverträglichen Kanalraum (Krausoperatoren je mit festem Torusgewicht, d.h. Gram = G_u (+) G_d, "
      "Spurerhaltung u_k ⟂ d_(k+1), native Raten |u_k|^2 = p_k, |d_k|^2 = q_k) ist Phi(sigma) diagonal und die "
      "Petz-Verschränkungsfidelität F = Tr M^2 mit M = Summe a_k u_k u_k^+ + b_k d_k d_k^+, a_k = sigma_k^(3/2) "
      "tau_(k+1)^(-1/2), b_k = sigma_k^(3/2) tau_k^(-1/2). Cauchy-Schwarz |<u_k,u_l>|^2 <= p_k p_l gibt F <= (Summe a p)^2 + "
      "(Summe b q)^2 = F_koh mit Gleichheit genau, wenn alle u_k und alle d_k parallel sind: der kohärente Kanal, bis auf "
      "diagonale Eichung. Geprüft m = 2..6 für gleichförmige, zufällige und thermische KMS-Gewichte (z80, beta = 4 pi und 2, "
      "mit den thermischen Raten): Formel = Petz numerisch, 150 zufällige Toruskanäle je Fall stets darunter, Eichphasen "
      "ändern nichts, eine Störung der Ausrichtung senkt F; gleichförmig F_koh = (m+1)/(2m). Beispiele: "
      + "; ".join(f"m = {m} {n}: F_koh {v['F_koh']:.6f}, max zufällig {v['max zufällig']:.6f}"
                  for (m, n), v in satz.items() if m == 3), ok_satz)

# ---------------------------------------------------------------- (b) ohne Torus
say("  (b) Ohne Torusverträglichkeit: Blockkanäle und freie Optimierung")


def partitionen(items):
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for P in partitionen(rest):
        yield [[first]] + P
        for i in range(len(P)):
            yield P[:i] + [[first] + P[i]] + P[i + 1:]


def verboten(m, zustand=True):
    f = {frozenset({("u", k), ("d", k + 1)}) for k in range(m - 1)}
    if zustand:
        f |= {frozenset({("u", k), ("d", k)}) for k in range(m)}
    return f


def blockkanal(P, p, q):
    m = len(p)
    Uv, Dv = np.zeros((m, len(P))), np.zeros((m, len(P)))
    for bi, Bl in enumerate(P):
        for typ, k in Bl:
            if typ == "u":
                Uv[k, bi] = np.sqrt(p[k])
            else:
                Dv[k, bi] = np.sqrt(q[k])
    return kraus_aus(Uv, Dv)


def spiegel(v, m):
    return ("d", m - 1 - v[1]) if v[0] == "u" else ("u", m - 1 - v[1])


def superop_petz(Ks, sig):
    m = Ks.shape[2]
    S = np.diag(sig)
    Phi = lambda X: np.einsum("iab,bc,idc->ad", Ks, X, Ks.conj())
    Phid = lambda Y: np.einsum("iba,bc,icd->ad", Ks.conj(), Y, Ks)
    ti = msqrt(Phi(S), inv=True)
    sq = msqrt(S)
    R = lambda Y: sq @ Phid(ti @ Y @ ti) @ sq
    Sup = np.zeros((m*m, m*m), complex)
    for a in range(m*m):
        E = np.zeros(m*m)
        E[a] = 1
        Sup[:, a] = R(Phi(E.reshape(m, m))).reshape(-1)
    return Sup, Phi


blk = {}
ok_blk = True
for m in (2, 3, 4, 5):
    p, q = vakuum(m)
    pf, qf = [Fr(k + 1, m + 1) for k in range(m)], [Fr(m - k, m + 1) for k in range(m)]
    items = [("u", k) for k in range(m)] + [("d", k) for k in range(m)]
    wt = {("u", k): pf[k] for k in range(m)}
    wt.update({("d", k): qf[k] for k in range(m)})
    vb = verboten(m)
    c2 = Fr(m + 1, m**3)
    best, bestP, bestS, bestSP = Fr(-1), None, Fr(-1), None
    for P in partitionen(items):
        if any(frozenset({x, y}) in vb for Bl in P for x in Bl for y in Bl if x != y):
            continue
        sc = sum(sum(wt[v] for v in Bl)**2 for Bl in P)
        if sc > best:
            best, bestP = sc, P
        if sorted(sorted(spiegel(v, m) for v in Bl) for Bl in P) == sorted(sorted(Bl) for Bl in P) and sc > bestS:
            bestS, bestSP = sc, P
    Ks = blockkanal(bestP, p, q)
    sig = np.ones(m)/m
    Fb = petz_F(Ks, sig)
    tp = np.abs(np.einsum("iba,ibc->ac", Ks.conj(), Ks) - np.eye(m)).max()
    Sup, Phi = superop_petz(Ks, sig)
    zust = np.abs(Phi(np.eye(m)/m) - np.eye(m + 1)/(m + 1)).max()
    rq = (np.abs(Ks)**2).sum(0)
    raten = max(np.abs(rq[np.arange(1, m + 1), np.arange(m)] - p).max(), np.abs(rq[np.arange(m), np.arange(m)] - q).max())
    petzsym = np.abs(Sup - Sup.conj().T).max()
    gemischt = any(any(v[0] == "u" for v in Bl) and any(v[0] == "d" for v in Bl) for Bl in bestP)
    KsS = blockkanal(bestSP, p, q)
    FbS = petz_F(KsS, sig)
    PhiS = superop_petz(KsS, sig)[1]
    Xr = rng.normal(size=(m, m)) + 1j*rng.normal(size=(m, m))
    Fm, Fm1 = np.eye(m)[::-1], np.eye(m + 1)[::-1]
    spiegel_fehl = np.abs(PhiS(Fm @ Xr @ Fm) - Fm1 @ PhiS(Xr) @ Fm1).max()
    j_fehl = np.abs(PhiS(Xr.conj()) - PhiS(Xr).conj()).max() + np.abs(Phi(Xr.conj()) - Phi(Xr).conj()).max()
    blk[m] = {"F_koh": (m + 1)/(2*m), "F_Block": Fb, "exakt": str(c2*best), "Delta": str(c2*best - c2*Fr(m*m, 2)),
              "Partition": [[f"{t}{k}" for t, k in Bl] for Bl in bestP], "F_Block_spiegel": FbS,
              "Partition_spiegel": [[f"{t}{k}" for t, k in Bl] for Bl in bestSP], "Spiegelfehler": float(spiegel_fehl),
              "J-Fehler": float(j_fehl)}
    ok_blk &= (abs(Fb - float(c2*best)) < 1e-12 and tp < 1e-12 and zust < 1e-12 and raten < 1e-12 and petzsym < 1e-12
               and (gemischt or m == 2) and spiegel_fehl < 1e-12 and j_fehl < 1e-12)
# thermische Referenz (z80, beta = 4 pi) und zufällige Referenz bei m = 2: Formel F = Summe_Block (Summe c_v w_v)^2
th = {}
for m, (lab, dat) in itertools.product((2, 3), (("KMS beta=4pi", "th"), ("zufällig", "zu"))):
    if dat == "th":
        p, q, sig = thermisch(m, BETA)
    else:
        p, q = vakuum(m)
        sig = np.array([0.62, 0.38]) if m == 2 else np.array([0.5, 0.3, 0.2])
    n_ = len(p)
    Fk, a, b = F_koh_formel(p, q, sig)
    coef = {("u", k): a[k]*p[k] for k in range(n_)}
    coef.update({("d", k): b[k]*q[k] for k in range(n_)})
    vb = verboten(n_)
    best, bestP = -1, None
    for P in partitionen(sorted(coef)):
        if any(frozenset({x, y}) in vb for Bl in P for x in Bl for y in Bl if x != y):
            continue
        sc = sum(sum(coef[v] for v in Bl)**2 for Bl in P)
        if sc > best:
            best, bestP = sc, P
    th[(m, lab)] = {"F_koh": Fk, "F_Block": float(best), "Petz numerisch": petz_F(blockkanal(bestP, p, q), sig)}
    ok_blk &= abs(th[(m, lab)]["Petz numerisch"] - best) < 1e-9*max(1, best)
# nur Spurerhaltung und Raten (ohne Zustandserhaltung): freie Optimierung (m = 2) und Blocksuche (m = 2..4)


def vek_aus_param(x, m, r, p, q):
    z = x[:2*m*r].reshape(2, m, r) + 1j*x[2*m*r:].reshape(2, m, r)
    Uv, Dv = z[0].copy(), z[1].copy()
    for k in range(m):
        Uv[k] *= np.sqrt(p[k])/np.linalg.norm(Uv[k])
    for k in range(m):
        if k >= 1:
            Dv[k] -= Uv[k - 1]*(np.vdot(Uv[k - 1], Dv[k])/np.vdot(Uv[k - 1], Uv[k - 1]))
        Dv[k] *= np.sqrt(q[k])/np.linalg.norm(Dv[k])
    return Uv, Dv


nur_tp = {}
for m in (2, 3, 4):
    p, q = vakuum(m)
    sig = np.ones(m)/m
    items = [("u", k) for k in range(m)] + [("d", k) for k in range(m)]
    vb = verboten(m, zustand=False)
    best = -1
    for P in partitionen(items):
        if any(frozenset({x, y}) in vb for Bl in P for x in Bl for y in Bl if x != y):
            continue
        Fv = petz_F(blockkanal(P, p, q), sig)
        if Fv > best:
            best, bP = Fv, P
    eintrag = {"F_koh": (m + 1)/(2*m), "Block": best, "Partition": [[f"{t}{k}" for t, k in Bl] for Bl in bP]}
    if m == 2:
        r = 2*m
        bo = -1
        for _ in range(8):
            res = minimize(lambda x: -petz_F(kraus_aus(*vek_aus_param(x, m, r, p, q)), sig), rng.normal(size=4*m*r),
                           method="L-BFGS-B", options={"maxiter": 500})
            bo = max(bo, -res.fun)
        eintrag["frei optimiert"] = bo
    nur_tp[m] = eintrag
ZAHLEN["2b_bloecke_zustandserhaltend"] = {str(k): v for k, v in blk.items()}
ZAHLEN["2b_bloecke_thermisch_zufaellig"] = {f"{m}|{n}": v for (m, n), v in th.items()}
ZAHLEN["2b_nur_spurerhaltend"] = {str(k): v for k, v in nur_tp.items()}
ok_neg = (all(Fr(blk[m]["Delta"]) > 0 for m in (3, 4, 5)) and Fr(blk[2]["Delta"]) == 0
          and blk[4]["F_Block_spiegel"] > blk[4]["F_koh"] + 1e-9
          and th[(3, "KMS beta=4pi")]["F_Block"] > th[(3, "KMS beta=4pi")]["F_koh"] + 1e-9
          and th[(2, "zufällig")]["F_Block"] > th[(2, "zufällig")]["F_koh"] + 1e-9
          and nur_tp[2]["frei optimiert"] > 0.85 and nur_tp[3]["Block"] > 0.82)
check(f"[E negativ] (b) Die Eindeutigkeitsbehauptung der Analyse 16 gilt ohne Torus nicht. Exakte Blockkanäle (je Block ein "
      f"Krausoperator mit den Wurzeln der nativen Raten; erlaubt sind alle Paare außer u_k-d_k und u_k-d_(k+1)) erfüllen "
      f"Raten, Spurerhaltung, Zustandserhaltung I/m -> I/(m+1), J-Realität und Petz-Symmetrie (alle Fehler < 1e-12) und "
      f"übertreffen F_koh: Delta F = { {m: blk[m]['Delta'] for m in blk} } (m = 2 Gleichstand), etwa m = 3: "
      f"{blk[3]['F_Block']:.6f} gegen {blk[3]['F_koh']:.6f} mit {blk[3]['Partition']}; spiegelverträglich (Partition "
      f"spiegelinvariant) noch m = 4: {blk[4]['F_Block_spiegel']:.6f} gegen 0.625 mit {blk[4]['Partition_spiegel']} (für "
      f"m = 3, 5 ist unter den spiegelinvarianten Blockkanälen der kohärente der beste: "
      f"{blk[3]['F_Block_spiegel']:.6f}, {blk[5]['F_Block_spiegel']:.6f}). Mit den "
      f"KMS-Gewichten z80 (beta = 4 pi, m = 3): {th[(3, 'KMS beta=4pi')]['F_Block']:.6f} gegen "
      f"{th[(3, 'KMS beta=4pi')]['F_koh']:.6f}; mit ungleichen Gewichten schon m = 2: {th[(2, 'zufällig')]['F_Block']:.6f} "
      f"gegen {th[(2, 'zufällig')]['F_koh']:.6f}. Nur mit Spurerhaltung und Raten: m = 2 frei optimiert "
      f"{nur_tp[2]['frei optimiert']:.6f} gegen 0.75, m = 3 Block {nur_tp[3]['Block']:.6f} gegen 0.666667 "
      f"({nur_tp[3]['Partition']}: die großen Raten beider Schrittarten in einem kohärenten Krausoperator). Die Maximierer "
      "mischen die Schrittarten kohärent und brechen den Torus", ok_blk and ok_neg)

# ---------------------------------------------------------------- (c) native Dilatation
say("  (c) Ist die native Dilatation des Schritts L^nat der kohärente Kanal?")


def Lnat(m):
    L = np.zeros((m + 1, m))
    for k in range(m):
        L[k + 1, k] = k + 1                                          # Teilchen hoch, Amplitude n_1 + 1
        L[k, k] = m - k                                              # Loch runter, Amplitude n_2 + 1
    return L


nat = {}
ok_nat, ok_pol = True, True
for m in range(2, 9):
    L = Lnat(m)
    Lp, Lh = np.diag(np.diag(L, -1), -1)[:, :m], np.eye(m + 1, m)*np.diag(L)
    LtL = L.T @ L
    off = np.abs(LtL - np.diag(np.diag(LtL))).max()
    raten_ad = (np.diag(L, -1)**2/(np.diag(L, -1)**2 + np.diag(L)[:m]**2))
    p, q = vakuum(m)
    V = L @ msqrt(LtL, inv=True)
    Fv = petz_F(V[None, :, :].astype(complex), np.ones(m)/m)
    band = np.zeros_like(V, dtype=bool)
    band[np.arange(1, m + 1), np.arange(m)] = True
    band[np.arange(m), np.arange(m)] = True
    aus_band = float((np.abs(V[~band])**2).sum()/m)
    polar = {}
    for t in (0.0, 0.25, 0.5, 0.75, 1.0):
        ap, ah = np.arange(1, m + 1)**t, np.arange(m, 0, -1)**t                    # |L^p|^t, |L^h|^t
        nrm = np.sqrt(ap**2 + ah**2)
        Kp = np.zeros((m + 1, m))
        Kh = np.zeros((m + 1, m))
        Kp[np.arange(1, m + 1), np.arange(m)] = ap/nrm
        Kh[np.arange(m), np.arange(m)] = ah/nrm
        polar[t] = float(np.abs((ap/nrm)**2 - p).max())
        if t == 0.5:
            Kk = koh_kraus(p, q)
            choi_ab = np.abs(np.kron(Kp, Kp) + np.kron(Kh, Kh) - sum(np.kron(K_, K_.conj()) for K_ in Kk)).max()
    nat[m] = {"TP-Nebendiagonale (k+1)(m-k-1)": off, "Raten Ad_L - p": float(np.abs(raten_ad - p).max()),
              "F_Petz polare Isometrie": Fv, "Gewicht außerhalb der Schritte": aus_band, "Polarfamilie |Rate - p|": polar,
              "Choi(t=1/2) - Phi_1/2": float(choi_ab)}
    ok_nat &= (np.abs(Lp).sum() > 0 and np.abs(Lh).sum() > 0 and off >= m - 1 - 1e-12 and np.abs(raten_ad - p).max() > 0.05
               and abs(Fv - 1) < 1e-10 and aus_band > 1e-3)
    ok_pol &= polar[0.5] < 1e-12 and all(v > 1e-3 for t, v in polar.items() if t != 0.5) and choi_ab < 1e-12
ZAHLEN["2c_native_dilatation"] = {str(k): v for k, v in nat.items()}
check(f"[E negativ] (c-i) Der native Schritt L^nat_-1 (z83: Amplituden n_1 + 1 Teilchen hoch, n_2 + 1 Loch runter, beide "
      f"positiv) ist ein einziger Fockoperator: Krausrang 1, ohne Register für die Schrittart, und er trägt beide "
      f"Torusgewichte (+1 und 0) kohärent. Er ist mit keiner Konstante spurerhaltend (L^T L hat die Nebendiagonale "
      f"(k+1)(m-k-1), max {nat[8]['TP-Nebendiagonale (k+1)(m-k-1)']:.0f} bei m = 8); Ad_L hat die Raten "
      f"(k+1)^2/((k+1)^2 + (m-k)^2) statt (k+1)/(m+1) (Abweichung {nat[4]['Raten Ad_L - p']:.3f} bei m = 4); seine polare "
      f"Isometrie ist perfekt rekonstruierbar (F = {nat[4]['F_Petz polare Isometrie']:.12f}), springt aber außerhalb der "
      f"nativen Schritte (Gewicht {nat[4]['Gewicht außerhalb der Schritte']:.4f} je Eingang bei m = 4). Das Original liefert "
      "also keinen Kanal der Klasse, insbesondere nicht den kohärenten", ok_nat)
check(f"[E] (c-ii) Gewichtsaufgelöste Dilatation (Register = Hardy-Hälfte, die gehoben hat; Krausoperatoren aus der "
      f"Polarfamilie V_a |L^a|^t, eingangsnormiert): Sie ist torusverträglich und spurerhaltend für jedes t, hat die nativen "
      f"Raten aber genau bei t = 1/2 (Abweichungen m = 5: { {t: round(v, 4) for t, v in nat[5]['Polarfamilie |Rate - p|'].items()} }), "
      f"und dort ist sie Phi_1/2 (Choi-Abstand {max(v['Choi(t=1/2) - Phi_1/2'] for v in nat.values()):.0e}, m = 2..8): die "
      "Kanalform von z91 (iii). Das Register 'nur die Schrittart' ist der Inhalt von P_koh, nicht eine Folge des Originals",
      ok_pol)
diag = {}
ok_diag = True
for m in range(2, 7):
    p, q = vakuum(m)
    sig = np.ones(m)/m
    Lo = Lnat(m - 1)
    Cn = m*(m - 1)*np.eye(m) - Lo @ Lo.T                              # Casimir der L^nat-sl(2,R) auf H_m (z91)
    idd = [a*m + a for a in range(m)]
    ido = [a*m + b for a in range(m) for b in range(m) if a != b]
    spek = {}
    dek = np.concatenate([kraus_aus(np.diag(np.sqrt(p)), np.zeros((m, m))), kraus_aus(np.zeros((m, m)), np.diag(np.sqrt(q)))])
    for name, Ks in (("kohärent", koh_kraus(p, q)), ("zufällig torus", zufalls_torus(p, q)), ("dekohärent", dek)):
        Sup, _ = superop_petz(Ks, sig)
        Cs = m*(m + 1)*(np.eye(m*m) - Sup)
        ok_diag &= (np.abs(Cs[np.ix_(idd, idd)] - Cn).max() < 1e-10 and np.abs(Cs[np.ix_(ido, idd)]).max() < 1e-10)
        ev = np.sort(np.linalg.eigvals(Cs[np.ix_(ido, ido)]).real)
        spek[name] = (round(float(ev.min()), 4), round(float(ev.max()), 4))
    diag[m] = spek
ZAHLEN["2c_petz_komposit_kohaerenzen"] = {str(k): v for k, v in diag.items()}
check(f"[E] (c-iii) Für jeden torusverträglichen Kanal mit nativen Raten ist der Petz-Komposit m(m+1)(1 - R Phi) auf den "
      f"Diagonalen exakt der native Casimir C = L_0(L_0 - 1) - L_-1 L_1 der L^nat-sl(2,R) (m = 2..6, kohärent, zufällig, "
      f"dekohärent): Die nativen Daten (Raten, Torus, Casimir, KMS-Gewichte) sind in der ganzen Klasse dieselben. Die Kanäle "
      f"unterscheiden sich nur auf den Kohärenzen (Spektrum dort, min/max, m = 4: {diag[4]}); genau diese Freiheit wählt "
      "das Petz-Prinzip, das Original wählt sie nicht", ok_diag)
befund("C: Das Petz-Prinzip ersetzt die Hälfte 'Krausrang 2' von P_koh, nicht P_koh ganz. Im torusverträglichen Raum mit "
       "nativen Raten ist maximale Petz-Rekonstruierbarkeit für jede treue diagonale Referenz (auch die KMS-Gewichte) genau "
       "der kohärente Kanal; ohne Torus wählt sie schrittartmischende Blockkanäle. Die Torusverträglichkeit des Kanals "
       "liefert das Original nicht: L^nat trägt beide Torusgewichte kohärent, seine diagonalen Daten sind in der ganzen "
       "Klasse gleich. Also P_koh = P_Torus + P_Petz (bedingt), nicht nativ", True)

# =================================================================== Teil 3
say("\nTeil 3: Bilanz")
bilanz = [
    ("E (negativ)", "KMS-Route: eingeschränkter KMS-Zustand (Ladungsenergie, jede additive Nahtdynamik) auf (1,1) maximal "
                    "gemischt bzw. Produkt; Kopplung 0, I(A:B) = 0 gegen log(3/(1 + 2 e^-4pi)) und 0,2 beim Twirl"),
    ("E (negativ)", "Eichroute: Kette W_12 und Tausch Gamma(pi) fixieren den KMS-Zustand, Kopplung 0; die kontinuierliche "
                    "Transfer-SU(2) bricht Q^2/2, der Twirl hebt P(N_A != 1) von 1,7e-6 auf 1/6: kein Eichrahmen"),
    ("E", "eichinvariant: KMS = Twirl auf invarianten Observablen, Kopplung -log(1 + 2 e^-4pi); log 3 = log dim(Isospin 1) "
          "ist der Vielfachheitsterm der Hilbertraumspur"),
    ("E (Vergleich)", "Hüpfdynamik (nicht nativ): Superaustausch mit entgegengesetztem Vorzeichen"),
    ("E", "N Nähte: K_H(S) = log(dim W/f) + O(e^-4pi), N = 3 rein paarweise, N = 4 mit Vierkörperterm, J_N = -1,099, "
          "-0,924, -0,786, ... ~ -4/N (alle-mit-allen, Mittelfeld)"),
    ("E (negativ)", "z105/z108: anzahlerhaltend, vertauscht mit allen E_ij, keine 1->3-Beimischung, gamma = 0; "
                    "nichtadditiv (I_AB|C > 0), aber diagonal in den Besetzungen"),
    ("E", "Vorzeichen = Statistik des Dubletts (fermionisch ferro, bosonisch antiferro); P_ident lässt nur das Singulett"),
    ("C", "Austausch der Analyse 16 = Prämisse P_twirl; Modulgewicht ohne 1->2-Amplitude"),
    ("E", "Petz-Satz im Torusraum: F = Tr M^2 <= F_koh, Gleichheit genau kohärent, jede treue diagonale Referenz "
          "(gleichförmig, KMS z80 bei 4 pi und 2, zufällig)"),
    ("E (negativ)", "ohne Torus: Blockkanäle mit Raten, Spur- und Zustandserhaltung, J, Petz-Symmetrie übertreffen "
                    "F_koh (Delta 1/54, 1/20, 22/375 für m = 3, 4, 5), spiegelinvariant noch bei m = 4 (0,669 gegen "
                    "0,625; bei m = 3, 5 unter Blockkanälen nicht); nur spurerhaltend 0,881 gegen 0,75"),
    ("E (negativ)", "native Dilatation: L^nat Rang 1, zwei Torusgewichte, nicht spurerhaltend, Ad_L mit quadratischen Raten; "
                    "gewichtsaufgelöst mit nativen Raten genau t = 1/2 = Phi_1/2; diagonale Daten klassengleich"),
    ("C", "P_koh = P_Torus + P_Petz: Petz ersetzt 'Krausrang 2', der Torus bleibt Prämisse"),
    ("O", "eine native Begründung der Torusverträglichkeit des Schrittkanals (etwa Lokalität der Hardy-Zerlegung: Teilchen- "
          "und Lochseite als getrennte Teilsysteme, Schritt als lokale Operation)"),
]
for tag, txt in bilanz:
    say(f"   [{tag}] {txt}")

urteil = ("U3 ist beantwortet, beide Vorschläge der sechzehnten Analyse sind bedingt, keiner ist nativ. (1) Austausch: Die "
          "Mittelung über die Transfer-SU(2) folgt weder aus dem KMS-Zustand noch aus einer Eichredundanz. Der KMS-Zustand der "
          "Ladungsenergie ist auf den festen Ladungen (1,1) maximal gemischt, der jeder additiven Nahtdynamik ein Produkt "
          "(z108); die diskreten Eichungen der gemeinsamen CAR (Kette W_12, Tausch Gamma(pi)) fixieren ihn und geben keine "
          "Kopplung. Die kontinuierliche Transfer-SU(2) ist durch die native Ladungsenergie gebrochen, und ihr Twirl verändert "
          "beobachtbare Nahtladungen (P(N_A != 1) von 1,7e-6 auf 1/6). Eichinvariant gelesen sind KMS und Twirl gleich und die "
          "Kopplung ist -log(1 + 2 e^-4pi), praktisch null; die z116-Kopplung log(3/(1 + 2 e^-4pi)) ist dieselbe Größe mit "
          "der Hilbertraumspur, und die Differenz log 3 = log dim(Isospin 1) ist ein Vielfachheitsterm, keine Energie. Der "
          "Austausch hängt also an der Prämisse P_twirl. Für N Nähte ist K = log(dim W/f) + O(e^-4pi) eine Funktion des "
          "Gesamtspins: bei N = 3 rein paarweise, ab N = 4 mit Vierkörperterm, die Paarkopplung fällt wie 1/N (-1,099, "
          "-0,924, -0,786, ...; alle-mit-allen), also keine feste Heisenbergkopplung. Das Vorzeichen folgt der Statistik "
          "des Dubletts (fermionisch ferromagnetisch, bosonisch antiferromagnetisch); eine Hüpfdynamik gäbe das "
          "entgegengesetzte. Der Term erhält jede Besetzungszahl und vertauscht mit allen Transferströmen: keine "
          "1->2-Amplitude, keine 1->3-Beimischung, gamma = 0, kein Box-Nachkomme (z105); er liegt im nichtadditiven Teil "
          "(I_AB|C > 0), korreliert aber nur wie der z94-Zustand. Für K15 und T3/T5 ändert sich nichts. In der "
          "Szenario-Suche gehört 'Transferaustausch (Twirl)' als bedingter Wert (Prämisse P_twirl, Modulgewicht, "
          "anzahlerhaltend, kein Vertex) auf die Wechselwirkungsachse und zählt nicht als Wechselwirkung für T3/T5. (2) Petz: "
          "Im torusverträglichen Raum mit den nativen Raten ist maximale Petz-Rekonstruierbarkeit exakt und eindeutig der "
          "kohärente Kanal, für jede treue diagonale Referenz einschließlich der KMS-Gewichte aus z80 (Satz: F = Tr M^2, "
          "Cauchy-Schwarz). Ohne Torus gilt die Eindeutigkeit nicht: Blockkanäle mit Raten, Spur- und Zustandserhaltung, "
          "J-Realität und Petz-Symmetrie übertreffen F_koh (Delta F = 1/54, 1/20, 22/375 für m = 3, 4, 5; spiegelinvariant "
          "noch bei m = 4), nur "
          "spurerhaltend erreicht man 0,881 gegen 0,75. Die native Dilatation ist nicht der kohärente Kanal: L^nat ist ein "
          "Fockoperator vom Rang 1 mit beiden Torusgewichten, nicht spurerhaltend, und die gewichtsaufgelöste Dilatation hat "
          "die nativen Raten genau an der Tomita-Wurzel t = 1/2, wo sie Phi_1/2 ist; die diagonalen nativen Daten (Raten, "
          "Casimir, KMS-Gewichte) sind in der ganzen Klasse gleich. Das Petz-Prinzip ersetzt also die Hälfte 'Krausrang 2' "
          "von P_koh durch ein Extremalprinzip, die Torusverträglichkeit des Schrittkanals bleibt Prämisse: P_koh = P_Torus + "
          "P_Petz, bedingt. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
say(f"\n   (Laufzeit {time.time() - T0:.1f} s)")
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z119_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "bilanz": bilanz, "urteil": urteil}, indent=1,
    ensure_ascii=False, default=str) + "\n")
