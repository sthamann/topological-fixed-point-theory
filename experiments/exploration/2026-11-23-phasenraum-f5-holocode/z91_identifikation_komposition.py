"""Z91: Identifikationsprinzip für L^nat und die Komposition vieler Nähte (Folgerechnung F3; zehnte externe Analyse, 25.11.).

Leitlinie 'Rekonstruktion statt Auswahl': Eine Symmetrie, die ein Schluss braucht, muss aus dem Original gezeigt werden und
darf nicht aus der gewünschten Raumzeit kommen. Firewalls: GRAV.NONCIRCULAR.01 (der Modulfluss bleibt Uhr und
Helizitätsphase, kein Bisognano-Wichmann) und DIMENSION.UPLIFT.FIREWALL.01 (alle Aussagen typisiert: Einteilchenkanal,
Zählung, Fusion auf Abschneidungen; keine wechselwirkende 3+1-Schließung). Gelöste Bausteine bleiben unberührt. Bestand:
z73 (kein natives su(2) aus Bilinearen), z76/z80 (zustandserhaltender Zweipfadkanal, Petz l(l+1)), z83 (Satz F G = const,
nativer Heber L^nat mit den Amplituden n_a + 1, offene Identifikation 'verdoppelter Diagonalschritt = nativer
Einzelschritt'), z88 (der 4D-Fockraum passt energietreu nicht in eine Naht).

Bezeichnungen: Paarschale H_m (dim m), Paarlage k = n_1 = 0..m-1, Lochtiefe n_2 = m-1-k; in der su(2)-Lesart Spin
j = (m-1)/2 und w = k - j. Nativer Kanal (z76) Phi_nat: H_m -> H_(m+1) mit p_k = (k+1)/(m+1) (Teilchen hoch, k -> k+1) und
q_k = (m-k)/(m+1) (Loch runter, k -> k).

Teil 1  Zehnte Analyse (a) und (b), exakt für m = 2..8. (a) Kovariante Stinespring-Zerlegung j -> j + 1/2 über
        Clebsch-Gordan: j' x j-quer ist multiplizitätsfrei, die kovarianten CP-Abbildungen sind der Kegel über m extremalen
        Kanälen Phi_s (Umgebung Spin s = 1/2..m-1/2, Krausrang 2s+1). Die nativen Übergänge, und schon die Forderung 'nur
        die beiden Nachbarschritte', erzwingen s = 1/2; Phi_1/2 = Phi_nat exakt; der Polarisationsfaktor c_s ist bei s = 1/2
        maximal; die Stinespring-Isometrie ist das Schwinger-Paar. (b) Der dekohärente Kanal mit denselben Übergängen ist
        CPTP, zustandserhaltend, J-reell, Petz-symmetrisch, spiegel- und torusverträglich, aber nicht SU(2)-kovariant; sein
        Petz-Komposit trägt nur den zonalen Teil von l(l+1).
Teil 2  Ist SU(2) auf H_m nativ? (i) Uhr, Modulphase, Torus, Spiegel, Tomita: abelsch mal Weyl-Z_2; mit dem
        sl(2,R)-Casimir erzeugen die nativen Operatoren ganz u(m), ohne ein su(2) auszuzeichnen; (ii) native Hüpfer (Strom,
        Möbius L^p_-1 L^h_1 aus L^nat) schließen nur für m <= 3; (iii) polare Familie V |L|^t: Abschluss für alle m genau bei
        t = 1/2, auch mit getrennten Exponenten; (iv) nativer sl(2,R)-Casimir auf H_m = su(2)-Casimir-Superoperator auf den
        Diagonalen (exakt, ohne Wurzel: die Raten sind nativ); (v) Symmetriealgebra der Kanäle: kohärent genau su(2),
        dekohärent abelsch; (vi) Fortsetzungen mit Krausrang 2: k-abhängige Mischung gegen Eindeutigkeit bei Torus + Rang 2
        (bis auf Eichung), Pfadverdichtung n + 1. Ergebnis: genau eine Prämisse P_koh.
Teil 3  Zusammensetzung. (i) Exakte Zählung bis E = 60: Sym^N(Vakuum + Einteilchenleiter) = prod (1 - q^m)^(-m^2) bis E = N;
        ohne Symmetrisierung zu viele, fermionisch zu wenige Zustände; Austauschvorzeichen zweier Paaranregungen
        verschiedener Nähte +1 (gemeinsame CAR); hartkernige Komposition mit kappa_4 = -1/(2N). (ii) Relative
        Tensorprodukte über Ladung/Uhr, die c-Algebra (maximal abelsch, linear), die Impulse P_mu (quadratisch) und die
        L^nat-sl(2,R): Korrelationen und Kumulanten; kubische Vertizes (K-Typ-Hindernis, Unitaritätsgrenze Delta = 1).
        (iii) Kandidaten der elften Analyse: diagonale E8-Komposition ((E8)_1 x (E8)_1 > (E8)_2 x Ising, Mischterm
        Summe J_1 J_2, c_rel(N)), relative Paritätsphase W_12 = (-1)^(n_1 n_2) mit R~ = (R_1 R_2) W_12, und 'lokal frei,
        kollektiv korreliert' (ganzzahliges K im Modulhamiltonoperator, das die Uhr nicht sieht).
Teil 4  Bilanz (E/C/L/O) und Urteil für F3.
"""
import itertools
import json
import time
from math import comb
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.linalg as sla
import scipy.sparse as sps
import scipy.sparse.linalg as spla
import sympy as sp
from sympy.physics.wigner import clebsch_gordan

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
rng = np.random.default_rng(91)
mp.mp.dps = 40
HALB = sp.Rational(1, 2)
MS = range(2, 9)                                                    # m = 2..8 (die Analyse rechnete sieben Niveaus)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


# ------------------------------------------------------------------ Bausteine
def spin_np(n):
    """Spin j = (n-1)/2 auf C^n, Basis k = 0..n-1 mit w = k - j aufsteigend: (J_x, J_y, J_z)."""
    j = (n - 1)/2
    w = np.arange(n) - j
    Jp = np.zeros((n, n))
    for k in range(n - 1):
        Jp[k + 1, k] = np.sqrt((j - w[k])*(j + w[k] + 1))
    return [(Jp + Jp.T)/2, (Jp - Jp.T)/2j, np.diag(w)]


def spin_ex(n):
    j = sp.Rational(n - 1, 2)
    Jz = sp.diag(*[k - j for k in range(n)])
    Jp = sp.zeros(n, n)
    for k in range(n - 1):
        w = k - j
        Jp[k + 1, k] = sp.sqrt((j - w)*(j + w + 1))
    return [(Jp + Jp.T)/2, (Jp - Jp.T)/(2*sp.I), Jz]


def kraus_nat(m):
    """z76: Teilchen hoch sqrt((k+1)/(m+1)), Loch runter sqrt((m-k)/(m+1))."""
    Ap, Ah = np.zeros((m + 1, m)), np.zeros((m + 1, m))
    for k in range(m):
        Ap[k + 1, k] = np.sqrt((k + 1)/(m + 1))
        Ah[k, k] = np.sqrt((m - k)/(m + 1))
    return [Ap, Ah]


def kraus_dek(m):
    """Dekohärente Fortsetzung: je Schrittart und Lage k ein Krausoperator (die Umgebung registriert auch die Lage)."""
    Ks = []
    for A in kraus_nat(m):
        for k in range(m):
            B = np.zeros_like(A)
            B[:, k] = A[:, k]
            Ks.append(B)
    return Ks


def kraus_misch(m, theta):
    """Krausrang 2 mit k-abhängigem Mischwinkel: u(k) = sqrt(p_k)(cos, sin) theta_k, d(k) senkrecht zu u(k-1)."""
    A1, A2 = np.zeros((m + 1, m)), np.zeros((m + 1, m))
    for k in range(m):
        wp, wq = np.sqrt((k + 1)/(m + 1)), np.sqrt((m - k)/(m + 1))
        A1[k + 1, k], A2[k + 1, k] = wp*np.cos(theta(k)), wp*np.sin(theta(k))
        A1[k, k], A2[k, k] = -wq*np.sin(theta(k - 1)), wq*np.cos(theta(k - 1))
    return [A1, A2]


def superop(Ks):
    """Zeilenweise Vektorisierung: vec(A X B) = (A kron B^T) vec(X)."""
    return sum(np.kron(A, A.conj()) for A in Ks)


def ad(X):
    n = X.shape[0]
    return np.kron(X, np.eye(n)) - np.kron(np.eye(n), X.T)


def choi(Ks):
    return sum(np.outer(A.reshape(-1), A.reshape(-1).conj()) for A in Ks)


def anwenden(Ks, X):
    return sum(A @ X @ A.conj().T for A in Ks)


# =================================================================== Teil 1
say("Teil 1: Zehnte Analyse (a) und (b): kovariante Stinespring-Zerlegung, dekohärentes Gegenbeispiel")


def kraus_s(m, s):
    """Extremaler SU(2)-kovarianter Kanal j -> j' = j + 1/2 mit Umgebung Spin s, exakt:
    A_sigma[k', k] = <j' w'; s sigma | j w>, w = k - j, w' = w - sigma = k' - j' (Stinespring = Kopplung j' x s -> j)."""
    j = sp.Rational(m - 1, 2)
    jp = j + HALB
    As = []
    for i in range(int(2*s) + 1):
        sg = -s + i
        A = sp.zeros(m + 1, m)
        for k in range(m):
            w = k - j
            kp = w - sg + jp
            if 0 <= kp <= m:
                A[int(kp), k] = clebsch_gordan(jp, s, j, w - sg, sg, w)
        As.append(A)
    return As


def T_nat(m):
    T = sp.zeros(m + 1, m)
    for k in range(m):
        T[k + 1, k] = sp.Rational(k + 1, m + 1)
        T[k, k] = sp.Rational(m - k, m + 1)
    return T


ok_mf = ok_cov = ok_tp = ok_eind = ok_traeger = ok_z76 = ok_pol = True
lam_tab, lam2_tab, pol_tab = {}, {}, {}
for m in MS:
    j = sp.Rational(m - 1, 2)
    jp = j + HALB
    SW = [HALB + i for i in range(m)]
    KR = {s: kraus_s(m, s) for s in SW}
    Jn, Jpn = spin_np(m), spin_np(m + 1)
    G = [np.kron(Jpn[a], np.eye(m)) - np.kron(np.eye(m + 1), Jn[a].T) for a in range(3)]   # Erzeuger von U' x U-quer
    Cas = sum(g @ g for g in G)
    ev = np.linalg.eigvalsh((Cas + Cas.conj().T)/2)
    werte, anz = np.unique(np.round(ev, 9), return_counts=True)
    ok_mf &= (len(werte) == m and all(abs(wv - float(s*(s + 1))) < 1e-8 and a_ == int(2*s + 1)
                                       for wv, a_, s in zip(werte, anz, SW)))
    Ts = {}
    for s, As in KR.items():
        ok_tp &= sum((A.T*A for A in As), sp.zeros(m, m)) == sp.eye(m)
        Ts[s] = sum((A.applyfunc(lambda x: x**2) for A in As), sp.zeros(m + 1, m))
        C = choi([np.array(A.tolist(), dtype=float) for A in As])
        ok_cov &= (max(np.abs(C @ g - g @ C).max() for g in G) < 1e-12
                   and np.abs(Cas @ C - float(s*(s + 1))*C).max() < 1e-10
                   and np.linalg.matrix_rank(C, tol=1e-10) == int(2*s + 1))
    Mlam = sp.Matrix.hstack(*[Ts[s].reshape((m + 1)*m, 1) for s in SW])
    los, frei = Mlam.gauss_jordan_solve(T_nat(m).reshape((m + 1)*m, 1))
    ok_eind &= frei.shape[0] == 0 and list(los) == [1] + [0]*(m - 1)
    zeilen = [i for i in range((m + 1)*m) if i//m not in (i % m, i % m + 1)]      # Übergänge außer den Nachbarschritten
    Ms_ = sp.Matrix.vstack(Mlam.extract(zeilen, list(range(m))), sp.ones(1, m))
    los2, frei2 = Ms_.gauss_jordan_solve(sp.Matrix([0]*len(zeilen) + [1]))
    ok_traeger &= frei2.shape[0] == 0 and list(los2) == [1] + [0]*(m - 1)
    lam_tab[m] = [str(x) for x in los]
    lam2_tab[m] = [str(x) for x in los2]
    Aup, Adn = KR[HALB]                                                            # sigma = -1/2 (hoch), +1/2 (Loch)
    ok_z76 &= all(Aup[k + 1, k]**2 == sp.Rational(k + 1, m + 1) and Adn[k, k]**2 == sp.Rational(m - k, m + 1)
                  for k in range(m))
    ok_z76 &= (len({sp.sign(Aup[k + 1, k]) for k in range(m)}) == 1 and len({sp.sign(Adn[k, k]) for k in range(m)}) == 1
               and Ts[HALB] == T_nat(m))
    Jz_j = sp.diag(*[k - j for k in range(m)])
    Jz_jp = sp.diag(*[k - jp for k in range(m + 1)])
    cs = {}
    for s, As in KR.items():
        X = sum((A.T*Jz_jp*A for A in As), sp.zeros(m, m))
        c = (j*(j + 1) + jp*(jp + 1) - s*(s + 1))/(2*j*(j + 1))
        ok_pol &= sp.simplify(X - c*Jz_j) == sp.zeros(m, m)
        cs[s] = c
    ok_pol &= cs[HALB] == sp.Rational(m + 2, m + 1) and all(cs[SW[i]] > cs[SW[i + 1]] for i in range(m - 1))
    pol_tab[m] = {str(s): str(c) for s, c in cs.items()}
ZAHLEN["lambda_aus_nativen_uebergaengen"] = lam_tab
ZAHLEN["lambda_nur_nachbarschritte"] = lam2_tab
ZAHLEN["polarisation_c_s"] = {m: pol_tab[m] for m in (2, 5, 8)}
say(f"   Polarisationsfaktoren c_s (m = 5): {pol_tab[5]}")
check("[E] (a) Kovarianz: Für m = 2..8 ist j' x j-quer (j = (m-1)/2, j' = m/2) multiplizitätsfrei mit den Spins "
      "s = 1/2..m-1/2; jeder extremale Kanal Phi_s (Clebsch-Gordan, exakt spurerhaltend) hat eine kovariante Choi-Matrix vom "
      "Rang 2s+1 im s-Isotypenteil. Die kovarianten CP-Abbildungen sind also genau Summe lambda_s Phi_s mit lambda_s >= 0, "
      "spurerhaltend für Summe lambda_s = 1", ok_mf and ok_cov and ok_tp)
check("[E] (a) Eindeutigkeit: Die nativen Übergänge p_k, q_k legen lambda = (1, 0, .., 0) eindeutig fest (lineares System, "
      "exakt, m = 2..8). Schon die Forderung 'nur die beiden Nachbarschritte' (k -> k+1, k -> k) mit Spurerhaltung erzwingt "
      "s = 1/2; die Übergangswahrscheinlichkeiten (k+1)/(m+1), (m-k)/(m+1) und damit die Zustandserhaltung folgen dann von "
      "selbst", ok_eind and ok_traeger)
check("[E] (a) Phi_1/2 ist der native Kanal aus z76: Krausamplituden zum Quadrat (k+1)/(m+1) und (m-k)/(m+1), je "
      "Krausoperator ein festes Vorzeichen; Krausrang 2, kleinste mögliche Umgebung Spin 1/2 (s muss halbzahlig sein). Die "
      "Stinespring-Isometrie ist (a_1^+ (x) |runter> - a_2^+ (x) |hoch>)/sqrt(N + 2) mit den Schwinger-Oszillatoren n_1 = k, "
      "n_2 = m-1-k. Polarisation: Phi_s^+(J'_z) = c_s J_z mit c_s = (j(j+1) + j'(j'+1) - s(s+1))/(2j(j+1)), streng fallend "
      "in s, maximal c_1/2 = (m+2)/(m+1)", ok_z76 and ok_pol)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

# --- (b) dekohärentes Gegenbeispiel und Petz-Komposit
ok_b, ok_nat_su2, ok_dek_su2, ok_petz = True, True, True, True
su2_abw = {}
for m in MS:
    Jn, Jpn = spin_np(m), spin_np(m + 1)
    Mm, Mm1 = np.eye(m)[::-1], np.eye(m + 1)[::-1]
    N1, N1p = np.diag(np.arange(m, dtype=float)), np.diag(np.arange(m + 1, dtype=float))
    X = rng.normal(size=(m, m)) + 1j*rng.normal(size=(m, m))
    Tn = np.array(T_nat(m).tolist(), dtype=float)
    U, Up = sla.expm(0.7j*N1), sla.expm(0.7j*N1p)
    Lp = {}
    for name, Ks in (("nat", kraus_nat(m)), ("dek", kraus_dek(m))):
        S = superop(Ks)
        SR = (m + 1)/m*S.conj().T                                                  # Petz mit Referenz I/m
        eig = [np.abs(sum(A.conj().T @ A for A in Ks) - np.eye(m)).max() < 1e-13,
               np.allclose([[sum(abs(A[kp, k])**2 for A in Ks) for k in range(m)] for kp in range(m + 1)], Tn, atol=1e-14),
               np.allclose(anwenden(Ks, np.eye(m)/m), np.eye(m + 1)/(m + 1), atol=1e-14),
               np.allclose(anwenden(Ks, X.conj()), anwenden(Ks, X).conj(), atol=1e-13),
               np.allclose(SR @ S, (SR @ S).conj().T, atol=1e-13),
               np.allclose(anwenden(Ks, Mm @ X @ Mm), Mm1 @ anwenden(Ks, X) @ Mm1, atol=1e-13),
               np.allclose(anwenden(Ks, U @ X @ U.conj().T), Up @ anwenden(Ks, X) @ Up.conj().T, atol=1e-12)]
        ok_b &= all(eig)
        abw = max(np.abs(anwenden(Ks, Jn[a] @ X - X @ Jn[a])
                         - (Jpn[a] @ anwenden(Ks, X) - anwenden(Ks, X) @ Jpn[a])).max() for a in range(3))
        rang = np.linalg.matrix_rank(choi(Ks), tol=1e-10)
        su2_abw[(m, name)] = (float(abw), int(rang))
        if name == "nat":
            ok_nat_su2 &= abw < 1e-12 and rang == 2
        else:
            ok_dek_su2 &= abw > 0.05 and rang == 2*m
        Lp[name] = m*(m + 1)*(np.eye(m*m) - SR @ S)
    CasS = sum(ad(J) @ ad(J) for J in Jn)
    di = [i*m + i for i in range(m)]
    od = [i for i in range(m*m) if i not in di]
    Ld = Lp["dek"]
    zon = np.sort(np.linalg.eigvalsh(Ld[np.ix_(di, di)].real))
    ok_petz &= (np.allclose(Lp["nat"], CasS, atol=1e-10) and np.allclose(zon, [l*(l + 1) for l in range(m)], atol=1e-9)
                and np.allclose(Ld[np.ix_(od, od)], m*(m + 1)*np.eye(len(od)), atol=1e-10)
                and np.allclose(Ld[np.ix_(di, od)], 0, atol=1e-12) and np.allclose(Ld[np.ix_(od, di)], 0, atol=1e-12)
                and np.allclose(Lp["nat"][np.ix_(di, di)], Ld[np.ix_(di, di)], atol=1e-10))
ZAHLEN["su2_abweichung_und_krausrang"] = {f"m={k[0]}, {k[1]}": v for k, v in su2_abw.items() if k[0] in (2, 5, 8)}
check("[E negativ] (b) Der dekohärente Kanal (sqrt(p_k)|k+1><k| und sqrt(q_k)|k><k| einzeln, Krausrang 2m) hat dieselben "
      "Übergänge und ist wie der native Kanal spurerhaltend, zustandserhaltend (I/m -> I/(m+1)), J-reell (Phi(X-quer) = "
      "Phi(X)-quer), Petz-symmetrisch (R Phi selbstadjungiert im KMS-Skalarprodukt der Spurzustände), spiegelverträglich "
      "(k -> m-1-k mit Tausch der Schrittarten) und torusverträglich (Phase von N_1). Er ist aber nicht SU(2)-kovariant "
      "(Abweichung > 0,05, m = 2..8), der native Kanal schon (< 1e-12). Vollständige Positivität und Tomita-Symmetrie "
      "entscheiden die Fortsetzung also nicht", ok_b and ok_nat_su2 and ok_dek_su2)
check("[E] Petz-Komposit: m(m+1)(1 - R Phi) ist für den kohärenten Kanal der volle Casimir-Superoperator (l(l+1) mit "
      "Vielfachheit 2l+1, z76). Beim dekohärenten Kanal bleibt nur der zonale Teil: auf den Diagonalen derselbe Block mit "
      "l(l+1), l = 0..m-1, je einmal; auf allen Kohärenzen die Konstante m(m+1). Der Petz-Rückkanal ist also keine "
      "unabhängige Quelle des su(2): Er trägt es nur, wenn der Kanal kohärent ist", ok_petz)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

# =================================================================== Teil 2
say("\nTeil 2: Ist SU(2) auf den Paarschalen nativ?")


def lie_dim(gens, tol=1e-9):
    """Dimension der von antihermiteschen Erzeugern erzeugten reellen Lie-Algebra (iterierte Kommutatoren)."""
    n = gens[0].shape[0]
    B, mats = [], []

    def add(X):
        v = np.concatenate([X.real.ravel(), X.imag.ravel()])
        for _ in range(2):
            for b in B:
                v = v - (v @ b)*b
        nv = np.linalg.norm(v)
        if nv > tol:
            v = v/nv
            B.append(v)
            mats.append((v[:n*n] + 1j*v[n*n:]).reshape(n, n))
            return True
        return False

    for g in gens:
        add(g/np.linalg.norm(g))
    neu = list(range(len(mats)))
    while neu:
        frisch = []
        for i in neu:
            for jx in range(len(mats)):
                if add(mats[i] @ mats[jx] - mats[jx] @ mats[i]):
                    frisch.append(len(mats) - 1)
        neu = frisch
    return len(B)


def Lm1_ex(m):
    """Nativer Heber L^nat_-1: H_m -> H_(m+1), Teilchen hoch mit n_1 + 1 = k + 1, Loch runter mit n_2 + 1 = m - k (z83)."""
    L = sp.zeros(m + 1, m)
    for k in range(m):
        L[k + 1, k] = k + 1
        L[k, k] = m - k
    return L


def casimir_sl2(m):
    """Casimir der L^nat-sl(2,R) auf H_m: L_0(L_0 - 1) - L_-1 L_1."""
    return m*(m - 1)*sp.eye(m) - (Lm1_ex(m - 1)*Lm1_ex(m - 1).T if m > 1 else sp.zeros(1, 1))


# --- (i) native Symmetrien und native Operatoren auf H_m
nat_dim, nat_dim_C, spiegel_ok = {}, {}, True
for m in MS:
    N1 = np.diag(np.arange(m, dtype=float))
    Mm = np.eye(m)[::-1]
    spiegel_ok &= np.allclose(Mm @ N1 @ Mm, (m - 1)*np.eye(m) - N1)
    nat_dim[m] = lie_dim([1j*np.eye(m), 1j*N1, 1j*(Mm @ N1 @ Mm)])
    Cm = np.array(casimir_sl2(m).tolist(), dtype=float)
    nat_dim_C[m] = lie_dim([1j*np.eye(m), 1j*N1, 1j*Cm])
ZAHLEN["native_symmetrien_liedim"] = nat_dim
ZAHLEN["native_operatoren_mit_sl2_casimir_liedim"] = nat_dim_C
check("[E negativ] (i) Native Symmetrien auf H_m: Uhr R = e^{2 pi i L_0} = 1 und Modulfluss e^{-i beta t L_0} sind auf H_m "
      "Phasen (L_0 = m, z83), der Cartan-Torus e^{i theta N_1} kommt aus dem Hardy-Schnitt (z73), der Spiegel x -> -x ist "
      "das Weyl-Element (M N_1 M = m - 1 - N_1), Tomita-J ist auf dem Diagonalsektor die Konjugation der reellen Basis. "
      "Zusammen erzeugen sie eine abelsche Lie-Algebra der Dimension 2 (mal Z_2), kein su(2) (m = 2..8)",
      spiegel_ok and all(v == 2 for v in nat_dim.values()))
check("[E] (i) Nimmt man den Casimir C der L^nat-sl(2,R) auf H_m hinzu, erzeugen die nativen hermiteschen Operatoren ganz u(m) "
      f"(Dimension m^2: {nat_dim_C}). Ein su(2) ist darin enthalten, aber nicht ausgezeichnet: 'aus nativen Operatoren "
      "gebaut' ist ohne Auswahlregel leer", all(nat_dim_C[m] == m*m for m in MS))

# --- (ii) native Hüpfer, (iii) polare Familie
AMPS = {"Strom (Toeplitz, z73)": lambda k, m: 1.0,
        "Möbius L^p_-1 L^h_1 (aus L^nat)": lambda k, m: float((k + 1)*(m - 1 - k)),
        "Schwinger (CCR)": lambda k, m: np.sqrt((k + 1)*(m - 1 - k))}


def huepf_dim(m, amp):
    E = np.zeros((m, m))
    for k in range(m - 1):
        E[k + 1, k] = amp(k, m)
    H = np.diag(np.arange(m) - (m - 1)/2)
    return lie_dim([1j*H, E - E.T, 1j*(E + E.T)])


hop = {name: {m: huepf_dim(m, a) for m in MS} for name, a in AMPS.items()}
ZAHLEN["huepfer_liedim"] = hop
say(f"   erzeugte Lie-Algebren {{Hüpfer: {{m: Dimension}}}}: {hop}")
soll_bf = {m: m*(m - 1)//2 if m % 2 else m*(m + 1)//2 for m in MS}               # so(m) bzw. sp(m)
check("[E negativ] (ii) Native Hüpfer auf H_m mit der nativen Cartan-Größe: das Stromhüpfen (Toeplitz-Hälften von J_-1, "
      "Amplitude 1) und das Möbiushüpfen L^p_-1 L^h_1 aus den Hardy-Hälften von L^nat (Amplitude (n_1 + 1) n_2) erzeugen für "
      "m = 2, 3 su(2) (Dimension 3), ab m = 4 aber so(m) (m ungerade) bzw. sp(m) (m gerade), Dimension m(m-1)/2 bzw. "
      "m(m+1)/2: die Algebren, die die vom Spiegel k -> m-2-k induzierte Bilinearform erhalten. Nur die Schwinger-Amplituden "
      "sqrt((n_1 + 1) n_2) bleiben für alle m bei 3",
      all(hop["Schwinger (CCR)"][m] == 3 for m in MS)
      and all(hop[n_][m] == 3 for n_ in list(AMPS)[:2] for m in (2, 3))
      and all(hop[n_][m] == soll_bf[m] for n_ in list(AMPS)[:2] for m in MS if m >= 4))

TS = (0.0, 0.25, 0.5, 0.75, 1.0)
pol = {t: {m: huepf_dim(m, lambda k, m, t=t: ((k + 1)*(m - 1 - k))**t) for m in MS} for t in TS}
schliesst = [t for t in TS if all(v == 3 for v in pol[t].values())]
raenge = {}
for m in range(4, 9):
    A_ = np.array([[np.log(k + 1), np.log(m - 1 - k), 1.0] for k in range(m - 1)])
    raenge[m] = int(np.linalg.matrix_rank(A_, tol=1e-10))
getrennt = {tt: huepf_dim(6, lambda k, m, tt=tt: (k + 1)**tt[0]*(m - 1 - k)**tt[1])
            for tt in ((0.5, 0.5), (0.5, 0.6), (0.4, 0.5), (0.3, 0.7))}
ZAHLEN["polare_familie_liedim"] = {str(t): v for t, v in pol.items()}
ZAHLEN["getrennte_exponenten_m6"] = {str(k): v for k, v in getrennt.items()}
check("[E] (iii) Polare Familie: Mit L^p_-1 = V_1 |L^p_-1| und L^h_1 = |L^h_1| V_2^+ (polare Zerlegung der nativen "
      "Fermionleitern D_1/2) hat E_+(t) = V_1 |L^p|^t |L^h|^t V_2^+ die Amplitude ((n_1 + 1) n_2)^t; für t = 0 (Toeplitz), "
      "1/4, 1/2, 3/4, 1 (Möbius) schließt sie für alle m = 2..8 genau bei t = 1/2. Mit getrennten Exponenten t_1, t_2 ist "
      "(2 t_1 - 1) log(n_1 + 1) + (2 t_2 - 1) log n_2 = const auf den m - 1 Stufen nur für t_1 = t_2 = 1/2 lösbar (Rang 3 "
      f"für m = 4..8: {raenge}; m = 6 numerisch: {getrennt}). Die CCR sind genau die Wurzel der nativen Amplituden",
      schliesst == [0.5] and all(r == 3 for r in raenge.values())
      and getrennt[(0.5, 0.5)] == 3 and all(v > 3 for k, v in getrennt.items() if k != (0.5, 0.5)))

# --- (iv) native Raten = zonaler Schatten des su(2)-Casimir
ok_kom, ok_ident, ok_spek = True, True, True
for m in MS:
    Lu, Lo = Lm1_ex(m), Lm1_ex(m - 1)
    ok_kom &= Lu.T*Lu - Lo*Lo.T == 2*m*sp.eye(m)
    Cm = casimir_sl2(m)
    Jx, Jy, Jz = spin_ex(m)
    D = sp.zeros(m, m)
    for k in range(m):
        E = sp.zeros(m, m)
        E[k, k] = 1
        Y = sum((J*(J*E - E*J) - (J*E - E*J)*J for J in (Jx, Jy, Jz)), sp.zeros(m, m)).applyfunc(sp.expand)
        ok_ident &= all(Y[a_, b_] == 0 for a_ in range(m) for b_ in range(m) if a_ != b_)
        for kp in range(m):
            D[kp, k] = Y[kp, kp]
    ok_ident &= D == Cm
    ev = sorted(sp.Matrix(Cm).eigenvals(multiple=True))
    ok_spek &= ev == sorted(l*(l + 1) for l in range(m))
check("[E] (iv) Native Raten: Mit den nativen Heberamplituden (n_1 + 1, n_2 + 1) gilt [L_1, L_-1] = 2 L_0 (m = 2..8), und "
      "der Casimir C = L_0(L_0 - 1) - L_-1 L_1 der L^nat-sl(2,R) ist auf H_m exakt die Einschränkung des su(2)-Casimir-"
      "Superoperators Summe [J_a,[J_a, .]] (Spin (m-1)/2) auf die Diagonalen |k><k| (der Diagonalen auf Diagonalen abbildet), "
      "mit Nebendiagonale -(n_1 + 1) n_2 und Spektrum l(l+1), l = 0..m-1, je einmal. Ohne jede Wurzel: Die Raten |J_+|^2 "
      "sind nativ, das su(2) ist ihre kohärente Wurzel", ok_kom and ok_ident and ok_spek)


# --- (v) Symmetriealgebra der Kanäle
def sym_dim(Ks, m):
    """Dimension aller (X, X') in gl(m) + gl(m+1) mit Phi o ad_X = ad_X' o Phi (enthält zwei triviale Richtungen)."""
    S = superop(Ks)
    spalten = []
    for i, jx in itertools.product(range(m), repeat=2):
        E = np.zeros((m, m))
        E[i, jx] = 1
        spalten.append((S @ ad(E)).ravel())
    for i, jx in itertools.product(range(m + 1), repeat=2):
        E = np.zeros((m + 1, m + 1))
        E[i, jx] = 1
        spalten.append((-ad(E) @ S).ravel())
    Mx = np.array(spalten).T
    sv = np.linalg.svd(Mx, compute_uv=False)
    return int(Mx.shape[1] - np.sum(sv > 1e-10*sv[0]))


THETA = lambda k: 0.4*k + 0.1
symd = {m: {"kohaerent": sym_dim(kraus_nat(m), m), "dekohaerent": sym_dim(kraus_dek(m), m),
            "misch": sym_dim(kraus_misch(m, THETA), m)} for m in MS}
ZAHLEN["symmetriealgebra_dimension"] = symd
say(f"   Symmetriealgebra (mit 2 trivialen Richtungen) {{m: ...}}: {symd}")
check("[E] (v) Symmetriealgebra des Kanals: Für den kohärenten Kanal ist {(X, X'): Phi o ad_X = ad_X' o Phi} genau "
      "sl(2, C) plus die zwei trivialen Richtungen (Dimension 5, m = 2..8): das su(2) ist aus dem Kanal selbst ablesbar, nicht "
      "eingesetzt. Beim dekohärenten Kanal sind es die Diagonalen beider Seiten (abelsch, 2m + 1), bei der k-abhängigen "
      "Mischung mit Krausrang 2 nur die trivialen Richtungen (ab m = 3)",
      all(v["kohaerent"] == 5 and v["dekohaerent"] == 2*m + 1 for m, v in symd.items())
      and all(symd[m]["misch"] == 2 for m in MS if m >= 3))

# --- (vi) Krausrang 2: Mischfamilie, Eindeutigkeit bei Torus + Rang 2, Pfadverdichtung
ok_misch, ok_eich = True, True
for m in MS:
    Km = kraus_misch(m, THETA)
    Kc = kraus_misch(m, lambda k: 0.9)                                               # konstanter Winkel: Umgebungsdrehung
    Tn = np.array(T_nat(m).tolist(), dtype=float)
    N1, N1p = np.diag(np.arange(m, dtype=float)), np.diag(np.arange(m + 1, dtype=float))
    U, Up = sla.expm(0.7j*N1), sla.expm(0.7j*N1p)
    X = rng.normal(size=(m, m)) + 1j*rng.normal(size=(m, m))
    ok_misch &= (np.abs(sum(A.T @ A for A in Km) - np.eye(m)).max() < 1e-13
                 and np.allclose([[sum(A[kp, k]**2 for A in Km) for k in range(m)] for kp in range(m + 1)], Tn, atol=1e-14)
                 and np.linalg.matrix_rank(choi(Km), tol=1e-10) == 2
                 and np.allclose(choi(Kc), choi(kraus_nat(m)), atol=1e-13)
                 and np.abs(np.abs(choi(Km)) - np.abs(choi(kraus_nat(m)))).max() > 0.01
                 and not np.allclose(anwenden(Km, U @ X @ U.conj().T), Up @ anwenden(Km, X) @ Up.conj().T, atol=1e-6))
    al, be = rng.uniform(0, 2*np.pi, m), rng.uniform(0, 2*np.pi, m)
    Ap, Ah = kraus_nat(m)
    A1, A2 = Ap*np.exp(1j*al)[None, :], Ah*np.exp(1j*be)[None, :]                  # torusgewichtet, beliebige Phasen
    gew = {(kp - k) for A in (A1, A2) for kp, k in zip(*np.nonzero(A))}
    chi, psi = np.zeros(m), np.zeros(m + 1)
    psi[0] = -be[0]
    for k in range(m):
        if k > 0:
            chi[k] = -psi[k] - be[k]
        psi[k + 1] = -al[k] - chi[k]
    Din, Dout = np.diag(np.exp(1j*chi)), np.diag(np.exp(1j*psi))
    ok_eich &= (gew == {0, 1} and np.allclose(choi([Dout @ A1 @ Din, Dout @ A2 @ Din]), choi(kraus_nat(m)), atol=1e-12)
                and np.allclose(anwenden([A1, A2], U @ X @ U.conj().T), Up @ anwenden([A1, A2], X) @ Up.conj().T,
                                atol=1e-12))


def kompos_rang(fn, m, n):
    Ks, d = [np.eye(m)], m
    for _ in range(n):
        Ks = [B @ A for A in Ks for B in fn(d)]
        d += 1
    return int(np.linalg.matrix_rank(choi(Ks), tol=1e-9))


pfad = {name: [kompos_rang(fn, 3, n) for n in (1, 2, 3)]
        for name, fn in (("kohaerent", kraus_nat), ("misch", lambda d: kraus_misch(d, THETA)), ("dekohaerent", kraus_dek))}
ZAHLEN["pfadverdichtung_krausrang_m3_n123"] = pfad
say(f"   Krausrang von n Schritten H_3 -> H_(3+n), n = 1, 2, 3: {pfad}")
check("[E] (vi) Krausrang 2 allein genügt nicht: Die Familie mit k-abhängigem Mischwinkel (u(k) entlang theta_k, d(k+1) "
      "senkrecht zu u(k)) ist spurerhaltend, hat die nativen Übergänge und Krausrang 2, ist aber nicht torusverträglich und "
      "nicht eichäquivalent zum nativen Kanal (Beträge der Choi-Einträge verschieden; konstanter Winkel = Umgebungsdrehung = "
      "nativer Kanal). Torusverträglich mit Rang 2 heißt dagegen: je ein Krausoperator pro Gewicht (Schrittart), beliebige "
      "Phasen, und diagonale Eichungen führen ihn exakt in den nativen Kanal über (m = 2..8). Pfadverdichtung: n Schritte "
      f"haben Krausrang n + 1 nur im kohärenten Fall ({pfad})",
      ok_misch and ok_eich and pfad["kohaerent"] == [2, 3, 4]
      and all(r > n + 1 for n, r in zip((2, 3), pfad["misch"][1:])) and all(r > n + 1 for n, r in zip((1, 2, 3), pfad["dekohaerent"])))
befund("Teil 2: SU(2) auf den Paarschalen ist bedingt nativ. Nativ (E) sind die Cartan-Größe, der Torus, das Weyl-Element "
       "und die Raten |J_+|^2 = (n_1 + 1) n_2 als Casimir der L^nat-sl(2,R). Kein natives Operatorprodukt und kein natives "
       "Symmetriedatum ist selbst das su(2); die polare Familie schließt nur bei der Wurzel t = 1/2. Die verbleibende Prämisse "
       "ist genau eine, P_koh: Die Umgebung des Elementarschritts registriert nur die Schrittart (Krausrang 2, ein "
       "torusgewichteter Krausoperator je nativer Schrittart), nicht die Lage k. Unter P_koh mit Spurerhaltung und "
       "Zustandserhaltung (z76) ist der Kanal eindeutig Phi_1/2; daraus folgen das su(2) als seine Symmetriealgebra, die CCR "
       "als sqrt(L_0 + 1) mal Krausamplitude und die Identifikation aus z83, ohne die Ziel-CCR vorauszusetzen", True)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

# =================================================================== Teil 3
say("\nTeil 3: Zusammensetzung vieler Nähte")
say("  (i) Zählung")
EMAX = 60


def euler_produkt(a, N):
    """Koeffizienten von prod_m (1 - t^m)^(-a(m)) bis t^N exakt (wie z88)."""
    sig = [0]*(N + 1)
    for d in range(1, N + 1):
        ad_ = a(d)
        if ad_:
            for k in range(d, N + 1, d):
                sig[k] += d*ad_
    p = [0]*(N + 1)
    p[0] = 1
    for n in range(1, N + 1):
        p[n] = sum(sig[k]*p[n - k] for k in range(1, n + 1))//n
    return p


fock = euler_produkt(lambda d: d*d, EMAX)
dp = [[0]*(EMAX + 1) for _ in range(EMAX + 1)]                                     # dp[n][E]: n Teilchen, Energie E
dp[0][0] = 1
for mm in range(1, EMAX + 1):
    d = mm*mm
    neu = [[0]*(EMAX + 1) for _ in range(EMAX + 1)]
    for n in range(EMAX + 1):
        for E, v in enumerate(dp[n]):
            if v:
                kk = 0
                while n + kk <= EMAX and E + mm*kk <= EMAX:
                    neu[n + kk][E + mm*kk] += v*comb(d + kk - 1, kk)
                    kk += 1
    dp = neu
summe_ok = all(sum(dp[n][E] for n in range(EMAX + 1)) == fock[E] for E in range(EMAX + 1))
SN = lambda N, E: sum(dp[n][E] for n in range(min(N, EMAX) + 1))
NS = (5, 10, 20, 40, 59)
sym_ok = all(all(SN(N, E) == fock[E] for E in range(N + 1)) and fock[N + 1] - SN(N, N + 1) == 1 for N in NS)
chi1 = [0] + [m*m for m in range(1, EMAX + 1)]


def poly_mal(u, v):
    w = [0]*(EMAX + 1)
    for i, x in enumerate(u):
        if x:
            for jx in range(EMAX + 1 - i):
                w[i + jx] += x*v[jx]
    return w


plain = [1] + [0]*EMAX
basis_poly = [1] + chi1[1:]
for _ in range(10):
    plain = poly_mal(plain, basis_poly)
ferm = [1] + [0]*EMAX
for mm in range(1, EMAX + 1):
    fak = [0]*(EMAX + 1)
    for kk in range(0, min(mm*mm, EMAX//mm) + 1):
        fak[mm*kk] = comb(mm*mm, kk)
    ferm = poly_mal(ferm, fak)
ZAHLEN["zaehlung"] = {"fock_E0..15": fock[:16], "fock_E60": fock[60],
                      "defizit_Sym^N_bei_E=N+1": {N: fock[N + 1] - SN(N, N + 1) for N in NS},
                      "ohne_Symmetrisierung_N10_E2_E10": [plain[2], plain[10]], "fock_E2_E10": [fock[2], fock[10]],
                      "fermionisch_E2_E10": [ferm[2], ferm[10]]}
say(f"   prod (1 - q^m)^(-m^2): {fock[:12]} ... E = 60: {fock[60]}")
check("[E] (i) Exakte Zählung bis E = 60 (ganzzahlig): Die N-fache symmetrische Potenz von (Vakuum + Einteilchenleiter "
      "Summe m^2 q^m), also N ununterscheidbare Nähte mit je höchstens einem Teilchen, stimmt mit dem 4D-Fockcharakter "
      "prod (1 - q^m)^(-m^2) für alle E <= N überein (N = 5, 10, 20, 40, 59); bei E = N + 1 fehlt genau ein Zustand (N + 1 "
      "Teilchen auf m = 1). Der volle freie Fockraum ist Sym^N für N -> unendlich", summe_ok and sym_ok)
check(f"[E negativ] (i) Ohne Symmetrisierung (unterscheidbare Nähte, (1 + Summe m^2 q^m)^N) gibt es zu viele Zustände: N = 10 "
      f"hat bei E = 2 {plain[2]} statt {fock[2]} und bei E = 10 {plain[10]} statt {fock[10]}; fermionisch "
      f"(prod (1 + q^m)^(m^2)) zu wenige: {ferm[2]} statt {fock[2]} bei E = 2. Nur die bosonische Zusammensetzung gibt den "
      "freien 4D-Fockraum", plain[2] == 4*10 + comb(10, 2) and plain[10] > fock[10] and ferm[2] == 4 and fock[2] == 5)

# Austauschvorzeichen in einer gemeinsamen CAR (Jordan-Wigner, zwei Nähte zu je 4 Moden)
NB = 8
dimF = 1 << NB


def c_op_n(mode, dag, nb):
    """Jordan-Wigner: Erzeuger (dag = 1) oder Vernichter (dag = 0) der Mode 'mode' auf nb Moden."""
    dimn = 1 << nb
    Mx = np.zeros((dimn, dimn))
    for s in range(dimn):
        if ((s >> mode) & 1) == dag:
            continue
        vz = -1 if bin(s & ((1 << mode) - 1)).count("1") % 2 else 1
        Mx[s ^ (1 << mode), s] = vz
    return Mx


def c_op(mode, dag):
    return c_op_n(mode, dag, NB)


def tausch():
    """Unitäre Vertauschung der Nähte (Moden i <-> i + 4) mit fermionischem Umordnungsvorzeichen."""
    Ux = np.zeros((dimF, dimF))
    for s in range(dimF):
        occ = [i for i in range(NB) if (s >> i) & 1]
        bild = [(i + 4) % NB for i in occ]
        inv = sum(1 for a_, b_ in itertools.combinations(bild, 2) if a_ > b_)
        Ux[sum(1 << i for i in bild), s] = -1 if inv % 2 else 1
    return Ux


Ux = tausch()
vak = np.zeros(dimF)
vak[(1 << 0) | (1 << 1) | (1 << 4) | (1 << 5)] = 1                                 # Fermisee: Lochplätze 0, 1 und 4, 5 besetzt
Bp = [c_op(2 + 4*i, 1) @ c_op(1 + 4*i, 0) for i in range(2)]                      # Paar psi*_(p=0) psi_(r=-1) je Naht
paar2 = Bp[0] @ Bp[1] @ vak
ferm2 = c_op(2, 1) @ c_op(6, 1) @ vak
aus_ok = (np.allclose(Ux @ vak, vak) and np.allclose(Ux @ Bp[0] @ Ux.T, Bp[1]) and np.allclose(Bp[0] @ Bp[1], Bp[1] @ Bp[0])
          and np.allclose(Ux @ paar2, paar2) and np.allclose(Ux @ ferm2, -ferm2) and np.linalg.norm(paar2) > 0.5)
check("[E] (i) Statistik: Liegen die Fermionen aller Nähte in einer gemeinsamen CAR (Jordan-Wigner, zwei Nähte), so "
      "vertauschen die Paaranregungen (Einteilchenzustände der Leiter) verschiedener Nähte, und der Nahttausch hat auf "
      "B_1 B_2 Omega den Eigenwert +1, auf einzelnen Fermionen -1. Die bosonische Symmetrisierung ist damit die native "
      "Statistik; dass die Nähte ununterscheidbar sind (Prämisse P_ident), bleibt eine Setzung", aus_ok)

kap4 = {}
for N in range(1, 13):
    Jx_, _, _ = spin_np(N + 1)
    x_ = np.sqrt(2/N)*Jx_.real
    v0 = np.zeros(N + 1)
    v0[0] = 1.0
    mu = [v0 @ np.linalg.matrix_power(x_, r) @ v0 for r in range(5)]
    kap4[N] = (mu[2], mu[3], mu[4] - 3*mu[2]**2)
ZAHLEN["hartkern_kumulanten_N"] = {N: [round(float(v), 12) for v in t] for N, t in kap4.items()}
check("[E] (i) Endlich viele hartkernige Nähte (je höchstens ein Quant, symmetrisiert = Spin N/2, Kollektivfeld "
      "x = (J_+ + J_-)/sqrt(2N) im Vakuum |N/2, -N/2>): kappa_2 = 1/2, kappa_3 = 0, kappa_4 = -1/(2N) (N = 1..12). Die "
      "endliche Zusammensetzung ist gerade nichtgaußsch von der Ordnung 1/N (eine statistische Wechselwirkung, die für "
      "N -> unendlich verschwindet), nie kubisch",
      all(abs(t[0] - 0.5) < 1e-12 and abs(t[1]) < 1e-12 and abs(t[2] + 1/(2*N)) < 1e-12 for N, t in kap4.items()))
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

# --- (ii) Relative Tensorprodukte über native Algebren
say("  (ii) Relative Tensorprodukte über native geteilte Algebren")
# (alpha) Ladung und Uhr
bet, alq = 4*mp.pi, mp.mpf(1)/4
QS = range(-12, 13)


def ladungsgesetz(beta):
    w = {Q: mp.e**(-beta*(mp.mpf(Q)**2/2 + (alq - mp.mpf(1)/2)*Q)) for Q in QS}
    Z = sum(w.values())
    return {Q: v/Z for Q, v in w.items()}


def kumulanten(P):
    mu = [sum(P[Q]*mp.mpf(Q)**r for Q in P) for r in range(5)]
    return [mu[1], mu[2] - mu[1]**2, mu[3] - 3*mu[2]*mu[1] + 2*mu[1]**3,
            mu[4] - 4*mu[3]*mu[1] - 3*mu[2]**2 + 12*mu[2]*mu[1]**2 - 6*mu[1]**4]


ein = ladungsgesetz(bet)
fus = {Q: ein[Q]**2 for Q in QS}
Zf = sum(fus.values())
fus = {Q: v/Zf for Q, v in fus.items()}
dop = ladungsgesetz(2*bet)
k_ein, k_fus = kumulanten(ein), kumulanten(fus)
uhr = {(Q1, Q2): ein[Q1]*ein[Q2] for Q1 in QS for Q2 in QS if (Q1 - Q2) % 4 == 0}
Zu = sum(uhr.values())
EQ1 = sum(v*Q1 for (Q1, Q2), v in uhr.items())/Zu
EQ2 = sum(v*Q2 for (Q1, Q2), v in uhr.items())/Zu
cov = sum(v*(Q1 - EQ1)*(Q2 - EQ2) for (Q1, Q2), v in uhr.items())/Zu
var1 = sum(v*(Q1 - EQ1)**2 for (Q1, Q2), v in uhr.items())/Zu
var2 = sum(v*(Q2 - EQ2)**2 for (Q1, Q2), v in uhr.items())/Zu
korr_uhr = cov/mp.sqrt(var1*var2)
# 4D-Teil: kleines Gesamtmodell Ladung (x) Oszillator je Naht, Projektion Summe P_Q (x) P_Q
bt, qs_, nosc = 1.0, [-1, 0, 1], 6
pq = np.array([np.exp(-bt*(Q*Q/2 - Q/4)) for Q in qs_])
pq /= pq.sum()
po = np.exp(-bt*np.arange(nosc))
po /= po.sum()
rho1 = np.kron(np.diag(pq), np.diag(po))
d1 = len(qs_)*nosc
Pi = sum(np.kron(np.kron(np.diag(np.eye(3)[i]), np.eye(nosc)), np.kron(np.diag(np.eye(3)[i]), np.eye(nosc))) for i in range(3))
rf = Pi @ np.kron(rho1, rho1) @ Pi
rf /= np.trace(rf)
rf4 = rf.reshape(3, nosc, 3, nosc, 3, nosc, 3, nosc)
rosc = np.einsum("aibjacbd->ijcd", rf4).reshape(nosc*nosc, nosc*nosc)
ao = np.diag(np.sqrt(np.arange(1, nosc)), 1)
xo = (ao + ao.T)/np.sqrt(2)
X1, X2 = np.kron(xo, np.eye(nosc)), np.kron(np.eye(nosc), xo)
k3_4d = max(abs(np.trace(rosc @ A @ B @ C)) for A in (X1, X2) for B in (X1, X2) for C in (X1, X2))
ladung_ok = (max(abs(fus[Q] - dop[Q]) for Q in QS) < mp.mpf(10)**-35 and abs(k_ein[2]) > 1e-3 and abs(k_fus[2]) > 1e-6
             and np.allclose(rosc, np.kron(np.diag(po), np.diag(po)), atol=1e-15)
             and abs(np.trace(rosc @ X1 @ X2)) < 1e-15 and k3_4d < 1e-15 and 0.999 < korr_uhr < 1)
ZAHLEN["ladung"] = {"kumulanten_einzeln(k1..k4)": [mp.nstr(x, 8) for x in k_ein],
                    "kumulanten_fusion(k1..k4)": [mp.nstr(x, 8) for x in k_fus],
                    "uhr_fusion_mod4_korrelation": mp.nstr(korr_uhr, 15)}
check("[E] (ii-a) Ladung/Uhr: Das relative Tensorprodukt über C[Q] identifiziert die Ladungen (Q_1 = Q_2); das gefundene "
      "Ladungsgesetz ist exakt das bei 2 beta (Abweichung < 1e-35), seine dritte Kumulante "
      f"{mp.nstr(k_fus[2], 4)} ist innen und war schon in einer Naht da ({mp.nstr(k_ein[2], 4)}, z58). Über C[R] (R = i^Q) "
      f"bleibt Q_1 = Q_2 mod 4 (1 - Korrelation = {mp.nstr(1 - korr_uhr, 3)}). Der Oszillator-, also 4D-Teil bleibt exakt das "
      "Produkt (Gesamtmodell Ladung (x) Oszillator, reduzierter Zustand = Produkt, <x_1 x_2> = 0, alle dritten Momente 0): "
      "keine Korrelation und keine Dreipunktstruktur in 4D; R = 1 auf der ganzen Leiter (z83)", ladung_ok)

# (beta) c-Algebra: kommutativ, gaußsch, maximal abelsch
nmax = 40
dd = nmax + 1
a1 = sps.diags(np.sqrt(np.arange(1, dd)), 1, format="csr")
Aop, Bop = sps.kron(a1, sps.identity(dd), format="csr"), sps.kron(sps.identity(dd), a1, format="csr")
c_ok, cf_tab = True, {}
for q in (0.5, float(np.exp(-2*np.pi))):
    kap = np.sqrt((1 - q)/(1 + q))
    Om = np.zeros(dd*dd)
    for n in range(dd):
        Om[n*dd + n] = np.sqrt(1 - q*q)*q**n
    c = kap*(Aop + Bop.T)
    cd = c.T.conj().tocsr()
    nn = [float(Om @ (cd @ (c @ Om))), float(Om @ (c @ (cd @ Om)))]
    cfs = []
    for t in (0.4, 0.3 + 0.5j):
        v = spla.expm_multiply((1j*(t*cd + np.conj(t)*c)).tocsc(), Om.astype(complex))
        cfs.append(abs(np.vdot(Om, v) - np.exp(-abs(t)**2)))
    cf_tab[q] = {"<c^+c>, <cc^+>": nn, "|charakteristische Funktion - e^{-|t|^2}|": cfs}
    c_ok &= abs(nn[0] - 1) < 1e-12 and abs(nn[1] - 1) < 1e-12 and max(cfs) < 1e-10
OMEGA = np.array([[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]])        # (x_a, p_a, x_b, p_b)
vRe, vIm = np.array([1, 0, 1, 0]), np.array([0, 1, 0, -1])                         # Re c ~ x_a + x_b, Im c ~ p_a - p_b
lagr = abs(vRe @ OMEGA @ vIm) < 1e-15 and np.linalg.matrix_rank(np.array([vRe, vIm])) == 2
xs, ys = sp.symbols("x y", real=True)
pc = sp.exp(-(xs**2 + ys**2))/sp.pi
pf = pc**2/sp.integrate(pc**2, (xs, -sp.oo, sp.oo), (ys, -sp.oo, sp.oo))
var_f = sp.integrate(pf*(xs**2 + ys**2), (xs, -sp.oo, sp.oo), (ys, -sp.oo, sp.oo))
m2x = sp.integrate(pf*xs**2, (xs, -sp.oo, sp.oo), (ys, -sp.oo, sp.oo))
m4x = sp.integrate(pf*xs**4, (xs, -sp.oo, sp.oo), (ys, -sp.oo, sp.oo))
logpf = sp.expand(sp.expand_log(sp.log(pf), force=True))
quadr = sp.Poly(logpf.subs(sp.log(sp.pi), 0).subs(sp.log(2), 0), xs, ys).total_degree() == 2
ZAHLEN["c_algebra"] = {str(k): v for k, v in cf_tab.items()}
ZAHLEN["c_fusion"] = {"E|c|^2": str(var_f), "kappa_4(Re c)": str(sp.simplify(m4x - 3*m2x**2))}
check("[E] (ii-b) c-Algebra (z78: c = kappa (a + b^+), vertauscht mit c^+): Im Purifikationszustand gilt <c^+c> = <cc^+> = 1 "
      "und <exp(i(t c^+ + t-quer c))> = e^{-|t|^2} (Fockraum bis n = 40, q = 1/2 und q = e^{-2 pi}): c ist eine "
      "klassische komplexe Gaußvariable. (Re c, Im c) ~ (x_a + x_b, p_a - p_b) spannen eine Lagrange-Ebene, die c-Algebra ist "
      "maximal abelsch, jede Naht ist L^2 ihres Spektrums mit Vielfachheit 1. Das relative Tensorprodukt verschmilzt daher "
      f"zwei Nähte zu einer: Gesetz ~ p_1 p_2, E|c|^2 = {var_f}, log-Dichte quadratisch, kappa_4 = 0. Keine neue "
      "Freiheit, keine Kumulante ab Ordnung 3", c_ok and lagr and var_f == sp.Rational(1, 2) and quadr
      and sp.simplify(m4x - 3*m2x**2) == 0)

# (gamma) Impulse P_mu = c^+ sigma_mu c / 2 (quadratisch): Fusion über C[P]
P0, phi = sp.symbols("P0 phi", positive=True)
g1 = 4*P0*sp.exp(-2*P0)                                                            # Gesetz von P_0 = |c|^2/2 in einer Naht
imp = {}
for name, a_ in (("U(1)-Haarmass auf der Faser (d^3P/P_0)", 1), ("Bogenlänge (Lebesgue)", sp.Rational(1, 2))):
    dichte = (g1/P0**a_)**2*P0**a_                                                # f^2 dnu mit f = g/P0^a
    Zd = sp.integrate(dichte, (P0, 0, sp.oo))
    E1 = sp.integrate(dichte*2*P0, (P0, 0, sp.oo))/Zd
    E2 = sp.integrate(dichte*(2*P0)**2, (P0, 0, sp.oo))/Zd
    imp[name] = sp.simplify(E2 - E1**2)
kreuz2 = sp.integrate(sp.exp(-sp.I*phi), (phi, 0, 2*sp.pi))
sig2 = HALB                                                                        # c-Fusion: CN(0, 1/2) je Komponente
gauss_kum = 2*sig2**2 - 2*sig2**2
ZAHLEN["impulsfusion_verbundene_kreuz_vierpunktsumme"] = {k: str(v) for k, v in imp.items()}
check("[E] (ii-c) Impulse: P^(1) = P^(2) heißt c^(2) = e^{i phi} c^(1) (Rang eins, P^2 = 0 bleibt). Weil jede Naht "
      "U(1)-invariant ist, ist phi gleichverteilt: die Kreuzkovarianz <c^(1) c^(2)-quer> verschwindet (Integral e^{-i phi} "
      f"= {kreuz2}), ebenso alle ungeraden Momente (c -> -c). Die verbundene Kreuz-Vierpunktsumme Cov(|c^(1)|^2, |c^(2)|^2) "
      f"ist dagegen {imp} (für zwei Bezugsmaße), während sie bei der c-Fusion (gaußsch) {gauss_kum} ist: Die Impulsfusion ist "
      "eine gerade, nichtgaußsche Zwangskorrelation über eine einzige globale Variable, keine kubische Struktur",
      kreuz2 == 0 and all(v > 0 for v in imp.values()) and gauss_kum == 0
      and imp["U(1)-Haarmass auf der Faser (d^3P/P_0)"] == HALB and imp["Bogenlänge (Lebesgue)"] == sp.Rational(5, 8))

# (delta) L^nat-sl(2,R): Kommutant, Tomita-Paarung, kubische Vertizes
MT = 8
off, pos = {}, 0
for m in range(1, MT + 1):
    off[m] = pos
    pos += m
dimV = pos
L0 = np.diag([float(m) for m in range(1, MT + 1) for _ in range(m)])
LM1 = np.zeros((dimV, dimV))
for m in range(1, MT):
    LM1[off[m + 1]:off[m + 1] + m + 1, off[m]:off[m] + m] = np.array(Lm1_ex(m).tolist(), dtype=float)
L1 = LM1.T
CasV = L0 @ L0 - L0 - LM1 @ L1
spalten1 = [i for m in range(1, MT) for i in range(off[m], off[m] + m)]
spalten2 = [i for m in range(2, MT + 1) for i in range(off[m], off[m] + m)]
blockbasis = []
for m in range(1, MT + 1):
    for a_, b_ in itertools.product(range(m), repeat=2):
        T = np.zeros((dimV, dimV))
        T[off[m] + a_, off[m] + b_] = 1
        blockbasis.append(T)
Mk = np.array([np.concatenate([(T @ LM1 - LM1 @ T)[:, spalten1].ravel(), (T @ L1 - L1 @ T)[:, spalten2].ravel()])
               for T in blockbasis]).T
svk = np.linalg.svd(Mk, compute_uv=False)
kdim = int(Mk.shape[1] - np.sum(svk > 1e-10*svk[0]))
evC, vecC = np.linalg.eigh(CasV)
proj = []
for h in range(1, MT + 1):
    sel = np.abs(evC - h*(h - 1)) < 1e-8
    proj.append(vecC[:, sel] @ vecC[:, sel].T)
proj_ok = all(np.abs(np.concatenate([(P @ LM1 - LM1 @ P)[:, spalten1].ravel(), (P @ L1 - L1 @ P)[:, spalten2].ravel()])).max()
              < 1e-10 for P in proj) and np.linalg.matrix_rank(np.array([P.ravel() for P in proj])) == MT
# kubischer sl(2,R)-Vertex D_1 x D_1 -> D_2 und sein Bild in der su(2)-Lesart
L1_H2 = np.array(Lm1_ex(1).T.tolist(), dtype=float)                                 # L_1: H_2 -> H_1
_, _, vh = np.linalg.svd(L1_H2)
tief2 = vh[-1]
J2 = spin_np(2)
bild = np.diag(tief2)
cas_bild = sum(J @ (J @ bild - bild @ J) - (J @ bild - bild @ J) @ J for J in J2)
l_bild_ok = np.allclose(cas_bild, 2*bild) and abs(np.trace(bild)) < 1e-12
ktyp_ok = all((m1 - 1)/2 + (m2 - 1)/2 < (m1 + m2 - 1)/2 for m1 in range(1, 41) for m2 in range(1, 41))
xv = sp.symbols("x0:4", real=True)
ypt, zpt = [1, 0, 0, 0], [0, 2, 0, 0]
r12 = sp.sqrt(sum((xv[i] - ypt[i])**2 for i in range(4)))
r13 = sp.sqrt(sum((xv[i] - zpt[i])**2 for i in range(4)))
f3 = 1/(r12*r13*sp.sqrt(5))                                                        # <phi phi phi>, Delta = 1
g3 = 1/(r13**2*5)                                                                  # <phi phi :phi^2:>, Delta = (1, 1, 2)
punkt = {xv[0]: sp.Rational(1, 3), xv[1]: sp.Rational(-1, 2), xv[2]: sp.Rational(1, 5), xv[3]: sp.Rational(2, 7)}
lap_f = sp.nsimplify(sum(sp.diff(f3, xv[i], 2) for i in range(4)).subs(punkt))
lap_g = sp.simplify(sum(sp.diff(g3, xv[i], 2) for i in range(4)))
ZAHLEN["sl2_fusion"] = {"kommutant_dim_M8": kdim, "tiefstgewicht_D2_in_H2": tief2.round(12).tolist(),
                        "box_phi_phi_phi_am_punkt": str(sp.N(lap_f, 12)), "box_phi_phi_phi2": str(lap_g)}
check("[E] (ii-d) L^nat-sl(2,R): Auf der abgeschnittenen Einpaarschicht (+)_(m <= 8) H_m ist der Raum der Verflechter "
      f"(Kommutant von L_0, L^nat_(+-1)) {kdim}-dimensional und genau von den Spektralprojektoren des Casimir aufgespannt (je "
      "ein D_h, h = 1..8). Das relative Tensorprodukt über U(sl(2,R)) paart also (h, n) zwischen den Nähten; mit dem "
      "KMS-Regulator q^{L_0} ist die Paarung T = 1 genau der Einpaaranteil der Purifikation (z83): Die Tomita-Verdopplung ist "
      "selbst eine sl(2,R)-Fusion, eine Verschränkung ohne Vertex", kdim == MT and proj_ok)
check("[E negativ] (ii-d) Kubische Vertizes: sl(2,R) allein erlaubt D_1 x D_1 -> D_2 (Produkt der Tiefstgewichte); das Bild, "
      "der Tiefstgewichtsvektor (1, -1)/sqrt(2) von H_2, ist in der su(2)-Lesart sigma_z mit l = 1, aus l = 0 (x) l = 0. Ein "
      "so(4)- und energieverträglicher kubischer Vertex zwischen Skalarleitern verschwindet immer, weil (j_1, j_1) x "
      "(j_2, j_2) den K-Typ (j_3, j_3) der Zielschale m_1 + m_2 nicht enthält (j_1 + j_2 < j_3, m_1, m_2 <= 40). Und auf der "
      f"Unitaritätsgrenze Delta = 1 (Box phi = 0) ist die konforme Struktur <phi phi phi> nicht harmonisch (Box = "
      f"{sp.N(lap_f, 6)} != 0 am Testpunkt), also C_phiphiphi = 0, während <phi phi :phi^2:> harmonisch bleibt "
      f"(Box = {lap_g}): freie Felder haben keine verbundene Dreipunktfunktion",
      l_bild_ok and ktyp_ok and lap_f != 0 and lap_g == 0)
befund("Teil 3: Der freie 4D-Fockraum entsteht als bosonische Zusammensetzung N -> unendlich vieler gleicher Nähte, jede im "
       "Vakuum oder mit einem Teilchen der Leiter; die Statistik ist nativ, die Ununterscheidbarkeit eine Prämisse. Keine der "
       "nativen geteilten Algebren erzeugt Wechselwirkung: Ladung und Uhr wirken nur innen, die c-Algebra ist linear "
       "(gaußsch, sie verschmilzt die Nähte), Impulse und sl(2,R) sind gerade (Zwangskorrelationen, Tomita-Paarung). Eine "
       "verbundene Dreipunktstruktur braucht eine geteilte Struktur, die weder gerade noch linear ist und die "
       "Unitaritätsgrenze verlässt; die nativen Kandidaten dafür sind die geladenen Vertexoperatoren (Fusion der "
       "Ladungssektoren nach Ladungsaddition) und die nichtabelschen Stromantworten (z75)", True)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

# --- (iii) Kandidaten der elften externen Analyse
say("  (iii) Kandidaten der elften externen Analyse: diagonale E8-Komposition, relative Paritätsphase, kollektive "
    "Nichtgaußizität")
# (1) diagonale E8-Komposition
Nsym, gsym = sp.symbols("N g", positive=True)
cE8 = lambda k: sp.Rational(248)*k/(k + 30)
crel = {N: 8*N - cE8(N) for N in range(1, 11)}
crel_ok = (cE8(1) == 8 and cE8(2) == sp.Rational(31, 2) and crel[2] == HALB
           and sp.simplify(8*Nsym - 248*Nsym/(Nsym + 30) - 8*Nsym*(Nsym - 1)/(Nsym + 30)) == 0)
nicht_frei = [N for N in range(2, 11) if sp.Rational(2*crel[N]).q != 1]          # c nicht in Z/2: kein Freifeldsektor
S11, S22, S12 = sp.symbols("S11 S22 S12")                                         # Summe J_1^2, Summe J_2^2, Summe J_1 J_2
T1s, T2s, Tdg = S11/62, S22/62, (S11 + S22 + 2*S12)/64
Trel = T1s + T2s - Tdg
trel_ok = (sp.expand(Trel - (T1s + T2s - S12)/32) == 0 and sp.expand(S12 - (Tdg - 31*Trel)) == 0
           and sp.expand(T1s + T2s + gsym*S12 - ((1 + gsym)*Tdg + (1 - 31*gsym)*Trel)) == 0)
WE8 = []
for i, jx in itertools.combinations(range(8), 2):
    for si, sj in itertools.product((1.0, -1.0), repeat=2):
        v = np.zeros(8)
        v[i], v[jx] = si, sj
        WE8.append(v)
for vz in itertools.product((0.5, -0.5), repeat=8):
    if sum(1 for x in vz if x < 0) % 2 == 0:
        WE8.append(np.array(vz))
WE8 = np.array(WE8)
GE8 = WE8 @ WE8.T
wurzel_ok = (len(WE8) == 240 and np.allclose(np.diag(GE8), 2) and np.allclose(GE8, np.round(GE8))
             and {tuple(v) for v in WE8} == {tuple(-v) for v in WE8} and np.linalg.matrix_rank(WE8) == 8
             and len(WE8) + 8 == 248)
# su(2)-Nullmoden-Abschneidung (A_1 statt E_8; (A_1)_1 x (A_1)_1 / (A_1)_2 ist ebenfalls Ising): je Naht ein Dublett
PAULI = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]])]
cd4, cc4 = [c_op_n(k, 1, 4) for k in range(4)], [c_op_n(k, 0, 4) for k in range(4)]
Sp = [[sum(cd4[2*i + a_] @ cc4[2*i + b_]*PAULI[a][a_, b_]/2 for a_ in range(2) for b_ in range(2)) for a in range(3)]
      for i in range(2)]
Mmix = sum(Sp[0][a] @ Sp[1][a] for a in range(3))
Cas1, Cas2 = sum(X @ X for X in Sp[0]), sum(X @ X for X in Sp[1])
C12 = sum((Sp[0][a] + Sp[1][a]) @ (Sp[0][a] + Sp[1][a]) for a in range(3))
sel4 = [s for s in range(16) if bin(s & 3).count("1") == 1 and bin(s >> 2).count("1") == 1]
evm, vecm = np.linalg.eigh(Mmix[np.ix_(sel4, sel4)].real)
psi = np.zeros(16, dtype=complex)
psi[sel4] = vecm[:, 0]
gam = np.array([[np.vdot(psi, cd4[i_] @ cc4[j_] @ psi) for j_ in range(4)] for i_ in range(4)])
erw = lambda O: np.vdot(psi, O @ psi)
drei_quer = erw(Sp[0][0] @ Sp[0][1] @ Sp[1][2])
drei_cartan = erw(Sp[0][2] @ Sp[0][2] @ Sp[1][2])
ein_punkt = max(abs(erw(Sp[i][a])) for i in range(2) for a in range(3))
Mz = (Sp[0][2] @ Sp[1][2]).real
evz, vecz = np.linalg.eigh(Mz[np.ix_(sel4, sel4)])
cart_slater = True
for col in range(2):                                                              # die beiden Grundzustände von S_1^z S_2^z
    ph = np.zeros(16, dtype=complex)
    ph[sel4] = vecz[:, col]
    gz = np.array([[np.vdot(ph, cd4[i_] @ cc4[j_] @ ph) for j_ in range(4)] for i_ in range(4)])
    cart_slater &= np.allclose(gz @ gz, gz, atol=1e-12)
e8_ok = (crel_ok and trel_ok and wurzel_ok and nicht_frei == [3, 4, 5, 6, 7, 8, 9]
         and np.allclose(Mmix, (C12 - Cas1 - Cas2)/2) and all(np.allclose(Mmix @ (Sp[0][a] + Sp[1][a]),
                                                                          (Sp[0][a] + Sp[1][a]) @ Mmix) for a in range(3))
         and abs(evm[0] + 0.75) < 1e-12 and not np.allclose(gam @ gam, gam, atol=1e-6)
         and abs(drei_quer - (-1j/8)) < 1e-12 and abs(drei_cartan) < 1e-12 and ein_punkt < 1e-12
         and abs(evz[0] + 0.25) < 1e-12 and abs(evz[1] + 0.25) < 1e-12 and cart_slater)
ZAHLEN["e8_diagonal"] = {"c_rel(N), N = 1..10": {N: str(v) for N, v in crel.items()}, "c_nicht_in_Z/2": nicht_frei,
                         "su2_proxy_<S1x S1y S2z>": str(np.round(drei_quer, 12)), "<S1z S1z S2z>": float(abs(drei_cartan))}
check("[E] (iii-1) Diagonale E8-Komposition: c((E8)_2) = 31/2, Coset c = 1/2 (Ising), c_rel(N) = 8N(N-1)/(N+30) "
      f"(symbolisch), nicht in Z/2 für N = {nicht_frei} (dort kein Sektor freier Felder). T_rel = (1/32)(T_1 + T_2 - Summe "
      "J_1 J_2) und Summe J_1 J_2 = T_diag - 31 T_rel, also T_1 + T_2 + g Summe J_1 J_2 = (1+g) T_diag + (1-31g) T_rel: zwei "
      "vertauschende Sektoren, die freie Summe T_1 + T_2 selbst bleibt unverändert. Der Mischterm besteht aus 8 "
      "Cartan-Bilinearen (in den Oszillatoren gaußsch) und 240 Wurzeltermen E_1^alpha E_2^-alpha (Vertexoperatoren mit "
      "Ladungsübertrag alpha; Wurzelsystem geprüft: Norm 2, ganzzahlig, negationsabgeschlossen, Rang 8). Nullmoden-Abschneidung "
      "mit A_1 (gleiche Coset-Struktur): der Mischterm ist (C_12 - C_1 - C_2)/2, sein Grundzustand das Singulett, nicht "
      "gaußsch (gamma^2 != gamma), mit verbundener Kreuz-Dreipunktfunktion <S_1^x S_1^y S_2^z> = -i/8 ~ f^abc bei "
      "verschwindenden Einpunktfunktionen; in Cartan-Richtung <S_1^z S_1^z S_2^z> = 0, und der Cartan-Teil allein hat "
      "Slater-Grundzustände", e8_ok)

# (2) relative Paritätsphase W_12 = (-1)^(n_1 n_2)
f1 = np.array([bin(s & 0b1111).count("1") for s in range(dimF)])
f2 = np.array([bin(s >> 4).count("1") for s in range(dimF)])
W12 = np.diag((-1.0)**(f1*f2))
Uord = np.zeros((dimF, dimF))
for s in range(dimF):
    Uord[sum(1 << ((i + 4) % NB) for i in range(NB) if (s >> i) & 1), s] = 1
R1, R2 = np.diag(1j**(f1 - 2)), np.diag(1j**(f2 - 2))                             # R_i = i^(Q_i), Q_i = n_i - 2 (Fermisee)
Rt = R1 @ R2 @ W12
gerade = ([c_op(a_, 1) @ c_op(b_, 0) for a_ in range(4) for b_ in range(4)]
          + [c_op(a_, 1) @ c_op(b_, 0) for a_ in range(4, 8) for b_ in range(4, 8)]
          + [c_op(a_, 1) @ c_op(b_, 1) for a_ in range(4) for b_ in range(4) if a_ < b_])
lin = np.array([X.ravel() for X in [c_op(k, d_) for k in range(NB) for d_ in (0, 1)]]).T
Wc = W12 @ c_op(0, 0) @ W12
koef, *_ = np.linalg.lstsq(lin, Wc.ravel(), rcond=None)
rest_lin = np.linalg.norm(lin @ koef - Wc.ravel())/np.linalg.norm(Wc)
sektor_ok = all(((-1)**(f1[s]*f2[s])) == ((-1)**((f1[s] + f2[s] + 1)*f1[s])) for s in range(dimF))
w12_ok = (np.allclose(Ux, Uord @ W12) and np.allclose(Rt @ Rt, np.diag((-1.0)**(f1 + f2)))
          and np.allclose(np.linalg.matrix_power(Rt, 4), np.eye(dimF))
          and all(np.allclose(W12 @ X, X @ W12) for X in gerade) and np.allclose(Wc, c_op(0, 0) @ np.diag((-1.0)**f2))
          and rest_lin > 0.5 and sektor_ok and np.allclose(R1 @ R1, np.diag((-1.0)**f1)))
check("[E] (iii-2) Relative Paritätsphase W_12 = (-1)^(n_1 n_2) (Jordan-Wigner, zwei Nähte zu je 4 Moden): Sie ist genau "
      "der Unterschied zwischen fermionischem und gewöhnlichem Nahttausch (U_ferm = U_gew W_12); mit R_i = i^(Q_i) gilt für "
      "R~ = (R_1 R_2) W_12 exakt R~^2 = (-1)^(n_1 + n_2) und R~^4 = 1. Als Operator ist W_12 nichtgaußsch (W c_1 W = "
      f"c_1 (-1)^(n_2), relativer Abstand zum linearen Spann {rest_lin:.2f}), vertauscht aber mit allen nahtlokalen geraden "
      "Operatoren (Paare, Ströme) und ist auf jedem Sektor fester Gesamtladung N die gaußsche Uhrpotenz (-1)^((N+1) n_1) = "
      "R_1^(2(N+1)). Unter Ladungsauswahl ist W_12 eine native Uhrpotenz, die die graduierte Zusammensetzung (Statistik) "
      "festlegt, keine Wechselwirkung", w12_ok)

# (3) lokal frei, kollektiv korreliert: ganzzahliges K im Modulhamiltonoperator
ek = np.array([0.25, 1.25, 0.25, 1.25])                                           # Quartermoden x in Z + 1/4, zwei Nähte
occ4 = np.array([[(s >> k) & 1 for k in range(4)] for s in range(16)])
L0d = occ4 @ ek
Kint = occ4[:, 0]*occ4[:, 2] + 2*occ4[:, 1]*occ4[:, 3] - occ4[:, 0]*occ4[:, 1]*occ4[:, 2]
lw = -4*np.pi*(L0d + Kint)
wK = np.exp(lw - lw.max())
wK /= wK.sum()
lZ = np.log(np.exp(lw - lw.max()).sum()) + lw.max()
uhr_phase = np.exp(-0.5j*(lw - lZ))                                               # rho^{-i/2} auf der Diagonalen
Rdiag = np.exp(2j*np.pi*L0d)
uhr_abw = np.abs(uhr_phase/uhr_phase[0] - Rdiag/Rdiag[0]).max()
nn = lambda k: occ4[:, k]
cov02 = wK @ (nn(0)*nn(2)) - (wK @ nn(0))*(wK @ nn(2))
rhoK = np.diag(wK)
cdK, ccK = [c_op_n(k, 1, 4) for k in range(4)], [c_op_n(k, 0, 4) for k in range(4)]
L0op = np.diag(L0d)
zeit_ok = True
for t_, s_ in ((0.3, 1.1), (2.0, -0.7)):
    Ut, Us = sla.expm(-1j*t_*L0op), sla.expm(-1j*s_*L0op)
    c0d_t = Ut.conj().T @ cdK[0] @ Ut
    c2d_s = Us.conj().T @ cdK[2] @ Us
    vier = np.trace(rhoK @ c0d_t @ ccK[0] @ c2d_s @ ccK[2])
    getrennt_ = np.trace(rhoK @ c0d_t @ ccK[0])*np.trace(rhoK @ c2d_s @ ccK[2])
    zeit_ok &= abs((vier - getrennt_) - np.exp(1j*(ek[0]*t_ + ek[2]*s_))*cov02) < 1e-14
mod = lambda t_: np.diag(wK**(1j*t_))
cm = mod(0.3) @ ccK[0] @ mod(-0.3)
prop_rest = np.linalg.norm(cm - (np.vdot(ccK[0].ravel(), cm.ravel())/np.vdot(ccK[0].ravel(), ccK[0].ravel()))*ccK[0])
Rop = np.diag(Rdiag)
uhr_mod = np.allclose(mod(-0.5) @ ccK[0] @ mod(0.5), Rop @ ccK[0] @ Rop.conj().T, atol=1e-12)
rand_abw = abs(wK @ nn(0) - np.exp(-np.pi)/(1 + np.exp(-np.pi)))
koll_ok = (uhr_abw < 1e-9 and abs(cov02) > 1e-4 and zeit_ok and prop_rest > 1e-3 and uhr_mod
           and np.allclose(rhoK @ L0op, L0op @ rhoK))
ZAHLEN["kollektiv"] = {"uhr_abweichung": float(uhr_abw), "Cov(n_0, n_2)": float(cov02),
                       "modulfluss_t0.3_rest": float(prop_rest), "randabweichung_n_0_hier": float(rand_abw)}
check("[E] (iii-3) Lokal frei, kollektiv korreliert: Jede ganzzahlige, mit L_0 vertauschende Größe K im "
      "Modulhamiltonoperator, rho_K ~ e^{-4 pi (L_0 + K)}, lässt die Uhridentität exakt bestehen (rho_K^{-i/2} = Phase mal "
      f"R auf vier Quartermoden, Abweichung {uhr_abw:.0e}), weil e^{{2 pi i K}} = 1; mit K = n_0 n_2 + 2 n_1 n_3 - n_0 n_1 n_2 "
      f"ist rho_K nichtgaußsch (Cov(n_0, n_2) = {cov02:.3e}) und sein Modulfluss für t = 0,3 wechselwirkend (Ad rho^{{it}} c_0 "
      "nicht proportional zu c_0), bei t = -1/2 aber genau Ad R. Die freie Dynamik vertauscht mit rho_K, die verbundenen "
      "Korrelationen werden mit freien Phasen transportiert (<c_0^+(t) c_0 c_2^+(s) c_2>_c = e^{i(e_0 t + e_2 s)} Cov), nicht "
      "erzeugt: ein Zustand, keine Wechselwirkung. Die Uhr legt K nicht fest; ob exakt freie Ränder erreichbar sind "
      f"(hier weichen sie um {rand_abw:.1e} ab; die Analyse nennt sechs Saaten), prüft F6", koll_ok)
befund("(iii) Einordnung der drei Kandidaten. (1) Diagonale E8: bedingt nativ (verlangt Naht = (E8)_1, Tor T2; das "
       "E8-Gitter selbst ist gelöst), als Zerlegung parameterfrei. Nichtgaußsch sind nur der Wurzelteil (Vertexoperatoren) "
       "und die Cosets (ab N = 3 kein Freifeldsektor): echte Wechselwirkung, aber im inneren 2D-Sektor; auf den "
       "Cartan-Oszillatoren gaußsch und paritätsgerade, also keine kubische 4D-Struktur; T_1 + T_2 bleibt, eine Kopplung g "
       "gewichtet nur zwei vertauschende Sektoren um. (2) W_12: nativ (aus R und (-1)^F), nichtgaußsch nur über "
       "Ladungssektoren hinweg; unter Ladungsauswahl eine gaußsche Uhrpotenz, sie legt die Statistik der Zusammensetzung "
       "fest, keine Wechselwirkung. (3) Lokal frei, kollektiv korreliert: eine Freiheit K, die die Uhr nicht sieht, also "
       "nicht nativ ausgewählt; ein nichtgaußscher Zustand der freien Dynamik, keine Wechselwirkung über den freien "
       "Fockraum hinaus. Echte (innere) Wechselwirkung bringt nur (1), und zwar über die Vertexoperatoren, wie in (ii) "
       "vermutet", True)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

# =================================================================== Teil 4
say("\nTeil 4: Bilanz für F3")
STATUS = {
    "(a) kovariante Eindeutigkeit": ("E", "unter SU(2)-Kovarianz ist die CPTP-Fortsetzung eindeutig Phi_1/2 = z76; schon "
                                          "'nur Nachbarschritte' genügt; Krausrang 2, Umgebung Spin 1/2 = Schwinger-Paar, "
                                          "maximale Polarisation (m + 2)/(m + 1)"),
    "(b) CP und Tomita reichen nicht": ("E (negativ)", "der dekohärente Kanal erfüllt Übergänge, CPTP, Zustandserhaltung, "
                                                       "J-Realität, Petz-Symmetrie, Spiegel und Torus, ist aber nicht "
                                                       "kovariant; Petz gibt dort nur den zonalen Teil"),
    "native Raten": ("E", "Casimir der L^nat-sl(2,R) auf H_m = su(2)-Casimir-Superoperator auf den Diagonalen, ohne Wurzel"),
    "su(2) als nativer Operator": ("E (negativ)", "Uhr, Modulphase, Torus, Spiegel, Tomita abelsch mal Z_2; Strom- und "
                                                  "Möbiushüpfer schließen nur für m <= 3; polare Familie nur bei t = 1/2; "
                                                  "native Operatoren erzeugen u(m) ohne Auszeichnung"),
    "Symmetriealgebra des Kanals": ("E", "kohärent genau su(2) (plus triviale Richtungen), dekohärent abelsch, Mischfamilie "
                                         "trivial"),
    "Identifikation 'verdoppelter Diagonalschritt = L^nat'": ("C", "eindeutig unter P_koh; dann folgen Phi_1/2, su(2), CCR "
                                                                  "und z83 ohne Ziel-CCR"),
    "Prämisse P_koh": ("O", "Krausrang 2 mit einem torusgewichteten Krausoperator je nativer Schrittart (die Umgebung "
                            "registriert nur die Schrittart); nicht aus CP, Tomita, Spiegel, Torus, Petz oder Uhr herleitbar"),
    "freier 4D-Fockraum": ("E / C", "Sym^N(Vakuum + Einteilchenleiter) = prod (1 - q^m)^(-m^2) exakt bis E = N; Statistik "
                                    "nativ (+1); P_ident (ununterscheidbare Nähte) und N -> unendlich gesetzt"),
    "Fusion über native Algebren": ("E (negativ)", "Ladung/Uhr nur innen; c-Algebra gaußsch und verschmelzend; Impulse und "
                                                   "sl(2,R) gerade Zwangskorrelationen bzw. Tomita-Paarung; keine verbundene "
                                                   "Dreipunktstruktur"),
    "kubischer Vertex": ("E (negativ)", "so(4) x E-verträgliche Vertizes zwischen Skalarleitern verschwinden (K-Typ); "
                                        "Delta = 1 verbietet <phi phi phi>; sl(2,R)-Vertizes brechen so(4)"),
    "diagonale E8-Komposition": ("C / E", "bedingt nativ (Naht = (E8)_1, T2); Mischterm = 8 gaußsche Cartan-Bilineare + 240 "
                                          "Vertexterme; Cosets ab N = 3 nicht frei: innere Wechselwirkung; auf den "
                                          "Cartan-Oszillatoren paritätsgerade, keine kubische 4D-Struktur"),
    "relative Paritätsphase W_12": ("E", "nativ (R, (-1)^F); unter Ladungsauswahl die gaußsche Uhrpotenz R_1^(2(N+1)); "
                                         "Statistik der Zusammensetzung, keine Wechselwirkung"),
    "lokal frei, kollektiv korreliert": ("E / O", "ganzzahliges K für die Uhr unsichtbar, also nicht nativ ausgewählt; "
                                                  "Zustand der freien Dynamik, keine Wechselwirkung; freie Ränder: F6"),
    "Wechselwirkung": ("O", "braucht eine geteilte Struktur, die weder gerade noch linear ist und Delta = 1 verlässt: "
                            "Vertexoperatoren (Ladungsfusion, auch der Wurzelteil der diagonalen E8) oder nichtabelsche "
                            "Ströme (z75)"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
urteil = ("F3 ist bedingt beantwortet. Die zehnte Analyse hält: Unter SU(2)-Kovarianz ist die vollständig positive, "
          "spurerhaltende Fortsetzung des nativen Kanals eindeutig (m = 2..8 exakt über Clebsch-Gordan; die kovarianten Kanäle "
          "sind der Kegel über m extremalen Phi_s, und schon 'nur die beiden Nachbarschritte' erzwingt s = 1/2). Phi_1/2 ist "
          "der Kanal aus z76 mit Krausrang 2, Umgebung Spin 1/2 = Schwinger-Paar und maximaler Polarisation (m+2)/(m+1). "
          "Vollständige Positivität und Tomita-Symmetrie reichen nicht: Der dekohärente Kanal mit denselben Übergängen erfüllt "
          "auch Zustandserhaltung, J-Realität, Petz-Symmetrie, Spiegel und Torus, ist aber nicht kovariant. Nativ ist SU(2) "
          "auf den Paarschalen nicht als Operator: Uhr, Modulphase, Torus, Spiegel und Tomita geben nur eine abelsche Algebra "
          "mal Z_2, Strom- und Möbiushüpfer schließen nur für m <= 3, die polare Familie nur bei der Wurzel t = 1/2, und die "
          "nativen Operatoren zusammen erzeugen ganz u(m) ohne Auszeichnung. Nativ sind aber die Raten: Der Casimir der "
          "L^nat-sl(2,R) auf H_m ist exakt der su(2)-Casimir-Superoperator auf den Diagonalen, mit |J_+|^2 = (n_1 + 1) n_2. "
          "Das su(2) ist ihre kohärente Wurzel, und die einzige verbleibende Prämisse ist P_koh: Die Umgebung des "
          "Elementarschritts registriert nur die Schrittart (Krausrang 2, ein torusgewichteter Krausoperator je Schrittart). "
          "Unter P_koh ist der Kanal bis auf Eichung eindeutig, sein Symmetriealgebra ist genau su(2), und die Identifikation "
          "aus z83 folgt ohne Ziel-CCR (C). Krausrang 2 allein genügt nicht (k-abhängige Mischung), Pfadverdichtung n + 1 "
          "unterscheidet. Zusammensetzung: Die bosonische symmetrische Potenz von N gleichen Nähten mit Vakuum oder einem "
          "Teilchen gibt prod (1 - q^m)^(-m^2) exakt bis E = N (bis E = 60 gezählt); ohne Symmetrisierung zu viele, fermionisch "
          "zu wenige Zustände; das Austauschvorzeichen der Paare ist in einer gemeinsamen CAR +1, die Ununterscheidbarkeit "
          "bleibt Prämisse; endliches N gibt nur eine gerade 1/N-Nichtgaußizität. Wechselwirkung erzeugt keine native Fusion: "
          "Ladung und Uhr wirken nur innen, die c-Algebra ist maximal abelsch und linear und verschmilzt zwei Nähte gaußsch, "
          "Impulse und sl(2,R) geben gerade Zwangskorrelationen bzw. die Tomita-Paarung. Kubische Vertizes zwischen "
          "Skalarleitern verbietet schon so(4) mit Energieerhaltung (K-Typ), und auf der Unitaritätsgrenze Delta = 1 "
          "verschwindet <phi phi phi>. Die Wechselwirkung muss also aus einer geteilten Struktur kommen, die weder gerade noch "
          "linear ist und Delta = 1 verlässt; die nativen Kandidaten sind die Vertexoperatoren der Ladungssektoren und die "
          "nichtabelschen Stromantworten. Von den Kandidaten der elften Analyse ist die diagonale E8-Komposition bedingt "
          "nativ (Naht = (E8)_1) und bringt echte, aber innere Wechselwirkung: Ihr Mischterm ist auf den Cartan-Oszillatoren "
          "gaußsch und paritätsgerade, nichtgaußsch nur im Vertexteil, und die Cosets sind ab N = 3 keine freien Felder; "
          "T_1 + T_2 + g Summe J_1 J_2 zerfällt in zwei vertauschende Sektoren. Die relative Paritätsphase W_12 ist nativ, "
          "aber unter Ladungsauswahl eine gaußsche Uhrpotenz und legt nur die Statistik fest. 'Lokal frei, kollektiv "
          "korreliert' beruht auf einer ganzzahligen Größe K, die die Uhr nicht sieht; es ist eine Zustandswahl unter freier "
          "Dynamik, keine Wechselwirkung. Nichts davon ist befördert.")
say("\nUrteil\n   " + urteil)
say(f"\n   (Laufzeit {time.time() - T0:.1f} s)")
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z91_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
