"""Z131: Ist I_step herleitbar? Der Tomita-verdoppelte Diagonalschritt gegen den nativen Einzelschritt (z83).

Leitlinie 'Rekonstruktion statt Auswahl'. H6 und die gelösten Bausteine bleiben unberührt, nichts wird befördert. Firewalls:
GRAV.NONCIRCULAR.01 (die radiale Zeit ist der parabolische Erzeuger aus den Stufenhebern, nicht der Modulfluss) und
DIMENSION.UPLIFT.FIREWALL.01 (alles ist Kinematik der Einpaarschicht). Statuskürzel: E exakt, C bedingt, L Lesart, O offen,
N numerisch. Bestand: z83 (radiale Standardform, Satz über Schalenfaktoren F G = const, nativer Möbiusheber L^nat_-1 =
Summe (k+1) psi*_(k+1) psi_k, Identifikation I_step), z76 (Zweipfadkanal, eindeutig durch Erhaltung des Originalzustands),
z91 (I_step bedingt unter P_koh 'die Umgebung registriert nur die Schrittart'; polare Familie V |L|^t schließt genau bei
t = 1/2; P_koh nicht aus CP, Tomita, Spiegel, Torus, Petz oder Uhr), z73 (su(2) braucht Oszillatorgewichte),
z92 (E^nat = [L^nat_1, L^nat_-1]/2 = L_0^osc + Q^2/2, Einzelfermion D_1/2), z117 (native Cl(2,4) auf W = T + T*,
Spinerzeuger = Fock-su(2,2) der Quarteroszillatoren), z79 (Twistorindex = Quarter), z72 (parabolischer Erzeuger),
z125 (Delta^(-i/2) = Paaruhr).

Rahmen F0 (aus z83 übernommen, nicht neu gewählt): zwei native Schrittarten S_1 (Teilchen hoch, k -> k+1) und S_2 (Loch
runter, k -> k) auf der Einpaarschicht H = (+)_m H_m mit Amplituden s_a(m, k) (n_1 = k, n_2 = m-1-k); Pfadunabhängigkeit
[S_1, S_2] = 0 (nativ für L^nat); Tomita-Verdopplung P_ab = S_a (x) J S_b J auf der gleichstufigen Standardform
D = (+)_m H_m (x) H-quer_m, K = P^+, E = (L_0 + L-quer_0)/2; Diagonaleinbettung delta|m,k> = |m,k> (x) J|m,k>. F0 enthält
z91s P_koh (ein Krausoperator je Schrittart); P_koh legt mit z76 nur die k-Form fest, nicht den Schalenfaktor.
Möbiusprämisse P_Moeb: (E, P^0, K^0 = (P^0)^+) mit P^0 = P_11 + P_22 ist eine sl(2,R) mit [K^0, P^0] = 2E, auf der
Diagonale delta (P_Moeb(delta)) bzw. auf ganz D (P_Moeb(D)).

Teil 1  I_step exakt: P^0 delta = delta L^nat_-1 (P^0 = P_11 + P_22); gleichwertig |s_a|^2 = n_a + 1, Vektorisierung
        (P^0 vec X = vec Summe_a S_a X S_a^+, auf Diagonaldichten der Kern |s_a|^2), Polarform S_a = U_a |S^nat_a|^(1/2).
        Native Daten im Fockraum (Jordan-Wigner, exakte Brüche).
Teil 2  Herleitungsversuche: (a) Tomita-Kovarianz, Modulfluss, natürlicher Kegel, KMS-Dualität (Petz), Purifikation als
        Kohärenzzustand; (b) Eindeutigkeit: Satz A (E_ext = c L_0 <=> |s_a|^2 = sqrt(c)(n_a + 1)), Cartan-Teil des
        Abschlusses, Klassifikation der erzeugten Lie-Algebren; (c) native radiale Energie E^nat: Schalenfaktorfamilie,
        P_Moeb(delta) (eine freie Größe je Schale), Satz C (P_Moeb(delta) + Spiegel C), Faktorisierung und
        Nullpunktsaufteilung, Satz D (P_Moeb(D) => I_step: Gram-Rang, Gröbner-Elimination, exakte Schalen 2..6);
        (d) Klein/Clifford (z117): Gamma^2 = Pf, Spinerzeuger, Homomorphiedefekt der Leiter;
        (e) Modulkovarianz und Parabolizität der Minkowskienergie. Negativkontrollen: reskalierte Schritte, nackte
        Hardy-Hälften, Petz-Kanal, naive Verdopplung, gekrümmte Phasen.
Teil 3  Bilanz (Äquivalenztafel über alle Familien), Prämissen, Urteil.
"""
import itertools
import json
import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
rng = np.random.default_rng(131)
BETA = 4*np.pi
TOL = 1e-9


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


def zeit(text="Laufzeit bisher"):
    print(f"   ({text} {time.time() - T0:.1f} s)", file=sys.stderr, flush=True)


def klein(x, schwelle=1e-12):
    """Residuen deterministisch anzeigen (unterhalb der Schwelle nur die Schranke)."""
    return f"< {schwelle:.0e}" if x < schwelle else f"{x:.3e}"


# =================================================================== Bausteine
SCHRITT = (lambda k: k + 1, lambda k: k)                  # S_1 Teilchen hoch, S_2 Loch runter
MMAX = 12                                                 # Familien definiert für Übergänge m -> m+1 mit m <= MMAX - 1
MP = 9                                                    # Leiterhöhe der Prädikattafel


def n_a(a, m, k):
    return k if a == 0 else m - 1 - k                     # Teilchenlage n_1 = k, Lochtiefe n_2 = m-1-k


def t_nat(a, m, k):
    return n_a(a, m, k) + 1                               # native Möbiusamplitude (z83): n_a + 1


def einseitig(M):
    lab1 = [(m, k) for m in range(1, M + 1) for k in range(m)]
    return lab1, {l: i for i, l in enumerate(lab1)}


def einzelschritte(M, s):
    lab1, ix1 = einseitig(M)
    S = []
    for a in range(2):
        X = np.zeros((len(lab1), len(lab1)), dtype=complex)
        for (m, k), i in ix1.items():
            if m < M:
                X[ix1[(m + 1, SCHRITT[a](k))], i] = s(a, m, k)
        S.append(X)
    return S


def kern(M, s):
    """Kern des Diagonalschritts auf delta: T[(m+1, k'), (m, k)] = |s_a(m, k)|^2."""
    lab1, ix1 = einseitig(M)
    T = np.zeros((len(lab1), len(lab1)))
    for (m, k), i in ix1.items():
        if m < M:
            for a in range(2):
                T[ix1[(m + 1, SCHRITT[a](k))], i] = abs(s(a, m, k))**2
    return T


def verdoppelt(M, s, bilinear=False):
    """Gleichstufige Standardform D: Labels (m, k, kb); P_ab = S_a (x) J S_b J, Amplitude s_a(m, k) conj(s_b(m, kb))."""
    lab = [(m, k, kb) for m in range(1, M + 1) for k in range(m) for kb in range(m)]
    ix = {l: i for i, l in enumerate(lab)}
    n = len(lab)
    P = {}
    for a in range(2):
        for b in range(2):
            X = np.zeros((n, n), dtype=complex)
            for (m, k, kb), i in ix.items():
                if m < M:
                    X[ix[(m + 1, SCHRITT[a](k), SCHRITT[b](kb))], i] = s(a, m, k)*np.conj(s(b, m, kb))
            P[(a, b)] = X
    E = np.diag([float(l[0]) for l in lab])
    if not bilinear:
        return lab, ix, P, E

    def bil(c_, a_, seite):                               # S_c S_a^+ auf der Original- bzw. Partnerseite
        X = np.zeros((n, n), dtype=complex)
        for (m, k, kb), i in ix.items():
            kk = k if seite == 0 else kb
            k1 = kk - (1 if a_ == 0 else 0)
            if m == 1 or not 0 <= k1 <= m - 2:
                continue
            k2 = SCHRITT[c_](k1)
            ziel = (m, k2, kb) if seite == 0 else (m, k, k2)
            amp = s(c_, m - 1, k1)*np.conj(s(a_, m - 1, k1))
            X[ix[ziel], i] = amp if seite == 0 else np.conj(amp)
        return X
    A = {(c_, a_): bil(c_, a_, 0) for c_ in range(2) for a_ in range(2)}
    B = {(c_, a_): bil(c_, a_, 1) for c_ in range(2) for a_ in range(2)}
    return lab, ix, P, E, A, B


def E_ext_von(P, K=None):
    K = {k: v.conj().T for k, v in P.items()} if K is None else K
    return sum(K[x] @ P[x] - P[x] @ K[x] for x in P)/4


def lie_dims(P, K, spalten):
    """Wie z83: Rang von {P, K}, + [K, P], + [[K, P], P oder K], + [[K, P], [K, P]] auf Spalten weg vom Rand."""
    nrm = lambda X: np.linalg.norm(X[:, spalten])
    gens = [(X, nrm(X)) for X in list(P.values()) + list(K.values())]
    s2 = [(K[x] @ P[y] - P[y] @ K[x], nrm(K[x])*nrm(P[y])) for x in K for y in P]
    s3 = [(X @ Y - Y @ X, a_*b_) for X, a_ in s2 for Y, b_ in gens]
    s4 = [(X @ Y - Y @ X, a_*b_) for (X, a_), (Y, b_) in itertools.combinations(s2, 2)]
    out, acc = [], []
    for L in (gens, s2, s3, s4):
        for X, skala in L:
            v = X[:, spalten].ravel()
            nv = np.linalg.norm(v)
            if nv > 1e-11*skala and nv > 0:
                acc.append(v/nv)
        sv = np.linalg.svd(np.array(acc), compute_uv=False)
        out.append(int(np.sum(sv > 1e-9*sv[0])))
    return out


def dims(s, M, Kfak=None):
    lab, _, P, _ = verdoppelt(M, s)
    K = {k: v.conj().T for k, v in P.items()}
    if Kfak is not None:                                  # K|m+1> = G(m) K_0|m+1> (Faktor an der Zielschale m)
        D = np.diag([Kfak(l[0]) for l in lab])
        K = {k: D @ v for k, v in K.items()}
    return lie_dims(P, K, np.array([l[0] for l in lab]) <= M - 3)


# --- Kettenlösung der Diagonal-sl(2) [K^0, P^0] = 2 E auf delta (Schale m, gegeben die Übergänge m-1 -> m)
def kette(tvor, m, u):
    tp = (lambda a, k: tvor.get((a, k), 0.0)) if tvor else (lambda a, k: 0.0)
    r = [2.0*m + tp(0, k - 1)**2 + tp(1, k)**2 for k in range(m)]
    p = [tp(0, k)*tp(1, k) for k in range(m - 1)]
    t = {(0, 0): u}
    x = r[0] - u*u
    if x <= 0:
        return None
    t[(1, 0)] = np.sqrt(x)
    for k in range(1, m):
        t2 = p[k - 1]/t[(0, k - 1)]
        x = r[k] - t2*t2
        if x <= 0:
            return None
        t[(1, k)] = t2
        t[(0, k)] = np.sqrt(x)
    return t


def gueltig(tvor, m, n=4000):
    """Intervall der Parameter u = t_1(m, 0), für die die ganze Kette positiv bleibt (Gitter + Bisektion)."""
    r0 = 2.0*m + (tvor.get((1, 0), 0.0)**2 if tvor else 0.0)
    us = np.sqrt(r0)*np.arange(1, n)/n
    idx = [i for i, u in enumerate(us) if kette(tvor, m, u) is not None]
    zusammen = idx == list(range(idx[0], idx[-1] + 1))
    lo_a, lo_b = (us[idx[0] - 1] if idx[0] > 0 else 0.0), us[idx[0]]
    hi_a, hi_b = us[idx[-1]], (us[idx[-1] + 1] if idx[-1] + 1 < len(us) else np.sqrt(r0))
    for _ in range(60):
        mid = 0.5*(lo_a + lo_b)
        lo_a, lo_b = (mid, lo_b) if kette(tvor, m, mid) is None else (lo_a, mid)
        mid = 0.5*(hi_a + hi_b)
        hi_a, hi_b = (mid, hi_b) if kette(tvor, m, mid) is not None else (hi_a, mid)
    return lo_b, hi_a, zusammen


# --- Familien von Amplituden s(a, m, k)
def aus_t(tf, zeichen=(1.0, 1.0), phase=None):
    def s(a, m, k):
        v = zeichen[a]*np.sqrt(tf(a, m, k))
        return complex(v*(phase(a, m, k) if phase is not None else 1.0))
    return s


t_ccr = lambda a, m, k: float(t_nat(a, m, k))
PHI = {(m, k): rng.uniform(0, 2*np.pi) for m in range(1, MMAX + 1) for k in range(m)}
THETA = {(a, m, k): rng.uniform(0, 2*np.pi) for a in range(2) for m in range(1, MMAX) for k in range(m)}
ZUF = {(a, m, k): rng.uniform(0.5, 2.0) for a in range(2) for m in range(1, MMAX) for k in range(m)}
SZ = {(a, m, k): rng.normal() + 1j*rng.normal() for a in range(2) for m in range(1, MMAX) for k in range(m)}
D_AFFIN = 0.2
G_AFFIN = {0: 1.0}
for _m in range(1, MMAX):                                 # g(m)^2 (m+1)^2 = 4 (m + d) + g(m-1)^2 (m-1)^2
    G_AFFIN[_m] = np.sqrt((4*(_m + D_AFFIN) + G_AFFIN[_m - 1]**2*(_m - 1)**2)/(_m + 1)**2)
KETTE, _tv, KETTE_BEREICHE = {}, None, {}
for _m in range(1, MMAX):                                 # Kettenfamilie: CCR auf m = 1, 2, u = 1.2 auf m = 3, dann Intervallmitte
    if _m < 3:
        _u = 1.0
    elif _m == 3:
        _u = 1.2
    else:
        _lo, _hi, _zs = gueltig(_tv, _m)
        KETTE_BEREICHE[_m] = (round(float(_lo), 6), round(float(_hi), 6), _zs)
        _u = 0.5*(_lo + _hi)
    _tv = kette(_tv, _m, _u)
    for (_a, _k), _v in _tv.items():
        KETTE[(_a, _m, _k)] = _v


def t_kohaerent(a, m, k):                                 # Zuflusssummen m an jeder Zielstelle, Abflusssummen nicht konstant
    if a == 0:
        return float(m) if k == m - 1 else (k + 1) + 0.25*(-1)**k
    return float(m) if k == 0 else m - ((k) + 0.25*(-1)**(k - 1))


FAM = {
    "CCR (Tomita-Wurzel)": aus_t(t_ccr),
    "flache Phasen": aus_t(t_ccr, phase=lambda a, m, k: np.exp(1j*(PHI[(m + 1, SCHRITT[a](k))] - PHI[(m, k)]))),
    "gekrümmte Phasen": aus_t(t_ccr, phase=lambda a, m, k: np.exp(1j*THETA[(a, m, k)])),
    "CCR mal 1.3 (global reskaliert)": aus_t(lambda a, m, k: 1.3*t_ccr(a, m, k)),
    "Schalenfaktor F = 1 + 1/m": aus_t(lambda a, m, k: (1 + 1/m)*t_ccr(a, m, k)),
    "Spurteil affin (d = 0.2)": aus_t(lambda a, m, k: G_AFFIN[m]*t_ccr(a, m, k)),
    "nackt (Hardy-Hälften)": aus_t(lambda a, m, k: 1.0, zeichen=(1.0, -1.0)),
    "Kanal (z76)": aus_t(lambda a, m, k: t_ccr(a, m, k)/(m + 1)),
    "Petz = Adjungierte (g = 1/m)": aus_t(lambda a, m, k: t_ccr(a, m, k)/m),
    "naiv verdoppelt (s = n_a + 1)": aus_t(lambda a, m, k: t_ccr(a, m, k)**2),
    "Nullpunkt c = 0.3": aus_t(lambda a, m, k: np.sqrt(t_ccr(a, m, k)*(t_ccr(a, m, k) + (0.3 if a == 0 else -0.3)))),
    "Nullpunkt c = -0.5": aus_t(lambda a, m, k: np.sqrt(t_ccr(a, m, k)*(t_ccr(a, m, k) + (-0.5 if a == 0 else 0.5)))),
    "Kette (P_Moeb, Abweichung Schale 3)": aus_t(lambda a, m, k: KETTE[(a, m, k)]),
    "Kohärenz ohne Zonalität": aus_t(t_kohaerent),
    "zufällig": aus_t(lambda a, m, k: ZUF[(a, m, k)]),
}
LAB1, IX1 = einseitig(MP)
MV1 = np.array([m for m, k in LAB1])
LNAT = kern(MP, FAM["CCR (Tomita-Wurzel)"])               # Matrix von L^nat_-1 auf H (Amplituden n_a + 1)


def praedikate(s, M=MP):
    """Alle Prädikate einer Familie auf der Leiter der Höhe M (Schalen weg vom Rand)."""
    T = kern(M, s)
    mv1 = np.array([m for m, k in einseitig(M)[0]])
    ln = LNAT if M == MP else kern(M, FAM["CCR (Tomita-Wurzel)"])
    d = {}
    d["r_Istep"] = float(np.abs(T - ln).max())
    maske = ln > 0
    q = T[maske]/ln[maske]
    d["lambda_fit"] = float(q.mean())
    d["r_Istep_skala"] = float(np.abs(q - q.mean()).max())
    inn1 = mv1 <= M - 1
    kom = T.T @ T - T @ T.T - 2*np.diag(mv1.astype(float))
    d["r_Moeb"] = float(np.abs(kom[np.ix_(inn1, inn1)]).max())
    labD_, _, PD_, ED_ = verdoppelt(M, s)
    P0_ = PD_[(0, 0)] + PD_[(1, 1)]
    KP_ = P0_.conj().T @ P0_ - P0_ @ P0_.conj().T                       # [K^0, P^0] auf der ganzen gleichstufigen Form D
    innD = np.array([l[0] for l in labD_]) <= M - 1
    blkKP, blkE = KP_[np.ix_(innD, innD)], ED_[np.ix_(innD, innD)]
    d["r_MoebD"] = float(np.abs(blkKP - 2*blkE).max())
    lam2 = float(np.real(np.vdot(blkE, blkKP)/np.vdot(blkE, blkE))/2)
    d["r_MoebD_skala"] = float(np.abs(blkKP - 2*lam2*blkE).max())
    d["r_C"] = float(max(abs(abs(s(0, m, k))**2 - abs(s(1, m, m - 1 - k))**2) for m in range(1, M) for k in range(m)))
    S = einzelschritte(M, s)
    d["r_comm"] = float(np.abs((S[0] @ S[1] - S[1] @ S[0])[:, mv1 <= M - 2]).max())
    aus = T.sum(axis=0)                                   # Abflusssummen x(m, k)
    ein = T.sum(axis=1)                                   # Zuflusssummen y(m, k)
    d["r_zonal"] = float(max(np.ptp(aus[mv1 == m]) for m in range(1, M)))
    d["r_petz_adj"] = float(max(abs(ein[mv1 == m] - 1).max() for m in range(2, M + 1)))
    w = np.zeros(len(mv1))
    w[0] = 1.0
    rk = 0.0
    for n in range(1, M):                                 # (P^0)^n |1>/n! gegen die gleichverteilte Schale (Purifikation)
        w = T @ w/n
        rk = max(rk, float(np.abs(w[mv1 == n + 1] - 1.0).max()))
    d["r_koh"] = rk
    lab, ix, P, E = verdoppelt(M, s)
    Ee = E_ext_von(P)
    mv = np.array([l[0] for l in lab])
    inn = mv <= M - 1
    d["r_PE"] = float(np.abs((Ee - E)[np.ix_(inn, inn)]).max())
    kb, eschale = 0.0, {}
    for m in range(1, M):
        sel = np.where(mv == m)[0]
        blk = Ee[np.ix_(sel, sel)]
        eschale[m] = round(float(blk[0, 0].real), 6)
        kb = max(kb, float(np.abs(blk - blk[0, 0]*np.eye(len(sel))).max()))
    d["r_kblind"] = kb
    d["E_ext_je_Schale"] = eschale
    rc, cs = 0.0, []
    spal = mv <= M - 2
    for x in P:
        X = (Ee @ P[x] - P[x] @ Ee)[:, spal]
        Y = P[x][:, spal]
        c = np.vdot(Y, X)/np.vdot(Y, Y)
        rc = max(rc, float(np.abs(X - c*Y).max()))
        cs.append(complex(c))
    d["r_cartan"] = rc
    d["cartan_c"] = round(float(np.mean(cs).real), 6)
    return d


# =================================================================== Teil 1
say("Teil 1: I_step exakt")
say("   I_step (z83): Auf der gleichstufigen Standardform D = (+)_m H_m (x) H-quer_m mit P^0 = P_11 + P_22, P_ab = S_a (x) J S_b J,")
say("   gilt P^0 o delta = delta o L^nat_-1 (mit K^0 o delta = delta o L^nat_1, E o delta = delta o L_0), delta|m,k> = |m,k> (x) J|m,k>.")

# --- 1.1 Native Daten im Fockraum (Jordan-Wigner, exakte Brüche)
KMIN, KMAX = -8, 7
MODEN = range(KMIN, KMAX + 1)


def bit(k):
    return 1 << (k - KMIN)


VAK = sum(bit(k) for k in range(KMIN, 0))


def nunter(s_, k):
    return bin(s_ & ((1 << (k - KMIN)) - 1)).count("1")


def bil(a, b, vec, coef=1):
    """coef psi*_a psi_b auf {Bitmaske: Koeffizient}."""
    out = {}
    for s_, cc in vec.items():
        if not (s_ >> (b - KMIN)) & 1:
            continue
        sg = -1 if nunter(s_, b) % 2 else 1
        s1 = s_ ^ bit(b)
        if (s1 >> (a - KMIN)) & 1:
            continue
        if nunter(s1, a) % 2:
            sg = -sg
        s2 = s1 | bit(a)
        out[s2] = out.get(s2, 0) + coef*sg*cc
    return {k: v for k, v in out.items() if v != 0}


def vadd(u, v, f=1):
    w = dict(u)
    for k, cc in v.items():
        w[k] = w.get(k, 0) + f*cc
    return {k: cc for k, cc in w.items() if cc != 0}


def opsum(terme):
    def f(vec):
        out = {}
        for coef, a, b in terme:
            if KMIN <= a <= KMAX and KMIN <= b <= KMAX and coef != 0:
                out = vadd(out, bil(a, b, vec, coef))
        return out
    return f


def skalar(u, v):
    return sum(cc*v.get(s_, 0) for s_, cc in u.items())


def kommutator(A_, B_, v):
    return vadd(A_(B_(v)), B_(A_(v)), -1)


LN_F = {-1: opsum([(Fr(k + 1), k + 1, k) for k in MODEN]), 1: opsum([(Fr(k), k - 1, k) for k in MODEN])}
TEIL = opsum([(Fr(k + 1), k + 1, k) for k in MODEN if k >= 0])    # Teilchenanteil S^nat_1
LOCH = opsum([(Fr(k + 1), k + 1, k) for k in MODEN if k < 0])     # Lochanteil S^nat_2


def paar(p, r):
    return bil(p, r, {VAK: 1})


ok_amp, ok_pfad, ok_sl2, ok_spiegel = True, True, True, True
for p in range(0, 5):
    for r in range(-5, 0):
        v = paar(p, r)
        w_ = LN_F[-1](v)
        hoch, runter = skalar(paar(p + 1, r), w_), skalar(paar(p, r - 1), w_)
        rest = vadd(vadd(w_, paar(p + 1, r), -hoch), paar(p, r - 1), -runter)
        ok_amp &= hoch == p + 1 and runter == -r and not rest                     # n_1 + 1, n_2 + 1
        ok_amp &= TEIL(v) == {s_: (p + 1)*c for s_, c in paar(p + 1, r).items()}
        ok_amp &= LOCH(v) == {s_: (-r)*c for s_, c in paar(p, r - 1).items()}
        ok_pfad &= TEIL(LOCH(v)) == LOCH(TEIL(v))
        ok_sl2 &= kommutator(LN_F[1], LN_F[-1], v) == {s_: 2*(p - r)*c for s_, c in v.items()}
        pb, rb = -r - 1, -p - 1                                                   # Spiegel C: (n_1, n_2) -> (n_2, n_1)
        wb = LN_F[-1](paar(pb, rb))
        ok_spiegel &= skalar(paar(pb + 1, rb), wb) == -r and skalar(paar(pb, rb - 1), wb) == p + 1
ok_d12 = True
for p in range(0, 5):
    v = {VAK | bit(p): 1}                                                         # Einzelteilchen psi*_p Omega (Q = 1)
    ok_d12 &= kommutator(LN_F[1], LN_F[-1], v) == {s_: (2*p + 1)*c for s_, c in v.items()}
for r in range(-5, 0):
    v = {VAK ^ bit(r): 1}                                                         # Einzelloch psi_r Omega (Q = -1)
    ok_d12 &= kommutator(LN_F[1], LN_F[-1], v) == {s_: (2*(-r - 1) + 1)*c for s_, c in v.items()}
ok_tief = LN_F[1](paar(0, -1)) == {} and LN_F[-1]({VAK: 1}) == {}
check("[E] Native Daten (Jordan-Wigner, 16 Moden, exakte Brüche, 25 Paarzustände): L^nat_-1 = S^nat_1 + S^nat_2 mit "
      "Teilchenschritt (n_1 + 1) und Lochschritt (n_2 + 1); die beiden nativen Schrittarten vertauschen (Pfadunabhängigkeit); "
      "[L^nat_1, L^nat_-1] = 2 L_0 = 2m auf der Einpaarschicht und = 2n + 1 auf Einzelteilchen psi*_p Omega und Einzellöchern "
      "psi_r Omega (jede Schrittart eine Einzelfermion-Leiter D_1/2, E^nat = n + 1/2 = n + Q^2/2, z92); Tiefstgewicht "
      "L^nat_1 psi*_0 psi_-1 Omega = 0 und L^nat_-1 Omega = 0; der Spiegel C (n_1, n_2) -> (n_2, n_1) vertauscht die "
      "Schrittamplituden", ok_amp and ok_pfad and ok_sl2 and ok_d12 and ok_tief and ok_spiegel)

# --- 1.2 I_step als Operatorgleichung und gleichwertige Formen
labD, ixD, PD, ED = verdoppelt(MP, FAM["CCR (Tomita-Wurzel)"])
DELTA = np.zeros((len(labD), len(LAB1)))
for (m, k), i in IX1.items():
    DELTA[ixD[(m, k, k)], i] = 1.0
P0 = PD[(0, 0)] + PD[(1, 1)]
r_p = float(np.abs(P0 @ DELTA - DELTA @ LNAT).max())
r_k = float(np.abs(P0.conj().T @ DELTA - DELTA @ LNAT.T).max())
r_e = float(np.abs(ED @ DELTA - DELTA @ np.diag(MV1.astype(float))).max())
s_z = lambda a, m, k: SZ[(a, m, k)]                                   # zufällige komplexe Amplituden
SZm = einzelschritte(MP, s_z)
labZ, ixZ, PZ, _ = verdoppelt(MP, s_z)
X = np.zeros((len(LAB1), len(LAB1)), dtype=complex)                   # blockdiagonale Matrix (+)_m X_m
for m in range(1, MP + 1):
    sel = np.where(MV1 == m)[0]
    X[np.ix_(sel, sel)] = rng.normal(size=(m, m)) + 1j*rng.normal(size=(m, m))


def vec(Y):
    return np.array([Y[IX1[(m, k)], IX1[(m, kb)]] for (m, k, kb) in labZ])


r_vec = max(float(np.abs(PZ[(a, b)] @ vec(X) - vec(SZm[a] @ X @ SZm[b].conj().T)).max()) for a in range(2) for b in range(2))
pdiag = rng.uniform(0.1, 1.0, len(LAB1))
TZ = kern(MP, s_z)
r_dicht = float(np.abs((PZ[(0, 0)] + PZ[(1, 1)]) @ vec(np.diag(pdiag)) - vec(np.diag(TZ @ pdiag))).max())
SC = einzelschritte(MP, FAM["CCR (Tomita-Wurzel)"])
U = [np.where(np.abs(Sa) > 0, 1.0, 0.0) for Sa in SC]                # nackte Verschiebungen (Hardy-Hälften, Betrag 1)
r_polar = 0.0
for a in range(2):
    mod2 = (SC[a].conj().T @ SC[a]).real                              # |S_a|^2
    Snat = U[a] @ np.diag([float(t_nat(a, m, k)) if m < MP else 0.0 for (m, k) in LAB1])
    r_polar = max(r_polar, float(np.abs(U[a] @ mod2 - Snat).max()))
TAB_I = {}
for name, s in FAM.items():
    TAB_I[name] = round(float(np.abs(kern(MP, s) - LNAT).max()), 6)
ZAHLEN["I_step_residuum_je_Familie"] = TAB_I
say(f"   |P^0 delta - delta L^nat_-1| je Familie: {TAB_I}")
erfuellt = sorted(n_ for n_, v in TAB_I.items() if v < TOL)
check(f"[E] I_step exakt formuliert (M = {MP}): für die CCR-Schritte P^0 delta = delta L^nat_-1, K^0 delta = delta L^nat_1, "
      f"E delta = delta L_0 (Residuen {klein(r_p)}, {klein(r_k)}, {klein(r_e)}); gleichwertige Formen: (i) Vektorisierung "
      f"P_ab vec X = vec(S_a X S_b^+) für zufällige komplexe Amplituden und Matrizen ({klein(r_vec)}), auf Diagonaldichten also "
      f"der klassische Kern |s_a|^2 ({klein(r_dicht)}); (ii) Polarform: CCR = U_a |S^nat_a|^(1/2), d. h. U_a |S_a|^2 = S^nat_a "
      f"({klein(r_polar)}; die CCR-Amplitude ist das geometrische Mittel der nackten Hardy-Hälfte und des nativen Schritts, "
      f"wie z83s Tomita-Wurzel und z91s polare Familie bei t = 1/2). "
      f"I_step sieht nur |s_a|^2: es gilt genau für {erfuellt}, für alle anderen Familien nicht",
      max(r_p, r_k, r_e, r_vec, r_dicht, r_polar) < 1e-12
      and erfuellt == ["CCR (Tomita-Wurzel)", "flache Phasen", "gekrümmte Phasen"])
befund("Teil 1: I_step ist die Aussage |s_a(m, k)|^2 = n_a + 1, also: der Kraus-Lift Summe_a S_a . S_a^+ des Diagonalschritts "
       "wirkt auf den Paarlagendichten genau als nativer Möbiusheber. Die Phasen der Einzelschritte sieht I_step nicht; der "
       "Rahmen F0 braucht sie flach (Pfadunabhängigkeit, nativ)", True)
zeit()

# =================================================================== Teil 2
say("\nTeil 2: Herleitungsversuche")
PRAED = {name: praedikate(s) for name, s in FAM.items()}
say(f"   (Prädikattafel für {len(PRAED)} Familien bei M = {MP})")
zeit()

# ----------------------------------------------------------- (a) Tomita, Modulfluss, Kegel, KMS
say("\n(a) Tomita-Kovarianz, Modulfluss, natürlicher Kegel, KMS-Dualität, Purifikation")
MV = 5


def standardform(M, s):
    """Volle Standardform (+)_(m, mb) H_m (x) H-quer_mb: Labels (m, mb, k, kb); P_ab hebt beide Seiten."""
    lab = [(m, mb, k, kb) for m in range(1, M + 1) for mb in range(1, M + 1) for k in range(m) for kb in range(mb)]
    ix = {l: i for i, l in enumerate(lab)}
    n = len(lab)
    P = {}
    for a in range(2):
        for b in range(2):
            X = np.zeros((n, n), dtype=complex)
            for (m, mb, k, kb), i in ix.items():
                if m < M and mb < M:
                    X[ix[(m + 1, mb + 1, SCHRITT[a](k), SCHRITT[b](kb))], i] = s(a, m, k)*np.conj(s(b, mb, kb))
            P[(a, b)] = X
    Pi = np.zeros((n, n))
    for (m, mb, k, kb), i in ix.items():
        Pi[ix[(mb, m, kb, k)], i] = 1.0
    E = np.diag([0.5*(l[0] + l[1]) for l in lab])
    logD = np.diag([-BETA*(l[0] - l[1]) for l in lab])
    return lab, P, Pi, E, logD


def J_(Pi, Xm):
    return Pi @ Xm.conj() @ Pi.T


TESTFAM = {"CCR (Tomita-Wurzel)": FAM["CCR (Tomita-Wurzel)"], "nackt (Hardy-Hälften)": FAM["nackt (Hardy-Hälften)"],
           "Kanal (z76)": FAM["Kanal (z76)"], "zufällig komplex": s_z, "gekrümmte Phasen": FAM["gekrümmte Phasen"]}
TOMITA = {}
for name, s in TESTFAM.items():
    labS, PS, Pi, ES, logD = standardform(MV, s)
    P0S = PS[(0, 0)] + PS[(1, 1)]
    K0S = P0S.conj().T
    JP0 = J_(Pi, P0S)
    TOMITA[name] = {
        "JPJ=P_ba": max(float(np.abs(J_(Pi, PS[(a, b)]) - PS[(b, a)]).max()) for a in range(2) for b in range(2)),
        "JEJ=E": float(np.abs(J_(Pi, ES) - ES).max()),
        "[E,JP0J]=+JP0J": float(np.abs(ES @ JP0 - JP0 @ ES - JP0).max()),
        "|JP0J-K0|": round(float(np.linalg.norm(JP0 - K0S)), 4),
        "[logDelta,P]": max(float(np.abs(logD @ PS[x] - PS[x] @ logD).max()) for x in PS),
    }
kegel = []
for name in ("CCR (Tomita-Wurzel)", "nackt (Hardy-Hälften)", "zufällig komplex"):
    Sk = einzelschritte(MP, TESTFAM[name])
    mins = []
    for _ in range(3):
        G = rng.normal(size=(len(LAB1), 4)) + 1j*rng.normal(size=(len(LAB1), 4))
        Xp = G @ G.conj().T                                                 # X >= 0
        mins.append(float(np.linalg.eigvalsh(sum(Sa @ Xp @ Sa.conj().T for Sa in Sk)).min()))
    kegel.append(min(mins))
ZAHLEN["tomita_modular"] = {k: {kk: (klein(vv) if kk != "|JP0J-K0|" else vv) for kk, vv in v.items()} for k, v in TOMITA.items()}
ZAHLEN["kegel_min_eigenwert"] = ["> -1e-12" if x > -1e-12 else x for x in kegel]
check("[E, negativ] Tomita und Modulfluss sind normierungsblind (volle Standardform (+)_(m, mb) H_m (x) H-quer_mb, M = 5, "
      "für CCR, nackt, Kanal, zufällige komplexe und gekrümmte Phasen): J P_ab J = P_ba, J E J = E, [log Delta, P_ab] = 0 "
      "für jede Normierung; weil J die Graduierung erhält, hat J P^0 J den Grad +1 und ist nie K^0 (|J P^0 J - K^0| >= "
      f"{min(v['|JP0J-K0|'] for v in TOMITA.values())}); P^0 = Summe_a S_a J S_a J bildet den natürlichen Kegel für jede "
      "Normierung in sich ab (kleinster Eigenwert von Summe S_a X S_a^+ für zufällige X >= 0: "
      f"{'> -1e-12' if min(kegel) > -1e-12 else min(kegel)})",
      all(v["JPJ=P_ba"] < 1e-12 and v["JEJ=E"] < 1e-12 and v["[E,JP0J]=+JP0J"] < 1e-12 and v["|JP0J-K0|"] > 0.5
          and v["[logDelta,P]"] < 1e-12 for v in TOMITA.values()) and min(kegel) > -1e-12)

# --- Petz-Dualität im KMS-Zustand (klassisch auf delta: R_(x <- y) = T_(y <- x) p_x / (T p)_y)
petz = {}
for name in ("CCR (Tomita-Wurzel)", "Petz = Adjungierte (g = 1/m)", "Kanal (z76)", "Spurteil affin (d = 0.2)"):
    T = kern(MP, FAM[name])
    res_beta, verh = [], {}
    for beta in (BETA, 1.0, 9.0):
        p = np.exp(-beta*(MV1 - 1.0))
        Tp = T @ p
        R = np.zeros_like(T)
        nz = Tp > 0
        R[:, nz] = (T.T*p[:, None])[:, nz]/Tp[nz]
        res_beta.append(R)
    beta_blind = max(float(np.abs(res_beta[0] - res_beta[i]).max()) for i in (1, 2))
    R = res_beta[0]
    for m in range(1, MP):                                                # Verhältnis Petz / Adjungierte je Schale
        sel = np.where(MV1 == m)[0]
        tgt = np.where(MV1 == m + 1)[0]
        blk_R, blk_A = R[np.ix_(sel, tgt)], T.T[np.ix_(sel, tgt)]
        q = blk_R[blk_A > 0]/blk_A[blk_A > 0]
        verh[m] = round(float(q.mean()), 6) if np.ptp(q) < 1e-12 else None
    petz[name] = {"beta_blind": klein(beta_blind), "Petz/Adjungierte je Schale": verh,
                  "dims": None}
for name in petz:
    petz[name]["dims"] = dims(FAM[name], 8)
ZAHLEN["petz"] = petz
say(f"   Petz/Adjungierte je Schale: { {k: v['Petz/Adjungierte je Schale'] for k, v in petz.items()} }")
check("[E, negativ] KMS-Dualität (Petz-Rückkanal des Diagonalkerns im KMS-Zustand, beta = 4 pi, 1, 9): der Petz-Dual ist "
      "beta-unabhängig und gleich K^0 / tau_m (tau_m = Zuflusssumme); für CCR ist Petz = K^0 / m, nicht die Adjungierte. "
      "Petz = Adjungierte verlangt tau_m = 1, in der bistochastischen Familie also g = 1/m: diese Normierung erzeugt die "
      f"Lie-Algebra {petz['Petz = Adjungierte (g = 1/m)']['dims']} (kein su(2,2)). Die KMS-Dualität wählt, wenn überhaupt, die "
      "falsche Normierung; die Temperatur trägt keine Normierungsinformation",
      all(v["beta_blind"].startswith("<") for v in petz.values())
      and petz["CCR (Tomita-Wurzel)"]["Petz/Adjungierte je Schale"] == {m: round(1/m, 6) for m in range(1, MP)}
      and petz["Petz = Adjungierte (g = 1/m)"]["Petz/Adjungierte je Schale"] == {m: 1.0 for m in range(1, MP)}
      and petz["Petz = Adjungierte (g = 1/m)"]["dims"][-1] > 15 and petz["CCR (Tomita-Wurzel)"]["dims"] == [8, 15, 15, 15])

# --- Purifikation als Kohärenzzustand
koh = {n_: (PRAED[n_]["r_koh"] < TOL, PRAED[n_]["r_zonal"] < TOL) for n_ in FAM}
koh_ja = sorted(n_ for n_, v in koh.items() if v[0])
koh_zon = sorted(n_ for n_, v in koh.items() if v[0] and v[1])
ZAHLEN["purifikation_kohaerent_(koh, zonal)"] = koh
check(f"[E] Purifikation: Der native Araki-Woods-Vektor Omega_beta = Summe_m q^(m-1) Summe_k |m,k,k> (q = e^(-beta/2), z83) "
      f"ist genau dann exp(q P^0)|1>, wenn (P^0)^n |1>/n! auf jeder Schale gleichverteilt mit Gewicht 1 ist, d. h. alle "
      f"Zuflusssummen tau_m = m (z83s Kohärenzbeobachtung als Bedingung). Das gilt für {koh_ja}: allein formabschwach (die "
      f"Familie 'Kohärenz ohne Zonalität' erfüllt es und verletzt I_step). Zusammen mit Zonalität von K^0 (Abflusssummen "
      f"konstant, K^0 erhält die KMS-Familie) bleibt genau {koh_zon}: P_KMS ist gleichwertig zu I_step (Beweis: konstante "
      f"Ab- und Zuflusssummen geben die Form g(m)(n_a + 1) (z76, Satz A3), tau_m = m g(m) = m gibt g = 1)",
      "Kohärenz ohne Zonalität" in koh_ja and "Kanal (z76)" not in koh_ja and "CCR mal 1.3 (global reskaliert)" not in koh_ja
      and koh_zon == ["CCR (Tomita-Wurzel)", "flache Phasen", "gekrümmte Phasen"])
befund("(a): Tomita-Kovarianz, Modulfluss und natürlicher Kegel nehmen jede Normierung an; J kann P nie in K überführen "
       "(Grad); die KMS-Dualität (Petz) wählt g = 1/L_0, die nicht schließt. Aus Tomita/KMS allein folgt I_step nicht. "
       "Gleichwertig zu I_step ist die Kohärenzform P_KMS (Purifikation = exp(q P^0)|1> und K^0 erhält die KMS-Familie)", True)
zeit()

# ----------------------------------------------------------- (b) Eindeutigkeit
say("\n(b) Eindeutigkeit: Satz A und die erzeugten Lie-Algebren")
# A1: E_ext diagonal mit (x x-quer - y y-quer)/4 (zufällige komplexe Amplituden)
EeZ = E_ext_von(PZ)
TZk = kern(MP, s_z)
xa, ya = TZk.sum(axis=0), TZk.sum(axis=1)
formel = np.array([(xa[IX1[(m, k)]]*xa[IX1[(m, kb)]] - ya[IX1[(m, k)]]*ya[IX1[(m, kb)]])/4 for (m, k, kb) in labZ])
mvZ = np.array([l[0] for l in labZ])
innZ = mvZ <= MP - 1
r_A1 = float(max(np.abs(EeZ - np.diag(np.diag(EeZ)))[np.ix_(innZ, innZ)].max(), np.abs(np.diag(EeZ).real - formel)[innZ].max()))
# A2: Lagrange-Identität der Lorentzform und Eindeutigkeit
x0, x1, y0, y1, cc = sp.symbols("x0 x1 y0 y1 c", positive=True)
lagrange = sp.expand((x0*x1 - y0*y1)**2 - (x0**2 - y0**2)*(x1**2 - y1**2) - (x0*y1 - x1*y0)**2) == 0
X1, Y1 = sp.symbols("X1 Y1", real=True)
los2 = sp.solve([x0*Y1 - X1*y0, X1**2 - Y1**2 - (x0**2 - y0**2)], [X1, Y1], dict=True)
pos2 = [l_ for l_ in los2 if sp.simplify(l_[X1] - x0) == 0 and sp.simplify(l_[Y1] - y0) == 0]
neg2 = [l_ for l_ in los2 if sp.simplify(l_[X1] + x0) == 0]
ok_A2 = lagrange and len(los2) == 2 and len(pos2) == 1 and len(neg2) == 1
# A3: bistochastische Form (z76): konstante Ab- und Zuflusssummen erzwingen t_a = g(m) (n_a + 1)
sig, tau = sp.symbols("sigma tau", positive=True)
ok_A3 = True
for m in range(1, 9):
    t1s = sp.symbols(f"a0:{m}")
    t2s = sp.symbols(f"b0:{m}")
    gl = [t1s[k] + t2s[k] - sig for k in range(m)]
    gl += [(t1s[j - 1] if j >= 1 else 0) + (t2s[j] if j <= m - 1 else 0) - tau for j in range(m + 1)]
    los = sp.solve(gl, list(t1s) + list(t2s) + [sig], dict=True)
    ok_A3 &= len(los) == 1 and all(sp.simplify(los[0][t1s[k]] - (k + 1)*tau/m) == 0
                                   and sp.simplify(los[0][t2s[k]] - (m - k)*tau/m) == 0 for k in range(m)) \
        and sp.simplify(los[0][sig] - (m + 1)*tau/m) == 0
# A4: E_ext = c m erzwingt g = sqrt(c): g(m)^2 (m+1)^2 = 4 c m + g(m-1)^2 (m-1)^2, exakt
ok_A4 = True
for c_ in (Fr(1), Fr(4, 9), Fr(169, 100)):
    g2 = {0: Fr(7, 3)}                                                    # beliebiger Startwert, Faktor (m-1)^2 = 0 bei m = 1
    for m in range(1, 61):
        g2[m] = (4*c_*m + g2[m - 1]*(m - 1)**2)/Fr((m + 1)**2)
    ok_A4 &= all(g2[m] == c_ for m in range(1, 61))
g2d = {0: Fr(0)}
for m in range(1, 8):                                                     # affiner Spurteil m + d: Schalenfaktor nicht konstant
    g2d[m] = (4*(m + Fr(1, 5)) + g2d[m - 1]*(m - 1)**2)/Fr((m + 1)**2)
ZAHLEN["satz_A"] = {"A1_residuum": klein(r_A1), "A4_affin_d=1/5_g2": {m: str(g2d[m]) for m in range(1, 8)}}
PE_ja = sorted(n_ for n_, v in PRAED.items() if v["r_PE"] < TOL)
kb_ja = sorted(n_ for n_, v in PRAED.items() if v["r_kblind"] < TOL)
check(f"[E] Satz A (P_E <=> I_step, beliebige komplexe Amplituden, nicht nur Schalenfunktionen): (A1) E_ext = (1/4) Summe_ab "
      f"[K_ab, P_ab] ist diagonal mit Eigenwert (x_k x_kb - y_k y_kb)/4 (x = Abfluss-, y = Zuflusssumme der |s|^2; Residuum "
      f"{klein(r_A1)}); (A2) ist E_ext auf einer Schale eine Konstante c > 0, so folgt aus der Lagrange-Identität "
      f"(x0 x1 - y0 y1)^2 - (x0^2 - y0^2)(x1^2 - y1^2) = (x0 y1 - x1 y0)^2 ({ok_A2}) x_k = x_kb, y_k = y_kb; (A3) konstante Ab- "
      f"und Zuflusssummen erzwingen die native Form t_a = g(m)(n_a + 1) (sympy, m = 1..8: {ok_A3}; das ist z76s Eindeutigkeit "
      f"des zustandserhaltenden Kanals); (A4) E_ext = c L_0 erzwingt g = sqrt(c) vom Tiefstgewicht an (exakt, m <= 60, drei "
      f"Werte c: {ok_A4}), ein affiner Spurteil m + d dagegen einen nicht konstanten Schalenfaktor. Also E_ext = L_0 <=> "
      f"P^0 delta = delta L^nat_-1. Numerisch: E_ext = L_0 genau für {PE_ja}; schalenskalar (k-blind) für {kb_ja}",
      r_A1 < 1e-10 and ok_A2 and ok_A3 and ok_A4 and len(set(g2d.values())) > 3
      and PE_ja == ["CCR (Tomita-Wurzel)", "flache Phasen", "gekrümmte Phasen"]
      and set(kb_ja) == set(PE_ja) | {"CCR mal 1.3 (global reskaliert)", "Schalenfaktor F = 1 + 1/m",
                                      "Spurteil affin (d = 0.2)", "Kanal (z76)", "Petz = Adjungierte (g = 1/m)"})

# --- Cartan-Teil des Abschlusses
cartan_ja = sorted(n_ for n_, v in PRAED.items() if v["r_cartan"] < 1e-8)
ZAHLEN["cartan"] = {n_: {"residuum": round(v["r_cartan"], 6), "c": v["cartan_c"], "E_ext_je_Schale": v["E_ext_je_Schale"]}
                    for n_, v in PRAED.items()}
check(f"[E] Cartan-Teil des Abschlusses: In jeder standardförmigen su(2,2)-Schließung ist E_ext das Zentrum von g_0 und "
      f"wirkt auf g_1 = span P als Skalar, [E_ext, P_ab] = c P_ab. Weil je zwei Ziele derselben Quelle gleiches E_ext tragen, "
      f"ist E_ext dann schalenskalar und affin (c m + d). Erfüllt von {cartan_ja} (c = "
      f"{[PRAED[n_]['cartan_c'] for n_ in cartan_ja]}); die affine Familie (d = 0.2) zeigt, dass der Cartan-Teil allein nicht "
      "genügt: erst der Quartikterm aus z83 (P K mit H(m) != H(m-1)) erzwingt d = 0, also g konstant",
      cartan_ja == sorted(["CCR (Tomita-Wurzel)", "flache Phasen", "gekrümmte Phasen", "CCR mal 1.3 (global reskaliert)",
                           "Spurteil affin (d = 0.2)"])
      and abs(PRAED["CCR mal 1.3 (global reskaliert)"]["cartan_c"] - 1.69) < 1e-6)

# --- Klassifikation der erzeugten Lie-Algebren
DIMS = {}
for name, s in FAM.items():
    DIMS[name] = {8: dims(s, 8), 10: dims(s, 10)}
DIMS["Petz-Paar (P = Kanal, K = Petz, z83)"] = {M: dims(FAM["Kanal (z76)"], M, Kfak=lambda m: (m + 1.0)/m) for M in (8, 10)}
ZAHLEN["lie_dimensionen"] = {k: {str(M): v for M, v in d.items()} for k, d in DIMS.items()}
for k, d in DIMS.items():
    say(f"   {k}: {d[8]} / {d[10]}")
schliesst = sorted(n_ for n_, d in DIMS.items() if d[8] == d[10] == [8, 15, 15, 15])
stabil = all(d[8] == d[10] for d in DIMS.values())
check(f"[N] Erzeugte Lie-Algebren (Rang von P, K und drei Kommutatorstufen, M = 8 und 10, abschneidestabil: {stabil}): "
      f"su(2,2) (8 -> 15 -> 15 -> 15) genau für {schliesst}; alle anderen wachsen, darunter die gekrümmten Phasen (I_step "
      "erfüllt, aber [P_11, P_22] != 0), die Nullpunkt- und Kettenfamilien (Diagonal-sl(2) erfüllt), der affine Spurteil "
      "(Cartan-Teil erfüllt), Kanal, Petz-Paar und die nackten Hardy-Hälften",
      stabil and schliesst == sorted(["CCR (Tomita-Wurzel)", "flache Phasen", "CCR mal 1.3 (global reskaliert)"]))
lam_out = {n_: round(PRAED[n_]["lambda_fit"], 6) for n_ in schliesst}
check(f"[E] Der native Einzelschritt ist Ausgabe des Abschlusses: für jede schließende Familie gilt P^0 delta = lambda delta "
      f"L^nat_-1 mit konstantem lambda ({lam_out}; Skalenresiduen "
      f"{[klein(PRAED[n_]['r_Istep_skala']) for n_ in schliesst]}). Bis auf die globale Konstante ist L^nat_-1 = L_-1 - "
      "(alpha - 1/2) J_-1 der einzige Diagonalschritt, dessen Tomita-Wurzel su(2,2) schließt (Satz A + Cartan-Teil + z83-Satz; "
      "numerisch bestätigt)", all(PRAED[n_]["r_Istep_skala"] < 1e-12 for n_ in schliesst)
      and lam_out == {"CCR (Tomita-Wurzel)": 1.0, "flache Phasen": 1.0, "CCR mal 1.3 (global reskaliert)": 1.3})
befund("(b): I_step ist äquivalent zu P_E (E_ext = L_0) und, bis auf die Konvention lambda, zum su(2,2)-Abschluss. Neu "
       "gegenüber z83: der Satz gilt für beliebige k-abhängige Amplituden; die Phasen müssen zusätzlich flach sein", True)
zeit()

# ----------------------------------------------------------- (c) native radiale Energie
say("\n(c) Native radiale Energie E^nat = [L^nat_1, L^nat_-1]/2 auf der Diagonale")
# c1: Schalenfunktionen P = F(E) P^osc, K = G(E) K^osc auf delta
MS = 6
Fs = sp.symbols(f"F1:{MS}")
Gs = sp.symbols(f"G1:{MS}")
lab6, ix6 = einseitig(MS)
Ts = sp.zeros(len(lab6), len(lab6))
Ks = sp.zeros(len(lab6), len(lab6))
for (m, k), i in ix6.items():
    if m < MS:
        for a in range(2):
            j = ix6[(m + 1, SCHRITT[a](k))]
            Ts[j, i] = Fs[m - 1]*t_nat(a, m, k)
            Ks[i, j] = Gs[m - 1]*t_nat(a, m, k)
komS = Ks*Ts - Ts*Ks - 2*sp.diag(*[m for (m, k) in lab6])
gls = [sp.expand(komS[i, j]) for i in range(len(lab6)) for j in range(len(lab6)) if lab6[i][0] <= MS - 1 and lab6[j][0] <= MS - 1]
gls = [g for g in gls if g != 0]
Hs = sp.symbols(f"H1:{MS}")
gls_H = [sp.expand(g.subs({Fs[i]*Gs[i]: Hs[i] for i in range(MS - 1)})) for g in gls]
losH = sp.solve(gls_H, Hs, dict=True)
ok_c1 = len(losH) == 1 and all(losH[0][h] == 1 for h in Hs)
Fz = {m: rng.uniform(0.4, 2.5) for m in range(1, MP)}
Tf = np.zeros((len(LAB1), len(LAB1)))
for (m, k), i in IX1.items():
    if m < MP:
        for a in range(2):
            Tf[IX1[(m + 1, SCHRITT[a](k))], i] = Fz[m]*t_nat(a, m, k)
Kf = np.zeros_like(Tf)
for (m, k), i in IX1.items():
    if m < MP:
        for a in range(2):
            Kf[i, IX1[(m + 1, SCHRITT[a](k))]] = t_nat(a, m, k)/Fz[m]
inn1 = MV1 <= MP - 1
r_FG1 = float(np.abs((Kf @ Tf - Tf @ Kf - 2*np.diag(MV1.astype(float)))[np.ix_(inn1, inn1)]).max())
phi_m = {m: rng.uniform(0, 2*np.pi) for m in range(1, MP + 1)}
_, _, Pcc, _ = verdoppelt(MP, FAM["CCR (Tomita-Wurzel)"])
Fd = np.diag([np.exp(1j*phi_m[l[0]]) for l in labD])                  # P|m> = F(m) P^osc|m> mit |F| = 1
theta = {1: 0.0}
for m in range(1, MP):
    theta[m + 1] = theta[m] + phi_m[m]
Ug = np.diag([np.exp(1j*theta[l[0]]) for l in labD])                   # U = e^(i theta(E)), theta(m+1) - theta(m) = phi(m)
r_eich = max(float(np.abs(Ug @ Pcc[x] @ Ug.conj().T - Pcc[x] @ Fd).max()) for x in Pcc)
check(f"[E] Schalenfunktionen (z83-Familie, sympy, M = {MS}): auf delta gilt [K^0, P^0] = 2 E^nat genau für H(m) = F(m) G(m) = 1 "
      f"für alle m (eindeutig: {ok_c1}; die k-Abhängigkeit der Diagonaleinträge erzwingt H(m) = H(m-1), das Tiefstgewicht "
      f"H(1) = 1); numerisch mit zufälligem F und G = 1/F erfüllt ({klein(r_FG1)}). Mit K = P^+ (nativ: L^nat_1 = "
      f"(L^nat_-1)^+) folgt |F| = 1, und eine Phase F(m) = e^(i phi(m)) ist reine Eichung: e^(i theta(E)) P^osc "
      f"e^(-i theta(E)) = F(E) P^osc ({klein(r_eich)}). Antwort auf (c) in dieser Familie: ja, E^nat erzwingt F G = 1 und "
      "F = 1 bis auf Eichung",
      ok_c1 and r_FG1 < 1e-10 and r_eich < 1e-12)

# c2: allgemeine Amplituden: die Diagonal-sl(2) P_Moeb lässt eine freie Größe je Schale
moeb_ja = sorted(n_ for n_, v in PRAED.items() if v["r_Moeb"] < 1e-8)
nicht_istep = sorted(n_ for n_ in moeb_ja if PRAED[n_]["r_Istep"] > TOL)
ZAHLEN["kettenfamilie_intervalle_(u_lo, u_hi, zusammenhaengend)"] = KETTE_BEREICHE
ZAHLEN["P_Moeb_residuen"] = {n_: klein(PRAED[n_]["r_Moeb"]) for n_ in FAM}
check(f"[E, negativ] Allgemeine Amplituden: P_Moeb(delta) ([K^0, P^0] = 2E auf delta) ist je Schale m eine Kette von m Kreis- und "
      f"m - 1 Hyperbelgleichungen für 2m Unbekannte, also eine freie Größe u = t_1(m, 0) je Schale. Explizite positive Lösungen "
      f"außer I_step: {nicht_istep} (die Kettenfamilie bleibt bis Schale {MMAX - 1} positiv, gültige u-Intervalle zusammenhängend "
      f"und breit, z. B. Schale 4: {KETTE_BEREICHE[4][:2]}). Sie verletzen I_step, P_E und den Abschluss; die Diagonale allein "
      "beantwortet (c) also nicht",
      nicht_istep == sorted(["Kette (P_Moeb, Abweichung Schale 3)", "Nullpunkt c = -0.5", "Nullpunkt c = 0.3"])
      and all(v[2] for v in KETTE_BEREICHE.values()) and KETTE_BEREICHE[4][1] - KETTE_BEREICHE[4][0] > 1.0
      and all(PRAED[n_]["r_PE"] > 0.01 and DIMS[n_][10][-1] > 15 for n_ in nicht_istep))

# c3: Satz C: P_Moeb + Spiegel C erzwingt I_step (monotone Kette + Induktion)
mono, f_null, mono_tab = True, True, {}
tvor = None
for m in range(1, 13):
    lo, hi, zs = gueltig(tvor, m)
    us = np.linspace(lo + 1e-9*(hi - lo), hi - 1e-9*(hi - lo), 400)
    kk = [kette(tvor, m, u) for u in us]
    t1 = np.array([[t[(0, k)] for k in range(m)] for t in kk])
    t2 = np.array([[t[(1, k)] for k in range(m)] for t in kk])
    f = t1[:, m - 1] - t2[:, 0]
    mono &= bool(np.all(np.diff(t1, axis=0) > 0) and np.all(np.diff(t2, axis=0) < 0) and np.all(np.diff(f) > 0)) and zs
    t_ccr_m = kette(tvor, m, 1.0)
    f_null &= t_ccr_m is not None and abs(t_ccr_m[(0, m - 1)] - t_ccr_m[(1, 0)]) < 1e-9 and lo < 1.0 < hi
    mono_tab[m] = (round(float(lo), 6), round(float(hi), 4))
    tvor = t_ccr_m
prev_t = None
sym_C = {}
for m in range(1, 5):                                                    # exakt: Schalen 1..4 mit Spiegel C
    xs = sp.symbols(f"x0:{m}", positive=True)
    t1d = {k: xs[k] for k in range(m)}
    t2d = {k: xs[m - 1 - k] for k in range(m)}
    tp = (lambda a, k, pv=prev_t: pv.get((a, k), sp.Integer(0))) if prev_t else (lambda a, k: sp.Integer(0))
    gl = [t1d[k]**2 + t2d[k]**2 - 2*m - tp(0, k - 1)**2 - tp(1, k)**2 for k in range(m)]
    gl += [t2d[k + 1]*t1d[k] - tp(0, k)*tp(1, k) for k in range(m - 1)]
    los = sp.solve(gl, xs, dict=True)
    sym_C[m] = [[str(l_[x]) for x in xs] for l_ in los]
    prev_t = {(0, k): sp.Integer(k + 1) for k in range(m)}
    prev_t.update({(1, k): sp.Integer(m - k) for k in range(m)})
ok_symC = all(sym_C[m] == [[str(k + 1) for k in range(m)]] for m in range(1, 5))
ZAHLEN["satz_C"] = {"monotone_kette_u_intervall_je_schale": mono_tab, "sympy_C_loesungen": sym_C}
check(f"[E] Satz C (P_Moeb(delta) + C => I_step): auf jeder Schale ist die Kette monoton (t_1(m, k) steigt, t_2(m, k) fällt mit u; "
      f"geprüft m = 1..12 auf je 400 Punkten des gültigen Intervalls: {mono}), also hat die Spiegelbedingung t_1(m, m-1) = "
      f"t_2(m, 0) höchstens eine Lösung; CCR ist eine ({f_null}). Induktion von der Schale 1 (C: t_1 = t_2 = 1) an gibt "
      f"t_a = n_a + 1 auf allen Schalen; exakt (sympy, C, Schalen 1..4) ist die einzige positive Lösung {sym_C[4][0]}. Der "
      "Spiegel C ist nativ (exakte Symmetrie des Originals, L^nat ist alpha-unabhängig, z83)",
      mono and f_null and ok_symC
      and sorted(n_ for n_ in moeb_ja if PRAED[n_]["r_C"] < TOL) == ["CCR (Tomita-Wurzel)", "flache Phasen", "gekrümmte Phasen"])

# c4: Faktorisierung über Pfadunabhängigkeit und Nullpunktsaufteilung
n1, n2 = sp.symbols("n1 n2", integer=True)
ph = sp.Function("phi")
A_ = lambda u, v: sp.exp(ph(u + 1, v) - ph(u, v))
B_ = lambda u, v: sp.exp(ph(u, v + 1) - ph(u, v))
flach = sp.simplify(B_(n1, n2)*A_(n1, n2 + 1) - A_(n1, n2)*B_(n1 + 1, n2)) == 0
kp = sp.simplify(A_(n1, n2)*B_(n1 + 1, n2 - 1)/(A_(n1, n2 - 1)*B_(n1, n2 - 1))
                 - sp.exp(2*(ph(n1 + 1, n2) - ph(n1, n2) - ph(n1 + 1, n2 - 1) + ph(n1, n2 - 1)))) == 0
cS, nS = sp.symbols("c n")
a2 = {-1: sp.Integer(0)}
for j in range(0, 8):                                                    # a(n)^2 - a(n-1)^2 = 2n + 1 + c
    a2[j] = sp.expand(a2[j - 1] + 2*j + 1 + cS)
ok_null = all(sp.expand(a2[j] - (j + 1)*(j + 1 + cS)) == 0 for j in range(8))
phi_zuf = {(u, v): 0.3*rng.normal() for u in range(MMAX + 1) for v in range(MMAX + 1)}
s_flach = lambda a, m, k: complex(np.exp(phi_zuf[(k + 1, m - 1 - k)] - phi_zuf[(k, m - 1 - k)]) if a == 0 else
                                  np.exp(phi_zuf[(k, m - k)] - phi_zuf[(k, m - 1 - k)]))
pr_flach = praedikate(s_flach)
comm_moeb = sorted(n_ for n_ in moeb_ja if PRAED[n_]["r_comm"] < 1e-8)
ZAHLEN["faktorisierung"] = {"flache_nicht_faktorisierte_Familie": {"r_comm": klein(pr_flach["r_comm"]),
                                                                    "r_Moeb": round(pr_flach["r_Moeb"], 4)},
                            "P_Moeb_und_P_comm": comm_moeb}
check(f"[E] Faktorisierung: Pfadunabhängigkeit [S_1, S_2] = 0 heißt log t = Gradient eines Potentials phi ({flach}); der "
      f"Nebendiagonalteil von P_Moeb(delta) ist dann exp(2 Delta_1 Delta_2 phi) = 1 ({kp}), also phi = f(n_1) + g(n_2) und t_1 = a(n_1), "
      f"t_2 = b(n_2). Der Diagonalteil trennt: a(n)^2 = (n+1)(n+1+c), b(n)^2 = (n+1)(n+1-c) ({ok_null}), d. h. "
      f"D_h1 (x) D_h2 mit h_1 + h_2 = 1: die einzige Restfreiheit ist die Aufteilung der Nullpunktsenergie auf die Schrittarten. "
      f"Kontrolle: eine flache, nicht faktorisierte Familie verletzt P_Moeb(delta) (Residuum {pr_flach['r_Moeb']:.3f}); P_Moeb(delta) und "
      f"Pfadunabhängigkeit zusammen erfüllen genau {comm_moeb}. Nativ legt z92 h_a = 1/2 fest (jede Schrittart ist eine "
      "Einzelfermion D_1/2, E^nat = n + Q^2/2, Teil 1), ebenso der Spiegel C: c = 0, I_step",
      flach and kp and ok_null and pr_flach["r_comm"] < 1e-12 and pr_flach["r_Moeb"] > 0.01
      and comm_moeb == ["CCR (Tomita-Wurzel)", "Nullpunkt c = -0.5", "Nullpunkt c = 0.3", "flache Phasen"])
# c5: Satz D: P_Moeb auf der ganzen gleichstufigen Form D erzwingt I_step ohne Spiegel und ohne Pfadunabhängigkeit
P0z = PZ[(0, 0)] + PZ[(1, 1)]
KPz = P0z.conj().T @ P0z - P0z @ P0z.conj().T
tz = lambda a, m, k: abs(s_z(a, m, k))**2 if (1 <= m < MMAX and 0 <= k < m) else 0.0
r_D1 = 0.0
for i, (m, k, kb) in enumerate(labZ):
    if m > MP - 1:
        continue
    gram = sum(tz(a, m, k)*tz(a, m, kb) for a in range(2)) - tz(0, m - 1, k - 1)*tz(0, m - 1, kb - 1) \
        - tz(1, m - 1, k)*tz(1, m - 1, kb)
    r_D1 = max(r_D1, abs(KPz[i, i] - gram))
    for j, (m2, k2, kb2) in enumerate(labZ):
        if j != i and abs(KPz[j, i]) > 1e-12 and not (m2 == m and k2 - k == kb2 - kb and abs(k2 - k) == 1):
            r_D1 = max(r_D1, abs(KPz[j, i]))
xs4 = sp.symbols("x0:4")
V3 = sp.Matrix([[0, xs4[0]], [xs4[1], xs4[2]], [xs4[3], 0]])
sS = sp.symbols("s")
rang_lemma = sp.expand((V3*V3.T + sS*sp.ones(3, 3)).det() - sS*(xs4[0]*xs4[1] - xs4[0]*xs4[3] + xs4[2]*xs4[3])**2) == 0
a0_, a1_, b0_, b1_, c_s = sp.symbols("a0 a1 b0 b1 c")
GB = sp.groebner([a0_**2 + b0_**2 - (5 - c_s), a1_**2 + b1_**2 - (5 + c_s), a0_*a1_ + b0_*b1_ - 4,
                  a0_**2*b1_**2 - (1 - c_s**2), a1_*b0_ - a1_*b1_ - a0_*b0_], a0_, a1_, b0_, b1_, c_s, order="lex")
elim_c = [sp.factor(g) for g in GB.exprs if g.free_symbols <= {c_s}]
ok_elim = len(elim_c) == 1 and sp.expand(elim_c[0] - c_s**3*(c_s - 1)*(c_s + 1)) == 0
sym_D = {}
for m in range(2, 7):                                                    # gegeben CCR auf Schale m-1: alle positiven Lösungen
    p_ = sp.symbols(f"p0:{m}", positive=True)
    q_ = sp.symbols(f"q0:{m}", positive=True)
    gl = [p_[k]*p_[kb] + q_[k]*q_[kb] - (k if k >= 1 else 0)*(kb if kb >= 1 else 0) - (m - 1 - k)*(m - 1 - kb) - 2*m
          for k in range(m) for kb in range(k, m)]
    gl += [q_[k + 1]*p_[k] - (k + 1)*(m - 1 - k) for k in range(m - 1)]
    los = sp.solve(gl, list(p_) + list(q_), dict=True)
    sym_D[m] = [[str(l_[x]) for x in list(p_) + list(q_)] for l_ in los]
ok_symD = all(sym_D[m] == [[str(k + 1) for k in range(m)] + [str(m - k) for k in range(m)]] for m in range(2, 7))
moebD_ja = sorted(n_ for n_, v in PRAED.items() if v["r_MoebD"] < 1e-8)
moebD_skala = sorted(n_ for n_, v in PRAED.items() if v["r_MoebD_skala"] < 1e-8)
ZAHLEN["satz_D"] = {"D1_residuum": klein(r_D1), "eliminationsideal_c": [str(g) for g in elim_c],
                    "sympy_loesungen_schale_2_bis_6": sym_D, "P_Moeb_D_erfuellt": moebD_ja,
                    "P_Moeb_D_bis_Skala": moebD_skala,
                    "P_Moeb_D_residuen": {n_: round(v["r_MoebD"], 6) for n_, v in PRAED.items()}}
check(f"[E] Satz D (P_Moeb auf der ganzen gleichstufigen Form D => I_step, ohne Spiegel C und ohne Pfadunabhängigkeit): "
      f"(D1) auf D ist [K^0, P^0] diagonal gleich der Gram-Form tau(k).tau(kb) - tau'(k).tau'(kb) (tau = Ab-, tau' = "
      f"Zuflussvektor (t_1, t_2)) plus Verschiebungen (k +- 1, kb +- 1), die genau die Plaketten der Diagonale sind (zufällige "
      f"komplexe Amplituden, {klein(r_D1)}); P_Moeb auf D heißt also: die Gram-Matrix der Abflussvektoren ist 2m 11^T + "
      f"Gram der Zuflussvektoren, vom Rang <= 2. (D2) Ranglemma: det(V V^T + s 11^T) = s (1 . (V_1 x V_2))^2 ({rang_lemma}), "
      f"also verlangt Schale m+1 >= 3, dass 1 im Spaltenraum der Zuflussvektoren liegt. (D3) Die Nullpunktsaufteilung c der "
      f"Schale 1 (t_1 = sqrt(1+c), t_2 = sqrt(1-c)) mit Schale 2 und der Rangbedingung von Schale 3: Eliminationsideal "
      f"{elim_c} (Gröbner), also c = 0 (c = +-1 hieße ein verschwindender Schritt). (D4) Gegeben CCR darunter, ist auf den "
      f"Schalen 2..6 die einzige positive Lösung CCR (sympy, {ok_symD}). Tafel: P_Moeb auf D gilt genau für {moebD_ja}, bis "
      f"auf die Skala für {moebD_skala} (= die su(2,2)-Familien); die Nullpunkt- und Kettenfamilien, die P_Moeb auf der "
      "Diagonale erfüllen, scheitern auf D, ebenso die gekrümmten Phasen (Phasenausrichtung der Plaketten)",
      r_D1 < 1e-10 and rang_lemma and ok_elim and ok_symD
      and moebD_ja == ["CCR (Tomita-Wurzel)", "flache Phasen"]
      and moebD_skala == sorted(["CCR (Tomita-Wurzel)", "CCR mal 1.3 (global reskaliert)", "flache Phasen"])
      and all(PRAED[n_]["r_MoebD"] > 0.01 for n_ in ("Nullpunkt c = 0.3", "Nullpunkt c = -0.5",
                                                      "Kette (P_Moeb, Abweichung Schale 3)", "gekrümmte Phasen")))
befund("(c): E^nat erzwingt in der Schalenfunktionsfamilie F G = 1 und F = 1 bis auf Eichung. Allgemein reicht die Diagonale "
       "nicht (eine freie Größe je Schale); mit Pfadunabhängigkeit bleibt genau die Nullpunktsaufteilung h_1 + h_2 = 1, die "
       "die native Einzelfermion-Energie n + 1/2 (z92, Q^2/2 = 1/2) bzw. der Spiegel C festlegt. Verlangt man die "
       "Möbiusstruktur der radialen Zeit auf der ganzen gleichstufigen Standardform statt nur auf der Diagonale, erzwingt sie "
       "I_step allein (Satz D)", True)
zeit()

# ----------------------------------------------------------- (d) Klein/Clifford
say("\n(d) Klein/Clifford-Struktur aus z117")
LC = np.zeros((4,)*4)
for _p in itertools.permutations(range(4)):
    LC[_p] = round(np.linalg.det(np.eye(4)[list(_p)]))


def stern(Q):
    return 0.5*np.einsum("abcd,cd->ab", LC, Q)


def pf(Q):
    return Q[0, 1]*Q[2, 3] - Q[0, 2]*Q[1, 3] + Q[0, 3]*Q[1, 2]


Z44 = np.zeros((4, 4))


def Gam(Q):
    return np.block([[Z44, -stern(Q)], [Q, Z44]])


ok_quad = True
for _ in range(6):
    Q = rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4))
    Q = Q - Q.T
    ok_quad &= np.allclose(Gam(Q) @ Gam(Q), pf(Q)*np.eye(8), atol=1e-12)
BIV = []
for i, j in itertools.combinations(range(4), 2):
    Q = np.zeros((4, 4), complex)
    Q[i, j], Q[j, i] = 1, -1
    BIV.append(Q)
spin, blockform = [], True
for Q1, Q2 in itertools.combinations(BIV, 2):
    C = (Gam(Q1) @ Gam(Q2) - Gam(Q2) @ Gam(Q1))/4
    blockform &= np.allclose(C[:4, 4:], 0) and np.allclose(C[4:, :4], 0) and np.allclose(C[:4, :4], -C[4:, 4:].T)
    spin.append(C[4:, 4:])
Vsp = np.array([x.ravel() for x in spin]).T
rang_spin = int(np.linalg.matrix_rank(Vsp))
ETA = np.diag([1.0, 1.0, -1.0, -1.0])


def Eij(i, j):
    Z = np.zeros((4, 4))
    Z[i, j] = 1
    return Z


phiA = lambda c_, a_: Eij(c_, a_) @ ETA
phiB = lambda d_, b_: Eij(2 + b_, 2 + d_) @ ETA
phiP = lambda a_, b_: Eij(a_, 2 + b_) @ ETA
phiK = lambda a_, b_: Eij(2 + b_, a_) @ ETA
Y15 = [ETA/2, phiA(0, 1), phiA(1, 0), (phiA(0, 0) - phiA(1, 1))/2, phiB(0, 1), phiB(1, 0), (phiB(0, 0) - phiB(1, 1))/2,
       phiP(0, 0), phiP(0, 1), phiP(1, 0), phiP(1, 1), phiK(0, 0), phiK(0, 1), phiK(1, 0), phiK(1, 1)]
r_span = max(float(np.linalg.norm(Vsp @ np.linalg.lstsq(Vsp, y.ravel().astype(complex), rcond=None)[0] - y.ravel()))
             for y in Y15)
VY = np.array([y.ravel() for y in Y15]).T
g15 = np.zeros((15, 15, 15), dtype=complex)
for i in range(15):
    for j in range(15):
        g15[i, j] = np.linalg.lstsq(VY, (Y15[i] @ Y15[j] - Y15[j] @ Y15[i]).ravel(), rcond=None)[0]
jac = 0.0
for i, j, k in itertools.combinations(range(15), 3):
    jac = max(jac, float(np.abs(g15[i, j] @ g15[:, k, :] + g15[j, k] @ g15[:, i, :] + g15[k, i] @ g15[:, j, :]).max()))
adj_ok = all(np.allclose(ETA @ phiP(a, b).conj().T @ ETA, phiK(a, b)) for a in range(2) for b in range(2))
check(f"[E] Clifford (z117): Gamma(Q) = [[0, -*Q], [Q, 0]] erfüllt Gamma^2 = Pf(Q) ({ok_quad}); die Spinerzeuger "
      f"[Gamma, Gamma]/4 sind blockdiagonal diag(-M^T, M) ({blockform}) und spannen sl(4) (Rang {rang_spin}); die 15 "
      f"Twistorbilder (E = eta/2, A0, B0, P, K) liegen darin ({klein(r_span)}), K = P^(G4-adjungiert) ({adj_ok}), Jacobi "
      f"{klein(jac)}: auf dem Einteilchenraum W sind die Strukturkonstanten und damit die Translationsnormierung relativ zu E "
      "festgelegt", ok_quad and blockform and rang_spin == 15 and r_span < 1e-12 and adj_ok and jac < 1e-12)

HOMO = {}
for name in ("CCR (Tomita-Wurzel)", "CCR mal 1.3 (global reskaliert)", "Kanal (z76)", "Nullpunkt c = 0.3",
             "Schalenfaktor F = 1 + 1/m"):
    M = 8
    lab, ix, P, E, A, B = verdoppelt(M, FAM[name], bilinear=True)
    K = {k: v.conj().T for k, v in P.items()}
    X15 = [E, A[(0, 1)], A[(1, 0)], (A[(0, 0)] - A[(1, 1)])/2, B[(0, 1)], B[(1, 0)], (B[(0, 0)] - B[(1, 1)])/2,
           P[(0, 0)], P[(0, 1)], P[(1, 0)], P[(1, 1)], K[(0, 0)], K[(0, 1)], K[(1, 0)], K[(1, 1)]]
    mv = np.array([l[0] for l in lab])
    spal = mv <= M - 2
    dfk = 0.0
    for i in range(15):
        for j in range(15):
            R = X15[i] @ X15[j] - X15[j] @ X15[i] - sum(g15[i, j, k]*X15[k] for k in range(15) if abs(g15[i, j, k]) > 1e-14)
            dfk = max(dfk, float(np.abs(R[:, spal]).max()))
    HOMO[name] = round(dfk, 6)
ZAHLEN["clifford_homomorphiedefekt"] = HOMO
say(f"   Homomorphiedefekt Twistorbild -> Leiter (E -> L_0): {HOMO}")
check(f"[E, bedingt] Übertragung auf die Leiter: die Abbildung der 15 Twistorbilder auf die Leiteroperatoren (E -> L_0, P -> P_ab, "
      f"K -> P_ab^+, A, B -> Schrittbilineare) ist genau für CCR ein Lie-Homomorphismus (Defekt {klein(HOMO['CCR (Tomita-Wurzel)'], 1e-10)}); "
      f"für die global reskalierte CCR ist der Defekt {HOMO['CCR mal 1.3 (global reskaliert)']} (die Clifford-"
      "Strukturkonstanten legen mit E -> L_0 und K = P^+ lambda = 1 fest), für Kanal, Nullpunkt und Schalenfaktor gibt es "
      "keinen Homomorphismus. Die Gammas leben auf W; "
      "die Form der Leiteramplituden sehen sie nur über die Fock-Identifikation 'Leiter = Fockraum von W mit CCR-Form' "
      "(P_Fock), und die ist I_step in anderer Gestalt. Zudem ist der Twistorindex von z117 der Quarterindex (z79), nicht "
      "die Schrittart: (d) ist keine Herleitung", HOMO["CCR (Tomita-Wurzel)"] < 1e-10
      and HOMO["CCR mal 1.3 (global reskaliert)"] > 0.5 and min(HOMO["Kanal (z76)"], HOMO["Nullpunkt c = 0.3"],
                                                              HOMO["Schalenfaktor F = 1 + 1/m"]) > 0.1)
zeit()

# ----------------------------------------------------------- (e) Minkowskienergie
say("\n(e) Modulkovarianz und Parabolizität der Minkowskienergie H_M = E + (P^0 + K^0)/2")
modkov = {}
for name in ("CCR (Tomita-Wurzel)", "nackt (Hardy-Hälften)", "Kanal (z76)", "zufällig komplex"):
    labS, PS, Pi, ES, logD = standardform(MV, TESTFAM[name])
    P0S = PS[(0, 0)] + PS[(1, 1)]
    HM = ES + 0.5*(P0S + P0S.conj().T)
    kom = float(np.abs(logD @ HM - HM @ logD).max())
    t_ = 0.3
    Ut = np.diag(np.exp(1j*t_*np.diag(logD)))
    borchers = float(np.linalg.norm(Ut @ HM @ Ut.conj().T - np.exp(-2*np.pi*t_)*HM))
    modkov[name] = (klein(kom), round(borchers, 4))
ZAHLEN["modulkovarianz_([logDelta, H_M], Borchers-Defekt)"] = modkov
check("[E, negativ] Modulkovarianz ist leer: auf der vollen Standardform vertauscht log Delta = -beta (L_0 - L-quer_0) "
      "(die Helizität, z83) mit H_M für jede Normierung (CCR, nackt, Kanal, zufällig); eine Borchers/Hislop-Longo-Beziehung "
      "Delta^(it) H_M Delta^(-it) = e^(-2 pi t) H_M ist daher für jede Normierung verletzt (Defekte "
      f"{[v[1] for v in modkov.values()]}). Der native Modulfluss kann die radiale Zeit weder dilatieren noch normieren "
      "(GRAV.NONCIRCULAR.01)", all(v[0].startswith("<") and v[1] > 0.1 for v in modkov.values()))

spek = {}
for lam in (0.8, 1.0, 1.25):
    spek[lam] = {}
    for M in (15, 30, 60):
        lab1, _ = einseitig(M)
        Ln = kern(M, FAM["CCR (Tomita-Wurzel)"])
        H = np.diag([float(m) for m, k in lab1]) + lam*(Ln + Ln.T)/2
        spek[lam][M] = round(float(np.linalg.eigvalsh(H)[0]), 6)
spek_null = {}
for M in (15, 30, 60):
    lab1, _ = einseitig(M)
    Ln = kern(M, FAM["Nullpunkt c = 0.3"])
    H = np.diag([float(m) for m, k in lab1]) + (Ln + Ln.T)/2
    spek_null[M] = round(float(np.linalg.eigvalsh(H)[0]), 6)
ZAHLEN["minkowski_kleinster_eigenwert_auf_delta"] = {str(k): v for k, v in spek.items()}
ZAHLEN["minkowski_nullpunkt_c=0.3"] = spek_null
say(f"   kleinster Eigenwert von E + lambda (L^nat_-1 + L^nat_1)/2 auf delta: {spek}; Nullpunktfamilie: {spek_null}")
check(f"[E] Parabolizität legt nur die Skala fest: mit P^0 = lambda L^nat_-1 auf delta ist H_M für lambda = 0.8 elliptisch "
      f"(Lücke {spek[0.8][60]} = sqrt(1 - lambda^2) h, h = 1), für lambda = 1 parabolisch (kleinster Eigenwert "
      f"{spek[1.0][15]}, {spek[1.0][30]}, {spek[1.0][60]} ~ 1,8/M -> 0, Spektrum [0, inf) wie z83) und für lambda = 1.25 "
      f"hyperbolisch (nach unten unbeschränkt: {spek[1.25][15]}, {spek[1.25][30]}, {spek[1.25][60]}). Die Forderung "
      "'Minkowskizeit = parabolisch und positiv' fixiert lambda = 1, also genau den Normierungsteil von I_step",
      all(abs(spek[0.8][M] - 0.6) < 1e-9 for M in (15, 30, 60)) and spek[1.0][15] > spek[1.0][30] > spek[1.0][60] > 0
      and 1.6 < 60*spek[1.0][60] < 1.84 and spek[1.25][15] > spek[1.25][30] > spek[1.25][60] and spek[1.25][15] < 0)
check(f"[E, negativ] Formblindheit auf delta: Die Nullpunktfamilie c = 0.3 (verletzt I_step) gibt auf delta exakt dasselbe "
      f"Minkowskispektrum wie CCR ({spek_null} gegen {spek[1.0]}), weil D_h1 (x) D_h2 mit h_1 + h_2 = 1 als sl(2)-Modul "
      "dasselbe ist wie D_1/2 (x) D_1/2 = (+)_(h >= 1) D_h. Positivität und Parabolizität der Zeit sehen die Form nicht",
      all(abs(spek_null[M] - spek[1.0][M]) < 1e-9 for M in (15, 30, 60)))
befund("(e): Der Modulfluss ist die Helizität und vertauscht mit der radialen Zeit; er trägt keine Normierung. Die "
       "Parabolizität der Minkowskienergie legt die globale Skala lambda fest, nicht die Form", True)

# ----------------------------------------------------------- Negativkontrollen
say("\nNegativkontrollen")
KONTR = ["CCR mal 1.3 (global reskaliert)", "Schalenfaktor F = 1 + 1/m", "nackt (Hardy-Hälften)", "Kanal (z76)",
         "naiv verdoppelt (s = n_a + 1)", "gekrümmte Phasen"]
ktab = {n_: {"I_step": PRAED[n_]["r_Istep"] < TOL, "P_E": PRAED[n_]["r_PE"] < TOL, "su(2,2)": DIMS[n_][10] == [8, 15, 15, 15],
             "E_ext_je_Schale": PRAED[n_]["E_ext_je_Schale"]} for n_ in KONTR}
ktab["Petz-Paar (P = Kanal, K = Petz, z83)"] = {"I_step": False, "P_E": False,
                                                "su(2,2)": DIMS["Petz-Paar (P = Kanal, K = Petz, z83)"][10] == [8, 15, 15, 15]}
ZAHLEN["negativkontrollen"] = ktab
check("[E, negativ] Kontrollen: global reskalierte Schritte schließen, verletzen aber I_step um genau lambda und P_E um "
      f"lambda^2 (E_ext = {PRAED['CCR mal 1.3 (global reskaliert)']['E_ext_je_Schale'][3]} auf Schale 3); Schalenfaktor "
      "F = 1 + 1/m, nackte Hardy-Hälften (E_ext nicht schalenskalar, Toeplitz), Kanal (E_ext = (2m-1)/(4m^2)), Petz-Paar und "
      "naive Verdopplung des nativen Schritts verletzen I_step, P_E und den Abschluss; gekrümmte Phasen erfüllen I_step und "
      "P_E, schließen aber nicht ([P, P] != 0)",
      ktab["CCR mal 1.3 (global reskaliert)"]["su(2,2)"] and not ktab["CCR mal 1.3 (global reskaliert)"]["I_step"]
      and abs(PRAED["CCR mal 1.3 (global reskaliert)"]["E_ext_je_Schale"][3] - 1.69*3) < 1e-6
      and all(not ktab[n_]["I_step"] and not ktab[n_]["P_E"] and not ktab[n_]["su(2,2)"]
              for n_ in ("Schalenfaktor F = 1 + 1/m", "nackt (Hardy-Hälften)", "Kanal (z76)", "naiv verdoppelt (s = n_a + 1)"))
      and not ktab["Petz-Paar (P = Kanal, K = Petz, z83)"]["su(2,2)"]
      and ktab["gekrümmte Phasen"]["I_step"] and ktab["gekrümmte Phasen"]["P_E"] and not ktab["gekrümmte Phasen"]["su(2,2)"]
      and abs(PRAED["Kanal (z76)"]["E_ext_je_Schale"][4] - 7/64) < 1e-6)
zeit()

# =================================================================== Teil 3
say("\nTeil 3: Bilanz")
TAFEL = {}
for n_, v in PRAED.items():
    TAFEL[n_] = {
        "I_step": v["r_Istep"] < TOL,
        "P_E": v["r_PE"] < TOL,
        "P_Moeb(delta)": v["r_Moeb"] < 1e-8,
        "P_Moeb(D)": v["r_MoebD"] < 1e-8,
        "C": v["r_C"] < TOL,
        "flach": v["r_comm"] < 1e-8,
        "P_KMS": v["r_koh"] < TOL and v["r_zonal"] < TOL,
        "su(2,2)": DIMS[n_][10] == [8, 15, 15, 15],
        "I_step bis lambda": v["r_Istep_skala"] < TOL,
    }
ZAHLEN["aequivalenztafel"] = TAFEL
for n_, v in TAFEL.items():
    say(f"   {n_}: " + ", ".join(f"{k}={'ja' if b else 'nein'}" for k, b in v.items()))
aeq = all(v["I_step"] == v["P_E"] == (v["P_Moeb(delta)"] and v["C"]) == v["P_KMS"] for v in TAFEL.values())
aeqD = all(v["P_Moeb(D)"] == (v["I_step"] and v["flach"]) for v in TAFEL.values())
absch = all(v["su(2,2)"] == (v["I_step bis lambda"] and v["flach"]) for v in TAFEL.values())
moeb_flach = sorted(n_ for n_, v in TAFEL.items() if v["P_Moeb(delta)"] and v["flach"])
check(f"[C] Äquivalenztafel über alle {len(TAFEL)} Familien: I_step <=> P_E <=> (P_Moeb(delta) und C) <=> P_KMS ({aeq}); "
      f"P_Moeb(D) <=> (I_step und flache Phasen) ({aeqD}); su(2,2)-Abschluss <=> (I_step bis auf lambda und flache Phasen) "
      f"({absch}); P_Moeb(delta) mit Pfadunabhängigkeit lässt genau die Nullpunktsaufteilung offen ({moeb_flach}). Tomita, "
      "Modulfluss und Kegel nehmen jede Familie an (Teil 2a)",
      aeq and aeqD and absch
      and moeb_flach == ["CCR (Tomita-Wurzel)", "Nullpunkt c = -0.5", "Nullpunkt c = 0.3", "flache Phasen"])

STATUS = {
    "I_step als Operatorgleichung": ("E", "P^0 delta = delta L^nat_-1 <=> |s_a|^2 = n_a + 1 <=> Kraus-Lift auf Paarlagendichten "
                                          "= nativer Möbiusheber; Polarform CCR = U |S^nat|^(1/2) (wie z83, z91)"),
    "Tomita, Modulfluss, Kegel": ("E (negativ)", "normierungsblind; J erhält den Grad, also nie J P J = K"),
    "KMS-Dualität (Petz)": ("E (negativ)", "beta-blind, Petz = K^0/tau; Petz = Adjungierte wählt g = 1/m, schließt nicht"),
    "Satz A": ("E", "E_ext = c L_0 <=> |s_a|^2 = sqrt(c)(n_a + 1), beliebige Amplituden (verallgemeinert z83, nutzt z76)"),
    "Abschluss": ("E/N", "su(2,2) <=> I_step bis lambda und flache Phasen; L^nat_-1 ist Ausgabe"),
    "E^nat (c)": ("E", "Schalenfunktionen: F G = 1, F = 1 bis Eichung; allgemein: Diagonale allein eine Größe je Schale frei"),
    "Satz C": ("E", "P_Moeb(delta) + Spiegel C => I_step (monotone Kette, Induktion)"),
    "Nullpunktsaufteilung": ("E", "P_Moeb(delta) + Pfadunabhängigkeit => D_h1 (x) D_h2, h_1 + h_2 = 1; nativ h_a = 1/2 "
                                  "(z92, C)"),
    "Satz D": ("E", "P_Moeb(D) => I_step ohne C: Gram-Rang <= 2, Eliminationsideal c^3 (c^2 - 1), exakt bis Schale 6"),
    "Verhältnis zu z91 (P_koh)": ("E", "P_koh ist Teil von F0 und legt mit z76 nur die k-Form fest; den Schalenfaktor legt erst "
                                       "P_Moeb (bzw. P_E, P_KMS) fest"),
    "Clifford (d)": ("C", "fixiert lambda relativ zu E; Form nur über P_Fock (= I_step), Twistorindex = Quarter"),
    "Minkowskienergie (e)": ("E", "Modulfluss vertauscht (leer); Parabolizität fixiert lambda, nicht die Form"),
    "P_Moeb selbst": ("O", "nicht aus Tomita/KMS herleitbar; minimale kinematische Form der radialen Zeit"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
PRAEMISSEN = {
    "P_Moeb (neu benannt)": "Die radiale Zeitleiter (E, P^0, K^0 = (P^0)^+), P^0 = P_11 + P_22, ist mit der nativen "
                            "Graduierung eine sl(2,R): [K^0, P^0] = 2E. Auf ganz D erzwingt sie I_step allein (Satz D); auf "
                            "der Diagonale ist sie echt schwächer als I_step und braucht den nativen Spiegel C (Satz C)",
    "C (nativ)": "Spiegel x -> -x mit Schrittartentausch, exakte Symmetrie des Originals (z83), mit derselben "
                 "Normierungsfunktion in beiden Quartersektoren (L^nat ist alpha-unabhängig); bei Pfadunabhängigkeit "
                 "gleichwertig die native Einzelfermion-Energie n + 1/2 je Schrittart (z92)",
    "Rahmen F0 (nativ, aus z83)": "zwei Schrittarten mit nativem Muster, Pfadunabhängigkeit, Tomita-Verdopplung mit einem "
                                  "Einzelschritt je Schrittart (= z91s P_koh), K = P^+",
    "gleichwertige Fassungen": "P_E (E_ext = L_0 auf D), P_KMS (Purifikation = exp(q P^0)|1>, K^0 erhält die KMS-Familie)",
    "Konvention": "lambda = 1 ist die Normierung [K^0, P^0] = 2E bzw. 'H_M parabolisch'; physikalisch leer",
}
ZAHLEN["praemissen"] = PRAEMISSEN

urteil = ("I_step ist nicht nativ im engen Sinn und nicht irreduzibel: Es ist bedingt auf die schwächere Möbiusprämisse P_Moeb "
          "zurückgeführt. (1) Exakt: I_step heißt P^0 delta = delta L^nat_-1, gleichwertig |s_a(m,k)|^2 = n_a + 1: der "
          "Kraus-Lift Summe_a S_a . S_a^+ des Tomita-verdoppelten Diagonalschritts wirkt auf den Paarlagendichten als nativer "
          "Möbiusheber; die CCR-Schritte sind U_a |S^nat_a|^(1/2) (Tomita-Wurzel wie z83, polare Wurzel t = 1/2 wie z91). Die "
          "Phasen der Einzelschritte sieht I_step nicht; der Abschluss braucht sie flach (gekrümmte Phasen erfüllen I_step, "
          "schließen aber nicht). (2) Tomita-Kovarianz, Modulfluss und natürlicher Kegel nehmen jede Normierung an; J erhält "
          "die Graduierung, also ist J P^0 J = K^0 unmöglich; die KMS-Dualität (Petz) ist beta-blind und wählt mit Petz = "
          "Adjungierte die Normierung g = 1/L_0, die nicht schließt; der Modulfluss ist die Helizität und vertauscht mit der "
          "Minkowskienergie. Aus Tomita/KMS allein folgt I_step nicht (vgl. z91: P_koh nicht aus CP, Tomita, Petz, Uhr). "
          "(3) Satz A (neu, für beliebige k-abhängige Amplituden, nicht nur Schalenfunktionen wie z83): E_ext = (1/4) Summe "
          "[K_ab, P_ab] = c L_0 genau dann, wenn |s_a|^2 = sqrt(c)(n_a + 1); Beweis über die Lagrange-Identität der "
          "Lorentzform (E_ext schalenskalar <=> bistochastisch), z76s Eindeutigkeit des zustandserhaltenden Kanals und die "
          "Tiefstgewichtsrekursion. Mit dem Cartan-Teil des Abschlusses und z83s Quartikterm schließt su(2,2) genau für "
          "I_step bis auf lambda (und flache Phasen); der native Möbiusheber L^nat_-1 ist Ausgabe des Abschlusses, nicht "
          "Eingabe. (4) Zu E^nat: In der Schalenfunktionsfamilie erzwingt [K^0, P^0] = 2 E^nat auf der Diagonale F G = 1 und "
          "mit K = P^+ F = 1 bis auf Eichung. Für allgemeine Amplituden lässt P_Moeb(delta) eine freie Größe je Schale "
          "(explizite positive Kettenlösungen); mit Pfadunabhängigkeit bleibt genau die Nullpunktsaufteilung h_1 + h_2 = 1 "
          "(D_h1 (x) D_h2), die die native Einzelfermion-Energie n + 1/2 = n + Q^2/2 (z92) und der native Spiegel C beide auf "
          "h_a = 1/2 festlegen; Satz C: P_Moeb(delta) + C => I_step (monotone Kette, Induktion). Satz D: Verlangt man die "
          "Möbiusstruktur der radialen Zeit auf der ganzen gleichstufigen Standardform D, so ist die Gram-Matrix der "
          "Abflussvektoren 2m 11^T + Gram der Zuflussvektoren vom Rang <= 2; die Nullpunktsaufteilung c stirbt an der "
          "Rangbedingung der dritten Schale (Eliminationsideal c^3 (c^2 - 1)), und gegeben CCR darunter ist CCR auf den "
          "Schalen 2..6 die einzige positive Lösung: P_Moeb(D) => I_step ohne C und ohne Pfadunabhängigkeit. z91s P_koh ist "
          "im Rahmen enthalten und legt nur die k-Form fest; den Schalenfaktor legt erst P_Moeb fest. (5) Clifford (z117) "
          "fixiert die Skala lambda relativ zu E, sieht die Leiterform aber nur über die Fock-Identifikation, die I_step "
          "selbst ist; der Twistorindex dort ist der Quarterindex. Die Parabolizität der Minkowskienergie fixiert ebenfalls "
          "nur lambda; auf der Diagonale ist sie formblind. Urteil: I_step = P_Moeb(D), gleichwertig P_Moeb(delta) + C "
          "(nativ), mit der Konvention lambda = 1. Die schwächere Prämisse heißt P_Moeb: Die radiale Zeitleiter (E, P^0, "
          "K^0) ist mit der nativen Graduierung eine sl(2,R). Gleichwertige Fassungen sind P_E (E_ext = L_0) und P_KMS "
          "(Purifikation = Kohärenzzustand von P^0, K^0 erhält die KMS-Familie). P_Moeb selbst ist aus Tomita/KMS nicht "
          "herleitbar (O). Folge für die Achsen: Spinorindex = Schrittart, Moden = Energiestufen und Zeit = radial konform "
          "hängen nicht mehr an einer Schrittidentifikation, sondern an P_Moeb, also an der minimalen kinematischen Form der "
          "Zeitachse selbst ('die radiale Zeit ist ein Möbiuserzeuger mit der nativen Energie'); die Normierung, L^nat_-1 "
          "als Diagonalschritt, die 4D-Konformgruppe, die Schrittart als Spinorindex und die Energiestufen folgen. Die "
          "Prämissenzahl des besten Kandidaten sinkt nicht von selbst (I_step -> P_Moeb, eins zu eins), die Prämisse wird aber "
          "schwächer und ist Teil der Lesart 'Zeit = radial konform'. Nichts befördert, Ledger unverändert.")
ZAHLEN["laufzeit_s"] = "< 180"
say("\nUrteil\n   " + urteil)
zeit("Laufzeit")
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z131_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z131.txt").write_text("\n".join(OUT) + "\n", encoding="utf-8")
