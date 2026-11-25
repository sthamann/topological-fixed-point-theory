"""Z107: Ursprung schließen -- P_tn aus <T> = 0 und Starrheit S(n) (Folgeauftrag N8).

Auftrag N8 (H6h46, Nächste Schritte): (1) die Zeitneutralität des Gesetzes P_tn aus der Bewegungsgleichung der
hintergrundfreien Naht herleiten (H6h4 A*: <T> = 0, für E8 gleich E6(tau) = 0); dann bliebe eine eigene Ursprungsprämisse.
(2) Satz S(n) für verschiedene m: S(3) für allgemeine ungerade m_1 < m_2 < m_3 exakt oder über eine vollständige endliche
Reduktion; allgemeines Argument (Induktion mit Nichtentartung und Kern 2n - 3) oder Hindernis. (3) Bilanz: eine oder zwei
Prämissen, Geltungsbereich, K19.
Leitlinie 'Rekonstruktion statt Auswahl'; gelöste Bausteine werden nicht geöffnet, nichts wird befördert. Status E/C/L/O.
Bestand: z94 (Gesetzestorus tau_Gesetz = i/T_a = i beta/(4 pi), P_clock <=> P_tn, S(2) allgemein, Kern 2n - 3, Fenster),
z87 (Saat rho_1, Quartermoden x = m/4 = NS-Moden der Überlagerung; die Rahmenkriterien K, H lassen jedes beta zu),
z3/z4 (A*: d log chi_E8/d tau = -(2 pi i/3) E6/E4, tau = i einzige stationäre Form, Sattel = Legendre-Struktur), n1.

Teil 1  P_tn aus <T> = 0. Unter dem Rahmen (Modulfluss = Drehung der Hardy-Naht) ist der Gesetzestorus aus z94 (Raumzyklus =
        Antiperiode des Modulflusses, Zeitzyklus = KMS) der euklidische Torus der Naht mit Spinstruktur (A, A); seine
        Zustandssumme ist das NS-Zeichen theta_3/eta eines komplexen Fermions (Jacobi). A*(a) (große Diffeomorphismen, hier S)
        legt die Casimirenergie -1/24 fest, A*(b) <T> = 0 heißt d log(theta_3/eta)/d tau = (pi i/12)(theta_2^4 - theta_4^4)
        = 0, der (A, A)-Faktor von E6 = (1/2)(theta_4^4 - theta_2^4)(theta_2^4 + theta_3^4)(theta_3^4 + theta_4^4).
        Nullstellen: auf der Achse genau y = 1, in der Halbebene genau Gamma_theta i. Nichtzirkularität: anderer Inhalt
        verschiebt oder tilgt den stationären Punkt; A* lässt P_sym frei und enthält die Minimalität; Basis- gegen
        Überlagerungstorus.
Teil 2  Starrheit S(n). Normalform (Teilchen-Loch: mu_j = eps_j m_j = 4 k_j + 1; Nichtentartung = verschiedene k_j),
        Kernzahl 2n - 3 aus den Restklassen; Partialbruchlemma (Zusatzkern = Polteil, oberster Block unberührt; n = 3: jede
        Paarprojektion invertierbar, also S(3) <=> ein Paar spaltet ab; n = 4: zwei Paare); ganzzahlige Starrheit aller 56
        Tripel m_3 <= 15 im Fenster; reelle Relaxation für n = 2 starr, für n = 3 nicht: S(n) ist arithmetisch (Hindernis).
Teil 3  Bilanz: Prämissenzählung, Status, K19, Urteil.
"""
import itertools
import json
import time
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import flint
import mpmath as mp
import numpy as np
import sympy as sp
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.optimize import root as sroot
from scipy.sparse import csr_matrix, lil_matrix, vstack

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
mp.mp.dps = 40
rng = np.random.default_rng(107)
T_START = time.time()
yS, sS = sp.symbols("y s")
MSER = 80                                                               # Ordnung der exakten q-Reihen (in p = q^(1/2))


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


# =================================================================== Teil 1
say("Teil 1: P_tn aus <T> = 0 (A*)")


# --- exakte ganzzahlige Potenzreihen in p = q^(1/2)
def smul(a, b):
    c = [0]*MSER
    for i, x in enumerate(a):
        if x:
            for j in range(MSER - i):
                if b[j]:
                    c[i + j] += x*b[j]
    return c


def spow(a, k):
    r = [1] + [0]*(MSER - 1)
    for _ in range(k):
        r = smul(r, a)
    return r


TH3, TH4, T2 = [0]*MSER, [0]*MSER, [0]*MSER
for n_ in range(-20, 21):
    if n_*n_ < MSER:
        TH3[n_*n_] += 1
        TH4[n_*n_] += (-1)**abs(n_)
for n_ in range(20):
    if n_*(n_ + 1) < MSER:
        T2[n_*(n_ + 1)] += 1
TH2_4 = [0] + [16*c for c in spow(T2, 4)][:MSER - 1]                   # theta_2^4 = 16 p (sum p^(n(n+1)))^4
TH3_4, TH4_4 = spow(TH3, 4), spow(TH4, 4)
GETA = [1] + [0]*(MSER - 1)                                             # prod (1 - p^(2n)) = q^(-1/24) eta
for n_ in range(1, MSER):
    if 2*n_ < MSER:
        GETA = [GETA[i] - (GETA[i - 2*n_] if i >= 2*n_ else 0) for i in range(MSER)]
PNS = [1] + [0]*(MSER - 1)                                              # prod (1 + p^(2n-1))^2 = Tr e^(-K_beta)
for n_ in range(1, MSER):
    if 2*n_ - 1 < MSER:
        for _ in range(2):
            PNS = [PNS[i] + (PNS[i - (2*n_ - 1)] if i >= 2*n_ - 1 else 0) for i in range(MSER)]


# --- mpmath-Funktionen
def th(k, tau):
    return mp.jtheta(k, 0, mp.exp(1j*mp.pi*tau))


def eta(tau):
    return mp.exp(2j*mp.pi*tau/24)*mp.qp(mp.exp(2j*mp.pi*tau))


def eis(k, tau, N=150):
    q = mp.exp(2j*mp.pi*tau)
    c = {4: 240, 6: -504}[k]
    return 1 + c*mp.fsum(n**(k - 1)*q**n/(1 - q**n) for n in range(1, N))


def logZ_ns(y):
    """log Z_NS(iy) = log theta_3(iy)/eta(iy), reell auf der Achse."""
    return mp.log(mp.re(th(3, 1j*y))) - (-mp.pi*y/12 + mp.log(mp.qp(mp.exp(-2*mp.pi*y))))


def E_ns(y):
    """Energie des (A, A)-Gesetzestorus mit Casimir: E(y) = -d_y log Z(iy) (Einheiten: d tau = i dy)."""
    return -mp.diff(logZ_ns, y)


# --- 1.1 Formulierung: Gesetzestorus = NS-Torus (A, A)
xs_o = np.array(sorted([m/4 for m in range(1, 16, 2)] + [-m/4 for m in range(1, 16, 2)]))


def antiperiode(beta):
    """wie z94: kleinste T > 0 mit e^(i beta x T) = -1 auf allen Quartermoden."""
    Ta = 4*np.pi/abs(beta)
    alle = np.abs(np.exp(1j*beta*xs_o*Ta) + 1).max() < 1e-9
    ts = np.linspace(1e-3, Ta*(1 - 1e-3), 4000)
    kleiner = bool((np.abs(np.exp(1j*np.outer(ts, beta*xs_o)) + 1).max(axis=1) < 1e-6).any())
    return Ta if alle and not kleiner else None


tau_tab = {bf: 1.0/antiperiode(bf*np.pi) for bf in (2, 4, 8, 20)}
jtp = smul(PNS, GETA) == TH3
spur = mp.fprod((1 + mp.exp(-mp.pi*m))**2 for m in range(1, 400, 2))   # Tr e^(-K) bei beta = 4 pi, x = m/4
q_i = mp.exp(-2*mp.pi)
ns_i = q_i**(mp.mpf(1)/24)*mp.re(th(3, 1j))/mp.re(eta(1j))
ZAHLEN["gesetzestorus_Im_tau_bei_beta/pi"] = tau_tab
ZAHLEN["Tr_e^-K_bei_beta_4pi"] = {"Produkt": mp.nstr(spur, 25), "q^(1/24) theta_3/eta (i)": mp.nstr(ns_i, 25)}
check("Formulierung: der Gesetzestorus aus z94 (kleinste Antiperiode T_a des Modulflusses K_beta = beta sum x n_x auf den "
      f"Quartermoden, KMS-Antiperiode 1) hat tau_Gesetz = i/T_a = i beta/(4 pi) (Im tau = {list(tau_tab.values())} bei beta = "
      "2, 4, 8, 20 pi); mit q = e^(2 pi i tau_Gesetz) = e^(-beta/2) ist das Boltzmanngewicht der Mode x = m/4 genau q^(m/2), "
      "also Tr e^(-K_beta) = prod_(m ungerade) (1 + q^(m/2))^2 = q^(1/24) theta_3/eta (Jacobi-Tripelprodukt exakt bis p^80, "
      f"p = q^(1/2); bei beta = 4 pi auf 40 Stellen: Abweichung {mp.nstr(abs(spur - ns_i), 2)}). Raumzyklus = Antiperiode "
      "(Drehung um 4 pi, NS), Zeitzyklus = KMS (antiperiodisch): Spinstruktur (A, A), Inhalt ein komplexes NS-Fermion (aus der "
      "Viertelzählung des Rahmens, ungerade m; weder Uhr noch P_tn gehen ein, beta bleibt frei)",
      all(abs(tau_tab[b] - b/4) < 1e-12 for b in tau_tab) and jtp and abs(spur - ns_i) < mp.mpf(10)**-35)

# --- 1.2 Casimirenergie aus A*(a): S-Invarianz der (A, A)-Zustandssumme


def logP(y):
    """log prod (1 + q^(n - 1/2))^2 auf der Achse, q = e^(-2 pi y): Gesetz ohne Casimir (normierter Operator)."""
    return 2*mp.log(mp.qp(-mp.exp(-mp.pi*y), mp.exp(-2*mp.pi*y)))


E0_aus_S = {float(y): (logP(1/mp.mpf(y)) - logP(mp.mpf(y)))/(2*mp.pi*(1/mp.mpf(y) - mp.mpf(y))) for y in (0.5, 0.8, 1.25, 2, 3)}
ohne_cas = [-mp.diff(logP, mp.mpf(y)) for y in np.linspace(0.2, 6, 30)]
ZAHLEN["Casimir_aus_S_invarianz"] = {k: mp.nstr(v, 30) for k, v in E0_aus_S.items()}
check("Casimir aus A*(a): der große Diffeomorphismus S (Zyklentausch) erhält (A, A); die Zustandssumme q^E0 prod (1 + q^(n-1/2))^2 "
      "ist genau dann S-invariant, wenn E0 = -1/24: aus Z(iy) = Z(i/y) folgt E0 = [log P(i/y) - log P(iy)]/(2 pi (1/y - y)) = "
      f"-1/24 für jedes y (Abweichung <= {mp.nstr(max(abs(v + mp.mpf(1)/24) for v in E0_aus_S.values()), 2)} bei y = 0,5 .. 3). "
      "Ohne Casimir (das normierte Gesetz rho = e^(-K)/Z) ist -d_y log P > 0 überall (min "
      f"{mp.nstr(min(ohne_cas), 4)}): das Gesetz allein sieht A* nicht, der Casimirterm kommt aus A*(a), nicht aus einer Wahl",
      max(abs(v + mp.mpf(1)/24) for v in E0_aus_S.values()) < mp.mpf(10)**-35 and min(ohne_cas) > 0)

# --- 1.3 Spannungstensor des Gesetzes und E6
pf_ = [i*c for i, c in enumerate(TH3)]
pg_ = [i*c for i, c in enumerate(GETA)]
links = [12*(a - b) - c for a, b, c in zip(smul(pf_, GETA), smul(pg_, TH3), smul(TH3, GETA))]
rechts = smul([a - b for a, b in zip(TH2_4, TH4_4)], smul(TH3, GETA))
E6S = [0]*MSER
E6S[0] = 1
for n_ in range(1, MSER):
    if 2*n_ < MSER:
        E6S[2*n_] = -504*sum(d**5 for d in range(1, n_ + 1) if n_ % d == 0)
fakt = smul(smul([a - b for a, b in zip(TH4_4, TH2_4)], [a + b for a, b in zip(TH2_4, TH3_4)]),
            [a + b for a, b in zip(TH3_4, TH4_4)])
t0_ = mp.mpc(0.13, 1.21)
dlog = {k: mp.diff(lambda t: mp.log(th(k, t)/eta(t)), t0_) for k in (3, 4, 2)}
a2, a3, a4 = th(2, t0_)**4, th(3, t0_)**4, th(4, t0_)**4
soll = {3: (mp.pi*1j/12)*(a2 - a4), 4: -(mp.pi*1j/12)*(a2 + a3), 2: (mp.pi*1j/12)*(a3 + a4)}
spin_ok = max(abs(dlog[k] - soll[k]) for k in dlog)
prod_ok = abs(dlog[3]*dlog[4]*dlog[2] - 2*(mp.pi*1j/12)**3*eis(6, t0_))
ZAHLEN["spannung_spinstrukturen_abweichung"] = {"log-Ableitungen": mp.nstr(spin_ok, 3), "Produkt gegen E6": mp.nstr(prod_ok, 3)}
check("Spannungstensor des Gesetzes (exakt bis p^80): q d/dq log(theta_3/eta) = (theta_2^4 - theta_4^4)/24, also <T> ~ d log "
      "Z/d tau = (pi i/12)(theta_2^4 - theta_4^4); und E6 = (1/2)(theta_4^4 - theta_2^4)(theta_2^4 + theta_3^4)(theta_3^4 + "
      "theta_4^4) (exakt bis p^80). Die drei Faktoren sind die Spannungen der drei geraden Spinstrukturen (A, A), (A, P), "
      f"(P, A) (40 Stellen: {mp.nstr(spin_ok, 2)}), ihr Produkt ist 2 (pi i/12)^3 E6 ({mp.nstr(prod_ok, 2)}): die Gesetzes-"
      "Entsprechung von E6(tau) = 0 ist der (A, A)-Faktor theta_4^4 - theta_2^4",
      links == rechts and [2*e for e in E6S] == fakt and spin_ok < mp.mpf(10)**-30 and prod_ok < mp.mpf(10)**-30)

# --- 1.4 Nullstellen
ys_ax = (0.5, 0.8, 0.9, 1.0, 1.1, 1.25, 2.0, 5.0)
E_ax = {y: E_ns(mp.mpf(y)) for y in ys_ax}
fallend = all(E_ax[a_] > E_ax[b_] for a_, b_ in zip(ys_ax, ys_ax[1:]))
konvex = min(mp.diff(logZ_ns, mp.mpf(y), 2) for y in np.linspace(0.25, 6, 24))
wurzel = mp.findroot(E_ns, mp.mpf("0.9"))
s_sym = abs(E_ns(mp.mpf(2)) + E_ns(mp.mpf("0.5"))/4)
pos_f = min(min(mp.re(th(2, 1j*y)**4 + th(3, 1j*y)**4), mp.re(th(3, 1j*y)**4 + th(4, 1j*y)**4)) for y in np.linspace(0.2, 6, 30))
e6_achse = {float(y): float(mp.re(eis(6, 1j*mp.mpf(y)))) for y in (0.5, 0.8, 1.25, 2.0)}
ZAHLEN["achse_E(y)"] = {str(k): mp.nstr(v, 12) for k, v in E_ax.items()}
nn_ = np.arange(-30, 31)


def th_np(k, t):
    t = np.asarray(t, complex)[..., None]
    if k == 3:
        return np.sum(np.exp(1j*np.pi*nn_**2*t), axis=-1)
    if k == 4:
        return np.sum((-1.0)**np.abs(nn_)*np.exp(1j*np.pi*nn_**2*t), axis=-1)
    return np.sum(np.exp(1j*np.pi*(nn_ + 0.5)**2*t), axis=-1)


f_aa = lambda t: th_np(2, t)**4 - th_np(4, t)**4
in_theta = lambda z: z.imag > 0.45 and abs(z.real) <= 1 + 1e-9 and abs(z) >= 1 - 1e-9
ng = 201
Xg, Yg = np.meshgrid(np.linspace(-1.05, 1.05, ng), np.linspace(0.45, 3.0, ng))
Ag = np.abs(f_aa(Xg + 1j*Yg))
nullst = []
for i in range(1, ng - 1):
    for j in range(1, ng - 1):
        if Ag[i, j] <= Ag[i - 1:i + 2, j - 1:j + 2].min():
            sol = sroot(lambda v: [f_aa(v[0] + 1j*v[1]).real, f_aa(v[0] + 1j*v[1]).imag], [Xg[i, j], Yg[i, j]], tol=1e-14)
            z = sol.x[0] + 1j*sol.x[1]
            if abs(f_aa(z)) < 1e-9 and in_theta(z) and not any(abs(z - w_) < 1e-6 for w_ in nullst):
                nullst.append(z)
# Spitze 1: f(1 - 1/tau') = tau'^2 (theta_3^4 + theta_4^4)(tau'); der abgeschnittene Spitzenbereich liegt bei Im tau' >= 2.1
spitze_trafo = abs(f_aa(1 - 1/(2.5j)) - (2.5j)**2*(th_np(3, 2.5j)**4 + th_np(4, 2.5j)**4))
tp = np.linspace(0, 1, 41)[:, None] + 1j*np.linspace(2.1, 12, 41)[None, :]
spitze_min = float(np.abs(th_np(3, tp)**4 + th_np(4, tp)**4).min())
lam_i = mp.re(th(2, 1j)**4/th(3, 1j)**4)
bahn, bahn_achse = 0, []
for a_, b_, c_, d_ in itertools.product(range(-12, 13), repeat=4):
    if a_*d_ - b_*c_ == 1:
        bahn += 1
        if a_*c_ + b_*d_ == 0:
            bahn_achse.append(c_*c_ + d_*d_)
ZAHLEN["nullstellen"] = {"Fundamentalbereich Thetagruppe": [str(np.round(z, 12)) for z in nullst], "lambda(i)": mp.nstr(lam_i, 30),
                         "Spitze 1 min |theta_3^4 + theta_4^4|": spitze_min, "SL2Z-Matrizen geprüft": bahn,
                         "davon i auf der Achse (c^2 + d^2)": dict(Counter(bahn_achse)), "E6(iy)": e6_achse}
check(f"Nullstellen der Gesetzesspannung: auf der Achse tau = iy (ein Zustand hat reelles beta, der Torus ist rechteckig) ist "
      f"E(y) = -d_y log Z streng fallend ({', '.join(f'{y}: {float(E_ax[y]):+.4f}' for y in ys_ax)}), d_y^2 log Z = (2 pi)^2 "
      f"Var > 0 (min {mp.nstr(konvex, 4)}), E(y) = -E(1/y)/y^2 (S, {mp.nstr(s_sym, 2)}); einzige Nullstelle y = "
      f"{mp.nstr(wurzel, 30)}. In der Halbebene ist theta_2^4 - theta_4^4 = theta_3^4 (2 lambda - 1) mit dem Hauptmodul lambda "
      f"von Gamma(2) (lambda(i) = {mp.nstr(lam_i, 20)}) und theta_3 ohne Nullstelle: Nullstellen genau Gamma(2) i = Gamma_theta i "
      f"(Gitterzensus im Fundamentalbereich der Thetagruppe: {[str(np.round(z, 9)) for z in nullst]}; Spitze 1 frei, "
      f"|theta_3^4 + theta_4^4| >= {spitze_min:.3f}); die SL(2,Z)-Bahn von i trifft die Achse nur in i ({bahn} Matrizen, alle "
      f"Achsentreffer mit c^2 + d^2 = 1). Die zwei übrigen E6-Faktoren sind auf der Achse positiv (min {mp.nstr(pos_f, 4)}): "
      "auch E6(tau_Gesetz) = 0 wählt genau tau_Gesetz = i (der E8-Inhalt von H6h4 gibt denselben Punkt)",
      fallend and konvex > 0 and abs(wurzel - 1) < mp.mpf(10)**-30 and s_sym < mp.mpf(10)**-30 and len(nullst) == 1
      and abs(nullst[0] - 1j) < 1e-10 and abs(lam_i - mp.mpf(1)/2) < mp.mpf(10)**-35 and spitze_trafo < 1e-9
      and spitze_min > 1.5 and set(bahn_achse) == {1} and pos_f > 0 and e6_achse[0.5] < 0 < e6_achse[1.25])
befund("P_tn folgt aus A*: Unter dem Rahmen ist der Gesetzestorus der Nahttorus mit (A, A); A*(a) macht seinen Inhalt "
       "S-invariant (Casimir -1/24), A*(b) verlangt <T> = 0, und die einzige Lösung auf den Zustandstori ist tau_Gesetz = i, "
       "also beta = 4 pi: genau die Selbstdualität des Gesetzestorus (z94). Mechanismus: S-Symmetrie der Funktion legt die "
       "Nullstelle auf den Fixpunkt, strenge Konvexität (Var > 0) macht sie eindeutig",
       all(CHECKS[k] for k in CHECKS if k.startswith(("Formulierung", "Casimir", "Spannungstensor", "Nullstellen"))))

# --- 1.5 Nichtzirkularität


def logZ_zwist(y, a):
    s = mp.nsum(lambda n: mp.exp(-mp.pi*y*(n + a)**2), [-mp.inf, mp.inf])
    return mp.log(s) - (-mp.pi*y/12 + mp.log(mp.qp(mp.exp(-2*mp.pi*y))))


def stationaer(fun, start=1.0):
    """Nullstelle von fun auf der Achse (0,05 < y < 50) per findroot, sonst None."""
    try:
        r = mp.findroot(fun, mp.mpf(start))
    except Exception:
        return None
    if not (mp.im(r) == 0 and 0.05 < mp.re(r) < 50):
        return None
    return r if abs(fun(r)) < mp.mpf(10)**-25 else None


zwist = {}
for a_ in (mp.mpf(1)/16, mp.mpf(1)/8, mp.mpf(3)/8, mp.mpf(1)/2):
    Ea = lambda y, a_=a_: -mp.diff(lambda t: logZ_zwist(t, a_), y)
    grid = [Ea(mp.mpf(y)) for y in np.linspace(0.3, 6, 20)]
    zwist[float(a_)] = {"y*": stationaer(Ea), "E(1)": Ea(mp.mpf(1)), "Vorzeichenwechsel": min(grid) < 0 < max(grid),
                        "S-Abweichung": abs(logZ_zwist(mp.mpf(2), a_) - logZ_zwist(mp.mpf("0.5"), a_))}
ap_vz = [mp.re(mp.diff(lambda t: mp.log(mp.re(th(4, 1j*t))) + mp.pi*t/12 - mp.log(mp.qp(mp.exp(-2*mp.pi*t))), mp.mpf(y)))
         for y in np.linspace(0.3, 6, 20)]
cas = {}
for E0 in (mp.mpf(0), -mp.mpf(1)/48, -mp.mpf(1)/24, -mp.mpf(1)/12):
    Ec = lambda y, E0=E0: 2*mp.pi*E0 - mp.diff(logP, y)
    cas[float(E0)] = stationaer(Ec)
ZAHLEN["nichtzirkular_zwist"] = {k: {kk: (mp.nstr(vv, 10) if isinstance(vv, mp.mpf) else vv) for kk, vv in v.items()}
                                 for k, v in zwist.items()}
ZAHLEN["nichtzirkular_casimir_y*"] = {k: (None if v is None else mp.nstr(v, 10)) for k, v in cas.items()}
z16, z8 = zwist[1/16], zwist[1/8]
check("Nichtzirkularität (Inhalt): 'spannungsfrei' und 'zeitneutral' sind verschiedene Prädikate. Verdrillter Inhalt "
      f"theta[a, 0]/eta (nicht S-invariant, Abweichung {mp.nstr(z8['S-Abweichung'], 3)} bei a = 1/8) ist bei y* = "
      f"{mp.nstr(z16['y*'], 8)} (a = 1/16) bzw. {mp.nstr(z8['y*'], 8)} (a = 1/8) spannungsfrei, also nicht am selbstdualen "
      f"Torus, und dort (y = 1) nicht spannungsfrei (E(1) = {mp.nstr(z16['E(1)'], 4)}, {mp.nstr(z8['E(1)'], 4)}); a = 3/8, "
      "Ramond (P, A) und (A, P) haben gar keinen stationären Punkt auf der Achse. Andere Casimirenergie (nicht S-invariante "
      f"Vervollständigung): E0 = 0 keiner, E0 = -1/48: y* = {mp.nstr(cas[-1/48], 8)}, E0 = -1/12: y* = {mp.nstr(cas[-1/12], 8)}; "
      f"nur E0 = -1/24 gibt y* = {mp.nstr(cas[-1/24], 8)}. A*(b) allein enthält tau = i also nicht; erst A*(a) (S-invarianter "
      "Inhalt) legt den Punkt auf den Fixpunkt",
      z16["y*"] is not None and z8["y*"] is not None and abs(z16["y*"] - 1) > 0.01 and abs(z8["y*"] - 1) > 0.05
      and abs(z8["E(1)"]) > 0.01 and z8["S-Abweichung"] > 0.01 and zwist[3/8]["y*"] is None
      and not zwist[3/8]["Vorzeichenwechsel"] and zwist[0.5]["y*"] is None and not zwist[0.5]["Vorzeichenwechsel"]
      and (min(ap_vz) > 0 or max(ap_vz) < 0) and cas[0.0] is None and abs(cas[-1/24] - 1) < mp.mpf(10)**-25
      and abs(cas[-1/48] - 1) > 0.1 and abs(cas[-1/12] - 1) > 0.1)

# Unabhängigkeit von P_sym und Minimalität (Qubit wie z94)
sx, sy_ = np.array([[0, 1], [1, 0]], dtype=complex), np.array([[0, -1j], [1j, 0]])


def j_marken(th_):
    z = [1, -1, np.exp(1j*th_), -np.exp(1j*th_)]
    lam = (z[0] - z[2])*(z[1] - z[3])/((z[0] - z[3])*(z[1] - z[2]))
    return float((256*(1 - lam + lam**2)**3/(lam**2*(1 - lam)**2)).real)


A3, B3 = sx, np.cos(np.pi/3)*sx + np.sin(np.pi/3)*sy_
G3 = (A3 @ B3 - B3 @ A3)/(2j*np.sin(np.pi/3))
w3 = np.linalg.eigvalsh(np.pi/2*G3)                                    # rho = e^(-pi G/2)/Z: Lücke g = pi
tau3 = (w3.max() - w3.min())/np.pi
jm3 = j_marken(np.pi/3)
E5, E3, E2 = E_ns(mp.mpf(5)), E_ns(mp.mpf(3)), E_ns(mp.mpf(2))
ZAHLEN["psym_minimalitaet"] = {"theta = pi/3: Im tau_Gesetz bei g = pi": tau3, "j(Marken) theta = pi/3": jm3,
                               "E(tau = 5i) (beta = 20 pi)": mp.nstr(E5, 10), "E(tau = 3i) (beta = 12 pi)": mp.nstr(E3, 10),
                               "E(tau = 2i) (P_sym mit g = 2 pi)": mp.nstr(E2, 10)}
check("Nichtzirkularität (Prämissen): A* legt beta = 4 pi bzw. die Saatlücke g = pi fest, unabhängig vom Winkel theta der "
      f"Binäroperationen; bei theta = pi/3 ist das Gesetz mit g = pi selbstdual (Im tau_Gesetz = {tau3:.6f}) und spannungsfrei, "
      f"die Marken haben aber j = {jm3:.1f} != 1728: A* impliziert P_sym nicht. Umgekehrt ist P_sym mit g = 2 pi (tau = 2i) nicht "
      f"spannungsfrei (E = {mp.nstr(E2, 6)}). A* enthält die Minimalität: beta = 20 pi (P_clock ohne Uhrenquotient, z94) gibt "
      f"tau = 5i mit E = {mp.nstr(E5, 6)}, beta = 12 pi (Delta^(-i/2) = R^-1) tau = 3i mit E = {mp.nstr(E3, 6)}; der Rahmen "
      "(z87: K, H) lässt jedes beta zu, enthält P_tn also nicht",
      abs(tau3 - 1) < 1e-9 and abs(jm3 - 1728) > 1 and abs(E2) > 0.1 and abs(E5) > 0.1 and abs(E3) > 0.1)

# Basis- gegen Überlagerungstorus
q_b = mp.exp(-2*mp.pi*2)
zb_prod = q_b**(-mp.mpf(1)/48)*mp.fprod((1 + q_b**(mp.mpf(m)/4))**2 for m in range(1, 400, 2))
zb_sekt = mp.exp(2*logZ_zwist(mp.mpf(2), mp.mpf(1)/4))
z_ns_i = mp.exp(logZ_ns(mp.mpf(1)))                                    # Z_NS(i) mit Casimir
cas_sektor = (mp.mpf(1)/4)**2/2 - mp.mpf(1)/24
Eb = lambda y: -mp.diff(lambda t: logZ_zwist(t, mp.mpf(1)/4), y)
yb = stationaer(Eb, 1.7)
sb_abw = abs(logZ_zwist(mp.mpf(4), mp.mpf(1)/4) - logZ_zwist(mp.mpf("0.25"), mp.mpf(1)/4))
ZAHLEN["basistorus"] = {"Casimir je Viertelsektor": mp.nstr(cas_sektor, 10), "stationaer tau_b/i": mp.nstr(yb, 20),
                        "S-Abweichung Basisinhalt (y = 4 gegen 1/4)": mp.nstr(sb_abw, 6)}
check("Basis- gegen Überlagerungstorus: auf dem Basistorus (Raumzyklus 2 pi in x-Einheiten, q_b = e^(-beta)) besteht der Inhalt "
      f"aus den zwei Viertelsektoren (Zwist 1/4, 3/4, Casimir je (1/4)^2/2 - 1/24 = {mp.nstr(cas_sektor, 6)}, zusammen -1/48); "
      f"Z_b(tau_b) = (theta[1/4, 0]/eta)^2(tau_b) = Z_NS(tau_b/2) (bei tau_b = 2i: {mp.nstr(abs(zb_prod - zb_sekt), 2)} und "
      f"{mp.nstr(abs(zb_prod - z_ns_i), 2)}). Stationär bei tau_b = {mp.nstr(yb, 15)} i, also wieder beta = 4 pi: A* wählt "
      "dasselbe beta unabhängig von der Normierung. Dort ist der Basistorus nicht selbstdual, und sein Inhalt ist nicht "
      f"S-invariant ({mp.nstr(sb_abw, 4)}): A*(a) zeichnet den Überlagerungstorus aus, genau die Antiperioden-Definition von "
      "tau_Gesetz",
      abs(zb_prod - zb_sekt) < mp.mpf(10)**-30 and abs(zb_prod - z_ns_i) < mp.mpf(10)**-30 and yb is not None
      and abs(yb - 2) < mp.mpf(10)**-25 and sb_abw > 0.1)
befund("A* ist nicht P_tn umbenannt: P_tn ist eine Aussage über die Geometrie des Gesetzestorus (tau = -1/tau), A*(b) über die "
       "Antwort des Inhalts (Casimir plus thermische Energie). Beide fallen genau für S-invarianten Inhalt zusammen, und den "
       "liefert A*(a); der Schluss benutzt weder P_tn noch P_clock noch P_sym. Er ist dieselbe Logik wie H6h4 für den "
       "E8-Nahttorus, jetzt auf dem Gesetzestorus von z94 mit dessen eigenem Inhalt geführt",
       all(CHECKS[k] for k in CHECKS if k.startswith(("Nichtzirkularität", "Basis-"))))
say(f"   (Laufzeit Teil 1: {time.time() - T_START:.0f} s)")

# =================================================================== Teil 2
say("\nTeil 2: Starrheit S(n)")


def eps_von(ms):
    return [1 if m % 4 == 1 else -1 for m in ms]


def paarbloecke(ms):
    """wie z94: Q_- -Randblöcke (Exponenten je Ladung) und mu = N/D gekürzt."""
    eps = eps_von(ms)
    rand = defaultdict(list)
    for st in itertools.product(((0, 0), (1, 1), (1, -1), (2, 0)), repeat=len(ms)):
        E = sum(m*nu for m, (nu, c) in zip(ms, st))
        qm = sum(e*nu for e, (nu, c) in zip(eps, st))
        rand[qm].append(E)
    mu = sum(e*2*yS**m/(1 + yS**m) for e, m in zip(eps, ms))
    N, D = sp.fraction(sp.cancel(sp.together(mu)))
    cf = lambda P: np.array([int(c) for c in reversed(sp.Poly(sp.expand(P), yS).all_coeffs())], dtype=np.int64)
    return {k: sorted(v) for k, v in rand.items()}, cf(N), cf(D), sp.expand(N), sp.expand(D)


def restklassen_kern(ms, Yval):
    """wie z94: Restklassensystem über Q bei Y = Yval; Kern = Unbekannte - Rang."""
    rand, _, _, Ne, De = paarbloecke(ms)
    qs = sorted(rand)
    zeilen = [[Fraction(0)]*len(qs) for _ in range(4)]
    for iq, q in enumerate(qs):
        for (e,), c in sp.Poly(sp.expand((q*De - Ne)*yS**(q % 4)), yS).terms():
            zeilen[e % 4][iq] += Fraction(int(c))*Fraction(Yval)**(e//4)
    M = flint.fmpq_mat(4, len(qs), [flint.fmpq(x.numerator, x.denominator) for r in zeilen for x in r])
    return len(qs) - M.rank()


def kern_formel(n):
    """Kern aus den Restklassen: Ladungen 2n + 1 aufeinanderfolgend; Klassen mod 4 mit k_r Werten; Z' hat 4 Restklassen-
    komponenten, jede Einzelklasse (k_r = 1) gibt eine Bedingung, jede Klasse mit k_r >= 3 trägt k_r - 2 verdeckte Richtungen."""
    k = Counter(q % 4 for q in range(2*n + 1))
    return 4 - sum(1 for r in range(4) if k[r] == 1) + sum(k[r] - 2 for r in range(4) if k[r] >= 3)


# --- 2.1 Normalform und Kernzahl
nf_ok = True
for ms in ((1, 3, 5), (3, 5, 7, 9), (1, 7, 9, 11)):
    eps = eps_von(ms)
    mus = [e*m for e, m in zip(eps, ms)]
    orig = sp.prod([(1 + sS**e*yS**m)**2 for e, m in zip(eps, ms)])
    nm = sum(1 for e in eps if e < 0)
    norm = sS**(-2*nm)*yS**(2*sum(m for e, m in zip(eps, ms) if e < 0))*sp.prod([(1 + sS*yS**mu)**2 for mu in mus])
    nf_ok &= sp.simplify(sp.expand(orig - norm)) == 0 and all((mu - 1) % 4 == 0 for mu in mus)
tupel = {2: list(itertools.combinations(range(1, 14, 2), 2)), 3: list(itertools.combinations(range(1, 14, 2), 3)),
         4: list(itertools.combinations(range(1, 12, 2), 4)), 5: [(1, 3, 5, 7, 9), (1, 3, 5, 7, 11), (3, 5, 7, 9, 11)]}
kern_ok, kern_tab = True, {}
for n, liste in tupel.items():
    ks = Counter(restklassen_kern(ms, Fraction(3, 7)) for ms in liste)
    kern_tab[n] = {"Tupel": len(liste), "Kerne": dict(ks), "Formel": kern_formel(n), "2n - 3": 2*n - 3}
    kern_ok &= set(ks) == {kern_formel(n)} and kern_formel(n) == 2*n - 3
ZAHLEN["kernzahl"] = kern_tab
formel_allg = all(kern_formel(n) == 2*n - 3 for n in range(2, 40))
check("Normalform (symbolisch für (1,3,5), (3,5,7,9), (1,7,9,11)): Teilchen-Loch auf den Paaren mit m = 3 mod 4 gibt Phi(y, s) "
      "= s^(-2 n_-) y^(2 sum_(eps = -1) m) prod (1 + s y^mu_j)^2 mit mu_j = eps_j m_j = 1 mod 4; mit X = ys, Y = y^4 ist Phi = "
      "prod (1 + X Y^k_j)^2, k_j = (mu_j - 1)/4, und die Uhr (Exponent = Ladung mod 4) ist eingebaut. Nichtentartung heißt "
      "genau: die k_j sind verschieden. Kernzahl aus den Restklassen (2n + 1 aufeinanderfolgende Ladungen, vier Klassen): "
      "4 - #Einzelklassen + sum_(k_r >= 3) (k_r - 2) = 2n - 3 für n = 2..39; exakt (flint) für "
      f"{', '.join(f'n = {n}: {v['Tupel']} Tupel Kern {list(v['Kerne'])}' for n, v in kern_tab.items())}",
      nf_ok and kern_ok and formel_allg)
befund("Lesart der Kernzahl: für n = 3 liegen alle drei freien Richtungen im Gesamtspektrum Z' selbst (keine Klasse mit drei "
       "Ladungen, eine Einzelklasse: die mittlere Ladung); ab n = 4 kommen verdeckte Umverteilungen innerhalb der Blöcke dazu, "
       "die Z' nicht sieht", kern_ok and Counter(q % 4 for q in range(7)) == Counter({0: 2, 1: 2, 2: 2, 3: 1}))

# --- 2.2 Partialbruchlemma
H1 = [lambda w: -4*w/(1 - w)**2, lambda w: (1 + 3*w)/(1 - w)**2, lambda w: -2*(1 + w)/(1 - w)**2, lambda w: (3 + w)/(1 - w)**2]
H2 = [lambda w: -2*w*(3*w + 5)/(1 - w)**3, lambda w: (3*w*w + 12*w + 1)/(1 - w)**3, lambda w: -(w*w + 12*w + 3)/(1 - w)**3,
      lambda w: 2*(5*w + 3)/(1 - w)**3]
zerl_ok = all(sum(u**r*H1[r](u**4) for r in range(4)) == u/(1 + u)**2 and sum(u**r*H2[r](u**4) for r in range(4)) == u/(1 + u)**3
              for u in (Fraction(1, 3), Fraction(5, 7), Fraction(-2, 9)))


def pb_system(ms, Y0):
    """4 x 2n-System für die Polkoeffizienten (A_j, B_j) bei Y = Y0 (exakt): Restklasse rho der Gleichung
    sum_j X^mu_j [A_j/(1 + X^mu_j)^2 + 2 B_j/(1 + X^mu_j)^3] = 0 auf der Kurve Y = X^4."""
    eps = eps_von(ms)
    mus = [e*m for e, m in zip(eps, ms)]
    ks = [(mu - 1)//4 for mu in mus]
    n = len(ms)
    M = sp.zeros(4, 2*n)
    for r in range(4):
        for j in range(n):
            w = Y0**mus[j]
            M[r, j] = sp.Rational(Y0**(r*ks[j])*H1[r](w))
            M[r, n + j] = sp.Rational(2*Y0**(r*ks[j])*H2[r](w))
    return M


pb_tab, pb_ok = {}, True
for ms in itertools.combinations(range(1, 16, 2), 3):
    eintrag = []
    for Y0 in (Fraction(3, 7), Fraction(5, 2)):
        M = pb_system(ms, Y0)
        ns = M.nullspace()
        V = sp.Matrix.hstack(*ns) if ns else sp.zeros(6, 0)
        minors = [V.extract([j, 3 + j], [0, 1]).det() if len(ns) == 2 else 0 for j in range(3)]
        eintrag.append((M.rank(), len(ns), all(mi != 0 for mi in minors)))
    pb_tab[str(ms)] = eintrag
    pb_ok &= all(e == (4, 2, True) for e in eintrag)
vier = {}
for ms in ((1, 3, 5, 7), (1, 3, 5, 9), (3, 5, 7, 9)):
    M = pb_system(ms, Fraction(3, 7))
    ns = M.nullspace()
    V = sp.Matrix.hstack(*ns)
    einzel = [V.extract([j, 4 + j], list(range(V.shape[1]))).rank() for j in range(4)]
    paare = [V.extract([i, j, 4 + i, 4 + j], [0, 1, 2, 3]).det() != 0 for i, j in itertools.combinations(range(4), 2)]
    vier[str(ms)] = {"Rang": M.rank(), "Polraum": len(ns), "Rang Einzelprojektion": einzel, "Paarpaare invertierbar": paare}
ZAHLEN["partialbruch_n3"] = {"Tripel": len(pb_tab), "alle (Rang 4, Polraum 2, drei Projektionen invertierbar)": pb_ok}
ZAHLEN["partialbruch_n4"] = vier
check("Partialbruchlemma (exakt, Restklassenzerlegung von u/(1+u)^2, u/(1+u)^3 geprüft): mit verschiedenen k_j hat F = W'/Phi "
      "die Form c(Y) + sum_j [A_j/(1 + X Y^k_j) + B_j/(1 + X Y^k_j)^2] (2n + 1 Unbekannte für 2n + 1 Blöcke), und die "
      "Ladungsidentität mit Uhr ist d_X F = 0 auf der Kurve Y = X^4: vier lineare Gleichungen nur für (A, B). Für alle 56 "
      "Tripel m_3 <= 15 bei Y = 3/7 und 5/2: Rang 4, Polraum der Dimension 2, und jede der drei Paarprojektionen (A_j, B_j) ist "
      f"invertierbar ({pb_ok}). Vier Paare (1,3,5,7), (1,3,5,9), (3,5,7,9): Polraum 4, Einzelprojektion Rang 2 (nicht "
      "injektiv), jede Projektion auf zwei Paare invertierbar "
      f"({all(all(v['Paarpaare invertierbar']) and v['Polraum'] == 4 for v in vier.values())})",
      zerl_ok and pb_ok and all(v["Rang"] == 4 and v["Polraum"] == 4 and all(v["Paarpaare invertierbar"])
                                and v["Rang Einzelprojektion"] == [2, 2, 2, 2] for v in vier.values()))
befund("Folgerung (Struktur für jedes Tripel verschiedener m, Invertierbarkeit exakt für die 56 Tripel m_3 <= 15): der Polteil "
       "hat in X höchstens den Grad 2n - 1, der oberste Block (Dimension 1) legt c = Y^c0 als Monom fest; spaltet EIN Paar ab "
       "((1 + X Y^k_j)^2 teilt W', also A_j = B_j = 0), verschwindet wegen der invertierbaren Projektion der ganze Polteil, "
       "W' = Y^c0 Phi, und S(3) gilt (ohne Positivität). Also S(3) <=> in jeder Lösung spaltet ein Paar ab; allgemein S(n) <=> "
       "Phi teilt W', und für n = 4 (drei Quadrupel geprüft) genügt das Abspalten zweier Paare. Die Nichtentartung trägt genau "
       "hier: bei gleichen k_j fallen die Pole zusammen, es gibt keine Paarblöcke",
       CHECKS[[k for k in CHECKS if k.startswith("Partialbruchlemma")][0]])

# --- 2.3 ganzzahlige Starrheit und reelle Relaxation


def aufbau(ms, lo, hi):
    """Q_- -Randproblem wie z94 (Relaxation des vollen Problems): w_(q,e) >= 0, e = q mod 4, Blocksummen, D M' = N Z'."""
    rand, Nc, Dc, _, _ = paarbloecke(ms)
    keys = sorted(rand)
    var = [(q, e) for q in keys for e in range(lo, hi + 1) if (e - q) % 4 == 0]
    idx = {v: i for i, v in enumerate(var)}
    G = max(len(Nc), len(Dc)) - 1
    A = lil_matrix((len(keys) + (hi + G - lo + 1), len(var)))
    lb, ub, r = [], [], 0
    for q in keys:
        for e in range(lo, hi + 1):
            if (q, e) in idx:
                A[r, idx[(q, e)]] = 1
        lb.append(len(rand[q]))
        ub.append(len(rand[q]))
        r += 1
    for e2 in range(lo, hi + G + 1):
        for (q, e), i in idx.items():
            d, val = e2 - e, 0
            if 0 <= d < len(Dc):
                val += q*Dc[d]
            if 0 <= d < len(Nc):
                val -= Nc[d]
            if val:
                A[r, i] = val
        lb.append(0)
        ub.append(0)
        r += 1
    rand_q = keys[0]
    fix = lil_matrix((1, len(var)))
    fix[0, idx[(rand_q, rand[rand_q][0])]] = 1                         # Randblock (Dimension 1) am Originalexponenten
    return vstack([csr_matrix(A[:r]), csr_matrix(fix)]), lb + [1], ub + [1], var, rand


def l1_ganz(ms, lo, hi, tl=120):
    """max sum_e |Z'_e - Z_e| über alle ganzzahligen Lösungen mit festem Randblock (Binärhilfsvariablen wie z94)."""
    A, lb, ub, var, rand = aufbau(ms, lo, hi)
    nv = len(var)
    Z0 = Counter(e for v in rand.values() for e in v)
    lev = list(range(lo, hi + 1))
    L, U = len(lev), 4**len(ms)
    B = lil_matrix((2*L, nv + 2*L))
    lb2, ub2 = [], []
    for li, e in enumerate(lev):
        for i, (q, ee) in enumerate(var):
            if ee == e:
                B[2*li, i] = -1
                B[2*li + 1, i] = 1
        B[2*li, nv + li], B[2*li, nv + L + li] = 1, -U
        B[2*li + 1, nv + li], B[2*li + 1, nv + L + li] = 1, U
        lb2 += [-np.inf, -np.inf]
        ub2 += [-Z0.get(e, 0), Z0.get(e, 0) + U]
    AA = vstack([csr_matrix(np.hstack([A.toarray(), np.zeros((A.shape[0], 2*L))])), csr_matrix(B)])
    c = np.concatenate([np.zeros(nv), -np.ones(L), np.zeros(L)])
    integ = np.concatenate([np.ones(nv), np.zeros(L), np.ones(L)])
    hib = np.concatenate([np.full(nv, np.inf), np.full(L, np.inf), np.ones(L)])
    res = milp(c, constraints=LinearConstraint(AA, lb + lb2, ub + ub2), integrality=integ, bounds=Bounds(np.zeros(len(c)), hib),
               options={"time_limit": tl})
    return (None if res.x is None else -res.fun), int(res.status)


def lp_reell(ms, lo, hi, versuche=4):
    """Reelle Relaxation mit festem Randblock: größte L1-Abweichung vom Original über zufällige lineare Ziele (Ecken)."""
    A, lb, ub, var, rand = aufbau(ms, lo, hi)
    Z0 = Counter((q, e) for q, Es in rand.items() for e in Es)
    w0 = np.array([Z0.get(v, 0) for v in var], dtype=float)
    best = 0.0
    for _ in range(versuche):
        cc = rng.normal(size=len(var))
        res = milp(-cc, constraints=LinearConstraint(A, lb, ub), integrality=np.zeros(len(var)),
                   bounds=Bounds(np.zeros(len(var)), np.full(len(var), np.inf)))
        if res.x is not None:
            best = max(best, float(np.abs(res.x - w0).sum()))
    return best


def lp_aussen(ms, lo, hi):
    """Reelle Relaxation: größte Masse außerhalb des Originalträgers und deren Reichweite."""
    A, lb, ub, var, rand = aufbau(ms, lo, hi)
    Z0 = Counter((q, e) for q, Es in rand.items() for e in Es)
    c = np.array([-(0.0 if v in Z0 else 1.0) for v in var])
    res = milp(c, constraints=LinearConstraint(A, lb, ub), integrality=np.zeros(len(var)),
               bounds=Bounds(np.zeros(len(var)), np.full(len(var), np.inf)))
    tr = sorted({e for (q, e), x in zip(var, res.x) if x > 1e-9})
    return -res.fun, (tr[0], tr[-1])


t2 = time.time()
ganz, reell3 = {}, {}
for ms in itertools.combinations(range(1, 16, 2), 3):
    lo, hi = -24, 2*sum(ms) + 24
    ganz[ms] = l1_ganz(ms, lo, hi)
    reell3[ms] = lp_reell(ms, lo, hi, versuche=3)
reell2 = {ms: lp_reell(ms, -24, 2*sum(ms) + 24, versuche=3) for ms in itertools.combinations(range(1, 16, 2), 2)}
aussen = {pad: lp_aussen((1, 3, 5), -pad, 18 + pad) for pad in (24, 36)}
ganz_ok = all(st == 0 and v is not None and v < 0.5 for v, st in ganz.values())
nicht_starr3 = sorted(ms for ms, v in reell3.items() if v > 1e-6)


def blockdim_symmetrisch(ms):
    """Volle Blöcke (Q_-, Q_+): dim(q, p) = dim(q, -p)? Dann hebt die proportionale Aufteilung jede reelle Q_- -Lösung auf
    das volle Blockproblem (Blockdimensionen und <Q_+> = 0 bleiben erfüllt)."""
    eps = eps_von(ms)
    dims = Counter((sum(e*nu for e, (nu, c) in zip(eps, st)), sum(c for nu, c in st))
                   for st in itertools.product(((0, 0), (1, 1), (1, -1), (2, 0)), repeat=len(ms)))
    return all(dims[(q, p)] == dims[(q, -p)] for (q, p) in dims)


hebung = all(blockdim_symmetrisch(ms) for ms in itertools.combinations(range(1, 16, 2), 3))
ZAHLEN["ganzzahlig_56_tripel"] = {"alle max L1 = 0": ganz_ok, "Fenster": "[-24, 2 sum m + 24]", "s": round(time.time() - t2, 1)}
ZAHLEN["reell_n2_max_abweichung"] = max(reell2.values())
ZAHLEN["reell_n3_nicht_starr"] = {"Anzahl": len(nicht_starr3), "von": len(reell3), "Tripel": [str(m) for m in nicht_starr3]}
ZAHLEN["reell_(1,3,5)_masse_ausserhalb"] = {f"Fenster +-{k}": {"Masse": round(v[0], 3), "Träger": v[1]} for k, v in aussen.items()}
check(f"ganzzahlige Starrheit, alle 56 ungeraden Tripel m_1 < m_2 < m_3 <= 15 (Fenster [-24, 2 sum m + 24], Randblock fest, HiGHS "
      f"mit Optimalitätsbeweis): max L1(Z' - Z) = 0 für jedes ({ganz_ok}); z94 hatte sechs Tripel", ganz_ok)
check(f"reelle Relaxation (gleiche Bedingungen, nur nichtnegativ reell): für n = 2 starr, alle {len(reell2)} Paare m <= 15 (größte "
      f"Abweichung {max(reell2.values()):.1e}, wie der lineare Beweis von S(2) verlangt); für n = 3 nicht: {len(nicht_starr3)} "
      f"von 56 Tripeln haben reelle nichtnegative Lösungen mit richtigen Blockdimensionen, die keine Spektren-Mischung sind "
      f"(etwa (1,3,5): {reell3[(1, 3, 5)]:.1f}); bei (1,3,5) erreicht die Masse außerhalb des Originalträgers jeden Fensterrand "
      f"({', '.join(f'+-{k}: {v[0]:.1f} auf {v[1]}' for k, v in aussen.items())}): reell gibt es keine Spannenschranke. Die "
      f"Blockdimensionen sind in Q_+ symmetrisch (alle 56 Tripel: {hebung}), also hebt die proportionale Aufteilung jede reelle "
      "Q_- -Lösung auf das volle Blockproblem mit <Q_+> = 0: das Hindernis gilt für den ganzen Blockweg, nicht nur für die "
      "Relaxation",
      max(reell2.values()) < 1e-6 and len(nicht_starr3) >= 10 and (1, 3, 5) in nicht_starr3 and hebung
      and aussen[24][1] == (-24, 42) and aussen[36][1][0] <= -30 and aussen[36][1][1] >= 50)
befund("Hindernis für den allgemeinen Satz (gescopt negativ für den linearen Blockweg der z94-Reduktion): die Starrheit ab n = 3 "
       "ist arithmetisch. Ein Beweis aus Kern 2n - 3, Nichtentartung und Positivität allein (also jede Induktion über den "
       "linearen Kern, auch über das Abspalten eines Paares) kann nicht gelingen, weil schon für (1,3,5) reelle Blocklösungen "
       "ohne Spannenschranke existieren; erst die Ganzzahligkeit der Eigenwertvielfachheiten schließt sie aus. Das entartete "
       "Gegenbeispiel (z94) ist der Fall, in dem auch die Ganzzahligkeit nicht mehr rettet. Offen bleibt ein Weg über die "
       "feineren Einteilchendaten (Eigenvektoren), den die Blockreduktion nicht sieht",
       (1, 3, 5) in nicht_starr3 and ganz_ok and hebung)
befund("S(3) für allgemeine ungerade m_1 < m_2 < m_3 bleibt offen: bewiesen ist die Reduktion 'S(3) <=> ein Paar spaltet ab' "
       "(exakt für alle 56 Tripel m_3 <= 15, als Struktur für jedes Tripel verschiedener m), ganzzahlig vollständig sind alle "
       "56 Tripel im Fenster; eine vollständige endliche Reduktion fehlt, weil ganzzahlig keine Spannenschranke bewiesen ist und "
       "reell keine gilt", True)
say(f"   (Laufzeit bis Teil 2: {time.time() - T_START:.0f} s)")

# =================================================================== Teil 3
say("\nTeil 3: Bilanz für N8")
STATUS = {
    "P_tn aus A*": ("E unter A* und Rahmen", "Gesetzestorus = Nahttorus (A, A), Inhalt theta_3/eta; A*(a): Casimir -1/24; A*(b): "
                    "(theta_2^4 - theta_4^4)(tau_Gesetz) = 0 <=> tau_Gesetz = i auf allen Zustandstori"),
    "Gesetzes-Entsprechung von E6": ("E", "E6 = (1/2)(theta_4^4 - theta_2^4)(theta_2^4 + theta_3^4)(theta_3^4 + theta_4^4): "
                                         "Produkt der drei geraden Spinstruktur-Spannungen; (A, A)-Faktor = Gesetz"),
    "Nichtzirkularität": ("E", "spannungsfrei != zeitneutral für nicht S-invarianten Inhalt (Zwist, Casimir); A* lässt P_sym "
                               "frei; Rahmen lässt jedes beta zu"),
    "Minimalität": ("E", "beta = 12 pi, 20 pi nicht spannungsfrei; der Uhrenquotient wird entbehrlich"),
    "S(3) <=> ein Paar spaltet ab": ("E (56 Tripel exakt, Struktur allgemein)", "Polraum 2, Paarprojektionen invertierbar"),
    "S(3) ganzzahlig": ("E (im Fenster)", "alle 56 Tripel m_3 <= 15, max L1 = 0"),
    "S(n) allgemein": ("O", "Hindernis: reelle Relaxation ab n = 3 nicht starr, Beweis muss arithmetisch sein"),
    "K19 Existenz der Naht": ("L", "A* setzt die Naht voraus; neu hergeleitet ist nur ihre Form: Raumzyklus = KMS-Zyklus"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
PRAEMISSEN = {
    "P_sym": "zwei Binäroperationen mit Vorzeichensymmetrie in derselben irreduziblen Darstellung -> Cl_2 = CAR, Marken als "
             "Quadrat, tau_Marken = i (z80, z87); unabhängig von A*",
    "A* (geteilt mit H6h4)": "hintergrundfreie Naht: große Diffeomorphismen sind Symmetrien (a), die Form erfüllt <T> = 0 (b); "
                            "trägt schon E8/GSO und tau = i der Naht, jetzt auch P_tn (tau_Gesetz = i, beta = 4 pi, g = pi)",
    "Rahmen": "Hardy/Szegő-Naht (Modulfluss = Drehung, Quartermoden x = m/4 = NS der Überlagerung), Ladungssuperauswahl, "
              "Nichtentartung (verschiedene k_j)",
    "P_tn, P_clock": "Folgen (P_clock <=> P_tn bis auf die Benennung Original/Partner, z94)",
}
ZAHLEN["praemissen"] = PRAEMISSEN
ZAHLEN["praemissen_anzahl"] = {"ursprungseigen": 1, "geteilt mit der Naht (A*)": 1, "Rahmen": 1,
                               "unabhaengig": "P_sym und A* (theta = pi/3 gegen g = 2 pi)"}
check("N8-Kriterien: P_tn folgt aus A* nichtzirkulär (Inhalt, Casimir und Prämissen geprüft); die Ursprungsprämissen sinken "
      "von zwei (P_sym, P_tn) auf eine eigene (P_sym) plus die schon vorhandene Nahtprämisse A*; S(3) ist auf das Abspalten "
      "eines Paares reduziert und im Fenster für 56 Tripel ganzzahlig starr, der allgemeine Satz hat ein benanntes Hindernis",
      CHECKS[[k for k in CHECKS if k.startswith("Nullstellen")][0]]
      and CHECKS[[k for k in CHECKS if k.startswith("Nichtzirkularität (Inhalt)")][0]]
      and CHECKS[[k for k in CHECKS if k.startswith("Nichtzirkularität (Prämissen)")][0]]
      and CHECKS[[k for k in CHECKS if k.startswith("Partialbruchlemma")][0]]
      and CHECKS[[k for k in CHECKS if k.startswith("ganzzahlige Starrheit")][0]])

urteil = ("N8 trägt in gescopter Form. (1) P_tn folgt aus A*: Unter dem Rahmen (Modulfluss = Drehung der Hardy-Naht) ist der "
          "Gesetzestorus aus z94, Raumzyklus = Antiperiode 4 pi/beta, Zeitzyklus = KMS, der euklidische Nahttorus mit "
          "Spinstruktur (A, A) und tau_Gesetz = i beta/(4 pi); seine Zustandssumme ist exakt das NS-Zeichen theta_3/eta eines "
          "komplexen Fermions. A*(a) macht diesen Inhalt S-invariant und legt damit die Casimirenergie -1/24 fest (das normierte "
          "Gesetz allein sieht sie nicht). A*(b) <T> = 0 ist d log(theta_3/eta)/d tau = (pi i/12)(theta_2^4 - theta_4^4) = 0, "
          "der (A, A)-Faktor von E6 = (1/2)(theta_4^4 - theta_2^4)(theta_2^4 + theta_3^4)(theta_3^4 + theta_4^4), dem Produkt "
          "der Spannungen der drei geraden Spinstrukturen. Auf den Zustandstori tau = iy ist die einzige Nullstelle y = 1 "
          "(S-Symmetrie legt sie auf den Fixpunkt, Var > 0 macht sie eindeutig), in der Halbebene genau Gamma_theta i "
          "(Hauptmodul lambda, lambda(i) = 1/2); auch E6(tau_Gesetz) = 0 wählt genau i. Also tau_Gesetz = i, beta = 4 pi, "
          "Saatlücke pi: P_tn und mit z94 P_clock samt Minimalität (beta = 12 pi und 20 pi sind nicht spannungsfrei). Der "
          "Schluss ist nicht zirkulär: Spannungsfreiheit und Zeitneutralität sind verschiedene Prädikate (verdrillter Inhalt ist "
          f"bei y* = {float(z8['y*']):.3f} spannungsfrei, am selbstdualen Torus nicht; falscher Casimir verschiebt oder tilgt den "
          "Punkt), sie "
          "fallen genau für S-invarianten Inhalt zusammen, und den liefert A*(a); der Rahmen lässt jedes beta zu; der "
          "Basistorus gibt dasselbe beta, ist aber nicht selbstdual, A*(a) zeichnet den Überlagerungstorus aus. Neu gegenüber "
          "H6h4, wo derselbe Schluss für den E8-Nahttorus stand, ist die Führung auf dem Gesetzestorus von z94 mit dessen "
          "eigenem Inhalt, die Gesetzes-Entsprechung von E6 und der Nachweis der Nichtzirkularität. A* impliziert P_sym nicht "
          f"(theta = pi/3: spannungsfrei, j = {jm3:.1f}). (2) S(n): Normalform mit mu_j = eps_j m_j = 4 k_j + 1, Nichtentartung = "
          f"verschiedene k_j; Kernzahl 2n - 3 aus den Restklassen (exakt für {sum(v['Tupel'] for v in kern_tab.values())} Tupel "
          "n = 2..5). Neu ist das Partialbruchlemma: "
          "Der Zusatzkern ist der Polteil von F = W'/Phi, der oberste Block macht c zum Monom, und für n = 3 ist jede "
          "Paarprojektion des zweidimensionalen Polraums invertierbar (alle 56 Tripel m_3 <= 15): S(3) gilt genau dann, wenn in "
          "jeder Lösung ein Paar abspaltet, dann schließt S(2) ohne Positivität; für n = 4 genügen zwei Paare (drei Quadrupel "
          "geprüft). Ganzzahlig sind "
          "alle 56 Tripel im Fenster starr. Ein allgemeiner Beweis von S(3) fehlt, und es gibt ein benanntes Hindernis: Die "
          f"reelle Relaxation ist für n = 2 starr, für n = 3 nicht ({len(nicht_starr3)} von 56 Tripeln im Fenster +-24, bei (1,3,5) ohne "
          "Spannenschranke; die Blockdimensionen sind in Q_+ symmetrisch, das gilt also für den ganzen Blockweg). Die Starrheit "
          "ab n = 3 ist arithmetisch, ein Beweis auf dem Blockweg muss die Ganzzahligkeit der "
          "Eigenwertvielfachheiten benutzen; eine Induktion über den linearen Kern mit Nichtentartung und Positivität allein "
          "kann nicht schließen. (3) Der Ursprung hat eine eigene Prämisse, P_sym, dazu die Nahtprämisse A*, die H6h4 schon "
          "für E8/GSO und tau = i braucht, und den Rahmen (Hardy-Naht, Ladungssuperauswahl, Nichtentartung). P_sym und A* "
          "sind voneinander unabhängig; formal bleiben zwei Prämissen, ursprungseigen ist nur noch eine. P_free ist Folge für "
          "n = 2 allgemein und für 56 Tripel im Fenster. K19 bleibt L: A* setzt die Naht voraus und leitet nur ihre Form "
          "her (Raumzyklus so lang wie der KMS-Zyklus), nicht ihre Existenz.")
say("\nUrteil\n   " + urteil)
say(f"\nLaufzeit {time.time() - T_START:.0f} s")
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z107_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
