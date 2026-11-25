"""Z94: Starrheit allgemein und P_clock aus P_sym (Folgeauftrag F6).

Auftrag F6 (GESAMTSTAND-UND-STRATEGIEN.md §10 Punkt 4, §11 Punkt 6): (1) den Gaußzustand aus der Uhrenidentität für beliebig
viele Fermipaare beweisen oder ein Gegenbeispiel finden; (2) prüfen, ob P_clock aus P_sym über den Tausch A <-> B folgt
(tau -> -1/tau bei tau = i, Zeitneutralität aus H6 K2); (3) zählen, wie viele Ursprungsprämissen bleiben.
Leitlinie 'Rekonstruktion statt Auswahl'; gelöste Bausteine werden nicht geöffnet, nichts wird befördert. Status E/C/L/O.
Bestand: z87 (Saatfamilie, CAR, Starrheit in Abschneidungen, P_sym, P_clock), z80 Teil 8 und z71 Teil 3 (zwei
Binäroperationen, Cl_2, Spurzustand), n1 (zeitneutrale Naht: beide OS-Rekonstruktionen geben dasselbe KMS-System), z88
(Leiter aus der Hardy-Invarianz). Neu geprüft: die Behauptung der elften Außenanalyse (gleiche Saaten, Bruch bei N = 6).

Teil 1  Starrheit allgemein. Satz S(n): Für Fermipaare m_1 < .. < m_n (ungerade) erzwingen Uhrenidentität, Ladungs-
        superauswahl je Quarter und Einteilchendaten den quasifreien Zustand. Allgemeine Reduktion (jedes n): Gitterspektrum
        je Ladungsblock, Ladungsidentitäten über Z[y] (e^pi transzendent), Max-Entropie/Subadditivität; S(n) folgt, sobald
        die Ladungsidentitäten globale Isospektralität erzwingen. Erzeugende Funktion Phi(y, s) = prod (1 + s^eps y^m)^2,
        Uhr als Phi(iy, s) = Phi(y, is). Neu: zwei Paare für alle ungeraden m_1 < m_2 und alle Exponenten bewiesen
        (Restklassen über Q(Y), Y = y^4, universelle Matrix bei Y = 1 vom Rang 4, Randblock der Dimension 1). Für n >= 3 hat
        das Restklassensystem einen (2n - 3)-dimensionalen Kern; ganzzahlige Suche (HiGHS) für drei und vier Paare;
        Lückenlemma; Kontrolle ohne Viertelverfeinerung. Außenbehauptung (gleiche Saaten): exakt geprüft.
Teil 2  P_clock aus P_sym: Tausch A <-> B auf dem Saat-Qubit (Algebra und Gesetz) und auf dem Original (Partner); der
        komplexe Modulzeit-Torus des Gesetzes (reelle Antiperiode gegen KMS-Antiperiode), Zeitneutralität im Sinn von n1;
        was folgt, was Prämisse bleibt, Unabhängigkeit von P_sym und Zeitneutralität.
Teil 3  Bilanz für F6: Prämissenzählung, Status, K19, Urteil.
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
import scipy.linalg as sla
import sympy as sp
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import csr_matrix, lil_matrix, vstack

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
mp.mp.dps = 50
rng = np.random.default_rng(94)
PI4 = 4*np.pi
YN = np.exp(-np.pi)
T_START = time.time()
yS, sS, tS = sp.symbols("y s t")


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def eps_von(ms):
    return [1 if m % 4 == 1 else -1 for m in ms]


def paarbloecke(ms):
    """Anregungen je Paar: (nu, c) mit nu = Teilchen + Loch, c = Teilchen - Loch; Q_- = Q_a - Q_b = sum eps nu,
    Q_+ = Q_a + Q_b = sum c, Exponent E = sum m nu (Gewicht y^E, y = e^-pi)."""
    eps = eps_von(ms)
    voll, rand = defaultdict(list), defaultdict(list)
    for st in itertools.product(((0, 0), (1, 1), (1, -1), (2, 0)), repeat=len(ms)):
        E = sum(m*nu for m, (nu, c) in zip(ms, st))
        qm = sum(e*nu for e, (nu, c) in zip(eps, st))
        qp = sum(c for nu, c in st)
        voll[(qm, qp)].append(E)
        rand[qm].append(E)
    mu = sum(e*2*yS**m/(1 + yS**m) for e, m in zip(eps, ms))
    N, D = sp.fraction(sp.cancel(sp.together(mu)))
    cf = lambda P: np.array([int(c) for c in reversed(sp.Poly(sp.expand(P), yS).all_coeffs())], dtype=np.int64)
    return ({k: sorted(v) for k, v in voll.items()}, {k: sorted(v) for k, v in rand.items()}, cf(N), cf(D), sp.expand(N),
            sp.expand(D))


# =================================================================== Teil 1
say("Teil 1: Starrheit für beliebig viele Fermipaare")

# --- 1.1 Formulierung: erzeugende Funktion, Uhr, Ladungsmittel
ms4 = (1, 3, 5, 7)
voll4, rand4, _, _, _, _ = paarbloecke(ms4)
Phi_enum = sp.expand(sum(sS**q*yS**E for q, Es in rand4.items() for E in Es))
Phi_prod = sp.expand(sp.prod([(1 + sS**e*yS**m)**2 for e, m in zip(eps_von(ms4), ms4)]))
Psi = sp.prod([(1 + sS**e*tS*yS**m)*(1 + sS**e*yS**m/tS) for e, m in zip(eps_von(ms4), ms4)])
mu4 = sum(e*2*yS**m/(1 + yS**m) for e, m in zip(eps_von(ms4), ms4))
uhr_fkt = sp.expand(Phi_prod.subs(yS, sp.I*yS) - Phi_prod.subs(sS, sp.I*sS)) == 0
mu_ok = sp.simplify(sp.diff(Phi_prod, sS).subs(sS, 1)/Phi_prod.subs(sS, 1) - mu4) == 0
qp_ok = sp.simplify(sp.diff(Psi, tS).subs({tS: 1, sS: 1})) == 0
gitter_ok = all((E - q) % 4 == 0 for q, Es in rand4.items() for E in Es)
check("Formulierung (vier Paare m = 1, 3, 5, 7, 256 Fockzustände): die Blockgewichte sind Phi(y, s) = prod_j (1 + s^eps_j y^m_j)^2 "
      "mit s^(Q_a - Q_b), eps_j = +-1 für m_j = +-1 mod 4; die Uhr R = e^(iK_0/2) ist die Funktionalgleichung Phi(iy, s) = "
      "Phi(y, is) (Exponent = Q_- mod 4 auf jedem Block); die Einteilchendaten geben <Q_a - Q_b> = d_s log Phi|_(s=1) = "
      "mu(y) = sum_j eps_j 2y^m_j/(1 + y^m_j) und <Q_a + Q_b> = 0 (symbolisch)",
      sp.expand(Phi_enum - Phi_prod) == 0 and uhr_fkt and mu_ok and qp_ok and gitter_ok)


# Prüfsystem zwei Paare (wie z87) für die Richtung Uhrenidentität <-> Gitterspektrum
def jw(n):
    I2_, Zp_, sm_ = np.eye(2), np.diag([1.0, -1.0]), np.array([[0.0, 1.0], [0.0, 0.0]])
    ops = []
    for j in range(n):
        mats = [Zp_]*j + [sm_] + [I2_]*(n - j - 1)
        op = mats[0]
        for M in mats[1:]:
            op = np.kron(op, M)
        ops.append(op)
    return ops


X4 = [0.25, -0.75, 0.75, -0.25]
SEK4 = [0, 0, 1, 1]
C4 = jw(4)
occ4 = np.array([np.diag(c.T @ c) for c in C4]).T.round().astype(int)
sea4 = (np.array(X4) < 0).astype(int)
E4 = ((occ4 - sea4) @ (4*np.array(X4))).round().astype(int)
Qa4 = ((occ4 - sea4)*(np.array(SEK4) == 0)).sum(axis=1)
Qb4 = ((occ4 - sea4)*(np.array(SEK4) == 1)).sum(axis=1)
R4 = np.diag(1j**((Qa4 - Qb4) % 4))
rho_qf = np.diag(YN**E4.astype(float))
rho_qf /= np.trace(rho_qf)
bl4 = defaultdict(list)
for s_ in range(16):
    bl4[(int(Qa4[s_]), int(Qb4[s_]))].append(s_)


def gitterzustand(verschiebung=0.0, mpd=False):
    """Ladungsblockdiagonal, Eigenwerte y^e mit e = Q_a - Q_b + 4k je Block, zufällige reelle Eigenvektoren im Block;
    ein Exponent im Block (0, 0) wird um 'verschiebung' verstellt. mpd: 50-stellige Matrix (Eigenwerte bis y^16)."""
    rho = mp.zeros(16, 16) if mpd else np.zeros((16, 16))
    for key, idx in bl4.items():
        qm = key[0] - key[1]
        ex = [qm + 4*int(rng.integers(-1, 2)) for _ in idx]
        ex = [e + (verschiebung if k == 0 and key == (0, 0) else 0.0) for k, e in enumerate(ex)]
        Qm, _ = np.linalg.qr(rng.normal(size=(len(idx), len(idx))))
        for a, ia in enumerate(idx):
            for b, ib in enumerate(idx):
                if mpd:
                    rho[ia, ib] = mp.fsum(mp.mpf(Qm[a, k])*mp.mpf(Qm[b, k])*mp.exp(-mp.pi*mp.mpf(ex[k])) for k in range(len(idx)))
                else:
                    rho[ia, ib] = sum(Qm[a, k]*Qm[b, k]*YN**ex[k] for k in range(len(idx)))
    return rho


def uhr_rest_mp(rho):
    """rho^(-i/2) R^-1 über die 50-stellige Spektralzerlegung; Rest = Abstand von einem Vielfachen der Eins."""
    w, V = mp.eigsy(rho)
    U = mp.zeros(16, 16)
    for k in range(16):
        ph = mp.exp(-0.5j*mp.log(w[k]))
        for a in range(16):
            for b in range(16):
                U[a, b] += ph*V[a, k]*V[b, k]
    Rinv = [mp.mpc(1)/mp.mpc(complex(R4[a, a])) for a in range(16)]
    U = mp.matrix([[U[a, b]*Rinv[b] for b in range(16)] for a in range(16)])
    return float(max(abs(U[a, b] - (U[0, 0] if a == b else 0)) for a in range(16) for b in range(16)))


reste = {v: max(uhr_rest_mp(gitterzustand(v, mpd=True)) for _ in range(3)) for v in (0.0, 4.0, 2.0, 0.5)}
ZAHLEN["uhr_rest_gitterzustaende"] = reste
check("Uhrenidentität <-> Gitterspektrum (jedes n; hier zwei Paare): ein ladungsblockdiagonaler Zustand mit Eigenwerten y^e, "
      "e = Q_a - Q_b mod 4 je Block und beliebigen Eigenvektoren erfüllt rho^(-i/2) R^-1 ~ 1 (Rest "
      f"{reste[0.0]:.0e}, auch nach Verschiebung eines Exponenten um 4: {reste[4.0]:.0e}); Verschiebung um 2 oder 1/2 "
      f"verletzt sie (Rest {reste[2.0]:.2f}, {reste[0.5]:.2f})", reste[0.0] < 1e-9 and reste[4.0] < 1e-9 and reste[2.0] > 0.1
      and reste[0.5] > 0.1)

# Max-Entropie-Form: K_0 ist Einteilchenoperator, also S(rho'||rho_qf) = S(rho_qf) - S(rho') bei gleichen Besetzungen
N4 = [np.diag(c.T @ c) for c in C4]
c4 = np.array([n @ np.diag(rho_qf).real for n in N4])


def entropie(rho):
    w = np.linalg.eigvalsh((rho + rho.conj().T)/2)
    w = w[w > 1e-300]
    return float(-np.sum(w*np.log(w)))


def relent_diag(p, q):
    return float(np.sum(p*(np.log(p) - np.log(q))))


Dd = (N4[0] - c4[0])*(N4[3] - c4[3])
p_qf = np.diag(rho_qf).real
p_l = p_qf*(1 + 0.5*Dd)
rho_l = np.diag(p_l)
c_x = np.array([0.1, 0.8, 0.2, 0.9])
p_x = np.prod(np.where(occ4 == 1, c_x, 1 - c_x), axis=1)
besetz_x = np.array([n @ p_x for n in N4])
d_l = relent_diag(p_l, p_qf) - (entropie(rho_qf) - entropie(rho_l))
d_x = relent_diag(p_x, p_qf) - (entropie(rho_qf) - entropie(np.diag(p_x)))
check("Max-Entropie-Form der Subadditivität (jedes n): K_0 = dGamma(4 pi x) ist ein Einteilchenoperator, also gilt bei gleichen "
      f"Besetzungen S(rho'||rho_qf) = S(rho_qf) - S(rho') >= 0 (z87-Deformation: Abweichung {abs(d_l):.0e}); ein Zustand mit "
      f"anderen Besetzungen verletzt die Identität (Abweichung {abs(d_x):.2f}). Globale Isospektralität gibt S(rho') = S(rho_qf), "
      "also rho' = rho_qf", abs(d_l) < 1e-12 and abs(d_x) > 1e-3 and np.abs(besetz_x - c4).max() > 1e-3)
befund("Satz S(n) ist für jedes n auf eine kombinatorische Aussage reduziert: Uhrenidentität + Superauswahl geben Blockgewichte "
       "W'_B(y) in Z_>=0[y, 1/y] mit Exponenten = Q_- mod 4; die Ladungsmittel werden wegen der Transzendenz von e^pi zu "
       "Identitäten D(y) sum_B Q_-(B) W'_B = N(y) sum_B W'_B und sum_B Q_+(B) W'_B = 0 (mu = N/D gekürzt); S(n) gilt, sobald "
       "diese Identitäten Z'(y) = y^(4k) Z(y) erzwingen (dann Max-Entropie). Lückenlemma: Eine Lücke > G = max(deg N, deg D) "
       "zwischen besetzten Exponenten zerlegt jede Lösung in zwei Teillösungen", True)

# --- 1.2 zwei Paare: Beweis für alle Exponenten
say("\n   1.2 Zwei Paare: Restklassen über Q(Y), Y = y^4")
I_ = sp.I
universal = {}
for n in (1, 2, 3, 4):
    qs_ = range(-n, n + 1)
    M = sp.Matrix([[q for q in qs_], [(-1)**q for q in qs_], [sp.expand((q - n*I_)*I_**q) for q in qs_],
                   [sp.expand((q + n*I_)*(-I_)**q) for q in qs_]])
    universal[n] = int(M.rank())
ZAHLEN["universelle_Y1_matrix_rang"] = universal


def restklassen(ms, Yval):
    """Gleichung sum_q (q D(y) - N(y)) y^(q mod 4) P_q(Y) = 0 nach Restklassen von y mod 4 zerlegt, exakt bei Y = Yval."""
    _, rand, _, _, Ne, De = paarbloecke(ms)
    qs = sorted(rand)
    zeilen = [[Fraction(0)]*len(qs) for _ in range(4)]
    for iq, q in enumerate(qs):
        for (e,), c in sp.Poly(sp.expand((q*De - Ne)*yS**(q % 4)), yS).terms():
            zeilen[e % 4][iq] += Fraction(int(c))*Fraction(Yval)**(e//4)
    M = flint.fmpq_mat(4, len(qs), [flint.fmpq(x.numerator, x.denominator) for r in zeilen for x in r])
    P = [sum(Fraction(Yval)**((E - q % 4)//4) for E in rand[q]) for q in qs]
    Pv = flint.fmpq_mat(len(qs), 1, [flint.fmpq(x.numerator, x.denominator) for x in P])
    return M.rank(), len(qs), (M*Pv) == flint.fmpq_mat(4, 1), [len(rand[q]) for q in qs]


zwei_ok, zwei_n = True, 0
for m1 in range(1, 32, 2):
    for m2 in range(m1 + 2, 32, 2):
        r1, nq, kern1, dims = restklassen((m1, m2), 1)
        r2, _, kern2, _ = restklassen((m1, m2), Fraction(3, 7))
        zwei_ok &= r1 == 4 and r2 == 4 and nq == 5 and kern1 and kern2 and dims[0] == 1 and dims[-1] == 1
        zwei_n += 1
ZAHLEN["zwei_paare_gepruefte_abschneidungen"] = zwei_n
check(f"universelle Matrix bei Y = 1 (Zeilen q', (-1)^q', (q' - ni) i^q', (q' + ni)(-i)^q', q' = q - sum eps in [-n, n]; sie "
      f"hängt nur von n ab, da mu(1) = sum eps, mu(+-i) = sum eps +- ni, y = -1 Pol): Rang {universal} für n = 1..4; für n = 2 "
      "also Rang 4 bei 5 Unbekannten für jede Zwei-Paar-Abschneidung", universal[2] == 4 and universal[3] == 4 and universal[4] == 4)
check(f"zwei Paare, alle {zwei_n} Abschneidungen 1 <= m_1 < m_2 <= 31 ungerade, exakt: Rang 4 bei Y = 1 und Y = 3/7, das "
      "Original P_q(Y) liegt im Kern, die Randblöcke Q_- = extrem haben Dimension 1. Also ist jede Lösung P'_q = R(Y) P_q mit "
      "einem R aus Q(Y), der Randblock (ein Eigenwert) erzwingt R = Y^k: blockisospektral bis auf die globale Verschiebung, "
      "für alle Exponenten (z87: vier Abschneidungen im Fenster -8..24)", zwei_ok)
befund("Zwei Paare sind allgemein starr (Satz, kein Fenster): Restklassenrang 4 universell, Kern über Q(Y) vom Original "
       "aufgespannt, Randblock der Dimension 1, dann Max-Entropie. Die Viertelverfeinerung (vier Restklassen) ist dabei "
       "wesentlich", zwei_ok and universal[2] == 4)

# --- 1.3 drei und vier Paare
say("\n   1.3 Drei und vier Paare")
kerne = {}
for ms in ((1, 3, 5), (1, 3, 7), (1, 5, 7), (3, 5, 7), (1, 3, 5, 7)):
    r_, nq, kern_, _ = restklassen(ms, Fraction(3, 7))
    kerne[str(ms)] = {"Rang": r_, "Unbekannte": nq, "Kern": nq - r_, "Original_im_Kern": bool(kern_)}
ZAHLEN["restklassen_kern"] = kerne
check(f"für n >= 3 reicht die lineare Algebra nicht: Rang 4 bei 2n + 1 Unbekannten, Kern der Dimension 2n - 3 "
      f"({', '.join(f'{k}: {v['Kern']}' for k, v in kerne.items())}); Positivität und Ganzzahligkeit müssen tragen",
      all(v["Kern"] == v["Unbekannte"] - 4 and v["Original_im_Kern"] for v in kerne.values())
      and kerne["(1, 3, 5, 7)"]["Kern"] == 5)


def milp_aufbau(ms, lo, hi, schritt=4, normierung="fenster"):
    """Q_- -Randproblem (Relaxation des vollen Problems): w_(q,e) >= 0 ganz, e = q mod schritt, Blocksummen, Identität
    D M' = N Z'. Normierung der globalen Verschiebung über das erste Moment."""
    _, rand, Nc, Dc, _, _ = paarbloecke(ms)
    keys = sorted(rand)
    var = [(q, e) for q in keys for e in range(lo, hi + 1) if (e - q) % schritt == 0]
    idx = {v: i for i, v in enumerate(var)}
    G = max(len(Nc), len(Dc)) - 1
    Z0 = Counter(e for v in rand.values() for e in v)
    F0 = sum(e*c for e, c in Z0.items())
    nz = 4**len(ms)
    A = lil_matrix((len(keys) + (hi + G - lo + 1) + 1, len(var)))
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
    for (q, e), i in idx.items():
        A[r, i] = e
    if normierung == "fenster":
        lb.append(F0 - 2*nz)
        ub.append(F0 + 2*nz - 1)
    else:
        lb.append(F0)
        ub.append(F0)
    r += 1
    return csr_matrix(A[:r]), lb, ub, var, Z0, G, rand


def l1_abweichung(ms, lo, hi, tl=240):
    """max sum_e |Z'_e - Z_e| über alle normierten Lösungen (Binärhilfsvariablen); 0 heißt global isospektral."""
    A, lb, ub, var, Z0, G, _ = milp_aufbau(ms, lo, hi)
    nv, lev = len(var), list(range(lo, hi + 1))
    L = len(lev)
    U = 4**len(ms)
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
    t0 = time.time()
    res = milp(c, constraints=LinearConstraint(AA, lb + lb2, ub + ub2), integrality=integ, bounds=Bounds(np.zeros(len(c)), hib),
               options={"time_limit": tl})
    return (None if res.x is None else -res.fun), res.status, time.time() - t0, nv, G


fenster3 = {}
for ms in ((1, 3, 5), (1, 3, 7), (1, 5, 7), (3, 5, 7), (1, 3, 9), (1, 5, 9)):
    lo, hi = -24, 2*sum(ms) + 24
    val, st, dt, nv, G = l1_abweichung(ms, lo, hi)
    fenster3[str(ms)] = {"Fenster": [lo, hi], "Variablen": nv, "G": G, "max_L1": val, "status": st, "s": round(dt, 1)}
    say(f"      Paare {ms}: Fenster [{lo}, {hi}], {nv} ganzzahlige Variablen, G = {G}, max L1(Z' - Z) = {val:.1e}, "
        f"Status {st}, {dt:.1f} s")
ZAHLEN["drei_paare_milp"] = fenster3
check("drei Paare, sechs Abschneidungen (1,3,5), (1,3,7), (1,5,7), (3,5,7), (1,3,9), (1,5,9), Exponentenfenster [-24, 2 sum m "
      "+ 24], vollständig (HiGHS, Optimalität bewiesen; Q_- -Randproblem, eine Relaxation des vollen Problems): jede Lösung der "
      "Ladungsidentitäten hat Z' = y^(4k) Z, max L1-Abweichung 0. Das geht über z87 hinaus (dort bis vier Eigenwertzüge)",
      all(v["status"] == 0 and v["max_L1"] is not None and v["max_L1"] < 0.5 for v in fenster3.values()))


DEF = {}


def defizit_niveau(e):
    A, lb, ub, var, Z0, tl = DEF["A"], DEF["lb"], DEF["ub"], DEF["var"], DEF["Z0"], DEF["tl"]
    nv = len(var)
    row = csr_matrix(np.array([[1.0 if ee == e else 0.0 for (q, ee) in var]]))
    res = milp(np.zeros(nv), constraints=LinearConstraint(vstack([A, row]), lb + [-np.inf], ub + [Z0[e] - 1]),
               integrality=np.ones(nv), bounds=Bounds(np.zeros(nv), np.full(nv, np.inf)), options={"time_limit": tl})
    return e, int(res.status)


def defizit_suche(ms, lo, hi, tl=120):
    """Für jedes besetzte Niveau e: gibt es eine normierte Lösung mit Z'_e <= Z_e - 1? Alle unlösbar = global isospektral."""
    A, lb, ub, var, Z0, G, _ = milp_aufbau(ms, lo, hi)
    DEF.update(A=A, lb=lb, ub=ub, var=var, Z0=Z0, tl=tl)
    t0 = time.time()
    st = dict(defizit_niveau(e) for e in sorted(Z0))
    return st, time.time() - t0, len(var), G


ms4f, lo4, hi4 = (1, 3, 5, 7), -4, 36
st4, dt4, nv4, G4 = defizit_suche(ms4f, lo4, hi4)
ZAHLEN["vier_paare_defizit"] = {"Fenster": [lo4, hi4], "Variablen": nv4, "G": G4, "status_je_Niveau": st4, "s": round(dt4, 1)}
say(f"      vier Paare {ms4f}: Fenster [{lo4}, {hi4}], {nv4} Variablen, G = {G4}, Defizitprobleme je Niveau: "
    f"{dict(Counter(st4.values()))} (2 = unlösbar), {dt4:.0f} s")
check(f"vier Paare (1, 3, 5, 7), Fenster [{lo4}, {hi4}]: für jedes der {len(st4)} besetzten Niveaus ist 'Z'_e <= Z_e - 1' "
      "unlösbar (HiGHS); da sum Z' = sum Z, hat jede Lösung Z' = Z bis auf die globale Verschiebung",
      all(v == 2 for v in st4.values()))


def teilloesung(ms):
    """Echte Teilmenge des Originalspektrums mit denselben Ladungsidentitäten? (Sie gäbe sofort ein Gegenbeispiel: um y^4
    verschieben.)"""
    _, rand, Nc, Dc, _, _ = paarbloecke(ms)
    w0 = Counter((q, E) for q, Es in rand.items() for E in Es)
    var = sorted(w0)
    G = max(len(Nc), len(Dc)) - 1
    los = [e for (q, e) in var]
    A = lil_matrix((max(los) + G + 2, len(var)))
    for e2 in range(0, max(los) + G + 1):
        for i, (q, e) in enumerate(var):
            d, val = e2 - e, 0
            if 0 <= d < len(Dc):
                val += q*Dc[d]
            if 0 <= d < len(Nc):
                val -= Nc[d]
            if val:
                A[e2, i] = val
    A[max(los) + G + 1, :] = 1
    lb = [0]*(max(los) + G + 1) + [1]
    ub = [0]*(max(los) + G + 1) + [4**len(ms) - 1]
    res = milp(np.zeros(len(var)), constraints=LinearConstraint(csr_matrix(A), lb, ub), integrality=np.ones(len(var)),
               bounds=Bounds(np.zeros(len(var)), np.array([w0[v] for v in var], dtype=float)), options={"time_limit": 120})
    return int(res.status)


teil = {str(ms): teilloesung(ms) for ms in ((1, 3), (1, 3, 5), (1, 3, 7), (1, 3, 5, 7))}
ZAHLEN["echte_teilloesungen_status"] = teil
check(f"notwendige Bedingung (zum Lückenlemma): keine echte Teilmenge des Originalspektrums erfüllt die Ladungsidentitäten "
      f"(Status {teil}, 2 = unlösbar); eine solche Teilmenge wäre um y^4 verschiebbar und ein Gegenbeispiel",
      all(v == 2 for v in teil.values()))


def moment2(ms, lo, hi, schritt, tl=60):
    A, lb, ub, var, Z0, G, _ = milp_aufbau(ms, lo, hi, schritt=schritt, normierung="fest")
    c = -np.array([e*e for (q, e) in var], dtype=float)
    res = milp(c, constraints=LinearConstraint(A, lb, ub), integrality=np.ones(len(var)),
               bounds=Bounds(np.zeros(len(var)), np.full(len(var), np.inf)), options={"time_limit": tl})
    return (None if res.x is None else -res.fun), sum(e*e*c_ for e, c_ in Z0.items())


m2_uhr, m2_0 = moment2((1, 3), -16, 24, 4)
m2_ohne, _ = moment2((1, 3), -16, 24, 1)
ZAHLEN["kontrolle_zweites_moment_(1,3)"] = {"Original": m2_0, "mit Uhr (mod 4)": m2_uhr, "ohne Uhr (nur Gitter Z)": m2_ohne}
check(f"Kontrolle (kann scheitern): ohne Viertelverfeinerung (Eigenwerte nur in y^Z, also bloße Periodizität Delta^(-2i) ~ 1) "
      f"haben die Ladungsidentitäten schon für zwei Paare nicht isospektrale Lösungen (größtes zweites Moment bei festem "
      f"ersten {m2_ohne:.0f} gegen {m2_0}); mit der Uhr R bleibt es {m2_uhr:.0f}: die Starrheit trägt die Uhr, nicht die "
      "Periodizität", m2_ohne is not None and m2_ohne > m2_0 + 0.5 and m2_uhr is not None and abs(m2_uhr - m2_0) < 0.5)

# --- 1.4 Außenbehauptung: gleiche Saaten, diagonale Deformation
say("\n   1.4 Außenbehauptung (elfte Analyse): N gleiche Saaten rho_1, diagonaler Zustand mit Gewichten y^(|n| + 4 f(n))")


def diag_suche(eps, J0=2, schritt=4, tl=120, Lambda=None):
    """Diagonale Zustände y^(E(b) + schritt f(b)), |f| <= J0, f(0) = 0; Randverteilung jedes Modus exakt (Multimengen
    {e(b): b_k = 1} = {e(b) + eps_k: b_k = 0}, Niveau für Niveau); Lambda: nur Niveaus < Lambda (Genauigkeit ~ y^Lambda)."""
    K = len(eps)
    confs = list(itertools.product((0, 1), repeat=K))
    E = [sum(e*b for e, b in zip(eps, c)) for c in confs]
    shifts = list(range(-J0, J0 + 1))
    var = [(ci, j) for ci in range(len(confs)) for j in shifts]
    idx = {v: i for i, v in enumerate(var)}
    Lmin, Lmax = min(E) - schritt*J0, max(E) + schritt*J0 + max(eps)
    nlev = Lmax - Lmin + 1
    A = lil_matrix((len(confs) + K*nlev + 2, len(var)))
    lb, ub, r = [], [], 0
    for ci in range(len(confs)):
        for j in shifts:
            A[r, idx[(ci, j)]] = 1
        lb.append(1)
        ub.append(1)
        r += 1
    for k in range(K):
        for (ci, j), i in idx.items():
            if confs[ci][k] == 1:
                A[r + E[ci] + schritt*j - Lmin, i] += 1
            else:
                A[r + E[ci] + schritt*j + eps[k] - Lmin, i] -= 1
        for L in range(Lmin, Lmax + 1):
            frei = Lambda is not None and L >= Lambda
            lb.append(-np.inf if frei else 0)
            ub.append(np.inf if frei else 0)
        r += nlev
    A[r, idx[(0, 0)]] = 1
    lb.append(1)
    ub.append(1)
    r += 1
    for ci in range(len(confs)):
        A[r, idx[(ci, 0)]] = 1
    lb.append(0)
    ub.append(len(confs) - 1)
    r += 1
    res = milp(np.array([abs(j) for (ci, j) in var], dtype=float), constraints=LinearConstraint(csr_matrix(A[:r]), lb, ub),
               integrality=np.ones(len(var)), bounds=Bounds(np.zeros(len(var)), np.ones(len(var))), options={"time_limit": tl})
    sol = None
    if res.x is not None:
        x = np.round(res.x).astype(int)
        sol = {confs[ci]: j for (ci, j), xi in zip(var, x) if xi and j != 0}
    return int(res.status), sol


saat = {N: diag_suche([1]*N, J0=8)[0] for N in range(2, 6)}
st6, sol6 = diag_suche([1]*6, J0=5, tl=200)
saat_ohne = {N: diag_suche([1]*N, J0=3, schritt=1)[0] for N in (3, 4, 5)}
confs6 = list(itertools.product((0, 1), repeat=6))
f6 = {c_: (sol6 or {}).get(c_, 0) for c_ in confs6}


def pruefe_diag(f, eps, faktor=1):
    """Exakt: Multimengenidentität je Modus; 50 Stellen: Randverteilungen, Uhr (Exponent = E mod 4), Paarkorrelationen,
    Entropie gegen das Produkt."""
    K = len(eps)
    cf = list(itertools.product((0, 1), repeat=K))
    e = {c_: faktor*(sum(a*b for a, b in zip(eps, c_)) + 4*f[c_]) for c_ in cf}
    eps_s = [faktor*x for x in eps]
    multi = all(Counter(e[c_] for c_ in cf if c_[k]) == Counter(e[c_] + eps_s[k] for c_ in cf if not c_[k]) for k in range(K))
    uhr = all((e[c_] - sum(a*b for a, b in zip(eps_s, c_))) % 4 == 0 for c_ in cf)
    yM = mp.exp(-mp.pi)
    gew = {c_: yM**e[c_] for c_ in cf}
    Zg = mp.fsum(gew.values())
    p = {c_: gew[c_]/Zg for c_ in cf}
    rand = max(abs(mp.fsum(v for c_, v in p.items() if c_[k]) - yM**eps_s[k]/(1 + yM**eps_s[k])) for k in range(K))
    paark = []
    for i_, j_ in itertools.combinations(range(K), 2):
        ni = mp.fsum(v for c_, v in p.items() if c_[i_])
        nj = mp.fsum(v for c_, v in p.items() if c_[j_])
        paark.append(mp.fsum(v for c_, v in p.items() if c_[i_] and c_[j_]) - ni*nj)
    S1 = -mp.fsum(v*mp.log(v) for v in p.values())
    S0 = -mp.fsum(c*mp.log(c) + (1 - c)*mp.log(1 - c) for c in (yM**x/(1 + yM**x) for x in eps_s))
    return {"Multimengen_exakt": multi, "Uhr": uhr, "max_Randfehler_50_Stellen": float(rand),
            "Paarkorrelationen_ungleich_0": sum(abs(v) > mp.mpf(10)**-45 for v in paark), "Paare": len(paark),
            "min_|Paarkorrelation|": float(min(abs(v) for v in paark)), "max_|Paarkorrelation|": float(max(abs(v) for v in paark)),
            "S": float(S1), "S_Produkt": float(S0), "f_Bereich": [min(f.values()), max(f.values())],
            "verschobene_Konfigurationen": sum(v != 0 for v in f.values())}


bsp6 = pruefe_diag(f6, [1]*6) if sol6 else {}
bsp6_3 = pruefe_diag(f6, [1]*6, faktor=3) if sol6 else {}            # Exponenten 3(|n| + 4f) = 3|n| + 4(3f)
ZAHLEN["aussenbehauptung_gleiche_saaten"] = {"N = 2..5, |f| <= 8 (Status)": saat, "N = 6, |f| <= 5 (Status)": st6,
                                             "Beispiel N = 6 (drei Paare m = 1)": bsp6, "skaliert: drei Paare m = 3": bsp6_3,
                                             "ohne Uhr, Schritt 1 (Status)": saat_ohne,
                                             "f(n) für N = 6": {"".join(map(str, c_)): v for c_, v in f6.items() if v}}
check(f"die Außenbehauptung stimmt (exakt): N = 2..5 gleiche Saaten sind diagonal starr (|f| <= 8, Status {saat}, 2 = "
      f"unlösbar); bei N = 6 gibt es einen nichtgaußschen diagonalen Zustand y^(|n| + 4 f(n)) mit f in {bsp6.get('f_Bereich')}, "
      f"{bsp6.get('verschobene_Konfigurationen')} verschobenen Konfigurationen: Randverteilungen exakt (Multimengen je Modus "
      f"{bsp6.get('Multimengen_exakt')}, Fehler {bsp6.get('max_Randfehler_50_Stellen'):.0e} bei 50 Stellen), Uhr exakt "
      f"({bsp6.get('Uhr')}), {bsp6.get('Paarkorrelationen_ungleich_0')} von 15 Paarkorrelationen ungleich 0 (|C| zwischen "
      f"{bsp6.get('min_|Paarkorrelation|'):.1e} und {bsp6.get('max_|Paarkorrelation|'):.1e}), S = {bsp6.get('S'):.4f} gegen "
      f"S_Produkt = {bsp6.get('S_Produkt'):.4f} (nicht isospektral, Subadditivität greift nicht)",
      all(v == 2 for v in saat.values()) and st6 == 0 and bsp6.get("Multimengen_exakt") and bsp6.get("Uhr")
      and bsp6.get("max_Randfehler_50_Stellen", 1) < 1e-40 and bsp6.get("Paarkorrelationen_ungleich_0") == 15
      and bsp6.get("S", 9) < bsp6.get("S_Produkt", 0) - 0.1)

entartet = {str(ms): diag_suche([m for m in ms for _ in range(2)], J0=6)[0] for ms in ((1, 1), (1, 1, 3), (1, 1, 5), (1, 3, 3))}
orig_diag = {str(ms): diag_suche([m for m in ms for _ in range(2)], J0=6)[0] for ms in ((1, 3, 5), (1, 3, 7), (1, 5, 9))}
orig_diag["(1, 3, 5, 7) |f|<=2"] = diag_suche([1, 1, 3, 3, 5, 5, 7, 7], J0=2, tl=200)[0]
ZAHLEN["entartete_paare_diagonal_status"] = entartet
ZAHLEN["originalpaare_diagonal_status"] = orig_diag
check(f"welche Prämisse rettet: drei Paare mit gleichem m sind dasselbe Gegenbeispiel (m = 1: die sechs Anregungsmoden; m = 3: "
      f"f -> 3f, Multimengen {bsp6_3.get('Multimengen_exakt')}, Uhr {bsp6_3.get('Uhr')}); Ladungssuperauswahl erfüllt jeder "
      f"diagonale Zustand von selbst. Teilweise entartet (1,1), (1,1,3), (1,1,5), (1,3,3) und die Paarstruktur des Originals "
      f"mit verschiedenen m, (1,3,5), (1,3,7), (1,5,9) mit |f| <= 6 und (1,3,5,7) mit |f| <= 2, haben keinen solchen Zustand "
      f"(Status {entartet}, {orig_diag}; 2 = unlösbar); ohne Uhr gibt es ihn schon bei N = 3 Saaten ({saat_ohne}, 0 = gefunden)",
      bsp6_3.get("Multimengen_exakt") and bsp6_3.get("Uhr") and all(v == 2 for v in entartet.values())
      and all(v == 2 for v in orig_diag.values()) and all(v == 0 for v in saat_ohne.values()))
befund("Die Außenbehauptung ist richtig, widerspricht z87 aber nicht: Uhrenidentität + Ladungssuperauswahl + Einteilchendaten "
       "erzwingen den Gaußzustand NICHT allgemein (Gegenbeispiel: sechs Anregungsmoden gleicher Energie, also drei Paare mit "
       "gleichem m). Die Superauswahl rettet nichts; was in allen gerechneten Fällen rettet, ist die nicht entartete "
       "Paarstruktur des Originals: jedes ungerade m genau einmal, Teilchen und Loch in verschiedenen Quartersektoren (Rahmen "
       "der Viertelzählung). Diese Nichtentartung muss als Teil des Rahmens ausdrücklich genannt werden, und ein allgemeiner "
       "Beweis von S(n) muss sie benutzen", bool(sol6) and all(v == 2 for v in orig_diag.values()))
befund("Allgemein (beliebig viele Paare mit verschiedenen m) bleibt S(n) offen: bewiesen für n = 2, vollständig im Fenster "
       "für sechs Drei-Paar-Abschneidungen (Fenster +-24) und (1, 3, 5, 7) (enges Fenster [-4, 36]), diagonal kollisionsfrei "
       "bis vier Paare; der lineare Weg endet bei "
       "Kern 2n - 3, das Gegenbeispiel mit entarteten Paaren zeigt, dass die Positivität allein nicht reicht", True)
say(f"   (Laufzeit Teil 1: {time.time() - T_START:.0f} s)")

# =================================================================== Teil 2
say("\nTeil 2: Folgt P_clock aus P_sym über den Tausch A <-> B?")
sx, sy_, sz = np.array([[0, 1], [1, 0]], dtype=complex), np.array([[0, -1j], [1j, 0]]), np.diag([1.0 + 0j, -1.0])
I2 = np.eye(2)


def implementierbar(bild_A, bild_B, A, B):
    zeilen = []
    for i in range(2):
        for j in range(2):
            U = np.zeros((2, 2), dtype=complex)
            U[i, j] = 1
            zeilen.append(np.concatenate([(U @ A - bild_A @ U).ravel(), (U @ B - bild_B @ U).ravel()]))
    sv = np.linalg.svd(np.array(zeilen).T, compute_uv=False)
    return bool(np.sum(sv < 1e-10) > 0)


winkel = np.unique(np.round(np.append(np.linspace(0.05, np.pi - 0.05, 61), np.pi/2), 12))
impl = {"Tausch A<->B": [], "Vorzeichen A->-A": [], "Viertel A->B->-A": []}
for th in winkel:
    A_, B_ = sx, np.cos(th)*sx + np.sin(th)*sy_
    for name, (bA, bB) in (("Tausch A<->B", (B_, A_)), ("Vorzeichen A->-A", (-A_, B_)), ("Viertel A->B->-A", (B_, -A_))):
        if implementierbar(bA, bB, A_, B_):
            impl[name].append(round(float(th/np.pi), 4))
ZAHLEN["algebra_implementierbar_bei_theta/pi"] = {k: (v if len(v) < 5 else f"{len(v)} von {len(winkel)}") for k, v in impl.items()}
check(f"Algebra-Seite ({len(winkel)} Winkel theta zwischen A und B in der irreduziblen 2D-Darstellung): der Tausch A <-> B ist "
      f"immer innerlich ({len(impl['Tausch A<->B'])} von {len(winkel)}, U = (A + B)/|A + B|), die Vorzeichenspiegelung und die "
      f"Vierteldrehung A -> B -> -A nur bei theta = pi/2 ({impl['Vorzeichen A->-A']}, {impl['Viertel A->B->-A']}): P_sym ist "
      "dieselbe Aussage wie 'Vierteldrehung (komplexe Multiplikation mit i) innerlich'; der Tausch selbst trägt keine Information",
      len(impl["Tausch A<->B"]) == len(winkel) and impl["Vorzeichen A->-A"] == [0.5] and impl["Viertel A->B->-A"] == [0.5])

A, B = sx, sy_
Gam = -1j*A @ B
Ux = (A + B)/np.sqrt(2)
Qv = sla.expm(-0.25j*np.pi*Gam)


def rho_g(g):
    """(1 + tanh(g/2) Gamma)/2 = diag(1, e^-g)/(1 + e^-g) in der Gamma-Basis (Gamma = -iAB = diag(1, -1)), ohne Auslöschung."""
    assert np.allclose(Gam, np.diag([1, -1]))
    return np.diag([1.0, np.exp(-g)])/(1 + np.exp(-g)) + 0j


def logd(rho):
    assert np.allclose(rho, np.diag(np.diag(rho)))
    return np.diag(np.log(np.diag(rho).real)).astype(complex)


def sigma(rho, X, t):
    L = logd(rho)
    return sla.expm(1j*t*L) @ X @ sla.expm(-1j*t*L)


gs = np.array([-3, -1, -0.5, 0.25, 0.5, 1, 2, 3, 5, 9])*np.pi
cov_ok = all(np.allclose(Ux @ rho_g(g) @ Ux.conj().T, rho_g(-g)) for g in gs)
inv_b0 = np.allclose(Ux @ (I2/2) @ Ux.conj().T, I2/2) and all(not np.allclose(Ux @ rho_g(g) @ Ux.conj().T, rho_g(g)) for g in gs)
rueck_ok = all(np.allclose(sigma(rho_g(-g), A, 0.3), sigma(rho_g(g), A, -0.3)) for g in gs)
viertel_inv = all(np.allclose(Qv @ rho_g(g) @ Qv.conj().T, rho_g(g)) for g in gs)
tQ_ok = all(np.allclose(sigma(rho_g(g), A, -np.pi/(2*g)), B) and np.allclose(sigma(rho_g(g), B, -np.pi/(2*g)), -A) for g in gs)
check("Gesetz-Seite auf dem Qubit (gerade Zustände rho_g = (1 + tanh(g/2) Gamma)/2, zehn g): der Tausch bildet rho_g auf "
      "rho_-g ab (Zeitumkehr, sigma^-g_t = sigma^g_-t), invariant nur der Spurzustand g = 0 ohne Zeit (z80); die Vierteldrehung "
      "lässt jedes rho_g invariant und liegt für jedes g != 0 in der Modulgruppe, bei t_Q = -pi/(2g): Tausch-Kovarianz, "
      "Tausch-Invarianz und Vierteldrehung legen g nicht fest", cov_ok and inv_b0 and rueck_ok and viertel_inv and tQ_ok)

tab2 = {}
for g in gs:
    rho = rho_g(g)
    w = np.sort(-np.log(np.diag(rho).real))
    gap = float(w.max() - w.min())                                       # Spektrum von log Delta auf ungeraden Elementen: +-g
    Ta = np.pi/gap
    ts = np.linspace(1e-3, Ta*(1 - 1e-3), 400)
    kleiner = any(np.allclose(sigma(rho, A, t), -A, atol=1e-8) for t in ts)
    anti = np.allclose(sigma(rho, A, Ta), -A) and np.allclose(sigma(rho, B, Ta), -B)
    X_, Y_ = A + 0.3*B, B - 0.2j*A
    L_ = logd(rho)
    kms = np.isclose(np.trace(rho @ X_ @ (sla.expm(L_) @ Y_ @ sla.expm(-L_))), np.trace(rho @ Y_ @ X_))
    tau = 1j*1.0/Ta
    tau_s = -1/tau
    tn = abs(tau.imag - tau_s.imag) < 1e-12 and abs(tau.real - tau_s.real) < 1e-12
    pc = np.allclose(sigma(rho, A, -0.5), B) and np.allclose(sigma(rho, B, -0.5), -A)
    tab2[round(float(g/np.pi), 3)] = {"T_a": round(Ta, 6), "tau_Gesetz": f"{tau.imag:.4g} i", "zeitneutral": bool(tn),
                                       "P_clock": bool(pc), "KMS_beta1": bool(kms), "antiperiodisch": bool(anti and not kleiner)}
ZAHLEN["qubit_modultorus"] = tab2
for k_, v_ in tab2.items():
    say(f"      g = {k_:6} pi: T_a = {v_['T_a']:.4f}, tau_Gesetz = {v_['tau_Gesetz']:>8}, zeitneutral {v_['zeitneutral']!s:5}, "
        f"P_clock {v_['P_clock']}")
tn_set = sorted(k for k, v in tab2.items() if v["zeitneutral"])
pc_set = sorted(k for k, v in tab2.items() if v["P_clock"])
check(f"komplexer Modulzeit-Torus des Gesetzes: reelle Antiperiode T_a = pi/|g| (kleinste, sigma_Ta = (-1)^F auf A, B), "
      f"KMS-Antiperiode 1 (KMS bei beta = 1 geprüft), Spinstruktur (A, A), tau_Gesetz = i|g|/pi. Zeitneutralität im Sinn von n1 "
      f"(beide Rekonstruktionen tau und -1/tau geben dieselben Daten) gilt genau bei g/pi in {tn_set}; P_clock (sigma_-1/2 = "
      f"Vierteldrehung) bei {pc_set}; Schnitt {sorted(set(tn_set) & set(pc_set))}; g = -pi ist das Tauschbild (Partner, "
      "Orientierung)", tn_set == [-1.0, 1.0] and pc_set == [-3.0, 1.0, 5.0, 9.0] and sorted(set(tn_set) & set(pc_set)) == [1.0]
      and all(v["KMS_beta1"] and v["antiperiodisch"] for v in tab2.values()))

# Original: KMS-Familie K_beta = beta sum x n_x
xs_o = np.array(sorted([m/4 for m in range(1, 16, 2)] + [-m/4 for m in range(1, 16, 2)]))


def antiperiode(beta):
    """Kleinste T > 0 mit e^(i beta x T) = -1 auf allen Moden: Kandidat aus dem kleinsten Modus |x| = 1/4, dann geprüft
    (alle Moden antiperiodisch, keine kleinere Antiperiode auf einem Raster)."""
    Ta = 4*np.pi/abs(beta)
    alle = np.abs(np.exp(1j*beta*xs_o*Ta) + 1).max() < 1e-9
    ts = np.linspace(1e-3, Ta*(1 - 1e-3), 4000)
    kleiner = bool((np.abs(np.exp(1j*np.outer(ts, beta*xs_o)) + 1).max(axis=1) < 1e-6).any())
    return Ta if alle and not kleiner else None


def uhr_beta(beta):
    """rho_beta = e^(-K_beta)/Z, K_beta = beta sum x (n_x - Meer) = beta E/4 auf zwei Paaren; rho^(-i/2) = e^(iK/2) bis auf
    eine Phase. Rest gegen R und gegen R^-1."""
    ph = np.exp(0.5j*beta*E4/4.0)
    r = np.diag(R4)
    a, b = ph/r, ph*r
    return float(np.abs(a - a[0]).max()), float(np.abs(b - b[0]).max())


tab_o = {}
for bf in (2, 4, 8, 12, 20, 4.04, -4):
    beta = bf*np.pi
    Ta = antiperiode(beta)
    tau = 1j*1.0/Ta
    tn = abs(tau.imag - (-1/tau).imag) < 1e-4
    rest, rest_inv = uhr_beta(beta)
    tab_o[bf] = {"T_a": round(Ta, 5), "tau_Gesetz": round(tau.imag, 5), "zeitneutral": bool(tn), "P_clock": bool(rest < 1e-8),
                 "Delta^(-i/2)=R^-1": bool(rest_inv < 1e-8)}
ZAHLEN["original_beta_scan"] = {f"beta = {k} pi": v for k, v in tab_o.items()}
for k_, v_ in tab_o.items():
    say(f"      beta = {k_:5} pi: T_a = {v_['T_a']:.5f}, tau_Gesetz = {v_['tau_Gesetz']:.4f} i, zeitneutral {v_['zeitneutral']!s:5}, "
        f"P_clock {v_['P_clock']!s:5}, Delta^(-i/2) = R^-1 {v_['Delta^(-i/2)=R^-1']}")
tn_o = sorted(k for k, v in tab_o.items() if v["zeitneutral"])
pc_o = sorted(k for k, v in tab_o.items() if v["P_clock"])
check(f"Original (16 Quartermoden +-m/4, m <= 15; Uhr auf dem Fockraum zweier Paare): tau_Gesetz = i beta/(4 pi); zeitneutral "
      f"genau bei beta/pi in {tn_o}, P_clock bei {pc_o}, Schnitt {sorted(set(tn_o) & set(pc_o))}; beta = -4 pi ist das Bild "
      "unter der Modenspiegelung x -> -x (Partner, R -> R^-1). Die Zeitneutralität wählt beta = 4 pi allein, ohne den "
      "Uhrenquotienten (z74), den z87 gegen beta = 20 pi brauchte",
      tn_o == [-4, 4] and pc_o == [4, 20] and sorted(set(tn_o) & set(pc_o)) == [4] and tab_o[-4]["Delta^(-i/2)=R^-1"])


def z_ns(tau):
    q12 = mp.exp(1j*mp.pi*tau)
    eta = mp.exp(2j*mp.pi*tau/24)*mp.qp(mp.exp(2j*mp.pi*tau))
    return abs(mp.jtheta(3, 0, q12)/eta)**2


s_inv = {float(yv): float(abs(z_ns(1j*mp.mpf(yv)) - z_ns(1j/mp.mpf(yv)))/z_ns(1j*mp.mpf(yv))) for yv in (0.5, 0.8, 1.25, 2.0)}
ZAHLEN["Z_NS_S_invarianz_rel"] = s_inv
check(f"S-Transformation: die NS-Zustandssumme |theta_3(tau)/eta(tau)|^2 ist für jedes tau = iy S-invariant (rel. Abweichung "
      f"<= {max(s_inv.values()):.0e} bei y = 0,5, 0,8, 1,25, 2): die Invarianz der Funktion wählt nichts aus, nur die "
      "Selbstdualität des einzelnen Gesetzestorus (tau = -1/tau) wählt tau = i", max(s_inv.values()) < 1e-40)


def j_marken(th):
    z = [1, -1, np.exp(1j*th), -np.exp(1j*th)]
    lam = (z[0] - z[2])*(z[1] - z[3])/((z[0] - z[3])*(z[1] - z[2]))
    return float((256*(1 - lam + lam**2)**3/(lam**2*(1 - lam)**2)).real)


jm = {"pi/2": j_marken(np.pi/2), "pi/3": j_marken(np.pi/3)}
gesetz_2pi = tab2[2.0]["tau_Gesetz"]
unab_1 = abs(jm["pi/2"] - 1728) < 1e-8 and not tab2[2.0]["zeitneutral"]
unab_2 = abs(jm["pi/3"] - 1728) > 1 and not implementierbar(-sx, np.cos(np.pi/3)*sx + np.sin(np.pi/3)*sy_, sx,
                                                               np.cos(np.pi/3)*sx + np.sin(np.pi/3)*sy_)
A3, B3 = sx, np.cos(np.pi/3)*sx + np.sin(np.pi/3)*sy_
G3 = (A3 @ B3 - B3 @ A3)/(2j*np.sin(np.pi/3))
rho3 = sla.expm(-np.pi/2*G3)
rho3 /= np.trace(rho3)
w3 = np.linalg.eigvalsh(-sla.logm(rho3).real)
tau3 = 1j*(w3.max() - w3.min())/np.pi
ZAHLEN["unabhaengigkeit"] = {"j(Marken) theta = pi/2": jm["pi/2"], "j(Marken) theta = pi/3": jm["pi/3"],
                             "tau_Gesetz bei P_sym und g = 2 pi": gesetz_2pi, "tau_Gesetz bei theta = pi/3, Lücke pi": str(tau3)}
check(f"Unabhängigkeit: P_sym ohne Zeitneutralität ist möglich (theta = pi/2, Marken j = {jm['pi/2']:.0f}, aber g = 2 pi gibt "
      f"tau_Gesetz = {gesetz_2pi}); Zeitneutralität ohne P_sym ebenso (theta = pi/3: Lücke pi gibt tau_Gesetz = "
      f"{tau3.imag:.4f} i, aber Marken mit j = {jm['pi/3']:.1f} != 1728, Vorzeichenspiegelung nicht innerlich). Unter P_sym "
      "ist Zeitneutralität dasselbe wie tau_Gesetz = tau_Marken ('Zeit = zweiter Zyklus', z87 Kette Schritt 4)",
      unab_1 and unab_2 and abs(tau3 - 1j) < 1e-9)
befund("P_clock folgt nicht aus P_sym: P_sym ist eine Aussage über die Algebra (Vierteldrehung innerlich, Marken als Quadrat, "
       "tau_Marken = i) und lässt jede Lücke g zu; der Tausch A <-> B ist auf der Algebra immer innerlich und auf dem Gesetz "
       "die Zeitumkehr (rho_g -> rho_-g), legt also nichts fest", True)
befund("P_clock ist gleichwertig zur Zeitneutralität des Gesetzes P_tn (n1 auf dem komplexen Modulzeit-Torus: reelle "
       "Antiperiode = KMS-Antiperiode, tau_Gesetz = i), bis auf die Orientierung, die der Tausch A <-> B als Benennung "
       "(Original/Partner) festlegt. P_tn ist schärfer als P_clock: es enthält die Minimalität (g = pi statt pi(4Z+1), "
       "beta = 4 pi statt 4 pi(4Z+1)), für die z87 den Uhrenquotienten brauchte", True)
say(f"   (Laufzeit bis Teil 2: {time.time() - T_START:.0f} s)")

# =================================================================== Teil 3
say("\nTeil 3: Bilanz für F6")
STATUS = {
    "Satz S(n), Reduktion (jedes n)": ("E", "Gitterspektrum je Ladungsblock, Identitäten über Z[y] (Gelfond), Max-Entropie: "
                                           "S(n) <=> Ladungsidentitäten erzwingen Z' = y^(4k) Z; Lückenlemma"),
    "S(2) zwei Paare": ("E (neu, allgemein)", "alle ungeraden m_1 < m_2, alle Exponenten: Restklassenrang 4, Randblock"),
    "S(3), S(4)": ("E (im Fenster)", "sechs Drei-Paar-Abschneidungen (Fenster +-24) und (1,3,5,7) (eng, [-4, 36]) "
                                      "vollständig; Kern 2n - 3"),
    "S(n) allgemein (verschiedene m)": ("O", "Beweis muss die Nichtentartung benutzen; kein Gegenbeispiel"),
    "Außenbehauptung gleiche Saaten": ("E (bestätigt)", "N = 6 gleiche Saaten = drei Paare mit gleichem m: exakter "
                                                        "nichtgaußscher Zustand mit Uhr; N <= 5 starr; Originalstruktur frei"),
    "P_clock aus P_sym": ("E (negativ, gescopt)", "P_sym lässt jede Lücke zu; Tausch = Zeitumkehr"),
    "P_clock <=> P_tn": ("E", "Qubit und Original; P_tn enthält die Minimalität (Uhrenquotient)"),
    "K19 Existenz der Naht": ("L", "unter P_tn ist die Naht der reelle Zyklus des Modulzeit-Torus, so lang wie der KMS-Zyklus; "
                                   "Existenz nicht hergeleitet"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
PRAEMISSEN = {
    "P_sym": "zwei Binäroperationen, Vorzeichensymmetrie bzw. Vierteldrehung in derselben irreduziblen Darstellung (z80) -> "
             "Cl_2 = CAR, Marken als Quadrat, tau_Marken = i",
    "P_tn (ersetzt P_clock)": "Zeitneutralität des Gesetzes: der komplexe Modulzeit-Torus (reelle Antiperiode x KMS-Antiperiode) "
                              "ist selbstdual unter tau -> -1/tau, tau_Gesetz = i; gleichwertig P_clock + Uhrenquotient; unter "
                              "P_sym gleich 'Zeit = zweiter Zyklus des Markentorus'",
    "Rahmen": "Hardy/Szegő-Naht mit Drehung L_0 und Viertelzählung (jedes ungerade m genau einmal, nicht entartet: für "
              "die Starrheit nötig), Ladungssuperauswahl je Quarter",
    "P_free": "Folge: allgemein für zwei Paare, im Fenster für drei und vier Paare; allgemein offen",
}
ZAHLEN["praemissen"] = PRAEMISSEN
ZAHLEN["praemissen_anzahl"] = {"benannt": 2, "davon neu gefasst": 1, "Rahmen": 1, "bedingt (P_free)": "Folge bis n = 2 "
                                                                                                      "allgemein, n = 3, 4 im Fenster"}
check("F6-Kriterien: Starrheit für die Originalstruktur ohne Gegenbeispiel (Satz für n = 2, Fenster für n = 3, 4), "
      "Außenbehauptung bestätigt, aber nur für entartete Paare; "
      "P_clock nicht aus P_sym, aber gleichwertig zu P_tn, das die Minimalität mitbringt; P_sym und P_tn voneinander "
      "unabhängig: zwei Prämissen in schärferer Form",
      CHECKS[[k for k in CHECKS if k.startswith("zwei Paare, alle")][0]]
      and CHECKS[[k for k in CHECKS if k.startswith("drei Paare")][0]]
      and CHECKS[[k for k in CHECKS if k.startswith("die Außenbehauptung")][0]]
      and CHECKS[[k for k in CHECKS if k.startswith("komplexer Modulzeit")][0]]
      and CHECKS[[k for k in CHECKS if k.startswith("Unabhängigkeit")][0]])

urteil = ("F6 trägt in gescopter Form. (1) Starrheit: Satz S(n) ist für jedes n auf eine Aussage über Blockgewichte reduziert "
          "(Uhrenidentität = Funktionalgleichung Phi(iy, s) = Phi(y, is), Ladungsidentitäten über Z[y] wegen der Transzendenz "
          "von e^pi, dann Max-Entropie: K_0 ist Einteilchenoperator, S(rho'||rho_qf) = S(rho_qf) - S(rho')). Neu bewiesen ist "
          "S(2) für alle ungeraden m_1 < m_2 und alle Exponenten: Die vier Restklassen von y mod 4 geben ein lineares System "
          "über Q(Y), dessen universelle Matrix bei Y = 1 Rang 4 hat; der Kern ist vom Original aufgespannt, der Randblock "
          "(ein Eigenwert) erzwingt R = Y^k. Für n >= 3 hat das System einen Kern der Dimension 2n - 3; ganzzahlig vollständig "
          "im Exponentenfenster sind sechs Drei-Paar-Abschneidungen (+-24) und (1, 3, 5, 7) (eng, [-4, 36]) starr, und keine "
          "echte Teilmenge des "
          "Originalspektrums erfüllt die Identitäten. Ohne Viertelverfeinerung (bloße Periodizität) bricht die Starrheit "
          "schon bei zwei Paaren: Die Uhr R trägt sie. Die Behauptung der elften Außenanalyse stimmt exakt: Für N = 6 gleiche "
          "Saaten (drei Paare mit gleichem m) gibt es einen diagonalen nichtgaußschen Zustand y^(|n| + 4 f(n)) mit exakten "
          "Randverteilungen, exakter Uhr und allen 15 Paarkorrelationen ungleich 0; N <= 5 ist starr. Das widerspricht z87 "
          "nicht, denn dessen Abschneidungen haben verschiedene m; die Ladungssuperauswahl rettet nichts (diagonale Zustände "
          "erfüllen sie), es rettet die nicht entartete Paarstruktur des Originals (jedes m einmal, Teilchen und Loch in "
          "verschiedenen Quartern): teilweise entartete Strukturen und die Originalstruktur bis vier Paare haben keinen "
          "solchen Zustand. Die Uhrenidentität mit Superauswahl und Einteilchendaten genügt also nicht allgemein; S(n) für "
          "verschiedene m bleibt offen, ein Beweis muss die Nichtentartung benutzen. (2) P_clock folgt nicht aus P_sym: P_sym (Vierteldrehung innerlich, tau_Marken "
          "= i) lässt jede Lücke g zu; der Tausch A <-> B ist auf der Algebra stets innerlich und auf dem Gesetz die "
          "Zeitumkehr rho_g -> rho_-g. Gleichwertig zu P_clock ist die Zeitneutralität des Gesetzes P_tn im Sinn von n1: "
          "Der komplexe Modulzeit-Torus (reelle Antiperiode pi/|g|, KMS-Antiperiode 1) ist genau bei |g| = pi selbstdual, "
          "auf dem Original genau bei beta = 4 pi. P_tn ist schärfer als P_clock, denn es enthält die Minimalität, für die "
          "z87 den Uhrenquotienten brauchte; die Orientierung ist die Benennung Original/Partner. (3) Es bleiben zwei "
          "Ursprungsprämissen, P_sym und P_tn, die beiden Hälften derselben Symmetrie tau -> -1/tau (komplexe Multiplikation "
          "mit i) für Algebra und Gesetz, nachweislich unabhängig, dazu der Rahmen der Hardy-Naht mit Ladungssuperauswahl, "
          "in dem die Nichtentartung der Viertelzählung jetzt ausdrücklich gebraucht wird. P_free ist Folge für zwei Paare "
          "allgemein und für drei und vier Paare im Fenster. K19 wird nicht hergeleitet; "
          "unter P_tn ist die Naht der reelle Zyklus des Modulzeit-Torus, gleich lang wie der KMS-Zyklus.")
say("\nUrteil\n   " + urteil)
say(f"\nLaufzeit {time.time() - T_START:.0f} s")
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z94_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
