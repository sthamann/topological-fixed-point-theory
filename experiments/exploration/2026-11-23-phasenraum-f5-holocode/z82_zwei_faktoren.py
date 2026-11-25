"""Z82: Strategie S1 'Zwei-Faktoren-Architektur' -- innen das Ladungsgitter (Uhr R, Marken, 16), außen die Oszillatoren.

Leitlinie 'Rekonstruktion statt Auswahl'. Firewalls: GRAV.NONCIRCULAR.01 (kein Schluss modular => Bisognano-Wichmann; hier
wird keine Raumzeitsymmetrie aus dem Modulfluss gefolgert), DIMENSION.UPLIFT.FIREWALL.01 (jede 3+1-Aussage ist als Kinematik
oder Lesart typisiert). Gelöste Bausteine (16, E8, Flavour, Massen, alpha) werden nur benutzt. Bestand: z58 (Bosonisierung,
Ladung x Oszillatoren), z76 (K/pi = Q_a - Q_b mod 4), z79 (Quarterpaar = eta-Spinor, Lorentzgruppe auf den Oszillatoren), z81
(R = i^(Q_a - Q_b), R^2 = (-1)^F), Träger e_m = z^m/(z^4 - 1) mit Y = (I - 5R^2)/12 (Atlas: rr-w-fixed-hypercharge,
quadratic-fock-five-carrier), V4-Schnitt (z34, z36), Z4 der 16: z_X = exp(i pi X/2), X = 5 - 2p (a3, s4, z23, z53).

Teil 1  Ist die Uhr R die Z4 der 16? R auf den fünf Trägermoden (Viertelcharaktere, R^2 teilt 3+2, Y); det R = -1 und
        R = (-1_5) g mit g in G_SM (beide 3+2-Schnitte); Spinorhebung auf Lambda^gerade C^5 = 16 gleich z_X mal
        Eichelement, naive Hebung eichtrivial; Ladungen der 16 und mod-16-Bilanz; Gegenproben (N = 4, verschobene
        Charaktere, Abtastung über N). Abbildung des Originals auf den Träger über die Hardy-Marken k = 4x
        (Residuenklassen der Laurentreihen): äquivariant, trifft nur die 2. Wo die 3 (R^2 = +1) herkäme: Ad_R auf der
        CAR-Algebra, Determinantenschranke, niedrigste Sektoren, Statistik der 16 bei exakten Trägercharakteren.
Teil 2  Spin-Statistik: geometrische Spinorhebung des eta-Dubletts (exakt); auf allen Ladungssektoren |Q_a|, |Q_b| <= 10 mit
        Oszillatorinhalt ist das Produkt der Spinoreigenwerte gleich R und sein Quadrat (2 pi-Drehung) gleich (-1)^F;
        alle 16 Zuordnungen der Quartergrade durchgemustert. Innere U(1): J_z = (Q_a - Q_b)/2. Implementierbarkeit der
        vollen eta-SU(2) im KMS-Zustand. Gegenprobe Außenbild: die 2 pi-Drehung der Lorentzgruppe aus z79 auf Strom-
        quanten und auf Fermimoden psi*_x|0>.
Teil 3  Universalität: Stromalgebra, Oszillatorzustand und T_osc sind in jedem Ladungssektor gleich (fermionisch in der
        Maya-Basis nachgerechnet: Matrixelemente, CCR, Zwei- und verbundene Vierpunktantworten bei beta = 4 pi und
        beta = 2, Kreuzkovarianz innen-außen); Sugawara L_n = L_n^osc + (Q + alpha - 1/2) J_n; die Tensorantwort
        <L_-n L_n>_Q hängt nur über die innere Energie E(Q) ab, mit gleichem Koeffizienten für beide Quarter.
Teil 4  Bilanz und Urteil für S1 (Erfolg und Abbruch aus docs/research/GESAMTSTAND-UND-STRATEGIEN.md, Abschnitt 7).
"""
import itertools
import json
import math
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.linalg import expm

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
rng = np.random.default_rng(82)
BETA = 4*np.pi


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


# ------------------------------------------------------------------ Hilfen: ein Quarter in der Maya-Basis
def partitionen(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitionen(n - k, k):
            yield (k,) + rest


class Quarter:
    """Fermionischer Fockraum eines Quarters (Moden x = y + alpha, y ganz), Sektoren Q, Oszillatorstufe <= N.

    Zustand (Q, lambda): besetzt sind y_i = lambda_i - i + Q (i = 1, 2, ...), also L_0 = E_alpha(Q) + |lambda|. Bitmaske über
    das Fenster [lo, hi); unterhalb lo ist alles besetzt, oberhalb hi alles leer. Vorzeichen der CAR: Anzahl besetzter
    Moden oberhalb (endlich, da die See nach unten unendlich ist)."""

    def __init__(self, alpha, Qs, N, M):
        self.alpha, self.N = alpha, N
        rand = max(abs(q) for q in Qs) + N + M + 2
        self.lo, self.hi = -rand, rand
        self.sek = {}
        for Q in Qs:
            zust, stufen = [], []
            for n in range(N + 1):
                for lam in partitionen(n):
                    zust.append(self.maya(Q, lam))
                    stufen.append(n)
            self.sek[Q] = (zust, {s: j for j, s in enumerate(zust)}, np.array(stufen))

    def maya(self, Q, lam):
        L, mask = len(lam), 0
        for i in range(L):
            mask |= 1 << (lam[i] - (i + 1) + Q - self.lo)
        for y in range(self.lo, Q - L):
            mask |= 1 << (y - self.lo)
        return mask

    @staticmethod
    def bilinear(mask, j, k):
        """psi*_j psi_k |mask> mit Vorzeichen (Bitpositionen)."""
        if not (mask >> k) & 1:
            return None, 0
        s = (mask >> (k + 1)).bit_count()
        m2 = mask ^ (1 << k)
        if (m2 >> j) & 1:
            return None, 0
        s += (m2 >> (j + 1)).bit_count()
        return m2 | (1 << j), (-1)**s

    def op(self, Q, koeff, m):
        """sum_y koeff(y) psi*_y psi_(y+m) auf dem Sektor Q (dicht; Bilder außerhalb der Abschneidung entfallen)."""
        zust, idx, _ = self.sek[Q]
        A = np.zeros((len(zust), len(zust)))
        W = self.hi - self.lo
        for c, s in enumerate(zust):
            for k in range(W):
                j = k - m
                if j < 0 or j >= W or not (s >> k) & 1:
                    continue
                neu, sg = self.bilinear(s, j, k)
                if neu is None:
                    continue
                r = idx.get(neu)
                if r is not None:
                    A[r, c] += sg*koeff(j + self.lo)
        return A

    def teilchen_loecher(self, mask):
        nP = nH = 0
        L0 = Fr(0)
        for b in range(self.hi - self.lo):
            y = b + self.lo
            besetzt = (mask >> b) & 1
            if y >= 0 and besetzt:
                nP += 1
                L0 += y + self.alpha
            elif y < 0 and not besetzt:
                nH += 1
                L0 -= y + self.alpha
        return nP, nH, L0


def E_ladung(alpha, Q):
    return Fr(Q*Q, 2) + (alpha - Fr(1, 2))*Q


# ------------------------------------------------------------------ Hilfen: äußere Algebra des Trägers C^5
SUB = [S for p in range(6) for S in itertools.combinations(range(5), p)]
IDX = {S: j for j, S in enumerate(SUB)}
GERADE = [IDX[S] for S in SUB if len(S) % 2 == 0]
X_LAD = np.array([5 - 2*len(S) for S in SUB])                    # X = 5 - 2p (a3, s4, z53)
Z_X = np.diag(1j**(X_LAD % 4))                                    # z_X = exp(i pi X/2), Zentrum von Spin(10)


def Lam(g):
    """Lambda(g) auf Lambda^* C^5: Matrixelemente sind die Minoren det g[T, S]."""
    M = np.zeros((32, 32), complex)
    for S in SUB:
        for T in SUB:
            if len(S) == len(T):
                M[IDX[T], IDX[S]] = np.linalg.det(g[np.ix_(T, S)]) if S else 1.0
    return M


def auf16(A):
    return A[np.ix_(GERADE, GERADE)]


# =================================================================== Teil 1
say("Teil 1: Ist die Uhr R die Z4 der 16?")

# --- 1.1 R auf den Trägermoden
z = sp.symbols("z")
E_MOD = [z**m/(z**4 - 1) for m in range(5)]
R_sym = sp.zeros(5, 5)
for k in range(5):
    poly = sp.Poly(sp.cancel(E_MOD[k].subs(z, sp.I*z)*(z**4 - 1)), z)
    for m in range(5):
        R_sym[m, k] = poly.coeff_monomial(z**m)
R2_sym = R_sym**2
Y_sym = (sp.eye(5) - 5*R2_sym)/12
soll_Y = sp.diag(*[sp.Rational(-1, 3) if m % 2 == 0 else sp.Rational(1, 2) for m in range(5)])
r2_diag = [R2_sym[m, m] for m in range(5)]
ok11 = (R_sym == sp.diag(*[sp.I**m for m in range(5)]) and Y_sym == soll_Y and Y_sym.trace() == 0
        and (Y_sym**2).trace() == sp.Rational(5, 6) and r2_diag.count(1) == 3 and r2_diag.count(-1) == 2)
ZAHLEN["traeger_charaktere"] = [str(R_sym[m, m]) for m in range(5)]
check("R: z -> i z auf dem Träger H0(O(D)), e_m = z^m/(z^4 - 1), exakt: R e_m = i^m e_m (Viertelcharaktere 1, i, -1, -i, 1); "
      "R^2 hat die Eigenräume 3 (+1: e_0, e_2, e_4) und 2 (-1: e_1, e_3); Y = (I - 5R^2)/12 = diag(-1/3, 1/2, -1/3, 1/2, -1/3), "
      "tr Y = 0, tr Y^2 = 5/6", ok11)
R = np.diag([1j**m for m in range(5)])
Rw, Rc = R[np.ix_([1, 3], [1, 3])], R[np.ix_([0, 2, 4], [0, 2, 4])]
check("auf der 2 (R^2 = -1) ist R = diag(i, -i) mit det 1 und Quadrat -1: das Spinorelement von SU(2), derselbe mu4-Modul wie "
      "das fermionische Quarterpaar und das eta-Dublett (z79); auf der 3 ist R = diag(1, -1, 1) mit det -1",
      np.allclose(Rw @ Rw, -np.eye(2)) and abs(np.linalg.det(Rw) - 1) < 1e-12 and abs(np.linalg.det(Rc) + 1) < 1e-12)

# --- 1.2 det R = -1, R = (-1_5) g mit g in G_SM, in beiden 3+2-Schnitten
E5 = np.eye(5)
SCHNITTE = {
    "R^2-Schnitt (Y = (I - 5R^2)/12)": (np.array([E5[1], E5[3]]).T, np.array([E5[0], E5[2], E5[4]]).T),
    "V4-Schnitt (z34, z36)": (np.array([(E5[0] - E5[4])/np.sqrt(2), E5[2]]).T,
                              np.array([E5[1], E5[3], (E5[0] + E5[4])/np.sqrt(2)]).T),
}
G_EICH = -R
ok12, schnitt_info = abs(np.linalg.det(R) + 1) < 1e-12, {}
for name, (Wb, Cb) in SCHNITTE.items():
    B = np.hstack([Wb, Cb])
    gb = B.conj().T @ G_EICH @ B
    blockdiag = np.allclose(gb[:2, 2:], 0) and np.allclose(gb[2:, :2], 0)
    dw, dc = np.linalg.det(gb[:2, :2]), np.linalg.det(gb[2:, 2:])
    Pw = Wb @ Wb.conj().T
    Yc = np.eye(5)/12 + 5/12*(2*Pw - np.eye(5))
    k_y = 0 if abs(dw - 1) < 1e-9 else 1
    h = G_EICH @ expm(-1j*np.pi*k_y*Yc)
    hb = B.conj().T @ h @ B
    ok12 &= (blockdiag and abs(dw*dc - 1) < 1e-9 and np.allclose(R @ Yc, Yc @ R)
             and abs(np.linalg.det(hb[:2, :2]) - 1) < 1e-9 and abs(np.linalg.det(hb[2:, 2:]) - 1) < 1e-9
             and np.allclose(np.sort(np.linalg.eigvalsh(Yc)), [-1/3]*3 + [1/2]*2))
    schnitt_info[name] = {"Blockdet(-R) schwach, Farbe": [complex(np.round(dw, 9)), complex(np.round(dc, 9))],
                          "e^{i pi Y} nötig": bool(k_y)}
ZAHLEN["minus_R_in_G_SM"] = schnitt_info
check("det R = i^(0+1+2+3+4) = -1, und R erhält beide 3+2-Schnitte: -R liegt in S(U(3) x U(2)) = G_SM (Bild in U(5)), im "
      "R^2-Schnitt sogar in SU(3) x SU(2) (Blockdeterminanten 1, 1), im V4-Schnitt mit Blockdeterminanten (-1, -1), also "
      "-R = h e^{i pi Y} mit h in SU(3) x SU(2). Modulo Eichgruppe ist R = -1_5", ok12)

# --- 1.3 Spinorhebung auf Lambda^gerade C^5 = 16
A1, A2 = rng.normal(size=(5, 5)) + 1j*rng.normal(size=(5, 5)), rng.normal(size=(5, 5)) + 1j*rng.normal(size=(5, 5))
ok_cb = np.allclose(Lam(A1 @ A2), Lam(A1) @ Lam(A2))
LamR, Lamg = Lam(R), Lam(G_EICH)
ok13, zweige = ok_cb, {}
for s in (1j, -1j):                                               # s^2 = det R: die beiden Hebungen
    Rt = LamR/s
    M = auf16(Rt @ np.linalg.inv(Lamg))
    skalar = np.allclose(M, M[0, 0]*np.eye(16))
    ist_zX = np.allclose(M, auf16(Z_X)) or np.allclose(M, auf16(np.linalg.inv(Z_X)))
    ok13 &= skalar and ist_zX and abs(abs(M[0, 0]) - 1) < 1e-12
    zweige[str(s)] = complex(np.round(M[0, 0], 12))
ZAHLEN["R_tilde_mal_Eichinverse_auf_16"] = zweige
check("Spinorhebung R~ = det(R)^(-1/2) Lambda(R) auf Lambda^gerade C^5 (beide Zweige; Cauchy-Binet geprüft): "
      "R~ Lambda(-R)^(-1) ist auf allen 16 der Skalar -+i = z_X^(-+1), z_X = exp(i pi X/2) mit X = 5 - 2p = 1 mod 4, das "
      "Zentrum von Spin(10): die Vierteldrehung der Marken ist auf der 16 die anomalierelevante Spin-Z4, modulo G_SM", ok13)
nu_c = IDX[()]
naiv_eichtrivial = np.allclose(auf16(LamR), auf16(Lamg))
check("die naive Hebung Lambda(R) ist auf der 16 gleich Lambda(-R), also eichtrivial; eichinvariant sichtbar am G_SM-Singulett "
      "nu^c = Lambda^0: dort wirkt Lambda(R) als 1, R~ als +-i (z_X), jedes Lambda(g) mit g in SU(5) als 1",
      naiv_eichtrivial and abs(LamR[nu_c, nu_c] - 1) < 1e-12 and abs(abs((LamR/1j)[nu_c, nu_c]) - 1) < 1e-12
      and abs((LamR/1j)[nu_c, nu_c] - 1) > 1 and abs(Lamg[nu_c, nu_c] - 1) < 1e-12)

# --- 1.4 Ladungen der 16 und mod-16-Bilanz
Ydiag = [float(Y_sym[m, m]) for m in range(5)]
NAMEN = {0.0: "nu^c", 1.0: "e^c", round(1/6, 9): "Q", round(-2/3, 9): "u^c", -0.5: "L", round(1/3, 9): "d^c"}
tab = {}
for S in SUB:
    if len(S) % 2:
        continue
    nm = NAMEN[round(sum(Ydiag[i] for i in S), 9)]
    r = sum(S) % 4
    d = tab.setdefault(nm, {"Lambda(R)": [], "R~ (Zweig i^X)": [], "X mod 4": []})
    d["Lambda(R)"].append(r)
    d["R~ (Zweig i^X)"].append((r + 1) % 4)
    d["X mod 4"].append((5 - 2*len(S)) % 4)
for d in tab.values():
    for k in d:
        d[k] = sorted(d[k])
ZAHLEN["ladungen_der_16"] = tab
xs = [(5 - 2*len(S)) % 4 for S in SUB if len(S) % 2 == 0]
nu = sum(1 if x == 1 else (-1 if x == 3 else 0) for x in xs)
nu_ohne = nu - 1                                                  # ohne nu^c (X = 5)
ZAHLEN["nu_mod16"] = {"eine Familie": nu, "ohne nu^c": nu_ohne, "drei Familien": 3*nu}
check("mod-16-Bilanz der Klasse von R (= z_X): alle 16 haben X = 1 mod 4, nu = 16 = 0 mod 16 je Familie, drei Familien "
      "48 = 0 mod 16; ohne nu_R 15 (Anomalie): die Kette mu_4 -> 16 (a3, s4, z53) gilt für die Vierteldrehung der Marken",
      nu == 16 and nu % 16 == 0 and (3*nu) % 16 == 0 and nu_ohne % 16 == 15)
Yc_R2 = np.diag(Ydiag)
LamR2 = Lam(R @ R)
weak_par = np.diag([(-1.0)**len(set(S) & {1, 3}) for S in SUB])
Rt_i = LamR/(-1j)                                                 # Zweig mit R~ Lambda(-R)^(-1) = z_X
Wv, Cv = SCHNITTE["V4-Schnitt (z34, z36)"]
Bv = np.hstack([Wv, Cv])
R2v = Bv.conj().T @ (R @ R) @ Bv
ok14 = (np.allclose(R @ R, expm(6j*np.pi*Yc_R2)) and np.allclose(LamR2, weak_par)
        and np.allclose(auf16(Rt_i @ Rt_i), -auf16(LamR2)) and np.allclose(auf16(LamR @ LamR), auf16(LamR2))
        and np.allclose(R2v[:2, :], np.eye(5)[:2, :]) and abs(np.linalg.det(R2v[2:, 2:]) - 1) < 1e-12)
gerade_Rt = sum(1 for d in tab.values() for c in d["R~ (Zweig i^X)"] if c % 2 == 0)
for nm, d in tab.items():
    say(f"   {nm:5s}: Lambda(R) {d['Lambda(R)']}, R~ {d['R~ (Zweig i^X)']}, X mod 4 {d['X mod 4']}")
check("R^2 = exp(6 pi i Y) ist im R^2-Schnitt das Zentrum von SU(2)_L ((-1) auf den Dubletts), im V4-Schnitt ein SU(3)-Element: "
      "in beiden ein Eichelement. R~^2 = (-1)^F Lambda(R^2) und Lambda(R)^2 = Lambda(R^2): modulo Eichung ist R~^2 = (-1)^F "
      "(Spin-Z4), Lambda(R)^2 = +1. Überträgt man die Regel R^2 = (-1)^F des Originals auf die 16, wählt sie die "
      "Spinorhebung; R~ selbst hat auf Q und L gerade Ladungen (8 von 16), ist also erst als Klasse modulo G_SM eine Spin-Z4",
      ok14 and gerade_Rt == 8)


def spin_z4_klasse(chars):
    """Gibt es c mit c^N = det R (c 1_N in der Klasse von R modulo SU(N)), dessen Spinorhebung auf Lambda^gerade C^N ein
    Skalar mit Quadrat -1 ist? Hebung von c 1_N auf Lambda^p: c^p det^(-1/2), skalar über p gerade genau für c^2 = 1."""
    N = len(chars)
    det = np.prod(chars)
    for k in range(N):
        c = np.exp((np.log(det + 0j) + 2j*np.pi*k)/N)
        skalar = N < 2 or abs(c*c - 1) < 1e-9
        if skalar and abs(1/c**N + 1) < 1e-9:                      # Quadrat der Hebung = det(c 1)^(-1) = c^(-N)
            return True
    return False


N_ok = [N for N in range(1, 17) if spin_z4_klasse([1j**m for m in range(N)])]
ZAHLEN["N_mit_Spin_Z4_Klasse_bis_16"] = N_ok
check("Gegenproben: N Trägermoden mit Viertelcharakteren i^0..i^(N-1) tragen modulo SU(N) genau dann eine Spin-Z4-Klasse auf "
      f"Lambda^gerade C^N, wenn det = -1 und N ungerade (bis 16: N = {N_ok}, also N = 5 mod 8); N = 4 (det -1, aber N gerade) "
      "und verschobene Charaktere i^1..i^5 (det -i) nicht: das kleinste ist N = 5 = h0(O(D)) der Viermarkennaht",
      N_ok[:1] == [5] and 4 not in N_ok and not spin_z4_klasse([1j**m for m in range(1, 6)]))

# --- 1.5 Hardy-Marken k = 4x des Originals -> Trägermoden (Residuenklassen)
w_ = sp.symbols("w")
ok_res = True
for m in range(5):
    s0 = sp.series(E_MOD[m], z, 0, 26).removeO()
    k0 = {k: s0.coeff(z, k) for k in range(26)}
    sI = sp.series(E_MOD[m].subs(z, 1/w_), w_, 0, 26).removeO()
    kI = {-j: sI.coeff(w_, j) for j in range(26)}
    ok_res &= all(c == 0 or (k - m) % 4 == 0 for k, c in k0.items())
    ok_res &= all(c == 0 or (k - m) % 4 == 0 for k, c in kI.items())
    ok_res &= any(c != 0 for c in k0.values()) and any(c != 0 for c in kI.values())
check("die Laurentreihen von e_m = z^m/(z^4 - 1) um 0 (Hardy-Hälfte k >= 0) und um unendlich (k < 0) tragen nur Hardy-Marken "
      "k = m mod 4: jede Trägermode ist eine Residuenklasse k mod 4 (e_0 und e_4 = 1 + e_0 beide in der Klasse 0)", ok_res)
ok_eq, bild = True, set()
for alpha in (Fr(1, 4), Fr(3, 4)):
    for r in range(-20, 21):
        x = r + alpha
        k = int(4*x)
        klasse = k % 4
        ok_eq &= abs(np.exp(2j*np.pi*float(x)) - complex(R_sym[klasse, klasse])) < 1e-12
        bild.add(klasse)
bild_r2 = sorted(int(R2_sym[k, k]) for k in bild)
check("die Einteilchenmoden des Originals (x in Z + 1/4, Z + 3/4, |x| <= 20) haben ungerade Hardy-Marken k = 4x; die Abbildung "
      "x -> e_(4x mod 4) ist R-äquivariant (e^{2 pi i x} = i^k = Charakter von e_(k mod 4)) und trifft genau e_1, e_3, also "
      "den R^2 = -1-Eigenraum (die 2); der dreidimensionale R^2 = +1-Eigenraum (e_0, e_2, e_4; Klassen k = 0, 2 mod 4, "
      "alpha = 0, 1/2) bleibt ohne Einteilchenquelle", ok_eq and bild == {1, 3} and bild_r2 == [-1, -1])
befund("Im R^2-Schnitt ist das Bild der Originalmoden das schwache Dublett (Y = 1/2); im V4-Schnitt (z34) liegt {z, z^3} in "
       "der Farbe. In beiden Schnitten fehlt der R^2 = +1-Eigenraum; was fehlt, sind die Quarterklassen alpha = 0 und 1/2", True)

# --- 1.6 Warum die 3 im Original keine fermionische Quelle hat
MODEN8 = [(Fr(1, 4), x) for x in (Fr(-7, 4), Fr(-3, 4), Fr(1, 4), Fr(5, 4))] + \
         [(Fr(3, 4), x) for x in (Fr(-5, 4), Fr(-1, 4), Fr(3, 4), Fr(7, 4))]
rexp, par = np.zeros(256, int), np.zeros(256, int)
for s in range(256):
    L0 = Fr(0)
    for b, (_, x) in enumerate(MODEN8):
        bes = (s >> b) & 1
        if x > 0 and bes:
            L0 += x
        elif x < 0 and not bes:
            L0 -= x
    rexp[s] = int(4*L0) % 4                                       # R = e^{2 pi i L_0} = i^(4 L_0)
    par[s] = bin(s).count("1") % 2
char = (rexp[:, None] - rexp[None, :]) % 4                        # Ad_R auf |s><t|
ungerade = (par[:, None] + par[None, :]) % 2 == 1
ok16a = set(np.unique(char[ungerade])) == {1, 3} and set(np.unique(char[~ungerade])) == {0, 2}
check("CAR-Algebra des Originals (8 Moden aus beiden Quartern, alle 65536 Matrixeinheiten): Ad_R hat auf ungeraden "
      "(fermionischen) Operatoren genau die Charaktere {i, -i}, auf geraden genau {1, -1}. Keine fermionische Mode des "
      "Originals hat einen Charakter +-1 der 3", ok16a)
min_gerade, min_farbe, nie_null = 9, 9, True
for d in itertools.product(range(4), repeat=5):
    if sum(d) % 4:
        continue                                                  # diagonales mu4-Element mit det 1 (in G_SM)
    ex = [(m + d[m]) % 4 for m in range(5)]
    n_ger = sum(1 for e in ex if e % 2 == 0)
    n_far = sum(1 for m in (0, 2, 4) if ex[m] % 2 == 0)
    min_gerade, min_farbe = min(min_gerade, n_ger), min(min_farbe, n_far)
    nie_null &= n_ger > 0
ZAHLEN["determinantenschranke"] = {"min. Charaktere mit Quadrat +1": min_gerade, "davon min. im Farbblock": min_farbe}
check("Determinantenschranke: det(R g) = -1 für jedes g in SU(5), fünf Charaktere mit Quadrat -1 gäben aber det^2 = -1. Über "
      f"alle 256 diagonalen mu4-Eichelemente (in G_SM) bleibt mindestens {min_gerade} Charakter mit Quadrat +1 (nie 0); er "
      f"kann in den Farb- oder in den Schwachblock gelegt werden (Minimum im Farbblock {min_farbe}). Modulo Eichung sind also "
      "höchstens vier der fünf Trägermoden fermionisch, nie alle", nie_null and min_gerade == 1 and min_farbe == 0)
sekt = []
for Qa in range(-3, 4):
    for Qb in range(-3, 4):
        L0 = E_ladung(Fr(1, 4), Qa) + E_ladung(Fr(3, 4), Qb)
        sekt.append((L0, Qa, Qb, (Qa - Qb) % 4, (Qa + Qb) % 2))
sekt.sort()
bos = [(str(L0), Qa, Qb, ["1", "i", "-1", "-i"][c]) for L0, Qa, Qb, c, f in sekt if f == 0][:5]
fer = [(str(L0), Qa, Qb, ["1", "i", "-1", "-i"][c]) for L0, Qa, Qb, c, f in sekt if f == 1][:5]
ZAHLEN["niedrigste_Sektoren (L0, Q_a, Q_b, R)"] = {"bosonisch": bos, "fermionisch": fer}
say(f"   niedrigste bosonische Sektoren (L0, Q_a, Q_b, R): {bos}")
say(f"   niedrigste fermionische Sektoren: {fer}")
befund("Die Charaktere der 3 (1, -1, 1) tragen im Original nur bosonische Sektoren (Q_a + Q_b gerade): Paare psi*_a psi*_b "
       "und die Ströme (k = 0), psi*_a psi_b, psi*_a psi*_a, psi*_b psi*_b (k = 2); die drei niedrigsten bosonischen Sektoren "
       "(0, 0), (1, -1), (1, 1) bei L_0 = 0, 1/2, 1 haben gerade die Charaktere 1, -1, 1 (Beobachtung, keine Zuordnung). Die 3 "
       "käme also aus Fermipaaren (Klassen alpha = 0 und 1/2), nicht aus Einteilchenmoden", True)
split_ok = True
for rho0 in (1, 1j, -1, -1j):
    r2 = [(rho0*np.prod([1j**i for i in S]))**2 for S in SUB if len(S) % 2 == 0]
    n_ferm = sum(1 for v in r2 if abs(v + 1) < 1e-9)
    split_ok &= n_ferm == 8
r2X = np.diag(auf16(Z_X @ Z_X))
check("Statistik der 16 bei exakten Trägercharakteren: Realisiert man die fünf Moden im Original als R-Eigenoperatoren mit "
      "den Charakteren 1, i, -1, -i, 1, dann ist für jede Referenzphase rho_0 des Zustands nu^c genau die Hälfte der 16 "
      "R^2-ungerade: im R^2-Schnitt haben die Dubletts Q, L die andere Statistik als nu^c, e^c, u^c, d^c (gescoptes Negativ); "
      "nur die Klasse z_X gibt allen 16 R^2 = -1", split_ok and np.allclose(r2X, -1))

# =================================================================== Teil 2
say("\nTeil 2: Spin-Statistik")

# --- 2.1 geometrische Spinorhebung des eta-Dubletts (exakt)
u = sp.symbols("u")


def eta_matrix(vorz):
    """Transport auf H0(E, eta (x) pi^*O(1)) = span{1, f}, f = w/(u^2 - 1), Kurve w^2 = 1 - u^4 (wie z79):
    (T b)(u, w) = f(u, w) b(-i u, vorz w). Elemente als A(u) + B(u) w, Reduktion w^2 = 1 - u^4."""
    F = 1/(u**2 - 1)
    basis = [(sp.Integer(1), sp.Integer(0)), (sp.Integer(0), F)]
    M = sp.zeros(2, 2)
    for j, (A, Bw) in enumerate(basis):
        A2, B2 = A.subs(u, -sp.I*u), vorz*Bw.subs(u, -sp.I*u)
        prod = (sp.cancel(F*B2*(1 - u**4)), sp.cancel(F*A2))      # (0 + F w)(A2 + B2 w)
        a, c = sp.cancel(prod[0]), sp.cancel(prod[1]*(u**2 - 1))
        assert a.free_symbols == set() and c.free_symbols == set()
        M[0, j], M[1, j] = a, c
    return M


G_ETA, G_ETA0 = eta_matrix(-1), eta_matrix(1)
ew_eta = sorted(G_ETA.eigenvals().keys(), key=lambda v: float(sp.arg(v)))
check("eta-Dublett span{1, w/(u^2 - 1)} exakt: die Spinorhebung (u, w) -> (-iu, -w) ist [[0, -1], [1, 0]] mit Eigenwerten "
      "{-i, i} und Quadrat -1 (2 pi-Drehung = -1); die andere Hebung (w -> w) ist [[0, 1], [1, 0]] mit Quadrat +1",
      G_ETA == sp.Matrix([[0, -1], [1, 0]]) and G_ETA**2 == -sp.eye(2) and ew_eta == [-sp.I, sp.I]
      and G_ETA0 == sp.Matrix([[0, 1], [1, 0]]) and G_ETA0**2 == sp.eye(2))
lam_a, lam_b = complex(ew_eta[1]), complex(ew_eta[0])              # Quarter a (alpha = 1/4) <-> i, b <-> -i (z79)
exp_von = lambda v: int(round(np.angle(v)/(np.pi/2))) % 4

# --- 2.2 alle Sektoren |Q_a|, |Q_b| <= 10 mit Oszillatorinhalt bis Stufe 5 je Quarter
QMAX, NLEV = 10, 5
daten = {}
for alpha in (Fr(1, 4), Fr(3, 4)):
    qq = Quarter(alpha, range(-QMAX, QMAX + 1), NLEV, 0)
    rows = []
    for Q in range(-QMAX, QMAX + 1):
        zust, _, stufen = qq.sek[Q]
        for s, n in zip(zust, stufen):
            nP, nH, L0 = qq.teilchen_loecher(s)
            assert L0 == E_ladung(alpha, Q) + n and nP - nH == Q
            rows.append((Q, nP, nH, int(4*L0) % 4, (nP + nH) % 2))
    daten[alpha] = np.array(rows)
Da, Db = daten[Fr(1, 4)], daten[Fr(3, 4)]
QA, QB = Da[:, 0][:, None], Db[:, 0][None, :]
R_exp = (Da[:, 3][:, None] + Db[:, 3][None, :]) % 4
F_par = (Da[:, 4][:, None] + Db[:, 4][None, :]) % 2


def spin_exp(ea, eb):
    return ((Da[:, 1] - Da[:, 2])[:, None]*ea + (Db[:, 1] - Db[:, 2])[None, :]*eb) % 4


ea, eb = exp_von(lam_a), exp_von(lam_b)
G_exp = spin_exp(ea, eb)
anz = R_exp.size
ok22 = (np.array_equal(R_exp, G_exp) and np.array_equal(R_exp, (QA - QB) % 4)
        and np.array_equal((2*R_exp) % 4, 2*F_par))
ZAHLEN["spin_statistik_zustaende"] = int(anz)
check(f"auf allen {anz} Zuständen (Q_a, Q_b in [-10, 10], Oszillatorstufe <= 5 je Quarter; L_0, Teilchen und Löcher aus den "
      "Besetzungen): das Produkt der eta-Spinoreigenwerte über alle Quanten (Teilchen lambda, Loch lambda-quer; a <-> i, "
      "b <-> -i) ist gleich R = e^{2 pi i L_0} = i^(Q_a - Q_b), und sein Quadrat, die 2 pi-Drehung der 4D-Lesart, ist gleich "
      "(-1)^F: Fermionen sind genau die eta-Spinoren", ok22)
treu_R, treu_F = [], []
for e1, e2 in itertools.product(range(4), repeat=2):
    Ge = spin_exp(e1, e2)
    if np.array_equal(Ge, R_exp):
        treu_R.append((e1, e2))
    if np.array_equal((2*Ge) % 4, 2*F_par):
        treu_F.append((e1, e2))
naht = spin_exp(1, 0)                                              # Nahtdublett V_1: Eigenwerte (i, 1) (z79)
skal = spin_exp(1, 1)
anteil_naht = float(np.mean((2*naht) % 4 != 2*F_par))
anteil_skal = float(np.mean(skal != R_exp))
ZAHLEN["zuordnungen"] = {"treu zu R": treu_R, "treu zur Statistik": treu_F, "Nahtdublett: Statistik falsch (Anteil)": anteil_naht,
                         "skalar (i, i): R falsch (Anteil)": anteil_skal}
check(f"alle 16 Zuordnungen (i^e_a, i^e_b) der Quartergrade durchgemustert: R treu nur {treu_R} (das eta-Spinordublett "
      f"{{1, 3}}, eindeutig wie in z79), die Statistik treu nur für die vier ungeraden Paare {treu_F}; Gegenproben: das "
      f"Nahtdublett V_1 (i, 1) verfehlt (-1)^F auf {anteil_naht:.1%}, die skalare Zuordnung (i, i) verfehlt R auf "
      f"{anteil_skal:.1%} der Zustände", treu_R == [(1, 3)] and sorted(treu_F) == [(1, 1), (1, 3), (3, 1), (3, 3)]
      and anteil_naht > 0.4 and anteil_skal > 0.4)

# --- 2.3 innere U(1): J_z = (Q_a - Q_b)/2; volle eta-SU(2) im KMS-Zustand
Jz2 = QA - QB                                                     # 2 J_z
ok23 = np.array_equal(Jz2 % 2 == 1, F_par == 1) and np.array_equal((Jz2 % 4), R_exp)
niedrig = [(L0, Qa, Qb, Fr(Qa - Qb, 2)) for L0, Qa, Qb, c, f in sekt if f == 1][:4]
ZAHLEN["innere_Helizitaet_niedrigste_Fermisektoren (L0, Q_a, Q_b, J_z)"] = [tuple(map(str, t)) for t in niedrig]
check("R = e^{i pi J_z} mit J_z = (Q_a - Q_b)/2: J_z ist genau auf den fermionischen Zuständen halbzahlig, e^{2 pi i J_z} = (-1)^F; "
      "die niedrigsten Fermisektoren haben J_z = +1/2 bei L_0 = 1/4 ((1, 0), (0, -1)) und J_z = -1/2 bei L_0 = 3/4 (ihre "
      "Konjugierten)", ok23 and [str(t[3]) for t in niedrig] == ["1/2", "1/2", "-1/2", "-1/2"]
      and [str(t[0]) for t in niedrig] == ["1/4", "1/4", "3/4", "3/4"])


def hs_summe(theta, X):
    """Powers-Stoermer: sum über Paare (x, -x), x in Z + 1/4, |x| <= X, von ||sqrt C - sqrt C'||^2 + ||sqrt(1-C) - sqrt(1-C')||^2
    für C' = U C U^+, U = Drehung um theta im eta-Dublett (endlich <=> implementierbar)."""
    U = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])

    def wurzel(A):
        ev, V = np.linalg.eigh(A)
        return V @ np.diag(np.sqrt(np.clip(ev, 0, None))) @ V.T
    tot = 0.0
    for r in range(-X - 1, X + 1):
        x = r + 0.25
        if abs(x) > X:
            continue
        C = np.diag([0.5*(1 - np.tanh(BETA*x/2)), 0.5*(1 + np.tanh(BETA*x/2))])
        C2 = U @ C @ U.T
        tot += np.sum((wurzel(C) - wurzel(C2))**2) + np.sum((wurzel(np.eye(2) - C) - wurzel(np.eye(2) - C2))**2)
    return tot


hs = {X: hs_summe(np.pi/4, X) for X in (10, 20, 40)}
torus = hs_summe(0.0, 40)
ZAHLEN["eta_SU2_HS_Summen"] = {str(k): v for k, v in hs.items()}
check("die volle eta-SU(2) (Drehung um pi/4 auf allen Paaren (x, -x)) ist im KMS-Zustand nicht implementierbar: die "
      "Powers-Stoermer-Summe wächst linear (X = 10, 20, 40: Verhältnisse ~2), nur der Torus (erzeugt von Q_a - Q_b, darin R) "
      "erhält den Zustand (Summe 0). Innen wirken also U(1) und Z4, keine Drehgruppe",
      abs(hs[20]/hs[10] - 2) < 0.1 and abs(hs[40]/hs[20] - 2) < 0.1 and torus < 1e-20)

# --- 2.4 Gegenprobe Außenbild: die 2 pi-Drehung der Lorentzgruppe aus z79
# z79: g in SL(2,C) wirkt als pi -> g pi, omega -> (g^+)^-1 omega auf T = span{a_i, b_i^+}; g = -1 ist a -> -a, b -> -b,
# implementiert durch U = (-1)^(N_a + N_b). Auf den Strömen des Originals: U_alle = (-1)^(sum_m N_m), U_1 = (-1)^(N_1).
qa = Quarter(Fr(1, 4), [0, 1], 9, 9)
Jq = {Q: {m: qa.op(Q, lambda y: 1.0, m) for m in list(range(-9, 0)) + list(range(1, 10))} for Q in (0, 1)}


def paritaet(Nop):
    ew, V = np.linalg.eigh(Nop)
    assert np.allclose(ew, np.round(ew), atol=1e-9)
    return V @ np.diag((-1.0)**np.round(ew)) @ V.T


U_alle = {Q: paritaet(sum(Jq[Q][-m] @ Jq[Q][m]/m for m in range(1, 10))) for Q in (0, 1)}
U_eins = {Q: paritaet(Jq[Q][-1] @ Jq[Q][1]) for Q in (0, 1)}
vak = qa.sek[0][0][0]
erw_alle, erw_eins, ok24 = [], [], True
for n in range(10):
    j = qa.sek[1][1][vak | (1 << (n - qa.lo))]                    # psi*_(n + 1/4) |0> = |Q = 1, lambda = (n)>
    v = np.zeros(len(qa.sek[1][0]))
    v[j] = 1
    ea_, e1_ = v @ U_alle[1] @ v, v @ U_eins[1] @ v
    erw_alle.append(round(float(ea_), 12))
    erw_eins.append(round(float(e1_), 12))
    soll_alle = 1.0 if n == 0 else (-1.0 if n == 1 else 0.0)
    soll_eins = sum((-2)**k/math.factorial(k) for k in range(n + 1))
    ok24 &= abs(ea_ - soll_alle) < 1e-10 and abs(e1_ - soll_eins) < 1e-10
ZAHLEN["aeussere_2pi_auf_psi*_(n+1/4)|0>"] = {"<U_alle>": erw_alle, "<U_1>": erw_eins}
jm1 = Jq[0][-1][:, 0]
osz_vz = float(jm1 @ U_alle[0] @ jm1/(jm1 @ jm1))
ungerade_bosonen = int(sum(1 for e in np.linalg.eigvalsh(sum(Jq[0][-m] @ Jq[0][m]/m for m in range(1, 10)))
                           if int(round(e)) % 2 == 1))
ZAHLEN["bosonische Oszillatorzustaende mit aeusserem Vorzeichen -1 (Stufe <= 9)"] = ungerade_bosonen
say(f"   <U_alle> auf psi*_(n+1/4)|0>, n = 0..9: {erw_alle}")
say(f"   <U_1>    auf psi*_(n+1/4)|0>, n = 0..9: {erw_eins}")
check("Gegenprobe Außenbild (z79): die 2 pi-Drehung U = (-1)^(Zahl der Stromquanten) gibt dem Bosonzustand J_-1|0> das "
      f"Vorzeichen {osz_vz:+.0f} (R^2 = +1; {ungerade_bosonen} solche Bosonzustände bis Stufe 9), und auf den Fermimoden "
      "psi*_(n+1/4)|0> ist sie nicht definit: <U_alle> = 1, -1, 0, 0, ... (exakt 0 für n >= 2: gleich viele Permutationen mit "
      "gerader und ungerader Zyklenzahl), <U_1> = sum_(k<=n) (-2)^k/k! -> e^(-2). Da U|0> = |0>, ist U psi*_x U^-1 nicht "
      "+-psi*_x. Gescoptes Negativ: die Hülle eines Vertexoperators mischt beide Vorzeichen, die äußere Lorentzgruppe gibt "
      "Fermionen keine konsistente Helizität", ok24 and abs(osz_vz + 1) < 1e-12 and ungerade_bosonen > 0)
befund("Spin-Statistik gilt innen exakt (eta-Spinor, R^2 = (-1)^F, J_z = (Q_a - Q_b)/2 halbzahlig genau für Fermionen), außen "
       "nicht: die Lorentzgruppe aus z79 wirkt auf der Stromhülle, die R nicht sieht, und gibt Strombosonen halbzahlige "
       "Twistorhelizität. Eine Spin x_Z2 Z4-Verklebung (äußeres -1 = R^2) verlangte (-1)^(Zahl der Quanten) = (-1)^F auf allen "
       "Zuständen; die Vertexoperatoren des Originals verletzen sie", True)

# =================================================================== Teil 3
say("\nTeil 3: Universalität -- dieselbe Oszillatorhülle in jedem Ladungssektor")
NU = 14
QS = range(-3, 4)
ok31, ok_ccr, ok32_q, ok32_w, ok33, ok33b, ok34 = True, True, True, True, True, True, True
antw, kreuz_max, tensor = {}, 0.0, {}
Jref = None
for alpha in (Fr(1, 4), Fr(3, 4)):
    qq = Quarter(alpha, QS, NU, NU)
    J = {Q: {m: qq.op(Q, lambda y: 1.0, m) for m in list(range(-NU, 0)) + list(range(1, NU + 1))} for Q in QS}
    if Jref is None:
        Jref = J[0]
    ok31 &= all(np.array_equal(J[Q][m], Jref[m]) for Q in QS for m in (1, 2, 3, -1, -2, -3))
    stufen = qq.sek[0][2]
    for m in (1, 2, 3):
        C = J[0][m] @ J[0][-m] - J[0][-m] @ J[0][m]
        tief = stufen <= NU - m
        ok_ccr &= np.allclose(C[np.ix_(tief, tief)], m*np.eye(int(tief.sum())))
    Losc = {n: sum(J[0][-k] @ J[0][n + k] for k in range(1, NU - n + 1)) + 0.5*sum(J[0][n - k] @ J[0][k] for k in range(1, n))
            for n in (1, 2, 3)}
    for Q in QS:
        mu = Q + float(alpha) - 0.5
        for n in (1, 2, 3):
            Lf = qq.op(Q, lambda y, n=n: float(y + alpha + Fr(n, 2)), n)
            Lo = sum(J[Q][-k] @ J[Q][n + k] for k in range(1, NU - n + 1)) + 0.5*sum(J[Q][n - k] @ J[Q][k] for k in range(1, n))
            ok33 &= np.abs(Lf - Lo - mu*J[Q][n]).max() < 1e-12
            ok33b &= np.array_equal(Lo, Losc[n])
            for beta in (BETA, 2.0):
                wgt = np.exp(-beta*stufen)
                wgt /= wgt.sum()
                ev = lambda A, B=None: float(np.sum(wgt*np.sum(A*(A if B is None else B), axis=0)))
                T = ev(Lf)
                A_n, B_n = ev(Lo), ev(J[Q][n])
                tensor.setdefault((n, beta), []).append((float(alpha), Q, T - 2*B_n*float(E_ladung(alpha, Q)), A_n + B_n/16))
                ok34 &= abs(ev(Lo, J[Q][n])) < 1e-12*np.sqrt(A_n*B_n)
                ok34 &= abs(ev(J[Q][n], Lf) - mu*B_n) < 1e-12*abs(mu*B_n)
    for beta in (BETA, 2.0):
        wgt = np.exp(-beta*stufen)
        wgt /= wgt.sum()
        ev = lambda A: float(np.sum(wgt*np.sum(A*A, axis=0)))
        n1, n2, n3 = ev(J[0][1]), ev(J[0][2]), ev(J[0][3])
        vier = ev(J[0][1] @ J[0][1]) - 2*n1*n1
        vier_m = ev(J[0][1] @ J[0][2]) - n1*n2
        for Q in QS:
            werte = (ev(J[Q][1]), ev(J[Q][2]), ev(J[Q][3]), ev(J[Q][1] @ J[Q][1]) - 2*ev(J[Q][1])**2,
                     ev(J[Q][1] @ J[Q][2]) - ev(J[Q][1])*ev(J[Q][2]))
            ok32_q &= np.allclose(werte, (n1, n2, n3, vier, vier_m), rtol=1e-13, atol=1e-300)
        tol = 1e-12 if beta == BETA else 1e-6
        soll = [m/(np.exp(beta*m) - 1) for m in (1, 2, 3)]
        ok32_w &= all(abs(v/s_ - 1) < tol for v, s_ in zip((n1, n2, n3), soll))
        ok32_w &= abs(vier) < tol*2*n1*n1 and abs(vier_m) < tol*n1*n2
        antw[(float(alpha), round(beta, 6))] = {"<J_-1 J_1>": n1, "<J_-2 J_2>": n2, "verbunden 4pt (1,1)": vier,
                                                "verbunden 4pt (1,2)": vier_m}
    # innen-außen: volle KMS-Verteilung über die Sektoren (Produkt), Kovarianz von Q^2 und J_-1 J_1
    for beta in (BETA, 2.0):
        gew, q2, nn = [], [], []
        for Q in QS:
            wQ = np.exp(-beta*(float(E_ladung(alpha, Q)) + qq.sek[Q][2]))
            gew.append(wQ)
            q2.append(np.full(len(wQ), Q*Q, float))
            nn.append(np.sum(J[Q][1]*J[Q][1], axis=0))
        gew = np.concatenate(gew)
        gew /= gew.sum()
        q2, nn = np.concatenate(q2), np.concatenate(nn)
        cov = gew @ (q2*nn) - (gew @ q2)*(gew @ nn)
        kreuz_max = max(kreuz_max, abs(cov)/max(gew @ q2, 1e-300)/max(gew @ nn, 1e-300))
ZAHLEN["stromantworten_je_Quarter_und_beta"] = {f"alpha={a}, beta={b}": v for (a, b), v in antw.items()}
check("in jedem Ladungssektor Q = -3..3 und für beide Quarter sind die fermionischen Matrixelemente von J_+-1, J_+-2, J_+-3 in "
      "der Maya-Basis exakt dieselben, und [J_m, J_-m] = m (m = 1, 2, 3): dieselbe Stromalgebra und derselbe Fockraum", ok31 and ok_ccr)
check("bedingte KMS-Antworten je Sektor (beta = 4 pi und Kontrolle beta = 2, Stufe <= 14): <J_-m J_m>_Q, die verbundenen "
      "Vierpunktantworten <J_-1^2 J_1^2> - 2<J_-1 J_1>^2 und <J_-2 J_-1 J_1 J_2> - <J_-1 J_1><J_-2 J_2> sind in allen sieben "
      "Sektoren gleich (relativ 1e-13)", ok32_q)
check("und stimmen mit der gaußschen Hülle überein: <J_-m J_m> = m/(e^{beta m} - 1), verbundene Vierpunktantworten 0 "
      "(beta = 4 pi relativ 1e-12, beta = 2 Abschneidefehler < 1e-6); die Kovarianz von Q^2 (innen) und J_-1 J_1 (außen) im "
      f"vollen KMS-Zustand verschwindet (relativ {kreuz_max:.0e}): der Zustand ist ein Produkt, die Materie sourct die Hülle "
      "nicht", ok32_w and kreuz_max < 1e-12)
check("Sugawara je Sektor, fermionisch: L_n = sum_x (x + n/2) psi*_x psi_(x+n) ist auf allen Zuständen gleich L_n^osc + "
      "(Q + alpha - 1/2) J_n (n = 1, 2, 3, beide Quarter, Fehler < 1e-12), und L_n^osc = sum J_-k J_(n+k) + ... ist in jedem "
      "Sektor dieselbe Matrix: T_osc ist universell, die Kopplung an den Sektor ist der lineare Term mu_Q J_n", ok33 and ok33b)
konst_ok, konst = True, {}
for (n, beta), rows in tensor.items():
    vals = np.array([r[2] for r in rows])
    pred = np.array([r[3] for r in rows])
    skala = np.abs(vals).max()
    konst_ok &= (vals.max() - vals.min()) < 1e-12*skala and np.allclose(vals, pred, rtol=1e-12)
    konst[f"n={n}, beta={round(beta, 6)}"] = float(vals.mean())
ZAHLEN["tensorantwort_konstante (<L_-n L_n>_Q - 2 B_n E(Q))"] = konst
check("Tensorantwort: <L_-n L_n>_Q = A_n + B_n/16 + 2 B_n E_alpha(Q) mit B_n = <J_-n J_n> (n = 1, 2, 3; beta = 4 pi und 2): nach "
      "Abzug von 2 B_n E(Q) bleibt in allen 14 Sektoren beider Quarter dieselbe Konstante (relativ 1e-12; (alpha - 1/2)^2 = "
      "1/16 für beide Quarter); die lineare Antwort <J_-n L_n>_Q = mu_Q B_n ist ladungsproportional, der Kreuzterm "
      "<L_-n^osc J_n> verschwindet", konst_ok and ok34)
befund("Universalität im genauen Sinn: Jeder Ladungssektor sieht denselben Oszillator-Fockraum, denselben Zustand und dasselbe "
       "T_osc; die Tensorantwort hängt vom inneren Inhalt nur über die innere Energie E(Q) ab, mit einem für beide Quarter "
       "gleichen Koeffizienten 2 B_n (Teleskopsumme, z79). Das ist eine äquivalenzprinzipartige Kopplung im 2D-Tensorkanal. "
       "Grenzen: Im KMS-Zustand ist die Kreuzkorrelation null (keine Quelle, G nicht bestimmt), die lineare Kopplung mu_Q J_n "
       "ist ladungs-, nicht energieproportional, und T(z) ist kein 4D-Spin 2 (z75)", True)

# =================================================================== Teil 4
say("\nTeil 4: Bilanz für S1")
STATUS = {
    "R auf dem Träger = Z4 der 16": ("E", "R = (-1_5) g, g in G_SM; Spinorhebung = z_X mal Eichelement; nu = 16 = 0 mod 16; "
                                     "N = 5 ist das kleinste N mit dieser Eigenschaft (Teil 1.2-1.4)"),
    "Hebung durch R^2 = (-1)^F gewählt": ("E", "nur die Spinorhebung erfüllt R^2 = (-1)^F modulo Eichung; die naive Hebung ist "
                                                "eichtrivial"),
    "Original -> Träger (Hardy-Residuen)": ("E", "äquivariant, Bild = R^2 = -1-Eigenraum (die 2)"),
    "die 3 aus dem Original": ("E (negativ, gescopt)", "keine fermionische Einteilchenquelle (Klassen alpha = 0, 1/2 fehlen); "
                                                       "modulo Eichung höchstens 4 von 5 Moden fermionisch; mit exakten "
                                                       "Trägercharakteren zerfällt die 16 statistisch 8/8"),
    "Spin-Statistik innen": ("E", "eta-Spinor, R^2 = (-1)^F auf allen Sektoren |Q| <= 10; Zuordnung {1, 3} eindeutig"),
    "4D-Helizität der Vertexoperatoren (Außenbild z79)": ("E (negativ, gescopt)", "äußeres 2 pi-Vorzeichen auf psi*_x|0> nicht "
                                                                                 "definit, Strombosonen mit Vorzeichen -1"),
    "Universalität der Hülle": ("E", "gleiche Stromalgebra, gleicher Zustand, gleiches T_osc; Tensorantwort nur über E(Q)"),
    "universelle Kopplung als Gravitation": ("O", "keine Quelle im KMS-Zustand, G offen, 2D-Kanal statt 4D-Spin 2"),
    "Verklebung innen/außen (Spin x_Z2 Z4)": ("O", "eine Lorentzwirkung mit -1 = R^2, die Vertexoperatoren kovariant lässt, "
                                                   "fehlt"),
}
for k, (s, t) in STATUS.items():
    say(f"   [{s}] {k}: {t}")
ZAHLEN["status"] = STATUS
ZAHLEN["S1_kriterien"] = {
    "Abbruch 1 (R ist nicht die Z4 der 16)": "nicht ausgelöst auf dem Träger (Klasse modulo G_SM = z_X); im Original teilweise: "
                                            "R erreicht nur die 2",
    "Abbruch 2 (keine konsistente 4D-Helizität der Vertexoperatoren)": "im Außenbild von z79 ausgelöst (gescopt); innen nur "
                                                                       "Vorzeichen und J_z, keine Drehgruppe",
    "Erfolg 'Materie = Ladungsgitter x Spinor'": "teilweise: Spinorvorzeichen innen exakt, Helizität außen nicht",
    "Erfolg 'Wechselwirkung = Vertexoperatoren mit universeller Hülle'": "Hülle universell (Teil 3), Kopplung nur über E(Q) "
                                                                         "bzw. Ladung, ohne Quelle",
}
urteil = ("S1 hält teilweise. Innen hält sie: Die Vierteldrehung der Marken auf dem Träger, R = diag(1, i, -1, -i, 1) mit "
          "det R = -1, ist R = (-1_5) g mit g in G_SM, ihre Spinorhebung auf Lambda^gerade C^5 ist z_X = exp(i pi X/2) mal "
          "Eichelement, also modulo Eichung genau die Spin-Z4 der 16 mit nu = 16 = 0 mod 16 (N = 5 ist das kleinste N mit "
          "dieser Eigenschaft); die naive Hebung wäre eichtrivial, und die Regel R^2 = (-1)^F des Originals wählt die "
          "Spinorhebung. Spin-Statistik gilt innen exakt: Auf allen Sektoren |Q_a|, |Q_b| <= 10 ist das Produkt der "
          "eta-Spinoreigenwerte R und sein Quadrat (-1)^F, mit eindeutiger Zuordnung {1, 3}; J_z = (Q_a - Q_b)/2 ist genau auf "
          "Fermionen halbzahlig. Die Hülle ist universell: dieselbe Stromalgebra, derselbe Zustand und dasselbe T_osc in jedem "
          "Sektor, und die Tensorantwort hängt nur über die innere Energie E(Q) mit quarterunabhängigem Koeffizienten ab. Nicht "
          "gehalten hat zweierlei. Erstens erreicht das Original über die Hardy-Residuen nur die 2 (R^2 = -1); die 3 hat keine "
          "fermionische Einteilchenquelle (es fehlen die Klassen alpha = 0, 1/2; modulo Eichung bleibt stets mindestens eine "
          "Trägermode R^2-gerade), und mit exakten Trägercharakteren zerfiele die 16 statistisch 8/8, R ist im Original "
          "also nur als Klasse modulo G_SM die Z4 der 16. Zweitens tragen die Vertexoperatoren im Außenbild von z79 keine "
          "konsistente Helizität: die äußere 2 pi-Drehung ist auf psi*_x|0> nicht definit (<U> = 0 ab n = 2) und gibt "
          "Strombosonen das Vorzeichen -1; innen gibt es nur Z4 und U(1), die volle eta-SU(2) ist im KMS-Zustand nicht "
          "implementierbar. Das Abbruchkriterium 'keine konsistente 4D-Helizität' ist damit für die wörtliche Fassung "
          "(Spin aus der äußeren Lorentzgruppe) erfüllt; offen bleibt eine Verklebung Spin x_Z2 Z4, die R^2 mit der äußeren "
          "2 pi-Drehung identifiziert und die Vertexoperatoren kovariant lässt. Universelle Kopplung im Sinn der Gravitation "
          "(Quelle, G) folgt nicht. Nichts ist befördert; 3+1-Aussagen bleiben Kinematik bzw. Lesart.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z82_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
