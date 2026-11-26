"""Z159: Die geeichte Schrittkette. Ist 'Transferschicht mal e^(-i (pi/2) H_el)' das Destri-de-Vega-Gitter aus z121,
und ist die Naht unter ihr stationär?

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart, N numerisch, O offen.

Frage (Welle C nach H6h50): z155 zeigt, dass der elektrische Term der geeichten Uhr nur als Schritt der Komposition
gelesen werden kann (P_schritt), nicht als Zustandsdynamik. z121 zeigt, dass die alternierende Kette aus Transfer und
Gatter das Destri-de-Vega-Lichtkegelgitter bei Theta = 1/2 ist. Hier: (1) Ist die Schrittkette mit H_el aus dem
Gauß-Gesetz auf den Schnittlinks exakt der Floquet-Operator von z121, und in welcher Linkgeometrie? (2) Ist das Produkt
nativer Nahtzustände (KMS bei beta = 4 pi, Leiter pi (4Z + 1)) unter einem Schritt stationär? (3) Ist der Grundzustand
der z152-Kette unter dem Schritt am Schnitt stationär? (4) Zeichnet Stationarität den Winkel pi/2 aus? (5) Bilanz.

Bestand: z120/z121 (Gatter U = B^+ e^(-i phi H_Q) B, Kette F = (prod_gerade U)(prod_ungerade U), XXX-Familie
R(u) = (1 + i u P)/(1 + i u) mit u = tan(phi/2), Destri-de-Vega bei Theta = 1/2; Transfer = Nullmode J_0, modendiagonal),
z125/z130 (P_univ: Paaruhr i^(Q_A - Q_B) = Delta^(-i/2) des Paars; P_TP), z133 (ungerade Gesamtladung: Umverteilung ab
Transferspin 3/2), z144 (P_uhr, P_Q, Leck nach (Q+2, Q-2) außer bei phi in (pi/2) Z; Treppenordnung vertauscht nicht
mit den Transfermatrizen), z145 (zirkulare Moden), z149 (N ungerade: einkörperige Tauschwurzel),
z152 (Naht = halbe Kette, t = 1, m = -+2), z155 (Gauß-Gesetz, P_eich, P_schnitt, P_schritt, F0 = -+1/4).

Modell der Nahtkette: L Nähte, je M Moden mit Modulniveaus eps_k = 4 pi x_k, x_k in Z + 1/4; Ladung Q_s gegen den
Dirac-See (x < 0 besetzt). Transfer auf der Bindung (s, s+1): J_x = (1/2) sum_k (c+_(s,k) c_(s+1,k) + h.c.), Strahlteiler
B = e^(-i (pi/2) J_x) (Vierteldrehung, P_sym). Elektrische Schicht: E = e^(-i phi H_el).

Literatur [L] (aus dem Gedächtnis zitiert, hier nicht nachgeschlagen; alle benutzten Identitäten sind unten selbst
geprüft, die Massenformel ist als [L] markiert): Destri/de Vega, Nucl. Phys. B 290 (1987) 363 (Lichtkegelgitter,
massives Thirring-Modell); Destri/de Vega, J. Phys. A 22 (1989) 1329 (integrable Feldtheorien aus dem Lichtkegelgitter);
Faddeev/Reshetikhin, Ann. Phys. 167 (1986) 227 (SU(2)-Fall); Vanicat/Zadnik/Prosen, Phys. Rev. Lett. 121 (2018) 030606
(integrable Trotterisierung); Rigol/Dunjko/Yurovsky/Olshanii, Phys. Rev. Lett. 98 (2007) 050405 (verallgemeinertes
Gibbs-Ensemble).

Teil A  Identifikation: Stern- und Kettengeometrie, Schichtform, Kontrollen.
Teil B  Stationarität der Naht: gleiche Nähte, gespiegelte Nähte, Tomita-Paar, Modulfluss, stationäre Zustände.
Teil C  Grundzustand der z152-Kette gegen den Schritt am Schnitt.
Teil D  Der Winkel.
Teil E  Bilanz und Urteil.
"""
import json
import os
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # deterministische Rundung (einfädig)

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
TOL = 1e-9
PI = float(np.pi)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val, text):
    BEFUNDE[name] = val
    say(f"   [Befund: {val}] {text}")


def r(x, n=5):
    return round(float(x), n)


def dag(a):
    return a.conj().T


def nrm(a):
    return float(np.linalg.norm(a))


def phasenabstand(A, B):
    """min_chi ||A - e^(i chi) B|| und die Phase."""
    ph = np.vdot(B.ravel(), A.ravel())
    ph = ph / abs(ph) if abs(ph) > 1e-14 else 1.0
    return nrm(A - ph * B), complex(ph)


def spurabstand(A):
    """(1/2) ||A||_1 für hermitesches A."""
    return 0.5 * float(np.sum(np.abs(np.linalg.eigvalsh((A + dag(A)) / 2))))


def spurgleich(A, B, nmax=6):
    """Spektralvergleich über die Spuren der Potenzen (unabhängig von der Sortierung auf dem Kreis)."""
    return max(abs(np.trace(np.linalg.matrix_power(A, n)) - np.trace(np.linalg.matrix_power(B, n))) for n in range(1, nmax + 1))


# ---------------------------------------------------------------------------------------------------------------------
# Werkzeuge: Fockraum (ganz oder Sektor) mit Jordan-Wigner-Vorzeichen; Nahtkette; XXX-Schaltkreis
class Fock:
    def __init__(self, n, states=None):
        self.n = n
        self.S = list(range(1 << n)) if states is None else sorted(states)
        self.idx = {b: i for i, b in enumerate(self.S)}
        self.d = len(self.S)
        self.occ = np.array([[(b >> i) & 1 for i in range(n)] for b in self.S], float)

    def eink(self, h):
        """Matrix von sum_ij h_ij c+_i c_j auf dem (unter Teilchenzahl abgeschlossenen) Zustandsraum."""
        M = np.zeros((self.d, self.d), complex)
        for i in range(self.n):
            for j in range(self.n):
                if abs(h[i, j]) < 1e-15:
                    continue
                lo, hi = min(i, j), max(i, j)
                mask = ((1 << hi) - 1) & ~((1 << (lo + 1)) - 1)
                for a, b in enumerate(self.S):
                    if not (b >> j) & 1 or (i != j and (b >> i) & 1):
                        continue
                    sg = -1.0 if bin(b & mask).count("1") % 2 else 1.0
                    M[self.idx[(b ^ (1 << j)) | (1 << i)], a] += h[i, j] * sg
        return M


class Nahtkette:
    """L Nähte mit je M Moden; Mode (s, k) hat den Index s M + k. vak = Besetzung des Dirac-Sees je Mode."""

    def __init__(self, L, M, ring, vak, states=None, versatz=0):
        self.L, self.M, self.ring = L, M, ring
        self.F = Fock(L * M, states)
        self.vak = np.array(vak, float)
        self.n = [[self.F.occ[:, s * M + k] for k in range(M)] for s in range(L)]
        self.Q = [sum(self.n[s][k] - self.vak[k] for k in range(M)) for s in range(L)]
        self.bonds = [(x, (x + 1) % L) for x in range(L if ring else L - 1)]
        self.Jx = {}
        for a, b in self.bonds:
            h = np.zeros((L * M, L * M))
            for k in range(M - versatz):
                h[a * M + k + versatz, b * M + k] = h[b * M + k, a * M + k + versatz] = 0.5
            self.Jx[(a, b)] = self.F.eink(h)

    def B(self, par, theta=PI / 2):
        """Transferschicht: Strahlteiler auf allen Bindungen mit linkem Platz der Parität par."""
        return expm(-1j * theta * sum(self.Jx[bd] for bd in self.bonds if bd[0] % 2 == par))

    def E_stern(self, phi, F0=0.0):
        return np.exp(-1j * phi * sum((F0 + q) ** 2 for q in self.Q) / 2)

    def E_kette(self, tau, F0=0.0):
        """Gauß-Gesetz der offenen Kette: Link l zwischen den Nähten l-1 und l trägt F0 + Ladung links davon."""
        acc, cum = np.zeros(self.F.d), np.zeros(self.F.d) + F0
        for l in range(1, self.L):
            cum = cum + self.Q[l - 1]
            acc = acc + cum ** 2 / 2
        return np.exp(-1j * tau * acc)

    def floquet(self, E):
        Be, Bo = self.B(0), self.B(1)
        return dag(Be) @ (E[:, None] * Be) @ dag(Bo) @ (E[:, None] * Bo)

    def neutral(self):
        """Zustände mit genau einem Fermion je Naht (M = 2) und ihr Spinindex (Bit s = besetzte Mode der Naht s)."""
        M = self.M
        ids = [i for i in range(self.F.d)
               if all(sum(self.F.occ[i, s * M + k] for k in range(M)) == 1 for s in range(self.L))]
        return ids, [sum(int(self.F.occ[i, M * s + 1]) << s for s in range(self.L)) for i in ids]

    def verschiebung(self):
        """Fermionische Translation um eine Naht, T c_(s,k) T^+ = c_(s+1,k)."""
        n, M, d = self.F.n, self.M, self.F.d
        T = np.zeros((d, d))
        for a, b in enumerate(self.F.S):
            img = [(i + M) % n for i in range(n) if (b >> i) & 1]
            inv = sum(1 for x in range(len(img)) for y in range(x + 1, len(img)) if img[x] > img[y])
            T[self.F.idx[sum(1 << i for i in img)], a] = (-1) ** inv
        return T


def swap_idx(L, x, y):
    idx = np.arange(2 ** L)
    bx, by = (idx >> x) & 1, (idx >> y) & 1
    return np.where(bx != by, idx ^ ((1 << x) | (1 << y)), idx)


def R_platz(L, x, y, u):
    P = np.zeros((2 ** L, 2 ** L))
    P[swap_idx(L, x, y), np.arange(2 ** L)] = 1
    return (np.eye(2 ** L) + 1j * u * P) / (1 + 1j * u)


def F_spin(L, u, ring=True, treppe=False):
    """z121: F = (prod_gerade R(u)) (prod_ungerade R(u)) auf C^2^L; treppe: R_01 R_12 R_23 ..."""
    bonds = [(x, (x + 1) % L) for x in range(L if ring else L - 1)]
    if treppe:
        out = np.eye(2 ** L, dtype=complex)
        for a, b in bonds:
            out = out @ R_platz(L, a, b, u)
        return out
    Ue, Uo = np.eye(2 ** L, dtype=complex), np.eye(2 ** L, dtype=complex)
    for a, b in bonds:
        if a % 2 == 0:
            Ue = R_platz(L, a, b, u) @ Ue
        else:
            Uo = R_platz(L, a, b, u) @ Uo
    return Ue @ Uo


def leck(Fm, ids):
    """Frobenius-Norm des Anteils, der den Unterraum ids verlässt."""
    rest = np.setdiff1d(np.arange(Fm.shape[0]), ids)
    return nrm(Fm[np.ix_(rest, ids)])


GITTER = [k * PI / 8 for k in range(17)] + [0.7, 2.0]              # Winkelraster: k pi/8 und zwei irrationale Kontrollen


def nullstellen(f, tol=TOL):
    """Winkel des Rasters (in Einheiten von pi), an denen f verschwindet."""
    return [r(g / PI, 4) for g in GITTER if f(g) < tol]


ALLE = [r(g / PI, 4) for g in GITTER]
VIERTEL = [0.0, 0.5, 1.0, 1.5, 2.0]
HALB = [0.0, 1.0, 2.0]
VOLL = [0.0, 2.0]

# ---------------------------------------------------------------------------------------------------------------------
say("Z159: Die geeichte Schrittkette gegen das Destri-de-Vega-Gitter und die Stationarität der Naht")
say("=" * 100)
say("\nTeil A  Identifikation, exakt")

# A1: zwei Nähte, M = 3
K2 = Nahtkette(2, 3, False, vak=[1, 0, 0])
NA, NB = sum(K2.n[0]), sum(K2.n[1])
hJ = np.zeros((6, 6), complex)
for k in range(3):
    hJ[k, 3 + k], hJ[3 + k, k] = -0.5j, 0.5j
Jx2, Jy2, Jz2 = K2.Jx[(0, 1)], K2.F.eink(hJ), np.diag((NA - NB) / 2).astype(complex)
B2 = expm(-1j * PI / 2 * Jx2)
a1 = max(nrm(dag(B2) @ np.diag(np.exp(-1j * ph * (NA ** 2 + NB ** 2) / 2)) @ B2
             - np.diag(np.exp(-1j * ph * (NA + NB) ** 2 / 4)) @ expm(-1j * ph * Jy2 @ Jy2)) for ph in (PI / 2, 0.7))
a1 += nrm(Jx2 @ Jy2 - Jy2 @ Jx2 - 1j * Jz2) + nrm(dag(B2) @ Jz2 @ B2 - Jy2)
check("A1 [E] Zwei Nähte (je 3 Moden): B^+ e^(-i phi (Q_A^2 + Q_B^2)/2) B = e^(-i phi N^2/4) e^(-i phi J_y^2) mit "
      "B^+ J_z B = J_y (z145, z155 A4 reproduziert; phi = pi/2 und 0,7)", a1 < TOL)

# A2: ein gemeinsamer Link
a2 = 0.0
for f in (0.0, 1.0, -0.25, 0.37):
    for tau in (PI, 0.9):
        for bg in (0, 1):
            X = np.diag(f + (NA + NB - 2 * bg) / 2) + Jy2
            a2 = max(a2, nrm(dag(B2) @ np.diag(np.exp(-1j * tau * (f + NA - bg) ** 2 / 2)) @ B2 - expm(-1j * tau * X @ X / 2)))
KP = Nahtkette(2, 2, False, vak=[1, 0])
idP, spP = KP.neutral()
BP = expm(-1j * PI / 2 * KP.Jx[(0, 1)])
link = {}
for f in (0.0, 1.0, 2.0, -0.25, 0.25, 0.5):
    G = dag(BP) @ np.diag(np.exp(-1j * PI * (f + KP.Q[0]) ** 2 / 2)) @ BP
    blk = G[np.ix_(idP, idP)]
    link[f] = (leck(G, idP), phasenabstand(blk, F_spin(2, 1.0, ring=False)[np.ix_(spP, spP)])[0],
               phasenabstand(blk, F_spin(2, -1.0, ring=False)[np.ix_(spP, spP)])[0])
ZAHLEN["A2 gemeinsamer Link, tau = pi: mittleres Feld f -> [Leck, Abstand zu R(1), Abstand zu R(-1)]"] = {
    str(f): [r(x) for x in v] for f, v in link.items()}
check("A2 [E] Ein gemeinsamer Link zwischen zwei Nähten (Feld f + q_A, Paarladung n): "
      "B^+ e^(-i tau (f + q_A)^2/2) B = e^(-i tau (f + n/2 + J_y)^2/2). Das Gatter ist die elektrische Phase des "
      "Linkfeldes im zirkularen Rahmen; maßgeblich ist das mittlere Feld f + n/2 des Links", a2 < TOL)
check("A2 [E] Neutrales Paar (ein Fermion je Naht), tau = pi: mittleres Feld 0 oder 2 gibt R(1) (das z120-Gatter), "
      "mittleres Feld 1 gibt R(-1) = R(1)^(-1), jeweils ohne Leck. Nichtganzes Feld leckt aus dem Sektor: "
      f"-1/4 mit Norm {r(link[-0.25][0], 4)} (das Singulett verlässt ihn mit Wahrscheinlichkeit "
      f"{r(link[-0.25][0] ** 2, 4)}), +1/4 ebenso, 1/2 mit Norm {r(link[0.5][0], 4)}",
      link[0.0][0] < TOL and link[0.0][1] < TOL and link[2.0][1] < TOL and link[1.0][0] < TOL and link[1.0][2] < TOL
      and link[-0.25][0] > 0.1 and link[0.25][0] > 0.1 and link[0.5][0] > 0.1)

# A3: Sterngeometrie, L = 4 Ring, ganzer Fockraum
K4 = Nahtkette(4, 2, True, vak=[1, 0])
id4, sp4 = K4.neutral()
Nk4 = [sum(K4.n[s][k] for s in range(4)) for k in range(2)]
Ntot4 = Nk4[0] + Nk4[1]
F4 = K4.floquet(K4.E_stern(PI / 2))
a3 = []
for F0 in (-0.25, 0.37):
    Fx = K4.floquet(K4.E_stern(PI / 2, F0))
    D = np.diag(Fx @ dag(F4))
    a3.append(nrm(Fx - D[:, None] * F4) + max(float(np.ptp(np.angle(D[Ntot4 == v] / D[Ntot4 == v][0]))) for v in range(9)))
d43, ph43 = phasenabstand(F4[np.ix_(id4, id4)], F_spin(4, 1.0)[np.ix_(sp4, sp4)])
F4b = K4.floquet(K4.E_stern(0.9))
d43b, _ = phasenabstand(F4b[np.ix_(id4, id4)], F_spin(4, float(np.tan(0.45)))[np.ix_(sp4, sp4)])
check("A3 [E] Sterngeometrie (jede Naht ein eigener Schnittlink mit Feld F0 + Q_s, H_el = sum_s (F0 + Q_s)^2/2), "
      "4 Nähte im Ring, ganzer Fockraum (256): Der Schrittoperator F = B_g^+ E B_g B_u^+ E B_u lässt den neutralen Sektor "
      "(ein Fermion je Naht) invariant und ist dort gleich dem Floquet-Operator von z121 mit R(1), ohne Basiswechsel "
      "und mit Phase 1; für phi = 0,9 gleich dem mit R(tan(phi/2))",
      leck(F4, id4) < TOL and d43 < TOL and abs(ph43 - 1) < TOL and d43b < TOL and leck(F4b, id4) < TOL)
check("A3 [E] Das Hintergrundfeld ist in der Sterngeometrie zentral: F(F0) = D F(0) mit einer Phase D, die nur von der "
      "Gesamtladung abhängt (F0 = -1/4 und 0,37). Der native Versatz -+1/4 aus z155 stört den Stern nicht",
      max(a3) < TOL)

# A4: L = 6 Ring, 12 Fermionmoden, alle Flavour-Sektoren mit sechs Fermionen
F6s = F_spin(6, 1.0)
a4, dims6 = [], []
for N0 in range(7):
    st = []
    for b in range(1 << 12):
        o = [(b >> i) & 1 for i in range(12)]
        if sum(o[0::2]) == N0 and sum(o[1::2]) == 6 - N0:
            st.append(b)
    K6 = Nahtkette(6, 2, True, vak=[1, 0], states=st)
    ids, sp = K6.neutral()
    F6 = K6.floquet(K6.E_stern(PI / 2))
    dd, ph = phasenabstand(F6[np.ix_(ids, ids)], F6s[np.ix_(sp, sp)])
    a4.append(max(dd, abs(ph - 1), leck(F6, ids)))
    dims6.append(K6.F.d)
ZAHLEN["A4 Sektorgrößen (6 Nähte, 12 Moden)"] = dims6
check("A4 [E] 6 Nähte im Ring (12 Fermionmoden, sieben Flavour-Sektoren bis Dimension 400): Im neutralen Sektor ist "
      "der Schrittoperator der Sterngeometrie exakt der 64-dimensionale Floquet-Operator von z121 (Phase 1, kein "
      "Leck)", max(a4) < TOL and sum(dims6) == 924)

# A5: Schichtform und Verschiebung
Be4, Bo4, E4 = K4.B(0), K4.B(1), np.diag(K4.E_stern(PI / 2))
W4, T4 = Be4 @ dag(Bo4), K4.verschiebung()
V4 = W4 @ E4 @ dag(T4)
a5 = (nrm(Bo4 @ F4 @ dag(Bo4) - dag(W4) @ E4 @ W4 @ E4) + nrm(T4 @ Be4 @ dag(T4) - Bo4) + nrm(T4 @ W4 @ dag(T4) - dag(W4))
      + nrm(dag(W4) @ E4 @ W4 @ E4 - T4 @ V4 @ V4 @ T4) + nrm(T4 @ T4 @ V4 - V4 @ T4 @ T4)
      + nrm(K4.floquet(K4.E_stern(0.0)) - np.eye(256)))
check("A5 [E] Schichtform: B_u F B_u^+ = W^+ E W E mit der freien Transferschicht W = B_g B_u^+ (zwei Lagen "
      "Strahlteiler, quasifrei). Mit der Verschiebung T um eine Naht ist T W T^+ = W^+, also "
      "W^+ E W E = T (W E T^+)^2 T: Der Floquet-Operator ist bis auf Konjugation T^2 mal das Quadrat des Einzelschritts "
      "'Transferschicht mal elektrische Schicht mal Verschiebung'. Ohne elektrische Schicht ist F = 1: Die Kette hat "
      "keine freie Dynamik", a5 < TOL)

# A6: Kontrollen
F1 = Be4 @ E4 @ Bo4 @ E4
G4 = []
for a, b in K4.bonds:
    Bq = expm(-1j * PI / 2 * K4.Jx[(a, b)])
    G4.append(dag(Bq) @ np.diag(np.exp(-1j * PI / 2 * (K4.Q[a] ** 2 + K4.Q[b] ** 2) / 2)) @ Bq)
F_tr = G4[0] @ G4[1] @ G4[2] @ G4[3]
a6 = (nrm(G4[0] @ G4[2] @ G4[1] @ G4[3] - F4), leck(F1, id4) ** 2 / 16, spurgleich(F1, F4),
      spurgleich(F_tr[np.ix_(id4, id4)], F4[np.ix_(id4, id4)]),
      phasenabstand(F_tr[np.ix_(id4, id4)], F_spin(4, 1.0, treppe=True)[np.ix_(sp4, sp4)])[0])
ZAHLEN["A6 Kontrollen: [F aus Gattern, Leckwahrscheinlichkeit einer einzelnen Transferschicht, Spurdifferenz dazu, "
       "Spurdifferenz Treppe gegen Ziegelmauer, Treppe gegen Spin-Treppe]"] = [r(x) for x in a6]
check("A6 [E] Kontrollen: (a) Mit nur einer Strahlteilerlage je Halbschritt, B_g E B_u E, verlässt die Kette den "
      f"neutralen Sektor (mittlere Leckwahrscheinlichkeit {r(a6[1], 4)}) und hat ein anderes Spektrum; nötig ist das "
      "Sandwich B^+ E B. (b) Die Treppenordnung U_01 U_12 U_23 U_30 hat ein anderes Spektrum als die Ziegelmauer "
      "(wie z144: sie vertauscht nicht mit den Transfermatrizen). Die Ziegelmauer ist eine Wahl (P_ziegel)",
      a6[0] < TOL and a6[1] > 0.1 and a6[2] > 0.1 and a6[3] > 0.1 and a6[4] < TOL)

# A7: Kettengeometrie mit dem Gauß-Gesetz der offenen Kette
K4o = Nahtkette(4, 2, False, vak=[1, 0])
id4o, sp4o = K4o.neutral()
kette = {}
for tau, F0 in ((PI, 0.0), (PI, 1.0), (PI, 2.0), (PI, -0.25), (PI, 0.5), (PI / 2, 0.0)):
    Fm = K4o.floquet(K4o.E_kette(tau, F0))
    blk = Fm[np.ix_(id4o, id4o)]
    kette[(r(tau / PI, 2), F0)] = [leck(Fm, id4o)] + [
        phasenabstand(blk, F_spin(4, u, ring=False)[np.ix_(sp4o, sp4o)])[0] for u in (1.0, -1.0, float(np.tan(PI / 8)))]
F_st_o = K4o.floquet(K4o.E_stern(PI / 2))
voll = phasenabstand(K4o.floquet(K4o.E_kette(PI, 0.0)), F_st_o)[0]
ZAHLEN["A7 Kettengeometrie (tau/pi, F0) -> [Leck, Abstand zu R(1), R(-1), R(tan(pi/8))]"] = {
    str(k): [r(x) for x in v] for k, v in kette.items()}
ZAHLEN["A7 Abstand Kette (tau = pi, F0 = 0) gegen Stern (phi = pi/2) auf dem ganzen Fockraum"] = r(voll)
check("A7 [E] Kettengeometrie (ein Link zwischen Nachbarnähten, jeder Link sieht die Gesamtladung links von sich; "
      "4 Nähte offen): Im neutralen Sektor ist jeder Link Schnittlink genau eines Paars, die übrigen Links sind "
      "Zuschauer mit konstanter Phase. Mit tau = g^2 T = pi und F0 = 0 oder 2 ist der Schrittoperator exakt der "
      "offene z121-Schaltkreis; F0 = 1 gibt R(-1) auf allen Bindungen (die zeitgespiegelte Kette)",
      kette[(1.0, 0.0)][0] < TOL and kette[(1.0, 0.0)][1] < TOL and kette[(1.0, 2.0)][1] < TOL
      and kette[(1.0, 1.0)][0] < TOL and kette[(1.0, 1.0)][2] < TOL)
check("A7 [E] Grenzen der Kettengeometrie: (a) Mit dem nativen Versatz F0 = -1/4 (z155) leckt die Kette aus dem "
      f"neutralen Sektor (Norm {r(kette[(1.0, -0.25)][0], 4)}), ebenso mit F0 = 1/2. (b) Bei tau = pi/2 (dem Wert "
      "des Sterns) entsteht R(tan(pi/8)), nicht R(1): Ein gemeinsamer Link zählt halb, P_viertel heißt hier "
      f"g^2 T = pi. (c) Außerhalb des neutralen Sektors ist die Kette nicht der Stern (Abstand {r(voll, 3)} auf dem "
      "Fockraum): Das Gatter trägt dort die Potenz P_sym^(f + n/2) mit der Ladung links vom Paar, eine Gauß-Schnur",
      kette[(1.0, -0.25)][0] > 0.1 and kette[(1.0, 0.5)][0] > 0.1 and kette[(0.5, 0.0)][3] < TOL
      and kette[(0.5, 0.0)][1] > 0.1 and voll > 0.1)

# A8: kleine Winkel, positive elektrische Energie
HQ4 = np.diag(sum(q ** 2 for q in K4.Q) / 2).astype(complex)
Heff = dag(Be4) @ HQ4 @ Be4 + dag(Bo4) @ HQ4 @ Bo4
dl = 1e-4
a8a = nrm((K4.floquet(K4.E_stern(dl)) - np.eye(256)) / (-1j * dl) - Heff)
PS4 = sum((np.eye(16) - np.eye(16)[swap_idx(4, a, b)]) / 2 for a, b in K4.bonds)
a8b = nrm(Heff[np.ix_(id4, id4)] - PS4[np.ix_(sp4, sp4)])
vak4 = np.zeros(256)
vak4[K4.F.idx[sum(1 << (2 * s) for s in range(4))]] = 1.0
a8c = nrm(Heff @ vak4) + abs(min(np.linalg.eigvalsh(Heff)))
check("A8 [E] Kleine Winkel: F(phi) = 1 - i phi H_eff + O(phi^2) mit H_eff = B_g^+ H_Q B_g + B_u^+ H_Q B_u >= 0 "
      f"(Rest {r(a8a, 4)} bei phi = 1e-4). Im neutralen Sektor ist H_eff = sum_Bindungen P_Singulett, die "
      "ferromagnetische Heisenberg-Kette; das Produkt der Dirac-Seen ist ein frustrationsfreier Grundzustand mit "
      "Energie 0", a8a < 1e-2 and a8b < TOL and a8c < TOL)


def t_matrix(L, u, xi):
    """Inhomogene Transfermatrix t(u) = Tr_0 R_(0,L-1)(u - xi_(L-1)) ... R_(0,0)(u - xi_0) mit Yangs R(w) = w + i P."""
    D = 2 ** L
    idx = np.arange(2 * D)
    a, sres = idx >> L, idx & (D - 1)
    M = np.eye(2 * D, dtype=complex)
    for j in range(L):
        bj = (sres >> j) & 1
        P = np.zeros((2 * D, 2 * D))
        P[np.where(a != bj, (bj << L) | (sres ^ (1 << j)), idx), idx] = 1
        M = ((u - xi[j]) * np.eye(2 * D) + 1j * P) @ M
    return M[:D, :D] + M[D:, D:]


a9 = []
for L in (4, 6):
    xi = np.array([0.5 if j % 2 == 0 else -0.5 for j in range(L)])
    Fs = F_spin(L, 1.0)
    X = t_matrix(L, -0.5, xi) @ np.linalg.inv(t_matrix(L, 0.5, xi))
    t1, t2 = t_matrix(L, 0.3 + 0.2j, xi), t_matrix(L, -1.1 + 0.7j, xi)
    a9.append(max(phasenabstand(Fs, X * nrm(Fs) / nrm(X))[0], nrm(t1 @ t2 - t2 @ t1) / nrm(t1 @ t2), nrm(t1 @ Fs - Fs @ t1) / nrm(t1)))
check("A9 [E] Gegenprobe zu z121 (L = 4, 6): Der Spin-Schaltkreis, dem die Schrittkette gleicht, ist das "
      "Lichtkegelgitter: F = t(-1/2) t(1/2)^(-1) bis auf Phase mit Yangs R(w) = w + i P und alternierenden "
      "Inhomogenitäten +-1/2 (Theta = 1/2), und [t(u), t(v)] = [t(u), F] = 0 für komplexe u, v", max(a9) < TOL)
befund("A Identifikation", "ja, in zwei Geometrien",
       "Teil A: Die geeichte Schrittkette ist exakt der Floquet-Operator von z121, ohne Basiswechsel. (S) Stern: jede "
       "Naht hat einen eigenen Schnittlink zu einem gemeinsamen Partner, H_el = sum_s (F0 + Q_s)^2/2; das ist das "
       "Gauß-Gesetz auf einem Sterngraphen (der Transfer ändert die Felder zweier Speichen um -+1). Dann gilt die "
       "Gleichheit im neutralen Sektor für jedes F0 mit g^2 T = pi/2, und die Fock-Fortsetzung ist die von z121/z144. "
       "(K) Kette: ein Link je Nachbarpaar mit dem Gauß-Gesetz der 1+1-dimensionalen Kette. Dann gilt sie nur im "
       "neutralen Sektor, nur für ganzes gerades mittleres Feld (ungerade: Zeitumkehr) und mit g^2 T = pi; der native "
       "Versatz -+1/4 ist dort unverträglich. In beiden Fällen ist das Sandwich B^+ E B und die Ziegelmauer nötig.")

# ---------------------------------------------------------------------------------------------------------------------
say("\nTeil B  Stationarität der Naht")
EPS2 = np.array([-3 * PI, PI])                                    # Moden x = -3/4 (See), +1/4


def produktzustand(K, leitern):
    w = np.zeros(K.F.d)
    for s in range(K.L):
        for k in range(K.M):
            w = w + leitern[s][k] * K.n[s][k]
    p = np.exp(-(w - w.min()))
    return p / p.sum(), w


p4, Kt4 = produktzustand(K4, [EPS2] * 4)
hS = np.zeros((8, 8))
for s in range(4):
    hS[2 * s, 2 * s + 1] = 1.0
Sp4 = K4.F.eink(hS)                                               # Flavour-Aufsteiger sum_s c+_(s,0) c_(s,1)
b1 = 0.0
for ph in GITTER:
    Fp = K4.floquet(K4.E_stern(ph))
    b1 = max(b1, max(nrm(Fp @ np.diag(x) - np.diag(x) @ Fp) for x in Nk4), nrm(Fp @ Sp4 - Sp4 @ Fp))
b1 += max(nrm(Be4 @ np.diag(x) - np.diag(x) @ Be4) for x in Nk4) + nrm(Be4 @ Sp4 - Sp4 @ Be4)
b1 += nrm(Kt4 - (EPS2[0] * Nk4[0] + EPS2[1] * Nk4[1]))
rho4 = np.diag(p4).astype(complex)
b1s = max(spurabstand(U @ rho4 @ dag(U) - rho4) for U in (F4, Be4, Bo4, E4, G4[0]))
check("B1 [E] Gleiche Nähte, modendiagonaler Transfer: Der Schrittoperator vertauscht für jedes phi des Rasters mit "
      "den Flavour-Ladungen N_k = sum_s n_(s,k) und mit der Flavour-SU(2); schon die Strahlteilerlage allein tut das. "
      "Der Modulgenerator des Produktzustands ist K = sum_k eps_k N_k, ein Cartan-Element der Flavour-Algebra", b1 < TOL)
check("B1 [E] Das Produkt nativer Nahtzustände (beta = 4 pi, Niveaus -3 pi und pi) ist deshalb exakt invariant unter "
      "dem Schritt, unter der Transferschicht allein, unter der elektrischen Schicht allein und unter dem "
      "Zweinahtgatter allein (Spurabstand < 1e-12, 4 Nähte im Ring)", b1s < 1e-12)

# B2: reduzierte Naht nach einem Schritt
r4 = (F4 @ rho4 @ dag(F4)).reshape(64, 4, 64, 4)
red = np.einsum("aibj->ij", r4).real
pr = np.diag(red)                                                 # Index = n_0 + 2 n_1 der Naht 0
niv = (float(np.log(pr[1] / pr[3]) / PI), float(np.log(pr[1] / pr[0]) / PI))
ZAHLEN["B2 Niveaus der Naht nach einem Schritt / pi (Teilchen, Loch)"] = [r(niv[0], 9), r(niv[1], 9)]
ZAHLEN["B2 Ladungsverteilung der Naht (Q = 0, +1, -1)"] = [r(pr[1], 8), r(pr[3], 8), r(pr[0], 8)]
check("B2 [E] Die reduzierte Dichtematrix einer Naht ist nach einem Schritt e^(-K)/Z mit dem nativen K: "
      f"Niveaus {r(niv[0], 9)} pi (Teilchen) und {r(niv[1], 9)} pi (Loch), diagonal",
      abs(niv[0] - 1) < 1e-9 and abs(niv[1] - 3) < 1e-9 and nrm(red - np.diag(pr)) < 1e-12)

# B3: Kontrollen (zwei Nähte, vier Moden)
XS = np.array([-7 / 4, -3 / 4, 1 / 4, 5 / 4])
EPS4 = 4 * PI * XS
KM = Nahtkette(2, 4, False, vak=(XS < 0).astype(float))
QA, QB = KM.Q
KA = sum(EPS4[k] * KM.n[0][k] for k in range(4))
KB = sum(EPS4[k] * KM.n[1][k] for k in range(4))
BM = expm(-1j * PI / 2 * KM.Jx[(0, 1)])


def UM(phi):
    return dag(BM) @ (np.exp(-1j * phi * (QA ** 2 + QB ** 2) / 2)[:, None] * BM)


def UMk(phi, kette_):
    Bq = expm(-1j * PI / 2 * kette_.Jx[(0, 1)])
    return dag(Bq) @ (np.exp(-1j * phi * (kette_.Q[0] ** 2 + kette_.Q[1] ** 2) / 2)[:, None] * Bq)


def gibbs(w):
    p = np.exp(-(w - w.min()))
    return np.diag(p / p.sum()).astype(complex)


rhoM = gibbs(KA + KB)
KV = Nahtkette(2, 4, False, vak=(XS < 0).astype(float), versatz=1)
b3 = (spurabstand(UM(PI / 2) @ rhoM @ dag(UM(PI / 2)) - rhoM),
      spurabstand(UMk(PI / 2, KV) @ rhoM @ dag(UMk(PI / 2, KV)) - rhoM),
      spurabstand(UM(PI / 2) @ gibbs(KA + 0.5 * KB) @ dag(UM(PI / 2)) - gibbs(KA + 0.5 * KB)))
ZAHLEN["B3 Spurabstand nach einem Gatter: [gleiche Nähte, verschobener Transfer J_1, Leiter der zweiten Naht halbiert]"] = [
    r(x, 6) for x in b3]
check("B3 [E] Kontrollen (zwei Nähte, Moden x = -7/4 .. 5/4): Mit dem um eine Stufe verschobenen Transfer J_1 ist das "
      f"Produkt nicht invariant (Spurabstand {r(b3[1], 4)}), ebenso mit ungleichen Leitern ({r(b3[2], 4)}). Nötig sind "
      "gleiche Leitern und der modendiagonale Transfer J_0 (P_gleich)", b3[0] < 1e-12 and b3[1] > 1e-3 and b3[2] > 1e-3)

# B4: unverschränktes Produkt Naht mal gespiegelte Naht
gesp = {}
for name, epsb in (("Schmidt-Paarung (eps <-> -eps)", np.array([3 * PI, -PI])),
                   ("Nachbarniveaus (-3 pi <-> -pi, pi <-> 3 pi)", np.array([-PI, 3 * PI]))):
    pg, Kg = produktzustand(K4, [EPS2, epsb, EPS2, epsb])
    rg = np.diag(pg).astype(complex)
    r1 = F4 @ rg @ dag(F4)
    rb = Be4 @ rg @ dag(Be4)
    occ, rr = [], rg.copy()
    for _ in range(200):
        rr = F4 @ rr @ dag(F4)
        occ.append([float(np.real(np.trace(rr * K4.n[0][k][None, :]))) for k in range(2)])
    occ = np.array(occ)
    n_a = [float(pg @ K4.n[0][k]) for k in range(2)]
    n_b = [float(pg @ K4.n[1][k]) for k in range(2)]
    scan = [spurabstand(Fp @ rg @ dag(Fp) - rg) for Fp in (K4.floquet(K4.E_stern(g)) for g in GITTER)]
    gesp[name] = dict(spur_F=spurabstand(r1 - rg), relent_F=float(np.real(np.trace((r1 - rg) * Kg[None, :]))),
                      spur_B=spurabstand(rb - rg), relent_B=float(np.real(np.trace((rb - rg) * Kg[None, :]))),
                      n_vor=n_a, n_partner=n_b, n_nach=occ[0].tolist(), n_mittel=occ.mean(0).tolist(),
                      n_gge=[(x + y) / 2 for x, y in zip(n_a, n_b)],
                      null=[ALLE[i] for i, x in enumerate(scan) if x < TOL])
    say(f"   {name}: Schritt: Spurabstand {r(gesp[name]['spur_F'], 5)}, Relativentropie {r(gesp[name]['relent_F'], 4)}; "
        f"Transferschicht allein: {r(gesp[name]['spur_B'], 5)}, {r(gesp[name]['relent_B'], 4)}; Besetzung Naht 0 "
        f"vorher {[r(x, 5) for x in n_a]}, nach 1 Schritt {[r(x, 5) for x in occ[0]]}, Mittel über 200 Schritte "
        f"{[r(x, 4) for x in occ.mean(0)]}, Flavour-Mittel {[r(x, 5) for x in gesp[name]['n_gge']]}")
ZAHLEN["B4 gespiegelte Nähte, unverschränkt"] = {k: {kk: ([r(x, 6) for x in vv] if isinstance(vv, list) else r(vv, 6))
                                                    for kk, vv in v.items()} for k, v in gesp.items()}
g1, g2 = gesp["Schmidt-Paarung (eps <-> -eps)"], gesp["Nachbarniveaus (-3 pi <-> -pi, pi <-> 3 pi)"]
check("B4 [E/N] Naht neben gespiegelter Naht (Leiter pi (4Z + 3)), unverschränkt, 4 Plätze im Ring: nicht stationär. "
      f"Mit der Schmidt-Paarung der Moden ändert ein Schritt den Zustand um Spurabstand {r(g1['spur_F'], 4)} "
      f"(Relativentropie {r(g1['relent_F'], 3)}); die Besetzungen der Naht liegen danach bei "
      f"{r(g1['n_nach'][0], 3)} und {r(g1['n_nach'][1], 3)}, im Zeitmittel bei {r(g1['n_mittel'][0], 3)} und "
      f"{r(g1['n_mittel'][1], 3)}, nahe dem Flavour-Mittel 1/2 (Niveau 0, Leiter weg). "
      f"Mit der Paarung benachbarter Niveaus sind es {r(g2['spur_F'], 5)} und {r(g2['relent_F'], 4)}. Stationär ist "
      "das Produkt im Raster nur für phi in 2 pi Z; die elektrische Schicht allein ändert nichts",
      g1["spur_F"] > 0.5 and g2["spur_F"] > 1e-3 and g1["null"] == VOLL and g2["null"] == VOLL
      and abs(g1["n_gge"][0] - 0.5) < 1e-9 and abs(g1["n_mittel"][0] - 0.5) < 0.1)

# B5: Automorphismus und Modulfluss (zwei Nähte, vier Moden)
def komm(U, d):
    return nrm(U * d[None, :] - d[:, None] * U)


zeiten = (0.13, 0.5, 1.7)
flu = {
    "gemeinsamer Modulfluss e^(i t (K_A + K_B))": nullstellen(lambda g: max(komm(UM(g), np.exp(1j * t * (KA + KB))) for t in zeiten)),
    "Paaruhr i^(Q_A - Q_B) = Delta^(-i/2) des Paars": nullstellen(lambda g: komm(UM(g), 1j ** (QA - QB))),
    "Nahtparität (-1)^Q_A": nullstellen(lambda g: komm(UM(g), (-1.0) ** QA)),
    "Nahtuhr R_A = i^Q_A": nullstellen(lambda g: komm(UM(g), 1j ** QA)),
    "Modulfluss einer Naht e^(i t K_A)": nullstellen(lambda g: max(komm(UM(g), np.exp(1j * t * KA)) for t in zeiten)),
}
Ug = UM(PI / 2)
ts = [k / 8 for k in range(17)] + [0.3]
t_null = [t for t in ts if komm(Ug, np.exp(1j * t * (KA - KB))) < TOL]
uhr_q = np.exp(0.5j * KA) / 1j ** QA
uhr = nrm(np.exp(0.5j * (KA - KB)) - 1j ** (QA - QB)) + float(np.max(np.abs(uhr_q - uhr_q[0])))
ZAHLEN["B5 Winkel (in pi), bei denen das Gatter vertauscht"] = flu
ZAHLEN["B5 Zeiten t, bei denen das Viertelgatter mit dem relativen Modulfluss e^(i t (K_A - K_B)) vertauscht"] = t_null
check("B5 [E] Der Schritt ist ein Automorphismus, der für jedes phi mit dem gemeinsamen Modulfluss, mit der Paaruhr "
      "(P_univ in der Form von z130) und mit der Nahtparität vertauscht. Mit der Uhr R_A einer einzelnen Naht und mit "
      "ihrem Modulfluss vertauscht das Gatter im Raster nur für phi in 2 pi Z, wo es trivial ist",
      all(flu[k] == ALLE for k in list(flu)[:3]) and flu["Nahtuhr R_A = i^Q_A"] == VOLL
      and flu["Modulfluss einer Naht e^(i t K_A)"] == VOLL)
check("B5 [E] Der relative Modulfluss e^(i t (K_A - K_B)) vertauscht mit dem Viertelgatter genau zu den Zeiten "
      "t in Z/2 (Raster t = k/8 und 0,3); bei t = 1/2 ist er die Paaruhr, und e^(i K_A/2) ist bis auf eine Konstante "
      "die Nahtuhr. Der Schritt vertauscht also mit der Uhrgruppe des Paars, nicht mit dem Fluss",
      t_null == [0.0, 0.5, 1.0, 1.5, 2.0] and uhr < 1e-8)

# B6: Tomita-Paar (Thermofeldform): Naht A und Partner mit den gepaarten Moden
nk = 1 / (np.exp(EPS4) + 1)
Om = np.zeros(256)
Om[0] = 1.0
for k in range(4):
    Om_neu = np.sqrt(1 - nk[k]) * Om
    for a, b in enumerate(KM.F.S):                                # at+_k a+_k auf Zustände ohne die Moden k und 4 + k
        if abs(Om[a]) > 0 and not (b >> k) & 1 and not (b >> (4 + k)) & 1:
            sg = (-1) ** (bin(b & ((1 << k) - 1)).count("1") + bin((b | (1 << k)) & ((1 << (4 + k)) - 1)).count("1"))
            Om_neu[KM.F.idx[b | (1 << k) | (1 << (4 + k))]] += np.sqrt(nk[k]) * sg * Om[a]
    Om = Om_neu


def rho_A(v):
    m = v.reshape(16, 16)                                         # Zeile = Moden des Partners, Spalte = Moden von A
    return m.T @ m.conj()


rA0 = rho_A(Om)
pQ = {q: float(np.sum(Om[QA == q] ** 2)) for q in (-2, -1, 0, 1, 2)}
par_soll = float(np.prod(np.tanh(np.abs(EPS4) / 2)))
b6 = (abs(Om @ Om - 1) + nrm(KM.Jx[(0, 1)] @ Om) + nrm(QA * Om - QB * Om)
      + max(nrm(UM(g) @ Om - np.exp(-1j * g * QA ** 2) * Om) for g in GITTER))
tom_eig = nullstellen(lambda g: 1 - abs(np.vdot(Om, UM(g) @ Om)))
tom_red = nullstellen(lambda g: nrm(rho_A(UM(g) @ Om) - rA0))
fid_tom = [abs(np.vdot(Om, UM(g) @ Om)) ** 2 for g in (PI / 2, PI)]
fid_formel = [abs(sum(p * np.exp(-1j * g * q * q) for q, p in pQ.items())) ** 2 for g in (PI / 2, PI)]
ZAHLEN["B6 Tomita-Paar: Ladungsverteilung p_Q (Q = -1, 0, 1, 2)"] = [r(pQ[q], 8) for q in (-1, 0, 1, 2)]
ZAHLEN["B6 Tomita-Paar: <(-1)^Q> = prod tanh(|eps_k|/2)"] = r(par_soll, 6)
ZAHLEN["B6 Tomita-Paar: Treue |<Omega|U Omega>|^2 bei phi = pi/2, pi"] = [r(x, 6) for x in fid_tom]
check("B6 [E] Tomita-Paar (Naht und Partner in Thermofeldform, vier Modenpaare): Der Tomita-Vektor Omega ist ein "
      "Transfer-Singulett (J_x Omega = 0, Q_A = Q_B auf Omega). Das Gatter wirkt auf ihm als die lokale Ladungsphase, "
      "U(phi) Omega = e^(-i phi Q_A^2) Omega: Der Transfer fällt heraus, die reduzierte Naht bleibt für jedes phi "
      "exakt e^(-K)", b6 < TOL and tom_red == ALLE)
check(f"B6 [E] Omega selbst ist im Raster nur für phi in 2 pi Z Eigenvektor. Treue bei phi = pi/2: {r(fid_tom[0], 5)} "
      f"= |sum_Q p_Q e^(-i pi Q^2/2)|^2, bei phi = pi: {r(fid_tom[1], 5)} = <(-1)^Q>^2 mit <(-1)^Q> = "
      f"prod_k tanh(|eps_k|/2) = {r(par_soll, 4)} (z130 hat mit zwei Moden auf dem Niveau pi tanh(pi/2)^2 = "
      f"{r(np.tanh(PI / 2) ** 2, 4)}). Der Schritt dreht die Purifikation, nicht die Naht",
      tom_eig == VOLL and max(abs(a - b) for a, b in zip(fid_tom, fid_formel)) < 1e-10
      and abs(np.sqrt(fid_tom[1]) - par_soll) < 1e-10 and abs(np.tanh(PI / 2) ** 2 - 0.841168) < 1e-6)

# B7: Zentralisator
sq = np.sqrt(p4)
fid_zentr = abs(np.sum(p4 * np.diag(F4))) ** 2
b7 = nrm(F4 @ (sq[:, None] * dag(F4)) - np.diag(sq))
ZAHLEN["B7 Treue der kanonischen Purifikation unter dem Schritt nur auf den Nähten (4 Nähte)"] = r(fid_zentr, 6)
check("B7 [E] Zentralisator: F liegt im Zentralisator des Produktzustands, F rho^(1/2) F^+ = rho^(1/2). Die "
      "kanonische Purifikation (jede Naht mit ihrem eigenen Partner) ist deshalb unter dem doppelten Schritt F mal "
      "J F J invariant; unter dem Schritt nur auf den Nähten bleibt die Naht, nicht die Purifikation "
      f"(Treue |tr(rho F)|^2 = {r(fid_zentr, 5)})", b7 < 1e-12 and 0.5 < fid_zentr < 1 - 1e-6)

# B8: ladungsscharfe Nähte (P_Q)
Nt, Jz2m = QA + QB, ((QA - QB) / 2) ** 2
ger = np.where(np.mod(Nt, 2) == 0)[0]
ung = np.where(np.mod(Nt, 2) == 1)[0]
pq_ger = nullstellen(lambda g: nrm((UM(g) * Jz2m[None, :] - Jz2m[:, None] * UM(g))[np.ix_(ger, ger)]))
pq_ung = nullstellen(lambda g: nrm((UM(g) * Jz2m[None, :] - Jz2m[:, None] * UM(g))[np.ix_(ung, ung)]))
klasse = np.array([min(a, b) * 10 + max(a, b) for a, b in zip(QA, QB)])
lw = -(KA + KB)
lw = lw - lw.max()
Zm = float(np.sum(np.exp(lw)))


def leckwahrscheinlichkeit(phi):
    """Wahrscheinlichkeit im nativen Produktzustand, dass ein Gatter das ungeordnete Ladungspaar ändert."""
    A = np.abs(UM(phi)) ** 2
    weg = np.array([float(np.sum(A[klasse != klasse[i], i])) for i in range(256)])
    return float(np.sum(np.exp(lw) * weg) / Zm)


lp = {str(r(g / PI, 4)): leckwahrscheinlichkeit(g) for g in (PI / 2, 0.7, PI / 4)}
ZAHLEN["B8 Leckfreiheit [U, J_z^2] = 0: Winkel (in pi) für gerade / ungerade Paarladung"] = [pq_ger, pq_ung]
ZAHLEN["B8 log10 der Umverteilungswahrscheinlichkeit im nativen Zustand (phi/pi)"] = {k: r(np.log10(v), 3) for k, v in lp.items()}
check("B8 [E] Ladungsscharfe Nähte (P_Q aus z144: Nahtladungen werden höchstens getauscht): Als Operatoraussage gilt "
      "[U(phi), (Q_A - Q_B)^2] = 0 bei gerader Paarladung genau für phi in (pi/2) Z, bei ungerader nur für phi in "
      "pi Z (Raster; Transferspin bis 2). Das Viertelgatter verteilt bei ungerader Paarladung um (wie z133, z149). Erst mit "
      "Niveauabgleich (gerade Paarladung) macht Stationarität ladungsscharfer Nähte den Winkel diskret: das ist "
      "P_Q <=> P_uhr aus z144, keine neue Herleitung", pq_ger == VIERTEL and pq_ung == HALB)
check("B8 [N] Im nativen Zustand ist diese Bedingung unsichtbar: Die Wahrscheinlichkeit, dass ein Gatter das "
      f"Ladungspaar umverteilt, ist 10^{r(np.log10(lp['0.5']), 2)} bei phi = pi/2 (ungerade Paare) und "
      f"10^{r(np.log10(lp[str(r(0.7 / PI, 4))]), 2)} bei phi = 0,7; der thermische Zustand ist fast ganz Dirac-See",
      all(v < 1e-8 for v in lp.values()) and lp["0.5"] > 0)

# B9: stationäre Zustände der integrablen Kette
kommutant = {}
for L in (4, 6, 8):
    ang = np.sort(np.mod(np.angle(np.linalg.eigvals(F_spin(L, 1.0))), 2 * PI))
    cl = [1]
    for g in np.diff(ang):
        if g < 1e-7:
            cl[-1] += 1
        else:
            cl.append(1)
    if ang[0] + 2 * PI - ang[-1] < 1e-7 and len(cl) > 1:
        cl[0] += cl.pop()
    kommutant[str(L)] = [len(cl), int(sum(c * c for c in cl)), 4 ** L]
ZAHLEN["B9 Spin-Sektor: L -> [verschiedene Eigenwerte, Dimension des Kommutanten, 4^L]"] = kommutant
R1 = (np.eye(4) + 1j * np.eye(4)[[0, 2, 1, 3]]) / (1 + 1j)
ra, rb2 = np.diag([0.9, 0.1]).astype(complex), np.array([[0.6, 0.2], [0.2, 0.4]], complex)
hom = (nrm(R1 @ np.kron(ra, ra) @ dag(R1) - np.kron(ra, ra)), nrm(R1 @ np.kron(rb2, rb2) @ dag(R1) - np.kron(rb2, rb2)),
       nrm(R1 @ np.kron(ra, rb2) @ dag(R1) - np.kron(ra, rb2)))
p_neutral = float(np.sum(p4[id4]))
p_ref = float(p4[K4.F.idx[sum(1 << (2 * s) for s in range(4))]])
Sz = np.array([sum(0.5 - ((x >> s) & 1) for s in range(4)) for x in sp4])
b9 = nrm(p4[id4] / p4[id4].sum() - np.exp(4 * PI * Sz) / np.sum(np.exp(4 * PI * Sz)))
ZAHLEN["B9 nativer Produktzustand (4 Nähte): [Gewicht des neutralen Sektors, Gewicht des Referenzzustands, "
       "Magnon-Fugazität e^(-4 pi)]"] = [r(p_neutral, 6), r(p_ref, 6), float(f"{np.exp(-4 * PI):.4e}")]
check("B9 [E/N] Stationär unter dem integrablen F sind genau die Funktionen seiner Erhaltungsgrößen. Der Kommutant "
      f"hat im Spin-Sektor Dimension {kommutant['4'][1]}, {kommutant['6'][1]}, {kommutant['8'][1]} für L = 4, 6, 8 "
      "(von 4^L). Unter den Produktzuständen lässt das Gatter R(1) genau die homogenen fest (rho_a x rho_b nur für "
      "rho_a = rho_b; Stichprobe)", hom[0] < TOL and hom[1] < TOL and hom[2] > 1e-2
      and kommutant["4"][1] < 4 ** 4 and kommutant["8"][1] < 4 ** 8)
check("B9 [E] Der native Zustand gehört dazu, aber nur über die ultralokalen Flavour-Ladungen: Im neutralen Sektor "
      "ist er e^(4 pi S^z)/Z, der unendlich heiße Zustand mit Magnetfeld 4 pi; alle Destri-de-Vega-Ladungen haben den "
      f"Multiplikator null. Er liegt zu {r(p_ref, 4)} auf dem Bethe-Referenzzustand (alle Nähte im Dirac-See; "
      f"neutraler Sektor {r(p_neutral, 4)}), also am ferromagnetischen Ende, nicht am antiferromagnetischen Vakuum "
      "der Destri-de-Vega-Lesart", b9 < 1e-9 and p_ref > 0.8 and abs(p_ref - pr[1] ** 4) < 1e-9)

# B10: Einmagnon-Dispersion über dem Referenzzustand
disp = 0.0
for L in (8, 12):
    for ph in (PI / 2, 0.9, PI / 4):
        u = float(np.tan(ph / 2))
        ca, cb = 1 / (1 + 1j * u), 1j * u / (1 + 1j * u)
        Ue, Uo = np.eye(L, dtype=complex), np.eye(L, dtype=complex)
        for x in range(L):
            Rm = np.eye(L, dtype=complex)
            y = (x + 1) % L
            Rm[x, x] = Rm[y, y] = ca
            Rm[x, y] = Rm[y, x] = cb
            if x % 2 == 0:
                Ue = Rm @ Ue
            else:
                Uo = Rm @ Uo
        ev = np.linalg.eigvals(Ue @ Uo)
        ks = 2 * PI * np.arange(L // 2) / (L // 2)
        eps_k = np.arccos(np.clip(np.cos(ph / 2) ** 2 - np.sin(ph / 2) ** 2 * np.cos(ks), -1, 1))
        soll = np.concatenate([np.exp(-1j * (ph - eps_k)), np.exp(-1j * (ph + eps_k))])
        disp = max(disp, max(float(np.min(np.abs(ev - s))) for s in soll), max(float(np.min(np.abs(soll - e))) for e in ev))
check("B10 [E/N] Über dem Referenzzustand hat ein Magnon die Quasienergie omega = phi -+ eps(k) mit "
      "cos eps = cos^2(phi/2) - sin^2(phi/2) cos k (k = Impuls der Zweinahtzelle; L = 8, 12; phi = pi/2, 0,9, pi/4). "
      "Bei k = 0 ist omega = (tan(phi/2)/4) k^2: quadratisch, nichtrelativistisch. Die lineare Berührung aus z121 "
      "liegt bei omega = phi (k = pi), Steigung sin(phi/2)", disp < 1e-8)
befund("B Stationarität", "ja für gleiche Nähte (für jedes phi), nein für Naht mit gespiegeltem Partner",
       "Teil B: Das Produkt gleicher nativer Nähte ist unter dem Schritt exakt stationär, die reduzierte Naht bleibt "
       "e^(-K) mit der Leiter pi (4Z + 1). Der Grund ist nicht die Integrabilität und nicht der Winkel: Der "
       "Modulgenerator ist ein Cartan-Element der Flavour-Symmetrie, mit der jede modendiagonale, ladungsdiagonale "
       "Komposition vertauscht. Der native Zustand hat Gibbs-Form in den Flavour-Ladungen allein, ist also der "
       "entropiemaximale stationäre Zustand bei gegebenen Flavour-Ladungen (unendlich heiß mit Feld 4 pi) und liegt "
       "fast ganz auf dem Bethe-Referenzzustand. Eine Naht neben ihrer "
       "gespiegelten Partnerin ist nicht stationär; im Tomita-Paar fällt der Transfer heraus, und der Schritt dreht "
       "nur die Purifikation.")

# ---------------------------------------------------------------------------------------------------------------------
say("\nTeil C  Der Grundzustand der Kette aus z152 gegen den Schritt am Schnitt")
NG, CUT = 12, 6
GK = Fock(NG, [b for b in range(1 << NG) if bin(b).count("1") == NG // 2])
h1 = np.zeros((NG, NG))
for n in range(NG - 1):
    h1[n, n + 1] = h1[n + 1, n] = -1.0
for n in range(NG):
    h1[n, n] = -2.0 if n % 2 == 0 else 2.0
HG = GK.eink(h1).real
wG, vG = np.linalg.eigh(HG)
psi0 = vG[:, 0] * np.sign(vG[np.argmax(np.abs(vG[:, 0])), 0])
la = np.array([b & ((1 << CUT) - 1) for b in GK.S])
ra_ = np.array([b >> CUT for b in GK.S])
nl = np.array([bin(x).count("1") for x in range(1 << CUT)])
QL = nl[la] - CUT / 2


def rho_L(ps):
    m = np.zeros((1 << CUT, 1 << CUT), complex)
    m[la, ra_] = ps
    return m @ dag(m)


def niveaus(ps):
    rr, e = rho_L(ps), {}
    for q in range(CUT + 1):
        ids = np.where(nl == q)[0]
        e[q] = -np.log(np.maximum(np.sort(np.linalg.eigvalsh(rr[np.ix_(ids, ids)]))[::-1], 1e-300))
    n0 = min(sorted(e), key=lambda k: e[k][0])
    x0 = e[n0][0]
    return dict(p=float((e[n0 + 1][0] - x0) / PI), h=float((e[n0 - 1][0] - x0) / PI), ph=float((e[n0][1] - x0) / PI),
                Z=float(np.exp(x0)))


nv0 = niveaus(psi0)
rL0 = rho_L(psi0)
pQL = {q: float(np.sum(psi0[QL == q] ** 2)) for q in (-3, -2, -1, 0, 1, 2, 3)}
ZAHLEN["C1 halbe Kette (12 Plätze): Niveaus / pi und Z"] = {k: r(v, 6) for k, v in nv0.items()}
ZAHLEN["C1 Ladungsverteilung der halben Kette (Q = -1, 0, +1)"] = [r(pQL[q], 8) for q in (-1, 0, 1)]
check(f"C1 [N] Grundzustand der freien Kette (12 Plätze, t = 1, m = -2, halbe Füllung, Lücke {r(wG[1] - wG[0], 4)}): "
      f"Die halbe Kette ist die Naht, Niveaus {r(nv0['p'], 5)} pi, {r(nv0['h'], 5)} pi, {r(nv0['ph'], 5)} pi für "
      f"Teilchen, Loch, Paar und Z = {r(nv0['Z'], 6)} (z152, z155 D1). Die Ladungsverteilung stimmt mit der des "
      "Tomita-Paars aus B6 auf 1e-5",
      abs(nv0["p"] - 1) < 2e-3 and abs(nv0["h"] - 3) < 2e-3 and abs(nv0["ph"] - 4) < 2e-3 and abs(nv0["Z"] - 1.043298) < 1e-5
      and max(abs(pQL[q] - pQ[q]) for q in (-1, 0, 1)) < 1e-5)


def E_schnitt(tau, F0=0.0):
    return np.exp(-1j * tau * (F0 + QL) ** 2 / 2)


TAUS = [k * PI / 4 for k in range(17)]
c2 = max(spurabstand(rho_L(E_schnitt(t, F0) * psi0) - rL0) for t in TAUS for F0 in (0.0, -0.25))
fid_E = {F0: abs(np.vdot(psi0, E_schnitt(PI, F0) * psi0)) ** 2 for F0 in (0.0, -0.25)}
dE_E = float(np.real(np.vdot(E_schnitt(PI) * psi0, HG @ (E_schnitt(PI) * psi0))) - wG[0])
c2f = abs(fid_E[0.0] - abs(sum(p * np.exp(-1j * PI * q * q / 2) for q, p in pQL.items())) ** 2)
kommHE = nrm(HG * E_schnitt(PI)[None, :] - E_schnitt(PI)[:, None] * HG) / nrm(HG)
ZAHLEN["C2 elektrischer Schritt allein (tau = pi): Treue für F0 = 0, -1/4; Energiezuwachs; ||[H, E]||/||H||"] = [
    r(fid_E[0.0], 6), r(fid_E[-0.25], 6), r(dE_E, 5), r(kommHE, 5)]
check("C2 [E] Elektrischer Schritt allein, e^(-i tau (F0 + Q_L)^2/2) auf dem Schnittlink: Die reduzierte Halbkette "
      "bleibt für jedes tau des Rasters und beide F0 exakt die Naht (Spurabstand < 1e-12); der Schritt liegt im "
      "Zentralisator. Der Grundzustand ist aber kein Eigenzustand: Treue bei tau = pi (phi = pi/2) "
      f"{r(fid_E[0.0], 5)} = |sum_Q p_Q e^(-i pi Q^2/2)|^2 wie im Tomita-Paar, Energiezuwachs {r(dE_E, 4)} t, und "
      f"[H, E] ist nicht null (relativ {r(kommHE, 4)})",
      c2 < 1e-12 and c2f < 1e-10 and abs(fid_E[0.0] - fid_tom[0]) < 1e-4 and dE_E > 0.1 and kommHE > 1e-2)

hc = np.zeros((NG, NG))
hc[CUT - 1, CUT] = hc[CUT, CUT - 1] = 1.0
B_alpha = expm(-1j * PI / 4 * GK.eink(hc))                         # Hüpfen über den Schnittlink, Vierteldrehung
_, vv1 = np.linalg.eigh(h1)
Ckor = vv1[:, :NG // 2] @ vv1[:, :NG // 2].T
Us, ss, Vts = np.linalg.svd(Ckor[:CUT, CUT:])
hb = np.zeros((NG, NG))
for k in range(CUT):
    if ss[k] > 1e-7:
        uu, wv = np.zeros(NG), np.zeros(NG)
        uu[:CUT], wv[CUT:] = Us[:, k], Vts[k, :]
        hb += 0.5 * (np.outer(uu, wv) + np.outer(wv, uu))
B_beta = expm(-1j * PI / 2 * GK.eink(hb))                          # modendiagonaler Transfer zwischen Schmidt-Partnern
hby = np.zeros((NG, NG), complex)
for k in range(CUT):
    if ss[k] > 1e-7:
        uu, wv = np.zeros(NG), np.zeros(NG)
        uu[:CUT], wv[CUT:] = Us[:, k], Vts[k, :]
        hby += -0.5j * np.outer(uu, wv) + 0.5j * np.outer(wv, uu)
Py_beta = expm(1j * PI * GK.eink(hby))                             # Tausch der Schmidt-Partner
ps_beta = dag(B_beta) @ (np.exp(-1j * PI * QL ** 2 / 2) * (B_beta @ psi0))
katze = min(nrm(ps_beta - c * ((1 - 1j) / 2 * psi0 + (1 + 1j) / 2 * sgn * (Py_beta @ psi0))) for c in (1, -1, 1j, -1j)
            for sgn in (1, -1, 1j, -1j))
ZAHLEN["C4 Abstand von U_beta|GS> zur Überlagerung ((1 - i)|GS> + (1 + i) z Tausch|GS>)/2, z^4 = 1; Überlapp "
       "<GS|Tausch|GS>"] = [r(katze, 6), r(abs(np.vdot(psi0, Py_beta @ psi0)), 8)]
schmidt_niv = [float(np.log((1 + np.sqrt(1 - 4 * x * x)) / (1 - np.sqrt(1 - 4 * x * x))) / PI) for x in ss[:3]]
trans = {}
for name, Bc in (("alpha: Hüpfen über den Schnittlink", B_alpha), ("beta: Schmidt-Moden", B_beta)):
    scan = []
    for t in TAUS:
        ps = dag(Bc) @ (E_schnitt(t) * (Bc @ psi0))
        scan.append((abs(np.vdot(psi0, ps)) ** 2, spurabstand(rho_L(ps) - rL0)))
    Uc = dag(Bc) @ (E_schnitt(PI)[:, None] * Bc)
    ps = Uc @ psi0
    trans[name] = dict(treue_B=abs(np.vdot(psi0, Bc @ psi0)) ** 2, spur_B=spurabstand(rho_L(Bc @ psi0) - rL0),
                       treue=scan[4][0], spur=scan[4][1], niveaus=niveaus(ps),
                       dE=float(np.real(np.vdot(ps, HG @ ps)) - wG[0]), komm=nrm(HG @ Uc - Uc @ HG) / nrm(HG),
                       null=[r(t / PI, 2) for t, s in zip(TAUS, scan) if s[1] < TOL],
                       min_spur=min(s[1] for t, s in zip(TAUS, scan) if 0 < t < 4 * PI - 1e-9))
    say(f"   {name}: Transfer allein: Treue {r(trans[name]['treue_B'], 5)}, Spurabstand der Halbkette "
        f"{r(trans[name]['spur_B'], 5)}; Schritt B^+ E(pi) B: Treue {r(trans[name]['treue'], 5)}, Spurabstand "
        f"{r(trans[name]['spur'], 5)}, Niveaus danach {[r(trans[name]['niveaus'][k], 4) for k in ('p', 'h', 'ph')]} pi, "
        f"Energiezuwachs {r(trans[name]['dE'], 4)} t, ||[H, U]||/||H|| = {r(trans[name]['komm'], 4)}")
ZAHLEN["C3 Schritt am Schnitt"] = {k: {kk: ({a: r(b, 5) for a, b in vv.items()} if isinstance(vv, dict) else
                                           (vv if isinstance(vv, list) else r(vv, 6))) for kk, vv in v.items()}
                                   for k, v in trans.items()}
ZAHLEN["C3 Schmidt-Niveaus der Halbkette / pi (erste drei)"] = [r(x, 4) for x in schmidt_niv]
ta, tb = trans["alpha: Hüpfen über den Schnittlink"], trans["beta: Schmidt-Moden"]
check("C3 [N] Schritt am Schnitt mit dem Transfer der Kette (Hüpfen über den Schnittlink, Vierteldrehung) und "
      f"E bei tau = pi, F0 = 0: Treue des Grundzustands {r(ta['treue'], 4)}, Spurabstand der Halbkette "
      f"{r(ta['spur'], 4)}, Niveaus danach {r(ta['niveaus']['p'], 3)} pi und {r(ta['niveaus']['h'], 3)} pi statt pi "
      f"und 3 pi, Energiezuwachs {r(ta['dE'], 3)} t bei Lücke {r(wG[1] - wG[0], 3)} t. Im Raster tau = k pi/4 bleibt "
      f"die Halbkette nur bei tau in 4 pi Z die Naht (kleinster Spurabstand sonst {r(ta['min_spur'], 4)})",
      ta["treue"] < 0.9 and ta["spur"] > 0.05 and abs(ta["niveaus"]["p"] - 1) > 0.1 and ta["null"] == [0.0, 4.0]
      and ta["min_spur"] > 1e-2)
check("C4 [N] Dasselbe mit dem modendiagonalen Transfer zwischen Schmidt-Partnern (im Thermofeldbild ist das "
      f"Paarerzeugung): Treue {r(tb['treue'], 4)}, Spurabstand {r(tb['spur'], 4)}; schon der Transfer allein lässt "
      f"nur Treue {r(tb['treue_B'], 4)}. Bei tau = pi ist der Zustand die gleichgewichtige Überlagerung aus "
      "Grundzustand und seinem Bild unter dem Tausch der Schmidt-Partner (Abstand < 1e-6, Bild orthogonal): Die "
      "Halbkette ist halb Naht, halb gespiegelte Naht. Auch hier bleibt die Naht im Raster nur bei tau in 4 pi Z",
      tb["treue"] < 0.9 and tb["spur"] > 0.05 and tb["null"] == [0.0, 4.0] and tb["treue_B"] < 0.5 and katze < 1e-6
      and abs(np.vdot(psi0, Py_beta @ psi0)) < 1e-9)

# C5: freier Transfer zwischen Nähten gegen den Schrittoperator
Hfrei = sum(2 * K4.Jx[bd] for bd in K4.bonds)
c5 = (nrm(F4 @ Hfrei - Hfrei @ F4) / nrm(Hfrei), komm(Hfrei, Kt4) / nrm(Hfrei), komm(F4, Kt4), ta["komm"], tb["komm"])
ZAHLEN["C5 [||[F, H_frei]||/||H_frei||, ||[H_frei, K]||, ||[F, K]||, ||[H, U_alpha]||/||H||, ||[H, U_beta]||/||H||]"] = [
    r(x, 5) for x in c5]
check("C5 [E] Vertauschen: Der freie Hamiltonoperator vertauscht nicht mit dem Schritt, weder in der Gitterlesart "
      f"(||[H, U]||/||H|| = {r(c5[3], 3)} bzw. {r(c5[4], 3)}) noch als freies Hüpfen zwischen gleichen Nähten "
      f"(||[F, H_frei]||/||H_frei|| = {r(c5[0], 3)}). Gemeinsam ist beiden nur die Flavour-Symmetrie: H_frei und F "
      "vertauschen beide mit dem Modulgenerator K = sum_k eps_k N_k. Das ist die genaue Form von 'Zustand und Gatter "
      "sind verträglich, obwohl sie zwei Enden sind'", c5[0] > 1e-2 and c5[1] < TOL and c5[2] < TOL and c5[3] > 1e-2)
befund("C Grundzustand", "nein (mit Transfer über den Schnitt); die Naht allein bleibt unter dem elektrischen Schritt",
       "Teil C: Der Grundzustand der z152-Kette ist kein Eigenzustand des Schritts am Schnitt. Der elektrische Schritt "
       "allein lässt die Naht (reduzierte Halbkette) für jedes tau exakt fest und dreht nur die Phasen zwischen den "
       "Ladungssektoren der Purifikation. Mit einem Transfer zwischen Naht und Partner davor und danach wird die Leiter "
       "zerstört, für jeden Winkel außer dem trivialen. Die Komposition mit P_schritt ist also eine Komposition "
       "gleicher Nähte untereinander, nicht der Naht mit ihrem Gitterpartner.")

# ---------------------------------------------------------------------------------------------------------------------
say("\nTeil D  Der Winkel")
nat_null = nullstellen(lambda g: spurabstand(K4.floquet(K4.E_stern(g)) @ rho4 @ dag(K4.floquet(K4.E_stern(g))) - rho4), tol=1e-12)
tabelle = [
    ("Produkt gleicher nativer Nähte stationär (B1)", nat_null),
    ("reduzierte Naht im Tomita-Paar invariant (B6)", tom_red),
    ("Schritt vertauscht mit gemeinsamem Modulfluss, Paaruhr, Nahtparität (B5)", flu["Paaruhr i^(Q_A - Q_B) = Delta^(-i/2) des Paars"]),
    ("Schritt vertauscht mit Uhr oder Modulfluss einer Naht (B5)", flu["Nahtuhr R_A = i^Q_A"]),
    ("Tomita-Vektor Eigenvektor (B6)", tom_eig),
    ("Naht neben gespiegelter Naht stationär (B4)", g1["null"]),
    ("Leckfreiheit P_Q, gerade Paarladung (B8)", pq_ger),
    ("Leckfreiheit P_Q, ungerade Paarladung (B8)", pq_ung),
]
for was, null in tabelle:
    txt = "alle Winkel des Rasters" if null == ALLE else f"phi/pi in {null}"
    say(f"   {was}: {txt}")
ZAHLEN["D1 Winkel (in pi), die das Kriterium erfüllen"] = {was: ("alle" if null == ALLE else null) for was, null in tabelle}
check("D1 [E] Kein Stationaritäts- oder Vertauschungskriterium zeichnet phi = pi/2 aus. Entweder gilt es für jeden "
      "Winkel (gleiche Nähte, reduzierte Naht, gemeinsamer Modulfluss, Paaruhr), oder nur für phi in 2 pi Z, wo das "
      "Gatter trivial ist (Uhr und Fluss einer Naht, Tomita-Vektor, gespiegelte Nähte). Diskret mit pi/2 wird der "
      "Winkel allein durch P_Q bei gerader Paarladung, und das ist P_uhr selbst",
      tabelle[0][1] == ALLE and tabelle[1][1] == ALLE and tabelle[2][1] == ALLE and tabelle[3][1] == VOLL
      and tabelle[4][1] == VOLL and tabelle[5][1] == VOLL and tabelle[6][1] == VIERTEL and tabelle[7][1] == HALB)

# D2: Wurzel des Nahttauschs
hy = np.zeros((8, 8), complex)
for k in range(4):
    hy[k, 4 + k], hy[4 + k, k] = -0.5j, 0.5j
JyM = KM.F.eink(hy)
Upi = UM(PI)
Py = expm(1j * PI * JyM)
d2 = (phasenabstand(Upi[np.ix_(ger, ger)] * np.exp(1j * PI * Nt[ger] ** 2 / 4)[None, :], Py[np.ix_(ger, ger)])[0],
      phasenabstand(Upi[np.ix_(ung, ung)] * np.exp(1j * PI * Nt[ung] ** 2 / 4)[None, :], np.eye(len(ung)))[0],
      nrm(Py @ np.diag(QA).astype(complex) @ dag(Py) - np.diag(QB)))
check("D2 [E] U(pi) ist bei gerader Paarladung der einkörperige Nahttausch P_sym = e^(i pi J_y) mal einer Phase der "
      "Paarladung, bei ungerader eine Phase. 'Eine Potenz des Schritts vertauscht die Nähte samt ihren Modulflüssen' "
      "gilt für U(phi)^n = U(n phi) bei phi = pi/n mit jedem n. P_viertel ist die Wahl n = 2 (Wurzel des Tauschs, "
      "z121); die Bedingung macht den Winkel rational, nicht fest", max(d2) < TOL)
befund("D Winkel", "nein",
       "Teil D: Stationarität wählt den Winkel nicht. Für die konsistente Lesart (gleiche Nähte) schränkt sie phi gar "
       "nicht ein; wo sie einschränkt, erzwingt sie das triviale Gatter. P_viertel folgt nicht aus P_univ: Paaruhr und "
       "Nahtparität vertauschen mit U(phi) für jedes phi. Diskret wird der Winkel nur über P_Q (gleichwertig P_uhr) "
       "oder über 'U^n ist der Nahttausch' mit der Wahl n = 2.")

# ---------------------------------------------------------------------------------------------------------------------
say("\nTeil E  Bilanz")
theta = {str(r(g / PI, 4)): (0.5 / np.tan(g / 2), 4 * np.exp(-PI * 0.5 / np.tan(g / 2))) for g in (PI / 2, PI / 4, PI / 8, PI / 16)}
ZAHLEN["E Destri-de-Vega [L]: phi/pi -> [Theta = cot(phi/2)/2, m a = 4 e^(-pi Theta)]"] = {
    k: [r(v[0], 5), float(f"{v[1]:.4e}")] for k, v in theta.items()}
check("E1 [E/L] Einordnung in das Lichtkegelgitter: Das Gatter bei Winkel phi ist R(u) mit u = tan(phi/2) (A3), in "
      "Yangs Form w + i P mit w = 1/u = 2 Theta, also Theta = cot(phi/2)/2 und Theta = 1/2 bei phi = pi/2 (z121). "
      "Nach der Literaturformel m a = 4 e^(-pi Theta) [L] ist dort m a = "
      f"{r(theta['0.5'][1], 4)}: keine Skalentrennung. Der massive Kontinuumslimes verlangt Theta -> unendlich, also "
      "phi -> 0, das dem Viertelwinkel entgegengesetzte Ende", abs(theta["0.5"][0] - 0.5) < 1e-12
      and abs(theta["0.5"][1] - 4 * np.exp(-PI / 2)) < 1e-12 and theta["0.0625"][1] < 1e-3)
bilanz = [
    ("E", "Identifikation", "Die geeichte Schrittkette ist exakt der Floquet-Operator von z121: in der Sterngeometrie "
     "(eigener Schnittlink je Naht, g^2 T = pi/2, jedes F0) im neutralen Sektor ohne Basiswechsel, in der "
     "Kettengeometrie (Gauß-Gesetz der 1+1-dimensionalen Kette, g^2 T = pi) nur im neutralen Sektor bei ganzem "
     "mittlerem Feld (gerade: z121, ungerade: Zeitumkehr). Mit dem nativen Versatz -+1/4 ist nur der Stern "
     "verträglich."),
    ("E", "Schichtform", "F ist bis auf Konjugation T^2 (W E T^+)^2 mit der quasifreien Transferschicht W = B_g B_u^+. "
     "Ohne elektrische Schicht ist F = 1; die Kette hat keine freie Dynamik. Für kleine phi ist sie die "
     "ferromagnetische Heisenberg-Kette mit H_eff >= 0."),
    ("E", "Stationarität", "Das Produkt gleicher nativer Nähte ist exakt stationär, die reduzierte Naht bleibt "
     "e^(-K) mit der Leiter pi (4Z + 1), für jedes phi. Grund: K ist ein Cartan-Element der Flavour-Symmetrie der "
     "Kette. P_schritt zerstört die Naht nicht."),
    ("E negativ", "Auszeichnung", "Die Stationarität benutzt weder die Integrabilität noch den Winkel und legt die "
     "Leiter nicht fest: Jedes homogene Produkt mit beliebigen Modenniveaus ist stationär. Der native Zustand ist "
     "unendlich heiß mit Feld 4 pi, fast ganz der Bethe-Referenzzustand."),
    ("E/N negativ", "Partnerlesart", "Naht neben gespiegelter Naht ist nicht stationär; der Grundzustand der "
     "z152-Kette ist kein Eigenzustand des Schritts am Schnitt, und mit Transfer über den Schnitt geht die Leiter "
     "verloren. Der elektrische Schritt allein lässt die Naht fest und dreht die Purifikation (Treue 0,920)."),
    ("E negativ", "Winkel", "Kein Stationaritätskriterium zeichnet pi/2 aus. P_viertel folgt nicht aus P_univ. "
     "Diskret wird der Winkel nur über P_Q (= P_uhr, bei gerader Paarladung) oder die Wahl n = 2 in "
     "'U^n = Nahttausch'."),
    ("C", "P_uhr", "P_null ist unter P_eich + P_schnitt ein Gauß-Gesetz (z155) und als Schritt mit der Naht "
     "verträglich (hier). P_viertel bleibt die Zahl g^2 T mal Zahl der bewegten Links = pi. P_uhr würde ersetzt durch "
     "P_eich + P_stern (oder P_kette0) + P_schritt + P_viertel; nötig bleiben P_sym, P_ord in der Form P_ziegel und "
     "P_gleich."),
    ("C", "P_univ, P_ord, P_sym", "P_univ gilt für den Schritt bei jedem Winkel und wählt nichts. P_ord wird zur "
     "Ziegelmauer (die Treppe ist ein anderer Operator). P_sym ist das Sandwich B^+ E B; eine einzelne Transferlage "
     "leckt."),
    ("E/C", "P_gleich gegen P_TP", "Die stationäre Lesart verlangt gleichchirale Nähte mit dem modendiagonalen Transfer "
     "J_0 (P_gleich). Im Tomita-Paar (P_TP), wo P_univ nach z130 nativ ist, ist der Tomita-Vektor ein "
     "Transfer-Singulett und das Gatter auf ihm nur die lokale Phase e^(-i phi Q^2). Nativität von P_univ und "
     "Wechselwirkung am Zustand liegen also in verschiedenen Lesarten; P_univ bleibt für den Kandidaten Prämisse."),
    ("L", "K15", "Die Kopplung bleibt g^2 T. Neu ist nur: Die Zahl hängt von der Linkgeometrie ab (Stern pi/2, Kette "
     "pi). Für kleine Kopplung ist die Schrittkette die Trotter-Form der Heisenberg-Kette."),
    ("L", "K16", "Bei festem Theta = 1/2 gibt es keinen massiven Kontinuumslimes (m a = 0,83 [L]); der Limes "
     "liegt bei phi -> 0. Um den antiferromagnetischen Zustand ist die Kette masselos (su(2)_1, z121 [N]). Der native "
     "Zustand liegt am anderen Ende: über dem Referenzzustand ist das Magnon quadratisch, omega = k^2/4 bei "
     "phi = pi/2."),
    ("O", "Nächster Test", "Die Linkgeometrie aus dem Bestand entscheiden. Der Stern verlangt einen gemeinsamen "
     "Partner aller Nähte, über den die Wilson-Linie des Transfers läuft; hat jede Naht ihren eigenen Partner "
     "(P_TP, P_CTM), so ist die Wilson-Linie ein gemeinsamer Link, g^2 T = pi, und der native Versatz -+1/4 lässt "
     "das Singulett mit Wahrscheinlichkeit 1/2 lecken. Rechnung: drei Nähte mit explizitem Linkfeld auf beiden "
     "Graphen, Transfer als eichinvarianter Operator, gegen den nativen Transfer J_0 (z116, z119) und den "
     "Viertelversatz (z113, z155)."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
urteil = ("Bedingt ja. (1) Identifikation [E]: Die geeichte Schrittkette ist exakt das Destri-de-Vega-Gitter aus z121, "
          "ohne Basiswechsel und mit Phase 1 im neutralen Sektor (4 und 6 Nähte, bis 12 Fermionmoden). Dafür sind zu "
          "wählen: das Sandwich B^+ E B (P_sym), die Ziegelmauer (P_ziegel, die Kettenform von P_ord) und eine "
          "Linkgeometrie. Im Stern (jede Naht ein eigener Schnittlink zu einem gemeinsamen Partner, P_stern) gilt die "
          "Gleichheit für jedes Hintergrundfeld mit g^2 T = pi/2. Mit dem Gauß-Gesetz einer 1+1-dimensionalen Kette "
          "(ein Link je Nachbarpaar, P_kette0) gilt sie nur im neutralen Sektor, nur bei ganzem mittlerem Feld "
          "(gerade: z121, ungerade: die zeitgespiegelte Kette) und mit g^2 T = pi; der native Versatz -+1/4 aus z155 "
          "erzeugt dort ein Leck (Singulett mit Wahrscheinlichkeit 1/2). (2) Stationarität [E]: Das "
          "Produkt gleicher nativer Nähte ist unter dem Schritt, unter der Transferschicht und unter dem Gatter exakt "
          "invariant, und die reduzierte Naht bleibt e^(-K) mit der Leiter pi (4Z + 1). P_schritt zerstört die Naht "
          "also nicht; 'P_eich + P_schnitt + P_schritt' ist eine konsistente Lesart, unter der weiteren Prämisse "
          "P_gleich (gleiche Leitern, modendiagonaler Transfer). (3) Der Grund ist trivial [E negativ]: Der "
          "Modulgenerator ist ein Cartan-Element der Flavour-Symmetrie, mit der jede solche Komposition vertauscht. "
          "Integrabilität und Winkel gehen nicht ein, die Leiter wird nicht festgelegt, und der native Zustand ist "
          "der unendlich heiße Zustand mit Feld 4 pi nahe am Bethe-Referenzzustand, nicht das Vakuum der "
          "Destri-de-Vega-Lesart. (4) Gitterlesart [E/N negativ]: Der Grundzustand der z152-Kette ist kein "
          "Eigenzustand des Schritts am Schnitt. Der elektrische Schritt allein lässt die Naht fest (Treue der "
          "Purifikation 0,920 bei phi = pi/2); mit Transfer zwischen Naht und Gitterpartner geht die Leiter verloren "
          "(Niveaus 0,32 pi und 2,25 pi statt pi und 3 pi). Die Komposition ist eine zwischen gleichen Nähten, nicht "
          "zwischen Naht und Partner. (5) Winkel [E negativ]: Stationarität und Vertauschen mit dem Modulfluss "
          "zeichnen pi/2 nicht aus; P_viertel folgt nicht aus P_univ. (6) Für den Kandidaten: P_uhr bleibt Prämisse, "
          "die neun Prämissen sind unverändert, keine entfällt. P_uhr wäre durch P_eich + P_stern + P_schritt und die "
          "Zahl g^2 T = pi/2 ersetzbar, also durch mehr. K13 ist nicht berührt. K15: die Kopplung bleibt offen und "
          "hängt von der Linkgeometrie ab. K16: bei Theta = 1/2 kein massiver Kontinuumslimes [L]. Nichts ist "
          "befördert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z159_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z159.txt").write_text("\n".join(OUT) + "\n")
