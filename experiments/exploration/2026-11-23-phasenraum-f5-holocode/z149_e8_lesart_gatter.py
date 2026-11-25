"""Z149: Trägt die E8-Lesart des z120-Gatters auf Operatorebene? (D5 und A3 als zirkulare Moden, Klebestrom als Transfer)

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart, N numerisch, O offen.

Frage: z145 Teil D hat auf Gitterebene vorgeschlagen, das z120-Gatter als 4 pi-Drehung des D5-Faktors der E8-Verklebung
zu lesen (D5 um 2 pi gegen A3 um -2 pi, Transfer = Klebestrom (16,4)). Hier wird das wörtlich als Operatoraussage geprüft:
mit den Nullmoden der Stufe-1-Darstellung von E8, der su(2) einer Klebewurzel gamma und B_gamma = e^(-i pi J_x/2).

Bestand: z120/z121/z133/z144 (U = B^+ e^(-i phi H_Q) B = e^(-i phi N^2/4) e^(-i phi J_y^2), P_uhr, P_sym, P_ord), z145
(Paarphase; Formen; E8-Ketten), z144 Teil A (zwei Lagrangesche), z130 (P_univ = Niveauabgleich, P_TP), z125 (Leckformel).

Teil A  Modell: Stufe-1-Darstellung von E8 auf dem Gitter, Niveau L_0 = 1 (240 Impulszustände + 8 Oszillatorzustände),
        Kozykel bilinear; die Nullmoden erfüllen E8 (Stichprobe von 9600 Kommutatoren). su(2) einer Klebewurzel.
Teil B  Die Drehungen als Operatoren: W = e^(4 pi i L_0^D5), T5 = e^(2 pi i L_0^D5), T3; Clifford-Form auf Klebeströmen.
Teil C  Wörtliche Prüfung: B_gamma^+ W B_gamma gegen das z120-Gatter, je Multiplett; höhere Spins in der Kette
        p + m gamma mit Oszillatoren längs gamma (exaktes su(2)_1-Teilmodell); die 2 pi-Drehung allein; zwei gleiche Nähte.
Teil D  Prämissen einzeln: P_uhr, P_null, P_univ, P_sym, P_ord, P_IF.
Teil E  Dynamik: erzeugte Gruppen (endlich oder dicht), Automorphismusfrage, Netzstruktur.
Teil F  Einordnung (zwei gleiche Nähte oder Partnerlesart P_TP), Bilanz und Urteil.
"""
import itertools
import json
import os
import time
from fractions import Fraction as Fr
from math import factorial, comb
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # deterministische Rundung (einfädig)

import numpy as np
import scipy.sparse as sps
from scipy.linalg import expm

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
TOL = 1e-9


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val, text):
    BEFUNDE[name] = val
    say(f"   [Befund: {val}] {text}")


def dag(a):
    return a.conj().T


def nrm(a):
    return float(np.linalg.norm(a))


def rd(x, n=6):
    return round(float(x), n) + 0.0


# ------------------------------------------------------------------------------------------------------------------
# Modellbau: E8-Gitter in verdoppelten Koordinaten (v = 2p ganzzahlig), Kozykel, Nullmoden auf dem Niveau L_0 = 1
# ------------------------------------------------------------------------------------------------------------------
def wurzeln(par):
    """240 Wurzeln von D8 plus Spinorklasse; par = 0: gerade Zahl von Minuszeichen (E8), par = 1: ungerade (E8')."""
    out = []
    for i, j in itertools.combinations(range(8), 2):
        for si, sj in itertools.product((2, -2), repeat=2):
            v = [0] * 8
            v[i], v[j] = si, sj
            out.append(tuple(v))
    for sg in itertools.product((1, -1), repeat=8):
        if sum(1 for s in sg if s < 0) % 2 == par:
            out.append(tuple(sg))
    return sorted(out)


def einfache(par):
    """Einfache Wurzeln (Bourbaki), verdoppelt; für par = 1 mit gespiegelter achter Koordinate."""
    s = [(1, -1, -1, -1, -1, -1, -1, 1), (2, 2, 0, 0, 0, 0, 0, 0)]
    for i in range(3, 9):
        v = [0] * 8
        v[i - 2], v[i - 3] = 2, -2
        s.append(tuple(v))
    if par:
        s = [tuple(v[:7]) + (-v[7],) for v in s]
    return s


def ip(a, b):
    """Skalarprodukt zweier verdoppelter Vektoren (ganzzahlig für Gittervektoren)."""
    s = sum(x * y for x, y in zip(a, b))
    assert s % 4 == 0
    return s // 4


class Modell:
    def __init__(self, par, kozykel=True):
        self.par = par
        self.kozykel = kozykel
        self.R = wurzeln(par)
        self.idx = {r: i for i, r in enumerate(self.R)}
        S = np.array(einfache(par), float) / 2
        Sinv = np.linalg.inv(S)
        self.cart = np.rint(S @ S.T).astype(int)
        self.koord = {}
        for r in self.R:
            c = np.rint(np.array(r, float) / 2 @ Sinv).astype(int)
            assert np.allclose(c @ S, np.array(r) / 2)
            self.koord[r] = c
        self.eta = np.array([[self.cart[i, j] if i > j else (1 if i == j else 0) for j in range(8)]
                             for i in range(8)], int)
        self.dim = 248
        self.H = [np.diag([r[i] / 2 for r in self.R] + [0.0] * 8) for i in range(8)]
        self._E = {}

    def eps(self, a, b):
        if not self.kozykel:
            return 1
        return -1 if int(self.koord[a] @ self.eta @ self.koord[b]) % 2 else 1

    def E(self, a):
        """Nullmode des Vertexoperators zur Wurzel a auf dem Niveau L_0 = 1 (reelle 248 x 248-Matrix)."""
        if a in self._E:
            return self._E[a]
        m = np.zeros((248, 248))
        for j, b in enumerate(self.R):
            s = ip(a, b)
            if s == -1:
                m[self.idx[tuple(x + y for x, y in zip(a, b))], j] = self.eps(a, b)
            elif s == -2:                              # b = -a: Koeffizient von z^-1 ist a . alpha_(-1) |0>
                for i in range(8):
                    m[240 + i, j] = self.eps(a, b) * a[i] / 2
        for i in range(8):                             # E_a alpha^i_(-1)|0> = -a^i |a>
            m[self.idx[a], 240 + i] = -a[i] / 2
        self._E[a] = m
        return m

    def cartan(self, w2):
        """w . H für einen verdoppelten Vektor w2 (w = w2/2)."""
        return sum(w2[i] / 2 * self.H[i] for i in range(8))

    def su2(self, g):
        Jp = self.E(g)
        Jm = Jp.T
        return (Jp + Jm) / 2, (Jp - Jm) / 2j, self.cartan(g) / 2


def fpoly(J, f, knoten=(-1.0, -0.5, 0.0, 0.5, 1.0), pw=None):
    """f(J) für hermitesches J mit Eigenwerten in den Knoten, als Polynom (ohne Diagonalisierung)."""
    kn = np.array(knoten)
    c = np.linalg.solve(np.vander(kn, len(kn), increasing=True).astype(complex),
                        np.array([f(x) for x in kn], complex))
    if pw is None:
        pw = potenzen(J, len(kn))
    return sum(c[n] * pw[n] for n in range(len(kn)))


def potenzen(J, n):
    pw = [np.eye(J.shape[0], dtype=complex), J.astype(complex)]
    for _ in range(n - 2):
        pw.append(pw[-1] @ pw[1])
    return pw


I2 = np.eye(2, dtype=complex)
Z2 = np.diag([1, -1]).astype(complex)


def kron(*ms):
    out = np.eye(1, dtype=complex)
    for m in ms:
        out = np.kron(out, m)
    return out


def fock(n):
    low = np.array([[0, 1], [0, 0]], complex)
    return [kron(*([Z2] * k + [low] + [I2] * (n - k - 1))) for k in range(n)]


# ------------------------------------------------------------------------------------------------------------------
say("Z149: Trägt die E8-Lesart des z120-Gatters auf Operatorebene?")
say("=" * 100)
say("\nTeil A: Modell. Stufe-1-Darstellung von E8 auf dem Gitter, Niveau L_0 = 1 (248 Zustände)")

M8 = Modell(0)
R = M8.R
typ = {}
for r in R:
    k = (sum(x * x for x in r[:5]), sum(x * x for x in r[5:]))
    typ[k] = typ.get(k, 0) + 1
kleb = [r for r in R if r[0] % 2]                      # Spinorklasse = Klebewurzeln (16,4) + (16-quer,4-quer)
gerade = [r for r in R if r[0] % 2 == 0]
check("[E] E8 = D8 + Spinorklasse, D5 auf den Koordinaten 1 bis 5, A3 = D3 auf 6 bis 8: 240 Wurzeln = 40 (D5) + 12 (A3) "
      "+ 60 ((10,6)) + 128 Klebewurzeln; Cartan-Matrix der einfachen Wurzeln hat Determinante 1",
      len(R) == 240 and typ == {(8, 0): 40, (0, 8): 12, (4, 4): 60, (5, 3): 128} and len(kleb) == 128
      and round(float(np.linalg.det(M8.cart))) == 1)
ok = all(M8.eps(a, a) == -1 for a in R) and all(M8.eps(a, b) * M8.eps(b, a) == (-1) ** ip(a, b) for a in R[::7] for b in R)
check("[E] Kozykel eps(a,b) = (-1)^(a^T eta b), bilinear in den Koordinaten der einfachen Wurzeln: eps(a,a) = -1 und "
      "eps(a,b) eps(b,a) = (-1)^(a.b)", ok)

rng = np.random.default_rng(149)
probe = [R[i] for i in sorted(rng.choice(240, 40, replace=False).tolist())]
spE = {a: sps.csr_matrix(M8.E(a)) for a in R}
spH = [sps.csr_matrix(h) for h in M8.H]
err_h, err_c, err_adj, n_kom, n_typ = 0.0, 0.0, 0.0, 0, {"Summe ist Wurzel": 0, "b = -a": 0, "sonst null": 0}
for a in R:
    err_adj = max(err_adj, abs(spE[a].T + spE[tuple(-x for x in a)]).max())
    for i in range(8):
        d = spH[i] @ spE[a] - spE[a] @ spH[i] - a[i] / 2 * spE[a]
        err_h = max(err_h, abs(d).max() if d.nnz else 0.0)
for a in probe:
    for b in R:
        kom = spE[a] @ spE[b] - spE[b] @ spE[a]
        s = ip(a, b)
        if s == -1:
            soll = M8.eps(a, b) * spE[tuple(x + y for x, y in zip(a, b))]
            n_typ["Summe ist Wurzel"] += 1
        elif s == -2:
            soll = M8.eps(a, b) * sps.csr_matrix(M8.cartan(a))
            n_typ["b = -a"] += 1
        else:
            soll = sps.csr_matrix((248, 248))
            n_typ["sonst null"] += 1
        d = kom - soll
        err_c = max(err_c, abs(d).max() if d.nnz else 0.0)
        n_kom += 1
cas = sum(h @ h for h in M8.H) + sum(M8.E(a) @ M8.E(a).T for a in R)
ZAHLEN["Kommutatoren der Stichprobe"] = {"gesamt": n_kom, **n_typ}
check("[E] Nullmoden erfüllen E8: [H_i, E_a] = a_i E_a für alle 240 Wurzeln; [E_a, E_b] = eps(a,b) E_(a+b), "
      "eps(a,-a) a.H oder 0 für 40 Wurzeln a gegen alle 240 b (9600 Kommutatoren, davon 2240 mit Wurzelsumme); "
      "E_a^+ = -E_(-a); Casimir = 60 = 2 h^vee auf allen 248 Zuständen",
      err_h < TOL and err_c < TOL and err_adj < TOL and n_kom == 9600 and n_typ["Summe ist Wurzel"] == 2240
      and nrm(cas - 60 * np.eye(248)) < TOL)

Mk = Modell(0, kozykel=False)                          # Kontrolle: ohne Kozykel
err_k = 0.0
for a in probe[:6]:
    for b in R[::4]:
        if ip(a, b) == -1:
            kom = Mk.E(a) @ Mk.E(b) - Mk.E(b) @ Mk.E(a)
            ziel = Mk.E(tuple(x + y for x, y in zip(a, b)))
            err_k = max(err_k, min(nrm(kom - ziel), nrm(kom + ziel)))
check("[E] Kontrolle: Ohne Kozykel (eps = 1) erfüllen dieselben Formeln E8 nicht ([E_a, E_b] ist dann nicht +-E_(a+b))",
      err_k > 0.5)

# Referenz-Klebewurzel: in der 16 mit 4 lambda_s . gamma = 1 und gamma_1 = +1/2
g0 = (1, 1, 1, -1, -1, 1, 1, 1)
assert g0 in M8.idx and sum(g0[:5]) == 1
Jx, Jy, Jz = M8.su2(g0)
C2 = Jx @ Jx + Jy @ Jy + Jz @ Jz
P0 = (C2 - 0.75 * np.eye(248)) @ (C2 - 2 * np.eye(248)) / 1.5
Ph = C2 @ (C2 - 2 * np.eye(248)) / (0.75 * -1.25)
P1 = C2 @ (C2 - 0.75 * np.eye(248)) / (2 * 1.25)
dims = [int(round(np.trace(P).real)) for P in (P0, Ph, P1)]
check("[E] su(2) der Klebewurzel gamma: J_+ = E_gamma, J_- = E_gamma^+ = -E_(-gamma), J_z = gamma.H/2 mit [J_x, J_y] = "
      "i J_z; 248 = 133 Singuletts + 56 Dubletts + 1 Triplett (E7 x SU(2): (133,1) + (56,2) + (1,3))",
      nrm(Jx @ Jy - Jy @ Jx - 1j * Jz) < TOL and nrm(Jy @ Jz - Jz @ Jy - 1j * Jx) < TOL and dims == [133, 112, 3])
befund("Modell", "exakt",
       "Teil A: Das Niveau L_0 = 1 der Stufe-1-Darstellung ist die adjungierte Darstellung, gebaut aus den "
       "Vertexoperatorformeln mit Kozykel. Spins 0, 1/2 und 1 der Klebe-su(2) kommen vor; höhere Spins liefert Teil C4.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil B: Die Drehungen der Faktoren als Operatoren")
L5 = np.array([sum(x * x for x in r[:5]) / 8 for r in R] + [1.0] * 5 + [0.0] * 3)      # L_0^D5 (Impuls + Oszillator)
L3 = np.array([sum(x * x for x in r[5:]) / 8 for r in R] + [0.0] * 5 + [1.0] * 3)      # L_0^A3
kl5 = np.array([sum(r[:5]) % 4 for r in R] + [0] * 8)                                  # Klasse a = 4 lambda_s . p mod 4
kl3 = np.array([sum(r[5:]) % 4 for r in R] + [0] * 8)
A1 = np.array([r[0] for r in R] + [0] * 8)                                             # A = 2 p_1
sig = np.diag((-1.0) ** A1).astype(complex)                                            # sigma = e^(2 pi i H_1)
zc = np.diag(1j ** kl5)                                                                # z = e^(2 pi i lambda_s . H)
W = np.diag(np.exp(4j * np.pi * L5))
T5 = np.diag(np.exp(2j * np.pi * L5))
T3 = np.diag(np.exp(2j * np.pi * L3))
Id = np.eye(248, dtype=complex)
check("[E] Klassen: jede der 248 Zustände hat L_0^D5 = 5 a^2/8 mod 1 und L_0^A3 = 3 a'^2/8 mod 1; in diesem Gitter gilt "
      "a + a' = 0 mod 4 (Lagrangesche <(1,3)>), und L_0^D5 + L_0^A3 = 1: T5 T3 = e^(2 pi i L_0) = 1",
      all(abs((L5[i] - 5 * kl5[i] ** 2 / 8 + 0.5) % 1 - 0.5) < TOL and abs((L3[i] - 3 * kl3[i] ** 2 / 8 + 0.5) % 1 - 0.5)
          < TOL and (kl5[i] + kl3[i]) % 4 == 0 for i in range(248)) and nrm(T5 @ T3 - Id) < TOL
      and np.allclose(L5 + L3, 1))
check("[E] W = e^(4 pi i L_0^D5) = i^(a^2) = (1+i)/2 + (1-i)/2 sigma mit sigma = e^(2 pi i H_1) = (-1)^a: W ist die "
      "Phasenwurzel der Involution sigma (W^2 = sigma, W^4 = 1). W = T5 T3^(-1) = T5^2 (relative Drehung); T5^4 = sigma, "
      "T5 hat Ordnung 8",
      nrm(W - np.diag(1j ** (kl5 ** 2 % 4))) < TOL and nrm(W - ((1 + 1j) / 2 * Id + (1 - 1j) / 2 * sig)) < TOL
      and nrm(W @ W - sig) < TOL and nrm(W - T5 @ dag(T3)) < TOL and nrm(W - T5 @ T5) < TOL
      and nrm(np.linalg.matrix_power(T5, 4) - sig) < TOL and nrm(np.linalg.matrix_power(T5, 8) - Id) < TOL
      and nrm(np.linalg.matrix_power(T5, 4) - Id) > 1)
fix_sig = sum(1 for r in R if r[0] % 2 == 0) + 8
fix_z = sum(1 for r in R if sum(r[:5]) % 4 == 0) + 8
ZAHLEN["Fixalgebren"] = {"sigma (Ordnung 2)": fix_sig, "z (Ordnung 4)": fix_z}
ok_inn = all(nrm(sig @ M8.E(a) @ sig - (-1) ** a[0] * M8.E(a)) < TOL for a in R[::5])
check("[E] sigma ist ein innerer Automorphismus von E8 der Ordnung 2 (Toruselement e^(2 pi i e_1.H)) mit Fixalgebra der "
      "Dimension 120 = so(16); z = e^(2 pi i lambda_s.H) = i^a hat Ordnung 4 und Fixalgebra der Dimension 60 = "
      "so(10) + su(4) (45 + 15)", ok_inn and fix_sig == 120 and fix_z == 60)
err_ev, err_cl, err_cl5, n16 = 0.0, 0.0, 0.0, 0
for a in gerade:
    err_ev = max(err_ev, nrm(W @ M8.E(a) - M8.E(a) @ W))
for g in kleb:
    Eg = M8.E(g)
    err_cl = max(err_cl, nrm(W @ Eg @ dag(W) - 1j * Eg @ sig), nrm(dag(W) @ Eg @ W + 1j * Eg @ sig))
    up = sum(g[:5]) % 4 == 1                          # gamma in der 16 (a -> a + 1) oder in der 16-quer (a -> a - 1)
    n16 += up
    for k in range(8):
        Tk = np.linalg.matrix_power(T5, k)
        zk = np.linalg.matrix_power(zc, k if up else (-k) % 4)
        err_cl5 = max(err_cl5, nrm(Tk @ Eg @ dag(Tk) - np.exp(5j * np.pi * k / 4) * Eg @ zk))
check("[E] Clifford-Form auf Operatorebene: W vertauscht mit allen 112 Strömen von D8 und bildet jeden der 128 "
      "Klebeströme auf Transfer mal Uhrwert ab, W E_gamma W^(-1) = i E_gamma sigma; W^(-1) E_gamma W = -i E_gamma sigma "
      "ist wörtlich z145 B3 (U X U^+ = -i X Gamma^(-1)) mit der Paaruhr sigma", err_ev < TOL and err_cl < TOL)
check("[E] Die acht vollen Drehungen T5^k erfüllen alle die Clifford-Form: T5^k E_gamma T5^(-k) = e^(5 pi i k/4) E_gamma "
      "z^(+-k) (64 Klebewurzeln in der 16, 64 in der 16-quer). Die 2 pi-Drehung von D5 allein trägt die volle Uhr z "
      "(Ordnung 4), die 4 pi-Drehung ihr Quadrat sigma", err_cl5 < TOL and n16 == 64)
# Kein Automorphismus: Bild von E_gamma liegt nicht in e8
gen = [M8.E(a) for a in R] + M8.H
gnorm = [float(np.sum(x * x)) for x in gen]
Xb = W @ M8.E(g0) @ dag(W)
koef = [np.sum(x * Xb) / n for x, n in zip(gen, gnorm)]
rest = np.sqrt(max(0.0, nrm(Xb) ** 2 - sum(abs(c) ** 2 * n for c, n in zip(koef, gnorm)))) / nrm(Xb)
Kf = np.zeros((248, 248))
for r in R:
    Kf[M8.idx[tuple(-x for x in r)], M8.idx[r]] = 1.0
for i in range(8):
    Kf[240 + i, 240 + i] = -1.0
reell = lambda X: nrm(Kf @ X.conj() @ Kf - X)          # Realitätsbedingung von so(248) bezüglich der invarianten Form
antih = [M8.E(g0) - M8.E(g0).T, 1j * (M8.E(g0) + M8.E(g0).T), 1j * M8.H[0], M8.E(R[0]) - M8.E(R[0]).T]
ZAHLEN["Anteil von W E_gamma W^-1 außerhalb von e8 (relative Norm)"] = rd(rest)
check("[E] W ist kein Automorphismus der E8-Ströme: W E_gamma W^(-1) = i E_gamma sigma liegt nicht in e8 (in der "
      "Spurform orthogonal zu allen 248 Erzeugern, relativer Rest 1). Ad(W) bildet e8 auf die verdrehte Kopie e8_gerade + i e8_ungerade sigma "
      "ab, die die Realitätsbedingung von so(248) verletzt (e8 selbst erfüllt sie)",
      abs(rest - 1) < 1e-8 and all(reell(x) < TOL for x in antih) and reell(W @ antih[0] @ dag(W)) > 1)
befund("4 pi-Drehung als Symmetrie", "nein",
       "Teil B: e^(4 pi i L_0^D5) ist die Phasenwurzel (1+i)/2 + (1-i)/2 sigma des inneren Automorphismus sigma "
       "(Ordnung 2, Fixalgebra so(16)); selbst ist sie kein Automorphismus. Auf Klebeströmen wirkt sie als Transfer mal "
       "Uhrwert, wie P_uhr in der Clifford-Form es verlangt.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil C: Wörtliche Prüfung gegen das z120-Gatter")
# C1: z120 im Fockmodell, Jordan-Schwinger-Form und das Gegenstück der E8-Formel
resF = {}
for K in (2, 3):
    c = fock(2 * K)
    cA, cB = c[:K], c[K:]
    NA = sum(dag(x) @ x for x in cA)
    NB = sum(dag(x) @ x for x in cB)
    nA, nB = np.rint(np.diag(NA).real).astype(int), np.rint(np.diag(NB).real).astype(int)
    Jp = sum(dag(a) @ b for a, b in zip(cA, cB))
    Fx, Fy = (Jp + dag(Jp)) / 2, (Jp - dag(Jp)) / 2j
    Bq = expm(-1j * np.pi / 2 * Fx)
    Nt = NA + NB
    fn = lambda f: dag(Bq) @ np.diag(f(nA, nB).astype(complex)) @ Bq
    Uv = fn(lambda a, b: np.exp(-1j * np.pi * (a * a + b * b) / 4))                    # U(pi/2)
    Ur = fn(lambda a, b: np.exp(1j * np.pi * (a * a + b * b) / 4))                     # U(-pi/2)
    r = {}
    r["JS"] = nrm(Uv - expm(-1j * np.pi / 8 * Nt @ Nt) @ expm(-1j * np.pi / 2 * Fy @ Fy)) \
        + nrm(dag(Bq) @ (NA - NB) / 2 @ Bq - Fy)
    r["eine Naht"] = nrm(fn(lambda a, b: 1j ** (a * a % 4)) - Ur @ expm(1j * np.pi / 2 * Nt @ Fy))
    th5 = lambda q: np.exp(2j * np.pi * 5 * q * q / 8)
    r["zwei Nähte"] = nrm(fn(lambda a, b: th5(a) * th5(b)) - expm(1j * np.pi * Nt) @ Ur)
    resF[K] = r
check("[E] z120 im Fockmodell (K = 2, 3): U(phi) = B^+ e^(-i phi H_Q) B = e^(-i phi N^2/4) e^(-i phi J_y^2) mit B^+ J_z B "
      "= J_y. Dort gilt auch: die gehobene 4 pi-Drehung einer Naht allein ist B^+ i^(Q_A^2) B = U(-pi/2) e^(i pi N J_y/2), "
      "und die Formfunktion beider Nähte B^+ theta_D5(Q_A) theta_D5(Q_B) B = (-1)^N U(-pi/2) (z145 C)",
      max(v for r in resF.values() for v in r.values()) < 1e-8)


def gatter(mod, g, Wd, A_diag, pw_cache=None):
    """B_g^+ W B_g und die Vergleichsoperatoren für die Klebewurzel g."""
    jx, jy, jz = mod.su2(g)
    s = g[0]                                           # s = 2 gamma_1 = +-1
    pwx, pwy = potenzen(jx, 5), potenzen(jy, 5)
    B = fpoly(jx, lambda x: np.exp(-1j * np.pi * x / 2), pw=pwx)
    Md = A_diag - s * np.diag(jz).real                 # M = 2 H_1 - s J_z, vertauscht mit der su(2)
    Md = (Md + 0.5) % 2 - 0.5                          # W kennt nur N = 2M mod 4: Vertreter M in {0, 1, 1/2, -1/2}
    lhs = dag(B) @ Wd @ B
    ref = np.zeros((248, 248), complex)
    rhs = np.zeros((248, 248), complex)
    for m in sorted(set(np.round(Md, 6).tolist())):
        Pm = np.diag((np.abs(Md - m) < 1e-6).astype(complex))
        ref += np.exp(1j * np.pi * m * m / 2) * Pm @ fpoly(jy, lambda x: np.exp(1j * np.pi * x * x / 2), pw=pwy)
        rhs += Pm @ fpoly(jy, lambda x, m=m: np.exp(1j * np.pi * (m + s * x) ** 2 / 2), pw=pwy)
    return dict(B=B, lhs=lhs, ref=ref, rhs=rhs, M=Md, jx=jx, jy=jy, jz=jz, s=s, pwx=pwx, pwy=pwy)


Wd = W
G0 = gatter(M8, g0, Wd, A1.astype(float))
Gam = fpoly(G0["jy"], lambda x: np.exp(1j * np.pi * x), pw=G0["pwy"])                  # Nahttausch e^(i pi J_y)
cM = np.diag(np.exp(1j * np.pi * G0["M"]))                                             # c = e^(i pi M) = i^N
af, bf = (1 - 1j) / 2, (1 + 1j) / 2
check("[E] Vierteldrehung: B_gamma^+ J_z B_gamma = J_y; die konjugierte Drehung ist B_gamma^+ W B_gamma = (1+i)/2 + "
      "(1-i)/2 c Gamma_y mit Gamma_y = e^(i pi J_y) und c = sigma e^(-i pi J_z), das mit der su(2) vertauscht "
      "(Referenzwurzel, s = 2 gamma_1 = +1)",
      nrm(dag(G0["B"]) @ G0["jz"] @ G0["B"] - G0["jy"]) < TOL
      and nrm(G0["lhs"] - (bf * Id + af * cM @ Gam)) < TOL and nrm(cM @ G0["jx"] - G0["jx"] @ cM) < TOL
      and nrm(cM - sig @ np.diag(np.exp(-1j * np.pi * np.diag(G0["jz"]).real))) < TOL)

# alle 128 Klebewurzeln
err_id, err_s1, err_s0, dev_h, cz = 0.0, 0.0, 0.0, set(), {}
for g in kleb:
    G = gatter(M8, g, Wd, A1.astype(float)) if g != g0 else G0
    err_id = max(err_id, nrm(G["lhs"] - G["rhs"]))
    cc = G["jx"] @ G["jx"] + G["jy"] @ G["jy"] + G["jz"] @ G["jz"]
    p1 = cc @ (cc - 0.75 * Id) / 2.5
    p0 = (cc - 0.75 * Id) @ (cc - 2 * Id) / 1.5
    ph = cc @ (cc - 2 * Id) / (0.75 * -1.25)
    D = G["lhs"] - G["ref"]
    err_s1 = max(err_s1, nrm(D @ p1))
    err_s0 = max(err_s0, nrm(D @ p0))
    dev_h.add(rd(nrm(D @ ph), 6))
    cdiag = np.round(np.exp(1j * np.pi * G["M"]), 6)
    key = tuple(sorted((str(complex(v)), int(round(np.sum((np.abs(cdiag - v) < 1e-6) * np.diag(P).real))))
                       for v in (1, -1, 1j, -1j) for P in [ph]))
    cz[key] = cz.get(key, 0) + 1
N4 = np.rint(2 * G0["M"]).astype(int) % 4
zaehl = {}
for name, P in (("Spin 0", P0), ("Spin 1/2", Ph), ("Spin 1", P1)):
    d = np.diag(P).real
    zaehl[name] = {f"N = {n} mod 4": int(round(float(np.sum(d[N4 == n])))) for n in range(4)}
ZAHLEN["Zustände je Spin und Gesamtladung N = 2M mod 4 (Referenzwurzel)"] = zaehl
ZAHLEN["Abweichung vom z120-Gatter auf den 56 Dubletts (Frobenius; Spektralnorm 2 sin(pi/8))"] = sorted(dev_h)
check("[E] Operatoridentität für alle 128 Klebewurzeln: B_gamma^+ e^(4 pi i L_0^D5) B_gamma = e^(i pi (M + s J_y)^2/2) mit "
      "M = 2 H_1 - s J_z (nur modulo 2 bestimmt; Vertreter 0, 1, +-1/2). Das z120-Gatter in derselben Schreibweise ist "
      "U(-pi/2) = e^(i pi (M^2 + J_y^2)/2) mit N = 2 M; der Quotient ist e^(i pi s M J_y) = Gamma_y^(s N/2), die "
      "(s N/2)-te Potenz des einkörperigen Nahttauschs",
      err_id < 1e-8)
check("[E] Auf dem Niveau L_0 = 1 stimmt die konjugierte Drehung mit dem z120-Gatter überein auf dem Triplett (E_gamma, "
      "gamma.H, E_(-gamma)) und den 133 Singuletts, für alle 128 Klebewurzeln; auf den 56 Dubletts (N ungerade) weicht "
      "sie um die einkörperige Wurzel des Tauschs e^(+-i pi J_y/2) ab (Spektralnorm 2 sin(pi/8) = 0,765367)",
      err_s1 < 1e-8 and err_s0 < 1e-8 and sorted(dev_h) == [rd(np.sqrt(112) * 2 * np.sin(np.pi / 8), 6)]
      and zaehl["Spin 1"] == {"N = 0 mod 4": 3, "N = 1 mod 4": 0, "N = 2 mod 4": 0, "N = 3 mod 4": 0}
      and zaehl["Spin 1/2"]["N = 0 mod 4"] == 0 and zaehl["Spin 1/2"]["N = 2 mod 4"] == 0)
trip = dag(G0["B"]) @ dag(Wd) @ G0["B"]
check("[E] Form a + b mal Tausch: Auf dem Triplett ist B_gamma^+ e^(-4 pi i L_0^D5) B_gamma = a + b Gamma_y mit a = "
      "(1-i)/2, b = (1+i)/2 und Sektorphase 1 (das Vorwärtsgatter U(+pi/2) bei N = 0 mod 4); mit e^(+4 pi i L_0^D5) "
      "steht dort (1+i)/2 + (1-i)/2 Gamma_y, das inverse Gatter",
      nrm((trip - (af * Id + bf * Gam)) @ P1) < TOL and nrm((G0["lhs"] - (bf * Id + af * Gam)) @ P1) < TOL)

# C4: höhere Spins. Kette p + m gamma mit Oszillatoren längs gamma (su(2)_1-Teilmodell, exakt)
LMAX = 9


def partitionen(n, L=LMAX):
    out = []

    def rec(rest, teil, acc):
        if teil == 0:
            if rest == 0:
                out.append(tuple(acc[::-1]))
            return
        for m in range(rest // teil, -1, -1):
            rec(rest - m * teil, teil - 1, acc + [m])
    rec(n, min(n, L), [])
    return sorted(tuple(list(t) + [0] * (L - len(t))) for t in out)


def grad(e):
    return sum((n + 1) * m for n, m in enumerate(e))


def pmul(p, q):
    out = {}
    for e1, c1 in p.items():
        for e2, c2 in q.items():
            e = tuple(x + y for x, y in zip(e1, e2))
            out[e] = out.get(e, 0) + c1 * c2
    return out


def verschiebe(e, cst):
    """Monom prod y_n^(m_n) mit y_n -> y_n + cst."""
    out = {tuple([0] * LMAX): Fr(1)}
    for n, m in enumerate(e):
        if m:
            fak = {}
            for j in range(m + 1):
                ee = [0] * LMAX
                ee[n] = j
                fak[tuple(ee)] = Fr(comb(m, j)) * Fr(cst) ** (m - j)
            out = pmul(out, fak)
    return out


_schur = {}


def schur(vz, d):
    """Homogener Teil vom Grad d von exp(vz sum_n y_n/n)."""
    if (vz, d) not in _schur:
        p = {}
        for e in partitionen(d):
            c = Fr(1)
            for n, m in enumerate(e):
                c *= Fr(vz, n + 1) ** m / factorial(m)
            p[e] = c
        _schur[(vz, d)] = p
    return _schur[(vz, d)]


def norm2(e):
    out = 1
    for n, m in enumerate(e):
        out *= (2 * (n + 1)) ** m * factorial(m)
    return out


def kettensektor(par, lev):
    """Zustände |k; lambda> mit k = gamma.p = par mod 2 und L_0 = k^2/4 + |lambda| = par/4 + lev; J_+, J_-, k."""
    st = []
    kmax = int(np.floor(2 * np.sqrt(par / 4 + lev))) + 2
    for k in range(-kmax, kmax + 1):
        if (k - par) % 2 == 0:
            d = Fr(par, 4) + lev - Fr(k * k, 4)
            if d >= 0 and d.denominator == 1:
                st += [(k, e) for e in partitionen(int(d))]
    pos = {s: i for i, s in enumerate(st)}
    Ep, Em = np.zeros((len(st), len(st))), np.zeros((len(st), len(st)))
    for (k, e), j in pos.items():
        d = grad(e)
        eta = (-1) ** ((k - par) // 2)
        for vz, Mat, ziel in ((1, Ep, d - k - 1), (-1, Em, d + k - 1)):
            if ziel < 0:
                continue
            bild = {}
            for e1, c1 in verschiebe(e, -2 * vz).items():
                if ziel - grad(e1) >= 0:
                    for e2, c2 in pmul({e1: c1}, schur(vz, ziel - grad(e1))).items():
                        bild[e2] = bild.get(e2, 0) + c2
            for e2, c2 in bild.items():
                if c2:
                    Mat[pos[(k + 2 * vz, e2)], j] = eta * float(c2) * np.sqrt(norm2(e2) / norm2(e))
    return st, Ep, -Em, np.array([k for k, _ in st], float)


sektoren = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 0), (1, 1), (1, 2), (1, 3), (1, 6)]
tab, ok_su2, ok_idk, ok_rev, spins_ges, dimsum = {}, True, True, True, {}, 0
for par, lev in sektoren:
    st, Jp_, Jm_, kk = kettensektor(par, lev)
    dimsum += len(st)
    ok_su2 &= nrm(Jm_ - Jp_.T) < 1e-9 and nrm(Jp_ @ Jm_ - Jm_ @ Jp_ - np.diag(kk)) < 1e-9
    jx_, jy_, jz_ = (Jp_ + Jm_) / 2, (Jp_ - Jm_) / 2j, np.diag(kk) / 2
    ew, ev = np.linalg.eigh(jx_)
    B_ = (ev * np.exp(-1j * np.pi * np.round(2 * ew) / 4)) @ ev.T
    ewy, evy = np.linalg.eigh(jy_)
    ewy = np.round(2 * ewy) / 2
    ok_su2 &= nrm(dag(B_) @ jz_ @ B_ - jy_) < 1e-8
    cas_ = jx_ @ jx_ + jy_ @ jy_ + jz_ @ jz_
    ewc, evc = np.linalg.eigh(cas_)
    jw = np.round(np.sqrt(4 * ewc + 1) - 1) / 2        # Spin j aus j(j+1)
    for j in sorted(set(jw.tolist())):
        spins_ges[(par, lev, j)] = int(np.sum(jw == j) / (2 * j + 1) + 0.5)
    for Mv in ((0.0, 1.0) if par == 0 else (0.5, -0.5)):
        Aw = np.rint(Mv + kk / 2).astype(int)
        W_ = np.diag(1j ** (Aw * Aw % 4))
        lhs_ = dag(B_) @ W_ @ B_
        rhs_ = (evy * np.exp(1j * np.pi * (Mv + ewy) ** 2 / 2)) @ dag(evy)
        ref_ = (evy * np.exp(1j * np.pi * (Mv * Mv + ewy ** 2) / 2)) @ dag(evy)
        ok_idk &= nrm(lhs_ - rhs_) < 1e-8
        if Mv == 1.0:
            ok_rev &= nrm(lhs_ + dag(ref_)) < 1e-8      # N = 2 mod 4: -U(+pi/2), umgekehrte Zeitrichtung
        cname = {0.0: "c = +1 (N = 0 mod 4)", 1.0: "c = -1 (N = 2 mod 4)", 0.5: "c = +i (N = 1 mod 4)",
                 -0.5: "c = -i (N = 3 mod 4)"}[Mv]
        for j in sorted(set(jw.tolist())):
            Pj = evc[:, jw == j]
            dv = rd(np.linalg.norm(dag(Pj) @ (lhs_ - ref_) @ Pj, 2), 6)
            tab.setdefault((cname, j), set()).add(dv)
ok_tab = all(len(v) == 1 for v in tab.values())
tabelle = {c: {f"Spin {j + 0.0:g}": sorted(v)[0] for (cc, j), v in sorted(tab.items()) if cc == c}
           for c in sorted({k[0] for k in tab})}
ZAHLEN["Abweichung ||B^+ W B - U(-pi/2)|| je Sektor und Spin (Kettenmodell)"] = tabelle
hoch = {f"{'gerade' if p == 0 else 'ungerade'} k, Niveau {l}": {f"Spin {j + 0.0:g}": n for (pp, ll, j), n in
        sorted(spins_ges.items()) if (pp, ll) == (p, l)} for p, l in sektoren}
ZAHLEN["Spininhalt der Kettensektoren"] = hoch
s2, s3 = rd(2 * np.sin(np.pi / 8)), rd(2 * np.sin(3 * np.pi / 8))
soll = {"c = +1 (N = 0 mod 4)": {0.0: 0.0, 1.0: 0.0, 2.0: 0.0}, "c = -1 (N = 2 mod 4)": {0.0: 0.0, 1.0: 2.0, 2.0: 2.0},
        "c = +i (N = 1 mod 4)": {0.5: s2, 1.5: s3, 2.5: s3}, "c = -i (N = 3 mod 4)": {0.5: s2, 1.5: s3, 2.5: s3}}
ok_soll = all(sorted(tab[(c, j)])[0] == v for c, d in soll.items() for j, v in d.items())
check("[E] Kettenmodell p + m gamma mit Oszillatoren längs gamma (neun Sektoren, 79 Zustände): J_- = J_+^+ und "
      "[J_+, J_-] = gamma.p; es enthält die Spins 0, 1, 2 (gerades gamma.p) und 1/2, 3/2, 5/2 (ungerades). Die "
      "Operatoridentität aus C2 gilt in allen vier Sektoren c = 1, -1, i, -i",
      ok_su2 and ok_idk and dimsum == 79 and spins_ges[(0, 4, 2.0)] == 1 and spins_ges[(1, 2, 1.5)] == 1
      and spins_ges[(1, 6, 2.5)] == 1)
check("[E] Wo es abweicht: Für c = +1 (N = 0 mod 4) ist die konjugierte Drehung auf allen Spins genau das z120-Gatter. "
      "Für c = -1 (N = 2 mod 4) ist sie -U(+pi/2), also das Gatter in umgekehrter Zeitrichtung (Abstand 2 ab Spin 1). "
      "Für c = +-i (N ungerade) trägt sie zusätzlich e^(+-i pi J_y/2) (Abstand 0,765367 bei Spin 1/2, 1,847759 darüber)",
      ok_tab and ok_soll and ok_rev)

# alle vier Sektoren kommen im E8-Gitter vor
vek = []
for v in itertools.product(range(-2, 3), repeat=8):
    if sum(x * x for x in v) == 4 and sum(v) % 2 == 0:
        vek.append(tuple(2 * x for x in v))
for v in itertools.product((-3, -1, 1, 3), repeat=8):
    if sum(x * x for x in v) == 16 and sum(v) % 4 == 0:
        vek.append(v)
vork = {}
for nm, liste in (("Norm 2", R), ("Norm 4", vek)):
    z = {}
    for v in liste:
        k = sum(x * y for x, y in zip(g0, v)) // 4
        key = f"N = {(2 * v[0] - g0[0] * k) % 4} mod 4, |gamma.p| = {abs(k)}"
        z[key] = z.get(key, 0) + 1
    vork[nm] = dict(sorted(z.items()))
ZAHLEN["Gittervektoren nach N = 4 omega_senkrecht . p mod 4 und gamma.p"] = vork
check("[E] Der Niveauabgleich N = 0 mod 4 der Klebe-su(2) ist in E8 nicht eingebaut: Unter den 240 Wurzeln und den 2160 "
      "Vektoren der Norm 4 kommen alle vier Werte von N mod 4 vor, darunter Tripletts mit N = 2 mod 4 (|gamma.p| = 2 "
      "bei Norm 4)", len(vek) == 2160 and all(any(k.startswith(f"N = {n} mod 4") for k in vork["Norm 4"]) for n in range(4))
      and vork["Norm 4"].get("N = 2 mod 4, |gamma.p| = 2", 0) > 0)

# C5: die 2 pi-Drehung allein und die Fassung mit zwei gleichartigen Nähten
h5 = dag(G0["B"]) @ T5 @ G0["B"]
w8 = np.exp(5j * np.pi / 4)
a5, b5 = (1 + w8) / 2, (1 - w8) / 2
check("[E] Die 2 pi-Drehung von D5 allein: B_gamma^+ T5 B_gamma ist auf dem Triplett e^(i 5 pi J_y^2/4) = a' + b' Gamma_y "
      "mit |a'|^2 = sin^2(pi/8) = 0,146447 und |b'|^2 = 0,853553: ein ungleich gewichteter Teiltausch mit achten "
      "Wurzeln (U(-5 pi/4) = U(-pi/4) mal Tausch), dessen Quadrat das Gatter ist. Das ist z145s 'halbes Gatter'",
      nrm((h5 - (a5 * Id + b5 * Gam)) @ P1) < TOL and abs(abs(a5) ** 2 - np.sin(np.pi / 8) ** 2) < TOL
      and nrm((h5 @ h5 - G0["lhs"]) @ P1) < TOL)
ZAHLEN["Gewichte des Teiltauschs der 2 pi-Drehung"] = [rd(abs(a5) ** 2), rd(abs(b5) ** 2)]
lam4 = (2, 2, 2, 2, 2, 0, 0, 0)                        # 4 lambda_s, verdoppelt
err_zw, n_zw, schnitt = 0.0, 0, set()
for g in kleb:
    kq = sum(g[:5])                                    # k = 4 lambda_s . gamma
    if kq != 1:
        continue
    n_zw += 1
    G = gatter(M8, g, Wd, A1.astype(float)) if g != g0 else G0
    QA = np.array([sum(r[:5]) for r in R] + [0] * 8, float)                            # 4 lambda_s . p
    QB = QA - 2 * np.diag(G["jz"]).real                                                # 4 s_gamma(lambda_s) . p
    prod = np.diag(np.exp(2j * np.pi * 5 * (QA ** 2 + QB ** 2) / 8))                   # T5 mal gespiegeltes T5
    Mq = (QA + QB) / 2
    refq = np.zeros((248, 248), complex)
    for m in sorted(set(np.round(Mq, 6).tolist())):
        Pm = np.diag((np.abs(Mq - m) < 1e-6).astype(complex))
        refq += np.exp(1j * np.pi * (2 * m + m * m / 2)) * Pm @ fpoly(G["jy"], lambda x: np.exp(1j * np.pi * x * x / 2),
                                                                      pw=G["pwy"])
    n2 = dag(G["B"]) @ dag(G["B"])
    err_zw = max(err_zw, nrm(dag(G["B"]) @ prod @ G["B"] - refq), nrm(n2 @ T5 @ dag(n2) - np.diag(
        np.exp(2j * np.pi * 5 * QB ** 2 / 8))))
    schnitt.add(sum(1 for r in R if r[0] % 2 == 0 and not any(r[5:]) and ip(r, g) == 0))
check("[E] Fassung mit zwei gleichartigen Ladungen: Q_A = 4 lambda_s.H (Klasse von D5) und Q_B = Q_A - gamma.H (Klasse "
      "des gespiegelten D5' = s_gamma(D5)). Dann ist B_gamma^+ T5 T5' B_gamma = (-1)^N U(-pi/2) auf allen 248 Zuständen, "
      "ohne Abweichung (40 Klebewurzeln mit 4 lambda_s.gamma = 1). D5 und D5' sind nicht unabhängig: Sie teilen die 20 "
      "Wurzeln eines su(5)", err_zw < 1e-8 and n_zw == 40 and schnitt == {20})
kontr = {}
for j in (1, 2, 3):
    th = j * np.pi / 8
    Bt = fpoly(G0["jx"], lambda x, th=th: np.exp(-1j * th * x), pw=G0["pwx"])
    Ut = (dag(Bt) @ Wd @ Bt) @ P1
    cols = np.array([P1.ravel(), (Gam @ P1).ravel()]).T
    kontr[j] = rd(np.linalg.norm(Ut.ravel() - cols @ np.linalg.lstsq(cols, Ut.ravel(), rcond=None)[0]))
ZAHLEN["Kontrolle: Abstand von B(theta)^+ W B(theta) zu span{1, Tausch} auf dem Triplett, theta = j pi/8"] = kontr
check("[E] Kontrolle: Mit einem anderen Transferwinkel (theta = pi/8, pi/4, 3 pi/8 statt pi/2) hat die konjugierte "
      "Drehung auf dem Triplett nicht die Form a + b mal Tausch", all(v > 0.1 for v in kontr.values()))
befund("wörtliche Gleichheit", "bedingt (N = 0 mod 4)",
       "Teil C: B_gamma^+ e^(4 pi i L_0^D5) B_gamma = U(-pi/2) Gamma_y^(s N/2). Gleichheit mit dem z120-Gatter gilt auf "
       "Multipletts mit N = 0 mod 4 (und auf allen Singuletts); bei N = 2 mod 4 läuft die Zeit umgekehrt, bei ungeradem N kommt die einkörperige "
       "Wurzel des Tauschs hinzu. Die 4 pi-Drehung ist die Paarphase einer einzigen zirkularen Mode, i^(Q_+^2).")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil D: Die Prämissen des Kandidaten einzeln")
# D1 P_uhr / P_null: Familie e^(i (t5 L_0^D5 + t3 L_0^A3)), Klassenfunktion?
dat = sorted({(int(kl5[i]), round(float(L5[i]), 6), round(float(L3[i]), 6)) for i in range(248)} | {(0, 0.0, 0.0)})
klf, ops = [], set()
for j5 in range(32):
    for j3 in range(32):
        werte = {}
        for a, l5, l3 in dat:
            werte.setdefault(a, set()).add(np.round(np.exp(1j * np.pi / 2 * (j5 * l5 + j3 * l3)), 6).item())
        if all(len(v) == 1 for v in werte.values()):
            klf.append((j5, j3))
            ops.add(tuple(str(sorted(werte[a])[0]) if False else complex(list(werte[a])[0]) for a in range(4)))
check("[E] P_uhr und P_null im E8-Bild: In der Familie e^(i (t5 L_0^D5 + t3 L_0^A3)) (Gitter (pi/2) Z, mit Vakuum und "
      "Oszillatoren) ist der Operator genau dann eine Funktion der Klassen, wenn t5 und t3 in 2 pi Z liegen. Das gibt "
      "die Z_8 der vollen Drehungen T5^k; der Oszillatoranteil fällt dann von selbst weg (P_null wie in z144)",
      sorted(klf) == [(a, b) for a in range(0, 32, 4) for b in range(0, 32, 4)] and len(ops) == 8)
ZAHLEN["volle Drehungen T5^k"] = {"k = 0": "1", "k = 4": "sigma (Automorphismus, einkörperig)",
                                  "k = 2, 6": "W, W^-1 (das Gatter, Ordnung 4)",
                                  "k = 1, 3, 5, 7": "achte Wurzeln (Formen D5, D7, D1, A3)"}
check("[E] Der Winkel ist damit nicht festgelegt: Sechs der acht vollen Drehungen sind keine Automorphismen "
      "(k = 1, 2, 3, 5, 6, 7), im Fockmodell bis auf die zentrale Parität zwei (phi = +-pi/2). Das Gatter ist k = +-2, die relative "
      "Drehung T5 T3^(-1); D5 allein (k = 1) ist ebenso eine Funktion der Uhrwerte. Die Wahl k = +-2 ist die Prämisse P_rel",
      all((nrm(np.linalg.matrix_power(T5, k) @ M8.E(g0) @ dag(np.linalg.matrix_power(T5, k))
               - np.linalg.matrix_power(sig, k // 4) @ M8.E(g0) @ np.linalg.matrix_power(sig, k // 4)) < TOL) == (k % 4 == 0)
          for k in range(8)))
# D2 P_univ
kl3_b = np.array([(sum(r[5:7]) - r[7]) % 4 for r in R] + [0] * 8)                      # A3-Klasse in der zweiten Verklebung
check("[E] P_univ: In der Faktorlesart (zirkulare Ladungen = Klassen von D5 und A3) ist der Abgleich die Verklebung "
      "selbst, und die gemeinsame Drehung T5 T3 ist 1. In der Verklebung <(1,1)> heißt er a = a' (Paaruhr 1, die native "
      "Form von P_univ aus z130), in <(1,3)> heißt er a + a' = 0 mod 4 (Gesamtladung, z145); der Wechsel der Verklebung "
      "tauscht die beiden Formen. In der su(2)-Lesart mit B_gamma ist die Gesamtladung N dagegen frei (C4): Die "
      "wörtliche Gleichheit braucht dort die Zusatzbedingung P_N4 (N = 0 mod 4)",
      all((kl5[i] + kl3[i]) % 4 == 0 and (kl5[i] - kl3_b[i]) % 4 == 0 for i in range(248))
      and zaehl["Spin 1/2"]["N = 1 mod 4"] > 0)
# D3 P_sym
B0 = G0["B"]
n0 = B0 @ B0
t0 = np.diag([(-1.0) ** ip(g0, r) for r in R] + [1.0] * 8)
cart_rest = max(nrm(n0 @ M8.H[i] @ dag(n0) - sum((np.eye(8)[i, j] - g0[i] * g0[j] / 4) * M8.H[j] for j in range(8)))
                for i in range(8))
nondiag = dag(B0) @ M8.cartan(g0) @ B0
stab = {}
for j in range(9):
    th = j * np.pi / 8
    Bt = fpoly(G0["jx"], lambda x, th=th: np.exp(-1j * th * x), pw=G0["pwx"])
    X = dag(Bt) @ G0["jz"] @ Bt
    Y = sig @ X @ sig
    stab[j] = rd(nrm(Y - np.vdot(X, Y) / np.vdot(X, X) * X) / nrm(X))
    stab[j] = (stab[j], rd(abs(np.sin(2 * th))))
antik = max(nrm((sig @ (dag(Bt) @ sig @ Bt) + (dag(Bt) @ sig @ Bt) @ sig - 2 * np.cos(th) * Id) @ Ph)
            for th in (0.3, 1.1, np.pi / 2) for Bt in [fpoly(G0["jx"], lambda x, th=th: np.exp(-1j * th * x),
                                                              pw=G0["pwx"])])
ZAHLEN["Defekt der sigma-Stabilität des gedrehten Rahmens über theta = j pi/8 (gerechnet, |sin 2 theta|)"] = stab
check("[E] P_sym: B_gamma liegt nicht in der Hebung der Weyl-Gruppe (B^+ gamma.H B = 2 J_y ist nicht diagonal). Erst "
      "B_gamma^2 = e^(-i pi J_x) normalisiert den Torus und hebt die Spiegelung s_gamma; B_gamma^4 = e^(i pi gamma.H) "
      "liegt im Torus, B_gamma^8 = 1", nrm(nondiag - np.diag(np.diag(nondiag))) > 1 and cart_rest < TOL
      and nrm(np.linalg.matrix_power(B0, 4) - t0) < TOL and nrm(np.linalg.matrix_power(B0, 8) - Id) < TOL)
check("[E] B_gamma ist die Cayley-Transformierte der sigma-ungeraden Wurzel gamma: In der Familie e^(-i theta J_x) ist "
      "der gedrehte Rahmen B^+ J_z B genau für theta in (pi/2) Z wieder sigma-stabil (Defekt |sin 2 theta|), und auf den "
      "Dubletts gilt {sigma, B^+ sigma B} = 2 cos theta. Das ist die Leckformel aus z125; E8 enthält jedes theta und "
      "erzwingt die Vierteldrehung nicht", all(abs(a - b) < 1e-6 for a, b in stab.values()) and antik < 1e-8)
# D4 P_ord
paar = {}
for g in kleb:
    s = ip(g0, g)
    if s in (-1, 0, 1) and s not in paar and g != g0:
        paar[s] = g
U0 = G0["lhs"]
komm = {}
for s in sorted(paar):
    Gs = gatter(M8, paar[s], Wd, A1.astype(float))
    komm[f"gamma.gamma' = {s}"] = rd(np.linalg.norm(U0 @ Gs["lhs"] - Gs["lhs"] @ U0, 2))
bahn, rand = {g0}, [g0]
while rand:
    neu = []
    for v in rand:
        for a in gerade:
            w = tuple(x - ip(v, a) * y for x, y in zip(v, a))
            if w not in bahn:
                bahn.add(w)
                neu.append(w)
    rand = sorted(neu)
ZAHLEN["Spektralnorm von [U_gamma, U_gamma'] (Niveau L_0 = 1)"] = komm
check("[E] P_ord: Die Gatter U_gamma = B_gamma^+ W B_gamma zu verschiedenen Klebewurzeln vertauschen nicht (Spektralnorm "
      "des Kommutators sqrt(3)/2 für gamma.gamma' = +-1), auch nicht für orthogonale Wurzeln (Norm 1; dort vertauschen "
      "die Involutionen nur bis auf das Toruselement e^(i pi (gamma + gamma').H)). "
      "Die Weyl-Gruppe von D8 wirkt transitiv auf den 128 Klebewurzeln: E8 zeichnet keine Reihenfolge aus",
      komm == {"gamma.gamma' = -1": rd(np.sqrt(3) / 2), "gamma.gamma' = 0": 1.0, "gamma.gamma' = 1": rd(np.sqrt(3) / 2)}
      and len(bahn) == 128)
# D5 P_IF: zweite Verklebung
M8b = Modell(1)
spieg = lambda v: tuple(v[:7]) + (-v[7],)
perm = [M8b.idx[spieg(r)] for r in R] + list(range(240, 248))
kl3b = np.array([sum(r[5:]) % 4 for r in M8b.R] + [0] * 8)
kl5b = np.array([sum(r[:5]) % 4 for r in M8b.R] + [0] * 8)
A1b = np.array([r[0] for r in M8b.R] + [0] * 8)
Wb = np.diag(1j ** (A1b ** 2 % 4))
okb = all(sum(1 for x in r if x < 0) % 2 == 1 for r in M8b.R if r[0] % 2) and all(
    (kl5b[i] - kl3b[i]) % 4 == 0 for i in range(248))
errb = 0.0
for a in probe[:8]:
    for b in M8b.R[::3]:
        a2 = spieg(a)
        kom = M8b.E(a2) @ M8b.E(b) - M8b.E(b) @ M8b.E(a2)
        s = ip(a2, b)
        soll_ = (M8b.eps(a2, b) * M8b.E(tuple(x + y for x, y in zip(a2, b))) if s == -1 else
                 M8b.eps(a2, b) * M8b.cartan(a2) if s == -2 else 0)
        errb = max(errb, nrm(kom - soll_))
Gb = gatter(M8b, spieg(g0), Wb, A1b.astype(float))
Pm_ = np.zeros((248, 248))
for i, j in enumerate(perm):
    Pm_[j, i] = 1.0 if i < 247 else -1.0               # Spiegelung: alpha^8 -> -alpha^8
nichtlokal = {sum(x * y for x, y in zip(g, spieg(h))) % 4 for g in kleb[::9] for h in kleb[::11]}
check("[E] P_IF: Die zweite Verklebung (Spinorklasse mit ungerader Zahl von Minuszeichen, a = a', Lagrangesche <(1,1)>) "
      "erfüllt E8 ebenso (Stichprobe). Die Spiegelung der achten Koordinate (Ladungskonjugation von A3) überführt die "
      "Modelle ineinander und lässt W fest: Die Gattermatrizen beider Verklebungen sind identisch. Klebewurzeln der "
      "beiden Verklebungen haben halbzahlige Skalarprodukte (gegenseitig nichtlokal)",
      okb and errb < TOL and nrm(Pm_ @ W @ Pm_.T - Wb) < TOL and nrm(Pm_ @ G0["lhs"] @ Pm_.T - Gb["lhs"]) < 1e-8
      and nrm(Pm_ @ G0["ref"] @ Pm_.T - Gb["ref"]) < 1e-8 and nichtlokal == {2})
check("[E] Zeitrichtung: Sie ist der Drehsinn, W gegen W^(-1) = e^(4 pi i L_0^A3), also die Wahl, welcher Faktor positiv "
      "dreht. Der Wechsel der Verklebung und die Ladungskonjugation p -> -p lassen W fest; die Zeitrichtung ist ein "
      "eigenes Bit neben P_IF",
      nrm(dag(W) - np.diag(np.exp(4j * np.pi * L3))) < TOL and nrm(W - np.diag(1j ** ((-A1) ** 2 % 4))) < TOL)
befund("ersetzte Prämissen", "keine",
       "Teil D: P_uhr wird zu P_rel (Wahl von k = +-2 in der Z_8 der vollen Drehungen) und bleibt eine Wahl; P_null "
       "folgt wie bisher aus der Klassenfunktion; P_univ ist in der Faktorlesart eingebaut (wie unter P_TP), die "
       "su(2)-Lesart braucht zusätzlich P_N4; P_sym ist als Cayley-Transformierte benannt, nicht erzwungen; P_ord und P_IF unverändert.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil E: Dynamik. Welche Gruppe wird erzeugt?")
# E1: Halbdrehungen und W: endlich
mono, werte_n = True, set()
for g in kleb[::8]:
    G = gatter(M8, g, Wd, A1.astype(float)) if g != g0 else G0
    n_ = G["B"] @ G["B"]
    blk = n_[:240, :240]
    mono &= all(int(np.sum(np.abs(blk[:, j]) > 1e-9)) == 1 for j in range(240)) and nrm(n_[240:, :240]) < TOL
    werte_n |= {complex(np.round(v, 6)) for v in blk[np.abs(blk) > 1e-9]}
    mono &= nrm(n_[240:, 240:] - (np.eye(8) - np.outer(g, g) / 4)) < TOL
    kk_ = W @ n_ @ dag(W) @ dag(n_)
    mono &= nrm(kk_ - np.diag(np.diag(kk_))) < TOL and all(abs(v ** 4 - 1) < 1e-8 for v in np.diag(kk_))
check("[E] Halbdrehungen sind vom Clifford-Typ: B_gamma^2 wirkt auf den 240 Wurzelvektoren monomial mit Einträgen in "
      "{+-1, +-i} und auf dem Cartan-Teil als Spiegelung s_gamma; W B_gamma^2 W^(-1) B_gamma^(-2) ist diagonal mit "
      "vierten Einheitswurzeln. Die von den B_gamma^2, dem 2-Torsionstorus und W erzeugte Gruppe ist endlich",
      mono and werte_n <= {1, -1, 1j, -1j})
# E2: Vierteldrehungen: unendlich
g1 = paar[-1]
G1 = gatter(M8, g1, Wd, A1.astype(float))


def zerlege(jx):
    """B = Rr + sqrt(2) Ss mit Einträgen in Q(i)."""
    pw = potenzen(jx, 5)
    Rr = fpoly(jx, lambda x: {0.0: 1, 1.0: -1j, -1.0: 1j}.get(x, 0), pw=pw)
    Ss = fpoly(jx, lambda x: {0.5: (1 - 1j) / 2, -0.5: (1 + 1j) / 2}.get(x, 0), pw=pw)
    return Rr, Ss


Ra, Sa = zerlege(G0["jx"])
Rb, Sb = zerlege(G1["jx"])
tr0 = np.trace(Ra @ Rb) + 2 * np.trace(Sa @ Sb)
tr1 = np.trace(Ra @ Sb) + np.trace(Sa @ Rb)
f0, f1 = Fr(float(tr0.real)).limit_denominator(1000), Fr(float(tr1.real)).limit_denominator(1000)
ZAHLEN["Spur von B_gamma B_gamma' auf den 248 Zuständen"] = f"{f0} + {f1} sqrt 2"
check("[E] Vierteldrehungen erzeugen eine unendliche Gruppe: Für zwei Klebewurzeln mit gamma.gamma' = -1 ist "
      "Spur(B_gamma B_gamma') = 425/4 + 55 sqrt 2 (Zerlegung B = R + sqrt 2 S mit Einträgen in Q(i)). Das ist keine ganze "
      "algebraische Zahl, also hat B_gamma B_gamma' unendliche Ordnung",
      nrm(G0["B"] - Ra - np.sqrt(2) * Sa) < TOL and abs(tr0.imag) + abs(tr1.imag) < 1e-8 and f0 == Fr(425, 4)
      and f1 == 55 and abs(float(f0) - tr0.real) < 1e-8 and abs(np.trace(G0["B"] @ G1["B"]) - (float(f0) + float(f1)
                                                                                                * np.sqrt(2))) < 1e-8)
# E3: Kommutant, invariante Zerlegung, Lie-Hülle
halb = [g for g in kleb if g[0] > 0]
etik = {}
for i, r in enumerate(R):
    etik.setdefault(tuple(ip(g, r) % 2 for g in halb), []).append(i)
etik.setdefault(tuple([0] * len(halb)), []).extend(range(240, 248))
bloecke = sorted(etik.values())
groessen = sorted(len(b) for b in bloecke)
unb = [(a, b) for blk in bloecke for a in blk for b in blk]
upos = {u: i for i, u in enumerate(unb)}
Gram = np.zeros((len(unb), len(unb)), complex)
ok_t, Bs = True, {}
for g in halb:
    jx_, _, _ = M8.su2(g)
    Bg = fpoly(jx_, lambda x: np.exp(-1j * np.pi * x / 2))
    Bg[np.abs(Bg) < 1e-12] = 0
    Bs[g] = Bg
    ok_t &= nrm(np.linalg.matrix_power(Bg, 4) - np.diag([(-1.0) ** ip(g, r) for r in R] + [1.0] * 8)) < 1e-8
    zeil, spal, wert = [], [], []
    nzr = [np.nonzero(Bg[b, :])[0] for b in range(248)]
    nzc = [np.nonzero(Bg[:, a])[0] for a in range(248)]
    for (a, b), u in upos.items():
        for c_ in nzr[b]:
            zeil.append(a * 248 + c_)
            spal.append(u)
            wert.append(Bg[b, c_])
        for r_ in nzc[a]:
            zeil.append(r_ * 248 + b)
            spal.append(u)
            wert.append(-Bg[r_, a])
    Mg = sps.csr_matrix((wert, (zeil, spal)), shape=(248 * 248, len(unb)))
    Gram += (Mg.conj().T @ Mg).toarray()
ewg = np.linalg.eigvalsh(Gram)
kern = int(np.sum(ewg < 1e-8))


def zust(r, vz):
    """Zustand |r> + vz |-r> (Bild von E_r + vz E_(-r) unter Operator -> Zustand)."""
    v = np.zeros(248, complex)
    v[M8.idx[r]] = 1
    v[M8.idx[tuple(-x for x in r)]] = vz
    return v




def rang_gf2(zeilen):
    m = [list(z) for z in zeilen]
    rang, sp = 0, 0
    while rang < len(m) and sp < len(m[0]):
        piv = next((i for i in range(rang, len(m)) if m[i][sp]), None)
        if piv is None:
            sp += 1
            continue
        m[rang], m[piv] = m[piv], m[rang]
        for i in range(len(m)):
            if i != rang and m[i][sp]:
                m[i] = [(x + y) % 2 for x, y in zip(m[i], m[rang])]
        rang += 1
        sp += 1
    return rang


tors = [[ip(g, r) % 2 for r in R] for g in halb]       # Exponenten von B_gamma^4 auf den Wurzelzuständen
rg_t = rang_gf2(tors)
rg_s = rang_gf2(tors + [[r[0] % 2 for r in R]])        # mit sigma
tH = np.array([1j ** r[0] for r in R] + [1.0] * 8)                                     # e^(i pi H_1), Quadrat sigma
vk = [zust(r, 1) for r in gerade if r > tuple(-x for x in r)] + [1j * zust(g, -1) for g in halb]
Qk, _ = np.linalg.qr(np.array(vk).T)
P120 = Qk @ dag(Qk)
inv120 = max(nrm(Bs[g] @ P120 - P120 @ Bs[g]) for g in halb)
Yg = {g: 1j * (M8.E(g) - M8.E(tuple(-x for x in g))) for g in halb}                    # log B_gamma bis auf Faktor
vg = {g: 1j * zust(g, -1) for g in halb}
ok_theta = all(nrm(sig @ Yg[g].conj() @ sig - Yg[g]) < TOL and nrm(Yg[g] + dag(Yg[g])) < TOL
               and nrm(Bs[g] - expm(-np.pi / 4 * Yg[g])) < 1e-8 for g in halb[::7])
Ymat = np.array([Yg[g] for g in halb])
Vmat = np.array([vg[g] for g in halb]).T


def reeller_rang(vs):
    m = np.array(vs)
    return int(np.linalg.matrix_rank(np.concatenate([m.real, m.imag], axis=1), tol=1e-8))


kl1 = [Ymat[i] @ Vmat for i in range(len(halb))]                                       # Zustände der Klammern [Y, Y']
alle = list(Vmat.T) + [k_[:, j] for k_ in kl1 for j in range(len(halb))]
rang_k = reeller_rang(alle)
in_k = max(nrm(v - P120 @ v) for v in alle)
Yz = np.array([zc @ Yg[g] @ dag(zc) for g in halb])                                    # mit der Uhr z gedreht
Vz = zc @ Vmat
V2 = np.concatenate([Vmat, Vz], axis=1)
alle_z = list(V2.T) + [(Y_ @ V2)[:, j] for Y_ in list(Ymat) + list(Yz) for j in range(V2.shape[1])]
rang_z = reeller_rang(alle_z)
ZAHLEN["Kommutant der B_gamma auf den 248 Zuständen"] = {"Unbekannte nach 2-Torsion": len(unb), "Dimension": kern,
                                                         "drittkleinster Eigenwert": rd(ewg[2], 4)}
ZAHLEN["Lie-Hülle"] = {"Vierteldrehungen allein": rang_k, "mit der Uhr z": rang_z}
check("[E] Die B_gamma wirken auf den 248 Zuständen nicht irreduzibel: B_gamma^4 = e^(i pi gamma.H) erzeugt den "
      "2-Torsionstorus (ein Achterblock, 120 Zweierblöcke, 544 Unbekannte), und der Kommutant der 64 Vierteldrehungen "
      "hat darin Dimension 2. Die Zerlegung ist 248 = 120 + 128; der 120er Teil ist k' = span{E_b + E_(-b) (b gerade), "
      "i (E_g - E_(-g)) (g Klebewurzel)}", ok_t and groessen == [2] * 120 + [8] and len(unb) == 544 and kern == 2
      and ewg[2] > 1e-3 and inv120 < 1e-8 and int(round(np.trace(P120).real)) == 120 and len(vk) == 120)
check("[E] Die Gruppe der Vierteldrehungen liegt dicht in einer Gruppe vom Typ Spin(16)/Z_2 (Dimension 120): Die "
      "Erzeuger log B_gamma = -(pi/4) i (E_g - E_(-g)) liegen in der Fixalgebra k' der Involution X -> sigma X-quer "
      "sigma, ihre Lie-Hülle hat Rang 120, und k' ist unter der Gruppe irreduzibel. k' ist die maximal kompakte "
      "Unteralgebra so(16) der zerfallenden Form E8(8), gedreht mit der Toruswurzel e^(i pi H_1) von sigma; sie ist nicht "
      "die Fixalgebra von sigma (Schnitt: 56 = so(8) + so(8))",
      ok_theta and rang_k == 120 and in_k < 1e-8 and kern == 2
      and sum(1 for r in gerade if r > tuple(-x for x in r)) == 56
      and max(nrm(tH * zust(r, 1) - P120 @ (tH * zust(r, 1))) for r in R if r > tuple(-x for x in r)) < 1e-8)
check("[E] Nimmt man die Uhr z = e^(2 pi i lambda_s.H) (Ordnung 4, Zentrum von Spin(10)) hinzu, hat die Lie-Hülle Rang "
      "248: Vierteldrehungen und Uhr erzeugen eine dichte Untergruppe der kompakten E8. Die Involution sigma liegt "
      "schon in der Gruppe der Vierteldrehungen (im 2-Torsionstorus der B_gamma^4, Rang 8 über F_2), W = (1+i)/2 + (1-i)/2 sigma erhält deshalb die Zerlegung 120 + 128",
      rang_z == 248 and nrm(W @ P120 - P120 @ W) < 1e-8 and rg_t == 8 and rg_s == 8)
# E4: Gatter selbst
mu = next(r for r in R if ip(r, g0) == 1 and ip(r, g1) == 0)
tri = [M8.idx[mu], M8.idx[tuple(x - y for x, y in zip(mu, g0))],
       M8.idx[tuple(x - y - z for x, y, z in zip(mu, g0, g1))]]
ander = [i for i in range(248) if i not in tri]
Ua, Ub = G0["lhs"], G1["lhs"]
inv3 = nrm(Ua[np.ix_(ander, tri)]) + nrm(Ub[np.ix_(ander, tri)])
ua, ub = Ua[np.ix_(tri, tri)], Ub[np.ix_(tri, tri)]
tr3 = np.trace(ua @ ub)
zopf = rd(np.linalg.norm(ua @ ub @ ua - ub @ ua @ ub, 2))
ZAHLEN["Dreierraum: Spur von U_gamma U_gamma', Zopfdefekt"] = [str(np.round(tr3, 6)), zopf]
check("[E] Auch die Gatter erzeugen eine unendliche Gruppe: Auf einem gemeinsamen invarianten Dreierraum (Triplett der "
      "su(3) von gamma und gamma') ist Spur(U_gamma U_gamma') = sigma(mu) + 3i/2, keine ganze Gaußsche Zahl. Der "
      "Zopfdefekt dort ist sqrt(3/8) = 0,612372, derselbe Wert wie für das z120-Gatter in z144 (Modell A)",
      inv3 < 1e-8 and abs(tr3.imag - 1.5) < 1e-8 and abs(abs(tr3.real) - 1) < 1e-8
      and abs(zopf - np.sqrt(3 / 8)) < 1e-6)
# E5: Netzstruktur
summen = {}
for g in kleb[::4]:
    for h in kleb:
        if ip(g, h) == -1:
            w = tuple(x + y for x, y in zip(g, h))
            summen[w[0] % 2] = summen.get(w[0] % 2, 0) + 1
komm_g = max(nrm(W @ fpoly(M8.su2(a)[0], lambda x: np.exp(-1j * np.pi * x / 2))
                 - fpoly(M8.su2(a)[0], lambda x: np.exp(-1j * np.pi * x / 2)) @ W) for a in gerade[::10])
check("[E] Netzstruktur: Die Summe zweier Klebewurzeln ist nie eine Klebewurzel (sigma-Gradierung; das Transfernetz ist "
      "zweigeteilt), und für gerade Wurzeln ist B^+ W B = W. Alle 64 Gatter wirken auf demselben Sektor (Zentralladung 8); E8 "
      "liefert keine Tensorstruktur, keine Kette und keine Reihenfolge", set(summen) == {0} and komm_g < TOL)
Yw = W @ Yg[g0] @ dag(W)
kw = [np.sum(x * Yw) / n for x, n in zip(gen, gnorm)]
rest_w = np.sqrt(max(0.0, nrm(Yw) ** 2 - sum(abs(c) ** 2 * n for c, n in zip(kw, gnorm)))) / nrm(Yw)
check("[E/O] Mit W verlässt die erzeugte Gruppe E8: W log(B_gamma) W^(-1) liegt nicht in e8 (orthogonal zu e8 in der Spurform). Die "
      "Hülle bleibt blockdiagonal in 120 + 128; welche Gruppe es genau ist, ist hier nicht bestimmt (offen). In jedem "
      "Fall ist sie nicht endlich und zeichnet keine Dynamik aus", abs(rest_w - 1) < 1e-8 and nrm(W @ P120 - P120 @ W) < 1e-8)
befund("Dynamik aus der Lesart", "nein",
       "Teil E: Mit Halbdrehungen ist alles vom Clifford-Typ (endlich). Mit den Vierteldrehungen, die das Gatter "
       "braucht, liegt die Transfergruppe dicht in einer Gruppe vom Typ Spin(16)/Z_2, mit der Uhr in E8, und mit W "
       "verlässt sie E8. Das ist eine dichte Gruppe auf einem einzigen Sektor, keine Kette und kein Netz. K15 wird "
       "nicht berührt.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil F: Einordnung, Bilanz und Urteil")


def lagrange(c1, c2):
    q = lambda x, y: Fr(c1 * x * x + c2 * y * y, 8)
    iso = [(x, y) for x in range(4) for y in range(4) if q(x, y).denominator == 1]
    lag = 0
    for gx, gy in iso:
        orb = sorted({((k * gx) % 4, (k * gy) % 4) for k in range(4)})
        if len(orb) == 4 and all(o in iso for o in orb) and (gx, gy) == orb[1]:
            lag += 1
    unter = {tuple(sorted({((k * gx) % 4, (k * gy) % 4) for k in range(4)})) for gx, gy in iso}
    return len(iso), len({u for u in unter if len(u) == 4 and all(o in iso for o in u)})


lag = {"D5 + A3": lagrange(5, 3), "D5 + D5": lagrange(5, 5), "A3 + A3": lagrange(3, 3)}
ZAHLEN["isotrope Elemente und Lagrangesche Untergruppen"] = {k: list(v) for k, v in lag.items()}
check("[E] Zwei gleichartige Faktoren verkleben nicht: D5 + D5 (und A3 + A3) hat nur zwei isotrope Elemente und keine "
      "Lagrangesche Untergruppe; D5 + A3 hat sechs und zwei. Die Verklebung braucht entgegengesetzte Formen, wie Naht "
      "und Tomita-Partner in P_TP (z130), nicht zwei gleiche Nähte (z99, z120)",
      lag == {"D5 + A3": (6, 2), "D5 + D5": (2, 0), "A3 + A3": (2, 0)})
check("[E/L] Entsprechung zu z130: Die gemeinsame Drehung T5 T3 = e^(2 pi i (L_0^D5 + L_0^A3)) = 1 entspricht der "
      "nativen 2 pi-Drehung der verdoppelten Naht (R = R-quer, Niveauabgleich); das Gatter ist die Gegendrehung "
      "T5 T3^(-1). Nativ ist die Komposition also trivial (wie z121), das Gatter verlangt den umgekehrten Drehsinn eines "
      "Faktors", nrm(T5 @ T3 - Id) < TOL and nrm(T5 @ dag(T3) - W) < TOL)

bilanz = [
    ("E", "Modell",
     "Niveau L_0 = 1 der Stufe-1-Darstellung (248 Zustände) aus Vertexoperatoren mit Kozykel; E8-Relationen an 9600 "
     "Kommutatoren; Kettenmodell mit Oszillatoren für Spins bis 5/2."),
    ("E", "Was die 4 pi-Drehung ist",
     "W = e^(4 pi i L_0^D5) = (1+i)/2 + (1-i)/2 sigma: Phasenwurzel des inneren Automorphismus sigma = e^(2 pi i H_1) "
     "(Ordnung 2, Fixalgebra so(16)). Kein Automorphismus; W E_gamma W^(-1) = i E_gamma sigma für alle 128 Klebeströme."),
    ("E", "Wörtliche Prüfung",
     "B_gamma^+ W B_gamma = U(-pi/2) Gamma_y^(s N/2). Gleich dem z120-Gatter (a + b Tausch, Sektorphase 1) genau für "
     "N = 0 mod 4; N = 2 mod 4: umgekehrte Zeitrichtung; N ungerade: zusätzlich die einkörperige Wurzel des Tauschs."),
    ("E", "2 pi-Drehung allein",
     "T5 = e^(2 pi i L_0^D5) ist die achte Wurzel: Teiltausch mit Gewichten sin^2(pi/8) und cos^2(pi/8), Uhrwert z der "
     "Ordnung 4. Sie ist die gehobene Uhr; das Gatter ist ihr Quadrat."),
    ("E", "Zwei gleichartige Ladungen",
     "Das z120-Gatter ohne Abweichung ist B_gamma^+ T5 T5' B_gamma mit dem gespiegelten D5' = s_gamma(D5). D5 und D5' "
     "teilen ein su(5); es sind keine zwei unabhängigen Nähte."),
    ("E", "Prämissen",
     "P_uhr -> P_rel (Wahl in Z_8, sechs wechselwirkende statt zwei): unverändert als Wahl. P_null: folgt aus der "
     "Klassenfunktion (wie z144). P_univ: in der Faktorlesart eingebaut (wie unter P_TP); die su(2)-Lesart braucht P_N4. P_sym: Cayley-Transformierte, "
     "nicht erzwungen. P_ord: unverändert. P_IF: unverändert, die Zeitrichtung ist ein eigenes Bit."),
    ("E/O", "Dynamik",
     "Halbdrehungen: endlich (Clifford-Typ). Vierteldrehungen: dicht in Spin(16)/Z_2 (so(16) der zerfallenden Form "
     "E8(8), Dimension 120); mit der Uhr z dicht in E8; mit W außerhalb von E8 (Hülle offen). Kein fester "
     "Automorphismus, aber auch keine Kette: ein Sektor, zweigeteiltes Transfernetz, keine Reihenfolge."),
    ("L", "Einordnung",
     "Nativ: E8 mit D5 + A3, die Formen, sigma, die Drehungen als Operatoren. Lesart: P_E8 (die komponierenden Objekte "
     "sind die Verklebungsfaktoren, der Transfer ist ein Klebestrom) und P_rel. Sie passt zur Partnerlesart P_TP, nicht "
     "zu zwei gleichen Nähten."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
urteil = ("Trägt bedingt (auf Multipletts mit N = 0 mod 4, unter den Prämissen P_E8 und P_rel); als Gewinn an "
          "Prämissen oder Dynamik gescopt negativ. (1) Auf Operatorebene ist die 4 pi-Drehung des D5-Faktors W = e^(4 pi i L_0^D5) "
          "= (1+i)/2 + (1-i)/2 sigma, die Phasenwurzel des inneren Automorphismus sigma = e^(2 pi i H_1) von E8 "
          "(Ordnung 2, Fixalgebra so(16)). Sie ist kein Automorphismus der Ströme, erfüllt aber die Clifford-Form von "
          "P_uhr für alle 128 Klebeströme: W^(-1) E_gamma W = -i E_gamma sigma. (2) Konjugiert mit der Vierteldrehung "
          "der Klebe-su(2) gilt exakt B_gamma^+ W B_gamma = e^(i pi (M + s J_y)^2/2) = U(-pi/2) Gamma_y^(s N/2) mit "
          "N = 2M = 4 omega_senkrecht.H. Das ist das z120-Gatter in der Form a + b mal Tausch mit Sektorphase 1 auf "
          "Multipletts mit N = 0 mod 4 (das Triplett des Niveaus L_0 = 1, alle Spins im Kettenmodell; auf Singuletts "
          "gilt die Gleichheit immer). "
          "Bei N = 2 mod 4 ist es das Gatter in umgekehrter Zeitrichtung, bei ungeradem N (die 56 Dubletts) kommt die "
          "einkörperige Wurzel des Tauschs hinzu. Alle vier Werte von N mod 4 kommen in E8 vor. Die 4 pi-Drehung ist "
          "die Paarphase einer einzigen zirkularen Mode, i^(Q_+^2), nicht beider. (3) Die 2 pi-Drehung von D5 allein "
          "ist die gehobene Uhr: ein ungleich gewichteter Teiltausch mit achten Wurzeln, dessen Quadrat das Gatter "
          "ist; auch sie ist eine Funktion der Uhrwerte. (4) Prämissen: Keine entfällt. P_uhr wird zur Wahl P_rel in "
          "der Z_8 der vollen Drehungen, mit sechs statt zwei wechselwirkenden Elementen; P_null folgt wie bisher; "
          "P_univ ist in der Faktorlesart die Verklebung selbst (wie unter P_TP, kein neuer Gewinn), die su(2)-Lesart "
          "braucht zusätzlich P_N4 (Gesamtladung N = 0 mod 4); P_sym ist die Cayley-Transformierte und gleichwertig zur Leckformel, aber von E8 nicht erzwungen; "
          "P_ord bleibt; P_IF ändert das Gatter nicht, die Zeitrichtung ist ein eigenes Bit (Drehsinn). Neu nötig "
          "sind P_E8 (die komponierenden Objekte sind die Verklebungsfaktoren, der Transfer ist ein Klebestrom), P_rel "
          "und für die wörtliche Gleichheit P_N4. (5) Dynamik: Mit "
          "Halbdrehungen ist die erzeugte Gruppe endlich (Clifford-Typ). Die Vierteldrehungen erzeugen eine dichte "
          "Untergruppe einer Gruppe vom Typ Spin(16)/Z_2 (die so(16) der zerfallenden Form E8(8), nicht die Fixalgebra "
          "von sigma); mit der Uhr z wird daraus E8, mit W verlässt die Gruppe E8. Das ist eine dichte Gruppe auf "
          "einem einzigen Sektor mit zweigeteiltem Transfernetz, keine Kette und kein Netz. "
          "(6) Die Lesart gehört zur Partnerlesart P_TP: Zwei gleichartige Faktoren haben keine Lagrangesche "
          "Untergruppe, die Verklebung verlangt entgegengesetzte Formen. Das Gatter zweier gleicher Ladungen gibt es "
          "in E8 nur als Produkt der 2 pi-Drehungen von D5 und seinem Spiegelbild, die ein su(5) teilen. Korrektur zu "
          "z145 Teil D: 'Niveauabgleich N = 0 mod 4 von selbst' gilt für die Klassen der Faktoren, nicht für die "
          "Multipletts der Klebe-su(2). P_uhr, P_ord, P_IF, K13 und K15 bleiben, wie sie sind. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z149_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z149.txt").write_text("\n".join(OUT) + "\n")
