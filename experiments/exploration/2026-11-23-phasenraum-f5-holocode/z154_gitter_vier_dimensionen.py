"""Z154: Ist die 16 der Naht die Kähler-Dirac-Zelle eines vierdimensionalen gestaffelten Gitters? Folgt daraus d = 4?

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart oder Literatur, N numerisch, O offen.

Anlass: z152 liest die Naht als Hälfte einer 1+1-dimensionalen Kogut-Susskind-Kette (Prämisse P_CTM). Offen ist K12:
Warum vier Dimensionen? Gestaffelte Fermionen in d Dimensionen sind Kähler-Dirac-Felder mit 2^d Komponenten je Zelle,
in d = 4 also 16. Frage dieser Rechnung: Ist die 16 von Spin(10) der Naht diese Zelle, und folgt d = 4 nativ?

Bestand: Die 16 ist Lambda^gerade(C^5) auf fünf Trägerplätzen (3 Farbe, 2 schwach; Y = (1 - 5 R^2)/12, z82). Die Uhr
ist das Zentrum Z_4 von Spin(10), i auf der ganzen 16 (X = 5(B-L) - 4Y = 1 mod 4). Die Spin-Z_4-Anomalie zählt
nu = n_1 - n_3 mod 16 über linkshändige Weyl-Fermionen (H6h5/H6h9). Träger mit Lorentz-Index aus z150: T16 = Twistor x 16
und der Indexkern K = (2_L x 16) + (2_R x 16-quer). Native Einteilchenzählung der 4D-Leiter: Summe m^2 q^m (z73).

Literatur [L], nur über Kurzfassung und Auszug gelesen, Gleichungsnummern nicht einzeln gegengeprüft:
  Butt/Catterall/Pradhan/Toga arXiv:2101.01026 (U(1) der Kähler-Dirac-Felder, gemischte Gravitationsanomalie ~ Euler-Zahl,
  Rest Z_4; zwei Kähler-Dirac-Felder = 8 Dirac = 16 Majorana); Catterall arXiv:2209.03828 (dasselbe für gestaffelte
  Felder); Catterall arXiv:2311.02487 (reduzierte Kähler-Dirac-Felder, Bedingung Summe q = 0 mod 4, vier reduzierte
  Felder, sechzehn linkshändige Weyl-Fermionen mit Pati-Salam-Inhalt); Peschel/Chung und Eisler/Peschel für Halbketten.

Teil 1  Darstellungstheorie, exakt: (a) 16 = Lambda*(C^4) unter SU(4) in SU(5), Formeninhalt unter SO(4);
        (b) Pati-Salam (2,1,4) + (1,2,4-quer) und der Bispinor unter der diagonalen Spin(4);
        (c) Verträglichkeit mit dem Standardmodell; welche Dimension d eine Verdrillung überhaupt zulässt.
Teil 2  Die Z_4: Kähler-Dirac-Z_4 als Element von Spin(10), Vergleich mit der Uhr, Anomaliezahlen.
Teil 3  Halbraum eines hyperkubischen gestaffelten Gitters: Kontinuum von Leiterabständen für D > 1.
Teil 4  Die native 4D-Leiter Summe m^2 q^m als Gitter-Ecke? Zählidentität, Symmetrie, drei Kandidaten.
Teil 5  Bilanz und Urteil.
"""
import itertools
import json
import os
import time
from collections import Counter
from fractions import Fraction
from math import comb, factorial
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # deterministische Rundung (einfädig)

import mpmath as mp
import numpy as np
import sympy as sp
from flint import fmpz_mat
from scipy.special import ellipk

mp.mp.dps = 40
HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val, text):
    BEFUNDE[name] = val
    say(f"   [Befund: {val}] {text}")


def rd(x, n=6):
    return round(float(x), n) + 0.0


# =====================================================================================================================
say("Z154: Ist die 16 der Naht die Kähler-Dirac-Zelle eines vierdimensionalen gestaffelten Gitters?")
say("=" * 100)
say("\nTeil 1  Darstellungstheorie (exakt, ganzzahlige Matrizen auf dem Fockraum der fünf Trägerplätze)")
say("   Plätze 1, 2, 3 = Farbe, Plätze 4, 5 = schwach; 16 = gerader Teil des Fockraums")

NM = 5


def popc(s):
    return bin(s).count("1")


def a_dag(i):
    M = np.zeros((32, 32), dtype=np.int64)
    for s in range(32):
        if not (s >> i) & 1:
            M[s | (1 << i), s] = (-1) ** popc(s & ((1 << i) - 1))
    return M


AD = [a_dag(i) for i in range(NM)]
AN = [m.T.copy() for m in AD]
I32 = np.eye(32, dtype=np.int64)
EVEN = [s for s in range(32) if popc(s) % 2 == 0]
ODD = [s for s in range(32) if popc(s) % 2 == 1]
I16 = np.eye(16, dtype=np.int64)


def r16(M):
    return M[np.ix_(EVEN, EVEN)]


def E(i, j):
    return AD[i] @ AN[j]


def NN(i):
    return E(i, i)


def comm(A, B):
    return A @ B - B @ A


def rank(M):
    M = np.asarray(M, dtype=np.int64)
    if M.size == 0:
        return 0
    M = np.unique(M, axis=0)
    M = M[np.any(M != 0, axis=1)]
    if M.shape[0] == 0:
        return 0
    return fmpz_mat(M.shape[0], M.shape[1], [int(v) for v in M.ravel()]).rank()


def joint_dim(pairs):
    """Dimension des gemeinsamen Eigenraums {O = v} auf der 16 (exakt über den Rang)."""
    return 16 - rank(np.vstack([r16(O) - v * I16 for O, v in pairs]))


def commutant_dim(gens):
    rows = [np.kron(r16(G), I16) - np.kron(I16, r16(G).T) for G in gens]
    return 256 - rank(np.vstack(rows))


def span_dim(gens):
    return rank(np.array([r16(G).ravel() for G in gens]))


def centraliser_dim(algebra, twist):
    """Dimension von {X in span(algebra): [X, T] = 0 für alle T in twist} (auf der 16, die treu ist)."""
    cols = [np.concatenate([r16(comm(S, T)).ravel() for T in twist]) for S in algebra]
    return len(algebra) - rank(np.array(cols).T)


def intersection_dim(g1, g2):
    return span_dim(g1) + span_dim(g2) - span_dim(list(g1) + list(g2))


ok = all(np.array_equal(AN[i] @ AD[j] + AD[j] @ AN[i], (i == j) * I32) for i in range(NM) for j in range(NM))
ok &= all(not np.any(AD[i] @ AD[j] + AD[j] @ AD[i]) for i in range(NM) for j in range(NM))
check("1.0 Fockraum: {a_i, a_j^+} = delta_ij, {a_i^+, a_j^+} = 0 (fünf Plätze, 32 Zustände, 16 gerade)", ok)

# --- Algebren (komplexifizierte Erzeuger, ganzzahlig) ------------------------------------------------------------------
SO10 = ([E(i, j) for i in range(5) for j in range(5) if i != j] + [2 * NN(i) - I32 for i in range(5)]
        + [AD[i] @ AD[j] for i in range(5) for j in range(i + 1, 5)]
        + [AN[i] @ AN[j] for i in range(5) for j in range(i + 1, 5)])
SU5 = [E(i, j) for i in range(5) for j in range(5) if i != j] + [NN(i) - NN(i + 1) for i in range(4)]
SU3 = [E(i, j) for i in range(3) for j in range(3) if i != j] + [NN(0) - NN(1), NN(1) - NN(2)]
SL_L = (E(3, 4), E(4, 3), NN(3) - NN(4))                                   # SU(2)_L auf den schwachen Plätzen
SL_R = (AD[3] @ AD[4], AN[4] @ AN[3], NN(3) + NN(4) - I32)                 # SU(2)_R: Paarerzeugung der schwachen Plätze
Y6 = -2 * (NN(0) + NN(1) + NN(2)) + 3 * (NN(3) + NN(4))                    # 6 Y
SM = SU3 + list(SL_L) + [Y6]
SO6 = ([E(i, j) for i in range(3) for j in range(3) if i != j] + [2 * NN(i) - I32 for i in range(3)]
       + [AD[i] @ AD[j] for i in range(3) for j in range(i + 1, 3)]
       + [AN[i] @ AN[j] for i in range(3) for j in range(i + 1, 3)])
SL_F1 = (E(0, 1), E(1, 0), NN(0) - NN(1))                                  # Spin(4) in Spin(6) auf den Farbplätzen 1, 2
SL_F2 = (AD[0] @ AD[1], AN[1] @ AN[0], NN(0) + NN(1) - I32)
NC = NN(0) + NN(1) + NN(2)
NW = NN(3) + NN(4)
NT = NC + NW


def is_sl2(t):
    e, f, h = t
    return (np.array_equal(comm(e, f), h) and np.array_equal(comm(h, e), 2 * e) and np.array_equal(comm(h, f), -2 * f))


def cas(t):
    """4 j (j + 1) für ein sl(2)-Tripel (e, f, h)."""
    e, f, h = t
    return h @ h + 2 * (e @ f + f @ e)


def tsum(t1, t2):
    return tuple(x + y for x, y in zip(t1, t2))


check("1.0 sl(2)-Tripel: SU(2)_L, SU(2)_R (schwach) und die zwei Faktoren von Spin(4) in Spin(6) erfüllen [e,f] = h, "
      "[h,e] = 2e; alle vier vertauschen paarweise",
      all(is_sl2(t) for t in (SL_L, SL_R, SL_F1, SL_F2))
      and all(not np.any(comm(x, y)) for t1, t2 in itertools.combinations((SL_L, SL_R, SL_F1, SL_F2), 2)
              for x in t1 for y in t2))
kd = {"so(10)": commutant_dim(SO10), "su(5)": commutant_dim(SU5), "SM": commutant_dim(SM),
      "Pati-Salam": commutant_dim(list(SL_L) + list(SL_R) + SO6)}
ZAHLEN["Kommutantendimension auf der 16"] = kd
say(f"   Kommutantendimension auf der 16 (= Summe der Quadrate der Vielfachheiten): {kd}")
check("1.0 Die 16 ist irreduzibel unter so(10); unter su(5) 3 Teile (1 + 10 + 5-quer), unter dem Standardmodell 6 "
      "(Q, u^c, d^c, L, e^c, nu^c), unter Pati-Salam 2",
      kd == {"so(10)": 1, "su(5)": 3, "SM": 6, "Pati-Salam": 2})

# --- (a) SU(4) in SU(5) auf vier Plätzen; SO(4) als reelle Drehungen dieser Plätze -------------------------------------
say("\n   (a) Spin(10) > SU(5) > SU(4) auf vier der fünf Plätze; SO(4) = reelle Drehungen dieser vier Plätze")


def so4_real(S):
    """Reelle Drehungen L_ab = E_ab - E_ba der vier Plätze S; Rückgabe: sechs Erzeuger und die zwei su(2)-Kasimire."""
    L = lambda a, b: E(S[a], S[b]) - E(S[b], S[a])
    J = [L(1, 2), L(2, 0), L(0, 1)]
    K = [L(0, 3), L(1, 3), L(2, 3)]
    A2 = [J[k] + K[k] for k in range(3)]                    # 2 A_k
    B2 = [J[k] - K[k] for k in range(3)]                    # 2 B_k
    zwei = all(not np.any(comm(a, b)) for a in A2 for b in B2)
    for X in (A2, B2):
        c = comm(X[0], X[1])
        zwei &= np.array_equal(c, 2 * X[2]) or np.array_equal(c, -2 * X[2])
    CA = -sum(x @ x for x in A2)                            # 4 j_A (j_A + 1)
    CB = -sum(x @ x for x in B2)
    return J + K, CA, CB, zwei


def inhalt(CA, CB, extra=()):
    res = {}
    for a, b in itertools.product((0, 3, 8), repeat=2):
        d = joint_dim([(CA, a), (CB, b)] + list(extra))
        if d:
            da, db = {0: 1, 3: 2, 8: 3}[a], {0: 1, 3: 2, 8: 3}[b]
            assert d % (da * db) == 0
            res[f"({da},{db})"] = d // (da * db)
    return res


S_A1 = (0, 1, 2, 3)                                         # drei Farbplätze + ein schwacher
S_A2 = (0, 1, 3, 4)                                         # zwei Farbplätze + zwei schwache
N4 = sum(NN(i) for i in S_A1)
pvert = [joint_dim([(N4, p)]) for p in range(5)]
check("1a Formengrad p = Besetzung der vier Plätze: Vielfachheiten 1, 4, 6, 4, 1; der fünfte Platz ist genau bei "
      "ungeradem p besetzt, also 16 = Lambda*(C^4) = 1 + 4 + 6 + 4-quer + 1 unter SU(4)",
      pvert == [1, 4, 6, 4, 1] and all(popc(s & 0b01111) % 2 == (s >> 4) & 1 for s in EVEN))
SU4A = [E(i, j) for i in S_A1 for j in S_A1 if i != j] + [NN(S_A1[k]) - NN(S_A1[k + 1]) for k in range(3)]
gew = {p: sorted(tuple(4 * ((s >> i) & 1) - p for i in S_A1) for s in EVEN if popc(s & 0b01111) == p) for p in (1, 3)}
check("1a Lambda^3 ist konjugiert zu Lambda^1 (spurfreie Gewichte genau entgegengesetzt); Kommutant von su(4) auf der "
      "16 hat Dimension 7 = 2^2 + 1 + 1 + 1 (zwei Singuletts, 4, 6, 4-quer)",
      gew[3] == sorted(tuple(-x for x in w) for w in gew[1]) and commutant_dim(SU4A) == 7)
SO4A, CA, CB, zwei = so4_real(S_A1)
check("1a so(4) der reellen Drehungen zerfällt in zwei vertauschende su(2) (Relationen geprüft)", zwei)
tab_a = {p: inhalt(CA, CB, [(N4, p)]) for p in range(5)}
ZAHLEN["(a) SO(4)-Inhalt je Formengrad p"] = tab_a
for p in range(5):
    say(f"      p = {p}: {tab_a[p]}")
SOLL_A = {0: {"(1,1)": 1}, 1: {"(2,2)": 1}, 2: {"(1,3)": 1, "(3,1)": 1}, 3: {"(2,2)": 1}, 4: {"(1,1)": 1}}
check("1a 16 unter SO(4) in SU(4): 1 + (2,2) + (3,1) + (1,3) + (2,2) + 1, nach Formengrad sortiert: das ist der "
      "Formeninhalt Lambda^p(R^4) eines Kähler-Dirac-Felds", tab_a == SOLL_A)
ges_a = inhalt(CA, CB)
SO4A2, CA2, CB2, zwei2 = so4_real(S_A2)
check("1a Dasselbe für die andere Platzwahl (zwei Farb- und zwei schwache Plätze)",
      zwei2 and inhalt(CA2, CB2) == ges_a == {"(1,1)": 2, "(1,3)": 1, "(2,2)": 2, "(3,1)": 1})


def su2x(a, c):
    return list(range(abs(a - c) + 1, a + c, 2))


def tens(r1, r2):
    """Tensorprodukt von Darstellungen von SU(2) x SU(2); Etiketten werden als Tupel verkettet."""
    res = Counter()
    for (a, b, l1), m1 in r1.items():
        for (c, d, l2), m2 in r2.items():
            for x in su2x(a, c):
                for y in su2x(b, d):
                    res[(x, y, l1 + l2)] += m1 * m2
    return res


def diag(rep):
    """Einschränkung von (Lorentz) x (innere Spin(4)) auf die Diagonale: rep hat Schlüssel (a, b, c, d, etikett)."""
    res = Counter()
    for (a, b, c, d, l), m in rep.items():
        for x in su2x(a, c):
            for y in su2x(b, d):
                res[(x, y, l)] += m
    return res


def ohne_etikett(rep):
    res = Counter()
    for (a, b, _), m in rep.items():
        res[(a, b)] += m
    return dict(sorted(res.items()))


S4 = Counter({(2, 1, ()): 1, (1, 2, ()): 1})               # Dirac-Spinor von Spin(4)
FORMEN = {(1, 1): 2, (1, 3): 1, (2, 2): 2, (3, 1): 1}
GERADE = {(1, 1): 2, (1, 3): 1, (3, 1): 1}
UNGERADE = {(2, 2): 2}
check("1a Bispinor: S x S mit S = (2,1) + (1,2) hat denselben Inhalt wie Lambda*(R^4)",
      ohne_etikett(tens(S4, S4)) == FORMEN
      and {f"({a},{b})": m for (a, b), m in FORMEN.items()} == ges_a)

# --- (b) Pati-Salam und die diagonale Spin(4) --------------------------------------------------------------------------
say("\n   (b) Spin(10) > Spin(4) x Spin(6): schwache Plätze tragen Spin(4) = SU(2)_L x SU(2)_R, Farbplätze Spin(6) = SU(4)")
CL, CR = cas(SL_L), cas(SL_R)
PC = I32 - 2 * (NC % 2)                                     # (-1)^(N_c): Chiralität von Spin(6)
ps = {}
for (cl, nl), (cr, nr), (pc, nm) in itertools.product(((0, 1), (3, 2)), ((0, 1), (3, 2)), ((-1, "4"), (1, "4-quer"))):
    d = joint_dim([(CL, cl), (CR, cr), (PC, pc)])
    if d:
        ps[f"({nl},{nr},{nm})"] = d
ZAHLEN["(b) Pati-Salam-Sektoren der 16 (Dimension)"] = ps
say(f"      Sektoren (Dimension): {ps}")
check("1b 16 = (2,1,4) + (1,2,4-quer) unter SU(2)_L x SU(2)_R x SU(4): zwei Sektoren der Dimension 8, irreduzibel "
      "(Kommutant 2); '4' = ungerade Farbbesetzung (Farbtriplett + Singulett, also Q und L)",
      ps == {"(2,1,4)": 8, "(1,2,4-quer)": 8})
DL, DR = tsum(SL_L, SL_F1), tsum(SL_R, SL_F2)
check("1b Spin(4) in Spin(6) über die Spinoreinbettung (Farbplätze 1, 2): 4 -> (2,1) + (1,2); die Diagonale mit der "
      "schwachen Spin(4) ist wieder ein Paar vertauschender sl(2)",
      is_sl2(DL) and is_sl2(DR) and all(not np.any(comm(x, y)) for x in DL for y in DR)
      and inhalt(cas(SL_F1), cas(SL_F2), [(PC, -1)]) == {"(1,2)": 2, "(2,1)": 2})
ges_b = inhalt(cas(DL), cas(DR))
ZAHLEN["(b) Inhalt der 16 unter der diagonalen Spin(4)"] = ges_b
say(f"      16 unter der diagonalen Spin(4): {ges_b}")
check("1b Unter der diagonalen Spin(4) ist die 16 wieder 2 (1,1) + 2 (2,2) + (3,1) + (1,3) = Lambda*(R^4) = S x S*: "
      "derselbe Inhalt wie in (a)", ges_b == ges_a)
DIAGB = list(DL) + list(DR)
check("1b (a) und (b) sind dieselbe Art Untergruppe: Kommutant auf der 16 jeweils Dimension 10 = 2^2 + 2^2 + 1 + 1",
      commutant_dim(SO4A) == 10 and commutant_dim(SO4A2) == 10 and commutant_dim(DIAGB) == 10)

# --- (c) Verträglichkeit mit dem Bestand -------------------------------------------------------------------------------
say("\n   (c) Verträglichkeit mit den Trägerplätzen 3 + 2 und dem Standardmodell")
SO4PS = list(SL_L) + list(SL_R)
vertr = {}
for name, tw in (("A1: SO(4) in SU(4) in SU(5), Plätze 1-4 (3 Farbe + 1 schwach)", SO4A),
                 ("A2: SO(4) in SU(4) in SU(5), Plätze 1, 2, 4, 5 (2 Farbe + 2 schwach)", SO4A2),
                 ("A3: diagonale Spin(4) aus (b), schwache Plätze mit zwei Farbplätzen", DIAGB),
                 ("B: schwache Spin(4) = SU(2)_L x SU(2)_R allein (Verdrillung mit äußerem Lorentz)", SO4PS)):
    vertr[name] = {"Dimension": span_dim(tw), "Schnitt mit SM": intersection_dim(SM, tw),
                   "SM-Erzeuger, die vertauschen": centraliser_dim(SM, tw),
                   "davon in su(3)": centraliser_dim(SU3, tw)}
    say(f"      {name}: {vertr[name]}")
ZAHLEN["(c) Verdrillungsgruppe gegen Standardmodell (Dimensionen)"] = vertr
check("1c Als innere Zelle (A1, A2, A3) verträgt sich die verdrillte Drehgruppe nicht mit dem Standardmodell: Höchstens "
      "einer der zwölf Eichgeneratoren vertauscht mit ihr, von der Farb-su(3) keiner; sie dreht Farbplätze in "
      "schwache Plätze (A1, A2) oder koppelt eine Farb-SU(2) starr an SU(2)_L (A3)",
      all(vertr[k]["SM-Erzeuger, die vertauschen"] <= 1 and vertr[k]["davon in su(3)"] == 0
          for k in vertr if k[0] == "A"))
kB = "B: schwache Spin(4) = SU(2)_L x SU(2)_R allein (Verdrillung mit äußerem Lorentz)"
check("1c Lesart B lässt die Farbe unberührt (alle acht su(3)-Erzeuger vertauschen), enthält aber die ganze schwache "
      "su(2)_L (Schnitt 3); die Hyperladung vertauscht nicht (sie enthält T_3R)",
      vertr[kB] == {"Dimension": 6, "Schnitt mit SM": 3, "SM-Erzeuger, die vertauschen": 8, "davon in su(3)": 8})

# Träger mit Lorentz-Index: psi in (2,1) x 16, psi-quer in (1,2) x 16-quer (euklidisch). Verdrillt = Diagonale.
R16 = Counter({(2, 1, ("4",)): 1, (1, 2, ("4q",)): 1})
R16Q = Counter({(2, 1, ("4q",)): 1, (1, 2, ("4",)): 1})


def lor_x(lor, rep):
    return Counter({(lor[0], lor[1], c, d, l): m for (c, d, l), m in rep.items()})


def je_etikett(rep):
    res = {}
    for (a, b, l), m in rep.items():
        res.setdefault(l[0], Counter())[(a, b)] += m
    return {k: dict(sorted(v.items())) for k, v in sorted(res.items())}


K_tw = diag(lor_x((2, 1), R16)) + diag(lor_x((1, 2), R16Q))
T_tw = diag(lor_x((2, 1), R16)) + diag(lor_x((1, 2), R16))
ZAHLEN["Indexkern K unter der verdrillten Gruppe, je SU(4)-Darstellung"] = {k: {str(a): b for a, b in v.items()}
                                                                           for k, v in je_etikett(K_tw).items()}
ZAHLEN["T16 unter der verdrillten Gruppe, je SU(4)-Darstellung"] = {k: {str(a): b for a, b in v.items()}
                                                                    for k, v in je_etikett(T_tw).items()}
say(f"      Indexkern K = (2_L x 16) + (2_R x 16-quer), verdrillt, je SU(4)-Darstellung: {je_etikett(K_tw)}")
say(f"      T16 = (2_L + 2_R) x 16, verdrillt, je SU(4)-Darstellung: {je_etikett(T_tw)}")
check("1c Lesart B, exakt: Der Indexkern K aus z150 ist unter diag(Lorentz x schwache Spin(4)) x SU(4) genau "
      "Lambda^gerade(R^4) x 4 + Lambda^ungerade(R^4) x 4-quer, also vier reduzierte Kähler-Dirac-Felder (64 = 4 mal 16)",
      je_etikett(K_tw) == {"4": GERADE, "4q": UNGERADE})
check("1c T16 = Twistor x 16 ist dagegen nicht nach Formengrad sortiert (je SU(4)-Darstellung (1,1) + (3,1) + (2,2)): "
      "kein Satz reduzierter Kähler-Dirac-Felder",
      je_etikett(T_tw)["4"] == {(1, 1): 1, (2, 2): 1, (3, 1): 1} and je_etikett(T_tw)["4"] != GERADE)
fam = ohne_etikett(tens(Counter({(2, 1, ()): 1}), S4))
check("1c Nebenrechnung K18: Wäre die Familien-SU(4) = A_3 (die 4 der Klebung (16,4)) die Geschmacksgruppe, so trüge "
      "jede linkshändige 16 den Inhalt (1,1) + (3,1) + (2,2): acht Komponenten, eine lorentzchirale Hälfte der Zelle, "
      "vier Geschmäcker statt drei Familien", fam == {(1, 1): 1, (2, 2): 1, (3, 1): 1})


# Welche Dimension d lässt eine Verdrillung zu? Spin(d) x Spin(10 - d) in Spin(10), Weyl-Spinor in d Dimensionen.
def spinor(n, chir):
    return [w for w in itertools.product((1, -1), repeat=n) if int(np.prod(w)) == chir]


W16 = spinor(5, 1)


def verdrillt(d):
    """Euklidischer Feldinhalt {psi, psi-quer} = S+ x 16 + (S-)* x 16* unter der diagonalen Spin(d); Gewichte."""
    n = d // 2
    ges, je = Counter(), {1: Counter(), -1: Counter()}
    for l in spinor(n, 1):
        for w in W16:
            v = tuple((l[i] + w[i]) // 2 for i in range(n))
            ges[v] += 1
            je[int(np.prod(w[n:]))][v] += 1
    for l in spinor(n, -1):
        for w in W16:
            v = tuple((-l[i] - w[i]) // 2 for i in range(n))
            ges[v] += 1
            je[-int(np.prod(w[n:]))][v] += 1                # 10 - d = 2 mod 4: Konjugation tauscht die Chiralität
    c = 2 ** (4 - n)
    formen, ger, ung = Counter(), Counter(), Counter()
    for v in itertools.product((-1, 0, 1), repeat=n):
        m = c * 2 ** sum(1 for x in v if x == 0)
        formen[v] = m
        (ger if sum(1 for x in v if x) % 2 == 0 else ung)[v] = m
    return c, ges == formen, (je[1] == ger and je[-1] == ung)


dtab = {}
for d in (2, 4, 6, 8):
    c, voll, sortiert = verdrillt(d)
    dtab[d] = {"Kopien c = 2^(4 - d/2)": c, "Inhalt = c mal Lambda*(R^d)": voll,
               "gerade x M, ungerade x M*": sortiert if d % 4 == 0 else False,
               "Komponenten 2^(d/2 + 4)": 2 ** (d // 2 + 4)}
    say(f"      d = {d}: {dtab[d]}")
ZAHLEN["Verdrillung von Weyl_d x 16 je Dimension d"] = dtab
check("1c Weyl-Spinor in d Dimensionen mal 16, verdrillt mit Spin(d) in Spin(10): Der Feldinhalt ist genau für d = 4 "
      "und d = 8 ein Satz reduzierter Kähler-Dirac-Felder (4 bzw. 1 Kopie), für d = 2 und d = 6 nicht (dort bleibt "
      "nur eine Lorentz-Chiralität)",
      [dtab[d]["gerade x M, ungerade x M*"] for d in (2, 4, 6, 8)] == [False, True, False, True]
      and [dtab[d]["Inhalt = c mal Lambda*(R^d)"] for d in (2, 4, 6, 8)] == [False, True, False, True])
check("1c Die native Platzteilung 3 + 2 (R^2, z82) zeichnet Spin(6) x Spin(4) aus: von beiden Faktoren kann nur "
      "Spin(4) verdrillen (d = 4, Kopiengruppe SU(4)); Spin(6) nicht (d = 6), und d = 8 verlangte die Teilung 4 + 1",
      dtab[4]["gerade x M, ungerade x M*"] and not dtab[6]["Inhalt = c mal Lambda*(R^d)"]
      and commutant_dim(list(SL_L) + list(SL_R) + SO6) == 2)
befund("1", "bedingt",
       "E: Die 16 ist auf zwei Wegen ein Bispinor Lambda*(R^4) (über SU(5) > SU(4) > SO(4) und über die Diagonale in "
       "Spin(4) x Spin(6)); beide Wege sind dieselbe Art Untergruppe. Als innere Zelle verträgt sich das nicht mit dem "
       "Bestand: Die verdrillte Drehgruppe läge in Spin(10), mischte Farb- und schwache Plätze und ließe vom "
       "Standardmodell höchstens einen Erzeuger übrig. Verträglich ist von den geprüften nur Lesart B: Der Lorentz-Index bleibt außen, "
       "die Geschmacksgruppe ist die schwache Spin(4) = SU(2)_L x SU(2)_R, die Kopiengruppe die Farb-Spin(6) = SU(4). "
       "Dann ist der Indexkern K aus z150 exakt vier reduzierte Kähler-Dirac-Felder. Das ist die Pati-Salam-Aussage "
       "der Literatur [L], hier für den Bestand nachgerechnet. Die 2^4 = 16 Komponenten einer Zelle und die 16 von "
       "Spin(10) sind dabei verschiedene Sechzehnen: Eine Zelle sind 8 Weyl-Fermionen, K sind vier Zellen (Feld und "
       "konjugiertes Feld zusammen).")

# =====================================================================================================================
say("\nTeil 2  Die Z_4: Kähler-Dirac-Z_4 = i (-1)^p gegen die Uhr (Zentrum von Spin(10))")
say("   Toruselement exp(i pi Summe t_j J_j), J_j = n_j - 1/2: auf einem Fockzustand die Phase i^(Summe t_j (2 n_j - 1))")


def bits(s):
    return [(s >> i) & 1 for i in range(5)]


def torus(t, s):
    return sum(t[j] * (2 * b - 1) for j, b in enumerate(bits(s))) % 4      # Exponent von i


T_UHR = (-1, -1, -1, -1, -1)                                # z = i^X, X = 5 - 2N
T_A = (2, 0, 0, 0, 1)                                       # Zelle (a): i (-1)^p, p = Besetzung der Plätze 1-4
T_B = (1, 1, 1, 2, 0)                                       # Lesart B: +i auf (2,1,4), -i auf (1,2,4-quer)
X = {s: 5 - 2 * popc(s) for s in EVEN}
YY = {s: Fraction(-1, 3) * popc(s & 0b00111) + Fraction(1, 2) * popc(s & 0b11000) for s in EVEN}
BL = {s: (X[s] + 4 * YY[s]) / 5 for s in EVEN}
check("2.0 Uhr: X = 5 - 2N = 1 mod 4 auf allen 16 Zuständen, also z = i^X = i; X = 5(B-L) - 4Y mit dem Spektrum "
      "B-L = 1/3 (6), -1/3 (6), -1 (2), +1 (2) und Hyperladungen Summe Y = 0",
      all(X[s] % 4 == 1 for s in EVEN) and all(torus(T_UHR, s) == 1 for s in EVEN)
      and Counter(BL.values()) == Counter({Fraction(1, 3): 6, Fraction(-1, 3): 6, Fraction(-1): 2, Fraction(1): 2})
      and sum(YY.values()) == 0)
pA = {s: popc(s & 0b01111) for s in EVEN}
check("2.1 Zelle (a): i (-1)^p ist auf der 16 das Toruselement t = (2,0,0,0,1) von Spin(10): Drehung um pi in der "
      "Ebene des fünften Platzes (mal -1)", all(torus(T_A, s) == (1 + 2 * pA[s]) % 4 for s in EVEN))
check("2.1 Lesart B: +i auf (2,1,4), -i auf (1,2,4-quer) ist das Toruselement t = (1,1,1,2,0)",
      all(torus(T_B, s) == (1 if popc(s & 0b11000) % 2 == 1 else 3) for s in EVEN))


def spur16(t):
    return sum(1j ** torus(t, s) for s in EVEN)


def spur10(t):
    return int(round(sum(2 * np.cos(np.pi * x) for x in t)))


sp_tab = {n: {"Spur auf der 16": str(complex(spur16(t))), "Spur auf der 10": spur10(t),
              "Eigenwerte i, -i": [sum(1 for s in EVEN if torus(t, s) == 1), sum(1 for s in EVEN if torus(t, s) == 3)]}
          for n, t in (("Uhr z", T_UHR), ("g_A (Zelle a)", T_A), ("g_B (Lesart B)", T_B))}
ZAHLEN["Z_4-Elemente: Spuren und Eigenwerte"] = sp_tab
for n, v in sp_tab.items():
    say(f"      {n}: {v}")
check("2.2 Weder g_A noch g_B ist die Uhr oder zu ihr konjugiert (Spur 0 gegen 16 i; die Uhr ist zentral); g_A und "
      "g_B sind auch untereinander nicht konjugiert (Spur auf der 10: 6 gegen -2)",
      sp_tab["Uhr z"]["Eigenwerte i, -i"] == [16, 0] and sp_tab["g_A (Zelle a)"]["Eigenwerte i, -i"] == [8, 8]
      and sp_tab["g_B (Lesart B)"]["Eigenwerte i, -i"] == [8, 8] and spur10(T_A) == 6 and spur10(T_B) == -2
      and spur10(T_UHR) == -10)
uA = all((torus(T_A, s) - torus(T_UHR, s)) % 4 == 2 * ((bits(s)[0] + bits(s)[1] + bits(s)[2] - 3 * bits(s)[3]) % 2)
         for s in EVEN)
uB = all((torus(T_B, s) + torus(T_UHR, s)) % 4 == 2 * ((bits(s)[3] - bits(s)[4]) % 2) for s in EVEN)
check("2.3 g_A = z mal exp(i pi (n_1 + n_2 + n_3 - 3 n_4)): Uhr mal ein Element des Standardmodell-Torus (spurfreie "
      "Diagonale)", uA)
check("2.3 g_B = z^-1 mal exp(i pi (n_4 - n_5)) = inverse Uhr mal Zentrum von SU(2)_L (gleichwertig z mal Zentrum von "
      "SU(2)_R): Als Gruppe ist die Kähler-Dirac-Z_4 die Uhr-Z_4 modulo einer Eichtransformation", uB)


# Kähler-Dirac-U(1) exp(i alpha Gamma) gegen Spin(10): Ordnung von e_alpha im Quotienten Z^6 / Gitter
def echelon(rows):
    M = [list(r) for r in rows if any(r)]
    ncol = len(M[0])
    basis = []
    for col in range(ncol):
        while True:
            nz = sorted((r for r in M if r[col] != 0), key=lambda r: (abs(r[col]), r))
            if len(nz) <= 1:
                break
            p = nz[0]
            for r in nz[1:]:
                f = r[col] // p[col]
                for c in range(ncol):
                    r[c] -= f * p[c]
        nz = [r for r in M if r[col] != 0]
        if nz:
            p = nz[0]
            M = [r for r in M if r is not p and any(r)]
            basis.append((col, p if p[col] > 0 else [-v for v in p]))
    return basis


def im_gitter(basis, v):
    v = list(v)
    for col, p in basis:
        if v[col] % p[col]:
            return False
        f = v[col] // p[col]
        v = [a - f * b for a, b in zip(v, p)]
    return not any(v)


def u1_ordnung(vorz):
    """Kleinste k mit k e_alpha im Gitter der Zeilen (2 n_j - 1, -2 Gamma(n)); None = unendlich (ganze U(1) liegt in Spin(10))."""
    rows = [[2 * b - 1 for b in bits(s)] + [-2 * vorz(s)] for s in EVEN]
    B = echelon(rows)
    for k in range(1, 65):
        if im_gitter(B, [0, 0, 0, 0, 0, k]):
            return k
    return None


kA = u1_ordnung(lambda s: 1 - 2 * (pA[s] % 2))
starr = all(-2 * (1 - 2 * (pA[s] % 2)) == 2 * (2 * bits(s)[4] - 1) for s in EVEN)   # letzte Spalte = 2 mal fünfte
kB_ = u1_ordnung(lambda s: 1 if popc(s & 0b11000) % 2 == 1 else -1)
ZAHLEN["Kähler-Dirac-U(1) in Spin(10): Ordnung k (alpha in 4 pi Z / k)"] = {"Zelle a": kA, "Lesart B": kB_}
check("2.4 Zelle (a): Die ganze Kähler-Dirac-U(1) exp(i alpha (-1)^p) liegt in Spin(10) (Drehung des fünften Platzes: "
      "(-1)^p = -2 J_5 auf der 16)", kA is None and starr)
check("2.4 Lesart B: Von der Kähler-Dirac-U(1) liegt genau die Z_4 in Spin(10) (alpha in (pi/2) Z): die Untergruppe, "
      "die nach [L] die Anomalie übrig lässt, ist genau der Teil, der in Spin(10) liegt", kB_ == 8)


def nu(t, traeger):
    """n_1 - n_3 über linkshändige Weyl-Fermionen. traeger: Liste (Chiralität +1 = links / -1 = rechts, Fockzustände)."""
    n1 = n3 = 0
    for chir, zust in traeger:
        for s in zust:
            q = torus(t, s) if chir == 1 else (-torus(t, s)) % 4
            assert q in (1, 3)
            n1 += q == 1
            n3 += q == 3
    return n1 - n3


nutab = {}
for n, t in (("Uhr z", T_UHR), ("g_A (Zelle a)", T_A), ("g_B (Lesart B)", T_B)):
    nutab[n] = {"2_L x 16 (eine Generation)": nu(t, [(1, EVEN)]),
                "T16 = (2_L + 2_R) x 16": nu(t, [(1, EVEN), (-1, EVEN)]),
                "K = 2_L x 16 + 2_R x 16-quer, als 16 Weyl-Felder": nu(t, [(1, EVEN)]),
                "Gegenprobe über den 2_R x 16-quer-Teil von K": nu(t, [(-1, ODD)])}
    say(f"      nu für {n}: {nutab[n]}")
ZAHLEN["Anomaliezahl nu = n_1 - n_3 (mod 16)"] = nutab
check("2.5 Uhr: nu = 16 auf einer Generation (alle 16 mit Ladung 1), 0 auf T16 (vektorartig); g_A und g_B: nu = 8 - 8 "
      "= 0 auf jedem Träger. Alle Werte sind 0 mod 16, die Uhr aber als 16 mal 1, die Kähler-Dirac-Z_4 als 8 - 8",
      nutab["Uhr z"]["2_L x 16 (eine Generation)"] == 16 and nutab["Uhr z"]["T16 = (2_L + 2_R) x 16"] == 0
      and all(v == 0 for n in ("g_A (Zelle a)", "g_B (Lesart B)") for v in nutab[n].values())
      and nutab["Uhr z"]["Gegenprobe über den 2_R x 16-quer-Teil von K"] == 16)
# Zelle mit Spin innen (16 Komponenten = 8 Weyl): Komponenten (s, f), s, f in {L, R} je zweifach
zelle = [(s, f) for s in (1, -1) for f in (1, -1) for _ in range(2)]
nu_zelle_g = sum((1 if (s * f == 1) == (s == 1) else -1) for s, f in zelle)       # i (-1)^p, p gerade <=> s f = 1
nu_zelle_z = sum((1 if s == 1 else -1) for s, f in zelle)                         # i auf allen Komponenten
check("2.5 Zelle mit dem Spin innen (16 Komponenten = 8 Weyl-Felder): nu = 4 - 4 = 0 für i (-1)^p und für i auf allen "
      "Komponenten; eine Zelle allein trägt keine 16", nu_zelle_g == 0 and nu_zelle_z == 0 and len(zelle) == 8)
ktab = []
for k in range(1, 9):
    nu_z, nu_g = 4 * k, 2 * k - 2 * k                        # (2,1) x k und (1,2) x k-quer, je 2k linkshändige Weyl-Felder
    inst = {"g_B, SU(2)_L": k % 4, "g_B, SU(2)_R": (3 * k) % 4, "Uhr, SU(2)_L": k % 4, "Uhr, SU(2)_R": k % 4}
    ktab.append({"k": k, "nu(Uhr)": nu_z, "nu(g_B)": nu_g, "Ladung des 't-Hooft-Vertex mod 4": inst,
                 "frei": nu_z % 16 == 0 and all(v == 0 for v in inst.values())})
ZAHLEN["k reduzierte Kähler-Dirac-Felder (2,1,k) + (1,2,k-quer)"] = ktab
say("      k reduzierte Felder: nu(Uhr) = 4k, nu(g_B) = 0, Instanton-Vertex von SU(2)_L trägt die Z_4-Ladung k mod 4")
check("2.6 Für k reduzierte Kähler-Dirac-Felder ist nu(Uhr) = 4k; die Kähler-Dirac-Z_4 hat nu = 0, aber der "
      "Instanton-Vertex der Geschmacks-SU(2) trägt die Ladung k mod 4. Beide Bedingungen sind dieselbe: k = 0 mod 4, "
      "also 16 Weyl-Felder. Das ist die Bedingung Summe q = 0 mod 4 der Literatur [L]",
      [r["k"] for r in ktab if r["frei"]] == [4, 8] and all((r["nu(Uhr)"] % 16 == 0) ==
                                                           (r["Ladung des 't-Hooft-Vertex mod 4"]["g_B, SU(2)_L"] == 0)
                                                           for r in ktab))
befund("2", "bedingt",
       "E: Als Elemente sind Kähler-Dirac-Z_4 und Uhr verschieden und nicht konjugiert (Eigenwerte 8 mal i und 8 mal -i "
       "gegen 16 mal i). In beiden Zuordnungen ist die Kähler-Dirac-Z_4 aber Uhr hoch +-1 mal ein Element des "
       "Standardmodell-Torus; in Lesart B ist dieses Element das Zentrum der schwachen SU(2)_L, also der "
       "Geschmacksgruppe selbst. Die Anomaliezahlen stimmen mod 16 immer überein (0), aber aus verschiedenem Grund: "
       "Uhr 16 mal 1, Kähler-Dirac 8 - 8. Inhaltlich gleich werden sie genau in Lesart B, wenn die Geschmacks-SU(2) "
       "geeicht (Bestand: schwache Eichgruppe) oder mit dem Spinzusammenhang verdrillt ist (Kähler-Dirac): Dann "
       "verlangen beide k = 0 mod 4 reduzierte Felder, also 16 Weyl-Fermionen.")

# =====================================================================================================================
say("\nTeil 3  Halbraum eines hyperkubischen gestaffelten Gitters")
say("   H = Summe i t eta_mu(x) (c^+_x c_(x+mu) - h.c.) + m eps(x) c^+_x c_x, eta_1 = 1, eta_2 = (-1)^x1, "
    "eta_3 = (-1)^(x1+x2), eps = (-1)^(Summe x)")


def gitter(L, m, t=1.0, tq=1.0, offen=True):
    """Gestaffeltes Fermion auf L[0] x L[1] x ...; Richtung 0 mit Hüpfen t, die übrigen mit tq."""
    dims = len(L)
    N = int(np.prod(L))
    H = np.zeros((N, N), dtype=complex)
    koord = list(itertools.product(*[range(l) for l in L]))
    pos = {x: i for i, x in enumerate(koord)}
    for x, i in pos.items():
        H[i, i] = m * (-1) ** sum(x)
        for mu in range(dims):
            y = list(x)
            y[mu] += 1
            if y[mu] >= L[mu]:
                if offen if mu == 0 else (offen is True):
                    continue
                y[mu] = 0
            j = pos[tuple(y)]
            amp = 1j * (t if mu == 0 else tq) * (-1) ** sum(x[:mu])
            H[i, j] += amp
            H[j, i] += np.conj(amp)
    return H, koord


def verschr(H, region):
    N = H.shape[0]
    _, v = np.linalg.eigh(H)
    occ = v[:, :N // 2]
    C = (occ @ occ.conj().T)[np.ix_(region, region)]
    z, U = np.linalg.eigh(C)
    z = np.clip(z, 1e-300, 1 - 1e-16)
    return np.log((1 - z) / z), U


def abstand(M, t=1.0):
    """Leiterabstand / pi der Halbkette mit gleichmäßigem Hüpfen t und Wechselpotential M (z152): 4 K(kappa')/K(kappa)."""
    kp2 = M * M / (4 * t * t + M * M)
    return 4 * float(ellipk(kp2) / ellipk(1 - kp2))


# 3.1 Dispersionsrelation in D = 3 (periodisch, 6^3): H^2 = m^2 + 4 t^2 Summe sin^2 k
H3, _ = gitter((6, 6, 6), 2.0, offen=False)
e2 = np.sort(np.linalg.eigvalsh(H3) ** 2)
soll = np.sort([4.0 + 4 * sum(np.sin(2 * np.pi * k / 6) ** 2 for k in kk) for kk in itertools.product(range(6), repeat=3)])
check("3.1 Gestaffeltes Fermion in D = 3 (6^3, periodisch): E^2 = m^2 + 4 t^2 Summe sin^2 k_mu (216 Werte, 1e-10)",
      np.allclose(e2, soll, atol=1e-10))
check("3.1 H ist hermitesch und das Spektrum symmetrisch mit Lücke 2m (halbe Füllung eindeutig)",
      np.allclose(H3, H3.conj().T) and np.allclose(np.sort(np.linalg.eigvalsh(H3)), -np.sort(np.linalg.eigvalsh(H3))[::-1])
      and abs(np.min(np.abs(np.linalg.eigvalsh(H3))) - 2.0) < 1e-10)
# 3.2 Monotonie des Abstands in der effektiven Masse (exakt über die Legendre-Relation)
ok = True
for kp in ("0.3", "0.5", "0.70710678", "0.8", "0.9"):
    kp = mp.mpf(kp)
    f = lambda x: mp.ellipk(x * x) / mp.ellipk(1 - x * x)
    abl = mp.diff(f, kp)
    geschl = mp.pi / (2 * kp * (1 - kp * kp) * mp.ellipk(1 - kp * kp) ** 2)
    ok &= abs(abl - geschl) < mp.mpf(10) ** -25 and geschl > 0
check("3.2 d/d kappa' [K(kappa')/K(kappa)] = pi / (2 kappa' kappa^2 K(kappa)^2) > 0 (fünf Punkte, 25 Stellen): Der "
      "Leiterabstand wächst streng mit der effektiven Masse, verschiedene Massen geben verschiedene Abstände", ok)
check("3.2 Kontrolle z152: M = 2t gibt den Abstand 4 pi (kappa' = 1/sqrt 2)", abs(abstand(2.0) - 4.0) < 1e-12)
# 3.3 2D-Streifen: offen in x, periodisch in y; Halbraum x < Lx/2
LX, LY, GR = 48, 8, 20.0


def streifen(m, tq, Ly=LY, Lx=LX):
    H, koord = gitter((Lx, Ly), m, tq=tq, offen="x")
    region = [i for i, x in enumerate(koord) if x[0] < Lx // 2]
    eps, U = verschr(H, region)
    return eps, U, [koord[i] for i in region]


def vorhersage(m, tq, Ly=LY, grenze=GR):
    lev, masses = [], []
    for n in range(Ly // 2):
        M = np.sqrt(m * m + 4 * tq * tq * np.sin(2 * np.pi * n / Ly) ** 2)
        masses.append(M)
        D = abstand(M) * np.pi
        j = 0
        while D * (2 * j + 1) / 4 < grenze:
            lev += [D * (2 * j + 1) / 4, -D * (2 * j + 1) / 4]
            j += 1
    return np.sort(lev), masses


str_tab = {}
ok = True
for tq in (1.0, 0.3, 0.0):
    eps, _, _ = streifen(2.0, tq)
    num = np.sort(eps[np.abs(eps) < GR])
    vor, masses = vorhersage(2.0, tq)
    gleich = len(num) == len(vor) and np.allclose(num, vor, atol=2e-5)
    ok &= gleich
    absts = sorted({rd(abstand(M), 6) for M in masses})
    str_tab[f"t_quer = {tq}"] = {"Abstände / pi": absts, "Zahl der Niveaus mit |eps| < 20": int(len(num)),
                                 "Abweichung von der Kettenformel < 2e-5": bool(gleich),
                                 "kleinste |eps| / pi": [rd(x / np.pi, 5) for x in np.sort(np.abs(num))[:8:2]]}
    say(f"      t_quer = {tq}: {str_tab[f't_quer = {tq}']}")
ZAHLEN["2D-Streifen 48 x 8, m = 2t: Halbraumspektrum"] = str_tab
check("3.3 2D-Streifen (48 x 8, m = 2t): Das Halbraumspektrum ist genau die Vereinigung der Kettenleitern mit "
      "M(k)^2 = m^2 + 4 t_quer^2 sin^2 k und Abstand 4 pi K'/K je Querimpuls (alle Niveaus |eps| < 20, 2e-5)", ok)
check("3.3 Bei t_quer = t gibt es drei verschiedene Abstände (4, 4,39 und 4,68 mal pi), bei t_quer = 0 genau einen "
      "(4 pi, achtfach): Die native Leiter erscheint nur ohne Querhüpfen",
      len(str_tab["t_quer = 1.0"]["Abstände / pi"]) == 3 and str_tab["t_quer = 0.0"]["Abstände / pi"] == [4.0]
      and len(str_tab["t_quer = 0.3"]["Abstände / pi"]) == 3)
ber = {f"D = {D}": [4.0, rd(abstand(2.0 * np.sqrt(D)), 6)] for D in (1, 2, 3, 4)}
ZAHLEN["Bereich der Abstände / pi im Halbraum, isotrop, m = 2t, D Raumdimensionen"] = ber
say(f"      Bereich der Abstände / pi (isotrop, m = 2t, unendliche Querausdehnung): {ber}")
check("3.3 Für D > 1 füllen die Abstände ein Intervall: D = 2 von 4 pi bis 4,68 pi, D = 3 bis 5,12 pi, D = 4 bis "
      "5,44 pi; für D = 1 ist es der Punkt 4 pi",
      ber["D = 1"][1] == 4.0 and all(ber[f"D = {D}"][1] > ber[f"D = {D - 1}"][1] + 0.3 for D in (2, 3, 4)))
# 3.4 Rohr mit einer Zelle Querschnitt und die modulare Uhr
eps, U, kr = streifen(2.0, 1.0, Ly=2)
sel = np.abs(eps) < GR
uhr = np.exp(0.5j * eps[sel])
par = np.array([np.sum(np.abs(U[:, i]) ** 2 * np.array([(-1) ** x[1] for x in kr])) for i in np.where(sel)[0]])
vz = np.sign(np.real(uhr[0] / 1j * par[0]))
check("3.4 Rohr mit einer Zelle Querschnitt (L_quer = 2, periodisch): Das Querhüpfen hebt sich weg, das Spektrum ist "
      "die native Leiter pi (4Z + 1) und pi (4Z + 3), und die modulare Uhr exp(i K/2) ist +-i (-1)^y: i (-1)^p auf "
      "der Querzelle, die Form der Kähler-Dirac-Z_4",
      np.allclose(np.abs(par), 1, atol=1e-8) and np.allclose(uhr, vz * 1j * par, atol=1e-4)
      and np.allclose(np.sort(np.abs(eps[sel]) / np.pi)[:4], [1, 1, 3, 3], atol=1e-5))
befund("3", "nein (gescopt)",
       "E/N: Im Halbraum eines hyperkubischen gestaffelten Gitters zerfällt das Feld je Querimpuls in Ketten mit "
       "M^2 = m^2 + 4 t_quer^2 Summe sin^2 k. Der Leiterabstand 4 pi K'/K wächst streng mit M. Für D > 1 und "
       "t_quer != 0 trägt der Halbraum daher ein Kontinuum von Abständen (D = 3: 4 pi bis 5,12 pi), nie die eine "
       "native Leiter. Sie erscheint nur für D = 1, für t_quer = 0 oder im Rohr von einer Zelle Querschnitt, also "
       "ohne Dynamik quer zum Schnitt. Dort ist die modulare Uhr i (-1)^p mit beiden Vorzeichen; die 16 der Naht "
       "trägt aber i auf allen Zuständen. Gilt für: ebener Schnitt längs einer Gitterebene, Nächste-Nachbar-Hüpfen, "
       "Grundzustand bei halber Füllung.")

# =====================================================================================================================
say("\nTeil 4  Die native 4D-Leiter Summe m^2 q^m = q (1 + q)/(1 - q)^3 als Gitter-Ecke?")
q = sp.symbols("q")
reihe = sp.Poly(sp.series(q * (1 + q) / (1 - q) ** 3, q, 0, 41).removeO(), q).all_coeffs()[::-1]
check("4.1 Zählidentität: q (1 + q)/(1 - q)^3 = Summe m^2 q^m (bis m = 40); Schale n von N^3 hat C(n+2,2) Punkte, und "
      "C(n+2,2) + C(n+1,2) = (n+1)^2",
      [int(c) for c in reihe] == [m * m for m in range(41)]
      and all(comb(n + 2, 2) + comb(n + 1, 2) == (n + 1) ** 2 for n in range(60)))


def monome(n, nv=4):
    return [e for e in itertools.product(range(n + 1), repeat=nv) if sum(e) == n]


ok = True
for n in range(2, 9):
    Pn = monome(n)
    pos = {e: i for i, e in enumerate(Pn)}
    rows = []
    for e in monome(n - 2):                                 # r^2 mal Monom vom Grad n - 2
        r = [0] * len(Pn)
        for i in range(4):
            f = list(e)
            f[i] += 2
            r[pos[tuple(f)]] += 1
        rows.append(r)
    kl = [e for e in Pn if e[3] <= 1]                       # Monome mit w-Exponent 0 oder 1
    for e in kl:
        r = [0] * len(Pn)
        r[pos[e]] = 1
        rows.append(r)
    ok &= len(kl) == (n + 1) ** 2 and len(rows) == len(Pn) and rank(np.array(rows)) == len(Pn)
check("4.2 Woher die Identität kommt: Die Monome x^a y^b z^c w^d mit d in {0, 1} sind eine Basis der harmonischen "
      "Polynome vom Grad n in vier Variablen modulo r^2 (Rang exakt, n = 2..8). Der 'Oktant' N^3 ist das Gitter der "
      "Exponenten (a, b, c), die 'Zweiplatzzelle' der Exponent d von w", ok)


def reihe_div(zaehler, nenner, N):
    """Potenzreihe zaehler/nenner bis q^N (ganzzahlig, nenner[0] = 1)."""
    out = []
    z = list(zaehler) + [0] * (N + 1)
    for n in range(N + 1):
        c = z[n] - sum(nenner[k] * out[n - k] for k in range(1, min(n, len(nenner) - 1) + 1))
        out.append(c)
    return out


NMAX = 12
gruppe = []
for perm in itertools.permutations(range(3)):
    for sg in itertools.product((1, -1), repeat=4):
        gruppe.append((perm, sg))
ok = True
for perm, sg in gruppe:
    g = sp.zeros(4, 4)
    for i in range(3):
        g[perm[i], i] = sg[i]
    g[3, 3] = sg[3]
    nen = [int(c) for c in sp.Poly((sp.eye(4) - q * g).det(), q).all_coeffs()[::-1]]
    harm = reihe_div([1, 0, -1], nen, NMAX)                 # (1 - q^2)/det(1 - q g): Charakter auf den Harmonischen
    for n in range(NMAX + 1):
        ch = 0
        for d in (0, 1):
            for e in monome(n - d, 3) if n - d >= 0 else []:
                if all(e[perm[i]] == e[i] for i in range(3)):
                    ch += int(np.prod([sg[i] ** e[i] for i in range(3)])) * sg[3] ** d
        ok &= ch == harm[n]
check("4.2 Die Identität ist äquivariant: Für alle 96 Elemente der Gruppe (vorzeichenbehaftete Permutationen von "
      "x, y, z; Vorzeichen von w) ist der Charakter auf den zwei Schalen gleich dem auf den Harmonischen der S^3 "
      "(n bis 12, exakt)", ok)


def s3_zerl(n):
    e, tau, c3 = (n + 1) ** 2, n + 1, 1 if n % 3 in (0, 1) else 0
    return {"trivial": (e + 3 * tau + 2 * c3) // 6, "Vorzeichen": (e - 3 * tau + 2 * c3) // 6, "zweidim.": (e - c3) // 3}


zerl = {f"m = {n + 1}": s3_zerl(n) for n in range(5)}
ZAHLEN["Niveau m (Dimension m^2) unter der Eckensymmetrie S_3: Vielfachheiten"] = zerl
say(f"      Niveau m unter der Eckensymmetrie S_3: {zerl}")


def haken(lam):
    lam = list(lam)
    n = sum(lam)
    konj = [sum(1 for x in lam if x > j) for j in range(lam[0])] if lam else []
    h = 1
    for i, li in enumerate(lam):
        for j in range(li):
            h *= li - j + konj[j] - i - 1
    return factorial(n) // h


def partitionen(n, mx=None):
    mx = n if mx is None else mx
    if n == 0:
        return [()]
    return [(k,) + p for k in range(min(n, mx), 0, -1) for p in partitionen(n - k, k)]


dimB4 = sorted(comb(4, k) * haken(l) * haken(mu) for k in range(5) for l in partitionen(k) for mu in partitionen(4 - k))
dimS4 = sorted(haken(l) for l in partitionen(4))
ZAHLEN["größte irreduzible Dimension: S_3, Würfelgruppe O_h, hyperoktaedrische Gruppe B_4"] = [2, max(dimS4), max(dimB4)]
check("4.3 Symmetrie: Niveau m ist unter SO(4) irreduzibel, (j,j) mit Dimension m^2. Die Eckengruppe S_3 hat höchstens "
      "Dimension 2, die Würfelgruppe 3, die hyperoktaedrische Gruppe des 4D-Gitters 8 (20 Darstellungen, Summe der "
      "Quadrate 384). Ab m = 2 (Ecke, Würfel) bzw. m = 3 (4D-Gitter) kann keine Gittersymmetrie die Entartung m^2 "
      "erzwingen; unter S_3 zerfällt schon m = 2 in 2 + 0 + 1",
      max(dimS4) == 3 and max(dimB4) == 8 and len(dimB4) == 20 and sum(d * d for d in dimB4) == 384
      and zerl["m = 2"] == {"trivial": 2, "Vorzeichen": 0, "zweidim.": 1}
      and all(v["trivial"] + v["Vorzeichen"] + 2 * v["zweidim."] == (i + 1) ** 2 for i, v in enumerate(zerl.values())))
ZIEL = np.array([1] + [2] * 4 + [3] * 9, dtype=float)


def passt(lev):
    lev = np.sort(np.abs(lev))[:14]
    return float(np.max(np.abs(lev / lev[0] - ZIEL))), [rd(x, 4) for x in lev / lev[0]]


# 4.4 Kandidat A: gestaffeltes Fermion in 3D am nativen Punkt m = 2t, Oktant einer offenen 12^3-Box
LB = 12
HA, koA = gitter((LB, LB, LB), 2.0)
okt = [i for i, x in enumerate(koA) if all(c < LB // 2 for c in x)]
epsA, _ = verschr(HA, okt)
abwA, verhA = passt(epsA)
nA = int(np.sum(np.abs(epsA) < 2 * np.min(np.abs(epsA))))
ZAHLEN["Kandidat A (gestaffeltes Fermion 12^3, Oktant 6^3)"] = {
    "kleinstes |eps| / pi": rd(np.min(np.abs(epsA)) / np.pi, 5), "erste 14 Verhältnisse": verhA,
    "Zahl der Niveaus unter dem Doppelten des kleinsten": nA, "größte Abweichung vom Ziel 1, 2 (4x), 3 (9x)": rd(abwA, 4)}
say(f"      Kandidat A: {ZAHLEN['Kandidat A (gestaffeltes Fermion 12^3, Oktant 6^3)']}")
check("4.4 Kandidat A (gestaffeltes Fermion, 3D, m = 2t, Oktant): passt nicht. Statt 1, 2 (4x), 3 (9x) beginnt bei "
      "0,75 pi ein dichtes Band von Ecken-, Kanten- und Flächenmoden (91 Niveaus unter dem Doppelten des kleinsten, "
      "das Ziel hätte dort eines)",
      abwA > 0.5 and nA > 20)
# 4.5 Kandidat B: harmonisches Gitter (massiver Skalar), Oktant


def boson_niveaus(Kmat, region):
    w, v = np.linalg.eigh(Kmat)
    Xh = (v * (0.5 / np.sqrt(w))) @ v.T
    Ph = (v * (0.5 * np.sqrt(w))) @ v.T
    XA, PA = Xh[np.ix_(region, region)], Ph[np.ix_(region, region)]
    wx, vx = np.linalg.eigh(XA)
    Xs = (vx * np.sqrt(wx)) @ vx.T
    nu2 = np.linalg.eigvalsh(Xs @ PA @ Xs)
    nuu = np.sqrt(np.clip(nu2, 0.25 + 1e-16, None))
    return np.sort(np.log((nuu + 0.5) / (nuu - 0.5))), np.sort(nu2)


def kette_K(N, w2):
    Kc = (w2 + 2.0) * np.eye(N)
    for n in range(N - 1):
        Kc[n, n + 1] = Kc[n + 1, n] = -1.0
    return Kc


def summe3(K1):
    N = K1.shape[0]
    I = np.eye(N)
    return np.kron(np.kron(K1, I), I) + np.kron(np.kron(I, K1), I) + np.kron(np.kron(I, I), K1)


LBo = 10
K1 = kette_K(LBo, 1.0 / 3)                                  # je Richtung omega0^2/3: zusammen (omega0^2 + 6) - Nachbarn
koB = list(itertools.product(range(LBo), repeat=3))
oktB = [i for i, x in enumerate(koB) if all(c < LBo // 2 for c in x)]
epsB, _ = boson_niveaus(summe3(K1), oktB)
abwB, verhB = passt(epsB)
ZAHLEN["Kandidat B (harmonisches Gitter 10^3, omega0 = 1, Oktant 5^3)"] = {
    "kleinstes eps": rd(epsB[0], 5), "erste 14 Verhältnisse": verhB,
    "größte Abweichung vom Ziel 1, 2 (4x), 3 (9x)": rd(abwB, 4)}
say(f"      Kandidat B: {ZAHLEN['Kandidat B (harmonisches Gitter 10^3, omega0 = 1, Oktant 5^3)']}")
check("4.5 Kandidat B (massiver Skalar auf dem kubischen Gitter, Nächste-Nachbar-Kopplung, Oktant): passt nicht "
      "(wieder ein Band von Flächenmoden, keine Entartungen 1, 4, 9)", abwB > 0.5)
# 4.6 Kandidat C: Produktstrukturen. Bosonen mit Produktkopplung K = K1 x K1 x K1; Fermionen mit Vorzeichenprodukt.
LP = 6
K1p = kette_K(LP, 1.0)
koP = list(itertools.product(range(LP), repeat=3))
oktP = [i for i, x in enumerate(koP) if all(c < LP // 2 for c in x)]
_, nu2_dir = boson_niveaus(np.kron(np.kron(K1p, K1p), K1p), oktP)
_, nu2_1 = boson_niveaus(K1p, list(range(LP // 2)))
nu2_prod = np.sort([16 * a * b * c for a in nu2_1 for b in nu2_1 for c in nu2_1])
e1 = np.log((np.sqrt(nu2_1) + 0.5) / (np.sqrt(nu2_1) - 0.5))
eprod = 2 * np.arctanh(1 / np.sqrt(4 * np.min(nu2_dir)))
ZAHLEN["Kandidat C (Produktkopplung 6^3)"] = {"Niveaus der Halbkette": [rd(x, 5) for x in np.sort(e1)],
                                             "kleinstes Niveau des Oktanten": rd(eprod, 5),
                                             "additiv wäre": rd(3 * np.min(e1), 5)}
say(f"      Kandidat C (Bosonen, Produktkopplung): {ZAHLEN['Kandidat C (Produktkopplung 6^3)']}")
check("4.6 Kandidat C, Bosonen mit Produktkopplung K1 x K1 x K1: coth(eps/2) = Produkt der coth(eps_i/2) (27 Werte, "
      "1e-9 relativ). Das Gesetz ist multiplikativ, nicht additiv",
      np.allclose(nu2_dir, nu2_prod, rtol=1e-9) and abs(eprod - 3 * np.min(e1)) > 0.5)
nat = [mp.pi * (2 * l + 1) for l in range(3)]
ferm = sorted(float(2 * mp.atanh(mp.tanh(a / 2) * mp.tanh(b / 2) * mp.tanh(c / 2)) / mp.pi)
              for a in nat for b in nat for c in nat)
ZAHLEN["Kandidat C (Fermionen, Vorzeichenprodukt der nativen Kette): |eps| / pi"] = [rd(x, 5) for x in ferm[:10]]
say(f"      Kandidat C (Fermionen, sgn(h) = sgn(h1) x sgn(h2) x sgn(h3), native Halbketten): |eps| / pi = "
    f"{[rd(x, 5) for x in ferm[:10]]}")
check("4.6 Kandidat C, Fermionen mit Vorzeichenprodukt: tanh(eps/2) = Produkt der tanh(eps_i/2). Aus der nativen Leiter "
      "entsteht 0,652 pi, dann ein Häufungspunkt bei 0,780 pi (alle Tripel (1, 1, l)), dann 0,999 pi: keine Leiter "
      "a + b + c (additiv wäre das unterste Niveau 3 pi)",
      abs(ferm[0] - 0.652) < 1e-3 and abs(ferm[1] - 0.780) < 1e-3 and abs(ferm[6] - 0.780) < 1e-3
      and abs(ferm[7] - 0.999) < 1e-3)
tw = Counter()
for na in itertools.product(range(9), repeat=4):            # Twistoroszillatoren a_1, a_2, b_1, b_2 mit N_a = N_b
    if na[0] + na[1] == na[2] + na[3] <= 8:
        tw[1 + na[0] + na[1]] += 1                          # E = 1 + (N_a + N_b)/2
osz = Counter()
for e in itertools.product(range(9), repeat=3):             # drei Oszillatoren und ein Zweiniveausystem
    for d in (0, 1):
        if sum(e) + d <= 8:
            osz[1 + sum(e) + d] += 1
zw = tw == osz == Counter({m: m * m for m in range(1, 10)})
check("4.7 Was die Zählung wirklich ist: drei Oszillatoren plus ein Zweiniveausystem (Exponenten von x, y, z und w), "
      "gleichwertig die vier Twistoroszillatoren des Bestands mit N_a = N_b (Zustände (n+1)^2, z73/z79). Das ist ein "
      "Fockraum, kein Ortsgitter", zw)
befund("4", "nein (gescopt)",
       "E/N: q (1 + q)/(1 - q)^3 = Summe m^2 q^m ist die Monombasis der Harmonischen auf S^3: Der 'Oktant' ist das "
       "Exponentengitter dreier Oszillatoren, die 'Zelle' der Exponent der vierten Variablen. Die Identität ist unter "
       "den 96 vorzeichenbehafteten Permutationen äquivariant, mehr Symmetrie hat keine Gitterecke; SO(4) verlangt "
       "aber irreduzible Niveaus der Dimension m^2. Drei Kandidaten für ein Ortsgitter sind gerechnet und passen "
       "nicht: gestaffeltes Fermion am nativen Punkt (Oktant, Band von Rand- und Flächenmoden), massiver Skalar "
       "(ebenso), Produktstrukturen (multiplikatives statt additives Gesetz). Kein Unmöglichkeitsbeweis für alle "
       "Gittermodelle. Die 4D-Leiter ist die Leiter von Oszillatoren (Twistoren), nicht die Verschränkungsleiter "
       "einer Ecke.")

# =====================================================================================================================
say("\nTeil 5  Bilanz")
bilanz = [
    ("E", "Zelle", "16 von Spin(10) = Lambda*(C^4) unter SU(4) in SU(5) und = 1 + (2,2) + (3,1) + (1,3) + (2,2) + 1 "
     "unter SO(4): der Formeninhalt einer Kähler-Dirac-Zelle in d = 4. Dasselbe über die Diagonale in Spin(4) x Spin(6)."),
    ("E negativ", "Zelle gegen Bestand", "Als innere Zelle scheidet das aus: Die verdrillte Drehgruppe läge in Spin(10), "
     "mischte Farb- und schwache Plätze und vertauschte mit höchstens einem der zwölf Eichgeneratoren. Eine Zelle sind "
     "außerdem nur 8 Weyl-Fermionen. 2^4 = 16 und dim 16 = 2^(5-1) sind verschiedene Sechzehnen."),
    ("E", "Lesart B", "Mit Lorentz außen: Der Indexkern K = (2_L x 16) + (2_R x 16-quer) aus z150 ist unter "
     "diag(Lorentz x SU(2)_L x SU(2)_R) x SU(4) genau Lambda^gerade x 4 + Lambda^ungerade x 4-quer, vier reduzierte "
     "Kähler-Dirac-Felder. T16 ist das nicht."),
    ("E", "Z_4", "Kähler-Dirac-Z_4 = i (-1)^p ist in Spin(10) das Element z^-1 mal Zentrum von SU(2)_L (Lesart B) "
     "bzw. z mal SM-Toruselement (Zelle a): nicht die Uhr, aber die Uhr modulo Eichtransformation. In Lesart B liegt "
     "von der Kähler-Dirac-U(1) genau diese Z_4 in Spin(10)."),
    ("E", "Anomalie", "nu(Uhr) = 16, nu(Kähler-Dirac) = 8 - 8 = 0 auf einer Generation; für k reduzierte Felder 4k "
     "gegen Instantonladung k mod 4. Beide verlangen k = 0 mod 4. Das deckt sich mit Summe q = 0 mod 4 aus [L]."),
    ("C", "Dimension vier", "Unter P_KD (Lorentz ist mit einem Spin(d)-Faktor von Spin(10) verdrillt, Träger ist ein "
     "Weyl-Feld mal 16) sind nur d = 4 und d = 8 möglich; die native Teilung 3 + 2 lässt d = 4. Die Vier ist dann die "
     "reelle Dimension der zwei schwachen Plätze. Das ist die kinematische Lesart von K12 in neuer Form, keine "
     "Herleitung: P_KD ist nicht nativ (z150: so(3,2) vertauscht mit allen Eichgeneratoren)."),
    ("E/N negativ", "Halbraum", "Ein hyperkubisches gestaffeltes Gitter in D > 1 hat als Halbraum ein Kontinuum von "
     "Leiterabständen (D = 3: 4 pi bis 5,12 pi). Die native Leiter gibt es nur ohne Querhüpfen. P_CTM lässt sich so "
     "nicht auf vier Dimensionen heben."),
    ("E/N negativ", "4D-Leiter", "Summe m^2 q^m ist die Monombasis der S^3-Harmonischen (Exponentengitter, "
     "Oszillatoren). Drei Ortsgitter-Kandidaten für eine Ecke passen nicht; Gittersymmetrie kann m^2 nicht erzwingen."),
    ("L", "Literatur", "Anomalie ~ Euler-Zahl, Rest-Z_4, zwei Kähler-Dirac-Felder = 16 Majorana, vier reduzierte Felder "
     "= Pati-Salam-Inhalt, Spiegelfermionen und ihr Entkoppeln: arXiv:2101.01026, 2209.03828, 2311.02487. Nur über "
     "Kurzfassung und Auszug gelesen; die Gitterdynamik dort ist Vermutung der Autoren."),
    ("C", "Neue Prämisse", "P_KD: Der Lorentz-Spin der Naht ist mit der schwachen Spin(4) verdrillt, und K lebt als "
     "reduziertes Kähler-Dirac-Feld auf einem vierdimensionalen Gitter. Folge: SU(2)_L wäre die Geschmackssymmetrie "
     "des Gitters, exakt erst im Kontinuum [L]; die Hyperladung enthielte T_3R aus der verdrillten Gruppe."),
    ("O", "Offen", "Ob ein natives Datum die Verdrillung erzwingt; ob es ein 4D-Gitter mit der nativen Leiter gibt "
     "(nicht hyperkubisch-gestaffelt mit ebenem Schnitt); Minkowski-Fassung der Verdrillung (Spin(4) ist kompakt, "
     "nativ ist so(3,2))."),
    ("-", "K12, K16, K18", "K12 bleibt L (kinematisch, jetzt als P_KD benannt). K16 unberührt: kein Kontinuumslimes, "
     "der Halbraum zeigt eher, dass die 1D-Kette nicht Teil eines isotropen 4D-Gitters ist. K18 unberührt: Die "
     "Kopiengruppe SU(4) ist die Farb-Spin(6) von Pati-Salam, nicht die Familien-A_3; A_3 als Geschmacksgruppe gäbe "
     "vier Geschmäcker und halbe Zellen."),
    ("-", "Neun Prämissen", "A*, I_step, P_lift, P_IF, P_univ, P_uhr, P_ord, P_curv, P_G unverändert. P_KD käme "
     "hinzu und ersetzt keine. P_CTM (z152) bleibt auf 1+1 Dimensionen beschränkt."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
urteil = ("Nein in der wörtlichen Form, bedingt in einer anderen. (1) Die 16 der Naht hat exakt den Formeninhalt einer "
          "vierdimensionalen Kähler-Dirac-Zelle, 1 + (2,2) + (3,1) + (1,3) + (2,2) + 1, auf zwei Wegen. Als Zelle "
          "gelesen müsste aber die verdrillte Drehgruppe in Spin(10) liegen; sie mischt dann Farb- und schwache Plätze "
          "und lässt vom Standardmodell höchstens einen Erzeuger. Eine Zelle sind zudem nur 8 Weyl-Fermionen. Die 16 "
          "der Literatur sind 16 Weyl-Fermionen aus zwei Zellen oder vier reduzierten Feldern, nicht 2^4. Die Aussage "
          "im Ausblick von z152 ist in diesem Punkt zu berichtigen. (2) Verträglich mit dem Bestand ist die Lesart "
          "mit Lorentz außen: Der Indexkern K aus z150 ist unter der Diagonale von Lorentz und schwacher Spin(4), mal "
          "der Farb-SU(4), exakt vier reduzierte Kähler-Dirac-Felder. Die Kähler-Dirac-Z_4 ist dann die inverse Uhr "
          "mal dem Zentrum von SU(2)_L, und beide Anomaliebedingungen sind dieselbe: vier reduzierte Felder, 16 "
          "Weyl-Fermionen. Das rechnet die Pati-Salam-Aussage der Literatur für den Bestand nach. (3) Die Dimension "
          "vier folgt daraus nicht nativ. Unter der neuen Prämisse P_KD lässt die Darstellungstheorie d = 4 oder 8 "
          "zu, und die native Platzteilung 3 + 2 wählt d = 4; die Vier ist dann die reelle Dimension der zwei "
          "schwachen Plätze. Die Verdrillung selbst ist nicht nativ: Im Bestand vertauscht die Raumzeitalgebra mit "
          "allen Eichgeneratoren. K12 bleibt Lesart. (4) Die Gitterlesart der Naht (P_CTM) lässt sich nicht "
          "hyperkubisch auf vier Dimensionen heben: Der Halbraum hat für D > 1 ein Kontinuum von Leiterabständen, "
          "die native Leiter gibt es nur ohne Querhüpfen. (5) Die native 4D-Leiter Summe m^2 q^m ist die Monombasis "
          "der S^3-Harmonischen, ein Exponentengitter von Oszillatoren; als Verschränkungsleiter einer Gitterecke "
          "passt keiner der drei gerechneten Kandidaten, und keine Gittersymmetrie kann die Entartung m^2 erzwingen. "
          "(6) Für den Kandidaten ändert sich nichts: neun Prämissen, P_KD käme als zehnte hinzu und ersetzt keine. "
          "K13, K15, K16 und K18 sind nicht berührt. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z154_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z154.txt").write_text("\n".join(OUT) + "\n")
