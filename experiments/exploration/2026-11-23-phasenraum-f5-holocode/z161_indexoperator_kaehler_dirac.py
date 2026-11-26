"""Z161: Hat der native Indexoperator auf dem Indexkern K die Kähler-Dirac-Form d - delta? Test der Prämisse P_KD aus z154.

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart, N numerisch, O offen.

Frage (Welle C nach H6h50, Folgetest aus z154/B1): z154 fand als einzige mit dem Bestand verträgliche Kähler-Dirac-Lesart
'Lesart B': Lorentz außen, Geschmacksgruppe = schwache Spin(4) = SU(2)_L x SU(2)_R, Kopiengruppe = Farb-SU(4); unter der
Diagonale ist der Indexkern K = S^+ cap P_- (z118/z150) genau Lambda^gerade x 4 + Lambda^ungerade x 4-quer. Das braucht die
Prämisse P_KD (Lorentz-Spin mit schwachem Spin verdrillt). Test: Hat der native Indexoperator die Form d - delta, also
(i) antikommutiert er mit Gamma_KD = gamma_5 x schwache Chiralität, (ii) ändert er den Formengrad um +-1, (iii) ist er
unter der verdrillten Gruppe ein Vektor? Streng: Das zählt nur, wenn es ohne Wahl der Gruppe folgt.

Bestand: z118 (256 = Cl(4,2) x Cl(10), S^+, Nullpaar gamma_+-, P_+-, Gamma_2, R = -i Gamma_2 gamma_5 auf S^+,
D_perp = gamma_- A + gamma_+ A^+ mit A_q = S^+(1 - q S^4), J), z150 (Träger T16 und K, so(3,2), Eichgeneratoren, direktes
Produkt), z154 (Lesart B, Kähler-Dirac-Z_4 = Uhr^-1 mal Zentrum von SU(2)_L), z129 (Kern von D_perp), z127 (Gewichte).

Rechenweise: Alle Matrizen haben Einträge in Z[i, 1/2] (dyadische Gauß-Zahlen); numpy rechnet darauf ohne Rundung, alle
Gleichungen werden mit exakter Gleichheit geprüft, Ränge exakt über flint (ganzzahlige reelle Einbettung), Vielfachheiten
über exakte Spuren mit Brüchen. Numerisch (N) ist nur der Markenteil (A_q, Wilson-Ecken).

Teil 1  Aufbau: Clifford-Daten, Träger, Eichgeneratoren, die verdrillte Gruppe und der Formengrad N_form (Lötung phi).
Teil 2  Der Test: (a) alle Lorentz-Vektor-Eich-Singulett-Operatoren je Träger, mit Uhrladung; (b) D_perp und der
        induzierte Weyl-Operator gegen Gamma_KD, N_form und die verdrillte Gruppe, Tautologieprobe; (c) welches native
        Datum die verdrillte von der unverdrillten Gruppe unterscheidet (reelle Form, J, so(3,2), Uhr, Eichgruppe,
        Gewichte der Mischspur).
Teil 3  Gitterseite: native Vierteldrehung, Wilson-Fortsetzung, Eckengrad.
Teil 4  Folgen für die Z_4-Anomalie, K12, K15, K16.
Teil 5  Bilanz und Urteil.
"""
import itertools
import json
import os
import time
from collections import Counter
from fractions import Fraction as Fr
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # deterministische Rundung (einfädig)

import numpy as np  # noqa: E402
from flint import fmpz_mat  # noqa: E402

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
    return bool(cond)


def befund(name, val, text):
    BEFUNDE[name] = val
    say(f"   [Befund: {val}] {text}")


# ---------------------------------------------------------------------------------------------------- Hilfsfunktionen
s0 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]], complex)


def kr(*m):
    out = m[0]
    for x in m[1:]:
        out = np.kron(out, x)
    return out


def euklid(n):
    G = []
    for j in range(n):
        G.append(kr(*([sz]*j + [sx] + [s0]*(n - j - 1))))
        G.append(kr(*([sz]*j + [sy] + [s0]*(n - j - 1))))
    return G


def eq(A, B):
    return bool(np.array_equal(A, B))


def null(A):
    return not np.any(A)


def com(A, B):
    return A @ B - B @ A


def acom(A, B):
    return A @ B + B @ A


def prod(L, d):
    M = np.eye(d, dtype=complex)
    for g in L:
        M = M @ g
    return M


def chiral(G):
    M = prod(G, G[0].shape[0])
    for ph in (1, 1j, -1, -1j):
        if eq((ph*M) @ (ph*M), np.eye(M.shape[0])):
            return ph*M
    raise ValueError


def frac(z):
    """Exakter Bruch aus einer dyadischen Zahl (Imaginärteil muss null sein)."""
    z = complex(z)
    assert z.imag == 0, z
    return Fr(z.real)


def ganz(W):
    sc = 1
    while not np.array_equal(W*sc, np.round(W*sc)):
        sc *= 2
        assert sc < 2**20
    Wi = np.round(W*sc).astype(np.int64)
    return fmpz_mat(Wi.shape[0], Wi.shape[1], [int(x) for x in Wi.ravel()])


def crank(mats):
    """Komplexe Dimension des Spanns einer Liste komplexer Matrizen (exakt)."""
    mats = [m for m in mats if np.any(m)]
    if not mats:
        return 0
    V = np.array([m.ravel() for m in mats])
    V = V[:, np.any(V != 0, axis=0)]
    V = np.unique(V, axis=1)
    r = ganz(np.block([[V.real, V.imag], [-V.imag, V.real]])).rank()
    assert r % 2 == 0
    return r // 2


def rrank(mats):
    """Reelle Dimension des reellen Spanns (exakt)."""
    V = np.array([np.concatenate([m.ravel().real, m.ravel().imag]) for m in mats])
    V = V[:, np.any(V != 0, axis=0)]
    return ganz(np.unique(V, axis=1)).rank()


def crank_mat(M):
    """Komplexer Rang einer Koeffizientenmatrix (exakt)."""
    M = M[np.any(M != 0, axis=1)]
    M = np.unique(M, axis=0)
    return ganz(np.block([[M.real, -M.imag], [M.imag, M.real]])).rank() // 2


def lagrange(vals, v):
    """Koeffizienten (aufsteigend) des Lagrange-Polynoms, das auf vals bei v eins und sonst null ist."""
    co = [Fr(1)]
    for w in vals:
        if w == v:
            continue
        neu = [Fr(0)]*(len(co) + 1)
        for k, c in enumerate(co):
            neu[k] += c*Fr(-w, v - w)
            neu[k + 1] += c*Fr(1, v - w)
        co = neu
    return co


def vielfachheiten(P, X, vals):
    """Vielfachheiten der Eigenwerte vals von X auf dem Bild des Projektors P (exakt über Potenzspuren)."""
    tr, M = [], P.copy()
    for _ in vals:
        tr.append(frac(np.trace(M)))
        M = M @ X
    res = []
    for v in vals:
        co = lagrange(vals, v)
        d = sum(c*t for c, t in zip(co, tr))
        assert d.denominator == 1
        res.append(int(d))
    return res


def is_sl2(t):
    e, f, h = t
    return eq(com(e, f), h) and eq(com(h, e), 2*e) and eq(com(h, f), -2*f)


def cas(t):
    """4 j (j + 1) für ein sl(2)-Tripel."""
    e, f, h = t
    return h @ h + 2*(e @ f + f @ e)


def inhalt(P, CL, CR, d=256):
    """Inhalt des Bilds von P unter sl(2) x sl(2) mit Casimiren CL, CR: {(dim_L, dim_R): Vielfachheit}."""
    powL, powR = [np.eye(d), CL, CL @ CL], [np.eye(d), CR, CR @ CR]
    T = {(i, j): frac(np.trace(P @ powL[i] @ powR[j])) for i in range(3) for j in range(3)}
    dims, res = {0: 1, 3: 2, 8: 3}, {}
    for a in (0, 3, 8):
        ca = lagrange((0, 3, 8), a)
        for b in (0, 3, 8):
            cb = lagrange((0, 3, 8), b)
            dd = sum(ca[i]*cb[j]*T[i, j] for i in range(3) for j in range(3))
            if dd != 0:
                m = dd/(dims[a]*dims[b])
                assert m.denominator == 1
                res[f"({dims[a]},{dims[b]})"] = int(m)
    return res


def triple(S, P):
    return ((S[0] + 1j*S[1]) @ P, (S[0] - 1j*S[1]) @ P, 2*S[2] @ P)


def gruppenordnung(gens, grenze=200):
    start = gens[0] @ np.linalg.inv(gens[0]) + 0.0
    start = np.round(start.real) + 1j*np.round(start.imag) + 0.0
    els = {start.tobytes()}
    front = [start]
    while front:
        neu = []
        for M in front:
            for g in gens:
                Pm_ = M @ g + 0.0
                k = Pm_.tobytes()
                if k not in els:
                    els.add(k)
                    neu.append(Pm_)
        front = neu
        if len(els) > grenze:
            return None
    return len(els)


# ============================================================================================================ Teil 1
say("Z161: Hat der native Indexoperator auf dem Indexkern K die Kähler-Dirac-Form d - delta? (Test von P_KD aus z154)")
say("=" * 100)
say("\nTeil 1  Aufbau (exakt: Einträge in Z[i, 1/2], Gleichheit ohne Toleranz, Ränge über flint)")

Gi = euklid(5)
G11 = chiral(Gi)
I8, I32, I256 = np.eye(8, dtype=complex), np.eye(32, dtype=complex), np.eye(256, dtype=complex)
Z256 = np.zeros((256, 256), complex)
E6 = euklid(3)
hh = [E6[0], E6[1], E6[2], E6[3], 1j*E6[4], 1j*E6[5]]                # Quadrate (+, +, +, +, -, -): R^(4,2)
h7 = chiral(hh)
G16 = [kr(x, I32) for x in hh] + [1j*kr(h7, gi) for gi in Gi]         # 16 Erzeuger, Signatur (4, 12)
Gt = chiral(G16)
Sp, Sm = (I256 + Gt)/2, (I256 - Gt)/2
gp, gm = (G16[3] + G16[5])/2, (G16[3] - G16[5])/2                     # Nullpaar der Dilatation (z118)
Pplus, Pminus = gm @ gp, gp @ gm
Gam2 = Pplus - Pminus
LI = (0, 1, 2, 4)                                                    # Lorentzrichtungen quer zum Nullpaar
ETA = [1, 1, 1, -1]
ga = [G16[a] for a in LI]
g5 = chiral(ga)
R = kr(I8, 1j*G11)                                                   # Uhr = i Gamma_11 (Zentrum von Spin(10))
Rinv = -R
H7 = kr(h7, I32)
PK, PKs = Sp @ Pminus, Sp @ Pplus                                    # Indexkern K und sein Komplement K' in S^+
PKm, PKsm = Sm @ Pminus, Sm @ Pplus                                  # dieselben Sektoren in S^-
PT16 = kr(np.diag((np.diag(h7).real > 0).astype(float)), np.diag((np.diag(G11).real > 0).astype(float))).astype(complex)
quad = [int((g @ g)[0, 0].real) for g in G16]
clifford = all(eq(acom(G16[i], G16[j]), 2*(i == j)*quad[i]*I256) for i in range(16) for j in range(16))
rang = {"S+": int(np.trace(Sp).real), "K": int(np.trace(PK).real), "K'": int(np.trace(PKs).real),
        "T16": int(np.trace(PT16).real), "T16 cap K": int(np.trace(PT16 @ PK).real)}
ZAHLEN["1.1 Clifford"] = {"Quadrate der 16 Erzeuger": quad, "Ränge": rang}
check(f"1.1 [E] 256 = Cl(4,2) x Cl(10): 16 antikommutierende Erzeuger, Signatur (4, 12) (Quadrate {quad[:6]} und zehnmal "
      f"-1); Gamma_ges = h7 x Gamma_11; Ränge {rang}",
      clifford and quad == [1, 1, 1, 1, -1, -1] + [-1]*10 and eq(Gt, kr(h7, G11))
      and rang == {"S+": 128, "K": 64, "K'": 64, "T16": 64, "T16 cap K": 32} and eq(R @ R, -I256))
check("1.2 [E] Bestand z118/z150: auf S^+ ist R = -i Gamma_2 gamma_5 = i (h7 x 1); auf K ist Gamma_2 = -1, R = i gamma_5 "
      "und gamma_5 = Gamma_11; auf T16 ist R = i und Gamma_2 = -gamma_5",
      eq(R @ Sp, -1j*Gam2 @ g5 @ Sp) and eq(R @ Sp, 1j*H7 @ Sp) and eq(Gam2 @ PK, -PK) and eq(R @ PK, 1j*g5 @ PK)
      and eq(g5 @ PK, kr(I8, G11) @ PK) and eq(R @ PT16, 1j*PT16) and eq(Gam2 @ PT16, -g5 @ PT16))
EPS4 = np.zeros((4, 4, 4, 4))
for p_ in itertools.permutations(range(4)):
    EPS4[p_] = round(np.linalg.det(np.eye(4)[list(p_)]))
spur_ok = True
for a, b, c_, d in itertools.product(range(4), repeat=4):
    if a < b and c_ < d:
        X = ga[a] @ ga[b] @ ga[c_] @ ga[d]
        spur_ok &= np.trace(R @ PK @ X) == 64*EPS4[a, b, c_, d]
        spur_ok &= np.trace(Gam2 @ PK @ X) == -64*(ETA[a]*(a == d)*ETA[b]*(b == c_) - ETA[a]*(a == c_)*ETA[b]*(b == d))
check("1.3 [E] Kontrolle z118: Tr_K(R g_a g_b g_c g_d) = 64 eps_abcd und Tr_K(Gamma_2 ...) = -64 (eta_ad eta_bc - "
      "eta_ac eta_bd), alle 36 Paare exakt", spur_ok)

# so(3,2) (z150): Richtungen (0, 1, 2, 4, 5)
IDX5 = (0, 1, 2, 4, 5)
ETA5 = [quad[a] for a in IDX5]
S5 = [[G16[IDX5[i]] @ G16[IDX5[j]]/2 if i != j else Z256 for j in range(5)] for i in range(5)]
PAIRS5 = list(itertools.combinations(range(5), 2))
PAIRS4 = list(itertools.combinations(range(4), 2))
SIG32 = [S5[i][j] for i, j in PAIRS5]
LOR = [S5[i][j] for i, j in PAIRS4]                                  # = gamma_a gamma_b / 2
PI = [S5[a][4] for a in range(4)]                                    # Translationen Pi_a = Sigma_a5
so32_ok = all(eq(com(S5[i][j], S5[k][l]), ETA5[j]*(j == k)*S5[i][l] - ETA5[i]*(i == k)*S5[j][l]
                 - ETA5[j]*(j == l)*S5[i][k] + ETA5[i]*(i == l)*S5[j][k])
              for i, j in PAIRS5 for k, l in PAIRS5)
check("1.4 [E] so(3,2) aus den Richtungen (0, 1, 2, 4, 5): Relationen exakt, [Pi_a, Pi_b] = Sigma_ab; alle zehn "
      "Erzeuger vertauschen mit der Uhr; mit Gamma_2 vertauschen genau die sechs Lorentz-Erzeuger (auf K wirkt nur "
      "Lorentz, wie z150)",
      so32_ok and all(eq(com(PI[a], PI[b]), LOR[PAIRS4.index((a, b))]) for a, b in PAIRS4)
      and all(null(com(X, R)) for X in SIG32) and all(null(com(X, Gam2)) for X in LOR)
      and all(np.any(com(X, Gam2)) for X in PI))

# Eichgeneratoren auf dem Fockraum der fünf Plätze (0, 1, 2 Farbe; 3, 4 schwach)
A5 = [(Gi[2*j] + 1j*Gi[2*j + 1])/2 for j in range(5)]
idx16 = np.where(np.diag(G11).real > 0)[0]
YQ = sorted([Fr(0)] + [Fr(1)] + [Fr(-2, 3)]*3 + [Fr(1, 6)]*6 + [Fr(1, 3)]*3 + [Fr(-1, 2)]*2)
OPS = None
for kand in (A5, [x.conj().T for x in A5]):
    nn = [k_.conj().T @ k_ for k_ in kand]
    y6 = -2*(nn[0] + nn[1] + nn[2]) + 3*(nn[3] + nn[4])              # 6 Y
    if sorted(Fr(int(round(v)), 6) for v in np.diag(y6).real[idx16]) == YQ:
        OPS = kand
        break
an = OPS
ad = [x.conj().T for x in an]
NN = [ad[j] @ an[j] for j in range(5)]
Y6 = -2*(NN[0] + NN[1] + NN[2]) + 3*(NN[3] + NN[4])


def E_(i, j):
    return ad[i] @ an[j]


SL_L = (E_(3, 4), E_(4, 3), NN[3] - NN[4])                           # SU(2)_L
SL_R = (ad[3] @ ad[4], an[4] @ an[3], NN[3] + NN[4] - I32)           # SU(2)_R
SU3 = [E_(i, j) for i in range(3) for j in range(3) if i != j] + [NN[0] - NN[1], NN[1] - NN[2]]
SM12 = SU3 + list(SL_L) + [Y6]
SO10 = [Gi[i] @ Gi[j]/2 for i, j in itertools.combinations(range(10), 2)]
SO6 = [Gi[i] @ Gi[j]/2 for i, j in itertools.combinations(range(6), 2)]
SO4W = [Gi[i] @ Gi[j]/2 for i, j in itertools.combinations(range(6, 10), 2)]
chiw = -(I32 - 2*NN[3]) @ (I32 - 2*NN[4])                            # schwache Chiralität: +1 auf SU(2)_L-Dubletts
chic = G11 @ chiw                                                    # Farbchiralität: +1 auf der 4 von SU(4)
CW, CC = kr(I8, chiw), kr(I8, chic)
cL = kr(I8, I32 - 2*(SL_L[2] @ SL_L[2]))                             # Zentrum von SU(2)_L: exp(i pi h_L)
car = all(eq(acom(an[i], ad[j]), (i == j)*I32) and null(acom(an[i], an[j])) for i in range(5) for j in range(5))


def kommutant_dim(gens, cart):
    """Dimension des Kommutanten auf der 32 (exakt); cart = diagonale Erzeuger (Gewichte), gens = übrige."""
    w = np.array([np.diag(c_).real for c_ in cart]).T
    pos = [(s_, t_) for s_ in range(32) for t_ in range(32) if np.array_equal(w[s_], w[t_])]
    rows = []
    for T in gens:
        M = np.zeros((1024, len(pos)), complex)
        for k, (s_, t_) in enumerate(pos):
            M[np.arange(32)*32 + t_, k] += T[:, s_]
            M[s_*32 + np.arange(32), k] -= T[t_, :]
        rows.append(M)
    return len(pos) - crank_mat(np.vstack(rows))


SU3g, SU3c = [E_(0, 1), E_(1, 0), E_(1, 2), E_(2, 1)], [NN[0] - NN[1], NN[1] - NN[2]]
SO6g = SU3g + [ad[0] @ ad[1], an[1] @ an[0]]
kd = {"so(10)": kommutant_dim([E_(i, i + 1) for i in range(4)] + [E_(i + 1, i) for i in range(4)]
                              + [ad[3] @ ad[4], an[4] @ an[3]], NN),
      "Pati-Salam": kommutant_dim(SO6g + [SL_L[0], SL_L[1], SL_R[0], SL_R[1]], NN),
      "Standardmodell": kommutant_dim(SU3g + [SL_L[0], SL_L[1]], SU3c + [SL_L[2], Y6]),
      "Farb-SU(4) allein": kommutant_dim(SO6g, NN[:3])}
Nc, Nw = NN[0] + NN[1] + NN[2], NN[3] + NN[4]


def dproj(nc, nw):
    d = ((np.round(np.diag(Nc).real) == nc) & (np.round(np.diag(Nw).real) == nw)).astype(float)
    return np.diag(d).astype(complex)


leer = int(np.where(np.diag(dproj(0, 0)).real > 0)[0][0])
voll = int(np.where(np.diag(dproj(3, 2)).real > 0)[0][0])
Eu = np.zeros((32, 32), complex)
Eu[voll, leer] = 1
C_G = {"so(10)": [I32, G11], "Pati-Salam": [I32, G11, chiw, G11 @ chiw],
       "Standardmodell": [dproj(nc, nw) for nc in range(4) for nw in range(3)] + [Eu, Eu.T.copy()]}
GENS_G = {"so(10)": SO10, "Pati-Salam": SO6 + list(SL_L) + list(SL_R), "Standardmodell": SM12}
basis_ok = all(all(null(com(c_, g)) for c_ in C_G[k] for g in GENS_G[k]) and crank(C_G[k]) == kd[k] for k in C_G)
ZAHLEN["1.5 Kommutanten auf der 32"] = kd
check(f"1.5 [E] Eichgeneratoren auf dem Fockraum (CAR exakt, Hyperladungsspektrum der 16 wie z82/z150, sl(2)-Tripel "
      f"SU(2)_L und SU(2)_R); Kommutantendimensionen auf der 32: {kd}; explizite Basen (so(10): 1, Gamma_11; "
      f"Pati-Salam: dazu die schwache Chiralität; Standardmodell: 12 Projektoren und zwei Abbildungen zwischen nu^c und "
      f"seinem Konjugierten) vertauschen und haben diesen Rang",
      OPS is not None and car and is_sl2(SL_L) and is_sl2(SL_R) and basis_ok
      and kd == {"so(10)": 2, "Pati-Salam": 4, "Standardmodell": 14, "Farb-SU(4) allein": 32}
      and all(null(com(kr(I8, T), X)) for T in SM12 for X in SIG32))

# Verdrillte Gruppe und Formengrad. Lötung phi: gamma_a -> c_a = i kappa_a s_a g^w_a (Clifford-verträglich, c_a^2 = eta_a)
GW = [G16[12 + i] for i in range(4)]                                 # die vier schwachen Erzeuger (Quadrat -1)
KAP = [1, 1, 1, 1j]
SROT = [-0.5j*ga[1] @ ga[2], -0.5j*ga[2] @ ga[0], -0.5j*ga[0] @ ga[1]]
PL, PR = (I256 + g5)/2, (I256 - g5)/2
LL, LR = triple(SROT, PL), triple(SROT, PR)
ph5 = [p for p in (1, 1j, -1, -1j) if eq(g5, p*prod(ga, 256))][0]


def loetung(perm, signs, gam=None, eta=None, kap=None, innen=None):
    gam, eta, kap, innen = gam or ga, eta or ETA, kap or KAP, innen or GW
    g5l = chiral(gam)
    ph = [p for p in (1, 1j, -1, -1j) if eq(g5l, p*prod(gam, 256))][0]
    c = [1j*kap[i]*signs[i]*innen[perm[i]] for i in range(4)]
    g5w = ph*prod(c, 256)
    eps = [(gam[i] + 1j*c[i])/2 for i in range(4)]
    iot = [eta[i]*(gam[i] - 1j*c[i])/2 for i in range(4)]
    Nf = sum(eps[i] @ iot[i] for i in range(4))
    srot = [-0.5j*gam[1] @ gam[2], -0.5j*gam[2] @ gam[0], -0.5j*gam[0] @ gam[1]]
    sw = [-0.5j*c[1] @ c[2], -0.5j*c[2] @ c[0], -0.5j*c[0] @ c[1]]
    lorL, lorR = triple(srot, (I256 + g5l)/2), triple(srot, (I256 - g5l)/2)
    wL, wR = triple(sw, (I256 + g5w)/2), triple(sw, (I256 - g5w)/2)
    twL, twR = tuple(x + y for x, y in zip(lorL, wL)), tuple(x + y for x, y in zip(lorR, wR))
    sig = [gam[i] @ gam[j]/2 for i, j in PAIRS4]
    sigc = [c[i] @ c[j]/2 for i, j in PAIRS4]
    tw = [x + y for x, y in zip(sig, sigc)]
    par = prod([I256 - 2*eps[i] @ iot[i] for i in range(4)], 256)
    return {"c": c, "g5w": g5w, "eps": eps, "iot": iot, "N": Nf, "twL": twL, "twR": twR, "tw": tw, "sig": sig,
            "sigc": sigc, "wL": wL, "wR": wR, "par": par, "gam": gam, "eta": eta, "g5": g5l}


def vektor(X, rot, eta):
    """[rot_bc, X_a] = eta_c delta_ca X_b - eta_b delta_ba X_c für alle sechs Erzeuger?"""
    for n, (b, c_) in enumerate(PAIRS4):
        for a in range(4):
            if not eq(com(rot[n], X[a]), eta[c_]*(c_ == a)*X[b] - eta[b]*(b == a)*X[c_]):
                return False
    return True


L0 = loetung((0, 1, 2, 3), (1, 1, 1, 1))
LB = L0 if eq(L0["g5w"], CW) else loetung((0, 1, 2, 3), (1, 1, 1, -1))    # gleichsinnig: Lorentz-L mit SU(2)_L
LB_um = loetung((0, 1, 2, 3), (1, 1, 1, 1)) if LB is not L0 else loetung((0, 1, 2, 3), (1, 1, 1, -1))
SGN_REF = [1, 1, 1, 1] if LB is L0 else [1, 1, 1, -1]
ZAHLEN["1.6 Vorzeichen s_a der Referenzlötung"] = SGN_REF
GKD = g5 @ CW                                                        # Gamma_KD = gamma_5 x schwache Chiralität
SLL256, SLR256 = [kr(I8, x) for x in SL_L], [kr(I8, x) for x in SL_R]
eps, iot, Nf, TW = LB["eps"], LB["iot"], LB["N"], LB["tw"]
fock = (all(eq(acom(eps[i], iot[j]), (i == j)*I256) and null(acom(eps[i], eps[j])) and null(acom(iot[i], iot[j]))
            for i in range(4) for j in range(4)))
check(f"1.6 [E] Verdrillte Gruppe (Lesart B): Lötung phi: gamma_a -> c_a = i kappa_a s_a g^w_a mit c_a^2 = eta_a (kappa = 1 "
      f"räumlich, i zeitlich; Vorzeichen s = {SGN_REF}, so gewählt, dass Lorentz-L mit SU(2)_L gepaart ist; das andere "
      f"Vorzeichen paart L mit SU(2)_R); Sigma^tw_ab = Sigma_ab + c_a c_b / 2 bildet zwei vertauschende sl(2)-Tripel; der schwache "
      "Anteil des L-Tripels spannt SU(2)_L, der des R-Tripels SU(2)_R auf; chirale Lorentz-Tripel und die sechs "
      "Sigma_ab spannen dieselbe komplexe Algebra auf",
      eq(LB["g5w"], CW) and eq(LB_um["g5w"], -CW) and is_sl2(LB["twL"]) and is_sl2(LB["twR"])
      and all(null(com(x, y)) for x in LB["twL"] for y in LB["twR"])
      and crank(list(LB["wL"]) + SLL256) == 3 and crank(list(LB["wR"]) + SLR256) == 3
      and crank(list(LL) + list(LR) + LOR) == 6 and crank(list(LB["twL"]) + list(LB["twR"]) + TW) == 6)
CLt, CRt = cas(LB["twL"]), cas(LB["twR"])
P4, P4q = (I256 + CC)/2, (I256 - CC)/2
inh = {"K x 4": inhalt(PK @ P4, CLt, CRt), "K x 4-quer": inhalt(PK @ P4q, CLt, CRt),
       "K' x 4": inhalt(PKs @ P4, CLt, CRt), "K' x 4-quer": inhalt(PKs @ P4q, CLt, CRt),
       "T16 x 4": inhalt(PT16 @ P4, CLt, CRt), "T16 x 4-quer": inhalt(PT16 @ P4q, CLt, CRt)}
ZAHLEN["1.7 Inhalt unter verdrillter Spin(4) (Vielfachheit einschließlich der SU(4)-Dimension 4)"] = inh
say(f"      Inhalt unter der verdrillten Spin(4), Vielfachheiten einschließlich Farbfaktor 4: {inh}")
GERADE, UNGERADE = {"(1,1)": 8, "(1,3)": 4, "(3,1)": 4}, {"(2,2)": 8}
check("1.7 [E] z154 mit Matrizen nachgerechnet: K = Lambda^gerade x 4 + Lambda^ungerade x 4-quer (vier reduzierte "
      "Kähler-Dirac-Felder); das Komplement K' = S^+ cap P_+ trägt die anderen Hälften (Lambda^ungerade x 4 + "
      "Lambda^gerade x 4-quer), S^+ also Lambda* x (4 + 4-quer); T16 ist nicht nach Formengrad sortiert",
      inh["K x 4"] == GERADE and inh["K x 4-quer"] == UNGERADE and inh["K' x 4"] == UNGERADE
      and inh["K' x 4-quer"] == GERADE and inh["T16 x 4"] == {"(1,1)": 4, "(2,2)": 4, "(3,1)": 4}
      and inh["T16 x 4-quer"] == {"(1,1)": 4, "(1,3)": 4, "(2,2)": 4})
nK = {"K": vielfachheiten(PK, Nf, (0, 1, 2, 3, 4)), "K x 4": vielfachheiten(PK @ P4, Nf, (0, 1, 2, 3, 4)),
      "K x 4-quer": vielfachheiten(PK @ P4q, Nf, (0, 1, 2, 3, 4)), "S+": vielfachheiten(Sp, Nf, (0, 1, 2, 3, 4))}
ZAHLEN["1.8 Formengrad N_form: Vielfachheiten der Grade 0..4"] = nK
check(f"1.8 [E] Formengrad-Operator: eps_a = (gamma_a + i c_a)/2, iota_a = eta_a (gamma_a - i c_a)/2 erfüllen die "
      f"Fock-Relationen, N_form = Summe eps_a iota_a = 2 - (i/2) Summe eta_a gamma_a c_a vertauscht mit der verdrillten "
      f"Gruppe; Grade 0..4 auf K x 4: {nK['K x 4']}, auf K x 4-quer: {nK['K x 4-quer']}, auf S^+: {nK['S+']}; "
      f"(-1)^N_form = Gamma_KD",
      fock and eq(Nf, 2*I256 - 0.5j*sum(ETA[i]*ga[i] @ LB["c"][i] for i in range(4)))
      and all(null(com(Nf, x)) for x in TW) and nK["K x 4"] == [4, 0, 24, 0, 4] and nK["K x 4-quer"] == [0, 16, 0, 16, 0]
      and nK["S+"] == [8, 32, 48, 32, 8] and eq(LB["par"], GKD) and eq(LB_um["par"], -GKD))
check("1.9 [E] Gamma_KD = -Gamma_2 x Farbchiralität x Gamma_ges auf allen 256; auf S^+ also Gamma_KD = -Gamma_2 chi_c, "
      "auf K gleich chi_c: Der Reduktionsprojektor (1 + Gamma_KD chi_c)/2 des Kähler-Dirac-Felds ist auf S^+ genau der "
      "Indexprojektor P_-",
      eq(GKD, -Gam2 @ CC @ Gt) and eq(GKD @ PK, CC @ PK) and eq((I256 + GKD @ CC)/2 @ Sp, Pminus @ Sp))
befund("1", "ja", "E: Der Bestand ist exakt reproduziert (Ränge 128/64/64/32, 64 eps, so(3,2), Kommutanten 2/4/14). "
       "Lesart B aus z154 ist mit Matrizen bestätigt und um den Operator N_form ergänzt. N_form liegt nicht in "
       "so(4,2) + so(10): Es ist eine Summe gemischter Bilinearer gamma_a g^w_b, also die Lötung von Lorentz- und "
       "schwachen Richtungen selbst.")

# ============================================================================================================ Teil 2
say("\nTeil 2  Der Test")
say("   (a) Alle Operatoren, die Lorentz-Vektor und Eich-Singulett sind")
g8 = [hh[a] for a in LI]
g58 = chiral(g8)
gp8, gm8 = (hh[3] + hh[5])/2, (hh[3] - hh[5])/2
G28 = gm8 @ gp8 - gp8 @ gm8
I8c = np.eye(8, dtype=complex)
U8 = [I8c, g58, gp8, gm8, gp8 @ g58, gm8 @ g58, G28, G28 @ g58]
UN = ["1", "gamma_5", "gamma_+", "gamma_-", "gamma_+ gamma_5", "gamma_- gamma_5", "Gamma_2", "Gamma_2 gamma_5"]
rot8 = [g8[b] @ g8[c_]/2 for b, c_ in PAIRS4]
rows = []
for n, (b, c_) in enumerate(PAIRS4):
    adm = np.kron(rot8[n], np.eye(8)) - np.kron(np.eye(8), rot8[n].T)
    for a in range(4):
        blk = [np.zeros((64, 64), complex) for _ in range(4)]
        blk[a] = blk[a] + adm
        if c_ == a:
            blk[b] = blk[b] - ETA[c_]*np.eye(64)
        if b == a:
            blk[c_] = blk[c_] + ETA[b]*np.eye(64)
        rows.append(np.hstack(blk))
dim_vec8 = 256 - crank_mat(np.vstack(rows))
check(f"2.1 [E] Lorentz-Vektoren in End(8) (Spinor von Cl(4,2)): Lösungsraum der Vektorbedingung hat Dimension "
      f"{dim_vec8}; Basis gamma_a u mit u in {{1, gamma_5}} x {{1, gamma_+, gamma_-, Gamma_2}}",
      dim_vec8 == 8 and crank(U8) == 8 and all(vektor([g8[a] @ u for a in range(4)], rot8, ETA) for u in U8))
TRAEGER = {"K": PK, "K'": PKs, "S+": Sp, "K-": PKm, "K'-": PKsm, "S-": Sm, "T16": PT16}
PAARE = [("K", "K"), ("K", "K'"), ("K'", "K"), ("S+", "S+"), ("K", "K-"), ("K", "K'-"), ("S+", "S-"), ("T16", "T16")]
tab = {}
for nm in ("so(10)", "Pati-Salam", "Standardmodell"):
    ops = [kr(g8[0] @ u, c_) for u in U8 for c_ in C_G[nm]]
    tab[nm] = {"alle auf 256": crank(ops)}
    for src, tgt in PAARE:
        pr = [TRAEGER[tgt] @ X @ TRAEGER[src] for X in ops]
        neutral = [p for p in pr if eq(R @ p @ Rinv, p)]
        ungerade = [p for p in pr if eq(R @ p @ Rinv, -p) and np.any(p)]
        assert all(eq(R @ p @ Rinv, p) or eq(R @ p @ Rinv, -p) for p in pr)
        tab[nm][f"{src} -> {tgt}"] = [crank(pr), crank(neutral), crank(ungerade)]
    say(f"      {nm}: Zahl der Vektor-Multipletts [gesamt, uhrneutral, Uhrladung -1]: {tab[nm]}")
ZAHLEN["2.2 Lorentz-Vektor x Eich-Singulett je Trägerpaar [gesamt, uhrneutral, Uhrladung -1]"] = tab
check("2.2 [E] Auf K selbst gibt es keinen Operator, der Lorentz-Vektor und Singulett unter so(10) oder Pati-Salam "
      "ist (K -> K: 0). Alle solchen Operatoren auf S^+ (4 bzw. 8 Multipletts) sind unter Gamma_2 ungerade: Sie bilden "
      "K auf K' ab (gamma_a gamma_-, 2 bzw. 4) und zurück (gamma_a gamma_+). Von K nach S^- führt gamma_a selbst nach "
      "K^- = S^- cap P_- (2 bzw. 4), nichts nach S^- cap P_+",
      all(tab[k]["K -> K"][0] == 0 and tab[k]["K -> K'-"][0] == 0 for k in ("so(10)", "Pati-Salam"))
      and [tab["so(10)"][f"{a} -> {b}"][0] for a, b in PAARE] == [0, 2, 2, 4, 2, 0, 4, 2]
      and [tab["Pati-Salam"][f"{a} -> {b}"][0] for a, b in PAARE] == [0, 4, 4, 8, 4, 0, 8, 4]
      and tab["so(10)"]["alle auf 256"] == 16 and tab["Pati-Salam"]["alle auf 256"] == 32)
check("2.3 [E] Nur unter dem Standardmodell allein gibt es zwei Vektor-Multipletts K -> K: Sie sitzen auf dem "
      "rechtshändigen Neutrino (nu^c und sein Konjugiertes sind dieselbe Standardmodell-Darstellung) und tragen die "
      "Uhrladung -1; die Uhr verbietet sie",
      tab["Standardmodell"]["K -> K"] == [2, 0, 2] and tab["Standardmodell"]["alle auf 256"] == 112
      and np.diag(G11).real[leer]*np.diag(G11).real[voll] == -1)
check("2.4 [E] Uhrladungen auf der Spinorseite: R^2 = -1, also hat Ad R nur die Eigenwerte +1 und -1. Kein Operator auf "
      "den 256 trägt die Uhrladung -i des Viertelschritts; alle Vektor-Singuletts unter so(10) und Pati-Salam sind "
      "uhrneutral. Die Ladung -i sitzt allein im Markenfaktor A_q (z118)",
      eq(R @ R, -I256) and all(tab[k][f"{a} -> {b}"][2] == 0 for k in ("so(10)", "Pati-Salam") for a, b in PAARE))
befund("2a", "nein (gescopt)", "E: Einen nativen Operator 'auf K' mit Lorentz-Index gibt es nicht. Die Vektor-Singuletts "
       "wechseln alle den Indexsektor (K <-> K', das sind die Translationen und speziellen konformen Erzeuger gamma_a "
       "gamma_+-, aus denen Pi_a besteht) oder die Gesamtchiralität (gamma_a: K -> K^-). Keiner trägt die Uhrladung -i. "
       "Gilt für Singuletts unter so(10) oder Pati-Salam; unter dem Standardmodell allein kommen nur zwei uhrgeladene "
       "Neutrinoterme hinzu.")

say("   (b) D_perp, der induzierte Weyl-Operator und die Kähler-Dirac-Bedingungen")
QQ = float(np.exp(-2*np.pi))


def A_rechteck(N, q):
    S = np.zeros((N + 9, N + 9))
    for n in range(N + 8):
        S[n + 1, n] = 1
    A = S.T @ (np.eye(N + 9) - q*np.linalg.matrix_power(S, 4))
    return A[:N - 1, :N]


NM = 41
Aq = A_rechteck(NM, QQ)
zm = np.zeros(NM)
zm[::4] = [QQ**n for n in range(len(zm[::4]))]
zm /= np.linalg.norm(zm)
svA = np.linalg.svd(Aq, compute_uv=False)
ker_A, ker_Ad = NM - int(np.sum(svA > 1e-12)), (NM - 1) - int(np.sum(svA > 1e-12))
Rs, Rt = np.diag([1j**k for k in range(NM)]), np.diag([1j**k for k in range(NM - 1)])
ZAHLEN["2.5 Marken"] = {"dim ker A_q": ker_A, "dim ker A_q^+": ker_Ad, "|A_q zm|": float(np.round(np.abs(Aq @ zm).max(), 12)),
                        "Nullmode auf Marken 0 mod 4": bool(np.all(zm[np.arange(NM) % 4 != 0] == 0))}
check(f"2.5 [E/N] D_perp = gamma_- A_q + gamma_+ A_q^+ (z118): gamma_+ P_- = 0 und gamma_- P_+ = 0 exakt, A_q hat Index "
      f"{ker_A - ker_Ad} mit der Nullmode Summe q^n |4n> (Rest {np.abs(Aq @ zm).max():.0e}), R A_q R^-1 = -i A_q; der Kern "
      f"von D_perp auf S^+ x Marken ist K x Nullmode, sein Kokern K^- x Nullmode (z129). D_perp verschwindet auf dem "
      f"Indexkern",
      null(gp @ Pminus) and null(gm @ Pplus) and eq(gm @ Pminus, Pplus @ gm) and eq(gp @ Sp, Sm @ gp)
      and ker_A == 1 and ker_Ad == 0 and np.abs(Aq @ zm).max() < 1e-14 and np.allclose(Rt @ Aq @ np.linalg.inv(Rs), -1j*Aq))
dperp = {"[gamma_+-, Gamma_KD] = 0": null(com(gp, GKD)) and null(com(gm, GKD)),
         "[gamma_+-, N_form] = 0": null(com(gp, Nf)) and null(com(gm, Nf)),
         "[gamma_+-, Sigma^tw] = 0": all(null(com(g, x)) for g in (gp, gm) for x in TW),
         "{gamma_+-, Gamma_2} = 0": null(acom(gp, Gam2)) and null(acom(gm, Gam2)),
         "[gamma_+-, gamma_5] = 0": null(com(gp, g5)) and null(com(gm, g5))}
ZAHLEN["2.6 D_perp"] = dperp
check("2.6 [E] D_perp gegen die drei Bedingungen: gamma_+- vertauschen mit Gamma_KD (D_perp ist Kähler-Dirac-gerade, "
      "nicht ungerade), mit N_form (Formengrad 0, nicht +-1) und mit der verdrillten Gruppe (Skalar, kein Vektor). "
      "D_perp antikommutiert nur mit der Indexgradierung Gamma_2. Das gilt für jede Lötung, weil gamma_+- mit jedem "
      "gamma_a g^w_b vertauschen", all(dperp.values())
      and all(null(com(g, ga[i] @ GW[j])) for g in (gp, gm) for i in range(4) for j in range(4)))
weyl = {"{gamma_a, Gamma_KD} = 0": all(null(acom(g, GKD)) for g in ga),
        "gamma_a = eps_a + eta_a iota_a": all(eq(ga[i], eps[i] + ETA[i]*iot[i]) for i in range(4)),
        "[N_form, eps_a] = eps_a, [N_form, iota_a] = -iota_a": all(eq(com(Nf, eps[i]), eps[i]) and eq(com(Nf, iot[i]), -iot[i])
                                                                 for i in range(4)),
        "gamma_a Vektor unter Sigma^tw": vektor(ga, TW, ETA),
        "gamma_a: K -> K^-": all(eq(PKm @ g @ PK, g @ PK) for g in ga)}
ZAHLEN["2.7 induzierter Weyl-Operator gamma^a p_a"] = weyl
check("2.7 [E] Der induzierte vierdimensionale Operator gamma^a p_a: K x Nullmode -> K^- x Nullmode (Kern nach Kokern "
      "von D_perp) erfüllt alle drei Bedingungen: Er antikommutiert mit Gamma_KD, ist eps^a p_a + eta iota^a p_a mit "
      "Formengrad +1 und -1 (also d - delta bis auf Vorzeichenkonvention), und gamma_a ist ein Vektor der "
      "verdrillten Gruppe", all(weyl.values()))
# Tautologieprobe: andere Lötungen
proben = {"umgekehrte Orientierung (Lorentz-L mit SU(2)_R)": LB_um,
          "andere Zuordnung der Richtungen": loetung((2, 0, 3, 1), (1, -1, 1, 1)),
          "Lötung an vier Farbrichtungen (innere Zelle, Art A aus z154)": loetung((0, 1, 2, 3), (1, 1, 1, 1),
                                                                                innen=[G16[6 + i] for i in range(4)])}
taut = {}
for nm, Lx in proben.items():
    taut[nm] = {"Fock": all(eq(acom(Lx["eps"][i], Lx["iot"][j]), (i == j)*I256) for i in range(4) for j in range(4)),
                "Grad +-1": all(eq(ga[i], Lx["eps"][i] + ETA[i]*Lx["iot"][i]) and eq(com(Lx["N"], Lx["eps"][i]), Lx["eps"][i])
                                for i in range(4)),
                "antikommutiert mit (-1)^N": all(null(acom(g, Lx["par"])) for g in ga),
                "Vektor der verdrillten Gruppe": vektor(ga, Lx["tw"], ETA),
                "D_perp Grad 0": null(com(gp, Lx["N"])) and null(com(gm, Lx["N"])),
                "N auf K (Grade 0..4)": vielfachheiten(PK, Lx["N"], (0, 1, 2, 3, 4)),
                "N vertauscht mit dem Standardmodell (Zahl der 12 Erzeuger)":
                    sum(1 for T in SM12 if null(com(kr(I8, T), Lx["N"])))}
ZAHLEN["2.8 Tautologieprobe: andere Lötungen"] = taut
for nm in proben:
    say(f"      {nm}: {taut[nm]}")
check("2.8 [E] Tautologie: Dieselben drei Eigenschaften gelten wörtlich für jede Lötung, auch für die umgekehrte "
      "Orientierung, für eine andere Zuordnung der Richtungen und für eine Lötung an vier Farbrichtungen (die in z154 "
      "ausgeschlossene innere Zelle). Sie folgen aus {gamma_a, g_b} = 0 in Cl(4,12), nicht aus dem Operator. Jeder "
      "Dirac-Operator mit Eich-Singulett-Symbol ist unter jeder verdrillten Gruppe d - delta",
      all(t["Fock"] and t["Grad +-1"] and t["antikommutiert mit (-1)^N"] and t["Vektor der verdrillten Gruppe"]
          and t["D_perp Grad 0"] and t["N auf K (Grade 0..4)"] == [4, 16, 24, 16, 4] for t in taut.values()))
adN2 = {}
for nm, u in (("gamma_a gamma_-", gm), ("gamma_a gamma_+", gp), ("gamma_a gamma_- gamma_5", gm @ g5),
              ("gamma_a gamma_+ gamma_5", gp @ g5), ("Pi_a", None)):
    X = PI[0] if u is None else ga[0] @ u
    adN2[nm] = {"Gamma_KD-ungerade": null(acom(X, GKD)), "reiner Grad +-1": eq(com(Nf, com(Nf, X)), X),
                "Vektor verdrillt": vektor([PI[a] if u is None else ga[a] @ u for a in range(4)], TW, ETA)}
ZAHLEN["2.9 Vektor-Singuletts auf S^+ gegen N_form"] = adN2
check("2.9 [E] Dasselbe für die Vektor-Singuletts auf S^+: gamma_a gamma_+- und die so(3,2)-Translationen Pi_a sind "
      "Gamma_KD-ungerade mit reinem Formengrad +-1 und Vektoren der verdrillten Gruppe; die Varianten mit gamma_5 sind "
      "Gamma_KD-ungerade, aber mit Hodge-Stern (Grad p -> 4 - p +- 1). Auch das gilt für jede Lötung",
      all(v["Gamma_KD-ungerade"] and v["Vektor verdrillt"] for v in adN2.values())
      and [v["reiner Grad +-1"] for v in adN2.values()] == [True, True, False, False, True])
# Was die verdrillte Gruppe zusätzlich zuließe: Zählung über Charaktere
FORMEN = Counter({(1, 1): 2, (3, 1): 1, (1, 3): 1, (2, 2): 2})


def su2x(a, c_):
    return list(range(abs(a - c_) + 1, a + c_, 2))


def tens(r1, r2):
    res = Counter()
    for (a, b), m1 in r1.items():
        for (c_, d), m2 in r2.items():
            for x in su2x(a, c_):
                for y in su2x(b, d):
                    res[(x, y)] += m1*m2
    return res


n_tw = 4*2*tens(FORMEN, FORMEN)[(2, 2)]                               # End(8) = 4 Lambda*, End(32)^SU(4) = 2 Lambda*
n_prod = 4*FORMEN[(2, 2)]*2*FORMEN[(1, 1)]
echt = {"c_a Vektor verdrillt": vektor(LB["c"], TW, ETA), "eps_a Vektor verdrillt": vektor(eps, TW, ETA),
        "eta_a iota_a Vektor verdrillt": vektor([ETA[i]*iot[i] for i in range(4)], TW, ETA),
        "eps_a Lorentz-Vektor": vektor(eps, LB["sig"], ETA),
        "c_a schwaches Singulett": all(null(com(x, kr(I8, T))) for x in LB["c"] for T in SL_L)}
ZAHLEN["2.10 Vektoren der verdrillten Gruppe mit SU(4)-Invarianz auf 256"] = {"verdrillt": n_tw, "davon Lorentz-Vektor x schwaches Singulett": n_prod, **echt}
check(f"2.10 [E] Was wirklich entscheidet: Unter verdrillter Gruppe x SU(4) gibt es auf den 256 genau {n_tw} "
      f"Vektor-Multipletts (Charaktere), davon sind {n_prod} Lorentz-Vektor x schwaches Singulett (= Zählung 2.2 für "
      f"Pati-Salam). Die übrigen {n_tw - n_prod} enthalten c_a, eps_a und eta_a iota_a (Index unten) einzeln: d und delta getrennt sind "
      f"Vektoren der verdrillten Gruppe, aber keine Lorentz-Vektoren. Eine Stütze für P_KD wäre ein nativer Operator in "
      f"diesen {n_tw - n_prod}; die nativen (2.2) liegen alle in den {n_prod}",
      n_tw == 128 and n_prod == 32 == tab["Pati-Salam"]["alle auf 256"] and echt["c_a Vektor verdrillt"]
      and echt["eps_a Vektor verdrillt"] and echt["eta_a iota_a Vektor verdrillt"] and not echt["eps_a Lorentz-Vektor"]
      and not echt["c_a schwaches Singulett"])
befund("2b", "nein (gescopt) / tautologisch", "E: Der native Indexoperator D_perp ist nicht von der Form d - delta: Er "
       "ist Kähler-Dirac-gerade, hat Formengrad 0 und ist ein Skalar der verdrillten Gruppe, für jede Lötung. Der "
       "induzierte Weyl-Operator gamma^a p_a und die Translationen Pi_a erfüllen alle drei Bedingungen, aber für jede "
       "Lötung gleich, auch für die in z154 ausgeschlossene. Die Aussage 'der Operator ist d - delta' trägt keine "
       "Information über die Gruppe. Entscheiden würde ein natives Datum außerhalb von Lorentz-Vektor x Singulett "
       "(d und delta getrennt, die zweite Clifford-Wirkung c_a oder N_form selbst).")

say("   (c) Unterscheidet ein natives Datum die verdrillte von der unverdrillten Gruppe?")
# Zerlegung so(16) = so(4,2) + so(10) + gemischt; Uhrladung
BIL_A = [G16[i] @ G16[j]/2 for i, j in itertools.combinations(range(6), 2)]
BIL_B = [G16[i] @ G16[j]/2 for i, j in itertools.combinations(range(6, 16), 2)]
BIL_M = [G16[i] @ G16[j]/2 for i in range(6) for j in range(6, 16)]
dims16 = [crank(BIL_A), crank(BIL_B), crank(BIL_M)]
uhr_bil = (all(eq(R @ X @ Rinv, X) for X in BIL_A + BIL_B), all(eq(R @ X @ Rinv, -X) for X in BIL_M))
check(f"2.11 [E] so(16) = so(4,2) + so(10) + gemischt mit Dimensionen {dims16}; die Uhr wirkt auf so(4,2) + so(10) "
      f"trivial und auf allen 60 gemischten Erzeugern mit -1. N_form - 2 liegt ganz im gemischten Teil: "
      f"R N_form R^-1 = 4 - N_form. Die Uhr tauscht eps_a und eta_a iota_a (d und delta); auf K ist sie i gamma_5, "
      f"der Hodge-Stern p -> 4 - p mal i. Der Formengrad ist keine uhrinvariante Größe, nur (-1)^p und |p - 2|",
      dims16 == [15, 45, 60] and all(uhr_bil) and crank(BIL_M + [Nf - 2*I256]) == 60
      and eq(R @ Nf @ Rinv, 4*I256 - Nf) and all(eq(R @ eps[i] @ Rinv, ETA[i]*iot[i]) for i in range(4))
      and eq(g5 @ Nf @ g5, 4*I256 - Nf) and null(com(R, GKD)))
check("2.12 [E] Unter der Lötung überlebt von Uhr und Zentrum von SU(2)_L nur das Produkt: [R, N_form] != 0, "
      "[c_L, N_form] != 0, [R^-1 c_L, N_form] = 0, und auf K ist i Gamma_KD = R^-1 c_L (z154: Kähler-Dirac-Z_4 = "
      "Uhr^-1 mal Zentrum von SU(2)_L). P_KD als Lötung bricht also die Uhr selbst auf diese Kombination",
      np.any(com(R, Nf)) and np.any(com(cL, Nf)) and null(com(Rinv @ cL, Nf)) and eq(1j*GKD @ PK, Rinv @ cL @ PK))


def zentralisator(alg, objs):
    return len(alg) - crank([np.concatenate([com(x, o) for o in objs], axis=1) for x in alg])


SM256, SO10_256 = [kr(I8, T) for T in SM12], [kr(I8, T) for T in SO10]
zen = {"von N_form im Standardmodell": zentralisator(SM256, [Nf]), "der verdrillten Gruppe im Standardmodell": zentralisator(SM256, TW),
       "von N_form in so(10)": zentralisator(SO10_256, [Nf]), "der verdrillten Gruppe in so(10)": zentralisator(SO10_256, TW),
       "von N_form in so(3,2)": zentralisator(SIG32, [Nf])}
ZAHLEN["2.13 Zentralisatoren (Dimension)"] = zen
check(f"2.13 [E] Eichgruppe: Mit N_form (und mit der verdrillten Gruppe) vertauschen vom Standardmodell genau die acht "
      f"Erzeuger der Farbe, von so(10) die 15 von so(6), von so(3,2) keiner: {zen}. SU(2)_L und die Hyperladung wären "
      f"keine inneren Symmetrien mehr",
      zen == {"von N_form im Standardmodell": 8, "der verdrillten Gruppe im Standardmodell": 8,
              "von N_form in so(10)": 15, "der verdrillten Gruppe in so(10)": 15, "von N_form in so(3,2)": 0})
# reelle Form
SIGW = [GW[i] @ GW[j]/2 for i, j in PAIRS4]
dR = TW + [1j*x for x in TW]
r_lor = 12 + 12 - rrank(dR + LOR + SIGW)
antiherm = [eq(x.conj().T, -x) for x in TW]
# euklidisches Paar (Uhrpaar E aus z118: Richtungen 4, 5; Querraum 0, 1, 2, 3)
gpE, gmE = (G16[4] + 1j*G16[5])/2, -(G16[4] - 1j*G16[5])/2
PplusE, PminusE = gmE @ gpE, gpE @ gmE
Gam2E = PplusE - PminusE
gaE = [G16[i] for i in range(4)]
PKE = Sp @ PminusE
LE0 = loetung((0, 1, 2, 3), (1, 1, 1, 1), gam=gaE, eta=[1, 1, 1, 1], kap=[1, 1, 1, 1])
LE = LE0 if eq(LE0["g5w"], CW) else loetung((0, 1, 2, 3), (1, 1, 1, -1), gam=gaE, eta=[1, 1, 1, 1], kap=[1, 1, 1, 1])
dRE = LE["tw"] + [1j*x for x in LE["tw"]]
r_eu = 12 + 12 - rrank(dRE + LE["sig"] + SIGW)
inhE = {"K_E x 4": inhalt(PKE @ P4, cas(LE["twL"]), cas(LE["twR"])), "K_E x 4-quer": inhalt(PKE @ P4q, cas(LE["twL"]), cas(LE["twR"]))}
schnitt = 6 + 10 - crank(LE["sig"] + SIG32)
ZAHLEN["2.14 reelle Form"] = {"reelle Dimension verdrillt cap (so(3,1) + so(4)), Dilatationspaar": r_lor,
                              "antihermitesche verdrillte Erzeuger (lorentzsch)": antiherm,
                              "reelle Dimension, euklidisches Uhrpaar": r_eu, "Inhalt K_E": inhE,
                              "dim so(4)_E cap so(3,2)": schnitt}
check(f"2.14 [E] Reelle Form: Die verdrillte Algebra ist im nativen lorentzschen Rahmen (Nullpaar der Dilatation, "
      f"Querraum (3,1)) nur komplex definiert. Ihr Schnitt mit so(3,1) + so(4)_schwach hat reelle Dimension {r_lor}: "
      f"nur die Drehungen so(3), diagonal mit SU(2)_(L+R); die drei Boost-Erzeuger sind Boost + i mal kompakt. (so(3,1) "
      f"ist einfach und nicht kompakt, es gibt keinen Homomorphismus in su(2) + su(2).) Im euklidischen Rahmen des "
      f"komplexen Uhrpaars E (Querraum (4,0)) ist die Dimension {r_eu}, und K_E trägt denselben Inhalt; diese so(4) "
      f"schneidet die native so(3,2) aber wieder nur in so(3) (Dimension {schnitt})",
      r_lor == 3 and antiherm == [True, True, False, True, False, False] and r_eu == 6
      and all(eq(x.conj().T, -x) for x in LE["tw"]) and inhE == {"K_E x 4": GERADE, "K_E x 4-quer": UNGERADE}
      and int(np.trace(PKE).real) == 64 and schnitt == 3 and eq(LE["g5w"], CW))
# Realstruktur J (z116/z118), explizit: Produkt der reellen Erzeuger
Ca = prod([g for g in hh if eq(g.conj(), g)], 8)
Cb = prod([g for g in Gi if eq(g.conj(), g)], 32)
C16 = np.kron(Ca, Cb)
C16i = np.linalg.inv(C16)
C16i = np.round(C16i.real) + 1j*np.round(C16i.imag)


def Jc(X):
    return C16 @ X.conj() @ C16i


J_ok = (eq(C16 @ C16i, I256) and all(eq(Jc(g), g) for g in G16) and eq(C16 @ C16.conj(), I256) and eq(Jc(Gt), Gt)
        and eq(Jc(Gam2), Gam2) and eq(Jc(R), R) and eq(Jc(PK), PK))
n_t = eps[3] @ iot[3]
jtab = {"J Sigma J^-1 = Sigma (Lorentz reell)": all(eq(Jc(x), x) for x in LOR),
        "J T J^-1 = T (schwach reell)": all(eq(Jc(x), x) for x in SIGW),
        "dim Spann(verdrillt, J verdrillt J^-1)": crank(TW + [Jc(x) for x in TW]),
        "J gamma_5 J^-1 = -gamma_5": eq(Jc(g5), -g5), "J Gamma_KD J^-1 = -Gamma_KD": eq(Jc(GKD), -GKD),
        "J N_form J^-1 = N_form + 1 - 2 n_t": eq(Jc(Nf), Nf + I256 - 2*n_t),
        "J Gamma_2^E J^-1 = -Gamma_2^E": eq(Jc(Gam2E), -Gam2E), "J K_E = K'_E": eq(Jc(PKE), Sp @ PplusE),
        "dim Spann(verdrillt_E, J verdrillt_E J^-1)": crank(LE["tw"] + [Jc(x) for x in LE["tw"]])}
ZAHLEN["2.15 Realstruktur J"] = jtab
check(f"2.15 [E] Realstruktur J (J^2 = 1, vertauscht mit allen 16 Erzeugern, erhält S^+, K, Gamma_2 und die Uhr): J "
      f"lässt die reellen Lorentz- und schwachen Erzeuger fest, bildet die verdrillte Algebra aber auf die konjugierte "
      f"Verdrillung ab (gemeinsamer Spann Dimension {jtab['dim Spann(verdrillt, J verdrillt J^-1)']} statt 6; J-stabil "
      f"sind nur die drei Drehungen). J tauscht gerade und ungerade Formen (J Gamma_KD J^-1 = -Gamma_KD, Dirac-"
      f"Konjugation mit der Zeitrichtung). Im euklidischen Rahmen ist die verdrillte Gruppe J-stabil, aber der Träger "
      f"K_E nicht (J tauscht K_E mit K'_E): In keinem der beiden Rahmen sind Träger und volle verdrillte Gruppe "
      f"zugleich reell",
      J_ok and jtab["J Sigma J^-1 = Sigma (Lorentz reell)"] and jtab["J T J^-1 = T (schwach reell)"]
      and jtab["dim Spann(verdrillt, J verdrillt J^-1)"] == 9 and jtab["J gamma_5 J^-1 = -gamma_5"]
      and jtab["J Gamma_KD J^-1 = -Gamma_KD"] and jtab["J N_form J^-1 = N_form + 1 - 2 n_t"]
      and jtab["J Gamma_2^E J^-1 = -Gamma_2^E"] and jtab["J K_E = K'_E"]
      and jtab["dim Spann(verdrillt_E, J verdrillt_E J^-1)"] == 6)
# so(3,2): schließt auf die unverdrillte Lorentz-Algebra
erz = TW + PI
komm = [com(x, y) for x, y in itertools.combinations(erz, 2)]
d_gen = crank(erz + komm)
d_abg = crank(erz + komm + [com(x, y) for x in SIG32 + SIGW for y in SIG32 + SIGW])
gem_wf = [GW[i] @ G16[6 + j]/2 for i in range(4) for j in range(6)]   # (4,6): schwach x Farbe
su3_256 = [kr(I8, T) for T in SU3]
inv_46 = zentralisator(gem_wf, su3_256)
farb6 = [G16[6 + j] for j in range(6)]
inv_6 = len(farb6) - crank([np.concatenate([com(x, o) for o in su3_256], axis=1) for x in farb6])
ZAHLEN["2.16 so(3,2) und Verdrillung"] = {"dim der von Sigma^tw und Pi_a erzeugten Algebra": d_gen,
                                          "dim so(3,2) + so(4)_schwach": d_abg, "su(3)-Invarianten in (4,6)": inv_46,
                                          "su(3)-Invarianten im Farbvektor 6": inv_6}
check(f"2.16 [E] Native Gravitationsalgebra: Die Translationen Pi_a sind Vektoren unter Sigma und unter Sigma^tw, "
      f"aber [Pi_a, Pi_b] = Sigma_ab ist die unverdrillte Lorentz-Algebra. Sigma^tw und Pi_a erzeugen so(3,2) + "
      f"so(4)_schwach (Dimension {d_gen}), keine verdrillte so(3,2). Verdrillte Translationen Pi_a + x_a im direkten "
      f"Produkt so(4,2) + so(10) bräuchten x_a im Teil (4,6) von so(10), also einen Farbvektor; su(3) hat dort keine Invariante ({inv_46} in (4,6), "
      f"{inv_6} in der 6). Die so(3,2) lässt sich nicht farbneutral verdrillen",
      vektor(PI, LOR, ETA) and vektor(PI, TW, ETA) and d_gen == 16 and d_abg == 16 and inv_46 == 0 and inv_6 == 0
      and crank(gem_wf) == 24)
# Gewichte: Mischspur vierter Ordnung auf K
tests = {"delta delta x delta delta": (0, 0, 0, 0), "eps x delta delta": (0, 5, 0, 0),
         "delta delta x eps": (0, 0, 0, 5), "eps x eps": (0, 5, 0, 5)}
GEW = {"1": I256, "Gamma_2": Gam2, "R": R, "chi_w (= Zentrum von SU(2)_R)": CW, "Gamma_KD": GKD}
misch = {}
for wn, Wm in GEW.items():
    misch[wn] = {}
    for tn, (p1, p2, p3, p4) in tests.items():
        v = complex(np.trace(PK @ Wm @ LB["sig"][p1] @ LB["sig"][p2] @ LB["sigc"][p3] @ LB["sigc"][p4]))
        misch[wn][tn] = (str(int(v.real)) if v.imag == 0 else f"{int(v.imag)}i")
    say(f"      Tr_K(W Sigma Sigma Sigma^w Sigma^w), W = {wn}: {misch[wn]}")
zweite = {"Tr_K(R Sigma_01 Sigma_24)": int(np.trace(PK @ R @ LB["sig"][0] @ LB["sig"][5]).real),
          "Tr_K(R Sigma^tw_01 Sigma^tw_24)": int(np.trace(PK @ R @ TW[0] @ TW[5]).real),
          "Tr_K(Gamma_KD Sigma^tw_01 Sigma^tw_24)": int(np.trace(PK @ GKD @ TW[0] @ TW[5]).real),
          "Tr_K(Gamma_KD)": int(np.trace(PK @ GKD).real), "Tr_K(-Gamma_2)": int(np.trace(-PK @ Gam2).real)}
ZAHLEN["2.17 Gewichte"] = {"Mischspur vierter Ordnung": misch, "zweite Ordnung": zweite}
check(f"2.17 [E] Gewichte: Die Euler-Struktur eps x eps (die Krümmungsform der Kähler-Dirac-Anomalie [L]) erscheint in "
      f"der Mischspur Tr_K(W Sigma Sigma Sigma^w Sigma^w) nur mit dem Gewicht Gamma_KD; die nativen Gewichte 1, "
      f"Gamma_2 (z118/z127) geben delta delta x delta delta, R gibt eps x delta delta. In zweiter Ordnung ändert die "
      f"Verdrillung die Uhrspur nicht (16 = 64/4 mit und ohne Verdrillung), und mit Gamma_KD verschwindet sie: "
      f"{zweite}. Der Index von D_perp bezüglich Gamma_2 ist 64, die Gamma_KD-Spur auf K ist 32 - 32 = 0",
      misch["1"] == {"delta delta x delta delta": "4", "eps x delta delta": "0", "delta delta x eps": "0", "eps x eps": "0"}
      and misch["Gamma_2"]["delta delta x delta delta"] == "-4" and misch["Gamma_2"]["eps x eps"] == "0"
      and misch["R"] == {"delta delta x delta delta": "0", "eps x delta delta": "-4", "delta delta x eps": "0", "eps x eps": "0"}
      and misch["Gamma_KD"] == {"delta delta x delta delta": "0", "eps x delta delta": "0", "delta delta x eps": "0", "eps x eps": "-4"}
      and misch["chi_w (= Zentrum von SU(2)_R)"]["delta delta x eps"] == "4i"
      and zweite == {"Tr_K(R Sigma_01 Sigma_24)": 16, "Tr_K(R Sigma^tw_01 Sigma^tw_24)": 16,
                     "Tr_K(Gamma_KD Sigma^tw_01 Sigma^tw_24)": 0, "Tr_K(Gamma_KD)": 0, "Tr_K(-Gamma_2)": 64}
      and eq(CW, -cL))
# Skalare im Indexkern nur unter der verdrillten Gruppe; die Uhr ist kein verdrilltes Gruppenelement
l0 = lagrange((0, 3, 8), 0)
powL, powR = [I256, CLt, CLt @ CLt], [I256, CRt, CRt @ CRt]
Psing = sum(l0[i].numerator/l0[i].denominator*powL[i] for i in range(3)) @ sum(l0[j].numerator/l0[j].denominator*powR[j] for j in range(3))
sing = {"verdrillte Singuletts in K": inhalt(PK, CLt, CRt)["(1,1)"],
        "Lorentz-Singuletts in K": inhalt(PK, cas(LL), cas(LR)).get("(1,1)", 0),
        "Spur der Uhr auf den verdrillten Singuletts": complex(np.round(np.trace(Psing @ PK @ R), 9)) == 0,
        "Spur von R^2": int(np.round(np.trace(Psing @ PK @ R @ R).real))}
ZAHLEN["2.18 Skalare im Indexkern"] = sing
check(f"2.18 [E/N] Nullmode und Indexkern: Die Nullmode von A_q ist ein Markenvektor und unter beiden Gruppen ein "
      f"Singulett. K enthält unter Lorentz keinen Skalar, unter der verdrillten Gruppe {sing['verdrillte Singuletts in K']} "
      f"(0- und 4-Formen mal 4). Auf ihnen hat die Uhr die Eigenwerte +i und -i je viermal (Spur 0, Spur von R^2 gleich "
      f"{sing['Spur von R^2']}): Die Uhr ist kein Element der verdrillten Gruppe mal Zentrum der Farbe, denn jedes "
      f"solche Element wirkt dort als ein Skalar",
      sing == {"verdrillte Singuletts in K": 8, "Lorentz-Singuletts in K": 0,
               "Spur der Uhr auf den verdrillten Singuletts": True, "Spur von R^2": -8})
befund("2c", "nein (gescopt)", "E: Kein natives Datum zeichnet die verdrillte Gruppe aus. Uhr, Gamma_2, gamma_5, Index, "
       "Nullmode und alle Vektor-Singuletts sind unter Lorentz und schwacher Spin(4) getrennt kovariant und damit blind. "
       "Vier native Daten sprechen gegen die verdrillte Gruppe als Drehgruppe der Naht: (1) Die Uhr tauscht d und delta "
       "und kehrt den Formengrad um; die Lötung bricht R auf R^-1 c_L. (2) Lorentzsch existiert die Verdrillung reell "
       "nur für so(3); J bildet sie auf die konjugierte Verdrillung ab. (3) Die Translationen der nativen so(3,2) "
       "schließen auf die unverdrillte Lorentz-Algebra, und eine verdrillte so(3,2) im direkten Produkt verlangt einen Farbvektor. (4) Die "
       "Lötung lässt vom Standardmodell nur die Farbe. Die Indexgradierung ist nicht (-1)^Formengrad, sondern "
       "Gamma_2 = -Gamma_KD chi_c auf S^+. Gilt für konstante Lötung im Kontinuum; ein Gitter mit nur diskreter "
       "Drehgruppe ist davon nicht erfasst (Teil 3).")

# ============================================================================================================ Teil 3
say("\nTeil 3  Gitterseite")
uhr_zentral = all(null(com(R, X)) for X in LOR + SO10_256 + SIG32)
g01 = (I256 + ga[0] @ ga[1]) @ (I256 + LB["c"][0] @ LB["c"][1])/2    # verdrillte Vierteldrehung in der Ebene (0, 1)
g01sq = g01 @ g01
spuren = [int(np.trace(PK @ M).real) for M in (I256, g01, g01sq, g01sq @ g01)]
n1 = (spuren[0] + 2*spuren[1] + spuren[2])//4
nm1 = (spuren[0] - 2*spuren[1] + spuren[2])//4
ni = (spuren[0] - spuren[2])//4
eig_g = {"1": n1, "-1": nm1, "i": ni, "-i": ni}
eig_R = {"i": int((64 + (np.trace(PK @ R)/1j).real)//2), "-i": int((64 - (np.trace(PK @ R)/1j).real)//2)}
lor_roh = I256 + ga[0] @ ga[1]
ZAHLEN["3.1 Vierteldrehungen auf K"] = {"verdrillte Gitterdrehung: Eigenwerte": eig_g, "Uhr: Eigenwerte": eig_R}
check(f"3.1 [E] Native Vierteldrehung: Die Uhr R = i Gamma_11 vertauscht mit allen Lorentz-, so(3,2)- und "
      f"so(10)-Erzeugern; sie dreht weder Lorentz- noch schwache Indizes. Eine Gitter-Vierteldrehung der verdrillten "
      f"Gruppe, g = exp((pi/2) Sigma^tw_01), dreht gamma_0 -> gamma_1 und zugleich c_0 -> c_1, hat g^4 = 1 (die "
      f"unverdrillte hat -1) und auf K die Eigenwerte {eig_g}; die Uhr hat {eig_R}. Die Vierteldrehung P_sym des "
      f"Transfers wirkt auf den Nahtindex, die Marken auf den Markenfaktor: beide vertauschen mit der Spinorseite "
      f"(direktes Produkt, z150)",
      uhr_zentral and eq(g01 @ ga[0], -ga[1] @ g01) and eq(g01 @ LB["c"][0], -LB["c"][1] @ g01)
      and eq(g01sq @ g01sq, I256) and eq(np.linalg.matrix_power(lor_roh, 4), -4*I256)
      and eig_g == {"1": 32, "-1": 0, "i": 16, "-i": 16} and eig_R == {"i": 32, "-i": 32} and null(com(g01, Nf)))
ord_c, ord_R = gruppenordnung(GW), gruppenordnung([R])
check(f"3.2 [E] Gruppenstruktur: Auf dem gestaffelten Gitter ist die Verschiebung um einen Platz die zweite "
      f"Clifford-Wirkung (Geschmacksmatrix) [L]; vier antikommutierende Einheitsschritte erzeugen eine nichtabelsche "
      f"Gruppe der Ordnung {ord_c}. Die native diskrete Gruppe einer Naht ist die zyklische Uhr der Ordnung {ord_R}; ihr "
      f"Erzeuger ist zentral. Vier antikommutierende native Viertelschritte gibt es nicht (ein Viertelschritt S^+, "
      f"Uhrladung -i)",
      ord_c == 32 and ord_R == 4 and all(null(acom(GW[i], GW[j])) for i in range(4) for j in range(i)))


def A_wilson(N, w):
    S = np.zeros((N, N))
    for n in range(N - 1):
        S[n + 1, n] = 1
    return (S.T - w*np.eye(N)) @ (np.eye(N) - QQ*np.linalg.matrix_power(S, 4))


def windung(w, M=4096):
    z = np.exp(2j*np.pi*np.arange(M)/M)
    a = (1/z - w)*(1 - QQ*z**4)
    ph = np.unwrap(np.angle(np.append(a, a[0])))
    return int(round((ph[-1] - ph[0])/(2*np.pi)))


ecken = sorted(itertools.product((0, 1), repeat=3))
w_ecke = {e: 2*sum(e) for e in ecken}                                 # w = Summe (1 - cos k_i), k_i in {0, pi}
grad = Counter(sum(e) for e in ecken)
wil = {}
for w in (0, 2, 4, 6):
    sv = np.linalg.svd(A_wilson(60, w), compute_uv=False)
    wil[w] = {"Ecken": sum(1 for e in ecken if w_ecke[e] == w), "Windung": windung(w), "Index": -windung(w),
              "kleinster Singulärwert": float(np.round(sv[-1], 6)), "Schranke (w-1)(1-q)": float(np.round((w - 1)*(1 - QQ), 6))}
    say(f"      w = {w}: {wil[w]}")
ZAHLEN["3.3 Wilson-Fortsetzung an den acht Ecken"] = {str(k): v for k, v in wil.items()}
moden_wilson = sum(v["Ecken"]*v["Index"] for v in wil.values())
check(f"3.3 [E/N] Wilson-Fortsetzung A_q(k) = (S^+ - w(k))(1 - q S^4), w = Summe (1 - cos k_i): An den acht Ecken ist "
      f"w = 2 p mit dem Eckengrad p = Zahl der Komponenten pi (Vielfachheiten {[grad[p] for p in range(4)]} wie "
      f"Lambda^p(R^3)); der Index ist 1 für |w| < 1 und 0 sonst (Windung des Symbols), also genau {moden_wilson} "
      f"Nullmode bei p = 0; die sieben anderen Ecken haben Singulärwerte >= (w - 1)(1 - q). Der Wilson-Term ist an den "
      f"Ecken das Doppelte des Formengrads im Geschmacksraum: Er hebt alles außer der 0-Form",
      [grad[p] for p in range(4)] == [1, 3, 3, 1] and moden_wilson == 1 and wil[0]["Index"] == 1
      and all(wil[w]["Index"] == 0 and wil[w]["kleinster Singulärwert"] >= wil[w]["Schranke (w-1)(1-q)"] for w in (2, 4, 6))
      and wil[0]["kleinster Singulärwert"] == 0 and [windung(x) for x in (0.5, 0.99, 1.01)] == [-1, -1, 0])
check("3.4 [E] Gegenprobe: Ohne Wilson-Term (w = 0 an allen Ecken, naive oder gestaffelte Fortsetzung) trüge jede der "
      "acht Ecken eine Nullmode; diese Entartung ist der Geschmacksraum, auf dem Gitterdrehungen die Ecken permutieren "
      "und so Spin und Geschmack verriegeln [L]. Die Verschiebung k_i -> k_i + pi (Geschmacksmatrix im Impulsraum) "
      "ändert w um 2 cos k_i und ist in der Wilson-Fortsetzung keine Symmetrie; w ist unter der Würfelgruppe "
      "invariant, die masselose Ecke k = 0 ist ihr Fixpunkt und trägt keinen Geschmacksindex",
      8*(-windung(0)) == 8 and all(w_ecke[tuple(1 - x if i == j else x for j, x in enumerate(e))] - w_ecke[e]
                                    == 2*(1 - 2*e[i]) for e in ecken for i in range(3))
      and all(w_ecke[tuple(e[p] for p in perm)] == w_ecke[e] for e in ecken for perm in itertools.permutations(range(3))))
befund("3", "nein (gescopt)", "E/N: Die native diskrete Struktur enthält keine Kopplung von Drehung und schwachem Spin. "
       "Die Uhr ist zentral und vertauscht mit beiden; sie ist auf K der Hodge-Stern mal i, kein Gruppenelement der "
       "verdrillten Drehungen. P_sym und die Marken wirken auf andere Faktoren. Die Wilson-Fortsetzung des Bestands ist "
       "das Gegenteil eines Kähler-Dirac-Gitters: Ihr Massenterm ist der Formengrad der Ecken, es bleibt eine Ecke "
       "ohne Geschmacksraum. Gilt für die eingesetzte Wilson-Schablone (die k-Abhängigkeit ist im Bestand nicht "
       "rekonstruiert) und für eine Naht; kein Unmöglichkeitsbeweis für ein anderes natives Gitter.")

# ============================================================================================================ Teil 4
say("\nTeil 4  Folgen")
PKL = PK @ PL                                                        # linkshändige Weyl-Felder in K: (2_L, 16)
n_weyl = int(np.trace(PKL).real)//2
nu_uhr = int((np.trace(PKL @ R)/1j).real)//2                         # Summe der Z_4-Ladungen (R = i^q)
kd_plus = int(np.trace(PKL @ (I256 + GKD)/2).real)//2
kd_minus = int(np.trace(PKL @ (I256 - GKD)/2).real)//2
ZAHLEN["4.1 Anomaliezahlen auf den linkshändigen Feldern von K"] = {"Weyl-Felder": n_weyl, "nu(Uhr)": nu_uhr,
                                                                    "Kähler-Dirac-Z_4: Ladung +1": kd_plus, "Ladung -1": kd_minus}
check(f"4.1 [E] Z_4-Anomalie: K enthält {n_weyl} linkshändige Weyl-Felder (2_L, 16), alle mit Uhrladung 1: nu(Uhr) = "
      f"{nu_uhr} steht im Bestand ohne P_KD (z154, Prüfung 2.5). Die Kähler-Dirac-Z_4 i Gamma_KD hat darauf {kd_plus} mal "
      f"+1 und {kd_minus} mal -1, also nu = 0; ihre Bedingung (vier reduzierte Felder) ist dieselbe Zahl 16 nur unter "
      f"P_KD. Mit dem Schluss von P_KD entfällt eine zweite Lesart derselben 16, keine Aussage des Bestands",
      n_weyl == 16 and nu_uhr == 16 and kd_plus == 8 and kd_minus == 8)
check("4.2 [E/L] K12, K15, K16: Die Vier als reelle Dimension der zwei schwachen Plätze (z154) verlangt die Lötung; "
      "sie ist nicht nativ und bricht die Uhr (2.12), also bleibt K12 kinematische Lesart. Für K15 hieße P_KD, dass "
      "SU(2)_L und Hyperladung keine inneren Symmetrien sind (2.13); ohne P_KD ist K15 nicht berührt. K16: Der "
      "Wilson-Halbraum des Bestands hat keinen Geschmacksraum (3.3); ein Kähler-Dirac-Kontinuumslimes wäre ein "
      "anderes Gitter als das im Bestand", zen["von N_form im Standardmodell"] == 8 and moden_wilson == 1
      and np.any(com(R, Nf)))
befund("4", "bedingt", "E/L: Für den Bestand ändert sich nichts. 'Sechzehn Weyl-Fermionen wegen der Z_4-Anomalie' ist "
       "die Uhranomalie nu = 16 und braucht P_KD nicht. K12 bleibt Lesart (Twistor-Kinematik), K15 und K16 sind nicht "
       "berührt; für K13 folgt nur, dass die Uhrspur 64 eps unter der Verdrillung gleich bleibt und die native so(3,2) "
       "sich nicht farbneutral verdrillen lässt.")

# ============================================================================================================ Teil 5
say("\nTeil 5  Bilanz")
bilanz = [
    ("E", "Aufbau", "256 = Cl(4,2) x Cl(10) mit 16 Erzeugern; S^+, K, T16, R, Gamma_2, so(3,2), Eichgeneratoren exakt "
     "wie z118/z150. Lesart B aus z154 mit Matrizen bestätigt: K = Lambda^gerade x 4 + Lambda^ungerade x 4-quer, "
     "K' die anderen Hälften. Neu: N_form = 2 - (i/2) Summe eta_a gamma_a c_a, ein gemischter Erzeuger von so(16), und "
     "Gamma_KD = -Gamma_2 chi_c Gamma_ges; der Reduktionsprojektor ist auf S^+ der Indexprojektor P_-."),
    ("E negativ", "D_perp", "Der native Indexoperator D_perp = gamma_- A_q + gamma_+ A_q^+ verschwindet auf K (K x "
     "Nullmode ist sein Kern), vertauscht mit Gamma_KD, hat Formengrad 0 und ist ein Skalar der verdrillten Gruppe. "
     "Er ist nicht d - delta, für keine Lötung."),
    ("E, tautologisch", "induzierter Operator", "gamma^a p_a: K -> K^- und die Translationen Pi_a: K <-> K' "
     "antikommutieren mit Gamma_KD, ändern den Formengrad um +-1 und sind Vektoren der verdrillten Gruppe. Dasselbe "
     "gilt für jede Lötung, auch für die in z154 ausgeschlossene innere Zelle. Die Eigenschaft folgt aus der "
     "Clifford-Algebra, nicht aus dem Operator."),
    ("E", "Operatorraum", "Auf K gibt es keinen Lorentz-Vektor, der Singulett unter so(10) oder Pati-Salam ist. Alle "
     "Vektor-Singuletts (16 bzw. 32 Multipletts auf den 256) sind uhrneutral; Uhrladung -i gibt es auf der Spinorseite "
     "nicht. Unter verdrillter Gruppe x SU(4) gäbe es 128 Vektor-Multipletts; die nativen liegen in den 32, die auch "
     "ohne Verdrillung Vektoren sind."),
    ("E negativ", "native Daten gegen die Verdrillung", "Uhr: tauscht d und delta, R N_form R^-1 = 4 - N_form, die "
     "Lötung bricht R auf R^-1 c_L. Reelle Form: lorentzsch nur so(3) reell, J bildet auf die konjugierte Verdrillung "
     "ab; euklidisch (Uhrpaar E) ist die Gruppe reell, der Träger aber nicht J-stabil. so(3,2): [Pi_a, Pi_b] ist die "
     "unverdrillte Sigma_ab, eine verdrillte so(3,2) im direkten Produkt braucht einen Farbvektor. Eichgruppe: Mit der Lötung vertauscht "
     "vom Standardmodell nur su(3)."),
    ("E", "Gewichte", "Die Euler-Struktur eps x eps erscheint in der Mischspur nur mit Gewicht Gamma_KD (= -i R^-1 "
     "c_L auf K), nicht mit den nativen Gewichten Gamma_2 und R. Die Uhrspur 64 eps ist unter der Verdrillung "
     "unverändert."),
    ("E/N negativ", "Gitterseite", "Die Uhr ist zentral und dreht nichts. Die Wilson-Fortsetzung hat w = 2 p an den "
     "Ecken: eine Nullmode bei p = 0, kein Geschmacksraum, das Gegenteil eines Kähler-Dirac-Gitters. Gescopt: "
     "eingesetzte Schablone, eine Naht."),
    ("C", "Was von P_KD bleibt", "Als Aussage über den Operator ist P_KD leer (tautologisch) oder falsch (D_perp). "
     "Mit Inhalt ist sie die Prämisse P_loet: Es gibt ein natives Datum im gemischten Teil (4, (2,2,1)) von so(16), "
     "eine Lötform Phi_am, die Lorentz x schwache Spin(4) auf die Diagonale bricht. P_loet hat Uhrladung -1, lässt "
     "vom Standardmodell nur die Farbe und verträgt sich reell nur mit dem euklidischen Uhrpaar (P_eukl)."),
    ("L", "Literatur", "Kähler 1962; Graf 1978; Becher-Joos 1982 und Banks-Dothan-Horn 1982 (Dirac-Kähler-Gleichung "
     "und Gitterfermionen); Golterman-Smit 1984, Kilcup-Sharpe 1987 (Symmetrien gestaffelter Fermionen, Verschiebung "
     "= Geschmacksmatrix); Callan-Harvey 1985, Kaplan 1992 (Nullmoden an der Wand); Witten 1988 (topologische "
     "Verdrillung); Butt-Catterall-Pradhan-Toga arXiv:2101.01026, Catterall arXiv:2209.03828, 2311.02487 (Anomalie "
     "proportional zur Euler-Zahl, reduzierte Felder, Pati-Salam). Aus dem Gedächtnis zitiert, nicht neu gelesen."),
    ("-", "Prämissen des Kandidaten", "A*, I_step, P_lift, P_IF, P_univ, P_uhr, P_ord, P_curv, P_G unverändert. P_KD "
     "kommt nicht als zehnte hinzu: In der Form von z154 ist sie geschlossen. P_loet und P_eukl sind benannt, aber "
     "keine Prämissen des Kandidaten."),
    ("O", "Offen", "Ob eine native Antwort (Zwei-Naht-Komposition, Klebeströme) einen Kanal im gemischten Teil "
     "(6,10) von so(16) mit nicht verschwindendem Erwartungswert trägt. Ob die aus Antworten rekonstruierte "
     "k-Abhängigkeit (statt der Wilson-Schablone) Ecken entartet lässt. Unter P_A3 trüge die Lötform die "
     "Quantenzahlen der E8-Ströme (10,6); nicht gerechnet."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
ZAHLEN["Bilanz"] = [{"Status": st, "Thema": was, "Text": text} for st, was, text in bilanz]
urteil = ("Nein, gescopt: P_KD ist in der Form von z154 geschlossen. (1) Der native Indexoperator D_perp = gamma_- A_q + "
          "gamma_+ A_q^+ hat nicht die Form d - delta. Er verschwindet auf dem Indexkern, vertauscht mit Gamma_KD = "
          "gamma_5 x schwache Chiralität, hat Formengrad 0 und ist unter der verdrillten Gruppe ein Skalar. Das gilt "
          "für jede Lötung von Lorentz- und inneren Richtungen. Auf K selbst gibt es überhaupt keinen Operator, der "
          "Lorentz-Vektor und Eich-Singulett ist. (2) Der induzierte vierdimensionale Operator gamma^a p_a von K nach "
          "K^- und die Translationen der so(3,2) erfüllen alle drei Bedingungen exakt, aber tautologisch: Sie gelten "
          "wörtlich für jede Lötung, auch für die in z154 ausgeschlossene innere Zelle, weil sie nur {gamma_a, g_b} = 0 "
          "benutzen. Die Kähler-Dirac-Struktur liegt ganz in der Wahl der Gruppe. (3) Was entscheidet, ist ein Datum "
          "außerhalb von Lorentz-Vektor x Eich-Singulett: die zweite Clifford-Wirkung c_a, d und delta getrennt, oder "
          "der Formengrad N_form = 2 - (i/2) Summe eta_a gamma_a c_a. N_form ist ein gemischter Erzeuger von so(16) "
          "und liegt nicht in so(4,2) + so(10), dem direkten Produkt des Bestands. Kein natives Datum liegt dort. "
          "(4) Vier native Daten sprechen gegen die verdrillte Gruppe als Drehgruppe: Die Uhr tauscht d und delta und "
          "kehrt den Formengrad um (auf K ist sie i mal Hodge-Stern); die Lötung bricht die Uhr auf R^-1 c_L. "
          "Lorentzsch ist die Verdrillung nur für so(3) reell, und J bildet sie auf die konjugierte ab; euklidisch ist "
          "der Träger nicht J-stabil. Die Translationen der so(3,2) schließen auf die unverdrillte Lorentz-Algebra, "
          "und eine verdrillte so(3,2) im direkten Produkt verlangt einen Farbvektor. Mit der Lötung vertauscht vom Standardmodell nur die "
          "Farbe. (5) Gitterseite: Die native Vierteldrehung ist zentral und dreht weder Lorentz- noch schwache "
          "Indizes. Die Wilson-Fortsetzung des Bestands hat an den Ecken w = 2 p, also den Formengrad als Masse, und "
          "behält eine Ecke ohne Geschmacksraum: das Gegenteil eines Kähler-Dirac-Gitters (gescopt auf die eingesetzte "
          "Schablone). (6) Folgen: Die Uhranomalie nu = 16 steht ohne P_KD; die Kähler-Dirac-Lesart derselben 16 "
          "entfällt. K12 bleibt kinematische Lesart, K15, K16 und K13 sind nicht berührt. Die neun Prämissen bleiben, "
          "P_KD wird keine zehnte. Was mit Inhalt bliebe, ist P_loet (eine native Lötform im gemischten Teil von "
          "so(16)), nicht unentscheidbar wie P_IF, sondern mit Folgen, die dem Bestand widersprechen (Uhrbruch, "
          "Verlust von SU(2)_L x U(1)_Y). Nächster Test: ob irgendeine native Antwort einen nicht verschwindenden "
          "Kanal (6,10) trägt; verschwindet er aus Symmetriegründen, ist auch P_loet geschlossen. Nichts ist "
          "befördert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z161_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z161.txt").write_text("\n".join(OUT) + "\n")
