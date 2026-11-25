"""Z116: Prüfung der sechzehnten und siebzehnten externen Analyse (25.11.): gemeinsame Clifford-Struktur, Chiralität und
Gravitation aus Nahtindex und Uhr, virtuelle Ladungstransfers, Petz-Rekonstruierbarkeit.

Teil 1  Cl(13,1) und Cl(3,11) auf 128 Dimensionen (Lorentz x Spin(10)): Cliffordrelationen, Realstruktur J (aus den
        Faktoren gebaut), J^2 und Majorana-Weyl; S^+ = (S_L x 16) + (S_R x 16-quer), J verbindet die Summanden.
Teil 2  Teilchen-Loch-Konjugation J_F auf fünf Moden (J_F a J_F^-1 = a^+, J_F^2 = 1, tauscht die Parität); Uhr
        R = i Gamma_11 als Zentrum von Spin(10) (z82); J_ges^2, [J_ges, Gamma_ges], J R J^-1 = R.
Teil 3  Spuridentität: Tr_(S^+-)(R gamma_ab gamma_cd) = +-64 eps_abcd (36 Paare), ohne Uhr nur die Metrik.
Teil 4  Spin(4,12) auf 256 Dimensionen (konform x Spin(10)): Realstruktur mit J^2 = 1, lichtartiger Schnitt P_+- =
        gamma_-+ gamma_+-, Lorentz- und J-verträglich, 2 x 16 komplexe Komponenten.
Teil 5  Nahtindex: A_q = S^+(1 - q S^4), q = e^-2pi, Index +1, Nullmode sqrt(1 - q^2) sum q^n |4n>; D_perp = gamma_- A +
        gamma_+ A^+, D^2 = P_- A^+A + P_+ AA^+; Tr e^(-s A^+A) - Tr e^(-s AA^+) = 1; relative Indexspur mit Uhr = eps.
Teil 6  Komplex balancierte Verklebung T (x) V_16: reelles Tensorprodukt = (T (x)_C V) + (T (x)_C V-quer), I_T (x) 1 =
        1 (x) I_F wählt einen Summanden; Hyperladungen: sum Y = 0, sum Y^3 = 0, sum Y^2 = 10/3, b_Y = 41/6.
Teil 7  Flächenregel (z57): (1/16) eps (B_p - B_k)^2 = det e_p (3/2 - e_2/2 + 3 e_4/2), also t = 1; gemischte Antwort
        d/d lambda eps F ^ F |_0 = 2 eps R ^ e ^ e; G = v_geo^-2 aus 2 c_3/4; Gewichte 1, 1 - q^2, (1 - q^2)^2.
Teil 8  Virtuelle Ladungstransfers (Analyse 16): Transfer-SU(2) der Ströme aus z99, Spin-Singulett = Isospin 1 (Pauli),
        Mittelung über die Transfer-SU(2): w_Singulett = (1 + 2 e^-4pi)/3, Austauschterm log(3/(1 + 2 e^-4pi)) P_Singulett.
Teil 9  Petz-Rekonstruierbarkeit des nativen Übergangskanals (z91): F_koh = (m+1)/(2m), F_deko = (2m+1)/(3 m^2).
Teil 10 Überlappkern det(I - Z W^+)^-1: Positivität bleibt unter konformer Reskalierung erhalten.
"""
import itertools
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
rng = np.random.default_rng(116)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


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


def chiral(G):
    M = np.eye(G[0].shape[0], dtype=complex)
    for g in G:
        M = M @ g
    for ph in (1, 1j, -1, -1j):
        if np.allclose((ph*M) @ (ph*M), np.eye(M.shape[0])):
            return ph*M
    raise ValueError


def konj(G, s):
    """C mit C g* C^-1 = s g für alle g (kleiner Faktor, SVD); None, wenn es keine gibt."""
    d = G[0].shape[0]
    rows = [np.kron(g.conj().T, np.eye(d)) - s*np.kron(np.eye(d), g) for g in G]
    _, sv, vh = np.linalg.svd(np.vstack(rows))
    if sv[-1] > 1e-9:
        return None
    C = vh[-1].conj().reshape(d, d, order="F")
    C = C/abs(np.linalg.det(C))**(1/d)
    return C


def clifford_ok(G):
    eta = [np.round((g @ g)[0, 0].real, 9) for g in G]
    d = G[0].shape[0]
    ok = all(np.allclose(G[m] @ G[n] + G[n] @ G[m], 2*(m == n)*eta[m]*np.eye(d)) for m in range(len(G)) for n in range(len(G)))
    return ok, sum(e > 0 for e in eta), sum(e < 0 for e in eta)


# ------------------------------------------------------------------ Teil 1
say("Teil 1: Cl(13,1) und Cl(3,11) auf 128 Dimensionen")
Gi = euklid(5)                                     # Cl(10,0) auf 32 = Fockraum von fünf Moden
G11 = chiral(Gi)
E4 = euklid(2)
g4 = [1j*E4[0]] + E4[1:4]                          # gamma_0^2 = -1, gamma_i^2 = +1
g5 = g4[0] @ g4[1] @ g4[2] @ g4[3]                 # g5^2 = -1
I32 = np.eye(32)
erg1 = {}
for sig, f in (("13,1", 1j), ("3,11", 1.0)):
    G = [kr(g, I32) for g in g4] + [kr(f*g5, gi) for gi in Gi]
    ok, npos, nneg = clifford_ok(G)
    Gt = chiral(G)
    beste = None
    for s4 in (1, -1):
        for s10 in (1, -1):
            C4, C10 = konj(g4, s4), konj(Gi, s10)
            if C4 is None or C10 is None:
                continue
            C = np.kron(C4, C10)
            alle = [np.allclose(C @ g.conj() @ np.linalg.inv(C), s*g) for s in (1, -1) for g in G]
            for s in (1, -1):
                if all(np.allclose(C @ g.conj() @ np.linalg.inv(C), s*g) for g in G):
                    JJ = C @ C.conj()
                    J2 = np.round(JJ[0, 0].real, 9)
                    kommut = np.allclose(C @ Gt.conj(), Gt @ C)
                    kand = (J2, kommut, s, C)
                    if beste is None or (kommut and not beste[1]) or (kommut == beste[1] and J2 > beste[0]):
                        beste = kand
    erg1[sig] = (ok, npos, nneg, beste[0], beste[1])
    if sig == "3,11":
        G311, C311, Gt311 = G, beste[3], Gt
ZAHLEN["clifford_14"] = {k: {"Clifford": v[0], "Signatur": [int(v[1]), int(v[2])], "J^2": v[3], "J vertauscht mit Chiralität": v[4]}
                         for k, v in erg1.items()}
Pp = (np.eye(128) + Gt311)/2
L4 = np.kron((np.eye(4) + 1j*g5)/2, I32)           # 4D-Chiralität (i g5 hat Eigenwerte +-1)
F16 = np.kron(np.eye(4), (np.eye(32) + G11)/2)
dimLF = np.linalg.matrix_rank(Pp @ L4 @ F16)
dimRF = np.linalg.matrix_rank(Pp @ (np.eye(128) - L4) @ (np.eye(128) - F16))
dimLB = np.linalg.matrix_rank(Pp @ L4 @ (np.eye(128) - F16))
JLF = C311 @ (Pp @ L4 @ F16).conj()
verbindet = np.linalg.matrix_rank(JLF) == 32 and np.allclose(Pp @ (np.eye(128) - L4) @ (np.eye(128) - F16) @ JLF, JLF)
check(f"Cliffordrelationen exakt für beide Signaturen ({ZAHLEN['clifford_14']}); Cl(3,11) hat eine Realstruktur mit J^2 = +1, "
      f"die mit der Gesamtchiralität vertauscht (Majorana-Weyl, 3 - 11 = -8), Cl(13,1) nur J^2 = -1. S^+ zerfällt in "
      f"(S_L x 16) mit Dimension {dimLF} und (S_R x 16-quer) mit {dimRF} (gemischt: {dimLB}); J bildet den ersten Summanden "
      f"auf den zweiten ab ({verbindet}): unabhängige Daten 2 x 16 komplex, eine Generation linkshändiger Weyl-16",
      all(v[0] for v in erg1.values()) and erg1["3,11"][3] == 1 and erg1["3,11"][4] and erg1["13,1"][3] == -1
      and dimLF == 32 and dimRF == 32 and dimLB == 0 and verbindet)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Teilchen-Loch-Konjugation, Uhr und kombinierte Realstruktur")
a5 = [(Gi[2*j] + 1j*Gi[2*j + 1])/2 for j in range(5)]          # Fermimoden aus den Majoranas der Cl(10)
car = all(np.allclose(a5[i] @ a5[j].conj().T + a5[j].conj().T @ a5[i], (i == j)*np.eye(32)) and
          np.allclose(a5[i] @ a5[j] + a5[j] @ a5[i], 0) for i in range(5) for j in range(5))
JF = np.eye(32, dtype=complex)
for j in range(5):
    JF = JF @ Gi[2*j]
JF_ok = all(np.allclose(JF @ a5[j] @ np.linalg.inv(JF), a5[j].conj().T) for j in range(5)) and np.allclose(JF @ JF, np.eye(32))
par = np.eye(32, dtype=complex)
for j in range(5):
    par = par @ (np.eye(32) - 2*a5[j].conj().T @ a5[j])
tauscht = np.allclose(JF @ par + par @ JF, 0)
gleich = np.allclose(G11, par) or np.allclose(G11, -par)
R = 1j*G11
zentrum = np.allclose(np.linalg.matrix_power(R, 4), np.eye(32)) and all(np.allclose(R @ (Gi[A] @ Gi[B]), (Gi[A] @ Gi[B]) @ R)
                                                                      for A in range(10) for B in range(10))
C4m = konj(g4, -1)
Jges = np.kron(C4m, JF)
J2ges = np.round((Jges @ Jges.conj())[0, 0].real, 9)
Gt4 = chiral(g4)
Gtot = np.kron(Gt4, par)
vert = np.allclose(Jges @ Gtot.conj(), Gtot @ Jges) or np.allclose(Jges @ Gtot.conj(), -Gtot @ Jges)
vert_plus = np.allclose(Jges @ Gtot.conj(), Gtot @ Jges)
Rges = np.kron(np.eye(4), R)
uhr = np.allclose(Jges @ Rges.conj(), Rges @ Jges)
check(f"fünf Moden aus den Majoranas von Cl(10) (CAR {car}); J_F = Produkt der fünf Majoranas a + a^+: J_F a_j J_F^-1 = a_j^+, "
      f"J_F^2 = 1 ({JF_ok}), tauscht gerade und ungerade Parität ({tauscht}), Gamma_11 = +-(-1)^F ({gleich}), also 16 = "
      f"Lambda^gerade C^5. R = i Gamma_11 hat Ordnung 4 und vertauscht mit Spin(10) ({zentrum}): das Zentrum von Spin(10), die "
      f"Uhr der 16 (z82). Kombinierte Realstruktur (Lorentz-Konjugation x J_F, antilinear): J^2 = {J2ges}, vertauscht mit der "
      f"Gesamtchiralität ({vert_plus}), J R J^-1 = R ({uhr})",
      car and JF_ok and tauscht and gleich and zentrum and vert_plus and uhr and J2ges in (1.0, -1.0))
ZAHLEN["J_ges^2"] = J2ges

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Spur mit Uhr = Levi-Civita")
G = G311
Pm = (np.eye(128) - Gt311)/2
Rg = np.kron(np.eye(4), 1j*G11)
lorentz = [(a, b) for a in range(4) for b in range(a + 1, 4)]
eps = np.zeros((4, 4, 4, 4))
for p in itertools.permutations(range(4)):
    eps[p] = round(np.linalg.det(np.eye(4)[list(p)]))
werte_p, werte_m, metr = {}, {}, {}
for (a, b) in lorentz:
    for (c, d) in lorentz:
        X = G[a] @ G[b] @ G[c] @ G[d]
        werte_p[(a, b, c, d)] = np.trace(Pp @ Rg @ X)
        werte_m[(a, b, c, d)] = np.trace(Pm @ Rg @ X)
        metr[(a, b, c, d)] = np.trace(Pp @ X)
k_p = {k: v/eps[k] for k, v in werte_p.items() if eps[k] != 0}
k_m = {k: v/eps[k] for k, v in werte_m.items() if eps[k] != 0}
konst_p = set(np.round(list(k_p.values()), 9))
konst_m = set(np.round(list(k_m.values()), 9))
null_sonst = all(abs(v) < 1e-9 for k, v in werte_p.items() if eps[k] == 0)
metr_nur = all(abs(v) < 1e-9 for k, v in metr.items() if {k[0], k[1]} != {k[2], k[3]})
ZAHLEN["spur_mit_uhr"] = {"S+ Konstante": [str(x) for x in konst_p], "S- Konstante": [str(x) for x in konst_m],
                          "ohne Uhr (a,b)=(c,d)": str(np.round(metr[(0, 1, 0, 1)], 6))}
check(f"Tr_(S^+)(R gamma_ab gamma_cd) = c eps_abcd für alle 36 Paare mit c = {konst_p} (und 0, wo eps = 0: {null_sonst}); auf S^- "
      f"c = {konst_m}; ohne Uhr nur die metrische Kontraktion (nur für {{a,b}} = {{c,d}} ungleich 0: {metr_nur}, "
      f"Tr(gamma_01 gamma_01) = {np.round(metr[(0, 1, 0, 1)].real, 6)})",
      len(konst_p) == 1 and len(konst_m) == 1 and abs(abs(list(konst_p)[0]) - 64) < 1e-9
      and abs(list(konst_p)[0] + list(konst_m)[0]) < 1e-9 and null_sonst and metr_nur)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Spin(4,12) auf 256 Dimensionen und der lichtartige Schnitt")
E6 = euklid(3)
h = [E6[0], E6[1], E6[2], E6[3], 1j*E6[4], 1j*E6[5]]    # vier positiv (drei Raum + e_4), zwei negativ (Zeit, e_-1)
h7 = chiral(h)
G16 = [kr(x, I32) for x in h] + [kr(h7, gi)*(1j if np.allclose(h7 @ h7, np.eye(8)) else 1.0) for gi in Gi]
ok16, np16, nn16 = clifford_ok(G16)
Gt16 = chiral(G16)
C16 = None
for s6 in (1, -1):
    for s10 in (1, -1):
        C8, C10 = konj(h, s6), konj(Gi, s10)
        if C8 is None or C10 is None:
            continue
        C = np.kron(C8, C10)
        for s in (1, -1):
            if all(np.allclose(C @ g.conj() @ np.linalg.inv(C), s*g) for g in G16):
                JJ = C @ C.conj()
                if np.allclose(JJ, np.eye(256)) and np.allclose(C @ Gt16.conj(), Gt16 @ C):
                    C16, s16 = C, s
gp = (G16[3] + G16[5])/2
gm = (G16[3] - G16[5])/2
null_ok = np.allclose(gp @ gp, 0) and np.allclose(gm @ gm, 0) and np.allclose(gp @ gm + gm @ gp, np.eye(256))
Pplus, Pminus = gm @ gp, gp @ gm
Sp16 = (np.eye(256) + Gt16)/2
dimP = np.linalg.matrix_rank(Pplus @ Sp16)
lor16 = [G16[a] @ G16[b] for a, b in itertools.combinations((0, 1, 2, 4), 2)]
lor_ok = all(np.allclose(L @ Pplus, Pplus @ L) for L in lor16)
J_ok = C16 is not None and np.allclose(C16 @ Pplus.conj(), Pplus @ C16) and np.allclose(C16 @ Sp16.conj(), Sp16 @ C16)
check(f"Cl(4,12) auf 256 Dimensionen: Cliffordrelationen exakt ({ok16}, Signatur ({np16}, {nn16})); Realstruktur aus "
      f"Cl(4,2)-Faktor x Cl(10)-Faktor mit J^2 = 1, verträglich mit der Gesamtchiralität ({C16 is not None}); lichtartige "
      f"gamma_+- = (h_4 +- h_-1)/2 mit gamma^2 = 0 und {{gamma_+, gamma_-}} = 1 ({null_ok}); P_+ = gamma_- gamma_+ halbiert S^+ "
      f"auf {dimP} komplexe Dimensionen, ist Lorentz-invariant ({lor_ok}) und J-verträglich ({J_ok}): mit J bleiben 64 reelle "
      "= 2 x 16 komplexe Komponenten",
      ok16 and np16 == 4 and nn16 == 12 and C16 is not None and null_ok and dimP == 64 and lor_ok and J_ok)

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: Nahtindex als Quelle des lichtartigen Schnitts")
mp.mp.dps = 30
qv = mp.exp(-2*mp.pi)
q = float(qv)
idx_ok, null_vek, spur_ok = True, [], []
for N in (13, 25, 41, 65):
    S = np.zeros((N + 1, N))
    for n in range(N):
        S[n + 1, n] = 1
    Sd = S.T[:N - 1, :N]                                            # S^+ abgeschnitten: C^N -> C^(N-1)
    S3 = np.zeros((N - 1, N))
    for n in range(N):
        if n + 3 < N - 1:
            S3[n + 3, n] = 1
    A = Sd - q*S3
    ev = np.linalg.svd(A, compute_uv=False)
    ker = N - int(np.sum(ev > 1e-14))
    cok = (N - 1) - int(np.sum(ev > 1e-14))
    idx_ok &= ker - cok == 1 and ker == 1
    _, _, vh = np.linalg.svd(A)
    v = vh[-1]*np.sign(vh[-1][0])
    soll = np.array([np.sqrt(1 - q**2)*q**(n//4) if n % 4 == 0 else 0.0 for n in range(N)])
    null_vek.append(float(np.max(np.abs(v - soll))))
    for s_ in (0.3, 1.0, 4.0):
        e1 = np.sum(np.exp(-s_*np.linalg.eigvalsh(A.T @ A)))
        e2 = np.sum(np.exp(-s_*np.linalg.eigvalsh(A @ A.T)))
        spur_ok.append(abs(e1 - e2 - 1))
D2_ok = True
Ntest = 13
S = np.zeros((Ntest, Ntest))
for n in range(Ntest - 1):
    S[n + 1, n] = 1
Afull = S.T - q*np.linalg.matrix_power(S, 3)
gpk, gmk = np.array([[0, 1], [0, 0]], complex), np.array([[0, 0], [1, 0]], complex)
Dperp = np.kron(gmk, Afull) + np.kron(gpk, Afull.T)
D2 = Dperp @ Dperp
D2_ok = np.allclose(D2, np.kron(gpk @ gmk, Afull.T @ Afull) + np.kron(gmk @ gpk, Afull @ Afull.T))
Tp = complex(np.trace(Pplus @ Sp16 @ np.kron(np.eye(8), 1j*G11) @ G16[0] @ G16[1] @ G16[2] @ G16[4]))
Tm = complex(np.trace(Pminus @ Sp16 @ np.kron(np.eye(8), 1j*G11) @ G16[0] @ G16[1] @ G16[2] @ G16[4]))
rel = {s_: Tm*np.sum(np.exp(-s_*np.linalg.eigvalsh(A.T @ A))) + Tp*np.sum(np.exp(-s_*np.linalg.eigvalsh(A @ A.T)))
       for s_ in (0.3, 1.0, 4.0)}
ZAHLEN["nahtindex"] = {"Nullvektor-Abweichung": null_vek, "max |Spurdifferenz - 1|": float(max(spur_ok)),
                       "Sektorspuren (P_-, P_+) mit Uhr, eps_0123": [str(np.round(Tm, 6)), str(np.round(Tp, 6))],
                       "relative Spur": {str(k): str(np.round(v, 6)) for k, v in rel.items()}}
gleich_rel = len(set(np.round(list(rel.values()), 6))) == 1 and abs(abs(list(rel.values())[0]) - abs(Tm)) < 1e-6
check(f"A_q = S^+(1 - q S^4), q = e^-2pi: auf rechteckigen Abschneidungen (N = 13 .. 65) Kern 1, Kokern 0 ({idx_ok}); der "
      f"Nullvektor ist sqrt(1 - q^2) sum q^n |4n> (Abweichung max {max(null_vek):.0e}); Tr e^(-s A^+A) - Tr e^(-s AA^+) = 1 für "
      f"s = 0,3, 1, 4 (Fehler {max(spur_ok):.0e}); D_perp = gamma_- A + gamma_+ A^+ gibt D^2 = P_- A^+A + P_+ AA^+ ({D2_ok}). Die "
      f"lichtartigen Sektoren tragen die Uhrspur mit entgegengesetztem Vorzeichen ({np.round(Tm, 3)} und {np.round(Tp, 3)} "
      f"bei eps_0123); die relative Spur ist für jedes s gleich ({gleich_rel}): der Index +1 macht aus der Aufhebung genau "
      "einen Levi-Civita-Beitrag",
      idx_ok and max(null_vek) < 1e-10 and max(spur_ok) < 1e-9 and D2_ok and abs(Tp + Tm) < 1e-9 and abs(Tp) > 1 and gleich_rel)

# ------------------------------------------------------------------ Teil 6
say("\nTeil 6: Komplex balancierte Verklebung und Hyperladung")
nT, nV = 4, 16
IT = np.kron(np.eye(nT), np.array([[0, -1], [1, 0]]))
IF = np.kron(np.eye(nV), np.array([[0, -1], [1, 0]]))
M = np.kron(IT, np.eye(2*nV)) - np.kron(np.eye(2*nT), IF)
Mm = np.kron(IT, np.eye(2*nV)) + np.kron(np.eye(2*nT), IF)
dim_bal = 2*nT*2*nV - np.linalg.matrix_rank(M)
dim_anti = 2*nT*2*nV - np.linalg.matrix_rank(Mm)
y = [sp.Rational(-1, 3)]*3 + [sp.Rational(1, 2)]*2
Y = [sum(y[i] for i in S_) for k in (0, 2, 4) for S_ in itertools.combinations(range(5), k)]
bY = sp.Rational(2, 3)*3*sum(v**2 for v in Y) + sp.Rational(1, 3)*2*sp.Rational(1, 4)
check(f"reelles Tensorprodukt von T (C^4) und V_16 (C^16): reelle Dimension {2*nT*2*nV} = {dim_bal} (I_T (x) 1 = 1 (x) I_F, "
      f"also T (x)_C V_16, komplex 64) + {dim_anti} (I_T (x) 1 = -1 (x) I_F, also T (x)_C V-quer): die balancierte Verklebung "
      f"entfernt genau einen "
      f"Summanden. Hyperladungen der 16 = Lambda^gerade C^5: sum Y = {sum(Y)}, sum Y^3 = {sum(v**3 for v in Y)}, sum Y^2 = "
      f"{sum(v**2 for v in Y)}; drei Familien und ein Higgs-Dublett: b_Y = {bY}",
      dim_bal == 128 and dim_anti == 128 and sum(Y) == 0 and sum(v**3 for v in Y) == 0
      and sum(v**2 for v in Y) == sp.Rational(10, 3) and bY == sp.Rational(41, 6))

# ------------------------------------------------------------------ Teil 7
say("\nTeil 7: Flächenregel, gemischte Krümmungsantwort, G")
ok_fl = True
for _ in range(64):
    ep = rng.integers(-3, 4, size=(4, 4))
    while round(np.linalg.det(ep)) == 0:
        ep = rng.integers(-3, 4, size=(4, 4))
    ek = rng.integers(-3, 4, size=(4, 4))
    epS, ekS = sp.Matrix(ep), sp.Matrix(ek)
    X = epS.inv()*ekS
    lam = sp.symbols("lam")
    cp = sp.Poly(sp.expand((sp.eye(4) + lam*X).det()), lam)
    e2, e4 = cp.coeff_monomial(lam**2), cp.coeff_monomial(lam**4)
    t_ = sp.symbols("t")
    gem = sp.expand((epS + t_*ekS).det()).coeff(t_, 2)            # gemischtes Volumen
    eps_BB = 24*epS.det() - 2*4*gem + 24*ekS.det()
    ok_fl &= sp.simplify(eps_BB/16 - epS.det()*(sp.Rational(3, 2) - e2/2 + sp.Rational(3, 2)*e4)) == 0
lam = sp.symbols("lambda")
Rs, ees = sp.symbols("R ee")
FF = Rs**2 + 2*lam*Rs*ees + lam**2*ees**2
gem_antw = sp.diff(FF, lam).subs(lam, 0)
c3 = 1/(8*sp.pi)
v = sp.symbols("v", positive=True)
Gn = sp.symbols("G", positive=True)
loes = sp.solve(sp.Eq(2*c3/4*v**2, 1/(16*sp.pi*Gn)), Gn)
qq = mp.exp(-4*mp.pi)
gew = [1, 1 - qq, (1 - qq)**2]
check(f"Flächenregel: (1/16) eps (B_p - B_k)^2 = det e_p (3/2 - e_2/2 + 3 e_4/2) für 64 ganzzahlige Vielbeinpaare exakt "
      f"({ok_fl}), bis auf das Vorzeichen das z57-Potential (-3/2) + e_2/2 - (3/2) e_4, also t = 1; gemischte Antwort d/d "
      f"lambda eps F^F bei lambda = 0 = {gem_antw} (eps R ^ e ^ e, ohne Euler- und Volumenterm); (2 c_3/4) v^2 = v^2/(16 pi) "
      f"gibt G = {loes[0]}; Kopplungsgewichte an den Nahtzustand: ganzer Nullzustand 1, erster Randplatz 1 - e^-4pi = "
      f"{mp.nstr(gew[1], 11)}, quartischer Überlapp {mp.nstr(gew[2], 11)}",
      ok_fl and gem_antw == 2*Rs*ees and loes[0] == v**-2 and mp.nstr(gew[1], 10) == "0.9999965127"
      and mp.nstr(gew[2], 10) == "0.9999930253")

# ------------------------------------------------------------------ Teil 8
say("\nTeil 8: Virtuelle Ladungstransfers zwischen zwei Nähten")


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


c = jw(4)                                                         # A_hoch, A_runter, B_hoch, B_runter
num = [x.T @ x for x in c]
NA, NB = num[0] + num[1], num[2] + num[3]
Ep = c[0].T @ c[2] + c[1].T @ c[3]
Em = Ep.T
Hiso = (NA - NB)/2
Siso = Ep @ Em + Hiso @ Hiso - Hiso                                  # Casimir der Transfer-SU(2)
Sx = (c[0].T @ c[1] + c[1].T @ c[0] + c[2].T @ c[3] + c[3].T @ c[2])/2
Sy = (-1j*(c[0].T @ c[1] - c[1].T @ c[0]) - 1j*(c[2].T @ c[3] - c[3].T @ c[2]))/2
Sz = (num[0] - num[1] + num[2] - num[3])/2
Sspin = Sx @ Sx + Sy @ Sy + Sz @ Sz
zwei = np.isclose(np.diag(NA + NB), 2)
idx2 = np.where(zwei)[0]
sub = lambda Mx: Mx[np.ix_(idx2, idx2)]
w_spin, V = np.linalg.eigh(sub(Sspin))
iso_s = sorted(np.round(np.diag(V.conj().T @ sub(Siso) @ V).real, 9))
paar = sorted(zip(np.round(w_spin, 9), np.round(np.diag(V.conj().T @ sub(Siso) @ V).real, 9)))
pauli = all((sw == 0 and si == 2) or (sw == 2 and si == 0) for sw, si in paar)
Eld = (NA @ NA + NB @ NB)/2
dE = np.diag(sub(Eld)) - 1
beta = 4*np.pi
Dg = np.diag(np.exp(-beta*dE))
haar_s, haar_t = [], []
singulett = [i for i, (sw, si) in enumerate(paar) if sw == 0]
P_s = V[:, np.isclose(w_spin, 0)]
P_t = V[:, np.isclose(w_spin, 2)]
ws_exakt = np.trace(P_s.conj().T @ Dg @ P_s).real/P_s.shape[1]
for _ in range(400):
    Zm = (rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3)))/np.sqrt(2)
    Q, Rr = np.linalg.qr(Zm)
    Q = Q @ np.diag(np.diag(Rr)/abs(np.diag(Rr)))
    U = P_s @ Q @ P_s.conj().T
    psi = P_s[:, [list(np.round(np.diag(P_s.conj().T @ sub(Hiso) @ P_s).real, 9)).index(0.0)]]
    haar_s.append((psi.conj().T @ U.conj().T @ Dg @ U @ psi)[0, 0].real)
ws = (1 + 2*np.exp(-4*np.pi))/3
kopplung = np.log(3/(1 + 2*np.exp(-4*np.pi)))
korr = (1 - ws)/(3 + ws)
ZAHLEN["austausch"] = {"w_Singulett exakt": ws_exakt, "Haar-Mittel (400)": float(np.mean(haar_s)),
                       "log(3/(1+2e^-4pi))": kopplung, "durch 4 pi": kopplung/(4*np.pi), "Korrelation": korr}
check(f"zwei Nähte mit je zwei Moden, Sektor zweier Teilchen: die Transferströme E^+- = sum_s a_(A s)^+ a_(B s) bilden eine "
      f"SU(2) (z99); Spin und Transfer-Isospin sind gekoppelt (Pauli): Spin-Singulett = Isospin 1, Spin-Triplett = Isospin 0 "
      f"({pauli}). Ladungsenergie Q^2/2: Doppelbesetzung kostet 1. Mittelung über die Transfer-SU(2) (Schur: Spur/Dimension) "
      f"gibt w_Singulett = {ws_exakt:.10f} = (1 + 2 e^-4pi)/3 (Haar-Stichprobe {np.mean(haar_s):.4f}), w_Triplett = 1; "
      f"Austauschterm log(3/(1 + 2 e^-4pi)) P_Singulett = {kopplung:.7f}, durch 4 pi {kopplung/(4*np.pi):.8f}; lokale Zustände "
      f"maximal gemischt, <s_A^i s_B^j> = (1 - w)/(3 + w) delta = {korr:.6f} delta",
      pauli and abs(ws_exakt - ws) < 1e-15 and abs(np.mean(haar_s) - ws) < 0.05 and abs(kopplung - 1.0986053) < 1e-7
      and abs(kopplung/(4*np.pi) - 0.08742423) < 1e-8 and abs(korr - 0.2) < 1e-5)
befund("C: Die Austauschkopplung ist im Kern log 3 = log(dim Isospin 1): Sie kommt aus der Mittelung über die native "
       "Transfer-SU(2); die Ladungsenergie trägt nur die Korrektur e^-4pi (ohne Ladungsenergie wäre w = 1 und keine "
       "Wechselwirkung, bei unendlicher w = 1/3). Prämisse: die relative Transferorientierung zweier Nähte ist nicht "
       "beobachtbar (Mittelung). Eine 4D-Eichdynamik folgt daraus nicht", True)

# ------------------------------------------------------------------ Teil 9
say("\nTeil 9: Petz-Rekonstruierbarkeit des nativen Übergangskanals")


def kanal_kraus(m, koh=True):
    K1 = np.zeros((m + 1, m))
    K2 = np.zeros((m + 1, m))
    for k in range(m):
        K1[k + 1, k] = np.sqrt((k + 1)/(m + 1))
        K2[k, k] = np.sqrt((m - k)/(m + 1))
    if koh:
        return [K1, K2]
    out = []
    for k in range(m):
        A1 = np.zeros((m + 1, m))
        A1[k + 1, k] = np.sqrt((k + 1)/(m + 1))
        A2 = np.zeros((m + 1, m))
        A2[k, k] = np.sqrt((m - k)/(m + 1))
        out += [A1, A2]
    return out


def msqrt(Mx, inv=False):
    w, U = np.linalg.eigh(Mx)
    w = np.clip(w, 0, None)
    d = np.array([1/np.sqrt(x) if (inv and x > 1e-14) else (np.sqrt(x) if not inv else 0.0) for x in w])
    return U @ np.diag(d) @ U.conj().T


def petz_fidelitaet(Ks, m):
    sigma = np.eye(m)/m
    Phis = sum(K @ sigma @ K.T for K in Ks)
    Pinv = msqrt(Phis, inv=True)
    Rk = [msqrt(sigma) @ K.T @ Pinv for K in Ks]
    psi = np.eye(m).reshape(-1)/np.sqrt(m)
    rho = np.outer(psi, psi)
    out = np.zeros((m*m, m*m))
    for Kx in Ks:
        for Rx in Rk:
            Op = np.kron(np.eye(m), Rx @ Kx)
            out += Op @ rho @ Op.T
    return float(psi @ out @ psi)


fk, fd = {}, {}
for m in (2, 3, 5, 8, 12):
    fk[m] = petz_fidelitaet(kanal_kraus(m, True), m)
    fd[m] = petz_fidelitaet(kanal_kraus(m, False), m)
ZAHLEN["petz"] = {"kohärent": fk, "dekohärent": fd}
check(f"Petz-Rückkanal (Referenz maximal gemischt), Verschränkungsfidelität: kohärenter Kanal (Krausrang 2, z91) "
      f"{ {m: round(v, 6) for m, v in fk.items()} } = (m+1)/(2m); dekohärenter Kanal mit denselben Übergängen "
      f"{ {m: round(v, 6) for m, v in fd.items()} } = (2m+1)/(3 m^2)",
      all(abs(fk[m] - (m + 1)/(2*m)) < 1e-10 for m in fk) and all(abs(fd[m] - (2*m + 1)/(3*m**2)) < 1e-10 for m in fd))
befund("C: Die Zahlen der Analyse 16 stimmen; der kohärente Kanal erhält mehr rekonstruierbare Information als der "
       "dekohärente. Als Prinzip ('der Übergang vernichtet keine Information, die die Raten nicht erzwingen') könnte das "
       "P_koh ersetzen; dass die Quelle dieses Extremalprinzip erfüllt, ist offen", True)

# ------------------------------------------------------------------ Teil 10
say("\nTeil 10: Überlappkern und konforme Reskalierung")
pts = []
while len(pts) < 8:
    Zc = (rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2)))*0.4
    if np.max(np.linalg.svd(Zc, compute_uv=False)) < 0.95:
        pts.append(Zc)
Kmat = np.array([[1/np.linalg.det(np.eye(2) - Zi @ Zj.conj().T) for Zj in pts] for Zi in pts])
Om = np.diag(rng.uniform(-2, 2, size=8))
ev0 = np.linalg.eigvalsh((Kmat + Kmat.conj().T)/2).min()
ev1 = np.linalg.eigvalsh(Om @ ((Kmat + Kmat.conj().T)/2) @ Om).min()
check(f"Überlappkern det(I - Z W^+)^-1 (nu = 1) auf acht Punkten der Einheitskugel: hermitesch und positiv (kleinster Eigenwert "
      f"{ev0:.2e}); nach Reskalierung mit reellen Omega(Z) bleibt er positiv ({ev1:.2e}): Die Positivität legt das Vorzeichen "
      "der Krümmung nicht fest; radiale AdS-Referenz und physikalisches Lambda sind zu trennen (vgl. z98, z104)",
      ev0 > -1e-12 and ev1 > -1e-12 and np.allclose(Kmat, Kmat.conj().T))

urteil = ("Die sechzehnte und die siebzehnte externe Analyse halten in allen gerechneten Punkten. (1) Cl(3,11) auf 128 "
          "Dimensionen hat eine Realstruktur mit J^2 = +1, die mit der Gesamtchiralität vertauscht (Majorana-Weyl), Cl(13,1) "
          "nur J^2 = -1; S^+ = (S_L x 16) + (S_R x 16-quer), J verbindet die Summanden, also eine Generation linkshändiger "
          "Weyl-16. (2) Die Teilchen-Loch-Konjugation J_F (Produkt der fünf Majoranas) tauscht die Parität; R = i Gamma_11 ist "
          "das Zentrum von Spin(10), also die Uhr der 16 (z82), und J R J^-1 = R. (3) Tr_(S^+-)(R gamma_ab gamma_cd) = "
          "+-64 eps_abcd, ohne Uhr nur die Metrik. (4) Cl(4,12) auf 256 Dimensionen mit J^2 = 1, lichtartiger Schnitt "
          "Lorentz- und J-verträglich, 2 x 16 komplexe Komponenten. (5) Der Nahtindex A_q = S^+(1 - q S^4) (Index +1, "
          "Nullmode sqrt(1 - q^2) sum q^n |4n>) gibt mit D_perp = gamma_- A + gamma_+ A^+ eine relative Spur, die die "
          "entgegengesetzten Uhrspuren der lichtartigen Sektoren zu genau einem Levi-Civita-Beitrag macht. (6) Die komplex "
          "balancierte Verklebung T (x)_C V_16 entfernt einen Summanden; Hyperladungen und b_Y = 41/6 unverändert. (7) Die "
          "Flächenregel ist die aus z57 (t = 1); die gemischte Antwort ist 2 eps R e e ohne Euler- und Volumenterm; G = "
          "v_geo^-2 folgt aus dem Gewicht 2 c_3, das offen bleibt. (8) Die Austauschkopplung log(3/(1 + 2 e^-4pi)) folgt aus "
          "der Mittelung über die native Transfer-SU(2) (Spin-Singulett = Isospin 1); sie ist im Kern log 3. (9) Petz: "
          "(m+1)/(2m) gegen (2m+1)/(3 m^2). (10) Die Positivität des Überlappkerns legt die Krümmung nicht fest. Offen: dass "
          "der H6-Tomita-Operator die gemeinsame Realstruktur erzeugt, dass die native Quelle D_perp und die Spur mit Uhr als "
          "Wirkung liefert, die Verklebung I_T = I_F, das Gewicht 2 c_3 und die Mittelungsprämisse.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z116_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
