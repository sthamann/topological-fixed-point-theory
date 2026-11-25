"""Z83: Radiale Quantisierung (Strategie S2): Sind die Moden m Energiestufen eines konformen 4D-Felds statt Spezies?

Leitlinie 'Rekonstruktion statt Auswahl'. Firewalls: GRAV.NONCIRCULAR.01 (die Minkowskizeit ist der parabolische Erzeuger
aus den Stufenhebern, nicht der Modulfluss; der Modulfluss bleibt Helizitätsphase) und DIMENSION.UPLIFT.FIREWALL.01 (alles
ist Kinematik der Einpaarschicht; keine wechselwirkende 3+1-Schließung). Bestand: z72 (parabolischer Erzeuger, Laguerre),
z73 (Paarschalen H_m = span{psi*_(r+m) psi_r Omega}, dim m; Toeplitz-Hürde; Summe m^2 q^m), z74 (Kern nu = 1), z76/z80
(zustandserhaltender Zweipfadkanal, Petz l(l+1), Mehrpaaranteil (m-2)/(3m-2)), z78/z79 (Araki-Woods, Twistorraum je Mode),
z75 (Speziesproblem), z81 (R = i^(Q_a - Q_b)).

Bezeichnungen: Paarzustand |m, k> = psi*_p psi_r Omega mit Teilchenlage n_1 = p = k und Lochtiefe n_2 = -r-1 = m-1-k;
L_0 = x_p - x_r = m in beiden Quartersektoren. Zwei native Elementarschritte H_m -> H_(m+1): S_1 'Teilchen hoch'
(k -> k+1) und S_2 'Loch runter' (k -> k). Tomita-Partner H-quer_m mit den J-Bildern der Schritte. Verdoppelter Raum
(+)_(m <= M) H_m (x) H-quer_m (gleiche Stufe: Standardform der KMS-Kompression je Schale).

Teil 1  Heber P_(a adot) = S_a (x) S-quer_adot, Senker K = P^+, Drehungen aus den Bilinearen S_c S_a^+, E = L_0 = m, für drei
        Normierungen: (i) nackte Hardy-Hälften von J_-1 (+1, -1; z73), (ii) Kraus-Amplituden des zustandserhaltenden
        Kanals (z76), (iii) Schwinger/CCR sqrt(n_1 + 1), sqrt(n_2 + 1). [E, P] = P, [P, P] = 0, [K, P] im Spann von E und
        so(4), erzeugte Lie-Algebra, Strukturkonstanten gegen gl(4), Jacobi, Killing-Signatur.
Teil 2  Native Suche: Satz über Schalenfaktoren (P = F(E) P^osc, K = G(E) K^osc schließt genau für F G = const);
        Kandidaten (Kanal mal sqrt(m+1), modulare Gewichtung, Petz als K, geometrisches Mittel, thermische Abbildung); der
        native Möbiusheber L^nat_-1 = L_-1 - (alpha - 1/2) J_-1 = Summe (k+1) psi*_(k+1) psi_k im Fockraum (paarzahlerhaltend,
        Amplituden n_a + 1, sl(2,R) auf der Einpaarschicht); Kanal = Möbius/L_0, CCR = Tomita-Wurzel; z76-Petz als Casimir
        der radialen sl(2,R).
Teil 3  Darstellung: K-Typen (j, j) bei E = 2j+1, Charakter, Kohärenzkern det(I - Z^+ W)^-1, Casimir (Skalar; volle
        Standardform mit jeder Helizität einmal), Minkowskienergie E + (P^0 + K^0)/2 (parabolisch, positiv, Spektrum [0, inf)
        mit Gamma-Maßen).
Teil 4  Thermische Gewichte: einseitige Kompression, Produkt mit dem Partner, Purifikation (Kohärenzzustand), Modulfluss.
Teil 5  Uhr R = e^{2 pi i L_0} auf Schritten und Schalen; Spiegel x -> -x zwischen den Quartersektoren.
Teil 6  Bilanz (E/C/L/O) und Urteil gegen Erfolgs- und Abbruchkriterien von S2.
"""
import itertools
import json
import time
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla
import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
rng = np.random.default_rng(83)
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


# ------------------------------------------------------------------ Schalenleiter, Schritte, Verdopplung
AMP = {
    "nackt": (lambda m, k: 1.0, lambda m, k: -1.0),                                        # Hardy-Hälften von J_-1 (z73, z80)
    "Kanal": (lambda m, k: np.sqrt((k + 1)/(m + 1)), lambda m, k: np.sqrt((m - k)/(m + 1))),  # z76
    "CCR": (lambda m, k: np.sqrt(k + 1), lambda m, k: np.sqrt(m - k)),                     # Schwinger: sqrt(n_1+1), sqrt(n_2+1)
}
SCHRITT = (lambda k: k + 1, lambda k: k)                                                   # S_1 Teilchen hoch, S_2 Loch runter


def verdoppelt(M, amp, voll=False):
    """Operatoren auf (+) H_m (x) H-quer_mb (voll=False: nur mb = m). Labels (m, mb, k, kb).
    P[a,b] = S_a (x) S-quer_b, A[c,a] = (S_c S_a^+) (x) 1, B[d,b] = 1 (x) (S-quer_d S-quer_b^+), E = (m + mb)/2."""
    if voll:
        lab = [(m, mb, k, kb) for m in range(1, M + 1) for mb in range(1, M + 1) for k in range(m) for kb in range(mb)]
    else:
        lab = [(m, m, k, kb) for m in range(1, M + 1) for k in range(m) for kb in range(m)]
    ix = {s: i for i, s in enumerate(lab)}
    n = len(lab)
    P = {}
    for a in range(2):
        for b in range(2):
            r, c, v = [], [], []
            for (m, mb, k, kb), i in ix.items():
                if m < M and mb < M:
                    r.append(ix[(m + 1, mb + 1, SCHRITT[a](k), SCHRITT[b](kb))])
                    c.append(i)
                    v.append(amp[a](m, k)*amp[b](mb, kb))
            P[(a, b)] = sps.csr_matrix((v, (r, c)), shape=(n, n))

    def bilinear(c_, a_, seite):
        r, cc, v = [], [], []
        for (m, mb, k, kb), i in ix.items():
            mm, kk = (m, k) if seite == 0 else (mb, kb)
            k1 = kk - (1 if a_ == 0 else 0)                          # S_a^+: zurück nach H_(mm-1)
            if mm == 1 or not 0 <= k1 <= mm - 2:
                continue
            k2 = SCHRITT[c_](k1)
            ziel = (m, mb, k2, kb) if seite == 0 else (m, mb, k, k2)
            r.append(ix[ziel])
            cc.append(i)
            v.append(amp[a_](mm - 1, k1)*amp[c_](mm - 1, k1))
        return sps.csr_matrix((v, (r, cc)), shape=(n, n))

    A = {(c_, a_): bilinear(c_, a_, 0) for c_ in range(2) for a_ in range(2)}
    B = {(c_, a_): bilinear(c_, a_, 1) for c_ in range(2) for a_ in range(2)}
    E = sps.diags([0.5*(s[0] + s[1]) for s in lab]).tocsr()
    return lab, P, A, B, E


def dicht(d):
    return {k: v.toarray() for k, v in d.items()}


def lie_dims(P, K, spalten):
    """Dimension des Spanns von {P, K}, + [K, P], + [[K, P], P oder K], + [[K, P], [K, P]] (auf Spalten weg vom Rand).
    Jeder Vektor wird normiert (Rang skaleninvariant); ein Kommutator gilt als null, wenn seine Norm unter 1e-11 mal dem
    Produkt der Normen seiner Faktoren liegt (Rundung, nicht Struktur)."""
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


def dims_fuer(amp, M):
    lab, P, _, _, _ = verdoppelt(M, amp)
    P = dicht(P)
    return lie_dims(P, {k: v.T for k, v in P.items()}, np.array([s[0] for s in lab]) <= M - 3)


def E_ext(P, K):
    return sum(K[x] @ P[x] - P[x] @ K[x] for x in P)/4


def schalenwerte(X, mv, ms):
    """Diagonalwert je Schale und Abweichung von einem Vielfachen der Eins."""
    out = {}
    for m in ms:
        sel = np.where(mv == m)[0]
        blk = X[np.ix_(sel, sel)]
        out[m] = (float(blk[0, 0]), float(np.abs(blk - blk[0, 0]*np.eye(len(sel))).max()))
    return out


# =================================================================== Teil 1
say("Teil 1: Heber und Senker auf der verdoppelten Schalenleiter, drei Normierungen")
M1 = 9
T1 = {}
for name, amp in AMP.items():
    lab, P, A, B, E = verdoppelt(M1, amp)
    P, A, B, E = dicht(P), dicht(A), dicht(B), E.toarray()
    K = {k: v.T for k, v in P.items()}
    mv = np.array([s[0] for s in lab])
    innen = mv <= M1 - 2
    eP = max(np.abs((E @ X - X @ E - X)[:, innen]).max() for X in P.values())
    pp = max(np.abs((X @ Y - Y @ X)[:, innen]).max() for X in P.values() for Y in P.values())
    T1[name] = dict(lab=lab, P=P, K=K, A=A, B=B, E=E, mv=mv, eP=eP, pp=pp)
check("für alle drei Normierungen gilt [E, P_(a adot)] = P_(a adot) mit E = L_0 = m und [P, P] = 0 (Pfadunabhängigkeit: erst "
      "Loch, dann Teilchen = erst Teilchen, dann Loch; M = 9, weg vom Rand): Graduierung und Vertauschbarkeit der Heber sind "
      "normierungsunabhängig", all(v["eP"] < 1e-12 and v["pp"] < 1e-12 for v in T1.values()))

# --- (i) nackt
t = T1["nackt"]
Ee = E_ext(t["P"], t["K"])
sw = schalenwerte(Ee, t["mv"], range(1, M1 - 1))
Q0 = np.diag([1.0 if s[2] == 0 else 0.0 for s in t["lab"]])                  # Projektor Teilchenlage k = 0 (links)
Q0b = np.diag([1.0 if s[3] == 0 else 0.0 for s in t["lab"]])
rand = t["K"][(0, 0)] @ t["P"][(0, 0)] - t["P"][(0, 0)] @ t["K"][(0, 0)] - (Q0 + Q0b - Q0 @ Q0b)
dims_n = {M: dims_fuer(AMP["nackt"], M) for M in (8, 10)}
ZAHLEN["nackt"] = {"E_ext_je_Schale(Wert, Abweichung von Skalar)": sw, "erzeugte_dimensionen_M8_M10": dims_n}
say(f"   nackt: E_ext je Schale {{m: (Diagonalwert, Abweichung)}} = {sw}; erzeugte Dimensionen {dims_n}")
check("[negativ, E] (i) nackte Schritte (Toeplitz, Amplituden +1/-1): [K_11, P_11] = Q_0 (x) 1 + 1 (x) Q_0 - Q_0 (x) Q_0 "
      "(Randprojektoren der Teilchenlage 0, wie z73), E_ext = (1/4) Summe [K, P] ist ab m = 3 kein Vielfaches der Eins, und "
      "die von P, K erzeugte Lie-Algebra wächst 8 -> 20 -> 68 -> 94 (gleich für M = 8 und 10, also kein Abschneideeffekt): "
      "kein su(2,2)", np.abs(rand[:, t["mv"] <= M1 - 1]).max() < 1e-12
      and all(sw[m][1] > 0.1 for m in range(3, M1 - 1)) and dims_n[8] == dims_n[10] == [8, 20, 68, 94])

# --- (ii) Kanal
t = T1["Kanal"]
Ee = E_ext(t["P"], t["K"])
sw = schalenwerte(Ee, t["mv"], range(1, M1 - 1))
soll = {m: (2*m - 1)/(4*m**2) for m in sw}
zweite_diff = [sw[m + 1][0] - 2*sw[m][0] + sw[m - 1][0] for m in range(2, M1 - 2)]
dims_k = {M: dims_fuer(AMP["Kanal"], M) for M in (8, 10)}
ZAHLEN["Kanal"] = {"E_ext_je_Schale": {m: v[0] for m, v in sw.items()}, "soll_(2m-1)/(4m^2)": soll,
                   "erzeugte_dimensionen_M8_M10": dims_k}
say(f"   Kanal: E_ext je Schale {{m: Wert}} = { {m: round(v[0], 6) for m, v in sw.items()} }; erzeugte Dimensionen {dims_k}")
check("[negativ, E] (ii) Kanalamplituden (z76): E_ext ist auf jeder Schale ein Vielfaches der Eins (das k-Profil stimmt), aber "
      "E_ext = (2m - 1)/(4m^2) ist nicht affin in m (zweite Differenzen != 0), also [E_ext, P] != lambda^2 P für jedes lambda; "
      "die erzeugte Lie-Algebra wächst 8 -> 24 -> 104 -> 155 (M = 8 und 10): kein su(2,2)",
      all(sw[m][1] < 1e-12 and abs(sw[m][0] - soll[m]) < 1e-12 for m in sw) and min(abs(z) for z in zweite_diff) > 1e-3
      and dims_k[8] == dims_k[10] == [8, 24, 104, 155])

# --- (iii) CCR: exakt (sympy, M = 5) und numerisch (M = 9)
MS = 5
labs = [(m, k, kb) for m in range(1, MS + 1) for k in range(m) for kb in range(m)]
ixs = {s: i for i, s in enumerate(labs)}
ns = len(labs)
hs = [lambda m, k: sp.sqrt(k + 1), lambda m, k: sp.sqrt(m - k)]


def P_sym(a, b):
    d = {}
    for (m, k, kb), i in ixs.items():
        if m < MS:
            d[(ixs[(m + 1, SCHRITT[a](k), SCHRITT[b](kb))], i)] = hs[a](m, k)*hs[b](m, kb)
    return sp.SparseMatrix(ns, ns, d)


def bil_sym(c_, a_, seite):
    d = {}
    for (m, k, kb), i in ixs.items():
        kk = k if seite == 0 else kb
        k1 = kk - (1 if a_ == 0 else 0)
        if m == 1 or not 0 <= k1 <= m - 2:
            continue
        k2 = SCHRITT[c_](k1)
        ziel = (m, k2, kb) if seite == 0 else (m, k, k2)
        d[(ixs[ziel], i)] = hs[a_](m - 1, k1)*hs[c_](m - 1, k1)
    return sp.SparseMatrix(ns, ns, d)


Ps = {(a, b): P_sym(a, b) for a in range(2) for b in range(2)}
Ks = {k: v.T for k, v in Ps.items()}
As = {(c_, a_): bil_sym(c_, a_, 0) for c_ in range(2) for a_ in range(2)}
Bs = {(c_, a_): bil_sym(c_, a_, 1) for c_ in range(2) for a_ in range(2)}
Es = sp.SparseMatrix(ns, ns, {(i, i): s[0] for i, s in enumerate(labs)})
A0s = {k: (As[k] - (As[(0, 0)] + As[(1, 1)])/2 if k[0] == k[1] else As[k]) for k in As}
B0s = {k: (Bs[k] - (Bs[(0, 0)] + Bs[(1, 1)])/2 if k[0] == k[1] else Bs[k]) for k in Bs}
innen_s = [i for i, s in enumerate(labs) if s[0] <= MS - 1]
ok_exakt = True
for a, b, c_, d_ in itertools.product(range(2), repeat=4):
    C = Ks[(a, b)]*Ps[(c_, d_)] - Ps[(c_, d_)]*Ks[(a, b)]
    T = sp.SparseMatrix(ns, ns, {})
    if a == c_ and b == d_:
        T = T + Es
    if a == c_:
        T = T + B0s[(d_, b)]
    if b == d_:
        T = T + A0s[(c_, a)]
    ok_exakt &= all(sp.simplify(x) == 0 for x in (C - T)[:, innen_s].values())
t = T1["CCR"]
A0 = {k: (t["A"][k] - (t["A"][(0, 0)] + t["A"][(1, 1)])/2 if k[0] == k[1] else t["A"][k]) for k in t["A"]}
B0 = {k: (t["B"][k] - (t["B"][(0, 0)] + t["B"][(1, 1)])/2 if k[0] == k[1] else t["B"][k]) for k in t["B"]}
res_kp = 0.0
for a, b, c_, d_ in itertools.product(range(2), repeat=4):
    C = t["K"][(a, b)] @ t["P"][(c_, d_)] - t["P"][(c_, d_)] @ t["K"][(a, b)]
    T = (t["E"] if (a == c_ and b == d_) else 0) + (B0[(d_, b)] if a == c_ else 0) + (A0[(c_, a)] if b == d_ else 0)
    res_kp = max(res_kp, np.abs((C - T)[:, t["mv"] <= M1 - 1]).max())
dims_c = {M: dims_fuer(AMP["CCR"], M) for M in (8, 10)}
ZAHLEN["CCR"] = {"residuum_[K,P]": res_kp, "erzeugte_dimensionen_M8_M10": dims_c}
check("[E] (iii) CCR-Amplituden sqrt(n_1+1), sqrt(n_2+1): [K_ab, P_cd] = delta_ac delta_bd E + delta_ac B0_db + delta_bd A0_ca "
      "mit E = L_0 = m und den spurfreien Bilinearen A0, B0 der Schritte, exakt (sympy, alle 16 Paare, M = 5, weg vom Rand) "
      f"und numerisch bei M = 9 (Residuum {res_kp:.1e}); die erzeugte Lie-Algebra hat Dimension 8 -> 15 -> 15 -> 15",
      ok_exakt and res_kp < 1e-12 and dims_c[8] == dims_c[10] == [8, 15, 15, 15])

# --- Strukturkonstanten, gl(4)-Bild, Jacobi, Killing
NAMEN = ["E", "A12", "A21", "HA", "B12", "B21", "HB", "P11", "P12", "P21", "P22", "K11", "K12", "K21", "K22"]


def basis15(t):
    A, B, P, K = t["A"], t["B"], t["P"], t["K"]
    return [t["E"], A[(0, 1)], A[(1, 0)], (A[(0, 0)] - A[(1, 1)])/2, B[(0, 1)], B[(1, 0)], (B[(0, 0)] - B[(1, 1)])/2,
            P[(0, 0)], P[(0, 1)], P[(1, 0)], P[(1, 1)], K[(0, 0)], K[(0, 1)], K[(1, 0)], K[(1, 1)]]


def strukturkonstanten(X, spalten):
    V = np.array([x[:, spalten].ravel() for x in X]).T
    f = np.zeros((len(X), len(X), len(X)), dtype=complex)
    res = 0.0
    for i in range(len(X)):
        for j in range(len(X)):
            C = (X[i] @ X[j] - X[j] @ X[i])[:, spalten].ravel()
            c, *_ = np.linalg.lstsq(V, C, rcond=None)
            res = max(res, np.abs(V @ c - C).max())
            f[i, j] = c
    return f, res


X15 = basis15(t)
f15, res15 = strukturkonstanten(X15, t["mv"] <= M1 - 2)
ETA = np.diag([1.0, 1.0, -1.0, -1.0])


def Eij(i, j):
    Z = np.zeros((4, 4))
    Z[i, j] = 1
    return Z


# Oszillatorbild psi = (a_1, a_2, b_1^+, b_2^+): Q(M) = psi^+ M psi -> M eta ist ein Lie-Homomorphismus nach gl(4)
phiA = lambda c_, a_: Eij(c_, a_) @ ETA
phiB = lambda d_, b_: Eij(2 + b_, 2 + d_) @ ETA
phiP = lambda a_, b_: Eij(a_, 2 + b_) @ ETA
phiK = lambda a_, b_: Eij(2 + b_, a_) @ ETA
Y15 = [ETA/2, phiA(0, 1), phiA(1, 0), (phiA(0, 0) - phiA(1, 1))/2, phiB(0, 1), phiB(1, 0), (phiB(0, 0) - phiB(1, 1))/2,
       phiP(0, 0), phiP(0, 1), phiP(1, 0), phiP(1, 1), phiK(0, 0), phiK(0, 1), phiK(1, 0), phiK(1, 1)]
g15, resY = strukturkonstanten(Y15, np.ones(4, dtype=bool))
jac = 0.0
for i, j, k in itertools.combinations(range(15), 3):
    s_ = f15[i, j] @ f15[:, k, :] + f15[j, k] @ f15[:, i, :] + f15[k, i] @ f15[:, j, :]
    jac = max(jac, np.abs(s_).max())
AD = [f15[i].T for i in range(15)]                                          # ad_i[k, j] = f_ij^k
KILL = np.array([[np.trace(AD[i] @ AD[j]) for j in range(15)] for i in range(15)])


def vek(koeff):
    v = np.zeros(15, dtype=complex)
    for nm, c in koeff.items():
        v[NAMEN.index(nm)] += c
    return v


HERM = [vek({"E": 1}), vek({"A12": 1, "A21": 1}), vek({"A12": -1j, "A21": 1j}), vek({"HA": 1}),
        vek({"B12": 1, "B21": 1}), vek({"B12": -1j, "B21": 1j}), vek({"HB": 1})]
for p_ in ("11", "12", "21", "22"):
    HERM += [vek({"P" + p_: 1, "K" + p_: 1}), vek({"P" + p_: -1j, "K" + p_: 1j})]
REELL = [1j*h for h in HERM]                                                 # reelle Form: i mal hermitesch
BR = np.array([[u @ KILL @ w for w in REELL] for u in REELL])
ewB = np.linalg.eigvalsh(BR.real)
sig = (int(np.sum(ewB > 1e-9)), int(np.sum(ewB < -1e-9)))
Ee = E_ext(t["P"], t["K"])
sw = schalenwerte(Ee, t["mv"], range(1, M1 - 1))
ZAHLEN["su22"] = {"abschluss_residuum": res15, "strukturkonstanten_gegen_gl4": float(np.abs(f15 - g15).max()),
                  "jacobi": jac, "killing_rang": int(np.linalg.matrix_rank(KILL)), "signatur(+,-)": sig,
                  "killing_diagonal_kompakt": np.round(np.diag(BR.real)[:7], 6).tolist(),
                  "killing_diagonal_nichtkompakt": np.round(np.diag(BR.real)[7:], 6).tolist(),
                  "E_ext_je_Schale": {m: v[0] for m, v in sw.items()}}
check(f"[E] CCR: die 15 Erzeuger (E, A0, B0, P, K) schließen (Residuum {res15:.1e}); ihre Strukturkonstanten sind exakt die "
      f"der gl(4)-Bilder psi^+ M psi -> M eta (Abweichung {np.abs(f15 - g15).max():.1e}), Jacobi {jac:.1e}, Killing-Form "
      f"nicht ausgeartet, Signatur der reellen Form {sig} (7 kompakt: E und so(4); 8 nichtkompakt: P + K, i(P - K)), also "
      "su(2,2) = so(4,2); E_ext = (1/4) Summe [K, P] = m = L_0 mit Faktor lambda = 1",
      res15 < 1e-12 and np.abs(f15 - g15).max() < 1e-12 and jac < 1e-12 and np.linalg.matrix_rank(KILL) == 15
      and sig == (8, 7) and all(abs(v[0] - m) < 1e-12 and v[1] < 1e-12 for m, v in sw.items()))
befund("Teil 1: Von den drei nativen Normierungen schließt nur die CCR-Normierung, mit Faktor 1 (E = L_0). Relativ zu den "
       "nackten Hardy-Hälften ist das der Schwinger-Faktor sqrt(n_a + 1) je Schritt, relativ zum zustandserhaltenden Kanal "
       "der Schalenfaktor sqrt(m + 1) = sqrt(L_0 der Zielschale) je Faktor, verdoppelt (m + 1). Das k-Profil des Kanals ist "
       "schon richtig; was ihm fehlt, ist allein der Schalenfaktor", True)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

# =================================================================== Teil 2
say("\nTeil 2: Native Suche nach der schließenden Normierung")
t = T1["CCR"]
POSC, KOSC, MV = t["P"], t["K"], t["mv"]
INNEN = MV <= M1 - 2


def mit_faktoren(P0_, K0_, F, G):
    """P|m> = F(m) P0_|m>, K|m+1> = G(m) K0_|m+1> (Schalenfunktionen)."""
    Fd = np.diag([F(m) for m in MV])
    Gd = np.diag([G(m) for m in MV])
    return {k: v @ Fd for k, v in P0_.items()}, {k: Gd @ v for k, v in K0_.items()}


# --- 2.1 Satz: P = F(E) P^osc, K = G(E) K^osc schließt genau für F G = const
Fz = dict(zip(range(1, M1 + 1), rng.uniform(0.3, 2.0, M1)))
Gz = dict(zip(range(1, M1 + 1), rng.uniform(0.3, 2.0, M1)))
Pz, Kz = mit_faktoren(POSC, KOSC, Fz.get, Gz.get)
Hm = np.diag([Fz[m]*Gz[m] for m in MV])
Hm1 = np.diag([Fz[m - 1]*Gz[m - 1] if m > 1 else 0.0 for m in MV])
res_satz = 0.0
for x in Kz:
    for y in Pz:
        lhs = Kz[x] @ Pz[y] - Pz[y] @ Kz[x]
        rhs = (KOSC[x] @ POSC[y] - POSC[y] @ KOSC[x]) @ Hm + (POSC[y] @ KOSC[x]) @ (Hm - Hm1)
        res_satz = max(res_satz, np.abs((lhs - rhs)[:, INNEN]).max())
A0o = {k: (t["A"][k] - (t["A"][(0, 0)] + t["A"][(1, 1)])/2 if k[0] == k[1] else t["A"][k]) for k in t["A"]}
B0o = {k: (t["B"][k] - (t["B"][(0, 0)] + t["B"][(1, 1)])/2 if k[0] == k[1] else t["B"][k]) for k in t["B"]}
SPANN7 = [t["E"]] + list(A0o.values()) + list(B0o.values())
quer = POSC[(0, 0)] @ KOSC[(1, 1)]                                          # (a_1^+ a_2) (x) (b_1^+ b_2): Produkt, keine Summe
unabh = []
for m in range(2, M1 - 1):
    sel = np.where(MV == m)[0]
    V = np.array([Xs[np.ix_(sel, sel)].ravel() for Xs in SPANN7]).T
    w = quer[np.ix_(sel, sel)].ravel()
    c, *_ = np.linalg.lstsq(V, w, rcond=None)
    unabh.append(float(np.linalg.norm(V @ c - w)/np.linalg.norm(w)))
H1 = [Fr(1) + Fr(1, 10**6)]                                                 # Spurteil mit Tiefstgewicht: H(1) = lambda^2 + eps
for m in range(2, 31):
    H1.append((4*m + H1[-1]*(m - 1)**2)/Fr((m + 1)**2))
homogen = all(H1[m - 1] - 1 == Fr(4, 10**6)/(m*(m + 1))**2 for m in range(1, 31))
ZAHLEN["satz_schalenfaktoren"] = {"identitaet_residuum": res_satz, "P11_K22_nicht_im_spann_rel_residuum_je_schale": unabh}
check("[E] Satz über Schalenfaktoren: für P = F(E) P^osc, K = G(E) K^osc (beliebige Schalenfunktionen, hier zufällig) gilt "
      f"[K_ab, P_cd] = H(m) [K_ab, P_cd]^osc + (H(m) - H(m-1)) P_cd^osc K_ab^osc mit H = F G (Residuum {res_satz:.1e}); der "
      f"Quartikterm, z. B. P_11 K_22 = (a_1^+ a_2) (x) (b_1^+ b_2), liegt auf keiner Schale m >= 2 im Spann von E, A0, B0 "
      f"(rel. Residuum >= {min(unabh):.2f}); der Spurteil erzwingt mit dem Tiefstgewicht (P K = 0 auf m = 1) H(1) = "
      "lambda^2, jede Abweichung eps läuft als 4 eps/(m(m+1))^2 (exakt, Brüche). Also schließt die Familie genau für "
      "F(m) G(m) = lambda^2; mit K = P^+ heißt das F = lambda: die CCR-Normierung bis auf eine globale Konstante",
      res_satz < 1e-12 and min(unabh) > 0.1 and homogen)

# --- 2.2 Kandidaten
KAND = {}
tk = T1["Kanal"]
tn = T1["nackt"]
# Kanal mal sqrt(m+1) je Faktor: verdoppelt Faktor (m+1)
KAND["Kanal mal sqrt(m+1) je Faktor"] = mit_faktoren(tk["P"], tk["K"], lambda m: m + 1.0, lambda m: m + 1.0)
# modulare Gewichtung im KMS-Zustand: einseitig Delta = e^{-beta L_0}, Delta^s S Delta^-s
lab1 = [(m, k) for m in range(1, M1 + 1) for k in range(m)]
ix1 = {l: i for i, l in enumerate(lab1)}
S_nackt = []
for a in range(2):
    S = np.zeros((len(lab1), len(lab1)))
    for (m, k), i in ix1.items():
        if m < M1:
            S[ix1[(m + 1, SCHRITT[a](k))], i] = AMP["nackt"][a](m, k)
    S_nackt.append(S)
mod_ok = True
for s_ in (-0.5, 0.25, 0.5):
    Dl = np.diag([np.exp(-s_*BETA*m) for (m, k) in lab1])
    Dr = np.diag([np.exp(s_*BETA*m) for (m, k) in lab1])
    for S in S_nackt:
        Sw = Dl @ S @ Dr
        mod_ok &= np.allclose(Sw, np.exp(-s_*BETA)*S, rtol=1e-12, atol=0)
Delta2 = np.array([np.exp(-BETA*s[0])*np.exp(BETA*s[1]) for s in tn["lab"]])  # Delta (x) Delta-quer^-1 auf gleichen Stufen
mod_ok &= np.allclose(Delta2, 1.0)
KAND["modulare Gewichtung Delta^s S Delta^-s"] = ({k: np.exp(-BETA*0.5)*v for k, v in tn["P"].items()},
                                                  {k: np.exp(-BETA*0.5)*v for k, v in tn["K"].items()})
# Petz als K: P = Kanal, K = Petz-Kraus ((m+1)/m) A^+ (x) A-quer^+
KAND["Petz als K (P = Kanal)"] = (tk["P"], mit_faktoren(tk["P"], tk["K"], lambda m: 1.0, lambda m: (m + 1.0)/m)[1])
# geometrisches Mittel von Kanal und Petz: je Faktor ((m+1)/m)^(1/4) A, K = P^+
gm = mit_faktoren(tk["P"], tk["K"], lambda m: np.sqrt((m + 1.0)/m), lambda m: np.sqrt((m + 1.0)/m))
KAND["geometrisches Mittel Kanal/Petz"] = gm
# thermische Abbildung: Phi(e^{-beta m} I_m) = e^{-beta(m+1)} I_(m+1) mit Kraus c(m) a^+, c^2 = e^{-beta}/m
KAND["thermisch zustandserhaltend (unnormiert)"] = mit_faktoren(POSC, KOSC, lambda m: np.exp(-BETA)/m,
                                                               lambda m: np.exp(-BETA)/m)
kand_tab = {}
for name, (Pk, Kk) in KAND.items():
    d = lie_dims(Pk, Kk, MV <= M1 - 3)
    Ee = E_ext(Pk, Kk)
    swk = schalenwerte(Ee, MV, range(1, 6))
    kand_tab[name] = {"erzeugte_dimensionen": d, "E_ext_je_Schale": {m: float(f"{v[0]:.6g}") for m, v in swk.items()},
                      "skalar_je_Schale": all(v[1] < 1e-9*max(1.0, abs(v[0])) for v in swk.values())}
    say(f"   {name}: erzeugte Dimensionen {d}, E_ext je Schale {kand_tab[name]['E_ext_je_Schale']}")
ZAHLEN["kandidaten"] = kand_tab
schliessen = sorted(n_ for n_, v in kand_tab.items() if v["erzeugte_dimensionen"] == [8, 15, 15, 15])
check("[E] Kandidaten: (a) Kanal mal sqrt(m+1) je Faktor ist identisch die CCR-Normierung und schließt (8 -> 15 -> 15 -> 15, "
      "E_ext = m); (b) modulare Gewichtungen Delta^s J_-1 Delta^-s im KMS-Zustand sind k-blind: einseitig nur der konstante "
      "Faktor e^{-s beta} (s = -1/2, 1/4, 1/2, exakt), verdoppelt ist Delta (x) Delta-quer^-1 = 1 auf gleichen Stufen; das "
      "Toeplitz-Profil bleibt, kein Abschluss; (c) Petz als K, (d) geometrisches Mittel von Kanal und Petz, (e) die "
      "thermisch zustandserhaltende Abbildung haben H = F G = 1/(m(m+1)), 1/(m(m+1)), e^{-2 beta}/m^2 und schließen nicht",
      schliessen == ["Kanal mal sqrt(m+1) je Faktor"] and mod_ok
      and kand_tab["modulare Gewichtung Delta^s S Delta^-s"]["erzeugte_dimensionen"] == [8, 20, 68, 94]
      and all(kand_tab[n_]["erzeugte_dimensionen"][-1] > 15 for n_ in kand_tab if n_ not in schliessen))

# --- 2.3 Der native Möbiusheber im Fockraum (Jordan-Wigner wie z73, Fenster k = -12..11)
KMIN, KMAX = -12, 11
MODEN = range(KMIN, KMAX + 1)


def bit(k):
    return 1 << (k - KMIN)


VAK = sum(bit(k) for k in range(KMIN, 0))


def nunter(s, k):
    return bin(s & ((1 << (k - KMIN)) - 1)).count("1")


def bil(a, b, vec, coef=1):
    """coef psi*_a psi_b auf {Bitmaske: Koeffizient}."""
    out = {}
    for s, cc in vec.items():
        if not (s >> (b - KMIN)) & 1:
            continue
        sg = -1 if nunter(s, b) % 2 else 1
        s1 = s ^ bit(b)
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
    return sum(cc*v.get(s, 0) for s, cc in u.items())


def paarzahl(s):
    return sum(1 for k in range(0, KMAX + 1) if (s >> (k - KMIN)) & 1)


def paar(p, r):
    return bil(p, r, {VAK: 1})


Lnat = {-1: opsum([(Fr(k + 1), k + 1, k) for k in MODEN]), 1: opsum([(Fr(k), k - 1, k) for k in MODEN])}


def L0_op(alpha):
    return lambda vec: {s: cc*sum((alpha + k)*(((s >> (k - KMIN)) & 1) - (1 if k < 0 else 0)) for k in MODEN)
                        for s, cc in vec.items()}


def L_alpha(n, alpha):                                                      # L_n e_x = (x - n/2) e_(x-n), x = k + alpha
    return opsum([(alpha + k - Fr(n, 2), k - n, k) for k in MODEN])


def J_op(n):
    return opsum([(1, k - n, k) for k in MODEN])


ALPHAS = (Fr(1, 4), Fr(3, 4))
ein = [(p, r) for p in range(0, 6) for r in range(-6, 0)]
ok_id, ok_pz = True, True
for al in ALPHAS:
    for p, r in ein:
        v = paar(p, r)
        for n_ in (-1, 1):
            ok_id &= vadd(L_alpha(n_, al)(v), J_op(n_)(v), -(al - Fr(1, 2))) == Lnat[n_](v)
for p, r in ein:
    for n_ in (-1, 1):
        ok_pz &= all(paarzahl(s) == 1 for s in Lnat[n_](paar(p, r)))
for (p1, r1), (p2, r2) in itertools.combinations(ein, 2):
    if p1 != p2 and r1 != r2:
        v = bil(p2, r2, paar(p1, r1))
        for n_ in (-1, 1):
            ok_pz &= all(paarzahl(s) == 2 for s in Lnat[n_](v))
roh = {str(al): L_alpha(-1, al)({VAK: 1}) for al in ALPHAS}
roh_ok = all(w_ == {next(iter(paar(0, -1))): (al - Fr(1, 2))*next(iter(paar(0, -1).values()))}
             for al, w_ in zip(ALPHAS, roh.values()))
anteil = {}
for m in range(2, 8):
    tot = multi = 0
    for kr in range(-m, 0):
        w_ = J_op(-1)(paar(kr + m, kr))
        tot += sum(c*c for c in w_.values())
        multi += sum(c*c for s, c in w_.items() if paarzahl(s) != 1)
    anteil[m] = Fr(multi, tot)
ZAHLEN["mehrpaaranteil_J_-1"] = {m: str(v) for m, v in anteil.items()}
check("[E] Fockraum (Jordan-Wigner, 24 Moden, exakte Brüche): L^nat_(+-1) := L^(alpha)_(+-1) - (alpha - 1/2) J_(+-1) ist für "
      "alpha = 1/4 und 3/4 derselbe Operator Summe (k+1) psi*_(k+1) psi_k bzw. Summe k psi*_(k-1) psi_k (der Möbiusheber des "
      "Oszillatorfaktors, z58: Original = Ladung (x) Oszillatoren) und erhält die Paarzahl auf allen Ein- und Zweipaarzuständen; "
      "der rohe Twist-Heber erzeugt an der Fermikante ein Paar, L^(alpha)_-1 Omega = (alpha - 1/2) psi*_0 psi_-1 Omega, und "
      "J_-1 hat den Mehrpaaranteil (m-2)/(3m-2) (wie z80), L^nat den Anteil 0",
      ok_id and ok_pz and roh_ok and all(anteil[m] == Fr(m - 2, 3*m - 2) for m in anteil))

amp_ok, sl2_ok, tief = True, True, {}
for al in ALPHAS:
    for m in range(1, 8):
        for kr in range(-m, 0):
            p, r = kr + m, kr
            v = paar(p, r)
            w_ = Lnat[-1](v)
            hoch, runter = skalar(paar(p + 1, r), w_), skalar(paar(p, r - 1), w_)
            rest = vadd(vadd(w_, paar(p + 1, r), -hoch), paar(p, r - 1), -runter)
            amp_ok &= hoch == p + 1 and runter == -r and not rest                  # n_1 + 1 und n_2 + 1
            amp_ok &= L0_op(al)(v) == {s: m*c for s, c in v.items()}               # L_0 = m in beiden Sektoren
for m in range(1, 7):
    for kr in range(-m, 0):
        v = paar(kr + m, kr)
        c_ = vadd(Lnat[1](Lnat[-1](v)), Lnat[-1](Lnat[1](v)), -1)
        sl2_ok &= c_ == {s: 2*m*cc for s, cc in v.items()}
for m in range(1, 8):                                                       # Tiefstgewichte: Kern von L^nat_1 auf H_m
    Mx = sp.Matrix([[skalar(paar(kq + m - 1, kq), Lnat[1](paar(kr + m, kr))) for kr in range(-m, 0)]
                    for kq in range(-(m - 1), 0)]) if m > 1 else sp.zeros(0, 1)
    tief[m] = m - (Mx.rank() if m > 1 else 0)
ZAHLEN["tiefstgewichte_je_schale"] = tief
check("[E] auf der Einpaarschicht hat L^nat_-1 die Hüpfamplituden n_1 + 1 (Teilchen hoch) und n_2 + 1 (Loch runter), beide "
      "positiv, L_0 = m in beiden Sektoren, und [L^nat_1, L^nat_-1] = 2 L_0 exakt (m <= 6): eine sl(2,R) ohne Leck; je "
      "Schale genau ein Tiefstgewichtsvektor (m = 1..7), also Einpaarschicht = D_1/2 (x) D_1/2 = (+)_(h >= 1) D_h (Teilchen- "
      "und Lochleiter je eine Fermionleiter D_1/2 mit Amplitude n + 1 = sqrt((n+1)(n + 2 h)), h = 1/2)",
      amp_ok and sl2_ok and all(v == 1 for v in tief.values()))

# --- 2.4 Kanal = Möbius / L_0 und CCR = Tomita-Wurzel; Diagonalsektor
w_mob = {(m, k): (Fr(k + 1), Fr(m - k)) for m in range(1, 12) for k in range(m)}  # (n_1 + 1, n_2 + 1), n_2 = m-1-k
kanal_ok = all(w_mob[(m, k)][0]/(m + 1) == Fr(k + 1, m + 1) and w_mob[(m, k)][1]/(m + 1) == Fr(m - k, m + 1)
               and sum(w_mob[(m, k)]) == m + 1 for (m, k) in w_mob)
wurzel_ok = all(abs(AMP["CCR"][a](m, k)**2 - float(w_mob[(m, k)][a])) < 1e-12 for (m, k) in w_mob for a in range(2))
L1mat = np.zeros((len(lab1), len(lab1)))
for (m, k), i in ix1.items():
    if m < M1:
        L1mat[ix1[(m + 1, k + 1)], i] = float(w_mob[(m, k)][0])
        L1mat[ix1[(m + 1, k)], i] = float(w_mob[(m, k)][1])
delta = np.zeros((len(t["lab"]), len(lab1)))
ixd = {s: i for i, s in enumerate(t["lab"])}
for (m, k), i in ix1.items():
    delta[ixd[(m, m, k, k)], i] = 1.0
P0o = POSC[(0, 0)] + POSC[(1, 1)]
K0o = KOSC[(0, 0)] + KOSC[(1, 1)]
m1 = np.array([m for (m, k) in lab1])
diag_ok = (np.abs(P0o @ delta - delta @ L1mat).max() < 1e-12 and np.abs(K0o @ delta - delta @ L1mat.T).max() < 1e-12
           and np.abs(t["E"] @ delta - delta @ np.diag(m1.astype(float))).max() < 1e-12)
gemischt = np.abs((POSC[(0, 1)] @ delta)).max()                            # P_12 verlässt den Diagonalsektor
check("[E] (a) Die Gewichte des zustandserhaltenden Kanals (z76) sind die nativen Möbiusamplituden geteilt durch L_0 der "
      "Zielschale: p_k = (n_1 + 1)/(m+1), q_k = (n_2 + 1)/(m+1), und (n_1 + 1) + (n_2 + 1) = m + 1 macht genau daraus einen "
      "stochastischen Kanal (exakt, m <= 11). (b) Die CCR-Amplituden sind die Tomita-Wurzeln der Möbiusamplituden, "
      "|s_a|^2 = n_a + 1. (c) Auf dem Diagonalsektor delta(H) = span{|m,k> (x) |m,k>-quer} (GNS-Bild der nativen "
      "Paarlagenobservablen aus der Purifikation) gilt P^0 delta = delta L^nat_-1, K^0 delta = delta L^nat_1, E delta = "
      "delta L_0 mit P^0 = P_11 + P_22 (M = 9); P_12 führt aus dem Diagonalsektor heraus (die echt vierdimensionalen Heber)",
      kanal_ok and wurzel_ok and diag_ok and gemischt > 0.5)

# --- 2.5 z76-Petz als Casimir der radialen sl(2,R)
tk_P0 = tk["P"][(0, 0)] + tk["P"][(1, 1)]                                  # Kanal als Superoperator: Summe_a A_a (x) A-quer_a
petzK0 = mit_faktoren(tk["P"], tk["K"], lambda m: 1.0, lambda m: (m + 1.0)/m)[1]
R0 = petzK0[(0, 0)] + petzK0[(1, 1)]                                       # Petz: ((m+1)/m) Summe_a A_a^+ (x) A-quer_a^+
Emat = t["E"]
lap = np.diag([s[0]*(s[0] + 1.0) for s in t["lab"]]) @ (np.eye(len(t["lab"])) - R0 @ tk_P0)
csl2 = Emat @ Emat + Emat - K0o @ P0o
Dg = {(c_, a_): t["A"][(c_, a_)] - t["B"][(a_, c_)] for c_ in range(2) for a_ in range(2)}
Jx, Jy, Jz = (Dg[(0, 1)] + Dg[(1, 0)])/2, (Dg[(0, 1)] - Dg[(1, 0)])/2j, (Dg[(0, 0)] - Dg[(1, 1)])/2
cso3 = Jx @ Jx + Jy @ Jy + Jz @ Jz
inn = MV <= M1 - 1
spek_ok = True
for m in range(1, M1):
    sel = np.where(MV == m)[0]
    ev = np.sort(np.linalg.eigvalsh(cso3[np.ix_(sel, sel)].real))
    spek_ok &= np.allclose(ev, np.sort(np.concatenate([[l*(l + 1)]*(2*l + 1) for l in range(m)])), atol=1e-9)
check("[E] Der Petz-Rückkanal aus z76 ist die radiale sl(2,R): verdoppelt ist der Kanal Phi = Summe_a A_a (x) A-quer_a der "
      "Spurteil P^0 der Kanalheber und R Phi = K^0_osc P^0_osc / (E(E+1)); daher m(m+1)(1 - R Phi) = E(E+1) - K^0 P^0 = "
      "Casimir der sl(2,R) {E, P^0, K^0} = Casimir der diagonalen so(3) (Spektrum l(l+1), l = 0..m-1, Vielfachheit 2l+1, "
      "M = 9): z76s 'Laplace der quantisierten Sphäre' ist die Dualität radiale sl(2,R) <-> räumliche so(3) mit h = l + 1",
      np.abs((lap - csl2)[:, inn]).max() < 1e-10 and np.abs((csl2 - cso3)[:, inn]).max() < 1e-10 and spek_ok)
befund("Teil 2: Die schließende Normierung ist bedingt nativ. Nativ (E) sind der paarzahlerhaltende Möbiusheber L^nat des "
       "Oszillatorfaktors mit den Amplituden n_a + 1 und der Kanal von z76 als dessen stochastische Normierung durch L_0. "
       "Die CCR-Schritte sind die Tomita-Wurzeln dieser Amplituden; unter der einen Identifikation 'der verdoppelte "
       "Diagonalschritt ist der native Einzelschritt' (P^0 = L^nat_-1 auf dem Diagonalsektor) ist die Normierung erzwungen. "
       "Die nackten Operatoren schließen nicht (die Toeplitz-Hürde aus z73 bleibt für sie bestehen), und keine modulare "
       "Gewichtung, kein Petz-Paar und kein geometrisches Mittel ersetzt die Identifikation", True)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

# =================================================================== Teil 3
say("\nTeil 3: Die Darstellung -- masseloser Skalar (nu = 1), Minkowskienergie")
t = T1["CCR"]
A_, B_ = t["A"], t["B"]
Lx, Ly, Lz = (A_[(0, 1)] + A_[(1, 0)])/2, (A_[(0, 1)] - A_[(1, 0)])/2j, (A_[(0, 0)] - A_[(1, 1)])/2
Rx, Ry, Rz = (B_[(0, 1)] + B_[(1, 0)])/2, (B_[(0, 1)] - B_[(1, 0)])/2j, (B_[(0, 0)] - B_[(1, 1)])/2
CL, CR = Lx @ Lx + Ly @ Ly + Lz @ Lz, Rx @ Rx + Ry @ Ry + Rz @ Rz
kt_ok = True
ktab = {}
for m in range(1, M1 + 1):
    sel = np.where(MV == m)[0]
    j = (m - 1)/2
    okm = (len(sel) == m**2 and np.allclose(CL[np.ix_(sel, sel)], j*(j + 1)*np.eye(m**2))
           and np.allclose(CR[np.ix_(sel, sel)], j*(j + 1)*np.eye(m**2)) and np.allclose(np.diag(t["E"])[sel], 2*j + 1))
    ktab[m] = {"dim": len(sel), "(j,j)": f"({j}, {j})", "E": 2*j + 1, "ok": bool(okm)}
    kt_ok &= okm
qs = sp.symbols("q")
reihe = sp.series(qs*(1 + qs)/(1 - qs)**3, qs, 0, 41).removeO()
char_ok = all(reihe.coeff(qs, m) == m**2 for m in range(1, 41))
ZAHLEN["K_typen"] = ktab
check("[E] K-Typen: auf der Schale m ist der verdoppelte Raum m^2-dimensional, mit Links- und Rechts-Casimir j(j+1), "
      "j = (m-1)/2, und E = 2j + 1 = m (m = 1..9): K-Typ (j, j); der Charakter Summe m^2 q^m = q(1+q)/(1-q)^3 (Koeffizienten bis "
      "q^40) ist der Einteilchencharakter des konformen Skalars auf R x S^3 (z73), jetzt als Spur von q^E über eine einzige "
      "Darstellung", kt_ok and char_ok)

# --- Kohärenzkern
MK = 30
labK, PK_, AK_, BK_, EK_ = verdoppelt(MK, AMP["CCR"])
ixK = {s: i for i, s in enumerate(labK)}
v0 = np.zeros(len(labK), dtype=complex)
v0[ixK[(1, 1, 0, 0)]] = 1.0


def kohaerent(Z):
    G = sum(Z[a, b]*PK_[(a, b)] for a in range(2) for b in range(2))
    out, term = v0.copy(), v0.copy()
    for k in range(1, MK):
        term = G @ term/k
        out = out + term
    return out


def zufall_Z(r):
    Z = rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2))
    return Z/np.linalg.norm(Z, 2)*r


ZS = [zufall_Z(0.4), zufall_Z(0.3), zufall_Z(0.45), 0.2*np.eye(2)]
VS = [kohaerent(Z) for Z in ZS]
kern_fehler = 0.0
for (Z, vz), (W, vw) in itertools.combinations_with_replacement(list(zip(ZS, VS)), 2):
    dt = np.linalg.det(np.eye(2) - Z.conj().T @ W)
    kern_fehler = max(kern_fehler, abs(np.vdot(vz, vw)*dt - 1))
ZAHLEN["kohaerenzkern_rel_fehler"] = kern_fehler
check(f"[E] Kohärenzzustände |Z> = exp(Summe Z_ab P_ab) |m = 1> (M = 30, vier Punkte der Matrixkugel, zehn Paare): "
      f"<Z|W> = det(I - Z^+ W)^-1 (rel. Fehler {kern_fehler:.1e}), der nu = 1-Kern aus z74/z79: die Darstellung ist die "
      "Leiterdarstellung des masselosen Skalars", kern_fehler < 1e-12)

# --- Casimir: Skalar und volle Standardform (alle Helizitäten)
KINV = np.linalg.inv(KILL)
Ccas = sum(KINV[i, j]*X15[i] @ X15[j] for i in range(15) for j in range(15))
cblk = Ccas[np.ix_(MV <= M1 - 1, MV <= M1 - 1)]
cas_skalar = (float(np.real(cblk.diagonal().mean())), float(np.abs(cblk - np.diag(cblk.diagonal())).max()),
              float(np.ptp(cblk.diagonal().real)))
Cfund = sum(KINV[i, j]*Y15[i] @ Y15[j] for i in range(15) for j in range(15))
MV9 = 9
labV, PV, AV, BV, EV = verdoppelt(MV9, AMP["CCR"], voll=True)              # dünn besetzt; volle Standardform
KV = {k: v.T.tocsr() for k, v in PV.items()}
XV = basis15(dict(P=PV, K=KV, A=AV, B=BV, E=EV))
mV = np.array([s[0] for s in labV])
mbV = np.array([s[1] for s in labV])
spV = np.where((mV <= MV9 - 2) & (mbV <= MV9 - 2))[0]
resV = 0.0                                                                  # dieselben Strukturkonstanten wie in Teil 1
for i in range(15):
    for j in range(15):
        R_ = (XV[i] @ XV[j] - XV[j] @ XV[i] - sum(f15[i, j, k]*XV[k] for k in range(15) if abs(f15[i, j, k]) > 1e-14))
        resV = max(resV, float(abs(R_[:, spV]).max()) if R_[:, spV].nnz else 0.0)
innenV = np.where((mV <= MV9 - 1) & (mbV <= MV9 - 1))[0]
CV = sum(KINV[i, j]*(XV[i] @ XV[j][:, innenV]) for i in range(15) for j in range(15) if abs(KINV[i, j]) > 1e-14)
CV = CV.toarray()                                                           # Spalten innenV
Kalle = sps.vstack([KV[k] for k in KV]).tocsc()
EVd = EV.diagonal()
hel = {}
hel_ok = True
for zwei_h in range(-4, 5):
    h = zwei_h/2
    pos = np.where((mV[innenV] - mbV[innenV]) == zwei_h)[0]
    blk = CV[innenV[pos]][:, pos]
    wert = float(np.real(blk.diagonal().mean()))
    sektor = np.where(mV - mbV == zwei_h)[0]
    u, s_, vh = np.linalg.svd(Kalle[:, sektor].toarray())
    kern = vh[int(np.sum(s_ > 1e-10)):]
    E_kern = sorted({float(EVd[sektor][i]) for v_ in kern for i in np.where(np.abs(v_) > 1e-10)[0]})
    hel[h] = {"casimir": wert, "soll_3(h^2-1)/8": 3*(h**2 - 1)/8, "tiefstgewichte": len(kern), "E_tiefst": E_kern}
    hel_ok &= (abs(wert - 3*(h**2 - 1)/8) < 1e-10 and np.abs(blk - wert*np.eye(len(pos))).max() < 1e-10
               and len(kern) == 2*abs(h) + 1 and E_kern == [abs(h) + 1])
ZAHLEN["casimir"] = {"skalar_(mittel, nichtdiagonal, spannweite)": cas_skalar, "fundamental": float(np.real(Cfund[0, 0])),
                     "helizitaeten": hel}
say(f"   Helizitätssektoren h = (m - mquer)/2: {{h: (Casimir, Tiefstgewichte, E_tiefst)}} = "
    f"{ {h: (round(v['casimir'], 6), v['tiefstgewichte'], v['E_tiefst']) for h, v in hel.items()} }")
check(f"[E] Casimir (Killing-Normierung): auf der verdoppelten Leiter überall {cas_skalar[0]:.6f} = -3/8 = Delta(Delta-4)/8 mit "
      f"Delta = 1 (Skalar), auf der Fundamentaldarstellung {np.real(Cfund[0, 0]):.6f} = 15/32; auf der vollen Standardform "
      "(+)_(m, mquer) H_m (x) H-quer_mquer (M = 9) zerfällt der Raum nach h = (m - mquer)/2 in Sektoren mit Casimir "
      "3(h^2 - 1)/8 und genau einem Tiefstgewichtstyp der Dimension 2|h| + 1 bei E = |h| + 1 (h = -2..2): jede masselose "
      "Helizität genau einmal", abs(cas_skalar[0] + 0.375) < 1e-12 and cas_skalar[1] < 1e-12 and cas_skalar[2] < 1e-12
      and abs(np.real(Cfund[0, 0]) - 15/32) < 1e-12 and np.allclose(Cfund, Cfund[0, 0]*np.eye(4)) and resV < 1e-12 and hel_ok)

# --- Minkowskienergie: parabolisch, positiv, Gamma-Maße
hvek = 1j*vek({"E": 1, "P11": 0.5, "P22": 0.5, "K11": 0.5, "K22": 0.5})        # reelle Form: i mal hermitesch
ad_h = sum(hvek[i]*AD[i] for i in range(15))
evek = 1j*vek({"E": 1})
B_hh, B_ee = float(np.real(hvek @ KILL @ hvek)), float(np.real(evek @ KILL @ evek))
nilp = float(np.abs(ad_h @ ad_h @ ad_h).max())
nilp2 = float(np.abs(ad_h @ ad_h).max())
MM = 40
labM, PM, AM, BM, EM = verdoppelt(MM, AMP["CCR"])
mM = np.array([s[0] for s in labM])
P0M = PM[(0, 0)] + PM[(1, 1)]
HM = (EM + 0.5*(P0M + P0M.T)).tocsr()
DgM = [AM[(0, 1)] - BM[(1, 0)], AM[(1, 0)] - BM[(0, 1)], (AM[(0, 0)] - AM[(1, 1)]) - (BM[(0, 0)] - BM[(1, 1)])]
innenM = np.where(mM <= MM - 1)[0]
kom_rot = max(abs((HM @ D - D @ HM)[:, innenM]).max() for D in DgM)
minev = {}
for Mt in (10, 20, 40):
    labt, Pt, _, _, Et = verdoppelt(Mt, AMP["CCR"])
    P0t = Pt[(0, 0)] + Pt[(1, 1)]
    Ht = (Et + 0.5*(P0t + P0t.T)).tocsr()
    if Ht.shape[0] < 3000:
        ev = np.linalg.eigvalsh(Ht.toarray())
        minev[Mt] = (float(ev[0]), int(np.sum(ev < 1)), float(Mt*ev[0]))
    else:
        ev = np.sort(spla.eigsh(Ht, k=4, sigma=-0.01, which="LM", return_eigenvectors=False))
        minev[Mt] = (float(ev[0]), None, float(Mt*ev[0]))
K0M = P0M.T.tocsr()
gamma_fehler = {}
for l in range(3):
    sel = np.where(mM == l + 1)[0]
    if l == 0:
        vl = np.zeros(len(labM))
        vl[sel[0]] = 1.0
    else:
        u, s_, vh = np.linalg.svd(K0M[:, sel].toarray())
        vl = np.zeros(len(labM))
        vl[sel] = vh[int(np.sum(s_ > 1e-10))]
    gamma_fehler[l] = max(abs(np.vdot(vl, spla.expm_multiply(1j*tt*HM, vl.astype(complex)))
                              - (1 - 1j*tt/2)**(-2*(l + 1))) for tt in (0.5, 1.0, 2.0))
ZAHLEN["minkowskienergie"] = {"B(H,H)": B_hh, "B(E,E)": B_ee, "max|(ad H)^3|": nilp, "max|(ad H)^2|": nilp2,
                              "[H, Drehungen]": float(kom_rot), "min_eig_(wert, anzahl<1, M*min)": minev,
                              "gamma_fehler_l": gamma_fehler}
say(f"   Minkowskienergie: kleinste Eigenwerte {minev}; Gamma-Fehler {gamma_fehler}")
check(f"[E] Minkowskienergie P^M_0 = E + (P^0 + K^0)/2 (P^0 = P_11 + P_22 = L_-1, K^0 = L_1, [K^0, P^0] = 2E wie z72): "
      f"parabolisch (Killing-Form der reellen Form: B(i P^M_0, i P^M_0) = {B_hh:.1e}, (ad P^M_0)^3 = 0, (ad P^M_0)^2 != 0; "
      f"dagegen B(iE, iE) = {B_ee:.0f} < 0: E elliptisch), vertauscht mit den räumlichen Drehungen A_ca - B_ac ({kom_rot:.1e}); alle Trunkierungen positiv, "
      f"der kleinste Eigenwert fällt wie 1/M (M min = {minev[10][2]:.3f}, {minev[20][2]:.3f}, {minev[40][2]:.3f} gegen "
      "j_11^2/8 = 1,8353), und die Spektralmaße der tiefsten Zustände mit l = 0, 1, 2 sind (1 - it/2)^(-2(l+1)) (Gamma, "
      f"absolut stetig auf (0, inf); max. Fehler {max(gamma_fehler.values()):.1e} bei M = 40): Spektrum [0, inf), positiv, "
      "auf jedem Drehimpulssektor die Laguerre-Leiter D_(l+1)",
      abs(B_hh) < 1e-10 and nilp < 1e-10 and nilp2 > 0.1 and B_ee < -1 and kom_rot < 1e-10
      and all(v[0] > 0 for v in minev.values()) and minev[10][0] > minev[20][0] > minev[40][0]
      and 1.6 < minev[40][2] < 1.84 and max(gamma_fehler.values()) < 1e-10)
befund("Teil 3: Mit der CCR-Normierung ist die ganze Schalenleiter EINE Darstellung: der masselose konforme Skalar "
       "(K-Typen (j, j) bei E = 2j + 1, Charakter Summe m^2 q^m, Kern nu = 1, Casimir -3/8). Die Moden m sind seine "
       "Energiestufen E = L_0 = m, keine Spezies. Die volle Standardform der Einpaarschicht (ungleiche Stufen m != mquer) "
       "enthält zusätzlich jede masselose Helizität h = (m - mquer)/2 genau einmal; der Modulgenerator des Originals ist "
       "dort beta (L_0 - L-quer_0) = 2 beta h, also die Helizität (wie z79: Zentrum = Modulfluss). Die Minkowskizeit ist der "
       "parabolische Erzeuger E + (P^0 + K^0)/2 mit Spektrum [0, inf)", True)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

# =================================================================== Teil 4
say("\nTeil 4: Thermische Gewichte des Originals gegen den thermischen Zustand der CFT auf S^3")
bs, bp, cs, xs_ = sp.symbols("beta beta_p c x", positive=True)
cfun = lambda x: 1/(1 + sp.exp(bs*x))
gew_ok = True
for al in (sp.Rational(1, 4), sp.Rational(3, 4)):
    for m in range(1, 7):
        for kr in range(-m, 0):
            xp, xr = kr + m + al, kr + al
            g_ = sp.simplify(cfun(xp)*(1 - cfun(xr))/((1 - cfun(xp))*cfun(xr)))
            gew_ok &= sp.simplify(g_ - sp.exp(-bs*m)) == 0
w1 = lambda m: m*sp.exp(-bs*m)                                              # einseitig, je Schale
wcft = lambda m: cs*m**2*sp.exp(-bp*m)
los = sp.solve([sp.Eq(w1(1), wcft(1)), sp.Eq(w1(2), wcft(2))], [cs, bp], dict=True)
los_ok = len(los) == 1 and sp.simplify(los[0][bp] - bs - sp.log(2)) == 0
abweichung3 = sp.simplify(wcft(3).subs(los[0])/w1(3)) if los_ok else None
check("[negativ, E] einseitig: Die Kompression des KMS-Zustands (beta = 4 pi, symbolisch in beta, beide Sektoren, m <= 6) "
      "gibt jedem Paarzustand das Gewicht e^{-beta m}, je Schale m e^{-beta m} (Stromcharakter h = 1, z72). Ein CFT-Gewicht "
      f"c m^2 e^{{-beta' m}} passt an m = 1, 2 nur mit beta' = beta + ln 2, c = 2 und verfehlt m = 3 um den Faktor "
      f"{abweichung3}: keine Temperatur, keine Normierung (Vielfachheit m statt m^2)",
      gew_ok and los_ok and abweichung3 == sp.Rational(3, 4))
prod_ok = all(sp.simplify(m**2*sp.exp(-bs*m)*sp.exp(-bs*m) - (m**2*sp.exp(-bp*m)).subs(bp, 2*bs)) == 0 for m in range(1, 31))
ZAHLEN["thermisch"] = {"einseitig_fit_m12": {str(k): str(v) for k, v in los[0].items()}, "verfehlung_m3": str(abweichung3),
                       "produkt_beta_strich": "2 beta = 8 pi"}
check("[E] Produktlesart: Original und Tomita-Partner je im KMS-Zustand bei beta = 4 pi, auf gleiche Stufen beschränkt, geben "
      "jedem der m^2 Zustände der Schale das Gewicht e^{-2 beta m}: genau der thermische Einteilchenzustand der CFT auf S^3, "
      "m^2 e^{-beta' m}, mit beta' = 2 beta = 8 pi und Normierung 1 relativ zur Fermisee (m <= 30, symbolisch); gleichwertig: "
      "KMS bezüglich L_0 + L-quer_0 bei 4 pi ist KMS bezüglich E = (L_0 + L-quer_0)/2 bei 8 pi", prod_ok)

q_ = np.exp(-BETA/2)
tfd = np.array([q_**(s[0] - 1) if s[2] == s[3] else 0.0 for s in labK], dtype=complex)
koh_q = kohaerent(q_*np.eye(2))
DgK = [AK_[(0, 1)] - BK_[(1, 0)], AK_[(1, 0)] - BK_[(0, 1)], (AK_[(0, 0)] - AK_[(1, 1)]) - (BK_[(0, 0)] - BK_[(1, 1)])]
zonal = max(np.abs(D @ tfd).max() for D in DgK)
schale_tfd = {m: float(sum(abs(tfd[i])**2 for i, s in enumerate(labK) if s[0] == m)/q_**(2*(m - 1))) for m in range(1, 6)}
kap2 = (1 - q_)/(1 + q_)
Delta_voll = sps.diags([np.exp(-BETA*(s[0] - s[1])) for s in labV]).tocsr()
mod_zentral = max(float(abs(Delta_voll @ X - X @ Delta_voll).max()) for X in XV)
ZAHLEN["purifikation"] = {"q": q_, "kappa^2 = (1-q)/(1+q)": kap2, "tanh(pi)": float(np.tanh(np.pi)),
                          "schalengewicht/q^(2(m-1))": schale_tfd, "norm": float(np.vdot(koh_q, koh_q).real)}
check("[E] Purifikation: Der Einpaaranteil des Araki-Woods-Vektors ist Summe_m q^(m-1) Summe_k |m,k> (x) |m,k>-quer mit "
      f"q = e^(-beta/2) = e^(-2 pi), exakt der Kohärenzzustand exp(q P^0)|m = 1> am Punkt Z = q I der Matrixkugel (Abweichung "
      f"{np.abs(tfd - koh_q).max():.0e}), Norm det(I - q^2 I)^-1 = (1 - q^2)^-2, rein zonal (l = 0, von allen räumlichen "
      f"Drehungen vernichtet: {zonal:.0e}), Schalengewicht m q^(2(m-1)); derselbe Parameter q wie bei den "
      f"Araki-Woods-Oszillatoren aus z78/z79 (kappa^2 = (1-q)/(1+q) = tanh(pi) = {kap2:.6f}). Der Modulgenerator ist auf der vollen Standardform beta (m - mquer) = 2 beta h, "
      f"vertauscht mit allen 15 Erzeugern ({mod_zentral:.0e}) und ist auf dem Skalar (h = 0) null: Modulfluss = "
      "Helizitätsphase, keine Raumzeitbewegung (GRAV.NONCIRCULAR.01)",
      np.abs(tfd - koh_q).max() < 1e-15 and abs(np.vdot(koh_q, koh_q).real - (1 - q_**2)**-2) < 1e-14 and zonal < 1e-15
      and all(abs(v - m) < 1e-12 for m, v in schale_tfd.items()) and abs(kap2 - np.tanh(np.pi)) < 1e-15
      and mod_zentral < 1e-12)
befund("Teil 4: Die thermischen Gewichte des Originals passen auf genau eine Weise zum thermischen Zustand der CFT auf S^3: "
       "als Produkt von Original und Partner (je beta = 4 pi) mit beta' = 8 pi. Einseitig gibt es keine Übereinstimmung, und "
       "die Purifikation ist kein thermischer S^3-Zustand, sondern der zonale Kohärenzzustand am Punkt Z = e^(-2 pi) I", True)

# =================================================================== Teil 5
say("\nTeil 5: Die Uhr R = e^{2 pi i L_0} und der Spiegel zwischen den Quartersektoren")
R_ok = True
for al in ALPHAS:
    for m in range(1, 8):
        for kr in range(-m, 0):
            v = paar(kr + m, kr)
            L0v = L0_op(al)(v)
            R_ok &= all(c_/v[s] == m for s, c_ in L0v.items())              # ganzzahlig: R = 1 relativ zu Omega
eins = {}
for al in ALPHAS:
    for p in range(0, 4):
        v = {VAK | bit(p): 1}                                               # psi*_p Omega (geladen)
        L0v = L0_op(al)(v)
        eins[(str(al), p)] = complex(np.exp(2j*np.pi*float(list(L0v.values())[0])))
R_eins_ok = all(abs(z - (1j if k[0] == "1/4" else -1j)) < 1e-12 for k, z in eins.items())
schritt_ok = True
for al in ALPHAS:
    for n_, op in ((-1, Lnat[-1]), (-1, J_op(-1))):
        for p, r in ein[:20]:
            v = paar(p, r)
            L0v = L0_op(al)(v)[next(iter(v))]/v[next(iter(v))]
            for s, c_ in op(v).items():
                L0w = L0_op(al)({s: 1})[s]
                schritt_ok &= L0w - L0v == -n_                              # jeder Schritt hebt L_0 um genau 1
R_dublett = np.diag([np.exp(2j*np.pi*1), np.exp(2j*np.pi*1)])              # Phase von Ad_R auf (S_1, S_2)
ordnung = next(k for k in range(1, 9) if np.allclose(np.linalg.matrix_power(R_dublett, k), np.eye(2)))
ZAHLEN["uhr"] = {"R_auf_geladenen_einzelzustaenden": {f"{k[0]}, p={k[1]}": str(np.round(z, 12)) for k, z in eins.items()},
                 "ordnung_auf_schrittdublett": ordnung, "z79_vergleich": "Quarterpaar Ordnung 2, V_1 Ordnung 4"}
check("[E] Uhr: R = e^{2 pi i L_0} ist auf allen Paarzuständen (beide Sektoren, m <= 7) relativ zu Omega gleich 1 (L_0 = m "
      "ganzzahlig, R^2 = (-1)^F = 1 bei gerader Fermionzahl); auf geladenen Einzelanregungen psi*_p Omega ist R = i "
      "(alpha = 1/4) bzw. -i (alpha = 3/4). Beide Elementarschritte (L^nat_-1 und J_-1, also auch ihre Hardy-Hälften) heben L_0 "
      f"um genau 1, also Ad_R = 1: auf dem Schrittdublett (Spinorindex a von S2) hat R die Ordnung {ordnung}, gegen 2 auf dem "
      "fermionischen Quarterpaar und 4 auf V_1 (z79). R ist auf der ganzen verdoppelten Leiter die Eins und vertauscht mit "
      "su(2,2): die Raumzeit von S2 ist R-neutral, verträglich mit S1 (R im Ladungsteil)",
      R_ok and R_eins_ok and schritt_ok and ordnung == 1)

xsym = sp.symbols("x", real=True)
spiegel_c = sp.simplify(cfun(-xsym) - (1 - cfun(xsym))) == 0
spiegel_ok = spiegel_c
for m in range(1, 8):
    for kr in range(-m, 0):
        p, r = kr + m, kr
        pb, rb = -r - 1, -p - 1                                             # Bild: Teilchen bei -r-1, Loch bei -p-1
        n1, n2 = p, -r - 1
        n1b, n2b = pb, -rb - 1
        w_ = Lnat[-1](paar(pb, rb))
        spiegel_ok &= (n1b, n2b) == (n2, n1)
        spiegel_ok &= skalar(paar(pb + 1, rb), w_) == n2 + 1 and skalar(paar(pb, rb - 1), w_) == n1 + 1
        spiegel_ok &= (pb + Fr(3, 4)) - (rb + Fr(3, 4)) == (p + Fr(1, 4)) - (r + Fr(1, 4)) == m
check("[E] Spiegel x -> -x (Teilchen-Loch, Moden k -> -k-1): er ist eine Symmetrie des Originals (c(-x) = 1 - c(x), symbolisch), "
      "führt den Sektor alpha = 1/4 in alpha = 3/4 über, erhält L_0 = m und die KMS-Gewichte und vertauscht die Schrittarten: "
      "(n_1, n_2) -> (n_2, n_1), die Möbiusamplituden gehen in die vertauschten über (m <= 7). Die beiden Quarterkopien sind "
      "also dieselbe su(2,2)-Darstellung, verbunden durch C mal Schrittvertauschung (ein Weyl-Element)", spiegel_ok)
befund("Teil 5: Der Spinorindex von S2 ist die Schrittart (Teilchen hoch / Loch runter), nicht das Quarter. Er trägt keine "
       "R-Ladung. Das Quarter bleibt als zweifache Kopie, die der native Spiegel C verbindet: Lesart als Teilchen und "
       "Antiteilchen eines komplexen Skalars (L)", True)

# =================================================================== Teil 6
say("\nTeil 6: Bilanz")
STATUS = {
    "Abschluss der Heber": ("E", "su(2,2) schließt auf der verdoppelten Schalenleiter genau für die CCR-Normierung (Faktor 1, "
                                 "E = L_0); innerhalb aller Schalenfunktionen F(E) eindeutig bis auf eine Konstante"),
    "nackte Hardy-Hälften, Kanal, modulare Gewichtung, Petz, geometrisches Mittel": (
        "E (negativ)", "schließen nicht (Toeplitz-Randprojektoren bzw. falscher Schalenfaktor; Delta ist innerhalb einer "
                       "Schale k-blind)"),
    "nativer Möbiusheber L^nat": ("E", "L_-1 - (alpha - 1/2) J_-1 = Summe (k+1) psi*_(k+1) psi_k, paarzahlerhaltend, "
                                       "Amplituden n_a + 1, sl(2,R) auf der Einpaarschicht = D_1/2 (x) D_1/2; Kanal (z76) = "
                                       "L^nat_-1 / L_0"),
    "CCR = Tomita-Wurzel": ("C", "unter der Identifikation P^0 = L^nat_-1 auf dem Diagonalsektor ist die schließende "
                                 "Normierung erzwungen"),
    "ein Feld statt Spezies": ("C", "die Moden m sind die Energiestufen E = m eines masselosen Skalars je Quartersektor; "
                                    "die zwei Quarter sind ein C-Paar (L)"),
    "Darstellung": ("E", "K-Typen (j, j), Summe m^2 q^m, Kern nu = 1, Casimir -3/8 (gegeben die CCR-Normierung)"),
    "Minkowskizeit": ("E", "parabolisch E + (P^0 + K^0)/2, Spektrum [0, inf), Gamma-Maße D_(l+1) (gegeben CCR)"),
    "z76-Petz": ("E", "m(m+1)(1 - R Phi) = Casimir der radialen sl(2,R) = l(l+1)"),
    "thermisch": ("E", "Produkt Original x Partner = CFT-Einteilchenzustand bei beta' = 8 pi; einseitig keine "
                       "Übereinstimmung; Purifikation = Kohärenzzustand bei Z = e^(-2 pi) I"),
    "Uhr": ("E", "R trivial auf Schritten und Schalen; Spinorindex = Schrittart, R-neutral"),
    "Helizitätsturm": ("E / L", "volle Standardform: jede masselose Helizität einmal, Modulfluss = Helizität; physikalische "
                                "Lesart (Fermionen, Spin 2) offen"),
    "Herleitung der Identifikation": ("O", "ein Prinzip, das den verdoppelten Diagonalschritt mit dem nativen Einzelschritt "
                                           "gleichsetzt (vollständig positive Fortsetzung, Standardform)"),
    "Dynamik": ("O", "alles Einteilchenkinematik; Wechselwirkung und zweite Quantisierung der Mehrpaarschichten offen "
                     "(DIMENSION.UPLIFT.FIREWALL.01)"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
urteil = ("S2 hält bedingt. Auf der Tomita-verdoppelten Schalenleiter (+) H_m (x) H-quer_m sind [E, P] = P und [P, P] = 0 für jede "
          "Normierung der beiden Elementarschritte richtig, aber [K, P] schließt genau dann zu su(2,2), wenn die Schritte die "
          "CCR-Amplituden sqrt(n_a + 1) tragen; dann mit Faktor 1 und E = L_0. Die nackten Hardy-Hälften von J_-1 scheitern an "
          "Randprojektoren (Toeplitz-Hürde aus z73, erzeugte Algebra 8 -> 20 -> 68 -> 94), der zustandserhaltende Kanal am "
          "Schalenfaktor (E_ext = (2m-1)/(4m^2), 8 -> 24 -> 104 -> 155); modulare Gewichtungen sind innerhalb einer Schale "
          "k-blind, Petz als K und das geometrische Mittel scheitern ebenfalls. Allgemein schließt P = F(E) P^osc, K = G(E) "
          "K^osc genau für F G = const. Neu ist ein nativer Ursprung der Normierung: Der Möbiusheber des Oszillatorfaktors, "
          "L^nat_-1 = L_-1 - (alpha - 1/2) J_-1 = Summe (k+1) psi*_(k+1) psi_k, erhält die Paarzahl exakt (J_-1 hat dagegen den "
          "Mehrpaaranteil (m-2)/(3m-2)), hat auf der Einpaarschicht die Amplituden n_a + 1 und bildet mit L_0 eine sl(2,R) = "
          "D_1/2 (x) D_1/2. Der Kanal aus z76 ist genau L^nat_-1 geteilt durch L_0 der Zielschale, und die CCR-Schritte sind die "
          "Tomita-Wurzeln von L^nat: Auf dem Diagonalsektor gilt P^0 = L^nat_-1, K^0 = L^nat_1, E = L_0. Unter dieser einen "
          "Identifikation ist die Normierung erzwungen (C). Dann ist die Schalenleiter eine einzige Darstellung, der masselose "
          "konforme Skalar (K-Typen (j, j) bei E = 2j + 1, Summe m^2 q^m, Kern det(I - Z^+ W)^-1, Casimir -3/8); die Moden sind "
          "Energiestufen, das Speziesproblem (eine Kopie je m) verschwindet, übrig bleibt das Quarterpaar als C-Paar. Die "
          "Minkowskienergie E + (P^0 + K^0)/2 ist parabolisch und positiv mit Spektrum [0, inf). z76s Petz-Laplace ist der "
          "Casimir der radialen sl(2,R). Thermisch passt nur die Produktlesart (Original und Partner je bei 4 pi) zur CFT auf "
          "S^3, bei beta' = 8 pi; die Purifikation ist der zonale Kohärenzzustand bei Z = e^(-2 pi) I, ihr Modulfluss ist die "
          "Helizitätsphase. R ist auf den Schritten trivial, die Raumzeit von S2 ist R-neutral. Für K12 heißt das: Die Gruppe "
          "über alle Moden ist von einer Lesart (L) zu einer bedingten Herleitung (C) geworden; offen bleiben die Herleitung der "
          "Identifikation aus einem Prinzip und die ganze Dynamik.")
say("\nUrteil\n   " + urteil)
say(f"\n   (Laufzeit {time.time() - T0:.1f} s)")
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z83_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
