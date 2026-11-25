"""Z122: Folgeschritt zu T5 'Ist die Linienregel P_RJ nativ? Die falschen gekleideten Vertizes und die native Paarung'.

Leitlinie 'Rekonstruktion statt Auswahl': eine Paarung zählt nur als nativ, wenn vorhandene TFPT-Struktur die andere ausschließt.
Firewalls: DIMENSION.UPLIFT.FIREWALL.01 (jede 3+1-Aussage ist Kinematik der radialen Einpaarschicht oder Lesart),
GRAV.NONCIRCULAR.01. H6 und die gelösten Bausteine (Flavour, Massen, alpha, die 16, E8) bleiben fest; nichts wird befördert.
Bestand: z97 (gekleidete Vertizes Psi_a = V^+ (x) S_a mit R = i und h + 1/2, Psi-quer_b = V^- (x) S-quer_b mit R = -i und
h - 1/2; G = U_2pi R^2, P_GSO; S_a ist (1/2, 0)-, S-quer_b (0, 1/2)-Spinor), z106 (P_CAR: gemeinsame CAR der Nähte, gekleidete
Vertizes verschiedener Nähte antikommutieren), z114 (Linienregeln: 256 -> 6 -> 2 J-kovariante, P_RJ und Gegenorientierung;
unter P_RJ + P_GSO + K2 Netto-Index 3 je 16, Spiegelpaar als nahtsterile Leiter im Spektrum; Empfehlung: Nativität von P_RJ
über die gekleideten Vertizes prüfen), z82 (R = i^(Q_a - Q_b) = e^(2 pi i L_0), Ladungsenergie Q^2/2 + (alpha - 1/2) Q,
eta-Dublett, niedrigste Sektoren), z83 (Standardform, Schrittart = Teilchen hoch / Loch runter, Modulgenerator 2 beta h,
Purifikation, nativer Heber L^nat_-1, Minkowskizeit P_0), z85/z89 (Tomita J = Twistorkonjugation, Modulladung D = N_a - N_b
= 2h dreht T), z116 (R = i Gamma_11, S^+ = (S_L x 16) + (S_R x 16-quer), Verklebung I_T = I_F offen), z102 (jeder
Z4-Helizitätscharakter, der h = +-1/2 trennt, trennt h = 0 von h = 1), z103 (native Orientierungsdaten, ein Bit).

Bezeichnungen: T = Q_a - Q_b (Transferladung, R = i^T), 2h = m - mquer, X = R i^(-2h) = i^N mit N = T - 2h, G = (-1)^N.
Seiten: o = Original (linkes Blatt m, Gerade P(omega), sieht h >= 0), p = Partner (rechtes Blatt mquer, P(pi), sieht h <= 0).
Gekleidete Vertizes: nativ Psi = V^+ (x) S (Klasse (i, +1/2)), Psi-quer = V^- (x) S-quer ((-i, -1/2)); falsch Phi = V^+ (x)
S-quer ((i, -1/2)), Phi-quer = V^- (x) S ((-i, +1/2)). V^+ = R-Schritt i (psi*_a oder psi_b), V^- = R-Schritt -i.

Teil 1  Sind die falschen Vertizes ausgeschlossen? Modell A (eine Naht, Ladung x Standardform): Klassen, G, X. CAR: zwei
        Nähte in der gemeinsamen CAR (z106), Blatttausch Theta. Modeninhalt: Maya-Fockraum des Originals (beide Quarter),
        welche Seite welchen Spinor trägt. Tomita: Modell B (Original und Partner mit Ladungen, Standardform, KMS-Vektor),
        J und S auf allen acht Elementarvertizes. Energiepositivität und Lokalität: E, P_0, Modulgewicht, Seitenlokalität,
        Ladungsenergie. Was die native Paarung ist (Z4-Verklebung, I_T = I_F, eta_+-Kleidung).
Teil 2  Wählt die native Paarung P_RJ? z114s Linienregeln nachgerechnet (6 -> 2), native Paarung als Bedingung (-> 1),
        Vorzeichen, Kontrollen (Spiegelpaarung, vektorartig, Modeninhalt, Modulpositivität).
Teil 3  Folgen für K16 unter T + P_GSO + K2 (+ P_RJ aus der Paarung): Superauswahl N = T - 2h, was vom Spiegel im Spektrum
        bleibt, Preis einer spektralen Z4-Projektion; Bilanz, Urteil, Folgen für K16 und die Achse chiral_matter.
"""
import itertools
import json
import math
from collections import Counter, deque
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
import scipy.sparse as sps

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
rng = np.random.default_rng(122)
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


def nrm(X):
    """Größter Betrag eines dünnen oder dichten Operators (0 für den Nulloperator)."""
    if sps.issparse(X):
        X = X.tocsr()
        X.eliminate_zeros()
        return float(abs(X).max()) if X.nnz else 0.0
    return float(np.abs(X).max()) if np.size(X) else 0.0


def kr(*ops):
    out = sps.csr_matrix(ops[0])
    for o in ops[1:]:
        out = sps.kron(out, sps.csr_matrix(o), format="csr")
    return out


# ------------------------------------------------------------------ radiale Standardform (wie z83/z97, CCR-Schritte)
AMP_CCR = (lambda m, k: np.sqrt(k + 1), lambda m, k: np.sqrt(m - k))
SCHRITT = (lambda k: k + 1, lambda k: k)                          # S_0 Teilchen hoch, S_1 Loch runter (z83)


def standardform(M):
    """(+)_(m, mquer <= M) H_m (x) H-quer_mquer wie z97 (dünn): P[a,b] = S_a (x) S-quer_b, A, B Bilineare links/rechts, E,
    einseitige Schritte SL[a] = S_a (x) 1, SR[b] = 1 (x) S-quer_b und der native Heber L^nat_-1 des linken Blatts (z83:
    Amplituden n_1 + 1 Teilchen hoch, n_2 + 1 Loch runter)."""
    lab = [(m, mb, k, kb) for m in range(1, M + 1) for mb in range(1, M + 1) for k in range(m) for kb in range(mb)]
    ix = {s: i for i, s in enumerate(lab)}
    n = len(lab)

    def dn(r, c, w):
        return sps.csr_matrix((w, (r, c)), shape=(n, n))
    P = {}
    for a in range(2):
        for b in range(2):
            r, c, w = [], [], []
            for (m, mb, k, kb), i in ix.items():
                if m < M and mb < M:
                    r.append(ix[(m + 1, mb + 1, SCHRITT[a](k), SCHRITT[b](kb))])
                    c.append(i)
                    w.append(AMP_CCR[a](m, k)*AMP_CCR[b](mb, kb))
            P[(a, b)] = dn(r, c, w)

    def bilinear(c_, a_, seite):
        r, cc, w = [], [], []
        for (m, mb, k, kb), i in ix.items():
            mm, kk = (m, k) if seite == 0 else (mb, kb)
            k1 = kk - (1 if a_ == 0 else 0)
            if mm == 1 or not 0 <= k1 <= mm - 2:
                continue
            k2 = SCHRITT[c_](k1)
            r.append(ix[(m, mb, k2, kb) if seite == 0 else (m, mb, k, k2)])
            cc.append(i)
            w.append(AMP_CCR[a_](mm - 1, k1)*AMP_CCR[c_](mm - 1, k1))
        return dn(r, cc, w)
    A = {(c_, a_): bilinear(c_, a_, 0) for c_ in range(2) for a_ in range(2)}
    B = {(c_, a_): bilinear(c_, a_, 1) for c_ in range(2) for a_ in range(2)}
    SL, SR = {}, {}
    for a in range(2):
        r, c, w, r2, c2, w2 = [], [], [], [], [], []
        for (m, mb, k, kb), i in ix.items():
            if m < M:
                r.append(ix[(m + 1, mb, SCHRITT[a](k), kb)])
                c.append(i)
                w.append(AMP_CCR[a](m, k))
            if mb < M:
                r2.append(ix[(m, mb + 1, k, SCHRITT[a](kb))])
                c2.append(i)
                w2.append(AMP_CCR[a](mb, kb))
        SL[a], SR[a] = dn(r, c, w), dn(r2, c2, w2)
    r, c, w = [], [], []
    for (m, mb, k, kb), i in ix.items():
        if m < M:
            for k2, amp in ((k + 1, k + 1), (k, m - k)):
                r.append(ix[(m + 1, mb, k2, kb)])
                c.append(i)
                w.append(float(amp))
    LNAT = dn(r, c, w)
    E = sps.diags([0.5*(s[0] + s[1]) for s in lab]).tocsr()
    return lab, ix, P, A, B, E, SL, SR, LNAT


def su22(PS, AS, BS, ES):
    gen = {}
    for a in range(2):
        for b in range(2):
            gen[f"P_{a}{b}"], gen[f"K_{a}{b}"] = PS[(a, b)], PS[(a, b)].T.tocsr()
            gen[f"A_{a}{b}"], gen[f"B_{a}{b}"] = AS[(a, b)], BS[(a, b)]
    gen["E"] = ES
    return gen


Elad = lambda al, Q: Fr(Q*Q, 2) + (al - Fr(1, 2))*Q                # Ladungsenergie eines Quarters (z82)
ALPHA = {"a": Fr(1, 4), "b": Fr(3, 4)}
RSTR = {0: "1", 1: "i", 2: "-1", 3: "-i"}                          # R = i^(T mod 4)

# =================================================================== Teil 1
say("Teil 1: Sind die falschen gekleideten Vertizes V^+ (x) S-quer und V^- (x) S nativ ausgeschlossen?")

# --- Modell A: eine Naht, Ladung (Q_a, Q_b) (x) radiale Standardform
MA, QW = 5, 2
LAB, IXA, PS, AS, BS, ES, SL, SR, LNAT = standardform(MA)
NR = len(LAB)
LAD = [(qa, qb) for qa in range(-QW, QW + 1) for qb in range(-QW, QW + 1)]
LIX = {q: i for i, q in enumerate(LAD)}
NL = len(LAD)


def lschritt(dqa, dqb):
    r, c = [], []
    for (qa, qb), i in LIX.items():
        if (qa + dqa, qb + dqb) in LIX:
            r.append(LIX[(qa + dqa, qb + dqb)])
            c.append(i)
    return sps.csr_matrix((np.ones(len(r)), (r, c)), shape=(NL, NL))


IL, IRd = sps.identity(NL, format="csr"), sps.identity(NR, format="csr")
VPL = {"psi*_a": lschritt(1, 0), "psi_b": lschritt(0, -1)}         # R -> i R
VMI = {"psi_a": lschritt(-1, 0), "psi*_b": lschritt(0, 1)}         # R -> -i R
ST = [(q, s) for q in LAD for s in LAB]
SIX = {x: i for i, x in enumerate(ST)}
T_ = np.array([q[0] - q[1] for q, s in ST])
H2 = np.array([s[0] - s[1] for q, s in ST])
EVA = np.array([0.5*(s[0] + s[1]) for q, s in ST])
NN = T_ - H2
Gd = sps.diags(((-1.0)**T_)*((-1.0)**H2)).tocsr()
Xv = 1j**(NN % 4)
VAK = SIX[((0, 0), (1, 1, 0, 0))]
TYPEN = {}
for nm, V in VPL.items():
    for a in range(2):
        TYPEN[("Psi", nm, a)] = kr(V, SL[a])                        # nativ
        TYPEN[("Phi", nm, a)] = kr(V, SR[a])                        # falsch
for nm, V in VMI.items():
    for a in range(2):
        TYPEN[("Psi-quer", nm, a)] = kr(V, SR[a])                   # nativ
        TYPEN[("Phi-quer", nm, a)] = kr(V, SL[a])                   # falsch


def klassen(O):
    O = sps.coo_matrix(O)
    return {(RSTR[int((T_[r] - T_[c]) % 4)], Fr(int(H2[r] - H2[c]), 2)) for r, c in zip(O.row, O.col)}


def ad_char(D, O):
    """Charakter von O unter Ad_D (D diagonal, als Vektor): +-1, +-i oder None."""
    O = sps.coo_matrix(O)
    ch = {complex(np.round(D[r]/D[c], 12)) for r, c in zip(O.row, O.col)}
    return ch.pop() if len(ch) == 1 else None


def kt(c):
    """Klasse (Delta R, Delta h) als Text, z. B. (i, +1/2)."""
    return f"({c[0]}, {'+' if c[1] > 0 else ''}{c[1]})"


def zs(x):
    """Vierte Einheitswurzel als Text."""
    return next((t for t, v in (("1", 1), ("i", 1j), ("-1", -1), ("-i", -1j)) if x is not None and np.isclose(x, v)), str(x))


KL_A = {}
for (typ, nm, a), O in TYPEN.items():
    KL_A.setdefault(typ, set()).update(klassen(O))
g_gerade = all(nrm(Gd @ O @ Gd - O) == 0.0 for O in TYPEN.values())
nackt_V = kr(VPL["psi*_a"], IRd)
nackt_S = kr(IL, SL[0])
g_nackt = nrm(Gd @ nackt_V @ Gd + nackt_V) == 0.0 and nrm(Gd @ nackt_S @ Gd + nackt_S) == 0.0
X_char = {typ: {ad_char(Xv, O) for (t_, nm, a), O in TYPEN.items() if t_ == typ} for typ in KL_A}
ZAHLEN["Modell A"] = {"Zustände": len(ST), "Klassen (Delta R, Delta h)": {k: sorted(map(kt, v)) for k, v in KL_A.items()},
                      "Ad_X je Typ": {k: [zs(x) for x in v] for k, v in X_char.items()}}
check("[E] Gekleidete Vertizes und ihre Klassen (Modell A: Ladung (Q_a, Q_b) in [-2, 2]^2 (x) radiale Standardform M = 5, "
      f"{len(ST)} Zustände). Die vier Typen Psi_a = V^+ (x) S_a, Psi-quer_b = V^- (x) S-quer_b (nativ, z97), Phi_b = V^+ (x) "
      "S-quer_b und Phi-quer_a = V^- (x) S_a (falsch), je mit beiden Realisierungen V^+ aus psi*_a oder psi_b und V^- aus psi_a "
      f"oder psi*_b (16 Operatoren), sind alle G-gerade (G X G = X exakt: {g_gerade}); nackter Ladungsschritt und nackter Schritt "
      f"sind G-ungerade (Kontrolle: {g_nackt}). Klassen (Delta R, Delta h): Psi {sorted(map(kt, KL_A['Psi']))}, Psi-quer "
      f"{sorted(map(kt, KL_A['Psi-quer']))}, Phi {sorted(map(kt, KL_A['Phi']))}, Phi-quer "
      f"{sorted(map(kt, KL_A['Phi-quer']))}. P_GSO lässt alle vier zu; sie unterscheiden sich genau im Z4-Charakter "
      f"X = R i^(-2h): Ad_X = {[zs(x) for x in X_char['Psi']]} auf Psi, {[zs(x) for x in X_char['Psi-quer']]} auf Psi-quer, "
      f"{[zs(x) for x in X_char['Phi']]} auf Phi, {[zs(x) for x in X_char['Phi-quer']]} auf Phi-quer",
      g_gerade and g_nackt and KL_A["Psi"] == {("i", Fr(1, 2))} and KL_A["Psi-quer"] == {("-i", Fr(-1, 2))}
      and KL_A["Phi"] == {("i", Fr(-1, 2))} and KL_A["Phi-quer"] == {("-i", Fr(1, 2))}
      and X_char["Psi"] == {1} and X_char["Psi-quer"] == {1} and X_char["Phi"] == {-1} and X_char["Phi-quer"] == {-1})

# --- 1.2 CAR: zwei Nähte in der gemeinsamen CAR (z106), Blatttausch Theta
LABc, IXc, PSc, ASc, BSc, ESc, SLc, SRc, _ = standardform(2)
NRc = len(LABc)
NB = 4                                                            # Naht 0: Moden (See, alpha), Naht 1: dito; Quarter a
dimF = 1 << NB
s_ = np.arange(dimF)
occ = (s_[:, None] >> np.arange(NB)[None, :]) & 1
cF = []
for j in range(NB):
    par = occ[:, :j].sum(1) % 2 if j else np.zeros(dimF, dtype=np.int64)
    src = s_[occ[:, j] == 1]
    cF.append(sps.csr_matrix(((1 - 2*par[src]).astype(float), (src ^ (1 << j), src)), shape=(dimF, dimF)))
cdF = [x.T.tocsr() for x in cF]
nF = [sps.diags(occ[:, j].astype(float)).tocsr() for j in range(NB)]
IF, I9 = sps.identity(dimF, format="csr"), sps.identity(NRc, format="csr")


def U_naht(i):
    """Ladungsschritt V^+ der Naht i: Maya-Verschiebung |Q> -> |Q+1> aus den Moden des Quarters (z106), ungerade."""
    lo, hi = 2*i, 2*i + 1
    return (cdF[hi] @ nF[lo] + cdF[lo] @ (IF - nF[hi])).tocsr()


auf = lambda i, X: kr(IF, X, I9) if i == 0 else kr(IF, I9, X)
ferm = lambda X: kr(X, I9, I9)
VT = {}
for i in range(2):
    Ui, Udi = ferm(U_naht(i)), ferm(U_naht(i).T.tocsr())
    for a in range(2):
        VT[("Psi", i, a)] = Ui @ auf(i, SLc[a])
        VT[("Psi-quer", i, a)] = Udi @ auf(i, SRc[a])
        VT[("Phi", i, a)] = Ui @ auf(i, SRc[a])
        VT[("Phi-quer", i, a)] = Udi @ auf(i, SLc[a])
anti, prod = 0.0, []
for k0, Xo in VT.items():
    for k1, Yo in VT.items():
        if k0[1] == 0 and k1[1] == 1:
            anti = max(anti, nrm(Xo @ Yo + Yo @ Xo))
            prod.append(nrm(Xo @ Yo))
perm_c = np.array([IXc[(x[1], x[0], x[3], x[2])] for x in LABc])
Th1 = sps.csr_matrix((np.ones(NRc), (perm_c, np.arange(NRc))), shape=(NRc, NRc))
Th = auf(0, Th1)
moden = [ferm(x) for x in cF + cdF]
th_moden = max(nrm(Th @ x @ Th.T - x) for x in moden)
th_vert = max(max(nrm(Th @ VT[("Psi", 0, a)] @ Th.T - VT[("Phi", 0, a)]),
                  nrm(Th @ VT[("Psi-quer", 0, a)] @ Th.T - VT[("Phi-quer", 0, a)])) for a in range(2))
MVc, MBc = np.array([x[0] for x in LABc]), np.array([x[1] for x in LABc])
Q0 = occ[:, 1] - (1 - occ[:, 0])
G0 = kr(sps.diags((-1.0)**Q0), sps.diags(np.where((MVc - MBc) % 2 == 0, 1.0, -1.0)), I9)
R0 = kr(sps.diags(1j**(Q0 % 4)), I9, I9)
h0 = auf(0, sps.diags((MVc - MBc)/2))
th_G, th_R, th_E = nrm(Th @ G0 @ Th.T - G0), nrm(Th @ R0 @ Th.T - R0), nrm(Th @ auf(0, ESc) @ Th.T - auf(0, ESc))
th_h, th_h_nt = nrm(Th @ h0 @ Th.T + h0), nrm(Th @ h0 @ Th.T - h0)     # Theta h Theta = -h, also Theta != 1 auf h
genc = su22(PSc, ASc, BSc, ESc)
bild = {}
for nm, Xg in genc.items():
    Yg = Th1 @ Xg @ Th1.T
    bild[nm] = [n2 for n2, X2 in genc.items() if nrm(Yg - X2) < 1e-12]
th_su = all(len(v) == 1 for v in bild.values()) and sorted(v[0] for v in bild.values()) == sorted(genc)
S0n = auf(0, SLc[0])
kom_nackt = max(nrm(S0n @ VT[("Psi", 1, a)] - VT[("Psi", 1, a)] @ S0n) for a in range(2))
akom_nackt = max(nrm(S0n @ VT[("Psi", 1, a)] + VT[("Psi", 1, a)] @ S0n) for a in range(2))
ZAHLEN["CAR zwei Nähte"] = {"Dimension": dimF*NRc*NRc, "max |{X^(0), Y^(1)}|, 64 Paare": anti,
                            "min max|X^(0) Y^(1)|": min(prod), "Theta: Moden, Vertizes, G, R, E, h": [th_moden, th_vert, th_G,
                                                                                                  th_R, th_E, th_h],
                            "Theta auf su(2,2)": bild, "nackter Schritt: [ , ] und { , } mit Psi^(1)": [kom_nackt, akom_nackt]}
check("[E, negativ] CAR: Zwei Nähte in der gemeinsamen CAR (z106: Quarter a mit See- und alpha-Mode je Naht, Ladungsschritt = "
      f"Maya-Verschiebung, Standardform M = 2 je Naht, Dimension {dimF*NRc*NRc}). Alle 64 Antikommutatoren {{X^(0), Y^(1)}} "
      f"der acht Vertizes (nativ und falsch, beide Spinorindizes) verschwinden (max {anti:.0e}, Produkte selbst vom Betrag "
      f"{min(prod):.1f}): Die CAR macht jeden gekleideten Vertex ungerade, gleich welches Blatt er hebt. Der Blatttausch Theta "
      "(m <-> mquer der Naht 0) lässt alle 8 Moden fest (Abweichung "
      f"{th_moden:.0e}), bildet Psi auf Phi und Psi-quer auf Phi-quer ab ({th_vert:.0e}), vertauscht mit G, R und E "
      f"({th_G:.0e}, {th_R:.0e}, {th_E:.0e}) und permutiert die 17 Erzeuger der radialen su(2,2) (P_ab -> P_ba, A <-> B, E "
      f"fest: {th_su}). Jede Bedingung aus CAR, G, R, E und su(2,2) ist also blind für die Paarung. Kontrollen: Theta ist "
      f"nicht trivial (Theta h Theta + h = {th_h:.0e}, Theta h Theta - h = {th_h_nt:.1f}), und die CAR sieht die Statistik: "
      f"Der nackte Schritt vertauscht mit den Vertizes der anderen Naht ([ , ] = {kom_nackt:.0e}, {{ , }} = "
      f"{akom_nackt:.1f})",
      anti == 0.0 and min(prod) > 0.5 and th_moden == 0.0 and th_vert == 0.0 and th_G == 0.0 and th_R == 0.0
      and th_E == 0.0 and th_su and th_h == 0.0 and th_h_nt > 0.5 and kom_nackt == 0.0 and akom_nackt > 1.0)

# --- 1.3 Modeninhalt: Maya-Fockraum des Originals (beide Quarter)
KMIN, KMAX = -9, 8
vakuum = lambda Q: frozenset(range(KMIN, Q))


def L0_maya(s, al):
    return sum((k + al) for k in s if k >= 0) - sum((k + al) for k in range(KMIN, 0) if k not in s)


def ladung(s):
    return sum(1 for k in s if k >= 0) - sum(1 for k in range(KMIN, 0) if k not in s)


def paarstufen(s):
    """(Teilchenstufe n_1, Lochstufe n_2) eines Einpaarzustands relativ zum Vakuum seiner Ladung (z83: Schrittart)."""
    Q = ladung(s)
    v = vakuum(Q)
    t, l = sorted(k for k in s if k not in v), sorted(k for k in v if k not in s)
    return (t[0] - Q, Q - 1 - l[0]) if len(t) == 1 and len(l) == 1 else None


moden_tab, ok_pos, n_paare = [], True, 0
for q, al in ALPHA.items():
    for Q in range(-2, 3):
        s0 = vakuum(Q)
        for art, rr, dq in (("psi*", range(Q, Q + 5), 1), ("psi", range(Q - 5, Q), -1)):
            for r in rr:
                s1 = (s0 | {r}) if art == "psi*" else (s0 - {r})
                dL = L0_maya(s1, al) - L0_maya(s0, al)
                dosc = dL - (Elad(al, Q + dq) - Elad(al, Q))
                ok_pos &= dosc >= 0 and ladung(s1) == Q + dq
                n_paare += 1
                dT = dq if q == "a" else -dq
                moden_tab.append({"Quarter": q, "Mode": art, "Q": Q, "r": r, "Delta L_0": dL, "Delta 2h": dosc,
                                  "Delta R": RSTR[dT % 4], "G-gerade": dosc % 2 == 1, "Stufen": paarstufen(s1)})
tief = {}
for t in moden_tab:
    if t["Q"] == 0 and t["Delta 2h"] == 1:
        tief.setdefault(t["Delta R"], []).append((f"{t['Mode']}_{t['Quarter']}", t["Delta L_0"]))
TIEF_TXT = {k: [f"{a} bei L_0 = {b}" for a, b in sorted(v)] for k, v in tief.items()}
art3 = {(t["Mode"], t["Quarter"]): (t["Delta R"], t["Stufen"]) for t in moden_tab if t["Q"] == 0 and t["Delta 2h"] == 3}
ART3_TXT = {f"{a}_{b}": f"R -> {v[0]} R, (n_1, n_2) = {v[1]}" for (a, b), v in art3.items()}


def partitionen(n, mx=None):
    mx = n if mx is None else mx
    if n == 0:
        yield ()
        return
    for k in range(min(n, mx), 0, -1):
        for rest in partitionen(n - k, k):
            yield (k,) + rest


def maya(Q, lam):
    s = set(range(KMIN, Q - len(lam)))
    for i, l in enumerate(lam):
        s.add(l - (i + 1) + Q)
    return frozenset(s)


ok_L = all(L0_maya(maya(Q, lam), al) == Elad(al, Q) + n and ladung(maya(Q, lam)) == Q
           for al in ALPHA.values() for Q in range(-2, 3) for n in range(5) for lam in partitionen(n))
npart = [sum(1 for _ in partitionen(n)) for n in range(7)]
spek = []
for Qa in range(-3, 4):
    for Qb in range(-3, 4):
        for na in range(7):
            for nb in range(7 - na):
                L = Elad(ALPHA["a"], Qa) + na + Elad(ALPHA["b"], Qb) + nb
                spek.append((RSTR[(Qa - Qb) % 4], na + nb, L, (-1)**(Qa + Qb + na + nb), npart[na]*npart[nb]))
tiefst = {}
for R in ("i", "-i"):
    ger = [(L, h2, m) for (RR, h2, L, g, m) in spek if RR == R and g == 1]
    Lmin = min(L for L, _, _ in ger)
    tiefst[R] = {"L_0": Lmin, "2h": sorted({h2 for L, h2, m in ger if L == Lmin}),
                 "Zustände": sum(m for L, h2, m in ger if L == Lmin)}
neg_h = sum(1 for (RR, h2, L, g, m) in spek if RR in ("i", "-i") and h2 < 0)
ZAHLEN["Modeninhalt"] = {"(Zustand, Mode)-Paare": n_paare, "tiefste G-gerade Einzelmoden auf |0, 0> (Delta 2h = 1)":
                         {k: [(a, str(b)) for a, b in v] for k, v in tief.items()},
                         "Schrittart bei Delta 2h = 3": ART3_TXT,
                         "tiefste G-gerade Fermionstufe des Originals": {k: {kk: str(vv) for kk, vv in v.items()}
                                                                         for k, v in tiefst.items()}}
check("[E, negativ] Modeninhalt (Maya-Fockraum, beide Quarter, alpha = 1/4, 3/4; L_0 = E_lad + |lambda| auf allen "
      f"Maya-Zuständen bestätigt: {ok_L}): Jede native Einzelmode psi*_r, psi_r auf jedem Ladungsvakuum hebt nur L_0^osc des "
      f"Originals (Delta 2h >= 0 auf {n_paare} Paaren: {ok_pos}); das Original trägt also für Teilchen und Löcher sein eigenes "
      f"Blatt. G-gerade mit Delta 2h = 1 sind auf |0, 0>: R -> i R durch {TIEF_TXT['i']} und R -> -i R durch {TIEF_TXT['-i']}. Der "
      "falsche Vertex V^- (x) S (Klasse (-i, +1/2)) ist damit eine native Mode des Originals (Loch in a, Teilchen in b), der "
      "native Psi-quer (Klasse (-i, -1/2)) dagegen keine. Bei Delta 2h = 3 liegen Teilchen rein in der Schrittart 0 "
      f"(Teilchen hoch) und Löcher rein in 1 (Loch runter): {ART3_TXT}; jedes Dublett ist nativ realisiert, das R = i-Dublett "
      "durch (psi*_a, psi_b), das R = -i-Dublett durch (psi*_b, psi_a). Volles Spektrum (|Q| <= 3, Oszillatorstufe <= 6, "
      f"Partner im Vakuum): tiefste G-gerade Fermionstufe R = i bei L_0 = {tiefst['i']['L_0']} mit 2h = {tiefst['i']['2h']} "
      f"({tiefst['i']['Zustände']} Zustände), R = -i bei {tiefst['-i']['L_0']} mit 2h = {tiefst['-i']['2h']} "
      f"({tiefst['-i']['Zustände']}); Zustände mit R = +-i und h < 0: {neg_h}. Der Modeninhalt schließt den falschen Vertex "
      "nicht aus, er realisiert ihn",
      ok_L and ok_pos and sorted(tief["i"]) == [("psi*_a", Fr(5, 4)), ("psi_b", Fr(5, 4))]
      and sorted(tief["-i"]) == [("psi*_b", Fr(7, 4)), ("psi_a", Fr(7, 4))]
      and art3[("psi*", "a")] == ("i", (2, 0)) and art3[("psi", "b")] == ("i", (0, 2))
      and art3[("psi*", "b")] == ("-i", (2, 0)) and art3[("psi", "a")] == ("-i", (0, 2))
      and tiefst["i"] == {"L_0": Fr(5, 4), "2h": [1], "Zustände": 4}
      and tiefst["-i"] == {"L_0": Fr(7, 4), "2h": [1], "Zustände": 4} and neg_h == 0)

# --- Modell B: Original und Partner mit Ladung (Quarter a je Seite), Standardform, KMS-Vektor
MB, QB = 4, 2
LABb, IXb, PSb, ASb, BSb, ESb, SLb, SRb, _ = standardform(MB)
NRb = len(LABb)
LAD2 = [(qo, qp) for qo in range(-QB, QB + 1) for qp in range(-QB, QB + 1)]
L2 = {q: i for i, q in enumerate(LAD2)}
NL2 = len(LAD2)


def sch2(do, dp):
    r, c = [], []
    for (qo, qp), i in L2.items():
        if (qo + do, qp + dp) in L2:
            r.append(L2[(qo + do, qp + dp)])
            c.append(i)
    return sps.csr_matrix((np.ones(len(r)), (r, c)), shape=(NL2, NL2))


CH = {("o", 1): sch2(1, 0), ("o", -1): sch2(-1, 0), ("p", 1): sch2(0, 1), ("p", -1): sch2(0, -1)}   # c = +1: R_Seite -> i R
BLATT = {"o": SLb, "p": SRb}
STb = [(q, s) for q in LAD2 for s in LABb]
SIXb = {x: i for i, x in enumerate(STb)}
Tb = np.array([q[0] - q[1] for q, s in STb])                     # R_phys = R R-quer^-1 = i^(Q_o - Q_p)
H2b = np.array([s[0] - s[1] for q, s in STb])
EL = {(sc, c, sl): kr(CH[(sc, c)], BLATT[sl][0]) for sc in "op" for c in (1, -1) for sl in "op"}


def klasse_b(O):
    O = sps.coo_matrix(O)
    k = {(RSTR[int((Tb[r] - Tb[c]) % 4)], Fr(int(H2b[r] - H2b[c]), 2)) for r, c in zip(O.row, O.col)}
    return k.pop() if len(k) == 1 else None


KLB = {k: klasse_b(O) for k, O in EL.items()}
NATIV = {("i", Fr(1, 2)), ("-i", Fr(-1, 2))}
FALSCH = {("i", Fr(-1, 2)), ("-i", Fr(1, 2))}
perm_b = np.array([SIXb[((q[1], q[0]), (s[1], s[0], s[3], s[2]))] for q, s in STb])
PiB = sps.csr_matrix((np.ones(len(STb)), (perm_b, np.arange(len(STb)))), shape=(len(STb), len(STb)))
Jop = lambda O: (PiB @ O.conj() @ PiB.T).tocsr()                  # J O J^-1, J antiunitär (Tausch + Konjugation)
J_bild = {}
for k, O in EL.items():
    JO = Jop(O)
    J_bild[k] = [k2 for k2, O2 in EL.items() if nrm(JO - O2) < 1e-12]
J_klasse = {KLB[k]: KLB[v[0]] for k, v in J_bild.items() if len(v) == 1}
J_ok = all(len(v) == 1 for v in J_bild.values()) and all(v[0] == (("p" if k[0] == "o" else "o"), k[1],
                                                                  ("p" if k[2] == "o" else "o")) for k, v in J_bild.items())
regeln_kl = []
for si in itertools.product([(), ("+",), ("-",), ("+", "-")], repeat=2):
    regel = {("i", Fr(1, 2) if x == "+" else Fr(-1, 2)) for x in si[0]} | {("-i", Fr(1, 2) if x == "+" else Fr(-1, 2))
                                                                           for x in si[1]}
    regeln_kl.append(frozenset(regel))
j_inv = [r for r in regeln_kl if {J_klasse[c] for c in r} == set(r)]
j_chiral = [r for r in j_inv if len(r) == 2 and {c[0] for c in r} == {"i", "-i"}]
halb = {("i", Fr(1, 2))}
# KMS-Vektor (Purifikation): Ladung sqrt(e^(-beta E(Q))) |Q>|Q>, radial sum q^(m-1) |m,k>|m,k>-quer (z83), q = e^(-beta/2)
qq = np.exp(-BETA/2)
om_ch = np.zeros(NL2)
for Q in range(-QB, QB + 1):
    om_ch[L2[(Q, Q)]] = np.sqrt(np.exp(-BETA*float(Elad(ALPHA["a"], Q))))
om_ch /= np.linalg.norm(om_ch)
om_rad = np.zeros(NRb)
for m in range(1, MB + 1):
    for k in range(m):
        om_rad[IXb[(m, m, k, k)]] = qq**(m - 1)
om_rad /= np.linalg.norm(om_rad)
OM = np.kron(om_ch, om_rad)
J_OM = float(np.abs(PiB @ OM.conj() - OM).max())
innen = np.array([s[0] < MB and s[1] < MB for s in LABb])
tom_rad = max(float(np.abs(SLb[a].T @ om_rad - qq*(SRb[a] @ om_rad))[innen].max()) for a in range(2))
v1, v2 = CH[("o", -1)] @ om_ch, CH[("p", 1)] @ om_ch
gleich_traeger = set(np.nonzero(v1)[0]) == set(np.nonzero(v2)[0])
tom_lad = max(abs(v1[i]/v2[i] - np.exp(-BETA*float(Elad(ALPHA["a"], qp) - Elad(ALPHA["a"], qp - 1))/2))
              / np.exp(-BETA*float(Elad(ALPHA["a"], qp) - Elad(ALPHA["a"], qp - 1))/2)
              for (qo, qp), i in L2.items() if v2[i] != 0)


def klasse_vektor(v):
    tr = np.nonzero(np.abs(v) > 1e-14*np.abs(v).max())[0]
    return {(RSTR[int(Tb[i] % 4)], Fr(int(H2b[i]), 2)) for i in tr}


S_nat = klasse_vektor(EL[("o", 1, "o")].T @ OM)                 # A* Omega für A = V^+_o (x) S (nativ, seitenlokal)
S_falsch = klasse_vektor(EL[("o", -1, "o")].T @ OM)             # A* Omega für A = V^-_o (x) S (falsch, seitenlokal)
ZAHLEN["Tomita"] = {"J auf den Elementarvertizes (Ladungsseite, c, Blattseite)": {f"({k[0]}, {'+i' if k[1] == 1 else '-i'}, {k[2]})": f"({v[0][0]}, {'+i' if v[0][1] == 1 else '-i'}, {v[0][2]})" for k, v in J_bild.items()},
                    "J auf den Klassen": {kt(k): kt(v) for k, v in J_klasse.items()},
                    "J-invariante Klassenregeln": len(j_inv), "davon chiral": [sorted(map(kt, r)) for r in j_chiral],
                    "J Omega - Omega": J_OM, "a_c Omega - q b_c^+ Omega (innen)": tom_rad,
                    "V^-_o Omega / V^+_p Omega gegen Delta-Gewicht (relativ)": tom_lad,
                    "Klasse von A* Omega": {"nativ (o,+i,o)": sorted(map(kt, S_nat)),
                                            "falsch (o,-i,o)": sorted(map(kt, S_falsch))}}
check("[E, negativ] Tomita (Modell B: Ladung Q_o, Q_p in [-2, 2] je Seite, Standardform M = 4, "
      f"{len(STb)} Zustände, R_phys = R R-quer^-1 = i^(Q_o - Q_p), J antiunitär mit Seitentausch und m <-> mquer, J R_phys J = "
      "R_phys). J bildet jeden der acht Elementarvertizes (Ladungsseite, Charakter c, Blattseite) auf (Gegenseite, c, "
      f"Gegenseite) ab ({J_ok}) und auf den Klassen {{nativ}} -> {{nativ}}, {{falsch}} -> {{falsch}}: "
      f"{ {kt(k): kt(v) for k, v in J_klasse.items()} }. Von den 16 Klassenregeln sind {len(j_inv)} J-invariant, chiral "
      f"davon genau {len(j_chiral)} (nativ und Spiegel); eine halbe Paarung {{(i, +1/2)}} ist nicht J-abgeschlossen (Kontrolle: "
      f"{'nein' if {J_klasse[c] for c in halb} != halb else 'ja'}). KMS-Vektor: J Omega = Omega ({J_OM:.0e}); Tomita-"
      f"Relationen a_c Omega = q b_c^+ Omega (radial, q = e^(-2 pi), {tom_rad:.0e}) und V^-_o Omega = Delta-Gewicht mal "
      f"V^+_p Omega (gleicher Träger: {gleich_traeger}, relativ {tom_lad:.0e}): Das Loch des Originals ist auf Omega das "
      "Teilchen des Partners. S = J Delta^(1/2) bildet A Omega auf A* Omega ab: für den nativen seitenlokalen Vertex V^+_o (x) "
      f"S liegt A* Omega in {sorted(map(kt, S_nat))} (nativ), für den falschen V^-_o (x) S in {sorted(map(kt, S_falsch))} "
      "(falsch). J und S erhalten beide Paarungen; sie wählen nicht",
      J_ok and set(J_klasse) == NATIV | FALSCH and all((J_klasse[c] in NATIV) == (c in NATIV) for c in J_klasse)
      and len(j_inv) == 4 and len(j_chiral) == 2 and {J_klasse[c] for c in halb} != halb and J_OM == 0.0
      and tom_rad < 1e-15 and gleich_traeger and tom_lad < 1e-10 and S_nat == {("-i", Fr(-1, 2))}
      and S_falsch == {("i", Fr(-1, 2))})

# --- 1.5 Energiepositivität und Lokalität
P0 = (ES + (PS[(0, 0)] + PS[(1, 1)] + PS[(0, 0)].T + PS[(1, 1)].T)/2).toarray()      # Minkowskizeit (z83)
HV = np.array([s[0] - s[1] for s in LAB])
p0_min = {hh: float(np.linalg.eigvalsh(P0[np.ix_(HV == hh, HV == hh)]).min()) for hh in (1, -1)}
e_hub = {typ: sorted({float(EVA[r] - EVA[c]) for O in [TYPEN[k] for k in TYPEN if k[0] == typ]
                      for r, c in zip(*sps.coo_matrix(O).nonzero())}) for typ in KL_A}
KDIAG = {}
for (q, s), i in SIXb.items():
    KDIAG[i] = (Elad(ALPHA["a"], q[0]) - Elad(ALPHA["a"], q[1])) + (s[0] - s[1])     # K = L_0 - L-quer_0 (Ladung + Blatt)
start_b = SIXb[((0, 0), (1, 1, 0, 0))]
dK = {}
for k, O in EL.items():
    col = sps.csc_matrix(O).getcol(start_b)
    ziel = col.indices[0]
    dK[k] = KDIAG[ziel] - KDIAG[start_b]
modpos = sorted({KLB[k] for k in EL if dK[k] > 0}, key=str)
ek = lambda k: f"({k[0]}, {'+i' if k[1] == 1 else '-i'}, {k[2]})"
DK_TXT = {ek(k): str(v) for k, v in dK.items()}
signum = all((dK[k] > 0) == (KLB[k][1] > 0) for k in EL)
lokal = {c: [k for k in EL if KLB[k] == c and k[0] == k[2]] for c in sorted(NATIV | FALSCH, key=kt)}   # sortiert: reproduzierbare Ausgabe
lad_lokal = {kt(c): str(Elad(ALPHA["a"], v[0][1])) for c, v in lokal.items()}
lad_gem = {kt(c): str(Elad(ALPHA["a"], k[1])) for c in sorted(NATIV | FALSCH, key=kt) for k in EL if KLB[k] == c and k[0] != k[2]}
boltz = math.exp(-BETA*float(Fr(7, 4) - Fr(5, 4)))
ZAHLEN["Energie und Lokalität"] = {"min P_0 auf h = +1/2 / -1/2 (M = 5)": p0_min, "Delta E je Typ": e_hub,
                                   "Modulgewicht Delta K (Ladung + Blatt)": DK_TXT,
                                   "Modulpositiv": [kt(c) for c in modpos], "seitenlokale Realisierung": {kt(c): str(v) for c, v in lokal.items()},
                                   "Ladungsenergie seitenlokal": lad_lokal, "Ladungsenergie gemischt": lad_gem,
                                   "Boltzmannverhältnis falsch/nativ (beta = 4 pi)": boltz}
check("[E, negativ] Energiepositivität und Lokalität. (a) Alle vier Typen heben E um genau 1/2 "
      f"({ {k: v for k, v in e_hub.items()} }); die Minkowskizeit P_0 = E + (P^0 + K^0)/2 hat auf den Sektoren h = +1/2 und "
      f"-1/2 dasselbe positive Spektrum (min {p0_min[1]:.6f} und {p0_min[-1]:.6f}, M = 5): Positivität ist helizitätsblind "
      "(P_0 ist Theta-invariant). (b) Modulgewicht Delta K (K = L_0 - L-quer_0 auf Ladung und Blatt, Modell B): "
      f"{DK_TXT}; sgn Delta K = sgn h für alle acht ({signum}), also wählt "
      f"Modulpositivität die Klassen {[kt(c) for c in modpos]}: beide R mit der Helizität des Originals (universelle "
      "Gerade) und verwirft den nativen Psi-quer, nicht den falschen Vertex. (c) Seitenlokalität (Element von M oder M'): "
      f"jede Klasse hat genau eine seitenlokale Realisierung ({ {kt(c): len(v) for c, v in lokal.items()} }). (d) Was die "
      f"Paarungen trennt, ist allein die Ladungsenergie der seitenlokalen Realisierung: nativ {lad_lokal['(i, +1/2)']} "
      f"(tiefste G-gerade Fermionstufe je Seite, L_0 = 5/4), falsch {lad_lokal['(i, -1/2)']} (L_0 = 7/4), "
      f"Boltzmannverhältnis e^(-2 pi) = {boltz:.3e}; gemischt realisiert kehrt sich das um ({lad_gem}). Energie ordnet, "
      "schließt aber nichts aus, und auch das nur zusammen mit Seitenlokalität",
      all(v == [0.5] for v in e_hub.values()) and abs(p0_min[1] - p0_min[-1]) < 1e-12 and p0_min[1] > 0 and signum
      and set(modpos) == {("i", Fr(1, 2)), ("-i", Fr(1, 2))} and all(len(v) == 1 for v in lokal.values())
      and all(lad_lokal[kt(c)] == "1/4" for c in NATIV) and all(lad_lokal[kt(c)] == "3/4" for c in FALSCH)
      and all(lad_gem[kt(c)] == "3/4" for c in NATIV) and all(lad_gem[kt(c)] == "1/4" for c in FALSCH)
      and abs(boltz - math.exp(-2*math.pi)) < 1e-15)

# --- 1.6 Was die native Paarung ist
# I_T = Ad Z_SU: Zentrum von SU(2,2) aus den nativen Erzeugern, Z_SU = exp(i pi (E + 2 J3_R)) (z97), Homogenitätscharakter
ZSU_rad = np.exp(1j*np.pi*(ES.diagonal() + (BS[(0, 0)] - BS[(1, 1)]).diagonal()))
zsu_formel = float(np.abs(ZSU_rad - 1j**(np.array([s[0] - s[1] - 2 for s in LAB]) % 4)).max())
zsu_zentral = max(nrm(sps.diags(ZSU_rad) @ X - X @ sps.diags(ZSU_rad)) for X in su22(PS, AS, BS, ES).values())
UTH = np.tile(ZSU_rad, NL)                                         # Z_SU auf Ladung (x) Standardform
RV = 1j**(T_ % 4)                                                 # I_F = Ad_R (z116: R = i Gamma_11 = i auf der 16)
phasen = {}
for (typ, nm, a), O in TYPEN.items():
    phasen.setdefault(typ, set()).add((ad_char(UTH, O), ad_char(RV, O)))
balanciert = {typ: {pT == pF for pT, pF in v} for typ, v in phasen.items()}
N_erh = {typ: {int(NN[r] - NN[c]) for (t_, nm, a), O in TYPEN.items() if t_ == typ
               for r, c in zip(*sps.coo_matrix(O).nonzero())} for typ in KL_A}
gen_A = su22(PS, AS, BS, ES)
N_su = all(np.all(NN[sps.coo_matrix(kr(IL, X)).row] == NN[sps.coo_matrix(kr(IL, X)).col]) for X in gen_A.values())
W_rand = kr(VPL["psi*_a"] @ VPL["psi_b"], IRd)                      # Randpaar psi*_a psi_b: R -> -R, h fest (z82: L_0 = 1/2)
D_nat = kr(IL, LNAT @ LNAT)                                        # (L^nat_-1)^2: h + 1, R fest (z83)
ung = {nm: (nrm(Gd @ O @ Gd - O) == 0.0, ad_char(Xv, O)) for nm, O in (("Randpaar psi*_a psi_b", W_rand),
                                                                       ("(L^nat_-1)^2", D_nat))}
z97 = json.loads((HERE / "z97_ergebnis.json").read_text())["zahlen"]["E8_kleber"]["Regel aus Kleber + P_A3: erlaubte h je R"]
tab_if, tab_a3 = {}, {}
for k in range(4):
    for zh in range(-3, 4):
        if (zh - k) % 2:
            continue                                              # P_GSO
        X_ = 1j**((k - zh) % 4)
        Y_ = 1j**((k + zh) % 4)
        if np.isclose(X_, 1):
            tab_if.setdefault(RSTR[k], []).append(str(Fr(zh, 2)))
        if np.isclose(Y_, -1):
            tab_a3.setdefault(RSTR[k], []).append(str(Fr(zh, 2)))
gleich_ferm = all(sorted(tab_if[R]) == sorted(tab_a3[R]) for R in ("i", "-i"))
ungleich_bos = all(set(tab_if[R]).isdisjoint(tab_a3[R]) for R in ("1", "-1"))
a3_z97 = all(sorted(tab_a3[R]) == sorted(z97[R]) for R in z97)
PH_TXT = {k: sorted((zs(a), zs(b)) for a, b in v) for k, v in phasen.items()}
UNG_TXT = {k: (v[0], zs(v[1])) for k, v in ung.items()}
ZAHLEN["native Paarung"] = {"(Twistorphase, Uhrphase) je Typ": PH_TXT,
                            "Delta N je Typ": {k: sorted(v) for k, v in N_erh.items()},
                            "X-ungerade native Bosonen (G-gerade, Ad_X)": UNG_TXT,
                            "P_IF (X = 1) erlaubte h je R": tab_if, "P_A3 (aus Kleber) erlaubte h je R": tab_a3}
check("[E] Was die native Paarung ist. Als Twistorphase dient das Zentrum von SU(2,2) aus den nativen Erzeugern, "
      f"Z_SU = exp(i pi (E + 2 J3_R)) = i^(2h - 2) (Homogenitätscharakter, z97/z114; Abweichung {zsu_formel:.0e}, zentral "
      f"{zsu_zentral:.0e}), I_T = Ad Z_SU, gegen die Uhr I_F = Ad_R (z116). Auf den gekleideten Vertizes ist (I_T, I_F) = "
      f"{PH_TXT}, gleich genau auf den nativen ({balanciert}). Die native Paarung ist also eine der beiden Verklebungen "
      "I_T = +-I_F aus z116, die mit Ad Z_SU = Ad R (R Z_SU^-1 erhalten, auf dem Vakuumsektor R = -Z_SU), gleichwertig X = "
      "R i^(-2h) = 1 auf jedem Einzelquant und, je Seite, die eta_+-Kleidung (das Blatt einer Seite mit ihren "
      "Ad_R = +i-Moden, Teil 1.3). Ihre ganzzahlige Hebung N = T - 2h ändert sich unter Psi, Psi-quer um "
      f"{N_erh['Psi']}, {N_erh['Psi-quer']}, unter Phi, Phi-quer um {N_erh['Phi']}, {N_erh['Phi-quer']}, unter su(2,2) gar "
      f"nicht ({N_su}); G = (-1)^N ist P_GSO. X ist aber keine Symmetrie der nativen geraden Algebra: das Randpaar psi*_a "
      f"psi_b (R -> -R, z82, L_0 = 1/2) und (L^nat_-1)^2 (z83) sind G-gerade und X-ungerade ({UNG_TXT}). Als Spektralregel fällt "
      f"P_IF auf den Fermionen mit P_A3 zusammen ({ {R: tab_if[R] for R in ('i', '-i')} }; P_A3 aus z97 reproduziert: "
      f"{a3_z97}) und ist auf den Bosonen genau entgegengesetzt (P_IF {[(R, tab_if[R]) for R in ('1', '-1')]}, P_A3 "
      f"{[(R, tab_a3[R]) for R in ('1', '-1')]}): P_IF behält die neutrale Skalarleiter und verliert die R-neutralen "
      "Vektoren, P_A3 umgekehrt (z102)",
      balanciert == {"Psi": {True}, "Psi-quer": {True}, "Phi": {False}, "Phi-quer": {False}} and zsu_formel < 1e-12
      and zsu_zentral < 1e-12
      and N_erh["Psi"] == {0} and N_erh["Psi-quer"] == {0} and N_erh["Phi"] == {2} and N_erh["Phi-quer"] == {-2}
      and N_su and all(v == (True, -1) for v in ung.values()) and gleich_ferm and ungleich_bos and a3_z97
      and "0" in tab_if["1"] and "1" not in tab_if["1"] and "1" in tab_a3["1"] and "0" not in tab_a3["1"])
befund("Teil 1: Die falschen Vertizes sind nicht nativ ausgeschlossen. CAR, G, R, E, P_0 und su(2,2) sind blind (der "
       "Blatttausch Theta fixiert sie und tauscht die Paarungen), J und Tomita-S erhalten beide Paarungen, Modulpositivität "
       "wählt die universelle Gerade, und der Modeninhalt realisiert den falschen Vertex V^- (x) S sogar als native Mode des "
       "Originals (psi_a, psi*_b bei L_0 = 7/4). Was nativ trennt, ist allein die Ladungsenergie der seitenlokalen "
       "Realisierung (nativ 5/4, falsch 7/4). Die native Paarung ist die Prämisse P_IF: Verklebung der Twistorphase mit der Uhr, "
       "Ad Z_SU = Ad R (eine der beiden Verklebungen I_T = +-I_F aus z116), gleich X = R i^(-2h) = 1 auf den Einzelquanten, "
       "gleich eta_+-Kleidung je Seite", True)

# =================================================================== Teil 2
say("\nTeil 2: Wählt die native Paarung P_RJ? z114s Linienregeln mit der Paarung als Bedingung")


def cech_O(n, L=None):
    """h^0, h^1 von O(n) auf CP^1 (Cech, Laurent-Monome, wie z114)."""
    L = abs(n) + 3 if L is None else L
    C1 = list(range(-L, L + 1))
    ix = {p: i for i, p in enumerate(C1)}
    sp_ = [(ix[p0], 1.0) for p0 in C1 if n - p0 >= 0] + [(ix[p0], -1.0) for p0 in C1 if p0 >= 0]
    D = np.zeros((len(C1), len(sp_)))
    for j, (i, v) in enumerate(sp_):
        D[i, j] = v
    r = int(np.linalg.matrix_rank(D))
    return D.shape[1] - r, D.shape[0] - r


def coh(F, n):
    h = [cech_O(a + n) for a in F]
    return sum(x for x, _ in h), sum(y for _, y in h)


def kuenneth(F, n):
    f0, f1 = coh(F, 0)
    l0, l1 = cech_O(n)
    return f0*l0, f0*l1 + f1*l0, f1*l1


def restriktion(n, x, Kmax=5):
    """Rang von H^1(PT', O(n)) -> H^1(L_x, O(n)) auf Cech-Klassen mit omega-Grad <= Kmax (wie z114)."""
    zeilen = list(range(n + 1, 0))
    zi = {p: i for i, p in enumerate(zeilen)}
    spalten = []
    for k in range(Kmax + 1):
        s2 = n - k
        for a0 in range(k + 1):
            for p0 in range(s2 + 1, 0):
                spalten.append((a0, k - a0, p0, s2 - p0))
    R = np.zeros((len(zeilen), len(spalten)), complex)
    for c, (a0, a1, p0, p1) in enumerate(spalten):
        for i in range(a0 + 1):
            for j in range(a1 + 1):
                e0, e1 = p0 + i + j, p1 + (a0 - i) + (a1 - j)
                if e0 <= -1 and e1 <= -1:
                    R[zi[e0], c] += (math.comb(a0, i)*x[0, 0]**i*x[0, 1]**(a0 - i)
                                     * math.comb(a1, j)*x[1, 0]**j*x[1, 1]**(a1 - j))
    if R.size == 0:
        return 0
    sv = np.linalg.svd(R, compute_uv=False)
    return int(np.sum(sv > 1e-9*max(1.0, sv[0])))


ind = lambda hh: hh[0] - hh[1]
n_z97 = lambda h: int(2*h - 2)
n_pen = lambda h: int(-2*h - 2)
EM = [-1, 0, 0, 0]
xs = [np.zeros((2, 2), complex), rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2))]
HEL3 = [Fr(k, 2) for k in range(-3, 4)]
SICHT = {"p": {h: max(restriktion(n_z97(h), x) for x in xs) > 0 for h in HEL3},
         "o": {h: max(restriktion(n_z97(-h), x) for x in xs) > 0 for h in HEL3}}
chiF = {"E_-": ind(coh(EM, 0)), "E_-^v (x) K": ind(coh([-1, -2, -2, -2], 0)), "4 O": ind(coh([0]*4, 0)),
        "Spiegelfamilie O(-4) + 3 O": ind(coh([-4, 0, 0, 0], 0))}
GSO = lambda R, h: (h.denominator == 1) == (R in (1, -1))
MULT = {1: 1, 1j: chiF["E_-"], -1: coh([-1, -1, -1, 0, 0, 0], 0)[0], -1j: -chiF["E_-^v (x) K"]}
TEILM = [frozenset(s) for s in ((), ("o",), ("p",), ("o", "p"))]
Jseite = lambda S: frozenset({"o": "p", "p": "o"}[x] for x in S)
RS = (1, 1j, -1, -1j)
RNAME = {1: "R = 1", 1j: "R = i", -1: "R = -1", -1j: "R = -i"}
Rinv = {1: 1, 1j: -1j, -1: -1, -1j: 1j}


def bewerte(regel, m16=None):
    m16 = MULT[1j] if m16 is None else m16
    wz = lambda R, h, m_: m_ if GSO(R, h) and any(SICHT[s][h] for s in regel[R]) else 0
    N16 = wz(1j, Fr(1, 2), m16) - wz(1j, Fr(-1, 2), m16)
    N16q = wz(-1j, Fr(-1, 2), m16) - wz(-1j, Fr(1, 2), m16)
    spiegel = min(wz(1j, Fr(1, 2), m16), wz(1j, Fr(-1, 2), m16))
    return {"N16": N16, "N16q": N16q, "chiral": N16 != 0 and N16 == N16q and spiegel == 0,
            "skalar": wz(1, Fr(0), 1) > 0, "vektor": wz(1, Fr(1), 1) > 0 and wz(1, Fr(-1), 1) > 0,
            "alle": all(len(regel[R]) > 0 for R in RS), "J": all(regel[Rinv[R]] == Jseite(regel[R]) for R in RS)}


REGELN = [dict(zip(RS, t)) for t in itertools.product(TEILM, repeat=4)]
EV_ = [(r, bewerte(r)) for r in REGELN]
VOLL = [(r, e) for r, e in EV_ if e["chiral"] and e["skalar"] and e["vektor"] and e["alle"]]
VOLL_J = [(r, e) for r, e in VOLL if e["J"]]
fmt = lambda r: {RNAME[R]: sorted(r[R]) for R in RS}
z114 = json.loads((HERE / "z114_ergebnis.json").read_text())["zahlen"]
z114_lsg = z114["Linienregeln"]["Restriktion"]["J-kovariante Lösungen"]
reprod = (len(VOLL) == z114["Linienregeln"]["Restriktion"]["chiral, Skalar, beide Vektoren, alle koppeln"]
          and len(VOLL_J) == z114["Linienregeln"]["Restriktion"]["davon J-kovariant"]
          and sorted((fmt(r) | {"N16": e["N16"]} for r, e in VOLL_J), key=str) == sorted(z114_lsg, key=str))
ZAHLEN["Linienregeln (Nachrechnung z114)"] = {"Sichtbarkeit": {s: {str(h): v for h, v in d.items()} for s, d in SICHT.items()},
                                              "voll": [fmt(r) | {"N16": e["N16"], "J": e["J"]} for r, e in VOLL]}
check("[E] z114 nachgerechnet: Linienrestriktion (x = 0 und ein zufälliges komplexes x) gibt Felddaten auf der Partnergeraden "
      f"P(pi) genau für h <= 0 ({[str(h) for h in HEL3 if SICHT['p'][h]]}), auf der Originalgeraden P(omega) für h >= 0 "
      f"({[str(h) for h in HEL3 if SICHT['o'][h]]}); von 256 Linienregeln sind {len(VOLL)} chiral mit Skalarleiter, beiden "
      f"Vektoren und koppeln überall, J-kovariant davon {len(VOLL_J)}: P_RJ (N16 = "
      f"{[e['N16'] for r, e in VOLL_J if r[1j] == frozenset({'o'})]}) und die Gegenorientierung (N16 = "
      f"{[e['N16'] for r, e in VOLL_J if r[1j] == frozenset({'p'})]}); identisch mit z114_ergebnis.json ({reprod}). Die sechs "
      "unterscheiden sich nur in der Fermionseite (R = i über o oder p, R = -i gegenüber) und in R = -1 (o, p oder beide)",
      reprod and len(VOLL) == 6 and len(VOLL_J) == 2 and all(SICHT["p"][h] == (h <= 0) for h in HEL3)
      and all(SICHT["o"][h] == (h >= 0) for h in HEL3) and all(r[1] == frozenset({"o", "p"}) for r, e in VOLL))

# --- 2.2 native Paarung als Bedingung: Geraden aus den Blättern, die die Erzeuger eines Sektors heben (Modell A)
MAPR = {"1": 1, "i": 1j, "-1": -1, "-i": -1j}
nat_ops = {k: O for k, O in TYPEN.items() if k[0] in ("Psi", "Psi-quer") and k[1] in ("psi*_a", "psi_a")}
spi_ops = {k: O for k, O in TYPEN.items() if k[0] in ("Phi", "Phi-quer") and k[1] in ("psi*_a", "psi_a")}
MVA, MBA = np.array([s[1][0] for s in ST]), np.array([s[1][1] for s in ST])


def geraden(ops, produkte=True):
    """Sektor R -> Menge der Seiten, deren Blatt die Erzeuger (einzeln und paarweise) heben: o bei Delta m > 0, p bei
    Delta mquer > 0 (ein Quant mit h > 0 hat Felddaten nur auf o, mit h < 0 nur auf p)."""
    lis = list(ops.values())
    kand = lis + ([X @ Y for X in lis for Y in lis] if produkte else [])
    out = {}
    for O in kand:
        O = sps.coo_matrix(O)
        if O.nnz == 0:
            continue
        for r, c in zip(O.row, O.col):
            R = MAPR[RSTR[int((T_[r] - T_[c]) % 4)]]
            seiten = ({"o"} if MVA[r] > MVA[c] else set()) | ({"p"} if MBA[r] > MBA[c] else set())
            out.setdefault(R, set()).update(seiten)
    return {R: frozenset(v) for R, v in out.items()}


G_nat, G_spi = geraden(nat_ops), geraden(spi_ops)
G_vek = {R: G_nat[R] | G_spi[R] for R in (1j, -1j)}
G_mod = {1j: frozenset({"o", "p"}), -1j: frozenset({"o", "p"})}  # Modeninhalt: R = i aus Original (o) und Partner (p) usw.
G_pos = {1j: frozenset({"o"}), -1j: frozenset({"o"})}            # Modulpositivität: nur h > 0


def mit_bedingung(lin):
    f = [(r, e) for r, e in EV_ if all(r[R] == lin[R] for R in (1j, -1j))]
    fv = [(r, e) for r, e in f if e["chiral"] and e["skalar"] and e["vektor"] and e["alle"]]
    fj = [(r, e) for r, e in fv if e["J"]]
    return {"Regeln": len(f), "voll": len(fv), "J-kovariant": len(fj), "N16": [e["N16"] for _, e in fj],
            "chiral überhaupt": sum(e["chiral"] for _, e in f), "Lösung": [fmt(r) for r, _ in fj]}


B_nat, B_spi, B_vek = mit_bedingung(G_nat), mit_bedingung(G_spi), mit_bedingung(G_vek)
B_mod, B_pos = mit_bedingung(G_mod), mit_bedingung(G_pos)
J_dann_nat = [e["N16"] for r, e in VOLL_J if r[1j] == G_nat[1j] and r[-1j] == G_nat[-1j]]
abgeleitet = {R: G_nat[R] for R in RS}
P_RJ = [r for r, e in VOLL_J if e["N16"] == 3][0]
gleich_PRJ = all(abgeleitet[R] == P_RJ[R] for R in RS)
ZAHLEN["native Paarung als Linienbedingung"] = {"Geraden nativ (abgeleitet)": fmt(abgeleitet),
                                                "Geraden Spiegelpaarung": {RNAME[R]: sorted(v) for R, v in G_spi.items()},
                                                "nativ": B_nat, "Spiegel": B_spi, "vektorartig (P_GSO allein)": B_vek,
                                                "Modeninhalt (seitenlokal)": B_mod, "Modulpositivität": B_pos}
check("[E, bedingt] Native Paarung als Bedingung. Aus Modell A abgeleitet (Seite = Blatt, das die nativen Erzeuger eines Sektors "
      f"heben, einzeln und paarweise): {fmt(abgeleitet)}. Mit der Fermionbedingung (R = i und R = -i) bleiben von 256 Regeln "
      f"{B_nat['Regeln']}, davon voll {B_nat['voll']}, J-kovariant {B_nat['J-kovariant']} mit N16 = {B_nat['N16']}; in "
      f"z114s Reihenfolge 6 -> 2 (J) -> {len(J_dann_nat)} (nativ), N16 = {J_dann_nat}. Die bosonischen Sektoren der "
      "nativen Algebra (Psi Psi-quer hebt beide Blätter, Psi Psi nur o, Psi-quer Psi-quer nur p) geben R = +-1 über beide "
      f"Geraden, die abgeleitete Regel ist P_RJ selbst ({gleich_PRJ}), ohne J. Kontrollen: Spiegelpaarung -> "
      f"{B_spi['J-kovariant']} Regel mit N16 = {B_spi['N16']}; P_GSO allein (beide Paarungen) -> chiral "
      f"{B_vek['chiral überhaupt']}; Modeninhalt (Original trägt beide R mit h > 0, Partner beide mit h < 0) -> chiral "
      f"{B_mod['chiral überhaupt']}; Modulpositivität (nur h > 0) -> chiral {B_pos['chiral überhaupt']} (universelle "
      "Gerade, Eichchiralität 0)",
      B_nat["Regeln"] == 16 and B_nat["voll"] == 3 and B_nat["J-kovariant"] == 1 and B_nat["N16"] == [3]
      and J_dann_nat == [3] and gleich_PRJ and B_spi["J-kovariant"] == 1 and B_spi["N16"] == [-3]
      and B_vek["chiral überhaupt"] == 0 and B_mod["chiral überhaupt"] == 0 and B_pos["chiral überhaupt"] == 0
      and G_nat[1j] == frozenset({"o"}) and G_nat[-1j] == frozenset({"p"}))

# --- 2.3 Vorzeichen und Zählung unter P_RJ aus der Paarung
o_line = {h: kuenneth(EM, n_pen(h)) for h in (Fr(1, 2), Fr(-1, 2))}
e_prj = bewerte(P_RJ)
kontr = {name: bewerte(P_RJ, m16=chiF[name])["N16"] for name in ("4 O", "Spiegelfamilie O(-4) + 3 O", "E_-^v (x) K")}
ZAHLEN["Vorzeichen"] = {"P_RJ": fmt(P_RJ), "Bewertung": e_prj, "16 auf P(omega), K2": {str(h): v for h, v in o_line.items()},
                        "Kontrollen N16": kontr}
check("[E unter P_IF] Vorzeichen: Die native Paarung legt R = i (die 16, z_X = i) auf die Originalgerade P(omega) "
      f"(Penrose-Konvention O(-2h - 2), sieht h >= 0) und R = -i (16q) auf die Partnergerade P(pi). Unter P_GSO und K2: "
      f"chiraler Slot h = +1/2 (H^0, H^1, H^2) = {o_line[Fr(1, 2)]}, Spiegelslot h = -1/2 {o_line[Fr(-1, 2)]}; N16 = "
      f"{e_prj['N16']}, N16q = {e_prj['N16q']} (CPT-gleich), Skalarleiter {e_prj['skalar']}, beide Vektoren "
      f"{e_prj['vektor']}. Die 16 ist damit Weyl mit h = +1/2 in der Konvention h = (m - mquer)/2 (Blatt des Originals). "
      f"Kontrollen: triviales 4 O {kontr['4 O']}, Spiegelfamilie {kontr['Spiegelfamilie O(-4) + 3 O']}, Serre-Dual "
      f"{kontr['E_-^v (x) K']}",
      o_line[Fr(1, 2)] == (0, 6, 0) and o_line[Fr(-1, 2)] == (0, 0, 0) and e_prj["N16"] == 3 and e_prj["N16q"] == 3
      and e_prj["chiral"] and e_prj["skalar"] and e_prj["vektor"]
      and kontr == {"4 O": 4, "Spiegelfamilie O(-4) + 3 O": 0, "E_-^v (x) K": -3})
befund("Teil 2: Unter der Prämisse P_IF (native Paarung) folgt P_RJ eindeutig, ohne J und ohne Auswahl unter den sechs: "
       "256 -> 16 -> 3 -> 1 bzw. 6 -> 2 -> 1, abgeleitet sogar die volle Regel. Das Vorzeichen ist fest: R = i auf die "
       "Originalgerade, N16 = +3. Die Spiegelpaarung gibt die Gegenorientierung; Modeninhalt, Modulpositivität und P_GSO "
       "allein geben keine Chiralität", True)

# =================================================================== Teil 3
say("\nTeil 3: Folgen für K16 unter T + P_GSO + K2 + P_IF: entfernt die Paarung den Spiegel auch aus dem Spektrum?")


def bfs(start, ops):
    adj = [sps.csc_matrix(O) for O in ops]
    seen, dq = {start}, deque([start])
    while dq:
        j = dq.popleft()
        for O in adj:
            col = O.getcol(j)
            for i in col.indices[np.abs(col.data) > 1e-14]:
                if i not in seen:
                    seen.add(i)
                    dq.append(i)
    return seen


nativ_alg = [O for k, O in TYPEN.items() if k[0] in ("Psi", "Psi-quer")]
nativ_alg = nativ_alg + [O.T.tocsr() for O in nativ_alg] + [kr(IL, X) for X in gen_A.values()]
reach0 = bfs(VAK, nativ_alg)
kl0 = Counter((RSTR[int(T_[i] % 4)], str(Fr(int(H2[i]), 2))) for i in reach0)
N0_alle = int(np.sum(NN == 0))
seed = SIX[((1, -1), (1, 1, 0, 0))]                               # tiefster R = -1-Sektor (z82, L_0 = 1/2), skalar
reach2 = bfs(seed, nativ_alg)
kl2 = Counter((RSTR[int(T_[i] % 4)], str(Fr(int(H2[i]), 2))) for i in reach2)
mit_W = bfs(VAK, nativ_alg + [W_rand, W_rand.T.tocsr()])
mit_D = bfs(VAK, nativ_alg + [D_nat, D_nat.T.tocsr()])
klW = Counter((RSTR[int(T_[i] % 4)], str(Fr(int(H2[i]), 2))) for i in mit_W)
klD = Counter((RSTR[int(T_[i] % 4)], str(Fr(int(H2[i]), 2))) for i in mit_D)
gso_alle = int(np.sum(NN % 2 == 0))
SPIEGEL, VEKT, SKAL2 = ("i", "-1/2"), (("1", "1"), ("1", "-1")), ("-1", "0")
ZAHLEN["Superauswahl"] = {"vom Vakuum erreicht (native Algebra)": len(reach0), "Zustände mit N = 0": N0_alle,
                          "N-Werte erreicht": sorted({int(NN[i]) for i in reach0}),
                          "Klassen vom Vakuum (R, h)": {f"({a}, {b})": v for (a, b), v in sorted(kl0.items())},
                          "vom R = -1-Skalar erreicht": {f"({a}, {b})": v for (a, b), v in sorted(kl2.items())},
                          "mit Randpaar W": len(mit_W), "mit (L^nat_-1)^2": len(mit_D), "P_GSO-Zustände (N gerade)": gso_alle}
check("[E] Superauswahl (Modell A, Breitensuche über die Basis, alle Erzeuger monomial). Die native Algebra (Psi, Psi-quer "
      "in beiden Realisierungen, ihre Adjungierten und die 17 su(2,2)-Erzeuger) erreicht vom Vakuum (Q = 0, (E, h) = (1, 0)) "
      f"genau {len(reach0)} Zustände, das sind alle {N0_alle} Zustände mit N = T - 2h = 0 im Fenster (N erreicht: "
      f"{sorted({int(NN[i]) for i in reach0})}): native Leitern (i, +1/2) {kl0[('i', '1/2')]}, (-i, -1/2) "
      f"{kl0[('-i', '-1/2')]}, neutrale Skalarleiter (1, 0) {kl0[('1', '0')]}, R = -1 mit h = +-1 "
      f"{kl0[('-1', '1')]} + {kl0[('-1', '-1')]}; Spiegel (i, -1/2) {kl0[SPIEGEL]}, R-neutrale Vektoren (1, +-1) "
      f"{kl0[VEKT[0]] + kl0[VEKT[1]]}, R = -1-Skalare {kl0[SKAL2]}. Ein einziger nativer X-ungerader Boson-Erzeuger hebt die "
      f"Superauswahl auf: mit dem Randpaar psi*_a psi_b erreicht man {len(mit_W)}, mit (L^nat_-1)^2 {len(mit_D)} Zustände, "
      f"also alle {gso_alle} P_GSO-Zustände (Spiegel {klW[SPIEGEL]} bzw. {klD[SPIEGEL]}, R-neutrale Vektoren "
      f"{klW[VEKT[0]] + klW[VEKT[1]]})",
      len(reach0) == N0_alle and {int(NN[i]) for i in reach0} == {0} and kl0[("i", "1/2")] > 0 and kl0[("-i", "-1/2")] > 0
      and kl0[("1", "0")] > 0 and kl0[SPIEGEL] == 0 and kl0[VEKT[0]] + kl0[VEKT[1]] == 0 and kl0[SKAL2] == 0
      and len(mit_W) == gso_alle and len(mit_D) == gso_alle and klW[SPIEGEL] > 0 and klD[SPIEGEL] > 0)

# --- 3.2 der Spiegel im Sektor X = -1: native Kleidung der R-neutralen Vektoren, Preis der Z4-Projektion
vek_tief = SIX[((0, 0), (1, 3, 0, 0))]                            # R = 1, h = -1, E = 2
psi0 = TYPEN[("Psi", "psi*_a", 0)]
e_vek = np.zeros(len(ST))
e_vek[vek_tief] = 1.0
bild_v = psi0 @ e_vek
ziel_kl = {(RSTR[int(T_[i] % 4)], str(Fr(int(H2[i]), 2))) for i in np.nonzero(np.abs(bild_v) > 1e-14)[0]}
ZIEL_TXT = sorted(kt((a, Fr(b))) for a, b in ziel_kl)
weyl = {}
for lab_, (Tq, h2) in (("nativ 16 (i, +1/2)", (1, 1)), ("Spiegel 16 (i, -1/2)", (1, -1)), ("nativ 16q (-i, -1/2)", (-1, -1)),
                       ("Spiegel 16q (-i, +1/2)", (-1, 1))):
    sel = [i for i, (q, s) in enumerate(ST) if q == ((1, 0) if Tq == 1 else (-1, 0)) and H2[i] == h2]
    weyl[lab_] = [sum(1 for i in sel if EVA[i] == n + 0.5) for n in range(1, 5)]
PX = sps.diags((1 + Xv + Xv**2 + Xv**3)/4)
kom_PX = max(nrm(PX @ O - O @ PX) for O in nativ_alg)
kom_PX_W, kom_PX_D = nrm(PX @ W_rand - W_rand @ PX), nrm(PX @ D_nat - D_nat @ PX)
pxd = np.real(PX.diagonal())
verlust = {name: (int(np.sum((NN % 2 == 0) & sel)), int(np.sum((pxd > 0.5) & sel)))
           for name, sel in (("Spiegel (i, -1/2)", (T_ % 4 == 1) & (H2 == -1)), ("R-neutrale Vektoren", (T_ % 4 == 0)
                                                                                 & (np.abs(H2) == 2)),
                             ("R = -1-Skalare", (T_ % 4 == 2) & (H2 == 0)), ("native 16", (T_ % 4 == 1) & (H2 == 1)),
                             ("neutrale Skalarleiter", (T_ % 4 == 0) & (H2 == 0)))}
ZAHLEN["Spiegel im Spektrum"] = {"Psi_0 auf (R = 1, h = -1, E = 2)": ZIEL_TXT,
                                 "Weyl-Zählung je Ladungssektor bei E = n + 1/2": weyl,
                                 "[P_X, native Algebra], [P_X, W], [P_X, (L^nat)^2]": [kom_PX, kom_PX_W, kom_PX_D],
                                 "P_GSO -> P_X (Zustände im Fenster)": verlust}
check("[E] Der Spiegel liegt im Sektor X = -1 zusammen mit den R-neutralen Vektoren und den R = -1-Skalaren, und die nativen "
      f"Vertizes verbinden sie: vom tiefsten R = -1-Skalar (Ladung (1, -1), z82) erreicht die native Algebra den Spiegel "
      f"(i, -1/2) {kl2[SPIEGEL]}-mal, die R-neutralen Vektoren {kl2[VEKT[0]] + kl2[VEKT[1]]}-mal, die native Leiter "
      f"(i, +1/2) {kl2[('i', '1/2')]}-mal; Psi_0 bildet den tiefsten Vektor (R = 1, h = -1, E = 2) auf {ZIEL_TXT} ab: "
      "die Spiegel-16 ist die native Kleidung der R-neutralen Vektorleiter. Beide Weyl-Leitern sind voll (je Ladungssektor "
      f"bei E = n + 1/2: {weyl}). Die Projektion P_X = (1 + X + X^2 + X^3)/4 vertauscht mit der nativen Algebra "
      f"({kom_PX:.0e}), nicht mit dem Randpaar ({kom_PX_W:.1f}) und (L^nat_-1)^2 ({kom_PX_D:.1f}); sie nimmt (P_GSO -> P_X, "
      f"Zustände im Fenster) {verlust}: Den Spiegel spektral zu entfernen kostet genau die R-neutralen Vektoren und die "
      "R = -1-Skalare (z102 in Vertexform)",
      kl2[SPIEGEL] > 0 and kl2[VEKT[0]] + kl2[VEKT[1]] > 0 and kl2[("i", "1/2")] == 0 and ziel_kl == {("i", "-1/2")}
      and weyl["nativ 16 (i, +1/2)"] == weyl["Spiegel 16 (i, -1/2)"] == [2, 6, 12, 20]
      and weyl["nativ 16q (-i, -1/2)"] == weyl["Spiegel 16q (-i, +1/2)"] == [2, 6, 12, 20]
      and kom_PX < 1e-12 and kom_PX_W > 0.5 and kom_PX_D > 0.5
      and all(v[0] > 0 and v[1] == 0 for k, v in verlust.items() if k in ("Spiegel (i, -1/2)", "R-neutrale Vektoren",
                                                                         "R = -1-Skalare"))
      and all(v[0] == v[1] > 0 for k, v in verlust.items() if k in ("native 16", "neutrale Skalarleiter")))

# --- 3.3 Linie gegen Spektrum unter T + P_GSO + K2 + P_IF
h0E, h1E = coh(EM, 0)                                            # Familiennullmoden der Naht (K2-Faktor)
artinhalt = {"native Leiter": [h0E*x for x in weyl["nativ 16 (i, +1/2)"]],
             "Spiegelleiter": [h0E*x for x in weyl["Spiegel 16 (i, -1/2)"]]}
linie = {"chiraler Slot": kuenneth(EM, n_pen(Fr(1, 2))), "Spiegelslot": kuenneth(EM, n_pen(Fr(-1, 2)))}
p_univ = bewerte({R: frozenset({"p"}) for R in RS})
ZAHLEN["Linie gegen Spektrum"] = {"Artinhalt K2 je 16 bei E = n + 1/2": artinhalt,
                                  "Linienkohomologie (Original, K2)": {k: v for k, v in linie.items()},
                                  "Netto-Index je 16": e_prj["N16"]}
check("[E/C] Linie gegen Spektrum unter T + P_GSO + K2 + P_IF (=> P_RJ): Auf der Linie hat jede 16 Netto-Index "
      f"{e_prj['N16']} (Familienindex chi(E_-) = {chiF['E_-']} mal Twistorindex; chiraler Slot {linie['chiraler Slot']}, "
      f"Spiegelslot {linie['Spiegelslot']}, nahtsteril). Im Spektrum ist der Artinhalt in K2 H^0(E_-) (x) Standardform "
      f"(h^0 = {h0E}, h^1 = {h1E}; z114), je 16 und E = n + 1/2 also {artinhalt}: drei native und drei Spiegelleitern (16 "
      "bei h = -1/2, 16q bei h = +1/2) mit gleicher Zählung. "
      "Die Paarung entfernt den Spiegel aus den Liniendaten und aus allen Einzelvertex-Anregungen des Vakuums (N = 0), nicht "
      "aus dem Spektrum: Dort ist er der Sektor X = -1, der mit den R-neutralen Vektoren steht und fällt. Gegenprobe: ohne "
      f"die Paarung (universelle Gerade) wäre der Linienindex N16 = {p_univ['N16']}, N16q = {p_univ['N16q']}, nicht chiral",
      e_prj["N16"] == 3 and linie["chiraler Slot"] == (0, 6, 0) and linie["Spiegelslot"] == (0, 0, 0)
      and (h0E, h1E) == (3, 0) and artinhalt["native Leiter"] == artinhalt["Spiegelleiter"] == [6, 18, 36, 60]
      and not p_univ["chiral"] and p_univ["N16"] == -3)
befund("Teil 3: Die native Paarung entfernt die Spiegelleiter nicht aus dem Spektrum. Sie macht N = T - 2h zur "
       "Superauswahlladung der Vertexalgebra: Der Vakuumsektor N = 0 enthält die nativen Leitern und die Skalarleiter, aber "
       "weder Spiegel noch R-neutrale Vektoren noch R = -1-Skalare. Der Spiegel lebt im Sektor X = -1, ist dort die native "
       "Kleidung der R-neutralen Vektoren und R = -1-Skalare und wird von nativen X-ungeraden Bosonen (Randpaar psi*_a "
       "psi_b, (L^nat_-1)^2) erreicht. Wer ihn spektral entfernt (P_X), entfernt die R-neutralen Vektoren mit", True)

# =================================================================== Bilanz
say("\nBilanz und Urteil")
STATUS = {
    "Klassen der gekleideten Vertizes": ("E", "P_GSO lässt nativ (i, +1/2), (-i, -1/2) und falsch (i, -1/2), (-i, +1/2) zu; "
                                              "Unterschied genau X = R i^(-2h)"),
    "CAR (gemeinsam, zwei Nähte)": ("E (negativ)", "64 Antikommutatoren 0 für beide Paarungen; Blatttausch Theta fixiert CAR, "
                                                   "G, R, E, su(2,2) und tauscht die Paarungen"),
    "Modeninhalt": ("E (negativ)", "Original trägt beide R mit h = +1/2: R = i bei L_0 = 5/4 (psi*_a, psi_b), R = -i bei "
                                   "7/4 (psi_a, psi*_b); V^- (x) S ist native Mode"),
    "Tomita J und S": ("E (negativ)", "beide Paarungen J-abgeschlossen, 4 J-invariante Klassenregeln, 2 chirale; S erhält "
                                      "nativ und falsch"),
    "Energiepositivität": ("E (negativ)", "E, P_0 helizitätsblind; Modulpositivität wählt h > 0 für beide R (universelle "
                                          "Gerade)"),
    "Lokalität": ("E (negativ)", "jede Klasse genau eine seitenlokale Realisierung"),
    "Energieordnung": ("E", "seitenlokal nativ 1/4 (L_0 = 5/4), falsch 3/4 (7/4), e^(-2 pi); gemischt umgekehrt"),
    "native Paarung": ("C (Prämisse P_IF)", "Ad Z_SU = Ad R (eine Verklebung I_T = +-I_F, z116) = X = 1 je Einzelquant = "
                                            "eta_+-Kleidung; ganzzahlig N = T - 2h; auf Fermionen = P_A3, auf Bosonen "
                                            "entgegengesetzt"),
    "P_IF => P_RJ": ("E unter P_IF", "256 -> 16 -> 3 -> 1, 6 -> 2 -> 1, abgeleitet = P_RJ; R = i -> Originalgerade, N16 = +3"),
    "Spiegel im Spektrum": ("E (negativ)", "Sektor X = -1 mit R-neutralen Vektoren und R = -1-Skalaren; native X-ungerade "
                                           "Bosonen verbinden; P_X kostet die Vektoren"),
    "Herkunft von P_IF": ("O", "kein geprüftes natives Datum wählt Ad Z_SU = +Ad R statt -Ad R (ein Bit); Kandidaten: "
                               "Nahtindex +1 von A_q (z116), Energieordnung alpha = 1/4"),
    "K16": ("L", "Linienindex 3 je 16 unter T, P_GSO, K2, P_IF; Spiegel spektral vorhanden, nahtsteril"),
    "gesamt": ("bedingt", "P_RJ auf P_IF zurückgeführt (ein Vorzeichenbit), Ausschluss der falschen Vertizes nicht nativ, "
                          "spektral gescoptes Negativ"),
}
for k_, (s2, t2) in STATUS.items():
    say(f"   [{s2}] {k_}: {t2}")
ZAHLEN["status"] = STATUS
ZAHLEN["z122_urteil_kurz"] = {
    "Urteil": "bedingt (Prämisse P_IF); Ausschluss der falschen Vertizes nicht nativ (gescoptes Negativ); spektral gescoptes "
              "Negativ",
    "K16": "bleibt L: Spiegel im Sektor X = -1 mit den R-neutralen Vektoren",
    "chiral_matter": "Twistorlinienindex mit P_RJ jetzt aus P_IF (Vorzeichen fest, N16 = +3); Modeninhalt und "
                     "Modulpositivität: nicht chiral; P_IF = P_A3 auf Fermionen"}
urteil = ("Die falschen gekleideten Vertizes V^+ (x) S-quer und V^- (x) S sind nicht nativ ausgeschlossen; die native Paarung "
          "ist eine benannte Prämisse, aus der P_RJ dann eindeutig folgt. (1) Alle vier Typen sind G-gerade, P_GSO lässt "
          "sie zu; sie unterscheiden sich genau im Z4-Charakter X = R i^(-2h). (2) CAR: In der gemeinsamen CAR zweier Nähte "
          "antikommutieren native und falsche Vertizes gleich (alle 64 Paare); der Blatttausch Theta (m <-> mquer) fixiert "
          "alle Moden, G, R, E und die su(2,2) und tauscht die Paarungen, also ist jede Bedingung aus diesen Daten blind. "
          "(3) Modeninhalt: Jede native Mode des Originals hebt nur dessen Blatt (Delta 2h >= 0). Teilchen liegen in der "
          "Schrittart Teilchen hoch, Löcher in Loch runter; das R = i-Dublett ist (psi*_a, psi_b) bei L_0 = 5/4, das "
          "R = -i-Dublett (psi*_b, psi_a) bei 7/4, beide mit h = +1/2. Der falsche Vertex V^- (x) S ist also eine native "
          "Mode, der native Psi-quer keine; das Original allein hat nur h >= 0 (universelle Gerade). (4) Tomita: J bildet "
          "jeden Elementarvertex auf sein Gegenstück der anderen Seite ab und erhält beide Paarungen (4 J-invariante "
          "Klassenregeln, 2 chirale); auf dem KMS-Vektor gilt a_c Omega = q b_c^+ Omega und V^-_o Omega = Delta-Gewicht mal "
          "V^+_p Omega, und S bildet nativ auf nativ, falsch auf falsch ab. (5) Energie und Lokalität: E und P_0 sind "
          "helizitätsblind; Modulpositivität wählt h > 0 für beide R und verwirft den nativen Psi-quer; jede Klasse hat "
          "genau eine seitenlokale Realisierung. Was trennt, ist allein die Ladungsenergie: seitenlokal liegt die native "
          "Klasse auf der tiefsten G-geraden Fermionstufe (1/4, L_0 = 5/4), die falsche 1/2 höher (e^(-2 pi)); gemischt "
          "realisiert ist es umgekehrt. Das ordnet, schließt aber nichts aus. (6) Die native Paarung ist damit die Prämisse "
          "P_IF: Verklebung der Twistorphase mit der Uhr, Ad Z_SU = Ad R mit dem Zentrum Z_SU = i^(2h - 2) von SU(2,2) (eine der "
          "beiden Verklebungen I_T = +-I_F aus z116), gleichwertig "
          "X = 1 auf jedem Einzelquant oder eta_+-Kleidung je Seite; ihre ganzzahlige Hebung N = T - 2h ist eine "
          "Erhaltungsgröße der Vertexalgebra, G = (-1)^N ist P_GSO. Als Spektralregel fällt P_IF auf den Fermionen mit P_A3 "
          "zusammen und ist auf den Bosonen entgegengesetzt. (7) Teil 2: Unter P_IF folgt P_RJ eindeutig (256 -> 16 -> 3 "
          "-> 1, in z114s Reihenfolge 6 -> 2 -> 1; die aus der Vertexalgebra abgeleitete volle Regel ist P_RJ selbst, ohne "
          "J), mit festem Vorzeichen: R = i (die 16) auf die Originalgerade P(omega), R = -i auf die Partnergerade, N16 = "
          "+3 unter P_GSO und K2. Die Spiegelpaarung gibt die Gegenorientierung (-3); Modeninhalt, Modulpositivität und "
          "P_GSO allein geben keine Chiralität. (8) Teil 3: Die Paarung entfernt die Spiegelleiter nicht aus dem Spektrum. "
          "Vom Vakuum erreicht die native Algebra genau den Sektor N = 0 (native Leitern, neutrale Skalarleiter, R = -1 mit "
          "h = +-1), ohne Spiegel, ohne R-neutrale Vektoren und ohne R = -1-Skalare. Der Spiegel lebt im Sektor X = -1; "
          "dort ist er die native Kleidung der R-neutralen Vektoren und R = -1-Skalare, und native X-ungerade Bosonen "
          "(Randpaar psi*_a psi_b mit L_0 = 1/2, (L^nat_-1)^2) verbinden die Sektoren. Unter T + P_GSO + K2 + P_IF bleiben "
          "drei Spiegelleitern (16 bei h = -1/2, 16q bei h = +1/2, je n(n+1)) im Spektrum, nahtsteril auf der Linie. Wer "
          "sie spektral entfernt, verliert die R-neutralen Vektoren (z102 in Vertexform). Folge für K16: bleibt L; der "
          "Spiegel ist jetzt genau verortet. Folge für chiral_matter: Der Wert 'Twistorlinienindex mit P_RJ' ruht nun auf "
          "P_IF (ein Vorzeichenbit wie P_vol in z103) statt auf einer freien Linienregel; Modeninhalt und Modulpositivität "
          "sind als Werte nicht chiral. Nichts ist befördert; 3+1-Aussagen bleiben Kinematik der Einpaarschicht.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z122_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
