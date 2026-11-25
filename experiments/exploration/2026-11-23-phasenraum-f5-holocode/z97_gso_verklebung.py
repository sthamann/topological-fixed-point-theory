"""Z97: Folgefrage G1 'Verklebung von Spin und Statistik: GSO-artige Regel (-1)^(2h) = R^2'.

Leitlinie 'Rekonstruktion statt Auswahl': eine Regel zählt nur als nativ, wenn sie aus vorhandener TFPT-Struktur folgt.
Firewalls: DIMENSION.UPLIFT.FIREWALL.01 (jede 3+1-Aussage ist Kinematik der radialen Einpaarschicht oder Lesart),
GRAV.NONCIRCULAR.01 (kein Schluss modular => Raumzeit). Gelöste Bausteine (16, E8, Flavour) werden benutzt, nicht neu
geöffnet. Bestand: z82 (R = i^(Q_a - Q_b), R^2 = (-1)^F, Klasse von R = z_X modulo G_SM, eta-Dublett, eta-SU(2) im
KMS-Zustand nicht implementierbar), z83 (radiale Standardform mit CCR-Schritten S_a, su(2,2), E = (m + mquer)/2,
h = (m - mquer)/2, Modulgenerator 2 beta h, Ad_R = 1 auf den Schritten), z91 (bosonische Komposition der Nähte), z92
(U_2pi = (-1)^(2h), Spinoren der Leiter sind R-Bosonen, E^nat = L_0^osc + Q^2/2, Familie U^osc c(Q) scheitert), z95
(J R J = R^-1), TFPT-STAND (E8 = D5 + A3 über Z4, A3 = Familiengitter mit q(A3) = 3/4).

Teil 1  Die Regel und ihr Abschluss. Kombinierter Raum Ladungsgitter (Q_a, Q_b) (x) radiale Standardform (z83, CCR, M = 7),
        G = U_2pi (x) R^2, P = (1 + G)/2 (gleichwertig: die 2 pi-Drehung mit der Uhr verdrehen und G = 1 verlangen).
        Abschluss unter su(2,2) (alle 15 Erzeuger), unter Vertexoperatoren (nackte Ladungsschritte, native Fermimoden
        psi*_x und der native Heber, gekleidete Vertizes V (x) S_a mit einseitigem Schritt), unter Fusion (Gradierung und
        Austauschvorzeichen). Was überlebt: (E, h, R), Weyl-Zählung, Index je R-Klasse.
Teil 2  Ist die Regel nativ? (a) Tomita J und Modulgruppe bei beta = 4 pi auf Original (x) Partner (Ladungen beider Seiten,
        KMS-Vektor), (b) Z4-Kleber von E8 = D5 + A3: Zentrum von SU(2,2) auf der radialen Leiter, Prämisse P_A3, Gehalt
        und Konflikte, (c) Spinstruktur der Doppelüberlagerung (eta-Dublett): zentrale involutive Hebungen, Ad_R-Charaktere.
Teil 3  Bilanz, Urteil für G1, Folgen für chirale Materie (T4).
"""
import itertools
import json
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
import scipy.sparse as sps
from scipy.linalg import expm

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
rng = np.random.default_rng(97)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


# ------------------------------------------------------------------ radiale Standardform (wie z83/z92, CCR-Schritte)
AMP_CCR = (lambda m, k: np.sqrt(k + 1), lambda m, k: np.sqrt(m - k))
SCHRITT = (lambda k: k + 1, lambda k: k)                          # S_1 Teilchen hoch, S_2 Loch runter (z83)


def standardform(M):
    """(+)_(m, mquer <= M) H_m (x) H-quer_mquer: P[a,b] = S_a (x) S-quer_b, A = Bilineare links, B rechts, E, und die
    einseitigen Schritte SL[a] = S_a (x) 1, SR[b] = 1 (x) S-quer_b."""
    lab = [(m, mb, k, kb) for m in range(1, M + 1) for mb in range(1, M + 1) for k in range(m) for kb in range(mb)]
    ix = {s: i for i, s in enumerate(lab)}
    n = len(lab)

    def dicht(r, c, w):
        return sps.csr_matrix((w, (r, c)), shape=(n, n)).toarray()

    P = {}
    for a in range(2):
        for b in range(2):
            r, c, w = [], [], []
            for (m, mb, k, kb), i in ix.items():
                if m < M and mb < M:
                    r.append(ix[(m + 1, mb + 1, SCHRITT[a](k), SCHRITT[b](kb))])
                    c.append(i)
                    w.append(AMP_CCR[a](m, k)*AMP_CCR[b](mb, kb))
            P[(a, b)] = dicht(r, c, w)

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
        return dicht(r, cc, w)

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
        SL[a], SR[a] = dicht(r, c, w), dicht(r2, c2, w2)
    return lab, ix, P, A, B, np.diag([0.5*(s[0] + s[1]) for s in lab]), SL, SR


def car(n):
    """Jordan-Wigner: Vernichter a_0..a_(n-1) auf (C^2)^n (wie z92)."""
    dim = 2**n
    ops = []
    for j in range(n):
        A_ = np.zeros((dim, dim))
        for s in range(dim):
            if (s >> j) & 1:
                A_[s ^ (1 << j), s] = (-1)**bin(s & ((1 << j) - 1)).count("1")
        ops.append(A_)
    return ops


E_lad = lambda al, Q: Fr(Q*Q, 2) + (al - Fr(1, 2))*Q              # Ladungsenergie eines Quarters (z82/z92)
E_nat = lambda al, Q: Fr(Q*Q, 2)                                   # native radiale Ladungsenergie (z92)
QUARTER = {"a": Fr(1, 4), "b": Fr(3, 4)}
R_von = lambda qa, qb: 1j**((qa - qb) % 4)                         # R = i^(Q_a - Q_b) (z81, z82)

# =================================================================== Teil 1
say("Teil 1: Die Regel (-1)^(2h) = R^2 auf Ladung (x) radialer Standardform und ihr Abschluss")

M = 7
LAB, IXS, PS, AS, BS, ES, SL, SR = standardform(M)
NR = len(LAB)
MV, MBV = np.array([s[0] for s in LAB]), np.array([s[1] for s in LAB])
ZH = MV - MBV                                                     # 2h
SG = np.where(ZH % 2 == 0, 1, -1)                                 # (-1)^(2h)
par1 = lambda x: 1 if x % 2 == 0 else -1
zstr = lambda x: next(s for s, v in (("1", 1), ("i", 1j), ("-1", -1), ("-i", -1j)) if np.isclose(x, v))
EV = 0.5*(MV + MBV)
HH = np.diag(ZH/2)

# --- 1.1 G, P und su(2,2)
DG = {(c_, a_): AS[(c_, a_)] - BS[(a_, c_)] for c_ in range(2) for a_ in range(2)}   # räumliche Drehungen (z83)
JX, JY, JZ = (DG[(0, 1)] + DG[(1, 0)])/2, (DG[(0, 1)] - DG[(1, 0)])/2j, (DG[(0, 0)] - DG[(1, 1)])/2
achse = rng.normal(size=3)
achse /= np.linalg.norm(achse)
JN = achse[0]*JX + achse[1]*JY + achse[2]*JZ
U2 = np.zeros((NR, NR), complex)
for m in range(1, M + 1):
    for mb in range(1, M + 1):
        sel = np.where((MV == m) & (MBV == mb))[0]
        U2[np.ix_(sel, sel)] = expm(2j*np.pi*JN[np.ix_(sel, sel)])
u2_fehler = float(np.abs(U2 - np.diag(SG)).max())
U2 = np.diag(SG.astype(float))
GEN = {f"P_{a}{b}": PS[(a, b)] for a in range(2) for b in range(2)}
GEN.update({f"K_{a}{b}": PS[(a, b)].T for a in range(2) for b in range(2)})
GEN.update({f"A_{c}{a}": AS[(c, a)] for c in range(2) for a in range(2)})
GEN.update({f"B_{c}{a}": BS[(c, a)] for c in range(2) for a in range(2)})
GEN["E"] = ES
komm_U = max(float(np.abs(U2 @ X - X @ U2).max()) for X in GEN.values())
komm_H = max(float(np.abs(HH @ X - X @ HH).max()) for X in GEN.values())
nichttrivial = sum(1 for X in GEN.values() if np.abs(X - np.diag(np.diag(X))).max() > 1e-9)
QWIN = range(-3, 4)
LADUNG = [(qa, qb) for qa in QWIN for qb in QWIN]
R2_lad = np.array([par1(qa + qb) for qa, qb in LADUNG])
G_komb = np.outer(R2_lad, SG)                               # G auf Ladung (x) Paar, diagonal
anteil_behalten = float((G_komb == 1).mean())
klassen = {}
for zh in range(-4, 5):
    for k in range(4):
        klassen[(zh, k)] = par1(zh) == par1(k)                    # R = i^k, R^2 = (-1)^k
erlaubt = sorted((str(Fr(zh, 2)), ["1", "i", "-1", "-i"][k]) for (zh, k), ok in klassen.items() if ok)
weyl = {}
for zh in (1, -1):
    weyl[str(Fr(zh, 2))] = [int(np.sum((ZH == zh) & np.isclose(EV, nn + 0.5))) for nn in range(1, 5)]
index_je_R = {}
for k in (1, 3):
    plus = sum(int(np.sum((ZH == 1) & (G_komb[j] == 1))) for j, (qa, qb) in enumerate(LADUNG) if (qa - qb) % 4 == k)
    minus = sum(int(np.sum((ZH == -1) & (G_komb[j] == 1))) for j, (qa, qb) in enumerate(LADUNG) if (qa - qb) % 4 == k)
    index_je_R[["1", "i", "-1", "-i"][k]] = {"h = +1/2": plus, "h = -1/2": minus, "Index": plus - minus}
ZAHLEN["regel"] = {"U_2pi - (-1)^(2h) (expm, zufällige Achse)": u2_fehler, "[U_2pi, su(2,2)] max": komm_U,
                   "[h, su(2,2)] max": komm_H, "nichtdiagonale Erzeuger": nichttrivial,
                   "Anteil behalten (|Q| <= 3, M = 7)": anteil_behalten, "erlaubte Klassen (h, R), |h| <= 2": erlaubt,
                   "Weyl-Zählung h = +-1/2 bei E = n + 1/2": weyl, "Index je R-Klasse (h = +1/2 minus h = -1/2)": index_je_R}
check("[E] Definition: auf Ladungsgitter (x) radialer Standardform (z83, CCR, M = 7, 784 Paarzustände, |Q_a|, |Q_b| <= 3) ist "
      f"U_2pi = exp(2 pi i n.J) blockweise (-1)^(2h) (Fehler {u2_fehler:.0e}), G = U_2pi (x) R^2 und P = (1 + G)/2. U_2pi und h "
      f"vertauschen mit allen Erzeugern P, K, A, B, E der radialen Leiter (17 Operatoren, Spann u(2,2) mit su(2,2), "
      f"{nichttrivial} davon nichtdiagonal; Residuen {komm_U:.0e}, "
      f"{komm_H:.0e}), R mit allem auf dem Paarfaktor: P schließt unter su(2,2) exakt (die Regel verdreht die 2 pi-Drehung mit "
      "der Uhr, G = U_2pi R^2 ist zentral)",
      u2_fehler < 1e-12 and komm_U == 0.0 and komm_H == 0.0 and nichttrivial == 12)
check("[E] Inhalt: es überleben genau die ganzzahligen Helizitäten mit R = +-1 und die halbzahligen mit R = +-i; der Sektor "
      "h = +-1/2 hat bei E = n + 1/2 die Weyl-Zählung n(n+1) = 2, 6, 12, 20 und bleibt genau mit ungerader Ladung. Aber für "
      f"jede Klasse R = +-i überleben h = +1/2 und h = -1/2 gleich oft (Index {[d['Index'] for d in index_je_R.values()]}): "
      "die Regel gibt Spin-Statistik, ist aber vektorartig (Dirac-Paare, keine Chiralität)",
      all(klassen[(zh, k)] == ((zh - k) % 2 == 0) for zh, k in klassen) and weyl["1/2"] == [2, 6, 12, 20]
      and weyl["-1/2"] == [2, 6, 12, 20] and all(d["Index"] == 0 and d["h = +1/2"] > 0 for d in index_je_R.values()))

# --- 1.2 Vertexoperatoren
# (a) nackte Ladungsschritte V^+ (x) 1
nackt = []
for j, (qa, qb) in enumerate(LADUNG):
    for dq in ((1, 0), (0, 1), (-1, 0), (0, -1)):
        ziel = (qa + dq[0], qb + dq[1])
        if ziel in LADUNG:
            jz = LADUNG.index(ziel)
            beh = G_komb[j] == 1
            nackt.append(float(np.mean(G_komb[jz][beh] == 1)))
nackt_max = max(nackt)

# (b) native Fermimoden psi*_x und der native Heber auf Maya-Zuständen (exakte Brüche, wie z92)
KMIN, KMAX = -9, 9


def bit(k):
    return 1 << (k - KMIN)


def nunter(s, k):
    return bin(s & ((1 << (k - KMIN)) - 1)).count("1")


def erzeuge(a, vec):
    out = {}
    for s, cc in vec.items():
        if (s >> (a - KMIN)) & 1:
            continue
        out[s | bit(a)] = out.get(s | bit(a), 0) + (-1 if nunter(s, a) % 2 else 1)*cc
    return out


def bil(a, b, vec, coef=1):
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
    return {k: x for k, x in out.items() if x != 0}


def partitionen(n, maxpart=None):
    maxpart = n if maxpart is None else maxpart
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitionen(n - k, k):
            yield (k,) + rest


def maya(Q, lam):
    s = 0
    for i in range(len(lam)):
        s |= bit(lam[i] - (i + 1) + Q)
    for y in range(KMIN, Q - len(lam)):
        s |= bit(y)
    return s


def L0(s, al):
    return sum((k + al) for k in range(0, KMAX + 1) if (s >> (k - KMIN)) & 1) - \
        sum((k + al) for k in range(KMIN, 0) if not (s >> (k - KMIN)) & 1)


def lnat_m1(vec):
    out = {}
    for k in range(KMIN, KMAX):
        for s2, x in bil(k + 1, k, vec, k + 1).items():
            out[s2] = out.get(s2, 0) + x
    return {k: x for k, x in out.items() if x != 0}


ZUST = [(Q, lam, maya(Q, lam)) for Q in range(-2, 3) for nn in range(4) for lam in partitionen(nn)]
ok_dosc, n_paare, erhalt = True, 0, {}
for al in QUARTER.values():
    for Q, lam, s in ZUST:
        for r in range(-3, 4):
            w = erzeuge(r, {s: 1})
            if not w:
                continue
            t = next(iter(w))
            d_osc = (L0(t, al) - E_lad(al, Q + 1)) - (L0(s, al) - E_lad(al, Q))
            ok_dosc &= d_osc == r - Q
            n_paare += 1
            erhalt.setdefault(r, {}).setdefault(Q, set()).add((d_osc + 1) % 2 == 0)   # Ad_G = (-1)^(d_osc) (-1)
eindeutig = all(len(v) == 1 for d in erhalt.values() for v in d.values())
sektoren = {r: sorted(Q for Q, v in d.items() if True in v) for r, d in erhalt.items()}
alle_sektoren = [r for r, qs in sektoren.items() if qs == sorted(erhalt[r])]
paritaet = all((True in erhalt[r][Q]) == ((r - Q) % 2 == 1) for r in erhalt for Q in erhalt[r])
kante = all(erhalt[Q][Q] == {False} for Q in range(-2, 3) if Q in erhalt and Q in erhalt[Q])
ok_heber = True
for Q, lam, s in ZUST:
    for al in QUARTER.values():
        for t in lnat_m1({s: 1}):
            ok_heber &= L0(t, al) - L0(s, al) == 1                # gleiches Q, L_0^osc + 1: G-ungerade
ZAHLEN["vertex"] = {"nackte Ladungsschritte: Anteil behalten -> behalten (max)": nackt_max,
                    "native Moden: (Zustand, Mode)-Paare": n_paare,
                    "Sektoren Q, in denen psi*_r die Regel erhält": {str(r): q for r, q in sorted(sektoren.items())}}
check(f"[negativ, E] Vertexoperatoren: nackte Ladungsschritte V^+- (x) 1 bilden behaltene Zustände vollständig auf "
      f"verworfene ab (Anteil behalten -> behalten {nackt_max}): R^2 kippt, h nicht. Mit 2h = L_0^osc - L-quer_0^osc auf allen "
      "Oszillatorzuständen (auf der Einpaarschicht m - mquer, wie U^osc in z92): Die nativen Fermimoden psi*_r verschieben "
      f"L_0^osc um r - Q ({n_paare} Paare, beide Quarter, exakt), also Ad_G psi*_r = -(-1)^(r - Q) psi*_r: die Regel bleibt "
      "genau in den Sektoren mit r - Q ungerade erhalten, keine Mode in allen Sektoren, das Feld psi* als Ganzes nicht. Die "
      "Kantenmode psi*_Q (|Q> -> |Q+1>, Oszillatorvakuum, das tiefste Fermion aus z92) und der native Heber L^nat_-1 "
      "(L_0^osc + 1 bei festem Q) sind G-ungerade: das tiefste geladene Fermion des Originals (h = 0, R^2 = -1) wird verworfen",
      nackt_max == 0.0 and ok_dosc and eindeutig and paritaet and not alle_sektoren and kante and ok_heber)

# (c) gekleidete Vertizes: Ladungsschritt (x) einseitiger Schritt
inner = MV <= M - 1
ok_kov = all(np.abs(U2 @ S @ U2 + S).max() == 0.0 for S in list(SL.values()) + list(SR.values()))
ok_spinor = True
for c_ in range(2):
    for d_ in range(2):
        for b_ in range(2):
            ok_spinor &= np.abs(AS[(c_, d_)] @ SL[b_] - SL[b_] @ AS[(c_, d_)] - (d_ == b_)*SL[c_]).max() < 1e-12
            ok_spinor &= np.abs(BS[(c_, d_)] @ SL[b_] - SL[b_] @ BS[(c_, d_)]).max() < 1e-12
            ok_spinor &= np.abs(BS[(c_, d_)] @ SR[b_] - SR[b_] @ BS[(c_, d_)] - (d_ == b_)*SR[c_]).max() < 1e-12
ok_E = all(np.abs(ES @ S - S @ ES - 0.5*S).max() < 1e-12 for S in list(SL.values()) + list(SR.values()))
ok_h = all(np.abs(HH @ SL[a] - SL[a] @ HH - 0.5*SL[a]).max() < 1e-12 and
           np.abs(HH @ SR[a] - SR[a] @ HH + 0.5*SR[a]).max() < 1e-12 for a in range(2))
ok_P = all(np.abs(SL[a] @ SR[b] - PS[(a, b)]).max() < 1e-12 and np.abs(SR[b] @ SL[a] - PS[(a, b)]).max() < 1e-12
           for a in range(2) for b in range(2))
grund = np.zeros(NR)
grund[IXS[(1, 1, 0, 0)]] = 1
bild = np.array([SL[a] @ grund for a in range(2)])
traeger = sorted({(int(MV[i]), int(MBV[i])) for v in bild for i in np.nonzero(np.abs(v) > 1e-12)[0]})
R_dressed = R_von(1, 0)
check("[E] gekleidete Vertizes Psi_a = V^+_a (x) S_a (Ladungsschritt mal einseitiger Schritt, S_a = native Schrittart aus "
      "z83) und Psi-quer_b = V^- (x) S-quer_b: U_2pi S U_2pi = -S exakt, also erhalten sie G. S_a ist ein "
      "(1/2, 0)-Spinor ([A_cd, S_b] = delta_db S_c, [B, S] = 0; S-quer ein (0, 1/2)-Spinor), hebt E um 1/2 und h um +1/2 "
      "(S-quer: -1/2) und liegt nicht in der Hülle von su(2,2) (h ist zentral). Aus dem skalaren Grundzustand (E, h) = (1, 0) "
      "macht Psi_a den Weyl-Tiefstgewichtsdoppelpunkt (3/2, 1/2) mit R = i. Ihre Produkte geben die radialen Heber: "
      "(V^+ (x) S_a)(V^- (x) S-quer_b) = 1 (x) P_ab. Die Schrittart ist der Spinorindex, die Ladung trägt die Statistik",
      ok_kov and ok_spinor and ok_E and ok_h and ok_P and traeger == [(2, 1)] and R_dressed == 1j)

# --- 1.3 Fusion und Austausch
dub = [IXS[(2, 1, 0, 0)], IXS[(2, 1, 1, 0)]]
Jd = [X[np.ix_(dub, dub)] for X in (JX, JY, JZ)]
cas_d = sum(X @ X for X in Jd)
I2 = np.eye(2)
Jt = [np.kron(X, I2) + np.kron(I2, X) for X in Jd]
Ut = expm(2j*np.pi*sum(achse[i]*Jt[i] for i in range(3)))
SW = np.eye(4)[[0, 2, 1, 3]]
cas_t = sum(X @ X for X in Jt)
sym, anti = (np.eye(4) + SW)/2, (np.eye(4) - SW)/2
vs, va = np.linalg.eigh(sym)[1][:, 1:], np.linalg.eigh(anti)[1][:, 3:]
cas_sym = np.round(np.linalg.eigvalsh(vs.conj().T @ cas_t @ vs), 9)
cas_anti = np.round(np.linalg.eigvalsh(va.conj().T @ cas_t @ va), 9)
c4 = car(4)
cd = [X.T for X in c4]
anti_ok = np.abs(cd[0] @ cd[2] + cd[2] @ cd[0]).max() < 1e-12 and np.abs(cd[0] @ cd[1] + cd[1] @ cd[0]).max() < 1e-12
R_paar = R_von(1, 0)*R_von(1, 0)
G_weyl_weyl = float(np.round(np.real(Ut[0, 0]), 12))*(R_paar**2).real
G_weyl_skalar = -1*(R_von(1, 0)**2).real
ZAHLEN["fusion"] = {"Casimir des Weyl-Dubletts": float(np.round(cas_d[0, 0].real, 12)),
                    "U_2pi auf zwei Weyl-Quanten": complex(np.round(Ut[0, 0], 12)),
                    "Casimir Sym^2 (bosonische Komposition, z91)": cas_sym.tolist(),
                    "Casimir Lambda^2 (gradierte Komposition)": cas_anti.tolist()}
check("[E] Fusion: G ist multiplikativ. Auf zwei Weyl-Quanten (h = 1/2, R = i) ist die 2 pi-Drehung aus der summierten "
      "Drehung J (x) 1 + 1 (x) J gleich +1 und R^2 = (+1): behalten; Weyl (x) Skalar hat U = -1, R^2 = -1: behalten; "
      "behalten (x) behalten liegt also in behalten. Austausch: im Original antikommutieren Fermimoden beider Quarter (CAR, "
      "Vorzeichen -1 = (-1)^(2h_1 2h_2) auf behaltenen Quanten). Die native Komposition der Nähte (z91, bosonisch, Sym^2) "
      f"gibt zwei Weyl-Quanten an der Schwelle E = 3 aber den Drehimpuls j = 1 (Casimir {cas_sym.tolist()}) statt j = 0 "
      f"(Lambda^2, Casimir {cas_anti.tolist()}): Spin-Statistik über Nähte hinweg braucht eine gradierte Komposition",
      np.allclose(cas_d, 0.75*np.eye(2)) and np.allclose(Ut, np.eye(4)) and G_weyl_weyl == 1.0 and G_weyl_skalar == 1.0
      and anti_ok and cas_sym.tolist() == [2.0, 2.0, 2.0] and cas_anti.tolist() == [0.0])
befund("Teil 1: Die Regel ist widerspruchsfrei und schließt unter su(2,2) und unter Fusion; behalten werden ganze Helizitäten "
       "mit R = +-1 und halbe mit R = +-i, mit Weyl-Zählung, aber vektorartig (beide Helizitäten je R). Unter den nativen "
       "Vertexoperatoren schließt sie nicht: nackte Ladungsschritte, die Kantenmode (tiefstes Fermion) und der einseitige "
       "Heber sind G-ungerade, psi*_r nur in den Sektoren r - Q ungerade. Geschlossen ist sie unter gekleideten Vertizes "
       "Ladung (x) einseitiger Schritt, deren Schrittart den Spinorindex liefert; die Verklebung ist genau diese Paarung", True)

# =================================================================== Teil 2
say("\nTeil 2: Ist die Regel nativ? Tomita J, E8-Kleber, Spinstruktur")

# --- 2.1 Tomita J und Modulgruppe bei beta = 4 pi auf Original (x) Partner
BETA = 4*np.pi
QW = range(-2, 3)
DOPP = [(qa, qb, pa, pb) for qa in QW for qb in QW for pa in QW for pb in QW]
dix = {d: i for i, d in enumerate(DOPP)}
perm_l = np.array([dix[(pa, pb, qa, qb)] for qa, qb, pa, pb in DOPP])
perm_r = np.array([IXS[(s[1], s[0], s[3], s[2])] for s in LAB])
par = lambda d, seite: par1(d[0] + d[1]) if seite == 0 else par1(d[2] + d[3])
G_orig = np.outer([par(d, 0) for d in DOPP], SG)
G_part = np.outer([par(d, 1) for d in DOPP], SG)
G_tot = G_orig*np.outer([par(d, 1) for d in DOPP], np.ones(NR))
J_von = lambda D: np.conj(D)[np.ix_(perm_l, perm_r)]               # J D J^-1 für diagonale D (antilinear, Tausch)
R_orig = np.outer([R_von(d[0], d[1]) for d in DOPP], np.ones(NR))
R_part = np.outer([R_von(d[2], d[3]) for d in DOPP], np.ones(NR))
j_orig = bool(np.array_equal(J_von(G_orig), G_part)) and not np.array_equal(G_orig, G_part)
j_tot = bool(np.array_equal(J_von(G_tot), G_tot))
j_R = bool(np.allclose(J_von(R_orig), 1/R_part))
anteil_jbruch = float(np.mean(G_orig != G_part))
beide = int(np.sum(G_tot == 1)) > 0 and int(np.sum(G_tot == -1)) > 0
QK = range(-8, 9)
gew = {(qa, qb): np.exp(-BETA*float(E_lad(QUARTER["a"], qa) + E_lad(QUARTER["b"], qb))) for qa in QK for qb in QK}
Zk = sum(gew.values())
p_ungerade = sum(w for (qa, qb), w in gew.items() if (qa + qb) % 2)/Zk
diag_l = [dix[(qa, qb, qa, qb)] for qa in QW for qb in QW]       # KMS-Vektor: Q = Q-quer, zonal |m,k> (x) |m,k>-quer
diag_r = [IXS[(m, m, k, k)] for m in range(1, M + 1) for k in range(m)]
ungerade_l = np.array([(qa + qb) % 2 == 1 for qa in QW for qb in QW])
kms_tot = bool(np.all(G_tot[np.ix_(diag_l, diag_r)] == 1))
kms_orig = bool(np.all(G_orig[np.ix_(diag_l, diag_r)][ungerade_l] == -1)
                and np.all(G_orig[np.ix_(diag_l, diag_r)][~ungerade_l] == 1))


# exakte Suche: Modulfluss e^{i theta (L_0 - L-quer_0)} mal U(1)^4-Ladungsdrehung gleich G? Phasen in Einheiten von 2 pi
def modulsuche(energie):
    treffer, zeiten = [], []
    for j8 in range(32):                                          # theta = j8 pi/8 in [0, 4 pi)
        t = Fr(j8, 16)                                            # theta/(2 pi)
        if any((t*d - Fr(d, 2)) % 1 for d in range(-(M - 1), M)):
            continue                                              # verfehlt U_2pi auf den neutralen Paaren
        zeiten.append(t)
        for al in QUARTER.values():
            for s_ in (1, -1):                                    # Original (+) oder Partner (-)
                phi = (Fr(1, 2) - t*s_*energie(al, 1)) % 1        # aus Q = 1 erzwungen
                rest = [(t*s_*energie(al, q) + phi*q - Fr(q % 2, 2)) % 1 for q in QW]
                if not any(rest):
                    treffer.append((t, al, s_))
    return zeiten, treffer


zeiten_lad, treffer_lad = modulsuche(E_lad)
zeiten_nat, treffer_nat = modulsuche(E_nat)
ordnung_lad = max((Fr(1, 2)*E_lad(al, q)).denominator for al in QUARTER.values() for q in QW)
ordnung_nat = max((Fr(1, 2)*E_nat(al, q)).denominator for al in QUARTER.values() for q in QW)
vierte = all((4*Fr(1, 2)*E_lad(al, q) - Fr(q % 2, 2)) % 1 == 0 for al in QUARTER.values() for q in QW)
zweite_nat = all((2*Fr(1, 2)*E_nat(al, q) - Fr(q % 2, 2)) % 1 == 0 for al in QUARTER.values() for q in QW)
zweite_diff = {str(al): str((Fr(1, 2)*(E_lad(al, 2) - 2*E_lad(al, 1) + E_lad(al, 0))) % 1) for al in QUARTER.values()}
ZAHLEN["tomita"] = {"Anteil G_orig != J G_orig J": anteil_jbruch, "KMS-Gewicht ungerader Ladung (beta = 4 pi)": p_ungerade,
                    "Modulwinkel theta/(2 pi), die U_2pi auf den Paaren geben": [str(t) for t in zeiten_lad],
                    "Treffer G (L_0)": len(treffer_lad), "Treffer G (E^nat)": len(treffer_nat),
                    "Ordnung von e^{i pi K} auf den Ladungen (L_0 / E^nat)": [ordnung_lad, ordnung_nat],
                    "zweite Ladungsdifferenz der Phase bei theta = pi (Einheiten 2 pi)": zweite_diff}
check("[E] Tomita J (Tausch Original <-> Partner mit Konjugation, Blätterwechsel m <-> mquer; J R J = R-quer^-1, vgl. z95): "
      f"J bildet G_orig = U_2pi R^2 auf G_Partner ab (verschieden auf {anteil_jbruch:.1%} der Zustände), J-kovariant ist nur "
      "G_tot = U_2pi (-1)^(F + F-quer). Diese Fassung hält den KMS-Vektor fest (Komponenten Q = Q-quer, zonal h = 0), die "
      f"einseitige Fassung schneidet sein Gewicht ungerader Ladung ab ({p_ungerade:.4f} bei beta = 4 pi). J vertauscht mit "
      "G_tot und lässt beide Eigenräume invariant: J ist verträglich, wählt aber nichts aus",
      j_orig and j_tot and j_R and beide and kms_tot and kms_orig and 0.07 < p_ungerade < 0.09)
check("[negativ, E] Modulgruppe bei beta = 4 pi (exakte Suche über theta = j pi/8 in [0, 4 pi) und beliebige U(1)^4-"
      "Ladungsdrehungen, Original und Partner, beide Quarter): U_2pi auf den neutralen Paaren erzwingt theta = pi mod 2 pi "
      f"(theta/2 pi = {[str(t) for t in zeiten_lad]}); dann bleibt in der Ladung die zweite Differenz e^(i theta) = -1 des "
      f"Terms Q^2/2 ({zweite_diff}), und kein Element ist gleich G (Treffer {len(treffer_lad)} mit L_0, {len(treffer_nat)} mit "
      f"E^nat). Das zu U_2pi passende Element hat auf den Ladungen die Ordnung {ordnung_lad} (vierte Potenz = (-1)^F: "
      f"Spin x_Z2 Z8, wie z92 an den Marken), mit E^nat die Ordnung {ordnung_nat} (Quadrat = (-1)^F, z92 kappa = i). Auf den "
      "Paaren ist U_2pi der Modulwinkel pi, auf den Ladungen ist R^2 der Modulwinkel 4 pi",
      zeiten_lad == [Fr(1, 4)*2, Fr(3, 4)*2] and zeiten_nat == zeiten_lad and not treffer_lad and not treffer_nat
      and ordnung_lad == 8 and ordnung_nat == 4 and vierte and zweite_nat and set(zweite_diff.values()) == {"1/2"})

# --- 2.2 E8 = D5 + A3 über Z4: Kleber, Zentrum von SU(2,2) auf der Leiter, Prämisse P_A3
wurzeln = []
for i, j in itertools.combinations(range(8), 2):
    for si, sj in itertools.product((1, -1), repeat=2):
        w = np.zeros(8)
        w[i], w[j] = si, sj
        wurzeln.append(w)
for s in itertools.product((0.5, -0.5), repeat=8):
    if sum(1 for x in s if x < 0) % 2 == 0:
        wurzeln.append(np.array(s))
WZ = np.array(wurzeln)
zX = np.exp(1j*np.pi*WZ[:, :5].sum(1))                            # Zentrum Spin(10): z_X = i auf der 16, -1 auf der 10
zA = np.exp(1j*np.pi*WZ[:, 5:].sum(1))                            # Zentrum Spin(6) = SU(4)
n5 = np.round((WZ[:, :5]**2).sum(1), 9)
n3 = np.round((WZ[:, 5:]**2).sum(1), 9)
kleber = float(np.abs(zX*zA - 1).max())
kl = {}
for a5, a3 in sorted(set(zip(n5, n3))):
    sel = (n5 == a5) & (n3 == a3)
    kl[f"D5-Norm {a5}, A3-Norm {a3}"] = {"Anzahl": int(sel.sum()),
                                         "(z_X, z_A3)": sorted({(complex(np.round(x, 9)), complex(np.round(y, 9)))
                                                                for x, y in zip(zX[sel], zA[sel])}, key=str)}
q_A3 = float(n3[(n5 == 1.25)][0])
# radiales Zentrum: Bild von i 1_4 aus nativen Erzeugern, zwei Darstellanten
J3L, J3R = (AS[(0, 0)] - AS[(1, 1)])/2, (BS[(0, 0)] - BS[(1, 1)])/2
Z_SU = expm(1j*np.pi*(ES + 2*J3R))
Z_SU2 = expm(-1j*np.pi*(ES + 2*J3L))
E4 = np.diag([0.5, 0.5, -0.5, -0.5])
tw = [expm(1j*np.pi*(E4 + 2*sg*np.diag([0, 0, 0.5, -0.5]))) for sg in (1, -1)] + \
     [expm(-1j*np.pi*(E4 + 2*sg*np.diag([0.5, -0.5, 0, 0]))) for sg in (1, -1)]
twistor_ok = all(np.allclose(T, 1j*np.eye(4)) for T in tw)
zentral = max(float(np.abs(Z_SU @ X - X @ Z_SU).max()) for X in GEN.values())
z_formel = float(np.abs(np.diag(Z_SU) - (-(1j**(ZH % 4)))).max())
Z_h = np.diag(1j**(ZH % 4))                                       # e^{i pi h}: Modulfluss der Standardform beim Viertelwinkel
regel_su = {}
for zh in range(-4, 5):
    zsu = -(1j**(zh % 4))
    for k in range(4):
        regel_su[(zh, k)] = np.isclose(1j**k*zsu, 1)              # Kleber z_X z_A3 = 1 mit R = z_X, z_A3 -> Z_SU
teil_gso = all(klassen[key] for key, ok in regel_su.items() if ok)
halb = sum(regel_su.values())*2 == sum(klassen.values())
chiral = {["1", "i", "-1", "-i"][k]: sorted(str(Fr(zh, 2)) for zh in range(-3, 4) if regel_su[(zh, k)]) for k in range(4)}
skalar_verboten = not regel_su[(0, 0)]
inhalt = {}
for name, wert in (("z_X = 1 (45 + 15)", 1), ("z_X = -1 (10, 6)", -1), ("z_X = i (16, 4)", 1j), ("z_X = -i (16q, 4q)", -1j)):
    k = [1, 1j, -1, -1j].index(wert)
    inhalt[name] = {"Wurzeln": int(np.sum(np.isclose(zX, wert))),
                    "Helizitäten (|h| <= 3/2)": sorted(str(Fr(zh, 2)) for zh in range(-3, 4) if regel_su[(zh, k)])}
ZAHLEN["E8_kleber"] = {"max |z_X z_A3 - 1| auf 240 Wurzeln": kleber, "Klassen": kl, "A3-Norm des Kleberanteils": q_A3,
                       "Z_SU zentral (max Kommutator)": zentral, "Z_SU - (-i^(2h))": z_formel,
                       "Regel aus Kleber + P_A3: erlaubte h je R": chiral, "4D-Lesart der E8-Wurzeln": inhalt}
check("[E] E8-Kleber (240 Wurzeln, D5 auf den ersten 5, A3 = D3 auf den letzten 3 Koordinaten): z_X = exp(i pi Summe w_D5) "
      "und z_A3 = exp(i pi Summe w_A3) erfüllen z_X z_A3 = 1 auf allen Wurzeln (Abweichung "
      f"{kleber:.0e}): (45,1), (1,15) haben (1, 1), (10, 6) hat (-1, -1), (16, 4) hat (i, -i); der Kleberanteil von A3 hat die "
      f"Norm {q_A3} = q(A3), die Gluenorm des gelösten Flavourblocks (A3 = Familiengitter, TFPT-STAND)",
      kleber < 1e-12 and q_A3 == 0.75 and inhalt["z_X = i (16, 4)"]["Wurzeln"] == 64
      and inhalt["z_X = -1 (10, 6)"]["Wurzeln"] == 60 and inhalt["z_X = 1 (45 + 15)"]["Wurzeln"] == 52)
check("[E] Zentrum von SU(2,2) auf der radialen Leiter: aus den nativen Erzeugern ist Z_SU = exp(i pi (E + 2 J3_R)) = "
      "exp(-i pi (E + 2 J3_L)) (zwei Darstellanten, im Twistor 4 beide gleich i 1_4, für beide Orientierungen), zentral "
      f"(Kommutator {zentral:.0e} mit allen Erzeugern) und gleich i^(2h - 2) = -i^(2h) (Abweichung {z_formel:.0e}), der "
      "Homogenitätscharakter des Twistortransformats (Helizität h <-> O(2h - 2)). Z_SU^2 = U_2pi, Z_SU^4 = 1. Der Viertelwinkel "
      "des Modulflusses der Standardform, e^(i pi h), ist -Z_SU und kein Bild des Zentrums",
      twistor_ok and np.allclose(Z_SU, Z_SU2) and zentral < 1e-12 and z_formel < 1e-12
      and np.allclose(Z_SU @ Z_SU, U2) and np.allclose(np.linalg.matrix_power(Z_SU, 4), np.eye(NR))
      and np.allclose(Z_h, -Z_SU))
check("[C, mit Konflikt] Liest man A3 als komplexifizierte radiale su(2,2) (Prämisse P_A3: z_A3 -> Z_SU), so gibt der Kleber "
      "R = Z_SU^-1, also quadriert genau die Regel (jede erlaubte Klasse ist GSO-behalten, die Hälfte davon) und dazu eine "
      f"chirale Z4-Regel: erlaubte h je R {chiral}. Die E8-Wurzeln lesen sich dann als SO(10)-Inhalt: 45 Vektoren (h = +-1), "
      "10 Skalare (h = 0), 16 eine Helizität (h = +1/2), 16-quer die andere. Aber: die neutrale radiale Skalarleiter aus "
      "z83 (h = 0, R = 1) ist kleberverboten (Z_SU = -1), ebenso neutrales h = +-2, und P_A3 deutet das Familiengitter A3 "
      "(q(A3) = 3/4 im gelösten Flavourblock) in die Raumzeit um",
      teil_gso and halb and skalar_verboten and chiral["i"] == ["-3/2", "1/2"] and chiral["1"] == ["-1", "1"]
      and chiral["-1"] == ["0"] and inhalt["z_X = 1 (45 + 15)"]["Helizitäten (|h| <= 3/2)"] == ["-1", "1"]
      and inhalt["z_X = -1 (10, 6)"]["Helizitäten (|h| <= 3/2)"] == ["0"])

# --- 2.3 Spinstruktur der Doppelüberlagerung (eta-Dublett)
S_eta = np.array([[0, -1], [1, 0]])                               # Spinorhebung (u, w) -> (-iu, -w) (z82)
S_triv = np.array([[0, 1], [1, 0]])
ew_eta = sorted(np.round(np.linalg.eigvals(S_eta), 12), key=lambda x: x.imag)
lam = {"a": 1j, "b": -1j}                                         # eta-Eigenwert eines Teilchens je Quarter (z82)
eta_paar = {q: lam[q]*np.conj(lam[q]) for q in lam}               # Teilchen und Loch desselben Quarters: Paar neutral
eta_lift = np.outer(R2_lad, np.ones(NR))                          # 2 pi-Hebung der eta-Struktur: 1 auf Paaren, R^2 auf Ladung
weg = float(np.mean(eta_lift*np.outer(np.ones(len(LADUNG)), SG) != G_komb))
G_unterschied = float(np.mean(eta_lift != np.outer(np.ones(len(LADUNG)), SG)))
hebungen = [k for k in range(4) if all(np.isclose((1j**k)**(2*q), 1) for q in range(4))]   # (U (x) R^k)^2 = 1 für alle R
Rd = np.diag([np.prod([[1j, 1j, -1j, -1j][j] if (s >> j) & 1 else 1 for j in range(4)]) for s in range(16)])
adR = lambda X: next((ch for ch in (1, 1j, -1, -1j) if np.allclose(Rd @ X @ np.linalg.inv(Rd), ch*X)), None)
eta_ch = [adR(cd[0]), adR(cd[2])]
schritt_ch = [1, 1]                                               # R = 1 auf der ganzen Leiter (z83): Ad_R S_a = S_a
hom_R = sum(1 for x in eta_ch for y in schritt_ch if np.isclose(x, y))
hom_gekleidet = sum(1 for x in eta_ch for y in schritt_ch if np.isclose(x, y*1j))
ZAHLEN["eta_spinstruktur"] = {"Eigenwerte der Spinorhebung": [str(x) for x in ew_eta],
                              "eta-Produkt auf neutralen Paaren": {q: complex(v) for q, v in eta_paar.items()},
                              "Anteil eta-Hebung != U_2pi (= verworfener Anteil)": G_unterschied,
                              "involutive zentrale Hebungen U (x) R^k": hebungen,
                              "Ad_R eta-Dublett": [zstr(x) for x in eta_ch], "dim Hom_R(eta, Schritte)": hom_R,
                              "dim Hom_R nach Ladungsschritt (Charakter i)": hom_gekleidet}
check("[negativ, E] Spinstruktur der Doppelüberlagerung: die eta-Spinorhebung hat die Eigenwerte {-i, i} und Quadrat -1; "
      "ihre 2 pi-Hebung ist auf den Ladungen R^2 = (-1)^F (z82), auf jedem neutralen Paar aber 1 (Teilchen und Loch "
      "desselben Quarters, eta-Produkt lambda lambda-quer = 1), für jede Helizität. Sie weicht von U_2pi genau auf den "
      f"ungeraden 2h ab ({G_unterschied:.1%} des kombinierten Raums), also genau dort, wo die Regel etwas fordert. Zentrale "
      f"involutive Hebungen der 2 pi-Drehung auf Ladung (x) Paar sind genau U (x) R^k mit k in {hebungen}: zwei Spin"
      "strukturen (unverdreht, verdreht), die die eta-Daten (R, (-1)^F, Ad_R auf den Schritten) nicht unterscheiden. "
      f"Ad_R ist {[zstr(x) for x in eta_ch]} auf dem eta-Dublett und 1 auf dem Schrittdublett: Hom_R = {hom_R}; erst nach "
      f"einem Ladungsschritt (Charakter i) gibt es {hom_gekleidet} äquivariante Zuordnungen, und das ist die Kleidung aus 1.2",
      np.allclose(sorted(ew_eta, key=lambda x: x.imag), [-1j, 1j]) and np.allclose(S_eta @ S_eta, -np.eye(2))
      and np.allclose(S_triv @ S_triv, np.eye(2)) and all(np.isclose(v, 1) for v in eta_paar.values())
      and hebungen == [0, 2] and hom_R == 0 and hom_gekleidet == 2 and abs(G_unterschied - (1 - anteil_behalten)) < 1e-12
      and eta_ch == [1j, -1j] and weg == 0.0)
befund("Teil 2: Keine native Struktur erzwingt die Regel. Tomita J ist verträglich (nur in der Fassung mit der Gesamtparität "
       "von Original und Partner) und wählt nichts aus; die Modulgruppe gibt wegen Q^2/2 nie G, ihr zu U_2pi passendes "
       "Element hat auf den Ladungen Ordnung 8. Die eta-Spinstruktur ist rein innen und lässt zwei Hebungen offen. Der "
       "E8-Kleber gibt die Regel (sogar chiral) nur unter P_A3, und P_A3 verbietet die neutrale radiale Skalarleiter und "
       "deutet das Familiengitter um", True)

# =================================================================== Teil 3
say("\nTeil 3: Bilanz und Urteil für G1")
STATUS = {
    "Regel G = U_2pi R^2 auf Ladung (x) Paar": ("E", "zentral, schließt unter su(2,2) und Fusion; behalten: ganze h mit R = +-1, "
                                                     "halbe h mit R = +-i; Weyl-Zählung n(n+1)"),
    "Chiralität unter der Regel": ("E (negativ)", "Index 0 je R-Klasse: beide Helizitäten je R, vektorartig"),
    "native Vertexoperatoren": ("E (negativ)", "nackte Ladungsschritte, Kantenmode psi*_Q (tiefstes Fermion) und L^nat_-1 "
                                               "G-ungerade; psi*_r nur in Sektoren r - Q ungerade"),
    "gekleidete Vertizes V (x) S_a": ("E", "erhalten G; S_a ist (1/2, 0)-Spinor, E + 1/2, h + 1/2; Produkte = P_ab"),
    "Austausch": ("E / O", "CAR im Original richtig; native Nahtkomposition (bosonisch) gibt j = 1 statt j = 0: gradierte "
                           "Komposition als Prämisse"),
    "Tomita J bei beta = 4 pi": ("E (negativ als Herleitung)", "verträglich nur als G_tot mit (-1)^(F + F-quer), hält den "
                                                              "KMS-Vektor; wählt nicht aus"),
    "Modulgruppe": ("E (negativ)", "kein Element gleich G; passendes Element Ordnung 8 (Z8) bzw. 4 mit E^nat"),
    "E8-Z4-Kleber": ("C mit P_A3, im Konflikt", "R = Z_SU^-1 gibt die Regel und chirale Z4; verbietet die neutrale "
                                                 "Skalarleiter; A3 ist das Familiengitter"),
    "eta-Spinstruktur": ("E (negativ)", "innen; zwei zentrale Hebungen; Hom_R(eta, Schritte) = 0"),
    "G1 gesamt": ("bedingt", "Prämisse P_GSO: die physikalische 2 pi-Drehung ist die mit der Uhr verdrehte Hebung "
                             "(Spin x_Z2 Z4 als Symmetriegruppe), nicht nativ"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
ZAHLEN["G1_urteil_kurz"] = {"Regel": "widerspruchsfrei, nicht nativ", "Prämisse": "P_GSO (verdrehte Spinstruktur)",
                            "Chiralität": "nein unter Z2; Z4-Verfeinerung nur mit P_A3 und Verlust der Skalarleiter"}
urteil = ("G1 ist bedingt, nicht nativ. Die Regel (-1)^(2h) = R^2 auf Ladungsgitter (x) radialer Standardform ist "
          "widerspruchsfrei. G = U_2pi R^2 ist zentral; die Projektion schließt exakt unter allen 15 Erzeugern der radialen "
          "su(2,2) und unter Fusion, denn G ist multiplikativ. Es überleben die ganzzahligen Helizitäten mit R = +-1 und die "
          "halbzahligen mit R = +-i. Der Sektor h = 1/2 trägt die Weyl-Zählung n(n+1) und bleibt genau mit ungerader Ladung. "
          "Unter den nativen Vertexoperatoren schließt die Regel nicht: Nackte Ladungsschritte sind G-ungerade, die "
          "Fermimoden psi*_r erhalten sie nur in den Sektoren mit r - Q ungerade, und das tiefste Fermion des Originals "
          "(Kantenmode, h = 0, R^2 = -1) wird verworfen, ebenso der einseitige native Heber. Geschlossen ist sie unter "
          "gekleideten Vertizes, Ladungsschritt mal einseitigem Schritt S_a. Die Schrittart aus z83 ist dabei der Spinorindex, "
          "(1/2, 0) mit E + 1/2 und h + 1/2, und die Produkte zweier gekleideter Vertizes sind die radialen Heber P_ab. Die "
          "Verklebung ist genau diese Paarung. Austausch: Im Original stimmt das Vorzeichen (CAR). Die native bosonische "
          "Komposition der Nähte gäbe zwei Weyl-Quanten an der Schwelle j = 1 statt j = 0; über Nähte hinweg braucht es eine "
          "gradierte Komposition. Herkunft: Tomita J ist verträglich, aber nur in der Fassung mit der Gesamtparität von "
          f"Original und Partner, die den KMS-Vektor festhält; die einseitige Fassung schnitte {p_ungerade:.1%} seines "
          "Gewichts ab. J wählt nichts aus. Kein Element der Modulgruppe mal Ladungsdrehungen ist G: U_2pi auf den Paaren "
          "erzwingt den Modulwinkel pi, und dann verdirbt der Term Q^2/2 die Ladungen; das passende Element hat dort "
          "Ordnung 8. Die eta-Spinstruktur ist rein innen und lässt zwei zentrale Hebungen U (x) R^k (k = 0, 2) offen; "
          "Hom_R(eta, Schritte) = 0. Der E8-Kleber z_X z_A3 = 1 gibt die Regel nur unter der Prämisse P_A3, dass A3 die "
          "komplexifizierte radiale su(2,2) ist. Das Zentrum von SU(2,2) wirkt auf der Leiter als i^(2h - 2), und dann folgt "
          "sogar eine chirale Z4-Regel mit dem SO(10)-Inhalt (45 Vektoren, 10 Skalare, 16 eine Helizität). P_A3 verbietet "
          "aber die neutrale radiale Skalarleiter aus z83 und deutet das Familiengitter A3 (q(A3) = 3/4) um. Urteil: bedingt "
          "mit der benannten Prämisse P_GSO, die physikalische 2 pi-Drehung sei die mit der Uhr verdrehte Hebung (Spin x_Z2 Z4 "
          "als Symmetriegruppe). Folge für T4: Die Z2-Regel liefert Spin-Statistik, aber keine Chiralität (Index 0 je "
          "R-Klasse). Chiralität bräuchte die Z4-Verfeinerung R = Z_SU^-1, und die ist nur mit P_A3 zu haben, auf Kosten der "
          "neutralen Skalarleiter. Nichts ist befördert; 3+1-Aussagen bleiben Kinematik der Einpaarschicht.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z97_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
