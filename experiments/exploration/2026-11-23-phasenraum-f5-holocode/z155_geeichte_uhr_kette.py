"""Z155: Die geeichte Uhr auf der nativen Kette (Gitter-Schwinger-Modell). Folgt P_uhr aus der Eichung?

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart, N numerisch, O offen.

Frage (Welle B nach H6h50): z152 liest die Naht als Hälfte der Kette H = -t sum (c+_n c_(n+1) + h.c.) + m sum (-1)^n n_n
mit m = -+2t bei halber Füllung (Kogut-Susskind-Fermion, frei). Das wechselwirkende Gatter des Kandidaten ist
U(phi) = B^+ e^(-i phi (Q_A^2 + Q_B^2)/2) B bei phi = pi/2 (z120, z145; Prämisse P_uhr = P_null + P_viertel). Eicht man
die Ladung auf der Kette, entsteht das Gitter-Schwinger-Modell; in 1+1 Dimensionen gibt das Gauß-Gesetz
H_el = (g^2/2) sum_Links E_l^2 mit E_l = Hintergrundfeld + Ladung links davon. Folgt P_uhr (Form und Viertelwinkel) aus
'die Uhr ist auf der nativen Kette geeicht' (P_eich)?

Bestand: z113 (L_0 = (Q -+ 1/4)^2/2 - 1/32 + Oszillatoren, e^(2 pi i L_0) = R), z120/z121 (Gatter, nativ trivial), z133
(P_null <=> P_sl2; inneres Gatter diagonal in (Q_A, Q_B) ist einkörperig), z135 (für Yang-Mills fehlen P_link, P_E), z144
(P_uhr), z145 (drei Formen von P_uhr, Z_16-Sektoren n = 4Q - 1, acht Uhrgatter), z152 (Kette, K = c mal Gitter-Boost).
Literatur (aus dem Gedächtnis zitiert, in dieser Rechnung nicht nachgeschlagen; alle benutzten Identitäten sind unten
selbst geprüft): Kogut/Susskind 1975; Banks/Susskind/Kogut 1976; Hamer/Weihong/Oitmaa 1997 (Gauß-Gesetz eliminiert);
Dempsey/Klebanov/Pufu/Zan, arXiv:2206.05308 (Massenverschiebung g^2 a/8 und Verschiebung um einen Platz als diskrete
chirale Drehung); Villain-Form und Selbstdualität der Z_N-Modelle (Elitzur/Pearson/Shigemitsu 1979,
Ukawa/Windey/Guth 1980).

Konvention: Plätze n = 0..N-1, Link l = 1..N-1 liegt zwischen den Plätzen l-1 und l, also ist jeder Link ein Schnitt
und Q_l die Ladung der Halbkette links davon. t = 1, m = -2 (gerade Plätze bevorzugt; z152: Leiter pi (4Z + 1)).
Zwei Hintergründe: 'gleichmäßig' q_n = n_n - 1/2 (Feld F0 am linken Rand) und 'gestaffelt' q_n = n_n - (1 auf geraden
Plätzen), die übliche Kogut-Susskind-Form mit theta = 2 pi F.

Teil A  Form, exakt: H_el als quadratische Form in den Ladungen; Kreuzterm; Modenbasis; Rest.
Teil B  Hintergrundfeld und das Viertel: z113-Formel, Konstante, Uhr; ist F = -+1/4 nativ; reine Z_4-Eichtheorie.
Teil C  Der Winkel: Clifford-Bedingung, Modulperiode und Boost-Gewicht, Z_8 der Uhrgatter, Selbstdualität.
Teil D  Exakte Diagonalisierung (N <= 12): Bleibt die Naht die Naht?
Teil E  Bilanz und Urteil.
"""
import json
import os
import time
from fractions import Fraction as Fr
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # deterministische Rundung (einfädig)

import mpmath as mp
import numpy as np
import scipy.sparse as sps
import sympy as sp
from scipy.linalg import expm
from scipy.optimize import brentq, least_squares
from scipy.sparse.linalg import eigsh

mp.mp.dps = 60
HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
TOL = 1e-10
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


# ---------------------------------------------------------------------------------------------------------------------
# Werkzeuge: Kette mit eliminiertem Gauß-Gesetz (halbe Füllung), Verschränkungsspektrum der Halbkette
class Kette:
    def __init__(self, N):
        self.N = N
        self.B = [b for b in range(1 << N) if bin(b).count("1") == N // 2]
        self.idx = {b: i for i, b in enumerate(self.B)}
        self.occ = np.array([[(b >> n) & 1 for n in range(N)] for b in self.B], float)
        rr, cc, r2, c2 = [], [], [], []
        for i, b in enumerate(self.B):
            for n in range(N - 1):
                if ((b >> n) & 1) != ((b >> (n + 1)) & 1):
                    rr.append(self.idx[b ^ (3 << n)])
                    cc.append(i)
                    if n % 2 == 1:
                        r2.append(rr[-1])
                        c2.append(i)
        self.hop = sps.csr_matrix((-np.ones(len(rr)), (rr, cc)), shape=(len(self.B),) * 2)
        self.hop2 = sps.csr_matrix((-np.ones(len(r2)), (r2, c2)), shape=(len(self.B),) * 2)   # Bindungen (n, n+1), n ungerade
        self.stag = (self.occ * np.array([(-1.0) ** n for n in range(N)])).sum(1)
        self.Qg = np.cumsum(self.occ - 0.5, axis=1)[:, :-1]                     # gleichmäßiger Hintergrund
        self.Qs = np.cumsum(self.occ - np.array([1.0 if n % 2 == 0 else 0.0 for n in range(N)]), axis=1)[:, :-1]
        self.v0 = np.ones(len(self.B)) / np.sqrt(len(self.B))
        self.blk = {}

    def diag(self, m, g2, F0, bg):
        Q = self.Qg if bg == "gleich" else self.Qs
        return m * self.stag + 0.5 * g2 * ((F0 + Q) ** 2).sum(1)

    def grund(self, m, g2, F0, bg="gleich", dicht=False, t2=1.0):
        d = self.diag(m, g2, F0, bg)
        Hm = self.hop + sps.diags(d) + (t2 - 1.0) * self.hop2
        if dicht:
            w, v = np.linalg.eigh(Hm.toarray())
            return w[:2], v[:, 0]
        w, v = eigsh(Hm, k=2, which="SA", v0=self.v0, tol=0, ncv=40, maxiter=20000)
        o = np.argsort(w)
        return w[o], v[:, o[0]]

    def bloecke(self, cut):
        if cut not in self.blk:
            mask = (1 << cut) - 1
            bl = {}
            for i, b in enumerate(self.B):
                bl.setdefault(bin(b & mask).count("1"), []).append(i)
            out = {}
            for na in sorted(bl):
                ids = bl[na]
                As = sorted(set(self.B[i] & mask for i in ids))
                Rs = sorted(set(self.B[i] >> cut for i in ids))
                ia = {a: k for k, a in enumerate(As)}
                ir = {x: k for k, x in enumerate(Rs)}
                out[na] = (np.array(ids), np.array([ia[self.B[i] & mask] for i in ids]),
                           np.array([ir[self.B[i] >> cut] for i in ids]), As, len(Rs))
            self.blk[cut] = out
        return self.blk[cut]

    def spektrum(self, psi, cut, vektoren=False):
        """Verschränkungsenergien xi = -2 ln(Schmidt-Wert) je Teilchenzahl der linken Hälfte (aufsteigend)."""
        res = {}
        for na, (ids, ia, ir, As, lr) in self.bloecke(cut).items():
            M = np.zeros((len(As), lr))
            M[ia, ir] = psi[ids]
            if vektoren:
                U, s, _ = np.linalg.svd(M, full_matrices=False)
                res[na] = (-2 * np.log(np.maximum(s, 1e-300)), U, As)
            else:
                s = np.linalg.svd(M, compute_uv=False)
                res[na] = -2 * np.log(np.maximum(s, 1e-300))
        return res

    def niveaus(self, m, g2, F0, bg="gleich", cut=None, dicht=False, t2=1.0):
        """Niedrige Vielteilchenniveaus in Einheiten von pi, relativ zum tiefsten: p (Q = +1), h (Q = -1),
        ph (zweites Niveau bei Q = 0), p1b (zweites bei Q = +1), p2 (Q = +2), h1b, Z = e^(xi_0)."""
        cut = self.N // 2 if cut is None else cut
        w, psi = self.grund(m, g2, F0, bg, dicht, t2)
        e = self.spektrum(psi, cut)
        n0 = min(sorted(e), key=lambda k: e[k][0])
        x0 = e[n0][0]
        return dict(p=(e[n0 + 1][0] - x0) / PI, h=(e[n0 - 1][0] - x0) / PI, ph=(e[n0][1] - x0) / PI,
                    p1b=(e[n0 + 1][1] - x0) / PI, p2=(e[n0 + 2][0] - x0) / PI, h1b=(e[n0 - 1][1] - x0) / PI,
                    Z=float(np.exp(x0)), luecke=float(w[1] - w[0]), n0=n0)


def freie_leiter(N, cut, m=-2.0):
    """Freie Kette: Einteilchen-Verschränkungsenergien / pi (nach Betrag sortiert) und wahrscheinlichste Teilchenzahl."""
    H = np.zeros((N, N))
    for n in range(N - 1):
        H[n, n + 1] = H[n + 1, n] = -1.0
    for n in range(N):
        H[n, n] = m if n % 2 == 0 else -m
    _, v = np.linalg.eigh(H)
    occ = v[:, :N // 2]
    z = np.clip(np.linalg.eigvalsh((occ @ occ.T)[:cut, :cut]), 1e-300, 1 - 1e-16)
    e = np.log((1 - z) / z) / PI
    e = e[np.abs(e) < 9.5]
    return sorted(e.tolist(), key=abs), int(np.sum(z > 0.5))


# ---------------------------------------------------------------------------------------------------------------------
say("Z155: Die geeichte Uhr auf der nativen Kette (Gitter-Schwinger-Modell)")
say("=" * 100)
say("\nTeil A  Form, exakt: die elektrische Energie als quadratische Form in den Ladungen")

# A1: Gauß-Gesetz eliminiert, symbolisch (N = 7 Plätze, 6 Links)
NS = 7
F = sp.Symbol("F")
q = sp.symbols(f"q0:{NS}")
Hel = sp.expand(sum((F + sum(q[:l])) ** 2 for l in range(1, NS)) / 2)
ok = sp.expand(Hel.coeff(F, 2)) == sp.Rational(NS - 1, 2)
for k in range(NS):
    ok &= sp.expand(Hel.coeff(q[k], 2)) == sp.Rational(NS - 1 - k, 2)
    ok &= sp.expand(Hel.coeff(F, 1).coeff(q[k], 1)) == NS - 1 - k
    for k2 in range(k + 1, NS):
        ok &= sp.expand(Hel.coeff(q[k], 1).coeff(q[k2], 1)) == NS - 1 - k2
check("A1 [E] Gauß-Gesetz eliminiert (symbolisch, 7 Plätze): H_el/g^2 = (1/2) sum_l (F + Q_l)^2 mit Q_l = Ladung links von "
      "Link l. Der Koeffizient von q_k q_k' ist die Zahl der Links rechts von beiden Plätzen (lineares Coulomb-Potential), "
      "der von F q_k die Zahl der Links rechts von k", ok)

# A2: zwei Gebiete A = {0,1,2}, Lücke {3}, B = {4,5}, Rest {6}: Kreuzterm
A_s, B_s = [0, 1, 2], [4, 5]
QA, QB = sum(q[k] for k in A_s), sum(q[k] for k in B_s)
kreuz = sum(Hel.coeff(q[a], 1).coeff(q[b], 1) * q[a] * q[b] for a in A_s for b in B_s)
n_aus = NS - 1 - max(B_s)
DB = sum((max(B_s) - k) * q[k] for k in B_s)
check("A2 [E] Zwei Gebiete A links von B: Der gemischte Teil ist Q_A (n_aus Q_B + D_B). Der Kreuzterm Q_A Q_B hat den "
      "Koeffizienten g^2 mal n_aus = Zahl der Links rechts von beiden Gebieten (Außenlinks); D_B ist das Dipolmoment von B "
      "um sein rechtes Ende, keine Funktion der Ladung", sp.expand(kreuz - QA * (n_aus * QB + DB)) == 0 and n_aus == 1)
zw = sp.expand(sum((F + sum(q[:l])) ** 2 for l in (3, 4)) / 2)          # die zwei Links zwischen A und B
check("A2 [E] Die Links zwischen A und B tragen (F + Q_A + Lückenladung)^2: nur Q_A, kein Q_B und kein Kreuzterm. Endet die "
      "offene Kette an B (n_aus = 0), gibt es keinen Term Q_A Q_B",
      all(sp.expand(zw.coeff(q[b], 1)) == 0 for b in B_s) and
      sp.expand(zw.subs(q[3], 0) - (F + QA) ** 2) == 0)

# A3: Ring mit Fluss (exakt, Brüche): Gauß lösbar nur bei Gesamtladung 0; dann E = E_0 + Q_l
okr = True
for n1, n2 in ((1, 1), (2, 3), (1, 4)):
    for qa in range(-3, 4):
        for e0 in (Fr(0), Fr(1, 4), Fr(-1, 4), Fr(2)):
            E = [e0 + qa] * n1 + [e0] * n2                        # Feld auf den Links von A nach B und von B zurück nach A
            okr &= (E[0] - E[-1] == qa) and (E[-1] - E[n1 - 1] == -qa)   # Gauß an A und an B = -A
            okr &= sum(x * x for x in E) / 2 == (n1 * (e0 + qa) ** 2 + n2 * e0 ** 2) / 2
check("A3 [E] Ring mit Fluss: Gauß verlangt Q_A + Q_B = 0; H_el/g^2 = (n_1 (E_0 + Q_A)^2 + n_2 E_0^2)/2 mit der "
      "Windungszahl E_0 als einziger neuer Marke. Auch hier nur eine Variable je Link, kein unabhängiger Term Q_A Q_B", okr)
befund("A Geometrie", "E",
       "In jeder der drei Geometrien (Link zwischen den Gebieten, Außenlinks, Ring) ist die Ladungsform bei fester "
       "Gesamtladung eine Funktion von Q_A allein. Der Kreuzterm Q_A Q_B kommt nur über Links, die beide Ladungen sehen, "
       "und ist dann (Q_A + Q_B)^2, also zentral für jeden ladungserhaltenden Transfer.")

# A4: Fockmodell zweier Nähte aus z145 (K = 3 Moden je Naht): In welcher Basis ist die elektrische Entwicklung das Gatter?
I2 = np.eye(2, dtype=complex)
Z2 = np.diag([1, -1]).astype(complex)


def kron(*ms):
    out = np.eye(1, dtype=complex)
    for x in ms:
        out = np.kron(out, x)
    return out


def fock(n):
    low = np.array([[0, 1], [0, 0]], complex)
    return [kron(*([Z2] * k + [low] + [I2] * (n - k - 1))) for k in range(n)]


def schmidtrang(Op, d):
    T = Op.reshape(d, d, d, d).transpose(0, 2, 1, 3).reshape(d * d, d * d)
    return int(np.sum(np.linalg.svd(T, compute_uv=False) > 1e-9))


KM = 3
cf = fock(2 * KM)
cA, cB = cf[:KM], cf[KM:]
NA = sum(dag(x) @ x for x in cA)
NB = sum(dag(x) @ x for x in cB)
nA, nB = np.rint(np.diag(NA).real), np.rint(np.diag(NB).real)
nT = nA + nB
Jp = sum(dag(a) @ b for a, b in zip(cA, cB))
Jx, Jy, Jz = (Jp + dag(Jp)) / 2, (Jp - dag(Jp)) / 2j, (NA - NB) / 2
Bq = expm(-1j * np.pi / 2 * Jx)
dim = len(nA)


def Uz(phi):
    """z120/z145-Gatter U(phi) = B^+ e^(-i phi (Q_A^2 + Q_B^2)/2) B."""
    return dag(Bq) @ np.diag(np.exp(-1j * phi * (nA ** 2 + nB ** 2) / 2)) @ Bq


res = {}
res["zwei Links"] = max(
    nrm(dag(Bq) @ np.diag(np.exp(-1j * t * ((Fv + nA) ** 2 + (Fv + nB) ** 2) / 2)) @ Bq
        - np.diag(np.exp(-1j * t * (Fv * nT + Fv * Fv))) @ Uz(t))
    for Fv in (0.0, -0.25, 0.25, 0.37) for t in (np.pi / 2, 0.7))
check("A4 [E] Zwei Nähte, jede mit eigenem Schnittlink (Feld F + Q_A bzw. F + Q_B): B^+ e^(-i t H_el) B = "
      "e^(-i g^2 t (F N + F^2)) U(phi = g^2 t) für jedes F. Die elektrische Entwicklung ist die Phase des z120-Gatters; "
      "das Hintergrundfeld geht nur in eine zentrale Phase der Gesamtladung", res["zwei Links"] < TOL)
rang_naht = schmidtrang(np.diag(np.exp(-1j * np.pi / 2 * (nA ** 2 + nB ** 2) / 2)), 2 ** KM)
rang_gatter = schmidtrang(Uz(np.pi / 2), 2 ** KM)
ZAHLEN["Operator-Schmidt-Rang: elektrische Phase in Nahtbasis, Gatter"] = [rang_naht, rang_gatter]
check("A4 [E] In der Nahtbasis ist e^(-i phi H_Q) ein Produkt von Einnahtphasen (Operator-Schmidt-Rang 1): Die geeichte "
      "Nahtladung allein gibt das triviale Gatter (wie z121, z145). Rang des Gatters mit B-Sandwich: "
      f"{rang_gatter}", rang_naht == 1 and rang_gatter > 1)
Gy, Gz = expm(1j * np.pi * Jy), expm(1j * np.pi * Jz)
gate = Uz(np.pi / 2)
res["zirkular"] = nrm(gate - expm(-1j * np.pi / 8 * (NA + NB) @ (NA + NB)) @ expm(-1j * np.pi / 2 * Jy @ Jy))
res["Paaruhr"] = nrm(gate @ Gz - Gz @ gate) + nrm(gate @ Gy - Gy @ gate)
check("A4 [E] Das Gatter ist elektrisch in den zirkularen Moden: U(phi) = e^(-i phi N^2/4) e^(-i phi J_y^2), "
      "J_y = (Q_+ - Q_-)/2. Geeicht ist dabei die Ladung der Vierteldrehung P_sym = e^(i pi J_y), nicht die Nahtladung; "
      "das Gatter vertauscht mit P_sym und mit der Paaruhr e^(i pi J_z) der Nähte (P_univ)",
      res["zirkular"] < TOL and res["Paaruhr"] < TOL)
# ein gemeinsamer Link: symmetrisches Feld F = -N/2 gegen F = -1/4
W_sym = np.diag(np.exp(-1j * np.pi * (nA - nT / 2) ** 2 / 2))                # t g^2 = pi
res["ein Link sym"] = nrm(dag(Bq) @ W_sym @ Bq - np.diag(np.exp(1j * np.pi / 2 * nT ** 2 / 4)) @ gate)
W_q = np.diag(np.exp(-1j * np.pi * (nA - 0.25) ** 2 / 2))
Rr = dag(Bq) @ W_q @ Bq @ dag(gate)
sek = np.diag((nT == 2).astype(complex))
abst = min(1 - abs(np.trace(sek @ np.linalg.matrix_power(dag(Gy), k) @ Rr)) / np.trace(sek).real for k in range(4))
ZAHLEN["ein gemeinsamer Link, F = -1/4: Abstand der Restdrehung von der Uhrgruppe (N = 2)"] = r(abst, 6)
check("A4 [E] Ein gemeinsamer Link zwischen A und B: B^+ e^(-i t g^2 (F + Q_A)^2/2) B ist das Gatter (phi = g^2 t/2, Faktor 2 "
      "gegenüber zwei Links) genau bei symmetrischem Feld F = -N/2. Mit F = -1/4 bleibt die einkörperige Drehung "
      "e^(-i phi (N - 1/2) J_y), die nicht in der Uhrgruppe liegt",
      res["ein Link sym"] < TOL and abst > 0.05)
O_aus = np.diag(np.exp(-1j * 0.7 * (0.25 + nT) ** 2 / 2))
check("A4 [E] Außenlink (F + Q_A + Q_B)^2: vertauscht mit dem Transfer B und ist in jedem Sektor fester Gesamtladung eine "
      "Zahl. Der Kreuzterm Q_A Q_B aus dem Gauß-Gesetz wechselwirkt also nicht",
      nrm(O_aus @ Bq - Bq @ O_aus) < TOL and schmidtrang(O_aus, 2 ** KM) > 1)

# A5: Rest. Kette N = 12 bei g = 0, Naht = linke Hälfte
K12 = Kette(12)
w0, psi0 = K12.grund(-2.0, 0.0, -0.25, dicht=True)
cut = 6
mask = (1 << cut) - 1
QA_c = K12.Qg[:, cut - 1]
R_in = 0.5 * ((-0.25 + K12.Qg[:, :cut - 1]) ** 2).sum(1)                   # innere Links der Naht
E_cut = 0.5 * (-0.25 + QA_c) ** 2
streu = {}
for qv in (-1.0, 0.0, 1.0):
    sel = np.abs(QA_c - qv) < 1e-9
    pw = psi0[sel] ** 2 / np.sum(psi0[sel] ** 2)
    mw = float(np.sum(pw * R_in[sel]))
    streu[int(qv)] = [r(mw, 5), r(np.sqrt(np.sum(pw * (R_in[sel] - mw) ** 2)), 5),
                      r(np.ptp(E_cut[sel]), 12), r(np.ptp(R_in[sel]), 5)]
ZAHLEN["Rest: Q_A -> [Mittel, Streuung, Spanne Schnittlink, Spanne innere Links] der elektrischen Energie / g^2"] = streu
say(f"   Naht = linke Hälfte (6 Plätze), g = 0: Q_A -> [Mittel, Streuung der inneren Links; Spanne Schnittlink; Spanne innen] = {streu}")
check("A5 [E] Schnittlink: (F + Q_A)^2/2 ist bei fester Nahtladung eine Zahl (Spanne 0). Das ist P_null für den Link "
      "zwischen Naht und Partner, eine Folge des Gauß-Gesetzes", all(v[2] == 0.0 for v in streu.values()))
check("A5 [N] Innere Links der Naht: Ihre elektrische Energie ist bei fester Nahtladung keine Zahl (Spanne und Streuung "
      "im nativen Zustand > 0). Der Rest ist nur null, wenn das Feld allein auf den Schnittlinks lebt (P_schnitt)",
      all(v[1] > 1e-3 and v[3] > 0.5 for v in streu.values()))
befund("A Form", "bedingt",
       "Teil A: Die elektrische Energie des Schnittlinks ist exakt eine Funktion der Nahtladung ohne Oszillatoranteil "
       "(P_null folgt aus P_eich, wenn das Feld nur auf den Schnittlinks lebt: P_schnitt). Zwei Nähte mit je einem "
       "Schnittlink geben genau H_Q = (Q_A^2 + Q_B^2)/2; mit dem Transfer B davor und danach ist das das z120-Gatter mit "
       "phi = g^2 t. Die Wechselwirkung kommt nicht aus dem Kreuzterm des Gauß-Gesetzes, sondern aus dem Basiswechsel B "
       "(P_sym); den liefert die Eichung nicht.")

# ---------------------------------------------------------------------------------------------------------------------
say("\nTeil B  Hintergrundfeld und das Viertel")
# B1: native Ladungsparabel je Schnitt (freie Kette), gerade und ungerade Schnitte
lad = {}
for c_ in (99, 100, 101):
    e, nstar = freie_leiter(200, c_)
    lad[c_] = (e[:4], nstar)
    say(f"   freie Kette (200 Plätze), Schnitt nach {c_} Plätzen: Niveaus / pi = {[r(x, 5) for x in e[:4]]}, "
        f"wahrscheinlichste Teilchenzahl {nstar}")
soll_g = [1, -3, 5, -7]
okl = all(len(lad[c_][0]) == 4 for c_ in (99, 100, 101))
okl &= max(abs(a - b) for a, b in zip(lad[100][0], soll_g)) < 1e-5
okl &= all(max(abs(a + b) for a, b in zip(lad[c_][0], soll_g)) < 1e-5 for c_ in (99, 101))
okl &= lad[100][1] == 50 and lad[99][1] == 50 and lad[101][1] == 51
check("B1 [N] Gerader Schnitt: Leiter pi (4Z + 1) (Naht), wahrscheinlichste Teilchenzahl j/2. Ungerader Schnitt: "
      "gespiegelte Leiter pi (4Z + 3) (Partner), Teilchenzahl (j + 1)/2 (vier Niveaus, 1e-5)", okl)


def parabel(teil, loch, Q):
    """Tiefste Energie bei Ladung Q in Einheiten 4 pi: Pauli-Summe über die Leiter."""
    return sum(teil + j for j in range(Q)) if Q >= 0 else sum(loch + j for j in range(-Q))


okp = True
for Q in range(-8, 9):
    xg = Fr(Q)                                                   # gerader Schnitt: Q_gleich = Q
    xu = Fr(Q) + Fr(1, 2)                                        # ungerader Schnitt: Q_gleich = Q + 1/2
    okp &= parabel(Fr(1, 4), Fr(3, 4), Q) == (xg - Fr(1, 4)) ** 2 / 2 - Fr(1, 32)
    okp &= parabel(Fr(3, 4), Fr(1, 4), Q) == (xu - Fr(1, 4)) ** 2 / 2 - Fr(1, 32)
    okp &= (Fr(Q) - Fr(1, 4)) ** 2 / 2 - Fr(1, 32) == Fr(Q * Q, 2) - Fr(Q, 4)
    okp &= (Fr(Q) + Fr(1, 4)) ** 2 / 2 - Fr(1, 32) == Fr(Q * Q, 2) + Fr(Q, 4)
check("B1 [E] Pauli-Summe über die Leiter (|Q| <= 8): tiefstes Niveau bei Ladung Q ist 4 pi ((Q_gl - 1/4)^2/2 - 1/32) mit "
      "Q_gl = Ladung gegen den gleichmäßigen Hintergrund (n - 1/2), an geraden und an ungeraden Schnitten. Das sind die "
      "zwei z113-Formeln (Q - 1/4)^2/2 - 1/32 und (Q + 1/4)^2/2 - 1/32 in einer Variablen", okp)
okd = True
tab = {}
for c_ in (5, 6, 7):
    e = K12.spektrum(psi0, c_)
    n0 = min(sorted(e), key=lambda k: e[k][0])
    row = []
    for dn in (-1, 0, 1):
        x = n0 + dn - c_ / 2
        soll = ((x - 0.25) ** 2 - (n0 - c_ / 2 - 0.25) ** 2) / 2
        ist = (e[n0 + dn][0] - e[n0][0]) / (4 * PI)
        row.append([r(x, 2), r(ist, 5), r(soll, 5)])
        okd &= abs(ist - soll) < 5e-4
    tab[c_] = row
ZAHLEN["ED N = 12, g = 0: Schnitt -> [Q_gl, xi_min/4pi gemessen, Parabel]"] = tab
check("B1 [N] Exakte Diagonalisierung (12 Plätze, g = 0, Schnitte nach 5, 6, 7 Plätzen): tiefstes Niveau je Ladung auf der "
      "Parabel mit Scheitel Q_gl = +1/4 (5e-4)", okd)
_, psi_p = K12.grund(2.0, 0.0, 0.25, dicht=True)
e = K12.spektrum(psi_p, 6)
check("B1 [N] Umgekehrtes Vorzeichen des Wechselpotentials (m = +2t): Scheitel bei Q_gl = -1/4. Das Vorzeichen des "
      "Viertels folgt dem von m (z152: m = -+2t)",
      abs((e[4][0] - e[3][0]) / PI - 3) < 2e-3 and abs((e[2][0] - e[3][0]) / PI - 1) < 2e-3)

# B2: z113-Formel, Konstante, Uhr
okc = True
for Q in range(-8, 9):
    hel = (Fr(Q) - Fr(1, 4)) ** 2 / 2                              # H_el / g^2 auf dem Schnittlink mit F0 = -1/4
    L0 = Fr(Q * Q, 2) - Fr(Q, 4)
    okc &= hel - L0 == Fr(1, 32)
    okc &= hel == Fr((4 * Q - 1) ** 2, 32)                         # Gewicht q(n) = n^2/32 des Sektors n = 4Q - 1 (z145)
    okc &= (L0 - Fr(Q, 4)).denominator == 1                        # e^(2 pi i L_0) = i^Q
check("B2 [E] Mit F0 = -1/4 ist H_el/g^2 auf dem Schnittlink (Q - 1/4)^2/2 = L_0^Ladung + 1/32 = (4Q - 1)^2/32, das "
      "Gewicht des Z_16-Sektors n = 4Q - 1 aus z145. Die Form der z113-Formel ist reproduziert, die Konstante nicht: "
      "H_el >= g^2/32 (Feldenergie F^2/2 des Hintergrunds), z113 hat 0 bei Q = 0. Bei g^2 t = 2 pi gilt "
      "e^(-i t H_el) = e^(-2 pi i/32) R^-1: die Uhr bis auf die Twistkonstante aus z145", okc)

# B3: gleichmäßiger Hintergrund mit F0 = gestaffelter Hintergrund mit theta/2pi = F0 + 1/4 plus Massenverschiebung g^2/8
okq = True
for N in (8, 10):
    for b in [b for b in range(1 << N) if bin(b).count("1") == N // 2]:
        o = [(b >> n) & 1 for n in range(N)]
        for F0 in (Fr(-1, 4), Fr(0), Fr(1, 4), Fr(1, 3)):
            Qg = Qs = Fr(0)
            hg = hs = Fr(0)
            for n in range(N - 1):
                Qg += o[n] - Fr(1, 2)
                Qs += o[n] - (1 if n % 2 == 0 else 0)
                hg += (F0 + Qg) ** 2 / 2
                hs += (F0 + Fr(1, 4) + Qs) ** 2 / 2
            st = sum((-1) ** n * o[n] for n in range(N))
            okq &= hg == hs + Fr(st, 8) + (F0 + Fr(1, 4)) / 4 - Fr(N, 16) + Fr(N - 1, 32)
check("B3 [E] Identität auf allen neutralen Besetzungen (8 und 10 Plätze, vier Werte F0): gleichmäßiger Hintergrund mit "
      "Randfeld F0 = Kogut-Susskind-Hintergrund mit theta/2pi = F0 + 1/4, plus (g^2/8) sum (-1)^n n_n, plus Konstante. "
      "F0 = -1/4 ist also theta = 0 mit der Massenverschiebung g^2 a/8, F0 = +1/4 ist theta = pi; ein gleichmäßiges "
      "Feld theta = pi/2 wäre F0 = 0", okq)
K10 = Kette(10)
wa = np.linalg.eigvalsh((K10.hop + sps.diags(K10.diag(-2.0, 1.3, -0.25, "gleich"))).toarray())
wb = np.linalg.eigvalsh((K10.hop + sps.diags(K10.diag(-2.0 + 1.3 / 8, 1.3, 0.0, "gestaffelt"))).toarray())
check("B3 [N] Gegenprobe im Spektrum (10 Plätze, g^2 = 1,3): H(gleichmäßig, F0 = -1/4, m) und H(gestaffelt, theta = 0, "
      "m + g^2/8) haben dasselbe Spektrum bis auf eine Konstante (1e-10)",
      float(np.max(np.abs((wa - wa[0]) - (wb - wb[0])))) < 1e-10)

# B4: Welche F0 sind mit Ladungskonjugation mal Verschiebung um einen Platz verträglich?
ct = []
for k in range(-8, 9):
    F0 = Fr(k, 16)
    gerade = F0 % 1                                 # Feldwerte auf geraden Links: F0 + Z, auf ungeraden: F0 + 1/2 + Z
    ungerade = (F0 + Fr(1, 2)) % 1
    if (-gerade) % 1 == ungerade:                   # C: E -> -E, T: gerade <-> ungerade
        ct.append(str(F0))
ZAHLEN["CT-symmetrische Randfelder F0 in (1/16) Z, |F0| <= 1/2"] = ct
check("B4 [E] Das Gauß-Gesetz mit gleichmäßigem Hintergrund ist unter Ladungskonjugation mal Verschiebung um einen Platz "
      "(E_l -> -E_(l+1)) genau für F0 = +-1/4 mod 1 invariant: zwei Klassen (theta = 0 und theta = pi), kein Kontinuum",
      ct == ["-1/4", "1/4"])
feld = sorted(set(np.round(-0.25 + K12.Qg[np.argmax(np.abs(psi0))], 6).tolist()))
check("B4 [E] Native Klasse: Im Zustand größter Amplitude (bevorzugte Plätze besetzt) trägt mit F0 = -1/4 jeder Link das "
      "Feld -+1/4 (kleinste Feldenergie, theta = 0). Der Scheitel +1/4 der nativen Parabel (B1) ist dieses F0",
      feld == [-0.25, 0.25])
befund("B Viertel", "bedingt",
       "Teil B: F = -+1/4 ist auf der Kette kein gleichmäßiges Hintergrundfeld theta = pi/2. Es ist der Versatz zwischen "
       "halber Füllung (gleichmäßiger Hintergrund) und Dirac-See, dieselbe Zahl wie der Leiterversatz 1/4 aus z152, und "
       "in der Standardlesart theta = 0 mit der bekannten Massenverschiebung g^2 a/8. Als Parameter des Gauß-Gesetzes ist "
       "F0 trotzdem eine Eingabe; die Symmetrie CT legt sie bis auf ein Bit (theta = 0 oder pi) fest (P_CT).")

# B5: reine Z_4-Eichtheorie: Energiefunktionen auf Z_4
okz = all(((Fr(Q) - Fr(1, 4)) ** 2 / 2) % 1 == (Fr(1, 32) + Fr(Q, 4)) % 1 for Q in range(-12, 13))
check("B5 [E] Modulo 1 ist (Q - 1/4)^2/2 = 1/32 + Q/4 eine Funktion auf Z_4 und linear: Auf Z_4 wohldefiniert ist von der "
      "quadratischen Energie mit Viertelversatz nur die Entwicklung über g^2 t = 2 pi, und die ist die Uhr", okz)
wohl = []
for j in range(0, 17):                                           # tau = j pi/8
    ph = [(-Fr(j, 16) * ((Fr(Q) - Fr(1, 4)) ** 2 / 2)) % 1 for Q in range(-6, 7)]     # Phase / 2 pi
    d = set((ph[i + 4] - ph[i]) % 1 for i in range(len(ph) - 4))
    ph0 = [(-Fr(j, 16) * Fr(Q * Q, 2)) % 1 for Q in range(-6, 7)]
    d0 = set((ph0[i + 4] - ph0[i]) % 1 for i in range(len(ph0) - 4))
    if len(d) == 1:
        wohl.append((j, str(next(iter(d))), d0 == {Fr(0)}))
ZAHLEN["quadratische Phase e^(-i tau E^2/2) hängt nur von E mod 4 ab: tau = j pi/8 -> (j, Phase je Umlauf / 2 pi, ohne Versatz strikt)"] = wohl
check("B5 [E] Die quadratische Phase e^(-i tau E^2/2) hängt genau für tau in (pi/2) Z nur von E mod 4 ab (ohne Versatz "
      "strikt, mit Versatz 1/4 bis auf eine Konstante i^k je Umlauf). Für eine Z_4-Eichtheorie mit quadratischem "
      "Feldgesetz ist die Schrittphase also auf acht Werte quantisiert: die acht Uhrgatter aus z145",
      [w[0] for w in wohl] == [0, 4, 8, 12, 16] and all(w[2] for w in wohl))
# alle Funktionen g: Z_4 -> 16. Einheitswurzeln mit g(0) = 1
quad, nichtentartet, clifford, uhr4 = 0, 0, 0, 0
gleich = True
for a1 in range(16):
    for a2 in range(16):
        for a3 in range(16):
            g = [0, a1, a2, a3]                                   # Exponenten: g(Q) = e^(2 pi i g[Q]/16)
            b = [[(g[(x + y) % 4] - g[x] - g[y]) % 16 for y in range(4)] for x in range(4)]
            bich = all(b[(x + x2) % 4][y] == (b[x][y] + b[x2][y]) % 16 for x in range(4) for x2 in range(4) for y in range(4))
            h = [(g[(E - 1) % 4] - g[E]) % 16 for E in range(4)]   # Transfer X -> X mal h(E)
            cl = all((h[(E + 1) % 4] - h[E]) % 16 == (h[1] - h[0]) % 16 for E in range(4)) and (4 * (h[1] - h[0])) % 16 == 0
            gleich &= bich == cl
            if bich:
                quad += 1
                ne = b[1][1] % 8 == 4                              # b(1,1) = +-i: Paarung nichtentartet
                nichtentartet += ne
                vier = [(4 * x) % 16 for x in g]                   # vierte Potenz
                if ne and any(all((vier[Q] - s * 4 * Q) % 16 == 0 for Q in range(4)) for s in (1, -1)):
                    uhr4 += 1
ZAHLEN["Funktionen auf Z_4 (16. Einheitswurzeln): alle, quadratisch, nichtentartet"] = [4096, quad, nichtentartet]
check("B5 [E] Von 4096 Phasenfunktionen auf Z_4 sind 16 quadratisch und 8 nichtentartet (vier Formen mal Charakter). "
      "Die Clifford-Bedingung (Transfer -> Transfer mal ein Uhrelement, z145) gilt genau für die quadratischen",
      quad == 16 and nichtentartet == 8 and gleich)
check("B5 [E] Die vierte Potenz eines nichtentarteten Gatters ist die Parität R^2, nie die Uhr R: Gatter (quadratisch) und "
      "Uhr (linear, vom Viertelversatz) sind nicht Viertel- und Vollperiode derselben Energiefunktion auf Z_4 "
      "(wie z133: Gauß-Summe gegen lineare Uhr)", uhr4 == 0)
befund("B Z_4", "E",
       "Reine Z_4-Eichtheorie: Jede elektrische Phase ist von selbst eine Funktion der Uhrwerte ohne Oszillatoranteil, "
       "denn E_l mod 4 ist der Uhrwert des Schnitts. Die Energiefunktion hat aber drei freie Phasen (zwei mit "
       "Ladungskonjugation) statt einer; die wechselwirkenden Uhrgatter sind darin acht isolierte Punkte. Quadratisch ist "
       "eine eigene Prämisse (P_quad).")

# ---------------------------------------------------------------------------------------------------------------------
say("\nTeil C  Der Winkel phi = g^2 mal Zeit")
# C1: Clifford-Bedingung auf der Kette (N = 8): Hüpfen über Link l wird mit der Uhr des Schnitts l gekleidet
K8 = Kette(8)
nb = len(K8.B)
E8 = -0.25 + K8.Qg


def huepf(l):
    X = np.zeros((nb, nb), complex)
    for i, b in enumerate(K8.B):
        if (b >> (l - 1)) & 1 and not (b >> l) & 1:
            X[K8.idx[b ^ (3 << (l - 1))], i] = 1.0                  # Teilchen von Platz l-1 nach l
    return X


def uhr(l):
    qi = np.rint(K8.Qg[:, l - 1] - (0.5 if l % 2 else 0.0)).astype(int)    # ganzzahlige Ladung des Schnitts
    return np.diag(1j ** (qi % 4))


okcl, rest_min = True, 1.0
for l in (3, 4):
    X, Rl = huepf(l), uhr(l)
    kand = [(X @ np.linalg.matrix_power(Rl, k)).ravel() for k in range(4)]
    for k in range(8):
        U = np.diag(np.exp(-1j * (k * np.pi / 2) * (E8 ** 2).sum(1) / 2))
        M = U @ X @ dag(U)
        Y = X @ np.linalg.matrix_power(Rl, k % 4)
        cph = np.vdot(Y.ravel(), M.ravel()) / np.vdot(Y.ravel(), Y.ravel())
        okcl &= nrm(M - cph * Y) < TOL and abs(abs(cph) - 1) < TOL
    for tau in (0.3, 0.7, np.pi / 4, 1.0, 2.0):
        U = np.diag(np.exp(-1j * tau * (E8 ** 2).sum(1) / 2))
        v = (U @ X @ dag(U)).ravel()
        for y in kand:                                              # bester Abstand zu Konstante mal X R^j
            cj = np.vdot(y, v) / np.vdot(y, y)
            rest_min = min(rest_min, float(np.linalg.norm(v - cj * y) / np.linalg.norm(v)))
ZAHLEN["Clifford auf der Kette: kleinster relativer Abstand zu einem Uhrelement außerhalb (pi/2) Z"] = r(rest_min, 5)
check("C1 [E] Auf der geeichten Kette (8 Plätze, F0 = -1/4) gilt e^(-i t H_el) X_l e^(i t H_el) = Konstante mal X_l R_l^k "
      "(ein einzelnes Uhrelement) genau für g^2 t = k pi/2 (X_l = Hüpfen über Link l, R_l = Uhr des Schnitts l; Kontrollen "
      "0,3; 0,7; pi/4; 1; 2). "
      "P_viertel heißt in Eichsprache g^2 t = pi/2 mod pi: 'Ladung um Fluss' aus z145", okcl and rest_min > 0.05)

# C2: Modulperiode und Boost-Gewicht (N = 12, erste Ordnung in g^2)
eps = 1e-3


def stufen(g2):
    _, ps = K12.grund(-2.0, g2, -0.25, dicht=True)
    return K12.spektrum(ps, 6, vektoren=True)


s0, sp_, sm_ = stufen(0.0), stufen(eps), stufen(-eps)
marken = [(4, 0, "pi, Q=+1"), (2, 0, "3pi, Q=-1"), (3, 1, "4pi, Q=0"), (4, 1, "5pi, Q=+1"), (5, 0, "6pi, Q=+2"),
          (2, 1, "7pi, Q=-1")]
c_lem = float(mp.gamma(mp.mpf(1) / 4) ** 2 / (2 * mp.sqrt(2 * mp.pi)))
occA = {a: np.array([(a >> k) & 1 for k in range(6)], float) for a in range(64)}


def erwartung(na, i, f):
    xi, U, As = s0[na]
    return float(sum(U[j, i] ** 2 * f(occA[a]) for j, a in enumerate(As)))


def f_boost(o):
    Q = np.cumsum(o - 0.5)
    return sum(c_lem * (6 - l) * 0.5 * (-0.25 + Q[l - 1]) ** 2 for l in range(1, 6))


def f_schnitt(o):
    return 0.5 * (-0.25 + np.sum(o - 0.5)) ** 2


steig, boost, schn, lev0 = [], [], [], []
b0, s_0 = erwartung(3, 0, f_boost), erwartung(3, 0, f_schnitt)
d0 = (sp_[3][0][0] - sm_[3][0][0]) / (2 * eps)
for na, i, nm in marken:
    steig.append((sp_[na][0][i] - sm_[na][0][i]) / (2 * eps) - d0)
    boost.append(erwartung(na, i, f_boost) - b0)
    schn.append(erwartung(na, i, f_schnitt) - s_0)
    lev0.append((s0[na][0][i] - s0[3][0][0]) / PI)
steig, boost, schn = np.array(steig), np.array(boost), np.array(schn)
ZAHLEN["erste Ordnung: Niveau/pi bei g = 0"] = [r(x, 4) for x in lev0]
ZAHLEN["erste Ordnung: d xi / d g^2 (exakt)"] = [r(x, 4) for x in steig]
ZAHLEN["erste Ordnung: Erwartung c sum d E_d^2/2 (Boost-Gewicht, Schnittlink mit Gewicht 0)"] = [r(x, 4) for x in boost]
ZAHLEN["erste Ordnung: Erwartung (F + Q_A)^2/2 (nur Schnittlink)"] = [r(x, 4) for x in schn]
say(f"   Niveaus / pi bei g = 0:              {[r(x, 4) for x in lev0]}  ({', '.join(m[2] for m in marken)})")
say(f"   d xi / d g^2, exakt (12 Plätze):     {[r(x, 4) for x in steig]}")
say(f"   Boost-gewichteter elektrischer Term: {[r(x, 4) for x in boost]}")
say(f"   nur Schnittlink (F + Q_A)^2/2:       {[r(x, 4) for x in schn]}")
komm = 0.0
for g2 in (0.0, 1.0, 5.0):
    _, ps = K12.grund(-2.0, g2, -0.25, dicht=True)
    Mfull = np.zeros((64, 64))
    for i, b in enumerate(K12.B):
        Mfull[b & 63, b >> 6] = ps[i]
    rhoA = Mfull @ Mfull.T
    RA = np.diag(1j ** np.array([bin(a).count("1") % 4 for a in range(64)]))
    komm = max(komm, nrm(rhoA @ RA - RA @ rhoA))
check("C2 [E] Der Modulfluss vertauscht für jedes g^2 mit der Uhr (||[rho_A, R_A]|| < 1e-12 bei g^2 = 0, 1, 5), weil die "
      "Nahtladung erhalten ist. Das quantisiert g^2 nicht", komm < 1e-12)
ng1 = steig[2] - steig[0] - steig[1]
ZAHLEN["erste Ordnung: d/dg^2 von (xi_ph - xi_p - xi_h) und von (xi_h - 3 xi_p)"] = [r(ng1, 4), r(steig[1] - 3 * steig[0], 4)]
check("C2 [N] Die native Gleichung 'Modulfluss über die halbe Periode = Uhr' (e^(iK/2) = Konstante mal R, also Niveaus "
      "= pi Q mod 4 pi) bricht in erster Ordnung in g^2: d xi/d g^2 = 0,093; 0,696; 0,578 für die Niveaus pi, 3 pi, 4 pi, "
      "nicht proportional zu den Niveaus und nicht additiv (Teilchen plus Loch ungleich Paar)",
      abs(steig[0] - 0.0926) < 2e-3 and abs(steig[1] - 0.696) < 3e-3 and abs(steig[2] - 0.578) < 3e-3 and ng1 < -0.15)
Afit = np.array([boost, schn]).T
xf = np.linalg.lstsq(Afit, steig, rcond=None)[0]
rest_bs = float(np.linalg.norm(steig - Afit @ xf) / np.linalg.norm(steig))
xs = np.linalg.lstsq(schn[:, None], steig, rcond=None)[0]
rest_s = float(np.linalg.norm(steig - schn * xs[0]) / np.linalg.norm(steig))
korr = float(np.corrcoef(steig, boost)[0, 1])
ZAHLEN["Anpassung d xi/d g^2 = alpha Boost + gamma Schnitt: alpha, gamma, relativer Rest"] = [r(xf[0], 3), r(xf[1], 3), r(rest_bs, 3)]
ZAHLEN["Anpassung nur Schnittlink: relativer Rest; Korrelation mit Boost-Term"] = [r(rest_s, 3), r(korr, 3)]
say(f"   Anpassung alpha Boost + gamma Schnitt: alpha = {xf[0]:.3f}, gamma = {xf[1]:.3f}, Rest {rest_bs:.3f}; "
    f"nur Schnittlink: Rest {rest_s:.3f}; Korrelation mit Boost-Term {korr:.3f}")
check("C2 [N] Boost-Gewicht: Die Antwort der Niveaus folgt dem mit dem Abstand gewichteten elektrischen Term (Korrelation "
      "> 0,95), aber nicht exakt (Faktor alpha zwischen 0,5 und 0,8). Der Schnittlink selbst trägt wenig Gewicht "
      "(gamma unter einem Fünftel von alpha c, dem Gewicht des ersten inneren Links); 'nur Schnittlink' scheitert "
      "(Rest > 0,3). Die Ladungsenergie (Q - 1/4)^2 des Schnitts geht also kaum in den eigenen Modulgenerator ein; die "
      "Modulperiode sieht g^2 vor allem über das Innere und quantisiert es nicht",
      korr > 0.95 and 0.5 < xf[0] < 0.8 and 0 <= xf[1] < 0.2 * xf[0] * c_lem and rest_s > 0.3)

# C3: Z_8 der Uhrgatter und eine native Zeit
verh = mp.sqrt(8) / 2                                             # E_max / E_min der freien Kette
check("C3 [E] Die freie Kette hat keine eigene Periode: E_max/E_min = sqrt(2) ist irrational (z152: kappa = 1/sqrt(2)). "
      "Die Laufzeit T des elektrischen Schritts ist eine zweite neue Größe; fest ist durch P_viertel nur das Produkt "
      "g^2 T = pi/2 (ein gemeinsamer Link: pi). Lesart ohne Stütze: gleiche Schrittdauer für Transfer (2 t T = pi/2) und "
      "Feld gäbe g^2 = 2t, also g a = 1 neben m a = 1", abs(verh ** 2 - 2) < mp.mpf(10) ** -50)


# C4: Selbstdualität: Villain-Gewicht auf Z_4
def w_villain(a, k):
    return mp.nsum(lambda n: mp.exp(-a * (k - 4 * n) ** 2 / 2), [-mp.inf, mp.inf])


def dual_defekt(a):
    w = [w_villain(a, k) for k in range(4)]
    wd = [mp.fsum(w[k] * mp.cos(mp.pi * j * k / 2) for k in range(4)) for j in range(4)]
    return wd[1] / wd[0] - w[1] / w[0], wd[2] / wd[0] - w[2] / w[0], w


a_sd = mp.findroot(lambda a: dual_defekt(a)[0], mp.mpf("1.3"))
d1, d2, wsd = dual_defekt(mp.pi / 2)
vorz = [mp.sign(dual_defekt(mp.mpf(x) / 10)[0]) for x in range(3, 60, 3)]
ZAHLEN["selbstduales Villain-Gewicht auf Z_4: a"] = mp.nstr(a_sd, 30)
ZAHLEN["dort w(1)/w(0), w(2)/w(0)"] = [mp.nstr(wsd[1] / wsd[0], 15), mp.nstr(wsd[2] / wsd[0], 15)]
check("C4 [E] Selbstdualität: Das uhrperiodisierte Gewicht w_a(E) = sum_n e^(-a (E - 4n)^2/2) ist genau für a = 2 pi/4 = "
      "pi/2 seine eigene Z_4-Fourier-Transformierte (40 Stellen; ein Vorzeichenwechsel auf 0,3 <= a <= 5,7). "
      "pi/2 = beta/8 mit beta = 4 pi: das KMS-Gewicht des Ladungsgitters über ein Achtel der Periode",
      abs(a_sd - mp.pi / 2) < mp.mpf(10) ** -40 and abs(d1) < mp.mpf(10) ** -40 and abs(d2) < mp.mpf(10) ** -40
      and sum(1 for i in range(len(vorz) - 1) if vorz[i] != vorz[i + 1]) == 1)


def theta2(k, tau):
    """Thetafunktion der Stufe 2: sum_n e^(2 pi i tau 2 (n - k/4)^2), Gewichte k^2/8 auf Z_4."""
    return mp.fsum(mp.exp(2j * mp.pi * tau * 2 * (n - mp.mpf(k) / 4) ** 2) for n in range(-40, 41))


tau0 = mp.mpc("0.3", "1.1")
okth = all(abs(theta2(k, 1j) - wsd[k]) < mp.mpf(10) ** -40 for k in range(4))
okth &= all(abs(theta2(k, tau0 + 1) - mp.exp(2j * mp.pi * mp.mpf(k * k) / 8) * theta2(k, tau0)) < mp.mpf(10) ** -40
            for k in range(4))
okth &= all(abs(theta2(k, -1 / tau0) - mp.sqrt(-1j * tau0) / 2 * mp.fsum(mp.exp(-2j * mp.pi * j * k / 4) * theta2(j, tau0)
                                                                       for j in range(4))) < mp.mpf(10) ** -40 for k in range(4))
okth &= all(abs(mp.exp(-1j * mp.pi / 2 * mp.mpf(k * k) / 2) * mp.exp(2j * mp.pi * mp.mpf(k * k) / 8) - 1) < mp.mpf(10) ** -50
            for k in range(-6, 7))
check("C4 [E] Einordnung: w_(pi/2)(k) ist die Thetafunktion der Stufe 2 mit Gewichten k^2/8 bei tau = i. S wirkt auf diesen "
      "vier Funktionen als Z_4-Fourier-Transformation, tau = i ist ihr Fixpunkt. Das Gatter e^(-i (pi/2) Q^2/2) ist die "
      "inverse T-Matrix (tau -> tau + 1) derselben Funktionen (40 Stellen). Die native Selbstdualität tau = i legt also "
      "den euklidischen Schritt fest; dass der reelle Schritt ein ganzer Dehn-Twist ist, ist P_viertel", okth)
x1, x2 = wsd[1] / wsd[0], wsd[2] / wsd[0]
okat = abs(2 * x1 + x2 - 1) < mp.mpf(10) ** -40
for x1t in (mp.mpf("0.2"), mp.mpf("0.35"), mp.sqrt(2) - 1):
    x2t = 1 - 2 * x1t
    okat &= abs((1 - x2t) / (1 + 2 * x1t + x2t) - x1t) < mp.mpf(10) ** -40
    okat &= abs((1 - 2 * x1t + x2t) / (1 + 2 * x1t + x2t) - x2t) < mp.mpf(10) ** -40
okat &= abs((mp.sqrt(2) - 1) ** 2 - (1 - 2 * (mp.sqrt(2) - 1))) < mp.mpf(10) ** -40
check("C4 [E] Für eine allgemeine ladungskonjugationssymmetrische Energie auf Z_4 ist die selbstduale Menge eine Linie, "
      "2 w(1) + w(2) = w(0). Das Villain-Gewicht trifft sie bei a = pi/2, die Kosinus-Energie bei w(1)/w(0) = sqrt(2) - 1 "
      "(zwei Ising-Modelle; der Hadamard-Fixpunkt aus z144). Erst quadratisch plus selbstdual legt den Winkel fest", okat)
befund("C Winkel", "nein (nativ), bedingt (diskret unter benannten Prämissen)",
       "Teil C: Kein natives Datum der Kette legt g^2 fest. (a) Der Modulfluss vertauscht immer mit der Uhr; gleich der "
       "Uhr bleibt er nur bei g^2 = 0, und weil sich die Niveaus nicht proportional verschieben, geht die Periode "
       "generisch verloren. Der elektrische Term geht näherungsweise mit dem Boost-Gewicht ein, der Schnittlink mit "
       "wenig Gewicht. (b) P_viertel ist g^2 T = pi/2 mit zwei neuen Größen. (c) Selbstdualität "
       "legt das euklidische Gewicht eindeutig auf pi/2 = beta/8 (tau = i der Uhr-Thetafunktionen), aber nur für die "
       "quadratische Form; das reelle Gatter ist die T-Matrix derselben Funktionen, und 'Schritt = T' ist P_viertel. "
       "Diskret wird der Winkel unter P_Z4quad (Schrittphase hängt nur von E mod 4 ab, Feldgesetz quadratisch): acht "
       "Werte, davon vier wechselwirkend = zwei Zeitrichtungen mal zentrale Parität.")

# ---------------------------------------------------------------------------------------------------------------------
say("\nTeil D  Exakte Diagonalisierung: Bleibt die Naht die Naht?")
Zex = float(mp.mpf(2) ** mp.mpf("0.25") * mp.exp(-mp.pi / 24))
fs = {}
for N, c_ in ((8, 4), (10, 4), (12, 6)):
    KN = K12 if N == 12 else (K10 if N == 10 else K8)
    o0, o1 = KN.niveaus(-2.0, 0.0, -0.25, cut=c_), KN.niveaus(-2.0, 1.0, -0.25, cut=c_)
    fs[N] = [[r(o0[k], 5) for k in ("p", "h", "ph")], [r(o1[k], 5) for k in ("p", "h", "ph")]]
    say(f"   {N} Plätze, Schnitt nach {c_}: g^2 = 0: p, h, ph = {fs[N][0]}; g^2 = 1: {fs[N][1]}")
ZAHLEN["endliche Größe: N -> [p, h, ph bei g^2 = 0], [bei g^2 = 1] (Einheit pi)"] = fs
o0 = K12.niveaus(-2.0, 0.0, -0.25)
check("D1 [N] g = 0, 12 Plätze: Vielteilchenniveaus pi (1, 3, 4, 5, 6) für (Q=+1, Q=-1, Paar, zweites Q=+1, Q=+2) auf 2e-3 "
      "und Zustandssumme 2^(1/4) e^(-pi/24) = 1,043298 auf 1e-5: Die halbe freie Kette ist die Naht (z152)",
      max(abs(o0["p"] - 1), abs(o0["h"] - 3), abs(o0["ph"] - 4), abs(o0["p1b"] - 5), abs(o0["p2"] - 6)) < 2e-3
      and abs(o0["Z"] - Zex) < 1e-5)
check("D1 [N] Endliche Größe: Zwischen 8, 10 und 12 Plätzen ändern sich p, h, ph bei g^2 = 1 um weniger als 3e-3 pi "
      "(Korrelationslänge 0,567 Plätze); die Deformation durch g^2 ist zehn- bis hundertmal größer",
      max(abs(fs[8][1][i] - fs[12][1][i]) for i in range(3)) < 3e-3 and fs[12][1][1] - 3 > 0.15)
G2 = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0)
LES = (("gleich", -0.25, "CT-symmetrisch, theta = 0 (nativ)"), ("gleich", 0.25, "CT-symmetrisch, theta = pi"),
       ("gleich", 0.0, "gleichmäßiges Feld theta = pi/2"), ("gestaffelt", 0.0, "Kogut-Susskind theta = 0 ohne Massenkorrektur"))
scan = {}
for bg, F0, name in LES:
    say(f"   Hintergrund {bg}, F0 = {F0} ({name}):")
    rows = []
    for g2 in G2:
        o = K12.niveaus(-2.0, g2, F0, bg)
        row = [g2, r(o["p"], 4), r(o["h"], 4), r(o["p"] / (o["p"] + o["h"]), 4), r(o["ph"] - o["p"] - o["h"], 4),
               r(o["p1b"], 4), r(o["p2"], 4), r(o["Z"], 5)]
        rows.append(row)
        say(f"      g^2 = {g2}: p = {row[1]}, h = {row[2]}, Versatz p/(p+h) = {row[3]}, Paar - (p + h) = {row[4]}, "
            f"zweites Q=+1: {row[5]}, Q=+2: {row[6]}, Z = {row[7]}")
    scan[f"{bg}, F0 = {F0}"] = rows
ZAHLEN["Raster 12 Plätze: [g^2, p, h, Versatz, Paar - (p+h), zweites Q=+1, Q=+2, Z] (Niveaus in pi)"] = scan
nat = scan["gleich, F0 = -0.25"]
check("D2 [N] Native Lesart (F0 = -1/4): Mit g^2 wachsen p und h monoton, der Versatz p/(p+h) fällt monoton unter 1/4 "
      "(0,2432 bei g^2 = 1), die Zustandssumme fällt. Die Leiter pi (4Z + 1) deformiert sich glatt",
      all(nat[i + 1][1] > nat[i][1] and nat[i + 1][2] > nat[i][2] and nat[i + 1][3] < nat[i][3] and nat[i + 1][7] < nat[i][7]
          for i in range(len(nat) - 1)) and abs(nat[3][3] - 0.2432) < 5e-4)
ng = []
for bg, F0, name in LES:
    oa, ob = K12.niveaus(-2.0, 1e-3, F0, bg), K12.niveaus(-2.0, -1e-3, F0, bg)
    ng.append(((oa["ph"] - oa["p"] - oa["h"]) - (ob["ph"] - ob["p"] - ob["h"])) / 2e-3)
ZAHLEN["d/dg^2 von Paar - (p + h) bei g = 0, vier Lesarten (Einheit pi)"] = [r(x, 4) for x in ng]
say(f"   Steigung von Paar - (p + h) bei g = 0 in den vier Lesarten: {[r(x, 4) for x in ng]} pi je g^2")
check("D2 [N] Nicht quasifrei: Paar - (p + h) ist für jedes g^2 > 0 negativ und fällt monoton, in allen vier Lesarten. Die "
      "Steigung bei g = 0 ist in allen dieselbe (-0,067 pi je g^2): Der Zweikörperanteil hängt in erster Ordnung nicht vom "
      "Hintergrundfeld ab",
      all(all(rows[i + 1][4] < rows[i][4] - 1e-3 for i in range(len(rows) - 1)) for _, rows in sorted(scan.items()))
      and max(ng) - min(ng) < 2e-4 and abs(ng[0] + 0.0672) < 1e-3)
check("D2 [N] 'Modulfluss = Uhr' verlangt p = 1 mod 4. In allen vier Lesarten liegt p für 0 < g^2 <= 8 strikt zwischen "
      "1,005 und 1,7: kein g^2 > 0 im Raster stellt die Gleichung wieder her",
      all(1.005 < row[1] < 1.7 for rows in scan.values() for row in rows[1:]))
# Masse nachstellen: p = pi erzwingen
nach = []
for g2 in (0.25, 0.5, 1.0, 2.0, 4.0, 8.0):
    ms = brentq(lambda m: K12.niveaus(m, g2, -0.25)["p"] - 1.0, -4.0, -0.3, xtol=1e-11)
    o = K12.niveaus(ms, g2, -0.25)
    nach.append([g2, r(ms, 4), r(o["h"], 4), r(o["ph"], 4), r(o["p1b"], 4), r(o["p2"], 4)])
    say(f"   Masse nachgestellt (p = pi): g^2 = {g2}: m = {nach[-1][1]}, h = {nach[-1][2]}, Paar = {nach[-1][3]}, "
        f"zweites Q=+1: {nach[-1][4]}, Q=+2: {nach[-1][5]}  (Naht: 3, 4, 5, 6)")
ZAHLEN["Masse nachgestellt: [g^2, m, h, Paar, zweites Q=+1, Q=+2] (Niveaus in pi)"] = nach
check("D3 [N] Stellt man die Masse nach, so dass das erste Niveau pi bleibt, weichen die anderen ab: h = 3,03 bis 3,62 "
      "statt 3, monoton in g^2. Mit m = -1,872 bei g^2 = 1",
      all(nach[i + 1][2] > nach[i][2] for i in range(len(nach) - 1)) and nach[0][2] > 3.03 and abs(nach[2][1] + 1.8718) < 1e-3)
# drei Stellgrößen (m, F0, Dimerisierung t2) gegen fünf Niveaus bei festem g^2
rest3 = []
for g2 in (0.0, 0.25, 0.5, 1.0, 2.0, 4.0):
    def f5(x, g2=g2):
        o = K12.niveaus(x[0], g2, x[1], t2=x[2])
        return [o["p"] - 1, o["h"] - 3, o["ph"] - 4, o["p1b"] - 5, o["p2"] - 6]
    sol = least_squares(f5, [-2.0, -0.25, 1.0], xtol=1e-12, ftol=1e-12, gtol=1e-12)
    o = K12.niveaus(sol.x[0], g2, sol.x[1], t2=sol.x[2])
    rest3.append([g2, r(np.max(np.abs(sol.fun)), 4), r(sol.x[0], 2), r(sol.x[1], 2) if g2 > 0 else None, r(sol.x[2], 2),
                  r(o["ph"] - o["p"] - o["h"], 3)])
    say(f"   m, F0 und Dimerisierung frei (lokale Suche ab dem nativen Punkt), g^2 = {g2}: kleinster Rest der fünf Niveaus (pi, 3pi, 4pi, 5pi, 6pi) = "
        f"{rest3[-1][1]} pi (m = {rest3[-1][2]}, F0 = {rest3[-1][3]}, t2/t1 = {rest3[-1][4]}); Paar - (p + h) = {rest3[-1][5]}")
ZAHLEN["m, F0, t2 frei: [g^2, kleinster Rest von fünf Niveaus in pi, m, F0, t2/t1, Paar - (p+h)]"] = rest3
check("D4 [N] Auch mit Masse, Hintergrundfeld und Dimerisierung als freien Stellgrößen lassen sich die fünf ersten Niveaus "
      "für g^2 > 0 nicht treffen (kleinste Quadrate, lokale Suche ab dem nativen Punkt): Der kleinste Rest wächst monoton mit g^2 (0,010 pi bei g^2 = 0,25) und ist nur bei "
      "g^2 = 0 auf dem Niveau der endlichen Größe (< 1e-3). Die Nichtadditivität Paar - (p + h) < 0 bleibt, weil sie ein "
      "Zweikörpereffekt ist",
      rest3[0][1] < 1e-3 and all(rest3[i + 1][1] > rest3[i][1] for i in range(len(rest3) - 1)) and rest3[1][1] > 5e-3
      and all(x[5] < -0.01 for x in rest3[1:]))
befund("D Naht", "nein",
       "Teil D: Die Naht bleibt nicht die Naht. Mit g^2 > 0 im Hamiltonoperator, dessen Grundzustand die Naht definiert, "
       "deformiert sich die Leiter glatt und schon in erster Ordnung (Versatz unter 1/4, Niveaus nicht mehr additiv). Kein "
       "g^2 im Raster bis 8 stellt sie wieder her, auch nicht mit nachgestellter Masse; mit freier Masse, freiem "
       "Hintergrundfeld und Dimerisierung ist das bis g^2 = 4 geprüft (lokale Suche, 12 Plätze). In der geprüften Familie "
       "(Dimerisierung, Wechselpotential, Hintergrundfeld, g^2) liegt die Naht bei g^2 = 0; das passt zu z152, wo zwei "
       "verschachtelte Nähte den freien Hamiltonoperator festlegen. Der elektrische Term kann deshalb nur als Schritt der "
       "Komposition gelesen werden (P_schritt), nicht als Teil der Dynamik, die den Nahtzustand festlegt.")

# ---------------------------------------------------------------------------------------------------------------------
say("\nTeil E  Bilanz")
bilanz = [
    ("E", "Form (P_null)", "Gauß-Gesetz: Das Feld auf einem Schnittlink ist F + Nahtladung. Seine Energie ist eine Funktion "
     "der Ladung ohne Oszillatoranteil; zwei Nähte mit je einem Schnittlink geben exakt H_Q. P_null folgt aus P_eich mit "
     "P_schnitt (Feld nur auf den Schnittlinks). Mit Feld auf den inneren Links ist der Rest nicht null."),
    ("E negativ", "Kreuzterm", "Der Term Q_A Q_B des Gauß-Gesetzes kommt nur von Außenlinks, als (Q_A + Q_B)^2, und ist "
     "zentral. In der Nahtbasis ist die geeichte Uhr ein Produkt von Einnahtphasen: triviales Gatter, wie z121 und z145."),
    ("E", "Basis", "Das Gatter ist die elektrische Phase in den zirkularen Moden, gleichwertig die elektrische Phase der "
     "Ladung J_y der Vierteldrehung P_sym. Die Eichung liefert den Basiswechsel B nicht; P_sym bleibt."),
    ("E/N", "Viertel im Feld", "F0 = -+1/4 reproduziert (Q -+ 1/4)^2/2 an jedem Schnitt, bis auf die Konstante 1/32 je Link "
     "(Twistkonstante e^(2 pi i/32) aus z145). Es ist der Versatz halbe Füllung gegen Dirac-See, in Standardlesart "
     "theta = 0 mit Massenverschiebung g^2 a/8, kein Feld theta = pi/2. Fest ist es unter P_CT bis auf ein Bit."),
    ("E", "Z_4 statt U(1)", "In einer Z_4-Eichtheorie ist jede elektrische Phase von selbst eine Funktion der Uhrwerte. "
     "Frei sind dann drei Phasen statt einer. Mit quadratischem Feldgesetz ist die Schrittphase auf (pi/2) Z "
     "quantisiert: die acht Uhrgatter, vier davon wechselwirkend (zwei Zeitrichtungen mal zentrale Parität)."),
    ("E/N negativ", "Winkel nativ", "Modulperiode, Boost-Gewicht und die Zeiten der freien Kette legen g^2 nicht fest. Der "
     "Modulfluss bleibt nur bei g^2 = 0 gleich der Uhr."),
    ("E/C", "Selbstdualität", "Das euklidische Villain-Gewicht auf Z_4 ist genau bei a = pi/2 = beta/8 selbstdual; das "
     "ist tau = i für die Uhr-Thetafunktionen (Gewichte k^2/8). Das Gatter ist die T-Matrix derselben Funktionen. tau = i "
     "legt den euklidischen Schritt fest, nicht den reellen: 'Schritt = ganzer Dehn-Twist' bleibt P_viertel (unter "
     "P_quad). Für allgemeine Z_4-Energien ist die selbstduale Menge eine Linie. Dieselbe Fourier-Bedingung ist der "
     "selbstduale Punkt der Z_4-Eichtheorie in 3+1 Dimensionen [L, nicht gerechnet]."),
    ("N negativ", "Naht unter g^2", "Im Grundzustand der geeichten Kette ist die halbe Kette für g^2 > 0 nicht die Naht "
     "(12 Plätze, g^2 bis 8; mit freier Masse, freiem F und Dimerisierung bis g^2 = 4). P_eich ist nur als "
     "Kompositionsschritt verträglich (P_schritt)."),
    ("C", "Prämissenbilanz", "P_uhr (eine Prämisse) würde ersetzt durch P_eich (= P_E aus z135; P_link ist in 1+1 leer) + "
     "P_schnitt + P_Z4quad (oder g^2 T = pi/2 als Zahl) + P_schritt; dazu P_CT für F. P_sym, P_ord und P_univ werden "
     "nicht ersetzt. Die Zahl der Prämissen sinkt nicht."),
    ("C", "K15", "In 1+1 Dimensionen hat das Eichfeld keine eigenen Freiheitsgrade. P_eich liefert genau den elektrischen "
     "Term, der in z135 fehlte, und der ist H_Q. Die absolute Kopplung bleibt offen: P_viertel und die offene Kopplung "
     "von K15 sind dieselbe Zahl g^2 T. Photonen und Plaketten gibt es erst ab 2+1."),
    ("O", "Nächster Test", "Die geeichte Schrittkette U = (Transferschicht) mal e^(-i (pi/2) H_el) auf den Schnittlinks: Ist "
     "sie das Destri-de-Vega-Lichtkegelgitter aus z121, und ist die halbe Kette ihres stationären Zustands die Naht? "
     "Danach die Z_4-Eichtheorie in 3+1 am selbstdualen Punkt (dort ist die Kopplung keine freie Zahl)."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
urteil = ("Bedingt, und ohne Gewinn in der Zählung. (1) Die Form folgt: Das Gauß-Gesetz macht das Feld auf dem Schnittlink "
          "zur Nahtladung plus Hintergrund, seine Energie ist eine Funktion der Ladung ohne Oszillatoranteil, und zwei "
          "Schnittlinks geben exakt H_Q = (Q_A^2 + Q_B^2)/2. Das ist P_null, unter P_eich mit der Geometrie P_schnitt. "
          "(2) Die Wechselwirkung folgt nicht aus der Eichung: Der Kreuzterm Q_A Q_B des Gauß-Gesetzes ist zentral, in "
          "der Nahtbasis ist die geeichte Uhr ein Produkt von Einnahtphasen. Das Gatter ist elektrisch in den zirkularen "
          "Moden, also bezüglich der Ladung der Vierteldrehung P_sym; der Basiswechsel B bleibt Prämisse. (3) Das Viertel "
          "im Hintergrund ist auf der Kette nativ bis auf ein Bit: F0 = -+1/4 ist der Versatz zwischen halber Füllung und "
          "Dirac-See, dieselbe Zahl wie der Leiterversatz aus z152, in Standardlesart theta = 0 mit der Massenverschiebung "
          "g^2 a/8. Die Idee 'theta = pi/2' ist damit korrigiert. Die z113-Formel wird bis auf 1/32 je Link reproduziert, "
          "und bei g^2 t = 2 pi ist die elektrische Entwicklung die Uhr mal der Twistkonstante aus z145. (4) Der Winkel "
          "folgt nicht nativ. g^2 T = pi/2 ist eine Bedingung an zwei neue Größen; weder Modulperiode noch Boost-Gewicht "
          "noch eine Zeit der freien Kette legt sie fest. Diskret wird der Winkel erst unter der benannten Prämisse "
          "P_Z4quad (Z_4-Eichgruppe mit quadratischem Feldgesetz: acht Werte, vier wechselwirkend). Die Selbstdualität "
          "tau = i legt nur den euklidischen Schritt fest (a = pi/2 = beta/8); das Gatter ist die T-Matrix derselben "
          "Uhr-Thetafunktionen. Beides sind Umformulierungen von P_viertel, keine Herleitung. (5) Als Dynamik ist die Eichung mit der Naht unverträglich: Im Grundzustand der geeichten Kette "
          "deformiert sich die Leiter pi (4Z + 1) schon in erster Ordnung in g^2, die Niveaus sind nicht mehr additiv, und "
          "kein g^2 im Raster bis 8 stellt sie wieder her (12 Plätze; mit freier Masse, freiem Feld und Dimerisierung bis "
          "g^2 = 4 geprüft, lokale Suche). In der geprüften Familie liegt die Naht bei g^2 = 0; der elektrische Term kann "
          "nur ein Schritt der Komposition sein. (6) Für den Kandidaten: P_uhr bleibt Prämisse, die neun Prämissen sind unverändert. Gewonnen "
          "ist eine physikalische Lesart (P_null ist das Gauß-Gesetz, P_viertel ist die Eichkopplung mal Schrittdauer und "
          "damit dieselbe offene Zahl wie die Kopplung von K15) und ein scharfes Negativ für die geeichte Kette als "
          "Hamiltonoperator der Naht. K13 ist nicht berührt. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z155_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z155.txt").write_text("\n".join(OUT) + "\n")
