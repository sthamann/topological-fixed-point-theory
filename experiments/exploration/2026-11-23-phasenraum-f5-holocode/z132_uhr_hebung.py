"""Z132: Folgeschritt 'Ist P_lift nativ? Die Hebung der Uhr auf die 16, ihr Zweig und die fehlende Drei'.

Leitlinie 'Rekonstruktion statt Auswahl': ein Zweig zählt nur als nativ, wenn vorhandene TFPT-Struktur den anderen
ausschließt. Firewalls: DIMENSION.UPLIFT.FIREWALL.01 (keine 3+1-Aussage; Helizität bleibt Kinematik der Einpaarschicht),
GRAV.NONCIRCULAR.01 (die Uhr tritt nur als Symmetrie bzw. Gewicht auf, nie als Zeit). H6 und die gelösten Bausteine (die 16
mit ihren Hyperladungen, E8, Familien, Massen, alpha) sind feste Eingaben und werden nicht neu geöffnet; nichts wird
befördert. Bestand: z82 (Träger e_m = z^m/(z^4 - 1), R = diag(1, i, -1, -i, 1), det R = -1, R = (-1_5) g mit g in G_SM,
Spinorhebung R~ = det(R)^(-1/2) Lambda(R) mit zwei Zweigen, R~ Lambda(-R)^-1 = z_X^(-+1); naive Hebung eichtrivial; die Regel
R^2 = (-1)^F wählt die Spinorhebung; das Original erreicht über die Hardy-Marken nur die 2; Statistik 8/8;
Determinantenschranke), z92 (die Drei = Sym^2(2)/(uv)^2 mit erzwungenem Teiler; SL(2)-Hebung; Z8 an den Marken), z97
(E8 = D5 + A3 über Z4, z_X z_A3 = 1; e^(i pi K) hat auf den Ladungen Ordnung 8), z113 (vier chirale Nähte: (E8)_1, Uhr auf
R = -i eingefroren), z116 (R = i Gamma_11, Nahtindex A_q = S^+(1 - q S^4) mit Index +1, J_F), z117 (P_glue: I_T = I_F,
Uhrspur c = 64 x Chiralität x Schnitt x Uhrrichtung, 'das Vorzeichen hängt am Hebungszweig'), z118 ('Uhrrichtung durch
die Marken fest, J R J^-1 = R'), z122 (P_IF: native Paarung, 'R = i (die 16) auf die Originalgerade', Blatttausch Theta).

Bezeichnungen: S = Lambda^* C^5 (x) det^(-1/2) (32), 16 = Lambda^gerade, 16q = Lambda^ungerade, nu^c = Lambda^0 (das
G_SM-Singulett der 16), X = 5 - 2p, z_X = exp(i pi X/2) (= i auf der 16). Zweig eps = +-1: R~_eps = eps i Lambda(R), also
R~_eps(nu^c) = eps i und R~_eps = z_X^eps Lambda(-R). eps = +1 ist die Konvention von z82 ('Zweig i^X'), z116-z118
(R = i Gamma_11) und z122 ('R = i, die 16'). sigma_IF: Vorzeichen der nativen Paarung (z122), sigma_glue: Vorzeichen der
Spinormodul-Verklebung I_T = sigma z_X (z117). Statuskürzel: E exakt, C bedingt, L Lesart, O offen.

Teil 1  P_lift als endliche Aussage: Bestandsabgleich; die zwei Hebungen (Ordnung, Spinregel modulo G_SM, Klasse z_X^eps,
        gleiche Wirkung auf den zehn Majoranas); Gruppe Z4 x Z2 mit genau zwei Schnitten; eichinvariante Lesung über nu^c
        und G_SM-Singuletts; der Zweig als mod-8-Datum ganzzahliger Logarithmen, als Hauptwinkel der R = -1-Mode e_2, als
        Determinante des Ordnung-4-Lifts auf die homogenen Koordinaten (u, v) und als Richtung der U(1)_X-Halbdrehung.
Teil 2  Wählt ein natives Datum den Zweig? X-Ladung und 3+2-Schnitte; Hardy-Orientierung und Index von A_q (Filtrationen an
        0 und unendlich, Torusfortsetzung der leichtesten Energien des Originals, Gegenregel); Spiegelungen tau, s, K, J_F;
        verklebter Spinormodul (P_glue, P_IF, Blatttausch Theta); z122s Forderung, Uhrspur und z118; die eingefrorene Uhr
        der Viernaht-(E8)_1 gegen D5 + A3.
Teil 3  Die fehlende Drei: Original und Partner, Kompositionen über bis zu drei Nähte und Quarter mit alpha = 0, 1/2,
        Determinantenschranke, die bosonische Drei aus z92 und ihre Gradierungen (Produkt- gegen Einsrelation, Z8-Daten).
Teil 4  Bilanz und Urteil, Folgen für die Prämissenliste des besten Kandidaten und für das Chiralitätsvorzeichen.
"""
import itertools
import json
from collections import Counter
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.linalg import expm

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
rng = np.random.default_rng(132)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def zs(x):
    """Vierte Einheitswurzel als Text."""
    for t, v in (("1", 1), ("i", 1j), ("-1", -1), ("-i", -1j)):
        if abs(x - v) < 1e-9:
            return t
    return str(np.round(x, 6))


def zt(e):
    """Zweig als Text mit Vorzeichen."""
    return f"{e:+d}" if e is not None else "None"


# ------------------------------------------------------------------ Träger, Spinormodul, Fockoperatoren (wie z82, z116)
SUB = [S for p in range(6) for S in itertools.combinations(range(5), p)]
IDX = {S: j for j, S in enumerate(SUB)}
GER = [IDX[S] for S in SUB if len(S) % 2 == 0]
UNG = [IDX[S] for S in SUB if len(S) % 2 == 1]
PGRAD = np.array([len(S) for S in SUB])
XL = 5 - 2*PGRAD
ZX = np.diag(1j**(XL % 4))                                         # Zentrum von Spin(10)
ZXP = {1: ZX, -1: np.linalg.inv(ZX)}
NU = IDX[()]
I32 = np.eye(32)
EPS = (1, -1)


def Lam(g):
    """Lambda(g) auf Lambda^* C^5: Matrixelemente sind die Minoren det g[T, S] (z82)."""
    M = np.zeros((32, 32), complex)
    for S in SUB:
        for T in SUB:
            if len(S) == len(T):
                M[IDX[T], IDX[S]] = np.linalg.det(g[np.ix_(T, S)]) if S else 1.0
    return M


def erzeuger(j):
    """a_j^+ auf Lambda^* C^5 (e_j ^ e_S mit Vorzeichen)."""
    A = np.zeros((32, 32), complex)
    for S in SUB:
        if j not in S:
            A[IDX[tuple(sorted(S + (j,)))], IDX[S]] = (-1)**sum(1 for k in S if k < j)
    return A


AD = [erzeuger(j) for j in range(5)]
AN = [a.conj().T for a in AD]
assert all(np.allclose(AN[j] @ AD[k] + AD[k] @ AN[j], (j == k)*I32) for j in range(5) for k in range(5))
MAJ = []
for _j in range(5):
    MAJ += [AN[_j] + AD[_j], 1j*(AD[_j] - AN[_j])]


def dLam(X):
    """Lie-Algebra-Wirkung sum_jk X_jk a_j^+ a_k."""
    return sum(X[j, k]*AD[j] @ AN[k] for j in range(5) for k in range(5))


_Xt = rng.normal(size=(5, 5)) + 1j*rng.normal(size=(5, 5))
_Xt = 0.3*(_Xt - _Xt.conj().T)
assert np.allclose(Lam(expm(_Xt)), expm(dLam(_Xt)))                # Lambda und dLambda sind verträglich


def auf16(A):
    return A[np.ix_(GER, GER)]


def vektorwirkung(U):
    """U gamma_b U^-1 = sum_a M_ab gamma_a auf den zehn Majoranas (Bild in SO(10))."""
    Ui = np.linalg.inv(U)
    return np.array([[np.trace(MAJ[a] @ U @ MAJ[b] @ Ui).real/32 for b in range(10)] for a in range(10)])


def hebung(W):
    """Spinorhebung von R = exp(i pi W/2) längs des Weges exp(i t pi W/2): auf Lambda^S der Faktor
    exp(i pi/2 (sum_S w - sum w/2)) (diagonal)."""
    return np.array([np.exp(1j*np.pi/2*(sum(W[i] for i in S) - sum(W)/2)) for S in SUB])


R = np.diag([1j**m for m in range(5)])
LR, LmR = Lam(R), Lam(-R)
RT = {e: e*1j*LR for e in EPS}                                     # R~_eps
DRT = {e: np.diag(RT[e]) for e in EPS}


def zweig(h):
    """Welcher Zweig ist die diagonale Hebung h? (+1, -1 oder None)"""
    z = [e for e in EPS if np.allclose(h, DRT[e])]
    return z[0] if len(z) == 1 else None


# =================================================================== Teil 1
say("Teil 1: P_lift als endliche Aussage")

# --- 1.0 Bestandsabgleich (nur lesen)
d82 = json.loads((HERE / "z82_ergebnis.json").read_text())
d113 = json.loads((HERE / "z113_ergebnis.json").read_text())
d116 = json.loads((HERE / "z116_ergebnis.json").read_text())
d117 = json.loads((HERE / "z117_ergebnis.json").read_text())
d118 = json.loads((HERE / "z118_ergebnis.json").read_text())
d122 = json.loads((HERE / "z122_ergebnis.json").read_text())
zw82 = d82["zahlen"]["R_tilde_mal_Eichinverse_auf_16"]
zw82_num = {complex(k): complex(v) for k, v in zw82.items()}
eigene = {}
for s_ in (1j, -1j):                                                # z82: R~ = Lambda(R)/s mit s^2 = det R
    M = auf16((LR/s_) @ np.linalg.inv(LmR))
    eigene[s_] = complex(np.round(M[0, 0], 12))
bestand82 = len(zw82_num) == 2 and all(any(abs(k - s_) < 1e-12 and abs(v - eigene[s_]) < 1e-12
                                           for k, v in zw82_num.items()) for s_ in eigene)
eigene_txt = {zs(k): zs(v) for k, v in eigene.items()}
treffer = lambda d, frag: [v for k, v in d["checks"].items() if frag in k]
konv = {"z116: R = i Gamma_11": treffer(d116, "R = i Gamma_11 hat Ordnung 4"),
        "z117: I_F = R = i Gamma_11": treffer(d117, "I_F = R = i Gamma_11"),
        "z118: Uhrrichtung durch die Marken": treffer(d118, "Uhrrichtung durch die"),
        "z122: R = i (die 16) auf die Originalgerade": treffer(d122, "R = i (die 16, z_X = i)")}
hebungszweig_117 = any("Hebungszweig" in k for k in d117["befunde"])
koeff113 = d113["zahlen"]["chirales_original"]["relative Koeffizienten (Phase/2pi) von R, R^2, R^3 bei N = 4"]
ZAHLEN["Bestand"] = {"z82 R~ Lambda(-R)^-1 je s (Datei)": zw82, "nachgerechnet": eigene_txt,
                     "Konventionsstellen (Prüfungen bestanden)": {k: v for k, v in konv.items()},
                     "z113 Koeffizienten R, R^2, R^3": koeff113}
check("[E] Bestandsabgleich (nur gelesen): z82s Zweigtabelle R~ Lambda(-R)^-1 für R~ = Lambda(R)/s ist nachgerechnet "
      f"({eigene_txt}: {bestand82}); z116 (R = i Gamma_11), z117 (I_F = R = i Gamma_11), z118 (Uhrrichtung durch die Marken) "
      "und z122 (R = i, die 16, auf die Originalgerade) benutzen alle den Zweig eps = +1 (s = -i), z117 nennt den "
      f"Hebungszweig als offenes Vorzeichen ({hebungszweig_117}); z113s Uhrkoeffizienten {koeff113}",
      bestand82 and eigene_txt == {"i": "-i", "-i": "i"} and all(len(v) == 1 and v[0] for v in konv.values())
      and hebungszweig_117 and koeff113 == {"R": "1/4", "R2": "1/2", "R3": "3/4"})

# --- 1.1 die beiden Hebungen
ok_ord = all(np.allclose(np.linalg.matrix_power(RT[e], 4), I32) and not np.allclose(RT[e] @ RT[e], I32)
             and not np.allclose(RT[e] @ RT[e], -I32) for e in EPS)
ok_quad = all(np.allclose(RT[e] @ RT[e], -Lam(R @ R)) for e in EPS)
R2_eich = np.allclose(R @ R, expm(6j*np.pi*np.diag([-1/3, 1/2, -1/3, 1/2, -1/3])))   # R^2 = exp(6 pi i Y) in G_SM
ok_zx = all(np.allclose(RT[e] @ np.linalg.inv(LmR), ZXP[e]) for e in EPS)
ok_kern = np.allclose(RT[-1], ZX @ ZX @ RT[1]) and np.allclose(ZX @ ZX, -I32)
MV = {"eps = +1": vektorwirkung(RT[1]), "eps = -1": vektorwirkung(RT[-1]), "naiv Lambda(R)": vektorwirkung(LR)}
mv0 = MV["eps = +1"]
ok_vek = (all(np.allclose(M, mv0) for M in MV.values()) and np.allclose(mv0 @ mv0.T, np.eye(10))
          and abs(np.linalg.det(mv0) - 1) < 1e-12 and not np.allclose(mv0, np.eye(10)))
winkel = [float(np.round(np.angle(mv0[2*j, 2*j] + 1j*mv0[2*j + 1, 2*j])/(np.pi/2), 9)) % 4 for j in range(5)]
naiv = np.allclose(auf16(LR), auf16(LmR)) and np.allclose(LR @ LR, Lam(R @ R))
diff_nu = RT[1][NU, NU] - RT[-1][NU, NU]
ZAHLEN["Hebungen"] = {"R~_eps(nu^c)": {e: zs(RT[e][NU, NU]) for e in EPS},
                      "Drehwinkel der Vektorwirkung je Ebene (Einheiten pi/2)": winkel}
check("[E] Die zwei Hebungen R~_eps = eps i Lambda(R) auf S (32): beide haben genau Ordnung 4, R~^2 = -Lambda(R^2) mit "
      f"Lambda(R^2) = Lambda(e^(6 pi i Y)) eichartig ({R2_eich}), also R~^2 = (-1)^F modulo G_SM; R~_eps Lambda(-R)^-1 = "
      f"z_X^eps auf allen 32 ({ok_zx}); R~_-eps = z_X^2 R~_eps = -R~_eps ({ok_kern}): der Zweigwechsel ist das zentrale -1 von "
      f"Spin(10), der Kern von Spin(10) -> SO(10). Beide Zweige und die naive Hebung wirken auf den zehn Majoranas gleich "
      f"({ok_vek}; Drehwinkel je Ebene {winkel} x pi/2 = die Viertelcharaktere von R): jedes Datum, das über die "
      f"Vektordarstellung läuft, ist zweigblind. Naive Hebung eichtrivial auf der 16 und ohne Spinregel ({naiv}). "
      f"Kontrolle: R~_+(nu^c) - R~_-(nu^c) = {zs(diff_nu/2)} x 2",
      ok_ord and ok_quad and R2_eich and ok_zx and ok_kern and ok_vek and naiv and abs(diff_nu - 2j) < 1e-12
      and winkel == [0.0, 1.0, 2.0, 3.0, 0.0])

# --- 1.2 Gruppe: Urbild von <R> in {(g, s): s^2 = det g}
ELEM = [(k, s) for k in range(4) for s in ((1, -1) if k % 2 == 0 else (1j, -1j))]
mul = lambda x, y: ((x[0] + y[0]) % 4, x[1]*y[1])
gleich = lambda x, y: x[0] == y[0] and abs(x[1] - y[1]) < 1e-12


def ordnung(x):
    n, y = 1, x
    while not gleich(y, (0, 1)):
        y, n = mul(y, x), n + 1
    return n


abgeschl = all(any(gleich(mul(x, y), z) for z in ELEM) for x in ELEM for y in ELEM)
abelsch = all(gleich(mul(x, y), mul(y, x)) for x in ELEM for y in ELEM)
ordn = Counter(ordnung(x) for x in ELEM)
schnitte = []
for s_ in itertools.product((1, -1), (1j, -1j), (1, -1), (1j, -1j)):
    sec = {k: (k, s_[k]) for k in range(4)}
    if all(gleich(mul(sec[k], sec[l]), sec[(k + l) % 4]) for k in range(4) for l in range(4)):
        schnitte.append(tuple(zs(s_[k]) for k in range(4)))
schnitte.sort()
ZAHLEN["Gruppe"] = {"Ordnungen im Urbild": dict(sorted(ordn.items())), "Schnitte k -> (R^k, s_k)": schnitte}
check("[E] Das Urbild von <R> in der Doppelüberlagerung {(g, s): s^2 = det g} hat 8 Elemente, ist abgeschlossen und abelsch "
      f"({abgeschl}, {abelsch}) mit Ordnungen {dict(sorted(ordn.items()))}, also Z4 x Z2 (nicht Z8, nicht Q8); genau "
      f"{len(schnitte)} Schnitte sind Homomorphismen ({schnitte}), s_1 = -eps i gehört zu R~_eps: Die Gruppentheorie bietet "
      "zwei gleichwertige Z4-Hebungen und bevorzugt keine",
      len(ELEM) == 8 and abgeschl and abelsch and ordn == Counter({1: 1, 2: 3, 4: 4}) and len(schnitte) == 2
      and schnitte == [("1", "-i", "-1", "i"), ("1", "i", "-1", "-i")])

# --- 1.3 eichinvariante Lesung
FARBE, SCHWACH = [0, 2, 4], [1, 3]
E5 = lambda a, b: np.outer(np.eye(5)[a], np.eye(5)[b])
GENS = []
for blk in (FARBE, SCHWACH):
    for a, b in itertools.combinations(blk, 2):
        GENS += [E5(a, b) + E5(b, a), -1j*E5(a, b) + 1j*E5(b, a)]
GENS += [E5(0, 0) - E5(2, 2), E5(0, 0) + E5(2, 2) - 2*E5(4, 4), E5(1, 1) - E5(3, 3)]
YHY = np.diag([-1/3, 1/2, -1/3, 1/2, -1/3])
GENS.append(YHY)
G16 = [auf16(dLam(X)) for X in GENS]


def nullraum(M, tol=1e-9):
    _, s, vh = np.linalg.svd(M)
    r = int(np.sum(s > tol*max(1.0, s[0])))
    return vh[r:].conj()


N16 = nullraum(np.vstack(G16))
nu16 = np.zeros(16)
nu16[GER.index(NU)] = 1
sing16 = N16.shape[0] == 1 and abs(abs(N16[0] @ nu16) - 1) < 1e-12
I16 = np.eye(16)
G2 = [np.kron(A, I16) + np.kron(I16, A) for A in G16]
N2 = nullraum(np.vstack(G2))
nunu = N2.shape[0] == 1 and abs(abs(N2[0] @ np.kron(nu16, nu16)) - 1) < 1e-12      # das Singulett ist nu^c x nu^c
RT16 = {e: auf16(RT[e]) for e in EPS}
paar = {e: np.allclose((np.kron(RT16[e], RT16[e]) @ N2.T), -N2.T) for e in EPS}
eich2 = np.allclose(np.kron(auf16(LmR), auf16(LmR)) @ N2.T, N2.T)
zufall_ok = True
for _ in range(6):
    c_ = rng.normal(size=len(GENS))
    g_ = expm(1j*sum(c*X for c, X in zip(c_, GENS)))
    zufall_ok &= abs(np.linalg.det(g_) - 1) < 1e-10 and np.allclose(Lam(g_)[:, NU], I32[:, NU])
drei = {e: zs(RT[e][NU, NU]**3) for e in EPS}
ZAHLEN["Eichinvariant"] = {"G_SM-Singuletts in 16": N16.shape[0], "in 16 x 16": N2.shape[0],
                           "R~ x R~ auf 16 x 16-Singuletts": {e: "-1" if paar[e] else "?" for e in EPS},
                           "R~^3 auf nu^c x nu^c x nu^c": drei}
check("[E] Eichinvariante Lesung: Unter su(3) + su(2) + u(1)_Y (Träger, R^2-Schnitt) hat die 16 genau ein Singulett, nu^c "
      f"({sing16}); jedes g in G_SM (zufällige Gruppenelemente) lässt nu^c fest ({zufall_ok}). Also ist eps i = R~_eps(nu^c) "
      f"die Uhrladung eines G_SM-Singuletts und eichunabhängig. 16 x 16 hat genau ein Singulett, nu^c x nu^c ({nunu}); dort "
      f"wirkt R~ x R~ für beide Zweige als -1 ({paar}, Eichteil trivial: {eich2}); erst ungerade Produkte sehen den Zweig "
      f"(nu^c x nu^c x nu^c: {drei})",
      sing16 and zufall_ok and nunu and all(paar.values()) and eich2 and drei == {1: "-i", -1: "i"})

# --- 1.4 der Zweig als mod-8-Datum, Hauptwinkel
LOGS = list(itertools.product(*[[m + 4*k for k in range(-2, 3)] for m in range(5)]))
H_ALL = np.array([hebung(W) for W in LOGS])
ok_log, zaehl = True, Counter()
for W, h in zip(LOGS, H_ALL):
    e = zweig(h)
    ok_log &= e is not None and sum(W) % 4 == 2 and e == (-1 if sum(W) % 8 == 2 else 1)
    zaehl[e] += 1
HAUPT = [W for W in LOGS if all(abs(w) <= 2 for w in W)]
haupt_zweig = {str(W): zweig(hebung(W)) for W in HAUPT}
ok_haupt = len(HAUPT) == 2 and all(zweig(hebung(W)) == -int(np.sign(W[2])) for W in HAUPT) \
    and all(W[0] == 0 and W[1] == 1 and W[3] == -1 and W[4] == 0 for W in HAUPT)
ZAHLEN["mod 8"] = {"Logarithmen (w_m in m + 4{-2..2})": len(LOGS), "Zweige": {e: zaehl[e] for e in EPS},
                   "Hauptlogarithmen -> Zweig": haupt_zweig}
check(f"[E] Der Zweig als mod-8-Datum: Für alle {len(LOGS)} ganzzahligen Logarithmen R = exp(i pi W/2), w_m in m + 4{{-2..2}}, "
      "ist die Spinorhebung längs exp(i t pi W/2) einer der beiden Zweige; tr W = 2 mod 4 immer (die Charaktere legen tr W nur "
      f"mod 4 fest), und eps = -1 genau für tr W = 2 mod 8, eps = +1 für tr W = 6 mod 8 ({ok_log}; Anteile "
      f"{ {e: zaehl[e] for e in EPS} }). Mit Hauptwinkeln |theta| <= pi sind e_0, e_1, e_3, e_4 eindeutig (0, pi/2, -pi/2, 0), "
      f"nur die R = -1-Mode e_2 hat zwei (+-pi): eps = -sgn theta_2 ({haupt_zweig}). Der Zweig sitzt also auf der "
      "Halbdrehung von e_2, der Mitte der Drei und dem erzwungenen Teiler (uv)^2 von z92",
      ok_log and ok_haupt and zaehl[1] > 0 and zaehl[-1] > 0)

# --- 1.5 GL(2)-Lift auf (u, v), SL(2)-Lift, U(1)_X-Halbdrehung
om = np.exp(1j*np.pi/4)
torus = {}
for mu in (1, 1j, -1, -1j):
    la = 1j*mu
    wirk = np.array([la**m*mu**(4 - m) for m in range(5)])
    torus[zs(mu)] = {"wirkt als R auf Sym^4": bool(np.allclose(wirk, np.diag(R))), "det": zs(la*mu),
                     "Zweig": zweig(np.array([(la*mu)**(-5)*np.prod([la**m*mu**(4 - m) for m in S]) for S in SUB]))}
ok_torus = all(v["wirkt als R auf Sym^4"] for v in torus.values()) and \
    all(v["Zweig"] == (-1 if v["det"] == "i" else 1) for v in torus.values())
sl2 = np.array([om**m*om**(-(4 - m)) for m in range(5)])
ok_sl2 = np.allclose(sl2, -np.diag(R)) and np.allclose(auf16(LmR), auf16(LR)) and np.allclose(LmR @ LmR, Lam(R @ R))
ux = lambda th: np.diag(np.exp(1j*th*(PGRAD - 2.5)))                # Hebung von e^(i th) 1_5 längs des Weges
ok_ux = np.allclose(ux(-np.pi), ZX) and np.allclose(ux(np.pi), np.linalg.inv(ZX)) and \
    all(np.allclose(RT[e], ux(-e*np.pi) @ LmR) for e in EPS)
ZAHLEN["GL2"] = {"Torus diag(i mu, mu)": torus, "SL(2) diag(om, om^-1) auf Sym^4": [zs(x) for x in sl2]}
check("[E] Lift auf die homogenen Koordinaten: Die vier Ordnung-4-Elemente g = diag(i mu, mu) (mu^4 = 1) wirken auf "
      f"Sym^4 = H^0(O(4)) alle als R; ihre kanonische Torushebung (Charakter (lambda mu)^-5) ist eps = -1 genau für det g = i "
      f"und eps = +1 für det g = -i ({torus}). Der SL(2)-Lift diag(om, om^-1) (Ordnung 8, z92) wirkt als -R mit kanonischer "
      f"Hebung Lambda(-R): eichtrivial, ohne Spin-Z4 ({ok_sl2}). U(1)_X: R = e^(i theta_X) 1_5 (-R) mit theta_X = +-pi; die "
      f"Hebung längs theta_X ist z_X für -pi und z_X^-1 für +pi, also eps = -theta_X/pi ({ok_ux}); z82s 'Zweig i^X' ist die "
      "Halbdrehung im Uhrzeigersinn der komplexen Struktur des Trägers",
      ok_torus and ok_sl2 and ok_ux)

befund("P_lift als endliche Aussage: 'Die Uhr des Originals wirkt auf S = Lambda^* C^5 (x) det^(-1/2) als R~_eps = eps i "
       "Lambda(R) mit eps = +1.' Gleichwertig: (a) die Uhrladung des G_SM-Singuletts nu^c ist eps i; (b) tr W = 6 mod 8 für "
       "jeden ganzzahligen Logarithmus von R; (c) die R = -1-Mode e_2 dreht im Hauptwinkel um -pi; (d) die Vierteldrehung "
       "wirkt auf (u, v) mit det = -eps i; (e) ihr U(1)_X-Anteil ist die Halbdrehung -pi, R~ = z_X Lambda(-R); (f) die "
       "Identifikation der Original-Z4 mit dem Zentrum Z4 von Spin(10) ist R -> z_X^eps. Zweigfrei und exakt ist der "
       "Rest: beide Zweige geben modulo G_SM die anomalierelevante Spin-Z4 (z82), die Regel R^2 = (-1)^F wählt die "
       "Spinorhebung; auf der Drei ist die Wirkung die des Trägers, nicht die des Originals (Teil 3)", True)

# =================================================================== Teil 2
say("\nTeil 2: Wählt ein natives Datum den Zweig?")

# --- 2.1 X-Ladung, 3+2-Schnitte, Hyperladung
spek_X = dict(sorted(Counter(XL[GER].tolist()).items()))
Y16 = np.real(np.diag(auf16(dLam(YHY))))
spek_Y = dict(sorted(Counter(np.round(Y16, 9).tolist()).items()))
klasse = {e: np.diag(auf16(RT[e] @ np.linalg.inv(LmR))) for e in EPS}
nu_bilanz = {e: int(round(sum(1 if abs(x - 1j) < 1e-9 else -1 for x in klasse[e]).real)) for e in EPS}
quad_gleich = np.allclose(RT[1] @ RT[1], RT[-1] @ RT[-1])
PW = np.zeros((10, 10))
for j in SCHWACH:
    PW[2*j, 2*j] = PW[2*j + 1, 2*j + 1] = 1
Wv = (np.array([np.eye(5)[0] - np.eye(5)[4], np.eye(5)[2]])/np.array([[np.sqrt(2)], [1]])).T
schnitt_V4 = np.zeros((10, 10))
for col in Wv.T:
    re, im = np.zeros(10), np.zeros(10)
    re[0::2], im[1::2] = col.real, col.real
    schnitt_V4 += np.outer(re, re) + np.outer(im, im)
haelt = all(np.allclose(MV[k] @ P_ @ MV[k].T, P_) for k in MV for P_ in (PW, schnitt_V4))
u1_regel = {"X = 5 - 2p": zweig(np.diag(ZX)*np.diag(LmR)),
            "X' = 2p - 5": zweig(np.diag(np.linalg.inv(ZX))*np.diag(LmR))}
u1_x, u1_xs = u1_regel["X = 5 - 2p"], u1_regel["X' = 2p - 5"]
spek_Xs = dict(sorted(Counter((-XL[GER]).tolist()).items()))
ZAHLEN["X und 3+2"] = {"X auf der 16": spek_X, "Y auf der 16": spek_Y, "mod-16-Bilanz je Zweig": nu_bilanz,
                       "U1-Regel R~ = exp(i pi X/2) Lambda(-R) je Orientierung": u1_regel}
check("[E, negativ] X-Ladung und 3+2-Schnitte: X = 5 - 2p ist ohne R definiert (auf der 16: "
      f"{spek_X}), Y = (I - 5R^2)/12 hängt nur an R^2, und R~_+^2 = R~_-^2 ({quad_gleich}); beide Zweige halten den R^2- und den "
      f"V4-Schnitt (gleiche Vektorwirkung, {haelt}); die mod-16-Bilanz ist für beide 16 eps = 0 mod 16 (nu = {nu_bilanz}, drei "
      "Familien 48 eps, ohne nu^c 15 eps: für beide anomal). U1s Regel 'Zweig aus der X-Ladung', R~ = exp(i pi X/2) "
      f"Lambda(-R), gibt eps = {zt(u1_x)}; dieselbe Regel mit der anderen Orientierung X' = -X (Spektrum "
      f"{spek_Xs}, gleiche Beträge, gleiche Anomaliezählung) gibt eps = {zt(u1_xs)}: Die "
      "Regel ist die Wahl der Orientierung von U(1)_X relativ zur Uhr, kein natives Datum. Kontrolle: nu^c sieht den Zweig "
      "(1.1)",
      spek_X == {-3: 5, 1: 10, 5: 1} and quad_gleich and haelt and nu_bilanz == {1: 16, -1: -16}
      and all(v % 16 == 0 for v in nu_bilanz.values()) and u1_regel == {"X = 5 - 2p": 1, "X' = 2p - 5": -1})

# --- 2.2 Hardy-Orientierung, Index von A_q, Filtrationen, Torusfortsetzung
z, w = sp.symbols("z w")
E_MOD = [z**m/(z**4 - 1) for m in range(5)]
KO = 9
C0 = sp.Matrix([[sp.series(E_MOD[m], z, 0, KO).removeO().coeff(z, k) for k in range(KO)] for m in range(5)])
CI = sp.Matrix([[sp.series(E_MOD[m].subs(z, 1/w), w, 0, KO).removeO().coeff(w, k) for k in range(KO)] for m in range(5)])


def ordnungen(C):
    dims = [5 - C[:, :k].rank() for k in range(KO)]                # dim {f : ord f >= k}
    return sorted(k for k in range(KO - 1) for _ in range(dims[k] - dims[k + 1]))


def ord_je_mode(C):
    return [next(k for k in range(KO) if C[m, k] != 0) for m in range(5)]


ord0, ordI = ordnungen(C0), ordnungen(CI)
om0, omI = ord_je_mode(C0), ord_je_mode(CI)
char_ok = all(sp.I**om0[m] == sp.I**m and sp.I**(-omI[m]) == sp.I**m for m in range(5))
W_null, W_unend = tuple(om0), tuple(-k for k in omI)
zw_null, zw_unend = zweig(hebung(W_null)), zweig(hebung(W_unend))
qn = np.exp(-2*np.pi)
index = {}
for N in (13, 25):
    Sd = np.zeros((N - 1, N))
    for n in range(N - 1):
        Sd[n, n + 1] = 1                                            # S^+ abgeschnitten: C^N -> C^(N-1)
    S3 = np.zeros((N - 1, N))
    for n in range(N):
        if n + 3 < N - 1:
            S3[n + 3, n] = 1
    A = Sd - qn*S3
    for nm, Op in (("H^2 an 0 (A_q)", A), ("getauscht (A_q^T)", A.T)):
        sv = np.linalg.svd(Op, compute_uv=False)
        rk = int(np.sum(sv > 1e-12))
        index.setdefault(nm, []).append((Op.shape[1] - rk, Op.shape[0] - rk))     # (Kern, Kokern)
    _, _, vh_ = np.linalg.svd(A)
    nullv = vh_[-1]*np.sign(vh_[-1][0])
    traeger_null = sorted({int(n) % 4 for n in np.nonzero(np.abs(nullv) > 1e-12)[0]})
LEICHT = {}
for al, q in ((Fr(1, 4), "a"), (Fr(3, 4), "b")):
    for n in range(-3, 4):
        x = al + n
        en, ch = (x, sp.I**int(4*x)) if x > 0 else (-x, sp.I**int(-4*x % 4))   # Teilchen bzw. Loch
        key = str(sp.nsimplify(ch))
        LEICHT[key] = min(LEICHT.get(key, en), en)
w1, w3 = int(4*LEICHT["I"]), int(4*LEICHT["-I"])
al_, be_ = sp.symbols("alpha beta")
los = sp.solve([al_ + 3*be_ - w1, 3*al_ + be_ - w3], [al_, be_], dict=True)
W_torus = tuple(int(m*los[0][al_] + (4 - m)*los[0][be_]) for m in range(5)) if len(los) == 1 else None
W_min = tuple(m % 4 if m < 4 else 0 for m in range(5))
LEICHT_TXT = {{"I": "i", "-I": "-i"}[k]: str(v) for k, v in sorted(LEICHT.items())}
ZAHLEN["Hardy"] = {"Index (H^2 an 0, getauscht)": index, "Verschwindungsordnungen an 0 / unendlich": [om0, omI],
                   "Gewichte an 0 / unendlich": [W_null, W_unend], "Zweige an 0 / unendlich": [zw_null, zw_unend],
                   "leichteste Energien je Uhrklasse des Originals": LEICHT_TXT,
                   "Torusfortsetzung (alpha, beta)": {str(k): str(v) for k, v in los[0].items()} if los else None,
                   "Gewichte Torus / minimal je Mode": [W_torus, W_min],
                   "Zweig Torus / minimal": [zweig(hebung(W_torus)), zweig(hebung(W_min))]}
check("[C] Hardy-Orientierung: A_q = S^+(1 - q S^4) hat auf H^2 am Punkt 0 (Kern, Kokern) = (1, 0), Index +1, mit der "
      f"Nullmode auf den Marken = {traeger_null} mod 4; ausgetauscht (0, 1), Index -1 ({index}); die "
      f"Verschwindungsordnungen von H^0(O(D)) sind an 0 und an unendlich je {ord0} und {ordI}, mit "
      f"Charakteren i^k bzw. i^-k ({char_ok}). Die Filtration am Punkt 0, wo die Uhr um +pi/2 dreht, gibt die Gewichte "
      f"{W_null} und eps = {zt(zw_null)}; die an unendlich {W_unend} und eps = {zt(zw_unend)}. Das Original liefert die "
      f"leichtesten Energien je Uhrklasse ({LEICHT_TXT}: R = i bei 1/4, R = -i bei 3/4, Hardy-Marken "
      f"{w1}, {w3}); ihre Fortsetzung auf einen Torusfluss der homogenen Koordinaten ist eindeutig ({los}) und gibt "
      f"{W_torus}, also eps = {zt(zweig(hebung(W_torus)))}, entgegen der Bestandskonvention. Gegenregel: leichteste Energie "
      f"je Mode, {W_min}, gibt eps = {zt(zweig(hebung(W_min)))}. Die Hardy-Orientierung wählt den Zweig nur zusammen mit einer "
      "Gradierungsregel; beide Regeln sind 'leichteste positive Energie' und unterscheiden sich allein in der zweiten "
      "R = 1-Mode e_4 der Drei (Energie 1 oder 0)",
      index == {"H^2 an 0 (A_q)": [(1, 0), (1, 0)], "getauscht (A_q^T)": [(0, 1), (0, 1)]} and traeger_null == [0]
      and ord0 == ordI == [0, 1, 2, 3, 4] and char_ok
      and zw_null == -1 and zw_unend == 1 and LEICHT == {"I": Fr(1, 4), "-I": Fr(3, 4)} and W_torus == W_null
      and zweig(hebung(W_torus)) == -1 and W_min == (0, 1, 2, 3, 0) and zweig(hebung(W_min)) == 1)

# --- 2.3 Spiegelungen: tau (z -> 1/z-quer), s (z -> 1/z), K, J_F
s_ok = all(sp.simplify(E_MOD[m].subs(z, 1/z) + E_MOD[4 - m]) == 0 for m in range(5))
P_s = np.zeros((5, 5))
for m in range(5):
    P_s[4 - m, m] = -1                                                # s^* e_m = -e_(4-m)
LPs = Lam(P_s)
marken = [1, 1j, -1, -1j]
tau_marken = all(abs(1/np.conj(p) - p) < 1e-12 for p in marken)
s_RR = np.allclose(P_s @ R @ np.linalg.inv(P_s), np.linalg.inv(R))
tau_R = np.allclose(P_s @ np.conj(R) @ np.linalg.inv(P_s), R)
tau_Y = np.allclose(P_s @ YHY @ np.linalg.inv(P_s), YHY)
tau_X = np.allclose(LPs @ np.diag(XL) @ np.linalg.inv(LPs), np.diag(XL))
tau_16 = np.allclose(LPs[np.ix_(UNG, GER)], 0)
JFl = np.eye(32, dtype=complex)
for j in range(5):
    JFl = JFl @ MAJ[2*j]
jf_bild = {e: bool(np.allclose(JFl @ np.conj(RT[e]) @ np.linalg.inv(JFl), RT[e])) for e in EPS}
jf_zx = {e: bool(np.allclose(JFl @ np.conj(ZXP[e]) @ np.linalg.inv(JFl), ZXP[e])) for e in EPS}
jf_tausch = np.allclose(JFl[np.ix_(GER, GER)], 0)
k_bild = {e: bool(np.allclose(np.conj(RT[e]), np.linalg.inv(RT[e]))) for e in EPS}
Rinv = np.linalg.inv(R)


def rel(Hdiag, C):
    """Relativer Zweig einer diagonalen Hebung H der Uhr C: H = rel i Lambda(C) (rel = +-1, sonst None)."""
    q = Hdiag/(1j*np.diag(Lam(C)))
    return next((e for e in EPS if np.allclose(q, e)), None)


def drehwinkel(lam_, punkt):
    """Drehwinkel der Uhr z -> lam z am Fixpunkt 0 bzw. unendlich (w = 1/z: w -> w/lam), in Einheiten pi/2."""
    return int(round(np.angle(lam_ if punkt == "0" else 1/lam_)/(np.pi/2)))


# Operation: (Bild der Hebung als Diagonale, Bilduhr, ihr Drehfaktor, Bild des Hardy-Punkts 0)
OPS = {"J_F (Teilchen-Loch, z116)": (lambda e: np.diag(JFl @ np.conj(RT[e]) @ np.linalg.inv(JFl)), R, 1j, "0"),
       "K (z -> z-quer, R -> R^-1)": (lambda e: np.diag(np.conj(RT[e])), Rinv, -1j, "0"),
       "s (z -> 1/z)": (lambda e: np.diag(LPs @ RT[e] @ np.linalg.inv(LPs)), Rinv, -1j, "unendlich"),
       "tau (z -> 1/z-quer)": (lambda e: np.diag(LPs @ np.conj(RT[e]) @ np.linalg.inv(LPs)), R, 1j, "unendlich")}
op_tab, ok_ops = {}, True
for nm, (bild, uhr, lam_, punkt) in OPS.items():
    relz = {e: rel(bild(e), uhr) for e in EPS}
    korr = drehwinkel(lam_, punkt) == 1                             # Uhr dreht am (Bild des) Hardy-Punkt(s) um +pi/2
    erh = relz == {1: 1, -1: -1}
    op_tab[nm] = {"Zweig relativ zur Bilduhr": relz, "Hardy-Korrelation erhalten": korr}
    ok_ops &= (erh == korr) and (relz in ({1: 1, -1: -1}, {1: -1, -1: 1}))
tau_kippt = op_tab["tau (z -> 1/z-quer)"]["Zweig relativ zur Bilduhr"] == {1: -1, -1: 1}
ZAHLEN["Spiegelungen"] = {"Operationen": op_tab, "J_F R~_eps J_F^-1 = R~_eps": jf_bild,
                          "K: conj R~_eps = R~_eps^-1": k_bild}
check("[E] Spiegelungen. tau: z -> 1/z-quer (Spiegelung am Kreis der Marken) hält jede Marke fest "
      f"({tau_marken}), vertauscht mit der Vierteldrehung (tau R tau^-1 = R: {tau_R}), hält Y, X und die 16 "
      f"({tau_Y}, {tau_X}, {tau_16}) und kippt den Zweig: tau R~_eps tau^-1 = R~_(-eps) ({tau_kippt}). Tabelle (Zweig relativ "
      f"zur Bilduhr, und ob die Uhr am Bild des Hardy-Punkts 0 weiter um +pi/2 dreht): {op_tab}. Genau die Operationen, die "
      f"diese native Korrelation von Uhrrichtung und Hardy-Orientierung erhalten (J_F, s: z -> 1/z mit s^* e_m = -e_(4-m): "
      f"{s_ok}, R -> R^-1: {s_RR}), erhalten eps; die sie brechen (K, tau), kippen eps ({ok_ops}). Der Zweig hängt also an "
      "dieser Korrelation, aber mit einem Vorzeichen, das keine Operation festlegt. J_F hält beide Zweige "
      f"({jf_bild}, z_X^eps: {jf_zx}, tauscht 16 und 16q: {jf_tausch}); die Tomita-artige Konjugation bildet R~_eps auf "
      f"R~_eps^-1 ab ({k_bild}), für beide Zweige kovariant; auch z118s Begründung 'J R J^-1 = R' gilt für beide",
      tau_marken and tau_R and tau_Y and tau_X and tau_16 and tau_kippt and s_ok and s_RR and ok_ops
      and all(jf_bild.values()) and all(jf_zx.values()) and jf_tausch and all(k_bild.values()))

# --- 2.4 verklebter Spinormodul: P_glue, P_IF, Theta, Zweigwechsel B
I8 = np.eye(8)
IT = np.diag([1j]*4 + [-1j]*4)                                      # I_T = Z_SU^-1 auf W = T + T* (z117)
TH = np.block([[np.zeros((4, 4)), np.eye(4)], [np.eye(4), np.zeros((4, 4))]])   # Theta: T <-> T*


def kern_proj(sig, e):
    """Projektor auf ker(I_T (x) 1 - sig 1 (x) z_X^e) in W (x) S (diagonal)."""
    d = np.kron(np.diag(IT), np.ones(32)) - sig*np.kron(np.ones(8), np.diag(ZXP[e]))
    return np.diag((np.abs(d) < 1e-12).astype(float))


P_T16 = np.kron(np.diag([1.0]*4 + [0.0]*4), np.diag([1.0 if len(S) % 2 == 0 else 0.0 for S in SUB]))
mod_tab, ok_mod = {}, True
for sig in (1, -1):
    for e in EPS:
        Pk = kern_proj(sig, e)
        glue = kern_proj(sig*e, 1)
        t16 = int(round(np.trace(Pk @ P_T16)))
        th = np.kron(TH, I32) @ Pk @ np.kron(TH, I32).T
        mod_tab[f"sigma_IF = {sig:+d}, eps = {e:+d}"] = {"= P_glue mit sigma_glue": sig*e, "dim (T x 16) darin": t16,
                                                         "Dimension": int(round(np.trace(Pk)))}
        ok_mod &= (np.allclose(Pk, glue) and int(round(np.trace(Pk))) == 128 and t16 == (64 if sig*e == 1 else 0)
                   and np.allclose(th, kern_proj(-sig, e)) and np.allclose(kern_proj(sig, -e), kern_proj(-sig, e))
                   and np.allclose(np.kron(TH, I32) @ kern_proj(sig, -e) @ np.kron(TH, I32).T, Pk))
ZAHLEN["verklebter Modul"] = mod_tab
check("[E] Verklebter Spinormodul (W (x) S, 256): Die native Paarung in Spinorform, R_Original = sigma_IF Z_SU^-1 mit "
      "R_Original = z_X^eps (P_lift), ist exakt die Verklebung ker(I_T - sigma_glue z_X) von z117 mit sigma_glue = sigma_IF "
      f"eps, also P_glue = P_IF x P_lift ({ok_mod}); die 16 sitzt bei T genau für sigma_IF eps = +1 ({mod_tab}). Der "
      "Blatttausch Theta (T <-> T*) kippt sigma bei festem eps, der Zweigwechsel B kippt eps bei festem sigma, auf dem Modul "
      "wirken beide gleich und B Theta hält alle vier Moduln fest: Die Buchung 'P_IF = P_glue, ein Bit' (z117, z122, "
      "OriginLab) enthält stillschweigend eps = +1",
      ok_mod)

# --- 2.5 z122s Forderung, Uhrspur, z118
forderung = {f"sigma_IF = {sig:+d}": [e for e in EPS if int(round(np.trace(kern_proj(sig, e) @ P_T16))) == 64]
             for sig in (1, -1)}
spur = {e: (zs(np.trace(auf16(ZXP[e]))/16), zs(np.trace(ZXP[e][np.ix_(UNG, UNG)])/16)) for e in EPS}
inv_mod = {e: bool(np.allclose(RT[-e] @ RT[e], Lam(R @ R))) for e in EPS}
ZAHLEN["Forderung und Uhrspur"] = {"16 bei T (Originalseite) für eps": forderung,
                                   "Tr_16 z_X^eps / 16, Tr_16q / 16": {e: list(v) for e, v in spur.items()}}
check("[L] z122s Forderung 'die 16, nicht die 16q, trägt R = +i auf der Originalgeraden' ist zirkulär: unter P_IF "
      f"(sigma_IF = +1) gilt sie genau für eps = +1, unter der Spiegelpaarung (sigma_IF = -1) genau für eps = -1 "
      f"({forderung}); sie ist mit P_IF "
      "gleichwertig zu P_glue (sigma_glue = +1), also das Chiralitätsvorzeichen selbst. Die Uhrspur des inneren Faktors kippt "
      f"mit dem Zweig: Tr_16 z_X^eps = 16 eps i, Tr_16q = -16 eps i ({spur}); z117s c = 64 x Chiralität x Schnitt x "
      "Uhrrichtung hat auf der 16 die Uhrrichtung eps, denn R~_(-eps) R~_eps = Lambda(R^2) ist eichartig, also "
      f"R~_(-eps) = R~_eps^-1 modulo G_SM ({inv_mod}): z118s Umkehrung 'R -> R^-1' ist auf der 16 der Zweigwechsel; die "
      "Marken legen R nur auf C^5 fest (1.1)",
      forderung == {"sigma_IF = +1": [1], "sigma_IF = -1": [-1]} and spur == {1: ("i", "-i"), -1: ("-i", "i")}
      and all(inv_mod.values()))

# --- 2.6 die eingefrorene Uhr der Viernaht-(E8)_1 gegen D5 + A3
WZ = []
for i, j in itertools.combinations(range(8), 2):
    for si, sj in itertools.product((1, -1), repeat=2):
        v = np.zeros(8)
        v[i], v[j] = si, sj
        WZ.append(v)
for s8 in itertools.product((0.5, -0.5), repeat=8):
    if sum(1 for x in s8 if x < 0) % 2 == 0:
        WZ.append(np.array(s8))
WZ = np.array(WZ)
zX8, zA8 = np.exp(1j*np.pi*WZ[:, :5].sum(1)), np.exp(1j*np.pi*WZ[:, 5:].sum(1))
kleber = float(np.abs(zX8*zA8 - 1).max())
ns = np.array(list(itertools.product(range(-2, 2), repeat=8)))
VS = []
for c in range(4):
    y4 = c + 4*ns                                                    # 4 y, y in (c/4) 1 + Z^8
    sel = ((y4**2).sum(1) == 32) & (((y4.sum(1)//4) + c) % 4 == 0)
    VS += [(c, v/4) for v in y4[sel]]
YV = np.array([v for _, v in VS])
CS = np.array([c for c, _ in VS])
gram = YV @ YV.T
schluessel = lambda v: tuple(np.round(4*v).astype(int))
RS8, RSV = {schluessel(v) for v in WZ}, {schluessel(v) for v in YV}
spiegel_zu = all(schluessel(b - (a @ b)*a) in RSV for a in YV[::7] for b in YV)
idx_c = [int(np.sum(CS == c)) for c in range(4)]


def mult(eig):
    """Vielfachheiten eines inneren Automorphismus auf der 248 (240 Wurzeln + 8 Cartan); andere Eigenwerte gesammelt."""
    lab = Counter()
    for e in eig:
        lab[zs(e)] += 1
    lab["1"] += 8
    out = {k: lab[k] for k in ("1", "i", "-1", "-i")}
    rest = sum(v for k, v in lab.items() if k not in out)
    if rest:
        out["sonst"] = rest
    return out


m_uhr = mult(np.exp(2j*np.pi*YV[:, 0]))                              # Sektorgradierung i^c (Original: c = 3, R = -i)
m_zx = mult(zX8)
H_EPS = {-1: np.array([-j/4 for j in range(5)] + [0, 0, 0]), 1: np.array([1 - 0/4] + [-j/4 for j in range(1, 5)] + [0, 0, 0])}
wD = lambda S: np.array([0.5 - (j in S) for j in range(5)])
d5_ok = all(abs(np.exp(2j*np.pi*H_EPS[e][:5] @ wD(S)) - RT[e][IDX[S], IDX[S]]) < 1e-12 for e in EPS for S in SUB)
m_rt = {e: mult(np.exp(2j*np.pi*WZ @ H_EPS[e])) for e in EPS}


def klammern(S, alle, plus):
    out = set()
    for a in S:
        for b in S:
            k = schluessel(a + b if plus else a - b)
            if k in alle:
                out.add(k)
    return len(out)


V_uhr = YV[CS == 3]
V16q = WZ[np.isclose(zX8, -1j)]
kl = {"Viernaht, R = -i": (klammern(V_uhr, RSV, True), klammern(V_uhr, RSV, False) + int(np.linalg.matrix_rank(V_uhr))),
      "(16q, 4q)": (klammern(V16q, RS8, True), klammern(V16q, RS8, False) + int(np.linalg.matrix_rank(V16q)))}
koeff = [Fr(0)] + [Fr(koeff113[k]) for k in ("R", "R2", "R3")]
proj = {zs(r): zs(sum(np.exp(2j*np.pi*float(c))*r**k for k, c in enumerate(koeff))/4) if abs(
    sum(np.exp(2j*np.pi*float(c))*r**k for k, c in enumerate(koeff))/4) > 1e-9 else "0" for r in (1, 1j, -1, -1j)}
min_norm = float(np.min(((3 + 4*ns)**2).sum(1))/16)                  # Originalsektor c = 3: y in -1/4 + Z^8
ZAHLEN["eingefrorene Uhr"] = {"Viernaht-Gitter: Wurzeln je Sektor c": idx_c, "ganzzahlig, Spiegelungen": [
    bool(np.allclose(gram, np.round(gram))), bool(spiegel_zu)], "Vielfachheiten (1, i, -1, -i)": {
    "Viernaht-Uhr": m_uhr, "z_X (D5 + A3)": m_zx, "R~_+ in E8": m_rt[1], "R~_- in E8": m_rt[-1]},
    "[V,V] und [V,V-quer]": {k: [int(x) for x in v] for k, v in kl.items()}, "z113-Projektor auf R": proj,
    "Minimalnorm des Originalsektors": min_norm}
check("[E, negativ] Eingefrorene Uhr (z113): Die Viernaht ist exakt E8 im Rahmen y in (c/4) 1 + Z^8, sum y = -c mod 4 "
      f"({len(YV)} Wurzeln, je Sektor {idx_c}, ganzzahlig und spiegelungsabgeschlossen); z113s Koeffizienten projizieren auf "
      f"R = -i ({proj}), das ist der Originalsektor c = 3 (Minimalnorm der Originalimpulse {min_norm}, kein E8-Vektor; erst "
      f"die Stufe 1 hat 64 Zustände). Die Viernaht-Uhr hat die Vielfachheiten {m_uhr} (Zentralisator u(8)), z_X von D5 + A3 "
      f"{m_zx} (so(10) + su(4), Kleber z_X z_A3 = 1 bis {kleber:.0e}), die gehobene Uhr R~_eps in E8 für beide Zweige "
      f"dieselben {m_rt[1]} (D5-Treue {d5_ok}): drei paarweise nicht konjugierte Elemente, die Zweige in E8 ununterscheidbar. "
      f"Die eingefrorene Schicht ist auch als Unterraum weder (16q, 4q) noch (16, 4): [V, V] und [V, V-quer] haben "
      f"{kl['Viernaht, R = -i']} gegen {kl['(16q, 4q)']} Dimensionen. Der Wert R = -i lässt sich auf die 16 nicht übertragen",
      len(YV) == 240 and idx_c == [56, 64, 56, 64] and np.allclose(gram, np.round(gram)) and spiegel_zu
      and proj == {"1": "0", "i": "0", "-1": "0", "-i": "1"} and min_norm == 0.5
      and m_uhr == {"1": 64, "i": 64, "-1": 56, "-i": 64} and m_zx == {"1": 60, "i": 64, "-1": 60, "-i": 64}
      and m_rt[1] == m_rt[-1] and m_rt[1] != m_uhr and m_rt[1] != m_zx and d5_ok and kleber < 1e-12
      and kl["Viernaht, R = -i"] == (56, 64) and kl["(16q, 4q)"] == (60, 60))

befund("Teil 2: Kein natives Datum wählt den Zweig. Alles, was über die Vektordarstellung läuft (Träger, Marken, R auf "
       "C^5, X, Y, beide 3+2-Schnitte, mod-16-Bilanz), ist blind; J_F und die Tomita-artige Konjugation sind für beide "
       "Zweige kovariant; die Spiegelung tau am Markenkreis kippt den Zweig bei festem R, festen Marken, Y, X und 16 und "
       "kehrt dabei nur die Hardy- und Energieorientierung um. Die Hardy-Orientierung (Index +1) wählt deshalb erst mit einer "
       "Gradierungsregel: Torusfortsetzung der leichtesten Energien des Originals eps = -1, leichteste Energie je Mode "
       "eps = +1. Die eingefrorene Uhr ist nicht übertragbar, z122s Forderung ist P_glue. Exakt gilt P_glue = P_IF x P_lift",
       True)

# =================================================================== Teil 3
say("\nTeil 3: Die fehlende Drei")

# --- 3.1 Original und Partner
CH8 = [1j, 1j, -1j, -1j, -1j, -1j, 1j, 1j]                          # a_o, a_o, b_o, b_o, a_p, a_p, b_p, b_p (R_phys)
st = np.arange(256)
bes = (st[:, None] >> np.arange(8)[None, :]) & 1
val = np.prod(np.where(bes == 1, np.array(CH8)[None, :], 1), axis=1)
par = bes.sum(1) % 2
chm = val[:, None]/val[None, :]
ung = (par[:, None] + par[None, :]) % 2 == 1
c_ung = sorted({zs(x) for x in chm[ung]})
c_ger = sorted({zs(x) for x in chm[~ung]})
spin_ok = np.allclose(val**2, (-1.0)**par)
check("[E, negativ] Original und Partner (je zwei Moden pro Quarter, Partner mit R_phys = R R-quer^-1, 256 Zustände, alle "
      f"65536 Matrixeinheiten): Ad_R hat auf ungeraden Operatoren genau {c_ung}, auf geraden {c_ger}, R^2 = (-1)^F überall "
      f"({spin_ok}). Der Partner liefert wie das Original nur die 2; eine fermionische Mode mit Charakter +-1 gibt es nicht",
      c_ung == ["-i", "i"] and c_ger == ["-1", "1"] and spin_ok)

# --- 3.2 Kompositionen über Nähte, Quarter mit alpha = 0, 1/2
komp, ok_komp = [], True
for N in (1, 2, 3):
    for k in itertools.product(range(4), repeat=N):
        chars = []
        for s_ in range(N):
            chars += [1j**k[s_], (-1j)**k[s_]]
        n_ = len(chars)
        st_ = np.arange(1 << n_)
        bes_ = (st_[:, None] >> np.arange(n_)[None, :]) & 1
        val_ = np.prod(np.where(bes_ == 1, np.array(chars)[None, :], 1), axis=1)
        par_ = bes_.sum(1) % 2
        ch_ = val_[:, None]/val_[None, :]
        odd_ = (par_[:, None] + par_[None, :]) % 2 == 1
        quelle = bool(np.any(np.isclose(ch_[odd_], 1) | np.isclose(ch_[odd_], -1)))
        bruch = not np.allclose(val_**2, (-1.0)**par_)
        ok_komp &= quelle == bruch and (not quelle) == all(x % 2 for x in k)
        komp.append(quelle)
alpha_tab = {}
for al in (Fr(0), Fr(1, 4), Fr(1, 2), Fr(3, 4)):
    ch = np.exp(2j*np.pi*float(al))
    alpha_tab[str(al)] = {"Charakter": zs(ch), "R^2 auf dem Fermion": zs(ch**2)}
ok_alpha = all((v["Charakter"] in ("1", "-1")) == (v["R^2 auf dem Fermion"] == "1") for v in alpha_tab.values())
ZAHLEN["Kompositionen"] = {"Uhren": len(komp), "ohne fermionische +-1": komp.count(False), "Quarter alpha": alpha_tab}
check(f"[E, negativ] Kompositionen: Über N = 1, 2, 3 Nähte und alle {len(komp)} Uhren R_k = prod R_s^(k_s) gibt es eine "
      "fermionische (ungerade) Mode mit Charakter +-1 genau dann, wenn R_k^2 = (-1)^F verletzt ist, und das genau, wenn ein "
      f"k_s gerade ist ({ok_komp}; ohne Quelle bleiben {komp.count(False)} = 2 + 4 + 8 Uhren, alle k_s ungerade). Ein Quarter "
      f"mit alpha = 0 oder 1/2 hätte Charakter +-1 und R^2 = +1 auf seinem Fermion ({alpha_tab}). Jede fermionische Quelle der "
      "Drei bricht also genau die Regel R^2 = (-1)^F, die z82 die Spinorhebung wählen lässt",
      ok_komp and komp.count(False) == 14 and ok_alpha)

# --- 3.3 Determinantenschranke (z82)
min_ger = min(sum(1 for m in range(5) if (m + d[m]) % 2 == 0)
              for d in itertools.product(range(4), repeat=5) if sum(d) % 4 == 0)
check(f"[E] Determinantenschranke (z82 nachgerechnet): Über alle 256 diagonalen mu4-Eichelemente mit det 1 bleibt mindestens "
      f"{min_ger} Charakter von R g mit Quadrat +1: Auch modulo G_SM sind nie alle fünf Trägermoden fermionisch, weil det R = -1",
      min_ger == 1)

# --- 3.4 die bosonische Drei (z92) und ihre Gradierungen
rel_prod = (sp.cancel(E_MOD[1]**2/E_MOD[2] - E_MOD[0]) == 0 and sp.cancel(E_MOD[1]*E_MOD[3]/E_MOD[2] - E_MOD[2]) == 0
            and sp.cancel(E_MOD[3]**2/E_MOD[2] - E_MOD[4]) == 0)
rel_eins = sp.cancel(E_MOD[4] - E_MOD[0] - 1) == 0
FEN = [[m + 4*k for k in range(-3, 4)] for m in range(5)]
prod_ok = lambda W: 2*W[2] == W[1] + W[3] and W[0] + W[2] == 2*W[1] and W[4] + W[2] == 2*W[3]
eins_ok = lambda W: W[0] == W[4]
GP, GE = [], []
for W in itertools.product(*FEN):
    if prod_ok(W):
        GP.append(W)
    if eins_ok(W):
        GE.append(W)
beide = [W for W in GP if eins_ok(W)]
zw = lambda W: zweig(hebung(W))
zp, ze = Counter(zw(W) for W in GP), Counter(zw(W) for W in GE)
z8 = {}
for W in GP:
    z8.setdefault(f"(w_1, w_3) mod 8 = ({W[1] % 8}, {W[3] % 8})", set()).add(zw(W))
z8 = {k: sorted(v) for k, v in sorted(z8.items())}
leicht_p = [W for W in GP if min(W) >= 0 and (W[1], W[3]) == (w1, w3)]
leicht_e = min((W for W in GE if min(W) >= 0 and (W[1], W[3]) == (w1, w3)), key=lambda W: (sum(W), W))
ZAHLEN["Gradierungen der Drei"] = {"produktverträglich": len(GP), "einsverträglich": len(GE), "beides": len(beide),
                                    "Zweige produkt / eins": [{e: zp[e] for e in EPS}, {e: ze[e] for e in EPS}],
                                    "Z8-Daten der 2 -> Zweige": z8, "leichteste produkt / eins": [leicht_p, leicht_e]}
check("[E] Die bosonische Drei (z92): e_1^2/e_2 = e_0, e_1 e_3/e_2 = e_2, e_3^2/e_2 = e_4 exakt, und e_4 - e_0 = 1 (die "
      f"Einsrelation des Markenschnitts) exakt ({rel_prod}, {rel_eins}). Ganzzahlige Gradierungen, die die Charaktere heben "
      f"(Fenster m + 4{{-3..3}}): produktverträglich {len(GP)} (Zweige {dict(sorted(zp.items()))}), einsverträglich {len(GE)} "
      f"({dict(sorted(ze.items()))}), beides {len(beide)}: die Produktrelationen erzwingen w_4 - w_0 = 2(w_3 - w_1) = 4 mod 8. "
      f"Mit den leichtesten Energien des Originals auf der 2 ist die produktverträgliche Fortsetzung eindeutig {leicht_p} "
      f"(eps = {zt(zw(leicht_p[0])) if leicht_p else None}, Torus, Filtration an 0), die leichteste einsverträgliche "
      f"{leicht_e} (eps = {zt(zw(leicht_e))}). Die Z8-Daten der 2 (Halbuhr e^(i pi L_0) der Ordnung 8, z97) lassen beide "
      f"Zweige offen ({z8}); produktverträglich hängt eps an w_1 + w_3 mod 16. Der Zweig liegt jenseits der nativen "
      "Z8-Verfeinerung: erst ganzzahlige Energien und eine Gradierungsregel für die Drei legen ihn fest",
      rel_prod and rel_eins and len(beide) == 0 and zp[1] > 0 and zp[-1] > 0 and ze[1] > 0 and ze[-1] > 0
      and leicht_p == [(0, 1, 2, 3, 4)] and zw(leicht_p[0]) == -1 and leicht_e == (0, 1, 2, 3, 0) and zw(leicht_e) == 1
      and all(v == [-1, 1] for v in z8.values()))

befund("Teil 3: Die Drei bekommt keine native fermionische Quelle. Original und Partner liefern auf ungeraden Operatoren nur "
       "+-i; in jeder Komposition über Nähte ist eine fermionische +-1-Mode gleichbedeutend mit dem Bruch von "
       "R^2 = (-1)^F, Quarter mit alpha = 0, 1/2 ebenso, und det R = -1 lässt modulo G_SM stets eine R^2-gerade Mode. Bosonisch "
       "existiert die Drei (z92: Sym^2(2)/(uv)^2), aber keine Gradierung respektiert zugleich ihre Produktrelationen und die "
       "Einsrelation des Markenschnitts, und genau dort sitzt der Zweig: produktverträglich und leichteste Energien eps = -1, "
       "einsverträglich und leichteste Energien eps = +1", True)

# =================================================================== Teil 4
say("\nTeil 4: Bilanz und Urteil")
DATEN = {
    "X-Ladung und 3+2-Schnitt": "blind (Vektordarstellung); U1-Regel = Orientierung von U(1)_X",
    "Hardy-Orientierung / Index von A_q": "wählt nur mit Gradierungsregel: Torus eps = -1, je Mode eps = +1",
    "Tomita J, J_F, K": "kovariant für beide Zweige",
    "Spiegelung tau am Markenkreis": "kippt eps bei festem R, Marken, Y, X, 16; kehrt Hardy- und Energieorientierung um",
    "eingefrorene Uhr R = -i (z113)": "nicht übertragbar (u(8) gegen so(10) + su(4); R~_+ und R~_- in E8 gleich)",
    "z122-Forderung (16 mit R = +i auf der Originalgeraden)": "zirkulär: mit P_IF gleich P_glue",
}
ohne_regel = [k for k, v in DATEN.items() if "wählt nur" not in v and "zirkulär" not in v]
kritik = (CHECKS[[k for k in CHECKS if k.startswith("[E, negativ] X-Ladung")][0]]
          and CHECKS[[k for k in CHECKS if k.startswith("[C] Hardy")][0]]
          and CHECKS[[k for k in CHECKS if k.startswith("[E] Verklebter")][0]]
          and CHECKS[[k for k in CHECKS if k.startswith("[E, negativ] Eingefrorene")][0]]
          and CHECKS[[k for k in CHECKS if k.startswith("[E] Die bosonische Drei")][0]])
check(f"[C] Auftragskriterien: {len(DATEN)} native Daten geprüft, jedes mit Gegenprobe im anderen Zweig; keines wählt den "
      f"Zweig ohne Zusatzidentifikation ({len(ohne_regel)} blind, kovariant, kippend oder nicht übertragbar; zwei wählen nur "
      "mit Regel bzw. zirkulär), die beiden natürlichen Regeln geben entgegengesetzte Zweige, und der Zweig ist exakt das "
      "Relativvorzeichen sigma_glue sigma_IF",
      kritik and len(ohne_regel) == 4 and len(DATEN) == 6)
STATUS = {
    "P_lift als endliche Aussage": ("E", "R~_eps = eps i Lambda(R), eps = +1; gleichwertig nu^c-Ladung, tr W mod 8, "
                                         "Hauptwinkel von e_2, det des (u, v)-Lifts, U(1)_X-Halbdrehung"),
    "Existenz der Spin-Z4 (zweigfrei)": ("E", "beide Zweige modulo G_SM = z_X^(+-1); R^2 = (-1)^F wählt die Spinorhebung; "
                                              "Vektordarstellung zweigblind"),
    "X-Ladung, 3+2-Schnitt": ("E (negativ)", "blind; U1-Regel ist Orientierung von U(1)_X"),
    "Hardy-Orientierung": ("C", "Torusfortsetzung der leichtesten Energien eps = -1 (gegen die Konvention), leichteste "
                                "Energie je Mode eps = +1"),
    "Tomita, J_F, K": ("E (negativ)", "kovariant für beide Zweige; z118s 'J R J^-1 = R' für beide"),
    "Spiegelung tau": ("E", "kippt eps, hält R, Marken, Y, X, 16, kehrt Hardy- und Energieorientierung um"),
    "eingefrorene Uhr (z113)": ("E (negativ)", "nicht auf die 16 übertragbar; Zweige in E8 gleich"),
    "z122-Forderung": ("L", "zirkulär, unter P_IF gleich P_glue"),
    "P_glue = P_IF x P_lift": ("E", "sigma_glue = sigma_IF eps; Theta und Zweigwechsel wirken auf dem Modul gleich"),
    "fehlende Drei": ("E (negativ)", "keine fermionische Quelle ohne Bruch von R^2 = (-1)^F; bosonisch ohne verträgliche "
                                     "Gradierung"),
    "Herkunft des Zweigs": ("O", "Datum der Drei jenseits der nativen Z8-Verfeinerung; kein natives Datum; Kandidat: Z8 der "
                                 "Verzweigungsfelder an den Marken (z92)"),
    "gesamt": ("irreduzibel", "nicht nativ, nicht gleich P_IF; exakt das Relativbit von P_glue gegen P_IF"),
}
for k_, (s2, t2) in STATUS.items():
    say(f"   [{s2}] {k_}: {t2}")
ZAHLEN["status"] = STATUS
ZAHLEN["native_Daten"] = DATEN
ZAHLEN["z132_urteil_kurz"] = {
    "Urteil": "P_lift irreduzibel: nicht nativ, nicht dasselbe Bit wie P_IF; exakt P_glue = P_IF x P_lift",
    "Prämissenliste": "P_lift bleibt; Inhalt geschärft: (i) zweigfrei die Uhr auf der unquelligen Drei, (ii) Zweig eps = "
                      "sigma_glue sigma_IF; die Buchung P_IF = P_glue setzt eps = +1 voraus",
    "Chiralität": "Existenz (eine Weyl-16, N16 = +-3) zweigunabhängig; Vorzeichen (Seite der 16, N16, eps-Term der Uhrspur) "
                  "= sigma_IF eps",
    "Konvention": "Bestand eps = +1; Torusfortsetzung der leichtesten Energien gäbe eps = -1"}
urteil = ("P_lift ist nicht nativ und nicht auf P_IF reduzierbar; es ist ein eigenes Bit, das sich exakt als Relativvorzeichen "
          "schreiben lässt: eps = sigma_glue sigma_IF, also P_glue = P_IF x P_lift. (1) Endliche Aussage [E]: Die Uhr des "
          "Originals wirkt auf S = Lambda^* C^5 (x) det^(-1/2) als R~_eps = eps i Lambda(R), eps = +1. Beide Zweige haben "
          "Ordnung 4, erfüllen R~^2 = -Lambda(R^2) (Spinregel modulo G_SM), sind z_X^eps mal Eichelement und wirken auf den "
          "zehn Majoranas gleich; der Zweigwechsel ist das zentrale -1 von Spin(10). Das Urbild von <R> ist Z4 x Z2 mit genau "
          "zwei Schnitten. Gleichwertig ist eps i die Uhrladung des G_SM-Singuletts nu^c, tr W = 6 mod 8 für jeden ganzzahligen "
          "Logarithmus, die Halbdrehung -pi der R = -1-Mode e_2 im Hauptwinkel, det = -i des Ordnung-4-Lifts auf (u, v) und "
          "die U(1)_X-Halbdrehung -pi. Zweigfrei und exakt ist die Spin-Z4 selbst (z82). (2) Native Daten: Alles, was über die "
          "Vektordarstellung läuft (Träger, Marken, R auf C^5, X, Y, beide 3+2-Schnitte, mod-16-Bilanz), ist blind; U1s "
          "Regel 'Zweig aus X' ist die Orientierung von U(1)_X. J_F und die Tomita-artige Konjugation sind für beide Zweige "
          "kovariant, auch z118s Begründung 'J R J^-1 = R'. Die Spiegelung tau: z -> 1/z-quer am Kreis der Marken hält R, die "
          "Marken, Y, X und die 16 fest und kippt eps; sie kehrt nur die Hardy-Orientierung (Index von A_q) und das "
          "Energievorzeichen um. Darum wählt die Hardy-Orientierung den Zweig erst mit einer Gradierungsregel, und die beiden "
          "natürlichen Regeln widersprechen sich: Die eindeutige Torusfortsetzung der leichtesten Energien des Originals "
          "(1/4 für R = i, 3/4 für R = -i) gibt die Gewichte (0, 1, 2, 3, 4) und eps = -1, entgegen der Bestandskonvention; "
          "die leichteste Energie je Mode gibt (0, 1, 2, 3, 0) und eps = +1. Die eingefrorene Uhr R = -i der "
          "Viernaht-(E8)_1 ist nicht übertragbar: Viernaht-Uhr (Zentralisator u(8)), z_X (so(10) + su(4)) und die gehobene "
          "Uhr (in E8 für beide Zweige gleich) sind paarweise nicht konjugiert, und die eingefrorene Stufe-1-Schicht ist "
          "weder (16, 4) noch (16q, 4q). z122s Forderung '16 mit R = +i auf der Originalgeraden' ist unter P_IF gleich "
          "P_glue, also zirkulär. (3) Die fehlende Drei hat keine native fermionische Quelle: Partner und Original liefern auf "
          "ungeraden Operatoren nur +-i; in allen 84 Uhren über bis zu drei Nähte gibt es eine fermionische +-1-Mode genau "
          "dann, wenn R^2 = (-1)^F verletzt ist; Quarter mit alpha = 0, 1/2 ebenso; det R = -1 lässt modulo G_SM stets eine "
          "R^2-gerade Mode. Bosonisch existiert sie (z92), aber keine ganzzahlige Gradierung respektiert zugleich ihre "
          "Produktrelationen und die Einsrelation e_4 - e_0 = 1 des Markenschnitts, und die Z8-Daten der 2 (Halbuhr, z97) "
          "lassen beide Zweige offen: Der Zweig ist ein Datum der Drei jenseits der nativen Z8-Verfeinerung. (4) Urteil: "
          "irreduzibel. Folgen: Die Existenz der Chiralität "
          "(eine Weyl-16, N16 = +-3) hängt nicht am Zweig. Ihr Vorzeichen (16 auf der Original- oder Partnergeraden, N16 = "
          "+3 oder -3, Vorzeichen des eps-Terms der Uhrspur) ist sigma_IF eps; z122s 'N16 = +3 auf der Originalgeraden' und "
          "die Buchung 'P_IF = P_glue, ein Bit' setzen eps = +1 stillschweigend voraus. In der Prämissenliste bleibt P_lift "
          "stehen, mit geschärftem Inhalt: zweigfrei die Uhr auf der unquelligen Drei, dazu der Zweig als Relativbit von "
          "P_glue gegen P_IF; wer P_IF und P_glue als ein Bit zählt, muss P_lift behalten. Nichts ist befördert; "
          "3+1-Aussagen bleiben Kinematik der Einpaarschicht.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z132_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
