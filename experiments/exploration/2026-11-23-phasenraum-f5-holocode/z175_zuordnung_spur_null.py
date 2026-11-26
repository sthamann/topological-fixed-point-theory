"""Z175: Schritt 2 der Kernstrategie aus z173. Zuordnung: P_lift gegen Sp R = 0, Sp R = 0 gegen A*, und I_step.

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Bausteine bleiben unberührt, nichts wird befördert.
Statuskürzel: E exakt, N numerisch, C bedingt, L Lesart, O offen.
Bestand: z82 (P_lift: Vierteldrehung auf dem Träger H0(O(D)), R = diag(1, i, -1, -i, 1) = (-1_5) g mit g in G_SM,
Spinorhebung = z_X mal Eichelement, zwei Zweige), z83 (I_step: der verdoppelte Diagonalschritt ist der native
Einzelschritt; su(2,2) schließt nur für die CCR-Normierung; die Uhr ist auf Schritten und Schalen trivial), Z34/Z35
(die acht Schnitte mit den Deckhälften 5 + 3 sind reell die 16 Nahteltern (10,1) + (1,6)), z165 (Volumenaufhebung
Sp_248 R = 0), z169 (Uhr R = exp(2 pi i v.P), v = (1/2^5, 0^3)), z170 (138 Klassen der Ordnung 2 bis 8 in E8; 'Spur null
und Ordnung höchstens 4' wählt die Uhr), z168 (<K>(tau) = -E_6/(3 E_4), <T> = 0 heißt tau = i).

Teil A  P_lift gegen Sp R = 0: was 'Spur null' auswählt, was daran Zusatz ist, die fermionische Normalform, der
        Zentralisatorsatz, die Markenzahl, die exakte Vierteldrehung auf den acht Schnitten, und was P_lift genau annimmt.
Teil B  Folgt Sp R = 0 aus <T> = 0? Uhr-eingesetzte und Uhr-verdrillte Charaktere aller 138 Klassen bei tau = i
        (mpmath), S-Dualität, lineare Antwort, Grundniveau-Bilanz.
Teil C  I_step: Blindheit der Uhr und der Sektoren für die Schrittnormierung; der Kandidat 'Spinorindex = A3-Spinor'.
Teil D  Bilanz: was ersetzt wird und was bleibt.
"""
import json
import sys
import time
from collections import Counter
from fractions import Fraction as Fr
from itertools import product
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.linalg
import sympy as sp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parents[2] / "research" / "originlab"))
import clock_sector_search as cs  # noqa: E402  (nur lesend)

OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T_START = time.time()
mp.mp.dps = 30


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def zeit():
    return f"(Laufzeit bisher {time.time() - T_START:.0f} s)"


ROOTS = cs.ROOTS
R2 = np.array([[int(2 * x) for x in r] for r in ROOTS], dtype=np.int64)           # doppelte Wurzeln, ganzzahlig
GANZ = np.array([all(x.denominator == 1 for x in r) for r in ROOTS])              # die 112 Wurzeln von so(16)
H, Q4, NULL = Fr(1, 2), Fr(1, 4), Fr(0)


def mult(v, N):
    """Vielfachheiten von exp(2 pi i v.) auf den 248 (None, wenn die Ordnung N nicht teilt)."""
    m = [8] + [0] * (N - 1)
    for r in ROOTS:
        x = cs.dot(v, r) * N
        if x.denominator != 1:
            return None
        m[int(x) % N] += 1
    return tuple(m)


def spur(m):
    N = len(m)
    return sum(c * sp.exp(2 * sp.pi * sp.I * k / N) for k, c in enumerate(m))


def spur_null(m):
    return cs._trace_zero(list(m))


def plus(a, b):
    return tuple(x + y for x, y in zip(a, b))


def hebung2(v):
    """Die zweite Hebung nach Spin(16): v_1 -> v_1 + 1 (multipliziert die Spinorwurzeln mit -1)."""
    return (v[0] + 1,) + tuple(v[1:])


# =================================================================== Teil A
say("Teil A: P_lift gegen Sp R = 0")
KLASSEN = cs.clock_candidates()
W = cs.coweights()
for c in KLASSEN:
    N, s = c["order"], c["kac"][1:]
    c["v"] = tuple(sum(Fr(s[i], N) * W[i][k] for i in range(8)) for k in range(8))
    c["spur0"] = "trace_zero" in c["clauses"]
SPURLOS = [c for c in KLASSEN if c["spur0"]]
UHRKLASSE = next(c for c in KLASSEN if c["order"] == 4 and c["spur0"])
UHRM = tuple(UHRKLASSE["multiplicities"])
NAME = {tuple(c["multiplicities"]): c["fixed"] for c in KLASSEN if c["order"] == 4}
ZAHLEN["A1_spurlos"] = [{"ordnung": c["order"], "fest": c["fixed"], "vielfachheiten": c["multiplicities"],
                         "grundniveau": str(c["twist_ground"])} for c in SPURLOS]
check("A1 [E] Bestand reproduziert (OriginLab clock_sector_search): 138 Klassen der Ordnung 2 bis 8 in E8; Sp_248 R = 0 "
      "haben genau fünf: " + "; ".join(f"Z{c['order']} {c['fixed']}" for c in SPURLOS) + ". Nur die Uhr (Z4, A3 + D5) hat "
      "keinen U(1)-Rest; die vier anderen haben Ordnung 6 oder 8 und einen U(1)-Rest",
      len(KLASSEN) == 138 and len(SPURLOS) == 5 and [c["order"] for c in SPURLOS] == [4, 6, 6, 8, 8]
      and ["semisimple_fixed" in c["clauses"] for c in SPURLOS] == [True, False, False, False, False]
      and UHRM == (60, 64, 60, 64))

ORD4 = [c for c in KLASSEN if c["order"] == 4]
FP = [c for c in ORD4 if "square_is_fermion_parity" in c["clauses"]]
ZAHLEN["A2_ordnung4"] = [{"fest": c["fixed"], "vielfachheiten": c["multiplicities"],
                          "spur": int(c["multiplicities"][0] - c["multiplicities"][2]),
                          "quadrat_ist_fermionparitaet": c in FP} for c in ORD4]
check("A2 [E] Was an 'Spur null und Ordnung höchstens 4' Zusatz ist: Die Ordnung 4 ist nativ (mu_4 stabilisiert tau = i, "
      "H6 unter A*), ebenso R^2 = (-1)^F (z82). Ordnung 4 haben 7 Klassen, davon mit R^2 in der Klasse der Fermionparität "
      "(120 | 128) genau drei: " + ", ".join(f"{c['fixed']} (Sp = {c['multiplicities'][0] - c['multiplicities'][2]})"
                                              for c in FP) + ". 'Spur null' trifft also eine Auswahl unter drei; die "
      "Ordnungsschranke und 'kein U(1)-Rest' sind daneben nicht nötig",
      len(ORD4) == 7 and len(FP) == 3 and sorted(c["multiplicities"][0] - c["multiplicities"][2] for c in FP) == [0, 8, 64]
      and sum(c["spur0"] for c in ORD4) == 1)

# fermionische Normalform: Involution r auf den 16 reellen Eltern mit p Eigenwerten -1, und die komplexe Struktur
zeilen, ok = [], True
for p in range(0, 17, 2):
    v = tuple([H] * (p // 2) + [NULL] * (8 - p // 2))
    m = mult(v, 4) if (p // 2) % 2 else mult(v, 2)
    ordn = 4 if (p // 2) % 2 else (2 if p not in (0, 16) else 1)
    if (p // 2) % 2:
        sp_ = m[0] - m[2]
        quad = all((2 * cs.dot(v, r)) % 1 == (NULL if g else H) for r, g in zip(ROOTS, GANZ))
    else:
        sp_ = m[0] - m[1]
        quad = False
    formel = Fr((16 - 2 * p) ** 2 - 16, 2) if 0 < p < 16 else Fr(248)
    ok &= (sp_ == formel) and (quad == ((p // 2) % 2 == 1))
    zeilen.append({"p": p, "ordnung": ordn, "spur_248": int(sp_), "quadrat_ist_(-1)^F": quad})
vJ = (Q4,) * 7 + (-Q4,)
mJ = mult(vJ, 4)
ZAHLEN["A3_normalform"] = {"involutionen": zeilen, "komplexe_struktur": {"vielfachheiten": mJ, "spur": mJ[0] - mJ[2]}}
check("A3 [E] Fermionische Normalform: Für die Spin-Hebung einer Involution r der 16 reellen Eltern mit p Eigenwerten -1 "
      "ist Sp_248 R = ((16 - 2p)^2 - 16)/2 (p = 2, ..., 14; der Spinoranteil verschwindet), und R^2 = (-1)^F genau für "
      "p = 2 mod 4. Sp R = 0 heißt (n_+ - n_-)^2 = n_+ + n_- = 16, also die Teilung 10 | 6. p = 2 gibt Sp = 64 "
      "(D7 + U(1)), die komplexe Struktur Sp = 8 (A7 + U(1)); das sind die drei Klassen aus A2",
      ok and mJ == (64, 64, 56, 64) and [z["spur_248"] for z in zeilen if z["p"] in (2, 6, 10, 14)] == [64, 0, 0, 64])

# die gewählte Klasse: (16, 4), (10, 6), Zentrum
vU = (H,) * 5 + (NULL,) * 3
ph = [int(4 * cs.dot(vU, r)) % 4 for r in ROOTS]
fest = [r for r, k in zip(ROOTS, ph) if k == 0]
d5 = [r for r in fest if all(x == 0 for x in r[5:])]
a3 = [r for r in fest if all(x == 0 for x in r[:5])]
eins = [r for r, k in zip(ROOTS, ph) if k == 1]
pr5 = Counter(r[:5] for r in eins)
pr3 = Counter(r[5:] for r in eins)
zwei = [r for r, k in zip(ROOTS, ph) if k == 2]
lam_gerade = [w for w in product((H, -H), repeat=5) if sum(1 for x in w if x < 0) % 2 == 0]
zX = {sp.exp(2 * sp.pi * sp.I * sum(H * x for x in w)) for w in lam_gerade}       # exp(2 pi i v.w), v = (1/2)^5
zX_formel = {sp.I ** (5 - 2 * sum(1 for x in w if x < 0)) for w in lam_gerade}
check("A4 [E] Die gewählte Klasse liefert die Aussage von P_lift: Fest sind D5 (40 Wurzeln) + A3 (12); der Eigenraum "
      "R = i ist (16, 4): 16 verschiedene D5-Gewichte der Norm 5/4, je viermal, mit halbzahligen Koordinaten (Spinor), "
      "und 4 A3-Gewichte der Norm 3/4, je sechzehnmal; R = -1 ist (10, 6). R ist auf (16, 4) der Skalar i, liegt also im "
      "Zentrum Z4 von Spin(10) (und von SU(4)); auf Lambda^gerade C^5 ist R = i^(5-2p) = i, das z_X aus z82. Die zweite "
      "Hebung ist R^-1 (Skalar -i) in derselben Klasse: Der Zweig ist eine Konvention",
      len(fest) == 52 and len(d5) == 40 and len(a3) == 12 and len(eins) == 64 and len(zwei) == 60
      and len(pr5) == 16 and set(pr5.values()) == {4} and all(cs.dot(w, w) == Fr(5, 4) for w in pr5)
      and all(all(x.denominator == 2 for x in w) for w in pr5)
      and len(pr3) == 4 and set(pr3.values()) == {16} and all(cs.dot(w, w) == Fr(3, 4) for w in pr3)
      and zX == {sp.I} and zX_formel == {sp.I}
      and mult(hebung2(vU), 4) == UHRM
      and all((int(4 * cs.dot(hebung2(vU), r)) % 4) == (-k) % 4 for r, k in zip(ROOTS, ph)))

# Zentralisatorsatz: alle mu_4-Elemente des Torus von Spin(16) (beide Hebungen) gegen die 240 Wurzeln
A8 = np.array(list(product(range(4), repeat=8)), dtype=np.int64)
A8 = np.vstack([A8, A8 + np.array([4, 0, 0, 0, 0, 0, 0, 0])])                      # zweite Hebung
PH = A8 @ R2.T                                                                    # 8 v.r in Einheiten von 1/8 ... a.r2 = 8 v.r
gueltig = (PH % 2 == 0).all(axis=1)
A8, PH = A8[gueltig], PH[gueltig]
IDX = (PH // 2) % 4
n = np.stack([(IDX == k).sum(axis=1) for k in range(4)], axis=1)
n[:, 0] += 8
spur0 = (n[:, 0] == n[:, 2]) & (n[:, 1] == n[:, 3])
D5B = np.array([all(x == 0 for x in r[5:]) and g for r, g in zip(ROOTS, GANZ)])
A3B = np.array([all(x == 0 for x in r[:5]) and g for r, g in zip(ROOTS, GANZ)])
zentral = (IDX[:, D5B] == 0).all(axis=1) & (IDX[:, A3B] == 0).all(axis=1)
elemente = {bytes(row.astype(np.uint8)) for row in IDX[zentral]}
ordnung = lambda row: 1 if not row.any() else (2 if not (row % 2).any() else 4)   # noqa: E731
zent = sorted({(ordnung(row), tuple(int(x) for x in nn)) for row, nn in zip(IDX[zentral], n[zentral])})
fp_exakt = ((2 * IDX[:, GANZ]) % 4 == 0).all(axis=1) & ((2 * IDX[:, ~GANZ]) % 4 == 2).all(axis=1)
typ = lambda a: (int((a % 4 == 0).sum()), int((a % 2 == 1).sum()), int((a % 4 == 2).sum()))   # noqa: E731
typen0 = Counter(typ(a) for a in A8[spur0])
typen0_fp = Counter(typ(a) for a in A8[spur0 & fp_exakt])
ZAHLEN["A5_torus"] = {"mu4_elemente_gueltig": int(len(A8)), "verschiedene_zentrale": len(elemente), "zentrale": zent,
                      "spurlos_typen_(n_+1, n_+-i, n_-1)": {str(k): v for k, v in sorted(typen0.items())},
                      "spurlos_und_quadrat_(-1)^F_typen": {str(k): v for k, v in sorted(typen0_fp.items())}}
check("A5 [E] Zentralisatorsatz: Unter allen mu_4-Elementen des Torus (65536 mal zwei Hebungen, davon "
      f"{len(A8)} mit Ordnung höchstens 4 auf den 248) vertauschen mit beiden Deckhälften, also mit so(10) + so(6) der "
      "Eltern, genau vier verschiedene Elemente: 1, (-1)^F und die Uhr c mit c^-1 (Vielfachheiten "
      f"{[z[1] for z in zent]}). Ein Element, das mit einer Unteralgebra vollen Rangs vertauscht, liegt im Torus; der "
      "Zentralisator ist also die Gruppe Z4, die die Uhr erzeugt. Jedes Element der Ordnung 4, das mit den Basiswechseln "
      "beider Hälften vertauscht, ist die Uhr oder ihr Inverses, mit Sp R = 0 und R^2 = (-1)^F als Folge",
      len(elemente) == 4 and zent == [(1, (248, 0, 0, 0)), (2, (120, 0, 128, 0)), (4, (60, 64, 60, 64))]
      and sum(1 for row in IDX[zentral] if ordnung(row) == 4) > 0
      and all(spur0[i] and fp_exakt[i] for i in np.where(zentral)[0] if ordnung(IDX[i]) == 4)
      and len({bytes(row.astype(np.uint8)) for row in IDX[zentral] if ordnung(row) == 4}) == 2)

# Markenzahl: Deckhälften (n + 1 | n - 1) bei n Marken
zeilen, ok = [], True
for nm in range(1, 9):
    N, p = 4 * nm, 2 * (nm - 1)                                                   # reelle Eltern, Eigenwerte -1
    # Spur der Involution auf so(N): Paare (i < j) mit Vorzeichen r_i r_j
    r = [1] * (N - p) + [-1] * p
    tr = sum(r[i] * r[j] for i in range(N) for j in range(i + 1, N))
    ok &= tr == 8 - 2 * nm
    zeilen.append({"marken": nm, "haelften": (nm + 1, nm - 1), "spur_so(4n)": tr, "quadrat_ist_(-1)^F": (p // 2) % 2 == 1})
ZAHLEN["A6_marken"] = zeilen
check("A6 [E] Markenzahl: Bei n Marken auf der Sphäre sind die Deckhälften h0(O(D)) = n + 1 und h0(K(D)) = n - 1, also "
      "4n reelle Eltern. Die Deckinvolution (+1 | -1) hat auf so(4n) die Spur 8 - 2n. Sie verschwindet genau für n = 4, "
      "die Markenzahl aus Gauß-Bonnet, bei der auch c = 2n = 8 und E8 vorliegen; ihre Spin-Hebung hat das Quadrat (-1)^F "
      "genau für gerades n. Für die Uhr als gehobene Deckinvolution ist Sp R = 0 also keine eigene Forderung",
      ok and [z["marken"] for z in zeilen if z["spur_so(4n)"] == 0] == [4]
      and [z["marken"] for z in zeilen if z["quadrat_ist_(-1)^F"]] == [2, 4, 6, 8])

# die exakte Vierteldrehung auf den acht Schnitten
z = sp.Symbol("z")


def charaktere(nm):
    zeta = sp.exp(2 * sp.pi * sp.I / nm)
    traeger = [sp.simplify((z ** m / (z ** nm - 1)).subs(z, zeta * z) / (z ** m / (z ** nm - 1))) for m in range(nm + 1)]
    # Schnitte von K(D): z^m dz/(z^n - 1), dz -> zeta dz
    index = [sp.simplify(zeta * (z ** m / (z ** nm - 1)).subs(z, zeta * z) / (z ** m / (z ** nm - 1))) for m in range(nm - 1)]
    return zeta, traeger, index


regulaer = True
for nm in range(1, 9):
    zeta, tr_, ix_ = charaktere(nm)
    zahl = Counter(int(sp.nsimplify(sp.log(c) / sp.log(zeta))) % nm if nm > 1 else 0 for c in tr_ + ix_)
    regulaer &= all(zahl[k] == 2 for k in range(nm))
_, tr4, ix4 = charaktere(4)
vGEO = (NULL, Q4, H, -Q4, NULL, H, Q4, -Q4)                                       # (1, i, -1, -i, 1 | -1, i, -i)
vTR = vGEO[:5] + (NULL,) * 3
mGEO, mGEO2, mTR = mult(vGEO, 4), mult(hebung2(vGEO), 4), mult(vTR, 4)
quad_geo = all((2 * cs.dot(vGEO, r)) % 1 == (NULL if g else H) for r, g in zip(ROOTS, GANZ))
fest_geo_d5 = sum(1 for r in d5 if cs.dot(vGEO, r).denominator == 1)
ZAHLEN["A7_vierteldrehung"] = {"traeger": [str(c) for c in tr4], "index": [str(c) for c in ix4],
                               "acht_schnitte": list(mGEO), "nur_traeger": list(mTR),
                               "klasse_nur_traeger": NAME.get(mTR), "spur_nur_traeger": mTR[0] - mTR[2],
                               "feste_wurzeln_im_traeger_D5": fest_geo_d5,
                               "spurlose_torustypen": {str(k): v for k, v in sorted(typen0.items())}}
check("A7 [E] Die exakte Vierteldrehung z -> i z hat auf dem Träger die Charaktere (1, i, -1, -i, 1) und auf dem "
      "Familienindex (i, -1, -i); zusammen ist das zweimal die reguläre Darstellung von mu_4 (für n = 1, ..., 8 Marken "
      "allgemein zweimal die reguläre von mu_n). Ihre Hebung auf die 248 hat die Vielfachheiten (60, 64, 60, 64) in "
      "beiden Zweigen: Sie liegt in der Uhrklasse, Sp R = 0 gilt nativ. Der Träger allein gibt "
      f"{mTR} (Klasse {NAME.get(mTR)}, Sp = {mTR[0] - mTR[2]}): Erst mit der Familienhälfte verschwindet die Spur",
      regulaer and tr4 == [1, sp.I, -1, -sp.I, 1] and ix4 == [sp.I, -1, -sp.I]
      and mGEO == UHRM and mGEO2 == UHRM and mTR == (82, 56, 54, 56) and not spur_null(mTR))
check("A8 [E, negativ] Die exakte Vierteldrehung ist aber nicht die Uhr von H6: Ihr Quadrat ist nicht die Fermionparität "
      "der Eltern (es dreht acht der sechzehn Eltern), und ihr Zentralisator enthält vom Träger-so(10) nur "
      f"{fest_geo_d5} Wurzeln, also nicht G_SM (8 Wurzeln). In der Uhrklasse liegen drei Torustypen (Zahl der komplexen "
      "Eigenwerte +1, +-i, -1): (5, 0, 3) und (3, 0, 5) (Deckinvolution, 10 | 6), (2, 4, 2) (Vierteldrehung, regulär), "
      "(4, 4, 0) und (0, 4, 4). Nur der Decktyp hat R^2 = (-1)^F der Eltern. Sp R = 0 legt also die Klasse fest, nicht "
      "den Darsteller",
      (not quad_geo) and fest_geo_d5 == 6
      and set(typen0) == {(5, 0, 3), (3, 0, 5), (2, 4, 2), (4, 4, 0), (0, 4, 4)}
      and set(typen0_fp) == {(5, 0, 3), (3, 0, 5)})

# Gegenprobe: die anderen spurlosen Klassen und die reguläre Darstellung
v8 = tuple(Fr(k, 8) for k in range(8))                                            # reguläre Darstellung von mu_8 auf C^8
m8 = mult(v8, 8)
kl8 = [c for c in KLASSEN if c["order"] == 8 and tuple(c["multiplicities"]) == m8]
ZAHLEN["A8b_mu8"] = {"vielfachheiten": list(m8), "klasse": [c["fixed"] for c in kl8]}
check("A8b [E] Gegenprobe zu den vier anderen spurlosen Klassen: Die reguläre Darstellung von mu_8 auf den acht "
      "komplexen Eltern liegt in genau einer von ihnen ("
      + ", ".join(f"Z8 {c['fixed']}" for c in kl8) + "); sie verlangte acht Marken, also 16 komplexe Eltern statt 8. "
      "Die beiden Klassen der Ordnung 6 sind kein Vielfaches einer regulären Darstellung (16 ist nicht durch 6 "
      "teilbar). Aus der Markengeometrie mit c = 8 entsteht unter den fünf nur die Uhr",
      len(kl8) == 1 and kl8[0]["spur0"] and m8 == (30, 32, 30, 32, 30, 32, 30, 32) and 16 % 6 != 0)

# was P_lift genau annimmt: Nebenklasse der Vierteldrehung unter Basiswechseln der Hälften
zahl_f1, zahl_fn = Counter(), Counter()
for a in product(range(4), repeat=5):
    if sum(a) % 4:
        continue                                                                  # SU(5): Determinante 1
    vv = tuple(H + Fr(x, 4) for x in a)
    zahl_f1[mult(vv + (NULL,) * 3, 4)] += 1
    zahl_fn[mult(vv + (H, Q4, -Q4), 4)] += 1
neben_u3, neben_su3 = set(), set()
for a in product(range(4), repeat=5):
    if sum(a) % 4:
        continue
    for b in product(range(4), repeat=3):
        vv = plus(vGEO, tuple(Fr(x, 4) for x in a) + tuple(Fr(x, 4) for x in b))
        if not all(cs.dot(vv, r).denominator == 1 for r in d5 + a3):
            continue                                                              # nicht zentral für die Hälften
        for vl in (vv, hebung2(vv)):
            m = mult(vl, 4)
            neben_u3.add(m)
            if sum(b) % 4 == 0:
                neben_su3.add(m)
ZAHLEN["A9_nebenklasse"] = {"eichtorus_256_familie_trivial": {str(k): v for k, v in sorted(zahl_f1.items())},
                            "eichtorus_256_familie_nativ": {str(k): v for k, v in sorted(zahl_fn.items())},
                            "zentral_mod_SU5xU3": sorted(neben_u3), "zentral_mod_SU5xSU3": sorted(neben_su3)}
check("A9 [E] Was P_lift annimmt: In der Nebenklasse der Vierteldrehung unter den Basiswechseln der Hälften (Torus von "
      "SU(5) auf dem Träger, Torus von U(3) auf den Familien) liegen als zentrale Elemente genau 1, (-1)^F, c und c^-1; "
      "mit SU(3) statt U(3) bleiben nur 1 und (-1)^F, also keine Uhr. Die Uhr von P_lift ist das Element der Ordnung 4 "
      "dieser Nebenklasse, das mit beiden Hälften vertauscht; die Familienzahl-U(1) muss dabei zu den Basiswechseln "
      "zählen. Sp R = 0 ersetzt das nicht: In der Eich-Nebenklasse von c sind "
      f"{zahl_f1[UHRM]} von 256 Elementen spurlos (mit nativen Familiencharakteren {zahl_fn[UHRM]} von 256)",
      neben_u3 == {(248, 0, 0, 0), (120, 0, 128, 0), (60, 64, 60, 64)} and neben_su3 == {(248, 0, 0, 0), (120, 0, 128, 0)}
      and zahl_f1[UHRM] == 51 and zahl_fn[UHRM] == 60 and sum(zahl_f1.values()) == 256)
befund("Teil A: Sp R = 0 folgt nativ (Deckinvolution bei vier Marken; Vierteldrehung als zweimal reguläre Darstellung) "
       "und legt mit der Ordnung 4 die Klasse der Uhr fest, in der R das Zentrum Z4 von Spin(10) ist. P_lift fällt damit "
       "nicht weg: Es wird zur Aussage P_zentral 'die Uhr vertauscht mit den Basiswechseln beider Deckhälften'", True)
say(f"   {zeit()}")

# =================================================================== Teil B
say("\nTeil B: Folgt Sp R = 0 aus <T> = 0?")
NMAX = 14


def q_hoch(tau, x):
    return mp.exp(2j * mp.pi * tau * x)


def stuecke(zc, tau, verdrillt):
    """Vier eindimensionale Summen (Wert, q d/dq) für eine Koordinate: ganz, ganz mit (-1)^n, halb, halb mit e^{i pi n}."""
    out = []
    for halb in (0, 1):
        for vz in (0, 1):
            w = d = mp.mpc(0)
            for k in range(-NMAX, NMAX + 1):
                nn = mp.mpf(k) + (mp.mpf(1) / 2 if halb else 0)
                e = (nn + zc) ** 2 / 2 if verdrillt else nn ** 2 / 2
                t = q_hoch(tau, e)
                if not verdrillt:
                    t *= mp.exp(2j * mp.pi * nn * zc)
                if vz:
                    t *= mp.exp(1j * mp.pi * nn)
                w += t
                d += e * t
            out.append((w, d))
    return out


def theta(v, tau, verdrillt=False):
    """Theta-Summe über E8 mit eingesetzter Uhr (oder über E8 + v) und ihre Ableitung q d/dq."""
    st = [stuecke(mp.mpf(x.numerator) / x.denominator, tau, verdrillt) for x in v]
    wert = abl = mp.mpc(0)
    for s in range(4):
        f = [st[j][s][0] for j in range(8)]
        g = [st[j][s][1] for j in range(8)]
        p = mp.mpf(1)
        for x in f:
            p *= x
        dsum = mp.mpc(0)
        for j in range(8):
            t = g[j]
            for i in range(8):
                if i != j:
                    t *= f[i]
            dsum += t
        wert += p / 2
        abl += dsum / 2
    return wert, abl


def eta(tau):
    q = q_hoch(tau, 1)
    return q_hoch(tau, mp.mpf(1) / 24) * mp.nprod(lambda k: 1 - q ** k, [1, mp.inf])


def E2(tau):
    q = q_hoch(tau, 1)
    return 1 - 24 * mp.nsum(lambda k: k * q ** k / (1 - q ** k), [1, mp.inf])


def K_wert(v, tau, verdrillt=False):
    """<K> = q d/dq log (Theta/eta^8), K = L_0 - 1/3."""
    w, d = theta(v, tau, verdrillt)
    return d / w - E2(tau) / 3, w / eta(tau) ** 8


TI = mp.mpc(0, 1)
K0, Z0 = K_wert((NULL,) * 8, TI)
tau2 = mp.mpc("0.3", "0.8")
zU = [K_wert(vU, t)[1] for t in (TI, tau2)]
zEta = [(eta(2 * t) / eta(4 * t)) ** 4 for t in (TI, tau2)]
fehl = max(abs(a - b) / abs(b) for a, b in zip(zU, zEta))
check("B1 [N, 25 Stellen] Kontrolle: Ohne Uhr ist <K>(i) = 0 (E_6(i) = 0, also <T> = 0 bei tau = i) und Z(i) = "
      f"E_4/eta^8 = {mp.nstr(Z0.real, 12)}; mit der Uhr ist Sp(R q^(L_0 - 1/3)) = (eta(2 tau)/eta(4 tau))^4 an zwei "
      f"Punkten (relativer Fehler unter 1e-24)",
      abs(K0) < mp.mpf(10) ** -25 and abs(Z0 - mp.mpf("12") ) < 1e-20 and fehl < mp.mpf(10) ** -24)

# lineare Antwort und <K> bei tau = i für alle Klassen
T4 = mp.mpc(0, 4)
q4, qi = q_hoch(T4, 1), q_hoch(TI, 1)
zeilen, ok_lin, ok_dual, kmin, dual_max = [], True, True, mp.inf, mp.mpf(0)
for c in KLASSEN:
    m = c["multiplicities"]
    N = c["order"]
    sp_exakt = sum(mp.mpf(x) * mp.cos(2 * mp.pi * k / N) for k, x in enumerate(m))
    k4, _ = K_wert(c["v"], T4)
    lin = (k4 + mp.mpf(1) / 3) / q4
    ok_lin &= abs(lin - sp_exakt) < mp.mpf(10) ** -6
    ki, zi = K_wert(c["v"], TI)
    kt, zt = K_wert(c["v"], TI, verdrillt=True)
    dual_max = max(dual_max, abs(ki + kt), abs(zi - zt) / abs(zi))
    kmin = min(kmin, abs(ki))
    c["K_i"], c["sp"] = ki.real, sp_exakt
    zeilen.append((c["order"], c["fixed"], int(mp.nint(sp_exakt)), ki.real))
ok_dual = dual_max < mp.mpf(10) ** -20
ZAHLEN["B2_K_bei_i"] = {f"Z{c['order']} {c['fixed']} {c['kac']}": {"Sp_R": int(mp.nint(c["sp"])),
                                                                    "<K>_R(i)": mp.nstr(c["K_i"], 15)} for c in SPURLOS}
check("B2 [N] Lineare Antwort: <K>_R(q) = -1/3 + (Sp R) q + O(q^2) für alle 138 Klassen (bei tau = 4i gegen die exakte "
      "Spur, Abweichung < 1e-6). Sp R = 0 ist also gleichwertig damit, dass die Uhr-eingesetzte Energie keine Antwort "
      "erster Ordnung hat. Das ist eine Umformulierung, keine Herleitung", ok_lin)
check("B3 [N] <T> = 0 sektorweise gilt für keine Klasse: |<K>_R(i)| ist mindestens "
      f"{mp.nstr(kmin, 6)} über alle 138 Klassen; für die fünf spurlosen liegt es bei 1/3 ("
      + ", ".join(mp.nstr(c["K_i"], 9) for c in SPURLOS) + ")",
      kmin > mp.mpf("0.02") and all(abs(c["K_i"] + mp.mpf(1) / 3) < mp.mpf("1e-4") for c in SPURLOS))
check("B4 [N] <T> = 0 für das Paar gilt für jede Klasse: Der Uhr-verdrillte Sektor (Summe über E8 + v) hat bei tau = i "
      "denselben Wert und die entgegengesetzte Energie, <K>^R(i) = -<K>_R(i) (S-Dualität; größte Abweichung über 138 "
      f"Klassen unter 1e-20). Jede modulinvariante Summe ist bei tau = i stationär. Die Bedingung E_6(i) = 0 "
      "unterscheidet also keine Klasse und wählt Sp R = 0 nicht aus", ok_dual)

# neutraler Sektor (Summe über die Potenzen der Uhr) bei tau = i
zeilen_n, nmin = {}, mp.inf
for c in SPURLOS + [x for x in ORD4 if not x["spur0"]]:
    N = c["order"]
    zs = ds = mp.mpc(0)
    for k in range(N):
        vk = tuple(k * x for x in c["v"])
        w, d = theta(vk, TI)
        zs += w
        ds += d
    kn = (ds / zs - E2(TI) / 3).real
    nmin = min(nmin, abs(kn))
    zeilen_n[f"Z{N} {c['fixed']} {c['kac']}"] = mp.nstr(kn, 12)
ZAHLEN["B5_neutraler_sektor"] = zeilen_n
check("B5 [N] Auch der uhrneutrale Sektor (Summe über alle Potenzen von R, das Vakuummodul der festen Algebra) erfüllt "
      f"<K>(i) = 0 für keine der elf geprüften Klassen (fünf spurlose und die sechs übrigen der Ordnung 4; kleinster "
      f"Betrag {mp.nstr(nmin, 6)}); für die Uhr ist <K>_0(i) = {zeilen_n[f'Z4 {UHRKLASSE['fixed']} {UHRKLASSE['kac']}']}",
      nmin > mp.mpf("0.01") and len(zeilen_n) == 11)


# Grundniveau-Bilanz der verdrillten Sektoren
def naechster(x):
    """Kleinster Abstand^2 von x zum Gitter E8 (D8 und D8 + (1/2)^8), exakt."""
    best = None
    for sh in (NULL, H):
        y = [xi - sh for xi in x]
        f = [Fr(round(yi)) for yi in y]
        if sum(f) % 2:
            j = max(range(8), key=lambda i: abs(y[i] - f[i]))
            f[j] += 1 if y[j] > f[j] else -1
        d = sum((yi - fi) ** 2 for yi, fi in zip(y, f))
        best = d if best is None or d < best else best
    return best


bil = {}
for c in KLASSEN:
    N = c["order"]
    hs = [naechster(tuple(k * x for x in c["v"])) / 2 for k in range(N)]
    c["bilanz"] = sum(h - Fr(1, 3) for h in hs)
    c["hs"] = hs
bil_null = [c for c in KLASSEN if c["bilanz"] == 0]
ZAHLEN["B6_grundniveau"] = {f"Z{c['order']} {c['fixed']} {c['kac']}": {"h_k": [str(h) for h in c["hs"]],
                                                                        "summe(h_k - 1/3)": str(c["bilanz"])}
                            for c in SPURLOS}
check("B6 [E] Grundniveau-Bilanz: Für die Uhr sind die Grundniveaus der vier Sektoren h = 0, 3/8, 1/2, 3/8, also "
      "L_0 - 1/3 = -1/3, 1/24, 1/6, 1/24 mit Summe -1/12. Die Summe verschwindet für keine der fünf spurlosen Klassen ("
      + ", ".join(str(c["bilanz"]) for c in SPURLOS) + f") und für {len(bil_null)} der 138 Klassen. Eine Casimir-Bilanz "
      "der Uhrsektoren gibt Sp R = 0 nicht",
      UHRKLASSE["hs"] == [0, Fr(3, 8), Fr(1, 2), Fr(3, 8)] and UHRKLASSE["bilanz"] == Fr(-1, 12)
      and all(c["bilanz"] != 0 for c in SPURLOS))
t0 = mp.mpc("0.25", "0.25")
k_t0 = K_wert(vU, t0)[0]
ZAHLEN["B7_nullstelle"] = {"tau": "(1 + i)/4", "<K>_R": mp.nstr(abs(k_t0), 5)}
check("B7 [N; Lesart L] Die Uhr-eingesetzte Energie verschwindet bei tau = (1 + i)/4, wo 2 tau zum quadratischen Punkt "
      "äquivalent ist (|<K>_R| unter 1e-20; J(2 tau) = 1, J(tau) = "
      f"{mp.nstr(mp.kleinj(t0).real, 8)} mit J(i) = 1); dieser Punkt ist nicht zu tau = i äquivalent",
      abs(k_t0) < 1e-20 and abs(mp.kleinj(2 * t0) - 1) < 1e-20 and abs(mp.kleinj(t0) - 1) > 0.5)
befund("Teil B: Sp R = 0 folgt nicht aus <T> = 0. Sektorweise erfüllt keine Klasse die Bedingung, paarweise jede; die "
       "Bilanz der Grundniveaus verschwindet für keine spurlose Klasse. Die Spur ist die Antwort erster Ordnung, also "
       "eine eigene Aussage. Herleitbar ist Sp R = 0 aus der Markengeometrie (Teil A), nicht aus der Stationarität", True)
say(f"   {zeit()}")

# =================================================================== Teil C
say("\nTeil C: I_step")
M = 8
basis = [(m, n1, b1) for m in range(1, M + 1) for n1 in range(m) for b1 in range(m)]   # n2 = m-1-n1, nb2 = m-1-b1
pos = {s: i for i, s in enumerate(basis)}
DIM = len(basis)
schale = np.array([s[0] for s in basis])
NIEDRIG = schale <= M - 4


def schritte(F):
    """P_ab = F(Zielschale) a_a^+ b_b^+ auf gleichen Stufen (z83), K = P^T."""
    P = {}
    for a in (0, 1):
        for b in (0, 1):
            X = np.zeros((DIM, DIM))
            for (m, n1, b1), i in pos.items():
                if m == M:
                    continue
                na = n1 if a == 0 else m - 1 - n1
                nb = b1 if b == 0 else m - 1 - b1
                ziel = (m + 1, n1 + (a == 0), b1 + (b == 0))
                X[pos[ziel], i] = F(m + 1) * np.sqrt((na + 1) * (nb + 1))
            P[(a, b)] = X
    return P


def rang(mats):
    V = np.array([x[np.ix_(NIEDRIG, NIEDRIG)].ravel() for x in mats])
    s = np.linalg.svd(V, compute_uv=False)
    return int((s > 1e-9 * s[0]).sum())


def lie_wachstum(F, stufen=3):
    P = schritte(F)
    gen = list(P.values()) + [x.T for x in P.values()]
    alle, out = list(gen), [rang(gen)]
    for _ in range(stufen):
        neu = [x @ g - g @ x for x in alle for g in gen]
        kand = alle + neu
        # unabhängige Auswahl auf den niedrigen Schalen
        V = np.array([x[np.ix_(NIEDRIG, NIEDRIG)].ravel() for x in kand])
        _, r_, piv = scipy.linalg.qr(V.T, pivoting=True, mode="economic")
        k = int((np.abs(np.diag(r_)) > 1e-9 * abs(r_[0, 0])).sum())
        alle = [kand[i] for i in sorted(piv[:k])]
        out.append(k)
    return out, P


FAMILIEN = {"CCR (F = 1)": lambda m: 1.0, "global reskaliert (F = 1,3)": lambda m: 1.3,
            "Schalenfaktor F = 1 + 1/m": lambda m: 1.0 + 1.0 / m, "Schalenfaktor F = m^(-1/2)": lambda m: m ** -0.5}
wachs, eext, uhrspur = {}, {}, {}
for name, F in FAMILIEN.items():
    w, P = lie_wachstum(F)
    wachs[name] = w
    E_ext = sum(P[k].T @ P[k] - P[k] @ P[k].T for k in P) / 4
    dg = np.diag(E_ext)
    eext[name] = [round(float(dg[schale == m].mean()), 6) for m in range(1, M - 1)]
    eext[name + " (zweite Differenzen)"] = [round(eext[name][i + 1] - 2 * eext[name][i] + eext[name][i - 1], 6)
                                            for i in range(1, len(eext[name]) - 1)]
    # Uhr R = exp(2 pi i E) mit E = m: auf jeder Schale die Eins, für jede Normierung
    uhrspur[name] = [[int(round((np.exp(2j * np.pi * k * schale[schale == m] / 1.0)).sum().real)) for m in range(1, 5)]
                     for k in range(4)]
ZAHLEN["C1_lie"] = wachs
ZAHLEN["C1_E_ext"] = eext
check("C1 [E/N] Verdoppelte Schalenleiter (z83, zwei Oszillatoren je Seite, Schalen m <= 8, Auswertung auf m <= 4): Die "
      "von den acht Schritten P_ab, K_ab erzeugte Lie-Algebra hat die Dimension 15 (su(2,2)) für die CCR-Normierung und "
      f"ihre globale Reskalierung ({wachs['CCR (F = 1)']}, {wachs['global reskaliert (F = 1,3)']}) und wächst für "
      f"Schalenfaktoren ({wachs['Schalenfaktor F = 1 + 1/m']}, {wachs['Schalenfaktor F = m^(-1/2)']}); E_ext ist nur für "
      "konstantes F proportional zu m (CCR: E_ext = m). I_step ist eine Aussage über diesen Schalenfaktor",
      wachs["CCR (F = 1)"] == [8, 15, 15, 15] and wachs["global reskaliert (F = 1,3)"] == [8, 15, 15, 15]
      and wachs["Schalenfaktor F = 1 + 1/m"][-1] > 15 and wachs["Schalenfaktor F = m^(-1/2)"][-1] > 15
      and eext["CCR (F = 1)"] == [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
      and not any(eext["global reskaliert (F = 1,3) (zweite Differenzen)"])
      and any(eext["Schalenfaktor F = 1 + 1/m (zweite Differenzen)"])
      and any(eext["Schalenfaktor F = m^(-1/2) (zweite Differenzen)"]))
z83 = json.loads((HERE / "z83_ergebnis.json").read_text())["checks"]
uhr83 = [k for k, v in z83.items() if k.startswith("[E] Uhr: R = e^{2 pi i L_0}") and v]
satz83 = [k for k, v in z83.items() if k.startswith("[E] Satz über Schalenfaktoren") and v]
check("C2 [E, negativ] Die Uhr sieht den Schalenfaktor nicht: Auf der verdoppelten Leiter ist R = exp(2 pi i L_0) die "
      "Eins (Bestand z83, dort geprüft: Ad_R = 1 auf beiden Elementarschritten), und die Uhrspuren Sp(R^k) je Schale "
      "sind für alle vier Normierungen dieselben (m^2 für jedes k). Jede Bedingung, die nur Uhrwerte, Uhrspuren oder "
      "E8-Klassen enthält (Sp R = 0, Ordnung, Zentralität), gilt für alle Normierungen gleich und kann I_step nicht "
      "festlegen",
      len(uhr83) == 1 and len(satz83) == 1
      and all(u == [[1, 4, 9, 16]] * 4 for u in uhrspur.values()))
# Kandidat: Spinorindex der Schritte = Spinor von A3 (su(4) als kompakte Form von su(2,2))
ok_q = True
for qd in range(1, 5):
    for w in product(range(4), repeat=qd):                                        # Gewichte von 4^(x qd)
        m1 = sum(1 for x in w if x == 0) - sum(1 for x in w if x == 1)            # 2 j_1-Komponente
        m2 = sum(1 for x in w if x == 2) - sum(1 for x in w if x == 3)
        ok_q &= (m1 + m2) % 2 == qd % 2
check("C3 [E, negativ] Kandidat 'Spinorindex der Schritte = Spinor von A3' (su(4) als kompakte Form von su(2,2)): In "
      "jeder Darstellung von SU(4) mit Quadralität q gilt 2 j_1 + 2 j_2 = q mod 2 unter SU(2) x SU(2); ein einzelner "
      "Spinorindex (1/2, 0) liegt also nur in Darstellungen ungerader Quadralität, auf denen die Uhr +-i ist. Die "
      "Schritte sind uhrneutral (C2). Der Schrittindex ist damit kein A3-Spinor; uhrneutral wären nur Paare (1/2, 1/2) "
      "in der 15. Die Sektorstruktur D5 x A3 liefert I_step nicht",
      ok_q)
befund("Teil C: I_step folgt nicht aus derselben Quelle. Geltungsbereich: alle Bedingungen in Uhrwerten, Uhrspuren und "
       "E8-Klassen; Schalenfaktorfamilie der verdoppelten Einpaarleiter. I_step bleibt eine eigene Annahme, gleichwertig "
       "zum Abschluss von su(2,2) (z83)", True)
say(f"   {zeit()}")

# =================================================================== Teil D
say("\nTeil D: Bilanz")
bilanz = [
    "[E] Sp_248 R = 0 und Ordnung 4 legen die Klasse der Uhr fest; unter den drei Klassen mit R^2 = (-1)^F ist es eine "
    "Auswahl aus drei. Ordnungsschranke und 'kein U(1)-Rest' sind neben der nativen Ordnung 4 nicht nötig (A1, A2)",
    "[E] In dieser Klasse ist R das Zentrum Z4 von Spin(10), Skalar i auf (16, 4), z_X auf Lambda^gerade C^5; der Zweig "
    "ist die Wahl R oder R^-1 (A4)",
    "[E] Sp R = 0 ist nativ: Für die gehobene Deckinvolution der Hälften (n + 1 | n - 1) ist die Spur 8 - 2n, null genau "
    "bei vier Marken; die exakte Vierteldrehung ist zweimal regulär und liegt in der Uhrklasse (A3, A6, A7)",
    "[E] Zentralisatorsatz: Ordnung 4 und Vertauschen mit so(10) + so(6) der Hälften geben c oder c^-1 (A5)",
    "[E, negativ] Sp R = 0 legt den Darsteller nicht fest (drei Torustypen; 51 von 256 in der Eich-Nebenklasse); die "
    "exakte Vierteldrehung ist nicht die Uhr von H6 (A8, A9)",
    "[C] P_lift wird durch P_zentral ersetzt: 'die Uhr vertauscht mit den Basiswechseln beider Deckhälften, die "
    "Familienzahl eingeschlossen'. Hebungszweig und 'modulo G_SM' entfallen als eigene Wahlen; die Zahl der Annahmen "
    "bleibt gleich",
    "[N/E, negativ] Sp R = 0 folgt nicht aus <T> = 0: sektorweise falsch für alle 138 Klassen, paarweise wahr für alle; "
    "Grundniveau-Bilanz ungleich null (B3 bis B6)",
    "[E, negativ] I_step folgt weder aus Sp R = 0 noch aus der Sektorstruktur D5 x A3 (C2, C3)",
]
for b in bilanz:
    say("   " + b)
ZAHLEN["D_bilanz"] = bilanz
ZAHLEN["D_praemissen"] = {
    "P_lift": "ersetzbar durch P_zentral (gleiche Zahl; schärfer und ohne Zweigwahl)",
    "P_zentral": "die Uhr hat Ordnung 4 und vertauscht mit den Basiswechseln der Schnitte beider Deckhälften "
                 "(so(10) + so(6) der 16 Nahteltern, Familienzahl-U(1) eingeschlossen); gibt R = c^(+-1), Sp R = 0, "
                 "R^2 = (-1)^F und die Aussage von P_lift; nicht aus A* hergeleitet",
    "I_step": "bleibt; nicht aus Uhr- oder Sektordaten herleitbar",
    "Sp R = 0": "keine eigene Annahme: Folge von P_zentral bzw. der Markengeometrie; nicht Folge von <T> = 0"}
check("D1 [E] Bilanz der Zuordnung: Keine der beiden Annahmen entfällt. P_lift lässt sich gegen P_zentral tauschen "
      "(A5, A9); Sp R = 0 ist dabei Folge und keine Annahme (A6, A7); I_step bleibt (C2, C3). Der Abstand zur "
      "Komplettlösung ändert sich durch Schritt 2 nicht",
      all(CHECKS.values()))

urteil = (
    "Schritt 2 ist gerechnet. (1) P_lift gegen Sp R = 0: Die Spurbedingung legt zusammen mit der nativen Ordnung 4 die "
    "E8-Klasse der Uhr fest, und in dieser Klasse ist R genau das Zentrum Z4 von Spin(10) mit dem Skalar i auf der 16; der "
    "Hebungszweig ist die Wahl zwischen R und R^-1. Sp R = 0 ist dabei keine neue Forderung: Für die gehobene "
    "Deckinvolution der beiden Hälften ist die Spur 8 - 2n und verschwindet genau bei vier Marken, und die exakte "
    "Vierteldrehung wirkt auf den acht Schnitten als zweimal reguläre Darstellung und liegt in derselben Klasse. Die "
    "Spurbedingung legt aber nur die Klasse fest, nicht den Darsteller; die exakte Vierteldrehung ist nicht die Uhr von H6. "
    "Was P_lift tatsächlich annimmt, ist die Zentralität: Die Uhr vertauscht mit den Basiswechseln beider Deckhälften. "
    "Unter dieser einen Aussage (P_zentral) folgen c oder c^-1, Sp R = 0, R^2 = (-1)^F und die Wirkung auf der 16 exakt. "
    "P_lift wird also schärfer gefasst, entfällt aber nicht. (2) Sp R = 0 aus A*: nein. Die Bedingung <T> = 0 ist bei "
    "tau = i für jedes Paar aus Uhr-eingesetztem und Uhr-verdrilltem Sektor automatisch erfüllt und für keinen einzelnen "
    "Sektor, bei allen 138 Klassen; die Bilanz der Grundniveaus verschwindet für keine spurlose Klasse. Die Spur ist die "
    "Antwort erster Ordnung der eingesetzten Energie. (3) I_step: Die Uhr ist auf der verdoppelten Leiter die Eins und "
    "sieht den Schalenfaktor nicht; der Schrittindex ist kein A3-Spinor. I_step bleibt eigenständig. Der Abstand zur "
    "Komplettlösung ändert sich nicht. Keine Komplettlösung.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z175_ergebnis.json").write_text(json.dumps({"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil},
                                                    indent=1, ensure_ascii=False, default=str))
(HERE / "ausgabe_z175.txt").write_text("\n".join(x for x in OUT if "Laufzeit bisher" not in x) + "\n")
say(f"{zeit()}")
