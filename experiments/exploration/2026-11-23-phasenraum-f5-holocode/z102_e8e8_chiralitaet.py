"""Z102: Folgefrage N3 'Chiralität aus der Gegenstruktur: E8 + E8 mit Tomita-Kleber' (Hypothese H3).

Leitlinie 'Rekonstruktion statt Auswahl': eine Lesart zählt nur als nativ, wenn sie aus vorhandener TFPT-Struktur folgt.
Firewalls: DIMENSION.UPLIFT.FIREWALL.01 (jede 3+1-Aussage ist Kinematik der radialen Einpaarschicht oder Lesart),
GRAV.NONCIRCULAR.01 (kein Schluss modular => Raumzeit). Gelöste Bausteine (E8 = D5 + A3 über Z4, 16 je Generation, Flavour
mit dem Familiengitter A3, q(A3) = 3/4, Familiencharaktere {-1, i, -i}) werden benutzt, nicht neu geöffnet. Bestand: z97
(Kleber z_X z_A3 = 1 auf 240 Wurzeln, Zentrum der radialen SU(2,2) Z_SU = i^(2h-2), unter P_A3 R = Z_SU^-1 mit chiraler
Z4-Regel, aber Konflikt mit der neutralen Skalarleiter h = 0, R = 1 und dem Familiengitter), z95 (E8 + E8(-1) mit isotroper
Diagonale vom Rang 8, J R J = R^-1), z91/z89 (diagonale E8-Komposition), z83 (radiale Standardform), z23 ((16,4) =
(16,3) + (16,1); Familien 3 x 16 als direkte Summe), z12 (Familiencharaktere aus der freien mu4-Bahn der Marken).

Hypothese H3: das Original E8 = D5 + A3_fam bleibt; der Tomita-Partner E8' = D5' + A3' liest A3' als konforme Algebra
(Prämisse P_A3' nur auf dem Partner); der Tomita-Kleber paart die 16 des Originals mit dem Twistor 4 des Partners. Erhofft:
eine Helizität, Familien intakt, neutrale Skalarleiter erlaubt.

Teil 1  Konstruktion. Blöcke D5, A3, D5', A3' mit Z4-Diskriminanten aus expliziten Vektoren; alle isotropen Kleber-Codes
        in (Z4)^4 für E8 + E8 (definit) und E8 + E8(-1) (Signatur (8,8)); explizite Gitter (Basis, Gram-Determinante,
        Wurzelsysteme); Tomita J als Blocktausch mit und ohne Konjugation; die isotrope Diagonale aus z95 und ihre Klassen
        (= Tomita-gepaarter Code); Ganzzahligkeit der Paarung 16 <-> 4' mit und ohne Familie.
Teil 2  Die Prüfungen (a) Chiralität, (b) Familienklasse, (c) neutrale Skalarleiter, (d) GSO, (e) drei Familien, für den
        Tomita-gepaarten Code, erschöpfend für alle Codes und an Negativkontrollen (ein E8 mit P_A3 im Original, ein E8 ohne
        Kleber, direkte Summe, A3-Tausch, Diagonalcode, D16+); Kern auf der radialen Leiter (z83) für jede Z4-wertige
        Helizitätsfunktion.
Teil 3  Bilanz, Urteil für H3, Folgen für T4 und K16.
"""
import itertools
import json
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
import scipy.sparse as sps
import sympy as sp
from scipy.linalg import expm

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


# ------------------------------------------------------------------ Blöcke (Koordinaten doppelt skaliert, ganzzahlig)
ART = ("D", "A", "D", "A")                                         # D5, A3 = D3, D5', A3'
DIM = (5, 3, 5, 3)
OFF = (0, 5, 8, 13)
NAME = ("D5", "A3_fam", "D5'", "A3'_konf")


def kurzvektoren(n):
    """Alle Vektoren von D_n^* mit Norm <= 2 (doppelt skaliert): ganz (gerade Einträge) oder halbganz (ungerade)."""
    out = []
    for v in itertools.product(range(-2, 3), repeat=n):
        if all(x % 2 == 0 for x in v) or all(x % 2 == 1 for x in v):
            nn = Fr(sum(x*x for x in v), 4)
            if nn <= 2:
                out.append(v)
    return out


def klasse(v, art):
    """Z4-Klasse über exp(i pi Summe w) = i^(Summe v) (v doppelt skaliert): D5: z_X = i^a, A3: z_A3 = i^(-b)."""
    k = sum(v) % 4
    return k if art == "D" else (-k) % 4


norm = lambda v: Fr(sum(x*x for x in v), 4)
zs = lambda z: next(t for t, w in (("1", 1), ("i", 1j), ("-1", -1), ("-i", -1j)) if np.isclose(z, w))
KURZ = {}
for art, n in (("D", 5), ("A", 3)):
    tab = {k: {} for k in range(4)}
    for v in kurzvektoren(n):
        tab[klasse(v, art)].setdefault(norm(v), []).append(v)
    KURZ[art] = tab
qmin = {art: [min(KURZ[art][k]) for k in range(4)] for art in "DA"}
REP = {art: [sorted(KURZ[art][k][qmin[art][k]])[-1] for k in range(4)] for art in "DA"}
dot = lambda u, w: Fr(sum(x*y for x, y in zip(u, w)), 4)

# =================================================================== Teil 1
say("Teil 1: E8 + E8' mit Tomita-Kleber: Blöcke, alle Kleber-Codes, explizite Gitter, isotrope Diagonale")

# --- 1.1 Blöcke und der Kleber eines E8 (Reproduktion z97)
E8 = []                                                            # Standardkonstruktion wie z97, doppelt skaliert
for i, j in itertools.combinations(range(8), 2):
    for si, sj in itertools.product((2, -2), repeat=2):
        v = [0]*8
        v[i], v[j] = si, sj
        E8.append(tuple(v))
for vz in itertools.product((1, -1), repeat=8):
    if vz.count(-1) % 2 == 0:
        E8.append(vz)
E8 = [(v, (klasse(v[:5], "D"), klasse(v[5:], "A"))) for v in E8]  # Klassen (a, b) von D5 und A3 = D3
klassen_e8 = sorted({kl for _, kl in E8})
zahl_e8 = {f"{kl}": sum(1 for _, k in E8 if k == kl) for kl in klassen_e8}
kleber_e8 = all((a - b) % 4 == 0 for a, b in klassen_e8)
einfach = [(2, -2, 0, 0, 0, 0, 0, 0), (0, 2, -2, 0, 0, 0, 0, 0), (0, 0, 2, -2, 0, 0, 0, 0), (0, 0, 0, 2, -2, 0, 0, 0),
           (0, 0, 0, 2, 2, 0, 0, 0), (0, 0, 0, 0, 0, 2, -2, 0), (0, 0, 0, 0, 0, 0, 2, -2), (0, 0, 0, 0, 0, 0, 2, 2)]
klebvek = REP["D"][1] + REP["A"][1]


def zbasis(zeilen):
    """Ganzzahlige Zeilenstufenbasis des Z-Spanns (Hermite-artig, exakt)."""
    zeilen = [list(r) for r in zeilen if any(r)]
    basis = []
    for sp_ in range(len(zeilen[0]) if zeilen else 0):
        piv = [r for r in zeilen if r[sp_] != 0]
        rest = [r for r in zeilen if r[sp_] == 0]
        while len(piv) > 1:
            piv.sort(key=lambda r: abs(r[sp_]))
            p, neu = piv[0], [piv[0]]
            for r in piv[1:]:
                f = r[sp_] // p[sp_]
                r2 = [x - f*y for x, y in zip(r, p)]
                (neu if r2[sp_] != 0 else rest).append(r2)
            piv = neu
        if piv:
            basis.append(piv[0])
        zeilen = [r for r in rest if any(r)]
    return basis


def gram(basis, metrik):
    G = [[sum(Fr(x*y*m, 4) for x, y, m in zip(u, w, metrik)) for w in basis] for u in basis]
    return sp.Matrix(G)


B_e8 = zbasis(einfach + [klebvek])
G_e8 = gram(B_e8, [1]*8)
det_e8 = G_e8.det()
gerade_e8 = all(G_e8[i, i] % 2 == 0 for i in range(8)) and all(x.q == 1 for x in G_e8)
ZAHLEN["block"] = {"Diskriminantennorm D5 je Klasse": [str(x) for x in qmin["D"]],
                   "Diskriminantennorm A3 je Klasse": [str(x) for x in qmin["A"]],
                   "E8-Wurzeln je Klasse (a, b)": zahl_e8, "det E8 (Basis aus D5 + A3 + Klebevektor)": str(det_e8)}
check("[E] Blöcke aus expliziten Vektoren: die Z4-Klassen (über exp(i pi Summe w)) haben die Minimalnormen D5 "
      f"{[str(x) for x in qmin['D']]} und A3 = D3 {[str(x) for x in qmin['A']]}, also q(A3) = 3/4 für die Familienklasse. "
      f"E8 = D5 + A3: genau {len(E8)} Wurzeln der Norm 2 in den Klassen {zahl_e8} (Kleber a = b, z_X z_A3 = 1, wie z97), "
      f"Basis aus 8 einfachen Wurzeln von D5 + A3 und einem Klebevektor (16, 4): det = {det_e8}, gerade, unimodular",
      qmin["D"] == [0, Fr(5, 4), 1, Fr(5, 4)] and qmin["A"] == [0, Fr(3, 4), 1, Fr(3, 4)] and len(E8) == 240
      and all(norm(v) == 2 for v, _ in E8) and kleber_e8
      and zahl_e8 == {"(0, 0)": 52, "(1, 1)": 64, "(2, 2)": 60, "(3, 3)": 64} and det_e8 == 1 and gerade_e8)

# --- 1.2 alle isotropen Kleber-Codes in (Z4)^4 für E8 + E8 (s = +1) und E8 + E8(-1) (s = -1)
ALLE = list(itertools.product(range(4), repeat=4))
plus4 = lambda x, y: tuple((u + w) % 4 for u, w in zip(x, y))
SIG = lambda s: (1, 1, s, s)


def qform(x, s):
    return sum(sg*norm(REP[ART[i]][x[i]]) for i, sg in enumerate(SIG(s)))


def bform(x, y, s):
    return sum(sg*dot(REP[ART[i]][x[i]], REP[ART[i]][y[i]]) for i, sg in enumerate(SIG(s))) % 1


def spann(gens):
    S, rand = {(0, 0, 0, 0)}, [(0, 0, 0, 0)]
    while rand:
        neu = []
        for x in rand:
            for g in gens:
                y = plus4(x, g)
                if y not in S:
                    S.add(y)
                    neu.append(y)
        rand = neu
    return frozenset(S)


def codes(s):
    iso = [x for x in ALLE if qform(x, s) % 2 == 0]
    alle, rand = {frozenset({(0, 0, 0, 0)})}, [frozenset({(0, 0, 0, 0)})]
    while rand:
        neu = []
        for C in rand:
            for x in iso:
                if x not in C and all(bform(x, c, s) == 0 for c in C):
                    D = spann(list(C) + [x])
                    if D not in alle and all(qform(y, s) % 2 == 0 for y in D):
                        alle.add(D)
                        neu.append(D)
        rand = neu
    return alle


CODES = {s: codes(s) for s in (1, -1)}
MAXI = {s: [C for C in CODES[s] if len(C) == 16] for s in (1, -1)}
form_ok = all(qform(tuple(k*np.array(e)), 1) % 2 == (k*k*qform(e, 1)) % 2 for k in range(4)
              for e in ((1, 0, 0, 0), (0, 1, 0, 0)))


def gitterbasis(C, s):
    gen = []
    for i in range(4):
        n_ = DIM[i]
        einf = [[0]*n_ for _ in range(n_)]
        for j in range(n_ - 1):
            einf[j][j], einf[j][j + 1] = 2, -2
        einf[n_ - 1][n_ - 2], einf[n_ - 1][n_ - 1] = 2, 2
        for r in einf:
            v = [0]*16
            v[OFF[i]:OFF[i] + n_] = r
            gen.append(v)
    for x in C:
        gen.append([c for i in range(4) for c in REP[ART[i]][x[i]]])
    return zbasis(gen)


def wurzeln(C, s, nur=None):
    """Norm-2-Vektoren des Übergitters (definit) bzw. einseitige Norm-+-2-Vektoren (nur = 'orig' / 'partner')."""
    out = []
    for x in C:
        if nur == "orig" and x[2:] != (0, 0) or nur == "partner" and x[:2] != (0, 0):
            continue
        listen = [KURZ[ART[i]][x[i]] for i in range(4)]
        for ns in itertools.product(*[sorted(l_) for l_ in listen]):
            if nur is None and sum(ns) != 2 or nur == "orig" and (ns[0] + ns[1] != 2 or ns[2] + ns[3] != 0) \
                    or nur == "partner" and (ns[2] + ns[3] != 2 or ns[0] + ns[1] != 0):
                continue
            for teile in itertools.product(*[listen[i][ns[i]] for i in range(4)]):
                out.append(sum(teile, ()))
    return out


def komponenten(W):
    A = np.array(W, dtype=float)
    adj = np.abs(A @ A.T) > 1e-9
    frei, groessen = set(range(len(W))), []
    while frei:
        st = [frei.pop()]
        n_ = 1
        while st:
            i = st.pop()
            for j in np.nonzero(adj[i])[0]:
                if j in frei:
                    frei.remove(j)
                    st.append(j)
                    n_ += 1
        groessen.append(n_)
    return sorted(groessen)


gitter = {}
for s in (1, -1):
    for C in MAXI[s]:
        Bm = gitterbasis(sorted(C), s)
        Gm = gram(Bm, [1]*8 + [s]*8)
        eintrag = {"Ordnung": len(C), "det": int(Gm.det()), "gerade": all(Gm[i, i] % 2 == 0 for i in range(16))
                   and all(x.q == 1 for x in Gm), "Rang": len(Bm)}
        if s == 1:
            W = wurzeln(C, s)
            eintrag["Wurzeln"] = len(W)
            eintrag["Komponenten"] = komponenten(W)
        else:
            eintrag["einseitige Wurzeln (Original, Partner)"] = [len(wurzeln(C, s, "orig")), len(wurzeln(C, s, "partner"))]
        gitter[(s, C)] = eintrag
typ_def = sorted(str(gitter[(1, C)]["Komponenten"]) for C in MAXI[1])
alle_unimodular = all(abs(e["det"]) == 1 and e["gerade"] and e["Rang"] == 16 for e in gitter.values())
n_e8e8 = sum(1 for C in MAXI[1] if gitter[(1, C)]["Komponenten"] == [240, 240])
n_d16 = sum(1 for C in MAXI[1] if gitter[(1, C)]["Komponenten"] == [480])
mit_orig = {s: [C for C in MAXI[s] if (1, 1, 0, 0) in C or (1, 3, 0, 0) in C] for s in (1, -1)}
partner_e8 = all(frozenset(x for x in C if x[:2] == (0, 0)) in (spann([(0, 0, 1, 1)]), spann([(0, 0, 1, 3)]))
                 for s in (1, -1) for C in mit_orig[s])
groessen = {s: {k: sum(1 for C in CODES[s] if len(C) == k) for k in (1, 2, 4, 8, 16)} for s in (1, -1)}
ZAHLEN["codes"] = {"isotrope Codes je Ordnung (definit)": groessen[1], "isotrope Codes je Ordnung (Signatur (8,8))": groessen[-1],
                   "maximale Codes, Wurzelkomponenten (definit)": typ_def,
                   "maximale Codes mit intaktem Original-E8 (definit, (8,8))": [len(mit_orig[1]), len(mit_orig[-1])]}
check(f"[E] Alle isotropen Kleber-Codes von D5 + A3 + D5' + A3' (Diskriminante (Z4)^4, q aus den Blöcken, q(kx) = k^2 q(x) "
      f"geprüft): definit {len(CODES[1])} Codes ({groessen[1]}), Signatur (8,8) {len(CODES[-1])} Codes ({groessen[-1]}). "
      f"Die maximalen (Ordnung 16) geben explizit gerade unimodulare Gitter (Basis aus 16 einfachen Wurzeln und Klebevektoren, "
      f"|det| = 1, Rang 16, alle {len(gitter)} geprüft): definit {n_e8e8} mal E8 + E8 (Wurzelkomponenten 240 + 240) und "
      f"{n_d16} mal D16+ (eine Komponente 480), in Signatur (8,8) zehn Fassungen von II_(8,8). Das Original-E8 ist in "
      f"{len(mit_orig[1])} bzw. {len(mit_orig[-1])} davon intakt; sein Partner ist dann ein eigenes E8' (Kleber (1, +-1))",
      form_ok and alle_unimodular and n_e8e8 == 8 and n_d16 == 2 and len(MAXI[1]) == 10 and len(MAXI[-1]) == 10
      and partner_e8 and len(mit_orig[1]) == 4 and len(mit_orig[-1]) == 4)

# --- 1.3 Tomita J, isotrope Diagonale (z95) und der gepaarte Code
J_plus = lambda x: (x[2], x[3], x[0], x[1])                       # Tausch Original <-> Partner
J_minus = lambda x: tuple((-y) % 4 for y in (x[2], x[3], x[0], x[1]))   # Tausch mit Konjugation
PAAR = spann([(1, 1, 1, 1)])
PAAR_K = spann([(1, 1, 3, 3)])
kl4 = lambda w: tuple(klasse(w[OFF[i]:OFF[i] + DIM[i]], ART[i]) for i in range(4))
metrik_m = [1]*8 + [-1]*8
diag = [(v + v, kl4(v + v)) for v, _ in E8]                          # (v, v) in E8 + E8(-1)
anti = [(v + tuple(-x for x in v), kl4(v + tuple(-x for x in v))) for v, _ in E8]   # (v, -v)
D_arr = np.array([v for v, _ in diag], dtype=float)
A_arr = np.array([v for v, _ in anti], dtype=float)
rang_diag = [int(np.linalg.matrix_rank(D_arr)), int(np.linalg.matrix_rank(A_arr))]
iso_diag = max(float(np.abs(X @ np.diag(metrik_m) @ X.T).max()) for X in (D_arr, A_arr))
klassen_diag = sorted({kl for _, kl in diag})
klassen_anti = sorted({kl for _, kl in anti})
det_paar = {}
for s in (1, -1):
    Bp = gitterbasis(sorted(PAAR), s)
    det_paar[s] = int(gram(Bp, [1]*8 + [s]*8).det())
wurzeln_paar = len(wurzeln(PAAR, 1))
einseitig_paar = [len(wurzeln(PAAR, -1, "orig")), len(wurzeln(PAAR, -1, "partner"))]
paar_min = [v for v in itertools.product(*[KURZ[ART[i]][1][qmin[ART[i]][1]] for i in range(4)])]
zerlegt = all(norm(v[0] + v[1]) == 2 and norm(v[2] + v[3]) == 2 for v in paar_min)
voll_paar = {s: [C for C in MAXI[s] if PAAR <= C] for s in (1, -1)}
J_inv_paar = all(J_plus(x) in C and J_minus(x) in C for C in (PAAR, PAAR_K) for x in C)
ZAHLEN["tomita"] = {"Rang der Diagonale (v, v) und der Antidiagonale (v, -v)": rang_diag,
                    "max |Skalarprodukt| auf beiden in E8 + E8(-1)": iso_diag,
                    "Klassen der Diagonale": [str(k) for k in klassen_diag],
                    "Klassen der Antidiagonale": [str(k) for k in klassen_anti],
                    "det des gepaarten Gitters (definit, (8,8))":
                        [det_paar[1], det_paar[-1]], "Wurzeln des gepaarten Gitters (definit)": wurzeln_paar,
                    "einseitige Wurzeln (Original, Partner), (8,8)": einseitig_paar,
                    "Minimalvektoren der Klasse (1,1,1,1)": len(paar_min),
                    "maximale Codes, die den gepaarten Code enthalten (definit, (8,8))": [len(voll_paar[1]), len(voll_paar[-1])]}
check("[E] Tomita-Paarung = isotrope Diagonale: in E8 + E8(-1) spannen die 240 Vektoren (v, v) wie die 240 Vektoren "
      f"(v, -v) je einen isotropen Teilraum vom Rang {rang_diag} (max |Skalarprodukt| {iso_diag:.0e}, wie z95). Die "
      f"Kleberklassen der Diagonale sind genau {[''.join(map(str, k)) for k in klassen_diag]}, die der Antidiagonale "
      f"{[''.join(map(str, k)) for k in klassen_anti]}: die 16 des Originals (mit 4_fam) ist an die 16' und den 4' des "
      "Partners geklebt (bzw. an 16q', 4q'). Dieser gepaarte Code P = <(1,1,1,1)> (konjugiert <(1,1,3,3)>) ist unter J+ "
      "(Tausch) und J- (Tausch mit Konjugation) invariant, aber "
      f"nicht maximal: das Gitter D5 + A3 + D5' + A3' + P hat |det| = {abs(det_paar[1])} (definit) bzw. {abs(det_paar[-1])} "
      f"((8,8)), ist also nicht selbstdual. Es hat nur die {wurzeln_paar} Wurzeln der vier Blöcke (definit; in (8,8) "
      f"einseitig {einseitig_paar}); die 16 kommt nur in den {len(paar_min)} Minimalvektoren der Klasse (1,1,1,1) vor, und "
      "jeder davon ist Summe einer Original-Wurzel (16, 4_fam) und einer Partner-Wurzel (16', 4'). Selbstduale Hüllen von P: "
      f"{len(voll_paar[1])} definit, {len(voll_paar[-1])} in (8,8)",
      rang_diag == [8, 8] and iso_diag == 0.0 and set(klassen_diag) == set(PAAR) and set(klassen_anti) == set(PAAR_K)
      and J_inv_paar and abs(det_paar[1]) == 16
      and abs(det_paar[-1]) == 16 and wurzeln_paar == 104 and einseitig_paar == [52, 52] and len(paar_min) == 4096
      and zerlegt and len(voll_paar[1]) == 3 and len(voll_paar[-1]) == 3)

# --- 1.4 Ganzzahligkeit: 16 (x) 4_fam (x) 4'_konf als ein Teilchen?
tripel = {s: [str(qform(x, s)) for x in ((1, 1, 0, 1), (1, 1, 0, 3))] for s in (1, -1)}
kreuz = {s: [str(qform(x, s)) for x in ((1, 0, 0, 1), (1, 0, 0, 3))] for s in (1, -1)}
kreuz_orig = [str(bform((1, 1, 0, 0), x, s)) for s in (1, -1) for x in ((1, 0, 0, 1), (1, 0, 0, 3))]
codes_kreuz = {s: [C for C in CODES[s] if (1, 0, 0, 1) in C or (1, 0, 0, 3) in C] for s in (1, -1)}
kreuz_mit_orig = sum(1 for s in (1, -1) for C in codes_kreuz[s] if (1, 1, 0, 0) in C or (1, 3, 0, 0) in C)
C0 = spann([(1, 1, 0, 0), (0, 0, 1, 1)])
W_C0 = wurzeln(C0, 1)
orig_d = sorted({klasse(w[13:16], "A") for w in W_C0 if any(w[0:8])})   # Partnerklasse d einseitiger Originalwurzeln
ZAHLEN["ganzzahligkeit"] = {"Norm der Klasse (16, 4_fam, 4') bzw. (16, 4_fam, 4q') (definit)": tripel[1],
                            "dieselbe in (8,8)": tripel[-1], "Norm der Kreuzklasse (16, 4') / (16, 4q') (definit, (8,8))":
                                [kreuz[1], kreuz[-1]], "b((1,1,0,0), Kreuzklasse)": kreuz_orig,
                            "Codes mit Kreuzklasse (definit, (8,8))": [len(codes_kreuz[1]), len(codes_kreuz[-1])],
                            "davon mit intaktem Original-E8": kreuz_mit_orig}
check("[negativ, E] Ein Teilchen 16 (x) 4_fam (x) 4'_konf gibt es nicht: die Klasse (1,1,0,+-1) hat die Norm "
      f"{tripel[1][0]} (definit) bzw. {tripel[-1][0]} ((8,8)), ist also in keinem ganzzahligen Gitter. Die Kreuzpaarung "
      f"16 <-> 4' allein, (1,0,0,+-1), hat die Norm {kreuz[1][0]} (definit, erlaubt) bzw. {kreuz[-1][0]} ((8,8), verboten); "
      f"mit dem Kleber des Originals hat sie die Paarung b = {kreuz_orig[0]} != 0: unter den {len(codes_kreuz[1])} "
      f"definiten Codes mit Kreuzklasse ist keiner mit intaktem Original-E8 ({kreuz_mit_orig}). Einseitige Anregungen des "
      "Originals verschieben die Klasse um (a, b, 0, 0): alle 240 Originalwurzeln von E8 + E8' haben die Partnerklasse "
      f"d = {orig_d}, das Zentrum Z_SU' der Kommutante ist auf ihnen 1; der Partner sieht die Helizität der Originalmaterie "
      "nicht",
      tripel[1] == ["11/4", "11/4"] and tripel[-1] == ["5/4", "5/4"] and kreuz[1] == ["2", "2"]
      and kreuz[-1] == ["1/2", "1/2"] and set(kreuz_orig) == {"1/4"} and kreuz_mit_orig == 0
      and len(codes_kreuz[-1]) == 0 and orig_d == [0] and len(W_C0) == 480)
befund("Teil 1: Der Tomita-Kleber von H3 ist als Gitter genau die Klassenmenge der isotropen Diagonale (z95), P = <(1,1,1,1)>: "
       "J-invariant, aber nicht selbstdual (|det| = 16), ohne Ein-Teilchen-16; die 16 tritt nur als zweiseitiges Paar "
       "Original-Wurzel + Partner-Wurzel auf. Ein Teilchen 16 (x) 4_fam (x) 4' ist nicht ganzzahlig (Norm 11/4 bzw. 5/4); die "
       "Kreuzpaarung 16 <-> 4' allein geht nur definit und nur ohne intaktes Original-E8", True)

# =================================================================== Teil 2
say("\nTeil 2: Die Prüfungen (a) bis (e), erschöpfend über alle Codes, Kern auf der radialen Leiter, Negativkontrollen")

# 4D-Lesart mit P_A3' (nur Partner): R = i^a (Klasse von R = z_X, z82), Z_SU = z_A3' = i^(-d) = i^(2h - 2)
H2 = range(-3, 4)                                                  # 2h, |h| <= 3/2
erlaubt_h = lambda d: [str(Fr(z, 2)) for z in H2 if (z - 2 + d) % 4 == 0]
RS = ["1", "i", "-1", "-i"]


def eigenschaften(C):
    a1 = sorted({x[3] for x in C if x[0] == 1})
    D0 = sorted({x[3] for x in C if x[0] == 0})
    return {"chiral": len(a1) == 1 and a1[0] % 2 == 1,
            "b1": bool(a1) and (all((x[0] - x[1]) % 4 == 0 for x in C) or all((x[0] + x[1]) % 4 == 0 for x in C)),
            "b2": (1, 1, 0, 0) in C or (1, 3, 0, 0) in C,
            "c_streng": (0, 0, 0, 2) in C, "c_R": 2 in D0,
            "gso": all((x[0] - x[3]) % 2 == 0 for x in C), "d bei a = 1": a1, "D0": D0}


def tabelle(C):
    return {RS[r]: sorted({h for x in C if x[0] == r for h in erlaubt_h(x[3])}, key=Fr) for r in range(4)
            if any(x[0] == r for x in C)}


# --- 2.1 der Tomita-gepaarte Code P
eP = eigenschaften(PAAR)
tab_P = tabelle(PAAR)
eK = eigenschaften(PAAR_K)
tab_K = tabelle(PAAR_K)
z97_tab = {}                                                       # Kontrolle K1: ein E8, P_A3 im Original (d := b)
K1 = spann([(1, 1, 0, 1)])
K1 = frozenset(x for x in K1 if x[2] == 0)
z97_tab = tabelle(K1)
zentrum_gleich = all((x[1] - x[3]) % 4 == 0 for x in PAAR)        # z_A3 = Z_SU' auf allen Klassen
inhalt_P = {}
for k in range(4):
    inhalt_P[str((k, k, k, k))] = {"SO(10)": ["45+1", "16", "10", "16q"][k], "A3_fam": ["15+1", "4", "6", "4q"][k],
                                   "SO(10)'": ["45+1", "16", "10", "16q"][k], "A3'_konf": ["15+1", "4", "6", "4q"][k],
                                   "R": RS[k], "Z_SU": RS[(-k) % 4], "h": erlaubt_h(k)}
ZAHLEN["gepaart"] = {"Eigenschaften": dict(eP), "erlaubte h je R": tab_P, "Inhalt": inhalt_P,
                     "z97-Tabelle (ein E8, P_A3)": z97_tab}
inhalt_txt = "; ".join(f"({v['SO(10)']}, {v['A3_fam']}, {v['SO(10)' + chr(39)]}, {v['A3' + chr(39) + '_konf']}) mit "
                       f"R = {v['R']}, h in {v['h']}" for v in inhalt_P.values())
check("[E] Der gepaarte Code P unter P_A3' (nur Partner): (a) chiral, die 16 hat nur die Helizitätsklasse Z_SU = -i "
      f"(h in {erlaubt_h(1)}); (b) als Kleberrelation intakt, a = b auf allen Klassen, die Familienklasse hat q = 3/4; "
      f"(d) GSO R^2 = (-1)^(2h) auf allen Klassen. Inhalt: {inhalt_txt}. Die erlaubten h je R sind {tab_P}, Klasse für "
      f"Klasse identisch mit der z97-Regel des einen E8 unter P_A3 ({z97_tab}), denn auf P gilt z_A3 = Z_SU' (b = d auf "
      "allen Klassen): der Tomita-Kleber ist P_A3 auf Zentrumsebene. Die Skalare (h = 0) liegen nur im Sektor R = -1, als "
      f"(10, 6_fam, 10', 6'). Die Antidiagonale <(1,1,3,3)> gibt dasselbe mit der anderen Helizität ({tab_K})",
      eP["chiral"] and eP["b1"] and eP["gso"] and not eP["c_R"] and tab_P == z97_tab and zentrum_gleich
      and tab_P == {"1": ["-1", "1"], "i": ["-3/2", "1/2"], "-1": ["0"], "-i": ["-1/2", "3/2"]}
      and eK["chiral"] and eK["b1"] and eK["gso"] and not eK["c_R"] and eK["d bei a = 1"] == [3]
      and tab_K["i"] == ["-1/2", "3/2"])

# --- 2.2 Ein-Teilchen-Inhalt und diagonale SO(10)
spin_D5 = KURZ["D"][1][Fr(5, 4)]
diag_klassen = {"J+": sorted({klasse(tuple(x + y for x, y in zip(u, w)), "D") for u in spin_D5 for w in spin_D5}),
                "J-": sorted({klasse(tuple(x - y for x, y in zip(u, w)), "D") for u in spin_D5 for w in spin_D5})}
kopien = len(paar_min)//(len(spin_D5)*len(KURZ["A"][1][Fr(3, 4)]))
ZAHLEN["ein_teilchen"] = {"Klassen der diagonalen SO(10)-Gewichte der Paare (J+, J-)": diag_klassen,
                          "Kopien der Original-16 je Twistorkomponente (SO(10) x SO(10)')": kopien}
check("[negativ, E] Ein-Teilchen-Inhalt von P: keine 16 (definit nur die 104 Blockwurzeln, in (8,8) einseitig 52 + 52; die "
      "16 nur als Paar Original-Wurzel + Partner-Wurzel, Norm 2 + 2 bzw. 2 - 2 = 0). Identifiziert J die beiden SO(10) "
      f"(Gewicht w + w' für J+, w - w' für J-), so liegen die 16 x 16 Paargewichte in den Klassen {diag_klassen} (0 = "
      "Wurzelklasse, 2 = Vektorklasse): 16 x 16' = 10 + 120 + 126 bzw. 16 x 16q' = 1 + 45 + 210, keine 16. Ohne "
      f"Identifikation ist die Original-16 {kopien}-fach (4 Familien mal 16' des Partners) da",
      diag_klassen == {"J+": [2], "J-": [0]} and kopien == 64 and wurzeln_paar == 104 and einseitig_paar == [52, 52])

# --- 2.3 Kern auf der radialen Leiter: Weyl-Trennung = Skalar/Vektor-Trennung (z83-Standardform, CCR)
AMP_CCR = (lambda m, k: np.sqrt(k + 1), lambda m, k: np.sqrt(m - k))
SCHRITT = (lambda k: k + 1, lambda k: k)


def standardform(M):
    lab = [(m, mb, k, kb) for m in range(1, M + 1) for mb in range(1, M + 1) for k in range(m) for kb in range(mb)]
    ix = {s_: i for i, s_ in enumerate(lab)}
    n_ = len(lab)

    def dicht(r, c, w):
        return sps.csr_matrix((w, (r, c)), shape=(n_, n_)).toarray()

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

    B_ = {(c_, a_): bilinear(c_, a_, 1) for c_ in range(2) for a_ in range(2)}
    return lab, B_, np.diag([0.5*(s_[0] + s_[1]) for s_ in lab])


LAB, BS, ES = standardform(5)
ZH = np.array([s_[0] - s_[1] for s_ in LAB])
Z_SU = expm(1j*np.pi*(ES + 2*(BS[(0, 0)] - BS[(1, 1)])/2))
nichtdiag = float(np.abs(Z_SU - np.diag(np.diag(Z_SU))).max())
zw = {}
for zh in (-2, -1, 0, 1, 2):
    werte = np.diag(Z_SU)[ZH == zh]
    zw[zh] = complex(np.round(werte.mean(), 12)) if np.allclose(werte, werte[0]) else None
viertel = {zh: complex(np.round(np.exp(1j*np.pi*zh/2), 12)) for zh in (-2, -1, 0, 1, 2)}
q_weyl, q_sv = zw[1]/zw[-1], zw[2]/zw[0]
q_weyl_v, q_sv_v = viertel[1]/viertel[-1], viertel[2]/viertel[0]
skalarleiter = int(np.sum(ZH == 0))
alpha_tab = {}
for j in range(8):                                                 # chi(h) = e^(i alpha h), alpha = j pi/2 in [0, 4 pi)
    chi = lambda h: sp.exp(sp.I*sp.pi*sp.Rational(j, 2)*h)
    werte = [sp.nsimplify(sp.simplify(chi(sp.Rational(z, 2)))) for z in range(-3, 4)]
    z4 = all(sp.simplify(w**4 - 1) == 0 for w in werte)
    weyl = sp.simplify(chi(sp.Rational(1, 2)) - chi(-sp.Rational(1, 2))) != 0
    sv = sp.simplify(chi(1) - chi(0)) != 0
    alpha_tab[f"{sp.Rational(j, 2)} pi"] = {"Z4-wertig": z4, "trennt h = +-1/2": weyl, "trennt h = 0, 1": sv}
kern_ok = all((v["trennt h = +-1/2"] == v["trennt h = 0, 1"]) for v in alpha_tab.values() if v["Z4-wertig"])
z4_alphas = [k for k, v in alpha_tab.items() if v["Z4-wertig"]]
ZAHLEN["kern"] = {"Z_SU je 2h (radiale Leiter, M = 5)": {str(k): zs(v) for k, v in zw.items()},
                  "Z_SU(1/2)/Z_SU(-1/2), Z_SU(1)/Z_SU(0)": [zs(q_weyl), zs(q_sv)],
                  "Viertelwinkel e^(i pi h): dieselben Quotienten": [zs(q_weyl_v), zs(q_sv_v)],
                  "Zustände der neutralen Skalarleiter (h = 0)": skalarleiter, "Helizitätscharaktere e^(i alpha h)": alpha_tab}
check("[E] Kern auf der radialen Leiter (Standardform z83, CCR, M = 5): Z_SU = exp(i pi (E + 2 J3_R)) ist diagonal "
      f"(Rest {nichtdiag:.0e}) und konstant je Helizität: {', '.join(f'2h = {k}: {zs(v)}' for k, v in zw.items())}; die "
      f"{skalarleiter} Zustände der neutralen Skalarleiter haben Z_SU = -1. Der Quotient, der die Weyl-Helizitäten trennt, "
      f"Z_SU(1/2)/Z_SU(-1/2) = {zs(q_weyl)}, ist derselbe, der Skalar und Vektor trennt, Z_SU(1)/Z_SU(0) = {zs(q_sv)}; ebenso "
      f"für den Viertelwinkel e^(i pi h) = -Z_SU ({zs(q_weyl_v)}, {zs(q_sv_v)}). Für jeden Helizitätscharakter e^(i alpha h) gilt "
      "chi(1/2)/chi(-1/2) = chi(1)/chi(0) = e^(i alpha); Z4-wertig sind genau alpha in "
      f"{z4_alphas}, und jeder davon trennt die Weyl-Helizitäten genau dann, wenn er Skalar und Vektor trennt",
      nichtdiag < 1e-12 and zw == {-2: 1, -1: 1j, 0: -1, 1: -1j, 2: 1} and q_weyl == -1 and q_sv == -1
      and q_weyl_v == -1 and q_sv_v == -1 and kern_ok and z4_alphas == ["0 pi", "1 pi", "2 pi", "3 pi"]
      and alpha_tab["1 pi"]["trennt h = +-1/2"] and not alpha_tab["2 pi"]["trennt h = +-1/2"])

# --- 2.4 erschöpfend: Nebenklassensatz und die Verträglichkeit von (a), (b), (c), (d)
neben_ok, zaehl, liste_abd, max_chiral = True, {}, [], {}
for s in (1, -1):
    for C in CODES[s]:
        e = eigenschaften(C)
        for r in range(4):
            ds = {x[3] for x in C if x[0] == r}
            if ds:
                d_r = min(ds)
                neben_ok &= ds == {(d_r + d0) % 4 for d0 in e["D0"]}
        key = tuple(k for k in ("chiral", "b1", "b2", "c_R", "gso") if e[k])
        zaehl[(s, key)] = zaehl.get((s, key), 0) + 1
        if e["chiral"] and e["b1"] and e["gso"]:
            liste_abd.append((s, sorted(C), len(C), all(J_plus(x) in C for x in C), all(J_minus(x) in C for x in C)))
        if e["chiral"] and len(C) == 16:
            max_chiral.setdefault(s, []).append((sorted(C), e["b1"], e["b2"], e["c_R"]))
n_ac = sum(v for (s, key), v in zaehl.items() if "chiral" in key and "c_R" in key)
n_ab2 = sum(v for (s, key), v in zaehl.items() if "chiral" in key and "b2" in key)
n_ac_streng = sum(1 for s in (1, -1) for C in CODES[s] if (0, 0, 0, 2) in C)
abd_max = [x for x in liste_abd if x[2] == 16]
abd_J = [x for x in liste_abd if x[3] and x[4]]
zyklisch = all(len(x[1]) == 4 and any(spann([g]) == frozenset(x[1]) and g[0] == 1 and all(c % 2 for c in g[1:])
                                      for g in x[1]) for x in liste_abd)
erz_J = sorted({next(g for g in x[1] if g[0] == 1) for x in abd_J})
q_streng = [str(qform((0, 0, 0, 2), s)) for s in (1, -1)]
huellen = [(s, eigenschaften(C)) for s in (1, -1) for C in voll_paar[s]]
huellen_txt = [f"{'definit' if s == 1 else '(8,8)'}: chiral {e['chiral']}, GSO {e['gso']}, Familie {e['b1']}"
               for s, e in huellen]
ZAHLEN["huellen_von_P"] = huellen_txt
ZAHLEN["erschoepfend"] = {"Kombinationen (Signatur, erfüllte Eigenschaften): Anzahl Codes":
                              {f"{s}, {'+'.join(key) or '-'}": v for (s, key), v in sorted(zaehl.items(), key=str)},
                          "(a) und (b1) und (d)": [{"Signatur": x[0], "Code": str(x[1]), "Ordnung": x[2], "J+": x[3],
                                                    "J-": x[4]} for x in liste_abd],
                          "maximale chirale Codes (definit)": [str(x) for x in max_chiral.get(1, [])],
                          "maximale chirale Codes ((8,8))": [str(x) for x in max_chiral.get(-1, [])]}
check("[negativ, E] Erschöpfend über alle 140 Codes (beide Signaturen): in jedem Code ist die Menge der Helizitätsklassen d "
      "jeder R-Klasse eine Nebenklasse derselben Untergruppe D0 = {d : R = 1} (Nebenklassensatz, alle Codes geprüft). "
      "Chiralität der 16 verlangt D0 = {0}, die neutrale Skalarleiter (R = 1, h = 0, also d = 2) verlangt 2 in D0 (unter dem "
      "Viertelwinkel ist d = 2 die Vektorklasse: dann fehlen statt der Skalare die R-neutralen Vektoren): (a) und (c) "
      f"zusammen in {n_ac} Codes, auch mit beliebiger Familien- und Partnerklasse; die strenge Klasse (0,0,0,2) (nur der "
      f"Vektor 6' des Partners) hat die ungerade Minimalnorm {q_streng} (definit, (8,8); die freien Fermionen der "
      f"Vektorklasse) und liegt in keinem geraden Code ({n_ac_streng}). Die sechs "
      f"selbstdualen Hüllen von P: {huellen_txt}; keine ist chiral und hält die Familie. (a) mit intaktem "
      f"Original-E8 (b2) in {n_ab2} Codes, denn dessen eigene 16 hat d = 0. (a), (b1) und (d) zusammen genau in den "
      f"{len(liste_abd)} zyklischen Codes <(1, +-1, +-1, +-1)> der Ordnung 4 (je Signatur 8, keiner selbstdual); unter J+ und "
      f"J- invariant davon ({len(abd_J)}) nur die mit den Erzeugern {[''.join(map(str, g)) for g in erz_J]}: P = Klassen der "
      "Diagonale (v, v), <(1,1,3,3)> = Klassen der Antidiagonale (v, -v), und beide mit umgekehrter A3-Orientierung. "
      f"Maximal und chiral sind nur {len(max_chiral.get(1, []))} definite Codes, die A3-Tausche (D5 + A3' und D5' + A3_fam, "
      f"beide Orientierungen), alle ohne (b1), (b2) und (c); in (8,8) {len(max_chiral.get(-1, []))}",
      neben_ok and n_ac == 0 and n_ab2 == 0 and n_ac_streng == 0 and q_streng == ["1", "-1"] and len(huellen) == 6
      and all(not (e["chiral"] and e["b1"]) for _, e in huellen) and sum(e["chiral"] for _, e in huellen) == 1
      and len(liste_abd) == 16 and not abd_max and zyklisch
      and sum(1 for x in liste_abd if x[0] == 1) == 8 and len(abd_J) == 8
      and erz_J == [(1, 1, 1, 1), (1, 1, 3, 3), (1, 3, 1, 3), (1, 3, 3, 1)]
      and len(max_chiral.get(1, [])) == 4 and all(not b1 and not b2 and not cR for _, b1, b2, cR in max_chiral[1])
      and len(max_chiral.get(-1, [])) == 0)

# --- 2.5 (e) drei Familien: Markenzyklus auf den vier Kleberzuständen von A3_fam
A4 = sorted(KURZ["A"][1][Fr(3, 4)])
W_D3 = []
for perm in itertools.permutations(range(3)):
    for vz in itertools.product((1, -1), repeat=3):
        if vz.count(-1) % 2 == 0:
            g = np.zeros((3, 3), int)
            for i, p in enumerate(perm):
                g[p, i] = vz[i]
            W_D3.append(g)
zyklen = []
for g in W_D3:
    bild = [A4.index(tuple(int(x) for x in g @ np.array(v))) for v in A4]
    Pm = np.zeros((4, 4))
    for i, j in enumerate(bild):
        Pm[j, i] = 1
    if np.allclose(np.linalg.matrix_power(Pm, 4), np.eye(4)) and not np.allclose(np.linalg.matrix_power(Pm, 2), np.eye(4)):
        zyklen.append((g, Pm))
g4, P4 = zyklen[0]
ew_zust = sorted({complex(np.round(x, 9)) for x in np.linalg.eigvals(P4)}, key=lambda z: (z.real, z.imag))
ew_gitter = sorted({complex(np.round(x, 9)) for x in np.linalg.eigvals(g4.astype(float))}, key=lambda z: (z.real, z.imag))
summe_null = all(x == 0 for x in np.sum(np.array(A4), axis=0))
zentrum_A = sorted({complex(np.round(np.exp(1j*np.pi*sum(v)/2), 9)) for v in A4}, key=str)
ZAHLEN["familien"] = {"Kleberzustände 4_fam (doppelt skaliert)": [str(v) for v in A4],
                      "4-Zyklen in W(D3) = S4 auf den vier Zuständen": len(zyklen),
                      "Eigenwerte auf den vier Zuständen": [zs(z) for z in ew_zust],
                      "Eigenwerte auf dem Familiengitter (Rang 3)": [zs(z) for z in ew_gitter],
                      "Zentrum z_A3 auf den vier Zuständen": [zs(z) for z in zentrum_A],
                      "chirale Kopien der 16 in P": kopien}
check("[negativ, E] (e) drei Familien: die vier Kleberzustände der Familienklasse (q = 3/4) summieren sich zu null "
      f"(Summe-null-Gitter der vier Marken, Rang 3). Ein Markenzyklus (4-Zyklus in W(D3) = S4, {len(zyklen)} Stück) hat auf "
      f"den vier Zuständen die Eigenwerte {[zs(z) for z in ew_zust]}, auf dem Familiengitter genau "
      f"{[zs(z) for z in ew_gitter]}: das sind die drei gelösten Familiencharaktere {{-1, i, -i}}. Das Zentrum, das allein "
      f"der Kleber sieht, ist auf allen vier Zuständen gleich ({[zs(z) for z in zentrum_A]}, Schur): der Kleber gibt allen "
      f"vier dieselbe Helizität, auch der mu4-trivialen Komponente. Chirale Zählung in P: {kopien} Kopien der Original-16 "
      "(4 x 16'), unter der diagonalen SO(10) 0; nie 3",
      summe_null and len(zyklen) == 6 and ew_zust == [-1, -1j, 1j, 1] and ew_gitter == [-1, -1j, 1j]
      and len(zentrum_A) == 1 and np.isclose(zentrum_A[0], -1j) and kopien == 64)

# --- 2.6 Negativkontrollen
C1 = spann([(1, 0, 0, 1), (0, 1, 1, 0)])
CD = spann([(1, 0, 1, 0), (0, 1, 0, 1)])
D16 = spann([(0, 0, 2, 2), (0, 2, 0, 2), (0, 2, 2, 0), (1, 1, 1, 1)])
K2_a1 = sorted({d for d in range(4) if (1 - d) % 2 == 0})         # ein E8, A3 = Familie, Z_SU frei (nur GSO)
kontrollen = {}
for name, C, s in (("K3 direkte Summe E8 + E8', P_A3'", C0, 1), ("K4 A3-Tausch (D5 + A3', D5' + A3_fam)", C1, 1),
                   ("K5 Diagonalcode (a, b, a, b) in (8,8)", CD, -1), ("K6 D16+", D16, 1)):
    kontrollen[name] = dict(eigenschaften(C))
    kontrollen[name]["selbstdual"] = len(C) == 16 and C in MAXI[s]
    kontrollen[name]["GSO-verletzende Klassen"] = sum(1 for x in C if (x[0] - x[3]) % 2)
eK1 = eigenschaften(K1)
kontrollen["K1 ein E8, P_A3 im Original (z97)"] = {**eK1, "Familiengitter": "als konforme Algebra umgedeutet"}
kontrollen["K2 ein E8, A3 = Familie, kein Kleber an Z_SU"] = {"d bei a = 1": K2_a1, "chiral": len(K2_a1) == 1}
tab_C1 = tabelle(C1)
komp_C1 = gitter[(1, C1)]["Komponenten"]
komp_D16 = gitter[(1, D16)]["Komponenten"]
e8A = [w for w in wurzeln(C1, 1) if not any(w[5:13])]              # Wurzeln ohne A3_fam- und D5'-Anteil
ZAHLEN["kontrollen"] = {k: {kk: (str(vv) if not isinstance(vv, (bool, int, str)) else vv) for kk, vv in v.items()}
                        for k, v in kontrollen.items()}
ZAHLEN["kontrollen"]["K4 erlaubte h je R"] = tab_C1
check("[E] Negativkontrollen: K1 ein E8 mit P_A3 im Original (z97) ist chiral und GSO-treu, deutet aber das Familiengitter "
      f"um und verbietet die Skalarleiter (c_R = {eK1['c_R']}); K2 ein E8 mit A3 als Familie ohne Kleber an Z_SU lässt für "
      f"R = i beide Klassen d = {K2_a1} zu (vektorartig, wie z97 Teil 1); K3 die direkte Summe E8 + E8' unter P_A3' ist "
      f"nicht chiral (d bei a = 1: {kontrollen[list(kontrollen)[0]]['d bei a = 1']}) und hat "
      f"{kontrollen[list(kontrollen)[0]]['GSO-verletzende Klassen']} GSO-verletzende Klassen (die eigene 16 des Originals hat "
      "h ganz), erlaubt aber die Skalarleiter als (10', 6'); K4 der A3-Tausch ist selbstdual, chiral und GSO-treu mit "
      f"genau der z97-Tabelle {tab_C1}, aber die Familie klebt an D5' (b1, b2 falsch) und die Skalarleiter ist verboten; "
      f"Wurzelkomponenten {komp_C1}, eine davon D5 + A3' ({len(e8A)} Wurzeln); K5 der Diagonalcode in (8,8) (d bei a = 1: "
      f"{kontrollen['K5 Diagonalcode (a, b, a, b) in (8,8)']['d bei a = 1']}) und K6 D16+ (Wurzelkomponenten {komp_D16}) "
      "sind nicht chiral",
      eK1["chiral"] and eK1["gso"] and not eK1["c_R"] and len(K2_a1) == 2
      and not kontrollen["K3 direkte Summe E8 + E8', P_A3'"]["chiral"]
      and kontrollen["K3 direkte Summe E8 + E8', P_A3'"]["GSO-verletzende Klassen"] == 8
      and kontrollen["K3 direkte Summe E8 + E8', P_A3'"]["c_R"] and kontrollen["K3 direkte Summe E8 + E8', P_A3'"]["b2"]
      and kontrollen["K4 A3-Tausch (D5 + A3', D5' + A3_fam)"]["chiral"]
      and kontrollen["K4 A3-Tausch (D5 + A3', D5' + A3_fam)"]["gso"]
      and not kontrollen["K4 A3-Tausch (D5 + A3', D5' + A3_fam)"]["b1"]
      and not kontrollen["K4 A3-Tausch (D5 + A3', D5' + A3_fam)"]["c_R"] and tab_C1 == z97_tab
      and komp_C1 == [240, 240] and len(e8A) == 240 and komp_D16 == [480]
      and not kontrollen["K5 Diagonalcode (a, b, a, b) in (8,8)"]["chiral"] and not kontrollen["K6 D16+"]["chiral"]
      and all(kontrollen[k]["selbstdual"] for k in list(kontrollen)[:4]))

# --- 2.7 J-Kovarianz der Prämisse P_A3'
lesart = lambda w, seite: ("fam" if seite == 0 else "konf") if any(w[5:8] if seite == 0 else w[13:16]) else "-"
wechsel = 0
for w in W_C0:
    seite = 0 if any(w[0:8]) else 1
    wJ = w[8:16] + w[0:8]
    if lesart(w, seite) != lesart(wJ, 1 - seite):
        wechsel += 1
ZAHLEN["J_kovarianz"] = {"Wurzeln von E8 + E8'": len(W_C0), "Wurzeln, deren Lesart J wechselt": wechsel}
check(f"[E] J-Kovarianz: J (Tausch Original <-> Partner) bildet die {len(W_C0)} Wurzeln von E8 + E8' aufeinander ab; unter "
      f"der asymmetrischen Lesart von H3 (A3 Familie, A3' konform) wechseln {wechsel} davon die Lesart (Familienindex <-> "
      "Twistorindex auf (1,15), (10,6) und (16,4)), nur die 80 SO(10)-Wurzeln nicht. J-kovariant sind nur (Familie, Familie), "
      "ohne Helizität, und (konform, konform), also P_A3 auf beiden Seiten mit dem z97-Konflikt doppelt: P_A3' ist eine "
      "Auswahl, keine Rekonstruktion",
      len(W_C0) == 480 and wechsel == 400)
befund("Teil 2: H3 scheitert an (c) und (e), und zwar strukturell. Jede gruppenwertige Kleberregel macht die Helizitätsklassen "
       "jeder R-Klasse zu einer Nebenklasse von D0; Chiralität verlangt D0 = {0}, die neutrale Skalarleiter 2 in D0. Der "
       "Grund liegt auf der Leiter: der Quotient, der h = +-1/2 trennt, ist derselbe, der h = 0 von h = 1 trennt. Der "
       "gepaarte Code P ist chiral, GSO-treu und hält die Kleberrelation der Familie, ist aber P_A3 auf Zentrumsebene, nicht "
       "selbstdual und ohne Ein-Teilchen-16; er zählt 64 oder 0 Kopien der 16, nie 3", True)

# =================================================================== Teil 3
say("\nTeil 3: Bilanz und Urteil für H3")
STATUS = {
    "Tomita-Kleber als Gitter": ("E", "Klassen der isotropen Diagonale = P = <(1,1,1,1)>, J-invariant, |det| = 16, nicht "
                                      "selbstdual; selbstduale Hüllen: 3 definit, 3 in (8,8)"),
    "Teilchen 16 (x) 4_fam (x) 4'": ("E (negativ)", "Norm 11/4 bzw. 5/4, in keinem ganzzahligen Gitter"),
    "(a) Chiralität": ("E auf Klassen, C mit P_A3'", "P: 16 nur mit Z_SU = -i; aber keine Ein-Teilchen-16, unter der "
                                                     "diagonalen SO(10) keine 16"),
    "(b) Familienklasse": ("E", "Kleberrelation a = b und q = 3/4 bleiben; das Original-E8 als Teilgitter nicht"),
    "(c) neutrale Skalarleiter": ("E (negativ)", "in keinem chiralen Code; Nebenklassensatz und Leiterquotient"),
    "(d) GSO": ("E", "R^2 = (-1)^(2h) auf P"),
    "(e) drei Familien": ("E (negativ)", "Kleber familienblind: 64 oder 0 Kopien, mu4-triviale Komponente mitchiral"),
    "P_A3' (nur Partner)": ("Auswahl", "nicht J-kovariant; auf P gleich P_A3 auf Zentrumsebene"),
    "bedingte Fassung P_A3' + Paarung P": ("nein", "(a), (b1), (d) ja; (c), (e) nein; keine Ein-Teilchen-16"),
    "H3 gesamt": ("nein (gescoptes Negativ)", "Chiralität und neutrale Skalarleiter schließen sich in jeder Z4-Kleberregel aus"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
ZAHLEN["H3_urteil_kurz"] = {"Urteil": "nein, gescoptes Negativ",
                            "T4": "Kleberroute geschlossen: Chiralität aus einer Z4-Kleberregel verbietet R-neutrale "
                                  "Skalare (Zentrum) bzw. R-neutrale Vektoren (Viertelwinkel)",
                            "K16": "chiraler Sektor nur als nicht selbstduales Gitter (|det| = 16); selbstduale Hüllen sind "
                                   "vektorartig, GSO-verletzend oder verlegen die Familie"}
urteil = ("H3 ist nein, als gescoptes Negativ. Der Tomita-Kleber der Hypothese existiert als Gitterobjekt: Er ist genau "
          "die Klassenmenge der isotropen Diagonale von E8 + E8(-1) aus z95, der Code P = <(1,1,1,1)>; die Antidiagonale "
          "(v, -v) gibt <(1,1,3,3)> mit der anderen Helizität. P ist unter J invariant (Tausch mit und ohne Konjugation). "
          "Unter P_A3' (A3' des Partners als konforme Algebra) ist P chiral: "
          "Die 16 hat nur die Helizitätsklasse Z_SU = -i, also h = 1/2 oder -3/2. P ist GSO-treu und hält die Kleberrelation "
          "der Familie (a = b, q(A3) = 3/4). Aber: (1) P ist nicht selbstdual (|det| = 16) und enthält keine Ein-Teilchen-16. "
          "Ein Teilchen 16 (x) 4_fam (x) 4' hat die Norm 11/4 (definit) bzw. 5/4 (Signatur (8,8)) und liegt in keinem "
          "ganzzahligen Gitter. Die 16 kommt nur als Paar Original-Wurzel + Partner-Wurzel vor, 64-fach (4 Familien mal 16'), "
          "und unter der diagonalen SO(10) gar nicht (16 x 16' bzw. 16 x 16q' enthalten keine 16). (2) Auf P gilt "
          "z_A3 = Z_SU', die erlaubten h je R sind Klasse für Klasse die z97-Regel: Der Tomita-Kleber ist P_A3 auf "
          "Zentrumsebene. P_A3' selbst ist nicht J-kovariant, J wechselt auf 400 von 480 Wurzeln Familien- und "
          "Twistorindex. (3) Die neutrale Skalarleiter bleibt verboten, und zwar in jedem Code: Erschöpfend über alle 140 "
          "isotropen Codes beider Signaturen sind die Helizitätsklassen jeder R-Klasse eine Nebenklasse von D0. Chiralität "
          "verlangt D0 = {0}, die Skalarleiter 2 in D0. Der Grund steht auf der radialen Leiter: Z_SU(1/2)/Z_SU(-1/2) = "
          "Z_SU(1)/Z_SU(0) = -1, und für jeden Z4-wertigen Helizitätscharakter trennt er die Weyl-Helizitäten genau dann, "
          "wenn er Skalar und Vektor trennt. Mit dem Viertelwinkel statt des Zentrums verliert man statt der R-neutralen "
          "Skalare die R-neutralen Vektoren (die 45). (4) Mit intaktem Original-E8 ist kein Code chiral, denn dessen eigene "
          "16 hat d = 0; die Kommutante sieht die Helizität einseitiger Originalmaterie nicht. Selbstdual und chiral sind nur "
          "die vier A3-Tausche (D5 + A3', D5' + A3_fam). Sie geben genau die z97-Tabelle und legen die Familie an das SO(10)' "
          "des Partners. (5) Der Kleber ist familienblind: Das Zentrum ist auf den vier Kleberzuständen gleich, der "
          "Markenzyklus hat dort die Eigenwerte 1, i, -1, -i und auf dem Familiengitter die gelösten {-1, i, -i}. Die "
          "mu4-triviale Komponente wird mitchiral, drei Familien folgen nicht. Folge für T4: Die Kleberroute zur Chiralität "
          "ist geschlossen, unabhängig davon, welches E8 die konforme A3 trägt und ob Tomita paart. Chiralität aus einer "
          "Z4-Kleberregel und R-neutrale Skalare und Vektoren zugleich gibt es nicht; T4 bleibt offen, nur über einen "
          "nichtzentralen Mechanismus (Index oder Fluss, gelöster Familienindex). Folge für K16: Der chirale Sektor ist nur "
          "ein nicht selbstduales Gitter; seine selbstdualen Hüllen sind vektorartig, GSO-verletzend oder verlegen die "
          "Familie. Das ist die Gitterform des Spiegelproblems: Chiralität braucht einen relativen Sektor mit Zufluss. "
          "K16 bleibt L. Nichts ist befördert; 3+1-Aussagen bleiben Kinematik der Einpaarschicht.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z102_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
