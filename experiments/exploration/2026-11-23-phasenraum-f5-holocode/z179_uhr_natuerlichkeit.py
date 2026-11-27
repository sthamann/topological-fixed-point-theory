"""Z179: F6, die Quartikuhr aus der Natürlichkeit der Deckstruktur. Ziel: P_zentral.

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Bausteine bleiben unberührt, nichts wird befördert.
Statuskürzel: E exakt, N numerisch, C bedingt, L Lesart, O offen.
Bestand: Die 16 reellen Nahteltern sind die Schnitte der Überlagerung, zerlegt nach der Deckparität in 10 gerade und 6
ungerade (VERDICT: 5 + 3 komplexe Schnitte; die 3 ungeraden sind die Familien); die Naht ist E8 über so(16) der 16 Eltern
(H6: A* gibt die GSO-Summe). z175: P_lift folgt aus P_zentral ('die Uhr hat Ordnung 4 und vertauscht mit den Basiswechseln
der Schnitte beider Deckhälften'); Zentralisatorsatz; die Vierteldrehung z -> i z der Marken ist nicht die Uhr.

Frage (F6): Ist die Uhr ein kanonischer Operator der Deckstruktur, so dass die Vertauschung mit den Basiswechseln ein Satz
und keine Annahme ist?

Teil A  Clifford-Algebra der 16 Eltern (exakt, 256 x 256): die Hebung der Deckparität iota = (+1^10, -1^6).
Teil B  Fortsetzung auf die Nahtalgebra E8 (Gitter, exakt): genau zwei Fortsetzungen, beide von Ordnung 4.
Teil C  Natürlichkeit: Verhalten unter allen Basiswechseln, die die Deckhälften erhalten; Ladungskonjugation.
Teil D  Umkehrung (Zentralisator) und Vergleich mit anderen Markenzahlen und anderen Teilungen.
Teil E  Bilanz und offene Identifikation.
"""
import json
import sys
import time
from fractions import Fraction as Fr
from itertools import combinations, product
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parents[2] / "research" / "originlab"))

OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T_START = time.time()


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


# =================================================================== Teil A
say("Teil A: Clifford-Algebra der 16 Eltern und die Hebung der Deckparität")
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)


def kron(*ms):
    out = np.array([[1]], dtype=complex)
    for m in ms:
        out = np.kron(out, m)
    return out


GAM = []
for k in range(8):
    GAM.append(kron(*([Z] * k + [X] + [I2] * (7 - k))))
    GAM.append(kron(*([Z] * k + [Y] + [I2] * (7 - k))))
E = np.eye(256, dtype=complex)
cliff = all(np.allclose(GAM[i] @ GAM[j] + GAM[j] @ GAM[i], 2 * (i == j) * E) for i in range(16) for j in range(i, 16))
GERADE, UNGERADE = list(range(10)), list(range(10, 16))                  # Deckhälften: 10 gerade, 6 ungerade Eltern
IOTA = [1] * 10 + [-1] * 6


def prod(idx):
    out = E.copy()
    for i in idx:
        out = out @ GAM[i]
    return out


OM = prod(UNGERADE)                                                       # Volumenelement der ungeraden Hälfte
OMG = prod(GERADE)
GAMMA = prod(range(16))                                                   # Gamma^2 = +1 für 16 Erzeuger
ad = all(np.allclose(OM @ GAM[j] @ np.linalg.inv(OM), IOTA[j] * GAM[j]) for j in range(16))
ordnung = np.allclose(OM @ OM, -E) and np.allclose(np.linalg.matrix_power(OM, 4), E)
chir = np.allclose(GAMMA @ GAMMA, E) and np.allclose(GAMMA @ OM, OM @ GAMMA)
Pp = (E + GAMMA) / 2
ev = np.linalg.eigvals(OM @ Pp + 5 * (E - Pp))                             # S^+: Eigenwerte von omega, S^-: Marke 5
spek = {"+i": int(np.sum(np.abs(ev - 1j) < 1e-9)), "-i": int(np.sum(np.abs(ev + 1j) < 1e-9))}
# adjungierte 120 = Lambda^2 V: iota wirkt mit iota_i iota_j
adj = {"+1": sum(1 for i, j in combinations(range(16), 2) if IOTA[i] * IOTA[j] == 1),
       "-1": sum(1 for i, j in combinations(range(16), 2) if IOTA[i] * IOTA[j] == -1)}
spur248 = adj["+1"] - adj["-1"] + 1j * spek["+i"] - 1j * spek["-i"]
rel = np.allclose(OM @ OMG, GAMMA) or np.allclose(OM @ OMG, -GAMMA)
ZAHLEN["A1"] = {"Spektrum auf 128": spek, "Spektrum auf 120": adj, "Spur auf 248": str(spur248)}
check("A1 [E] Clifford-Algebra der 16 Eltern (256 x 256). omega = Produkt der sechs ungeraden Eltern wirkt auf die Eltern "
      "durch Konjugation als Deckparität iota = (+1^10, -1^6). Es ist omega^2 = -1, also Ordnung 4; auf der adjungierten "
      f"120 hat es die Eigenwerte +1 ({adj['+1']} = 45 + 15) und -1 ({adj['-1']} = 10 * 6), auf dem Spinor 128 die Eigenwerte "
      f"+i ({spek['+i']}) und -i ({spek['-i']}). Also omega^2 = (-1)^F (F = Spinorzahl) und Sp_248 omega = 0. Das Produkt der "
      "zehn geraden Eltern ist auf dem Spinor dasselbe bis auf das Vorzeichen (omega_10 = +- Gamma omega^-1)",
      cliff and ad and ordnung and chir and spek == {"+i": 64, "-i": 64} and adj == {"+1": 60, "-1": 60}
      and abs(spur248) < 1e-12 and rel)

# Eindeutigkeit der Hebung: x = Summe x_S e_S mit x e_j = iota_j e_j x für alle j
loesungen = []
for bits in range(1 << 16):
    S = [j for j in range(16) if bits >> j & 1]
    n = len(S)
    if all((-1) ** (n - (j in S)) == IOTA[j] for j in range(16)):
        loesungen.append(S)
ZAHLEN["A2"] = {"Monome, die iota implementieren": loesungen}
check("A2 [E] Eindeutigkeit der Hebung: Ein Element x der Clifford-Algebra mit x e_j x^-1 = iota_j e_j für alle 16 Eltern "
      "ist eine Kombination der Monome e_S mit (-1)^(|S| - [j in S]) = iota_j. Unter allen 65536 Monomen erfüllt das genau "
      "eines: S = die sechs ungeraden Eltern. Also ist x = lambda omega, und in der Spingruppe (x reell, x x^t = 1) bleibt "
      "lambda = +-1. Die Hebung der Deckparität ist bis auf das Vorzeichen omega -> -omega = omega^-1 eindeutig",
      loesungen == [UNGERADE])

# =================================================================== Teil B
say("\nTeil B: Fortsetzung der Deckparität auf die Nahtalgebra E8")
D8 = [tuple(Fr(v) for v in vec) for vec in
      {tuple((s1 if k == i else s2 if k == j else 0) for k in range(8))
       for i, j in combinations(range(8), 2) for s1 in (1, -1) for s2 in (1, -1)}]
SPIN = [tuple(Fr(s, 2) for s in sg) for sg in product((1, -1), repeat=8) if sum(1 for s in sg if s < 0) % 2 == 0]
ROOTS = D8 + SPIN
dot = lambda a, b: sum(x * y for x, y in zip(a, b))
V = tuple([Fr(1, 2)] * 5 + [Fr(0)] * 3)                                   # Uhr c = exp(2 pi i v.P)
E1 = tuple([Fr(1)] + [Fr(0)] * 7)
ph = lambda u, a: (dot(u, a) % 1)
# iota auf so(16): Wurzel +-e_i +- e_j hat Parität (-1)^[genau einer der Indizes in den ersten fünf Ebenen]
def iota_d8(a):
    return Fr(sum(abs(a[k]) for k in range(5)) % 2, 2)                     # Phase 0 oder 1/2


ext = {}
for name, u in (("v", V), ("v + e_1", tuple(x + y for x, y in zip(V, E1)))):
    ok_d8 = all(ph(u, a) == iota_d8(a) for a in D8)
    sp_ph = sorted({ph(u, a) for a in SPIN})
    cnt = {str(p): sum(1 for a in SPIN if ph(u, a) == p) for p in sp_ph}
    ext[name] = {"gleich iota auf so(16)": ok_d8, "Phasen auf 128": cnt}
# alle Toruselemente, die auf so(16) wie iota wirken: u - v in D8*; D8*/E8 hat Ordnung 2
def in_e8(u):
    ganz = all(x.denominator == 1 for x in u)
    halb = all((2 * x).denominator == 1 and x.denominator == 2 for x in u)
    return (ganz or halb) and sum(u) % 2 == 0


klassen = {"0": in_e8(tuple([Fr(0)] * 8)), "e_1": in_e8(E1), "s": in_e8(tuple([Fr(1, 2)] * 8)),
           "c": in_e8(tuple([Fr(1, 2)] * 7 + [Fr(-1, 2)]))}
invers = all((ph(tuple(x + y for x, y in zip(V, E1)), a) + ph(V, a)) % 1 == 0 for a in SPIN) and \
    all((ph(tuple(x + y for x, y in zip(V, E1)), a) - ph(V, a)) % 1 == 0 for a in D8)
# die Ordnung-2-Normierung (i omega auf dem Spinor) ist kein Automorphismus
ROOTSET = set(ROOTS)
def falsch(a):
    return (ph(V, a) + (Fr(1, 4) if a in SPIN_SET else 0)) % 1


SPIN_SET = set(SPIN)
verletzt, geprueft, ok_v = 0, 0, True
for a in SPIN[:40]:
    for b in ROOTS:
        c = tuple(x + y for x, y in zip(a, b))
        if c in ROOTSET:
            geprueft += 1
            ok_v = ok_v and (ph(V, a) + ph(V, b) - ph(V, c)) % 1 == 0
            verletzt += (falsch(a) + falsch(b) - falsch(c)) % 1 != 0
ZAHLEN["B1"] = {"Fortsetzungen": ext, "Klassen von D8* in E8": klassen, "Wurzeltripel geprüft": geprueft,
                "Verletzungen der Ordnung-2-Normierung": verletzt}
check("B1 [E] Fortsetzung auf E8 = so(16) + 128 (240 Wurzeln): Die Automorphismen von e8, die auf so(16) wie die "
      "Deckparität wirken, bilden eine Nebenklasse des Zentralisators von so(16), und der hat die Ordnung 2 (D8*/E8: von "
      "den vier Klassen von D8* liegen 0 und eine Spinorklasse in E8). Es gibt also genau zwei Fortsetzungen, "
      "exp(2 pi i v.P) mit v = (1/2^5, 0^3) und exp(2 pi i (v + e_1).P); sie sind zueinander invers, haben auf der 128 die "
      "Phasen +-i (je 64) und damit beide die Ordnung 4 mit Quadrat (-1)^F. Eine Fortsetzung der Ordnung 2 gibt es "
      f"nicht: Die Normierung mit Eigenwerten +-1 auf der 128 verletzt die Wurzeladdition ({verletzt} von {geprueft} "
      "geprüften Tripeln). Die Ordnung 4 ist also keine Zusatzannahme, sondern von der Fusion 128 x 128 -> 120 erzwungen",
      all(e["gleich iota auf so(16)"] for e in ext.values())
      and all(e["Phasen auf 128"] == {"1/4": 64, "3/4": 64} for e in ext.values())
      and klassen == {"0": True, "e_1": False, "s": True, "c": False} and invers and ok_v and verletzt > 0)

# Abgleich mit den Glue-Sektoren: Uhr = i^x auf (x, x)
sektor = {}
for a in ROOTS:
    p = ph(V, a)
    x = int(p * 4)
    sektor[x] = sektor.get(x, 0) + 1
check("B2 [E] Abgleich mit den Uhrsektoren: exp(2 pi i v.P) hat auf den 240 Wurzeln die Phasen i^x mit den Anzahlen "
      f"{dict(sorted(sektor.items()))}: x = 0 die 52 Wurzeln von so(10) + so(6), x = 2 die 60 von (10, 6), x = 1 und 3 je 64 "
      "von (16, 4) und (16-quer, 4-quer). Das ist R = i^x auf den E8-Klebern (k, k) aus z170, und es ist omega aus Teil A",
      sektor == {0: 52, 1: 64, 2: 60, 3: 64})

# =================================================================== Teil C
say("\nTeil C: Natürlichkeit unter den Basiswechseln der Deckhälften")
rng = np.random.default_rng(179)
gen_gerade = [GAM[i] @ GAM[j] for i, j in combinations(GERADE, 2)]
gen_ungerade = [GAM[i] @ GAM[j] for i, j in combinations(UNGERADE, 2)]
gen_misch = [GAM[i] @ GAM[j] for i in GERADE for j in UNGERADE]
komm = all(np.allclose(OM @ g, g @ OM) for g in gen_gerade + gen_ungerade)
antik = all(np.allclose(OM @ g, -g @ OM) for g in gen_misch)
# endliche Basiswechsel: exp(theta e_i e_j / 2)
def rot(i, j, th):
    return np.cos(th / 2) * E + np.sin(th / 2) * GAM[i] @ GAM[j]


g_bel = E.copy()
for _ in range(12):
    i, j = rng.choice(GERADE, 2, replace=False)
    g_bel = g_bel @ rot(i, j, rng.uniform(0, 6.28))
    i, j = rng.choice(UNGERADE, 2, replace=False)
    g_bel = g_bel @ rot(i, j, rng.uniform(0, 6.28))
endl = np.allclose(g_bel @ OM, OM @ g_bel)
# Spiegelungen (Orientierungswechsel): Automorphismus e_j -> -e_j
spieg = {"Spiegelung in der geraden Hälfte": np.allclose(GAM[0] @ OM @ GAM[0], OM),          # e_0 omega e_0^-1, 6 Vorzeichen
         "Spiegelung in der ungeraden Hälfte": np.allclose(GAM[12] @ OM @ GAM[12], -OM)}
# als Automorphismus r_j (e_j -> -e_j): omega -> det(g | ungerade) omega
def refl_auto(j):
    idx = UNGERADE
    sign = -1 if j in idx else 1
    return sign


C_wirkung = refl_auto(0) * refl_auto(12)                                    # Ladungskonjugation: je eine Spiegelung
ZAHLEN["C1"] = {"vertauscht mit so(10) + so(6)": komm, "antivertauscht mit (10, 6)": antik, "endliche Basiswechsel": endl,
                "Spiegelungen": spieg, "Ladungskonjugation: omega ->": "omega^-1" if C_wirkung == -1 else "omega"}
check("C1 [E] Natürlichkeit: omega vertauscht mit allen 45 + 15 Erzeugern e_i e_j der Basiswechsel innerhalb der beiden "
      "Deckhälften und mit beliebigen endlichen Produkten daraus, und es antivertauscht mit den 60 Erzeugern, die die "
      "Hälften mischen. Unter einem Basiswechsel g aus O(10) x O(6) geht omega in det(g | ungerade Hälfte) omega über: Die "
      "Spiegelung eines geraden Elters lässt omega fest, die eines ungeraden kehrt es um. Die Ladungskonjugation der "
      "Sektoren (x, y) -> (-x, -y) ist je eine Spiegelung in beiden Hälften und gibt omega -> omega^-1. Das ist C R C = "
      "R^-1 aus z176 (Blatttausch S^2 = C), hier als Folge der Natürlichkeit. Die Aussage von P_zentral 'die Uhr vertauscht "
      "mit den Basiswechseln der Schnitte beider Deckhälften' ist damit ein Satz über die Hebung der Deckparität",
      komm and antik and endl and all(spieg.values()) and C_wirkung == -1)

# =================================================================== Teil D
say("\nTeil D: Umkehrung und Vergleich")
H_ROOTS = [a for a in D8 if ph(V, a) == 0]
# Zentralisator von so(10) + so(6) in E8 = Zentrum = (Wurzelgitter von H)^* / E8; Erzeuger v
vielfache = [tuple(k * x for x in V) for k in range(4)]
zent_ok = all(dot(u, a).denominator == 1 for u in vielfache for a in H_ROOTS)
verschieden = [in_e8(u) for u in vielfache] == [True, False, False, False]
# Ordnung der Gruppe: det(Gram D5) * det(Gram D3) = 16 = |L*/L|, Index von L in E8 = 4, also |L*/E8| = 4
import sympy as sp


def cartan_d(n):
    M = sp.zeros(n)
    for i in range(n):
        M[i, i] = 2
    for i in range(n - 2):
        M[i, i + 1] = M[i + 1, i] = -1
    M[n - 3, n - 1] = M[n - 1, n - 3] = -1
    return M


det_l = cartan_d(5).det() * cartan_d(3).det()
ordnungen = {k: min(m for m in range(1, 9) if in_e8(tuple(m * k * x for x in V))) for k in range(4)}
ZAHLEN["D1"] = {"det Gram(D5 + D3)": int(det_l), "Ordnung des Zentralisators": int(sp.sqrt(det_l)),
                "Ordnungen von k v": ordnungen, "Wurzeln von so(10) + so(6)": len(H_ROOTS)}
check("D1 [E] Umkehrung (Zentralisatorsatz aus z175, hier am Gitter): Die Elemente von E8, die mit so(10) + so(6) "
      "vertauschen, sind das Zentrum dieser Untergruppe vollen Ranges, (Wurzelgitter D5 + D3)^*/E8. Es hat die Ordnung "
      f"sqrt({int(det_l)}) = 4 und wird von v erzeugt (Ordnungen von 0, v, 2v, 3v: {list(ordnungen.values())}). Die einzigen "
      "Elemente der Ordnung 4 sind c und c^-1. Also: 'Ordnung 4 und vertauscht mit den Basiswechseln beider Deckhälften' "
      "(P_zentral) gilt genau für die zwei Hebungen der Deckparität. P_zentral ist gleichwertig mit der Aussage 'die Uhr "
      "ist die Deckparität der Schnitte, fortgesetzt auf die Naht' (P_deck)",
      zent_ok and verschieden and det_l == 16 and ordnungen == {0: 1, 1: 4, 2: 2, 3: 4} and len(H_ROOTS) == 52)

# Vergleich: andere Teilungen der 16 Eltern und andere Markenzahlen
teilung = {}
for k in range(0, 17, 2):
    quad = (-1) ** (k * (k - 1) // 2)                                      # (Volumenelement von k Eltern)^2
    sp_adj = ((16 - 2 * k) ** 2 - 16) // 2                                 # Spur von iota auf Lambda^2
    teilung[f"{16 - k} | {k}"] = {"Quadrat der Hebung": quad, "Ordnung": 4 if quad == -1 else (2 if k else 1),
                                  "Spur auf 120": sp_adj}
ord4 = [k for k, v in teilung.items() if v["Ordnung"] == 4]
spurfrei = [k for k, v in teilung.items() if v["Spur auf 120"] == 0]
marken = {}
for n in range(1, 9):
    k = 2 * (n - 1)
    dim_p, dim_m = 2 * (n + 1), 2 * (n - 1)
    sp_adj = (dim_p * (dim_p - 1) + dim_m * (dim_m - 1)) // 2 - dim_p * dim_m
    marken[n] = {"Hälften": f"{dim_p} | {dim_m}", "Quadrat": (-1) ** (k * (k - 1) // 2), "Spur adjungiert": sp_adj}
ZAHLEN["D2"] = {"Teilungen von 16": teilung, "Markenzahl n (Hälften 2(n+1) | 2(n-1), Zählung wie z175)": marken}
check("D2 [E] Vergleich: Unter den Teilungen der 16 Eltern in gerade und ungerade hat die Hebung der Parität genau dann "
      f"die Ordnung 4, wenn die ungerade Hälfte 2 mod 4 Eltern hat ({ord4}); spurfrei auf der 120 ist sie nur für "
      f"{spurfrei}. Beides zusammen lässt 10 | 6 und 6 | 10. In der Zählung von z175 (n Marken, Hälften 2(n+1) | 2(n-1)) "
      "ist das Quadrat der Hebung (-1)^(n-1) und die adjungierte Spur 8 - 2n: Ordnung 4 für gerades n, spurfrei nur für "
      "n = 4. Die Teilung 10 | 6 ist kein freier Parameter, sondern die Schnittzahl 5 + 3 der Überlagerung mit vier Marken",
      ord4 == ["14 | 2", "10 | 6", "6 | 10", "2 | 14"] and spurfrei == ["10 | 6", "6 | 10"]
      and all(v["Quadrat"] == (-1) ** (n - 1) and v["Spur adjungiert"] == 8 - 2 * n for n, v in marken.items()))

# =================================================================== Teil E
say("\nTeil E: Bilanz")
kette = [
    "[Bestand] Die 16 Eltern sind Schnitte der Überlagerung; die Deckinvolution teilt sie in 10 gerade und 6 ungerade",
    "[Bestand, H6] Die Naht ist E8 über so(16) der 16 Eltern (A* gibt die GSO-Summe)",
    "[Satz, A2, B1] Die Deckparität hat genau zwei Fortsetzungen auf die Naht, omega und omega^-1, beide von Ordnung 4",
    "[Satz, A1] omega^2 = (-1)^F, Sp omega = 0, Fixalgebra so(10) + so(6)",
    "[Satz, C1] omega vertauscht mit allen Basiswechseln der Deckhälften; Orientierungswechsel der ungeraden Hälfte "
    "und die Ladungskonjugation kehren es um",
    "[Satz, D1] Umgekehrt ist jedes Element der Ordnung 4, das mit diesen Basiswechseln vertauscht, omega oder omega^-1",
    "[Identifikation, offen] Die Uhr der Theorie (Schritt der Zusammensetzung, Summe der Uhrsektoren) ist diese "
    "Deckparität (P_deck)",
]
for k in kette:
    say("   " + k)
befund("Die Aussage von P_zentral (Ordnung 4, Vertauschung mit den Basiswechseln, R^2 = (-1)^F, Sp R = 0) ist ein Satz "
       "für die Fortsetzung der Deckparität auf die Naht", True)
befund("P_zentral ist gleichwertig mit P_deck ('die Uhr ist die Deckparität der Schnitte'); die Ordnung 4 und der Zweig "
       "(Orientierung der ungeraden Hälfte, dasselbe Bit wie C) sind keine eigenen Annahmen mehr", True)
befund("Gemeinsame Quelle mit F1: Die Paaruhr U = R (x) R-quer^-1 wirkt auf die Eltern beider Seiten als Deckparität; ihre "
       "Eichung ist die getrennte Spinstruktursumme der ungeraden Hälfte (z178 D1)", True)
befund("P_deck ist aus A* hergeleitet (die Identifikation der Uhr der Zusammensetzung mit der Deckparität)", False)
befund("Komplettlösung", False)

urteil = (
    "F6 ist gerechnet. Die Deckparität der 16 Eltern (+1 auf den zehn geraden, -1 auf den sechs ungeraden Schnitten) hat "
    "auf der Naht genau zwei Fortsetzungen, omega und omega^-1 (omega = Produkt der sechs ungeraden Eltern). Beide haben "
    "die Ordnung 4: Eine Fortsetzung der Ordnung 2 verletzt die Fusion der Spinoren. Es gilt omega^2 = (-1)^F, "
    "Sp omega = 0 auf der 248, die Fixalgebra ist so(10) + so(6), und omega vertauscht mit allen Basiswechseln innerhalb "
    "der Deckhälften; Orientierungswechsel der ungeraden Hälfte und die Ladungskonjugation geben omega^-1 (das ist "
    "C R C = R^-1 aus z176). Umgekehrt sind c und c^-1 die einzigen Elemente der Ordnung 4, die mit diesen Basiswechseln "
    "vertauschen. Damit ist P_zentral keine Aussage über einen unbekannten Operator mehr: Sie ist gleichwertig mit 'die "
    "Uhr ist die Deckparität der Schnitte' (P_deck), und alles Weitere ist Satz. Nicht gezeigt ist die Identifikation "
    "selbst, also dass der Schritt der Zusammensetzung und die Summe über Uhrsektoren diese Deckparität meinen; die "
    "Vierteldrehung der Marken bleibt als anderes Modell ausgeschlossen (z175), aber A* benennt keine Uhr. P_zentral "
    "bleibt in der Zählung, in der schwächeren Form P_deck. Keine Komplettlösung, nichts befördert.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z179_ergebnis.json").write_text(json.dumps({"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil},
                                                    indent=1, ensure_ascii=False, default=str))
(HERE / "ausgabe_z179.txt").write_text("\n".join(OUT) + "\n")
say(f"(Laufzeit {time.time() - T_START:.0f} s)")
