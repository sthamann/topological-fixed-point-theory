"""Z183: Angriff auf P_std. Ist die Deckparität der Naht eine Symmetrie oder eine Redundanz?

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Bausteine bleiben unberührt, nichts wird befördert.
Statuskürzel: E exakt, N numerisch, C bedingt, L Lesart, O offen.
Bestand: z170 (27 Spiegelvervollständigungen der 16 Uhrsektoren, A = Z4 x Z4, q = (5 x^2 + 3 y^2)/8; Antwortpotenz 5/4
mit der 16 als nativem Träger), z178 (K19 folgt aus P_std, der GNS-Regel und A*; Gegenmodelle E8 (x) E8-quer und eine
Spin-Form; die Eichsumme ist die getrennte Spinstruktursumme der Deckhälften), z179 (die Uhr ist die Deckparität der
Schnitte, fortgesetzt auf die Naht: P_zentral <=> P_deck), z180 (Nahttorus = Doppelüberlagerung der vier Marken).

Frage. K19 und K20 ruhen auf P_std, und zwar auf seinem ersten Teil: Die Paaruhr ist superselektiert. Lässt sich das
aus A*, aus der Zeitneutralität oder aus der Überlagerungsgeometrie des Bestands herleiten, oder gibt es Modelle, die
alles erfüllen und es verletzen?

Teil A  Exakt auf der Sektorebene: Deckparität als Monodromie der Eltern, Markenfelder, ein Lemma (Markenfeld <=> kein
        einseitiger Strom mischt die Deckhälften) und der Satz 'A* und ein Markenfeld geben die Eichsumme'.
Teil B  Die Gegenmodelle unter A*: Klassifikation, Zustandssummen, und welche Strukturbedingung welchen Typ ausschließt
        (Orientierungswechsel der Schnittbasen, Zeitumkehr, Mischströme).
Teil C  Liefert der Bestand die Aussage? A allein (Bahnen der relativen Spinstruktur), die Überlagerungsgeometrie
        (Windungszahlen der Zyklen um die Marken), der gelöste alpha-Block (Träger der Antwortpotenz 5/4).
Teil D  Bilanz: was aus P_std wird, und der Abstand zur Komplettlösung.
"""
import json
import sys
import time
from fractions import Fraction as Fr
from itertools import product
from pathlib import Path

import mpmath as mp
import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parents[2] / "research" / "originlab"))
import clock_sector_search as cs  # noqa: E402  (nur lesend)

OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T_START = time.time()
mp.mp.dps = 40


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
say("Teil A: Deckparität als Monodromie, Markenfelder und die Eichsumme (exakt)")
A = cs.A
COMP = [frozenset(s) for s in cs.isotropic_subgroups() if len(s) == 16]
ZERO = (0, 0)
M_L, M_R = ((2, 2), ZERO), (ZERO, (2, 2))          # der Mischstrom (10, 6), links und rechts
H2 = Fr(1, 2)


def vertreter(n, k):
    """Vertreter der Klasse k von D_n^*/D_n = Z4 (n ungerade): 0, Spinor, Vektor, konjugierter Spinor."""
    if k == 0:
        return [Fr(0)] * n
    if k == 2:
        return [Fr(1)] + [Fr(0)] * (n - 1)
    v = [H2] * n
    if k == 3:
        v[-1] = -H2
    return v


def klasse_block(v):
    """Klasse eines Vektors von D_n^* in Z4."""
    if all(c.denominator == 1 for c in v):
        return 0 if sum(v) % 2 == 0 else 2
    d = [c - H2 for c in v]
    return 1 if sum(d) % 2 == 0 else 3


def klasse(lam):
    return klasse_block(lam[:5]), klasse_block(lam[5:])


def e8_wurzeln():
    out = []
    for i in range(8):
        for j in range(i + 1, 8):
            for si, sj in product((1, -1), repeat=2):
                v = [Fr(0)] * 8
                v[i], v[j] = Fr(si), Fr(sj)
                out.append(tuple(v))
    for signs in product((1, -1), repeat=8):
        if sum(1 for s in signs if s < 0) % 2 == 0:
            out.append(tuple(Fr(s, 2) for s in signs))
    return out


WURZELN = e8_wurzeln()
WSET = set(WURZELN)
rng = np.random.default_rng(183)
mono_ok, gew_ok = True, True
mono_tab = {}
for (x, y) in A:
    lam = vertreter(5, x) + vertreter(3, y)
    for _ in range(20):                                        # derselbe Sektor, verschoben um Gittervektoren
        sh = [int(c) for c in rng.integers(-3, 4, size=8)]
        if sum(sh[:5]) % 2:
            sh[0] += 1
        if sum(sh[5:]) % 2:
            sh[5] += 1
        mu = [a + b for a, b in zip(lam, sh)]
        m = [1 if c.denominator == 1 else -1 for c in mu]      # exp(2 pi i mu.e_j) mit dem Elter e_j
        mono_ok &= m[:5] == [(-1) ** x] * 5 and m[5:] == [(-1) ** y] * 3 and klasse(mu) == (x, y)
        gew_ok &= (sum(c * c for c in mu) / 2) % 1 == cs.q_form((x, y))
    mono_tab[(x, y)] = ((-1) ** x, (-1) ** y)
par = lambda a: (a[0] + a[1]) % 2                              # 1: die Eltern sehen die Deckparität
ungerade = [a for a in A if par(a)]
gew_ungerade = sorted({str(cs.q_form(a)) for a in ungerade})
par_ist_achtel = all(par(a) == (8 * cs.q_form(a)) % 2 for a in A)
par_ist_misch = all(Fr(par(a), 2) == cs.b_form(a, (2, 2)) for a in A)
fluss_ok = all(par(a) == ((a[0] - s * a[1]) % 2) for a in A for s in (1, 3))
kl = {}
for r in WURZELN:
    kl[klasse(r)] = kl.get(klasse(r), 0) + 1
wurzel_mono = {tuple(1 if c.denominator == 1 else -1 for c in r) for r in WURZELN}
ZAHLEN["A1_klassen_der_E8_wurzeln"] = {str(k): v for k, v in sorted(kl.items())}
ZAHLEN["A1_gewichte_mit_deckparitaet"] = gew_ungerade
check("A1 [E] Eltern und Sektoren im Gitter D5 + D3: Die 16 Eltern sind die Vektorgewichte +-e_j (zehn in D5, sechs in D3). "
      "Ein Feld im Sektor (x, y) hat mit ihnen die Monodromie (-1)^x auf den zehn und (-1)^y auf den sechs (16 Klassen, je 20 "
      "Vertreter; h = (5 x^2 + 3 y^2)/8 mod 1). Die Eltern sehen also genau dann die Deckparität (+1 auf einer Hälfte, -1 auf "
      "der anderen), wenn x + y ungerade ist. Das sind acht Sektoren, die mit ungeradem Uhrfluss für beide Kleber, mit den "
      f"Gewichten {gew_ungerade} mod 1: genau die Sektoren mit h in ungeraden Achteln, und genau die, die mit dem Mischstrom "
      "(10, 6) = (2, 2) die Monodromie -1 haben. Die 240 Wurzeln von E8 liegen in den Klassen "
      f"{ZAHLEN['A1_klassen_der_E8_wurzeln']} (alle gerade): Jedes Feld von E8 hat mit allen 16 Eltern dieselbe Monodromie "
      "(+1 für so(16), -1 für die 128, also (-1)^F), nie die Deckparität",
      mono_ok and gew_ok and len(ungerade) == 8 and gew_ungerade == ["1/8", "3/8", "5/8", "7/8"] and par_ist_achtel
      and par_ist_misch and fluss_ok and kl == {(0, 0): 52, (2, 2): 60, (1, 1): 64, (3, 3): 64}
      and wurzel_mono == {(1,) * 8, (-1,) * 8})

marke = lambda L: any(par(g[0]) and par(g[1]) for g in L)
einseitig = lambda L: any(par(g[0]) != par(g[1]) for g in L)
ohne_misch = lambda L: M_L not in L and M_R not in L
t_inv = lambda L: all(cs.Q_lr(g) == 0 for g in L)
links = lambda L: frozenset(g[0] for g in L if g[1] == ZERO)
rechts = lambda L: frozenset(g[1] for g in L if g[0] == ZERO)
graph = lambda L: len(links(L)) == 1 and len(rechts(L)) == 1
spektrum = lambda L: tuple(sorted((str(cs.q_form(g[0])), str(cs.q_form(g[1]))) for g in L))
lemma = all(marke(L) == ohne_misch(L) for L in COMP)
beidseitig = all((M_L in L) == (M_R in L) for L in COMP)
nie_einseitig = not any(einseitig(L) for L in COMP)
dual = all(all(cs.B_lr(g, M_L) == 0 for g in L) == (M_L in L) for L in COMP)        # L = L^perp
ZAHLEN["A2_mit_markenfeld (alle, mit T)"] = [sum(marke(L) for L in COMP), sum(marke(L) and t_inv(L) for L in COMP)]
check("A2 [E] Lemma (alle 27 Vervollständigungen): Eine Vervollständigung enthält genau dann ein Markenfeld, also ein "
      "lokales Feld, um das die Eltern beider Seiten die Deckparität als Monodromie haben, wenn der Mischstrom (10, 6) "
      "weder links noch rechts zur chiralen Algebra gehört. Grund: Die Parität ist die Monodromie mit dem Mischstrom, und "
      "jede Vervollständigung ist ihr eigenes Orthogonal, also fehlt links ein ungerades Feld genau dann, wenn der "
      "Mischstrom dort lokal ist. Linke und rechte Parität stimmen in jedem Feld überein (aus Spin in Z/2), einseitige "
      "Markenfelder gibt es in keiner der 27, und der Mischstrom ist immer auf beiden Seiten zugleich lokal oder nicht. "
      f"Markenfelder tragen {ZAHLEN['A2_mit_markenfeld (alle, mit T)'][0]} der 27",
      lemma and beidseitig and nie_einseitig and dual and ZAHLEN["A2_mit_markenfeld (alle, mit T)"] == [12, 4])

T10 = [L for L in COMP if t_inv(L)]
mit_marke = [L for L in T10 if marke(L)]
sp_marke = {spektrum(L) for L in mit_marke}
iso = sorted(tuple(sorted(dict(L)[a] for a in ((1, 0), (0, 1)))) for L in mit_marke)
check("A3 [E] Satz: Unter A* (S und T, 10 der 27) tragen genau vier Vervollständigungen ein Markenfeld, und alle vier "
      "sind Graphen: jeder der 16 Sektoren genau einmal, gepaart durch (x, y) -> (+-x, +-y). Sie haben dasselbe Spektrum; "
      "es sind die Eichsummen der Paaruhr aus z178 (für die vier Wahlen von Spiegelbit und Ladungskonjugation rechts). "
      "A* und ein einziges Markenfeld geben also die diagonale Vervollständigung mit allen 16 Uhrsektoren, ohne "
      "Superauswahl und ohne GNS-Regel als Eingabe",
      len(T10) == 10 and len(mit_marke) == 4 and all(graph(L) for L in mit_marke) and len(sp_marke) == 1
      and iso == sorted(tuple(sorted(((sx, 0), (0, sy)))) for sx in (1, 3) for sy in (1, 3)))


def eichsumme(L):
    """Gibt es eine Paaruhr (U = i^(k - k') oder U' = i^(k + k'), Kleber s, t), deren Superauswahl mit GNS-Regel L gibt?"""
    for s, t, sig in product((1, 3), (1, 3), (1, 3)):
        gl = [(k % 4, (s * k) % 4) for k in range(4)]
        gr = [(k % 4, (t * k) % 4) for k in range(4)]
        untw = frozenset((a, b) for a in gl for b in gr)
        neutral = frozenset(((k % 4, (s * k) % 4), ((sig * k) % 4, (sig * t * k) % 4)) for k in range(4))
        if L & untw == neutral:
            return True
    return False


lokal16 = lambda L: any(g[0] in ((1, 0), (3, 0)) and g[1] in ((1, 0), (3, 0)) for g in L)
lokal4 = lambda L: any(g[0] in ((0, 1), (0, 3)) and g[1] in ((0, 1), (0, 3)) for g in L)
tab = [(eichsumme(L), marke(L), ohne_misch(L), graph(L), lokal16(L), lokal4(L), spektrum(L) in sp_marke) for L in T10]
gleichwertig = all(len(set(row)) == 1 for row in tab)
alle27 = [(eichsumme(L), marke(L), graph(L), lokal16(L), lokal4(L)) for L in COMP]
ohne_T = {"Eichsumme": sum(r[0] for r in alle27), "Markenfeld": sum(r[1] for r in alle27), "Graph": sum(r[2] for r in alle27),
          "16 lokal": sum(r[3] for r in alle27), "4 lokal": sum(r[4] for r in alle27)}
ZAHLEN["A4_ohne_T (Anzahl unter den 27)"] = ohne_T
check("A4 [E] Gleichwertig unter A* (jede der zehn Vervollständigungen erfüllt alle oder keine): (i) Superauswahl einer "
      "Paaruhr mit GNS-Regel (z178), (ii) ein Markenfeld, (iii) kein einseitiger Strom mischt die Deckhälften, (iv) Graph, "
      "(v) die 16 ist ein lokales Skalarfeld des Paars (h = h-quer = 5/8, Skalendimension 5/4), (vi) die 4 ist eines "
      "(Dimension 3/4), (vii) das Spektrum der Eichsumme. Ohne T fallen die Aussagen auseinander "
      f"({ohne_T}); das Lemma (ii) <=> (iii) gilt immer",
      gleichwertig and sum(r[0] for r in tab) == 4)
befund("Teil A: Der für K19 nötige Teil von P_std (Superauswahl der Paaruhr, mit der GNS-Regel) ist unter A* gleichwertig "
       "mit einer Lokalitätsaussage ohne Bezug auf die Standardform: Kein einseitiger Strom mischt die Deckhälften, "
       "gleichbedeutend: die 16 ist ein lokales Feld des Paars", gleichwertig and lemma)

# =================================================================== Teil B
say("\nTeil B: die Gegenmodelle unter A* und was sie ausschließt")


def typ(L):
    n = len(links(L))
    return {1: "Graph (so(10) + su(4), 60 Ströme je Seite)", 2: "so(16) (120 Ströme je Seite)",
            4: "E8 (x) E8-quer (248 Ströme je Seite)"}[n]


typen = {}
for L in T10:
    typen.setdefault(typ(L), []).append(L)
spektren = {}
for L in T10:
    spektren.setdefault(spektrum(L), []).append(L)
ZAHLEN["B1_typen"] = {k: len(v) for k, v in typen.items()}
check(f"B1 [E] Die zehn Vervollständigungen unter A* nach chiraler Algebra: {ZAHLEN['B1_typen']}. Es sind drei Spektren. "
      "Bosonische einseitige Erweiterungen gibt es nur drei (q = 0 auf der Untergruppe): keine, um den Mischstrom (2, 2) zu "
      "so(16), um einen Kleber zu E8. Die sechs Modelle ohne Eichsumme sind die Gegenmodelle zu P_std: Sie erfüllen A*, und "
      "in ihnen ist die Deckparität eine Symmetrie, keine Redundanz",
      ZAHLEN["B1_typen"] == {"Graph (so(10) + su(4), 60 Ströme je Seite)": 4, "so(16) (120 Ströme je Seite)": 2,
                             "E8 (x) E8-quer (248 Ströme je Seite)": 4} and len(spektren) == 3
      and sorted(len(v) for v in spektren.values()) == [2, 4, 4])


def charaktere(tau):
    qn = mp.exp(mp.j * mp.pi * tau)
    t2, t3, t4 = mp.jtheta(2, 0, qn), mp.jtheta(3, 0, qn), mp.jtheta(4, 0, qn)
    eta = mp.exp(mp.j * mp.pi * tau / 12) * mp.qp(mp.exp(2 * mp.j * mp.pi * tau))

    def so(n):
        return {0: (t3 ** n + t4 ** n) / (2 * eta ** n), 2: (t3 ** n - t4 ** n) / (2 * eta ** n),
                1: t2 ** n / (2 * eta ** n), 3: t2 ** n / (2 * eta ** n)}
    c5, c3 = so(5), so(3)
    return {(x, y): c5[x] * c3[y] for (x, y) in A}, (t2, t3, t4, eta)


def zsumme(L, tau):
    ch, _ = charaktere(tau)
    return sum(ch[g[0]] * mp.conj(ch[g[1]]) for g in L)


TAU = mp.mpf("0.23") + mp.mpf("0.81") * mp.j
inv_fehler = max(max(abs(zsumme(L, TAU) - zsumme(L, -1 / TAU)), abs(zsumme(L, TAU) - zsumme(L, TAU + 1))) / abs(zsumme(L, TAU))
                 for L in T10)
werte = {}
for name, Ls in typen.items():
    zi = {mp.nstr(mp.re(zsumme(L, mp.j)), 14) for L in Ls}
    z2 = {mp.nstr(mp.re(zsumme(L, 2 * mp.j)), 14) for L in Ls}
    zt = {mp.nstr(mp.re(zsumme(L, TAU)), 14) for L in Ls}
    werte[name] = (sorted(zi), sorted(z2), sorted(zt))
G_, S_, E_ = (typen[k][0] for k in ("Graph (so(10) + su(4), 60 Ströme je Seite)", "so(16) (120 Ströme je Seite)",
                                    "E8 (x) E8-quer (248 Ströme je Seite)"))
ZG = zsumme(G_, mp.j)
ch_i, _ = charaktere(mp.j)
gew_marke = sum(abs(ch_i[a]) ** 2 for a in ungerade) / ZG
achse = max(abs(zsumme(S_, y * mp.j) - zsumme(E_, y * mp.j)) / abs(zsumme(E_, y * mp.j))
            for y in (mp.mpf("0.4"), mp.mpf("0.7"), 1, mp.mpf("1.3"), 2, 3))
abseits = abs(zsumme(S_, TAU) - zsumme(E_, TAU)) / abs(zsumme(E_, TAU))
u_, v_ = mp.mpf("0.37"), mp.mpf("1.91")                        # (a + b + c)^2 = 2 (a^2 + b^2 + c^2) mit Jacobi
jacobi = abs((u_ ** 2 + (u_ + v_) ** 2 + v_ ** 2) ** 2 - 2 * (u_ ** 4 + (u_ + v_) ** 4 + v_ ** 4))
ZAHLEN["B2_Z(i), Z(2i), Z(0,23 + 0,81 i) je Typ"] = werte
ZAHLEN["B2_gewicht_markensektoren"] = mp.nstr(mp.re(gew_marke), 15)
ZAHLEN["B2_so16_gegen_E8E8 (imaginäre Achse, abseits)"] = [mp.nstr(achse, 3), mp.nstr(abseits, 6)]
check("B2 [N, 40 Stellen] Zustandssummen: Alle zehn sind S- und T-invariant (rel. Fehler "
      f"{mp.nstr(inv_fehler, 3)} bei tau = 0,23 + 0,81 i). Werte bei tau = i, 2 i und 0,23 + 0,81 i: {werte}. Die Eichsumme "
      f"hat bei tau = i den Wert 80 + 48 sqrt(2), die beiden anderen Spektren 144. Die so(16)-Form und E8 (x) E8-quer haben "
      f"auf der ganzen imaginären Achse dieselbe Zustandssumme (rel. Abweichung {mp.nstr(achse, 3)} an sechs Punkten; mit "
      "theta_3^4 = theta_2^4 + theta_4^4 ist (a + b + c)^2 = 2 (a^2 + b^2 + c^2) für die achten Potenzen) und unterscheiden "
      f"sich nur abseits davon (rel. {mp.nstr(abseits, 4)}), also im Spin, nicht in der Energie: 496 Zustände der "
      "Gesamtenergie 1 sind einmal 248 + 248 Ströme, einmal 120 + 120 Ströme und 256 Felder (1/2, 1/2). Nur in der Eichsumme "
      f"tragen die Sektoren mit Deckparität Gewicht: (12 sqrt(2) - 13)/14 = {ZAHLEN['B2_gewicht_markensektoren']} (z176, z178)",
      inv_fehler < mp.mpf(10) ** -30 and abs(ZG - (80 + 48 * mp.sqrt(2))) < mp.mpf(10) ** -30
      and abs(gew_marke - (12 * mp.sqrt(2) - 13) / 14) < mp.mpf(10) ** -30
      and all(len(x) == 1 for v in werte.values() for x in v)
      and len({v[2][0] for v in werte.values()}) == 3 and achse < mp.mpf(10) ** -30 and abseits > mp.mpf(10) ** -4
      and jacobi < mp.mpf(10) ** -30 and abs(zsumme(E_, mp.j) - 144) < mp.mpf(10) ** -30
      and 2 * 248 == 2 * 120 + 16 * 16)


def bild(L, fx, fy):
    f = lambda a: ((fx * a[0]) % 4, (fy * a[1]) % 4)
    return frozenset((f(g[0]), f(g[1])) for g in L)


sp_wurzel = [r for r in WURZELN if r[0].denominator == 2]
gesp = [tuple([-r[0]] + list(r[1:])) for r in WURZELN]
bleibt = sum(g in WSET for g in gesp)
gesp8 = [tuple(list(r[:7]) + [-r[7]]) for r in WURZELN]
refl = {name: (sum(bild(L, 3, 1) == L for L in Ls), sum(bild(L, 1, 3) == L for L in Ls), len(Ls)) for name, Ls in typen.items()}
ZAHLEN["B3_invariant unter Spiegelung (gerade Hälfte, ungerade Hälfte, Anzahl)"] = refl
check("B3 [E] Orientierungswechsel einer Schnittbasis (Spiegelung eines Elters, auf beiden Seiten zugleich): Im Gitter "
      f"bleiben von den 240 Wurzeln {bleibt} in E8 (die von so(16)), die 128 Spinorwurzeln gehen in die andere Klasse, also "
      "E8 in die zweite, isomorphe Verklebung (das Bit P_IF); dasselbe für einen Elter der ungeraden Hälfte. Auf den "
      f"Sektoren ist das x -> -x bzw. y -> -y. Invariant sind {refl}: alle Graphen und beide so(16)-Formen, keine der vier "
      "faktorisierten. Gilt 'Basiswechsel der Schnitte sind Symmetrien' einschließlich der Orientierung, fällt "
      "E8 (x) E8-quer; die so(16)-Formen bleiben",
      bleibt == 112 and sum(g in WSET for g in gesp8) == 112 and len(sp_wurzel) == 128
      and refl == {"Graph (so(10) + su(4), 60 Ströme je Seite)": (4, 4, 4), "so(16) (120 Ströme je Seite)": (2, 2, 2),
                   "E8 (x) E8-quer (248 Ströme je Seite)": (0, 0, 4)})

iota = np.diag([1.0] * 10 + [-1.0] * 6)
komm, anti = 0, 0
for i in range(16):
    for j in range(i + 1, 16):
        X = np.zeros((16, 16))
        X[i, j], X[j, i] = 1, -1
        komm += np.allclose(iota @ X - X @ iota, 0)
        anti += np.allclose(iota @ X + X @ iota, 0)
v_uhr = tuple([Fr(0)] * 5 + [H2] * 3)
phasen = {}
for r in WURZELN:
    p = sum(a * b for a, b in zip(r, v_uhr)) % 1
    phasen[str(p)] = phasen.get(str(p), 0) + 1
ZAHLEN["B4_phasen_der_wurzeln_unter_der_deckparitaet"] = phasen
check("B4 [E] Was der Mischstrom bedeutet: In so(16) vertauschen 45 + 15 = 60 Erzeuger mit der Deckparität, die 60 "
      f"Erzeuger (10, 6) antivertauschen mit ihr ({komm}, {anti}). Auf E8 hat die Fortsetzung der Deckparität die Phasen "
      f"{phasen} (als Bruchteil eines Umlaufs): fest sind 52 Wurzeln und der Rang, die Fixalgebra ist so(10) + so(6) mit 60 "
      "Erzeugern (z179). 'Kein einseitiger Strom mischt die Deckhälften' heißt: Die Teilung 10 | 6 wird von jeder stetigen "
      "Symmetrie des Paars erhalten. In den so(16)-Formen und in E8 (x) E8-quer dreht eine stetige Symmetrie die Teilung in "
      "jede andere; dort zeichnet nichts die Teilung des Bestands aus",
      komm == 60 and anti == 60 and phasen == {"0": 52, "1/2": 60, "1/4": 64, "3/4": 64})

tausch = lambda L: frozenset((g[1], g[0]) for g in L)
tau_inv = {name: (sum(tausch(L) == L for L in Ls), len(Ls)) for name, Ls in typen.items()}
ZAHLEN["B5_invariant unter Seitentausch"] = tau_inv
check("B5 [E] Zeitneutralität (Orientierungswechsel der Zyklen, tauscht Naht und Partner): invariant sind "
      f"{tau_inv}. Sie schließt nur die zwei faktorisierten Formen mit verschiedenem Kleber links und rechts aus; "
      "E8 (x) E8-quer mit gleichem Kleber und die so(16)-Formen bleiben",
      tau_inv == {"Graph (so(10) + su(4), 60 Ströme je Seite)": (4, 4), "so(16) (120 Ströme je Seite)": (2, 2),
                  "E8 (x) E8-quer (248 Ströme je Seite)": (2, 4)})
befund("Teil B: Die Gegenmodelle sind vollständig: sechs unter A*, in zwei Spektren. Den Typ E8 (x) E8-quer schließt der "
       "Orientierungswechsel der Schnittbasen aus, den Typ so(16) nur die Aussage selbst (kein Mischstrom)", True)

# =================================================================== Teil C
say("\nTeil C: Liefert der Bestand die Aussage?")
H1 = [(0, 0), (1, 0), (0, 1), (1, 1)]
S_mod = lambda d: (d[1] % 2, d[0] % 2)                         # S = [[0, -1], [1, 0]] modulo 2
T_mod = lambda d: ((d[0] + d[1]) % 2, d[1] % 2)                # T = [[1, 1], [0, 1]] modulo 2
bahnen = []
rest = set(H1)
while rest:
    o, todo = set(), [min(rest)]
    while todo:
        d = todo.pop()
        if d not in o:
            o.add(d)
            todo += [S_mod(d), T_mod(d)]
    bahnen.append(sorted(o))
    rest -= o
_, (t2, t3, t4, eta) = charaktere(TAU)
th = [abs(t2), abs(t3), abs(t4)]
z_gemeinsam = sum(t ** 16 for t in th) / (2 * abs(eta) ** 16)
z_getrennt = sum(t ** 10 for t in th) * sum(t ** 6 for t in th) / (4 * abs(eta) ** 16)
z_so16 = zsumme(typen["so(16) (120 Ströme je Seite)"][0], TAU)
z_graph = zsumme(typen["Graph (so(10) + su(4), 60 Ströme je Seite)"][0], TAU)
ZAHLEN["C1_bahnen_der_relativen_spinstruktur"] = [str(b) for b in bahnen]
check("C1 [E/N] A allein: Die relative Spinstruktur delta der beiden Deckhälften liegt in H^1(Torus, Z2), und SL(2, Z) hat "
      f"darauf die Bahnen {ZAHLEN['C1_bahnen_der_relativen_spinstruktur']}. delta = 0 ist für sich invariant: 'Keine "
      "ausgezeichneten Zyklen' verlangt die Summe über die Spinstruktur, aber nicht die Summe über delta. Die gemeinsame "
      "Summe (delta = 0) ist die so(16)-Form, die Summe über alle delta die Eichsumme (rel. Fehler "
      f"{mp.nstr(max(abs(z_gemeinsam - z_so16) / abs(z_so16), abs(z_getrennt - z_graph) / abs(z_graph)), 3)}). "
      "Beide erfüllen A*",
      sorted(len(b) for b in bahnen) == [1, 3] and [(0, 0)] in bahnen
      and abs(z_gemeinsam - z_so16) / abs(z_so16) < mp.mpf(10) ** -30
      and abs(z_getrennt - z_graph) / abs(z_graph) < mp.mpf(10) ** -30)

# Überlagerung: p(z) - e_k = F_k(z)^2 mit F_k = pi theta_a theta_b theta_(k+1)(pi z)/theta_1(pi z), tau = i
QI = mp.exp(-mp.pi)
T2, T3, T4 = mp.jtheta(2, 0, QI), mp.jtheta(3, 0, QI), mp.jtheta(4, 0, QI)
E_K = [mp.pi ** 2 * (T2 ** 4 + 2 * T4 ** 4) / 3, mp.pi ** 2 * (T2 ** 4 - T4 ** 4) / 3, -mp.pi ** 2 * (2 * T2 ** 4 + T4 ** 4) / 3]


def F(k, z):
    pref = [T3 * T4, T2 * T4, T2 * T3][k]
    return mp.pi * pref * mp.jtheta(k + 2, mp.pi * z, QI) / mp.jtheta(1, mp.pi * z, QI)


def wp(z):
    return E_K[0] + F(0, z) ** 2


def wp_strich(z):
    return -2 * F(0, z) * F(1, z) * F(2, z)


z0 = mp.mpf("0.137") + mp.mpf("0.291") * mp.j
konsistent = max(abs(E_K[k] + F(k, z0) ** 2 - wp(z0)) for k in range(3))
abl = abs(mp.diff(wp, z0) - wp_strich(z0)) / abs(wp_strich(z0))
ungerade_fkt = abs(wp_strich(-z0) + wp_strich(z0)) / abs(wp_strich(z0))
period = max(abs(wp_strich(z0 + 1) - wp_strich(z0)), abs(wp_strich(z0 + mp.j) - wp_strich(z0))) / abs(wp_strich(z0))


def windung(k, weg, n=1200):
    """Windungszahl von F_k längs des geschlossenen Wegs; p - e_k windet doppelt so oft."""
    tot, prev = mp.mpf(0), mp.arg(F(k, weg(mp.mpf(0))))
    for i in range(1, n + 1):
        cur = mp.arg(F(k, weg(mp.mpf(i) / n)))
        d = cur - prev
        while d > mp.pi:
            d -= 2 * mp.pi
        while d < -mp.pi:
            d += 2 * mp.pi
        tot += d
        prev = cur
    return tot / (2 * mp.pi)


mp.mp.dps = 20
wege = {"a (waagerecht)": lambda t: t + mp.mpf("0.2") * mp.j, "b (senkrecht)": lambda t: mp.mpf("0.2") + t * mp.j,
        "a + b (diagonal)": lambda t: mp.mpf("0.2") + t * (1 + mp.j)}
wind = {}
for name, weg in wege.items():
    w = [2 * windung(k, weg) for k in range(3)]                # Windung des Bildes um e1, 0, -e1
    wind[name] = [int(mp.nint(x)) for x in w]
    assert max(abs(x - mp.nint(x)) for x in w) < mp.mpf("1e-6")
mp.mp.dps = 40
ZAHLEN["C2_windung des Bildes um die Marken (e1, 0, -e1)"] = wind
ZAHLEN["C2_e1"] = mp.nstr(E_K[0], 12)
check("C2 [N] Überlagerungsgeometrie (tau = i, p - e_k = F_k^2 geprüft, Fehler "
      f"{mp.nstr(konsistent, 3)}; e1 = {ZAHLEN['C2_e1']}, e2 = 0): Das Bild eines Zyklus des Nahttorus umläuft die endlichen "
      f"Marken (e1, 0, -e1) mit den Windungszahlen {wind}. Jede der drei Zyklusklassen umläuft zwei Marken ungerade oft, "
      "die Gesamtwindung ist gerade; die drei Klassen sind die drei Paarungen der Marken (z180 A2). Ein unter der "
      "Deckinvolution ungerader Schnitt ist deshalb auf dem Nahttorus einwertig (p' ist ungerade und periodisch, Fehler "
      f"{mp.nstr(max(ungerade_fkt, period), 3)}). Die Überlagerung liefert also keine relative Verdrillung der beiden "
      "Deckhälften: Aus der Geometrie folgt delta = 0, die gemeinsame Summe, nicht die Eichsumme",
      konsistent < mp.mpf(10) ** -30 and abl < mp.mpf(10) ** -15 and max(ungerade_fkt, period) < mp.mpf(10) ** -30
      and abs(E_K[1]) < mp.mpf(10) ** -35
      and all(sum(abs(x) % 2 for x in w) == 2 and sum(w) % 2 == 0 for w in wind.values())
      and len({tuple(abs(x) % 2 for x in w) for w in wind.values()}) == 3)

z170 = json.loads((HERE / "z170_ergebnis.json").read_text())
z170_ok = any(k.startswith("C2") and v for k, v in z170["checks"].items()) \
    and any(k.startswith("B4") and v for k, v in z170["checks"].items())
nur_graph = all(lokal16(L) == graph(L) and lokal4(L) == graph(L) for L in T10)
kein_beides = not any(lokal16(L) and M_L in L for L in COMP)
check("C3 [E] Gelöster alpha-Block: Die Antwortpotenz 5/4 hat als einzigen nativen positiven Träger den Kohärenzkern der "
      "16 (h = 5/8, z170 C2), und 3/4 gehört zur 4 (z170 B4). Diese Felder haben h in ungeraden Achteln, sind also "
      "Markenfelder: Unter A* sind sie genau in den vier Graphen lokal. Keine der 27 Vervollständigungen enthält die 16 als "
      "lokales Feld und zugleich den Mischstrom. In den sechs Gegenmodellen hat die Antwortpotenz daher keinen lokalen "
      "Träger; die 16 ist dort das Ende einer Defektlinie. Die Gleichsetzung der Potenz mit dem Feld ist in z170 eine "
      "Lesart, kein Satz: Sie stützt P_std, leitet es aber nicht her",
      z170_ok and nur_graph and kein_beides)
befund("Teil C: A*, die Zeitneutralität und die Überlagerungsgeometrie legen die Redundanz der Deckparität nicht fest; "
       "die Geometrie allein gibt sogar die gemeinsame Summe", True)
befund("Der gelöste alpha-Block verlangt P_std (als Satz)", False)

# =================================================================== Teil D
say("\nTeil D: Bilanz")
import reconstruction_search as rs  # noqa: E402  (nur lesend)
dist = rs.complete_solution_distance()
ZAHLEN["D1_abstand"] = dist["distance"]
ZAHLEN["D1_offen"] = dist["open_total"]
kette = [
    "[Bestand] A*: Die Zustandssumme des Paars ist S- und T-invariant; zehn Vervollständigungen, drei Spektren (B1)",
    "[Lemma, A2] Markenfeld <=> kein einseitiger Mischstrom (10, 6); exakt für alle 27",
    "[Satz, A3, A4] A* und eine der gleichwertigen Aussagen (Superauswahl der Paaruhr mit GNS-Regel; Markenfeld; kein "
    "Mischstrom; die 16 als lokales Feld) geben die Eichsumme mit allen 16 Uhrsektoren (K19)",
    "[Gegenmodelle, B1 bis B5] vier Formen E8 (x) E8-quer, zwei Formen so(16): A* erfüllt, Deckparität nur Symmetrie",
    "[Sperre, C1, C2] A lässt delta = 0 für sich zu, und die Überlagerung gibt delta = 0",
    "[Stütze, C3] Der Träger der Antwortpotenz 5/4 ist nur in der Eichsumme ein lokales Feld (Lesart aus z170)",
    "[offen] P_std bleibt Annahme; sein zweiter Teil (P_uhr: die Phase der Komposition ist Funktion der Uhrwerte) ist hier "
    "nicht berührt",
]
for k in kette:
    say("   " + k)
check("D1 [E] Abstand zur Komplettlösung nach der Runde (OriginLab): "
      f"{dist['distance']}, {dist['open_total']} unabhängige offene Punkte; K19 und K20 ruhen weiter auf P_std (und "
      "P_zentral). Kein Punkt ist geschlossen",
      dist["open_total"] == 8 and set(dist["bench_derived"]) == {"K19", "K20"})
befund("P_std ist aus A*, der Zeitneutralität oder der Überlagerungsgeometrie hergeleitet", False)
befund("P_std (Teil Paaruhr) ist unter A* gleichwertig mit 'kein einseitiger Strom mischt die Deckhälften' und mit 'die 16 "
       "ist ein lokales Feld des Paars'", gleichwertig and lemma)
befund("Die Gegenmodelle zu P_std sind vollständig aufgezählt (sechs, zwei Spektren) und physikalisch verschieden von der "
       "Eichsumme", len(spektren) == 3)
befund("Komplettlösung", False)

urteil = (
    "P_std ist angegriffen und nicht gefallen, aber eingegrenzt. (1) Neue exakte Form. Die 16 Eltern sind die "
    "Vektorgewichte von D5 + D3; ein Feld im Sektor (x, y) hat mit ihnen die Monodromie (-1)^x und (-1)^y. Die Eltern sehen "
    "die Deckparität genau um die acht Sektoren mit h in ungeraden Achteln. Lemma (alle 27): Eine Vervollständigung enthält "
    "ein solches Markenfeld genau dann, wenn der Mischstrom (10, 6) auf keiner Seite lokal ist. Satz: Unter A* sind das "
    "genau die vier Graphen mit einem Spektrum, die Eichsummen aus z178. Der für K19 nötige Teil von P_std (Superauswahl "
    "der Paaruhr, mit der GNS-Regel) ist damit unter A* gleichwertig mit jeder der Aussagen: kein einseitiger Strom mischt "
    "die Deckhälften; die 16 ist ein lokales Feld des Paars; die 4 ist eines. Die Standardform wird dafür nicht gebraucht. "
    "(2) Gegenmodelle, vollständig. Unter A* gibt es zehn Vervollständigungen in drei Spektren: vier Eichsummen, zwei "
    "so(16)-Formen, vier Formen E8 (x) E8-quer. Bei tau = i haben die Eichsummen Z = 80 + 48 sqrt(2), die anderen 144; "
    "die so(16)-Form und E8 (x) E8-quer stimmen auf der ganzen imaginären Achse überein und unterscheiden sich nur im Spin. "
    "Die sechs Gegenmodelle erfüllen A*; in ihnen ist die Deckparität eine Symmetrie. E8 (x) E8-quer fällt, wenn Basiswechsel "
    "der Schnitte einschließlich des Orientierungswechsels Symmetrien des Paars sind (eine Spiegelung führt E8 in die "
    "andere Verklebung). Die so(16)-Formen fallen nur durch die Aussage selbst. Die Zeitneutralität schließt zwei der zehn "
    "aus, keines der beiden Spektren. (3) Der Bestand entscheidet nicht. 'Keine ausgezeichneten Zyklen' lässt die "
    "relative Spinstruktur delta = 0 für sich zu. Die Überlagerung der vier Marken gibt sogar delta = 0: Jeder Zyklus des "
    "Nahttorus umläuft zwei Marken, ungerade Schnitte sind dort einwertig. Die getrennte Summe kommt also nicht aus der "
    "Geometrie; sie verlangt, dass die Deckparität eine Redundanz ist. (4) Stütze aus einem gelösten Block: Der native "
    "Träger der Antwortpotenz 5/4 (die 16, z170) ist nur in der Eichsumme ein lokales Feld; in allen Gegenmodellen ist er "
    "das Ende einer Defektlinie. Das ist eine Lesart, kein Satz. P_std bleibt Annahme, in schärferer Form: 'Die Teilung "
    "10 | 6 wird von jeder stetigen Symmetrie des Paars erhalten.' Der zweite Teil von P_std (P_uhr) ist nicht berührt. "
    "Abstand unverändert, acht offene Punkte, nichts befördert.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z183_ergebnis.json").write_text(json.dumps({"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil},
                                                    indent=1, ensure_ascii=False, default=str))
(HERE / "ausgabe_z183.txt").write_text("\n".join(OUT) + "\n")
print(f"(Laufzeit {time.time() - T_START:.0f} s)")
