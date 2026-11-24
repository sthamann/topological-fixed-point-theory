"""Z64: Das zweite O(1) -- folgt Lesart T (Naht = Twistorgerade) aus den TFPT-Daten? (Lücke B der Vogelperspektive)

Leitlinie 'Rekonstruktion statt Auswahl': keine neuen Auswahlprinzipien, nur Konstruktionen aus dem Bestand. Eine Twistorgerade
braucht das Normalenbündel N = O(1) + O(1): Rang 2, Grad 2, balanciert (h0(N(-2)) = 0), also h0 = 4, h1 = 0. z61 fand für die
Doppelüberlagerung O(1) (x) pi_* O_E = O(1) + O(-1). Geprüft wird systematisch, welche kanonischen Konstruktionen aus dem Bestand
(Naht P(V) mit O(1), Euler-Folge, vier Marken mit D4 > V4, E: w^2 = 1 - z^4, (Z2)^3-Überlagerung mit r_* r^* O(1) =
O(1) + 6 O + O(-1), Deckhälften, Spinsektoren, Standardform) O(1) + O(1) liefern, welche kanonisch sind, und was das für T heißt.

Werkzeuge (exakt): Riemann-Roch-Räume L(D) auf E (Funktionen f(z) + w g(z), Rang über Q(i) mit flint), Spaltungstyp aus h0
der Twists; Bündel auf P^1 über Übergangsmatrizen (Jets, Euler-Erweiterung) und über Hecke-Punktbedingungen.

Teil 0  Werkzeuge: O(a) + O(b), Riemann-Roch auf E an 80 Divisoren, lineare Äquivalenzen der Sonderpunkte.
Teil 1  Doppelüberlagerung: Satz 'Sprung genau beim Pullback'; die z61-Klasse (Blätter mit zurückgezogenem Bündel) ist nie
        balanciert; Spinstrukturen: T_C (x) pi_* theta = O(1) + O(1) genau für gerade theta (Szegő-Kern), O(2) + O für theta = O.
Teil 2  (Z2)^3-Überlagerung: Sektoren von r_* L^k mit L = r^* O(1); kein O(1) + O(1) in r_* L; in r_* L^2 = r_* K_{E/B} genau ein
        D4-invarianter balancierter Sektor (antipodal); Nicht-Pullback-Spinlifts sind balanciert, haben aber kein Dublett.
Teil 3  Hecke-Modifikationen an den Marken: O(1) + O(1) genau an zwei Marken; keine Zweiermenge ist V4- oder D4-invariant.
Teil 4  Partner und Dublett: J^1(O(2)) = O(1) (x) H0(O(1)); Euler-Kanal; Grad des Partners; Lorentzkriterium (nur W = V-bar^*);
        Schur (h durch die Marken erzwungen); Blattblindheit: kein Überlagerungssektor trägt eine lorentzkovariante Faser.
Teil 5  Bilanz und Status von T.
"""
import itertools
import json
import random
from fractions import Fraction as Fr
from math import comb
from pathlib import Path

import flint
import mpmath as mp
import numpy as np

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN, STATUS = {}, {}, {}, {}
rng = np.random.default_rng(64)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def status(key, lab, text):
    STATUS[key] = {"status": lab, "aussage": text}
    say(f"   [{lab}] {key}: {text}")


# ------------------------------------------------------------------ Werkzeuge: Q(i) und exakter Rang
class Q:
    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a, self.b = Fr(a), Fr(b)

    def __add__(s, o):
        o = q(o)
        return Q(s.a + o.a, s.b + o.b)
    __radd__ = __add__

    def __sub__(s, o):
        o = q(o)
        return Q(s.a - o.a, s.b - o.b)

    def __mul__(s, o):
        o = q(o)
        return Q(s.a * o.a - s.b * o.b, s.a * o.b + s.b * o.a)
    __rmul__ = __mul__

    def __neg__(s):
        return Q(-s.a, -s.b)

    def __pow__(s, n):
        r = Q(1)
        for _ in range(n):
            r = r * s
        return r


def q(x):
    return x if isinstance(x, Q) else Q(x)


def fq(f):
    return flint.fmpq(f.numerator, f.denominator)


def rank_gauss(rows, ncols):
    """Rang über Q(i) = halber Rang der reellen Blockmatrix [[A, -B], [B, A]] (exakt, flint)."""
    if not rows:
        return 0
    ent = []
    for r in rows:
        ent += [fq(x.a) for x in r] + [fq(-x.b) for x in r]
    for r in rows:
        ent += [fq(x.b) for x in r] + [fq(x.a) for x in r]
    rk = flint.fmpq_mat(2 * len(rows), 2 * ncols, ent).rank()
    assert rk % 2 == 0
    return rk // 2


# ------------------------------------------------------------------ Werkzeuge: L(D) auf E: w^2 = 1 - z^4
MARKS = {"p1": Q(1), "pi": Q(0, 1), "pm1": Q(-1), "pmi": Q(0, -1)}       # Verzweigungspunkte über den vier Marken
PTS = list(MARKS) + ["inf+", "inf-", "o+", "o-"]                           # inf+-: w/z^2 -> +-i; o+- = (0, +-1)
PB1 = {"inf+": 1, "inf-": 1}                                                # pi^* O(1)
R4 = {a: 1 for a in MARKS}                                                  # Verzweigungsdivisor R


def sqrt_series(n):
    c = [Fr(0)] * (n + 1)
    b, k = Fr(1), 0
    while 4 * k <= n:
        c[4 * k] = b * (-1) ** k
        b = b * (Fr(1, 2) - k) / (k + 1)
        k += 1
    return c


def ell(D, extra=0, part=None):
    """l(D) für D mit Träger in den Sonderpunkten. Ansatz phi = (P + w Q)/(prod (z - a)^e_a z^e0); an den Marken trennen
    sich gerade (P) und ungerade (w Q) Ordnungen exakt, an inf+- und o+- Reihen in t = 1/z bzw. z. part = 'gerade'/'ungerade'
    beschränkt auf den Deck-Eigenteil (nur für iota-invariante D sinnvoll)."""
    n = {p: D.get(p, 0) for p in PTS}
    e = {a: max(0, (n[a] + 1) // 2) for a in MARKS}
    e0 = max(0, n["o+"], n["o-"])
    dd = sum(e.values()) + e0
    N = dd + max(0, n["inf+"], n["inf-"]) + 3 + extra
    nun = 2 * (N + 1)
    s = sqrt_series(2 * N + dd + 20 + max(0, -min(n.values())))
    rows = []
    for a, val in MARKS.items():
        for off, mreq in ((0, e[a] - n[a] // 2), (N + 1, e[a] - (n[a] + 1) // 2)):
            for j in range(max(0, mreq)):
                row = [Q(0)] * nun
                for k in range(j, N + 1):
                    row[off + k] = Q(comb(k, j)) * val ** (k - j)
                rows.append(row)
    for eps, key in ((1, "o+"), (-1, "o-")):
        for j in range(max(0, e0 - n[key])):
            row = [Q(0)] * nun
            if j <= N:
                row[j] = row[j] + 1
            for l_ in range(j + 1):
                if j - l_ <= N and s[l_] != 0:
                    row[N + 1 + j - l_] = row[N + 1 + j - l_] + Q(eps) * s[l_]
            rows.append(row)
    for eps, key in ((1, "inf+"), (-1, "inf-")):
        for j in range(max(0, N + 2 - dd - n[key])):
            row = [Q(0)] * nun
            if 2 <= j <= N + 2:
                row[N - (j - 2)] = row[N - (j - 2)] + 1
            for l_ in range(j + 1):
                idx = N - (j - l_)
                if 0 <= idx <= N and s[l_] != 0:
                    row[N + 1 + idx] = row[N + 1 + idx] + Q(0, eps) * s[l_]
            rows.append(row)
    if part in ("gerade", "ungerade"):
        off = N + 1 if part == "gerade" else 0
        for k in range(N + 1):
            row = [Q(0)] * nun
            row[off + k] = Q(1)
            rows.append(row)
    return nun - rank_gauss(rows, nun)


def dadd(*Ds):
    out = {}
    for D in Ds:
        for k, v in D.items():
            out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v != 0}


def dsc(D, c):
    return {k: c * v for k, v in D.items()}


def deg(D):
    return sum(D.values())


def split_from_h0(h, rank):
    ms = sorted(h)
    cnt = {}
    for m in ms[2:]:
        c = h[m] - 2 * h[m - 1] + h[m - 2]
        if c:
            cnt[-m] = c
    if sum(cnt.values()) != rank or h[ms[0]] != 0 or h[ms[1]] != 0:
        return None
    return tuple(sorted([a for a, c in cnt.items() for _ in range(c)], reverse=True))


def split_push(D, extra=0):
    """Spaltungstyp von pi_* O_E(D) (Rang 2) aus h0 = l(D + m pi^* O(1))."""
    c = -(deg(D) - 2) // 2
    return split_from_h0({m: ell(dadd(D, dsc(PB1, m)), extra) for m in range(c - 4, c + 4)}, 2)


def split_part(D, part):
    """Grad des Deck-Eigenteils (Linienbündel) von pi_* O_E(D)."""
    c = -(deg(D) - 2) // 2
    return split_from_h0({m: ell(dadd(D, dsc(PB1, m)), part=part) for m in range(c - 5, c + 5)}, 1)


# ------------------------------------------------------------------ Werkzeuge: Bündel auf P^1
def h0_trans(G, m, Nmax=14):
    """h0 des Bündels mit s_inf = z^-m G(z) s_0 (G: Laurentpolynome als {Potenz: Koeffizient})."""
    r = len(G)
    nun = r * (Nmax + 1)
    rows = {}
    for i in range(r):
        for j in range(r):
            for pw, cf in G[i][j].items():
                for k in range(Nmax + 1):
                    p = pw + k - m
                    if p > 0:
                        rows.setdefault((i, p), [Q(0)] * nun)
                        rows[(i, p)][j * (Nmax + 1) + k] = rows[(i, p)][j * (Nmax + 1) + k] + cf
    return nun - rank_gauss(list(rows.values()), nun)


def split_trans(G, Nmax=14):
    return split_from_h0({m: h0_trans(G, m, Nmax) for m in range(-6, 7)}, len(G))


def h0_hecke(k1, k2, S, lines, m):
    """W = O(k1) e1 + O(k2) e2, Aufmodifikation an den Punkten S längs Geraden l_p: Schnitte sigma/prod(z - p) von W(m)."""
    d1, d2 = k1 + m + len(S), k2 + m + len(S)
    n1, n2 = max(d1 + 1, 0), max(d2 + 1, 0)
    rows = []
    for p, (l1, l2) in zip(S, lines):
        rows.append([q(l2) * p ** k for k in range(n1)] + [-q(l1) * p ** k for k in range(n2)])
    return n1 + n2 - rank_gauss(rows, n1 + n2)


def split_hecke(k1, k2, S, lines):
    c = -(k1 + k2 + len(S)) // 2
    return split_from_h0({m: h0_hecke(k1, k2, S, lines, m) for m in range(c - 5, c + 5)}, 2)


def bal2(t):
    return t is not None and sum(t) == 2 and t == (1, 1)


# ================================================================== Teil 0
say("Teil 0: Werkzeuge")
okt = True
for a_, b_ in [(1, 1), (2, 0), (1, -1), (0, -2), (3, -1), (0, 0)]:
    G = [[{-a_: Q(1)}, {}], [{}, {-b_: Q(1)}]]
    okt &= split_trans(G) == tuple(sorted((a_, b_), reverse=True)) == split_trans(G, Nmax=18)
check("Übergangsmatrix-Werkzeug: Spaltungstyp von O(a) + O(b) für sechs Paare richtig und stabil (Nmax 14 -> 18)", okt)
random.seed(64)
okrr, nrr = True, 0
for _ in range(80):
    D = {p: random.randint(-2, 3) for p in PTS}
    D = {k: v for k, v in D.items() if v}
    l_, d_ = ell(D), deg(D)
    okrr &= (l_ == d_) if d_ >= 1 else ((l_ == 0) if d_ < 0 else l_ in (0, 1))
    okrr &= l_ == ell(D, extra=3)
    nrr += 1
okrr &= ell({}) == 1 and ell({"p1": 1, "pm1": 1}) == 2 and ell(R4) == 4
check(f"Riemann-Roch auf E (Geschlecht 1): l(D) = deg D für deg >= 1, 0 für deg < 0, in {{0, 1}} für deg 0, an {nrr} "
      "Zufallsdivisoren auf den Sonderpunkten; stabil unter Ansatzvergrößerung", okrr)
eq = lambda D1, D2: ell(dadd(D1, dsc(D2, -1))) == 1                         # lineare Äquivalenz bei gleichem Grad
pairs = {"antipodal": ({"p1": 1, "pm1": 1}, {"pi": 1, "pmi": 1}), "seite1": ({"p1": 1, "pi": 1}, {"pm1": 1, "pmi": 1}),
         "seite2": ({"p1": 1, "pmi": 1}, {"pi": 1, "pm1": 1})}
ok_eq = all(eq({a: 2}, PB1) for a in MARKS) and eq({"o+": 1, "o-": 1}, PB1)
ok_eq &= all(eq(S_, Sc) for S_, Sc in pairs.values())
ok_eq &= not eq(pairs["antipodal"][0], pairs["seite1"][0]) and not eq(pairs["antipodal"][0], pairs["seite2"][0]) \
    and not eq(pairs["seite1"][0], pairs["seite2"][0]) and not eq(pairs["antipodal"][0], PB1)
zwei_inf = [k for k, (S_, _) in pairs.items() if eq({"inf+": 2}, S_)]
zwei_o = [k for k, (S_, _) in pairs.items() if eq({"o+": 2}, S_)]
ok_eq &= zwei_inf == ["antipodal"] and zwei_o == ["antipodal"] and eq({"inf+": 1, "p1": 1}, {"inf-": 1, "pm1": 1})
ZAHLEN["aequivalenzen"] = {"2 inf+ ~": zwei_inf, "2 o+ ~": zwei_o}
check("lineare Äquivalenzen: 2 p_a ~ inf+ + inf- ~ o+ + o- (= pi^* O(1)); Komplementärpaare äquivalent, die drei Paarklassen "
      "verschieden; 2 inf+ ~ 2 o+ ~ p_1 + p_-1 (antipodal), also inf+ - inf- ~ p_-1 - p_1", ok_eq)

# ================================================================== Teil 1
say("\nTeil 1: Doppelüberlagerung pi: E -> P^1 (z61-Klasse, Satz 'Sprung genau beim Pullback', Spinstrukturen)")
tab_pb = {k: split_push(dsc(PB1, k)) for k in range(-1, 4)}
ZAHLEN["pi_pi^*O(k)"] = {str(k): v for k, v in tab_pb.items()}
check("z61-Klasse: pi_* pi^* O(k) = O(k) + O(k-2) für k = -1..3 (Projektionsformel); jeder Twist mit O(j) der Naht verschiebt "
      "beide Summanden, die Lücke 2 bleibt: kein Bündel 'Blätter mit zurückgezogenem Bündel' ist O(1) + O(1)",
      all(tab_pb[k] == (k, k - 2) for k in tab_pb) and tab_pb[1] == (1, -1))
lists = {
    0: [{}, {"p1": 1, "pm1": -1}, {"p1": 1, "pi": -1}, {"p1": 1, "pmi": -1}, {"inf+": 1, "inf-": -1}, {"inf+": 1, "p1": -1},
        {"o+": 1, "p1": -1}, {"o+": 1, "o-": -1}, {"inf+": 1, "o+": -1}, {"p1": 2, "inf+": -1, "inf-": -1},
        {"o+": 1, "o-": 1, "inf+": -1, "inf-": -1}, {"p1": 1, "pm1": 1, "pi": -1, "pmi": -1}],
    2: [PB1, {"p1": 2}, {"o+": 1, "o-": 1}, {"p1": 1, "pm1": 1}, {"p1": 1, "pi": 1}, {"p1": 1, "pmi": 1}, {"inf+": 2},
        {"inf+": 1, "o+": 1}, {"o+": 1, "p1": 1}, {"o+": 2}, {"inf+": 1, "p1": 1}],
    4: [dsc(PB1, 2), R4, dadd(PB1, {"p1": 1, "pm1": 1}), dadd(PB1, {"p1": 1, "pi": 1}), {"inf+": 4}, {"inf+": 2, "o+": 2},
        {"p1": 3, "o+": 1}, dadd(PB1, {"o+": 2})]}
ok_s, cnt_pb, cnt_np, tab_s = True, 0, 0, {}
for d_, Ds in lists.items():
    k = d_ // 2
    for D in Ds:
        pull = ell(dadd(D, dsc(PB1, -k))) == 1
        sp_ = split_push(D)
        ok_s &= sp_ == ((k, k - 2) if pull else (k - 1, k - 1))
        cnt_pb += pull
        cnt_np += not pull
        tab_s[str(D)] = {"grad": d_, "pullback": pull, "spaltung": sp_}
ok_s &= all(split_push(D) == (0, -1) for D in ({"p1": 1}, {"inf+": 1}, {"o+": 1})) and \
    all(split_push(D) == (1, 0) for D in ({"p1": 1, "pm1": 1, "pi": 1}, {"inf+": 2, "o+": 1}))
ZAHLEN["satz_sprung_pullback"] = tab_s
check(f"Satz: für M auf E vom Grad 2k ist pi_* M = O(k-1) + O(k-1) genau dann, wenn M nicht pi^* O(k) ist, sonst O(k) + O(k-2) "
      f"({cnt_pb} Pullbacks, {cnt_np} andere, Grad 0/2/4; ungerader Grad stets O(k) + O(k-1)): der einzige Sprungpunkt jeder "
      "Gradklasse ist das zurückgezogene Bündel", ok_s and cnt_pb >= 9 and cnt_np >= 20)
thetas = {"O (ungerade, P,P)": {}, "eta_antipodal": {"p1": 1, "pm1": -1}, "eta_seite1": {"p1": 1, "pi": -1},
          "eta_seite2": {"p1": 1, "pmi": -1}}
tab_t = {k: {"h0(theta)": ell(v), "T_C (x) pi_* theta": split_push(dadd(v, dsc(PB1, 2)))} for k, v in thetas.items()}
ZAHLEN["spinstrukturen_E0"] = tab_t
check("Spinstrukturen von E (K_E = O): theta = O hat die Nullmode h0 = 1 und T_C (x) pi_* theta = O(2) + O; die drei geraden "
      "theta (h0 = 0, Szegő-Kern existiert) geben T_C (x) pi_* theta = O(1) + O(1)",
      tab_t["O (ungerade, P,P)"] == {"h0(theta)": 1, "T_C (x) pi_* theta": (2, 0)} and
      all(tab_t[k] == {"h0(theta)": 0, "T_C (x) pi_* theta": (1, 1)} for k in list(thetas)[1:]))
mp.mp.dps = 30
qn = mp.exp(-mp.pi)
th = {k: mp.jtheta(k, 0, qn) for k in (1, 2, 3, 4)}
ok_th = abs(th[1]) < mp.mpf(10) ** -25 and abs(th[2] - th[4]) < mp.mpf(10) ** -25 and abs(th[3] - th[2]) > 0.1 and \
    abs(th[3] ** 4 - th[2] ** 4 - th[4] ** 4) < mp.mpf(10) ** -25
M_i = np.array([[0, 1], [1, 0]])                                            # i auf E[2] = (Z/2)^2 (Gitter Z[i])
fix2 = [v for v in [(1, 0), (0, 1), (1, 1)] if tuple((M_i @ np.array(v)) % 2) == v]
ok_th &= fix2 == [(1, 1)] and eq({"pi": 1, "pmi": 1}, {"p1": 1, "pm1": 1}) and not eq({"pi": 1, "pm1": 1}, {"p1": 1, "pi": 1})
ZAHLEN["theta_tau_i"] = {f"theta{k}(0)": mp.nstr(v, 12) for k, v in th.items()}
check("Szegő bei tau = i: theta_1(0) = 0 (ungerade, kein Kern), theta_2(0) = theta_4(0) != theta_3(0) (Jacobi erfüllt); die "
      "Vierteldrehung hält genau einen nichttrivialen 2-Torsionspunkt fest, und auf E bildet sie die Seitenklassen aufeinander "
      "ab und hält die antipodale fest: eta_* = O(p_1 - p_-1) = O(inf+ - inf-) ist die vierteldrehungsfeste gerade Struktur",
      ok_th)
sp_R = split_push(dadd(PB1, R4))
deg2_sums = {sum(split_push(D)) for D in lists[2]}
check("Kandidat (3) mit Verzweigungsdivisor: pi_*(pi^* O(1) (x) O(R)) = O(3) + O(1) (Grad 4); pi_* L mit deg L = 2 hat stets "
      "Grad 0: beide verfehlen Grad 2, erst die halbe Verzweigung (relative Spinstruktur) mit T_C trifft O(1) + O(1)",
      sp_R == (3, 1) and deg2_sums == {0})
parts = {nm: (split_part(D, "gerade"), split_part(D, "ungerade")) for nm, D in
         {"O": {}, "eta_antipodal": thetas["eta_antipodal"], "eta_seite1": thetas["eta_seite1"]}.items()}
ZAHLEN["deck_eigenteile"] = parts
check("Deck-Eigenteile: pi_* O = O (gerade) + O(-2) (ungerade); pi_* eta = O(-1) + O(-1) mit je einem Rang-1-Teil pro "
      "Deckparität: die Decktransformation wirkt auf jeder Faser eines geraden Sektors als diag(1, -1), nicht skalar",
      parts["O"] == ((0,), (-2,)) and parts["eta_antipodal"] == ((-1,), (-1,)) and parts["eta_seite1"] == ((-1,), (-1,)))

# ================================================================== Teil 2
say("\nTeil 2: (Z2)^3-Überlagerung r: E' -> P^1 (Bestand z34) und ihre Sektoren")
marks_c = [1, 1j, -1, -1j]
names = ["p1", "pi", "pm1", "pmi"]


def perm(f):
    return tuple(min(range(4), key=lambda j: abs(marks_c[j] - f(marks_c[i]))) for i in range(4))


def comp(p, q_):                                                            # erst q_, dann p
    return tuple(p[q_[i]] for i in range(4))


rot, ref = perm(lambda z: 1j * z), perm(lambda z: 1 / z)
ide = (0, 1, 2, 3)
D4 = {ide}
while True:
    new = {comp(g, h) for g in D4 | {rot, ref} for h in D4 | {rot, ref}} | D4
    if new == D4:
        break
    D4 = new
r2 = comp(rot, rot)
V4m = {ide, r2, comp(rot, ref), comp(r2, comp(rot, ref))}
check("D4 = <z -> iz, z -> 1/z> hat Ordnung 8 auf den Marken; V4 = {1, r^2, rs, r^3 s} besteht aus den drei "
      "Doppeltranspositionen (Kleingruppe)", len(D4) == 8 and len(V4m) == 4 and
      all(sum(1 for i in range(4) if g[i] == i) == 0 for g in V4m - {ide}))
elem = {0: (0, 0), 1: (1, 0), 2: (0, 1), 3: (1, 1)}                         # Marke -> E[2] (p_1 als Ursprung)
chars = {}
for dsg in (1, -1):
    for c1, c2 in itertools.product((1, -1), repeat=2):
        Sx = frozenset(i for i in range(4) if dsg * (c1 ** elem[i][0]) * (c2 ** elem[i][1]) == -1)
        chars[(dsg, (c1, c2))] = Sx
even_sets = {frozenset(s) for k in (0, 2, 4) for s in itertools.combinations(range(4), k)}
blocks = {c: (chars[(1, c)], chars[(-1, c)]) for c in itertools.product((1, -1), repeat=2)}
ok_ch = set(chars.values()) == even_sets and len(set(chars.values())) == 8 and \
    all(Sp | Sm == frozenset(range(4)) and not (Sp & Sm) for Sp, Sm in blocks.values())
mono = sorted([1 - len(S_) // 2 for S_ in chars.values()], reverse=True)
e0side = {}
for c, (Sp, Sm) in blocks.items():
    if not Sp:
        e0side[c] = split_push({})
    else:
        x_, y_ = sorted(Sp)
        e0side[c] = split_push({names[x_]: 1, names[y_]: -1})
ok_ch &= all(e0side[c] == tuple(sorted((-(len(Sp) // 2), -(len(Sm) // 2)), reverse=True)) for c, (Sp, Sm) in blocks.items())
ok_ch &= mono == [1, 0, 0, 0, 0, 0, 0, -1]
ZAHLEN["sektoren_r_O"] = {str(c): {"S+": sorted(Sp), "S-": sorted(Sm), "E0": e0side[c]} for c, (Sp, Sm) in blocks.items()}
check("Monodromie (z34) und E-Rechnung stimmen überein: Charaktere <-> gerade Markenmengen S, Linienbündel O(-|S|/2); "
      "V4-Sektoren = Zerlegungen {S, S^c}: trivial O + O(-2) = pi_* O, sonst O(-1) + O(-1) = pi_* eta_c; "
      "r_* r^* O(1) = O(1) + 6 O + O(-1)", ok_ch)
hom_o1 = sum(1 for S_ in chars.values() if 1 - len(S_) // 2 >= 1)
check("kein O(1) + O(1) in r_* L (L = r^* O(1)): Hom(O(1), r_* L) = H0(r_* L (-1)) ist eindimensional, und nur ein Charakter "
      "trägt O(1): weder eine Einbettung noch eine deckäquivariante Surjektion r_* L -> O(1) + O(1) existiert", hom_o1 == 1)


def type_of(S1, S2, k):
    return tuple(sorted((k - len(S1) // 2, k - len(S2) // 2), reverse=True))


def d4_inv(pair):
    return all({frozenset(g[i] for i in S_) for S_ in pair} == set(pair) for g in D4)


tab_k = {}
for k in (1, 2, 3):
    prs = list(itertools.combinations(sorted(even_sets, key=sorted), 2))
    bal = [p for p in prs if type_of(*p, k) == (1, 1)]
    inv = [p for p in prs if d4_inv(p)]
    tab_k[k] = {"paare": len(prs), "balanciert_grad2": len(bal), "D4_invariant": [(sorted(map(sorted, p)), type_of(*p, k))
                                                                                  for p in inv],
                "D4_inv_balanciert": [sorted(map(sorted, p)) for p in inv if type_of(*p, k) == (1, 1)]}
ZAHLEN["klassifikation_r_L^k"] = tab_k
anti = [[0, 2], [1, 3]]
check("Klassifikation der 28 deckinvarianten Rang-2-Unterbündel von r_* L^k: O(1) + O(1) nur für k = 2 (L^2 = r^* O(2) = "
      "K_{E/B}) und genau für zwei Paarcharaktere (15); D4-invariant sind nur zwei: der triviale Sektor O(2) + O und der "
      "antipodale Sektor {p_1, p_-1 | p_i, p_-i} = O(1) + O(1)",
      tab_k[1]["balanciert_grad2"] == 0 and tab_k[3]["balanciert_grad2"] == 0 and tab_k[2]["balanciert_grad2"] == 15 and
      len(tab_k[2]["D4_invariant"]) == 2 and tab_k[2]["D4_inv_balanciert"] == [anti])
h0_blocks = sorted(sum(max(t + 1, 0) for t in type_of(Sp, Sm, 2)) for Sp, Sm in blocks.values())
check("h0(L^2) = 16 = 4 + 4 + 4 + 4: jeder V4-Sektor von r_* K_{E/B} hat vier Schnitte; nur die drei nichttrivialen sind "
      "balanciert (Twistortyp), der triviale ist O(2) + O", h0_blocks == [4, 4, 4, 4])
xi = {"inf+": 1, "p1": -1}
ok_np = ell(dsc(xi, 2)) == 0 and ell(dsc(xi, 4)) == 1
np_split = {nm: split_push(dadd(xi, eta, PB1)) for nm, eta in
            {"O": {}, "pi-p1": {"pi": 1, "p1": -1}, "pm1-p1": {"pm1": 1, "p1": -1}, "pmi-p1": {"pmi": 1, "p1": -1}}.items()}
ok_np &= all(v == (0, 0) for v in np_split.values())
ZAHLEN["nicht_pullback_spinlift"] = np_split
check("relative Selbstdualität L^2 = K_{E/B} ohne die Annahme L = r^* M: drei weitere Lösungen L' = r^* O(1) (x) eta; über die "
      "Isogenie ist r_* L' = O(1) (x) Summe_c pi_*(xi eta_c) mit xi = O(inf+ - p_1) der Ordnung 4: alle vier Blöcke O + O, also "
      "r_* L' = 8 O (balanciert, kein O(1)-Summand, kein schwaches Dublett; passt zu 4 + 4 außerhalb (P,P) in z38)", ok_np)

# ================================================================== Teil 3
say("\nTeil 3: Hecke-Modifikationen von O(1) + O(-1) = O(1) (x) pi_* O_E an den Marken")
pts_m = {"p1": Q(1), "pi": Q(0, 1), "pm1": Q(-1), "pmi": Q(0, -1)}
ok_h, tab_h = True, {}
for k in range(5):
    for S in itertools.combinations(names, k):
        base = split_hecke(1, -1, [pts_m[a] for a in S], [(0, 1)] * k)
        eside = split_push(dadd(PB1, {a: 1 for a in S}))
        hb = [h0_hecke(1, -1, [pts_m[a] for a in S], [(0, 1)] * k, m) for m in range(-3, 3)]
        he = [ell(dadd(PB1, {a: 1 for a in S}, dsc(PB1, m))) for m in range(-3, 3)]
        ok_h &= base == eside and hb == he and (bal2(base) == (k == 2))
        tab_h["+".join(S) or "leer"] = base
ZAHLEN["hecke_ungerade"] = tab_h
check("Hecke längs der ungeraden (zweiten Blatt-)Geraden an S: Basisrechnung = E-Rechnung O(1) (x) pi_* O_E(S) für alle 16 "
      "Teilmengen (h0-Folgen gleich); Ergebnis O(1) + O(|S| - 1), also O(1) + O(1) genau für |S| = 2", ok_h)
lns = {"ungerade": (0, 1), "generisch": (1, 2), "gerade": (1, 0)}
tab_l = {}
for (n1, l1), (n2, l2) in itertools.product(lns.items(), repeat=2):
    tab_l[f"{n1}/{n2}"] = split_hecke(1, -1, [Q(1), Q(-1)], [l1, l2])
ZAHLEN["hecke_geraden_zwei_marken"] = tab_l
check("an zwei Marken hängt das Ergebnis von den Geraden ab: O(1) + O(1), solange keine die gerade Richtung (Fasern des "
      "destabilisierenden O(1)) ist; sonst O(2) + O bzw. O(3) + O(-1): auch die kanonische gerade Wahl springt",
      all(bal2(v) == ("gerade" not in key.split("/")) for key, v in tab_l.items()) and tab_l["gerade/ungerade"] == (2, 0)
      and tab_l["gerade/gerade"] == (3, -1))
two = [frozenset(s) for s in itertools.combinations(range(4), 2)]
fixed_V4 = [s for s in two if all(frozenset(g[i] for i in s) == s for g in V4m)]
fixed_D4 = [s for s in two if all(frozenset(g[i] for i in s) == s for g in D4)]
parts_ = {frozenset({s, frozenset(range(4)) - s}) for s in two}
pf_D4 = [p for p in parts_ if all({frozenset(g[i] for i in s) for s in p} == set(p) for g in D4)]
pf_V4 = [p for p in parts_ if all({frozenset(g[i] for i in s) for s in p} == set(p) for g in V4m)]
ZAHLEN["hecke_symmetrie"] = {"2-Mengen fest unter V4": len(fixed_V4), "unter D4": len(fixed_D4),
                             "Zerlegungen fest unter D4": [sorted(map(sorted, p)) for p in pf_D4], "unter V4": len(pf_V4)}
check("keine der sechs Zweiermengen ist unter V4 oder D4 fest (die Modifikation selbst ist nie symmetrisch); fest sind nur "
      "Zerlegungen: alle drei unter V4, unter D4 genau die antipodale. Da S und S^c dasselbe Linienbündel auf E geben, ist die "
      "Isomorphieklasse O(1) + O(1) kanonisch, die Einbettung eine binäre Wahl",
      len(fixed_V4) == 0 and len(fixed_D4) == 0 and len(pf_V4) == 3 and
      [sorted(map(sorted, p)) for p in pf_D4] == [anti])

# ================================================================== Teil 4
say("\nTeil 4: Partner, Dublett und Lorentzkovarianz")
jet = {}
for k in range(5):
    G = [[{-k: Q(1)}, {}], [{1 - k: Q(k)} if k else {}, {2 - k: Q(-1)}]]                # (s, ds) -> (s_inf, ds_inf/dw)
    jet[k] = split_trans(G)
ZAHLEN["jets"] = jet
check("Jetbündel aus Übergangsfunktionen: J^1(O) = O + O(-2), J^1(O(k)) = O(k-1) + O(k-1) für k >= 1; also J^1(O(2)) = "
      "J^1(T_C) = J^1(K(D)) = O(1) + O(1) = O(1) (x) H0(O(1)) (Kandidat 4, kanonisch, nur aus P1)",
      jet[0] == (0, -2) and all(jet[k] == (k - 1, k - 1) for k in range(1, 5)) and jet[2] == (1, 1))
eul = {}
for c in (0, 1, 3):
    G = [[{1: Q(1)}, {0: Q(c)} if c else {}], [{}, {-1: Q(1)}]]
    Gt = [[{0: Q(1)}, {-1: Q(c)} if c else {}], [{}, {-2: Q(1)}]]                     # O(1) (x) Erweiterung
    eul[c] = (split_trans(G), split_trans(Gt))
ZAHLEN["euler_kanal"] = {str(c): v for c, v in eul.items()}
check("Euler-Kanal (Rang-1-Kanal z35/z37 = Erweiterungsklasse in H1(O(-2))): ohne Kanal O(1) + O(-1) (z61-Sektor), mit Kanal "
      "O + O = V (x) O (volles Dublett trivial); O(1) (x) (Dublett) = O(1) + O(1), ohne Kanal O(2) + O",
      eul[0] == ((1, -1), (2, 0)) and eul[1] == ((0, 0), (1, 1)) and eul[3] == ((0, 0), (1, 1)))
th_ = np.linspace(0, 2 * np.pi, 4001)
zc = np.exp(1j * th_)


def winding(e0, einf):
    f = np.array([np.vdot(einf(z), e0(z)) / np.vdot(einf(z), einf(z)) for z in zc])
    par = max(abs(e0(z)[0] * einf(z)[1] - e0(z)[1] * einf(z)[0]) for z in zc[::50])
    return int(round((np.unwrap(np.angle(f))[-1] - np.angle(f[0])) / (2 * np.pi))), par


frames = {"tautologisch l": (lambda z: np.array([1, z]), lambda z: np.array([1 / z, 1])),
          "konjugiert l-bar (glatt)": (lambda z: np.array([1, np.conj(z)]), lambda z: np.array([1 / np.conj(z), 1])),
          "sigma^* l-bar (holomorph)": (lambda z: np.array([-z, 1]), lambda z: np.array([-1, 1 / z])),
          "V-bar (x) O": (lambda z: np.array([1, 0]), lambda z: np.array([1, 0]))}
wd = {k: winding(*v) for k, v in frames.items()}
degs = {k: -w for k, (w, _) in wd.items()}
ok_sig = True
for _ in range(50):
    z = complex(*rng.normal(size=2))
    sz = -1 / np.conj(z)
    v = np.conj(np.array([1, sz]))
    ok_sig &= abs(v[0] * 1 - v[1] * (-z)) < 1e-12
ZAHLEN["grade_partner"] = degs
o1Vbar = split_trans([[{-1: Q(1)}, {}], [{}, {-1: Q(1)}]])
check("Grad des Partners (Windungszahlen der Rahmen auf |z| = 1): tautologische Linie -1, ihre Konjugierte +1 (nur glatt), "
      "über die Antipodenabbildung holomorph gemacht wieder -1 (Rahmen (-z, 1)), das konjugierte Dublett V-bar (x) O hat Grad 0; "
      "O(1) (x) V-bar = O(1) + O(1): der Partner liefert Vielfachheit 2, keine zweite Windung",
      degs == {"tautologisch l": -1, "konjugiert l-bar (glatt)": 1, "sigma^* l-bar (holomorph)": -1, "V-bar (x) O": 0}
      and all(p_ < 1e-9 for _, p_ in wd.values()) and ok_sig and o1Vbar == (1, 1))


def rsl2():
    g = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    return g / np.sqrt(np.linalg.det(g))


rho = {"C^2 (trivial)": lambda g: np.eye(2), "V": lambda g: g, "V-bar^*": lambda g: np.linalg.inv(g.conj().T)}
res = {}
gs = [rsl2() for _ in range(200)]
for nm, rf in rho.items():
    herm_ok, char_real = True, True
    for g in gs:
        A = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        x = A + A.conj().T
        y = rf(g) @ x @ np.linalg.inv(g)
        herm_ok &= np.allclose(y, y.conj().T, atol=1e-9)
        char_real &= abs((np.trace(rf(g)) * np.trace(np.linalg.inv(g))).imag) < 1e-9
    Ms = [np.kron(np.linalg.inv(g).T, rf(g)) - np.eye(4) for g in gs[:4]]
    inv_dim = 4 - np.linalg.matrix_rank(np.vstack(Ms), tol=1e-9)
    Cs = [np.kron(np.eye(2), rf(g)) - np.kron(rf(g).T, np.eye(2)) for g in gs[:4]]
    comm_dim = 4 - np.linalg.matrix_rank(np.vstack(Cs), tol=1e-9)
    res[nm] = {"Herm(2) erhalten": bool(herm_ok), "Charakter reell": bool(char_real), "invariante Vektoren": int(inv_dim),
               "Kommutant": int(comm_dim)}
ZAHLEN["lorentzkriterium"] = res
check("Lorentzkriterium für N = O(1) (x) W, H0 = Hom(V, W), x -> rho_W(g) x g^-1, 200 zufällige g in SL(2,C): nur W = V-bar^* "
      "erhält die reelle Form Herm(2) und hat reellen Charakter |tr g|^2 (Lorentzvektor (1/2,1/2)); W = V gibt 3 + 1 (invarianter "
      "Vektor, nur Drehungen), W = C^2 zwei chirale Dubletts",
      res["V-bar^*"]["Herm(2) erhalten"] and res["V-bar^*"]["Charakter reell"] and res["V-bar^*"]["invariante Vektoren"] == 0
      and not res["V"]["Herm(2) erhalten"] and res["V"]["invariante Vektoren"] == 1 and not res["V"]["Charakter reell"]
      and not res["C^2 (trivial)"]["Herm(2) erhalten"] and not res["C^2 (trivial)"]["Charakter reell"])
w8 = np.exp(1j * np.pi / 4)
rt = np.diag([1 / w8, w8])                                                  # Hebung von z -> iz (z = v2/v1)
st = np.array([[0, 1j], [1j, 0]])                                           # Hebung von z -> 1/z


def closure(gens):
    G = [np.eye(2, dtype=complex)]
    changed = True
    while changed:
        changed = False
        for a in list(G):
            for b in gens:
                c = a @ b
                if not any(np.allclose(c, d) for d in G):
                    G.append(c)
                    changed = True
    return G


def herm_inv_dim(G):
    rows = []
    for g in G:
        M = []
        basis = [np.array([[1, 0], [0, 0]]), np.array([[0, 0], [0, 1]]), np.array([[0, 1], [1, 0]]),
                 np.array([[0, -1j], [1j, 0]])]
        for B in basis:
            D_ = g.conj().T @ B @ g - B
            M.append(np.concatenate([D_.real.ravel(), D_.imag.ravel()]))
        rows.append(np.array(M).T)
    return int(4 - np.linalg.matrix_rank(np.vstack(rows), tol=1e-9))


Q8 = closure([rt @ rt, rt @ st])
B16 = closure([rt, st])
hardy = np.array([[np.mean(np.conj(zc[:-1] ** j) * zc[:-1] ** k) for k in range(2)] for j in range(2)])
dims = {"Q8 (V4 gehoben)": herm_inv_dim(Q8), "2D4 (D4 gehoben)": herm_inv_dim(B16), "<r> allein": herm_inv_dim(closure([rt]))}
inv_ok = all(np.allclose(g.conj().T @ hardy @ g, hardy) for g in B16)
ZAHLEN["schur"] = {**dims, "|Q8|": len(Q8), "|2D4|": len(B16)}
check("Schur: die gehobene Markensymmetrie (Q8 aus V4, Ordnung 8; binär-dihedral aus D4, Ordnung 16) wirkt irreduzibel auf V, "
      "der Raum invarianter Hermitescher Formen ist eindimensional; es ist die Hardy-Form der Kreislinie durch die Marken "
      "(Gram = 1 auf {1, z}); Kontrolle: die Drehung allein lässt zwei Formen zu",
      len(Q8) == 8 and len(B16) == 16 and dims["Q8 (V4 gehoben)"] == 1 and dims["2D4 (D4 gehoben)"] == 1
      and dims["<r> allein"] == 2 and np.allclose(hardy, np.eye(2), atol=1e-12) and inv_ok)
dk = np.diag([1, -1])
noncomm = all(not np.allclose(dk @ rho["V-bar^*"](g), rho["V-bar^*"](g) @ dk) for g in gs[:50])
check("Blattblindheit: Lorentzkovarianz mit dem Deck als innerer, mit SL(2,C) vertauschender Symmetrie verlangt nach Schur "
      "(Kommutant von V-bar^* ist "
      f"{res['V-bar^*']['Kommutant']}-dimensional) skalare Deckwirkung auf N; in r_* L^k ist jeder der acht Charaktere einfach "
      "(Teil 2) und in jedem geraden Sektor wirkt das Deck als diag(1, -1) (Teil 1), das mit keinem Boost-Element vertauscht: "
      "kein Überlagerungssektor ist ein lorentzkovariantes Normalenbündel",
      res["V-bar^*"]["Kommutant"] == 1 and res["C^2 (trivial)"]["Kommutant"] == 4 and noncomm and len(set(chars.values())) == 8
      and parts["eta_antipodal"] == ((-1,), (-1,)))
befund("Für T' braucht es weder Zustand noch Hermitesche Form: Hom(V, V-bar^*) sind die Sesquilinearformen auf V, ihre reelle "
       "Form die Hermiteschen Formen, der Kegel die positiven; das alles ist kanonisch. Die Markensymmetrie erzwingt nur die "
       "Hermitesche Form h (Hardy-Form), also die Identifikation mit End(V) aus z62 und ein Ruhesystem der Marken. Welchen "
       "Punkt des Kegels der Dublettzustand (Lücke A) wählt, legt sie nicht fest, sobald der Zustand sie bricht", True)

# ================================================================== Teil 5
say("\nTeil 5: Bilanz -- welche Konstruktion gibt O(1) + O(1)?")
tafel = [
    ("z61: O(1) (x) pi_* O_E = [r_* L]_trivial", str(tab_pb[1]), "kanonisch", "nein (Grad 0, Lücke 2)"),
    ("O(1) (x) [r_* L]_trivial = [r_* K_{E/B}]_trivial", str(tab_t["O (ungerade, P,P)"]["T_C (x) pi_* theta"]), "kanonisch",
     "nein (Nullmodensektor (P,P), schwaches Dublett)"),
    ("(3) T_C (x) pi_* theta, theta gerade (Szegő)", "(1, 1)", "Klasse kanonisch, D4 hält genau eta_* fest", "holomorph ja"),
    ("(5) [r_* K_{E/B}]_antipodal = Farbe_c* + Familie_c*, (x) O(1)", "(1, 1)",
     "einziger D4-invarianter balancierter Sektor", "holomorph ja, nicht lorentzkovariant (Deck nicht skalar)"),
    ("(5) Unter-/Quotientenbündel von r_* L, E+, E-", "-", "-", "nein (kein O(1)+O(1), keine äquivariante Surjektion)"),
    ("(2) Hecke an zwei Marken (ungerade Gerade)", "(1, 1)", "nur Isomorphieklasse (Zerlegung)", "holomorph ja"),
    ("(2) Hecke an 0, 1, 3, 4 Marken oder längs gerader Geraden", "O(1)+O(|S|-1) bzw. springend", "-", "nein"),
    ("(4) O(1) (x) H0(O(1)) = J^1(O(2)) = O(1) (x) [Euler-Dublett]", "(1, 1)", "kanonisch aus P1",
     "holomorph ja, aber 3 + 1 (keine Boosts)"),
    ("(1) O(1) (x) V-bar^* (Purifikationspartner)", "(1, 1)", "kanonisch, keine Wahl",
     "einziger lorentzkovarianter Kandidat; als Normalenbündel Prämisse T'"),
]
for row in tafel:
    say("   " + " | ".join(row))
ZAHLEN["tafel"] = [dict(zip(["konstruktion", "spaltung", "kanonisch", "ergebnis"], r_)) for r_ in tafel]
say("")
status("Spaltungstyp", "E", "O(1) + O(1) ist kein Hindernis: es entsteht kanonisch in den Szegő-Sektoren (gerade Spin-"
       "strukturen), aus J^1(O(2)) = O(1) (x) Euler-Dublett und aus O(1) (x) V-bar^*")
status("z61-Klasse", "E", "No-go: Blätter mit zurückgezogenem Bündel geben stets O(k) + O(k-2); der Pullback ist der einzige "
       "Sprungpunkt jeder Gradklasse, und genau er trägt das schwache Dublett (5 + 3)")
status("Kanonischer Überlagerungssektor", "E", "in r_* K_{E/B} ist der antipodale Sektor der einzige D4-invariante "
       "balancierte; seine Faser ist intern (Farbe_c*, Familie_c*)")
status("Lorentzkovarianz", "E", "für N = O(1) (x) W ist H0(N) genau dann der Lorentzvektor, wenn W = V-bar^* (Charakter "
       "|tr g|^2); SL(2)-äquivariante holomorphe Konstruktionen geben 3 + 1 oder chirale Dubletts; jeder Überlagerungssektor "
       "scheitert an der Blattblindheit (Schur), sofern das Deck als innere Symmetrie mit der Lorentzwirkung vertauscht")
status("Hermitesche Form", "E", "durch die Markensymmetrie bis auf Skala erzwungen (Q8 bzw. 2D4 irreduzibel): Hardy-Form")
status("Lesart T", "C", "folgt vollständig unter der einen Prämisse T': die Naht ist P(V) in P(V + V-bar^*), also "
       "N = O(1) (x) V-bar^*; T' ist kanonisch (keine Wahl, einziger lorentzkovarianter Kandidat der Form O(1) (x) W), aber "
       "nicht aus den Daten abgeleitet")
status("Faser S_* = (Farbe_c*, Familie_c*) als zweiter Spinor", "L", "nur holomorphe Koinzidenz; als Raumzeitspinor "
       "ausgeschlossen (Deck nicht skalar)")
status("Nahtlagen als Antworten", "O", "dass die Naht überhaupt deformiert (eine Einbettung in eine 3-Mannigfaltigkeit), ist "
       "aus keiner Antwort des Originals rekonstruiert")
befund("Das zweite O(1) ist kein Gradproblem und kein Problem der Blätter: Es muss blattblind sein, und der einzige "
       "blattblinde, lorentzkovariante Träger ist das konjugierte Dublett V-bar^* des Partners", True)
befund("Spannung im Bestand: Der Sektor, der 5 + 3 und das schwache Dublett trägt (Pullback, (P,P)), ist genau der, in dem T "
       "scheitert; die Sektoren mit Szegő-Kern sind balanciert, haben aber kein Dublett (4 + 4)", True)

urteil = ("Das zweite O(1) der Lesart T ist kein Grad- oder Spaltungsproblem: O(1) + O(1) entsteht kanonisch aus den TFPT-Daten, "
          "nämlich als T_C (x) pi_* theta für jede gerade Spinstruktur von E (Szegő-Sektor), als einziger D4-invarianter "
          "balancierter Sektor von r_* K_{E/B} (antipodal: Farbe_c* + Familie_c*), als J^1(O(2)) = O(1) (x) Euler-Dublett und "
          "als O(1) (x) V-bar^*. Exakt ausgeschlossen ist die z61-Klasse: Blätter mit zurückgezogenem Bündel geben stets "
          "O(k) + O(k-2), weil der Pullback der einzige Sprungpunkt ist; derselbe Sprung trägt das schwache Dublett (5 + 3). "
          "Entscheidend ist die Lorentzkovarianz: Unter N = O(1) (x) W ist H0(N) nur für W = V-bar^* der Lorentzvektor; "
          "holomorphe Konstruktionen geben 3 + 1 ohne Boosts, und kein Überlagerungssektor ist blattblind (Deck nicht skalar, "
          "Schur, mit dem Deck als innerer Symmetrie). Die "
          "Hermitesche Form ist durch die Markensymmetrie erzwungen. Damit bleibt T bedingt (C) mit genau einer kanonischen "
          "Prämisse T': Die Naht ist P(V) in P(V + V-bar^*), ihre Normalenrichtungen sind das konjugierte Dublett des Partners. "
          "T' folgt nicht aus den Daten; ausgeschlossen ist es auch nicht.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z64_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "status": STATUS, "zahlen": ZAHLEN, "urteil": urteil}, indent=1,
    ensure_ascii=False, default=str) + "\n")
