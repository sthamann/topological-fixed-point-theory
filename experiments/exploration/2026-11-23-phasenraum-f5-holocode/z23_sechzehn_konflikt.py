"""Z23: Die zwei Sechzehnen. Naht-16 (Vektor von SO(16)) gegen Generations-16 (Spinor von Spin(10)).

Frage. Gibt es eine konsistente Zuordnung, in der die 16 Majorana-Weyl-Moden der Naht und die drei
leichten 16 der Generationen aus einer Struktur kommen, oder eine Obstruktion?

Gerechnet wird (alles aus Wurzel- und Gewichtssystemen, Freudenthal, keine zitierten Tabellen):
 1. Verzweigungen: E8 > SO(16) > SO(10) x SO(6) > SO(10) x SU(3) x U(1)_X; E8 > E6 x SU(3);
    E8 > SU(5) x SU(5); E8 > E7 x SU(2); SM-Inhalt der Generation.
 2. Gitter und Verklebung D5 + A3 (R102): wo Naht-Eltern und Generation im E8-Gitter liegen; Zentrum.
 3. Automorphismen, Trialität, Realität; der J-Satz für ein reelles Multiplett (Z18-Lesart).
 4. Ströme zu Fermionen: Clifford-Wörterbuch, Statistik (h mod 1), J auf E8, E8-Klammer auf 128 x 128,
    Chiralität mit Fluss und im Konstantmoden-Orbifold T^2/Z4.
 5. Getrennte Zuordnung: 4D-Felder und Anomalien (SO(10), SM, Z16 für X, gemeinsames Spin-Z4).
 6. Urteil.
"""
import itertools
import json
from collections import Counter
from fractions import Fraction as Fr
from math import comb, gcd
from pathlib import Path

import numpy as np
import sympy as sp

HERE = Path(__file__).parent
LINES, CHECKS, BEFUNDE, DATA = [], {}, {}, {}
RNG = np.random.default_rng(20261123)


def say(s=""):
    print(s)
    LINES.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    return "ja" if cond else "NEIN"


def befund(name, val):
    BEFUNDE[name] = bool(val)
    return val


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def vec(*xs):
    return tuple(Fr(x) for x in xs)


# ---------------------------------------------------------------------------------------------
# Lie-Algebra-Werkzeug: Cartan-Matrix, positive Wurzeln, Freudenthal, Zerlegung
# ---------------------------------------------------------------------------------------------
class Lie:
    def __init__(self, A):
        self.A = [list(map(int, r)) for r in A]
        self.n = n = len(A)
        Ainv = sp.Matrix(self.A).inv()
        D = 1
        for x in Ainv:
            D = sp.ilcm(D, sp.fraction(sp.nsimplify(x))[1])
        self.D = int(D)
        self.F = [[int(Ainv[i, j] * self.D) for j in range(n)] for i in range(n)]
        self.pos = self._posroots()
        self.Fpos = [[sum(self.F[i][j] * a[j] for j in range(n)) for i in range(n)] for a in self.pos]
        self.cache = {}

    def ip(self, x, y):
        n = self.n
        return sum(x[i] * self.F[i][j] * y[j] for i in range(n) for j in range(n))

    def height(self, mu):
        return sum(self.F[i][j] * mu[j] for i in range(self.n) for j in range(self.n))

    def _posroots(self):
        n, A = self.n, self.A
        simple = [tuple(1 if j == i else 0 for j in range(n)) for i in range(n)]
        roots, layer = set(simple), list(simple)
        while layer:
            new = []
            for b in layer:
                lab = [sum(b[k] * A[k][j] for k in range(n)) for j in range(n)]
                for j in range(n):
                    r = 0
                    while True:
                        c = list(b); c[j] -= r + 1
                        if tuple(c) in roots:
                            r += 1
                        else:
                            break
                    if r - lab[j] > 0:
                        c = list(b); c[j] += 1; c = tuple(c)
                        if c not in roots:
                            roots.add(c); new.append(c)
            layer = new
        return [tuple(sum(b[k] * A[k][j] for k in range(n)) for j in range(n)) for b in roots]

    def char(self, lam):
        lam = tuple(lam)
        if lam in self.cache:
            return self.cache[lam]
        n, A = self.n, self.A
        lr = tuple(l + 1 for l in lam)
        N = self.ip(lr, lr)
        twoD = 2 * self.D
        mult, layer = {lam: 1}, [lam]
        while layer:
            cand = set()
            for mu in layer:
                for i in range(n):
                    cand.add(tuple(mu[j] - A[i][j] for j in range(n)))
            new = []
            for mu in cand:
                if mu in mult:
                    continue
                mr = tuple(m + 1 for m in mu)
                den = N - self.ip(mr, mr)
                if den <= 0:
                    continue
                s = 0
                for a, fa in zip(self.pos, self.Fpos):
                    base = sum(mu[i] * fa[i] for i in range(n))
                    k, nu = 1, tuple(mu[i] + a[i] for i in range(n))
                    while nu in mult:
                        s += mult[nu] * (base + k * twoD)
                        k += 1
                        nu = tuple(nu[i] + a[i] for i in range(n))
                assert (2 * s) % den == 0
                m = 2 * s // den
                if m > 0:
                    mult[mu] = m
                    new.append(mu)
            layer = new
        self.cache[lam] = mult
        return mult

    def dim(self, lam):
        return sum(self.char(lam).values())

    def decompose(self, counter):
        rem = {k: v for k, v in counter.items() if v}
        res = {}
        while rem:
            cands = [mu for mu, v in rem.items() if v > 0 and all(x >= 0 for x in mu)]
            mu = max(cands, key=self.height)
            c = rem[mu]
            res[mu] = res.get(mu, 0) + c
            for w, m in self.char(mu).items():
                rem[w] = rem.get(w, 0) - c * m
                if rem[w] == 0:
                    del rem[w]
            assert all(v > 0 for v in rem.values()), "negative Restmultiplizität"
        return res


def simple_roots(roots):
    g = [1.0 / (np.pi ** (k + 1)) + 1e-3 * k for k in range(len(roots[0]))]
    pos = [r for r in roots if sum(float(x) * y for x, y in zip(r, g)) > 0]
    pset = set(pos)
    simple = []
    for r in pos:
        if not any(tuple(a - b for a, b in zip(r, s)) in pset for s in pos if s != r):
            simple.append(r)
    return simple


def cartan(simple):
    return [[int(dot(a, b)) for b in simple] for a in simple]


def labels(w, simple):
    out = []
    for s in simple:
        x = dot(w, s)
        assert x.denominator == 1
        out.append(int(x))
    return tuple(out)


def components(A):
    n = len(A)
    seen, comps = set(), []
    for i in range(n):
        if i in seen:
            continue
        st, comp = [i], []
        seen.add(i)
        while st:
            k = st.pop(); comp.append(k)
            for j in range(n):
                if j not in seen and A[k][j] != 0:
                    seen.add(j); st.append(j)
        comps.append(sorted(comp))
    return comps


class Sub:
    """Halbeinfache Unteralgebra (gegeben durch einfache Wurzeln in R^8) plus U(1)-Richtungen."""

    def __init__(self, simple, u1=(), names=None):
        self.simple = simple
        self.A = cartan(simple)
        self.lie = Lie(self.A)
        self.comps = components(self.A)
        self.u1 = u1
        self.names = names or {}

    def key(self, w):
        return labels(w, self.simple), tuple(dot(w, u) for u in self.u1)

    def branch(self, weights):
        by_q = {}
        for w in weights:
            lab_, q = self.key(w)
            by_q.setdefault(q, Counter())[lab_] += 1
        out = Counter()
        for q, cnt in by_q.items():
            for hw, m in self.lie.decompose(cnt).items():
                out[(hw, q)] += m
        return out

    def factor_dims(self, hw):
        res = []
        for comp in self.comps:
            sub = [[self.A[i][j] for j in comp] for i in comp]
            L = Lie(sub)
            res.append(L.dim(tuple(hw[i] for i in comp)))
        return res

    def name(self, hw, q=()):
        parts = []
        for ci, comp in enumerate(self.comps):
            part = tuple(hw[i] for i in comp)
            nm = self.names.get((ci, part))
            if nm is None:
                sub = [[self.A[i][j] for j in comp] for i in comp]
                nm = str(Lie(sub).dim(part))
            parts.append(nm)
        s = "(" + ",".join(parts) + ")"
        if q:
            s += "_" + ",".join(str(x) for x in q)
        return s

    def table(self, weights):
        br = self.branch(weights)
        items = sorted(br.items(), key=lambda kv: (-self.lie.dim(kv[0][0]), str(kv[0])))
        return [(self.name(hw, q), m, self.lie.dim(hw)) for (hw, q), m in items], br


def comp_hw(sub, ci, weights):
    """Höchstgewicht (Labels des Faktors ci) einer irreduziblen Gewichtsmenge."""
    comp = sub.comps[ci]
    subA = [[sub.A[i][j] for j in comp] for i in comp]
    L = Lie(subA)
    cnt = Counter(tuple(labels(w, sub.simple)[i] for i in comp) for w in weights)
    dec = L.decompose(cnt)
    assert len(dec) == 1 and list(dec.values())[0] == 1, dec
    return list(dec.keys())[0]


def fmt_table(rows):
    return " + ".join(f"{m}x{n}" if m > 1 else n for n, m, d in rows)


# ---------------------------------------------------------------------------------------------
# 0. E8 explizit
# ---------------------------------------------------------------------------------------------
def e8_roots():
    R = []
    for i in range(8):
        for j in range(i + 1, 8):
            for si in (1, -1):
                for sj in (1, -1):
                    v = [Fr(0)] * 8; v[i] = Fr(si); v[j] = Fr(sj); R.append(tuple(v))
    for signs in itertools.product((1, -1), repeat=8):
        if signs.count(-1) % 2 == 0:
            R.append(tuple(Fr(s, 2) for s in signs))
    return R


ROOTS = e8_roots()
ZERO = tuple(Fr(0) for _ in range(8))
ADJ = ROOTS + [ZERO] * 8
D8R = [r for r in ROOTS if all(x.denominator == 1 for x in r)]
SPIN128 = [r for r in ROOTS if r[0].denominator == 2]
VEC16 = [tuple(Fr(s) if k == i else Fr(0) for k in range(8)) for i in range(8) for s in (1, -1)]


def in_e8_lattice(v):
    ints = all(x.denominator == 1 for x in v)
    halfs = all((2 * x).denominator == 1 and x.denominator == 2 for x in v)
    return (ints or halfs) and sum(v) % 2 == 0


def odd_minus5(w):
    """SO(10)-Teil (Koordinaten 1..5) in der Klasse der Generation: halbzahlig, ungerade viele Minus."""
    return all(x.denominator == 2 for x in w[:5]) and sum(1 for x in w[:5] if x < 0) % 2 == 1


def even_minus5(w):
    return all(x.denominator == 2 for x in w[:5]) and sum(1 for x in w[:5] if x < 0) % 2 == 0


say("0. E8-Wurzelsystem in R^8")
say("")
E8S = simple_roots(ROOTS)
E8L = Lie(cartan(E8S))
detE8 = sp.Matrix(cartan(E8S)).det()
adj_dec = E8L.decompose(Counter(labels(w, E8S) for w in ADJ))
say(f"   240 Wurzeln, 8 einfache; det Cartan = {detE8}; 248 als E8-Irrep: {list(adj_dec.values())} mal "
    f"Dimension {E8L.dim(list(adj_dec)[0])}")
# Diagrammautomorphismen
Ae8 = cartan(E8S)
nauto = sum(1 for perm in itertools.permutations(range(8))
            if all(Ae8[perm[i]][perm[j]] == Ae8[i][j] for i in range(8) for j in range(8)))
say(f"   Automorphismen des E8-Dynkindiagramms: {nauto} (keine äußeren Automorphismen)")
check("E8: det Cartan = 1 (Gewichtsgitter = Wurzelgitter, Zentrum trivial), Diagramm ohne Symmetrie, 248 irreduzibel",
      detE8 == 1 and nauto == 1 and len(adj_dec) == 1 and E8L.dim(list(adj_dec)[0]) == 248)

# ---------------------------------------------------------------------------------------------
# 1. Verzweigungen
# ---------------------------------------------------------------------------------------------
say("\n1. Verzweigungen (Freudenthal auf den Gewichten)\n")

# 1a: E8 > D8
D8S = simple_roots(D8R)
d8 = Sub(D8S)
br = d8.branch(ADJ)
d8_rows = sorted([(d8.lie.dim(hw), m) for (hw, q), m in br.items()])
say(f"   1a. E8 > SO(16): 248 = " + " + ".join(str(d) for d, m in d8_rows))
check("E8 > SO(16): 248 = 120 + 128", d8_rows == [(120, 1), (128, 1)])
vec_in_248 = any(d8.lie.dim(hw) == 16 for (hw, q) in br)
say(f"       Kommt der Vektor 16 von SO(16) in der 248 vor? {vec_in_248}")

# 1b: D8 > D5 x D3
D5R = [r for r in D8R if all(x == 0 for x in r[5:])]
D3R = [r for r in D8R if all(x == 0 for x in r[:5])]
D5S, D3S = simple_roots(D5R), simple_roots(D3R)
# Namen festlegen
s5 = Sub(D5S)
c16 = [w[:5] + (Fr(0),) * 3 for w in SPIN128 if odd_minus5(w)]
c16b = [w[:5] + (Fr(0),) * 3 for w in SPIN128 if even_minus5(w)]
hw16 = comp_hw(s5, 0, list(set(c16)))
hw16b = comp_hw(s5, 0, list(set(c16b)))
hw10 = comp_hw(s5, 0, [v for v in VEC16 if any(v[:5])])
# Partner der Generation in SU(4): Klasse der Koordinaten 6..8 bei ungerader Minuszahl in 1..5
four = sorted(set(w[5:] for w in SPIN128 if odd_minus5(w)))
fourb = sorted(set(w[5:] for w in SPIN128 if even_minus5(w)))
s3 = Sub(D3S)
hw4 = comp_hw(s3, 0, [(Fr(0),) * 5 + f for f in four])
hw4b = comp_hw(s3, 0, [(Fr(0),) * 5 + f for f in fourb])
hw6 = comp_hw(s3, 0, [v for v in VEC16 if any(v[5:])])
NAMES_D5 = {hw16: "16", hw16b: "16b", hw10: "10", (0,) * 5: "1"}
NAMES_D3 = {hw4: "4", hw4b: "4b", hw6: "6", (0,) * 3: "1"}
s53 = Sub(D5S + D3S)
assert s53.comps == [[0, 1, 2, 3, 4], [5, 6, 7]]
for hw, nm in list(NAMES_D5.items()):
    s53.names[(0, hw)] = nm
for hw, nm in list(NAMES_D3.items()):
    s53.names[(1, hw)] = nm
adjD5 = comp_hw(s5, 0, D5R + [ZERO] * 5)
adjD3 = comp_hw(s3, 0, D3R + [ZERO] * 3)
s53.names[(0, adjD5)] = "45"
s53.names[(1, adjD3)] = "15"
t16, _ = s53.table(VEC16)
t120, _ = s53.table(D8R + [ZERO] * 8)
t128, _ = s53.table(SPIN128)
say(f"   1b. SO(16) > SO(10) x SO(6):  16 = {fmt_table(t16)}")
say(f"                                120 = {fmt_table(t120)}")
say(f"                                128 = {fmt_table(t128)}")
ok1b = (sorted(n for n, m, d in t16) == ["(1,6)", "(10,1)"]
        and sorted(n for n, m, d in t120) == ["(1,15)", "(10,6)", "(45,1)"]
        and sorted(n for n, m, d in t128) == ["(16,4)", "(16b,4b)"])
check("SO(16) > SO(10) x SO(6): 16 = (10,1)+(1,6), 120 = (45,1)+(1,15)+(10,6), 128 = (16,4)+(16b,4b)", ok1b)

# 1c: SU(4) > SU(3) x U(1)_X'
A2S = [vec(0, 0, 0, 0, 0, 1, -1, 0), vec(0, 0, 0, 0, 0, 0, 1, -1)]
uX = vec(0, 0, 0, 0, 0, 2, 2, 2)   # X' = 2(x6+x7+x8): 4 = 3_1 + 1_-3
s532 = Sub(D5S + A2S, u1=(uX,))
for hw, nm in NAMES_D5.items():
    s532.names[(0, hw)] = nm
s532.names[(0, adjD5)] = "45"
trip = [(Fr(0),) * 5 + f for f in four if sum(f) == Fr(1, 2)]
hw3 = comp_hw(Sub(A2S), 0, trip)
hw3b = tuple(reversed(hw3))
s532.names[(1, hw3)] = "3"
s532.names[(1, hw3b)] = "3b"
s532.names[(1, (1, 1))] = "8"
s532.names[(1, (0, 0))] = "1"
t248c, br248c = s532.table(ADJ)
say(f"   1c. SO(10) x SU(3) x U(1)_X' (SU(4) > SU(3) x U(1), 4 = 3_1 + 1_-3):")
say(f"       248 = {fmt_table(t248c)}")
names_c = Counter()
for n, m, d in t248c:
    names_c[n] += m
fam_ok = names_c.get("(16,3)_1", 0) == 1 and names_c.get("(16,1)_-3", 0) == 1
check("(16,4) = (16,3)_1 + (16,1)_-3; 248 unter SO(10) x SU(3) x U(1) vollständig", fam_ok and
      sum(m * d for n, m, d in t248c) == 248)

# 1d: E6 x SU(3)
E6R = [r for r in ROOTS if r[5] == r[6] == r[7]]
E6S = simple_roots(E6R)
A2R = [r for r in ROOTS if all(x == 0 for x in r[:5]) and sum(r[5:]) == 0]
s63 = Sub(E6S + A2S)
assert len(E6R) == 72 and len(A2R) == 6
e6 = Sub(E6S)
# 27: E6-Irrep der Gewichte von 248 mit A2-Klasse wie die 3 (Labels der 3 -> Klasse)
def a2class(w):
    a = labels(w, A2S)
    return (a[0] + 2 * a[1]) % 3
cl3 = a2class(trip[0])
w27x3 = [w for w in ADJ if a2class(w) == cl3]
cnt27 = Counter(labels(w, E6S) for w in w27x3)
dec27 = e6.lie.decompose(cnt27)
hw27 = list(dec27.keys())[0]
hw27b = list(e6.lie.decompose(Counter(labels(w, E6S) for w in ADJ if a2class(w) == (3 - cl3) % 3)).keys())[0]
hw78 = comp_hw(e6, 0, E6R + [ZERO] * 6)
ci_e6 = 0 if len(s63.comps[0]) == 6 else 1
ci_a2 = 1 - ci_e6
for hw, nm in ((hw27, "27"), (hw27b, "27b"), (hw78, "78"), ((0,) * 6, "1")):
    s63.names[(ci_e6, hw)] = nm
for hw, nm in ((hw3, "3"), (hw3b, "3b"), ((1, 1), "8"), ((0, 0), "1")):
    s63.names[(ci_a2, hw)] = nm
t63, _ = s63.table(ADJ)
say(f"   1d. E8 > E6 x SU(3): 248 = {fmt_table(t63)}")
ok63 = sorted(n for n, m, d in t63) == sorted(["(78,1)", "(1,8)", "(27,3)", "(27b,3b)"]) or \
       sorted(n for n, m, d in t63) == sorted(["(1,78)", "(8,1)", "(3,27)", "(3b,27b)"])
# E6 > SO(10) x U(1)
s5u = Sub(D5S, u1=(vec(0, 0, 0, 0, 0, 2, 2, 2),))
for hw, nm in NAMES_D5.items():
    s5u.names[(0, hw)] = nm
s5u.names[(0, adjD5)] = "45"
# Gewichte einer 27: die mit A2-Gewicht = höchstes Gewicht der 3 (eine Komponente)
top3 = max(trip, key=lambda w: labels(w, A2S))
w27 = [w for w in ADJ if a2class(w) == cl3 and labels(w, A2S) == labels(top3, A2S)]
t27, _ = s5u.table(w27)
w78 = E6R + [ZERO] * 6
t78, _ = s5u.table(w78)
say(f"       E6 > SO(10) x U(1): 27 = {fmt_table(t27)};  78 = {fmt_table(t78)}")
ok27 = sorted(n for n, m, d in t27) == sorted(["(16)_1", "(10)_-2", "(1)_4"])
check("E8 > E6 x SU(3): 248 = (78,1)+(1,8)+(27,3)+(27b,3b); 27 = 16_1 + 10_-2 + 1_4", ok63 and ok27)
# Abgleich der beiden Ketten
fam_SO16 = set(w for w in SPIN128 if odd_minus5(w) and sum(w[5:]) == Fr(1, 2))
fam_E6 = set(w for w in ADJ if a2class(w) == cl3 and odd_minus5(w))
fourth = set(w for w in SPIN128 if odd_minus5(w) and sum(w[5:]) == Fr(-3, 2))
fourth_in_78 = all(w[5] == w[6] == w[7] for w in fourth)
say(f"       (16,3) aus der SO(16)-Kette und aus der E6-Kette: {len(fam_SO16)} bzw. {len(fam_E6)} Gewichte, "
    f"gleich: {fam_SO16 == fam_E6}")
say(f"       (16,1) (die vierte Familie) liegt in (78,1): {fourth_in_78} (16 = Spinoranteil von E6 > SO(10) x U(1))")
check("(16,3) aus E6 x SU(3) und aus SO(16) sind dieselben 48 Gewichte; (16,1) liegt in (78,1)",
      fam_SO16 == fam_E6 and len(fam_SO16) == 48 and fourth_in_78 and len(fourth) == 16)

# 1e: SU(5) x SU(5)
A4R1 = [r for r in D8R if all(x == 0 for x in r[5:]) and sum(r[:5]) == 0]
A4R2 = [r for r in ROOTS if r[0] == r[1] == r[2] == r[3] == r[4]]
assert len(A4R1) == 20 and len(A4R2) == 20


def chain_order(simple):
    """A_n: einfache Wurzeln entlang der Kette ordnen (für Klassen und Konjugation durch Umkehr)."""
    A = cartan(simple)
    n = len(simple)
    deg = [sum(1 for j in range(n) if j != i and A[i][j] != 0) for i in range(n)]
    start = deg.index(1)
    order, prev = [start], None
    while len(order) < n:
        cur = order[-1]
        nxt = [j for j in range(n) if j != cur and A[cur][j] != 0 and j not in order]
        order.append(nxt[0])
    return [simple[i] for i in order]


A4S1, A4S2 = chain_order(simple_roots(A4R1)), chain_order(simple_roots(A4R2))
s55 = Sub(A4S1 + A4S2)
su5a, su5b = Sub(A4S1), Sub(A4S2)
hw5 = comp_hw(su5a, 0, [v for v in VEC16 if any(x > 0 for x in v[:5])])
hw5b = tuple(reversed(hw5))
hw10a = comp_hw(su5a, 0, [tuple(VEC16[2 * i][k] + VEC16[2 * j][k] for k in range(8))
                          for i in range(5) for j in range(i + 1, 5)])
hw10ab = tuple(reversed(hw10a))
for hw, nm in ((hw5, "5"), (hw5b, "5b"), (hw10a, "10"), (hw10ab, "10b"), ((1, 0, 0, 1), "24"), ((0,) * 4, "1")):
    s55.names[(0, hw)] = nm
# SU(5)': 5' := die Partnerdarstellung der 10 in einem (10, .)-Term
br55 = s55.branch(ADJ)
for (hw, q), m in br55.items():
    if tuple(hw[:4]) == hw10a:
        hw5p = tuple(hw[4:])
hw5pb = tuple(reversed(hw5p))
for (hw, q), m in br55.items():
    if tuple(hw[:4]) == hw5b:
        hw10p = tuple(hw[4:])
hw10pb = tuple(reversed(hw10p))
for hw, nm in ((hw5p, "5'"), (hw5pb, "5'b"), (hw10p, "10'"), (hw10pb, "10'b"), ((1, 0, 0, 1), "24'"),
               ((0,) * 4, "1")):
    s55.names[(1, hw)] = nm
t55, _ = s55.table(ADJ)
say(f"   1e. E8 > SU(5) x SU(5)': 248 = {fmt_table(t55)}")
ok55 = sorted(n for n, m, d in t55) == sorted(["(24,1)", "(1,24')", "(10,5')", "(10b,5'b)", "(5b,10')",
                                                "(5,10'b)"])
# wo die Teile einer Generation liegen (Familie f0 = eine Komponente der (16,3))
gen1 = [w for w in fam_SO16 if w[5:] == max(w2[5:] for w2 in fam_SO16)]
where = Counter()
for w in gen1:
    lab_ = labels(w, A4S1 + A4S2)
    # Zuordnung über Klassen der beiden SU(5)
    c1 = sum((i + 1) * a for i, a in enumerate(lab_[:4])) % 5
    c2 = sum((i + 1) * a for i, a in enumerate(lab_[4:])) % 5
    if (c1, c2) == (0, 0):
        # (24,1) und (1,24') haben dieselbe Klasse; unterscheiden über den SU(5)-Anteil
        where["(1,24')" if all(x == 0 for x in lab_[:4]) else "(24,1)"] += 1
    else:
        where[(c1, c2)] += 1
cls = {}
for (hw, q), m in br55.items():
    c1 = sum((i + 1) * a for i, a in enumerate(hw[:4])) % 5
    c2 = sum((i + 1) * a for i, a in enumerate(hw[4:])) % 5
    cls.setdefault((c1, c2), []).append(s55.name(hw))
gen_parts = {("/".join(cls[k]) if isinstance(k, tuple) else k): v for k, v in where.items()}
say(f"       Eine Generation (16 Gewichte einer Familienkomponente) verteilt sich auf: {gen_parts}")
say("       10 der Generation in (10,5'), 5b in (5b,10'), nu^c in (1,24'): drei verschiedene SU(5)'-")
say("       Darstellungen, kein gemeinsamer Familienindex. SU(5) x SU(5)' enthält SO(10) nicht.")
check("E8 > SU(5) x SU(5)': 248 = (24,1)+(1,24')+(10,5')+(10b,5'b)+(5b,10')+(5,10'b); Generation auf drei Terme verteilt",
      ok55 and gen_parts == {"(10,5')": 10, "(5b,10')": 5, "(1,24')": 1})

# 1f: E7 x SU(2)
E7R = [r for r in ROOTS if r[6] + r[7] == 0]
A1R = [vec(0, 0, 0, 0, 0, 0, 1, 1), vec(0, 0, 0, 0, 0, 0, -1, -1)]
assert len(E7R) == 126
E7S = simple_roots(E7R)
s72 = Sub(E7S + [A1R[0]])
t72, br72 = s72.table(ADJ)
say(f"   1f. E8 > E7 x SU(2): 248 = {fmt_table(t72)}")
# SO(10)-Spinoren in (56,2) und (133,1)
in56 = [w for w in SPIN128 if odd_minus5(w) and dot(w, A1R[0]) != 0]
in133 = [w for w in SPIN128 if odd_minus5(w) and dot(w, A1R[0]) == 0]
say(f"       Gewichte der (16,4) mit SU(2)-Ladung != 0 (in (56,2)): {len(in56)}, mit Ladung 0 (in (133,1)): "
    f"{len(in133)}")
ok72 = sorted((d, m) for n, m, d in t72) == [(3, 1), (112, 1), (133, 1)]
check("E8 > E7 x SU(2): 248 = (133,1)+(1,3)+(56,2); die (16,4) verteilt sich auf beide Teile",
      ok72 and len(in56) == 32 and len(in133) == 32)

# 1g: SM-Inhalt der Generation (Fock-Raum Lambda^gerade C^5)
say("\n   1g. Eine Generation = Lambda^gerade C^5 (x_i = n_i - 1/2, gerade Besetzung = ungerade Minuszahl)")
SU3c = [vec(1, -1, 0, 0, 0, 0, 0, 0), vec(0, 1, -1, 0, 0, 0, 0, 0)]
SU2L = [vec(0, 0, 0, 1, -1, 0, 0, 0)]


def Y(w):
    return -Fr(1, 3) * (w[0] + w[1] + w[2]) + Fr(1, 2) * (w[3] + w[4])


def X(w):
    return -2 * sum(w[:5])


gen = []
for n in itertools.product((0, 1), repeat=5):
    if sum(n) % 2 == 0:
        gen.append(tuple(Fr(k) - Fr(1, 2) for k in n) + (Fr(0),) * 3)
cnt_sm = {}
for w in gen:
    cnt_sm.setdefault(Y(w), Counter())[labels(w, SU3c + SU2L)] += 1
Lsm = Lie(cartan(SU3c + SU2L))
sm_rows = []
nm3 = {(1, 0): "3", (0, 1): "3b", (0, 0): "1"}
nm2 = {(1,): "2", (0,): "1"}
for y, cnt in sorted(cnt_sm.items()):
    for hw, m in Lsm.decompose(cnt).items():
        sm_rows.append((nm3[hw[:2]], nm2[hw[2:]], y, m))
xs = sorted(set(X(w) for w in gen))
say("       " + ", ".join(f"({a},{b},{y})" for a, b, y, m in sm_rows))
expected = {("3", "2", Fr(1, 6)), ("3b", "1", Fr(-2, 3)), ("3b", "1", Fr(1, 3)), ("1", "2", Fr(-1, 2)),
            ("1", "1", Fr(1)), ("1", "1", Fr(0))}
say(f"       X = 5(B-L) - 4Y = -2 sum x_i: Werte {xs}, alle = 1 mod 4: {all(x % 4 == 1 for x in xs)}")
check("Generation = Lambda^gerade C^5: genau Q, u^c, d^c, L, e^c, nu^c; X = 1 mod 4 auf allen 16",
      set((a, b, y) for a, b, y, m in sm_rows) == expected and all(x % 4 == 1 for x in xs)
      and set(gen) == set(w[:5] + (Fr(0),) * 3 for w in fam_SO16))

# ---------------------------------------------------------------------------------------------
# 2. Gitter, Verklebung D5 + A3 (R102), Zentrum
# ---------------------------------------------------------------------------------------------
say("\n2. Gitter und Verklebung D5 + A3 (R102)\n")


def c5(w):
    return int(2 * sum(w[:5])) % 4


def c3(w):
    return int(2 * sum(w[5:])) % 4


glue = Counter((c5(w), c3(w)) for w in ADJ)
say(f"   Klassen (Z4 von Spin(10), Z4 von SU(4)) aller 248 Gewichte: {dict(sorted(glue.items()))}")
diag = all((a + b) % 4 == 0 for a, b in glue)
say(f"   E8 = D5 + A3 mit Verklebung {{(0,0), (1,3), (2,2), (3,1)}} (c5 + c3 = 0 mod 4): {diag}")
par = sorted(set((c5(v), c3(v)) for v in VEC16))
famcl = sorted(set((c5(w), c3(w)) for w in fam_SO16))
say(f"   Naht-Eltern (10,1)+(1,6): Klassen {par} -> c5 + c3 = 2: außerhalb der Verklebung")
say(f"   Generation (16,3): Klasse {famcl} -> in der Verklebung (als (16,4) in der 128)")
vec_in = [in_e8_lattice(v) for v in VEC16]
say(f"   Vektorgewichte +-e_i im E8-Gitter: {any(vec_in)}; Generationsgewichte im E8-Gitter: "
    f"{all(in_e8_lattice(w) for w in fam_SO16)}")
check("Verklebung: alle 248 Gewichte in der Diagonale c5 + c3 = 0; Naht-Eltern (2,0),(0,2) außerhalb, Generation innerhalb",
      diag and par == [(0, 2), (2, 0)] and all((a + b) % 4 == 0 for a, b in famcl))
# zentrales Element
lam = tuple(Fr(1, 2) for _ in range(8))
ph_roots = set(dot(lam, r) % 1 for r in ROOTS)
ph_vec = set(dot(lam, v) % 1 for v in VEC16)
say(f"   g = exp(2 pi i H), H = (1/2)^8 (im E8-Gitter: {in_e8_lattice(lam)}): Phasen auf allen 240 Wurzeln "
    f"{sorted(ph_roots)} -> g = 1 in E8; auf 16_v: {sorted(ph_vec)} -> g = -1.")
check("exp(2 pi i (1/2)^8) ist trivial in E8 und -1 auf dem Vektor 16: 16_v ist keine Darstellung der E8-Untergruppe Spin(16)/Z2",
      ph_roots == {0} and ph_vec == {Fr(1, 2)} and in_e8_lattice(lam))
befund("Naht-16 (Vektor von SO(16)) ist in keiner E8-Darstellung enthalten", not any(vec_in))
say("   Folgerung: Die 16 Nahtfermionen sind Bausteine (freie Fermionen vor der GSO-Projektion), nie Komponenten")
say("   eines E8-Multipletts. 'Verschiedene Komponenten einer E8' ist für die Naht-16 unmöglich.")

# ---------------------------------------------------------------------------------------------
# 3. Automorphismen, Trialität, Realität, J-Satz
# ---------------------------------------------------------------------------------------------
say("\n3. Automorphismen, Trialität, Realität\n")
T4 = [[Fr(1, 2)] * 4, [Fr(1, 2), Fr(1, 2), Fr(-1, 2), Fr(-1, 2)], [Fr(1, 2), Fr(-1, 2), Fr(1, 2), Fr(-1, 2)],
      [Fr(1, 2), Fr(-1, 2), Fr(-1, 2), Fr(1, 2)]]


def tri(w):
    a = tuple(sum(T4[i][j] * w[j] for j in range(4)) for i in range(4))
    b = tuple(sum(T4[i][j] * w[4 + j] for j in range(4)) for i in range(4))
    return a + b


rootset = set(ROOTS)
tri_ok = set(tri(r) for r in ROOTS) == rootset
tri_vec = [tri(v) for v in VEC16]
say(f"   Trialität diag(T,T) auf D4 x D4 (T vertauscht 8v und 8s): erhält die 240 Wurzeln: {tri_ok}")
say(f"   Bild von 16_v: z.B. {tuple(str(x) for x in tri_vec[0])}, Norm {dot(tri_vec[0], tri_vec[0])}, im E8-Gitter: "
    f"{any(in_e8_lattice(v) for v in tri_vec)}")
d5_img = [tri(r) for r in D5R]
d5_in_d8 = all(all(x.denominator == 1 for x in r) for r in d5_img)
say(f"   Bild der SO(10)-Wurzeln bleibt in SO(16): {d5_in_d8} (die SO(10) wird mitbewegt)")
check("Trialität ist E8-Automorphismus, bildet 16_v auf (8s,1)+(1,8s) ab (wieder nicht im E8-Gitter, Norm 1)",
      tri_ok and not any(in_e8_lattice(v) for v in tri_vec) and all(dot(v, v) == 1 for v in tri_vec))

# Realität: Irreps von Spin(10) mit Dimension <= 16
L5 = Lie(cartan(D5S))
small = []
for lab_ in itertools.product(range(3), repeat=5):
    if sum(lab_) > 3:
        continue
    # schnelle Dimension über Weyl-Formel
    num, den = Fr(1), Fr(1)
    rho = (1,) * 5
    for a in L5.pos:
        Fa = [sum(L5.F[i][j] * a[j] for j in range(5)) for i in range(5)]
        num *= sum((lab_[i] + 1) * Fa[i] for i in range(5))
        den *= sum(Fa[i] for i in range(5))
    d = num / den
    if d <= 16:
        ch = L5.char(lab_)
        selfc = set(ch) == set(tuple(-x for x in k) for k in ch) and \
            all(ch[k] == ch[tuple(-x for x in k)] for k in ch)
        small.append((int(d), lab_, selfc))
small.sort()
two_dims = []
for lab_ in itertools.product(range(3), repeat=5):
    if sum(lab_) == 2:
        num, den = Fr(1), Fr(1)
        for a in L5.pos:
            Fa = [sum(L5.F[i][j] * a[j] for j in range(5)) for i in range(5)]
            num *= sum((lab_[i] + 1) * Fa[i] for i in range(5))
            den *= sum(Fa[i] for i in range(5))
        two_dims.append(num / den)
say("   Irreps von Spin(10) mit Dimension <= 16: " + ", ".join(
    f"{d} ({'selbstkonj.' if s else 'komplex'})" for d, l_, s in small))
say(f"   (Labelsumme 2 gibt schon Dimension >= {min(two_dims)}; die Dimension wächst mit jedem Label.)")
# FS-Indikator selbstkonjugierter: (-1)^(2rho . lambda)
two_rho = [sum(L5.F[i][j] * a[j] for a in L5.pos for j in range(5)) for i in range(5)]  # scaled
fs = {}
for d, l_, s in small:
    if s:
        val = sum(l_[i] * two_rho[i] for i in range(5))
        assert val % L5.D == 0
        fs[d] = "reell" if (val // L5.D) % 2 == 0 else "pseudoreell"
say(f"   Typ der selbstkonjugierten: {fs}")
# reelle 16-dim Darstellungen: Summen, in denen komplexe Irreps paarweise auftreten
irr = [(d, s) for d, l_, s in small]
real16 = []
for k1 in range(17):
    for k10 in range(2):
        for k16 in range(2):
            for k16b in range(2):
                if k1 + 10 * k10 + 16 * (k16 + k16b) == 16 and k16 == k16b:
                    real16.append({"1": k1, "10": k10, "16": k16, "16b": k16b})
say(f"   Reelle 16-dimensionale Darstellungen von Spin(10): {real16}")
check("Reelle 16-dim Darstellungen von Spin(10) sind genau 16 x 1 und 10 + 6 x 1; die 16 ist komplex",
      sorted(tuple(sorted(r.items())) for r in real16) == sorted(
          [tuple(sorted({"1": 16, "10": 0, "16": 0, "16b": 0}.items())),
           tuple(sorted({"1": 6, "10": 1, "16": 0, "16b": 0}.items()))])
      and min(two_dims) > 16
      and any(d == 16 and not s for d, l_, s in small) and fs.get(10) == "reell")
smset = Counter((a, b, y) for a, b, y, m in sm_rows for _ in range(1))
conj = {"3": "3b", "3b": "3", "1": "1", "2": "2"}
smconj = Counter((conj[a], b, -y) for (a, b, y) in smset)
check("SM-Generation nicht selbstkonjugiert (chiral): keine reelle 16 kann sie tragen", smset != smconj)
say("   Satz R: Einschränkung einer reellen Darstellung ist reell, und jeder Automorphismus erhält den Typ.")
say("   Also trägt der Vektor 16 von SO(16) unter keiner Einbettung von Spin(10) oder G_SM eine Generation,")
say("   und kein Automorphismus (E8 hat keine äußeren; Trialität; Wahl einer anderen maximalen Untergruppe)")
say("   bildet die Naht-16 auf eine Generations-16 ab.")
befund("kein Automorphismus, keine Trialität, keine Untergruppenwahl verbindet Naht-16 und Generations-16",
       nauto == 1 and tri_ok and not any(in_e8_lattice(v) for v in tri_vec) and smset != smconj)

# J-Satz: reelles Multiplett W, komplexe Struktur J (J^2 = -1, reell), mit SO(10) vertauschend
say("\n   J-Satz (Z18-Lesart). Fermionen in Majorana-Basis: reeller Raum W mit reeller SO(10)-Wirkung; J reell,")
say("   J^2 = -1, [J, SO(10)] = 0. Dann ist xi- = konj(xi+) als SO(10)-Darstellung. Probe mit Zufalls-J:")


def gammas(n):
    s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.diag([1, -1]).astype(complex)
    I2 = np.eye(2)
    G = []
    for k in range(n // 2):
        for s in (s1, s2):
            m = np.array([[1]], complex)
            for j in range(n // 2):
                m = np.kron(m, s3 if j < k else (s if j == k else I2))
            G.append(m)
    return G


G10 = gammas(10)
Gam = np.eye(32, dtype=complex)
for g in G10:
    Gam = Gam @ g
Gam = Gam * (1j ** 5)   # Chiralität, Quadrat 1
evals, evecs = np.linalg.eigh(Gam)
Pp = evecs[:, evals > 0]
gens_sp = {}
for a in range(10):
    for b in range(a + 1, 10):
        S = 0.25 * (G10[a] @ G10[b] - G10[b] @ G10[a])
        gens_sp[(a, b)] = Pp.conj().T @ S @ Pp          # 16 x 16, antihermitesch
cart_sp = [gens_sp[(2 * k, 2 * k + 1)] for k in range(5)]
spin_w = np.round(np.real(np.linalg.eigvals(sum((k + 1.3) ** 1 * (-1j) * c for k, c in enumerate(cart_sp)))), 6)


def realify(M):
    return np.block([[M.real, -M.imag], [M.imag, M.real]])


def vec_gen(a, b, n=16):
    M = np.zeros((n, n)); M[a, b] = 1; M[b, a] = -1
    return M


def build_W(kind):
    """Liefert die 45 SO(10)-Generatoren auf W (reell) und die 5 Cartan-Generatoren."""
    gens, carts = [], []
    for (a, b), S in gens_sp.items():
        blocks = []
        if "vv" in kind:
            blocks += [vec_gen(a, b), vec_gen(a, b)]
        if "s" in kind:
            blocks += [realify(S)]
        n = sum(x.shape[0] for x in blocks)
        M = np.zeros((n, n)); o = 0
        for x in blocks:
            k = x.shape[0]; M[o:o + k, o:o + k] = x; o += k
        gens.append(M)
        if b == a + 1 and a % 2 == 0:
            carts.append(M)
    return gens, carts


def random_J(kind):
    """Zufällige reelle komplexe Struktur im Kommutanten (Schur-Blöcke)."""
    blocks = []
    if "vv" in kind:
        # 10-isotyp: 1_10 (x) j2 auf (Kopie1, Kopie2); 12 Singuletts: zufällige komplexe Struktur
        A = RNG.normal(size=(2, 2)); j2 = A - A.T
        j2 = j2 / np.sqrt(abs(np.linalg.det(j2)))
        A12 = RNG.normal(size=(12, 12)); A12 = A12 - A12.T
        # Polarzerlegung: j12 = A (A^T A)^(-1/2)
        w_, V_ = np.linalg.eigh(A12.T @ A12)
        j12 = A12 @ V_ @ np.diag(w_ ** -0.5) @ V_.T
        Jv = np.zeros((32, 32))
        idx1 = list(range(10)); idx2 = list(range(16, 26))
        for i, k in zip(idx1, idx2):
            Jv[np.ix_([i, k], [i, k])] = j2
        sing = list(range(10, 16)) + list(range(26, 32))
        Jv[np.ix_(sing, sing)] = j12
        blocks.append(Jv)
    if "s" in kind:
        sgn = RNG.choice([1, -1])
        blocks.append(sgn * np.block([[np.zeros((16, 16)), -np.eye(16)], [np.eye(16), np.zeros((16, 16))]]))
    n = sum(x.shape[0] for x in blocks)
    J = np.zeros((n, n)); o = 0
    for x in blocks:
        k = x.shape[0]; J[o:o + k, o:o + k] = x; o += k
    return J


def weights_on(J, carts, sign):
    w_, V_ = np.linalg.eig(J)
    idx = np.where(np.abs(w_ - sign * 1j) < 1e-8)[0]
    B = V_[:, idx]
    Q, _ = np.linalg.qr(B)
    mats = [Q.conj().T @ C @ Q for C in carts]
    coeff = [1.0, np.sqrt(2), np.sqrt(3), np.sqrt(5), np.sqrt(7)]
    M = sum(c * m for c, m in zip(coeff, mats))
    ev, U = np.linalg.eig(M)
    ws = []
    for k in range(U.shape[1]):
        u = U[:, k] / np.linalg.norm(U[:, k])
        ws.append(tuple(round(float(np.real(-1j * (u.conj() @ m @ u))), 6) for m in mats))
    return Counter(ws)


jsatz_ok, jsatz_comm = True, True
contents = {}
for kind in ("vv", "s", "vvs"):
    gens, carts = build_W(kind)
    for trial in range(7):
        J = random_J(kind)
        ok_alg = np.allclose(J @ J, -np.eye(J.shape[0]), atol=1e-9) and \
            max(np.abs(J @ g - g @ J).max() for g in gens) < 1e-9
        jsatz_comm &= ok_alg
        wp, wm = weights_on(J, carts, +1), weights_on(J, carts, -1)
        neg = Counter({tuple(-x for x in k): v for k, v in wp.items()})
        jsatz_ok &= (wm == +neg) or (wm == neg)
        if trial == 0:
            nspin = sum(v for k, v in wp.items() if all(abs(abs(x) - 0.5) < 1e-6 for x in k))
            nvec = sum(v for k, v in wp.items() if sum(abs(x) for x in k) == 1)
            contents[kind] = (sum(wp.values()), nvec, nspin)
say(f"   W = 16_v + 16_v (Paare), W = (16)_R, W = beides: Inhalt von xi+ (Dim, Vektorgewichte, Spinorgewichte): "
    f"{contents}")
say(f"   J^2 = -1 und [J, so(10)] = 0 in allen Proben: {jsatz_comm}; Gewichte von xi- = -(Gewichte von xi+): {jsatz_ok}")
check("J-Satz: für reelles W und reelles J im Kommutanten ist xi- = konj(xi+) (21 Zufallsproben)", jsatz_comm and jsatz_ok)
say("   Z18 nimmt Naht = xi+ (x) sqrt(q) und Familien = xi- (x) f_k, k = 0,1,2. Als SO(10)-Darstellung also")
say("   Familien = 3 x konj(Naht-Multiplett). Mit Naht = 16_v = (10,1)+(1,6) folgt Familien = 3 x [(10,1)+(1,6)]:")
say("   reell, vektorartig, keine Generation. Umgekehrt verlangt Familien = 3 x 16 die Naht in 16b.")
befund("ein reelles Multiplett mit Paarstruktur J liefert Familien = 3 x konj(Naht), nie Naht reell und Familien chiral",
       jsatz_ok and contents["vv"][2] == 0)


# GSO für die Gegenrichtung (Naht in 16b: 16 komplexe = 32 MW)
def gso_counts(N, TT=6):
    def pmul(a, b):
        out = [0] * (TT + 1)
        for i_, x in enumerate(a):
            if x:
                for j_, y in enumerate(b[:TT + 1 - i_]):
                    out[i_ + j_] += x * y
        return out
    nsp, nsm, rr = [1] + [0] * TT, [1] + [0] * TT, [1] + [0] * TT
    for n in range(1, TT + 1):
        if 2 * n - 1 <= TT:
            for _ in range(N):
                f = [0] * (TT + 1); f[0] = 1; f[2 * n - 1] = 1; nsp = pmul(nsp, f)
                f = [0] * (TT + 1); f[0] = 1; f[2 * n - 1] = -1; nsm = pmul(nsm, f)
        if 2 * n <= TT:
            for _ in range(N):
                f = [0] * (TT + 1); f[0] = 1; f[2 * n] = 1; rr = pmul(rr, f)
    shift = N // 8          # h_R = N/16 in Einheiten t = q^(1/2): Exponent N/8
    R = [0] * (TT + 1)
    for i_ in range(TT + 1 - shift):
        R[i_ + shift] = 2 ** (N // 2) * rr[i_] // 2
    return [(a + b) // 2 + c for a, b, c in zip(nsp, nsm, R)], Fr(N, 16)


g16c, h16 = gso_counts(16)
g32c, h32 = gso_counts(32)
say(f"   GSO: 16 MW -> Gewicht-1-Zustände {g16c[2]}, h_R = {h16};  32 MW (Naht in 16b) -> {g32c[2]}, h_R = {h32}")
check("Familien = 3 x 16 erzwingt Naht in 16b: 32 MW, c = 16, 496 Ströme = so(32), h_R = 2, kein E8",
      g16c[2] == 248 and g32c[2] == 496 and h32 == 2)

# ---------------------------------------------------------------------------------------------
# 4. Ströme zu Fermionen
# ---------------------------------------------------------------------------------------------
say("\n4. Ströme zu Fermionen: das Wörterbuch 128 > (16,4) > (16,3) + (16,1)\n")
# 4a Clifford
fock = []
for n in itertools.product((0, 1), repeat=8):
    fock.append((tuple(Fr(k) - Fr(1, 2) for k in n), sum(n) % 2))
even_fock = set(w for w, p in fock if p == 0)
odd_fock = set(w for w, p in fock if p == 1)
spinset = set(SPIN128)
cliff_ok = even_fock == spinset or odd_fock == spinset
four_fock = sorted(set(w[5:] for w in spinset if odd_minus5(w)))
occ = Counter(sum(1 for x in f if x > 0) for f in four_fock)
say(f"   16 Nullmoden psi_a (a = 1..16) = 8 komplexe Moden; Fockraum 256 = 128 + 128'. Die GSO-behaltene Hälfte")
say(f"   ist genau die Spinorhälfte von E8: {cliff_ok}. Die 10 Moden aus (10,1) geben Lambda^* C^5 = 16 + 16b;")
say(f"   die 6 aus (1,6) geben Lambda^* C^3 = 4 + 4b. Die 4 neben der Generation hat Besetzungszahlen {dict(occ)}")
say(f"   von C^3: 4 = Lambda^0 + Lambda^2 = 1 + 3, also (16,1) = leere, (16,3) = doppelt besetzte (1,6)-Moden.")
check("Clifford: die 128 ist der gerade Fockraum der 16 Naht-Nullmoden; Generation = Lambda^gerade C^5 der (10,1)-Moden",
      cliff_ok and occ == Counter({0: 1, 2: 3}))
say("   Das ist die genaue Beziehung der beiden 16: dim 16 der Generation = 2^(10/2 - 1) aus den zehn (10,1)-")
n_unique = [n for n in range(2, 17, 2) if 2 ** (n // 2 - 1) == 16]
say(f"   Moden; unter den Untergruppen SO(n) der 16 Moden hat nur n = {n_unique} einen Spinor der Dimension 16.")

# 4b Statistik
say("\n   4b. Statistik der Nebenklassen von D8 (bosonisiert: h = |lambda|^2/2 + ganze Oszillatoren), c = 8")
grid_i = np.array(list(itertools.product(range(-2, 3), repeat=8)), dtype=np.int64)
grid_h = np.array(list(itertools.product((-3, -1, 1, 3), repeat=8)), dtype=np.int64)  # 2x
res_cos = {}
for nm, G, half in (("ganz", grid_i, False), ("halb", grid_h, True)):
    s = G.sum(axis=1)
    nrm2 = (G ** 2).sum(axis=1)       # bei halb: 4|lambda|^2
    if half:
        par_ = (s // 2) % 2                # Summe gerade/ungerade (in Einheiten 1)
        for p_ in (0, 1):
            sel = par_ == p_
            hs = set((Fr(int(x), 8)) % 1 for x in nrm2[sel])
            minus5 = (G[sel][:, :5] < 0).sum(axis=1) % 2
            res_cos["s" if p_ == 0 else "c"] = (hs, set(minus5.tolist()))
    else:
        for p_ in (0, 1):
            sel = (s % 2) == p_
            hs = set(Fr(int(x), 2) % 1 for x in nrm2[sel])
            res_cos["0" if p_ == 0 else "v"] = (hs, None)
for k in ("0", "v", "s", "c"):
    hs = res_cos[k][0]
    spin = sorted(set((h - Fr(1, 3)) % 1 for h in hs))
    say(f"     Klasse {k}: h mod 1 = {sorted(hs)}; Drehphase eines Rings (h - c/24) mod 1 = {spin}")
stat_ok = (res_cos["0"][0] == {0} and res_cos["v"][0] == {Fr(1, 2)} and res_cos["s"][0] == {0}
           and res_cos["c"][0] == {0})
say("   SO(10)-Spinorquantenzahlen (halbzahlige Koordinaten 1..5) treten genau in s und c auf, und dort ist h ganz.")
say("   Relativ zum Grundzustand des Rings haben die 128 Ströme dieselbe Drehphase (h = 1 ganz): bosonisch.")
say("   Fermionisch (h halbzahlig) ist nur die Klasse v, und sie trägt (10,1)+(1,6), keine Generation.")
say("   In der fermionischen Fassung ohne GSO-Summe gibt es die 128 auf einem gewöhnlichen Ring (NS, er berandet")
say("   eine Scheibe) gar nicht; in der bosonischen E8-Fassung existiert sie mit h = 1. In keinem Fall ist sie")
say("   fermionisch.")
check("Nebenklassen: Klasse v hat h in Z+1/2, Klassen 0, s, c haben h ganz; SO(10)-Spinoren nur in s, c", stat_ok)
befund("Spinorströme (128 > (16,3)) sind relativ zum Ringgrundzustand bosonisch (h ganz)", res_cos["s"][0] == {0})

# 4c J auf E8 und die 3+1-Teilung
say("\n   4c. Elemente von E8, die mit SO(10) vertauschen (H = (c,c,c,c,c+d, h6,h7,h8), h in Z/16)")
grid = []
for c in (Fr(0), Fr(1, 2)):
    for d_ in (0, 1):
        for h in itertools.product(range(32), repeat=3):
            grid.append((c, c, c, c, c + d_, Fr(h[0], 16), Fr(h[1], 16), Fr(h[2], 16)))
Hn = np.array([[int(x * 32) for x in H] for H in grid], dtype=np.int64)          # 32*H
R2 = np.array([[int(2 * x) for x in r] for r in ROOTS], dtype=np.int64)           # 2*w
S2 = np.array([[int(2 * x) for x in r] for r in SPIN128], dtype=np.int64)
P_roots = (Hn @ R2.T // 2) % 32       # Exponent in Einheiten 2 pi i / 32 : 32 H.w
P_128 = (Hn @ S2.T // 2) % 32
assert ((Hn @ R2.T) % 2 == 0).all() and ((Hn @ S2.T) % 2 == 0).all()
D5n = np.array([[int(2 * x) for x in r] for r in D5R], dtype=np.int64)
comm = (((Hn @ D5n.T) // 2) % 32 == 0).all(axis=1)
say(f"     Gitterpunkte: {len(grid)}, alle mit SO(10) vertauschend (nachgeprüft: {bool(comm.all())}); das Gitter")
say("     enthält alle Elemente, deren Phasen auf der 128 Vielfache von 2 pi/32 sind.")


def order_of(row):
    g = 32
    for x in row:
        g = gcd(g, int(x))
    return 32 // g


ord248 = np.array([order_of(r) for r in P_roots])
# 3+1-Teilung mit Vorzeichen auf der 4 (bei festem SO(10)-Gewicht)
w0 = [i for i, w in enumerate(SPIN128) if w[:5] == (Fr(-1, 2), Fr(1, 2), Fr(1, 2), Fr(1, 2), Fr(1, 2))]
assert len(w0) == 4
split31 = []
for idx in range(len(grid)):
    ph = sorted(P_128[idx, w0])
    cnt = Counter(ph)
    if len(cnt) == 2:
        (a, na), (b, nb) = sorted(cnt.items(), key=lambda t: -t[1])
        if na == 3 and (a - b) % 32 == 16:
            split31.append(idx)
min_ord31 = min(ord248[split31]) if split31 else None
inv31 = [i for i in split31 if ord248[i] <= 2]
say(f"     Elemente mit 3 + 1 und relativem Vorzeichen auf der 4: {len(split31)}; davon Involutionen: {len(inv31)}; "
    f"kleinste Ordnung auf der 248: {min_ord31}")
# alle Teilungen, die Involutionen auf der 4 leisten
inv_idx = [i for i in range(len(grid)) if ord248[i] <= 2]
inv_splits = sorted(set(tuple(sorted(Counter(P_128[i, w0]).values())) for i in inv_idx))
say(f"     Teilungen der 4 durch Involutionen (Ordnung <= 2 auf der 248): {inv_splits}")
check("Kein SO(10)-vertauschendes Element der Ordnung <= 2 teilt die 4 als 3 + 1 mit Vorzeichen; kleinste Ordnung 8",
      len(inv31) == 0 and min_ord31 == 8 and (3, 1) not in [tuple(sorted(s, reverse=True)) for s in inv_splits])
say("     E8 hat triviales Zentrum: J mit J^2 = -1 kann auf keiner E8-Darstellung ein E8-Element sein; auf den")
say("     248 ist J eine Involution (Z12). Involutionen geben auf mu und -mu dieselbe Phase, trennen also 16 und")
say("     16b nicht, und sie teilen die 4 nur als 4 + 0 oder 2 + 2. Die Z18-Sortierung (Familien gegen")
say("     Nahtmode mit Vorzeichen) ist auf der (16,4) kein Element der Ordnung 2, sondern mindestens 8.")

# 4d E8-Klammer auf 128 x 128
say("\n   4d. Symmetrie der E8-Klammer 128 x 128 -> 120 (Sym^2 und Lambda^2 der 128 unter SO(16))")
lab128 = [labels(w, D8S) for w in SPIN128]
sym2, alt2 = Counter(), Counter()
for i in range(128):
    for j in range(i, 128):
        s_ = tuple(a + b for a, b in zip(lab128[i], lab128[j]))
        sym2[s_] += 1
        if j > i:
            alt2[s_] += 1
dec_sym = {d8.lie.dim(k): v for k, v in d8.lie.decompose(sym2).items()}
dec_alt = {d8.lie.dim(k): v for k, v in d8.lie.decompose(alt2).items()}
say(f"     Sym^2(128) = {dict(sorted(dec_sym.items()))};  Lambda^2(128) = {dict(sorted(dec_alt.items()))}")
check("120 liegt in Lambda^2(128), nicht in Sym^2(128): die Klammer [128,128] -> 120 ist antisymmetrisch",
      dec_alt.get(120) == 1 and 120 not in dec_sym and sum(k * v for k, v in dec_sym.items()) == 8256)
say("     Eine fermionische 128 bräuchte eine symmetrische Klammer (Superalgebra so(16) + 128). Die gibt es nicht:")
say("     Die 128 als Fermionen ist mit der E8-Struktur unverträglich.")
befund("E8-Klammer verbietet eine fermionische 128 (keine Superalgebra so(16) + 128)", 120 not in dec_sym)

# 4e Chiralität mit Fluss
say("\n   4e. Chiralität mit Fluss auf einer geschlossenen inneren Fläche (SO(10) ungebrochen)")
four_vecs = [w[5:] for w in SPIN128 if w[:5] == (Fr(-1, 2), Fr(1, 2), Fr(1, 2), Fr(1, 2), Fr(1, 2))]
trace_zero = all(sum(v[k] for v in four_vecs) == 0 for k in range(3))
s128 = Counter(SPIN128)
selfconj128 = s128 == Counter(tuple(-x for x in w) for w in SPIN128)
n16 = sum(1 for w in SPIN128 if odd_minus5(w)) // 16
n16b = sum(1 for w in SPIN128 if even_minus5(w)) // 16
say(f"     128 selbstkonjugiert: {selfconj128}; als 4D-Weylfelder #16 = {n16}, #16b = {n16b}: vektorartig")
qX = sorted(int(2 * sum(v)) for v in four_vecs)
say(f"     Gewichte der 4 summieren zu null (jede Cartan-Richtung von SU(4)): {trace_zero}; U(1)_X'-Ladungen {qX}")
say("     Index auf einer Fläche = Fluss x Ladung. Netto #16 - #16b = n x (Spur der Ladung über die 4) = 0.")
say("     Mit U(1)_X'-Fluss n: (16,3)_1 gibt 3n Familien, (16,1)_-3 gibt -3n, also 3n Anti-Familien.")
net_flux = sum(qX)
# C''-Kopplung: eine U(1), die die ganze (16,4) gleich lädt, liegt nicht in E8
M = np.array([[float(x) for x in w] for w in SPIN128 if odd_minus5(w)])
rk = np.linalg.matrix_rank(np.hstack([M, -np.ones((64, 1))]))
rk0 = np.linalg.matrix_rank(M)
say(f"     Gibt es H in R^8 mit H.w = 1 auf allen 64 Gewichten der (16,4)? Rang [W | 1] = {rk}, Rang W = {rk0}: "
    f"{'ja' if rk == rk0 else 'nein'}")
check("128 vektorartig (4 x 16 + 4 x 16b); Fluss in SU(4): netto 0, die (16,1) trägt genau die Gegenchiralität (3n - 3n); gleiche Ladung auf ganzer (16,4) liegt nicht in E8",
      trace_zero and net_flux == 0 and qX.count(1) == 3 and qX.count(-3) == 1 and rk != rk0
      and selfconj128 and n16 == n16b == 4)
befund("glatter Fluss in E8 mit ungebrochener SO(10): die vierte Familie (16,1) hebt die drei (16,3) auf",
       net_flux == 0 and qX.count(-3) == 1)
say("     Der C''-Index 3 (Grad 3 = deg K^(1/2) + deg O(D)) koppelt alle vier Komponenten gleich und braucht eine")
say("     U(1) außerhalb von E8 (oben: keine solche Richtung in R^8).")

# 4f Orbifold T^2/Z4, Konstantmoden
say("\n   4f. T^2/Z4-Orbifold ohne Fluss, komplexes 6D-Weylfeld in der 128 (Spinorphase e^(+-i pi/4)),")
say("       Eichzwirn P (vertauscht mit SO(10)), Fermionzahlphase eta = e^(2 pi i m/32).")
nets = {"E8": Counter(), "Spin16": Counter()}
example3 = None
F_mask = np.array([odd_minus5(w) for w in SPIN128])
Fb_mask = np.array([even_minus5(w) for w in SPIN128])
d8_mask = np.array([all(x.denominator == 1 for x in r) for r in ROOTS])
for idx in np.where(comm)[0]:
    pr = P_roots[idx]
    p128 = P_128[idx]
    r8 = set((p128 % 8).tolist())
    if len(r8) != 1:
        continue                       # (eta P)^4 = -1 muss auf der ganzen 128 gelten
    r = r8.pop()
    bos_spin16 = ((4 * pr[d8_mask]) % 32 == 0).all()
    bos_e8 = ((4 * pr) % 32 == 0).all()
    if not bos_spin16:
        continue
    for t in range(4):
        m = (4 - r + 8 * t) % 32
        Lsv = (4 + m + p128) % 32 == 0
        Rsv = (-4 + m + p128) % 32 == 0
        net = (int((Lsv & F_mask).sum()) - int((Lsv & Fb_mask).sum()) - int((Rsv & F_mask).sum())
               + int((Rsv & Fb_mask).sum()))
        assert net % 16 == 0
        net //= 16
        nets["Spin16"][net] += 1
        if bos_e8:
            nets["E8"][net] += 1
        if net == 4 and bos_e8 and example3 is None:
            example3 = (grid[idx], m, order_of(np.concatenate([pr, p128])), bos_e8)
say(f"     Volumen mit E8-Eichbosonen (P^4 = 1 auf 248): erreichbare Netto-#16 = {sorted(nets['E8'])}")
say(f"     Volumen nur mit Spin(16)-Eichbosonen (P^4 = 1 auf 120): erreichbare Netto-#16 = {sorted(nets['Spin16'])}")
if example3:
    H3, m3, o3, _ = example3
    say(f"     Beispiel für 4 mit E8: H = {tuple(str(x) for x in H3)}, eta-Exponent m = {m3}/32, Ordnung von P = {o3}")
all_even = all(k % 2 == 0 for k in list(nets["E8"]) + list(nets["Spin16"]))
check("Orbifold T^2/Z4 (Konstantmoden, alle mit SO(10) vertauschenden Zwirne): netto #16 stets gerade, nie 3, auch ohne 128-Bosonen",
      all_even and 3 not in nets["Spin16"] and 4 in nets["E8"] and len(nets["E8"]) > 1)
befund("Z4-Orbifold der 128 gibt nie drei Familien (netto stets gerade)", all_even)
say("     Grund: Das 6D-Weylfeld ist komplex, die 128 aber reell. (eta P)^4 = -1 muss auf (16,4) und (16b,4b)")
say("     zugleich gelten; das erzwingt m = 0 mod 4 und Zwirnphasen p_j = 0 mod 4 auf der 4. Mit det = 1 auf der")
say("     4 (sum p_j = 0 mod 32) ist die Netto-Zahl dann gerade.")
say("     Mit E8-Bosonen im Volumen ist 4 erreichbar (P = Zentrum von SU(4) auf der 4): die ganze (16,4),")
say("     also vier Familien mit gleichem Charakter, nicht drei mit {-1, i, -i}. Drei bräuchten eine halbe")
say("     Zählung (Majorana-Weyl in 6D gibt es für eine reelle Darstellung nicht) oder Fixpunktfelder.")

# ---------------------------------------------------------------------------------------------
# 5. Getrennte Zuordnung: 4D-Felder und Anomalien
# ---------------------------------------------------------------------------------------------
say("\n5. Getrennte Zuordnung (Option 3): 4D-Inhalt und Anomalien\n")
say("   Linkshändige Weylfelder in 4D:")
say("     F_k in 16 von SO(10), k = 0,1,2 (leicht, ungepaart; Charaktere {-1, i, -i} der Vierteldrehung, C'')")
say("     P_a in (10,1) + (1,6) von SO(10) x SO(6), a = 1..16 (schwer; Masse Phi P P, Phi komplex, Singulett)")
say("   Bosonen: SO(10)-Eichfelder (45), Phi (J-Ladung 2), Higgs 10_H und 126_H für die SM-Brechung (nicht Thema).")
hgen = (Fr(1), Fr(2), Fr(3), Fr(5), Fr(7))


def cubic(ws, h):
    return sum(dot(w[:5], h) ** 3 for w in ws)


h5 = hgen
a16 = cubic(gen, h5)
a10 = cubic([v for v in VEC16 if any(v[:5])], h5)
a5 = cubic([v for v in VEC16 if any(x > 0 for x in v[:5])], h5)
say(f"   SO(10)^3 (Summe (h.w)^3, generisches h): 16 -> {a16}, 10 -> {a10}; Gegenprobe 5 von SU(5) -> {a5}")
check("SO(10)^3-Anomalie verschwindet für 16 und 10 (Gegenprobe: die 5 hat sie)", a16 == 0 and a10 == 0 and a5 != 0)
# SM-Anomalien je Generation
hc = (Fr(2), Fr(-1), Fr(-1), Fr(0), Fr(0))
hw2 = (Fr(0), Fr(0), Fr(0), Fr(1), Fr(-1))
anom = {
    "Y^3": sum(Y(w) ** 3 for w in gen),
    "Y-grav": sum(Y(w) for w in gen),
    "SU3^2 Y": sum(Y(w) * dot(w[:5], hc) ** 2 for w in gen),
    "SU2^2 Y": sum(Y(w) * dot(w[:5], hw2) ** 2 for w in gen),
    "SU3^3": sum(dot(w[:5], (Fr(1), Fr(2), Fr(-3), Fr(0), Fr(0))) ** 3 for w in gen),
}
ndoub = sum(1 for w in gen if dot(w[:5], hw2) != 0) // 2
say(f"   SM je Generation: {anom}; SU(2)-Dubletts: {ndoub} (gerade: kein Witten-Defekt)")
trX = sum(X(w) for w in gen); trX3 = sum(X(w) ** 3 for w in gen)
nuX = 3 * sum(1 if X(w) % 4 == 1 else (-1 if X(w) % 4 == 3 else 0) for w in gen)
say(f"   X = 5(B-L) - 4Y: Spur {trX}, Spur X^3 {trX3}; nu_X = #(X=1) - #(X=3) mod 4 über drei Generationen = {nuX}"
    f" = {nuX % 16} mod 16")
check("Leichter Inhalt anomaliefrei: SM-Spuren 0, Dubletts gerade, tr X = tr X^3 = 0, nu_X = 48 = 0 mod 16",
      all(v == 0 for v in anom.values()) and ndoub % 2 == 0 and trX == 0 and trX3 == 0 and nuX % 16 == 0)
# Eltern
h3 = (Fr(1), Fr(3), Fr(7))
a6 = sum(dot(v[5:], h3) ** 3 for v in VEC16 if any(v[5:]))
a4 = sum(dot(f, h3) ** 3 for f in four)
say(f"   Eltern: reell. SO(6)^3 für die 6: {a6} (Gegenprobe 4: {a4}); SO(10)^3 für die 10: {a10}.")
check("Eltern (10,1)+(1,6): alle kubischen Anomalien 0 (reell), auch bei geeichter SO(6) (Gegenprobe 4 != 0)",
      a6 == 0 and a10 == 0 and a4 != 0)
# SU(3)_fam geeicht?
hA2 = (Fr(1), Fr(2), Fr(-3))
a3fam = 16 * sum(dot(f, hA2) ** 3 for f in four if sum(f) == Fr(1, 2))
say(f"   Wäre SU(3)_fam aus SU(4) geeicht: kubische Anomalie der (16,3) = 16 x A(3) ~ {a3fam} != 0.")
befund("SU(3)_fam ist mit 3 x 16 allein nicht eichbar (16 x A(3) != 0)", a3fam != 0)
# Gemeinsames Spin-Z4: J = z_X^{-1} (x) R_Paar
say("\n   Gemeinsames Spin-Z4: J = z_X^(-1) (x) R_Paar, z_X = exp(i pi X/2) (Zentrum von Spin(10)),")
say("   R_Paar dreht (psi, chi) der Eltern (Orientierung je Block so, dass xi+ auf allen 16 J = +i hat).")


def zX(Xval):
    return 1j ** (int(Xval) % 4)


fermions = []
for w in gen:
    for k in range(3):
        fermions.append(("F", zX(X(w)) ** -1))
for v in VEC16:
    xv = -2 * sum(v[:5])
    zz = zX(xv) ** -1
    eps = 1 if zz == 1 else -1           # Orientierung von R_Paar je Block
    fermions.append(("xi+", zz * (1j * eps)))
    fermions.append(("xi-", zz * (-1j * eps)))
bosons = [("45", 1), ("10_H", zX(2) ** -1), ("126_H", zX(-10) ** -1), ("Phi", -1 + 0j)]
J2f = all(abs(j ** 2 + 1) < 1e-12 for _, j in fermions)
J2b = all(abs(j ** 2 - 1) < 1e-12 for _, j in bosons)
xiplus = [j for nm, j in fermions if nm == "xi+"]
all_plus = all(abs(j - 1j) < 1e-12 for j in xiplus)
fam_ch = set(np.round(j, 6) for nm, j in fermions if nm == "F")


def nu(js):
    return sum(1 if abs(j - 1j) < 1e-9 else (-1 if abs(j + 1j) < 1e-9 else 0) for j in js)


nu_pair = nu([j for _, j in fermions])
nu_free = nu([j for nm, j in fermions if nm != "xi-"])
z16h = zX(1) ** -1
say(f"   J^2 = -1 auf allen Fermionen: {J2f}; J^2 = +1 auf 45, 10_H, 126_H, Phi: {J2b}; xi+ hat J = +i auf allen 16: "
    f"{all_plus}; Familien J = {sorted(str(complex(x)) for x in fam_ch)}")
say(f"   nu_J = n(+i) - n(-i): Paarmodell {nu_pair} = {nu_pair % 16} mod 16; partnerfreie Fassung {nu_free} = "
    f"{nu_free % 16} mod 16")
say(f"   Eine skalare 16_H hätte J^2 = {np.round(z16h ** 2, 6)}: mit Spin-Z4 unverträglich (B-L-Brechung über 126_H).")
check("Gemeinsames Spin-Z4 J = z_X^-1 (x) R_Paar: J^2 = (-1)^F auf allen Feldern, nu_J = 0 mod 16 (beide Fassungen)",
      J2f and J2b and all_plus and nu_pair % 16 == 0 and nu_free % 16 == 0 and fam_ch == {-1j})
befund("skalare 16_H ist mit dem gemeinsamen Spin-Z4 unverträglich", abs(z16h ** 2 + 1) < 1e-12)
spinZ4X_broken = any(abs(zX(-2 * sum(v[:5])) ** 2 + 1) > 1e-9 for v in VEC16)
befund("Spin x_Z2 Spin(10) ((-1)^F = z_X^2) ist durch die Eltern verletzt (X gerade auf Fermionen)", spinZ4X_broken)
# Z4_fam x SO(10)^2 Instanton
alpha = vec(1, -1, 0, 0, 0, 0, 0, 0)
n0_16 = sum(1 for w in gen if dot(w, alpha) != 0) // 2
n0_10 = sum(1 for v in VEC16 if any(v[:5]) and dot(v, alpha) != 0) // 2
qfam = [2, 1, 3]   # Charaktere -1, i, -i = i^q
inst = sum(q * n0_16 for q in qfam) % 4
grav = sum(q * 16 for q in qfam) % 4
say(f"   Nullmoden im minimalen SO(10)-Instanton: 16 -> {n0_16}, 10 -> {n0_10}. Z4_fam (Charaktere i^q, q = 2,1,3):")
say(f"   sum q n_0 = {sum(q * n0_16 for q in qfam)} = {inst} mod 4; Z4_fam-Gravitation sum q dim = {sum(q * 16 for q in qfam)}"
    f" = {grav} mod 4")
check("Z4_fam der Familiencharaktere ist von SO(10)-Instantonen und Gravitation nicht gebrochen (Ibanez-Ross)",
      n0_16 == 4 and n0_10 == 2 and inst == 0 and grav == 0)
befund("Option 3 bleibt eine direkte Summe zweier Multipletts (Naht-Eltern und Familien)",
       jsatz_ok and contents["vv"][2] == 0 and not any(in_e8_lattice(v) for v in VEC16))

# ---------------------------------------------------------------------------------------------
# 6. Urteil
# ---------------------------------------------------------------------------------------------
say("\n6. Urteil\n")
URTEIL = (
    "Obstruktion für jede Vereinigung, konsistente Zuordnung nur getrennt. (i) Die Naht-16 (Vektor von SO(16), "
    "Klassen (2,0)+(0,2) der D5+A3-Verklebung) liegt in keinem E8-Multiplett; kein Automorphismus, keine "
    "Trialität, keine andere maximale Untergruppe bildet sie auf eine Generations-16 ab (Satz R: reell bleibt "
    "reell). (ii) Ein reelles Multiplett mit Paarstruktur J gibt xi- = konj(xi+): Familien = 3 x konj(Naht); "
    "mit Naht 16_v sind die Z18-Familien 3 x [(10,1)+(1,6)], vektorartig; Familien 3 x 16 erzwingen Naht 16b "
    "mit so(32)_1 statt E8. (iii) Das Wörterbuch Ströme -> Fermionen stimmt nur in den Quantenzahlen "
    "(128 = gerader Fockraum der 16 Nullmoden, (16,4) = (16,3)+(16,1)); es scheitert an der Statistik (h ganz, "
    "J^2 = +1, E8-Klammer antisymmetrisch) und an der Chiralität (Fluss: (16,1) hebt (16,3) auf; Z4-Orbifold "
    "netto stets gerade, also nie 3). (iv) Konsistent ist die getrennte Zuordnung: 3 x 16 leicht, 16 Eltern in "
    "(10,1)+(1,6) schwer, anomaliefrei, mit gemeinsamem Spin-Z4 J = z_X^-1 R_Paar (nu_J = 0 mod 16). Die Zahl 16 "
    "steht zweimal für Z16 (zwei Spin-Z4-Zählungen) und einmal für 2^(10/2-1); kein Objekt wird doppelt benutzt, "
    "aber keines folgt aus dem anderen.")
say("   " + URTEIL)
say("\n   Verbleibende Annahmen der getrennten Zuordnung: zwei Multipletts als direkte Summe; Phi koppelt nur an")
say("   die Eltern (ersetzt den J-Projektor aus Z18); Familien koppeln an K(D) (C''); 16 Eltern aus der Forderung")
say("   E8 auf der Naht; SO(10) als Eichgruppe; Orientierung von R_Paar je Block (10 gegen 6) gewählt.")

DATA.update({
    "verzweigung_SO10xSO6": {"16": fmt_table(t16), "120": fmt_table(t120), "128": fmt_table(t128)},
    "verzweigung_SO10xSU3xU1": fmt_table(t248c),
    "verzweigung_E6xSU3": fmt_table(t63), "E6>SO10xU1": {"27": fmt_table(t27), "78": fmt_table(t78)},
    "verzweigung_SU5xSU5": fmt_table(t55), "verzweigung_E7xSU2": fmt_table(t72),
    "generation_in_SU5xSU5": gen_parts,
    "verklebung_klassen": {str(k): v for k, v in glue.items()},
    "naht_eltern_klassen": [list(p_) for p_ in par],
    "orbifold_netto": {k: sorted(v) for k, v in nets.items()},
    "orbifold_beispiel_4_E8": [str(x) for x in example3[0]] + [example3[1], example3[2]] if example3 else None,
    "sym2_128": dec_sym, "alt2_128": dec_alt,
    "j_satz_inhalt_xi_plus": {k: list(v) for k, v in contents.items()},
    "nu_J": {"paar": nu_pair, "partnerfrei": nu_free}, "nu_X": nuX,
})

say("\nZusammenfassung der Prüfungen")
for n_, c_ in CHECKS.items():
    say(f"   {'bestanden' if c_ else 'NICHT bestanden'}: {n_}")
say("Befunde (negativ oder begrenzend, nicht erzwungen)")
for n_, c_ in BEFUNDE.items():
    say(f"   {'ja' if c_ else 'nein'}: {n_}")
(HERE / "ausgabe_z23.txt").write_text("\n".join(LINES) + "\n")
(HERE / "z23_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "urteil": URTEIL,
     "annahmen": ["zwei Multipletts als direkte Summe (Naht-Eltern 16_v, Familien 3 x 16)",
                  "Phi koppelt nur an die Eltern",
                  "Familien koppeln an K(D) (C'')",
                  "16 Eltern aus der Forderung E8 auf der Naht",
                  "SO(10) als Eichgruppe",
                  "Orientierung von R_Paar je Block gewählt (J = z_X^-1 R_Paar)"],
     "data": DATA}, ensure_ascii=False, indent=1, default=str))
assert all(CHECKS.values())
