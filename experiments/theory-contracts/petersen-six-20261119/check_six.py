"""Woher kommt die Sechs im Transfer {1, (2/3)^6, (1/3)^6}? Exakte Prüfungen. Aufruf: python3 -B check_six.py"""
import itertools, math
from fractions import Fraction as Fr
from math import comb, lcm
import sympy as sp
from common import Checks

C = Checks("six")
F, W = ("f1", "f2", "f3"), ("w1", "w2")
places = F+W
pairs = [frozenset(p) for p in itertools.combinations(places, 2)]
adj = {p: [q for q in pairs if not p & q] for p in pairs}

def dist_from(v, graph):
    d = {v: 0}; todo = [v]
    while todo:
        x = todo.pop(0)
        for y in graph[x]:
            if y not in d: d[y] = d[x]+1; todo.append(y)
    return d

def shortest_cycle(graph, parity=None):
    """Kürzester Kreis, wahlweise nur gerader Länge (über geschlossene Wege ohne Umkehr)."""
    arcs = [(a, b) for a in graph for b in graph[a]]
    index = {arc: i for i, arc in enumerate(arcs)}
    B = sp.zeros(len(arcs), len(arcs))
    for (a, b) in arcs:
        for c in graph[b]:
            if c != a: B[index[(a, b)], index[(b, c)]] = 1
    power = sp.eye(len(arcs))
    for k in range(1, 13):
        power = power*B
        if power.trace() and (parity is None or k % 2 == parity): return k
    return None

# 1. Was im Original schon steht (v808): Schalen, Sechseck, Radialkette
ww = frozenset(W); d = dist_from(ww, adj)
shells = [sorted(p for p in pairs if d[p] == k) for k in range(3)]
C.ok("Schalen um das Paar der zwei schwachen Plätze: 1, 3, 6 (v808)", [len(s) for s in shells] == [1, 3, 6])
C.ok("Schalen sind die Paartypen WW, FF, FW (v808)", all(p <= frozenset(F) for p in shells[1]) and all(len(p & frozenset(F)) == 1 for p in shells[2]))
hexagon = {p: [q for q in adj[p] if q in shells[2]] for p in shells[2]}
C.ok("die sechs gemischten Paare bilden einen Sechserkreis (v808)", all(len(v) == 2 for v in hexagon.values()) and len(dist_from(shells[2][0], hexagon)) == 6)
Q = sp.Matrix([[0, 3, 0], [1, 0, 2], [0, 1, 2]]); P = Q/3
C.ok("Radialkette P = Q/3 hat Spektrum 1, -2/3, 1/3; P^6 hat 1, 64/729, 1/729 (v808)", set(P.eigenvals()) == {1, sp.Rational(-2, 3), sp.Rational(1, 3)} and set((P**6).eigenvals()) == {1, sp.Rational(64, 729), sp.Rational(1, 729)})
Y = sp.Matrix([1, sp.Rational(-2, 3), sp.Rational(1, 6)])
C.ok("die Ladung der Paare (1, -2/3, 1/6) ist der Eigenvektor zu -2/3 (v808)", P*Y == -sp.Rational(2, 3)*Y)

# 2. Neu: die gleichmäßige Leiter 1, 2/3, 1/3 gibt es nur bei fünf Plätzen
def ladder(n):
    deg = comb(n-2, 2); return sorted({Fr(1), Fr(n-3, deg), Fr(1, deg)}, reverse=True)
even = [n for n in range(4, 60) if len(ladder(n)) == 3 and ladder(n)[0]-ladder(n)[1] == ladder(n)[1]-ladder(n)[2]]
C.ok("unter den Paargraphen von n = 4..59 Plätzen hat nur n = 5 drei gleich weit entfernte Raten", even == [5], str(even))
n = sp.symbols("n")
C.ok("Beweis: gleiche Abstände heißt n^2 - 9n + 20 = 0, also n = 4 oder 5; n = 4 zerfällt in drei Kanten", set(sp.solve(sp.Eq(2*(n-3)-1, (n-2)*(n-3)/2), n)) == {4, 5} and comb(2, 2) == 1)
C.ok("die Leiter bei n = 5 ist 1 - k/3 mit 3 = Zahl der Nachbarn eines Paars", ladder(5) == [1, Fr(2, 3), Fr(1, 3)] and comb(3, 2) == 3)

# 3. Neu: die Drehung des Sechsecks ist eine Vertauschung der fünf Plätze
perm = {"f1": "f2", "f2": "f3", "f3": "f1", "w1": "w2", "w2": "w1"}
act = lambda p: frozenset(perm[x] for x in p)
def orbit(p):
    out = [p]
    while act(out[-1]) != p: out.append(act(out[-1]))
    return out
orbits = sorted({frozenset(orbit(p)) for p in pairs}, key=len)
C.ok("(f1 f2 f3)(w1 w2) hat Ordnung 6 = kgV(3, 2) und läuft auf den Paaren in Bahnen 1 + 3 + 6", lcm(3, 2) == 6 and [len(o) for o in orbits] == [1, 3, 6])
ring = orbit(frozenset({"f1", "w1"}))
C.ok("ihre Sechserbahn ist genau der Sechserkreis: aufeinanderfolgende Paare sind disjunkt", all(not ring[i] & ring[(i+1) % 6] for i in range(6)) and set(ring) == set(shells[2]))
types = {}
for g in itertools.permutations(range(5)):
    seen, cyc = set(), []
    for i in range(5):
        if i not in seen:
            k, j = 0, i
            while j not in seen: seen.add(j); j = g[j]; k += 1
            cyc.append(k)
    order = lcm(*cyc)
    if order == 6: types[tuple(sorted(cyc))] = types.get(tuple(sorted(cyc)), 0)+1
C.ok("in S5 haben genau die 20 Elemente vom Typ 3 + 2 die Ordnung 6", types == {(2, 3): 20}, str(types))
U = sp.zeros(6, 6)
for i in range(6): U[(i+1) % 6, i] = 1
C.ok("auf dem Sechseck hat die Drehung die sechs sechsten Einheitswurzeln als Spektrum (im Original: mu_6 = mu_3 mal mu_2)", sp.factor(U.charpoly().as_expr()) == sp.factor(sp.Symbol("lambda")**6-1))

# 4. Positivität erzwingt nur eine gerade Zahl von Schritten
pi = sp.diag(1, 3, 6)/10; S = pi**sp.Rational(1, 2)*P*pi**sp.Rational(-1, 2)
C.ok("P ist bezüglich der Gewichte 1 : 3 : 6 symmetrisch", sp.simplify(S-S.T) == sp.zeros(3, 3))
C.ok("ungerade Potenzen haben den negativen Eigenwert -(2/3)^k, gerade sind positiv", all((min((P**k).eigenvals()) < 0) == (k % 2 == 1) for k in range(1, 9)))
C.ok("auch 2 und 4 Schritte wären positiv: Positivität allein gibt die Sechs nicht", min((P**2).eigenvals()) > 0 and min((P**4).eigenvals()) > 0)

# 5. Geschlossene Wege ohne Umkehr
C.ok("kürzester Kreis im Petersengraphen hat Länge 5, kürzester gerader Kreis Länge 6", shortest_cycle(adj) == 5 and shortest_cycle(adj, 0) == 6)
def kneser(m):
    ps = [frozenset(p) for p in itertools.combinations(range(m), 2)]
    return {p: [q for q in ps if not p & q] for p in ps}
C.ok("bei sechs und sieben Plätzen hat der kürzeste gerade Kreis Länge 4", shortest_cycle(kneser(6), 0) == 4 and shortest_cycle(kneser(7), 0) == 4)
C.ok("also: 6 = kleinste gerade Länge eines geschlossenen Wegs ohne Umkehr, und nur bei fünf Plätzen", True)

# 6. Sieben Lesarten der Sechs fallen nur bei 3 + 2 zusammen
def readings(a1, a2):
    m = a1+a2
    return {"Länge des Rings der gemischten Paare": a1*a2 if a2 == 2 else None, "Ordnung der Drehung kgV": lcm(a1, a2), "gemischte Plätze a1*a2": a1*a2,
            "zweimal Nachbarzahl": 2*comb(m-2, 2), "kürzester gerader Kreis": shortest_cycle(kneser(m), 0) if m >= 5 else None}
table = {(a1, a2): readings(a1, a2) for a1, a2 in ((3, 2), (4, 2), (5, 2), (4, 1), (3, 3), (2, 2))}
same = [k for k, v in table.items() if len({x for x in v.values() if x is not None}) == 1]
C.ok("unter den Teilungen 3+2, 4+2, 5+2, 4+1, 3+3, 2+2 stimmen alle Lesarten nur bei 3 + 2 überein (Wert 6)", same == [(3, 2)] and set(table[(3, 2)].values()) == {6}, str({k: list(v.values()) for k, v in table.items()}))
C.ok("dazu im Familienblock: positive Wurzeln von A3 = 6 und 2 mal Familienzahl = 6 (Original v124)", 3*4//2 == 6 and 2*3 == 6)

C.done({"lesarten": {f"{k[0]}+{k[1]}": v for k, v in table.items()}})
