"""NEUSING20261120: Entsteht E8 aus der Symmetriegruppe der Punkte der Periode vier (SL2(F5) als Einheitsquaternionen)?
Exakt in Q(sqrt5): Zahl p + q*sqrt5 als Paar (p, q) von Bruechen."""
import itertools, json
from fractions import Fraction as Fr
import sympy as sp
rows = []
def check(name, ok, detail=""): rows.append({"name": name, "ok": bool(ok), "detail": str(detail)})
Z0 = (Fr(0), Fr(0))
def add(a, b): return (a[0] + b[0], a[1] + b[1])
def mul(a, b): return (a[0] * b[0] + 5 * a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def neg(a): return (-a[0], -a[1])
def qmul(x, y):
    a1, b1, c1, d1 = x; a2, b2, c2, d2 = y
    s = lambda *t: (sum((u[0] for u in t), Fr(0)), sum((u[1] for u in t), Fr(0)))
    return (s(mul(a1, a2), neg(mul(b1, b2)), neg(mul(c1, c2)), neg(mul(d1, d2))),
            s(mul(a1, b2), mul(b1, a2), mul(c1, d2), neg(mul(d1, c2))),
            s(mul(a1, c2), neg(mul(b1, d2)), mul(c1, a2), mul(d1, b2)),
            s(mul(a1, d2), mul(b1, c2), neg(mul(c1, b2)), mul(d1, a2)))
h = Fr(1, 2); one = (Fr(1), Fr(0))
gen1 = ((h, Fr(0)),) * 4                                       # (1+i+j+k)/2, Ordnung 6
phi = (h, h); phin = (h, -h)                                   # (1+sqrt5)/2 und (1-sqrt5)/2
gen2 = ((phi[0] / 2, phi[1] / 2), (phin[0] / 2, phin[1] / 2), (h, Fr(0)), Z0)   # (phi + phi' i + j)/2
E = (one, Z0, Z0, Z0)
G = {E}; frontier = [E]
while frontier:
    nxt = []
    for x in frontier:
        for g in (gen1, gen2):
            y = qmul(x, g)
            if y not in G: G.add(y); nxt.append(y)
    frontier = nxt
check("erzeugte Gruppe von Einheitsquaternionen hat 120 Elemente", len(G) == 120, len(G))
def order(x):
    y, n = x, 1
    while y != E: y = qmul(y, x); n += 1
    return n
import collections
stat = dict(sorted(collections.Counter(order(x) for x in G).items()))
check("Ordnungsstatistik wie SL2(F5): {1:1,2:1,3:20,4:30,5:24,6:20,10:24}", stat == {1: 1, 2: 1, 3: 20, 4: 30, 5: 24, 6: 20, 10: 24}, stat)
# Form: |x|^2 = a + b sqrt5  ->  euklidische Norm 2(a+b)?  Konvention: Q(x) = 2*(a + b) mit der Spurform tr((1/sqrt5)*phi-Gewicht)
def norm2(x):
    s = Z0
    for c in x: s = add(s, mul(c, c))
    return s
def bil(x, y):
    s = Z0
    for c, d in zip(x, y): s = add(s, mul(c, d))
    return s
def lam(a, sign):       # rationale Linearform auf Q(sqrt5): a = p + q sqrt5 -> 2(p + sign*5q/... ) ; getestet werden zwei Formen
    return 2 * (a[0] + sign * a[1])
vec = lambda x: [c for comp in x for c in comp]
units = sorted(G)
res = {}
for sign in (1, -1):
    # Z-Basis aus Einheiten greedy waehlen
    basis = []
    for x in units:
        M = sp.Matrix([vec(b) for b in basis + [x]])
        if M.rank() == len(basis) + 1: basis.append(x)
        if len(basis) == 8: break
    def gram(bs): return sp.Matrix(8, 8, lambda i, j: lam(bil(bs[i], bs[j]), sign))
    # Basis verbessern: Gitter aller Einheiten = Z-Spann; nehme Hermite Normalform der Koordinaten bzgl. erster Basis
    B0 = sp.Matrix([vec(b) for b in basis]).T
    coords = [list(B0.solve(sp.Matrix(vec(x)))) for x in units]
    den = sp.ilcm(*[sp.Rational(c).q for row in coords for c in row])
    from sympy.matrices.normalforms import hermite_normal_form
    Mint = sp.Matrix([[int(c * den) for c in row] for row in coords])
    Hn = hermite_normal_form(Mint.T).T if False else None
    # einfacher: Determinante des Gitters ueber den Index: Smith Normalform der ganzzahligen Koordinatenmatrix
    from sympy.matrices.normalforms import smith_normal_form
    snf = smith_normal_form(Mint, domain=sp.ZZ)
    diag = [abs(snf[i, i]) for i in range(8)]
    index = sp.Rational(den ** 8, sp.prod(diag))          # Index des Einheitenspanns ueber dem Spann der 8 gewaehlten
    detG = gram(basis).det() / index ** 2
    even = all(lam(norm2(x), sign) % 2 == 0 for x in units)
    integral = all(sp.Rational(lam(bil(x, y), sign)).q == 1 for x in units[:40] for y in units)
    res[sign] = {"det": str(detG), "gerade": even, "ganzzahlig": integral, "norm der Einheiten": str(lam(norm2(units[0]), sign))}
check("Z-Spann der 120 Einheiten mit Form 2(a+b), |x|^2 = a + b sqrt5: gerade, ganzzahlig, Determinante 1 (also E8)",
      any(r["det"] == "1" and r["gerade"] and r["ganzzahlig"] for r in res.values()), res)
out = {"result_id": "NEUSING20261120", "part": "icosians", "passed": sum(r["ok"] for r in rows), "total": len(rows), "rows": rows}
json.dump(out, open("check_icosians.json", "w"), indent=1, ensure_ascii=False)
print(out["passed"], out["total"]); print(res)
