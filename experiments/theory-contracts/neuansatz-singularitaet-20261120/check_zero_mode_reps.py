"""NEUSING20261120: Darstellungsinhalt der Nullmoden unter der fixierten Algebra (exakt)."""
import json, collections, itertools, sympy as sp
from fractions import Fraction as Fr
import check_orbifold_zero_modes as Z

R, h = Z.R, Z.h
T = {v: t for t, v in R}
g0 = [v for _, v in R if Z.grade(v, h, 8) == 0]
comps = sorted(Z.components(g0), key=lambda c: -len(c))
names = []
for c in comps:
    src = collections.Counter(T[v] for v in c)
    names.append({"typ": Z.typ(c), "wurzeln": len(c), "herkunft": dict(src)})
add = lambda u, v: tuple(x + y for x, y in zip(u, v))

def irreps(k):
    ws = set(v for _, v in R if Z.grade(v, h, 8) == k)
    out = []; left = set(ws)
    while left:
        w = left.pop(); comp = {w}; st = [w]
        while st:
            x = st.pop()
            for r in g0:
                y = add(x, r)
                if y in ws and y not in comp: comp.add(y); st.append(y)
        left -= comp
        dims = []
        for c in comps:
            sub = {w}; st = [w]
            while st:
                x = st.pop()
                for r in c:
                    y = add(x, r)
                    if y in comp and y not in sub: sub.add(y); st.append(y)
            dims.append(len(sub))
        n16 = sum(1 for v in comp if T[v] == "16x4")
        out.append({"dim": len(comp), "dims": dims, "aus16": n16, "aus16bar": len(comp) - n16, "w": comp})
    return out

# Projektion auf den halbeinfachen Teil: Vergleich einer Darstellung mit ihrer Konjugierten
B = sp.Matrix([list(v) for v in g0])
ns = B.nullspace()          # Richtung des u(1)
u1 = ns[0] if ns else None
def proj(v):
    v = sp.Matrix(v)
    if u1 is None: return tuple(v)
    return tuple(v - (v.dot(u1) / u1.dot(u1)) * u1)
def charge(v):
    return sp.Matrix(v).dot(u1) / sp.sqrt(u1.dot(u1)) if u1 is not None else 0

res = {}
for k in (7, 1, 3, 5):
    lst = irreps(k)
    res[k] = [{"dim": x["dim"], "dims(%s)" % ",".join(n["typ"] for n in names): x["dims"],
               "aus16": x["aus16"], "aus16bar": x["aus16bar"],
               "u1": str(sp.nsimplify(charge(next(iter(x["w"])))))} for x in lst]
# chiral unter dem halbeinfachen Teil?  (Gewichtsmenge gleich ihrer negativen)
def ss_set(k):
    return collections.Counter(proj(v) for _, v in R if Z.grade(v, h, 8) == k)
L = ss_set(7)
Lneg = collections.Counter(tuple(-x for x in w) for w in L.elements())
real_ss = (L == Lneg)
# kubische Anomalie je einfachem Faktor: Summe (w.x)^3 fuer x im Span des Faktors
def anomaly(c):
    basis = list(c)[:6]; tot = []
    for x in basis:
        for y in basis:
            tot.append(sum((sum(p*q for p, q in zip(v, x)) ** 2) * sum(p*q for p, q in zip(v, y))
                           for _, v in R if Z.grade(v, h, 8) == 7))
    return any(t != 0 for t in tot)
anom = {n["typ"] + "_%d" % i: anomaly(c) for i, (n, c) in enumerate(zip(names, comps))}
out = {"faktoren": names, "nullmoden": res, "reell unter halbeinfachem Teil": real_ss, "kubische Anomalie je Faktor": anom}
json.dump(out, open("check_zero_mode_reps.json", "w"), indent=1, ensure_ascii=False, default=str)
print(json.dumps(names, ensure_ascii=False))
for k in (7, 3): print(k, json.dumps(res[k], ensure_ascii=False))
print("reell unter halbeinfachem Teil:", real_ss, anom)
