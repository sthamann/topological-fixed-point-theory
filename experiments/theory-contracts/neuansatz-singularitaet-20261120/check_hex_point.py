"""NEUSING20261120: dieselbe Rechnung am hexagonalen Punkt, und Vorzeichen der Darstellungen am quadratischen.

Dreierdrehung z -> omega z der Marken {1, omega, omega^2, unendlich}: Funktionen z^m/(z^3-1), m=0..4,
Pullbackwinkel d3=(0,1,2,0,1)/3 (gleiche Regel wie im Bestand fuer die Viertelwirkung, dort m/4).
Familie: Dreierzyklus auf drei Gewichten, kanonischer Lift a3=(1,-1,0,0)/3.  ANNAHME: kanonischer Lift;
anders als beim Viererzyklus ist er nicht eindeutig (eine freie Phase).
"""
import json, collections, sympy as sp
from fractions import Fraction as Fr
import check_orbifold_zero_modes as Z
R = Z.R; T = {v: t for t, v in R}
add = lambda u, v: tuple(x + y for x, y in zip(u, v))
dot = lambda u, v: sum(p * q for p, q in zip(u, v))

def analyse(h, n, grades):
    g = collections.Counter(Z.grade(v, h, n) for _, v in R)
    dims = [g[k] + (8 if k == 0 else 0) for k in range(n)]
    g0 = [v for _, v in R if Z.grade(v, h, n) == 0]
    comps = sorted(Z.components(g0), key=lambda c: -len(c))
    types = [Z.typ(c) for c in comps]
    res = {}
    for k in grades:
        ws = set(v for _, v in R if Z.grade(v, h, n) == k); left = set(ws); lst = []
        while left:
            w = left.pop(); comp = {w}; st = [w]
            while st:
                x = st.pop()
                for r in g0:
                    y = add(x, r)
                    if y in ws and y not in comp: comp.add(y); st.append(y)
            left -= comp
            dd = []; sign = []
            for c in comps:
                sub = {w}; st = [w]
                while st:
                    x = st.pop()
                    for r in c:
                        y = add(x, r)
                        if y in comp and y not in sub: sub.add(y); st.append(y)
                dd.append(len(sub))
                # Vorzeichen der kubischen Form entlang einer festen Richtung des Faktors: trennt N von N-quer
                cl = sorted(c); x = tuple(sum(3 ** i * r[j] for i, r in enumerate(cl[:5])) for j in range(9))
                s = sum(dot(v, x) ** 3 for v in sub)
                sign.append(0 if s == 0 else (1 if s > 0 else -1))
            n16 = sum(1 for v in comp if T[v] == "16x4")
            lst.append({"dim": len(comp), "dims": dd, "kubisch": sign, "aus16": n16, "aus16bar": len(comp) - n16})
        res[k] = sorted(lst, key=lambda x: -x["dim"])
    return {"graddimensionen": dims, "fixierte_algebra": types, "u1": 8 - sum(Z.rank(c) for c in comps), "grade": res}

quad = analyse(Z.h, 8, (7, 3))
d3 = [Fr(0), Fr(1, 3), Fr(2, 3), Fr(0), Fr(1, 3)]
a3 = [Fr(1, 3), Fr(-1, 3), Fr(0), Fr(0)]
hexa = analyse(d3 + a3, 6, (1, 5, 3, 2, 4))
# Kontrolle: nur Familienzyklus (ohne Traegerwirkung) -> Ordnung 3, kein -1 auf der Verklebung
fam_only = analyse([Fr(0)] * 5 + a3, 3, (1, 2))
out = {"quadratisch": quad, "hexagonal": hexa, "nur_familie": fam_only}
json.dump(out, open("check_hex_point.json", "w"), indent=1, ensure_ascii=False, default=str)
for name, o in out.items():
    print(name, o["graddimensionen"], o["fixierte_algebra"], "u1:", o["u1"])
    for k, lst in o["grade"].items():
        print("  Grad", k, [(x["dims"], x["kubisch"], x["aus16"], x["aus16bar"]) for x in lst])
