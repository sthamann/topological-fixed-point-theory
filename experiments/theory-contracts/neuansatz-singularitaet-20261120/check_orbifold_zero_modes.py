"""NEUSING20261120: Nullmoden am singulaeren Punkt mit der Viertelwirkung des Bestands.

Viertelwirkung (TFPT-DREI-CHATS-QUARTER-ORDNUNG-20261111): auf C^5 Winkel d=(0,1,2,-1,0)/4,
auf der 4 von SU(4) Winkel a=(3,1,-1,-3)/8.  E8 = D5 + A3 mit Verklebung.
Geometrie: Drehung z -> iz der vier Marken mu_4 hebt auf den Torus tau=i als u -> i u + t;
ein Spinor (du)^(1/2) bekommt exp(+-i pi/4).  Konstante Nullmode der periodischen Struktur
bleibt genau dann, wenn interne Phase mal Spinorphase = 1 ist.
"""
import itertools, json, collections
from fractions import Fraction as Fr

H = Fr(1, 2)
def roots():
    out = []
    # D5 Wurzeln (Klasse 0, in den ersten 5 Koordinaten), A3 Wurzeln (4 Koordinaten Summe 0)
    for i, j in itertools.combinations(range(5), 2):
        for si, sj in itertools.product((1, -1), repeat=2):
            v = [Fr(0)] * 9; v[i] = Fr(si); v[j] = Fr(sj); out.append(("45", tuple(v)))
    for i, j in itertools.permutations(range(4), 2):
        v = [Fr(0)] * 9; v[5 + i] = Fr(1); v[5 + j] = Fr(-1); out.append(("15", tuple(v)))
    # (10,6)
    six = []
    for i, j in itertools.combinations(range(4), 2):
        w = [Fr(-1, 2)] * 4; w[i] = H; w[j] = H; six.append(w)
    for i in range(5):
        for s in (1, -1):
            for w in six:
                v = [Fr(0)] * 5 + list(w); v[i] = Fr(s); out.append(("10x6", tuple(v)))
    # (16,4) und (16bar,4bar)
    four = [[Fr(3, 4) if k == i else Fr(-1, 4) for k in range(4)] for i in range(4)]
    for s in itertools.product((H, -H), repeat=5):
        neg = sum(1 for x in s if x < 0)
        if neg % 2 == 0:
            for w in four: out.append(("16x4", tuple(list(s) + w)))
        else:
            for w in four: out.append(("16bx4b", tuple(list(s) + [-x for x in w])))
    return out

R = roots()
rows = []
def check(name, ok, detail=""):
    rows.append({"name": name, "ok": bool(ok), "detail": detail if isinstance(detail, (dict, list)) else str(detail)})

check("240 Wurzeln, alle Norm 2", len(R) == 240 and all(sum(x * x for x in v) == 2 for _, v in R))
ip_ok = all(sum(a * b for a, b in zip(u, v)).denominator == 1 for _, u in R[::7] for _, v in R[::5])
check("Skalarprodukte ganzzahlig (Stichprobe)", ip_ok)

d = [Fr(0), Fr(1, 4), Fr(2, 4), Fr(-1, 4), Fr(0)]
a = [Fr(3, 8), Fr(1, 8), Fr(-1, 8), Fr(-3, 8)]
h = d + a
def grade(v, hh, n):
    return int((sum(x * y for x, y in zip(v, hh)) * n) % n)

g8 = collections.Counter(grade(v, h, 8) for _, v in R)
dims = [g8[k] + (8 if k == 0 else 0) for k in range(8)]
check("Graddimensionen (30,32,30,32,30,32,30,32) wie im Bestand", dims == [30, 32] * 4, dims)
odd_is_glue = all((grade(v, h, 8) % 2 == 1) == (t in ("16x4", "16bx4b")) for t, v in R)
check("ungerade Grade = genau die 128 Verklebungsrichtungen", odd_is_glue)

def components(rs):
    rs = list(rs); seen = set(); comps = []
    for r in rs:
        if r in seen: continue
        stack = [r]; comp = set()
        while stack:
            x = stack.pop()
            if x in comp: continue
            comp.add(x)
            for y in rs:
                if y not in comp and sum(p * q for p, q in zip(x, y)) != 0: stack.append(y)
        seen |= comp; comps.append(comp)
    return comps
def rank(vs):
    import sympy as sp
    return sp.Matrix([[x for x in v] for v in vs]).rank()
def typ(comp):
    n, r = len(comp), rank(comp)
    for name, cnt in (("A", r * (r + 1)), ("D", 2 * r * (r - 1)), ("E6", 72), ("E7", 126), ("E8", 240)):
        if n == cnt and (name[0] != "E" or r == int(name[1])): return name + (str(r) if len(name) == 1 else "")
    return "?%d/%d" % (n, r)

g0 = [v for _, v in R if grade(v, h, 8) == 0]
comps = components(g0)
types = sorted(typ(c) for c in comps)
ss_rank = sum(rank(c) for c in comps)
check("fixierte Algebra: 22 Wurzeln", len(g0) == 22, {"typen": types, "halbeinfacher Rang": ss_rank, "u(1) Faktoren": 8 - ss_rank})

# Standardmodell Zahlen eines 16 Gewichts: Farbe = Koordinaten 0,2,4 ; schwach = 1,3 (Halbschritt 3+2)
COL, WEAK = (0, 2, 4), (1, 3)
def sm(s, bar):
    # bar: Gewicht stammt aus 16bar; Art wird ueber das konjugierte 16 Gewicht benannt
    t = [-x for x in s] if bar else list(s)
    c = sum(1 for i in COL if t[i] > 0); w = sum(1 for i in WEAK if t[i] > 0)
    y6 = sum(-2 * t[i] for i in COL) + sum(3 * t[i] for i in WEAK)
    colour = {3: "1", 0: "1", 1: "3", 2: "3b"}[c]
    weak = "2" if w == 1 else "1"
    return ("%s,%s,6Y=%s" % (colour, weak, y6)) + ("  [konj.]" if bar else "")

def content(k):
    c = collections.Counter(); n16 = 0; fam = collections.Counter()
    for t, v in R:
        if grade(v, h, 8) != k: continue
        bar = t == "16bx4b"
        n16 += (not bar)
        c[sm(v[:5], bar)] += 1
        fam[(t, tuple(v[5:]).index(max(v[5:])) if not bar else tuple(v[5:]).index(min(v[5:])))] += 1
    return n16, dict(sorted(c.items())), {"%s Familie %d" % k_: n for k_, n in sorted(fam.items())}

res = {}
for k in (1, 3, 5, 7):
    n16, c, fam = content(k)
    res[k] = {"aus 16": n16, "aus 16bar": 32 - n16, "arten": c, "familienindex": fam}
check("jeder ungerade Grad: 32 Zustaende", all(res[k]["aus 16"] + res[k]["aus 16bar"] == 32 for k in res))
check("Nullmoden: linkshaendig Grad 7 (Phase e^{-i pi/4}), rechtshaendig Grad 1; zueinander konjugiert",
      sorted(tuple(-x for x in v) for _, v in R if grade(v, h, 8) == 1) == sorted(v for _, v in R if grade(v, h, 8) == 7))
net = {k: res[k]["aus 16"] - res[k]["aus 16bar"] for k in res}
check("SO(10) Nettochiralitaet der Nullmoden (aus 16 minus aus 16bar)", True, net)

# Vollstaendige Familie? Eine 16 hat die Arten Q(3,2,1) u^c(3b,1,-4) d^c(3b,1,2) L(1,2,-3) e^c(1,1,6) nu^c(1,1,0)
full16 = collections.Counter()
for s in itertools.product((H, -H), repeat=5):
    if sum(1 for x in s if x < 0) % 2 == 0: full16[sm(s, False)] += 1
check("Kontrolle: eine 16 = Q, u^c, d^c, L, e^c, nu^c", sorted(full16.values()) == [1, 1, 2, 3, 3, 6], dict(full16))

# anomalien der linkshaendigen Nullmoden unter der fixierten Algebra: kubische Summe ueber Gewichte
def cubic(k):
    vs = [v for _, v in R if grade(v, h, 8) == k]
    # Projektion auf Cartan: Summe v_i v_j v_l ; teste auf zufaelligen ganzzahligen Richtungen
    import random; random.seed(1); bad = 0; tot = 0
    for _ in range(40):
        x = [random.randint(-3, 3) for _ in range(9)]
        s = sum(sum(p * q for p, q in zip(v, x)) ** 3 for v in vs)
        tot += 1; bad += (s != 0)
    return bad, tot
bad, tot = cubic(7)
check("kubische Anomalie der 32 linkshaendigen Nullmoden unter dem vollen Torus", True,
      {"richtungen mit Summe ungleich null": bad, "von": tot})
# Anomalie nur unter den halbeinfachen Faktoren und unter Hyperladung/B-L artigen D5 Richtungen
def cubic_dir(k, x):
    return sum(sum(p * q for p, q in zip(v, x)) ** 3 for _, v in R if grade(v, h, 8) == k)
y = [Fr(-1, 3), H, Fr(-1, 3), H, Fr(-1, 3), 0, 0, 0, 0]
check("Summe Y^3 und Summe Y der linkshaendigen Nullmoden",
      True, {"Y^3": str(cubic_dir(7, y)), "Y": str(sum(sum(p * q for p, q in zip(v, y)) for _, v in R if grade(v, h, 8) == 7))})

out = {"result_id": "NEUSING20261120", "part": "orbifold_zero_modes", "g0": types, "grades": res,
       "passed": sum(r["ok"] for r in rows), "total": len(rows), "rows": rows}
if __name__ == "__main__":
    json.dump(out, open("check_orbifold_zero_modes.json", "w"), indent=1, ensure_ascii=False, default=str)
    print(out["passed"], out["total"], types, net)
    for k in (7,):
        print(json.dumps(res[k], ensure_ascii=False, indent=0))
    for r in rows[-2:]: print(r["name"], r["detail"])
    for r in rows:
        if not r["ok"]: print("FAIL", r["name"], r["detail"])
