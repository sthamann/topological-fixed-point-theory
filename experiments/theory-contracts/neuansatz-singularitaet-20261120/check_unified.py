"""NEUSING20261120, Nachtrag: eine Regel fuer Netz und Familie? Exakt (sympy)."""
import itertools, json
import sympy as sp

rows = []
def check(name, ok, detail=""):
    rows.append({"name": name, "ok": bool(ok), "detail": str(detail)})

# --- Netz: erweitertes E8 Diagramm, Marken 1..6,4,2,3 --------------------
# Kette 1-2-3-4-5-6-4-2 und Abzweig 3 am Knoten mit Marke 6
marks = [1, 2, 3, 4, 5, 6, 4, 2, 3]
edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (5, 8)]
A = sp.zeros(9, 9)
for a, b in edges:
    A[a, b] = A[b, a] = 1
d = sp.Matrix(marks)
check("A d = 2 d (Marken sind Dimensionen)", A * d == 2 * d)
T_net = (A + 2 * sp.eye(9)) / 4          # Definition aus v395

# Klassen von SL2(F5): Ordnung -> Drehwinkel theta mit Spur 2cos(theta)
x = sp.symbols("x")
thetas = [0, sp.pi, 2*sp.pi/3, sp.pi/2, 2*sp.pi/5, 4*sp.pi/5, sp.pi/3, sp.pi/5, 3*sp.pi/5]
target = sorted([sp.nsimplify(sp.cos(t/2)**2) for t in thetas], key=lambda v: float(v))
cp = T_net.charpoly(x).as_expr()
ok = all(sp.simplify(cp.subs(x, v)) == 0 for v in target)
check("Spektrum von T_net = cos^2(theta/2) ueber die neun Klassen", ok,
      [str(sp.simplify(v)) for v in target])
# cos^2(theta/2) ist die Norm von (1+g)/2 fuer g in SU(2):
th = sp.symbols("th", real=True)
g = sp.Matrix([[sp.exp(sp.I*th), 0], [0, sp.exp(-sp.I*th)]])
H = (sp.eye(2) + g) / 2
check("(1+g)/2 mal Adjungierte = cos^2(theta/2) * 1",
      sp.simplify(H * H.H - sp.cos(th/2)**2 * sp.eye(2)) == sp.zeros(2, 2))
# Markov-Form: P_ij = d_j T_ij / d_i ist stochastisch, stationaer d_i^2/120
P = sp.Matrix(9, 9, lambda i, j: T_net[i, j] * d[j] / d[i])
check("d-konjugierter Netzschritt ist stochastisch", all(sum(P.row(i)) == 1 for i in range(9)))
pi_ = sp.Matrix([[m*m for m in marks]]) / 120
check("stationaer: Marke^2/120 (Plancherel)", pi_ * P == pi_)

# --- Familie: gleiche Form (1+tau)/2, gemittelt ueber die drei Schritte --
def perm(p):
    M = sp.zeros(4, 4)
    for i, j in enumerate(p):
        M[j, i] = 1
    return M
taus = [perm(p) for p in [(1, 0, 3, 2), (2, 3, 0, 1), (3, 2, 1, 0)]]
F = sum(((sp.eye(4) + t) / 2 for t in taus), sp.zeros(4, 4)) / 3
check("Familienregel = (1/3)I + (1/6)J", F == sp.eye(4)/3 + sp.ones(4, 4)/6)
blk = F[1:3, 1:3]
check("Paarblock 1/2, 1/6 wie v395", blk == sp.Matrix([[sp.Rational(1, 2), sp.Rational(1, 6)]]*2).applyfunc(lambda v: v) or
      blk == sp.Matrix([[sp.Rational(1, 2), sp.Rational(1, 6)], [sp.Rational(1, 6), sp.Rational(1, 2)]]))
check("Paarblock^6 hat 64/729 und 1/729",
      sorted((blk**6).eigenvals().keys()) == [sp.Rational(1, 729), sp.Rational(64, 729)])

# --- Sechs: Kandidaten, nichts davon erzwingt die Potenz ------------------
sl2f2 = [m for m in itertools.product(range(2), repeat=4) if (m[0]*m[3] - m[1]*m[2]) % 2 == 1]
check("Symmetriegruppe der 2-Torsion hat 6 Elemente", len(sl2f2) == 6)
evs = {k: sorted((blk**k).eigenvals().keys()) for k in (2, 4, 6, 8)}
check("jede gerade Potenz ist positiv; 6 ist dadurch nicht ausgezeichnet",
      all(min(v) > 0 for v in evs.values()), {k: [str(e) for e in v] for k, v in evs.items()})

# --- c3 und phi_base als Klassenanteile? ----------------------------------
sizes = [1, 1, 20, 30, 12, 12, 20, 12, 12]
sums = set()
for r in range(1, 10):
    for c in itertools.combinations(sizes, r):
        sums.add(sum(c))
check("Anteil 1/6 (phi_base = 1/(6 pi)) ist ein Klassenanteil (Ordnung 3)", 20 in sums)
check("Anteil 1/8 (c3 = 1/(8 pi)) ist KEINE Vereinigung von Klassen", 15 not in sums)

out = {"result_id": "NEUSING20261120", "part": "unified", "passed": sum(r["ok"] for r in rows),
       "total": len(rows), "rows": rows}
if __name__ == "__main__":
    json.dump(out, open("check_unified.json", "w"), indent=1, ensure_ascii=False)
    print(out["passed"], out["total"])
    for r in rows:
        if not r["ok"]:
            print("FAIL", r["name"], r["detail"])
