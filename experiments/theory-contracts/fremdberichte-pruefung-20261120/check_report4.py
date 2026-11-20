"""FREMD20261120, vierter Bericht: Familienkanal aus W3 und D3, V4 Zerlegung von End_0(F). Exakt."""
import json, itertools, sympy as sp
rows = []
def check(name, ok, detail=""): rows.append({"name": name, "ok": bool(ok), "detail": str(detail)})
I3 = sp.eye(3); h = sp.Rational
def Phi(rho):
    W = (rho.trace() * I3 - rho.T) / 2
    D = rho.trace() * I3 / 3
    return rho / 2 + W / 3 + D / 6
# Clockbasis der Vierteldrehung u = [[0,-1,0],[1,0,0],[0,0,1]]
u = sp.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
vs = [sp.Matrix([1, -sp.I, 0]) / sp.sqrt(2), sp.Matrix([1, sp.I, 0]) / sp.sqrt(2), sp.Matrix([0, 0, 1])]
check("Clockbasis: Eigenvektoren von u", all(sp.simplify(u * v - lam * v) == sp.zeros(3, 1) for v, lam in zip(vs, (sp.I, -sp.I, 1))))
B = sp.Matrix(3, 3, lambda i, j: sp.simplify((vs[i].H * Phi(vs[j] * vs[j].H) * vs[i])[0]))
check("Populationswirkung in der Clockbasis = (1/18)[[13,1,4],[1,13,4],[4,4,10]]", 18 * B == sp.Matrix([[13, 1, 4], [1, 13, 4], [4, 4, 10]]), (18 * B).tolist())
Bstd = sp.Matrix(3, 3, lambda i, j: Phi(sp.Matrix(3, 3, lambda a, b: 1 if a == b == j else 0))[i, i])
check("in der reellen Standardbasis dagegen (1/3)I + (2/9)J (Raten 1, 1/3, 1/3, nicht B)", Bstd == I3 / 3 + sp.ones(3, 3) * h(2, 9))
S = sp.zeros(9, 9)
for j, (a, b) in enumerate(itertools.product(range(3), repeat=2)):
    E = sp.zeros(3, 3); E[a, b] = 1; S[:, j] = Phi(E).reshape(9, 1)
ev = S.eigenvals()
check("volles Spektrum 1, (2/3) x3, (1/3) x5", ev == {sp.Integer(1): 1, h(2, 3): 3, h(1, 3): 5}, ev)
# vollstaendig positiv: Choi Matrix
C = sp.zeros(9, 9)
for a, b in itertools.product(range(3), repeat=2):
    E = sp.zeros(3, 3); E[a, b] = 1
    C += sp.kronecker_product(E, Phi(E))
check("Choi Matrix positiv semidefinit, spurtreu", min(C.eigenvals()) >= 0 and all(Phi(sp.Matrix(3, 3, lambda i, j: 1 if (i, j) == (a, b) else 0)).trace() == (1 if a == b else 0) for a in range(3) for b in range(3)), sorted(C.eigenvals().items()))
# V4 = diag Vorzeichen (det 1) auf F; Zerlegung der spurfreien Operatoren
V4 = [sp.diag(1, 1, 1), sp.diag(1, -1, -1), sp.diag(-1, 1, -1), sp.diag(-1, -1, 1)]
def chars(basis):
    out = []
    for X in basis:
        out.append(tuple(1 if g * X * g == X else -1 for g in V4))
    return sorted(out)
sym = [sp.diag(1, -1, 0), sp.diag(1, 1, -2)] + [sp.Matrix(3, 3, lambda i, j: 1 if {i, j} == {a, b} else 0) for a, b in ((0, 1), (0, 2), (1, 2))]
asym = [sp.Matrix(3, 3, lambda i, j: (1 if (i, j) == (a, b) else -1 if (i, j) == (b, a) else 0)) for a, b in ((0, 1), (0, 2), (1, 2))]
triv = (1, 1, 1, 1)
check("Quadrupole: 2 mal trivial + drei verschiedene Charaktere", chars(sym).count(triv) == 2 and len(set(chars(sym))) == 4)
check("Dipole: drei verschiedene nichttriviale Charaktere", triv not in chars(asym) and len(set(chars(asym))) == 3)
comm = lambda X, Y: X * Y - Y * X
is_sym = lambda X: X.T == X; is_asym = lambda X: X.T == -X
check("[5,5] in 3, [3,5] in 5, [3,3] in 3", all(is_asym(comm(X, Y)) for X in sym for Y in sym) and all(is_sym(comm(X, Y)) and comm(X, Y).trace() == 0 for X in asym for Y in sym) and all(is_asym(comm(X, Y)) for X in asym for Y in asym))
out = {"result_id": "FREMD20261120", "part": "report4", "passed": sum(r["ok"] for r in rows), "total": len(rows), "rows": rows}
json.dump(out, open("experiments/theory-contracts/fremdberichte-pruefung-20261120/check_report4.json", "w"), indent=1, ensure_ascii=False)
print(out["passed"], out["total"]); [print("FAIL", r["name"], r["detail"]) for r in rows if not r["ok"]]
