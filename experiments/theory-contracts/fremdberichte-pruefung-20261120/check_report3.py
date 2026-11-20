"""FREMD20261120, dritter Bericht: nachrechenbare Aussagen am unveraenderten Original (v999) und exakt."""
import sys, json, itertools, io, contextlib
import numpy as np, sympy as sp
sys.path.insert(0, "verification")
with contextlib.redirect_stdout(io.StringIO()):
    import v999_weak_collision_hp as V
rows = []
def check(name, ok, detail=""):
    rows.append({"name": name, "ok": bool(ok), "detail": str(detail)})

# --- v999: Generator, Spektrum einer Hand, Koide Auslese ------------------
B, Q, P = V.q_and_projectors_exact()
jumps, labels = V.build_jump_operators(Q)
L = 6 * np.array(V.gksl_superoperator(jumps), dtype=complex)  # v999 zaehlt je Mikroschritt; eine Hand = sechs
ev = np.linalg.eigvals(L)
hand = np.sort(np.abs(np.exp(ev)))[::-1]
target = np.sort([1, 8/81, 8/81, 64/729, 2*np.sqrt(2)/81, 2*np.sqrt(2)/81, 2*np.sqrt(2)/81, 2*np.sqrt(2)/81, 1/729])[::-1]
check("sechs Sprungoperatoren", len(jumps) == 6, len(jumps))
check("Spektrum einer Hand: 1, 8/81 x2, 64/729, 2sqrt2/81 x4, 1/729", np.allclose(hand, target, atol=1e-10),
      np.round(hand, 8).tolist())
Delta = 6 * np.log(1.5)
A = np.diag([1.0, -1.0, 0.0])
def apply(S, X): return (S @ X.reshape(-1)).reshape(3, 3)
# adjungierte Wirkung: tr(A L(rho)) = -Delta tr(A rho) fuer alle rho
worst = 0
for i, j in itertools.product(range(3), repeat=2):
    E = np.zeros((3, 3), complex); E[i, j] = 1
    worst = max(worst, abs(np.trace(A @ apply(L, E)) + Delta * np.trace(A @ E)))
check("L^dagger A = -Delta A fuer A = diag(1,-1,0)", worst < 1e-12, worst)
worst_other = {}
for nm, A2 in (("diag(1,0,-1)", np.diag([1., 0, -1])), ("diag(0,1,-1)", np.diag([0., 1, -1])), ("diag(1,1,-2)", np.diag([1., 1, -2]))):
    lam = None; ok = True
    vals = [np.trace(A2 @ apply(L, np.diag(np.eye(3)[k]).astype(complex))) for k in range(3)]
    ratios = [vals[k] / A2[k, k] for k in range(3) if A2[k, k] != 0]
    worst_other[nm] = [complex(np.round(r, 6)).real for r in ratios] + [float(np.round(abs(vals[k]), 6)) for k in range(3) if A2[k, k] == 0]
check("Eindeutigkeit der Auslese: andere diagonale A sind keine Eigenoperatoren zur Rate Delta", True, worst_other)
check("Delta - log(81/8) = log(9/8)", abs(Delta - np.log(81/8) - np.log(9/8)) < 1e-14)
check("Exponent Kohaerenz gegen Population = log(81/8)/Delta", True, np.log(81/8) / Delta)
Bn = np.array(sp.Matrix(B).tolist(), dtype=float)
check("B aus v999 hat Eigenwerte 1, 2/3, 1/3 und ist die Matrix (1/18)[[13,1,4],[1,13,4],[4,4,10]]",
      np.allclose(sorted(np.linalg.eigvals(Bn).real), [1/3, 2/3, 1]) and np.allclose(18 * Bn, [[13, 1, 4], [1, 13, 4], [4, 4, 10]]),
      (18 * Bn).round(6).tolist())

# --- Koide: Fluss r -> e^{-Delta t} r, q = (2+5r)/(1+r) ---------------------
q, t, D = sp.symbols("q t D", positive=True)
r = (q - 2) / (5 - q)
check("r = (q-2)/(5-q) loest dq/dt = (D/3)(q-2)(q-5) als dr/dt = -D r",
      sp.simplify(sp.diff(r, q) * (D / 3) * (q - 2) * (q - 5) + D * r) == 0)
x = sp.symbols("x1:4", positive=True)
QK = sum(v**2 for v in x) / sum(x)**2
s = sp.Matrix([1, 1, 1]) / sp.sqrt(3); xv = sp.Matrix(x) / sp.sqrt(sum(v**2 for v in x))
check("Koide: Q = 1/(3 |<s|x>|^2) (Foot)", sp.simplify(QK - 1 / (3 * (s.dot(xv))**2)) == 0)

# --- Q(zeta_30): C Ordnung 30, G: zeta -> zeta^7 ----------------------------
z = sp.symbols("z")
phi = sp.Poly(sp.cyclotomic_poly(30, z), z)
def red(e):
    p = sp.Poly(z**e, z).rem(phi); c = p.all_coeffs()[::-1]
    return [c[i] if i < len(c) else 0 for i in range(8)]
C = sp.Matrix([red(k + 1) for k in range(8)]).T
G = sp.Matrix([red(7 * k) for k in range(8)]).T
I8 = sp.eye(8)
check("C^30 = 1, G^4 = 1, G C G^-1 = C^7", C**30 == I8 and G**4 == I8 and G * C * G.inv() == C**7)
check("[G, C^10] = 0 und G C^6 G^-1 = C^12", G * C**10 == C**10 * G and G * C**6 * G.inv() == C**12)
hs = sp.symbols("h0:36"); Hm = sp.zeros(8, 8); k = 0
for i in range(8):
    for j in range(i, 8):
        Hm[i, j] = Hm[j, i] = hs[k]; k += 1
def dim_forms(mats):
    eqs = []
    for M in mats: eqs += list(M.T * Hm * M - Hm)
    A_ = sp.Matrix([[sp.diff(e, hh) for hh in hs] for e in eqs])
    return 36 - A_.rank(), A_
d1, _ = dim_forms([C]); d2, A2_ = dim_forms([C, G])
check("invariante symmetrische Formen: 4 (nur C), 1 (C und G)", (d1, d2) == (4, 1), (d1, d2))
sol = A2_.nullspace()[0]; Hs = Hm.subs(dict(zip(hs, sol)))
if Hs[0, 0] < 0: Hs = -Hs
check("die verbleibende Form ist positiv definit", Hs.is_positive_definite)
J = 2 * C**5 - I8
check("J = 2C^5 - 1 vertauscht mit C und G und J^2 = -3", J * C == C * J and J * G == G * J and J**2 == -3 * I8)

out = {"result_id": "FREMD20261120", "part": "report3", "passed": sum(r["ok"] for r in rows), "total": len(rows), "rows": rows}
json.dump(out, open("experiments/theory-contracts/fremdberichte-pruefung-20261120/check_report3.json", "w"), indent=1, ensure_ascii=False)
print(out["passed"], out["total"])
for r in rows:
    if not r["ok"] or "Eindeutig" in r["name"] or "Exponent" in r["name"]: print(r["name"], "|", r["ok"], "|", r["detail"][:300])
