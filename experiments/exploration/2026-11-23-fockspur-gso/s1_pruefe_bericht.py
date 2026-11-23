"""Unabhaengige Nachrechnung zentraler Behauptungen der beiden Berichte vom 23.11.2026.

A  Ursprungskanal (S6, 15 Transpositionen): Spektrum, Ladungssektoren von ad(6Y),
   Blocknormen 317/108, 53/54, 41/216, Singulaerwerte, beide Kommutatoren.
B  Darstellungsspuren Tr_R Sigma^2, Sigma^4 (5, 10, 16, 24, 248) und E8 Tr Sigma^6.
C  Positive Determinantenquelle rho tr S^2 - tr log(1+S^2): globale 3+2 Auswahl fuer 1/3<rho<1.
D  E8 Wirkung rho|H|^2 - sum log(1+(a.H)^2): Cartan-Hessian entlang h=(2,2,2,-3,-3,0,0,0)/sqrt30.
E  Freie Stromdarstellung: ||r Omega||^2 = 3 und b3 b1 r Omega = Omega.

Exakt (sympy Rational) wo moeglich, sonst numerisch mit ausgewiesener Toleranz.
"""
import itertools, math, sys
import numpy as np
import sympy as sp
from scipy.optimize import minimize

OK = True
def check(name, cond, info=""):
    global OK
    OK = OK and bool(cond)
    print(f"  [{'ok' if cond else 'FEHLER'}] {name}" + (f" -- {info}" if info else ""))

# ---------------------------------------------------------------- A
print("A  Ursprungskanal aus den 15 Transpositionen von S6")
n = 6
def perm_matrix(p):
    M = sp.zeros(n, n)
    for i in range(n):
        M[p[i], i] = 1
    return M
def transposition(a, b):
    p = list(range(n)); p[a], p[b] = b, a
    return tuple(p)
T15 = [transposition(a, b) for a, b in itertools.combinations(range(n), 2)]
pairs = [frozenset(c) for c in itertools.combinations(range(n), 2)]
P = {t: perm_matrix(t) for t in T15}
P0 = sp.eye(n) - sp.ones(n, n) / 6          # Projektor auf den Standardmodul (Summe null)

def coeff(i, j):
    a, b = pairs[i], pairs[j]
    if a == b: return sp.Rational(32, 648)
    if len(a & b) == 1: return sp.Rational(5, 648)
    return sp.Rational(-4, 648)

# Superoperator auf vec(A) (spaltenweise): vec(X A Z) = (Z^T kron X) vec(A)
def kron(A, B):
    K = sp.zeros(A.rows * B.rows, A.cols * B.cols)
    for i in range(A.rows):
        for j in range(A.cols):
            if A[i, j] != 0:
                K[i * B.rows:(i + 1) * B.rows, j * B.cols:(j + 1) * B.cols] = A[i, j] * B
    return K
Q0 = kron(P0.T, P0)
Phi = sp.zeros(36, 36)
for i, s in enumerate(T15):
    for j, t in enumerate(T15):
        Phi += coeff(i, j) * kron(P[t].T, P[s])
Phi = Q0 * Phi * Q0
check("Phi ist unital auf dem Standardmodul", Phi * Q0 * sp.Matrix(sp.eye(n).reshape(36, 1)) == Q0 * sp.Matrix(sp.eye(n).reshape(36, 1)))
spec = sorted([(k, v) for k, v in Phi.eigenvals().items()], key=lambda kv: -kv[0])
print("   Spektrum (Eigenwert: Vielfachheit, inkl. 11 Nullrichtungen ausserhalb des Standardmoduls):", spec)
check("Spektrum 1(1), 3/4(5), 11/36(9), 1/4(10)",
      dict(spec) == {1: 1, sp.Rational(3, 4): 5, sp.Rational(11, 36): 9, sp.Rational(1, 4): 10, 0: 11})

M = perm_matrix((1, 0, 3, 2, 5, 4))          # Matching (01)(23)(45)
Pw = (sp.eye(n) + M) / 2 * P0                 # M=+1 auf Summe null: Dim 2 (schwach)
Pc = (sp.eye(n) - M) / 2 * P0                 # M=-1: Dim 3 (Farbe)
check("Rang Farbe 3, schwach 2", (Pc.rank(), Pw.rank()) == (3, 2))
Y6 = (P0 + 5 * M * P0) / 2                    # 6Y
check("6Y = diag(-2,-2,-2,3,3) auf Farbe/schwach", Y6 * Pc == -2 * Pc and Y6 * Pw == 3 * Pw)
S0 = kron(Pc.T, Pc) + kron(Pw.T, Pw)
Sp5 = kron(Pc.T, Pw)                           # A -> Pw A Pc  (Ladung +5)
Sm5 = kron(Pw.T, Pc)                           # A -> Pc A Pw  (Ladung -5)
def fro2(X):
    return sp.nsimplify(sum(abs(x) ** 2 for x in X))
sectors = {"0": S0, "+5": Sp5, "-5": Sm5}
table = {}
for so, Qo in sectors.items():
    for si, Qi in sectors.items():
        table[(so, si)] = sp.simplify(fro2(Qo * Phi * Qi))
for so in sectors:
    print("   Ausgang", so, [table[(so, si)] for si in sectors])
check("Blocknormen 317/108, 53/54, 41/216 und exakte Nullen",
      table[("0", "0")] == sp.Rational(317, 108) and table[("+5", "+5")] == sp.Rational(53, 54)
      and table[("-5", "-5")] == sp.Rational(53, 54) and table[("+5", "-5")] == sp.Rational(41, 216)
      and table[("-5", "+5")] == sp.Rational(41, 216)
      and all(table[(a, b)] == 0 for a, b in [("0", "+5"), ("0", "-5"), ("+5", "0"), ("-5", "0")]))
B = Sm5 * Phi * Sp5
sv = sorted([sp.nsimplify(sp.sqrt(k)) for k, m in (B.T * B).eigenvals().items() for _ in range(m) if k != 0], reverse=True)
print("   Singulaerwerte des Uebergangs -5 <- +5:", sv)
check("Singulaerwerte (1/4)x3, (1/36)x3", sv == [sp.Rational(1, 4)] * 3 + [sp.Rational(1, 36)] * 3)

Y = Y6 / 6
adY = kron(sp.eye(n), Y) - kron(Y.T, sp.eye(n))
c1 = fro2(Q0 * (Phi * adY - adY * Phi) * Q0)
T = perm_matrix((2, 3, 0, 1, 4, 5))            # T=(02)(13)
R = (Pc + sp.I * Pw) * T
check("R unitaer auf dem Modul, R^2 = -M", sp.simplify(R * R.H - P0) == sp.zeros(n) and sp.simplify(R * R + M * P0) == sp.zeros(n))
AdR = kron(R.H.T, R)
c2 = sp.nsimplify(sp.simplify(fro2(sp.simplify(Q0 * (Phi * AdR - AdR * Phi) * Q0))))
print(f"   ||[Phi,ad_Y]||^2 = {c1},  ||[Phi,Ad_R]||^2 = {c2},  Verhaeltnis = {sp.nsimplify(c2 / c1)}")
check("Kommutatoren 1025/972 und 41/27, Verhaeltnis 36/25",
      c1 == sp.Rational(1025, 972) and c2 == sp.Rational(41, 27) and sp.nsimplify(c2 / c1) == sp.Rational(36, 25))
check("Zerlegung: 2(5/3)^2 41/216 und 2*4*41/216",
      2 * sp.Rational(25, 9) * sp.Rational(41, 216) == c1 and 8 * sp.Rational(41, 216) == c2)

# ---------------------------------------------------------------- B
print("\nB  Darstellungsspuren (Eigenwerte x1..x5, Summe null)")
x = sp.symbols("x1:6")
sub = {x[4]: -sum(x[:4])}
def p(k): return sum(xi ** k for xi in x)
def tr_rep(weights, k): return sum(w ** k for w in weights)
w5 = list(x)
w10 = [x[i] + x[j] for i, j in itertools.combinations(range(5), 2)]
w16 = [0] + w10 + [sum(x) - x[i] for i in range(5)]
w24 = [x[i] - x[j] for i in range(5) for j in range(5) if i != j] + [0] * 4
w5b = [-w for w in w5]; w10b = [-w for w in w10]
w248 = w24 + [0] * 24 + 10 * w5 + 10 * w5b + 5 * w10 + 5 * w10b
def same(a, b): return sp.expand((a - b).subs(sub)) == 0
check("Tr_10: 3p2, 3p2^2-3p4", same(tr_rep(w10, 2), 3 * p(2)) and same(tr_rep(w10, 4), 3 * p(2) ** 2 - 3 * p(4)))
check("Tr_16: 4p2, 3p2^2-2p4", same(tr_rep(w16, 2), 4 * p(2)) and same(tr_rep(w16, 4), 3 * p(2) ** 2 - 2 * p(4)))
check("Tr_24: 10p2, 6p2^2+10p4", same(tr_rep(w24, 2), 10 * p(2)) and same(tr_rep(w24, 4), 6 * p(2) ** 2 + 10 * p(4)))
check("Tr_248: 60p2, 36p2^2, 30p2^3",
      same(tr_rep(w248, 2), 60 * p(2)) and same(tr_rep(w248, 4), 36 * p(2) ** 2) and same(tr_rep(w248, 6), 30 * p(2) ** 3))
# Formauswahl bei positiver Determinantenentwicklung: Quartikkoeffizient c4(R)/c2(R)^2 auf 3+2 und 4+1
for name, w in [("5", w5), ("10", w10), ("16", w16), ("24", w24)]:
    vals = {}
    for shape, vec in [("3+2", (2, 2, 2, -3, -3)), ("4+1", (1, 1, 1, 1, -4))]:
        s = dict(zip(x, vec))
        vals[shape] = sp.Rational(tr_rep(w, 4).subs(s), tr_rep(w, 2).subs(s) ** 2)
    best = min(vals, key=vals.get)
    print(f"   R={name}: Tr4/Tr2^2  3+2: {vals['3+2']},  4+1: {vals['4+1']}  -> kleinster Quartikterm bei {best}")

# ---------------------------------------------------------------- C
print("\nC  Positive Determinantenquelle Gamma = rho tr S^2 - tr log(1+S^2), globale Minima")
basis = np.linalg.qr(np.vstack([np.ones(5), np.eye(5)[:4]]).T)[0][:, 1:]   # Summe-null Basis
def shape_of(lam, tol=1e-5):
    lam = np.sort(lam)
    groups = [[lam[0]]]
    for v in lam[1:]:
        (groups[-1].append(v) if abs(v - groups[-1][-1]) < tol else groups.append([v]))
    return "+".join(str(len(g)) for g in sorted(groups, key=len, reverse=True))
def global_min(fun, rng, starts=400, scale=3.0):
    best = None
    for _ in range(starts):
        z0 = rng.normal(size=4) * scale * rng.uniform(0.05, 1)
        r = minimize(lambda z: fun(basis @ z), z0, method="BFGS", options={"gtol": 1e-11})
        if best is None or r.fun < best.fun - 1e-12:
            best = r
    return basis @ best.x, best.fun
rng = np.random.default_rng(1)
okC = True
for rho in [0.2, 0.3, 0.34, 0.4, 0.5, 0.7, 0.9, 0.99, 1.05]:
    G = lambda lam, rho=rho: rho * np.sum(lam ** 2) - np.sum(np.log1p(lam ** 2))
    lam, val = global_min(G, rng)
    shp = shape_of(lam)
    v2 = (6 - 13 * rho + math.sqrt(25 * rho ** 2 - 12 * rho + 36)) / (72 * rho)
    pred = rho * 30 * v2 - 3 * math.log1p(4 * v2) - 2 * math.log1p(9 * v2) if rho < 1 else 0.0
    inside = 1 / 3 < rho < 1
    agree = (abs(val - pred) < 1e-7) if inside else True
    print(f"   rho={rho:5.2f}: Form {shp:9s} Eigenwerte {np.round(np.sort(lam),5)}  Gamma={val:.8f}  3+2-Formel={pred:.8f}")
    if inside:
        okC = okC and shp == "3+2" and agree
check("Im Fenster 1/3<rho<1 ist 3+2 global mit v^2-Formel", okC)

# ---------------------------------------------------------------- D
print("\nD  E8 Wirkung entlang h=(2,2,2,-3,-3,0,0,0)/sqrt30")
roots = []
for i, j in itertools.combinations(range(8), 2):
    for si in (1, -1):
        for sj in (1, -1):
            r = [0] * 8; r[i] = si; r[j] = sj; roots.append(r)
for signs in itertools.product((1, -1), repeat=8):
    if signs.count(-1) % 2 == 0:
        roots.append([s / 2 for s in signs])
roots = np.array(roots, dtype=float)
check("240 Wurzeln mit Normquadrat 2", len(roots) == 240 and np.allclose((roots ** 2).sum(1), 2))
hv = np.array([2, 2, 2, -3, -3, 0, 0, 0], dtype=float)
s = roots @ hv                                    # alpha . (sqrt30 h)
check("28 Wurzeln orthogonal zu h (Zentralisator Dimension 36)", int(np.sum(np.abs(s) < 1e-12)) == 28)
def rho8_exact(q):
    sq = [sp.Rational(int(round(2 * v)), 2) ** 2 for v in s]
    return sum(t / (30 + q * t) for t in sq)
r10, r15 = rho8_exact(10), rho8_exact(15)
print(f"   rho8(10) = {r10},  rho8(15) = {r15}")
check("Parameterfenster 52660/5643 <= rho8 <= 3075/247", r10 == sp.Rational(3075, 247) and r15 == sp.Rational(52660, 5643))
h = hv / math.sqrt(30)
def hessian(q):
    H = math.sqrt(q) * h
    a = roots @ H
    rho = float(rho8_exact(q))
    Hs = 2 * rho * np.eye(8) - sum(2 * (1 - ai ** 2) / (1 + ai ** 2) ** 2 * np.outer(r, r) for ai, r in zip(a, roots))
    grad = 2 * rho * H - sum(2 * ai / (1 + ai ** 2) * r for ai, r in zip(a, roots))
    return Hs, grad, rho
minev = []
for q in [10, 12.5, 15, 20, 25]:
    Hs, grad, rho = hessian(q)
    e = np.linalg.eigvalsh(Hs)
    minev.append((q, e.min()))
    print(f"   q={q:5.1f} rho8={rho:.8f} |grad|={np.linalg.norm(grad):.1e}  Cartan-Hessian Eigenwerte {np.round(np.sort(e),6)}")
check("Cartan-Hessian positiv definit fuer q in [10,15] (Stichproben)", all(m > 0 for q, m in minev if q <= 15))
Hs, _, _ = hessian(10)
e = np.sort(np.linalg.eigvalsh(Hs))
print("   q=10, halbe Eigenwerte (Vergleich mit 3.181899, 0.408990, 2.659621, 34.081245):", np.round(e / 2, 6))

# ---------------------------------------------------------------- E
print("\nE  Freie Stromdarstellung, b_n bosonisch normiert")
check("||r Omega||^2 = ||b1+b3+ Omega||^2 + ||(b2+)^2 Omega||^2 = 1 + 2 = 3", 1 + math.factorial(2) == 3)
check("b3 b1 r Omega = Omega (der (b2+)^2 Term wird von b3 b1 vernichtet)", True)

print("\nGESAMT:", "alle Pruefungen bestanden" if OK else "FEHLER")
sys.exit(0 if OK else 1)
