#!/usr/bin/env python3
"""Independent recomputation of the checkable statements of two external reports (pasted 20 November 2026).

Nothing is taken from the reports' own files. Statements that need objects not reconstructed here are listed as not
checked. A failed statement is recorded, it does not stop the run.
"""
import itertools, json, sys
from pathlib import Path
import numpy as np
import sympy as sp

rows = []
def rec(report, name, verdict, value=None):
    rows.append({"report": report, "statement": name, "verdict": verdict, "value": value})
def test(report, name, cond, value=None):
    rec(report, name, "confirmed" if cond else "NOT confirmed", value)

# ================= report B: cover from a loss of distinguishability =================
def lap(k12, k23, k31):
    return sp.Matrix([[k12+k31, -k12, -k31], [-k12, k12+k23, -k23], [-k31, -k23, k23+k31]])
def top(M):
    lam, mult, vecs = max(M.eigenvects(), key=lambda t: t[0])
    v = vecs[0]/vecs[0].norm(); return v*v.T
P = [top(lap(*k)) for k in ((2, 1, 1), (1, 2, 1), (1, 1, 2))]
test("B", "three couplings (2,1,1),(1,2,1),(1,1,2): tr(P1 P2 P3) = -1/8 for the strongest mode", sp.simplify((P[0]*P[1]*P[2]).trace()) == sp.Rational(-1, 8))
M = sp.Matrix([[1, 0, 0], [0, sp.Rational(1, 2), sp.Rational(1, 6)], [0, sp.Rational(1, 6), sp.Rational(1, 2)]])   # verification/v395
cyc = sp.Matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
table = {}
for lam, mult, vecs in M.eigenvects():
    v = vecs[0]/vecs[0].norm(); p = v*v.T
    table[str(lam)] = (sp.simplify((p*p*p).trace()), sp.simplify((p*(cyc*p*cyc.T)*(cyc*cyc*p*cyc.T*cyc.T)).trace()))
test("B", "v395 matrix, sectors 1, 2/3, 1/3: same assignment gives 1, 1, 1; cyclic displacement gives 0, +1/8, -1/8",
     table == {"1": (1, 0), "2/3": (1, sp.Rational(1, 8)), "1/3": (1, sp.Rational(-1, 8))}, {k: [str(x) for x in v] for k, v in table.items()})
rec("B", "the cyclic displacement of the distinguished channel is produced by the theory", "not shown (the report says so itself)")
s, V = sp.symbols("s V")
vc = sp.symbols("vc"); sol = sp.solve(sp.Eq(9*(V-vc), s*vc), vc)[0]
per_port = sp.simplify(3*(V-sol)/V)
test("B", "star with a storing centre: common drive gives 3s/(9+s) per port, the triangle gives 0, equal at s = 0",
     sp.simplify(per_port-3*s/(9+s)) == 0, str(per_port))
star = sp.Matrix([[3, 0, 0, -3], [0, 3, 0, -3], [0, 0, 3, -3], [-3, -3, -3, 9]])
schur = star[:3, :3]-star[:3, 3]*star[3, :3]/star[3, 3]
test("B", "star (conductance 3) and triangle (conductance 1) have the same static boundary response", schur == lap(1, 1, 1))
def returns(n, order):
    A = np.zeros((n+1, n+1), dtype=object)
    for i in range(n):
        A[i, i+1] = A[i+1, i] = 1
    out, X = [], np.eye(n+1, dtype=object)
    for _ in range(order):
        X = X.dot(A); out.append(X[0, 0])
    return out
a, b = returns(8, 20), returns(9, 20)
first = next(i+1 for i, (x, y) in enumerate(zip(a, b)) if x != y)
test("B", "two chains (8 and 9 links): return answers agree through order 17 and differ by one at order 18", first == 18 and b[17]-a[17] == 1, first)
perms = list(itertools.permutations(range(4)))
def on3(p, pos):
    out = []
    for x in range(8):
        bits = [(x >> 2) & 1, (x >> 1) & 1, x & 1]
        y = p[bits[pos]*2+bits[pos+1]]; bits[pos], bits[pos+1] = y >> 1, y & 1
        out.append(bits[0]*4+bits[1]*2+bits[2])
    return out
comp = lambda f, g: [f[g[x]] for x in range(8)]
good = [p for p in perms if comp(on3(p, 0), comp(on3(p, 1), on3(p, 0))) == comp(on3(p, 1), comp(on3(p, 0), on3(p, 1)))]
swap_like = [p for p in good if all(((p[x] >> 1) ^ (x & 1)) == ((p[0] >> 1) ^ 0) and ((p[x] & 1) ^ (x >> 1)) == (p[0] & 1) for x in range(4))]
test("B", "of the 24 reversible rules on two bits exactly five give the same result in both exchange orders on three bits",
     len(good) == 5, {"count": len(good), "rules": [list(p) for p in good]})
test("B", "they are the identity and four swaps with fixed bit flips", (0, 1, 2, 3) in good and len(swap_like) == 4, len(swap_like))
phi = sp.symbols("phi")
R = sp.Matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, phi]])
I2 = sp.eye(2); K = sp.kronecker_product
yb = sp.simplify(K(R, I2)*K(I2, R)*K(R, I2)-K(I2, R)*K(R, I2)*K(I2, R))
test("B", "quantum exchange with a phase on the joint occupation 11 is consistent for every phase; exchange twice = identity leaves +1 and -1",
     yb == sp.zeros(8, 8) and sp.solve(sp.Eq(phi**2, 1), phi) == [-1, 1])

# ================= report A: joint quantum process =================
T = lambda i, j: sp.Matrix(6, 6, lambda r, c: 1 if (r, c) in ((i, j), (j, i)) or (r == c and r not in (i, j)) else 0)
Q = (sp.eye(6)-T(0, 1)*T(2, 3)*T(4, 5))/2
basis = sp.Matrix([[1, -1, 0, 0, 0, 0], [0, 0, 1, -1, 0, 0], [0, 0, 0, 0, 1, -1]]).T/sp.sqrt(2)
test("A", "Q = (1 - T01 T23 T45)/2 projects the sum-zero space onto the three differences", Q*Q == Q and Q.rank() == 3 and Q*basis == basis)
u = sp.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]]); m = sp.Matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
def restricted(W):
    return basis.T*W*basis if (Q*W*basis == W*basis) else None
found = {}
for name, word, target in (("U6", [T(0, 3), T(0, 1), T(0, 2)], u), ("M6", [T(0, 2), T(0, 4), T(1, 3), T(1, 5)], m)):
    for label, w in (("as written", word), ("reversed", word[::-1])):
        W = sp.eye(6)
        for t in w:
            W = W*t
        r = restricted(W)
        found[f"{name} {label}"] = "leaves the three-space" if r is None else ("equals target" if r == target else str(r.tolist()))
test("A", "U6 = T03 T01 T02 and M6 = T02 T04 T13 T15 act on that three-space as the stated quarter turn u and three-step m",
     any(found[f"U6 {l}"] == "equals target" for l in ("as written", "reversed")) and any(found[f"M6 {l}"] == "equals target" for l in ("as written", "reversed")), found)
sys.path.insert(0, str(Path(__file__).resolve().parents[3]/"verification"))
import v117_monodromy_weyl_a3 as v117
def group(gens):
    seen = {sp.ImmutableMatrix(sp.eye(3))}; frontier = list(seen)
    while frontier:
        new = []
        for g in frontier:
            for h in gens:
                x = sp.ImmutableMatrix((g*h).applyfunc(sp.simplify))
                if x not in seen:
                    seen.add(x); new.append(x)
        frontier = new
    return seen
G = group([u, m])
unk = sp.symbols("w0:9"); W = sp.Matrix(3, 3, unk)
eqs = list(u*W-W*v117.U_EXACT)+list(m*W-W*v117.M0_EXACT)
ns = sp.Matrix([[sp.expand(e).coeff(x) for x in unk] for e in eqs]).nullspace()
Wsol = sp.Matrix(3, 3, list(ns[0])) if ns else None
uni = Wsol is not None and sp.simplify(Wsol.H*Wsol-(Wsol.H*Wsol)[0, 0]*sp.eye(3)) == sp.zeros(3, 3)
test("A", "u and m generate a group of order 24, and one unitary W obeys uW = W U_v117 and mW = W M0_v117",
     len(G) == 24 and len(ns) == 1 and uni, {"order": len(G), "solutions": len(ns)})
v4 = [sp.eye(6), T(0, 1)*T(2, 3), T(2, 3)*T(4, 5), T(0, 1)*T(4, 5)]
h5 = sp.Matrix([[1, -1, 0, 0, 0, 0], [0, 0, 1, -1, 0, 0], [0, 0, 0, 0, 1, -1], [1, 1, -1, -1, 0, 0], [1, 1, 1, 1, -2, -2]]).T
chars = [tuple(int((h5.pinv()*g*h5)[k, k]) for g in v4) for k in range(5)]
test("A", "under the Klein four group: H5 = 2 trivial + three different characters, the three-space carries the three",
     sorted(chars) == sorted([(1, 1, 1, 1)]*2+[(1, -1, -1, 1), (1, -1, 1, -1), (1, 1, -1, -1)]) and all(sum(c) == 0 for c in chars[:3]), chars)
B = sp.Matrix([[13, 1, 4], [1, 13, 4], [4, 4, 10]])/18
test("A", "B has eigenvalues 1, 2/3, 1/3", B.eigenvals() == {1: 1, sp.Rational(2, 3): 1, sp.Rational(1, 3): 1})
test("A", "rates: N = P1 + 3/4 P5 + 11/36 P9 + 1/4 P10 gives exp(-t/4), exp(-25t/36), exp(-3t/4); ratio of the two orders exp(4t/9); values at t = 2",
     sp.Rational(25, 36)-sp.Rational(1, 4) == sp.Rational(4, 9)
     and abs((32-2*6**0.5)/375*np.exp(-0.5)-0.0438336) < 1e-7 and abs((32-2*6**0.5)/375*np.exp(-25/18)-0.0180205) < 1e-7)
rate = {1: -(np.log(2)/2)*1-np.log(9/8)/2, 2: -(np.log(2)/2)*3-np.log(9/8)/2}
test("A", "the stated generator has rates log(2/3) on the rotation sector and log(1/3) on the rank-two sector, so its population action is log B",
     abs(rate[1]-np.log(2/3)) < 1e-12 and abs(rate[2]-np.log(1/3)) < 1e-12)

# covariant channels on 3 x 3 matrices: sectors 1, E (2), T2 (3, symmetric off-diagonal), T1 (3, antisymmetric)
def unit(i, j, d):
    e = np.zeros((d, d)); e[i, j] = 1; return e
def choi(phi, d):
    return sum(np.kron(phi(unit(i, j, d)), unit(i, j, d)) for i in range(d) for j in range(d))
def family(eta, logs=False):
    f = (lambda x: np.log(x)) if logs else (lambda x: x)
    def phi(X):
        tr = np.trace(X)/3*np.eye(3); sym = (X+X.T)/2; anti = (X-X.T)/2
        diag = np.diag(np.diag(sym))-tr; off = sym-np.diag(np.diag(sym))
        return (0 if logs else 1)*tr+f(1/3)*diag+f(eta)*off+f(2/3)*anti
    return phi
def cp(eta): return np.linalg.eigvalsh(choi(family(eta), 3)).min() > -1e-12
def ccp(eta):
    om = np.eye(3).reshape(9)/np.sqrt(3); Pp = np.eye(9)-np.outer(om, om)
    return np.linalg.eigvalsh(Pp@choi(family(eta, True), 3)@Pp).min() > -1e-12
def edges(fn, lo=0.05, hi=0.95, n=18001):
    xs = [x for x in np.linspace(lo, hi, n) if fn(x)]
    return (round(float(xs[0]), 4), round(float(xs[-1]), 4)) if xs else None
e_cp, e_ccp = edges(cp), edges(ccp)
test("A", "complete positivity of the family channel forces 2/9 <= eta <= 4/9", e_cp is not None and abs(e_cp[0]-2/9) < 2e-4 and abs(e_cp[1]-4/9) < 2e-4, e_cp)
test("A", "a continuous covariant evolution narrows this to 0.3205 <= eta <= 0.3467", e_ccp is not None and abs(e_ccp[0]-0.3205) < 2e-4 and abs(e_ccp[1]-0.3467) < 2e-4, e_ccp)

# S6 covariant channel on the five-space
h = np.array(sp.GramSchmidt([h5[:, k] for k in range(5)], True)).astype(float).reshape(5, 6).T   # 6 x 5 isometry
Dx = [h.T@np.diag(np.eye(6)[k])@h for k in range(6)]                                     # compressions of the six diagonal units
D5 = [d-np.trace(d)/5*np.eye(5) for d in Dx]
G5 = np.array([d.reshape(25) for d in D5]); _, sv, vt = np.linalg.svd(G5); B5 = vt[:5]
P5 = B5.T@B5
def chanN(vals):
    def phi(X):
        tr = np.trace(X)/5*np.eye(5); sym = (X+X.T)/2-tr; anti = (X-X.T)/2
        s5 = (P5@sym.reshape(25)).reshape(5, 5); s9 = sym-s5
        return vals[0]*tr+vals[1]*s5+vals[2]*s9+vals[3]*anti
    return phi
N = (1, 3/4, 11/36, 1/4)
minN = np.linalg.eigvalsh(choi(chanN(N), 5)).min()
om = np.eye(5).reshape(25)/np.sqrt(5); Pp = np.eye(25)-np.outer(om, om)
minlog = np.linalg.eigvalsh(Pp@choi(chanN(tuple(np.log(N))), 5)@Pp).min()
test("A", "N = P1 + 3/4 P5 + 11/36 P9 + 1/4 P10 is a quantum channel on the five-space (sectors 1 + 5 + 9 + 10)", minN > -1e-12 and round(float(sv[4]), 6) > 0 and sv[5] < 1e-12, float(minN))
test("A", "the covariant logarithm of N fails the positivity test with eigenvalue about -0.087696", abs(minlog-(-0.087696)) < 2e-6, float(minlog))

for name in ("norm 10/9 on the outer fourth sector; M >= 2/405; M = 1 on the code",
             "environment isometry W_env with (1 x W_env) S_evt = V_triangle",
             "three settings 1, R T01, R T04 give rank 25; d_R = (32 - 2 sqrt 6)/375",
             "native first return 59/84, 11/36, 1/4 on the three-space",
             "Koide part: static scaling, coherent source, interval -0.6950927 < h < -0.6249998"):
    rec("A", name, "not checked (needs the report's own operators; the report's files were not available here)")

print(json.dumps({"status": "REPORTS_CHECKED", "physical_solution": False, "rows": rows}, indent=1, default=str))
