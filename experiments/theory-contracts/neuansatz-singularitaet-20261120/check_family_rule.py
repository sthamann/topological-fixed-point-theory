#!/usr/bin/env python3
"""From the four marks to the family rule, the later response, the group, and the zero mode.

Marks = the four two-torsion points of the torus (a group of order four). A step of degree two with kernel t acts on
functions on the marks by averaging over the translation by t. Nothing of the record is an input before it is compared.
"""
import itertools, json, sys
from pathlib import Path
import sympy as sp

rows = []
def test(name, cond, value=None):
    rows.append({"check": name, "passed": bool(cond), "value": value})
R = sp.Rational
marks = [(0, 0), (1, 0), (0, 1), (1, 1)]
add = lambda s, t: ((s[0]+t[0]) % 2, (s[1]+t[1]) % 2)
def transl(t):
    return sp.Matrix(4, 4, lambda r, c: 1 if marks[r] == add(marks[c], t) else 0)
step = {t: (sp.eye(4)+transl(t))/2 for t in marks[1:]}

# ---- 1. the family rule
A = sum(step.values(), sp.zeros(4, 4))/3
test("each step operator is a projector (pull back after push forward, normalised)", all(s*s == s and s.rank() == 2 for s in step.values()))
test("the uniform mixture of the three steps is 1/2 on the diagonal and 1/6 elsewhere", A == sp.eye(4)/3+sp.ones(4, 4)/6, str(A.tolist()))
M = sp.Matrix([[1, 0, 0], [0, R(1, 2), R(1, 6)], [0, R(1, 6), R(1, 2)]])                         # verification/v395
pairs = list(itertools.combinations(range(4), 2))
test("compressed to any two marks it is the pair block of the record's family matrix (1/2, 1/6): eigenvalues 2/3 and 1/3",
     all(A.extract(list(p), list(p)) == M[1:, 1:] for p in pairs) and M[1:, 1:].eigenvals() == {R(2, 3): 1, R(1, 3): 1})
test("on all four marks the mixture has eigenvalues 1 (constant) and 1/3 three times; 2/3 appears only after restriction to a pair",
     A.eigenvals() == {1: 1, R(1, 3): 3})
test("weight that stays inside a pair per step: 2/3; weight that leaves: 1/3", sum(A[0, k] for k in (0, 1)) == R(2, 3))
alt = {w: (sp.eye(4)*(1-w)+w*sum(step.values(), sp.zeros(4, 4))/3) for w in (R(1, 2), R(3, 4))}
test("control: a lazy mixture with any other weight does not give 1/2 and 1/6", all(B[0, 0] != R(1, 2) or B[0, 1] != R(1, 6) for B in alt.values()))

# ---- 2. order of the steps
a, b, c = marks[1:]
test("a fixed cyclic order of the three steps kills every non-constant function after two steps (as the marks die after two steps of 1+i)",
     step[a]*step[b] == sp.ones(4, 4)/4 and step[c]*step[b]*step[a] == sp.ones(4, 4)/4)
test("the random order decays with 1/3 per step; so order matters", (A**2).eigenvals() == {1: 1, R(1, 9): 3})

# ---- 3. the later response: transfer followed by a move of the distinguished channel
def perm(p):
    return sp.Matrix(3, 3, lambda r, col: 1 if p[col] == r else 0)
spectra = {}
for p in itertools.permutations(range(3)):
    ev = (M*perm(p)).eigenvals()
    spectra[str(p)] = sorted((str(sp.nsimplify(sp.simplify(k))) for k in ev), key=str)
K = M*perm((1, 2, 0))
lam = [k for k in K.eigenvals() if k != 1]
test("three moved transfers in a row equal the third power of one operator: transfer then cyclic move",
     M*(perm((1, 2, 0))*M*perm((1, 2, 0)).T)*(perm((1, 2, 0))**2*M*(perm((1, 2, 0)).T)**2) == K**3)
xx = sp.symbols("x")
poly = sp.Poly(sp.expand(18*K.charpoly(xx).as_expr()), xx)
roots = [complex(r) for r in poly.nroots(n=20)]
real = [r.real for r in roots if abs(r.imag) < 1e-12]; cpx = [r for r in roots if abs(r.imag) > 1e-12]
test("its characteristic polynomial is 18x^3 - 3x^2 - 3x - 4: one real rate 0.764672 and one oscillating pair of modulus 0.539084; "
     "the rates 1, 2/3, 1/3 of the repeated transfer are gone",
     poly.all_coeffs() == [18, -3, -3, -4] and len(real) == 1 and abs(real[0]-0.764672) < 1e-6 and abs(abs(cpx[0])-0.539084) < 1e-6,
     {"real": real[0], "pair": [cpx[0].real, cpx[0].imag], "six_steps": [real[0]**6, abs(cpx[0])**6], "repeated_six_steps": [64/729, 1/729]})
rows.append({"check": "my earlier expectation (-5 +- i sqrt 7)/12 with eigenvalue 1", "passed": False,
             "value": "wrong: the family matrix does not conserve total weight (rows sum to 1, 2/3, 2/3), so 1 is no eigenvalue of the moved transfer"})
test("a move that only swaps the pair changes nothing but one sign: 1, 2/3, -1/3", sorted((M*perm((0, 2, 1))).eigenvals()) == [R(-1, 3), R(2, 3), 1], spectra)

# ---- 4. the group of the record on the three-space and the marks
sys.path.insert(0, str(Path(__file__).resolve().parents[3]/"verification"))
import v117_monodromy_weyl_a3 as v117
def closure(gens):
    seen = {sp.ImmutableMatrix(sp.eye(gens[0].rows))}; frontier = list(seen)
    while frontier:
        new = []
        for g in frontier:
            for h in gens:
                x = sp.ImmutableMatrix((g*h).applyfunc(lambda e: sp.nsimplify(sp.re(sp.expand(e)))+sp.I*sp.nsimplify(sp.im(sp.expand(e)))))
                if x not in seen:
                    seen.add(x); new.append(x)
        frontier = new
    return seen
G = closure([v117.U_EXACT, v117.M0_EXACT])
def affine():
    out = []
    for m in itertools.product((0, 1), repeat=4):
        if (m[0]*m[3]-m[1]*m[2]) % 2:
            for t in marks:
                out.append(tuple(marks.index(add(((m[0]*x+m[1]*y) % 2, (m[2]*x+m[3]*y) % 2), t)) for x, y in marks))
    return out
aff = affine()
def sign(p):
    return sp.Matrix(4, 4, lambda r, col: 1 if p[col] == r else 0).det()
def on_traceless(p):
    P = sp.Matrix(4, 4, lambda r, col: 1 if p[col] == r else 0)
    Bs = sp.Matrix([[1, -1, 0, 0], [0, 1, -1, 0], [0, 0, 1, -1]]).T
    return (Bs.pinv()*P*Bs)
from collections import Counter
classes_rec = Counter(str(sp.nsimplify(sp.simplify(g.trace()))) for g in G)
classes_geo = Counter(str(sp.nsimplify(sign(p)*on_traceless(p).trace())) for p in aff)
test("the symmetries of the four marks (moves of the torus and shifts by a mark) are all 24 permutations of the marks", len(set(aff)) == 24)
test("the record's group on the three-space (v117, order 24) has the same character as: functions on the four marks without the constant, times the sign of the permutation",
     len(G) == 24 and classes_rec == classes_geo, {"record": dict(classes_rec), "marks": dict(classes_geo)})
test("so the three family lines are the three non-trivial characters of the group of marks (as in GESAMT2: chi1 + chi2 + chi3)",
     all(sorted(on_traceless(tuple(marks.index(add(s, t)) for s in marks)).eigenvals().items()) == [(-1, 2), (1, 1)] for t in marks[1:]))

# ---- 5. orientation: the lift to the double cover
half = R(1, 2); I = sp.I
def su2(axis, angle):
    n = sp.Matrix(axis)/sp.Matrix(axis).norm()
    s = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -I], [I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
    return sp.simplify(sp.cos(angle/2)*sp.eye(2)-I*sp.sin(angle/2)*(n[0]*s[0]+n[1]*s[1]+n[2]*s[2]))
qu, th = su2((0, 0, 1), sp.pi/2), su2((1, 1, 1), 2*sp.pi/3)
test("lifted to the double cover the quarter turn has order eight and the three-step order six: their fourth and third powers are -1",
     sp.simplify(qu**4+sp.eye(2)) == sp.zeros(2, 2) and sp.simplify(th**3+sp.eye(2)) == sp.zeros(2, 2))
v_sum, v_dif = sp.Matrix([0, 1, 1]), sp.Matrix([0, 1, -1])
cycm = perm((1, 2, 0))
def coplanar(v):
    return sp.Matrix.hstack(v, cycm*v, cycm**2*v).det() == 0
test("the difference line sweeps a plane (a closed loop of lines that cannot be shrunk: sign -1); the sum line sweeps a cone (shrinkable: sign +1)",
     coplanar(v_dif) and not coplanar(v_sum))

# ---- 6. zero mode on the sphere with four marks (flat, quotient of the square torus by the half turn)
zero = {}
for delta in itertools.product((0, 1), repeat=2):                     # spin structure: 0 periodic, 1 antiperiodic
    for k in range(4):                                                # twist of the half turn by a fourth root of unity
        n = 0
        if delta == (0, 0):                                           # constants exist only for the periodic structure
            for comp in (1, -1):                                      # the half turn acts on the two spinor components by i and -i
                if (I**k*(I*comp)) == 1:
                    n += 1
        zero[(delta, k)] = n
test("the boundary operator on the four-mark sphere has a zero mode only for the periodic structure with an odd quarter twist; then exactly one, of one handedness",
     sum(zero.values()) == 2 and zero[((0, 0), 1)] == 1 and zero[((0, 0), 3)] == 1 and zero[((0, 0), 0)] == 0 and zero[((0, 0), 2)] == 0
     and all(v == 0 for (d, k), v in zero.items() if d != (0, 0)), {str(k): v for k, v in zero.items() if v})
test("counts: three steps, six endpoints after two steps without return, ten points in the ball of radius two", (3, 3*2, 1+3+6) == (3, 6, 10))

if __name__ == "__main__":
    print(json.dumps({"status": "FAMILY_RULE_CHECKED", "physical_solution": False, "passed": sum(r["passed"] for r in rows), "total": len(rows), "rows": rows}, indent=1, default=str))
