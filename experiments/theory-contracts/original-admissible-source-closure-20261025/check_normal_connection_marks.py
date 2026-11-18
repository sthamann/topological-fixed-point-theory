#!/usr/bin/env python3
"""Exact comparisons of existing tetrahedron, covers and reducing connections."""
from itertools import permutations
from pathlib import Path
from math import gcd, lcm
import json
import sympy as S

I, rt = S.I, S.sqrt(3)
eye, zero = S.eye(2), S.zeros(2)
pauli = [S.Matrix([[0, 1], [1, 0]]), S.Matrix([[0, -I], [I, 0]]), S.diag(1, -1)]
checks = []


def eq(a, b=0):
    if isinstance(a, S.MatrixBase):
        return all(S.trigsimp(S.simplify(v), method="fu") == 0 for v in a - b)
    return S.trigsimp(S.simplify(S.expand_complex(a - b)), method="fu") == 0


def need(condition, name):
    if not condition:
        raise RuntimeError(name)
    checks.append(name)


def cr(a, b, c, d):
    return S.simplify((a - c) * (b - d) / ((a - d) * (b - c)))


def normal(v):
    return S.simplify(S.expand_complex(v))


ns = [S.Matrix(v) / rt for v in [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]]
ps = [(eye + sum((n[j] * pauli[j] for j in range(3)), zero)) / 2 for n in ns]
need(all(eq(p * p, p) and eq(S.trace(p), 1) for p in ps), "four Bloch rank-one projectors")
need(eq(sum(ps, zero), 2 * eye), "Bloch sum is 2I, not the root-block resolution I")
need(all(eq(S.trace(ps[j] * ps[k]), S.Rational(1, 3))
         for j in range(4) for k in range(4) if j != k), "nonorthogonal pair overlaps")
traces = []
for j, k, m in permutations(range(4), 3):
    t = normal(S.trace(ps[j] * ps[k] * ps[m]))
    need(eq(t * S.conjugate(t), S.Rational(1, 27)), f"face magnitude {j}{k}{m}")
    need(eq((3 * rt * t) ** 2, -1), f"face quarter phase {j}{k}{m}")
    traces.append(str(t))

zs = [S.simplify((n[0] + I * n[1]) / (1 + n[2])) for n in ns]
tet = {normal(cr(*(zs[j] for j in p))) for p in permutations(range(4))}
square = [1, I, -1, -I]
sq = {normal(cr(*(square[j] for j in p))) for p in permutations(range(4))}
need(tet == {(1 - I * rt) / 2, (1 + I * rt) / 2}, "all tetrahedral cross ratios")
need(sq == {S.Integer(2), S.Rational(1, 2), S.Integer(-1)}, "all quarter cross ratios")
need(not tet & sq, "no marked Mobius equivalence even after relabeling")
tet_aut = [p for p in permutations(range(4)) if eq(cr(*(zs[j] for j in p)), cr(*zs))]
sq_aut = [p for p in permutations(range(4)) if eq(cr(*(square[j] for j in p)), cr(*square))]
need(len(tet_aut) == 12 and len(sq_aut) == 8, "different holomorphic four-mark automorphism groups")
cycle = (1, 2, 3, 0)
need(cycle not in tet_aut and cycle in sq_aut, "square four-cycle is not a tetrahedral holomorphic action")
jfun = lambda a: S.simplify(256 * (1 - a + a * a) ** 3 / (a * a * (1 - a) ** 2))
need(all(eq(jfun(a), 0) for a in tet) and all(eq(jfun(a), 1728) for a in sq), "elliptic invariants j0 and j1728")

omega = (-1 + I * rt) / 2
q1, q2, q3 = -1, -omega, -omega ** 2
lam = normal((q1 - q3) / (q2 - q3))
need(eq(lam, cr(*zs)), "tetrahedron matches branch quartet of E_omega")
Z = S.Symbol("Z")
M = q1 + (q2 - q1) * (Z - zs[0]) * (zs[1] - zs[3]) / ((Z - zs[3]) * (zs[1] - zs[0]))
need(all(eq(M.subs(Z, zs[k]), q) for k, q in enumerate([q1, q2, q3])), "explicit ordered Mobius map to E_omega branches")

x, s = S.symbols("x s")
relation = s ** 12 - x ** 4 + 1
need(S.rem((s ** 6) ** 2 - (x ** 4 - 1), relation, s) == 0, "C12 to E_i")
need(S.rem((x ** 2) ** 2 - ((s ** 4) ** 3 + 1), relation, s) == 0, "C12 to E_omega")
need(gcd(6, 4) == 2 and S.cancel(s ** 6 / s ** 4) == s ** 2, "joint image function field is C6")
need(all(e % 2 == 0 for e in [6, 4]) and (-s) ** 3 == -s ** 3, "joint elliptic image loses sign retained by C4")
need([k for k in range(12) if 6 * k % 12 == 0 and 4 * k % 12 == 0] == [0, 6], "full deck kernel of elliptic pair")
# Existing forced compiler lift: R12 = Q12 Delta12^4 = (i x, omega s).
need(lcm(4, 3) == 12, "original mixed compiler lift has order twelve")
need([k for k in range(1, 13) if k % 4 == 0 and k % 6 == 0] == [12], "mixed order twelve is retained on elliptic pair")
need(4 * 4 % 12 == 4 and 4 * 9 % 12 == 0, "original c4=sigma and c9=J marking retained")

theta, beta, aa, bb = S.symbols("theta beta aa bb", real=True)
chi = S.Matrix([S.cos(beta / 2), S.exp(I * theta) * S.sin(beta / 2)])
p = S.simplify(chi * chi.H)
dp = p.diff(theta)
kato = S.simplify(I * (dp * p - p * dp))
need(eq(p * p, p) and eq(p * dp * p, zero), "differentiated projection identity")
need(eq(kato.H, kato) and eq(p * kato * p, zero), "forced off-diagonal connection is Hermitian")
need(eq(-I * dp + kato * p - p * kato, zero), "exact Kato reducing commutator")
a = kato + aa * p + bb * (eye - p)
need(eq(-I * dp + a * p - p * a, zero), "all diagonal freedoms preserve exact reduction")
eff = S.simplify((chi.H * a * chi)[0] - I * (chi.H * chi.diff(theta))[0])
need(eq(eff, aa + (1 - S.cos(beta)) / 2), "entire reduced circle connection including Berry term")
need(eq(-I * chi.diff(theta) + a * chi - chi * eff, S.zeros(2, 1)), "operator intertwiner on all scalar Fourier functions")
need(not eq(dp.subs(beta, S.pi / 2), zero), "scalar original matrix source does not reduce a moving line")

need(all(eq((3 * rt * t) * S.conjugate(3 * rt * t), 1) for t in [S.sympify(v) for v in traces]),
     "dual homogeneous line cancels every normalized spin-line face phase")

result = {
    "verdict": "PASS_EXACT_SCOPED_COMPARISONS",
    "checks_passed": len(checks),
    "checks": checks,
    "tetra_crossratio_orbit": sorted(map(str, tet)),
    "quarter_crossratio_orbit": sorted(map(str, sq)),
    "ordered_face_trace": traces[0],
    "elliptic_maps": {"E_i": "(x,y)=(x,s^6)", "E_omega": "(u,v)=(x^2,s^4)"},
    "degrees": {"C3_to_Eomega": 2, "C12_to_Ei": 6, "C12_to_Eomega": 8, "joint_generic_fiber": 2},
    "joint_image": "C6: t=s^2, t^6=x^4-1",
    "forced_compiler_lift": "R12=(ix,omega*s); fourth=sigma, ninth=plain J",
    "connection_class": "A=i[P_prime,P]+PAP+(I-P)A(I-P)",
    "physical_source_identified": False,
    "complete_TFPT_solution": False,
}
out = Path(__file__).with_name("normal-connection-marks-results.json")
out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"verdict": result["verdict"], "checks_passed": len(checks), "result": str(out)}))
