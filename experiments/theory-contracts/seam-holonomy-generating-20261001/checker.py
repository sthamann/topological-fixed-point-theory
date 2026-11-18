#!/usr/bin/env python3
"""seam-holonomy-generating-20261001 -- exact checker (experiments firewall).

Question: are T, R, B_Y, R_Y and a family (flavor) operator Taylor coefficients
of ONE object, the response of the original local pair products
B_beta = 1/2 :(beta.h)^2: to a Cartan holonomy z on the seam,

    W_S(z) = sum_{beta in S} exp(i z.beta) B_beta   <->   A_S(z) = sum exp(i z.beta) beta beta^T ?

and what does the E8 design property say about which operators the closed
source can produce?  All roots are the ORIGINAL Chevalley roots of
tfpt_explorer (doubled integer coordinates; beta_std = r/2).

Not claimed: a physical seam holonomy, the P1 origin of the selection S=48,
a 4D field map, a mass or mixing prediction.  Run from the repository root:
    PYTHONPATH=. python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/checker.py
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, product
from pathlib import Path

import numpy as np
import sympy as sp

from tfpt_explorer.neutral_source_response import _actual_current_source, _native_neutral_branch
from tfpt_explorer.source_charge_response import hypercharge_vector
from tfpt_explorer.source_majorana_pairs import _native_pairs

HERE = Path(__file__).resolve().parent
RESULTS: list[dict] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append({"check": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


# ---------------------------------------------------------------- roots
src = _actual_current_source()
branch = _native_neutral_branch(0)
R240 = np.array([list(map(int, r)) for r in src["chevalley"].roots])
R64 = np.array([list(map(int, r)) for r in branch["spinor"]])
R48 = np.array([list(map(int, r)) for r in branch["retained"]])
check("C0.root counts 240/64/48", (len(R240), len(R64), len(R48)) == (240, 64, 48))


def pairings(idx):
    if not idx:
        yield []
        return
    for j in range(1, len(idx)):
        rest = idx[1:j] + idx[j + 1:]
        for p in pairings(rest):
            yield [(idx[0], idx[j])] + p


def isotropic(rank: int) -> np.ndarray:
    d = np.eye(8, dtype=np.int64)
    out = np.zeros((8,) * rank, dtype=np.int64)
    letters = "abcdefgh"[:rank]
    for p in pairings(list(range(rank))):
        out += np.einsum(",".join(letters[i] + letters[j] for i, j in p) + "->" + letters, *([d] * len(p)))
    return out


def moment(S: np.ndarray, rank: int) -> np.ndarray:
    letters = "abcdefgh"[:rank]
    return np.einsum(",".join("z" + l for l in letters) + "->" + letters, *([S] * rank))


def is_isotropic(S: np.ndarray, rank: int) -> bool:
    M, I = moment(S, rank), isotropic(rank)
    c = F(int(M.flat[0]), int(I.flat[0]))
    return bool(np.array_equal(M * c.denominator, I * c.numerator))


# ------------------------------------------- C1: design theorem on the closure
# Full tensor identities (a proof, not a sample): the even moments of rank 2,4,6
# of the 240 original roots are isotropic and all odd moments vanish.
iso240 = {k: is_isotropic(R240, k) for k in (2, 4, 6)}
odd240 = {k: not moment(R240, k).any() for k in (1, 3, 5)}
check("C1.240 even moments rank2/4/6 isotropic (7-design)", all(iso240.values()), str(iso240))
check("C1.240 odd moments rank1/3/5 vanish", all(odd240.values()), str(odd240))
# rank 8 is NOT isotropic: first W(E8)-specific (degree-8 Casimir) order.
rng = np.random.default_rng(20261001)
z, y8 = rng.integers(-5, 6, 8), rng.integers(-5, 6, 8)
lhs = sum(int(z @ b) ** 6 * int(y8 @ b) ** 2 for b in R240)
zz, yy, zy = int(z @ z), int(y8 @ y8), int(z @ y8)
pred = F(240 * 8 ** 4, 8 * 10 * 12 * 14) * (15 * zz ** 3 * yy + 90 * zz ** 2 * zy ** 2)
check("C1.240 rank8 NOT isotropic (exact counterexample point)", lhs != pred, f"{lhs} vs {pred}")
# Consequence in std coordinates: sum_240 (z.b)^2 b b^T = 12(|z|^2 I + 2 z z^T), hence
#   W_240 = 60 T  - (1/2)(10 T + 12 :J_z J_z:) |_{z=y} eps^2 + O(eps^4, isotropic) ...
zs = sp.Matrix(sp.symbols("z0:8"))
B = [sp.Matrix(r) / 2 for r in R240]
c0 = sum((b * b.T for b in B), sp.zeros(8))
yv = hypercharge_vector()
c2y = sum(((yv.T * b)[0] ** 2 * b * b.T for b in B), sp.zeros(8))
check("C1.W_240(0) = 60 T (sum b b^T = 60 I)", c0 == 60 * sp.eye(8))
check("C1.closure hypercharge Hessian = 10 T + 12 :J_Y J_Y:",
      c2y == 12 * ((yv.T * yv)[0] * sp.eye(8) + 2 * yv * yv.T),
      "A = 12(|y|^2 I + 2 y y^T), |y|^2=5/6")
check("C1.selection hierarchy: 64 is 2-design only, 48 is not even 2-design",
      is_isotropic(R64, 2) and not is_isotropic(R64, 4) and not is_isotropic(R48, 2))

# ------------------------------------- C2: one generating function, Taylor terms
B48 = [sp.Matrix(r) / 2 for r in R48]
A0 = sum((b * b.T for b in B48), sp.zeros(8))
Rm = A0 - 12 * sp.eye(8)
check("C2.order0: A_48(0) = 12 I + R, <R,R>=48 (project values)",
      Rm == branch["R_matrix"] and A0 == branch["A48"])
A2y = sum(((yv.T * b)[0] ** 2 * b * b.T for b in B48), sp.zeros(8))
ARY = A2y - sp.Rational(5, 2) * sp.eye(8)
gram = lambda X, Y: (X.T * Y).trace() / 2  # vacuum pairing of 1/2:h^T X h: states, <T,T>=4
check("C2.order2 along y: A = 5/2 I + A_RY, tr A_RY = 0", ARY.trace() == 0)
check("C2.Gram(R,R_Y) = [[48,10],[10,35/3]] reproduced from the Taylor coefficients",
      [gram(Rm, Rm), gram(Rm, ARY), gram(ARY, ARY)] == [48, 10, sp.Rational(35, 3)])
c1y = sum(((yv.T * b)[0] * b * b.T for b in B48), sp.zeros(8))
# Mixed anomaly content of the first-order (odd) term on SM x family Cartan directions.
colour = [sp.Matrix([1, -1, 0, 0, 0, 0, 0, 0]), sp.Matrix([0, 1, -1, 0, 0, 0, 0, 0])]
weak = [sp.Matrix([0, 0, 0, 1, -1, 0, 0, 0])]
fam_dir = sp.Matrix([0] * 5 + list(branch["family"][0])) / 2
check("C2.order1 along y: Tr Y H H' = 0 for H,H' in colour, weak, Y (anomaly-free)",
      all((H.T * c1y * K)[0] == 0 for H in colour + weak + [yv] for K in colour + weak + [yv]),
      "first-order seam-holonomy response = mixed anomaly tensor")
check("C2.order1 along y is NOT the zero matrix (chiral selection)", c1y != sp.zeros(8),
      "nonzero only off the gauge Cartan block")

# ------------------------------------------ C3: commutant and hypercharge line
pairs = _native_pairs()
alpha = [sp.Matrix(pairs["algebra"].roots[i]) / 2 for i, _ in pairs["currents"]]
fam_roots = [alpha[a] - alpha[b] for a, b in combinations((1, 2, 3), 2)]
r240set = {tuple(int(x) for x in r) for r in R240}  # already doubled
check("C3.family transitions alpha_a-alpha_b are original E8 roots",
      all(tuple(int(2 * x) for x in g) in r240set for g in fam_roots))
constraints = sp.Matrix.vstack(*[g.T for g in colour + weak + fam_roots[:2]])
comm = constraints.nullspace()
X = sp.Matrix([1] * 5 + [0] * 3)
U = sp.Matrix(list(sum(alpha[1:4], sp.zeros(8, 1))))  # family-number direction
check("C3.Cartan commutant of su3_c+su2_w+su3_fam has dim 3", len(comm) == 3)
basis = sp.Matrix.hstack(*comm)
line = (sp.Matrix.vstack(X.T, U.T) * basis).nullspace()
yline = basis * line[0] if len(line) == 1 else None
check("C3.commutant orthogonal to carrier number X and family number U is ONE line = y",
      yline is not None and sp.Matrix.hstack(yline, yv).rank() == 1)

# ---------------------- C4: what each layer of W can do on the six neutrino pairs
def pair_values(A: sp.Matrix) -> list:
    """Q_A zero mode on the six symmetric retained pair states (exact lattice formula).

    a=b: e^{2 alpha_a};  a<b: h_{-1}^{v} e^{alpha_a+alpha_b}, v=alpha_a-alpha_b.
    Returns diagonal values and whether each state is an eigenvector.
    """
    out = []
    for a, b in combinations_with_replacement((1, 2, 3), 2):
        p = alpha[a] + alpha[b]
        val = (p.T * A * p)[0] / 2
        eig = True
        if a != b:
            v = alpha[a] - alpha[b]
            Av = A * v
            osc = (v.T * Av)[0] / (v.T * v)[0]
            eig = sp.Matrix.hstack(Av, v).rank() <= 1
            val += osc
        out.append((sp.nsimplify(val), eig))
    return out


vals_R = pair_values(Rm)
check("C4.neutral layer: R_0 = 4 on all six pairs (project value reproduced)",
      all(v == 4 and e for v, e in vals_R), str([str(v) for v, _ in vals_R]))
rY = pair_values(ARY)
check("C4.neutral layer: R_Y,0 = -5/3 on all six pairs (project value reproduced)",
      all(v == sp.Rational(-5, 3) and e for v, e in rY), str([str(v) for v, _ in rY]))
plane = sp.Matrix.hstack(*fam_roots[:2])
P2 = plane * (plane.T * plane).inv() * plane.T
check("C4.Schur reason: R and A_RY restricted to the family A2 plane are scalar",
      all(sp.simplify(P2 * M * P2 - ((P2 * M * P2).trace() / 2) * P2) == sp.zeros(8)
          and sp.simplify(P2 * M - M * P2) == sp.zeros(8) for M in (Rm, ARY)),
      "=> commutes with SU3_fam zero modes => scalar on the irreducible 6")
# principal (Kac) twist of the family A2: f in the A2 plane with f.gamma_i = 1/3
g1, g2 = alpha[1] - alpha[2], alpha[2] - alpha[3]
cf = sp.symbols("cf0:2")
fv = cf[0] * g1 + cf[1] * g2
sol = sp.solve([(fv.T * g1)[0] - sp.Rational(1, 3), (fv.T * g2)[0] - sp.Rational(1, 3)], cf)
f = fv.subs(sol)
fam_charges = [(f.T * alpha[a])[0] for a in (1, 2, 3)]
diffs = sorted(sp.nsimplify(x - min(fam_charges)) for x in fam_charges)
check("C4.principal family twist rho/3 gives family charges with residue pattern (0,1/3,2/3)",
      diffs == [0, sp.Rational(1, 3), sp.Rational(2, 3)], str(diffs))
regular = [(m1, m2) for m1 in range(3) for m2 in range(3)
           if m1 % 3 and m2 % 3 and (m1 + m2) % 3]
check("C4.Kac: regular order-3 twists of A2 are only +-rho/3 (Weyl-conjugate)",
      regular == [(1, 1), (2, 2)], str(regular))
Af = sum(((f.T * b)[0] ** 2 * b * b.T for b in B48), sp.zeros(8))
vf = pair_values(Af)
check("C4.Cartan-twist layer: order2 along f is non-scalar on the six pairs (hierarchy)",
      len({v for v, _ in vf}) > 1, str([(str(v), e) for v, e in vf]))
momenta = [tuple(alpha[a] + alpha[b]) for a, b in combinations_with_replacement((1, 2, 3), 2)]
check("C4.Cartan layer cannot mix: six pair momenta are pairwise distinct",
      len(set(momenta)) == 6, "every Q_A (any A) preserves momentum -> diagonal, no angles/phases")
moves = sum(tuple(alpha[a] + alpha[b] + s * g) in set(momenta)
            for g in fam_roots for s in (1, -1)
            for a, b in combinations_with_replacement((1, 2, 3), 2))
check("C4.root layer: E_{+-(alpha_a-alpha_b)} connect pair momenta (mixing needs directed roots)",
      moves > 0, f"{moves} momentum-preserving pair-to-pair moves")
c1f = sum(((f.T * b)[0] * b * b.T for b in B48), sp.zeros(8))
check("C4.order1 along f nonzero: response distinguishes clock orientation f vs -f",
      c1f != sp.zeros(8))

# ------------- C5: the 3+1 family selection from "port norm = topological count"
# Class: all nonempty unions of the 16 D5 x A3 blocks of the original 240 roots.
# Condition: <R_S,R_S> = |S|, i.e. the norm of the stress-orthogonal state equals
# the component count used in the original delta_top = 48 c3^4.
roots_int = [tuple(map(int, r)) for r in src["chevalley"].roots]
blk: dict = {}
for r in roots_int:
    d5, a3 = r[:5], r[5:]
    if all(abs(x) == 1 for x in r):
        key = ("S+" if sum(x == -1 for x in d5) % 2 else "S-", a3)
    elif not any(a3):
        key = ("45",)
    elif not any(d5):
        key = ("15",)
    else:
        key = ("10", a3)
    blk.setdefault(key, []).append(np.array(r))
keys = sorted(blk, key=str)
Mb = {k: sum(np.outer(x, x) for x in blk[k]) for k in keys}  # doubled coords, std = /4
setb = {k: frozenset(map(tuple, blk[k])) for k in keys}
negb = {k: next(j for j in keys if setb[j] == frozenset(tuple(-x) for x in blk[k])) for k in keys}
n_all = n_nopm = 0
carrier_solutions = []
for mask in range(1, 1 << len(keys)):
    sel = [keys[i] for i in range(len(keys)) if mask >> i & 1]
    A = sum(Mb[k] for k in sel)
    cnt = sum(len(blk[k]) for k in sel)
    Rs = 4 * A - np.trace(A) / 2 * np.eye(8)        # = 16 R_std, so <R,R> = sum(Rs^2)/512
    if int((Rs * Rs).sum()) == 512 * cnt:
        n_all += 1
        if all(negb[k] not in sel for k in sel):
            n_nopm += 1
        if all(k[0] == "S+" for k in sel):
            carrier_solutions.append(tuple(k[1] for k in sel))
check("C5.class has 16 blocks / 65535 nonempty unions", len(keys) == 16)
check("C5.norm=count: 65 unions, 32 without +-beta pairs (condition alone not enough)",
      (n_all, n_nopm) == (65, 32), f"{n_all}, {n_nopm}")
check("C5.with the P2 carrier chirality (S+ blocks only): exactly the four 3-of-4 family selections",
      len(carrier_solutions) == 4 and all(len(c) == 3 for c in carrier_solutions),
      str(len(carrier_solutions)))
kvals = {}
for k in range(1, 5):
    for keep in combinations(range(4), k):
        A = sum(Mb[("S+", tuple(branch["family"][a]))] for a in keep)
        Rs = 4 * A - np.trace(A) / 2 * np.eye(8)
        kvals.setdefault(k, set()).add(int((Rs * Rs).sum()) // 512)
check("C5.<R_S,R_S> = 16 k (4-k) for k kept family blocks; equals |S|=16k only at k=3",
      kvals == {k: {16 * k * (4 - k)} for k in range(1, 5)}, str(kvals))
# the four solutions are conjugate under an E8 Weyl element fixing the D5 coordinates
flip = lambda r, i, j: tuple(-x if t in (i, j) else x for t, x in enumerate(r))
fam_t = [tuple(f) for f in branch["family"]]
orbit = {fam_t[0]} | {tuple(flip((0,) * 5 + fam_t[0], 5 + i, 5 + j)[5:]) for i, j in combinations(range(3), 2)}
check("C5.even family sign flips (in W(D8), preserve the 240) map the omitted weight transitively",
      orbit == set(fam_t) and all(frozenset(flip(r, 5 + i, 5 + j) for r in roots_int) == frozenset(roots_int)
                                   for i, j in combinations(range(3), 2)))

# -------- C6: the charge Hessian = (mean Y^2 over a 16) x order-0 term + carrier part
B48op = A0  # A-matrix of B_48 = 12T + R
SY = A2y - (yv.T * yv)[0] / 4 * B48op
d5 = sp.diag(*([1] * 5 + [0] * 3))
check("C6.B_Y = (|y|^2/4) B_48 + S_Y with |y|^2/4 = 5/24 = mean Y^2 on a 16",
      sp.Rational(5, 24) == (yv.T * yv)[0] / 4 ==
      sum(((yv.T * sp.Matrix(r) / 2)[0]) ** 2 for r in blk[("S+", fam_t[1])]) / 16)
check("C6.S_Y lives only in the carrier D5 block and is orthogonal to T and R",
      d5 * SY * d5 == SY and SY.trace() == 0 and gram(SY, Rm) == 0 and gram(SY, SY) == sp.Rational(115, 12))

# ------- C7: the selection as a joint eigenspace of finite-order Cartan holonomies
# z4 = -(1/2)(1,1,1,1,1;0,0,0): order-4 inner element (centre of Spin(10)), its
# eigenvalue-i root space is the chiral (16,4) = the S+ blocks (the mu4 glue sector).
std = [tuple(F(x, 2) for x in r) for r in roots_int]
dot = lambda u, v: sum((a * b for a, b in zip(u, v)), F(0))
frac = lambda q: q - (q.numerator // q.denominator)
z4 = tuple([F(-1, 2)] * 5 + [F(0)] * 3)
grade4 = {}
for r, rs in zip(roots_int, std):
    grade4.setdefault(frac(dot(z4, rs)), set()).add(r)
splus = {r for k in keys if k[0] == "S+" for r in map(tuple, blk[k])}
check("C7.z4 grades the 240 as (45+15) | (16,4) | (10,6) | (16b,4b) with grade 1/4 = S+ chiral sector",
      {g: len(v) for g, v in grade4.items()} == {F(0): 52, F(1, 4): 64, F(1, 2): 60, F(3, 4): 64}
      and grade4[F(1, 4)] == splus)
# all family holonomies z = (0^5; t), t in (1/n)Z^3, n = 2..6: which eigenspaces S inside
# the grade-1/4 sector satisfy <R_S,R_S> = |S| ?
hits = set()
for nord in range(2, 7):
    for t in product(range(2 * nord), repeat=3):
        zt = tuple([F(0)] * 5 + [F(x, nord) for x in t])
        classes = {}
        for r in grade4[F(1, 4)]:
            classes.setdefault(frac(dot(zt, tuple(F(x, 2) for x in r))), []).append(r)
        for cls in classes.values():
            A = sum(np.outer(x, x) for x in map(np.array, cls))
            Rs = 4 * A - np.trace(A) / 2 * np.eye(8)
            if int((Rs * Rs).sum()) == 512 * len(cls):
                hits.add(frozenset(cls))
retained = {tuple(map(int, r)) for r in branch["retained"]}
check("C7.over all family holonomies of order 2..6 the norm=count eigenspaces are exactly the four 3-of-4 sets",
      len(hits) == 4 and all(len(h) == 48 for h in hits) and frozenset(retained) in hits, str(len(hits)))
z3 = tuple([F(0)] * 5 + [F(2, 3)] * 3)
joint = {r for r, rs in zip(roots_int, std)
         if frac(dot(z4, rs)) == F(1, 4) and frac(dot(z3, rs)) == F(1, 3)}
check("C7.S=48 = joint eigenspace exp(2 pi i(z4+z3)) = i*omega (order 12); z3 has centraliser SU3xU1",
      joint == retained)

# ---- C8: directed family permutation and diagonal principal twist are one element
refl = lambda g: (lambda v: v - (v.T * g)[0] * g)  # roots have norm 2
s1, s2 = refl(g1), refl(g2)
cyc = [s1(s2(alpha[a])) for a in (1, 2, 3)]
perm = [next(b for b in (1, 2, 3) if cyc[i] == alpha[b]) for i in range(3)]
check("C8.Coxeter element s1 s2 of the family A2 permutes the three families cyclically (order 3)",
      sorted(perm) == [1, 2, 3] and perm != [1, 2, 3]
      and all(s1(s2(s1(s2(s1(s2(alpha[a])))))) == alpha[a] for a in (1, 2, 3)), str(perm))
w = sp.exp(2 * sp.pi * sp.I / 3)
P = sp.zeros(3)
for i, b in enumerate(perm):
    P[b - 1, i] = 1
Fm = sp.Matrix(3, 3, lambda j, k: w ** (j * k)) / sp.sqrt(3)
D = (Fm.H * P * Fm).applyfunc(lambda u: sp.nsimplify(sp.simplify(sp.expand_complex(u))))
twist = [sp.exp(2 * sp.pi * sp.I * (x - fam_charges[0])) for x in fam_charges]
same = lambda xs, ys: all(any(sp.simplify(sp.expand_complex(x - y)) == 0 for y in ys) for x in xs)
roots3 = [sp.Integer(1), w, w ** 2]
check("C8.(Kostant, A2) permutation = Fourier-conjugate of the principal twist exp(2 pi i f)",
      D.is_diagonal() and same(list(D.diagonal()), roots3) and same(twist, roots3) and same(roots3, twist),
      "conjugator = Z3 Fourier matrix, |entries|^2 = 1/3; a skeleton, not a PMNS prediction")

# ---- C9: the H^1 mu4 character action on families is NOT a Spin(10)-commuting holonomy
# QGEO.MARKS.02/COHOM.01: H^1 = chi_1+chi_2+chi_3 (trivial character absent). Realising
# the four family blocks with characters {1,i,-1,-i} by z = (z_D5; z_A3) requires block
# values {0,1/4,1/2,3/4} mod 1.  Proof of impossibility for ALL real z:
#   (a) Spin(10)-commuting z_D5 lies in Z^5 or (1/2+Z)^5 and adds the SAME quarter value c
#       on every S+ block; (b) the four family weights sum to zero, so the A3 parts sum to 0.
#   Hence the four block values sum to 4c = 0 mod 1, while 0+1/4+1/2+3/4 = 3/2 = 1/2 mod 1.
# An executed 1/8-grid census (262144 z, 0 solutions) agrees; it is not rerun here.
d5_opts = [tuple(F(x) for x in v) for v in product((0, 1), repeat=5)] + \
          [tuple(F(2 * x + 1, 2) for x in v) for v in product((0, 1), repeat=5)]
s_plus_d5 = {tuple(r[:5]) for r in splus}
uniform_quarter = True
for zd in d5_opts:
    vals = {frac(sum((F(x, 2) * y for x, y in zip(d, zd)), F(0))) for d in s_plus_d5}
    uniform_quarter &= len(vals) == 1 and (4 * vals.pop()).denominator == 1
fam_sum = [sum(f[i] for f in fam_t) for i in range(3)]
check("C9.mu4 family characters {1,i,-1,-i} are not realisable by any Spin(10)-commuting Cartan holonomy",
      uniform_quarter and fam_sum == [0, 0, 0] and frac(F(3, 2)) == F(1, 2),
      "obstruction class 1/2 mod 1 (a Z2); grid census 262144/0")

# ---- C10: the deck mu4 PERMUTES the family lines; its linear lift is a double cover
# Weyl element w: (x6,x7,x8) -> (x7,-x6,-x8), identity on the D5 coordinates (in W(D8) < W(E8)).
wmap = lambda r: tuple(r[:5]) + (r[6], -r[5], -r[7])
root_set = frozenset(roots_int)
cyc_fam = {f: wmap((0,) * 5 + f)[5:] for f in fam_t}
def cycle_len(start):
    n, f = 1, cyc_fam[start]
    while f != start:
        f, n = cyc_fam[f], n + 1
    return n
check("C10.w preserves the 240, fixes the carrier coordinates and 4-cycles the four family blocks",
      frozenset(map(wmap, roots_int)) == root_set and cycle_len(fam_t[0]) == 4
      and all(frozenset(map(wmap, blk[("S+", f)])) == frozenset(map(tuple, blk[("S+", cyc_fam[f])])) for f in fam_t))
# Every SU(4) lift of a 4-cycle of lines is g = D P with det g = 1, and then g^4 = -I.
d = sp.symbols("d0:4")
Pc = sp.Matrix(4, 4, lambda i, j: 1 if i == (j + 1) % 4 else 0)
g = sp.diag(*d) * Pc
g4 = sp.expand(g ** 4)
check("C10.every SU(4) lift g of the 4-cycle has g^4 = det(D) I = -I (sgn of a 4-cycle is -1)",
      g4 == sp.prod(d) * sp.eye(4) and sp.expand(g.det()) == -sp.prod(d),
      "det g = 1 forces prod d = -1, so the deck has order 8 on the family 4: mu4 lifts to mu8")
# g^4 = -I on the family 4 is the torus element z_F = (0^5;1,0,0): it is -1 exactly on the
# 128 half-integer (spinor) roots and +1 on the 112 integer roots: the D8 spinor parity (-1)^F.
zF = tuple([F(0)] * 5 + [F(1), F(0), F(0)])
minus = {r for r, rs in zip(roots_int, std) if frac(dot(zF, rs)) == F(1, 2)}
plus = {r for r, rs in zip(roots_int, std) if frac(dot(zF, rs)) == 0}
half_int = {r for r in roots_int if all(abs(x) == 1 for x in r)}
check("C10.g^4 = (-1)^F: -1 on the 128 spinor roots (16,4)+(16b,4b), +1 on the 112 integer roots (120 = 112 roots + 8 Cartan = 45+15+(10,6))",
      minus == half_int and len(minus) == 128 and len(plus) == 112 and minus | plus == set(roots_int)
      and all(frac(dot(zF, tuple(F(x, 2) for x in r))) == F(1, 2) for f in fam_t for r in map(tuple, blk[("S+", f)])))

summary = {
    "contract": "seam-holonomy-generating-20261001",
    "checks": len(RESULTS), "passed": sum(r["pass"] for r in RESULTS),
    "pair_values": {"R": [str(v) for v, _ in vals_R], "R_Y": [str(v) for v, _ in rY],
                    "A_f": [[str(v), bool(e)] for v, e in vf]},
    "family_twist_f": [str(x) for x in f], "family_charges": [str(x) for x in fam_charges],
    "results": RESULTS,
}
(HERE / "results.json").write_text(json.dumps(summary, indent=2))
print(f"{summary['passed']}/{summary['checks']} checks passed")
sys.exit(0 if summary["passed"] == summary["checks"] else 1)
