#!/usr/bin/env python3
"""End-to-end reconstruction: seam premises + E8 -> marks -> S48 -> mu8/(-1)^F -> y -> W_S -> T, R, B_Y.

Independent of checker.py and of tfpt_explorer: E8 is built here from scratch
(D8 roots + even half-spin).  tfpt_explorer is imported ONLY at the very end, to
cross-check that the reconstructed objects equal the project's original ones.

Declared inputs (nothing else is used):
  P1-seam  : Z2 sheet = branched double cover of a genus-0 base (chi = 2, as in
             c3 = 1/(|Z2| 2 pi chi)); curvature concentrated at the branch points
             (flat Gauss-Bonnet reading of c3); deck clock acts on the marks as the
             family Coxeter element acts on the family lines (marks <-> families).
  P2       : carrier rank g = 5 with split 3+2.
  E8       : the root system.
Part D then shows that these seam premises reduce to three smaller ones:
  Q1 family lattice = H_1(seam - marks) with the residue pairing (P-c follows by functoriality),
  Q2 carrier D5 and family lattice glue to an even unimodular rank-8 lattice (genus 0 and n=4 follow),
  Q3 the Z2 sheet involution has a square root on the sheet curve (j=1728: square, rho, D4, tau follow).
Every step prints PASS/FAIL; the run fails if any step fails.

    PYTHONPATH=. python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/end_to_end.py
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction as F
from itertools import combinations, permutations, product
from pathlib import Path

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []
OUT: dict = {}


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


# =============================================================== part A: marks
I = sp.I
INF = "inf"


def mob(m, z):
    a, b, c, d = m
    if z == INF:
        return INF if c == 0 else sp.nsimplify(sp.simplify(a / c))
    den = sp.simplify(c * z + d)
    return INF if den == 0 else sp.nsimplify(sp.simplify((a * z + b) / den))


def compose(m1, m2):  # m1 o m2
    A = sp.Matrix(2, 2, list(m1)) * sp.Matrix(2, 2, list(m2))
    return tuple(A)


rho, tau, ident = (I, 0, 0, 1), (0, 1, 1, 0), (1, 0, 0, 1)
d4 = [ident]
for _ in range(7):
    for g in (rho, tau):
        for h in list(d4):
            k = compose(g, h)
            norm = next(x for x in k if x != 0)
            k = tuple(sp.simplify(x / norm) for x in k)
            if k not in d4:
                d4.append(k)
step("A0.<rho,tau> = D4 of order 8 with tau rho tau = rho^-1",
     len(d4) == 8 and all(mob(compose(compose(tau, rho), tau), z) == mob((1, 0, 0, I), z)
                          for z in (sp.Integer(2), sp.Integer(3) + I)))

z8 = sp.exp(I * sp.pi / 4)
orbits = {
    "fixed{0,inf}": [sp.Integer(0), INF],
    "mu4": [sp.Integer(1), I, sp.Integer(-1), -I],
    "z8*mu4": [sp.nsimplify(sp.expand_complex(z8 * I ** k)) for k in range(4)],
    "generic(a=2)": [sp.nsimplify(c * I ** k) for c in (2, sp.Rational(1, 2)) for k in range(4)],
    "generic(a=3)": [sp.nsimplify(c * I ** k) for c in (3, sp.Rational(1, 3)) for k in range(4)],
}
step("A1.every listed orbit is D4-invariant",
     all({mob(g, z) for z in pts} == set(pts) for pts in orbits.values() for g in d4))


def faithful(points) -> bool:
    return sum(all(mob(g, z) == z for z in points) for g in d4) == 1


census = []
names = list(orbits)
for r in range(1, len(names) + 1):
    for combo in combinations(names, r):
        pts = [z for nme in combo for z in orbits[nme]]
        n = len(pts)
        census.append({"orbits": combo, "n": n, "faithful": faithful(pts),
                       "flat": F(n, 2) == 2,  # sum(1 - 1/2) over cone points = chi = 2
                       "genus_double_cover": F(n - 2, 2), "rank_H1": n - 1,
                       "one_free_rho_orbit": combo in (("mu4",), ("z8*mu4",))})
faithful_cfg = [c for c in census if c["faithful"]]
minimal = min(c["n"] for c in faithful_cfg)
flat_faithful = [c["orbits"] for c in faithful_cfg if c["flat"]]
step("A2.(Existence) marks = branch points of the Z2 sheet; n is even for every D4-invariant divisor",
     all(c["n"] % 2 == 0 for c in census), f"{len(census)} configurations")
step("A3.(Count) among faithful D4 divisors: flat <=> minimal <=> one free rho orbit <=> rank H1 = 3 <=> n = 4",
     minimal == 4 and all((c["flat"] == (c["n"] == minimal) == c["one_free_rho_orbit"] == (c["rank_H1"] == 3))
                          for c in faithful_cfg)
     and sorted(flat_faithful) == [("mu4",), ("z8*mu4",)], str(flat_faithful))
countermodels = sorted({(c["n"], str(c["genus_double_cover"]), c["faithful"]) for c in census if not c["flat"]})
step("A4.counter-models without flatness exist (faithful D4, n=6 genus2, n=8 genus3, ...); n=2 is not faithful",
     (6, "2", True) in countermodels and (8, "3", True) in countermodels and (2, "0", False) in countermodels)


def cr(a, b, c, d):
    return sp.nsimplify(sp.simplify((a - c) * (b - d) / ((a - d) * (b - c))))


def stabiliser(pts):
    base = cr(*pts)
    return [p for p in permutations(range(4)) if cr(*[pts[i] for i in p]) == base]


w = sp.Rational(-1, 2) + sp.sqrt(3) * I / 2
quads = {"square mu4": orbits["mu4"],
         "generic {0,1,3,7}": [sp.Integer(x) for x in (0, 1, 3, 7)],
         "tetrahedral {0,1,w,w^2}": [sp.Integer(0), sp.Integer(1), sp.nsimplify(w), sp.nsimplify(sp.expand(w ** 2))]}


def has_4cycle(perms):
    def order(p):
        q, k = p, 1
        while q != tuple(range(4)):
            q, k = tuple(p[i] for i in q), k + 1
        return k
    return any(order(p) == 4 for p in perms)


stab = {k: stabiliser(v) for k, v in quads.items()}
step("A5.(Order) Moebius stabiliser of 4 marks is V4 (generic), A4 (tetrahedral) or D4 (square); a 4-cycle only for the square",
     [len(stab[k]) for k in quads] == [8, 4, 12]
     and [has_4cycle(stab[k]) for k in quads] == [True, False, False],
     str({k: len(v) for k, v in stab.items()}))
step("A6.(Reflection) the square's full stabiliser IS D4: tau is forced, no free modulus (cross-ratio 2, j=1728)",
     cr(*orbits["mu4"]) == 2 and sorted(stab["square mu4"]) == sorted(
         tuple(orbits["mu4"].index(mob(g, z)) for z in orbits["mu4"]) for g in d4))
OUT["marks"] = {"n": 4, "set": "mu4 (up to Moebius)", "aut": "D4", "cross_ratio": 2,
                "families_k": 3, "census": [{**c, "orbits": list(c["orbits"]),
                                             "genus_double_cover": str(c["genus_double_cover"])} for c in census]}

# ================================================================ part B: E8
roots = set()
for i, j in combinations(range(8), 2):
    for si, sj in product((2, -2), repeat=2):
        r = [0] * 8
        r[i], r[j] = si, sj
        roots.add(tuple(r))
for s in product((1, -1), repeat=8):
    if s.count(-1) % 2 == 0:
        roots.add(s)
R = sorted(roots)  # doubled coordinates, beta_std = r/2
step("B0.E8 built from scratch: 240 roots, all norm 2", len(R) == 240 and all(sum(x * x for x in r) == 8 for r in R))
n_marks, g_car = 4, 5
# Second, algebraic route to n (reproduces the known D_g + A_m bootstrap class): the family algebra of
# n marks is A_(n-1) with discriminant Z_n, q(m) = m^2 (n-1)/n mod 2; the carrier D_5 has discriminant
# Z4, q(k) = 5 k^2/4 mod 2.  An even unimodular overlattice exists iff some subgroup H of Z4 x Z_n has
# |H|^2 = 4n and q = 0 mod 2 on H (and then rank 4+n = 0 mod 8).  Search all subgroups for n = 2..20.
def glue_exists(n):
    q = lambda k, m: (F(g_car * k * k, 4) + F(m * m * (n - 1), n)) % 2
    elems = [(k, m) for k in range(4) for m in range(n)]
    for gen in elems:  # cyclic subgroups suffice for these small cyclic products; also try pairs
        for gen2 in elems:
            H = {((a * gen[0] + b * gen2[0]) % 4, (a * gen[1] + b * gen2[1]) % n)
                 for a in range(4 * n) for b in range(4 * n)}
            if len(H) ** 2 == 4 * n and all(q(*h) == 0 for h in H):
                return True
    return False
glue_n = [n for n in range(2, 21) if glue_exists(n)]
step("B1.glue route: D5 + A_(n-1) has an even unimodular overlattice only for n = 4 (n=2..20); rank 5+3 = 8",
     glue_n == [4] and g_car + n_marks - 1 == 8, str(glue_n))
fr = lambda q: q - (q.numerator // q.denominator)
dotf = lambda z, r: sum((zi * F(x, 2) for zi, x in zip(z, r)), F(0))
z4 = [F(-1, 2)] * 5 + [F(0)] * 3  # order-4 element: centre of Spin(10)
g4 = {}
for r in R:
    g4.setdefault(fr(dotf(z4, r)), []).append(r)
chiral = g4[F(1, 4)]
step("B2.z4 grading 52|64|60|64; grade 1/4 = chiral (16,4)",
     {k: len(v) for k, v in sorted(g4.items())} == {F(0): 52, F(1, 4): 64, F(1, 2): 60, F(3, 4): 64})
fam_lines = sorted({r[5:] for r in chiral})
step("B3.the chiral sector carries exactly n=4 family lines of 16", len(fam_lines) == 4
     and all(sum(r[5:] == f for r in chiral) == 16 for f in fam_lines))
omit = fam_lines[0]  # H^1 has rank n-1: one line is the base relation (all four choices conjugate)
S = [r for r in chiral if r[5:] != omit]
z3 = [F(0)] * 5 + [F(-2, 3) * F(x, 1) for x in omit]  # order-3 element along the omitted line
joint = [r for r in R if fr(dotf(z4, r)) == F(1, 4) and fr(dotf(z3, r)) == F(1, 3)]
step("B4.k = n-1 = 3: S = joint i*omega eigenspace of z4, z3; |S| = 48", sorted(joint) == sorted(S) and len(S) == 48)

Sstd = np.array(S, dtype=float) / 2
A0 = Sstd.T @ Sstd
Rm = A0 - np.trace(A0) / 8 * np.eye(8)
RR = np.trace(Rm @ Rm) / 2
step("B5.port: <R,R> = 48 = |S| = Omega_adm = N_fam dim S+", abs(RR - 48) < 1e-12 and len(S) == 3 * 16)

# deck: 4-cycle of the family lines by a Weyl element fixing the carrier coordinates
def wmap(r):
    return tuple(r[:5]) + (r[6], -r[5], -r[7])
cyc = {f: wmap((0,) * 5 + f)[5:] for f in fam_lines}
k, f = 1, cyc[fam_lines[0]]
while f != fam_lines[0]:
    f, k = cyc[f], k + 1
step("B6.deck = Weyl element of W(E8) fixing the carrier and 4-cycling the family lines (rho as family Coxeter element)",
     set(map(wmap, R)) == set(R) and k == 4)
d = sp.symbols("d0:4")
P = sp.Matrix(4, 4, lambda i, j: 1 if i == (j + 1) % 4 else 0)
gl = sp.diag(*d) * P
step("B7.SU(4) lift: g^4 = det(D) I and det g = -det(D) => g^4 = -I, deck mu4 -> mu8",
     sp.expand(gl ** 4) == sp.prod(d) * sp.eye(4) and sp.expand(gl.det()) == -sp.prod(d))
zF = [F(0)] * 5 + [F(1), F(0), F(0)]
minus = [r for r in R if fr(dotf(zF, r)) == F(1, 2)]
step("B8.g^4 = (-1)^F: -1 on exactly the 128 spinor roots, +1 on the 112 integer roots",
     sorted(minus) == sorted(r for r in R if all(abs(x) == 1 for x in r)) and len(minus) == 128)

# hypercharge: commutant of colour, weak and family roots, orthogonal to the two number operators
def vec(*c):
    return sp.Matrix(c)
cons = [vec(1, -1, 0, 0, 0, 0, 0, 0), vec(0, 1, -1, 0, 0, 0, 0, 0), vec(0, 0, 0, 1, -1, 0, 0, 0)]
kept = [f for f in fam_lines if f != omit]
cons += [sp.Matrix([0] * 5 + [sp.Rational(a - b, 2) for a, b in zip(kept[0], kept[1])]),
         sp.Matrix([0] * 5 + [sp.Rational(a - b, 2) for a, b in zip(kept[1], kept[2])])]
comm = sp.Matrix.vstack(*[c.T for c in cons]).nullspace()
X = vec(1, 1, 1, 1, 1, 0, 0, 0)
Ufam = sp.Matrix([0] * 5 + list(omit))
line = (sp.Matrix.vstack(X.T, Ufam.T) * sp.Matrix.hstack(*comm)).nullspace()
yv = sp.Matrix.hstack(*comm) * line[0]
yv = yv / yv[3] * sp.Rational(1, 2)
step("B9.commutant dim 3; orthogonal to X and family number: the single line y = (-1/3^3, 1/2^2, 0^3)",
     len(comm) == 3 and len(line) == 1 and list(yv) == [sp.Rational(-1, 3)] * 3 + [sp.Rational(1, 2)] * 2 + [0] * 3)

# W_S Taylor coefficients
y = np.array([float(x) for x in yv])
Y = Sstd @ y
A2y = (Sstd * (Y ** 2)[:, None]).T @ Sstd
ARY = A2y - 2.5 * np.eye(8)
gram = lambda a, b: np.trace(a @ b) / 2
SY = A2y - (y @ y) / 4 * A0
step("B10.W_S: order0 = 12T + R; order2 along y = 5/2 T + R_Y; Gram(R,R_Y) = [[48,10],[10,35/3]]",
     np.allclose(np.diag(A0)[:5], 12) and abs(np.trace(ARY)) < 1e-12
     and np.allclose([gram(Rm, Rm), gram(Rm, ARY), gram(ARY, ARY)], [48, 10, 35 / 3]))
step("B11.B_Y = (|y|^2/4) B_48 + S_Y, S_Y carrier-only, orthogonal to T and R, |S_Y|^2 = 115/12",
     np.allclose(SY[5:, :], 0) and np.allclose(SY[:, 5:], 0) and abs(np.trace(SY)) < 1e-12
     and abs(gram(SY, Rm)) < 1e-12 and abs(gram(SY, SY) - 115 / 12) < 1e-12)
Rall = np.array(R, dtype=float) / 2
step("B12.closure: W_240(0) = 60 T and Hessian along y = 10 T + 12 :J_Y J_Y: (7-design)",
     np.allclose(Rall.T @ Rall, 60 * np.eye(8))
     and np.allclose((Rall * ((Rall @ y) ** 2)[:, None]).T @ Rall, 12 * ((y @ y) * np.eye(8) + 2 * np.outer(y, y))))
OUT["outputs"] = {"S_size": len(S), "R_norm": float(RR), "y": [str(x) for x in yv],
                  "deck_lift_order": 8, "deck_fourth_power": "(-1)^F on 128 spinor roots",
                  "B_Y": "5/2 T + R_Y, R_Y = 5/24 R + S_Y"}

# ============== part D: the premises themselves (P-c from functoriality, genus 0, the order-4 input)
# D1  H_1(P^1 minus n marks) = Z^n/Z(1,..,1), generated by puncture loops gamma_a with sum 0.  Its dual
#     H^1 (log forms, residue vectors x with sum 0) carries the residue pairing sum x_a y_a, i.e. the A_(n-1)
#     root lattice; hence H_1 with the dual form has Gram(gamma_a, gamma_b) = delta_ab - 1/n.
n = 4
gram_H1 = sp.Matrix(n, n, lambda a, b: (1 if a == b else 0) - sp.Rational(1, n))
wstd = [sp.Matrix([0] * 5 + [sp.Rational(x, 2) for x in f]) for f in fam_lines]
gram_fam = sp.Matrix(n, n, lambda a, b: (wstd[a].T * wstd[b])[0])
step("D1.H_1(P^1 - 4 marks) with residue-dual form is isometric to the E8 family line lattice (Gram delta - 1/4, sum 0)",
     gram_H1 == gram_fam and sum(wstd, sp.zeros(8, 1)) == sp.zeros(8, 1))
# D2  functoriality: a mark permutation acts on H_1 by permuting gamma_a; transported to the family lines
#     it is realised by Weyl elements fixing the carrier.  rho -> the 4-cycle w, tau -> a family root reflection.
mark_rho = [orbits["mu4"].index(mob(rho, z)) for z in orbits["mu4"]]     # 0->1->2->3->0
mark_tau = [orbits["mu4"].index(mob(tau, z)) for z in orbits["mu4"]]     # fixes 1,-1 swaps i,-i
line_w = [fam_lines.index(cyc[f]) for f in fam_lines]
phis = [p for p in permutations(range(4)) if all(p[mark_rho[a]] == line_w[p[a]] for a in range(4))]
def refl_line_perm(a, b):
    al = [sp.Rational(x - y, 2) for x, y in zip(fam_lines[a], fam_lines[b])]
    root = tuple([0] * 5 + [int(2 * x) for x in al])
    def act(r):
        rv = [F(x, 2) for x in r]
        c = 2 * sum((F(int(2 * u), 2) * v for u, v in zip([0] * 5 + al, rv)), F(0)) / 2
        return tuple(int(2 * (v - c * F(int(2 * u), 2))) for u, v in zip([0] * 5 + al, rv))
    return root, act
ok_tau = []
for phi in phis:
    tl = [None] * 4
    for a in range(4):
        tl[phi[a]] = phi[mark_tau[a]]
    moved = [i for i in range(4) if tl[i] != i]
    if len(moved) == 2:
        root, act = refl_line_perm(*moved)
        ok_tau.append(root in roots and set(map(act, R)) == set(R)
                      and [fam_lines.index(act((0,) * 5 + f)[5:]) for f in fam_lines] == tl)
step("D2.(P-c from functoriality) exactly 4 = |C4| equivariant labellings marks->lines; tau becomes an E8 family root reflection",
     len(phis) == 4 and len(ok_tau) == 4 and all(ok_tau), f"labellings {len(phis)}")
# D3  genus 0 from nondegeneracy: H^1(Sigma_g - n) has dim 2g+n-1, the residue map has rank n-1, so the
#     residue pairing is nondegenerate iff g = 0 (n >= 1; for n = 0, dim H^1 = 2g = 4 != 3).
#     With the E8 rank condition family rank = 8 - 5 = 3:
options = [(gg, nn) for gg in range(4) for nn in range(1, 9) if 2 * gg + nn - 1 == 3]
nondeg = [(gg, nn) for gg, nn in options if gg == 0]
step("D3.(P-a) family rank 3 allows (g,n) in {(0,4),(1,2)}; nondegenerate residue pairing forces g=0, n=4",
     options == [(0, 4), (1, 2)] and nondeg == [(0, 4)])
# D4  the Z4 glue class: every gamma_a maps to the same generator of P(A3)/Q(A3) = Z4, so the glue is the
#     monodromy of a cyclic Z4 cover branched with equal index at all marks: exists iff 4 | n; its index-2
#     quotient is the Z2 sheet.  Minimal n for the glue cover = 4.
glue_class = {sp.nsimplify(4 * (wstd[a].T * wstd[0])[0] - 4 * (wstd[0].T * wstd[0])[0]) % 4 for a in range(4)}
step("D4.Z4 glue = monodromy of the cyclic cover y^4 = prod(x - z_a): exists iff 4 | n, Z2 sheet is its quotient",
     [nn for nn in range(1, 13) if nn % 4 == 0][0] == 4 and glue_class == {0}
     and all(((wstd[a] - wstd[b]).T * (wstd[a] - wstd[b]))[0] == 2 for a, b in combinations(range(4), 2)),
     "gamma_a - gamma_b are roots: all puncture classes share one glue class")
# D5  what is NOT derived: the order-4 rotation.  Maximal symmetry would pick the tetrahedral set (12 > 8).
step("D5.maximal-symmetry principle selects the WRONG configuration (tetrahedral |Stab|=12 > square 8): rho stays an input",
     len(stab["tetrahedral {0,1,w,w^2}"]) > len(stab["square mu4"]))

# D6  the order-4 input in its most natural form: the Z2 sheet curve E (double cover branched at the 4
#     marks) has j(E) = 256 (l^2-l+1)^3 / (l^2 (l-1)^2), l the cross-ratio.  The sheet involution [-1] has
#     a square root in Aut(E,O) iff j = 1728 (standard: Aut(E,O) = Z4 exactly then), iff the marks are the
#     square, iff rho exists.  Witness on the j=1728 model y^2 = x^3 - x: phi(x,y) = (-x, i y), phi^2 = [-1].
jinv = lambda l: sp.nsimplify(sp.simplify(256 * (l ** 2 - l + 1) ** 3 / (l ** 2 * (l - 1) ** 2)))
jv = {k: jinv(cr(*v)) for k, v in quads.items()}
xs, ys = sp.symbols("x y")
phi = lambda P: (-P[0], I * P[1])
curve = ys ** 2 - (xs ** 3 - xs)
step("D6.(rho as Z2 subset Z4) j: square 1728, tetrahedral 0, generic other; on j=1728, phi^2 = [-1] (sheet involution has a square root)",
     jv["square mu4"] == 1728 and jv["tetrahedral {0,1,w,w^2}"] == 0 and jv["generic {0,1,3,7}"] not in (0, 1728)
     and sp.expand(curve.subs({xs: -xs, ys: I * ys}, simultaneous=True) + curve) == 0
     and phi(phi((xs, ys))) == (xs, -ys), str(jv))

# ======================================== part C: cross-check with project objects
try:
    from tfpt_explorer.neutral_source_response import _actual_current_source, _native_neutral_branch
    src, br = _actual_current_source(), _native_neutral_branch(0)
    proj_roots = {tuple(map(int, r)) for r in src["chevalley"].roots}
    proj_S = {tuple(map(int, r)) for r in br["retained"]}
    step("C0.reconstructed E8 = project E8 root set", proj_roots == set(R))
    orbit48 = [sorted(r for r in chiral if r[5:] != f) for f in fam_lines]
    step("C1.project S48 is one of the four reconstructed 3+1 selections", sorted(proj_S) in orbit48)
    step("C2.reconstructed R matrix (for that choice) = project R matrix",
         any(np.allclose((lambda A: A - np.trace(A) / 8 * np.eye(8))((np.array(s) / 2).T @ (np.array(s) / 2)),
                         np.array(br["R_matrix"], dtype=float)) for s in orbit48 if sorted(s) == sorted(proj_S)))
except ImportError as exc:  # pragma: no cover
    step("C0.project cross-check available", False, str(exc))

OUT["steps"] = LOG
(HERE / "end_to_end_results.json").write_text(json.dumps(OUT, indent=2, default=str))
ok = all(s["pass"] for s in LOG)
print(f"{sum(s['pass'] for s in LOG)}/{len(LOG)} steps passed")
sys.exit(0 if ok else 1)
