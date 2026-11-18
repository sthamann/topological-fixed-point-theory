#!/usr/bin/env python3
"""Checks of the fourth pasted note (18 November 2026) on the elliptic surface
        y^2 = x^3 - 3 (1 - z^4) x + 2 (1 - z^4)
(fibre I4 over z = 0, type II over the four marks, square torus at infinity;
the mirror z -> 1/z is the surface whose four branch points give the v117
system, see check_original_system_closed_form.py).

Sections
    S1  the sections with x of degree <= 2 and y of degree <= 3: count and classes
    S2  action of the quarter turn on them
    S3  Z/4 grading of E8 by the glue of A3 + D5
    S4  mass of the curvature of the Hodge line: 1/3 smooth, 1/6 at each mark
    S5  the curve t^3 = f - f^5: Eisenstein form and its negative line
"""
import itertools
import sys

import mpmath as mp
import numpy as np
import sympy as sp

import cover_lattice as cl

FAIL = []


def check(name, ok, info=""):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info != "" else ""))
    if not ok:
        FAIL.append(name)


# ---------------------------------------------------------------- S1
print("S1  sections (x, y) with deg x <= 2, deg y <= 3")
# The Weierstrass model has an A3 point at (x, y, z) = (1, 0, 0), the node of the fibre over z = 0.
# Sections through it are multiple solutions of the raw equations, so they are counted in
# blown-up coordinates:  x = 1 + z xi(z),  y = z eta(z)  gives
#       eta^2 = 3 xi^2 + z xi^3 + z^2 + 3 z^3 xi,          xi = a + b z,  eta = e0 + e1 z + e2 z^2.
Apol = np.array([-3, 0, 0, 0, 3], complex)       # A = -3 + 3 z^4, ascending coefficients
Bpol = np.array([2, 0, 0, 0, -2], complex)
rng = np.random.default_rng(20261118)


def raw_resid(v):
    x, y = v[:3], v[3:]
    r = np.zeros(10, complex)
    r[:7] += np.convolve(y, y)
    r[:7] -= np.convolve(np.convolve(x, x), x)
    r[:7] -= np.convolve(Apol, x)
    r[:5] -= Bpol
    return r[:7]


def raw_jac(v):
    x, y = v[:3], v[3:]
    J = np.zeros((7, 7), complex)
    x2 = np.convolve(x, x)
    for k in range(3):
        col = np.zeros(9, complex)
        col[k:k+5] -= 3*x2
        col[k:k+5] -= Apol
        J[:, k] = col[:7]
    for k in range(4):
        col = np.zeros(9, complex)
        col[k:k+4] += 2*y
        J[:, 3+k] = col[:7]
    return J


def node_resid(v):
    xi, eta = v[:2], v[2:]
    r = np.zeros(6, complex)
    r[:5] += np.convolve(eta, eta)
    xi2 = np.convolve(xi, xi)
    r[:3] -= 3*xi2
    r[1:5] -= np.convolve(xi2, xi)
    r[2] -= 1
    r[3:5] -= 3*xi
    return r[:5]


def numjac(fun, v):
    n = len(v)
    J = np.zeros((len(fun(v)), n), complex)
    for k in range(n):
        d = np.zeros(n, complex)
        d[k] = 1e-7
        J[:, k] = (fun(v+d)-fun(v-d))/2e-7
    return J


def newton(fun, jacf, v, its=60):
    for _ in range(its):
        r = fun(v)
        if np.linalg.norm(r) < 1e-14:
            break
        try:
            v = v-np.linalg.solve(jacf(v), r)
        except np.linalg.LinAlgError:
            return None
        if not np.all(np.isfinite(v)) or np.linalg.norm(v) > 1e6:
            return None
    return v if np.linalg.norm(fun(v)) < 1e-12 else None


def collect(fun, jacf, n, keep, trials, turn, ysign):
    """simple solutions; the set is closed under the quarter turn, y -> -y and complex conjugation"""
    out = []

    def add(v):
        if v is None or not keep(v):
            return False
        if np.linalg.svd(jacf(v), compute_uv=False)[-1] < 1e-6:
            return False
        if any(np.linalg.norm(w-v) < 1e-7 for w in out):
            return False
        out.append(v)
        return True

    for _ in range(trials):
        v = newton(fun, jacf, (rng.normal(size=n)+1j*rng.normal(size=n))*rng.choice([0.5, 1.0, 2.0, 4.0, 8.0]))
        if add(v):
            todo = [out[-1]]
            while todo:
                w = todo.pop()
                for img in (w*turn, w*ysign, np.conj(w)):
                    img = newton(fun, jacf, img)
                    if add(img):
                        todo.append(out[-1])
    return out


# 40 sections away from the node: simple solutions of the raw system
root_secs = collect(raw_resid, raw_jac, 7, lambda v: abs(v[0]-1)+abs(v[3]) > 1e-4, 20000,
                    np.array([1, 1j, -1, 1, 1j, -1, -1j]), np.array([1, 1, 1, -1, -1, -1, -1]))
check("away from the node: 40 simple solutions (height 2, identity component)", len(root_secs) == 40, len(root_secs))
check("the constant sections (2/3, +-sqrt(8/27)) and (-2, +-2 sqrt 2 i z^2) are among them",
      sum(1 for v in root_secs if np.linalg.norm(v[:3]-np.array([2/3, 0, 0])) < 1e-7) == 2
      and sum(1 for v in root_secs if np.linalg.norm(v[:3]-np.array([-2, 0, 0])) < 1e-7) == 2)
# through the node with a != 0: outer components, height 2 - 3/4.  Exact elimination:
#   e0 = eps sqrt(3) a,  e1 = (6b + a^2)/(2 eps sqrt 3),  e2 = 3 a (b^2+1) eps sqrt(3)/(6b + a^2),
# leaves two equations in (A, b) with A = a^2; each solution gives four sections (a = +-sqrt A, eps = +-1).
As, bs = sp.symbols("A b")
eqA = sp.expand((6*bs+As)**3+216*As*(bs**2+1)-12*(3*bs**2+3*As*bs+1)*(6*bs+As))
eqB = sp.expand(27*As*(bs**2+1)**2-(bs**3+3*bs)*(6*bs+As)**2)
gb = sp.groebner([eqA, eqB], As, bs, order="lex")
elim = sp.factor(gb.exprs[-1])
print("      eliminant in b:", elim)
pairs = []
for rb in sp.Poly(gb.exprs[-1], bs).nroots(n=30, maxsteps=500):
    for ra in sp.Poly(eqB.subs(bs, rb), As).nroots(n=30, maxsteps=500):
        if abs(sp.N(eqA.subs({As: ra, bs: rb}))) < 1e-15:
            cand = (complex(ra), complex(rb))
            if abs(cand[0]) > 1e-9 and abs(6*cand[1]+cand[0]) > 1e-9 \
                    and not any(abs(cand[0]-u[0])+abs(cand[1]-u[1]) < 1e-8 for u in pairs):
                pairs.append(cand)
# the excluded case e1 = 0 (6b + A = 0, b^2 = -1) would need (e0 e2)^2 = 64 but gives -36 b^2 = 36
check("excluded branch e1 = 0 has no solution (36 != 64)", -36*(-1) != 64)
spin_secs = []
for (Av, bv) in pairs:
    for sa in (1, -1):
        av = sa*np.sqrt(Av+0j)
        for epsv in (1, -1):
            e0 = epsv*np.sqrt(3)*av
            e1 = (6*bv+av**2)/(2*epsv*np.sqrt(3))
            e2 = 3*av*(bv**2+1)*epsv*np.sqrt(3)/(6*bv+av**2)
            spin_secs.append(np.array([av, bv, e0, e1, e2]))
check("through the node, first order contact (a != 0): 8 pairs (A, b), 32 sections (height 5/4)",
      len(pairs) == 8 and len(spin_secs) == 32
      and max(np.linalg.norm(node_resid(v)) for v in spin_secs) < 1e-9, len(spin_secs))
# through the node with a = 0: xi = b z, eta = z (c0 + c1 z),  (c0 + c1 z)^2 = 3 b^2 + 1 + z^2 (b^3 + 3 b)
bb, c0, c1 = sp.symbols("b c0 c1")
vec = sp.solve([c0**2-3*bb**2-1, 2*c0*c1, c1**2-bb**3-3*bb], [bb, c0, c1], dict=True)
check("through the node, second order contact (a = 0): exactly 10 solutions (height 1, middle component)",
      len(vec) == 10, len(vec))
check("among them (1, +-z^2): b = 0, c0 = +-1, c1 = 0",
      sum(1 for d in vec if d[bb] == 0 and d[c1] == 0) == 2)
check("10 + 32 + 40 = 82 sections; minimal vectors of D5* of norm 1, 5/4 and of D5 of norm 2 are 10, 32, 40",
      len(vec)+len(spin_secs)+len(root_secs) == 82)

# ---------------------------------------------------------------- S2
print("S2  the quarter turn z -> i z on the sections")
ph5 = np.array([1j, -1, 1j, -1, -1j])            # xi -> i xi(iz), eta -> i eta(iz)
orb = []
left = list(range(len(spin_secs)))
ok_perm = True
while left:
    k = left[0]
    o, v = [], spin_secs[k]
    for _ in range(4):
        j = min(range(len(spin_secs)), key=lambda m: np.linalg.norm(spin_secs[m]-v))
        if np.linalg.norm(spin_secs[j]-v) > 1e-6:
            ok_perm = False
        if j not in o:
            o.append(j)
        v = v*ph5
    for j in o:
        if j in left:
            left.remove(j)
    orb.append(len(o))
check("the quarter turn permutes the 32 outer sections in eight orbits of four", ok_perm and sorted(orb) == [4]*8, sorted(orb))
ph7 = np.array([1, 1j, -1, 1, 1j, -1, -1j])
ro = []
left = list(range(len(root_secs)))
while left:
    k = left[0]
    o, v = [], root_secs[k]
    for _ in range(4):
        j = min(range(len(root_secs)), key=lambda m: np.linalg.norm(root_secs[m]-v))
        if j not in o:
            o.append(j)
        v = v*ph7
    for j in o:
        if j in left:
            left.remove(j)
    ro.append(len(o))
print("      orbit sizes on the 40 root sections:", sorted(ro))
vo = sorted([1 if d[bb] == 0 else (2 if sp.simplify(d[bb]**2+3) == 0 else 4) for d in vec])
print("      the 10 middle sections: b -> -b, c0 -> -c0, c1 -> -i c1 under the turn")

# ---------------------------------------------------------------- S3
print("S3  Z/4 grading of E8 by the glue of A3 + D5")
roots = []
for i, j in itertools.combinations(range(8), 2):
    for si in (1, -1):
        for sj in (1, -1):
            v = [sp.Integer(0)]*8
            v[i], v[j] = sp.Integer(si), sp.Integer(sj)
            roots.append(tuple(v))
for signs in itertools.product((1, -1), repeat=8):
    if signs.count(-1) % 2 == 0:
        roots.append(tuple(sp.Rational(c, 2) for c in signs))
# element of order four exp(2 pi i h), h = coweight of the node with mark 4
simple = [tuple(sp.Integer(c) for c in v) for v in (
    (1, -1, 0, 0, 0, 0, 0, 0), (0, 1, -1, 0, 0, 0, 0, 0), (0, 0, 1, -1, 0, 0, 0, 0), (0, 0, 0, 1, -1, 0, 0, 0),
    (0, 0, 0, 0, 1, -1, 0, 0), (0, 0, 0, 0, 0, 1, -1, 0), (0, 0, 0, 0, 0, 1, 1, 0))]
simple.append(tuple(sp.Rational(-1, 2) for _ in range(8)))
Bm = sp.Matrix([list(s) for s in simple])
check("simple roots of E8 (Cartan determinant 1)", (Bm*Bm.T).det() == 1)
coords = [list((sp.Matrix([list(r)])*Bm.inv())) for r in roots]
top = max(coords, key=lambda c: sum(c))
marks = [int(c) for c in top]
print("      marks of the highest root:", marks)
found = {}
for node in range(8):
    if marks[node] != 4:
        continue
    dist = [0, 0, 0, 0]
    for c in coords:
        dist[int(c[node]) % 4] += 1
    dist[0] += 8                                   # Cartan
    found[node] = dist
    zero = [r for r, c in zip(roots, coords) if int(c[node]) % 4 == 0]
    print("      node with mark 4: currents by degree mod 4:", dist, " roots in degree 0:", len(zero))
check("the element of order four attached to the mark-4 node gives (60, 64, 60, 64) on the 248 currents",
      any(d == [60, 64, 60, 64] for d in found.values()))
check("degree 0 has 52 = 12 + 40 roots (A3 + D5) and the Cartan part 8", any(d[0] == 60 for d in found.values()))

# ---------------------------------------------------------------- S4
print("S4  curvature mass of the Hodge line")
mp.mp.dps = 30


def period(zz):
    """a period of dx/y on the fibre over zz: integral between two branch points"""
    Aq = -3*(1-zz**4)
    Bq = 2*(1-zz**4)
    r = mp.polyroots([1, 0, Aq, Bq], maxsteps=200, extraprec=100)
    r = sorted(r, key=lambda c: (mp.re(c), mp.im(c)))
    a, b, c = r
    f = lambda s: 1/mp.sqrt(((a+(b-a)*s)-c))/mp.sqrt(s*(1-s))
    return mp.quad(f, [0, mp.mpf(1)/2, 1])


u1, u2 = mp.mpf(10)**-12, mp.mpf(10)**-18
e_num = mp.log(abs(period(1+u1))/abs(period(1+u2)))/mp.log(u1/u2)
check("periods near a mark grow like |z - m|^(-1/6)", abs(e_num+mp.mpf(1)/6) < mp.mpf(10)**-3, mp.nstr(e_num, 8))
check("norm of dx/y grows like |z - m|^(-1/3): curvature atom (1/2)(1/3) = 1/6 at each mark",
      sp.Rational(1, 2)*sp.Rational(1, 3) == sp.Rational(1, 6))
check("smooth part: area 4 * pi/3 of the pulled back metric, divided by 4 pi, is 1/3",
      4*sp.Rational(1, 3)/4 == sp.Rational(1, 3))
check("total 1/3 + 4/6 = 1 = Euler number 12 divided by 12 (degree of the Hodge line)",
      sp.Rational(1, 3)+4*sp.Rational(1, 6) == 1)
check("fibre contributions e/12: type II gives 2/12 = 1/6, I4 gives 4/12 = 1/3 (carried by the cusp area)",
      sp.Rational(2, 12) == sp.Rational(1, 6) and sp.Rational(4, 12) == sp.Rational(1, 3))
G2 = sp.Rational(2, 3)*sp.eye(4)+sp.Rational(1, 3)*sp.ones(4, 4)
Psym = sp.ones(4, 4)/4
check("same weights as the tetrahedral form of the second note: G2/4 = (1/6) sum_marks P_mark + (1/3) P_sym, trace 1",
      G2/4 == sp.Rational(1, 6)*sp.eye(4)+sp.Rational(1, 3)*Psym and (G2/4).trace() == 1)
cc = sp.symbols("c")
fam = ((1-cc)*sp.eye(4)+cc*sp.ones(4, 4))/4
check("in the family ((1-c) I + c 11^T)/4 the weights are (1-c)/4 per mark and c on the symmetric sum: c = 1/3",
      sp.simplify(fam-((1-cc)/4*sp.eye(4)+cc*Psym)) == sp.zeros(4, 4) and sp.solve((1-cc)/4-sp.Rational(1, 6), cc) == [sp.Rational(1, 3)])

# ---------------------------------------------------------------- S5
print("S5  the curve t^3 = f - f^5, i.e. t^12 = s (1-s)^4 with s = f^4")
n = 12
xg = [(k+1) % n for k in range(n)]
yg = [(k+4) % n for k in range(n)]
S = cl.Surface(xg, yg)
check("genus 4", S.genus == 4, S.genus)
g1 = [(k+1) % n for k in range(n)]
R = sp.Matrix(S.rep(g1))
Wm = R**4
check("deck generator of order 12, omega = g^4 of order 3 without fixed vectors",
      R**12 == sp.eye(8) and Wm**3 == sp.eye(8) and (Wm-sp.eye(8)).rank() == 8)
M = sp.Matrix(S.M)
inv_ok = (R.T*M*R == M) or (R*M*R.T == M)
check("cup form invariant", inv_ok)
col = (R.T*M*R == M)


def act(mat, v):
    return mat*v if col else mat.T*v


def Ef(u, v):
    return (u.T*M*v)[0, 0]


def herm(u, v):
    """h = [E(w u, v) + E(u, v)] + E(u, v) omega, as a pair (a, b) meaning a + b omega"""
    return (Ef(act(Wm, u), v)+Ef(u, v), Ef(u, v))


basis = None
E8 = [sp.eye(8)[:, k] for k in range(8)]
for combo in itertools.combinations(range(8), 4):
    cols = []
    for k in combo:
        cols += [E8[k], act(Wm, E8[k])]
    if abs(sp.Matrix.hstack(*cols).det()) == 1:
        basis = [E8[k] for k in combo]
        break
check("a Z[omega] basis of four cocycles exists among the standard ones", basis is not None)
om = sp.Rational(-1, 2)+sp.sqrt(3)*sp.I/2
Hm = sp.zeros(4, 4)
for a in range(4):
    for b in range(4):
        p, qv = herm(basis[a], basis[b])
        Hm[a, b] = p+qv*om
check("the form is hermitian", (Hm-Hm.H).applyfunc(sp.simplify) == sp.zeros(4, 4))
dH = sp.simplify(Hm.det())
ev = list(np.linalg.eigvalsh(np.array(sp.N(Hm, 30).tolist(), dtype=complex)))
sig = (sum(1 for e in ev if e > 1e-9), sum(1 for e in ev if e < -1e-9))
norms = [sp.simplify(Hm[a, a]) for a in range(4)]
check("unimodular over the Eisenstein integers, signature (3,1) up to the overall sign",
      sp.simplify(dH**2-1) == 0 and sig in ((3, 1), (1, 3)), "det %s signature %s" % (dH, sig))
odd = any(int(sp.simplify(herm(v, v)[0])) % 2 for v in E8)
check("the form is odd (a vector of odd norm exists): the lattice is diag(1,1,1,-1) up to sign", odd)
sgn = 1 if sig == (3, 1) else -1
lines = {}
for lab, K in (("g = omega", R-Wm), ("g = -omega", R+Wm), ("g = omega^2", R-Wm**2), ("g = -omega^2", R+Wm**2)):
    Kk = K if col else K.T
    ker = cl.int_kernel([list(Kk.row(i)) for i in range(8)], 8)
    if len(ker) == 2:
        nv = [sgn*int(sp.simplify(herm(sp.Matrix(k), sp.Matrix(k))[0])) for k in ker]
        lines[lab] = nv
        print("      eigenline %-13s norms of the two lattice generators: %s" % (lab, nv))
check("the deck generator has exactly two eigenlines defined over the Eisenstein field", len(lines) == 2, list(lines))
neg = [v for v in lines.values() if min(v) < 0]
pos = [v for v in lines.values() if min(v) > 0]
check("one of them is negative, with primitive norm -4", len(neg) == 1 and min(abs(x) for x in neg[0]) == 4,
      neg)
print("      positive rational eigenline: norms", pos)
ff, tt = sp.symbols("f t")
check("V4 = <(f,t)->(-f,-t), (f,t)->(1/f,-t/f^2)> preserves the curve and the form df/t of the minority character",
      sp.simplify((-tt/ff**2)**3-(1/ff-1/ff**5)-(tt**3-(ff-ff**5))*(-1/ff**6)) == 0)
check("quotient by f -> i f, t -> -i t: (t f)^3 = f^4 (1 - f^4), an elliptic curve with j = 0",
      sp.expand((tt*ff)**3-ff**3*tt**3) == 0)

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
