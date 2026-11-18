#!/usr/bin/env python3
"""Session check of 18 November 2026 (third pasted note).

Original data (contract text of 10 November 2026, section 3): two local
images v = (-1, sqrt 2), w = (1, sqrt 2) in the normalised basis of the
oscillator directions (2) and (1,1); preparations v, w, v+w, v+iw; directed
local return  -w* P v / 2  for a positive contraction P; the two recorded
rank-one selectors P_+- give -1/2.

Q1  optimal binary discrimination: directions and strengths
Q2  weighting that equalises the strengths (diagonal and general real)
Q3  directed return under those weightings, all three pairings
Q4  v, w as level-two parts of opposite root vertex operators; Gram formula
Q5  two-point kernels of the free field as weightings
Q6  tetrahedron of qubit states, the two invariant forms, the degeneracy
Q7  reflections that preserve the four marks generate a finite group
"""
import itertools
from fractions import Fraction as Fr
import sympy as sp

I = sp.I
s2 = sp.sqrt(2)
v = sp.Matrix([-1, s2])
w = sp.Matrix([1, s2])
a = v+w
b = v+I*w
states = {"v": v, "w": w, "v+w": a, "v+iw": b}


def ip(x, y, G):
    return sp.simplify((x.H*G*y)[0, 0])


def fid(x, y, G):
    return sp.simplify(sp.Abs(ip(x, y, G))**2/(ip(x, x, G)*ip(y, y, G)))


def helstrom(pair1, pair2, G):
    """Squared strengths of the two difference operators and their overlap.
    For unit trace projectors: Tr D^2 = 2 - 2F, strength^2 = Tr D^2 / 2."""
    (x1, y1), (x2, y2) = pair1, pair2
    st1 = 1-fid(x1, y1, G)
    st2 = 1-fid(x2, y2, G)
    cross = fid(x1, x2, G)-fid(x1, y2, G)-fid(y1, x2, G)+fid(y1, y2, G)
    return sp.simplify(st1), sp.simplify(st2), sp.simplify(cross)


E2 = sp.eye(2)
pairings = {"(v,w | v+w,v+iw)": ((v, w), (a, b)),
            "(v,v+w | w,v+iw)": ((v, a), (w, b)),
            "(v,v+iw | w,v+w)": ((v, b), (w, a))}
print("Q1. optimal discrimination in the original two-mode space")
for name, (p1, p2) in pairings.items():
    st1, st2, cross = helstrom(p1, p2, E2)
    cosang = sp.simplify(cross/(2*sp.sqrt(st1*st2)))
    print("    %-20s squared strengths %s, %s | cosine between the two axes %s" % (name, st1, st2, cosang))
st1, st2, cross = helstrom(*pairings["(v,w | v+w,v+iw)"], E2)
assert (st1, st2, cross) == (sp.Rational(8, 9), sp.Rational(1, 3), 0)

print("Q2. weighting that makes the two full operators a quarter turn apart")
g, c = sp.symbols("g c", real=True)
Gd = sp.diag(1, g)
st1, st2, cross = helstrom(*pairings["(v,w | v+w,v+iw)"], Gd)
sol = sp.solve(sp.Eq(st1, st2), g)
print("    diagonal weighting diag(1, g): axes stay orthogonal:", sp.simplify(cross) == 0,
      "| equal strengths for g =", sol, "| common squared strength", sp.simplify(st1.subs(g, sol[0])))
assert [s for s in sol if s > 0] == [sp.Rational(1, 6)]
G6 = sp.diag(1, sp.Rational(1, 6))
tri = [v, w, -a]
print("    with diag(1, 1/6): norms of v, w, -(v+w):", [ip(x, x, G6) for x in tri],
      "| cosines:", [sp.simplify(ip(x, y, G6)/sp.sqrt(ip(x, x, G6)*ip(y, y, G6))) for x, y in itertools.combinations(tri, 2)])
# general real symmetric weighting [[1, c], [c, g]]: both conditions as polynomials, the factor
# c^2 - g (singular weighting) removed, common zeros through the resultant in c
Gs = sp.Matrix([[1, c], [c, g]])
def fidp(x_, y_):
    num = sp.expand((x_.H*Gs*y_)[0, 0])
    return sp.simplify((sp.re(num)**2+sp.im(num)**2)/(sp.expand((x_.H*Gs*x_)[0, 0])*sp.expand((y_.H*Gs*y_)[0, 0])))
res = {}
for name, ((x1, y1), (x2, y2)) in pairings.items():
    e1 = sp.numer(sp.together(fidp(x1, y1)-fidp(x2, y2)))
    e2 = sp.numer(sp.together(fidp(x1, x2)-fidp(x1, y2)-fidp(y1, x2)+fidp(y1, y2)))
    e1 = sp.cancel(e1/(c**2-g))
    e2 = sp.cancel(e2/(c**2-g))
    R = sp.factor(sp.resultant(e1, e2, c))
    good = []
    for g0 in sp.real_roots(sp.Poly(R, g)) if R != 0 else []:
        for c0 in sp.solve(e1.subs(g, g0), c):
            c0 = sp.nsimplify(sp.simplify(c0))
            if c0.is_real and sp.simplify(e2.subs({g: g0, c: c0})) == 0 and sp.simplify(g0-c0**2) > 0:
                good.append((c0, sp.nsimplify(g0)))
    res[name] = sorted(set(good), key=str)
    print("    general real weighting [[1,c],[c,g]], %-20s nonsingular positive solutions (c, g): %s" % (name, res[name]))
assert res["(v,w | v+w,v+iw)"] == [(0, sp.Rational(1, 6))]

print("Q3. directed local return  -w* P v / 2")
Pp = sp.Rational(1, 3)*sp.Matrix([[1, s2], [s2, 2]])
Pm = sp.Rational(1, 3)*sp.Matrix([[1, -s2], [-s2, 2]])
d = lambda P: sp.simplify(-(w.H*P*v)[0, 0]/2)
print("    recorded selectors P+, P-:", d(Pp), d(Pm), "| identity:", d(E2), "| diag(1, 1/6):", d(G6))
assert d(Pp) == d(Pm) == sp.Rational(-1, 2) and d(G6) == sp.Rational(1, 3)
for name, good in res.items():
    for cv, gv in good:
        G0 = sp.Matrix([[1, cv], [cv, gv]])
        top = max(G0.eigenvals().keys(), key=lambda e: sp.N(e))
        dmax = sp.simplify(d(G0)/top)                   # largest contraction in this ray
        print("    %-20s weighting (c, g) = (%s, %s): return at the largest contraction %s  (%.6f)"
              % (name, cv, gv, sp.nsimplify(dmax), float(dmax)))

print("Q4. v and w as level-two parts of root vertex operators")
# one boson, [a_m, a_n] = m delta: |a_-2|^2 = 2, |a_-1^2|^2 = 2; S_2(alpha) = alpha a_-2/2 + alpha^2 a_-1^2/2
al = sp.sqrt(2)
coords = lambda x: sp.Matrix([x/2*sp.sqrt(2), x**2/2*sp.sqrt(2)])      # normalised basis
assert coords(al) == w and coords(-al) == v
print("    S_2(+alpha) = w, S_2(-alpha) = v for alpha^2 = 2 (a root): exact")
# multi-direction check of <S_n(alpha), S_n(beta)> = binom(alpha.beta + n - 1, n) by generating function
z, x = sp.symbols("z x")
for cc in (2, 1, 0, -1, -2):
    ser = sp.series((1-z)**(-cc), z, 0, 5).removeO()
    vals = [ser.coeff(z, n) for n in range(1, 4)]
    norm = [sp.binomial(n+1, n) for n in range(1, 4)]
    print("    alpha.beta = %2d: Gram at levels 1,2,3 = %s | normalised %s" % (cc, vals, [sp.Rational(a_, b_) for a_, b_ in zip(vals, norm)]))
# explicit Fock check at level two in rank 2: alpha = (sqrt2, 0), beta at product c
def s2vec(alpha):
    """Level-two state in the basis a1_-2, a2_-2, a1_-1^2, a1_-1 a2_-1, a2_-1^2 with its Gram."""
    a1, a2 = alpha
    return sp.Matrix([a1/2, a2/2, a1**2/2, a1*a2, a2**2/2])
GramF = sp.diag(2, 2, 2, 1, 2)
for cc, beta in ((1, (1/sp.sqrt(2), sp.sqrt(sp.Rational(3, 2)))), (0, (0, sp.sqrt(2))),
                 (-1, (-1/sp.sqrt(2), sp.sqrt(sp.Rational(3, 2)))), (-2, (-sp.sqrt(2), 0))):
    A, B = s2vec((sp.sqrt(2), 0)), s2vec(beta)
    val = sp.simplify((A.T*GramF*B)[0, 0])
    assert val == sp.binomial(cc+1, 2), (cc, val)
print("    explicit two-direction Fock check of the level-two values 1, 0, 0, 1 for products 1, 0, -1, -2: passed")
# E8 roots: values of the normalised level-two Gram, symmetric four-root configurations
roots = []
for i, j in itertools.combinations(range(8), 2):
    for s_i, s_j in itertools.product((2, -2), repeat=2):
        r = [0]*8; r[i], r[j] = s_i, s_j; roots.append(tuple(r))
roots += [s for s in itertools.product((1, -1), repeat=8) if s.count(-1) % 2 == 0]
dot = lambda p, q: sum(x_*y_ for x_, y_ in zip(p, q))//4
lvl2 = lambda cc: Fr(cc*(cc+1)//2, 3)
print("    E8: normalised level-two Gram values over all root pairs:", sorted({str(lvl2(dot(p, q))) for p in roots[:40] for q in roots}))
r0 = roots[0]
nb = [q for q in roots if dot(r0, q) == 1]
simplex = next(t for t in itertools.combinations(nb, 3) if all(dot(p, q) == 1 for p, q in itertools.combinations(t, 2)))
orth = [q for q in roots if dot(r0, q) == 0]
frame = next(t for t in itertools.combinations(orth, 3) if all(dot(p, q) == 0 for p, q in itertools.combinations(t, 2)))
for nm, cfg in (("four roots at 60 degrees", (r0,)+simplex), ("four orthogonal roots", (r0,)+frame)):
    Gm = sp.Matrix(4, 4, lambda i_, j_: sp.Rational(lvl2(dot(cfg[i_], cfg[j_])).numerator, lvl2(dot(cfg[i_], cfg[j_])).denominator))
    print("    %-26s level-two Gram = %s" % (nm, "I" if Gm == sp.eye(4) else ("(2/3) I + (1/3) 11^T" if Gm == sp.Rational(2, 3)*sp.eye(4)+sp.ones(4)/3 else Gm)))
# neutrality against full symmetry: four roots with equal pairwise product and zero sum need 2 + 3c = 0
print("    four roots with zero total charge and equal pairwise products need product -2/3: impossible; ")
P1 = (1,)*8
P2 = (-1, 1, 1, 1, 1, -1, 1, 1)
GA = tuple(x_-y_ for x_, y_ in zip(P1, P2))
word = [P2, GA, tuple(-x_ for x_ in P1)]
print("    recorded three-letter word (p2, gamma, -p1): pairwise products", [dot(p, q) for p, q in itertools.combinations(word, 2)],
      "| level-one normalised", [str(Fr(dot(p, q), 2)) for p, q in itertools.combinations(word, 2)],
      "| level-two", [str(lvl2(dot(p, q))) for p, q in itertools.combinations(word, 2)])

print("Q5. two-point kernels of the free field as weightings of the modes (2) and (1,1)")
rho, th = sp.symbols("rho theta", positive=True)
zz, zb = sp.symbols("zz zb")
K = (1-zz*zb)**(-2)                                   # <J(z)^dagger J(z)> for the weight-one current
dJ = sp.simplify(sp.diff(K, zz, zb).subs(zb, rho/zz))  # derivative field = mode (2)
JJ = 2*(1-rho)**(-4)                                   # weight-two quasi-primary of norm 2 = mode (1,1)
ratio = sp.simplify(dJ/JJ)
print("    reflected pair (z, 1/conj z), rho = |z|^2: weight ratio (2):(1,1) =", sp.factor(ratio))
f = -1/(4*sp.sin(th/2)**2)                             # same-circle kernel of the current
ratio_c = sp.simplify(-sp.diff(f, th, 2)/(2*f**2))
print("    two points on the mark circle at angle theta: ratio =", sp.simplify(sp.expand_trig(ratio_c)), "= 2 + cos(theta)")
assert sp.simplify(ratio_c-(2+sp.cos(th))) == 0 and sp.simplify(ratio-(1+2*rho)) == 0
eta = sp.symbols("eta")
dk = sp.simplify(-(w.T*sp.diag(1+2*eta, 1)*v)[0, 0]/2)
print("    ratio 1 + 2 eta with 0 <= eta <= 1 (never 6); directed return relative to the (1,1) weight:", sp.factor(dk))
print("    eta = 0 (unweighted, or antipodal marks): ", dk.subs(eta, 0), "| quarter separation (eta = 1/2):", dk.subs(eta, sp.Rational(1, 2)),
      "| coincidence (eta = 1):", dk.subs(eta, 1))

print("Q6. tetrahedron of qubit states")
n = [sp.Matrix(t)/sp.sqrt(3) for t in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))]
pauli = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -I], [I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
rhos = [(sp.eye(2)+sum((nv[k]*pauli[k] for k in range(3)), sp.zeros(2)))/2 for nv in n]
GramS = sp.Matrix(4, 4, lambda i_, j_: sp.simplify((rhos[i_]*rhos[j_]).trace()))
GramB = sp.Matrix(4, 4, lambda i_, j_: sp.simplify((n[i_].T*n[j_])[0, 0]))
assert GramS == sp.Rational(2, 3)*sp.eye(4)+sp.ones(4)/3
assert GramB == sp.Rational(4, 3)*sp.eye(4)-sp.ones(4)/3 and GramB.rank() == 3
assert GramS == sp.ones(4)/2+GramB/2
print("    Tr(rho_i rho_j) = (2/3) I + (1/3) 11^T (the second form of the earlier note)")
print("    = (1/2) 11^T + (1/2) Bloch Gram; Bloch Gram = (4/3) I - (1/3) 11^T has rank 3 (cross entry -1/3)")
Dm = sp.Matrix.hstack(*[n[k]-n[0] for k in (1, 2, 3)])
print("    weights -> Bloch vector is affine with determinant", sp.simplify(Dm.det()),
      "(nonzero): one isolated degeneracy at the centre, winding number of the direction map = sign = %+d" % sp.sign(Dm.det()))
quartet = [fid(x_, y_, E2) for x_, y_ in itertools.combinations([v, w, a, b], 2)]
print("    original quartet v, w, v+w, v+iw: pairwise overlaps squared", quartet, "(not equal: not isotropic)")

print("Q7. maps of the sphere that preserve the four marks")
# Moebius and anti-Moebius maps permuting {1, i, -1, -i}: z -> u z^(+-1), with or without conjugation, u^4 = 1
# written on exponents k of the marks i^k: (u, e, conj) sends k -> u + e*(+-k) mod 4
maps = [(u, e, cj) for u in range(4) for e in (1, -1) for cj in (0, 1)]
act = lambda m, k: (m[0]+m[1]*(-k if m[2] else k)) % 4
assert all(sorted(act(m, k) for k in range(4)) == [0, 1, 2, 3] for m in maps)
def compose(m1, m2):                                   # m1 after m2
    u1, e1, c1 = m1; u2, e2, c2 = m2
    return ((u1+e1*(-u2 if c1 else u2)) % 4, e1*e2, (c1+c2) % 2)
grp = set(maps)
assert all(compose(m1, m2) in grp for m1 in maps for m2 in maps)
refl = [m for m in maps if m[2] == 1 and compose(m, m) == (0, 1, 0)]
prods = {compose(m1, m2) for m1 in refl for m2 in refl}
def order(m):
    k, cur = 1, m
    while cur != (0, 1, 0):
        cur, k = compose(m, cur), k+1
    return k
print("    group order", len(grp), "| reflections", len(refl), "| orders of products of two reflections:", sorted({order(m) for m in prods}))
print("    every product has finite order: no dilation z -> lambda z with |lambda| != 1 arises")

print("Q8. recorded sphere family: ratio 2 C_u(a, -3a) / C_u(a, -a)")
import mpmath as mp
mp.mp.dps = 30
def C_same(u):                                         # R * C_u(x, -x), recorded closed form
    rest = mp.nsum(lambda l: (-1)**int(l)*((l+mp.mpf(1)/2)/mp.sqrt((l+mp.mpf(1)/2)**2+u)-1), [0, mp.inf])
    return 1/(8*mp.pi)+rest/(4*mp.pi)
def C_far(u, k=3):                                     # R * C_u(x, -k x): radial modes (R/r')^nu / (2 nu sqrt(R r'))
    tot = mp.nsum(lambda l: (-1)**int(l)*(2*l+1)*mp.mpf(k)**(-mp.sqrt((l+mp.mpf(1)/2)**2+u))/(2*mp.sqrt((l+mp.mpf(1)/2)**2+u)), [0, mp.inf])
    return tot/(4*mp.pi*mp.sqrt(k))
assert abs(C_same(0)-1/(8*mp.pi)) < mp.mpf(10)**-20 and abs(C_far(0)-1/(16*mp.pi)) < mp.mpf(10)**-20
vals = [(u, 2*C_far(u)/C_same(u)) for u in (0, mp.mpf(1)/100, mp.mpf(1)/4, 1, 4, 25)]
print("    undeformed: C(a,-a) = 1/(8 pi |a|), C(a,-3a) = 1/(16 pi |a|) (Coulomb kernel 1/(4 pi distance)): ratio 1")
print("    ratio at u = 0, 1/100, 1/4, 1, 4, 25:", [mp.nstr(r, 10) for _, r in vals])
mono = all(vals[k][1] < vals[k+1][1] for k in range(len(vals)-1)) or all(vals[k][1] > vals[k+1][1] for k in range(len(vals)-1))
print("    monotone on the sampled values:", mono, "(equality with 1 only at u = 0 on this sample; not a proof for all u)")

print("Q9. returns of the three-family matter space through the omitted family pair (root adjacency only)")
is_int = lambda r: all(x_ % 2 == 0 for x_ in r)
MATTER = [r for r in roots if not is_int(r)]
RS = set(roots)
addr = lambda p, q: tuple(x_+y_ for x_, y_ in zip(p, q))
fam = lambda r: tuple(r[5:]) if list(r[:5]).count(-1) % 2 == 0 else tuple(-x_ for x_ in r[5:])
labels = sorted(set(map(fam, MATTER)))
omit = labels[0]
P96 = [r for r in MATTER if fam(r) != omit]
O32 = set(r for r in MATTER if fam(r) == omit)
assert len(P96) == 96 and len(O32) == 32
classes = {
    "weak carrier x family (24)": [r for r in roots if is_int(r) and sum(1 for i_ in (3, 4) if r[i_]) == 1 and sum(1 for i_ in (5, 6, 7) if r[i_]) == 1],
    "colour carrier x family (36)": [r for r in roots if is_int(r) and sum(1 for i_ in (0, 1, 2) if r[i_]) == 1 and sum(1 for i_ in (5, 6, 7) if r[i_]) == 1],
    "family x family (12)": [r for r in roots if is_int(r) and sum(1 for i_ in (5, 6, 7) if r[i_]) == 2],
    "weak pair (4)": [r for r in roots if is_int(r) and r[3] and r[4]],
}
for nm, cls in classes.items():
    counts = sorted({sum(1 for gmm in cls if addr(a_, gmm) in O32) for a_ in P96})
    print("    %-30s number of returns through the omitted pair per matter root: %s%s"
          % (nm, counts, "  -> scalar" if len(counts) == 1 else ""))
print("ALL CHECKS PASSED")
