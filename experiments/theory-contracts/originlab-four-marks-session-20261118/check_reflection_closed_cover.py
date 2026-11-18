#!/usr/bin/env python3
"""Session check of 18 November 2026: hypothesis H-PQ.

The reflection in the circle through the four marks, z -> 1/conj(z),
exchanges the two symplectic kernel covers P and Q.  Their compositum Yhat
is the full two-torsion cover of the base curve C (group (Z/2)^4, genus 17).

H1. Lifts of the reflections of the z-sphere to Yhat: which are involutions,
    how many circles they fix, whether they divide the curve, how they sit
    relative to the quarter lifts.
H2. Integral lattice of the six square Prym tori (three of line P, three of
    line Q) inside H^1(Yhat, Z); the other nine Prym tori; action of the
    circle reflection.

Run:  python -B check_reflection_closed_cover.py
"""
import cmath
import itertools
import os
import sys
import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cover_lattice import (build_cover_s, Surface, automorphisms, antiholomorphic_maps, mul, inv,
                           perm_order, lattice, coords, smith)

pts, x, y = build_cover_s(("P", "Q"))
S = Surface(x, y)
n, N = len(pts), 2*S.genus
auts = automorphisms(x, y)
ident = list(range(n))
print("H0. compositum of P and Q over the s-sphere: degree", n, "| orders", perm_order(x), perm_order(y),
      perm_order(mul(x, y)), "| genus", S.genus, "| deck group", len(auts))

b = pts[0]
unit = lambda v: abs(abs(v)-1) < 1e-6
quarter_of = lambda v: round(cmath.phase(v)/(cmath.pi/2)) % 4


def hol_label(j):
    a = pts[j]
    r = a[1]/b[1]
    if abs(r**4-1) < 1e-6:                              # z -> i^k z
        kz, typ = quarter_of(r), "rot"
    else:                                               # z -> i^k / z
        kz, typ = quarter_of(a[1]*b[1]), "inv"
    lab = dict(typ=typ, kz=kz, kw=None, sg=None)
    if typ == "rot" and kz == 0:
        lab["kw"] = round(cmath.phase(a[2]/b[2])/(2*cmath.pi/3)) % 3
        if lab["kw"] == 0:
            lab["sg"] = tuple(round((a[k]/b[k]).real) for k in range(3, 7))
    return lab


lab = [hol_label(j) for j in range(n)]
N16 = [j for j in range(n) if lab[j]["sg"] is not None]
assert len(N16) == 16
G48 = [j for j in range(n) if lab[j]["typ"] == "rot" and lab[j]["kz"] == 0]
quarter = [j for j in range(n) if lab[j]["typ"] == "rot" and lab[j]["kz"] == 1]
assert len(G48) == 48 and len(quarter) == 48
rep = lambda a: S.C @ S.A[:, np.concatenate([np.array(inv(a)), n+np.array(inv(a)), 2*n+np.array(inv(a))])].T
Mn = np.array(S.M.tolist(), dtype=np.int64)
rho = {j: rep(auts[j]) for j in set(N16+quarter+G48)}
for r in rho.values():
    assert (r.T @ Mn @ r == Mn).all()

# ======================================================== H1. reflections
anti = antiholomorphic_maps(x, y)
assert all(a is not None for a in anti)


def anti_label(j):
    a = pts[j]                                          # value of z after the map at the base point
    r = a[1]/b[1]
    if abs(r**4-1) < 1e-6:                              # z -> i^k / conj(z)  (|z0| = 1)
        return ("circle", quarter_of(r))
    return ("conj", quarter_of(a[1]*b[1]))              # z -> i^k conj(z)


def ovals(Rp):
    fixed = [e for e in range(n) if Rp[e] == e]
    types = [len(fixed)]
    fixed += [n+e for e in range(n) if Rp[e] == S.x[e]]
    types.append(len(fixed)-types[0])
    fixed += [2*n+e for e in range(n) if S.y[Rp[e]] == e]
    types.append(len(fixed)-types[0]-types[1])
    par = list(range(S.nV))
    def find(a):
        while par[a] != a:
            par[a] = par[par[a]]; a = par[a]
        return a
    verts = set()
    for ed in fixed:
        verts |= {S.tail[ed], S.head[ed]}
        par[find(S.tail[ed])] = find(S.head[ed])
    comps = len({find(v) for v in verts})
    fs = set(fixed)
    tp = list(range(2*n))
    def tf(a):
        while tp[a] != a:
            tp[a] = tp[tp[a]]; a = tp[a]
        return a
    for ed in range(3*n):
        if ed not in fs:
            u, v = S.duals[ed]
            tp[tf(u)] = tf(v)
    return comps, len({tf(k) for k in range(2*n)}), tuple(types)


print("H1. antiholomorphic maps of Yhat (384), by the reflection of the z-sphere they cover")
table = {}
circle_invol = []
for j, Rp in enumerate(anti):
    key = anti_label(j)
    inv_ = mul(Rp, Rp) == ident
    entry = table.setdefault(key, dict(lifts=0, involutions=0, kinds={}))
    entry["lifts"] += 1
    if inv_:
        entry["involutions"] += 1
        c, pcs, types = ovals(Rp)
        entry["kinds"][(c, pcs)] = entry["kinds"].get((c, pcs), 0)+1
        if key == ("circle", 0):
            circle_invol.append((j, c, pcs, types))
for key in sorted(table):
    name = {"circle": "z -> i^%d / conj(z)", "conj": "z -> i^%d conj(z)"}[key[0]] % key[1]
    print("    %-22s lifts %2d | involutions %2d | (fixed circles, pieces): count %s"
          % (name, table[key]["lifts"], table[key]["involutions"], table[key]["kinds"]))
# black vertices are the points over the four marks
Bid = S.tail[:n]
print("    mark-circle reflection z -> 1/conj(z): involutive lifts", len(circle_invol))
for j, c, pcs, types in circle_invol[:1]:
    Rp = anti[j]
    fixedB = len({Bid[e] for e in range(n) if Bid[Rp[e]] == Bid[e] and Rp[e] in
                  [e2 for e2 in range(n) if Bid[e2] == Bid[e]]})
print("    kinds among them:", sorted({(c, pcs, types) for _, c, pcs, types in circle_invol}),
      "(types = fixed edges over |z|=1, over the mark rays, over the diagonal rays)")
# relation to the quarter lifts
comm = inverts = 0
for j, c, pcs, types in circle_invol:
    Rp = anti[j]
    for q in quarter:
        U = auts[q]
        V = mul(Rp, mul(U, Rp))
        comm += V == U
        inverts += V == inv(U)
print("    pairs (involutive circle lift, quarter lift): commuting", comm, "| inverting", inverts,
      "| of", len(circle_invol)*len(quarter))
conj_inv = [j for j, Rp in enumerate(anti) if anti_label(j) == ("conj", 0) and mul(Rp, Rp) == ident]
comm = inverts = 0
for j in conj_inv:
    Rp = anti[j]
    for q in quarter:
        V = mul(Rp, mul(auts[q], Rp))
        comm += V == auts[q]
        inverts += V == inv(auts[q])
print("    pairs (involutive lift of z -> conj(z), quarter lift): commuting", comm, "| inverting", inverts,
      "| of", len(conj_inv)*len(quarter))

# ============================================== H2. lattice of the Prym tori
print("H2. Prym tori inside H^1(Yhat, Z)")
gens = []
for j in N16:                                           # four generators of (Z/2)^4
    if sum(v == -1 for v in lab[j]["sg"]) == 1:
        gens.append(j)
assert len(gens) == 4
one = np.eye(N, dtype=np.int64)
chars = [c for c in itertools.product((0, 1), repeat=4) if any(c)]
chi = lambda c, j: int(np.prod([lab[j]["sg"][k] for k in range(4) if c[k]]))
kind = lambda c: "P" if c[2] == c[3] == 0 else ("Q" if c[0] == c[1] == 0 else "mixed")
L, forms, Jx = {}, {}, {}
toM = lambda a: sp.Matrix(a.tolist())
for c in chars:
    L[c] = lattice([toM(rho[j]-chi(c, j)*one) for j in gens], N)
    assert L[c].shape[1] == 2
    forms[c] = smith(L[c].T*S.M*L[c])
    Jx[c] = None
    for q in quarter:
        img = toM(rho[q])*L[c]
        try:
            X = coords(L[c], img)
        except AssertionError:
            continue
        if X*X == -sp.eye(2):
            Qf = (L[c].T*S.M*L[c])*X
            Jx[c] = X if Qf[0, 0] > 0 else -X
            break
for k in ("P", "Q", "mixed"):
    cs = [c for c in chars if kind(c) == k]
    print("    %-5s tori: %d | forms %s | square (quarter lift with square -1): %s"
          % (k, len(cs), sorted({tuple(forms[c]) for c in cs}), [Jx[c] is not None for c in cs]))


def hull(cs, tag):
    Bsum = sp.Matrix.hstack(*[L[c] for c in cs])
    coef = sum((sum(chi(c, j) for c in cs)*rho[j] for j in N16), np.zeros((N, N), dtype=np.int64))
    Lam = lattice([toM(16*one-coef)], N)
    T = coords(Lam, Bsum)
    res = dict(rank=Lam.shape[1], index=abs(T.det()), glue=[d for d in smith(T) if d != 1],
               form=smith(Lam.T*S.M*Lam))
    if all(Jx[c] is not None for c in cs):
        Jsum = sp.diag(*[Jx[c] for c in cs])
        JL = T*Jsum*T.inv()
        res["J_integral"] = all(v.is_integer for v in JL)
        Q = T.inv().T*((Bsum.T*S.M*Bsum)*Jsum)*T.inv()
        res["Q_integral"] = all(v.is_integer for v in Q)
        if res["Q_integral"]:
            res["Q_divisors"] = smith(Q)
            res["Q_even"] = all(Q[i, i] % 2 == 0 for i in range(Q.shape[0]))
    print("    %-22s %s" % (tag, res))
    return Lam, T


Pc = [c for c in chars if kind(c) == "P"]
Qc = [c for c in chars if kind(c) == "Q"]
LamP, _ = hull(Pc, "three P tori:")
LamQ, _ = hull(Qc, "three Q tori:")
Lam6, T6 = hull(Pc+Qc, "six square tori:")
both = sp.Matrix.hstack(LamP, LamQ)
X6 = coords(Lam6, both)
print("    six-torus lattice over (P lattice + Q lattice): index", abs(X6.det()),
      "| quotient", [d for d in smith(X6) if d != 1])

# action of an involutive lift of the mark-circle reflection
j0 = circle_invol[0][0]
Rm = S.rep_real(anti[j0])
assert Rm*Rm == sp.eye(N) and Rm.T*S.M*Rm == -S.M
swap = {}
for c in Pc+Qc:
    img = Rm*L[c]
    for c2 in Pc+Qc:
        try:
            coords(L[c2], img); swap[c] = c2
        except AssertionError:
            pass
print("    circle reflection on the six tori:", {kind(c)+str(c): kind(d)+str(d) for c, d in swap.items()})
fix = lattice([Rm-sp.eye(N)], N)
anti_ = lattice([Rm+sp.eye(N)], N)
tot = sp.Matrix.hstack(fix, anti_)
print("    on H^1(Yhat, Z): fixed rank", fix.shape[1], "| reversed rank", anti_.shape[1],
      "| index of their sum: 2^%d" % sp.log(abs(tot.det()), 2))
R6 = coords(Lam6, Rm*Lam6)
f6 = lattice([R6-sp.eye(12)], 12)
a6 = lattice([R6+sp.eye(12)], 12)
print("    on the six-torus lattice: fixed rank", f6.shape[1], "| reversed rank", a6.shape[1],
      "| index of their sum: 2^%d" % sp.log(abs(sp.Matrix.hstack(f6, a6).det()), 2))

# ---------------------------------------- H3. the glue code of the six tori
print("H3. glue code of the six square tori in orthonormal torus coordinates")
def adapted(c):
    """Basis (b, Jb) of L[c] with E(b, Jb) minimal (= 4)."""
    G = (L[c].T*S.M*L[c])*Jx[c]
    best = None
    for u in itertools.product(range(-6, 7), repeat=2):
        if u != (0, 0):
            v = sp.Matrix(u)
            q = (v.T*G*v)[0, 0]
            if best is None or q < best[0]:
                best = (q, v)
    assert best[0] == 4
    return sp.Matrix.hstack(L[c]*best[1], L[c]*Jx[c]*best[1])
order6 = Pc+Qc
Bad = sp.Matrix.hstack(*[adapted(c) for c in order6])
Qad = sp.Matrix(12, 12, lambda i, j: 0)
Tad = coords(Lam6, Bad)
half = 2*Tad.inv()
assert all(v.is_integer for v in half)
words = set()
cols = [tuple(int(v) % 2 for v in half[:, k]) for k in range(12)]
span = {tuple([0]*12)}
for cw in cols:
    span |= {tuple((a+b) % 2 for a, b in zip(w, cw)) for w in span}
assert len(span) == 16
wt = {}
for w in span:
    wt[sum(w)] = wt.get(sum(w), 0)+1
print("    code length 12, dimension 4, weight distribution", dict(sorted(wt.items())))
supp = [k for k in range(12) if any(w[k] for w in span)]
print("    coordinates used:", len(supp), "| pairs (torus) fully used:", sum(1 for t in range(6) if 2*t in supp and 2*t+1 in supp))
pairsum = {tuple((w[2*t]+w[2*t+1]) % 2 for t in range(6)) for w in span}
print("    J-swap invariant:", all(tuple(w[k ^ 1] for k in range(12)) in span for w in span),
      "| words with equal entries on every torus:", sum(1 for w in span if all(w[2*t] == w[2*t+1] for t in range(6))))
for w in sorted(span, key=lambda w: (sum(w), w))[:9]:
    print("     ", "".join(map(str, w[:6])), "|", "".join(map(str, w[6:])))
# the lattice with the form E(x, Jx)/2: minimal vectors
Qfull = Tad.inv().T*((Bad.T*S.M*Bad)*sp.diag(*[sp.Matrix([[0, -1], [1, 0]])]*6))*Tad.inv()
if Qfull[0, 0] < 0:
    Qfull = -Qfull
assert Qfull == Qfull.T and all(v % 2 == 0 for v in Qfull)
Qh = Qfull/2
# count vectors of norm 2 in {x/sqrt2 : x in Z^12, x mod 2 in code}: weight-4 words give 2^4 each, 2e_k give 2 each
roots = 2*12+sum(16 for w in span if sum(w) == 4)
print("    lattice (E(x,Jx)/2): determinant", Qh.det(), "| even:", all(Qh[i, i] % 2 == 0 for i in range(12)),
      "| vectors of norm 2:", roots)

# ------------------------------------------------ H4. the dividing reflections
print("H4. dividing reflections of Yhat")
div = []
for j, Rp in enumerate(anti):
    if mul(Rp, Rp) == ident:
        c, pcs, types = ovals(Rp)
        if pcs == 2:
            div.append((j, c, types))
print("    count", len(div), "| over", sorted({anti_label(j) for j, _, _ in div}), "| circles", sorted({c for _, c, _ in div}),
      "| fixed edges (mark circle, mark rays, diagonal rays)", sorted({t for _, _, t in div}))
print("    half surface: genus", (S.genus+1-div[0][1])//2, "with", div[0][1], "boundary circles")
jd = div[0][0]
Rd = S.rep_real(anti[jd])
assert Rd*Rd == sp.eye(N) and Rd.T*S.M*Rd == -S.M
fd, ad = lattice([Rd-sp.eye(N)], N), lattice([Rd+sp.eye(N)], N)
lam = sp.log(abs(sp.Matrix.hstack(fd, ad).det()), 2)
assert lam == S.genus+1-div[0][1]                      # Comessatti characteristic g + 1 - (number of circles)
assert sp.log(abs(tot.det()), 2) == S.genus+1-circle_invol[0][1]
print("    control: index of fixed + reversed lattice = 2^(g+1-circles) for the dividing and for the mark-circle lift")
kinds = {}
for c in Pc+Qc:
    img = Rd*L[c]
    try:
        Xc = coords(L[c], img)
    except AssertionError:
        kinds[kind(c)+str(c)] = "moved"
        continue
    assert Xc*Xc == sp.eye(2) and Xc*Jx[c] == -Jx[c]*Xc
    f2, a2 = lattice([Xc-sp.eye(2)], 2), lattice([Xc+sp.eye(2)], 2)
    kinds[kind(c)+str(c)] = "index %d" % abs(sp.Matrix.hstack(f2, a2).det())
print("    action on the six square tori (real structure; index 1 = two real circles, 2 = one):", kinds)
comm = sum(mul(anti[jd], mul(auts[q], anti[jd])) == auts[q] for q in quarter)
invq = sum(mul(anti[jd], mul(auts[q], anti[jd])) == inv(auts[q]) for q in quarter)
print("    quarter lifts: commuting", comm, "| inverted", invq, "| of", len(quarter))
# does a dividing reflection exist on Y_q alone?  (answered in check_gluing_and_reflection.py: none of 20)
print("ALL CHECKS PASSED")
