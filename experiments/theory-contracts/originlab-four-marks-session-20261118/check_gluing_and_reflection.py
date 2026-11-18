#!/usr/bin/env python3
"""Session check of 18 November 2026.

G. Gluing test for H-K': the four square tori (three Prym curves of the sign
   lines of the kernel cover Y_q, and the mark torus y^2 = z^4 - 1) inside
   the integral cohomology of Y_q (genus 5) and of the common cover X17
   (genus 17).  Compared with the recorded Gaussian E8 (four coordinate tori
   with twice the unit form, glue = Hamming code, index 16).
R. Second pasted note: the two S4-invariant positive forms, deck map against
   boundary reflection, the real structure R with R U R = U^-1, and the
   flavour transition matrix of v117 on the holomorphic three-block.

Run:  python -B check_gluing_and_reflection.py
"""
import cmath
import itertools
import sys
import os
import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cover_lattice import (build_cover, Surface, automorphisms, mul, inv, perm_order,
                           lattice, coords, smith, match, radicands)

I2 = sp.eye(2)


def labels(pts, use_v):
    out = []
    b = pts[0]
    for a in pts:
        kz = round(cmath.phase(a[0]/b[0])/(cmath.pi/2)) % 4
        kw = round(cmath.phase(a[1]/b[1])/(2*cmath.pi/3)) % 3
        s1 = round((a[2]/b[2]).real) if abs((a[2]/b[2]).imag) < 1e-6 else 0
        s2 = round((a[3]/b[3]).real) if abs((a[3]/b[3]).imag) < 1e-6 else 0
        s3 = (round((a[4]/b[4]).real) if abs((a[4]/b[4]).imag) < 1e-6 else 0) if use_v else 1
        sy = round(((a[1]*a[4])/(b[1]*b[4])).real) if use_v and kz == 0 else 1
        out.append(dict(kz=kz, kw=kw, sp=s1, sq=s2, sv=s3, sy=sy))
    return out


def reduce_gram(G):
    """Pairwise reduction of a positive Gram matrix (returns G', unimodular T)."""
    n = G.shape[0]
    T = sp.eye(n)
    G = G.copy()
    changed = True
    while changed:
        changed = False
        for i in range(n):
            for j in range(n):
                if i != j:
                    f = sp.floor(sp.Rational(G[i, j], G[j, j])+sp.Rational(1, 2))
                    if f != 0:
                        E = sp.eye(n); E[j, i] = -f
                        G2 = E.T*G*E
                        if G2[i, i] < G[i, i]:
                            G, T, changed = G2, T*E, True
    return G, T


def short_vectors(G, bound, box=2):
    Gr, _ = reduce_gram(G)
    Gn = np.array(Gr.tolist(), dtype=np.int64)
    n = Gn.shape[0]
    rng = np.arange(-box, box+1)
    grid = np.array(np.meshgrid(*([rng]*n), indexing="ij")).reshape(n, -1).T
    q = np.einsum("ki,ij,kj->k", grid, Gn, grid)
    res = {}
    for val in sorted(set(q.tolist())):
        if 0 < val <= bound:
            res[int(val)] = int((q == val).sum())
    return res


def analyse(line, use_v, verbose=True):
    pts, x, y = build_cover(line, use_v)
    S = Surface(x, y)
    auts = automorphisms(x, y)
    n, N = len(pts), 2*S.genus
    lab = labels(pts, use_v)
    M = S.M
    rho = [S.rep(a) for a in auts]
    for r in rho:
        assert r.T*M*r == M
    assert S.rep(mul(auts[5], auts[7])) == rho[5]*rho[7]
    V4 = [j for j in range(n) if lab[j]["kz"] == 0 and lab[j]["kw"] == 0 and lab[j]["sv"] == 1]
    assert len(V4) == 4
    eps = [j for j in range(n) if lab[j]["kz"] == 0 and lab[j]["kw"] == 0 and lab[j]["sv"] == -1
           and lab[j]["sp"] == 1 and lab[j]["sq"] == 1]
    A4 = [j for j in range(n) if lab[j]["kz"] == 0 and lab[j]["sy"] == 1]
    assert len(A4) == 12
    sig = next(j for j in A4 if lab[j]["kw"] == 1)
    assert perm_order(auts[sig]) == 3
    quarter = [j for j in range(n) if lab[j]["kz"] == 1]
    one = sp.eye(N)
    epsm = [rho[eps[0]]-one] if use_v else []
    chars = {"p": lambda j: lab[j]["sp"], "q": lambda j: lab[j]["sq"], "pq": lambda j: lab[j]["sp"]*lab[j]["sq"]}
    L = {}
    for name, ch in chars.items():
        L[name] = lattice(epsm+[rho[j]-ch(j)*one for j in V4], N)
    LC = lattice(epsm+[rho[j]-one for j in V4], N)
    sumV4 = sum((rho[j] for j in V4), sp.zeros(N))
    LamP = lattice(epsm+[sumV4], N)
    out = dict(line=line, use_v=use_v, genus=S.genus, group=n)
    names = ["p", "q", "pq"]
    if use_v:
        L["m"] = lattice([rho[j]-one for j in A4], N)
        sumA4 = sum((rho[j] for j in A4), sp.zeros(N))
        LamA = lattice([24*one-3*(one+rho[eps[0]])*(4*one-sumV4)-2*sumA4], N)
        names.append("m")
    form = lambda B: B.T*M*B
    out["rank"] = {k: L[k].shape[1] for k in names}
    out["form_divisors"] = {k: smith(form(L[k])) for k in names}
    out["base_curve_part"] = dict(rank=LC.shape[1], divisors=smith(form(LC)))
    # sigma permutes the three Prym lattices
    perm = {}
    for k in ("p", "q", "pq"):
        img = rho[sig]*L[k]
        for k2 in ("p", "q", "pq"):
            try:
                coords(L[k2], img); perm[k] = k2
            except AssertionError:
                pass
    out["sigma_on_pryms"] = perm
    # complex structure from automorphisms acting with square -1
    J = {}
    for k in names:
        found = None
        for j in range(n):
            img = rho[j]*L[k]
            try:
                X = coords(L[k], img)
            except AssertionError:
                continue
            if X*X == -I2:
                found = (j, X); break
        if found is None:
            J[k] = None
            continue
        j, X = found
        Q = form(L[k])*X
        assert Q == Q.T
        if Q[0, 0] < 0:
            X = -X
        J[k] = (j, X)
    out["J_from_quarter_class"] = {k: (None if J[k] is None else lab[J[k][0]]["kz"]) for k in names}

    def glue(Lam, ks, tag):
        Bsum = sp.Matrix.hstack(*[L[k] for k in ks])
        T = coords(Lam, Bsum)
        res = dict(rank=Lam.shape[1], index=abs(T.det()), glue_group=[d for d in smith(T) if d != 1],
                   form_divisors=smith(form(Lam)))
        if all(J[k] is not None for k in ks):
            Jsum = sp.diag(*[J[k][1] for k in ks])
            JL = T*Jsum*T.inv()
            res["J_integral"] = all(v.is_integer for v in JL)
            if res["J_integral"]:
                Q = form(Lam)*JL
                assert Q == Q.T and Q.is_positive_definite
                res["Q_divisors"] = smith(Q)
                res["Q_even"] = all(Q[i, i] % 2 == 0 for i in range(Q.shape[0]))
                res["Q_short_vectors"] = short_vectors(Q, 4*min(Q[i, i] for i in range(Q.shape[0])))
            # glue vectors: coordinates of the lattice basis in the torus bases, mod 1
            Ti = T.inv()
            res["glue_vectors"] = sorted({tuple(str(sp.nsimplify(v % 1)) for v in Ti[:, c]) for c in range(Ti.shape[1])})
        out[tag] = res

    glue(LamP, ["p", "q", "pq"], "three_pryms")
    if use_v:
        glue(LamA, names, "four_tori")
    out["_obj"] = dict(S=S, pts=pts, x=x, y=y, auts=auts, lab=lab, rho=rho, L=L, J=J, sig=sig,
                       quarter=quarter, V4=V4, eps=eps, A4=A4, LamP=LamP)
    return out


def show(o):
    for k, v in o.items():
        if not k.startswith("_"):
            print("   ", k, ":", v)


# =============================================================== G. gluing test
print("G1. recorded kernel cover Y_q (line P), genus 5")
yq = analyse("P", False); show(yq)
print("G2. common cover X17 of Y_q and the mark torus, genus 17")
x17 = analyse("P", True); show(x17)
print("G3. comparison lines (same group data, other sign lines)")
for ln in ("Q", "L0"):
    o = analyse(ln, False)
    print("   line", ln, ": Prym forms", o["form_divisors"], "| J from automorphism:",
          o["J_from_quarter_class"], "| three-Prym glue", o["three_pryms"]["glue_group"])

# reference: recorded Gaussian E8 = Z^8 with Hamming code, complex pairs (0,1),(2,3),(4,5),(6,7)
print("G4. reference (recorded Gaussian E8): each coordinate torus has form 2 x unit, index 16,")
print("    glue group (Z/2)^4 = Hamming code of minimum weight 4, quadratic form even unimodular")

# verdict of the literal test
ft = x17["four_tori"]
d = x17["form_divisors"]
equal_forms = len({tuple(v) for v in d.values()}) == 1
print("G5. literal test: equal forms on the four tori:", equal_forms,
      "| index", ft["index"], "| glue", ft["glue_group"])


# ================================================== R. second pasted note
print("R1. the two invariant positive forms")
G1 = sp.eye(4)
G2 = sp.Rational(2, 3)*sp.eye(4)+sp.Rational(1, 3)*sp.ones(4)
for p in itertools.permutations(range(4)):
    P = sp.zeros(4)
    for a in range(4):
        P[p[a], a] = 1
    assert P.T*G1*P == G1 and P.T*G2*P == G2
assert G2.is_positive_definite and all(G2[k, k] == 1 for k in range(4))
c = sp.symbols("c")
Gc = (1-c)*sp.eye(4)+c*sp.ones(4)
print("    G2 eigenvalues", G2.eigenvals(), "| family (1-c)I + c 11^T eigenvalues", Gc.eigenvals())
print("    positive semidefinite exactly for -1/3 <= c <= 1; rank 3 at c = -1/3 (the four sources sum to zero)")
# four vectors with full symmetry in a three-dimensional space: e_k minus mean
E4 = sp.eye(4)-sp.ones(4)/4
Gram3 = E4.T*E4
Gram3n = Gram3/Gram3[0, 0]
assert Gram3n == Gc.subs(c, -sp.Rational(1, 3)) and Gram3.rank() == 3
print("    sum-zero (residue) space: normalised cross entry", Gram3n[0, 1], "squared", Gram3n[0, 1]**2)

print("R2. deck map against boundary reflection")
Dk, Rf = sp.diag(-1, -1), sp.diag(1, -1)
lam = sp.symbols("lam")
print("    char. polynomials", sp.factor(Dk.charpoly(lam).as_expr()), "|", sp.factor(Rf.charpoly(lam).as_expr()),
      "| det", Dk.det(), Rf.det())


def real_structure(o, line):
    ob = o["_obj"]
    pts, x, y, auts, S = ob["pts"], ob["x"], ob["y"], ob["auts"], ob["S"]
    n = len(pts)
    R = []
    for pt in pts:
        f = radicands(line, pt[1])
        img = [pt[0].conjugate(), pt[1].conjugate(), (pt[2]*pt[3]/f[1]).conjugate(), pt[3].conjugate()]
        if len(pt) == 5:
            img.append(pt[4].conjugate())
        R.append(match(tuple(img), pts))
    ident = list(range(n))
    assert mul(R, R) == ident
    assert mul(R, mul(x, R)) == inv(x) and mul(R, mul(y, R)) == inv(y)
    return R


def ovals(S, Rp):
    """Fixed circles of an antiholomorphic involution Rp and whether it separates."""
    n = S.n
    fixed = [e for e in range(n) if Rp[e] == e]
    fixed += [n+e for e in range(n) if Rp[e] == S.x[e]]
    fixed += [2*n+e for e in range(n) if S.y[Rp[e]] == e]
    par = list(range(S.nV))
    def find(a):
        while par[a] != a:
            par[a] = par[par[a]]; a = par[a]
        return a
    verts = set()
    for ed in fixed:
        a, b = find(S.tail[ed]), find(S.head[ed])
        verts |= {S.tail[ed], S.head[ed]}
        par[a] = b
    comps = len({find(v) for v in verts})
    fs = set(fixed)
    tp = list(range(2*n))
    def tf(a):
        while tp[a] != a:
            tp[a] = tp[tp[a]]; a = tp[a]
        return a
    for ed in range(3*n):
        if ed not in fs:
            a, b = S.duals[ed]
            tp[tf(a)] = tf(b)
    pieces = len({tf(k) for k in range(2*n)})
    return comps, pieces


print("R3. real structure R(w,a,b,c) = (conj w, conj a, conj c, conj b) on the genus 5 curve")
for ln, o in (("P", yq),):
    ob = o["_obj"]
    R = real_structure(o, ln)
    auts, quarter, S = ob["auts"], ob["quarter"], ob["S"]
    good = [j for j in quarter if mul(R, mul(auts[j], R)) == inv(auts[j])]
    print("    line", ln, ": R^2 = 1, R x R = x^-1, R y R = y^-1 hold; quarter lifts with R U R = U^-1:",
          len(good), "of", len(quarter), "| orders of those lifts", sorted({perm_order(auts[j]) for j in good}))
    kinds = {}
    for j in range(len(auts)):
        Rp = mul(R, auts[j])
        if mul(Rp, Rp) == list(range(S.n)):
            k = ovals(S, Rp)
            kinds[k] = kinds.get(k, 0)+1
    print("    R itself: (fixed circles, pieces after cutting) =", ovals(S, R))
    print("    all antiholomorphic involutions R.a of the curve, (circles, pieces): count =", kinds)
    Rm = S.rep_real(R)
    assert Rm*Rm == sp.eye(2*S.genus) and Rm.T*S.M*Rm == -S.M
    print("    on H^1: R is an involution and reverses the intersection form (orientation reversing)")

print("R4. flavour transition matrix on the holomorphic three-block of Y_q")
ob = yq["_obj"]
Lk, Jk, rho, S = ob["L"], ob["J"], ob["rho"], ob["S"]
Bsum = sp.Matrix.hstack(Lk["p"], Lk["q"], Lk["pq"])
Jsum = sp.diag(Jk["p"][1], Jk["q"][1], Jk["pq"][1])
Gsum = np.array((Bsum.T*S.M*Bsum).tolist(), dtype=float)
Jn = np.array(Jsum.tolist(), dtype=float)
real = lambda j: np.array(coords(Bsum, rho[j]*Bsum).tolist(), dtype=float)
for j in ob["A4"]+ob["quarter"]:
    assert np.allclose(real(j) @ Jn, Jn @ real(j)), "automorphism must be holomorphic"
Pp = (np.eye(6)-1j*Jn)/2                                # +i eigenspace of J
basis = Pp[:, [0, 2, 4]]
herm = lambda u, v: (np.conj(u) @ (Gsum @ Jn) @ v)       # positive on the block
def block(j):
    return np.linalg.lstsq(basis, real(j) @ basis, rcond=None)[0]
Hm = np.array([[herm(basis[:, a], basis[:, b]) for b in range(3)] for a in range(3)])
assert np.all(np.linalg.eigvalsh((Hm+Hm.conj().T)/2) > 0)
sg = block(ob["sig"])
ev, vec = np.linalg.eig(sg)
order = np.argsort(np.round(np.angle(ev)/(2*np.pi/3)) % 3)
vec = vec[:, order]
tables = {}
for j in ob["quarter"]:
    U = block(j)
    T = np.zeros((3, 3))
    for a in range(3):
        for b in range(3):
            num = abs(np.conj(vec[:, a]) @ Hm @ U @ vec[:, b])**2
            den = (np.conj(vec[:, a]) @ Hm @ vec[:, a]).real*(np.conj(vec[:, b]) @ Hm @ vec[:, b]).real
            T[a, b] = num/den
    key = tuple(tuple(sp.nsimplify(round(v, 9), rational=True) for v in row) for row in T)
    tables[key] = tables.get(key, 0)+1
print("    eigenvalues of the three-cycle on the block:", np.round(ev[order], 6))
for t in sorted(tables, key=str):
    print("    transition table of", tables[t], "of the", len(ob["quarter"]), "quarter lifts:", [[str(v) for v in row] for row in t])

# ============================== G6. fine structure of the three-Prym lattice
print("G6. fine structure inside H^1(Y_q, Z)")
from sympy.matrices.normalforms import hermite_normal_form as _hnf
hermite_normal_form = lambda Mx: _hnf(sp.Matrix(Mx.shape[0], Mx.shape[1], [int(v) for v in Mx]))
ob = yq["_obj"]
Lk, Jk, rho, S, LamP = ob["L"], ob["J"], ob["rho"], ob["S"], ob["LamP"]
Bsum = sp.Matrix.hstack(Lk["p"], Lk["q"], Lk["pq"])
Jsum = sp.diag(Jk["p"][1], Jk["q"][1], Jk["pq"][1])
T = coords(LamP, Bsum)
Ti2 = 2*T.inv()                                   # lattice basis in torus coordinates, times 2
assert all(v.is_integer for v in Ti2)
hull = hermite_normal_form(sp.Matrix.hstack(Ti2, Jsum*Ti2))
idx_hull = abs(Ti2.det())/abs(hull.det())
print("    [Lambda + J Lambda : Lambda] =", idx_hull, "= [Lambda : Lambda cap J Lambda];",
      "hull equals half the torus sum:", abs(hull.det()) == 1)
Qsum = (Bsum.T*S.M*Bsum)*Jsum
print("    quadratic form on the torus sum:", [Qsum[i, i] for i in range(6)], "diagonal:", Qsum.is_diagonal() or smith(Qsum))
QL = (T.inv().T*Qsum*T.inv())
assert all(v.is_integer for v in QL)
print("    E(x, Jx) on Lambda is integral although J is not: divisors", smith(QL),
      "| short vectors", short_vectors(QL, 4), "| even:", all(QL[i, i] % 2 == 0 for i in range(6)))
# four mark tori: invariants of the four subgroups of order three
A4 = ob["A4"]
threes = [j for j in A4 if perm_order(ob["auts"][j]) == 3]
subs = {}
for j in threes:
    key = frozenset({j, next(k for k in threes if mul(ob["auts"][j], ob["auts"][j]) == ob["auts"][k])})
    subs[key] = j
assert len(subs) == 4
N = 2*S.genus
Lm = [lattice([rho[j]-sp.eye(N)], N) for j in subs.values()]
print("    four tori Y_q / Z3 (one per subgroup of order three): ranks", [B.shape[1] for B in Lm],
      "| forms", [smith(B.T*S.M*B) for B in Lm])
Bm = sp.Matrix.hstack(*Lm)
Xm = sp.Matrix.hstack(*[coords(LamP, B) for B in Lm])       # they lie in the Prym lattice
gens = hermite_normal_form(Xm)
print("    their sum inside the three-Prym lattice: rank", Xm.rank(), "| index", abs(gens.det()),
      "| quotient", [d for d in smith(gens) if d != 1])
# complex structure on each mark torus, and integrality on their sum
Jreal = T*Jsum*T.inv()                                      # J in the basis of Lambda (rational)
stable = [all(v.is_integer for v in (coords(LamP, B).T*coords(LamP, B)).inv()*coords(LamP, B).T*Jreal*coords(LamP, B)) and
          coords(LamP, B)*((coords(LamP, B).T*coords(LamP, B)).inv()*coords(LamP, B).T*Jreal*coords(LamP, B)) == Jreal*coords(LamP, B)
          for B in Lm]
print("    each of these four tori is stable under J (so it is the square torus):", stable)
Xg = gens
Jg = Xg.inv()*Jreal*Xg
print("    J integral on the sum of the four tori:", all(v.is_integer for v in Jg))
if all(v.is_integer for v in Jg):
    Qg = Xg.T*QL*Xg
    assert Qg == Qg.T and Qg.is_positive_definite
    print("    quadratic form on that sum: divisors", smith(Qg), "| even:", all(Qg[i, i] % 2 == 0 for i in range(6)),
          "| short vectors", short_vectors(Qg, 12))
# normalised overlaps of the four mark lines in the holomorphic block
Tin = np.array(T.inv().tolist(), dtype=float)
vecs = []
for B in Lm:
    u = Tin @ np.array(coords(LamP, B).tolist(), dtype=float)[:, 0]   # torus-sum coordinates
    vecs.append(Pp @ u)
ov = [[abs(herm(a, b))**2/(herm(a, a).real*herm(b, b).real) for b in vecs] for a in vecs]
print("    squared overlaps of the four tetrahedral lines:", [[str(sp.nsimplify(round(v, 9), rational=True)) for v in row] for row in ov])

# pure S4 statement: the table depends only on the induced permutation of the four marks
print("R5. the same tables from the three-dimensional representation of S4 alone")
def pm(p):
    P = np.zeros((4, 4))
    for a in range(4):
        P[p[a], a] = 1
    return P
sgn = lambda p: round(np.linalg.det(pm(p)))
c3 = pm((1, 2, 0, 3))
w, v4 = np.linalg.eig(c3)
keep = [k for k in range(4) if abs(v4[:, k].sum()) < 1e-9 or abs(w[k]-1) > 1e-9]
z0 = np.ones(4)/2
e1 = np.array([1, 1, 1, -3])/np.sqrt(12)
lines = [e1]+[v4[:, k]/np.linalg.norm(v4[:, k]) for k in range(4) if abs(w[k]-1) > 1e-9]
res = {}
for p in itertools.permutations(range(4)):
    if sgn(p) == -1:
        Tt = tuple(sorted(round(abs(np.conj(a) @ pm(p) @ b)**2, 6) for a in lines for b in lines))
        res[Tt] = res.get(Tt, 0)+1
for k, cnt in res.items():
    print("    odd permutations:", cnt, "give entries", sorted(set(k)))

print("R6. the reflection in the circle through the four marks")
OMc = cmath.exp(2j*cmath.pi/3)
for zz in (0.3+0.4j, -1.1+0.2j, 0.7-0.9j):
    W = (zz*zz-1)/(zz*zz+1)
    zr = 1/zz.conjugate()
    Wr = (zr*zr-1)/(zr*zr+1)
    assert abs(Wr+W.conjugate()) < 1e-12              # w^3 -> -conj(w^3), so w -> -conj(w) up to a cube root of one
for mk in (1, 1j, -1, -1j):
    assert abs(1/complex(mk).conjugate()-mk) < 1e-15   # fixes all four marks
    assert sum(abs(complex(mk).conjugate()*u-mk) < 1e-15 for u in (1, 1j, -1, -1j)) == 1
wt = 0.37+0.21j
P_img = sorted([(-wt.conjugate()+OMc**k) for k in range(3)], key=lambda v: (round(v.real, 9), round(v.imag, 9)))
Q_ref = sorted([-(wt-OMc**k).conjugate() for k in range(3)], key=lambda v: (round(v.real, 9), round(v.imag, 9)))
assert all(abs(a-b) < 1e-12 for a, b in zip(P_img, Q_ref))
print("    z -> 1/conj(z) fixes all four marks and sends w -> -conj(w): the radicands of line P go to those of line Q")
print("    z -> i^k conj(z) (the reflections that act on Y_q) fix two marks (k even) or none (k odd)")

print("ALL CHECKS PASSED")
