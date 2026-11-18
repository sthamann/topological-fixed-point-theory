#!/usr/bin/env python3
"""Exact integral cohomology of the covers used in the session checks of
18 November 2026.

Curves (all over the t-sphere, t = z^4, branched only over 0, 1, infinity):
    z^4 = t,  w^3 = (z^2-1)/(z^2+1)                 base curve C, genus 2
    p^2 = f0 f1,  q^2 = f1 f2                       V4 cover given by a sign line
    v^2 = w                                         optional: adjoins y^2 = z^4-1
with f_k = w + omega^k (line 'P', the recorded kernel cover Y_q),
f_k = w - omega^k (line 'Q'), f_k = w^2 - omega^(2k) (line 'L0').

The monodromy permutations x (around 0) and y (around 1) are obtained by
numerical continuation of radicals; everything after that is exact integer
arithmetic on the triangulated surface (two triangles per edge of the
dessin, vertices ordered black < white < face centre):
    H^1(X, Z) by tree/cotree, cup form by Alexander-Whitney,
    representation of the deck group on H^1(X, Z).
"""
import cmath
from collections import deque
import numpy as np
import sympy as sp

OM = cmath.exp(2j*cmath.pi/3)


def radicands(line, w):
    if line == "P":
        return [w+OM**k for k in range(3)]
    if line == "Q":
        return [w-OM**k for k in range(3)]
    if line == "L0":
        return [w*w-OM**(2*k) for k in range(3)]
    raise ValueError(line)


def fibre(t, line, use_v):
    pts = []
    for kz in range(4):
        z = abs(t)**0.25*cmath.exp(1j*(cmath.phase(t)+2*cmath.pi*kz)/4)
        W = (z*z-1)/(z*z+1)
        for kw in range(3):
            w = abs(W)**(1/3)*cmath.exp(1j*(cmath.phase(W)+2*cmath.pi*kw)/3)
            f = radicands(line, w)
            p0, q0 = cmath.sqrt(f[0]*f[1]), cmath.sqrt(f[1]*f[2])
            for s1 in (1, -1):
                for s2 in (1, -1):
                    if use_v:
                        for s3 in (1, -1):
                            pts.append((z, w, s1*p0, s2*q0, s3*cmath.sqrt(w)))
                    else:
                        pts.append((z, w, s1*p0, s2*q0))
    return pts


def follow(pt, path, line, use_v):
    cur = list(pt)
    for t in path:
        z = cur[0]*(t/cur[0]**4)**0.25
        W = (z*z-1)/(z*z+1)
        w = cur[1]*(W/cur[1]**3)**(1/3)
        f = radicands(line, w)
        p = cur[2]*cmath.sqrt(f[0]*f[1]/cur[2]**2)
        q = cur[3]*cmath.sqrt(f[1]*f[2]/cur[3]**2)
        nxt = [z, w, p, q]
        if use_v:
            nxt.append(cur[4]*cmath.sqrt(w/cur[4]**2))
        cur = nxt
    return tuple(cur)


def match(pt, pts):
    d = [sum(abs(a-b) for a, b in zip(pt, q)) for q in pts]
    j = int(np.argmin(d))
    assert d[j] < 1e-6, d[j]
    return j


def build_cover(line="P", use_v=False, steps=4000):
    t0 = 0.5
    pts = fibre(t0, line, use_v)
    th = [2*cmath.pi*k/steps for k in range(1, steps+1)]
    loop0 = [0.5*cmath.exp(1j*a) for a in th]
    loop1 = [1-0.5*cmath.exp(1j*a) for a in th]
    x = [match(follow(p, loop0, line, use_v), pts) for p in pts]
    y = [match(follow(p, loop1, line, use_v), pts) for p in pts]
    assert sorted(x) == list(range(len(pts))) and sorted(y) == list(range(len(pts)))
    return pts, x, y


def inv(p):
    r = [0]*len(p)
    for i, j in enumerate(p):
        r[j] = i
    return r


def mul(a, b):                       # (a after b)
    return [a[k] for k in b]


def perm_order(p):
    n, q = 1, list(p)
    ident = list(range(len(p)))
    while q != ident:
        q, n = mul(p, q), n+1
    return n


def cycles_id(p):
    ident = [-1]*len(p)
    c = 0
    for s in range(len(p)):
        if ident[s] < 0:
            k = s
            while ident[k] < 0:
                ident[k] = c
                k = p[k]
            c += 1
    return ident, c


def automorphisms(x, y):
    n = len(x)
    auts = []
    for j in range(n):
        a = [-1]*n
        a[0] = j
        dq = deque([0])
        while dq:
            e = dq.popleft()
            for g in (x, y):
                e2, img = g[e], g[a[e]]
                if a[e2] < 0:
                    a[e2] = img
                    dq.append(e2)
                else:
                    assert a[e2] == img, "cover is not Galois"
        assert sorted(a) == list(range(n))
        auts.append(a)
    return auts


class Surface:
    """Triangulated surface of a dessin (x, y) with exact H^1 and cup form."""

    def __init__(self, x, y):
        n = self.n = len(x)
        self.x, self.y = list(x), list(y)
        self.xi, self.yi = inv(x), inv(y)
        phi = [self.xi[self.yi[e]] for e in range(n)]
        B, nB = cycles_id(self.x)
        W, nW = cycles_id(self.y)
        F, nF = cycles_id(phi)
        self.nV = nB+nW+nF
        self.genus = (2-(self.nV-3*n+2*n))//2
        assert 2-2*self.genus == self.nV-n
        E = 3*n
        tail = [B[e] for e in range(n)]+[B[e] for e in range(n)]+[nB+W[e] for e in range(n)]
        head = [nB+W[e] for e in range(n)]+[nB+nW+F[e] for e in range(n)]+[nB+nW+F[e] for e in range(n)]
        self.tail, self.head = tail, head
        # cocycle rows: triangle L(e) = e, R(e) = n + e
        D = np.zeros((2*n, E), dtype=np.int64)
        for e in range(n):
            D[e, 2*n+e] += 1; D[e, n+e] -= 1; D[e, e] += 1
            D[n+e, 2*n+self.y[e]] += 1; D[n+e, n+self.xi[e]] -= 1; D[n+e, e] += 1
        self.D = D
        duals = [None]*E
        for e in range(n):
            duals[e] = (e, n+e)
            duals[n+e] = (e, n+self.x[e])
            duals[2*n+e] = (e, n+self.yi[e])
        for ed in range(E):
            for nd in duals[ed]:
                assert D[nd, ed] != 0
        self.duals = duals
        # primal spanning tree
        adj = [[] for _ in range(self.nV)]
        for ed in range(E):
            adj[tail[ed]].append((ed, head[ed], +1))
            adj[head[ed]].append((ed, tail[ed], -1))
        par = {0: None}
        depth = {0: 0}
        dq = deque([0])
        intree = [False]*E
        while dq:
            u = dq.popleft()
            for ed, v, s in adj[u]:
                if v not in par:
                    par[v] = (u, ed, s)          # going u -> v traverses ed with sign s
                    depth[v] = depth[u]+1
                    intree[ed] = True
                    dq.append(v)
        assert len(par) == self.nV
        # dual spanning tree among non-tree edges
        dadj = [[] for _ in range(2*n)]
        for ed in range(E):
            if not intree[ed]:
                a, b = duals[ed]
                dadj[a].append((ed, b)); dadj[b].append((ed, a))
        dpar = {0: None}
        ddepth = {0: 0}
        dq = deque([0])
        incot = [False]*E
        while dq:
            u = dq.popleft()
            for ed, v in dadj[u]:
                if v not in dpar:
                    dpar[v] = (u, ed)
                    ddepth[v] = ddepth[u]+1
                    incot[ed] = True
                    dq.append(v)
        assert len(dpar) == 2*n
        left = [ed for ed in range(E) if not intree[ed] and not incot[ed]]
        assert len(left) == 2*self.genus, (len(left), self.genus)
        sgn = lambda ed, nd: int(D[nd, ed])*(1 if nd < n else -1)   # incidence form
        A = np.zeros((len(left), E), dtype=np.int64)
        C = np.zeros((len(left), E), dtype=np.int64)
        for j, ed in enumerate(left):
            n1, n2 = duals[ed]
            A[j, ed] = sgn(ed, n1)                  # leave n1 through ed
            a, b = n2, n1                           # dual tree path a -> b
            up_a, up_b = [], []
            while a != b:
                if ddepth[a] >= ddepth[b]:
                    p, e2 = dpar[a]; up_a.append((a, e2)); a = p
                else:
                    p, e2 = dpar[b]; up_b.append((b, e2)); b = p
            for node, e2 in up_a:                   # leave `node` upward through e2
                A[j, e2] = sgn(e2, node)
            for node, e2 in up_b:                   # arrive at `node` from its parent
                A[j, e2] = -sgn(e2, node)
            # primal cycle: ed from tail to head, then tree path head -> tail
            C[j, ed] += 1
            a, b = head[ed], tail[ed]
            while a != b:
                if depth[a] >= depth[b]:
                    p, e2, s = par[a]; C[j, e2] -= s; a = p     # a -> parent reverses (p -> a)
                else:
                    p, e2, s = par[b]; C[j, e2] += s; b = p     # eventually parent -> b
        assert not (D @ A.T).any(), "not cocycles"
        P = A @ C.T
        assert (np.abs(P) == np.eye(len(left), dtype=np.int64)).all(), "dual bases"
        A = (np.diag(P)[:, None])*A
        self.A, self.C = A, C
        yy = np.array(self.y)
        wf = A[:, 2*n:]
        M = A[:, :n] @ (wf-wf[:, yy]).T
        assert (M == -M.T).all()
        self.M = sp.Matrix(M.tolist())
        assert abs(self.M.det()) == 1, "cup form must be unimodular"

    def rep(self, a):
        """Matrix of the automorphism a (a permutation of the dessin edges
        commuting with x and y) on H^1(X, Z): columns = images of the basis."""
        n = self.n
        ai = np.array(inv(a))
        idx = np.concatenate([ai, n+ai, 2*n+ai])
        return sp.Matrix((self.C @ self.A[:, idx].T).tolist())

    def rep_real(self, R):
        """Matrix on H^1 of an antiholomorphic involution given on the fibre
        (R x R = x^-1, R y R = y^-1)."""
        n = self.n
        R = np.array(R)
        idx = np.concatenate([R, n+np.array(self.xi)[R], 2*n+np.array(self.y)[R]])
        return sp.Matrix((self.C @ self.A[:, idx].T).tolist())


# ------------------------------------------------------------ integer lattices
def int_kernel(rows, ncols):
    """Saturated Z-basis of {v in Z^ncols : rows . v = 0} (list of columns)."""
    m = len(rows)
    cols = [[int(rows[i][j]) for i in range(m)]+[int(j == k) for k in range(ncols)] for j in range(ncols)]
    r = 0
    for i in range(m):
        while True:
            nz = [j for j in range(r, ncols) if cols[j][i] != 0]
            if len(nz) <= 1:
                break
            piv = min(nz, key=lambda j: abs(cols[j][i]))
            for j in nz:
                if j != piv:
                    f = cols[j][i]//cols[piv][i]
                    if f:
                        cols[j] = [a-f*b for a, b in zip(cols[j], cols[piv])]
        nz = [j for j in range(r, ncols) if cols[j][i] != 0]
        if nz:
            cols[r], cols[nz[0]] = cols[nz[0]], cols[r]
            r += 1
    return [c[m:] for c in cols[r:]]


def lattice(mats, ncols):
    """Z^ncols intersected with the common kernel of the integer matrices."""
    rows = []
    for Mx in mats:
        rows += Mx.tolist()
    ker = int_kernel(rows, ncols)
    return sp.Matrix(ker).T if ker else sp.zeros(ncols, 0)


def coords(Bbig, Bsmall):
    """Integer matrix X with Bbig X = Bsmall (asserts integrality)."""
    X = (Bbig.T*Bbig).inv()*Bbig.T*Bsmall
    assert Bbig*X == Bsmall
    assert all(v.is_integer for v in X), "not a sublattice"
    return X


def smith(Mx):
    from sympy.matrices.normalforms import smith_normal_form
    S = smith_normal_form(Mx, domain=sp.ZZ)
    return sorted(abs(int(S[i, i])) for i in range(min(S.shape)))


# ------------------------------------------------- covers over the s-sphere
# s = -(t-1)^2/(4t) is invariant under t -> 1/t.  Unit circle |t| = 1 maps to
# the segment [0, 1], positive real t to s < 0, negative real t to s > 1, so
# both t -> conj(t) and t -> 1/conj(t) lie over s -> conj(s) and every fixed
# curve of an antiholomorphic involution lies in the one-skeleton.
import itertools as _it


def _troots(s):
    d = cmath.sqrt(s*s-s)
    return (1-2*s)+2*d, (1-2*s)-2*d


def fibre_s(s, lines):
    pts = []
    for t in _troots(s):
        for kz in range(4):
            z = abs(t)**0.25*cmath.exp(1j*(cmath.phase(t)+2*cmath.pi*kz)/4)
            W = (z*z-1)/(z*z+1)
            for kw in range(3):
                w = abs(W)**(1/3)*cmath.exp(1j*(cmath.phase(W)+2*cmath.pi*kw)/3)
                rads = []
                for ln in lines:
                    f = radicands(ln, w)
                    rads += [cmath.sqrt(f[0]*f[1]), cmath.sqrt(f[1]*f[2])]
                for sg in _it.product((1, -1), repeat=len(rads)):
                    pts.append((t, z, w)+tuple(a*b for a, b in zip(sg, rads)))
    return pts


def follow_s(pt, path, lines):
    cur = list(pt)
    for s in path:
        c1, c2 = _troots(s)
        t = c1 if abs(c1-cur[0]) < abs(c2-cur[0]) else c2
        z = cur[1]*(t/cur[1]**4)**0.25
        W = (z*z-1)/(z*z+1)
        w = cur[2]*(W/cur[2]**3)**(1/3)
        nxt = [t, z, w]
        k = 3
        for ln in lines:
            f = radicands(ln, w)
            nxt.append(cur[k]*cmath.sqrt(f[0]*f[1]/cur[k]**2))
            nxt.append(cur[k+1]*cmath.sqrt(f[1]*f[2]/cur[k+1]**2))
            k += 2
        cur = nxt
    return tuple(cur)


def build_cover_s(lines=("P", "Q"), steps=6000):
    s0 = 0.5
    pts = fibre_s(s0, lines)
    th = [2*cmath.pi*k/steps for k in range(1, steps+1)]
    loop0 = [0.5*cmath.exp(1j*a) for a in th]
    loop1 = [1-0.5*cmath.exp(1j*a) for a in th]
    arr = np.array(pts)
    def m(pt):
        d = np.abs(arr-np.array(pt)).sum(axis=1)
        j = int(np.argmin(d))
        assert d[j] < 1e-6, d[j]
        return j
    x = [m(follow_s(p, loop0, lines)) for p in pts]
    y = [m(follow_s(p, loop1, lines)) for p in pts]
    n = len(pts)
    assert sorted(x) == list(range(n)) and sorted(y) == list(range(n))
    return pts, x, y


def antiholomorphic_maps(x, y):
    """All maps R of the fibre with R x = x^-1 R and R y = y^-1 R."""
    n = len(x)
    xi, yi = inv(x), inv(y)
    out = []
    for j in range(n):
        a = [-1]*n
        a[0] = j
        dq = deque([0])
        ok = True
        while dq and ok:
            e = dq.popleft()
            for g, gi in ((x, xi), (y, yi)):
                e2, img = g[e], gi[a[e]]
                if a[e2] < 0:
                    a[e2] = img
                    dq.append(e2)
                elif a[e2] != img:
                    ok = False
                    break
        out.append(a if ok and sorted(a) == list(range(n)) else None)
    return out
