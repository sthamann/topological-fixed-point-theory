"""Computation 2: which Weyl-group orbits of weight-lattice points are 'chiral anomaly-free families'?
A finite set S of charge vectors (weights of the Cartan torus) is
  anomaly-free  iff  sum v = 0 and sum v(x)v(x)v = 0   (all linear and cubic abelian anomalies; = odd part of a spherical 3-design)
  chiral        iff  S and -S are disjoint (no state has its conjugate in the set).
Second moment is automatically isotropic for an irreducible reflection group (=> universal coupling normalisation).
Class enumerated: ALL irreducible root systems of rank<=8, ALL dominant weights with Dynkin labels <=CMAX, orbit size <=BOUND.
Exact integer arithmetic in Dynkin-label coordinates (vanishing of a tensor is basis independent).
Controls: (a) groups containing -1 must give nothing; (b) drop transitivity: unions of two orbits;
(c) drop lattice-rationality: non-crystallographic I2(m) and an irrational A4 orbit."""
import itertools, time
import numpy as np

def cartan(typ, n):
    A = 2 * np.eye(n, dtype=int)
    def link(i, j, a=-1, b=-1): A[i, j] = a; A[j, i] = b
    if typ in 'ABCF':
        for i in range(n - 1): link(i, i + 1)
        if typ == 'B': A[n - 2, n - 1] = -2
        if typ == 'C': A[n - 1, n - 2] = -2
        if typ == 'F': A[1, 2] = -2
    elif typ == 'D':
        for i in range(n - 2): link(i, i + 1)
        link(n - 3, n - 1)
    elif typ == 'E':
        for i in range(n - 2): link(i, i + 1)
        link(2, n - 1)
    elif typ == 'G':
        A = np.array([[2, -1], [-3, 2]])
    return A

def orbit(A, lam, bound):
    n = len(A); start = tuple(lam); seen = {start}; stack = [start]
    while stack:
        v = stack.pop()
        for i in range(n):
            if v[i]:
                w = tuple(v[j] - v[i] * A[i, j] for j in range(n))
                if w not in seen:
                    seen.add(w); stack.append(w)
                    if len(seen) > bound: return None
    return seen

def moments(S):
    O = np.array(sorted(S), dtype=np.int64)
    return O.sum(0), np.einsum('ni,nj,nk->ijk', O, O, O)

def systems(maxrank):
    out = []
    for n in range(1, maxrank + 1):
        out.append(('A', n))
        if n >= 2: out.append(('B', n))
        if n >= 3: out.append(('C', n))
        if n >= 4: out.append(('D', n))
        if n in (6, 7, 8): out.append(('E', n))
        if n == 4: out.append(('F', 4))
        if n == 2: out.append(('G', 2))
    return out

def dominant_weights(A, bound, cmax):
    n = len(A)
    for r in range(1, n + 1):
        for supp in itertools.combinations(range(n), r):
            lam = [0] * n
            for i in supp: lam[i] = 1
            if orbit(A, lam, bound) is None: continue
            for vals in itertools.product(range(1, cmax + 1), repeat=r):
                lam = [0] * n
                for i, v in zip(supp, vals): lam[i] = v
                yield tuple(lam)

if __name__ == '__main__':
    t0 = time.time(); BOUND = 300; CMAX = 6
    print('=== single orbits: rank<=8, labels<=%d, orbit size<=%d ===' % (CMAX, BOUND))
    hits = []; store = {}
    for typ, n in systems(8):
        A = cartan(typ, n); tot = chir = free = 0; mincubfree = None
        for lam in dominant_weights(A, BOUND, CMAX):
            S = orbit(A, lam, BOUND); tot += 1
            anti = tuple(-x for x in lam) in S
            m1, T = moments(S)
            assert not m1.any()
            cub = not T.any()
            if cub and (mincubfree is None or len(S) < mincubfree): mincubfree = len(S)
            if not anti: chir += 1
            if cub: free += 1
            if cub and not anti: hits.append((len(S), n, typ, lam))
            if n <= 6 and len(S) <= 40 and max(lam) <= 4: store.setdefault((typ, n), []).append((lam, S, T, anti))
        print('%s%d: orbits tested=%5d  chiral(antipode-free)=%5d  cubic-anomaly-free=%5d  smallest cubic-free orbit=%s' % (typ, n, tot, chir, free, mincubfree))
    hits.sort()
    print('\nchiral AND anomaly-free single orbits (sorted by size), total %d:' % len(hits))
    for h in hits[:25]: print('   size=%3d  rank=%d  %s%d  dominant weight (Dynkin labels)=%s' % (h[0], h[1], h[2], h[1], h[3]))

    print('\n=== control (b): unions of TWO distinct orbits of the same group, rank<=6, each <=40 points, labels<=4 ===')
    u = []
    for (typ, n), L in store.items():
        for i in range(len(L)):
            for j in range(i + 1, len(L)):
                l1, S1, T1, a1 = L[i]; l2, S2, T2, a2 = L[j]
                if a1 or a2 or not T1.any(): continue
                if (T1 + T2).any(): continue
                if any(tuple(-x for x in v) in S2 for v in S1): continue
                u.append((len(S1) + len(S2), n, typ, l1, l2, len(S1), len(S2)))
    u.sort()
    print('two-orbit chiral anomaly-free unions found: %d' % len(u))
    for x in u[:12]: print('   size=%3d rank=%d %s%d  weights %s (%d pts) + %s (%d pts)' % (x[0], x[1], x[2], x[1], x[3], x[5], x[4], x[6]))

    print('\n=== control (c): drop lattice rationality ===')
    for m in range(3, 10):
        th = 2 * np.pi * np.arange(m) / m; P = np.stack([np.cos(th), np.sin(th)], 1)
        T = np.einsum('ni,nj,nk->ijk', P, P, P)
        anti = any(np.allclose(-p, q) for p in P for q in P)
        print('   I2(%d) regular %d-gon in R^2: antipode-free=%s  cubic moment max=%.2e  => chiral anomaly-free=%s' % (m, m, not anti, abs(T).max(), (not anti) and abs(T).max() < 1e-9))
    r = (-3 + 5 ** 0.5) / 2; lam = np.array([r, r, 1.0, 1.0, -2 * (r + 1)])
    P = np.array(sorted(set(itertools.permutations(lam))))
    T = np.einsum('ni,nj,nk->ijk', P, P, P)
    print('   irrational A4 orbit of (r,r,1,1,-2r-2), r=(-3+sqrt5)/2: %d points, sum=%.1e, cubic max=%.1e, antipode-free=%s' % (len(P), abs(P.sum(0)).max(), abs(T).max(), not any(np.allclose(-p, q) for p in P for q in P)))

    print('\n=== check of the winner: D5 half-spinor as 16 SM-like states ===')
    S = [v for v in itertools.product([1, -1], repeat=5) if np.prod(v) == 1]
    P = np.array(S) / 2.0
    Y = lambda v: -(v[0] + v[1] + v[2]) / 3.0 + (v[3] + v[4]) / 2.0
    T3 = lambda v: (v[3] - v[4]) / 2.0
    Ys = sorted(Y(v) for v in P)
    print('   hypercharges Y=-(c1+c2+c3)/3+(w1+w2)/2 over the 16:', sorted(set(round(y, 4) for y in Ys)), ' multiplicities', [Ys.count(y) for y in sorted(set(Ys))])
    print('   Tr T3^2 = %.4f, Tr Q^2 = %.4f, ratio sin^2 = %.4f (3/8=0.375)' % (sum(T3(v) ** 2 for v in P), sum((T3(v) + Y(v)) ** 2 for v in P), sum(T3(v) ** 2 for v in P) / sum((T3(v) + Y(v)) ** 2 for v in P)))
    print('   second moment matrix = %s * identity' % (P.T @ P)[0, 0], ' offdiag max', abs(P.T @ P - np.diag(np.diag(P.T @ P))).max())
    print('total time %.1fs' % (time.time() - t0))
