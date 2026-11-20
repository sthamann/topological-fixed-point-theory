"""Computation 1: exhaustive enumeration of chiral anomaly-free 'ledgers'.
Class: gauge group SU(N)xSU(M)xU(1) (SU(1)=trivial; also the U(1)-less variant),
Weyl multiplets (rN,rM)_q with rN,rM in {singlet, fundamental, antifundamental}, integer q, |q|<=Q.
Conditions: SU(n)^3, Witten SU(2), SU(n)^2 U(1), grav^2 U(1), U(1)^3; chirality (no gauge-invariant
mass bilinear); irreducibility (no proper anomaly-free sub-ledger); every nonabelian factor acts.
Part A: minimal number of Weyl states per group (any number of multiplets), S<=SMAX.
Part B: all RIGID ledgers (anomaly variety has an isolated point up to scale => charges forced);
        rigid => multiplets k<=5 (Jacobian rank<=4), so enumerate ALL k<=5 for N<=7.
The SM is never put in: it must show up (or not) as an output."""
import itertools, sys, time
import numpy as np
from math import gcd
from functools import reduce, lru_cache

def types_for(N, M):
    rn = [0] if N == 1 else ([0, 1] if N == 2 else [0, 1, -1])
    rm = [0] if M == 1 else ([0, 1] if M == 2 else [0, 1, -1])
    ts = [(a, b) for a in rn for b in rm]
    key = lambda t: (0 if (t[0] and t[1]) else 1 if t[0] else 2 if t[1] else 3)
    return sorted(ts, key=key)

def tdim(t, N, M): return (N if t[0] else 1) * (M if t[1] else 1)
def conj(t, N, M): return ((-t[0] if N >= 3 else t[0]), (-t[1] if M >= 3 else t[1]))
def ndoublets(t, N, M): return (1 if (t[0] and N == 2) else 0) + (1 if (t[1] and M == 2) else 0)
def name(t, N, M):
    f = lambda r, n: '1' if r == 0 else (str(n) if r == 1 else str(n) + 'b')
    return '(%s,%s)' % (f(t[0], N), f(t[1], M))

def type_ok(counts, ts, N, M):
    for ax, n, m in ((0, N, M), (1, M, N)):
        if n == 1: continue
        other = lambda t: (m if t[1 - ax] else 1)
        if sum(c for c, t in zip(counts, ts) if t[ax]) == 0: return False
        if n >= 3 and sum(c * t[ax] * other(t) for c, t in zip(counts, ts)) != 0: return False
        if n == 2 and sum(c * other(t) for c, t in zip(counts, ts) if t[ax]) % 2: return False
    return True

@lru_cache(None)
def multisets(Q, m, selfconj, real):
    a = np.array(list(itertools.combinations_with_replacement(range(-Q, Q + 1), m)), dtype=np.int32).reshape(-1, m)
    if selfconj:
        keep = np.ones(len(a), bool)
        for i in range(m):
            for j in range(i + 1, m):
                keep &= (a[:, i] + a[:, j] != 0)
        if real: keep &= (a != 0).all(axis=1)
        a = a[keep]
    return a

def solve(N, M, counts, ts, Q):
    """return canonical primitive chiral solutions for this multiplet content"""
    mult = []  # expanded multiplet types
    for c, t in zip(counts, ts): mult += [t] * c
    k = len(mult)
    d = np.array([tdim(t, N, M) for t in mult], dtype=np.int64)
    cA = np.array([(M if t[1] else 1) if (t[0] and N >= 2) else 0 for t in mult], dtype=np.int64)
    cB = np.array([(N if t[0] else 1) if (t[1] and M >= 2) else 0 for t in mult], dtype=np.int64)
    rows = np.zeros((1, 0), dtype=np.int32)
    pos = 0
    nons = [(c, t) for c, t in zip(counts, ts) if c and t != (0, 0)]
    ns = sum(c for c, t in zip(counts, ts) if t == (0, 0))
    lastA = max([i for i, t in enumerate(mult) if t[0]], default=-1)
    lastB = max([i for i, t in enumerate(mult) if t[1]], default=-1)
    doneA = doneB = False
    for c, t in nons:
        arr = multisets(Q, c, conj(t, N, M) == t, ndoublets(t, N, M) % 2 == 0)
        out = []
        chunk = max(1, 4000000 // max(1, len(arr)))
        newpos = pos + c
        for s in range(0, len(rows), chunk):
            r = rows[s:s + chunk]
            x = np.concatenate([np.repeat(r, len(arr), axis=0), np.tile(arr, (len(r), 1))], axis=1)
            keep = np.ones(len(x), bool)
            for i in range(pos):            # cross-type chirality
                if conj(mult[i], N, M) == t and mult[i] != t:
                    for j in range(pos, newpos): keep &= (x[:, i] + x[:, j] != 0)
            if lastA < newpos and not doneA and lastA >= 0: keep &= (x @ cA[:newpos] == 0)
            if lastB < newpos and not doneB and lastB >= 0: keep &= (x @ cB[:newpos] == 0)
            out.append(x[keep])
        rows = np.concatenate(out) if out else rows[:0]
        pos = newpos
        if lastA < pos: doneA = True
        if lastB < pos: doneB = True
        if len(rows) == 0: return [], mult
    g = rows.astype(np.int64) @ d[:pos]
    cu = (rows.astype(np.int64) ** 3) @ d[:pos]
    if ns == 0:
        sol = rows[(g == 0) & (cu == 0)]
    else:
        S = multisets(Q, ns, True, True).astype(np.int64)
        s1 = S.sum(1); s3 = (S ** 3).sum(1)
        BIG = 10 ** 9
        key = s1 * BIG + s3
        o = np.argsort(key, kind='stable'); key = key[o]; S = S[o]
        need = (-g) * BIG + (-cu)
        lo = np.searchsorted(key, need, 'left'); hi = np.searchsorted(key, need, 'right')
        idx = np.nonzero(hi > lo)[0]
        parts = []
        for i in idx:
            blk = S[lo[i]:hi[i]]
            parts.append(np.concatenate([np.repeat(rows[i:i + 1], len(blk), axis=0), blk.astype(np.int32)], axis=1))
        sol = np.concatenate(parts) if parts else np.zeros((0, k), dtype=np.int32)
    # canonicalise
    bounds = []; p = 0
    for c, t in zip(counts, ts):
        if c: bounds.append((p, p + c)); p += c
    res = set()
    for q in sol.tolist():
        gg = reduce(gcd, [abs(v) for v in q])
        if gg == 0:
            if Q != 0: continue
        else:
            q = [v // gg for v in q]
        a = tuple(v for lo_, hi_ in bounds for v in sorted(q[lo_:hi_]))
        b = tuple(v for lo_, hi_ in bounds for v in sorted(-u for u in q[lo_:hi_]))
        res.add(max(a, b))
    return sorted(res), mult

def subset_free(q, mult, N, M, idx):
    for ax, n, m in ((0, N, M), (1, M, N)):
        if n == 1: continue
        oth = lambda t: (m if t[1 - ax] else 1)
        if n >= 3 and sum(mult[i][ax] * oth(mult[i]) for i in idx) != 0: return False
        if n == 2 and sum(oth(mult[i]) for i in idx if mult[i][ax]) % 2: return False
        if sum(oth(mult[i]) * q[i] for i in idx if mult[i][ax]) != 0: return False
    if sum(tdim(mult[i], N, M) * q[i] for i in idx) != 0: return False
    if sum(tdim(mult[i], N, M) * q[i] ** 3 for i in idx) != 0: return False
    return True

def irreducible(q, mult, N, M):
    k = len(mult)
    for r in range(1, k):
        for idx in itertools.combinations(range(k), r):
            if subset_free(q, mult, N, M, idx): return False
    return True

def local_dim(q, mult, N, M):
    k = len(mult)
    d = [tdim(t, N, M) for t in mult]
    J = [[(M if t[1] else 1) if (t[0] and N >= 2) else 0 for t in mult],
         [(N if t[0] else 1) if (t[1] and M >= 2) else 0 for t in mult],
         d, [3 * di * qi * qi for di, qi in zip(d, q)]]
    return k - np.linalg.matrix_rank(np.array(J, dtype=float))

def fmt(q, mult, N, M):
    return ' + '.join('%s_%d' % (name(t, N, M), v) for t, v in zip(mult, q))

def count_vectors_states(ts, N, M, S, maxsing):
    dims = [tdim(t, N, M) for t in ts]
    def rec(i, rem):
        if i == len(ts):
            if rem == 0: yield ()
            return
        top = rem // dims[i]
        if ts[i] == (0, 0): top = min(top, maxsing)
        for c in range(top + 1):
            for rest in rec(i + 1, rem - c * dims[i]): yield (c,) + rest
    return rec(0, S)

def count_vectors_k(ts, k):
    def rec(i, rem):
        if i == len(ts):
            if rem == 0: yield ()
            return
        for c in range(rem + 1):
            for rest in rec(i + 1, rem - c): yield (c,) + rest
    return rec(0, k)

def gname(N, M, Q): return 'SU(%d)xSU(%d)%s' % (N, M, 'xU(1)' if Q else ' [no U(1)]')

def partA(Q, SMAX, NMAX, maxsing, extra=0):
    print('\n=== PART A: minimal chiral anomaly-free irreducible ledger per group; |q|<=%d, states<=%d, singlets<=%d ===' % (Q, SMAX, maxsing))
    for N in range(1, NMAX + 1):
        for M in range(1, N + 1):
            if Q == 0 and N == 1: continue
            ts = types_for(N, M); found = None; t0 = time.time()
            for S in range(1, SMAX + 1):
                sols = []
                for cv in count_vectors_states(ts, N, M, S, maxsing):
                    if not type_ok(cv, ts, N, M): continue
                    res, mult = solve(N, M, cv, ts, Q)
                    for q in res:
                        if irreducible(q, mult, N, M):
                            sols.append((q, mult, local_dim(q, mult, N, M)))
                if sols:
                    rig = [s for s in sols if s[2] == 1]
                    print('%-22s min states = %2d   #solutions(primitive, in box) = %d   #rigid = %d   [%.1fs]' % (gname(N, M, Q), S, len(sols), len(rig), time.time() - t0))
                    for q, mult, ld in (rig[:3] + [s for s in sols if s[2] != 1][:3]):
                        print('      %s   {multiplets=%d, local dim of solution set=%d%s}' % (fmt(q, mult, N, M), len(mult), ld - 1 if Q else 0, ', RIGID' if ld == 1 else ''))
                    if found is None: found = S
                    if S >= found + extra: break
            if found is None:
                print('%-22s none with <= %d states   [%.1fs]' % (gname(N, M, Q), SMAX, time.time() - t0))
            sys.stdout.flush()

def partB(Q, NMAX, KMAX):
    print('\n=== PART B: ALL rigid chiral irreducible ledgers with <= %d multiplets, N<=%d, |q|<=%d ===' % (KMAX, NMAX, Q))
    print('(rigid = solution set of the anomaly equations is locally just the scaling line; rank(Jacobian)<=4 => rigid needs k<=5)')
    allr = []; nonrig = 0
    for N in range(1, NMAX + 1):
        for M in range(1, N + 1):
            ts = types_for(N, M)
            for k in range(1, KMAX + 1):
                for cv in count_vectors_k(ts, k):
                    if not type_ok(cv, ts, N, M): continue
                    res, mult = solve(N, M, cv, ts, Q)
                    for q in res:
                        if not irreducible(q, mult, N, M): continue
                        if local_dim(q, mult, N, M) == 1:
                            allr.append((sum(tdim(t, N, M) for t in mult), N, M, fmt(q, mult, N, M)))
                        else: nonrig += 1
    allr.sort()
    print('rigid chiral ledgers found: %d   (non-rigid chiral ledgers with k<=%d in the box: %d)' % (len(allr), KMAX, nonrig))
    for S, N, M, s in allr: print('  states=%2d  %s : %s' % (S, gname(N, M, 1), s))

if __name__ == '__main__':
    t = time.time()
    partA(Q=10, SMAX=16, NMAX=5, maxsing=6, extra=0)
    print('\n--- control: same search with the U(1) removed (all charges 0) ---')
    partA(Q=0, SMAX=30, NMAX=5, maxsing=0)
    partB(Q=14, NMAX=7, KMAX=5)
    print('\ntotal time %.1fs' % (time.time() - t))
