"""Consistent-world generators.  A world is an assignment of symbols 0..3 to cells;
it is valid when every active four-cell ("four-mark") constraint sums to zero in a
group G.  G = Z2xZ2 ('xor') is native to the quartic code; G = Z4 ('z4') has the same
shape but the wrong group and serves as the mismatched control.
"""
import itertools, os
import numpy as np
import sympy

# ---- five partitions of the 4x4 grid into quads (rule families)
def _families():
    c = lambda r, k: 4 * r + k
    R = [[c(r, k) for k in range(4)] for r in range(4)]
    C = [[c(r, k) for r in range(4)] for k in range(4)]
    B = [[c(2 * i + a, 2 * j + b) for a in (0, 1) for b in (0, 1)] for i in (0, 1) for j in (0, 1)]
    D = [[c(r, (r + k) % 4) for r in range(4)] for k in range(4)]
    A = [[c(r, (k - r) % 4) for r in range(4)] for k in range(4)]
    return np.array([R, C, B, D, A])

FAM = _families()                       # (5,4,4)
QUADS16 = FAM.reshape(20, 4)
NAMES = "RCBDA"
SINGLES = [(i,) for i in range(5)]
PAIRS = list(itertools.combinations(range(5), 2))
TRAIN_PAIRS = [(0, 1), (0, 2), (1, 3), (2, 4), (3, 4)]
TRAIN_COMBOS = SINGLES + TRAIN_PAIRS
OOD_COMBOS = [p for p in PAIRS if p not in TRAIN_PAIRS] + list(itertools.combinations(range(5), 3))


def gsum(vals, group, axis=-1):
    if group == "xor":
        return np.bitwise_xor.reduce(vals, axis=axis)
    return vals.sum(axis=axis) % 4


def _null_f2(M):
    M = M.copy() % 2
    m, n = M.shape
    piv, r = [], 0
    for c in range(n):
        rows = [i for i in range(r, m) if M[i, c]]
        if not rows:
            continue
        M[[r, rows[0]]] = M[[rows[0], r]]
        for i in range(m):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        piv.append(c); r += 1
        if r == m:
            break
    free = [c for c in range(n) if c not in piv]
    basis = []
    for f in free:
        v = np.zeros(n, dtype=np.int64); v[f] = 1
        for i, pc in enumerate(piv):
            v[pc] = M[i, f]
        basis.append(v)
    return np.array(basis).reshape(-1, n)


_cache = {}
def _kernel(combo, group, quads, n):
    key = (combo, group, n)
    if key in _cache:
        return _cache[key]
    M = np.zeros((len(quads), n), dtype=np.int64)
    for i, q in enumerate(quads):
        M[i, q] = 1
    f2 = _null_f2(M)
    zb = None
    if group == "z4":
        ns = sympy.Matrix(M.tolist()).nullspace()
        zb = []
        for v in ns:
            l = sympy.ilcm(*[sympy.Rational(x).q for x in v])
            zb.append(np.array([int(x * l) for x in v], dtype=np.int64))
        zb = np.array(zb).reshape(-1, n)
    _cache[key] = (f2, zb)
    return _cache[key]


def sample_world(combo, group, rng, quads=None, n=16):
    quads = FAM[list(combo)].reshape(-1, 4) if quads is None else quads
    f2, zb = _kernel(tuple(combo), group, quads, n)
    if group == "xor":   # uniform over the solution space, bit-plane wise
        bits = [(rng.integers(0, 2, len(f2)) @ f2) % 2 for _ in range(2)]
        return bits[0] + 2 * bits[1]
    w = rng.integers(0, 4, len(zb)) @ zb + 2 * (rng.integers(0, 2, len(f2)) @ f2)
    return w % 4


def n_solutions(world, masked, quads, group):
    k = len(masked)
    cand = np.array(list(itertools.product(range(4), repeat=k)), dtype=np.int64)
    W = np.repeat(world[None], len(cand), 0)
    W[:, masked] = cand
    return int((gsum(W[:, quads], group) == 0).all(1).sum())


def make_grid_set(n, combos, ks, group, seed):
    """4x4 rule-combination worlds.  Returns dict of arrays; every puzzle has a unique completion."""
    rng = np.random.default_rng(seed)
    X, K, R = [], [], []
    while len(X) < n:
        combo = combos[rng.integers(len(combos))]
        quads = FAM[list(combo)].reshape(-1, 4)
        w = sample_world(combo, group, rng)
        assert (gsum(w[quads], group) == 0).all()
        k = int(rng.choice(ks))
        for _ in range(6):
            masked = rng.choice(16, k, replace=False)
            if n_solutions(w, masked, quads, group) == 1:
                known = np.ones(16, bool); known[masked] = False
                rule = np.zeros(5, bool); rule[list(combo)] = True
                X.append(w); K.append(known); R.append(rule)
                break
    R = np.array(R)
    return dict(x=np.array(X), known=np.array(K), rule=R, active=np.repeat(R, 4, axis=1))


def cached(name, fn):
    os.makedirs("out/cache", exist_ok=True)
    path = f"out/cache/{name}.npz"
    if os.path.exists(path):
        return dict(np.load(path))
    d = fn()
    np.savez_compressed(path, **d)
    return d


# ---- plaquette lattice: every 2x2 window of an n x n grid is a four-mark cell
def lattice_quads(n):
    return np.array([[r * n + c, r * n + c + 1, (r + 1) * n + c, (r + 1) * n + c + 1]
                     for r in range(n - 1) for c in range(n - 1)])


def lattice_world(n, group, rng):
    a, b = rng.integers(0, 4, n), rng.integers(0, 4, n)
    if group == "xor":
        return (a[:, None] ^ b[None, :]).reshape(-1)
    sgn = (-1) ** (np.add.outer(np.arange(n), np.arange(n)))
    return ((sgn * (a[:, None] + b[None, :])) % 4).reshape(-1)


def propagate_depth(known, quads):
    """Rounds of local unit propagation (a plaquette with 3 known cells fixes the 4th).
    Returns (#rounds, solved?)."""
    known = known.copy()
    d = 0
    while not known.all():
        kq = known[quads]
        hit = kq.sum(1) == 3
        if not hit.any():
            return d, False
        known[quads[hit][~kq[hit]]] = True
        d += 1
    return d, True


def connected(known, n):
    """Is the completion unique up to nothing (bipartite row/column graph connected)?"""
    par = list(range(2 * n))
    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for i in np.flatnonzero(known):
        par[find(i // n)] = find(n + i % n)
    return len({find(x) for x in range(2 * n)}) == 1


def make_lattice_set(m, n, group, seed, kind="mix"):
    """kind: 'random' (random givens, locally solvable), 'chain' (one row + one column given:
    inference must travel across the whole grid), 'mix', or 'global' (uniquely determined
    but NOT solvable by local unit propagation)."""
    rng = np.random.default_rng(seed)
    quads = lattice_quads(n)
    X, K, D = [], [], []
    while len(X) < m:
        w = lattice_world(n, group, rng)
        kd = kind if kind != "mix" else ("chain" if rng.random() < 0.5 else "random")
        known = np.zeros(n * n, bool)
        if kd == "chain":
            r0, c0 = rng.integers(n), rng.integers(n)
            g = known.reshape(n, n); g[r0, :] = True; g[:, c0] = True
            known |= rng.random(n * n) < rng.uniform(0, 0.15)
        else:
            known = rng.random(n * n) < rng.uniform(0.25, 0.7)
        d, ok = propagate_depth(known, quads)
        if kd == "global":
            if ok or not connected(known, n):
                continue
        elif not ok or d == 0:
            continue
        X.append(w); K.append(known); D.append(d)
    return dict(x=np.array(X), known=np.array(K), depth=np.array(D))
