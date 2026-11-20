"""H2: exhaustive search over ALL symmetric +-1 relations (Seidel matrices S, zero diagonal) on n<=8 points,
modulo the gauge 'switching' (S -> D S D, D diagonal +-1; fixed by first row = +1).
Rule tested:  S^2 = (n-1) I   ("perfectly balanced relation").
Wider control class: S^2 = a S + b I (exactly two eigenvalues), and the eigenvalue-count histogram."""
import numpy as np, itertools, math
from collections import Counter

def all_normalized(n):
    m = n-1; pairs = [(i, j) for i in range(1, n) for j in range(i+1, n)]
    P = len(pairs); N = 1 << P
    for start in range(0, N, 1 << 18):
        idx = np.arange(start, min(N, start + (1 << 18)), dtype=np.int64)
        S = np.zeros((len(idx), n, n), dtype=np.int16)
        S[:, 0, 1:] = 1; S[:, 1:, 0] = 1
        for b, (i, j) in enumerate(pairs):
            v = (1 - 2*((idx >> b) & 1)).astype(np.int16)
            S[:, i, j] = v; S[:, j, i] = v
        yield S

def canon(S):
    """canonical form of a switching class under permutations (brute force), n<=8 only for few matrices"""
    n = len(S); best = None
    for p in itertools.permutations(range(n)):
        T = S[np.ix_(p, p)]; d = T[0].copy(); d[0] = 1
        T = T * d[:, None] * d[None, :]
        key = T[np.triu_indices(n, 1)].tobytes()
        if best is None or key < best: best = key
    return best

print("n | #switching-gauge-fixed relations | S^2=(n-1)I | two-eigenvalue (S^2=aS+bI) | values of a | iso classes of conference solutions")
store = {}
for n in range(2, 9):
    tot = conf = two = 0; avals = Counter(); sols = []
    eye = np.eye(n, dtype=bool)
    for S in all_normalized(n):
        tot += len(S)
        S2 = np.matmul(S, S)
        off = S2 * S                      # off-diagonal: S2_ij * S_ij (= a if S2 = aS + bI)
        o = off[:, ~eye]
        is_two = (o == o[:, :1]).all(axis=1)
        is_conf = (S2[:, ~eye] == 0).all(axis=1)
        two += int(is_two.sum()); conf += int(is_conf.sum())
        for a in o[is_two, 0]: avals[int(a)] += 1
        for s in S[is_conf]: sols.append(s.copy())
    classes = len({canon(s) for s in sols}) if sols and n <= 8 and len(sols) < 400 else ('-' if not sols else '?')
    store[n] = sols
    print(n, '|', tot, '|', conf, '|', two, '|', dict(avals), '|', classes)

print("\nhistogram: number of distinct eigenvalues over ALL gauge-fixed relations, n=6 (control)")
h = Counter()
for S in all_normalized(6):
    for s in S:
        ev = np.round(np.linalg.eigvalsh(s.astype(float)), 6); h[len(set(ev))] += 1
print(dict(sorted(h.items())))

S = store[6][0].astype(int)
print("\nthe n=6 solution (one representative):\n", S)
ev, U = np.linalg.eigh(S.astype(float))
print("eigenvalues:", np.round(ev, 6), " (+-sqrt5 =", round(math.sqrt(5), 6), ")")
P = (np.eye(6) + S/math.sqrt(5))/2
print("rank of projector P+ :", int(round(np.trace(P))), " -> six vectors P e_i in a 3-dim space")
G = P/ P[0, 0]
print("normalized Gram (cosines):\n", np.round(G, 6))
print("all |cos| off-diagonal equal 1/sqrt5 =", np.allclose(np.abs(G[~np.eye(6, dtype=bool)]), 1/math.sqrt(5)))
# automorphisms: permutations p such that S[p,p] is switching-equivalent to S
def sw_equiv(T, S):
    d = (T[0]*S[0]).copy(); d[0] = 1     # need D with D T D = S ; first rows fix D up to sign
    # general: find d by propagating from vertex 0
    d = np.ones(len(S), dtype=int)
    for j in range(1, len(S)): d[j] = T[0, j]*S[0, j]
    return np.array_equal(T * d[:, None] * d[None, :], S)
aut = [p for p in itertools.permutations(range(6)) if sw_equiv(S[np.ix_(p, p)], S)]
anti = [p for p in itertools.permutations(range(6)) if sw_equiv(-S[np.ix_(p, p)], S)]
print("permutations preserving the switching class:", len(aut), " (of 720);  mapping S -> -S (self-complementary):", len(anti))
orders = Counter()
for p in aut:
    q = tuple(range(6)); k = 0
    while True:
        q = tuple(p[i] for i in q); k += 1
        if q == tuple(range(6)): break
    orders[k] += 1
print("element orders in that group:", dict(sorted(orders.items())), " (A5: 1:1, 2:15, 3:20, 5:24)")
# 2-transitivity
print("2-transitive on the 6 points:", len({(p[0], p[1]) for p in aut}) == 30)

print("\nNON-exhaustive constructive controls for larger n (Paley construction, q prime = 1 mod 4):")
for q in (5, 13, 17, 29):
    sq = {(x*x) % q for x in range(1, q)}
    n = q+1; S = np.zeros((n, n), dtype=int); S[0, 1:] = 1; S[1:, 0] = 1
    for i in range(q):
        for j in range(q):
            if i != j: S[i+1, j+1] = 1 if (i-j) % q in sq else -1
    ok = np.array_equal(S @ S, (n-1)*np.eye(n, dtype=int))
    print("  n=%d: S^2=(n-1)I %s -> eigenspace dims %d + %d" % (n, ok, n//2, n//2))
print("  => the rule alone gives dimension n/2 for every admissible n; only 'smallest non-trivial' gives 3.")
