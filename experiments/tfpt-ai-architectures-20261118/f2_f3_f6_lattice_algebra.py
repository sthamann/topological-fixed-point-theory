"""F2, F3, F6 (algebraic parts) - exact integer / rational computations on root systems.

F2  Does one spinor (doublet) current complete D8 to E8?  Same question for D5 + A3 with one
    Z4 glue vector.  Closure under root addition, generation paths, Frenkel-Kac cocycle, Jacobi.
F3  Why four marks?  Scan the mark number m in the class D_{m+1} + A_{m-1} (function space of
    dimension m+1, difference lattice of rank m-1): for which m can the two be glued to an even
    unimodular lattice?  Plus: words in two involutions with (AB)^2 = +1 or -1 up to depth 6.
F6  Chirality count: how do the E8 roots split under D5 + A3, and do mirror families appear?
    Standard Model charges, anomaly sums and the coupling normalisation from the half-spinor.
"""
import itertools, json, random
from fractions import Fraction as Fr
import numpy as np

res = {}

# ------------------------------------------------------------------ helpers
def d_roots(n):
    R = []
    for i, j in itertools.combinations(range(n), 2):
        for a in (1, -1):
            for b in (1, -1):
                v = [Fr(0)] * n; v[i], v[j] = Fr(a), Fr(b); R.append(tuple(v))
    return R

def spinors(n, parity):
    return [tuple(Fr(s, 2) for s in sg) for sg in itertools.product((1, -1), repeat=n)
            if sum(x < 0 for x in sg) % 2 == parity]

def ip(a, b):
    return sum(x * y for x, y in zip(a, b))

def close(gens):
    """Smallest set of norm-2 vectors containing gens and closed under addition and negation.
    Returns the set and, for each element, the step at which it first appeared."""
    S = {}
    for g in gens:
        S[g] = 0; S[tuple(-x for x in g)] = 0
    step, frontier = 0, list(S)
    while frontier:
        step += 1
        new = []
        cur = list(S)
        for a in frontier:
            for b in cur:
                if ip(a, b) == -1:
                    c = tuple(x + y for x, y in zip(a, b))
                    if c not in S:
                        S[c] = step; new.append(c)
        frontier = new
    return S

# ------------------------------------------------------------------ F2
D8 = d_roots(8)
cand = spinors(8, 0)[0]
full = close(D8 + [cand])
f2 = {"D8_roots": len(D8), "norm_of_candidate": str(ip(cand, cand)),
      "closure_with_one_spinor_current": len(full), "generation_steps_needed": max(full.values()),
      "currents_incl_8_Cartan": len(full) + 8,
      "control_no_candidate": len(close(D8)),
      "control_conjugate_spinor": len(close(D8 + [spinors(8, 1)[0]])),
      "control_both_spinors_norm2_violation": None,
      "control_vector_weight_norm": str(ip((Fr(1),) + (Fr(0),) * 7, (Fr(1),) + (Fr(0),) * 7))}
# S+ and S- together: inner products become half-integral -> not an integral lattice
sp, sm = spinors(8, 0)[0], spinors(8, 1)[0]
f2["control_both_spinors_norm2_violation"] = f"inner product {ip(sp, sm)} is not an integer"
# D5 + A3 with one glue vector (half-spinor of D5 times fundamental weight of A3)
D5 = [r + (Fr(0),) * 4 for r in d_roots(5)]
A3 = [(Fr(0),) * 5 + tuple(Fr(int(k == i) - int(k == j)) for k in range(4))
      for i in range(4) for j in range(4) if i != j]
glue = spinors(5, 0)[0] + (Fr(3, 4), Fr(-1, 4), Fr(-1, 4), Fr(-1, 4))
g = close(D5 + A3 + [glue])
f2["D5_A3"] = dict(roots=len(D5) + len(A3), glue_norm=str(ip(glue, glue)), closure_with_one_glue_vector=len(g),
                   steps=max(g.values()), dimension=len(g) + 8)
# Lie algebra level: Frenkel-Kac cocycle on the closed root set, Jacobi identity on random triples
roots = [np.array([float(x) for x in r]) for r in full]
key = {tuple(np.round(r * 2).astype(int)): i for i, r in enumerate(roots)}
simple = [np.array(v, float) for v in (
    [1, -1, 0, 0, 0, 0, 0, 0], [0, 1, -1, 0, 0, 0, 0, 0], [0, 0, 1, -1, 0, 0, 0, 0], [0, 0, 0, 1, -1, 0, 0, 0],
    [0, 0, 0, 0, 1, -1, 0, 0], [0, 0, 0, 0, 0, 1, -1, 0], [0, 0, 0, 0, 0, 1, 1, 0], [-.5] * 8)]
Bm = np.array(simple)
coords = [np.round(np.linalg.solve(Bm.T, r)).astype(int) for r in roots]
assert all(np.allclose(c @ Bm, r) for c, r in zip(coords, roots)), "simple roots do not span"
Gm = np.round(Bm @ Bm.T).astype(int)
def eps(a, b):                       # bilinear cocycle: eps(a_i,a_j) = (-1)^{G_ij} for i>j, -1 for i=j
    e = 0
    for i in range(8):
        for j in range(8):
            if a[i] % 2 and b[j] % 2 and (i > j and Gm[i, j] % 2 or i == j):
                e ^= 1
    return -1 if e else 1
def bracket(i, j):
    """[e_i, e_j] for root vectors: returns (coefficient, index) or None; Cartan part handled separately."""
    s = roots[i] + roots[j]
    k = key.get(tuple(np.round(s * 2).astype(int)))
    return None if k is None or abs(s @ s - 2) > 1e-9 else (eps(coords[i], coords[j]), k)
random.seed(0)
viol = tested = 0
for _ in range(20000):
    i, j, k = random.sample(range(len(roots)), 3)
    # Jacobi on triples whose pairwise sums are never zero (pure root-vector part)
    if any(np.allclose(roots[a] + roots[b], 0) for a, b in ((i, j), (j, k), (i, k))):
        continue
    if np.allclose(roots[i] + roots[j] + roots[k], 0):
        continue
    tot = {}
    for a, b, c in ((i, j, k), (j, k, i), (k, i, j)):
        x = bracket(b, c)
        if x:
            y = bracket(a, x[1])
            if y:
                tot[y[1]] = tot.get(y[1], 0) + x[0] * y[0]
    tested += 1
    viol += any(v != 0 for v in tot.values())
f2["jacobi_root_triples_tested"] = tested
f2["jacobi_violations"] = viol
res["F2"] = f2

# ------------------------------------------------------------------ F3
def disc_group(kind, n):
    """discriminant group as list of (order) cyclic factors and the quadratic form values mod 2."""
    if kind == "A":       # A_n: Z_{n+1}, q(k) = k^2 n/(n+1)
        return [n + 1], {k: Fr(k * k * n, n + 1) % 2 for k in range(n + 1)}
    if n % 2:             # D_n odd: Z4, q = n/4, 1, n/4
        return [4], {0: Fr(0), 1: Fr(n, 4) % 2, 2: Fr(1), 3: Fr(n, 4) % 2}
    return [2, 2], {0: Fr(0), 1: Fr(n, 4) % 2, 2: Fr(1), 3: Fr(n, 4) % 2}
scan = {}
for m in range(2, 41):
    gA, qA = disc_group("A", m - 1)
    gD, qD = disc_group("D", m + 1)
    same_group = sorted(gA) == sorted(gD)
    rank = 2 * m
    ok = False
    if same_group and gA == [4]:
        # glue with an isomorphism k -> +-k : need qA(k) + qD(+-k) = 0 mod 2 for all k
        ok = any(all((qA[k] + qD[(s * k) % 4]) % 2 == 0 for k in range(4)) for s in (1, 3))
    scan[m] = dict(rank=rank, disc_A=gA, disc_D=gD, same_group=same_group, rank_multiple_of_8=rank % 8 == 0,
                   even_unimodular_gluing=ok)
res["F3"] = dict(
    mark_numbers_with_even_unimodular_gluing=[m for m, v in scan.items() if v["even_unimodular_gluing"]],
    mark_numbers_with_matching_discriminant_group=[m for m, v in scan.items() if v["same_group"]],
    binary_depths={d: scan[2 ** d] for d in (1, 2, 3, 4, 5)},
    scope="class D_{m+1} + A_{m-1} only; the weighting of transitions is not derived here")

# words in two involutions A, B with (AB)^2 = sign; normal forms and number of distinct elements
def words(sign, depth):
    # elements as (phase, a, b, order) normal form  phase * A^a B^b  using BA = sign * AB
    seen = {(1, 0, 0)}
    frontier = [(1, 0, 0)]
    per_depth = []
    for _ in range(depth):
        new = []
        for ph, a, b in frontier:
            for gen in "AB":
                if gen == "A":              # (A^a B^b) A = sign^b A^{a+1} B^b
                    e = (ph * (sign ** b), (a + 1) % 2, b)
                else:
                    e = (ph, a, (b + 1) % 2)
                if e not in seen:
                    seen.add(e); new.append(e)
        frontier = new
        per_depth.append(len(seen))
    return per_depth
res["F3"]["two_involution_words"] = {
    "(AB)^2=+1": dict(elements_by_depth=words(1, 6), irreps="four 1-dimensional characters (Klein four group)"),
    "(AB)^2=-1": dict(elements_by_depth=words(-1, 6), irreps="one 2-dimensional representation (Pauli pair), Re w(ABAB) = -1"),
    "note": "normal form phase*A^a*B^b is unique; both rewrite systems are confluent. The sign of ABAB is "
            "the single datum that separates four commuting marks from one anticommuting pair."}

# ------------------------------------------------------------------ F6
# split the 240 roots built from D5 + A3 + glue by (D5 class, A3 class)
def d5_class(v):
    w = v[:5]
    if all(x.denominator == 1 for x in w):
        return "adj/0" if ip(w, w) in (0, 2) else "vector"
    return "16" if sum(x < 0 for x in w) % 2 == 0 else "16bar"
cnt = {}
for r in g:
    k = (d5_class(r), str(ip(r[5:], r[5:])))
    cnt[k] = cnt.get(k, 0) + 1
split = {f"D5 {a}, A3 norm {b}": n for (a, b), n in sorted(cnt.items())}
n16 = sum(n for (a, b), n in cnt.items() if a == "16")
n16b = sum(n for (a, b), n in cnt.items() if a == "16bar")
# one family: the 16 half-spinor weights, 3+2 marking X = (-2,-2,-2,3,3)
fam = spinors(5, 0)
Y, T3 = [], []
for w in fam:
    s = [int(2 * x) for x in w]                       # +-1
    Y.append(Fr(-1, 3) * sum(s[:3]) / 2 + Fr(1, 2) * sum(s[3:]) / 2)   # hypercharge from the 3+2 marking
    T3.append(Fr(s[3] - s[4], 4))
Q = [y + t for y, t in zip(Y, T3)]
col = [tuple(int(2 * x) for x in w[:3]) for w in fam]
triplet = [c.count(1) in (1, 2) for c in col]
anom = dict(sum_Y=str(sum(Y)), sum_Y3=str(sum(y ** 3 for y in Y)),
            sum_Y_on_weak_doublets=str(sum(y for y, t in zip(Y, T3) if t != 0)),
            sum_Y_on_colour_triplets=str(sum(y for y, t in zip(Y, triplet) if t)),
            sum_Q=str(sum(Q)))
res["F6"] = dict(root_split=split, family_copies_16=n16 // 16, mirror_copies_16bar=n16b // 16,
                 net_chirality_in_adjoint=(n16 - n16b) // 16,
                 charges_in_one_16=sorted({str(q) for q in Q}), anomaly_sums=anom,
                 sin2_thetaW_at_common_normalisation=str(sum(t * t for t in T3) / sum(q * q for q in Q)),
                 reading="E8 -> D5 + A3 contains four 16 and four 16bar with equal multiplicity: the algebra alone "
                         "is vector-like.  Three chiral families need the selection (which of the four weights, "
                         "and removal of the mirrors) from the actual weighted source; that part is not computed.")
json.dump(res, open("out/f2_f3_f6_lattice_algebra.json", "w"), indent=1, default=str)
print(json.dumps(res, indent=1, default=str))
