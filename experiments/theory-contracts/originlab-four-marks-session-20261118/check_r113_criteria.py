#!/usr/bin/env python3
"""Session test of 17 November 2026: the algebraic R113 criteria on the E8
matter space with nested family projectors.

This is NOT the native R113 call (that call stays blocked at Q_0[a,b]). It
evaluates the recorded finite criteria
    A0 = PHP, CQ = QHQ, B = QHP, G = B^T B, K = B^T CQ B, L = B^T CQ^2 B,
    D = L - K G^+ K,   quarter class: D = 0, K = lam G, A0 = lam, G = lam^2
    golden class:  D = 0, K = lam G, A0 = 2 lam          (on the support of G)
for H = c J^2 + J with J = X_gamma + X_-gamma, every middle current gamma and
every set of retained families. Integer arithmetic. Prints JSON.
"""
import json
from fractions import Fraction as Fr
from itertools import combinations, product
import numpy as np

matter = [s for s in product((1, -1), repeat=8) if s.count(-1) % 2 == 0]
index = {w: n for n, w in enumerate(matter)}
def family(w):
    f = w[5:]
    return tuple(-x for x in f) if w[:5].count(-1) % 2 else f
fams = sorted({family(w) for w in matter})
def sign(w, i, j):                      # Jordan-Wigner sign of a_i^# a_j^#, i < j
    return -1 if sum(1 for k in range(i+1, j) if w[k] == 1) % 2 else 1
def current(i, s, j, t):
    X = np.zeros((128, 128), dtype=object); X[:] = 0
    for w in matter:
        if w[i] == -s and w[j] == -t:
            v = list(w); v[i], v[j] = s, t
            X[index[tuple(v)], index[w]] = sign(w, i, j)
    return X
def pinv_diag_support(G):
    return [n for n in range(G.shape[0]) if any(G[n, m] != 0 for m in range(G.shape[0]))]

def classify(H, keep):
    p = [n for n, w in enumerate(matter) if family(w) in keep]
    q = [n for n in range(128) if n not in p]
    A0, CQ, B = H[np.ix_(p, p)], H[np.ix_(q, q)], H[np.ix_(q, p)]
    G, K, L = B.T@B, B.T@CQ@B, B.T@CQ@CQ@B
    sup = pinv_diag_support(G)
    if not sup:
        return "closed"                                   # no return through the complement
    Gs, Ks, Ls, As = (M[np.ix_(sup, sup)] for M in (G, K, L, A0))
    g = Gs[0, 0]
    if not (Gs == g*np.eye(len(sup), dtype=object)).all():
        return "G_not_scalar"
    lam = Fr(int(Ks[0, 0]), int(g))
    if not all(Fr(int(Ks[a, b])) == lam*int(Gs[a, b]) for a in range(len(sup)) for b in range(len(sup))):
        return "K_not_lam_G"
    D = np.array([[Fr(int(Ls[a, b]))-sum(Fr(int(Ks[a, c]))*Fr(int(Ks[c, b]))/int(g) for c in range(len(sup)))
                   for b in range(len(sup))] for a in range(len(sup))], dtype=object)
    if any(x != 0 for x in D.flat):
        return "D_nonzero"
    diag = all(Fr(int(As[a, b])) == (lam if a == b else 0) for a in range(len(sup)) for b in range(len(sup)))
    gold = all(Fr(int(As[a, b])) == (2*lam if a == b else 0) for a in range(len(sup)) for b in range(len(sup)))
    if diag and g == lam*lam: return f"quarter(lam={lam})"
    if gold and lam != 0: return f"golden(lam={lam})"
    return f"D0_KlamG_only(lam={lam},G={g},A0_scalar={diag})"

out = {}
for c in (0, 1, 2):
    tally = {}
    for i in range(5):
        for j in (5, 6, 7):
            for s, t in ((1, 1), (1, -1)):
                X = current(i, s, j, t); J = X+X.T
                assert ((J@J@J) == J).all()
                H = c*(J@J)+J
                for k in (1, 2, 3):
                    for keep in combinations(fams, k):
                        key = (k, classify(H, keep))
                        tally[key] = tally.get(key, 0)+1
    out[f"H = {c} J^2 + J"] = {f"keep_{k}: {v}": n for (k, v), n in sorted(tally.items())}
print(json.dumps(out, indent=1))
