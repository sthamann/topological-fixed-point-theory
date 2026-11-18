#!/usr/bin/env python3
"""Session scenarios of 17 November 2026: the fourth family as return channel.

Exact checks on the D8 half-spinor (128 matter roots of E8) with Cl(16) in
Jordan-Wigner form. No recorded value is used as input. Prints JSON.
"""
import json
from fractions import Fraction as Fr
from itertools import product
import numpy as np

N = 8
I2, Z = np.eye(2), np.diag([1.0, -1.0])
LOW = np.array([[0.0, 1.0], [0.0, 0.0]])          # lowers the occupation
def kron(ms):
    out = np.eye(1)
    for m in ms:
        out = np.kron(out, m)
    return out
A = [kron([Z]*k+[LOW]+[I2]*(N-k-1)) for k in range(N)]     # annihilators
G1 = [a+a.T for a in A]                                    # real directions
G2 = [1j*(a-a.T) for a in A]
STATES = list(product((1, 0), repeat=N))                   # bit 1 = empty in this basis order
weight = lambda s: tuple(1 if b == 0 else -1 for b in s)   # doubled weight
W = [weight(s) for s in STATES]
EVEN = [n for n, w in enumerate(W) if w.count(-1) % 2 == 0]
assert len(EVEN) == 128

def family(w):
    f = w[5:]
    return tuple(-x for x in f) if w[:5].count(-1) % 2 else f
FAMS = sorted({family(W[n]) for n in EVEN})
assert len(FAMS) == 4
YH = (Fr(-1, 3), Fr(-1, 3), Fr(-1, 3), Fr(1, 2), Fr(1, 2), 0, 0, 0)

out = {}
# 1. hypercharge census of the 60 middle currents (sigma e_i + tau e_j, i<5<=j)
weak = [(i, s, j, t) for i in (3, 4) for s in (1, -1) for j in (5, 6, 7) for t in (1, -1)]
colour = [(i, s, j, t) for i in (0, 1, 2) for s in (1, -1) for j in (5, 6, 7) for t in (1, -1)]
assert {abs(YH[i]) for i, *_ in weak} == {Fr(1, 2)} and {abs(YH[i]) for i, *_ in colour} == {Fr(1, 3)}
out["middle_currents"] = {"weak_Y_half_centre_neutral": len(weak), "colour_Y_third_centre_charged": len(colour)}

# 2. returns through the omitted family: number operator sum_gamma P X_-gamma Q X_gamma P
def steps(w, cur):
    i, s, j, t = cur
    if w[i] == -s and w[j] == -t:
        v = list(w); v[i], v[j] = s, t
        return tuple(v)
census = {}
for omit in FAMS:
    for name, S in (("all", weak+colour), ("weak", weak), ("colour", colour)):
        vals = set()
        for n in EVEN:
            w = W[n]
            if family(w) == omit:
                continue
            through = sum(1 for c in S if (v := steps(w, c)) and family(v) == omit)
            total = sum(1 for c in S if steps(w, c))
            vals.add((through, total))
        census.setdefault(name, set()).update(vals)
out["returns_per_state_(through_fourth,total)"] = {k: sorted(v) for k, v in census.items()}

# 3. real middle current H = i Gamma_v Gamma_omega on the half-spinor
idx = np.array(EVEN)
def restrict(M):
    return M[np.ix_(idx, idx)]
rng = np.random.default_rng(7)
def spectrum(M):
    vals = np.round(np.linalg.eigvalsh(M), 9)+0.0
    return {str(float(v)): int(c) for v, c in zip(*np.unique(vals, return_counts=True))}
fam_dirs = [G1[5], G1[6], G1[7], G2[5], G2[6], G2[7]]
report = []
for trial in range(6):
    c = rng.normal(size=6); c /= np.linalg.norm(c)
    if trial < 3:
        c = np.eye(6)[trial]
    v = rng.normal(size=4); v /= np.linalg.norm(v)
    Gv = v[0]*G1[3]+v[1]*G2[3]+v[2]*G1[4]+v[3]*G2[4]          # the four weak (Higgs) directions
    Gw = sum(x*g for x, g in zip(c, fam_dirs))
    H = restrict(1j*Gv@Gw)
    assert np.allclose(H, H.conj().T) and np.allclose(H@H, np.eye(128))
    for omit in FAMS:
        keep = np.array([family(W[n]) != omit for n in EVEN])
        P, Q = np.diag(keep.astype(float)), np.diag((~keep).astype(float))
        A0 = (P@H@P)[np.ix_(keep, keep)]
        Gm = (P@H@Q@H@P)[np.ix_(keep, keep)]
        CQ = Q@H@Q
        assert np.allclose(CQ, 0)                               # no direct term in the fourth family
        assert np.allclose(A0@A0+Gm, np.eye(96))                # A0^2 + G = 1
        report.append((spectrum(A0), spectrum(Gm)))
assert all(r == report[0] for r in report)
out["real_weak_current_on_three_families"] = {"A0=PHP": report[0][0], "G=PHQHP": report[0][1],
                                              "identity": "A0^2 + G = 1, QHQ = 0"}

# 4. two real currents with one carrier direction and orthogonal family directions
K1, K2 = restrict(G1[4]@G1[5]), restrict(G1[4]@G1[6])
out["ABAB"] = "-1" if np.allclose(K1@K2@K1@K2, -np.eye(128)) else "other"
out["anticommute"] = bool(np.allclose(K1@K2+K2@K1, 0))

# 5. the three coordinate family directions pair the four families (V4 torsor)
pairs = {}
for j in (5, 6, 7):
    M = restrict(G1[4]@G1[j])
    link = set()
    for a, n in enumerate(EVEN):
        for b, m in enumerate(EVEN):
            if abs(M[b, a]) > 1e-12:
                link.add(frozenset((FAMS.index(family(W[n])), FAMS.index(family(W[m])))))
    pairs[j] = sorted(sorted(p) for p in link)
out["family_pairings_by_direction"] = pairs
print(json.dumps(out, indent=1))
