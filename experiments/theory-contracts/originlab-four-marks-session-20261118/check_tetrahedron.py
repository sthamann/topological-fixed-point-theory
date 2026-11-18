#!/usr/bin/env python3
"""Session patterns of 17 November 2026: the family tetrahedron.

1. Casimir split of so(16) on the 128 half-spinor, resolved by family:
   stay / hop rates of one normalized step.
2. Which D_p x D_q in D_8 (and in other ranks) has stay = 1/2.
3. Compression of that step onto one, two, three families and its sixth power.
4. Wick pairings of a weight-one four-point function at the marks 1, i, -1, -i.
Exact rational arithmetic. No recorded value enters. Prints JSON.
"""
import json
from fractions import Fraction as Fr
from itertools import combinations, product
import sympy as sp

def spinor(n):
    return [w for w in product((1, -1), repeat=n) if w.count(-1) % 2 == 0]   # doubled weights

def roots(n):
    out = []
    for i, j in combinations(range(n), 2):
        for a, b in product((2, -2), repeat=2):
            v = [0]*n; v[i], v[j] = a, b
            out.append(tuple(v))
    return out

def split(p, q):
    """Casimir of so(2(p+q)) on its half-spinor, by class of generator."""
    n = p+q
    w = spinor(n)[0]
    own = Fr(0); mixed = Fr(0)
    for r in roots(n):
        if all(w[k]*r[k] <= 0 for k in range(n)) and sum(1 for k in range(n) if r[k]) == 2 \
           and all(w[k] == -r[k]//2 for k in range(n) if r[k]):
            inside = sum(1 for k in range(n) if r[k] and k < p)
            if inside == 1: mixed += 1
            else: own += 1
    own += Fr(n, 4)                      # Cartan part: |weight|^2 with roots of norm 2
    return own, mixed

out = {}
own, mixed = split(5, 3)
out["casimir_128"] = {"own_D5_plus_A3": str(own), "mixed_10_6": str(mixed), "total": str(own+mixed),
                      "stay": str(own/(own+mixed))}
assert own == 15 and mixed == 15

# family-resolved step on the 128: stay or move to each other family
def family(w):
    f = w[5:]
    return tuple(-x for x in f) if w[:5].count(-1) % 2 else f
S = spinor(8); FAM = sorted({family(w) for w in S})
kernel = {}
for w in S:
    row = {f: Fr(0) for f in FAM}
    row[family(w)] += Fr(8, 4)
    for r in roots(8):
        if all(w[k] == -r[k]//2 for k in range(8) if r[k]):
            v = tuple(w[k]+r[k] for k in range(8))
            row[family(v)] += 1
    key = (family(w), tuple(str(row[f]/30) for f in FAM))
    kernel[key] = kernel.get(key, 0)+1
rows = {FAM.index(f): r for (f, r) in kernel}
assert len(kernel) == 4
out["family_step"] = rows

M = sp.Matrix(4, 4, lambda a, b: sp.Rational(1, 2) if a == b else sp.Rational(1, 6))
for k in (4, 3, 2, 1):
    sub = M[:k, :k]
    out[f"keep_{k}_families"] = {"step": {str(e): m for e, m in sub.eigenvals().items()},
                                  "six_steps": {str(e**6): m for e, m in sub.eigenvals().items()}}

# uniqueness of stay = 1/2 among D_p x D_q
hits = []
for n in range(2, 41):
    for q in range(1, n//2+1):
        p = n-q
        o = Fr(p*(2*p-1)+q*(2*q-1), 4); m = Fr(p*q)
        if o == m: hits.append((p, q))
out["stay_half_solutions_up_to_rank_40"] = hits

# Wick pairings at the square marks
marks = [sp.Integer(1), sp.I, sp.Integer(-1), -sp.I]
K = lambda a, b: sp.simplify(a*b/(a-b)**2)
pairings = {"diagonal": [(0, 2), (1, 3)], "side_a": [(0, 1), (2, 3)], "side_b": [(1, 2), (3, 0)]}
wt = {k: sp.simplify(K(marks[a], marks[b])*K(marks[c], marks[d])) for k, ((a, b), (c, d)) in pairings.items()}
tot = sum(wt.values())
out["wick_pairings"] = {k: [str(v), str(sp.nsimplify(v/tot))] for k, v in wt.items()}
print(json.dumps(out, indent=1))
