#!/usr/bin/env python3
"""Session check of 17 November 2026: the recorded flavour monodromy (v117)
against the family coordinates of the D3 = A3 factor of E8.

1. An explicit unitary S carries the recorded pair (U, M0) to signed
   permutation matrices on the three family coordinates x5, x6, x7.
2. The images are det(g) g for elements g of W(D3); the four lines fixed by
   the order-three elements are the four family weights (+-1,+-1,+-1).
3. Clock ladder: middle-current returns of a matter state that stay inside
   k retained families, for every choice of the families.
Exact arithmetic. Prints JSON. Argument: repository root.
"""
import json, sys
from fractions import Fraction as Fr
from itertools import combinations, permutations, product
import sympy as sp

sys.path.insert(0, sys.argv[1]+"/verification")
from v117_monodromy_weyl_a3 import M0_EXACT as M0, U_EXACT as U

out = {}
# all 24 rotations of the cube as signed permutation matrices with det +1
rot = []
for perm in permutations(range(3)):
    for signs in product((1, -1), repeat=3):
        m = sp.zeros(3)
        for r, c in enumerate(perm):
            m[r, c] = signs[r]
        if m.det() == 1:
            rot.append(m)
assert len(rot) == 24
order = lambda m: next(k for k in range(1, 13) if m**k == sp.eye(3))
found = None
for R4 in (m for m in rot if order(m) == 4):
    for C3 in (m for m in rot if order(m) == 3):
        if (C3*R4).trace() != (M0*U).trace() or order(C3*R4) != order(sp.simplify(M0*U)):
            continue
        # intertwiner by averaging over the generated group words
        words = [(sp.eye(3), sp.eye(3))]
        seen = {tuple(sp.eye(3))}
        frontier = list(words)
        while frontier:
            a, b = frontier.pop()
            for ga, gb in ((U, R4), (M0, C3)):
                na, nb = sp.simplify(ga*a), gb*b
                if tuple(nb) not in seen:
                    seen.add(tuple(nb)); words.append((na, nb)); frontier.append((na, nb))
        if len(words) != 24:
            continue
        S = sum((b*sp.Matrix(3, 3, [1, 2, 3, 5, 7, 11, 13, 17, 19])*a.H for a, b in words), sp.zeros(3))
        S = sp.simplify(S)
        if S.rank() == 3 and sp.simplify(S*U-R4*S) == sp.zeros(3) and sp.simplify(S*M0-C3*S) == sp.zeros(3):
            found = (R4, C3, S); break
    if found: break
assert found, "no intertwiner"
R4, C3, S = found
G = sp.simplify(S.H*S)
out["intertwiner"] = {"quarter_turn_image": R4.tolist(), "order_three_image": C3.tolist(),
                      "S_dagger_S_is_scalar": bool(sp.simplify(G-G[0, 0]*sp.eye(3)) == sp.zeros(3))}
# det(g) g with g in W(D3): signed permutations with an even number of minus signs
def in_wd3(m):
    return sum(1 for x in m if x == -1) % 2 == 0
out["images_are_det_times_W(D3)"] = bool(all(in_wd3(m.det()*m) or in_wd3(-m) for m in (R4, C3))
                                         and in_wd3(-R4) and in_wd3(C3))
axes = set()
for m in rot:
    if order(m) == 3:
        v = (m-sp.eye(3)).nullspace()[0]
        v = v/abs(v[0])
        axes.add(tuple(v) if v[0] > 0 else tuple(-v))
out["order_three_axes"] = sorted(str(a) for a in axes)

# clock ladder from the lattice
R = []
for i, j in combinations(range(8), 2):
    for a, b in product((2, -2), repeat=2):
        v = [0]*8; v[i], v[j] = a, b; R.append(tuple(v))
matter = [s for s in product((1, -1), repeat=8) if s.count(-1) % 2 == 0]
RS = set(R) | set(matter)
bif = [r for r in R if sum(1 for i in range(5) if r[i]) == 1]
def family(w):
    f = w[5:]
    return tuple(-x for x in f) if w[:5].count(-1) % 2 else f
fams = sorted({family(w) for w in matter})
ladder = {}
for k in (4, 3, 2, 1):
    vals = set()
    for keep in combinations(fams, k):
        for a in matter:
            if family(a) not in keep: continue
            tot = [tuple(x+y for x, y in zip(a, g)) for g in bif]
            tot = [b for b in tot if b in RS]
            vals.add(Fr(sum(1 for b in tot if family(b) in keep), len(tot)))
    ladder[f"keep_{k}"] = sorted(str(v) for v in vals)
out["middle_current_return_fraction"] = ladder
print(json.dumps(out, indent=1, default=str))
