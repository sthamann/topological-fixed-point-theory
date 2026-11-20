"""RMFAMILY20261119 (idea 1): are carrier 5 and E8 two members of the self-dual
Reed-Muller family RM(r,2r+1)?  Exact finite arithmetic over F2."""
from common import Checks
from itertools import combinations, product, permutations
from fractions import Fraction as Fr
from math import comb
from collections import Counter
c = Checks("reed_muller")

def rm(r, m):
    pts = list(product([0, 1], repeat=m))
    rows = []
    for d in range(r + 1):
        for S in combinations(range(m), d):
            rows.append((S, [int(all(p[i] for i in S)) for p in pts]))
    return pts, rows
def rank_f2(vs):
    vs = [int("".join(map(str, v)), 2) for v in vs]; r = 0
    while vs:
        p = max(vs); vs.remove(p)
        if p == 0: continue
        r += 1; hb = p.bit_length() - 1
        vs = [v ^ p if (v >> hb) & 1 else v for v in vs]
    return r
def span(vs):
    words = {tuple([0] * len(vs[0]))}
    for v in vs: words |= {tuple(a ^ b for a, b in zip(w, v)) for w in words}
    return words

# 1. Pascal condition == self-duality of RM(2,g)
c.ok("2^g = g^2+g+2 has the single solution g=5 for g<200", [g for g in range(1, 200) if 2**g == g * g + g + 2] == [5])
c.ok("that equation is dim RM(2,g) = 2^(g-1), i.e. RM(2,g) has rate 1/2", all((2**g == g*g+g+2) == (sum(comb(g, k) for k in range(3)) == 2**(g-1)) for g in range(1, 60)))
pts5, rows5 = rm(2, 5); G5 = [v for _, v in rows5]
c.ok("RM(2,5): length 32, dimension 16", len(pts5) == 32 and rank_f2(G5) == 16)
c.ok("RM(2,5) is self-orthogonal, hence self-dual", all(sum(a * b for a, b in zip(u, v)) % 2 == 0 for u in G5 for v in G5))
W5 = Counter(sum(w) for w in span(G5))
c.ok("RM(2,5) is doubly even with minimum distance 8 (Type II [32,16,8])", all(k % 4 == 0 for k in W5) and min(k for k in W5 if k) == 8, dict(sorted(W5.items())))
pts3, rows3 = rm(1, 3); G3 = [v for _, v in rows3]
W3 = Counter(sum(w) for w in span(G3))
c.ok("RM(1,3) is the self-dual [8,4,4] Hamming code (E8 by Construction A)", rank_f2(G3) == 4 and dict(W3) == {0: 1, 4: 14, 8: 1})
c.ok("self-dual members RM(r,2r+1): m = 1,3,5,7 have dimension 2^(m-1)", all(sum(comb(2*r+1, k) for k in range(r+1)) == 2**(2*r) for r in range(6)))

# 2. does RM(2,5) carry the 3+2 split and hypercharge of the 16?
y = [Fr(-1, 3)] * 3 + [Fr(1, 2)] * 2
SM = Counter({Fr(1, 6): 6, Fr(-2, 3): 3, Fr(1, 3): 3, Fr(-1, 2): 2, Fr(1): 1, Fr(0): 1})
mono = Counter(sum((y[i] for i in S), Fr(0)) for S, _ in rows5)
even = Counter(sum((y[i] for i in S), Fr(0)) for d in (0, 2, 4) for S in combinations(range(5), d))
hodge = Counter(sum((y[i] for i in (S if len(S) != 1 else tuple(j for j in range(5) if j not in S))), Fr(0)) for S, _ in rows5)
c.ok("even exterior algebra (the recorded carrier) gives the Standard Model hypercharges", even == SM)
c.ok("RM(2,5) monomials of degree <= 2 do NOT give them directly: the 5 comes out conjugate", mono != SM, dict(mono))
c.ok("they do after the complement map on degree one (1+5+10 -> 1+10+5bar)", hodge == SM)

# 3. Hamming code, a family triple and the group of the four marks
def affine_maps():
    import itertools
    mats = [M for M in product(product([0, 1], repeat=3), repeat=3)
            if rank_f2([list(r) for r in M]) == 3]
    for M in mats:
        for t in product([0, 1], repeat=3):
            yield lambda p, M=M, t=t: tuple((sum(M[i][j] * p[j] for j in range(3)) + t[i]) % 2 for i in range(3))
maps = list(affine_maps())
c.ok("Aut of the Hamming code as AGL(3,2) has order 1344", len(maps) == 1344)
T = [(0, 0, 0), (1, 0, 0), (0, 1, 0)]
stab = [f for f in maps if {f(p) for p in T} == set(T)]
perms = {tuple(pts3.index(f(p)) for p in pts3) for f in stab}
orbs = []; seen = set()
for i in range(8):
    if i in seen: continue
    o = {p[i] for p in perms}; orbs.append(sorted(o)); seen |= o
c.ok("stabiliser of a 3-set of coordinates has order 24", len(perms) == 24)
c.ok("its orbits on the 8 coordinates have sizes 3, 1, 4", sorted(len(o) for o in orbs) == [1, 3, 4], orbs)
four = [o for o in orbs if len(o) == 4][0]
c.ok("it acts on the 4-orbit as the full S4 = W(A3), the group permuting four marks", len({tuple(p[i] for i in four) for p in perms}) == 24)
c.ok("all 56 three-sets are equivalent (one orbit), so the split 8 = 3+1+4 needs one choice only", len({frozenset(f(p) for p in T) for f in maps}) == 56)
# 4. glue norms from the same 5+3 split of the 8 coordinates
s = [Fr(1, 2)] * 8
c.ok("the all-halves E8 root splits over 5+3 coordinates with norms 5/4 and 3/4", sum(x * x for x in s[:5]) == Fr(5, 4) and sum(x * x for x in s[5:]) == Fr(3, 4))
# 5. doubling tower
d = lambda r, m: sum(comb(m, k) for k in range(r + 1)) if r >= 0 else 0
c.ok("two doubling steps: RM(2,5) = RM(2,3)+2*RM(1,3)+RM(0,3), 16 = 7+4+4+1", d(2,3) + 2*d(1,3) + d(0,3) == 16 and (d(2,3), d(1,3), d(0,3)) == (7, 4, 1))
c.ok("the doubling does NOT reproduce the recorded 1+10+5 or 3+2 content by dimension", sorted([7, 4, 4, 1]) != sorted([1, 10, 5]))
c.done({"weights_RM25": dict(W5), "orbits": orbs})
