"""Follow-up to computation 1: (i) all chiral irreducible ledgers of SU(3)xSU(2)xU(1) with <=16 states;
(ii) for every rigid ledger: does every multiplet carry U(1) charge, is the U(1) charge multiset itself chiral,
and what is the order of the subgroup of the centre Z_N x Z_M x U(1) acting trivially (global form)."""
import itertools
from math import gcd
from functools import reduce
from fractions import Fraction
import comp1_anomaly_ledgers as c

def kernel_order(q, mult, N, M):
    D = N * M
    cnt = 0
    for a in range(N):
        for b in range(M):
            for th in range(D):   # theta = th/D ; primitive charges => theta in (1/D)Z suffices
                if all((Fraction(a * t[0], N) + Fraction(b * t[1], M) + Fraction(th * v, D)).denominator == 1 for t, v in zip(mult, q)):
                    cnt += 1
    return cnt

def u1_chiral(q, mult, N, M):
    ch = []
    for t, v in zip(mult, q): ch += [v] * c.tdim(t, N, M)
    from collections import Counter
    cn = Counter(ch)
    return cn[0] == 0 and all(cn[v] != cn[-v] or True for v in cn) and all(cn[-v] == 0 for v in cn)

print('--- (i) SU(3)xSU(2)xU(1): every chiral irreducible ledger with <=16 states, |q|<=12, singlets<=6 ---')
N, M = 3, 2; ts = c.types_for(N, M)
for S in range(1, 17):
    for cv in c.count_vectors_states(ts, N, M, S, 6):
        if not c.type_ok(cv, ts, N, M): continue
        res, mult = c.solve(N, M, cv, ts, 12)
        res = [q for q in res if c.irreducible(q, mult, N, M)]
        if res:
            rig = [q for q in res if c.local_dim(q, mult, N, M) == 1]
            print('states=%d multiplets=%d: %d solutions in box, %d rigid' % (S, len(mult), len(res), len(rig)))
            for q in (rig + [q for q in res if q not in rig])[:4]:
                print('     %s  rigid=%s allcharged=%s' % (c.fmt(q, mult, N, M), c.local_dim(q, mult, N, M) == 1, all(v != 0 for v in q)))

print('\n--- (ii) all rigid ledgers (k<=5, N<=7, |q|<=14) with extra flags ---')
rows = []
for N in range(1, 8):
    for M in range(1, N + 1):
        ts = c.types_for(N, M)
        for k in range(1, 6):
            for cv in c.count_vectors_k(ts, k):
                if not c.type_ok(cv, ts, N, M): continue
                res, mult = c.solve(N, M, cv, ts, 14)
                for q in res:
                    if c.irreducible(q, mult, N, M) and c.local_dim(q, mult, N, M) == 1:
                        rows.append((sum(c.tdim(t, N, M) for t in mult), N, M, c.fmt(q, mult, N, M), all(v != 0 for v in q), u1_chiral(q, mult, N, M), kernel_order(q, mult, N, M)))
rows.sort()
for r in rows:
    print('states=%2d SU(%d)xSU(%d)xU(1): %-55s every multiplet charged=%-5s U(1)-set chiral=%-5s trivially-acting centre order=%d' % r)
print('\nminimal rigid ledger overall:', rows[0][:4])
full = [r for r in rows if r[4]]
print('minimal rigid ledger with every multiplet U(1)-charged:', full[0][:4], ' centre kernel order', full[0][6])
