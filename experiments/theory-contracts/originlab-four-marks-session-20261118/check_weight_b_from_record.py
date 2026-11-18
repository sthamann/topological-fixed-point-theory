#!/usr/bin/env python3
"""The weight b of the side pairings against the record (18 November 2026).

Recorded, unchanged: ledger row SEAM.CODE.TYPEII.01 / verification/v799_seam_code_typeii.py:
    eight sides = four marks x two sides, mark k owning the coordinates {2k, 2k+1};
    A1 mutual locality  <=> the code is self-orthogonal,
    A2 integer conformal spin (trivial T phases i^wt) <=> doubly even,
    A3 holomorphy / index one <=> self-dual;
    the survivors of all three are the 30 copies of the Hamming code (enumerator 1 + 14 x^4 + x^8).

Sections
    B1  all self-dual codes of length eight: 105 + 30
    B2  their theta functions: theta3^8 (one spin structure) and (theta2^8 + theta3^8 + theta4^8)/2 (E8)
    B3  the modular maps on the three determinants; which axiom fixes b
    B4  eleventh pasted note: overlap on the sixteen positions; positivity of the mark matrix
"""
import itertools
import sys

import mpmath as mp
import numpy as np

mp.mp.dps = 30
PI = mp.pi
FAIL = []


def check(name, ok, info=None):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info is not None else ""))
    if not ok:
        FAIL.append(name)


# ---------------------------------------------------------------- B1
print("B1  self-dual codes of length eight")
wt = lambda w: bin(w).count("1")
even = [w for w in range(1, 256) if wt(w) % 2 == 0 and w != 255]
orth = lambda a, b: wt(a & b) % 2 == 0
codes = set()
for v2 in even:
    for v3 in even:
        if v3 <= v2 or not orth(v2, v3):
            continue
        for v4 in even:
            if v4 <= v3 or not orth(v2, v4) or not orth(v3, v4):
                continue
            span = {0}
            for g in (255, v2, v3, v4):
                span |= {s ^ g for s in span}
            if len(span) == 16:
                codes.add(frozenset(span))
check("there are 135 self-dual codes (A1 and A3 together)", len(codes) == 135, len(codes))


def enum(code):
    out = [0]*9
    for w in code:
        out[wt(w)] += 1
    return tuple(out)


kinds = {}
for c in codes:
    kinds.setdefault(enum(c), []).append(c)
pair_enum = (1, 0, 4, 0, 6, 0, 4, 0, 1)          # (1 + x^2)^4
ham_enum = (1, 0, 0, 0, 14, 0, 0, 0, 1)
check("exactly two kinds: 105 with enumerator (1 + x^2)^4 and 30 with 1 + 14 x^4 + x^8 (recorded count 30)",
      set(kinds) == {pair_enum, ham_enum} and len(kinds[pair_enum]) == 105 and len(kinds[ham_enum]) == 30)
mark_pairs = [0b11 << (2*k) for k in range(4)]
span = {0}
for g in mark_pairs:
    span |= {s ^ g for s in span}
check("the code 'one word per mark, both sides together' (recorded bookkeeping) is one of the 105",
      frozenset(span) in kinds[pair_enum])
check("it has words of weight two: T phase i^2 = -1, half-integer spin; A2 fails for all 105 and holds for the 30",
      all(any(wt(w) % 4 for w in c) for c in kinds[pair_enum]) and all(all(wt(w) % 4 == 0 for w in c) for c in kinds[ham_enum]))
check("105 = 7 * 5 * 3: the ways of grouping eight sides into four pairs", 7*5*3 == 105)

# ---------------------------------------------------------------- B2
print("B2  theta functions")


def thetas(tau):
    N = 80
    q = lambda e: mp.exp(1j*PI*tau*e)
    return (mp.nsum(lambda n: q((n+mp.mpf(1)/2)**2), [-N, N]), mp.nsum(lambda n: q(n*n), [-N, N]),
            mp.nsum(lambda n: (-1)**int(n)*q(n*n), [-N, N]))


tau = mp.mpc("0.17", "0.93")
t2, t3, t4 = thetas(tau)
a_, b_, _ = thetas(2*tau)[1], thetas(2*tau)[0], None
check("pair code: (a^2 + b^2)^4 = theta3(tau)^8, the lattice Z^8, sixteen free fermions with one spin structure",
      abs((a_**2+b_**2)**4-t3**8) < mp.mpf(10)**(-22))
check("Hamming code: a^8 + 14 a^4 b^4 + b^8 = (theta2^8 + theta3^8 + theta4^8)/2, the lattice E8",
      abs(a_**8+14*a_**4*b_**4+b_**8-(t2**8+t3**8+t4**8)/2) < mp.mpf(10)**(-22))
print("      so given A1 and A3 the weight b can only be 0 (pair code) or 1 (Hamming code); nothing in between.")

# ---------------------------------------------------------------- B3
print("B3  which axiom fixes b")
eta = lambda t: mp.exp(1j*PI*t/12)*mp.nprod(lambda n: 1-mp.exp(2j*PI*t*n), [1, mp.inf])
Z = lambda t: [(x/eta(t))**8 for x in thetas(t)]
z0, zT, zS = Z(tau), Z(tau+1), Z(-1/tau)
om = mp.exp(-2j*PI/3)
check("tau -> tau + 1: the determinants of theta3 and theta4 are exchanged, theta2 stays, common phase exp(-2 pi i/3)",
      abs(zT[1]-om*z0[2]) < mp.mpf(10)**(-20) and abs(zT[2]-om*z0[1]) < mp.mpf(10)**(-20) and abs(zT[0]-om*z0[0]) < mp.mpf(10)**(-20))
check("tau -> -1/tau: theta2 and theta4 exchanged, theta3 fixed, no phase",
      abs(zS[0]-z0[2]) < mp.mpf(10)**(-20) and abs(zS[2]-z0[0]) < mp.mpf(10)**(-20) and abs(zS[1]-z0[1]) < mp.mpf(10)**(-20))
fam = lambda zz, bb: zz[1]+bb*(zz[0]+zz[2])
dev_S = max(abs(fam(zS, bb)-fam(z0, bb)) for bb in (0, mp.mpf("0.4"), 1))
dev_T = {str(bb): abs(fam(zT, bb)-om*fam(z0, bb)) for bb in (0, mp.mpf("0.4"), 1)}
check("the family theta3 + b (theta2 + theta4) is invariant under tau -> -1/tau for every b (the quarter turn)",
      dev_S < mp.mpf(10)**(-20))
check("it is invariant under tau -> tau + 1 up to the common phase only for b = 1: this is A2, integer spin",
      dev_T["1"] < mp.mpf(10)**(-20) and dev_T["0"] > mp.mpf("0.1") and dev_T["0.4"] > mp.mpf("0.1"),
      {k: mp.nstr(v, 3) for k, v in dev_T.items()})
print("      b = 0: the free seam, fermionic words allowed.  b = 1: only words of integer spin.")

# ---------------------------------------------------------------- B4
print("B4  eleventh pasted note")
N = 16
ks = [(n+0.5)*2*np.pi/N for n in range(-N//2, 0)]          # occupied antiperiodic momenta
sites = np.arange(N)
g = np.where((sites//4) % 2 == 0, 1.0+0j, 1j)              # four sites per quadrant, marks on bond midpoints
M = np.array([[np.sum(g*np.exp(-1j*(k-kp)*sites))/N for kp in ks] for k in ks])
ov = abs(np.linalg.det(M))
check("quarter-phase probe on the sixteen positions of v519: overlap 0.461663 (the note's correction of 0.338394)",
      abs(ov-0.461663) < 5e-6, round(float(ov), 6))
C = np.array([[0, 1, 2, 1], [1, 0, 1, 2], [2, 1, 0, 1], [1, 2, 1, 0]], float)
ev = np.linalg.eigvalsh(-C)
check("the recorded mark matrix with zero diagonal has eigenvalues of both signs: it is no positive Gram form",
      ev.min() < 0 < ev.max(), np.round(ev, 6))

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
