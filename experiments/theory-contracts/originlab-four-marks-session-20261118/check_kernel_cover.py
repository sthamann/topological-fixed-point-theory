#!/usr/bin/env python3
"""Session check of 17 November 2026: the A4 kernel cover Y_q of the
four-mark sphere (recorded genus 5) as a candidate origin of the carrier.

1. Representation of A4 on the holomorphic differentials of Y_q
   (Chevalley-Weil, cross-checked by the holomorphic Lefschetz trace).
2. Operators on that five-dimensional space commuting with A4 and with the
   quarter covariance (which exchanges the two classes of three-cycles).
3. Even exterior powers of the five-space under the resulting charge.
4. The gamma-invariant unramified V4 covers of the genus-two curve
   C: eta^2 = 1 - w^6, with the Weil pairing on J(C)[2].
5. Exact theta values at tau = i.
Exact arithmetic (sympy, Fractions) except item 5 (mpmath, 40 digits).
"""
import json
from fractions import Fraction as Fr
from itertools import combinations, permutations
import sympy as sp
import mpmath as mp

out = {}
w = sp.Rational(-1, 2)+sp.sqrt(3)*sp.I/2

# ---- A4 as even permutations of four marks --------------------------------
def sign(p):
    s = 1
    for a in range(4):
        for b in range(a+1, 4):
            if p[a] > p[b]: s = -s
    return s
A4 = [p for p in permutations(range(4)) if sign(p) == 1]
mul = lambda p, q: tuple(p[q[k]] for k in range(4))
inv = lambda p: tuple(sorted(range(4), key=lambda k: p[k]))
order = lambda p: next(n for n in range(1, 5) if all((lambda x: x == k)(k) for k in [0]) and _pow(p, n) == (0, 1, 2, 3))
def _pow(p, n):
    r = (0, 1, 2, 3)
    for _ in range(n): r = mul(p, r)
    return r
fix = lambda p: sum(1 for k in range(4) if p[k] == k)
c3 = (1, 2, 0, 3)                                   # a three-cycle fixing mark 3
def cls(p):                                         # 0 identity, 'v' double transposition, +1/-1 the two 3-cycle classes
    if p == (0, 1, 2, 3): return "e"
    if fix(p) == 0: return "v"
    conj = {mul(mul(g, c3), inv(g)) for g in A4}
    return "+" if p in conj else "-"
chars = {"1": {"e": 1, "v": 1, "+": 1, "-": 1}, "w": {"e": 1, "v": 1, "+": w, "-": w**2},
         "wbar": {"e": 1, "v": 1, "+": w**2, "-": w}, "3": {"e": 3, "v": -1, "+": 0, "-": 0}}
# local monodromies: classes alternate under the quarter turn (recorded U conjugation)
local = ["+", "-", "+", "-"]
# Chevalley-Weil for quotient genus 0: n_V = -dim V + sum_j sum_alpha N_{j,alpha} <-alpha/3>  (+1 for trivial)
def eig_mult(rep, c):                               # multiplicities of 1, w, w^2 for a 3-cycle of class c
    if rep == "3": return [1, 1, 1]
    val = chars[rep][c]
    return [int(val == 1), int(sp.simplify(val-w) == 0), int(sp.simplify(val-w**2) == 0)]
cw = {}
for rep in chars:
    n = -chars[rep]["e"]+(1 if rep == "1" else 0)
    for c in local:
        m = eig_mult(rep, c)
        n += m[1]*Fr(2, 3)+m[2]*Fr(1, 3)
    cw[rep] = str(n)
out["chevalley_weil_multiplicities"] = cw
# Lefschetz: character of H^0(K): identity 5; free elements 1; three-cycles 1 - sum_P 1/(1-zeta_P) conj
genus = 1+12*Fr(-2+4*Fr(2, 3))/2
lef = {"e": genus, "v": 1}
for c in ("+", "-"):
    # each cyclic subgroup of order three fixes one point over every mark; rotation number by local class
    zs = [w if (l == c) else w**2 for l in local]
    lef[c] = sp.simplify(1-sum(1/(1-z) for z in zs))
def mult(rep):
    size = {"e": 1, "v": 3, "+": 4, "-": 4}
    return sp.simplify(sum(size[k]*sp.conjugate(chars[rep][k])*lef[k] for k in size)/12)
out["genus"] = str(genus)
out["lefschetz_multiplicities"] = {rep: str(mult(rep)) for rep in chars}

# ---- invariant operators and the charge ------------------------------------
a, b, c = sp.symbols("a b c")
# A4 commutant on 3 + w + wbar is diag(a,a,a,b,c); the quarter covariance swaps w and wbar: b = c
sol = sp.solve([3*a+2*b], [a], dict=True)[0]
Y = [sp.Rational(-1, 3)]*3+[sp.Rational(1, 2)]*2
out["traceless_covariant_operator"] = {"shape": "diag(a,a,a,b,b)", "a_over_b": str(sol[a]/b),
                                       "normalised": [str(y) for y in Y], "trace": str(sum(Y))}
table = {}
for k in (0, 2, 4):
    for sub in combinations(range(5), k):
        col = sum(1 for s in sub if s < 3); weak = k-col
        y = sum(Y[s] for s in sub)
        key = (k, col, weak, str(y))
        table[key] = table.get(key, 0)+1
out["even_exterior_powers[k,colour_slots,weak_slots,charge,count]"] = [list(k)+[v] for k, v in sorted(table.items())]
out["even_exterior_dimension"] = sum(table.values())

# ---- 2-torsion of J(C), C: eta^2 = 1 - w^6 ---------------------------------
pts = ["P0", "P1", "P2", "Q0", "Q1", "Q2"]            # w = -omega^a (z = 0 side), w = omega^a (z = infinity side)
norm = lambda s: min(frozenset(s), frozenset(pts)-frozenset(s), key=lambda x: (len(x), sorted(x)))
tors = {norm(s) for k in (0, 2) for s in combinations(pts, k)}
assert len(tors) == 16
add2 = lambda s, t: norm(set(s) ^ set(t))
weil = lambda s, t: len(set(s) & set(t)) % 2
rot = lambda s: norm({p[0]+str((int(p[1])+1) % 3) for p in s})
zero = norm(())
lines = set()
for s in tors:
    if s == zero: continue
    grp = frozenset({zero, s, rot(s), rot(rot(s))})
    if len(grp) == 4 and all(add2(x, y) in grp for x in grp for y in grp):
        lines.add(grp)
kinds = []
for g in sorted(lines, key=lambda g: sorted(map(sorted, g))):
    nz = [x for x in g if x != zero]
    iso = all(weil(x, y) == 0 for x in nz for y in nz)
    kinds.append({"elements": [sorted(x) for x in nz], "type": "Lagrangian" if iso else "symplectic"})
out["gamma_invariant_V4_covers"] = {"count": len(lines), "Lagrangian": sum(k["type"] == "Lagrangian" for k in kinds),
                                    "symplectic": sum(k["type"] == "symplectic" for k in kinds), "lines": kinds}
recorded = frozenset({zero, norm({"P0", "P1"}), norm({"P1", "P2"}), norm({"P0", "P2"})})
out["recorded_Y_q_line"] = {"in_list": recorded in lines,
                            "type": "symplectic" if weil(norm({"P0", "P1"}), norm({"P1", "P2"})) else "Lagrangian"}

# ---- theta values at tau = i ------------------------------------------------
mp.mp.dps = 40
n = mp.e**(-mp.pi)
t2, t3, t4 = (mp.jtheta(k, 0, n) for k in (2, 3, 4))
eta = mp.e**(-mp.pi/12)*mp.qp(mp.e**(-2*mp.pi))
out["tau_i"] = {"(theta3/eta)^8": mp.nstr((t3/eta)**8, 20), "(theta2/eta)^8": mp.nstr((t2/eta)**8, 20),
                "(theta4/eta)^8": mp.nstr((t4/eta)**8, 20),
                "D8_vacuum_character": mp.nstr(((t3/eta)**8+(t4/eta)**8)/2, 20),
                "D8_spinor_character": mp.nstr((t2/eta)**8/2, 20)}
print(json.dumps(out, indent=1))
