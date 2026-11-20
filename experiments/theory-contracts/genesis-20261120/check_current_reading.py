#!/usr/bin/env python3
"""Do the (10,6) currents of the record measure a charge on the family label, or do they only move states?

Uses the record's own exact response functions (seam-collar-multifield-20261016/e8_boundary_response.py, read only)
and the mark formulas of the E8 word contract (Z0, ZU, ZY, ZUY as quadratic traces). The label f is +1 on the three
family words of one chirality, -1 on their mirrors, 0 on the fourth family.
"""
import json, sys
from fractions import Fraction as Fr
from itertools import product
from pathlib import Path
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"seam-collar-multifield-20261016"))
import e8_boundary_response as rec

H = sp.Rational(1, 2)
words = [w for w in product((1, -1), repeat=8) if sp.prod(w) == 1]
three = {w: 0 if w[5] == w[6] == w[7] else w[5]*w[6]*w[7] for w in words}
four = {w: w[5]*w[6]*w[7] for w in words}
vec = lambda w: tuple(H*x for x in w)
unit = lambda i, s=1: tuple(sp.Integer(s) if k == i else sp.S.Zero for k in range(8))
checks = []
def ok(name, cond, value=None):
    checks.append({"check": name, "passed": bool(cond), "value": value}); assert cond, (name, value)

# 1. Cartan current of the record's net on a matter pair: residue of <H_d(u) V_w(1) V_-w(0)> at u = 1 is d.w
u = sp.symbols("u")
def cartan_residue(d, w):
    p = vec(w); n = rec.adjoint_sign(p)*rec.vertex_response((p, rec.neg(p)), (1, 0))
    val = rec.adjoint_sign(p)*sum(rec.dot(d, a)/(u-z) for a, z in ((p, 1), (rec.neg(p), 0)))*rec.vertex_response((p, rec.neg(p)), (1, 0))
    return sp.residue(val/n, u, 1)
sample = words[5]
ok("record Ward rule: Cartan current in direction d has residue d.w on the matter pair", cartan_residue(unit(6), sample) == H*sample[6]
   and rec.current_response((vec(sample), rec.neg(vec(sample))), (1, 0), ((unit(6), 3),)) /
       rec.vertex_response((vec(sample), rec.neg(vec(sample))), (1, 0)) == H*sample[6]*(sp.Rational(1, 2)-sp.Rational(1, 3)))
first = {i: sum(three[w]*H*w[i] for w in words) for i in range(8)}
ok("three-family label under the eight Cartan currents of the net: carrier 0, partner -16 each",
   [first[i] for i in range(8)] == [0]*5+[-16]*3, [str(first[i]) for i in range(8)])

# 2. a (10,6) current pair: coefficient of 1/(x-y) in <V_a(x) V_-a(y) V_w(1) V_-w(0)> / (leading term) is the charge a.w
x, y = sp.symbols("x y")
a = tuple(sp.Integer(v) for v in (1, 0, 0, 0, 0, 1, 0, 0))
def pair_charge(w):
    p = vec(w); k = int(rec.dot(a, p))
    # exact four-point function of the record: (x-y)^-2 * ((x-1)/(y-1))^k * (y/x)^k * const
    ratio = ((x-1)/(y-1))**k*(y/x)**k
    lead = sp.series(ratio.subs(x, y+sp.Symbol("e")), sp.Symbol("e"), 0, 2).removeO().coeff(sp.Symbol("e"), 1)
    return sp.simplify(lead/(1/(y-1)-1/y))
num = rec.vertex_response((a, rec.neg(a), vec(sample), rec.neg(vec(sample))), (5, 4, 1, 0))
k = int(rec.dot(a, vec(sample)))
ok("record four-point function has the form used", sp.simplify(abs(num) - sp.Rational(1)*abs((5-4)**-2*((5-1)*sp.Rational(1, 4-1))**k*(sp.Rational(4, 5))**k)) == 0, str(num))
ok("the pair V_a V_-a of a (10,6) current measures the charge a.w at first order", all(pair_charge(w) == rec.dot(a, vec(w)) for w in words[:16]))
ok("summed with the three-family label this charge is -16; with the four-family label 0",
   sum(three[w]*rec.dot(a, vec(w)) for w in words) == -16 and sum(four[w]*rec.dot(a, vec(w)) for w in words) == 0)

# 3. the mark contract of the record: Z0, ZU, ZY, ZUY are traces of even order with a carrier marker
X6 = (-2, -2, -2, 3, 3, 0, 0, 0)                                    # 6 Y, marker of the recorded word
T3 = (0, 0, 0, 1, -1, 0, 0, 0)
currents = [tuple((si if k == i else sj if k == j else 0) for k in range(8)) for i in range(5) for j in range(5, 8)
            for si in (1, -1) for sj in (1, -1)]
def move(w, r):                                                      # w + r if it is a glue vector
    return tuple(-s if c else s for s, c in zip(w, r)) if sum(c*s for c, s in zip(r, w)) == -2 else None
def odd_marks(f, r, h):
    hv = lambda w: Fr(sum(c*s for c, s in zip(h, w)), 2)
    z0 = sum(f.values()); zu = sum(f[w]*hv(w)**2 for w in words)
    zy = sum(f[w] for w in words if move(w, r)); zuy = sum(f[w]*hv(w)*hv(move(w, r)) for w in words if move(w, r))
    return z0, zu, zy, zuy
single = {odd_marks(three, r, X6) for r in currents}
herm = {tuple(p+q for p, q in zip(odd_marks(three, r, X6), odd_marks(three, tuple(-c for c in r), X6))) for r in currents}
ok("odd part of Z0 and ZU on the three-family label vanishes for the carrier markers Y and T3",
   all(odd_marks(three, r, h)[:2] == (0, 0) for r in currents for h in (X6, T3)))
ok("a single root current gives odd ZY = +-8 on the label (it is not self-adjoint)", {m[2] for m in single} == {8, -8}, sorted(map(str, single)))
ok("a current together with its adjoint gives odd ZY = 0 and odd ZUY = 0 for the carrier marker", herm == {(0, 0, 0, 0)}, sorted(map(str, herm)))
fam = (0, 0, 0, 0, 0, 1, 1, 1)
ok("a family marker (1,1,1) in the partner sees the label at first order: -48",
   sum(three[w]*Fr(sum(c*s for c, s in zip(fam, w)), 2) for w in words) == -48)

out = {"status": "CURRENT_READING_CHECKED", "physical_solution": False, "checks": checks,
       "verdict": {"record_mark_contract": "markers are Cartan currents of the carrier (6Y, T3); the (10,6) current enters only quadratically, as Tr P X_-g X_g over both chiralities; the mirror-odd part is never formed. Formed separately it is +-8 for one current and cancels between g and -g; no quantity of order one or three in a partner direction occurs",
                   "stipulated_E8_net": "all eight Cartan currents and the pairs V_a V_-a are local fields with exact Ward rule; partner Cartan currents and every (10,6) pair see the three-family label at first order (-16)",
                   "consequence": "three blind chiral families hold in the mark contract and fail in the full local net; the record already states that both are different algebras and that the compression to 96 states is no homomorphism"}}
print(json.dumps(out, indent=1, ensure_ascii=False))
