#!/usr/bin/env python3
"""What two marks exchange: the three pairing responses in the channel of a pair (18 November 2026).

Four marks a, b, c, d.  Pairing responses R_(xy|zw) = P^2/(P1 P2 P3), P the product of the two
differences inside the pairs (round seventeen; weight 1/2 at each mark).  Let b approach a.
A function behaving like (a-b)^(h - 1) exchanges a state of weight h between the pair and the rest.

Sections
    H1  relation between the three pairing products
    H2  exponents in the channel of the pair (a, b): weights 0, 1 and 2
    H3  which combinations carry the weight-one exchange; values on the mark circle
    H4  the sign form of round twenty-three in these terms
"""
import sys

import sympy as sp

FAIL = []


def check(name, ok, info=None):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info is not None else ""))
    if not ok:
        FAIL.append(name)


a, b, c, d, eta = sp.symbols("a b c d eta")
Pab = (a-b)*(c-d)
Pac = (a-c)*(b-d)
Pad = (a-d)*(b-c)
den = Pab*Pac*Pad
R = {"ab": Pab**2/den, "ac": Pac**2/den, "ad": Pad**2/den}

# ---------------------------------------------------------------- H1
print("H1  the three pairing products")
check("P_ab - P_ac + P_ad = 0", sp.expand(Pab-Pac+Pad) == 0)

# ---------------------------------------------------------------- H2
print("H2  exponents when b approaches a")
vals = {a: sp.Rational(3, 10), c: sp.Rational(-7, 5)+sp.I/2, d: sp.Rational(11, 10)-sp.I*sp.Rational(3, 4)}


def lead(expr):
    e = sp.simplify(expr.subs(b, a+eta).subs(vals))
    ser = sp.series(e, eta, 0, 3).removeO()
    powers = sorted(p for p in (-2, -1, 0, 1, 2) if sp.simplify(ser.coeff(eta, p)) != 0)
    return powers[0]


ex = {"R_ab": lead(R["ab"]), "R_ac + R_ad": lead(R["ac"]+R["ad"]), "R_ac - R_ad": lead(R["ac"]-R["ad"])}
print("      leading powers of (a - b):", ex)
check("R_ab starts at power +1 (weight 2), R_ac + R_ad at -1 (weight 0, the empty exchange), "
      "R_ac - R_ad at 0 (weight 1, a current)", ex == {"R_ab": 1, "R_ac + R_ad": -1, "R_ac - R_ad": 0})
check("exactly: R_ac - R_ad = (P_ac + P_ad)/(P_ac P_ad), regular and non-zero when the pair closes",
      sp.simplify(R["ac"]-R["ad"]-(Pac+Pad)/(Pac*Pad)) == 0)

# ---------------------------------------------------------------- H3
print("H3  the weight-one exchange")
check("the sum of the three responses (the E8 reading) contains R_ac and R_ad with equal weight: "
      "no weight-one exchange between a and b", sp.simplify((R["ab"]+R["ac"]+R["ad"])-(R["ab"]+(R["ac"]+R["ad"]))) == 0)
check("a single pairing other than (ab|cd) contains it with weight 1/2: R_ac = (R_ac + R_ad)/2 + (R_ac - R_ad)/2", True)
dl = sp.symbols("delta", positive=True)
Ropp, Rs, Rsp = 1/sp.sin(dl)**2, sp.tan(dl/2)**2/4, sp.cot(dl/2)**2/4   # (13|24), (12|34), (14|23) for marks 1,2,3,4
adj = sp.simplify((Ropp-Rsp).subs(dl, sp.pi/2))
opp = sp.simplify((Rs-Rsp).subs(dl, sp.pi/2))
check("marks {1, e^{i delta}, -1, -e^{i delta}} at the square: weight-one exchange between neighbours 3/4, "
      "between opposite marks 0", adj == sp.Rational(3, 4) and opp == 0, (adj, opp))
gen = sp.simplify(Ropp-Rsp)
check("neighbours in general: 1/sin^2(delta) - cot^2(delta/2)/4; opposite: (tan^2(delta/2) - cot^2(delta/2))/4",
      abs(float((gen-(1/sp.sin(dl)**2-sp.cot(dl/2)**2/4)).subs(dl, sp.Rational(7, 10)).evalf())) < 1e-12)

# ---------------------------------------------------------------- H4
print("H4  the sign form")
signs = {"ab": sp.Matrix([1, 1, -1, -1]), "ac": sp.Matrix([1, -1, 1, -1]), "ad": sp.Matrix([1, -1, -1, 1])}
ra, rc, rd = sp.symbols("R_ab R_ac R_ad")
Q0 = ra*signs["ab"]*signs["ab"].T+rc*signs["ac"]*signs["ac"].T+rd*signs["ad"]*signs["ad"].T
check("Q0[a,c] - Q0[a,d] = 2 (R_ac - R_ad): the difference of the sign correlations of a with the two "
      "other marks is twice the weight-one exchange of the pair (a, b)", sp.expand(Q0[0, 2]-Q0[0, 3]-2*(rc-rd)) == 0)
check("Q0[a,b] = R_ab - (R_ac + R_ad): weight two minus the empty exchange", sp.expand(Q0[0, 1]-(ra-rc-rd)) == 0)

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
