#!/usr/bin/env python3
"""A positive form on the four marks with its own reference, built from the recorded mark matrix
(18 November 2026).  Hypothesis, not an original producer.

Each mark carries a sign.  A configuration s of the four signs has the weight
        w(s) = exp(-(pi/2) s^T G s),
G the recorded seam matrix at the marks {+-1, +-e^{i delta}} with the recorded zero diagonal
(v484, v485, v512).  Only neutral configurations (sum of the signs zero) are kept: these are the
three pairings, each twice.  Then
        Z_ref   = sum over the three pairings of w,
        Q0[a,b] = sum over the three pairings of w * s_a s_b,
        C       = Q0 / Z_ref  (correlation of the signs).

Sections
    S1  weights of all sixteen configurations at the square
    S2  the form: positive, rank three, neutral; closed forms; reference = E8 character
    S3  values at the square and at the counterwitness; the free seam for comparison
    S4  equal weights (the tetrahedral case): (4/3) I - (1/3) 11^T
"""
import itertools
import sys

import sympy as sp

FAIL = []


def check(name, ok, info=None):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info is not None else ""))
    if not ok:
        FAIL.append(name)


d = sp.symbols("delta", positive=True)
Gk = lambda th: -sp.log(2*sp.sin(th/2))/sp.pi


def G4(delta):
    a, b, o = Gk(delta), Gk(sp.pi-delta), Gk(sp.pi)
    return sp.Matrix([[0, a, o, b], [a, 0, b, o], [o, b, 0, a], [b, o, a, 0]])


def weight(s, delta):
    v = sp.Matrix(s)
    return sp.simplify(sp.exp(-(sp.pi/2)*(v.T*G4(delta)*v)[0]))


# ---------------------------------------------------------------- S1
print("S1  the sixteen sign configurations at the square")
tab = {}
for s in itertools.product((1, -1), repeat=4):
    tab[s] = sp.nsimplify(weight(s, sp.pi/2))
by_total = {}
for s, w in tab.items():
    by_total.setdefault(sum(s), []).append(w)
print("      total sign -> weights:", {k: sorted(set(v)) for k, v in sorted(by_total.items())})
check("neutral: 1 for the alternating signs, 1/4 for the two side pairings; one sign different: 1; all equal: 16",
      sorted(set(by_total[0])) == [sp.Rational(1, 4), 1] and set(by_total[2]) == {1} and set(by_total[4]) == {16})

# ---------------------------------------------------------------- S2
print("S2  the form")
pairings = {"opposite": (1, -1, 1, -1), "side": (1, 1, -1, -1), "side'": (1, -1, -1, 1)}
R = {"opposite": 1/sp.sin(d)**2, "side": sp.tan(d/2)**2/4, "side'": sp.cot(d/2)**2/4}
ok = all(abs(float((weight(s, sp.Rational(9, 10))-R[k].subs(d, sp.Rational(9, 10))).evalf())) < 1e-12 for k, s in pairings.items())
check("weights of the three pairings: 1/sin^2(delta), tan^2(delta/2)/4, cot^2(delta/2)/4", ok)
Zref = sp.simplify(sum(R.values()))
check("reference Z_ref = (4 - sin^2 delta)/(2 sin^2 delta)",
      all(abs(float((Zref-(4-sp.sin(d)**2)/(2*sp.sin(d)**2)).subs(d, dv).evalf())) < 1e-12 for dv in (sp.Rational(2, 5), sp.Rational(9, 10), sp.Rational(13, 10))))
chi = sp.Rational(2)**sp.Rational(5, 3)*(16*sp.sin(d)**2)**sp.Rational(1, 3)*Zref
check("Z_ref times 2^(5/3) |Disc|^(1/6) is the E8 character 4 (4 - sin^2 delta)/sin^(4/3) delta (round nineteen)",
      abs(float((chi-4*(4-sp.sin(d)**2)/sp.sin(d)**sp.Rational(4, 3)).subs(d, sp.Rational(7, 10)).evalf())) < 1e-12)
Q0 = sp.zeros(4, 4)
for k, s in pairings.items():
    v = sp.Matrix(s)
    Q0 += R[k]*v*v.T
C = (Q0/Zref).applyfunc(sp.simplify)
Cn = C.subs(d, sp.Rational(9, 10)).evalf()
ev = sorted(float(sp.re(e)) for e in Cn.eigenvals(multiple=True))
check("C is positive of rank three, annihilates the symmetric sum, and has unit diagonal",
      abs(ev[0]) < 1e-12 and ev[1] > 0 and all(abs(float(Cn[i, i])-1) < 1e-12 for i in range(4)), [round(e, 6) for e in ev])
c_opp = sp.simplify(C[0, 2])
check("sign correlation of opposite marks: sin^2(delta)/(4 - sin^2 delta)",
      abs(float((c_opp-sp.sin(d)**2/(4-sp.sin(d)**2)).subs(d, sp.Rational(7, 10)).evalf())) < 1e-12)

# ---------------------------------------------------------------- S3
print("S3  values")
Csq = C.subs(d, sp.pi/2).applyfunc(sp.nsimplify)
check("square: C = circ(1, -2/3, 1/3, -2/3); weights of the pairings 2/3, 1/6, 1/6",
      Csq == sp.Matrix([[1, sp.Rational(-2, 3), sp.Rational(1, 3), sp.Rational(-2, 3)],
                        [sp.Rational(-2, 3), 1, sp.Rational(-2, 3), sp.Rational(1, 3)],
                        [sp.Rational(1, 3), sp.Rational(-2, 3), 1, sp.Rational(-2, 3)],
                        [sp.Rational(-2, 3), sp.Rational(1, 3), sp.Rational(-2, 3), 1]]))
dcw = sp.atan(sp.Rational(4, 3))
Ccw = C.subs(d, dcw).applyfunc(lambda e: sp.nsimplify(sp.simplify(e)))
print("      counterwitness: first row of C =", list(Ccw.row(0)))
check("counterwitness: opposite marks 4/21", Ccw[0, 2] == sp.Rational(4, 21), Ccw[0, 2])
ratio = sp.nsimplify((Zref.subs(d, dcw)/Zref.subs(d, sp.pi/2)).evalf(30), rational=True)
check("reference at the counterwitness against the square: 7/4", ratio == sp.Rational(7, 4), ratio)
free = sp.Matrix(pairings["opposite"])*sp.Matrix(pairings["opposite"]).T
check("free seam (only the alternating signs): C = (+1, -1, +1, -1) pattern of rank one, opposite marks +1",
      free.rank() == 1 and free[0, 2] == 1 and free[0, 1] == -1)

# ---------------------------------------------------------------- S4
print("S4  equal weights")
Ceq = sp.zeros(4, 4)
for s in pairings.values():
    v = sp.Matrix(s)
    Ceq += v*v.T/3
check("three equal weights: C = (4/3) I - (1/3) 11^T, the member c = -1/3 of the family (1-c) I + c 11^T",
      Ceq == sp.Rational(4, 3)*sp.eye(4)-sp.ones(4, 4)/3)

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
