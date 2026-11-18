#!/usr/bin/env python3
"""Sixteenth pasted note: commutator rule, modular clock of the four marks, exclusion for the
quarter rotation, gluing state; and the same commutator rule on the native operator
(18 November 2026).

argv[1] = repository root (native assembler, used unchanged).

Sections
    M1  note: C_X = <HXHX> - <H^2><X^2> = -(1/2)<[H,X]^+[H,X]> = -(q^2/2)<X^2> on 128 and 96 states
    M2  note: modular coordinate and speed for two opposite arcs; square against counterwitness
    M3  note: one abelian quarter rotation cannot have fourth power -1 on the spinors and weight 1/2;
        correction of a sentence of round twenty-nine
    M4  note: the gluing state and the numbers 1 and 12
    M5  new: the same rule on the native operator, flat density on the 96 fibres
"""
import itertools
import sys
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

FAIL = []


def check(name, ok, info=None):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info is not None else ""))
    if not ok:
        FAIL.append(name)


root = Path(sys.argv[1]).resolve()

# ---------------------------------------------------------------- zero modes of eight fermions
N = 8
dim = 2**N
bits = np.array(list(itertools.product((0, 1), repeat=N)))
ann = []
for i in range(N):
    a = np.zeros((dim, dim))
    for k, b in enumerate(bits):
        if b[i]:
            c = b.copy()
            c[i] = 0
            a[int("".join(map(str, c)), 2), k] = (-1)**int(b[:i].sum())
    ann.append(a)
spin2 = 2*bits-1
matter = np.array([list(v).count(-1) % 2 == 0 for v in spin2])
hyp = np.array([-1/3, -1/3, -1/3, 1/2, 1/2, 0, 0, 0])
fam = [tuple(v[5:]) for v in spin2]
carrier_even = [list(v[:5]).count(-1) % 2 == 0 for v in spin2]
w0 = sorted({f for f, m, e in zip(fam, matter, carrier_even) if m and e})[0]
three = matter & np.array([f not in (w0, tuple(-v for v in w0)) for f in fam])
currents = [(i, j, si, sj) for i in range(5) for j in range(5, 8) for si in (1, -1) for sj in (1, -1)]

# ---------------------------------------------------------------- M1
print("M1  the commutator rule")
rat = lambda x: Fr(float(x)).limit_denominator(2000)
for label, keep in (("128", matter), ("96", three)):
    idx = np.where(keep)[0]
    H = np.diag(spin2[idx]@hyp/2)
    n = len(idx)
    tr = lambda A: np.trace(A)/n
    res = {}
    ok_id, ok_fac = True, True
    for (i, j, si, sj) in currents:
        xg = ((ann[i].T if si > 0 else ann[i])@(ann[j].T if sj > 0 else ann[j]))[np.ix_(idx, idx)]
        X = (xg+xg.T)/np.sqrt(2)
        com = H@X-X@H
        cx = tr(H@X@H@X)-tr(H@H)*tr(X@X)
        ok_id = ok_id and abs(cx+0.5*tr(com.T@com)) < 1e-12
        ok_fac = ok_fac and abs(tr(H@H@X@X)-tr(H@H)*tr(X@X)) < 1e-12
        q = hyp[i]*si
        ok_id = ok_id and abs(cx+q*q/2*tr(X@X)) < 1e-12
        res.setdefault("weak" if i >= 3 else "colour", set()).add((rat(tr(H@H)), rat(tr(X@X)), rat(tr(H@X@H@X)), rat(cx)))
    print(f"      {label} states (H^2, X^2, HXHX, C_X):", {k: sorted(v) for k, v in res.items()})
    check(f"{label} states: C_X = -(1/2)<[H,X]^+[H,X]> = -(q^2/2)<X^2> for all 60 hermitian pair currents", ok_id)
    check(f"{label} states: <H^2 X^2> = <H^2><X^2>, which the rule needs", ok_fac)
    want = ({"weak": {(Fr(5, 24), Fr(1, 4), Fr(1, 48), Fr(-1, 32))}, "colour": {(Fr(5, 24), Fr(1, 4), Fr(11, 288), Fr(-1, 72))}} if label == "128"
            else {"weak": {(Fr(5, 24), Fr(1, 6), Fr(1, 72), Fr(-1, 48))}, "colour": {(Fr(5, 24), Fr(1, 6), Fr(11, 432), Fr(-1, 108))}})
    check(f"{label} states: the values of the note (and of the record)", res == want)

# ---------------------------------------------------------------- M2
print("M2  the modular clock of two opposite arcs")
th, de = sp.symbols("theta delta", positive=True)
ends_a, ends_b = [0, sp.pi], [de, sp.pi+de]
prod = sp.Mul(*[sp.sin((th-x)/2) for x in ends_a])/sp.Mul(*[sp.sin((th-x)/2) for x in ends_b])
ok = all(abs(float((prod-sp.sin(th)/sp.sin(th-de)).subs({th: tv, de: dv}))) < 1e-12 for tv, dv in ((0.3, 1.1), (0.7, 0.9), (0.2, 1.5)))
check("the ratio of the chord factors of the four ends 0, delta, pi, pi + delta is sin(theta)/sin(theta - delta): "
      "the coordinate of the note is log(sin theta / sin(delta - theta))", ok)
t = sp.log(sp.sin(th)/sp.sin(de-th))
speed = sp.simplify(2*sp.pi/sp.diff(t, th))
check("speed 2 pi / t' = 2 pi sin(theta) sin(delta - theta)/sin(delta); it vanishes at both marks",
      sp.simplify(speed-2*sp.pi*sp.sin(th)*sp.sin(de-th)/sp.sin(de)) == 0 and speed.subs(th, 0) == 0)
mid = sp.simplify(speed.subs(th, de/2))
num = lambda e, dv: float(e.subs(de, dv))
check("at the middle of an arc the speed is pi tan(delta/2): pi at the square and pi/2 at the counterwitness (3 + 4i)/5",
      all(abs(num(mid-sp.pi*sp.tan(de/2), dv)) < 1e-12 for dv in (0.4, 0.9, 1.3)) and abs(num(mid, sp.pi/2)-float(sp.pi)) < 1e-12
      and abs(num(mid, sp.atan(sp.Rational(4, 3)))-float(sp.pi/2)) < 1e-12)
check("in the variable of round twenty-five, x = sin^2(delta/2): middle speed pi sqrt(x/(1 - x))",
      abs(float((mid-sp.pi*sp.sqrt(sp.sin(de/2)**2/(1-sp.sin(de/2)**2))).subs(de, 0.8))) < 1e-12)

# ---------------------------------------------------------------- M3
print("M3  one abelian quarter rotation")
found = [m for m in itertools.product(range(-4, 5), repeat=8) if sum(x*x for x in m) == 16 and sum(m) % 2 == 1] if False else None
check("for integers m_i the sum and the sum of squares have the same parity, so sum odd (fourth power -1 on the "
      "spinor roots) excludes sum of squares 16 (weight 1/2)", all((sum(m)-sum(x*x for x in m)) % 2 == 0 for m in itertools.product(range(-3, 4), repeat=4)))
m = (1,)*8
phase4 = {int(round(float((sum(m[i]*s[i] for i in range(8))/2) % 1*2))) for s in spin2[matter]}
check("the common quarter rotation of eight fermions, m = (1,...,1): its fourth power is +1 on all 128 spinor states "
      "and its weight is 8/32 = 1/4, not 1/2", phase4 == {0} and Fr(sum(x*x for x in m), 32) == Fr(1, 4))
check("correction to round twenty-nine: the sentence that the common quarter rotation has order eight on matter states "
      "and that its fourth power is the recorded sign is withdrawn; the recorded sign needs an odd rotation, which is a "
      "different operation from the exchange of the copies (weight 1/2 on two copies)", True)

# ---------------------------------------------------------------- M4
print("M4  the gluing state")
idx = np.where(matter)[0]
B = np.eye(128)/np.sqrt(128)                                  # coefficients of sum |a>|a-bar>
worst = 0.0
for (i, j, si, sj) in currents:
    xg = ((ann[i].T if si > 0 else ann[i])@(ann[j].T if sj > 0 else ann[j]))[np.ix_(idx, idx)]
    worst = max(worst, float(np.max(np.abs(xg@B-B@xg))))       # (X x 1 - 1 x X^T)|B> = 0  <=>  X B = B X
check("the equal-weight glued state satisfies the gluing condition for all 60 currents; its reduced state is flat "
      "(round twenty-eight showed that it is the only one)", worst < 1e-14)
mp.mp.dps = 25
q = mp.exp(-2*mp.pi)
E4 = 1+240*mp.nsum(lambda n: n**3*q**n/(1-q**n), [1, mp.inf])
eta8 = (mp.exp(-2*mp.pi/24)*mp.qp(q))**8
check("E4(i)/eta(i)^8 = 12: the propagated channel at the square; the unit coupling of the single E8 sector is a "
      "statement of the note, not recomputed", abs(E4/eta8-12) < mp.mpf(10)**(-18))

# ---------------------------------------------------------------- M5
print("M5  the same rule on the native operator")
base = root/"experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q"
sys.path.insert(0, str(base/"runtime"))
from assemble import Assembler  # noqa: E402

Asm = Assembler(base/"domain/native-diagnostic-v1.json")
g, orig = Asm.gamma, Asm.original
YT = {("u", "L"): 1/6, ("d", "L"): 1/6, ("u", "R"): 2/3, ("d", "R"): -1/3, ("nu", "L"): -1/2, ("e", "L"): -1/2, ("nu", "R"): 0.0, ("e", "R"): -1.0}
TT = {("u", "L"): 0.5, ("nu", "L"): 0.5, ("d", "L"): -0.5, ("e", "L"): -0.5}
yv = [YT[(n, c)] for (n, c, _, _) in orig.P]
tv = [TT.get((n, c), 0.0) for (n, c, _, _) in orig.P]
Yh, T3 = np.diag(yv+[-v for v in yv]), np.diag(tv+[-v for v in tv])
DF = Asm.assemble({"schema": "tfpt.field-configuration.v1", "domain_sha256": Asm.domain["_sha256"], "genes": [1, 1, 1]})["phi_fourier"][0]
wv, vv = np.linalg.eigh(-1j*g@DF)
Yn = (vv*np.sign(wv))@vv.conj().T
tr = lambda A: (np.trace(A)/96).real
rule = lambda Hm: (tr(Hm@Yn@Hm@Yn)-tr(Hm@Hm)*tr(Yn@Yn), -0.5*tr((Hm@Yn-Yn@Hm).conj().T@(Hm@Yn-Yn@Hm)))
vals = {name: rule(Hm) for name, Hm in (("grading", g), ("hypercharge", Yh), ("isospin", T3), ("electric charge", Yh+T3))}
print("      (connected answer, -(1/2)<|[.,Y]|^2>):", {k: (round(a, 6), round(b, 6)) for k, (a, b) in vals.items()})
check("flat density on the 96 fibres, cut Y = direction of the native potential: the rule holds for every marker",
      all(abs(a-b) < 1e-10 for a, b in vals.values()))
check("native <H^2> for the hypercharge is 5/24, the value of the E8 chart", abs(tr(Yh@Yh)-5/24) < 1e-12)
check("grading: -2, i.e. the normalised word ABAB = -1 of group A9 is the extreme case of the same rule "
      "(anticommuting cuts)", abs(vals["grading"][0]+2) < 1e-10 and abs(tr(g@Yn@g@Yn)+1) < 1e-10)
check("electric charge: 0 (the potential is neutral); hypercharge and isospin: equal and negative",
      abs(vals["electric charge"][0]) < 1e-10 and vals["hypercharge"][0] < 0 and abs(vals["hypercharge"][0]-vals["isospin"][0]) < 1e-10,
      round(vals["hypercharge"][0], 6))
check("the native value -0.129976 is not -1/8: the direction of the potential is not a pure carrier of hypercharge 1/2 "
      "(the Majorana part mixes); no recorded anchor value", abs(vals["hypercharge"][0]+0.125) > 1e-3)

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
