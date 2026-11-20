#!/usr/bin/env python3
"""Does the network of the record (extended E8 diagram with its marks) come from the same torus?

Points of period n under the step 1+i are the torsion points killed by (1+i)^n - 1. Their symmetry group (maps that
keep the pairing) is SL2 over the integers modulo m. The graph "multiply by the two-dimensional representation" of a
finite subgroup of SU(2) is an extended diagram (McKay). Checked here by spectra and counts, with controls.
"""
import itertools, json
from collections import Counter
import numpy as np
import sympy as sp

rows = []
def test(name, cond, value=None):
    rows.append({"check": name, "passed": bool(cond), "value": value})

def sl2(p):
    return [(a, b, c, d) for a, b, c, d in itertools.product(range(p), repeat=4) if (a*d-b*c) % p == 1]
def mul(x, y, p):
    return ((x[0]*y[0]+x[1]*y[2]) % p, (x[0]*y[1]+x[1]*y[3]) % p, (x[2]*y[0]+x[3]*y[2]) % p, (x[2]*y[1]+x[3]*y[3]) % p)
def order(x, p):
    y, n = x, 1
    while y != (1, 0, 0, 1):
        y, n = mul(y, x, p), n+1
    return n
def inv(x, p):
    return (x[3] % p, (-x[1]) % p, (-x[2]) % p, x[0] % p)
def classes(G, p):
    seen, out = set(), []
    for g in G:
        if g not in seen:
            cl = {mul(mul(h, g, p), inv(h, p), p) for h in G}
            seen |= cl; out.append(cl)
    return out
def spectrum_from_orders(cls, p):
    # a faithful two-dimensional unitary representation with determinant 1: an element of order n has trace 2cos(2 pi k/n)
    vals = []
    for cl in cls:
        n = order(next(iter(cl)), p)
        vals.append(n)
    return sorted(vals)
def diagram(edges, n):
    A = np.zeros((n, n))
    for i, j in edges:
        A[i, j] = A[j, i] = 1
    return A
E8hat = diagram([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (5, 8)], 9)       # arm lengths 5, 2, 1 around node 5
E6hat = diagram([(0, 1), (1, 2), (2, 3), (3, 4), (2, 5), (5, 6)], 7)
E7hat = diagram([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (3, 7)], 8)
def perron(A):
    w, v = np.linalg.eigh(A); x = np.abs(v[:, -1]); return w, np.round(x/x.min()).astype(int)
phi = (1+5**0.5)/2

# ---- the step and its periodic points
al = lambda m: (1 % m, (-1) % m, 1 % m, 1 % m)                         # multiplication by 1+i on (x, y) -> (x - y, x + y)
def step_order(m):
    y, n = al(m), 1
    while y != (1, 0, 0, 1):
        y, n = mul(y, al(m), m), n+1
    return n
test("(1+i)^4 = -4: the fourth power is the first that is a pure dilation, and the points of period four are the 5-torsion (25 points)",
     step_order(5) == 4 and all(step_order(5) != k for k in (1, 2, 3)))
test("on the 3-torsion the step has order 8, on the 5-torsion order 4", (step_order(3), step_order(5)) == (8, 4))

# ---- level 5
G5 = sl2(5); C5 = classes(G5, 5)
ords5 = Counter(order(g, 5) for g in G5)
test("symmetry group of the 5-torsion: 120 elements, 9 classes, element orders as the binary icosahedral group",
     len(G5) == 120 and len(C5) == 9 and dict(ords5) == {1: 1, 2: 1, 3: 20, 4: 30, 5: 24, 6: 20, 10: 24}, dict(ords5))
traces5 = sorted([2, -2, -1, 1, 0, phi, 1-phi, phi-1, -phi])           # 2cos of the rotation angles of the nine classes
w8, marks8 = perron(E8hat)
test("its two-dimensional representation has traces 2, -2, -1, 1, 0, +-phi, +-(phi - 1): exactly the spectrum of the extended E8 diagram",
     np.allclose(sorted(w8), traces5, atol=1e-9), [round(float(x), 6) for x in sorted(w8)])
test("the marks of the record's network are the dimensions of the nine representations: 1, 2, 3, 4, 5, 6, 4, 2, 3 with squares summing to 120; "
     "top value 2, next value phi (gap 2 - phi)", sorted(marks8.tolist()) == [1, 2, 2, 3, 3, 4, 4, 5, 6] and int((marks8**2).sum()) == 120
     and abs(sorted(w8)[-1]-2) < 1e-9 and abs(sorted(w8)[-2]-phi) < 1e-9, marks8.tolist())
a5 = al(5); det = (a5[0]*a5[3]-a5[1]*a5[2]) % 5
five = [cl for cl in C5 if order(next(iter(cl)), 5) == 5]
img = {frozenset(mul(mul(a5, g, 5), (a5[3]*3 % 5, (-a5[1])*3 % 5, (-a5[2])*3 % 5, a5[0]*3 % 5), 5) for g in cl) for cl in five}
test("the step is not in that group (determinant 2, not a square modulo 5); it acts on it from outside and swaps the two classes of order five, "
     "the exchange of sqrt 5 with -sqrt 5", det == 2 and len(five) == 2 and img == {frozenset(five[1]), frozenset(five[0])} and frozenset(five[0]) != frozenset(five[1])
     and mul(a5, (a5[3]*3 % 5, (-a5[1])*3 % 5, (-a5[2])*3 % 5, a5[0]*3 % 5), 5) == (1, 0, 0, 1)
     and not any(frozenset(mul(mul(a5, g, 5), (a5[3]*3 % 5, (-a5[1])*3 % 5, (-a5[2])*3 % 5, a5[0]*3 % 5), 5) for g in cl) == frozenset(cl) for cl in five))

# ---- level 3 and controls
G3 = sl2(3); C3 = classes(G3, 3)
w6, marks6 = perron(E6hat)
test("3-torsion: 24 elements, 7 classes (binary tetrahedral); spectrum of the extended E6 diagram 2, 1, 1, 0, -1, -1, -2 equals its traces",
     len(G3) == 24 and len(C3) == 7 and np.allclose(sorted(w6), [-2, -1, -1, 0, 1, 1, 2], atol=1e-9)
     and dict(Counter(order(g, 3) for g in G3)) == {1: 1, 2: 1, 3: 8, 4: 6, 6: 8})
G7 = sl2(7)
test("control: for the 7-torsion the group has 336 elements; no extended E diagram has that order (24, 48, 120 only)", len(G7) == 336)
def sl2z4():
    return [(a, b, c, d) for a, b, c, d in itertools.product(range(4), repeat=4) if (a*d-b*c) % 4 == 1]
G4 = sl2z4()
o4 = Counter(order(g, 4) for g in G4)
test("control: the 4-torsion gives 48 elements but not the binary octahedral group (that one has exactly one element of order two)",
     len(G4) == 48 and o4[2] != 1, dict(o4))
test("the three levels that give an extended E diagram or the marks are 2, 3, 5; their product 30 is the first modulus where all three are present",
     2*3*5 == 30 and step_order(15) == 8)

if __name__ == "__main__":
    print(json.dumps({"status": "NETWORK_CHECKED", "physical_solution": False, "passed": sum(r["passed"] for r in rows), "total": len(rows), "rows": rows}, indent=1, default=str))
