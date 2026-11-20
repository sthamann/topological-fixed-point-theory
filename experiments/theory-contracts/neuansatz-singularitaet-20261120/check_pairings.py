#!/usr/bin/env python3
"""Are the three steps of degree two out of the torus the three pairings of the external report B, and does the
torus itself move the distinguished one cyclically?

A step of degree two has a mark t as kernel. It sends t to the origin and merges the other two marks: one
distinguished channel and a coupled pair. The symmetries of the torus permute the three marks.
"""
import itertools, json
import sympy as sp

checks = []
def ok(name, cond, value=None):
    checks.append({"check": name, "passed": bool(cond), "value": value}); assert cond, (name, value)

marks = [(1, 0), (0, 1), (1, 1)]                       # a = 1/2, b = tau/2, c = (1+tau)/2 in the basis (1, tau), mod 2
def act(mat):                                          # action of an integer matrix on the three marks
    img = [((mat[0][0]*x+mat[0][1]*y) % 2, (mat[1][0]*x+mat[1][1]*y) % 2) for x, y in marks]
    return tuple(marks.index(p) for p in img)
# multiplication by i on Z + iZ: 1 -> i, i -> -1 ; multiplication by omega on Z + omega Z: 1 -> omega, omega -> -1 - omega
quarter = [[0, -1], [1, 0]]; third = [[0, -1], [1, -1]]
ok("each step of degree two: kernel mark goes to the origin, the other two are merged (one channel and a pair)",
   all(len({((x+t[0]) % 2, (y+t[1]) % 2) for x, y in marks if (x, y) != t} | {s for s in marks if s != t}) == 2 for t in marks))
ok("square torus: the symmetry of order four swaps two marks and fixes the third", act(quarter) == (1, 0, 2), act(quarter))
ok("hexagonal torus: the symmetry of order three moves the three marks cyclically", sorted(act(third)) == [0, 1, 2]
   and all(act(third)[k] != k for k in range(3)), act(third))
S = [[0, -1], [1, 0]]; ST = [[0, -1], [1, 1]]; Tm = [[1, 1], [0, 1]]
ok("modular group on the marks: S and T give swaps, ST gives the cycle of length three; all six permutations occur",
   len(set(act(S))) == 3 and sum(act(S)[k] == k for k in range(3)) == 1 and sum(act(Tm)[k] == k for k in range(3)) == 1
   and all(act(ST)[k] != k for k in range(3)))
lam = sp.symbols("lam")
j = 256*(lam**2-lam+1)**3/(lam**2*(lam-1)**2)
l0 = sp.Rational(1, 2)+sp.sqrt(3)*sp.I/2
ok("in the position of the fourth mark: the cycle fixes exp(+-i pi/3) with j = 0; the swap fixes -1, 2, 1/2 with j = 1728",
   sp.expand(l0**2-l0+1) == 0 and all(sp.simplify(j.subs(lam, v)) == 1728 for v in (-1, 2, sp.Rational(1, 2)))
   and sp.expand((1-l0)*l0-1) == 0 and all(sp.simplify(w-v) == 0 for v, w in ((-1, 1/sp.Integer(-1)), (2, sp.Integer(2)/(2-1)), (sp.Rational(1, 2), 1-sp.Rational(1, 2)))))

# the family matrix of the record (verification/v395): distinguished channel 0, pair (1, 2)
M = sp.Matrix([[1, 0, 0], [0, sp.Rational(1, 2), sp.Rational(1, 6)], [0, sp.Rational(1, 6), sp.Rational(1, 2)]])
swap = sp.Matrix([[1, 0, 0], [0, 0, 1], [0, 1, 0]])
ok("the matrix commutes with the swap of its pair, so its eigenlines are fixed by the symmetry alone: channel, sum, difference",
   M*swap == swap*M)
def perm(p):
    return sp.Matrix(3, 3, lambda r, c: 1 if p[c] == r else 0)
def invariants(p, steps):
    out = {}
    for val, mult, vecs in M.eigenvects():
        v = vecs[0]/vecs[0].norm(); prod, g = sp.eye(3), sp.eye(3)
        for _ in range(steps):
            prod = prod*(g*v*v.T*g.T); g = perm(p)*g
        out[str(val)] = sp.simplify(prod.trace())
    return out
hexa, square = invariants(act(third), 3), invariants((0, 2, 1), 2)
ok("hexagonal torus: its own symmetry moves the distinguished channel through all three places; closed comparison 0, +1/8, -1/8",
   hexa == {"1": 0, "2/3": sp.Rational(1, 8), "1/3": sp.Rational(-1, 8)}, {k: str(v) for k, v in hexa.items()})
sq = {}
for val, mult, vecs in M.eigenvects():
    v = vecs[0]/vecs[0].norm(); g = perm((1, 0, 2))
    sq[str(val)] = sp.simplify(((v*v.T)*(g*v*v.T*g.T)).trace())
ok("square torus: its symmetry only swaps; the closed comparison of two places is 0, 1/4, 1/4 and carries no sign",
   sq == {"1": 0, "2/3": sp.Rational(1, 4), "1/3": sp.Rational(1, 4)}, {k: str(v) for k, v in sq.items()})
ok("order of the symmetry on the torus itself: six at the hexagon (three on the marks), four at the square (two on the marks)",
   sp.Matrix(third)**3 == sp.eye(2) and (-sp.Matrix(third))**6 == sp.eye(2) and (-sp.Matrix(third))**3 == -sp.eye(2)
   and sp.Matrix(quarter)**2 == -sp.eye(2))

print(json.dumps({"status": "PAIRINGS_CHECKED", "physical_solution": False, "checks": checks,
  "verdict": "the three steps are the three pairings; the cyclic displacement is produced by the torus itself at the hexagonal point, not at the square point"}, indent=1, default=str))
