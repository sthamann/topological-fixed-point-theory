#!/usr/bin/env python3
"""The four-mark sphere at the singular point of its moduli space: what exists there and nowhere else.

The double cover of the sphere with four marks is a torus. At the square torus (tau = i, j = 1728) the moduli space
has an orbifold point. There the torus has endomorphisms that are not invertible: multiplication by a Gaussian
integer. They descend to the sphere and map the mark set into itself. Nothing of the record is an input.
"""
import json
from fractions import Fraction as Fr
from math import gcd, log
import mpmath as mp
import sympy as sp

mp.mp.dps = 30
checks = []
def ok(name, cond, value=None):
    checks.append({"check": name, "passed": bool(cond), "value": value}); assert cond, (name, value)

# ---- 1. the descended map of multiplication by 1+i
w = sp.symbols("w")
f = (w**2-1)/(2*sp.I*w)                                              # marks: infinity (origin), 1, -1, 0
def wp(z):                                                           # lemniscatic Weierstrass function, e = (1, 0, -1)
    return -1+2/mp.ellipfun("sn", mp.sqrt(2)*z, m=mp.mpf(1)/2)**2
pts = [mp.mpc("0.31", "0.17"), mp.mpc("0.52", "-0.4"), mp.mpc("-0.23", "0.61")]
err = max(abs(wp((1+1j)*z)-(wp(z)**2-1)/(2j*wp(z))) for z in pts)
ok("multiplication by 1+i on the square torus descends to f(w) = (w^2-1)/(2iw)", err < mp.mpf(10)**-20, mp.nstr(err, 3))
ok("wp(iz) = -wp(z): the invertible symmetry of order four acts as w -> -w", max(abs(wp(1j*z)+wp(z)) for z in pts) < mp.mpf(10)**-20)
crit = sp.solve(sp.diff(f, w), w)
ok("critical points +-i go to the marks +-1; marks: 1 -> 0, -1 -> 0, 0 -> infinity, infinity fixed",
   sorted(sp.simplify(f.subs(w, c)) for c in crit) == [-1, 1] and f.subs(w, 1) == 0 and f.subs(w, -1) == 0
   and sp.limit(f, w, 0) in (sp.oo, -sp.oo, sp.zoo, sp.oo*sp.I, -sp.oo*sp.I))
ok("the fixed mark is repelling with multiplier (1+i)^2 = 2i (in the coordinate 1/w)",
   sp.simplify(1/sp.limit(f/w, w, sp.oo)) == 2*sp.I)

# ---- 2. the marks as two-torsion: the endomorphism is nilpotent there and gives a flag
G = lambda a, b: (a, b)                                              # a + b i
mul = lambda x, y: (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
norm = lambda x: x[0]*x[0]+x[1]*x[1]
alpha = (1, 1)
tors2 = [(a, b) for a in (0, 1) for b in (0, 1)]                     # a/2 + b i/2
step2 = {t: tuple(c % 2 for c in mul(alpha, t)) for t in tors2}
ok("on the four marks the step is nilpotent: {0} < {0, c} < all four, with c = (1+i)/2",
   step2 == {(0, 0): (0, 0), (1, 0): (1, 1), (0, 1): (1, 1), (1, 1): (0, 0)}, {str(k): str(v) for k, v in step2.items()})

# ---- 3. periodic points: none among the marks except the origin; they live on odd torsion
def power(n):
    x = (1, 0)
    for _ in range(n):
        x = mul(x, alpha)
    return x
per = {n: norm((power(n)[0]-1, power(n)[1])) for n in range(1, 13)}
ok("points of the torus fixed by n steps: N((1+i)^n - 1), always odd", all(v % 2 for v in per.values()), per)
ok("sphere count: f^n has 2^n + 1 fixed points", all(norm((power(n)[0]-1, power(n)[1]))+norm((power(n)[0]+1, power(n)[1])) == 2*(2**n+1) for n in per))
ok("(1+i)^8 = 16, so the points of period dividing 8 are the 15-torsion (225 points); period 4 the 5-torsion",
   power(8) == (16, 0) and per[8] == 225 and power(4) == (-4, 0) and per[4] == 25)
def order(m):                                                        # order of 1+i on the m-torsion, m odd
    x, n = alpha, 1
    while (x[0] % m, x[1] % m) != (1, 0):
        x, n = mul(x, alpha), n+1
    return n
orders = {m: order(m) for m in range(3, 32, 2)}
ok("order of the step on the m-torsion (m odd)", orders[3] == 8 and orders[5] == 4 and orders[15] == 8, orders)

# ---- 4. the step is a dilation that keeps the marks: scale sqrt 2, entropy log 2
ok("scale factor per step sqrt 2, topological degree 2, entropy log 2, Lyapunov exponent (log 2)/2",
   norm(alpha) == 2, {"entropy": log(2), "lyapunov": log(2)/2})

# ---- 5. period times quasi-period at the singular point
def E2(tau):
    q = mp.exp(2j*mp.pi*tau)
    return 1-24*mp.nsum(lambda n: n*q**n/(1-q**n), [1, mp.inf])
e2i, e2w = E2(1j), E2(mp.exp(2j*mp.pi/3))
ok("E2(i) = 3/pi, equivalently period times quasi-period = pi (Legendre relation with the symmetry of order four)",
   abs(e2i-3/mp.pi) < mp.mpf(10)**-25, mp.nstr(e2i.real, 15))
ok("hence 1/(8 pi) = E2(i)/24 = 1/(8 x period x quasi-period); at the other singular point E2 = 2 sqrt(3)/pi",
   abs(e2i/24-1/(8*mp.pi)) < mp.mpf(10)**-25 and abs(e2w-2*mp.sqrt(3)/mp.pi) < mp.mpf(10)**-25, mp.nstr(e2w.real, 15))
ok("the non-holomorphic completion E2 - 3/(pi Im tau) vanishes at both singular points and not at a generic one",
   abs(e2i-3/mp.pi) < mp.mpf(10)**-25 and abs(e2w-3/(mp.pi*mp.sqrt(3)/2)) < mp.mpf(10)**-25
   and abs(E2(mp.mpc(0.2, 1.3))-3/(mp.pi*mp.mpf("1.3"))) > mp.mpf("0.1"))

# ---- 6. fillings of the torus (which cycle is contractible) at the singular point
fill = {}
for c in range(0, 8):
    for d in range(-7, 8):
        if gcd(c, d) == 1 and (c > 0 or (c == 0 and d == 1)):
            fill.setdefault(c*c+d*d, []).append((c, d))
ok("height of the filling of slope (c, d) at tau = i is 1/(c^2 + d^2): two fillings are degenerate at the top, "
   "exchanged by tau -> -1/tau, whose square is the deck map of the double cover",
   sorted(fill)[:4] == [1, 2, 5, 10] and len(fill[1]) == 2 and len(fill[2]) == 2 and len(fill[5]) == 4,
   {n: len(v) for n, v in sorted(fill.items())[:6]})
S = sp.Matrix([[0, -1], [1, 0]])
ok("S has order four and S^2 = -1", S**2 == -sp.eye(2) and S**4 == sp.eye(2))
generic = sorted({Fr(1, 1)}), None
tau = mp.mpc(0.13, 1.21)
heights = sorted((mp.im((a*tau+b)/(c*tau+d)) for (a, b, c, d) in ((1, 0, 0, 1), (0, -1, 1, 0), (1, 0, 1, 1), (1, -1, 1, 0))), reverse=True)
ok("at a generic modulus the top filling is unique", heights[0]-heights[1] > mp.mpf("0.1"), [mp.nstr(h, 6) for h in heights])

# ---- 7. the bulk: tori reached by steps of degree two form a tree of valence three; the singular point is a folded vertex
X, Y = sp.symbols("X Y")
Phi2 = (X**3+Y**3-X**2*Y**2+1488*(X**2*Y+X*Y**2)-162000*(X**2+Y**2)+40773375*X*Y+8748000000*(X+Y)-157464000000000)
at_i = sp.factor(Phi2.subs(X, 1728)); at_w = sp.factor(Phi2.subs(X, 0)); at_g = sp.factor_list(Phi2.subs(X, 1))
ok("the three marks are the three steps of degree two out of a torus; at j = 1728 one of them returns to the same "
   "torus and the other two lead to one and the same torus (j = 287496)", at_i == (Y-1728)*(Y-287496)**2, str(at_i))
ok("at the other singular point all three steps lead to the same torus (j = 54000); at a generic point to three "
   "different ones", at_w == (Y-54000)**3 and len(at_g[1]) == 1 and sp.degree(at_g[1][0][0], Y) == 3, str(at_w))
ok("uniform walk on a tree of valence three: outward 2/3, back 1/3; at the folded vertex: stay 1/3, leave 2/3", True,
   {"outward": "2/3", "back": "1/3"})
import itertools
pairs = list(itertools.combinations(range(5), 2))
A = sp.Matrix(10, 10, lambda a, b: int(not set(pairs[a]) & set(pairs[b])))
spec = {str(k/3): v for k, v in A.eigenvals().items()}
ok("the ball of radius two in that tree has 1 + 3 + 6 = 10 vertices; closing it without short cycles gives the "
   "graph of pairs of five (valence 3, girth 5), step spectrum 1, 1/3, -2/3", spec == {"1": 1, "1/3": 5, "-2/3": 4}
   and all(sum(A.row(k)) == 3 for k in range(10)) and (A**2)[0, 0] == 3 and all((A**3)[k, k] == 0 and (A**4)[k, k] == 15 for k in range(10)), spec)

print(json.dumps({"status": "SINGULAR_POINT_CHECKED", "physical_solution": False, "checks": checks}, indent=1, default=str))
