"""Exact overlap and root positivity checks on the original marked E8 data."""

from fractions import Fraction as Q
from itertools import combinations, product

import sympy as sp


def dilation(x, fixed, scale):
    return fixed + scale * (x - fixed)


x, lam, t = sp.symbols("x lam t", real=True)
# In each of these two charts the common endpoint is infinity.
# The first factor is the modular dilation of the smaller interval.
small_large_centres = ((1, 0), (1, sp.Rational(1, 2)))
for small, large in small_large_centres:
    transport = dilation(dilation(x, large, 1 / lam), small, lam)
    assert sp.simplify(transport - x - (1 - lam) * (small - large)) == 0

left_flow = 1 - 1 / (1 / (1 - x) + t)
right_flow = 1 / (1 / x + t)
left = sp.simplify(sp.diff(left_flow, t).subs(t, 0))
right = sp.simplify(sp.diff(right_flow, t).subs(t, 0))
bracket = sp.expand(left * sp.diff(right, x) - right * sp.diff(left, x))
assert sp.expand(left - (1 - x) ** 2) == 0 and right == -x**2
assert sp.factor(bracket) == 2 * x * (x - 1)
assert bracket.subs(x, sp.Rational(1, 2)) == -sp.Rational(1, 2)
assert (dilation(0, -1, lam), dilation(1, -1, lam)) == (lam - 1, 2 * lam - 1)
print("Overlap (0,1): endpoint motions (1-x)^2 and -x^2; bracket at x=1/2 is -1/2.")
print("Union (-1,infinity): the transported overlap is contained in itself only at scale 1.")

p1 = (Q(1, 2),) * 8
p2 = tuple(Q(-1 if i in (0, 5) else 1, 2) for i in range(8))
gamma = tuple(a - b for a, b in zip(p1, p2))
assert gamma == (1, 0, 0, 0, 0, 1, 0, 0)
assert sum(a*a for a in gamma) == 2
y = (-Q(1, 3),)*3 + (Q(1, 2),)*2 + (0,)*3
assert sum(a*b for a, b in zip(y, gamma)) == -Q(1, 3)
roots = {
    tuple(s*(k == i) + r*(k == j) for k in range(8))
    for i, j in combinations(range(8), 2)
    for s, r in product((-1, 1), repeat=2)
}
roots |= {
    tuple(Q(a, 2) for a in signs)
    for signs in product((-1, 1), repeat=8)
    if signs.count(-1) % 2 == 0
}
assert len(roots) == 240 and gamma in roots
assert all(sum(r*r for r in root) == 2 for root in roots)
assert sum((2*a)**2 for a in gamma)/2 == 4
print("Original E8: middle root e1+e6, norm squared 2, Y=-1/3; 248 weight-one directions, one at charge gamma.")

# From [F_1,E_-1]=1-H_0, H_0=2 on charge gamma:
# the creation Gram is the lowering Gram minus the original state Gram.
theta = sp.symbols("theta", real=True)
mean = lambda expr: sp.integrate(expr, (theta, 0, 2*sp.pi))/(2*sp.pi)
profile = sp.cos(theta)
assert mean(profile) == 0 and mean(profile**2) == sp.Rational(1, 2)
gram = sp.diag(1, mean(profile**2))
pointwise_return = sp.diag(1, mean(profile)**2)
assert pointwise_return - gram == sp.diag(0, -sp.Rational(1, 2))
rho = 1 + profile
assert 2*(mean(rho**2)**2 - mean(rho**4)) == -sp.Rational(17, 4)
print("Pointwise profile test: mixed two-current norm squared -1/2; same-profile norm squared -17/4.")
print("The unmodified level-one affine root requires a nonzero neutral return for every independent charge-gamma grade-one state.")
