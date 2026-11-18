"""Exact checks of the original marked current probe under spatial motion."""
from fractions import Fraction
from itertools import product

import sympy as s

# Two-point function of the weight-four isotropic Wick word on the circle.
theta = s.symbols("theta", real=True)
kernel = s.sin(theta/2)**-8
response = kernel/kernel.subs(theta, s.pi/2)
assert response.subs(theta, s.pi) == s.Rational(1,16)
assert -response.subs(theta, s.pi)/5 == -s.Rational(1,80)
log_response = s.log(kernel)
assert s.diff(log_response, theta).subs(theta, s.pi/2) == -4
assert s.diff(log_response, theta, 2).subs(theta, s.pi/2) == 4
assert s.diff(response, theta, 2).subs(theta, s.pi/2) == 20
print("Move second original quarter point from i to -1: A=1/16, B=-1/80 relative to the initial unclocked value.")
print("Moving normalization: B/A=-1/5; common log slope -4 and log curvature 4 at the initial quarter separation.")

# The two Ward insertions are compared in the same plane coordinates.
x,y,w = s.symbols("x y w")
F = (x-y)**-8
ward = (4/(w-x)**2+4/(w-y)**2)*F + s.diff(F,x)/(w-x) + s.diff(F,y)/(w-y)
normalized = 4*(x-y)**2/((w-x)**2*(w-y)**2)
assert s.cancel(ward/F-normalized) == 0
values = [s.simplify(normalized.subs({x:1,y:z,w:-s.I})) for z in (s.I,-1)]
assert values == [1,4]
print("Plane normalized stress response at w=-i: 1 and 4 before and after movement, identical for A and B.")

# Project B s_4 onto the orbits of the mark-preserving clock.
# The five states in its orbit are different simplex vertices.
orbit_gram = [[Fraction(1) if i==j else Fraction(-1,5) for j in range(5)] for i in range(5)]
# Coefficients of the cyclic autocorrelation diagonalized by the five DFT modes.
# For a nonzero Fourier index the sum of all fifth roots vanishes.
unweighted = sum(orbit_gram[0])/5
nontrivial = (orbit_gram[0][0]-orbit_gram[0][1])/5
assert unweighted == Fraction(1,25)
assert nontrivial == Fraction(6,25)
assert unweighted + 4*nontrivial == 1
print("Marker-moving state: squared norm 1/25 in the invariant line and 6/25 in each of four nontrivial characters.")

# Real representation check: a squared generator permutes all four
# nontrivial characters in one orbit.
nontrivial_orbit = tuple(pow(2,i,5) for i in range(4))
assert sorted(nontrivial_orbit)==[1,2,3,4]
for multiplicity in product(range(4), repeat=5):
    if sum(multiplicity) <= 3 and all(multiplicity[k] == multiplicity[(2*k)%5] for k in range(5)):
        assert not any(multiplicity[1:])
print("The quarter conjugation requires a full orbit of four nontrivial momentum characters; a three-dimensional space can only carry trivial C5.")
