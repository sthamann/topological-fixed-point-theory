"""Neutral plane stress on the actual fixed three-region E8 word."""
from itertools import combinations

import sympy as sp

z0, z1, z2, z_stress = sp.symbols("z0 z1 z2 z_stress")
positions = (z0, z1, z2)
form = 1 / ((z0 - z1) * (z1 - z2) * (z0 - z2))
ward = sum(1 / (z_stress - z)**2
           + sp.diff(form, z) / form / (z_stress - z) for z in positions)
pair_expression = sum(
    (left - right)**2 / (z_stress - left)**2 / (z_stress - right)**2
    for left, right in combinations(positions, 2)
) / 2
assert sp.cancel(ward - pair_expression) == 0
points = {z0: 1, z1: sp.I, z2: -1, z_stress: -sp.I}
normalized_response = sp.simplify(ward.subs(points))
assert normalized_response == sp.Rational(3, 4)
forms = [-sp.I / 4, sp.I / 4, -sp.I / 4]
marked_stress = [sp.simplify(value * normalized_response) for value in forms]
assert marked_stress == [-3 * sp.I / 16, 3 * sp.I / 16, -3 * sp.I / 16]
print("same original marks 1,i,-1,-i; normalized plane Ward response", normalized_response)
print("0,1,2 inner half-phase insertions", forms)
print("same marked forms with the neutral plane stress", marked_stress)
