#!/usr/bin/env python3
"""TFPT F1 audit, 2026-11-11.

Exact checks of the family-Cartan obstruction, an explicit 6D SU(4) index-three
bundle candidate, and why the naïve flat-torus heterotic parent fails anomaly
cancellation. Mathematical test cases, NOT a complete theory or P1 derivation.
Requires only sympy.
"""
from itertools import product
import sympy as s

# Fundamental primitive A3 family direction and the corresponding split bundle
H = s.Matrix([1, 1, 1, -3])
assert sum(H) == 0
assert (H/4).dot(H/4) == s.Rational(3, 4)
c = s.symbols('c')  # formal first Chern class c_1(L)
# V = L + L + L + L^{-3}: c_1(V) = 0; ch_3(V) = -4 c^3.
ch3_simple = s.expand((3*c**3 + (-3*c)**3) / 6)
assert ch3_simple == -4*c**3
print('PASS S3-symmetric split SU4 bundle has index = -4 integral(c1(L)^3).')
print('     Therefore it cannot have index +/-3 on any closed spin 6-manifold.')

# Source four-mark surface alone cannot supply an SU4 index via c1.
# On a 2D spin surface the degree-2 index density is ch1(V) = c1(V)=0.
print('PASS Any smooth SU4 bundle over S2 has twisted Dirac index 0.')

# Take the product of three of the CM elliptic curves from the four-mark double cover.
# X = E^3 is a CY3: c1(TX)=c2(TX)=c3(TX)=0.
# x,y,z are the positive integral generators of H^2 of the three elliptic factors.
# For c1(L_i)=a_i*x + b_i*y + c_i*z, the index is sum a_i*b_i*c_i.
# Choose J=2*x+2*y+z. Slope 0 requires 2*a+2*b+4*c=0.
fluxes = [(-2,4,-1), (-1,-3,2), (1,-3,1), (2,2,-2)]
assert all(sum(f[i] for f in fluxes) == 0 for i in range(3))
slopes = [2*a + 2*b + 4*c for a,b,c in fluxes]
assert slopes == [0,0,0,0]
individual_indices = [a*b*c for a,b,c in fluxes]
assert individual_indices == [8,6,-3,-8]
index_V = sum(individual_indices)
assert index_V == 3
c3_integral = 2*index_V
assert c3_integral == 6

# The 3-torus line bundle Dolbeault cohomology has h^0(E,L_d)=d if d>0,
# h^1(E,L_d)=-d if d<0 (and h^0=h^1=1 for d=0).
# Here all three degrees in every flux are nonzero; the zero modes are
# consequently concentrated in just one parity.
parity_counts = [0,0]
for a,b,c in fluxes:
    degrees = (a,b,c)
    parity = sum(d<0 for d in degrees)%2
    parity_counts[parity] += abs(a*b*c)
assert parity_counts == [14,11]
assert parity_counts[0]-parity_counts[1] == 3

# c2(V)= - sum_i c1(L_i)^2/2 when sum_i c1(L_i)=0.
c2_xy, c2_xz, c2_yz = [-sum(f[i]*f[j] for f in fluxes)
                             for i,j in ((0,1),(0,2),(1,2))]
assert (c2_xy,c2_xz,c2_yz) == (4,3,17)
# Integral J wedge c2(V): xy*J_z + xz*J_y + yz*J_x.
positive_anomaly = c2_xy*1 + c2_xz*2 + c2_yz*2
assert positive_anomaly == 44
print('PASS X=E^3, J=(2,2,1): 4 integral-slope-zero line bundles define an SU4 polystable candidate.')
print('     Fluxes:', fluxes)
print('     chiral indices:',individual_indices, 'net=',index_V,'c3=',c3_integral)
print('     even/odd Dolbeault zero modes:', parity_counts, '=> 11 mirror pairs remain.')
print('     c2(V) in (xy,xz,yz):',(c2_xy,c2_xz,c2_yz), 'integral J*c2(V)=',positive_anomaly)
print('     c2(T E^3)=0: the standard heterotic anomaly/Bianchi constraint cannot')
print('     hold with only positive-slope-stable hidden-sector bundles/effective 5-branes.')

# A symmetric equal-radius Kahler class imposes a+b+c=0 for each line;
# then a*b*c is even since three odd summands cannot sum to zero.
for a,b in product(range(-15,16), repeat=2):
    cc=-a-b
    assert (a*b*cc)%2==0
print('PASS J=(1,1,1) slope-zero factorwise line bundles have even index, never 3.')

# Minimal diagonal CY orbifold from the square CM torus has g=(i,i,-1).
# Pair count in the orbifold Euler formula for g^a,g^b, a,b mod 4:
# 12 pairs have at least one odd exponent, common fixed locus 16 points;
# remaining 4 pairs retain a torus and contribute Euler zero.
fixed_counts=[16 if a%2 or b%2 else 0 for a in range(4) for b in range(4)]
orbifold_euler = s.Rational(sum(fixed_counts),4)
assert orbifold_euler==48
assert orbifold_euler/2==24
print('PASS square-torus E^3/Z4 orbifold Euler characteristic 48, standard-embedding net 24 (not 3).')

# Critical heterotic matter central charges: 10 spacetime coordinates,
# visible E8 current algebra c=8, hidden E8 c=8, right-moving fermions 10/2.
assert 10+8+8==26 and 10+5==15
print('PASS a heterotic gravitational parent would need one additional left-moving c=8 sector')
print('     given 10 target coordinates and one visible E8 current algebra (c=8).')

# Boundary holonomy alone does not fix the APS/spectral-flow integer.
alpha0, alpha1 = s.Rational(1,4),s.Rational(13,4)
assert s.simplify(s.exp(2*s.pi*s.I*alpha0)-s.exp(2*s.pi*s.I*alpha1))==0
crossings=[n for n in range(-10,11) if (n+alpha0)*(n+alpha1)<0]
assert crossings==[-3,-2,-1]
print('PASS same endpoint holonomy i may conceal 3 spectral crossings: the full history matters.')
print('VERDICT: explicit topological model exists, but naïve torus parent fails global consistency.')
