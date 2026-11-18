"""Exact implication test, not a new TFPT field model or physical evaluation."""
import json
import sympy as s

a, b = s.symbols('a b', positive=True, real=True)
x = s.Matrix([[0, 1], [1, 0]])
z = s.diag(1, -1)
I = s.eye(2)
tau = s.kronecker_product(z, I)
grading = s.kronecker_product(I, z)
D = s.kronecker_product(a*I + b*x, x)
reference = s.kronecker_product(a*I, x)
assert (D + tau*D*tau)/2 == reference
assert D*grading + grading*D == s.zeros(4)
assert reference*grading + grading*reference == s.zeros(4)
ratio = s.factor(D.det()/reference.det())
assert s.simplify(ratio - (1-b*b/(a*a))**2) == 0

q = s.symbols('q')
relative_zeta = 2*((a+b)**(-2*q) + (a-b)**(-2*q) - 2*a**(-2*q))
assert s.residue(relative_zeta, q, 1) == 0
assert relative_zeta.subs(q, 0) == 0
derivative = s.diff(relative_zeta, q).subs(q, 0)
probe_derivative = s.simplify(derivative.subs({a:2,b:1}))
assert s.simplify(s.exp(-probe_derivative/2) - s.Rational(9,16)) == 0

v0,v2,v4,J,g = s.symbols('v0 v2 v4 J g')
# Laurent coefficients of the relative Mellin numerator; no heat data fitted.
numerator = v4/q + (J-v0/s.Integer(2)-v2)
germ = s.expand((q+s.EulerGamma*q*q)*numerator)
assert s.simplify(germ.coeff(q,0)-v4) == 0
assert s.simplify(germ.coeff(q,1)-(J-v0/2-v2+s.EulerGamma*v4)) == 0

print(json.dumps({
 'scope':'Spectral implication test only, not an actual TFPT background or full source',
 'bound_reference':'(D+tau D tau)/2 = a I tensor sigma_x, independent of b',
 'domain':'0 < b < a; same finite domain, grading and no zero crossing',
 'relative_residue_at_1':0,
 'relative_finite_part_at_0':0,
 'relative_eta':0,
 'determinant_ratio':str(ratio),
 'probe':{'a':2,'b':1,'minimum_absolute_eigenvalue':1,
          'residue_exponent_weight':1,
          'relative_determinant_modulus':'9/16',
          'relative_zeta_derivative':str(probe_derivative)},
 'mellin_derivative':'J - v0/2 - v2 + EulerGamma*v4',
 'verdict':'RESIDUE_PAIR_AND_BOUND_REFERENCE_ALONE_DO_NOT_RECOVER_FULL_DETERMINANT'
},indent=2))
