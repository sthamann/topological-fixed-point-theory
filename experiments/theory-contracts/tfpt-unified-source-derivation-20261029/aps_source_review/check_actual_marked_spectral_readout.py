"""Exact full marked reader versus its linear scalar centralization.

No source, metric, clock or state is selected by this calculation.
M0 is the already verified original v118/common Fuchs matrix.
"""
import json
from pathlib import Path
import sympy as s

i = s.I
M = s.Matrix([[0, -(1+i)/2, (1-i)/2],
              [-(1+i)/2, -i/2, -s.Rational(1, 2)],
              [(1-i)/2, -s.Rational(1, 2), i/2]])
P = s.diag(1, 0, 0)
Dsign = s.diag(1, -1, 1)
a, b, y, delta, u, v = s.symbols('alpha beta y delta u v')
F = a*M*P - b*P*M
N = -Dsign*F*Dsign
c = (a-b)/3
centralize = lambda X: s.simplify(sum(
    (M**k*X*M**(2*k) for k in range(3)), s.zeros(3))/3)
full = s.expand((y*s.eye(3)-delta*M-u*F-v*N).det())
expected_full = (y**3-delta**3-(a-b)*delta**2*u
                 -i*(a-b)*delta*y*v+a*b*delta*u**2+2*i*a*b*y*u*v)
central = s.expand((y*s.eye(3)-(delta+c*u)*M-i*c*v*M**2).det())
expected_central = y**3-(delta+c*u)**3-(i*c*v)**3-3*y*(delta+c*u)*(i*c*v)
expected_difference = (delta*u**2*(a*a+a*b+b*b)/3
    +i*y*u*v*(a*a+4*a*b+b*b)/3+c**3*u**3-i*c**3*v**3)
rho = s.Rational(1, 2)-i*s.sqrt(3)/2
C = (s.sqrt(6)-s.sqrt(2))/4-i*(s.sqrt(6)+s.sqrt(2))/4
A, B = s.symbols('A B', real=True)
E = lambda r, q: s.eye(3)[:, r]*s.eye(3)[q, :]
Fc = C*(A*E(0, 2)+B*E(2, 1))
Nc = -i*C*(A*E(1, 2)+B*E(2, 0))
checks = {
    'original_cube': s.simplify(M**3)==s.eye(3),
    'F_centralization': s.simplify(centralize(F)-c*M)==s.zeros(3),
    'N_centralization': s.simplify(centralize(N)-i*c*M**2)==s.zeros(3),
    'full_marked_determinant': s.simplify(full-expected_full)==0,
    'central_determinant': s.simplify(central-expected_central)==0,
    'nonlinear_difference': s.simplify(full-central-expected_difference)==0,
    'linear_agreement': all(s.simplify(s.diff(full-central, z).subs({u:0,v:0}))==0 for z in (u,v)),
    'actual_phase': s.simplify(s.expand_complex(rho*s.conjugate(C)-i*C))==0,
    'full_adjoint_residual': s.simplify(s.expand_complex(
        Nc-rho*Fc.conjugate().T+i*C*(A+B)*(E(1,2)+E(2,0))))==s.zeros(3),
    'quadratic_not_identically_zero_for_real_ratio': s.discriminant(1+s.Symbol('r')+s.Symbol('r')**2,s.Symbol('r'))==-3,
}
result = {'checks':checks,'passed':all(checks.values()),
    'full_determinant':str(expected_full),
    'central_determinant':str(expected_central),
    'full_minus_central':str(expected_difference),
    'scope':'Full original marked coefficient family. Centralization is a linear scalar readout, not a whole nonlinear source replacement. Original real delta contract retained.'}
path=Path(__file__).with_name('ROOT-ACTUAL-MARKED-SPECTRAL-READOUT-CHECK.json')
path.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
if not result['passed']:
    raise SystemExit(1)
