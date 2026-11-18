"""Exact original Q-twin transport; no physical source selection is assumed."""
from pathlib import Path
import sys
import sympy as s

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / 'verification'))
from v568_bit_selector_ladder import Qc, UV

I = s.eye(3)
G = s.Matrix([[1, 0, 0], [0, 1, 0], [0, 2, -1]])
Sigma = s.diag(1, -1, -1)
U0, V0 = UV(Qc(0)); U2, V2 = UV(Qc(2))
if not (G*G == I and G.det() == -1):
    raise AssertionError('control at original line 13')
if not (G*U0 == U2*G and G*V0 == V2*G):
    raise AssertionError('control at original line 14')
if not (G*Sigma == Sigma*G):
    raise AssertionError('control at original line 15')
if not (G*s.Matrix([1, 0, 0]) == s.Matrix([1, 0, 0])):
    raise AssertionError('control at original line 16')
if not (G*s.ones(3, 1) == s.ones(3, 1)):
    raise AssertionError('control at original line 17')
if not (G*s.Matrix([1, 2, 0]) == s.Matrix([1, 2, 4])):
    raise AssertionError('control at original line 18')

# The simultaneous intertwiner space is exactly one-dimensional, not merely
# an example. Its nonzero members cannot preserve the stipulated Sigma form.
xx = s.symbols('x0:9'); X = s.Matrix(3, 3, xx)
solutions = s.linsolve(list(X*U0-U2*X)+list(X*V0-V2*X), xx)
if not (solutions == s.FiniteSet(tuple(-xx[8]*G))):
    raise AssertionError('control at original line 24')
B = G.T*Sigma*G
if not (B == s.Matrix([[1, 0, 0], [0, -5, 2], [0, 2, -1]])):
    raise AssertionError('control at original line 26')
if not (B != Sigma):
    raise AssertionError('control at original line 27')
theta = lambda A: Sigma*A.T*Sigma
theta_transport = lambda A: G*theta(G*A*G)*G
if not (s.trace(theta(V0)*V0) == 4):
    raise AssertionError('control at original line 30')
if not (s.trace(theta(V2)*V2) == 8):
    raise AssertionError('control at original line 31')
if not (s.trace(theta_transport(V2)*V2) == 4):
    raise AssertionError('control at original line 32')
if not (theta_transport(U2) == G*theta(U0)*G):
    raise AssertionError('control at original line 33')

M = s.Matrix([[1, 0, 0], [0, s.Rational(1, 2), s.Rational(1, 6)],
              [0, s.Rational(1, 6), s.Rational(1, 2)]])
c = s.symbols('c', real=True)
Gc = s.Matrix([[1,0,0],[0,s.Rational(1,2),s.Rational(1,2)],
               [0,(c+1)/2,(c-1)/2]])
Qplus = (Qc(c)+Sigma*Qc(c)*Sigma)/2
if not (s.simplify(Gc*M-Qplus*Gc/3) == s.zeros(3)):
    raise AssertionError('control at original line 41')
if not (Gc.det() == -s.Rational(1,2)):
    raise AssertionError('control at original line 42')

print('PASS: complete original pairs are GL3(Z)-conjugate; winding and grading retained.')
print('Fixed reflected word values: 4 and 8. Fully transported reflected value: 4 and 4.')
print('The original three-rate clock admits an intertwiner for both twins.')
print('No physical P1 equivalence, selected source metric, state or normalization is inferred.')
