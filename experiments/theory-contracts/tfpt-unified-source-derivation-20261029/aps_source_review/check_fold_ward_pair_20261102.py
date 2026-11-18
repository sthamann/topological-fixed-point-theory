"""Exact algebraic Ward target: not a primitive P1 correlation calculation."""
import json
from pathlib import Path
import sympy as s

root = Path(__file__).resolve().parents[2]
p = json.loads((root / 'v16-operator-pullback-20261014/fold-results.json').read_text())
M = s.Matrix([[s.Rational(x) for x in row] for row in p['explicit_Weyl_matrix']])
q = s.ones(8, 1)
lam_source = q / 2
lam_target = M * lam_source
expected = s.Matrix([0, 0, 1, 0, 0, 1, 0, 0])
assert M.T * M == s.eye(8) and lam_target == expected
assert (lam_source.T * lam_source)[0] == 2
assert M * q == 2 * expected

z, w, u = s.Integer(0), s.Integer(1), s.Integer(3)
pair = 1 / (w - u)**2
ward_factor = 1 / (z - w) - 1 / (z - u)
source = [s.factor(x * pair * ward_factor) for x in lam_source]
target = [s.factor(x * pair * ward_factor) for x in lam_target]
assert pair == s.Rational(1, 4)
assert source == [-s.Rational(1, 12)] * 8
assert target == [0, 0, -s.Rational(1, 6), 0, 0, -s.Rational(1, 6), 0, 0]
assert (q.T * lam_source)[0] * pair * ward_factor == -s.Rational(2, 3)
assert (2 * lam_target.T * lam_target)[0] * pair * ward_factor == -s.Rational(2, 3)
X = s.Matrix([[0, 1], [1, 0]])
Y = s.Matrix([[0, -s.I], [s.I, 0]])
Z = s.diag(1, -1)
t1, t5 = s.kronecker_product(X, s.eye(2)) / 2, s.kronecker_product(Z, s.eye(2)) / 2
assert s.trace(t1*t1) == s.trace(t5*t5) == 1 and s.trace(t1*t5) == 0
assert t1*t5-t5*t1 == -s.I*s.kronecker_product(Y, s.eye(2))/2
print('PASS: exact existing fold matrix, root charges and Ward target')
print('source eight:', source, 'target eight:', target, 'total:', -s.Rational(2,3))
print('PASS: two native current directions are independent and do not commute')
print('SCOPE: Ward target uses normalized lattice-vacuum CFT; no P1 source correlation computed')
