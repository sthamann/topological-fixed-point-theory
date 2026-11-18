"""Exact algebraic checks accompanying the analytic source compatibility proofs.

The symbolic polynomial check covers arbitrary charges and real a. Local net
membership, convergence and clustering are proved in the companion texts;
the finite matrices below are not substitutes for those arguments.
"""
import argparse
import json
from pathlib import Path
import sympy as s

parser = argparse.ArgumentParser()
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
checks = {}

def require(name, condition):
    if not condition:
        raise RuntimeError(name)
    checks[name] = True

a = s.symbols('a', real=True)
p = s.Matrix(s.symbols('p0:8', real=True))
delta = s.Matrix([(j+a)/8-s.Rational(1, 2) for j in range(8)])
star = s.Matrix([s.Rational(1, 2)]*4+[-s.Rational(1, 2)]*4)
q = s.ones(8, 1)
A = s.zeros(8)
for j in range(8):
    A[j, (j+4) % 8] = 1

def energy(v):
    return s.expand(8*(v.dot(v)/2+delta.dot(v)))

Tp = A*p+star
require('A_orthogonal_involution', A.T*A == s.eye(8) and A*A == s.eye(8))
require('affine_delta_identity', A*delta == delta+star)
require('star_reverses', A*star == -star)
require('T_involution_all_charges', A*Tp+star == p)
require('H_polynomial_identity_all_charges', s.simplify(energy(Tp)-energy(p)) == 0)
require('Q_polynomial_identity_all_charges', s.simplify(q.dot(Tp-p)) == 0)
require('source_star_zero', energy(star) == 0)
require('conformal_star_weight_one', star.dot(star)/2 == 1)
e0, e7 = s.eye(8)[:, 0], s.eye(8)[:, 7]
require('vector_minimum_candidates', energy(e0) == a and energy(-e7) == 1-a)
require('spinor_minimum_candidates', energy(A*e0+star) == a and energy(-A*e7+star) == 1-a)
require('Cartan_discriminator', e0[0] == 1 and (A*e0+star)[0] == s.Rational(1, 2))
require('L0_not_preserved', (Tp.dot(Tp)-p.dot(p)).expand() != 0)
u = s.Matrix(s.symbols('u0:8'))
require('refined_fugacity_transport', s.expand(u.dot(Tp)-u.dot(star)-(A*u).dot(p)) == 0)
require('twice_root_marking_in_s', int(sum(q/2)) % 2 == 0 and (q/2).dot(q/2) == 2)

# These matrices describe only the already proved two-dimensional H kernel.
Omega = s.Matrix([1, 0])
psi = s.Matrix([0, 1])
Cartan = s.diag(0, 2)
X = s.Matrix([[s.sqrt(2)/2, s.I*s.sqrt(2)/2],
              [s.I*s.sqrt(2)/2, s.sqrt(2)/2]])
require('doublet_unitary_control', X.H*X == s.eye(2))
connected = (Omega.H*X.H*X*Omega)[0] - s.Abs((Omega.H*X*Omega)[0])**2
require('constant_connected_control', connected == s.Rational(1, 2))
require('ward_distinguishes_without_transition',
        (psi.H*Cartan*psi)[0] == 2 and (psi.H*Cartan*Omega)[0] == 0)

result = {'verdict': 'EXACT_ALGEBRA_PASS', 'checks': checks,
          'scope': 'Symbolic H,Q identity and ground controls. Full lattice minima, bounded local witness, simple-current locality and P1 typing are established or audited in the proof texts, not inferred from finite checks.'}
args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
print(json.dumps({'verdict': result['verdict'], 'checks': len(checks)}))
