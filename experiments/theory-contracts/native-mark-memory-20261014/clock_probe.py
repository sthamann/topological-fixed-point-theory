#!/usr/bin/env python3
"""Exact small-matrix check of the original C/J action on a fixed P5."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import sympy as sp

ROOT = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
CERT_REL = 'experiments/theory-contracts/compiler-integral-triality-20260918/certificate.json'
CERT_PATH = ROOT / CERT_REL
CERT_SHA256 = hashlib.sha256(CERT_PATH.read_bytes()).hexdigest()
cert = json.loads(CERT_PATH.read_text(encoding='utf-8'))
C = sp.Matrix(cert['clock_matrices']['C']['vector'])
J = sp.Matrix(cert['clock_matrices']['J']['vector'])
P = sp.diag(1, 1, 1, 1, 1, 0, 0, 0)
I = sp.eye(8)

# Row-major vectorization. vec(XP-PX) = (P.T kron I - I kron P) vec(X).
def constraints(projectors: list[sp.Matrix]) -> sp.Matrix:
    blocks = [(sp.kronecker_product(Q.T, I) - sp.kronecker_product(I, Q)) for Q in projectors]
    return sp.Matrix.vstack(*blocks)


def rank_mod(matrix: sp.Matrix, prime: int) -> int:
    """Rank over F_p after reducing rational entries (all denominators invertible)."""
    rows, cols = matrix.shape
    a = []
    for row in range(rows):
        vals = []
        for col in range(cols):
            x = matrix[row, col]
            den = int(x.q) % prime
            if den == 0:
                raise ValueError('denominator vanishes modulo prime')
            vals.append((int(x.p) * pow(den, -1, prime)) % prime)
        a.append(vals)
    r = 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][c], -1, prime)
        a[r] = [(v * inv) % prime for v in a[r]]
        for i in range(r + 1, rows):
            fac = a[i][c]
            if fac:
                a[i] = [(x - fac * y) % prime for x, y in zip(a[i], a[r])]
        r += 1
        if r == rows:
            break
    return r


projectors = [P, C * P * C.T, J * P * J.T, C**2 * P * (C**2).T]
checks = {
    'C_orthogonal': C.T * C == I,
    'J_orthogonal': J.T * J == I,
    'P_is_rank5_orthogonal_projector': P * P == P and P.T == P and sp.trace(P) == 5,
    'all_four_are_rank5_orthogonal_projectors': all(Q*Q == Q and Q.T == Q and sp.trace(Q) == 5 for Q in projectors),
}
# Sequential joint-commutant constraints, without imposing symmetry on X.
ranks_mod101 = [rank_mod(constraints(projectors[:k]), 101) for k in range(1, 5)]
full = constraints(projectors)
rank_q = full.rank()
null = full.nullspace()
kernel_matrix = sp.Matrix(8, 8, list(null[0])) if len(null) == 1 else sp.zeros(8)
identity_kernel = len(null) == 1 and sp.simplify(kernel_matrix - sp.trace(kernel_matrix) * I / 8) == sp.zeros(8)
checks['joint_commutant_rank_over_Q_is_63'] = rank_q == 63
checks['joint_commutant_kernel_is_exactly_scalars'] = identity_kernel
checks['commutator_norm_matches_95_over_128'] = sp.trace((P*(C*P*C.T) - (C*P*C.T)*P).T * (P*(C*P*C.T) - (C*P*C.T)*P)) == sp.Rational(95,128)

result = {
    'certificate_path': CERT_REL,
    'certificate_sha256': CERT_SHA256,
    'clock_vectors': {'C': [[str(x) for x in row] for row in C.tolist()], 'J': [[str(x) for x in row] for row in J.tolist()]},
    'projector_definition': 'P5=diag(1,1,1,1,1,0,0,0); P, C P C^T, J P J^T, C^2 P (C^2)^T',
    'unknown_matrix_space': 'all complex 8x8 matrices X; no symmetry constraint imposed',
    'incremental_constraint_ranks_mod_101': ranks_mod101,
    'joint_constraint_rank_over_Q': rank_q,
    'joint_kernel_dimension_over_Q': len(null),
    'joint_kernel_is_scalar_identity': identity_kernel,
    'commutator_hs_norm_squared': str(sp.trace((P*(C*P*C.T) - (C*P*C.T)*P).T * (P*(C*P*C.T) - (C*P*C.T)*P))),
    'checks': checks,
    'all_checks_pass': all(checks.values()),
    'scope': 'The eight-dimensional index representation from the original C/J certificate only. It does not identify these indices with eight states or with physical time, and does not establish that C or J is a symmetry of the Round20–26 Hamiltonian.'
}
out = Path(__file__).with_name('clock-result.json')
out.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
print(json.dumps(result, indent=2))
if not result['all_checks_pass']:
    raise SystemExit(1)
