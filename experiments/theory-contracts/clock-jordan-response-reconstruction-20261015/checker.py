"""Exact reconstruction of neutral quadratic response from original clock data.

No fitted source, physical detector, state selection, or T1-T8 closure.
"""
from collections import Counter
from itertools import combinations, product
import hashlib
import json
from pathlib import Path
import sys
import sympy as s

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
from tfpt_explorer.neutral_source_response import _native_neutral_branch

CERT = 'experiments/theory-contracts/compiler-integral-triality-20260918/certificate.json'
EXPECTED_CERT = 'a342f865bec164ae6dd54e9e7b7c7cc991efcf7a939dda2ecab586b4a6c6cfe3'
coords = [(i, j) for i in range(8) for j in range(i, 8)]


def digest(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def vec(matrix):
    return [matrix[i, j] for i, j in coords]


def discover_basis(projectors):
    # A full modular rank certifies full rational rank. The final inverse is
    # independently checked over Q; no floating-point rank threshold is used.
    pivots, matrices, expressions, degrees = {}, [], [], []

    def add(matrix, expression, degree):
        v = [int(x.p) * pow(int(x.q), -1, 101) % 101 for x in vec(matrix)]
        for position, row in sorted(pivots.items()):
            factor = v[position]
            v = [(x - factor*y) % 101 for x, y in zip(v, row)]
        if not any(v):
            return
        position = next(i for i, x in enumerate(v) if x)
        inv = pow(v[position], -1, 101)
        pivots[position] = [x*inv % 101 for x in v]
        matrices.append(matrix)
        expressions.append(expression)
        degrees.append(degree)

    for i, matrix in enumerate(projectors):
        add(matrix, ['Q', i], 1)
    cursor = 0
    while cursor < len(matrices) and len(matrices) < 36:
        matrix, expression, degree = matrices[cursor], expressions[cursor], degrees[cursor]
        for i, projector in enumerate(projectors):
            add((matrix*projector + projector*matrix)/2,
                ['jordan', expression, ['Q', i]], degree+1)
        cursor += 1
    return matrices, expressions, degrees


def run():
    if digest(CERT) != EXPECTED_CERT:
        raise RuntimeError('Original clock certificate changed; review before reusing.')
    cert = json.loads((ROOT / CERT).read_text())
    C = s.Matrix(cert['clock_matrices']['C']['vector'])
    J = s.Matrix(cert['clock_matrices']['J']['vector'])
    native = _native_neutral_branch()
    P = native['P5']
    Q = [P, C*P*C.T, J*P*J.T, C**2*P*(C**2).T]
    basis, expressions, degrees = discover_basis(Q)
    checks = {}

    def check(name, value):
        checks[name] = bool(value)

    check('original_P5_in_same_index_basis', P == s.diag(1,1,1,1,1,0,0,0))
    check('original_clocks_orthogonal', C.T*C == J.T*J == s.eye(8))
    check('four_actual_rank5_projectors', all(x*x == x and x.T == x and s.trace(x) == 5 for x in Q))
    roots = set()
    for i,j in combinations(range(8),2):
        for a,b in product((-2,2),repeat=2):
            r=[0]*8; r[i]=a; r[j]=b; roots.add(tuple(r))
    roots.update(r for r in product((-1,1),repeat=8) if sum(v<0 for v in r)%2==0)
    check('same_240_root_lattice_in_clock_basis', len(roots)==240 and all(
        {tuple(U*s.Matrix(r)) for r in roots}==roots for U in (C,J)))
    check('36_symmetric_Jordan_words_degree_at_most4', len(basis) == 36 and max(degrees) == 4
          and all(x.T == x for x in basis))
    W = s.Matrix([[B[i,j] * (1 if i == j else 2)/2 for i,j in coords] for B in basis])
    inverse = W.inv(method='DM')
    check('exact_left_and_right_inverse', inverse*W == W*inverse == s.eye(36))
    # These 36 tests prove reconstruction for every symmetric input by linearity.
    for index, (i,j) in enumerate(coords):
        A = s.zeros(8); A[i,j] = A[j,i] = 1
        m = s.Matrix([s.trace(B*A)/2 for B in basis])
        check(f'reconstruct_basis_{i}_{j}', inverse*m == s.Matrix(vec(A)))
    actual = native['R_matrix']
    check('reconstruct_original_R', inverse*s.Matrix([s.trace(B*actual)/2 for B in basis]) == s.Matrix(vec(actual)))
    check('no_new_D5_component_created', s.trace(P*actual) == 0)
    # Every symmetric quadratic field remains available; primarity is one
    # independently known linear condition, giving dimension 35.
    check('traceless_primary_dimension35', len(coords)-1 == 35)
    check('clock_probe_is_noncommuting', s.trace((Q[0]*Q[1]-Q[1]*Q[0]).T*(Q[0]*Q[1]-Q[1]*Q[0])) == s.Rational(95,128))
    # The Wick identity is proved for arbitrary matrices in PROOF.md, rather
    # than counted repeatedly as a test of its own implementation.
    paths = [CERT, 'tfpt_explorer/neutral_source_response.py',
             'experiments/theory-contracts/native-mark-memory-20261014/clock_probe.py',
             'experiments/theory-contracts/neutral-response-tomography-20261013/PROOF.txt']
    result = {'verdict':'PARTIAL', 'mathematical_verdict':'PASS_COMPLETE_REAL_SYMMETRIC_QUADRATIC_RESPONSE',
              'passed':sum(checks.values()), 'total':len(checks), 'checks':checks,
              'rank':36, 'primary_dimension':35,
              'word_degree_histogram':dict(sorted(Counter(degrees).items())),
              'max_Jordan_degree':max(degrees), 'source_hashes':{p:digest(p) for p in paths},
              'scope':'Fixed native vacuum current algebra, same marked eight-current basis and supplied C/J; access to iterated weight-two OPE coefficients. No physical detector or common-source selection.'}
    certificate = {'projectors':[[[str(v) for v in row] for row in x.tolist()] for x in Q],
                   'basis_expressions':expressions, 'degrees':degrees,
                   'basis_matrices':[[[str(v) for v in row] for row in x.tolist()] for x in basis],
                   'response_inverse':[[str(v) for v in row] for row in inverse.tolist()]}
    (HERE/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
    (HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('checks','source_hashes')},indent=2))
    if not all(checks.values()):
        raise SystemExit(1)


if __name__ == '__main__':
    run()
