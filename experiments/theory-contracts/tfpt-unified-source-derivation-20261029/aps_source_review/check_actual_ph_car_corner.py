"""Independent CAR replay of the actual finite E/Eh source summand."""
from pathlib import Path
import json
import sympy as s

n, size = 6, 64
ann = []
for j in range(n):
    entries = {}
    for k in range(size):
        if k & (1 << j):
            entries[(k ^ (1 << j), k)] = (-1) ** (
                (k & ((1 << j) - 1)).bit_count())
    ann.append(s.SparseMatrix(size, size, entries))


def lift(matrix):
    out = s.zeros(size)
    for i in range(n):
        for j in range(n):
            if matrix[i, j] != 0:
                out += matrix[i, j] * (ann[i].T * ann[j])
    return out


def unit(i, j):
    return s.eye(3)[:, i] * s.eye(3)[j, :]


def zero(value):
    return s.simplify(s.expand_complex(value)) == 0


a, b = s.symbols('a b', real=True, nonzero=True)
sech, tanh = s.symbols('sech tanh', real=True)
phase = (s.sqrt(6) - s.sqrt(2)
         - s.I * (s.sqrt(6) + s.sqrt(2))) / 4
forward = phase * (a * unit(0, 2) + b * unit(2, 1))
reverse = -s.I * phase * (a * unit(1, 2) + b * unit(2, 0))
face = s.Matrix([[sech, tanh], [tanh, -sech]])
of = lift(s.kronecker_product(face, forward))
on = lift(s.kronecker_product(face, reverse))
omega = s.zeros(size, 1)
omega[7] = 1  # The three original occupied modes; complements are empty.


def expectation(operator):
    return s.simplify((omega.T * operator * omega)[0])


checks = {
    'exact_CAR': all(
        ann[i] * ann[j].T + ann[j].T * ann[i]
        == (s.eye(size) if i == j else s.zeros(size))
        for i in range(n) for j in range(n)),
    'filled_defect_Pauli':
        lift(s.diag(forward, s.zeros(3))) * omega == s.zeros(size, 1),
    'FF_unadjointed': zero(expectation(of * of)),
    'NN_unadjointed': zero(expectation(on * on)),
    'FN_connected': zero(expectation(of * on)
                         - tanh ** 2 * s.trace(forward * reverse)),
    'NF_connected': zero(expectation(on * of)
                         - tanh ** 2 * s.trace(reverse * forward)),
    'positive_F_adjoint': zero(expectation(of.conjugate().T * of)
                               - tanh ** 2 * (a * a + b * b)),
    'positive_N_adjoint': zero(expectation(on.conjugate().T * on)
                               - tanh ** 2 * (a * a + b * b)),
    'positive_FN_cross_orthogonal':
        zero(expectation(of.conjugate().T * on)),
}
quadratures = [
    (of + of.conjugate().T) / s.sqrt(2),
    (of - of.conjugate().T) / (s.I * s.sqrt(2)),
    (on + on.conjugate().T) / s.sqrt(2),
    (on - on.conjugate().T) / (s.I * s.sqrt(2)),
]
covariance = s.Matrix(4, 4, lambda i, j: expectation(
    (quadratures[i] * quadratures[j]
     + quadratures[j] * quadratures[i]) / 2))
z = s.simplify(s.trace(forward * reverse))
h = a * a + b * b
zr, zi = s.re(s.expand_complex(z)), s.im(s.expand_complex(z))
real_cov = s.Matrix([[h, 0, zr, zi], [0, h, zi, -zr],
                     [zr, zi, h, 0], [zi, -zr, 0, h]])
checks['four_selfadjoint_quadratures'] = all(
    zero(covariance[i, j] - tanh ** 2 * real_cov[i, j])
    for i in range(4) for j in range(4))
q = s.symbols('q')
checks['positive_covariance_polynomial'] = s.factor(
    (covariance - q * s.eye(4)).det()
    - ((tanh ** 2 * h - q) ** 2 - tanh ** 4 * 4 * a * a * b * b) ** 2
) == 0
out = {
    'scope': 'Actual canonical finite E/Eh summand, one Pgamma+ carrier. '
             'Ordinary full retained trace multiplies by16 and (2/pi)^2. '
             'Not the whole unsmoothed Z response or a full JF-real observable.',
    'checks': checks, 'passed': all(checks.values()),
    'directed_FN_trace': str(z),
    'hermitian_covariance_without_common_factor': str(real_cov),
    'covariance_eigenvalues_without_common_factor': [
        '(abs(a)+abs(b))**2, multiplicity2',
        '(abs(a)-abs(b))**2, multiplicity2'],
    'finite_CAR_dimension': size,
    'face_parameters': 'symbolic sech and tanh; all expectations verified '
                       'before substituting ell_NS',
}
if not out['passed']:
    raise RuntimeError('Independent CAR replay failed')
Path(__file__).with_name('ROOT-ACTUAL-PH-CAR-AND-ADJOINT-REPLAY.json').write_text(
    json.dumps(out, indent=2) + '\n')
print(json.dumps(out))
