"""Constant R081 singular pair: proposed transported diagnostic clock.

No claim identifies either the original P1 clock or its reflection with these
operators. Truncations check analytic formulae; they do not prove tail bounds.
"""
import json
from pathlib import Path

import numpy as np

rng = np.random.default_rng(252)
yu = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
mass = .3 * np.linalg.svd(yu, compute_uv=False)[0]
ell, z = .2, .3
I = np.eye(2)
P = np.diag([1., 0.])
Q = I - P
Z = np.diag([1., -1.])
X = np.array([[0., 1.], [1., 0.]])
checks = {}

def check(name, condition):
    checks[name] = bool(condition)
    if not condition:
        raise AssertionError(name)

def w(p):
    t = np.arctan(mass / p) / 2
    return np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])

def c(p):
    e = np.hypot(p, mass)
    return .5 * (I + np.tanh(ell * e) / e * (p * Z + mass * X))

def occ(p):
    return .5 * (1 + np.sign(p) * np.tanh(ell * np.hypot(p, mass)))

for alpha in (.25, .75):
    for n in (-20, -1, 0, 20):
        p = n + alpha
        W = w(p)
        h = p * Z + mass * X
        diagonal_h = np.sign(p) * np.hypot(p, mass) * Z
        check(f'rotation_{alpha}_{n}', np.allclose(W.T @ h @ W, diagonal_h, atol=1e-12))
        check(f'covariance_{alpha}_{n}', np.allclose(W.T @ c(p) @ W, np.diag([occ(p), 1-occ(p)]), atol=1e-12))
        check(f'quarter_reflection_{alpha}_{n}', np.allclose(X @ W @ X, w(-p), atol=1e-12))
    N = 5000
    ns = np.arange(-N, N + 1, dtype=float)
    p = ns + alpha
    o, op = occ(p), occ(p - 1)
    trace_norm = np.sum(np.abs(op - o))
    determinant = (1 - z*occ(-N-1+alpha))/(1 - z*occ(N+alpha))
    check(f'trace_norm_telescope_{alpha}', abs(trace_norm-1) < 1e-12)
    check(f'determinant_telescope_{alpha}', abs(determinant - 1/(1-z)) < 1e-12)

# The n,n block is W_n Q W_n*, and n,n-1 is W_n P W_{n-1}*.
# An x-space multiplication operator has the same n,n block for every n.
q0, q1 = w(.25) @ Q @ w(.25).T, w(1.25) @ Q @ w(1.25).T
check('nonconstant_fourier_diagonal', np.linalg.norm(q0-q1) > .01)
check('offdiagonal_exact', abs(q0[0,1] + mass/(2*np.hypot(.25,mass))) < 1e-12)

summaries = []
for N in (100, 1000, 10000):
    ns = np.arange(-N, N+1)
    # One offdiagonal diagonal matrix element is a trace-norm lower bound.
    diag_lb = sum(abs((w(n+.25) @ Q @ w(n+.25).T)[0,1]) for n in ns)
    hs2 = sum(np.linalg.norm(w(n+.25) @ Q @ w(n+.25).T-Q, 'fro')**2
              + np.linalg.norm(w(n+.25) @ P @ w(n-1+.25).T-P, 'fro')**2
              for n in ns)
    summaries.append({'N': N, 'trace_norm_diagonal_lower_bound': diag_lb,
                      'HS_norm': np.sqrt(hs2)})
check('diagonal_lower_bound_grows', summaries[0]['trace_norm_diagonal_lower_bound'] < summaries[1]['trace_norm_diagonal_lower_bound'] < summaries[2]['trace_norm_diagonal_lower_bound'])
check('HS_sums_bounded_in_control', summaries[-1]['HS_norm'] < 1.1*summaries[-2]['HS_norm'])

out = {'scope': 'R081 constant 2x2 singular pair; transported conditional diagnostic shift',
       'mass': mass, 'ell': ell, 'z': z, 'quarter_phases': [.25, .75],
       'fixed_shift_claimed_original': False, 'count_checks': len(checks),
       'checks': checks, 'locality_diagonal_difference': float(np.linalg.norm(q0-q1)),
       'truncation_samples': summaries, 'relative_determinant': 1/(1-z)}
Path(__file__).with_name('COMOVING-CLOCK-CONTROL.json').write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps(out, indent=2))
