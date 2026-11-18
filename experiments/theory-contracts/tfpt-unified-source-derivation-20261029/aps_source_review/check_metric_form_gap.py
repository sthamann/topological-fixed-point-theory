"""R082: diagnostic form-bound counterexamples on the existing v311 transfer.

The auxiliary N and V are chosen tests, not the source's geometric operators.
"""
from pathlib import Path
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'verification'))
from v311_gap_decoupled_measure import transfer

T, lam2, lam3 = transfer()
values, vectors = np.linalg.eigh(T)
order = np.argsort(values)[::-1]
values, vectors = values[order], vectors[:, order]
projectors = [np.outer(vectors[:, i], vectors[:, i]) for i in range(3)]
p0, p1, p2 = projectors
H = sum(-np.log(values[i]) * projectors[i] for i in range(3))
I = np.eye(3)
delta = 6 * np.log(1.5)
v = 31 / (8 * np.pi ** 2)

def gap(A):
    eig = np.linalg.eigvalsh(A)
    return float(eig[1] - eig[0])

def form_norm(N, V):
    w, U = np.linalg.eigh(N + I)
    inv = (U * (w ** -0.5)) @ U.T
    return float(np.linalg.norm(inv @ V @ inv, 2))

checks = {}
def check(name, condition):
    checks[name] = bool(condition)
    if not condition:
        raise RuntimeError('failed: ' + name)

check('native_transfer_spectrum', np.allclose(values, [1, lam2, lam3], atol=1e-13))
check('native_H_gap', abs(gap(H) - delta) < 1e-12)
check('native_projectors', np.linalg.norm(sum(projectors) - I) < 1e-12)

# Stronger-than-stated diagnostic assumption N=H, no physical identification.
V = v * p0 - v * (1 + delta) * p1
bound_gap = gap(H + V)
check('same_generator_graph_norm', abs(form_norm(H, V) - v) < 1e-12)
check('same_generator_gap_exact', abs(bound_gap - (delta-v*(delta+2))) < 1e-12)
check('same_generator_old_bound_fails', bound_gap < delta - 2*v)
check('same_generator_gap_positive', bound_gap > 0)

# The ordinary norm bound really does imply the advertised numerical bound.
Vordinary = v * (p0 - p1)
check('ordinary_norm', abs(np.linalg.norm(Vordinary,2) - v) < 1e-12)
check('ordinary_bound_sharp', abs(gap(H + Vordinary) - (delta - 2*v)) < 1e-12)

small_gaps = []
for eps in [0.1, 0.01, 0.001]:
    t = (delta - eps) / v - 1
    N = t * p1
    Veps = -(delta-eps) * p1
    check('positive_N_' + str(eps), np.linalg.eigvalsh(N)[0] > -1e-12)
    check('form_norm_' + str(eps), abs(form_norm(N, Veps) - v) < 1e-12)
    check('actual_gap_' + str(eps), abs(gap(H + Veps) - eps) < 1e-12)
    small_gaps.append({'epsilon': eps, 'actual_gap': gap(H + Veps),
                       'ordinary_norm_V': float(np.linalg.norm(Veps,2)),
                       'relative_norm_V': form_norm(N, Veps)})

result = {'scope': 'v311 transfer representative; auxiliary N and V are diagnostic, not original TFPT geometry',
          'delta': float(delta), 'v': float(v), 'advertised_delta_minus_2v': float(delta-2*v),
          'same_H_form_bound_gap': bound_gap, 'same_H_formula': float(delta-v*(delta+2)),
          'small_gaps': small_gaps, 'checks': checks}
out = Path(__file__).with_name('metric-form-gap-results.json')
out.write_text(json.dumps(result, indent=2, ensure_ascii=False)+'\n')
print(json.dumps(result, indent=2, ensure_ascii=False))
