"""Use the reconstructed R063 control law to predict an independent paired quarter word."""

import ast
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import sys

import numpy as np


root = Path(sys.argv[1]).resolve()
source = (root / 'experiments/theory-contracts/tfpt-unified-source-derivation-20261029'
          / 'aps_source_review/check_quarter_tower_higgs_inverse.py')
tree = ast.parse(source.read_text())
assert sum(isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)
           and isinstance(node.value.func, ast.Attribute) and node.value.func.attr == 'write_text'
           for node in tree.body) == 1
tree.body = [node for node in tree.body if not (isinstance(node, ast.Expr)
             and isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Attribute)
             and node.value.func.attr == 'write_text')]
old = {'__file__': str(source)}
with redirect_stdout(io.StringIO()):
    exec(compile(tree, str(source), 'exec'), old)
assert old['out']['all_passed']
native, mp = old['v'], old['mp']


def vectors():
    _, _, right_basis = np.linalg.svd(old['Ys'][0])
    right = right_basis[0].conj()
    left = old['Ys'][0] @ right / np.sqrt(old['y2'])
    vector = np.zeros((96, 2), complex)
    for j in range(3):
        vector[native.IDX[('u', 'R', 0, j)], 0] = right[j]
        vector[native.IDX[('u', 'L', 0, j)], 1] = left[j]
    return vector


def direct_word(probabilities, z):
    # Diagonal trace of a_R1 Q O_R a_R1^* in the actual six-mode paired marginal.
    value = 0j
    for state in range(0, 64, 2):
        prob = np.prod([probabilities[(state >> (2 * j)) & 3] for j in range(3)])
        filled = state | 1
        n_right = sum((filled >> j) & 1 for j in (0, 2, 4))
        value += prob * z**filled.bit_count() * 4 * n_right * (n_right - 1)
    return value


def pair_response(p, z):
    paired = vectors()
    model, direct, right_only = 0j, 0j, 0j
    for radial, weight in old['recovered'].items():
        if not weight:
            continue
        u = float((mp.pi / old['ell'])**2 * mp.mpf(radial.numerator) / radial.denominator / old['yq'])
        h = p * old['gamma'] - 1j * old['gamma'] @ old['phi']((0j, np.sqrt(u) + 0j))
        eigenvalues, basis = np.linalg.eigh(h)
        covariance = (basis * (.5 * (1 - np.tanh(float(old['ell']) * eigenvalues)))) @ basis.conj().T
        local = paired.conj().T @ covariance @ paired
        det = np.linalg.det(local).real
        occupation = local.diagonal().real
        assert abs(occupation.sum() - 1) < 1e-12
        probabilities = [1 - occupation.sum() + det, occupation[0] - det, occupation[1] - det, det]
        direct += float(weight) * direct_word(probabilities, z)

        energy = np.sqrt(p*p + old['y2']*u)
        t = .5 * (1 + p * np.tanh(float(old['ell']) * energy) / energy)
        d = 1 / (4 * np.cosh(float(old['ell']) * energy)**2)
        w0, w1 = d + z*(1 - t - d), t - d + z*d
        model += float(weight) * (16*z**2*w0**2*w1 + 24*z**3*w0*w1**2)
        right_only += float(weight) * (16*z**2*t*(1 - t)**2 + 24*z**3*t*t*(1 - t))
    assert abs(model - direct) < 1e-12
    assert abs(model - right_only) > 1e-4
    return {'p': p, 'z': str(z), 'prediction': str(model), 'direct_six_mode': str(direct),
            'right_marginal_only': str(right_only), 'error': abs(model - direct),
            'error_if_partner_omitted': abs(model - right_only)}


print(json.dumps({'old_control_count': old['out']['control_count'],
                  'recovered_R063_control_atoms': {str(k): str(v) for k, v in old['recovered'].items() if v},
                  'ordered_paired_holdout': [pair_response(.25, 1j), pair_response(.75, -1j)]}, indent=2))
