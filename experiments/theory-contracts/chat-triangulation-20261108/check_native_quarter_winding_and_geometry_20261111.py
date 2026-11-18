"""Check the original winding on marked quarter arms and audit their geometry inputs."""

from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from scipy.linalg import eigh, solve
import sympy as sp


root = Path(sys.argv[1]).resolve()
base = root / 'experiments/theory-contracts/chat-triangulation-20261108'
with redirect_stdout(io.StringIO()):
    shared = runpy.run_path(str(base / 'check_joint_reflected_responses_20261111.py'))

source = shared['source']
frames = source['frames']
grading = source['grading']
identity = shared['identity']
idx = source['native'].IDX
weak, hypercharge = shared['gauge_operators']()
photon_charge = weak[2] + hypercharge
raising = weak[0] + 1j * weak[1]


def maxnorm(matrix):
    return float(np.max(np.abs(matrix)))


assert maxnorm(photon_charge @ raising - raising @ photon_charge - raising) < 1e-13
u = [idx[('u', 'L', 0, family)] for family in range(3)]
d = [idx[('d', 'L', 0, family)] for family in range(3)]
assert maxnorm(raising[np.ix_(u, d)] - np.eye(3)) < 1e-13


def arm(momentum, reference):
    steps, sea = [], 0.
    for profile, duration in zip(source['history'], source['durations'], strict=True):
        h = momentum * grading
        if not reference:
            h = h - 1j * grading @ (source['inner_form'](profile) + source['majorana'])
        assert maxnorm(h - h.conj().T) < 1e-13
        eigenvalues, basis = eigh(h)
        steps.append((basis * np.exp(-duration * eigenvalues)) @ basis.conj().T)
        sea += duration * eigenvalues[:len(grading)//2].sum()
    return steps, sea


def direct(steps):
    result = identity.copy()
    for step in steps:
        result = step @ result
    return result


def sew(steps, cuts):
    pieces = []
    for begin, end in zip((0, *cuts), (*cuts, len(steps)), strict=True):
        piece = identity.copy()
        for i in range(begin, end):
            transport = frames[i].conj().T @ frames[i - 1] if i else identity
            piece = (frames[i].conj().T @ steps[i] @ frames[i]) @ transport @ piece
        pieces.append(piece)
    return direct(pieces)


def relative_cross(log_radius=0.):
    momenta = (.25 * np.exp(-log_radius), .75 * np.exp(-log_radius))
    a, b = (arm(p, False) for p in momenta)
    s_a, s_b = (arm(p, True) for p in momenta)
    sea_shift = 2 * (a[1] + b[1] - s_a[1] - s_b[1])

    def gram(transfer):
        return identity + transfer.conj().T @ transfer

    def logtrace(transfer):
        sign, logdet = np.linalg.slogdet(gram(transfer))
        assert abs(sign - 1) < 1e-12
        return logdet

    reference = logtrace(direct(s_a[0])) + logtrace(direct(s_b[0]))

    def response(first, second, mark):
        cross = second.conj().T @ mark @ first
        dressed = solve(gram(second), cross, assume_a='pos')
        dressed = solve(gram(first), dressed.conj().T, assume_a='pos').conj().T
        log_empty = sea_shift + logtrace(first) + logtrace(second) - reference
        return dressed, log_empty

    # u_cov: W0 -> W1 on this Fourier mode. Its matrix coefficient is one.
    # The one-way mark is off-diagonal. The derivative of the closed trace
    # vanishes, while the charged open leg gives the resolvent-dressed cross.
    direct_pair, log_empty = response(direct(a[0]), direct(b[0]), raising)
    endpoint_mark = frames[-1].conj().T @ raising @ frames[-1]
    errors, weight_errors = [], []
    for cuts in ((2,), (4,), (2, 4)):
        sewn_a, sewn_b = sew(a[0], cuts), sew(b[0], cuts)
        sewn_pair, sewn_weight = response(sewn_a, sewn_b, endpoint_mark)
        errors.append(maxnorm(direct_pair - frames[0] @ sewn_pair @ frames[0].conj().T))
        weight_errors.append(abs(log_empty - sewn_weight))
    assert max(errors) < 1e-12 and max(weight_errors) < 1e-11
    block = direct_pair[np.ix_(u, d)]
    assert np.linalg.matrix_rank(block, tol=1e-11) == 3
    return {'reference': 'common free reflected S Fock trace, same p and durations per arm',
            'relative_log_empty': log_empty,
            'uLdL_family_singular_values': np.linalg.svd(block, compute_uv=False).tolist(),
            'uLdL_max_imaginary_entry': maxnorm(block.imag), 'cut_max_errors': errors,
            'cut_log_weight_errors': weight_errors,
            'unnormalized_uLdL_max_entry_over_S': float(np.exp(log_empty)) * maxnorm(block)}, block


def exact_checks():
    n, width, length, k = sp.symbols('n W L k', positive=True)
    generator = 2 * sp.pi * (k + sp.Rational(1, 2)) / length
    log_odds = 2 * width * generator
    uniform_scaling = sp.cancel((width * sp.diff(log_odds, width) + length * sp.diff(log_odds, length)))
    assert uniform_scaling == 0
    # The old annulus only sees W/L. Length stretch with W fixed and angular
    # period stretch with W fixed are the same tangent in that formula.
    assert sp.simplify(sp.diff(log_odds.subs(length, n*length), n).subs(n, 1)
                       - length*sp.diff(log_odds, length)) == 0

    A, G, f_r, grav_shift, photon_shift, j = sp.symbols('A G f_R s c j', positive=True)
    i_ref = sp.Function('I_ref')(n, j)
    i_rel = sp.Function('I_rel')(n, j)
    i_adm = i_rel + i_ref
    def angle(action):
        return (n * sp.diff(action, n) - action).subs({n: 1, j: 0})
    assert sp.simplify(angle(i_adm) - angle(i_rel) - angle(i_ref)) == 0
    change = (n-1)*grav_shift*A + photon_shift*j*j/2
    assert change.subs({n: 1, j: 0}) == 0
    assert angle(change) == grav_shift*A
    assert sp.diff(change, j, 2) == photon_shift
    assert sp.diff(change, j, 2, n) == 0
    assert angle((n-1)*A*f_r/(4*G))/A == f_r/(4*G)

    euler, a1, a2, sd = sp.symbols('eta A_1 A_2 S_dS', positive=True)
    chi_bulk, chi_horizon = sp.symbols('chi_bulk chi_horizon')
    euler_replica = -euler * (n*chi_bulk + (1-n)*chi_horizon)
    assert sp.simplify(angle(euler_replica) - euler * chi_horizon) == 0
    s1 = a1*f_r/(4*G) + 2*euler
    s2 = a2*f_r/(4*G) + 2*euler
    assert sp.cancel((s2-s1)/(a2-a1)) == f_r/(4*G)
    ratio = (sp.Rational(2, 3)*sd + 4*euler)/(sd + 2*euler)
    assert sp.cancel(ratio - sp.Rational(2, 3) - 8*euler/(3*(sd + 2*euler))) == 0

    V = sp.diag(sp.I, -sp.I)
    carry = sp.Symbol('v', nonzero=True)
    winding = sp.Matrix([[0, carry], [1, 0]])
    assert V*winding*V.inv() == -winding and winding*winding == carry*sp.eye(2)
    return {'V_conjugates_original_winding_to': '-u_cov',
            'one_turn_folds': '[[0,u_base],[1,0]]', 'common_annulus_scale_jet': str(uniform_scaling),
            'circle_stretch_and_angle_at_fixed_width': 'same tangent in the original angular log odds',
            'absolute_conical_response': 'relative conical response + reference conical response',
            'separate_local_response_shifts': ['s A', 'c'], 'wald_target_per_area': 'f_R/(4G)',
            'euler_entropy_nariai_and_ds': ['4 eta', '2 eta'],
            'same_topology_fixed_curvature_area_difference': 'f_R/(4G)'}


baseline, _ = relative_cross()
metric = []
for step in (.001, .0005, .00025):
    up, upper = relative_cross(step)
    down, lower = relative_cross(-step)
    metric.append({'step': step, 'd_log_empty': (up['relative_log_empty'] - down['relative_log_empty']) / (2*step),
                   'd_cross_block_max': maxnorm((upper - lower) / (2*step))})
assert abs(metric[-1]['d_log_empty'] - metric[-2]['d_log_empty']) < 1e-6
assert abs(metric[-1]['d_cross_block_max'] - metric[-2]['d_cross_block_max']) < 1e-6
print(json.dumps({'conditional_R032_quarter_cross': baseline, 'same_branch_smooth_radius_response': metric,
                  'exact_original_structure': exact_checks()}, indent=2))
