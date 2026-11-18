"""Test reflected sewing on R032 and the constant-curvature horizon comparison."""

from contextlib import redirect_stdout
from functools import cache
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from scipy.linalg import eigh, solve
import sympy as sp


repo = Path(sys.argv[1]).resolve()
base = repo / 'experiments/theory-contracts/chat-triangulation-20261108'
with redirect_stdout(io.StringIO()):
    source = runpy.run_path(str(base / 'check_joint_reflected_responses_20261111.py'))

identity = source['identity']
frames = source['frames']
charges = source['charges']
index = source['source']['native'].IDX
families = {kind: [[index[(kind, side, 0, j)] for j in range(3)] for side in ('L', 'R')]
            for kind in ('u', 'd', 'e')}


def maxnorm(matrix):
    return float(np.max(np.abs(matrix)))


def logdet(matrix):
    sign, value = np.linalg.slogdet(matrix)
    assert abs(sign - 1) < 1e-12
    return float(value)


def direct(steps, twist):
    kernel = identity.copy()
    for i, step in enumerate(steps):
        kernel = step @ (twist if i == source['link_before_step'] else identity) @ kernel
    gram = kernel.conj().T @ kernel
    return logdet(identity + gram), identity - solve(identity + gram, identity, assume_a='pos')


def strips(steps, twist, cuts):
    out = []
    for begin, end in zip((0, *cuts), (*cuts, len(steps))):
        kernel = identity.copy()
        for i in range(begin, end):
            link = identity
            if i:
                link = frames[i].conj().T @ (twist if i == source['link_before_step'] else identity) @ frames[i-1]
            kernel = (frames[i].conj().T @ steps[i] @ frames[i]) @ link @ kernel
        out.append(kernel)
    return out


def join_at_cut(early, late):
    # Compute the Fock determinant and return kernel on the cut, without BA.
    p, q = early @ early.conj().T, late.conj().T @ late
    cp = identity - solve(identity + p, identity, assume_a='pos')
    cq = identity - solve(identity + q, identity, assume_a='pos')
    overlap = (identity - cp) @ (identity - cq) + cp @ cq
    log_empty = logdet(identity + p) + logdet(identity + q) + logdet(overlap)
    density = early.conj().T @ q @ solve(identity + p @ q, early)
    return log_empty, density, logdet(overlap)


@cache
def evaluate(higgs, angle, geometry, reference):
    steps_and_sea = [source['heat_step'](i, higgs if i == source['higgs_step'] else 0.,
                                       geometry if i == source['geometry_step'] else 0., reference)
                     for i in range(len(frames))]
    steps = [step for step, _ in steps_and_sea]
    sea = 2 * sum(sea for _, sea in steps_and_sea)
    twist = np.diag(np.exp(1j * angle * charges))
    exact, density = direct(steps, twist)
    routes = {}
    for cuts in ((2,), (4,), (2, 4)):
        pieces = strips(steps, twist, cuts)
        if len(pieces) == 3:
            middle = pieces[1]
            a = join_at_cut(middle @ pieces[0], pieces[2])
            b = join_at_cut(pieces[0], pieces[2] @ middle)
            assert abs(a[0] - b[0]) < 1e-11 and maxnorm(a[1] - b[1]) < 1e-11
        else:
            a = join_at_cut(*pieces)
        other_density = frames[0] @ a[1] @ frames[0].conj().T
        errors = {kind: maxnorm(density[np.ix_(right, left)] - other_density[np.ix_(right, left)])
                  for kind, (left, right) in families.items()}
        assert maxnorm(other_density - density) < 1e-11 and abs(exact - a[0]) < 1e-11
        routes[str(cuts)] = {'log': a[0] + sea, 'full_density_error': maxnorm(other_density - density),
                             'flavor_block_errors': errors, 'cut_log_overlap': a[2]}
    blocks = {kind: maxnorm(density[np.ix_(right, left)]) for kind, (left, right) in families.items()}
    return {'direct_log': exact + sea, 'routes': routes, 'flavor_block_maxima': blocks}


def difference(step):
    def value(h, a, g):
        full, empty = evaluate(h, a, g, False), evaluate(0., a, g, True)
        return {'direct': full['direct_log'] - empty['direct_log'],
                **{key: row['log'] - empty['routes'][key]['log'] for key, row in full['routes'].items()}}

    rows = {}
    for h in (-step, step):
        for g in (-step, step):
            low, mid, high = (value(h, a, g) for a in (-step, 0., step))
            for key in low:
                rows[key] = rows.get(key, 0.) + np.sign(h * g) * (low[key] + high[key] - 2 * mid[key])
    return {key: val / (4 * step**4) for key, val in rows.items()}


def endpoint_phase():
    kernel = identity.copy()
    for i in range(len(frames)):
        heat, _ = source['heat_step'](i, 0., 0., False)
        kernel = heat @ kernel
    twist = np.diag(np.exp(1j * np.pi * charges))
    other = twist @ kernel
    gram = kernel.conj().T @ kernel
    diagonal_error = maxnorm(gram - other.conj().T @ other)
    mixed_error = maxnorm(gram - kernel.conj().T @ other)
    transported_error = maxnorm(gram - kernel.conj().T @ twist.conj().T @ other)
    assert diagonal_error < 1e-12 and transported_error < 1e-12 and mixed_error > 0.1
    return {'diagonal_error': diagonal_error, 'mixed_difference_fixed_terminal_frame': mixed_error,
            'mixed_error_with_terminal_frame_transported': transported_error}


def horizon():
    n, area, bulk, shift, lam, mass, grav = sp.symbols('n A I B Lambda Ms G', positive=True)
    r = sp.Symbol('R')
    action = r - 2 * lam + r**2 / (6 * mass**2)
    f_r = sp.diff(action, r).subs(r, 4 * lam)
    replica = n * (bulk + shift) + (n - 1) * area * f_r / (4 * grav)
    entropy = sp.simplify((n * sp.diff(replica, n) - replica).subs(n, 1))
    assert sp.simplify(entropy - area * f_r / (4 * grav)) == 0
    nariai_total = entropy.subs(area, 8 * sp.pi / lam)
    ds = entropy.subs(area, 12 * sp.pi / lam)
    assert sp.cancel(nariai_total / ds) == sp.Rational(2, 3)
    nariai_action = -action.subs(r, 4 * lam) * 16 * sp.pi**2 / lam**2 / (16 * sp.pi * grav)
    assert sp.simplify(nariai_action + nariai_total) == 0
    x = sp.Symbol('x', positive=True)
    temperature_ratio = (2*x + 1) / (x*(x + 2))
    ratio_gap = (1-x)*(x+1)/(x*(x+2))
    assert sp.cancel(temperature_ratio - 1 - ratio_gap) == 0
    return {'entropy_per_area': str(sp.factor(entropy / area)), 'smooth_bulk_shift_response': str(sp.diff(entropy, shift)),
            'nariai_each': str(sp.factor(nariai_total / 2)), 'nariai_over_ds': '2/3',
            'sds_temperature_ratio_minus_one': str(ratio_gap)}


full, empty = evaluate(0., 0., 0., False), evaluate(0., 0., 0., True)
contrast = difference(.06)
assert max(abs(v - contrast['direct']) for v in contrast.values()) < 1e-7
print(json.dumps({'R032': {'A': full, 'S': empty, 'log_relative': full['direct_log'] - empty['direct_log'],
                           'mixed_relative_contrast': contrast, 'relative_endpoint_control': endpoint_phase()},
                  'constant_curvature_replica': horizon()}, indent=2))
