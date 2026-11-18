"""Decisive tests of what recursive SOM dynamics selects by itself.

No claim about all native TFPT processes: U(35) is a symmetry of this numerical
update, not generally a symmetry of the native quartic product or full field net.
"""
import json
from dataclasses import asdict

import numpy as np

from machine import HERE, Config, Machine, unit_rows


def need(value, message):
    if not value:
        raise RuntimeError(message)


def fresh_history(machine):
    machine.history = []
    machine.measure()


def haar_unitary(rng, size=35):
    a = rng.normal(size=(size, size)) + 1j*rng.normal(size=(size, size))
    q, r = np.linalg.qr(a)
    diagonal = np.diag(r)
    return q * (diagonal/abs(diagonal))[None, :]


def code_weight(x, readout):
    return float(np.linalg.norm(x@readout)**2/np.linalg.norm(x)**2)


def response_at(machine):
    """Reconstruct one frozen response, then verify it predicts real steps."""
    c = machine.config
    x, w = machine.x.copy(), machine.memory.copy()
    machine.advance(1)
    raw = (1-c.coupling)*x+c.coupling*(machine.coupler@x)
    R = ((1-c.coupling)*np.eye(60)+c.coupling*machine.coupler)/np.linalg.norm(raw, axis=1)[:, None]
    S = np.eye(60)[machine.last_winners]
    weights = machine.neighbors[machine.last_winners]
    mass = weights.sum(axis=0)
    B = weights.T/np.where(mass[:, None] > 0, mass[:, None], 1)
    L = np.diag(1-c.learning*(mass > 0))
    W_x = c.learning*B@R
    pre = (1-c.feedback)*(R@x)+c.feedback*S@(W_x@x+L@w)
    nx = np.linalg.norm(pre, axis=1)[:, None]
    A = np.block([[((1-c.feedback)*R+c.feedback*S@W_x)/nx, c.feedback*S@L/nx], [W_x, L]])
    initial = np.concatenate([x, w])
    eig = np.linalg.eigvals(A)
    periods = sorted([float(2*np.pi/abs(np.angle(z))) for z in eig
                      if abs(abs(z)-1) < 1e-8 and z.imag > 1e-7])
    machine.advance(999)
    prediction = np.linalg.matrix_power(A, 1000)@initial
    error = float(np.linalg.norm(prediction-np.concatenate([machine.x, machine.memory]))/np.linalg.norm(prediction))
    return dict(periods_steps=periods, frozen_prediction_relative_error=error,
                prediction_valid_in_test=error < 1e-9)


def equivariance():
    rng = np.random.default_rng(20261026)
    u = haar_unitary(rng)
    a, b = Machine(), Machine()
    b.x = b.x@u
    b.memory = b.memory@u
    fresh_history(b)
    a.advance(5000, measure_every=100)
    b.advance(5000, measure_every=100)
    error_x = float(np.linalg.norm(b.x-a.x@u)/np.linalg.norm(a.x))
    error_w = float(np.linalg.norm(b.memory-a.memory@u)/np.linalg.norm(a.memory))
    need(error_x < 1e-10 and error_w < 1e-10, 'Unitary covariance failed')
    need(np.array_equal(a.last_winners, b.last_winners), 'Winner covariance failed')
    code = a.native['code']
    fixed = code_weight(b.x, code)
    transported = code_weight(b.x, u.conj().T@code)
    original = code_weight(a.x, code)
    need(abs(transported-original) < 1e-10, 'Passive readout transport failed')
    # A larger orientation population needs no duplicate run: exact equivariance
    # transports the already tested trajectory. Only one trajectory is rerun above.
    orientations = [code_weight(a.x@haar_unitary(rng), code) for _ in range(128)]
    return dict(original=a.history[-1], rotated=b.history[-1], signal_relative_error=error_x,
                memory_relative_error=error_w, fixed_code_weight=fixed,
                transported_code_weight=transported, original_code_weight=original,
                orientation_samples=128, orientation_code_mean=float(np.mean(orientations)),
                orientation_code_min=float(min(orientations)), orientation_code_max=float(max(orientations)),
                isotropic_expected_code_weight=5/35,
                scope='U(35) covariance of the declared update; generic U does not preserve native product inputs')


def generic_starts():
    results = []
    for seed in [0, 7, 19]:
        m = Machine(Config(seed=seed))
        rng = np.random.default_rng(1000+seed)
        m.x = unit_rows(rng.normal(size=(60, 35))+1j*rng.normal(size=(60, 35)))
        m.memory = m.x[rng.permutation(60)].copy()
        fresh_history(m)
        m.advance(5000, measure_every=100)
        results.append(dict(seed=seed, **m.history[-1]))
    return results


def rates():
    results = []
    # Factors of two around the implemented rates; no tuning to physical values.
    for feedback in [.125, .25, .5]:
        for learning in [.075, .15, .3]:
            m = Machine(Config(feedback=feedback, learning=learning))
            m.advance(5000, measure_every=100)
            stats = dict(m.history[-1])
            response = response_at(m)
            results.append(dict(config=asdict(m.config), at_5000=stats, response=response))
    return results


def main():
    result = dict(equivariance=equivariance(), generic_starts=generic_starts(), rates=rates())
    path = HERE/'runs'/'fundamental_audit.json'
    path.write_text(json.dumps(result, indent=2)+'\n')
    e = result['equivariance']
    print('U35 covariance:', e['signal_relative_error'], e['memory_relative_error'])
    print('Original/fixed rotated/passively transported code:',e['original_code_weight'],e['fixed_code_weight'],e['transported_code_weight'])
    print('128 orientations mean/min/max:',e['orientation_code_mean'],e['orientation_code_min'],e['orientation_code_max'])
    for row in result['generic_starts']:
        print('Generic:', row['seed'], row['effective_rank'], row['code_weight'])
    for row in result['rates']:
        c,s,r = row['config'],row['at_5000'],row['response']
        print('Rates:',c['feedback'],c['learning'],'rank',s['effective_rank'],'code',s['code_weight'], 'change',s['change'],'periods',r['periods_steps'],'prediction',r['frozen_prediction_relative_error'])


if __name__ == '__main__':
    main()
