"""Explain ongoing motion in a saved run by a frozen response operator.

This is not the full Jacobian. Winners and row normalization factors are frozen
at the observed state, then its predicted trajectory is checked against the
actual nonlinear machine. No frequencies are inputs to the machine.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from machine import HERE, Machine


def inspect(checkpoint, steps):
    machine = Machine.load(checkpoint)
    start = machine.step_number
    x, w, config = machine.x.copy(), machine.memory.copy(), machine.config
    machine.advance(1)
    raw = (1-config.coupling)*x + config.coupling*(machine.coupler@x)
    response = ((1-config.coupling)*np.eye(60) + config.coupling*machine.coupler)
    response /= np.linalg.norm(raw, axis=1)[:, None]
    select = np.eye(60)[machine.last_winners]
    weights = machine.neighbors[machine.last_winners]
    mass = weights.sum(axis=0)
    batch = weights.T / np.where(mass[:, None] > 0, mass[:, None], 1)
    retain = np.diag(1-config.learning*(mass > 0))
    learned_response = config.learning*batch@response
    pre = (1-config.feedback)*(response@x) + config.feedback*select@(learned_response@x+retain@w)
    norms = np.linalg.norm(pre, axis=1)[:, None]
    operator = np.block([
        [((1-config.feedback)*response + config.feedback*select@learned_response)/norms,
         (config.feedback*select@retain)/norms],
        [learned_response, retain]])
    state = np.concatenate([x, w])
    eig = np.linalg.eigvals(operator)
    modes = [dict(real=float(z.real), imag=float(z.imag), modulus=float(abs(z)),
                  frozen_response_period_steps=float(2*np.pi/abs(np.angle(z))))
             for z in eig if abs(abs(z)-1) < 1e-8 and z.imag > 1e-7]
    machine.advance(999)
    predicted = np.linalg.matrix_power(operator, 1000)@state
    error = float(np.linalg.norm(predicted-np.concatenate([machine.x, machine.memory]))/np.linalg.norm(predicted))
    if error > 1e-9:
        raise RuntimeError(f'Frozen response does not reproduce this trajectory: {error}')
    report = dict(frozen_response_modes=modes, frozen_1000_step_relative_prediction_error=error, history=[])
    endpoints = sorted({start+1000, *[start+n for n in [10000, 40000, steps] if n <= steps]})
    for endpoint in endpoints:
        machine.advance(endpoint-machine.step_number, measure_every=100)
        report['history'].append(dict(**machine.history[-1],
            signal_distance_from_start=float(np.linalg.norm(machine.x-x)/np.sqrt(60))))
    report['full_history'] = machine.history
    machine.save(HERE/'runs'/'long_run.npz')
    (HERE/'runs'/'long_run.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({key: value for key, value in report.items() if key != 'full_history'}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checkpoint', type=Path, default=HERE/'runs'/'recursive.npz')
    parser.add_argument('--steps', type=int, default=90000)
    args = parser.parse_args()
    if args.steps < 1000:
        parser.error('--steps must be at least 1000')
    inspect(args.checkpoint, args.steps)
