"""Run the machine, compare causal controls, and verify actual save/restore."""
import json
from pathlib import Path
import tempfile

import numpy as np

from machine import Config, Machine

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def run():
    checks = []
    m = Machine()
    eig = np.linalg.eigvalsh(m.native['kernel'])
    require(np.max(abs(eig - np.array([0]*25 + [1.5]*30 + [3]*5))) < 1e-11, 'Native spectrum')
    checks.append('Original signed Gram spectrum: 0^25, 1.5^30, 3^5')
    m.advance(137)
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / 'state.npz'
        m.save(path)
        restored = Machine.load(path)
        m.perturb(.1)
        restored.perturb(.1)
        m.advance(123)
        restored.advance(123)
        require(np.array_equal(m.x, restored.x) and np.array_equal(m.memory, restored.memory), 'Restore trajectory')
        require(m.history == restored.history and m.events == restored.events, 'Restore metadata')
    checks.append('Save, cold load, random perturbation and continuation reproduce the exact state')
    a = Machine(Config(feedback=0))
    b = Machine(Config(feedback=0, learning=0))
    a.advance(200)
    b.advance(200)
    require(np.array_equal(a.x, b.x), 'Memory leaks into disabled feedback')
    checks.append('With feedback=0, learning cannot affect the signal trajectory')
    require(np.max(abs(np.linalg.norm(m.x, axis=1)-1)) < 1e-12, 'Signal normalization')
    checks.append('All 60 signal vectors retain unit norm after perturbation and recursion')
    cases = [('recursive', {}), ('no_feedback', {'feedback': 0}),
             ('frozen_memory', {'learning': 0}), ('no_interaction', {'kernel': 'none'}),
             ('absolute', {'kernel': 'absolute'}), ('exact_start', {'initial_noise': 0})]
    comparisons = []
    for name, options in cases:
        machine = Machine(Config(**options))
        machine.advance(10000, measure_every=100)
        machine.save(HERE / 'runs' / (name + '.npz'))
        comparisons.append(dict(name=name, config=machine.public_state()['config'],
                                history=machine.history, final=machine.history[-1]))
        print(name, json.dumps(machine.history[-1]), flush=True)
    seeds = []
    for seed in range(20):
        machine = Machine(Config(seed=seed))
        machine.advance(5000, measure_every=250)
        seeds.append(dict(seed=seed, **machine.history[-1]))
    # Check microscopic recurrence rather than infer fixed points from stable plots.
    machine = Machine.load(HERE / 'runs' / 'recursive.npz')
    states = []
    for _ in range(128):
        machine.advance(1, measure_every=1000)
        states.append(machine.x.copy())
    lags = [dict(lag=lag, rms_distance=float(np.linalg.norm(states[-1]-states[-1-lag])/np.sqrt(60)))
            for lag in range(1, 65)]
    baseline = Machine.load(HERE / 'runs' / 'recursive.npz')
    perturbed = Machine.load(HERE / 'runs' / 'recursive.npz')
    perturbed.perturb(.1)
    recovery = []
    for steps in [0, 1, 10, 100, 1000, 5000]:
        baseline.advance(steps - (baseline.step_number - 10000), measure_every=100)
        perturbed.advance(steps - (perturbed.step_number - 10000), measure_every=100)
        rho_a = baseline.x.T @ baseline.x.conj() / 60
        rho_b = perturbed.x.T @ perturbed.x.conj() / 60
        recovery.append(dict(after_perturbation=steps,
            population_distance=float(np.linalg.norm(rho_a-rho_b)),
            signal_rms_distance=float(np.linalg.norm(baseline.x-perturbed.x)/np.sqrt(60)),
            code_weight=perturbed.history[-1]['code_weight'],
            effective_rank=perturbed.history[-1]['effective_rank']))
    report = dict(checks=checks, comparisons=comparisons, seeds=seeds, lag_distances=lags,
                  perturbation=recovery, scope='Finite experiment of the declared feedback rule')
    (HERE / 'runs' / 'experiment_results.json').write_text(json.dumps(report, indent=2)+'\n')
    ranks = [r['effective_rank'] for r in seeds]
    codes = [r['code_weight'] for r in seeds]
    print('20 seeds: ranks', min(ranks), max(ranks), 'code weights', min(codes), max(codes), flush=True)
    print('Best lag:', min(lags, key=lambda item: item['rms_distance']), flush=True)
    print('Perturbation:', json.dumps(recovery), flush=True)
    print('Checks:', json.dumps(checks), flush=True)


if __name__ == '__main__':
    run()
