"""Test separable scalar compensation on a reflected original R032 cut."""
import contextlib
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from scipy.linalg import eigh

repo = Path(__file__).resolve().parents[3]
checker = repo / "experiments/theory-contracts/chat-triangulation-20261108/check_original_whole_two_strips_20261109.py"
sys.argv = [str(checker), str(repo)]
with contextlib.redirect_stdout(io.StringIO()):
    original = runpy.run_path(str(checker))
history = original["history"][:3]
widths = original["durations"][:3]
p = original["source"]["p"]
inner_form = original["source"]["inner_form"]
grade = original["grading"]
M = original["majorana"]
I = original["identity"]


def half(h, momentum):
    steps, sea = [], 0.
    for field, width in zip(history, widths):
        K = momentum * grade - 1j * grade @ (inner_form((1+h)*np.asarray(field)) + M)
        values, vectors = eigh(K)
        steps.append((vectors * np.exp(-width*values)) @ vectors.conj().T)
        sea += width * values[:len(values)//2].sum()
    return steps, sea


def ordered(steps):
    R = I.copy()
    for step in steps:
        R = step @ R
    return R


scales = [0., .05, -.05]
log_s = original["direct_log"]("S", 0.)  # Free S is invariant under p -> -p.


def reflected_gram(momentum):
    halves = [half(h, momentum) for h in scales]
    logs = np.zeros((3, 3))
    phases = np.zeros((3, 3), dtype=complex)
    errors = []
    for i, (a, s_a) in enumerate(halves):
        for j, (b, s_b) in enumerate(halves):
            glued = ordered(a).conj().T @ ordered(b)
            direct = ordered(b + a[::-1])
            errors.append(float(np.max(np.abs(glued-direct))))
            phases[i, j], lg = np.linalg.slogdet(I+glued)
            logs[i, j] = lg+s_a+s_b-log_s
    return logs, phases, max(errors)


logs, phases, error = reflected_gram(p)
mirror_logs, mirror_phases, mirror_error = reflected_gram(-p)
mirror_residual = max(np.max(np.abs(logs-mirror_logs)),
                      np.max(np.abs(phases-mirror_phases.conj())))
assert mirror_residual < 5e-12

assert max(error, mirror_error) < 5e-12
assert abs(logs[0,0]-original["relative_s"][0]) < 1e-11
assert np.max(np.abs(phases-phases.conj().T)) < 5e-12
assert np.max(np.abs(logs-logs.T)) < 5e-12
log_cross_ratio = logs[0,0]+logs[1,1]-logs[0,1]-logs[1,0]
assert log_cross_ratio > 1e-6
rho = np.exp(-log_cross_ratio / 2)
normalized_gram = phases * np.exp(logs-(logs.diagonal()[:,None]+logs.diagonal()[None,:])/2)
gram_eigenvalues = np.linalg.eigvalsh(normalized_gram)
assert min(gram_eigenvalues) > 0
cycle_phase = np.angle(phases[0,1]*phases[1,2]*phases[2,0])
assert abs(cycle_phase) > 1e-9
print(json.dumps({"higgs_amplitude_parameters": scales,
                  "S_relative_reflected_fermion_gram_log_magnitudes": logs.tolist(),
                  "max_abs_determinant_phase": float(np.max(np.abs(np.angle(phases)))),
                  "log_cross_ratio": float(log_cross_ratio),
                  "normalized_off_diagonal_overlap": float(rho),
                  "normalized_gram_determinant": float(-np.expm1(-log_cross_ratio)),
                  "three_history_normalized_gram_eigenvalues": gram_eigenvalues.tolist(),
                  "three_history_rephasing_invariant_cycle_phase": float(cycle_phase),
                  "mirror_conjugacy_residual": float(mirror_residual),
                  "mirror_paired_normalized_gram_determinant": float(-np.expm1(-2*log_cross_ratio)),
                  "max_direct_sewing_error": max(error, mirror_error)}, indent=2))
