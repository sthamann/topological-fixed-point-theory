"""Mixed metric/Higgs response of an actual R032 matter mark on a fixed history."""
import contextlib
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from scipy.linalg import eigh, expm, solve

repo = Path(__file__).resolve().parents[3]
old_checker = repo / "experiments/theory-contracts/chat-triangulation-20261108/check_original_whole_two_strips_20261109.py"
sys.argv = [str(old_checker), str(repo)]
with contextlib.redirect_stdout(io.StringIO()):
    original = runpy.run_path(str(old_checker))

history = original["history"]
durations = original["durations"]
inner_form = original["source"]["inner_form"]
p0 = original["source"]["p"]
M = original["majorana"]
grade = original["grading"]
I = original["identity"]
n = len(I)
right = original["charged_right"]
reference = original["relative_s"]
assert len(right) == 3

def insertion_jets():
    out = []
    for field, width in zip(history, durations):
        higgs = inner_form(field)
        a = -width * (p0 * grade - 1j * grade @ (higgs + M))
        a_metric = width * p0 * grade
        a_higgs = 1j * width * grade @ higgs
        z = np.zeros_like(a)
        out.append(expm(np.block([[a, a_metric, a_higgs, z],
                                   [z, a, z, a_higgs],
                                   [z, z, a, a_metric],
                                   [z, z, z, a]])))
    return out

def product(items):
    result = np.eye(4 * n, dtype=complex)
    for item in items:
        result = item @ result
    return result

def read_jet(jet):
    value = jet[:n, :n]
    metric = jet[:n, n:2*n]
    higgs = jet[:n, 2*n:3*n]
    mixed = jet[:n, 3*n:]
    G = solve(I + value, I)
    cov = I - G
    cov_metric = G @ metric @ G
    cov_higgs = G @ higgs @ G
    cov_mixed = (G @ mixed @ G - G @ metric @ G @ higgs @ G
                 - G @ higgs @ G @ metric @ G)
    return np.array([np.trace(x[np.ix_(right, right)]) for x in
                     (cov, cov_metric, cov_higgs, cov_mixed)])

def direct(metric, higgs):
    transition = I.copy()
    sea = 0.0
    for field, width in zip(history, durations):
        p = p0 * np.exp(-metric)
        phi = inner_form((1 + higgs) * np.asarray(field)) + M
        h = p * grade - 1j * grade @ phi
        vals, vectors = eigh(h)
        sea += width * vals[:n // 2].sum()
        transition = ((vectors * np.exp(-width * vals)) @ vectors.conj().T) @ transition
    C = solve(I + transition, transition)
    return np.trace(C[np.ix_(right, right)]), sea + np.linalg.slogdet(I + transition)[1]

factors = insertion_jets()
full = product(factors)
jet = read_jet(full)
split_errors = [np.max(np.abs(read_jet(product(factors[k:]) @ product(factors[:k])) - jet))
                for k in range(1, len(factors))]
assert max(split_errors) < 4e-11
assert np.max(np.abs(jet.imag)) < 2e-11
step = 2e-4
def central_jet(f):
    return np.array([
        f(0, 0),
        (f(step, 0) - f(-step, 0)) / (2 * step),
        (f(0, step) - f(0, -step)) / (2 * step),
        (f(step, step) - f(step, -step) - f(-step, step)
         + f(-step, -step)) / (4 * step**2),
    ])

central = central_jet(lambda sigma, h: direct(sigma, h)[0])
assert np.max(np.abs(central - jet)) < 2e-7
assert abs(jet[3]) > 1e-4

# Keep the actual filled-sea factor and the same fixed S reference when
# converting the conditional matter mark into an unnormalized contraction.
sea_h = sea_mh = 0.0
for field, width in zip(history, durations):
    H = inner_form(field)
    phi = H + M
    Q = p0*p0*I + phi @ phi
    eig, vec = eigh(Q)
    diag = np.diag(vec.conj().T @ (phi @ H + H @ phi) @ vec).real
    sea_h -= width / 4 * np.sum(diag / np.sqrt(eig))
    sea_mh -= width * p0**2 / 4 * np.sum(diag / eig**1.5)

G = solve(I + full[:n, :n], I)
R_m = full[:n, n:2*n]
R_h = full[:n, 2*n:3*n]
R_mh = full[:n, 3*n:]
L_h = sea_h + np.trace(G @ R_h).real
L_mh = sea_mh + np.trace(G @ R_mh - G @ R_h @ G @ R_m).real
L = np.array([reference[0], reference[1], L_h, L_mh])
def direct_unnormalized(sigma, h):
    value, log_a = direct(sigma, h)
    return value * np.exp(log_a - original["direct_log"]("S", sigma))

N = np.exp(L[0]) * np.array([
    jet[0], jet[1] + L[1] * jet[0], jet[2] + L[2] * jet[0],
    jet[3] + L[1] * jet[2] + L[2] * jet[1] + (L[3] + L[1] * L[2]) * jet[0],
])
N_direct = central_jet(direct_unnormalized)
assert np.max(np.abs(N - N_direct)) < 2e-11
print(json.dumps({
    "scope": "actual fixed reflected v252 history; all one-particle intermediate directions retained",
    "mark": "sum of right charged-lepton particle occupations in the actual three-family chart",
    "empty_relative_metric_jet_S": reference.tolist(),
    "normalized_matter_value_metric_Higgs_mixed": jet.real.tolist(),
    "independent_spectral_finite_difference": central.real.tolist(),
    "log_S_relative_trace_metric_Higgs_mixed": L.tolist(),
    "S_relative_unnormalized_matter_metric_Higgs_mixed": N.real.tolist(),
    "independent_unnormalized_finite_difference": N_direct.real.tolist(),
    "whole_vs_split_errors": [float(x) for x in split_errors],
}, indent=2))
