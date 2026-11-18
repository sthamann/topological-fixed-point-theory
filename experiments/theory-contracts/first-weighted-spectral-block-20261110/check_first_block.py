"""Replay the first spectral shell using the unmodified R113 full96 builder."""

from contextlib import redirect_stdout
from pathlib import Path
import io
import json
import runpy
import sys

import numpy as np
import sympy as sp


repo = Path(sys.argv[1]).resolve()
original = repo / "experiments/theory-contracts/chat-triangulation-20261108"
sys.path.insert(0, str(repo / "verification"))
entrypoint = "check_r113_dressed_open_20261109.py"
sys.argv = [entrypoint, str(original)]
with redirect_stdout(io.StringIO()):
    source = runpy.run_path(str(original / entrypoint))

gamma = source["grading"]
majorana = source["majorana"]
higgs = source["inner_form"](np.array([0., 1.]))
assert np.count_nonzero(gamma @ higgs + higgs @ gamma) == 0
assert np.count_nonzero(gamma @ majorana + majorana @ gamma) == 0
graded_plus = np.flatnonzero(gamma.diagonal() == 1)
graded_minus = np.flatnonzero(gamma.diagonal() == -1)
assert len(graded_plus) == len(graded_minus) == 48
off_diagonal = np.ix_(graded_plus, graded_minus)
b_y = higgs[off_diagonal]
b_m = majorana[off_diagonal]

def exact(matrix):
    return sp.Matrix([[sp.Rational(float(x.real)) + sp.I * sp.Rational(float(x.imag))
                       for x in row] for row in matrix])

yu, yd, ye, ynu = source["qref"].yukawa
assert all(exact(y).det() != 0 for y in (yu, yd, ye, ynu))
def positive_definite(matrix):
    return all(sp.expand(matrix[:n, :n].det()) > 0 for n in range(1, matrix.rows + 1))

def column_blocks(matrix):
    # Connected column components in the exact Gram sparsity.
    touching = np.abs(matrix) != 0
    adjacency = touching.T @ touching
    remaining = set(range(matrix.shape[1]))
    while remaining:
        block = {remaining.pop()}
        frontier = set(block)
        while frontier:
            column = frontier.pop()
            neighbors = remaining.intersection(np.flatnonzero(adjacency[column]))
            remaining.difference_update(neighbors)
            block.update(neighbors)
            frontier.update(neighbors)
        yield sorted(block)

yukawa_floor, majorana_ceiling = sp.Rational(33, 250), sp.Rational(21, 5)
for columns in column_blocks(b_y):
    y = exact(b_y[:, columns])
    assert positive_definite(y.H * y - yukawa_floor**2 * sp.eye(len(columns)))
for columns in column_blocks(b_m):
    m = exact(b_m[:, columns])
    assert positive_definite(majorana_ceiling**2 * sp.eye(len(columns)) - m.H * m)

s_y = np.linalg.svd(b_y, compute_uv=False)
s_m = np.linalg.svd(b_m, compute_uv=False)
radius = 64
assert radius * yukawa_floor > 2 * majorana_ceiling
far = radius * higgs + majorana
s_far = np.linalg.svd(far[off_diagonal], compute_uv=False)
assert radius * s_y.min() - s_m.max() > s_m.max()
assert s_far.min() > s_m.max()

t = float(source["qref"].time)

def shell(p, phi, singular):
    h = p*gamma - 1j*gamma@phi
    eig2 = np.linalg.eigvalsh(h@h)
    squared = p*p + singular*singular
    np.testing.assert_allclose(eig2, np.sort(np.repeat(squared, 2)), atol=5e-9)
    heat_log_ratio = 2*np.sum(np.logaddexp(0., -t*np.sqrt(squared))
                              - np.logaddexp(0., -t*abs(p)))
    spectral_trace_s1 = 2*np.sum(np.exp(-squared) - np.exp(-p*p))
    spectral_trace_thermal = np.sum(np.logaddexp(0., -t*np.sqrt(np.maximum(0., eig2)))
                                    - np.logaddexp(0., -t*abs(p)))
    np.testing.assert_allclose(heat_log_ratio, spectral_trace_thermal, atol=5e-14)
    return float(spectral_trace_s1), float(heat_log_ratio)

rows = {}
for quarter in (.25, .75):
    # Cut off |p|, so reflection (alpha,n)->(1-alpha,-n-1)
    # carries both the retained and newly added modes exactly.
    initial_n = 0 if quarter == .25 else -1
    added_n = -1 if quarter == .25 else 0
    retained_p, added_p = initial_n + quarter, added_n + quarter
    assert abs(retained_p) == .25 and abs(added_p) == .75
    start = shell(added_p, majorana, s_m)
    end = shell(added_p, far, s_far)
    assert end[0] < start[0] and end[1] < start[1]
    rows[str(quarter)] = {
        "retained_momentum": retained_p,
        "added_momentum": added_p,
        "relative_spectral_action_slice_s1_H0": start[0],
        "relative_spectral_action_slice_s1_Hfar": end[0],
        "bare_log_relative_reweighting_Hfar_vs_H0_s1": start[0]-end[0],
        "fermion_sea_corrected_log_shell_H0": start[1],
        "fermion_sea_corrected_log_shell_Hfar": end[1],
        "net_log_reweighting_slice_s1_plus_fermion": start[0]-end[0]+end[1]-start[1],
    }
assert all(rows["0.25"][key] == rows["0.75"][key]
           for key in rows["0.25"] if "momentum" not in key)

# The existing Q plugin supplies this spatial first harmonic and all 384 entries.
runtime = repo / ("experiments/theory-contracts/original-paper-source-reconciliation-20261102/"
                  "contribution_20261106Q/runtime")
sys.path.insert(0, str(runtime))
from assemble import Assembler

assembler = Assembler(runtime.parent / "domain/native-diagnostic-v1.json")
candidate = {"schema": "tfpt.field-configuration.v1",
             "domain_sha256": assembler.domain["_sha256"], "genes": [2, 1, 1]}
assembled = assembler.assemble(candidate)
assert assembled["modes"] == [-2, -1, 0, 1]
assert all(error < 1e-14 for error in assembled["kinematic_residuals"].values())
h = assembled["h_N"]
low = np.arange(96, 288)
high = np.r_[0:96, 288:384]
retained, eliminated = h[np.ix_(low, low)], h[np.ix_(high, high)]
cross = h[np.ix_(low, high)]
np.testing.assert_array_equal(cross[:96, 96:], 0)
np.testing.assert_array_equal(cross[96:, :96], 0)
np.testing.assert_array_equal(eliminated[:96, 96:], 0)
np.testing.assert_allclose(cross[:96,:96], (-1j*assembler.gamma@assembler.V[0])*13/200,
                           atol=1e-15)
assert np.linalg.svd(cross, compute_uv=False).min() > 0
z = 1j
low_z = z*np.eye(len(low))-retained
high_z = z*np.eye(len(high))-eliminated
self_energy = cross @ np.linalg.solve(high_z, cross.conj().T)
effective_z = low_z-self_energy
assert np.linalg.eigvalsh((self_energy-self_energy.conj().T)/(2j)).max() < 0
full_z = z*np.eye(len(h))-h
resolvent = np.linalg.inv(full_z)[np.ix_(low,low)]
coarse_resolvent = np.linalg.inv(low_z)
np.testing.assert_allclose(resolvent, np.linalg.inv(effective_z), atol=1e-14)
full_sign, full_logabs = np.linalg.slogdet(full_z)
shell_sign, shell_logabs = np.linalg.slogdet(high_z)
effective_sign, effective_logabs = np.linalg.slogdet(effective_z)
np.testing.assert_allclose(full_logabs, shell_logabs+effective_logabs, atol=1e-12)
np.testing.assert_allclose(full_sign, shell_sign*effective_sign, atol=1e-12)
mode_result = {
    "role": assembled["h_N_role"],
    "actual_proposal_genes": candidate["genes"],
    "cross_mode_rank": int(np.linalg.matrix_rank(cross)),
    "return_gram_frobenius_norm": float(np.linalg.norm(cross@cross.conj().T)),
    "bare_compression_resolvent_error_frobenius": float(np.linalg.norm(resolvent-coarse_resolvent)),
    "schur_resolvent_error_max_entry": float(np.max(np.abs(resolvent-np.linalg.inv(effective_z)))),
    "schur_determinant_logabs_error": float(abs(full_logabs-shell_logabs-effective_logabs)),
}

print(json.dumps({
    "scope": "Original constant full96 circle coefficient family, not weighted field integration.",
    "thermal_time_already_given": t,
    "exact_full_yukawa_singular_lower_bound": str(yukawa_floor),
    "exact_majorana_norm_upper_bound": str(majorana_ceiling),
    "s_min_full_yukawa": float(s_y.min()),
    "majorana_operator_norm": float(s_m.max()),
    "comparison_H_radius": radius,
    "s_min_full_phi_at_radius": float(s_far.min()),
    "original_quarter_shells": rows,
    "nonconstant_existing_Q_first_harmonic": mode_result,
    "thermal_profile_f0": float(np.log(2.)),
    "thermal_profile_scale_invariant_moment_ratio": float(5 * sp.log(2) * sp.zeta(5) / sp.zeta(3)**2),
}, indent=2))
