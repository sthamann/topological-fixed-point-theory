"""Apply the paper's covariance reader on the unchanged reflected R032 chart."""
import contextlib
import io
import json
import runpy
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

repo = Path(sys.argv[1]).resolve()
directory = repo / "experiments/theory-contracts/chat-triangulation-20261108"
previous = directory / "check_original_whole_two_strips_20261109.py"
argv = sys.argv
try:
    sys.argv = [str(previous), str(repo)]
    with contextlib.redirect_stdout(io.StringIO()):
        source = runpy.run_path(str(previous))
finally:
    sys.argv = argv

sys.path.insert(0, str(repo / "verification"))
from v258_dirac_covariance_induction import induce

gamma = source["grading"]
duration = sum(source["durations"])
identity = source["identity"]
left = source["charged_left"][0]
right = source["charged_right"][0]


def maximum(a):
    return float(np.max(np.abs(a)))


def read(matrix):
    values, vectors = eigh(matrix)
    assert values.min() > 0
    covariance = (vectors * (values / (1 + values))) @ vectors.conj().T
    full_log = -(vectors * np.log(values)) @ vectors.conj().T
    implementation_error = maximum(induce(covariance) - full_log)
    assert implementation_error < 1e-11
    # Per existing strip time. The unknown common conversion mu_geo is not fitted.
    odd = (full_log - gamma @ full_log @ gamma) / (2 * duration)
    return covariance, full_log / duration, odd, implementation_error


def snapshot(odd, full, midpoint, midpoint_odd):
    singular = np.abs(np.linalg.eigvalsh(odd))
    singular_full = np.abs(np.linalg.eigvalsh(full))
    singular_mid = np.abs(np.linalg.eigvalsh(midpoint))
    singular_mid_odd = np.abs(np.linalg.eigvalsh(midpoint_odd))
    # This shape invariant ignores any positive common time or length scaling.
    def shape(values):
        return float(np.sum(values**4) / np.sum(values**2)**2)

    return {
        "ordered_odd_norm": float(np.linalg.norm(odd, 2)),
        "ordered_odd_smallest_singular": float(singular.min()),
        "ordered_odd_spectral_shape": shape(singular),
        "ordered_full_smallest_singular": float(singular_full.min()),
        "ordered_full_spectral_shape": shape(singular_full),
        "midpoint_full_smallest_singular": float(singular_mid.min()),
        "midpoint_full_spectral_shape": shape(singular_mid),
        "midpoint_odd_spectral_shape": shape(singular_mid_odd),
    }


def scalar_pair(z):
    return [float(z.real), float(z.imag)]


rows = {}
ordered_odds = []
for metric in (-.001, .001):
    transfer, sea = source["direct_whole"]("A", metric)
    reference, reference_sea = source["direct_whole"]("S", metric)
    covariance, full, odd, implementation_error = read(transfer)
    _, _, reference_odd, _ = read(reference)
    mid, phi_mid, p = source["generator"](2, "A", metric)
    h_mid_odd = -1j * gamma @ phi_mid
    assert maximum((mid - gamma @ mid @ gamma) / 2 - h_mid_odd) < 1e-13
    assert maximum(reference_odd) < 1e-13
    # The original reflected product is K* K. Reconstruct it independently
    # from the existing first half, with the original unchanged time letters.
    wing = identity.copy()
    for i, dt in enumerate(source["durations"][:3]):
        h, _, _ = source["generator"](i, "A", metric)
        e, v = eigh(h)
        wing = ((v * np.exp(-dt * e)) @ v.conj().T) @ wing
    split = wing.conj().T @ wing
    _, _, odd_split, _ = read(split)
    assert maximum(transfer - split) < 1e-11
    assert maximum(odd - odd_split) < 1e-11
    spectrum, eigenvectors = eigh(split)
    polar = wing @ ((eigenvectors / np.sqrt(spectrum)) @ eigenvectors.conj().T)
    assert maximum(polar.conj().T @ polar - identity) < 1e-11
    _, middle_full, middle_fixed_odd, _ = read(wing @ wing.conj().T)
    assert maximum(middle_full - polar @ full @ polar.conj().T) < 1e-11
    transported_grading = polar @ gamma @ polar.conj().T
    middle_transported_odd = (middle_full - transported_grading @ middle_full @ transported_grading) / 2
    assert maximum(middle_transported_odd - polar @ odd @ polar.conj().T) < 1e-11
    fixed_grading_error = float(np.linalg.norm(middle_fixed_odd - middle_transported_odd, 2))
    assert fixed_grading_error > .01
    relative = (sea-reference_sea+np.linalg.slogdet(identity+transfer)[1]
                - np.linalg.slogdet(identity+reference)[1])
    rows[str(metric)] = {
        **snapshot(odd, full, mid, h_mid_odd),
        "log_radius": metric,
        "radius": float(np.exp(metric)),
        "relative_empty_A_over_S": float(np.exp(relative)),
        "normalized_untwisted_L0_R0": scalar_pair(covariance[left, right]),
        "relative_weighted_untwisted_L0_R0": scalar_pair(
            np.exp(relative)*covariance[left, right]),
        "ordered_odd_L0_R0": scalar_pair(odd[left, right]),
        "midpoint_odd_L0_R0": scalar_pair(h_mid_odd[left, right]),
        "midpoint_to_ordered_odd_operator_difference": float(np.linalg.norm(odd-h_mid_odd, 2)),
        "polar_transport_grading_change": float(np.linalg.norm(transported_grading-gamma, 2)),
        "middle_odd_reader_error_if_grading_not_transported": fixed_grading_error,
        "middle_odd_reader_covariance_error": maximum(middle_transported_odd-polar@odd@polar.conj().T),
        "whole_vs_reflected_split_error": maximum(odd-odd_split),
        "paper_reader_vs_direct_log_error": implementation_error,
        "S_reference_odd_norm": maximum(reference_odd),
    }
    ordered_odds.append(odd)

# No new field, state or weight. Reversing only the tangential term is a
# mathematical check of the exact symmetry, not a third physical metric.
def reflected_transfer_with_p(p):
    wing = identity.copy()
    for i, dt in enumerate(source["durations"][:3]):
        _, phi, _ = source["generator"](i, "A", 0)
        e, v = eigh(p * gamma - 1j * gamma @ phi)
        wing = ((v * np.exp(-dt * e)) @ v.conj().T) @ wing
    return wing.conj().T @ wing


positive = reflected_transfer_with_p(.25)
negative = reflected_transfer_with_p(-.25)
parity_error = maximum(negative - gamma @ np.linalg.solve(positive, gamma))
_, positive_full, positive_odd, _ = read(positive)
_, negative_full, negative_odd, _ = read(negative)
assert parity_error < 1e-11
assert maximum(positive_odd-negative_odd) < 1e-11
assert maximum(negative_full+gamma@positive_full@gamma) < 1e-11
difference = float(np.linalg.norm(ordered_odds[1] - ordered_odds[0], 2))
assert difference > 1e-7
assert rows["-0.001"]["midpoint_odd_L0_R0"] == rows["0.001"]["midpoint_odd_L0_R0"]
assert all(abs(row["ordered_odd_spectral_shape"]-row["midpoint_odd_spectral_shape"]) > .01
           for row in rows.values())
print(json.dumps({
    "scope": "native fixed R032 chart, no full P_adm or original dynamical measure",
    "duration": duration,
    "same_original_Higgs_history_and_reference_S": rows,
    "ordered_odd_change_across_metrics_operator_norm": difference,
    "opposite_p_reflection_identity_max_error": parity_error,
    "opposite_p_ordered_odd_max_error": maximum(positive_odd-negative_odd),
    "actual_Q_Theta_and_D_adm_applied": False,
}, indent=2))
