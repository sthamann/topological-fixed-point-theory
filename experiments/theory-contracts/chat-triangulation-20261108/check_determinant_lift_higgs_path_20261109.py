"""Test the determinant-sign gauge lift on the existing R032 Higgs histories."""
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
path = directory / "check_original_whole_two_strips_20261109.py"
argv = sys.argv
try:
    sys.argv = [str(path), str(repo)]
    with contextlib.redirect_stdout(io.StringIO()):
        original = runpy.run_path(str(path))
finally:
    sys.argv = argv

native = original["native"]
data = original["source"]
gamma = original["grading"]
identity = original["identity"]
swap = data["swap"]
history = [np.array(h) for h in original["history"]]
durations = original["durations"]
majorana = original["majorana"]
inner_form = data["inner_form"]
weak_identity = np.eye(2)


def norm(a):
    return float(np.linalg.norm(a, 2))


def physical_lift(weak):
    # Higgs doublet transforms by lambda*q. The matched central color choice
    # at lambda=i is the determinant-one 3+2 carrier element with U3=-I3.
    rep = native.rep(1j, -1j * weak, 1j * np.eye(3))
    return rep @ swap @ rep.conj() @ swap


def fixed_higgs_lift(h):
    return 2 * np.outer(h, h.conj()) / np.vdot(h, h).real - weak_identity


weak = [fixed_higgs_lift(h) for h in history]
lifts = [physical_lift(a) for a in weak]
for h, a, lift in zip(history, weak, lifts):
    assert norm(a @ a - weak_identity) < 1e-14
    assert abs(np.linalg.det(a) + 1) < 1e-14
    assert np.linalg.norm(a @ h - h) < 1e-14
    assert norm(lift @ lift - identity) < 1e-14
    assert norm(lift - lift.conj().T) < 1e-14
    assert norm(lift @ majorana - majorana @ lift) < 1e-14
    assert norm(lift @ gamma - gamma @ lift) < 1e-14

# The existing undressed links are identities. Use the actual dressed
# links independently to check the same transport failure in another chart.
frames = original["frames"]
dressed_links = original["links"]
seams = []
for i, (a, b) in enumerate(zip(weak, weak[1:])):
    before, after = lifts[i], lifts[i+1]
    dress_before = frames[i].conj().T @ before @ frames[i]
    dress_after = frames[i+1].conj().T @ after @ frames[i+1]
    link = dressed_links[i]
    squared_overlap = (abs(np.vdot(history[i], history[i+1]))**2
                       / (np.vdot(history[i], history[i]).real
                          * np.vdot(history[i+1], history[i+1]).real))
    undressed_defect = norm(after - before)
    dressed_defect = norm(dress_after @ link - link @ dress_before)
    angle_formula = 2 * np.sqrt(max(0, 1 - squared_overlap))
    assert abs(undressed_defect - angle_formula) < 1e-12
    assert abs(dressed_defect - undressed_defect) < 1e-12
    seams.append({"time": sum(durations[:i+1]), "weak_line_overlap_squared": float(squared_overlap),
                  "parallel_intertwining_defect": undressed_defect})

# Candidate global order-two transformation, held at the actual first
# endpoint. It transforms the Higgs histories themselves and is not being
# identified with the separate determinant-sector operator of the paper.
anchor = lifts[0]
anchor_weak = weak[0]
idx_r = native.IDX[("e", "R", 0, 0)]
idx_l = native.IDX[("e", "L", 0, 0)]
assert norm((anchor + identity)[:, idx_r:idx_r+1]) < 1e-14
neutral = physical_lift(1j * weak_identity)
assert abs(np.linalg.det(1j * weak_identity) + 1) < 1e-14
assert abs((neutral @ neutral)[idx_l, idx_l] + 1) < 1e-14
assert abs((neutral @ neutral)[idx_r, idx_r] - 1) < 1e-14

# The involutive weak representatives form the conjugacy class of a
# traceless reflection. Its class mean on this representation can be
# evaluated exactly by the six opposite Pauli directions: all doublet
# terms are linear in the weak matrix and average to zero.
pauli = (np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]),
         np.diag([1, -1]))
class_mean = sum(physical_lift(s * a) for a in pauli for s in (1, -1)) / 6
assert norm(class_mean[:, idx_l:idx_l+1]) < 1e-14
assert abs(class_mean[idx_r, idx_r] + 1) < 1e-14
assert norm(class_mean @ class_mean - identity) > .9


def strip(fields, p, times=durations):
    wing = identity.copy()
    sea = 0.0
    for h, dt in zip(fields, times):
        op = p * gamma - 1j * gamma @ (inner_form(h) + majorana)
        e, v = eigh(op)
        wing = ((v * np.exp(-dt * e)) @ v.conj().T) @ wing
        sea += dt * e[:len(e)//2].sum()
    return wing, sea


metrics = {}
plus = (identity + anchor) / 2
minus = (identity - anchor) / 2
for m in (-.001, .001):
    p = .25 * np.exp(-m)
    local_commutators = []
    anchored_commutators = []
    for i in range(len(history)):
        operator, _, _ = original["generator"](i, "A", m)
        local_commutators.append(norm(lifts[i] @ operator - operator @ lifts[i]))
        anchored_commutators.append(norm(anchor @ operator - operator @ anchor))
    assert max(local_commutators) < 1e-12
    assert max(anchored_commutators) > 1
    original_strip, sea = strip(history, p)
    transformed_strip, transformed_sea = strip([anchor_weak @ h for h in history], p)
    covariance_error = norm(anchor @ original_strip @ anchor.conj().T - transformed_strip)
    assert covariance_error < 1e-11
    assert abs(sea - transformed_sea) < 1e-11
    original_trace = sea + np.linalg.slogdet(identity + original_strip)[1]
    transformed_trace = transformed_sea + np.linalg.slogdet(identity + transformed_strip)[1]
    half, _ = strip(history[:3], p, durations[:3])
    assert norm(half.conj().T @ half - original_strip) < 1e-11
    averaged_empty = (original_strip + transformed_strip) / 2
    averaged_empty_commutator = norm(averaged_empty @ anchor - anchor @ averaged_empty)
    empty_leakage = norm(minus @ half @ plus)
    charged_leakage = float(np.linalg.norm((plus @ half)[:, idx_r]))
    assert averaged_empty_commutator < 1e-11
    assert min(empty_leakage, charged_leakage) > .001
    # These are lower bounds on the corresponding trace diagnostics on the
    # no-particle/one-particle portions of the actual exterior representation.
    # The reflected gauge-related half has exactly the same norms.
    transformed_half = anchor @ half @ anchor
    assert abs(empty_leakage-norm(minus @ transformed_half @ plus)) < 1e-12
    assert abs(charged_leakage-np.linalg.norm((plus @ transformed_half)[:, idx_r])) < 1e-12
    metrics[str(m)] = {
        "max_local_generator_commutator": max(local_commutators),
        "max_parallel_anchor_generator_commutator": max(anchored_commutators),
        "whole_transfer_covariance_error": covariance_error,
        "fermion_log_weight_gauge_pair_error": abs(original_trace - transformed_trace),
        "gauge_pair_averaged_whole_empty_commutator": averaged_empty_commutator,
        "half_leakage_from_even_to_odd_norm": empty_leakage,
        "charged_eR_half_leakage_into_even_norm": charged_leakage,
    }

def pair(z):
    return [float(z.real), float(z.imag)]


print(json.dumps({
    "scope": "actual R032 fixed field data and original v252 gauge action; determinant gauge candidate, not QTheta identification",
    "link_checks": seams,
    "two_metric_checks": metrics,
    "anchored_global_lift_eR": pair(anchor[idx_r, idx_r]),
    "order_four_weak_central_lift_eL_squared": pair((neutral@neutral)[idx_l, idx_l]),
    "order_four_weak_central_lift_eR_squared": pair((neutral@neutral)[idx_r, idx_r]),
    "orientation_averaged_gauge_lift_eL": pair(class_mean[idx_l, idx_l]),
    "orientation_averaged_gauge_lift_eR": pair(class_mean[idx_r, idx_r]),
    "orientation_averaged_gauge_lift_is_not_involution_norm": norm(class_mean@class_mean-identity),
    "full_physical_source_and_QTheta_gap_evaluated": False,
}, indent=2))
