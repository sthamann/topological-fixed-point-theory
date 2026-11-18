"""Separate the C58 circle state limit from the absolute cap scalar."""

import contextlib
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from scipy.special import zeta


repo = Path(sys.argv[1]).resolve()
contract = repo / "experiments/theory-contracts/chat-triangulation-20261108"
with contextlib.redirect_stdout(io.StringIO()):
    source = runpy.run_path(str(contract / "check_tfpt_conditioned_r113_weights_20261110.py"))
gamma = source["gamma"]
ell = .125


def singular_values(radius, sign):
    potential = source["clifford_majorana"] + radius * source["clifford_higgs"]
    np.testing.assert_allclose(gamma @ potential + potential @ gamma, 0, atol=1e-14)
    occupied = np.flatnonzero(sign * np.diag(gamma) < 0)
    empty = np.flatnonzero(sign * np.diag(gamma) > 0)
    return np.linalg.svd(potential[np.ix_(occupied, empty)], compute_uv=False)


def cap_data(momentum, singular):
    energy = np.hypot(momentum, singular)
    a = np.tanh(ell * energy)
    pairs = singular * a / (energy + momentum * a)
    # Stable evaluation of log(A/e^(ell p)) in each singular direction.
    ratio = momentum / energy
    relative_log = (ell * singular**2 / (energy + momentum)
                    + np.log1p((ratio-1)*(1-np.exp(-2*ell*energy))/2))
    return pairs, relative_log


result = {}
for radius in (0., .55):
    row = {}
    for cutoff in (2048, 4096):
        log_norm, missing_bound, log_cap, relative_majorana_cap = (0.,)*4
        coefficient = 0.
        for sign, shift in ((1, .25), (-1, .75)):
            momentum = (np.arange(cutoff) + shift)[:,None]
            singular = singular_values(radius, sign)
            majorana = singular_values(0., sign)
            pairs, relative_log = cap_data(momentum, singular)
            assert np.max(pairs-singular/(2*momentum)) < 1e-14
            log_norm += np.log1p(pairs**2).sum()
            missing_bound += (singular**2).sum()/4*zeta(2, cutoff+shift)
            log_cap += relative_log.sum()
            relative_majorana_cap += (relative_log-cap_data(momentum, majorana)[1]).sum()
            coefficient += ell*((singular**2).sum()-(majorana**2).sum())/2
        row[str(cutoff)] = {
            "old_vacuum_squared_overlap_lower": float(np.exp(-log_norm-missing_bound)),
            "old_vacuum_squared_overlap_upper": float(np.exp(-log_norm)),
            "missing_log_norm_upper_bound": float(missing_bound),
            "free_relative_log_cap": float(log_cap),
            "majorana_relative_log_cap": float(relative_majorana_cap),
            "majorana_relative_log_cutoff_coefficient": float(coefficient),
        }
    lo, hi = row["2048"], row["4096"]
    assert lo["old_vacuum_squared_overlap_lower"] <= hi["old_vacuum_squared_overlap_lower"]
    assert lo["old_vacuum_squared_overlap_upper"] >= hi["old_vacuum_squared_overlap_upper"]
    doubling = hi["majorana_relative_log_cap"] - lo["majorana_relative_log_cap"]
    assert abs(doubling-hi["majorana_relative_log_cutoff_coefficient"]*np.log(2)) < .01
    row["majorana_relative_log_cap_cutoff_doubling"] = float(doubling)
    result[str(radius)] = row

print(json.dumps({"ell":ell, "quarter":.25, "circle_radius":1,
                  "scope":"Old C58, fixed Higgs and Majorana on the quarter circle, using the old occupied free sea. No dynamical boson integration, actual R113 reference, or projected q inferred.",
                  "results":result}, indent=2))
