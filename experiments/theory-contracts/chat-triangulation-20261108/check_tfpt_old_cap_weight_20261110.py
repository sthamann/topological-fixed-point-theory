"""Audit the old C57/C58 excitation cap against its own free reference."""

import contextlib
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
import sympy as sp


repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
sys.path.insert(0, str(repo / "verification"))
from v330_qgamb_admissible_measure import symmetric_transfer
from v987_os_dilation_package import B_SYM


pi = sp.Matrix([[sp.Rational(1, 2), sp.Rational(3, 10), sp.Rational(1, 5)]])
v240_diagonal_metric = sp.diag(*pi)
macro = B_SYM**6
assert np.max(abs(np.array(macro, dtype=float) - symmetric_transfer()[0])) < 1e-15
assert pi * B_SYM == sp.Matrix([[sp.Rational(19, 45), sp.Rational(13, 45), sp.Rational(13, 45)]])
v987_balance_defect = (v240_diagonal_metric * B_SYM - B_SYM.T * v240_diagonal_metric)[0, 1]
v330_balance_defect = (v240_diagonal_metric * macro - macro.T * v240_diagonal_metric)[0, 1]
assert v987_balance_defect == sp.Rational(1, 90)
assert v330_balance_defect == sp.Rational(1267, 21870)
uniform = sp.ones(1, 3) / 3
assert uniform * macro == uniform

sys.argv = ["check_r113_dressed_open_20261109.py", str(base)]
with contextlib.redirect_stdout(io.StringIO()):
    source = runpy.run_path(str(base / sys.argv[0]))
gamma = source["grading"]
majorana = source["majorana"]
inner_form = source["inner_form"]
unit_higgs = inner_form(np.array([0., 1.]))
assert np.max(abs(unit_higgs - source["yukawa"])) < 1e-14
assert np.max(abs(gamma @ unit_higgs + unit_higgs @ gamma)) < 1e-14
assert abs(np.trace(unit_higgs @ majorana)) < 1e-14
late = 1 / 8  # Existing C58 width, not a fit or new reference.
q = (3 + 1 / np.sqrt(2)) / 4
data = {}
for radius in (0., .55, 1., 4., 8., 16.):
    phi = inner_form(np.array([0., radius])) + majorana
    eigen = np.linalg.eigvalsh(phi)
    assert np.max(abs(eigen + eigen[::-1])) < 1e-12
    masses = eigen[len(eigen) // 2:]
    assert masses[-1] >= radius * np.linalg.norm(unit_higgs, 2) - np.linalg.norm(majorana, 2) - 1e-12
    results = {}
    for cutoff in (128, 512, 2048):
        p = abs(np.arange(-cutoff, cutoff) + .25)[:, None]
        energy = np.hypot(p, masses)
        vacancy = (1 - p / energy) / 2
        cap_terms = np.log1p(vacancy * np.expm1(-2 * late * energy))
        bulk = late * np.sum(masses**2 / (energy + p))
        cap = cap_terms.sum()
        selection = p[:, 0] <= masses[-1]
        if np.exp(-2 * late * masses[-1]) <= 1 / 2:
            assert np.max(cap_terms[selection, -1]) < np.log(q) + 1e-14
            assert cap <= selection.sum() * np.log(q) + 1e-12
        results[str(cutoff)] = {"log_excitation_cap": float(cap),
                                "log_free_relative_cap": float(cap + bulk),
                                "free_reference_bulk": float(bulk)}
    assert results["2048"]["log_excitation_cap"] <= results["128"]["log_excitation_cap"] + 1e-12
    log_ratio = np.log(4)
    observed_slope = ((results["2048"]["free_reference_bulk"] - results["512"]["free_reference_bulk"])
                      / log_ratio)
    predicted_slope = late * np.sum(masses**2)
    assert abs(observed_slope / predicted_slope - 1) < .03
    data[str(radius)] = dict(largest_mass=float(masses[-1]),
                             logarithmic_free_bulk_slope=float(predicted_slope),
                             observed_bulk_slope=float(observed_slope), cutoffs=results)

print(json.dumps({"scope": "C57/C58 constant R032 quarter-circle cap with the existing old Slater state, no physical R113 boson measure or selector",
                  "three_state_transfer": {"v987_six_step_equals_v330": True,
                                           "v240_sample_under_v987_B": list(map(str, pi * B_SYM)),
                                           "v987_detailed_balance_defect_in_v240_sample": str(v987_balance_defect),
                                           "v330_detailed_balance_defect_in_v240_sample": str(v330_balance_defect),
                                           "v330_uniform_state_stationary": True},
                  "quarter": .25, "cap_time": late,
                  "upper_bound_per_large_mass_mode": q,
                  "higgs_response_norm": float(np.linalg.norm(unit_higgs, 2)),
                  "majorana_norm": float(np.linalg.norm(majorana, 2)),
                  "higgs_magnitudes": data}, indent=2))
