"""Check the installed source gate and the exact sector algebra for its requested word."""
import hashlib
import json
import sys
from pathlib import Path

import sympy as sp

repo = Path(sys.argv[1]).resolve()
caller = (repo / "experiments/theory-contracts"
          / "original-paper-source-reconciliation-20261102/contribution_20261106Q")
sys.path.insert(0, str(caller / "runtime"))

from assemble import load_domain
from source_adapter import evaluate_contribution, evaluate_source, whole_capability_report

path = caller / "domain/native-diagnostic-v1.json"
domain, _, _ = load_domain(path)
assembled = {
    "domain_id": domain["domain_id"],
    "native_inputs_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    "h_N_role": "not_used_for_original_endpoint_probe",
}
whole = evaluate_source({}, assembled, domain)
fiber = evaluate_contribution({"schema": "tfpt.original-source-unit.v1",
                               "original_configuration_id": None}, assembled, domain)
capability = whole_capability_report(domain)
assert capability["missing_operations"] == ["evaluate_original_full_blocked_endpoint_source"]
assert whole["status"] == "original_full_source_operation_not_supplied"
assert not capability["numeric_bundle_ready"]
assert all(value is None for value in whole["values"].values())
assert fiber["missing_operations"] == ["original_blocked_fiber_packet"]
assert fiber["original_measure_mass"] is None

# Pure projector algebra in the four formal nonoverlapping sectors, no physical
# matrices, weights, state or gap are assigned. Columns are P, Pstr-P, Eplus-Pstr, Eminus.
readouts = sp.Matrix([[1, 1, 1, 1], [1, 1, 1, -1],
                      [1, 1, 0, 0], [1, 0, 0, 0]])
assert readouts.inv() == sp.Matrix([[0, 0, 0, 1], [0, 0, 1, -1],
                                   [sp.Rational(1, 2), sp.Rational(1, 2), -1, 0],
                                   [sp.Rational(1, 2), -sp.Rational(1, 2), 0, 0]])

# Actual printed zero-parameter C6 pole and all 18 transport singular values.
# This only evaluates D_f's published comparison scale, never D_adm or H_OS.
z = (794 - 7 * sp.sqrt(9961)) / 2187
delta = z ** sp.Rational(1, 6)
assert sp.simplify(1594323*z**2 - 1157652*z + 47449) == 0
assert (sp.Rational(1, 2) < delta) and (delta < sp.Rational(2, 3))
cusps = [sp.Integer(1), sp.Rational(1, 3), sp.Rational(2, 3)]
gap_squared = {(str(y), n): y**2 + delta**2 - 2*y*delta*sp.cos(n*sp.pi/3)
               for y in cusps for n in range(6)}
first = min(gap_squared, key=lambda key: gap_squared[key].evalf(40))
assert first == ("2/3", 0)
epsilon_f = sp.Rational(2, 3) - delta
assert sp.simplify(gap_squared[first] - epsilon_f**2) == 0

print(json.dumps({
    "native_endpoint_call": {
        "domain_id": domain["domain_id"], "domain_sha256": assembled["native_inputs_sha256"],
        "whole_source_bundle": domain["whole_source_bundle"],
        "status": whole["status"], "missing_operations": whole["missing_operations"],
        "registered_whole_values": whole["values"],
        "fiber_status": fiber["status"], "fiber_missing_operations": fiber["missing_operations"],
        "original_measure_mass": fiber["original_measure_mass"],
        "evaluated_physical_metric_pairs": 0,
        "qtheta_marked_charged_word": None, "same_cap_empty_value": None,
    },
    "exact_sector_algebra": {
        "readout_order": ["I", "QTheta", "Pstr", "Padm"],
        "sector_order": ["retained", "low_gap_in_Pstr", "other_structural_reject_in_Qplus", "Qminus"],
        "rank": readouts.rank(),
        "twice_inverse": (2 * readouts.inv()).tolist(),
    },
    "original_C6_gap_comparison_only": {
        "delta_ph": str(delta.evalf(20)), "minimizing_cusp_and_phase_index": first,
        "epsilon_f": str(epsilon_f.evalf(20)),
        "half_epsilon_f": str((epsilon_f/2).evalf(20)),
        "epsilon_f_squared_over_four": str((epsilon_f**2/4).evalf(20)),
        "computed_comparison_singular_values": len(gap_squared),
        "identification_with_original_Dadm_established": False,
    },
}, indent=2, default=int))
