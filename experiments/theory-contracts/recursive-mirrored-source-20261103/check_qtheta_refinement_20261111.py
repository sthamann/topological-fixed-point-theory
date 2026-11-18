"""Compare the installed Q assembler with a nested spatial Fourier refinement."""
from copy import deepcopy
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np
from scipy.linalg import eigvalsh, solve


def main(repository):
    package = repository / (
        "experiments/theory-contracts/original-paper-source-reconciliation-20261102"
        "/contribution_20261106Q"
    )
    sys.path.insert(0, str(package / "runtime"))
    from assemble import Assembler
    from source_adapter import evaluate_contribution, evaluate_source

    path = package / "domain/native-diagnostic-v1.json"
    coarse = Assembler(path)
    genes = [1, 1, 1]
    original = coarse.assemble({"genes": genes})
    whole = evaluate_source({}, original, coarse.domain)
    fiber = evaluate_contribution({}, original, coarse.domain)
    assert all(value is None for value in whole["values"].values())
    assert fiber["original_measure_mass"] is None

    amplitude, harmonic, unit = 0.13, 8, "13/100"
    fine_domain = deepcopy(json.loads(path.read_text()))
    fine_domain["modes"] = list(range(-7, 7))
    fine_domain["field_proposal"]["slots"].append(
        {"kind": "cos", "harmonic": harmonic, "component": 0, "unit": unit}
    )
    with TemporaryDirectory() as directory:
        fine_path = Path(directory) / "refinement.json"
        fine_path.write_text(json.dumps(fine_domain))
        fine = Assembler(fine_path)
        zero = fine.assemble({"genes": genes + [1]})["h_N"]
        moved = fine.assemble({"genes": genes + [2]})["h_N"]
        opposite = fine.assemble({"genes": genes + [0]})["h_N"]

    modes = fine_domain["modes"]
    indices = [96*modes.index(n) + i for n in coarse.domain["modes"] for i in range(96)]
    h0 = original["h_N"]
    assert np.array_equal(zero[np.ix_(indices, indices)], h0)
    assert np.array_equal(moved[np.ix_(indices, indices)], h0)
    assert not np.array_equal(zero, moved)
    perturbation = -1j*fine.gamma @ fine.V[0]
    norm_square = np.vdot(perturbation, perturbation).real
    square_trace = np.vdot(moved, moved).real - np.vdot(zero, zero).real
    expected_trace = (len(modes)-harmonic)*amplitude**2*norm_square/2
    assert abs(square_trace-expected_trace) < 1e-9
    returned = moved[:, indices].conj().T @ moved[:, indices] - h0 @ h0
    expected_return = (amplitude/2)**2*(perturbation @ perturbation)
    for i, n in enumerate(coarse.domain["modes"]):
        block = slice(96*i, 96*(i+1))
        expected = expected_return if n in (-2, 1) else 0
        np.testing.assert_allclose(returned[block, block], expected, rtol=0, atol=1e-13)
        remainder = returned[block].copy()
        remainder[:, block] = 0
        assert not np.any(remainder)

    # In this refinement each nonzero harmonic couples just two Fourier
    # modes, so the full trace is the sum over these invariant blocks.
    groups = [[n for n in modes if n % harmonic == r] for r in range(harmonic)]
    signs = np.ones(len(zero))
    for group in groups:
        assert len(set(group) & set(coarse.domain["modes"])) <= 1
        if len(group) == 2:
            partner = next(n for n in group if n not in coarse.domain["modes"])
            start = 96*modes.index(partner)
            signs[start:start+96] = -1
    assert np.array_equal(opposite, signs[:, None]*moved*signs[None, :])
    assert np.all(signs[indices] == 1)
    group_indices = [
        [96*modes.index(n)+i for n in group for i in range(96)] for group in groups
    ]
    proper_time = 0.001
    heat_difference = 0.0
    for ix in group_indices:
        old_energy = eigvalsh(zero[np.ix_(ix, ix)], check_finite=False)
        new_energy = eigvalsh(moved[np.ix_(ix, ix)], check_finite=False)
        heat_difference += np.sum(np.exp(-proper_time*new_energy**2)
                                  - np.exp(-proper_time*old_energy**2))
    assert heat_difference < -1e-4

    def block(n):
        i = 96*modes.index(n)
        return slice(i, i+96)

    left, right = block(-2), block(6)
    eye = np.eye(96)
    X = zero[left, left] - 1j*eye
    Y = zero[right, right] - 1j*eye
    B = moved[left, right]
    parent = moved[np.r_[left, right], :][:, np.r_[left, right]] - 1j*np.eye(192)
    direct = solve(parent, np.eye(192))[:96, :96]
    schur = solve(X-B@solve(Y, B.conj().T), eye)
    schur_error = np.linalg.norm(direct-schur)
    resolvent_change = np.linalg.norm(direct-solve(X, eye))
    assert schur_error < 1e-12 and resolvent_change > 1e-3

    # The exact first coefficients of the conditional scalar moment reader.
    pair_indices = np.r_[left, right]
    start_pair = zero[np.ix_(pair_indices, pair_indices)]
    active_pair = moved[np.ix_(pair_indices, pair_indices)]
    other_pair = opposite[np.ix_(pair_indices, pair_indices)]
    assert np.allclose((active_pair @ active_pair)[:96, :96]
                       - (start_pair @ start_pair)[:96, :96], expected_return)
    signed_mark = np.trace(perturbation @ active_pair[96:, :96]).real
    assert np.isclose(signed_mark, amplitude*norm_square/2)
    assert np.isclose(np.trace(perturbation @ other_pair[96:, :96]).real, -signed_mark)

    print(json.dumps({
        "unmodified_domain_id": original["domain_id"],
        "native_endpoint_missing": whole["missing_operations"],
        "native_fiber_missing": fiber["missing_operations"],
        "native_values": whole["values"],
        "native_measure_mass": fiber["original_measure_mass"],
        "coarse_exactly_equal": True,
        "fine_modes": modes,
        "new_spatial_harmonic": harmonic,
        "delta_fine_square_trace": square_trace,
        "coarse_return_square_trace": float(np.trace(returned).real),
        "delta_fine_heat_trace_at_s_0_001": heat_difference,
        "coarse_return_resolvent_frobenius_change_at_i": resolvent_change,
        "schur_cross_check_frobenius_error": schur_error,
        "opposite_sign_conjugacy_fixes_entire_coarse_space": True,
        "odd_cross_mode_mark_at_positive_amplitude": signed_mark,
        "physical_integration_computed": False,
    }, indent=2))


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
