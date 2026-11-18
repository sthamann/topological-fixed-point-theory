"""Evaluate the native Q spatial factor and inspect the original source call."""
from itertools import product
import json
from pathlib import Path
import sys

import numpy as np


def main(repo):
    package = repo / (
        "experiments/theory-contracts/original-paper-source-reconciliation-20261102"
        "/contribution_20261106Q"
    )
    sys.path.insert(0, str(package / "runtime"))
    from assemble import Assembler
    from original_source_adapter import capability_report, integration_plan

    native = Assembler(package / "domain/native-diagnostic-v1.json")
    profiles = {
        genes: native.assemble({"genes": list(genes)})
        for genes in product(range(3), repeat=3)
    }
    # A Fourier coefficient of H*H is a convolution of H modes.
    norm_harmonics = set()
    for profile in profiles.values():
        field = {
            n: v[[0, 2]] + 1j * v[[1, 3]]
            for n, v in profile["higgs_fourier"].items()
        }
        for k in range(17):
            coefficient = sum(
                np.vdot(v, field[n+k]) for n, v in field.items() if n+k in field
            )
            if coefficient != 0:
                norm_harmonics.add(k)
    assert max(norm_harmonics) == 4

    # Vary the actual first and third proposal slots. Keep the second at +1.
    operators = {
        (i-1, j-1): profiles[i, 2, j]["h_N"]
        for i, j in product(range(3), repeat=2)
    }
    h = operators[0, 0]
    b = (operators[1, 0] - operators[-1, 0]) / 2
    c = (operators[0, 1] - operators[0, -1]) / 2
    for (i, j), operator in operators.items():
        np.testing.assert_allclose(operator, h+i*b+j*c, atol=1e-14, rtol=0)

    bc = b @ c
    overlap = float(np.vdot(bc, bc).real)
    coefficient = 4*overlap + 2*np.einsum("ij,ji->", bc, bc).real
    assert overlap > 0 and coefficient >= 2*overlap
    weights = {-1: 1, 0: -2, 1: 1}

    def cross_difference(values):
        return sum(weights[i] * weights[j] * value
                   for (i, j), value in values.items())

    eigenvalues = {key: np.linalg.eigvalsh(value) for key, value in operators.items()}
    fourth_difference = cross_difference(
        {key: np.sum(e**4) for key, e in eigenvalues.items()}
    )
    np.testing.assert_allclose(fourth_difference, 4*coefficient, atol=1e-9, rtol=0)

    heat_difference = {}
    for t in (0.02, 0.01, 0.001, 1.0):
        # Remove the constant and linear Taylor terms before summing.
        heat = {key: np.sum(np.expm1(-t*e**2) + t*e**2)
                for key, e in eigenvalues.items()}
        heat_difference[str(t)] = float(cross_difference(heat))
    assert all(heat_difference[t] > 0 for t in ("0.02", "0.01", "0.001"))
    direct_heat = cross_difference({key: np.sum(np.exp(-e**2))
                                    for key, e in eigenvalues.items()})
    np.testing.assert_allclose(direct_heat, heat_difference["1.0"], atol=1e-10, rtol=0)
    assert direct_heat < -0.0008
    assert abs(heat_difference["0.001"] / 0.001**2 - 2*coefficient) < 0.05*coefficient

    status = capability_report(native.domain)
    assert not status["whole_source_bundle_supplied"]
    plan = integration_plan(native.domain)
    assert not plan["ready"] and plan["missing_operations"] == ["original_blocked_fiber_packet"]
    print(json.dumps({
        "native_domain_id": native.domain["domain_id"],
        "native_profiles": len(profiles),
        "highest_native_higgs_norm_harmonic": max(norm_harmonics),
        "U7_and_U16_zero_on_proposal_catalogue": True,
        "spatial_operator_size": len(h),
        "trace_h4_coefficient_x_squared_y_squared": coefficient,
        "trace_h4_mixed_difference": float(fourth_difference),
        "heat_mixed_difference_by_proper_time": heat_difference,
        "whole_source_bundle_supplied": status["whole_source_bundle_supplied"],
        "missing_direct_operation": status["missing_operations"],
        "missing_fibre_operation": plan["missing_operations"],
    }, indent=2))


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
