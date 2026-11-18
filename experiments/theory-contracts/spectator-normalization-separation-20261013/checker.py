#!/usr/bin/env python3
"""Diagnostic of the existing v1020 model, not a new physical source."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = ROOT / "verification/v1020_spectator_sector_obstruction.py"
sys.path.insert(0, str(ROOT / "verification"))
spec = importlib.util.spec_from_file_location("original_v1020", SOURCE)
original = importlib.util.module_from_spec(spec)
spec.loader.exec_module(original)


def complex_commutant_dimension(h, p):
    """Complex nullity of X -> ([X,h], [X,p]); column vectorization."""
    n = len(h)
    identity = np.eye(n)
    matrix = np.vstack([
        np.kron(h.T, identity) - np.kron(identity, h),
        np.kron(p.T, identity) - np.kron(identity, p),
    ])
    singular = np.linalg.svd(matrix, compute_uv=False)
    return n * n - int(np.sum(singular > 1e-9 * singular[0]))


def thermal(h, beta):
    values, vectors = np.linalg.eigh(h)
    weights = np.exp(-beta * values)
    z = float(weights.sum())
    return (vectors * (weights / z)) @ vectors.conj().T, z


def run():
    # These are diagnostic choices, not TFPT constants: a regular momentum,
    # the documented exceptional momentum, and a nonzero scalar/basis change.
    beta, shift = 0.7, 0.4
    samples = []
    for momentum in (0.0, 0.37):
        h = original.strip_hamiltonian(momentum, original.MASS_B)
        n = len(h)
        p = original.seam_projector(n)
        u = np.diag(np.exp(1j * 0.217 * np.arange(n)))
        changed = u @ h @ u.conj().T
        changed_p = u @ p @ u.conj().T
        translated = h + shift * np.eye(n)
        rho, z = thermal(h, beta)
        rho_shifted, z_shifted = thermal(translated, beta)
        old = original.commutant_dimension(h, p, n)
        old_changed = original.commutant_dimension(changed, changed_p, n)
        cc = complex_commutant_dimension(h, p)
        cc_changed = complex_commutant_dimension(changed, changed_p)
        cc_shifted = complex_commutant_dimension(translated, p)
        kd = original.krylov_dimension(h, original.SEAM_DIMS, n)
        kd_changed = original.krylov_dimension(changed, original.SEAM_DIMS, n)
        kd_shifted = original.krylov_dimension(translated, original.SEAM_DIMS, n)
        assert old == cc == cc_changed == cc_shifted
        assert kd == kd_changed == kd_shifted
        assert np.max(np.abs(rho - rho_shifted)) < 1e-12
        assert abs(z_shifted / z - np.exp(-beta * shift)) < 1e-12
        samples.append({
            "momentum": momentum, "dimension": n,
            "commutant_original_helper": old,
            "commutant_original_helper_after_basis_change": old_changed,
            "commutant_complex": cc, "commutant_complex_after_basis_change": cc_changed,
            "commutant_complex_after_scalar_shift": cc_shifted,
            "krylov": kd, "krylov_after_basis_change": kd_changed,
            "krylov_after_scalar_shift": kd_shifted,
            "normalized_thermal_max_error": float(np.max(np.abs(rho - rho_shifted))),
            "partition_ratio": z_shifted / z,
            "expected_partition_ratio": float(np.exp(-beta * shift)),
            "ground_energy_shift": float(np.linalg.eigvalsh(translated)[0] - np.linalg.eigvalsh(h)[0]),
        })
    assert samples[0]["commutant_original_helper_after_basis_change"] != samples[0]["commutant_complex"]
    assert samples[1]["krylov"] == samples[1]["dimension"]
    result = {
        "verdict": "PARTIAL", "assertions_passed": True,
        "source": str(SOURCE.relative_to(ROOT)),
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "parameters": {"beta": beta, "scalar_shift": shift, "basis_phase_step": 0.217},
        "samples": samples,
        "scope": "Finite original v1020 one-particle Hamiltonian. Thermal Gibbs normalization at positive beta; no preservation claim for v1020's zero-chemical-potential Fermi projector. Scalar-shift admissibility in the full original TFPT source is not established.",
        "finding": "Seam cyclicity and joint commutant do not determine absolute normalization, even for the seam-cyclic sample. Original helper's general complex-domain annotation fails under a complex unitary basis change at the exceptional momentum; its original real-input values are not refuted.",
    }
    (HERE / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    run()
