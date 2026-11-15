"""Eliminate the native spatial full96 Fourier tails as matrix-valued returns."""
import json
from pathlib import Path
import subprocess
import sys

import sympy as sp

from higgs_moment_search import field_coefficients
from source import NATIVE, complex_pair

RUNTIME = NATIVE.parents[1] / "runtime"
INPUTS = (RUNTIME / "assemble.py", RUNTIME / "repo_paths.py", Path("verification/v252_full_finite_triple.py"),
          Path("experiments/theory-contracts/tfpt-unified-source-derivation-20261029/"
               "aps_source_review/check_original_higgs_fiber.py"))
GENES = [2, 1, 1]


def spatial_cut(z, diagonal, edge, layers):
    """Return the compressed resolvent and the hopping determinant relative to the diagonal."""
    import numpy as np

    size = edge.shape[0]
    identity = np.eye(size, dtype=complex)
    central = np.block([[z * identity - diagonal(-1), -edge],
                        [-edge.conj().T, z * identity - diagonal(0)]])
    phase, logabs = 1 + 0j, 0.0

    def include(block, reference):
        nonlocal phase, logabs
        sign, log = np.linalg.slogdet(block)
        reference_sign, reference_log = np.linalg.slogdet(reference)
        phase *= sign / reference_sign
        logabs += log - reference_log

    for first, direction, coupling, target in ((-1, -1, edge.conj().T, slice(0, size)),
                                               (0, 1, edge, slice(size, 2 * size))):
        contribution = np.zeros_like(identity)
        for distance in range(layers, 0, -1):
            bare = z * identity - diagonal(first + direction * distance)
            outside = bare - contribution
            include(outside, bare)
            contribution = coupling @ np.linalg.solve(outside, coupling.conj().T)
        central[target, target] -= contribution
    bare_retained = np.block([[z * identity - diagonal(-1), np.zeros_like(edge)],
                              [np.zeros_like(edge), z * identity - diagonal(0)]])
    include(central, bare_retained)
    return np.linalg.solve(central, np.eye(2 * size, dtype=complex)), phase, logabs


def analyze(domain):
    import numpy as np
    from assemble import Assembler

    assembler = Assembler(domain)
    fixture = json.loads(domain.read_text(), parse_float=str)
    coefficients = field_coefficients(fixture["field_proposal"], GENES)
    active = {mode: value for mode, value in coefficients.items() if value.is_zero_matrix is False}
    link = sp.Matrix([sp.Rational(13, 200), 0, 0, 0])
    if set(active) != {-1, 0, 1} or active[-1] != link or active[1] != link:
        raise ArithmeticError("native nearest-neighbour Higgs proposal changed")
    yukawas = [sp.Matrix([[complex_pair(entry) for entry in row] for row in fixture["couplings"][name]])
               for name in ("Yu", "Yd", "Ye", "Ynu")]
    # A unit quaternion rotates each left weak doublet. The graded lift
    # repeats each Yukawa singular value on the two chiral and two real-lift legs.
    determinants = [sp.expand(matrix.det()) for matrix in yukawas]
    if any(value == 0 for value in determinants):
        raise ArithmeticError("a native Yukawa block lost invertibility")
    family_norms = [sp.expand(link[0]**2 * sp.trace(matrix.H * matrix)) for matrix in yukawas]
    norm_squared = 4 * sum(multiplicity * value for multiplicity, value in zip((3, 3, 1, 1), family_norms))
    edge_norm_upper = sp.Rational(5, 16)
    if max(family_norms) >= edge_norm_upper**2:
        raise ArithmeticError("native Fourier edge no longer obeys the certified norm bound")
    eta = 2
    ratio = 2 * edge_norm_upper / eta
    z = 1j * eta
    candidate = {"schema": "tfpt.field-configuration.v1", "domain_sha256": assembler.domain["_sha256"],
                 "genes": GENES}
    model = assembler.assemble(candidate)
    h = model["h_N"]
    size = 96
    if model["modes"] != [-2, -1, 0, 1] or h.shape != (4 * size, 4 * size):
        raise ArithmeticError("unexpected assembly of the original four Fourier modes")
    gamma = model["gamma96"]
    potential = -1j * gamma @ model["phi_fourier"][0]
    edge = -1j * gamma @ model["phi_fourier"][-1]

    def diagonal(mode):
        return (mode + model["alpha"]) * gamma + potential

    expected = np.block([[diagonal(n) if n == m else edge if m == n + 1 else
                          edge.conj().T if m == n - 1 else np.zeros_like(edge)
                          for m in model["modes"]] for n in model["modes"]])
    np.testing.assert_allclose(h, expected, atol=1e-14, rtol=0)
    np.testing.assert_allclose(np.linalg.norm(edge, "fro")**2, float(norm_squared), atol=1e-14, rtol=0)
    full = z * np.eye(len(h)) - h
    direct = np.linalg.inv(full)[size:3 * size, size:3 * size]
    result, phase, logabs = spatial_cut(z, diagonal, edge, 1)
    np.testing.assert_allclose(result, direct, atol=1e-14, rtol=0)
    direct_phase, direct_logabs = np.linalg.slogdet(full)
    for mode in model["modes"]:
        bare_phase, bare_logabs = np.linalg.slogdet(z * np.eye(size) - diagonal(mode))
        direct_phase /= bare_phase
        direct_logabs -= bare_logabs
    np.testing.assert_allclose([phase, logabs], [direct_phase, direct_logabs], atol=1e-11, rtol=0)

    depths = (0, 1, 2, 4, 8, 16)
    reference, _, _ = spatial_cut(z, diagonal, edge, depths[-1])
    rows = []
    for depth in depths:
        approx, phase, logabs = spatial_cut(z, diagonal, edge, depth)
        bound = 2 * ratio**(2 * (depth + 1)) / (eta * (1 - ratio))
        rows.append({"eliminated_layers_per_side": depth, "krylov_rank_through_depth": 2 * size * (depth + 1),
                     "first_different_return_power": 2 * (depth + 1),
                     "infinite_resolvent_error_upper_bound": str(bound),
                     "numeric_max_entry_delta_to_depth_16": float(np.max(np.abs(approx - reference))),
                     "relative_spatial_determinant_phase": [float(phase.real), float(phase.imag)],
                     "relative_spatial_determinant_logabs": float(logabs)})
    return {"status": "EXACT_NATIVE_SPATIAL_MEMORY_OBSTRUCTION", "genes": GENES,
            "source_sha256": model["native_inputs_sha256"], "retained_modes": [-1, 0],
            "exact_edge_frobenius_squared": str(norm_squared), "exact_edge_rank": size,
            "edge_operator_norm_upper_bound": str(edge_norm_upper),
            "resolvent_z": f"{eta}*I", "depths": rows,
            "determinant_reference": "spatial operator with the same diagonal blocks and no Fourier hopping",
            "finite_static_exact_realization_exists": False, "physical_OS_time_identified": False,
            "absolute_weighted_R113_produced": False}


def replay(root):
    runtime = root / RUNTIME
    child = subprocess.run([sys.executable, "-B", str(Path(__file__).resolve()), str(root / NATIVE), str(runtime)],
                           cwd=root, capture_output=True, text=True, timeout=60, check=True)
    return json.loads(child.stdout)


if __name__ == "__main__":
    sys.path.insert(0, sys.argv[2])
    print(json.dumps(analyze(Path(sys.argv[1]))))
