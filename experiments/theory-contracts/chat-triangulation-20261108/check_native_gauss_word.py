"""Placement test using the native v989 Hamiltonian and both interior stars."""

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply


source = Path(sys.argv[1]) / "verification/v989_lattice_gate_battery.py"
sys.path.insert(0, str(source.parent))
spec = importlib.util.spec_from_file_location("v989", source)
native = importlib.util.module_from_spec(spec)
spec.loader.exec_module(native)

psi, number = native._build_ham_ops()
h = csr_matrix(native.build_H(psi, number))
g1, g2 = (csr_matrix(native.gauss_unitary(number, x)) for x in (1, 2))
one = csr_matrix(np.eye(native.D_H))


def group_average(generator):
    return (one + generator + generator**2 + generator**3) / 4


def max_entry(matrix):
    return float(np.max(np.abs(matrix.data), initial=0))


p = group_average(g1) @ group_average(g2)
v = csr_matrix(native._link_op(1, native.Z1_H**2))
checks = {
    "projector": max_entry(p @ p - p),
    "projector_adjoint": max_entry(p - p.getH()),
    "hamiltonian_hermitian": max_entry(h - h.getH()),
    "time_preserves_sector": max_entry(h @ p - p @ h),
    "time_commutes_with_first_star": max_entry(h @ g1 - g1 @ h),
    "time_commutes_with_second_star": max_entry(h @ g2 - g2 @ h),
    "link_first_star_character_minus": max_entry(g1 @ v + v @ g1),
    "link_second_star_character_minus": max_entry(g2 @ v + v @ g2),
    "single_link_on_projected_space": max_entry(p @ v @ p),
}
assert all(value == 0 for value in checks.values()), checks
assert max_entry(v.getH() @ v - one) == 0

# Original projected Gibbs prescription, evaluated on its whole support.
eigenvalues, vectors = np.linalg.eigh(p.toarray())
basis = vectors[:, eigenvalues > 0.5]
assert basis.shape == (native.D_H, 64)
energies, rotation = np.linalg.eigh(basis.conj().T @ (h @ basis))
physical = basis @ rotation
beta = 1.0
weights = np.exp(-beta * (energies - energies[0]))
weights /= weights.sum()
inserted = v @ physical


def trace_on_physical(images):
    return complex(np.dot(weights, np.einsum("ij,ij->j", physical.conj(), images)))


rows = []
for tau in (0, 0.4, 1.0):
    transfer = -(h - energies[0] * one) * tau
    moved = expm_multiply(transfer, inserted)
    closed = p @ (v.getH() @ moved)
    amplitudes = {
        "closed": trace_on_physical(closed),
        "phase_inside": trace_on_physical(p @ (v.getH() @ (g1 @ moved))),
        "phase_after_return": trace_on_physical(p @ (g1 @ closed)),
        "projection_inside": trace_on_physical(p @ (v.getH() @ (p @ moved))),
    }
    assert amplitudes["closed"].real > 0
    assert abs(amplitudes["phase_inside"] + amplitudes["closed"]) < 1e-12
    assert abs(amplitudes["phase_after_return"] - amplitudes["closed"]) < 1e-12
    assert abs(amplitudes["projection_inside"]) < 1e-12
    assert all(abs(value.imag) < 1e-12 for value in amplitudes.values())
    rows.append({"tau": tau, **{name: value.real for name, value in amplitudes.items()}})

print(json.dumps({
    "source": "verification/v989_lattice_gate_battery.py",
    "beta": beta, "time_shift": float(energies[0]),
    "dimension": native.D_H, "projected_rank": basis.shape[1],
    "max_exact_entry_residuals": checks, "amplitudes": rows,
}, indent=2))
