"""Test the literal fold Cartan and phase on the existing R032 source operator."""

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sparse
import sympy as sp

repo = Path(sys.argv[1]).resolve()
contract = repo / "experiments/theory-contracts/chat-triangulation-20261108"
sys.path.insert(0, str(repo / "verification"))
spec = importlib.util.spec_from_file_location("quarter", contract / "check_quarter_reference.py")
quarter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(quarter)

gamma = sparse.csr_matrix(quarter.grading)
majorana = sparse.csr_matrix(-1j * quarter.grading @ quarter.mass_matrix)
yukawa = sparse.csr_matrix(-1j * quarter.grading @ quarter.higgs_matrix)
modes = np.arange(-8, 8)
one = sparse.eye(len(modes), format="csr")
neighbors = sparse.diags([np.ones(len(modes) - 1)] * 2, [-1, 1], format="csr")
mode_identity = sparse.eye(gamma.shape[0], format="csr")
selected = (modes % 8 == 0).astype(float)
projector = sparse.kron(sparse.diags(selected), mode_identity, format="csr")
phase = sparse.kron(sparse.diags((1j)**(modes % 4)), mode_identity, format="csr")
cartan_rotation = sparse.kron(sparse.diags((1j)**selected), mode_identity, format="csr")
u = 0.23 / 2  # unchanged profile in check_fold_r032_bulk_20261109.py
yukawa_blocks = sparse.kron(u * neighbors, yukawa, format="csr")
translated_neighbors = sparse.diags(
    [-1j * np.ones(len(modes) - 1), 1j * np.ones(len(modes) - 1)], [1, -1], format="csr"
)
translated_yukawa = sparse.kron(u * translated_neighbors, yukawa, format="csr")

roots = sp.I ** (sp.Rational(1, 2))  # eighth root with positive imaginary part
for m in modes:
    filtered = sum(roots ** (int(m) * j) for j in range(8)) / 8
    assert sp.simplify(sp.expand_complex(filtered)) == int(m % 8 == 0)


def square_norm(matrix):
    return float(np.vdot(matrix.data, matrix.data).real)


rows = []
for alpha in (0.25, 0.75):
    arcs = np.arange(8)
    folded_vector = np.exp(2j * np.pi * alpha * arcs / 8) / np.sqrt(8)
    cartan_arc = np.outer(folded_vector, folded_vector.conj())
    np.testing.assert_allclose(cartan_arc @ cartan_arc, cartan_arc, atol=1e-15)
    np.testing.assert_allclose(np.diag(cartan_arc), np.ones(8) / 8, atol=1e-15)
    off_diagonal = cartan_arc - np.diag(np.diag(cartan_arc))
    np.testing.assert_allclose(np.linalg.norm(off_diagonal)**2, 7 / 8, atol=1e-15)

    base = sparse.kron(sparse.diags(modes + alpha), gamma) + sparse.kron(one, majorana)
    assert (base @ projector - projector @ base).nnz == 0
    assert (base @ phase - phase @ base).nnz == 0
    h = base + yukawa_blocks
    marker_norm = square_norm(h @ projector - projector @ h)
    phase_norm = square_norm(h @ phase - phase @ h)
    scale = u**2 * square_norm(yukawa)
    np.testing.assert_allclose(marker_norm, 6 * scale)
    np.testing.assert_allclose(phase_norm, 60 * scale)
    assert (phase @ h @ phase.getH() - base - translated_yukawa).nnz == 0
    cartan_moved = cartan_rotation @ h @ cartan_rotation.getH()
    block_size = gamma.shape[0]
    edge_from_minus8 = cartan_moved[block_size:2 * block_size, :block_size]
    edge_from_minus7 = cartan_moved[2 * block_size:3 * block_size, block_size:2 * block_size]
    assert (edge_from_minus8 + 1j * u * yukawa).nnz == 0
    assert (edge_from_minus7 - u * yukawa).nnz == 0

    rows.append({
        "quarter": alpha,
        "arc_off_diagonal_squared": float(np.linalg.norm(off_diagonal)**2),
        "commutator_cartan_squared": marker_norm,
        "commutator_phase_squared": phase_norm,
        "relative_cartan_squared": 6,
        "relative_phase_squared": 60,
    })

print(json.dumps({
    "unchanged_higgs_profile_amplitude": 0.23,
    "yukawa_rank": int(np.linalg.matrix_rank(yukawa.toarray())),
    "finite_signed_modes": modes.tolist(),
    "phase_rotates_whole_geometric_background": True,
    "cartan_rotation_equal_fourier_displacement_factors": ["-i", "1"],
    "checks": rows,
}, indent=2))
