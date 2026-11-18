"""Compare the old C58 Slater cap with the R113 exterior vacuum type."""

import contextlib
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from scipy.linalg import expm, solve


repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
sys.path.insert(0, str(repo / "verification"))
sys.argv = ["check_r113_dressed_open_20261109.py", str(base)]
with contextlib.redirect_stdout(io.StringIO()):
    source = runpy.run_path(str(base / sys.argv[0]))

gamma = source["grading"]
identity = np.eye(len(gamma))
ell = 1 / 8  # Actual old C58 width.
native, swap = source["native"], source["swap"]


def real_lie(a):
    return -1j * (a + swap @ a.conj() @ swap)


weak = real_lie(native.rep(0, 1j * np.array([[0., 1.], [1., 0.]]), np.zeros((3, 3))))
hyper = real_lie(native.rep(3j, np.zeros((2, 2)), -1j * np.eye(3))) / 6
weak_kinetic, hyper_kinetic = weak @ gamma, hyper @ gamma
assert np.array_equal(weak_kinetic @ hyper_kinetic, hyper_kinetic @ weak_kinetic)
rows = {}
for amplitude in (0., .55):
    phi = source["inner_form"](np.array([0., amplitude])) + source["majorana"]
    for momentum in (.25, -.25):
        h0 = momentum * gamma
        h = h0 - 1j * gamma @ phi
        assert np.max(abs(h-h.conj().T)) < 1e-14
        transfer = expm(-ell * h)
        graph00 = solve(identity + transfer.conj().T @ transfer, identity, assume_a="pos")
        graph10 = transfer @ graph00
        recovered = solve(graph00.T, graph10.T).T
        graph_error = np.max(abs(recovered-transfer))
        assert graph_error < 1e-13
        old = np.flatnonzero(np.diag(h0) < 0)
        new = np.flatnonzero(np.diag(h0) > 0)
        assert len(old) == len(new) == 48
        occupied = recovered[:, old]
        overlap = np.linalg.slogdet(occupied[old, :])[1]
        norm_square = np.linalg.slogdet(occupied.conj().T @ occupied)[1]
        # These are determinants in the original exterior norm of the Slater vector.
        old_ray_fraction = np.exp(2*overlap - norm_square)
        log_free_relative = overlap + ell * np.trace(h0[np.ix_(old, old)]).real
        log_excitation = overlap + ell * np.linalg.eigvalsh(h)[:len(old)].sum()
        off_diagonal = occupied[new, :]
        assert np.linalg.norm(off_diagonal, 2) > 1e-6 and old_ray_fraction < 1
        # In the true hole chart the Gauss factors retain pair creation and
        # annihilation, in addition to the old occupied determinant scalar.
        block_a = occupied[old, :]
        block_b = recovered[np.ix_(old, new)]
        block_d = recovered[np.ix_(new, new)]
        pair_graph = solve(block_a.T, off_diagonal.T).T
        schur = block_d - pair_graph @ block_b
        assert min(np.linalg.eigvalsh(schur)) > 0
        zero = np.zeros_like(block_a)
        small_identity = np.eye(len(old))
        reconstructed = (np.block([[small_identity, zero], [pair_graph, small_identity]])
                         @ np.block([[block_a, zero], [zero, schur]])
                         @ np.block([[small_identity, solve(block_a, block_b)], [zero, small_identity]]))
        expected = np.block([[block_a, block_b], [off_diagonal, block_d]])
        factor_error = np.max(abs(reconstructed-expected))
        assert factor_error < 1e-13
        pair_log_norm = np.linalg.slogdet(small_identity + pair_graph.conj().T @ pair_graph)[1]
        assert abs(np.exp(-pair_log_norm) - old_ray_fraction) < 1e-13
        null_mixed = None
        if amplitude == 0:
            for insertion in (weak_kinetic, hyper_kinetic):
                assert np.array_equal(insertion @ h, h @ insertion)
                assert np.array_equal(insertion[np.ix_(old, new)], np.zeros_like(block_a))
                assert abs(np.trace(insertion[np.ix_(old, old)])) < 1e-14
            # Commutativity and block preservation exponentiate exactly.
            # These finite values independently replay the resulting identity.
            values = []
            for a, b in ((0., 0.), (.02, 0.), (0., .03), (.02, .03)):
                perturbed = expm(-ell * (h+a*weak_kinetic+b*hyper_kinetic))
                values.append(np.linalg.slogdet(perturbed[np.ix_(old, old)])[1])
            assert max(abs(value-overlap) for value in values) < 1e-13
            null_mixed = {"commuting_constant_gauge_directions": True,
                          "occupied_weak_and_hypercharge_traces": [0, 0],
                          "deterministic_mixed_log_cap_jet": 0,
                          "largest_finite_log_cap_change": float(max(abs(value-overlap) for value in values))}
        rows[f"H={amplitude},p={momentum}"] = {
            "source_graph_transfer_error": float(graph_error),
            "particle_hole_gauss_factorization_error": float(factor_error),
            "particle_hole_pair_graph_norm": float(np.linalg.norm(pair_graph, 2)),
            "off_diagonal_particle_hole_transfer_norm": float(np.linalg.norm(off_diagonal, 2)),
            "off_diagonal_particle_hole_transfer_rank": int(np.linalg.matrix_rank(off_diagonal, tol=1e-10)),
            "slater_retained_norm_fraction": float(old_ray_fraction),
            "slater_excitation_log_cap": float(log_excitation),
            "slater_free_relative_log_cap": float(log_free_relative),
            "unweighted_exterior_grade_zero_matrix_element": 1,
            "Higgs_zero_fixed_reference_gauge_test": null_mixed,
        }

print(json.dumps({"rows": rows,
                  "scope": "Native fixed-H R032 and old Slater C58 fiber. Neither actual boson integration nor QTheta/gap specified."}, indent=2))
