"""Null-word check on the existing R032/96 two-face comparison."""

from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from scipy.special import expit


repo = Path(sys.argv[1]).resolve()
comparison = repo / "experiments/theory-contracts/chat-triangulation-20261108"
sys.path[:0] = [str(repo / "verification"), str(comparison)]
entry = comparison / "check_r113_dressed_open_20261109.py"
sys.argv = [str(entry), str(comparison)]
with redirect_stdout(io.StringIO()):
    original = runpy.run_path(str(entry))

from check_mixed_higgs_eich import gauge_operators

native, grading = original["native"], original["grading"]
mass = -1j * grading @ (
    original["majorana"] + original["inner_form"](np.array(original["half_higgs"][0]))
)
np.testing.assert_allclose(mass, mass.conj().T, atol=1e-14)
weak, hypercharge = gauge_operators()
incoming = native.IDX[("nu", "R", 0, 0)]
outgoing = native.IDX[("d", "R", 0, 0)]
operator = np.zeros_like(mass)
operator[outgoing, incoming] = 1
assert np.count_nonzero(6 * (hypercharge @ operator - operator @ hypercharge) + 2 * operator) == 0
assert all(np.count_nonzero(g @ operator - operator @ g) == 0 for g in weak)
assert np.count_nonzero(grading @ operator - operator @ grading) == 0


def sqrt_marginals(momentum):
    energies, vectors = np.linalg.eigh(momentum * grading + mass)
    c = (vectors * np.sqrt(expit(-4 * np.pi * energies))) @ vectors.conj().T
    one_minus_c = (vectors * np.sqrt(expit(4 * np.pi * energies))) @ vectors.conj().T
    return c, one_minus_c


for alpha in (0.25, 0.75):
    # Independent Fourier blocks. Each analytic block is rank one:
    # both square roots are injective at finite inverse temperature.
    squared_norms = []
    for k in range(-2, 2):
        p = k + alpha
        _, incoming_factor = sqrt_marginals(p)
        outgoing_factor, _ = sqrt_marginals(p + 1)
        block = outgoing_factor @ operator @ incoming_factor
        norm_squared = float(np.vdot(block, block).real)
        assert norm_squared > 0
        squared_norms.append(norm_squared)
    largest, second = sorted(squared_norms, reverse=True)[:2]
    lower_bound = 4 * largest * second
    assert lower_bound > 0
    print(json.dumps({
        "original_quarter": alpha,
        "probe": "native nu_R to d_R; matching Y and weak charge only",
        "particle_hole_block_norm_squared": squared_norms,
        "four_excitation_norm_squared_lower_bound": lower_bound,
    }))

# A sharp vacuum polarization of the actual free S-reference Hamiltonian.
# On the right-handed modes h_S(p)=-p, hence a +1 shift crosses exactly once.
for alpha in (0.25, 0.75):
    source_modes = [k + alpha for k in range(-3, 3)]
    crossings = [p for p in source_modes if p < 0 and p + 1 > 0]
    assert crossings == [alpha - 1]
assert np.count_nonzero(operator @ operator) == 0
assert np.array_equal(
    (operator @ operator.conj().T - operator.conj().T @ operator) @ operator
    - operator @ (operator @ operator.conj().T - operator.conj().T @ operator),
    2 * operator,
)
print("Sharp h_S vacuum comparison: one crossing, one-pair norm 1, neutral square 0.")
