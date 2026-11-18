"""Gauge-sign calibration in the original v989 Hamiltonian and its true Gauss sector."""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo / "verification"))
import v989_lattice_gate_battery as native

psi, num = native._build_ham_ops()
hamiltonian = native.build_H(psi, num)
g1, g2 = [native.gauss_unitary(num, site) for site in (1, 2)]
identity = np.eye(len(hamiltonian))
sign = g1 @ g1


def error(operator):
    return float(np.max(np.abs(operator)))


def neutral_projector(star):
    return sum(np.linalg.matrix_power(star, n) for n in range(4)) / 4


gauss = neutral_projector(g1) @ neutral_projector(g2)
assert error(gauss @ gauss-gauss) < 1e-13
assert error(hamiltonian @ gauss-gauss @ hamiltonian) < 1e-13
assert error(hamiltonian @ sign-sign @ hamiltonian) < 1e-13
assert error(sign @ gauss-gauss) < 1e-13

# Read the two real adjacent lattice ends from the native hopping word:
# psi_1^* Z_1 psi_2. The link dresses the left insertion to site 1's charge.
right = psi[1].conj().T
left = native._link_op(1, native.Z1_H) @ psi[2]
for field in (right, left):
    assert error(sign @ field @ sign + field) < 1e-13
    assert error(g2 @ field @ g2.conj().T-field) < 1e-13
assert error(g1 @ right @ g1.conj().T - 1j*right) < 1e-13
assert error(g1 @ left @ g1.conj().T + 1j*left) < 1e-13

# Same state as the existing Gauss receipt: the native projected Gibbs
# state at beta=1. Choose the strip times from the R032 question as a
# readout, not a fit or a change of the physical Hamiltonian.
beta, half_time = 1., .4
_, vectors = eigh(gauss, subset_by_value=(.5, 1.5))
assert vectors.shape[1] == 64
energies, rotation = eigh(vectors.conj().T @ hamiltonian @ vectors)
vectors = vectors @ rotation
full_energies, full_basis = eigh(hamiltonian)
def propagate(v, time):
    return full_basis @ (np.exp(-time*full_energies)[:,None]*(full_basis.conj().T @ v))


time = 2*half_time
density = np.exp(-(beta-time)*energies)
normalization = np.exp(-beta*energies).sum()
def value(field_l, field_r, marker):
    middle = propagate(field_r @ vectors, half_time)
    if marker is not None:
        middle = marker @ middle
    finish = field_l @ propagate(middle, half_time)
    return (density * np.sum(vectors.conj()*finish, axis=0)).sum()


empty = value(identity, identity, None)
empty_signed = value(identity, identity, sign)
charged = value(left, right, None)
charged_signed = value(left, right, sign)
assert abs(empty-normalization) < 1e-11*normalization
assert abs(empty_signed-empty) < 1e-11*normalization
assert abs(charged_signed+charged) < 1e-11*normalization
assert abs(charged) > .01*normalization

def pair(z):
    return [float(z.real), float(z.imag)]


print(json.dumps({
    "scope": "native v989 Gauss calibration, no TFPT QTheta or physical metric identification",
    "same_beta_as_existing_Gauss_receipt": beta,
    "half_times": [half_time,half_time],
    "actual_gauss_rank": vectors.shape[1],
    "raw_empty_partition_factor": float(normalization),
    "normalized_empty": pair(empty/normalization),
    "normalized_empty_with_G1_squared": pair(empty_signed/normalization),
    "normalized_two_field_ends": pair(charged/normalization),
    "normalized_two_field_ends_with_G1_squared": pair(charged_signed/normalization),
    "charged_odd_identity_residual": float(abs(charged+charged_signed)/normalization),
}, indent=2))
