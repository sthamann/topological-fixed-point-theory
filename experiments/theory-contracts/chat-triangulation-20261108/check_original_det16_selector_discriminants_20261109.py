"""Distinguish native gauge, number parity and rank-one DET signs on actual DET16."""
import json
import sys
from pathlib import Path

import numpy as np

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo / "verification"))
import v1002_spin10_mirror_projector as native

n = native.SPIN10_MODES
dimension = 1 << n
modes = tuple(range(n))
projector = native.determinant_projector(n, modes, native.PHI)
states = np.arange(dimension)
fock_vacuum = np.zeros(dimension, dtype=complex)
fock_vacuum[0] = 1
phi = projector @ fock_vacuum
phi /= np.linalg.norm(phi)
assert np.max(np.abs(projector@phi-phi)) < 1e-14

# Read the original positive Weyl basis and its Cartan generators.
gamma = native.gamma_matrices_10d()
gamma_product = np.linalg.multi_dot(gamma)
positive = np.flatnonzero(np.diag((-1j)**5*gamma_product).real > .5)
generators,pairs = native.positive_weyl_generators()
cartans = [generators[pairs.index((2*j,2*j+1))] for j in range(5)]
cartan_weights = np.array([np.diag(x).real for x in cartans]).T
expected_weights = np.array([[(bit >> (4-j) & 1)-.5 for j in range(5)]
                             for bit in positive])
assert np.array_equal(cartan_weights,expected_weights)
assert all(bit.bit_count()%2 == 0 for bit in positive.tolist())

# E=E3+E2, nu^c is the internal vacuum and e^c the top weak wedge.
# This fixes two native carrier coordinates, not a map to the physical net.
nu_index = int(np.flatnonzero(positive==0)[0])
e_index = int(np.flatnonzero(positive==3)[0])
angles = np.pi*np.array([-1,-1,1,1,0])
carrier_sign = np.rint(np.exp(1j*cartan_weights@angles).real).astype(int)
assert set(carrier_sign) == {-1,1}
assert carrier_sign[nu_index] == 1 and carrier_sign[e_index] == -1
assert np.prod(carrier_sign) == 1
gauge_diag = np.ones(dimension)
for i,sign in enumerate(carrier_sign):
    if sign == -1:
        gauge_diag[(states>>i&1) == 1] *= -1
parity_diag = 1-2*native.parity_table(n)

# Same flat DET term and coefficient in the original v1013 table.
def propagate(v,time):
    retained = projector@v
    return retained+np.exp(-2*time)*(v-retained)


def identity(v):
    return v


def gauge(v):
    return gauge_diag*v


def parity(v):
    return parity_diag*v


def rank_one_det_sign(v):
    return 2*(projector@v)-v


def number_contrast(v):
    return (2*(states>>e_index&1)-1)*v


markers = {"gauge":gauge,"fermion_parity":parity,"rank_one_DET":rank_one_det_sign}
endpoint_states = {
    "empty":phi,
    "charged_e_conjugate_address":native.mode_operator(n,e_index,True)@phi,
    "neutral_nu_conjugate_address":native.mode_operator(n,nu_index,True)@phi,
    "neutral_number_contrast":number_contrast(phi),
}
half_time = .4
rows = {}
for name,v in endpoint_states.items():
    half = propagate(v,half_time)
    plain = np.vdot(half,half).real
    rows[name] = {"unmarked":float(plain)}
    for label,mark in markers.items():
        marked = np.vdot(half,mark(half))
        assert abs(marked.imag) < 1e-14
        ratio = marked.real/plain
        assert abs(abs(ratio)-1) < 1e-14
        rows[name][label] = {"marked":float(marked.real),"ratio":float(np.rint(ratio))}

assert [rows["charged_e_conjugate_address"][label]["ratio"] for label in markers] == [-1,-1,-1]
assert [rows["neutral_nu_conjugate_address"][label]["ratio"] for label in markers] == [1,-1,-1]
assert [rows["neutral_number_contrast"][label]["ratio"] for label in markers] == [1,1,-1]
contrast = endpoint_states["neutral_number_contrast"]
assert np.allclose(np.flatnonzero(contrast),[0,dimension-1])
assert abs(np.vdot(phi,contrast)) < 1e-14

# An independently exact test of the neutral witness from only the original
# two-plane algebra. E_0F is the rotation-free exchange whose + ray is phi.
alpha = np.exp(1j*native.PHI)
two_plane = np.array([[0,alpha.conjugate()],[alpha,0]])
assert np.max(np.abs(two_plane@contrast[[0,dimension-1]]+contrast[[0,dimension-1]])) < 1e-14

print(json.dumps({
    "scope":"actual native DET16 mirror factor and native Weyl carrier; no original QTheta, joint metric or physical net identification",
    "native_mode_indices":{"e_conjugate":e_index,"nu_conjugate":nu_index},
    "native_internal_wedge_masks":{"e_conjugate":int(positive[e_index]),"nu_conjugate":int(positive[nu_index])},
    "state":"actual local phi ray, not full physical source state",
    "original_flat_energy":2,
    "half_time":half_time,
    "results":rows,
    "neutral_number_contrast_stays_in_empty_full_plane":True,
},indent=2))
