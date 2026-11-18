"""Original Pin Green reflection on the massive native pair and its NS cover."""

import contextlib
import io
from itertools import product
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss


repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
with contextlib.redirect_stdout(io.StringIO()):
    native = runpy.run_path(str(base / "check_tfpt_conditioned_r113_weights_20261110.py"))
phi_m = native["majorana"]
phi_h = native["inner"](np.array([0., 1.]))
assert np.linalg.matrix_rank(phi_h) == 96
z = 1j
seed = (z, 1)


def higgs_move(label):
    c, b = label
    return -c, -b


def reflected_move(label):
    c, b = label
    return -c, b


higgs_pair = {seed, higgs_move(seed)}
assert {c*b for c, b in higgs_pair} == {1j}
orbit = higgs_pair | {reflected_move(label) for label in higgs_pair}
assert len(orbit) == 4
assert {higgs_move(label) for label in orbit} == orbit
assert {reflected_move(label) for label in orbit} == orbit
labels = sorted(orbit,key=lambda label:(label[0].imag,label[1]))
scalar_signs = []
for signs in product((1,-1),repeat=len(labels)):
    q = dict(zip(labels,signs,strict=True))
    if all(q[label]==q[higgs_move(label)]==q[reflected_move(label)] for label in labels):
        scalar_signs.append(signs)
assert scalar_signs==[(1,)*4,(-1,)*4]


def spatial(epsilon, cutoff):
    # On b=+1 the original spatial H mode n has energy n*pi/2 and
    # J=sigma1 R_w character (-1)^n. The n=0 input lies in the even sector.
    # The transported bounded scalar rho and sigma1 preserve this sector.
    n = 2*np.arange(-cutoff,cutoff+1)
    energy = n*np.pi/2
    nodes, weights = leggauss(max(256,8*cutoff))
    w = (nodes+1)/2
    bump_coord = (w-.5)/.25
    bump = np.zeros_like(w)
    inside = abs(bump_coord)<1
    bump[inside] = np.exp(1-1/(1-bump_coord[inside]**2))
    rho_minus_one = np.expm1(epsilon*bump)
    frequencies = np.arange(-2*cutoff,2*cutoff+1)
    coeff = np.cos(np.pi*frequencies[:,None]*w) @ (weights*rho_minus_one/2)
    rho = np.eye(len(n))+coeff[(n[:,None]-n[None,:])//2+2*cutoff]
    # sigma1 is n -> -n on this b=+1 even sector.
    mark = rho[::-1,:]
    np.testing.assert_allclose(mark,mark.T,atol=1e-14)
    return energy,mark,mark[:,cutoff]


def multipliers(energy, u):
    denom = np.logaddexp(0.,2*energy)
    return np.exp(u*energy-denom),np.exp((u+1)*energy-denom)


def compress(epsilon, cutoff, radius, t, s):
    energy,mark,test = spatial(epsilon,cutoff)
    masses, eigenbasis = np.linalg.eigh(phi_m+radius*phi_h)
    phi_h_in_mass_basis = eigenbasis.conj().T@phi_h
    # tau2 diagonalization on the paired b=+/- domains yields
    # the two spatial blocks H +/- sigma1*rho*m for each native mass.
    normal_ns = np.zeros(len(masses))
    normal_im = np.zeros(len(masses))
    for mass in np.unique(masses.round(11)):
        mask = abs(masses-mass)<1e-10
        exact_mass = float(masses[mask].mean())
        for sign in (1,-1):
            eig,u = np.linalg.eigh(np.diag(energy)+sign*exact_mass*mark)
            overlaps = abs(u.conj().T@test)**2
            real,imag = multipliers(eig,t+s)
            normal_ns[mask] += (overlaps@real)/2
            normal_im[mask] += (overlaps@imag)/2
    ns = phi_h_in_mass_basis.conj().T @ (normal_ns[:,None]*phi_h_in_mass_basis)
    im = phi_h_in_mass_basis.conj().T @ (normal_im[:,None]*phi_h_in_mass_basis)
    assert np.linalg.eigvalsh(ns).min()>0 and np.linalg.eigvalsh(im).min()>0
    return ns,im


def direct_paired_spin_control():
    energy,mark,test = spatial(.1,8)
    mass = np.linalg.eigvalsh(phi_m+.55*phi_h)[-1]
    h = np.diag(energy)
    full = np.block([[h,1j*mass*mark],[-1j*mass*mark,h]])
    direct_energy, direct_vectors = np.linalg.eigh(full)
    vector = np.r_[np.zeros(len(energy)),test]
    direct = (abs(direct_vectors.conj().T@vector)**2
              @ np.stack(multipliers(direct_energy,.3),axis=1))
    reduced = np.zeros(2)
    for sign in (1,-1):
        eig,u = np.linalg.eigh(h+sign*mass*mark)
        reduced += abs(u.conj().T@test)**2@np.stack(multipliers(eig,.3),axis=1)/2
    np.testing.assert_allclose(direct,reduced,atol=1e-13)
    return float(np.max(abs(direct-reduced)))


# Check the actual twisted Green boundary and jump on sample masses.
for e in (-2.,0.,.7,3.):
    for q in (1j,-1j):
        cp = 1/(1-np.exp(-e)/q)
        assert abs(np.exp(-e)*cp-q*(cp-1))<1e-14
        real,imag = multipliers(np.array([e]),.5)
        direct = -np.exp(.5*e)*(cp-1)
        np.testing.assert_allclose(direct,real[0]+q*imag[0],atol=1e-14)

out = {}
for epsilon in (0.,.1):
    for radius in (0.,.55):
        values = {}
        for cutoff in (16,32,64):
            real,imag = compress(epsilon,cutoff,radius,.2,.3)
            values[str(cutoff)] = {
                "cover_return_trace_over_96":float(np.trace(real).real/96),
                "fixed_quarter_imaginary_trace_over_96":float(np.trace(imag).real/96),
                "cover_minimum_gram_eigenvalue":float(np.linalg.eigvalsh(real)[0]),
                "quarter_imaginary_minimum_gram_eigenvalue":float(np.linalg.eigvalsh(imag)[0]),
                "cover_rank":int(np.linalg.matrix_rank(real)),
            }
        np.testing.assert_allclose(list(values["32"].values()),list(values["64"].values()),
                                   atol=2e-9)
        out[f"epsilon={epsilon},H={radius}"] = values

# For a flat geometric Dirac the canonical paired square is D0²+Phi².
# Keep this gap distinct from the spatial generator and relative S/T terms.
gaps = {}
for radius in (0.,.55):
    phi = phi_m+radius*phi_h
    masses = np.linalg.eigvalsh(phi)
    relative_s = float(min(abs(masses)))
    zero = np.zeros_like(phi)
    d_rel_s = np.block([[zero,phi],[phi,zero]])
    first = d_rel_s[:96,:96]
    square = (d_rel_s.conj().T@d_rel_s)[:96,:96]
    np.testing.assert_allclose(first,zero,atol=1e-15)
    gaps[str(radius)] = {"absolute_flat_paired_dirac_gap":float(np.hypot(np.pi/2,relative_s)),
                         "flat_spatial_Higgs_generator_gap":relative_s,
                         "S_relative_block_minimum_singular_value":relative_s,
                         "T_relative_block_minimum_singular_value":0.,
                         "S_first_order_bare_seam_compression_rank":int(np.linalg.matrix_rank(first)),
                         "S_ambient_square_compressed_to_bare_seam_rank":int(np.linalg.matrix_rank(square))}
phase = {}
for t,s in ((.1,.2),(.2,.3),(.3,.4)):
    re,im = compress(.1,64,.55,t,s)
    phase[f"{t+s:.1f}"] = float(np.trace(im).real/np.trace(re).real)
assert max(phase.values())-min(phase.values())>.1
early,late = .15,.35
early_gram,_ = compress(.1,64,.55,early,early)
late_gram,_ = compress(.1,64,.55,late,late)
cross_gram,_ = compress(.1,64,.55,early,late)
two_times = np.block([[early_gram,cross_gram],[cross_gram.conj().T,late_gram]])
assert np.linalg.eigvalsh(two_times)[0]>-1e-12
print(json.dumps({
    "Higgs_pair_labels":sorted([[str(c),b] for c,b in higgs_pair]),
    "closure_under_original_Green_reflection_and_Higgs":sorted([[str(c),b] for c,b in orbit]),
    "same_native_Higgs_ends_both_quarters_and_NS":out,
    "fixed_quarter_trace_imaginary_over_real_same_state_by_time_sum":phase,
    "cover_two_time_gram":{
        "times":[early,late],"rank":int(np.linalg.matrix_rank(two_times)),
        "minimum_eigenvalue":float(np.linalg.eigvalsh(two_times)[0])
    },
    "independent_full_pair_spin_control_error":direct_paired_spin_control(),
    "label_scalar_signs_commuting_with_higgs_and_reflection":scalar_signs,
    "flat_distinct_gap_objects":gaps,
    "scope":"Static native paired Pin Green comparison in original x time and period-two NS cover. No QTheta or D_adm assignment, no R113 weight."
},indent=2))
