"""Native finite SU2/U2 link scans, covariant Gaussian and charged Gram tests."""

from fractions import Fraction as F
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import expm_frechet
from scipy.special import eval_chebyu


repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
old = runpy.run_path(str(base / "check_tfpt_conditioned_r113_weights_20261110.py"))
alpha = 1 / 137.0359992168407  # v3 result; independent implementation recomputes the root.
t2 = 16 * np.pi * alpha / 3
t1 = .6 * t2
radius = old["radius_max"]
width, duration = old["width"], old["duration"]
blocks, identity = old["blocks"], old["identity"]
gamma = old["gamma"]
native = old["native"]
he, det_y, ref_trace = (old[name] for name in ("hermitian_exponential", "determinant_and_yukawa", "ref_trace"))
majorana, higgs = old["clifford_majorana"], old["clifford_higgs"]
momentum = .25
h0 = momentum * gamma
late_h = h0 + majorana + radius * higgs
middle, middle_jet = expm_frechet(-2 * width * late_h, -2 * width * radius * higgs)
late_sea = np.linalg.eigvalsh(late_h)[:48].sum()
log_ref = ref_trace(h0)


def native_values(r, link_blocks):
    early, energies = he(h0 + majorana + r * higgs, width)
    sea = 2 * width * (energies[:48].sum() + late_sea)
    return det_y(early, middle, middle_jet, link_blocks, sea - log_ref)


def records(masses, f0, labels):
    z, zu, zy, zuy = masses.T
    c = zuy / z - (zu / z) * (zy / z)
    return {label: dict(Z0=float(z[i]), F0=float(f0[i]), C_UY=float(c[i]),
                        F0_over_Z0=float(f0[i] / z[i])) for i, label in enumerate(labels)}


def character_link_scan():
    links, _ = old["weak_arrays"](49, 25)
    theta = np.arange(1, 50) * np.pi / 50
    _, axis_weight = leggauss(25)
    order = np.arange(49)
    j = order / 2
    casimir = j * (j + 1)
    a_rate = casimir / 1000 + (j > 0)
    labels = ["Haar", "heat_0.10", "heat_alpha", "heat_0.15", "heat_0.20", "heat_0.25", "heat_2pi", "A_0.8", "A_2pi"]
    fourier = np.array([j == 0, *[np.exp(-t * casimir) for t in (.10, t2, .15, .20, .25, 2*np.pi)],
                       np.exp(-.8 * a_rate), np.exp(-2*np.pi * a_rate)])
    weights = (2 / 50 * np.sin(theta) * ((fourier * (order + 1)) @ np.sin(np.outer(order + 1, theta))))
    weights = weights[:, :, None] * axis_weight[None, None, :] / 2
    assert np.max(abs(weights.sum(axis=(1, 2)) - 1)) < 1e-12
    points, coeffs = leggauss(20)
    result = np.zeros((len(labels), 4))
    for r, mass in zip(radius * (points + 1) / 2, coeffs * (points + 1)**3 / 4):
        z, y = native_values(r, links)
        az, ay = (np.einsum("kij,ij->k", weights, v) for v in (z, y))
        result += mass * np.array([az, r*r*az, ay, r*r*ay]).T
    f0 = null_fermion() * 2 / (np.pi**2 * radius**4)
    return records(result, np.full(len(labels), f0), labels)


def null_fermion():
    k0, _ = native_values(0, [np.eye(len(i)) for i, _ in blocks])
    return float(k0)


def disk_heat_grid():
    # For G=e^(i phi/2) U and h_L,h_R along e2, the relative SU2
    # matrix is X=U exp(-i phi T3). Thus X22=U22 exp(i phi/2).
    # Haar U22 is uniform on the complex unit disk. The circle heat
    # multiplies Fourier angular mode n by exp[-t1 n²/4].
    points, masses = leggauss(25)
    psi = 2 * np.pi * np.arange(100) / 100
    rr = np.sqrt((1 + points) / 2)[:, None]
    a, z = rr * np.cos(psi), -rr * np.sin(psi)
    x = np.sqrt(1 - rr**2) * np.ones_like(psi)
    links = []
    for i, _ in blocks:
        fixed, doublet, ix, iz = (old[name][np.ix_(i, i)] for name in ("fixed", "doublet", "ix", "iz"))
        links.append(fixed + a[..., None, None] * doublet + x[..., None, None] * ix + z[..., None, None] * iz)
    density = sum((n + 1) * np.exp(-t2 * n * (n + 2) / 4) * eval_chebyu(n, a) for n in range(49))
    harmonics = np.fft.fftfreq(100, 1 / 100)
    density_spectrum = np.fft.fft(density, axis=1)
    return links, masses, density, harmonics, density_spectrum


def disk_joint_heat():
    links, masses, density, harmonics, density_spectrum = disk_heat_grid()
    joint = np.fft.ifft(density_spectrum * np.exp(-t1 * harmonics**2 / 4), axis=1).real
    weights = masses[:, None] * np.array([density, joint]) / 200
    assert np.max(abs(weights.sum(axis=(1, 2)) - 1)) < 1e-12
    r_nodes, r_mass = leggauss(20)
    result = np.zeros((2, 4))
    for r, mass in zip(radius * (1 + r_nodes) / 2, r_mass * (1 + r_nodes)**3 / 4):
        z0, y = native_values(r, links)
        az, ay = (np.einsum("kij,ij->k", weights, v) for v in (z0, y))
        result += mass * np.array([az, r*r*az, ay, r*r*ay]).T
    f0 = null_fermion() * 2 / (np.pi**2 * radius**4)
    # Twisting the X22 character by the open e_R determinant changes
    # the abelian Fourier exponent, so treating X and phi as independent
    # is not the central U2 heat law. This ratio is an exact class test.
    e = np.exp(-t1)
    moment_base = np.exp(-.75 * t2 - t1 / 4)
    moment_charged = np.exp(-.75 * t2 - t1 * (1 - .5)**2)
    moment_charged_conjugate_orientation = np.exp(-.75 * t2 - t1 * (1 + .5)**2)
    return {"rows": records(result, [f0, f0], ["SU2_disk_control", "central_U2_heat"]),
            "determinant_character_e_vs_nu": float(e),
            "X22_mean": float(moment_base),
            "X22_times_e_character_mean": float(moment_charged),
            "ratio_for_this_charged_orientation_moment": float(moment_charged / moment_base),
            "conjugate_X22_times_e_character_mean": float(moment_charged_conjugate_orientation),
            "ratio_for_conjugate_charged_orientation_moment": float(moment_charged_conjugate_orientation / moment_base),
            "factorized_product_for_comparison": float(e * moment_base)}


def gaussian_source():
    # Change variable to transported H_L. The weak link disappears from
    # the scalar half-word with just the two fixed end frames. For the
    # Gaussian law only the transported Higgs value remains.
    taus = np.array([.02, .05, .10])
    labels = [f"tau_{tau:.2f}" for tau in taus]
    theta = np.arange(1, 50) * np.pi / 50
    _, axis_weight = leggauss(25)
    haar = 2 / 50 * np.sin(theta)[:, None]**2 * axis_weight[None, :] / 2
    links, _ = old["weak_arrays"](49, 25)
    end = radius + 10 * np.sqrt(max(taus))
    points, masses = leggauss(38)
    result = np.zeros((len(taus), 4))
    measure_mass = np.zeros(len(taus))
    for r, dr in zip(end * (points + 1) / 2, masses * end / 2):
        z, y = native_values(r, links)
        exponent = -(r*r + radius**2 - 2*r*radius*np.cos(theta))
        density = (r**3 / (2 * taus**2))[:, None, None] * np.exp(exponent[None, :, None] / (2 * taus[:, None, None])) * haar
        measure_mass += dr * density.sum(axis=(1, 2))
        az, ay = (np.einsum("kij,ij->k", density, v) for v in (z, y))
        result += dr * np.array([az, r*r*az, ay, r*r*ay]).T
    assert np.max(abs(measure_mass - 1)) < 1e-8
    f0 = null_fermion() * np.exp(-radius**2 / (2 * taus)) / (2 * np.pi * taus)**2
    return {"rows": records(result, f0, labels), "integrated_gaussian_mass": measure_mass.tolist(),
            "upper_radial_bound_for_quadrature": float(end)}


def exact_and_reflection_controls():
    # An ordinary positive density can have a negative transfer channel.
    # chi_1 has Haar norm one; d_1=3. Min density >= 1/2.
    bad_fourier = -F(1, 18)
    # On (SU2 x circle)/Z2, χ_1 cos(2θ) descends. Both jump marginals
    # equal Haar; the mixed Fourier coefficient is ±1/36.
    mixed_plus = F(1, 6) / (3 * 2)
    assert mixed_plus == F(1, 36) and bad_fourier < 0
    j = np.array([.5, 1., 1.5, 2., 25.])
    casimir = j * (j+1)
    apparent = .8 * (casimir / 1000 + 1 - (j == 25) / 5202) / casimir

    # Actual native E K and W K N columns in the lepton block. All
    # weak, Higgs, Majorana and frame operators are from the same chart.
    lepton, _ = blocks[0]
    sub = np.ix_(lepton, lepton)
    e_idx = native.IDX[("e", "R", 0, 0)]
    n_idx = native.IDX[("nu", "R", 0, 0)]
    e = identity[:, e_idx][lepton]
    n = identity[:, n_idx][lepton]
    early, _ = he(h0 + majorana + .35 * higgs, width)
    late, _ = he(late_h, width)
    p_e, p_n = (np.diag([1 if native.P[i % native.NP][:2] == (tag, "R") else 0 for i in range(96)])
                for tag in ("e", "nu"))
    op_e, op_n = (radius * higgs @ p_e)[sub], (.35 * higgs @ p_n)[sub]
    q = native.rep(0, 1j * old["sigma_x"], np.zeros((3,3)))
    swap = old["swap"]
    weak_mark = (-1j * (q + swap @ q.conj() @ swap) @ gamma)[sub]
    links, _ = old["weak_arrays"](49,25)
    rel_link = links[0]
    def reflected_entries(link):
        k = late[sub] @ link @ early[sub]
        vectors = np.stack((k @ e, k @ n, op_e @ k @ e, weak_mark @ k @ op_n @ n), axis=-1)
        return vectors.swapaxes(-1, -2).conj() @ vectors

    local = reflected_entries(rel_link)
    theta = np.arange(1,50) * np.pi/50
    _, w = leggauss(25)
    nn = np.arange(1,50)
    jj = (nn-1)/2
    a = np.array([np.exp(-t2*jj*(jj+1)), np.exp(-.8*(jj*(jj+1)/1000 + (jj>0)))])
    theta_mass = 2/50*np.sin(theta)*((a*nn) @ np.sin(np.outer(nn,theta)))
    weights = theta_mass[:, :, None] * w[None,None,:] / 2
    gram = np.einsum("kij,ijab->kab", weights, local)
    eig = np.linalg.eigvalsh(gram)
    assert eig.min() > 1e-5
    # Product split extension: average electron phase multiplies E/N
    # off diagonal by exp(-t1), keeping every diagonal identical.
    mask = np.array([-1,0,-1,0])
    phase = np.exp(-t1*(mask[:,None]-mask[None,:])**2)
    gram_heat = gram * phase
    assert np.linalg.eigvalsh(gram_heat).min() > 1e-5
    assert np.array_equal(gram_heat.diagonal(axis1=-2, axis2=-1), gram.diagonal(axis1=-2, axis2=-1))
    disk_links, disk_masses, _, harmonics, density_spectrum = disk_heat_grid()
    entries = reflected_entries(disk_links[0])
    shifts = mask[None, :] - mask[:, None]
    joint_gram = np.zeros((4, 4), complex)
    for charge in (-1, 0, 1):
        joint_density = np.fft.ifft(density_spectrum * np.exp(-t1 * (charge - harmonics/2)**2), axis=1)
        contracted = np.einsum("i,ij,ijab->ab", disk_masses / 200, joint_density, entries)
        joint_gram[shifts == charge] = contracted[shifts == charge]
    assert np.max(abs(joint_gram - joint_gram.conj().T)) < 1e-12
    assert np.linalg.eigvalsh(joint_gram).min() > 1e-5

    # Native U2 stabilizer det phase and the Wilson return. The invariant
    # two-leg comparison needs a fixed reference return in the same frame.
    phi = .37
    g = np.diag([np.exp(1j * phi), 1.])
    lam = np.exp(.5j * phi)
    color = np.exp(-1j * phi / 6) * np.eye(3)
    p = native.rep(lam, g / lam, color)
    lift = p @ swap @ p.conj() @ swap
    assert abs(lift[e_idx,e_idx] - np.exp(-1j * phi)) < 1e-13
    assert abs(lift[n_idx,n_idx] - 1) < 1e-13
    assert np.max(abs(lift @ early - early @ lift)) < 1e-13
    return dict(A_effective_heat_times=apparent.tolist(), A_probe_spins=j.tolist(),
                strictly_positive_pointwise_weight_negative_channel=str(bad_fourier),
                positive_joint_laws_identical_marginals_mixed_jump_coefficients=[str(mixed_plus), str(-mixed_plus)],
                finite_e_n_gram_eigenvalues_identity_phase=eig.tolist(),
                finite_e_n_gram_eigenvalues_heat_phase=np.linalg.eigvalsh(gram_heat).tolist(),
                finite_e_n_gram_eigenvalues_central_U2=np.linalg.eigvalsh(joint_gram).tolist(),
                finite_e_n_gram_central_U2_real=joint_gram.real.tolist(),
                finite_e_n_gram_central_U2_imaginary=joint_gram.imag.tolist(),
                phase_transfer_e_over_n=float(np.exp(-t1)),
                scopes="4x4 native tagged half-word span at fixed r=.35 and Z3-trivial lepton block; no full OS/P_adm. Exact U2 mixed jump example has Haar marginals and positive selfadjoint exponential." )


if __name__ == "__main__":
    print(json.dumps(dict(t2=float(t2), t1_physical_Y=float(t1), alpha=alpha)), flush=True)
    print(json.dumps({"character_scan": character_link_scan()}), flush=True)
    print(json.dumps({"central_U2": disk_joint_heat()}), flush=True)
    print(json.dumps({"Gaussian": gaussian_source()}), flush=True)
    print(json.dumps({"reflection": exact_and_reflection_controls()}), flush=True)
