"""Three declared positive weights on one native R032 window comparison."""

import contextlib
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import expm_frechet


repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
sys.path.insert(0, str(repo / "verification"))
with contextlib.redirect_stdout(io.StringIO()):
    args = sys.argv
    try:
        sys.argv = ["check_r113_dressed_open_20261109.py", str(base)]
        source = runpy.run_path(str(base / sys.argv[0]))
    finally:
        sys.argv = args

native = source["native"]
gamma, majorana, swap = (source[k] for k in ("grading", "majorana", "swap"))
inner, weak_lift = (source[k] for k in ("inner_form", "real_action"))
identity = np.eye(len(gamma))
late = .55  # Same amplitude as the established separated-window comparison.
radius_max = late  # Explicit conditional compact radial base, no fitted parameter.
width = .2
duration = 4 * width
omega = np.exp(2j * np.pi / 3)
labels = native.P
leptons = [i for i in range(96) if labels[i % native.NP][0] in ("e", "nu")]
quark_groups = [[i + sign * native.NP for i, label in enumerate(labels)
                 if label[0] in ("u", "d") and label[2] == color]
                for sign in range(2) for color in range(3)]
blocks = [(np.array(leptons), 0)] + [(np.array(ix), 1 if i < 3 else -1)
                                        for i, ix in enumerate(quark_groups)]
assert sorted(i for ix, _ in blocks for i in ix) == list(range(96))
assert sorted(len(ix) for ix, _ in blocks) == [12] * 6 + [24]
sigma_x = np.array([[0., 1.], [1., 0.]])
sigma_z = np.diag([1., -1.])
doublet = (weak_lift(np.eye(2)) - weak_lift(-np.eye(2))) / 2
assert np.count_nonzero(np.diag(doublet)) == 48
fixed = identity - doublet
ix = (weak_lift(1j * sigma_x) - weak_lift(-1j * sigma_x)) / 2
iz = (weak_lift(1j * sigma_z) - weak_lift(-1j * sigma_z)) / 2
clifford_higgs = -1j * gamma @ inner(np.array([0., 1.]))
clifford_majorana = -1j * gamma @ majorana


def weak_arrays(order, axis_order):
    # The 24 native doublets give degree <=48 in the reflected word.
    # Residual diagonal conjugation removes the azimuth around sigma_z.
    theta = np.arange(1, order + 1) * np.pi / (order + 1)
    axis, weights = leggauss(axis_order)
    a = np.cos(theta)[:, None]
    x = np.sin(theta)[:, None] * np.sqrt(1 - axis**2)[None, :]
    z = np.sin(theta)[:, None] * axis[None, :]
    link_blocks = [(fixed[np.ix_(i, i)] + a[..., None, None] * doublet[np.ix_(i, i)] +
                    x[..., None, None] * ix[np.ix_(i, i)] +
                    z[..., None, None] * iz[np.ix_(i, i)]) for i, _ in blocks]
    j = np.arange(order) / 2
    casimir = j * (j + 1)
    a_rate = casimir / 1000 + (j > 0) - (j == 25) / 5202
    central = np.array([j == 0, np.exp(-2 * np.pi * casimir), np.exp(-duration * a_rate)])
    theta_weight = 2 / (order + 1) * np.sin(theta) * (
        (central * (2 * j + 1)) @ np.sin(np.outer(2 * j + 1, theta)))
    assert np.max(abs(theta_weight.sum(axis=1) - 1)) < 1e-12
    return link_blocks, theta_weight[:, :, None] * weights[None, None, :] / 2


def determinant_and_yukawa(early, middle, middle_jet, links, log_sea_minus_ref):
    shape = links[0].shape[:-2]
    det_product = np.full((3, *shape), np.exp(log_sea_minus_ref), complex)
    log_derivative = np.zeros_like(det_product)
    for (i, charge), g in zip(blocks, links, strict=True):
        start = early[np.ix_(i, i)]
        k = g @ start
        adjoint = k.swapaxes(-1, -2).conj()
        r = adjoint @ middle[np.ix_(i, i)] @ k
        dr = adjoint @ middle_jet[np.ix_(i, i)] @ k
        for power in range(1 if charge == 0 else 3):
            phase = omega ** (charge * power)
            matrix = np.eye(len(i)) + phase * r
            sign, logdet = np.linalg.slogdet(matrix)
            derivative = np.linalg.solve(matrix, phase * dr).trace(axis1=-2, axis2=-1)
            if charge == 0:
                det_product *= sign * np.exp(logdet)
                log_derivative += derivative
            else:
                det_product[power] *= sign * np.exp(logdet)
                log_derivative[power] += derivative
    empty = det_product.mean(axis=0)
    marked = -(det_product * log_derivative).mean(axis=0) / (2 * width)
    assert np.max(abs(empty.imag)) < 2e-13 * np.max(abs(empty.real))
    assert np.max(abs(marked.imag)) <= 2e-13 * np.max(abs(marked.real))
    return empty.real, marked.real


def hermitian_exponential(h, time):
    energies, vectors = np.linalg.eigh(h)
    return (vectors * np.exp(-time * energies)) @ vectors.conj().T, energies


def ref_trace(h):
    matrix, energies = hermitian_exponential(h, duration)
    zero = np.zeros_like(matrix)
    empty, _ = determinant_and_yukawa(identity, matrix, zero,
                                     [np.eye(len(i)) for i, _ in blocks],
                                     duration * energies[:48].sum())
    return float(np.log(empty))


def scope_checks(momentum, log_ref, middle, middle_jet, e_late):
    # With link -I only the early Higgs sign is reversed; switching both
    # signs together is an exact global weak gauge change.
    h0 = momentum * gamma
    early, energies = hermitian_exponential(h0 + clifford_majorana + .35 * clifford_higgs, width)
    early_minus, _ = hermitian_exponential(h0 + clifford_majorana - .35 * clifford_higgs, width)
    sea = 2 * width * (energies[:48].sum() + e_late[:48].sum()) - log_ref
    plus = [np.eye(len(i)) for i, _ in blocks]
    minus = [weak_lift(-np.eye(2))[np.ix_(i, i)] for i, _ in blocks]
    rplus, yplus = determinant_and_yukawa(early, middle, middle_jet, plus, sea)
    rminus, yminus = determinant_and_yukawa(early, middle, middle_jet, minus, sea)
    rflip, yflip = determinant_and_yukawa(early_minus, middle, middle_jet, plus, sea)
    assert np.allclose([rminus, yminus], [rflip, yflip], atol=1e-15, rtol=1e-11)

    def lift(g):
        lam = np.exp(.5j * np.angle(np.linalg.det(g)))
        color = np.exp(-1j * np.angle(lam) / 3) * np.eye(3)
        rep = native.rep(lam, g / lam, color)
        return rep @ swap @ rep.conj() @ swap

    gl = np.exp(.21j) * source["quaternion"](np.array([.6, .8j]))
    gr = np.exp(-.38j) * source["quaternion"](np.array([.8j, .6]))
    l, r = lift(gl), lift(gr)
    actual_g = source["quaternion"](np.array([.6j, .8]))
    actual_link = weak_lift(actual_g)
    moved_link = r @ actual_link @ l.conj().T
    moved_early, _ = hermitian_exponential(
        h0 + clifford_majorana - 1j * gamma @ inner(gl @ np.array([0., .35])), width)
    moved_y = -1j * gamma @ inner(gr @ np.array([0., late]))
    moved_middle, moved_jet = expm_frechet(-2 * width * (h0 + clifford_majorana + moved_y),
                                           -2 * width * moved_y)
    err = max(np.max(abs(moved_early - l @ early @ l.conj().T)),
              np.max(abs(moved_middle - r @ middle @ r.conj().T)),
              np.max(abs(moved_jet - r @ middle_jet @ r.conj().T)))
    assert err < 1e-12
    one = determinant_and_yukawa(early, middle, middle_jet,
                                 [actual_link[np.ix_(i, i)] for i, _ in blocks], sea)
    changed = determinant_and_yukawa(moved_early, moved_middle, moved_jet,
                                     [moved_link[np.ix_(i, i)] for i, _ in blocks], sea)
    assert np.allclose(one, changed, atol=1e-15, rtol=1e-11)
    # On the same fixed history, sew the two time pieces with all 96 channels.
    late_block, _ = hermitian_exponential(h0 + clifford_majorana + late * clifford_higgs, width)
    whole = early @ actual_link.conj().T @ middle @ actual_link @ early
    half = late_block @ actual_link @ early
    sew_error = np.max(abs(whole - half.conj().T @ half))
    assert sew_error < 1e-12
    return dict(Z0_center_plus=float(rplus), Z0_center_minus=float(rminus),
                late_Y_center_plus=float(yplus / rplus), late_Y_center_minus=float(yminus / rminus),
                passive_U2_frame_error=float(err), full_96_sewing_error=float(sew_error))


def solve_metric(m, radial_order, order=49, axis_order=25):
    p = .25 * np.exp(-m)
    h0 = p * gamma
    log_ref = ref_trace(h0)
    late_h = h0 + clifford_majorana + late * clifford_higgs
    middle, middle_jet = expm_frechet(-2 * width * late_h, -2 * width * late * clifford_higgs)
    e_late = np.linalg.eigvalsh(late_h)
    angular_links, angular_weights = weak_arrays(order, axis_order)
    radial_nodes, radial_weights = leggauss(radial_order)
    u_grid = (radius_max * (1 + radial_nodes) / 2)**2
    measure = radial_weights * (1 + radial_nodes)**3 / 4
    assert abs(measure.sum() - 1) < 1e-13
    n = np.zeros((3, 4))
    log_boson_ref = -duration * len(gamma) * np.exp(-p * p)

    for u, mass in zip(u_grid, measure, strict=True):
        early, e_early = hermitian_exponential(h0 + clifford_majorana + np.sqrt(u) * clifford_higgs, width)
        sea = 2 * width * (e_early[:48].sum() + e_late[:48].sum())
        z, y = determinant_and_yukawa(early, middle, middle_jet, angular_links, sea - log_ref)
        # Candidate 1: declared local finite trace heat profile f(x)=exp(-x),
        # original and S evaluated at the same p and cutoff scale 1.
        # The two other rows have normalized central link laws and no such action.
        log_boson = -2 * width * (np.exp(-e_early**2).sum() + np.exp(-e_late**2).sum())
        scalar = np.array([np.exp(log_boson - log_boson_ref), 1., 1.])
        az = np.einsum("kij,ij->k", angular_weights, z) * scalar
        ay = np.einsum("kij,ij->k", angular_weights, y) * scalar
        n += mass * np.array([az, u * az, ay, u * ay]).T
    # At U=0 the early transfer commutes with the weak group; the three
    # normalized angular integrations all give the same fermionic factor.
    early0, e0 = hermitian_exponential(h0 + clifford_majorana, width)
    identity_links = [np.eye(len(i)) for i, _ in blocks]
    sea0 = 2 * width * (e0[:48].sum() + e_late[:48].sum())
    k0, _ = determinant_and_yukawa(early0, middle, middle_jet, identity_links, sea0 - log_ref)
    log_boson0 = -2 * width * (np.exp(-e0**2).sum() + np.exp(-e_late**2).sum())
    f0 = 2 / (np.pi**2 * radius_max**4) * k0 * np.array([np.exp(log_boson0 - log_boson_ref), 1, 1])
    z0, zu, zy, zuy = n.T
    assert min(z0) > 0
    connected = zuy / z0 - zu * zy / z0**2
    controls = scope_checks(p, log_ref, middle, middle_jet, e_late)
    controls["spectral_boson_S_factor"] = float(np.exp(log_boson_ref))
    return np.array([z0, f0, connected, zu / z0, zy / z0]).T, controls


names = ("finite_96_heat_profile_with_Haar", "SU2_heat_at_beta_2pi", "prior_A_at_same_strip_time")
fields = ("Z0_original_over_weight_matched_S", "F0_cartesian_early_H", "C_UY_disjoint_windows", "U_mean", "Y_mean")


if __name__ == "__main__":
    rows = {}
    for m in (-.001, 0., .001):
        values, checks = solve_metric(m, 20)
        rows[m] = values
        print(json.dumps({"m": m, "rows": {name: dict(zip(fields, vals, strict=True))
                                           for name, vals in zip(names, values.tolist(), strict=True)},
                          "spectral_Z0_over_common_bare_S": values[0, 0] * checks["spectral_boson_S_factor"],
                          "spectral_F0_over_common_bare_S": values[0, 1] * checks["spectral_boson_S_factor"],
                          "controls": checks}), flush=True)
    more, _ = solve_metric(0., 28)
    residual = np.max(abs(more - rows[0.]) / np.maximum(abs(more), 1e-16))
    print(json.dumps({"radial_refinement_largest_relative_residual": float(residual),
                      "weight_matched_log_Z0_first_metric_jet": ((np.log(rows[.001][:, 0]) - np.log(rows[-.001][:, 0])) / .002).tolist(),
                      "weight_matched_log_Z0_second_metric_jet": ((np.log(rows[.001][:, 0]) - 2 * np.log(rows[0.][:, 0]) + np.log(rows[-.001][:, 0])) / 1e-6).tolist(),
                      "scope": "Conditional finite 96-mode R032 comparison with actual Z3 only; no physical Q_Theta/gap/P_adm or original R113 boson measure."}), flush=True)
    assert residual < 1e-7
