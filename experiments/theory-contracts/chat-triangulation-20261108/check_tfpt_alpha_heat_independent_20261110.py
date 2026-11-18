"""Independent full-matrix disk quadrature for the conditional R032 weight test."""

import argparse
import importlib.util
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import expm_frechet


parser = argparse.ArgumentParser()
parser.add_argument("repo", type=Path)
parser.add_argument("--metrics", nargs="+", type=float, default=[0.0])
parser.add_argument("--radial-order", type=int, default=20)
args = parser.parse_args()
repo = args.repo.resolve()
sys.path.insert(0, str(repo / "verification"))
from v3_em_alpha import make_F
from tfpt_constants import c3

spec = importlib.util.spec_from_file_location(
    "original_quarter", repo / "experiments/theory-contracts/chat-triangulation-20261108/check_quarter_reference.py"
)
original = importlib.util.module_from_spec(spec)
spec.loader.exec_module(original)
native = original.native
grading = original.grading
swap = native.Kswap()
identity = np.eye(native.NH, dtype=complex)
right = native.rep(1, np.zeros((2, 2)), np.zeros((3, 3)))
neutral_higgs = native.rep(0, np.eye(2), np.zeros((3, 3)))
directed = neutral_higgs @ (original.higgs_matrix @ right - right @ original.higgs_matrix)
hermitian_higgs = directed + directed.conj().T
higgs = -1j * grading @ (hermitian_higgs + swap @ hermitian_higgs.conj() @ swap)
majorana = -1j * grading @ original.mass_matrix
width, outer_radius = 0.2, 0.55
duration = 4 * width
omega = np.exp(2j * np.pi / 3)
color = np.array([(-1 if i < native.NP else 1) if native.P[i % native.NP][0] in ("u", "d") else 0
                  for i in range(native.NH)])
phase = omega ** color

with mp.workdps(60):
    alpha = mp.findroot(make_F(c3, 41), (mp.mpf(".007"), mp.mpf(".008")))
    g2 = 4 * mp.pi * alpha / mp.mpf(".375")
    heat_t = float(g2 / 2)
    print(json.dumps({"alpha_inverse": float(1 / alpha), "g_squared": float(g2), "hypothesis_heat_time": heat_t}), flush=True)


def lift(q):
    a = native.rep(1, q, np.eye(3))
    return a @ swap @ a.conj() @ swap


def disk_rule():
    # Haar SU2 induces uniform planar measure on z in the disk, and uniform
    # conditional phase of the other complex entry. The reflected words here
    # are invariant under that phase by the fixed Higgs stabilizer.
    disk_nodes, disk_weights = leggauss(25)
    azimuths = 2 * np.pi * np.arange(97) / 97
    z = (np.sqrt((1 + disk_nodes) / 2)[:, None] * np.exp(1j * azimuths)[None, :]).ravel()
    b = np.sqrt(1 - abs(z)**2)
    area_weights = np.repeat(disk_weights / (2 * len(azimuths)), len(azimuths))
    # Truncate only beyond the degree visible in a reflected full exterior word.
    # Each arm sees the one link once, and each of the 24 doublets has exterior
    # degree at most 1. Products including the central character have degree <=96.
    x = z.real
    characters = [np.ones_like(x), 2*x]
    while len(characters) < 49:
        characters.append(2 * x * characters[-1] - characters[-2])
    spin = np.arange(49) / 2
    casimir = spin * (spin + 1)
    transforms = np.array([
        np.exp(-heat_t * casimir),
        spin == 0,
        np.exp(-2 * np.pi * casimir),
        # A's jump character at spin 25 is above the observable bandwidth.
        np.exp(-duration * (casimir / 1000 + (spin > 0))),
    ])
    weights = ((transforms * (2 * spin + 1)) @ characters) * area_weights
    assert np.max(abs(weights.sum(axis=1) - 1)) < 2e-13
    # A real-linear expansion avoids thousands of calls into the native rep.
    unit_i = lift(np.diag([-1j, 1j]))
    unit_b = lift(np.array([[0, 1], [-1, 0]], complex))
    unit_minus = lift(-np.eye(2))
    doublets = (identity - unit_minus) / 2
    singlets = identity - doublets
    links = (singlets[None] + x[:, None, None] * doublets[None]
             + z.imag[:, None, None] * (unit_i - singlets)[None]
             + b[:, None, None] * (unit_b - singlets)[None])
    for k in (0, len(z) // 2, len(z) - 1):
        q = np.array([[z[k].conjugate(), b[k]], [-b[k], z[k]]])
        assert np.max(abs(links[k] - lift(q))) < 2e-15
    return links, weights


def heat(h, time):
    energy, basis = np.linalg.eigh(h)
    return (basis * np.exp(-time * energy)) @ basis.conj().T, energy[:native.NP].sum()


def exterior_trace(matrix, derivative=None):
    """Full 96-dimensional determinant, with one common color average."""
    values = []
    for twist in (np.ones_like(phase), phase):
        mat = identity + twist[:, None] * matrix
        sign, logmodulus = np.linalg.slogdet(mat)
        det = sign * np.exp(logmodulus)
        if derivative is not None:
            log_derivative = np.linalg.solve(mat, twist[:, None] * derivative).trace(axis1=-2, axis2=-1)
            det = np.stack((det, det * log_derivative), axis=-1)
        values.append(det)
    # The conjugate twist is the complex conjugate contribution by Hermiticity
    # and color conservation. The untwisted contribution is real.
    return (values[0].real + 2 * values[1].real) / 3


def solve(metric, links, angular, radial_order):
    p = .25 * np.exp(-metric)
    ref, ref_sea = heat(p * grading, duration)
    log_ref = np.log(exterior_trace(ref)) + duration * ref_sea
    h_late = p * grading + majorana + outer_radius * higgs
    h_early_zero = p * grading + majorana
    middle, middle_derivative = expm_frechet(-2 * width * h_late, -2 * width * outer_radius * higgs)
    late_sea = np.linalg.eigvalsh(h_late)[:native.NP].sum()
    nodes, weights = leggauss(radial_order)
    radial_mass = weights * (1 + nodes)**3 / 4
    assert abs(sum(radial_mass) - 1) < 1e-14
    moment = np.zeros((4, 4))
    for r, mass in zip(outer_radius * (1 + nodes) / 2, radial_mass, strict=True):
        early, sea = heat(h_early_zero + r * higgs, width)
        weighted = np.zeros((4, 2))
        for start in range(0, len(links), 48):
            end = min(start + 48, len(links))
            k = links[start:end] @ early
            adjoint = k.conj().swapaxes(-1, -2)
            value = exterior_trace(adjoint @ middle @ k, adjoint @ middle_derivative @ k)
            weighted += angular[:, start:end] @ value
        empty, marked = weighted.T * np.exp(2 * width * (sea + late_sea) - log_ref)
        marked = -marked / (2 * width)
        moment += mass * np.array([empty, r*r*empty, marked, r*r*marked]).T
    early0, sea0 = heat(h_early_zero, width)
    assert max(np.linalg.norm(early0 @ links[k] - links[k] @ early0) for k in (0, 1200, len(links)-1)) < 1e-13
    value0 = exterior_trace(early0 @ middle @ early0) * np.exp(2 * width * (sea0 + late_sea) - log_ref)
    density0 = 2 * value0 / (np.pi**2 * outer_radius**4)
    z0, zu, zy, zuy = moment.T
    rows = zip(("alpha_heat", "haar", "heat_2pi", "prior_A"),
               zip(z0, zu/z0, zy/z0, zuy/z0 - zu*zy/z0**2, strict=True), strict=True)
    return dict(metric=metric, radial_order=radial_order,
                same_F0=density0, same_log_S_reference=log_ref,
                rows={k: dict(zip(("Z0", "U_mean", "Y_mean", "C_UY"), row, strict=True)) for k, row in rows},
                scope="Conditional original R032 96-mode fixed-field chart and color projection. Heat time from alpha is a tested hypothesis.")


if __name__ == "__main__":
    links, angular = disk_rule()
    for metric in args.metrics:
        print(json.dumps(solve(metric, links, angular, args.radial_order)), flush=True)
