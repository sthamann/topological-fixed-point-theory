"""Native C58 background Hessian on two disjoint collar windows."""

import contextlib
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from scipy.integrate import quad
from scipy.linalg import expm


repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
with contextlib.redirect_stdout(io.StringIO()):
    old = runpy.run_path(str(base / "check_tfpt_conditioned_r113_weights_20261110.py"))
width = 1 / 8
higgs = old["clifford_higgs"]
majorana = old["clifford_majorana"]
gamma = old["gamma"]
left = (width / 8, width / 4)
mirror = (3 * width / 4, 7 * width / 8)
shifted = (width / 2, 5 * width / 8)
shifted_reflected = (3 * width / 8, width / 2)
off_block = np.ones_like(gamma, dtype=bool)
for indices, _ in old["blocks"]:
    off_block[np.ix_(indices, indices)] = False
for matrix in (gamma, higgs, majorana):
    assert np.max(abs(matrix[off_block])) < 1e-14


def geometry(index, p, mass):
    grade = gamma[np.ix_(index, index)]
    h = p * grade + mass[np.ix_(index, index)]
    eig, u = np.linalg.eigh(h)
    occupied = u.conj().T[:, p * np.diag(grade) < 0]
    v = u.conj().T @ higgs[np.ix_(index, index)] @ u
    return eig - eig[0], occupied, v


def transfer_derivatives(eig, v, duration):
    diff = np.abs(eig[:, None] - eig)
    a = duration * diff
    ratio = np.ones_like(a)
    np.divide(-np.expm1(-a), a, out=ratio, where=a != 0)
    exponential = np.exp(-duration * eig)
    derivative = -duration * v * np.maximum(exponential[:, None], exponential) * ratio
    return exponential, derivative


def cross(eig, occupied, v, earlier, later):
    a, b = earlier
    c, d = later
    e, early = transfer_derivatives(eig, v, b - a)
    f, late = transfer_derivatives(eig, v, d - c)
    before = np.exp(-a * eig)[:, None] * occupied
    after = occupied.conj().T * np.exp(-(width - d) * eig)
    between = np.exp(-(c - b) * eig)[:, None]
    cap = occupied.conj().T @ (np.exp(-width * eig)[:, None] * occupied)
    de = after @ (f[:, None] * between * (early @ before))
    df = after @ (late @ (between * e[:, None] * before))
    mixed = after @ (late @ (between * (early @ before)))
    derivatives = np.linalg.solve(cap, mixed) - np.linalg.solve(cap, de) @ np.linalg.solve(cap, df)
    return float(np.trace(derivatives).real)


def direct(eig, occupied, v, earlier, later, delta):
    a, b = earlier
    c, d = later
    return sum(
        x * y * np.linalg.slogdet(
            (occupied.conj().T * np.exp(-(width-d)*eig))
            @ expm(-(d-c)*(np.diag(eig)+y*delta*v))
            @ (np.exp(-(c-b)*eig)[:,None] * expm(-(b-a)*(np.diag(eig)+x*delta*v)))
            @ (np.exp(-a*eig)[:,None]*occupied)
        )[1]
        for x in (-1,1) for y in (-1,1)
    )/(4*delta**2)


def free_cross(p, earlier, later):
    a, b = earlier
    c, d = later
    weight = np.trace(higgs @ higgs).real / 2
    q = 2 * abs(p)
    return weight * np.exp(-q*(c-b)) * np.expm1(-q*(b-a)) * np.expm1(-q*(d-c)) / q**2


def total(first, second, cutoff):
    mixed = 0.
    reference = 0.
    for sign, shift in ((1, .25), (-1, .75)):
        for p in sign * (np.arange(cutoff) + shift):
            reference += free_cross(p, first, second)
            for index, _ in old["blocks"]:
                mixed += cross(*geometry(index, p, majorana), first, second)
    return mixed, reference


for p in (.25, -.75):
    for index, _ in old["blocks"]:
        system = geometry(index, p, majorana)
        for window in (mirror, shifted):
            np.testing.assert_allclose(cross(*system, left, window),
                                       direct(*system, left, window, .01), atol=4e-9)
    free = sum(cross(*geometry(index, p, 0*majorana), left, shifted) for index, _ in old["blocks"])
    np.testing.assert_allclose(free, free_cross(p, left, shifted), atol=1e-14)

cases = {}
for key, window in (("mirror", mirror), ("shifted", shifted), ("shifted_reflected", shifted_reflected)):
    cases[key] = {str(cutoff): total(left, window, cutoff) for cutoff in (256, 512, 1024)}
    np.testing.assert_allclose(cases[key]["512"], cases[key]["1024"], atol=1e-10)

# OLD transported S sees odd H and has no Majorana; T sees even H and M_R.
a, free = cases["shifted"]["1024"]
a_reverse, free_reverse = cases["shifted_reflected"]["1024"]
ref_s = (free - free_reverse) / 2
ref_t = (a + a_reverse) / 2
shifted_relative = {"A": a, "reference_S": ref_s, "reference_T": ref_t,
                    "A_minus_S": a-ref_s, "A_minus_T": a-ref_t}

with contextlib.redirect_stdout(io.StringIO()):
    radial = runpy.run_path(str(base / "check_tfpt_c58_relative_radial_measure_20261110.py"))
finite_part = radial["log_cap_finite_part"]
base_log = finite_part(radial["base_masses"], 1024)
test_log = finite_part(radial["mass_squared"](.55), 1024)
subtraction_c = -24.70819538394082
logs = {}
for c, body in ((0., 0.), (subtraction_c, 0.), (subtraction_c, -subtraction_c)):
    def density(r):
        fermion = radial["relative_log"](r, 512)-c*r*r
        return 2*np.pi**2*r**3*np.exp(fermion-body*r*r)

    mass, _ = quad(density, 0., 8., epsabs=1e-13)
    logs[f"fermion_c={c},body_c={body}"] = float(mass)
assert abs(logs["fermion_c=0.0,body_c=0.0"]-logs[
    f"fermion_c={subtraction_c},body_c={-subtraction_c}"]) < 1e-13

print(json.dumps({
    "scope":"Original 96 native fixed-background Slater cap and pulled-back old S/T background maps; collar width, not physical OS time.",
    "width":width, "windows":{"left":left, "mirror":mirror, "shifted":shifted,
                              "shifted_reflected":shifted_reflected},
    "raw_cross_columns":["original_Majorana","massless_control"],
    "raw_cross":cases,"shifted_relative_background_jets":shifted_relative,
    "seam_even_projected_disjoint_pairs": {"A_minus_S":ref_t,"A_minus_T":0.},
    "constant_even_higgs_modewise_finite_parts": {
        "A_minus_S_H_zero_log":base_log,"A_minus_S_H_zero_factor":float(np.exp(base_log)),
        "A_minus_S_H_0_55_log":test_log,"A_minus_S_H_0_55_factor":float(np.exp(test_log)),
        "A_minus_T_both_log":0.,"A_minus_T_both_factor":1.,
    },
    "radial_integration_scheme_reallocation":logs,
    "local_counterterm_hessians_per_lambda_on_mirror":{"raw":0., "relative_S":-(left[1]-left[0]),
                                                        "relative_T":left[1]-left[0]},
    "local_counterterm_hessians_per_lambda_on_shifted":{"raw":0., "relative_S":0., "relative_T":0.}
}, indent=2))
