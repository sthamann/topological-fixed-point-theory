"""Evaluate radial source probes on the existing full96 Calderón pullback."""
from contextlib import redirect_stdout
from pathlib import Path
import io
import json
import runpy
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss

repo = Path(sys.argv[1]).resolve()
original = repo / "experiments/theory-contracts/recursive-mirrored-source-20261103"
with redirect_stdout(io.StringIO()):
    source = runpy.run_path(str(original / "check_tfpt_calderon_volume_20261110.py"))

metric = source["metric"]
circle = source["circle"]
volume = source["volume"]
packed = source["packed"]
alpha, time = source["alpha"], source["t_ref"]

# Evaluate the same S reference on the retained momenta. exp(-x) is
# the existing diagnostic slice of a positive proper-time profile.
def radial_weight(radius, modes):
    masses2 = np.concatenate([
        np.linalg.eigvalsh(mass + radius * vertices[2]) ** 2
        for _, mass, vertices in packed
    ])
    squared = modes[:, None] ** 2 + masses2[None, :]
    heat_trace = np.exp(-squared).sum()
    thermal_trace = np.logaddexp(0, -time * np.sqrt(squared)).sum()
    return np.exp(thermal_trace - heat_trace)


nodes, weights = leggauss(40)
angles = (nodes + 1) * np.pi / 4
xs, xweights = np.tan(angles), weights * np.pi / 4 / np.cos(angles) ** 2
scale = 1.0
observable = np.exp(-xs**2)
regulator = np.exp(-scale * xs**2)

# On R=N*x the limiting metric uses the unmodified principal
# coefficients with p/R in [-1/x, 1/x]. The Majorana term scales as 1/N.
angles_inside = np.arctan(1 / xs)
inside_nodes, inside_weights = leggauss(56)
limiting = []
for x, angle in zip(xs, angles_inside, strict=True):
    positions = angle * inside_nodes
    positions_weights = angle * inside_weights / np.cos(positions) ** 2
    metric_limit = sum(
        weight * metric(p, 1, sign=True, mr=False)
        for p, weight in zip(np.tan(positions), positions_weights, strict=True)
    )
    radial = 1 / np.linalg.inv(metric_limit)[2, 2]
    limiting.append((x * volume(metric_limit), np.sqrt(radial / x)))


def laplace_ratios(radial_densities):
    mass = (xweights * regulator)[:, None] * np.asarray(radial_densities)
    return (observable[:, None] * mass).sum(axis=0) / mass.sum(axis=0)


def batched_circle(radius, momenta):
    # Same divided differences as the original metric, batched over p.
    gram = np.zeros((4, 4))
    length = source["ell"]
    for grading, mass, vertices in packed:
        values, vectors = np.linalg.eigh(
            momenta[:, None, None] * grading + radius * vertices[2] + mass
        )
        axes = vectors.swapaxes(-1, -2).conj()
        variations = axes[:, None] @ vertices[None] @ vectors[:, None]
        z = length * values
        exp_z = np.exp(-np.abs(z))
        sech = 2 * exp_z / (1 + exp_z**2)
        tanh = np.tanh(z)
        dx = values[:, :, None] - values[:, None, :]
        coincident = np.abs(dx) < 1e-11
        for function, derivative in (
            (tanh, length * sech**2), (sech, -length * sech * tanh)
        ):
            df = function[:, :, None] - function[:, None, :]
            average_derivative = (derivative[:, :, None] + derivative[:, None, :]) / 2
            divided = np.divide(df, dx, out=average_derivative, where=~coincident)
            terms = variations * divided[:, None]
            gram += 0.5 * np.einsum("naij,nbij->ab", terms.conj(), terms).real
    return gram


np.testing.assert_allclose(
    batched_circle(.55, np.arange(-3, 4) + alpha), circle(.55, 3, alpha),
    rtol=1e-10, atol=1e-10
)
limit = laplace_ratios(limiting)
rows = []
for cutoff in (7, 31, 127, 511):
    momenta = np.arange(-cutoff, cutoff + 1) + alpha
    masses = []
    unweighted_masses = []
    for x in xs:
        radius = cutoff * x
        graph_metric = batched_circle(radius, momenta)
        radial = 1 / np.linalg.inv(graph_metric)[2, 2]
        density = np.array((x**3 * cutoff**2 * volume(graph_metric), np.sqrt(cutoff * radial)))
        weight = radial_weight(radius, momenta)
        masses.append(weight * density)
        unweighted_masses.append(density)
    actual = laplace_ratios(masses)
    unweighted = laplace_ratios(unweighted_masses)
    rows.append({
        "cutoff": cutoff,
        "weighted_laplace": actual.tolist(),
        "bare_laplace": unweighted.tolist(),
        "difference_from_limit": (actual - limit).tolist(),
    })

errors = np.abs([row["difference_from_limit"] for row in rows])
assert np.all(errors[1:] < errors[:-1]), (rows, limit)
assert np.all(errors[-1] < 0.035), (rows, limit)
assert np.max(np.abs(np.asarray(rows[-1]["bare_laplace"]) - limit)) < 0.003
assert errors[1].min() > 0.1  # The spectral factor still changes small cutoffs substantially.
assert np.all((limit > 0) & (limit < 1))
assert abs(limit[0] - limit[1]) > 0.10
# The actual cutoff weight tends to a positive constant at large Higgs norm.
all_modes = np.arange(-31, 32) + alpha
assert abs(radial_weight(64, all_modes) - 1) < 1e-25
assert abs(radial_weight(0, all_modes) - 1) > 0.1
# In the full-tower-first regime the gamma law follows independently
# by integrating the predicted radial power against exp(-x²).
from scipy.integrate import quad
powers = (1, -0.5)
gamma_controls = []
for power in powers:
    mass = quad(lambda x: x**power * np.exp(-x*x), 0, np.inf)[0]
    mark = quad(lambda x: x**power * np.exp(-2*x*x), 0, np.inf)[0]
    exact = 2**(-(power + 1) / 2)
    assert abs(mark / mass - exact) < 1e-10
    gamma_controls.append(exact)

print(json.dumps({
    "input": "existing R032/v252 full96 constant chart and R032 free S reference",
    "scope": "calderon graph base and the known spectral and sea factors only",
    "metric_order": ["four_real_components", "intrinsic_weak_angle_quotient"],
    "kappa_times_cutoff_squared": scale,
    "probe": "exp(-abs(H)^2 / cutoff^2)",
    "joint_limit": limit.tolist(),
    "native_controls": rows,
    "tower_first_probe_exp_minus_kappa_U": gamma_controls,
    "primitive_F1_computed": False,
}, indent=2))
