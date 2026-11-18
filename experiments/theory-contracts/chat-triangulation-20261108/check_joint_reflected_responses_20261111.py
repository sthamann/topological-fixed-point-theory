"""Separated Higgs, relative hypercharge and metric response on the fixed R032 history."""

import contextlib
from functools import cache
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from scipy.linalg import eigh, eigvalsh


repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo / "verification"))
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
sys.path.insert(0, str(base))
from check_mixed_higgs_eich import gauge_operators

source_file = base / "check_r113_dressed_open_20261109.py"
sys.argv = [str(source_file), str(base)]
with contextlib.redirect_stdout(io.StringIO()):
    source = runpy.run_path(str(source_file))

grading = source["grading"]
majorana = source["majorana"]
history = source["history"]
durations = source["durations"]
frames = source["frames"]
identity = np.eye(len(grading))
hypercharge = gauge_operators()[1]
charges = hypercharge.diagonal().real
assert np.allclose(hypercharge, np.diag(charges))
assert np.allclose(grading @ hypercharge, hypercharge @ grading)
assert np.allclose(majorana @ hypercharge, hypercharge @ majorana)

higgs_step = 1
link_before_step = 3
geometry_step = 4


@cache
def heat_step(i, higgs, geometry, reference):
    momentum = source["p"] * np.exp(-geometry)
    generator = momentum * grading
    if not reference:
        phi = source["inner_form"](history[i]) * (1 + higgs) + majorana
        generator = generator - 1j * grading @ phi
    assert np.max(np.abs(generator - generator.conj().T)) < 1e-13
    values, vectors = eigh(generator)
    width = durations[i]
    matrix = (vectors * np.exp(-width * values)) @ vectors.conj().T
    return matrix, float(width * values[:len(values) // 2].sum())


def reflected_log(higgs, angle, geometry, *, reference=False, gauge_control=False,
                  framed=False, link_step=link_before_step):
    twist = np.diag(np.exp(1j * angle * charges))
    matrix = identity.copy()
    sea = 0.
    for i in range(len(history)):
        heat, weight = heat_step(i, higgs if i == higgs_step else 0.,
                                 geometry if i == geometry_step else 0., reference)
        if gauge_control and i >= link_step:
            heat = twist @ heat @ twist.conj().T
        link = twist if i == link_step else identity
        if framed:
            if i:
                link = frames[i].conj().T @ link @ frames[i - 1]
            heat = frames[i].conj().T @ heat @ frames[i]
        matrix = heat @ link @ matrix
        sea += weight
    gram = matrix.conj().T @ matrix
    if framed:
        sign, log_det = np.linalg.slogdet(identity + gram)
        assert abs(sign - 1) < 1e-12
    else:
        eigenvalues = eigvalsh(gram)
        assert eigenvalues[0] > 0
        log_det = np.log1p(eigenvalues).sum()
    return float(2 * sea + log_det)


def relative_log(higgs, angle, geometry, **options):
    return (reflected_log(higgs, angle, geometry, **options)
            - reflected_log(0., angle, geometry, reference=True, **options))


def contrast(step, **options):
    even_gauge = {
        (h, g): (relative_log(h, step, g, **options)
                 + relative_log(h, -step, g, **options)
                 - 2 * relative_log(h, 0., g, **options))
        for h in (-step, step) for g in (-step, step)
    }
    numerator = sum(np.sign(h * g) * value for (h, g), value in even_gauge.items())
    return float(numerator), float(numerator / (4 * step**4))


rows = []
for step in (0.12, 0.06, 0.03, 0.015):
    numerator, value = contrast(step)
    rows.append({"step": step, "raw_difference": numerator, "mixed_response": value})

step = 0.06
_, framed = contrast(step, framed=True)
_, gauge = contrast(step, gauge_control=True)
_, moved = contrast(step, link_step=2)
baseline_gauge = relative_log(step, step, -step, gauge_control=True)
baseline_untwisted = relative_log(step, 0., -step)
assert abs(baseline_gauge - baseline_untwisted) < 1e-11
assert abs(framed - rows[1]["mixed_response"]) < 1e-7
assert abs(gauge) < 1e-7
assert abs(moved - rows[1]["mixed_response"]) > 1e-5
assert abs(rows[2]["mixed_response"] - rows[1]["mixed_response"]) < 0.02 * abs(rows[2]["mixed_response"])

print(json.dumps({"fiber_dimension": len(grading), "durations": durations,
                  "higgs_step_0_based": higgs_step, "link_before_step_0_based": link_before_step,
                  "geometry_step_0_based": geometry_step,
                  "reference": "S on the same R032 history, relative reflected log trace",
                  "rows": rows, "framed_at_0.06": framed,
                  "pure_gauge_at_0.06": gauge,
                  "link_before_step_2_at_0.06": moved,
                  "pure_gauge_value_residual": abs(baseline_gauge - baseline_untwisted)}, indent=2))
