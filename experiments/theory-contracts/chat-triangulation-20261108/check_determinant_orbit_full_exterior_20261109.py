"""Exact finite determinant readout of the conditional R032 gauge-pair test."""
import contextlib
import io
import json
import runpy
import sys
from pathlib import Path

import numpy as np

repo = Path(sys.argv[1]).resolve()
own = Path(__file__).resolve().parent / "check_determinant_lift_higgs_path_20261109.py"
args = sys.argv
try:
    sys.argv = [str(own), str(repo)]
    with contextlib.redirect_stdout(io.StringIO()):
        src = runpy.run_path(str(own))
finally:
    sys.argv = args

strip, history, widths = (src[key] for key in ("strip", "history", "durations"))
original = src["original"]
gauge, identity, idx_r = (src[key] for key in ("anchor", "identity", "idx_r"))
assert np.linalg.norm(gauge[:, idx_r] + identity[:, idx_r]) < 1e-14

# Trace over the original exterior algebra, with one OUTER + projector.
# If B is the one-particle whole transfer, the charged diagonal coefficient
# of Tr Γ(B) is det(I+B)-det((I+B) with the right row/column omitted).
def determinant(x):
    sign, logabs = np.linalg.slogdet(x)
    return sign * np.exp(logabs)


def read(b):
    empty = []
    charged = []
    for x in (b, b @ gauge):
        z = determinant(identity + x)
        minor = determinant(np.delete(np.delete(identity + x, idx_r, axis=0), idx_r, axis=1))
        covariance = np.linalg.solve(identity+x, x)
        assert abs(z-minor-z*covariance[idx_r,idx_r]) < 1e-10*max(abs(z),abs(minor))
        empty.append(z)
        charged.append(z - minor)
    # G^s anticommutes with the odd right creation operator when s=1.
    return (empty[0]+empty[1])/2, (charged[0]-charged[1])/2


rows = {}
for m in (-.001, .001):
    k,half_sea = strip(history[:3], .25*np.exp(-m), widths[:3])
    kg = gauge @ k @ gauge
    results = []
    for half in (k, kg):
        b = half.conj().T @ half
        twisted = half.conj().T @ gauge @ half
        results.append(np.array([read(b), read(twisted)]))
    assert np.max(np.abs(results[0]-results[1])) < 1e-11 * np.max(np.abs(results[0]))
    paired = sum(results) / 2
    assert np.max(np.abs(paired.imag)) < 1e-11 * np.max(np.abs(paired.real))
    plain, marked = paired.real
    ratios = marked/plain
    assert np.all(np.abs(ratios) <= 1 + 1e-11)
    assert 1 - ratios[0] > .001
    assert 1 + ratios[1] > .001
    # Keep the same original S reference used for the other R032 tests.
    reference, ref_sea = original["direct_whole"]("S", m)
    assert np.max(np.abs(reference@gauge-gauge@reference)) < 1e-12
    ref_plain = np.array(read(reference))
    ref_marked = np.array(read(reference@gauge))
    assert np.max(np.abs(ref_marked - ref_plain * [1, -1])) < 1e-12*ref_plain[0].real
    relative = np.exp(2*half_sea-ref_sea)*paired.real/ref_plain[0].real
    rows[str(m)] = {"empty_marked_over_plain": float(ratios[0]),
                    "charged_eR_diagonal_marked_over_plain": float(ratios[1]),
                    "normalized_charged_eR_plus_sector_weight": float((1+ratios[1])/2),
                    "plain_relative_empty_and_charged": relative[0].tolist(),
                    "marked_relative_empty_and_charged": relative[1].tolist(),
                    "S_reference_marked_over_plain": (ref_marked/ref_plain).real.tolist()}

print(json.dumps({
    "scope": "conditional R032 gauge pair; whole finite exterior trace and outer candidate gauge-plus projection only",
    "actual_Q_Theta_and_Padm_used": False,
    "two_metrics": rows,
}, indent=2))
