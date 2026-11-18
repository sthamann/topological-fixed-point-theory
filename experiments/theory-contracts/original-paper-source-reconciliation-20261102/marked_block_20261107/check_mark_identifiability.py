"""Corroborate the analytic mark separation with the existing v252 fixture."""

import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np
from scipy.special import expit


def run(repo):
    sys.path.insert(0, str(repo / "verification"))
    import v252_full_finite_triple as v252

    # Same seed and arrays as v252.run: structural controls, not physical inputs.
    rng = np.random.default_rng(252)
    n = v252.NG
    ys = [rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)) for _ in range(4)]
    m0 = rng.normal(size=(n, n))
    mr = m0 + m0.T
    zero = np.zeros_like(mr)
    gamma = v252.grading()
    p, time = 0.25, 4 * np.pi
    h0 = p * gamma - 1j * gamma @ v252.full_D(*([zero] * 4), mr)
    k = -1j * gamma @ v252.full_D(*ys, zero)
    assert np.allclose(h0, h0.conj().T) and np.allclose(k, k.conj().T)
    assert np.linalg.norm(mr) > 0 and np.linalg.norm(k) > 0

    def mark(radius):
        ev, v = np.linalg.eigh(h0 + radius * k)
        diag = np.einsum("ij,ij->j", v.conj(), k @ v).real
        return radius * np.dot(expit(-time * ev), diag)

    def f(radius):
        ev = np.linalg.eigvalsh(h0 + radius * k)
        return np.logaddexp(0, -time * ev).sum()

    radii = np.array([0.25, 0.5, 1.0, 2.0])
    u = radii**2
    y = np.array([mark(r) for r in radii])
    delta = 1e-4
    other = -radii * np.array([(f(r + delta) - f(r - delta)) / (2 * delta) for r in radii]) / time
    assert np.allclose(y, other, rtol=1e-7, atol=1e-7)
    base = np.stack((np.ones_like(u), u, y))
    minors = np.array([(-1)**i * np.linalg.det(base[:, columns])
                       for i, columns in enumerate(itertools.combinations(range(4), 3))])
    # combinations omit the columns in reverse order, giving a common sign.
    direction = minors[::-1] / np.max(np.abs(minors))
    weights = np.full(4, 0.25)
    plus, minus = weights + direction / 8, weights - direction / 8
    assert min(plus.min(), minus.min()) >= 0.125
    assert np.linalg.norm(base @ (plus - minus), ord=np.inf) < 1e-12
    mixed_difference = float((u * y) @ (plus - minus))
    assert abs(mixed_difference) > 1e-4

    slope = -np.minimum(np.linalg.eigvalsh(k), 0).sum()
    slope_sv = 2 * sum(m * np.linalg.svd(a, compute_uv=False).sum()
                       for m, a in zip((3, 3, 1, 1), ys))
    assert np.isclose(slope, slope_sv, rtol=1e-13, atol=1e-13)
    large_radius = 1024.0
    slope_error = abs(-mark(large_radius) / large_radius - slope)
    bound = 4 * len(k) * (time * np.linalg.norm(h0, 2) + np.log(2)) / (time * large_radius)
    assert slope_error < bound
    eig, vec = np.linalg.eigh(k)
    negative = vec[:, eig < -1e-10]
    swap = v252.Kswap()
    color_checks = []
    for g in (np.diag(np.exp(1j * np.array([0.37, 0.91, -1.28]))),
              np.roll(np.eye(3), 1, axis=0)):
        pi = v252.rep(1, np.eye(2), g.conj())
        color = pi @ swap @ pi.conj() @ swap
        comm = max(np.linalg.norm(color @ h - h @ color, ord=np.inf) for h in (h0, k))
        restriction = negative.conj().T @ color @ negative
        error = abs(np.linalg.det(restriction) - 1)
        assert comm < 1e-12 and error < 1e-12
        color_checks.append({"max_commutator_norm": float(comm),
                             "negative_exterior_character_error": float(error)})
    return {
        "scope": "v252 structural fixture; float corroboration, no physical source evaluation",
        "p": p, "time": time, "radii": radii.tolist(), "raw_marks": y.tolist(),
        "max_log_trace_derivative_disagreement": float(np.max(np.abs(y - other))),
        "minimum_positive_weight": float(min(plus.min(), minus.min())),
        "three_moment_difference": (base @ (plus - minus)).tolist(),
        "mixed_difference": mixed_difference,
        "negative_spectral_slope": float(slope), "singular_value_slope": float(slope_sv),
        "slope_error_at_1024": float(slope_error), "slope_bound_at_1024": float(bound),
        "real_SU3_action_controls": color_checks,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("repo", type=Path)
    print(json.dumps(run(parser.parse_args().repo.resolve()), indent=2))
