"""Original R032/v252 empty scalar: source inputs, midpoint bounds and tails."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.special import digamma

REPO_CHECK = ("experiments/theory-contracts/chat-triangulation-20261108/"
              "check_quarter_reference.py")
EULER = float(np.euler_gamma)


def source_module(repo):
    spec = importlib.util.spec_from_file_location("r032_original", repo / REPO_CHECK)
    source = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(source)
    return source


def half_spectrum(hermitian):
    eig = np.linalg.eigvalsh(hermitian)
    np.testing.assert_allclose(eig[:48], -eig[48:][::-1], atol=1e-10)
    return np.maximum(eig[48:], 0.0)


def midpoint_leading(masses, circle):
    masses = np.asarray(masses)
    result = np.zeros_like(masses)
    pos = masses > 0
    m = masses[pos]
    result[pos] = circle * m**2 * (np.log(4 * circle * m) + EULER - 0.5)
    return result


def sum_bounds(masses, circle, count):
    """Evaluate positive modal remainder with a certified analytic tail."""
    m = np.asarray(masses)
    p = ((np.arange(count) + .5) / (2 * circle))[:, None]
    energy = np.hypot(p, m)
    finite = (m**4 / (2 * p * (energy + p)**2)).sum(axis=0)
    tail = m**4 * (2 * circle)**3 / (16 * count**2)
    return finite, finite + tail


def verify(repo):
    src = source_module(repo)
    m, y, gamma = src.mass_matrix, src.higgs_matrix, src.grading
    time = float(src.time)
    circle = 1.0
    a = np.vdot(y, y).real / 2
    mu = half_spectrum(m)
    y_eig = half_spectrum(y)
    assert abs(a - np.sum(y_eig**2)) < 1e-10
    assert abs(np.trace(m @ y)) < 1e-12
    assert np.linalg.norm(gamma @ m + m @ gamma) < 1e-12
    assert np.linalg.norm(gamma @ y + y @ gamma) < 1e-12

    # Compare the original 96x96 tangential operator with the energy pairs.
    for amplitude in (0, .55, 2, 8):
        phi = m + amplitude * y
        masses = half_spectrum(phi)
        np.testing.assert_allclose(
            np.linalg.eigvalsh(.25 * gamma - 1j * gamma @ phi)[48:],
            np.hypot(.25, masses), atol=1e-11,
        )
        assert abs(np.sum(masses**2) - np.sum(mu**2) - a * amplitude**2) < 1e-10

    # Independent direct modal sums and analytic midpoint error.
    count = 2**17
    test = np.r_[0., 1e-3, .25, 1., 4., 16., 64., half_spectrum(m + 8*y)]
    lo, hi = sum_bounds(test, circle, count)
    base = midpoint_leading(test, circle)
    assert np.all(lo <= base + 1 / (16 * circle) + 1e-12)
    assert np.all(hi + 1e-12 >= base)
    assert lo[0] == hi[0] == 0.0

    finite0, upper0 = sum_bounds(mu, circle, count)
    rows = []
    for amp in (.55, 2., 8., 32., 128., 512., 2048.):
        masses = half_spectrum(m + amp * y)
        low, high = sum_bounds(masses, circle, count)
        analytic = time * (midpoint_leading(masses, circle).sum()
                           - midpoint_leading(mu, circle).sum())
        error = time * 48 / (16 * circle)
        low_bound = max(analytic - error, time * (low.sum() - upper0.sum()))
        high_bound = min(analytic + error, time * (high.sum() - finite0.sum()))
        assert low_bound <= high_bound
        rows.append([amp, low_bound, high_bound])

    # Same bosonic endpoints and geometry. The smooth even loop uses the
    # original radial field direction, H(s)=amp*sin(pi*s/T)*e2.
    profile_mass = time / 2
    # Coefficient at high p, calculated from unchanged h instead of a counterterm.
    u = np.linspace(0, time, 2049)
    x = .55 * np.sin(np.pi * u / time)
    big_p = 1e4
    energy_diff = []
    for field in x:
        n = half_spectrum(m + field * y)
        # Stable subtraction of p in each spectrum.
        energy_diff.append(np.sum(n**2/(np.hypot(big_p, n) + big_p)
                                  - mu**2/(np.hypot(big_p, mu) + big_p)))
    mode_contrast = np.trapezoid(energy_diff, u)
    expected = a * .55**2 * profile_mass
    np.testing.assert_allclose(2 * big_p * mode_contrast, expected, rtol=2e-7)

    # Fixed cutoff removeability on bounded histories is uniform, verified
    # separately from the large amplitude and from the circle radius limit.
    harmonic = digamma(count + .5) - digamma(.5)
    np.testing.assert_allclose(
        harmonic, np.sum(1/(np.arange(count)+.5)), atol=1e-13,
    )
    return {
        "circle": circle, "time": time, "half_trace_Y_squared": a,
        "half_trace_M_squared": float(np.sum(mu**2)),
        "leading_fixed_circle_coefficient": time * circle * a,
        "C58_scalar_coefficient_same_fixture_if_time_is_1_over_8": -circle*a/8,
        "same_endpoint_loop_log_cutoff_coefficient": -circle * expected,
        "direct_high_momentum_loop_coefficient": float(-2 * big_p * circle * mode_contrast),
        "rows_amp_log_factor_low_high_relative_H_zero": rows,
        "pointwise_tail_bound_cutoff": count,
        "midpoint_error_upper_per_mass": 1/(16*circle),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("repo", type=Path)
    options = parser.parse_args()
    print(json.dumps(verify(options.repo), indent=2))
