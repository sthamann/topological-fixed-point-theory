"""Derive the weak first response by the original NS doubled Green inverse."""

import json
import sys
from pathlib import Path

import numpy as np

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo / "verification"))
import v252_full_finite_triple as v

D = v.full_D(np.diag([1, 2, 3]), np.diag([4, 5, 6]), np.diag([2, 3, 4]),
             np.diag([1, 3, 5]), np.array([[2, 1, 0], [1, 3, 1], [0, 1, 4]]))
T = v.rep(0, 1j * np.array([[0, 1], [1, 0]]), np.zeros((3, 3)))
J = v.Kswap()
Q = -1j * (T + J @ T.conj() @ J)
gamma = v.grading()
V = Q @ gamma
mass = -1j * gamma @ D
weak_neutral = np.eye(96) - Q @ Q
s, arc = 1 / 8, np.pi / 4
cutoff = 384
nodes, weights = np.polynomial.legendre.leggauss(64)
u = (nodes + 1) * s / 2
weights = weights * s / 2


def green_nonnegative(lambdas, normal_time):
    return np.exp(-normal_time[:, None] * lambdas - np.logaddexp(0, -2 * s * lambdas))


def divided_difference(lambdas, values, derivative):
    delta = lambdas[:, None] - lambdas[None, :]
    close = np.abs(delta) < 1e-8
    out = np.empty_like(delta)
    np.divide(values[:, None] - values[None, :], delta, out=out, where=~close)
    midpoint = (lambdas[:, None] + lambdas[None, :]) / 2
    out[close] = derivative(midpoint[close])
    return out


def norm(a):
    return float(np.linalg.norm(a))


results = []
for alpha in (0.25, 0.75):
    X = np.zeros_like(D)
    Xprime_green = np.zeros_like(D)
    Xprime_compare = np.zeros_like(D)
    max_mode_error = 0.0
    max_frozen_projection_norm = 0.0
    neutral_compression_defect = None
    for m in range(-cutoff, cutoff + 1):
        p = m + alpha
        lambdas, U = np.linalg.eigh(p * gamma + mass)
        moved = U.conj().T @ V @ U
        positive = green_nonnegative(lambdas, np.array([0, s, 2 * s]))
        A, cross, C = positive

        # -G V G on the actual anti-periodic double. When the first
        # factor crosses the normal cut its Green value has a minus sign.
        def contraction(left, right):
            return green_nonnegative(lambdas, left).T @ (
                weights[:, None] * green_nonnegative(lambdas, right)
            )

        gx = -contraction(s - u, u) + contraction(2 * s - u, s + u)
        ga = contraction(2 * s - u, u) + contraction(s - u, s + u)
        Xgreen = gx * moved
        Agreen = ga * moved

        def cross_function(z):
            e = np.exp(-s * np.abs(z))
            return e / (1 + e**2)

        compare = divided_difference(lambdas, cross,
                                     lambda z: -s * np.tanh(s * z) * cross_function(z))
        Xcompare = compare * moved
        Acompare = divided_difference(lambdas, A,
                                      lambda z: s / (2 * np.cosh(s * z)**2)) * moved
        max_mode_error = max(max_mode_error, norm(Xgreen - Xcompare), norm(Agreen - Acompare))
        if m == 0:
            physical = np.block([[np.diag(A), np.diag(cross)],
                                 [np.diag(cross), np.diag(C)]])
            tangent = np.block([[Agreen, Xgreen], [Xgreen, -Agreen]])
            max_frozen_projection_norm = norm(physical @ tangent @ physical)
            Cf = (U * C) @ U.conj().T
            Xf = (U * cross) @ U.conj().T
            P = weak_neutral
            removed = np.eye(96) - P
            compressed_C = P @ Cf @ P
            compressed_X = P @ Xf @ P
            defect = compressed_C - compressed_C @ compressed_C - compressed_X @ compressed_X
            leakage = P @ Cf @ removed @ Cf @ P + P @ Xf @ removed @ Xf @ P
            assert norm(defect - leakage) < 1e-13
            assert np.linalg.eigvalsh(defect).min() > -1e-13
            neutral_compression_defect = norm(defect)
            assert neutral_compression_defect > 0.1

        endpoint_phase = np.exp(-1j * p * arc) / (2 * np.pi)
        X += endpoint_phase * ((U * cross) @ U.conj().T)
        Xprime_green += endpoint_phase * (U @ Xgreen @ U.conj().T)
        Xprime_compare += endpoint_phase * (U @ Xcompare @ U.conj().T)

    dressed = Xprime_green - 1j * arc * Q @ X
    assert norm(Xprime_green - Xprime_compare) < 1e-12
    assert max_mode_error < 1e-12
    assert max_frozen_projection_norm < 1e-13
    assert abs(norm(dressed) - 0.0002582371501243) < 1e-12
    results.append({
        "quarter": alpha,
        "max_single_mode_derivative_comparison_error": max_mode_error,
        "endpoint_green_vs_frechet_norm": norm(Xprime_green - Xprime_compare),
        "endpoint_green_first_response_norm": norm(Xprime_green),
        "endpoint_dressed_first_response_norm": norm(dressed),
        "frozen_cauchy_compression_tangent_norm_at_m_0": max_frozen_projection_norm,
        "weak_neutral_reader_squared_cross_defect_norm_at_m_0": neutral_compression_defect,
    })

print(json.dumps({"normal_double": "NS antiperiodic", "arc": "pi/4",
                  "collar_width": "1/8", "absolute_mode_cutoff": cutoff,
                  "quadrature_nodes_per_interval": len(nodes), "checks": results}, indent=2))
