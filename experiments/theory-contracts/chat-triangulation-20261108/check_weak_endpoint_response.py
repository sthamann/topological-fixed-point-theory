"""First weak derivative of the existing v252 two-face kernel at the Fold endpoints."""

import json
import sys
from pathlib import Path

import numpy as np

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo / "verification"))
import v252_full_finite_triple as v

# The actual data used in ACTUAL-WEAK-CONNECTION-AND-JOINT-SOURCE.
D = v.full_D(np.diag([1, 2, 3]), np.diag([4, 5, 6]), np.diag([2, 3, 4]),
             np.diag([1, 3, 5]), np.array([[2, 1, 0], [1, 3, 1], [0, 1, 4]]))
T = v.rep(0, 1j * np.array([[0, 1], [1, 0]]), np.zeros((3, 3)))
J = v.Kswap()
Q = -1j * (T + J @ T.conj() @ J)
gamma = v.grading()
h1 = Q @ gamma
h0_mass = -1j * gamma @ D

# exp(3 pi i Q_em) in the actual particle/antiparticle order.
three_q = {"u": 2, "d": -1, "e": -3, "nu": 0}
signs = [1 if three_q[label] % 2 == 0 else -1 for label, *_ in v.P]
em_flip = np.diag(signs + signs)
assert np.array_equal(em_flip @ D @ em_flip, D)
assert np.array_equal(em_flip @ gamma @ em_flip, gamma)
assert np.array_equal(em_flip @ Q @ em_flip, -Q)

s = 1 / 8  # Existing C58 collar diagnostic width.
ell = np.pi / 4  # The two actual Fold endpoints.
cutoffs = (256, 384)
rows = []


def transfer(eigenvalues):
    exp = np.exp(-s * np.abs(eigenvalues))
    return exp / (1 + exp**2)


def mode(p):
    lambdas, vectors = np.linalg.eigh(p * gamma + h0_mass)
    values = transfer(lambdas)
    delta = lambdas[:, None] - lambdas[None, :]
    close = np.abs(delta) < 1e-8
    jacobian = np.empty_like(delta)
    np.divide(values[:, None] - values[None, :], delta, out=jacobian, where=~close)
    midpoint = (lambdas[:, None] + lambdas[None, :]) / 2
    jacobian[close] = -s * np.tanh(s * midpoint[close]) * transfer(midpoint[close])
    moved = vectors.conj().T @ h1 @ vectors
    return ((vectors * values) @ vectors.conj().T,
            vectors @ (jacobian * moved) @ vectors.conj().T)


def norm(m):
    return float(np.linalg.norm(m))


for alpha in (0.25, 0.75):
    X = np.zeros_like(D)
    Xprime = np.zeros_like(D)
    X_short = np.zeros_like(D)
    Xprime_short = np.zeros_like(D)
    for m in range(-cutoffs[1], cutoffs[1] + 1):
        p = m + alpha
        kernel, derivative = mode(p)
        phase = np.exp(-1j * p * ell) / (2 * np.pi)
        X += phase * kernel
        Xprime += phase * derivative
        if abs(m) <= cutoffs[0]:
            X_short += phase * kernel
            Xprime_short += phase * derivative

    Wprime_X = 1j * ell * Q @ X
    dressed_prime = Xprime - Wprime_X
    assert norm(em_flip @ X @ em_flip - X) < 1e-10
    assert norm(em_flip @ Xprime @ em_flip + Xprime) < 1e-10
    assert abs(np.trace(dressed_prime)) < 1e-10
    assert norm(dressed_prime) > 1e-9, (norm(Xprime), norm(Wprime_X), norm(dressed_prime))

    # Independent symmetric difference for the derivative at one actual mode.
    p = alpha
    eps = 1e-5
    calc = []
    for t in (eps, -eps):
        l, v0 = np.linalg.eigh(p * gamma + h0_mass + t * h1)
        calc.append((v0 * transfer(l)) @ v0.conj().T)
    difference_error = norm((calc[0] - calc[1]) / (2 * eps) - mode(p)[1])
    assert difference_error < 1e-8
    rows.append({
        "quarter": alpha,
        "raw_endpoint_kernel_norm": norm(X),
        "kernel_first_derivative_norm": norm(Xprime),
        "isolated_wilson_contact_norm": norm(Wprime_X),
        "dressed_first_derivative_norm": norm(dressed_prime),
        "dressed_scalar_trace_absolute": float(abs(np.trace(dressed_prime))),
        "cutoff_change_kernel_norm": norm(X - X_short),
        "cutoff_change_derivative_norm": norm(Xprime - Xprime_short),
        "one_mode_finite_difference_error": difference_error,
    })
gap_bound = cutoffs[1] - float(np.linalg.norm(D, 2))
geometric_tail = np.exp(-s * gap_bound) / (np.pi * (1 - np.exp(-s)))
print(json.dumps({"arc_length": "pi/4", "collar_width": "1/8", "cutoffs": cutoffs,
                  "em_flip_exact_on_all_96_slots": True,
                  "tail_bound_kernel_norm": float(np.sqrt(96) * geometric_tail),
                  "tail_bound_first_derivative_norm": float(
                      max(s, 1 / gap_bound) * norm(h1) * geometric_tail),
                  "checks": rows}, indent=2))
