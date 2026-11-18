"""Check endpoint assignments on the existing marked v252 quarter pair."""

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
three_q = {"u": 2, "d": -1, "e": -3, "nu": 0}
signs = [(-1)**(three_q[label] % 2) for label, *_ in v.P]
R = np.diag(signs + signs)
V = Q @ gamma
mass = -1j * gamma @ D
L = R @ J
S = R @ gamma @ J
s, arc = 1 / 8, np.pi / 4


def norm(a):
    return float(np.linalg.norm(a))


# Three exact comparisons of the polynomial coefficients of h_p(a).
for c in (gamma, V, mass):
    assert np.array_equal(L @ c @ L.T, -c)
    assert np.array_equal(S @ c.conj() @ S.T, -c)
for c, sign in ((gamma, 1), (V, -1), (mass, -1)):
    assert np.array_equal(R @ c.conj() @ R.T, sign * c)

# The geometric reference h_ref=p gamma is covered by the same checks.
assert np.array_equal(R @ Q.conj() @ R.T, -Q)
assert np.array_equal(S @ Q.conj() @ S.T, Q)
assert np.array_equal(L @ Q @ L.T, Q)

# Quarter conjugation at fixed spatial coordinates sends p to -p. On the
# existing pair this is m -> -m-1, so the original F0 mark is not fixed.
for m in range(-32, 32):
    assert (m % 8 == 0) == ((-m - 1) % 8 == 7)
    assert (-m - 1 + 0.75) == -(m + 0.25)

# A compensating one-mode shift restores the residue but changes h by gamma.
assert norm(gamma) == np.sqrt(96)
assert np.array_equal((-0.25 + 1) * gamma + mass - (-0.25 * gamma + mass), gamma)


def kernel_and_derivative(p):
    w, U = np.linalg.eigh(p * gamma + mass)
    e = np.exp(-s * np.abs(w))
    vals = e / (1 + e**2)
    dw = w[:, None] - w[None, :]
    close = np.abs(dw) < 1e-8
    dif = np.empty_like(dw)
    np.divide(vals[:, None] - vals[None, :], dw, out=dif, where=~close)
    mid = (w[:, None] + w[None, :]) / 2
    emid = np.exp(-s * np.abs(mid))
    dif[close] = (-s * np.tanh(s * mid) * emid / (1 + emid**2))[close]
    return ((U * vals) @ U.conj().T,
            U @ (dif * (U.conj().T @ V @ U)) @ U.conj().T)


def endpoints(alpha, direction):
    X = np.zeros_like(D)
    Xprime = np.zeros_like(D)
    # These paired cutoffs are exactly closed under m -> -m-1.
    for m in range(-384, 384):
        p = m + alpha
        block, derivative = kernel_and_derivative(p)
        phase = np.exp(-1j * p * arc * direction) / (2 * np.pi)
        X += phase * block
        Xprime += phase * derivative
    return X, Xprime - 1j * arc * direction * Q @ X


x, weak = endpoints(0.25, 1)
target, target_weak = endpoints(0.75, 1)
reflected, reflected_weak = endpoints(0.25, -1)
errors = {
    "paired_pointwise_antiunitary_cross_residual": norm(R @ x.conj() @ R.T - target),
    "paired_pointwise_antiunitary_weak_residual": norm(R @ weak.conj() @ R.T - target_weak),
    "fixed_quarter_reflected_antiunitary_cross_residual": norm(S @ x.conj() @ S.T - reflected),
    "fixed_quarter_reflected_antiunitary_weak_residual": norm(S @ weak.conj() @ S.T - reflected_weak),
    "linear_same_endpoints_cross_residual": norm(L @ x @ L.T - x),
    "linear_same_endpoints_weak_residual": norm(L @ weak @ L.T - weak),
}
assert max(errors.values()) < 1e-12, errors
print(json.dumps({
    "endpoints": [0, "pi/4"], "absolute_mode_bound": 384,
    "paired_pointwise_antiunitary_fixed_residue_zero_image": 7,
    "mode_repair_change_in_tangential_generator_norm": norm(gamma),
    "weak_derivative_norm": norm(weak),
    **errors,
}, indent=2))
