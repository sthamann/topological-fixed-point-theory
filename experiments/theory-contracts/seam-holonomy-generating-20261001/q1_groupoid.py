#!/usr/bin/env python3
"""Part J: Q1 on the level of the fundamental groupoid, not H_1.

Directed mark visits gamma_k (base point z = 0, the canonical rho-fixed point; out along the
ray to the mark i^k, once around it, back) are transported with the ORIGINAL v117 Fuchs
connection A(z) = sum_k U^k A0 U^-k / (z - i^k) by numerical path integration.  No monodromy
matrix is inserted; the exact v117 matrix M is used only for comparison at the end.

  J1 (Q1a) generators: every gamma_k has order-3 monodromy (exponents 0, 1/3, 2/3).
  J2 functor: transport of a concatenated path = product of transports; reversed path = inverse.
  J3 the quarter rotation acts on paths: T(i * gamma_k) = U T(gamma_k) U^-1 = T(gamma_k+1).
  J4 (Q1b) the surface relation: gamma_3 gamma_2 gamma_1 gamma_0 ~ big circle, trivial monodromy.
  J5 (Q1c) the orbifold relation (M U^-1)^2 = I is a PATH relation (loop around x = infinity
     of the quotient x = z^4) and reproduces diag M = (0, -i/2, i/2) of v117.
  J6 group images: four mark visits generate A4 (12), adding the rotation gives S4 (24).
  J7 (Q1d) non-abelian necessity: two words with the SAME H_1 class have NON-conjugate flavor action.
    PYTHONPATH=. python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/q1_groupoid.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from tfpt_explorer.flavor_path_transport import exact_clock_dictionary, original_fuchs_connection

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


def transport(segments):
    """segments: list of (z(t), dz(t)) on t in [0,1]; later segments act on the left."""
    T = np.eye(3, dtype=complex)
    for path in segments:
        def rhs(t, y, path=path):
            z, dz = path(t)
            return (dz * original_fuchs_connection(z) @ y.reshape(3, 3)).ravel()
        sol = solve_ivp(rhs, (0, 1), np.eye(3, dtype=complex).ravel(), method="DOP853", rtol=1e-12, atol=1e-14)
        T = sol.y[:, -1].reshape(3, 3) @ T
    return T


DELTA = 0.4
def ray(a, b):
    return lambda t: (a + (b - a) * t, b - a)
def circle(c, r, start_angle, turns=1):
    return lambda t: (c + r * np.exp(1j * (start_angle + 2 * np.pi * turns * t)),
                      2j * np.pi * turns * r * np.exp(1j * (start_angle + 2 * np.pi * turns * t)))
def gamma(k, sign=1):
    mk = 1j ** k
    p = (1 - DELTA) * mk
    segs = [ray(0, p), circle(mk, DELTA, np.angle(-mk), sign), ray(p, 0)]
    return segs
def reverse(segs):
    return [(lambda t, s=s: (s(1 - t)[0], -s(1 - t)[1])) for s in reversed(segs)]


M = [transport(gamma(k)) for k in range(4)]
w = np.exp(2j * np.pi / 3)
step("J1.(Q1a) each directed mark visit gamma_k has order-3 monodromy with eigenvalues {1, w, w^2}",
     all(np.allclose(np.linalg.matrix_power(Mk, 3), np.eye(3), atol=1e-8) for Mk in M)
     and all(np.allclose(sorted(np.angle(np.linalg.eigvals(Mk))), sorted(np.angle([1, w, w * w])), atol=1e-8) for Mk in M))
T01 = transport(gamma(0) + gamma(1))
Tinv = transport(reverse(gamma(0)))
step("J2.functor: T(gamma_0 then gamma_1) = T(gamma_1) T(gamma_0); T(reverse gamma_0) = T(gamma_0)^-1",
     np.allclose(T01, M[1] @ M[0], atol=1e-8) and np.allclose(Tinv, np.linalg.inv(M[0]), atol=1e-8))
U = np.diag([1, 1j, -1j])
step("J3.quarter rotation on paths: T(gamma_k+1) = U T(gamma_k) U^-1 (from i A(iz) = U A(z) U^-1)",
     all(np.allclose(M[(k + 1) % 4], U @ M[k] @ np.linalg.inv(U), atol=1e-8) for k in range(4)))
big = transport([circle(0, 2.0, 0.0)])
prod = M[3] @ M[2] @ M[1] @ M[0]
step("J4.(Q1b) surface relation: gamma_3 gamma_2 gamma_1 gamma_0 ~ big circle; both monodromies trivial (no log at infinity)",
     np.allclose(prod, big, atol=1e-7) and np.allclose(big, np.eye(3), atol=1e-7),
     f"|big - I| = {np.abs(big - np.eye(3)).max():.1e}")
Ui = np.linalg.inv(U)
rel = M[0] @ Ui
Mex = np.array(exact_clock_dictionary()["M"].evalf(), dtype=complex)
step("J5.(Q1c) orbifold path relation: (M_0 U^-1)^2 = I (loop around x = infinity of x = z^4, order 2); diag M_0 = (0, -i/2, i/2) as v117",
     np.allclose(rel @ rel, np.eye(3), atol=1e-8) and np.allclose(np.diag(M[0]), [0, -0.5j, 0.5j], atol=1e-8)
     and np.allclose(M[0] * M[0].T, Mex * Mex.T, atol=1e-8)
     and abs(M[0][0, 1] * M[0][1, 2] * M[0][2, 0] - Mex[0, 1] * Mex[1, 2] * Mex[2, 0]) < 1e-8,
     "equals the exact v117 M up to a diagonal similarity (raw frames are not orthonormal): M_ij M_ji and the 3-cycle product agree")
def closure(gens):
    key = lambda X: tuple(np.round(X, 6).ravel())
    seen = {key(np.eye(3)): np.eye(3, dtype=complex)}
    frontier = list(seen.values())
    while frontier:
        new = []
        for X in frontier:
            for G in gens:
                Y = G @ X
                if key(Y) not in seen:
                    seen[key(Y)] = Y
                    new.append(Y)
        frontier = new
        if len(seen) > 200:
            break
    return len(seen)
step("J6.images: the four mark visits generate A4 (order 12); adding the rotation U gives S4 (order 24)",
     closure(M) == 12 and closure(M + [U]) == 24, f"{closure(M)}, {closure(M + [U])}")
comm = M[1] @ M[0] @ np.linalg.inv(M[1]) @ np.linalg.inv(M[0])   # [gamma_0, gamma_1]: H_1 class 0
w1 = M[2] @ M[1] @ M[0]       # gamma_0 gamma_1 gamma_2
w2 = M[1] @ M[2] @ M[0]       # gamma_0 gamma_2 gamma_1  -- same H_1 class (1,1,1,0)
step("J7.(Q1d) non-abelian necessity: the commutator [g0,g1] has the H_1 class of the empty word but trace != 3; equal-class words g0g1g2 != g0g2g1",
     abs(np.trace(comm) - 3) > 1e-3 and not np.allclose(w1, w2, atol=1e-6),
     f"tr[g0,g1] = {np.round(np.trace(comm).real, 6)} (empty word: 3); |g0g1g2 - g0g2g1| = {np.abs(w1 - w2).max():.3f}")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG),
       "M0": np.round(M[0], 10).tolist().__repr__()}
(HERE / "q1_groupoid_results.json").write_text(json.dumps(out, indent=2))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
