#!/usr/bin/env python3
"""Part K: SEAM.ACCESS.01 -- are the four raw stations the four local accesses of ONE seam?

Instead of postulating "marks = stations", ask which implementations of the raw clock are
possible at all.  Ring vacuum C on N sites (the reconstructed NS kernel).
  K1 the vacuum's nearest-neighbour structure is the N-cycle: |C| is maximal exactly at
     distance 1, so every site permutation preserving |C| lies in Aut(C_N) = D_N.
  K2 elements of D_N of order 4 whose orbits all have size 4: only the rotations by +-N/4.
  K3 sign lifts: a rotation by N/4 with site signs preserving C exists only as +- the twisted
     quarter translation (signs fixed up to a global sign); its fourth power is -I.
  => any vacuum-preserving, locality-preserving (signed site permutation) implementation of a
     fixed-point-free order-4 clock is the quarter rotation, and every access orbit is a square.
  K4 control 'four copies': four disconnected M-site rings with the same block clock: U J = J Q holds,
     but all mixed responses J_i^dag C J_j vanish and the state is not homogeneous (violates K0).
  K5 control 'delocalised stations': raw labels spread over the circle (raw-side sign gauge): the
     clock is still the quarter rotation on sites, mixed responses are non-zero, but stations are
     not arcs.  So U J = J Q plus mixed responses does NOT make stations arcs -- and lambda = -1
     does not need it: it needs only K1-K3.
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/seam_access.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


def ns_kernel(N):
    d = np.subtract.outer(np.arange(N), np.arange(N))
    with np.errstate(divide="ignore", invalid="ignore"):
        C = np.sin(np.pi * d / 2) / (N * np.sin(np.pi * d / N))
    C[d == 0] = 0.5
    return C


N, M = 16, 4
C = ns_kernel(N)
A = np.abs(C - 0.5 * np.eye(N))
mx = A.max()
nn = {(i, j) for i in range(N) for j in range(N) if abs(A[i, j] - mx) < 1e-12}
step("K1.|C| is maximal exactly at circular distance 1: the vacuum's nearest-neighbour graph is the N-cycle (Aut = D_N)",
     nn == {(i, (i + s) % N) for i in range(N) for s in (1, -1)})
dihedral = [tuple((s * j + t) % N for j in range(N)) for t in range(N) for s in (1, -1)]
pres = [p for p in dihedral if np.allclose(A[np.ix_(p, p)], A)]
step("K1b.all 2N dihedral maps preserve |C|", len(pres) == 2 * N)


def orbit_sizes(p):
    seen, sizes = set(), []
    for j in range(N):
        if j in seen:
            continue
        k, n = j, 0
        while k not in seen:
            seen.add(k); k = p[k]; n += 1
        sizes.append(n)
    return sizes


ord4 = [p for p in dihedral if set(orbit_sizes(p)) == {4}]
rot = [tuple((j + N // 4) % N for j in range(N)), tuple((j - N // 4) % N for j in range(N))]
step("K2.order-4 elements of D_N with all orbits of size 4 are exactly the rotations by +-N/4", sorted(ord4) == sorted(rot))


def sign_lifts(p):
    """All s in {+-1}^N with diag(s) P C P^T diag(s) = C, P the permutation p."""
    P = np.zeros((N, N)); P[list(p), list(range(N))] = 1
    Cp = P @ C @ P.T
    sols = []
    for s0 in (1, -1):
        s = np.zeros(N); s[0] = s0
        ok = True
        for _ in range(N):
            for i in range(N):
                for j in range(N):
                    if abs(C[i, j]) > 1e-12 and i != j and s[i] != 0 and s[j] == 0:
                        s[j] = s[i] * np.sign(Cp[i, j]) * np.sign(C[i, j])
        if (s == 0).any():
            ok = False
        if ok and np.allclose(np.outer(s, s) * Cp, C):
            sols.append(s)
    return sols, P


lifts, P = sign_lifts(rot[0])
Tq = np.diag(lifts[0]) @ P if lifts else None
Ttw = np.zeros((N, N))
for j in range(N):
    Ttw[(j + M) % N, j] = -1.0 if j + M >= N else 1.0
step("K3.sign lift of the quarter rotation preserving C: exactly two (global sign), = +-twisted quarter translation, T^4 = -I",
     len(lifts) == 2 and (np.allclose(Tq, Ttw) or np.allclose(Tq, -Ttw))
     and np.allclose(np.linalg.matrix_power(Tq, 4), -np.eye(N)))
z = np.exp(2j * np.pi * np.arange(N) / N)
crs = [((z[p] - z[(p + 2 * M) % N]) * (z[(p + M) % N] - z[(p + 3 * M) % N])
        / ((z[p] - z[(p + 3 * M) % N]) * (z[(p + M) % N] - z[(p + 2 * M) % N]))) for p in range(N)]
step("K3b.every access orbit of this clock is a square (cross-ratio 2): lambda = -1 without 'marks = stations'",
     np.allclose(crs, 2))

# raw block clock and access maps
Rst = np.zeros((4, 4)); Rst[1, 0] = Rst[2, 1] = Rst[3, 2] = 1; Rst[0, 3] = -1
Rraw = np.kron(Rst, np.eye(M))
J = [np.kron(np.eye(4)[:, [a]], np.eye(M)) for a in range(4)]          # station a -> block a
Q = Rst
UJ_JQ = all(np.allclose(Rraw @ J[a], sum(Q[b, a] * J[b] for b in range(4))) for a in range(4))
mixed_ring = [np.abs(J[i].T @ C @ J[j]).max() for i in range(4) for j in range(4) if i != j]
Cc = np.kron(np.eye(4), ns_kernel(M))                                   # four disconnected copies
mixed_copy = [np.abs(J[i].T @ Cc @ J[j]).max() for i in range(4) for j in range(4) if i != j]
T1 = np.zeros((N, N))
for j in range(N):
    T1[(j + 1) % N, j] = -1.0 if j + 1 >= N else 1.0
step("K4.control four copies: U J = J Q and R-invariance hold, but every mixed response vanishes and K0 (homogeneity) fails",
     UJ_JQ and np.allclose(Rraw @ Cc @ Rraw.T, Cc) and max(mixed_copy) < 1e-14 and min(mixed_ring) > 1e-3
     and not np.allclose(T1 @ Cc @ T1.T, Cc),
     f"min mixed response ring {min(mixed_ring):.3f}, copies {max(mixed_copy):.1e}")
# delocalised stations: raw label (a,l) -> site l + M*((a + c_l) mod 4) with c_l not constant
c = [0, 1, 0, 1]
site = {(a, l): l + M * ((a + c[l]) % 4) for a in range(4) for l in range(M)}
Pd = np.zeros((N, N))
for (a, l), sidx in site.items():
    Pd[sidx, a * M + l] = 1
# raw-side sign gauge: move the single -1 of each 4-orbit to where the site orbit crosses N
Rsite = Pd @ Rraw @ Pd.T
# find a site-diagonal sign D with D Rsite D = +-Ttw (always possible: one -1 per 4-orbit)
D = np.ones(N)
for l in range(M):
    orbit = [l + M * k for k in range(4)]
    for k in range(3):
        j = orbit[k]; nxt = orbit[(k + 1) % 4]
        want = Ttw[nxt, j]; have = Rsite[nxt, j]
        D[nxt] = D[j] * want / have
Rg = np.diag(D) @ Rsite @ np.diag(D)
stations_arcs = all(sorted(site[(a, l)] for l in range(M)) == list(range(M * a, M * a + M)) for a in range(4))
mixed_deloc = min(np.abs((Pd @ J[i]).T @ C @ (Pd @ J[j])).max() for i in range(4) for j in range(4) if i != j)
step("K5.control delocalised stations: after raw-side sign gauge the clock IS the quarter rotation, mixed responses non-zero, yet stations are not arcs",
     np.allclose(Rg, Ttw) and not stations_arcs and mixed_deloc > 1e-3 and np.allclose(Rg @ C @ Rg.T, C),
     "U J = J Q + mixed responses do not localise stations; lambda = -1 needs only K1-K3")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG)}
(HERE / "seam_access_results.json").write_text(json.dumps(out, indent=2))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
