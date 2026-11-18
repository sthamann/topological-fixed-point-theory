#!/usr/bin/env python3
"""Part F: GRAV.FRAME.01 -- the frame field of the mark walk as vierbein, and what its vacuum does (not) induce.

  F1 frame = vierbein: replacing the axis Clifford operators by n_i.alpha with arbitrary unit vectors n_i keeps
     every controlled partial shift unitary; the IR dispersion is omega^2 = k^T G k with G_ij = n_i.n_j (the
     Gram matrix of the frame = inverse spatial metric), checked for random frames.
  F2 class restriction: G_ii = 1 for every such frame, hence det G <= 1 (Hadamard) with equality only for
     orthonormal frames: angles between the frame legs are free, a conformal factor or lapse is not.
  F3 sheet identity: every primitive step flips the sheet, so the sheet sign s anticommutes with the step
     operator and commutes with every on-site coin: s U s = -U for arbitrary on-site coins (position space,
     random coins).  Together with the +-omega pairing of the mass-coin walks the Floquet Dirac sea has
     quasi-energy exactly -pi per cell for every frame, uniform mass and sheet-odd mass.
  F4 CORRECTED (the first reading "no induced frame potential" was too strong): the same identity makes the full
     Floquet-sea energy exactly -pi per cell for EVERY background, including a magnetic field (Landau chain),
     although a continuum Dirac sea responds to B (vacuum polarisation).  Every state is paired with its
     temporal doubler (omega -> omega + pi); the full-sea measure is therefore blind and decides nothing about
     induced frame or gauge dynamics.  An IR-restricted measure (|omega| < pi/2) does respond to B.
  F5 (statement) Weinberg-Witten / spin 2: the IR walk is a free Dirac field with a Lorentz-covariant conserved
     stress tensor, so a massless spin-2 state cannot be a composite of the walk fields; the only graviton
     candidate here is the frame field itself (fundamental, the premise that is evaded).  Its dynamics is not
     induced at homogeneous order (F3/F4); the curvature response and the acceptance tests of
     GRAV.SPIN2.EMERGENCE.01 are not performed.  T7 stays open, no spin-2 marker moves.
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/frame_gravity.py
"""
from __future__ import annotations

import json
import sys
from itertools import product
from pathlib import Path

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


Zp = np.diag([1.0, -1.0]).astype(complex); I2 = np.eye(2, dtype=complex)
lower = np.array([[0, 1], [0, 0]], dtype=complex)
a0 = np.kron(lower, I2); a1 = np.kron(Zp, lower)
g = [a0 + a0.conj().T, -1j * (a0 - a0.conj().T), a1 + a1.conj().T, -1j * (a1 - a1.conj().T)]
w = np.array([[-.5, -.5, -.5], [-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]])
beta = sum(g) / 2
e = [sum(w[a, i] * g[a] for a in range(4)) for i in range(3)]
alpha = [1j * beta @ e[i] for i in range(3)]
I4 = np.eye(4)
Q = np.pi * np.ones(3)
rng = np.random.default_rng(23)


def unit_frame():
    v = rng.normal(size=(3, 3))
    return [x / np.linalg.norm(x) for x in v]


def Wn(k, n):
    out = I4.astype(complex)
    for i in range(3):
        a = sum(n[i][j] * alpha[j] for j in range(3))
        out = out @ (np.cos(k[i]) * I4 - 1j * np.sin(k[i]) * a)
    return out


frames = [unit_frame() for _ in range(3)] + [[np.array([1, 0, 0.]), np.array([np.cos(0.6), np.sin(0.6), 0.]), np.array([0, 0, 1.])]]
ks = rng.normal(size=(4, 3))
unit = all(np.allclose(Wn(k, n).conj().T @ Wn(k, n), I4) for n in frames for k in ks)
dirs = rng.normal(size=(12, 3)); dirs /= np.linalg.norm(dirs, axis=1)[:, None]
def ir_dev(n, kmag):
    G = np.array([[x @ y for y in n] for x in n])
    return max(abs(np.max(np.angle(np.linalg.eigvals(Wn(kmag * d, n)))) / np.sqrt(kmag ** 2 * d @ G @ d) - 1) for d in dirs)
devs = [(ir_dev(n, 2e-3), ir_dev(n, 1e-3)) for n in frames]
step("F1.frame = vierbein: with unit vectors n_i the controlled partial shifts stay unitary and the IR dispersion is omega^2 = k^T G k, G_ij = n_i.n_j (inverse spatial metric), for random frames",
     unit and all(d1 < 1e-2 and d1 / d2 > 1.6 for d1, d2 in devs), "max relative IR deviation " + ", ".join(f"{d1:.1e}->{d2:.1e}" for d1, d2 in devs))

Gs = [np.array([[x @ y for y in n] for x in n]) for n in frames]
diag_one = all(np.allclose(np.diag(G), 1) for G in Gs)
hadamard = all(np.linalg.det(G) <= 1 + 1e-12 for G in Gs)
step("F2.class restriction: G_ii = 1 for every unit frame, so det G <= 1 (Hadamard), equality only for orthonormal frames: angles are free, a conformal factor or lapse is not representable in this class",
     diag_one and hadamard, "det G = " + ", ".join(f"{np.linalg.det(G):.3f}" for G in Gs))

# F3: position-space sheet identity with random on-site coins -----------------------------------------
Lb = 4
sites = [x for x in product(range(Lb), repeat=3) if len({c % 2 for c in x}) == 1]
idx = {x: i for i, x in enumerate(sites)}
signs = [np.array(s) for s in product((-1, 1), repeat=3)]
A = {tuple(s): (I4 + s[0] * alpha[0]) / 2 @ (I4 + s[1] * alpha[1]) / 2 @ (I4 + s[2] * alpha[2]) / 2 for s in signs}
dim = 4 * len(sites)
Wr = np.zeros((dim, dim), complex)
for x, i in idx.items():
    for s_, As in A.items():
        y = tuple((x[c] + s_[c]) % Lb for c in range(3))
        Wr[4 * idx[y]:4 * idx[y] + 4, 4 * i:4 * i + 4] += As
def rand_u4():
    q, r = np.linalg.qr(rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)))
    return q
C1 = np.zeros_like(Wr); C2 = np.zeros_like(Wr); Sg = np.zeros_like(Wr)
for x, i in idx.items():
    C1[4 * i:4 * i + 4, 4 * i:4 * i + 4] = rand_u4(); C2[4 * i:4 * i + 4, 4 * i:4 * i + 4] = rand_u4()
    Sg[4 * i:4 * i + 4, 4 * i:4 * i + 4] = (1 if x[0] % 2 == 0 else -1) * I4
Ur = C1 @ Wr @ C2
sheet_id = np.allclose(Sg @ Ur @ Sg, -Ur) and np.allclose(Sg @ Wr @ Sg, -Wr)

sx = np.array([[0, 1], [1, 0]])
def U8(k, n, m, ms):
    Wk = np.zeros((8, 8), complex); Wk[:4, :4] = Wn(k, n); Wk[4:, 4:] = Wn(k + Q, n)
    Cc = expm(-0.5j * (m * np.kron(np.eye(2), beta) + ms * np.kron(sx, beta)))
    return Cc @ Wk @ Cc
Ng = 10
kg = 2 * np.pi * (np.arange(Ng) + 0.31) / Ng
grid = list(product(kg, repeat=3))
def evac(n, m, ms):
    tot, pair = 0.0, True
    for k in grid:
        q = np.sort(np.angle(np.linalg.eigvals(U8(np.array(k), n, m, ms))))
        tot += np.sum(q[q < 0])
        pair &= np.allclose(np.sort(q), np.sort(-q), atol=1e-9)
    return tot / len(grid) / 2, pair
E_tab, pair_all = [], True
for n in frames[:3]:
    for m, ms in ((0.0, 0.0), (0.3, 0.0), (0.2, 0.5)):
        Ev, pr = evac(n, m, ms)
        E_tab.append(Ev); pair_all &= pr
step("F3.sheet identity: every step flips the sheet, so s U s = -U for ARBITRARY on-site coins (position space, random coins); with the +-omega pairing of the mass-coin walks the Floquet Dirac sea has quasi-energy exactly -pi per cell for every frame, uniform mass and sheet-odd mass",
     sheet_id and pair_all and np.allclose(E_tab, -np.pi, atol=1e-10), f"E_vac over 3 frames x 3 coin choices: max |E + pi| = {max(abs(x + np.pi) for x in E_tab):.1e}")
Pp_ = [(I4 + alpha[i]) / 2 for i in range(3)]; Pm_ = [(I4 - alpha[i]) / 2 for i in range(3)]
def landau(Lx, p, ky, kz):
    Bf = 2 * np.pi * p / Lx
    V1 = np.zeros((4 * Lx, 4 * Lx), complex); V2 = np.zeros_like(V1); V3 = np.zeros_like(V1)
    for x in range(Lx):
        xp = (x + 1) % Lx
        V1[4 * xp:4 * xp + 4, 4 * x:4 * x + 4] += Pp_[0]; V1[4 * x:4 * x + 4, 4 * xp:4 * xp + 4] += Pm_[0]
        k2 = ky - Bf * x
        V2[4 * x:4 * x + 4, 4 * x:4 * x + 4] = np.exp(-1j * k2) * Pp_[1] + np.exp(1j * k2) * Pm_[1]
        V3[4 * x:4 * x + 4, 4 * x:4 * x + 4] = np.exp(-1j * kz) * Pp_[2] + np.exp(1j * kz) * Pm_[2]
    return V1 @ V2 @ V3
full_B, ir_B = [], []
for p_ in (3, 6, 12):
    q = np.angle(np.linalg.eigvals(landau(240, p_, 0.13, 0.21)))
    full_B.append(np.sum(q[q < 0]) / 240); ir_B.append(np.sum(q[(q < 0) & (q > -np.pi / 2)]) / 240)
step("F4.CORRECTED: the full Floquet-sea energy is -pi per cell for every frame, coin AND magnetic field (Landau chain, three fluxes) although a continuum Dirac sea responds to B: the measure is blind (temporal-doubler pairing) and does not decide induced frame or gauge dynamics; the IR-restricted measure (|omega| < pi/2) does respond to B",
     np.ptp(E_tab) < 1e-10 and np.allclose(full_B, -np.pi, atol=1e-10) and np.ptp(ir_B) > 1e-3,
     f"full sea per site at B = 2 pi p/240, p = 3, 6, 12: {[round(x, 10) for x in full_B]}; IR-restricted: {[round(x, 5) for x in ir_B]}")
step("F5.(statement) Weinberg-Witten / spin 2: the IR walk is a free Dirac field with covariant conserved T_munu, so massless spin 2 cannot be a composite of walk fields; the graviton candidate is the frame field itself (fundamental: the evaded premise); whether the walk induces its dynamics is undecided (F4); curvature response and GRAV.SPIN2.EMERGENCE.01 acceptance tests not performed; T7 open",
     True, "no spin-2 or TOE marker moves")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "E_vac": E_tab}
(HERE / "frame_gravity_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
