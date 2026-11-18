#!/usr/bin/env python3
"""Part H: CHIRAL.FLOQUET.01 -- can the free mark walk have a chiral IR at quasi-energy 0?

On one chirality block (Part K2: alpha_i = -Gamma S_i, a Pauli triple) the walk is W(k) = prod_i (cos k_i - i sin k_i
sigma_i).  Two free-level escape routes from the doubler of Part L5 are tested:
  H1 Floquet separation (positive part): adding one uncontrolled primitive step, U = T_eps W (eps a bcc step,
     i.e. one extra mark visit), moves the opposite-chirality nodes to quasi-energy -pi/2: all quasi-energy-0
     nodes then have the SAME chirality (net -2 per bcc zone).
  H2 ... but it is not a Weyl semimetal: (a) the cones are overtilted (tilt |eps| = sqrt3 > speed 1, type II);
     (b) exactly: det U(k) = e^{-2i n.k}, so along every axis with n_i != 0 the two bands wind and MUST cross
     quasi-energy 0 at least 2|n_i| times on every line -- a Floquet metal with quasi-energy-0 surfaces.
     Every nonzero tilt n (a lattice translation) gives such a winding.
  H3 without tilt: every circuit built from controlled single-axis shifts and on-site SU(2) coins is a pointwise
     product of maps T^3 -> SU(2) that factor through a circle, hence null-homotopic (pi_1(S^3) = 0): its degree
     is 0 and the net chirality at quasi-energy 0 (and at pi) vanishes.  Checked with the winding integral on the
     BDT walk and random circuits, validated on a (non-local) degree-1 control map.
  H4 (statement) free-level chirality no-go for the walk architecture of Parts K/L/N (controlled shifts, coins,
     translations): a chiral IR at quasi-energy 0 without quasi-energy Fermi surfaces is impossible; the chiral
     Standard Model needs interactions (a symmetric gap of the doubler sector, CHIRAL4D.NOMIRROR.01).
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/chiral_floquet.py
"""
from __future__ import annotations

import json
import sys
from itertools import product
from pathlib import Path

import numpy as np
from scipy.linalg import expm, logm
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1.0 + 0j, -1]); I2 = np.eye(2)
P = [sx, sy, sz]
Uf = lambda t, A: np.cos(t) * I2 - 1j * np.sin(t) * A
W = lambda k: Uf(k[0], sx) @ Uf(k[1], sy) @ Uf(k[2], sz)


def gap(U):
    ev = np.linalg.eigvals(U)
    return abs(ev[0] - ev[1])


def census(f, N=24):
    ks = 2 * np.pi * np.arange(N) / N - np.pi
    G = np.array([[[gap(f(np.array([a, b, c]))) for c in ks] for b in ks] for a in ks])
    pts = []
    for i, j, l in product(range(N), repeat=3):
        v = G[i, j, l]
        nb = [G[(i + a) % N, (j + b) % N, (l + c) % N] for a, b, c in product((-1, 0, 1), repeat=3) if (a, b, c) != (0, 0, 0)]
        if v < 0.8 and v <= min(nb):
            r = minimize(lambda k: gap(f(k)), np.array([ks[i], ks[j], ks[l]]), method="Nelder-Mead", options={"xatol": 1e-11, "fatol": 1e-14, "maxiter": 5000})
            if r.fun < 1e-6:
                pts.append(np.mod(r.x + np.pi, 2 * np.pi) - np.pi)
    uniq = []
    for p in pts:
        if not any(np.allclose(np.mod(p - q + np.pi, 2 * np.pi) - np.pi, 0, atol=1e-5) for q in uniq):
            uniq.append(p)
    out = []
    for k0 in uniq:
        U0 = f(k0)
        lam = np.linalg.eigvals(U0)[0]
        M = np.zeros((3, 3))
        for jj in range(3):
            dk = np.zeros(3); dk[jj] = 1e-5
            H = 1j * logm(f(k0 + dk) / lam)
            M[:, jj] = [np.real(np.trace(H @ P[a])) / 2 / 1e-5 for a in range(3)]
        out.append((k0, np.angle(lam), np.sign(np.linalg.det(M)), M))
    return out


eps = np.array([1.0, 1.0, 1.0])
Ut = lambda k: np.exp(-1j * eps @ k) * W(k)
nodes = census(Ut)
z = [n for n in nodes if abs(n[1]) < 1e-6]
others = [n for n in nodes if abs(n[1]) >= 1e-6]
net0 = sum(n[2] for n in z)
phases_other = sorted({round(n[1], 4) for n in others})
step("H1.Floquet separation: U = T_eps W (one extra uncontrolled primitive step) puts all quasi-energy-0 nodes at the SAME chirality (net -2 per bcc zone), the opposite charges at quasi-energy -pi/2",
     len(z) == 8 and abs(net0) == 8 and all(n[2] == z[0][2] for n in z) and phases_other == [round(-np.pi / 2, 4)],
     f"{len(z)} qe-0 nodes in the cubic zone (= 2 per bcc zone), net {net0:+.0f}; other nodes at {phases_other}")

tilt = max(np.linalg.norm(np.linalg.solve(n[3].T, eps)) for n in z)
rng = np.random.default_rng(31)
cross_counts = []
for _ in range(12):
    kx, ky = rng.uniform(-np.pi, np.pi, 2)
    kz = np.linspace(-np.pi, np.pi, 2001)
    ph = np.array([np.sort(np.angle(np.linalg.eigvals(Ut(np.array([kx, ky, t]))))) for t in kz])
    c = 0
    for b in range(2):
        s = np.sign(ph[:, b])
        jumps = np.where(np.diff(s) != 0)[0]
        c += sum(1 for j in jumps if abs(ph[j, b]) < 0.5 and abs(ph[j + 1, b]) < 0.5)       # genuine zero crossings, not +-pi wraps
    cross_counts.append(c)
detw = all(np.isclose(np.linalg.det(Ut(k)), np.exp(-2j * eps @ k)) for k in rng.normal(size=(5, 3)))
step("H2.not a Weyl semimetal: the qe-0 cones are overtilted (|M^-T eps| > 1, type II) and det U = e^{-2i n.k} winds, so on EVERY k_z line the bands cross quasi-energy 0 at least 2|n_z| = 2 times: quasi-energy-0 surfaces (Floquet metal)",
     tilt > 1 and detw and min(cross_counts) >= 2, f"tilt ratio {tilt:.3f}; zero crossings per random k_z line {cross_counts}")


def rand_coin(theta=0.35):
    v = rng.normal(size=3); v /= np.linalg.norm(v)
    return expm(-1j * theta * (v[0] * sx + v[1] * sy + v[2] * sz))


Ng = 48
kk = 2 * np.pi * (np.arange(Ng) + 0.5) / Ng - np.pi
KX, KY, KZ = [x.ravel() for x in np.meshgrid(kk, kk, kk, indexing="ij")]
fac = lambda t, A: np.cos(t)[:, None, None] * I2 - 1j * np.sin(t)[:, None, None] * A
dfac = lambda t, A: -np.sin(t)[:, None, None] * I2 - 1j * np.cos(t)[:, None, None] * A


def degree(Us, dUs):
    """(1/24 pi^2) int eps^{ijk} tr(U^-1 d_i U U^-1 d_j U U^-1 d_k U) d^3k on the grid (exact derivatives)."""
    Ui = np.linalg.inv(Us)
    Ad = [Ui @ d for d in dUs]
    tot = 0
    for (i, j, k), sg in (((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((0, 2, 1), -1), ((2, 1, 0), -1), ((1, 0, 2), -1)):
        tot += sg * np.trace(Ad[i] @ Ad[j] @ Ad[k], axis1=1, axis2=2).sum()
    return (tot * (2 * np.pi / Ng) ** 3 / (24 * np.pi ** 2)).real


def circuit_degree(C):
    X, Y, Zf = fac(KX, sx), fac(KY, sy), fac(KZ, sz)
    dX, dY, dZ = dfac(KX, sx), dfac(KY, sy), dfac(KZ, sz)
    U = C[0] @ X @ C[1] @ Y @ C[2] @ Zf @ C[3]
    dU = [C[0] @ dX @ C[1] @ Y @ C[2] @ Zf @ C[3], C[0] @ X @ C[1] @ dY @ C[2] @ Zf @ C[3], C[0] @ X @ C[1] @ Y @ C[2] @ dZ @ C[3]]
    return degree(U, dU)


def control_degree(m=2.0):
    """Degree-1 control (not strictly local): U = (n0 - i n.sigma)/|n|, n = (m - sum cos k, sin k)."""
    h = 1e-5
    def U_of(kx_, ky_, kz_):
        n0 = m - (np.cos(kx_) + np.cos(ky_) + np.cos(kz_))
        n = np.stack([np.sin(kx_), np.sin(ky_), np.sin(kz_)], axis=1)
        nn = np.sqrt(n0 ** 2 + (n ** 2).sum(1))
        return (n0[:, None, None] * I2 - 1j * (n[:, 0, None, None] * sx + n[:, 1, None, None] * sy + n[:, 2, None, None] * sz)) / nn[:, None, None]
    U = U_of(KX, KY, KZ)
    dU = [(U_of(KX + h, KY, KZ) - U_of(KX - h, KY, KZ)) / (2 * h), (U_of(KX, KY + h, KZ) - U_of(KX, KY - h, KZ)) / (2 * h), (U_of(KX, KY, KZ + h) - U_of(KX, KY, KZ - h)) / (2 * h)]
    return degree(U, dU)


ctrl = control_degree()
degs = [circuit_degree([rand_coin() for _ in range(4)]) for _ in range(4)] + [circuit_degree([rand_coin(1.2) for _ in range(4)]) for _ in range(2)]
bdt = circuit_degree([I2, I2, I2, I2])
step("H3.without tilt: every circuit of controlled single-axis shifts and on-site SU(2) coins is null-homotopic (each factor factors through a circle, pi_1(S^3) = 0), hence degree 0 and net chirality 0 at quasi-energy 0 and at pi; winding integral: BDT walk and six random circuits give 0, a (non-local) degree-1 control gives 1",
     abs(abs(ctrl) - 1) < 0.02 and abs(bdt) < 1e-6 and max(abs(d) for d in degs) < 1e-6,
     f"control {ctrl:.4f}; BDT {bdt:.1e}; random circuits {[f'{d:.1e}' for d in degs]}")
step("H4.(statement) free-level chirality no-go for the walk architecture of Parts K/L/N: a chiral IR at quasi-energy 0 without quasi-energy Fermi surfaces is impossible (H2, H3, L5); the chiral Standard Model needs interactions -- a symmetric gap of the doubler sector (CHIRAL4D.NOMIRROR.01)",
     True, "no chirality or TOE marker moves")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "tilt_ratio": tilt, "crossings": cross_counts, "random_circuit_degrees": degs, "control_degree": ctrl}
(HERE / "chiral_floquet_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
