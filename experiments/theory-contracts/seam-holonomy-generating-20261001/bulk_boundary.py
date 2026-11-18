#!/usr/bin/env python3
"""Part X: TFPT.BULK.BOUNDARY.RECONSTRUCTION -- is the seam Weyl function a LOCAL bulk problem?

Bulk = the Euclidean disc whose reflected double is the chi = 2 sphere of P1 (radial quantisation:
the seam circle is a time slice, the disc path integral prepares the vacuum, disc + reflected disc
= S^2 is the reflection-positive norm).
  X1 (Q1) M(z) is the ROBIN-TO-DIRICHLET map of the LOCAL disc Laplacian: Delta u = 0, d_r u - z u = f
     on the seam, u|_seam = M(z) f.  Computed with a sparse polar finite-difference Laplacian, a local
     point source and a local trace (no Fourier ansatz): K(-2) - K(-1) at the square marks converges
     to the geometric values (-1/2pi, 0.019473, 0.036194, 0.019473).
  X2 (Q1 uniqueness, Calderon) a massive bulk (Delta - m^2) or a radial conductivity changes the
     Dirichlet-to-Neumann spectrum away from |n|; in 2D the DtN map determines an isotropic
     conductivity (Astala-Paivarinta, cited) -- so DtN = |D| selects the conformal (massless,
     constant-conductivity) disc; any conformal deformation keeps it.
  X3 (Q2) the bulk Dirichlet energy of the harmonic extension is the H^{1/2} form 2 pi sum |n||f_n|^2:
     the canonical 1/sqrt|n| access of Part W is the trace map H^1(disc) -> H^{1/2}(S^1).
  X4 (Q4) radial quantisation: DtN = r d_r at r = 1 on harmonic functions = the dilation generator
     L0 + Lbar0 built from the bulk stress tensor; M(z) is the resolvent of that generator.
  X5 (limit) the bulk is a 2D Euclidean disc (compiler geometry); it does not produce 3+1D (Theorem D).
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/bulk_boundary.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla
import sympy as sp
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


def disc_operator(Nr, Nt, z, sigma=lambda r: 1.0, mass=0.0):
    """Sparse polar FD for div(sigma grad u) - mass^2 u = 0 with Robin d_r u - z u = f at r = 1."""
    h, dt = 1.0 / Nr, 2 * np.pi / Nt
    idx = lambda i, k: 1 + (i - 1) * Nt + (k % Nt)            # 0 = centre
    n = 1 + Nr * Nt
    rows, cols, vals = [], [], []
    def add(a, b, v):
        rows.append(a); cols.append(b); vals.append(v)
    add(0, 0, -4 / h ** 2 - mass ** 2)
    for k in range(Nt):
        add(0, idx(1, k), 4 / (h ** 2 * Nt))
    for i in range(1, Nr + 1):
        r = i * h
        sp_, sm_ = sigma(r + h / 2), sigma(r - h / 2)
        for k in range(Nt):
            a = idx(i, k)
            cr_out = (r + h / 2) * sp_ / (r * h ** 2)
            cr_in = (r - h / 2) * sm_ / (r * h ** 2)
            ct = sigma(r) / (r ** 2 * dt ** 2)
            add(a, a, -(cr_out + cr_in) - 2 * ct - mass ** 2)
            add(a, idx(i, k + 1), ct); add(a, idx(i, k - 1), ct)
            add(a, idx(i - 1, k) if i > 1 else 0, cr_in)
            if i < Nr:
                add(a, idx(i + 1, k), cr_out)
            else:  # ghost: u_{N+1} = u_{N-1} + 2h (f + z u_N)  (sigma(1) = 1 at the seam)
                add(a, idx(i - 1, k), cr_out)
                add(a, a, cr_out * 2 * h * z)
    A = sps.csr_matrix((vals, (rows, cols)), shape=(n, n))
    return A, idx, h


def robin_response(Nr, Nt, z, marks_k, src_k, **kw):
    A, idx, h = disc_operator(Nr, Nt, z, **kw)
    dt = 2 * np.pi / Nt
    f = np.full(Nt, -1 / (2 * np.pi)); f[src_k] += 1 / dt     # mean-zero point source (density per angle)
    b = np.zeros(A.shape[0])
    r = 1.0
    cr_out = (r + h / 2) / (r * h ** 2)
    for k in range(Nt):
        b[idx(Nr, k)] = -cr_out * 2 * h * f[k]
    u = spla.spsolve(A.tocsc(), b)
    return np.array([u[idx(Nr, k)] for k in marks_k])


geo_row = [sum(np.cos(n * th) * (1 / (n + 2) - 1 / (n + 1)) for n in range(1, 200000)) / np.pi
           for th in (0.0, np.pi / 2, np.pi, 3 * np.pi / 2)]
devs, rows_ = [], {}
for Nr, Nt in ((32, 128), (64, 256), (128, 512)):
    marks = [k * Nt // 4 for k in range(4)]
    row = robin_response(Nr, Nt, -2.0, marks, 0) - robin_response(Nr, Nt, -1.0, marks, 0)
    rows_[(Nr, Nt)] = row
    devs.append(np.abs(row - geo_row).max())
step("X1.(Q1) M(z) = Robin-to-Dirichlet map of the LOCAL sparse disc Laplacian: K(-2)-K(-1) at the square marks -> geometric row, converging",
     devs[2] < devs[1] < devs[0] and devs[2] < 2e-3,
     "max dev " + ", ".join(f"{d:.1e}" for d in devs) + "; finest " + str(np.round(rows_[(128, 512)], 5).tolist())
     + " vs " + str(np.round(geo_row, 5).tolist()))


def dtn_eigen(n, mass=0.0, sigma=lambda r: 1.0, dsigma=lambda r: 0.0):
    """DtN eigenvalue for mode n: regular radial solution of (sigma u')' + sigma u'/r - (sigma n^2/r^2 + m^2) u = 0."""
    r0 = 1e-4
    u0, du0 = r0 ** n, n * r0 ** (n - 1)
    def rhs(r, y):
        u, du = y
        s, ds = sigma(r), dsigma(r)
        return [du, -(ds / s + 1 / r) * du + (n ** 2 / r ** 2 + mass ** 2 / s) * u]
    sol = solve_ivp(rhs, (r0, 1.0), [u0, du0], rtol=1e-11, atol=1e-14)
    u, du = sol.y[:, -1]
    return du / u                                               # sigma(1) = 1 at the seam
flat = [dtn_eigen(n) for n in (1, 2, 3, 4)]
massive = [dtn_eigen(n, mass=1.0) for n in (1, 2, 3, 4)]
cond = [dtn_eigen(n, sigma=lambda r: 1 + r * r * (1 - r) ** 2 * 4, dsigma=lambda r: 4 * (2 * r * (1 - r) ** 2 - 2 * r * r * (1 - r)))
        for n in (1, 2, 3, 4)]
step("X2.(uniqueness) flat disc DtN = |n|; massive bulk and a radial conductivity (sigma = 1 at the seam) give DtN != |n| -> different M(z) (Calderon: DtN fixes an isotropic conductivity)",
     np.allclose(flat, [1, 2, 3, 4], atol=1e-6) and min(abs(a - b) for a, b in zip(massive, flat)) > 1e-2
     and min(abs(a - b) for a, b in zip(cond, flat)) > 1e-3,
     f"flat {np.round(flat, 5).tolist()}, massive {np.round(massive, 4).tolist()}, conductivity {np.round(cond, 4).tolist()}")

rr, th = sp.symbols("r theta", positive=True)
nn = sp.symbols("n", positive=True, integer=True)
u = rr ** nn * sp.cos(nn * th)
energy = sp.integrate(sp.integrate(sp.simplify((sp.diff(u, rr) ** 2 + sp.diff(u, th) ** 2 / rr ** 2) * rr), (th, 0, 2 * sp.pi)), (rr, 0, 1))
step("X3.(Q2) Dirichlet energy of the harmonic extension of cos(n theta) = pi n = 2 pi |n| |f_n|^2: the H^{1/2} trace form, i.e. the canonical 1/sqrt|n| access",
     sp.simplify(energy - sp.pi * nn) == 0)
dil = sp.simplify(rr * sp.diff(u, rr) / u)
step("X4.(Q4) on harmonic functions DtN = r d_r = dilation generator (radial quantisation: L0 + Lbar0 from the bulk stress tensor); M(z) is its resolvent",
     dil == nn)
step("X5.(statement) the reconstructing bulk is the 2D Euclidean disc (compiler geometry; its reflected double is the chi = 2 sphere of P1): no 3+1D geometry follows",
     True, "Theorem D stays open; the disc is the state-preparing Euclidean half of the RP double")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG),
       "X1_devs": devs, "dtn": {"flat": flat, "massive": massive, "conductivity": cond}}
(HERE / "bulk_boundary_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
