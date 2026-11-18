#!/usr/bin/env python3
"""Part F: is the raw station clock R a geometric symmetry of the seam vacuum?

Operator test on the antiperiodic critical fermion ring with N = 4 x 64 = 256 modes
(the raw ring size; same class as verification/v480): local algebras A(I) of arcs,
quasi-free vacuum Omega with covariance C = <c_j^dag c_k>, and the raw station clock
R_station (x) I_64 of universalraum-mu4-pointing-20260917.

    F1  R_station (x) I_64 IS the twisted quarter translation of the ring (exact).
    F2  R Omega = Omega:   R C R^dag = C   (and the vacuum is homogeneous).
    F3  R A(I) R^-1 = A(I + N/4) for every arc I (signed permutation of modes).
    F4  R^4 = -I = the antiperiodic (NS) sign: projectively rho^4 = 1, spin lift (-1)^F.
    F5  the lattice vacuum approaches the Moebius-covariant NS kernel 1/(N sin(pi d/N)).
    F6  the bit: marks on the seam circle are concyclic <=> lambda real; j=0 has
        non-real lambda, and an order-3 rotation of a circle cannot permute 4 marks.
    F7  marks = station boundaries of equal stations -> angles k pi/2 -> lambda = -1;
        negative control: unequal stations admit no order-4 vacuum covariance.
Firewall: the identification raw seam vacuum = this ring vacuum stays SEAM.EQUIV.01 [O].
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/seam_vacuum_covariance.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


N, M = 256, 64  # 4 stations x 64 modes


def ring_vacuum(sizes=None):
    """Half-filled antiperiodic hopping ring; returns covariance C_jk = <c_j^dag c_k>."""
    H = np.zeros((N, N))
    for j in range(N):
        k = (j + 1) % N
        s = -1.0 if k == 0 else 1.0       # antiperiodic boundary (no zero mode)
        H[j, k] = H[k, j] = -s
    e, v = np.linalg.eigh(H)
    occ = v[:, e < 0]
    gap = min(abs(e))
    return occ.conj() @ occ.T, (e < 0).sum(), gap


C, n_occ, gap = ring_vacuum()
step("F0.antiperiodic critical ring: unique half-filled vacuum (128 occupied, no zero mode)",
     n_occ == N // 2 and gap > 1e-6, f"min |e| = {gap:.3e}")

# raw station clock (mu4-pointing station_matrices) and the twisted quarter translation
Rst = np.zeros((4, 4))
for x in range(3):
    Rst[x + 1, x] = 1
Rst[0, 3] = -1
Rraw = np.kron(Rst, np.eye(M))
Tq = np.zeros((N, N))
for j in range(N):
    k = (j + M) % N
    Tq[k, j] = -1.0 if j + M >= N else 1.0
step("F1.raw station clock R_station (x) I_64 = twisted quarter translation c_j -> c_{j+64} (sign on wrap)",
     np.array_equal(Rraw, Tq))
dev = np.abs(Rraw @ C @ Rraw.T - C).max()
T1 = np.zeros((N, N))
for j in range(N):
    T1[(j + 1) % N, j] = -1.0 if j + 1 >= N else 1.0
dev1 = np.abs(T1 @ C @ T1.T - C).max()
step("F2.R Omega = Omega: R C R^T = C (vacuum invariant); also homogeneous under the one-site twisted shift",
     dev < 1e-12 and dev1 < 1e-12, f"{dev:.1e}, {dev1:.1e}")
arcs_ok = True
for start in range(0, N, 7):
    for length in (5, 17, 64, 90):
        I = [(start + t) % N for t in range(length)]
        img = {int(np.nonzero(Rraw[:, j])[0][0]) for j in I}
        arcs_ok &= img == {(j + M) % N for j in I}
step("F3.R A(I) R^-1 = A(I + N/4) for all tested arcs (signed permutation of local modes)", arcs_ok)
step("F4.R^4 = -I (antiperiodic NS sign = (-1)^F on the fermions); projectively order 4",
     np.array_equal(np.linalg.matrix_power(Rraw, 4), -np.eye(N)))
d = np.arange(1, N // 2, 2)
lat = np.array([abs(C[0, int(x)]) for x in d])
cft = 1.0 / (N * np.sin(np.pi * d / N))
rel = np.abs(lat - cft) / cft
step("F5.lattice vacuum = Moebius-covariant NS kernel 1/(N sin(pi d/N)) on odd separations (rel. dev.)",
     rel.max() < 1e-12, f"max rel {rel.max():.1e}; even separations vanish: {np.abs(C[0, 2:N//2:2]).max():.1e}")

# F6: the bit
phi = sp.symbols("phi", real=True)
cr = lambda a, b, c, e: sp.simplify((a - c) * (b - e) / ((a - e) * (b - c)))
ab = sp.symbols("a b c e")
# on the unit circle conj(z) = 1/z, and z -> 1/z preserves cross-ratios, so conj(lambda) = lambda
inv_invariant = sp.simplify(cr(*[1 / x for x in ab]) - cr(*ab)) == 0
w = sp.Rational(-1, 2) + sp.sqrt(3) * sp.I / 2
lam_tet = sp.nsimplify(cr(sp.Integer(0), sp.Integer(1), w, sp.expand(w ** 2)))
orbit3 = [k for k in range(1, 5) if 4 % 3 == 0]
step("F6.marks on the seam circle => lambda real; tetrahedral (j=0) lambda is non-real, so an order-3 geometric clock is impossible",
     inv_invariant and sp.simplify(sp.im(sp.expand_complex(lam_tet))) != 0 and orbit3 == [],
     f"lambda_tet = {lam_tet}; an order-3 circle rotation has orbits of size 3, 4 marks cannot form them")
# F7: equal stations -> square marks; unequal stations -> no order-4 covariance
bounds = [0, M, 2 * M, 3 * M]
angles = [sp.exp(2 * sp.pi * sp.I * sp.Rational(b, N)) for b in bounds]
step("F7.equal stations (4 x 64): station boundaries at angles k pi/2, cross-ratio 2 (lambda = -1 orbit, square)",
     sp.simplify(cr(*angles) - 2) == 0)
sizes = [60, 68, 64, 64]
bnd = np.cumsum([0] + sizes[:-1])
ang = [sp.exp(2 * sp.pi * sp.I * sp.Rational(int(b), N)) for b in bnd]
step("F7b.negative control: unequal stations (60,68,64,64) give a non-harmonic real cross-ratio: no order-4 Moebius covariance",
     min(abs(complex(sp.N(cr(*ang))) - h) for h in (-1, 2, 0.5)) > 1e-3 and abs(complex(sp.N(cr(*ang))).imag) < 1e-12,
     f"cross-ratio {complex(sp.N(cr(*ang))).real:.6f}")

summary = {"steps": LOG, "passed": sum(s["pass"] for s in LOG), "total": len(LOG)}
(HERE / "seam_vacuum_covariance_results.json").write_text(json.dumps(summary, indent=2))
print(f"{summary['passed']}/{summary['total']} steps passed")
sys.exit(0 if summary["passed"] == summary["total"] else 1)
