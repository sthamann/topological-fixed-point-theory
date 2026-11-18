#!/usr/bin/env python3
"""Part P: AX.ACCESS from the existing microscopic CAR field (QWZ edge, actual filled sea).

The access is NOT postulated: Psi_N(g) = a(R_N g) of microscopic-charged-car-limit (exactly local on
the top boundary row, original r=1 seam holonomy) in the actual negative-energy sea of the original
width-8, mass-1 QWZ cylinder.  Single-site insertions at the clock-orbit points (D.LAMBDA).
  P1 log|<Psi(x) Psi(y)^*>| obeys the chordal law: after removing its circle mean (the zero mode of
     |D|, i.e. the canonical reference) it converges to -ln|2 sin(pi d/N)|.
  P2 the charge-neutral four-point function of CAR insertions is a Wick determinant whose modulus is
     exp(sum_{i<j} e_i e_j G_ij): the neutral port kernel is the Hessian of log charged correlators.
  P3 the port kernel M = (1/pi) G at the four orbit points with zero diagonal: E_+ M E_+ -> rank one
     along the uniform direction, M -> -(ln2/2pi) circ(0,1,2,1): the uniform port is an OUTPUT of
     the microscopic access + zero-mode reference.
  P4 control: a bulk row has exponentially small correlations (no log law, no port).
  P5 control: the bottom row (opposite chirality) gives the same modulus law.
    PYTHONPATH=. python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/car_access.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import scipy.linalg

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


spec = importlib.util.spec_from_file_location(
    "car_limit", ROOT / "experiments/theory-contracts/microscopic-charged-car-limit/checker.py")
car = importlib.util.module_from_spec(spec)
spec.loader.exec_module(car)
_, source = car.inherited()                       # validates the source pins


def sea(n):
    h = np.ascontiguousarray(np.asarray(source.qwz_cylinder(n, 8, 1, 1), dtype=complex))
    ev, v = scipy.linalg.eigh(h, driver="evr")            # full source; numpy's zheevd fails at 16N = 4096
    occ = ev < 0
    assert int(occ.sum()) == 8 * n
    return v[:, occ] @ v[:, occ].conj().T


def site_vector(n, x, row="top"):
    samples = np.zeros(n, complex)
    samples[x] = 1.0
    if row == "top":
        return car.local_smearing(n, samples)              # original R_N g, top row, r = 1 holonomy
    gauge = np.exp(-.5j * np.pi * np.arange(n) / n)
    spin = np.zeros(16, complex)
    if row == "bottom":
        spin[:2] = np.array([1, -1]) / np.sqrt(2)
    else:                                                   # a bulk row y = 3
        spin[6:8] = np.array([1, 1]) / np.sqrt(2)
    return np.kron(gauge * samples / np.sqrt(n), spin)


def two_point(n, P, row="top"):
    V = np.column_stack([site_vector(n, x, row) for x in range(n)])
    return V.conj().T @ P @ V                               # <a^*(v_x) a(v_y)>-type hole covariance


res = {}
for n in (64, 128, 256):
    P = sea(n)
    Ct = two_point(n, P, "top")
    d = np.arange(1, n)
    L = np.log(np.abs(Ct[0, 1:]))
    G = L - L.mean() - np.log(n) / (n - 1)                  # zero-mode (circle-mean) reference
    exact = -np.log(np.abs(2 * np.sin(np.pi * d / n)))
    far = (d >= n // 8) & (d <= n - n // 8)
    res[n] = dict(P=P, C=Ct, G=G, dev=np.abs(G - exact)[far].max())
devs = [res[n]["dev"] for n in (64, 128, 256)]
step("P1.microscopic CAR two-point on the top row: log|C| minus its circle mean -> -ln|2 sin(pi d/N)| (chordal law)",
     devs[-1] < 2e-2 and devs[-1] < devs[0], "max dev (d in [N/8, 7N/8]) " + ", ".join(f"{x:.1e}" for x in devs))

n = 256
Ct, G = res[n]["C"], res[n]["G"]
pos = [5 + k * n // 4 for k in range(4)]
ket, bra = [pos[1], pos[3]], [pos[0], pos[2]]               # Psi^* at pos0, pos2; Psi at pos1, pos3
W = np.array([[Ct[a, b] for b in ket] for a in bra])
wick = abs(np.linalg.det(W))
eps = [1, -1, 1, -1]
Gfun = lambda a, b: G[(b - a) % n - 1] if a != b else 0.0
quad = sum(eps[i] * eps[j] * Gfun(pos[i], pos[j]) for i in range(4) for j in range(i + 1, 4))
norm = abs(Ct[pos[0], pos[1]]) / np.exp(Gfun(pos[0], pos[1]))   # common normalisation of one contraction
step("P2.neutral four-point of CAR insertions: |Wick det| = norm^2 exp(-sum_{i<j} e_i e_j G_ij) (Cauchy identity on microscopic data)",
     abs(np.log(wick) - (2 * np.log(norm) - quad)) < 3e-2,
     f"log|det| = {np.log(wick):.5f}, from quadratic form {2 * np.log(norm) - quad:.5f}")
M = np.array([[Gfun(a, b) / np.pi for b in pos] for a in pos])
s = np.ones(4) / 2
v = np.array([1, -1, 1, -1]) / 2
Q4 = np.roll(np.eye(4), 1, axis=0)
Eplus = (np.eye(4) + Q4 @ Q4) / 2
Cm = np.array([[0, 1, 2, 1], [1, 0, 1, 2], [2, 1, 0, 1], [1, 2, 1, 0]])
ratio = abs((v @ M @ v) / (s @ M @ s))
step("P3.port from the CAR access: E_+ M E_+ rank one along s (lam_v/lam_s -> 0), M -> -(ln2/2pi) circ(0,1,2,1): uniform port is an OUTPUT",
     ratio < 2e-2 and np.allclose(M, -(np.log(2) / (2 * np.pi)) * Cm, atol=5e-3) and np.abs(Q4 @ M - M @ Q4).max() < 5e-3,
     f"lam_v/lam_s = {ratio:.2e}; M/(-ln2/2pi) row0 = {np.round(M[0] / (-(np.log(2) / (2 * np.pi))), 3).tolist()}")
Cb = two_point(n, res[n]["P"], "bulk")
Cd = two_point(n, res[n]["P"], "bottom")
far_bulk = np.abs(Cb[0, n // 8: n - n // 8]).max()
Ld = np.log(np.abs(Cd[0, 1:]))
Gd = Ld - Ld.mean() - np.log(n) / (n - 1)
dd = np.arange(1, n)
far = (dd >= n // 8) & (dd <= n - n // 8)
step("P4.control bulk row y=3: correlations beyond N/8 are exponentially small (no log law, no port)",
     far_bulk < 1e-3 * np.abs(Ct[0, n // 8: n - n // 8]).min() or far_bulk < 1e-8, f"max |C_bulk| far = {far_bulk:.1e}")
step("P5.control bottom row (opposite chirality): same chordal modulus law",
     np.abs(Gd - (-np.log(np.abs(2 * np.sin(np.pi * dd / n)))))[far].max() < 2e-2)

# P6: the declared source itself satisfies the abstract state axioms used in Parts H/K/M:
#     homogeneous (|C(x,x+d)| independent of x), reflection symmetric (|C(d)| = |C(-d)|),
#     half filled (rank 8N of 16N), quarter-translation invariant up to the r=1 gauge phase;
#     positive energy Hcal_N >= 0 is proven in microscopic-charged-car-limit (README section 1).
absC = np.abs(Ct)
homog = max(np.abs(absC[x, (x + dd) % n] - absC[0, dd % n]).max() for x in range(0, n, 17) for dd in [np.arange(1, n)])
refl = np.abs(absC[0, 1:] - absC[0, 1:][::-1]).max()
Tq = np.roll(np.eye(n), n // 4, axis=0)
phase = np.abs(np.abs(Tq @ Ct @ Tq.T) - absC).max()
step("P6.the declared QWZ source satisfies K0 (homogeneous), K2 (reflection), K1 (half filling 8N/16N) and quarter-translation covariance; positive energy is proven in the CAR contract",
     homog < 1e-10 and refl < 1e-10 and phase < 1e-10 and np.isclose(np.trace(res[n]["P"]).real, 8 * n),
     f"homogeneity {homog:.1e}, reflection {refl:.1e}, quarter shift {phase:.1e}")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "P1_devs": devs}
(HERE / "car_access_results.json").write_text(json.dumps(out, indent=2))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
