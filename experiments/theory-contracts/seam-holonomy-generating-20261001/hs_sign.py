#!/usr/bin/env python3
"""Part HS: CHIRAL4D.HS_SIGN.01 -- is the Hubbard-Stratonovich decoupling of Theta_1 sign-free by a Kramers symmetry?

Theta_1 (frozen): Weyl block, H_walk = (i/2)(W - W^dag) (local, same eight bcc steps, Hermitian, nodes where W = +-1),
frozen filter F (real, point-group invariant), H_int = -g sum O_A^2 decoupled with ten real fields phi_A(x, tau):
e^{dtau g O_A^2} ~ int dphi exp(-phi^2/(4 g dtau) + phi O_A).  Kramers criterion: an antiunitary Theta = u K with
Theta^2 = -1 that maps EVERY decoupled single-particle operator to itself (walk, all ten phi_A O_A, filter) makes the
Bogoliubov transfer matrix Kramers-paired, det(1 + B) = prod |1 + lambda|^2 >= 0 and the Pfaffian sign constant on the
connected complement of codimension-2 zeros.
  HS1 Hamiltonian walk: h_eps = (i/8) eps.sigma - (1/8) nu_eps.sigma with nu_eps = (eps2 eps3, -eps1 eps3, eps1 eps2) (the
      ordering terms of the product walk); spectrum zero exactly at W = +-1: k = 0 (light), +P, H, -P (all filtered,
      f = 1, 2, 1); the time-reversal-friendly part alone (nu dropped) has nodal LINES.
  HS2 internal part: a Spin(10)-compatible antiunitary must map psi to psi^dag (16 is not self-conjugate; Schur), the
      pairing condition det(b_s) Y^T B_A^* Y = -B_A^* is solvable with Y = 1 iff det b_s = -1 (B_A real/imaginary
      pattern checked), and number-preserving Theta would need an intertwiner 16 -> 16bar (none).
  HS3 search: for all 48 point-group elements g (signed axis permutations, acting on positions) and both Nambu types,
      the spin conditions are linear; NO solution with Theta^2 = -1 exists (exact rank computation): the Kramers route
      fails for the BDT Hamiltonian walk; the obstruction is the ordering term nu_eps, the very term that removes the
      nodal lines.
  NOTE (Part SR, sign_remedy.py): lattice maps g are NOT admissible for a term-wise HS argument (they send phi(x) to
      phi(gx)); the search above is a superset, its negative result stands and is tightened to the whole class there.
  HS4 (statement) no sign-free guarantee for Theta_1 from Kramers symmetry; options: a different local H_walk with the
      same light node and filtered doublers that admits Theta, Majorana-class criteria (two anticommuting Majorana
      time reversals), or measuring the average sign; QMC for Gate B is NOT started.
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/hs_sign.py
"""
from __future__ import annotations

import json
import sys
from itertools import permutations, product
from pathlib import Path

import numpy as np
from scipy.linalg import null_space

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    sys.stdout.flush()


sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1.0 + 0j, -1]); I2 = np.eye(2)
S = [sx, sy, sz]
signs = [np.array(s) for s in product((-1, 1), repeat=3)]
def A_eps(e):
    return (I2 + e[0] * sx) / 2 @ (I2 + e[1] * sy) / 2 @ (I2 + e[2] * sz) / 2
h = {tuple(e): 0.5j * (A_eps(e) - A_eps(-e).conj().T) for e in signs}
nu = lambda e: np.array([e[1] * e[2], -e[0] * e[2], e[0] * e[1]])
form_ok = all(np.allclose(h[tuple(e)], sum(((1j / 8) * e[i] - (1 / 8) * nu(e)[i]) * S[i] for i in range(3))) for e in signs)
herm_ok = all(np.allclose(h[tuple(-e)], h[tuple(e)].conj().T) for e in signs)
Hk = lambda k: sum(np.exp(-1j * e @ k) * h[tuple(e)] for e in signs)
Wk = lambda k: (np.cos(k[0]) * I2 - 1j * np.sin(k[0]) * sx) @ (np.cos(k[1]) * I2 - 1j * np.sin(k[1]) * sy) @ (np.cos(k[2]) * I2 - 1j * np.sin(k[2]) * sz)
rng = np.random.default_rng(61)
match = all(np.allclose(Hk(k), 0.5j * (Wk(k) - Wk(k).conj().T)) for k in rng.normal(size=(5, 3)))
P = np.pi / 2 * np.ones(3); H = np.array([np.pi, 0, 0])
nodes = all(np.allclose(Hk(k), 0) for k in (np.zeros(3), P, -P, H))
f = lambda k: 1 - np.cos(k[0]) * np.cos(k[1]) * np.cos(k[2])
Hsym = lambda k: sum(np.exp(-1j * e @ k) * sum((1j / 8) * e[i] * S[i] for i in range(3)) for e in signs)
line = all(np.allclose(Hsym(np.array([np.pi / 2, np.pi / 2, t])), 0) for t in np.linspace(-np.pi, np.pi, 9))
line_bdt = max(np.abs(Hk(np.array([np.pi / 2, np.pi / 2, t]))).max() for t in np.linspace(-np.pi, np.pi, 9))
step("HS1.Hamiltonian walk: h_eps = (i/8) eps.sigma - (1/8) nu_eps.sigma (nu = ordering terms), Hermitian, = (i/2)(W - W^dag); zeros at k = 0 and the filtered doublers +P, -P, H (f = 1, 1, 2); the time-reversal-friendly part alone vanishes on the whole line (pi/2, pi/2, k_z)",
     form_ok and herm_ok and match and nodes and line and line_bdt > 1e-3 and np.allclose([f(P), f(-P), f(H)], [1, 1, 2]),
     f"max |H_BDT| on that line {line_bdt:.3f} (the ordering terms lift the nodal line)")

# HS2 ----------------------------------------------------------------------------------------------------
Zq = np.diag([1.0, -1.0]).astype(complex); lowc = np.array([[0, 1], [0, 0]], dtype=complex)
c5 = []
for j in range(5):
    mats = [Zq] * j + [lowc] + [I2.astype(complex)] * (4 - j)
    out = mats[0]
    for m_ in mats[1:]:
        out = np.kron(out, m_)
    c5.append(out)
Gm = []
for j in range(5):
    Gm += [c5[j] + c5[j].conj().T, -1j * (c5[j] - c5[j].conj().T)]
occ = np.array([[(b >> (4 - j)) & 1 for j in range(5)] for b in range(32)])
even_idx = [b for b in range(32) if occ[b].sum() % 2 == 0]
Pe = np.eye(32)[:, even_idx]
Cc = Gm[0] @ Gm[2] @ Gm[4] @ Gm[6] @ Gm[8]
BA = [Pe.T @ (Cc @ Gm[A]) @ Pe for A in range(10)]
Sig = [Pe.T @ (0.25j * (Gm[A] @ Gm[B] - Gm[B] @ Gm[A])) @ Pe for A in range(10) for B in range(A + 1, 10)]
real_pattern = [("real" if np.allclose(b.imag, 0) else "imag" if np.allclose(b.real, 0) else "complex") for b in BA]
ct_Y1 = all(np.allclose(b.conj(), b) or np.allclose(b.conj(), -b) for b in BA)
# CT type with Y = 1: det(b_s) B_A^* = -B_A^*  -> det b_s = -1 for every A (independent of the real/imag pattern)
# number-preserving: Y^T B^* Y = c B with Y commuting with Spin(10) up to conjugation: 16 vs 16bar test
gens16 = Sig
inter = null_space(np.array([np.concatenate([(Yb @ g.conj() - g @ Yb).ravel() for g in gens16]) for Yb in [np.eye(16)[:, [i]] @ np.eye(16)[[j], :] for i in range(16) for j in range(16)]]).T, rcond=1e-9)
no_16_16bar = inter.shape[1] == 0
step("HS2.internal part: no intertwiner between the 16 and its conjugate (16 is not self-conjugate), so a Spin(10)-compatible antiunitary must be particle-hole type psi -> b psi^dag; with Y = 1 its pairing condition det(b_s) Y^T B_A^* Y = -B_A^* holds for all ten channels iff det b_s = -1",
     no_16_16bar and ct_Y1, f"B_A pattern {real_pattern}")

# HS3 ----------------------------------------------------------------------------------------------------
pg = []
for p in permutations(range(3)):
    for sg in product((-1, 1), repeat=3):
        R = np.zeros((3, 3))
        for i in range(3):
            R[i, p[i]] = sg[i]
        pg.append(R)
def solve(kind, R):
    """linear conditions on the 2x2 spin matrix X (vectorised, complex): returns basis of solutions."""
    rows = []
    for e in signs:
        he = h[tuple(e)]
        hg = h[tuple((R @ e).astype(int))]
        for i, j in product(range(2), range(2)):
            E = np.zeros((2, 2), complex); E[i, j] = 1
        if kind == "N":       # X^dag h_eps^* X = h_{g eps}  <=>  h_eps^* X - X h_{g eps} = 0
            M = np.kron(he.conj(), I2) - np.kron(I2, hg.T)
        else:                 # CT: X^dag h_eps^* X = - h_{g eps}^*  <=>  h_eps^* X + X h_{g eps}^* = 0
            M = np.kron(he.conj(), I2) + np.kron(I2, hg.conj().T)
        rows.append(M)
    return null_space(np.vstack(rows), rcond=1e-10)
found, report = [], {}
for kind in ("N", "CT"):
    cnt = 0
    for R in pg:
        ns = solve(kind, R)
        if ns.shape[1] == 0:
            continue
        # look for X in the solution space with X unitary, and Theta^2 = -1 (N: X^* X = -1 requires g^2 = 1; CT: X^2 = -1), det X = -1 for CT
        for c in range(ns.shape[1]):
            X = ns[:, c].reshape(2, 2)
            X = X / np.sqrt(abs(np.linalg.det(X))) if abs(np.linalg.det(X)) > 1e-12 else X
            unitary = np.allclose(X.conj().T @ X, np.eye(2), atol=1e-8)
            if kind == "N":
                sq = np.allclose(X.conj() @ X, -np.eye(2), atol=1e-8) and np.allclose(R @ R, np.eye(3))
                ok = unitary and sq and no_16_16bar is False
            else:
                sq = np.allclose(X @ X, -np.eye(2), atol=1e-8) and np.allclose(R @ R, np.eye(3))
                ok = unitary and sq and abs(np.linalg.det(X) + 1) < 1e-8
            if ok:
                found.append((kind, R.tolist()))
        cnt += ns.shape[1]
    report[kind] = cnt
# exact reduction: N type is excluded by HS2; CT type with Y in {1, i 1} (the only internal unitaries preserving all
# ten B_A up to a common factor) needs X^2 = -det... -> the spin part X must be a scalar (eigenvalue argument), and then
# the walk condition reads h_{g eps} = -h_eps for every eps:
scalar_ok = [R.tolist() for R in pg if all(np.allclose(h[tuple((R @ e).astype(int))], -h[tuple(e)]) for e in signs)]
eig_arg = True
for Y2, detc in ((1, -1), (-1, +1)):          # Y = 1: det X = -1, X^2 = -1 ;  Y = i 1: det X = +1, X^2 = +1 (beta^2 = X^2 (x) Y^2 = -1)
    # X^2 = -Y2 and det X = detc force X = scalar: eigenvalues x1, x2 with x^2 = -Y2 and x1 x2 = detc
    roots = np.roots([1, 0, Y2])
    pairs = [(a_, b_) for a_ in roots for b_ in roots if abs(a_ * b_ - detc) < 1e-12]
    eig_arg &= all(abs(a_ - b_) < 1e-12 for a_, b_ in pairs) and len(pairs) > 0
step("HS3.search over all 48 point-group elements g: particle-number-preserving Theta is excluded (HS2); for the particle-hole type the conditions Theta^2 = -1 and the pairing determinant force a SCALAR spin part (eigenvalue argument), and then the walk would need h_{g eps} = -h_eps for all eps -- satisfied by no g (inversion fails on the ordering term nu_eps); the linear null spaces confirm it: the Kramers route to a sign-free decoupling fails for the BDT Hamiltonian walk",
     len(found) == 0 and len(scalar_ok) == 0 and eig_arg, f"linear solution dimensions {report}; g with h_(g eps) = -h_eps: {len(scalar_ok)}; admissible Theta: {len(found)}")
step("HS4.(statement) no sign-free guarantee for Theta_1 from a Kramers symmetry; options: a different local walk Hamiltonian with the same light node and filtered doublers that admits Theta (the time-reversal-friendly part alone has nodal lines), Majorana-class criteria, or measuring the average sign; QMC for Gate B is not started",
     True, "Gate B stays open")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "solution_dims": report, "B_pattern": real_pattern}
(HERE / "hs_sign_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
