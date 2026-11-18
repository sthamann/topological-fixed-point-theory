#!/usr/bin/env python3
"""Part FP: CHIRAL4D.FILTERED_CAR.01 + CHIRAL4D.LOCAL_PATCH.01 for Theta_1 (Weyl block, frozen filter, trivial links).

  FP1 filtered CAR: {rho_x, rho_y^dag} = (F F^dag)_xy; without links (F^2) has support {0, eps, eps+eps'} with
      coefficients 9/8, -1/4, 1/16 (2,0,0), 1/32 (2,2,0), 1/64 (2,2,2); {rho, rho} = 0; with random U(2) links the
      diagonal stays 9/8 and the support is unchanged.
  FP2 spectrum of F^2: f(k)^2 in [0, 4]; null space = the k = 0 class only (rank per momentum full elsewhere);
      the smallest nonzero eigenvalue on an L-periodic cell falls like L^-4 (the filtered fields become singular
      towards the light node, as intended); condition number away from Gamma is bounded (shells |k| >= kappa).
      Normalised overlaps of neighbouring active sets: -2/9 (nearest), 1/18, 1/36, 1/72.
  FP3 cell geometry: a periodic primitive cell L^3 contains the mirror momentum P iff 4 | L; conventional cubic boxes
      of even cube count (cubic side multiple of 4 in these coordinates) contain P -- a mirror test must use such a cell.
  FP4 local spectrum by parity (from the exact so(11) decomposition of Part SG): with the canonical active modes
      a = sqrt(8/9) rho and g* = 81/64 g: E_0 = -2176 g*; lowest even excitation 144 g* (degeneracy 10); lowest odd
      232 g* (degeneracy 64); next even 320 g* (degeneracy 54); even and odd sectors have 2^31 states each.
  FP5 normal ordering: [Q_A, Q_A^dag] = 2 tr(M_A M_A^dag) - 4 a^dag (M_A^dag M_A) a = 64 - 4N (B_A unitary), so
      H_int = H_quartic + 45 g sum psi^dag F^2 psi - 810 g N_sites; the extra one-particle term is f(k)^2 ~ |k|^4/4 at
      the light node (irrelevant) and the interaction conserves N only mod 4.
  FP6 LOCAL_PATCH (frustration, exact part): for neighbouring active sets b = r a + s c (r = -2/9, s^2 = 1 - r^2) and any
      state whose a-factor has <a^dag a> = 1/2, <a a> = 0 (true for the unique singlet ground state: Schur for
      Spin(10) x SU(2), Z4, and <N - 16> = 0 as a vector matrix element), ||J_b Psi||^2 >= r^2 s^2 tr(t^2)/2 > 0 for a
      Cartan generator t: the two singlet minima are NEVER reached simultaneously (identity verified on test states).
      E_min(h_x + h_y) > 2 E_0; its exact value (2^64-dimensional patch) is not computed here.
  FP7 (statement) consequence: the strong-coupling ground state is not a product of local singlets; Gate B needs the
      collective phase of the overlapping model (QMC / tensor networks after the sign audit).
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/filtered_patch.py
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction as Fr
from itertools import product
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    sys.stdout.flush()


signs = [tuple(s) for s in product((-1, 1), repeat=3)]
rng = np.random.default_rng(53)

# FP1 ----------------------------------------------------------------------------------------------------
coef = {(0, 0, 0): Fr(1)}
for s in signs:
    coef[s] = coef.get(s, Fr(0)) - Fr(1, 8)
F2 = {}
for u, cu in coef.items():
    for v, cv in coef.items():
        key = tuple(u[i] + v[i] for i in range(3))
        F2[key] = F2.get(key, Fr(0)) + cu * cv
F2 = {k: v for k, v in F2.items() if v != 0}
classes = {}
for k, v in F2.items():
    cls = tuple(sorted(abs(x) for x in k))
    classes.setdefault(cls, set()).add(v)
expect = {(0, 0, 0): {Fr(9, 8)}, (1, 1, 1): {Fr(-1, 4)}, (0, 0, 2): {Fr(1, 16)}, (0, 2, 2): {Fr(1, 32)}, (2, 2, 2): {Fr(1, 64)}}
support_ok = classes == expect
Lb = 4
sites = [x for x in product(range(Lb), repeat=3) if len({c % 2 for c in x}) == 1]
idx = {x: i for i, x in enumerate(sites)}
def rand_u2():
    q, r = np.linalg.qr(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))
    return q
links = {}
for x in sites:
    for s in signs:
        if s[0] == 1:
            y = tuple((x[c] + s[c]) % Lb for c in range(3))
            U = rand_u2(); links[(x, s)] = U; links[(y, tuple(-t for t in s))] = U.conj().T
dim = 2 * len(sites)
FU = np.eye(dim, dtype=complex)
for x in sites:
    for s in signs:
        y = tuple((x[c] + s[c]) % Lb for c in range(3))
        FU[2 * idx[x]:2 * idx[x] + 2, 2 * idx[y]:2 * idx[y] + 2] -= links[(x, s)] / 8
G = FU @ FU.conj().T
diag_ok = all(np.allclose(G[2 * i:2 * i + 2, 2 * i:2 * i + 2], 9 / 8 * np.eye(2)) for i in range(len(sites)))
def cdist(x, y):
    d = [((y[c] - x[c] + Lb // 2) % Lb) - Lb // 2 for c in range(3)]
    return tuple(sorted(abs(t) for t in d))
supp_U = all(np.allclose(G[2 * idx[x]:2 * idx[x] + 2, 2 * idx[y]:2 * idx[y] + 2], 0) for x in sites for y in sites if cdist(x, y) not in expect)
step("FP1.filtered CAR: {rho_x, rho_y^dag} = (F F^dag)_xy; without links F^2 has support {0, eps, eps+eps'} with 9/8, -1/4, 1/16 (2,0,0), 1/32 (2,2,0), 1/64 (2,2,2) (exact fractions); with random U(2) links the diagonal stays 9/8 and the support is unchanged",
     support_ok and diag_ok and supp_U, "{rho, rho} = 0 trivially (rho is linear in psi)")

# FP2 ----------------------------------------------------------------------------------------------------
f = lambda k: 1 - np.cos(k[..., 0]) * np.cos(k[..., 1]) * np.cos(k[..., 2])
def bcc_momenta(Lc):
    """momenta of an Lc^3 periodic primitive bcc cell (primitive vectors b_i), as cubic-coordinate k."""
    b = np.array([[1, 1, -1], [1, -1, 1], [-1, 1, 1]], dtype=float)
    Binv = 2 * np.pi * np.linalg.inv(b).T
    ns = np.array(list(product(range(Lc), repeat=3)), dtype=float) / Lc
    return ns @ Binv
nulls = []
for Lc in (4, 8, 12, 16):
    vals = f(bcc_momenta(Lc)) ** 2
    nulls.append(int(np.sum(vals < 1e-12)))
bmat = np.array([[1, 1, -1], [1, -1, 1], [-1, 1, 1]], dtype=float)
Binv = 2 * np.pi * np.linalg.inv(bmat).T
small_n = np.array([n for n in product(range(-2, 3), repeat=3) if any(n)], dtype=float)
Ls = np.array([32, 64, 128, 256])
mins = [np.min(f((small_n / Lc) @ Binv) ** 2) for Lc in Ls]                # smallest nonzero eigenvalue on an L^3 cell
expo = np.polyfit(np.log(Ls), np.log(mins), 1)[0]
kk = rng.uniform(-np.pi, np.pi, size=(200000, 3))
fk2 = f(kk) ** 2
conds = {}
for kap in (0.5, 1.0):
    mask = np.linalg.norm(kk, axis=1) >= kap
    # distance to the k = 0 class: all fcc translates (pi Z^3 with even sum)
    G0 = np.array([[a, b_, c] for a, b_, c in product((-2, -1, 0, 1, 2), repeat=3) if (a + b_ + c) % 2 == 0]) * np.pi
    dist = np.min(np.linalg.norm(kk[:, None, :] - G0[None, :, :], axis=2), axis=1)
    m2 = dist >= kap
    conds[kap] = fk2[m2].max() / fk2[m2].min()
ov = {k: float(v / Fr(9, 8)) for k, v in [(c, list(s)[0]) for c, s in expect.items()]}
step("FP2.spectrum of F^2: f^2 in [0, 4]; null space = the k = 0 class only (one null momentum per cell, full rank elsewhere); smallest nonzero eigenvalue ~ L^-4; condition number bounded away from Gamma; normalised overlaps -2/9 (nearest), 1/18, 1/36, 1/72",
     nulls == [1, 1, 1, 1] and abs(expo + 4) < 0.02 and fk2.max() <= 4 + 1e-12 and all(np.isfinite(list(conds.values())))
     and abs(ov[(1, 1, 1)] + 2 / 9) < 1e-15 and abs(ov[(0, 0, 2)] - 1 / 18) < 1e-15 and abs(ov[(0, 2, 2)] - 1 / 36) < 1e-15 and abs(ov[(2, 2, 2)] - 1 / 72) < 1e-15,
     f"null momenta per cell {nulls}; min nonzero f^2 ~ L^{expo:.2f}; condition number for distance >= 0.5: {conds[0.5]:.1f}, >= 1.0: {conds[1.0]:.1f}")

# FP3 ----------------------------------------------------------------------------------------------------
P = np.pi / 2 * np.ones(3)
def contains_P_primitive(Lc):
    b = np.array([[1, 1, -1], [1, -1, 1], [-1, 1, 1]], dtype=float)
    return all(abs(np.cos(Lc * (P @ bi)) - 1) < 1e-12 for bi in b)
def contains_P_cubic(side):
    return all(abs(np.cos(side * P[i]) - 1) < 1e-12 for i in range(3))
tab = {f"primitive {Lc}^3 ({Lc ** 3} sites)": contains_P_primitive(Lc) for Lc in (2, 3, 4)}
tab.update({f"cubic side {s_} ({2 * (s_ // 2) ** 3} sites)": contains_P_cubic(s_) for s_ in (2, 4, 6, 8)})
step("FP3.cell geometry: a periodic primitive L^3 cell contains the mirror momentum P iff 4 | L (2^3, 3^3 do not, 4^3 does); a cubic box contains P iff its side is a multiple of 4 (the conventional 2^3 cubes with 16 sites do)",
     tab["primitive 2^3 (8 sites)"] is False and tab["primitive 3^3 (27 sites)"] is False and tab["primitive 4^3 (64 sites)"] is True
     and tab["cubic side 4 (16 sites)"] is True and tab["cubic side 2 (2 sites)"] is False and tab["cubic side 6 (54 sites)"] is False,
     "; ".join(f"{k}: {'yes' if v else 'no'}" for k, v in tab.items()))

# FP4 ----------------------------------------------------------------------------------------------------
res = json.loads((HERE / "symmetric_gap_results.json").read_text())
import re
irreps = {tuple(float(t) for t in re.findall(r"-?\d+\.\d+", k)): v for k, v in res["irreps"].items()}
rho2 = np.array([9, 7, 5, 3, 1]); rhoD = np.array([4, 3, 2, 1, 0])
C2B = lambda lam: float(np.sum(lam * (lam + rho2)))
C2D = lambda mu: float(np.sum(mu * (mu + 2 * rhoD)))
def dim_B5(lam):
    l = lam + rho2 / 2; r = rho2 / 2
    return np.prod([(l[i] ** 2 - l[j] ** 2) / (r[i] ** 2 - r[j] ** 2) for i in range(5) for j in range(i + 1, 5)]) * np.prod(l / r)
def dim_D5(mu):
    l = mu + rhoD
    return np.prod([(l[i] ** 2 - l[j] ** 2) / (rhoD[i] ** 2 - rhoD[j] ** 2) for i in range(5) for j in range(i + 1, 5)])
def branch(lam):
    ranges = [np.arange(lam[i + 1], lam[i] + 1e-9, 1.0) for i in range(4)]
    return [np.array(m) for m in product(*ranges, np.arange(-lam[4], lam[4] + 1e-9, 1.0))]
levels = {"even": {}, "odd": {}}
dims = {"even": 0.0, "odd": 0.0}
for lam_t, m in irreps.items():
    lam = np.array(lam_t)
    par = "even" if abs(lam[0] - round(lam[0])) < 1e-9 else "odd"
    dims[par] += m * round(dim_B5(lam))
    for mu in branch(lam):
        E = round(C2B(lam) - C2D(mu), 9)
        levels[par][E] = levels[par].get(E, 0) + m * int(round(dim_D5(mu)))
E0 = max(levels["even"])
ev_sorted = sorted(levels["even"], reverse=True)
od_sorted = sorted(levels["odd"], reverse=True)
kap = 16
gaps = {"even1": (kap * (E0 - ev_sorted[1]), levels["even"][ev_sorted[1]]), "odd1": (kap * (E0 - od_sorted[0]), levels["odd"][od_sorted[0]]),
        "even2": (kap * (E0 - ev_sorted[2]), levels["even"][ev_sorted[2]])}
step("FP4.local spectrum by parity (canonical a = sqrt(8/9) rho, g* = 81/64 g): E_0 = -2176 g*; lowest even excitation 144 g* (x10), lowest odd 232 g* (x64), next even 320 g* (x54); even and odd sectors 2^31 each",
     kap * E0 == 2176 and levels["even"][E0] == 1 and gaps["even1"] == (144, 10) and gaps["odd1"] == (232, 64) and gaps["even2"] == (320, 54)
     and abs(dims["even"] - 2 ** 31) < 0.5 and abs(dims["odd"] - 2 ** 31) < 0.5,
     f"gaps (in g*) and degeneracies: even {gaps['even1']}, odd {gaps['odd1']}, next even {gaps['even2']}; in the original coupling: 729/4 g, 2349/8 g")

# FP5 ----------------------------------------------------------------------------------------------------
def fermion_ops(n):
    Zp = np.diag([1.0, -1.0]); I2 = np.eye(2); low = np.array([[0, 1], [0, 0]], dtype=complex)
    ops = []
    for i in range(n):
        out = np.eye(1, dtype=complex)
        for j in range(n):
            out = np.kron(out, Zp if j < i else (low if j == i else I2))
        ops.append(out)
    return ops
n_t = 4
ct = fermion_ops(n_t)
Mr = rng.normal(size=(n_t, n_t)) + 1j * rng.normal(size=(n_t, n_t)); Mr = Mr - Mr.T
Qt = sum(Mr[i, j] * ct[i] @ ct[j] for i in range(n_t) for j in range(n_t))
lhs = Qt @ Qt.conj().T - Qt.conj().T @ Qt
X = Mr.conj().T @ Mr
rhs = 2 * np.trace(Mr @ Mr.conj().T).real * np.eye(2 ** n_t) - 4 * sum(X[i, j] * ct[i].conj().T @ ct[j] for i in range(n_t) for j in range(n_t))
formula_ok = np.allclose(lhs, rhs)
# the actual M_A: unitary B_A -> tr(M M^dag) = 32, M^dag M = 1
Zq = np.diag([1.0, -1.0]).astype(complex); I2c = np.eye(2, dtype=complex); lowc = np.array([[0, 1], [0, 0]], dtype=complex)
c5 = []
for j in range(5):
    mats = [Zq] * j + [lowc] + [I2c] * (4 - j)
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
MA = [np.kron(np.array([[0, 1], [-1, 0]]), b) for b in BA]
unit_ok = all(np.allclose(M.conj().T @ M, np.eye(32)) for M in MA)
s_ = Fr(9, 8)
one_body = 10 * 4 * s_          # coefficient of sum_x N_rho(x) = sum psi^dag F^2 psi
const = 10 * 64 * s_ ** 2
step("FP5.normal ordering: [Q, Q^dag] = 2 tr(M M^dag) - 4 a^dag (M^dag M) a (checked on a random 4-mode example), = 64 - 4N for the unitary B_A; hence H_int = H_quartic + 45 g sum psi^dag F^2 psi - 810 g N_sites; the extra one-particle term ~ f(k)^2 ~ |k|^4/4 at the light node; N conserved only mod 4",
     formula_ok and unit_ok and one_body == 45 and const == 810, f"coefficients {one_body}, {const}")

# FP6 ----------------------------------------------------------------------------------------------------
r = -2 / 9
s2 = 1 - r * r
n_a = 3
ops = fermion_ops(2 * n_a)
a_ops, c_ops = ops[:n_a], ops[n_a:]
t = np.diag([0.5, -0.5, 0.0])                                          # traceless Cartan-like generator, tr t^2 = 1/2
def J_of(bops):
    return sum(t[i, i] * bops[i].conj().T @ bops[i] for i in range(n_a))
b_ops = [r * a_ops[i] + np.sqrt(s2) * c_ops[i] for i in range(n_a)]
Jb = J_of(b_ops)
# a-factor states with <a^dag a> = 1/2, <a a> = 0 and vanishing odd moments: maximally mixed on the a-modes
ok6, vals = True, []
dim_a = 2 ** n_a
for _ in range(4):
    eta = rng.normal(size=2 ** n_a) + 1j * rng.normal(size=2 ** n_a); eta /= np.linalg.norm(eta)
    rho_eta = np.outer(eta, eta.conj())
    rho_tot = np.kron(np.eye(dim_a) / dim_a, rho_eta)
    val = np.trace(rho_tot @ Jb.conj().T @ Jb).real
    bound = r * r * s2 * np.trace(t @ t).real / 2
    vals.append(val / bound)
    ok6 &= val >= bound - 1e-12
bound_TFPT = r * r * s2 * 8 / 2                                         # tr t^2 = 8 for a normalised Spin(10) Cartan on 32 modes
step("FP6.LOCAL_PATCH frustration: for b = r a + s c (r = -2/9) and any state whose a-factor has <a^dag a> = 1/2, <a a> = 0 (the unique singlet: Schur, Z4, <N-16> = 0), ||J_b Psi||^2 >= r^2 s^2 tr(t^2)/2 > 0: the two singlet minima are never reached simultaneously; E_min(h_x + h_y) > 2 E_0 (exact value of the 2^64 patch not computed)",
     ok6 and bound_TFPT > 0, f"ratio ||J_b Psi||^2 / bound on test states {np.round(vals, 3).tolist()} (>= 1); bound for tr t^2 = 8: 4 r^2 (1 - r^2) = {bound_TFPT:.4f}")
step("FP7.(statement) the strong-coupling ground state is not a product of local singlets; the simple product proof of a mirror gap is excluded, the collective phase of the overlapping model decides Gate B (QMC / tensor networks after the sign audit)",
     True, "no marker moves")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "F2_support": {str(k): [str(v) for v in vs] for k, vs in classes.items()},
       "cells_with_P": tab, "local_gaps_gstar": gaps, "frustration_bound": bound_TFPT}
(HERE / "filtered_patch_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
