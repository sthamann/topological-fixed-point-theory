#!/usr/bin/env python3
"""THETA1.UNIVERSALITY.CLASS.01 -- preregistered universality test of the Gate-B phase over (b, theta, gamma); exact prechecks only.

Question: if the frozen point P0 = (b, theta, gamma) = (0, 120 deg, 0) passes Gate B, does the symmetric-gapped chiral phase persist
in an open neighbourhood, so that P1/P2 need only select a region R_chiral instead of one microscopic formula?
Nothing is run with QMC here: Gate B for P0 is undecided and the run condition below is not met.  This file freezes the design
and computes the exact free / local / IR invariants of the five design points.
  UC1 design points P0, P_b+- = (+-0.25, 120, 0), P_th+- = (0, 100, 0), (0, 140, 0): all minimally doubled (Gamma +1, P -1).
  UC2 the mirror stays at P for every point, so the existing P-compatible cells (side divisible by 4) remain valid.
  UC3 the filter f_b: f(Gamma) = 0 with Hessian I/2 for all b; f(+P) = 1 + b, f(-P) = 1 - b; PSD for |b| <= 1.
  UC4 canonical normalisation: diag(F F^dag) = 9/8 + b^2/8, so g* = g (9/8 + b^2/8)^2 (b = 0: 81/64 g): the exact one-site
      so(11) spectrum is unchanged up to this rescaling -> the local necessary condition holds at every point.
  UC5 the mirror sees the filtered quartic with weight f(P)^4 = (1 + b)^4: b = +-0.25 changes the effective mirror coupling by
      2.44 / 0.32 at fixed g; the analysis must report g_mirror = g (1 + b)^4 next to g.
  UC6 IR identity: every point has the same unit isotropic Gamma cone (|E - |k|| = O(k^2)) and f(Gamma) = 0 quadratically, so the IR
      protection Sigma = f^2 Sigma^ holds unchanged; differences are confined to the mirror sector.
  UC7 free reference: lowest positive single-particle energy on the primitive L = 4 and L = 8 cells for each point (pins for the
      free-limit validation any future QMC run must reproduce first).
  UC8 (statement) preregistration frozen (see contract_index.json): status PREREGISTERED_NOT_RUN.
  UC9 amendment 1 (before any data): arm A fixed g, arm B matched mirror coupling g(b) = g0/(1+b)^4; Gamma and mirror Weyl
      speeds are exactly 1 at every design point (V_Gamma = I, V_P orthogonal), so no velocity renormalisation is needed.
    python3 -B experiments/theory-contracts/theta1-universality-class-20261005/prechecks.py
"""
from __future__ import annotations

import sys
from itertools import product
from pathlib import Path

import numpy as np
from scipy.optimize import fsolve

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "minimal-runtime-reconstruction-20261005"))
from ur_lib import Log  # noqa: E402

log = Log()
sig = [np.array([[0, 1], [1, 0]], dtype=complex), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]).astype(complex)]
B = np.array([[1, 1, -1], [1, -1, 1], [-1, 1, 1]], dtype=float); Binv = np.linalg.inv(B)
Cp = np.roll(np.eye(3), 1, axis=1)
circ = lambda c: c[0] * np.eye(3) + c[1] * Cp + c[2] * Cp @ Cp
DESIGN = {"P0": (0.0, 120.0, 0.0), "P_b+": (0.25, 120.0, 0.0), "P_b-": (-0.25, 120.0, 0.0),
          "P_th-": (0.0, 100.0, 0.0), "P_th+": (0.0, 140.0, 0.0)}


def member(th_deg, gam):
    th = np.radians(th_deg); s_b = 1 - 3 * gam
    beta = circ([(s_b + 2 * np.cos(th)) / 3, (s_b - np.cos(th) + np.sqrt(3) * np.sin(th)) / 3,
                 (s_b - np.cos(th) - np.sqrt(3) * np.sin(th)) / 3])
    return dict(beta=beta, g1=-gam, g2=gam, alpha=np.eye(3) - gam * np.ones((3, 3)))


def d_p(p, beta, g1, g2, alpha):
    S = p.sum()
    return alpha @ np.sin(p) + beta @ (np.cos(p) - 1) + (g2 * np.sin(S) + g1 * (np.cos(S) - 1)) * np.ones(3)


def jac_p(p, beta, g1, g2, alpha):
    S = p.sum()
    return alpha * np.cos(p)[None, :] - beta * np.sin(p)[None, :] + (g2 * np.cos(S) - g1 * np.sin(S)) * np.ones((3, 3))


def census(m, n=30):
    grid = np.linspace(0, 2 * np.pi, n, endpoint=False)
    G = np.array(np.meshgrid(grid, grid, grid, indexing="ij")).reshape(3, -1).T
    val = np.array([np.linalg.norm(d_p(p, **m)) for p in G])
    roots = []
    for p0 in G[np.argsort(val)[:400]]:
        p, _, ier, _ = fsolve(lambda q: d_p(q, **m), p0, fprime=lambda q: jac_p(q, **m), full_output=True, xtol=1e-13)
        if ier != 1 or np.linalg.norm(d_p(p, **m)) > 1e-10:
            continue
        p = np.mod(p, 2 * np.pi)
        if not any(np.linalg.norm(np.angle(np.exp(1j * (p - q)))) < 1e-6 for q, _ in roots):
            roots.append((p, int(np.sign(np.linalg.det(jac_p(p, **m))))))
    return roots


Pp = np.pi / 2 * np.ones(3)
summary = {}
ok1 = ok2 = True
for name, (b, th, gam) in DESIGN.items():
    m = member(th, gam)
    r = census(m)
    labels = sorted(("Gamma" if np.allclose(np.angle(np.exp(1j * p)), 0, atol=1e-6) else
                     "P" if np.allclose(np.angle(np.exp(1j * (p - Pp))), 0, atol=1e-6) else "other", c) for p, c in r)
    ok1 &= labels == [("Gamma", 1), ("P", -1)]
    ok2 &= np.allclose(d_p(Pp, **m), 0)
    summary[name] = {"b": b, "theta": th, "gamma": gam, "nodes": labels}
log.step("UC1 all five design points are minimally doubled: exactly Gamma (+1) and P (-1)", ok1,
         "; ".join(f"{k}: {v['nodes']}" for k, v in summary.items()))
cells_ok = all(np.allclose(np.mod(Pp * L / (2 * np.pi), 1), 0) for L in (4, 8)) and not np.allclose(np.mod(Pp * 6 / (2 * np.pi), 1), 0)
log.step("UC2 the mirror stays at P = (pi/2)(1,1,1) for every point (p-coordinates); primitive cells with side divisible by 4 contain P, "
         "so the existing P-compatible cells remain valid", ok2 and cells_ok)

ccc = lambda k: np.prod(np.cos(k)); sss = lambda k: np.prod(np.sin(k))
f_b = lambda k, b: 1 - ccc(k) + b * sss(k)
kk = 1e-3 * np.array([0.3, -0.7, 0.5])
ok3 = True
for name, (b, th, gam) in DESIGN.items():
    hess_ok = abs(f_b(kk, b) - 0.5 * kk @ kk) < 1e-8
    ok3 &= hess_ok and np.isclose(f_b(Pp, b), 1 + b) and np.isclose(f_b(-Pp, b), 1 - b) and abs(f_b(np.zeros(3), b)) < 1e-15
    summary[name].update({"f_P": 1 + b, "f_minusP": 1 - b})
log.step("UC3 filter f_b: f(Gamma) = 0 with Hessian I/2 for every point; f(+P) = 1 + b, f(-P) = 1 - b (PSD for |b| <= 1, "
         "SOURCE.FILTER.UNIQUENESS.01)", ok3)

steps = [np.array(e) for e in product((-1, 1), repeat=3)]
ok4, g_ratio = True, {}
for name, (b, th, gam) in DESIGN.items():
    eta = -b / 8
    c = {tuple(e): (-1 / 8 + 1j * eta if np.prod(e) == -1 else -1 / 8 - 1j * eta) for e in steps}
    diag = 1 + sum(abs(v) ** 2 for v in c.values())
    ok4 &= np.isclose(diag, 9 / 8 + b ** 2 / 8)
    sym_ok = all(abs(sum(cv * np.exp(1j * np.array(e) @ k) for e, cv in c.items()) + 1 - f_b(k, b)) < 1e-12
                 for k in np.random.default_rng(1).normal(size=(5, 3)))
    ok4 &= sym_ok
    g_ratio[name] = diag ** 2
    summary[name].update({"diag_FFdag": diag, "g_star_over_g": diag ** 2})
log.step("UC4 canonical normalisation diag(F F^dag) = 9/8 + b^2/8: g* = g (9/8 + b^2/8)^2 (81/64 at b = 0); the exact one-site so(11) "
         "spectrum (ground state unique, gaps 144/232 g*) holds at every point up to this rescaling",
         ok4 and np.isclose(g_ratio["P0"], 81 / 64), "; ".join(f"{k}: g*/g = {v:.4f}" for k, v in g_ratio.items()))

mirror = {name: (1 + b) ** 4 for name, (b, th, gam) in DESIGN.items()}
for name in DESIGN:
    summary[name]["mirror_weight_(1+b)^4"] = mirror[name]
log.step("UC5 the mirror sees the filtered quartic with weight (1 + b)^4: P_b+ 2.441, P_b- 0.316 at fixed g -- report g_mirror = g (1 + b)^4 "
         "next to g; theta changes only the mirror frame", np.isclose(mirror["P_b+"], 1.25 ** 4) and np.isclose(mirror["P_b-"], 0.75 ** 4),
         "; ".join(f"{k}: {v:.3f}" for k, v in mirror.items()))


def h_k(k, m):
    d = Binv @ d_p(B @ k, **m)
    return sum(di * s for di, s in zip(d, sig))


ok6, spread = True, []
k_ir = 1e-3 * np.array([0.2, 0.5, -0.4])
for name, (b, th, gam) in DESIGN.items():
    e = np.linalg.eigvalsh(h_k(k_ir, member(th, gam)))
    spread.append(abs(e[1] - np.linalg.norm(k_ir)))
    VG = Binv @ jac_p(np.zeros(3), **member(th, gam)) @ B
    ok6 &= np.allclose(VG, np.eye(3))
log.step("UC6 IR identity: unit isotropic Gamma cone at every point (V_Gamma = I, |E - |k|| = O(k^2)) and f(Gamma) = 0 quadratically -> "
         "IR protection Sigma = f^2 Sigma^ unchanged; differences are confined to the mirror sector", ok6 and max(spread) < 1e-5,
         f"max |E - |k|| at |k| = {np.linalg.norm(k_ir):.1e}: {max(spread):.1e}")


def lowest_positive(m, L):
    best = np.inf
    for n in product(range(L), repeat=3):
        p = 2 * np.pi * np.array(n) / L
        if np.allclose(np.mod(p, 2 * np.pi), 0):
            continue                                   # Gamma: decoupled zero modes (Z = 2^32 Z_active)
        e = np.linalg.norm(d_p(p, **m))
        if e > 1e-12:
            best = min(best, e)
    return best


ref = {}
for name, (b, th, gam) in DESIGN.items():
    ref[name] = {f"L{L}": round(float(lowest_positive(member(th, gam), L)), 6) for L in (4, 8)}
    summary[name]["free_lowest_positive"] = ref[name]
log.step("UC7 free reference energies (primitive L = 4, 8 cells, Gamma excluded) recorded as pins for the free-limit validation of any "
         "future run", all(np.isfinite(v) for r in ref.values() for v in r.values()),
         "; ".join(f"{k}: {v}" for k, v in ref.items()))
armB = {name: 1.0 / (1 + b) ** 4 for name, (b, th, gam) in DESIGN.items()}
vel = {}
for name, (b, th, gam) in DESIGN.items():
    m = member(th, gam)
    VG = Binv @ jac_p(np.zeros(3), **m) @ B
    VP = Binv @ jac_p(Pp, **m) @ B
    vel[name] = (float(np.linalg.svd(VG)[1].max()), float(np.linalg.svd(VG)[1].min()), float(np.linalg.svd(VP)[1].max()),
                 float(np.linalg.svd(VP)[1].min()))
    summary[name].update({"armB_g_over_g0": armB[name], "velocity_Gamma_max_min": vel[name][:2], "velocity_P_max_min": vel[name][2:]})
log.step("UC9 amendment 1 (before any data): control arm B with matched mirror coupling g(b) = g0/(1+b)^4 (P_b+ 0.4096 g0, P_b- 3.1605 g0) "
         "next to arm A (fixed g); Weyl velocities need no renormalisation: unit isotropic speed at Gamma and at P for every design point",
         np.isclose(armB["P_b+"], 0.4096) and np.isclose(armB["P_b-"], 1 / 0.75 ** 4)
         and all(np.allclose(v, 1.0, atol=1e-9) for v in vel.values()),
         "; ".join(f"{k}: g/g0 = {armB[k]:.4f}" for k in DESIGN))
log.step("UC8 (statement) THETA1.UNIVERSALITY.CLASS.01 = PREREGISTERED_NOT_RUN: run condition (Gate B CONDITIONAL at P0) not met; no QMC "
         "launched", True)

log.dump(HERE / "prechecks_results.json", {"research_id": "THETA1.UNIVERSALITY.CLASS.01", "verdict": "PREREGISTERED_NOT_RUN",
                                           "design": summary})
