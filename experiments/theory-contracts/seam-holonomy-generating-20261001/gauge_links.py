#!/usr/bin/env python3
"""Part N: GAUGE.LINK.01 -- internal holonomies as link variables of the mark walk.

The mark walk factorises into three controlled partial shifts, W = V_1 V_2 V_3 with
V_i = T_i (x) Pi_i^+ + T_i^dag (x) Pi_i^-, Pi_i^+- = (1 +- alpha_i)/2 (alpha_i from Part K).  Replacing T_i by
covariant shifts with link variables U_i(x) on the cubic edges gives the gauged walk W_U.
  N1 factorisation (position space): V_1 V_2 V_3 equals the eight-step bcc walk sum_eps A_eps T_eps.
  N2 for random U(1) x SU(2) links W_U is unitary and exactly gauge covariant: W_{U^g} = G W_U G^dag.
  N3 a constant link phase is the flat family-Cartan twist of W_S (Part T): W_A(k) = W(k - A).
  N4 minimal coupling in the continuum: a uniform field B along z gives relativistic Landau levels
     omega_n = sqrt(2 n B) (n = 1, 2, 3) and a zero level, with lattice corrections that shrink as B -> 0.
  N5 Nielsen-Ninomiya in gauged form: in the field, one chiral block carries chiral n = 0 branches near
     k_z = 0 (light Weyl node) and near k_z = pi/2 (the doubler P) with OPPOSITE slopes d omega/d k_z = +-1: the
     chiral anomaly of the light node is cancelled by its doubler.  A link acting on one chirality block alone
     therefore couples vectorially (light node + doubler).
  N6 count: TFPT matter per family is Lambda^even C^5 = 16 (1 + 10 + 5 of SU(5), S48 = 3 x 16); the Z16
     index per family is 16 = 0 mod 16, the necessary condition for gapping a mirror/doubler copy without
     breaking the chiral symmetry (symmetric mass generation, cited) -- condition only, not a mechanism.
  N7 (statement) gauge fields enter as link variables on the partial shifts with exact covariance and the
     correct continuum coupling; the chiral gauge structure of the Standard Model is NOT obtained at the free
     level (doubler); it needs an interacting, symmetric gap for the doubler sector (CHIRAL4D.NOMIRROR.01).
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/gauge_links.py
"""
from __future__ import annotations

import json
import sys
from itertools import combinations, product
from math import comb
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
Gam = -g[0] @ g[1] @ g[2] @ g[3]
Pp_, Pm_ = [(I4 + s * alpha[i]) / 2 for i in range(3) for s in (1,)], [(I4 - alpha[i]) / 2 for i in range(3)]
signs = [np.array(s) for s in product((-1, 1), repeat=3)]
A = {tuple(s): (I4 + s[0] * alpha[0]) / 2 @ (I4 + s[1] * alpha[1]) / 2 @ (I4 + s[2] * alpha[2]) / 2 for s in signs}


def W(k):
    return (np.cos(k[0]) * I4 - 1j * np.sin(k[0]) * alpha[0]) @ (np.cos(k[1]) * I4 - 1j * np.sin(k[1]) * alpha[1]) @ (np.cos(k[2]) * I4 - 1j * np.sin(k[2]) * alpha[2])


# position-space machinery on an L^3 torus with an n-dim internal (gauge) space -----------------------
def walk_from_links(L, links, n):
    """W_U = V_1 V_2 V_3; links[i][x] = U_i(x) on the edge x -> x + e_i (n x n unitary)."""
    sites = list(product(range(L), repeat=3))
    idx = {x: j for j, x in enumerate(sites)}
    d = 4 * n
    def V(i):
        M = np.zeros((d * len(sites), d * len(sites)), dtype=complex)
        for x in sites:
            xp = list(x); xp[i] = (xp[i] + 1) % L; xp = tuple(xp)
            # (T_i psi)(x+e_i) = U_i(x) psi(x) ; (T_i^dag psi)(x) = U_i(x)^dag psi(x+e_i)
            r, c = idx[xp], idx[x]
            M[d * r:d * r + d, d * c:d * c + d] += np.kron(Pp_[i], links[i][x])
            M[d * c:d * c + d, d * r:d * r + d] += np.kron(Pm_[i], links[i][x].conj().T)
        return M
    return V(0) @ V(1) @ V(2), sites, idx


L, n = 4, 1
ones = [{x: np.eye(1, dtype=complex) for x in product(range(L), repeat=3)} for _ in range(3)]
W0, sites, idx = walk_from_links(L, ones, 1)
Wb = np.zeros_like(W0)
for x in sites:
    for s_, As in A.items():
        y = tuple((x[c] + s_[c]) % L for c in range(3))                 # amplitude A_eps moves x -> x + eps
        Wb[4 * idx[y]:4 * idx[y] + 4, 4 * idx[x]:4 * idx[x] + 4] += As
step("N1.factorisation: V_1 V_2 V_3 (controlled partial shifts along the cubic axes) equals the eight-step bcc walk sum_eps A_eps T_eps in position space",
     np.allclose(W0, Wb), f"max deviation {np.abs(W0 - Wb).max():.1e}")

rng = np.random.default_rng(17)
def rand_u2():
    q, r = np.linalg.qr(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))
    return q * np.exp(1j * rng.uniform(0, 2 * np.pi))
L2 = 3
links = [{x: rand_u2() for x in product(range(L2), repeat=3)} for _ in range(3)]
WU, sites2, idx2 = walk_from_links(L2, links, 2)
gx = {x: rand_u2() for x in sites2}
links_g = [{x: gx[tuple((x[c] + (c == i)) % L2 for c in range(3))] @ links[i][x] @ gx[x].conj().T for x in sites2} for i in range(3)]
WUg, _, _ = walk_from_links(L2, links_g, 2)
G = np.zeros_like(WU)
for x in sites2:
    G[8 * idx2[x]:8 * idx2[x] + 8, 8 * idx2[x]:8 * idx2[x] + 8] = np.kron(I4, gx[x])
unitary = np.allclose(WU.conj().T @ WU, np.eye(WU.shape[0]))
covariant = np.allclose(WUg, G @ WU @ G.conj().T)
step("N2.random U(2) = U(1) x SU(2) links: the gauged walk W_U is unitary and exactly gauge covariant, W_{U^g} = G W_U G^dag",
     unitary and covariant, f"unitarity {np.abs(WU.conj().T @ WU - np.eye(WU.shape[0])).max():.1e}, covariance {np.abs(WUg - G @ WU @ G.conj().T).max():.1e}")

Aconst = np.array([0.31, -0.17, 0.05])
phase_links = [{x: np.array([[np.exp(-1j * Aconst[i])]]) for x in product(range(L), repeat=3)} for i in range(3)]
WA, _, _ = walk_from_links(L, phase_links, 1)
kk = 2 * np.pi * np.array([1, 2, 3]) / L
psi_k = np.concatenate([np.exp(1j * kk @ np.array(x)) * np.eye(4) for x in sites], axis=0) / np.sqrt(len(sites))
red = psi_k.conj().T @ WA @ psi_k
flat = np.allclose(red, W(kk - Aconst)) or np.allclose(red, W(kk + Aconst))
step("N3.a constant link phase is the flat family-Cartan twist (Part T): the gauged walk on a plane wave is W(k -+ A)",
     flat, "flat background = minimal shift of the node")


# Landau problem: chain in x (period Lx), plane waves in y, z, Landau gauge A_y = B x ----------------------
def landau_walk(Lx, p, ky, kz, P=None):
    B = 2 * np.pi * p / Lx
    d = 4 if P is None else P.shape[1]
    proj = (lambda M: M) if P is None else (lambda M: P.conj().T @ M @ P)
    Pi = [(proj((I4 + alpha[i]) / 2), proj((I4 - alpha[i]) / 2)) for i in range(3)]
    V1 = np.zeros((d * Lx, d * Lx), dtype=complex)
    V2 = np.zeros_like(V1); V3 = np.zeros_like(V1)
    for x in range(Lx):
        xp = (x + 1) % Lx
        V1[d * xp:d * xp + d, d * x:d * x + d] += Pi[0][0]
        V1[d * x:d * x + d, d * xp:d * xp + d] += Pi[0][1]
        k2 = ky - B * x
        V2[d * x:d * x + d, d * x:d * x + d] = np.exp(-1j * k2) * Pi[1][0] + np.exp(1j * k2) * Pi[1][1]
        V3[d * x:d * x + d, d * x:d * x + d] = np.exp(-1j * kz) * Pi[2][0] + np.exp(1j * kz) * Pi[2][1]
    return V1 @ V2 @ V3, B


def levels(Lx, p):
    U, B = landau_walk(Lx, p, 0.0, 0.0)
    om = np.sort(np.abs(np.angle(np.linalg.eigvals(U))))
    return om, B
dev, info = [], []
for Lx, p in ((480, 4), (960, 4)):
    om, B = levels(Lx, p)
    zero = np.sum(om < 0.2 * np.sqrt(2 * B))
    pos = om[om > 0.5 * np.sqrt(2 * B)]
    found = []
    for nn in (1, 2, 3):
        target = np.sqrt(2 * nn * B)
        found.append(pos[np.argmin(np.abs(pos - target))] / target - 1)
    dev.append(max(abs(x) for x in found))
    info.append(f"B = {B:.4f}: zero-level states {zero}, omega_n/sqrt(2nB) - 1 = {np.round(found, 4).tolist()}")
step("N4.minimal coupling: a uniform B along z gives relativistic Landau levels omega_n = sqrt(2nB) (n = 1,2,3) plus a zero level; the lattice deviation shrinks as B -> 0",
     dev[0] < 0.05 and dev[1] < dev[0] and dev[1] < 0.03, "; ".join(info))

ev_, V_ = np.linalg.eigh(Gam)
Pblk = V_[:, ev_ > 0]
Lx, p = 240, 3
B = 2 * np.pi * p / Lx
def zero_branch_slopes(kz0, h=2e-3):
    sl = []
    for kz in (kz0,):
        Ua, _ = landau_walk(Lx, p, 0.0, kz - h, Pblk)
        Ub, _ = landau_walk(Lx, p, 0.0, kz + h, Pblk)
        oa = np.angle(np.linalg.eigvals(Ua)); ob = np.angle(np.linalg.eigvals(Ub))
        oa = np.sort(oa[np.abs(oa) < 0.3 * np.sqrt(2 * B)]); ob = np.sort(ob[np.abs(ob) < 0.3 * np.sqrt(2 * B)])
        if len(oa) == len(ob) and len(oa) > 0:
            sl = (ob - oa) / (2 * h)
    return np.array(sl)
s_light = zero_branch_slopes(0.0)
s_doub = zero_branch_slopes(np.pi / 2)
step("N5.Nielsen-Ninomiya in gauged form: in the field, one chirality block has chiral n = 0 branches near k_z = 0 (light node) and near k_z = pi/2 (doubler P) with OPPOSITE slopes +-1: the doubler cancels the chiral anomaly, a block-only link couples vectorially",
     len(s_light) > 0 and len(s_doub) == len(s_light) and np.allclose(np.abs(s_light), 1, atol=0.05) and np.allclose(np.abs(s_doub), 1, atol=0.05)
     and np.all(np.sign(s_light) == np.sign(s_light[0])) and np.all(np.sign(s_doub) == -np.sign(s_light[0])),
     f"{len(s_light)} light zero-branch states with slopes {np.round(np.unique(np.round(s_light, 3)), 3).tolist()}, {len(s_doub)} doubler states with slopes {np.round(np.unique(np.round(s_doub, 3)), 3).tolist()}")

per_family = sum(comb(5, r) for r in (0, 2, 4))
su5 = [comb(5, 0), comb(5, 2), comb(5, 4)]
step("N6.count: matter per family = Lambda^even C^5 = 16 (1 + 10 + 5 of SU(5)), S48 = 3 x 16; Z16 index per family 16 = 0 mod 16: the necessary anomaly condition for a symmetric (chirality-preserving) gap of a doubler copy holds -- condition only (mechanism cited, not tested)",
     per_family == 16 and su5 == [1, 10, 5] and 3 * per_family == 48 and per_family % 16 == 0, "symmetric mass generation: Wang-Wen 2018, Razamat-Tong 2021 (cited)")
step("N7.(statement) gauge fields = link variables on the partial shifts: exact covariance, correct continuum coupling (Landau levels); the chiral Standard-Model gauge structure is NOT obtained at the free level (doubler cancels the anomaly); it needs an interacting symmetric gap of the doubler sector (CHIRAL4D.NOMIRROR.01)",
     True, "no chirality or TOE marker moves")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "landau": info,
       "slopes_light": s_light.tolist(), "slopes_doubler": s_doub.tolist()}
(HERE / "gauge_links_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
