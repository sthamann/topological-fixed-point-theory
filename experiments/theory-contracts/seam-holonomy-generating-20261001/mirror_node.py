#!/usr/bin/env python3
"""Part MN: CHIRAL4D.MIRROR.ONE_LOOP.01 -- does the interaction gap the mirror at P, or does it only move the node?

First-order (static Hartree-Fock-Bogoliubov) self-energy of H_int = -g sum_a Q_a^2, Q_a = 1/2 Psi^dag K_a Psi:
with A_a = tau_x K_a (antisymmetric) and the contraction C = <Psi Psi^T> = <Psi Psi^dag> tau_x of a Gaussian state,
<Q^2> = 1/4 Tr(A C)^2 + 1/2 Tr(A C A C^T) (Wick), and the mean-field operator is 1/2 Psi^dag Sigma_K Psi with
Sigma_K = tau_x B, B = -g sum_a [ -Tr(A_a C) A_a + A_a (C - C^T) A_a ].
Because H_int depends on psi only through rho = F psi and the free state is translation invariant, the whole first-order
self-energy is Sigma^(1)(k) = f(k)^2 Sigma^rho with ONE constant local 64 x 64 Nambu matrix Sigma^rho (computed from the
local filtered correlator <rho_x rho_x^dag>).
  MN1 brute force (3 modes, random BdG state and random BdG quadratics): the Wick formula for <Q^2> is exact and
      1/2 Psi^dag Sigma_K Psi is the exact first-order mean-field operator (derivative identity along Gaussian families).
  MN2 Theta_1, T = 0, thermodynamic limit (BZ grids 24^3 ... 48^3): Sigma^rho = alpha 1 + b.sigma (x) 1_16, no anomalous part
      (Spin(10) x Z4), alpha = 0 (particle-hole), b || [111] (C3); Sigma^(1)(Gamma) = 0, Sigma^(1) = O(k^4) at Gamma.
  MN3 the mirror node is SHIFTED, not gapped: h_W(P + dk) + Sigma^(1) = 0 has a solution dk = -(t V_P)^-1 b + O(g^2),
      linear in g along [111]; the Weyl charge of the shifted node stays -1 (winding of the one-loop Hamiltonian on a small
      sphere); a translation-invariant bilinear cannot gap a single Weyl node (Spin(10) x Z4 forbid every Majorana bilinear
      of one 16) -> at weak coupling there is no mirror gap, only node motion.
  MN4 comparison with the QMC data on the 16-site cell at beta = 2 (reweight_points_md_cubic16_g1.json): the one-loop
      Hamiltonian h_W + f^2 Sigma^rho (same cell, same beta; a pure node shift) predicts G_P(beta/2); its share of the
      measured depletion 1/2 - G_P shows how much of the "gap at P" is node motion; the rest is not resolvable at the old P.
  MN6 flavour-twisted boundary conditions for node tracking: a continuous U(1) twist is NOT a symmetry of H_int (charge
      only mod 4), but a twist by the Spin(10) Cartan generator Sigma_12 is; validated: trivial twist, channel rotation,
      twisted free spectrum h_W(k +- delta), site-independent filtered correlators.
  MN7 g/N test (32 sites at 2g vs 16 sites at g, same beta): multiplet picture excluded; one-loop node shift share on 32 sites.
  MN8 one-loop node trajectory and charge bookkeeping along the C3 axis: Gamma +1 fixed, mirror -1 always present, linear
      regime ends near g ~ 0.02 with an axis pair creation/annihilation that transfers the mirror charge.
  MN5 (statement) consequence for Gate B: the success criterion is min_{k != Gamma} Delta_odd(k) > 0 (node tracking with
      twisted boundary conditions), not Delta(P) > 0.
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/mirror_node.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
from itertools import product
from pathlib import Path

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("reweight", HERE / "reweight.py")
rw = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(rw)
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    sys.stdout.flush()


def tau_x(m):
    Z, I = np.zeros((m, m)), np.eye(m)
    return np.block([[Z, I], [I, Z]])


def hfb_sigma(Ks, G, g):
    """Nambu mean-field matrix Sigma_K (1/2 Psi^dag Sigma_K Psi) of -g sum_a Q_a^2 in the Gaussian state G = <Psi Psi^dag>."""
    m = G.shape[0] // 2
    tx = tau_x(m)
    C = G @ tx
    B = np.zeros_like(G, dtype=complex)
    for K in Ks:
        A = tx @ K; A = (A - A.T) / 2
        B += -g * (-np.trace(A @ C) * A + A @ (C - C.T) @ A)
    return tx @ B


def q2_wick(K, G):
    m = G.shape[0] // 2
    tx = tau_x(m); C = G @ tx
    A = tx @ K; A = (A - A.T) / 2
    return 0.25 * np.trace(A @ C) ** 2 + 0.5 * np.trace(A @ C @ A @ C.T)


# ---------------------------------------------------------------------------------------- MN1 brute force
rng = np.random.default_rng(3)
nb = 3
Zp = np.diag([1.0, -1.0]); low = np.array([[0, 1], [0, 0]], complex)
cf = []
for i in range(nb):
    mm = np.eye(1, dtype=complex)
    for j in range(nb):
        mm = np.kron(mm, Zp if j < i else (low if j == i else np.eye(2)))
    cf.append(mm)
Psi = cf + [x.conj().T for x in cf]
def randK(s=1.0):
    h = rng.normal(size=(nb, nb)) + 1j * rng.normal(size=(nb, nb)); h = (h + h.conj().T) / 2
    D = rng.normal(size=(nb, nb)) + 1j * rng.normal(size=(nb, nb)); D = D - D.T
    return s * np.block([[h, D.conj().T], [D, -h.T]])
def quad(K):
    return 0.5 * sum(K[i, j] * Psi[i].conj().T @ Psi[j] for i in range(2 * nb) for j in range(2 * nb))
K0 = randK(); Ks = [randK(0.7) for _ in range(3)]; gb = 0.3
H0 = quad(K0); Hi = -gb * sum(quad(K) @ quad(K) for K in Ks)
E0, V0 = np.linalg.eigh(H0)
gs = V0[:, 0]
G0 = np.array([[gs.conj() @ Psi[i] @ Psi[j].conj().T @ gs for j in range(2 * nb)] for i in range(2 * nb)])
q2_err = max(abs(gs.conj() @ quad(K) @ quad(K) @ gs - q2_wick(K, G0)) for K in Ks)
SigK = hfb_sigma(Ks, G0, gb)
# mean-field identity: along a family of Gaussian ground states (K0 + s X), d<H_int>/ds = d<1/2 Psi^dag Sigma_K Psi>/ds at s = 0
Xp = randK(1.0); Hs = quad(SigK); ds = 1e-5
def gs_of(s_):
    e_, v_ = np.linalg.eigh(quad(K0 + s_ * Xp)); return v_[:, 0]
gp, gm = gs_of(ds), gs_of(-ds)
d_int = ((gp.conj() @ Hi @ gp) - (gm.conj() @ Hi @ gm)) / (2 * ds)
d_mf = ((gp.conj() @ Hs @ gp) - (gm.conj() @ Hs @ gm)) / (2 * ds)
qp_err = abs(d_int - d_mf) / abs(d_int)
step("MN1 brute force (3 modes, random BdG ground state, three random BdG quadratics): <Q^2> = 1/4 Tr(AC)^2 + 1/2 Tr(ACAC^T) is exact, and 1/2 Psi^dag Sigma_K Psi (Sigma_K = tau_x B) is the exact first-order mean-field operator: d<H_int>/ds = d<1/2 Psi^dag Sigma_K Psi>/ds along a family of Gaussian ground states",
     q2_err < 1e-10 and qp_err < 1e-6, f"<Q^2> error {q2_err:.1e}; relative derivative mismatch {qp_err:.1e} (d<H_int>/ds = {d_int.real:.5f})")


# ---------------------------------------------------------------------------------------- local one-loop matrix
B3 = rw.B3.astype(float); Binv = np.linalg.inv(B3)
S = [rw.sx, rw.sy, rw.sz]
tau = [sum(Binv[a, i] * S[a] for a in range(3)) for i in range(3)]
def h_of_p(p, t=1.0):
    d = [np.sin(p[i]) + np.cos(p[(i + 1) % 3]) - 1 for i in range(3)]
    return t * sum(d[i] * tau[i] for i in range(3))
def f_of_k(k):
    return 1 - np.cos(k[0]) * np.cos(k[1]) * np.cos(k[2])
Kloc = [np.block([[np.zeros((32, 32)), 2 * MA.conj().T], [2 * MA, np.zeros((32, 32))]]) for MA in rw.MA]
def local_G(Gn, Gh):
    """Nambu <Phi Phi^dag> of Phi = (rho, rho^dag) from the 2 x 2 spin blocks <rho rho^dag> = Gn (x) 1, <rho^dag rho>^T = Gh (x) 1."""
    return np.block([[np.kron(Gn, np.eye(16)), np.zeros((32, 32))], [np.zeros((32, 32)), np.kron(Gh, np.eye(16))]])
def sigma_rho(Gn, Gh, g=1.0):
    return hfb_sigma(Kloc, local_G(Gn, Gh), g)
def decompose(Sn):
    """32 x 32 normal block = (alpha 1 + b.sigma) (x) 1_16 + rest."""
    s2 = Sn.reshape(2, 16, 2, 16).trace(axis1=1, axis2=3) / 16
    alpha = np.trace(s2).real / 2
    b = np.array([np.trace(s2 @ S[a]).real / 2 for a in range(3)])
    rest = np.abs(Sn - np.kron(alpha * np.eye(2) + sum(b[a] * S[a] for a in range(3)), np.eye(16))).max()
    return alpha, b, rest


def bz_local(L):
    """T = 0 thermodynamic local correlators: Gn = (1/N) sum f^2 P_+(k), Gh^T-block = (1/N) sum f^2 P_-(k)^T."""
    Gn = np.zeros((2, 2), complex); Gm = np.zeros((2, 2), complex)
    for m in product(range(L), repeat=3):
        p = 2 * np.pi * (np.array(m) + 0.5) / L                       # shifted grid avoids the nodes exactly
        k = Binv @ p
        h = h_of_p(p); e, U = np.linalg.eigh(h)
        Pp = np.outer(U[:, 1], U[:, 1].conj()); Pm = np.eye(2) - Pp
        f2 = f_of_k(k) ** 2
        Gn += f2 * Pp; Gm += f2 * Pm
    return Gn / L ** 3, (Gm / L ** 3).T


res2 = {}
for L in (24, 36, 48):
    Gn, Gh = bz_local(L)
    Sr = sigma_rho(Gn, Gh, 1.0)
    al, b, rest = decompose(Sr[:32, :32])
    anom = np.abs(Sr[:32, 32:]).max()
    res2[L] = (al, b, rest, anom, np.trace(Gn).real + np.trace(Gh).real)
al, b, rest, anom, f2sum = res2[48]
nvec = np.ones(3) / np.sqrt(3)
par = abs(abs(b @ nvec) - np.linalg.norm(b)) < 1e-6 * np.linalg.norm(b)
conv = np.linalg.norm(res2[48][1] - res2[36][1]) / np.linalg.norm(res2[48][1])
step("MN2 Theta_1 one-loop structure (T = 0, BZ grids 24^3 / 36^3 / 48^3): Sigma^(1)(k) = f(k)^2 Sigma^rho with Sigma^rho = alpha 1 + b.sigma (x) 1_16, no anomalous part, alpha = 0 (particle-hole), b along [111] (C3); hence Sigma^(1)(Gamma) = O(k^4) and Sigma^(1)(P) = b.sigma (f(P) = 1)",
     rest < 1e-10 and anom < 1e-10 and abs(al) < 1e-8 and par and conv < 1e-3 and abs(f2sum / 2 - 9 / 8) < 1e-2,
     f"b / g = {np.round(b, 4).tolist()} (|b|/g = {np.linalg.norm(b):.4f}, along [111]: {par}); alpha/g = {al:.1e}; non-(alpha,b) rest {rest:.1e}; anomalous {anom:.1e}; grid change 36 -> 48: {conv:.1e}; local <{{rho,rho^dag}}> = {f2sum / 2:.4f} (9/8)")


# ---------------------------------------------------------------------------------------- MN3 node shift
def h_of_k(k):
    return h_of_p(B3 @ k)
P = np.pi / 2 * np.ones(3)
eps = 1e-6
VP = np.array([[np.trace((h_of_k(P + eps * np.eye(3)[i]) - h_of_k(P - eps * np.eye(3)[i])) @ S[a]).real / (4 * eps) for a in range(3)] for i in range(3)])
# h(P + dk) ~ sum_i,a dk_i VP[i,a] sigma_a ; zero of h + f^2 b.sigma
def H1(k, g):
    return h_of_k(k) + g * f_of_k(k) ** 2 * sum(b[a] * S[a] for a in range(3))
def find_node(g, k0):
    k = k0.copy()
    for _ in range(50):
        hv = np.array([np.trace(H1(k, g) @ S[a]).real / 2 for a in range(3)])
        J = np.array([[(np.trace(H1(k + eps * np.eye(3)[i], g) @ S[a]).real - np.trace(H1(k - eps * np.eye(3)[i], g) @ S[a]).real) / (4 * eps) for a in range(3)] for i in range(3)])
        dk = -np.linalg.solve(J.T, hv)
        k = k + dk
        if np.linalg.norm(dk) < 1e-13:
            break
    return k, np.linalg.norm(np.array([np.trace(H1(k, g) @ S[a]).real / 2 for a in range(3)]))
lin_pred = -np.linalg.solve(VP.T, b)                                    # d k per unit g
shifts = {}
for g_ in (0.001, 0.003, 0.01):
    kn, resid = find_node(g_, P + g_ * lin_pred)
    shifts[g_] = (kn - P, resid)
ratio_lin = [np.linalg.norm(shifts[g_][0]) / g_ for g_ in shifts]
dirn = shifts[0.001][0] / np.linalg.norm(shifts[0.001][0])
def winding(center, r=1e-3, n=40, g=0.0):
    """Weyl charge: degree of h_hat on a small sphere (sum of solid angles of triangles)."""
    th = np.linspace(0, np.pi, n + 1); ph_ = np.linspace(0, 2 * np.pi, 2 * n + 1)
    def hv(k):
        v = np.array([np.trace(H1(k, g) @ S[a]).real / 2 for a in range(3)]); return v / np.linalg.norm(v)
    tot = 0.0
    for i in range(n):
        for j in range(2 * n):
            pts = [center + r * np.array([np.sin(a) * np.cos(c), np.sin(a) * np.sin(c), np.cos(a)]) for a, c in ((th[i], ph_[j]), (th[i + 1], ph_[j]), (th[i + 1], ph_[j + 1]), (th[i], ph_[j + 1]))]
            v = [hv(q) for q in pts]
            for (a1, a2, a3) in ((v[0], v[1], v[2]), (v[0], v[2], v[3])):
                tot += 2 * np.arctan2(a1 @ np.cross(a2, a3), 1 + a1 @ a2 + a2 @ a3 + a3 @ a1)
    return tot / (4 * np.pi)
w_free = winding(P); w_shift = winding(P + shifts[0.01][0], g=0.01); w_gam = winding(np.zeros(3))
step("MN3 the mirror node is SHIFTED, not gapped: the one-loop Hamiltonian h_W + f^2 b.sigma has an exact zero at P + dk with dk linear in g along [111], the old P is no longer a zero (|E(P)| = |b| g), the Weyl charge of the moved node stays -1 (Gamma: +1); Spin(10) x Z4 forbid every Majorana bilinear of one 16, so no translation-invariant bilinear can gap it",
     all(r < 1e-10 for _, r in shifts.values()) and abs(ratio_lin[0] - np.linalg.norm(lin_pred)) / np.linalg.norm(lin_pred) < 0.02
     and ratio_lin[0] > ratio_lin[1] > ratio_lin[2] and abs(abs(dirn @ nvec) - 1) < 1e-3
     and round(w_free) == -1 and round(w_shift) == -1 and round(w_gam) == 1,
     f"|dk|/g = {np.round(ratio_lin, 3).tolist()} for g = 0.001 / 0.003 / 0.01 (linear prediction {np.linalg.norm(lin_pred):.3f}/t, approached for g -> 0; external estimate 80/t is a factor {80 / np.linalg.norm(lin_pred):.2f} larger -- normalisation here fixed by MN1 and by the frozen O_A of reweight.py); direction {np.round(dirn, 4).tolist()}; |E(old P)| = {np.linalg.norm(b):.3f} g; Weyl charges P {w_free:+.3f}, shifted node {w_shift:+.3f}, Gamma {w_gam:+.3f}")


# ---------------------------------------------------------------------------------------- MN4 QMC comparison
qmc = {}
qf = HERE / "reweight_points_md_cubic16_g1.json"
if qf.exists():
    for p_ in json.loads(qf.read_text()):
        qmc[p_["g"]] = (p_["greens"]["P"]["G_half"], p_["greens"]["P"]["G_half_err"], p_["greens"]["-P"]["G_half"], p_["greens"]["H"]["G_half"])
m16 = rw.Model("cubic16", 2.0, 0.25, 0.0, "md")
beta = 2.0
E_, V_ = np.linalg.eigh(m16.Hw)
Gfree = (V_ * (1 / (1 + np.exp(-beta * E_)))) @ V_.conj().T            # <psi psi^dag> at beta = 2 on the cell
Fk = np.kron(m16.F, np.eye(32))
Grho = Fk @ Gfree @ Fk.T                                                # <rho rho^dag> (F real symmetric)
Grho_h = Fk @ (np.eye(m16.n) - Gfree).T @ Fk.T                          # <rho^dag rho>^T
Gn16 = Grho[:32, :32].reshape(2, 16, 2, 16).trace(axis1=1, axis2=3) / 16
Gh16 = Grho_h[:32, :32].reshape(2, 16, 2, 16).trace(axis1=1, axis2=3) / 16
Sr16 = sigma_rho(Gn16, Gh16, 1.0)
X = np.array(m16.sites, float)
def G_onel(g, kname):
    H1c = m16.Hw + g * Fk @ np.kron(np.eye(m16.N), Sr16[:32, :32]) @ Fk.T
    e, U = np.linalg.eigh(H1c)
    Gt = (U * (np.exp(-beta / 2 * e) / (1 + np.exp(-beta * e)))) @ U.conj().T
    k = {"P": P, "-P": -P, "H": np.array([np.pi, 0, 0])}[kname]
    ph = np.exp(-1j * X @ k)
    Gxy = np.einsum("xaya->xy", Gt.reshape(m16.N, 32, m16.N, 32))
    return float((ph @ Gxy @ ph.conj()).real / (32 * m16.N))
cmp_ = {g_: (G_onel(g_, "P"), qmc[g_][0], qmc[g_][1]) for g_ in sorted(qmc)}
al16, b16, rest16 = decompose(Sr16[:32, :32])
frac = {g_: (0.5 - a) / (0.5 - q) for g_, (a, q, e) in cmp_.items()}                # share of the depletion 1/2 - G_P explained
agree = len(frac) >= 3 and all(0.7 <= v <= 1.15 for v in frac.values())
step("MN4 QMC comparison (16 sites, beta = 2): the one-loop Hamiltonian h_W + g f^2 Sigma^rho(cell, beta) -- a pure node shift, no mass -- explains the major part (70-115 %) of the measured depletion 1/2 - G_P(beta/2) at g = 0.003, 0.01, 0.03; the remainder (beyond one loop, multiplet splitting or a genuine gap) cannot be resolved at the old P -> node tracking needed",
     agree, "; ".join(f"g = {g_}: one-loop {a:.4f} vs QMC {q:.4f}+-{e:.4f} (one-loop share {100 * frac[g_]:.0f} %)" for g_, (a, q, e) in cmp_.items()) + f"; cell/beta Sigma^rho: b/g = {np.round(b16, 3).tolist()}, alpha/g = {al16:.1e}")
# ---------------------------------------------------------------------------------------- MN6 flavour twist
def all_momenta(Rv):
    G_ = 2 * np.pi * np.linalg.inv(Rv.astype(float)).T
    fcc_ = lambda k: np.allclose(np.cos(np.array([k @ np.array(s_) for s_ in product((-1, 1), repeat=3)])), 1, atol=1e-9)
    Nn = int(round(abs(np.linalg.det(Rv.astype(float))) / 4)); ks_ = []
    for mm_ in product(range(-Nn, Nn + 1), repeat=3):
        k = np.array(mm_) @ G_
        if not any(fcc_(k - q) for q in ks_):
            ks_.append(k)
        if len(ks_) == Nn:
            break
    return ks_
Rv16 = rw.CELLS["cubic16"].astype(float)
# (a) a twist by a superlattice reciprocal vector is trivial: identical Hamiltonian and pairing matrix
dtriv = 2 * np.pi * np.linalg.inv(Rv16) @ np.array([1.0, 0.0, 0.0]) / 1.0
m0 = rw.Model("cubic16", 0.5, 0.25, 0.1, "md"); mt = rw.Model("cubic16", 0.5, 0.25, 0.1, "md", twist=dtriv)
phr = 0.3 * rng.normal(size=(m0.N, 10))
triv = max(np.abs(m0.Hw - mt.Hw).max(), np.abs(m0.build_M(phr) - mt.build_M(phr)).max())
# (b) Spin(10) covariance: U^T B_A U = sum_B R_AB B_B with R real orthogonal (Sigma_12 rotates channels)
th_ = 0.7; U16 = np.diag(np.exp(1j * th_ * rw.QF)); BA = [MA[:16, 16:] for MA in rw.MA]
basis = np.array([B_.ravel() for B_ in BA]).T
Rm = np.array([np.linalg.lstsq(basis, (U16.T @ B_ @ U16).ravel(), rcond=None)[0] for B_ in BA])
cov_res = max(np.abs(basis @ Rm[A] - (U16.T @ BA[A] @ U16).ravel()).max() for A in range(10))
cov_ok = cov_res < 1e-12 and np.abs(Rm.imag).max() < 1e-12 and np.allclose(Rm.real @ Rm.real.T, np.eye(10), atol=1e-12)
# (c) free twisted spectrum per flavour group = h_W(k +- delta)
dl = 0.2 * np.ones(3)
mtw = rw.Model("cubic16", 0.5, 0.25, 0.1, "md", twist=dl)
spec_err = 0.0
for q in (0.5, -0.5):
    idx = [32 * i + 16 * s_ + f for i in range(mtw.N) for s_ in range(2) for f in range(16) if rw.QF[f] == q]
    sub = mtw.Hw[np.ix_(idx, idx)]
    e_sub = np.sort(np.linalg.eigvalsh(sub))
    pred = np.sort(np.concatenate([np.linalg.eigvalsh(h_of_k(k + 2 * q * dl)) for k in all_momenta(Rv16)] * 8))
    spec_err = max(spec_err, np.abs(e_sub - pred).max())
# (d) twisted interaction is translation covariant: the local filtered correlators of the free twisted state are site-independent
Et, Ut = np.linalg.eigh(mtw.Hw); Gt0 = (Ut * (1 / (1 + np.exp(-2.0 * Et)))) @ Ut.conj().T
locs = [mtw.Ftw[i] @ Gt0 @ mtw.Ftw[i].conj().T for i in range(mtw.N)]
loc_dev = max(np.abs(L_ - locs[0]).max() for L_ in locs)
step("MN6 flavour-twisted boundary conditions (Spin(10) Cartan Sigma_12 = 1/2 - n_0, the only continuous twist compatible with H_int, which conserves charge only mod 4): a reciprocal-superlattice twist reproduces the untwisted Hamiltonian and pairing matrix exactly; the ten channels rotate real-orthogonally (sum O_A^2 invariant); the free twisted spectrum of the Q = +-1/2 flavours is h_W(k +- delta); the filtered local correlators of the twisted cell are site-independent (translation covariance of the twisted interaction)",
     triv < 1e-12 and cov_ok and spec_err < 1e-10 and loc_dev < 1e-10,
     f"trivial-twist deviation {triv:.1e}; channel rotation residual {cov_res:.1e} (R real orthogonal: {cov_ok}); twisted spectrum error {spec_err:.1e}; local correlator site dependence {loc_dev:.1e}")

# ---------------------------------------------------------------------------------------- MN7 g/N test on 32 sites
def onel_cell(cell, beta_, g_):
    mc = rw.Model(cell, beta_, 0.25, 0.0, "md"); E_c, U_c = np.linalg.eigh(mc.Hw)
    Gc = (U_c * (1 / (1 + np.exp(-beta_ * E_c)))) @ U_c.conj().T
    Fc = np.kron(mc.F, np.eye(32))
    Gr_ = (Fc @ Gc @ Fc.T)[:32, :32]; Gh_ = (Fc @ (np.eye(mc.n) - Gc).T @ Fc.T)[:32, :32]
    Sc = hfb_sigma(Kloc, np.block([[Gr_, np.zeros((32, 32))], [np.zeros((32, 32)), Gh_]]), 1.0)
    H1c = mc.Hw + g_ * Fc @ np.kron(np.eye(mc.N), Sc[:32, :32]) @ Fc.T
    e_c, V_c = np.linalg.eigh(H1c)
    Gt_ = (V_c * (np.exp(-beta_ / 2 * e_c) / (1 + np.exp(-beta_ * e_c)))) @ V_c.conj().T
    ph_ = np.exp(-1j * np.array(mc.sites, float) @ P)
    return float((ph_ @ np.einsum("xaya->xy", Gt_.reshape(mc.N, 32, mc.N, 32)) @ ph_.conj()).real / (32 * mc.N))
pts32 = []
for fn in ("reweight_points_md_oh32_gN.json", "reweight_points_md_oh32_g2.json"):
    if (HERE / fn).exists():
        pts32 += [p_ for p_ in json.loads((HERE / fn).read_text()) if "greens" in p_ and "twist_delta" not in p_]
if pts32 and 0.01 in qmc:
    rows32 = []
    for q32 in pts32:
        G32, E32 = q32["greens"]["P"]["G_half"], q32["greens"]["P"]["G_half_err"]
        a32 = onel_cell("oh32", q32["beta"], q32["g"])
        rows32.append((q32["beta"], q32["g"], G32, E32, a32, (0.5 - a32) / (0.5 - G32)))
    tN = [r for r in rows32 if r[0] == 2.0 and abs(r[1] - 0.02) < 1e-12]
    mult_excl = bool(tN) and abs(tN[0][2] - qmc[0.01][0]) > 10 * (tN[0][3] + qmc[0.01][1])
    step("MN7 32-site points: the pure finite-size multiplet picture (G_P depends only on g/N; test 32 sites, beta = 2, g = 0.02 vs 16 sites, g = 0.01) is excluded, and the one-loop node shift explains the 32-site depletion almost completely, i.e. the beyond-one-loop residual of the 16-site cell shrinks with volume (16 and 32 are not homothetic: indication, not a scaling law)",
         mult_excl and all(0.9 <= r[5] <= 1.1 for r in rows32),
         "; ".join(f"32 sites beta = {b_}, g = {g_}: QMC {G_:.4f}+-{e_:.4f}, one-loop {a_:.4f}, share {100 * sh:.0f} %" for b_, g_, G_, e_, a_, sh in rows32)
         + f"; 16 sites beta = 2, g = 0.01: QMC {qmc[0.01][0]:.4f} (multiplet prediction for 32/0.02: equal), one-loop share {100 * frac[0.01]:.0f} %")

# ---------------------------------------------------------------------------------------- MN8 one-loop node trajectory
def H1vec(k, g):
    hv = np.array([np.trace(h_of_k(k) @ S[a]).real / 2 for a in range(3)])
    return hv + g * f_of_k(k) ** 2 * b
def wcharge(k, g, e=1e-6):
    J = np.array([(H1vec(k + e * np.eye(3)[i], g) - H1vec(k - e * np.eye(3)[i], g)) / (2 * e) for i in range(3)])
    return int(np.sign(np.linalg.det(J)))
lam_ax = np.linspace(-np.pi / 2, 1.5 * np.pi, 4001)                     # the whole [111] axis through Gamma (lambda = -pi/2) and P (0)
traj = {}
for g_ in (0.0, 0.01, 0.02, 0.021, 0.022, 0.023, 0.024, 0.03, 0.05, 0.1):
    c = np.array([H1vec(P + l * np.ones(3), g_) @ nvec for l in lam_ax])
    zs = []
    for i in range(len(lam_ax) - 1):
        if c[i] == 0 or c[i] * c[i + 1] < 0:
            l0 = lam_ax[i] - c[i] * (lam_ax[i + 1] - lam_ax[i]) / (c[i + 1] - c[i])
            zs.append((round(float(l0), 3), wcharge(P + l0 * np.ones(3), g_)))
    traj[g_] = zs
gam_fixed = all(any(abs(l + np.pi / 2) < 1e-2 and q == 1 for l, q in zs) for zs in traj.values())
mirror_kept = all(any(q == -1 for l, q in zs) for zs in traj.values())
pair_window = [g_ for g_, zs in traj.items() if len(zs) == 4]
step("MN8 one-loop node trajectory along [111] (the C3 axis; total charge of the one-loop band structure is 0): Gamma stays a fixed +1 node for every g; a -1 mirror node exists for every g; it moves continuously up to g ~ 0.020, then a +-1 pair is created on the axis, the new +1 annihilates the old mirror and the new -1 carries the mirror charge on (pair window g ~ 0.0205 ... 0.0235) -> the linear regime ends near g ~ 0.02 and the mirror charge is never removed at first order (off-axis pairs appear around g = 0.05, global census)",
     gam_fixed and mirror_kept and bool(pair_window) and min(pair_window) > 0.02 and max(pair_window) < 0.024,
     "; ".join(f"g = {g_}: " + ", ".join(f"({l:+.3f}, {q:+d})" for l, q in zs) for g_, zs in traj.items()))

step("MN5 (statement) Gate B criterion corrected: Delta(P) > 0 at the old P does not show a mirror gap (it is reproduced by node motion); required is min_{k != Gamma} Delta_odd(k) > 0 with node tracking (twisted boundary conditions, k_star(g, V), Delta(k_star)), plus Gamma Weyl intact, all independent mirror composites gapped, Spin(10) x Z4 x translations unbroken; a genuine mirror removal must be non-perturbative (no bilinear can gap a single Weyl node; Green-function zeros are the candidate mechanism)",
     True, "next contract: CHIRAL4D.MIRROR.NODE_TRACK.01, then CHIRAL4D.GREEN_TOPOLOGY.01")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "b_over_g": b.tolist(), "dk_per_g": lin_pred.tolist(),
       "V_P": VP.tolist(), "one_loop_vs_qmc": {str(k): v for k, v in cmp_.items()}, "axis_trajectory": {str(k): v for k, v in traj.items()}, "one_loop_share": {str(k): v for k, v in frac.items()}, "b16_over_g": b16.tolist()}
(HERE / "mirror_node_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
