#!/usr/bin/env python3
"""Part PF: QMC.PFAFFIAN.GENERAL.01 -- sign-exact HS weights of Theta_1 without any square root, and the HMC force.

W[phi] = Tr prod_l X_l, X_l = e^{1/2 Psi^dag K~_l Psi} the slice operator whose Nambu image is S_l = e^{K_w/2} e^{K_phi,l}
e^{K_w/2} (reweight.py).  Every X_l is a positive Hermitian Gaussian; its Grassmann representation is
omega(X_l) = c_l exp(1/2 theta^T T_l theta) with T_l = Om^dag (S_l - 1)(S_l + 1)^-1 Om (Majorana basis, Psi = Om gamma / sqrt2)
and c_l = 2^-n sqrt det(1 + S_l) > 0 (no branch: positive operator).  Products compose exactly by a Grassmann integral:
c_XY = c_X c_Y Pf([[T_X, 1], [-1, T_Y]]) / Pf([[0, 1], [-1, 0]]),  T_XY = A^-1_11 + A^-1_22 - A^-1_12 - A^-1_21,
and W = 2^n c_{0..L-1}.  The Pfaffian is a polynomial in the entries: sign and phase are exact, no square root and no
continuity tracking.  log W is accumulated as log|.| + i arg.  The phase-quenched HMC force is
d log W / d phi = 1/2 Tr[(1 + B)^-1 dB] (single valued although W is not), via the exact Frechet derivative of e^{K_phi}.
Unit tests (the seven required ones):
  PF1 brute force (Fock space, n = 3): Pfaffian routine, single-Gaussian representation (c, T), composition and traces of
      products of six random positive BdG Gaussians with complex pairing -- including phase; the principal square root of
      det(1 + B) gets the sign wrong on a sizeable fraction of these products.
  PF2 4-site comparison: on the debug4 cell the frozen h_W commutes with n.sigma, n = [111]/sqrt3; a particle-hole
      transformation of the down_n spins gives the number-conserving special formula
      W = 2^32 e^{-beta 64 sqrt3 t} det(1 + S_active), exact with sign/phase; the Pfaffian weight must reproduce it on
      thousands of random HS histories (several field strengths), and continuous tracking of det(1 + B) must agree.
  PF3 limits: g = 0 (W = det(1 + e^{-beta H_w}) on debug4, cell8, cubic16) and t = 0 (spin SU(2) exact -> the
      number-conserving formula holds on every cell; checked on debug4, cell8 and cubic16).
  PF4 sign crossing: on a real-weight family (n = 3, real BdG) W(s) changes sign; the Pfaffian follows the brute-force
      sign through the crossing, sqrt(det) with continuity tracking cannot see it (det = W^2 touches zero).
  PF5 HMC force: analytic d log W / d phi (real and imaginary part) against central finite differences of the Pfaffian
      log W on debug4.
  PF6 leapfrog: exact reversibility (forward, momentum flip, backward) and Delta H = O(eps^2).
  PF7 stability: debug4, beta up to 16 (64 slices): Pfaffian log|W| = 1/2 log det_UDV(1 + B) and 2 arg W = arg det mod
      2 pi, while the unstabilised product fails.
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/qmc_pfaffian.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("reweight", HERE / "reweight.py")
rw = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(rw)
LOG: list[dict] = []
wrap = lambda a: (a + np.pi) % (2 * np.pi) - np.pi


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    sys.stdout.flush()


# ---------------------------------------------------------------------------------------- Pfaffian Gaussian algebra
def pfaffian(A):
    """Parlett-Reid with pivoting (Wimmer's LTL form); returns (log|Pf|, arg Pf) with arg accumulated, not wrapped."""
    A = A.astype(complex).copy(); n = A.shape[0]; la, ph = 0.0, 0.0
    for k in range(0, n - 1, 2):
        kp = k + 1 + np.abs(A[k + 1:, k]).argmax()
        if kp != k + 1:
            A[[k + 1, kp], k:] = A[[kp, k + 1], k:]; A[k:, [k + 1, kp]] = A[k:, [kp, k + 1]]; ph += np.pi
        p = A[k, k + 1]
        if p == 0:
            return -np.inf, 0.0
        la += np.log(abs(p)); ph += np.angle(p)
        if k + 2 < n:
            tau = A[k, k + 2:] / p; col = A[k + 2:, k + 1].copy()
            A[k + 2:, k + 2:] += np.outer(tau, col) - np.outer(col, tau)
    return la, ph


_PF0: dict[int, complex] = {}
def omega_mat(n):
    I = np.eye(n)
    return np.block([[I, 1j * I], [I, -1j * I]]) / np.sqrt(2)


def gauss_from_S(S, Om):
    """positive Hermitian Nambu image S -> (log c, T)."""
    m = S.shape[0]; n = m // 2
    T = Om.conj().T @ (S - np.eye(m)) @ np.linalg.inv(S + np.eye(m)) @ Om
    return 0.5 * np.linalg.slogdet(np.eye(m) + S)[1] - n * np.log(2), (T - T.T) / 2


def compose(a, b, want_T=True):
    lca, Ta = a; lcb, Tb = b
    m = Ta.shape[0]; I = np.eye(m)
    A = np.block([[Ta, I], [-I, Tb]])
    if m not in _PF0:
        la0, ph0 = pfaffian(np.block([[0 * I, I], [-I, 0 * I]])); _PF0[m] = la0 + 1j * ph0
    la, ph = pfaffian(A)
    lc = lca + lcb + (la + 1j * ph) - _PF0[m]
    if not want_T:
        return lc, None
    Ai = np.linalg.inv(A)
    S = Ai[:m, :m] + Ai[m:, m:] - Ai[:m, m:] - Ai[m:, :m]
    return lc, (S - S.T) / 2


def logW_pf(Ss):
    """log W = log Tr prod_l X_l (complex: log|W| + i arg W, arg accumulated) from the Nambu slice images S_0 ... S_{L-1}."""
    n = Ss[0].shape[0] // 2
    Om = omega_mat(n)
    acc = gauss_from_S(Ss[0], Om)
    for l in range(1, len(Ss)):
        acc = compose(acc, gauss_from_S(Ss[l], Om), want_T=(l < len(Ss) - 1))
    return acc[0] + n * np.log(2)


# ---------------------------------------------------------------------------------------- PF1 brute force
rng = np.random.default_rng(11)
X6 = rng.normal(size=(6, 6)) + 1j * rng.normal(size=(6, 6)); X6 = X6 - X6.T
def pf_rec(A):
    if A.shape[0] == 0:
        return 1.0
    s = 0
    for j in range(1, A.shape[0]):
        idx = [k for k in range(A.shape[0]) if k not in (0, j)]
        s += (-1) ** (j + 1) * A[0, j] * pf_rec(A[np.ix_(idx, idx)])
    return s
la, ph = pfaffian(X6)
X40 = rng.normal(size=(40, 40)) + 1j * rng.normal(size=(40, 40)); X40 = X40 - X40.T
la40, ph40 = pfaffian(X40); sd40, ld40 = np.linalg.slogdet(X40)
pf_ok = abs(np.exp(la + 1j * ph) - pf_rec(X6)) < 1e-10 and abs(2 * la40 - ld40) < 1e-9 and abs(wrap(2 * ph40 - np.angle(sd40))) < 1e-9

nb = 3
Zp = np.diag([1.0, -1.0]); low = np.array([[0, 1], [0, 0]], complex)
cf = []
for i in range(nb):
    mm = np.eye(1, dtype=complex)
    for j in range(nb):
        mm = np.kron(mm, Zp if j < i else (low if j == i else np.eye(2)))
    cf.append(mm)
Psi = cf + [x.conj().T for x in cf]
gam = [x + x.conj().T for x in cf] + [-1j * (x - x.conj().T) for x in cf]
Om3 = omega_mat(nb)
def randK(s=1.0, real=False):
    if real:
        h = rng.normal(size=(nb, nb)); h = (h + h.T) / 2
        D = rng.normal(size=(nb, nb)); D = D - D.T
    else:
        h = rng.normal(size=(nb, nb)) + 1j * rng.normal(size=(nb, nb)); h = (h + h.conj().T) / 2
        D = rng.normal(size=(nb, nb)) + 1j * rng.normal(size=(nb, nb)); D = D - D.T
    return s * np.block([[h, D.conj().T], [D, -h.T]])
def op(K):
    return expm(0.5 * sum(K[i, j] * Psi[i].conj().T @ Psi[j] for i in range(2 * nb) for j in range(2 * nb)))
def omega_bf(X):
    c0 = np.trace(X) / 2 ** nb; T = np.zeros((2 * nb, 2 * nb), complex)
    for i in range(2 * nb):
        for j in range(i + 1, 2 * nb):
            T[i, j] = np.trace((gam[i] @ gam[j]).conj().T @ X) / 2 ** nb / c0; T[j, i] = -T[i, j]
    return c0, T
K0 = randK(); c0, Tb = omega_bf(op(K0)); lc0, T0 = gauss_from_S(expm(K0), Om3)
single_ok = abs(np.exp(lc0) - c0) < 1e-10 * abs(c0) and np.abs(T0 - Tb).max() < 1e-10
K1, K2 = randK(), randK(); cb, Tb12 = omega_bf(op(K1) @ op(K2))
lc12, T12 = compose(gauss_from_S(expm(K1), Om3), gauss_from_S(expm(K2), Om3))
comp_ok = abs(np.exp(lc12) - cb) < 1e-10 * abs(cb) and np.abs(T12 - Tb12).max() < 1e-9
errs, wrong_sqrt, neg = [], 0, 0
for trial in range(200):
    Ks = [randK(1.2) for _ in range(6)]
    Wbf = np.trace(np.linalg.multi_dot([op(K) for K in Ks]))
    lw = logW_pf([expm(K) for K in Ks])
    errs.append(abs(np.exp(lw) - Wbf) / abs(Wbf))
    dB = np.linalg.det(np.eye(2 * nb) + np.linalg.multi_dot([expm(K) for K in Ks]))
    wrong_sqrt += abs(np.sqrt(dB) - Wbf) > 1e-6 * abs(Wbf)
    neg += Wbf.real < 0
step("PF1 brute force (Fock space, 3 modes): Pfaffian routine exact (recursive expansion, Pf^2 = det); single positive Gaussian (c, T = Om^dag tanh(K/2) Om) and the composition rule reproduce the Fock-space Grassmann coefficients; 200 traces of six random BdG Gaussians with complex pairing agree including phase, although the principal sqrt det(1 + B) has the wrong sign on many of them",
     pf_ok and single_ok and comp_ok and max(errs) < 1e-9 and wrong_sqrt > 20,
     f"max rel. error {max(errs):.1e}; principal sqrt wrong in {wrong_sqrt}/200; Re W < 0 in {neg}/200")


# ---------------------------------------------------------------------------------------- Theta_1 helpers
def build_M(m, phi_l):
    M = np.zeros((m.n, m.n), complex)
    for i in range(m.N):
        M += np.kron(m.FF[i], np.tensordot(phi_l[i], rw.MA, axes=(0, 0)))
    return M


def slices(m, phi):
    return [m.slice_op(phi[l]) for l in range(m.Lt)]


def nc_projectors(m, nvec):
    """isometries onto up_n / down_n (site x spin x flavour ordering 32 i + 16 s + f)."""
    ns = sum(nvec[a] * [rw.sx, rw.sy, rw.sz][a] for a in range(3))
    ev, U = np.linalg.eigh(ns)
    up, dn = U[:, 1], U[:, 0]
    Cu = np.kron(np.eye(m.N), np.kron(up[:, None], np.eye(16)))
    Cd = np.kron(np.eye(m.N), np.kron(dn[:, None], np.eye(16)))
    Tchi = np.block([[Cu.conj().T, np.zeros((Cu.shape[1], m.n))], [np.zeros((Cd.shape[1], m.n)), Cd.T]])
    Tbar = np.block([[np.zeros((Cu.shape[1], m.n)), Cu.T], [Cd.conj().T, np.zeros((Cd.shape[1], m.n))]])
    return Tchi, Tbar, Cu, Cd


def logW_nc(m, Ss, Tchi, Tbar, Cu):
    """number-conserving special formula after the particle-hole transformation of down_n:
    W = det(1 + S_chi) / sqrt(det S_chi) with the closed-form prefactor sqrt(det S_chi) = prod_l det(T W_h T^dag)
    = exp(-beta Tr H_up) (the pairing slices are traceless in chi: det = 1; Tr H = 0) -- so W = e^{beta Tr H_up} det(1 + S_chi),
    no numerical determinant of an ill-conditioned slice.  Returns log W, the relative chi/chibar leak, the chi slices."""
    chi = [Tchi @ S @ Tchi.conj().T for S in Ss]
    leak = max(np.abs(Tbar @ S @ Tchi.conj().T).max() / np.abs(S).max() for S in Ss)
    la, ph, _ = m.logdet(chi)
    TrHu = np.trace(Cu.conj().T @ m.Hw @ Cu).real
    return la + 1j * ph + m.beta * TrHu, leak, chi


def logW_track(m, phi, nsub=1):
    """continuous tracking of arg det(1 + B) along t phi, t: 0 -> 1 (reweight.py procedure); returns 1/2 log det."""
    def ld(t):
        a, p, _ = m.logdet([m.slice_op(t * phi[l]) for l in range(m.Lt)])
        return a, p
    th, t_, (_, p_prev), dt = 0.0, 0.0, ld(0.0), 0.125
    while t_ < 1 - 1e-12:
        dt = min(dt, 1 - t_)
        a, p = ld(t_ + dt)
        inc = wrap(p - p_prev)
        if abs(inc) > np.pi / 4 and dt > 1e-4:
            dt /= 2; continue
        th += inc; p_prev = p; t_ += dt; dt *= 1.5
    return 0.5 * (a + 1j * th)


# ---------------------------------------------------------------------------------------- PF2 four-site special formula
t0 = time.time()
nvec = np.ones(3) / np.sqrt(3)
md4 = rw.Model("debug4", 1.0, 0.25, 0.1, "md")
NS = np.kron(np.eye(md4.N), np.kron(sum(nvec[a] * [rw.sx, rw.sy, rw.sz][a] for a in range(3)), np.eye(16)))
comm_n = np.abs(md4.Hw @ NS - NS @ md4.Hw).max()
Tchi, Tbar, Cu, Cd = nc_projectors(md4, nvec)
TrHu = np.trace(Cu.conj().T @ md4.Hw @ Cu).real
# Gamma subspace (f = 0, h = 0): projector on the uniform mode of the cell
Pg = np.kron(np.ones((md4.N, 1)) / np.sqrt(md4.N), np.eye(16))
Gam_chi = np.block([[np.kron(np.ones((md4.N, 1)) / np.sqrt(md4.N), np.eye(16)), np.zeros((16 * md4.N, 16))],
                    [np.zeros((16 * md4.N, 16)), np.kron(np.ones((md4.N, 1)) / np.sqrt(md4.N), np.eye(16))]])
res2 = []
gam_fix = 0.0
hist = 0
for g_, beta_, nh in ((0.1, 1.0, 600), (0.5, 1.0, 600), (2.0, 1.0, 600), (0.5, 2.0, 120)):
    mloc = rw.Model("debug4", beta_, 0.25, g_, "md")
    sig = np.sqrt(2 * g_ * 0.25)
    for h_ in range(nh):
        phi = sig * rng.normal(size=(mloc.Lt, mloc.N, 10))
        Ss = slices(mloc, phi)
        lp = logW_pf(Ss)
        ln, leak, chi = logW_nc(mloc, Ss, Tchi, Tbar, Cu)
        for c_ in chi:                                                # Gamma block of every chi slice is the identity and decoupled
            gam_fix = max(gam_fix, np.abs(Gam_chi.conj().T @ c_ @ Gam_chi - np.eye(32)).max() / np.abs(c_).max(),
                          np.abs(c_ @ Gam_chi - Gam_chi).max() / np.abs(c_).max())
        res2.append((g_, beta_, abs(lp.real - ln.real), abs(wrap(lp.imag - ln.imag)), wrap(ln.imag), leak, abs(ln.real)))
        hist += 1
res2 = np.array(res2)
# closed form of the prefactor: exp(-1/2 sum log det S_l) = exp(beta Tr H_u); |Tr H_u| = 64 sqrt3 t
ntrack, track_err = 0, 0.0
mtr = rw.Model("debug4", 1.0, 0.25, 2.0, "md")
for h_ in range(60):
    phi = np.sqrt(2 * 2.0 * 0.25) * rng.normal(size=(mtr.Lt, mtr.N, 10))
    Ss = slices(mtr, phi)
    lt = logW_track(mtr, phi); ln, _, _ = logW_nc(mtr, Ss, Tchi, Tbar, Cu)
    track_err = max(track_err, abs(lt.real - ln.real) + abs(wrap(lt.imag - ln.imag))); ntrack += 1
big_phase = float(np.mean(np.abs(res2[:, 4]) > np.pi / 2))
# outside the tested range (not gating): prior fields at g = 5 (sigma_phi = 1.6) reach the double-precision conditioning
# limit of both evaluations; the production chains stay below the g = 2 prior range (phi_max ~ 1.7 vs ~ 3.5 here)
mx = rw.Model("debug4", 1.0, 0.25, 5.0, "md"); ext = 0.0
for h_ in range(20):
    phi = np.sqrt(2 * 5.0 * 0.25) * rng.normal(size=(mx.Lt, mx.N, 10))
    Ss = slices(mx, phi); lp = logW_pf(Ss); ln, _, _ = logW_nc(mx, Ss, Tchi, Tbar, Cu)
    ext = max(ext, abs(lp.real - ln.real) + abs(wrap(lp.imag - ln.imag)))
per_g = {}
for g_ in sorted(set(res2[:, 0])):
    r_ = res2[res2[:, 0] == g_]
    per_g[float(g_)] = (float(r_[:, 2].max()), float(r_[:, 3].max()), float((r_[:, 2] / r_[:, 6]).max()))
crit = all((e1 < 1e-8 and e2 < 1e-8) if g_ <= 0.5 else (e1 < 1e-4 and e2 < 1e-4 and rel < 1e-7) for g_, (e1, e2, rel) in per_g.items())
step("PF2 4-site comparison (debug4, frozen h_W): h_W commutes with n.sigma, n = [111]/sqrt3; after the particle-hole transformation of down_n, W = 2^32 e^{-beta 64 sqrt3 t} det(1 + S_active) exactly (Gamma block of every slice = 1, Tr H_up = -64 sqrt3); the Pfaffian weight reproduces it on 1920 random HS histories (g = 0.1, 0.5: |d log W| < 1e-8; g = 2, log|W| ~ 2e3: relative < 1e-7, conditioning-limited), continuous tracking of det(1 + B) agrees; on this cell all weights are real positive (the non-trivial sign tests are PF1 and PF4)",
     comm_n < 1e-12 and res2[:, 5].max() < 1e-12 and abs(abs(TrHu) - 64 * np.sqrt(3)) < 1e-9 and gam_fix < 1e-12 and crit and track_err < 1e-4,
     f"{hist} histories; per g (max |dlog|W||, max |d arg|, max relative): { {k: tuple(f'{x:.1e}' for x in v) for k, v in per_g.items()} }; |arg W| > pi/2 in {100 * big_phase:.1f}%, max |arg W| {np.abs(res2[:, 4]).max():.1e}; "
     f"Tr H_up = {TrHu:.4f} (64 sqrt3 = {64 * np.sqrt(3):.4f}); relative Gamma-block deviation {gam_fix:.1e}; "
     f"tracking vs exact on {ntrack} histories (g = 2) {track_err:.1e}; outside the range, g = 5 prior (conditioning limit, not gating): {ext:.1e}; {time.time() - t0:.0f} s")


# ---------------------------------------------------------------------------------------- PF3 limits
t0 = time.time()
free = {}
for cell, beta_ in (("debug4", 2.0), ("cell8", 1.0), ("cubic16", 1.0)):
    mloc = rw.Model(cell, beta_, 0.25, 0.0, "md")
    Ss = [mloc.Wh @ mloc.Wh for _ in range(mloc.Lt)]
    lp = logW_pf(Ss)
    E = np.linalg.eigvalsh(mloc.Hw)
    exact = np.sum(np.log1p(np.exp(-beta_ * E)))
    free[cell] = (abs(lp.real - exact), abs(wrap(lp.imag)))
t0lim = {}
for cell, beta_, nh in (("debug4", 1.0, 100), ("cell8", 1.0, 6), ("cubic16", 0.5, 1)):
    mloc = rw.Model(cell, beta_, 0.25, 0.5, "md")
    mloc.Wh = np.eye(2 * mloc.n); mloc.Hw = np.zeros_like(mloc.Hw)                  # t = 0: no kinetic factor, no prefactor
    Tc, Tb_, Cu_, _ = nc_projectors(mloc, np.array([0.0, 0.0, 1.0]))
    e_ = []
    for h_ in range(nh):
        phi = np.sqrt(2 * 0.5 * 0.25) * rng.normal(size=(mloc.Lt, mloc.N, 10))
        Ss = slices(mloc, phi)
        lp = logW_pf(Ss); ln, leak, _ = logW_nc(mloc, Ss, Tc, Tb_, Cu_)
        e_.append(abs(lp.real - ln.real) + abs(wrap(lp.imag - ln.imag)) + leak)
    t0lim[cell] = max(e_)
step("PF3 limits: free g = 0 (W = det(1 + e^{-beta H_w}) on debug4, cell8, cubic16) and interaction-only t = 0 (spin SU(2) exact, number-conserving formula on every cell: debug4, cell8, cubic16) reproduced by the Pfaffian weight",
     all(a < 1e-8 and b < 1e-8 for a, b in free.values()) and all(v < 1e-8 for v in t0lim.values()),
     f"free: { {k: f'{a:.1e}' for k, (a, b) in free.items()} }; t = 0: { {k: f'{v:.1e}' for k, v in t0lim.items()} }; {time.time() - t0:.0f} s")


# ---------------------------------------------------------------------------------------- PF4 sign crossing
pos = neg_ = None
for trial in range(400):
    Ks = [randK(1.5, real=True) for _ in range(4)]
    W = np.trace(np.linalg.multi_dot([op(K) for K in Ks])).real
    if W > 0.5 and pos is None:
        pos = Ks
    if W < -0.5 and neg_ is None:
        neg_ = Ks
    if pos is not None and neg_ is not None:
        break
ss = np.linspace(0, 1, 81)
Wbf, Wpf, Wtr = [], [], []
th, p_prev = 0.0, None
for s in ss:
    Ks = [(1 - s) * a + s * b for a, b in zip(pos, neg_)]
    Wbf.append(np.trace(np.linalg.multi_dot([op(K) for K in Ks])).real)
    Wpf.append(np.exp(logW_pf([expm(K) for K in Ks])))
    d = np.linalg.det(np.eye(2 * nb) + np.linalg.multi_dot([expm(K) for K in Ks]))
    if p_prev is not None:
        th += wrap(np.angle(d) - p_prev)
    p_prev = np.angle(d)
    Wtr.append(np.sqrt(abs(d)) * np.exp(0.5j * th))
Wbf, Wpf, Wtr = np.array(Wbf), np.array(Wpf), np.array(Wtr)
cross = int(np.sum(np.diff(np.sign(Wbf)) != 0))
pf_follow = np.max(np.abs(Wpf - Wbf) / np.maximum(np.abs(Wbf), 1e-3)) < 1e-8
tr_fail = np.max(np.abs(Wtr - Wbf)) > 0.1
step("PF4 sign crossing: on a real-weight family (3 modes, real BdG, 4 slices) W(s) changes sign; the Pfaffian weight follows the brute-force sign through the crossing, sqrt det(1 + B) with continuity tracking stays on one sheet (det = W^2 only touches zero) and is wrong after the crossing",
     cross >= 1 and pf_follow and tr_fail, f"sign changes {cross}; Pfaffian max rel. error {np.max(np.abs(Wpf - Wbf) / np.maximum(np.abs(Wbf), 1e-3)):.1e}; tracked-sqrt max error {np.max(np.abs(Wtr - Wbf)):.2f}")


# ---------------------------------------------------------------------------------------- PF5 HMC force
MAd = np.array([x.conj().T for x in rw.MA])
def dlogW(m, phi):
    """d log W / d phi_(l,i,A) = 1/2 Tr[(1 + B)^-1 dB], evaluated in the symmetrised, stabilised form
    1/2 Tr[(1 - G~_l) E_l^-1/2 dE_l E_l^-1/2],  G~_l = (1 + E_l^1/2 W_h S_{l+1}..S_{L-1} S_0..S_{l-1} W_h E_l^1/2)^-1
    (cyclic UDV chain, pivoted QR), E^-1/2 dE E^-1/2 = U [dK^ o sinhc((lam_i - lam_j)/2)] U^dag: the amplification is
    sinh(Delta/2)/(Delta/2) of one slice instead of the condition number of B."""
    Ss = slices(m, phi); L = len(Ss); nn = 2 * m.n; I = np.eye(nn, dtype=complex)
    out = np.zeros((L, m.N, 10), complex)
    for l in range(L):
        M = build_M(m, phi[l]); K = np.block([[np.zeros_like(M), 2 * M.conj().T], [2 * M, np.zeros_like(M)]])
        lam, U = np.linalg.eigh(K)
        Eh = (U * np.exp(lam / 2)) @ U.conj().T
        udv = (I, np.ones(nn), I)
        for Mx in [Eh, m.Wh] + [Ss[k] for k in range(l - 1, -1, -1)] + [Ss[k] for k in range(L - 1, l, -1)] + [m.Wh, Eh]:
            udv = m._udv_left(Mx, udv)
        Uq, Dq, Vq = udv
        Db, Ds = np.maximum(Dq, 1.0), np.minimum(Dq, 1.0)
        Vi = np.linalg.inv(Vq)
        Gt = Vi @ np.linalg.solve((Uq.conj().T @ Vi) / Db[:, None] + np.diag(Ds), Uq.conj().T / Db[:, None])
        dl = (lam[:, None] - lam[None, :]) / 2
        with np.errstate(divide="ignore", invalid="ignore"):
            Sig = np.where(np.abs(dl) > 1e-10, np.sinh(dl) / np.where(np.abs(dl) > 1e-10, dl, 1), 1.0)
        Z = U @ ((U.conj().T @ (I - Gt) @ U) * Sig) @ U.conj().T
        n = m.n
        Z12 = Z[:n, n:].reshape(m.N, 32, m.N, 32); Z21 = Z[n:, :n].reshape(m.N, 32, m.N, 32)
        # Tr(Z dK) = Tr(Z_12 dK_21) + Tr(Z_21 dK_12) = 2 Tr(Z_12 dM) + 2 Tr(Z_21 dM^dag), dM = FF_i (x) M_A
        t12 = np.einsum("iy,ix,xyA->iA", m.F, m.F, np.einsum("xayb,Aba->xyA", Z12, rw.MA))
        t21 = np.einsum("iy,ix,xyA->iA", m.F, m.F, np.einsum("xayb,Aba->xyA", Z21, MAd))
        out[l] = t12 + t21
    return out


mf = rw.Model("debug4", 1.0, 0.25, 0.5, "md")
phi5 = np.sqrt(2 * 0.5 * 0.25) * rng.normal(size=(mf.Lt, mf.N, 10))
an = dlogW(mf, phi5)
fd_err = []
for _ in range(16):
    l, i, A = rng.integers(mf.Lt), rng.integers(mf.N), rng.integers(10)
    h = 1e-5
    pp = phi5.copy(); pp[l, i, A] += h; pm = phi5.copy(); pm[l, i, A] -= h
    lp, lm = logW_pf(slices(mf, pp)), logW_pf(slices(mf, pm))
    fd = ((lp.real - lm.real) + 1j * wrap(lp.imag - lm.imag)) / (2 * h)
    fd_err.append(abs(fd - an[l, i, A]) / max(abs(an[l, i, A]), 1e-3))
step("PF5 HMC force: analytic d log W / d phi = 1/2 Tr[(1 + B)^-1 dB] in the symmetrised stabilised form (cyclic UDV Green function, exact Frechet derivative of e^{K_phi}; real part = phase-quenched force, imaginary part = d theta) agrees with central finite differences of the Pfaffian log W on debug4 at g = 0.5 (cond B ~ 1e18; 16 random components)",
     max(fd_err) < 1e-6, f"max rel. deviation {max(fd_err):.1e}; |force| range {np.abs(an).min():.2e} .. {np.abs(an).max():.2e}")


# ---------------------------------------------------------------------------------------- PF6 leapfrog
def H_and_F(m, phi):
    la_, _, _ = m.logdet(slices(m, phi))
    S_ = np.sum(phi ** 2) / (4 * m.g * m.dtau) - 0.5 * la_
    Fo = -phi / (2 * m.g * m.dtau) + dlogW(m, phi).real
    return S_, Fo
def leapfrog(m, phi, pi, eps, nsteps):
    phi, pi = phi.copy(), pi.copy()
    S_, Fo = H_and_F(m, phi)
    H0 = 0.5 * np.sum(pi ** 2) + S_
    for _ in range(nsteps):
        pi += 0.5 * eps * Fo
        phi += eps * pi
        S_, Fo = H_and_F(m, phi)
        pi += 0.5 * eps * Fo
    return phi, pi, 0.5 * np.sum(pi ** 2) + S_ - H0
mh = rw.Model("debug4", 0.5, 0.25, 0.5, "md")
phi6 = np.sqrt(2 * 0.5 * 0.25) * rng.normal(size=(mh.Lt, mh.N, 10)); pi6 = rng.normal(size=phi6.shape)
tau_md = 0.2
p1, q1, dH1 = leapfrog(mh, phi6, pi6, tau_md / 8, 8)
p2, q2, _ = leapfrog(mh, p1, -q1, tau_md / 8, 8)
rev = max(np.abs(p2 - phi6).max(), np.abs(-q2 - pi6).max())
dHs = [abs(leapfrog(mh, phi6, pi6, tau_md / ns, ns)[2]) for ns in (4, 8, 16)]
ratios6 = [dHs[0] / dHs[1], dHs[1] / dHs[2]]
step("PF6 leapfrog: forward trajectory, momentum flip, backward trajectory returns to the start (reversibility), and |Delta H| scales like eps^2 (ratio 4 per halving)",
     rev < 1e-9 and all(2.8 < r < 5.5 for r in ratios6), f"reversibility error {rev:.1e}; |Delta H| {[f'{x:.2e}' for x in dHs]}; ratios {np.round(ratios6, 2).tolist()}")


# ---------------------------------------------------------------------------------------- PF7 stability
stab = []
for beta_ in (1.0, 4.0, 8.0, 16.0):
    ms = rw.Model("debug4", beta_, 0.25, 0.1, "md")
    phi = np.sqrt(2 * 0.1 * 0.25) * rng.normal(size=(ms.Lt, ms.N, 10))
    Ss = slices(ms, phi)
    lp = logW_pf(Ss)
    la_, ph_, _ = ms.logdet(Ss)
    sn, lnaive = np.linalg.slogdet(np.eye(2 * ms.n) + np.linalg.multi_dot(Ss) if len(Ss) > 1 else np.eye(2 * ms.n) + Ss[0])
    stab.append((beta_, abs(2 * lp.real - la_) / abs(la_), abs(wrap(2 * lp.imag - ph_)), abs(lnaive - la_) / abs(la_), lp.real))
stab_ok = all(a < 1e-9 and b < 1e-8 for _, a, b, _, _ in stab)
step("PF7 stability (debug4, g = 0.1, beta = 1 ... 16, up to 64 slices): Pfaffian log|W| equals 1/2 log det(1 + B) of the stabilised UDV chain and 2 arg W = arg det mod 2 pi; the unstabilised product drifts with beta",
     stab_ok, "; ".join(f"beta={b:g}: rel |2 log|W| - UDV| {a:.1e}, phase {p:.1e}, naive rel. error {nv:.1e}, log|W| {lw:.1f}" for b, a, p, nv, lw in stab))

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG),
       "pf2_histories": int(hist), "pf2_max_dlog": float(res2[:, 2].max()), "pf2_max_darg": float(res2[:, 3].max()), "TrH_up": TrHu,
       "stability": stab}
(HERE / "qmc_pfaffian_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
