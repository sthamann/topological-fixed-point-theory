#!/usr/bin/env python3
"""Part SPP: SECTOR.PROJECTED.GREEN.PRODUCTION.01 -- the start-sector projected Green function in the Theta_1 QMC.

    C_r(tau) = Tr[P_r e^{-(beta - tau) H} O e^{-tau H} O^dag] / Tr[P_r e^{-beta H}],   P_r = (1/64) sum_U chi_r(U)^* U,
U = Gamma(i^m u_a): the 16 (twisted) cell translations u_a times the four Z4 elements e^{i pi m N / 2}.  Only the start
projector is inserted: the HS-averaged integrator commutes with every U (ST4; twisted translations SPP3), and O = c_k (one
momentum, Z4 charge 1) maps sectors to sectors, so P_s is redundant in the average (it is not redundant per configuration).
Per HS configuration X = S_0 .. S_{L-1}:
    f_U = Tr[U X] / Tr[X]  (sign-exact): f_U^2 = det(u) det(1 + S_U S_X) / det(1 + S_X); the root is continued along
          u(t) = F e^{i t phi} F^-1 from t = 0 (f = 1) in the common eigenbasis F of all u (plane waves), so every step is one
          stabilised LU of the 2n x 2n matrix conj(U_X) / D_b + S_u(t) V_X^T D_s;
    g_U(tau) = Tr[U X_A O X_C O^dag] / Tr[U X] = tr Pi^dag S_C (1 + S_U S_A S_C)^-1 Pi / 16  (Wick for the Gaussian U X):
          the UDV chains of greens_k with V_A -> V_A S_U^T (a unitary on the stable side; S_U monomial).
Reweighted: D_r = <e^{i theta} sum_U chi_r^* f_U / 64> / <e^{i theta}>  (sector fraction Z_r / Z; reference sector = the
largest D_r at the coldest beta, from the data, not assumed), N_r(tau) likewise with f_U g_U(tau), C_r = N_r / D_r with its
ABSOLUTE normalisation (C_r(0) = <P_r O O^dag> / <P_r>; no per-curve renormalisation).  Fixed physical times TAUS are compared
across beta (first beta -> infinity at fixed tau, then the tau dependence).
  SPP1 brute force (3 modes, random BdG slices with pairing, random unitaries incl. an eigenvalue -1): f_U and f_U g_U from
       the production functions equal the Fock-space traces.
  SPP2 free limit of the production operator (twisted cubic16, g = 0): f_U and g_U(tau) for all 64 elements equal the closed
       forms det(1 + u e^{-beta h}) / det(1 + e^{-beta h}) and tr Pi^dag e^{-tau h} (1 + u e^{-beta h})^-1 Pi / 16.
  SPP3 symmetry of the twisted production operator: [u_a, H_w] = 0, u_a are a representation diagonal in F, and the HS
       composites map covariantly (u^T B_{iA} u in span{B_{pi(i) B}} with an orthogonal 10 x 10 matrix) -> the HS-averaged
       slice commutes with every U.
  SPP4 interacting configuration (g = 0.018): the U = 1 term reproduces greens_k, and f_U^2 from the path equals the
       independent determinant ratio det(u) det(mid_U) / det(mid_1) at several cuts l.
  SPP5 (statement) production: sector fractions, reference sector and C_r(tau) at fixed tau for beta = 2, 4, 6, 8, with the
       free closed form and the first-order HF quasi-particle model as references.
  SPP7 SPG budget test (SPG.BUDGET.01, criteria fixed 2026-10-03 before applying them; B2 sharpened after the first
       application, see code): B1 R_phase > 0.5; B2 projection active (D_ref >= 2/64 and D_ref - 1/64 > 3 sigma, else
       SPG_INACTIVE: C_r = G carries no extra information; sigma >= 1/64 counts as noise); B3 projection sign
       R_ref^proj = |<d_ref>| / <|d_ref|> >= 0.1 and resolved (> 3 sigma); B4 effective blocks N / (2 tau_int) >= 8;
       B5 physical sector fractions (no D_r < -3 sigma, none > 1 + 3 sigma); B6 >= 2 seeds compatible (D_ref and C_ref at
       the fixed tau within 2 sigma); B7 stable block errors (jackknife with 4 and 8 blocks within a factor 2).  All met ->
       SPG_RESOLVED; B2 failed -> SPG_INACTIVE; otherwise SPG_PROJECTION_UNRESOLVED (no ground-state statement).
       Unit test on synthetic per-configuration weights (resolved / cancelling).
  SPP8 (statement) the budget verdict of every production run; promotion only with SPG_RESOLVED.
  SPP9 (statement) POOLED budget per (g, beta, lambda) over independent chains (files spg_points_<g><b>[s<k>|c<k>]): chain-level
       delete-one-chain jackknife (chains are independent), N_eff = sum over chains of n / (2 tau_int), B6 as chi^2 / dof of
       the chain means of D_ref (< 2) with matching reference sector, B7 chain jackknife vs within-chain blocks within x2.
       Extended 2026-10-03 BEFORE the ten-chain pool was evaluated: B5 also requires Im D_r compatible with 0 (3 sigma, chain
       jackknife; nothing clipped), B6 also requires the chain values of the ratio C_ref = N / D at tau = 0.5, 1, 1.5 to be
       compatible (chi^2 / dof < 2); a failure is classified as MIXING (R^proj >= 0.5) or PROJECTION_PHASE (R^proj < 0.5).
  SPP6 (statement) promotion test at g = 0.018: SECTOR_LOW_GAP from the multiplicity-free beta-slope of the sector ratios
       ln(D_s / D_ref) between beta = 6 and 8 (Delta E = -d ln(D_s / D_ref) / d beta), SPG_LOW_GAP from the reference-sector
       curve at fixed tau (beta = 6 vs 8) -- with the noise diagnostics (negative D_r, per-configuration |f_U| >> 1).
    python3 -B spg_production.py                    (checks + evaluation of spg_points_*.json)
    python3 -B spg_production.py run <beta> <g> <sweeps> <seed> <tag> <dx,dy,dz> [hook_every]
    python3 -B spg_production.py runlong <beta> <g> <sweeps> <burn> <seed> <tag> <dx,dy,dz> <hook_every> <f_every>   (QMC.MIXING.MODE.01:
        long burn-in, per-sweep field features and every f_every sweeps the 64 sector traces streamed to spgm_monitor_<tag>.jsonl,
        final configuration saved to spgm_phi_<tag>.npy)
    python3 -B spg_production.py env <beta> <g> <sweeps> <seed> <tag> <dx,dy,dz> <alpha> <lambda>   (QMC.SECTOR_SAMPLING.PRODUCTION.01:
        envelope h = sqrt(|d_ref|^2 + alpha^2 + lambda^2 sum_{tau = 0.5, 1, 1.5} |y_ref(tau)|^2), reference sector (0, 0) fixed in advance)
"""
from __future__ import annotations

import importlib.util as ilu
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import expm, logm, lu_factor, lu_solve

HERE = Path(__file__).resolve().parent
_sp = ilu.spec_from_file_location("reweight", HERE / "reweight.py"); rw = ilu.module_from_spec(_sp); _sp.loader.exec_module(rw)
udv_left, inv_sum = rw.Model._udv_left, rw.Model._inv_sum
TAUS = (0.0, 0.5, 1.0, 1.5, 2.0, 3.0)                                  # predeclared fixed physical times (t = 1 units)
LOG: list[dict] = []
wrap = lambda a: (a + np.pi) % (2 * np.pi) - np.pi


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    sys.stdout.flush()


# ------------------------------------------------------------------------------------------------ core (generic in S)
def chains(S, Sinv):
    L = len(S); nn = S[0].shape[0]; I = np.eye(nn, dtype=complex); unit = (I, np.ones(nn), I)
    A = [unit]                                                         # A[l] = UDV of (S_0 .. S_{l-1})^T
    for l in range(L):
        A.append(udv_left(S[l].T, A[-1]))
    C = [unit] * (L + 1)                                               # C[l] = UDV of ((S_l .. S_{L-1})^-1)^T
    for l in range(L - 1, -1, -1):
        C[l] = udv_left(Sinv[l].T, C[l + 1])
    return A, C


def rho_prep(A_L, Fn):
    """C' = Fn^-1 conj(U) / D_b, B' = Fn^-1 V^T D_s with S_X^T = U D V (S_X = V^T D U^T): det(1 + S_u S_X) / det(1 + S_X) =
    det(C' + D_u B') / det(C' + B') for S_u = Fn D_u Fn^-1."""
    U, D, V = A_L
    Db, Ds = np.maximum(D, 1.0), np.minimum(D, 1.0)
    Fi = np.linalg.inv(Fn)
    return Fi @ (U.conj() / Db[None, :]), Fi @ (V.T * Ds[None, :])


def rho_track(Cp, Bp, ang):
    """f = Tr[Gamma(u) X] / Tr[X] for u with eigen-angles ang (n particle modes in the basis F), sign-exact by continuity:
    h(t) = sqrt(R(t)) from h(0) = 1, R(t) = det(C' + e^{i t ang_N} B') / det(C' + B'), f = h(1) e^{i sum(ang) / 2}."""
    dph = np.concatenate([ang, -ang])
    s0, l0 = np.linalg.slogdet(Cp + Bp)
    def logR(t):
        s, l = np.linalg.slogdet(Cp + np.exp(1j * t * dph)[:, None] * Bp)
        return l - l0, np.angle(s) - np.angle(s0)
    t, dt, acc, prev, nstep = 0.0, 0.25, 0.0, 0.0, 0
    lr = 0.0
    while t < 1 - 1e-12:
        dt = min(dt, 1 - t)
        lr_new, ph_new = logR(t + dt); nstep += 1
        inc = wrap(ph_new - prev)
        if abs(inc) > 0.4 and dt > 1e-6:
            dt /= 2; continue
        acc += inc; prev = ph_new; lr = lr_new; t += dt; dt *= 1.5
    return np.exp(0.5 * (lr + 1j * acc) + 0.5j * np.sum(ang)), nstep


def green_prep(A, C, l):
    U1, D1, V1 = C[l]; U2, D2, V2 = A[l]
    return {"U1": U1, "D1": D1, "V1": V1, "U2": U2, "D2": D2, "V2i": np.linalg.inv(V2), "U12": U1.conj().T @ U2}


def green_U(P, SUc, Pi, want_det=False):
    """(1/k) tr Pi^dag G_U Pi with G_U^T = U_C D_C V_C + U_A D_A V_A S_U^T (Loh-Gubernatis split); SUc = conj(S_U) as a
    callable acting on rows (left) / columns (right) or as a dense matrix."""
    Db1, Ds1 = np.maximum(P["D1"], 1.0), np.minimum(P["D1"], 1.0)
    Db2, Ds2 = np.maximum(P["D2"], 1.0), np.minimum(P["D2"], 1.0)
    V1c = SUc("right", P["V1"]) if callable(SUc) else P["V1"] @ SUc   # V1 conj(S_U)
    mid = Ds1[:, None] * (V1c @ P["V2i"]) / Db2[None, :] + P["U12"] * Ds2[None, :] / Db1[:, None]
    lu = lu_factor(mid)
    rhs = (P["U1"].conj().T @ Pi.conj()) / Db1[:, None]
    Y = P["V2i"] @ (lu_solve(lu, rhs) / Db2[:, None])
    Y = SUc("left", Y) if callable(SUc) else SUc @ Y                    # conj(S_U) V2^-1 ...
    val = np.trace(Pi.T @ Y) / Pi.shape[1]
    if want_det:
        d = np.diag(lu[0]); piv = lu[1]
        sgn = (-1) ** int(np.sum(piv != np.arange(len(piv))))
        return val, np.sum(np.log(np.abs(d))) + 1j * (np.sum(np.angle(d)) + (np.pi if sgn < 0 else 0.0))
    return val


# ------------------------------------------------------------------------------------------------ SPP1 brute force
rng = np.random.default_rng(37)
nb = 3
Zp = np.diag([1.0, -1.0]); low = np.array([[0, 1], [0, 0]], complex)
cf = []
for i in range(nb):
    mm = np.eye(1, dtype=complex)
    for j in range(nb):
        mm = np.kron(mm, Zp if j < i else (low if j == i else np.eye(2)))
    cf.append(mm)
Psi = cf + [x.conj().T for x in cf]
def randK(s_=0.7):
    h = rng.normal(size=(nb, nb)) + 1j * rng.normal(size=(nb, nb)); h = (h + h.conj().T) / 2
    D = rng.normal(size=(nb, nb)) + 1j * rng.normal(size=(nb, nb)); D = D - D.T
    return s_ * np.block([[h, D.conj().T], [D, -h.T]])
def fock(K):
    return expm(0.5 * sum(K[i, j] * Psi[i].conj().T @ Psi[j] for i in range(2 * nb) for j in range(2 * nb)))
def fock_unitary(u):
    A_ = -1j * logm(u)
    return expm(1j * sum(A_[i, j] * cf[i].conj().T @ cf[j] for i in range(nb) for j in range(nb)))
def nambu_u(u):
    return np.block([[u, np.zeros_like(u)], [np.zeros_like(u), u.conj()]])
errs1 = []
for trial in range(6):
    L = 5; Ks = [randK() for _ in range(L)]
    S = [expm(K) for K in Ks]; Sinv = [expm(-K) for K in Ks]
    Xf = [fock(K) for K in Ks]
    A, C = chains(S, Sinv)
    Q_, R_ = np.linalg.qr(rng.normal(size=(nb, nb)) + 1j * rng.normal(size=(nb, nb))); u = Q_ * (np.diag(R_) / np.abs(np.diag(R_)))
    if trial % 2 == 0:
        w_, V_ = np.linalg.eig(u); w_[0] = -1.0; Qv, _ = np.linalg.qr(V_); u = Qv @ np.diag(w_ / np.abs(w_)) @ Qv.conj().T
    w_, F_ = np.linalg.eig(u); Fq, _ = np.linalg.qr(F_); ang = np.angle(np.diag(Fq.conj().T @ u @ Fq))
    Fn = np.block([[Fq, np.zeros_like(Fq)], [np.zeros_like(Fq), Fq.conj()]])
    Cp, Bp = rho_prep(A[L], Fn)
    f, _ = rho_track(Cp, Bp, ang)
    Xtot = np.linalg.multi_dot(Xf); Uf = fock_unitary(u)
    f_bf = np.trace(Uf @ Xtot) / np.trace(Xtot)
    Pi = np.zeros((2 * nb, 2), complex); Pi[:nb] = rng.normal(size=(nb, 2)) + 1j * rng.normal(size=(nb, 2))
    SUc = nambu_u(u).conj()
    for l in (1, 3, 5):
        g = green_U(green_prep(A, C, l), SUc, Pi)
        XA = np.linalg.multi_dot([np.eye(2 ** nb), np.eye(2 ** nb)] + Xf[:l]); XC = np.linalg.multi_dot([np.eye(2 ** nb), np.eye(2 ** nb)] + Xf[l:])
        num = 0
        for col in range(2):
            O = sum(np.conj(Pi[i, col]) * cf[i] for i in range(nb))
            num += np.trace(Uf @ XA @ O @ XC @ O.conj().T) / 2
        errs1.append(abs(f * g - num / np.trace(Xtot)) / abs(num / np.trace(Xtot)))
    errs1.append(abs(f - f_bf) / abs(f_bf))
step("SPP1 brute force (3 modes, random BdG slices with pairing, random unitaries, three of six with an eigenvalue -1): the sign-exact path continuation f_U = Tr[Gamma(u) X] / Tr[X] and the stabilised f_U g_U(tau) = Tr[Gamma(u) X_A O X_C O^dag] / Tr[X] equal the Fock-space traces",
     max(errs1) < 1e-8, f"max rel. error {max(errs1):.1e} over 6 configurations x (1 trace + 3 cuts)")


# ------------------------------------------------------------------------------------------------ production operator
def canon_index(m):
    Rv = rw.CELLS[m.cell].astype(float); Rinv = np.linalg.inv(Rv)
    def canon(x):
        return tuple(np.round(np.mod(np.round(x @ Rinv, 9), 1.0), 6))
    keys = {canon(x): i for i, x in enumerate(m.sites)}
    return canon, keys


class Sector:
    """translations u_a = G^-1 (P_a x 1_32) G (twisted cell), Z4, common eigenbasis F, characters, momentum projector Pi."""
    def __init__(self, m, kname="P+d"):
        self.m = m; N, n = m.N, m.n
        canon, keys = canon_index(m)
        X = np.array(m.sites, float)
        self.perm = [np.array([keys[canon(m.sites[x] + m.sites[a])] for x in range(N)]) for a in range(N)]
        self.id_a = [a for a in range(N) if np.array_equal(self.perm[a], np.arange(N))][0]
        qf = np.tile(rw.QF, 2)                                         # index s * 16 + f
        tw = np.zeros(3) if m.twist is None else m.twist
        best = None
        for sg in (1.0, -1.0):
            gvec = np.exp(1j * sg * 2 * (X @ tw)[:, None] * qf[None, :]).reshape(n)
            us = [self._mono(a, gvec) for a in range(N)]
            res = max(np.abs(self._dense(u) @ m.Hw - m.Hw @ self._dense(u)).max() for u in us[:4])
            if best is None or res < best[0]:
                best = (res, sg, gvec, us)
        self.comm_res, self.sigma, self.gvec, self.us = best
        # common eigenbasis: F = G^-1 (plane waves x 1_32); plane-wave momenta of the supercell
        Rv = rw.CELLS[m.cell].astype(float)
        cand = [2 * np.pi * np.linalg.solve(Rv, np.array(nv, float)) for nv in np.ndindex(4, 4, 4)]
        chars, Ks = [], []
        for K in cand:
            ch = np.exp(-1j * X @ K)
            if not any(np.allclose(ch, c2, atol=1e-9) for c2 in chars):
                chars.append(ch); Ks.append(K)
        assert len(Ks) == N
        self.Ks = np.array(Ks)
        Fs = np.exp(1j * X @ self.Ks.T) / np.sqrt(N)                  # (x, K)
        self.F = (np.kron(Fs, np.eye(32)) / self.gvec[:, None])        # G^-1 (Fs x 1)
        self.Fn = np.block([[self.F, np.zeros((n, n))], [np.zeros((n, n)), self.F.conj()]])
        # elements (a, m): eigen-angles in F, characters chi_{K,q}(a, m) = lambda_K(a) i^{q m}
        self.elems, self.ang = [], []
        diag_res = 0.0
        for a in range(N):
            Ud = self._dense(self.us[a])
            Dg = self.F.conj().T @ Ud @ self.F
            diag_res = max(diag_res, np.abs(Dg - np.diag(np.diag(Dg))).max())
            lam = np.diag(Dg)
            for mz in range(4):
                self.elems.append((a, mz)); self.ang.append(np.angle(lam * 1j ** mz))
        self.diag_res = diag_res
        lamK = np.array([[np.exp(1j * np.angle(np.mean(np.diag(self.F.conj().T @ self._dense(self.us[a]) @ self.F)[32 * k:32 * k + 32]))) for a in range(N)] for k in range(N)])
        self.lamK = lamK                                                # (K, a): eigenvalue of u_a on plane wave K
        self.sectors = [(k, q) for k in range(N) for q in range(4)]
        self.chi = np.array([[lamK[k, a] * 1j ** (q * mz) for (a, mz) in self.elems] for (k, q) in self.sectors])
        # momentum projector of the measured mode (greens_k convention): Q = +1/2 flavours at k + delta
        ks = m.momenta(); base = kname[:-2]
        v = ks[base] + (tw if kname.endswith("+d") else -tw)
        ph = np.exp(-1j * X @ v)
        qsel = np.where(rw.QF == (0.5 if kname.endswith("+d") else -0.5))[0]
        cols = []
        for s_ in range(2):
            for f in qsel:
                c = np.zeros(2 * n, complex)
                c[np.arange(N) * 32 + s_ * 16 + f] = ph.conj() / np.sqrt(N)
                cols.append(c)
        self.Pi = np.array(cols).T

    def _mono(self, a, gvec):
        n = self.m.n
        src = np.arange(n); x = src // 32; c = src % 32
        dst = self.perm[a][x] * 32 + c
        return (dst, gvec[src] / gvec[dst])                            # u[dst, src] = phase

    def _dense(self, u):
        dst, ph = u; n = self.m.n
        M = np.zeros((n, n), complex); M[dst, np.arange(n)] = ph
        return M

    def SUc(self, a, mz):
        """conj(S_U) for U = i^m u_a, S_U = diag(u, u^*) (monomial), as a callable for green_U."""
        dst, ph = self.us[a]; n = self.m.n
        z = 1j ** mz
        dstN = np.concatenate([dst, dst + n]); phN = np.concatenate([ph * z, np.conj(ph * z)])
        cph = np.conj(phN)                                             # conj(S_U)[dstN[j], j] = cph[j]
        def op(side, Mx):
            if side == "right":                                        # Mx conj(S_U): column j <- Mx[:, dstN[j]] cph[j]
                return Mx[:, dstN] * cph[None, :]
            out = np.empty_like(Mx); out[dstN] = cph[:, None] * Mx    # conj(S_U) Mx: row dstN[j] <- cph[j] Mx[j]
            return out
        return op


def slices(m, phi):
    if not hasattr(m, "Wh_inv"):
        m.Wh_inv = np.linalg.inv(m.Wh)
    S = [m.slice_op(phi[l]) for l in range(m.Lt)]
    Sinv = [m.Wh_inv @ (m.Wh_inv @ m.slice_op(-phi[l]) @ m.Wh_inv) @ m.Wh_inv for l in range(m.Lt)]
    return S, Sinv


def measure(sec, phi):
    m = sec.m; L = m.Lt
    S, Sinv = slices(m, phi)
    A, C = chains(S, Sinv)
    Cp, Bp = rho_prep(A[L], sec.Fn)
    f = np.zeros(len(sec.elems), complex); nsteps = 0
    for e_, ang in enumerate(sec.ang):
        f[e_], ns = rho_track(Cp, Bp, ang); nsteps += ns
    ls = [L - int(round(t / m.dtau)) for t in TAUS if t <= m.beta / 2 + 1e-9]
    g = np.zeros((len(sec.elems), len(ls)), complex)
    for j, l in enumerate(ls):
        P = green_prep(A, C, l)
        for e_, (a, mz) in enumerate(sec.elems):
            g[e_, j] = green_U(P, sec.SUc(a, mz), sec.Pi)
    return {"f": f, "g": g, "taus": [t for t in TAUS if t <= m.beta / 2 + 1e-9], "path_steps": nsteps}


# ------------------------------------------------------------------------------------------------ CLI: production run
if len(sys.argv) > 1 and sys.argv[1] == "run":
    beta, g_, sweeps, seed, tag = float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
    twist = np.array([float(x) for x in sys.argv[7].split(",")])
    every = int(sys.argv[8]) if len(sys.argv) > 8 else 2
    holder = {}
    def hook(m, phi):
        if "sec" not in holder:
            holder["sec"] = Sector(m)
        t0 = time.time(); r = measure(holder["sec"], phi); r["seconds"] = time.time() - t0
        return r
    res = rw.run("cubic16", beta, 0.25, g_, sweeps, seed, "md", True, twist, None, hook=hook, hook_every=every)
    sec = holder["sec"]
    recs = res.pop("hook")
    out = {k: res[k] for k in ("cell", "beta", "dtau", "g", "sweeps", "seed", "R", "R_err", "acceptance", "seconds")}
    out["twist_delta"] = twist.tolist(); out["k"] = "P+d"; out["taus"] = recs[0][2]["taus"]
    out["G_tau_greens"] = res["greens"]["P+d"]["G_tau"]
    out["theta"] = [r_[0] for r_ in recs]
    out["f_re"] = [r_[2]["f"].real.tolist() for r_ in recs]; out["f_im"] = [r_[2]["f"].imag.tolist() for r_ in recs]
    out["g_re"] = [r_[2]["g"].real.tolist() for r_ in recs]; out["g_im"] = [r_[2]["g"].imag.tolist() for r_ in recs]
    out["G_greens_at_taus"] = [[float(np.real(r_[1]["P+d"][int(round(t / 0.25))])) for t in recs[0][2]["taus"]] for r_ in recs]
    out["identity_element"] = sec.elems.index((sec.id_a, 0))
    out["path_steps"] = [r_[2]["path_steps"] for r_ in recs]; out["meas_seconds"] = [r_[2]["seconds"] for r_ in recs]
    out["sectors"] = [[int(k), int(q)] for k, q in sec.sectors]; out["K"] = sec.Ks.tolist()
    out["chi_re"] = sec.chi.real.tolist(); out["chi_im"] = sec.chi.imag.tolist()
    (HERE / f"spg_points_{tag}.json").write_text(json.dumps(out, default=float))
    print(json.dumps({k: out[k] for k in ("beta", "g", "R", "seconds")}, default=float))
    sys.exit(0)


def field_features(sec, phi):
    """global HS-field observables: uniform channel components, tau-Fourier modes of the spatial mean, |phi_A(K, omega = 0)|^2
    at all 16 cell momenta, total norm."""
    m = sec.m; X = np.array(m.sites, float)
    mA = phi.mean(axis=(0, 1))
    mt = phi.mean(axis=1)                                              # (L, 10)
    Ft = np.fft.rfft(mt, axis=0)[:3] / mt.shape[0]
    Phi = phi.sum(axis=0)                                              # (N, 10)
    ph = np.exp(-1j * sec.Ks @ X.T)                                    # (16, N)
    PK = np.abs(ph @ Phi) ** 2 / (m.N * m.Lt)                          # (16, 10)
    stag = (phi * ((-1.0) ** np.arange(phi.shape[0]))[:, None, None]).mean(axis=0)   # alternating imaginary-time mode (N, 10)
    return {"mA": mA.tolist(), "Ft_abs": np.abs(Ft).tolist(), "PK": PK.tolist(), "norm2": float(np.sum(phi ** 2)),
            "stagA": stag.mean(axis=0).tolist(), "stagK": (np.abs(ph @ stag) ** 2 / m.N).sum(axis=1).tolist()}


def f_only(sec, phi):
    m = sec.m
    S = [m.slice_op(phi[l]) for l in range(m.Lt)]
    I = np.eye(2 * m.n, dtype=complex); A = (I, np.ones(2 * m.n), I)
    for l in range(m.Lt):
        A = udv_left(S[l].T, A)
    Cp, Bp = rho_prep(A, sec.Fn)
    return np.array([rho_track(Cp, Bp, ang)[0] for ang in sec.ang])


if len(sys.argv) > 1 and sys.argv[1] == "runlong":
    beta, g_, sweeps, burn, seed, tag = float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), sys.argv[7]
    twist = np.array([float(x) for x in sys.argv[8].split(",")]); every, f_every = int(sys.argv[9]), int(sys.argv[10])
    phi_start = np.load(sys.argv[11]) if len(sys.argv) > 11 else None    # continuation from a saved configuration
    holder = {}
    mon_file = HERE / f"spgm_monitor_{tag}.jsonl"; mon_file.write_text("")
    def ensure(m):
        if "sec" not in holder:
            holder["sec"] = Sector(m); holder["ref"] = holder["sec"].sectors.index((0, 0))
        return holder["sec"]
    def monitor(m, phi, state):
        sec = ensure(m)
        rec = {"sw": state["sw"], "theta": state["theta"], "logabsW": state["logabsW"], **field_features(sec, phi)}
        if state["sw"] % f_every == 0:
            f = f_only(sec, phi)
            rec["f_re"] = f.real.tolist(); rec["f_im"] = f.imag.tolist()
            rec["d_ref"] = complex(np.exp(1j * state["theta"]) * (f @ sec.chi[holder["ref"]].conj()) / 64).real
            rec["kappa_char"] = float(np.sum(np.abs(f)) / max(abs(f @ sec.chi[holder["ref"]].conj()), 1e-300))
        with mon_file.open("a") as fh:
            fh.write(json.dumps(rec) + "\n")
        return {"sw": state["sw"]}
    def hook(m, phi):
        sec = ensure(m); t0 = time.time(); r = measure(sec, phi); r["seconds"] = time.time() - t0
        return r
    res = rw.run("cubic16", beta, 0.25, g_, sweeps, seed, "md", True, twist, None, hook=hook, hook_every=every,
                 burn=burn, monitor=monitor, monitor_every=1, save_phi=str(HERE / f"spgm_phi_{tag}.npy"), phi0=phi_start)
    sec = holder["sec"]; recs = res.pop("hook"); res.pop("monitor", None)
    out = {k: res[k] for k in ("cell", "beta", "dtau", "g", "sweeps", "seed", "R", "R_err", "acceptance", "seconds", "burn")}
    out["continued_from"] = sys.argv[11] if len(sys.argv) > 11 else None
    out["twist_delta"] = twist.tolist(); out["k"] = "P+d"; out["taus"] = recs[0][2]["taus"]
    out["theta"] = [r_[0] for r_ in recs]
    out["f_re"] = [r_[2]["f"].real.tolist() for r_ in recs]; out["f_im"] = [r_[2]["f"].imag.tolist() for r_ in recs]
    out["g_re"] = [r_[2]["g"].real.tolist() for r_ in recs]; out["g_im"] = [r_[2]["g"].imag.tolist() for r_ in recs]
    out["G_greens_at_taus"] = [[float(np.real(r_[1]["P+d"][int(round(t / 0.25))])) for t in recs[0][2]["taus"]] for r_ in recs]
    out["identity_element"] = sec.elems.index((sec.id_a, 0)); out["meas_seconds"] = [r_[2]["seconds"] for r_ in recs]
    out["sectors"] = [[int(k), int(q)] for k, q in sec.sectors]; out["K"] = sec.Ks.tolist()
    out["chi_re"] = sec.chi.real.tolist(); out["chi_im"] = sec.chi.imag.tolist()
    (HERE / f"spgl_points_{tag}.json").write_text(json.dumps(out, default=float))
    print(json.dumps({k: out[k] for k in ("beta", "g", "R", "burn", "seconds")}, default=float))
    sys.exit(0)


if len(sys.argv) > 1 and sys.argv[1] == "env":
    beta, g_, sweeps, seed, tag = float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
    twist = np.array([float(x) for x in sys.argv[7].split(",")]); alpha, lam_e = float(sys.argv[8]), float(sys.argv[9])
    holder = {}
    def outer(m, phi):
        if "sec" not in holder:
            holder["sec"] = Sector(m); holder["ref"] = holder["sec"].sectors.index((0, 0))
        sec = holder["sec"]; t0 = time.time(); r = measure(sec, phi); r["seconds"] = time.time() - t0
        cref = sec.chi[holder["ref"]].conj()
        dref = abs(r["f"] @ cref) / 64
        jt = [j for j, t in enumerate(r["taus"]) if t in (0.5, 1.0, 1.5)]
        y = [abs((r["f"] * r["g"][:, j]) @ cref) / 64 for j in jt]
        return float(np.sqrt(dref ** 2 + alpha ** 2 + lam_e ** 2 * sum(v ** 2 for v in y))), r
    res = rw.run("cubic16", beta, 0.25, g_, sweeps, seed, "md", False, twist, None, outer=outer)
    sec = holder["sec"]; recs = res.pop("outer")
    out = {k: res[k] for k in ("cell", "beta", "dtau", "g", "sweeps", "seed", "R", "R_err", "acceptance", "seconds")}
    out.update({"sampler": "envelope", "alpha": alpha, "lambda": lam_e, "outer_acceptance": res["outer_acceptance"], "twist_delta": twist.tolist(), "k": "P+d",
                "taus": recs[0][2]["taus"], "theta": [r_[0] for r_ in recs], "h": [r_[1] for r_ in recs],
                "f_re": [r_[2]["f"].real.tolist() for r_ in recs], "f_im": [r_[2]["f"].imag.tolist() for r_ in recs],
                "g_re": [r_[2]["g"].real.tolist() for r_ in recs], "g_im": [r_[2]["g"].imag.tolist() for r_ in recs],
                "meas_seconds": [r_[2]["seconds"] for r_ in recs], "identity_element": sec.elems.index((sec.id_a, 0)),
                "sectors": [[int(k), int(q)] for k, q in sec.sectors], "K": sec.Ks.tolist(), "chi_re": sec.chi.real.tolist(), "chi_im": sec.chi.imag.tolist()})
    (HERE / f"spge_points_{tag}.json").write_text(json.dumps(out, default=float))
    print(json.dumps({k: out[k] for k in ("beta", "g", "R", "outer_acceptance", "seconds")}, default=float))
    sys.exit(0)


# ------------------------------------------------------------------------------------------------ SPP2 free closed forms
m0 = rw.Model("cubic16", 1.0, 0.25, 0.0, "md", np.array([0.39, 0.39, 0.39]))
sec0 = Sector(m0)
phi0 = np.zeros((m0.Lt, m0.N, 10))
r0 = measure(sec0, phi0)
Eb = expm(-m0.beta * m0.Hw); I_n = np.eye(m0.n)
Pn = sec0.Pi[:m0.n]
err_f, err_g = 0.0, 0.0
for e_, (a, mz) in enumerate(sec0.elems):
    u = sec0._dense(sec0.us[a]) * 1j ** mz
    f_ex = np.linalg.det(I_n + u @ Eb) / np.linalg.det(I_n + Eb)
    err_f = max(err_f, abs(r0["f"][e_] - f_ex) / max(abs(f_ex), 1e-300))
    for j, t in enumerate(r0["taus"]):
        g_ex = np.trace(Pn.conj().T @ expm(-t * m0.Hw) @ np.linalg.solve(I_n + u @ Eb, Pn)) / Pn.shape[1]
        err_g = max(err_g, abs(r0["g"][e_, j] - g_ex))                     # absolute (g = O(1/2); some exact values are ~0)
step("SPP2 free limit of the production operator (twisted cubic16, lambda = 0.39, beta = 1, g = 0): f_U and g_U(tau) from the production path for all 64 elements (16 twisted translations x Z4) equal the closed forms det(1 + u e^{-beta h}) / det(1 + e^{-beta h}) and tr Pi^dag e^{-tau h} (1 + u e^{-beta h})^-1 Pi / 16",
     err_f < 1e-8 and err_g < 1e-8, f"max rel. error f {err_f:.1e}, max abs. error g {err_g:.1e}; path steps per configuration {r0['path_steps']}")

# ------------------------------------------------------------------------------------------------ SPP3 symmetry of the twisted operator
m1 = rw.Model("cubic16", 1.0, 0.25, 0.018, "md", np.array([0.39, 0.39, 0.39]))
sec1 = Sector(m1)
rep_res = 0.0
for a in range(4):
    for b in range(4):
        cidx = [c for c in range(m1.N) if np.array_equal(sec1.perm[c], sec1.perm[a][sec1.perm[b]])][0]
        Uab = sec1._dense(sec1.us[a]) @ sec1._dense(sec1.us[b])
        rep_res = max(rep_res, np.abs(Uab - sec1._dense(sec1.us[cidx])).max())
cov_res, orth_res = 0.0, 0.0
rng3 = np.random.default_rng(5)
for a in (1, 2, 5):                                                    # three translations
    u = sec1._dense(sec1.us[a])
    for i in (0, 3, 7):                                                 # three sites
        Bs_i = [m1.Ftw[i].T @ rw.MA[A_] @ m1.Ftw[i] for A_ in range(10)]
        j = sec1.perm[a][i]
        Bs_j = np.array([(m1.Ftw[j].T @ rw.MA[A_] @ m1.Ftw[j]).ravel() for A_ in range(10)]).T
        Rm = np.zeros((10, 10), complex)
        for A_ in range(10):
            Y = (u.conj() @ Bs_i[A_] @ u.conj().T).ravel()             # psi -> u psi: psi^T B psi -> psi^T u^T B u psi; B -> u^-T B u^-1
            coef, *_ = np.linalg.lstsq(Bs_j, Y, rcond=None)
            cov_res = max(cov_res, np.abs(Bs_j @ coef - Y).max() / np.abs(Y).max())
            Rm[A_] = coef
        orth_res = max(orth_res, np.abs(Rm.conj().T @ Rm - np.eye(10)).max(), np.abs(Rm.imag).max())
step("SPP3 symmetry of the twisted production operator (cubic16, lambda = 0.39): the twisted translations u_a = G^-1 (P_a x 1) G commute with H_w, form a representation, are diagonal in the common plane-wave basis F, and map every HS composite B_{iA} into span{B_{pi(i) B}} with a real orthogonal 10 x 10 matrix (a Spin(10) Cartan rotation) -- with the isotropic Gaussian prior the HS-averaged slice commutes with all 64 sector elements (Z4: as in ST4)",
     sec1.comm_res < 1e-12 and rep_res < 1e-12 and sec1.diag_res < 1e-10 and cov_res < 1e-10 and orth_res < 1e-10,
     f"[u_a, H_w] {sec1.comm_res:.1e} (gauge sign {sec1.sigma:+.0f}), representation {rep_res:.1e}, diagonal in F {sec1.diag_res:.1e}, composite covariance {cov_res:.1e}, orthogonality / reality {orth_res:.1e}")

# ------------------------------------------------------------------------------------------------ SPP4 interacting configuration
rng4 = np.random.default_rng(8)
phi1 = np.sqrt(2 * 0.018 * 0.25) * rng4.normal(size=(m1.Lt, m1.N, 10)) * 3.0
r1 = measure(sec1, phi1)
Gk = m1.greens_k(phi1, m1.momenta())["P+d"]
idx_t = [int(round(t / m1.dtau)) for t in r1["taus"]]
e_id = sec1.elems.index((sec1.id_a, 0))
err_id = max(abs(r1["g"][e_id, j] - Gk[jt]) for j, jt in enumerate(idx_t)) + abs(r1["f"][e_id] - 1)
S1, Sinv1 = slices(m1, phi1); A1, C1 = chains(S1, Sinv1)
err_det = 0.0
for l in (1, 2, 4):
    P = green_prep(A1, C1, l)
    _, d1 = green_U(P, sec1.SUc(sec1.id_a, 0), sec1.Pi, want_det=True)
    for e_ in (1, 5, 22, 47, 63):
        a, mz = sec1.elems[e_]
        _, dU = green_U(P, sec1.SUc(a, mz), sec1.Pi, want_det=True)
        logf2 = 1j * np.sum(sec1.ang[e_]) + dU - d1                        # det(u) det(mid_U) / det(mid_1)
        err_det = max(err_det, abs(np.exp(logf2) - r1["f"][e_] ** 2) / abs(r1["f"][e_] ** 2))
step("SPP4 interacting configuration (g = 0.018, lambda = 0.39, beta = 1, HS field 3x the prior width): the U = 1 term reproduces greens_k at every measured tau, and f_U^2 from the sign-exact path equals the independent determinant ratio det(u) det(mid_U) / det(mid_1) at three cuts l (stabilisation and the det(u) factor)",
     err_id < 1e-10 and err_det < 1e-7, f"U = 1 vs greens_k {err_id:.1e}; f^2 vs det ratio {err_det:.1e}; |f_U| range {np.abs(r1['f']).min():.2e}..{np.abs(r1['f']).max():.2f}; path steps {r1['path_steps']}")


# ------------------------------------------------------------------------------------------------ SPP5 evaluation
def evaluate(d):
    th = np.array(d["theta"]); z = np.exp(1j * th)
    f = np.array(d["f_re"]) + 1j * np.array(d["f_im"]); g = np.array(d["g_re"]) + 1j * np.array(d["g_im"])
    chi = np.array(d["chi_re"]) + 1j * np.array(d["chi_im"])
    nm = len(z)
    def est(sel):
        zz = z[sel]
        Dr = (zz[:, None] * (f[sel] @ chi.conj().T / chi.shape[1])).mean(0) / zz.mean()                  # (n_sec,)
        Nr = np.einsum("m,me,met,se->st", zz, f[sel], g[sel], chi.conj()) / chi.shape[1] / len(zz) / zz.mean()
        return Dr, Nr
    Dr, Nr = est(np.arange(nm))
    nbk = min(nm, 8); blocks = np.array_split(np.arange(nm), nbk)
    jk = [est(np.setdiff1d(np.arange(nm), b)) for b in blocks]
    Dr_e = np.sqrt((nbk - 1) * np.var([j_[0].real for j_ in jk], axis=0))
    ref = int(np.argmax(Dr.real))
    Cr = Nr[ref] / Dr[ref]
    Cr_jk = np.array([j_[1][ref] / j_[0][ref] for j_ in jk])
    Cr_e = np.sqrt((nbk - 1) * np.var(Cr_jk.real, axis=0))
    tot = Nr.sum(0) / Dr.sum()
    return {"D": Dr, "D_err": Dr_e, "ref": ref, "C_ref": Cr, "C_ref_err": Cr_e, "C_total": tot, "n": nm}


def free_reference(beta, twist, taus):
    """free (g = 0) sector fractions D_r and the projected C_r(tau) of the largest free sector, same cell / twist / beta -- exact
    product formulas in the common eigenbasis (every u is a phase on each 32 x 32 plane-wave block of H_w; no large determinants)."""
    mf = rw.Model("cubic16", beta, 0.25, 0.0, "md", np.asarray(twist))
    sf = Sector(mf)
    Hb = sf.F.conj().T @ mf.Hw @ sf.F
    PiF = sf.F.conj().T @ sf.Pi[:mf.n]                                 # measured modes in the eigenbasis
    blocks = []
    for k in range(mf.N):
        e, V = np.linalg.eigh(Hb[32 * k:32 * k + 32, 32 * k:32 * k + 32]); blocks.append((e, V))
    fU, gU = [], []
    for (a, mz) in sf.elems:
        lf = 0.0; g_t = np.zeros(len(taus), complex)
        for k, (e, V) in enumerate(blocks):
            lam = sf.lamK[k, a] * 1j ** mz
            lf += np.sum(np.log1p(lam * np.exp(-beta * e)) - np.log1p(np.exp(-beta * e)))
            Pk = V.conj().T @ PiF[32 * k:32 * k + 32]
            if np.abs(Pk).max() > 1e-12:
                w = np.sum(np.abs(Pk) ** 2, axis=1)
                g_t += np.array([np.sum(w * np.exp(-t * e) / (1 + lam * np.exp(-beta * e))) for t in taus]) / sf.Pi.shape[1]
        fU.append(np.exp(lf)); gU.append(g_t)
    fU, gU = np.array(fU), np.array(gU)
    D = (sf.chi.conj() @ fU) / 64
    ref = int(np.argmax(D.real))
    C = (sf.chi.conj()[ref] @ (fU[:, None] * gU)) / 64 / D[ref]
    return D.real, ref, C.real


_mn = (HERE / "mirror_node.py").read_text()
exec(compile("def tau_x(m):" + _mn.split("def tau_x(m):")[1].split("# ---------------------------------------------------------------------------------------- MN1 brute force")[0], "hfb", "exec"))
Kloc = [np.block([[np.zeros((32, 32)), 2 * MA_.conj().T], [2 * MA_, np.zeros((32, 32))]]) for MA_ in rw.MA]


def hf_reference(beta, twist, g):
    """first-order HFB quasi-particle model (local Sigma of mirror_node.py from the free local density, filtered with the
    twisted filter): number of modes with |e| < 0.1, max |<U>| (U != 1), max D_r and the number of sectors with D_r > 0.02."""
    mu = rw.Model("cubic16", beta, 0.25, 0.0, "md")
    E_c, U_c = np.linalg.eigh(mu.Hw); Gc = (U_c * (1 / (1 + np.exp(-beta * E_c)))) @ U_c.conj().T
    Fc = np.kron(mu.F, np.eye(32))
    Gr_ = (Fc @ Gc @ Fc.T)[:32, :32]; Gh_ = (Fc @ (np.eye(mu.n) - Gc).T @ Fc.T)[:32, :32]
    Sc = hfb_sigma(Kloc, np.block([[Gr_, np.zeros((32, 32))], [np.zeros((32, 32)), Gh_]]), 1.0)[:32, :32]
    mt = rw.Model("cubic16", beta, 0.25, 0.0, "md", np.asarray(twist)); st = Sector(mt)
    H1 = mt.Hw + g * sum(mt.Ftw[i].conj().T @ Sc @ mt.Ftw[i] for i in range(mt.N))
    e1 = np.linalg.eigvalsh(H1); Eb = expm(-beta * H1); I_n = np.eye(mt.n)
    s0, l0 = np.linalg.slogdet(I_n + Eb); fs = []
    for (a, mz) in st.elems:
        s1, l1 = np.linalg.slogdet(I_n + st._dense(st.us[a]) * 1j ** mz @ Eb); fs.append(s1 / s0 * np.exp(l1 - l0))
    fs = np.array(fs); idn = st.elems.index((st.id_a, 0)); D = (st.chi.conj() @ fs).real / 64
    return int((np.abs(e1) < 0.1).sum()), float(np.abs(np.delete(fs, idn)).max()), float(D.max()), int((D > 0.02).sum())


rows = []
for fp in sorted(HERE.glob("spg_points_*.json")):
    d = json.loads(fp.read_text())
    ev = evaluate(d)
    Df, reff, Cf = free_reference(d["beta"], d["twist_delta"], d["taus"])
    hf = hf_reference(d["beta"], d["twist_delta"], d["g"])
    f_ = np.array(d["f_re"]) + 1j * np.array(d["f_im"]); fmax = float(np.abs(np.delete(f_, d["identity_element"], axis=1)).max())
    sec_lab = d["sectors"][ev["ref"]]
    top = np.argsort(-ev["D"].real)[:4]
    rows.append(f"{fp.name[11:-5]} g = {d['g']} beta = {d['beta']:g} R = {d['R']:.3f} n = {ev['n']}: top sectors (K, q): "
                + ", ".join(f"{tuple(d['sectors'][s_])} D = {ev['D'][s_].real:.3f}+-{ev['D_err'][s_]:.3f}" for s_ in top)
                + f"; reference {tuple(sec_lab)}: C_r(tau) at tau = {d['taus']}: " + ", ".join(f"{c.real:.3f}+-{e:.3f}" for c, e in zip(ev["C_ref"], ev["C_ref_err"]))
                + f"; unprojected (sum over sectors) " + ", ".join(f"{c.real:.3f}" for c in ev["C_total"])
                + f"; sectors with D > 0.02: {int((ev['D'].real > 0.02).sum())} | FREE reference (g = 0, same beta, twist): max D = {Df.max():.3f} at {tuple(d['sectors'][reff])}, sectors with D > 0.02: {int((Df > 0.02).sum())}, C_r(tau) = " + ", ".join(f"{c:.3f}" for c in Cf)
                + f" | QMC per-configuration max |f_U| (U != 1) {fmax:.2e} | first-order HF quasi-particles (g = {d['g']}): {hf[0]} modes with |e| < 0.1, max |<U>| {hf[1]:.2e}, max D {hf[2]:.4f}, sectors with D > 0.02: {hf[3]}")
step("SPP5 (statement) production at the GROUND_CONNECTED_3BETA_CANDIDATE branch (cubic16, Q+ at P + lambda (1,1,1)): sector fractions D_r = Z_r / Z (direct sector traces), reference sector (largest D_r), absolute start-sector projected C_r(tau) at the predeclared fixed tau, the unprojected sum (= ordinary G), the free (g = 0) closed-form reference at the same beta and twist (the sector distribution concentrates in the free system only for beta >~ 6-8: lowest free |eps| = 0.53) and the first-order HF quasi-particle model (16 mirror modes moved to |e| ~ 0.01 at g = 0.018) -- the direct sector traces test whether near-zero modes exist, independent of the operator",
     True, " | ".join(rows) if rows else "production runs not finished yet")

# ------------------------------------------------------------------------------------------------ SPP7 / SPP8 budget test
_sb = (HERE / "spectral_bounds.py").read_text()
exec(compile(_sb.split("# ---------------------------------------------------------------------------------------------- Gamma method")[1].split("# SB1 ---")[0], "gamma", "exec"))


def budget_stats(z, d_ref, D_all, D_err):
    """per-run budget quantities from the per-configuration projected weights d_ref = e^{i theta} sum chi^* f / 64."""
    n = len(z)
    R_phase = abs(z.mean())
    Rp = abs(d_ref.mean()) / np.mean(np.abs(d_ref))
    jk = []
    for nbk in (4, 8):
        bl = np.array_split(np.arange(n), nbk)
        vals = [abs(np.delete(d_ref, b).mean()) / np.mean(np.abs(np.delete(d_ref, b))) for b in bl]
        jk.append(np.sqrt((nbk - 1) * np.var(vals)))
    _, eD, tau = uwerr(np.vstack([d_ref.real, z.real]), f=lambda m: m[0] / m[1])
    return {"R_phase": R_phase, "R_proj": Rp, "R_proj_err": jk[1], "jk4": jk[0], "jk8": jk[1], "tau_int": tau, "N_eff": n / (2 * tau),
            "D_ref_gamma_err": eD, "neg": int(np.sum(D_all < -3 * D_err)), "over": int(np.sum(D_all > 1 + 3 * D_err))}
def imag_diag(z, f_, chi):
    """diagnostic (not a criterion): the exact sector fractions are real; max |Im D_r| shows the true noise level."""
    return float(np.abs((z[:, None] * (f_ @ chi.conj().T) / chi.shape[1]).mean(0).imag / abs(z.mean())).max())


def budget_verdict(st, D_ref, D_ref_err, seeds_ok):
    fails = []
    if st["R_phase"] <= 0.5: fails.append("B1")
    sig = max(D_ref_err, st["D_ref_gamma_err"])
    # B2 sharpened after the first application (2026-10-03): 'active' also needs materiality D_ref >= 2/64 (g018b4 had
    # D_ref - 1/64 = 2e-4 at > 3 sigma); a large error (sigma >= 1/64) is noise, not inactivity
    if sig >= 1 / 64: fails.append("B2-noise")
    elif not (D_ref >= 2 / 64 and D_ref - 1 / 64 > 3 * sig): return "SPG_INACTIVE", ["B2"]
    if not (st["R_proj"] >= 0.1 and st["R_proj"] > 3 * st["R_proj_err"]): fails.append("B3")
    if st["N_eff"] < 8: fails.append("B4")
    if st["neg"] or st["over"]: fails.append("B5")
    if seeds_ok is not True: fails.append("B6" + ("" if seeds_ok is False else " (one seed)"))
    if not (0.5 <= (st["jk4"] + 1e-12) / (st["jk8"] + 1e-12) <= 2): fails.append("B7")
    return ("SPG_RESOLVED" if not fails else "SPG_PROJECTION_UNRESOLVED"), fails


# SPP7 unit test: synthetic per-configuration weights
rngb = np.random.default_rng(71)
zs = np.ones(200, complex)
good = 0.2 + 0.05 * rngb.normal(size=200)                                # resolved projected weight
bad = 0.02 + np.exp(2j * np.pi * rngb.random(200))                        # tiny mean on top of unit-modulus cancelling terms
stg = budget_stats(zs, good.astype(complex), np.array([0.2]), np.array([0.01]))
stb = budget_stats(zs, bad, np.array([0.02]), np.array([0.07]))
vg, fg = budget_verdict(stg, 0.2, 0.004, True); vb, fb = budget_verdict(stb, abs(bad.mean()), 0.07, True)
step("SPP7 SPG budget test (SPG.BUDGET.01, B1-B7 fixed before application): synthetic resolved per-configuration weights pass (SPG_RESOLVED); a tiny mean on unit-modulus cancelling weights (R^proj ~ 0.02) is not interpreted (SPG_INACTIVE / SPG_PROJECTION_UNRESOLVED)",
     vg == "SPG_RESOLVED" and vb != "SPG_RESOLVED", f"resolved: R^proj = {stg['R_proj']:.2f}, tau_int {stg['tau_int']:.2f} -> {vg}; cancelling: R^proj = {stb['R_proj']:.3f} +- {stb['R_proj_err']:.3f} -> {vb} {fb}")

runs = {}
for fp in sorted(HERE.glob("spg_points_*.json")) + sorted(HERE.glob("spgl_points_*.json")):
    d = json.loads(fp.read_text()); ev = evaluate(d)
    z = np.exp(1j * np.array(d["theta"])); f_ = np.array(d["f_re"]) + 1j * np.array(d["f_im"])
    chi = np.array(d["chi_re"]) + 1j * np.array(d["chi_im"])
    d_ref = z * (f_ @ chi[ev["ref"]].conj()) / chi.shape[1]
    st = budget_stats(z, d_ref, ev["D"].real, ev["D_err"])
    st["imag_max"] = imag_diag(z, f_, chi)
    dr_all = (z[:, None] * (f_ @ chi.conj().T) / chi.shape[1]).real
    st["frozen_neg"] = int(np.sum(np.all(dr_all < 0, axis=0)))           # sectors negative in EVERY configuration of the run
    runs[fp.name[11:-5] if fp.name.startswith("spg_points_") else fp.name[12:-5]] = (d, ev, st)
def _key(name):
    return name.split("s")[0] if "s" in name[5:] else name                 # g018b6s2 -> g018b6
groups = {}
for name in runs:
    groups.setdefault(name[:7] if len(name) > 6 and name[6] in "LM" else name[:6], []).append(name)
lines = []
verdicts = {}
for gname, members in sorted(groups.items()):
    seeds_ok = None
    if len(members) >= 2:
        (dA, eA, _), (dB, eB, _) = runs[members[0]], runs[members[1]]
        same_ref = eA["ref"] == eB["ref"]
        okD = same_ref and abs(eA["D"][eA["ref"]].real - eB["D"][eB["ref"]].real) < 2 * np.hypot(eA["D_err"][eA["ref"]], eB["D_err"][eB["ref"]])
        okC = same_ref and all(abs(a.real - b.real) < 2 * np.hypot(ea, eb) for a, b, ea, eb in zip(eA["C_ref"], eB["C_ref"], eA["C_ref_err"], eB["C_ref_err"]))
        seeds_ok = bool(okD and okC)
    for name in members:
        d, ev, st = runs[name]
        v, fl = budget_verdict(st, ev["D"][ev["ref"]].real, ev["D_err"][ev["ref"]], seeds_ok)
        verdicts[name] = (v, fl)
        lines.append(f"{name} (g = {d['g']}, beta = {d['beta']:g}): {v} {fl if fl else ''} -- R_phase {st['R_phase']:.3f}, D_ref {ev['D'][ev['ref']].real:.3f} (1/64 = 0.016), R^proj {st['R_proj']:.3f} +- {st['R_proj_err']:.3f}, tau_int {st['tau_int']:.2f}, N_eff {st['N_eff']:.1f}, "
                     f"negative / >1 sectors beyond 3 sigma {st['neg']} / {st['over']} (sectors negative in every configuration: {st['frozen_neg']}), max |Im D_r| {st['imag_max']:.3f}, jackknife 4 / 8 blocks {st['jk4']:.3f} / {st['jk8']:.3f}, seeds {'compatible' if seeds_ok else ('incompatible' if seeds_ok is False else 'one')}; C_ref(tau = 0.5, 1, 1.5) = "
                     + ", ".join(f"{c.real:.3f}" for t, c in zip(d["taus"], ev["C_ref"]) if t in (0.5, 1.0, 1.5)))
step("SPP8 (statement) SPG budget verdict of every production run (promotion of a GROUND_CONNECTED candidate only with SPG_RESOLVED at two temperatures)",
     True, " | ".join(lines) if lines else "no runs")

# SPP9 pooled budget over chains -------------------------------------------------------------------------------------
def pooled(names):
    """chain-pooled sector fractions, reference-sector curve and budget quantities."""
    per = []
    for nm in names:
        d, ev, st = runs[nm]
        z = np.exp(1j * np.array(d["theta"])); f_ = np.array(d["f_re"]) + 1j * np.array(d["f_im"]); g_ = np.array(d["g_re"]) + 1j * np.array(d["g_im"])
        chi = np.array(d["chi_re"]) + 1j * np.array(d["chi_im"])
        per.append((z, f_, g_, chi, st, d))
    def est(sel):
        Z = np.concatenate([per[i][0] for i in sel]); F = np.concatenate([per[i][1] for i in sel]); G = np.concatenate([per[i][2] for i in sel])
        chi = per[0][3]
        Dr = (Z[:, None] * (F @ chi.conj().T) / 64).mean(0) / Z.mean()
        Nr = np.einsum("m,me,met,se->st", Z, F, G, chi.conj()) / 64 / len(Z) / Z.mean()
        return Dr, Nr
    allc = list(range(len(per)))
    Dr, Nr = est(allc)
    ref = int(np.argmax(Dr.real))
    nc = len(per)
    jk = [est([j for j in allc if j != i]) for i in allc] if nc > 1 else []
    D_err = np.sqrt((nc - 1) / nc * np.sum([(j_[0].real - np.mean([k_[0].real for k_ in jk], axis=0)) ** 2 for j_ in jk], axis=0)) if nc > 1 else np.full(len(Dr), np.nan)
    Di_err = np.sqrt((nc - 1) / nc * np.sum([(j_[0].imag - np.mean([k_[0].imag for k_ in jk], axis=0)) ** 2 for j_ in jk], axis=0)) if nc > 1 else np.full(len(Dr), np.nan)
    C = Nr[ref] / Dr[ref]
    C_jk = np.array([j_[1][ref] / j_[0][ref] for j_ in jk]) if nc > 1 else None
    C_err = np.sqrt((nc - 1) / nc * np.sum((C_jk.real - C_jk.real.mean(0)) ** 2, axis=0)) if nc > 1 else np.full(len(C), np.nan)
    Dimag = float(np.abs(Dr.imag).max() / abs(np.concatenate([q[0] for q in per]).mean()))
    chain_D = np.array([ (q[0] * (q[1] @ q[3][ref].conj()) / 64).mean().real / q[0].mean().real for q in per])
    chain_err = np.array([q[4]["D_ref_gamma_err"] for q in per])
    chi2 = float(np.sum((chain_D - chain_D.mean()) ** 2 / np.maximum(chain_err, 1e-9) ** 2) / max(nc - 1, 1)) if nc > 1 else np.nan
    within = float(np.sqrt(np.sum(chain_err ** 2)) / nc)
    N_eff = float(sum(q[4]["N_eff"] for q in per))
    R_phase = float(abs(np.concatenate([q[0] for q in per]).mean()))
    chain_C = np.array([runs[nm][1]["C_ref"].real for nm in names]); chain_Ce = np.array([runs[nm][1]["C_ref_err"] for nm in names])
    jt = [j for j, t in enumerate(per[0][5]["taus"]) if t in (0.5, 1.0, 1.5)]
    chi2C = float(np.mean([np.sum((chain_C[:, j] - chain_C[:, j].mean()) ** 2 / np.maximum(chain_Ce[:, j], 1e-9) ** 2) / max(nc - 1, 1) for j in jt])) if nc > 1 else np.nan
    imag_bad = int(np.sum(np.abs(Dr.imag) > 3 * Di_err)) if nc > 1 else -1
    R_proj_pool = float(abs(np.concatenate([q[0] * (q[1] @ q[3][ref].conj()) for q in per]).mean()) / np.mean(np.abs(np.concatenate([q[1] @ q[3][ref].conj() for q in per]))))
    return {"D": Dr, "D_err": D_err, "Di_err": Di_err, "imag_bad": imag_bad, "chi2C": chi2C, "R_proj": R_proj_pool, "ref": ref, "C": C, "C_err": C_err, "imag": Dimag, "chi2": chi2, "N_eff": N_eff, "R_phase": R_phase,
            "chain_D": chain_D, "jk_vs_within": float(D_err[ref] / max(within, 1e-12)) if nc > 1 else np.nan, "nc": nc, "sectors": per[0][5]["sectors"], "taus": per[0][5]["taus"]}


def pooled_verdict(P):
    fails = []
    if P["R_phase"] <= 0.5: fails.append("B1")
    if P["nc"] < 2: return "SPG_PROJECTION_UNRESOLVED", ["B6 (one chain)"]
    sig = P["D_err"][P["ref"]]; Dref = P["D"][P["ref"]].real
    if sig >= 1 / 64: fails.append("B2-noise")
    elif not (Dref >= 2 / 64 and Dref - 1 / 64 > 3 * sig): return "SPG_INACTIVE", ["B2"]
    if P["N_eff"] < 8: fails.append("B4")
    if np.any(P["D"].real < -3 * P["D_err"]) or np.any(P["D"].real > 1 + 3 * P["D_err"]) or P["imag_bad"] > 0: fails.append("B5")
    if not (P["chi2"] < 2 and P["chi2C"] < 2): fails.append("B6")
    if not (0.5 <= P["jk_vs_within"] <= 2): fails.append("B7")
    if fails:
        fails.append("cause: " + ("MIXING" if P["R_proj"] >= 0.5 else "PROJECTION_PHASE"))
    return ("SPG_RESOLVED" if not fails else "SPG_PROJECTION_UNRESOLVED"), fails


pool_groups = {}
for name in runs:
    pool_groups.setdefault(name[:7] if len(name) > 6 and name[6] in "LM" else name[:6], []).append(name)   # long-burn-in (..L<k>) and continued (..M<k>) chains pooled separately
plines, pooled_res = [], {}
for gk, mem in sorted(pool_groups.items()):
    P = pooled(sorted(mem)); v, fl = pooled_verdict(P); pooled_res[gk] = (P, v, fl)
    nneg = int(np.sum(P["D"].real < -3 * P["D_err"])) if P["nc"] > 1 else -1
    plines.append(f"{gk} ({P['nc']} chains): {v} {fl if fl else ''} -- D_ref {P['D'][P['ref']].real:.3f} +- {P['D_err'][P['ref']]:.3f} ({tuple(P['sectors'][P['ref']])}), chain means {np.round(P['chain_D'], 3).tolist()}, chi^2/dof {P['chi2']:.2f}, "
                  f"N_eff {P['N_eff']:.1f}, sectors below -3 sigma {nneg}, Im D beyond 3 sigma {P['imag_bad']}, max |Im D| {P['imag']:.3f}, chi^2/dof of C_ref chains {P['chi2C']:.2f}, pooled R^proj {P['R_proj']:.3f}, chain-jackknife / within-chain error {P['jk_vs_within']:.2f}; C_ref(tau) " + ", ".join(f"{t:g}: {c.real:.3f}+-{e:.3f}" for t, c, e in zip(P["taus"], P["C"], P["C_err"])))
step("SPP9 (statement) pooled SPG budget per (g, beta, lambda) over independent chains (chain-level jackknife, B6 from the chain means)",
     True, " | ".join(plines) if plines else "no runs")

fA, fB = HERE / "spg_points_g018b6.json", HERE / "spg_points_g018b8.json"
if fA.exists() and fB.exists():
    dA, dB = json.loads(fA.read_text()), json.loads(fB.read_text())
    eA, eB = evaluate(dA), evaluate(dB)
    ref = int(np.argmax(eB["D"].real)); txt = []
    for s_ in np.argsort(-eB["D"].real)[:4]:
        if s_ == ref:
            continue
        rA, rB = eA["D"][s_].real / eA["D"][ref].real, eB["D"][s_].real / eB["D"][ref].real
        sA = rA * np.hypot(eA["D_err"][s_] / eA["D"][s_].real, eA["D_err"][ref] / eA["D"][ref].real)
        sB = rB * np.hypot(eB["D_err"][s_] / eB["D"][s_].real, eB["D_err"][ref] / eB["D"][ref].real)
        dE = -(np.log(rB) - np.log(rA)) / (dB["beta"] - dA["beta"]); sdE = np.hypot(sA / rA, sB / rB) / (dB["beta"] - dA["beta"])
        txt.append(f"{tuple(dB['sectors'][s_])} vs {tuple(dB['sectors'][ref])}: D ratio {rA:.3f} (beta {dA['beta']:g}) -> {rB:.3f} (beta {dB['beta']:g}), slope Delta E = {dE:+.3f} +- {sdE:.3f}")
    fixed = ", ".join(f"tau {t:g}: {a.real:.3f} -> {b.real:.3f}" for t, a, b in zip(dA["taus"], eA["C_ref"], eB["C_ref"]))
    noise = f"min D_r {eA['D'].real.min():.3f} / {eB['D'].real.min():.3f} (negative = noise), jackknife errors from 16 measurements (tau_int ~ 5 sweeps: underestimated ~x1.6)"
    step("SPP6 (statement) promotion test at g = 0.018 (lambda = 0.39, 16 sites): SECTOR_LOW_GAP from the multiplicity-free beta-slope of the sector ratios between beta = 6 and 8 (free mirror addition energy 0.53), SPG_LOW_GAP from the reference-sector curve at fixed tau -- interpretable only if SPP8 gives SPG_RESOLVED for both temperatures",
         True, " | ".join(txt) + f" | reference-sector C_r(tau), beta {dA['beta']:g} -> {dB['beta']:g}: {fixed} | {noise} | budget: beta {dA['beta']:g} {verdicts.get('g018b6', ('?',))[0]}, beta {dB['beta']:g} {verdicts.get('g018b8', ('?',))[0]}")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG)}
(HERE / "spg_production_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
