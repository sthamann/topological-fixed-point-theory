#!/usr/bin/env python3
"""Part RW: CHIRAL4D.REWEIGHT.01 -- phase-quenched sampling of Theta_1 and the severity of its phase problem.

Theta_1 frozen: minimal-doubling C3 kinetics h_W (theta1_md.py, walk = "md", default); the legacy BDT Hamiltonian walk
H_w = (i/2)(W - W^dag) is available as walk = "bdt" (all runs before 2026-10-02 09:55 used it); frozen filter F, ten Spin(10)
channels, trivial links; the full squares O_A^2 are HS-decoupled (the normal-ordering one-body term is NOT added again).
Per slice: S_l = e^{K_w/2} e^{K_phi,l} e^{K_w/2}; W[phi]^2 = det(1 + prod S_l) exactly; phase-quenched measure
P_pq ~ exp(-sum phi^2 / (4 g dtau)) |W[phi]|; the phase theta of W is tracked continuously along the Markov chain
(start phi = 0 where W > 0), never by a principal square root; determinants via a pivoted-QR UDV chain.
Measured per point: R = |<e^{i theta}>_pq| (binned errors), Var theta, Delta f = -log R / (beta V),
N_req = (1 - R^2) / (eps^2 R^2) for eps = 0.1, integrated autocorrelation of cos theta, smallest singular value of the
stabilised inner matrix (zero detector), acceptance; with "meas": reweighted G_k(tau) at Gamma, +-P, H (mirror-gap estimator
G_P(beta/2) = 1 / (2 cosh(beta Delta / 2)), free reference from the same code path) and the uniform susceptibility chi_O of
the ten filtered pairing composites from the HS fields (symmetry-breaking diagnostic); 3 x 3 correlation matrix of independent
composite seeds at P (degree 1 psi_P, degree 3 (2 M^dag psi^dag)_P, degree 5 (4 M^dag M psi)_P) with a whitened GEVP.
Usage: python3 -B reweight.py <cell> <beta> <dtau> <g> <sweeps> <seed> [tag] [walk md|bdt] [meas] [dx,dy,dz flavour twist|none] [t1x,t1y,t1z;t2x,... target twists]   (cell: debug4 | cell8 | cubic16 | oh32 | prim64)
Results are appended to reweight_points.json; reweight_check.py evaluates them.
"""
from __future__ import annotations

import json
import sys
import time
from itertools import product
from pathlib import Path

import numpy as np
from scipy.linalg import expm, qr

HERE = Path(__file__).resolve().parent
B3 = np.array([[1, 1, -1], [1, -1, 1], [-1, 1, 1]])
CELLS = {"debug4": np.array([B3[0] - B3[1], B3[1] - B3[2], B3[0] + B3[1] + 2 * B3[2]]),
         "cell8": np.array([[4, 0, 0], [2, 2, 0], [2, 1, 1]]) @ B3,          # best index-8 P-compatible cell (supercell_search.py)
         "cubic16": np.array([[4, 0, 0], [0, 4, 0], [0, 0, 4]]),
         "oh32": np.array([[8, 0, 0], [6, 2, 0], [6, 0, 2]]) @ B3,           # full-O_h index-32 P-compatible cell (supercell_search.py)
         "prim64": 4 * B3}                                                  # primitive 4 x 4 x 4 cell (contains P); 16 / 32 / 64 is NOT homothetic
signs = [np.array(s) for s in product((-1, 1), repeat=3)]
I2 = np.eye(2)
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1.0 + 0j, -1])


def A_eps(e):
    return (I2 + e[0] * sx) / 2 @ (I2 + e[1] * sy) / 2 @ (I2 + e[2] * sz) / 2


h_eps = {tuple(e): 0.5j * (A_eps(e) - A_eps(-e).conj().T) for e in signs}                 # legacy BDT Hamiltonian walk (4 nodes)
# frozen Theta_1 kinetics (theta1_md.py): h_W(k) = t sum_i d_i(p) tau_i, d_i = sin p_i + cos p_{i+1} - 1, p = B3 k, tau_i = sum_a (B3^-1)_{ai} sigma_a
_Binv = np.linalg.inv(B3.astype(float))
_tau = [sum(_Binv[a, i] * [sx, sy, sz][a] for a in range(3)) for i in range(3)]
h_md = {}
for i in range(3):
    h_md[tuple(B3[i])] = 0.5j * _tau[i] + 0.5 * _tau[(i - 1) % 3]
    h_md[tuple(-B3[i])] = -0.5j * _tau[i] + 0.5 * _tau[(i - 1) % 3]
h_md_onsite = -sum(_tau)


def build_cell(Rv):
    Rinv = np.linalg.inv(Rv.astype(float))
    def canon(x):
        fr = np.mod(np.round(x @ Rinv, 9), 1.0)
        return tuple(np.round(fr, 6))
    sites, keys = [], {}
    for m in product(range(-6, 7), repeat=3):
        x = np.array(m) @ B3
        k = canon(x)
        if k not in keys:
            keys[k] = len(sites); sites.append(x)
    N = len(sites)
    F = np.eye(N)
    nb = {}
    for i, x in enumerate(sites):
        for e in signs:
            j = keys[canon(x + e)]
            F[i, j] -= 1 / 8
            nb[(i, tuple(e))] = j
    return sites, F, nb


def windings(Rv, sites, nb):
    """integer superlattice winding n of each bond: x_i + e = x_j + n Rv (rows of Rv = superlattice vectors)."""
    Rinv = np.linalg.inv(Rv.astype(float))
    return {(i, e): np.rint((sites[i] + np.array(e) - sites[j]) @ Rinv).astype(int) for (i, e), j in nb.items()}


def spin10():
    Zq = np.diag([1.0, -1.0]).astype(complex); lowc = np.array([[0, 1], [0, 0]], dtype=complex)
    c5 = []
    for j in range(5):
        mats = [Zq] * j + [lowc] + [np.eye(2, dtype=complex)] * (4 - j)
        out = mats[0]
        for m_ in mats[1:]:
            out = np.kron(out, m_)
        c5.append(out)
    Gm = []
    for j in range(5):
        Gm += [c5[j] + c5[j].conj().T, -1j * (c5[j] - c5[j].conj().T)]
    occ = np.array([[(b >> (4 - j)) & 1 for j in range(5)] for b in range(32)])
    Pe = np.eye(32)[:, [b for b in range(32) if occ[b].sum() % 2 == 0]]
    Cc = Gm[0] @ Gm[2] @ Gm[4] @ Gm[6] @ Gm[8]
    return [np.kron(np.array([[0, 1], [-1, 0]]), Pe.T @ (Cc @ Gm[A]) @ Pe) for A in range(10)]


MA = spin10()
_even = [b_ for b_ in range(32) if bin(b_).count("1") % 2 == 0]
QF = np.array([0.5 - ((b_ >> 4) & 1) for b_ in _even])                 # Sigma_12 = -(i/2) Gamma_1 Gamma_2 = 1/2 - n_0 on the 16
FAST = True                                                           # split-chain determinant updates (identical chains to the full QR chain, checked on debug4 and cubic16)
wrap = lambda a: (a + np.pi) % (2 * np.pi) - np.pi


class Model:
    def __init__(self, cell, beta, dtau, g, walk="md", twist=None):
        self.sites, self.F, self.nb = build_cell(CELLS[cell])
        self.N = len(self.sites)
        n = 32 * self.N
        self.twist = None if twist is None or not np.any(twist) else np.asarray(twist, float)
        # flavour-twisted boundary conditions psi(x + R) = exp(2 i (R.delta) Q) psi(x), Q = Sigma_12 in Spin(10): a symmetry
        # of H_int (the ten channels rotate orthogonally), so the twisted cell is the same theory with momenta k +- delta
        if self.twist is not None:
            wd = windings(CELLS[cell], self.sites, self.nb)
            Rv = CELLS[cell].astype(float)
            Uf = lambda nvec: np.exp(2j * (np.asarray(nvec) @ Rv @ self.twist) * QF)       # flavour phases of a winding
        Hw = np.zeros((n, n), complex)
        hop = h_md if walk == "md" else h_eps
        for i in range(self.N):
            for e in signs:
                if tuple(e) not in hop:
                    continue
                j = self.nb[(i, tuple(e))]
                fl = np.eye(16) if self.twist is None else np.diag(Uf(wd[(i, tuple(e))]).conj())
                Hw[32 * j:32 * j + 32, 32 * i:32 * i + 32] += np.kron(hop[tuple(e)], fl)
            if walk == "md":
                Hw[32 * i:32 * i + 32, 32 * i:32 * i + 32] += np.kron(h_md_onsite, np.eye(16))
        if self.twist is not None:                                     # rho_i = F_theta[i] psi (32 x n)
            self.Ftw = []
            for i in range(self.N):
                Fi = np.zeros((32, n), complex)
                Fi[:, 32 * i:32 * i + 32] += np.eye(32)
                for e in signs:
                    j = self.nb[(i, tuple(e))]
                    Fi[:, 32 * j:32 * j + 32] -= np.kron(np.eye(2), np.diag(Uf(wd[(i, tuple(e))]))) / 8
                self.Ftw.append(Fi)
        self.Hw, self.walk, self.cell = Hw, walk, cell
        assert np.allclose(Hw, Hw.conj().T)
        Kw = -dtau * np.block([[Hw, np.zeros_like(Hw)], [np.zeros_like(Hw), -Hw.T]])
        self.Wh = expm(Kw / 2)
        self.FF = [np.outer(self.F[i], self.F[i]) for i in range(self.N)]
        self.beta, self.dtau, self.g = beta, dtau, g
        self.Lt = int(round(beta / dtau))
        self.n = n

    def build_M(self, phi):
        M = np.zeros((self.n, self.n), complex)
        if self.twist is None:
            for i in range(self.N):
                M += np.kron(self.FF[i], np.tensordot(phi[i], MA, axes=(0, 0)))
        else:
            for i in range(self.N):
                M += self.Ftw[i].T @ np.tensordot(phi[i], MA, axes=(0, 0)) @ self.Ftw[i]
        return M

    def slice_op(self, phi):
        M = self.build_M(phi)
        # e^{K}, K = [[0, 2M^dag], [2M, 0]], exactly via the SVD M = U s V^dag (half the dimension of an eigh of K)
        try:
            U, sv, Vh = np.linalg.svd(M)
        except np.linalg.LinAlgError:                                  # LAPACK gesdd non-convergence: fall back to gesvd
            from scipy.linalg import svd as sp_svd
            U, sv, Vh = sp_svd(M, lapack_driver="gesvd")
        V = Vh.conj().T
        ch, sh = np.cosh(2 * sv), np.sinh(2 * sv)
        E = np.block([[(V * ch) @ Vh, (V * sh) @ U.conj().T], [(U * sh) @ Vh, (U * ch) @ U.conj().T]])
        return self.Wh @ E @ self.Wh

    def momenta(self):
        """Gamma, P, -P, H if allowed in the cell (e^{i k.R} = 1 on the superlattice)."""
        Rv = CELLS[self.cell].astype(float)
        cand = {"Gamma": np.zeros(3), "P": np.pi / 2 * np.ones(3), "-P": -np.pi / 2 * np.ones(3), "H": np.array([np.pi, 0, 0])}
        return {k: v for k, v in cand.items() if np.allclose(np.cos(Rv @ v), 1, atol=1e-9)}

    @staticmethod
    def _udv_left(M, udv):
        """M (U D V) -> U' D' V' (pivoted QR of (M U) D, column-graded: stable)."""
        U, D, V = udv
        Q, Rr, piv = qr((M @ U) * D[None, :], pivoting=True)
        d = np.abs(np.diag(Rr))
        return Q, d, (Rr / d[:, None])[:, np.argsort(piv)] @ V

    @staticmethod
    def _inv_sum(a, b):
        """(U1 D1 V1 + U2 D2 V2)^-1, U unitary, with the big/small split of both D (Loh-Gubernatis)."""
        U1, D1, V1 = a; U2, D2, V2 = b
        Db1, Ds1 = np.maximum(D1, 1.0), np.minimum(D1, 1.0)
        Db2, Ds2 = np.maximum(D2, 1.0), np.minimum(D2, 1.0)
        V2i = np.linalg.inv(V2)
        mid = Ds1[:, None] * (V1 @ V2i) / Db2[None, :] + (U1.conj().T @ U2) * Ds2[None, :] / Db1[:, None]
        return (V2i / Db2[None, :]) @ np.linalg.solve(mid, U1.conj().T / Db1[:, None])

    def greens_k(self, phi, ks):
        """reweightable estimator G_k(tau_j) = (1 / 32) tr <psi_k(tau_j) psi_k^dag(0)>, tau_j = j dtau, j = 0..Lt, normal block,
        G(tau_C) = (S_C^-1 + S_A)^-1 with S_A = S_0..S_{l-1}, S_C = S_l..S_{L-1} (validated by brute force).  Stabilised:
        S_A^T = S_{l-1}^T..S_0^T and (S_C^-1)^T = S_l^-T..S_{L-1}^-T are left-multiplication chains (column-graded QR), so
        G^T = (U_C D_C V_C + U_A D_A V_A)^-1 with unitary U; slice inverses exact (E(phi)^-1 = E(-phi))."""
        if not hasattr(self, "Wh_inv"):
            self.Wh_inv = np.linalg.inv(self.Wh)
        L = self.Lt; nn = 2 * self.n; I = np.eye(nn, dtype=complex)
        unit = (I, np.ones(nn), I)
        S = [self.slice_op(phi[l]) for l in range(L)]
        Sinv = [self.Wh_inv @ (self.Wh_inv @ self.slice_op(-phi[l]) @ self.Wh_inv) @ self.Wh_inv for l in range(L)]
        A = [unit]                                                   # A[l] = UDV of S_A^T = (S_0 .. S_{l-1})^T
        for l in range(L):
            A.append(self._udv_left(S[l].T, A[-1]))
        C = [unit] * (L + 1)                                          # C[l] = UDV of (S_C^-1)^T = ((S_l .. S_{L-1})^-1)^T
        for l in range(L - 1, -1, -1):
            C[l] = self._udv_left(Sinv[l].T, C[l + 1])
        X = np.array(self.sites, float)
        if self.twist is not None:                                     # Q = +1/2 flavours sit at k + delta, Q = -1/2 at k - delta
            ks = {**{k + "+d": v + self.twist for k, v in ks.items()}, **{k + "-d": v - self.twist for k, v in ks.items()}}
        ph = {k: np.exp(-1j * X @ v) for k, v in ks.items()}
        out = {k: np.zeros(L + 1, complex) for k in ks}
        spin = {k: np.zeros((L + 1, 2, 2), complex) for k in ks}          # spin-resolved C_st(tau, k), flavour-averaged per sector
        # composite seeds at P with the quantum numbers of psi (Nambu rows, normal component of the momentum projector):
        # degree 1: psi_P; degree 3: (2 M(phi_l)^dag psi^dag)_P = HS form of [psi, phi.O]; degree 5: (4 M^dag M psi)_P.
        # The operator at boundary l uses the next slice l (cyclically 0 for the end insertion); all vanish at Gamma.
        nN = self.n
        PiP = np.kron(ph["P"][None, :] / np.sqrt(self.N), np.eye(32)) if ("P" in ks and self.twist is None) else None
        def seeds(l):
            M = self.build_M(phi[l % L])
            Md = M.conj().T
            z = np.zeros((32, nN), complex)
            s1 = np.hstack([PiP, z])
            s3 = np.hstack([z, 2 * PiP @ Md]) / self.dtau
            s5 = np.hstack([4 * PiP @ Md @ M, z]) / self.dtau ** 2
            return [s1, s3, s5]
        Cmat = np.zeros((L + 1, 3, 3), complex)
        s0 = seeds(0) if PiP is not None else None
        for l in range(L + 1):
            G = self._inv_sum(C[l], A[l]).T
            Gr = G[:self.n, :self.n].reshape(self.N, 2, 16, self.N, 2, 16)
            j = L - l                                                  # tau = (L - l) dtau
            if self.twist is None:
                Gst = np.einsum("xsfytf->stxy", Gr)
                for k in ks:
                    spin[k][j] = np.einsum("x,stxy,y->st", ph[k], Gst, ph[k].conj()) / (16 * self.N)
                    out[k][j] = np.trace(spin[k][j]) / 2
            else:
                Gq = {q: np.einsum("xsfytf->stxy", Gr[:, :, QF == q][:, :, :, :, :, QF == q]) for q in (0.5, -0.5)}
                for k in ks:
                    q = 0.5 if k.endswith("+d") else -0.5
                    spin[k][j] = np.einsum("x,stxy,y->st", ph[k], Gq[q], ph[k].conj()) / (8 * self.N)
                    out[k][j] = np.trace(spin[k][j]) / 2
            if l == L and self.twist is None:
                # equal-time density structure factor S_n(q) = (1/N) sum_xy e^{iq(x-y)} <dn_x dn_y> by Wick with the Nambu G:
                # connected <n_a n_b> = <c_a^dag c_b><c_a c_b^dag> - <c_a^dag c_b^dag><c_a c_b>
                n_ = self.n
                G11, G12, G21 = G[:n_, :n_], G[:n_, n_:], G[n_:, :n_]
                Cc = (np.eye(n_) - G11.T) * G11 - G21 * G12
                Cxy = Cc.reshape(self.N, 32, self.N, 32).sum(axis=(1, 3))
                out["S_density"] = {k: complex(ph[k].conj() @ Cxy @ ph[k] / self.N) for k in ks}
            if PiP is not None:
                sl = seeds(l)
                for a_ in range(3):
                    for b_ in range(3):
                        Cmat[j, a_, b_] = np.trace(sl[a_] @ G @ s0[b_].conj().T) / 32
        if PiP is not None:
            out["seeds_P"] = Cmat
        out["spin"] = spin
        return out

    def sweep_cache(self, S):
        """Y[l] = UDV of S_{l+1} .. S_{L-1} (left chain from the end, Y[L-1] = 1)."""
        L = len(S); nn = S[0].shape[0]; I = np.eye(nn, dtype=complex)
        Y = [None] * L
        Y[L - 1] = (I, np.ones(nn), I)
        for l in range(L - 2, -1, -1):
            Y[l] = self._udv_left(S[l + 1], Y[l + 1])
        return Y

    @staticmethod
    def logdet_split(Sl, X, Y):
        """log det(1 + X_l S_l Y_l) = log det(1 + S_0 .. S_{L-1}) with X_l = S_0 .. S_{l-1} = V_X^T D_X U_X^T (UDV of X_l^T,
        a left chain) and Y_l = S_{l+1} .. S_{L-1} = U_Y D_Y V_Y:  1 + X S Y = V_X^T (V_X^-T V_Y^-1 + D_X M D_Y) V_Y with
        M = U_X^T S_l U_Y, big/small split of both D.  Cost: two products and three LU instead of a full QR chain."""
        Ux, Dx, Vx = X; Uy, Dy, Vy = Y
        Mid = Ux.T @ Sl @ Uy
        W0 = np.linalg.inv(Vx.T) @ np.linalg.inv(Vy)
        Xb, Xs = np.maximum(Dx, 1.0), np.minimum(Dx, 1.0)
        Yb, Ys = np.maximum(Dy, 1.0), np.minimum(Dy, 1.0)
        br = W0 / Xb[:, None] / Yb[None, :] + Xs[:, None] * Mid * Ys[None, :]
        s1, l1 = np.linalg.slogdet(Vx); s2, l2 = np.linalg.slogdet(Vy); s3, l3 = np.linalg.slogdet(br)
        return l1 + l2 + l3 + np.sum(np.log(Xb)) + np.sum(np.log(Yb)), np.angle(s1) + np.angle(s2) + np.angle(s3)

    def logdet(self, S, want_smin=False):
        nn = S[0].shape[0]
        U, D, V = np.eye(nn, dtype=complex), np.ones(nn), np.eye(nn, dtype=complex)
        for M in reversed(S):
            X = (M @ U) * D[None, :]
            Q, Rr, piv = qr(X, pivoting=True)
            d = np.abs(np.diag(Rr))
            V = (Rr / d[:, None])[:, np.argsort(piv)] @ V
            U, D = Q, d
        Db, Ds = np.maximum(D, 1.0), np.minimum(D, 1.0)
        inner = (U.conj().T @ np.linalg.inv(V)) / Db[:, None] + np.diag(Ds)
        s1, l1 = np.linalg.slogdet(U); s2, l2 = np.linalg.slogdet(inner); s3, l3 = np.linalg.slogdet(V)
        logabs = l1 + l2 + l3 + np.sum(np.log(Db))
        phase = np.angle(s1) + np.angle(s2) + np.angle(s3)
        smin = np.linalg.svd(inner, compute_uv=False)[-1] if want_smin else np.nan
        return logabs, phase, smin


def run(cell, beta, dtau, g, sweeps, seed, walk="md", measure=False, twist=None, targets=None, hook=None, hook_every=1, outer=None,
        burn=None, monitor=None, monitor_every=1, save_phi=None, phi0=None, global_moves=None, bias=None):
    """hook(m, phi) -> any: optional extra measurement on every hook_every-th measured configuration (raw records with the
    tracked phase are returned in res["hook"]); does not touch the integrator.
    outer(m, phi) -> (h, payload): QMC.SECTOR_SAMPLING envelope (surrogate-transition Metropolis): every sweep runs forward or
    backward with probability 1/2 (the mixture is reversible w.r.t. the phase-quenched measure), then the new configuration
    is accepted with min(1, h' / h) -> stationary measure e^{-S} |W| h; payloads and h of the current state are recorded in
    res["outer"] after thermalisation.  outer = None leaves the sampler and its random stream unchanged.
    burn: thermalisation sweeps (default sweeps // 5; measurements and R start after it, the step is adapted during it);
    monitor(m, phi, state) -> dict: cheap per-sweep record (every monitor_every sweeps, also during burn-in) in res["monitor"];
    save_phi: path for the final HS configuration (.npy); phi0: start configuration (continuation) instead of a prior draw.
    global_moves = {"scale_eps": e, "pcn_rho": r, "every": k} (QMC.MIXING.MODE.01 targeted update of the identified slow
    direction, the overall HS amplitude): after every k-th sweep (i) phi -> s phi, s = exp(e xi), accepted with
    min(1, s^D e^{-(s^2 - 1)|phi|^2 / (4 g dtau)} |W(s phi)| / |W(phi)|) (D = number of field components: Jacobian of the
    multiplicative move), (ii) pCN on the time-constant modes c(x, A) = sum_l phi_l(x, A) / sqrt(L): c' = sqrt(1 - r^2) c +
    r sigma xi (prior exact), accepted with min(1, |W'| / |W|); the phase is continued along the straight path; own random
    stream (global_moves = None leaves the sampler unchanged).
    bias = {"center": c, "kappa": k[, "kind": "norm2" | "uniform"]}: umbrella window in |phi|^2 or in the K = 0 pairing
    order parameter X = |M|^2 / N (extra action k/2 (coord - c)^2 in every acceptance; removed again by WHAM / umbrella
    integration); bias = None leaves the sampler unchanged."""

    rng = np.random.default_rng(seed)
    rec_H = []
    rng_o = np.random.default_rng(seed + 7919); rec_O = []; o_acc = 0
    burn = sweeps // 5 if burn is None else int(burn)
    rec_M = []
    rng_g = np.random.default_rng(seed + 104729); g_acc = {"scale": [0, 0], "pcn": [0, 0]}
    def bias_val(ph_):
        """umbrella bias kappa/2 (coord - centre)^2; coord = |phi|^2 (kind 'norm2', default) or X = |M|^2 / N with
        M_A = sum_{l,x} phi_A in prior units, X = sum_A M_A^2 / (L N 2 g dtau) / N (kind 'uniform', the K = 0 order parameter)."""
        if bias.get("kind", "norm2") == "uniform":
            Mv = ph_.sum(axis=(0, 1)); coord = float(np.sum(Mv ** 2)) / (ph_.shape[0] * ph_.shape[1] * 2 * g * dtau) / ph_.shape[1]
        else:
            coord = float(np.sum(ph_ ** 2))
        return 0.5 * bias["kappa"] * (coord - bias["center"]) ** 2
    m = Model(cell, beta, dtau, g, walk, twist)
    # TWIST.MULTIREWEIGHT.01: same HS configurations, target twists reweighted by W_t / |W_anchor| (phase continued from the
    # tracked anchor phase by the nearest branch; the jump is monitored)
    tmods = [Model(cell, beta, dtau, g, walk, t) for t in (targets or [])]
    rec_T = [[] for _ in tmods]
    ks = m.momenta() if measure else {}
    rec_G, rec_Phi, sentinel, phimax = [], [], 0.0, 0.0
    Lt, N = m.Lt, m.N
    sig_prior = np.sqrt(2 * g * dtau)
    phi = sig_prior * rng.normal(size=(Lt, N, 10))                    # start from the prior; initial phase tracked from t = 0
    if phi0 is not None:
        phi = np.array(phi0, float).reshape(Lt, N, 10).copy()          # continuation (the prior draw above keeps the random stream)
    def ph_at(t):
        return m.logdet([m.slice_op(t * phi[l]) for l in range(Lt)])[1]
    thetaD, t_, ph_prev, dt = 0.0, 0.0, ph_at(0.0), 0.125
    while t_ < 1 - 1e-12:
        dt = min(dt, 1 - t_)
        ph_new = ph_at(t_ + dt)
        inc = wrap(ph_new - ph_prev)
        if abs(inc) > np.pi / 4 and dt > 1e-4:
            dt /= 2; continue
        thetaD += inc; ph_prev = ph_new; t_ += dt; dt *= 1.5
    S = [m.slice_op(phi[l]) for l in range(Lt)]
    la, ph_D, smin = m.logdet(S)
    step_sig = 0.5 * sig_prior / np.sqrt(N)
    acc, tries = 0, 0
    acc_prev, tries_prev = 0, 0
    rec_theta, rec_smin = [], []
    t0 = time.time()
    nn_ = 2 * m.n; unit_ = (np.eye(nn_, dtype=complex), np.ones(nn_), np.eye(nn_, dtype=complex)); drift = 0.0
    def upd(l, Xl, Yl):
        nonlocal S, la, ph_D, smin, thetaD, acc, tries
        new = phi[l] + step_sig * rng.normal(size=(N, 10))
        dS = (np.sum(new ** 2) - np.sum(phi[l] ** 2)) / (4 * g * dtau)
        if bias is not None:
            ph_n = phi.copy(); ph_n[l] = new
            dS += bias_val(ph_n) - bias_val(phi)
        Sl_new = m.slice_op(new)
        if FAST:
            la2, ph2 = m.logdet_split(Sl_new, Xl, Yl); sm2 = np.nan
        else:
            S_new = S.copy(); S_new[l] = Sl_new
            la2, ph2, sm2 = m.logdet(S_new)
        inc = wrap(ph2 - ph_D)
        if abs(inc) > np.pi / 2:                                      # refine the path to keep the phase continuous
            inc_tot, cur, ph_cur = 0.0, phi[l].copy(), ph_D
            for k in range(1, 9):
                mid = phi[l] + (new - phi[l]) * k / 8
                if FAST:
                    phm = m.logdet_split(m.slice_op(mid), Xl, Yl)[1]
                else:
                    Sm = S.copy(); Sm[l] = m.slice_op(mid)
                    _, phm, _ = m.logdet(Sm)
                inc_tot += wrap(phm - ph_cur); ph_cur = phm
            inc = inc_tot
        tries += 1
        if np.log(rng.random()) < -dS + 0.5 * (la2 - la):
            phi[l] = new; S = S.copy(); S[l] = Sl_new; la, ph_D, smin = la2, ph2, sm2
            thetaD += inc; acc += 1
    if outer is not None:
        h_cur, pay_cur = outer(m, phi)
    for sw in range(sweeps):
        backward = False
        if outer is not None:
            saved = (phi.copy(), list(S), la, ph_D, thetaD, smin)
            backward = rng_o.random() < 0.5
        if FAST:
            Ych, Xch = m.sweep_cache(S), unit_                        # Y_l = S_{l+1}..S_{L-1}; X_l = S_0..S_{l-1} grows below
            la_s, ph_s = m.logdet_split(S[0], unit_, Ych[0])          # resynchronise with a fresh chain once per sweep
            drift = max(drift, abs(la_s - la))
            thetaD += wrap(ph_s - ph_D); la, ph_D = la_s, ph_s
        if not backward:
            for l in range(Lt):
                upd(l, Xch if FAST else None, Ych[l] if FAST else None)
                if FAST:
                    Xch = m._udv_left(S[l].T, Xch)
        else:                                                         # backward sweep: prefixes fixed, suffix grows
            Xc = [unit_]
            for l in range(Lt - 1):
                Xc.append(m._udv_left(S[l].T, Xc[-1]))
            Yb = unit_
            for l in range(Lt - 1, -1, -1):
                upd(l, Xc[l], Yb)
                Yb = m._udv_left(S[l], Yb)
        if outer is not None:
            h_new, pay_new = outer(m, phi)
            if rng_o.random() < h_new / max(h_cur, 1e-300):
                h_cur, pay_cur = h_new, pay_new; o_acc += 1
            else:
                phi[:] = saved[0]; S = saved[1]; la, ph_D, thetaD, smin = saved[2], saved[3], saved[4], saved[5]
            if sw >= burn:
                rec_O.append((thetaD / 2, h_cur, pay_cur))
        if global_moves is not None and sw % int(global_moves.get("every", 1)) == 0:
            sig2 = 2 * g * dtau; Dn = phi.size
            def global_try(phi_new, extra_logw):
                nonlocal phi, S, la, ph_D, smin, thetaD
                # phase along the straight path phi(t) = phi + t (phi_new - phi), adaptive like the initial tracking
                dphi = phi_new - phi
                def ph_lt(t):
                    return m.logdet([m.slice_op(phi[l] + t * dphi[l]) for l in range(Lt)])
                t_, dt, inc_tot, ph_prev = 0.0, 0.5, 0.0, ph_D
                la_new = None
                while t_ < 1 - 1e-12:
                    dt = min(dt, 1 - t_)
                    la_t, ph_t, _ = ph_lt(t_ + dt)
                    inc = wrap(ph_t - ph_prev)
                    if abs(inc) > np.pi / 4 and dt > 1e-4:
                        dt /= 2; continue
                    inc_tot += inc; ph_prev = ph_t; t_ += dt; dt *= 1.5; la_new = la_t
                logacc = extra_logw + 0.5 * (la_new - la)
                if np.log(rng_g.random()) < logacc:
                    phi = phi_new; S = [m.slice_op(phi[l]) for l in range(Lt)]; la, ph_D = la_new, ph_prev; thetaD += inc_tot
                    return True
                return False
            e_ = float(global_moves.get("scale_eps", 0.0))
            if e_ > 0:
                sfac = np.exp(e_ * rng_g.normal())
                n2 = float(np.sum(phi ** 2))
                lw = Dn * np.log(sfac) - (sfac ** 2 - 1) * n2 / (4 * g * dtau)
                if bias is not None:
                    lw -= bias_val(sfac * phi) - bias_val(phi)
                g_acc["scale"][1] += 1
                if global_try(sfac * phi, lw):
                    g_acc["scale"][0] += 1
            eu_ = float(global_moves.get("uniform_scale_eps", 0.0))
            if eu_ > 0:                                                   # scale only the 10 uniform channel coordinates c_u -> s c_u
                cu = phi.sum(axis=(0, 1)) / np.sqrt(Lt * N)
                su = np.exp(eu_ * rng_g.normal())
                phi_s = phi + (((su - 1) * cu) / np.sqrt(Lt * N))[None, None, :]
                lwu = cu.size * np.log(su) - (su ** 2 - 1) * float(np.sum(cu ** 2)) / (2 * sig2)   # Jacobian s^10, Gaussian prior of c_u
                if bias is not None:
                    lwu -= bias_val(phi_s) - bias_val(phi)
                g_acc.setdefault("uniform_scale", [0, 0]); g_acc["uniform_scale"][1] += 1
                if global_try(phi_s, lwu):
                    g_acc["uniform_scale"][0] += 1
            ru_ = float(global_moves.get("pcn_uniform_rho", 0.0))
            if ru_ > 0:                                                   # pCN on the 10 uniform (K = 0, omega = 0) channel modes only
                cu = phi.sum(axis=(0, 1)) / np.sqrt(Lt * N)
                cu_new = np.sqrt(1 - ru_ ** 2) * cu + ru_ * np.sqrt(sig2) * rng_g.normal(size=cu.shape)
                phi_u = phi + ((cu_new - cu) / np.sqrt(Lt * N))[None, None, :]
                lbu = 0.0
                if bias is not None:
                    lbu = -(bias_val(phi_u) - bias_val(phi))
                g_acc.setdefault("pcn_uniform", [0, 0]); g_acc["pcn_uniform"][1] += 1
                if global_try(phi_u, lbu):
                    g_acc["pcn_uniform"][0] += 1
            r_ = float(global_moves.get("pcn_rho", 0.0))
            if r_ > 0:
                c = phi.sum(axis=0) / np.sqrt(Lt)                             # (N, 10), orthonormal time-constant coordinates
                c_new = np.sqrt(1 - r_ ** 2) * c + r_ * np.sqrt(sig2) * rng_g.normal(size=c.shape)
                g_acc["pcn"][1] += 1
                phi_p = phi + ((c_new - c) / np.sqrt(Lt))[None, :, :]
                lb = 0.0
                if bias is not None:
                    lb = -(bias_val(phi_p) - bias_val(phi))
                if global_try(phi_p, lb):
                    g_acc["pcn"][0] += 1
        rec_theta.append(thetaD / 2)
        if monitor is not None and sw % monitor_every == 0:
            rec_M.append(monitor(m, phi, {"sw": sw, "theta": thetaD / 2, "logabsW": 0.5 * la, "burn": burn}))
        if measure and sw >= burn:
            Gm_ = m.greens_k(phi, ks)
            rec_G.append((thetaD / 2, Gm_))
            if hook is not None and (sw - burn) % hook_every == 0:
                rec_H.append((thetaD / 2, Gm_, hook(m, phi)))
            for it, mt in enumerate(tmods):
                St = [mt.slice_op(phi[l]) for l in range(Lt)]
                la_t, ph_t, _ = mt.logdet(St)
                la_a, ph_a, _ = m.logdet(S)
                jump = 0.5 * wrap(ph_t - ph_a)
                rec_T[it].append((0.5 * (la_t - la_a), thetaD / 2 + jump, abs(jump), mt.greens_k(phi, mt.momenta())))
            if "Gamma" in Gm_:
                sentinel = max(sentinel, float(np.abs(Gm_["Gamma"] - 0.5).max()))         # Gamma decouples exactly: G = 1/2
            rec_Phi.append(np.sum(phi, axis=(0, 1)) ** 2)                                 # uniform order-parameter channel
            phimax = max(phimax, float(np.abs(phi).max()))
        if sw % 5 == 0 or sw == sweeps - 1:
            rec_smin.append(m.logdet(S, want_smin=True)[2])
        if sw < burn:                                                 # adapt the step on the acceptance of the last sweep
            rate = (acc - acc_prev) / max(tries - tries_prev, 1)
            step_sig *= 1.3 if rate > 0.6 else (0.7 if rate < 0.35 else 1.0)
            if sw == burn - 1:
                acc, tries = 0, 0
        acc_prev, tries_prev = acc, tries
    th = np.array(rec_theta[burn:])
    z = np.exp(1j * th)
    nb_ = max(2, min(8, len(z) // 4))
    bins = np.array_split(z, nb_)
    means = np.array([b.mean() for b in bins])
    Rm = abs(z.mean()); Rerr = np.std(np.abs(means)) / np.sqrt(nb_)
    c = np.cos(th) - np.cos(th).mean()
    tau_int = 0.5
    if np.var(c) > 1e-14:
        for t in range(1, min(len(c) // 3, 50)):
            rho = np.mean(c[:-t] * c[t:]) / np.var(c)
            if rho < 0.05:
                break
            tau_int += rho
    R_safe = max(Rm, 1e-12)
    res = {"walk": walk, "cell": cell, "N": N, "beta": beta, "dtau": dtau, "g": g, "sweeps": sweeps, "seed": seed,
           "R": Rm, "R_err": Rerr, "minus_log_R": float(-np.log(R_safe)), "var_theta": float(np.var(th)), "delta_f": float(-np.log(R_safe) / (beta * N)),
           "L_eff": float(abs(np.linalg.det(CELLS[cell].astype(float))) ** (1 / 3)), "c_beta": float(beta / abs(np.linalg.det(CELLS[cell].astype(float))) ** (1 / 3)),
           "N_req": float((1 - R_safe ** 2) / (0.01 * R_safe ** 2)), "tau_int": tau_int, "acceptance": acc / max(tries, 1),
           "min_smin": float(np.min(rec_smin)), "seconds": time.time() - t0, "fast": FAST, "logdet_drift": drift}
    if measure and rec_G:
        # reweighted G_k(tau) = <G e^{i theta}>_pq / <e^{i theta}>_pq, binned errors; free reference from the same code path
        th_m = np.array([t for t, _ in rec_G]); z_m = np.exp(1j * th_m)
        Hfree = Model(cell, beta, dtau, 0.0, walk, twist)
        Gfree = Hfree.greens_k(np.zeros((m.Lt, m.N, 10)), ks)
        meas = {}
        nbm = max(2, min(8, len(z_m) // 3))
        if "seeds_P" in rec_G[0][1]:
            Cs = np.array([G["seeds_P"] for _, G in rec_G])                              # (n_meas, L+1, 3, 3)
            Cbar = (Cs * z_m[:, None, None, None]).mean(0) / z_m.mean()
            Cbar = 0.5 * (Cbar + np.conj(np.transpose(Cbar, (0, 2, 1))))                 # Hermitian part
            t0_ = 1
            # GEVP by whitening with the positive part of C(t0) (HS-dressed seeds carry contact terms; directions with
            # non-positive or tiny eigenvalue of C(t0) are dropped and reported)
            lam0, U0 = np.linalg.eigh(Cbar[t0_])
            keep = lam0 > 1e-8 * lam0.max()
            Wt = U0[:, keep] / np.sqrt(lam0[keep])[None, :]
            gev = [sorted(np.linalg.eigvalsh(Wt.conj().T @ Cbar[j_] @ Wt).real, reverse=True) for j_ in range(t0_ + 1, m.Lt // 2 + 1)]
            Eeff = []
            for j_ in range(len(gev) - 1):
                Eeff.append([float(-np.log(a_ / b_) / dtau) if a_ > 0 and b_ > 0 else None for a_, b_ in zip(gev[j_ + 1], gev[j_])])
            res["composites_P"] = {"seeds": ["psi_P (deg 1)", "(2 M^dag psi^dag)_P (deg 3)", "(4 M^dag M psi)_P (deg 5)"],
                                   "C_diag_tau": [[float(Cbar[j_, a_, a_].real) for a_ in range(3)] for j_ in range(m.Lt + 1)],
                                   "C_t0_eigs": lam0.tolist(), "kept_seeds": int(keep.sum()), "gevp_lambda": gev, "E_eff": Eeff}
        if "S_density" in rec_G[0][1]:
            Sd = {k: np.array([G["S_density"][k] for _, G in rec_G]) for k in rec_G[0][1]["S_density"]}
            res["S_density"] = {k: {"mean": float(((v * z_m).mean() / z_m.mean()).real),
                                    "err": float(np.std([((a_ * zz).mean() / zz.mean()).real for a_, zz in zip(np.array_split(v, nbm), np.array_split(z_m, nbm))]) / np.sqrt(nbm))}
                                for k, v in Sd.items()}
        for k in [k_ for k_ in rec_G[0][1] if k_ not in ("seeds_P", "spin", "S_density")]:
            arr = np.array([G[k] for _, G in rec_G])
            num = (arr * z_m[:, None]).mean(0) / z_m.mean()
            bins_k = [(a_ * zz[:, None]).mean(0) / zz.mean() for a_, zz in zip(np.array_split(arr, nbm), np.array_split(z_m, nbm))]
            err = np.std(np.array(bins_k).real, axis=0) / np.sqrt(nbm)
            half = m.Lt // 2
            Gh = num.real[half]
            gap = float(2 / beta * np.arccosh(min(1e6, 1 / (2 * Gh)))) if 0 < Gh < 0.5 else 0.0
            sp = np.array([G["spin"][k] for _, G in rec_G])                              # (n_meas, L+1, 2, 2)
            sp_rw = (sp * z_m[:, None, None, None]).mean(0) / z_m.mean()
            meas[k] = {"C_spin_re": sp_rw.real.tolist(), "C_spin_im": sp_rw.imag.tolist(), "G_tau": num.real.tolist(), "G_tau_imag_max": float(np.abs(num.imag).max()), "G_err": err.tolist(),
                       "G_free": Gfree[k].real.tolist(), "G_half": float(Gh), "G_half_err": float(err[half]),
                       "gap_estimate": gap, "G_half_free": float(Gfree[k].real[half])}
        # symmetry-breaking diagnostic from the HS identity <phi_i phi_j> = 2 g dtau delta_ij + (2 g dtau)^2 <O_i O_j>:
        # chi_A = (1 / (beta N)) int int sum_xy <O_A O_A> = (<Phibar_A^2> - 2 g beta N) / (4 g^2 beta N); finite in the
        # symmetric phase, ~ beta N in a broken phase (O_A = filtered Spin(10) 10-channel pairing composite)
        P2 = np.array(rec_Phi)
        P2m = (P2 * z_m[:, None]).mean(0).real / z_m.mean().real
        chiA = (P2m - 2 * g * beta * m.N) / (4 * g ** 2 * beta * m.N)
        bins_c = [((a_ * zz[:, None]).mean(0).real / zz.mean().real - 2 * g * beta * m.N).mean() / (4 * g ** 2 * beta * m.N)
                  for a_, zz in zip(np.array_split(P2, nbm), np.array_split(z_m, nbm))]
        res["chi_O"] = {"mean": float(chiA.mean()), "err": float(np.std(bins_c) / np.sqrt(nbm)), "per_channel": chiA.tolist()}
        if tmods:
            tres = []
            for it, mt in enumerate(tmods):
                lw = np.array([r_[0] for r_ in rec_T[it]]); th_t = np.array([r_[1] for r_ in rec_T[it]])
                w = np.exp(lw - lw.max()); zt = w * np.exp(1j * th_t)
                ess = float(w.sum() ** 2 / (w ** 2).sum())
                gt = {}
                for k in [k_ for k_ in rec_T[it][0][3] if k_ not in ("seeds_P", "spin", "S_density")]:
                    arr = np.array([r_[3][k] for r_ in rec_T[it]])
                    num_t = (arr * zt[:, None]).sum(0) / zt.sum()
                    bins_t = [(a_ * zz[:, None]).sum(0) / zz.sum() for a_, zz in zip(np.array_split(arr, nbm), np.array_split(zt, nbm))]
                    err_t = np.std(np.array(bins_t).real, axis=0) / np.sqrt(nbm)
                    gt[k] = {"G_tau": num_t.real.tolist(), "G_err": err_t.tolist(), "G_half": float(num_t.real[m.Lt // 2]), "G_half_err": float(err_t[m.Lt // 2])}
                tres.append({"twist": None if mt.twist is None else mt.twist.tolist(), "ESS": ess, "n": len(w), "R_target": float(abs(zt.sum()) / w.sum()),
                             "max_phase_jump": float(max(r_[2] for r_ in rec_T[it])), "greens": gt})
            res["multitwist"] = tres
        res["greens"] = meas
        res["n_meas"] = len(rec_G)
        # raw measurement series for autocorrelation-corrected errors (autocorr.py): phase and G(beta/2) per momentum
        half = m.Lt // 2
        res["series"] = {"theta": th_m.tolist(),
                         "G_half": {k: [float(np.real(G[k][half])) for _, G in rec_G] for k in rec_G[0][1] if k not in ("seeds_P", "spin", "S_density")},
                         "G_half_imag": {k: [float(np.imag(G[k][half])) for _, G in rec_G] for k in rec_G[0][1] if k not in ("seeds_P", "spin", "S_density")}}
        if tmods:
            res["series"]["targets"] = [{"log_ratio": [r_[0] for r_ in rec_T[it]], "theta": [r_[1] for r_ in rec_T[it]],
                                         "G_half": {k: [float(np.real(r_[3][k][half])) for r_ in rec_T[it]] for k in rec_T[it][0][3] if k not in ("seeds_P", "spin", "S_density")}}
                                        for it in range(len(tmods))]
        res["gamma_sentinel"] = sentinel                                               # precision monitor of the stabilised G
        res["phi_max"] = phimax
    if hook is not None:
        res["hook"] = rec_H
    if outer is not None:
        res["outer"] = rec_O; res["outer_acceptance"] = o_acc / max(sweeps, 1)
    if monitor is not None:
        res["monitor"] = rec_M
    if save_phi is not None:
        np.save(save_phi, phi)
    res["burn"] = burn
    if global_moves is not None:
        res["global_acceptance"] = {k: (v[0] / v[1] if v[1] else None) for k, v in g_acc.items()}
    return res


if __name__ == "__main__":
    cell, beta, dtau, g, sweeps, seed = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
    tag = sys.argv[7] if len(sys.argv) > 7 else ""
    walk = sys.argv[8] if len(sys.argv) > 8 else "md"
    measure = len(sys.argv) > 9 and sys.argv[9] == "meas"
    twist = np.array([float(x) for x in sys.argv[10].split(",")]) if len(sys.argv) > 10 and sys.argv[10] != "none" else None
    targets = [np.array([float(x) for x in t.split(",")]) for t in sys.argv[11].split(";")] if len(sys.argv) > 11 else None
    r = run(cell, beta, dtau, g, sweeps, seed, walk, measure, twist, targets)
    if twist is not None:
        r["twist_delta"] = twist.tolist()
    fpath = HERE / f"reweight_points_{walk}_{cell}{tag}.json"
    data = json.loads(fpath.read_text()) if fpath.exists() else []
    data.append(r)
    fpath.write_text(json.dumps(data, indent=2, default=float))
    print(json.dumps(r, default=float))
