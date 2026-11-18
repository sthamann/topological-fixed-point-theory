#!/usr/bin/env python3
"""Part SG: CHIRAL4D.SECTOR_GAP.01 (method) -- symmetry-projected thermal traces with exact signs.

Z_{k,q}(beta) = Tr[P_k P_q e^{-beta H}] finds the lowest state of a whole symmetry sector, whether or not any chosen
operator overlaps with it (dark composites included).  In the HS representation every term is Tr[U X_0 ... X_{L-1}] with
U a symmetry element: translations and Z4 rotations are Gaussian unitaries U = Gamma(u) (u the single-particle unitary).
Their Grassmann representation is c_U = Tr U / 2^n = det(1 + u) / 2^n (exact, complex, no square root) and
T_U = Om^dag (S_U - 1)(S_U + 1)^-1 Om with Nambu image S_U = diag(u, u^*); the Pfaffian composition of qmc_pfaffian.py then
gives the sign-exact trace.  If u has an eigenvalue -1 (c_U = 0, e.g. fermion parity or a translation with e^{ikR} = -1)
the identity U = (U e^{i eta N}) e^{-i eta N} splits it into two Gaussians with det(1 + u) != 0 -- exact, no limit; the two
factors are placed cyclically around the X's (Tr[U_1 U_2 X] = Tr[U_2 X U_1]) so that no partial product has c = 0.
  SG1 brute force (3 modes): Tr[Gamma(u) X_0 X_1 X_2] for random number-conserving unitaries u (generic, and with
      eigenvalues -1 via the eta split) and random positive BdG Gaussians X_l, including the phase.
  SG2 brute force: the Z4 sector projector P_q = 1/4 sum_m e^{-i q m pi/2} Z4^m (Z4 = e^{i pi N/2}) and a translation-type
      permutation projector inserted into products of BdG Gaussians reproduce the Fock-space projected traces.
  SG4 the 32 free Gamma modes: post-hoc deconvolution of total-charge sectors has rank 3 (character zero at m = 2) ->
      KILL; the active charge Z_a = e^{i pi (N - N_Gamma)/2} is inserted directly (brute-force validated).
  SG5 beta plateau, exactly solvable free limit on 16 sites: the slope estimator converges to the exact (enumerated)
      sector gap; the naive estimator is biased by log(degeneracy)/beta.
  SG3 (statement) application plan: zero twist, active charge, beta-plateau slope, comparison with the desired-Weyl
      continuum (SECTOR_EXCESS); cost N_cell x 4 compositions per configuration.
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/sector_gap.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
_src = (HERE / "qmc_pfaffian.py").read_text()
exec(compile(_src.split("# ---------------------------------------------------------------------------------------- Pfaffian Gaussian algebra")[1]
             .split("# ---------------------------------------------------------------------------------------- PF1 brute force")[0], "pfaffian_algebra", "exec"))
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    sys.stdout.flush()


def gauss_unitary(u, Om):
    """Gamma(u) for a number-conserving single-particle unitary u: (log c_U, T_U); requires det(1 + u) != 0."""
    n = u.shape[0]
    S_U = np.block([[u, np.zeros((n, n))], [np.zeros((n, n)), u.conj()]])
    T = Om.conj().T @ (S_U - np.eye(2 * n)) @ np.linalg.inv(S_U + np.eye(2 * n)) @ Om
    sgn, ld = np.linalg.slogdet(np.eye(n) + u)
    return ld + 1j * np.angle(sgn) - n * np.log(2), (T - T.T) / 2


def trace_with(u_list, Ss):
    """Tr[Gamma(u_1) ... Gamma(u_m) X_0 ... X_{L-1}] (complex log), exact."""
    n = Ss[0].shape[0] // 2
    Om = omega_mat(n)
    # cyclic order Tr[U_1 U_2 X] = Tr[U_2 X U_1]: with the eta split U_1 U_2 = Gamma(u) may have c = 0, so the two factors are
    # separated by the X's and the last composition needs only the Pfaffian (no inverse)
    facs = [gauss_unitary(u, Om) for u in u_list[1:]] + [gauss_from_S(S, Om) for S in Ss] + [gauss_unitary(u_list[0], Om)]
    acc = facs[0]
    for k in range(1, len(facs)):
        acc = compose(acc, facs[k], want_T=(k < len(facs) - 1))
    return acc[0] + n * np.log(2)


def eta_split(u, eta=0.37):
    """U = Gamma(u e^{i eta}) Gamma(e^{-i eta}) -- exact factorisation avoiding det(1 + u) = 0."""
    n = u.shape[0]
    return [u * np.exp(1j * eta), np.exp(-1j * eta) * np.eye(n)]


# brute-force Fock space -----------------------------------------------------------------------------------
rng = np.random.default_rng(23)
nb = 3
Zp = np.diag([1.0, -1.0]); low = np.array([[0, 1], [0, 0]], complex)
cf = []
for i in range(nb):
    mm = np.eye(1, dtype=complex)
    for j in range(nb):
        mm = np.kron(mm, Zp if j < i else (low if j == i else np.eye(2)))
    cf.append(mm)
Psi = cf + [x.conj().T for x in cf]
Nop = sum(c.conj().T @ c for c in cf)
def randK(s=1.0):
    h = rng.normal(size=(nb, nb)) + 1j * rng.normal(size=(nb, nb)); h = (h + h.conj().T) / 2
    D = rng.normal(size=(nb, nb)) + 1j * rng.normal(size=(nb, nb)); D = D - D.T
    return s * np.block([[h, D.conj().T], [D, -h.T]])
def op(K):
    return expm(0.5 * sum(K[i, j] * Psi[i].conj().T @ Psi[j] for i in range(2 * nb) for j in range(2 * nb)))
def fock_unitary(u):
    A = -1j * __import__("scipy").linalg.logm(u)                         # u = e^{iA}, A Hermitian
    return expm(1j * sum(A[i, j] * cf[i].conj().T @ cf[j] for i in range(nb) for j in range(nb)))
def rand_unitary():
    Q, R = np.linalg.qr(rng.normal(size=(nb, nb)) + 1j * rng.normal(size=(nb, nb)))
    return Q * (np.diag(R) / np.abs(np.diag(R)))

# SG1 ------------------------------------------------------------------------------------------------------
err1 = []
for trial in range(30):
    Ks = [randK(1.0) for _ in range(3)]
    u = rand_unitary()
    if trial % 3 == 0:                                                    # force an eigenvalue -1
        w, V = np.linalg.eig(u); w[0] = -1.0; u = V @ np.diag(w) @ np.linalg.inv(V)
        Q_, _ = np.linalg.qr(V); u = Q_ @ np.diag(w / np.abs(w)) @ Q_.conj().T
    bf = np.trace(fock_unitary(u) @ np.linalg.multi_dot([op(K) for K in Ks]))
    us = eta_split(u) if abs(np.linalg.det(np.eye(nb) + u)) < 1e-8 else [u]
    lw = trace_with(us, [expm(K) for K in Ks])
    err1.append(abs(np.exp(lw) - bf) / abs(bf))
step("SG1 brute force (3 modes): Tr[Gamma(u) X_0 X_1 X_2] with c_U = det(1 + u)/2^n and T_U from S_U = diag(u, u^*) composed by Pfaffians reproduces the Fock-space trace including phase, also for u with an eigenvalue -1 via the exact split U = Gamma(u e^{i eta}) Gamma(e^{-i eta})",
     max(err1) < 1e-9, f"30 traces (10 with det(1 + u) = 0), max rel. error {max(err1):.1e}")

# SG2 ------------------------------------------------------------------------------------------------------
Ks = [randK(0.8) for _ in range(4)]
Xop = np.linalg.multi_dot([op(K) for K in Ks]); Ss = [expm(K) for K in Ks]
err2 = []
for q in range(4):
    Pq = sum(np.exp(-1j * q * m_ * np.pi / 2) * expm(1j * np.pi / 2 * m_ * Nop) for m_ in range(4)) / 4
    bf = np.trace(Pq @ Xop)
    tot = 0
    for m_ in range(4):
        u = np.exp(1j * np.pi / 2 * m_) * np.eye(nb)
        us = eta_split(u) if abs(np.linalg.det(np.eye(nb) + u)) < 1e-8 else [u]
        tot += np.exp(-1j * q * m_ * np.pi / 2) * np.exp(trace_with(us, Ss)) / 4
    err2.append(abs(tot - bf) / max(abs(bf), 1e-12))
perm = np.eye(nb)[[1, 2, 0]].astype(complex)                              # a 3-cycle "translation"
err3 = []
for kq in range(3):
    Pk = sum(np.exp(-2j * np.pi * kq * r / 3) * fock_unitary(np.linalg.matrix_power(perm, r)) for r in range(3)) / 3
    bf = np.trace(Pk @ Xop)
    tot = sum(np.exp(-2j * np.pi * kq * r / 3) * np.exp(trace_with(eta_split(np.linalg.matrix_power(perm, r)) if abs(np.linalg.det(np.eye(nb) + np.linalg.matrix_power(perm, r))) < 1e-8 else [np.linalg.matrix_power(perm, r)], Ss)) for r in range(3)) / 3
    err3.append(abs(tot - bf) / max(abs(bf), 1e-12))
step("SG2 brute force: the Z4 sector projector P_q (q = 0..3; Z4 = e^{i pi N/2}; the element Z4^2 = (-1)^N has c = 0 and needs the eta split) and a 3-cycle translation projector P_k inserted into a product of four BdG Gaussians reproduce the Fock-space projected traces",
     max(err2) < 1e-9 and max(err3) < 1e-9, f"Z4 sectors max rel. error {max(err2):.1e}; translation sectors {max(err3):.1e}")
# SG4: the 32 free Gamma modes: their Z4 character chi_Gamma(m) = Tr e^{i pi m N_Gamma / 2} = (1 + i^m)^32 vanishes at m = 2,
# so a post-hoc deconvolution of measured total-charge sector traces has rank 3 (one sector combination is lost exactly);
# the active charge Z_a = e^{i pi (N - N_Gamma)/2} must be inserted directly.
chi = [(1 + 1j ** m_) ** 32 for m_ in range(4)]
conv = np.array([[np.exp(-1j * np.pi * q * m_ / 2) * chi[m_] for m_ in range(4)] for q in range(4)]) / 4
rank_conv = np.linalg.matrix_rank(conv @ np.array([[np.exp(1j * np.pi * q * m_ / 2) for q in range(4)] for m_ in range(4)]) / 4, tol=1e-6 * 2 ** 32)
# active-charge projector, brute force: mode 0 plays a decoupled free "Gamma" mode
err4 = []
Ks4 = []
for _ in range(3):
    K = randK(0.9)
    K[0, :] = 0; K[:, 0] = 0; K[nb, :] = 0; K[:, nb] = 0                  # mode 0 decoupled and free (zero energy)
    Ks4.append(K)
X4 = np.linalg.multi_dot([op(K) for K in Ks4]); S4 = [expm(K) for K in Ks4]
Nact = Nop - cf[0].conj().T @ cf[0]
for q in range(4):
    Pa = sum(np.exp(-1j * q * m_ * np.pi / 2) * expm(1j * np.pi / 2 * m_ * Nact) for m_ in range(4)) / 4
    bf = np.trace(Pa @ X4)
    tot = 0
    for m_ in range(4):
        u = np.diag([1.0] + [np.exp(1j * np.pi / 2 * m_)] * (nb - 1)).astype(complex)
        us = eta_split(u) if abs(np.linalg.det(np.eye(nb) + u)) < 1e-8 else [u]
        tot += np.exp(-1j * q * m_ * np.pi / 2) * np.exp(trace_with(us, S4)) / 4
    err4.append(abs(tot - bf) / abs(np.trace(X4)))                     # relative to the full trace (sectors can nearly cancel)
step("SG4 Gamma modes: the Z4 character of the 32 free Gamma modes, (1 + i^m)^32 = [2^32, 2^16, 0, 2^16], vanishes at m = 2, so deconvolving total-charge sector traces afterwards has rank 3 (KILL); the active charge Z_a = e^{i pi (N - N_Gamma)/2} inserted directly as the Gaussian unitary u = diag(1_Gamma, i^m 1_active) reproduces the brute-force active-sector traces",
     abs(chi[2]) == 0 and rank_conv == 3 and max(err4) < 1e-9, f"chi_Gamma = {[f'{abs(x):.3g}' for x in chi]}; deconvolution rank {rank_conv}; active-charge projector max rel. error {max(err4):.1e}")
# SG5: beta-plateau test of the sector-gap estimator in the exactly solvable free limit (cubic16, g = 0)
import importlib.util as _ilu
_sp = _ilu.spec_from_file_location("reweight", HERE / "reweight.py"); _rw = _ilu.module_from_spec(_sp); _sp.loader.exec_module(_rw)
m16 = _rw.Model("cubic16", 1.0, 0.25, 0.0, "md")
Xs = np.array(m16.sites, float); Rv = _rw.CELLS["cubic16"].astype(float)
Gr = 2 * np.pi * np.linalg.inv(Rv).T
def fccq(k):
    return np.allclose(np.cos(np.array([k @ np.array(s_) for s_ in __import__("itertools").product((-1, 1), repeat=3)])), 1, atol=1e-9)
ks = []
for mm_ in __import__("itertools").product(range(-4, 5), repeat=3):
    k = np.array(mm_) @ Gr
    if not any(fccq(k - q_) for q_ in ks):
        ks.append(k)
    if len(ks) == 16:
        break
def kidx(k):
    return next(i for i, q_ in enumerate(ks) if fccq(k - q_))
Bm = _rw.B3.astype(float); Binv = np.linalg.inv(Bm); Sg = [_rw.sx, _rw.sy, _rw.sz]
tauv = [sum(Binv[a_, i] * Sg[a_] for a_ in range(3)) for i in range(3)]
def hk(k):
    p_ = Bm @ k
    return sum((np.sin(p_[i]) + np.cos(p_[(i + 1) % 3]) - 1) * tauv[i] for i in range(3))
eps = {i: np.linalg.eigvalsh(hk(k)) for i, k in enumerate(ks)}            # two bands per k, x 16 flavours
gam = kidx(np.zeros(3))
# exact free sector gap: zero modes (Gamma: no active charge; P: 32 modes at zero energy) are excited for free -- every P
# particle (hole) shifts (K, q) by (+P, +1) ((-P, -1)), any number; the costly quasi-particles (cheapest particle / hole per
# momentum) are enumerated up to three
iP = kidx(np.pi / 2 * np.ones(3))
types = []
for i in range(16):
    if i in (gam, iP):
        continue
    e_ = eps[i]
    types.append((i, +1, float(e_[e_ > 1e-12].min())))
    types.append((i, -1, float(-e_[e_ < -1e-12].max())))
exact_gap = {}
from itertools import combinations_with_replacement as cwr
for nq in range(0, 4):
    for combo in cwr(range(len(types)), nq):
        K0 = np.zeros(3); q0 = 0; cost = 0.0
        for t in combo:
            i, ch, c_ = types[t]
            K0 = K0 + ch * ks[i]; q0 += ch; cost += c_
        for nP in range(4):                                            # free P zero-mode shifts (4P is a reciprocal vector)
            key = (kidx(K0 + nP * ks[iP]), (q0 + nP) % 4)
            exact_gap[key] = min(exact_gap.get(key, np.inf), cost)
import mpmath as _mp
_mp.mp.dps = 60
def logZ_sector(beta_, K_i, q, precise=True):
    """log Tr[P_K P_q^a e^{-beta H_0}] exactly: Tr[Gamma(u) e^{-beta H}] = prod_modes (1 + u e^{-beta eps}) in the momentum basis;
    the symmetry sum cancels to e^{-beta Delta}, so it is evaluated with 60 digits and with EXACT roots of unity for the phases
    (float phases carry 1e-14 errors and destroy the cancellation once e^{-beta Delta} < 1e-14)."""
    I_ = _mp.mpc(0, 1)
    def quarter(x):                                                    # e^{i x} with x in (pi/2) Z exactly (roots of unity)
        n_ = x / (np.pi / 2); assert abs(n_ - round(n_)) < 1e-9
        return I_ ** (int(round(n_)) % 4)
    tot = _mp.mpc(0)
    for R in Xs:
        for m_ in range(4):
            lt = _mp.mpc(0)
            for i, k in enumerate(ks):
                ph = quarter(float(k @ R)) * (I_ ** m_ if i != gam else 1)
                for e_ in eps[i]:
                    lt += 16 * _mp.log(1 + ph * _mp.exp(-beta_ * _mp.mpf(float(e_))))
            tot += _mp.exp(lt) * quarter(-float(ks[K_i] @ R)) * I_ ** ((-q * m_) % 4)
    return complex(_mp.log(tot / (len(Xs) * 4)))
sectors = [(kidx(np.pi / 2 * np.ones(3)), 1), (kidx(np.array([np.pi, 0, 0])), 1), (kidx(-np.pi / 2 * np.ones(3)), 1), (gam, 1), (kidx(np.array([np.pi, 0, 0])), 2)]
ref = (gam, 0)
rows5 = []
for K_i, q in sectors:
    lr = {b_: (logZ_sector(b_, K_i, q) - logZ_sector(b_, *ref)).real for b_ in (4.0, 8.0, 16.0)}
    naive = -lr[16.0] / 16.0
    slope = -(lr[16.0] - lr[8.0]) / 8.0
    rows5.append((K_i, q, exact_gap.get((K_i, q), np.inf), naive, slope, -(lr[8.0] - lr[4.0]) / 4.0))
ok5 = all(abs(sl - ex) < 0.02 for _, _, ex, _, sl, _ in rows5 if np.isfinite(ex))
step("SG5 beta plateau in the exactly solvable free limit (cubic16, g = 0, active Z4 charge; 60-digit arithmetic and exact roots of unity, because the symmetry sum cancels to e^{-beta Delta}): the slope estimator -[log(Z_s/Z_r)(2 beta) - log(Z_s/Z_r)(beta)] / beta converges to the exact sector gap (enumerated from up to three quasi-particles), while the naive -log(Z_s/Z_r)/beta is biased by log(degeneracy)/beta; with float phases the cancellation fails beyond e^{-beta Delta} ~ 1e-14, and in QMC the signal-to-noise falls like e^{-beta Delta}: the sector test is limited to moderate beta Delta",
     ok5, "; ".join(f"K#{K_i} q={q}: exact {ex:.4f}, naive(beta=16) {nv:.4f}, slope(8->16) {sl:.4f}, slope(4->8) {s4:.4f}" for K_i, q, ex, nv, sl, s4 in rows5))

step("SG3 (statement) application plan for Theta_1: zero twist only (the exact Gamma factorisation Z = 2^32 Z_active fails with a twist); the sector is defined by the ACTIVE Z4 charge inserted directly (post-hoc Gamma deconvolution killed, SG4); the HS slice of reweight.py preserves translations and Z4 exactly (symmetric_transfer.py), so no group-averaged integrator is needed; the gap must come from a beta plateau of the slope -d log(Z_s/Z_r)/d beta (degeneracy factor cancels), and be compared with the desired-Weyl multi-particle continuum (SECTOR_EXCESS); per-configuration cost N_cell x 4 Pfaffian compositions: feasible on debug4 / cell8, on 16+ sites only with a faster Pfaffian",
     True, "method validated; physics run pending")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG)}
(HERE / "sector_gap_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
