#!/usr/bin/env python3
"""Part SR: CHIRAL4D.SIGN_REMEDY.01 -- measured fermion phase of the HS decoupling, and a Kramers-symmetric walk Hamiltonian.

Cell: periodic cubic box of side 4 (16 bcc sites, contains the mirror momentum P), 32 modes per site (Weyl spinor x 16),
512 modes, Nambu dimension 1024; beta = 1, four symmetric Trotter slices; ten real HS fields per site and slice,
phi ~ N(0, 2 g dtau); fermion weight W[phi] = Tr prod_l e^{K_w/2} e^{K_phi,l} e^{K_w/2}, with W^2 = det(1 + B) exactly and
the phase of W tracked continuously along phi -> t phi (W(0) = Tr e^{-beta H_walk} > 0).
  SR1 method: W^2 = det(1 + prod e^{K_l}) and the continuous square-root tracking reproduce the exact Fock trace
      (including its phase) on a random 3-mode system with hopping and pairing.
  SR2 measured phases (stabilised UDV determinant, phase tracked along t): for Theta_1 (BDT Hamiltonian walk) and for
      the time-reversal-symmetric walk of SR3 the weights W[phi] carry small but nonzero phases at beta = 1 on the
      16-site cell: a phase problem exists; at this size it is mild (weighted <cos> close to 1).
  SR3 a time-reversal-symmetric walk exists: odd anti-Hermitian bonds h_eps = i b_eps.sigma, b_eps = (eps + eps1 eps2 eps3 n)/8,
      Hermitian, isotropic light node at k = 0, nodal lines lifted by the s1 s2 s3 term, all zeros isolated, every zero
      off the k = 0 class has f >= 1/2; the walk alone is invariant under on-site T = i sigma_y K.
  SR4 structural no-go (corrects the first reading and tightens Part HS): a sign-free Kramers argument needs an
      ON-SITE antiunitary that keeps every HS term O_A(x) individually (inversion or other lattice maps send phi(x) to
      phi(gx) and are not allowed).  Number-preserving Theta: complex conjugation acts on the ten channels of the complex
      16 by the pattern (+,-,+,-,...) with determinant -1, not realisable by Spin(10) x U(1) -> excluded for ANY walk.
      Particle-hole Theta: needs b^dag sigma_i^* b = -sigma_i^* for all three Pauli directions -> no 2 x 2 solution ->
      excluded for any walk whose hopping spans su(2) (every Weyl walk).  Hence no Kramers sign-free decoupling exists
      for the whole class (two-component Weyl walk + Spin(10) 10-channel pairing).
  SR5 (statement) the sign problem of Gate B is structural, not an artefact of the BDT ordering; QMC needs reweighting
      (feasible only while the phase stays mild), a different interaction channel or representation (Majorana-class
      criteria), or tensor networks.
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/sign_remedy.py
"""
from __future__ import annotations

import json
import sys
from itertools import product
from pathlib import Path

import numpy as np
from scipy.linalg import expm, null_space
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    sys.stdout.flush()


rng = np.random.default_rng(81)
wrap = lambda a: (a + np.pi) % (2 * np.pi) - np.pi


def track_phase(Dfun, tol=np.pi / 4, min_dt=1e-5):
    """continuous phase of D(t) on [0, 1] (adaptive); returns total phase and number of evaluations."""
    t, ph_prev, Phi, dt, nev = 0.0, np.angle(Dfun(0.0)), 0.0, 0.1, 1
    while t < 1 - 1e-12:
        dt = min(dt, 1 - t)
        ph_new = np.angle(Dfun(t + dt)); nev += 1
        inc = wrap(ph_new - ph_prev)
        if abs(inc) > tol and dt > min_dt:
            dt /= 2
            continue
        Phi += inc; ph_prev = ph_new; t += dt; dt *= 1.5
    return Phi, nev


# SR1 ------------------------------------------------------------------------------------------------------
Zp = np.diag([1.0, -1.0]); I2 = np.eye(2); low = np.array([[0, 1], [0, 0]], complex)
def fock_ops(n):
    out = []
    for i in range(n):
        m = np.eye(1, dtype=complex)
        for j in range(n):
            m = np.kron(m, Zp if j < i else (low if j == i else I2))
        out.append(m)
    return out
nt = 3
cops = fock_ops(nt)
Qs, Ks = [], []
for _ in range(4):
    hh = rng.normal(size=(nt, nt)) + 1j * rng.normal(size=(nt, nt)); hh = (hh + hh.conj().T) / 2; hh -= np.trace(hh) / nt * np.eye(nt)
    Mm = rng.normal(size=(nt, nt)) + 1j * rng.normal(size=(nt, nt)); Mm = Mm - Mm.T
    Q = sum(hh[i, j] * cops[i].conj().T @ cops[j] for i in range(nt) for j in range(nt))
    Pp = sum(Mm[i, j] * cops[i] @ cops[j] for i in range(nt) for j in range(nt))
    Qs.append(0.25 * (Q + Pp + Pp.conj().T)); Ks.append(0.25 * np.block([[hh, 2 * Mm.conj().T], [2 * Mm, -hh.T]]))
h0 = rng.normal(size=(nt, nt)) + 1j * rng.normal(size=(nt, nt)); h0 = (h0 + h0.conj().T) / 2; h0 -= np.trace(h0) / nt * np.eye(nt)
Q0 = sum(h0[i, j] * cops[i].conj().T @ cops[j] for i in range(nt) for j in range(nt)); K0 = np.block([[h0, 0 * h0], [0 * h0, -h0.T]])
def W_exact(t):
    return np.trace(expm(-Q0) @ np.linalg.multi_dot([expm(t * Q) for Q in Qs]))
def D_small(t):
    return np.linalg.det(np.eye(2 * nt) + expm(-K0) @ np.linalg.multi_dot([expm(t * K) for K in Ks]))
Phi_s, _ = track_phase(D_small)
W1 = W_exact(1.0)
Wt = np.sqrt(abs(D_small(1.0))) * np.exp(0.5j * Phi_s)
step("SR1.method: W^2 = det(1 + prod e^{K_l}) and the continuous square-root tracking from t = 0 reproduce the exact Fock trace including its phase (random 3-mode hopping + pairing products)",
     abs(W_exact(0.0).imag) < 1e-9 and W_exact(0.0).real > 0 and abs(Wt - W1) / abs(W1) < 1e-8, f"exact {W1:.6g}, tracked {Wt:.6g}")

# lattice ---------------------------------------------------------------------------------------------------
Lb = 4
sites = [x for x in product(range(Lb), repeat=3) if len({c % 2 for c in x}) == 1]
idx = {x: i for i, x in enumerate(sites)}
ns = len(sites)
signs = [np.array(s) for s in product((-1, 1), repeat=3)]
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1.0 + 0j, -1]); S = [sx, sy, sz]
def A_eps(e):
    return (I2 + e[0] * sx) / 2 @ (I2 + e[1] * sy) / 2 @ (I2 + e[2] * sz) / 2
h_bdt = {tuple(e): 0.5j * (A_eps(e) - A_eps(-e).conj().T) for e in signs}
nhat = np.ones(3) / np.sqrt(3)
def h_ct(eta):
    return {tuple(e): 1j * sum(((e + eta * e[0] * e[1] * e[2] * nhat) / 8)[i] * S[i] for i in range(3)) for e in signs}
F = np.eye(ns)
for x in sites:
    for e in signs:
        y = tuple((x[c] + e[c]) % Lb for c in range(3))
        F[idx[x], idx[y]] -= 1 / 8
def lattice_hop(hb):
    H = np.zeros((ns * 32, ns * 32), complex)
    for x in sites:
        for e in signs:
            y = tuple((x[c] + e[c]) % Lb for c in range(3))
            H[32 * idx[y]:32 * idx[y] + 32, 32 * idx[x]:32 * idx[x] + 32] += np.kron(hb[tuple(e)], np.eye(16))
    return H
# Spin(10) pairing matrices (as in Parts SG/FP)
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
MA = [np.kron(np.array([[0, 1], [-1, 0]]), Pe.T @ (Cc @ Gm[A]) @ Pe) for A in range(10)]
FF = [np.outer(F[i], F[i]) for i in range(ns)]                      # rho_x rho_x structure (F real symmetric)
def K_phi(phi):
    """Nambu matrix of sum_{x,A} phi_{xA} O_A(x)."""
    M = np.zeros((ns * 32, ns * 32), complex)
    for i in range(ns):
        M += np.kron(FF[i], sum(phi[i, A] * MA[A] for A in range(10)))
    Z = np.zeros_like(M)
    return np.block([[Z, 2 * M.conj().T], [2 * M, Z]])
dtau, Lt = 0.25, 4
from scipy.linalg import qr


def stable_phase(Ms):
    """phase of det(1 + M_0 M_1 ... M_{L-1}) via a pivoted-QR UDV chain (Loh-Gubernatis splitting)."""
    n = Ms[0].shape[0]
    U, D, V = np.eye(n, dtype=complex), np.ones(n), np.eye(n, dtype=complex)
    for M in reversed(Ms):
        X = (M @ U) * D[None, :]
        Q, R, piv = qr(X, pivoting=True)
        d = np.abs(np.diag(R))
        V = (R / d[:, None])[:, np.argsort(piv)] @ V
        U, D = Q, d
    Db, Ds = np.maximum(D, 1.0), np.minimum(D, 1.0)
    inner = (U.conj().T @ np.linalg.inv(V)) / Db[:, None] + np.diag(Ds)
    return np.angle(np.linalg.slogdet(U)[0]) + np.angle(np.linalg.slogdet(inner)[0]) + np.angle(np.linalg.slogdet(V)[0])


def measured(hb, g, ncfg):
    Hw = lattice_hop(hb)
    Kw = -dtau * np.block([[Hw, np.zeros_like(Hw)], [np.zeros_like(Hw), -Hw.T]])
    Wh = expm(Kw / 2)
    phases = []
    for _ in range(ncfg):
        facs = []
        for l in range(Lt):
            phi = rng.normal(scale=np.sqrt(2 * g * dtau), size=(ns, 10))
            lam, V = np.linalg.eigh(K_phi(phi))
            facs.append((Wh @ V, lam, V.conj().T @ Wh))
        Dph = lambda t: np.exp(1j * stable_phase([(L_ * np.exp(t * lam)) @ R_ for L_, lam, R_ in facs]))
        Phi, _ = track_phase(Dph)
        phases.append(wrap(Phi / 2))
    return np.array(phases)


eta = 1.0
hct = h_ct(eta)
res = {}
for name, hb in (("Theta_1 BDT", h_bdt), ("TR-symmetric walk", hct)):
    for g in (0.02, 0.1):
        ph = measured(hb, g, 3)
        res[(name, g)] = (float(np.max(np.abs(ph))), float(np.mean(np.cos(ph))))
step("SR2.measured phases (stabilised UDV determinant, tracked along t, beta = 1, 16-site cell with P): the weights carry nonzero phases for Theta_1 and for the time-reversal-symmetric walk -- a phase problem exists, mild at this size",
     all(v[0] > 1e-6 for v in res.values()) and all(v[1] > 0.5 for v in res.values()),
     "; ".join(f"{k[0]}, g = {k[1]}: max |phase| {v[0]:.2e}, <cos> {v[1]:.4f}" for k, v in res.items()))

# SR3 ------------------------------------------------------------------------------------------------------
herm = all(np.allclose(hct[tuple(-e)], hct[tuple(e)].conj().T) for e in signs)
odd = all(np.allclose(hct[tuple(-e)], -hct[tuple(e)]) for e in signs)
Hk = lambda k, hb: sum(np.exp(-1j * e @ k) * hb[tuple(e)] for e in signs)
dvec = lambda k: np.array([np.real(np.trace(Hk(k, hct) @ S[i])) / 2 for i in range(3)])
lin = all(np.allclose(dvec(1e-4 * d) / 1e-4, d, atol=1e-6) for d in np.eye(3))
f = lambda k: 1 - np.cos(k[0]) * np.cos(k[1]) * np.cos(k[2])
N = 24
grid = 2 * np.pi * np.arange(N) / N - np.pi
nrm = np.array([[[np.linalg.norm(dvec(np.array([a, b, c]))) for c in grid] for b in grid] for a in grid])
cands = []
for i, j, l in product(range(N), repeat=3):
    v = nrm[i, j, l]
    nb = [nrm[(i + a) % N, (j + b_) % N, (l + c) % N] for a, b_, c in product((-1, 0, 1), repeat=3) if (a, b_, c) != (0, 0, 0)]
    if v <= min(nb) and v < 0.3:
        r = minimize(lambda k: np.linalg.norm(dvec(k)), np.array([grid[i], grid[j], grid[l]]), method="Nelder-Mead", options={"xatol": 1e-12, "fatol": 1e-14, "maxiter": 6000})
        if r.fun < 1e-7:
            cands.append(np.mod(r.x + np.pi, 2 * np.pi) - np.pi)
zeros = []
for p_ in cands:
    if not any(np.allclose(np.mod(p_ - q + np.pi, 2 * np.pi) - np.pi, 0, atol=1e-5) for q in zeros):
        zeros.append(p_)
def chir(k0):
    J = np.array([(dvec(k0 + 1e-5 * e) - dvec(k0 - 1e-5 * e)) / 2e-5 for e in np.eye(3)]).T
    return np.sign(np.linalg.det(J)), abs(np.linalg.det(J))
chis = [chir(z) for z in zeros]
recip = lambda k: np.allclose(np.cos(np.array([k @ s_ for s_ in signs])), 1, atol=1e-7)
others = [z for z in zeros if not recip(z)]
fmin = min(f(z) for z in others) if others else np.inf
isolated = all(c[1] > 1e-6 for c in chis)
line_lifted = min(np.linalg.norm(dvec(np.array([np.pi / 2, np.pi / 2, t]))) for t in np.linspace(0.1, np.pi - 0.1, 9)) > 1e-3
HwC = lattice_hop(hct)
KwC = np.block([[HwC, np.zeros_like(HwC)], [np.zeros_like(HwC), -HwC.T]])
bT = np.kron(np.eye(ns), np.kron(1j * sy, np.eye(16)))                 # on-site T = i sigma_y K, number-preserving
UT = np.block([[bT, np.zeros_like(bT)], [np.zeros_like(bT), bT.conj()]])
walk_T = np.allclose(UT.conj().T @ KwC.conj() @ UT, KwC)
step("SR3.a time-reversal-symmetric walk exists: h_eps = i b_eps.sigma with b_eps = (eps + eps1 eps2 eps3 n)/8 is Hermitian and odd, light node isotropic (d = k + O(k^3)), nodal lines lifted, all zeros isolated, every zero off the k = 0 class has f >= 1/2, chiralities sum to zero; the walk alone is invariant under on-site T = i sigma_y K",
     herm and odd and lin and line_lifted and isolated and fmin >= 0.5 and abs(sum(c[0] for c in chis)) < 1e-9 and walk_T,
     f"{len(zeros)} zeros in the cubic zone, min f off the k = 0 class {fmin:.3f}")

# SR4 ------------------------------------------------------------------------------------------------------
Bm = [Pe.T @ (Cc @ Gm[A]) @ Pe for A in range(10)]
pattern = np.array([1 if np.allclose(b.conj(), b) else -1 if np.allclose(b.conj(), -b) else 0 for b in Bm])
det_pattern = int(np.prod(pattern))
# a unitary Y on the 16 preserving span{B_A} acts by Y^T B_A Y = c (R B)_A with R in SO(10) (Spin(10) x U(1)); number-preserving
# Theta needs Y^T B_A^* Y = c' B_A, i.e. R = c'' diag(pattern) with real orthogonal R -> c'' = +-1 and det R = det(pattern) = -1
T_action = []
for A in range(10):
    KA = K_phi(np.outer(np.eye(ns)[0], np.eye(10)[A])); KT = UT.conj().T @ KA.conj() @ UT
    T_action.append(1 if np.allclose(KT, KA) else -1 if np.allclose(KT, -KA) else 0)
pair_T_fails = any(t_ != 1 for t_ in T_action) and T_action == [int(x) for x in pattern]
ct_null = null_space(np.vstack([np.kron(sm.conj(), I2) + np.kron(I2, sm.conj().T) for sm in S]))
# inversion is not admissible: a particle-hole Theta combined with x -> -x keeps the walk but maps O_A(x) to O_A(-x)
inv = np.zeros((ns, ns))
for x in sites:
    inv[idx[tuple((-c) % Lb for c in x)], idx[x]] = 1
b_inv = np.kron(inv, np.kron(1j * np.eye(2), np.eye(16)))
U_inv = np.block([[np.zeros_like(b_inv), b_inv], [b_inv.conj(), np.zeros_like(b_inv)]])
x0 = [i for i, x in enumerate(sites) if x != tuple((-c) % Lb for c in x)][0]
xm = idx[tuple((-c) % Lb for c in sites[x0])]
Kx = K_phi(np.outer(np.eye(ns)[x0], np.eye(10)[0])); Kmx = K_phi(np.outer(np.eye(ns)[xm], np.eye(10)[0]))
inv_keeps_walk = np.allclose(U_inv.conj().T @ KwC.conj() @ U_inv, KwC)
inv_moves = np.allclose(U_inv.conj().T @ Kx.conj() @ U_inv, Kmx) and not np.allclose(U_inv.conj().T @ Kx.conj() @ U_inv, Kx)
step("SR4.structural no-go: a Kramers argument needs an ON-SITE antiunitary keeping every HS term O_A(x) (lattice maps send phi(x) to phi(gx)); number-preserving Theta fails because conjugation acts on the ten channels of the complex 16 with determinant -1 (not in Spin(10) x U(1)); particle-hole Theta needs b with s^* b = -b s^* for all three Pauli directions, which has no solution -> no Kramers sign-free decoupling for ANY two-component Weyl walk with the Spin(10) 10-channel pairing",
     det_pattern == -1 and pair_T_fails and ct_null.shape[1] == 0 and inv_keeps_walk and inv_moves,
     f"conjugation pattern on B_A {pattern.tolist()} (det {det_pattern}); T acts on the ten channels as {T_action} (= conjugation pattern); particle-hole spin solutions: {ct_null.shape[1]}; inversion-type Theta keeps the walk but sends O_A(x) to O_A(-x): {inv_moves}")
step("SR5.(statement) the sign problem of Gate B is structural for this class, not an artefact of the BDT ordering; QMC needs reweighting (only while the phase stays mild), a different channel / representation (Majorana-class criteria), or tensor networks; Theta_1 stays the candidate, no walk is swapped",
     True, "QMC not started")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "measured_phases": {str(k): v for k, v in res.items()},
       "tr_walk_zeros": len(zeros), "tr_walk_fmin": fmin, "conjugation_pattern": pattern.tolist()}
(HERE / "sign_remedy_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
