#!/usr/bin/env python3
"""Part SPG: SECTOR.PROJECTED.GREEN.01 (method) -- remove thermal excited -> excited lines by projecting the start sector.

The thermal Green function C(tau) = Tr[e^{-(beta - tau) H} O e^{-tau H} O^dag] / Z contains transitions from every thermally
populated state.  An excited -> excited transition can sit at omega = 0 while the ground-sector excitation costs a finite
gap.  The sector-projected correlator
    C_{r -> s}(tau) = Tr[P_r e^{-(beta - tau) H} O P_s e^{-tau H} O^dag] / Tr[P_r e^{-beta H}]
starts in the reference sector r (the sector of the ground state, determined from the sector traces, not guessed).
  SPG1 Fock-space toy (two modes, H = 2 n_a + eps_b n_b - 2 n_a n_b, ground state empty): the ordinary C_a(tau) has a
       zero-frequency line whose residue decays like e^{-beta eps_b} on cooling (a thermal line -- GROUNDLINE classifies it
       THERMAL), while the true ground-sector gap for adding a is 2; the projected C_{0 -> 1} has no zero-frequency weight
       and decays with the gap 2 at every beta.
       The zero-line residue series is classified THERMAL by the GROUNDLINE.RESIDUE rule (groundline_residue.py).
  SPG2 the reference sector is identified from the sector traces (lowest free energy at large beta), not assumed.
  SPG3 brute force (3 modes, BdG): the per-configuration Gaussian formula
           Tr[U_m X_2 O U_m' X_1 O^dag] = Tr[Z] a^T S_B (1 + S_Z)^-1 tau b,   Z = U_m X_2 U_m' X_1,  B = U_m' X_1,
       (O = a^T Psi a pure annihilation form, O^dag = b^T Psi, S the Nambu images, tau the c <-> c^dag swap; Tr[Z] by
       the sign-exact Pfaffian composition with the eta split for Z4^2 = (-1)^N), summed with the Z4 characters of P_r and
       P_s, reproduces the Fock-space C_{r -> s} numerator for random non-number-conserving slices (the per-configuration
       HS slices are not Z4-symmetric, so both projectors are needed).  Exact zeros Tr[Z] = 0 need the nilpotent form
       e^{O' O^dag} = 1 + (e^s - 1)/s O' O^dag (not needed generically; flagged).
  SPG4 (statement) QMC implementation on the Theta_1 cells pending: cost 16 compositions per (configuration, tau cut).
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/sector_projected_green.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
_src = (HERE / "qmc_pfaffian.py").read_text()
exec(compile(_src.split("# ---------------------------------------------------------------------------------------- Pfaffian Gaussian algebra")[1]
             .split("# ---------------------------------------------------------------------------------------- PF1 brute force")[0], "pfaffian_algebra", "exec"))
_sg = (HERE / "sector_gap.py").read_text()
exec(compile("def gauss_unitary" + _sg.split("def gauss_unitary")[1].split("# brute-force Fock space")[0], "sector_gap", "exec"))
_gl = (HERE / "groundline_residue.py").read_text()
exec(compile("def classify" + _gl.split("def classify")[1].split("b_ = np.array")[0], "groundline", "exec"))
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    sys.stdout.flush()


Zp = np.diag([1.0, -1.0]); low = np.array([[0, 1], [0, 0]], float)
a = np.kron(low, np.eye(2)); b = np.kron(Zp, low)
na, nb = a.T @ a, b.T @ b
eps_b = 1.5
H = 2.0 * na + eps_b * nb - 2.0 * na @ nb                               # adding a costs 2 from |0>, 0 if b is occupied
Nop = na + nb
def proj(n):
    w, V = np.linalg.eigh(Nop)
    return (V[:, np.abs(w - n) < 1e-9]) @ V[:, np.abs(w - n) < 1e-9].T
def curve(beta, projected, L=16):
    tau = np.linspace(0, beta, L + 1)
    out = []
    for t in tau:
        A = expm(-(beta - t) * H); Bm = expm(-t * H)
        if projected:
            num = np.trace(proj(0) @ A @ a @ proj(1) @ Bm @ a.T); den = np.trace(proj(0) @ expm(-beta * H))
        else:
            num = np.trace(A @ a @ Bm @ a.T); den = np.trace(expm(-beta * H))
        out.append(num / den)
    return tau, np.array(out)
def zero_line_weight(beta):
    """exact spectral decomposition of the ordinary C_a: weight of the omega = 0 line (b occupied -> a and b occupied)."""
    E, V = np.linalg.eigh(H); w = np.exp(-beta * (E - E.min())); Zs = w.sum()
    M = V.T @ a.T @ V                                                   # <n| a^dag |m>
    wt = 0.0
    for m in range(4):
        for n in range(4):
            if abs(E[n] - E[m]) < 1e-9 and abs(M[n, m]) > 1e-12:
                wt += (w[m] + w[n]) * M[n, m] ** 2 / Zs
    return wt
res = {}
for beta in (2.0, 4.0, 6.0):
    t_, Cord = curve(beta, False); t_, Cpr = curve(beta, True)
    # effective energy from the projected curve (single exponential e^{-E tau} for a ground-sector start)
    Epr = -np.log(Cpr[8] / Cpr[0]) / t_[8]
    res[beta] = (zero_line_weight(beta), Epr, float(Cord[8]), float(Cpr[8]))
dec = [res[b_][0] for b_ in (2.0, 4.0, 6.0)]
cl, _ = classify([2.0, 4.0, 6.0], dec, [1e-3 * d for d in dec])
step("SPG1 Fock-space toy: the ordinary thermal C_a carries an omega = 0 line whose residue decays like e^{-beta eps_b} on cooling (thermal excited -> excited line), while the ground-sector cost of adding a is 2; the start-projected C_{0 -> 1} has no zero-frequency weight and its effective energy is the true gap 2 at every beta",
     all(abs(res[b_][1] - 2.0) < 1e-9 for b_ in res) and dec[0] > dec[1] > dec[2] and abs(np.log(dec[1] / dec[2]) / 2.0 - eps_b) < 0.05 and cl == "THERMAL",
     "; ".join(f"beta = {b_:g}: zero-line residue {r[0]:.2e}, projected effective energy {r[1]:.4f}" for b_, r in res.items()) + f"; residue decay rate {np.log(dec[1] / dec[2]) / 2.0:.3f} (eps_b = {eps_b}); GROUNDLINE rule -> {cl}")
# SPG2: the reference sector from sector traces
Fs = {n: float(-np.log(np.trace(proj(n) @ expm(-20.0 * H))) / 20.0) for n in (0, 1, 2)}
ref = min(Fs, key=Fs.get)
step("SPG2 the reference sector is identified from the sector traces (lowest sector free energy at large beta), not assumed",
     ref == 0, f"sector free energies at beta = 20: { {k: round(v, 4) for k, v in Fs.items()} } -> reference N = {ref}")
# SPG3: per-configuration Gaussian formula for the projected Green function, brute force with BdG slices
rng = np.random.default_rng(29)
nb = 3
low_c = low.astype(complex)
cf = []
for i in range(nb):
    mm = np.eye(1, dtype=complex)
    for j in range(nb):
        mm = np.kron(mm, Zp if j < i else (low_c if j == i else np.eye(2)))
    cf.append(mm)
Psi = cf + [x.conj().T for x in cf]
Nf = sum(c.conj().T @ c for c in cf)
tau_sw = np.block([[np.zeros((nb, nb)), np.eye(nb)], [np.eye(nb), np.zeros((nb, nb))]])
def randK(s_=0.8):
    h = rng.normal(size=(nb, nb)) + 1j * rng.normal(size=(nb, nb)); h = (h + h.conj().T) / 2
    D = rng.normal(size=(nb, nb)) + 1j * rng.normal(size=(nb, nb)); D = D - D.T
    return s_ * np.block([[h, D.conj().T], [D, -h.T]])
def fock(K):
    return expm(0.5 * sum(K[i, j] * Psi[i].conj().T @ Psi[j] for i in range(2 * nb) for j in range(2 * nb)))
def trace_cyclic(factors):
    """sign-exact Tr of a cyclic product of Gaussians given as ('S', S) or ('U', u); unitaries with det(1 + u) = 0 are eta-split."""
    Om = omega_mat(nb); facs = []
    for kind, M in factors:
        if kind == "S":
            facs.append(gauss_from_S(M, Om))
        else:
            us = eta_split(M) if abs(np.linalg.det(np.eye(nb) + M)) < 1e-8 else [M]
            facs += [gauss_unitary(u_, Om) for u_ in us]
    acc = facs[0]
    for k in range(1, len(facs)):
        acc = compose(acc, facs[k], want_T=(k < len(facs) - 1))
    return np.exp(acc[0] + nb * np.log(2))
def S_unit(u):
    return np.block([[u, np.zeros((nb, nb))], [np.zeros((nb, nb)), u.conj()]])
err4, scale = [], []
for trial in range(6):
    K2 = [randK() for _ in range(2)]; K1 = [randK() for _ in range(2)]
    X2f, X1f = fock(K2[0]) @ fock(K2[1]), fock(K1[0]) @ fock(K1[1])
    S2, S1 = expm(K2[0]) @ expm(K2[1]), expm(K1[0]) @ expm(K1[1])
    a_ = rng.normal(size=nb) + 1j * rng.normal(size=nb)
    Of = sum(a_[i] * cf[i] for i in range(nb)); av = np.concatenate([a_, np.zeros(nb)]); bv = np.concatenate([np.zeros(nb), a_.conj()])
    for r in range(4):
        for s_ in range(4):
            Pr = sum(np.exp(-1j * r * m * np.pi / 2) * expm(1j * np.pi / 2 * m * Nf) for m in range(4)) / 4
            Ps = sum(np.exp(-1j * s_ * m * np.pi / 2) * expm(1j * np.pi / 2 * m * Nf) for m in range(4)) / 4
            bf = np.trace(Pr @ X2f @ Of @ Ps @ X1f @ Of.conj().T)
            tot = 0
            for m in range(4):
                for mp in range(4):
                    um, ump = np.exp(1j * np.pi / 2 * m) * np.eye(nb, dtype=complex), np.exp(1j * np.pi / 2 * mp) * np.eye(nb, dtype=complex)
                    trZ = trace_cyclic([("S", S2), ("U", ump), ("S", S1), ("U", um)])      # Tr[U_m X_2 U_m' X_1], cyclically
                    SB = S_unit(ump) @ S1; SZ = S_unit(um) @ S2 @ SB
                    tot += np.exp(-1j * (r * m + s_ * mp) * np.pi / 2) / 16 * trZ * (av @ SB @ np.linalg.solve(np.eye(2 * nb) + SZ, tau_sw @ bv))
            err4.append(abs(tot - bf)); scale.append(abs(np.trace(X2f @ Of @ X1f @ Of.conj().T)))
step("SPG3 brute force (3 modes, random non-number-conserving BdG slices, random annihilation form O): the per-configuration Gaussian formula Tr[U_m X_2 O U_m' X_1 O^dag] = Tr[Z] a^T S_B (1 + S_Z)^-1 tau b with sign-exact Pfaffian traces (eta split for (-1)^N), summed with the Z4 characters of P_r and P_s, reproduces the Fock-space projected numerator Tr[P_r X_2 O P_s X_1 O^dag] for all 16 (r, s)",
     max(e / sc for e, sc in zip(err4, scale)) < 1e-9, f"6 configurations x 16 sector pairs, max error relative to the unprojected numerator {max(e / sc for e, sc in zip(err4, scale)):.1e}")
step("SPG4 (statement) QMC implementation of C_{r -> s} on the Theta_1 cells (reference sector from the measured sector traces, 16 compositions per configuration and tau cut, active charge as in sector_gap.py) pending; exact zeros Tr[Z] = 0 would need the nilpotent form e^{O' O^dag} = 1 + (e^s - 1)/s O' O^dag",
     True, "method validated (SPG1-SPG3); production not started")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG)}
(HERE / "sector_projected_green_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
