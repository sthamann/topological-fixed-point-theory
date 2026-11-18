#!/usr/bin/env python3
"""Part M: MARK.PARITY.01 -- fermion parity of the mark Dirac algebra, the reference extension, and the joint sheet-mass test.

External review (2026-10-01): if the four gamma_a are PHYSICAL Majorana zero modes, their Fock space C^4 has parity
Pi_M = -gamma_0 gamma_1 gamma_2 gamma_3 (= Gamma of Part K); the kinetic alpha_i are even but every Dirac mass is odd.
  M1 R1 obstruction: all Hermitian M with {M, alpha_i} = 0 form the real 2-space {a beta + b i chi beta}
     (chi = -i alpha_1 alpha_2 alpha_3 = +-Pi_M), and every such M is parity-ODD: no even Dirac mass on the isolated
     4-Majorana Fock space (witness: e^{-i pi beta/4} sends an even state to odd with probability 1/2).
  M2 R1 repair: with a reference pair eta, zeta (six Majoranas) Gamma_a = i gamma_a eta satisfy Cl(4); an isometry
     V: C^4 -> fixed-total-parity sector intertwines Gamma_a V = V gamma_a, so every word is preserved and the even
     mass beta~ = i beta eta restricts EXACTLY to beta: all spectral results of Parts K/L/N/H transfer verbatim; the
     mass pulse changes total parity with probability 0.
  M3 price of R1: V^dag n_ref V = (1 +- chi)/2, so a reference splitting Delta adds the chiral potential (Delta/2) chi;
     H = alpha.k + m beta + mu chi has E^2 = m^2 + (|k| +- mu)^2 (not Lorentz unless mu = 0); freezing the reference
     kills the mass (P beta~ P = 0).  R1 therefore needs an exactly degenerate (zero-mode) reference pair.
  M4 mass classes on the walk: for EVERY uniform anticommuting mass pulse M (single or symmetric) U_M(+-P)^2 = I with
     spectrum {1,1,-1,-1}: the doubler stays (agrees with L6).  A uniform chiral phase e^{-i mu chi} removes +-1 at P
     but leaves exact zero modes elsewhere (U(mu,0,0) has eigenvalue 1 twice).  AX.SHEETMASS as implemented
     (causal_dirac L7): C(x) = exp(-(i/2)(m + m_s s(x)) beta), U = C W C, s(x) = bcc sublattice (sheet) sign, i.e. a
     coupling of k and k + Q.  JOINT test at the SAME parameters (m, m_s) = (0.2, 0.6): gapped at quasi-energy 0 and
     pi over the whole zone, Dirac limit omega^2 = k^2 + m_eff^2 isotropic near k = 0, and exponential clustering.
  M5 R1 origin condition: an odd zero mode gamma with [H, gamma] = 0 on the UNCHANGED seam space would make gamma Omega
     an odd ground state, contradicting the unique even vacuum; bound zero modes need a defect sector (they add no
     c = 1/2 propagating branch, unlike a seventeenth chiral species).
  M6 R2 (coefficient reading) has no parity problem: second-quantised bilinears psi^dag M psi commute with (-1)^F for
     ANY single-particle M (checked on the 4-mode Fock space), and the spin lift acts as -I on all 16 internal
     components of psi in C^4 (x) C^16 (direct sum): every species gets spin 1/2 -- the (-1)^16 cancellation of K7
     concerns only species-diagonal Fock lifts.
  M7 (statement) AX.MARKMAJORANA is replaced by AX.MARKCLIFFORD (R2): the single-particle coin space of the propagating
     field is an irreducible, internal-singlet Cl(4)-module whose generators are labelled by the four marks.  R1 is
     an alternative realisation that needs six zero modes, a degenerate reference and a defect sector.
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/mark_parity.py
"""
from __future__ import annotations

import json
import sys
from itertools import product
from pathlib import Path

import numpy as np
from scipy.linalg import eigh as sp_eigh, expm, null_space

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


Zp = np.diag([1.0, -1.0]).astype(complex); I2 = np.eye(2, dtype=complex)
lower = np.array([[0, 1], [0, 0]], dtype=complex)


def majoranas(nmodes):
    def op(i, m):
        out = np.eye(1, dtype=complex)
        for j in range(nmodes):
            out = np.kron(out, Zp if j < i else (m if j == i else I2))
        return out
    out = []
    for i in range(nmodes):
        a = op(i, lower)
        out += [a + a.conj().T, -1j * (a - a.conj().T)]
    return out


w = np.array([[-.5, -.5, -.5], [-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]])
g = majoranas(2)
I4 = np.eye(4)
beta = sum(g) / 2
e = [sum(w[a, i] * g[a] for a in range(4)) for i in range(3)]
alpha = [1j * beta @ e[i] for i in range(3)]
PiM = -g[0] @ g[1] @ g[2] @ g[3]
chi = -1j * alpha[0] @ alpha[1] @ alpha[2]
acomm = lambda x, y: x @ y + y @ x

# M1 ---------------------------------------------------------------------------------------------------
basis = []
for i, j in product(range(4), repeat=2):
    E = np.zeros((4, 4), complex); E[i, j] = 1; basis.append(E)
rows = []
for E in basis:
    rows.append(np.concatenate([acomm(E, a).ravel() for a in alpha]))
ns = null_space(np.array(rows).T)
herm = []
for c in range(ns.shape[1]):
    Mc = sum(ns[k, c] * basis[k] for k in range(16))
    for Mh in ((Mc + Mc.conj().T) / 2, (Mc - Mc.conj().T) / 2j):
        if np.linalg.norm(Mh) > 1e-9:
            herm.append(Mh.ravel())
H_span = np.array(herm).T
rank_real = np.linalg.matrix_rank(np.vstack([H_span.real, H_span.imag]), tol=1e-9)
in_span = np.linalg.matrix_rank(np.column_stack([np.vstack([H_span.real, H_span.imag]),
                                                np.concatenate([beta.ravel().real, beta.ravel().imag]),
                                                np.concatenate([(1j * chi @ beta).ravel().real, (1j * chi @ beta).ravel().imag])]), tol=1e-9) == rank_real
all_odd = all(np.allclose(acomm(PiM, h.reshape(4, 4)), 0) for h in herm)
kin_even = all(np.allclose(PiM @ a, a @ PiM) for a in alpha)
ev, V = np.linalg.eigh(PiM)
even_state = V[:, ev > 0][:, 0]
out_state = expm(-1j * np.pi * beta / 4) @ even_state
p_odd = np.linalg.norm(V[:, ev < 0].conj().T @ out_state) ** 2
step("M1.R1 obstruction: the Hermitian masses anticommuting with all alpha_i form the real 2-space {a beta + b i chi beta}; all are parity-odd while the kinetic alpha_i are even; witness e^{-i pi beta/4}: even -> odd with probability 1/2",
     rank_real == 2 and in_span and all_odd and kin_even and (np.allclose(chi, PiM) or np.allclose(chi, -PiM)) and abs(p_odd - 0.5) < 1e-12,
     f"real dimension of the mass space {rank_real}; chi = {'+' if np.allclose(chi, PiM) else '-'}Pi_M; odd probability {p_odd:.3f}")

# M2 ---------------------------------------------------------------------------------------------------
g6 = majoranas(3)
gam6, eta, zeta = g6[:4], g6[4], g6[5]
Gam6 = [1j * x @ eta for x in gam6]
Pi_tot = (-1j) ** 3 * np.linalg.multi_dot(g6)                   # (-i)^3 g1...g6 = total parity on 3 modes
Pi_tot = Pi_tot if np.allclose(Pi_tot @ Pi_tot, np.eye(8)) else 1j * Pi_tot
clif6 = all(np.allclose(acomm(Gam6[a], Gam6[b]), 2 * (a == b) * np.eye(8)) for a in range(4) for b in range(4))
bil = all(np.allclose(Gam6[a] @ Gam6[b], gam6[a] @ gam6[b]) for a in range(4) for b in range(4))
best = None
for sector in (-1, 1):
    Ps = (np.eye(8) + sector * Pi_tot) / 2
    eqs = []
    basV = []
    for i, j in product(range(8), range(4)):
        E = np.zeros((8, 4), complex); E[i, j] = 1; basV.append(E)
    for E in basV:
        eqs.append(np.concatenate([(Gam6[a] @ E - E @ g[a]).ravel() for a in range(4)] + [((np.eye(8) - Ps) @ E).ravel()]))
    nsV = null_space(np.array(eqs).T)
    if nsV.shape[1] >= 1:
        Vm = sum(nsV[k, 0] * basV[k] for k in range(32))
        Vm = Vm / np.sqrt(np.real(np.trace(Vm.conj().T @ Vm)) / 4)
        best = (sector, Vm, nsV.shape[1])
        break
sector, Vm, mult = best
iso = np.allclose(Vm.conj().T @ Vm, I4)
inter = all(np.allclose(Gam6[a] @ Vm, Vm @ g[a]) for a in range(4))
beta6 = sum(gam6) / 2
beta_t = 1j * beta6 @ eta
even_mass = np.allclose(beta_t @ Pi_tot, Pi_tot @ beta_t)
restr = np.allclose(Vm.conj().T @ beta_t @ Vm, beta)
psi_in = Vm @ even_state
psi_out = expm(-1j * np.pi * beta_t / 4) @ psi_in
p_change = np.linalg.norm((np.eye(8) - sector * Pi_tot) / 2 @ psi_out) ** 2
step("M2.R1 repair: with a reference pair (eta, zeta) the Gamma_a = i gamma_a eta satisfy Cl(4) with unchanged even bilinears; an isometry V onto a fixed-total-parity sector intertwines Gamma_a V = V gamma_a; the even mass i beta eta restricts EXACTLY to beta (all walk spectra transfer); the mass pulse changes total parity with probability 0",
     clif6 and bil and iso and inter and mult == 1 and even_mass and restr and p_change < 1e-14,
     f"sector Pi_tot = {sector:+d}, intertwiner unique ({mult}); parity-change probability {p_change:.1e}")

# M3 ---------------------------------------------------------------------------------------------------
Pi_ref = -1j * eta @ zeta
n_ref = (np.eye(8) - Pi_ref) / 2
red = Vm.conj().T @ n_ref @ Vm
chiral_pot = np.allclose(red, (I4 + chi) / 2) or np.allclose(red, (I4 - chi) / 2)
mu, m = 0.17, 0.11
disp_ok = True
for k in np.random.default_rng(2).normal(size=(5, 3)):
    Hc = sum(k[i] * alpha[i] for i in range(3)) + m * beta + mu * chi
    E = np.sort(np.linalg.eigvalsh(Hc))
    kk = np.linalg.norm(k)
    pred = np.sort([s * np.sqrt(m ** 2 + (kk + t * mu) ** 2) for s in (1, -1) for t in (1, -1)])
    disp_ok &= np.allclose(E, pred)
Pfro = (I4 + chi) / 2
frozen = np.allclose(Pfro @ beta @ Pfro, 0)
step("M3.price of R1: V^dag n_ref V = (1 +- chi)/2, so a reference splitting adds the chiral potential (Delta/2) chi; alpha.k + m beta + mu chi has E^2 = m^2 + (|k| +- mu)^2 (Lorentz only for mu = 0); a frozen reference kills the mass (P beta P = 0): R1 needs an exactly degenerate reference pair",
     chiral_pot and disp_ok and frozen, "condition: Delta = 0 (zero-mode reference)")

# M4 ---------------------------------------------------------------------------------------------------
W = lambda k: (np.cos(k[0]) * I4 - 1j * np.sin(k[0]) * alpha[0]) @ (np.cos(k[1]) * I4 - 1j * np.sin(k[1]) * alpha[1]) @ (np.cos(k[2]) * I4 - 1j * np.sin(k[2]) * alpha[2])
Pc = np.pi / 2 * np.ones(3)
rng = np.random.default_rng(9)
class_ok = True
for _ in range(6):
    a_, b_ = rng.normal(size=2)
    M = a_ * beta + b_ * 1j * chi @ beta
    for U in (expm(-1j * M) @ W(Pc), expm(-0.5j * M) @ W(Pc) @ expm(-0.5j * M), expm(-1j * M) @ W(-Pc)):
        class_ok &= np.allclose(U @ U, I4) and np.allclose(np.sort(np.linalg.eigvals(U).real), [-1, -1, 1, 1], atol=1e-9)
muc = np.pi / 4
Uchi = lambda k: expm(-1j * muc * chi) @ W(k)
corner_free = np.min(np.abs(np.abs(np.linalg.eigvals(Uchi(Pc)).real) - 1)) > 1e-3 or np.min(np.abs(np.angle(np.linalg.eigvals(Uchi(Pc))))) > 1e-3
zeros_elsewhere = np.sum(np.abs(np.angle(np.linalg.eigvals(Uchi(np.array([muc, 0, 0]))))) < 1e-9)
Q = np.pi * np.ones(3)
sx = np.array([[0, 1], [1, 0]])
def U8(k, mm, mms):
    Wk = np.zeros((8, 8), complex); Wk[:4, :4], Wk[4:, 4:] = W(k), W(k + Q)
    Cc = expm(-0.5j * (mm * np.kron(np.eye(2), beta) + mms * np.kron(sx, beta)))
    return Cc @ Wk @ Cc
mJ, msJ = 0.2, 0.6
g0 = gpi = np.inf
gr = 2 * np.pi * np.arange(20) / 20 - np.pi
for k in list(product(gr, repeat=3)) + [tuple(x) for x in (np.zeros(3), Pc, -Pc, np.array([np.pi, 0, 0]))]:
    q = np.angle(np.linalg.eigvals(U8(np.array(k), mJ, msJ)))
    g0 = min(g0, np.min(np.abs(q))); gpi = min(gpi, np.min(np.pi - np.abs(q)))
dirs = rng.normal(size=(10, 3)); dirs /= np.linalg.norm(dirs, axis=1)[:, None]
om0 = np.min(np.abs(np.angle(np.linalg.eigvals(U8(np.zeros(3), mJ, msJ)))))
def mass_fit(kmag):
    vals = []
    for d in dirs:
        q = np.abs(np.angle(np.linalg.eigvals(U8(kmag * d, mJ, msJ))))
        vals.append(np.min(q) ** 2 - kmag ** 2)
    return np.array(vals)
f1, f2 = mass_fit(0.02), mass_fit(0.01)
iso1, iso2 = np.ptp(f1) / om0 ** 2, np.ptp(f2) / om0 ** 2
Lb = 16
sites = [x for x in product(range(Lb), repeat=3) if len({c % 2 for c in x}) == 1]
idx = {x: i for i, x in enumerate(sites)}
signs = [np.array(s) for s in product((-1, 1), repeat=3)]
Aeps = {tuple(s): (I4 + s[0] * alpha[0]) / 2 @ (I4 + s[1] * alpha[1]) / 2 @ (I4 + s[2] * alpha[2]) / 2 for s in signs}
n4 = 4 * len(sites)
Wr = np.zeros((n4, n4), complex); Cd = np.zeros_like(Wr)
for x, i in idx.items():
    for s_, As in Aeps.items():
        y = tuple((x[c] - s_[c]) % Lb for c in range(3))
        Wr[4 * i:4 * i + 4, 4 * idx[y]:4 * idx[y] + 4] += As
    Cd[4 * i:4 * i + 4, 4 * i:4 * i + 4] = expm(-0.5j * (mJ + msJ * (1 if x[0] % 2 == 0 else -1)) * beta)
Ur = Cd @ Wr @ Cd
evh, Vh = sp_eigh((Ur - Ur.conj().T) / 2j, driver="evr")
Vn = Vh[:, evh < 0]
o = idx[(0, 0, 0)]
cvals = np.array([np.linalg.norm(Vn[4 * idx[(r, r, r)]:4 * idx[(r, r, r)] + 4] @ Vn[4 * o:4 * o + 4].conj().T) for r in (1, 3, 5)])
slope = np.polyfit([1, 3, 5], np.log(cvals), 1)[0]
step("M4.mass classes: every uniform anticommuting mass pulse (single or symmetric) leaves U(+-P)^2 = I with spectrum {1,1,-1,-1}; a uniform chiral phase e^{-i pi/4 chi} frees P but has two exact zero quasi-energies at (pi/4,0,0); AX.SHEETMASS (C(x) = exp(-(i/2)(m + m_s s(x)) beta), U = C W C) passes the JOINT test at (m, m_s) = (0.2, 0.6): gaps at 0 and pi, isotropic Dirac limit, exponential clustering",
     class_ok and corner_free and zeros_elsewhere == 2 and g0 > 0.05 and gpi > 0.05 and iso2 < iso1 and iso2 < 0.05 and slope < -0.2 and min(np.min(np.abs(evh)), 1) > 1e-3,
     f"gaps {g0:.3f} (0), {gpi:.3f} (pi); m_eff = {om0:.4f}; anisotropy of omega^2 - k^2 relative to m_eff^2: {iso1:.1e} -> {iso2:.1e}; |<psi^dag psi>| at r = 1,3,5: {np.array2string(cvals, precision=2)}, rate {-slope:.2f}")

step("M5.(statement) R1 origin condition: an odd zero mode commuting with H on the UNCHANGED seam space would give an odd ground state gamma Omega, contradicting the unique even vacuum; bound zero modes need a defect sector or changed boundary conditions; they add no c = 1/2 branch (refines K9)",
     True, "vacuum uniqueness of the reconstruction contract")

# M6 ---------------------------------------------------------------------------------------------------
def fermion_ops(n):
    def op(i, mm):
        out = np.eye(1, dtype=complex)
        for j in range(n):
            out = np.kron(out, Zp if j < i else (mm if j == i else I2))
        return out
    return [op(i, lower) for i in range(n)]
psi = fermion_ops(4)
Nop = sum(p.conj().T @ p for p in psi)
PF = np.diag(np.where(np.round(np.diag(Nop).real).astype(int) % 2 == 0, 1.0, -1.0))
bil_beta = sum(psi[i].conj().T @ psi[j] * beta[i, j] for i in range(4) for j in range(4))
bil_even = np.allclose(bil_beta @ PF, PF @ bil_beta)
rot = expm(-1j * 2 * np.pi * (-(1j / 2) * (e[1] @ e[2])))                 # 2 pi rotation on the coin spinor
center = np.allclose(rot, -I4)
full = np.kron(rot, np.eye(16))
direct_sum = np.allclose(full, -np.eye(64))
step("M6.R2 has no parity problem: the second-quantised bilinear psi^dag beta psi commutes with (-1)^F; the coin spin lift (2 pi rotation = -1) acts as -I on all 16 internal components of psi in C^4 (x) C^16: every species carries spin 1/2 (the (-1)^16 of K7 concerns only species-diagonal Fock lifts)",
     bil_even and center and direct_sum, "single-particle coin reading")
step("M7.(statement) AX.MARKMAJORANA -> AX.MARKCLIFFORD (R2): the single-particle coin space of the propagating field is an irreducible internal-singlet Cl(4)-module with generators labelled by the four marks; R1 (physical zero modes) is an alternative needing six zero modes, a degenerate reference and a defect sector",
     True, "AX.SHEETMASS defined explicitly in M4")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "joint": {"m": mJ, "m_s": msJ, "gap0": g0, "gappi": gpi, "m_eff": om0, "clustering": cvals.tolist()}}
(HERE / "mark_parity_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
