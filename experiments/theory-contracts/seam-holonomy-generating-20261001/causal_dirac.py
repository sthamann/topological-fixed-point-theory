#!/usr/bin/env python3
"""Part L: TIME.EVENT_DEPTH.01 / SPACETIME.MICROCAUSALITY.01 for the mark Dirac walk (Part K).

The walk U_D = e^{-i m beta/2} W e^{-i m beta/2}, W(k) = prod_i (cos k_i - i sin k_i alpha_i), with alpha_i, beta
the mark Clifford operators of Part K, acts on the bcc lattice (cubic coordinates, steps eps in {+-1}^3).
  L1 microcausality (exact): after t steps (event depth t) the propagator from the origin is supported in
     the cube |x_i| <= t and reaches its corners: the anticommutator {psi(x,t), psi^dag(0,0)} = G_t(x) vanishes
     identically outside the hard cone; no tails (strict Lieb-Robinson cone of a QCA).
  L2 emergent cone inside the hard cone: every group velocity over the Brillouin zone lies in the cube
     |v_i| <= 1; the IR cone (m = 0) is isotropic with speed 1: the round light cone is inscribed in the cube.
  L3 Lorentz is emergent, not microscopic: omega(k)/|k| - 1 is direction dependent and vanishes LINEARLY with
     |k| (the leading violation is O(k l)); there is no exact boost symmetry on the lattice.
  L4 proper time: event depth t is coordinate time; the phase of a massive excitation moving with group
     velocity v obeys omega - k.v = m sqrt(1 - v^2) in the continuum limit (time dilation emerges).
  L5 doubling census (CORRECTS the reading of Part Q/K as a single Weyl node): at quasi-energy 0 each chiral
     block has TWO nodes, k = 0 and a Brillouin-zone corner P = +-(pi/2)(1,1,1) (bcc "P point"), of opposite
     Berry charge: Nielsen-Ninomiya holds, the walk is vector-like in the IR; the pi-nodes (H = (pi,0,0) and
     -+P) cancel likewise.
  L6 the on-site mass coin gaps k = 0 but NOT the doubler pair: U_D(+-P, m) = +-Gamma-type with eigenvalue +1
     for every m (the partner at the same k sits at quasi-energy pi); gapping needs a coupling of P and -P,
     whose difference Q = (pi,pi,pi) is the bcc sublattice sign.
  L7 sheets: total winding parity is a well-defined Z2 grading of H_1 = Z^4/(1,1,1,1) (the relation has even
     sum), flipped by every primitive step: it is the deck index of the double cover branched at the four
     marks and equals the bcc sublattice sign e^{iQ.x}.  The sheet-odd mass coin exp(-i/2 (m + m_s s(x)) beta)
     gaps the walk at quasi-energies 0 AND pi everywhere, and its quasi-free Dirac sea clusters exponentially;
     without the sheet-odd term (m_s = 0) the vacuum does not cluster exponentially.
  L8 (statement) dimension typing (DIMENSION.UPLIFT.FIREWALL.01, F2): hypotheses = 3D bcc lattice + discrete step
     time (+ AX.MARKMAJORANA, + sheet-odd mass); conclusion = free 3+1D kinematics of a vector-like Dirac pair
     (light at k = 0, doubler pair at +-P with sheet-odd mass) with strict causality, emergent Lorentz cone,
     proper time and clustering vacuum.  NOT derived: a chiral IR, interactions, gauge fields, gravity.
     Weinberg-Witten (GRAV.WEINBERG.WITTEN.01) only cross-linked: no microscopic Lorentz symmetry (L3).
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/causal_dirac.py
"""
from __future__ import annotations

import json
import sys
from itertools import product
from pathlib import Path

import numpy as np
from scipy.linalg import eigh as sp_eigh, expm

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


# mark Clifford operators (Part K) -------------------------------------------------------------------
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
signs = [np.array(s) for s in product((-1, 1), repeat=3)]
A = {tuple(s): (I4 + s[0] * alpha[0]) / 2 @ (I4 + s[1] * alpha[1]) / 2 @ (I4 + s[2] * alpha[2]) / 2 for s in signs}


def W(k):
    return (np.cos(k[0]) * I4 - 1j * np.sin(k[0]) * alpha[0]) @ (np.cos(k[1]) * I4 - 1j * np.sin(k[1]) * alpha[1]) @ (np.cos(k[2]) * I4 - 1j * np.sin(k[2]) * alpha[2])


def UD(k, m):
    c = expm(-0.5j * m * beta)
    return c @ W(k) @ c


def quasi(k, m):
    return np.sort(np.angle(np.linalg.eigvals(UD(k, m))))


# L1 exact cone in position space ------------------------------------------------------------------
m, T = 0.3, 6
Lg = 2 * T + 4
c = expm(-0.5j * m * beta)
psi = np.zeros((Lg, Lg, Lg, 4), dtype=complex)
psi[0, 0, 0] = np.array([1, 0.3, -0.2j, 0.5]) / np.linalg.norm([1, 0.3, 0.2, 0.5])
for _ in range(T):
    psi = psi @ c.T
    new = np.zeros_like(psi)
    for s, As in A.items():
        new += np.roll(psi, shift=s, axis=(0, 1, 2)) @ As.T
    psi = new @ c.T
coords = np.array(np.meshgrid(*(np.where(np.arange(Lg) <= Lg // 2, np.arange(Lg), np.arange(Lg) - Lg) for _ in range(3)), indexing="ij"))
cheb = np.max(np.abs(coords), axis=0)
amp = np.linalg.norm(psi, axis=-1)
outside = amp[cheb > T].max()
corner = amp[(np.abs(coords[0]) == T) & (np.abs(coords[1]) == T) & (np.abs(coords[2]) == T)].max()
norm = np.sum(amp ** 2)
step("L1.microcausality: after t = 6 steps the propagator is supported EXACTLY in the cube |x_i| <= t (zero outside, nonzero at the corners), norm preserved: strict causal cone, {psi(x,t), psi^dag(0,0)} = 0 outside",
     outside == 0.0 and corner > 1e-6 and abs(norm - 1) < 1e-12, f"max outside {outside:.1e}, corner amplitude {corner:.2e}, norm {norm:.12f}")

# L2 group velocities --------------------------------------------------------------------------------
def band_velocities(k, m, h=1e-6):
    e0 = quasi(k, m)
    return np.array([(quasi(k + h * np.eye(3)[i], m) - e0) / h for i in range(3)]).T     # bands x 3
rng = np.random.default_rng(11)
vmax = 0.0
for k in rng.uniform(-np.pi, np.pi, size=(400, 3)):
    v = band_velocities(k, 0.0)
    vmax = max(vmax, np.abs(v).max())
dirs = rng.normal(size=(20, 3)); dirs /= np.linalg.norm(dirs, axis=1)[:, None]
vir = [np.linalg.norm(band_velocities(1e-3 * n, 0.0)[-1]) for n in dirs]
step("L2.emergent cone inside the hard cone: all group velocities over the Brillouin zone satisfy |v_i| <= 1 (cube); the IR cone (m = 0) is isotropic with speed 1, inscribed in the cube",
     vmax <= 1 + 1e-4 and np.allclose(vir, 1, atol=1e-3), f"max |v_i| over BZ {vmax:.4f}; IR speeds in [{min(vir):.5f}, {max(vir):.5f}]")

# L3 Lorentz violation scaling -----------------------------------------------------------------------
def lv(kmag):
    vals = [quasi(kmag * n, 0.0)[-1] / kmag - 1 for n in dirs]
    return max(abs(x) for x in vals), np.ptp(vals)
l1, l2 = lv(0.02), lv(0.01)
step("L3.Lorentz is emergent: omega/|k| - 1 is direction dependent and vanishes linearly with |k| (leading violation O(k l)); no exact lattice boosts",
     l1[1] > 1e-4 and abs(l1[0] / l2[0] - 2) < 0.1, f"max |omega/|k| - 1| = {l1[0]:.2e} (|k| = 0.02), {l2[0]:.2e} (0.01); spread over directions {l1[1]:.1e}")

# L4 proper time -------------------------------------------------------------------------------------
def dil(sc):
    mm, out = 0.05 * sc, []
    for n in dirs[:8]:
        for kmag in (0.03 * sc, 0.06 * sc):
            k = kmag * n
            v = band_velocities(k, mm, h=1e-7)[-1]
            om = quasi(k, mm)[-1]
            out.append(abs((om - k @ v) / (mm * np.sqrt(1 - v @ v)) - 1))
    return max(out)
d1, d2 = dil(1.0), dil(0.5)
step("L4.proper time: event depth t is coordinate time; a massive excitation with group velocity v accumulates phase (omega - k.v) t = m sqrt(1 - v^2) t in the continuum limit: emergent time dilation",
     d1 < 0.05 and d1 / d2 > 1.6, f"max relative deviation {d1:.2e} -> {d2:.2e} at half scale")

# L5 doubling census ---------------------------------------------------------------------------------
def berry_block(k0, P, target, r=0.03, n=30):
    th, ph = np.linspace(0, np.pi, n + 1), np.linspace(0, 2 * np.pi, n + 1)
    def vec(t, f):
        k = k0 + r * np.array([np.sin(t) * np.cos(f), np.sin(t) * np.sin(f), np.cos(t)])
        ev, V = np.linalg.eig(P.conj().T @ W(k) @ P)
        j = np.argmax(np.angle(ev * np.exp(-1j * target)))
        return V[:, j] / np.linalg.norm(V[:, j])
    flux = 0.0
    for i in range(n):
        for j in range(n):
            u = [vec(th[i], ph[j]), vec(th[i + 1], ph[j]), vec(th[i + 1], ph[j + 1]), vec(th[i], ph[j + 1])]
            flux += np.angle(np.vdot(u[0], u[1]) * np.vdot(u[1], u[2]) * np.vdot(u[2], u[3]) * np.vdot(u[3], u[0]))
    return flux / (2 * np.pi)
ev_, V_ = np.linalg.eigh(Gam)
blocks = {"+": V_[:, ev_ > 0], "-": V_[:, ev_ < 0]}
recip = lambda k: np.allclose(np.cos(np.array([k @ s for s in signs])), 1, atol=1e-9)     # k in the fcc reciprocal lattice
Nn = 24
grid = 2 * np.pi * np.arange(Nn) / Nn - np.pi
Pc = np.pi / 2 * np.ones(3)
census, ok5 = {}, True
for b, P in blocks.items():
    z, pi_ = [], []
    for k in product(grid, repeat=3):
        k = np.array(k)
        ev = np.angle(np.linalg.eigvals(P.conj().T @ W(k) @ P))
        if np.max(np.abs(ev)) < 1e-9:
            z.append(k)
        if np.max(np.pi - np.abs(ev)) < 1e-9:
            pi_.append(k)
    cls0 = {("G" if recip(k) else "+P" if recip(k - Pc) else "-P" if recip(k + Pc) else "other") for k in z}
    clspi = {("H" if recip(k - np.array([np.pi, 0, 0])) else "+P" if recip(k - Pc) else "-P" if recip(k + Pc) else "other") for k in pi_}
    reps0 = {"G": np.zeros(3), "+P": Pc, "-P": -Pc}
    repspi = {"H": np.array([np.pi, 0, 0]), "+P": Pc, "-P": -Pc}
    ch0 = {c: round(berry_block(reps0[c], P, 0.0), 4) for c in cls0 if c != "other"}
    chpi = {c: round(berry_block(repspi[c], P, np.pi), 4) for c in clspi if c != "other"}
    census[b] = (ch0, chpi)
    ok5 &= "other" not in cls0 | clspi and len(ch0) == 2 and "G" in ch0 and sum(ch0.values()) == 0 and sum(chpi.values()) == 0 and all(abs(abs(x) - 1) < 1e-6 for x in list(ch0.values()) + list(chpi.values()))
step("L5.doubling census (corrects the single-node reading of Parts Q/K): per chiral block the quasi-energy-0 nodes are k = 0 and one BZ corner +-P = +-(pi/2)(1,1,1) with OPPOSITE charge; the pi-nodes cancel likewise: Nielsen-Ninomiya holds, the walk is vector-like in the IR",
     ok5, f"block +: qe 0 {census['+'][0]}, qe pi {census['+'][1]}; block -: qe 0 {census['-'][0]}, qe pi {census['-'][1]}")

ms_list = (0.1, 0.3, 1.0)
gapless = all(np.min(np.abs(quasi(sgn * Pc, mm))) < 1e-12 for mm in ms_list for sgn in (1, -1))
gap_G = all(abs(np.min(np.abs(quasi(np.zeros(3), mm))) - mm) < 1e-12 for mm in ms_list)
Q = np.pi * np.ones(3)
step("L6.the on-site mass coin gaps k = 0 (quasi-energies +-m) but leaves the doubler pair EXACTLY massless: U_D(+-P, m) has eigenvalue 1 for every m; P and -P differ by Q = (pi,pi,pi), not a reciprocal vector",
     gapless and gap_G and not recip(2 * Pc) and np.allclose(2 * Pc, Q),
     f"m = {ms_list}: min |quasi-energy| at +-P = 0 (exact)")

# L7 sheets ------------------------------------------------------------------------------------------
par_ok = all((sum(s) % 2) == 1 for s in signs) and (1 + 1 + 1 + 1) % 2 == 0
sheet_sign = all(np.isclose(np.cos(Q @ np.array(x)), (-1) ** (int(sum(x)) % 2)) for x in [(0, 0, 0), (1, 1, 1), (2, 0, 0), (1, -1, 1), (3, 1, -1), (2, 2, 4)])
def U8(k, mm, mms):
    Wk = np.zeros((8, 8), dtype=complex)
    Wk[:4, :4], Wk[4:, 4:] = W(k), W(k + Q)
    sx = np.array([[0, 1], [1, 0]])
    Cc = expm(-0.5j * (mm * np.kron(np.eye(2), beta) + mms * np.kron(sx, beta)))
    return Cc @ Wk @ Cc
def gaps8(mm, mms, Ng=16):
    g0 = gpi = np.inf
    gr = 2 * np.pi * np.arange(Ng) / Ng - np.pi
    for k in product(gr, repeat=3):
        q = np.angle(np.linalg.eigvals(U8(np.array(k), mm, mms)))
        g0 = min(g0, np.min(np.abs(q))); gpi = min(gpi, np.min(np.pi - np.abs(q)))
    for k in (np.zeros(3), Pc, -Pc, np.array([np.pi, 0, 0])):
        q = np.angle(np.linalg.eigvals(U8(k, mm, mms)))
        g0 = min(g0, np.min(np.abs(q))); gpi = min(gpi, np.min(np.pi - np.abs(q)))
    return g0, gpi
g_stag = gaps8(0.5, 0.5)
g_unif = gaps8(0.5, 0.0)

Lb = 16
sites = [x for x in product(range(Lb), repeat=3) if len({c % 2 for c in x}) == 1]          # bcc = eee + ooo
index = {x: i for i, x in enumerate(sites)}
def U_real(mm, mms):
    n = len(sites)
    Wr = np.zeros((4 * n, 4 * n), dtype=complex)
    for x, i in index.items():
        for s_ in signs:
            y = tuple((x[c] - s_[c]) % Lb for c in range(3))
            Wr[4 * i:4 * i + 4, 4 * index[y]:4 * index[y] + 4] += A[tuple(s_)]
    Cd = np.zeros_like(Wr)
    for x, i in index.items():
        sh = 1 if x[0] % 2 == 0 else -1
        Cd[4 * i:4 * i + 4, 4 * i:4 * i + 4] = expm(-0.5j * (mm + mms * sh) * beta)
    return Cd @ Wr @ Cd
def vacuum_corr(mm, mms):
    U = U_real(mm, mms)
    Hs = (U - U.conj().T) / 2j
    ev, V = sp_eigh(Hs, driver="evr")
    Vn = V[:, ev < 0]
    o = index[(0, 0, 0)]
    vals = []
    for r in range(0, Lb // 2 + 1):
        j = index[(r % Lb, r % Lb, r % Lb)]
        vals.append(np.linalg.norm(Vn[4 * j:4 * j + 4] @ Vn[4 * o:4 * o + 4].conj().T))
    return np.array(vals), np.min(np.abs(ev))
cs_, hmin_s = vacuum_corr(0.5, 0.5)
cu_, hmin_u = vacuum_corr(0.5, 0.0)
rr = np.arange(len(cs_))
odd_r = rr[1:6:2]                                                # opposite sheet, r <= 5 < L/2 (periodic images); same-sheet entries vanish
same_sheet_zero = np.max(cs_[2::2]) < 1e-12
fit = np.polyfit(odd_r, np.log(cs_[odd_r]), 1)
lin_res = np.abs(np.polyval(fit, odd_r) - np.log(cs_[odd_r])).max()
step("L7.sheets: total winding parity is a Z2 grading of H_1 (every primitive step flips it, the relation (1,1,1,1) is even) = deck index of the double cover branched at the marks = bcc sublattice sign e^{iQ.x}; the sheet-odd mass gaps quasi-energies 0 AND pi everywhere and its Dirac sea clusters exponentially; with m_s = 0 the walk stays gapless at +-P",
     par_ok and sheet_sign and min(g_stag) > 0.1 and min(g_unif) < 1e-9 and hmin_s > 0.05 and fit[0] < -0.3 and lin_res < 0.3 and hmin_u < 1e-9 and same_sheet_zero,
     f"gaps (0, pi): sheet-odd {np.round(g_stag, 3).tolist()}, uniform only {np.round(g_unif, 6).tolist()}; |<psi^dag(r,r,r) psi(0)>| at r = 1,3,5: {np.array2string(cs_[odd_r], precision=2)} (same-sheet r even: 0); decay rate {-fit[0]:.3f}/step; uniform-only min |sin omega| {hmin_u:.1e} (no Dirac sea)")
step("L8.(statement) dimension typing: hypotheses = 3D bcc lattice + step time + AX.MARKMAJORANA + sheet-odd mass; conclusion = free 3+1D kinematics of a vector-like Dirac pair (light pair at k = 0, sheet-gapped doubler pair at +-P) with strict causality, emergent Lorentz cone, proper time and clustering vacuum; NOT derived: chiral IR, interactions, gauge fields, gravity; Weinberg-Witten only cross-linked",
     True, "no spin-2, chirality or TOE marker moves")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "vmax": vmax,
       "lorentz_violation": [l1[0], l2[0]], "census": {b: [list(c[0].items()), list(c[1].items())] for b, c in census.items()},
       "gaps_sheet_odd": list(g_stag), "gaps_uniform": list(g_unif), "clustering": cs_.tolist()}
(HERE / "causal_dirac_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
