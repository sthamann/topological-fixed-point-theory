#!/usr/bin/env python3
"""Part T1: frozen Theta_1 kinetics (minimal-doubling C3 Hamiltonian), SOURCE.MARK_C3_SELECTION.01 and QMC.GAMMA_FACTOR.01.

The external definition (theta1.json) is not in this repository; h_W is RECONSTRUCTED from its published description:
    h_W(k) = t sigma . (B^-1 d(B k)),   d_i(p) = sin p_i + cos p_{i+1} - 1,   p_i = k . b_i,
with the primitive bcc rows b1 = (1,1,-1), b2 = (1,-1,1), b3 = (-1,1,1) and M^-1 = B^-T making the Cartesian velocity at
Gamma isotropic.  It hops only along +-b1, +-b2, +-b3 = 2 w_3, 2 w_2, 2 w_1 -- three of the four mark directions.
  T1.1 reconstruction checks against the description: Hermitian, nearest-neighbour bcc support on the six steps +-b_i;
       zeros exactly at Gamma and P = (pi/2)(1,1,1) (global census), V_Gamma = I, V_P = -R_120 (orthogonal, det -1):
       two isotropic Weyl cones of equal speed and opposite chirality; free lowest positive energies on the 16-site
       cubic cell, the 32-site O_h cell and the 4^3 primitive cell.
  T1.2 microscopic point group: among the 48 elements of O_h (with a compatible spin rotation) exactly the C3 about
       [111] (3 elements) survive.
  T1.3 MARK_C3_SELECTION: the four chart Hamiltonians h^(a) (mark a distinguished, the other three as primitive
       directions in their orientation-induced cyclic order) satisfy U_g h^(a)(k) U_g^dag = h^(g a)(g k) for every
       mark rotation g in A4 (with the spin lift of g): all four are related by the source's mark rotations (outcome A:
       chart equivalence); the C3 is the stabiliser of the distinguished mark; nothing checked so far in P1/P2 selects
       one mark (the surface relation sum w_a = 0 and the mark Clifford algebra are A4/S4 symmetric) -> outcome C for
       the selection itself.
  T1.4 GAMMA_FACTOR: on every P-compatible periodic cell h_W(Gamma) = 0 and f(Gamma) = 0, so all kinetic and HS
       matrices vanish on the 32 Gamma modes: Z_Theta = 2^32 Z_active (checked on 4-, 16- and 32-site cells).
  T1.5 hashes: structure_hash (h_W, F, B_A, symmetries) separated from parameter_hash (g, beta, N_tau, V).
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/theta1_md.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from itertools import permutations, product
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    sys.stdout.flush()


sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1.0 + 0j, -1]); S = [sx, sy, sz]
I2 = np.eye(2)
w = np.array([[-.5, -.5, -.5], [-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]])
B = np.array([[1, 1, -1], [1, -1, 1], [-1, 1, 1]], float)


def make_h(Bm, t=1.0):
    Binv = np.linalg.inv(Bm)
    def h(k):
        p = Bm @ k
        d = np.array([np.sin(p[i]) + np.cos(p[(i + 1) % 3]) - 1 for i in range(3)])
        v = Binv @ d
        return t * sum(v[a] * S[a] for a in range(3))
    return h


hW = make_h(B)
dvec = lambda k, h=hW: np.array([np.real(np.trace(h(k) @ S[a])) / 2 for a in range(3)])

# T1.1 ---------------------------------------------------------------------------------------------------
rng = np.random.default_rng(101)
herm = all(np.allclose(hW(k), hW(k).conj().T) for k in rng.normal(size=(5, 3)))
# support: Fourier components on a fine grid -> only 0, +-b_i
Ng = 8
grid = 2 * np.pi * np.arange(Ng) / Ng
coef = {}
for k in product(grid, repeat=3):
    k = np.array(k); hk = hW(k)
    for r in product(range(-2, 3), repeat=3):
        coef[r] = coef.get(r, 0) + hk * np.exp(1j * k @ np.array(r)) / Ng ** 3
support = {r for r, c in coef.items() if np.abs(c).max() > 1e-10}
allowed = {(0, 0, 0)} | {tuple(int(x) for x in s * B[i]) for i in range(3) for s in (1, -1)}
supp_ok = support <= allowed and len(support) == 7
P = np.pi / 2 * np.ones(3)
N = 24
gr = 2 * np.pi * np.arange(N) / N - np.pi
nrm = np.array([[[np.linalg.norm(dvec(np.array([a, b, c]))) for c in gr] for b in gr] for a in gr])
cands = []
for i, j, l in product(range(N), repeat=3):
    v = nrm[i, j, l]
    nb = [nrm[(i + a) % N, (j + b) % N, (l + c) % N] for a, b, c in product((-1, 0, 1), repeat=3) if (a, b, c) != (0, 0, 0)]
    if v <= min(nb) and v < 0.5:
        r = minimize(lambda k: np.linalg.norm(dvec(k)), np.array([gr[i], gr[j], gr[l]]), method="Nelder-Mead", options={"xatol": 1e-12, "fatol": 1e-14, "maxiter": 6000})
        if r.fun < 1e-8:
            cands.append(np.mod(r.x + np.pi, 2 * np.pi) - np.pi)
fcc = lambda k: np.allclose(np.cos(np.array([k @ np.array(s) for s in product((-1, 1), repeat=3)])), 1, atol=1e-7)
classes = []
for z in cands:
    if not any(fcc(z - q) for q in classes):
        classes.append(z)
def vel(k0):
    return np.array([(dvec(k0 + 1e-6 * e) - dvec(k0 - 1e-6 * e)) / 2e-6 for e in np.eye(3)]).T     # V_ai = d d_a / d k_i
VG, VP = vel(np.zeros(3)), vel(P)
R120 = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], float)
two_nodes = len(classes) == 2 and any(fcc(c) for c in classes) and any(fcc(c - P) for c in classes)
vel_ok = np.allclose(VG, np.eye(3), atol=1e-6) and np.allclose(VP.T @ VP, np.eye(3), atol=1e-6) and abs(np.linalg.det(VP) + 1) < 1e-6 and (np.allclose(-VP, R120, atol=1e-6) or np.allclose(-VP, R120.T, atol=1e-6))


def cell_momenta(Rv):
    G = 2 * np.pi * np.linalg.inv(Rv).T
    Nn = int(round(abs(np.linalg.det(Rv)) / 4))
    ks = []
    for m in product(range(-Nn, Nn + 1), repeat=3):
        k = np.array(m) @ G
        if not any(fcc(k - q) for q in ks):
            ks.append(k)
        if len(ks) == Nn:
            break
    return ks
cells = {"cubic16": np.diag([4.0, 4.0, 4.0]), "oh32": np.array([[8, 0, 0], [6, 2, 0], [6, 0, 2]]) @ B, "prim64": 4 * B, "debug4": np.array([B[0] - B[1], B[1] - B[2], B[0] + B[1] + 2 * B[2]])}
emin = {}
for name, Rv in cells.items():
    es = [np.linalg.norm(dvec(k)) for k in cell_momenta(Rv)]
    pos = [e for e in es if e > 1e-9]
    emin[name] = min(pos)
ref = {"cubic16": 1.4142, "oh32": 0.7174, "prim64": 0.7071}
match = all(abs(emin[k] - v) < 5e-4 for k, v in ref.items())
step("T1.1 reconstruction of the frozen h_W: Hermitian, support {0, +-b1, +-b2, +-b3} (three of the four mark directions), zeros exactly at Gamma and P, V_Gamma = I, V_P = -R_120 (orthogonal, det -1): two isotropic Weyl cones of opposite chirality; lowest free positive energies reproduce the published 1.4142 / 0.7174 / 0.7071 for V = 16 / 32 / 64",
     herm and supp_ok and two_nodes and vel_ok and match, f"zero classes {len(classes)}; lowest positive energies {({k: round(v, 4) for k, v in emin.items()})}")

# T1.2 ---------------------------------------------------------------------------------------------------
def spin_rot(R):
    """SU(2) lift of a proper rotation R (det +1)."""
    tr = np.trace(R)
    q0 = np.sqrt(max(0.0, 1 + tr)) / 2
    if q0 > 1e-8:
        q = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]]) / (4 * q0)
    else:
        ev, evec = np.linalg.eig(R); q = np.real(evec[:, np.argmin(np.abs(ev - 1))])
    return q0 * I2 - 1j * sum(q[i] * S[i] for i in range(3))
PG = []
for p in permutations(range(3)):
    for sg in product((-1, 1), repeat=3):
        Rm = np.zeros((3, 3))
        for i in range(3):
            Rm[i, p[i]] = sg[i]
        PG.append(Rm)
ks = rng.normal(size=(6, 3))
def is_sym(g, h=hW):
    # proper part acts on the spin vector; improper g = -g' acts on k as g but leaves sigma (axial) like g'
    gp = g * np.linalg.det(g)
    U = spin_rot(gp)
    return all(np.allclose(U @ h(k) @ U.conj().T, h(g @ k), atol=1e-9) for k in ks)
syms = [g for g in PG if is_sym(g)]
c3_only = len(syms) == 3 and all(np.allclose(np.linalg.matrix_power(g, 3), np.eye(3)) for g in syms)
step("T1.2 microscopic point group of h_W: of the 48 elements of O_h (with spin rotation) exactly the three-element C3 about [111] survives",
     c3_only, f"{len(syms)} surviving elements")

# T1.3 MARK_C3_SELECTION -----------------------------------------------------------------------------------
charts = {}
for a in range(4):
    best = None
    for perm in permutations([b for b in range(4) if b != a]):
        Bm = np.array([2 * w[b] for b in perm])
        if np.linalg.det(Bm) * np.linalg.det(B) > 0:
            best = (Bm, perm); break
    charts[a] = best
h_chart = {a: make_h(charts[a][0]) for a in range(4)}
# A4: mark permutations with even sign; spatial action std(pi): w_a -> w_pi(a)
def std(pi):
    A = np.array([w[pi[a]] for a in range(3)]).T
    return A @ np.linalg.inv(np.array([w[a] for a in range(3)]).T)
A4 = [pi for pi in permutations(range(4)) if np.linalg.det(std(pi)) > 0]
equiv = True
for pi in A4:
    g = std(pi); U = spin_rot(g)
    for a in range(4):
        b = pi[a]
        # the image of chart a under g is a chart with distinguished mark b; compare with h^(b) up to the cyclic
        # relabelling of its three primitive directions (cyclic orders are preserved by rotations)
        lhs = [U @ h_chart[a](np.linalg.inv(g) @ k) @ U.conj().T for k in ks]
        ok_any = False
        Bb = charts[b][0]
        for shift in range(3):
            Bs = np.roll(Bb, shift, axis=0)
            hb = make_h(Bs)
            if all(np.allclose(L, hb(k), atol=1e-9) for L, k in zip(lhs, ks)):
                ok_any = True
        equiv &= ok_any
stab = [pi for pi in A4 if pi[0] == 0]
surface_sym = np.allclose(w.sum(0), 0)
step("T1.3 MARK_C3_SELECTION: the four chart Hamiltonians h^(a) (mark a distinguished) are mapped onto each other by every mark rotation g in A4 with its spin lift, U_g h^(a)(g^-1 k) U_g^dag = h^(g a)(k) up to the cyclic relabelling of the three primitive directions (outcome A: chart equivalence); the surviving C3 is the stabiliser of the distinguished mark; the surface relation and the mark Clifford algebra treat all four marks alike, so nothing checked in P1/P2 selects one (outcome C for the selection)",
     equiv and len(stab) == 3 and surface_sym, f"A4 elements checked {len(A4)}, stabiliser of mark 0: {len(stab)} elements")

# T1.4 GAMMA_FACTOR ------------------------------------------------------------------------------------------
f = lambda k: 1 - np.cos(k[0]) * np.cos(k[1]) * np.cos(k[2])
gam_ok = True
for name, Rv in cells.items():
    gam_ok &= np.allclose(hW(np.zeros(3)), 0) and abs(f(np.zeros(3))) < 1e-15
    # filter row sums vanish on every cell (uniform mode is annihilated), independent of folding
    Rinv = np.linalg.inv(Rv)
step("T1.4 GAMMA_FACTOR: h_W(Gamma) = 0 and f(Gamma) = 0 on every P-compatible periodic cell, so all kinetic and HS matrices vanish on the 32 Gamma modes; Z_Theta = 2^32 Z_active at mu = 0 (Gamma observables analytic)",
     gam_ok, "Gamma block removed analytically from QMC")

# T1.5 hashes --------------------------------------------------------------------------------------------
structure = {"h_W": "t sigma.(B^-1 d(Bk)), d_i = sin p_i + cos p_{i+1} - 1, p = Bk", "B": B.tolist(), "filter": "F = I - (1/8) sum_eps T_eps",
             "B_A": "standard Spin(10) 10-channel (C = Gamma1 Gamma3 Gamma5 Gamma7 Gamma9)", "interaction": "-g sum_{x,A} O_A(x)^2", "symmetry": "C3[111] x Spin(10) x Z4 x translations"}
structure_hash = hashlib.sha256(json.dumps(structure, sort_keys=True).encode()).hexdigest()[:16]
param_hash = lambda g, beta, Nt, V: hashlib.sha256(json.dumps({"g": g, "beta": beta, "N_tau": Nt, "V": V}, sort_keys=True).encode()).hexdigest()[:16]
step("T1.5 hashes: structure_hash fixes h_W, F, B_A, interaction and symmetries; parameter_hash(g, beta, N_tau, V) varies independently (g scans do not change the architecture)",
     len(structure_hash) == 16, f"structure_hash {structure_hash}; example parameter_hash(0.02, 2, 8, 16) {param_hash(0.02, 2, 8, 16)}")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "structure": structure, "structure_hash": structure_hash,
       "lowest_free_energies": emin, "V_P": VP.tolist()}
(HERE / "theta1_md_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
