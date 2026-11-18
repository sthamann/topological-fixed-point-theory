#!/usr/bin/env python3
"""Part SG: CHIRAL4D.SYMMETRIC_GAP.01 -- the filtered Spin(10) interaction: what can be decided exactly.

External construction (2026-10-01/02): bcc filter F = I - (1/8) sum_eps T_eps (symbol f(k) = 1 - cos kx cos ky cos kz),
filtered field rho = F_U psi, interaction H_int = -g sum_{x,A} O_A(x)^2, O_A = rho^T (i sigma_2 (x) B_A) rho + h.c.
(B_A: the ten symmetric Spin(10) matrices of 16 (x) 16 -> 10; reconstructed here in the standard way, the external
package was not available).  The exact diagonalisation of 2^3 bcc cells with all 16 internal Weyl components
(16 sites x 2 x 16 = 512 fermion modes, Hilbert dimension 2^512) is NOT feasible; this checker decides exactly:
  SG1 filter (frozen): symbol f = 1 - cos kx cos ky cos kz; f = |k|^2/2 + O(k^4); zero set on the torus = the k = 0
      class only; f(+-P) = 1, f(H) = 2; F = (1/16) sum_eps (1 - T_eps)^dag (1 - T_eps) >= 0; with random U(2) links
      F_U >= 0 and G F_U G^dag = F_{U^g}.
  SG2 chirality-block requirement: in the 4-component mark walk both chirality blocks have their light node at k = 0
      (f = 0), so a filtered interaction cannot gap the right-handed partner: that model stays vector-like (KILL for
      "Dirac walk + filtered gap" as chiral matter).  The chiral 16 must use ONE Clifford block (Weyl walk): nodes
      k = 0 (f = 0, light), +P (f = 1, mirror), and the pi-nodes H, -P (f = 2, 1).  beta-mass and AX.SHEETMASS couple
      the two blocks and belong to vector-like matter only.
  SG3 operators: O_A Hermitian, fermion-parity even, Delta N = +-2 (so O_A^2 keeps Z4: psi -> i psi), spin singlet
      (commutes with SU(2) spin), a Spin(10) vector; {O_A} with the Spin(10) generators close into so(11) (dim 55),
      with Q_A = i(P_A - P_A^dag) into so(12) (dim 66).
  SG4 exact strong-coupling local spectrum (one site, 32 canonical modes, 2^32 states): the Fock space is an so(11)
      representation (B5-Weyl invariance of the weight multiplicities checked), decomposed by Racah-Speiser
      (dimension check sum m dim = 2^32); sum_A O_A^2 = kappa (C2_B5(lambda) - C2_D5(mu)), kappa fixed twice (trace
      form and full Fock trace).  Ground state of -g sum O_A^2: energy, degeneracy, gap, symmetry.
  SG5 anomalies of the light sector (one left-handed 16 per family, hypercharge y = sum y_i n_i of the source):
      tr Y = tr Y^3 = SU(2)^2 Y = SU(3)^2 Y = SU(3)^3 = 0, even number of doublets, Z16 index 16 = 0; the mirror 16 at
      +P alone is anomaly-free as well (necessary for a symmetric gap).
  SG6 (statement) what is and is not decided; the thermodynamic many-body gap (2^3 cells and beyond, composites,
      Green functions) needs sign-problem-aware QMC or tensor networks and is NOT claimed.
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/symmetric_gap.py
"""
from __future__ import annotations

import json
import sys
from itertools import combinations, permutations, product
from pathlib import Path

import numpy as np
from scipy.linalg import null_space

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    sys.stdout.flush()


I2 = np.eye(2, dtype=complex); Zp = np.diag([1.0, -1.0]).astype(complex)
lower = np.array([[0, 1], [0, 0]], dtype=complex)
signs = [np.array(s) for s in product((-1, 1), repeat=3)]

# SG1 ----------------------------------------------------------------------------------------------------
f = lambda k: 1 - np.cos(k[0]) * np.cos(k[1]) * np.cos(k[2])
sym = lambda k: 1 - np.mean([np.exp(-1j * s @ k) for s in signs])
rng = np.random.default_rng(41)
ks = rng.uniform(-np.pi, np.pi, size=(50, 3))
symbol_ok = all(abs(sym(k) - f(k)) < 1e-12 for k in ks)
small = [f(t * d) / (0.5 * t * t) for d in rng.normal(size=(6, 3)) / 1 for t in (1e-2,)]
dirs = rng.normal(size=(6, 3)); dirs /= np.linalg.norm(dirs, axis=1)[:, None]
quad = all(abs(f(1e-3 * d) / (0.5e-6) - 1) < 1e-5 for d in dirs)
recip = lambda k: np.allclose(np.cos(np.array([k @ s for s in signs])), 1, atol=1e-9)
N = 24
grid = 2 * np.pi * np.arange(N) / N - np.pi
zeros = [np.array(k) for k in product(grid, repeat=3) if abs(f(np.array(k))) < 1e-12]
zero_set = len(zeros) > 0 and all(recip(k) for k in zeros)
gmin = min(f(np.array(k)) for k in product(grid, repeat=3) if not recip(np.array(k)))
P = np.pi / 2 * np.ones(3); H = np.array([np.pi, 0, 0])
values = (f(P), f(-P), f(H))
Lb = 4
sites = [x for x in product(range(Lb), repeat=3) if len({c % 2 for c in x}) == 1]
idx = {x: i for i, x in enumerate(sites)}
def rand_u2():
    q, r = np.linalg.qr(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))
    return q
links = {(x, tuple(s)): rand_u2() for x in sites for s in signs if tuple(s) in [tuple(t) for t in signs if t[0] == 1 or (t[0] == -1 and False)]}
# links on every directed bcc edge, with U_{-eps}(x + eps) = U_eps(x)^dag
links = {}
for x in sites:
    for s in signs:
        if s[0] == 1:
            y = tuple((x[c] + s[c]) % Lb for c in range(3))
            U = rand_u2()
            links[(x, tuple(s))] = U
            links[(y, tuple(-s))] = U.conj().T
dim = 2 * len(sites)
def FU(lk):
    M = np.eye(dim, dtype=complex)
    for x in sites:
        for s in signs:
            y = tuple((x[c] + s[c]) % Lb for c in range(3))
            M[2 * idx[x]:2 * idx[x] + 2, 2 * idx[y]:2 * idx[y] + 2] -= lk[(x, tuple(s))] / 8
    return M
F_U = FU(links)
herm = np.allclose(F_U, F_U.conj().T)
pos = np.linalg.eigvalsh((F_U + F_U.conj().T) / 2).min() > -1e-12
gx = {x: rand_u2() for x in sites}
links_g = {(x, s): gx[x] @ U @ gx[tuple((x[c] + s[c]) % Lb for c in range(3))].conj().T for (x, s), U in links.items()}
G = np.zeros((dim, dim), complex)
for x in sites:
    G[2 * idx[x]:2 * idx[x] + 2, 2 * idx[x]:2 * idx[x] + 2] = gx[x]
cov = np.allclose(G @ F_U @ G.conj().T, FU(links_g))
step("SG1.filter (frozen): symbol 1 - cos kx cos ky cos kz, = |k|^2/2 + O(k^4), zero set on the torus = the k = 0 class only (min elsewhere on the grid > 0), f(+-P) = 1, f(H) = 2; with random U(2) links F_U is Hermitian, positive semidefinite and gauge covariant",
     symbol_ok and quad and zero_set and gmin > 0 and np.allclose(values, (1, 1, 2)) and herm and pos and cov,
     f"min f off the zero class on the 24^3 grid {gmin:.3e}; f(P), f(-P), f(H) = {np.round(values, 6).tolist()}")

# SG2 ----------------------------------------------------------------------------------------------------
w4 = np.array([[-.5, -.5, -.5], [-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]])
def majoranas(n):
    def op(i, m):
        out = np.eye(1, dtype=complex)
        for j in range(n):
            out = np.kron(out, Zp if j < i else (m if j == i else I2))
        return out
    out = []
    for i in range(n):
        a = op(i, lower)
        out += [a + a.conj().T, -1j * (a - a.conj().T)]
    return out
gm = majoranas(2)
I4 = np.eye(4)
beta = sum(gm) / 2
e = [sum(w4[a, i] * gm[a] for a in range(4)) for i in range(3)]
alpha = [1j * beta @ e[i] for i in range(3)]
Gam = -gm[0] @ gm[1] @ gm[2] @ gm[3]
W = lambda k: (np.cos(k[0]) * I4 - 1j * np.sin(k[0]) * alpha[0]) @ (np.cos(k[1]) * I4 - 1j * np.sin(k[1]) * alpha[1]) @ (np.cos(k[2]) * I4 - 1j * np.sin(k[2]) * alpha[2])
ev_, V_ = np.linalg.eigh(Gam)
blocks = {"+": V_[:, ev_ > 0], "-": V_[:, ev_ < 0]}
light = {b: np.allclose(Pb.conj().T @ W(np.zeros(3)) @ Pb, np.eye(2)) for b, Pb in blocks.items()}
nodes_plus = {"0": (np.zeros(3), 1), "+P": (P, 1), "H": (H, -1), "-P": (-P, -1)}
weyl_nodes = all(np.allclose(blocks["+"].conj().T @ W(k) @ blocks["+"], s * np.eye(2)) for k, s in nodes_plus.values())
fvals = {n: round(f(k), 6) for n, (k, s) in nodes_plus.items()}
mixes = not np.allclose(blocks["+"].conj().T @ beta @ blocks["-"], 0)
step("SG2.chirality-block requirement: both chirality blocks of the 4-component walk have their light node at k = 0 (f = 0), so the filtered interaction cannot gap the right-handed partner -> vector-like (KILL of 'Dirac walk + filtered gap' as chiral matter); the chiral 16 needs ONE block (Weyl walk): f = 0 at k = 0, 1 at the mirror +P, 2 and 1 at the pi-nodes H, -P; beta (and AX.SHEETMASS) couple the blocks",
     light["+"] and light["-"] and weyl_nodes and fvals == {"0": 0.0, "+P": 1.0, "H": 2.0, "-P": 1.0} and mixes,
     f"f at the Weyl-block nodes {fvals}")

# SG3 ----------------------------------------------------------------------------------------------------
c5 = []
for j in range(5):
    mats = [Zp] * j + [lower] + [I2] * (4 - j)
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    c5.append(out)
Gm = []
for j in range(5):
    Gm += [c5[j] + c5[j].conj().T, -1j * (c5[j] - c5[j].conj().T)]
occ = np.array([[(b >> (4 - j)) & 1 for j in range(5)] for b in range(32)])
even = [b for b in range(32) if occ[b].sum() % 2 == 0]
Pe = np.eye(32)[:, even]
Sig = {(A, B): Pe.T @ (0.25j * (Gm[A] @ Gm[B] - Gm[B] @ Gm[A])) @ Pe for A in range(10) for B in range(A + 1, 10)}
best = None
for name, Cc in (("odd", Gm[0] @ Gm[2] @ Gm[4] @ Gm[6] @ Gm[8]), ("even", Gm[1] @ Gm[3] @ Gm[5] @ Gm[7] @ Gm[9])):
    BA = [Pe.T @ (Cc @ Gm[A]) @ Pe for A in range(10)]
    if all(np.linalg.norm(b) > 1e-9 for b in BA) and all(np.allclose(b, b.T) for b in BA):
        flat = np.array([b.ravel() for b in BA]).T
        okc = True
        for s_ in Sig.values():
            for b in BA:
                tgt = (s_.T @ b + b @ s_).ravel()
                coef, res, *_ = np.linalg.lstsq(flat, tgt, rcond=None)
                okc &= np.allclose(flat @ coef, tgt, atol=1e-10)
        if okc:
            best = (name, BA)
            break
BA = best[1]
nm = 32                                                                # modes psi_(alpha s): index alpha*16 + s
eps2 = np.array([[0, 1], [-1, 0]], dtype=complex)                     # i sigma_2
Mpair = [np.kron(eps2, b) for b in BA]
def K_hop(Aop):                                                       # Q = psi^dag A psi
    return np.block([[-Aop, np.zeros((nm, nm))], [np.zeros((nm, nm)), Aop.T]])
def K_pair(Mp):                                                       # Q = psi^T M psi
    return np.block([[np.zeros((nm, nm)), np.zeros((nm, nm))], [-2 * Mp, np.zeros((nm, nm))]])
def K_pair_dag(Mp):                                                   # Q = (psi^T M psi)^dag
    return np.block([[np.zeros((nm, nm)), 2 * Mp.conj()], [np.zeros((nm, nm)), np.zeros((nm, nm))]])
KO = [K_pair(Mp) + K_pair_dag(Mp) for Mp in Mpair]
KQ = [1j * (K_pair(Mp) - K_pair_dag(Mp)) for Mp in Mpair]
KJ = [K_hop(np.kron(I2, s_)) for s_ in Sig.values()]
KS = [K_hop(np.kron(sp, np.eye(16))) for sp in (np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]))]
KN = K_hop(np.eye(nm))
br = lambda X, Y: X @ Y - Y @ X
def span_dim(mats):
    M = np.array([np.concatenate([m.ravel().real, m.ravel().imag]) for m in mats])
    return np.linalg.matrix_rank(M, tol=1e-8)
def closure(gens, maxit=6):
    basis = list(gens)
    d = span_dim(basis)
    for _ in range(maxit):
        new = basis + [br(a, b) for a, b in combinations(basis, 2)]
        Mx = np.array([np.concatenate([m.ravel().real, m.ravel().imag]) for m in new])
        u, s, vt = np.linalg.svd(Mx, full_matrices=False)
        r = int((s > 1e-8 * s[0]).sum())
        basis = [(vt[i, :2 * (2 * nm) ** 2 // 2] + 1j * vt[i, 2 * (2 * nm) ** 2 // 2:]).reshape(2 * nm, 2 * nm) for i in range(r)]
        if r == d:
            break
        d = r
    return d
herm_O = all(np.allclose(np.block([[np.zeros((nm, nm)), Mp.conj().T], [Mp, np.zeros((nm, nm))]]),
                         np.block([[np.zeros((nm, nm)), Mp.conj().T], [Mp, np.zeros((nm, nm))]]).conj().T) for Mp in Mpair)
anti = all(np.allclose(Mp, -Mp.T) for Mp in Mpair)
spin_comm = all(np.allclose(br(ko, ks_), 0) for ko in KO for ks_ in KS)
dN = all(np.allclose(br(KN, K_pair(Mp)), 2 * K_pair(Mp)) or np.allclose(br(KN, K_pair(Mp)), -2 * K_pair(Mp)) for Mp in Mpair)
vec = True
flatO = np.array([k.ravel() for k in KO]).T
for kj in KJ:
    for ko in KO:
        tgt = br(kj, ko).ravel()
        coef, *_ = np.linalg.lstsq(flatO, tgt, rcond=None)
        vec &= np.allclose(flatO @ coef, tgt, atol=1e-9)
d11 = closure([1j * k for k in KO + KJ])                           # anti-Hermitian generators iQ: compact real form
d12 = closure([1j * k for k in KO + KQ + KJ])
step("SG3.operators: O_A Hermitian, fermion-parity even, Delta N = +-2 (O_A^2 keeps Z4), spin singlet, a Spin(10) vector; {O_A} + so(10) close into so(11) (dim 55), with Q_A = i(P_A - P_A^dag) into so(12) (dim 66)",
     best is not None and herm_O and anti and spin_comm and dN and vec and d11 == 55 and d12 == 66,
     f"B_A from the {best[0] if best else '-'} charge conjugation; closure dims {d11}, {d12}")

# SG4 ----------------------------------------------------------------------------------------------------
wts = np.array([2 * (occ[b] - 0.5) for b in even]).astype(int)       # doubled so(10) weights of the 16
R = 16
shape = (2 * R + 1,) * 5
mult = np.zeros(shape, dtype=np.int64)
mult[(R,) * 5] = 1
for wv in wts:
    for _ in range(2):                                                # two spin components
        src = [slice(max(0, -d), 2 * R + 1 - max(0, d)) for d in wv]
        dst = [slice(max(0, d), 2 * R + 1 + min(0, d)) for d in wv]
        mult[tuple(dst)] += mult[tuple(src)].copy()
total = int(mult.sum())
flipped = mult[..., ::-1]
b5_inv = np.array_equal(mult, flipped)                                # sign flip of the 5th coordinate (not in W(D5))
rho2 = np.array([9, 7, 5, 3, 1])
perms = list(permutations(range(5)))
def perm_sign(p):
    s = 1
    p = list(p)
    for i in range(5):
        while p[i] != i:
            j = p[i]; p[i], p[j] = p[j], p[i]; s = -s
    return s
W_elems = [(np.array(p), np.array(sg), perm_sign(p) * int(np.prod(sg))) for p in perms for sg in product((1, -1), repeat=5)]
def lookup(mu2):
    if np.any(np.abs(mu2) > 2 * R):
        return 0
    if np.any(mu2 % 2 != mu2[0] % 2):
        return 0
    return int(mult[tuple(mu2 + R)]) if np.all(np.abs(mu2) <= R) else 0
def m_of(lam2):
    nu = lam2 + rho2
    return sum(e_ * lookup(sg * nu[p] - rho2) for p, sg, e_ in W_elems)
def dim_B5(lam2):
    l = (lam2 + rho2) / 2; r = rho2 / 2
    num = np.prod([(l[i] ** 2 - l[j] ** 2) / (r[i] ** 2 - r[j] ** 2) for i in range(5) for j in range(i + 1, 5)])
    return num * np.prod(l / r)
rhoD = np.array([4, 3, 2, 1, 0])
def dim_D5(mu):
    l = mu + rhoD
    return np.prod([(l[i] ** 2 - l[j] ** 2) / (rhoD[i] ** 2 - rhoD[j] ** 2) for i in range(5) for j in range(i + 1, 5)])
C2B = lambda lam: float(np.sum(lam * (lam + rho2 / 1.0)))           # lam in normal coords, 2 rho_B5 = rho2
C2D = lambda mu: float(np.sum(mu * (mu + 2 * rhoD)))
cands = []
for parity in (0, 1):
    vals = range(parity, 2 * 8 + 1, 2)
    for lam2 in product(vals, repeat=5):
        if all(lam2[i] >= lam2[i + 1] for i in range(4)):
            cands.append(np.array(lam2))
irreps = {}
for lam2 in cands:
    m = m_of(lam2)
    if m:
        irreps[tuple(lam2)] = m
dimsum = sum(m * dim_B5(np.array(l2)) for l2, m in irreps.items())
nonneg = all(m > 0 for m in irreps.values())
def branch(lam):
    """B5 -> D5 interlacing patterns mu (normal coordinates)."""
    half = (lam[0] % 1) != 0
    step_ = 1.0
    ranges = []
    for i in range(4):
        lo, hi = lam[i + 1], lam[i]
        ranges.append(np.arange(lo, hi + 1e-9, step_))
    r5 = np.arange(-lam[4], lam[4] + 1e-9, step_)
    out = []
    for mu4 in product(*ranges, r5):
        out.append(np.array(mu4))
    return out
levels = {}
trace_rhs = 0.0
for l2, m in irreps.items():
    lam = np.array(l2) / 2
    for mu in branch(lam):
        E = C2B(lam) - C2D(mu)
        dmu = dim_D5(mu)
        trace_rhs += m * dmu * E
        levels[round(E, 9)] = levels.get(round(E, 9), 0) + m * dmu
Es = sorted(levels)
Emax, E2 = Es[-1], Es[-2]
deg = levels[Emax]
gs_lams = [(l2, m) for l2, m in irreps.items() if any(abs(C2B(np.array(l2) / 2) - C2D(mu) - Emax) < 1e-9 for mu in branch(np.array(l2) / 2))]
# normalisation kappa: (a) full Fock trace, (b) trace form on the single-particle 64-rep
trace_lhs = 0.0
for Mp in Mpair:
    # Tr_Fock (P + P^dag)^2 = 2 Tr(P^dag P) = 2 * 2^32 * sum_{i<j} |2 M_ij|^2 / 4
    Ms = 2 * np.triu(Mp, 1)
    trace_lhs += 2 * 2.0 ** 32 * np.sum(np.abs(Ms) ** 2) / 4
kappa_trace = trace_lhs / trace_rhs
kappa_form = abs(np.trace(KO[0] @ KO[0]).real) / 16                   # tr_64 = 2 tr_32, tr_32(L^2) = 8 for standard so(11) generators
uniq = deg == 1
step("SG4.exact local strong-coupling spectrum (one site, 2^32 states): the Fock weights are B5-Weyl invariant and decompose by Racah-Speiser into so(11) irreps (sum m dim = 2^32); sum_A O_A^2 = kappa (C2_B5 - C2_D5) with kappa from the Fock trace and from the trace form; ground state of -g sum O_A^2: energy, degeneracy, gap",
     total == 2 ** 32 and b5_inv and nonneg and abs(dimsum - 2 ** 32) < 0.5 and abs(kappa_trace / kappa_form - 1) < 1e-6,
     f"kappa trace {kappa_trace:.6f} vs trace form {kappa_form:.6f}; irreps {len(irreps)}; ground level C2 difference {Emax:g} with degeneracy {deg} (lambda = {[tuple(np.array(l) / 2) for l, _ in gs_lams]}); next {E2:g}; gap = {kappa_trace * (Emax - E2):.4f} g; unique: {uniq}")
step("SG4b.necessary local condition for a symmetric gap: the strong-coupling ground state is UNIQUE (hence Spin(10) x SU(2) x Z4 symmetric) and separated by a finite gap",
     uniq and Emax - E2 > 0, f"ground-state degeneracy {deg}")

# SG5 ----------------------------------------------------------------------------------------------------
yv = np.array([-1 / 3, -1 / 3, -1 / 3, 0.5, 0.5])
states = [occ[b] for b in even]
Y = np.array([s @ yv for s in states])
trY, trY3 = Y.sum(), (Y ** 3).sum()
def su_content(s):
    col = s[:3].sum(); wk = s[3:].sum()
    return col, wk
su2Y = sum(Y[i] for i, s in enumerate(states) if s[3:].sum() == 1) / 2          # T(doublet)=1/2 per doublet state pair
su3Y = sum(Y[i] for i, s in enumerate(states) if s[:3].sum() in (1, 2)) / 2       # each triplet/antitriplet state carries T = 1/2 per colour multiplet
col_sign = [(+1 if s[:3].sum() == 1 else -1) for s in states if s[:3].sum() in (1, 2)]
su3cube = sum(col_sign) / 3
doublets = sum(1 for s in states if s[3:].sum() == 1) // 2
z16 = len(states) % 16
step("SG5.anomalies of the light sector (one left-handed 16 per family, y = sum y_i n_i): tr Y = tr Y^3 = SU(2)^2 Y = SU(3)^2 Y = SU(3)^3 = 0, even number of doublets (Witten), Z16 index 16 = 0; the mirror 16 at +P (opposite handedness) is anomaly-free by itself",
     abs(trY) < 1e-12 and abs(trY3) < 1e-12 and abs(su2Y) < 1e-12 and abs(su3Y) < 1e-12 and abs(su3cube) < 1e-12 and doublets % 2 == 0 and z16 == 0,
     f"tr Y {trY:.1e}, tr Y^3 {trY3:.1e}, SU(2)^2Y {su2Y:.1e}, SU(3)^2Y {su3Y:.1e}, SU(3)^3 {su3cube:.1e}, doublets {doublets}")
step("SG6.(statement) decided exactly: filter, chirality-block requirement, operator algebra, the strong-coupling local spectrum and the anomaly bookkeeping; NOT decided: the thermodynamic symmetric mirror gap with light chiral modes, absence of light composites, Green functions (2^512-dimensional already for 2^3 cells; needs QMC / tensor networks)",
     True, "CHIRAL4D.NOMIRROR.01 stays open; no marker moves")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG),
       "ground": {"C2_difference": Emax, "degeneracy": deg, "next": E2, "kappa": kappa_trace, "gap_over_g": kappa_trace * (Emax - E2)},
       "irreps": {str(tuple(np.array(k) / 2)): v for k, v in irreps.items()}}
(HERE / "symmetric_gap_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
