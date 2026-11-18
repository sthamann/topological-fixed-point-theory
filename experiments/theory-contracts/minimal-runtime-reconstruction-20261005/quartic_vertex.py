#!/usr/bin/env python3
"""SOURCE.QUARTIC.VERTEX.01 -- does the source select the Theta_1 vertex -g sum_A O_A^2?

  QV1 premise (imported from Part H, not re-derived here): the source vacuum is a PURE quasi-free state (Hardy sea, a
      Slater state of fixed parity), so Wick factorisation holds and survives the GSO/parity restriction; its connected
      four-point function vanishes in every channel and the source selects NO channel.  The check is a consistency check
      of the Wick premise (number-conserving and pairing Gaussian states, interacting control nonzero); it also records
      that parity-projected MIXED Gaussian states do not factorise (not the source state).
  QV2 symmetry alone: per site V = 16 (x) 2 (Spin(10) x SU(2)_spin), Z4: psi -> i psi allows Delta N in {0, +-4}.
      Lambda^2 V = (120,3) + (10,1) + (126,1), multiplicity-free (Casimir spectrum, dims 360/10/126): three independent
      number-conserving invariants |P_120,3|^2, |P_10|^2, |P_126|^2.
  QV3 (Lambda^4 V)^G is one-dimensional (zero-weight Casimir kernel); QV3b it lies in the span of the wedge products of
      the (10,1) channel and is nonzero: the unique Delta N = 4 invariant is proportional to sum_A P_A P_A.
  QV4 real dimension of Hermitian, Z4-invariant, G-invariant quartic vertices = 3 + 2 = 5; the global U(1) phase of psi
      (a symmetry of h_W and F_U, assuming no Xi-map or Majorana constraint breaks it) rotates the Delta N = 4 coefficient,
      so 4 physical couplings remain (modulo one-body normal-ordering terms); Theta_1's
      sum_A O_A^2 = sum_A (P_A^2 + h.c.) + sum_A {P_A, P_A^dag} is ONE ray.  Every member keeps the IR protection
      (depends on psi only through rho = F psi), so IR protection does not select either.
  QV5 (statement) SOURCE.QUARTIC.VERTEX.01 = NOT_SELECTED.
    python3 -B experiments/theory-contracts/minimal-runtime-reconstruction-20261005/quartic_vertex.py
"""
from __future__ import annotations

import sys
from itertools import combinations, product
from pathlib import Path

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ur_lib import Log, jordan_wigner  # noqa: E402

log = Log()
rng = np.random.default_rng(3)

# ---------------------------------------------------------------- QV1 quasi-free -> connected 4-point = 0
n = 4
a = jordan_wigner(n)
ad = [x.conj().T for x in a]


def gibbs(H):
    w, v = np.linalg.eigh(H)
    r = v @ np.diag(np.exp(-(w - w.min()))) @ v.conj().T
    return r / np.trace(r)


def two(rho, X, Y): return np.trace(rho @ X @ Y)


def connected(rho):
    worst = 0.0
    ops = a + ad
    for i, j, k, l in product(range(2 * n), repeat=4):
        A_, B_, C_, D_ = ops[i], ops[j], ops[k], ops[l]
        full = np.trace(rho @ A_ @ B_ @ C_ @ D_)
        wick = (two(rho, A_, B_) * two(rho, C_, D_) - two(rho, A_, C_) * two(rho, B_, D_) + two(rho, A_, D_) * two(rho, B_, C_))
        worst = max(worst, abs(full - wick))
    return worst


h = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)); h = h + h.conj().T
D = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)); D = D - D.T
H_nc = sum(h[i, j] * ad[i] @ a[j] for i in range(n) for j in range(n))
H_bdg = H_nc + sum(0.5 * (D[i, j] * ad[i] @ ad[j] + np.conj(D[i, j]) * a[j] @ a[i]) for i in range(n) for j in range(n))
H_int = H_nc + 1.5 * (ad[0] @ ad[1] @ a[1] @ a[0] + ad[2] @ ad[3] @ a[3] @ a[2] + ad[0] @ ad[1] @ a[3] @ a[2]
                       + ad[2] @ ad[3] @ a[1] @ a[0])
c_nc, c_bdg, c_int = connected(gibbs(H_nc)), connected(gibbs(H_bdg)), connected(gibbs(H_int))
Ppar = np.diag([(-1) ** bin(i).count("1") for i in range(2 ** n)]).astype(complex)
rho_m = (np.eye(2 ** n) + Ppar) / 2 @ gibbs(H_nc) @ (np.eye(2 ** n) + Ppar) / 2
c_proj = connected(rho_m / np.trace(rho_m))
w_, v_ = np.linalg.eigh(H_nc)
slater = np.outer(v_[:, 0], v_[:, 0].conj())                                  # pure Gaussian ground state (fixed parity)
c_pure = connected(slater)
log.step("QV1 Wick premise: connected four-point function = 0 for number-conserving, pairing and pure Slater Gaussian states "
         "(all 4096 operator orderings); interacting control nonzero -> the source selects no quartic channel",
         c_nc < 1e-13 and c_bdg < 1e-13 and c_pure < 1e-12 and c_int > 1e-2 and c_proj > 1e-4,
         f"max connected: thermal {c_nc:.1e}, pairing {c_bdg:.1e}, pure Slater {c_pure:.1e}; parity-projected mixed {c_proj:.1e} "
         f"(not the source); interacting {c_int:.2f}")

# ---------------------------------------------------------------- Spin(10) x SU(2) on V = 16 (x) 2
X = np.array([[0, 1], [1, 0]], dtype=complex); Y = np.array([[0, -1j], [1j, 0]]); Zm = np.diag([1, -1]).astype(complex)
I2 = np.eye(2, dtype=complex)


def kron(*ms):
    out = np.array([[1.0 + 0j]])
    for m in ms:
        out = np.kron(out, m)
    return out


Gam = []
for j in range(5):
    pre, post = [Zm] * j, [I2] * (4 - j)
    Gam += [kron(*pre, X, *post), kron(*pre, Y, *post)]
chir = kron(*[Zm] * 5)
plus = [i for i in range(32) if chir[i, i].real > 0]
gens10 = [(0.5j * Gam[p] @ Gam[q])[np.ix_(plus, plus)] for p, q in combinations(range(10), 2)]
gens = [np.kron(g, I2) for g in gens10] + [np.kron(np.eye(16), s / 2) for s in (X, Y, Zm)]
d16 = len(plus)
closure = all(np.allclose(g, g.conj().T) for g in gens)
weights = np.array([[np.real(np.kron((0.5j * Gam[2 * j] @ Gam[2 * j + 1])[np.ix_(plus, plus)], I2)[i, i]) for j in range(5)]
                    + [np.real(np.kron(np.eye(16), Zm / 2)[i, i])] for i in range(32)])
mono = all(np.count_nonzero(np.abs(g) > 1e-12, axis=0).max() <= 2 for g in gens)


def wedge_index(k):
    basis = list(combinations(range(32), k))
    return basis, {b: i for i, b in enumerate(basis)}


def act(g, vec, index):
    """Derivation action of a one-particle generator on a sparse wedge vector {sorted tuple: coeff}."""
    out = {}
    cols = [np.nonzero(np.abs(g[:, c]) > 1e-12)[0] for c in range(32)]
    for t, cf in vec.items():
        for slot, c in enumerate(t):
            for m in cols[c]:
                if m in t and m != c:
                    continue
                new = list(t); new[slot] = m
                perm_sign = 1
                arr = new[:]
                for i in range(len(arr)):                 # bubble sort with sign
                    for j in range(len(arr) - 1 - i):
                        if arr[j] > arr[j + 1]:
                            arr[j], arr[j + 1] = arr[j + 1], arr[j]; perm_sign = -perm_sign
                key = tuple(arr)
                out[key] = out.get(key, 0) + cf * g[m, c] * perm_sign
    return {k: v for k, v in out.items() if abs(v) > 1e-13}


def casimir_on(vecs_basis, k):
    """Matrix of C = sum_a T_a^2 on a list of wedge basis tuples closed under C (zero-weight sectors)."""
    idx = {t: i for i, t in enumerate(vecs_basis)}
    Cm = np.zeros((len(vecs_basis), len(vecs_basis)), dtype=complex)
    for i, t in enumerate(vecs_basis):
        acc = {}
        for g in gens:
            v1 = act(g, {t: 1.0}, None)
            v2 = act(g, v1, None)
            for kk, vv in v2.items():
                acc[kk] = acc.get(kk, 0) + vv
        for kk, vv in acc.items():
            if kk in idx:
                Cm[idx[kk], i] += vv
            elif abs(vv) > 1e-10:                         # nonzero weight would contradict [C, H] = 0
                raise RuntimeError("Casimir left the zero-weight sector")
    return Cm


# QV2: Lambda^2 V decomposition by Casimirs (full 496-dim space)
b2, i2 = wedge_index(2)
def full_matrix(g, basis, index):
    M = np.zeros((len(basis), len(basis)), dtype=complex)
    for i, t in enumerate(basis):
        for kk, vv in act(g, {t: 1.0}, None).items():
            M[index[kk], i] += vv
    return M
G2 = [full_matrix(g, b2, i2) for g in gens]
C10 = sum(m @ m for m in G2[:45]); C2s = sum(m @ m for m in G2[45:])
ev10 = np.round(np.linalg.eigvalsh(C10), 6); ev2 = np.round(np.linalg.eigvalsh(C2s), 6)
Ctot = C10 + 1000 * C2s
evt = np.round(np.linalg.eigvalsh(Ctot), 4)
vals, counts = np.unique(evt, return_counts=True)
log.step("QV2 Lambda^2 V = (120,3) + (10,1) + (126,1), multiplicity-free (joint Casimir eigenspaces of dims 360, 10, 126): three "
         "independent number-conserving invariants |P_(120,3)|^2, |P_10|^2, |P_126|^2",
         closure and mono and len(b2) == 496 and sorted(counts.tolist()) == [10, 126, 360],
         f"joint eigenvalue multiplicities {dict(zip(vals.tolist(), counts.tolist()))}")

# QV3: (Lambda^4 V)^G via the zero-weight Casimir kernel
zero4 = [t for t in combinations(range(32), 4) if np.allclose(weights[list(t)].sum(0), 0)]
C4 = casimir_on(zero4, 4)
ev4 = np.linalg.eigvalsh((C4 + C4.conj().T) / 2)
n_sing = int((np.abs(ev4) < 1e-8).sum())
log.step("QV3 (Lambda^4 V)^G is one-dimensional: the zero-weight Casimir kernel of Lambda^4(16 x 2) has dimension 1 "
         "(the unique Delta N = 4 invariant, hence proportional to sum_A P_A P_A)", n_sing == 1,
         f"zero-weight subspace {len(zero4)} of C(32,4) = 35960; kernel {n_sing}; next eigenvalue {np.sort(np.abs(ev4))[1]:.3f}")

# QV3b: the invariant is built from the (10,1) channel and is nonzero
w2, v2 = np.linalg.eigh(Ctot)
ten = v2[:, np.isclose(w2, 9.0, atol=1e-6)].T
quad = np.array(list(combinations(range(32), 4)))
I_, J_, K_, L_ = quad.T
def as_matrix(vec):
    M_ = np.zeros((32, 32), dtype=complex)
    for (i, j), c in zip(b2, vec):
        M_[i, j], M_[j, i] = c, -c
    return M_
Vs = [as_matrix(x) for x in ten]
cols = []
for A_ in range(10):
    for B_ in range(A_, 10):
        V_, U_ = Vs[A_], Vs[B_]
        cols.append(V_[I_, J_] * U_[K_, L_] - V_[I_, K_] * U_[J_, L_] + V_[I_, L_] * U_[J_, K_]
                    + V_[J_, K_] * U_[I_, L_] - V_[J_, L_] * U_[I_, K_] + V_[K_, L_] * U_[I_, J_])
Wfull = np.array(cols).T
zvec = np.linalg.eigh((C4 + C4.conj().T) / 2)[1][:, 0]
zfull = np.zeros(len(quad), dtype=complex)
pos = {tuple(t): i for i, t in enumerate(map(tuple, quad))}
for t, c in zip(zero4, zvec):
    zfull[pos[tuple(t)]] = c
cfit, *_ = np.linalg.lstsq(Wfull, zfull, rcond=None)
resid = np.linalg.norm(Wfull @ cfit - zfull)
log.step("QV3b the unique Delta N = 4 invariant lies in span{P_A ^ P_B} of the (10,1) channel (all 35960 components) and is nonzero: "
         "it is proportional to sum_A P_A P_A", resid < 1e-8 and np.linalg.norm(zfull) > 0.5 and ten.shape[0] == 10,
         f"fit residual {resid:.1e}; nonzero invariant components {int((np.abs(zfull) > 1e-10).sum())}")

log.step("QV4 Hermitian Z4- and G-invariant quartic vertices: real dimension 3 + 2 = 5; the exact U(1) phase of psi (symmetry of h_W "
         "and F_U, no U(1)-breaking Xi-map assumed) removes one -> 4 physical couplings (mod one-body terms); sum_A O_A^2 is one ray; "
         "all keep the IR protection (functions of rho = F psi)",
         n_sing == 1 and sorted(counts.tolist()) == [10, 126, 360], "dim = #irreps(Lambda^2 V) + 2 dim(Lambda^4 V)^G = 3 + 2")
log.step("QV5 (statement) SOURCE.QUARTIC.VERTEX.01 = NOT_SELECTED: zero connected source response and a 4-parameter symmetric space", True)

log.dump(HERE / "quartic_vertex_results.json", {
    "research_id": "SOURCE.QUARTIC.VERTEX.01",
    "verdict": "NOT_SELECTED",
    "pins": {"lambda2_multiplicities": sorted(counts.tolist()), "lambda4_singlets": n_sing, "zero_weight_dim": len(zero4),
             "hermitian_invariant_dim": 5, "physical_couplings": 4},
})
