#!/usr/bin/env python3
"""SOURCE.FAMILY.SYMMETRY.AT_GATEB.01 -- which family group can be an exact symmetry when the Gate-B interaction acts?

Premise: the three-family Theta_1 keeps the Delta N = +-4 part of -g sum_A O_A^2 (U(1) -> Z4); a symmetric gap of the mirror 16
cannot keep U(1) exact (sum q = sum q^3 = 16).
  FS1 anomaly premise: U(1) sums of one 16 are 16 (nonzero); Z4 is the residual symmetry of Theta_1 (SG5 cited for its index).
  FS2 centre criterion: scalar elements of a family group act on the light families as bulk U(1) phases, so an exact family group G
      at Gate B must satisfy G cap U(1)*I in mu_4.  Necessary condition checked directly: invariants of 3^{x4} (the family part of
      any Delta N = 4 operator) vanish for U(3) and SU(3) (centre Z3), exist for SO(3), A4, the signed permutations B3 and
      Z4^3 x| S3.
  FS3 source layers: L0 seam: Spin(6) = SU(4) on the family Majoranas; L1 matter selection (family number 2): stabiliser u(3)
      (dimension 9, contains the overall U(1)); L2 flavor transport: the centraliser of the v117 monodromies in u(3) is one-dimensional
      (Schur), the transport is of real type, and a basis exists in which all twelve monodromies are signed permutation matrices
      (A4-monomial basis); the deck relabeling a -> a+1 is a combined mark (x) family symmetry with an intertwiner that is a real
      signed permutation up to phase in that basis.
  FS4 consequence: exact SU(4), U(3) or SU(3) family symmetry at Gate B is excluded; the family structure must already be real
      (G in O(3) x mu_4); in the source the real structure comes from the flavor transport.  Ordering: seam SU(4) -> matter
      selection U(3) -> real flavor structure (transport) -> Gate B -> flavor/Yukawa.
  FS5 (statement) SOURCE.FAMILY.SYMMETRY.AT_GATEB.01 = REAL_FAMILY_GROUP_REQUIRED (conditional on Theta_1's Delta N = 4 vertex for three
      light families).
    python3 -B experiments/theory-contracts/family-symmetry-gateb-20261005/family_symmetry.py
"""
from __future__ import annotations

import sys
from itertools import combinations, permutations, product
from pathlib import Path

import numpy as np
from scipy.linalg import expm, logm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent / "minimal-runtime-reconstruction-20261005"))
sys.path.insert(0, str(ROOT))
from ur_lib import Log, jordan_wigner  # noqa: E402
from tfpt_explorer.flavor_path_transport import exact_clock_dictionary  # noqa: E402

log = Log()
rng = np.random.default_rng(8)

# ---------------------------------------------------------------- FS1 anomaly premise
q = np.ones(16)
log.step("FS1 premise: one 16 with unit fermion number has sum q = sum q^3 = 16 != 0, so a symmetric mirror gap cannot keep U(1); Theta_1 "
         "keeps Z4 (psi -> i psi) through its Delta N = 4 part (QV3b: unique invariant proportional to sum_A P_A P_A; SG5: index 16 = 0 mod 16)",
         q.sum() == 16 and (q ** 3).sum() == 16)

# ---------------------------------------------------------------- FS2 centre criterion: invariants of 3^{x4}
def tensor4(gens=None, group=None, d=3):
    I = np.eye(d)
    if group is not None:
        return int(round(np.real(sum(np.trace(np.kron(np.kron(g, g), np.kron(g, g))) for g in group) / len(group))))
    ops = []
    for T in gens:
        terms = []
        for slot in range(4):
            mats = [I] * 4; mats[slot] = T
            terms.append(np.kron(np.kron(mats[0], mats[1]), np.kron(mats[2], mats[3])))
        ops.append(sum(terms))
    C = sum(o.conj().T @ o for o in ops)
    return int(np.sum(np.abs(np.linalg.eigvalsh((C + C.conj().T) / 2)) < 1e-8))


def su3_gens():
    E = lambda i, j: np.eye(3)[:, [i]] @ np.eye(3)[[j], :]
    g = []
    for i in range(3):
        for j in range(i + 1, 3):
            g += [E(i, j) + E(j, i), -1j * (E(i, j) - E(j, i))]
    g += [np.diag([1, -1, 0]).astype(complex), np.diag([1, 1, -2]).astype(complex)]
    return g


Lso3 = [np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]]), np.array([[0, 0, 1j], [0, 0, 0], [-1j, 0, 0]]),
        np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]])]
signed_perms = []
for perm in permutations(range(3)):
    for s in product((1, -1), repeat=3):
        M = np.zeros((3, 3)); M[list(perm), range(3)] = s
        signed_perms.append(M.astype(complex))
A4mono = [M for M in signed_perms if np.isclose(np.linalg.det(M), 1) and
          any(np.allclose(np.abs(M), np.roll(np.eye(3), k, 0)) for k in range(3))]
Z4S3 = []
for perm in permutations(range(3)):
    for k in product(range(4), repeat=3):
        M = np.zeros((3, 3), dtype=complex); M[list(perm), range(3)] = [1j ** kk for kk in k]
        Z4S3.append(M)
inv = {"U(3)": tensor4(gens=su3_gens() + [np.eye(3, dtype=complex)]), "SU(3)": tensor4(gens=su3_gens()),
       "SO(3)": tensor4(gens=Lso3), "A4 (monomial)": tensor4(group=A4mono), "B3 (signed perms)": tensor4(group=signed_perms),
       "Z4^3 x| S3": tensor4(group=Z4S3)}
scalars = {"U(3)": "U(1)", "SU(3)": "Z3", "SO(3)": "{1}", "A4 (monomial)": "{1}", "B3 (signed perms)": "{+-1}", "Z4^3 x| S3": "mu_4"}
log.step("FS2 centre criterion: scalar elements act as bulk U(1) phases, so G cap U(1)*I must lie in mu_4; invariants of 3^{x4} (family part "
         "of every Delta N = 4 operator) vanish for U(3) and SU(3), exist for SO(3), A4, B3 and Z4^3 x| S3",
         inv["U(3)"] == 0 and inv["SU(3)"] == 0 and all(inv[k] > 0 for k in ("SO(3)", "A4 (monomial)", "B3 (signed perms)", "Z4^3 x| S3"))
         and len(A4mono) == 12 and len(signed_perms) == 48 and len(Z4S3) == 384,
         "; ".join(f"{k}: inv {v}, scalars {scalars[k]}" for k, v in inv.items()))

# ---------------------------------------------------------------- FS3 source layers
f = jordan_wigner(3)
dag = lambda x: x.conj().T
maj = []
for o in f:
    maj += [o + dag(o), -1j * (o - dag(o))]
spin6 = [0.5j * maj[p] @ maj[q_] for p, q_ in combinations(range(6), 2)]
Nf = sum(dag(o) @ o for o in f)
P2 = np.diag([1.0 if bin(i).count("1") == 2 else 0.0 for i in range(8)])
A = np.array([(T @ P2 - P2 @ T).ravel() for T in spin6]).T
stab_dim = 15 - np.linalg.matrix_rank(A, tol=1e-9)
u, s_, vh = np.linalg.svd(A)
stab = [sum(c * T for c, T in zip(vh[k].conj(), spin6)) for k in range(np.linalg.matrix_rank(A, tol=1e-9), 15)]
contains_U1 = np.linalg.matrix_rank(np.array([X.ravel() for X in stab] + [(Nf - 1.5 * np.eye(8)).ravel()]), tol=1e-9) == len(stab)
M = np.array(exact_clock_dictionary()["M"].evalf(), dtype=complex)
Uc = np.diag([1, 1j, -1j])
Ma = [np.linalg.matrix_power(Uc, a) @ M @ np.linalg.matrix_power(Uc, -a) for a in range(4)]
GM = [expm(-1j * sum((1j * logm(m))[i, j] * dag(f[i]) @ f[j] for i in range(3) for j in range(3))) for m in Ma]
B = np.array([np.concatenate([(X @ G - G @ X).ravel() for G in GM]) for X in stab]).T
cent_dim = len(stab) - np.linalg.matrix_rank(B, tol=1e-8)
# group, V4, monomial basis
grp = [np.eye(3, dtype=complex)]
frontier = list(grp)
while frontier:
    new = []
    for g in frontier:
        for h in Ma:
            x = h @ g
            if not any(np.allclose(x, e, atol=1e-9) for e in grp):
                grp.append(x); new.append(x)
    frontier = new
inv2 = [g for g in grp if np.allclose(g @ g, np.eye(3), atol=1e-9) and not np.allclose(g, np.eye(3))]
eqs = np.array([np.concatenate([(m.T @ E @ m - E).ravel() for m in Ma]) for E in
                [np.eye(3)[:, [i]] @ np.eye(3)[[j], :] for i in range(3) for j in range(3)]]).T
uJ, sJ, vJ = np.linalg.svd(eqs)
Jsym = vJ[-1].conj().reshape(3, 3)
real_type = np.sum(sJ < 1e-9) == 1 and np.allclose(Jsym, Jsym.T, atol=1e-9)
w_, W = np.linalg.eig(inv2[0] + 2 * inv2[1])
W = W / np.linalg.norm(W, axis=0)
cyc = next(g for g in grp if not np.allclose(g @ g, np.eye(3), atol=1e-9) and np.allclose(np.linalg.matrix_power(g, 3), np.eye(3), atol=1e-9))
c_mono = np.linalg.inv(W) @ cyc @ W
perm_idx = [int(np.argmax(np.abs(c_mono[:, j]))) for j in range(3)]
phases = np.ones(3, dtype=complex)
j = 0
for _ in range(2):
    i = perm_idx[j]
    phases[i] = phases[j] * c_mono[i, j] / abs(c_mono[i, j])
    j = i
Wm = W @ np.diag(phases)
mono = [np.linalg.inv(Wm) @ g @ Wm for g in grp]
all_signed = all(np.allclose(np.abs(m), np.round(np.abs(m)), atol=1e-8) and np.allclose(m.imag, 0, atol=1e-8)
                 and np.allclose(np.round(m.real), m.real, atol=1e-8) for m in mono)
deck_V = np.linalg.inv(Wm) @ Uc @ Wm
ph = deck_V[np.unravel_index(np.argmax(np.abs(deck_V)), deck_V.shape)]
deck_real = np.allclose((deck_V / (ph / abs(ph))).imag, 0, atol=1e-8)
deck_ok = all(np.allclose(Uc @ Ma[a] @ np.linalg.inv(Uc), Ma[(a + 1) % 4]) for a in range(4))
log.step("FS3 source layers: L0 seam Spin(6) = SU(4) (15 generators); L1 matter selection (family number 2) leaves u(3) (dimension 9, contains "
         "the overall U(1)); L2 flavor transport: centraliser of the v117 monodromies in u(3) is one-dimensional (Schur), the transport is of "
         "real type, all twelve monodromies are signed permutation matrices in the A4-monomial basis, and the deck intertwiner U (a -> a+1) is "
         "a real matrix up to a phase there",
         len(spin6) == 15 and stab_dim == 9 and contains_U1 and cent_dim == 1 and len(grp) == 12 and len(inv2) == 3 and real_type
         and all_signed and deck_ok and deck_real,
         f"stabiliser dim {stab_dim}, centraliser dim {cent_dim}, |A4| {len(grp)}, real type {real_type}, monomial {all_signed}")
log.step("FS4 consequence: exact SU(4), U(3) or SU(3) family symmetry at Gate B is excluded (FS2); the family structure must already be real "
         "(G in O(3) x mu_4); in the source the real structure is supplied by the flavor transport (FS3) -> ordering: seam SU(4) -> matter "
         "selection U(3) -> real flavor structure -> Gate B -> flavor/Yukawa; the SU(3)-singlet sieve candidate (omitted line) refers to a "
         "group that cannot be exact at Gate B", inv["SU(3)"] == 0 and real_type)
log.step("FS5 (statement) SOURCE.FAMILY.SYMMETRY.AT_GATEB.01 = REAL_FAMILY_GROUP_REQUIRED (conditional on Theta_1's Delta N = 4 vertex for "
         "three light families)", True)

log.dump(HERE / "family_symmetry_results.json", {
    "research_id": "SOURCE.FAMILY.SYMMETRY.AT_GATEB.01", "verdict": "REAL_FAMILY_GROUP_REQUIRED",
    "pins": {"invariants_3x4": inv, "stabiliser_dim": int(stab_dim), "centraliser_dim": int(cent_dim)},
    "A4_monomial_basis_columns": [[[float(z.real), float(z.imag)] for z in Wm[:, k]] for k in range(3)],
})
np.save(HERE / "a4_monomial_basis.npy", Wm)
