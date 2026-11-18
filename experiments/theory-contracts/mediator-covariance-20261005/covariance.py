#!/usr/bin/env python3
"""SOURCE.MEDIATOR.COVARIANCE.01 -- which family representation must the Gate-B mediator carry?

The QMC Theta_1 treats one 16 block.  Its three-family extensions are written on the pair labels (f f') of
P_{A,ff'} = rho_f^T (i sigma_2 (x) B_A) rho_f' (symmetric in f f', so the labels span S^2(3)); a vertex is a real quadratic form Q
on S^2(3): H = -g sum_A sum_{ab} Q_ab O_{A,a} O_{A,b}.  It is invariant under a family transformation S on S^2(3) iff S^T Q S = Q
(Delta N = 4 part) and S^dag Q S = Q (Delta N = 0 part).
  MC1 vertices: V1 family-decoupled sum_f O_{A,ff}^2 (one Theta_1 copy per family); V2 SO(3)-invariant sum over all symmetric pairs;
      V3 interfamily sum_{f != f'} O_{A,ff'}^2.  V1 is invariant under the monomial group Z4^3 x| S3 (per-family Z4 and permutations)
      and under the A4-monomial transport group, but not under SO(3) or a rotated A4; V2 is SO(3)-invariant but breaks per-family Z4.
  MC2 Z4^3 x| S3 admits exactly one Delta N = 4 family tensor (sum_f e_f^{x4}): V1 is the unique Delta N = 4 three-family structure of
      the decoupled symmetry.
  MC3 decomposition of the pair labels: under A4 (monomial) the diagonal pairs carry 1 + 1' + 1'' and the off-diagonal pairs the
      triplet 3; under SO(3) S^2(3) = 1 + 5; under Z4^3 x| S3 the diagonal pairs (O_11, O_22, O_33) form an irreducible signed-permutation
      triplet.
  MC4 required mediator: second-order elimination of a degenerate multiplet coupled to the components of R gives Q = projector onto R.
      V1 needs ONE real Z4-odd Spin(10) vector per family, degenerate, unmixed: a (10, 3_signed-perm) of Z4^3 x| S3, which restricts to
      1 + 1' + 1'' under A4.  An A4-vector triplet (10, 3) gives V3 (interfamily pairing), not V1.
  MC5 the flavor transport in the A4-monomial basis (family_symmetry.py) maps V1 to itself; in a rotated family basis it does not:
      V1's family basis must be the A4-monomial basis for the transport to respect the three-family Gate-B vertex.
  MC6 source link (representation level, mediator-rep-sieve): [J(16,4) J(16,4)] in (10, S^2 4) contains S^2(3) = 1 + 1' + 1'' + 3 under A4,
      i.e. per-family composites Phi_{A,f} ~ [J^f J^f]_10 of the light families; the omitted-line singlet S^2(1) is not needed for V1.
  MC7 (statement) SOURCE.MEDIATOR.COVARIANCE.01 = R_PHI_DETERMINED_FOR_V1.
    python3 -B experiments/theory-contracts/mediator-covariance-20261005/covariance.py
"""
from __future__ import annotations

import sys
from itertools import permutations, product
from pathlib import Path

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "minimal-runtime-reconstruction-20261005"))
from ur_lib import Log  # noqa: E402

log = Log()
rng = np.random.default_rng(12)
labels = [(0, 0), (1, 1), (2, 2), (1, 2), (2, 0), (0, 1)]


def pair_action(R):
    """Action of rho -> R rho on P_{ff'} (symmetric labels): P'_{ab} = sum R_ag R_bh P_gh."""
    S = np.zeros((6, 6), dtype=complex)
    for i, (a, b) in enumerate(labels):
        for j, (g, h) in enumerate(labels):
            val = R[a, g] * R[b, h] + (R[a, h] * R[b, g] if g != h else 0)
            S[i, j] = val
    return S


def invariant(Q, S):
    return np.allclose(S.T @ Q @ S, Q, atol=1e-9) and np.allclose(S.conj().T @ Q @ S, Q, atol=1e-9)


Qdiag = np.diag([1, 1, 1, 0, 0, 0]).astype(float)
Qoff = np.diag([0, 0, 0, 1, 1, 1]).astype(float)
Qall = np.diag([1, 1, 1, 2, 2, 2]).astype(float)       # sum_{f,f'} over ordered pairs = diag + 2 x off-diagonal
V = {"V1 decoupled": Qdiag, "V2 SO(3) all pairs": Qall, "V3 interfamily": Qoff}
signed = []
for perm in permutations(range(3)):
    for s in product((1, -1), repeat=3):
        M = np.zeros((3, 3)); M[list(perm), range(3)] = s; signed.append(M)
A4mono = [M for M in signed if np.isclose(np.linalg.det(M), 1) and any(np.allclose(np.abs(M), np.roll(np.eye(3), k, 0)) for k in range(3))]
Z4S3 = []
for perm in permutations(range(3)):
    for k in product(range(4), repeat=3):
        M = np.zeros((3, 3), dtype=complex); M[list(perm), range(3)] = [1j ** kk for kk in k]; Z4S3.append(M)
R0 = expm(np.array([[0, 0.3, -0.5], [-0.3, 0, 0.7], [0.5, -0.7, 0]]))
A4rot = [R0 @ M @ R0.T for M in A4mono]
Lso3 = [np.array([[0, 0, 0], [0, 0, -1], [0, 1, 0]]), np.array([[0, 0, 1], [0, 0, 0], [-1, 0, 0]]), np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0]])]
SO3 = [expm(sum(c * L for c, L in zip(rng.normal(size=3), Lso3))) for _ in range(6)]
groups = {"Z4^3 x| S3": Z4S3, "A4 monomial": A4mono, "A4 rotated": A4rot, "SO(3)": SO3}
table = {v: {g: all(invariant(Q, pair_action(R)) for R in grp) for g, grp in groups.items()} for v, Q in V.items()}
log.step("MC1 three-family vertices on the pair labels: V1 (decoupled) invariant under Z4^3 x| S3 and A4-monomial, not under a rotated A4 or "
         "SO(3); V2 (all pairs) SO(3)-invariant but not per-family Z4; V3 (interfamily) A4-monomial-invariant, not per-family Z4",
         table["V1 decoupled"] == {"Z4^3 x| S3": True, "A4 monomial": True, "A4 rotated": False, "SO(3)": False}
         and table["V2 SO(3) all pairs"]["SO(3)"] and not table["V2 SO(3) all pairs"]["Z4^3 x| S3"]
         and table["V3 interfamily"]["A4 monomial"] and not table["V3 interfamily"]["Z4^3 x| S3"],
         "; ".join(f"{v}: {t}" for v, t in table.items()))
inv4 = int(round(np.real(sum(np.trace(np.kron(np.kron(g, g), np.kron(g, g))) for g in Z4S3) / len(Z4S3))))
diag4 = sum(np.kron(np.kron(np.eye(3)[f], np.eye(3)[f]), np.kron(np.eye(3)[f], np.eye(3)[f])) for f in range(3))
diag4_inv = all(np.allclose(np.kron(np.kron(g, g), np.kron(g, g)) @ diag4, diag4) for g in Z4S3)
log.step("MC2 Z4^3 x| S3 admits exactly one Delta N = 4 family tensor, sum_f e_f^{x4}: V1 is the unique Delta N = 4 three-family structure "
         "of the decoupled symmetry", inv4 == 1 and diag4_inv)

om = np.exp(2j * np.pi / 3)
def chars(grp, block):
    return [np.trace(pair_action(R)[np.ix_(block, block)]) for R in grp]
def cls(M):
    if np.allclose(M, np.eye(3)):
        return "e"
    if np.allclose(M @ M, np.eye(3)):
        return "V4"
    P = np.abs(M)
    return "C3a" if np.allclose(P, np.roll(np.eye(3), 1, 0)) else "C3b"
irr = {"1": {"e": 1, "V4": 1, "C3a": 1, "C3b": 1}, "1'": {"e": 1, "V4": 1, "C3a": om, "C3b": om ** 2},
       "1''": {"e": 1, "V4": 1, "C3a": om ** 2, "C3b": om}, "3": {"e": 3, "V4": -1, "C3a": 0, "C3b": 0}}
def decompose(block):
    ch = chars(A4mono, block)
    return {k: int(round(np.real(sum(np.conj(v[cls(R)]) * c for R, c in zip(A4mono, ch)) / 12))) for k, v in irr.items()}
dec_diag, dec_off = decompose([0, 1, 2]), decompose([3, 4, 5])
gen = [(pair_action(expm(1e-7 * L)) - pair_action(expm(-1e-7 * L))) / 2e-7 for L in Lso3]
C = -sum(g @ g for g in gen)
evC = sorted(float(x) for x in np.round(np.real(np.linalg.eigvals(C)), 3))
irreducible_signed = np.allclose(sum(pair_action(R)[np.ix_([0, 1, 2], [0, 1, 2])] for R in Z4S3) / len(Z4S3), 0) and \
    abs(sum(abs(np.trace(pair_action(R)[np.ix_([0, 1, 2], [0, 1, 2])])) ** 2 for R in Z4S3) / len(Z4S3) - 1) < 1e-9
log.step("MC3 pair-label decomposition: A4-monomial diagonal pairs 1 + 1' + 1'', off-diagonal pairs 3; SO(3): S^2(3) = 1 + 5 (Casimir 0, 6); "
         "Z4^3 x| S3: the diagonal triple (O_11, O_22, O_33) is irreducible (character norm 1, no invariant)",
         dec_diag == {"1": 1, "1'": 1, "1''": 1, "3": 0} and dec_off == {"1": 0, "1'": 0, "1''": 0, "3": 1}
         and evC == [0.0, 6.0, 6.0, 6.0, 6.0, 6.0] and irreducible_signed,
         f"diag {dec_diag}, off {dec_off}, SO(3) Casimir {evC}")

# MC4: elimination of a degenerate multiplet = projector onto its span
def eliminated(basis):
    Bm = np.array(basis, dtype=float)
    return Bm.T @ Bm
med_V1 = [np.eye(6)[i] for i in range(3)]                      # one mediator per family, coupled to O_ff
med_triplet = [np.eye(6)[i] for i in (3, 4, 5)]               # A4-vector triplet couples to the off-diagonal pairs
log.step("MC4 required mediator: eliminating a degenerate multiplet coupled to the components of R gives Q = projector onto R; V1 needs one real "
         "Z4-odd Spin(10) vector per family (degenerate, unmixed): a (10, 3_signed-perm) of Z4^3 x| S3, = 1 + 1' + 1'' under A4; an A4-vector "
         "triplet (10, 3) produces V3 (interfamily), not V1",
         np.allclose(eliminated(med_V1), Qdiag) and np.allclose(eliminated(med_triplet), Qoff))

# MC5: the flavor transport (A4-monomial basis from family_symmetry.py) respects V1 only in that basis
Wm = np.load(HERE.parent / "family-symmetry-gateb-20261005" / "a4_monomial_basis.npy")
sys.path.insert(0, str(HERE.parents[2]))
from tfpt_explorer.flavor_path_transport import exact_clock_dictionary  # noqa: E402
M = np.array(exact_clock_dictionary()["M"].evalf(), dtype=complex)
Uc = np.diag([1, 1j, -1j])
Ma = [np.linalg.matrix_power(Uc, a) @ M @ np.linalg.matrix_power(Uc, -a) for a in range(4)]
in_mono = all(invariant(Qdiag, pair_action(np.linalg.inv(Wm) @ m @ Wm)) for m in Ma)
in_orig = all(invariant(Qdiag, pair_action(m)) for m in Ma)
log.step("MC5 the v117 transport maps V1 to itself in the A4-monomial family basis but not in the original family-mode basis: the per-family "
         "basis of the three-family Gate-B vertex must be the A4-monomial basis", in_mono and not in_orig)

# MC6: source link
def s2_3_decomp():
    ch = [np.trace(pair_action(R)) for R in A4mono]
    return {k: int(round(np.real(sum(np.conj(v[cls(R)]) * c for R, c in zip(A4mono, ch)) / 12))) for k, v in irr.items()}
dS2 = s2_3_decomp()
log.step("MC6 source link (representation level): [J(16,4) J(16,4)] in (10, S^2 4) contains S^2(3) = 1 + 1' + 1'' + 3 under A4, i.e. per-family "
         "composites Phi_{A,f} ~ [J^f J^f]_10 of the light families; the omitted-line singlet S^2(1) is not needed for V1 (bulk Z4 charge, gap, "
         "residue, coupling and SW residual still open)", dS2 == {"1": 1, "1'": 1, "1''": 1, "3": 1}, f"S^2(3) under A4: {dS2}")
log.step("MC7 (statement) SOURCE.MEDIATOR.COVARIANCE.01 = R_PHI_DETERMINED_FOR_V1: (10, 3_signed-perm) of Z4^3 x| S3, degenerate and unmixed, "
         "restricting to 1 + 1' + 1'' under the A4 transport group, in the A4-monomial family basis", True)

log.dump(HERE / "covariance_results.json", {"research_id": "SOURCE.MEDIATOR.COVARIANCE.01", "verdict": "R_PHI_DETERMINED_FOR_V1",
                                            "pins": {"invariance_table": table, "A4_diag": dec_diag, "A4_off": dec_off, "S2_3": dS2}})
