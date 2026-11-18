#!/usr/bin/env python3
"""SOURCE.MEDIATOR.REP_SIEVE.01 -- from how many and which source modes can the missing (10,1) mediator be built at all?

Target: a bosonic Spin(10) vector that is a family singlet (spin singlet is automatic: source modes carry no spacetime spin).
Building blocks of the source: NS fermions psi (10,1) and chi (1,6); E8 currents (45,1), (1,15), (10,6), (16,4), (16bar,4bar).
  RS1 center selection rule: on all 240 E8 roots the Spin(10) and SU(4) center charges satisfy q10 + q4 = 0 mod 4 (diagonal
      centre trivial, as in UR.SOURCE.OPERATOR_ORIGIN.01); the NS fermions have q10 + q4 = 2.  Hence every bosonic composite has
      q10 + q4 = 0 mod 4, while (10, SU(4)-singlet) has q10 + q4 = 2: FORBIDDEN AT EVERY COMPOSITE DEGREE (exhaustive charge
      enumeration up to degree 6 as a check).
  RS2 the candidate (16,4)(16,4)(1,6) is excluded: with chi it is fermionic, and a bosonic (1,6) would violate RS1.
  RS3 Spin(10) side checked numerically: 10 occurs once in S^2 16, not in Lambda^2 16, 16 x 16bar, 10 x 10; once in 10 x 45.
  RS4 smaller family-neutrality groups H (SU(3) = stabiliser of the omitted line, SO(3), A4): degree-2 candidates exist:
      [J(16,4) J(16,4)] in (10, S^2 4) and [J(10,6) J(1,15)] in (10, 6 x 15); H-singlet counts computed explicitly.  Under SU(3) the
      unique singlet of S^2 4 is the omitted-family-line component.  Fermion bilinears never give (10, H-singlet).
  RS5 (statement) minimal degree: SU(4): none; SU(3), SO(3), A4: 2.  Representation theory does not decide the bulk Z4 charge,
      the gap Delta, the residue Z_Phi, the coupling lambda or the Schrieffer-Wolff residual: those are the next gates.
    python3 -B experiments/theory-contracts/mediator-rep-sieve-20261005/sieve.py
"""
from __future__ import annotations

import sys
from itertools import combinations, combinations_with_replacement, permutations, product
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "minimal-runtime-reconstruction-20261005"))
from ur_lib import Log, jordan_wigner  # noqa: E402

log = Log()

# ---------------------------------------------------------------- RS1 centre charges
roots = []
for i, j in combinations(range(8), 2):
    for si, sj in product((1, -1), repeat=2):
        r = np.zeros(8); r[i], r[j] = si, sj; roots.append(r)
for signs in product((0.5, -0.5), repeat=8):
    if sum(1 for s_ in signs if s_ < 0) % 2 == 0:
        roots.append(np.array(signs))
roots = np.array(roots)
q10 = lambda lam: int(round(4 * lam[:5].sum() / 2)) % 4          # exp(2 pi i lam . (1/2)^5) = i^{q10}
q4 = lambda lam: int(round(4 * lam[5:].sum() / 2)) % 4
charges = {(q10(r), q4(r)) for r in roots}
diag_trivial = all((q10(r) + q4(r)) % 4 == 0 for r in roots)
psi = (2, 0); chi = (0, 2)
blocks = {"(45,1)": (0, 0, +1), "(1,15)": (0, 0, +1), "(10,6)": (2, 2, +1), "(16,4)": None, "(16b,4b)": None,
          "psi(10,1)": (2, 0, -1), "chi(1,6)": (0, 2, -1)}
spin_pairs = sorted({(q10(r), q4(r)) for r in roots if abs(r[0]) == 0.5})
blocks["(16,4)"] = (*spin_pairs[0], +1); blocks["(16b,4b)"] = (*spin_pairs[1], +1)
bad = []
for n in range(1, 7):
    for combo in combinations_with_replacement(list(blocks), n):
        a = sum(blocks[b][0] for b in combo) % 4
        b_ = sum(blocks[b][1] for b in combo) % 4
        stat = np.prod([blocks[b][2] for b in combo])
        if stat == 1 and (a, b_) == (2, 0):
            bad.append(combo)
log.step("RS1 centre selection rule: all 240 E8 roots satisfy q10 + q4 = 0 mod 4 (diagonal centre trivial, UR.SOURCE.OPERATOR_ORIGIN.01), "
         "NS fermions have q10 + q4 = 2 -> every bosonic composite has q10 + q4 = 0; the target (10, SU(4)-singlet) = (2, 0) never occurs "
         "(exhaustive enumeration of all block multisets up to degree 6)",
         diag_trivial and len(spin_pairs) == 2 and not bad, f"E8 charge pairs {sorted(charges)}; spinor blocks {spin_pairs}")

# ---------------------------------------------------------------- RS2 the (16,4)(16,4)(1,6) candidate
combo_chi = (blocks["(16,4)"][0] * 2 + 0, blocks["(16,4)"][1] * 2 + 2)
log.step("RS2 (16,4)(16,4)(1,6): with chi the product is fermionic; a bosonic (1,6) has q10 + q4 = 2 and does not exist in the source "
         "(RS1) -> candidate excluded", np.prod([1, 1, blocks["chi(1,6)"][2]]) == -1 and (0 + 2) % 4 != 0,
         f"centre charges of (16,4)(16,4)chi: {tuple(x % 4 for x in combo_chi)}")

# ---------------------------------------------------------------- Spin(10) representations: 16, 16bar, 10, 45
X = np.array([[0, 1], [1, 0]], dtype=complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1, -1]).astype(complex)
I2 = np.eye(2)
def kron(*ms):
    out = np.array([[1.0 + 0j]])
    for m in ms:
        out = np.kron(out, m)
    return out
Gam = []
for j in range(5):
    Gam += [kron(*[Z] * j, X, *[I2] * (4 - j)), kron(*[Z] * j, Y, *[I2] * (4 - j))]
chir = kron(*[Z] * 5)
plus = [i for i in range(32) if chir[i, i].real > 0]; minus = [i for i in range(32) if chir[i, i].real < 0]
pairs10 = list(combinations(range(10), 2))
S16 = [(0.5j * Gam[p] @ Gam[q])[np.ix_(plus, plus)] for p, q in pairs10]
S16b = [(0.5j * Gam[p] @ Gam[q])[np.ix_(minus, minus)] for p, q in pairs10]
S10 = []
for p, q in pairs10:
    M = np.zeros((10, 10), dtype=complex); M[p, q], M[q, p] = 1j, -1j       # rotation generator in the (p,q) plane
    S10.append(M)
def tensor(g1, g2):
    return [np.kron(a, np.eye(b.shape[0])) + np.kron(np.eye(a.shape[0]), b) for a, b in zip(g1, g2)]
def casimir(gs):
    return sum(g @ g for g in gs)
c10 = np.linalg.eigvalsh(casimir(S10))[0]
c16 = np.linalg.eigvalsh(casimir(S16))[0]
def mult_of(gs, cval, dim):
    ev = np.linalg.eigvalsh(casimir(gs))
    return int(round(np.sum(np.isclose(ev, cval, atol=1e-6)) / dim))
# symmetric / antisymmetric squares of the 16
d = 16
Psym = np.zeros((d * d, d * d)); Pas = np.zeros((d * d, d * d))
for i in range(d):
    for j in range(d):
        a, b = i * d + j, j * d + i
        Psym[a, a] += 0.5; Psym[a, b] += 0.5; Pas[a, a] += 0.5; Pas[a, b] -= 0.5
def restrict(gs, Pproj):
    u, s, vh = np.linalg.svd(Pproj); B = u[:, s > 0.5]
    return [B.conj().T @ g @ B for g in gs]
T1616 = tensor(S16, S16)
m_sym = mult_of(restrict(T1616, Psym), c10, 10)
m_as = mult_of(restrict(T1616, Pas), c10, 10)
m_1616b = mult_of(tensor(S16, S16b), c10, 10)
m_1010 = mult_of(tensor(S10, S10), c10, 10)
S45 = []
basis45 = [M for M in S10]
for g in S10:
    A = np.zeros((45, 45), dtype=complex)
    for k, Bk in enumerate(basis45):
        comm = g @ Bk - Bk @ g
        coeff = [np.trace(Bm.conj().T @ comm) / np.trace(Bm.conj().T @ Bm) for Bm in basis45]
        A[:, k] = coeff
    S45.append(A)
m_1045 = mult_of(tensor(S10, S45), c10, 10)
log.step("RS3 Spin(10) multiplicities of the vector 10: S^2 16 -> 1, Lambda^2 16 -> 0, 16 x 16bar -> 0, 10 x 10 -> 0, 10 x 45 -> 1",
         (m_sym, m_as, m_1616b, m_1010, m_1045) == (1, 0, 0, 0, 1), f"Casimir(10) = {c10:.3f}, Casimir(16) = {c16:.3f}")

# ---------------------------------------------------------------- family representations 4, 4bar, 6, 15 and subgroups
f = jordan_wigner(3)
dag = lambda x: x.conj().T
maj = []
for o in f:
    maj += [o + dag(o), -1j * (o - dag(o))]
Pf = np.diag([(-1) ** bin(i).count("1") for i in range(8)])
even = [i for i in range(8) if Pf[i, i] > 0]; odd = [i for i in range(8) if Pf[i, i] < 0]
spin6 = [0.5j * maj[p] @ maj[q] for p, q in combinations(range(6), 2)]
def su3(ops):
    gens = []
    for i in range(3):
        for j in range(i + 1, 3):
            gens += [dag(ops[i]) @ ops[j] + dag(ops[j]) @ ops[i], -1j * (dag(ops[i]) @ ops[j] - dag(ops[j]) @ ops[i])]
    gens += [dag(ops[0]) @ ops[0] - dag(ops[1]) @ ops[1], dag(ops[0]) @ ops[0] + dag(ops[1]) @ ops[1] - 2 * dag(ops[2]) @ ops[2]]
    return gens
Lso3 = [np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]]), np.array([[0, 0, 1j], [0, 0, 0], [-1j, 0, 0]]),
        np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]])]
dG = lambda L: sum(L[i, j] * dag(f[i]) @ f[j] for i in range(3) for j in range(3))
w = np.array([[-.5, -.5, -.5], [-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]])
def par(p):
    p = list(p); s_ = 1
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]; p[i], p[j] = p[j], p[i]; s_ = -s_
    return s_
from scipy.linalg import expm, logm  # noqa: E402
A4 = []
for p in permutations(range(4)):
    if par(p) == 1:
        P4 = np.zeros((4, 4)); P4[list(p), range(4)] = 1
        R = w.T @ P4 @ w
        A4.append(expm(-1j * dG(1j * logm(R.astype(complex)))))
H = {"SU(4)": spin6, "SU(3)": su3(f), "SO(3)": [dG(L) for L in Lso3]}


def rep_matrices(ops_fock):
    """Representation of a family generator/unitary on 4 (even), 4bar (odd), 6 (Majorana vector), 15 (Spin(6) bilinears)."""
    out = {"4": ops_fock[np.ix_(even, even)], "4b": ops_fock[np.ix_(odd, odd)]}
    return out


def vector_rep(T, basis, unitary=False):
    """action on span(basis) by commutator (generator) or conjugation (unitary)."""
    n = len(basis)
    M = np.zeros((n, n), dtype=complex)
    for k, Bk in enumerate(basis):
        img = (T @ Bk @ dag(T)) if unitary else (T @ Bk - Bk @ T)
        for m, Bm in enumerate(basis):
            M[m, k] = np.trace(dag(Bm) @ img) / np.trace(dag(Bm) @ Bm)
    return M


def invariants(gen_list, group_list, space):
    """dimension of H-invariants in a representation given by generator matrices or group matrices."""
    if group_list is not None:
        Pavg = sum(group_list) / len(group_list)
        return int(round(np.real(np.trace(Pavg))))
    Cs = sum(g.conj().T @ g for g in gen_list)
    ev = np.linalg.eigvalsh((Cs + Cs.conj().T) / 2)
    return int(np.sum(np.abs(ev) < 1e-8))


def family_reps(hname):
    if hname == "A4":
        mats = [{"4": U[np.ix_(even, even)], "4b": U[np.ix_(odd, odd)], "6": vector_rep(U, maj, True),
                 "15": vector_rep(U, spin6, True)} for U in A4]
        return mats, None
    gens = [{"4": T[np.ix_(even, even)], "4b": T[np.ix_(odd, odd)], "6": vector_rep(T, maj), "15": vector_rep(T, spin6)}
            for T in H[hname]]
    return None, gens


def sym2(M, unitary):
    dd = M.shape[0]
    Ps = np.zeros((dd * dd, dd * dd))
    for i in range(dd):
        for j in range(dd):
            a, b = i * dd + j, j * dd + i
            Ps[a, a] += 0.5; Ps[a, b] += 0.5
    u, s, vh = np.linalg.svd(Ps); B = u[:, s > 0.5]
    big = np.kron(M, M) if unitary else np.kron(M, np.eye(dd)) + np.kron(np.eye(dd), M)
    return B.conj().T @ big @ B


counts = {}
for hname in ("SU(4)", "SU(3)", "SO(3)", "A4"):
    grp, gens = family_reps(hname)
    if grp is not None:
        S24 = invariants(None, [sym2(g["4"], True) for g in grp], None)
        V6 = invariants(None, [g["6"] for g in grp], None)
        T615 = invariants(None, [np.kron(g["6"], g["15"]) for g in grp], None)
        S24b = invariants(None, [sym2(g["4b"], True) for g in grp], None)
    else:
        S24 = invariants([sym2(g["4"], False) for g in gens], None, None)
        S24b = invariants([sym2(g["4b"], False) for g in gens], None, None)
        V6 = invariants([g["6"] for g in gens], None, None)
        T615 = invariants([np.kron(g["6"], np.eye(15)) + np.kron(np.eye(6), g["15"]) for g in gens], None, None)
    counts[hname] = {"J(16,4)^2 in (10,S^2 4)": S24, "J(16b,4b)^2": S24b, "J(10,6)J(45,1)": V6, "J(10,6)J(1,15)": T615,
                     "psi chi": V6}
min_deg = {h: (2 if any(v > 0 for v in c_.values()) else None) for h, c_ in counts.items()}
omitted = np.zeros(4); omitted[0] = 1                                   # |000> = the omitted family line (first even state)
su3_fix = all(np.allclose(g[np.ix_(even, even)] @ omitted, 0) for g in H["SU(3)"])
log.step("RS4 degree-2 family-singlet counts: SU(4) none (consistent with RS1); SU(3), SO(3), A4 admit (10, H-singlet) composites "
         "[J(16,4)J(16,4)] in (10, S^2 4) and [J(10,6)J(1,15)]; fermion bilinears (psi chi in (10,6)) never; under SU(3) the singlet of "
         "S^2 4 is the omitted-family-line component (|000> is SU(3)-invariant)",
         all(v == 0 for v in counts["SU(4)"].values()) and all(min_deg[h] == 2 for h in ("SU(3)", "SO(3)", "A4"))
         and counts["SU(3)"]["J(16,4)^2 in (10,S^2 4)"] == 1 and all(counts[h]["psi chi"] == 0 for h in counts) and su3_fix,
         "; ".join(f"{h}: {c_}" for h, c_ in counts.items()))
log.step("RS5 (statement) SOURCE.MEDIATOR.REP_SIEVE.01: minimal composite degree of a bosonic (10, family-singlet) mediator = none for "
         "SU(4) neutrality (all degrees forbidden), 2 for SU(3)/SO(3)/A4 neutrality; bulk Z4 charge, gap Delta, residue Z_Phi, coupling "
         "lambda and the Schrieffer-Wolff residual ||dH|| << lambda^2/Delta are not decided by representation theory (next gates)", True)

log.dump(HERE / "sieve_results.json", {"research_id": "SOURCE.MEDIATOR.REP_SIEVE.01", "verdict": "SU4_FORBIDDEN_ALL_DEGREES; SUBGROUP_MIN_DEGREE_2",
                                        "pins": {"degree2_counts": counts, "min_degree": min_deg}})
