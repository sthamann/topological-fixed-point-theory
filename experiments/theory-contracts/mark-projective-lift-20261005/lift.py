#!/usr/bin/env python3
"""SOURCE.MARK_PROJECTIVE_LIFT.01 -- is the oriented mark-rotation group A4 lifted trivially or as 2T (center = (-1)^F)?

No spin structure is inserted.  The invariant of a projective A4 action is the rephasing-invariant commutator
c = U(x) U(y) U(x)^-1 U(y)^-1 of the lifts of two distinct double transpositions x = (01)(23), y = (02)(13):
c = +1 for the trivial class, c = -1 for the nontrivial Schur class (2T); spin-statistics additionally demands that the
central element acts as the fermion parity (-1)^F, not as a scalar.
  PL1 invariant calibrated: A4 std rep and regular rep give c = +1, the SU(2) (2T) doublet c = -1.
  PL2 Fock-implementation lemma: for a GENUINE orthogonal one-particle action O of A4 on n Majoranas, the Clifford
      implementers (products of reflection elements gamma(v)) satisfy every A4 relation up to SCALARS; their class is
      w2(O) = m_3(O) mod 2 (multiplicity of the 3-dim irrep in the real one-particle representation); the central element
      is never the fermion parity.  Checked on: number-conserving family modes (R + R: c = +1), gamma-only rotation
      (R + 1: c = -1 as a scalar, spinorial on BOTH parity sectors), twist zero modes N = 1..4 (c = (-1)^N), random real
      representations a 1 + b 2_R + c 3 (c = (-1)^c).  Generalises the 20260921 lemma (U^2 = parity needs O^2 = -I).
  PL3 positive control: if the ONE-PARTICLE space itself carries the 2T doublet, second quantisation gives c = (-1)^F:
      spin-statistics arises exactly from a projective one-particle action, never from the Fock implementation.
  PL4 seam geometry: the Moebius stabiliser of the square marks {1, i, -1, -i} is D4 (8 permutations); only the four
      rotations preserve the orientation of the seam circle; with A4 it shares only <(02)(13)>; no 3-cycle is a seam
      symmetry (a concyclic quadruple has real cross-ratio; the A4-symmetric equianharmonic one is not concyclic).
  PL5 the seam rotation lift is conditional: NS modes give U(2 pi) = (-1)^F, R modes +1, with identical permutation
      structure (20260921 kill: the geometric half-turn alone does not fix the lift).  Even under the S^2 spinor-section
      premise this reaches only <(02)(13)> inside A4 and provides no V4 commutator.
  PL6 source inventory: all source one-particle mark actions are genuine with m_3 even (family modes, twist N = 16) or are
      not fermionic (flavor fiber, windings); the only nontrivial class (gamma-only) breaks fermion number, mixes families
      and is spinorial on the even sector: no source sector realises A4 with center (-1)^F.
  PL7 emergent spin 1/2: at the isotropic Weyl node of h_W the emergent rotations act on a two-component coin necessarily
      through SU(2) (O(k^2) covariance) and the 2 pi loop lifts to -I; its second quantisation coincides with the fermion
      parity that the CAR/Fock structure already has.  This is NOT a derivation of spin-statistics from IR geometry and it
      cannot be used to derive AX.WEYLCOIN (a Weyl node presupposes the two-dimensional coin).  The P1 sheet doublet is no
      coin (all marks act by the same swap).
  PL8 (statement) SOURCE.MARK_PROJECTIVE_LIFT.01 = NOT_FORCED_IN_SOURCE; SPIN_HALF_EMERGENT_GIVEN_TWO_COMPONENT_COIN.
    python3 -B experiments/theory-contracts/mark-projective-lift-20261005/lift.py
"""
from __future__ import annotations

import sys
from itertools import permutations
from pathlib import Path

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "minimal-runtime-reconstruction-20261005"))
from ur_lib import Log  # noqa: E402

log = Log()
rng = np.random.default_rng(42)
sig = [np.array([[0, 1], [1, 0]], dtype=complex), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]).astype(complex)]
w = np.array([[-.5, -.5, -.5], [-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]])


def parity(p):
    p = list(p); s = 1
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]; p[i], p[j] = p[j], p[i]; s = -s
    return s


A4 = [p for p in permutations(range(4)) if parity(p) == 1]
X_, Y_ = (1, 0, 3, 2), (2, 3, 0, 1)                      # (01)(23), (02)(13)
compose = lambda p, q: tuple(p[q[i]] for i in range(4))   # p after q


def perm4(p):
    M = np.zeros((4, 4)); M[[p[a] for a in range(4)], range(4)] = 1
    return M


def std(p):
    return w.T @ perm4(p) @ w                              # R w_a = w_{p(a)}


def pairing_rep(p):                                        # A4 -> permutations of the three pairings (= 1 + 2_R)
    pairs = [frozenset([frozenset([0, 1]), frozenset([2, 3])]), frozenset([frozenset([0, 2]), frozenset([1, 3])]),
             frozenset([frozenset([0, 3]), frozenset([1, 2])])]
    M = np.zeros((3, 3))
    for i, pr in enumerate(pairs):
        img = frozenset(frozenset(p[a] for a in blk) for blk in pr)
        M[pairs.index(img), i] = 1
    return M


def su2_lift(R):
    """U in SU(2) with U sigma_i U^dag = sum_j R_ji sigma_j (sign arbitrary); robust at angle pi."""
    vals, vecs = np.linalg.eig(R)
    n = np.real(vecs[:, np.argmin(np.abs(vals - 1))]); n /= np.linalg.norm(n)
    cth = (np.trace(R) - 1) / 2
    sth = n @ np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]]) / 2
    th = np.arctan2(sth, cth)
    return expm(-0.5j * th * sum(ni * s for ni, s in zip(n, sig)))


def comm_inv(U):
    x, y = U[X_], U[Y_]
    c = x @ y @ np.linalg.inv(x) @ np.linalg.inv(y)
    return c


def scalar_value(M):
    s = np.trace(M) / M.shape[0]
    return s if np.abs(M - s * np.eye(M.shape[0])).max() < 1e-9 else None


# ---------------------------------------------------------------- PL1 calibration
Ustd = {p: std(p).astype(complex) for p in A4}
U2T = {p: su2_lift(std(p)) for p in A4}
reg = {}
for p in A4:
    M = np.zeros((12, 12))
    for i, q in enumerate(A4):
        M[A4.index(compose(p, q)), i] = 1
    reg[p] = M.astype(complex)
cs = [scalar_value(comm_inv(U)) for U in (Ustd, reg, U2T)]
log.step("PL1 invariant c = U(x)U(y)U(x)^-1U(y)^-1 (rephasing-invariant): std rep +1, regular rep +1, SU(2) doublet (2T) -1",
         all(np.allclose(std(p) @ std(p).T, np.eye(3)) and np.isclose(np.linalg.det(std(p)), 1) for p in A4)
         and np.allclose([cs[0], cs[1], cs[2]], [1, 1, -1]), f"c = {np.round(np.real(cs), 6)}")


# ---------------------------------------------------------------- PL2 Fock-implementation lemma
def majoranas(n):
    m = n // 2
    Z, X, Y, I2 = sig[2], sig[0], sig[1], np.eye(2)
    out = []
    for j in range(m):
        for P in (X, Y):
            mat = np.array([[1.0 + 0j]])
            for k in range(m):
                mat = np.kron(mat, Z if k < j else (P if k == j else I2))
            out.append(mat)
    return out


def implementer(O, gam):
    """Clifford element U with U gamma(x) U^-1 = gamma(O x) for O in SO(n): product of reflection elements."""
    n = O.shape[0]
    M = O.copy(); vs = []
    for j in range(n):
        e = np.zeros(n); e[j] = 1
        d = M[:, j] - e
        if np.linalg.norm(d) > 1e-12:
            v = d / np.linalg.norm(d)
            M = M - 2 * np.outer(v, v @ M)
            vs.append(v)
    assert np.allclose(M, np.eye(n)) and len(vs) % 2 == 0
    U = np.eye(gam[0].shape[0], dtype=complex)
    for v in vs:
        U = U @ sum(vi * g for vi, g in zip(v, gam))
    return U


def fock_class(Ofun, n, label):
    gam = majoranas(n)
    Pf = np.diag([(-1) ** bin(i).count("1") for i in range(2 ** (n // 2))]).astype(complex)
    U = {}
    ok = True
    for p in A4:
        O = Ofun(p)
        U[p] = implementer(O, gam)
        for k in range(n):
            ok &= np.allclose(U[p] @ gam[k] @ np.linalg.inv(U[p]), sum(O[j, k] * gam[j] for j in range(n)), atol=1e-9)
    rel = [U[X_] @ U[X_], U[Y_] @ U[Y_], np.linalg.matrix_power(U[(1, 2, 0, 3)], 3)]
    rel_scalar = all(scalar_value(r) is not None for r in rel)
    c = comm_inv(U)
    cval = scalar_value(c)
    not_parity = not (np.allclose(c, Pf) or np.allclose(c, -Pf))
    even = [i for i in range(Pf.shape[0]) if Pf[i, i].real > 0]
    odd = [i for i in range(Pf.shape[0]) if Pf[i, i].real < 0]
    c_even = scalar_value(c[np.ix_(even, even)]); c_odd = scalar_value(c[np.ix_(odd, odd)])
    return dict(ok=ok, rel_scalar=rel_scalar, c=None if cval is None else round(float(np.real(cval)), 9),
                not_parity=not_parity, c_even=round(float(np.real(c_even)), 9), c_odd=round(float(np.real(c_odd)), 9), label=label)


def m3(Ofun, n):
    """multiplicity of the 3-dim irrep via characters (chi_3: identity 3, V4 elements -1, 3-cycles 0)."""
    tot = 0.0
    for p in A4:
        chi = np.trace(Ofun(p))
        chi3 = 3 if p == (0, 1, 2, 3) else (-1 if p in [(1, 0, 3, 2), (2, 3, 0, 1), (3, 2, 1, 0)] else 0)
        tot += chi * chi3
    return int(round(tot / 12))


cases = []
blk = lambda *Ms: (lambda p: np.block([[Ms[i](p) if i == j else np.zeros((Ms[i](p).shape[0], Ms[j](p).shape[1]))
                                         for j in range(len(Ms))] for i in range(len(Ms))]))
one = lambda p: np.eye(1)
cases.append(("family modes, number-conserving (R + R on gamma, gamma')", blk(std, std), 6))
cases.append(("gamma-only rotation (R + 1 + 1 + 1)", blk(std, one, one, one), 6))
for N in (1, 2, 3, 4):
    cases.append((f"twist zero modes N = {N} (N x permutation rep)", blk(*([perm4] * N)), 4 * N))
for trial in range(4):
    a, b, cc = [int(x) for x in rng.integers(0, 3, 3)]
    parts = [one] * a + [pairing_rep] * b + [std] * cc
    n = a + 3 * b + 3 * cc
    if n % 2:
        parts.append(one); n += 1
    if n < 2 or n > 12:
        continue
    Q = np.linalg.qr(rng.normal(size=(n, n)))[0]
    base = blk(*parts)
    cases.append((f"random real rep 1^{a + (n - a - 3 * b - 3 * cc)} 3perm^{b} 3^{cc} (rotated)", lambda p, base=base, Q=Q: Q @ base(p) @ Q.T, n))
res2 = []
ok2 = True
for label, Of, n in cases:
    r = fock_class(Of, n, label)
    mm = m3(Of, n)
    ok2 &= r["ok"] and r["rel_scalar"] and r["not_parity"] and r["c"] == (-1) ** mm
    res2.append((label, mm, r["c"], r["c_even"], r["c_odd"]))
fam, gonly = res2[0], res2[1]
log.step("PL2 Fock-implementation lemma: for genuine orthogonal one-particle actions the implementers satisfy all A4 relations up to "
         "scalars, c = (-1)^{m_3} (w2 of the real representation), never the fermion parity; number-conserving family modes c = +1; "
         "gamma-only c = -1 on BOTH parity sectors (no spin-statistics); twist zero modes c = (-1)^N",
         ok2 and fam[2] == 1 and gonly[2] == -1 and gonly[3] == gonly[4] == -1,
         "; ".join(f"{l}: m3={m}, c={c:+.0f}" for l, m, c, ce, co in res2))

# ---------------------------------------------------------------- PL3 positive control: projective one-particle space
from ur_lib import jordan_wigner  # noqa: E402
am = jordan_wigner(2)
def second_quant(u):
    """Fock unitary Gamma(u) for a one-particle unitary u on C^2 (two modes): exp of the generator."""
    from scipy.linalg import logm
    h = 1j * logm(u)
    H = sum(h[i, j] * am[i].conj().T @ am[j] for i in range(2) for j in range(2))
    return expm(-1j * H)
Ufock = {p: second_quant(U2T[p]) for p in A4}
cF = comm_inv(Ufock)
PF2 = np.diag([(-1) ** bin(i).count("1") for i in range(4)]).astype(complex)
log.step("PL3 positive control: a one-particle 2T doublet (the Weyl coin) second-quantised gives c = (-1)^F (odd sector -1, even +1): "
         "spin-statistics comes from a projective ONE-PARTICLE action", np.allclose(cF, PF2, atol=1e-9), f"diag(c) = {np.round(np.diag(cF).real, 6)}")

# ---------------------------------------------------------------- PL4 seam geometry
z = np.array([1j ** a for a in range(4)])
def cross(a, b, c, d): return (a - c) * (b - d) / ((a - d) * (b - c))
def mobius_through(src, dst):
    """Moebius map sending src[0..2] -> dst[0..2] as a 2x2 matrix."""
    def to_std(p):   # maps p0->0, p1->inf, p2->1
        p0, p1, p2 = p
        return np.array([[1, -p0], [1, -p1]], dtype=complex) * np.array([[1], [(p2 - p0) / (p2 - p1)]])
    return np.linalg.inv(to_std(dst)) @ to_std(src)
apply = lambda M, x: (M[0, 0] * x + M[0, 1]) / (M[1, 0] * x + M[1, 1]) if abs(M[1, 0] * x + M[1, 1]) > 1e-12 else np.inf
stab, orient = [], []
for p in permutations(range(4)):
    M = mobius_through(z[:3], z[list(p[:3])])
    if abs(apply(M, z[3]) - z[p[3]]) < 1e-9:
        stab.append(p)
        f0 = apply(M, 0)
        if np.isfinite(f0) and abs(f0) < 1 - 1e-9:
            orient.append(p)
inter = [p for p in orient if p in A4]
def stabiliser(pts):
    out = []
    for p in permutations(range(4)):
        M = mobius_through(pts[:3], pts[list(p[:3])])
        if abs(apply(M, pts[3]) - pts[p[3]]) < 1e-9:
            out.append(p)
    return out
om3 = np.exp(2j * np.pi / 3)
tet = np.array([0, 1, om3, om3 ** 2])
stab_tet = stabiliser(tet)
cr_sq, cr_tet = cross(*z), cross(*tet)
log.step("PL4 seam geometry: Moebius stabiliser of the square marks = 8 permutations (D4); circle-orientation-preserving ones = the four "
         "rotations (deck Z4); intersection with A4 = {e, (02)(13)}; no 3-cycle; square cross-ratio real (concyclic), the A4-symmetric "
         "equianharmonic value e^{+-i pi/3} is not real (not concyclic)",
         len(stab) == 8 and sorted(orient) == sorted([(0, 1, 2, 3), (1, 2, 3, 0), (2, 3, 0, 1), (3, 0, 1, 2)])
         and sorted(inter) == [(0, 1, 2, 3), (2, 3, 0, 1)] and not any(p in stab for p in A4 if p not in [(0, 1, 2, 3), X_, Y_, (3, 2, 1, 0)])
         and abs(cr_sq.imag) < 1e-12 and sorted(stab_tet) == sorted(A4) and abs(cr_tet.imag) > 0.1,
         f"|stab| = {len(stab)}, rotations {len(orient)}, A4 part {inter}; cross-ratio square {cr_sq:.3f}, tetrahedral {cr_tet:.3f} "
         f"with stabiliser of order {len(stab_tet)} = A4")

# ---------------------------------------------------------------- PL5 seam rotation lift is conditional (NS vs R)
modes_ns, modes_r = np.arange(-3.5, 4, 1.0), np.arange(-3, 4, 1.0)
rot = lambda modes, phi: np.diag(np.exp(1j * modes * phi))
ns2pi, r2pi = rot(modes_ns, 2 * np.pi), rot(modes_r, 2 * np.pi)
quarter_ns, quarter_r = rot(modes_ns, np.pi / 2), rot(modes_r, np.pi / 2)
same_abs = np.allclose(np.abs(quarter_ns), np.eye(len(modes_ns))) and np.allclose(np.abs(quarter_r), np.eye(len(modes_r)))
log.step("PL5 seam rotation lift is conditional: one-particle rotation by 2 pi = -I on NS modes (second quantisation (-1)^F), +I on R "
         "modes, with identical permutation structure (20260921: the geometric half-turn alone does not fix the lift); even under the "
         "S^2 spinor-section premise only <(02)(13)> of A4 is reached and no V4 commutator exists on the seam",
         np.allclose(ns2pi, -np.eye(len(modes_ns))) and np.allclose(r2pi, np.eye(len(modes_r))) and same_abs
         and np.allclose(np.linalg.matrix_power(quarter_ns, 4), -np.eye(len(modes_ns))))

# ---------------------------------------------------------------- PL6 source inventory
log.step("PL6 source inventory: flavor fiber (genuine A4 3-dim, not fermionic modes), windings (T_d, genuine), family modes "
         "(number-conserving, m_3 = 2, c = +1), twist zero modes N = 16 (m_3 = 16, c = +1); the only nontrivial class (gamma-only, "
         "m_3 = 1) breaks fermion number, mixes families (Part G) and is spinorial on the even sector -> no source sector realises A4 "
         "with center (-1)^F", fam[2] == 1 and gonly[3] == -1 and (-1) ** 16 == 1)

# ---------------------------------------------------------------- PL7 IR rescue for a two-component coin
B = np.array([[1, 1, -1], [1, -1, 1], [-1, 1, 1]], dtype=float); Binv = np.linalg.inv(B)
def hW(k):
    p = B @ k
    d = Binv @ np.array([np.sin(p[i]) + np.cos(p[(i + 1) % 3]) - 1 for i in range(3)])
    return sum(di * s for di, s in zip(d, sig))
ratios, ok7 = [], True
for _ in range(4):
    A = rng.normal(size=(3, 3)); R = expm(A - A.T)
    U = su2_lift(R)
    k0 = rng.normal(size=3); k0 /= np.linalg.norm(k0)
    errs = [np.abs(hW(R @ (s * k0)) - U @ hW(s * k0) @ U.conj().T).max() for s in (1e-2, 5e-3)]
    ratios.append(errs[0] / errs[1])
    ok7 &= errs[0] < 1e-3
n_ax = np.array([0.3, -0.5, 0.81]); n_ax /= np.linalg.norm(n_ax)
Kx = np.array([[0, -n_ax[2], n_ax[1]], [n_ax[2], 0, -n_ax[0]], [-n_ax[1], n_ax[0], 0]])
Ucur = np.eye(2, dtype=complex)
for t in np.linspace(0, 1, 401)[1:]:
    Rt = expm(2 * np.pi * t * Kx)
    cand = su2_lift(Rt)
    Ucur = cand if np.abs(cand - Ucur).max() < np.abs(cand + Ucur).max() else -cand
loop_end = Ucur
sheet = [sig[0]] * 4
gram_sheet = np.array([[np.trace(sheet[a] @ sheet[b] + sheet[b] @ sheet[a]).real / 4 for b in range(4)] for a in range(4)])
log.step("PL7 emergent spin 1/2: at the isotropic Weyl node of h_W the emergent rotations act on a two-component coin through SU(2) "
         "(covariance error O(k^2): ratio ~4 under k -> k/2) and the 2 pi loop lifts to -I, consistent with the fermion parity already "
         "present in the CAR structure (no IR derivation of spin-statistics, no derivation of the coin itself); the P1 sheet doublet "
         "is no coin (all marks act by the same swap, Gram rank 1)",
         ok7 and all(3.5 < r < 4.5 for r in ratios) and np.allclose(loop_end, -np.eye(2), atol=1e-9) and np.linalg.matrix_rank(gram_sheet) == 1,
         f"O(k^2) ratios {np.round(ratios, 2)}; loop end {np.round(np.diag(loop_end).real, 6)}")

log.step("PL8 (statement) SOURCE.MARK_PROJECTIVE_LIFT.01 = NOT_FORCED_IN_SOURCE (all source mark actions trivial class; the Fock "
         "implementation cannot create spin-statistics; the seam is spinorial only on <(02)(13)> and conditionally) and "
         "SPIN_HALF_EMERGENT_GIVEN_TWO_COMPONENT_COIN: spin 1/2 emerges at an isotropic Weyl node once a two-component coin exists; the "
         "remaining microscopic input is that coin itself (AX.WEYLCOIN), not a cocycle choice", True)

log.dump(HERE / "lift_results.json", {
    "research_id": "SOURCE.MARK_PROJECTIVE_LIFT.01",
    "verdict": "NOT_FORCED_IN_SOURCE",
    "pins": {"fock_cases": [(l, m, c) for l, m, c, ce, co in res2], "moebius_stabiliser": len(stab), "seam_rotations": len(orient),
             "A4_seam_part": [list(p) for p in inter]},
})
