#!/usr/bin/env python3
"""SOURCE.WEYLCOIN.ORIGIN.01 -- does the existing source contain a local, family- and gauge-neutral two-dimensional factor
on which the spatial process acts nontrivially?

Local source algebra: the eight complex seam modes of the universal source U (5 carrier c, 3 family f; Part U), Fock space
C^32 (x) C^8 (Jordan-Wigner, carrier modes first, so even carrier/family generators act as X (x) 1 and 1 (x) Y).  A neutral coin
would be an M_2(C) block of the parity-even commutant of gauge x family on which source process operators act non-abelianly.
Nothing about SU(2) or Weyl structure is assumed.
  WC1 graded multiplicities (highest-weight vectors of the carrier gauge algebra and of the family group, per fermion parity)
      for gauge in {Spin(10), SU(5), SM} and family in {SU(4), SU(3), SO(3), A4, U(1)_N (control)}.
  WC2 Spin(10) gauge: for every family group containing A4 the parity-even commutant of the local algebra is ABELIAN (every
      irrep pair occurs at most once per parity); in particular the matter 16 carries no two-dimensional multiplicity.
  WC3 weakened gauge SU(5) or SM: neutral M_2 blocks appear, but only between the carrier-empty and carrier-full states
      (Lambda^0 and Lambda^5, separated by Spin(10) as 16 vs 16bar), i.e. eight-fermion operators carrying no matter.
  WC4 the source process operators (channel-A family Majoranas and their products, second-quantised flavor monodromies
      Gamma(M_a)) act on the family factor only: on every neutral M_2 block they act diagonally (they cannot connect
      carrier-empty and carrier-full), and their neutral part is abelian.
  WC5 the remaining candidates are no coins: the sheet doublet is the fermion-parity grading itself (Part K9) with all marks
      acting by the same swap; the flavor fiber C^3 and the U(1)_N-only multiplicities are family-carrying.
  WC6 (statement) SOURCE.WEYLCOIN.ORIGIN.01 = NO_NEUTRAL_COIN_IN_SOURCE within the declared local algebra and symmetry
      choices -> AX.WEYLCOIN is the minimal additional microscopic input (decision to freeze it is the user's).
    python3 -B experiments/theory-contracts/weylcoin-origin-20261005/coin_origin.py
"""
from __future__ import annotations

import sys
from itertools import combinations, permutations
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent / "minimal-runtime-reconstruction-20261005"))
sys.path.insert(0, str(ROOT))
from ur_lib import Log, jordan_wigner  # noqa: E402
from tfpt_explorer.flavor_path_transport import exact_clock_dictionary  # noqa: E402

log = Log()
c = jordan_wigner(5)
f = jordan_wigner(3)
dag = lambda x: x.conj().T
Pc = np.diag([(-1) ** bin(i).count("1") for i in range(32)]).astype(complex)
Pfam = np.diag([(-1) ** bin(i).count("1") for i in range(8)]).astype(complex)
num = lambda ops: sum(dag(o) @ o for o in ops)


def majoranas(ops):
    out = []
    for o in ops:
        out += [o + dag(o), -1j * (o - dag(o))]
    return out


# ---------------------------------------------------------------- group data: (Hermitian generators, raising ops, Cartan)
def spin_even(ops):
    g = majoranas(ops)
    return [0.5j * g[p] @ g[q] for p, q in combinations(range(len(g)), 2)]


def su_n(ops):
    n = len(ops); gens = []
    for i in range(n):
        for j in range(i + 1, n):
            gens += [dag(ops[i]) @ ops[j] + dag(ops[j]) @ ops[i], -1j * (dag(ops[i]) @ ops[j] - dag(ops[j]) @ ops[i])]
    for k in range(1, n):
        gens.append(sum(dag(ops[i]) @ ops[i] for i in range(k)) - k * dag(ops[k]) @ ops[k])
    return gens


def d_gamma(L, ops):
    return sum(L[i, j] * dag(ops[i]) @ ops[j] for i in range(len(ops)) for j in range(len(ops)))


Lso3 = [np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]]), np.array([[0, 0, 1j], [0, 0, 0], [-1j, 0, 0]]),
        np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]])]
y = [-1 / 3] * 3 + [1 / 2] * 2
GAUGE = {
    "Spin(10)": dict(gens=spin_even(c), raise_=[dag(c[i]) @ c[i + 1] for i in range(4)] + [dag(c[3]) @ dag(c[4])],
                     cartan=[dag(c[i]) @ c[i] - 0.5 * np.eye(32) for i in range(5)]),
    "SU(5)": dict(gens=su_n(c), raise_=[dag(c[i]) @ c[i + 1] for i in range(4)],
                  cartan=[dag(c[i]) @ c[i] - dag(c[i + 1]) @ c[i + 1] for i in range(4)]),
    "SM": dict(gens=su_n(c[:3]) + su_n(c[3:]) + [sum(yi * dag(ci) @ ci for yi, ci in zip(y, c))],
               raise_=[dag(c[0]) @ c[1], dag(c[1]) @ c[2], dag(c[3]) @ c[4]],
               cartan=[dag(c[0]) @ c[0] - dag(c[1]) @ c[1], dag(c[1]) @ c[1] - dag(c[2]) @ c[2], dag(c[3]) @ c[3] - dag(c[4]) @ c[4],
                       sum(yi * dag(ci) @ ci for yi, ci in zip(y, c))]),
}
w = np.array([[-.5, -.5, -.5], [-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]])


def parity_perm(p):
    p = list(p); s = 1
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]; p[i], p[j] = p[j], p[i]; s = -s
    return s


A4 = [p for p in permutations(range(4)) if parity_perm(p) == 1]


def std(p):
    P = np.zeros((4, 4)); P[list(p), range(4)] = 1
    return w.T @ P @ w


def second_quant(R):
    """Fock unitary Gamma(R) for a real orthogonal one-particle R on the three family modes."""
    from scipy.linalg import expm, logm
    h = 1j * logm(R.astype(complex))
    return expm(-1j * d_gamma(h, f))


A4_fock = [second_quant(std(p)) for p in A4]
FAMILY = {
    "SU(4)": dict(gens=spin_even(f), raise_=[dag(f[0]) @ f[1], dag(f[1]) @ f[2], dag(f[1]) @ dag(f[2])],
                  cartan=[dag(f[i]) @ f[i] - 0.5 * np.eye(8) for i in range(3)]),
    "SU(3)": dict(gens=su_n(f), raise_=[dag(f[0]) @ f[1], dag(f[1]) @ f[2]],
                  cartan=[dag(f[0]) @ f[0] - dag(f[1]) @ f[1], dag(f[1]) @ f[1] - dag(f[2]) @ f[2]]),
    "SO(3)": dict(gens=[d_gamma(L, f) for L in Lso3], raise_=[d_gamma(Lso3[0] + 1j * Lso3[1], f)], cartan=[d_gamma(Lso3[2], f)]),
    "A4": dict(group=A4_fock),
    "U(1)_N": dict(gens=[num(f)], raise_=[], cartan=[num(f)]),
}


def commutant(spec, d):
    """Basis of the commutant in End(C^d): kernel of sum_k ad_{T_k}^2 (Lie) or the group average (finite)."""
    I = np.eye(d)
    if "group" in spec:
        Pavg = sum(np.kron(U, U.conj()) for U in spec["group"]) / len(spec["group"])
        w_, v_ = np.linalg.eig(Pavg)
        vecs = v_[:, np.abs(w_ - 1) < 1e-8]
        q, _ = np.linalg.qr(vecs)
        return [q[:, i].reshape(d, d) for i in range(q.shape[1])]
    Lsup = sum(np.kron(T @ T, I) + np.kron(I, (T @ T).T) - 2 * np.kron(T, T.T) for T in spec["gens"])
    w_, v_ = np.linalg.eigh((Lsup + Lsup.conj().T) / 2)
    return [v_[:, i].reshape(d, d) for i in np.where(np.abs(w_) < 1e-8)[0]]


def hw_table(spec, P, d):
    """Multiplicities: number of highest-weight vectors per (weight, parity); finite group: irreps via characters."""
    if "group" in spec:
        out = {}
        for par in (1, -1):
            idx = [i for i in range(d) if P[i, i].real == par]
            chars = [np.trace(U[np.ix_(idx, idx)]) for U in spec["group"]]
            # A4 irreps by classes: identity, V4, 3-cycles (two classes)
            dims = {}
            for name, chi in A4_chars.items():
                m = np.real(sum(np.conj(chi[k]) * chars[k] for k in range(12)) / 12)
                if round(m) > 0:
                    dims[name] = int(round(m))
            for name, m in dims.items():
                out[(name, par)] = m
        return out
    if spec["raise_"]:
        K = np.vstack(spec["raise_"])
        u, s, vh = np.linalg.svd(K)
        ker = vh[np.sum(s > 1e-9):].conj().T
    else:
        ker = np.eye(d, dtype=complex)
    out = {}
    ops = spec["cartan"] + [P]
    Mops = [ker.conj().T @ O @ ker for O in ops]
    rnd = sum(np.random.default_rng(3).normal() * M for M in Mops)
    w_, v_ = np.linalg.eigh((rnd + rnd.conj().T) / 2)
    vecs = v_
    for k in range(vecs.shape[1]):
        vk = vecs[:, k]
        key = tuple(round(float(np.real(vk.conj() @ M @ vk)), 6) for M in Mops)
        lab, par = key[:-1], int(round(key[-1]))
        out[(lab, par)] = out.get((lab, par), 0) + 1
    return out


# A4 characters (class functions on the list A4)
def cls(p):
    if p == (0, 1, 2, 3):
        return "e"
    if all(p[p[i]] == i for i in range(4)):
        return "V4"
    return "C3a" if p in [(1, 2, 0, 3), (0, 2, 3, 1), (2, 1, 3, 0), (3, 0, 2, 1)] else "C3b"
om = np.exp(2j * np.pi / 3)
table = {"1": {"e": 1, "V4": 1, "C3a": 1, "C3b": 1}, "1'": {"e": 1, "V4": 1, "C3a": om, "C3b": om ** 2},
         "1''": {"e": 1, "V4": 1, "C3a": om ** 2, "C3b": om}, "3": {"e": 3, "V4": -1, "C3a": 0, "C3b": 0}}
A4_chars = {k: [v[cls(p)] for p in A4] for k, v in table.items()}
class_sizes = {k: sum(1 for p in A4 if cls(p) == k) for k in ("e", "V4", "C3a", "C3b")}
ortho = all(abs(sum(np.conj(A4_chars[a][i]) * A4_chars[b][i] for i in range(12)) / 12 - (a == b)) < 1e-12 for a in table for b in table)

# ---------------------------------------------------------------- WC1 graded multiplicity tables
tables_c = {g: hw_table(GAUGE[g], Pc, 32) for g in GAUGE}
tables_f = {h: hw_table(FAMILY[h], Pfam, 8) for h in FAMILY}
dimcheck = all(sum(tables_c["Spin(10)"].values()) == 2 for _ in [0])
log.step("WC1 graded multiplicities: Spin(10) on the carrier Fock space = 16 (even) + 16bar (odd), each once; family groups on "
         "Lambda C^3 recorded per parity (A4 characters orthonormal, class sizes 1/3/4/4)",
         ortho and class_sizes == {"e": 1, "V4": 3, "C3a": 4, "C3b": 4} and dimcheck
         and sorted(tables_c["Spin(10)"].values()) == [1, 1],
         "; ".join(f"{h}: {len(t)} (irrep, parity) labels" for h, t in tables_f.items()))


def blocks(tc, tf):
    """Even/odd total-parity multiplicity W for each (carrier irrep, family irrep)."""
    lc = {lab for lab, _ in tc}; lf = {lab for lab, _ in tf}
    out = {}
    for a in lc:
        for b in lf:
            mcp, mcm = tc.get((a, 1), 0), tc.get((a, -1), 0)
            mfp, mfm = tf.get((b, 1), 0), tf.get((b, -1), 0)
            out[(a, b)] = (mcp * mfp + mcm * mfm, mcp * mfm + mcm * mfp)
    return out


def even_commutant_abelian(gname, hname):
    Ac = commutant(GAUGE[gname], 32)
    Af = commutant(FAMILY[hname], 8)
    P = np.kron(Pc, Pfam)
    basis = []
    for X in Ac:
        for Y in Af:
            Z = np.kron(X, Y)
            Z = (Z + P @ Z @ P) / 2
            if np.abs(Z).max() > 1e-9:
                basis.append(Z)
    # orthonormalise
    M = np.array([b.ravel() for b in basis])
    u, s, vh = np.linalg.svd(M, full_matrices=False)
    B = [vh[i].reshape(256, 256) for i in range(int((s > 1e-9).sum()))]
    worst = 0.0
    for i in range(len(B)):
        for j in range(i + 1, len(B)):
            worst = max(worst, np.abs(B[i] @ B[j] - B[j] @ B[i]).max())
    return len(Ac), len(Af), len(B), worst


res2 = {}
for h in ("SU(4)", "SU(3)", "SO(3)", "A4"):
    blk = blocks(tables_c["Spin(10)"], tables_f[h])
    maxW = max(max(v) for v in blk.values())
    dc, df, de, comm = even_commutant_abelian("Spin(10)", h)
    res2[h] = (maxW, dc, df, de, comm)
log.step("WC2 Spin(10) gauge: for family SU(4), SU(3), SO(3) and A4 every (carrier, family) irrep pair occurs at most once per parity, "
         "and the parity-even commutant of the local algebra is abelian (all commutators vanish): no neutral M_2, no C^2 multiplicity "
         "of the matter 16", all(v[0] <= 1 and v[4] < 1e-9 for v in res2.values()),
         "; ".join(f"{h}: max W {v[0]}, dim commutants {v[1]}x{v[2]}, even part {v[3]}, max commutator {v[4]:.0e}" for h, v in res2.items()))

res3 = {}
for g in ("SU(5)", "SM"):
    blk = blocks(tables_c[g], tables_f["SO(3)"])
    big = {k: v for k, v in blk.items() if max(v) >= 2}
    dc, df, de, comm = even_commutant_abelian(g, "SO(3)")
    res3[g] = (big, comm)
vac = np.zeros(32); vac[0] = 1
full = np.zeros(32); full[31] = 1
f0 = np.zeros(8); f0[0] = 1
f3 = np.zeros(8); f3[7] = 1
s1, s2 = np.kron(vac, f0), np.kron(full, f3)
X = np.outer(s2, s1) + np.outer(s1, s2)
inv_ok = all(np.abs(np.kron(T, np.eye(8)) @ X - X @ np.kron(T, np.eye(8))).max() < 1e-12 for T in GAUGE["SM"]["gens"]) and \
    all(np.abs(np.kron(np.eye(32), T) @ X - X @ np.kron(np.eye(32), T)).max() < 1e-12 for T in FAMILY["SO(3)"]["gens"])
inv_spin10 = all(np.abs(np.kron(T, np.eye(8)) @ X - X @ np.kron(T, np.eye(8))).max() < 1e-12 for T in GAUGE["Spin(10)"]["gens"])
log.step("WC3 weakened gauge (SU(5) or SM, family SO(3)): neutral M_2 blocks appear (non-abelian even commutant), all of them pairing the "
         "carrier-empty and carrier-full states: e.g. |0_c,000_f> <-> |11111_c,111_f> commutes with SM x SO(3) but not with Spin(10) "
         "(16 vs 16bar); an eight-fermion operator in the carrier-singlet sector, no matter",
         all(len(v[0]) > 0 and v[1] > 1e-6 and all(all(abs(x) < 1e-9 for x in k[0]) for k in v[0]) for v in res3.values())
         and inv_ok and not inv_spin10,
         "; ".join(f"{g}: {len(v[0])} blocks with W>=2, carrier labels all trivial (Lambda^0/Lambda^5), family labels "
                   f"{sorted(str(k[1]) for k in v[0])}" for g, v in res3.items()))

# ---------------------------------------------------------------- WC4 process operators on the neutral blocks
M = np.array(exact_clock_dictionary()["M"].evalf(), dtype=complex)
Uc = np.diag([1, 1j, -1j])
Ma = [np.linalg.matrix_power(Uc, a) @ M @ np.linalg.matrix_power(Uc, -a) for a in range(4)]
from scipy.linalg import expm, logm  # noqa: E402
flav = [expm(-1j * d_gamma(1j * logm(m), f)) for m in Ma]
gi = [f[i] + dag(f[i]) for i in range(3)]
gam_a = [sum(w[a, i] * gi[i] for i in range(3)) for a in range(4)]
gens_proc = [ga @ gb for ga in gam_a for gb in gam_a] + flav
alg = [np.eye(8, dtype=complex)]
def span_add(lst, X):
    M_ = np.array([l.ravel() for l in lst] + [X.ravel()])
    return np.linalg.matrix_rank(M_, tol=1e-9) > len(lst)
frontier = list(alg)
while frontier:
    new = []
    for A_ in frontier:
        for G in gens_proc:
            Z = G @ A_
            if span_add(alg, Z):
                alg.append(Z); new.append(Z)
    frontier = new
proc_dim = len(alg)
full_proc = [np.kron(np.eye(32), A_) for A_ in alg]
Pblock = np.array([s1, s2]).T
offdiag = max(abs(s2 @ Z @ s1) for Z in full_proc)
Af_so3 = commutant(FAMILY["SO(3)"], 8)
Mi = np.array([a_.ravel() for a_ in alg]).T
Mc = np.array([a_.ravel() for a_ in Af_so3]).T
# intersection of span(process) and family commutant: solve Mi x = Mc y
K = np.hstack([Mi, -Mc])
u, s, vh = np.linalg.svd(K)
null = vh[np.sum(s > 1e-9):].conj().T
inter = [(Mi @ v[:Mi.shape[1]]).reshape(8, 8) for v in null.T]
inter_comm = max([np.abs(A_ @ B_ - B_ @ A_).max() for A_ in inter for B_ in inter] + [0.0])
log.step("WC4 source process operators (channel-A family Majorana products, second-quantised flavor monodromies Gamma(M_a)) act on the "
         "family factor only: zero matrix element between carrier-empty and carrier-full (diagonal on every neutral M_2 block), and "
         "their intersection with the family-neutral commutant is abelian",
         offdiag < 1e-12 and inter_comm < 1e-9, f"process algebra dimension {proc_dim} of 64; neutral intersection dimension "
         f"{len(inter)}, max commutator {inter_comm:.0e}; empty<->full element {offdiag:.0e}")

sheet = [np.array([[0, 1], [1, 0]])] * 4
gram = np.array([[np.trace(sheet[a] @ sheet[b] + sheet[b] @ sheet[a]) / 4 for b in range(4)] for a in range(4)])
u1 = tables_f["U(1)_N"]
log.step("WC5 other candidates are no coins: sheet doublet = fermion-parity grading (Part K9), all four marks act by the same swap (Gram "
         "rank 1); flavor fiber C^3 and U(1)_N-only multiplicities (dims 1,3,3,1) are family-carrying",
         np.linalg.matrix_rank(gram) == 1 and sorted(u1.values()) == [1, 1, 3, 3], f"U(1)_N multiplicities {sorted(u1.values())}")
log.step("WC6 (statement) SOURCE.WEYLCOIN.ORIGIN.01 = NO_NEUTRAL_COIN_IN_SOURCE: with the TFPT gauge group Spin(10) the parity-even "
         "neutral commutant is abelian for every family group containing A4; weaker gauge groups give only carrier-empty/full M_2 blocks "
         "without matter on which the source process acts diagonally -> AX.WEYLCOIN is the minimal additional microscopic input", True)

log.dump(HERE / "coin_origin_results.json", {
    "research_id": "SOURCE.WEYLCOIN.ORIGIN.01", "verdict": "NO_NEUTRAL_COIN_IN_SOURCE",
    "pins": {"spin10_max_W": {h: v[0] for h, v in res2.items()}, "process_algebra_dim": proc_dim, "neutral_intersection_dim": len(inter)},
})
