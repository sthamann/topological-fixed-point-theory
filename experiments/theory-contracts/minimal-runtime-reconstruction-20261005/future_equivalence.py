#!/usr/bin/env python3
"""UR.FUTURE_EQUIVALENCE.FINITERANK.01 -- future indistinguishability as a closure problem, not a prediction problem.

Definition: two histories h, h' are equivalent iff every admissible continuation gives the same answers.
The checks show that this infinite condition is decided by a FINITE algebraic closure:
  FE1 button example: two states with equal present answers are separated exactly when an admissible
      operation X exists that maps their difference into the measured effects.
  FE2 Heisenberg closure O_{n+1} = O_n + sum_a Phi_a^dag(O_n) terminates; the dimension rises strictly until it
      stops, so there are at most d^2 - 1 enlargements (d = Hilbert dimension).
  FE3 soundness and completeness against brute force: a state difference orthogonal to O_* gives identical
      outcome probabilities for EVERY word up to length 6 (never used in the construction); a difference not
      orthogonal is separated by a word whose length is at most the closure index.
  FE4 history kernel K(u,v) = omega(C_u^dag C_v): ranks r_L from kernel values only; r_{L+1} = r_L implies
      closure for all longer words; rank <= d * rank(rho); the realization G_a = K_BB^-1 K_{B,aB} reproduces
      kernel entries of longer, unused words; adjoint letters are realized as K-adjoints.
  FE5 relations under the quotient: pi_omega is a unital *-homomorphism, so EVERY relation of the source operations
      survives (none is removed); additional relations come only from ker pi_omega = {x : omega(y^dag x z) = 0}, i.e.
      when the state lives on a sector where they already hold; with a faithful state the quotient is isomorphic to the
      generated algebra.  Premise: the kernel has the form omega(C_u^dag C_v) with single-operator letters (then the null
      space is a left ideal and the G_a are well defined); for multi-Kraus instruments the Heisenberg closure (FE2/FE3)
      is the relevant equivalence, not the GNS quotient.
  FE6 survival examples used later: the A4 relations survive on every quotient and the V4 commutator never becomes -I
      (an A4-specific fact, not a general class argument: genuine representations of Q8 or 2^(1+4) do anticommute);
      commuting generators stay commuting, and commuting operators with scalar anticommutators 2 g_ab, g > 0, do not
      exist for two or more generators.
  FE7 two equivalences: tests restricted to an IR window (low-energy projector) give a coarser quotient: two
      Hamiltonians that agree on the window are IR-equivalent but fundamentally distinguishable.
    python3 -B experiments/theory-contracts/minimal-runtime-reconstruction-20261005/future_equivalence.py
"""
from __future__ import annotations

import sys
from itertools import product
from pathlib import Path

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ur_lib import (Log, all_words, amplify, apply_channel, gns_quotient, kernel_from_ops, observability_closure,  # noqa: E402
                    purify, realized_kernel, word_op)

rng = np.random.default_rng(20261005)
log = Log()


def rand_unitary(d):
    z = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    q, r = np.linalg.qr(z)
    return q * (np.diag(r) / abs(np.diag(r)))


def rand_kraus(d, n=2):
    ks = [rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)) for _ in range(n)]
    s = sum(k.conj().T @ k for k in ks)
    w, v = np.linalg.eigh(s)
    isq = v @ np.diag(w ** -0.5) @ v.conj().T
    return [k @ isq for k in ks]


def rand_state(d, rank=None):
    rank = rank or d
    a = rng.normal(size=(d, rank)) + 1j * rng.normal(size=(d, rank))
    r = a @ a.conj().T
    return r / np.trace(r)


def word_probs(kraus_sets, effects, rho, L):
    out = []
    for w in all_words(len(kraus_sets), L):
        r = rho
        for a in w:
            r = apply_channel(kraus_sets[a], r)
        out.extend(np.trace(e @ r).real for e in effects)
    return np.array(out)


# ---------------------------------------------------------------- FE1 button example
I2 = np.eye(2, dtype=complex)
Z = np.diag([1, -1]).astype(complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Hd = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
P0 = np.diag([1, 0]).astype(complex)
plus = np.array([[1, 1], [1, 1]]) / 2
minus = np.array([[1, -1], [-1, 1]]) / 2
Rz = expm(-0.7j * Z)
effects = [I2, P0]
O_noX, d_noX = observability_closure([[Rz]], effects)
O_X, d_X = observability_closure([[Rz], [Hd]], effects)
sep = lambda O, d: max(abs(np.vdot(b.reshape(d, d).conj().T.ravel(), (plus - minus).ravel())) for b in O)
p_noX = word_probs([[Rz]], [P0], plus, 5) - word_probs([[Rz]], [P0], minus, 5)
p_X = word_probs([[Rz], [Hd]], [P0], plus, 2) - word_probs([[Rz], [Hd]], [P0], minus, 2)
log.step("FE1 button: |+> and |-> agree on the present test P0; with only Rz they stay equivalent for all words (closure span{I,Z}); "
         "the admissible Hadamard X separates them (closure grows to M_2)",
         abs(np.trace(P0 @ (plus - minus))) < 1e-15 and d_noX[-1] == 2 and sep(O_noX, 2) < 1e-12 and abs(p_noX).max() < 1e-12
         and d_X[-1] == 4 and sep(O_X, 2) > 0.5 and abs(p_X).max() > 0.9,
         f"closure dims without X {d_noX}, with X {d_X}; max prob difference {abs(p_noX).max():.1e} vs {abs(p_X).max():.2f}")

# ---------------------------------------------------------------- FE2 termination bound
ok, hist = True, []
for d in (2, 3, 4, 5):
    for trial in range(4):
        ks = [rand_kraus(d) for _ in range(2)]
        eff = [np.eye(d, dtype=complex), np.diag(rng.random(d)).astype(complex)]
        O, dims = observability_closure(ks, eff)
        strictly = all(dims[i + 1] > dims[i] for i in range(len(dims) - 2))
        ok &= strictly and dims[-1] <= d * d and (len(dims) - 2) <= d * d - 1
        hist.append((d, dims))
log.step("FE2 closure terminates: dimension rises strictly until the first repetition; <= d^2 - 1 enlargements",
         ok, "; ".join(f"d={d}:{dims}" for d, dims in hist[::4]))

# ---------------------------------------------------------------- FE3 brute-force soundness / completeness
# structured instrument: operations and effects act on the first qubit of C^2 (x) C^2 -> O_* = M_2 (x) I
d = 4
loc = lambda m: np.kron(m, I2)
ks = [[loc(rand_unitary(2))], [loc(k) for k in rand_kraus(2)]]
eff = [np.eye(4, dtype=complex), loc(P0)]
O, dims = observability_closure(ks, eff)
herm = lambda a: (a + a.conj().T) / 2
D = herm(np.kron(rng.normal(size=(2, 2)), rng.normal(size=(2, 2)) * 1j + rng.normal(size=(2, 2))))
proj = sum(np.vdot(b, D.ravel()) * b for b in O).reshape(4, 4)
Delta = herm(D - proj)
Delta -= np.trace(Delta) / 4 * np.eye(4)
rho1 = rand_state(4)
rho2 = rho1 + 0.5 * np.linalg.eigvalsh(rho1).min() * Delta / np.linalg.norm(Delta, 2)
pos = np.linalg.eigvalsh(rho2).min() > 0
dev_equiv = abs(word_probs(ks, eff, rho1, 6) - word_probs(ks, eff, rho2, 6)).max()
O0 = observability_closure([], eff)[0]
D3 = sum(np.vdot(b, D.ravel()) * b for b in O).reshape(4, 4)                 # inside O_* ...
D3 = herm(D3 - sum(np.vdot(b, D3.ravel()) * b for b in O0).reshape(4, 4))    # ... but invisible to the present effects
rho3 = rho1 + 0.5 * np.linalg.eigvalsh(rho1).min() * D3 / np.linalg.norm(D3, 2)
shortest = None
for L in range(0, 6):
    if abs(word_probs(ks, eff, rho1, L) - word_probs(ks, eff, rho3, L)).max() > 1e-8:
        shortest = L
        break
orth = max(abs(np.vdot(b, Delta.ravel())) for b in O)
log.step("FE3 brute force: difference orthogonal to O_* (dim 4 = M_2 (x) I) -> identical probabilities for all 127 words up to "
         "length 6; a difference inside O_* but invisible now is separated at length 1 <= closure index",
         dims[-1] == 4 and pos and orth < 1e-12 and dev_equiv < 1e-13 and shortest == 1 and shortest <= len(dims) - 2,
         f"closure dims {dims}; equivalent pair max deviation {dev_equiv:.1e}; present-invisible pair separated at length {shortest}")

# ---------------------------------------------------------------- FE4 history kernel, ranks, realization
ok4, info = True, []
for d, rk, nl in ((3, 1, 2), (3, 3, 2), (4, 2, 3)):
    ops = [rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)) for _ in range(nl)]
    ops = [o / np.linalg.norm(o, 2) for o in ops]
    alph = ops + [o.conj().T for o in ops]
    rho = rand_state(d, rk)
    A = amplify(alph)
    om = purify(rho)
    kern = kernel_from_ops(A, om)
    q = gns_quotient(kern, len(A), 2 * d * d)
    r = len(q["basis"])
    stable = q["closed_at"] is not None
    test = [w for w in all_words(len(A), 4) if len(w) >= 3][:: 7]
    err = max(abs(realized_kernel(q, u, v) - kern(u, v)) for u in test[:40] for v in test[:40])
    Kbb = q["gram"]
    adj = max(np.abs(q["G"][a + nl] - np.linalg.solve(Kbb, q["G"][a].conj().T @ Kbb)).max() for a in range(nl))
    ok4 &= stable and r <= d * rk and err < 1e-8 and adj < 1e-8
    info.append(f"d={d}, rank(rho)={rk}: r_L={q['ranks']}, err {err:.1e}, adj {adj:.1e}")
log.step("FE4 kernel-only quotient: rank stabilises (r_{L+1} = r_L), r <= d rank(rho), realization reproduces unused longer "
         "words, adjoint letters = K-adjoints", ok4, "; ".join(info))

# ---------------------------------------------------------------- FE5 no new relations
sx, sy, szm = X, np.array([[0, -1j], [1j, 0]]), Z
A_ops = [sx, sy, szm]                                            # Clifford on the first sector
B_ops = [herm(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))) for _ in range(3)]   # generic second sector
C = [np.block([[a, np.zeros((2, 2))], [np.zeros((2, 2)), b]]) for a, b in zip(A_ops, B_ops)]
def clifford_defect(G):
    n = len(G)
    dd = G[0].shape[0]
    worst = 0.0
    for a in range(n):
        for b in range(n):
            s = G[a] @ G[b] + G[b] @ G[a]
            worst = max(worst, np.abs(s - np.trace(s) / dd * np.eye(dd)).max())
    return worst
res = {}
for name, rho in (("support in Clifford sector", np.diag([0.5, 0.5, 0, 0]).astype(complex)),
                  ("faithful", np.diag([0.4, 0.3, 0.2, 0.1]).astype(complex))):
    q = gns_quotient(kernel_from_ops(amplify(C), purify(rho)), 3, 8)
    res[name] = (len(q["basis"]), clifford_defect(q["G"]))
log.step("FE5 quotient relations = source relations + ker pi_omega: {C_a, C_b} = 2 delta_ab appears on the quotient exactly when the "
         "state lives on the sector where it already holds; with a faithful state the generic sector survives and it fails",
         res["support in Clifford sector"][1] < 1e-9 and res["faithful"][1] > 1e-2,
         ", ".join(f"{k}: rank {v[0]}, defect {v[1]:.1e}" for k, v in res.items()))

# ---------------------------------------------------------------- FE6 class preservation, commuting -> no Clifford
a4_s = np.diag([1, -1, -1]).astype(complex)                      # V4 element
a4_t = np.roll(np.eye(3), 1, axis=0).astype(complex)             # order-3 cyclic permutation
gens = [a4_s, a4_t]
ok6 = True
for rho in (np.diag([1, 0, 0]).astype(complex), np.eye(3) / 3):
    q = gns_quotient(kernel_from_ops(amplify(gens + [g.conj().T for g in gens]), purify(rho)), 4, 12)
    G = q["G"]
    r = len(q["basis"])
    Ir = np.eye(r)
    rel = (np.abs(G[0] @ G[0] - Ir).max() < 1e-9 and np.abs(np.linalg.matrix_power(G[1], 3) - Ir).max() < 1e-9
           and np.abs(np.linalg.matrix_power(G[0] @ G[1], 3) - Ir).max() < 1e-9)
    comm = G[0] @ G[1] @ G[0] @ np.linalg.inv(G[1])              # s t s t^-1 lies in V4, never -I
    ok6 &= rel and np.abs(comm + Ir).max() > 0.5
# commuting family: best scalar fit of anticommutators leaves a residual whenever g would be positive definite
Uu = rand_unitary(6)
comm_ops = [Uu @ np.diag(rng.choice([-1, 1], 6)).astype(complex) @ Uu.conj().T for _ in range(4)]
ok6 &= clifford_defect(comm_ops) > 0.5
log.step("FE6 survival: A4 relations (s^2 = t^3 = (st)^3 = 1) hold on every quotient and the V4 commutator never becomes -I "
         "(A4-specific); commuting involutions never acquire scalar anticommutators", ok6,
         "pure and tracial states; commuting witness defect %.2f" % clifford_defect(comm_ops))

# ---------------------------------------------------------------- FE7 IR versus fundamental equivalence
d = 6
H_low = herm(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))
H1 = np.zeros((d, d), dtype=complex); H1[:2, :2] = H_low; H1[2:, 2:] = np.diag([5.0, 6.0, 7.0, 8.0])
H2 = H1.copy(); H2[2:, 2:] = np.diag([5.0, 6.5, 9.0, 8.0])
Pir = np.zeros((d, d), dtype=complex); Pir[:2, :2] = np.eye(2)
rho0 = np.zeros((d, d), dtype=complex); rho0[0, 0] = 1
U1, U2 = expm(-1j * 0.3 * H1), expm(-1j * 0.3 * H2)
ir_eff = [Pir @ e @ Pir for e in (np.eye(d), np.diag([1, 0, 0, 0, 0, 0]).astype(complex))]
dev_ir = max(abs(np.trace(e @ apply_channel([Pir @ U @ Pir], rho0)) - np.trace(e @ apply_channel([Pir @ V @ Pir], rho0))).real
             for e in ir_eff for U, V in ((U1, U2), (U1 @ U1, U2 @ U2)))
psi = np.zeros(d, dtype=complex); psi[1] = 1 / np.sqrt(2); psi[4] = 1 / np.sqrt(2)
phase_test = lambda U: abs(np.vdot(psi, U @ U @ psi)) ** 2
log.step("FE7 IR versus fundamental: H1, H2 agree on the IR window -> identical IR-restricted answers; an admissible "
         "high-energy intervention (window-crossing superposition) separates them",
         dev_ir < 1e-13 and abs(phase_test(U1) - phase_test(U2)) > 1e-2,
         f"IR deviation {dev_ir:.1e}; fundamental return-probability difference {abs(phase_test(U1) - phase_test(U2)):.3f}")

log.dump(HERE / "future_equivalence_results.json", {
    "research_id": "UR.FUTURE_EQUIVALENCE.FINITERANK.01",
    "verdict": "PASS_MATHEMATICAL",
    "statement": "Future indistinguishability is decided by a finite Heisenberg/Hankel closure (<= d^2 - 1 enlargements); "
                 "the quotient is the GNS homomorphism pi_omega: it removes no source relation and adds only those in "
                 "ker pi_omega (none for faithful states). Class obstructions must be shown channel by channel.",
})
