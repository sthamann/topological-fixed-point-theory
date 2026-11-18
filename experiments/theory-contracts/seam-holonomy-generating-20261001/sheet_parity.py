#!/usr/bin/env python3
"""Part X: SHEET.PARITY.01 -- can the existing source supply the spacetime spin factor at all?

  X1 internal no-go: the real commutant of su(5) on the carrier fields R^10 = C^5 is C (dimension 2,
     abelian; for the Standard-Model subalgebra it is C + C, still abelian).  Any SU(2) acting linearly on the sixteen elementary fields and commuting with the carrier gauge
     algebra therefore acts trivially on the ten carrier fields (SU(2) has no nontrivial map to U(1)): carrier
     fermions cannot get spin 1/2 from an internal field-space action.  With Part G (family spinor kills) and
     Part K9 (16 sheet-odd species) no internal route remains.
  X2 geometric no-go: with L_n^dag = L_-n the Hermitian combinations of L_0, L_+-1 span su(1,1) (indefinite
     Killing form), not su(2): the unitary seam carries only the compact rotations e^{i theta L_0} (abelian;
     rho_F^4 = (-1)^F), no non-abelian geometric spin.
  X3 target correction for the 2O test: the Pin(4) lifts of the 24 mark relabelings (Part K) form a group of
     order 48 with 13 involutions and 12 elements of order 8 -- the GL(2,3)-type double cover of S4, not 2O
     (which has one involution); its even part is the binary tetrahedral group 2T.  The oriented seam realises
     T_d on the windings (std, det = sgn, -I never occurs), i.e. the rotations T only; the octahedral rotations
     O (whose cover is 2O) need the loop-orientation reversal.  2O is the cover of the FLAVOR rotation action.
  X4 walk symmetry and intermediate spaces: the mark walk is exactly covariant under the Klein four (axis
     pi-rotations); under the three-fold rotation it is only conjugate (by one partial factor) to the rotated
     walk, i.e. isospectral, not equal.  Each controlled partial shift of Part N maps the bcc coset of Z^3 to a
     different coset (index 4): Part N already runs on the union of the four intermediate spaces.
  X5 (statement) SHEET.PARITY.01: the spin factor is necessarily an EXTENSION of the present source (X1, X2,
     Part G, K9, Part P); P1 in its current formalisation does not force it.  Decision on AX.MARKMAJORANA (or
     any other spin carrier) is a choice, not a derivation.
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/sheet_parity.py
"""
from __future__ import annotations

import json
import sys
from itertools import permutations, product
from pathlib import Path

import numpy as np
from scipy.linalg import null_space

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


# X1 ----------------------------------------------------------------------------------------------
gens = []
for i in range(5):
    for j in range(i + 1, 5):
        E = np.zeros((5, 5), complex); E[i, j] = 1; E[j, i] = -1; gens.append(E)
        F = np.zeros((5, 5), complex); F[i, j] = 1j; F[j, i] = 1j; gens.append(F)
for i in range(4):
    D = np.zeros((5, 5), complex); D[i, i] = 1j; D[i + 1, i + 1] = -1j; gens.append(D)
realify = lambda A: np.block([[A.real, -A.imag], [A.imag, A.real]])
Rg = [realify(A) for A in gens]
basis = [np.eye(10)[:, [i]] @ np.eye(10)[[j], :] for i in range(10) for j in range(10)]
eqs = np.array([np.concatenate([(X @ G - G @ X).ravel() for G in Rg]) for X in basis]).T
comm = null_space(eqs)
Cs = [sum(comm[k, c] * basis[k] for k in range(100)) for c in range(comm.shape[1])]
abelian = all(np.allclose(A @ B, B @ A) for A in Cs for B in Cs)
def is_sm(A):
    blk = np.abs(A[:3, 3:]).max() + np.abs(A[3:, :3]).max()
    return blk < 1e-12
sm = [G for G in gens if is_sm(G)]                                       # su(3) + su(2) + the diagonal Cartans
sm_ok = []
for G in sm:
    t3, t2 = np.trace(G[:3, :3]), np.trace(G[3:, 3:])
    if abs(t3) < 1e-12 and abs(t2) < 1e-12:
        sm_ok.append(G)
sm_ok.append(1j * np.diag([-1 / 3, -1 / 3, -1 / 3, 0.5, 0.5]))                # hypercharge
eqs_sm = np.array([np.concatenate([(X @ realify(G) - realify(G) @ X).ravel() for G in sm_ok]) for X in basis]).T
comm_sm = null_space(eqs_sm)
Csm = [sum(comm_sm[k, c] * basis[k] for k in range(100)) for c in range(comm_sm.shape[1])]
abelian_sm = all(np.allclose(A @ B, B @ A) for A in Csm for B in Csm)
step("X1.internal no-go: the real commutant of the carrier gauge algebra on the carrier fields R^10 is abelian (su(5): C, dim 2; Standard-Model su(3)+su(2)+u(1): C + C, dim 4); an internal SU(2) commuting with it acts trivially on the carrier fields -> no spin 1/2 for carrier fermions from field-space actions",
     comm.shape[1] == 2 and abelian and len(gens) == 24 and comm_sm.shape[1] == 4 and abelian_sm and len(sm_ok) == 12,
     f"dim su(5) = {len(gens)}, commutant {comm.shape[1]}; dim sm = {len(sm_ok)}, commutant {comm_sm.shape[1]}")

# X2 ----------------------------------------------------------------------------------------------
def bracket(a, b):
    """Witt structure constants on span{L_-1, L_0, L_1}: [L_m, L_n] = (m - n) L_{m+n}; vectors indexed by n+1."""
    out = np.zeros(3, complex)
    for m in (-1, 0, 1):
        for n in (-1, 0, 1):
            if abs(m + n) <= 1:
                out[m + n + 1] += a[m + 1] * b[n + 1] * (m - n)
    return out
Lm1, L0, L1 = np.eye(3, dtype=complex)
Xh, Yh, Zh = (L1 + Lm1) / 2, 1j * (L1 - Lm1) / 2, L0                      # Hermitian for L_n^dag = L_-n
realb = [1j * Xh, 1j * Yh, 1j * Zh]                                      # anti-Hermitian generators of the unitary group
Bm = np.array(realb).T
def coords(v):
    c = np.linalg.solve(Bm, v)
    return c
ad = [np.array([coords(bracket(x, y)) for y in realb]).T for x in realb]
real_ad = all(np.allclose(A.imag, 0, atol=1e-12) for A in ad)
kill = np.array([[np.trace(A @ Bk).real for Bk in ad] for A in ad])
ev = np.linalg.eigvalsh(kill)
step("X2.geometric no-go: the Hermitian combinations of L_0, L_+-1 (with L_n^dag = L_-n) span su(1,1): Killing form indefinite, not negative definite; only e^{i theta L_0} is a compact unitary seam rotation (abelian, rho_F^4 = (-1)^F)",
     real_ad and (ev > 1e-9).any() and (ev < -1e-9).any(), f"Killing eigenvalues {np.round(ev, 6).tolist()}")

# X3 ----------------------------------------------------------------------------------------------
Zp = np.diag([1.0, -1.0]).astype(complex); I2 = np.eye(2, dtype=complex)
lower = np.array([[0, 1], [0, 0]], dtype=complex)
a0 = np.kron(lower, I2); a1 = np.kron(Zp, lower)
g = [a0 + a0.conj().T, -1j * (a0 - a0.conj().T), a1 + a1.conj().T, -1j * (a1 - a1.conj().T)]
I4 = np.eye(4)
w = np.array([[-.5, -.5, -.5], [-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]])


def pin_lift(pi):
    perm, out = list(range(4)), I4.astype(complex)
    for pos in range(4):
        if perm[pos] != pi[pos]:
            q = perm.index(pi[pos])
            out = (g[perm[pos]] - g[perm[q]]) / np.sqrt(2) @ out
            perm[pos], perm[q] = perm[q], perm[pos]
    return out


def std(pi):
    A = np.array([w[pi[a]] for a in range(3)]).T
    return A @ np.linalg.inv(np.array([w[a] for a in range(3)]).T)


gens_pin = [pin_lift(p) for p in permutations(range(4))]
key = lambda X: tuple(np.round(X, 8).ravel())
group = {key(I4): I4.astype(complex)}
frontier = [I4.astype(complex)]
while frontier:
    new = []
    for X in frontier:
        for G in gens_pin + [-I4]:
            Y = G @ X
            if key(Y) not in group:
                group[key(Y)] = Y; new.append(Y)
    frontier = new
def order(X):
    Y = X.copy()
    for n in range(1, 25):
        if np.allclose(Y, I4):
            return n
        Y = Y @ X
    return -1
orders = [order(X) for X in group.values()]
stats = {n: orders.count(n) for n in sorted(set(orders))}
Gam = -g[0] @ g[1] @ g[2] @ g[3]
even = [X for X in group.values() if np.allclose(X @ Gam, Gam @ X)]
even_stats = {n: [order(X) for X in even].count(n) for n in sorted({order(X) for X in even})}
stds = [std(p) for p in permutations(range(4))]
no_inversion = not any(np.allclose(R, -np.eye(3)) for R in stds)
dets = [round(np.linalg.det(R)) for R in stds]
rotations = sum(1 for d in dets if d == 1)
step("X3.target correction: the Pin(4) lifts of the 24 mark relabelings form a group of order 48 with 13 involutions and 12 elements of order 8 (GL(2,3)-type cover, not 2O with 1 involution); even part = 2T (order 24); the oriented seam realises T_d on the windings (12 rotations, -I absent): 2O is the cover of the flavor rotations O, not of the spatial action",
     len(group) == 48 and stats.get(2, 0) == 13 and stats.get(8, 0) == 12 and len(even) == 24 and even_stats == {1: 1, 2: 1, 3: 8, 4: 6, 6: 8}
     and no_inversion and rotations == 12, f"element orders {stats}; even part {even_stats}")

# X4 ----------------------------------------------------------------------------------------------
beta = sum(g) / 2
e = [sum(w[a, i] * g[a] for a in range(4)) for i in range(3)]
alpha = [1j * beta @ e[i] for i in range(3)]
Uf = lambda t, A_: np.cos(t) * I4 - 1j * np.sin(t) * A_
W = lambda k: Uf(k[0], alpha[0]) @ Uf(k[1], alpha[1]) @ Uf(k[2], alpha[2])
def adj(p):
    return np.array([[np.trace(alpha[j] @ p @ alpha[i] @ np.linalg.inv(p)).real / 4 for i in range(3)] for j in range(3)])
k = np.array([0.37, -0.81, 0.52])
klein_ok, c3_info = True, None
for pi in [(1, 0, 3, 2), (2, 3, 0, 1), (3, 2, 1, 0)]:
    p = pin_lift(pi); Rm = adj(p)
    klein_ok &= np.allclose(p @ W(k) @ np.linalg.inv(p), W(Rm @ k))
p3 = pin_lift((0, 2, 3, 1)); R3 = adj(p3)
lhs = p3 @ W(k) @ np.linalg.inv(p3)
rk = R3 @ k
eq = np.allclose(lhs, W(rk))
iso = np.allclose(np.sort(np.angle(np.linalg.eigvals(lhs))), np.sort(np.angle(np.linalg.eigvals(W(rk)))))
conj_ok = False
for i in range(3):
    P = Uf(rk[i], alpha[i])
    if np.allclose(lhs, np.linalg.inv(P) @ W(rk) @ P) or np.allclose(lhs, P @ W(rk) @ np.linalg.inv(P)):
        conj_ok = True
    for j in range(3):
        P2 = Uf(rk[i], alpha[i]) @ Uf(rk[j], alpha[j])
        if np.allclose(lhs, np.linalg.inv(P2) @ W(rk) @ P2) or np.allclose(lhs, P2 @ W(rk) @ np.linalg.inv(P2)):
            conj_ok = True
coset = lambda x: tuple(np.mod(np.array(x) - x[0], 2))
moves = all(coset(np.array(x) + np.eye(3, dtype=int)[i]) != coset(x) for x in [(0, 0, 0), (1, 1, 1), (2, 0, 4), (3, 1, -1)] for i in range(3))
step("X4.walk symmetry: exactly covariant under the Klein four (axis pi-rotations); under the 3-fold rotation only conjugate to the rotated walk by partial factors (isospectral, not equal); every controlled partial shift of Part N leaves the bcc coset of Z^3 (index 4): Part N runs on the union of the four intermediate spaces",
     klein_ok and not eq and iso and conj_ok and moves)
step("X5.(statement) SHEET.PARITY.01: no internal (X1, Part G, K9), geometric (X2) or projective-transport (Part P) route supplies the spin factor inside the present source; it is necessarily an extension, P1 in its current form does not force it",
     True, "decision required: AX.MARKMAJORANA or another explicitly added spin carrier")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "pin_group_orders": stats, "killing": ev.tolist()}
(HERE / "sheet_parity_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
