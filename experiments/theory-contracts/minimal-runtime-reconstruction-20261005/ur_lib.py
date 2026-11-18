"""Shared tools for UR.MINIMAL_RUNTIME.RECONSTRUCTION.01.

Everything here works from KERNEL VALUES only, the way the Universalraum principle demands:
  * gns_quotient: the history kernel K(u, v) = omega(C_u^dag C_v) on words u, v; a basis of words is grown
    letter by letter, rank r_L = dim span{C_u Omega : |u| <= L}; once r_{L+1} = r_L the span is closed under
    every letter, hence under every longer word (no future oracle).  The induced operators
    G_a = K_BB^{-1} K_{B, aB} are computed from kernel entries alone and are checked against kernel entries
    of longer words that were never used in the construction.
  * observability_closure: Heisenberg-picture closure O_{n+1} = O_n + sum_a Phi_a^dag(O_n) of the effect
    space for CP instruments; it stops after at most d^2 - 1 enlargements.
  * Majorana monomials as (phase, bitmask) with exact signs, for the twist-zero-mode channel.
"""
from __future__ import annotations

import json
from itertools import product
from pathlib import Path

import numpy as np

TOL = 1e-10


class Log:
    def __init__(self) -> None:
        self.rows: list[dict] = []

    def step(self, name: str, ok: bool, detail: str = "") -> bool:
        self.rows.append({"step": name, "pass": bool(ok), "detail": detail})
        print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
        return bool(ok)

    def dump(self, path: Path, extra: dict | None = None) -> dict:
        out = {"steps": self.rows, "passed": sum(r["pass"] for r in self.rows), "total": len(self.rows)}
        if extra:
            out.update(extra)
        path.write_text(json.dumps(out, indent=2, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))
        print(f"{out['passed']}/{out['total']} steps passed -> {path.name}")
        return out


def word_op(ops, word):
    """C_u for u = (u_1, ..., u_n): u_1 acts first."""
    d = ops[0].shape[0]
    m = np.eye(d, dtype=complex)
    for a in word:
        m = ops[a] @ m
    return m


def gns_quotient(kernel, n_letters: int, L_max: int, tol: float = TOL):
    """Grow independent basis words from kernel values only.

    kernel(u, v) -> omega(C_u^dag C_v).  Returns basis words, ranks per length, K_BB and the induced G_a.
    """
    basis: list[tuple] = []
    gram = np.zeros((0, 0), dtype=complex)

    def try_add(w):
        nonlocal gram
        col = np.array([kernel(b, w) for b in basis], dtype=complex)
        kww = kernel(w, w)
        if basis:
            coeff = np.linalg.solve(gram, col)
            resid = kww - col.conj() @ coeff
        else:
            resid = kww
        if resid.real > tol * max(1.0, abs(kww)):
            basis.append(w)
            g = np.zeros((len(basis), len(basis)), dtype=complex)
            g[:-1, :-1] = gram
            g[:-1, -1] = col
            g[-1, :-1] = col.conj()
            g[-1, -1] = kww
            gram = g
            return True
        return False

    try_add(())
    ranks = [len(basis)]
    frontier = list(basis)
    closed_at = None
    for L in range(1, L_max + 1):
        new = []
        for w in frontier:
            for a in range(n_letters):
                if try_add(w + (a,)):
                    new.append(w + (a,))
        ranks.append(len(basis))
        frontier = new
        if not new:
            closed_at = L - 1
            break
    G = []
    if basis:
        for a in range(n_letters):
            kab = np.array([[kernel(u, v + (a,)) for v in basis] for u in basis], dtype=complex)
            G.append(np.linalg.solve(gram, kab))
    return {"basis": basis, "ranks": ranks, "closed_at": closed_at, "gram": gram, "G": G}


def realized_kernel(q, u, v):
    """Kernel predicted by the minimal realization: c_u^dag K_BB c_v with c_w = G_w e_empty."""
    r = len(q["basis"])
    e0 = np.zeros(r, dtype=complex)
    e0[q["basis"].index(())] = 1

    def coords(w):
        c = e0.copy()
        for a in w:
            c = q["G"][a] @ c
        return c

    return coords(u).conj() @ q["gram"] @ coords(v)


def all_words(n_letters: int, L: int):
    for n in range(L + 1):
        yield from product(range(n_letters), repeat=n)


def kernel_from_ops(ops, omega_vec):
    def k(u, v):
        return np.vdot(word_op(ops, u) @ omega_vec, word_op(ops, v) @ omega_vec)
    return k


def purify(rho):
    """Cyclic vector on H (x) H for the state rho: |Omega> = vec(sqrt(rho))."""
    w, v = np.linalg.eigh(rho)
    s = v @ np.diag(np.sqrt(np.clip(w, 0, None))) @ v.conj().T
    return s.reshape(-1)


def amplify(ops):
    """C -> C (x) I on the purification space."""
    d = ops[0].shape[0]
    return [np.kron(c, np.eye(d)) for c in ops]


def span_dim(mats, tol: float = 1e-9) -> int:
    if not mats:
        return 0
    a = np.array([np.asarray(m).ravel() for m in mats])
    s = np.linalg.svd(a, compute_uv=False)
    return int((s > tol * max(1.0, s[0])).sum())


def orth_basis(mats, tol: float = 1e-9):
    a = np.array([np.asarray(m).ravel() for m in mats])
    u, s, vh = np.linalg.svd(a, full_matrices=False)
    r = int((s > tol * max(1.0, s[0])).sum())
    return vh[:r]


def observability_closure(kraus_sets, effects, tol: float = 1e-9):
    """O_0 = span{effects}; O_{n+1} = O_n + sum_a Phi_a^dag(O_n).  Returns basis rows and the dimension history."""
    d = effects[0].shape[0]
    basis = orth_basis(effects, tol)
    dims = [basis.shape[0]]
    while True:
        cand = [b.reshape(d, d) for b in basis]
        for ks in kraus_sets:
            for b in basis:
                x = b.reshape(d, d)
                cand.append(sum(k.conj().T @ x @ k for k in ks))
        nb = orth_basis(cand, tol)
        dims.append(nb.shape[0])
        if nb.shape[0] == basis.shape[0]:
            return nb, dims
        basis = nb


def apply_channel(ks, rho):
    return sum(k @ rho @ k.conj().T for k in ks)


# ----------------------------------------------------------------- Majorana monomials (exact signs)
def mono_mul(p, q):
    """(phase_p, mask_p) * (phase_q, mask_q) with gamma_S ordered ascending."""
    (ph1, s), (ph2, t) = p, q
    sign = 0
    tt = t
    while tt:
        low = tt & -tt
        j = low.bit_length() - 1
        sign += bin(s >> (j + 1)).count("1")      # elements of S greater than j that gamma_j must pass
        tt ^= low
    return (ph1 * ph2 * (-1) ** sign, s ^ t)       # repeated elements cancel: gamma_j^2 = 1


def mono_dag(p):
    ph, s = p
    n = bin(s).count("1")
    return (np.conj(ph) * (-1) ** (n * (n - 1) // 2), s)


def jordan_wigner(n_modes: int):
    """Annihilators a_j on (C^2)^n, a_j = Z...Z sigma^- I...I."""
    z = np.diag([1.0, -1.0])
    sm = np.array([[0.0, 1.0], [0.0, 0.0]])
    ops = []
    for j in range(n_modes):
        m = np.array([[1.0]])
        for k in range(n_modes):
            m = np.kron(m, z if k < j else (sm if k == j else np.eye(2)))
        ops.append(m.astype(complex))
    return ops
