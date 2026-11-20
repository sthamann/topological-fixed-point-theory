"""Shared helpers for P1BIND20261120. Linear algebra in double precision with fixed seeds; no fitted inputs."""
import json, sys
import numpy as np

TOL = 1e-9

class Checks:
    def __init__(self, name): self.name, self.rows = name, []
    def ok(self, label, cond, detail=""):
        self.rows.append({"check": label, "pass": bool(cond), "detail": str(detail)})
        print(("PASS " if cond else "FAIL ") + label + ("  | " + str(detail) if detail != "" else ""))
    def done(self, extra=None):
        n = sum(r["pass"] for r in self.rows)
        print(f"{self.name}: {n}/{len(self.rows)} checks passed")
        json.dump({"name": self.name, "checks": self.rows, "data": extra or {}},
                  open(f"results_{self.name}.json", "w"), indent=1, ensure_ascii=False, default=str)
        sys.exit(0 if n == len(self.rows) else 1)

def odd(D):
    """Self-adjoint operator anticommuting with the normal grading sigma_z x 1, from its off-diagonal block D."""
    m = D.shape[0]
    Z = np.zeros((m, m), complex)
    return np.block([[Z, D.conj().T], [D, Z]])

def grading(m):
    return np.diag([1.0] * m + [-1.0] * m).astype(complex)

def spectral(B, sign):
    w, v = np.linalg.eigh(B)
    cols = v[:, w * sign > 0]
    return cols @ cols.conj().T

def expmh(H, t):
    """exp(t H) for self-adjoint H."""
    w, v = np.linalg.eigh(H)
    return (v * np.exp(t * w)) @ v.conj().T

def projector(cols):
    q, _ = np.linalg.qr(cols)
    return q @ q.conj().T

def randc(rng, *shape):
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


def fock(n):
    """Annihilation operators a_0..a_{n-1} on the 2^n dimensional Fock space (Jordan-Wigner)."""
    sz, lo, one = np.diag([1.0, -1.0]), np.array([[0.0, 1.0], [0.0, 0.0]]), np.eye(2)
    ops = []
    for i in range(n):
        m = np.array([[1.0]])
        for j in range(n):
            m = np.kron(m, sz if j < i else (lo if j == i else one))
        ops.append(m.astype(complex))
    return ops


def dgamma(A, ann):
    """Second quantisation sum_ij A_ij a_i^+ a_j."""
    n = len(ann)
    out = np.zeros_like(ann[0])
    for i in range(n):
        for j in range(n):
            if A[i, j] != 0:
                out = out + A[i, j] * ann[i].conj().T @ ann[j]
    return out


def slater(cols, ann):
    """Fock vector with the orthonormal columns of cols filled."""
    q, _ = np.linalg.qr(cols)
    v = np.zeros(ann[0].shape[0], complex)
    v[0] = 1.0                                   # all sites empty (first basis vector is |0...0>)
    for k in range(q.shape[1]):
        c = sum(q[i, k] * ann[i].conj().T for i in range(len(ann)))
        v = c @ v
    return v / np.linalg.norm(v)


def twopoint(v, ann):
    """C_ij = <a_j^+ a_i> ; for a Slater state this is the projector onto the filled space."""
    n = len(ann)
    return np.array([[v.conj() @ (ann[j].conj().T @ (ann[i] @ v)) for j in range(n)] for i in range(n)])
