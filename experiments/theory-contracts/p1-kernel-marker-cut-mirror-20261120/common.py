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
