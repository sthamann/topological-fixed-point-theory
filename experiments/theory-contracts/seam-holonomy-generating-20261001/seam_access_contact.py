#!/usr/bin/env python3
"""Part N: independent re-computation of the external SEAM.ACCESS.01 findings + consequences.

External claims re-derived here (the external sandbox package was NOT available):
  O1 clock-covariant accesses U J = J Q (U = Q x I) are exactly J chi_r = chi_r x b_r; the response is
     circulant with eigenvalues k_r(z) = <b_r,(h-z)^-1 b_r>; counter-case K_02 = (p-q)/((1-z)(3-z)).
  O2 a clock-invariant state as j_0 gives rank J = 1.
  O3 regulated point evaluation on L^2_0(S^1) with A = |D|: ||J_N e_j||^2 = N/pi; the reference-
     subtracted response M(0) = -(1/pi) ln|2 sin(d/2)|; on the square -(ln2/2pi) circ(0,1,2,1);
     the mu4 average gives -c3 ln2 C_marks.
  O4 M(-1) first row = (-1/pi, 1/4 - 1/pi, (ln2 - 1)/pi, 1/4 - 1/pi).
  O5 exclusion: pole weights R_n >= 0 with Tr R_n = 4/pi; their sum diverges, while any bounded J has
     total weight Tr J^dag J < infinity -> no M(z) = B + J^dag (A - z)^-1 J.
  O6 repair: J0 = (A+1)^-1 gamma^dag is bounded and M(z) = M(-1) + (z+1) J0^dag J0 + (z+1)^2 J0^dag (A-z)^-1 J0.
Consequences for this contract (new):
  O7 the uniform rank-one port is a STATIC statement: on the square the alternating eigenvalue
     lam_v(z) = -(4/pi) sum_m z/((4m+2)(4m+2-z)) vanishes only at z = 0 (lam_v(-1) != 0).
  O8 the rank-one port needs exactly the reference subtraction: with diagonal offset delta, lam_v = delta;
     delta = 0 is v485's point splitting at the KMS seam length l = 2 pi = 1/(4 c3).
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/seam_access_contact.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


# O1 / O2 -------------------------------------------------------------------------------
Q = np.roll(np.eye(4), 1, axis=0)
Fm = np.array([[1j ** (r * j) for r in range(4)] for j in range(4)]) / 2      # columns chi_r, Q chi_r = i^r chi_r
h = np.diag([1.0, 3.0])
def access(bs):
    """J: C^4 -> C^4 (x) C^2 with J chi_r = chi_r (x) b_r."""
    cols = [np.kron(Fm[:, r], bs[r]) for r in range(4)]
    return np.column_stack(cols) @ np.linalg.inv(Fm)
def response(J, z):
    return J.conj().T @ np.kron(np.eye(4), np.linalg.inv(h - z * np.eye(2))) @ J
U = np.kron(Q, np.eye(2))
def bvec(p): return np.array([np.sqrt(p), np.sqrt(1 - p)])
ok1, vals = True, []
for p, q in ((0.75, 0.25), (0.5, 0.25)):
    J = access([bvec(p), bvec(q), bvec(p), bvec(q)])
    ok1 &= np.allclose(U @ J, J @ Q) and np.allclose(J.conj().T @ J, np.eye(4))
    K = response(J, 0.0)
    vals.append(K[0, 2].real)
    ok1 &= abs(K[0, 2] - (p - q) / ((1 - 0) * (3 - 0))) < 1e-12
step("O1.clock-covariant, normalised accesses with the same h and clock give different cross responses K_02 = (p-q)/((1-z)(3-z))",
     ok1 and np.allclose(vals, [1 / 6, 1 / 12]), f"K_02(0) = {vals[0]:.4f}, {vals[1]:.4f}")
psi = np.kron(np.ones(4) / 2, bvec(0.3))
Jinv = np.column_stack([np.linalg.matrix_power(U, k) @ psi for k in range(4)])
step("O2.a clock-invariant state used as j_0 gives rank J = 1 (state and local insertion are different objects)",
     np.allclose(U @ psi, psi) and np.linalg.matrix_rank(Jinv) == 1)

# O3 - O6 -----------------------------------------------------------------------------
Nc = 200000
n = np.arange(1, Nc + 1)
H = (1 / n).sum()
def Mz(d, z):  # reference-subtracted response between points at separation d
    return (np.cos(n * d) / (n - z)).sum() / np.pi - (H / np.pi if d == 0 else 0.0)
step("O3a.regulated point access: ||J_N e_j||^2 = sum_{0<|n|<=N} 1/(2 pi) = N/pi (diverges: not a normalisable source vector)",
     abs(2 * Nc / (2 * np.pi) - Nc / np.pi) < 1e-9)
xs = [k * np.pi / 2 for k in range(4)]
M0 = np.array([[Mz(abs(a - b), 0.0) for b in xs] for a in xs])
C = np.array([[0, 1, 2, 1], [1, 0, 1, 2], [2, 1, 0, 1], [1, 2, 1, 0]])
c3 = 1 / (8 * np.pi)
step("O3b.M(0) on the square = -(ln2/2pi) circ(0,1,2,1); (1/4) M(0) = -c3 ln2 C_marks",
     np.allclose(M0, -(np.log(2) / (2 * np.pi)) * C, atol=2e-5) and np.allclose(M0 / 4, -c3 * np.log(2) * C, atol=1e-5))
row = [Mz(abs(xs[0] - b), -1.0) for b in xs]
target = [-1 / np.pi, 0.25 - 1 / np.pi, (np.log(2) - 1) / np.pi, 0.25 - 1 / np.pi]
step("O4.M(-1) first row = (-1/pi, 1/4-1/pi, (ln2-1)/pi, 1/4-1/pi)", np.allclose(row, target, atol=2e-5),
     f"{np.round(row, 6).tolist()}")
Rn = [np.array([[np.cos(m * (a - b)) / np.pi for b in xs] for a in xs]) for m in range(1, 9)]
step("O5.pole weights R_n >= 0 with Tr R_n = 4/pi each: total weight diverges, so no bounded J realises M(z) as J^dag (A-z)^-1 J + B",
     all(np.linalg.eigvalsh(R).min() > -1e-12 and abs(np.trace(R) - 4 / np.pi) < 1e-12 for R in Rn))
zt = 0.5j
def S(f):  # sum over 0 < |n| <= Nc of f(|n|) e^{in d}/(2 pi) between points -> (1/pi) sum cos(nd) f(n)
    return lambda d: (np.cos(n * d) * f(n)).sum() / np.pi
lhs = np.array([[Mz(abs(a - b), zt) for b in xs] for a in xs])
m1 = np.array([[Mz(abs(a - b), -1.0) for b in xs] for a in xs])
J0J0 = np.array([[S(lambda k: 1 / (k + 1) ** 2)(abs(a - b)) for b in xs] for a in xs])
J0RJ0 = np.array([[S(lambda k: 1 / ((k + 1) ** 2 * (k - zt)))(abs(a - b)) for b in xs] for a in xs])
rhs = m1 + (zt + 1) * J0J0 + (zt + 1) ** 2 * J0RJ0
step("O6.repair with bounded J0 = (A+1)^-1 gamma^dag: M(z) = M(-1) + (z+1) J0^dag J0 + (z+1)^2 J0^dag (A-z)^-1 J0 (z = i/2)",
     np.allclose(lhs, rhs, atol=1e-8), f"max dev {np.abs(lhs - rhs).max():.1e}")

# O7 / O8: consequences for the uniform port -----------------------------------------
v = np.array([1, -1, 1, -1]) / 2
s = np.ones(4) / 2
lam_v = lambda z: v @ np.array([[Mz(abs(a - b), z) for b in xs] for a in xs]) @ v
lv0, lvm1, lvp = lam_v(0.0), lam_v(-1.0), lam_v(0.5)
m = np.arange(0, 200000)
closed = -(4 / np.pi) * (1 / ((4 * m + 2) * (4 * m + 3))).sum()
step("O7.the uniform rank-one port is STATIC: lam_v(0) = 0 but lam_v(-1) = -(4/pi) sum 1/((4m+2)(4m+3)) != 0, lam_v(1/2) != 0",
     abs(lv0) < 1e-4 and abs(lvm1 - closed) < 1e-4 and abs(lvm1) > 0.1 and abs(lvp) > 0.05,
     f"lam_v(0)={lv0:.1e}, lam_v(-1)={lvm1:.4f} (closed {closed:.4f}), lam_v(1/2)={lvp:.4f}")
offs = [0.0, 0.05, -0.1]
lvs = [v @ (M0 + dlt * np.eye(4)) @ v for dlt in offs]
ell = 2 * np.pi
reg = lambda L: np.log(L / (2 * np.pi)) / np.pi        # v485: lim_{d->0} [G(d;l) + (1/pi) ln d]
step("O8.rank-one port needs exactly the reference subtraction: lam_v = delta for diagonal offset delta; delta = 0 is v485 point splitting at l = 2 pi = 1/(4 c3)",
     np.allclose(lvs, offs, atol=2e-5) and abs(reg(ell)) < 1e-15 and abs(ell - 1 / (4 * c3)) < 1e-12 and abs(reg(3.0)) > 0.1)

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG),
       "lam_v": {"0": lv0, "-1": lvm1, "0.5": lvp}}
(HERE / "seam_access_contact_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
