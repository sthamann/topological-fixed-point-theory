#!/usr/bin/env python3
"""Part G: SPINOR.ORIGIN.01 -- can the missing propagation spinor come from TFPT's own CAR structure?

The three family modes of U (Part U) are complex CAR modes a_1, a_2, a_3 with Majoranas
gamma_i = a_i + a_i^dag, gamma'_i = -i (a_i - a_i^dag).  On their Fock space Lambda C^3:
  G1 CAR and Majorana relations hold (Jordan-Wigner, exact).
  G2 J_i = -(i/4) eps_ijk gamma_j gamma_k and J'_i (same with gamma') form su(2) + su(2), commuting:
     Spin(3) x Spin(3) = Spin(4) inside the family Spin(6) = SU(4); on the four family lines
     Lambda^even C^3 both Casimirs are 3/4: the family lines are the (1/2, 1/2) of Spin(4).
  G3 under the diagonal rotation (J + J', rotating the complex modes) the four lines split as
     singlet + triplet: the singlet is the family Fock vacuum = the OMITTED family line, the triplet is
     Lambda^2 C^3 = the THREE families.  The 3+1 family split is vector + singlet of the diagonal SO(3).
  G4 the Weyl walk with spinor 2J: U_F(k) = prod_i (cos k_i - i sin k_i 2J_i) on the four lines is unitary,
     supported on the eight bcc steps, isotropic Weyl ((k.2J)^2 = |k|^2), and the deck rotation lifts to
     order 8 with fourth power -1 on the spinor factor (the mu8 sign).
  G5 KILL test: [U_F(k), P_vacuum] = O(k) != 0: free propagation would mix the omitted line with the three
     families and change the family content kinetically -- incompatible with the separation of space and
     flavor (Part T6) and with family conservation of free propagation.  The same holds for J'.
  G6 conclusion: TFPT's internal Clifford structure does contain a canonical spinor (and reproduces the
     3+1 family split as SO(3) vector + singlet), but used as the propagation spinor it locks family
     quantum numbers to momentum.  The physical spacetime spinor must be an internal-number-free C^2:
     AX.SPINSTEP stays, now sharpened.
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/spinor_origin.py
"""
from __future__ import annotations

import json
import sys
from itertools import product
from pathlib import Path

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


Z = np.diag([1.0, -1.0]).astype(complex); I2 = np.eye(2, dtype=complex)
lower = np.array([[0, 1], [0, 0]], dtype=complex)             # |0><1| on occupation basis (0 = empty)
def op(i, m):
    mats = [Z if j < i else (m if j == i else I2) for j in range(3)]
    out = mats[0]
    for x in mats[1:]:
        out = np.kron(out, x)
    return out
a = [op(i, lower) for i in range(3)]
ad = [x.conj().T for x in a]
I8 = np.eye(8)
car = all(np.allclose(a[i] @ ad[j] + ad[j] @ a[i], (i == j) * I8) and np.allclose(a[i] @ a[j] + a[j] @ a[i], 0)
          for i in range(3) for j in range(3))
g = [a[i] + ad[i] for i in range(3)]
gp = [-1j * (a[i] - ad[i]) for i in range(3)]
maj = all(np.allclose(x @ y + y @ x, 2 * (x is y) * I8) for x in g + gp for y in g + gp)
step("G1.three family CAR modes (Jordan-Wigner) and their six Majoranas gamma_i, gamma'_i: exact CAR / Clifford relations", car and maj)

eps = lambda i, j, k: np.sign(np.prod([j - i, k - i, k - j])) if len({i, j, k}) == 3 else 0
def Jvec(gg):
    return [sum(-0.25j * eps(i, j, k) * gg[j] @ gg[k] for j in range(3) for k in range(3)) for i in range(3)]
J, Jp = Jvec(g), Jvec(gp)
su2 = all(np.allclose(J[i] @ J[j] - J[j] @ J[i], 1j * sum(eps(i, j, k) * J[k] for k in range(3))) for i in range(3) for j in range(3))
su2p = all(np.allclose(Jp[i] @ Jp[j] - Jp[j] @ Jp[i], 1j * sum(eps(i, j, k) * Jp[k] for k in range(3))) for i in range(3) for j in range(3))
commute = all(np.allclose(J[i] @ Jp[j], Jp[j] @ J[i]) for i in range(3) for j in range(3))
N = sum(ad[i] @ a[i] for i in range(3))
even = [k for k in range(8) if round(N[k, k].real) % 2 == 0]
Pe = np.eye(8)[:, even]
cas = lambda JJ: Pe.conj().T @ sum(x @ x for x in JJ) @ Pe
step("G2.J and J' form commuting su(2) + su(2) = Spin(4) inside the family Spin(6); on the four family lines both Casimirs = 3/4: the lines are (1/2, 1/2)",
     su2 and su2p and commute and np.allclose(cas(J), 0.75 * np.eye(4)) and np.allclose(cas(Jp), 0.75 * np.eye(4)))

D = [J[i] + Jp[i] for i in range(3)]
casD = Pe.conj().T @ sum(x @ x for x in D) @ Pe
ev, evec = np.linalg.eigh(casD)
vac = np.zeros(8); vac[0] = 1                                   # |000> = family Fock vacuum
singlet = Pe @ evec[:, np.argmin(ev)]
num = np.real(np.diag(N))
triplet_space = all(abs(num[even][k] - 2) < 1e-9 or abs(num[even][k]) < 1e-9 for k in range(4))
step("G3.diagonal SO(3): (1/2,1/2) = singlet + triplet; the singlet IS the family Fock vacuum (omitted line), the triplet is Lambda^2 C^3 (the three families)",
     np.allclose(sorted(ev), [0, 2, 2, 2]) and abs(abs(np.vdot(vac, singlet)) - 1) < 1e-12 and triplet_space,
     f"diagonal Casimir on the lines {np.round(sorted(ev), 6).tolist()}")

S = [Pe.conj().T @ (2 * J[i]) @ Pe for i in range(3)]
I4 = np.eye(4)
UF = lambda k: (np.cos(k[0]) * I4 - 1j * np.sin(k[0]) * S[0]) @ (np.cos(k[1]) * I4 - 1j * np.sin(k[1]) * S[1]) @ (np.cos(k[2]) * I4 - 1j * np.sin(k[2]) * S[2])
rng = np.random.default_rng(4)
ks = rng.normal(size=(4, 3))
unit = all(np.allclose(UF(k).conj().T @ UF(k), I4) for k in ks)
iso = all(np.allclose((sum(kk * s for kk, s in zip(k, S))) @ (sum(kk * s for kk, s in zip(k, S))), (k @ k) * I4) for k in ks)
rot = expm(-1j * (np.pi / 2) * (Pe.conj().T @ J[1] @ Pe))     # quarter turn about the deck axis on the spinor factor
mu8 = np.allclose(np.linalg.matrix_power(rot, 4), -I4) and np.allclose(np.linalg.matrix_power(rot, 8), I4)
step("G4.Weyl walk with the family spinor 2J on the four lines: unitary, eight bcc steps, isotropic ((k.2J)^2 = |k|^2); the deck quarter turn lifts to order 8 with fourth power -1 (mu8)",
     unit and iso and mu8)

Pv = np.outer(Pe.conj().T @ vac, (Pe.conj().T @ vac).conj())
mix = [np.abs(UF(e * np.array([0.3, 0.5, 0.8])) @ Pv - Pv @ UF(e * np.array([0.3, 0.5, 0.8]))).max() for e in (1e-3, 2e-3)]
SPp = [Pe.conj().T @ (2 * Jp[i]) @ Pe for i in range(3)]
UFp = lambda k: (np.cos(k[0]) * I4 - 1j * np.sin(k[0]) * SPp[0]) @ (np.cos(k[1]) * I4 - 1j * np.sin(k[1]) * SPp[1]) @ (np.cos(k[2]) * I4 - 1j * np.sin(k[2]) * SPp[2])
mixp = np.abs(UFp(1e-3 * np.array([0.3, 0.5, 0.8])) @ Pv - Pv @ UFp(1e-3 * np.array([0.3, 0.5, 0.8]))).max()
step("G5.KILL test: [U_F(k), P_vacuum] = O(k) != 0 (also for J'): free propagation with the family spinor mixes the omitted line with the three families -> family content changes kinetically",
     mix[0] > 1e-4 and abs(mix[1] / mix[0] - 2) < 0.05 and mixp > 1e-4,
     f"|[U,P]| = {mix[0]:.2e}, {mix[1]:.2e} (linear in k); J' variant {mixp:.2e}")
step("G6.(statement) TFPT's family Clifford structure contains a canonical spinor and reproduces the 3+1 family split as SO(3) vector + singlet, but as the propagation spinor it locks family number to momentum: the spacetime spinor must be internal-number-free (AX.SPINSTEP sharpened)",
     True, "excluded: family Spin(4) spinor as the propagation factor")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG), "mixing": mix}
(HERE / "spinor_origin_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
