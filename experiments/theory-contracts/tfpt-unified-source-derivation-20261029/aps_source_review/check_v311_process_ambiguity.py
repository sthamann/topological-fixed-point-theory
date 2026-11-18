"""Exact source-relative process ambiguity, not a physical TFPT source.

The concrete 3-state matrix is reconstructed exactly from verification/v311's
spectral projectors and checked against that module's numerical implementation.
The CPTP channels below are diagnostic extensions, not claims that v311 or TFPT
physically chooses these channels or identifies their qutrit with E8 or CAR.
"""
from pathlib import Path
import sys

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "verification"))
from v311_gap_decoupled_measure import transfer

R = sp.Rational
I = sp.eye(3)
v, w = sp.Matrix([1, -1, 0]), sp.Matrix([1, 1, -2])
T = sp.ones(3) / 3 + R(64, 729) * v * v.T / 2 + R(1, 729) * w * w.T / 6
eps = sp.symbols("eps", real=True)
limit = min(T[i, i] for i in range(3))


def check(label, condition):
    if not condition:
        raise RuntimeError("FAIL: " + label)
    print("PASS", label)


def channel(x, e):
    """E(X)=e X+sum_ij(T_ij-e delta_ij) E_ij X E_ij*."""
    result = e * x
    for i in range(3):
        for j in range(3):
            unit = sp.zeros(3)
            unit[i, j] = 1
            result += (T[i, j] - e * int(i == j)) * unit * x * unit.T
    return sp.simplify(result)


def choi(e):
    result = sp.zeros(9)
    for a in range(3):
        for b in range(3):
            unit = sp.zeros(3)
            unit[a, b] = 1
            result += sp.kronecker_product(channel(unit, e), unit)
    return sp.simplify(result)


native, lam2, lam3 = transfer()
check("v311 numerical matrix matches exact spectral-projector expression",
      np.allclose(native, np.array(T.tolist(), dtype=float), rtol=0, atol=1e-15))
check("v311 spectral values", abs(lam2 - 64 / 729) < 1e-15 and abs(lam3 - 1 / 729) < 1e-15)
check("T columns sum exactly to one", all(sum(T[i, j] for i in range(3)) == 1 for j in range(3)))
check("all T entries positive", all(T[i, j] > 0 for i in range(3) for j in range(3)))
check("sufficient positivity interval endpoint", limit == R(731, 2187))
check("all diagonal coefficients nonnegative at both interval endpoints",
      all(T[i, i] >= 0 and T[i, i] - limit >= 0 for i in range(3)))

J = choi(eps)
omega = sp.zeros(9, 1)
for i in range(3):
    omega[3 * i + i, 0] = 1
diagonal = sp.zeros(9)
for i in range(3):
    for j in range(3):
        diagonal[3 * i + j, 3 * i + j] = T[i, j] - eps * int(i == j)
check("Choi decomposition exact", sp.simplify(J - eps * omega * omega.T - diagonal) == sp.zeros(9))
check("Choi partial output trace is identity for all eps",
      sp.Matrix(3, 3, lambda a, b: sp.simplify(sum(J[3*i+a, 3*i+b] for i in range(3)))) == I)
print("PROOF CPTP: Choi is eps*|Omega><Omega| plus diagonal nonnegative weights")
print("            for every 0 <= eps <=", limit, "; trace is identity.")

for a in range(3):
    for b in range(3):
        x = sp.zeros(3)
        x[a, b] = 1
        y = channel(x, eps)
        expected = T[:, a] if a == b else sp.zeros(3, 1)
        check(f"population intertwining E({a},{b})",
              sp.Matrix([sp.simplify(y[i, i]) for i in range(3)]) == expected)

print("PROOF HISTORY: diag(E(X))=T diag(X) for every X; induction gives")
print("               diag(E^n(X))=T^n diag(X) for all n >= 0.")
rho = sp.Matrix([[R(1, 2), R(1, 2), 0], [R(1, 2), R(1, 2), 0], [0, 0, 0]])
check("diagnostic density matrix normalized and rank one", rho.trace() == 1 and rho * rho == rho)
check("exact coherence response for every eps", channel(rho, eps)[0, 1] == eps / 2)
for e in (R(0), R(1, 5), R(1, 3)):
    check(f"epsilon {e} lies in proved CPTP interval", 0 <= e <= limit)
    print(f"  eps={e}, output coherence={channel(rho, e)[0, 1]}")
check("distinct exact outputs at same population history", channel(rho, 0) != channel(rho, R(1, 3)))
print("SCOPE: Only the v311 population reader is fixed. Additional original TFPT")
print("       source constraints may remove any or all diagnostic extensions.")
