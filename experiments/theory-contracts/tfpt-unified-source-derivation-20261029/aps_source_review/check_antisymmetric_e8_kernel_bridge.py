#!/usr/bin/env python3
"""Exact phase, equivariance and positive-compression audit for the 6 -> 3 bridge.

The six input vectors are the actual ordered mixed-spin products in the
natural q/p section frames.  This matters: in those frames the source spin
half-forms give R mu_ab = rho^-1 mu_(b+1,a+1), rather than a phase-free
permutation.  The orbit-normalized U6 frame is unitarily equivalent, but it
hides the raw R^3=-transpose comparison used by the E8 antisymmetrizer.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as S


HERE = Path(__file__).resolve().parent
OUT = HERE / "ANTISYMMETRIC-E8-KERNEL-BRIDGE.json"

SOURCE_NAMES = [
    "MIXED-SPIN-PRODUCT-U6-SOURCE-REVIEW.txt",
    "OFFDIAGONAL-K-E8-PRODUCT-TEST.txt",
    "ORIGINAL-FUCHS-ALGEBRAIC-FRAME.txt",
    "SPIN-PRODUCT-BOUNDARY-WEIGHT-REVIEW.txt",
    "SPIN-PRODUCT-CAR-U6-TRANSPORT-REVIEW.txt",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix_strings(matrix: S.Matrix) -> list[list[str]]:
    return [[str(S.simplify(matrix[i, j])) for j in range(matrix.cols)]
            for i in range(matrix.rows)]


checks: list[str] = []


def check(condition: bool, label: str) -> None:
    if not bool(condition):
        raise RuntimeError(f"FAILED: {label}")
    checks.append(label)


I = S.I
sqrt3 = S.sqrt(3)
rho = S.Rational(1, 2) + I * sqrt3 / 2
omega = -S.Rational(1, 2) + I * sqrt3 / 2
lam = sqrt3 / 2 + I / 2

check(S.simplify(rho**2 - omega) == 0, "rho^2=omega")
check(S.simplify(lam**2 - rho) == 0, "lambda^2=rho")
check(S.simplify((lam / rho)**2 - 1 / rho) == 0,
      "two R-pulled spin sections give rho^-1")
check(S.simplify((-I)**2 + 1) == 0,
      "two s-pulled spin sections give the scalar -1")

# Natural ordered-pair basis and the oriented wedge basis.
pairs = [(0, 1), (1, 0), (1, 2), (2, 1), (2, 0), (0, 2)]
pair_index = {pair: j for j, pair in enumerate(pairs)}
wedges = [(0, 1), (1, 2), (2, 0)]

# c e_a=e_(a+1), actual v117 quarter K_Q e0=-e1, e1=e0, e2=e2.
c = S.Matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
v4a = S.diag(1, -1, -1)
v4b = S.diag(-1, 1, -1)
kq = S.Matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])


def ordered_pair_action(g: S.Matrix) -> S.Matrix:
    """g tensor g on off-diagonal ordered pairs, using the fixed real B."""
    out = S.zeros(6, 6)
    for col, (a, b) in enumerate(pairs):
        for i in range(3):
            for j in range(3):
                if i != j:
                    out[pair_index[(i, j)], col] += g[i, a] * g[j, b]
    return out


def wedge_action(g: S.Matrix) -> S.Matrix:
    out = S.zeros(3, 3)
    for col, (a, b) in enumerate(wedges):
        for row, (i, j) in enumerate(wedges):
            out[row, col] = S.simplify(g[i, a] * g[j, b] - g[i, b] * g[j, a])
    return out


# A(mu_ab-mu_ba)/sqrt(2)=e_a wedge e_b.
A = S.zeros(3, 6)
for row, (a, b) in enumerate(wedges):
    A[row, pair_index[(a, b)]] = 1 / S.sqrt(2)
    A[row, pair_index[(b, a)]] = -1 / S.sqrt(2)

swap = S.zeros(6, 6)
perm = S.zeros(6, 6)
for col, (a, b) in enumerate(pairs):
    swap[pair_index[(b, a)], col] = 1
    perm[pair_index[((b + 1) % 3, (a + 1) % 3)], col] = 1

# Raw spin-product action.  The scalar follows directly from
# R*Phi=lambda Chi and rho*w+omega^a=rho(w-omega^(a+1)).
Rraw = S.simplify(perm / rho)
sraw = -swap

check(S.simplify(Rraw**3 - sraw) == S.zeros(6), "raw R^3=-transpose")
check(S.simplify(Rraw**6 - S.eye(6)) == S.zeros(6), "raw R^6=I")
check(S.simplify(A * sraw - A) == S.zeros(3, 6),
      "antisymmetric sector has raw s=+I")
check(S.simplify(A * Rraw - omega * wedge_action(c) * A) == S.zeros(3, 6),
      "A Rraw=omega Lambda2(c) A")

# Compare with the old positive U6 orbit frame.  In forward orbit order set
# e_j=rho^(-j) mu_(f^j(01)); then R e_j=e_(j+1).  The j+3 partner contains
# the forced minus rho^(-3)=-1, so raw antisymmetry becomes the +1 sheet of
# the phase-free shift, not its -1 sheet.
orbit = [(0, 1)]
for _ in range(5):
    a, b = orbit[-1]
    orbit.append(((b + 1) % 3, (a + 1) % 3))
orbit_frame = S.zeros(6, 6)
positive_shift = S.zeros(6, 6)
for j, pair in enumerate(orbit):
    orbit_frame[pair_index[pair], j] = rho**(-j)
    positive_shift[(j + 1) % 6, j] = 1
check(S.simplify(Rraw * orbit_frame - orbit_frame * positive_shift) == S.zeros(6),
      "raw R is the phase-free positive U6 shift in the orbit frame")
raw_antisym = (S.eye(6) - swap) / 2
positive_sheet = (S.eye(6) + positive_shift**3) / 2
check(S.simplify(orbit_frame.inv() * raw_antisym * orbit_frame - positive_sheet)
      == S.zeros(6), "raw antisymmetry is U6^3=+1 in the positive orbit frame")

check(S.simplify(A * A.H - S.eye(3)) == S.zeros(3), "A is a coisometry")
check(S.simplify(A.H * A - (S.eye(6) - swap) / 2) == S.zeros(6),
      "A^*A is the antisymmetric projector")

for name, g in (("c", c), ("V4a", v4a), ("V4b", v4b), ("KQ", kq)):
    check(g.T * g == S.eye(3), f"{name} is B-orthogonal")
    check(g.det() == 1, f"{name} preserves the chosen volume")
    check(S.simplify(A * ordered_pair_action(g) - wedge_action(g) * A)
          == S.zeros(3, 6), f"A intertwines {name}")

# With orientation e0 e1 e2=1, Hodge sends 01->e2, 12->e0, 20->e1.
hodge = S.Matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
W = hodge * A
for name, g in (("c", c), ("V4a", v4a), ("V4b", v4b), ("KQ", kq)):
    check(S.simplify(W * ordered_pair_action(g) - g * W) == S.zeros(3, 6),
          f"Hodge-B map intertwines {name}")
check(S.simplify(W * Rraw - omega * c * W) == S.zeros(3, 6),
      "Hodge-B target carries the fixed omega*c R action")

# The exact bare positive structural Gram, with a nonzero regulator included.
y = S.Rational(2, 3)
delta = S.Rational(1, 2)
epsilon = S.Rational(1, 7)
H6 = S.simplify((y * S.eye(6) - delta * Rraw).H
                * (y * S.eye(6) - delta * Rraw)
                + epsilon**2 * S.eye(6))
H3 = S.simplify(W * H6 * W.H)
H3_expected = S.simplify(
    (y**2 + delta**2 + epsilon**2) * S.eye(3)
    - y * delta * (omega * c + S.conjugate(omega) * c.T)
)
check(S.simplify(H3 - H3_expected) == S.zeros(3),
      "positive Gram compresses exactly to the three-sector Gram")

expected_eigenvalues = {
    S.simplify((y - delta)**2 + epsilon**2): 1,
    S.simplify(y**2 + y * delta + delta**2 + epsilon**2): 2,
}
check(H3.eigenvals() == expected_eigenvalues,
      "compressed Gram has the exact 1+2 positive spectrum")
check(all(value > 0 for value in expected_eigenvalues), "compressed Gram is positive definite")

# The inverse is an exact rational example of K3=A K6 A^*.  Functional
# calculus gives the same restriction identity for H^(-1/2).
K6_green = H6.inv()
K3_green = S.simplify(W * K6_green * W.H)
check(S.simplify(K3_green - H3.inv()) == S.zeros(3),
      "Green kernel compression equals the inverse compressed Gram")

# K_Q transports this source to its conjugate; it is not a fixed symmetry of
# the selected c-oriented Gram.
quarter_commutator = S.simplify(kq * H3 - H3 * kq)
check(quarter_commutator != S.zeros(3),
      "the actual quarter gives covariance, not fixed K3 invariance")

result = {
    "status": "PASS",
    "check_count": len(checks),
    "checks": checks,
    "basis": {
        "ordered_pairs": [f"{a}{b}" for a, b in pairs],
        "wedges": [f"{a}{b}" for a, b in wedges],
        "family": ["e0", "e1", "e2"],
    },
    "raw_spin_transport": {
        "R_scalar": "rho^-1",
        "R_label_map": "(a,b)->(b+1,a+1)",
        "s_action": "mu_ab->-mu_ba",
        "antisymmetric_R_action": "omega*Lambda2(c)",
        "antisymmetric_R_cubed": "+I",
        "positive_U6_orbit": [f"{a}{b}" for a, b in orbit],
        "positive_U6_sheet": "U6^3=+I",
    },
    "maps": {
        "A": matrix_strings(A),
        "Hodge_times_A": matrix_strings(W),
        "KQ": matrix_strings(kq),
    },
    "positive_kernel_example": {
        "y": str(y),
        "delta": str(delta),
        "epsilon": str(epsilon),
        "H3": matrix_strings(H3),
        "H3_eigenvalues": {str(k): v for k, v in expected_eigenvalues.items()},
        "quarter_commutator_nonzero": True,
    },
    "source_sha256": {name: sha256(HERE / name) for name in SOURCE_NAMES},
    "scope": [
        "exact coefficient-bundle and positive-kernel compression; the full map is id_KC tensor A",
        "the K_C one-form and full quarter spin/endpoint factors are retained",
        "no identification of v117 KQ with the native full-E8 kappa",
        "no scalar Yukawa or physical P1 response selection",
    ],
}

OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(json.dumps({"status": result["status"], "check_count": result["check_count"],
                  "output": OUT.name}, sort_keys=True))
