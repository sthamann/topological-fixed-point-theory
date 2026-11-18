#!/usr/bin/env python3
"""Part E: the two remaining selection problems as moduli problems.

  (1) Which condition fixes the cross-ratio lambda = -1 (square, tau = i)?
  (2) Why does the natural three-dimensional H^1 carry an IRREDUCIBLE family action
      with an order-3 element, although the conformal D4 acts reducibly (1+2)?

Tested candidates: extremal functionals (det Laplacian), clock consistency with the
raw station clock R of the antiperiodic seam ring, Moebius covariance.  Run:
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/selection.py
"""
from __future__ import annotations

import json
import sys
from itertools import permutations
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


I = sp.I
cr = lambda a, b, c, d: sp.simplify((a - c) * (b - d) / ((a - d) * (b - c)))

# ------------------------------------------------ E1: the modulus is a genuine parameter
th = sp.symbols("theta", real=True)
B = [sp.Integer(1), sp.exp(I * th), sp.Integer(-1), -sp.exp(I * th)]
lam = sp.simplify(sp.expand_complex(sp.simplify(B[0] - B[1]) * (B[2] - B[3]) / ((B[0] - B[3]) * (B[2] - B[1]))))
lam_target = -sp.tan(th / 2) ** 2
step("E1.B_theta = {1, e^{i th}, -1, -e^{i th}}: cross-ratio = -tan^2(th/2) (one real modulus); square <=> th = pi/2 <=> lambda = -1",
     sp.simplify(sp.expand_complex(lam - lam_target)) == 0 and sp.simplify(lam_target.subs(th, sp.pi / 2)) == -1)

# ------------------------- E2: anharmonic stabilisers: order 2 OR order 3, never both
z = sp.symbols("l")
anh = {"1-l": 1 - z, "1/l": 1 / z, "l/(l-1)": z / (z - 1),            # order 2
       "1/(1-l)": 1 / (1 - z), "(l-1)/l": (z - 1) / z}                 # order 3
fix = {k: set(sp.solve(sp.Eq(v, z), z)) for k, v in anh.items()}
ord2 = fix["1-l"] | fix["1/l"] | fix["l/(l-1)"]
ord3 = fix["1/(1-l)"] | fix["(l-1)/l"]
step("E2.conformal stabiliser of 4 marks = V4 x (stabiliser of lambda in S3): order-2 points {-1,1/2,2} (D4), order-3 points e^{+-i pi/3} (A4); disjoint",
     {sp.nsimplify(x) for x in ord2} - {0, 1} == {-1, sp.Rational(1, 2), 2}
     and all(sp.simplify(x ** 2 - x + 1) == 0 for x in ord3) and not (ord2 & ord3))
step("E3.NO-GO: no four-mark seam carries BOTH a conformal order-4 clock and a conformal order-3 family clock (S4 never conformal)",
     not (ord2 & ord3), "one of them must be a non-conformal (mapping-class / Fuchs) monodromy")

# ----------------------------- E4: H^1 representations: D4 reducible, S4 irreducible
def perm_matrix(p):
    M = sp.zeros(4)
    for i, j in enumerate(p):
        M[j, i] = 1
    return M
def rep_on_H1(p):  # H^1 = residue vectors with sum 0 = permutation rep minus trivial
    return perm_matrix(p)
def char_inner(group):
    chis = [rep_on_H1(p).trace() - 1 for p in group]  # standard = perm - trivial
    return sp.Rational(sum(c * c for c in chis), len(group))
S4 = list(permutations(range(4)))
rho_p, tau_p = (1, 2, 3, 0), (0, 3, 2, 1)
D4 = {tuple(range(4))}
while True:
    new = {tuple(g[h[i]] for i in range(4)) for g in D4 | {rho_p, tau_p} for h in D4 | {rho_p, tau_p}} | D4
    if new == D4:
        break
    D4 = new
step("E4.on H^1 (dim 3): D4 acts reducibly (<chi,chi>=2: 1+2), the full mapping-class image S4 = W(A3) irreducibly (<chi,chi>=1)",
     len(D4) == 8 and char_inner(sorted(D4)) == 2 and char_inner(S4) == 1
     and any(p != tuple(range(4)) and all(p[p[p[i]]] == i for i in range(4)) for p in S4)
     and not any(p != tuple(range(4)) and all(p[p[p[i]]] == i for i in range(4)) for p in D4),
     "Mod(S^2,4) -> S4 is onto for EVERY lambda; the order-3 element is a braid monodromy, conformal only at j=0")

# ------------------------------- E5: extremal principle on the sheet torus (det Laplacian)
mp.mp.dps = 30
def eta(t):
    q = mp.e ** (2j * mp.pi * t)
    return q ** (mp.mpf(1) / 24) * mp.nprod(lambda n: 1 - q ** n, [1, mp.inf])
def logdet(t):  # flat torus of fixed area: log det' Delta = log(Im t |eta(t)|^4) + const
    return mp.log(mp.im(t) * abs(eta(t)) ** 4)
F = lambda x, y: logdet(mp.mpc(x, y))
sq, hx = (0, 1), (0.5, mp.sqrt(3) / 2)
h = mp.mpf("1e-4")
def grad_hess(x, y):
    fx = (F(x + h, y) - F(x - h, y)) / (2 * h)
    fy = (F(x, y + h) - F(x, y - h)) / (2 * h)
    fxx = (F(x + h, y) - 2 * F(x, y) + F(x - h, y)) / h ** 2
    fyy = (F(x, y + h) - 2 * F(x, y) + F(x, y - h)) / h ** 2
    fxy = (F(x + h, y + h) - F(x + h, y - h) - F(x - h, y + h) + F(x - h, y - h)) / (4 * h ** 2)
    return (fx, fy), (fxx * fyy - fxy ** 2)
g_sq, det_sq = grad_hess(*sq)
g_hx, det_hx = grad_hess(*hx)
step("E5.extremal principle FAILS to select the square: log det' Delta is critical at tau=i but a SADDLE; its maximum is hexagonal (j=0)",
     abs(g_sq[0]) < 1e-8 and abs(g_sq[1]) < 1e-8 and det_sq < 0 and det_hx > 0 and F(*hx) > F(*sq),
     f"logdet(i)={mp.nstr(F(*sq), 10)}, logdet(hex)={mp.nstr(F(*hx), 10)}; Osgood-Phillips-Sarnak")

# ---------------------------- E6: the raw antiperiodic station clock is the det-1 lift
R = np.zeros((4, 4), dtype=np.int64)
for x in range(3):
    R[x + 1, x] = 1
R[0, 3] = -1                       # universalraum-mu4-pointing-20260917, station_matrices()
J = np.zeros((4, 4), dtype=np.int64)
J[0, 0] = 1
for x in (1, 2, 3):
    J[4 - x, x] = -1
R4 = np.linalg.matrix_power(R, 4)
proj = [int(np.nonzero(R[:, i])[0][0]) for i in range(4)]
step("E6.raw station clock R: projectively the 4-cycle of the stations, det R = 1, R^4 = -I, R^8 = I; J R J = R^-1",
     proj == [1, 2, 3, 0] and round(np.linalg.det(R)) == 1 and (R4 == -np.eye(4)).all()
     and (np.linalg.matrix_power(R, 8) == np.eye(4)).all() and (J @ R @ J == np.linalg.inv(R).round()).all(),
     "the SU(4) lift with g^4 = -I found in E8 (C10) is already the raw seam's antiperiodic station clock")

# ------------- E7: Moebius covariance turns the raw station clock into the square
phi = sp.symbols("phi", real=True)
a = sp.exp(I * phi)
orbit = [a * I ** k for k in range(4)]
step("E7.an order-4 Moebius map preserving the seam circle is a rotation; every orbit on the circle is harmonic (lambda = -1 / 2), for all phi",
     sp.simplify(cr(*orbit) - 2) == 0,
     "elliptic order-4 element of PSU(1,1) is conjugate to z->iz; orbit {a,ia,-a,-ia}")

# --------- E8: inside E8 the det-1 condition is automatic (no U(1) commutes with Spin(10) x SU(4))
from itertools import combinations, product
roots = set()
for i, j in combinations(range(8), 2):
    for si, sj in product((2, -2), repeat=2):
        r = [0] * 8; r[i], r[j] = si, sj; roots.add(tuple(r))
roots |= {sgn for sgn in product((1, -1), repeat=8) if sgn.count(-1) % 2 == 0}
d5 = [r for r in roots if not any(r[5:]) ]
cart = sp.Matrix([list(r) for r in d5]).nullspace()
cent_roots = [r for r in roots if all(sum(a * b for a, b in zip(r, q)) == 0 for q in d5)]
step("E8.centraliser of Spin(10) in E8 has Lie algebra su(4) exactly (Cartan 3 + 12 roots = 15, no extra u(1)): det-1 lift forced",
     len(cart) == 3 and len(cent_roots) == 12 and len(roots) == 240 and len(d5) == 40)

summary = {"steps": LOG, "passed": sum(s["pass"] for s in LOG), "total": len(LOG),
           "logdet": {"square": mp.nstr(F(*sq), 15), "hexagonal": mp.nstr(F(*hx), 15)}}
(HERE / "selection_results.json").write_text(json.dumps(summary, indent=2))
print(f"{summary['passed']}/{summary['total']} steps passed")
sys.exit(0 if summary["passed"] == summary["total"] else 1)
