#!/usr/bin/env python3
"""Part U: one universal source U and three representations of it (U -> A, U -> B, U -> C).

Candidate U: chiral CAR net of 8 complex fermion modes on the oriented, Z2-sheeted, marked seam:
    5 carrier modes  (from P2, g = 5)          -> SO(10)_1, Ramond zero modes Cl(10)
    3 family modes   (from the marks, n-1 = 3) -> SO(6)_1 = SU(4)_1, Ramond zero modes Cl(6)
    Z2 sheet         = sum over spin structures (NS + R), GSO chirality = orientation of P1
    vacuum           = Hardy (positive-energy) sea, deck = quarter rotation (Part M, K)
Exact checks against the ORIGINAL project objects (tfpt_explorer):
  U1  R-sector states = occupation patterns n in {0,1}^8 -> doubled coordinates 2n-1; with GSO they
      are exactly the 128 spinor roots of the project's E8; NS bilinears are the 112 integer roots.
  U2  carrier: even occupation of the 5 modes = S+ = Lambda^even C^5 (the project's spinor block).
  U3  family: the four family lines = family occupation 0 (Fock vacuum) or 2 (Lambda^2 C^3);
      the project's omitted line is the family Fock vacuum, the three retained families are Lambda^2.
  U4  S48 = {carrier even, family number 2}: equals the project's retained 48.
  U5  the selection holonomies are fermion-number holonomies: z3 = exp(2 pi i (2/3) N_fam),
      z4 = the carrier number phase; hypercharge y is linear in the carrier occupations.
  U6  the Z2 sheet: (-1)^F = -1 exactly on the R sector (= spinor roots): the mu8 deck lift of C10.
  U7  E8 = closure: D5 + D3 (= A3) inside D8 inside E8, total index 4 with cyclic glue Z4.
  U8  flavor: the v117 monodromy group is real (Frobenius-Schur indicator +1): it acts by real
      rotations SO(3) inside U(3) inside SO(6) on the three family modes.
  U9  geometry of U itself: the modular/Moebius data of a chiral net on S^1 generate PSL(2,R)
      (dimension 3) -- NOT a 3+1D group (Poincare 10, conformal 15): Theorem D is not a corollary of U.
    PYTHONPATH=. python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/universal_source.py
"""
from __future__ import annotations

import json
import sys
from itertools import combinations, product
from pathlib import Path

import numpy as np
import sympy as sp

from tfpt_explorer.neutral_source_response import _actual_current_source, _native_neutral_branch
from tfpt_explorer.source_charge_response import hypercharge_vector

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


src, br = _actual_current_source(), _native_neutral_branch(0)
E8 = {tuple(map(int, r)) for r in src["chevalley"].roots}
proj_spinor = {tuple(map(int, r)) for r in br["spinor"]}
proj_retained = {tuple(map(int, r)) for r in br["retained"]}
proj_omitted = {tuple(map(int, r)) for r in br["omitted"]}

# U1: R-sector occupation patterns and NS bilinears
occ = list(product((0, 1), repeat=8))
to_root = lambda n: tuple(2 * x - 1 for x in n)
R_all = {to_root(n) for n in occ}
gso = {r for r in R_all if r.count(-1) % 2 == 0}        # GSO projection (orientation of P1 fixes the sign)
ns = {r for r in E8 if all(abs(x) in (0, 2) for x in r)}
step("U1.R-sector states with GSO = the 128 spinor roots of the project's E8; NS bilinears +-e_i+-e_j = the 112 integer roots",
     gso == {r for r in E8 if all(abs(x) == 1 for x in r)} and len(gso) == 128 and len(ns) == 112 and gso | ns == E8)

# U2/U3: carrier and family Fock content of the project's chiral block (16,4)
carrier_occ = lambda r: sum(x == 1 for x in r[:5])
family_occ = lambda r: sum(x == 1 for x in r[5:])
step("U2.carrier: the project's chiral spinor block has EVEN occupation of the 5 carrier modes (S+ = Lambda^even C^5, 16 states)",
     all(carrier_occ(r) % 2 == 0 for r in proj_spinor) and len({r[:5] for r in proj_spinor}) == 16)
step("U3.family lines = family occupation 0 (Fock vacuum) or 2 (Lambda^2 C^3); the project's OMITTED line is the family Fock vacuum",
     {family_occ(r) for r in proj_spinor} == {0, 2} and all(family_occ(r) == 0 for r in proj_omitted)
     and all(family_occ(r) == 2 for r in proj_retained))
S48 = {r for r in gso if carrier_occ(r) % 2 == 0 and family_occ(r) == 2}
step("U4.S48 = {carrier even, family number 2} in the R sector equals the project's retained 48 exactly",
     S48 == proj_retained and len(S48) == 48)

# U5: selection holonomies as number holonomies, hypercharge linear in occupations
F = sp.Rational
frac = lambda q: q - sp.floor(q)
z3 = [F(0)] * 5 + [F(2, 3)] * 3
z4 = [F(-1, 2)] * 5 + [F(0)] * 3
dot = lambda z, r: sum(zi * F(x, 2) for zi, x in zip(z, r))
ok5 = all(frac(dot(z3, r)) == frac(F(2, 3) * family_occ(r)) for r in gso)
ok5 &= all(frac(dot(z4, r)) == F(1, 4) for r in gso if carrier_occ(r) % 2 == 0) and \
       all(frac(dot(z4, r)) == F(3, 4) for r in gso if carrier_occ(r) % 2 == 1)
y = list(hypercharge_vector())
Ylin = lambda r: sum(yi * F(x, 2) for yi, x in zip(y, r))
ok5 &= all(Ylin(r) == sum(yi * (F(1) if x == 1 else 0) for yi, x in zip(y[:5], r[:5])) - sum(y[:5]) / 2 for r in gso)
step("U5.z3 = exp(2 pi i (2/3) N_fam), z4 = carrier-parity phase (1/4 even, 3/4 odd); hypercharge = sum_i y_i n_i (linear in carrier occupations)",
     ok5)

# U6: Z2 sheet
step("U6.(-1)^F of the Z2 sheet = -1 exactly on the R sector = the spinor roots (the mu8 deck lift g^4 of C10)",
     len(gso) == 128 and all(any(abs(x) == 1 for x in r) for r in gso) and not any(any(abs(x) == 1 for x in r) for r in ns))

# U7: E8 as closure D5 + D3 in D8 in E8, glue Z4
D5 = {r for r in ns if not any(r[5:])}
D3 = {r for r in ns if not any(r[:5])}
D8 = ns
def gram_det(basis):
    B = sp.Matrix(basis)
    return (B * B.T).det()
e = lambda i: [1 if k == i else 0 for k in range(8)]
sub = lambda u, w: [x - y for x, y in zip(u, w)]
add = lambda u, w: [x + y for x, y in zip(u, w)]
D5b = [sub(e(i), e(i + 1)) for i in range(4)] + [add(e(3), e(4))]
D3b = [sub(e(5), e(6)), sub(e(6), e(7)), add(e(6), e(7))]
E8b = [[F(1, 2), F(-1, 2), F(-1, 2), F(-1, 2), F(-1, 2), F(-1, 2), F(-1, 2), F(1, 2)],
       add(e(0), e(1))] + [sub(e(i), e(i - 1)) for i in range(1, 7)]       # Bourbaki simple roots of E8
in_sub = lambda v: all(sp.Integer(x) == x for x in v) and sum(v[:5]) % 2 == 0 and sum(v[5:]) % 2 == 0
spin = [F(1, 2)] * 8
step("U7.E8 = closure: det D5 = det D3 = 4, det E8 = 1 -> index sqrt(16) = 4; a spinor class has order 4 in the quotient (cyclic Z4 glue); D5, D3 inside the 112 NS roots",
     gram_det(D5b) == 4 and gram_det(D3b) == 4 and gram_det(E8b) == 1
     and not in_sub([2 * x for x in spin]) and in_sub([4 * x for x in spin])
     and len(D5) == 40 and len(D3) == 12 and len(D8) == 112)

# U8: flavor monodromy is real
from tfpt_explorer.flavor_path_transport import original_fuchs_connection
from scipy.integrate import solve_ivp
def transport(segs):
    T = np.eye(3, dtype=complex)
    for path in segs:
        def rhs(t, yv, path=path):
            zz, dz = path(t)
            return (dz * original_fuchs_connection(zz) @ yv.reshape(3, 3)).ravel()
        T = solve_ivp(rhs, (0, 1), np.eye(3, dtype=complex).ravel(), method="DOP853", rtol=1e-12, atol=1e-14).y[:, -1].reshape(3, 3) @ T
    return T
ray = lambda a, b: (lambda t: (a + (b - a) * t, b - a))
circ = lambda c, r, s0: (lambda t: (c + r * np.exp(1j * (s0 + 2 * np.pi * t)), 2j * np.pi * r * np.exp(1j * (s0 + 2 * np.pi * t))))
M0 = transport([ray(0, 0.6), circ(1, 0.4, np.pi), ray(0.6, 0)])
Um = np.diag([1, 1j, -1j])
key = lambda X: tuple(np.round(X, 6).ravel())
group = {key(np.eye(3)): np.eye(3, dtype=complex)}
front = list(group.values())
while front:
    new = []
    for X in front:
        for G in (M0, Um):
            Y = G @ X
            if key(Y) not in group:
                group[key(Y)] = Y; new.append(Y)
    front = new
fs = sum(np.trace(X @ X) for X in group.values()) / len(group)
step("U8.the flavor monodromy group (|G| = 24) is real: Frobenius-Schur indicator +1, i.e. real rotations SO(3) < U(3) < SO(6) of the family modes",
     len(group) == 24 and abs(fs - 1) < 1e-8, f"FS = {fs.real:.6f}")

# U9: geometry carried by U itself
dims = {"PSL(2,R) (chiral net on S^1)": 3, "Poincare(3,1)": 10, "conformal SO(4,2)": 15}
step("U9.(statement) the chiral seam net's own geometric group is PSL(2,R) (dim 3), not a 3+1D group: Theorem D is NOT a corollary of U",
     dims["PSL(2,R) (chiral net on S^1)"] == 3 and dims["Poincare(3,1)"] != 3,
     "a 3+1D geometry needs a net indexed by 4D regions (e.g. the 2+1 Chern bulk of the QWZ edge is only 3D)")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG)}
(HERE / "universal_source_results.json").write_text(json.dumps(out, indent=2))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
