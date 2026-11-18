#!/usr/bin/env python3
"""Reconcile the user's unrestricted lift/cover tests with the fixed originals."""
import contextlib
import io
import json
from math import gcd, prod
from pathlib import Path
import runpy

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

HERE = Path(__file__).resolve().parent
with contextlib.redirect_stdout(io.StringIO()):
    original = runpy.run_path(str(HERE.parent / "e-series-joint-source-20261023/joint_clock.py"))
    relative = runpy.run_path(str(HERE.parent / "e-series-joint-source-20261023/relative_transport.py"))


def need(value, label):
    if not value:
        raise RuntimeError(label)


J, c = original["J"], original["c"]
power_phases = original["power_basis_phases"]
dual = original["dual_action"]
rank = original["gf2_rank"]
identity = sp.eye(8)
snf = smith_normal_form(identity-J, domain=sp.ZZ)
invariants = [abs(int(snf[i, i])) for i in range(8)]
need(invariants == [1, 1, 1, 1, 2, 2, 2, 2], "actual E8 lattice SNF(I-J)")

def norm(matrix, order):
    return sum((matrix**k for k in range(order)), sp.zeros(8))

need(norm(J, 4) == sp.zeros(8) and power_phases(J, 0, 4) == (1,)*8,
     "all 256 sign lifts of J have order four")
need(norm(c, 12) == sp.zeros(8) and power_phases(c, 0, 12) == (1,)*8,
     "all 256 sign lifts of c have order twelve")
j_gauge = {a ^ dual(J, a) for a in range(256)}
c_gauge = {a ^ dual(c, a) for a in range(256)}
need(len(j_gauge) == 16 and len(c_gauge) == 64, "unrestricted sign gauge orbit sizes")
need(rank([dual(J, 1 << i) ^ (1 << i) for i in range(8)]) == 4,
     "unrestricted J character coboundary rank")
need(rank([dual(c, 1 << i) ^ (1 << i) for i in range(8)]) == 6,
     "unrestricted c character coboundary rank")
fixed = original["result"]["fixed_Jhat_lifts"]
need(fixed["count"] == 4 and fixed["orbit_count"] == 1, "original constrained lift problem")
# SNF diagonal entries 1 and 2 imply every sign-valued character is an
# (I-J*) coboundary of a fourth-root-valued character. This is an exact
# lattice statement, not a change to the sign-gauge equivalence above.
need(all(4 % d == 0 and d in (1, 2) for d in invariants), "quarter-phase gauge suffices")
need(prod(gcd(d, 4) for d in invariants) == 16,
     "all mu4 candidate lifts modulo mu4 gauge instead have sixteen classes")

A = relative["A"]
winding = [1, 1, 1, 1, -4]
native, burau = [A**(9*k) for k in winding], [A**(4*k) for k in winding]
rel = [A**((7*k) % 12) for k in winding]
for g, b, r in zip(native, burau, rel):
    for H in relative["basis"]:
        need(b * H * g.inv() == r * H, "actual relative action of EACH pi1 generator")
need(sp.prod(rel) == sp.eye(4), "full punctured-sphere product relation")
need(len({tuple(A**(7*k)) for k in range(12)}) == 12,
     "global relative image has exactly twelve elements")
need({k for k in range(12) if A**(7*k) == sp.eye(4)} == {0}, "faithful winding factor")

# User's exact negative control: independent order-twelve peripheral images
# can generate C12 x C12 on a FOUR-punctured sphere. Local order is inadequate.
counter_generators = [(1, 0), (0, 1), (0, 0), (11, 11)]
counter_image = {(0, 0)}
while True:
    enlarged = counter_image | {((x+a) % 12, (y+b) % 12)
                               for x, y in counter_image for a, b in counter_generators}
    if enlarged == counter_image:
        break
    counter_image = enlarged
need(len(counter_image) == 144, "negative control distinguishes local orders from global image")
need(tuple(sum(a[i] for a in counter_generators) % 12 for i in range(2)) == (0, 0),
     "negative-control sphere relation")

result = {
    "research_id": "MATH.SEAM.ESERIES.INCIDENCE.07.CLOCK-AUDIT", "verdict": "PASS",
    "SNF_I_minus_J": invariants,
    "unrestricted_J_sign_lifts": 256, "J_sign_gauge_classes": 16,
    "J_classes_after_quarter_phase_gauge": 1,
    "quarter_phase_gauge_scope": "The 256 SIGN-lift candidates remain fixed; only their gauge equivalence is enlarged.",
    "all_mu4_J_lift_candidates_modulo_mu4_gauge_classes": 16,
    "unrestricted_c_sign_lifts": 256, "c_sign_gauge_classes": 4,
    "fixed_Jhat_c_lifts": fixed,
    "common_clock_matrix": [list(map(int, A.row(i))) for i in range(4)],
    "native_peripheral_exponents_mod12": [(9*k) % 12 for k in winding],
    "burau_peripheral_exponents_mod12": [(4*k) % 12 for k in winding],
    "relative_peripheral_exponents_mod12": [(7*k) % 12 for k in winding],
    "pi1": "<g1,g2,g3,g4,gInfinity | g1*g2*g3*g4*gInfinity=1>",
    "global_factorization": "rho_relative(word)=A^(7*nu(word)); nu(g1..g4)=1, nu(gInfinity)=-4",
    "global_image_order": 12, "minimal_cover": "C12: s^12=x^4-1",
    "local_order_only_negative_control_image_order": 144,
    "scope": "Original fixed sign cocycle and marked Jhat; Burau system induced from the specified cover and deck representation, not an independently identified Gauss-Manin family.",
}
print(json.dumps(result, indent=2))
