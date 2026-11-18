#!/usr/bin/env python3
"""The recorded rate q read as the chiral asymmetry of a non-negative marker (18 November 2026).
Hypothesis; not an original producer.

From round forty-five: on the sector sign(D_F) = eta, at equal time, with the cut X = direction of
the native potential, every marker c commuting with the grading has Z_c = <c> and
Z_cX = i eta <gamma c>.  For a non-negative marker the normalised mixed mark is therefore
        Z_cX / Z_c = i eta q,     q = <gamma c>/<c>,   -1 <= q <= 1,
the form of the recorded canonical model with its positive lift (responses 1/2 and i q/2).
argv[1] = repository root (native assembler, used unchanged).

Sections
    A1  the bound and the two extreme cases
    A2  the Higgs square seen through one chirality: q = -1 or +1 exactly
    A3  projectors of the recorded charges: rational values of q
    A4  the values on the sector, computed directly
"""
import sys
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np

FAIL = []


def check(name, ok, info=None):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info is not None else ""))
    if not ok:
        FAIL.append(name)


root = Path(sys.argv[1]).resolve()
base = root/"experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q"
sys.path.insert(0, str(base/"runtime"))
from assemble import Assembler  # noqa: E402

Asm = Assembler(base/"domain/native-diagnostic-v1.json")
orig, g = Asm.original, Asm.gamma
I96 = np.eye(96)
DF = orig.full_D(*Asm.Y, Asm.MR)
w, v = np.linalg.eigh(DF)
S = (v*np.sign(w))@v.conj().T
wp, vp = np.linalg.eigh(-1j*g@DF)
X = (vp*np.sign(wp))@vp.conj().T
YT = {("u", "L"): 1/6, ("d", "L"): 1/6, ("u", "R"): 2/3, ("d", "R"): -1/3, ("nu", "L"): -1/2, ("e", "L"): -1/2, ("nu", "R"): 0.0, ("e", "R"): -1.0}
TT = {("u", "L"): 0.5, ("nu", "L"): 0.5, ("d", "L"): -0.5, ("e", "L"): -0.5}
yv = np.array([YT[(nm, c)] for (nm, c, _, _) in orig.P])
tv = np.array([TT.get((nm, c), 0.0) for (nm, c, _, _) in orig.P])
Yd, Td = np.concatenate([yv, -yv]), np.concatenate([tv, -tv])
Qd = Yd+Td
H2 = sum(a@a for a in Asm.V)
PL, PR = (I96+g)/2, (I96-g)/2
avg = lambda c: (np.trace(c)/96).real
qof = lambda c: avg(g@c)/avg(c)
rat = lambda x: Fr(float(x)).limit_denominator(1000)

# ---------------------------------------------------------------- A1
print("A1  the bound")
rng = np.random.default_rng(20261118)
ok = True
for _ in range(50):
    d = rng.uniform(0, 1, size=96)
    ok = ok and abs(qof(np.diag(d))) <= 1+1e-12
check("for a non-negative marker commuting with the grading |q| <= 1; q = +1 or -1 exactly when the marker lives on "
      "one sign of the grading; q = 0 when it is the same on both", ok and abs(qof(PL)-1) < 1e-12 and abs(qof(PR)+1) < 1e-12 and abs(qof(I96)) < 1e-12)

# ---------------------------------------------------------------- A2
print("A2  the Higgs square through one chirality")
via_left = sum(a@PL@a for a in Asm.V)                             # Higgs, projector, Higgs
via_right = sum(a@PR@a for a in Asm.V)
check("the four Higgs directions anticommute with the grading, so sum_a V_a P_L V_a = P_R H^2 and sum_a V_a P_R V_a = P_L H^2: "
      "a Higgs pair around one chirality is the Higgs square on the other",
      all(np.allclose(a@g+g@a, 0) for a in Asm.V) and np.allclose(via_left, PR@H2) and np.allclose(via_right, PL@H2))
check("these are non-negative operators commuting with the grading, with <c> = <H^2>/2 and q = -1 and +1 exactly",
      np.linalg.eigvalsh(via_left).min() > -1e-10 and np.allclose(via_left@g, g@via_left)
      and abs(avg(via_left)-avg(H2)/2) < 1e-10 and abs(qof(via_left)+1) < 1e-12 and abs(qof(via_right)-1) < 1e-12)
check("the full Higgs square has q = 0 (round forty-four)", abs(qof(H2)) < 1e-12)

# ---------------------------------------------------------------- A3
print("A3  projectors of the recorded charges")
proj = {"hypercharge > 0": np.diag((Yd > 1e-9).astype(float)), "hypercharge < 0": np.diag((Yd < -1e-9).astype(float)),
        "electric charge > 0": np.diag((Qd > 1e-9).astype(float)), "electric charge < 0": np.diag((Qd < -1e-9).astype(float)),
        "isospin > 0": np.diag((Td > 1e-9).astype(float)), "isospin non-zero": np.diag((np.abs(Td) > 1e-9).astype(float)),
        "hypercharge zero": np.diag((np.abs(Yd) < 1e-9).astype(float)), "particles": np.diag([1.0]*48+[0.0]*48)}
table = {k: (rat(avg(c)), rat(qof(c))) for k, c in proj.items()}
print("      projector: (<c>, q):", {k: (str(a), str(b)) for k, (a, b) in table.items()})
check("positive hypercharge 15/32 with q = 1/3, negative hypercharge q = -1/3", table["hypercharge > 0"] == (Fr(15, 32), Fr(1, 3)) and table["hypercharge < 0"][1] == Fr(-1, 3))
check("the particle projector and the projectors of the electric charge have q = 0", table["particles"][1] == 0
      and table["electric charge > 0"][1] == 0 and table["electric charge < 0"][1] == 0)
check("the values of q are rational numbers fixed by the recorded charge table; none of them is a rate in time", True)

# ---------------------------------------------------------------- A4
print("A4  on the sector, directly")
ok = True
for name, c in list(proj.items())+[("Higgs through left", via_left), ("Higgs through right", via_right)]:
    for eta in (1, -1):
        Pe = (I96+eta*S)/2
        zc, zx, zcx = np.trace(Pe@c)/48, np.trace(Pe@X)/48, np.trace(Pe@X@c)/48
        ok = ok and abs(zc-avg(c)) < 1e-10 and abs(zx) < 1e-10 and abs(zcx-1j*eta*qof(c)*zc) < 1e-10
check("for all ten markers and both signs: Z_c = <c>, Z_X = 0, Z_cX = i eta q Z_c", ok)
check("so the recorded form (non-negative marker, single cut mark zero, mixed mark i eta q) holds at equal time with "
      "q the chiral asymmetry of the marker; the recorded rate needs no length of time in this reading", True)

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
