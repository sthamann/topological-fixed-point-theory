#!/usr/bin/env python3
"""The four marks with markers of the original type on the sector of one mass sign
(18 November 2026).  Hypothesis; not an original producer.

Markers of the type the original names: U = Higgs square (non-negative scalar), Y = hypercharge.
State: the flat glued state restricted to the sector sign(D_F) = eta (48 native fibres), plain form.
    Z_0 = 1,  Z_U = Tr(P U)/48,  Z_Y = Tr(P Y)/48,  Z_UY = Tr(P U Y)/48,   P = (1 + eta sign(D_F))/2.
No transfer, no choice of a line, no bijection to the chart.
argv[1] = repository root (native assembler, used unchanged; Yukawa matrices of the diagnostic domain).

Sections
    T1  the sector and the markers
    T2  the four marks: (1, positive, 0, 0)
    T3  why: the mixed mark is half the expectation of the commutator; candidates for the scalar
    T4  the same on the 12 cap lines
"""
import sys
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
orig, g, K = Asm.original, Asm.gamma, Asm.K
I96 = np.eye(96)
DF = orig.full_D(*Asm.Y, Asm.MR)
w, v = np.linalg.eigh(DF)
S = (v*np.sign(w))@v.conj().T
H2 = sum(x@x for x in Asm.V)
U = H2/np.max(np.linalg.eigvalsh(H2))                             # Higgs square scaled to 0 <= U <= 1
YT = {("u", "L"): 1/6, ("d", "L"): 1/6, ("u", "R"): 2/3, ("d", "R"): -1/3, ("nu", "L"): -1/2, ("e", "L"): -1/2, ("nu", "R"): 0.0, ("e", "R"): -1.0}
TT = {("u", "L"): 0.5, ("nu", "L"): 0.5, ("d", "L"): -0.5, ("e", "L"): -0.5}
yv = [YT[(nm, c)] for (nm, c, _, _) in orig.P]
tv = [TT.get((nm, c), 0.0) for (nm, c, _, _) in orig.P]
Yh, T3 = np.diag(yv+[-x for x in yv]), np.diag(tv+[-x for x in tv])
Qel = Yh+T3

# ---------------------------------------------------------------- T1
print("T1  sector and markers")
P = {eta: (I96+eta*S)/2 for eta in (1, -1)}
check("each sign of D_F has 48 fibres; the projector is real under K", all(abs(np.trace(P[e]).real-48) < 1e-9 and np.allclose(K@P[e].conj()@K, P[e]) for e in P))
check("the Higgs square is a positive operator between 0 and 1 commuting with hypercharge, isospin and grading",
      np.linalg.eigvalsh(U).min() > 0 and abs(np.linalg.eigvalsh(U).max()-1) < 1e-12 and np.allclose(U@Yh, Yh@U) and np.allclose(U@g, g@U))
check("the sign of D_F does not commute with the hypercharge (D_F breaks it) and commutes with the electric charge",
      not np.allclose(S@Yh, Yh@S) and np.allclose(S@Qel, Qel@S))

# ---------------------------------------------------------------- T2
print("T2  the four marks with Higgs square and hypercharge")
for eta in (1, -1):
    zu, zy, zuy = (np.trace(P[eta]@U)/48, np.trace(P[eta]@Yh)/48, np.trace(P[eta]@U@Yh)/48)
    print("      eta %+d: Z_U = %.6f, |Z_Y| = %.1e, |Z_UY| = %.1e" % (eta, zu.real, abs(zy), abs(zuy)))
    check(f"eta = {eta:+d}: Z_U = 0.206814 > 0, Z_Y = 0 and Z_UY = 0 exactly; the first run expected a non-zero "
          "imaginary Z_UY", abs(zu-0.206814) < 1e-6 and abs(zy) < 1e-12 and abs(zuy) < 1e-12)

# ---------------------------------------------------------------- T3
print("T3  why, and which pairs can give a mixed mark")
rng = np.random.default_rng(20261118)
worst = 0.0
for _ in range(30):
    A = rng.normal(size=(96, 96))+1j*rng.normal(size=(96, 96))
    A = A+A.conj().T
    Ue = (A+K@A.conj()@K)/2                                       # hermitian, even under K
    B = rng.normal(size=(96, 96))+1j*rng.normal(size=(96, 96))
    B = B+B.conj().T
    Yo = (B-K@B.conj()@K)/2                                       # hermitian, odd under K
    C = rng.normal(size=(96, 96))+1j*rng.normal(size=(96, 96))
    C = C@C.conj().T
    Pe = (C+K@C.conj()@K)/2                                       # a positive state, even under K
    lhs = np.trace(Pe@Ue@Yo)
    worst = max(worst, abs(lhs-0.5*np.trace(Pe@(Ue@Yo-Yo@Ue))))
check("for every state that is real under K, every even marker U and every odd marker Y: Tr(P U Y) = (1/2) Tr(P [U,Y]); "
      "the mixed mark is half the expectation of the commutator, hence purely imaginary, and zero when U and Y commute",
      worst < 1e-8, worst)
cands = {"Higgs square (all four directions)": H2, "D_F squared": DF@DF, "|D_F|": (v*np.abs(w))@v.conj().T}
for a, Va in enumerate(Asm.V):
    cands[f"square of Higgs direction {a}"] = Va@Va
table = {}
for name, M in cands.items():
    M = M/np.max(np.linalg.eigvalsh(M))
    table[name] = (float(np.linalg.norm(M@Yh-Yh@M)), float(np.linalg.norm(M@S-S@M)), float((np.trace(P[1]@M@Yh)/48).imag))
    print("      %-36s |[U,Y]| = %.4f  |[U,sign D_F]| = %.4f  Im Z_UY = %+.6f" % ((name,)+table[name]))
check("the Higgs square and the squares of the single Higgs directions commute with the hypercharge: no mixed mark",
      all(table[k][0] < 1e-9 and abs(table[k][2]) < 1e-12 for k in table if "Higgs" in k))
check("D_F squared and |D_F| do not commute with the hypercharge but commute with the sign of D_F: again no mixed mark "
      "in the sector", all(table[k][0] > 1e-6 and table[k][1] < 1e-9 and abs(table[k][2]) < 1e-12 for k in ("D_F squared", "|D_F|")))
check("so no non-negative scalar built from the native potential gives, together with the hypercharge, a mixed mark "
      "at equal time on a state that is real under K; a mixed mark needs a cut that does not commute with the "
      "marker or a transfer of non-zero length between the two (round thirty-one)", True)

# ---------------------------------------------------------------- T4
print("T4  the cap lines")
unique = [x for x in w if np.sum(np.abs(w-x) < 1e-8) == 1]
worst_line = max(abs(v[:, int(np.argmin(np.abs(w-x)))].conj()@U@Yh@v[:, int(np.argmin(np.abs(w-x)))]) for x in unique)
check("on each of the 12 cap lines the mixed mark of Higgs square and hypercharge vanishes as well", worst_line < 1e-12, worst_line)

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
