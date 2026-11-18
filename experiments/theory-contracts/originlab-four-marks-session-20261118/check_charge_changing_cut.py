#!/usr/bin/env python3
"""Second branch of round forty-four: a cut that changes the hypercharge instead of measuring it
(18 November 2026).

State: flat glued state on the sector sign(D_F) = eta (48 native fibres), plain form, equal time.
Marks for a marker c and a cut X:  Z_c = Tr(P c)/48,  Z_X = Tr(P X)/48,  Z_cX = Tr(P X c)/48.
argv[1] = repository root (native assembler, used unchanged).

Sections
    B1  cut = direction of the native potential: the marks are (1, <c>, 0, i eta <gamma c>)
    B2  the chiral averages <gamma c> of the recorded charges: zero for hypercharge, its cube,
        isospin and electric charge
    B3  cut = carried pair currents, marker = hypercharge
"""
import itertools
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
orig, g, K = Asm.original, Asm.gamma, Asm.K
I96 = np.eye(96)
DF = orig.full_D(*Asm.Y, Asm.MR)
w, v = np.linalg.eigh(DF)
S = (v*np.sign(w))@v.conj().T
wp, vp = np.linalg.eigh(-1j*g@DF)
Ypot = (vp*np.sign(wp))@vp.conj().T
YT = {("u", "L"): 1/6, ("d", "L"): 1/6, ("u", "R"): 2/3, ("d", "R"): -1/3, ("nu", "L"): -1/2, ("e", "L"): -1/2, ("nu", "R"): 0.0, ("e", "R"): -1.0}
TT = {("u", "L"): 0.5, ("nu", "L"): 0.5, ("d", "L"): -0.5, ("e", "L"): -0.5}
yv = [YT[(nm, c)] for (nm, c, _, _) in orig.P]
tv = [TT.get((nm, c), 0.0) for (nm, c, _, _) in orig.P]
Yh, T3 = np.diag(yv+[-x for x in yv]), np.diag(tv+[-x for x in tv])
Npart = np.diag([1.0]*48+[-1.0]*48)
P = {eta: (I96+eta*S)/2 for eta in (1, -1)}
rat = lambda x: Fr(float(x)).limit_denominator(1000)

# ---------------------------------------------------------------- B1
print("B1  cut = direction of the native potential")
check("the direction of the potential changes the hypercharge: it does not commute with it; it commutes with the "
      "electric charge", not np.allclose(Ypot@Yh, Yh@Ypot) and np.allclose(Ypot@(Yh+T3), (Yh+T3)@Ypot))
markers = {"grading": g, "hypercharge": Yh, "isospin": T3, "electric charge": Yh+T3, "particle minus conjugate": Npart,
           "hypercharge squared": Yh@Yh, "hypercharge cubed": Yh@Yh@Yh, "sign of hypercharge": np.sign(Yh),
           "lift (1 + grading)/2": (I96+g)/2, "projector hypercharge > 0": (np.sign(Yh)+np.abs(np.sign(Yh)))/2}
ok, table = True, {}
for name, c in markers.items():
    for eta in (1, -1):
        zc, zx, zcx = np.trace(P[eta]@c)/48, np.trace(P[eta]@Ypot)/48, np.trace(P[eta]@Ypot@c)/48
        ok = ok and abs(zc-np.trace(c)/96) < 1e-10 and abs(zx) < 1e-10 and abs(zcx-1j*eta*np.trace(g@c)/96) < 1e-10
    table[name] = (rat(np.trace(c).real/96), rat(np.trace(g@c).real/96))
print("      marker: (<c>, <gamma c>):", {k: (str(a), str(b)) for k, (a, b) in table.items()})
check("for every marker c that commutes with the grading: Z_c = <c>, Z_X = 0 and Z_cX = i eta <gamma c>, with <.> the "
      "average over the 96 fibres: the mixed mark is the chiral average of the marker", ok)

# ---------------------------------------------------------------- B2
print("B2  the chiral averages")
check("grading 1, lift 1/2 (the recorded canonical values)", table["grading"][1] == 1 and table["lift (1 + grading)/2"][1] == Fr(1, 2))
check("hypercharge, hypercharge cubed, isospin, electric charge, particle number: 0. For the hypercharge and its cube "
      "these are the conditions that the left and right sums agree (cancellation of the anomalies of the recorded "
      "charge table with a right-handed neutrino)",
      all(table[k][1] == 0 for k in ("hypercharge", "hypercharge cubed", "isospin", "electric charge", "particle minus conjugate")))
check("every marker that is the same for a field and its conjugate (hypercharge squared, and round forty-four: the "
      "Higgs square) has chiral average 0", table["hypercharge squared"][1] == 0)
check("non-zero only for markers correlated with the chirality: sign of the hypercharge 5/16, projector on positive "
      "hypercharge 5/32", table["sign of hypercharge"][1] == Fr(5, 16) and table["projector hypercharge > 0"][1] == Fr(5, 32))
check("so with the native cut that changes the hypercharge, the hypercharge marker gives no mixed mark at equal time: "
      "the recorded charge table forbids it", True)

# ---------------------------------------------------------------- B3
print("B3  cut = carried pair currents, marker = hypercharge")
bits = np.array(list(itertools.product((0, 1), repeat=8)))
ann = []
for i in range(8):
    a = np.zeros((256, 256))
    for k, b in enumerate(bits):
        if b[i]:
            c = b.copy()
            c[i] = 0
            a[int("".join(map(str, c)), 2), k] = (-1)**int(b[:i].sum())
    ann.append(a)
spin2 = 2*bits-1
matter = np.array([list(x).count(-1) % 2 == 0 for x in spin2])
even_fam = sorted({tuple(x[5:]) for x in spin2 if list(x[5:]).count(-1) % 2 == 0})
kept = even_fam[1:]
idx = np.where(matter & np.array([tuple(x[5:]) not in (even_fam[0], tuple(-y for y in even_fam[0])) for x in spin2]))[0]
pos = {int(k): n for n, k in enumerate(idx)}
Phi = np.zeros((96, 96))
for k in idx:
    s = spin2[k]
    particle = int(np.sum(s[:3] > 0)) in (1, 3)
    t = s if particle else -s
    quark = int(np.sum(t[:3] > 0)) == 1
    name = ("u" if t[3] > 0 else "d") if quark else ("nu" if t[3] > 0 else "e")
    fam = tuple(int(x) for x in t[5:])
    fam = fam if fam in kept else tuple(-x for x in fam)
    Phi[(0 if particle else 48)+orig.IDX[(name, "L" if t[3] != t[4] else "R", int(np.argmax(t[:3] > 0)) if quark else 0, kept.index(fam))], pos[int(k)]] = 1
rows = []
for i in range(5):
    for j in range(5, 8):
        for sj in (1, -1):
            X = Phi@((ann[i].T)@(ann[j].T if sj > 0 else ann[j]))[np.ix_(idx, idx)]@Phi.T
            for cut in ((X+X.T)/np.sqrt(2), 1j*(X-X.T)/np.sqrt(2)):
                rows.append(("weak" if i >= 3 else "colour", np.trace(P[1]@cut)/48, np.trace(P[1]@cut@Yh)/48, np.trace(cut@Yh)/96))
zx = max(abs(r[1]) for r in rows)
zcx = np.array([r[2] for r in rows])
check("on the full flat state every such mixed mark vanishes (the currents change the hypercharge)", max(abs(r[3]) for r in rows) < 1e-12)
print("      sector eta = +1: largest |Z_X| = %.6f; mixed marks: largest real part %.6f, largest imaginary part %.6f"
      % (zx, np.max(np.abs(zcx.real)), np.max(np.abs(zcx.imag))))
check("in the sector of one mass sign the single mark of a carried current and its mixed mark with the hypercharge are "
      "not zero in general; they are small and depend on the current (and on the order of the generations, round "
      "thirty-eight)", np.max(np.abs(zcx)) > 1e-6 and np.max(np.abs(zcx)) < 0.1)

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
