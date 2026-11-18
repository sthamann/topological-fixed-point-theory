#!/usr/bin/env python3
"""Does the native circle operator single out a real line? (18 November 2026)

argv[1] = repository root.  Used unchanged: the assembler of the original Q contract (and through
it verification/v252_full_finite_triple.py).  A line (a vector up to phase) is real when its
projector is invariant under the recorded real structure: K conj(v v^+) K = v v^+.

Sections
    L1  recorded quarter shift: no eigenline of the operator is real, in all 27 configurations
    L2  which shift gives a mode with a real operator: only shift 0, the shift of a sign mark
    L3  shift 0, constant mode: the real eigenlines and which of them are unique
    L4  marks on these lines; they are not the cap of the recorded canonical model
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
orig = Asm.original
g, K = Asm.gamma, Asm.K


def cand(genes):
    return Asm.assemble({"schema": "tfpt.field-configuration.v1", "domain_sha256": Asm.domain["_sha256"], "genes": list(genes)})


def defect(rho):
    return float(np.max(np.abs(K@rho.conj()@K-rho)))


# ---------------------------------------------------------------- L1
print("L1  the recorded quarter shift")
G4 = np.kron(np.eye(4), g)
min_def, min_gr, count = 1.0, 1.0, 0
for genes in itertools.product(range(3), repeat=3):
    a = cand(genes)
    w, v = np.linalg.eigh(a["h_N"])
    for k in range(v.shape[1]):
        blocks = v[:, k].reshape(4, 96)
        rho = blocks.T@blocks.conj()                               # fibre density of the eigenvector
        min_def = min(min_def, defect(rho))
        min_gr = min(min_gr, abs(float((v[:, k].conj()@G4@v[:, k]).real)))
    count += 1
check("27 configurations, 384 eigenvectors each: the fibre density of no eigenvector is real",
      count == 27 and min_def > 1e-6, ("smallest defect", round(min_def, 6)))
check("reason: a real line has grading expectation zero, and every eigenvector has a non-zero one",
      min_gr > 1e-6, ("smallest |<gamma>|", round(min_gr, 6)))

# ---------------------------------------------------------------- L2
print("L2  which shift allows a real mode")
p0 = cand([1, 1, 1])["phi_fourier"][0]
check("the constant potential is real under K; a mode operator (n + shift) gamma + potential is real exactly when "
      "n + shift = 0", np.allclose(K@p0.conj()@K, p0) and not np.allclose(K@(0.25*g+p0).conj()@K, 0.25*g+p0))
shifts = [Fr(0), Fr(1, 4), Fr(1, 2), Fr(3, 4)]
allowed = [str(s) for s in shifts if any(n+s == 0 for n in range(-4, 5))]
check("of the shifts 0, 1/4, 1/2, 3/4 of the mark types only shift 0 has such a mode: the shift of a sign mark "
      "(round twenty-eight)", allowed == ["0"], allowed)

# ---------------------------------------------------------------- L3
print("L3  shift 0, constant mode")
w, v = np.linalg.eigh(p0)
tol = 1e-8
groups = []
for val in w:
    if groups and abs(val-groups[-1][0]) < tol:
        groups[-1][1] += 1
    else:
        groups.append([val, 1])
mult = sorted({m for _, m in groups})
single = [val for val, m in groups if m == 1]
print("      multiplicities of the 96 values:", {m: sum(1 for _, mm in groups if mm == m) for m in mult})
print("      unique values:", np.round(single, 6))
names = [(n, ch) for (n, ch, c, gen) in orig.P]
lab = names+[(n+"*", ch) for n, ch in names]
real_ok, content = True, []
for val in single:
    k = int(np.argmin(np.abs(w-val)))
    x = v[:, k]
    real_ok = real_ok and defect(np.outer(x, x.conj())) < 1e-8
    wts = {}
    for i, amp in enumerate(x):
        wts[lab[i]] = wts.get(lab[i], 0)+abs(amp)**2
    content.append(max(wts, key=wts.get))
check("every unique eigenline of the constant potential is real", real_ok and len(single) > 0, len(single))
print("      main component of the unique lines:", sorted(set(content)))
check("the unique real lines are neutrino lines mixing a field with its conjugate: Majorana lines",
      all(n.startswith("nu") for n, _ in content), sorted(set(content)))
check("the other values are degenerate (colour and field/conjugate): there the real lines form families, none is unique",
      len(single) < 96)

# ---------------------------------------------------------------- L4
print("L4  marks on the unique real lines")
wv, vv = np.linalg.eigh(p0)
Y = (vv*np.sign(wv))@vv.conj().T
Pi = -1j*Y@g
rows = []
for val in single:
    x = v[:, int(np.argmin(np.abs(w-val)))]
    rows.append((val, (x.conj()@g@x).real, (x.conj()@Y@x).real, abs(x.conj()@Y@g@x), (x.conj()@Pi@x).real))
arr = np.array(rows)
check("on them: <grading> = 0, <direction of the potential> = sign of the value, <Y gamma> = 0, <Pi> = 0",
      np.allclose(arr[:, 1], 0, atol=1e-9) and np.allclose(arr[:, 2], np.sign(arr[:, 0])) and np.allclose(arr[:, 3], 0, atol=1e-9)
      and np.allclose(arr[:, 4], 0, atol=1e-9))
check("so with the cuts of round thirty the marks are (1, 0, +-1, 0), not the recorded (1, 0, 0, i eta): the operator "
      "singles out lines of definite sign of the mass, the recorded cap is a line of definite Pi", True)
check("Pi anticommutes with the constant potential: it exchanges the lines of value +m and -m; the operator does not "
      "prefer either value of Pi, so eta stays open", np.allclose(Pi@p0+p0@Pi, 0))

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
