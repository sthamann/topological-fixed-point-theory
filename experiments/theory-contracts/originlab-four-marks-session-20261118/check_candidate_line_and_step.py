#!/usr/bin/env python3
"""Do the line and the step of the candidate H-D follow from the record? (18 November 2026)

argv[1] = repository root.  The assembler of the original Q contract is used unchanged.
Candidate as in round thirty-four: U = grading, Y = direction of the native potential,
transfer exp(-s |h_N|) in three equal steps of length s.

Sections
    S1  the line: at the baseline the cap spans exactly the lowest level of |h_N|
    S2  zero step: the marks are (1, 0, 0, i eta) in all 27 configurations
    S3  finite step: the marks depend on the step; the recorded 4 pi does not single out a value
    S4  other modes and other lines
"""
import itertools
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
g = Asm.gamma
I4 = np.eye(4)


def sign_of(X):
    w, v = np.linalg.eigh(X)
    return (v*np.sign(w))@v.conj().T


class Conf:
    def __init__(self, genes):
        a = Asm.assemble({"schema": "tfpt.field-configuration.v1", "domain_sha256": Asm.domain["_sha256"], "genes": list(genes)})
        self.modes, self.h, self.DF = a["modes"], a["h_N"], a["phi_fourier"][0]
        self.w, self.v = np.linalg.eigh(self.DF)
        self.U, self.Y = np.kron(I4, g), np.kron(I4, sign_of(-1j*g@self.DF))
        self.wh, self.vh = np.linalg.eigh(self.h)
        self.unique = [x for x in self.w if np.sum(np.abs(self.w-x) < 1e-8) == 1]

    def cap(self, n, x):
        c = np.zeros(384, complex)
        i = self.modes.index(n)
        c[96*i:96*(i+1)] = self.v[:, int(np.argmin(np.abs(self.w-x)))]
        return c

    def marks(self, cap, s):
        T = (self.vh*np.exp(-s*np.abs(self.wh)))@self.vh.conj().T
        z0 = cap.conj()@T@T@T@cap
        return cap.conj()@T@T@self.U@T@cap/z0, cap.conj()@T@self.Y@T@T@cap/z0, cap.conj()@T@self.Y@T@self.U@T@cap/z0


# ---------------------------------------------------------------- S1
print("S1  the line")
b = Conf((1, 1, 1))
m = min(abs(x) for x in b.unique)
E = np.sort(np.abs(b.wh))
check("baseline: the smallest |value| of D_F overall is that of the lightest unique line", abs(np.min(np.abs(b.w))-m) < 1e-12, round(m, 6))
check("baseline: the lowest level of |h_N| is two-fold with value sqrt(1/16 + m^2) = 0.264108, separated from the next (0.282879)",
      abs(E[0]-np.sqrt(1/16+m*m)) < 1e-10 and abs(E[1]-E[0]) < 1e-10 and E[2]-E[1] > 0.01, np.round(E[:3], 6))
idx = np.argsort(np.abs(b.wh))[:2]
P = b.vh[:, idx]@b.vh[:, idx].conj().T
ov = [float((b.cap(0, s*m).conj()@P@b.cap(0, s*m)).real) for s in (1, -1)]
check("the two cap lines (mode 0, mass value +m and -m) lie in this level and span it: the line of the candidate is "
      "the lowest level of the native operator", all(abs(o-1) < 1e-10 for o in ov), ov)
c = Conf((0, 0, 0))
Ec = np.sort(np.abs(c.wh))
print("      configuration (0,0,0): lowest levels", np.round(Ec[:3], 6))
check("with a non-constant profile the lowest level splits (0.254195 and 0.273337); the cap keeps 0.9896 of its weight "
      "in the two lowest levels", abs(Ec[0]-0.254195) < 1e-6 and Ec[1]-Ec[0] > 0.01,
      round(float((c.cap(0, m).conj()@(c.vh[:, np.argsort(np.abs(c.wh))[:2]]@c.vh[:, np.argsort(np.abs(c.wh))[:2]].conj().T)@c.cap(0, m)).real), 6))

# ---------------------------------------------------------------- S2
print("S2  zero step")
worst = 0.0
for genes in itertools.product(range(3), repeat=3):
    cf = Conf(genes)
    mm = min(abs(x) for x in cf.unique)
    for sg in (1, -1):
        zu, zy, zuy = cf.marks(cf.cap(0, sg*mm), 0.0)
        worst = max(worst, abs(zu), abs(zy), abs(zuy-1j*sg))
check("at step zero the marks are exactly (1, 0, 0, i eta) in all 27 configurations, whatever the profile", worst < 1e-10, worst)

# ---------------------------------------------------------------- S3
print("S3  finite step, configuration (0,0,0), lightest line")
steps = [0.25, 0.5, 1.0, 4*np.pi/3, 2*np.pi, 4*np.pi, 40.0]
rows = [(s,)+c.marks(c.cap(0, m), s) for s in steps]
for s, zu, zy, zuy in rows:
    print("      step %7.3f: q = %.6f, |Z_U/Z_0| = %.4f, |Z_Y/Z_0| = %.4f" % (s, zuy.imag, abs(zu), abs(zy)))
qs = [r[3].imag for r in rows]
check("the rate falls monotonically with the step (0.9945 at 0.25, 0.9884 at 1, 0.612 at 40) and the single marks grow: "
      "there is no plateau", all(x > y for x, y in zip(qs, qs[1:])) and abs(rows[-1][1]) > 0.5)
check("the recorded inverse temperature 4 pi gives 0.9811 when read as the total of three steps and 0.9335 when read as "
      "one step: the record does not say which, and neither value is distinguished",
      abs(qs[3]-0.981118) < 1e-5 and abs(qs[5]-0.933452) < 1e-5)

# ---------------------------------------------------------------- S4
print("S4  other modes and lines, configuration (0,0,0), step 1")
by_mode = [(n, c.marks(c.cap(n, m), 1.0)[2].imag) for n in c.modes]
by_line = [(x, c.marks(c.cap(0, x), 1.0)[2].imag) for x in c.unique if x > 0]
print("      by mode:", [(n, round(float(q), 5)) for n, q in by_mode])
print("      by line:", [(round(float(x), 3), round(float(q), 5)) for x, q in by_line])
check("the four modes give rates between 0.988 and 0.992", min(q for _, q in by_mode) > 0.988 and max(q for _, q in by_mode) < 0.992)
check("the heavy lines give small rates (0.093 and 0.0097 for the values 4.174 and 5.568): the choice of the line matters",
      by_line[-1][1] < 0.02 and by_line[-2][1] < 0.1 and by_line[0][1] > 0.98)

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
