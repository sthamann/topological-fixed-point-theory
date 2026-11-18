#!/usr/bin/env python3
"""The four marks of hypothesis H-D on the full native operator, all 27 configurations
(18 November 2026).

argv[1] = repository root.  The assembler of the original Q contract is used unchanged.
Candidate (hypothesis, not an original producer):
    cap    = mode 0 of the circle times the unique eigenline of D_F of smallest |value| (both signs)
    U      = grading, Y = direction of the native potential -i gamma D_F
    transfer = exp(-|h_N|) for each of the three equal steps (the recorded order: U early, Y late)
    Z_0 = <cap|T^3|cap>, Z_U = <cap|T^2 U T|cap>, Z_Y = <cap|T Y T^2|cap>, Z_UY = <cap|T Y T U T|cap>

Sections
    J1  baseline: exact marks (1, 0, 0, i eta)
    J2  all 27 configurations: the rate q and the single marks
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


def marks(genes, step=1.0):
    a = Asm.assemble({"schema": "tfpt.field-configuration.v1", "domain_sha256": Asm.domain["_sha256"], "genes": list(genes)})
    h, phi = a["h_N"], a["phi_fourier"]
    DF = phi[0]
    mix = float(np.sqrt(sum(np.linalg.norm(v)**2 for k, v in phi.items() if k != 0)))
    Y = np.kron(I4, sign_of(-1j*g@DF))
    U = np.kron(I4, g)
    w, v = np.linalg.eigh(DF)
    uniq = [x for x in w if np.sum(np.abs(w-x) < 1e-8) == 1]
    small = min(abs(x) for x in uniq)
    wh, vh = np.linalg.eigh(h)
    T = (vh*np.exp(-step*np.abs(wh)))@vh.conj().T
    i0 = a["modes"].index(0)
    out = []
    for x in [x for x in uniq if abs(abs(x)-small) < 1e-9]:
        cap = np.zeros(384, complex)
        cap[96*i0:96*(i0+1)] = v[:, int(np.argmin(np.abs(w-x)))]
        z0 = cap.conj()@T@T@T@cap
        out.append((int(np.sign(x)), z0.real, cap.conj()@T@T@U@T@cap/z0, cap.conj()@T@Y@T@T@cap/z0, cap.conj()@T@Y@T@U@T@cap/z0))
    return mix, len(uniq), small, out


# ---------------------------------------------------------------- J1
print("J1  baseline [1,1,1]")
mix, nuniq, small, out = marks((1, 1, 1))
print("      profile beyond the constant part: %.3f; smallest unique |value| of D_F: %.6f" % (mix, small))
for eta, z0, zu, zy, zuy in out:
    print("      eta %+d: Z_0 = %.6f, Z_U/Z_0 = %.1e, Z_Y/Z_0 = %.1e, Z_UY/Z_0 = %.6f%+.6fi" % (eta, z0, abs(zu), abs(zy), zuy.real, zuy.imag))
check("the baseline has a constant profile; the marks are (1, 0, 0, i eta) to rounding, with eta the sign of the "
      "mass value: rate q = 1", mix < 1e-12 and all(abs(zu) < 1e-10 and abs(zy) < 1e-10 and abs(zuy-1j*eta) < 1e-10 for eta, _, zu, zy, zuy in out))
check("Z_0 = exp(-3 sqrt(1/16 + m^2)) for the line of mass value m", all(abs(z0-np.exp(-3*np.sqrt(1/16+small**2))) < 1e-10 for _, z0, _, _, _ in out))

# ---------------------------------------------------------------- J2
print("J2  all 27 configurations")
table = []
for genes in itertools.product(range(3), repeat=3):
    mix, nuniq, small, out = marks(genes)
    for eta, z0, zu, zy, zuy in out:
        table.append((genes, mix, nuniq, eta, z0, zu, zy, zuy))
q = np.array([t[7].imag*t[3] for t in table])
re = np.array([abs(t[7].real) for t in table])
su = np.array([abs(t[5]) for t in table])
sy = np.array([abs(t[6]) for t in table])
z0s = np.array([t[4] for t in table])
levels = sorted({round(float(x), 5) for x in q})
print("      rate q (imaginary part of Z_UY/Z_0 times eta):", levels)
print("      |Re Z_UY/Z_0| up to %.5f; |Z_U/Z_0| up to %.5f; |Z_Y/Z_0| up to %.5f; Z_0 from %.5f to %.5f"
      % (re.max(), su.max(), sy.max(), z0s.min(), z0s.max()))
check("54 evaluations (27 configurations, both signs); D_F keeps its 12 unique lines in every configuration",
      len(table) == 54 and all(t[2] == 12 for t in table))
check("in every configuration the mixed mark is i eta q with 0.988 < q <= 1; eta is the sign of the mass value",
      q.min() > 0.988 and q.max() <= 1+1e-9, (round(float(q.min()), 5), round(float(q.max()), 5)))
check("q is below one exactly when the profile is not constant: the rate comes from the mixing of the circle modes",
      all((t[1] < 1e-12) == (abs(t[7].imag*t[3]-1) < 1e-9) for t in table))
check("the single marks and the real part of the mixed mark stay small: below 0.03, 0.02 and 0.001",
      su.max() < 0.03 and sy.max() < 0.02 and re.max() < 0.001)
dev = max(abs(table[2*i][7]-np.conj(table[2*i+1][7])) for i in range(27))
check("the two signs of eta give complex conjugate mixed marks in each configuration, to five digits", dev < 1e-5, float("%.1e" % dev))

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
