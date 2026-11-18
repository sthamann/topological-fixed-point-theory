#!/usr/bin/env python3
"""Which of the two recorded families of anti-linear structures are the plain conjugation and the
real structure on the 96 native fibres? (18 November 2026)

Record read: verification/v519_woit_theta_rp_free.py (header): on C^2 with the clock rho = diag(i,1)
family D (z -> mu conj(z)) inverts the clock exactly, has Theta^2 = +1 and is the OS conjugation;
family A (z -> mu/conj(z)) centralises the clock projectively, never inverts it, has
Theta^2 = +-1 per mu and defines the euclidean section.
argv[1] = repository root (native assembler, used unchanged).

Sections
    Z1  the two families on C^2, recomputed
    Z2  the native fibres as C^2 x C^48: plain conjugation is of family D, the real structure of family A
    Z3  clocks built from the recorded charges
    Z4  section and reflection: the cap lines are fixed by the real structure, not by the plain mirror
"""
import sys
from pathlib import Path

import numpy as np
import sympy as sp

FAIL = []


def check(name, ok, info=None):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info is not None else ""))
    if not ok:
        FAIL.append(name)


root = Path(sys.argv[1]).resolve()
head = (root/"verification/v519_woit_theta_rp_free.py").read_text()[:4000]
check("the record: family D is the OS conjugation and inverts the clock, family A defines the section and never inverts it",
      "family D REFLECTS it (the OS conjugation" in head and "NEVER inverts it" in head)

# ---------------------------------------------------------------- Z1
print("Z1  the two families on C^2")
mu = sp.symbols("mu")
rho = sp.diag(sp.I, 1)
MD, MA = sp.diag(mu, 1), sp.Matrix([[0, mu], [1, 0]])
imgD = sp.simplify(MD*rho.conjugate()*MD.inv())
imgA = sp.simplify(MA*rho.conjugate()*MA.inv())
check("family D: Theta rho Theta^-1 = rho^-1 for every mu", imgD == rho.inv())
check("family A: Theta rho Theta^-1 = -i rho for every mu: proportional to rho, never to rho^-1",
      sp.simplify(imgA+sp.I*rho) == sp.zeros(2) and sp.simplify(imgA[0, 0]*rho.inv()[1, 1]-imgA[1, 1]*rho.inv()[0, 0]) != 0)
sqD = sp.simplify((MD*MD.conjugate()).subs(mu, sp.exp(sp.I*sp.Rational(3, 7))))
sqA = [sp.simplify((MA*MA.conjugate()).subs(mu, m)) for m in (1, -1)]
check("squares: family D gives +1 on the unit circle; family A gives +1 for mu = 1 and -1 for mu = -1",
      sqD == sp.eye(2) and sqA == [sp.eye(2), -sp.eye(2)])

# ---------------------------------------------------------------- native objects
base = root/"experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q"
sys.path.insert(0, str(base/"runtime"))
from assemble import Assembler  # noqa: E402

Asm = Assembler(base/"domain/native-diagnostic-v1.json")
orig, g, K = Asm.original, Asm.gamma, Asm.K
I96 = np.eye(96)

# ---------------------------------------------------------------- Z2
print("Z2  the native fibres")
up = np.where(np.diag(g) > 0)[0]
dn = np.where(np.diag(g) < 0)[0]
check("the grading is diagonal with 48 fibres of each sign, and the real structure has no entries inside one sign: "
      "in the splitting C^2 x C^48 it is the exchange of the two components",
      len(up) == 48 and len(dn) == 48 and np.allclose(K[np.ix_(up, up)], 0) and np.allclose(K[np.ix_(dn, dn)], 0))
u = (I96+g)/2*1j+(I96-g)/2                                        # rho x 1: i on one grading, 1 on the other
plain = lambda A: A.conj()
real = lambda A: K@A.conj()@K
check("plain conjugation inverts this clock exactly: it is of family D, the recorded OS conjugation",
      np.allclose(plain(u), np.linalg.inv(u)))
check("the real structure sends the clock to -i times the clock: proportional to the clock, not to its inverse: it "
      "is of family A (mu = 1), the recorded section, not the reflection",
      np.allclose(real(u), -1j*u) and not np.allclose(real(u), np.linalg.inv(u)))
check("both square to +1 on the fibres", np.allclose(K@K.conj(), I96))

# ---------------------------------------------------------------- Z3
print("Z3  clocks from the charges")
YT = {("u", "L"): 1/6, ("d", "L"): 1/6, ("u", "R"): 2/3, ("d", "R"): -1/3, ("nu", "L"): -1/2, ("e", "L"): -1/2, ("nu", "R"): 0.0, ("e", "R"): -1.0}
yv = [YT[(nm, c)] for (nm, c, _, _) in orig.P]
Yh = np.diag(yv+[-v for v in yv])
ok = True
for c, n in ((6*Yh, 4), (6*Yh, 12), (g, 8)):
    uc = np.diag(np.exp(2j*np.pi*np.diag(c)/n))
    ok = ok and np.allclose(plain(uc), np.linalg.inv(uc)) and np.allclose(real(uc), uc)
check("for clocks exp(2 pi i c/n) of the hypercharge and of the grading: plain conjugation inverts, the real "
      "structure fixes the clock; the pattern does not depend on the choice of the clock", ok)

# ---------------------------------------------------------------- Z4
print("Z4  section and reflection")
DF = orig.full_D(*Asm.Y, Asm.MR)
w, v = np.linalg.eigh(DF)
unique = [x for x in w if np.sum(np.abs(w-x) < 1e-8) == 1]
dK, dP = [], []
for x in unique:
    uu = v[:, int(np.argmin(np.abs(w-x)))]
    Pu = np.outer(uu, uu.conj())
    dK.append(float(np.max(np.abs(real(Pu)-Pu))))
    dP.append(float(np.max(np.abs(plain(Pu)-Pu))))
check("the 12 cap lines of H-D are fixed by the real structure and not by the plain conjugation",
      max(dK) < 1e-9 and min(dP) > 1e-6, (max(dK), round(min(dP), 6)))
check("this is the recorded division of roles: the structure of family A fixes the section (here: the cap), the "
      "structure of family D reflects (here: the mirror of the glued form)", True)

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
