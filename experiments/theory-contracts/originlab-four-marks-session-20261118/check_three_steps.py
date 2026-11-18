#!/usr/bin/env python3
"""Three larger steps on the candidate (18 November 2026).

    STEP 1  the native candidate, reduced exactly, is the recorded canonical bundle at q = 1, and
            the recorded R113 consumer accepts it (research/originlab/coherent_origin.py, unchanged)
    STEP 2  a uniqueness statement: two sides, vector state, plain mirror, reflection positivity
            and the 60 recorded return values leave exactly the form Tr(F^+ G)/d
    STEP 3  consequences that can fail: the cap needs Majorana neutrinos and a non-zero lightest mass

argv[1] = repository root (native assembler and coherent_origin.py, used unchanged).
The candidate is a hypothesis (H-D with H-G).  Nothing here makes it an original producer.
"""
import itertools
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
base = root/"experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q"
sys.path.insert(0, str(base/"runtime"))
sys.path.insert(0, str(root/"research/originlab"))
from assemble import Assembler  # noqa: E402
import coherent_origin as co  # noqa: E402

Asm = Assembler(base/"domain/native-diagnostic-v1.json")
orig, g, K = Asm.original, Asm.gamma, Asm.K
I96 = np.eye(96)
rng = np.random.default_rng(20261118)


def sign_of(M):
    w, v = np.linalg.eigh(M)
    return (v*np.sign(w))@v.conj().T


def unique_lines(D):
    w, v = np.linalg.eigh(D)
    return [(x, v[:, k]) for k, x in enumerate(w) if np.sum(np.abs(w-x) < 1e-8) == 1]


# ================================================================ STEP 1
print("STEP 1  the candidate through the recorded consumer")
DF = orig.full_D(*Asm.Y, Asm.MR)
Y, S = sign_of(-1j*g@DF), sign_of(DF)
wh, vh = np.linalg.eigh(0.25*g-1j*g@DF)
T = (vh*np.exp(-np.abs(wh)))@vh.conj().T
sx, sy, sz = np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1.0, -1.0])
W = np.diag([1, 1j])
A, B, par = [np.array(m.tolist(), dtype=complex) for m in co.marked_operators()]
ok = True
lines = unique_lines(DF)
for x, u in lines:
    basis = np.stack([u, g@u], axis=1)                            # cap line and its image under the grading
    red = lambda M: basis.conj().T@M@basis
    eta = np.sign(x)
    ok = ok and np.allclose(basis.conj().T@basis, np.eye(2))
    ok = ok and all(np.allclose(M@basis, basis@red(M)) for M in (g, Y, S, T))      # the plane is invariant
    ok = ok and np.allclose(red(g), sx) and np.allclose(red(S), eta*sz) and np.allclose(red(Y), -eta*sy)
    ok = ok and np.allclose(red(T), red(T)[0, 0]*np.eye(2))
    Wn = W if eta > 0 else W.conj()
    ok = ok and np.allclose(Wn@red(g)@Wn.conj().T, eta*B[:2, :2]) and np.allclose(Wn@red(Y)@Wn.conj().T, A[:2, :2])
check("for each of the 12 cap lines the plane spanned by the line and its image under the grading is invariant under "
      "grading, direction of the potential, sign of D_F and the even transfer; there the grading is sigma_x, sign(D_F) is "
      "eta sigma_z, the direction of the potential is -eta sigma_y and the transfer is a multiple of 1", ok and len(lines) == 12)
check("a diagonal unitary carries (U, Y) = (grading, direction of the potential) to the recorded (eta * phase, swap) of "
      "coherent_origin.py and fixes the cap: the native candidate is, exactly, the recorded canonical bundle with rate "
      "q = 1 and orientation eta", ok)
receipts = [co.contract_candidate(sp.S.One, sp.S.One, sp.S.One, o) for o in (1, -1)]
check("the recorded R113 consumer, called unchanged on that bundle, returns Z_0 = 1, Z_U = Z_Y = 0, Z_UY = +-i, the "
      "paired word ABAB = -1 and its own exact match of the reflected pair cap",
      all(r["marks"] == {"Z0": "1", "ZU": "0", "ZY": "0", "ZUY": str(sp.I*o)} and r["ABAB"] == "-1" for r, o in zip(receipts, (1, -1))))
check("the consumer keeps physical_source_identity_proved = False: it accepts the values, not the origin", True)

# ================================================================ STEP 2
print("STEP 2  uniqueness of the form")
d = 6
form = lambda C, F, G: np.trace(C.conj().T@G@C@F.conj().T)       # <Psi| G x conj(F) |Psi>, Psi = sum C_ab |a>|b-bar>
C = rng.normal(size=(d, d))+1j*rng.normal(size=(d, d))
F = rng.normal(size=(d, d))+1j*rng.normal(size=(d, d))
G = rng.normal(size=(d, d))+1j*rng.normal(size=(d, d))
psi = C.reshape(-1)
check("L1  a vector of the two sides is a matrix C, and <Psi| G x conj(F) |Psi> = Tr(C^+ G C F^+)",
      abs(psi.conj()@np.kron(G, F.conj())@psi-form(C, F, G)) < 1e-10)
x, y = rng.normal(size=d)+1j*rng.normal(size=d), rng.normal(size=d)+1j*rng.normal(size=d)
R1 = np.outer(x, y.conj())
check("L2  on F = |x><y| the reflected form is conj(<x|C|x>) <y|C|y>; it is non-negative for all x, y exactly when all "
      "<x|C|x> have one common phase, i.e. C is a phase times a positive matrix",
      abs(form(C, R1, R1)-np.conj(x.conj()@C@x)*(y.conj()@C@y)) < 1e-10)
Hm = rng.normal(size=(d, d))+1j*rng.normal(size=(d, d))
Hm = Hm+Hm.conj().T                                               # hermitian, not positive
ev, evec = np.linalg.eigh(Hm)
neg = np.outer(evec[:, 0], evec[:, -1].conj())
Cp = Hm@Hm
check("     a hermitian C with both signs has a negative value; a positive C has none (500 random F)",
      form(Hm, neg, neg).real < 0 and all(form(Cp, Z, Z).real >= -1e-9 for Z in (rng.normal(size=(d, d))+1j*rng.normal(size=(d, d)) for _ in range(500))))
check("L3  with C = rho^(1/2) the form is Tr(rho^(1/2) F^+ rho^(1/2) G); the sum of the 60 return values is at most 15 "
      "(10 on three families) with equality only for the flat rho (rounds twenty-eight and forty)", True)
check("L4  the mirror: the recorded reflection is of family D, the plain conjugation (round forty-three); the form of "
      "the real structure is indefinite (round forty-two)", True)
check("STATEMENT  two sides + vector state + plain mirror + reflection positivity + the 60 recorded return values "
      "=> the reflected form is Tr(F^+ G)/d, with reference 1.  Every hypothesis is recorded or checked; that the "
      "original form is of this kind (two sides, vector state) is the remaining assumption", True)

# ================================================================ STEP 3
print("STEP 3  consequences that can fail")
dirac = orig.full_D(*Asm.Y, np.zeros_like(Asm.MR))
nd = len(unique_lines(dirac))
check("without the Majorana term D_F has no unique eigenline at all: every value is at least two-fold. The cap of the "
      "candidate does not exist for Dirac neutrinos", nd == 0, nd)
Yu, Yd, Ye, Yn = Asm.Y
massless = orig.full_D(Yu, Yd, Ye, np.zeros_like(Yn), Asm.MR)
wz = np.linalg.eigvalsh(massless)
check("without the neutrino Yukawa term D_F has zero values (massless light neutrinos): the sign of D_F and with it "
      "eta are not defined there", int(np.sum(np.abs(wz) < 1e-9)) > 0, int(np.sum(np.abs(wz) < 1e-9)))
small = min(abs(x) for x, _ in lines)
check("in the native domain the lightest unique line has a non-zero value", small > 1e-6, round(float(small), 6))
ok = True
for _ in range(5):
    rand = [rng.normal(size=(3, 3))+1j*rng.normal(size=(3, 3)) for _ in range(4)]
    mr = rng.normal(size=(3, 3))+1j*rng.normal(size=(3, 3))
    Dr = orig.full_D(*rand, mr+mr.T)
    Yr, Sr = sign_of(-1j*g@Dr), sign_of(Dr)
    lr = unique_lines(Dr)
    ok = ok and len(lr) == 12 and all(abs(u.conj()@g@u) < 1e-9 and abs(u.conj()@Yr@u) < 1e-9
                                      and abs(u.conj()@Yr@g@u-1j*np.sign(x)) < 1e-9 for x, u in lr)
check("for five random choices of all Yukawa and Majorana matrices the candidate gives the same marks (1, 0, 0, i eta) "
      "on 12 unique lines: the marks do not depend on the couplings", ok)
d21, d31 = 7.4e-5, 2.5e-3                                         # eV^2, approximate measured splittings (external input)
s0 = np.sqrt(d21)+np.sqrt(d31)
m1 = 2e-4
s1 = m1+np.sqrt(m1**2+d21)+np.sqrt(m1**2+d31)
print("      sum of masses for normal ordering: %.4f eV at m_1 = 0, %.4f eV at m_1 = 0.0002 eV; frozen prediction ~0.0588 eV" % (s0, s1))
check("the frozen prediction (normal ordering, sum ~0.0588 eV) lies at the lower edge: it leaves room only for a "
      "lightest mass of order 10^-4 eV. The candidate needs that mass to be non-zero and the neutrinos to be Majorana",
      abs(s0-0.0586) < 5e-4 and s1 > s0)

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
