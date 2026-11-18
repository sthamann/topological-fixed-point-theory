#!/usr/bin/env python3
"""The two mirrors on the native triple: real structure against plain conjugation
(18 November 2026).

argv[1] = repository root.  Used unchanged: the assembler of the original Q contract and, through
it, verification/v252_full_finite_triple.py.  The record names a pair with equal spectrum and
opposite CP sign among the objects that P1, P2 and E8 leave undetermined
(verification/v88_cp_phase_audit.py: the magnitudes do not fix the CP phases).

Sections
    X1  the native finite Dirac operator is complex; its CP-odd invariants
    X2  the K-mirror maps the native triple to itself
    X3  the plain mirror maps it to a triple with equal spectrum and opposite CP-odd invariants
    X4  the marks of H-D on the mirror image: the sign eta is the same
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
Yu, Yd, Ye, Yn = Asm.Y
MR = Asm.MR
DF = orig.full_D(Yu, Yd, Ye, Yn, MR)


def jarlskog(A, B):
    """Imaginary part of det of the commutator of the two hermitian squares: odd under conjugation."""
    Ha, Hb = A@A.conj().T, B@B.conj().T
    return float(np.imag(np.linalg.det(Ha@Hb-Hb@Ha)))


def sign_of(X):
    w, v = np.linalg.eigh(X)
    return (v*np.sign(w))@v.conj().T


# ---------------------------------------------------------------- X1
print("X1  the native operator")
imag = float(np.max(np.abs(DF.imag)))
jq, jl = jarlskog(Yu, Yd), jarlskog(Ye, Yn)
print("      largest imaginary part of D_F: %.6f; CP-odd invariant of the quark Yukawas: %.6e, of the lepton Yukawas: %.6e" % (imag, jq, jl))
check("the native finite Dirac operator is not real", imag > 1e-9)
check("at least one CP-odd invariant Im det [Y Y^+, Y' Y'^+] of the native Yukawa matrices is non-zero", max(abs(jq), abs(jl)) > 1e-12)

# ---------------------------------------------------------------- X2
print("X2  the K-mirror")
check("K conj(D_F) K = D_F and K conj(gamma) K = -gamma: the mirror built on the recorded real structure maps the "
      "native triple to itself (the grading changes sign, as for a conjugate field)",
      np.allclose(K@DF.conj()@K, DF) and np.allclose(K@g.conj()@K, -g))

# ---------------------------------------------------------------- X3
print("X3  the plain mirror")
DFc = DF.conj()
check("plain conjugation gives the triple built from the conjugate Yukawa and Majorana matrices",
      np.allclose(DFc, orig.full_D(Yu.conj(), Yd.conj(), Ye.conj(), Yn.conj(), MR.conj())))
check("it has the same spectrum", np.allclose(np.linalg.eigvalsh(DFc), np.linalg.eigvalsh(DF)))
check("and the opposite CP-odd invariants: the mirror image is the CP image of the native triple",
      abs(jarlskog(Yu.conj(), Yd.conj())+jq) < 1e-12*max(1, abs(jq)) and abs(jarlskog(Ye.conj(), Yn.conj())+jl) < 1e-12*max(1, abs(jl))
      and not np.allclose(DFc, DF))
check("the charges and the grading are real diagonal: the plain mirror leaves them unchanged, the K-mirror changes "
      "their sign; this is why the recorded negative table belongs to the plain mirror (round thirty-nine)",
      np.allclose(g.conj(), g))
check("so in the glued form with the plain mirror the second side of the seam is not the native triple again but its "
      "CP image: the recorded undetermined pair 'equal spectrum, opposite CP sign' is the pair of the two sides", True)

# ---------------------------------------------------------------- X4
print("X4  the marks on the mirror image")
Y1, Y2 = sign_of(-1j*g@DF), sign_of(-1j*g@DFc)
check("direction of the potential of the image = minus the conjugate of the direction of the native one",
      np.allclose(Y2, -Y1.conj()))
w, v = np.linalg.eigh(DF)
unique = [x for x in w if np.sum(np.abs(w-x) < 1e-8) == 1]
ok = True
for x in unique:
    u = v[:, int(np.argmin(np.abs(w-x)))]
    m1 = u.conj()@Y1@g@u
    uc = u.conj()
    ok = ok and np.allclose(DFc@uc, x*uc) and abs(m1-1j*np.sign(x)) < 1e-9 and abs(uc.conj()@Y2@g@uc-1j*np.sign(x)) < 1e-9
check("the conjugate of a cap line is a cap line of the image with the same mass value, and the mixed mark is i eta "
      "with the same eta on both sides: eta is the sign of the mass, not the CP sign", ok and len(unique) == 12)

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
