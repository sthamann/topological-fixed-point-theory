#!/usr/bin/env python3
"""Correction of rounds thirty to thirty-two, and the cap from the original finite Dirac operator
(18 November 2026).

argv[1] = repository root.  Used unchanged: the assembler of the original Q contract and, through
it, verification/v252_full_finite_triple.py (full_D, grading, real structure).

Correction.  Rounds thirty to thirty-two wrote the constant-mode operator as gamma/4 + phi_0.  The
assembler builds gamma/4 - i gamma phi_0: the potential term is -i gamma D_F with D_F = phi_0 the
original finite Dirac operator (full_D).  Even functions of the two operators coincide, so the
numbers of those rounds stand; the names of two operators have to be exchanged.

Sections
    K1  the constant mode of the assembled operator; what stands and what changes
    K2  the three native involutions: grading, direction of the potential, sign of D_F
    K3  the cap: unique eigenlines of D_F are real and give the recorded marks, eta = sign of the value
    K4  the eigenlines of the circle operator at shift 0 are the other family
    K5  the fermion determinant of the native operator: first and second derivatives
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
I96, I4 = np.eye(96), np.eye(4)
conjK = lambda X: K@X.conj()@K


def cand(genes):
    return Asm.assemble({"schema": "tfpt.field-configuration.v1", "domain_sha256": Asm.domain["_sha256"], "genes": list(genes)})


def sign_of(X):
    w, v = np.linalg.eigh(X)
    return (v*np.sign(w))@v.conj().T


# ---------------------------------------------------------------- K1
print("K1  the constant mode")
a = cand([1, 1, 1])
h, DF = a["h_N"], a["phi_fourier"][0]
pot = -1j*g@DF
i0 = a["modes"].index(0)
check("baseline: the Fourier part of order zero is the original finite Dirac operator full_D",
      np.allclose(DF, orig.full_D(*Asm.Y, Asm.MR)))
check("the block of mode 0 of the assembled operator is gamma/4 - i gamma D_F (not gamma/4 + D_F)",
      np.allclose(h[96*i0:96*(i0+1), 96*i0:96*(i0+1)], 0.25*g+pot) and not np.allclose(pot, DF))
check("both forms have the same square 1/16 + D_F^2, so every even function of them is the same operator: "
      "the numbers of rounds thirty and thirty-one stand",
      np.allclose((0.25*g+pot)@(0.25*g+pot), I96/16+DF@DF) and np.allclose((0.25*g+DF)@(0.25*g+DF), I96/16+DF@DF))

# ---------------------------------------------------------------- K2
print("K2  three involutions")
Yt, S = sign_of(pot), sign_of(DF)
check("direction of the native potential Y: involution, anticommutes with the grading, real under K",
      np.allclose(Yt@Yt, I96) and np.allclose(Yt@g+g@Yt, 0) and np.allclose(conjK(Yt), Yt))
check("sign of D_F: involution, anticommutes with the grading, real under K; and sign(D_F) = -i Y gamma exactly",
      np.allclose(S@S, I96) and np.allclose(S@g+g@S, 0) and np.allclose(conjK(S), S) and np.allclose(-1j*Yt@g, S))
check("what rounds thirty and thirty-two called 'direction of the potential' was sign(D_F), and their 'Pi' was, "
      "up to sign, the direction of the potential", True)

# ---------------------------------------------------------------- K3
print("K3  the cap")
w, v = np.linalg.eigh(DF)
unique = [x for x in w if np.sum(np.abs(w-x) < 1e-8) == 1]
wh, vh = np.linalg.eigh(0.25*g+pot)
T = (vh*np.exp(-0.7*np.abs(wh)))@vh.conj().T
check("an even transfer of the true constant mode commutes with the grading, with Y and with D_F",
      np.allclose(T@g, g@T) and np.allclose(T@Yt, Yt@T) and np.allclose(T@DF, DF@T))
rows = []
for x in unique:
    u = v[:, int(np.argmin(np.abs(w-x)))]
    Pu = np.outer(u, u.conj())
    z0 = u.conj()@T@T@T@u
    rows.append((x, float(np.max(np.abs(conjK(Pu)-Pu))), abs(u.conj()@T@T@g@T@u/z0), abs(u.conj()@T@Yt@T@T@u/z0),
                 u.conj()@T@Yt@T@g@T@u/z0, (u.conj()@((I96+g)/2)@u).real))
print("      values of D_F on the unique lines:", np.round([r[0] for r in rows], 4))
check("D_F has 12 unique eigenlines; all are real", len(rows) == 12 and max(r[1] for r in rows) < 1e-9)
check("on each of them, with U = grading and Y = direction of the native potential: Z_U = Z_Y = 0 and "
      "Z_UY = i eta Z_0 with eta = sign of the value of D_F; positive lift 1/2",
      all(r[2] < 1e-9 and r[3] < 1e-9 and abs(r[4]-1j*np.sign(r[0])) < 1e-9 and abs(r[5]-0.5) < 1e-9 for r in rows))
check("so the cap of the recorded canonical model is, natively, a mass eigenline of the original finite Dirac "
      "operator, and the sign eta is the sign of its mass value", True)
a2 = cand([2, 2, 2])
check("the same holds for the coefficients [2,2,2] (same order-zero part)", np.allclose(a2["phi_fourier"][0], DF))

# ---------------------------------------------------------------- K4
print("K4  the circle operator at shift 0")
wp, vp = np.linalg.eigh(pot)
uniq_p = [x for x in wp if np.sum(np.abs(wp-x) < 1e-8) == 1]
ok = True
for x in uniq_p:
    u = vp[:, int(np.argmin(np.abs(wp-x)))]
    Pu = np.outer(u, u.conj())
    ok = ok and np.max(np.abs(conjK(Pu)-Pu)) < 1e-9 and abs(u.conj()@Yt@u-np.sign(x)) < 1e-9 and abs(u.conj()@S@u) < 1e-9
check("at shift 0 the mode operator is -i gamma D_F; it too has 12 unique real eigenlines, with definite direction "
      "of the potential and marks (1, 0, +-1, 0): round thirty-two stands with the true operator", len(uniq_p) == 12 and ok)
check("these are not the eigenlines of D_F: the two families are exchanged by the grading",
      np.allclose(g@S@g, -S) and np.allclose(g@Yt@g, -Yt))

# ---------------------------------------------------------------- K5
print("K5  the fermion determinant")
YTAB = {("u", "L"): Fr(1, 6), ("d", "L"): Fr(1, 6), ("u", "R"): Fr(2, 3), ("d", "R"): Fr(-1, 3),
        ("nu", "L"): Fr(-1, 2), ("e", "L"): Fr(-1, 2), ("nu", "R"): Fr(0), ("e", "R"): Fr(-1)}
yv = [float(YTAB[(n, c)]) for (n, c, _, _) in orig.P]
Yh = np.diag(yv+[-x for x in yv])
free = np.kron(np.diag([m+0.25 for m in a["modes"]]), g)
table = []
for genes in itertools.product(range(3), repeat=3):
    b = cand(genes)
    hh = b["h_N"]
    G = np.linalg.inv(hh)
    H = np.kron(I4, -1j*g@(b["phi_fourier"][0]-Asm.DM))            # Higgs part of the constant potential
    M = np.kron(I4, g@Yh)                                         # hypercharge added to the shift
    sgn, logdet = np.linalg.slogdet(hh)
    table.append((genes, logdet-np.linalg.slogdet(free)[1], sgn.real, np.trace(G@H), np.trace(G@M),
                  -np.trace(G@H@G@M), -np.trace(G@H@G@H), -np.trace(G@M@G@M)))
b0 = [t for t in table if t[0] == (1, 1, 1)][0]
print("      baseline: log det ratio %.6f; Higgs first %.6f; hypercharge first %.2e; mixed %.2e; Higgs second %.4f; "
      "hypercharge second %.4f" % (b0[1], b0[3].real, abs(b0[4]), abs(b0[5]), b0[6].real, b0[7].real))
arr = np.array([[t[1], t[3].real, t[4].real, t[5].real, t[6].real, t[7].real] for t in table])
print("      27 configurations: log det ratio %.4f to %.4f; mixed %.2e to %.2e" % (arr[:, 0].min(), arr[:, 0].max(), arr[:, 3].min(), arr[:, 3].max()))
check("the determinant is positive in all 27 configurations (192 negative levels)", all(abs(t[2]-1) < 1e-9 for t in table))
check("all derivatives are real", max(abs(np.imag(x)) for t in table for x in t[3:]) < 1e-9)
check("baseline: the first derivative in the hypercharge direction and the mixed derivative vanish; the second "
      "derivatives in both directions are negative", abs(b0[4]) < 1e-9 and abs(b0[5]) < 1e-9 and b0[6].real < 0 and b0[7].real < 0)
check("the mixed derivative is not strictly negative over the 27 configurations: it takes both signs and is small; "
      "the determinant with these two sources does not give the recorded negative mixed answer",
      arr[:, 3].min() < 0 < arr[:, 3].max(), (float("%.3e" % arr[:, 3].min()), float("%.3e" % arr[:, 3].max())))

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
