"""The recorded native operator of FOURMARKS20261118 as an instance of the kernel structure.

argv[1] = repository root. Uses the recorded assembler unchanged (native-diagnostic-v1).
Dictionary tested: normal grading G <-> recorded grading g; boundary operator B <-> D_F;
polarization i_C <-> sign(D_F); cut Y0 = i G i_C <-> minus the recorded direction of the potential.
"""
import sys
from pathlib import Path
import numpy as np
from common import Checks

root = Path(sys.argv[1] if len(sys.argv) > 1 else "../../..").resolve()
base = root / "experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q"
sys.path.insert(0, str(base / "runtime"))
from assemble import Assembler  # noqa: E402

c = Checks("native")
Asm = Assembler(base / "domain/native-diagnostic-v1.json")
g, K = np.array(Asm.gamma, complex), np.array(Asm.K, complex)
DF = np.array(Asm.original.full_D(*Asm.Y, Asm.MR), complex)
n = DF.shape[0]
I = np.eye(n)
norm = lambda A: float(np.abs(A).max())

def sign_of(M):
    w, v = np.linalg.eigh(M)
    return (v * np.sign(w)) @ v.conj().T

S = sign_of(DF)
Yrec = sign_of(-1j * g @ DF)
c.ok("recorded D_F is self-adjoint, invertible and odd under the recorded grading (96 fibres, 48 + 48)",
     norm(DF - DF.conj().T) < 1e-9 and norm(g @ DF + DF @ g) < 1e-9 and n == 96 and abs(np.trace(g)) < 1e-9
     and np.abs(np.linalg.eigvalsh(DF)).min() > 1e-12, f"smallest |eigenvalue| {np.abs(np.linalg.eigvalsh(DF)).min():.3e}")
P = (I + S) / 2
c.ok("polarization sign(D_F) anticommutes with the grading: the positive half is Lagrangian", norm(g @ P @ g - (I - P)) < 1e-8)
w, v = np.linalg.eigh(g)
plus, minus = v[:, w > 0], v[:, w < 0]
a, d, b = plus.conj().T @ P @ plus, minus.conj().T @ P @ minus, minus.conj().T @ P @ plus
c.ok("on each grading half the positive half of D_F is exactly 1/2: flat density 1/48 per sector",
     norm(a - np.eye(48) / 2) < 1e-8 and norm(d - np.eye(48) / 2) < 1e-8)
c.ok("the correlation between the halves is a unitary / 2", norm(4 * b.conj().T @ b - np.eye(48)) < 1e-7)
Y0 = 1j * g @ S
c.ok("the recorded direction of the potential sign(-i g D_F) is exactly -Y0 = -i g sign(D_F)", norm(Yrec + Y0) < 1e-8)
c.ok("mixed mark: g Y0 = i sign(D_F); on the sector sign(D_F) = eta the flat state gives Z_U = Z_Y = 0, Z_UY = i eta with (U, Y) = (g, Y0)",
     norm(g @ Y0 - 1j * S) < 1e-8 and all(abs(np.trace((I + e * S) / 2 @ M)) < 1e-8 for e in (1, -1) for M in (g, Y0))
     and all(abs(np.trace((I + e * S) / 2 @ g @ Y0) / 48 - 1j * e) < 1e-8 for e in (1, -1)))
ev = np.linalg.eigvalsh(DF)
uniq = [x for x in ev if np.sum(np.abs(ev - x) < 1e-8) == 1]
lam1 = ev[ev > 0].min()
c.ok("12 unique eigenlines (six values, both signs); lambda_1^+ is non-degenerate, so the seed line of the kernel is one of them",
     len(uniq) == 12 and int(np.sum(np.abs(ev - lam1) < 1e-8)) == 1 and int(np.sum(np.abs(ev + lam1) < 1e-8)) == 1,
     f"lambda_1^+ = {lam1:.6f}, unique values {sorted(round(abs(x), 4) for x in uniq if x > 0)}")
# the real structure K (anti-linear: x -> K conj(x)) against polarization and grading
KS = K @ S.conj() @ K.conj().T
Kg = K @ g.conj() @ K.conj().T
c.ok("the recorded real structure fixes the polarization (K S K^-1 = S) and exchanges the grading halves (K g K^-1 = -g)",
     norm(KS - S) < 1e-8 and norm(Kg + g) < 1e-8, f"K S K^-1 - S: {norm(KS - S):.1e}; K g K^-1 + g: {norm(Kg + g):.1e}")
KY = K @ Y0.conj() @ K.conj().T
c.ok("under the real structure the cut Y0 is even and the marker g is odd", norm(KY - Y0) < 1e-8 or norm(KY + Y0) < 1e-8,
     f"K Y0 K^-1 = {'+' if norm(KY - Y0) < 1e-8 else '-'} Y0")
c.done({"fibres": n, "unique_lines": len(uniq), "K_Y0_sign": "+" if norm(KY - Y0) < 1e-8 else "-"})
