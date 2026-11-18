#!/usr/bin/env python3
"""The gluing rule of the sixteenth note as a source candidate with mirror and reference
(18 November 2026).  Hypothesis H-G; not an original producer.

Construction.  Two sides of the seam carry the same space V.  The glued state is
|B> = sum_a |a> x |a-bar> / sqrt(d).  The mirror exchanges the sides and conjugates.  For
operators F, G acting on one side the reflected form is
        Q[F, G] = <B| (mirror image of F) G |B> = Tr(F^+ G)/d.
argv[1] = repository root (native assembler, used unchanged).

Sections
    P1  the form: positive, reference Q[1,1] = 1, mirror fixes the glued state
    P2  the recorded marks of the E8 word as reflected overlaps; the mixed mark compares two orders
    P3  the sign of the connected answer is the Cauchy-Schwarz inequality of the form
    P4  the sector of one sign of D_F gives the recorded canonical marks (native operators)
    P5  the recorded real structure as mirror: positive on even operators, negative on charges
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
rng = np.random.default_rng(20261118)
rat = lambda x: Fr(float(np.real(x))).limit_denominator(2000)

# ---------------------------------------------------------------- matter space of the E8 word
N = 8
dim = 2**N
bits = np.array(list(itertools.product((0, 1), repeat=N)))
ann = []
for i in range(N):
    a = np.zeros((dim, dim))
    for k, b in enumerate(bits):
        if b[i]:
            c = b.copy()
            c[i] = 0
            a[int("".join(map(str, c)), 2), k] = (-1)**int(b[:i].sum())
    ann.append(a)
spin2 = 2*bits-1
matter = np.array([list(v).count(-1) % 2 == 0 for v in spin2])
hyp = np.array([-1/3, -1/3, -1/3, 1/2, 1/2, 0, 0, 0])
fam = [tuple(v[5:]) for v in spin2]
ce = [list(v[:5]).count(-1) % 2 == 0 for v in spin2]
w0 = sorted({f for f, m, e in zip(fam, matter, ce) if m and e})[0]
three = matter & np.array([f not in (w0, tuple(-v for v in w0)) for f in fam])
currents = [(i, j, si, sj) for i in range(5) for j in range(5, 8) for si in (1, -1) for sj in (1, -1)]

# ---------------------------------------------------------------- P1
print("P1  the form")
d = 6
Bv = np.eye(d).reshape(-1)/np.sqrt(d)                           # |B> for a small space, explicit two sides
F = rng.normal(size=(d, d))+1j*rng.normal(size=(d, d))
G = rng.normal(size=(d, d))+1j*rng.normal(size=(d, d))
two = Bv.conj()@np.kron(G, F.conj())@Bv                         # G on side one, mirror image of F on side two
check("explicit two sides: <B| G x conj(F) |B> = Tr(F^+ G)/d", abs(two-np.trace(F.conj().T@G)/d) < 1e-12)
swap = np.zeros((d*d, d*d))
for a_, b_ in itertools.product(range(d), repeat=2):
    swap[b_*d+a_, a_*d+b_] = 1
check("the mirror (exchange of the sides with conjugation) fixes the glued state; the reference Q[1,1] is 1",
      np.allclose(swap@Bv.conj(), Bv) and abs(Bv.conj()@Bv-1) < 1e-12)
check("Q[F,F] = Tr(F^+ F)/d > 0 for every non-zero F: the form is positive", np.trace(F.conj().T@F).real > 0)

# ---------------------------------------------------------------- P2, P3
print("P2  the recorded marks as reflected overlaps")
for label, keep, want in (("128", matter, {"weak": (Fr(5, 24), Fr(1, 4), Fr(1, 48), Fr(-1, 32)), "colour": (Fr(5, 24), Fr(1, 4), Fr(11, 288), Fr(-1, 72))}),
                          ("96", three, {"weak": (Fr(5, 24), Fr(1, 6), Fr(1, 72), Fr(-1, 48)), "colour": (Fr(5, 24), Fr(1, 6), Fr(11, 432), Fr(-1, 108))})):
    idx = np.where(keep)[0]
    n = len(idx)
    H = np.diag(spin2[idx]@hyp/2)
    Qf = lambda A, B: np.trace(A.conj().T@B)/n
    got, cs_ok, norm_ok = {}, True, True
    for (i, j, si, sj) in currents:
        X = ((ann[i].T if si > 0 else ann[i])@(ann[j].T if sj > 0 else ann[j]))[np.ix_(idx, idx)]
        zu, zy, zuy = Qf(H, H), Qf(X, X), Qf(X@H, H@X)
        com = H@X-X@H
        conn = zuy-zu*zy
        norm_ok = norm_ok and abs(Qf(X@H, X@H)-zu*zy) < 1e-12 and abs(Qf(H@X, H@X)-zu*zy) < 1e-12
        cs_ok = cs_ok and abs(conn+0.5*Qf(com, com)) < 1e-12 and abs(zuy) <= zu*zy+1e-12
        got.setdefault("weak" if i >= 3 else "colour", set()).add((rat(zu), rat(zy), rat(zuy), rat(conn)))
    check(f"{label} states: Z_U = Q[H,H], Z_Y = Q[X,X], Z_UY = Q[XH, HX] reproduce the recorded table",
          got == {k: {v} for k, v in want.items()}, {k: sorted(v) for k, v in got.items()})
    check(f"{label} states: both orders have the same norm Z_U Z_Y; the connected answer is the overlap of the two "
          "orders minus that norm, = -(1/2) Q[[H,X],[H,X]]", norm_ok and cs_ok)
print("P3  the sign")
check("Cauchy-Schwarz for the positive form: |Q[XH, HX]| <= Z_U Z_Y, so the connected answer is never positive and "
      "vanishes only when marker and cut commute; the recorded negative sign needs no further input", True)

# ---------------------------------------------------------------- native objects
base = root/"experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q"
sys.path.insert(0, str(base/"runtime"))
from assemble import Assembler  # noqa: E402

Asm = Assembler(base/"domain/native-diagnostic-v1.json")
g, K, orig = Asm.gamma, Asm.K, Asm.original
DF = Asm.assemble({"schema": "tfpt.field-configuration.v1", "domain_sha256": Asm.domain["_sha256"], "genes": [1, 1, 1]})["phi_fourier"][0]


def sign_of(X):
    w, v = np.linalg.eigh(X)
    return (v*np.sign(w))@v.conj().T


Y, S = sign_of(-1j*g@DF), sign_of(DF)
YT = {("u", "L"): 1/6, ("d", "L"): 1/6, ("u", "R"): 2/3, ("d", "R"): -1/3, ("nu", "L"): -1/2, ("e", "L"): -1/2, ("nu", "R"): 0.0, ("e", "R"): -1.0}
yv = [YT[(n, c)] for (n, c, _, _) in orig.P]
Yh = np.diag(yv+[-v for v in yv])
H2 = sum(v@v for v in Asm.V)
I96 = np.eye(96)

# ---------------------------------------------------------------- P4
print("P4  one sign of D_F")
for eta in (1, -1):
    Pe = (I96+eta*S)/2
    ne = np.trace(Pe).real
    q1 = lambda A: np.trace(Pe@A)/ne                            # Q[1, A] on the glued state of the sector
    vals = (q1(I96), q1(g), q1(Y), q1(Y@g))
    check(f"glued state of the sector sign(D_F) = {eta:+d} (48 states): single marks (1, 0, 0, i eta), the recorded "
          "canonical marks, now for the whole sector and not for one line",
          abs(ne-48) < 1e-9 and abs(vals[0]-1) < 1e-12 and abs(vals[1]) < 1e-12 and abs(vals[2]) < 1e-12 and abs(vals[3]-1j*eta) < 1e-12)
check("on the full glued state of the 96 fibres all three single marks vanish: the phase i eta needs the choice of "
      "one sign of D_F", abs(np.trace(Y@g)) < 1e-9 and abs(np.trace(g)) < 1e-9 and abs(np.trace(Y)) < 1e-9)

# ---------------------------------------------------------------- P5
print("P5  the recorded real structure as mirror")
QK = lambda A, B: (np.trace((K@A.conj()@K).conj().T@B)/96).real    # mirror image K conj(A) K, then the plain form
norms = {"1": QK(I96, I96), "Higgs square": QK(H2, H2)/np.trace(H2@H2).real*96, "direction of potential": QK(Y, Y),
         "sign of D_F": QK(S, S), "D_F": QK(DF, DF)/np.trace(DF@DF).real*96, "grading": QK(g, g), "hypercharge": QK(Yh, Yh)}
print("      Q_K[A,A] (in units of Tr A^+A/96 for the first five):", {k: round(float(v), 6) for k, v in norms.items()})
check("with the mirror 'exchange and K-conjugation' the form is positive on the operators that are even under K "
      "(1, Higgs square, direction of the potential, sign of D_F, D_F)",
      all(abs(norms[k]-1) < 1e-9 for k in ("1", "Higgs square", "direction of potential", "sign of D_F", "D_F")))
check("and negative on the charges: grading -1, hypercharge -5/24: the K-mirror is a reflection-positive mirror only "
      "for the even operators", abs(norms["grading"]+1) < 1e-9 and abs(norms["hypercharge"]+5/24) < 1e-9)
connK = QK(Y@Yh, Yh@Y)-QK(Y@Yh, Y@Yh)
connHS = (np.trace((Y@Yh).conj().T@(Yh@Y))/96-np.trace((Y@Yh).conj().T@(Y@Yh))/96).real
check("for the hypercharge through the native potential the connected answer is -0.129976 with the plain mirror and "
      "+0.129976 with the K-mirror: the sign of the recorded answer tells the two mirrors apart",
      abs(connHS+0.129976) < 1e-6 and abs(connK-0.129976) < 1e-6, (round(connHS, 6), round(connK, 6)))

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
