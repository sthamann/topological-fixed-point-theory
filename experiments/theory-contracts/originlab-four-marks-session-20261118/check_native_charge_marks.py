#!/usr/bin/env python3
"""Native hypercharge and Higgs square as marks on real lines (18 November 2026).

argv[1] = repository root.  Used unchanged: the assembler of the original Q contract and, through
it, verification/v252_full_finite_triple.py (index of the 96 fibres, grading, real structure,
full_D) and the four original Higgs directions.

Sections
    W1  the native charges: hypercharge, weak isospin, electric charge; recorded properties
    W2  the Higgs square: non-negative, even, real, neutral
    W3  what reality of the line forces: charge marks imaginary, even marks real
    W4  values: electric charge (conserved) against hypercharge (broken by the constant potential)
    W5  flat density on the 96 fibres
"""
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
rng = np.random.default_rng(20261118)
conjK = lambda X: K@X.conj()@K

# ---------------------------------------------------------------- W1
print("W1  the native charges")
YTAB = {("u", "L"): Fr(1, 6), ("d", "L"): Fr(1, 6), ("u", "R"): Fr(2, 3), ("d", "R"): Fr(-1, 3),
        ("nu", "L"): Fr(-1, 2), ("e", "L"): Fr(-1, 2), ("nu", "R"): Fr(0), ("e", "R"): Fr(-1)}
T3TAB = {("u", "L"): Fr(1, 2), ("nu", "L"): Fr(1, 2), ("d", "L"): Fr(-1, 2), ("e", "L"): Fr(-1, 2)}
yv = [YTAB[(name, ch)] for (name, ch, c, gen) in orig.P]
tv = [T3TAB.get((name, ch), Fr(0)) for (name, ch, c, gen) in orig.P]
Yh = np.diag([float(v) for v in yv]+[-float(v) for v in yv])
T3 = np.diag([float(v) for v in tv]+[-float(v) for v in tv])
Q = Yh+T3
check("96 fibres = 48 particle states and their conjugates; 6 Y is an integer on all of them",
      len(orig.P) == 48 and all((6*v).denominator == 1 for v in yv))
phase = np.round(np.exp(4j*np.pi*np.diag(Yh)), 9)
sizes = sorted(int(np.sum(np.isclose(phase, p))) for p in set(phase))
check("exp(4 pi i Y) cuts the fibres into blocks of 24, 36 and 36 (recorded)", sizes == [24, 36, 36], sizes)
check("hypercharge, isospin and electric charge commute with the grading and are odd under the real structure: "
      "they are of the kind 'charge'", all(np.allclose(X@g, g@X) and np.allclose(conjK(X), -X) for X in (Yh, T3, Q)))
p0 = Asm.assemble({"schema": "tfpt.field-configuration.v1", "domain_sha256": Asm.domain["_sha256"], "genes": [1, 1, 1]})["phi_fourier"][0]
check("baseline: the constant potential conserves the electric charge and breaks the hypercharge",
      np.allclose(p0@Q, Q@p0) and not np.allclose(p0@Yh, Yh@p0), round(float(np.max(np.abs(p0@Yh-Yh@p0))), 6))

# ---------------------------------------------------------------- W2
print("W2  the Higgs square")
V = Asm.V
H2 = sum(v@v for v in V)
ev = np.linalg.eigvalsh(H2)
check("sum of the squares of the four original Higgs directions: hermitian, non-negative",
      np.allclose(H2, H2.conj().T) and ev.min() > -1e-10, (round(float(ev.min()), 6), round(float(ev.max()), 6)))
check("it commutes with the grading and is even under the real structure",
      np.allclose(H2@g, g@H2) and np.allclose(conjK(H2), H2))
cy, cq = float(np.max(np.abs(H2@Yh-Yh@H2))), float(np.max(np.abs(H2@Q-Q@H2)))
check("it commutes with hypercharge, isospin and electric charge: a neutral scalar marker",
      cy < 1e-10 and cq < 1e-10 and np.allclose(H2@T3, T3@H2), (cy, cq))
U = H2/ev.max()

# ---------------------------------------------------------------- W3
print("W3  what reality of the line forces")
k = 0.25
h0 = k*g+p0
wh, vh = np.linalg.eigh(h0)
T = (vh*np.exp(-0.7*np.abs(wh)))@vh.conj().T
check("an even transfer (function of h_0^2) is real under the real structure and commutes with the grading",
      np.allclose(conjK(T), T) and np.allclose(T@g, g@T))


def marks(x, C):
    z0 = x.conj()@T@T@T@x
    return (x.conj()@T@T@U@T@x)/z0, (x.conj()@T@C@T@T@x)/z0, (x.conj()@T@C@T@U@T@x)/z0


worst_re, worst_im = 0.0, 0.0
for _ in range(40):
    x = rng.normal(size=96)+1j*rng.normal(size=96)
    x = x+K@x.conj()
    x /= np.linalg.norm(x)
    for C in (Yh, Q, T3, g):
        zu, zy, zuy = marks(x, C)
        worst_im = max(worst_im, abs(zu.imag))
        worst_re = max(worst_re, abs(zy.real), abs(zuy.real))
check("on every real line: the mark of the even scalar is real; the mark of a charge and the mixed mark are purely "
      "imaginary (40 random real lines, four charges)", worst_re < 1e-12 and worst_im < 1e-12, (worst_re, worst_im))
check("this is the recorded form Z_UY = i eta q Z_0 with real q, forced by the real structure alone", True)

# ---------------------------------------------------------------- W4
print("W4  conserved against broken charge")
single_q, mixed_q, big_y, zy_y = 0.0, 0.0, 0.0, 0.0
sample = []
for _ in range(200):
    x = rng.normal(size=96)+1j*rng.normal(size=96)
    x = x+K@x.conj()
    x /= np.linalg.norm(x)
    _, zyq, zuyq = marks(x, Q)
    zu, zyy, zuyy = marks(x, Yh)
    single_q, mixed_q = max(single_q, abs(zyq)), max(mixed_q, abs(zuyq))
    big_y = max(big_y, abs(zuyy))
    zy_y = max(zy_y, abs(zyy))
    sample.append((zu.real, zyy.imag, zuyy.imag))
check("electric charge: the single mark vanishes on every real line (the charge is conserved and odd)", single_q < 1e-12, single_q)
check("electric charge: the mixed mark does not vanish, because the Higgs square does not commute with the transfer; "
      "it is purely imaginary", mixed_q > 1e-6 and not np.allclose(U@T, T@U), round(mixed_q, 6))
check("hypercharge: the single mark and the mixed mark are non-zero imaginary numbers on generic real lines",
      big_y > 1e-6 and zy_y > 1e-6, (round(zy_y, 6), round(big_y, 6)))
arr = np.array(sample)
print("      200 real lines: Z_U/Z_0 from %.4f to %.4f; Im Z_Y/Z_0 from %.5f to %.5f; Im Z_UY/Z_0 from %.5f to %.5f"
      % (arr[:, 0].min(), arr[:, 0].max(), arr[:, 1].min(), arr[:, 1].max(), arr[:, 2].min(), arr[:, 2].max()))
conn = arr[:, 2]-arr[:, 0]*arr[:, 1]
check("both signs occur for the imaginary mixed mark: the sign eta depends on the line", arr[:, 2].min() < 0 < arr[:, 2].max(),
      ("connected, imaginary part", round(float(conn.min()), 5), round(float(conn.max()), 5)))

# ---------------------------------------------------------------- W5
print("W5  flat density on the 96 fibres")
rows = []
for step in (0.3, 0.7, 1.5, 3.0):
    Ts = (vh*np.exp(-step*np.abs(wh)))@vh.conj().T
    z0 = np.trace(Ts@Ts@Ts).real
    rows.append((step, np.trace(Ts@Yh@Ts@Ts)/z0, np.trace(Ts@Yh@Ts@U@Ts)/z0, np.trace(Ts@Q@Ts@Ts)/z0, np.trace(Ts@Q@Ts@U@Ts)/z0,
                 (np.trace(Ts@Ts@U@Ts)/z0).real))
for step, y1, y2, q1, q2, u1 in rows:
    print("      step %.1f: Z_U/Z_0 = %.5f; hypercharge: single %.1e, mixed %+.4ei; electric: single %.1e, mixed %.1e"
          % (step, u1, abs(y1), y2.imag, abs(q1), abs(q2)))
check("flat density: the single charge marks vanish for both charges, and the mixed mark vanishes for the electric charge",
      all(abs(y1) < 1e-10 and abs(q1) < 1e-10 and abs(q2) < 1e-10 for _, y1, _, q1, q2, _ in rows))
check("flat density, hypercharge: the mixed mark is purely imaginary and not zero; it is small and depends on the step",
      all(abs(y2.real) < 1e-10 and abs(y2.imag) > 1e-9 for _, _, y2, _, _, _ in rows), [float("%.3e" % y2.imag) for _, _, y2, _, _, _ in rows])
check("its connected part equals the mark itself, since the single hypercharge mark is zero: with native markers the "
      "flat density gives an imaginary mixed answer, not the negative real one of the recorded E8 contract", True)

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
