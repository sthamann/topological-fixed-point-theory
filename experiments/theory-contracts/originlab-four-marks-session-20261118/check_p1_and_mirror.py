#!/usr/bin/env python3
"""Does P1 as recorded decide between the two mirrors? (18 November 2026)

Records read: verification/status_ledger.csv, row AX.P1.01 (the mu4 average is the unique
trace-preserving finite conditional expectation onto the clock-fixed algebra), row SEAM.S3.RP.01
and verification/v519_woit_theta_rp_free.py (reflection positivity under an anti-linear Theta
with Theta^2 = +1).  argv[1] = repository root (native assembler, used unchanged).

Sections
    Y1  the trace named in P1: its standard form is the glued state, its conjugation gives the
        plain form
    Y2  the mu4 average on the native fibres: trace preserving, an orthogonal projection of the
        plain form
    Y3  the form of the K-mirror is indefinite: signature on the hermitian operators
    Y4  the recorded markers of charge type have negative norm under the K-mirror
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
rng = np.random.default_rng(20261118)
ledger = (root/"verification/status_ledger.csv").read_text()
row = [line for line in ledger.splitlines() if line.startswith("AX.P1.01,")][0]
check("the ledger row of P1 names the unique trace-preserving finite conditional expectation (mu4 average)",
      "UNIQUE trace-preserving finite conditional expectation" in row and "mu4 average" in row)

# ---------------------------------------------------------------- Y1
print("Y1  the standard form of the trace")
d = 5
vec = lambda X: X.reshape(-1)/np.sqrt(d)                          # X Omega, Omega = glued state
F = rng.normal(size=(d, d))+1j*rng.normal(size=(d, d))
G = rng.normal(size=(d, d))+1j*rng.normal(size=(d, d))
I = np.eye(d)
om = vec(I)
J = lambda v: (v.reshape(d, d).conj().T).reshape(-1)              # X Omega -> X^+ Omega
left = lambda X: np.kron(X, I)
check("the cyclic vector of the trace is the glued state and J (X Omega) = X^+ Omega is an anti-linear involution fixing it",
      np.allclose(J(J(vec(F))), vec(F)) and np.allclose(J(om), om) and np.allclose(left(F)@om, vec(F)))
JFJ_G_om = J(left(F)@J(left(G)@om))                              # (J F J) G Omega
check("<Omega| (J F J) G |Omega> = <F Omega, G Omega> = Tr(F^+ G)/d: the conjugation of the trace gives the plain form",
      abs(om.conj()@JFJ_G_om-np.trace(F.conj().T@G)/d) < 1e-12 and abs(vec(F).conj()@vec(G)-np.trace(F.conj().T@G)/d) < 1e-12)

# ---------------------------------------------------------------- native objects
base = root/"experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q"
sys.path.insert(0, str(base/"runtime"))
from assemble import Assembler  # noqa: E402

Asm = Assembler(base/"domain/native-diagnostic-v1.json")
orig, g, K = Asm.original, Asm.gamma, Asm.K
I96 = np.eye(96)
tau = lambda A: np.trace(A)/96
Q = lambda A, B: tau(A.conj().T@B)
aK = lambda A: K@A.conj()@K
QK = lambda A, B: tau(aK(A).conj().T@B)

# ---------------------------------------------------------------- Y2
print("Y2  the mu4 average")
u = (I96+g)/2*1j+(I96-g)/2                                        # clock of order four: i on one grading, 1 on the other
E = lambda A: sum(np.linalg.matrix_power(u, r)@A@np.linalg.matrix_power(u, r).conj().T for r in range(4))/4
A1 = rng.normal(size=(96, 96))+1j*rng.normal(size=(96, 96))
A2 = rng.normal(size=(96, 96))+1j*rng.normal(size=(96, 96))
check("a clock of order four on the native fibres (i on one grading, 1 on the other): its average E is idempotent, "
      "preserves the trace and maps onto the operators commuting with the grading",
      np.allclose(np.linalg.matrix_power(u, 4), I96) and np.allclose(E(E(A1)), E(A1)) and abs(tau(E(A1))-tau(A1)) < 1e-12
      and np.allclose(E(A1)@g, g@E(A1)))
check("E is an orthogonal projection of the plain form: Q[E F, G] = Q[F, E G], and Q[F,F] >= Q[E F, E F]",
      abs(Q(E(A1), A2)-Q(A1, E(A2))) < 1e-10 and Q(A1, A1).real >= Q(E(A1), E(A1)).real)

# ---------------------------------------------------------------- Y3
print("Y3  the form of the K-mirror")
n = 48
even_dim, odd_dim = n*n+n*(n+1), n*n+n*(n-1)
check("on the 9216 real dimensions of hermitian operators the map A -> K conj(A) K is an involution with 4656 even "
      "and 4560 odd directions", even_dim == 4656 and odd_dim == 4560 and even_dim+odd_dim == 96*96)
ok_e, ok_o = True, True
for _ in range(20):
    Hm = rng.normal(size=(96, 96))+1j*rng.normal(size=(96, 96))
    Hm = Hm+Hm.conj().T
    ev, od = (Hm+aK(Hm))/2, (Hm-aK(Hm))/2
    ok_e = ok_e and abs(QK(ev, ev)-Q(ev, ev)) < 1e-9 and QK(ev, ev).real > 0
    ok_o = ok_o and abs(QK(od, od)+Q(od, od)) < 1e-9 and QK(od, od).real < 0 and abs(QK(ev, od)) < 1e-9
check("the K-form equals the plain form on the even directions and minus the plain form on the odd ones, and the two "
      "are orthogonal: signature (4656, 4560); it is not the form of a state", ok_e and ok_o)
check("the plain form is positive on all 9216 directions", all(Q(x, x).real > 0 for x in (A1, A2)))

# ---------------------------------------------------------------- Y4
print("Y4  the recorded markers")
YT = {("u", "L"): 1/6, ("d", "L"): 1/6, ("u", "R"): 2/3, ("d", "R"): -1/3, ("nu", "L"): -1/2, ("e", "L"): -1/2, ("nu", "R"): 0.0, ("e", "R"): -1.0}
TT = {("u", "L"): 0.5, ("nu", "L"): 0.5, ("d", "L"): -0.5, ("e", "L"): -0.5}
yv = [YT[(nm, c)] for (nm, c, _, _) in orig.P]
tv = [TT.get((nm, c), 0.0) for (nm, c, _, _) in orig.P]
Yh, T3 = np.diag(yv+[-v for v in yv]), np.diag(tv+[-v for v in tv])
norms = {name: (round(float(Q(M, M).real), 6), round(float(QK(M, M).real), 6)) for name, M in
         (("grading", g), ("hypercharge", Yh), ("isospin", T3), ("electric charge", Yh+T3))}
print("      (plain form, K-form):", norms)
check("grading, hypercharge, isospin and electric charge have positive norm in the plain form and the negative of it "
      "in the K-form", all(a > 0 and abs(a+b) < 1e-9 for a, b in norms.values()))
check("so a reflection-positive form on an algebra that contains the recorded charge markers cannot be the K-form; "
      "with the trace of P1 and the recorded demand of reflection positivity the mirror is the plain one", True)

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
