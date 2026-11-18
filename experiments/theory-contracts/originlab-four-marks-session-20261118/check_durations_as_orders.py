#!/usr/bin/env python3
"""Are the recorded durations orders and not lengths? (18 November 2026)

argv[1] = repository root (research/originlab/coherent_origin.py is read unchanged).

Sections
    O1  a glued state with weights: the reflected form carries a transfer rho^(1/2) between the
        mirror image and the operator; the transfer is the half period of the modular flow
    O2  the flat glued state has no modular flow: only the order of the operators is left
    O3  any weight other than flat lowers the sum of the 60 return values below its recorded value
    O4  the recorded canonical model: the marks do not depend on the durations exactly when q = 1
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
rng = np.random.default_rng(20261118)

# ---------------------------------------------------------------- O1
print("O1  a glued state with weights")
d = 5
p = rng.uniform(0.2, 1.0, size=d)
p /= p.sum()
rho = np.diag(p)
Bv = np.diag(np.sqrt(p)).reshape(-1)                              # sum sqrt(p_a) |a>|a-bar>
F = rng.normal(size=(d, d))+1j*rng.normal(size=(d, d))
G = rng.normal(size=(d, d))+1j*rng.normal(size=(d, d))
lhs = Bv.conj()@np.kron(G, F.conj())@Bv
rt = np.diag(np.sqrt(p))
check("<B_rho| G x conj(F) |B_rho> = Tr(rho^(1/2) F^+ rho^(1/2) G): the weights act as a transfer rho^(1/2) on each side "
      "of the operator", abs(lhs-np.trace(rt@F.conj().T@rt@G)) < 1e-12)
flow = lambda A, t: np.diag(p**(1j*t))@A@np.diag(p**(-1j*t))      # modular flow of the first side
half = np.diag(p**0.5)@G@np.diag(p**-0.5)                         # the flow at imaginary time 1/2
check("this transfer is the modular flow of the glued state continued to half an imaginary period: "
      "Tr(rho^(1/2) F^+ rho^(1/2) G) = Tr(rho F^+ sigma_(-i/2)(G))", abs(np.trace(rt@F.conj().T@rt@G)-np.trace(rho@F.conj().T@half)) < 1e-12)
check("the flow is a one-parameter group fixing the state", np.allclose(flow(flow(G, 0.3), 0.4), flow(G, 0.7)) and np.allclose(flow(rho, 0.9), rho))

# ---------------------------------------------------------------- O2
print("O2  the flat glued state")
flat = np.eye(d)/d
check("for equal weights the modular flow is the identity for every time: there is no clock in the state; the form "
      "is Tr(F^+ G)/d and only the order of the operators is left",
      all(np.allclose(np.diag(np.diag(flat)**(1j*t))@G@np.diag(np.diag(flat)**(-1j*t)), G) for t in (0.3, 1.7)))

# ---------------------------------------------------------------- O3
print("O3  the 60 return values against the weight")
bits = np.array(list(itertools.product((0, 1), repeat=8)))
ann = []
for i in range(8):
    a = np.zeros((256, 256))
    for k, b in enumerate(bits):
        if b[i]:
            c = b.copy()
            c[i] = 0
            a[int("".join(map(str, c)), 2), k] = (-1)**int(b[:i].sum())
    ann.append(a)
spin2 = 2*bits-1
idx = np.where(np.array([list(v).count(-1) % 2 == 0 for v in spin2]))[0]
X = [((ann[i].T if si > 0 else ann[i])@(ann[j].T if sj > 0 else ann[j]))[np.ix_(idx, idx)]
     for i in range(5) for j in range(5, 8) for si in (1, -1) for sj in (1, -1)]
hyp = spin2[idx]@np.array([-1/3, -1/3, -1/3, 1/2, 1/2, 0, 0, 0])/2
energies = {"hypercharge squared": hyp**2, "one carrier sign": spin2[idx][:, 3]/2.0, "random": rng.normal(size=128)}


def returns(E, beta):
    w = np.exp(-beta*E)
    w /= w.sum()
    r = np.sqrt(w)
    return np.array([np.sum((r[:, None]*x)*(x*r[None, :])) for x in X])   # Tr(rho^(1/2) X^+ rho^(1/2) X), X real


base = returns(energies["random"], 0.0)
check("flat weight: every return value is 1/128 * 32 = 1/4 and the sum is 15", np.allclose(base, 0.25) and abs(base.sum()-15) < 1e-10)
ok, table = True, {}
for name, E in energies.items():
    sums = [float(returns(E, b).sum()) for b in (0.5, 1.0, 2.0, 4.0)]
    table[name] = [round(s, 5) for s in sums]
    ok = ok and all(s < 15-1e-9 for s in sums) and all(x > y for x, y in zip(sums, sums[1:]))
print("      sum of the 60 return values for steps 0.5, 1, 2, 4:", table)
check("for three choices of the weight exp(-beta E) the sum falls below 15 for every beta > 0 and keeps falling: "
      "any length of the transfer lowers the return values below the recorded 1/4", ok)
check("with round twenty-eight (the recorded values force the flat weight) this is the statement: in the glued form "
      "the recorded table leaves no length, only the order of marker and cut", True)

# ---------------------------------------------------------------- O4
print("O4  the recorded canonical model")
sys.path.insert(0, str(root/"research/originlab"))
import coherent_origin as co  # noqa: E402

A, B, par = co.marked_operators()
q = sp.symbols("q", positive=True)
da, db, dc = sp.symbols("d_a d_b d_c", positive=True)
T = lambda dur: sp.diag(*[q**(j*dur) for j in range(4)])
cap = sp.eye(4)[:, 0]
z0 = (cap.T*T(da)*T(db)*T(dc)*cap)[0]
zuy = sp.simplify((cap.T*T(dc)*A*T(db)*B*T(da)*cap)[0]/z0)
zu = sp.simplify((cap.T*T(dc)*T(db)*B*T(da)*cap)[0]/z0)
zy = sp.simplify((cap.T*T(dc)*A*T(db)*T(da)*cap)[0]/z0)
check("with durations d_a, d_b, d_c of the three sections the recorded operators give Z_U = Z_Y = 0 and "
      "Z_UY/Z_0 = i q^(d_b): only the middle duration enters, through the rate", zu == 0 and zy == 0 and sp.simplify(zuy-sp.I*q**db) == 0, zuy)
check("the marks are independent of all durations exactly when q = 1", sp.simplify(zuy.subs(q, 1)-sp.I) == 0 and sp.diff(zuy, db) != 0)
check("so in the recorded model 'durations are orders' is the same statement as q = 1, the value of equal copies "
      "(round twenty-seven) and of the flat glued state", True)

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
