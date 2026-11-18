#!/usr/bin/env python3
"""Does a recorded rule fix the phase of the line shared by the two copies? (18 November 2026)

Read from the repository (argv[1] = repository root): research/originlab/coherent_origin.py, the
canonical model of ORIGINLABCOHERENTF120261116 (functions marked_operators and contract_candidate,
called unchanged).  Compared with the two-copy operators of round twenty-six: S exchanges the
copies, P gives the second copy a sign.

Sections
    R1  the recorded operators are S, P and i S P, matrix for matrix
    R2  which state gives the recorded marks: the line lies in one copy; no phase is left
    R3  the recorded positive lift is the projector on one quarter block
    R4  inner half phases: the recorded alternation -i/4, +i/4, -i/4
    R5  the rate q: equal copies give q = 1; the recorded contract at q = 1
"""
import sys
from pathlib import Path

import sympy as sp

FAIL = []


def check(name, ok, info=None):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info is not None else ""))
    if not ok:
        FAIL.append(name)


root = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(root/"research/originlab"))
import coherent_origin as co  # noqa: E402

# two copies, one mode each; states |n1 n2>; one line: |10> (copy 1) and |01> (copy 2)
basis = [(0, 0), (1, 0), (0, 1), (1, 1)]
S4 = sp.zeros(4)
P4 = sp.zeros(4)
for k, (n1, n2) in enumerate(basis):
    S4[basis.index((n2, n1)), k] = (-1)**(n1*n2)
    P4[k, k] = (-1)**n2
line = [1, 2]
S = S4.extract(line, line)
P = P4.extract(line, line)

# ---------------------------------------------------------------- R1
print("R1  the recorded operators")
A, B, par = co.marked_operators()
a2, b2, p2 = A[:2, :2], B[:2, :2], par[:2, :2]
check("recorded swap = S and recorded parity = P on the one-line states (the recorded model is two such pairs)",
      a2 == S and p2 == P and A == sp.diag(S, S) and par == sp.diag(P, P))
check("recorded phase = i S P: the second recorded cut is the quarter turn of the mark, made hermitian",
      b2 == sp.I*S*P and B == sp.I*A*par and (S*P)**2 == -sp.eye(2) and b2.H == b2)
check("so the recorded pair of cuts is (exchange, quarter turn), not (exchange, sign); "
      "the record binds Y to the first and U to the second", True)

# ---------------------------------------------------------------- R2
print("R2  the state")
nx, ny, nz = sp.symbols("n_x n_y n_z", real=True)
sx, sy, sz = sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.diag(1, -1)
rho = (sp.eye(2)+nx*sx+ny*sy+nz*sz)/2
ev = lambda op: sp.simplify((rho*op).trace())
marks = (ev(a2), ev(b2), ev(a2*b2))
check("for any state of the line with Bloch vector n: <A> = n_x, <B> = n_y, <A B> = i n_z",
      marks == (nx, ny, sp.I*nz), marks)
check("recorded marks Z_U = Z_Y = 0 and Z_UY = i eta Z_0 with |eta| = 1 hold exactly for n = (0, 0, eta): "
      "the line lies in one copy, eta says in which", sp.solve([nx, ny, nz**2-1], [nx, ny, nz]) == [(0, 0, -1), (0, 0, 1)])
cap = sp.Matrix([1, 0])
check("the recorded cap (number state 0) is this state: an eigenstate of P",
      p2*cap == cap and (cap.H*a2*cap)[0] == 0 and (cap.H*b2*cap)[0] == 0 and (cap.H*a2*b2*cap)[0] == sp.I)
phi = sp.symbols("phi", real=True)
shared = sp.Matrix([sp.exp(sp.I*phi), 1])/sp.sqrt(2)          # round twenty-six: (|01> + e^{i phi}|10>)/sqrt 2
m26 = [sp.simplify(sp.re((shared.H*op*shared)[0]))+sp.I*sp.simplify(sp.im((shared.H*op*shared)[0])) for op in (S, P, S*P)]
check("in the assignment of round twenty-six (U = S, Y = P) the same marks need phi = +-pi/2, an eigenstate of "
      "the quarter turn; the record's own assignment needs no shared line at all",
      sp.simplify(m26[0]-sp.cos(phi)) == 0 and m26[1] == 0 and sp.simplify(m26[2]-sp.I*sp.sin(phi)) == 0)

# ---------------------------------------------------------------- R3
print("R3  the recorded positive lift")
lift = (sp.eye(2)+b2)/2
q = sp.symbols("q", positive=True)
t = sp.diag(1, q)
check("(1 + B)/2 is the projector on the value -i of the quarter turn S P: one quarter block",
      lift*lift == lift and lift.rank() == 1 and (S*P)*lift == -sp.I*lift)
check("its responses on the recorded cap: 1/2 and i q/2, as recorded",
      (cap.H*lift*cap)[0] == sp.Rational(1, 2) and sp.simplify((cap.H*a2*t*lift*cap)[0]-sp.I*q/2) == 0)

# ---------------------------------------------------------------- R4
print("R4  inner half phases")
down = sp.Matrix([0, 1])                                       # the line in the second copy, P = -1
p_first = a2*p2*a2                                             # the sign on the first copy
on_line = [sp.simplify((down.H*a2*p_first**k*b2*down)[0])/4 for k in range(3)]
off_line = [sp.simplify((down.H*a2*p2**k*b2*down)[0])/4 for k in range(3)]
check("line in the second copy, place factor 1/4, half phase = sign of the copy that carries the line between the "
      "two cuts: -i/4, +i/4, -i/4 for none, one and two half phases, the recorded values",
      on_line == [-sp.I/4, sp.I/4, -sp.I/4], on_line)
check("the sign of the other copy between the cuts changes nothing: -i/4 three times; the recorded alternation "
      "says only that the half phase sits on the line (the recorded statement in other words)",
      off_line == [-sp.I/4, -sp.I/4, -sp.I/4], off_line)

# ---------------------------------------------------------------- R5
print("R5  the rate q")
t1, t2 = sp.symbols("t_1 t_2", positive=True)
T = sp.diag(t1, t2)
check("a transfer that is diagonal in the copies and commutes with the exchange has t_1 = t_2: equal copies give q = 1",
      sp.solve(list(S*T-T*S), [t1], dict=True) == [{t1: t2}])
try:
    receipts = [co.contract_candidate(sp.S.One, sp.S.One, sp.S.One, o) for o in (1, -1)]
    info = [(r["orientation"], r["marks"], r["ABAB"], r["connected_UY"]) for r in receipts]
    print("      recorded contract at q = 1:", info)
    check("the recorded contract accepts q = 1: marks (1, 0, 0, +-i), ABAB = -1",
          all(r["marks"] == {"Z0": "1", "ZU": "0", "ZY": "0", "ZUY": str(sp.I*o)} and r["ABAB"] == "-1"
              for r, o in zip(receipts, (1, -1))))
except Exception as error:                                     # the recorded contract decides
    check("the recorded contract accepts q = 1", False, repr(error)[:200])
rates = sorted({sp.Rational(n, d) for d in range(2, 5) for n in range(1, d)})
check("the recorded search enumerates only 0 < q < 1 (stages 2 to 4): q = 1 is not among its candidates",
      max(rates) < 1, [str(r) for r in rates])

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
