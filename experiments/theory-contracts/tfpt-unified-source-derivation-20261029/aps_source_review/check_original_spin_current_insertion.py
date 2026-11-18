"""R072: the actual arc/residue insertion, with the original fixed fold."""
from pathlib import Path
import json
import sympy as s

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
data = json.loads((REPO / "experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json").read_text())
M = s.Matrix([[s.Rational(x) for x in row] for row in data["explicit_Weyl_matrix"]])
zeta = (1 + s.I) / s.sqrt(2)
U = s.Matrix(8, 8, lambda r, k: s.expand_complex(zeta ** (r*k)) / s.sqrt(8))
checks = []


def require(name, condition):
    if not bool(condition):
        raise AssertionError(name)
    checks.append(name)


def simp(A):
    return A.applyfunc(s.simplify)


def hs2(A):
    return s.simplify(s.trace(A.conjugate().T * A))


require("actual fold is orthogonal", M.T*M == s.eye(8))
require("DFT core is unitary", simp(U.conjugate().T*U) == s.eye(8))
e = [s.eye(8)[:, j] for j in range(8)]
a, b = e[0]-e[1], e[1]-e[2]
rows = []
for index, root in enumerate([b, -b, a+b, -a-b]):
    q = M.T*root
    D = s.diag(*q)
    A = simp(U*D*U.conjugate().T)
    mean = sum(q)/8
    off = A-s.diag(*A.diagonal())
    require(f"{index}: all actual coordinates half signed", set(q) == {s.Rational(-1,2), s.Rational(1,2)})
    require(f"{index}: multiplicities2and6", sorted([list(q).count(s.Rational(-1,2)), list(q).count(s.Rational(1,2))]) == [2,6])
    require(f"{index}: level2", (q.T*q)[0] == 2)
    require(f"{index}: diagonal mean", A.diagonal() == s.ones(1,8)*mean)
    require(f"{index}: actual generator has offdiagonal norm3/2", hs2(off) == s.Rational(3,2))
    require(f"{index}: adjoint correct", A.conjugate().T == A)
    require(f"{index}: A squared onequarter", simp(A*A) == s.eye(8)/4)
    Vpi = -2*s.I*A
    require(f"{index}: pi source is unitary", simp(Vpi.conjugate().T*Vpi) == s.eye(8))
    require(f"{index}: pi offdiagonal norm6", hs2(Vpi-s.diag(*Vpi.diagonal())) == 6)
    # All alpha phases are diagonal row phases and preserve the entry moduli.
    for alpha in [s.Rational(1,4), s.Rational(3,4)]:
        phases = s.diag(*[s.exp(2*s.pi*s.I*alpha*r/8) for r in range(8)])
        Ap = phases*A*phases.conjugate().T
        require(f"{index},{alpha}: phased DFT diagonal mean", all(s.simplify(Ap[r,r]-mean) == 0 for r in range(8)))
        require(f"{index},{alpha}: phased offdiagonal norm3/2", hs2(Ap-s.diag(*Ap.diagonal())) == s.Rational(3,2))
    nonzero = [(r,t) for r in range(8) for t in range(8) if r != t and A[r,t] != 0]
    r,t = nonzero[0]
    P = e[t]*e[t].T
    require(f"{index}: actual cover singlearc support changes", simp(A*P-P*A) != s.zeros(8))
    R = s.eye(8)-2*P
    require(f"{index}: actual branch parity commutant test", simp(R*Vpi-Vpi*R) != s.zeros(8))
    rows.append({"root":list(map(str,root)), "q":list(map(str,q)), "sum_q":str(sum(q)),
                 "generator_offdiagonal_hs2":"3/2", "loop_pi_offdiagonal_hs2":"6",
                 "first_offdiagonal_entry":[r,t,str(A[r,t])], "nonzero_offdiagonal_entries":len(nonzero)})

# Exact all-profile formula, derived using q_j in {+1/2,-1/2}.
x = s.symbols("x", real=True)
for m in [s.Rational(-1,4),s.Rational(1,4)]:
    diagonal = s.cos(x/2)-2*s.I*m*s.sin(x/2)
    require(f"whole loop class at mean{m}", s.trigsimp(8*(1-s.expand_complex(diagonal*s.conjugate(diagonal)))-6*s.sin(x/2)**2) == 0)

# The scalar cover algebra and the residue Cartan algebra intersect in scalars.
# Independent symbolic Fourier reconstruction, without testing a sample profile.
d = s.symbols("d0:8")
B = U*s.diag(*d)*U.conjugate().T
linear_rows = [s.expand(B[r,t]) for r in range(8) for t in range(8) if r != t]
K,_ = s.linear_eq_to_matrix(linear_rows,d)
require("double diagonal equations rank7", K.rank() == 7)
require("double diagonal nullspace constant", K.nullspace() == [s.ones(8,1)])

result = {"result":"R072", "scope":"original fixed fold; smooth closed Cartan insertions, not sharp spinor vertex construction",
          "actual_roots":rows, "whole_class_offdiagonal_hs2":"6*sin(b/2)^2",
          "scalar_implementation_iff":"b in 2*pi*Z pointwise",
          "double_diagonal_rank":7, "passed":len(checks), "checks":checks}
(HERE/"ORIGINAL-SPIN-CURRENT-INSERTION.json").write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps({"passed":len(checks),"roots":4,"double_diagonal_rank":7,"result":"R072"}))
