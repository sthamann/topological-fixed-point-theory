"""Exact decomposition of the already computed original four-edge transport.

No physical sum prescription or source insertion is selected.
"""
from pathlib import Path
import hashlib
import json
import sympy as s

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / 'ACTUAL-DIAGONAL-SPINE-ENDPOINT-TRANSPORT-CHECK.json'
EXPECTED = '46787929dff0f7517949d225263fb4301bb4d3fe553b7bdfb4b937c6b6beae7d'
if hashlib.sha256(INPUT.read_bytes()).hexdigest() != EXPECTED:
    raise RuntimeError('Actual original endpoint input changed')
a = json.loads(INPUT.read_text())
def mat(x):
    return s.Matrix([[s.sympify(v) for v in row] for row in x])
def zero(x):
    return all(s.simplify(v) == 0 for v in x)
def check(ok, name):
    if not ok:
        raise RuntimeError(name)
T = [mat(x) for x in a['T_source_normalized_four']]
M = [mat(x) for x in a['M_four']]
H = mat(a['H_infinity'])
U = mat(a['U0'])
V = mat(a['U_infinity'])
P0 = mat(a['P0'])
E10 = s.zeros(3); E10[1,0] = 1
E01 = E10.T
S = sum(T,s.zeros(3))
check(zero(S + 6*s.sqrt(30)/5*E10), 'actual unweighted sum')
check(zero(S.adjoint()*H*S-8*P0), 'source-metric sum Gram')
check(all(zero(T[k]*M[k]-T[(k-1)%4]) for k in range(4)), 'absorption of actual based loop')
check(zero(sum((T[k]*M[k]*P0 for k in range(4)),s.zeros(3))-S), 'covariant word sum equals edge sum')
F = [s.simplify(sum((s.I**(-m*k)*T[k] for k in range(4)),s.zeros(3))/4) for m in range(4)]
ranks = [x.rank() for x in F]
weights = [s.simplify(s.trace(x.adjoint()*H*x)) for x in F]
check(ranks==[1,2,1,2], 'all four character channels survive')
check(weights==[s.Rational(1,2),1,s.Rational(1,2),1], 'physical source-metric character weights')
check(all(zero(T[k]-sum((s.I**(m*k)*F[m] for m in range(4)),s.zeros(3))) for k in range(4)), 'full Fourier reconstruction')
check(all(zero(V*F[m]*U.inv()-s.I**m*F[m]) for m in range(4)), 'typed quarter characters')
check(all(zero(F[m].adjoint()*H*F[n]) or s.simplify(s.trace(F[m].adjoint()*H*F[n]))==0 for m in range(4) for n in range(4) if m!=n), 'Hilbert-Schmidt character orthogonality')

# The first nontrivial inverse-square-root response need not be a holonomy
# tangent.  The Hermitian unitary-holonomy tangency test is independent of
# any later choice of a bounded additive V_adm field slot.
x,y=s.symbols('x y',real=True)
c=x+s.I*y
B=c*M[0]*P0
unitarity_tangent=s.simplify(M[0].adjoint()*B+B.adjoint()*M[0])
check(zero(unitarity_tangent-2*x*P0), 'raw Hardy word is not a generic unitary holonomy tangent')
check(s.simplify(s.trace(M[0].inv()*B)-c)==0, 'SU3 determinant tangent')
# A concrete original-loop packaging retains the loop after forgetting an
# unweighted edge sum.  All four returned loops have order3 and rank3.
check(all(zero(M[k]**3-s.eye(3)) and M[k].rank()==3 for k in range(4)), 'all actual loop channels retained')
result={
 'scope':'Actual four-edge algebra only. Edge average not promoted to original P1 consumer; V_adm not selected.',
 'input_sha256':EXPECTED,
 'unweighted_sum_rank':S.rank(),
 'unweighted_sum_source_metric_gram':'8 P0',
 'actual_loop_absorption':[True]*4,
 'fourier_channel_ranks':ranks,
 'fourier_source_metric_squared_norms':[str(w) for w in weights],
 'all_four_characters_reconstruct_actual_edges':True,
 'unitary_holonomy_tangent_defect':'2 Re(kappa) P0',
 'su3_determinant_tangent':'kappa; nonzero for any nonzero Hardy coefficient',
 'meaning':'Invariant scalar/linear compression loses channels. Whole covariant path/loop functor is not this compression.',
 'checks_passed':12,
}
output=ROOT/'actual-spine-fourier-consumers-results.json'
output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
