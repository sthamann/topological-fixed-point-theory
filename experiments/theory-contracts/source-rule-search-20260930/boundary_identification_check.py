"""Exact first-identification test for the original Dirac collar claim.

This is a principal-symbol/mode test, not a replacement TFPT model.
The exact flat Dirac operator is D=i*sigma1*partial_x+sigma2*K,
K=-i*partial_theta. Reflection (x,theta)->(-x,theta) has the
Dirac-preserving local lift R=sigma2. Antiperiodic spin structure
removes the zero mode. Positive-side L2 solutions select B=k*sigma3>0.
"""
import json
from pathlib import Path
import sympy as s

I=s.I
one=s.eye(2)
sx=s.Matrix([[0,1],[1,0]])
sy=s.Matrix([[0,-I],[I,0]])
sz=s.diag(1,-1)
cn=I*sx
R=sy
N=-I*cn
zero=s.zeros(2)

assert cn**2==-one and cn.H==-cn
assert R**2==one and R.H==R
assert R*cn*R==-cn
assert R*sy*R==sy

records=[]
for k in [s.Rational(-3,2),s.Rational(-1,2),s.Rational(1,2),s.Rational(3,2)]:
    B=k*sz
    C=s.sign(k)*sz
    P=(one+C)/2
    assert cn*B+B*cn==zero
    assert P**2==P and P.H==P
    assert s.exp(-s.Symbol('x',positive=True)*abs(k))*P == s.exp(-s.Symbol('x',positive=True)*B)*P
    assert R*P*R==one-P
    assert R*C+C*R==zero
    assert cn*P*cn.H==one-P
    # Correct quarter-order operation from normal Clifford and polarization.
    Q=N*C
    assert Q.H*Q==one and Q**2==-one and Q**4==one
    assert C != R
    records.append({'k':str(k),'C':str(C),'R_equals_C':False,
                    'reflection_complements_P':True,'quarter_square':'-I'})

# Locality obstruction, independent of the particular chosen reflection lift:
# a geometric lift over the same boundary point is independent of tangent xi,
# whereas the Dirac Calderon polarization has C(-xi)=-C(xi).
assert sz != -sz
assert (sz-(-sz)).det()!=0

# A finite C-invariant and R-invariant subspace has equal C eigenspace ranks.
# Hence a 3+2 carrier needs a non-invariant compression / relative-index map.
# Proof: CN=-NC maps N:E_+ -> E_-; N^2=I gives its inverse.
# Therefore the two restrictions are bijections and their dimensions agree.

out={
    'verdict':'ORIGINAL_GEOMETRIC_EQUALITY_NOT_DERIVED; LITERAL_LOCAL_EQUALITY_FAILS',
    'operator':'D=i sigma1 partial_x+sigma2(-i partial_theta)',
    'spin_structure':'antiperiodic; k in Z+1/2; no zero-mode convention',
    'exact_modes':records,
    'correct_general_relation':'gamma_n P_C gamma_n^* = I-P_C (orthogonal Dirac Calderon projection)',
    'local_symbol_obstruction':'C(-xi)=-C(xi); a local geometric reflection lift has no xi dependence',
    'quarter_operation':'N=-i gamma_n, Q=N(2P_C-I), Q^2=-I; derived boundary grading, not an independent TFPT clock',
    'finite_carrier_constraint':'if both involutions preserve E, dim E_+=dim E_-; 3+2 requires explicit compression/index morphism',
    'scope':'The archived theorem identifies R_tau with C without defining its transported action. The exact test excludes a literal local geometric lift. A nonlocal transported R_tau needs a separate definition and proof.',
    'not_proved':['source-state selection','canonical E8 lattice action','physical quarter time','full TFPT contradiction','fundamental full solution']
}
path=Path(__file__).with_name('boundary_identification_results.json')
path.write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
print(json.dumps(out,ensure_ascii=False))
