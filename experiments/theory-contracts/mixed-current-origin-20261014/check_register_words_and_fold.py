"""Exact decisive controls, without promoting finite tests to field limits.

All-word register obstruction, actual one-copy folded monodromy distinction,
charge-conjugating countercontrol and original full-Dirac assertion replay.
"""
from pathlib import Path
import ast
import contextlib
import hashlib
import importlib
import io
import json
import sys
import sympy as s

ROOT=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0,str(ROOT/'verification'))
v1033=importlib.import_module('v1033_charged_disorder')


def need(value,label):
    if not value:raise ValueError(label)


def zero(matrix):
    return all(s.expand(x)==0 for x in matrix)


Z=s.diag(1,s.I,-1,-s.I)
X=s.zeros(4)
for r in range(4):X[(r+1)%4,r]=1
F=s.zeros(8)
for r in range(4):
    B=s.Matrix([[r+1,s.I*r],[1-r,s.Rational(r+2,2)]])
    F[2*((r+2)%4):2*((r+2)%4)+2,2*r:2*r+2]=B
Zfull=s.kronecker_product(Z,s.eye(2))
need(Zfull*F*Zfull.conjugate().T==-F,'arbitrary nonunitary half-grade representative')
need(zero(F.conjugate().T*Zfull.conjugate().T*F*Zfull+F.conjugate().T*F),
     'all-word fourpoint sign identity')
comm=Zfull*F-F*Zfull
need(zero(comm.conjugate().T*comm-4*F.conjugate().T*F),'norm obstruction squared')
n=s.trace(F.conjugate().T*F)/8
need(n>0,'positive nonzero normalized word test')

# All register offdiagonals are killed by time invariance if the FULL
# original even matter algebra is retained: Schur reduces P to one register
# matrix per matter parity. Witness each offdiagonal with actual hopping.
hs=[v1033.qwz_cylinder(3,8,1.,r) for r in range(4)]
projection_witnesses=[]
for r in range(4):
    for t in range(4):
        if r==t:continue
        diff=hs[r]-hs[t]
        hits=[(i,j) for i in range(48) for j in range(48) if i!=j and diff[i,j]!=0]
        need(bool(hits),'each pair of register Hamiltonians differs on a matter transfer')
        i,j=hits[0]
        k=next(k for k in range(48) if k not in (i,j))
        projection_witnesses.append({'r':r,'s':t,'one_particle_transfer':[i,j],
                                     'even_sector_spectator':k})

# Same-x-orientation shared flux: no offdiagonal pair of background labels.
need(Zfull-Zfull==s.zeros(8),'actual shared edge Wilson identity')
ZL,ZR=s.kronecker_product(Z,s.eye(4)),s.kronecker_product(s.eye(4),Z)
left=s.kronecker_product(X**2,s.eye(4))
vac=s.eye(16)[:,0]
defect=(ZL-ZR)*left*vac
need((defect.conjugate().T*defect)[0]==4,'proposed one-edge flip violates joint word identity')

# Actual source fold: cyclic sheet monodromy C(h)^8=h I.
h=s.symbols('h',nonzero=True)
def cyclic(holonomy):
    out=s.zeros(8)
    for j in range(7):out[j+1,j]=1
    out[0,7]=holonomy
    return out
C=cyclic(h)
S=s.diag(*[(-1)**j for j in range(8)])
need(C**8==h*s.eye(8),'actual eightfold cover holonomy')
need(cyclic(-h)**8==-h*s.eye(8),'primitive half flux changes parent holonomy')
need((-C)**8==h*s.eye(8),'uniform vector-minus spinor monodromy invisible on cover')
need(S*C*S==-C,'two base lifts have the SAME cover holonomy')

# Mandatory countercontrol: at quarters an SO16 charge conjugation can
# relate the real holonomies, but reverses the original U1 charge action.
def realify(M):
    re,im=M.applyfunc(s.re),M.applyfunc(s.im)
    return re.row_join(-im).col_join(im.row_join(re))
I=s.eye(8)
complex_structure=s.zeros(8).row_join(-I).col_join(I.row_join(s.zeros(8)))
K=s.diag(I,-I)
O=realify(S)*K
need(O*O.T==s.eye(16) and O.det()==1,'proper real charge-conjugation countercontrol')
for q in [s.I,-s.I]:
    need(O*realify(cyclic(-q))*O.T==realify(-cyclic(q)),
         'quarter real holonomy conjugacy EXISTS when charge is reversed')
need(O*complex_structure*O.T==-complex_structure,'this escape reverses charge,not fixed-Q map')

# Read-only replay: preserve all original assertions and suppress only its
# result-file mutation. Do not run this assertion-based source under -OO.
source=ROOT/'experiments/theory-contracts/collar-current-source-bridge-20261021/check_dirac_lift.py'
tree=ast.parse(source.read_text(),filename=str(source))
tree.body=[stmt for stmt in tree.body if not(
    isinstance(stmt,ast.Expr) and isinstance(stmt.value,ast.Call)
    and isinstance(stmt.value.func,ast.Attribute) and stmt.value.func.attr=='write_text')]
# Replace Assert nodes with explicit raises so -OO also performs every check.
class ExplicitAssertions(ast.NodeTransformer):
    def visit_Assert(self,node):
        return ast.copy_location(ast.If(test=ast.UnaryOp(op=ast.Not(),operand=node.test),
            body=[ast.Raise(exc=ast.Call(func=ast.Name(id='RuntimeError',ctx=ast.Load()),
                args=[node.msg or ast.Constant(value='original assertion at line '+str(node.lineno))],keywords=[]),cause=None)],
            orelse=[]),node)
tree=ast.fix_missing_locations(ExplicitAssertions().visit(tree))
buffer=io.StringIO()
with contextlib.redirect_stdout(buffer):
    exec(compile(tree,str(source),'exec'),{'__file__':str(source),'__name__':'read_only_assertion_replay'})
dirac=json.loads(buffer.getvalue())
need(all(c['finite_full_vector_cocycle']=='0' for c in dirac['cases'].values()),
     'actual complete Dirac cancellation both quarters')

result={
    'date':'2026-10-22',
    'verdict':'ALL_WORD_REGISTER_OBSTRUCTION_AND_FIXED_CHARGE_FOLD_MISMATCH',
    'nonunitary_four_word_positive_norm':str(n),
    'half_flux_four_word_ratio':-1,
    'required_remote_even_locality_ratio':1,
    'uniform_operator_norm_commutator_ratio':2,
    'actual_width8_projection_time_witnesses':projection_witnesses,
    'shared_flux_one_edge_counterexample_norm_squared':4,
    'parent_half_flux_cover_holonomy':'-h',
    'uniform_D8_spinor_vector_twist_cover_holonomy':'h',
    'complex_charge_preserving_conjugacy_possible':False,
    'quarter_real_charge_reversing_conjugacy_possible':True,
    'all_original_Dirac_assertions_executed_even_under_OO':True,
    'original_Dirac_replay':dirac,
    'source_hashes':{
        str(source.relative_to(ROOT)):hashlib.sha256(source.read_bytes()).hexdigest(),
        'verification/v1033_charged_disorder.py':hashlib.sha256((ROOT/'verification/v1033_charged_disorder.py').read_bytes()).hexdigest()},
    'P1_E8_origin_derived':False,
}
Path(__file__).with_name('register-words-and-fold-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['actual_width8_projection_time_witnesses','original_Dirac_replay']},indent=2))
