"""Exact normal-fiber controls; read PROOF.txt for domains and full scope."""
import json
from pathlib import Path
import sympy as s

I=s.I
Z=s.diag(1,-1)
N=s.Matrix([[0,-1],[1,0]])
Ic=I*Z
C=I*N
b=s.symbols('b',real=True)
B=b*Z
q=s.symbols('q',positive=True)
a,c=q*q/(1+q*q),q/(1+q*q)
Ap=s.diag(a,1-a)
G=Ap.row_join(c*s.eye(2)).col_join((c*s.eye(2)).row_join(s.eye(2)-Ap))
Cf=s.kronecker_product(Z,C)
If=s.kronecker_product(s.eye(2),Ic)
T=s.zeros(4)
T[0,0]=T[2,1]=1
T[1,2]=I
T[3,3]=-I
scalar=s.Matrix([[a,c],[c,1-a]])
target=s.diag(1,1,1,1)
target[:2,:2]=scalar
target[2:,2:]=s.eye(2)-scalar
flip=s.zeros(4)
flip[:2,2:]=flip[2:,:2]=s.eye(2)
R=s.kronecker_product(s.Matrix([[0,1],[1,0]]),N*Z)

def eq(x,y):
    return all(s.cancel(e)==0 for e in x-y)

checks={
 'new_complex_structure':Ic*Ic==-s.eye(2),
 'real_structure_square':C*C==s.eye(2),
 'new_antilinearity':C*Ic==-Ic*C,
 'selfdual_B':C*B*C==-B,
 'field_embedding_is_fixed_by_real_structure':eq(C*s.Matrix([1,I]),s.Matrix([1,I])),
 'circle_generator_is_original_scalar':Ic*B==I*b*s.eye(2),
 'circle_action_preserves_real_form':C*Ic*B==Ic*B*C,
 'old_B_time_differs':C*I*B != I*B*C,
 'whole_graph_is_projection':eq(G*G,G),
 'oriented_graph_is_selfdual':eq(Cf*G*Cf,s.eye(4)-G),
 'charged_chart_is_unitary':T.conjugate().T*T==s.eye(4),
 'charged_chart_new_complex_linearity':If*T==T*I*s.diag(1,1,-1,-1),
 'charged_chart_real_flip':Cf*T==T*flip,
 'whole_graph_equals_scalar_charged_covariance':eq(T.conjugate().T*G*T,target),
 'reflection_is_involution':R*R==s.eye(4),
 'reflection_new_antilinearity':R*If*R==-If,
 'reflection_preserves_real_structure':R*Cf*R==Cf,
 'reflection_preserves_whole_graph':eq(R*G*R,G),
}
if not all(checks.values()):
    raise RuntimeError(checks)
result={'all_pass':True,'checks':checks,'scope':
 'Free original product collar in canonical graded Clifford typing; no original physical quantization selection or OS time theorem.'}
Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
