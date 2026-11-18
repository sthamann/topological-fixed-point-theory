import sympy as s
import json
M=s.Matrix([[1,0,0],[0,s.Rational(1,2),s.Rational(1,6)],[0,s.Rational(1,6),s.Rational(1,2)]])
I=s.eye(3); R=(I+M*M).inv(); C=M*M*R; X=M*R
P=(I-C).row_join(X).col_join(X.row_join(C))
perron=s.diag(1,0,0)
checks={
 'doubled_graph_is_projector':P*P==P,
 'doubled_graph_is_symmetric':P.T==P,
 'full_rank_three':P.rank()==3,
 'single_face_rank_three':C.rank()==3,
 'single_face_eigenvalues':set(C.eigenvals())=={s.Rational(1,2),s.Rational(4,13),s.Rational(1,10)},
 'single_face_not_perron':C!=perron and I-C!=perron,
 'calderon_reflection_squared_identity':(2*P-s.eye(6))**2==s.eye(6)
}
if not all(checks.values()):raise RuntimeError(checks)
print(json.dumps({'scope':'Given invertible finite family transfer used in the existing graph-Calderon formula; no physical source embedding inferred','checks':checks,'single_face_eigenvalues':['1/2','4/13','1/10']},indent=2))
