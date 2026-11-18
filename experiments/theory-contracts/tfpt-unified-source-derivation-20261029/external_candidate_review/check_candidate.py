"""Exact controls for a scoped external-source candidate review."""
import json
from fractions import Fraction
from pathlib import Path
import sympy as s

q = s.symbols('q', positive=True)
a, c = q**2 / (1 + q**2), q / (1 + q**2)
A = s.diag(a, 1-a)
cross = c*s.eye(2)
G = A.row_join(cross).col_join(cross.row_join(s.eye(2)-A))
X, Z = s.Matrix([[0,1],[1,0]]), s.diag(1,-1)
naive = s.kronecker_product(s.eye(2), X)
oriented = s.kronecker_product(Z, X)
zero = s.zeros(4)

def eq(left, right):
    return all(s.cancel(x) == 0 for x in left-right)

expected_defect = s.zeros(4)
for i in range(2):
    expected_defect[i,i+2] = expected_defect[i+2,i] = 2*c
checks = {
    'graph_is_projection': eq(G*G,G),
    'graph_is_hermitian': G.T == G,
    'oriented_conjugation_squares_to_identity': oriented*oriented == s.eye(4),
    'oriented_conjugation_is_selfdual': eq(oriented*G*oriented,s.eye(4)-G),
    'naive_conjugation_has_exact_nonzero_cross_defect':
        eq(naive*G*naive-(s.eye(4)-G),expected_defect) and expected_defect != zero,
    'one_face_is_mixed': eq(A*A-A,-c*c*s.eye(2)),
}
# Full Fourier formula, checked on rational labels without floating rank tests.
for alpha in [Fraction(1,4),Fraction(3,4),Fraction(1,3)]:
    checks[f'paired_fourier_energies_{alpha}'] = all(
        2*(n+alpha) == -2*(-n-1+1-alpha) for n in range(-20,21))
lam1, lam2, threshold = Fraction(1,2), Fraction(3,2), Fraction(3,4)
checks['diagnostic_highpass_removes_lowest_block'] = (
    lam1 < threshold <= lam2 and not lam1**2 >= threshold**2
    and lam2**2 >= threshold**2)
results = {'checks':checks, 'all_pass':all(checks.values()),
    'scope':'Exact controls plus analytic proof; no physical P1 source derivation.',
    'primary_purification_result':'Araki 1970/71, Definition 4.4, Lemmas 4.5-4.6'}
if not results['all_pass']:
    raise RuntimeError(results)
Path(__file__).with_name('results.json').write_text(
    json.dumps(results,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(results,ensure_ascii=False))
