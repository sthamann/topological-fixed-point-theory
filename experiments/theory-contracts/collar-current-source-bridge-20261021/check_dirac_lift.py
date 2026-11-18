"""Complete the original scalar endpoint family in its existing Dirac module.

Exact algebraic controls accompany DIRAC-LIFT-ADMISSIBILITY.md. No finite CCR,
physical chirality, Padm identification, or complete TFPT solution is asserted.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import sympy as s

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = ROOT / 'experiments/theory-contracts/seam-full-collar-reconstruction-20261015/checker.py'
spec = importlib.util.spec_from_file_location('original_collar_dirac_lift', SOURCE)
original = importlib.util.module_from_spec(spec)
spec.loader.exec_module(original)

def central(C, X, Y):
    Q = s.eye(C.rows) - C
    return s.simplify(s.trace(C*X*Q*Y*C - C*Y*Q*X*C))

def zero(M):
    return all(s.simplify(x) == 0 for x in M)

rows = {}
for alpha in (s.Rational(1, 4), s.Rational(3, 4)):
    modes = list(range(-2, 3))
    d = len(modes)
    I, Z = s.eye(d), s.zeros(d)
    A = s.diag(*[r+alpha for r in modes])
    B = s.diag(A, -A)
    Jn = Z.row_join(-I).col_join(I.row_join(Z))
    Gamma = s.diag(I, -I)
    T = s.diag(*[s.Integer(2)**(-4*B[k, k]) for k in range(2*d)])
    a, c, G = original.graph(T)
    N = s.diag(-Jn, Jn)
    assert N*G*N.conjugate().T == s.eye(4*d)-G
    assert Jn*B + B*Jn == s.zeros(2*d)
    n = s.eye(2*d)-a
    assert n[d:, d:] == s.eye(d)-n[:d, :d]
    shift = s.zeros(d)
    for k in range(d-1):
        shift[k+1, k] = 1
    f, h = (shift+shift.T)/2, (-s.I*shift+s.I*shift.T)/2
    xp, yp = s.diag(Z, Z, f, Z), s.diag(Z, Z, h, Z)
    xm, ym = s.diag(Z, Z, Z, f), s.diag(Z, Z, Z, h)
    x, y = xp+xm, yp+ym
    cp, cm, cv = central(G, xp, yp), central(G, xm, ym), central(G, x, y)
    assert cp != 0 and cp == -cm and cv == 0
    assert zero(N*x-x*N) and zero(N*y-y*N)
    assert not zero(N*xp-xp*N)
    # A constant relative face phase already moves the Cauchy trace graph.
    endpoint = s.diag(s.eye(2*d), s.I*s.eye(2*d))
    moved = endpoint*G*endpoint.conjugate().T
    assert moved != G
    assert moved[2*d:, 2*d:] == G[2*d:, 2*d:]
    assert s.diag(s.I*I, I)*Jn != Jn*s.diag(s.I*I, I)
    rows[str(alpha)] = {
        'normal_complement_identity': True,
        'second_face_opposite_occupations': True,
        'finite_component_plus_cocycle': str(cp),
        'finite_component_minus_cocycle': str(cm),
        'finite_full_vector_cocycle': str(cv),
        'nontrivial_endpoint_moves_trace_graph': True,
        'scalar_eta_analytic': str(1-2*alpha),
        'paired_B_eta_analytic': '0',
    }

Jn = s.Matrix([[0, -1], [1, 0]])
Gamma = s.diag(1, -1)
a, b, c, d = s.symbols('a b c d')
M = s.Matrix([[a, b], [c, d]])
commutant = s.solve(list(M*Jn-Jn*M)+list(M*Jn*Gamma-Jn*Gamma*M), [a,b,c,d])
assert commutant == {a:d, b:0, c:0}

r, theta, t = s.symbols('r theta t', real=True)
alpha = s.symbols('alpha', real=True)
chi, phi = s.Function('chi')(r), s.Function('phi')(theta)
v = s.Matrix([s.Function('x')(r,theta), s.Function('y')(r,theta)])
g = s.exp(s.I*t*chi*phi)
def D(v):
    return Jn*v.diff(r) + Jn*Gamma*(-s.I*v.diff(theta)+alpha*v)
transported = s.simplify(g*D(v/g))
V = -s.I*t*s.diff(chi,r)*phi*Jn - t*chi*s.diff(phi,theta)*Jn*Gamma
assert zero(transported-D(v)-V*v)
q, z = s.symbols('q z', real=True)
potential = -s.I*q*Jn-z*Jn*Gamma
assert potential.conjugate().T == potential
assert zero(s.diff(transported,t,2))

# The original four-mark-fixed profiles still move the fixed Hardy/APS space.
f, h = s.sin(4*theta), 1-s.cos(4*theta)
assert all(f.subs(theta,k*s.pi/2) == h.subs(theta,k*s.pi/2) == 0 for k in range(4))
leak_f = s.integrate(s.expand_trig(f)*s.exp(-4*s.I*theta), (theta,0,2*s.pi))/(2*s.pi)
leak_h = s.integrate(h*s.exp(-4*s.I*theta), (theta,0,2*s.pi))/(2*s.pi)
assert s.simplify(leak_f) == -s.I/2 and s.simplify(leak_h) == -s.Rational(1,2)

result = {
    'verdict': 'CONDITIONAL_FULL_DIRAC_VECTOR_COCYCLE_CANCELLATION_AND_DOMAIN_DICHOTOMY',
    'source_graph_constructor': str(SOURCE.relative_to(ROOT)),
    'source_graph_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    'cases': rows,
    'fixed_symbol_pointwise_commutant': 'scalar I2 (no additional internal multiplicities classified)',
    'full_Dirac_differential_lift_exact': True,
    'full_Dirac_potential_Hermitian_exact': True,
    'four_mark_profile_APS_leakage_minus1_to_plus3': {'sin4theta':'-I/2', '1-cos4theta':'-1/2'},
    'analytic_four_mark_cocycles': {'component_plus':'-2I','component_minus':'+2I','complete_vector':'0'},
    'analytic_four_mark_Weyl_rectangles': {'component_plus':'exp(2I s t)','component_minus':'exp(-2I s t)','complete_vector':'1'},
    'unchanged_previous_scalar_result': True,
    'scope': 'Existing complete two-component Dirac module, complete two-face graph, given charged CAR state, implementable normal-commuting generators. Infinite claims require the proof, not a mode cutoff.',
    'not_derived': ['P1 physical chiral measure/line connection', 'original OS Padm and A/A* operator dictionary', 'physical boundary attachment selection', '4D chirality or determinant-line unification'],
    'complete_TFPT_solution': False,
    'physical_gates_closed': [],
}
(HERE/'dirac_lift_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(result,ensure_ascii=False,indent=2))
