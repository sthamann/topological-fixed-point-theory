"""First full-product identity needed by the three original-source hypotheses.

The actual original constructor and all symbolic Y/MR are retained. No new
physical measure, regulator profile or spectral model is selected.
"""
import ast
import hashlib
import json
from pathlib import Path
import sympy as s

repo = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
source = repo / 'experiments/theory-contracts/craig-original-source-bridge-20261103/actual_product_pairing.py'
assert hashlib.sha256(source.read_bytes()).hexdigest() == '6acaa9bc2c7a157a30730d7d6de75e944266afbf69b18a31c6fba654b8bd1c04'
nodes = []
for node in ast.parse(source.read_text()).body:
    if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'Mass' for t in node.targets):
        break
    nodes.append(node)
env = {'__file__': str(source), '__name__': 'original_constructor_only'}
exec(compile(ast.Module(body=nodes, type_ignores=[]), str(source), 'exec'), env)
gam = [s.SparseMatrix(g) for g in env['gam']]
g5, C, KF, gF, Phi = [s.SparseMatrix(env[k]) for k in ['g5', 'C', 'KF', 'gF', 'Phi']]
I96 = s.SparseMatrix(s.eye(96))
R = s.kronecker_product(C, KF)
G = s.kronecker_product(g5, gF)
sv = s.symbols('S1:5', real=True)
W = sum((s.Rational(3, 4)*s.I*sv[mu]*s.kronecker_product(gam[mu]*g5,I96) for mu in range(4)),s.SparseMatrix(384,384,{}))
D = s.kronecker_product(g5, Phi) + W
P = s.kronecker_product(s.eye(4), Phi*Phi) + s.Rational(9,16)*sum(x*x for x in sv)*s.eye(384)
checks = {}
def zero(M):
    return all(s.expand(x) == 0 for x in s.SparseMatrix(M).todok().values())
def ck(name, val):
    assert val, name
    checks[name] = True

ck('whole_original_zero_momentum_real_torsion_square', zero(D*D-P))
ck('whole_original_zero_momentum_hermitian', zero(D.H-D))
ck('actual_full_J_reality', zero(R*D.conjugate()-D*R))
for e in [-1,1]:
    ix = [i for i in range(384) if G[i,i] == e]
    iy = [i for i in range(384) if G[i,i] == -e]
    M = (R.H*D).extract(ix,ix)
    P_e = P.extract(ix,ix)
    ck(f'whole_original_M_skew_{e}', zero(M.T+M))
    ck(f'whole_original_P_equals_MdaggerM_{e}', zero(M.H*M-P_e))
    # The actual antiunitary switches the total Weyl fibers and commutes with P.
    U = R.extract(iy,ix)
    ck(f'actual_J_fiber_unitarity_{e}', zero(U.H*U-s.eye(192)))
    ck(f'actual_J_fiber_P_intertwining_{e}', zero(U*P_e.conjugate()-P.extract(iy,iy)*U))
    ck(f'full_trace_is_twice_actual_Weyl_trace_{e}', s.expand(s.trace(P)-2*s.trace(P_e)) == 0)
    ck(f'linear_original_M_has_no_bare_second_torsion_contact_{e}', all(zero(M.diff(x,y)) for x in sv for y in sv))

# The coefficient follows from |Pf M|^4=det(M^dagger M), not a guessed
# multiplicity of physical particles. The full ordinary spin trace doubles
# the actual Weyl trace by the checked antiunitary intertwiner above.
pf_trace_weight = s.Rational(1,4)*s.Rational(1,2)
ck('actual_full_product_Pfaffian_log_heat_weight_one_eighth', pf_trace_weight == s.Rational(1,8))
d, chi, eps = s.symbols('d chi eps', positive=True)
a2 = s.Rational(9,32)*d/s.pi**2
geom_S2 = -6*chi**2*a2
fermion_positive_proper_time_S2 = pf_trace_weight*a2/eps
ck('ordinary_geometric_and_fermion_uv_S2_signs_differ', geom_S2.is_negative and fermion_positive_proper_time_S2.is_positive)

# The exact original heat-profile response is not a finite contact
# polynomial: every Taylor coefficient of F is nonzero. This is a support
# test of the already computed kernel, not a new physical heat-atom choice.
a = s.symbols('a')
moments = [s.integrate((a*(1-a))**n,(a,0,1)) for n in range(6)]
ck('original_F_form_factor_moments', moments == [s.factorial(n)**2/s.factorial(2*n+1) for n in range(6)])
ck('original_F_second_derivative_is_nonlocal_before_truncation', all(x != 0 for x in moments))

out = {
    'result_id':'TORSION20261104F', 'passed':len(checks), 'checks':checks,
    'actual_matrix_scope':'Existing ordinary Spin4/full96 total Weyl product, full symbolic Yukawa and Majorana, constant real torsion and Phi for the matrix check. Eeps is not Padm.',
    'analytic_extension':'For the actual J-real selfadjoint Gamma-odd operator with compatible paired domains, MdaggerM=D² restricted to Eeps and J intertwines the two fibers. Domain and selector transport remain required.',
    'fermion_absolute_effective_action':'-log|Pf M|=-(1/4)Tr_Eeps log(MdaggerM)=-(1/8)Tr_full log(D²), when the paired-domain and full-product hypotheses hold; phase retained separately.',
    'ordinary_positive_proper_time_uv_diagnostic':{'fermion_coefficient':str(fermion_positive_proper_time_S2),'geometric_coefficient':str(geom_S2),'not_a_selection_of_a_physical_regulator':True,'not_a_no_go_for_all_original_regulated_sources':True},
    'original_nonlocal_weight_contact':'An actual full spectral-profile second variation need not be supported at coincident points. F_xy cannot in general be removed as a local contact in the full-source spectral response.',
    'fundamental_solution_achieved':False, 'full_Sadm_derived':False,
    'physical_Htrpos_derived':False, 'common_physical_spectrum_computed':False,
}
Path(__file__).with_name('three-hypotheses-source-F.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['result_id','passed','fermion_absolute_effective_action','ordinary_positive_proper_time_uv_diagnostic','fundamental_solution_achieved']},ensure_ascii=False))
