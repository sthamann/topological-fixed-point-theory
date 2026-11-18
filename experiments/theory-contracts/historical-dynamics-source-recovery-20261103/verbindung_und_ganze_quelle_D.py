"""Exact operator and finite Berezin controls for the existing torsion source.

The exterior algebra below checks a general algebraic identity, not a physical
four-mode model. No new geometry, coupling, state or regulator is selected.
"""
import ast
import hashlib
import itertools
import json
from pathlib import Path
import sympy as s

REPO = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
original = REPO/'experiments/theory-contracts/craig-original-source-bridge-20261103/actual_product_pairing.py'
assert hashlib.sha256(original.read_bytes()).hexdigest() == '6acaa9bc2c7a157a30730d7d6de75e944266afbf69b18a31c6fba654b8bd1c04'
# Execute only the actual original Clifford assignments. Importing that research
# script would run unrelated work and overwrite its output.
env = {'s': s}
for node in ast.parse(original.read_text()).body:
    if isinstance(node, ast.Assign) and any(isinstance(t,ast.Name) and t.id in {'I2','I4','pauli','gam','g5','C'} for t in node.targets):
        exec(compile(ast.Module(body=[node],type_ignores=[]),str(original),'exec'),env)
gam,g5,C = env['gam'],env['g5'],env['C']
I4=s.eye(4); z4=s.zeros(4); checks={}
def zero(M): return all(s.expand(x)==0 for x in M)
def ck(name,value):
    assert value,name
    checks[name]=True

# Entire real orthogonal connection fiber K_muab=-K_muba: 24 coordinates.
labels=list(itertools.product(range(4),list(itertools.combinations(range(4),2))))
K=s.symbols('k:24',real=True)
A=[]
for mu in range(4):
    A.append(sum((K[n]*gam[a]*gam[b]/2 for n,(m,(a,b)) in enumerate(labels) if m==mu),z4))
V=s.expand(-s.I*sum((gam[mu]*A[mu] for mu in range(4)),z4))
columns=s.Matrix.hstack(*[s.Matrix(V.diff(k)).reshape(16,1) for k in K])
defect=s.Matrix.hstack(*[s.Matrix((V-V.H).diff(k)).reshape(16,1) for k in K])
null=defect.nullspace(); N=s.Matrix.hstack(*null)
ck('whole_metric_connection_dirac_image_rank8',columns.rank()==8)
ck('formal_symmetry_has_four_independent_conditions',defect.rank()==4)
ck('symmetric_connection_fiber_dimension20',N.shape==(24,20))
ck('symmetric_dirac_image_rank4', (columns*N).rank()==4)
ck('dirac_invisible_connection_kernel_dimension16',len(columns.nullspace())==16)
for mu in range(4):
    ck(f'orthogonal_spin_connection_skew_{mu}',zero(A[mu]+A[mu].H))
    for c in range(4):
        rhs=sum((K[n]*( (gam[a] if b==c else z4) - (gam[b] if a==c else z4)) for n,(m,(a,b)) in enumerate(labels) if m==mu),z4)
        ck(f'connection_clifford_compatibility_requires_moving_tangent_connection_{mu}_{c}',zero(A[mu]*gam[c]-gam[c]*A[mu]-rhs))

S=s.symbols('S1:5',real=True);c0=s.Rational(3,4)
W=s.I*c0*sum((S[mu]*gam[mu]*g5 for mu in range(4)),z4)
ak=[]
for mu in range(4):
    ak.append(sum((s.LeviCivita(mu,a,b,nu)*S[nu]*gam[a]*gam[b]/8 for a in range(4) for b in range(4) for nu in range(4)),z4))
ck('actual_axial_spin_connection_contraction',zero(-s.I*sum((gam[mu]*ak[mu] for mu in range(4)),z4)-W))
ck('real_axial_hermitian',zero(W-W.H))
ck('real_axial_J_real',zero(C*W.conjugate()-W*C))
ck('real_axial_spin_odd',zero(g5*W+W*g5))
for nu in range(4):
    ck(f'actual_axial_extraction_{nu}',s.expand(s.trace(g5*gam[nu]*W)/(3*s.I)-S[nu])==0)
    ck(f'axial_does_not_commute_with_full_Clifford_{nu}',zero(W*gam[nu]-gam[nu]*W+s.Rational(3,2)*s.I*S[nu]*g5))
elementary=[]
for i in range(4):
    for j in range(4):
        E=s.zeros(4);E[i,j]=1;elementary.append(E)
commutator_map=s.Matrix.hstack(*[s.Matrix.vstack(*[(E*g-g*E).reshape(16,1) for g in gam]) for E in elementary])
ck('whole_Clifford_commutant_is_scalar_dimension_one',commutator_map.rank()==15)
axial=s.Matrix.hstack(*[s.Matrix(W.diff(x)).reshape(16,1) for x in S])
ck('axial_image_is_whole_symmetric_metric_connection_dirac_image',s.Matrix.hstack(columns*N,axial).rank()==4)
# W acts as identity on the finite algebra. Thus [W,pi(a)]=0 on the
# original full product, preserving its first-order condition exactly.
ck('central_scalar_commutator_unchanged',zero(W*s.Symbol('f')-s.Symbol('f')*W))

# Normal r is the original fourth Clifford direction.
gamma_n=gam[3];normal=-c0*S[3]*g5
connection=-c0*sum((S[i]*gamma_n*gam[i]*g5 for i in range(3)),z4)
naive=s.I*gamma_n*W
ck('whole_collar_split',zero(naive-connection-normal))
ck('tangent_axial_part_is_skew_normal_connection',zero(connection.H+connection))
ck('normal_axial_part_is_hermitian_boundary_potential',zero(normal.H-normal))
ck('normal_connection_preserves_normal_clifford',zero(connection*gamma_n-gamma_n*connection))
ck('boundary_potential_anticommutes_with_normal',zero(normal*gamma_n+gamma_n*normal))
ck('operator_normal_connection_is_three_times_axial_representative_connection',zero(connection-3*ak[3]))
hermitian_defect=naive-naive.H
ck('naive_boundary_generator_hermiticity_defect',zero(hermitian_defect-2*connection))
ck('naive_hermitian_iff_all_three_tangent_components_zero',s.Matrix.hstack(*[s.Matrix(connection.diff(S[i])).reshape(16,1) for i in range(3)]).rank()==3)

# Full mass retained: arbitrary constant Hermitian Phi commutes with spin.
# In a single mass eigenspace m these are symbolic identities, so they lift
# to the entire original Phi by spectral calculus without fitting masses.
p=s.symbols('p1:4',real=True);m=s.Symbol('m',real=True)
B0=sum((s.I*gamma_n*gam[i]*p[i] for i in range(3)),z4)+s.I*gamma_n*g5*m
ck('original_whole_boundary_B0_hermitian',zero(B0-B0.H))
ck('original_whole_boundary_B0_normal_odd',zero(B0*gamma_n+gamma_n*B0))
ck('original_whole_boundary_square',zero(B0*B0-(sum(x*x for x in p)+m*m)*I4))
ck('new_boundary_potential_is_hermitian',zero(B0+normal-(B0+normal).H))

# Universal finite exterior algebra: four psi and four external eta.
# All six entries of M and all 24 current entries are independent symbols.
def product(a,b):
    if a&b:return None,0
    swaps=sum((b&((1<<i)-1)).bit_count() for i in range(8) if a>>i&1)
    return a|b,(-1)**swaps
def add(P,Q):
    R=P.copy()
    for k,v in Q.items():R[k]=R.get(k,0)+v
    return {k:s.expand(v) for k,v in R.items() if v!=0}
def mul(P,Q):
    R={}
    for a,x in P.items():
        for b,y in Q.items():
            k,sgn=product(a,b)
            if sgn:R[k]=R.get(k,0)+sgn*x*y
    return {k:s.expand(v) for k,v in R.items() if v!=0}
def scale(P,c):return {k:c*v for k,v in P.items()}
def exp_nil(P):
    R={0:s.Integer(1)};T=R
    for n in range(1,5):
        T=scale(mul(T,P),s.Rational(1,n));R=add(R,T)
    return R
pairs=list(itertools.combinations(range(4),2))
M=s.symbols('M:6');r=s.symbols('r:4');ell=s.symbols('ell:4')
sv=s.symbols('s:4');AA=s.Symbol('A',positive=True)
J=[{(1<<i)|(1<<j):s.Symbol(f'j{a}_{i}{j}') for i,j in pairs} for a in range(4)]
F={(1<<i)|(1<<j):-M[n] for n,(i,j) in enumerate(pairs)}
F.update({(1<<i)|(1<<(i+4)):-1 for i in range(4)})
raw=F
for a in range(4):raw=add(raw,scale(J[a],c0*sv[a]))
fermion_first=exp_nil(raw)
fermion_first={k:v for k,v in fermion_first.items() if k&15==15}
means=[(r[a]+ell[a])/(2*AA) for a in range(4)]
def gauss_poly(poly):
    total=0
    for powers,coef in s.Poly(poly,*sv).terms():
        term=coef
        for a,n in enumerate(powers):
            assert n<=2
            term*=1 if n==0 else means[a] if n==1 else means[a]**2+1/(2*AA)
        total+=term
    return s.expand(total)
direct={k:gauss_poly(v) for k,v in fermion_first.items()}
Feff=F
for a in range(4):
    Feff=add(Feff,scale(J[a],c0*(r[a]+ell[a])/(2*AA)))
    Feff=add(Feff,scale(mul(J[a],J[a]),c0*c0/(4*AA)))
aux_first={k:v for k,v in exp_nil(Feff).items() if k&15==15}
for k in sorted(set(direct)|set(aux_first)):
    ck(f'whole_fermion_source_coefficient_two_integration_orders_{k}',s.expand(direct.get(k,0)-aux_first.get(k,0))==0)
ck('all_external_eta_sectors_retained',set(direct)=={15,63,95,159,111,175,207,255})

kappa=s.Symbol('kappa_phys',positive=True)
Aactual=3/(4*kappa*kappa)
ck('exact_auxiliary_contact_2kappa_squared_over3',s.simplify(1/(2*Aactual)-2*kappa*kappa/3)==0)
ck('exact_current_response_kappa_squared_over2',s.simplify(c0/(2*Aactual)-kappa*kappa/2)==0)
ck('exact_transverse_current_spectral_weight_kappa_fourth_over4',s.simplify(c0*c0/(4*Aactual*Aactual)-kappa**4/4)==0)

# Check reference/extra-weight force AND second-derivative contacts without
# dividing by the weight. All coefficients are independent; zeros allowed.
x,y=s.symbols('x y');u,v=s.symbols('u v');variance=1/(2*AA)
mons=[1,x,y,x*x,x*y,y*y,x**3,x*x*y,x*y*y,y**3,x*x*y*y]
fc=s.symbols('f:11');extra=sum((a*b for a,b in zip(fc,mons)),s.Integer(0))
def moment(n,mean):
    row=[s.Integer(1),mean]
    for k in range(2,n+1):row.append(s.expand(mean*row[-1]+(k-1)*variance*row[-2]))
    return row[n]
def eg(poly):
    return s.expand(sum((coef*moment(i,u)*moment(j,v) for (i,j),coef in s.Poly(poly,x,y).terms()),s.Integer(0)))
for a,(xx,mm) in enumerate([(x,u),(y,v)]):
    ck(f'extra_weight_first_source_force_{a}',s.expand(eg(xx*extra)-mm*eg(extra)-variance*eg(s.diff(extra,xx)))==0)
    for b,(yy,nn) in enumerate([(x,u),(y,v)]):
        rhs=mm*nn*eg(extra)+(variance*eg(extra) if a==b else 0)+variance*(mm*eg(s.diff(extra,yy))+nn*eg(s.diff(extra,xx)))+variance**2*eg(s.diff(extra,xx,yy))
        ck(f'extra_weight_second_source_contacts_{a}_{b}',s.expand(eg(xx*yy*extra)-rhs)==0)

out={'result_id':'TORSION20261104D','passed':len(checks),'checks':checks,
 'connection_fiber':{'real_metric_dimension':24,'dirac_image_rank':8,'selfadjoint_constraints_rank':4,'symmetric_fiber_dimension':20,'symmetric_dirac_image_rank':4,'invisible_kernel_dimension':16},
 'collar':{'A_normal':'-3/4 sum_i S_i gamma4 gammai gamma5','B_selfadjoint':'B0-3/4 S4 gamma5','naive_shift':'A_normal-3/4 S4 gamma5','frame_ode':'partial_r U=-A_normal U, U(0)=I','frame_carry':'U* B_selfadjoint U including its tangential derivatives; both endpoints and original reference/state/domain'},
 'source_identity':{'Q':'r+c0 J5','gaussian_factor':'(pi/A)^(number_of_real_S_coordinates/2) exp(sum(ell+r)^2/(4A))','all_fermion_source_coefficients_checked':len(direct),'exact_relation':'G_SS,c=(2A)^-1 I+(2A)^-2 G_QQ,c','transverse_relation':'G_STST,c=(2kappa_phys^2/3) P_T+(kappa_phys^4/4) G_JTJT,c','extra_weight_unnormalized_identity':'4A^2 Z_SS=2A delta Z+Z[Q Q F+Q Fprime+Fprime Q+Fsecond] at ell0; connected subtraction retained; no division by F','scope':'actual regulated minimal historical P25 Gaussian parent; same original remaining fields, reference and sources; any extra S-dependent factor requires its own force/contact terms'},
 'primitive_original_torsion_field_image_selected':False,'full_Sadm_derived':False,'physical_Htrpos_derived':False,'common_physical_spectrum_computed':False}
Path(__file__).with_name('verbindung-und-ganze-quelle-D.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'result_id':out['result_id'],'passed':out['passed'],'connection_fiber':out['connection_fiber'],'whole_fermion_source_coefficients':len(direct),'full_solution':False},ensure_ascii=False))
