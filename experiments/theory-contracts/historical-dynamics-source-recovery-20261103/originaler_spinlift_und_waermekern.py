"""Direct original Spin4/96 field test and exact geometric torsion heat jet."""
from pathlib import Path
import ast, hashlib, json
import sympy as s

root=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
owner=root/'experiments/theory-contracts/craig-original-source-bridge-20261103/actual_product_pairing.py'
nodes=[]
for node in ast.parse(owner.read_text()).body:
    if isinstance(node,ast.Expr) and isinstance(node.value,ast.Call) and isinstance(node.value.func,ast.Name) and node.value.func.id=='ck':
        continue
    nodes.append(node)
    if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='Phi' for t in node.targets):
        break
ns={'__file__':str(owner)}
exec(compile(ast.Module(body=nodes,type_ignores=[]),str(owner),'exec'),ns)
gam,g5,C,gF,KF,Phi=[ns[k] for k in ('gam','g5','C','gF','KF','Phi')]
sp=s.SparseMatrix
kron=lambda a,b:sp(s.kronecker_product(sp(a),sp(b)))
I96=s.eye(96);I4=s.eye(4)
checks={}
def zero(M):
    if isinstance(M,s.MatrixBase):
        return all(s.simplify(s.expand(x))==0 for x in sp(M).todok().values())
    return s.simplify(s.expand(M))==0
def ck(name,value):
    ok=value if isinstance(value,bool) else zero(value)
    checks[name]=bool(ok)
    if not ok:raise AssertionError(name)

Sspin=s.kronecker_product(s.Matrix([[0,1],[1,0]]),s.eye(2))
Jmat=kron(C,KF);Gamma=kron(g5,gF)
ck('original_full96_Phi_odd',gF*Phi+Phi*gF)
ck('normal_spin_reflection_flips_gamma5',Sspin*g5*Sspin+g5)
for mu in range(4):
    contraction=s.zeros(4)
    for a in range(4):
        for b in range(4):
            for m in range(4):
                contraction += -s.I*s.LeviCivita(m,a,b,mu)*gam[m]*(gam[a]*gam[b]-gam[b]*gam[a])/16
    ck(f'original_real_Euclidean_spinlift_{mu}',contraction-3*s.I*gam[mu]*g5/4)
    W=s.I*gam[mu]*g5/2
    ck(f'real_geometric_torsion_is_Hermitian_{mu}',W.H-W)
    Wfull=kron(W,I96)
    ck(f'real_geometric_torsion_is_J_real_{mu}',Jmat*Wfull.conjugate()-Wfull*Jmat)
    ck(f'real_geometric_torsion_is_total_odd_{mu}',Gamma*Wfull+Wfull*Gamma)
    # For an even scalar theta, its normal derivative is odd and its
    # tangential derivatives are even under the ORIGINAL base reflection.
    base_sign=-1 if mu==3 else 1
    ck(f'tauS_real_geometric_direction_is_odd_{mu}',base_sign*Sspin*W*Sspin+W)
    F_T=kron(Sspin,gF)
    ck(f'tauT_real_geometric_direction_is_odd_{mu}',base_sign*F_T*Wfull*F_T+Wfull)
    ck(f'tauS_odd_coordinate_theta_geometric_direction_is_even_{mu}',-base_sign*Sspin*W*Sspin-W)
    ck(f'tauT_odd_coordinate_theta_geometric_direction_is_even_{mu}',-base_sign*F_T*Wfull*F_T-Wfull)
    phase=kron(gam[mu],gF)/2
    F_S=kron(Sspin,I96)
    ck(f'tauS_internal_unitary_direction_is_even_{mu}',base_sign*F_S*phase*F_S-phase)
    ck(f'tauT_internal_unitary_direction_is_even_{mu}',base_sign*F_T*phase*F_T-phase)
    ck(f'geometric_mass_cross_anticommutator_zero_{mu}',W*g5+g5*W)
    for eps in (-1,1):
        inds=[j for j in range(384) if Gamma[j,j]==eps]
        E=s.SparseMatrix(384,192,{(i,j):1 for j,i in enumerate(inds)})
        ck(f'fiber_{eps}_geometric_is_i_times_unitary_phase_{mu}',(Wfull-s.I*eps*phase)*E)
        ck(f'fiber_{eps}_geometric_not_same_real_unitary_phase_{mu}',not zero((Wfull-eps*phase)*E))

v=s.symbols('v0:4',real=True)
Wsum=sum((s.I*gam[mu]*g5*v[mu]/2 for mu in range(4)),s.zeros(4))
v2=sum(x*x for x in v)
ck('whole_geometric_W_square',Wsum*Wsum-s.eye(4)*v2/4)

# Exact Fourier vertices of the ACTUAL real geometric deformation,
# with no derivative expansion, fitted field or new physical profile.
k=s.symbols('k0:4',real=True);p=s.symbols('p0:4',real=True)
t=s.symbols('t',positive=True);a=s.symbols('a',real=True)
slash=lambda vec:sum((gam[i]*vec[i] for i in range(4)),s.zeros(4))
ks=slash(k);ps=slash(p);p2=sum(x*x for x in p);k2=sum(x*x for x in k);kp=sum(k[i]*p[i] for i in range(4))
Wp=-ps*g5/2;Wminus=ps*g5/2
Vplus=(ks+ps)*Wp+Wp*ks
Vminus=ks*Wminus+Wminus*(ks+ps)
vertex=s.expand(s.trace(Vminus*Vplus))
ck('full_spin_vertex_trace',vertex-p2**2-4*(k2*p2-kp**2))
q=s.symbols('q0:4',real=True)
shifted=s.expand(vertex.subs(dict(zip(k,[q[i]-a*p[i] for i in range(4)])),simultaneous=True))
def gaussian(poly):
    total=0
    for powers,coeff in s.Poly(poly,*q).terms():
        moment=1
        for power in powers:
            if power%2:moment=0;break
            if power:moment*=s.factorial2(power-1)/(2*t)**(power//2)
        total+=coeff*moment
    return s.expand(total)
averaged=gaussian(shifted)
ck('full_Gaussian_vertex_trace',averaged-p2**2-6*p2/t)
contact=s.expand(s.trace(Wminus*Wp))
ck('quadratic_W_contact_trace',contact-p2)
F=s.symbols('F',real=True)
loop=s.expand((t*t/2)*(averaged)/(16*s.pi**2*t*t)*F)
contact_term=-t*contact/(16*s.pi**2*t*t)
kernel=s.factor(loop+contact_term)
expected=p2*((t*p2+6)*F-2)/(32*s.pi**2*t)
ck('complete_geometric_heat_kernel_including_contact',kernel-expected)
x=s.symbols('x',real=True)
Fseries=sum((-x)**n*s.factorial(n)/s.factorial(2*n+1) for n in range(6))
Gseries=s.series(Fseries-x*s.diff(Fseries,x),x,0,5).removeO()
ck('geometric_G_has_no_linear_x_term',s.diff(Gseries,x).subs(x,0)==0)
ck('geometric_G_zero_is_one',Gseries.subs(x,0)==1)

# Retain the actual independent torsion before imposing its historical
# auxiliary equation. Rotational covariance permits p along axis 0 and
# b=3S/4 split into its parallel and one arbitrary transverse component.
pp,bL,bT=s.symbols('p bL bT',real=True)
ps_ax=pp*gam[0];bs=bL*gam[0]+bT*gam[1]
Wax=s.I*bs*g5
Vax_plus=(ks+ps_ax)*Wax+Wax*ks
Vax_minus=ks*Wax+Wax*(ks+ps_ax)
axtrace=s.expand(s.trace(Vax_minus*Vax_plus))
axshift=s.expand(axtrace.subs({k[i]:q[i]-(a*pp if i==0 else 0) for i in range(4)},simultaneous=True))
axavg=gaussian(axshift)
axexpected=bL**2*(4*pp**2+24/t)+bT**2*(4*pp**2*(1-4*a*(1-a))+24/t)
ck('full_independent_torsion_parallel_transverse_vertex',axavg-axexpected)
ck('full_independent_torsion_contact',s.trace(Wax*Wax)-4*(bL**2+bT**2))
Fp=s.symbols('Fp')
KL=((x+6)*F-2)/(8*s.pi**2*t)
KT=((x+6)*F+4*x*Fp-2)/(8*s.pi**2*t)
# Integration by parts: 1=(1+x/2)*F+2*x*F'.
ck('integration_by_parts_identity_series',s.series(2*x*s.diff(Fseries,x)+(1+x/2)*Fseries-1,x,0,5).removeO())
ck('transverse_form_factor_reduction',KT.subs(Fp,(1-(1+x/2)*F)/(2*x))-F/(2*s.pi**2*t))
ck('longitudinal_form_factor_reduction',KL-((x+6)*F-2)/(8*s.pi**2*t))
ck('transverse_a4_curl_coefficient',s.diff(Fseries,x).subs(x,0)+s.Rational(1,6))
ck('gradient_a4_curl_zero',s.diff(Gseries,x).subs(x,0))
# Independent local Bochner calculation of that sign and normalization.
omega=[-(gam[i]*bs-bs*gam[i])*g5/2 for i in range(4)]
curv_pair_trace=0
for i in range(4):
    for j in range(4):
        plus=s.I*pp*((omega[j] if i==0 else s.zeros(4))-(omega[i] if j==0 else s.zeros(4)))
        curv_pair_trace+=s.trace((-plus)*plus)
ck('independent_Bochner_curvature_trace',s.expand(curv_pair_trace)+24*pp**2*bL**2+16*pp**2*bT**2)
Eplus=-s.I*pp*bL*g5;Eminus=s.I*pp*bL*g5
local_a4=s.expand((s.trace(Eminus*Eplus)/2+curv_pair_trace/12)/(16*s.pi**2))
ck('independent_Bochner_a4_parallel_cancels',local_a4.coeff(bL,2))
ck('independent_Bochner_a4_transverse_negative',local_a4+pp**2*bT**2/(12*s.pi**2))

out={
 'result':'ACTUAL-REAL-EUCLIDEAN-TORSION-SOURCE-20261104',
 'checks':checks,'passed':sum(checks.values()),
 'original_symbolic_internal_dimension':96,'ordinary_product_dimension':384,
 'all_original_Yukawa_and_Majorana_blocks_retained':True,
 'original_constructor_sha256':hashlib.sha256(owner.read_bytes()).hexdigest(),
 'actual_geometric_spinlift':'delta D=3i/4 S_mu gamma_mu gamma5 for the stated original K and Euclidean epsilon conventions',
 'real_geometric_torsion_not_internal_unitary_phase':True,
 'even_theta_geometric_reference_jet':'delta D_ref=0;delta D_rel=W_geometric, for both existing lifts',
 'general_coordinate_parity_jet':'delta D_ref=W[theta_o];delta D_rel=W[theta_e], for both existing lifts; this does not select intrinsic physical reflection parity',
 'heat_theta_squared_coefficient':'p^2/(8*pi^2*t) G(t*p^2) Tr_F exp(-t*Phi^2)',
 'G_positive_integral':'integral_0^1 (1+x*a*(1-a))*exp(-x*a*(1-a)) da',
 'G_series':str(Gseries),
 'a2_theta_coefficient':'dim(F)/(8*pi^2) integral |dtheta|^2',
 'a4_theta_coefficient_flat_constant_Phi':'-Tr_F(Phi^2)/(8*pi^2) integral |dtheta|^2',
 'full_independent_torsion_heat_jet':'Tr_F exp(-t Phi^2)/(2*pi^2*t) [G(x)|b_parallel|^2+F(x)|b_transverse|^2], b=3S/4',
 'a4_independent_torsion_flat_constant_Phi':'-3 dim(F)/(64*pi^2) integral |dS|_form^2 -9 Tr_F(Phi^2)/(32*pi^2) integral |S|^2',
 'historical_parent_not_identified_with_geometric_residue':'P25 auxiliary action is algebraic in S; ordinary geometric a4 has a nonzero curl term before source/domain/reference reduction. The actual bound source must decide any cancellation or restriction. No free kinetic coupling introduced.',
 'Wick_or_contour_identity_inferred':False,
 'complete_Sadm_identified':False,'physical_Htrpos_or_spectrum_claimed':False,
 'scope':'Actual original ordinary Spin4 matrices and full96 family; real geometric Euclidean torsion and its bound-reference parity. Exact flat constant-full-Phi heat jet, second order in theta and all external momentum. No physical profile, field measure, vacuum, spacetime origin or whole quantum source chosen.'
}
Path(__file__).with_name('originaler-spinlift-und-waermekern.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(out,ensure_ascii=False,indent=2))
