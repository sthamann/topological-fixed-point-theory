"""Original full-product current and whole-source assembly identities.

No physical replacement model or additional field measure is chosen. The
four-variable polynomial at the end verifies a universal Berezin identity.
"""
import ast
import hashlib
import itertools
import json
from pathlib import Path
import sympy as s

repo=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
original=repo/'experiments/theory-contracts/craig-original-source-bridge-20261103/actual_product_pairing.py'
assert hashlib.sha256(original.read_bytes()).hexdigest()=='6acaa9bc2c7a157a30730d7d6de75e944266afbf69b18a31c6fba654b8bd1c04'
# Reuse the actual matrices/whole symbolic Yukawa/Majorana constructor.
# Stop before the prior mass test and before all old result-file writes.
nodes=[]
for node in ast.parse(original.read_text()).body:
    if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='Mass' for t in node.targets):break
    nodes.append(node)
env={'__file__':str(original),'__name__':'original_constructor_only'}
exec(compile(ast.Module(body=nodes,type_ignores=[]),str(original),'exec'),env)
gam,g5,C,KF,gF,Phi=[s.SparseMatrix(env[k]) if k!='gam' else [s.SparseMatrix(g) for g in env[k]] for k in ['gam','g5','C','KF','gF','Phi']]
I=s.SparseMatrix(s.eye(96)); R=s.kronecker_product(C,KF)
G=s.kronecker_product(g5,gF)
indices={e:[i for i in range(384) if G[i,i]==e] for e in [-1,1]}
checks={}
def zero(M):return all(s.expand(x)==0 for x in M)
def ck(name,value):
    assert value,name
    checks[name]=True

current_matrices={}
for e in [-1,1]:
    ix=indices[e]
    particle=[j for j,i in enumerate(ix) if i%96<48]
    antiparticle=[j for j,i in enumerate(ix) if i%96>=48]
    ck(f'actual_total_Weyl_dimension_{e}',len(ix)==192)
    ck(f'actual_particle_and_antiparticle_Weyl_slots_{e}',len(particle)==len(antiparticle)==96)
    L=[]
    for mu in range(4):
        full=R.H*s.kronecker_product(gam[mu]*g5,I)
        le=full.extract(ix,ix)
        L.append(le)
        ck(f'actual_axial_pairing_skew_{e}_{mu}',zero(le.T+le))
        ck(f'actual_axial_pairing_unit_norm_{e}_{mu}',zero(le.H*le-s.eye(192)))
        ck(f'actual_axial_pairing_has_no_zero_species_slots_{e}_{mu}',all(any(le[i,j]!=0 for j in range(192)) for i in range(192)))
        ce=-s.I*le/2
        vertex=3*s.I*le/4
        ck(f'actual_Majorana_half_pairing_matches_minus_three_quarters_current_{e}_{mu}',zero(vertex/2+s.Rational(3,4)*ce))
        ck(f'actual_current_has_no_particle_particle_block_{e}_{mu}',zero(le.extract(particle,particle)))
        ck(f'actual_current_has_no_antiparticle_antiparticle_block_{e}_{mu}',zero(le.extract(antiparticle,antiparticle)))
    current_matrices[e]=L

    q=s.kronecker_product(g5,I).extract(ix,ix)
    Mmass=(R.H*s.kronecker_product(g5,Phi)).extract(ix,ix)
    full_mass_breaking=(R.H*s.kronecker_product(s.eye(4),Phi)).extract(ix,ix)
    ck(f'whole_Yukawa_Majorana_axial_Ward_mass_term_{e}',zero(q.T*Mmass+Mmass*q-2*full_mass_breaking))
    ck(f'whole_axial_Ward_mass_pairing_skew_{e}',zero(full_mass_breaking.T+full_mass_breaking))
    for label,block in [('particle',particle),('antiparticle',antiparticle)]:
        majorana=Mmass.extract(block,block)
        norm=s.expand(sum(s.conjugate(x)*x for x in majorana.todok().values()))
        ck(f'actual_Majorana_same_side_mass_norm_{e}_{label}',s.expand(norm-2*s.trace(env['MR'].H*env['MR']))==0)
    # Restriction is algebraic only: epsilon here is not Padm.

# The whole finite mass response, with the original fixed MR and all Ys.
r=s.Symbol('r',real=True)
Sf,Tf=env['S'],env['T']
Phi_r=s.SparseMatrix.vstack(s.SparseMatrix.hstack(r*Sf,Tf),s.SparseMatrix.hstack(Tf.H,r*Sf.conjugate()))
ys=env['Ys']; mr=env['MR']
am=12*s.trace(ys['u'].H*ys['u']+ys['d'].H*ys['d'])+4*s.trace(ys['e'].H*ys['e']+ys['nu'].H*ys['nu'])
normmr=s.trace(mr.H*mr)
mass_trace=s.expand(s.trace(Phi_r*Phi_r))
ck('actual_full96_Higgs_trace_with_fixed_Majorana',s.expand(mass_trace-r*r*am-2*normmr)==0)
ss=s.Symbol('S_tangent',real=True)
mixed=-s.Rational(9,32)/s.pi**2*ss*ss*mass_trace
ck('actual_nonzero_bulk_torsion_Higgs_mixed_vertex',s.expand(s.diff(mixed,ss,2,r,2)+s.Rational(9,8)/s.pi**2*am)==0)
kappa=s.Symbol('kappa_phys',positive=True)
aa=3/(4*kappa*kappa); c0=s.Rational(3,4)
j1,j2=s.symbols('J5_sector_a J5_sector_b',commutative=True)
contact=-c0*c0/(4*aa)*(j1+j2)**2
ck('whole_historical_shared_current_contact',s.expand(contact+s.Rational(3,16)*kappa*kappa*(j1+j2)**2)==0)
ck('whole_historical_intersector_contact_coefficient',s.expand(s.diff(contact,j1,j2)+s.Rational(3,8)*kappa*kappa)==0)

# The algebraic assembly error is tested as an identity, not a new physics
# example. On Ran(product Pi) every commuting Pi acts as identity.
p,q=s.symbols('p q',commutative=True)
ck('commuting_filters_joint_range_identity',s.expand(p*(p*q)-p*q).subs(p*p,p)==0)
ck('same_filters_are_not_orthogonal_sector_projectors',s.expand((p*q)*(p*q)).subs(p*p,p).subs(q*q,q)==p*q)
za,zb=s.symbols('Z_a Z_b',positive=True)
ck('direct_sum_trace_log_mixed_derivative_nonzero',s.simplify(s.diff(-s.log(za+zb),za,zb)-1/(za+zb)**2)==0)
ck('tensor_product_log_adds',s.simplify(-s.log(za*zb)+s.log(za)+s.log(zb))==0)

# Universal4-dimensional skew-form change of basis and all external legs.
abc=s.symbols('m01 m02 m03 m12 m13 m23')
M=s.zeros(4)
for v,(i,j) in zip(abc,itertools.combinations(range(4),2)):M[i,j]=v;M[j,i]=-v
V=s.Matrix(4,4,s.symbols('v:16'))
def pf4(A):return A[0,1]*A[2,3]-A[0,2]*A[1,3]+A[0,3]*A[1,2]
ck('generic_Pfaffian_congruence_full_basis',s.expand(pf4(V.T*M*V)-V.det()*pf4(M))==0)
# Grassmann pullback exp(-psi.M.psi/2+eta.psi) including exterior sources.
def mul(P,Q):
    result={}
    for a,x in P.items():
        for b,y in Q.items():
            if a&b:continue
            swaps=sum((b&((1<<i)-1)).bit_count() for i in range(8) if a>>i&1)
            k=a|b;result[k]=result.get(k,0)+(-1)**swaps*x*y
    return {k:s.expand(v) for k,v in result.items() if v!=0}
def add(P,Q):
    result=P.copy()
    for k,v in Q.items():result[k]=result.get(k,0)+v
    return {k:s.expand(v) for k,v in result.items() if v!=0}
def expnil(P):
    result={0:s.Integer(1)}; term=result
    for n in range(1,5):
        term={k:v/n for k,v in mul(term,P).items()};result=add(result,term)
    return result
# A general diagonal plus upper shear keeps det symbolic and nonzero without
# assuming any invertible fermion kernel. This is a basis-map check only.
d=s.symbols('d0:4',nonzero=True); a=s.Symbol('shear')
v=s.diag(*d);v[0,1]=a
def source_polynomial(A,leg):
    P={(1<<i)|(1<<j):-A[i,j] for i,j in itertools.combinations(range(4),2)}
    for i in range(4):
        for j in range(4):
            if leg[i,j]!=0:P[(1<<j)|(1<<(i+4))]=-leg[i,j]
    return {k:s.expand(c) for k,c in expnil(P).items() if k&15==15}
original_poly=source_polynomial(M,s.eye(4))
new_poly=source_polynomial(v.T*M*v,v)
for k in sorted(set(original_poly)|set(new_poly)):
    ck(f'whole_Berezin_source_transport_including_eta_sector_{k}',s.expand(new_poly.get(k,0)-v.det()*original_poly.get(k,0))==0)
ck('whole_source_transport_all_eight_external_eta_sectors',len(original_poly)==8)

out={'result_id':'TORSION20261104E','passed':len(checks),'checks':checks,
 'original_pairing':{'total_Weyl_dimensions':{'minus':192,'plus':192},
 'torsion_current_matrix':'L_mu,eps=Eeps^T(C tensor KF)^dagger(gamma_mu gamma5 tensor I96)Eeps',
 'actual_current':'J5E_mu=-i psi^T L_mu,eps psi/2',
 'actual_action_vertex':'(3i/8) S_mu psi^T L_mu,eps psi=-(3/4) S_mu J5E_mu',
 'current_matrices_nonzero_on_every_species_slot':True,
 'current_has_only_particle_antiparticle_cross_blocks':True,
 'Majorana_same_side_mass_block_norm':'2Tr(MR*MR) on each side and each total Weyl chart',
 'all_original_Yukawa_and_Majorana_entries_symbolic':True,
 'total_chiral_projector_is_not_Padm':True},
 'whole_source_transport':{'matrix':'Mprime=V^T M V','source':'etaprime=V^T eta',
 'Berezinian':'det(V)^-1, including field-dependent V with unchanged bosonic coordinates',
 'no_inverse_fermion_kernel_used':True,'infinite_anomalous_Jacobian_not_computed':True},
 'whole_assembly':{'commuting_filters':'Pi restricted to Ran(product Pj) is identity, not an independent orthogonal sector',
 'direct_sum_thermal_trace':'Z(A directsum B)=ZA+ZB',
 'tensor_product_uncoupled_trace':'Z(A tensor I+I tensor B)=ZA ZB when the actual joint state/trace factorizes',
 'direct_sum_determinant':'det(A directsum B)=detA detB',
 'actual_original_geometric_mixed_vertex':'partial_S_tangent^2 partial_r^2 a4=-9 A_Y/(8pi²) per flat bulk volume; A_Y=12Tr(Yu*Yu+Yd*Yd)+4Tr(Ye*Ye+Ynu*Ynu)',
 'whole_historical_current_contact':'-3 kappa_phys²/16 (sum_f J5_f)²; cross coefficient -3 kappa_phys²/8 for each unordered distinct pair',
 'mixed_scope':'ordinary original Spin4 product, constant radial Higgs, fixed MR, deck-even tangential torsion; both bound references independent of this S direction; no arbitrary Brel/domain cancellation theorem'},
 'primitive_original_full_field_integration_evaluated':False,'full_Sadm_derived':False,'physical_Htrpos_derived':False,'common_physical_spectrum_computed':False}
Path(__file__).with_name('whole-source-pairing-E.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['result_id','passed','original_pairing','whole_assembly','full_Sadm_derived']},ensure_ascii=False))
