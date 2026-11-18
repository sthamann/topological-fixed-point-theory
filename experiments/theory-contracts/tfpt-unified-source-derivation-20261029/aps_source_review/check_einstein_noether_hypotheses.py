"""R081: decisive tests of four origin hypotheses on the existing source.

Native seed252 is a reproducible constructor control, not fitted physical Ys.
This does not construct the primitive measure, spacetime or a graviton.
"""
from pathlib import Path
import sys, json
import numpy as np
from scipy.linalg import expm
from scipy.special import expit

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'verification'))
import v252_full_finite_triple as original

checks = {}
metrics = {}
def need(name, condition):
    checks[name] = bool(condition)
    if not condition:
        raise RuntimeError(name)
def small(name, residual, tolerance=1e-8):
    metrics[name] = float(abs(residual))
    need(name, abs(residual) < tolerance)
def matrix_function(a, f):
    e,u = np.linalg.eigh(a)
    return (u*f(e)) @ u.conj().T
def divided_difference(e, f, df):
    delta=e[:,None]-e[None,:]
    numerator=f(e)[:,None]-f(e)[None,:]
    same=abs(delta)<1e-9
    out=np.zeros_like(delta)
    np.divide(numerator,delta,out=out,where=~same)
    out[same]=np.broadcast_to(df(e)[:,None],delta.shape)[same]
    return out
def frechet(a, b, f, df):
    e,u=np.linalg.eigh(a)
    return u@(divided_difference(e,f,df)*(u.conj().T@b@u))@u.conj().T

rng=np.random.default_rng(252)
ys=[rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)) for _ in range(4)]
mr=rng.normal(size=(3,3))
D=original.full_D(*ys,mr+mr.T)
gamma=original.grading(); swap=original.Kswap()
quaternions=[np.array([[0,1],[-1,0]],complex),
             np.array([[0,1j],[1j,0]],complex),np.eye(2),np.diag([-1j,1j])]
bp=original.rep(1,np.zeros((2,2)),np.zeros((3,3)))
V=[]
for q in quaternions:
    a=original.rep(0,q,np.zeros((3,3)))
    om=a@(D@bp-bp@D); sa=om+om.conj().T
    V.append(sa+swap@sa.conj()@swap)
ix=[i for i,f in enumerate(original.P) if f[0] in ('u','d')]
sub=lambda a:a[np.ix_(ix,ix)]
g=sub(gamma); T=[-1j*g@sub(v) for v in V]
need('actual_nondoubled_quark_sector_36',len(ix)==36)
need('original_force_tuple_hermitian_and_odd',all(
    np.linalg.norm(t-t.conj().T)<1e-10 and np.linalg.norm(t@g+g@t)<1e-10 for t in T))

# A. One given mode. Full CAR log partition, no classical force substitution.
ell=.2; p=.25; radius=.3; point=np.array([0.,0.,radius,0.])
B=[2*ell*t for t in T]
def generator(x):return 2*ell*p*g+sum((x[a]*B[a] for a in range(4)),np.zeros((36,36),complex))
def W(x):return float(np.logaddexp(0,np.linalg.eigvalsh(generator(x))).sum())
K=generator(point); C=matrix_function(K,expit)
derivs=[frechet(K,b,expit,lambda z:expit(z)*(1-expit(z))) for b in B]
gradient=np.array([np.trace(C@b).real for b in B])
hessian=np.array([[np.trace(derivs[a]@B[b]).real for b in range(4)] for a in range(4)])
ordinary=np.array([[np.trace(C@B[a]@(np.eye(36)-C)@B[b]).real for b in range(4)] for a in range(4)])
small('noether_tangent_first_ward',gradient[0])
small('noether_orbit_second_ward',radius**2*hessian[0,0]-radius*gradient[2])
small('gibbs_BKM_hessian_is_symmetric',np.linalg.norm(hessian-hessian.T))
need('gibbs_BKM_hessian_is_positive',np.linalg.eigvalsh(hessian).min()>0)
need('instantaneous_variance_is_not_BKM_hessian',np.linalg.norm(ordinary-hessian)>1e-3)
need('naive_variance_breaks_second_orbit_ward',abs(radius**2*ordinary[0,0]-radius*gradient[2])>1e-4)
small('whole_partition_is_constant_on_actual_higgs_rotation',W([radius*np.sin(.37),0,radius*np.cos(.37),0])-W(point))
eps=2e-4; unit=np.eye(4)[0]
small('BKM_hessian_matches_independent_partition_difference',
      (W(point+eps*unit)-2*W(point)+W(point-eps*unit))/eps**2-hessian[0,0],3e-6)
metrics['BKM_tangent']=float(hessian[0,0]);metrics['instantaneous_tangent_variance']=float(ordinary[0,0])
need('spectral_orbit_invariance_does_not_force_radial_stationarity',abs(gradient[2])>.01)
X=T[0]; orbit_tangent=1j*(X@K-K@X)
small('unitary_orbit_partition_derivative_is_identity',np.trace(C@orbit_tangent))
need('odd_force_as_unitary_generator_changes_principal_symbol',np.linalg.norm(X@g-g@X)>1)
metrics['radial_W_gradient']=float(gradient[2])
metrics['odd_force_principal_symbol_change_norm']=float(np.linalg.norm(X@g-g@X))

# B. Entanglement first law is exact at the chosen reference, not an equilibrium selector.
I=np.eye(36)
def entropy(c):
    e=np.linalg.eigvalsh(c)
    return float(np.sum(-e*np.log(e)-(1-e)*np.log1p(-e)))
dc=derivs[2]
predicted=-np.trace(K@dc).real
cp=matrix_function(generator(point+eps*np.eye(4)[2]),expit)
cm=matrix_function(generator(point-eps*np.eye(4)[2]),expit)
small('native_entropy_first_law', (entropy(cp)-entropy(cm))/(2*eps)-predicted,2e-6)
def relative_entropy(c,c0):
    lc=matrix_function(c,np.log);lc0=matrix_function(c0,np.log)
    lq=matrix_function(I-c,np.log);lq0=matrix_function(I-c0,np.log)
    return float(np.trace(c@(lc-lc0)+(I-c)@(lq-lq0)).real)
changed=matrix_function(generator(point+.05*np.eye(4)[2]),expit)
need('first_law_allows_distinct_positive_relative_entropy_state',relative_entropy(changed,C)>1e-5)
metrics['relative_entropy_to_distinct_radial_state']=relative_entropy(changed,C)

# C. Genuine native chiral marginal. P is diagnostic original grading+, not asserted P_F/Padm.
keep=np.where(np.diag(g).real>0)[0]
take=lambda a:a[np.ix_(keep,keep)]
logit=lambda z:np.log(z)-np.log1p(-z)
dlogit=lambda z:1/(z*(1-z))
CP=take(C); Keff=matrix_function(CP,logit); Kbare=take(K)
need('native_chiral_keep_space_18',len(keep)==18)
need('bare_higgs_vertex_compresses_to_zero',np.linalg.norm(take(B[2]))<1e-10)
need('native_marginal_differs_from_bare_compressed_gibbs',np.linalg.norm(CP-matrix_function(Kbare,expit))>1e-4)
dKeff=frechet(CP,take(dc),logit,dlogit)
need('integrating_original_partner_retains_higgs_response',np.linalg.norm(dKeff)>1e-3)
small('retained_vertex_reconstructs_full_marginal_derivative',np.linalg.norm(
    frechet(Keff,dKeff,expit,lambda z:expit(z)*(1-expit(z)))-take(dc)))
effective=lambda x:matrix_function(take(matrix_function(generator(x),expit)),logit)
small('effective_vertex_matches_independent_logit_difference',np.linalg.norm(
    (effective(point+eps*np.eye(4)[2])-effective(point-eps*np.eye(4)[2]))/(2*eps)-dKeff),2e-6)
metrics['native_marginal_covariance_defect']=float(np.linalg.norm(CP-matrix_function(Kbare,expit)))
metrics['native_effective_higgs_vertex_norm']=float(np.linalg.norm(dKeff))

# Independent 4-state Fock trace using an actual singular-pair subblock of this same source.
svals=np.linalg.svd(ys[0],compute_uv=False); s=float(svals[0])
kpair=2*ell*np.array([[p,radius*s],[radius*s,-p]])
a0=np.array([[0,1],[0,0]],complex); z=np.diag([1,-1])
ann=[np.kron(a0,np.eye(2)),np.kron(z,a0)]
dk=sum((kpair[i,j]*(ann[i].conj().T@ann[j]) for i in range(2) for j in range(2)),np.zeros((4,4),complex))
R=expm(dk); Z=float(np.trace(R).real)
marginal=R.reshape(2,2,2,2).trace(axis1=1,axis2=3)
cpair=matrix_function(kpair,expit); ce=float(cpair[0,0].real); ke=float(logit(ce))
Zkeep=1+np.exp(ke); scalar=Z/Zkeep
small('independent_fock_full_partition',Z-np.linalg.det(np.eye(2)+expm(kpair)).real)
small('independent_fock_full_unnormalized_marginal',np.linalg.norm(marginal-scalar*np.diag([1,np.exp(ke)])))
small('independent_fock_all_retained_words_via_density',np.linalg.norm(marginal/Z-np.diag([1-ce,ce])))
need('dropping_partition_ratio_changes_absolute_source',abs(scalar-1)>1e-2)
metrics['one_actual_singular_pair_absolute_trace_ratio']=scalar
# A requested independent review proposed sigmaX as a false counterexample.
false_counter=np.array([[0.,1.],[1.,0.]])
small('review_sigmaX_actually_has_half_diagonal',matrix_function(false_counter,expit)[0,0]-.5)
need('commutation_is_sufficient_not_necessary_for_compressed_logit_equality',
     abs(logit(matrix_function(false_counter,expit)[0,0].real))<1e-12
     and np.linalg.norm(false_counter@np.diag([1,0])-np.diag([1,0])@false_counter)>1)
E=np.sqrt(p*p+radius*radius*s*s)
small('actual_singular_pair_matches_exact_retained_occupation',ce-(1+(p/E)*np.tanh(ell*E))/2)
time=.7
physical_amplitude=expm(1j*time*kpair/(2*ell))[0,0]
need('exact_static_marginal_is_not_autonomous_original_physical_time',abs(physical_amplitude)<.999)
small('original_retained_time_amplitude_matches_full_partner_memory',
      physical_amplitude-(np.cos(E*time)+1j*p/E*np.sin(E*time)))
metrics['retained_physical_time_amplitude_modulus']=float(abs(physical_amplitude))

# D. Same original circle principal symbol and full finite internal fibre.
Phi=-1j*gamma@D
small('native_full96_Phi_is_hermitian',np.linalg.norm(Phi-Phi.conj().T))
small('native_full96_Phi_is_odd',np.linalg.norm(gamma@Phi+Phi@gamma))
mass2=np.linalg.eigvalsh(Phi@Phi)
times=[1e-3,2.5e-4,6.25e-5,1.5625e-5]
heats=[]; dimensions=[]
for t in times:
    ps=np.arange(-3000,3001)+.25
    theta=np.exp(-t*ps**2).sum()
    need('tail_control_at_t_'+str(t),np.exp(-t*3000**2)<1e-55)
    heat=float(theta*np.exp(-t*mass2).sum());heats.append(heat)
    poisson=np.sqrt(np.pi/t)*sum(np.exp(-np.pi**2*k*k/t)*np.cos(2*np.pi*k*.25) for k in range(-3,4))
    small('poisson_circle_heat_t_'+str(t),theta/poisson-1,2e-12)
    d=1+2*t*float(np.sum(mass2*np.exp(-t*mass2))/np.sum(np.exp(-t*mass2)))
    dimensions.append(d)
need('finite_internal_fibre_leaves_dimension_one',abs(dimensions[-1]-1)<.002)
need('finite_internal_multiplicity_does_not_give_dimension_three',abs(dimensions[-1]-3)>1.9)
metrics['heat_times']=times;metrics['spatial_spectral_dimensions']=dimensions
metrics['hypothesis_scope']='Given circle + fixed finite internal fibre only; no exclusion of other emergent or holographic spacetime maps.'

out={'date':'2026-11-02','result_id':'R081','checks':checks,'metrics':metrics,
     'passed':len(checks),'failed':0,'ledger_changed':False,'new_lean_run':False,
     'primitive_field_measure_derived':False,'physical_3plus1_or_spin2_derived':False}
Path(__file__).with_name('EINSTEIN-NOETHER-HYPOTHESES.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
