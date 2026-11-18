"""Original two-face graph -> face currents, with exact finite controls.

Infinite statements are proved in the companion proof. Finite matrices here
are controls of that same graph, never finite CCR representations or a new
physical model.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import sympy as s

ROOT=Path(__file__).resolve().parents[3]
SOURCE=ROOT/'experiments/theory-contracts/seam-full-collar-reconstruction-20261015/checker.py'
spec=importlib.util.spec_from_file_location('original_full_collar_checker',SOURCE)
original=importlib.util.module_from_spec(spec);spec.loader.exec_module(original)
OUT=Path(__file__).resolve().parent
OUT.mkdir(exist_ok=True)

def occupation(r,alpha):
    # Original ell=4 log2, beta=2ell; both quarter spectra are unchanged.
    return 1/(1+s.Integer(2)**(8*(r+alpha)))

def central(G,A,B):
    Q=s.eye(G.rows)-G
    return s.simplify(s.trace(G*A*Q*B*G-G*B*Q*A*G))

rows={}
for alpha in (s.Rational(1,4),s.Rational(3,4)):
    N=2;rs=list(range(-N,N+1));d=len(rs)
    T=s.diag(*[s.Integer(2)**(-4*(r+alpha)) for r in rs])
    A,C,G=original.graph(T)
    n=s.eye(d)-A
    assert G*G==G and G.T==G
    assert n==s.diag(*[occupation(r,alpha) for r in rs])
    S=s.zeros(d)
    for k in range(d-1):S[k+1,k]=1
    F=(S+S.T)/2;H=(-s.I*S+s.I*S.T)/2
    Z=s.zeros(d)
    a0=s.diag(F,Z);b0=s.diag(H,Z)
    a1=s.diag(Z,F);b1=s.diag(Z,H)
    c0=central(G,a0,b0);c1=central(G,a1,b1)
    total=central(G,a0+a1,b0+b1)
    assert total==0 and c0+c1==0
    assert central(G,a0,b1)==central(G,a1,b0)==0
    endpoint=n[0,0]-n[-1,-1]
    assert c1==-s.I*endpoint/2
    # Already on the actual fixed source, a common gauge and a one-face
    # unitary move different full covariances, even with equal face1 marginals.
    U=s.eye(d);U[0,0]=s.I
    joint=s.diag(U,U);face=s.diag(s.eye(d),U)
    Gj=joint*G*joint.conjugate().T
    Gf=face*G*face.conjugate().T
    assert Gj[d:,d:]==Gf[d:,d:] and Gj!=Gf
    assert Gj[:d,d:]==U*C*U.conjugate().T
    assert Gf[:d,d:]==C*U.conjugate().T
    # General orthogonal-graph identity for a relative endpoint transport.
    V=s.eye(d).col_join(U*T)
    relative=V*(V.conjugate().T*V).inv()*V.conjugate().T
    assert relative==Gf
    controls=[]
    for window in (2,4,8):
        for mode in (1,2,3):
            pp=sum(occupation(r,alpha)*(1-occupation(r+mode,alpha))
                   for r in range(-window,window-mode+1))
            mm=sum(occupation(r+mode,alpha)*(1-occupation(r,alpha))
                   for r in range(-window,window-mode+1))
            ratio=s.Rational(1,256)**mode
            assert s.factor(mm-ratio*pp)==0
            end=sum(occupation(r,alpha)-occupation(r+mode,alpha)
                    for r in range(-window,window-mode+1))
            assert s.factor(pp-mm-end)==0
            target=s.Integer(mode)/(1-ratio)
            controls.append({'window':window,'mode':mode,'exact_KMS_ratio':str(ratio),
                             'abs_error_to_infinite_two_point':float(abs(pp-target))})
    rows[str(alpha)]={'graph_projector_reused':True,
      'finite_source_face1_cocycle':str(c1),'common_two_face_cocycle':str(total),
      'infinite_face1_limit':'-I/2','infinite_face0_limit':'+I/2',
      'different_full_lifts_same_reduced_covariance':True,
      'relative_endpoint_graph_identity':True,'KMS_controls':controls}

# Actual scalar transport block, before attaching any full Dirac/real structure.
r,theta,z=s.symbols('r theta z',real=True)
alpha=s.symbols('alpha',real=True)
chi=s.Function('chi')(r);phi=s.Function('phi')(theta);psi=s.Function('psi')(r,theta)
g=s.exp(s.I*z*chi*phi)
transport=s.simplify(g*(s.diff(psi/g,r)-s.I*s.diff(psi/g,theta)+alpha*psi/g))
expected=s.diff(psi,r)-s.I*s.diff(psi,theta)+alpha*psi-z*(s.I*s.diff(chi,r)*phi+chi*s.diff(phi,theta))*psi
assert s.simplify(transport-expected)==0
assert s.diff(expected,z,2)==0

# A closed-loop diagnostic that fixes all four original circle marks.
f=s.sin(4*theta);h=1-s.cos(4*theta)
assert all(f.subs(theta,s.pi*k/2)==h.subs(theta,s.pi*k/2)==0 for k in range(4))
assert s.simplify(f.subs(theta,theta+s.pi/2)-f)==0
assert s.simplify(h.subs(theta,theta+s.pi/2)-h)==0
fk=-s.I/2;hk=-s.Rational(1,2)
marked_cocycle=s.simplify(4*(s.conjugate(fk)*hk-fk*s.conjugate(hk)))
assert marked_cocycle==-2*s.I

# Genuine source partition functions, retaining the entire charge carry in
# the infinite formulas. Finite products/sums below only check their identity.
q=s.Rational(1,256);partitions={}
for alpha in (s.Rational(1,4),s.Rational(3,4)):
    product=s.prod((1+q**(r+alpha))*(1+q**(r+1-alpha)) for r in range(12))
    charge=sum(q**(s.Rational(k*(k-1),2)+alpha*k) for k in range(-12,13))
    oscillator=s.prod(1/(1-q**m) for m in range(1,13))
    error=float(abs(product-charge*oscillator))
    assert error<1e-25
    current_mean=sum(occupation(r,alpha)-(1 if r<0 else 0) for r in range(-12,13))
    charge_mean=sum(k*q**(s.Rational(k*(k-1),2)+alpha*k) for k in range(-12,13))/charge
    mean_error=float(abs(current_mean-charge_mean))
    assert mean_error<1e-25
    partitions[str(alpha)]={'finite_product_charge_oscillator_error':error,
                            'original_sea_zero_mode_matches_charge_sum_error':mean_error}

result={'verdict':'EXACT_CONDITIONAL_TWO_FACE_TO_CURRENT_BRIDGE_WITH_DIFFERENT_SOURCE_LIFTS',
 'source_graph_constructor':str(SOURCE.relative_to(ROOT)),
 'source_graph_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
 'statistical_premise':'Given charged CAR and pure quasifree graph state; not derived from P1.',
 'analytic_scope':'All smooth real scalar source profiles, all Fourier modes, all ell>0. Infinite formulas require the companion proof; finite checks are controls only.',
 'cases':rows,'charge_partition_controls':partitions,
 'relative_endpoint_transport_block':'g_chi (partial_r+B_alpha) g_chi^* = partial_r+B_alpha-i chi_prime phi-chi phi_prime',
 'four_mark_loop':{'profiles':['sin(4 theta)','1-cos(4 theta)'],
                   'both_zero_at_all_four_marks':True,'both_quarter_periodic':True,
                   'face1_cocycle':str(marked_cocycle),'face1_Weyl_rectangle_phase':'exp(2 I s t)',
                   'common_two_face_rectangle_phase':'1'},
 'original_time_scale':'K_alpha from original source; [K_alpha,J_-m]=m J_-m/8. Normal length is not physical time.',
 'not_derived':['CAR/statistics from primitive P1','physical execution of a one-face intervention','original global P1 to E8 field identification','event menu/rates, pointer and record selection','3+1D chirality, coupling selection and gravity'],
 'complete_solution':False}
(OUT/'results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
print(json.dumps(result,ensure_ascii=False,indent=2))
