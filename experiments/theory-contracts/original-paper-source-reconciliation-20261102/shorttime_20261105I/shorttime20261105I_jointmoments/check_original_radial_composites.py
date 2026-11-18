"""Exact source coefficient tests on the EXISTING R048 quadratic radial caller.
The two labels below are marked source slots of its whole Green operator,
not a new two-mode model or a replacement for the original source.
"""
from pathlib import Path
import hashlib,json
import sympy as s
v,lam=s.symbols('v lam',positive=True)
t,r=s.symbols('t r',real=True)
gxx,gyy,gxy=s.symbols('Gxx Gyy Gxy',real=True)
m1,m2=s.symbols('m1 m2',real=True)
z=s.exp(m1*t+m2*r+(gxx*t*t+2*gxy*t*r+gyy*r*r)/2)
checks=[]
def need(name,result):
 checks.append({'name':name,'passed':bool(result)})
 if not result:raise AssertionError(name)
def moment(a,b):return s.diff(z,t,a,r,b).subs({t:0,r:0})
def Ux(fun):return v*v*fun+s.sqrt(2)*v*s.diff(fun,t)+s.diff(fun,t,2)/2
def Uy(fun):return v*v*fun+s.sqrt(2)*v*s.diff(fun,r)+s.diff(fun,r,2)/2
eu=s.simplify(Ux(z).subs({t:0,r:0}))
want=(v+m1/s.sqrt(2))**2+gxx/2
need('whole_R048_composite_first_moment',s.simplify(eu-want)==0)
euxuy=s.simplify(Ux(Uy(z)).subs({t:0,r:0}))
euy=eu.xreplace({m1:m2,gxx:gyy})
cov=s.simplify(euxuy-eu*euy)
need('whole_R048_composite_two_point',s.simplify(cov-(2*(v+m1/s.sqrt(2))*(v+m2/s.sqrt(2))*gxy+gxy*gxy/2))==0)
need('flat_stationary_background_is_not_quantum_u',s.simplify(eu.subs(m1,0)-v*v)==gxx/2)
need('flat_stationary_background_does_not_kill_u_covariance',s.simplify(cov.subs({m1:0,m2:0})-(2*v*v*gxy+gxy*gxy/2))==0)
eta=s.symbols('eta',real=True)
potential=lam*((v+eta/s.sqrt(2))**2-v*v)**2
need('actual_R048_potential_quadratic_cubic_quartic',s.expand(potential-(2*lam*v*v*eta**2+s.sqrt(2)*lam*v*eta**3+lam*eta**4/4))==0)
force_extra=s.diff(potential-2*lam*v*v*eta**2,eta)
gauss_force=s.expand(force_extra).coeff(eta,2)*moment(2,0)+s.expand(force_extra).coeff(eta,3)*moment(3,0)
need('actual_first_SD_residual_of_Gaussian_replacement',s.simplify(gauss_force.subs(m1,0)-3*s.sqrt(2)*lam*v*gxx)==0)
need('actual_second_SD_residual_of_Gaussian_replacement',s.simplify((3*s.sqrt(2)*lam*v*moment(2,1)+lam*moment(3,1)).subs({m1:0,m2:0})-3*lam*gxx*gxy)==0)

# Pointwise diagonal terms prove the source shift for the WHOLE original K;
# every offdiagonal entry of K is unchanged, no commutation with M_q assumed.
q,k,l=s.symbols('q k l',real=True)
need('whole_composite_source_exact_quadratic_completion',s.expand(-k*eta**2/2+l*eta+q*(v*v+s.sqrt(2)*v*eta+eta**2/2)-(-(k-q)*eta**2/2+(l+s.sqrt(2)*v*q)*eta+v*v*q))==0)
c,F0,R,alpha,theta=s.symbols('c F0 R alpha theta',real=True)
U=v*v+s.sqrt(2)*v*eta+eta**2/2
need('neutral_auxiliary_source_maps_to_composite_source',s.expand(q*U+theta*(F0+c*U+2*alpha*R)-alpha*theta**2-((q+c*theta)*U+theta*(F0+2*alpha*R)-alpha*theta**2))==0)
need('neutral_mixed_u_rho_contact',s.simplify(s.diff((q+c*theta)**2*cov/2-alpha*theta**2,q,theta)-c*cov)==0)
need('neutral_rho_second_contact',s.simplify(s.diff((q+c*theta)**2*cov/2-alpha*theta**2,theta,2)-(c*c*cov-2*alpha))==0)

# General exact logarithmic derivative of the SAME original integral.
o,oi,oj,si,sj,sij=s.symbols('O Oi Oj Si Sj Sij')
z0,zi,zj,zij=s.symbols('Z Zi Zj Zij',nonzero=True)
f=s.Function('Z')(t,r)
need('whole_log_source_second_derivative',s.simplify(s.diff(s.log(f),t,r)-(s.diff(f,t,r)/f-s.diff(f,t)*s.diff(f,r)/f**2))==0)

root=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sources=[
 ('experiments/theory-contracts/tfpt-unified-source-derivation-20261029/aps_source_review/COMMON-CURVATURE-RESPONSE-PROOF.txt',[[117,174]]),
 ('experiments/theory-contracts/tfpt-unified-source-derivation-20261029/aps_source_review/DYNAMICAL-HIGGS-SOURCE-PROOF.txt',[[34,78],[88,116]]),
 ('experiments/theory-contracts/tfpt-unified-source-derivation-20261029/aps_source_review/ORIGINAL-O1-HIGGS-LAW-PROOF.txt',[[1,124]]),
 ('verification/v361_grav_backreaction.py',[[1,40],[60,103]]),
 ('documents/archive--tfpt-45--source_extracts--04_qft_source/sections/007-subsection-renormalized-physical-and-scheme-readout-laye.tex',[[1,28],[119,145]]),
]
out={'passed':sum(c['passed']for c in checks),'total':len(checks),'checks':checks,
 'Atlas_revision_at_assignment':'d1ff953daecc63e28671b479',
 'existing_source':'R048 full fixed-background quadratic radial source, original K_R and original cubic/quartic vertices',
 'new_exact_source_coefficients':{'Eu':'(v+m/sqrt(2))**2+Gxx/2','Cov_u':'2*(v+mx/sqrt(2))*(v+my/sqrt(2))*Gxy+Gxy**2/2','first_SD_residual_if_quadratic_replaces_original':'3*sqrt(2)*lambda*v*Gxx','second_SD_residual_if_quadratic_replaces_original':'3*lambda*Gxx*Gxy'},
 'whole_composite_source':'Z_quad[0]*exp(v**2*<q,1>)*det(I-G*Mq)**(-1/2)*exp((<l+sqrt(2)*v*q,(K-Mq)^(-1)*(l+sqrt(2)*v*q)>-<l,G*l>)/2)',
 'whole_auxiliary_composite_source':'Z_u[q+c*theta]*exp(<theta,F0+2*alpha*R>-alpha*<theta,theta>)',
 'whole_joint_integration_evaluated':False,'no_new_source_family':True,'canonical_edits':False,
 'sources':[{'path':p,'sha256':hashlib.sha256((root/p).read_bytes()).hexdigest(),'lines':lines}for p,lines in sources]}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
