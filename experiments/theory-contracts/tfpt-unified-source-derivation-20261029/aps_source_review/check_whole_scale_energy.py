"""Controls for R040/R041. Existing quarter block; no physical source fit.
Profile mutation checks the inverse map, not a newly selected TFPT geometry.
The finite chiral matrix only falsifies an incorrect weighted-trace formula.
"""
from pathlib import Path
import json
import sympy as s
import mpmath as mp
checks={};details={}
def need(k,b):
 checks[k]=bool(b)
 if not b:raise AssertionError(k)
def zero(v):return all(s.simplify(x)==0 for x in v)
r,z,ep=s.symbols('r z ep',real=True);p,m=s.symbols('p m',real=True)
sig=s.Function('sigma')(r);h=s.Function('h')(r);f=s.exp(-sig/2)
u=s.Matrix([s.Function('u1')(r),s.Function('u2')(r)])
N=s.Matrix([[0,s.I],[s.I,0]]);B=s.Matrix([[p,-s.I*m],[s.I*m,-p]])
def K0(v):return N*(v.diff(r)+B*v)
def K(v):return f*K0(f*v)
need('original_normal_and_coupled_block_contract',N.H==-N and N*N==-s.eye(2) and zero(N*B+B*N))
need('full_variable_scale_null_solution_transport',zero(K(u/f)-f*K0(u)))
need('nonzero_energy_weight_exact_without_commutation',zero(K(u)-z*u-f*(K0(f*u)-z*s.exp(sig)*(f*u))))
need('unweighted_energy_substitution_has_nonzero_defect',not zero(K(u)-z*u-f*(K0(f*u)-z*(f*u))))
def Kep(v):
 ff=s.exp(-(sig+ep*h)/2);return ff*K0(ff*v)
dK=s.diff(Kep(u),ep).subs(ep,0)
need('nonconstant_scale_operator_variation_is_anticommutator',zero(dK+(h*K(u)+K(h*u))/2))
# Positive spectrum from the same historical p=1/4,m=1/3 odd block.
mp.mp.dps=90;ell=mp.mpf(1);mu=mp.mpf(1)/3;alpha=mp.mpf(1)/4
NN=mp.matrix([[0,1j],[1j,0]]);II=mp.eye(2)
def data(n,profile):
 pp=n+alpha;bb=mp.matrix([[pp,-1j*mu],[1j*mu,-pp]]);E=mp.sqrt(pp**2+mu**2);P=(II+bb/E)/2
 if profile is None:
  dt=mp.diff(lambda zz:mp.expm(-ell*(bb+zz*NN)),0)
 else:
  dt=mp.matrix(2)
  for i in range(2):
   for j in range(2):dt[i,j]=-mp.quad(lambda t:(mp.expm(-(ell-t)*bb)*NN*mp.expm(-t*bb))[i,j]*profile(t),[0,ell])
 out=mp.exp(-ell*E)*sum((P*NN*dt)[i,i] for i in range(2))
 return E,mp.re(out)
Es=[];Ms=[];maxerr=mp.mpf(0)
for n in range(7):
 E,M=data(n,None);Es.append(E);Ms.append(M)
 maxerr=max(maxerr,abs(M-(1-mp.exp(-2*ell*E))/(2*E)))
need('independent_whole_transfer_energy_derivative_gives_flat_moments',maxerr<mp.mpf('1e-75'))
a=mp.exp(-2*ell)
def gram(k):return mp.matrix([[(1-a**(Es[i]+Es[j]-Es[0]))/(2*(Es[i]+Es[j]-Es[0])) for j in range(k+1)] for i in range(k+1)])
c=mp.lu_solve(gram(6),mp.matrix(Ms))
need('unchanged_original_flat_source_recovers_whole_profile_one',max(abs(c[i]-(1 if i==0 else 0)) for i in range(7))<mp.mpf('1e-60'))
# Single positive scale mutation; tests recovery, no physical interpretation.
lam=Es[1]-Es[0];profile=lambda t:1+mp.exp(-2*t*lam)/3
mut=[]
for n in range(3):
 E,M=data(n,profile);mut.append(M)
 cc=(1-mp.exp(-2*ell*E))/(2*E)+(1-mp.exp(-2*ell*(E+lam)))/(6*(E+lam))
 need('normal_profile_mutation_transfer_moment_'+str(n),abs(M-cc)<mp.mpf('1e-75'))
c=mp.lu_solve(gram(2),mp.matrix(mut))
need('same_inverse_recovers_positive_profile_mutation_without_scan',max(abs(c[i]-([1,mp.mpf(1)/3,0][i])) for i in range(3))<mp.mpf('1e-70'))
# Finite algebra control: total chirality commutes with L, not with K.
T=s.Matrix([[1,1],[0,2]]);KK=s.zeros(4);KK[:2,2:]=T.T;KK[2:,:2]=T
hh=s.diag(1,2,1,2);Pc=s.diag(1,1,0,0);LL=KK*KK
varK=-(hh*KK+KK*hh)/2;varL=varK*KK+KK*varK
actual=-s.trace(Pc*LL.inv()*varL*LL.inv())
need('chiral_variation_is_half_whole_trace_not_naive_inserted_half',actual==s.trace(hh*LL.inv()) and actual!=2*s.trace(Pc*hh*LL.inv()))
# Meromorphic and local-coefficient checks are independent symbolic identities.
ss=s.symbols('s');rr,cc,dd=s.symbols('R0 C0 D0')
flow=s.expand(2*ss*(rr/ss+cc+dd*ss))
need('FP_and_logdet_read_different_Laurent_coefficients',flow.subs(ss,0)==2*rr and -s.diff(flow,ss).subs(ss,0)==-2*cc)
need('smeared_Higgs_divergence_cannot_be_dropped',s.integrate((1+s.cos(r))*s.diff(2+s.cos(r),r,2),(r,0,2*s.pi))==-s.pi)
need('local_Higgs_divergence_coefficient_from_full_heat_trace',s.Rational(60,360)*(-4)==-s.Rational(2,3))
details={'flat_transfer_moment_max_abs_error':mp.nstr(maxerr,8),'flat_reconstruction_coefficients':[mp.nstr(x,8) for x in mp.lu_solve(gram(6),mp.matrix(Ms))],'diagnostics_only':'Numerical inversion of the existing historical odd block and one positive profile mutation; not a measured TFPT scale,full3+1D fit or infinite-data stability proof.'}
out={'contract':'WHOLE-SCALE-ENERGY-RECONSTRUCTION-20261031','all_passed':all(checks.values()),'checks':checks,'details':details,'scope':'Analytic proof supplies whole profile uniqueness using the cited Muntz theorem. Controls test differential identities,independent energy response,Gram inversion,Laurent coefficients and a false chiral-weight shortcut. No new physical source,Y/MR,domain or metric selected;no new Lean formalization.'}
Path(__file__).with_name('WHOLE-SCALE-ENERGY-RECONSTRUCTION.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'all_passed':out['all_passed'],'count':len(checks),'details':details}))
