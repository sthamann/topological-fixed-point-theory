"""Exact nontruncated moment reduction for the GIVEN R048 quartic formula.
Numerical test parameters below are mathematical controls, not TFPT data.
"""
from pathlib import Path
import json,functools
import sympy as s
import mpmath as mp
checks=[]
def need(name,ok):
 if not ok:raise AssertionError(name)
 checks.append({'name':name,'passed':True})
J1,J2=s.symbols('J1 J2');b1,b2,a=s.symbols('beta1 beta2 a',nonzero=True);k11,k12,k22,l1,l2=s.symbols('K11 K12 K22 l1 l2');f=s.Function('Z')(J1,J2)
def E(expr,i):
 j=[J1,J2][i];beta=[b1,b2][i];kk=[[k11,k12],[k12,k22]][i];ll=[l1,l2][i]
 return beta*s.diff(expr,j,3)+3*a*beta*s.diff(expr,j,2)+sum(kk[t]*s.diff(expr,[J1,J2][t])for t in range(2))-(ll+j)*expr
need('whole_source_translation_operators_commute_for_actual_symmetric_K',s.expand(E(E(f,1),0)-E(E(f,0),1))==0)
base={idx:s.Symbol('M'+''.join(map(str,idx))) for idx in [(i,j)for i in range(3)for j in range(3)]}
def reducer(reverse):
 @functools.lru_cache(None)
 def red(m):
  if max(m)<3:return base[m]
  candidates=[i for i in range(2)if m[i]>=3];i=candidates[-1] if reverse else candidates[0]
  n=list(m);n[i]-=3;beta=[b1,b2][i];kk=[[k11,k12],[k12,k22]][i];ll=[l1,l2][i]
  rhs=ll*red(tuple(n))
  if n[i]:q=n.copy();q[i]-=1;rhs+=n[i]*red(tuple(q))
  for j in range(2):q=n.copy();q[j]+=1;rhs-=kk[j]*red(tuple(q))
  q=n.copy();q[i]+=2;rhs-=3*a*beta*red(tuple(q))
  return s.expand(rhs/beta)
 return red
r0,r1=reducer(False),reducer(True)
for n in [(3,3),(4,3),(3,4)]:need('exact_reduction_path_consistency_'+str(n),s.simplify(r0(n)-r1(n))==0)
# Whole non-Gaussian quadrature check of the original recurrence, not Wick closure.
mp.mp.dps=45;lam=mp.mpf('.1');vev=mp.mpf('.5');aa=mp.sqrt(2)*vev;K=mp.mpf(100);mean=mp.mpf('.15');sig=1/mp.sqrt(K);ell=K*mean
C=27*lam*aa**4/4
def integrand(z,n):
 eta=mean+sig*z;W=lam/4*(eta+3*aa)**2*((eta-aa)**2+2*aa*aa)
 return eta**n*mp.exp(-W-z*z/2)/mp.sqrt(2*mp.pi)
M=[mp.quad(lambda z:integrand(z,n),[-mp.inf,-8,0,8,mp.inf])for n in range(13)]
errs=[]
for n in range(10):
 residual=lam*M[n+3]-( (n*M[n-1] if n else 0)-K*M[n+1]-3*aa*lam*M[n+2]+ell*M[n])
 err=abs(residual);errs.append(float(err));need('nonGaussian_full_original_form_SD_recurrence_n'+str(n),err<mp.mpf('1e-30'))
result={'scope':'Actual R048 finite quartic formula; numerical controls are mathematical only, not native TFPT values','whole_source_moment_basis_size_upper_bound':'3^N','base_exponents_per_cell':[0,1,2],'no_moment_cumulant_truncation':True,'cubic_and_quartic_kept':True,'checks':checks,'passed':len(checks),'total':len(checks),'numerical_control_max_SD_residual':max(errs),'native_primitive_source_evaluated':False,'whole_primitive_S_extra_domain_ref_closure_proved':False}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
