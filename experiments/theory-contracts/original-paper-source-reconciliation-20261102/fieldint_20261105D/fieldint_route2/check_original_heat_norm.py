"""Targeted controls of the actual whole-circle sea-corrected heat class.
Original R055 fixtures only, unchanged v252 producer. Not a physical measure.
"""
from pathlib import Path
import sys,json
import numpy as np
import mpmath as mp
root=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0,str(root/'verification'))
import v252_full_finite_triple as v
mp.mp.dps=70
# Actual canonical elliptic branch ell=2pi; original heat convention t=2ell.
t=4*np.pi;tm=4*mp.pi
# Identical preexisting R055 diagnostic fixture, no new model or coupling.
rng=np.random.default_rng(252)
Ys=[rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)) for _ in range(4)]
mr=rng.normal(size=(3,3));MR=mr+mr.T
z3=np.zeros((3,3),complex);z2=np.zeros((2,2),complex)
DY=v.full_D(*Ys,z3);DM=v.full_D(z3,z3,z3,z3,MR)
P=v.rep(1,z2,z3);K=v.Kswap();gamma=v.grading()
def phi(H):
 h1,h2=H;q=np.array([[np.conj(h2),h1],[-np.conj(h1),h2]])
 a=v.rep(0,q,z3);B=a@(DY@P-P@DY);B=B+B.conj().T
 return B+K@B.conj()@K+DM
checks=[];native_errors=[];rows=[]
def check(name,ok):
 checks.append({'name':name,'passed':bool(ok)})
 if not ok:raise AssertionError(name)
def free_log(alpha,N):
 return 96*sum(mp.log1p(mp.exp(-tm*abs(mp.mpf(n)+alpha))) for n in range(-N,N+1))
def tail(alpha,N):
 return 96*(mp.exp(-tm*(N+1+alpha))+mp.exp(-tm*(N+1-alpha)))/(1-mp.exp(-tm))
ceilings=[]
for aa in [mp.mpf(1)/4,mp.mpf(3)/4]:
 a=float(aa);logC=free_log(aa,32);tr=tail(aa,32)
 base_sum=(mp.exp(-tm*aa)+mp.exp(-tm*(1-aa)))/(1-mp.exp(-tm))
 check(f'ceiling_geometric_majorant_alpha_{aa}',logC<=96*base_sum)
 # Independent positive products at wider existing mode partitions.
 for N in [0,1,3]:
  increment=free_log(aa,32)-free_log(aa,N)
  check(f'whole_log_tail_alpha_{aa}_N_{N}',0<=increment<=tail(aa,N))
 ceilings.append({'alpha':str(aa),'time':'4*pi = 2*(original ell=2*pi)',
  'log_C_truncated_approx':mp.nstr(logC,50),
  'C_truncated_approx':mp.nstr(mp.exp(logC),50),
  'numeric_arithmetic':'70 decimal digit floating arithmetic, not interval enclosure',
  'log_tail_N32_bound':mp.nstr(tr,18)})
 for radius in [0.,1.,64.]:
  ph=phi((0j,complex(radius)));evp=np.linalg.eigvalsh(ph)
  LF=sum(float(np.logaddexp(0.,-t*np.sqrt((n+a)**2+evp**2)).sum()) for n in range(-3,4))
  LC=float(free_log(aa,3))
  check(f'whole_finite_heat_bounds_a{a}_r{radius}',LF>=0 and LF<=LC+1e-10)
  rows.append({'alpha':a,'original_R055_radius':radius,'N':3,'log_F':LF,'log_free_ceiling':LC})
  for n in [-1,0,1]:
   p=n+a;h=p*gamma-1j*gamma@ph
   eig=np.linalg.eigvalsh(h)
   # Native full96 h spectrum vs complete Phi squared spectrum.
   native=float(np.logaddexp(0.,-t*np.abs(eig)).sum())
   reduced=float(np.logaddexp(0.,-t*np.sqrt(p*p+evp*evp)).sum())
   err=abs(native-reduced);native_errors.append(err)
   check(f'native96_heat_source_a{a}_r{radius}_n{n}',err<1e-9)
check('whole_alpha_quarter_symmetry',abs(free_log(mp.mpf(1)/4,32)-free_log(mp.mpf(3)/4,32))<mp.mpf('1e-65'))
out={'scope':'Whole actual constant-circle sea-corrected coordinate Fock heat class only. Actual alpha1/4,3/4 and original elliptic t=4pi; original R055 seed252 fixtures at original radii0,1,64. No new physical input, no original H integration or full aA/aR evaluation, no physical Padm replacement.',
 'new_premise_against_R055_R056':'Uses the already evaluated whole sea-corrected infinite fermion heat determinant, not only finite polynomial Pfaffian coefficients. The whole original additional aA/aR remains retained and unevaluated.',
 'checked':sum(c['passed'] for c in checks),'total':len(checks),'native_comparisons':len(native_errors),
 'max_absolute_error':max(native_errors),'checks':checks,'ceilings':ceilings,'whole_mode_controls':rows,
 'primitive_common_field_integral_evaluated':False}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['checked','total','native_comparisons','max_absolute_error','ceilings','scope']},indent=2))
