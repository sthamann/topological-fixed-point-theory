"""R084 original Higgs singular pair: asymmetrically shifted chirality.

Diagnostic U_plus is not asserted to be the original seam action. Test the
actual Schatten-class issue before taking an ordinary Fredholm determinant.
"""
from pathlib import Path
import json
import numpy as np
from scipy.special import expit

rng=np.random.default_rng(252)
ys=[rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)) for _ in range(4)]
s=float(np.linalg.svd(ys[0],compute_uv=False)[0]);m=.3*s;ell=.2;alpha=.25
def cross(n, mass=m):
 p=n+alpha; E=np.sqrt(p*p+mass*mass)
 return .5*mass/E*np.tanh(ell*E)
def full(n):
 p=n+alpha; E=np.sqrt(p*p+m*m);fac=.5*np.tanh(ell*E)/E
 return np.eye(2)/2+fac*np.array([[p,m],[m,-p]])
checks={}
def need(name,cond):
 checks[name]=bool(cond)
 if not cond: raise RuntimeError(name)

# Finite sections of S diag(b)-diag(b), with no cyclic wrap-around.
def partial_section(N,mass=m):
 n=np.arange(-N,N+1,dtype=float);b=cross(n,mass);A=-np.diag(b)
 A[1:,:-1]+=np.diag(b[:-1])
 return A,b
samples=[]
for N in [20,80,320]:
 A,b=partial_section(N)
 trace_norm=float(np.linalg.norm(A,ord='nuc'))
 diag_lower=float(np.sum(abs(b)))
 hs=float(np.linalg.norm(A,ord='fro'))
 need('pinching_lower_'+str(N),trace_norm+1e-10>=diag_lower)
 samples.append({'N':N,'trace_norm':trace_norm,'diagonal_lower_bound':diag_lower,'HS_norm':hs})
need('asymmetric_trace_lower_bound_grows',samples[-1]['diagonal_lower_bound']>samples[0]['diagonal_lower_bound'])
need('asymmetric_HS_bounded',samples[-1]['HS_norm']<2*samples[0]['HS_norm'])
need('massless_offdiagonal_is_zero',np.linalg.norm(partial_section(20,0)[0])==0)
small=cross(10000);need('asymptotic_coefficient',abs(10000*small-m/2)<1e-4)
full_diffs=[np.linalg.norm(full(n)-full(n-1),ord='nuc') for n in [80,160,320]]
need('common_shift_tail_summable',all(3.5<full_diffs[i]/full_diffs[i+1]<4.5 for i in range(2)))
out={'scope':'R081 Yu singular pair with hypothesized shift on gamma+ only; not an actual TFPT selector',
     'seed':252,'mass':m,'ell':ell,'alpha':alpha,'samples':samples,
     'common_shift_tail':full_diffs,'checks':checks,'passed':len(checks)}
Path(__file__).with_name('ASYMMETRIC-SHIFT-CLASS.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
