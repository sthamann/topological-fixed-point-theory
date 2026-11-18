"""R084: det2 of full selectively shifted original R081 singular pair.

Open Fourier cutoffs compress the infinite relative operator; no periodic
cutoff shift. This diagnostic is not the original TFPT seam or normalization.
"""
from pathlib import Path
import json
import numpy as np

rng=np.random.default_rng(252)
ys=[rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)) for _ in range(4)]
m=.3*float(np.linalg.svd(ys[0],compute_uv=False)[0]);ell=.2;alpha=.25;z=.3
def finite_section(N,mass):
 n=np.arange(-N,N+1);p=n+alpha;E=np.sqrt(p*p+mass*mass)
 f=.5*np.tanh(ell*E)/E;a=.5+p*f;d=.5-p*f;b=mass*f
 pm=p-1;em=np.sqrt(pm*pm+mass*mass);am=.5+pm*.5*np.tanh(ell*em)/em
 size=len(n);C=np.zeros((2*size,2*size));D=np.zeros_like(C)
 C[:size,:size]=np.diag(a);C[size:,size:]=np.diag(d)
 C[:size,size:]=np.diag(b);C[size:,:size]=np.diag(b)
 D[:size,:size]=np.diag(am-a)
 B=-np.diag(b);B[1:,:-1]+=np.diag(b[:-1])
 D[:size,size:]=B;D[size:,:size]=B.T
 X=-z*D@np.linalg.inv(np.eye(2*size)-z*C)
 sign,logdet=np.linalg.slogdet(np.eye(2*size)+X)
 return float(sign*np.exp(logdet-np.trace(X))),float(np.trace(X))

checks={}
def need(k,v):
 checks[k]=bool(v)
 if not v: raise RuntimeError(k)
rows=[]
for N in [8,24,64,128]:
 massive,tr=finite_section(N,m);massless,tr0=finite_section(N,0)
 rows.append({'N':N,'massive_det2':massive,'massless_det2':massless,
              'massive_finite_trace':tr,'massless_finite_trace':tr0})
 need('finite_nonzero_'+str(N),massive>0 and massless>0)
need('genuine_mass_dependence',abs(rows[-1]['massive_det2']-rows[-1]['massless_det2'])>.01)
need('massive_cutoff_difference_reduces',abs(rows[-1]['massive_det2']-rows[-2]['massive_det2'])<abs(rows[1]['massive_det2']-rows[0]['massive_det2']))
need('massless_cutoff_difference_reduces',abs(rows[-1]['massless_det2']-rows[-2]['massless_det2'])<1e-8)
out={'scope':'det2 of hypothesized chiral-only shift on actual R081 finite-H singular pair; no physical source or normalization identified',
     'z':z,'alpha':alpha,'ell':ell,'mass_control':m,'cutoffs':rows,'checks':checks,'passed':len(checks)}
Path(__file__).with_name('REGULARIZED-SELECTIVE-SOURCE.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
