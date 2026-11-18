"""Whole R081 retained quark space: winding multiplicity under scalar shift.

The shared scalar shift and the Fock determinant are conditional diagnostics;
this is not the physical Majorana Pfaffian or the full TFPT source.
"""
from pathlib import Path
import sys, json
import numpy as np
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'verification'))
import v252_full_finite_triple as original
rng=np.random.default_rng(252)
ys=[rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)) for _ in range(4)]
mr=rng.normal(size=(3,3));D=original.full_D(*ys,mr+mr.T)
gamma=original.grading();swap=original.Kswap()
bp=original.rep(1,np.zeros((2,2)),np.zeros((3,3)))
a=original.rep(0,np.eye(2),np.zeros((3,3)))
om=a@(D@bp-bp@D);sa=om+om.conj().T;V=sa+swap@sa.conj()@swap
ix=[i for i,f in enumerate(original.P) if f[0] in ('u','d')]
g=gamma[np.ix_(ix,ix)];t=-1j*g@V[np.ix_(ix,ix)]
keep=np.where(np.diag(g).real>0)[0];drop=np.where(np.diag(g).real<0)[0]
cross=t[np.ix_(keep,drop)];sing=np.linalg.svd(cross,compute_uv=False)
checks={}
def require(name,c):
 checks[name]=bool(c)
 if not c: raise RuntimeError(name)
require('actual_R081_retained_rank18',len(keep)==18 and len(drop)==18)
require('actual_radial_odd_force',np.linalg.norm(t@g+g@t)<1e-10)
require('actual_18_channel_singular_values',len(sing)==18)
mass=.3*sing; ell=.2;alpha=.25;z=.03
def occ(p):
 e=np.sqrt(p*p+mass*mass)
 return (1+p/e*np.tanh(ell*e))/2
values=[]
for n in [100,1000,10000]:
 logdet=float(np.sum(np.log1p(-z*occ(-n+alpha-1))-np.log1p(-z*occ(n+alpha))))
 values.append({'cutoff':n,'logdet':logdet,'det':float(np.exp(logdet))})
target=(1-z)**(-18)
require('retained_18_channel_winding_limit',abs(values[-1]['det']-target)<1e-7)
require('cannot_count_as_single_index_factor',abs(target-1/(1-z))>.1)
out={'scope':'full 18-mode diagnostic gamma+ quark component with hypothetical same scalar shift, fixed background; not physical Pfaffian',
 'channels':len(keep),'winding_index':-len(keep),'z':z,'target':target,'single_channel_target':1/(1-z),
 'samples':values,'checks':checks,'passed':len(checks)}
Path(__file__).with_name('WINDING-MULTIPLICITY.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
