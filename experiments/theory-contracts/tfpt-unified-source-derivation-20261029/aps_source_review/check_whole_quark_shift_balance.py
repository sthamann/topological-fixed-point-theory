"""R084: full original quark covariance under a conditional common shift.

The scalar shift on this carrier is a hypothesis; no Majorana Pfaffian,
complete E8 source, physical spacetime, or H-integrated state is computed.
"""
from pathlib import Path
import json,sys
import numpy as np
from scipy.special import expit

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
I=np.eye(len(ix));keep=np.where(np.diag(g).real>0)[0];drop=np.where(np.diag(g).real<0)[0]
mass=.3*np.linalg.svd(t[np.ix_(keep,drop)],compute_uv=False)
ell=.2;z=.3; alpha=.25;checks={};metrics={}
def need(name,condition):
 checks[name]=bool(condition)
 if not condition: raise RuntimeError(name)
def covariance(p,r=.3):
 h=p*g+r*t
 e,u=np.linalg.eigh(h)
 return (u*expit(2*ell*e))@u.conj().T
def logblock(p):
 e=np.sqrt(p*p+mass*mass)
 a=expit(2*ell*e)
 return float(np.log1p(-z*a).sum()+np.log1p(-z*(1-a)).sum())
need('actual_quark_and_grading_dimensions',len(ix)==36 and len(keep)==len(drop)==18)
need('actual_force_is_hermitian_odd',np.linalg.norm(t-t.conj().T)<1e-10 and np.linalg.norm(t@g+g@t)<1e-10)
for p in [.25,1.25,-.75]:
 C=covariance(p)
 e=np.linalg.eigvalsh(C)
 need('actual_full_block_spectrum_'+str(p),np.max(abs(e[:18]+e[18:][::-1]-1))<1e-12)
 sign,ld=np.linalg.slogdet(I-z*C)
 need('actual_full_block_formula_'+str(p),abs(ld-logblock(p))<1e-11 and abs(sign-1)<1e-12)
 Cplus=C[np.ix_(keep,keep)];Cminus=C[np.ix_(drop,drop)]
 blockprod=np.linalg.slogdet(np.eye(18)-z*Cplus)[1]+np.linalg.slogdet(np.eye(18)-z*Cminus)[1]
 metrics['nonfactorizing_offdiagonal_'+str(p)]=float(abs(ld-blockprod))
 need('full_covariance_cannot_be_replaced_by_both_marginals_'+str(p),abs(ld-blockprod)>1e-5)

# Show the local Higgs dependence remains present in the same full source.
dynamic=np.linalg.norm(covariance(.25)-covariance(.25,0))
need('actual_full_quark_covariance_changes_with_H',dynamic>1e-3)
metrics['actual_covariance_H_response_norm']=float(dynamic)

# Direct actual 36x36 tails: raw polarization is O(1/n), shift difference
# is O(1/n^2). These samples illustrate; the accompanying proof is analytic.
raw=[];shift=[]
for n in [80,160,320]:
 C=covariance(n+alpha);Cm=covariance(n+alpha-1)
 P=(I+g)/2
 raw.append(float(np.linalg.norm(C-P,ord='nuc')))
 shift.append(float(np.linalg.norm(C-Cm,ord='nuc')))
need('raw_polarization_has_one_over_n_tail',1.8<raw[0]/raw[1]<2.2 and 1.8<raw[1]/raw[2]<2.2)
need('shift_difference_has_one_over_n_squared_tail',3.6<shift[0]/shift[1]<4.4 and 3.6<shift[1]/shift[2]<4.4)
metrics['raw_trace_norms']=raw;metrics['shift_trace_norms']=shift

# Relative determinant of the 36x36 BLOCKS. Open Fourier restriction, not a
# fictitious cyclic shift at the cutoff. Telescope first, then take N→∞.
samples={}
for alpha in [.25,.75]:
 ratios=[]
 for N in [4,10,40]:
  blocksum=sum(logblock(n+alpha-1)-logblock(n+alpha) for n in range(-N,N+1))
  endpoint=logblock(-N+alpha-1)-logblock(N+alpha)
  need('full_telescope_'+str(alpha)+'_'+str(N),abs(blocksum-endpoint)<1e-10)
  ratios.append(float(np.exp(blocksum)))
 need('whole_relative_determinant_tends_one_'+str(alpha),abs(ratios[-1]-1)<1e-5)
 samples[str(alpha)]=ratios
metrics['whole_relative_determinant_samples']=samples

# A strictly diagonal chirality consumer has opposite index powers; the
# full interacting determinant was computed independently above.
need('two_chiral_ranks_match',len(keep)-len(drop)==0)
out={'scope':'actual finite 36-quark constant-H source with stipulated same scalar circle shift; ordinary relative determinant; no full E8/physical Pfaffian',
     'checks':checks,'passed':len(checks),'failed':0,'metrics':metrics,
     'source_seed':252,'ell_control':ell,'z_control':z,'ledger_changed':False}
Path(__file__).with_name('WHOLE-QUARK-SHIFT-BALANCE.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
