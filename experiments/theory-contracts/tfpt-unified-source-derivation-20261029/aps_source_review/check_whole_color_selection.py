"""R068: original whole graph, actual builder, independent exterior/Fock test."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import sys,json
import numpy as np
from scipy.linalg import eigh,expm,block_diag
r=Path(__file__).resolve().parents[4];sys.path.insert(0,str(r/'verification'))
import v252_full_finite_triple as v
checks=[]
def need(name,x,tol=2e-12):
 if abs(x)>=tol:raise RuntimeError((name,x,tol))
 checks.append(dict(name=name,error=float(abs(x)),tolerance=tol))
rng=np.random.default_rng(252)
ys=[rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)) for _ in range(4)]
m=rng.normal(size=(3,3));mr=m+m.T;z=np.zeros((3,3),complex)
dy=v.full_D(*ys,z);dm=v.full_D(z,z,z,z,mr);gamma=v.grading()
iq=[i for i,t in enumerate(v.P) if t[0] in ['u','d']]
gam=gamma[np.ix_(iq,iq)];d=36;eye=np.eye(d)
def lift(g):
 raw=v.rep(1,np.eye(2),g.conj());ks=v.Kswap()
 return (raw@ks@raw.conj()@ks)[np.ix_(iq,iq)]
ha=np.array([[.2,.3+.1j,0],[.3-.1j,-.1,.2j],[0,-.2j,-.1]])
hb=np.array([[.1,-.2j,.15],[.2j,.3,.1j],[.15,-.1j,-.4]])
ga=expm(1j*ha);gb=expm(1j*hb);ua=lift(ga);ub=lift(gb)
t=lift(np.diag([-1,-1,1]));gen=np.zeros(d)
for i,oi in enumerate(iq):gen[i]=[1,-1,0][v.P[oi][2]]
gen=np.diag(gen);zeta=np.exp(2j*np.pi/3)
need('original color matrices are not a commuting-only diagnostic',max(0,.01-np.linalg.norm(ga@gb-gb@ga)),1e-14)
for hh in [0.,.13,.6]:
 ph=(hh*dy+dm)[np.ix_(iq,iq)]
 for p in [.25,-.75]:
  h=p*gam-1j*gam@ph;ev,vec=eigh(h)
  cvals=np.exp(-np.logaddexp(0,4*np.pi*ev))
  sq=(vec*np.sqrt(cvals))@vec.conj().T
  sp=(vec*np.sqrt(1-cvals))@vec.conj().T
  c=(vec*cvals)@vec.conj().T
  e=np.vstack([sp,sq]);g=e@e.conj().T
  need(f'whole original graph projection H{hh} p{p}',np.linalg.norm(g@g-g),1e-11)
  need(f'original second-face marginal H{hh} p{p}',np.linalg.norm(g[d:,d:]-c))
  for label,ul,ur in [('diagonal',ua,ua),('relative',ua,ub),('oneface',eye,t)]:
   joint=np.linalg.det(e.conj().T@block_diag(ul,ur)@e)
   relative=np.linalg.det(eye-c+c@ul.conj().T@ur)
   need(f'whole joint-to-relative {label} H{hh} p{p}',joint-relative,2e-11)
   if label=='diagonal':need(f'whole diagonal color singlet H{hh} p{p}',joint-1,2e-11)
  udiag=block_diag(lift(zeta*np.eye(3)),lift(zeta*np.eye(3)))
  need(f'whole center commutes graph H{hh} p{p}',np.linalg.norm(udiag@g-g@udiag),1e-11)
  need(f'whole center Slater determinant H{hh} p{p}',np.linalg.det(e.conj().T@udiag@e)-1,3e-11)
  # Wick covariance of disjoint-face bilinears is -Tr(A G B G).
  al=block_diag(gen,np.zeros((d,d)));ar=block_diag(np.zeros((d,d)),gen)
  kl=np.trace(al@al@g-al@g@al@g).real
  kr=np.trace(ar@ar@g-ar@g@ar@g).real
  cross=-np.trace(al@g@ar@g).real
  k=np.trace(gen@gen@c@(eye-c)).real
  need(f'left current variance H{hh} p{p}',kl-k,4e-12)
  need(f'right current variance H{hh} p{p}',kr-k,4e-12)
  need(f'exact opposite-face connected current H{hh} p{p}',cross+k,4e-12)
  need(f'whole total color has zero variance H{hh} p{p}',kl+kr+2*cross,5e-12)

# Independent explicit 3-particle wedge in a 6-orbital/64-state Fock space.
# This is one existing color-degenerate eigenvalue triplet, not a new model.
c=F(1,5);ef=np.vstack([np.sqrt(float(1-c))*np.eye(3),np.sqrt(float(c))*np.eye(3)])
occ=list(combinations(range(6),3));amp={sum(1<<i for i in ii):np.linalg.det(ef[list(ii)]) for ii in occ}
need('independent64 Fock Slater normalization',sum(abs(a)**2 for a in amp.values())-1)
# Fock action on 3-particle sector constructed from minors, not graph determinant.
uf=block_diag(ga,ga)
vv=np.array([amp[sum(1<<i for i in ii)] for ii in occ])
wedge=np.array([[np.linalg.det(uf[np.ix_(ii,jj)]) for jj in occ] for ii in occ])
need('independent64 Fock all-SU3 diagonal singlet',np.linalg.norm(wedge@vv-vv))
probs={}
for state,a in amp.items():
 right=state>>3;probs[right]=probs.get(right,0)+abs(a)**2
one=0;center_weight=0;projected=0;global_center=0
for st,a in amp.items():
 nr=(st>>3).bit_count();nl=(st&7).bit_count()
 phase=(-1)**((st>>3&1)+(st>>4&1));prob=abs(a)**2
 one+=prob*phase
 if nr%3==0:center_weight+=prob;projected+=prob*phase
 global_center+=prob*zeta**(nl+nr)
need('whole singlet oneface nontrivial counting',one-float(F(9,25)))
need('whole total triality fixes state',global_center-1)
need('oneface postselection changes normalized counting',projected/center_weight-1)
# Exact rational verification, including normalization and conditioned state.
weights={st:c**st.bit_count()*(1-c)**(3-st.bit_count()) for st in range(8)}
raw=sum(pr*((-1)**((st&1)+((st>>1)&1))) for st,pr in weights.items())
den=sum(pr for st,pr in weights.items() if st.bit_count()%3==0)
num=sum(pr*((-1)**((st&1)+((st>>1)&1))) for st,pr in weights.items() if st.bit_count()%3==0)
assert raw==F(9,25) and den==F(13,25) and num/den==1
checks.append(dict(name='exact oneface center conditioning counterexample',raw=str(raw),weight=str(den),conditioned=str(num/den)))
# Restricting the state to center-invariant observables leaves this expectation.
# Group averaging a STATE by conjugation differs from charged-sector conditioning.
need('original Gstar commutes color center',np.linalg.norm(np.diag([-1,-1,1])@(zeta*np.eye(3))-(zeta*np.eye(3))@np.diag([-1,-1,1])))
need('quotient by fixed original reference cancels in normalized state',(F(7,11)/F(13,17))/(F(5,3)/F(13,17))-F(7,11)/F(5,3))
# First decisive representation test of a SAME-SPACE center compression.
raw=v.rep(1,np.eye(2),(zeta*np.eye(3)).conj());ks=v.Kswap()
zc=raw@ks@raw.conj()@ks
pc=(np.eye(v.NH)+zc+zc@zc)/3
need('original96 oneparticle center is projector',np.linalg.norm(pc@pc-pc))
need('same-space center retains24 lepton slots',np.trace(pc)-24)
need('same-space center kills36 particle quark slots',np.linalg.norm(pc[:,iq]))
need('same-space center kills36 antiquark slots',np.linalg.norm(pc[:,[i+v.NP for i in iq]]))
need('original quark center is fundamental',np.linalg.norm(zc[np.ix_(iq,iq)]-zeta*np.eye(36)))
lept=[i for i in range(v.NH) if i%v.NP not in iq]
need('same-space center retains all actual lepton slots',np.linalg.norm(pc[np.ix_(lept,lept)]-np.eye(24)))
out=dict(verdict='SCOPED_WHOLE_SOURCE_COLOR_SELECTION_AND_RELATIVE_READER',controls=len(checks),checks=checks,
 exact_counterexample=dict(single_face='9/25',total_center='1',conditioned_single_face='1',conditioning_probability='13/25'),
 scope='Existing conditional whole G quark source and actual pi/J color lift; original carrier-to-source admissibility map is NOT assumed; no new H measure, physical Gauss law or spacetime derivation.')
Path(__file__).with_name('WHOLE-COLOR-SELECTION.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(dict(verdict=out['verdict'],controls=len(checks))))
