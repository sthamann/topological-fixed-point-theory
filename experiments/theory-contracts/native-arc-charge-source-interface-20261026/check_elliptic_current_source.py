"""Current source functional of the existing elliptic Gaussian boundary.
No fitted temperature/state/coupling. Original Hardy source remains distinct.
"""
import argparse,hashlib,json
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.linalg import expm
from scipy.special import expit
import sympy as s
REPO=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
SRC=REPO/'experiments/theory-contracts/elliptic-moebius-source-bridge-20261026'
def need(ok,msg):
 if not bool(ok):raise RuntimeError(msg)
def multiply(a,b,cut):
 out={}
 for (e,q),v in a.items():
  for (f,r),w in b.items():
   if e+f<=cut:out[e+f,q+r]=out.get((e+f,q+r),0)+v*w
 return out

def charge_partition(sign,cut=80):
 fermion={(0,0):1}
 for r in range(cut//4+1):
  for power,charge in [(4*r+2-sign,1),(4*r+2+sign,-1)]:
   if power<=cut:fermion=multiply(fermion,{(0,0):1,(power,charge):1},cut)
 osc={(0,0):1}
 for n in range(1,cut//4+1):osc=multiply(osc,{(4*n*k,0):1 for k in range(cut//(4*n)+1)},cut)
 qs={}
 for q in range(-cut,cut+1):
  e=2*q*q-sign*q
  if e<=cut:qs[e,q]=1
 boson=multiply(osc,qs,cut)
 need(fermion==boson,'Exact charge-resolved elliptic fermion/boson partition polynomial')
 return len(fermion)

def symbolic_identities():
 x,y=s.symbols('x y',positive=True)
 c=1/(1+x);d=1/(1+x*y)
 need(s.factor(c*(1-d)-(c-d)/(1-1/y))==0,'Fermi telescoping positive mode identity')
 need(s.factor(d*(1-c)-(c-d)/(y-1))==0,'Fermi telescoping negative mode identity')
 r=s.symbols('r',integer=True)
 for epsilon in (0,1):
  alpha=s.Rational(2*epsilon+1,4)
  need(s.expand(2*s.pi*(2*r+epsilon+s.Rational(1,2))-4*s.pi*(r+alpha))==0,'Actual NS cover to quarter inverse temperature')

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,default=Path(__file__).with_name('elliptic_current_source_results.json'));args=parser.parse_args()
 mp.mp.dps=140;b=4*mp.pi
 symbolic_identities();formal=[charge_partition(sign) for sign in (1,-1)]
 rows=[]
 for alpha in (mp.mpf(1)/4,mp.mpf(3)/4):
  C=lambda r:1/(1+mp.exp(b*(r+alpha)))
  theta=sum(mp.exp(-b*(mp.mpf(q*q)/2+(alpha-mp.mpf('.5'))*q)) for q in range(-12,13))
  mean=sum(q*mp.exp(-b*(mp.mpf(q*q)/2+(alpha-mp.mpf('.5'))*q)) for q in range(-12,13))/theta
  fermionZ=mp.fprod(1+mp.exp(-b*(r+alpha)) for r in range(40))*mp.fprod(1+mp.exp(-b*(r-alpha)) for r in range(1,41))
  bosonZ=theta/mp.fprod(1-mp.exp(-b*n) for n in range(1,41))
  need(abs(fermionZ-bosonZ)<mp.mpf('1e-85'),'Actual geometric partition product')
  need(abs(sum(C(r)-(1 if r<0 else 0) for r in range(-40,41))-mean)<mp.mpf('1e-85'),'Actual sea-relative charge mean')
  for m in (1,2,5,10):
   positive=sum(C(r)*(1-C(r+m)) for r in range(-40,41))
   negative=sum(C(r+m)*(1-C(r)) for r in range(-40,41))
   target=mp.mpf(m)/(mp.exp(b*m)-1)
   need(abs(negative-target)<mp.mpf('1e-85'),'Neutral thermal current norm from original covariance')
   need(abs(positive-negative-m)<mp.mpf('1e-85'),'Unchanged level-one current commutator')
  rows.append({'alpha':str(alpha),'partition':mp.nstr(fermionZ,45),'mean_original_charge':mp.nstr(mean,45),'Jminus1_Jplus1':mp.nstr(1/(mp.exp(b)-1),45)})
 # Independent finite fermionic determinant for actual Weyl histories.
 # This is a convergence diagnostic, not a proof from a finite determinant.
 f={0:.37,1:.5+.2j,3:-.13j,5:.08};g={0:-.21,1:-.11+.33j,2:.17,3:.09j}
 profiles=[f,g];fg={k:f.get(k,0)+g.get(k,0) for k in set(f)|set(g)}
 def current_weyl(profile,alpha):
  zero=mp.mpf(str(profile.get(0,0)))
  weights=[mp.exp(-b*(mp.mpf(q*q)/2+(alpha-mp.mpf('.5'))*q)) for q in range(-12,13)]
  charge=sum(w*mp.exp(1j*zero*q) for w,q in zip(weights,range(-12,13)))/sum(weights)
  variance=sum(mp.mpf(n)*mp.coth(b*n/2)*mp.mpf(str(abs(v)**2)) for n,v in profile.items() if n)
  return complex(charge*mp.exp(-variance/2))
 def matrix(profile,rs):
  A=np.zeros((len(rs),len(rs)),complex)
  for i,r in enumerate(rs):
   for j,t in enumerate(rs):
    k=t-r
    A[i,j]=profile.get(0,0) if k==0 else np.conj(profile.get(k,0)) if k>0 else profile.get(-k,0)
  need(np.max(abs(A-A.conj().T))<1e-14,'Real source profile is Hermitian')
  return A
 cocycle=sum(n*(np.conj(f.get(n,0))*g.get(n,0)-f.get(n,0)*np.conj(g.get(n,0))) for n in set(f)|set(g) if n)
 dets=[]
 for alpha in (.25,.75):
  for cut in (16,32):
   rs=np.arange(-cut,cut);cov=expit(-float(b)*(rs+alpha));I=np.eye(len(rs))
   Af,Ag=[matrix(p,rs) for p in profiles]
   Uf,Ug=expm(1j*Af),expm(1j*Ag)
   def determinant(U,profile0):
    sign,log=np.linalg.slogdet(I-cov[:,None]*I+cov[:,None]*U)
    return sign*np.exp(log-1j*cut*profile0)
   single=determinant(Uf,f[0]);ordered=determinant(Uf@Ug,f[0]+g[0])
   expected=current_weyl(f,mp.mpf(str(alpha)))
   expected_order=np.exp(-cocycle/2)*current_weyl(fg,mp.mpf(str(alpha)))
   need(abs(single-expected)<1e-11,'Weyl functional matches independent fermion determinant')
   need(abs(ordered-expected_order)<1e-11,'Full ordered current source with original central phase')
   need(abs(ordered-current_weyl(fg,mp.mpf(str(alpha))))>1e-3,'Omitting current cocycle is detected')
   dets.append({'alpha':alpha,'mode_cut':cut,'single_error':float(abs(single-expected)),'ordered_error':float(abs(ordered-expected_order)),'ordered_value':[float(ordered.real),float(ordered.imag)]})
 # Independent charged Wick determinants versus the all-word current
 # vertex reconstruction. Global charge theta factor is essential.
 vertex=[]
 for alpha in (mp.mpf(1)/4,mp.mpf(3)/4):
  weights={q:mp.exp(-b*(mp.mpf(q*q)/2+(alpha-mp.mpf('.5'))*q)) for q in range(-12,13)}
  theta=sum(weights.values())
  def kernel(z,w):
   R=z/w
   # Sum the infinite occupied Hardy sea analytically. Only its
   # exponentially decaying thermal correction is numerically cut.
   def delta(r):
    x=mp.exp(b*(r+alpha))
    return -x/(1+x) if r<0 else 1/(1+x)
   return R**alpha*(1/(R-1)+sum(R**r*delta(r) for r in range(-80,81)))
  for pairs in (1,2,4,8):
   zs=[mp.mpf('.5')+mp.mpf(k)/32 for k in range(pairs)]
   ws=[mp.mpf(k+1)/64 for k in range(pairs)]
   R=mp.fprod(zs)/mp.fprod(ws)
   K=mp.matrix([[kernel(z,w) for w in ws] for z in zs])
   C0=mp.matrix([[w/(z-w) for w in ws] for z in zs])
   exponent=sum(sum(z**n-w**n for z,w in zip(zs,ws))*
                sum(w**(-n)-z**(-n) for z,w in zip(zs,ws))/
                (n*(mp.exp(b*n)-1)) for n in range(1,81))
   predicted=R**alpha*mp.det(C0)*sum(w*R**q for q,w in weights.items())/theta*mp.exp(exponent)
   actual=mp.det(K)
   relative=abs((actual-predicted)/predicted)
   need(relative<mp.mpf('1e-55'),f'Full neutral Wick word alpha={alpha},pairs={pairs},relative={mp.nstr(relative,15)}')
   vertex.append({'alpha':str(alpha),'pairs':pairs,'relative_error':mp.nstr(relative,8)})
 n1=1/(mp.exp(b)-1);need(n1>0,'Even Q0 current state differs from Hardy vacuum')
 need(abs(n1-1/(mp.exp(2*mp.pi)-1))>mp.mpf('.001'),'Wrong cover-to-quarter scale is detected')
 result={'verdict':'EXACT_SCOPED_ELLIPTIC_CURRENT_FUNCTIONAL; NEUTRAL_STATE_DIFFERS_FROM_HARDY; ORIGINAL_P1_MEASURE_OPEN','inverse_temperature_quarter':'4*pi','formal_partition_cut':80,'formal_partition_monomials':formal,'quarter_current_rows':rows,'finite_determinant_diagnostics':dets,'neutral_vertex_word_diagnostics':vertex,'neutral_Q0_Jminus1_Jplus1':mp.nstr(n1,45),'source_sha256':{str((SRC/name).relative_to(REPO)):hashlib.sha256((SRC/name).read_bytes()).hexdigest() for name in ['BOUNDARY-SOURCE.txt','ELLIPTIC-GREEN-RP.txt','UPDATE-REAL-SLOTS-AND-SOURCE.txt']},'scope':'Existing same-E free Gaussian quantization. All-profile formula proved separately; finite determinants are diagnostics. No physical temperature,time or primitive P1 measure selected. Old Hardy native source is preserved as a different state.'}
 args.out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
