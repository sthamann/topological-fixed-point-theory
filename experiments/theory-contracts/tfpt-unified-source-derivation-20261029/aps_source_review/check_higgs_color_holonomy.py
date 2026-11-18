"""R067: original builder and independent Fock diagnostics, not physical Y fits."""
from pathlib import Path
import sys,json
from fractions import Fraction
import numpy as np
import mpmath as mp
from scipy.linalg import expm,eigh,block_diag

repo=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(repo/'verification'))
import v252_full_finite_triple as v
mp.mp.dps=90
b=4*mp.pi
tests=[]
def need(name,error,tol):
    if abs(error)>=tol:raise RuntimeError(f'{name}: {error} >= {tol}')
    tests.append(dict(name=name,error=str(abs(error)),tolerance=str(tol)))

# Reuse the ORIGINAL constructor's seed252 diagnostic inputs, not measured
# or newly optimized Yukawas. The proof covers any fixed given matrices.
rng=np.random.default_rng(252)
ys=[rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)) for _ in range(4)]
mr0=rng.normal(size=(3,3));mr=mr0+mr0.T
zero=np.zeros((3,3),complex)
dy=v.full_D(*ys,zero);dm=v.full_D(zero,zero,zero,zero,mr)
gamma=v.grading()
quark=[i for i,(name,ch,c,g) in enumerate(v.P) if name in ['u','d']]
need('actual particle quark restriction dimension',len(quark)-36,1e-20)
tfull=np.zeros(v.NH)
for i,(name,ch,c,g) in enumerate(v.P):
    if name in ['u','d']:tfull[i]=[1,-1,0][c];tfull[v.NP+i]=-tfull[i]
t=np.diag(tfull[quark]);gam=gamma[np.ix_(quark,quark)]
need('actual quark Tr t2',np.trace(t@t)-24,1e-20)
need('actual chiral plus coefficient',np.trace((np.eye(36)+gam)/2@t@t)-12,1e-20)
need('actual chiral minus coefficient',np.trace((np.eye(36)-gam)/2@t@t)-12,1e-20)
sing=np.r_[np.linalg.svd(ys[0],compute_uv=False),np.linalg.svd(ys[1],compute_uv=False)]

def phi(h):
    h1,h2=h;rr=np.sqrt(abs(h1)**2+abs(h2)**2)
    if rr==0:return dm.copy()
    q=np.array([[np.conj(h2),h1],[-np.conj(h1),h2]],complex)/rr
    unit=np.eye(v.NH,dtype=complex)
    for g in range(3):
        for names,colors in [(('nu','e'),[0]),(('u','d'),range(3))]:
            for c in colors:
                ii=[v.IDX[(name,'L',c,g)] for name in names]
                unit[np.ix_(ii,ii)]=q
                jj=[i+v.NP for i in ii];unit[np.ix_(jj,jj)]=q.conj()
    return unit@(rr*dy+dm)@unit.conj().T

def group36(g):
    ug=np.zeros((36,36),complex)
    for i,oi in enumerate(quark):
        name,ch,c,gen=v.P[oi]
        for j,oj in enumerate(quark):
            name2,ch2,c2,gen2=v.P[oj]
            if (name,ch,gen)==(name2,ch2,gen2):ug[i,j]=g[c,c2]
    return ug

herm=np.array([[.2,.3+.1j,0],[.3-.1j,-.1,.2j],[0,-.2j,-.1]])
gg=expm(1j*herm)
cart=expm(1j*.7*np.diag([1,-1,0]))
for gi,g in enumerate([cart,gg,np.diag([-1,-1,1])]):
    raw=v.rep(1,np.eye(2),g.conj());ks=v.Kswap()
    actual=raw@ks@raw.conj()@ks
    need(f'original real color gauge lift G{gi}',np.linalg.norm(actual[np.ix_(quark,quark)]-group36(g)),1e-13)
cases=[(0,0),(.12,.23j),(.23j,.12),(.4+.1j,-.2+.3j)]
for k,higgs in enumerate(cases):
    ph=phi(higgs);u=sum(abs(z)**2 for z in higgs)
    need(f'original whole Phi odd case{k}',np.linalg.norm(gamma@ph+ph@gamma),2e-12)
    phq=ph[np.ix_(quark,quark)]
    need(f'original color commutes Phi case{k}',np.linalg.norm(t@phq-phq@t),1e-13)
    for p in [.25,-.75]:
        h=p*gam-1j*gam@phq
        ev,vec=eigh(h)
        occ=np.exp(-np.logaddexp(0,float(b)*ev))
        cov=(vec*occ)@vec.conj().T
        ee=np.sqrt(p*p+sing**2*u)
        ss=1/(4*np.cosh(float(b)*ee/2)**2)
        direct_var=np.trace(t@t@cov@(np.eye(36)-cov)).real
        need(f'actual36 forward variance H{k} p{p}',direct_var-4*ss.sum(),3e-13)
        for gi,g in enumerate([cart,gg]):
            ug=group36(g)
            direct=np.linalg.det(np.eye(36)-cov+cov@ug)
            expected=np.prod([np.linalg.det(np.eye(3)-s*(2*np.eye(3)-g-g.conj().T)) for s in ss])
            need(f'actual36 full SU3 determinant H{k} p{p} G{gi}',direct-expected,2e-13)

# Independent six-CAR/64-state Fock computation for a single flavor,
# with nondiagonal SU3 action, no one-particle determinant evaluator.
def annihilator(i,m=6):
    z=np.zeros((2**m,2**m),complex)
    for state in range(2**m):
        if state>>i&1:z[state^(1<<i),state]=(-1)**((state&((1<<i)-1)).bit_count())
    return z
aa=[annihilator(i) for i in range(6)]
def second(mat):
    out=np.zeros((64,64),complex)
    for i in range(6):
        for j in range(6):out+=mat[i,j]*(aa[i].conj().T@aa[j])
    return out
qop=second(block_diag(np.diag([1,-1,0]),np.diag([1,-1,0])))
ufock=expm(1j*second(block_diag(herm,herm)))
for mass in [0,.15,.6]:
    p=.25
    h=np.block([[p*np.eye(3),-1j*mass*np.eye(3)],[1j*mass*np.eye(3),-p*np.eye(3)]])
    val,vec=eigh(second(h));ww=np.exp(-float(b)*(val-val.min()))
    rho=(vec*ww)@vec.conj().T/ww.sum()
    ss=1/(4*np.cosh(float(b)*np.sqrt(p*p+mass*mass)/2)**2)
    expect=np.linalg.det(np.eye(3)-ss*(2*np.eye(3)-gg-gg.conj().T))
    need(f'independent64 Fock group trace mass{mass}',np.trace(rho@ufock)-expect,4e-13)
    need(f'independent64 zero color mean mass{mass}',np.trace(rho@qop),3e-13)
    need(f'independent64 color variance mass{mass}',np.trace(rho@qop@qop)-4*ss,4e-13)

# An exact rational Fock trace at a nontrivial group element.
f=Fraction(1,5);rational=Fraction()
for state in range(64):
    prob=Fraction(1);charge=0
    for i in range(6):
        occ=f if i<3 else 1-f
        bit=state>>i&1
        prob*=occ if bit else 1-occ
        charge+=[1,-1,0,1,-1,0][i]*bit
    rational+=prob*((-1)**abs(charge))
if rational!=Fraction(81,625):raise RuntimeError('rational sixCAR trace mismatch')
tests.append(dict(name='exact rational sixCAR trace',value=str(rational)))

# Whole original tower. These u values are constructor diagnostics only.
def sfun(e):return 1/(4*mp.cosh(b*e/2)**2)
def variance(a,u):
    return 4*sum(sfun(mp.sqrt((mp.mpf(n)+a)**2+mp.mpf(str(y*y))*u)) for n in range(-24,25) for y in sing)
def xi(a,u,theta):
    return mp.exp(2*sum(mp.log1p(-4*sfun(mp.sqrt((mp.mpf(n)+a)**2+mp.mpf(str(y*y))*u))*mp.sin(theta/2)**2)
                       for n in range(-24,25) for y in sing))
a1,a2=mp.mpf(1)/4,mp.mpf(3)/4
for a in [a1,a2]:
    need(f'free original tower variance a{a}',variance(a,0)-3/mp.pi,mp.mpf('1e-80'))
theta=mp.mpf('.7');vprobe=theta/(2*mp.pi)
norm=(mp.jtheta(3,mp.pi*vprobe,mp.exp(-mp.pi))/mp.jtheta(3,0,mp.exp(-mp.pi)))**2
need('whole free pair Quillen power12',xi(a1,0,theta)*xi(a2,0,theta)-norm**12,mp.mpf('1e-80'))
need('exact free quarter color parity1/8',xi(a1,0,mp.pi)-mp.mpf(1)/8,mp.mpf('1e-80'))
need('exact free pair color parity1/64',xi(a1,0,mp.pi)*xi(a2,0,mp.pi)-mp.mpf(1)/64,mp.mpf('1e-80'))
values=[]
for u in [mp.mpf(0),mp.mpf('.01'),mp.mpf('.1'),mp.mpf(1)]:
    values.append(dict(u=str(u),variance_pair=str(variance(a1,u)+variance(a2,u)),xi_pair=str(xi(a1,u,theta)*xi(a2,u,theta))))
for i in range(1,len(values)):
    if not mp.mpf(values[i]['variance_pair'])<mp.mpf(values[i-1]['variance_pair']):raise RuntimeError('variance monotonicity diagnostic')
    if not mp.mpf(values[i]['xi_pair'])>mp.mpf(values[i-1]['xi_pair']):raise RuntimeError('bounded group response monotonicity diagnostic')
    tests.append(dict(name=f'whole tower mixed source ordering {i}',passed=True))

out=Path(__file__).with_name('HIGGS-COLOR-HOLONOMY.json')
result=dict(verdict='EXACT_SCOPED_JOINT_HIGGS_COLOR_RESPONSE',controls=len(tests),tests=tests,diagnostic_values=values,
            original_quark_dimension=36,quark_t_squared_trace=24,chiral_coefficients=[12,12],
            diagnostic_inputs='Unchanged seed252 constructor controls, not physical Yukawa fits or H weights',
            scope='Given constant original conditional quarter CAR source; no primitive measure, local4D or EM identification')
out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(dict(verdict=result['verdict'],controls=len(tests),output=str(out))))
