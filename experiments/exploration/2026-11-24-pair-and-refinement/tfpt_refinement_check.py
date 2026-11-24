#!/usr/bin/env python3
"""Exact local edge refinement and native pair selection checks, H6h28. Standalone. No repo modifications."""
import sympy as s
from itertools import combinations
from fractions import Fraction
from collections import defaultdict, Counter
import numpy as np

checks=[]
def check(name,condition):
    assert condition,name
    checks.append(name); print('PASS',name)

# All diagonal tetrads. The homogeneous polynomial is the HR density.
t,e=s.symbols('t e')
d=s.symbols('d0:4'); v=s.symbols('v0:4')

def V(a,b,t):
    bcoef=[-2,s.Rational(1,2),0,s.Rational(1,2),-2]
    return s.expand(sum(bcoef[len(I)]*s.prod(b[i] if i in I else a[i] for i in range(4)) for n in range(5) for I in combinations(range(4),n))+t*s.prod(a[i]-b[i] for i in range(4))/2)
one=[1]*4; f=[1+e*x for x in d]; c=[1+e*x/2+e**2*y for x,y in zip(d,v)]
# expansion modulo e^5
raw=s.Poly(2*V(one,c,t)+2*V(c,f,t)-V(one,f,t),e)
check('exact affine refinement discrepancy begins only at quartic order',all(s.expand(raw.nth(i))==0 for i in range(4)))
q4=s.expand(raw.nth(4)); aff=s.factor(q4.subs(dict.fromkeys(v,0)))
check('affine half-edge additivity fixes t=0 for general diagonal perturbation',s.expand(aff+3*t*s.prod(d)/8)==0)
sol=s.solve([s.diff(q4,x) for x in v],v,dict=True)[0]
print('stationary intermediate second order:',sol)
stat=s.factor(q4.subs(sol));print('stationary two-edge discrepancy:',stat)
# Two decisive directions, conformal vs volume-preserving (linearized).
conf=s.factor(stat.subs(dict(zip(d,[1,1,1,1])))); shear=s.factor(stat.subs(dict(zip(d,[1,-1,0,0]))))
print('conformal:',conf,'shear:',shear)
check('stationary conformal composition requires t=1',s.simplify(conf-s.Rational(3,8)*(1-t))==0)
check('stationary shear discrepancy independent of t and nonzero',shear!=0 and t not in shear.free_symbols)
# general subdivision, in scalar polynomial notation, direct identity
z,l=s.symbols('z l'); aa,bb=s.symbols('a b'); ss=s.symbols('s',positive=True)
# V density at diagonal conformal e=aI and f=bI.
vc=lambda a,b: V([a]*4,[b]*4,t)
mid=(1-ss)*aa+ss*bb
err=s.factor(vc(aa,mid)/ss+vc(mid,bb)/(1-ss)-vc(aa,bb))
print('all affine conformal subdivisions:',err)
check('all affine subdivisions select t=0',s.simplify(err+3*t*ss*(1-ss)*(aa-bb)**4/2)==0)

# Native W exactly as native_source.py (64 internal labels, 60 channels)
ann=[]
for k in range(5):
    a=np.zeros((32,32),dtype=np.int64)
    for m in range(32):
        if (m>>k)&1:a[m^(1<<k),m]=(-1)**((m&((1<<k)-1)).bit_count())
    ann.append(a)
conj=np.eye(32,dtype=np.int64)
for a in ann:conj=conj@(a+a.T)
ix=[m for m in range(32) if m.bit_count()%2==0]
beta=[(conj@a)[np.ix_(ix,ix)] for a in ann+[a.T for a in ann]]
pairs=list(combinations(range(64),2)); fam=list(combinations(range(4),2)); W=[]
for B in beta:
    for a,b in fam:
        W.append({(i,j):int(B[i//4,j//4])*(int(i%4==a and j%4==b)-int(i%4==b and j%4==a)) for i,j in pairs if B[i//4,j//4] and {i%4,j%4}=={a,b}})
check('native sixty rows each contain eight supported pairs',len(W)==60 and all(len(x)==8 for x in W))
check('native W row orthogonality',all(sum(c*y.get(p,0) for p,c in x.items())==8*(i==j) for i,x in enumerate(W) for j,y in enumerate(W)))

def eps(a):
    if len(set(a))<len(a): return 0
    return (-1)**sum(a[i]>a[j] for i in range(len(a)) for j in range(i+1,len(a)))

def add(D,key,c):
    if c and len(set(key))==len(key):
        D[tuple(sorted(key))]+=c*eps(key)

def mul(D,E,out,coef=1):
    for p,a in D.items():
        for q,b in E.items():add(out,p+q,coef*a*b)
# unnormalized symmetric spin tensor: T00=2 sum W psi_i0 psi_j0, T11 likewise, T01 mixed
T=[]
for row in W:
    tup=[]
    for a,b in [(0,0),(0,1),(1,1)]:
        o=defaultdict(int)
        for (i,j),w in row.items():
            add(o,(2*i+a,2*j+b),w); add(o,(2*i+b,2*j+a),w)
        tup.append(dict(o))
    T.append(tup)
# Pair invariant is root-opposite spin vector and the SU4 epsilon on family two-forms.
Q=defaultdict(int)
for A in range(60):
    k,l=divmod(A,6);a,b=fam[l];c,d_=tuple(i for i in range(4) if i not in (a,b))
    bar=6*((k+5)%10)+fam.index((c,d_));eta=eps((a,b,c,d_))
    mul(T[A][0],T[bar][2],Q,eta)
    mul(T[A][2],T[bar][0],Q,eta)
    mul(T[A][1],T[bar][1],Q,-2*eta)
Q={k:v for k,v in Q.items() if v}
print('native holomorphic charge-four monomials:',len(Q),'coefficient census:',dict(sorted(Counter(Q.values()).items())))
check('native self-dual tensor exchange gives nonzero holomorphic quartic',len(Q)>0)
check('every nonzero holomorphic monomial uses four distinct family labels',all(set((i//2)%4 for i in k)==set(range(4)) for k in Q))
check('restricting to any three family labels kills this quartic exactly',all(not any(f not in {(i//2)%4 for i in k} for k in Q) for f in range(4)))
# Family determinant character invariance under swaps; spin SU2 infinitesimal raise.
for f in range(3):
    out=defaultdict(int)
    for key,c in Q.items():
        kk=[]
        for i in key:
            si,fl=divmod(i//2,4);sp=i%2
            kk.append(2*(4*si+(f+1 if fl==f else f if fl==f+1 else fl))+sp)
        add(out,tuple(kk),c)
    check(f'family transposition {f},{f+1} acts by determinant sign',{k:v for k,v in out.items() if v}=={k:-v for k,v in Q.items()})
# SU2 spin raising: replace any spin=1 by 0, check cancellation.
out=defaultdict(int)
for key,c in Q.items():
    for j,i in enumerate(key):
        if i%2:add(out,key[:j]+(i-1,)+key[j+1:],c)
check('Lorentz spinor raising generator annihilates holomorphic quartic',not any(out.values()))
# Spin10 torus charge: each basis Weyl weight (+/- 1/2)^5; invariant support sums zero.
check('all five Spin10 Cartan charges vanish on every monomial',all(all(sum(1-2*((ix[(i//2)//4]>>j)&1) for i in key)==0 for j in range(5)) for key in Q))
# Verify all 45 Spin10 compact Lie generators exactly via original Clifford generators.
gam=[a+a.T for a in ann]+[1j*(a.T-a) for a in ann]
for ia,ib in combinations(range(10),2):
    A=(gam[ia]@gam[ib])[np.ix_(ix,ix)]
    out=defaultdict(complex)
    for key,c in Q.items():
        for j,i in enumerate(key):
            spin,fl=divmod(i//2,4); sp=i%2
            for target in np.flatnonzero(A[:,spin]):
                ni=2*(4*int(target)+fl)+sp
                add(out,key[:j]+(ni,)+key[j+1:],c*A[target,spin])
    check(f'Spin10 Lie generator {ia},{ib} annihilates quartic',not any(out.values()))
print(f'{len(checks)} checks passed. No physical T4 or gravity UV closure claimed.')
# Exact geodesic world function of the diagonal DeWitt edge metric.
# In r=sqrt(prod e), beta_i=log e_i - mean(log e), the metric is
# (3/2)[-dr^2 + (r^2/3) sum d beta_i^2].
from math import cosh, sinh, atanh, exp, log, sqrt

def rg(a):
    u=np.log(np.array(a,float)); return np.exp(u.sum()/2),u-u.mean()
def world(a,b):
    r,ba=rg(a); s,bb=rg(b); ch=np.linalg.norm(bb-ba)/sqrt(3)
    return 1.5*(-(r-s)**2+4*r*s*sinh(ch/2)**2)
def intermediate(a,b,z):
    r,ba=rg(a); s,bb=rg(b); ch=np.linalg.norm(bb-ba)/sqrt(3)
    T=(1-z)*r+z*s*cosh(ch); X=z*s*sinh(ch)
    rm=sqrt(T*T-X*X); f=atanh(X/T)/ch if ch else z
    be=ba+f*(bb-ba)
    return np.exp(be)*sqrt(rm)
rng=np.random.default_rng(20261124); worst=0
for _ in range(80):
    a,b=np.exp(rng.normal(0,.3,(2,4))); z=rng.uniform(.05,.95)
    c=intermediate(a,b,z)
    worst=max(worst,abs(world(a,c)/z+world(c,b)/(1-z)-world(a,b)))
check('80 exact-geodesic interpolation identities for diagonal tetrads',worst<1e-12)
print('max numerical refinement residual:',worst)
# quartic series for two-direction shear, no HR t dependence.
q=s.symbols('q');l1=s.log(1+q).series(q,0,5).removeO();l2=s.log(1-q).series(q,0,5).removeO();sig=l1+l2
chisq=s.series((l1*l1+l2*l2-sig*sig/4)/3,q,0,5).removeO()
coshser=1+chisq/2+chisq**2/24
r=s.sqrt(1-q*q)
gser=s.series(s.Rational(3,2)*(-1-(1-q*q)+2*r*coshser),q,0,5).removeO(); hr=V([1]*4,[1+q,1-q,1,1],t)
print('shear exact-refinement world function:',gser,'HR:',hr)
check('geodesic kernel fixes missing shear quartic coefficient',s.simplify(gser-hr+q**4/36)==0)
# This fixed point is not itself a Hassan-Rosen density: unlike HR,
# it is generally nonlinear in the relative lapse with spatial tetrads fixed.
for b in [1.2,2.0]:
    aa=[1]*4
    vals=[world(aa,[N,b,1,1]) for N in [.7,1,1.3]]
    print('nonlinear lapse witness b=',b,'second difference=',vals[0]-2*vals[1]+vals[2])
    check('nonlinear lapse for spatial scale '+str(b),abs(vals[0]-2*vals[1]+vals[2])>1e-6)
print(f'FINAL: {len(checks)} checks passed.')
# Exact stress spectral measure from the Virasoro vacuum, not a chosen q-state.
from functools import lru_cache
from math import comb
cc=s.Integer(8)
x1,x2,x3,x4=s.symbols('x1 x2 x3 x4')
@lru_cache(None)
def stress(z):
    if not z:return s.Integer(1)
    if len(z)==1:return s.Integer(0)
    u,*vs=z; vs=tuple(vs); lower=stress(vs); ans=0
    for i,w in enumerate(vs):
        ans+=cc/2/(u-w)**4*stress(vs[:i]+vs[i+1:])+2/(u-w)**2*lower+s.diff(lower,w)/(u-w)
    return ans
z=(x1,x2,x3,x4); w={x1:0,x2:1,x3:2,x4:3}
g2=stress(z[:2]); g3=stress(z[:3]); g4=stress(z)
g4c=g4-sum(stress((x1,z[j]))*stress(tuple(z[i] for i in range(1,4) if i!=j)) for j in range(1,4))
print('E8 vacuum separated stress correlators:',g2.subs(w),g3.subs(w),g4.subs(w),'connected4:',s.simplify(g4c.subs(w)))
check('exact stress 3-point at separated points',g3.subs(w)==s.Integer(2))
check('exact non-Gaussian connected stress 4-point',s.simplify(g4c.subs(w))!=0)
# Pairwise spectral probabilities from Virasoro norm c(n^3-n)/12, n=m+2.
check('Virasoro stress weights are exactly NB4',all(cc*s.Integer((m+2)**3-(m+2))/12==cc/2*comb(m+3,3) for m in range(200)))
# Exact positive Jacobi recurrence reconstructed from negative-binomial generating function.
r=s.Rational(1,3); L=14; J=s.zeros(L)
for n in range(L):
    J[n,n]=2+(n+(n+4)*r)/(1-r)
    if n+1<L:
        J[n+1,n]=1;J[n,n+1]=r*(n+1)*(n+4)/(1-r)**2  # monic Jacobi, diagonally similar to symmetric
zv=s.symbols('z');gf=zv**2*((1-r)/(1-r*zv))**4
JJ=s.eye(L);cur=gf
for n in range(20):
    check(f'native stress Jacobi exact moment {n}',JJ[0,0]==s.factor(cur.subs(zv,1)))
    JJ=JJ*J; cur=s.factor(zv*s.diff(cur,zv))
print(f'FINAL: {len(checks)} checks passed.')
