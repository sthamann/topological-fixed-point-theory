"""Exact tests of one binary reflection law against original TFPT actions."""
from itertools import product, permutations
from collections import Counter
from pathlib import Path
import sys
import sympy as s

repo = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0,str(repo/'verification'))
from v418_cyclotomic_norm_triple import C5
from v117_monodromy_weyl_a3 import M0_EXACT, U_EXACT
from sympy.combinatorics import Permutation

checks=0
def check(test, label):
    global checks
    assert test,label
    checks+=1
    print('PASS',label)

# One operation on the free two-mark integer model. Check each identity
# symbolically; the universal property is proved in the accompanying text.
u,v,w=s.symbols('u v w')
r=lambda x,y: 2*y-x
d=lambda x:r(0,x)
check(r(u,u)==u and s.expand(r(r(u,v),v)-u)==0 and
      s.expand(r(r(u,v),w)-r(r(u,w),r(v,w)))==0,'reflection identities')
check(r(r(u,0),1)==u+2 and d(1)==2 and d(d(1))-r(1,0)==5,
      'translation and doubling from the same rule; anchor closure is a difference of five')
# Independent finite check of the least equivalence compatible with the
# rule, on periodic controls. The exact Z proof is independent of this scan.
def classes(N):
    parent=list(range(N))
    def root(i):
        while i!=parent[i]:
            i=parent[i]
        return i
    def merge(a,b):
        a,b=root(a),root(b)
        if a==b:return False
        parent[b]=a
        return True
    merge(4%N,-1%N)
    again=True
    while again:
        again=False
        for a,b in product(range(N),repeat=2):
            if root(a)!=root(b):continue
            for z in range(N):
                again=merge(r(a,z)%N,r(b,z)%N) or again
                again=merge(r(z,a)%N,r(z,b)%N) or again
    return len({root(i) for i in range(N)})
check(all(classes(N)==s.gcd(N,5) for N in range(2,36)),'anchor quotient checked on 34 cyclic controls')
check([N for N in range(3,200,2) if pow(2,4,N)==1 and pow(2,2,N)!=1]==[5,15],
      'order four alone admits the five and fifteen controls')

def matperm(n,func):
    m=s.zeros(n)
    for j in range(n):m[func(j),j]=1
    return m
C=matperm(5,lambda j:(j+1)%5)
G=matperm(5,lambda j:2*j%5)
R=matperm(5,lambda j:-j%5)
check(G*C*G.inv()==C**2 and G**2==R and G**4==s.eye(5),'affine clock and reflection on the derived quotient')
x=s.symbols('x')
phi=x**4+x**3+x**2+x+1
def col(poly):return s.Matrix([s.Poly(s.rem(poly,phi,x),x).nth(k) for k in range(4)])
B=s.Matrix.hstack(*[col(x**j) for j in range(5)])
nativeG=s.Matrix.hstack(*[col(x**(2*j)) for j in range(4)])
check(B*C==C5*B and B*G==nativeG*B,'exact integral simultaneous v418/v419 clock map')

# F_16 is Phi5 reduced mod 2, not an unrelated four-slot constructor.
def mul(a,b):
    t=0
    while b:
        if b&1:t^=a
        a<<=1
        if a&16:a^=31
        b>>=1
    return t
def pw(a,n):
    t=1
    for _ in range(n):t=mul(t,a)
    return t
points=[pw(2,j) for j in range(5)]
marks={pw(1^z,14) for z in points if z!=1}
trace_one={a for a in range(16) if a^pw(a,4)==1}
check(len(set(points))==5 and len(marks)==4 and marks==trace_one,'the nontrivial carrier phases are the relative trace-one torsor')
check(all(pw(pw(1^z,14),2)==pw(1^pw(z,2),14) for z in points[1:]),'field mark map intertwines actual squaring')
check(all(pw(a,4)==a^1 for a in marks),'second Frobenius is the nontrivial fixed-point-free translation')

# Compare literal pentagon positions and the clock's own 2D representation.
def fracmod(p,q):return s.rem(p*s.invert(q,phi,x),phi,x)
def cross(z):return s.cancel((z[0]-z[2])*(z[1]-z[3])/((z[0]-z[3])*(z[1]-z[2])))
z=[x**j for j in [1,2,4,3]]
n,de=s.fraction(cross(z))
lam=fracmod(n,de)
check(lam.free_symbols=={x} and fracmod(lam**2+lam-1,1)==0,
      'literal cyclotomic quartet has golden rather than harmonic cross-ratio')
G4=matperm(4,lambda j:(j+1)%4)
P=(s.eye(4)-G4**2)/2
a=P[:,0]; b=G4*a
square=[s.simplify((a.dot(P[:,j])+s.I*b.dot(P[:,j]))/a.dot(a)) for j in range(4)]
check(square==[1,s.I,-1,-s.I] and cross(square)==2,'orthogonal clock plane gives the original harmonic square')
# All 24 affine-semilinear maps on the canonical F4 torsor.
F4={z for z in range(16) if pw(z,4)==z}; base=min(marks)
pts=sorted(marks)
aff={tuple(pts.index(base^mul(a,pw(y^base,2**k))^b) for y in pts)
     for a,b,k in product(F4-{0},F4,range(2))}
check(len(aff)==24,'full affine-semilinear action on four marks is S4')
def finite_cross(p):
    a,b,c,d=p
    return mul(mul(a^c,b^d),pw(mul(a^d,b^c),14))
check({finite_cross(q) for q in permutations(pts)}==F4-{0,1},
      'four finite field positions are equianharmonic, not a good reduction of a harmonic square')


# Actual v117 group, including sign-twist, and an exact intertwiner.
def aug(p):
    return s.Matrix([[int(i==p[j])-int(i==p[3]) for j in range(3)] for i in range(3)])
def rep(p):return Permutation(list(p)).signature()*aug(p)
p4=(1,2,3,0)
X=s.Matrix(3,3,s.symbols('v:9'))
found=False
for p3 in permutations(range(4)):
    if Permutation(list(p3)).order()!=3:continue
    lin=s.linear_eq_to_matrix(list(U_EXACT*X-X*rep(p4))+list(M0_EXACT*X-X*rep(p3)),tuple(X))[0]
    if any(s.simplify(s.Matrix(3,3,k).det())!=0 for k in lin.nullspace()):
        found=True
        break
check(found,'sign-twisted augmentation matches actual v117 pair (3-cycle not uniquely marked)')
from itertools import combinations
pairs=list(combinations(range(3),2))
def exterior2(M):return s.Matrix(3,3,lambda i,j:M.extract(pairs[i],pairs[j]).det())
S=s.Matrix([[0,0,1],[0,-1,0],[1,0,0]])
check(all(S*exterior2(aug(p))*S.inv()==Permutation(list(p)).signature()*aug(p).inv().T
          for p in aff),
      'on all S4 actions the oriented area representation is sign-twisted dual')


# Glue counts using scaled integer coordinates (4 times the true vector).
roots=[]
for k in range(4):
    ds=[v for v in product(*[[z for z in range(-5,6) if z%4==2*k%4]]*5)
        if sum((z-2*k)//4 for z in v)%2==0 and sum(z*z for z in v)<=32]
    z0=[3*k,-k,-k,-k]
    es=[v for v in product(*[[z for z in range(-5,6) if z%4==j%4] for j in z0])
        if sum(v)==0 and sum(z*z for z in v)<=32]
    roots.extend((a+b,k) for a,b in product(ds,es) if sum(z*z for z in a+b)==32)
rootset={v for v,k in roots}
check(Counter(k for v,k in roots)=={0:52,1:64,2:60,3:64} and len(rootset)==240,
      'even index-four glue has 240 roots in sectors 52,64,60,64')
carrier_roots=[s.Matrix(v[:5]) for v,k in roots if k==0 and not any(v[5:])]
split=s.diag(1,1,1,-1,-1)
check(len(carrier_roots)==40 and R.trace()==split.trace()==1 and
      sum(R*v==v for v in carrier_roots)==4 and sum(split*v==v for v in carrier_roots)==12,
      'real 3+2 split exists but is not the marked D5 colour/weak split (4 vs 12 fixed roots)')

def on_root(m,v):return tuple(m*s.Matrix(v))
plain=s.diag(G,G4)
check(all(on_root(plain,v) in rootset for v in rootset),
      'native unsigned clock preserves the complete E8 glue on the ordinary A3 lattice')
direct=s.diag(G,-G4)
signed=s.diag(-G,-G4)
check(sum(on_root(direct,v) in rootset for v in rootset)==112,
      'native unsigned carrier clock and v117 sign twist preserve only the even glue sectors')
check(all(on_root(signed,v) in rootset for v in rootset) and
      signed*s.diag(C,s.eye(4))*signed.inv()==s.diag(C**2,s.eye(4)),
      'common sign twist restores the whole E8 action and the affine relation')
Y=s.Matrix(4,4,s.symbols('y:16'))
basis=s.linear_eq_to_matrix(list(Y*C5-C5*Y)+list(Y*nativeG+nativeG*Y),tuple(Y))[0].nullspace()
assert len(basis)==1
vv=basis[0]
vv=vv*s.ilcm(*[s.denom(q) for q in vv])
vv=vv/s.igcd(*[int(q) for q in vv])
W=s.Matrix(4,4,vv)
check(abs(W.det())==25,'signed repair is not a unimodular relabeling of native simultaneous cyclotomic clock')

# Parameter-free equal folding in the even subspace of the finite quotient.
P0=(s.eye(5)+R)/2
P1=C*P0*C.inv()
K=P0*P1*P0
mu=s.symbols('mu')
check(s.factor(K.charpoly(mu).as_expr())==mu**2*(mu-1)*(16*mu**2-12*mu+1)/16,
      'minimal reflection folding has positive nonzero spectrum 1,(3±sqrt(5))/8')
check(K*G!=G*K,'positive two-reflection transfer retains the extra chosen second mirror')
aver=sum((G**j*K*G**(-j) for j in range(4)),s.zeros(5))/4
check(aver*G==G*aver and aver.charpoly(mu).as_expr()==s.expand(mu**2*(mu-1)*(mu-s.Rational(3,8))**2),
      'clock average of this fold restores covariance but loses the split relaxation rates')
check(all(r(j,j)==j for j in range(5)),'diagonal iteration cannot select a state')
# Actual Pin normal-zero longitudinal twists. Work with 4*p/(2*pi)
# to keep the transport, signs and inverse integral.
def doubled_mode(n, quarter):
    lift=2*(4*n+quarter)
    return lift//4,lift%4
rows=[(n,q) for n in range(-100,101) for q in (1,3)]
check(all(doubled_mode(n,q)==(2*n+(q-1)//2,2) for n,q in rows),
      'native 1/4 and 3/4 towers are the two parity channels of NS under doubling')
check(all(doubled_mode(-n-1,4-q)==(-doubled_mode(n,q)[0]-1,2) for n,q in rows),
      'actual momentum reflection interchanges the two parity channels')
check(all(doubled_mode(*doubled_mode(n,q))==(4*n+q,0) for n,q in rows),
      'second doubling lands in odd Ramond modes, not the reflected quarter sector')
check(all(doubled_mode(n+1,q)[0]-doubled_mode(n,q)[0]==2 for n,q in rows),
      'a unit translation becomes a double translation on the combined target')
print(f'{checks} exact checks passed')
