"""Composite Picard differentials, v117 real form, and fixed-charge E8 glue."""
from itertools import permutations, product
from pathlib import Path
import sympy as s
from sympy.combinatorics import Permutation

I=s.I; rt=s.sqrt(3); omega=(-1+I*rt)/2
simple=lambda M:M.applyfunc(s.simplify)
def aug(p):
    return s.Matrix(3,3,lambda i,j:int(i==p[j])-int(i==p[3]))
def sym2(M):
    a,b,c,d=list(M)
    return s.Matrix([[a*a,a*b,b*b],[2*a*c,a*d+b*c,2*b*d],[c*c,c*d,d*d]])
V2s=s.diag(omega,s.conjugate(omega))
V2r=s.Matrix([[-1,-2],[-1,1]])/rt
S=simple(omega*sym2(V2s))
R=simple(-sym2(V2r))
assert simple((-1)*sym2(I*s.eye(2)))==s.eye(3)
pS=(0,2,3,1); pR=(1,0,3,2); pK=(0,1,3,2)
X=s.Matrix(3,3,s.symbols('x:9'))
eq=list(S*X-X*aug(pS))+list(R*X-X*aug(pR))
basis=s.linear_eq_to_matrix(eq,tuple(X))[0].nullspace()
assert len(basis)==1
W=s.Matrix(3,3,basis[0])
fixed=simple(W+W.conjugate()*aug(pK))
if fixed==s.zeros(3):
    fixed=simple(I*W+(I*W).conjugate()*aug(pK))
W=fixed
assert s.simplify(W.det())!=0
assert simple(W.conjugate()-W*aug(pK))==s.zeros(3)
assert simple(S*W-W*aug(pS))==simple(R*W-W*aug(pR))==s.zeros(3)
metric=s.diag(1,1,4)
q=simple(W.H*metric*W)
assert q==s.simplify(q[0,0]/2)*(s.eye(3)+s.ones(3)) and s.simplify(q[0,0])>0
print('Exact invertible Picard composite A4 intertwiner and anti-linear conjugation; positive metric is a scalar times the marked A3 metric.')
print('Intertwiner determinant in the chosen normalization:',s.simplify(W.det()))

# Check the actual v117 matrices against the signed and unsigned 4-cycle.
import sys
root=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0,str(root/'verification'))
from v117_monodromy_weyl_a3 import M0_EXACT as M0,U_EXACT as U
def signed(p):return int(Permutation(list(p)).signature())*aug(p)
b=(1,2,3,0)
assert s.trace(aug(b))==-1 and s.trace(U)==1
for a in permutations(range(4)):
    if Permutation(list(a)).order()!=3:continue
    eq=list(U*X-X*signed(b))+list(M0*X-X*signed(a))
    bs=s.linear_eq_to_matrix(eq,tuple(X))[0].nullspace()
    if bs:
        T=s.Matrix(3,3,bs[0])
        if s.simplify(T.det())!=0:break
assert s.simplify(T.det())!=0
print('Original v117 accepts the signed A3 permutation action. Untwisted linear 4-cycle traces are -1 versus +1.')
print('Picard conjugation is anti-linear; the verified relation is not a linear odd family generator.')

# Both A3 and D5 have discriminant Z/4. A3 spinorial weight differences
# are integral on permutations. An extra minus sends the generator to -1.
f=s.Matrix([3,-1,-1,-1])/4
d=s.ones(5,1)/2
assert (f.T*f)[0]==s.Rational(3,4) and (d.T*d)[0]==s.Rational(5,4)
assert sum(f)==0 and (2*f)[0].q==2 and (2*d).T*(2*d)==s.Matrix([[5]])
for p in permutations(range(4)):
    pf=s.Matrix([f[i] for i in p])
    assert all(t.q==1 for t in pf-f)
assert sum(2*d)==5
Y=(-s.Rational(1,3),)*3+(s.Rational(1,2),)*2
kept,conjugate=0,0
for p in permutations(range(5)):
    for z in product((-1,1),repeat=5):
        yp=tuple(z[j]*Y[p[j]] for j in range(5))
        flips=z.count(-1)
        if yp==Y:
            assert flips==0
            kept+=1
        if yp==tuple(-x for x in Y):
            assert flips==5
            conjugate+=1
assert kept==conjugate==12
print('Glue generator norms 3/4+5/4=2; signed odd A3 action negates its discriminant generator.')
print('D5 isometries preserving fixed diagonal hypercharge:',kept,'all even on sign flips.')
print('D5 signed permutations reversing hypercharge:',conjugate,'all odd on sign flips, exchange the two spinor classes.')
print('Consequently the odd signed v117 family action cannot extend with unchanged fixed D5 hypercharge in the original block-preserving lattice action.')

# Full curve invariance is stricter than the frozen odd v784 polarization.
h1,h2,h3,h4=s.symbols('h:4',real=True)
H=s.Matrix([[h1,h2+I*h3],[h2-I*h3,h4]])
conditions=[]
for M in [V2s,V2r]:
    conditions.extend(list(simple(M.H*H*M-H)))
sol=s.linsolve(conditions,(h1,h2,h3,h4))
assert sol==s.FiniteSet((h4/2,0,0,h4))
print('Full invariant Prym Hermitian forms form one ray, diag(1,2).')
print('The canonical integral Prym polarization is primitive of type (1,2); full-invariant integral scalings are integers and remain mod-2 degenerate.')
print('The odd (1,3) v784 form cannot be invariant under the full 48 automorphisms.')
