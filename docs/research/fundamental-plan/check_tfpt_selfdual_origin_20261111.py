"""Original E8 comparison and geometry of the original Q4 orientation."""
from pathlib import Path
import runpy
from contextlib import redirect_stdout
import io
from itertools import combinations, permutations
from collections import Counter
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy import QQ
field=QQ.algebraic_field(s.sqrt(5),s.sqrt(6))

root=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
base=root/'docs/research/fundamental-plan'
with redirect_stdout(io.StringIO()):
    old=runpy.run_path(str(base/'check_reflection_origin_eml_20261111.py'))
R=[s.Matrix(r)/4 for r,_ in old['roots']]
Cnat=s.diag(old['C'],s.eye(4)); Gnat=s.diag(old['G'],old['G4'])
fixed=[r for r in R if Cnat*r==r]
F=s.Matrix.hstack(*fixed)
moving=[r for r in R if F.T*r==s.zeros(len(fixed),1)]
E=s.eye(5); A=s.Matrix.hstack(*(E[:,i]-E[:,i+1] for i in range(4)))
def basis(roots):
    # Standard lexicographic chamber.
    pos={tuple(r) for r in roots if next(x for x in r if x)>0}
    sums={tuple(x+y for x,y in zip(a,b)) for a in pos for b in pos}
    simple=[s.Matrix(r) for r in sorted(pos-sums)]
    for seq in permutations(simple):
        B=s.Matrix.hstack(*seq)
        if B.T*B==A.T*A:return B
    raise ValueError('not an A4 root system')
E9=s.eye(9)
B1=s.Matrix.hstack(*(E9[:,i]-E9[:,i+1] for i in range(4)))
assert all(B1[:,i] in moving for i in range(4))
B2=basis(fixed)
X=A*(A.T*A).inv()
def to_pair(r):
    return X*B1.T*r,X*B2.T*r
def to_native(a,b):
    return B1*X.T*a+B2*X.T*b
pairs=[to_pair(r) for r in R]
def grade(x):
    k=5*x[4]
    assert k.is_Integer
    return int(k)%5
a_set={grade(y)*pow(grade(x),-1,5)%5 for x,y in pairs if grade(x)}
assert a_set in ({2},{3})
if a_set=={3}:
    B2=-B2
    pairs=[to_pair(r) for r in R]
assert all(grade(y)==2*grade(x)%5 for x,y in pairs)
pairset={(tuple(x),tuple(y)) for x,y in pairs}
assert len(pairset)==240
print('Original E8: explicit A4 x A4 coordinates, all 240 roots in the glue graph (k,2k).')
print('Root counts per glue grade:',dict(sorted(Counter(grade(x) for x,y in pairs).items())))
def pm(seq):
    return s.Matrix(5,5,lambda i,j:int(i==seq[j]))
C=pm([(i+1)%5 for i in range(5)]); G=pm([2*i%5 for i in range(5)])
def Jside(x,y):return -y,x
def gear(x,y):return -G*y,G*x
def shift(x,y):return C*x,C*y
def shift_x(x,y):return C*x,y
def shift_y(x,y):return x,C*y
def compose(f,g,r):return f(*g(*r))
for x,y in pairs:
    j=Jside(x,y); h=gear(x,y)
    assert (tuple(j[0]),tuple(j[1])) in pairset and (tuple(h[0]),tuple(h[1])) in pairset
    assert grade(j[0])==3*grade(x)%5 and grade(h[0])==3*grade(x)%5
    assert compose(Jside,Jside,(x,y))==(-x,-y)
    assert compose(gear,shift,(x,y))==shift(*shift(*gear(x,y)))
    assert compose(gear,shift_x,(x,y))==shift_y(*shift_y(*gear(x,y)))
    assert compose(gear,shift_y,(x,y))==shift_x(*shift_x(*gear(x,y)))
    assert gear(*gear(*gear(*gear(x,y))))==(x,y)
assert Cnat.charpoly().as_expr().factor()==(s.Symbol('lambda')-1)**5*s.cyclotomic_poly(5,s.Symbol('lambda'))
print('A common exchanging lattice lift exists: H=(-Gy,Gx), H^4=1, H(C,C)=(C,C)^2 H.')
print('The exchanging lift sends actual one-side C5 to the square of the other side clock, and conversely. Two crossings invert the original C5.')
print('Diagonal C acts on both sides, whereas original C5 acts only on one side.')
for r,(x,y) in zip(R,pairs):
    assert grade(to_pair(Gnat*r)[0])==grade(x)
    assert to_pair(Cnat*r)==shift_x(x,y)
print('Original G preserves the actual glue grades. Any exchanging map would move the fixed A4 of original C5, so cannot conjugate C5 to C5^2.')
# In the moving A4 choose A3 by omitting the fixed slot. Its actual root
# complement has the D5 certificate. This tests the supplied construction
# in the original E8 coordinates, not only an isomorphic root system.
a3=[to_native(E[:,i]-E[:,j],s.zeros(5,1)) for i in range(1,5) for j in range(1,5) if i!=j]
D=[r for r in R if all(r.dot(v)==0 for v in a3)]
nD=s.Matrix.hstack(*D)
assert len(a3)==12 and len(D)==40 and nD.rank()==5
assert all(Cnat*r==r for r in R if all(v==0 for v in r[:5]))
assert any(Cnat*r!=r for r in a3)
print('A3 picked in the moving A4 has 40 complementary rank-five roots, but is distinct from native fixed family A3.')

# Exact native orientation and its state cone.
Is=s.eye(5); Ps=s.ones(5)/5; P=Is-Ps; L=C-C.T+s.I*(C**2-(C.T)**2)
U=Ps+L/s.sqrt(5); K=G**2
sm=lambda M:M.applyfunc(s.simplify)
assert L*L.H==5*P and L.T==-L and G*L*G.T==-s.I*L
assert sm(U*U.conjugate())==Ps-P and sm(U*U.H)==Is
assert K*L*K.T==-L
sm=lambda M:M.applyfunc(s.simplify)
def trf(M):return sm(U*M.conjugate()*U.H)
def project(M):
    N=P*M*P
    E=(N+trf(N))/2
    return sm((E+K*E*K)/2)
assert K.trace()==1 and (K*P).trace()==0
assert K*U*P==-U*K*P
print('Native Q4 orientation yields Herm(2,H); the native half-turn fixes Herm(2,C), dimension 4 including trace.')
with redirect_stdout(io.StringIO()):
    q4=runpy.run_path(str(base/'check_tfpt_q4_picard_selection_20261111.py'))
S=q4['sv']
Q=[sm(S[:,j]*S[:,j].T) for j in range(10)]
assert len(Q)==10
Vplus=s.Matrix([[-2,-2],[1,0],[0,1],[0,1],[1,0]])
Gp=Vplus.T*Vplus
def chart_form(X):
    N=Vplus.T*X*Vplus+Vplus.T*U*X.conjugate()*U.H*Vplus
    return (N/2).applyfunc(s.expand)
bb=[]
for xi in Q:
    f=chart_form(xi)
    bb.append(s.Matrix([f[0,0],s.re(f[0,1]),s.im(f[0,1]),f[1,1]]).applyfunc(s.expand))
assert DomainMatrix.from_Matrix(s.Matrix.hstack(*bb)).convert_to(field).rank()==4
print('All four cone directions are spanned by projected native Q4 effects.')
qv=(U*S[:,0]).applyfunc(s.expand)
imaginary=(s.im(qv*qv.H)).applyfunc(s.expand)
imag_norm=s.simplify(s.trace(imaginary.T*imaginary))
assert imag_norm>0
print('The proposed antiunitary takes an original real Q4 effect outside the real response image. Imaginary squared norm:',imag_norm)

# Closed Q4 rate comes from the trace of the actual rank-one instrument.
assert s.simplify(sum(s.trace(q*q) for q in Q)/2)==5
assert q4['K'].trace()==5
print('Original Q4 supertrace equals 5 from ten unit projectors with weights 1/2. With fixed 1+9 support, this gives 4/9.')

# The real chart at each anchor is transported by the original C5.
# Check whether a fixed chart is preserved and what survives all anchors.
assert C*K*C.T*K==C**2
cyc=[P*(C+C.T),P*(C**2+(C.T)**2)]
assert s.Matrix.hstack(*(s.Matrix(list(x)) for x in cyc)).rank()==2
assert all(trf(X)==X and K*X*K==X for X in cyc)
assert G*cyc[0]*G.T==cyc[1] and G*cyc[1]*G.T==cyc[0]
assert cyc[0]+cyc[1]==-P
print('Two adjacent and all five anchored charts intersect in dimension two. Full affine invariance leaves only centered identity.')

# Six affine charts on the same five labels, two exact orientations.
def ray_key(M):
    a=next(x for x in M if x)
    return tuple(sm(M/a))
left={}; right={}
for perm in permutations(range(5)):
    Pp=pm(perm)
    Dm=Pp*L*Pp.T
    left[ray_key(Dm)]=Dm
    right[ray_key(Dm.conjugate())]=Dm.conjugate()
assert len(left)==len(right)==6 and not set(left)&set(right)
l=list(left.values()); r=[x.conjugate() for x in l]
for i,x in enumerate(l):
    for j,y in enumerate(l):
        assert s.simplify(abs(s.trace(x.H*y)/20)**2)==int(i==j)
        assert s.simplify(abs(s.trace(x.H*r[j])/20)**2)==(0 if i==j else s.Rational(1,5))
print('Exact six-clock census: each orientation is an orthonormal basis, opposite orientation overlap 0 or 1/5.')

# Original supported Weyl witness R059. These are the exact finite code
# blocks; the full source Weyl operators commute by disjoint support.
MA=s.diag(s.Rational(19,32),s.Rational(1,4),s.Rational(1,4),s.Rational(1,4),0)
MB=s.Matrix([[s.Rational(5,24),s.Rational(5,48),0,0,0],
             [s.Rational(5,4),s.Rational(5,8),0,0,0],
             [0,0,-s.Rational(1,2),s.Rational(1,2),0],
             [0,0,s.Rational(1,2),-s.Rational(1,2),0],
             [0,0,0,0,0]])
# Polynomial coefficients f_i -> symmetric tensor orbit coefficients,
# m_i=1,6,6,6,24. This isometry preserves the actual R059 and Q4 Grams.
Dq=s.diag(1,s.Rational(1,6),s.Rational(1,6),s.Rational(1,6),s.Rational(1,24))*s.sqrt(24)
assert Dq.T*q4['G']*Dq==s.diag(96,8,8,8,1)
Rslots=(q4['W'].T*q4['G']*Dq).applyfunc(s.expand)
Anew=sm(Gp.inv()*chart_form(sm(Rslots*MA*s.diag(96,8,8,8,1).inv()*Rslots.T)))
Bnew=sm(Gp.inv()*chart_form(sm(Rslots*MB*s.diag(96,8,8,8,1).inv()*Rslots.T)))
lie=sm(Anew*Bnew-Bnew*Anew)
print('Cone compression of R059 disjoint local Weyls: commutator squared norm',s.simplify(2*s.trace(Gp.inv()*lie.H*Gp*lie)))

# Root norm splitting for the proposed alternate family mark.
aBasis=s.Matrix.hstack(*(to_native(E[:,i]-E[:,i+1],s.zeros(5,1)) for i in range(1,4)))
pA=aBasis*(aBasis.T*aBasis).inv()*aBasis.T
counts=Counter((pA*r).dot(pA*r) for r in R)
assert counts=={0:40,s.Integer(2):12,s.Integer(1):60,s.Rational(3,4):128}
dpos={tuple(r) for r in D if next(x for x in r if x)>0}
dsums={tuple(x+y for x,y in zip(a,b)) for a in dpos for b in dpos}
DB=s.Matrix.hstack(*(s.Matrix(r) for r in sorted(dpos-dsums)))
assert DB.cols==5 and (DB.T*DB).det()==4
print('Alternate A3 complement certificate: determinant 4 and exact root partition 40,12,60,128.')
