"""Compare the six-mark selector with native TFPT operators and test propagation."""
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path

import sympy as s
from sympy import QQ
from sympy.polys.matrices import DomainMatrix

repo = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
original = repo/'_newest2/TFPT_Gesamtdokumentation2_20260927.md'
raw = original.read_text()
code = raw.split('Originalmitglied: `TFPT_Rekursive_Bindung_20260926/audit.py`',1)[1].split('```python\n',1)[1].split('```',1)[0]
expected = '4a19f8145d5ccaef08bfdc5893f64484feae630cf61e49639df52f89f56528c4'
assert expected in [sha256((code.rstrip()+'\n'*n).encode()).hexdigest() for n in range(4)]
native = {'__file__':str(original)}
exec(code.split('# Exact complete 35 x 35 two-module model.',1)[0],native)
G, K = native['Gs'], native['Recovery']
matrix = raw.split(r'W_6=\begin{pmatrix}',1)[1].split(r'\end{pmatrix}',1)[0]
W = s.Matrix([[int(value) for value in line.strip().split('&')]
              for line in matrix.split(r'\\') if line.strip()])
assert W.shape==(5,6) and W.T*G*W==48*s.eye(6)-8*s.ones(6)
assert W*W.T*G==48*s.eye(5)
Kron = s.kronecker_product
I = s.eye(25)
metric = Kron(G,G)
om = s.Matrix(list(G.inv()))
P1 = om*om.T*metric/5
P9 = s.Rational(9,4)*(K-P1)
swap = s.Matrix(25,25,lambda i,j:int(i==5*(j%5)+j//5))
anti = (I-swap)/2
square = s.Matrix.hstack(*(Kron(W[:,q],W[:,q])/40 for q in range(6)))
P = square*(square.T*metric*square).inv()*square.T*metric
assert P.trace()==6 and P*P==P and K*P==P1
assert P==(I+swap)/2-P9
assert (I-K)*(P-P1)==P-P1
print('Native Petz and W6: selector P=P1+P5=(I+swap)/2-P9; original h has energy 1 on P5.')

# Actual logical reflection quotient: each transposition occurs four times,
# with possible overall signs cancelling on the pair.
T = s.zeros(25)
for i,j in combinations(range(6),2):
    indices = list(range(6))
    indices[i],indices[j]=indices[j],indices[i]
    u = W[:,indices]*W.T*G/48
    assert u.T*G*u==G and u*u==s.eye(5)
    T += Kron(u,u)/15
assert T==P1+s.Rational(3,5)*(P-P1)+P9/3+anti/5
t = s.symbols('t')
f = s.interpolate([(s.Integer(1),1),(s.Rational(3,5),1),
                   (s.Rational(1,3),0),(s.Rational(1,5),0)],t)
assert sum((a*T**n for (n,),a in s.Poly(f,t).terms()),s.zeros(25))==P
print('Native 60-reflection quotient resolves P5; six-mark cutoff polynomial:',s.factor(f))

hmark=I-P
H3=Kron(hmark,s.eye(5))+Kron(s.eye(5),hmark)
want={0:6,s.Rational(3,4):19,s.Rational(19,20):5,
      s.Rational(21,20):5,s.Rational(5,4):19,2:71}
assert s.expand(H3.charpoly(t).as_expr()-s.prod((t-x)**m for x,m in want.items()))==0
assert all(H3*Kron(W[:,q],W[:,q],W[:,q])==s.zeros(125,1) for q in range(6))
print('Exact characteristic polynomial on the 125-dimensional original simplex chart:',want)

v=W[:,5]/s.sqrt(40)
w=W[:,0]+W[:,5]/5
w=w/s.sqrt((w.T*G*w)[0])
assert (v.T*G*w)[0]==0 and (w.T*G*w)[0]==1
a,b=Kron(w,v),Kron(v,w)
x=Kron(v,v)
local=w*v.T*G+v*w.T*G
A,B=Kron(local,s.eye(5)),Kron(s.eye(5),local)
comm=lambda a,b:a*b-b*a
assert s.simplify((a.T*metric*P*b)[0])==s.Rational(1,20)
assert s.simplify((x.T*metric*comm(comm(hmark,A),B)*x)[0])==s.Rational(1,10)
single=a*a.T*metric+b*b.T*metric
rem=(I-single)*hmark*a
assert s.simplify((rem.T*metric*rem)[0])==s.Rational(9,200)
print('An edge contributes 1/10 to the original-branch double commutator; single-excitation leakage is 9/200.')

swap23=Kron(s.eye(5),swap)
Htriangle=H3+swap23*Kron(hmark,s.eye(5))*swap23
V3=Kron(v,v,v)
A1=Kron(local,s.eye(5),s.eye(5))
A3=Kron(s.eye(5),s.eye(5),local)
G3=Kron(G,G,G)
assert all(Htriangle*Kron(W[:,q],W[:,q],W[:,q])==s.zeros(125,1) for q in range(6))
assert s.simplify((V3.T*G3*comm(comm(H3,A1),A3)*V3)[0])==0
assert s.simplify((V3.T*G3*comm(comm(Htriangle,A1),A3)*V3)[0])==s.Rational(1,10)
print('Same six ground branches on path and triangle; crossed 1-to-3 response is 0 vs 1/10.')


# The actual original polar frame, q*=5 and original sigma=(234).
q=5
U=(W[:,:5]-W[:,q]*s.ones(1,5))*(s.eye(5)+(1/s.sqrt(6)-1)*s.ones(5)/5)/s.sqrt(48)
assert (U.T*G*U).applyfunc(s.expand)==s.eye(5)
Y=s.diag(s.Rational(1,2),s.Rational(1,2),*([s.Rational(-1,3)]*3))
U2=Kron(U,U)
h=(I-U2.T*metric*K*U2).applyfunc(s.expand)
charges=[int(6*(Y[i,i]-Y[j,j])/5) for i,j in product(range(5),repeat=2)]
neutral=[i for i,q in enumerate(charges) if q==0]
charged=[i for i,q in enumerate(charges) if q!=0]
field=QQ.algebraic_field(s.sqrt(6))
B=h.extract(charged,neutral)
D=h.extract(charged,charged)
B1=DomainMatrix.from_Matrix(B).convert_to(field)
D1=DomainMatrix.from_Matrix(D).convert_to(field)
assert B1.rank()==4 and B1.hstack(D1*B1).rank()==4
leak=(B.T*B).applyfunc(s.expand)
print('Actual marked polar charge: neutral cyclic space 13+4=17; intermediate-return trace',s.expand(s.trace(leak)))
# Adding all charge projectors splits the coupled charged complement.
blocks=[B.extract([i for i,j in enumerate(charged) if charges[j]==q],list(range(13))) for q in (-1,1)]
assert [DomainMatrix.from_Matrix(M).convert_to(field).rank() for M in blocks]==[4,4]
rank4_projectors=[]
for q in (-1,1):
    projector=s.diag(*(int(charges[j]==q) for j in charged))
    rank4_projectors.append(DomainMatrix.from_Matrix(projector*B).convert_to(field))
closure=rank4_projectors[0].hstack(rank4_projectors[1])
assert closure.rank()==8 and closure.hstack(D1*closure).rank()==8
def perm(f):
    return s.Matrix(5,5,lambda i,j:int(i==f(j)))
X=perm(lambda i:(i+1)%5)
G5=perm(lambda i:2*i%5)
Gpair=Kron(G5,G5)
Uclock=Kron(X,X)
V0=s.eye(25)[:,neutral]
assert s.Matrix.hstack(*(Uclock**j*V0 for j in range(5))).rank()==25
print('Add charge queries: dimension 21; add affine shift: dimension 25.')

# Relative parity is affine invariant, but its candidate transfer and charge
# disagree with the source even with the original five slots held fixed.
Yrel=(I-5*swap)/12
Qpair=Kron(Y,s.eye(5))-Kron(s.eye(5),Y)
assert all(comm(Yrel,V)==s.zeros(25) for V in (Uclock,Gpair))
assert s.trace(Qpair)==0 and s.trace(Yrel)==0
assert Qpair.eigenvals()=={s.Rational(-5,6):6,0:13,s.Rational(5,6):6}
assert Yrel.eigenvals()=={s.Rational(-1,3):15,s.Rational(1,2):10}
omega=s.zeros(25,1)
for m in range(5):
    for r,c in [(1,1),(2,s.I),(4,-1),(3,-s.I)]:
        omega[5*((m-r)%5)+(m+r)%5]=c/s.sqrt(20)
assert swap*omega==-omega and Gpair*omega==-s.I*omega
Qrel=s.Rational(4,9)*anti+s.Rational(5,9)*omega*omega.H
Kslot=I-h
assert Qrel*Qrel.H==Qrel**2
assert (Kslot*Qrel).applyfunc(s.expand)==s.zeros(25)
assert s.simplify(s.trace((Kslot-Qrel)**2))==s.Rational(50,9)
print('Relative clock candidate: original and proposed active supports are orthogonal; squared HS distance 50/9.')

positions=list(product(range(5),repeat=2))
affine_actions=[
    tuple(positions.index(((a*i+b)%5,(a*j+b)%5)) for i,j in positions)
    for a,b in product(range(1,5),range(5))
]
pair_swap=tuple(positions.index((j,i)) for i,j in positions)
def centralizer_dimension(generators, require_charge=False):
    unused=set(product(range(25),repeat=2))
    kept=0
    while unused:
        first=unused.pop()
        orbit={first}
        queue=[first]
        for i,j in queue:
            for g in generators:
                point=(g[i],g[j])
                if point not in orbit:
                    orbit.add(point)
                    queue.append(point)
        unused-=orbit
        if not require_charge or all(charges[i]==charges[j] for i,j in orbit):
            kept+=1
    return kept
assert centralizer_dimension(affine_actions)==32
assert centralizer_dimension(affine_actions+[pair_swap])==20
assert centralizer_dimension(affine_actions,require_charge=True)==3
# The weak su(2) alone, conjugated by the same five-cycle, gives elementary
# generators along a connected five-cycle of color/weak slot positions.
edges={tuple(sorted((i,(i+1)%5))) for i in range(5)}
units={(i,j) for i,j in edges}|{(j,i) for i,j in edges}
while True:
    joined=units|{(i,l) for i,j in units for k,l in units if j==k and i!=l}
    if joined==units:
        break
    units=joined
assert len(units)==20
print('Pair commutant dimensions: affine 32; affine+swap 20; affine+fixed U(1) charge 3.')
print('Five translates of the actual weak root generate all 20 off-diagonal sl(5) root spaces.')

# The original Picard differential basis dx/y², x dx/y², dx/y.
from itertools import permutations
from sympy.combinatorics import Permutation
quarter=(-1,s.I,s.I)
omega3=(-1+s.I*s.sqrt(3))/2
deck=(omega3,omega3,omega3**2)
assert s.prod(quarter)==1 and s.simplify(s.prod(deck)-omega3)==0
fixed=s.simplify(s.prod(1-u for u in quarter)*s.prod(1-s.conjugate(u) for u in quarter))
full_pairs=sum(any(v%2 for v in powers) for powers in product(range(4),repeat=2))
assert fixed==16 and full_pairs==12 and fixed*full_pairs/4==48
branch=[(0,1),(1,1),((1+s.I*s.sqrt(3))/2,1),(1,0)]
def det(a,b):
    return a[0]*b[1]-a[1]*b[0]
def cross(a,b,c,d):
    return det(a,c)*det(b,d)/(det(a,d)*det(b,c))
target=cross(*branch)
stabilizer=[p for p in permutations(range(4))
            if s.simplify(cross(*(branch[i] for i in p))-target)==0]
assert len(stabilizer)==12 and all(Permutation(list(p)).signature()==1 for p in stabilizer)
print('Picard original quarter: 16 fixed points, Euler 48; quartic-cover branch group A4 of order 12.')
