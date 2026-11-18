from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path

import sympy as s
from sympy import QQ
from sympy.polys.matrices import DomainMatrix

repo = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
native = repo / 'docs/research/fundamental-plan/check_tfpt_native_marking_and_propagation_20261111.py'
exec(native.read_text().split('hmark=I-P')[0])
hQ = (I-T)*(3*I/5-T)
hB = I-K
assert hQ == s.Rational(8,45)*P9+s.Rational(8,25)*anti
assert hQ*hB == hB*hQ

def lift(h):
    return Kron(h,s.eye(5)) + Kron(s.eye(5),h)

HQ, HB = lift(hQ), lift(hB)
C = HQ*HB-HB*HQ
assert DomainMatrix.from_Matrix(C).convert_to(QQ).rank() == 48
assert -s.trace(C*C) == s.Rational(704,1125)

def partial_trace(matrix, axis, size=5):
    def index(pair_index, eliminated_index):
        indices = list(divmod(pair_index, size))
        indices.insert(axis, eliminated_index)
        first, middle, last = indices
        return first * size**2 + middle * size + last

    def entry(i, j):
        return sum(matrix[index(i, k), index(j, k)] for k in range(size))
    return s.Matrix(size**2,size**2,entry)

assert all(partial_trace(C,axis)==s.zeros(25) for axis in range(3))
# Real symmetric closed Gram of the two chronological orders.
nested = HQ*C-C*HQ + HB*C-C*HB
assert nested != s.zeros(125)
tau = s.Symbol('tau')
q, b = s.symbols('q b', commutative=False)
def exp_series(operator, scale):
    return 1 + sum((tau*scale*operator)**n/s.factorial(n) for n in range(1,4))

def gram_coefficients(outer, inner):
    a, b = exp_series(outer,-1), exp_series(inner,-2)
    series = s.expand(a*b*a)
    return [series.coeff(tau,n) for n in range(4)]

qb, bq = gram_coefficients(b,q), gram_coefficients(q,b)
c = q*b-b*q
assert qb[:3] == bq[:3]
assert s.expand(qb[3]-bq[3]+(q+b)*c-c*(q+b)) == 0
print('Native three-region commutator: rank 48, norm² 704/1125, all three partial traces zero.')
print('Reflected orders share their quadratic Gram coefficient; cubic difference has squared norm',s.trace(nested*nested))

M = s.zeros(25)
rho = [W[:,j]*W[:,j].T*G/40 for j in range(6)]
assert sum(rho,s.zeros(5)) == 6*s.eye(5)/5
for a in rho:
    M += s.Rational(5,6)*Kron(a,a)
assert M==P1+s.Rational(4,5)*(P-P1)
assert M*K==K*M==P1
assert all(K*s.Matrix(list(a*G.inv())) == s.Matrix(list(G.inv()/5)) for a in rho)
# Symmetric source and channel record alone versus original input distinction.
assert s.trace(sum((s.Rational(5,6)**2*a for a in rho),s.zeros(5)))==s.Rational(25,6)
print('Simplex instrument: Kraus coefficient sqrt(5/6); channel 1 + (4/5) P5; MK=KM=P1.')

def coord(size,i):
    return tuple(int(j==i) for j in range(size))
def dot(a,b):
    return sum(x*y for x,y in zip(a,b))
def norm(a):
    return dot(a,a)

# Exact D5+A3 construction with the fixed spinor-fundamental glue.
roots = []
for i,j in combinations(range(5),2):
    for a,b in product((-1,1),repeat=2):
        roots.append((0,tuple(a*x+b*y for x,y in zip(coord(5,i),coord(5,j)))+(0,)*4))
for i,j in combinations(range(4),2):
    for sign in (-1,1):
        roots.append((0,(0,)*5+tuple(sign*(x-y) for x,y in zip(coord(4,i),coord(4,j)))))
for signs in product((-1,1),repeat=5):
    glue_class=1 if signs.count(-1)%2==0 else 3
    for i in range(4):
        a=tuple(F(x,2) for x in signs)
        b=tuple((1 if glue_class==1 else -1)*(F(x)-F(1,4)) for x in coord(4,i))
        roots.append((glue_class,a+b))
for i,j in combinations(range(4),2):
    b=tuple(F(1,2) if k in (i,j) else F(-1,2) for k in range(4))
    for d,sign in product(range(5),(-1,1)):
        roots.append((2,tuple(sign*x for x in coord(5,d))+b))
assert len(roots)==len({r for _,r in roots})==240 and all(norm(r)==2 for _,r in roots)
assert Counter(c for c,_ in roots)=={0:52,1:64,2:60,3:64}
a=tuple(F(k,8) for k in (3,1,-1,-3))
# Use symmetric eigenangles (same U(5) operator as 0,1/4,1/2,3/4,0).
d=tuple(F(k,4) for k in (0,1,2,-1,0))
for k in range(4):
    shift=tuple(t+F(k,4) for t in d)+a
    data=[dot(shift,r) for _,r in roots]
    order=1
    while any((order*z).denominator!=1 for z in data):
        order+=1
    counts=Counter(int(order*z)%order for z in data)
    counts[0]+=8
    if k%2:
        assert order==4 and (16*norm(shift))%8==4
        assert all((4*z).denominator==1 for z in data)
    else:
        assert order==8 and all((4*dot(shift,r))%1==F(c%2,2) for c,r in roots)
    print('Common quarter scalar i^%d: order %d, grade dimensions %s, squared 4v norm %s' %
          (k,order,[counts[j] for j in range(order)],16*norm(shift)))

# Picard minimum via the ordinary D8/half-integer coordinates of E8.
picard_roots = []
for i,j in combinations(range(8),2):
    for a,b in product((-1,1),repeat=2):
        picard_roots.append(tuple(a*x+b*y for x,y in zip(coord(8,i),coord(8,j))))
picard_roots += [tuple(F(t,2) for t in signs) for signs in product((-1,1),repeat=8) if signs.count(-1)%2==0]
w=(2,1,1,0,0,0,0,0)
assert norm(w)==6
grades=Counter(int(dot(w,r))%4 for r in picard_roots)
grades[0]+=8
assert [grades[i] for i in range(4)]==[82,56,54,56]
print('Picard minimal shift: grade dimensions [82, 56, 54, 56], squared 4V norm 6 mod 8.')

# Square branch divisor of x versus tetrahedral branch divisor of y.
cross_square=(1-(-1))*(s.I-(-s.I))/((1-(-s.I))*(s.I-(-1)))
omega=(-1+s.sqrt(3)*s.I)/2
cross_picard=s.simplify((-1+omega**2)/(-1+omega))
assert s.simplify(cross_square)==2
assert s.simplify(cross_picard**2-cross_picard+1)==0
y=s.symbols('y')
def involution(z):
    return (-z+2)/(z+1)
assert s.cancel(involution(involution(y))-y)==0
assert s.cancel((y+1)*(involution(y)+1))==3
print('Branch cross ratios: square 2, Picard y divisor satisfies lambda²-lambda+1=0; curve lift rho²=1.')

# Invariant Hodge contributions in the torus, and character Euler indices.
h10=(2,1,1)
h01=tuple(-e%4 for e in h10)
h11=sum((a+b)%4==0 for a,b in product(h10,h01))
h21=sum((h10[i]+h10[j]+b)%4==0 for i,j in combinations(range(3),2) for b in h01)
assert (h11,h21)==(5,1)
indices=[]
for k in range(4):
    indices.append(sum((-1)**n for n in range(4)
                       for subset in combinations(h01,n) if (k+sum(subset))%4==0))
assert indices==[0,-2,0,2]
triplets=Counter(sum(indices[k] for k in ks) for ks in product(range(4),repeat=3) if sum(ks)%4==0)
assert triplets=={-4:3,0:10,4:3}
print('Picard invariant Hodge dimensions 5,1; flat equivariant SU(3) Euler indices -4,0,4 with counts 3,10,3.')
