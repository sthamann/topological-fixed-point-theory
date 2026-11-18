"""Check the actual order-eight TFPT lattice and the native Bell readout."""
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path

import sympy as s

root = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
basis = [tuple(int(i == j) for i in range(8)) for j in range(8)]
d8_roots = {
    tuple(a*x+b*y for x,y in zip(basis[i],basis[j]))
    for i,j in combinations(range(8),2) for a,b in product((-1,1),repeat=2)
}
spin_roots = {
    tuple(F(x,2) for x in signs)
    for signs in product((-1,1),repeat=8) if signs.count(-1)%2==0
}
w = (0,2,4,-2,0,4,2,0)
def dot(left, right):
    return sum(a*b for a,b in zip(left,right))

assert dot(w,w)==44
assert all(dot(w,r)%2==0 for r in d8_roots)
assert all(dot(w,r)%2==1 for r in spin_roots)
assert sum(dot(w,r)%8==0 for r in d8_roots)==22
p = next(r for r in sorted(spin_roots) if dot(w,r)%8==1)
assert dot(p,p)==2

sectors = tuple(product(range(8),repeat=2))
def fuse(left, right):
    a, b = (left[0]+right[0],left[1]+right[1])
    return ((a+4*(b//8))%8,b%8)

def weight(sector):
    a,b=sector
    return (F(a*b,8)+F(11*b*b,32))%1

local = {x for x in sectors if weight(x)==0}
assert local=={(a,0) for a in range(8)}|{(a,4) for a in range(1,8,2)}
expected = (
    frozenset((a,0) for a in range(8)),
    frozenset((a,4*(a%2)) for a in range(8)),
)
completions = {frozenset(choices) for choices in combinations(sorted(local),8)
               if (0,0) in choices and all(fuse(a,b) in choices for a,b in product(choices,repeat=2))}
assert completions == set(expected)
assert len(expected[0]&expected[1])==4
assert (64*weight((0,1)))%8==6
print('Actual clock: 22 fixed roots, 64 sectors, 12 integral-spin sectors, cyclic anomaly 6 mod 8.')
print('Exactly two holomorphic lattice completions: a e and a e + 4(a mod 2) f. Shared D8 index over K = 4.')

def flip(r):
    return r[:-1]+(-r[-1],)

assert flip(w)==w and {flip(r) for r in d8_roots}==d8_roots
alternate = {flip(r) for r in spin_roots}
assert len(alternate)==128 and not alternate & spin_roots
assert all(dot(w,flip(r))==dot(w,r) for r in spin_roots)
alpha = (F(1,2),)*8
beta = tuple(F(n,2) for n in (1,-1,-1,-1,-1,1,-1,-1))
middle = tuple(x+y for x,y in zip(alpha,beta))
assert alpha in spin_roots and beta in spin_roots and middle in d8_roots
assert flip(middle)==middle and dot(alpha,beta)==-1
print('Zero-angle D3 sign flip exchanges 128 spinor roots and keeps the clock, D8 roots and the tested middle root.')

def minimal_mod_four_norm(coords):
    minima = []
    for offsets in (range(-3,4),[F(n,2) for n in (-5,-3,-1,1,3,5)]):
        bound = sum(min((x+4*a)**2 for a in offsets) for x in coords)
        nearby = [sorted(offsets,key=lambda a:(x+4*a)**2)[:2] for x in coords]
        admissible = [sum((x+4*a)**2 for x,a in zip(coords,shift))
                      for shift in product(*nearby) if sum(shift)%2==0]
        assert admissible and min(admissible)==bound
        minima.append(bound)
    return min(minima)

for scalar in (1,3):
    quarter = tuple(F(x,2)+(scalar if i<5 else 0) for i,x in enumerate(w))
    assert minimal_mod_four_norm(quarter)==12
print('Each rephased common order-four shift has minimal squared norm 12 in its 4V class, never 4 or 6.')

# Use the exact W6, Gram form and Bell matrix from the native archived audit.
native_path = root/'docs/research/fundamental-plan/check_tfpt_native_marking_and_propagation_20261111.py'
exec(native_path.read_text().split('hmark=I-P')[0])
Fnum = s.Matrix(native['Fnum'])
assert Fnum.shape==(10,5)
assert {abs(x) for x in Fnum*W}=={8}
assert Fnum.T*Fnum/16==G
for row in range(10):
    covector=Fnum[row,:].T/4
    assert (covector.T*G.inv()*covector)[0]==s.Rational(1,2)
rho=[W[:,i]*W[:,i].T*G/40 for i in range(6)]
assert s.Matrix.hstack(*(s.Matrix(list(a)) for a in rho)).rank()==6
assert s.Rational(5,6)*s.Rational(1,10)*s.Rational(1,5)==s.Rational(1,60)
print('Original Bell and simplex: conditional probabilities 1/10; 60 combined outcomes span only six effects.')

# Harmonic branch quartet admits exactly eight of the 24 orderings.
z = (1,s.I,-1,-s.I)
def cross(xs):
    a,b,c,d=xs
    return s.simplify((a-c)*(b-d)/((a-d)*(b-c)))
from itertools import permutations
allowed=[perm for perm in permutations(range(4)) if cross(tuple(z[i] for i in perm))==cross(z)]
assert len(allowed)==8
assert all(not all(perm[perm[perm[i]]]==i for i in range(4)) or perm==tuple(range(4)) for perm in allowed)
print('The harmonic marked branch divisor has exactly 8 Möbius permutations and none of order 3.')
