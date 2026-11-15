"""Exact weight-one interface and conformal-collar checks on the actual E8 data."""
from itertools import product
from fractions import Fraction as Q

import sympy as s

def dot(a, b):
    return sum(x*y for x,y in zip(a,b))

def add(a, b):
    return tuple(x+y for x,y in zip(a,b))

def flip(a):
    return a[:-1]+(-a[-1],)

def negate(a):
    return tuple(-x for x in a)

clock = (0,2,4,-2,0,4,2,0)
vector = (0,0,0,0,0,0,0,-1)
alpha = (Q(1,2),)*8
beta = tuple(Q(x,2) for x in (1,-1,-1,-1,-1,1,-1,-1))
gamma = add(alpha,beta)
assert dot(clock,vector)==dot(gamma,vector)==0
assert dot(vector,vector)==1
assert add(alpha,vector)==flip(alpha)
assert add(negate(beta),vector)==negate(flip(beta))
assert all(dot(vector,p)==-Q(1,2) for p in (alpha,negate(beta)))
assert dot(clock,alpha)%8==dot(clock,flip(alpha))%8
assert dot(clock,beta)%8==dot(clock,flip(beta))%8

even = {tuple(Q(t,2) for t in v) for v in product((-1,1),repeat=8)
        if v.count(-1)%2==0}
odd = {flip(v) for v in even}
assert len(even)==len(odd)==128 and not even & odd
assert all(dot(p,p)==2 for p in even|odd)
assert {add(p,vector if p[-1]>0 else negate(vector)) for p in even}==odd
assert all(dot(clock,p)%8==dot(clock,flip(p))%8 for p in even)
assert dot(alpha,beta)==-1 and dot(gamma,gamma)==2
print('Crossed marked spinor spectra: 128 versus 128, zero overlap at fixed momentum.')
print('Weight-one transfer: a single unit vector preserves the clock and the middle root; both propagated momenta require v=-e8.')

# In D8*/D8, addition uses the integer and half-integer cosets.
cosets = {"0":(0,)*8,"v":negate(vector),"s":alpha,"c":flip(alpha)}
def lattice_class(p):
    for name, a in cosets.items():
        delta=add(p,negate(a))
        if all(Q(x).denominator==1 for x in delta) and sum(delta)%2==0:
            return name
    raise AssertionError(p)

assert lattice_class(add(cosets["s"],cosets["v"]))=="c"
assert lattice_class(add(cosets["s"],cosets["c"]))=="v"
assert [Q(dot(cosets[name],cosets[name]),2)%1 for name in ("0","v","s","c")]==[0,Q(1,2),0,0]
half_clock = tuple(Q(a,2) for a in clock)
assert lattice_class(add(vector,negate(half_clock)))=="0"
assert dot(clock,add(vector,negate(half_clock)))%8==2
assert (Q(2*4,8)+Q(11*4**2,32))%1==Q(1,2)
assert ((2+2+4)%8,(4+4)%8)==(0,0)
print('D8 fusion s+v=c; vector class 2e+4f has order two and spin 1/2; opposite pairs have integral spin.')

x,y=s.symbols('x y', positive=True)
z1,z2,z3=1,s.I*x,-x*y
W=s.factor(z1*z2*z3/((z1-z2)*(z2-z3)*(z1-z3)))
expected=-s.I*x*y/((1-s.I*x)*(s.I+y)*(1+x*y))
assert s.simplify(W-expected)==0
assert W.subs({x:1,y:1}).simplify()==-s.I/4
energies=[s.simplify((t*s.diff(W,t)).subs({x:1,y:1})) for t in (x,y)]
assert energies==[s.Rational(1,8),s.Rational(1,8)]
assert s.limit(s.limit(W/(x*y),x,0),y,0)==-1
print('Same-completion collar: -i/4 at zero radial separation, L0 insertions 1/8 and 1/8, leading radial coefficient -xy.')
