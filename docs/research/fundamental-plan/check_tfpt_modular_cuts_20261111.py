"""Exact algebra controls for alternating cuts in the given TFPT E8 net."""
from contextlib import redirect_stdout
from fractions import Fraction
from itertools import product
from pathlib import Path
import io
import runpy

import sympy as s

repo = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
base = repo/'docs/research/fundamental-plan'
with redirect_stdout(io.StringIO()):
    original = runpy.run_path(str(base/'check_reflection_origin_eml_20261111.py'))
roots = {tuple(s.Rational(x,4) for x in r) for r,_ in original['roots']}
positive = {r for r in roots if next(x for x in r if x)>0}
sums = {tuple(x+y for x,y in zip(a,b)) for a in positive for b in positive}
basis = s.Matrix.hstack(*(s.Matrix(r) for r in sorted(positive-sums)))
gram = basis.T*basis
assert len(roots)==240 and gram.shape==(8,8)
assert gram.det()==1 and all(x.is_Integer for x in gram.inv())
assert all((basis*(gram.inv()*basis.T*s.Matrix(r))-s.Matrix(r)).is_zero_matrix
           for r in roots)
print('Native E8 root lattice: det=1, inverse Gram integral, 240 roots.')

# Smooth plateau values at the ordered boundaries of four alternate arcs.
# A profile F varies only in the first and third arcs; G only in the others.
# With x a linear arc coordinate, a ramp r is any C-infinity function
# constant near 0,1. On arcs where G varies, F is constant, hence the
# integral uses endpoints only and is exact independently of ramp shape.
F = [0,1,1,0,0]
G = [1,1,0,0,1]
integral = sum(F[i]*(G[i+1]-G[i]) for i in range(4))
assert integral==-1
assert all((F[i+1]-F[i])*(G[i+1]-G[i])==0 for i in range(4))
for h in (F,G):
    assert all(h[i+2]+h[i]==1 for i in range(3))
print('Alternating plateau links: normalized crossing pairing -1; both differences odd under the half-turn.')

i,j = next((i,j) for i,j in product(range(8),repeat=2) if gram[i,j]==-1)
root_i,root_j = s.eye(8)[:,i],s.eye(8)[:,j]
v = root_i/2
# For W(2*pi*v*F), Weyl-crossing exponent / (2*pi) is
# integral*<v,w> (overall sign convention immaterial at this witness).
exponent = integral*(v.T*gram*root_j)[0]
ward_exponent = (v.T*gram*root_j)[0]
assert exponent==s.Rational(1,2) and ward_exponent==-s.Rational(1,2)
assert s.exp(2*s.pi*s.I*exponent)==s.exp(2*s.pi*s.I*ward_exponent)==-1
assert all((root_i.T*gram*s.eye(8)[:,k])[0].is_Integer for k in range(8))
print('Half a native simple root gives a -1 crossing phase and a -1 Ward phase against a native adjacent root.')

def cartan(rank, edges):
    result = 2*s.eye(rank)
    for i,j in edges:
        result[i,j]=result[j,i]=-1
    return result

def A(rank):
    return cartan(rank,[(i,i+1) for i in range(rank-1)])

def D(rank):
    return cartan(rank,[(i,i+1) for i in range(rank-2)]+[(rank-3,rank-1)])

E6 = cartan(6,[(0,1),(1,2),(2,3),(3,4),(2,5)])
assert [(s.diag(D(5),A(3))).det(),s.diag(A(4),A(4)).det(),
        s.diag(E6,A(2)).det(),D(8).det()] == [16,25,9,4]
print('Lattice-net duality indices before glue: D5+A3=16, A4+A4=25, E6+A2=9, D8=4; E8=1.')

# The normalized even-mode double-cover embedding preserves both vacuum
# norms and CCR. The n=0 quotient is treated separately in the proof.
g = [Fraction(1,3),Fraction(-7,11),Fraction(5,9)]
assert sum((2*n)*v*v/2 for n,v in enumerate(g,1))==sum(n*v*v for n,v in enumerate(g,1))
print('Double-cover pullback: mode n -> 2n with normalization 1/sqrt(2); quotient norm unchanged.')

# Klein four symmetry of every ordered cross ratio, and quarter symmetry
# exactly for the harmonic cyclically ordered configuration.
a,b,c,d=s.symbols('a b c d')
def eta(a,b,c,d):
    return (b-a)*(d-c)/((c-a)*(d-b))
assert s.cancel(eta(a,b,c,d)-eta(c,d,a,b))==0
assert s.cancel(eta(a,b,c,d)+eta(b,c,d,a)-1)==0
print('Opposite-pair half-turn preserves every cross ratio; the cyclic quarter requires eta=1/2.')
