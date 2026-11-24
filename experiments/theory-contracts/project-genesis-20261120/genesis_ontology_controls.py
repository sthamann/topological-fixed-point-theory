#!/usr/bin/env python3
"""Independent toy falsifiers for ontological starting assumptions.
Exact rational/complex arithmetic; not tests of the physical universe.
"""
from itertools import product
from sympy import Matrix, eye, zeros, I, sqrt, Rational, simplify, kronecker_product

# A. Tensor consistency alone does not single out Born's exponent 2.
a=[1,2]; b=[2,3]
for q in [1,2,4,8]:
    p=lambda v:[Rational(abs(x)**q,sum(abs(y)**q for y in v)) for x in v]
    joint=[x*y for x in a for y in b]
    expected=[x*y for x in p(a) for y in p(b)]
    assert p(joint)==expected
    print('Born candidate q=',q,'tensor-product normalized independently: YES')
# The same nontrivial linear 50-50 beamsplitter conserves 2-norm but not other p-norms.
H=Matrix([[1,1],[1,-1]])/sqrt(2)
x=Matrix([1,1]); y=H*x
for q in [1,2,4]:
    original=sum(abs(v)**q for v in x)
    mixed=simplify(sum(abs(v)**q for v in y))
    print('Interferometer q=',q,'norm before/after:',original,mixed)
    assert (original==mixed)==(q==2)

# B. Pairwise invertible observer relations need not be globally gluable.
valid=0
for e_ab,e_bc,e_ca in product([-1,1],repeat=3):
    hol=e_ab*e_bc*e_ca
    if hol==1: valid+=1
assert valid==4
print('Three-view gluing: 8 pairwise valid sign connections, only 4 admit a nonzero globally parallel section')

# C. One-qubit erasure correcting code without any spacetime geometry.
# Four physical qubits; [[4,1,2]] subcode.
n=4; dim=1<<n
def ket(k):
    z=zeros(dim,1);z[k,0]=1;return z
v0=(ket(0)+ket(15))/sqrt(2)
v1=(ket(3)+ket(12))/sqrt(2)
V=Matrix.hstack(v0,v1)
assert V.conjugate().T*V==eye(2)
Id=eye(2);X=Matrix([[0,1],[1,0]]);Y=Matrix([[0,-I],[I,0]]);Z=Matrix([[1,0],[0,-1]])
for i in range(n):
    for P in [X,Y,Z]:
        factors=[P if j==i else Id for j in range(n)]
        E=kronecker_product(*factors)
        result=simplify(V.conjugate().T*E*V)
        assert result==zeros(2)
print('QEC: [[4,1,2]] code detects all 12 single-site Pauli errors, corrects any one known erasure')
XX=kronecker_product(Id,Id,X,X)
assert V.conjugate().T*XX*V==X
print('QEC: a two-site Pauli acts as a nontrivial logical X; distance exactly 2')
print('QEC: no distance metric, physical dimension or Lorentz structure arises from these Knill-Laflamme checks')
print('ALL INDEPENDENT CONTROLS PASSED')
