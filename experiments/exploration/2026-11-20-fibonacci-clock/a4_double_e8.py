#!/usr/bin/env python3
"""Two A4 root systems + an EXTRA even unimodular Z5 glue yield E8.

This is an algebraic route from the root lattice of SU5; it does NOT mean
Fibonacci anyons supply a second physical SU5, a positive source lattice,
or a chiral E8_1 VOA. Invertible E8 stacking remains invisible to MTC data.
"""
from fractions import Fraction
from math import comb
from sympy import Matrix

def q(k):
    k%=5
    return Fraction(k*(5-k),5)
valid=[]
for multiplier in (1,2,3,4):
    norms=[q(k)+q(multiplier*k) for k in range(5)]
    good=all(x.denominator==1 and x.numerator%2==0 for x in norms)
    print("Z5 glue multiplier",multiplier,"even:",good,"norms:",list(map(str,norms)))
    if not good:continue
    valid.append(multiplier)
    # minimal vectors: 20 in each base A4, plus C(5,k)*C(5,mk) in each coset
    roots=40+sum(comb(5,k)*comb(5,(multiplier*k)%5)
                 for k in (1,2,3,4) if norms[k]==2)
    assert roots==240
    print(" root census: 20+20+50+50+50+50 =",roots)
assert valid==[2,3]
# Explicit integer basis for m=2 gluing:
# [alpha1 alpha2 alpha3 beta1 beta2 beta3 beta4 g=omega1+omega2].
# 5g=4alpha1+3alpha2+2alpha3+alpha4+3beta1+6beta2+4beta3+2beta4.
A=Matrix(4,4,lambda i,j:2 if i==j else (-1 if abs(i-j)==1 else 0))
G=Matrix.zeros(8)
for i in range(3):
    for j in range(3):G[i,j]=A[i,j]
for i in range(4):
    for j in range(4):G[i+3,j+3]=A[i,j]
G[7,7]=2; G[0,7]=G[7,0]=1; G[4,7]=G[7,4]=1
minors=[int(G[:k,:k].det()) for k in range(1,9)]
assert minors==[2,3,4,8,12,16,20,1]
assert G.det()==1 and all(G[i,i]%2==0 for i in range(8))
print("E8 Gram det=1, Sylvester minors=",minors)
print("EXACT: even positive unimodular rank 8; unique lattice E8.")
print("OPEN: physical doubling, lattice charge quantization, level 1, c3 and P1.")
