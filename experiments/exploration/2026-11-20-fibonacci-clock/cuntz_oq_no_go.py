#!/usr/bin/env python3
"""Exact Cuntz-Krieger K-theory and generic type-III_1 KMS null model.

The algebra isomorphism relies on published unital simple CK classification;
the quasi-free factor type on O2 relies on standard KMS factor theorems.
No claim of a native physical TFPT source follows from this test.
"""
from fractions import Fraction
from math import gcd
from sympy import Matrix
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

fib=Matrix([[0,1],[1,1]])
full=Matrix([[1,1],[1,1]])
for name,A in (("Fib",fib),("Full",full)):
    D=Matrix.eye(2)-A.T
    smith=smith_normal_form(D,domain=ZZ)
    assert D.det()==-1
    assert [abs(int(smith[i,i])) for i in range(2)] == [1,1]
    assert all(v>0 for v in A**2)
    print(name,"I-A^T=",D.tolist(),"Smith=",smith.tolist(),"K0=K1=0")

ps=[Fraction(a,b) for b in range(2,21) for a in range(1,b) if gcd(a,b)==1]
assert len(ps)==127 and sum(p==Fraction(1,2) for p in ps)==1
print("Unital simple CK classification: O_Fib isomorphic O_Full isomorphic O2.")
print("Rational quasi-free O2 weights (p,1-p) denominator up to20: 126/127 have")
print("incommensurate logs and, by the quasi-free factor theorem, type III_1.")
print("Proof: p^m=(1-p)^n for positive m,n => reduced coprime a=b-a=1 => p=1/2.")
print("Explicit control p=1/3 and 1-p=2/3 is III_1 without Fibonacci.")
print("WARNING: decorated gauge actions are distinct even if the C*-algebras are isomorphic.")
print("FAIL OF SPECIFICITY: III_1 is not a unique fingerprint of TFPT.")
