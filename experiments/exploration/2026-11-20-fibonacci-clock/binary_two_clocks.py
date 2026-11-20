#!/usr/bin/env python3
"""Complete 2-state primitive binary rule search and distinct modular clocks.

Assumptions: two state labels, all 0/1 local transition matrices,
strong connectivity + aperiodicity, and optional closure under
RELAXING the unique forbidden edge. Physical tensor independence
is a separate axiom, not forced by graph combinatorics.

Standard gauge KMS construction: O_A at beta=ln PF(A);
GNS factor for Fibonacci is III_{1/phi}, full binary is III_{1/2}
(the latter is O2). Assuming a tensor product, incommensurate
logs yield III_1. These facts rely on known operator-algebra
classification, while the finite enumeration is verified below.
"""
from itertools import product
from sympy import Matrix, sqrt, symbols, factor, simplify, log
x=symbols("x")
def primitive(A):
    B=Matrix.eye(2)
    for n in range(1,7):
        B=B*A
        if all(t>0 for t in B):return True
    return False
def label_class(A):
    S=Matrix([[0,1],[1,0]])
    return min(tuple(A),tuple(S*A*S))
C={}
for bits in product([0,1],repeat=4):
    A=Matrix(2,2,bits)
    if primitive(A):C.setdefault(label_class(A),[]).append(A)
assert len(C)==2 and sum(len(v) for v in C.values())==3
F=Matrix([[0,1],[1,1]])
U=Matrix([[1,1],[1,1]])
assert label_class(F) in C and label_class(U) in C
assert factor(F.charpoly(x).as_expr())==x*x-x-1
assert factor(U.charpoly(x).as_expr())==x*(x-2)
assert sum(int(t) for t in F)==3 and sum(int(t) for t in U)==4
assert U-F==Matrix([[1,0],[0,0]])
phi=(1+sqrt(5))/2
p=Matrix([phi**-2,phi**-1])
assert all(simplify(t)==0 for t in F*p-phi*p)
assert simplify(sum(p)-1)==0
f0,f1=0,1
for n in range(1,41):
    assert simplify(phi**n-(f1*phi+f0))==0
    f0,f1=f1,f0+f1
print("EXACT 2-state matrices 16, primitive labeled 3, unlabeled classes 2")
print("unique min-edge primitive: Fibonacci. Unique one-edge completion: full 2 shift.")
print("PF values phi and 2, canonical KMS betas log(phi), log(2).")
print("EXACT irrationality log(phi)/log(2) via phi^n=F_n phi+F_{n-1}.")
print("Conditional tensor product III_{1/phi} ⊗ III_{1/2} is type III_1.")
print("WARNING: product, physical time and spatial local net still not derived.")
