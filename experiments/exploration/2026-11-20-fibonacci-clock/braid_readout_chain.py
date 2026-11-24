#!/usr/bin/env python3
"""Exact fusion-tree braid measurement channel versus infinite-history KMS.

Five tau anyons of total tau; left-associated paths (0 vac, 1 tau):
010101, 010111, 011011, 011101, 011111. The off-diagonal
probability of a Fibonacci braid on a recoupled 2x2 pair is 1/phi.
Choose uniformly one of the four elementary braid locations and
measure again in the same tree basis. This CHOICE of measurement
defines a classical bistochastic channel, not an original TFPT source.
"""
from sympy import Matrix, symbols, sqrt, eye, simplify, factor, Rational
phi=(1+sqrt(5))/2
a=1/phi
edges={1:[],2:[(0,3),(1,4)],3:[(2,4)],4:[(0,1),(3,4)]}
L=Matrix.zeros(5)
for es in edges.values():
    for i,j in es:
        L[i,i]+=1; L[j,j]+=1; L[i,j]-=1;L[j,i]-=1
T=eye(5)-(a/4)*L
assert all(simplify(sum(T[i,j] for i in range(5))-1)==0 for j in range(5))
assert all(simplify(sum(T[i,j] for j in range(5))-1)==0 for i in range(5))
assert L.rank()==4
x=symbols("x")
assert factor(L.charpoly(x).as_expr())==x*(x-2)*(x**3-8*x*x+18*x-10)
assert all(simplify(v)==0 for v in T*Matrix([Rational(1,5)]*5)-Matrix([Rational(1,5)]*5))
assert simplify(phi**-2-Rational(2,5))!=0
assert simplify(phi**-1-Rational(3,5))!=0
print("Exact fusion-tree graph",edges)
print("Laplacian charpoly",factor(L.charpoly(x).as_expr()))
print("The unique stationary vector is (1/5)^5; 2/3 totals = 2/5, 3/5.")
print("The distinct infinite-history KMS weights are phi^-2, phi^-1.")
print("NO-GO: braid readout and history KMS are not the same canonical state.")
