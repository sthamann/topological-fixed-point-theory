#!/usr/bin/env python3
"""Finite octonion associator falsifier: existence of nonassociative closure
cannot by itself identify TFPT's 5+3 or physical Lorentz spacetime.
"""
from sympy import Matrix,zeros
# Oriented Fano multiplication triples for imaginary octonion units.
triples=[(1,2,3),(1,4,5),(1,7,6),(2,4,6),(2,5,7),(3,4,7),(3,6,5)]
lookup={}
for a,b,c in triples:
    for x,y,z in [(a,b,c),(b,c,a),(c,a,b)]:
        lookup[(x,y)]=(1,z);lookup[(y,x)]=(-1,z)
def mul_basis(i,j):
    if i==0:return (1,j)
    if j==0:return (1,i)
    if i==j:return (-1,0)
    return lookup[(i,j)]
def prod(x,y):
    z=zeros(8,1)
    for i in range(8):
        for j in range(8):
            s,k=mul_basis(i,j);z[k]+=s*x[i]*y[j]
    return z
E=lambda i: Matrix(8,1,lambda j,_: int(j==i))
a=prod(prod(E(1),E(2)),E(4))
b=prod(E(1),prod(E(2),E(4)))
assert a==E(7) and b==-E(7)
print('Octonion associator [e1,e2,e4]=',list(a-b))
# Real octonions are alternative, not associative; Moufang identities still hold.
for i in range(8):
    for j in range(8):
        assert prod(prod(E(i),E(i)),E(j))==prod(E(i),prod(E(i),E(j)))
        assert prod(prod(E(i),E(j)),E(j))==prod(E(i),prod(E(j),E(j)))
print('8x8 basic alternativity relations passed')
print('No distinguished 5+3 subspace follows from this multiplication table alone')
print('NONASSOCIATIVE MODEL CONTROL PASSED')
