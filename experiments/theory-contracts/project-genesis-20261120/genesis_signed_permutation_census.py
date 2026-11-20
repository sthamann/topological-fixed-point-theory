#!/usr/bin/env python3
"""Exhaust the 48 monomial signed-permutation matrices on the 3-family space.
Count the finite extension <V4,G> and its full group-algebra dimension exactly.
No physical inference about why G should exist or have a given energy.
"""
from itertools import permutations,product
from collections import deque,Counter
from sympy import Matrix,eye,diag
D=3
A=diag(-1,1,-1); B=diag(1,-1,-1); V4=[A,B]
key=lambda m:tuple(int(v) for v in m)

def closure(gens):
    seen={key(eye(D)):eye(D)}; todo=deque([eye(D)])
    while todo:
        g=todo.popleft()
        for h in gens:
            k=key(g*h)
            if k not in seen:
                seen[k]=g*h;todo.append(g*h)
    return list(seen.values())

def alg_dim(group):
    return Matrix.hstack(*[Matrix(9,1,list(g)) for g in group]).rank()

rows=[]
for perm in permutations(range(3)):
    for signs in product((-1,1),repeat=3):
        g=Matrix(3,3,lambda i,j: signs[j] if i==perm[j] else 0)
        gg=closure([A,B,g])
        rows.append((perm,signs,len(gg),alg_dim(gg)))
assert len(rows)==48
counts=Counter((order,rank) for _,_,order,rank in rows)
for (size,rank),count in sorted(counts.items()):
    print('generated group order=',size,' algebra dimension=',rank,' cases=',count)
full=[r for r in rows if r[-1]==9]
assert full and min(r[2] for r in full)==12
assert all(r[0] in ((1,2,0),(2,0,1)) for r in full)
print('FULL 3x3 MATRIX ALGEBRA: minimum order',min(r[2] for r in full),
      ' among',len(full),'out of 48 signed-permutation extensions of V4')
print('The extension must permute all three characters; a lone transposition is insufficient')
print('ALL 48 EXACT FINITE GROUP CENSUS CHECKS PASSED')
