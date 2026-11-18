"""Exact transport checks for the original E8 clock and Gaussian source."""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, permutations, product

import sympy as s

def dot(a, b):
    return sum(x*y for x,y in zip(a,b))

def add(a, b):
    return tuple(x+y for x,y in zip(a,b))

def negate(a):
    return tuple(-x for x in a)

def axis(index, sign=1):
    return tuple(sign*int(i==index) for i in range(8))

def norm_weight(a):
    return Q(dot(a,a),2)

w = (0,2,4,-2,0,4,2,0)
alpha = (Q(1,2),)*8
minus_beta = tuple(Q(x,2) for x in (-1,1,1,1,1,-1,1,1))
gamma = add(alpha,negate(minus_beta))
delta = tuple(Q(x,2) for x in w)
v = axis(7,-1)
dressing = add(v,negate(delta))
assert sum(delta)%2==1 and sum(dressing)%2==0
assert dot(w,delta)%8==6 and dot(w,dressing)%8==2 and dot(w,v)==0
assert dot(delta,gamma)==2 and dot(dressing,gamma)==-2 and dot(v,gamma)==0
assert [norm_weight(add(p,delta)) for p in (alpha,minus_beta)]==[9,7]
assert [norm_weight(add(p,v)) for p in (alpha,minus_beta)]==[1,1]
assert [dot(w,p)%8 for p in (alpha,minus_beta)]==[5,1]
assert not all((x/4).denominator in (1,2) for x in dressing)
print('Fourth affine step: correct D8 vector coset; conformal weights 9 and 7. Required endpoints have weights 1 and 1.')
print('Dressing has clock grade 2, cancels the clock grade 6 of the half-path and cancels the middle-current mode shift 2.')

# The minimum word length follows from the L1 norm. All minimum paths are
# enumerated using the coordinate signs of the required total shift.
assert sum(abs(x) for x in dressing)==8
candidates = [add(axis(i,-1 if dressing[i]<0 else 1),axis(j,-1 if dressing[j]<0 else 1))
              for i,j in combinations(range(8),2) if dressing[i] and dressing[j]]
count=0
for multiset in combinations_with_replacement(candidates,4):
    if tuple(map(sum,zip(*multiset)))!=dressing:
        continue
    for word in set(permutations(multiset)):
        momenta=[add(p,delta) for p in (alpha,minus_beta)]
        drops=[]
        for root in word:
            drops.append([1+dot(p,root) for p in momenta])
            momenta=[add(p,root) for p in momenta]
        assert momenta==[add(alpha,v),add(minus_beta,v)]
        if all(step<=0 for pair in drops for step in pair):
            count+=1
assert count==468
explicit=(add(axis(5,-1),axis(6,-1)),add(axis(5,-1),axis(7,-1)),
          add(axis(1,-1),axis(2,-1)),add(axis(2,-1),axis(3,1)))
momenta=[add(p,delta) for p in (alpha,minus_beta)]
weights=[tuple(map(norm_weight,momenta))]
for root in explicit:
    momenta=[add(p,root) for p in momenta]
    weights.append(tuple(map(norm_weight,momenta)))
assert weights==[(9,7),(6,5),(5,5),(2,2),(1,1)]
print('Minimal four-root repair: 468 ordered words nonincreasing on both legs. Example weights (9,7),(6,5),(5,5),(2,2),(1,1).')

# A real SU(2) x SU(2) block has dimension ab, or 2ab for quaternionic type.
# The complex trace index is halved in the orthogonal level-one normalization.
def index(d):
    return Q(d*(d*d-1),12)

def admissible_blocks(dimension, level1, level2):
    blocks=[]
    for a,b in product(range(1,dimension+1),repeat=2):
        copies=2 if (a%2==0)!=(b%2==0) else 1
        size=copies*a*b
        k1,k2=copies*b*index(a),copies*a*index(b)
        if size<=dimension and k1<=level1 and k2<=level2:
            blocks.append((a,b,copies,size,k1,k2))

    @lru_cache(None)
    def search(i, size, k1, k2):
        if not size:
            return ((),) if not k1 and not k2 else ()
        if i==len(blocks):
            return ()
        _,_,_,d,l1,l2=blocks[i]
        without=search(i+1,size,k1,k2)
        if d<=size and l1<=k1 and l2<=k2:
            return without+tuple((i,)+tail for tail in search(i,size-d,k1-l1,k2-l2))
        return without

    return [tuple(blocks[j][:3] for j in js) for js in search(0,dimension,level1,level2)]

cases=(admissible_blocks(16,1,39),admissible_blocks(16,3,13))
assert [len(x) for x in cases]==[1,2]
wanted=Counter([phase%8 for k in w for phase in (k,-k)])
assert wanted==Counter({0:6,2:3,4:4,6:3})
for blocks in cases[0]+cases[1]:
    assert all(copies==1 for _,_,copies in blocks) and len(set(blocks))==len(blocks)
    def weights(d):
        return range(1-d,d,2)

    for k,l in product(range(8),repeat=2):
        for signs in product((0,4),repeat=len(blocks)):
            if sum((sgn//4)*a*b for (a,b,_),sgn in zip(blocks,signs))%2:
                continue
            candidate=Counter((k*m+l*n+sgn)%8
                              for (a,b,_),sgn in zip(blocks,signs)
                              for m,n in product(weights(a),weights(b)))
            assert candidate!=wanted
print('SO(16) bilinear level pairs: 1 real block type for (1,39), 2 for (3,13). None admits the actual clock in its normalizer.')

# Native RM(1,3) roots in the published Gaussian source coordinates.
points=tuple(product((0,1),repeat=3))
binary_supports={tuple(j for j,p in enumerate(points) if (dot(linear,p)+offset)%2)
                 for linear in points for offset in (0,1)}
binary_supports={support for support in binary_supports if len(support)==4}
roots=[tuple(2*x for x in axis(i,sign)) for i,sign in product(range(8),(-1,1))]
for support in sorted(binary_supports):
    for signs in product((-1,1),repeat=4):
        root=[0]*8
        for j,value in zip(support,signs):
            root[j]=value
        roots.append(tuple(root))
assert len(set(roots))==240
rays={}
for root in roots:
    z=tuple(s.Integer(a)+s.I*b for a,b in zip(root[::2],root[1::2]))
    key=min(tuple((int(s.re(s.expand(u*x))),int(s.im(s.expand(u*x)))) for x in z) for u in (1,-1,s.I,-s.I))
    rays[key]=z
assert len(rays)==60
vectors=list(rays.values())
counts=Counter()
for a,b in product(vectors,repeat=2):
    overlap=s.expand(sum(x.conjugate()*y for x,y in zip(a,b)))
    counts[int(s.re(overlap))**2+int(s.im(overlap))**2]+=1
moments=[sum(n*Q(k,16)**power for k,n in counts.items())/3600 for power in range(1,5)]
assert moments==[Q(1,4),Q(1,10),Q(1,20),Q(1,32)]
print('Original Gaussian rays: 60; exact frame potentials 1/4,1/10,1/20,1/32. Their POVM spans all 16 source operators.')
