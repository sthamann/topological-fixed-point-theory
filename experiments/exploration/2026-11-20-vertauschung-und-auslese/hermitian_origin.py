"""Enumerate every sesquilinear form on the smallest finite conjugate plane.

No arguments or third-party packages are needed. 0,1,w,w+1 are coded 0,1,2,3;
the polynomial w*w+w+1 vanishes. Addition is bitwise XOR.
"""
from collections import Counter
from itertools import product, permutations
from fractions import Fraction

def times(a,b):
    term, result = a,0
    while b:
        if b&1:
            result ^= term
        term <<= 1
        if term&4:
            term ^= 7
        b >>= 1
    return result

def bar(a):
    return times(a,a)

def norm_on_ray(h,ray):
    a,b,c,d=h
    x,y=ray
    axby = times(a,x) ^ times(b,y)
    cxdy = times(c,x) ^ times(d,y)
    return times(bar(x),axby) ^ times(bar(y),cxdy)

def test():
    rays = [(0,1)] + [(1,t) for t in range(4)]
    invertible_forms, hermitian_forms, both = Counter(), Counter(), Counter()
    marks = Counter()
    for a,b,c,d in product(range(4),repeat=4):
        h = (a,b,c,d)
        null = tuple(i for i,v in enumerate(rays) if norm_on_ray(h,v)==0)
        invertible = (times(a,d) ^ times(b,c)) != 0
        hermitian = a==bar(a) and d==bar(d) and c==bar(b)
        if invertible:
            invertible_forms[len(null)] += 1
        if hermitian:
            hermitian_forms[len(null)] += 1
        if invertible and hermitian:
            both[len(null)] += 1
            marks[null] += 1
    print('null rays : count (nondegenerate sesquilinear):',sorted(invertible_forms.items()))
    print('null rays : count (all Hermitian):',sorted(hermitian_forms.items()))
    print('null rays : count (nondegenerate Hermitian):',sorted(both.items()))
    print('distinct marked triples / realisations each:',len(marks),set(marks.values()))
    assert both == {3:10} and len(marks)==10 and set(marks.values())=={1}
    h = (1,0,0,1)
    diagonal_null = [v for v in rays if norm_on_ray(h,v)==0]
    assert diagonal_null == [(1,1),(1,2),(1,3)]
    print('diagonal form null rays:',diagonal_null)
    # Controls with larger fixed subfield F_q; q²+1 lines, q+1 null lines.
    for q in (2,3,4,5):
        print('q, total, null, nonnull:',q,q*q+1,q+1,q*q-q)

def det_mod2(v):
    u,vv,a,d = v
    return (a*d + u*u - u*vv + vv*vv)%2

def eis_mul(x,y):
    # Coordinates in Z[w], w*w=-1-w.
    a,b=x
    c,d=y
    return (a*c-b*d,a*d+b*c-b*d)

def eis_conj(x):
    a,b=x
    return (a-b,-b)

def eis_neg(x):
    return tuple(-a for a in x)

def eis_norm(x):
    a,b=x
    return a*a-a*b+b*b

def diam(n):
    # Finite interval [0,nI] in the ordered integral lattice. Both endpoints
    # require 0<=a,d<=n and N(z)<=min(ad,(n-a)(n-d)).
    hist=Counter(eis_norm((u,v)) for u in range(-2*n,2*n+1)
                 for v in range(-2*n,2*n+1))
    return sum(sum(hist[k] for k in range(min(a*d,(n-a)*(n-d))+1))
               for a in range(n+1) for d in range(n+1))

def additional_tests():
    labels = list(product(range(2), repeat=4))
    zeros = {v for v in labels if not det_mod2(v) and any(v)}
    # This formula and order of labels are from lines 50-89 of the original
    # research/originlab/f1_binary_clifford_origin.py. Phases of Pauli
    # operators and the physical dictionary are not identified by this map.
    def native(v):
        u,w,a,d = v
        return (u+w+u*w+a*d)%2
    assert all(det_mod2(v)==native(v) for v in labels)
    native_order = [(0,0,1,0),(0,0,0,1),(1,0,1,1),(0,1,1,1),(1,1,1,1)]
    assert set(native_order) == zeros
    assert all((u==w==0) == ((a+d)%2==1) for u,w,a,d in zeros)
    assert all(det_mod2(tuple(a^b for a,b in zip(x,y)))==1
               for x in zeros for y in zeros if x!=y)
    ray_images=set()
    for x,y in product(range(4),repeat=2):
        if x==y==0:
            continue
        z=times(x,bar(y))
        ray_images.add((z&1,z>>1,times(x,bar(x)),times(y,bar(y))))
    assert ray_images==zeros
    print('native quadratic/marked agreement: 16/16 labels, 5/5 rays; Arf = 1')
    from itertools import combinations
    pairings = [((0,1),(2,3)), ((0,2),(1,3)), ((0,3),(1,2))]
    grammar_hits = 0
    for (i,j),(k,l) in pairings:
        for a,b in combinations(range(4),2):
            grammar_hits += all((x[i]*x[j] + x[k]*x[l] + x[a] + x[b]) % 2 == native(x)
                                for x in labels)
    assert grammar_hits == 1
    # Image-mask encoding transports the full input truth table.
    qvalues = [native(x) for x in labels]
    pair_mask = {labels.index(x) for x in native_order[:2]}
    gl_count = isometries = marked_isometries = 0
    for basis in permutations(range(1,16), 4):
        images=[0]
        for column in basis:
            images += [number ^ column for number in images]
        if len(set(images)) < 16:
            continue
        gl_count += 1
        if all(qvalues[i] == qvalues[j] for i,j in enumerate(images)):
            isometries += 1
            marked_isometries += {images[i] for i in pair_mask} == pair_mask
    assert (gl_count, isometries, marked_isometries)==(20160,120,12)
    print('same-cost fixed-coordinate quadratic formulas: 1/18; GL truth table + mark:',
          marked_isometries, '/', gl_count, '(quadratic alone:', isometries, ')')
    # The integer Gram matrix of the polarization of det is U + (-A2).
    # In coordinates (a,d,u,v) its two diagonal blocks have eigenvalues
    # +1,-1 and -1,-3, respectively.
    g1=((0,1),(1,0))
    g2=((-2,1),(1,-2))
    eig=[]
    for block, pair in ((g1,(1,-1)), (g2,(-1,-3))):
        assert sum(pair)==block[0][0]+block[1][1]
        assert pair[0]*pair[1]==block[0][0]*block[1][1]-block[0][1]**2
        eig.extend(pair)
    print('integral polarization eigenvalues:',sorted(eig),'discriminant:',-3)
    print('finite intervals [0,nI]:',[(n,diam(n)) for n in (1,2,4,8,16)])

    # Enumerate SU(2,Z[w]). A unitary integral column has unit norm, so
    # exactly one entry is an Eisenstein unit. Determinant one then fixes
    # the other column. Each off-diagonal matrix squares to -I.
    units=[(a,b) for a,b in product(range(-1,2),repeat=2)
           if eis_norm((a,b))==1]
    assert len(units)==6
    negidentity=(-1,0)
    site_actions=set()
    swap_actions=set()
    holdout=(1,0,3,2,4)
    for u in units:
        u2=eis_mul(u,u)
        for swap in (False,True):
            # diag(U,bar U) or [[0,U],[-bar U,0]]
            scale=eis_neg(u2) if swap else u2
            perm=[]
            for x,y,a,d in native_order:
                z=(x,y)
                if swap:
                    a,d=d,a
                    z=eis_conj(z)
                    assert eis_mul(u,eis_neg(eis_conj(u))) == negidentity
                z=eis_mul(scale,z)
                perm.append(native_order.index((z[0]%2,z[1]%2,a,d)))
            p=tuple(perm)
            site_actions.add(p)
            if swap:
                swap_actions.add(p)
                if u==(-1,-1):  # U=w²; the explicit witness g in the report
                    assert p==holdout
    assert len(site_actions)==6 and len(swap_actions)==3 and holdout in swap_actions
    print('integral SU(2) site image:',len(site_actions),'involutions:',len(swap_actions),
          'specified native permutation:',holdout)
    # The lifts of (1,t) and of (0,1) give positive rank-one observables.
    # Trace is one at the two axes and two on the three other directions.
    lift_vectors=[((1,0),(0,0)),((0,0),(1,0)),((1,0),(1,0)),
                  ((0,1),(1,0)),((-1,-1),(1,0))]
    tr=[eis_norm(x)+eis_norm(y) for x,y in lift_vectors]
    lift_labels=[]
    for x,y in lift_vectors:
        z=eis_mul(x,eis_conj(y))
        lift_labels.append((z[0]%2,z[1]%2,eis_norm(x)%2,eis_norm(y)%2))
    assert lift_labels==native_order
    assert [t%2 for t in tr]==[(a+d)%2 for _,_,a,d in native_order]
    print('primitive lifted traces on ordered five rays:',tr)

    ray_order=[(1,0),(0,1),(1,1),(2,1),(3,1)]
    def ray(x,y):
        if not y:
            return (1,0)
        return (times(x,bar(y)),1)
    def inner(x,y):
        return times(bar(x[0]),y[0]) ^ times(bar(x[1]),y[1])
    all_unitaries=0
    special_unitaries=0
    unitary_actions=set()
    special_actions=set()
    for a,b,c,d in product(range(4),repeat=4):
        col1,col2=(a,c),(b,d)
        if (inner(col1,col1),inner(col2,col2),inner(col1,col2))!=(1,1,0):
            continue
        all_unitaries+=1
        action=tuple(ray_order.index(ray(times(a,x)^times(b,y),times(c,x)^times(d,y)))
                     for x,y in ray_order)
        unitary_actions.add(action)
        if times(a,d)^times(b,c)==1:
            special_unitaries+=1
            special_actions.add(action)
    assert (all_unitaries,special_unitaries)==(18,6)
    assert unitary_actions==site_actions==special_actions
    identity=tuple(range(5))
    def compose(p,r):
        return tuple(p[i] for i in r)
    def order(p):
        a=p
        i=1
        while a!=identity:
            a=compose(a,p)
            i+=1
        return i
    char=Counter((order(p),sum(i==x for i,x in enumerate(p))) for p in site_actions)
    assert char=={(1,5):1,(2,1):3,(3,2):2}
    commute=sum(all(compose(p,g)==compose(g,p) for g in site_actions)
                for p in permutations(range(5)))
    assert commute==2
    print('|U(2,2)|, |SU(2,2)|, |projective site image|:',all_unitaries,
          special_unitaries,len(site_actions))
    print('site character on 1, transposition, 3-cycle: (5,1,2); equivariant self-bijections:',commute)
    print('native involution cycle type:',len(swap_actions),'/ 6; exact named action: 1 / 6')

    # Two distinct normalized RP kernels even with fixed integer trace time.
    # K(t+s) is a positive sum of outer products on every finite half-line.
    # The largest transfer eigenvalue is fixed, removing a common time rescaling.
    def k(t,r):
        return (Fraction(1,2)**t+r**t)/2
    r1,r2=Fraction(1,4),Fraction(1,8)
    for r in (r1,r2):
        assert k(0,r)==1 and k(2,r)-k(1,r)**2>0
    print('normalized RP controls at one tick:', k(1,r1),k(1,r2))

if __name__ == '__main__':
    test()
    additional_tests()
