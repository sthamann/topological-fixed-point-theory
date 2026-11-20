"""E8 in the complementary-pair frame of the eight subsets of three labels."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product

# Subsets are bit masks. Each complementary pair is oriented canonically
# by taking its member of smaller cardinality first.
universe=7
pairs=((0,7),(1,6),(2,5),(4,3))
family=frozenset((3,5,7))  # difference coordinates of the three active pairs
carrier=frozenset(range(8))-family


def affine_code():
    return {tuple(((v&x).bit_count()+b)%2 for x in range(8)) for v in range(8) for b in (0,1)}


def frame_roots():
    # Physical vector = z/sqrt(2) in the point basis; in the f frame = t/2.
    code=affine_code()
    roots=set()
    for i,sgn in product(range(8),(-1,1)):
        z=[0]*8
        z[i]=2*sgn
        roots.add(tuple(z))
    for w in code:
        support=[i for i in range(8) if w[i]]
        if len(support)!=4:
            continue
        for sign in product((-1,1),repeat=4):
            z=[0]*8
            for i,s in zip(support,sign):
                z[i]=s
            roots.add(tuple(z))
    return {tuple(y for i,j in pairs for y in (z[i]+z[j],z[i]-z[j])) for z in roots}


def subset_map(axes):
    return tuple(sum(1<<axes[i] for i in range(3) if x & 1<<i) for x in range(8))


def compose(p,q):
    return tuple(p[q[i]] for i in range(len(p)))


def permute_signs(p,mask):
    return sum(1<<p[i] for i in range(8) if mask & 1<<i)


def enumerate_lifts():
    # Store diagonal signs on the OUTPUT coordinates: U=(diagonal)*permutation.
    # a is a transposition, b a three-cycle. a^2=b^3=(ab)^2=1 presents S3.
    a=subset_map((1,0,2))
    b=subset_map((1,2,0))
    ab=compose(a,b)
    pmap={p:[permute_signs(p,m) for m in range(256)] for p in (a,b,compose(b,b),ab)}
    total=[]
    for s,t in product(range(256),repeat=2):
        if s^pmap[a][s] or t^pmap[b][t]^pmap[compose(b,b)][t]:
            continue
        u=s^pmap[a][t]
        if u^pmap[ab][u]:
            continue
        total.append((s,t))
    def preserves_active_differences(sign):
        return all(((sign>>i)^(sign>>j))&1==0 for i,j in pairs[1:])
    compatible=[(s,t) for s,t in total if all(map(preserves_active_differences,(s,t)))]
    swaps=[(s,t) for s,t in compatible if ((s>>0)^(s>>7))&1]
    parity_fixed=[(s,t) for s,t in swaps if not ((s|t)&1)]
    active_positive=[(s,t) for s,t in parity_fixed if all((v & 126)==0 for v in (s,t))]
    assert active_positive==[(128,0)]
    # e_0 is the extra Fano check coordinate, not a chosen Hecke vacuum.
    # Alternative: preserve positive point vectors on all seven original
    # coordinates, allowing signs ONLY on the newly appended e_0.
    old_positive=[(s,t) for s,t in total if ((s|t)&254)==0]
    old_positive_swap=[p for p in old_positive if p in swaps]
    assert set(old_positive)=={(0,0),(1,0)} and old_positive_swap==[(1,0)]
    # For ref orientations (e_7+e_0)/sqrt(2),(e_7-e_0)/sqrt(2),
    # sign(e_0) interchanges the actual vectors with positive coefficients.
    assert [(1,-1),(1,1)]==[(x,-y) for x,y in ((1,1),(1,-1))]
    assert all(s.bit_count()%2==1 for s,t in swaps)
    def gauge_classes(lifts, allowed_masks):
        # A gauge affects signs by the coboundary d XOR permutation(d).
        remaining=set(lifts)
        count=0
        while remaining:
            s,t=next(iter(remaining))
            orbit={(s^d^pmap[a][d],t^d^pmap[b][d]) for d in allowed_masks}
            remaining-=orbit
            count+=1
        return count
    gauges=[m for m in range(256) if preserves_active_differences(m)]
    print("S3 diagonal sign lifts, all and up to coordinate-sign gauge:",len(total),gauge_classes(total,range(256)))
    print("preserve active difference lines, all and up to pair-preserving gauge:",
          len(compatible),gauge_classes(compatible,gauges))
    print("also swap reference lines for odd permutations, all / gauge classes:",
          len(swaps),gauge_classes(swaps,gauges))
    print("also fix appended positive parity coordinate, all / gauge classes:",
          len(parity_fixed),gauge_classes(parity_fixed,gauges))
    print("also preserve positivity of the other six point axes: (sign masks for a,b)",active_positive)
    print("alternative, all 7 old point vectors stay positive; lifts / swap lifts:",old_positive,old_positive_swap)
    # A *signed integer* equivariant lift of the old parity function would
    # send the original all-ones word (fixed by S3) to a sign-fixed integer,
    # hence zero; its mod-two parity is 1. This isometry does not lift the encoder.
    assert (1,)*7 in {tuple(c[1:]) for c in affine_code()} and sum((1,)*7)%2==1


def hexagon_probe():
    # The order complex of nonempty proper subsets is a 6-cycle. Edge
    # orientations from singleton to doubleton are themselves S3-invariant.
    # Comparing a primitive cycle with its image determines H1 integrally.
    cycle=(1,3,2,6,4,5)
    oriented_edges={}
    for i,x in enumerate(cycle):
        y=cycle[(i+1)%6]
        edge=(x,y) if x.bit_count()==1 else (y,x)
        oriented_edges[edge]=1 if edge==(x,y) else -1
    assert len(oriented_edges)==6
    for p in permutations(range(3)):
        act=subset_map(p)
        inv=sum(p[i]>p[j] for i,j in combinations(range(3),2))
        sign=(-1)**inv
        assert all(oriented_edges[act[x],act[y]]==sign*v for (x,y),v in oriented_edges.items())
    print("hexagon H0 and H1 (transposition): +1, -1")


def main():
    roots=frame_roots()
    assert len(roots)==240
    int_roots={t for t in roots if all(v%2==0 for v in t)}
    spin={t for t in roots if all(abs(v)==1 for v in t)}
    assert len(int_roots)==112 and len(spin)==128
    assert all(sorted(abs(v) for v in t)==[0]*6+[2,2] for t in int_roots)
    spin_parities={sum(v<0 for v in t)%2 for t in spin}
    assert len(spin_parities)==1
    print("D8 vector roots / half-spin weights / parity of minus count:",len(int_roots),len(spin),spin_parities.pop())

    def supported(t,block):
        return all(t[k]==0 for k in range(8) if k not in block)
    fam={t for t in roots if supported(t,family)}
    car={t for t in roots if supported(t,carrier)}
    assert (len(fam),len(car))==(12,40)
    def chi(t):return sum(t[i] for i in family)%4
    counts=Counter(map(chi,roots))
    assert {t for t in roots if chi(t)==0}==fam|car
    assert [counts[i] for i in range(4)]==[52,64,60,64]
    print("D3(A3), D5 roots:",len(fam),len(car))
    print("glue degree 0,1,2,3 roots:",[counts[i] for i in range(4)])
    print("generator norms in D3*, D5*, sum:",Fraction(3,4),Fraction(5,4),Fraction(8,4))

    # Coordinate-root controls: every f-frame triple is D3 with D5 complement.
    hit=0
    for ids in combinations(range(8),3):
        block=frozenset(ids)
        if sum(supported(t,block) for t in roots)==12 and sum(supported(t,set(range(8))-block) for t in roots)==40:
            hit+=1
    assert hit==56
    symbols=[subset_map(p) for p in permutations(range(3))]
    def line_action(p):
        return tuple(2*next(j for j,(u,v) in enumerate(pairs) if {u,v}=={p[x],p[y]})+s
                     for x,y in pairs for s in (0,1))
    acts=list(map(line_action,symbols))
    triplets=[c for c in combinations(range(8),3) if all({p[i] for i in c}==set(c) for p in acts)]
    assert set(triplets)=={(2,4,6),(3,5,7)}
    print("D3+D5 among frame triples / S3-invariant triples:",hit,triplets)
    enumerate_lifts()
    hexagon_probe()


if __name__=="__main__":
    main()
