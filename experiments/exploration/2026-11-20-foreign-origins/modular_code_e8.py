#!/usr/bin/env python3
"""Independent exhaustive Type-II code control from modularly compatible bosonic lattice.
No E8 or Hamming code is inserted; enumerate all small doubly-even self-dual codes.
This *is* a known route already present in OriginLab; it is not a new TFPT origin.
"""
from itertools import permutations
from collections import Counter
from fractions import Fraction
import json


def admissible_codes(n):
    words=[x for x in range(1,1<<n) if x.bit_count()%4==0]
    spaces={(0,)}
    counts={0:1}
    for dim in range(1,n//2+1):
        next_spaces=set()
        for S in spaces:
            avail=[w for w in words if w not in S and all((w&s).bit_count()%2==0 for s in S)]
            for w in avail:
                T=tuple(sorted(set(S)|{s^w for s in S}))
                assert len(T)==2**dim
                next_spaces.add(T)
        spaces=next_spaces
        counts[dim]=len(spaces)
    return spaces,counts


def gaussian_binom_q2(n,k):
    f=Fraction(1)
    for i in range(k):f *= Fraction(2**n-2**i,2**k-2**i)
    assert f.denominator==1
    return f.numerator


def permute_word(w,perm):
    return sum(((w>>i)&1)<<perm[i] for i in range(len(perm)))


def orbit_size(code):
    n=8
    return len({tuple(sorted(permute_word(w,p) for w in code)) for p in permutations(range(n))})


def test():
    by_n={}
    for n in (2,4,6,8):
        codes, counts=admissible_codes(n)
        by_n[n]=(codes,counts)
    codes=by_n[8][0]
    assert all(len(by_n[n][0])==0 for n in (2,4,6))
    assert len(codes)==30
    ref=next(iter(codes))
    weight=dict(sorted(Counter(w.bit_count() for w in ref).items()))
    assert weight=={0:1,4:14,8:1}
    assert all(Counter(w.bit_count() for w in c)==Counter({0:1,4:14,8:1}) for c in codes)
    orb=orbit_size(ref)
    assert orb==len(codes)==30
    count_roots=2*8+16*weight[4]
    assert count_roots==240
    count_total=gaussian_binom_q2(8,4)
    assert count_total==200787
    return {'n2_n4_n6_type_II_codes':[len(by_n[n][0]) for n in (2,4,6)],'min_nonzero_type_II_code_length':8,
            'rank8_subspace_count_all':count_total,
            'rank8_type_II_labelled_codes':len(codes),
            'fraction_type_II_exact':f'{len(codes)}/{count_total}',
            'single_coordinate_permutation_orbit_size':orb,
            'unique_up_to_coordinate_permutation':True,
            'code_weight_enumerator':weight,
            'construction_A_covolume_exact':str(Fraction(2**(8-4),2**(8//2))),
            'construction_A_short_roots':{'coordinate_pm2':16,'four_pm1':224,'total':count_roots},
            'search_counts_by_code_dimension':by_n[8][1]}


if __name__=='__main__':
    print(json.dumps(test(),indent=2,ensure_ascii=False))
