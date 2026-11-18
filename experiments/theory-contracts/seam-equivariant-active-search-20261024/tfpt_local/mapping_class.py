"""Pure four-puncture Hurwitz action in the original quarter-clock stratum.

Declared extra equivalence: adjacent squared Hurwitz moves b0^2,b1^2.
No quotient of fixed reflection-marked objects is silently taken here.
"""
from collections import Counter
import json
from pathlib import Path

import numpy as np

from .engine import save_json
from .research import Transport, Task, act, linear


def hurwitz(group, vector, index):
    values=[(vector >> (8*j)&255,1,0) for j in range(4)]
    x,y=values[index:index+2]
    values[index:index+2]=[group.multiply(group.multiply(x,y),group.inverse(x)),x]
    return sum(c << (8*j) for j,(c,k,e) in enumerate(values))


def main():
    group=Transport(); task=Task(group)
    W=linear(task.rows[:8],32)['null_basis']
    N=lambda c:c ^ act(group.A,c)
    gauge=[sum(N(1 << i) << (8*j) for j in range(4)) for i in range(8)]
    pivots={}
    for x in gauge:
        while x:
            p=x.bit_length()-1
            if p in pivots: x ^= pivots[p]
            else: pivots[p]=x; break

    def reduce(values):
        out=np.asarray(values,dtype=np.uint32).copy()
        for p,v in sorted(pivots.items(),reverse=True):
            out ^= np.where((out >> p)&1,np.uint32(v),np.uint32(0))
        return out

    def apply_A(values):
        out=np.zeros_like(values,dtype=np.uint32)
        for j in range(4):
            for i,c in enumerate(group.A):
                out ^= np.where((values >> (8*j+i))&1,np.uint32(c << (8*j)),np.uint32(0))
        return out

    def pure(values,index):
        ci=(values >> (8*index))&255
        cj=(values >> (8*(index+1)))&255
        u=ci^cj
        au=np.zeros_like(values,dtype=np.uint32)
        for i,c in enumerate(group.A): au ^= np.where((u >> i)&1,np.uint32(c),np.uint32(0))
        u ^= au
        return values ^ (u << (8*index)) ^ (u << (8*(index+1)))

    # Check the proposed square action against exact original group products.
    for i in range(3):
        for v in W:
            expected=hurwitz(group,hurwitz(group,v,i),i)
            assert int(pure(np.uint32(v),i))==expected
            assert not any((r & expected).bit_count()%2 for r in task.rows[:8])
    Q=[int(reduce(np.uint32(v))) for v in W]
    independent={}
    for x in Q:
        while x:
            p=x.bit_length()-1
            if p in independent:x ^= independent[p]
            else:independent[p]=x;break
    assert len(independent)==20
    values=np.zeros(1 << 20,dtype=np.uint32)
    length=1
    for v in independent.values():
        values[length:2*length]=values[:length]^np.uint32(v);length*=2
    # Every character gauge quotient element is enumerated exactly once.
    assert len(np.unique(values))==1 << 20
    assert np.array_equal(reduce(values),values)
    alpha=reduce(apply_A(values))
    p0=reduce(pure(values,0));p1=reduce(pure(values,1))
    assert np.array_equal(reduce(pure(p0,0)),values)
    assert np.array_equal(reduce(pure(p1,1)),values)
    assert np.array_equal(reduce(pure(p0,1)),reduce(pure(p1,0)))
    assert np.array_equal(reduce(pure(alpha,0)),reduce(apply_A(p0)))
    assert np.array_equal(reduce(pure(alpha,1)),reduce(apply_A(p1)))
    # Complementary-pair twists agree modulo character gauge on this stratum.
    assert np.array_equal(p0,reduce(pure(values,2)))
    representatives=values.copy();fix_counts=[]
    for mask in range(8):
        transformed=values
        if mask&1:transformed=apply_A(transformed)
        if mask&2:transformed=pure(transformed,0)
        if mask&4:transformed=pure(transformed,1)
        transformed=reduce(transformed)
        fix_counts.append(int(np.sum(transformed==values)))
        representatives=np.minimum(representatives,transformed)
    reps,sizes=np.unique(representatives,return_counts=True)
    burnside=sum(fix_counts)//8
    assert len(reps)==burnside==157696
    alpha_fixed=reduce(apply_A(reps))==reps
    class_sizes=sizes//np.where(alpha_fixed,1,2)
    census=Counter(map(int,class_sizes))
    assert sum(size*count for size,count in census.items())==526336
    result={'status':'exact finite orbit classification with stated pure Hurwitz equivalence',
        'source_id':group.source.identity,'model_contribution':False,
        'coefficient_space':'V=Hom(E8,F2), native J action A, N=I+A, N^2=0',
        'squared_hurwitz_formula':'c_i -> c_i+N(c_i+c_(i+1)); c_(i+1) -> c_(i+1)+N(c_i+c_(i+1))',
        'generator_names':['b0^2','b1^2'],'action_image':'C2 x C2 on original Gamma moduli',
        'gauge_quotient_dimension':20,'all_quotient_points_checked':len(values),
        'fixed_points_of_A_P0_P1_combinations':fix_counts,
        'orbits_after_pure_mapping_class_equivalence':len(reps),
        'pure_mapping_class_orbits_by_number_of_original_Gamma_classes':dict(census),
        'normal_form':'After reducing by diagonal Im(I+A), choose the least nonnegative encoded integer among A^e P0^f P1^g images, e,f,g in F2',
        'scope':'All ordered four-J sign-lift tuples, first Gamma conjugation and then pure mapping classes; no extra reflection marking, no puncture permutations, no physical-origin claim',
        'topological_generator_reference':'https://archive.mpim-bonn.mpg.de/5138/1/mpim-preprint_2024-26.pdf, Section 2'}
    save_json('results/pure-mapping-class.json',result)
    np.savez_compressed('results/pure-mapping-class-normal-forms.npz',representatives=reps,
                        original_Gamma_class_count=class_sizes)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
