"""Exact grade-resolved lattice charges and bare chiral spectrum controls."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, product
from math import prod
from pathlib import Path
import json
import sys

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo / "verification"))
from v128_graded_hull import build_roots

def d_vectors(n):
    for i in range(n):
        for a in (-1, 1):
            yield tuple(a if j == i else 0 for j in range(n))

def d_roots(n):
    for i,j in combinations(range(n),2):
        for a,b in product((-1,1), repeat=2):
            yield tuple(a if k == i else b if k == j else 0 for k in range(n))

def a_roots():
    for i in range(4):
        for j in range(4):
            if i != j:
                yield tuple(1 if k == i else -1 if k == j else 0 for k in range(4))

def a_min(k):
    for occupied in combinations(range(4),k):
        yield tuple(F(int(i in occupied)) - F(k,4) for i in range(4))

def norm2(v):
    return sum(x*x for x in v)

def independent_roots(n):
    roots = {(x,(0,)*4):0 for x in d_roots(n)}
    roots.update({((0,)*n,f):0 for f in a_roots()})
    roots.update({(x,f):2 for x,f in product(d_vectors(n),a_min(2))})
    if n == 5:
        for signs in product((-1,1),repeat=n):
            x = tuple(F(s,2) for s in signs)
            g = 1 if signs.count(-1)%2 == 0 else 3
            roots.update({(x,f):g for f in a_min(g)})
    return roots

assert independent_roots(5) == build_roots()
Y = (-F(1,3),-F(1,3),-F(1,3),F(1,2),F(1,2))
A = (1,-1,0,0)

def spectrum(n):
    roots = independent_roots(n)
    grades = Counter(roots.values())
    charges = [(g,sum(x[i]*Y[i] for i in range(5)),sum(t*f for t,f in zip(A,z)))
               for (x,z),g in roots.items()]
    assert all(norm2(x)+norm2(z)==2 for x,z in roots)
    # Same A3 fundamental weight as the v983 generator. A Dn Weyl
    # reflection changing both weak signs supplies a charged component.
    omega_f = next(a_min(1))
    highest = (F(1,2),)*n
    weak_flipped = highest[:3] + (-F(1,2),)*2 + highest[5:]
    h = (norm2(weak_flipped)+norm2(omega_f))/2
    highest_y = sum(x*y for x,y in zip(highest,Y))
    charged_y = sum(x*y for x,y in zip(weak_flipped,Y))
    family_charge = sum(x*a for x,a in zip(omega_f,A))
    assert highest_y==0 and charged_y==-1 and family_charge==1
    # For <H(0)V(1)V*(3)>, with unit two-point normalization.
    ward_factor = -F(2,3)*F(1,2)**int(2*h)
    return {
        "norm2_roots_by_glue": [grades[g] for g in range(4)],
        "full_weight_one_by_glue": [grades[g]+(n+3 if g==0 else 0) for g in range(4)],
        "minimum_conformal_weights": ["0",str(F(n+3,8)),"1",str(F(n+3,8))],
        "minimal_odd_multiplicity_per_grade":2**(n-1)*4,
        "Y2_by_grade":[str(sum(y**2 for g,y,f in charges if g==k)) for k in range(4)],
        "A2_by_grade":[str(sum(f**2 for g,y,f in charges if g==k)) for k in range(4)],
        "YA_by_grade":[str(sum(y*f for g,y,f in charges if g==k)) for k in range(4)],
        "Y4_by_grade":[str(sum(y**4 for g,y,f in charges if g==k)) for k in range(4)],
        "Y2A2_by_grade":[str(sum(y*y*f*f for g,y,f in charges if g==k)) for k in range(4)],
        "joint_charge_moment_uniform_roots":str(sum(y*y*f*f for _,y,f in charges)/len(roots)),
        "odd_component_test":{
            "literal_highest_Y":str(highest_y),
            "weak_flipped_Y":str(charged_y),
            "family_coroot_charge":str(family_charge),
            "conformal_weight":str(h),
            "unit_normalized_Y_Ward_at_0_1_3":str(ward_factor*charged_y),
            "unit_normalized_A_Ward_at_0_1_3":str(ward_factor*family_charge),
        },
    }

def spinor(n):
    parity_pairs = Counter()
    Yhist = Counter()
    for signs in product((-1,1),repeat=n):
        # Match the even-occupation convention of the v1002 Weyl weights.
        if signs.count(1)%2:
            continue
        parity_pairs[(signs[:5].count(1)%2,signs[5:].count(1)%2)] += 1
        y = sum(F(x,2)*yi for x,yi in zip(signs,Y))
        Yhist[y] += 1
    return {
        "dimension":sum(Yhist.values()),
        "chirality_blocks":{str(a):b for a,b in sorted(parity_pairs.items())},
        "Y_by_multiplicity":{str(a):b for a,b in sorted(Yhist.items())},
        "trace_Y_power_1_2_3":[str(sum(y**p*k for y,k in Yhist.items())) for p in (1,2,3)],
    }

def d5_representation_dimensions():
    positive_roots = [v for v in d_roots(5) if next(x for x in v if x)!=-1]
    rho = (4,3,2,1,0)
    weights = [(1,0,0,0,0),(1,1,0,0,0),(1,1,1,0,0),
               (F(1,2),)*4+(-F(1,2),),(F(1,2),)*5]
    def dim(w):
        out = prod(F(sum((a+b)*v for a,b,v in zip(w,rho,r)),
                     sum(b*v for b,v in zip(rho,r))) for r in positive_roots)
        assert out.denominator == 1
        return int(out)
    fundamental = [dim(w) for w in weights]
    quadratic_min = min(dim(tuple(a+b for a,b in zip(x,y)))
                        for x,y in combinations_with_replacement(weights,2))
    assert fundamental == [10,45,120,16,16] and quadratic_min==54
    return {"fundamental":fundamental,"minimum_sum_of_two_fundamental_weights":quadratic_min}

s5,s13 = spectrum(5),spectrum(13)
assert [sum(map(F,s["Y2_by_grade"])) for s in (s5,s13)] == [50,50]
assert [sum(map(F,s["A2_by_grade"])) for s in (s5,s13)] == [120,120]
assert [sum(map(F,s["Y4_by_grade"])) for s in (s5,s13)] == [25,F(145,9)]
assert [sum(map(F,s["Y2A2_by_grade"])) for s in (s5,s13)] == [20,F(20,3)]
assert s5["Y2A2_by_grade"]==["0","20/3","20/3","20/3"]
assert s13["Y2A2_by_grade"]==["0","0","20/3","0"]
print(json.dumps({
    "D5_A3":s5,
    "D13_A3":s13,
    "Spin10_halfspinor":spinor(5),
    "Spin26_halfspinor":spinor(13),
    "D5_dimensions_for_all_orthogonal_26d_embeddings":d5_representation_dimensions(),
},indent=2))
