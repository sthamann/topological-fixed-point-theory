"""Exact free Ramond OPE data and isotropic gluing of D_(m+1)+A_(m-1)."""
from fractions import Fraction as F
from math import comb, gcd, lcm
import json

# Circle vacuum energy of one real Majorana in a common zeta convention.
def zeta_minus_one(a):
    return -(a*a-a+F(1,6))/2

e_ns = -zeta_minus_one(F(1,2))/2
e_r = -zeta_minus_one(F(1))/2
h_r = e_r-e_ns
assert (e_ns,e_r,h_r)==(F(-1,48),F(1,24),F(1,16))

def spin_ope(complex_count):
    # Fix the first Ramond weight s=(1/2)^N. At Hamming distance d
    # the second weight has pairing (N-2d)/4 and outgoing norm N-d.
    h = complex_count*h_r*2
    rows = []
    for d in range(complex_count+1):
        leading = F(complex_count-2*d,4)
        output_h = F(complex_count-d,2)
        assert output_h-2*h==leading
        if d%2==0:
            assert leading.denominator==1
        else:
            assert leading.denominator==2
        rows.append([d,comb(complex_count,d),str(leading),str(output_h)])
    return {"Majoranas":2*complex_count, "weight":str(h),
            "all_zero_mode_states":2**complex_count,
            "per_chirality":2**(complex_count-1),
            "S_lambda_S_mu_columns": ["Hamming distance","count from fixed weight",
                                      "leading power","output weight"],
            "S_lambda_S_mu":rows}

def qD(m,j):
    # m even => D_(m+1) has cyclic discriminant of order four.
    assert m%2==0
    return [F(0),F(m+1,4),F(1),F(m+1,4)][j%4]

def qA(m,j):
    j %= m
    return F(j*(m-j),m)

def certify_cyclic(m,d,a):
    order = lcm(4//gcd(4,d),m//gcd(m,a))
    group = [(t*d%4,t*a%m) for t in range(order)]
    assert len(set(group))==order
    assert all((qD(m,x)+qA(m,y))%2==0 for x,y in group)
    det = F(4*m,order**2)
    return {"marks":m,"generator":[d,a],"order":order,
            "total_rank":2*m,"determinant":str(det),
            "generator_h_min":str((qD(m,d)+qA(m,a))/2),
            "pure_D":sum(y==0 for x,y in group),
            "pure_A":sum(x==0 for x,y in group)}

# Both factors primitive: order-two graph exists iff m=4 mod 8;
# order-four graph exists iff m=4 mod 32.
for m in range(2,258,2):
    assert ((qD(m,2)+qA(m,m//2))%2==0)==(m%8==4)
    if m%4==0:
        assert ((qD(m,1)+qA(m,m//4))%2==0)==(m%32==4)

partial = [certify_cyclic(12,2,6),certify_cyclic(36,1,9)]
assert [p["determinant"] for p in partial]==["12","9"]
assert all(p["pure_A"]==p["pure_D"]==1 for p in partial)

unimodular = [certify_cyclic(4*k*k,1,k) for k in range(1,11)]
for k,row in enumerate(unimodular,1):
    assert row["determinant"]=="1" and row["order"]==4*k
    assert row["pure_D"]==1 and row["pure_A"]==k

r8 = spin_ope(8)
r16 = spin_ope(16)
# D8 s·c, s·v and c·v test independent from any cocycle.
s = (F(1,2),)*8
v = (0,)*7+(1,)
c = tuple(x-y for x,y in zip(s,v))
pair = lambda a,b: sum(x*y for x,y in zip(a,b))
assert (pair(s,c),pair(s,v),pair(c,v))==(F(3,2),F(1,2),F(-1,2))

print(json.dumps({
    "Majorana_NS_energy":str(e_ns), "Majorana_R_energy":str(e_r),
    "Majorana_R_relative_weight":str(h_r), "one_complex_channel_R_weight":str(2*h_r),
    "D8":r8, "D16":r16,
    "D8_pairings_s_c__s_v__c_v":[str(pair(s,c)),str(pair(s,v)),str(pair(c,v))],
    "partial_primitive_counterexamples":partial,
    "unimodular_nonprimitive_family":unimodular,
},indent=2))
