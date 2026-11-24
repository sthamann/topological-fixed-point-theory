#!/usr/bin/env python3
"""Exact route from the *minimal* perfect five-qubit code to the E8 isomorphism type.

The [[5,1,3]] code determines a 4-bit syndrome space and five projective
lines partitioning its fifteen nonzero elements. The *chosen* operation
"take any rank-3 syndrome hyperplane, build the binary Hamming kernel,
then append overall parity and apply Construction A" returns the E8 code
isomorphism type for ALL 15 choices. No hyperplane is preferred by C5.
This is a structural conditional bridge, not a derivation of the TFPT source.
"""
from itertools import combinations
from collections import Counter
from math import comb
import json
from sympy import Matrix, Symbol, Poly, factor
from five_qubit_qec import GENS, xz, syndrome, code_data


def rot5(mask):return ((mask<<1)&31) | (mask>>4)


def linear_map(v, columns):
    out=0
    for j in range(4):
        if v>>j&1:out ^= columns[j]
    return out


def gf2_polymul(a,b,poly=0b11111):
    out=0
    while b:
        if b&1:out ^= a
        b >>= 1
        a <<= 1
        if a&0b10000:a ^= poly
    return out


def gfpow(a,n):
    out=1
    for _ in range(n):out=gf2_polymul(out,a)
    return out


def hyperplane(a):
    return tuple(v for v in range(1,16) if (v&a).bit_count()%2==0)


def hamming_extension(H):
    assert len(H)==7
    C7=[]
    for w in range(1<<7):
        s=0
        for j,v in enumerate(H):
            if w>>j&1:s ^=v
        if s==0:C7.append(w)
    assert len(C7)==16
    C8=[w | ((w.bit_count()%2)<<7) for w in C7]
    dist7=dict(sorted(Counter(w.bit_count() for w in C7).items()))
    dist8=dict(sorted(Counter(w.bit_count() for w in C8).items()))
    assert dist7=={0:1,3:7,4:7,7:1}
    assert dist8=={0:1,4:14,8:1}
    assert all((a&b).bit_count()%2==0 for a in C8 for b in C8)
    return C8,dist7,dist8


def test():
    data=code_data(GENS)
    gs=data[0]
    one_error={}
    lines=[]
    for i in range(5):
        site=[]
        for p in 'XYZ':
            e=xz('I'*i+p+'I'*(4-i))
            s=syndrome(e,gs)
            assert 1<=s<=15 and s not in one_error
            one_error[s]=e
            site.append(s)
        assert site[0]^site[1]^site[2]==0
        lines.append(tuple(sorted(site)))
    assert len(one_error)==15 and len(set(lines))==5
    assert set(v for line in lines for v in line)==set(range(1,16))
    # Action induced by cyclically permuting the qubits, transported to syndromes.
    cols=[]
    for j in range(4):
        e=one_error[1<<j]
        cols.append(syndrome((rot5(e[0]),rot5(e[1])),gs))
    for s,e in one_error.items():
        assert linear_map(s,cols)==syndrome((rot5(e[0]),rot5(e[1])),gs)
    assert all(linear_map(linear_map(linear_map(linear_map(linear_map(v,cols),cols),cols),cols),cols)==v
               for v in range(16))
    assert all(linear_map(v,cols)!=v for v in range(1,16))
    x=Symbol('x')
    A=Matrix([[((cols[j]>>i)&1) for j in range(4)] for i in range(4)])
    cp=Poly(A.charpoly(x).as_expr(),x,modulus=2)
    assert cp.as_expr()==x**4+x**3+x**2+x+1 and cp.is_irreducible
    # An eigenvalue of the order-5 action is a nontrivial fifth root in F_16.
    t=0b0010
    assert gfpow(t,5)==1 and all(gfpow(t,k)!=1 for k in (1,2,3,4))
    frob=[gfpow(t,2**j) for j in range(5)]
    assert frob[0]==frob[4] and len(set(frob[:4]))==4
    # The multiplication-by-order-5 and fourth-order Frobenius maps yield 20
    # DISTINCT permutations: C5 semidirect C4, with Frobenius(t*v)=t^2*Frobenius(v).
    semilinear_actions={tuple(gf2_polymul(gfpow(t,k),gfpow(v,2**j))
                              for v in range(16)) for k in range(5) for j in range(4)}
    assert len(semilinear_actions)==20
    assert all(gfpow(gf2_polymul(t,v),2)==gf2_polymul(gfpow(t,2),gfpow(v,2))
               for v in range(16))
    # All 15 3D syndrome hyperplanes yield the same extended [8,4,4] Hamming type.
    profiles=[]
    hyperplanes=[hyperplane(a) for a in range(1,16)]
    assert len(set(hyperplanes))==15
    for H in hyperplanes:
        code,d7,d8=hamming_extension(H)
        profiles.append((tuple(d7.items()),tuple(d8.items())))
        assert 2*8+16*d8[4]==240
        assert {linear_map(v,cols) for v in H}!=set(H)
    assert len(set(profiles))==1
    # Exact matched-size null among arbitrary seven-point subsets of F_2^4 \ {0}.
    closed=0
    for S in combinations(range(1,16),7):
        U=set(S)|{0}
        if all((a^b) in U for a in U for b in U):closed+=1
    assert closed==15 and comb(15,7)==6435
    return {
        'code_single_error_syndromes':15,
        'syndrome_space':'F2^4',
        'five_site_lines_in_PG3F2':lines,
        'partition_15_into_5x3':True,
        'cyclic_site_action_columns_GF2':cols,
        'cyclic_site_action_characteristic_polynomial_GF2':str(cp.as_expr()),
        'irreducible_degree':4,
        'F16_frobenius_orbit_of_primitive_fifth_root':frob,
        'F16_frobenius_order':4,
        'semilinear_group_C5_semidirect_C4_order':len(semilinear_actions),
        'warning':'Frobenius C4 is not the same as the TFPT mu4 deck action',
        'syndrome_hyperplanes':len(hyperplanes),
        'all_hyperplanes_extend_to_Hamming_8_4_4_and_E8':True,
        'Hamming7_weight_distribution':dict(profiles[0][0]),
        'Hamming8_weight_distribution':dict(profiles[0][1]),
        'E8_norm2_vectors_all_choices':240,
        'C5_invariant_3D_hyperplanes':0,
        'null_same_size_7_of_15_subsets_that_are_linear_hyperplanes':f'{closed}/{comb(15,7)}',
        'new_indispensable_free_choice':'take a syndrome hyperplane and extend its Hamming kernel; no distinguished embedding',
        'not_derived':'TFPT marked 3+2 species, P1, continuum dynamics, Lorentz metric, physical c3'}


if __name__=='__main__': print(json.dumps(test(),indent=2,ensure_ascii=False))
