"""Exact conformal grading repair for 16 Majoranas via Ramond spin fields.

The standard free-fermion / even-lattice VOA construction assigns weight
one to 128 spinor root currents, fixing the M3 lifted number-clock
obstruction algebraically. This is NOT an original TFPT physical
spin-structure, chiral coordinate, OS source or 3+1D derivation.
"""
from itertools import combinations,product
from collections import Counter
import sympy as sp

SIMPLE_DOUBLED=(
 (1,-1,-1,-1,-1,-1,-1,1),
 (2,2,0,0,0,0,0,0),
 (-2,2,0,0,0,0,0,0),
 (0,-2,2,0,0,0,0,0),
 (0,0,-2,2,0,0,0,0),
 (0,0,0,-2,2,0,0,0),
 (0,0,0,0,-2,2,0,0),
 (0,0,0,0,0,-2,2,0))


def controls():
    # Doubled vectors: actual E8 charge alpha is vector /2.
    so=[]
    for i,j in combinations(range(8),2):
        for a,b in product((-2,2),repeat=2):
            v=[0]*8
            v[i],v[j]=a,b
            so.append(tuple(v))
    spin=[tuple(v) for v in product((-1,1),repeat=8) if v.count(-1)%2==0]
    roots=so+spin
    rootset=set(roots)
    assert len(so)==112 and len(spin)==128 and len(rootset)==240
    assert all(sum(v*v for v in root)==8 for root in roots)
    # 16 real Majoranas in eight complex pairs -> 8 bosons H_I;
    # spin field e^{i sum eps_i H_i/2} has h=(1/2)*sum(1/2)^2=1.
    # This Ramond ground energy is independent of zero-mode occupation.
    assert all(sp.Rational(sum(v*v for v in s),8)==1 for s in spin)
    assert all(sp.Rational(sum(v*v for v in s),8)==1 for s in so)
    pairings=Counter()
    ordered_adjacent=0
    for v in roots:
        for w in roots:
            numerator=sum(a*b for a,b in zip(v,w))
            assert numerator%4==0  # a dot b in the charge lattice
            inner=numerator//4
            pairings[inner]+=1
            assert -2<=inner<=2
            if inner==-1:
                summed=tuple(a+b for a,b in zip(v,w))
                assert summed in rootset
                ordered_adjacent+=1
            if inner==-2:
                assert tuple(-a for a in v)==w
    assert ordered_adjacent==13440
    gram=sp.Matrix(8,8,lambda i,j:sum(a*b for a,b in zip(
                 SIMPLE_DOUBLED[i],SIMPLE_DOUBLED[j]))//4)
    assert gram.det()==1 and list(gram.diagonal())==[2]*8
    assert all(gram[i,j]==int(gram[i,j]) for i in range(8) for j in range(8))
    # Rank-eight positive even unimodular lattice is E8, with c=8.
    # Its 240 roots all have norm2 / holomorphic conformal weight 1.
    h_ramond=sp.Rational(16,16)
    assert h_ramond==1
    # Clock/relative bosonization phase counts from eight complex qutrit
    # traceless modes. All four sectors are balanced algebraically.
    weights=(0,0,1,1,2,2,3,3)
    root_grades=Counter((sum(root[i]*weights[i] for i in range(8))//2)%4
                         for root in roots)
    root_grades[0]+=8
    assert [root_grades[k] for k in range(4)]==[60,64,60,64]
    # Ramond ground states have common weight h=1, while ordinary
    # finite-source Gamma(qI8) has occupation energies 0,2,4,6,8.
    # This is an actual algebraic conformal-weight repair, not an
    # identity of the original microscopic source's H.
    return {
     'status':'EXACT_RAMOND_SPIN16_GRADE_ONE_E8_LATTICE_BRIDGE',
     'Majoranas':16,
     'complex_fermion_pairs':8,
     'chiral_central_charge':8,
     'Ramond_ground_spin_field_weight':str(h_ramond),
     'Ramond_full_zero_mode_multiplicity':256,
     'Ramond_positive_chirality_spin_fields':128,
     'Ramond_even_spinor_weights_all_one':True,
     'SO16_NS_bilinear_weight_one_currents':120,
     'E8_spinor_Ramond_weight_one_currents':128,
     'E8_weight_one_currents':248,
     'D8_roots':112,
     'spinor_lattice_roots':128,
     'whole_E8_root_system_size':240,
     'simple_root_Gram_determinant':1,
     'root_pairs_checked':len(roots)**2,
     'ordered_dot_minus_one_pairs':ordered_adjacent,
     'quarter_grading': [root_grades[k] for k in range(4)],
     'finite_source_number_operator_automatically_equals_VOA_L0':False,
     'Ramond_spin_structure_GSO_from_original_P1_derived':False,
     'original_local_chiral_fermion_OPE_and_OS_source_constructed':False,
     'original_physical_3plus1D_spin2_parent_derived':False,
     'F1_solved':False,
    }


if __name__=='__main__':
 import json
 print(json.dumps(controls(),indent=2,sort_keys=True))
