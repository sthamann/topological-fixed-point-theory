"""Exact geometric mu4 (5+3) to native E8 order-four glue conjugacy.

This proves a FULL 240-root Weyl equivalence of the two C4 characters,
not merely a 248-dimensional phase census. Its physical assumption is
that the original meromorphic 5 and logarithmic differential 3 are
INDEPENDENT positive chiral complex fermions in the *original* source.
That source identification, BRST reduction, time, R113 and 3+1D F1
remain open.
"""
from collections import Counter, deque
from itertools import combinations
from affine_operators import ROOTS, charge_signature, components

# These eight modes are the DISJOINT, TYPED direct sum of the original
# meromorphic-section space (5) and logarithmic-differential space (3).
# Their geometric r:z->i*z degrees mod4:
GEOM8=(0,1,2,3,0,1,2,3)
# Coordinate reordering is harmless; use the grouped weights below for
# explicit root torus conjugacy, which matches the standard D8 chart.
GEOM_TORUS=(0,0,1,1,2,2,3,3)
NATIVE_GLUE_TORUS=(0,0,0,0,0,2,2,2)

# Eight doubled E8 simple roots, with Cartan determinant one and
# E8 branch-arm lengths 1,2,4. Their geometric E8 lattice is the
# ROOTS set in affine_operators.
SIMPLE_DOUBLED=(
 (1,-1,-1,-1,-1,-1,-1,1),
 (2,2,0,0,0,0,0,0),
 (-2,2,0,0,0,0,0,0),
 (0,-2,2,0,0,0,0,0),
 (0,0,-2,2,0,0,0,0),
 (0,0,0,-2,2,0,0,0),
 (0,0,0,0,-2,2,0,0),
 (0,0,0,0,0,-2,2,0),
)
# One explicit sequence of E8 Weyl simple-root reflections, 1-based.
# This maps 2*GEOM_TORUS to 2*NATIVE_GLUE_TORUS modulo 8*E8.
WEYL_WORD=(4,2,3,1,4,5,4,2,3,1,4,3,5,6,5,4,
           2,3,1,7,6,5,4,2,3,1,4,5,6,8,7)


def require(ok,why):
    if not ok:
        raise ArithmeticError(why)


def character(root, doubled_torus):
    """i^(root·torus) with both E8 root and torus in doubled coordinates."""
    s=sum(a*b for a,b in zip(root,doubled_torus))
    require(s%4==0, "torus must pair integrally with every E8 root")
    return (s//4)%4


def reflect(doubled_torus, doubled_root):
    scalar=sum(a*b for a,b in zip(doubled_torus,doubled_root))
    require(scalar%4==0,"Weyl E8 reflection loses integral lattice coordinates")
    return tuple(a-(scalar//4)*b
                 for a,b in zip(doubled_torus,doubled_root))


def controls():
    require(sum(GEOM8)%4==0 and sum(GEOM_TORUS)%4==0,
            "geometric 8 mode clock determinant not unity")
    require(set(SIMPLE_DOUBLED)<=set(ROOTS) and len(SIMPLE_DOUBLED)==8,
            "declared simple roots are not native E8 roots")
    # Even Fock spinor of 8 complex fermions: 2^7=128.
    even_phases=Counter(
      sum(GEOM8[j] for j in subset)%4
      for n in (0,2,4,6,8)
      for subset in combinations(range(8),n))
    require([even_phases[k] for k in range(4)]==[32,32,32,32],
            "geometric eight-mode positive half-spinor phase census changed")

    # Real 16-dimensional vector of Majorana fields: annihilators
    # and creators, hence phases q and -q for each complex mode.
    vector_phases=Counter(
      [k%4 for k in GEOM8]+[(-k)%4 for k in GEOM8])
    require([vector_phases[k] for k in range(4)]==[4,4,4,4],
            "geometric CAR realification does not yield 16 Majoranas")

    # so16=Lambda² of 16 real vector modes, dim120.
    vector_list=[k%4 for k in GEOM8]+[(-k)%4 for k in GEOM8]
    so16=Counter(
      (vector_list[i]+vector_list[j])%4
      for i,j in combinations(range(16),2))
    require([so16[k] for k in range(4)]==[28,32,28,32],
            "16 Majorana bilinear currents changed")

    # Exact ROOT-BASED E8 grading of the geometric torus.
    qgeom=tuple(2*k for k in GEOM_TORUS)
    qnative=tuple(2*k for k in NATIVE_GLUE_TORUS)
    geom=Counter(character(r,qgeom) for r in ROOTS)
    geom[0]+=8 # affine grade-one includes 8 Cartan currents
    native=Counter(charge_signature(r)[0] for r in ROOTS)
    native[0]+=8
    require([geom[k] for k in range(4)]==
            [native[k] for k in range(4)]==[60,64,60,64],
            "geometric C4 fixed algebra not native E8 current grade")
    fixed=[r for r in ROOTS if character(r,qgeom)==0]
    require(len(fixed)==52 and components(fixed)==["A3","D5"],
            "geometric fixed Lie algebra not A3+D5")
    # 120 currents come from D8 roots plus Cartan, 128 from
    # spinor roots: compare the two contributions separately.
    d8=Counter(character(r,qgeom) for r in ROOTS
               if all(x%2==0 for x in r))
    d8[0]+=8
    spinor=Counter(character(r,qgeom) for r in ROOTS
                   if all(x%2 for x in r))
    require([d8[k] for k in range(4)]==[28,32,28,32] and
            [spinor[k] for k in range(4)]==[32]*4,
            "E8 D8+spinor decomposition not reproduced by geometric character")

    transformed=qgeom
    for k in WEYL_WORD:
        transformed=reflect(transformed,SIMPLE_DOUBLED[k-1])
    # Explicit witness: w(qgeom) = qnative - 4*alpha_7.
    r7=SIMPLE_DOUBLED[6]
    expected=tuple(qnative[i]-4*r7[i] for i in range(8))
    require(transformed==expected,
            "stated 31-reflection Weyl witness is incorrect")
    require(all(character(r,transformed)==charge_signature(r)[0]
                for r in ROOTS),
            "Weyl-conjugated geometric C4 fails a native E8 root charge")

    return {
      "status":"EXACT_GEOMETRIC_C4_TO_NATIVE_E8_Weyl_CONJUGACY",
      "original_four_mark_mode_dimensions":[5,3],
      "geometric_complex_CAR_mode_number":8,
      "geometric_real_majorana_candidate_count":16,
      "geometric_complex_C4_degrees":list(GEOM8),
      "sixteen_real_majorana_C4_phase_counts":[4,4,4,4],
      "so16_bilinear_current_phase_counts":[28,32,28,32],
      "positive_128_spinor_phase_counts":[32,32,32,32],
      "full_E8_248_phase_counts":[60,64,60,64],
      "native_E8_248_phase_counts":[native[k] for k in range(4)],
      "geometric_fixed_root_system":["A3","D5"],
      "geometric_fixed_root_count":len(fixed),
      "explicit_Weyl_reflection_word_1based":list(WEYL_WORD),
      "weyl_word_length":len(WEYL_WORD),
      "geometric_torus_weight":list(GEOM_TORUS),
      "native_glue_torus_weight":list(NATIVE_GLUE_TORUS),
      "weyl_transformed_doubled_torus":list(transformed),
      "weyl_difference":"w(q_geom) = q_native - 4*alpha_7",
      "all_240_native_root_phase_equalities_proved":True,
      "geometry_to_fermion_CAR_source_identified":False,
      "five_plus_three_independent_physical_modes_proved":False,
      "original_Q0_and_R113_produced":False,
      "physical_local_3plus1D_parent_proved":False,
      "F1_solved":False,
    }


if __name__=="__main__":
    import json
    print(json.dumps(controls(),sort_keys=True,indent=2))
