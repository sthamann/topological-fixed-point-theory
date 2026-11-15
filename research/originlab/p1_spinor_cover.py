"""P1 meromorphic five-mode even spinor: native Z4 256/248 phase bridge.

The geometric C4 action on H0(P1,O(mu4)) gives five characters
(1,1,i,-1,-i). Their even exterior space is 16-dimensional.
End of this space is a *graded vector-space* spectral cover for
the affine E8 degree-one currents. It is not an E8 Lie action.
"""
from collections import Counter
from itertools import combinations
from affine_operators import ROOTS, charge_signature


def require(ok, why):
    if not ok:
        raise ArithmeticError(why)


def controls():
    # f0,..,f4=z^k/(z^4-1) (basis weight k mod 4).
    charges=(0,1,2,3,0)
    even_subsets=tuple(subset for n in (0,2,4)
                       for subset in combinations(range(5),n))
    spinor_degree=tuple(sum(charges[i] for i in subset)%4
                        for subset in even_subsets)
    spinor=Counter(spinor_degree)
    spinor_counts=[spinor[k] for k in range(4)]
    require(len(even_subsets)==16 and spinor_counts==[4,4,4,4],
            "original five-mode even spinor Z4 character is not uniform")

    adjoint=Counter((left-right)%4
                    for left in spinor_degree for right in spinor_degree)
    adjoint_counts=[adjoint[k] for k in range(4)]
    require(adjoint_counts==[64,64,64,64],
            "endomorphism of original 16 spinor not balanced as C4 module")

    native=Counter(charge_signature(root)[0] for root in ROOTS)
    native[0]+=8  # affine grade-one: 240 roots + 8 Cartan
    native_counts=[native[k] for k in range(4)]
    deficit=[a-b for a,b in zip(adjoint_counts,native_counts)]
    require(native_counts==[60,64,60,64] and deficit==[4,0,4,0],
            "native marked E8 current grading does not match spinor cover")

    # The natural geometric unitary has 4 equal eigenclasses in S+,
    # so its adjoint has 64 of every Z4 grade. This is a matching of
    # characters, NOT a representation of E8 by 16x16 matrices.
    trace_square_spinor=sum((-1)**k*spinor[k] for k in range(4))
    trace_square_end=sum((-1)**k*adjoint[k] for k in range(4))
    trace_square_native=sum((-1)**k*native[k] for k in range(4))
    require((trace_square_spinor,trace_square_end,trace_square_native)==(0,0,-8),
            "native parity-trace obstruction is not reproduced")
    return {
      "status":"EXACT_EVEN_SPINOR_Z4_VECTOR_COVER",
      "original_five_mode_C4_weights":[1,1,"i","-1","-i"],
      "even_exterior_16_dimension":len(even_subsets),
      "even_exterior_C4_phase_counts":spinor_counts,
      "endomorphism_256_dimension":len(even_subsets)**2,
      "endomorphism_C4_phase_counts":adjoint_counts,
      "native_affine_E8_grade_one_dimension":sum(native_counts),
      "native_affine_E8_grade_one_phase_counts":native_counts,
      "minimal_even_grade_complement_dimensions":deficit,
      "minimal_total_null_complement":sum(deficit),
      "native_current_parity_trace":trace_square_native,
      "cover_current_parity_trace":trace_square_end,
      "graded_vector_space_embedding_exists":True,
      "canonical_current_intertwiner_identified":False,
      "e8_lie_bracket_embedded_in_gl16":False,
      "e8_bracket_obstruction":"smallest nontrivial complex E8 representation has dimension 248, not 16",
      "spinor16_equals_16_Majorana_fields_proved":False,
      "physical_R113_source_proved":False,
      "F1_solved":False,
    }


if __name__=="__main__":
    import json
    print(json.dumps(controls(),indent=2,sort_keys=True))
