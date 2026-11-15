"""Select the count of marks from a weight-one local spinor extension.

For n punctures of P1, Riemann--Roch gives C_(n+1) and logarithmic
differentials F_(n-1): together 2n TYPED complex modes. If and ONLY
if the original physical source realizes ALL 2n as independent
positive chiral complex CAR fermions, it has 4n real Majoranas,
so(4n) current rank 2n and a Spin(4n) level-one spinor primary
with h=n/4. A local bosonic spinor extension requires n divisible
by four; a weight-one new current requires n=4 uniquely. This
is a conditional mathematical selection, not a TFPT source proof.
"""
from fractions import Fraction
from itertools import combinations,product


def require(ok,msg):
    if not ok:
        raise ArithmeticError(msg)


def puncture_class(n):
    require(type(n) is int and n>=2,"at least two distinct marked points")
    functions=n+1
    logforms=n-1
    complex_car=functions+logforms
    majoranas=2*complex_car
    d_rank=majoranas//2
    spinor_count=2**(d_rank-1)
    spinor_weight=Fraction(majoranas,16)
    spinor_squared=Fraction(d_rank,4)
    bilinear=majoranas*(majoranas-1)//2
    require(complex_car==2*n and majoranas==4*n and
            spinor_weight==spinor_squared/2==Fraction(n,4),
            "puncture/Spin chiral dimension dictionary changed")
    simple_local_spinor=spinor_weight.denominator==1
    new_weight_one=spinor_weight==1
    require(simple_local_spinor == (n%4==0) and
            new_weight_one == (n==4),
            "the local integral-spin extension gate changed")
    return {
      "punctures":n,
      "meromorphic_sections":functions,
      "logarithmic_differentials":logforms,
      "typed_complex_modes":complex_car,
      "candidate_real_majoranas":majoranas,
      "chiral_central_charge":2*n,
      "D_lattice_rank":d_rank,
      "D_lattice_root_count":2*d_rank*(d_rank-1),
      "SO_bilinear_current_dimension":bilinear,
      "spinor_ground_states":spinor_count,
      "spinor_squared_norm":str(spinor_squared),
      "spinor_conformal_weight":str(spinor_weight),
      "integral_spin_simple_current_allowed":simple_local_spinor,
      "new_spinor_weight_one_currents":new_weight_one,
      "spinor_overlattice_lattice_determinant_if_integral":1 if simple_local_spinor else None,
    }


def e8_weight_one_root_census():
    """Construct all D8 roots and the norm-two positive spinor coset."""
    d_rank=8
    d8=[]
    for i,j in combinations(range(d_rank),2):
        for a,b in product((-2,2),repeat=2):
            v=[0]*d_rank
            v[i],v[j]=a,b
            d8.append(tuple(v))
    spin=[v for v in product((-1,1),repeat=d_rank)
          if v.count(-1)%2==0]
    require(len(d8)==112 and len(spin)==128 and
            len(set(d8+spin))==240 and
            all(sum(x*x for x in v)==8 for v in d8+spin),
            "E8 index-two spinor gluing did not produce 240 norm-two roots")
    return len(d8),len(spin),len(d8)+len(spin)


def controls(max_marks=32):
    cases=[puncture_class(n) for n in range(2,max_marks+1)]
    integer=[c["punctures"] for c in cases
             if c["integral_spin_simple_current_allowed"]]
    weight_one=[c["punctures"] for c in cases
                if c["new_spinor_weight_one_currents"]]
    rootcounts=e8_weight_one_root_census()
    require(integer==list(range(4,max_marks+1,4)) and weight_one==[4] and
            rootcounts==(112,128,240),
            "primitive spinor weight-one selection or E8 roots changed")
    # At n=8 the even unimodular D16+ spinor has h=2:
    # a local higher-weight extension is available, but cannot
    # contribute any additional weight-one affine currents.
    c4=puncture_class(4)
    c8=puncture_class(8)
    require(c4["chiral_central_charge"]==8 and
            c4["SO_bilinear_current_dimension"]+c4["spinor_ground_states"]==248 and
            c8["spinor_conformal_weight"]=="2" and
            c8["SO_bilinear_current_dimension"]==496,
            "level-one minimal chiral mark count cannot select E8")
    return {
      "status":"EXACT_CONDITIONAL_MINIMAL_P1_SPINOR_MARK_COUNT",
      "mathematical_start":"n points on P1, sections n+1 and log forms n-1",
      "typed_complex_fermion_candidate_modes":"2n",
      "real_majorana_candidate_modes":"4n",
      "spinor_primary_conformal_weight":"h=n/4",
      "spinor_even_lattice_norm_squared":"n/2",
      "minimal_bosonic_local_spinor_extension_marks":4,
      "unique_n_with_new_spinor_weight_one_currents":4,
      "at_n4":{"meromorphic_sections":5,"log_differentials":3,
                "real_majoranas":16,"chiral_c":8,"D8_roots":112,
                "spinor_roots":128,"E8_roots":240,
                "full_E8_currents_including_Cartan":248},
      "at_n8":{"central_charge":16,"spinor_primary_weight":2,
                "spinor_adds_weight_one_current":False},
      "integer_local_extension_mark_counts_to_bound":integer,
      "candidate_original_source_independence_of_2n_CAR_modes_proved":False,
      "physical_spinor_extension_locality_and_GSO_proved":False,
      "raw_P1_square_modulus_selected":False,
      "original_weighted_R113_produced":False,
      "physical_3plus1D_parent_proved":False,
      "F1_solved":False
    }


if __name__=="__main__":
    import json
    print(json.dumps(controls(),sort_keys=True,indent=2))
