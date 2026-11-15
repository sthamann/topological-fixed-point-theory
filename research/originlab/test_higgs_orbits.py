import sympy as sp

from higgs_moment_search import field_coefficients, norm_coefficients, probe_plan
from higgs_orbits import even_word, gram_for_field, orbit_plan, phase_compatible
from r113_synthesis import from_input
from source import MARKS


PROPOSAL = {"base_constant":[0,0,1,0],"allowed_multipliers":[-1,0,1],
            "slots":[{"harmonic":1,"kind":"cos","component":0,"unit":"13/100"},
                     {"harmonic":1,"kind":"sin","component":1,"unit":"9/100"},
                     {"harmonic":2,"kind":"cos","component":3,"unit":"7/100"}]}


def test_second_harmonic_pair_has_identical_norm_and_distinct_orbits():
    negative, positive = [2,2,0],[2,2,2]
    assert norm_coefficients(PROPOSAL,negative)==norm_coefficients(PROPOSAL,positive)
    plan = orbit_plan(PROPOSAL,probe_plan(PROPOSAL)["norm_signature_classes"])
    assert plan["constant_gauge_and_rotation_orbits"] == 14
    assert plan["including_spatial_reflection_orbits"] == 11
    assert plan["norm_equal_distinct_orbit_pairs"] == 3
    witness = plan["second_harmonic_sign_witness"]
    assert not witness["same_gauge_spatial_orbit"]
    assert witness["values"] == ["77*I/1000000","-77*I/1000000"]


def test_phase_word_is_unchanged_by_rotating_and_reflecting_the_chart():
    g = gram_for_field(field_coefficients(PROPOSAL,[2,2,0]),range(-2,3))
    z = (3+4*sp.I)/5
    twin = {(n,m):sp.cancel(value*z**(m-n)) for (n,m),value in g.items()}
    reflected = {(n,m):twin[-n,-m] for n,m in g}
    assert phase_compatible(g,twin)
    assert phase_compatible(g,reflected,True)
    assert even_word(g,(-1,1),(0,2)) == even_word(reflected,(-1,1),(0,2))


def test_unit_modulus_per_difference_is_not_enough_for_one_translation():
    g = {(0,0):1,(0,2):1,(0,3):1}
    assert not phase_compatible(g,{**g,(0,2):-1})


def test_quotient_moment_inverse_and_positive_support_witness():
    fine = {**PROPOSAL,"slots":[*PROPOSAL["slots"],
                                {"harmonic":8,"kind":"cos","component":0,"unit":"13/100"}]}
    for proposal,count in ((PROPOSAL,11),(fine,22)):
        report = orbit_plan(proposal,probe_plan(proposal)["norm_signature_classes"])
        plan = report["weighted_orbit_measurement_plan"]
        assert plan["separates_all_orbits"]
        assert plan["independent_moments_including_empty"] == count
        assert max(len(row["word"]) for row in plan["moments"]) == 2
        measurement = sp.Matrix([[sp.Rational(x) for x in row["value_per_class"]] for row in plan["moments"]])
        inverse = sp.Matrix(plan["class_totals_from_ordered_measured_moments"]).applyfunc(sp.Rational)
        weights = sp.Matrix([sp.Rational(i+1,count*(count+1)//2) for i in range(count)])
        assert inverse*measurement*weights == weights

        test = plan["positive_support_test"]
        assert {"kind":"reflected_diagonal","mode":0,"component":"real","value":"1"} in test["constant_coordinates"]
        zeros = [list(map(sp.Rational,row)) for row in test["zeros"]]
        scales = list(map(sp.Rational,test["coordinate_scales"]))
        def polynomial(point):
            return sp.prod(sum(((x-y)/scale)**2 for x,y,scale in zip(point,zero,scales)) for zero in zeros)
        assert len({tuple(row) for row in zeros}) == count
        assert all(polynomial(point)==0 for point in zeros)
        outside = zeros[0][:]
        outside[0] = max(row[0] for row in zeros)+scales[0]
        assert polynomial(outside)>0


def test_original_format_replays_all_four_matrix_valued_moments():
    report = orbit_plan(PROPOSAL,probe_plan(PROPOSAL)["norm_signature_classes"])
    plan = report["weighted_orbit_measurement_plan"]
    truth = [sp.Matrix([[i,sp.I*j],[-sp.I*j,2*i]]) for i,j in zip(range(1,12),range(11))]
    def pairs(form):
        return [[[str(x) for x in v.as_real_imag()] for v in row] for row in form.tolist()]
    measured = []
    for row in plan["moments"]:
        form = sum((sp.Rational(a)*m for a,m in zip(row["value_per_class"],truth)),sp.zeros(2))
        measured.append({"word":row["word"],"ordered_source_forms":{
            mark:pairs(scale*form) for scale,mark in enumerate(MARKS,1)}})
    result = from_input({"kind":"orbit_source_moments","input_role":"declared_finite_oracle",
                         "field_proposal":PROPOSAL,"word_space_size":2,"measurements":measured})
    assert result["status"] == "EXACT_FORMS_ON_DECLARED_PROBE_CLASSES"
    assert result["reconstructed_marked_forms"] == [
        {mark:pairs(scale*form) for scale,mark in enumerate(MARKS,1)} for form in truth]
