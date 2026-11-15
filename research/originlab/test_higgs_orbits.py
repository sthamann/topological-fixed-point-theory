import sympy as sp

from higgs_moment_search import field_coefficients, norm_coefficients, probe_plan
from higgs_orbits import even_word, gram_for_field, orbit_plan, phase_compatible


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
