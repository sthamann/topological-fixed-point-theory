import sympy as sp

from higgs_moment_search import norm_coefficients, probe_plan


def test_invisible_harmonic_can_be_signed_by_the_norm_cross_term():
    proposal = {"base_constant":[0],"allowed_multipliers":[-1,0,1],
                "slots":[{"harmonic":1,"kind":"cos","component":0,"unit":"13/100"},
                         {"harmonic":8,"kind":"cos","component":0,"unit":"13/100"}]}
    negative = norm_coefficients(proposal,[2,0])
    positive = norm_coefficients(proposal,[2,2])
    assert positive[7] == -negative[7] == sp.Rational(169,20000)
    assert positive[16] == negative[16] == sp.Rational(169,40000)
    report = probe_plan(proposal)
    assert report["history_count"] == 9
    assert report["independent_moment_count_including_mass"] == 5
    measurement = sp.Matrix([[sp.Rational(v) for v in row["value_per_class"]]
                             for row in report["measurement_plan"]])
    reconstruction = sp.Matrix([[sp.Rational(v) for v in row]
                                for row in report["class_totals_from_ordered_measured_moments"]])
    assert reconstruction*measurement == sp.eye(5)
    assert not report["physical_source_derived"]


def test_moment_products_complete_an_incomplete_linear_sensor_set():
    proposal = {"base_constant":[0,0],"allowed_multipliers":[-1,0,1],
                "slots":[{"harmonic":1,"kind":"cos","component":0,"unit":1},
                         {"harmonic":2,"kind":"sin","component":1,"unit":1}]}
    report = probe_plan(proposal)
    assert report["independent_moment_count_including_mass"] == 4
    assert max(len(row["word"]) for row in report["measurement_plan"]) == 2
