import json
from itertools import product
from pathlib import Path

import pytest
import sympy as sp

from native_weak_w import replay as native_weak_w
from source import COMPONENTS, MARKS, evaluate
from verify import digest
from weak_clock import controls


ROOT = Path(__file__).resolve().parents[2]


def test_literal_fold_separates_the_root_lift_in_both_intervals():
    report = controls(ROOT)
    assert report["root_lift_mode_pairs_count"] == 96
    lift = report["joint_colour_lift"]
    assert lift["early_ranks_per_grade"] == [3,2]
    assert lift["map_space_rank"] == 6 and lift["insertion_multiplet"] == "(3,2)_(1/6)"
    assert [case["early_gap_K"] for case in report["cases"]] == ["3/16","1/16"]
    assert [case["late_gap_K"] for case in report["cases"]] == ["3/16","1/16"]
    assert {case["clock_preserving_maps_to_weak_endpoints"] for case in report["cases"]} == {0}
    assert set(report["marked_early_residues"]).isdisjoint(report["weak_early_residues"])
    assert report["noncommuting_two_quarter_norm_ratio_bound"] == "1024/18225"
    assert report["fold_factor_rectangular_ratio_on_each_axis_(1/4,1/16)"]["weak"] == "85/81"
    assert report["independent_power_seven_outside_character_residual"] == "144/125 + 144*I/125"


def test_whole_source_preserves_reference_and_separates_clock_ratios(tmp_path):
    def form(value):
        re,im = sp.expand(value).as_real_imag()
        return [[[str(re),str(im)]]]

    original = tmp_path/"original.json"
    evidence = tmp_path/"evidence.txt"
    evidence.write_text("Seven joint eigenvalues in an arithmetic control")
    source_ref = {"path":evidence.name,"sha256":digest(evidence)}
    native = {"domain_id":"weak-clock-control","source_originals":[source_ref]}
    original.write_text(json.dumps(native))
    binding = {"path":original.name,"sha256":digest(original)}
    source_id = native["domain_id"]+":"+binding["sha256"]
    z = (3+4*sp.I)/5
    sectors = [-sp.I,z,sp.conjugate(z),-sp.I*z,-sp.I*sp.conjugate(z),sp.I*z,sp.I*sp.conjugate(z)]
    samples = []
    widths = ("1/4","1/9","1/36","1/49","1/64")
    for x,y in [*product(widths,repeat=2),("1/16","1/25")]:
        ratio = sp.sqrt(sp.Rational(x)*sp.Rational(y))
        weights = [1+sp.I,2,2*ratio,3,3,7/(1+sp.Rational(x)*sp.Rational(y)),0]
        powers = [sum(w*v**k for w,v in zip(weights,sectors)) for k in range(8)]
        samples.append({"left_radial_factor":x,"right_radial_factor":y,"held_out":x=="1/16",
                        "power_forms":[form(2*value) for value in powers]})
    clock = {"source_id":source_id,"time_id":"T","samples":samples,
             **dict.fromkeys(("common_fields_and_frames_binding","fixed_joint_early_action_binding",
                              "same_histories_for_all_powers_binding","physical_to_radial_time_binding",
                              "absolute_reference_binding"),"arithmetic control")}
    packet = {"schema":"tfpt.original-blocked-source.v1","domain_id":native["domain_id"],
              "native_inputs_sha256":binding["sha256"],"input_role":"arithmetic_control_only",
              "matrix_convention":"sesquilinear_forms_in_original_word_basis",
              "ensemble_scope":"whole_original_ensemble",
              "representation":"original_exact_blocked_endpoint_source",
              "retained_words":[[]],"source_refs":[source_ref],
              "component_bindings":dict.fromkeys(COMPONENTS,"control"),
              "time":{"id":"T","kind":"original_euclidean_ordered","T":1,"t_U":0,"t_Y":1},
              "gram_form":form(1),"physical_selector_form":form(1),"boundary_pairing_form":form(1),
              "scalar_factor":[1,0],"reference_scalar":[2,0],
              "ordered_source_forms":dict.fromkeys(MARKS,form(2)),"weak_clock_sectors":clock}
    attached = tmp_path/"attached.json"

    def read():
        attached.write_text(json.dumps({**native,"original_domain_binding":binding,"whole_source_bundle":packet}))
        return evaluate(tmp_path,attached)["weak_clock_sectors"]

    result = read()
    assert result["status"] == "SUPPLIED_SEVEN_SECTOR_MODEL_REPLAYED"
    assert not result["physical_source_identity_proved"]
    holdout = result["samples"][-1]
    assert holdout["absolute_sector_amplitudes"]["spinor"] == "1 + I"
    assert holdout["absolute_sector_amplitudes"]["d8_upper"] == "2"
    assert holdout["absolute_sector_amplitudes"]["d8_lower"] == "1/10"
    assert holdout["weak_weight_comparison"]["d8"] == {
        "measured_ratio":"1/20","matches_literal_fold":True,"matches_equal_weights":False}
    assert holdout["weak_weight_comparison"]["weak_spinor_minus"]["matches_equal_weights"]
    assert holdout["power_seven_residual"] == "0"
    ranks = result["resolved_time_rank"]
    assert ranks["status"] == "FIXED_NATIVE_W_DRESSING_EXCLUDED"
    assert ranks["separating_sectors"] == ["weak_spinor_plus_upper"]
    assert ranks["complex_sector_ranks"]["weak_spinor_plus_upper"] == 5

    for value in samples[-1]["power_forms"]:
        re_im = value[0][0]
        re_im[0] = str(sp.Rational(re_im[0])+2)
    result = read()
    assert result["status"] == "SUPPLIED_SEVEN_SECTOR_MODEL_REFUTED"
    assert result["resolved_time_rank"]["status"] == "UNVERIFIED_SECTOR_AMPLITUDES"
    assert result["samples"][-1]["power_seven_residual"] == "144/125 + 144*I/125"
    clock["time_id"] = "other time"
    with pytest.raises(ValueError,match="same source|one original source"):
        read()


def test_native_w_word_has_exact_separation_rank_four():
    report = native_weak_w(ROOT)
    assert report["cross_chart_bracket_support_pairs"] == 3840
    assert report["nonzero_bracket_pairs"] == 960
    phase = report["phase_and_symmetry_replay"]
    assert phase["common_phase_chart"]["constraints"] == 1500
    assert phase["common_phase_chart"]["single_coefficient_mutant_rejected"]
    assert phase["weak_components_related_by_actual_W_symmetries"] == 12
    reduction = report["dressed_word_exact_time_reduction"]
    assert len(reduction["right"]["boson_numbers"]) == 4
    assert len(reduction["left"]["boson_numbers"]) == 13
    assert len(reduction["middle_form"]) == 13
    assert reduction["unnormalized_moment_minor_orders_zero_to_three"] == "22096281600*abs_g**12"
    assert reduction["mixed_time_rank_at_nonzero_g"] == 4
    assert reduction["mixed_time_rank_at_zero_g"] == 1
    assert report["signed_affine_cauchy_ranks_checked"] == [2,3,4,5]
    assert not report["grade_preserving_even_spinor_field_possible"]
    assert report["original_r113_values"] is None
    json.dumps(report)
