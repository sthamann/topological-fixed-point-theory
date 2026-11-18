import json

import sympy as sp

from four_mark_geometry import controls as geometry
from source import COMPONENTS, MARKS, evaluate
from verify import digest
from weak_return_completion import _operator, _return_gram, _words, controls, late_orbits


def test_return_support_and_tetrahedral_lift():
    report = controls()
    assert [row["combined_span_rank"] for row in report["support_by_grade"]] == [14,98]
    assert report["ordered_operator_span_by_maximum_length"] == [18,108,196]
    assert len(report["word_basis"]) == 196
    assert len(report["example_length_four_holdout"]) == 4
    action = geometry()
    assert action["both_holomorphic_moduli"] == []
    assert len(action["semilinear_moduli"]) == 2
    assert action["quarter_square_on_elliptic_differential"] == "+1"
    assert not action["common_invariant_reflection_circle_for_even_mark_action"]
    assert len(action["two_mark_reflection_orbit"]) == 6
    assert action["effective_mobius_translation_rank_bound_over_R"] == 2


def test_ordered_excess_reconstruction_uses_the_whole_absolute_form(tmp_path):
    roots = sum(late_orbits().values(),())
    gram = _return_gram(roots)
    metric,pivots,kets,bras,basis,_,holdout = _words(gram)
    psi = kets[:,holdout[0]]+sp.I*kets[:,holdout[-1]]
    density = psi*psi.H*metric/3+sp.eye(metric.rows)/metric.rows

    def form(scalar):
        return [[[str(sp.re(scalar)),str(sp.im(scalar))]]]

    def word_row(word):
        operator = _operator(word,gram,kets,bras)
        return {"word":list(word),"word_form":form((density*operator).trace())}

    (tmp_path/"evidence.txt").write_text("Arithmetic control")
    evidence = {"path":"evidence.txt","sha256":digest(tmp_path/"evidence.txt")}
    original = tmp_path/"original.json"
    native = {"domain_id":"ordered-excess-control","source_originals":[evidence]}
    original.write_text(json.dumps(native))
    binding = {"path":original.name,"sha256":digest(original)}
    observations = {
        "source_id":native["domain_id"]+":"+binding["sha256"],"time_id":"T",
        "charged_modes_and_adjoints_binding":"test","neutral_grade_and_excess_binding":"test",
        "absolute_reference_binding":"test","roots_doubled":[list(root) for root in roots],
        "basis_words":[word_row(word) for word in basis],"held_out_words":[word_row(holdout)],
    }
    bundle = {"schema":"tfpt.original-blocked-source.v1","domain_id":native["domain_id"],
              "native_inputs_sha256":binding["sha256"],"input_role":"arithmetic_control_only",
              "matrix_convention":"sesquilinear_forms_in_original_word_basis",
              "ensemble_scope":"whole_original_ensemble","representation":"original_exact_blocked_endpoint_source",
              "retained_words":[[]],"source_refs":[evidence],"component_bindings":dict.fromkeys(COMPONENTS,"test"),
              "time":{"id":"T","kind":"original_euclidean_ordered","T":1,"t_U":0,"t_Y":1},
              "gram_form":form(1),"physical_selector_form":form(1),"boundary_pairing_form":form(1),
              "scalar_factor":[1,0],"reference_scalar":[2,0],
              "ordered_source_forms":dict.fromkeys(MARKS,form(1)),"weak_return_grade_two":observations}
    endpoint = tmp_path/"endpoint.json"

    def read():
        endpoint.write_text(json.dumps({**native,"original_domain_binding":binding,"whole_source_bundle":bundle}))
        return evaluate(tmp_path,endpoint)["weak_return_grade_two"]

    result = read()
    assert result["status"] == "RECONSTRUCTED_POSITIVE_GRADE_TWO_EXCESS_FUNCTIONAL"
    assert result["absolute_excess_mass"] == "5/6"
    assert result["model_basis_root_indices"] == list(pivots)
    expected_form = (metric*density/2).applyfunc(sp.simplify)
    assert result["linear_candidate_form"] == [list(map(str,row)) for row in expected_form.tolist()]
    assert result["held_out"][0]["residual"] == "0"
    assert not result["physical_source_identity_proved"]
    observations["held_out_words"][0]["word_form"] = form(0)
    result = read()
    assert result["status"] == "REFUTED_SUPPLIED_GRADE_TWO_EXCESS_MODEL"
    assert result["positive_functional"]
    assert result["held_out"][0]["residual"] != "0"
