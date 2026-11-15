import json
from itertools import product

import pytest
import sympy as sp

from source import COMPONENTS, MARKS, evaluate
from verify import digest
from weak_source import (CHANNELS, PROBES, WINDING_PROBES, affine_controls, checked_parity_sector,
                         example_networks, joint_response_matrix, transport_character,
                         triangle_components, winding_characters)


def test_exact_modes_and_joint_probe_resolve_all_weak_completion_orbits():
    report = affine_controls(2)
    assert report["mode_checks"] == 84
    assert report["casimir_checks"] == 336
    assert report["colour_and_weak_invariance_checks"] == 24
    assert report["d8_root_string_checks"] == 24
    assert report["root_centre_counts"] == {0:38,1:36,2:45,3:40,4:45,5:36}
    assert list(report["character_controls"]["weak_first"].values()) == ["1","0","0","0"]
    assert set(report["character_controls"]["weak_inverse"].values()) == {"1"}
    assert report["original_coweight_controls"]["spinor"]["interval_winding"] == [-1,-1]
    assert report["original_coweight_controls"]["weak_spinor_minus"]["interval_winding"] == [-1,-1]
    fixed = report["fixed_middle_search"]
    assert fixed["examined_roots"] == 240
    assert fixed["paths"]["weak_spinor_minus"]["interval_glue"] == [3,1]
    assert fixed["paths"]["spinor"]["interval_glue"] == [1,3]
    assert fixed["joint_power_response"][1] == ["-I","3/5","-3*I/5","3*I/5"]
    assert fixed["powers_zero_to_three_determinant"] != "0"
    assert report["polar_parity_necessary_extension_control"]["vacuum_pair_order_coefficients"] == ["0","1"]


def test_weighted_word_reader_distinguishes_equal_untwisted_amplitudes(tmp_path):
    def pair(value):
        return [str(sp.re(value)),str(sp.im(value))]

    def form(value):
        return [[pair(value)]]

    original = tmp_path/"original.json"
    evidence = {"path":"evidence.txt"}
    (tmp_path/evidence["path"]).write_text("Arithmetic control")
    evidence["sha256"] = digest(tmp_path/evidence["path"])
    native = {"domain_id":"weak-source-control","source_originals":[evidence]}
    original.write_text(json.dumps(native))
    binding = {"path":original.name,"sha256":digest(original)}
    source_id = native["domain_id"]+":"+binding["sha256"]
    network = example_networks()
    triangles = triangle_components()
    response = joint_response_matrix(triangles)

    def link(edge):
        w,c,u = edge
        return {"weak":[[pair(value) for value in row] for row in w.tolist()],
                "colour":[[pair(value) for value in row] for row in c.tolist()],"u1":pair(u)}

    pairs = [*product(("1/2","1/3"),repeat=2), ("1/4","1/5")]
    samples = []
    for i,(x,y) in enumerate(pairs):
        bare = sp.Integer(2+i)
        words = {}
        for name,channels in CHANNELS.items():
            values = [c for channel in channels for c in channel.casimirs()]
            _,*signs = winding_characters(triangles[name])
            words[name] = {"bare_form":form(bare),
                           "casimir_forms":dict(zip(PROBES,(form(value*bare) for value in values))),
                           "winding_forms":dict(zip(WINDING_PROBES,(form(sign*bare) for sign in signs))),
                           "twisted_words":{key:{"links":[link(edge) for edge in edges],
                                                 "word_form":form(bare*transport_character(name,*edges))}
                                            for key,edges in network.items() if key!="identity"}}
        samples.append({"left_radial_factor":x,"right_radial_factor":y,
                        "held_out":i==len(pairs)-1,"words":words,
                        "joint_early_forms":[form(value) for value in response*sp.Matrix([bare,1,3,0])]})
    gauge = {"source_id":source_id,"time_id":"T","samples":samples,
             "comparison_paths":list(CHANNELS),
             **dict.fromkeys(("native_three_current_binding","same_dressing_and_inserted_actions_binding",
                              "casimir_insertions_binding","physical_to_radial_time_binding",
                              "absolute_reference_binding","original_coweight_binding",
                              "joint_early_insertion_binding"),"arithmetic control")}
    packet = {"schema":"tfpt.original-blocked-source.v1","domain_id":native["domain_id"],
              "native_inputs_sha256":binding["sha256"],"input_role":"arithmetic_control_only",
              "matrix_convention":"sesquilinear_forms_in_original_word_basis",
              "ensemble_scope":"whole_original_ensemble",
              "representation":"original_exact_blocked_endpoint_source",
              "retained_words":[[]],"source_refs":[evidence],
              "component_bindings":dict.fromkeys(COMPONENTS,"control"),
              "time":{"id":"T","kind":"original_euclidean_ordered","T":1,"t_U":0,"t_Y":1},
              "gram_form":form(1),"physical_selector_form":form(1),"boundary_pairing_form":form(1),
              "scalar_factor":[1,0],"reference_scalar":[2,0],
              "ordered_source_forms":dict.fromkeys(MARKS,form(1)),"weak_gauge_responses":gauge}
    attached = tmp_path/"attached.json"

    def read():
        attached.write_text(json.dumps({**native,"original_domain_binding":binding,"whole_source_bundle":packet}))
        return evaluate(tmp_path,attached)["weak_gauge_responses"]

    result = read()
    assert result["status"] == "MATCHES_SUPPLIED_AFFINE_RESPONSES"
    assert result["samples"][0]["spinor"]["absolute_bare"] == "1"
    assert result["samples"][0]["spinor"]["winding"]["early_half"]["absolute"] == "-1"
    assert result["samples"][0]["d8"]["absolute_casimir_words"]["early_weak"] == "3/4"
    assert result["samples"][0]["joint_early"]["absolute_path_coefficients"] == dict(zip(CHANNELS,["1","1/2","3/2","0"]))
    assert not result["physical_source_identity_proved"]
    holdout = samples[-1]["words"]["d8"]
    original_form = holdout["twisted_words"]["weak_first"]["word_form"]
    holdout["twisted_words"]["weak_first"]["word_form"] = holdout["bare_form"]
    result = read()
    assert result["status"] == "DISTINGUISHED_SUPPLIED_GAUGE_RESPONSES"
    assert result["samples"][-1]["d8"]["transport"]["weak_first"]["residual"] == "3"
    holdout["twisted_words"]["weak_first"]["word_form"] = original_form
    samples[-1]["joint_early_forms"][4] = form(0)
    result = read()
    assert result["status"] == "DISTINGUISHED_SUPPLIED_GAUGE_RESPONSES"
    assert result["samples"][-1]["joint_early"]["power_four_residual"] != "0"
    samples.pop()
    with pytest.raises(ValueError,match="holdout"):
        read()


def test_parity_bridge_rejects_an_odd_field_even_if_it_has_the_right_cartan_weight():
    basis = sp.eye(3)
    alpha = basis[:,2]*basis[1,:]
    delta = basis[:,2]*basis[0,:]
    data = {"parity_form":sp.diag(1,-1,-1),"alpha_cartan_form":sp.diag(0,-1,1),
            "alpha_raising_form":alpha,"alpha_lowering_form":alpha.H,
            "delta_form":delta,"delta_adjoint_form":delta.H,
            "neutral_projector_form":sp.diag(1,0,0),"early_projector_form":sp.zeros(3),
            "late_projector_form":sp.diag(0,0,1)}
    result = checked_parity_sector(data,basis,lambda value: value)
    assert result["status"]=="INCOMPATIBLE_SUPPLIED_ALPHA_WORDS"
    assert result["k_minus_rank"]==0
    assert result["residual_ranks"]["delta_parity"]==1
    assert result["residual_ranks"]["delta_cartan"]==0
