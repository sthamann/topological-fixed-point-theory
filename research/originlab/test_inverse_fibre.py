from r113_synthesis import from_input
from source import MARKS


def atom(name, mixed, mass=None):
    outputs = {key:[[[value,0]]] for key,value in zip(MARKS,(1,1,1,mixed))}
    record = {"original_configuration_id":name,"whole_fiber_binding":name,
              "original_measure_binding":"control-measure","unnormalized_scalar_weight":[1,0],
              "ordered_source_forms":outputs}
    if mass is not None:
        record["original_measure_mass"] = [mass,0]
    return record


def snapshot(fine, coarse, mapping, time="t1"):
    return {"fine_size":1,"fine_atoms":fine,"fine_scalar_factor":[1,0],"fine_reference_scalar":[1,0],
            "time":{"id":time,"t_U":0,"t_Y":1},
            "coarse_graining":{"time_id":time,"configuration_map_binding":"map",
                               "word_embedding_binding":"embedding","reference_binding":"empty",
                               "word_embedding":[[[1,0]]],"coarse_fiber_records":coarse,
                               "coarse_scalar_factor":[1,0],"coarse_reference_scalar":[1,0],
                               "fine_to_coarse":mapping}}


def search(*levels, **extra):
    return from_input({"kind":"inverse_fibres","input_role":"declared_finite_oracle",
                       "snapshots":list(levels),**extra})


def test_inverse_source_separates_nonunique_fit_from_unique_four_form_replay():
    level = snapshot([atom("a",1),atom("b",3)],[atom("c",2,1)],{"a":"c","b":"c"})
    result = search(level)
    assert result["fit"]["mass_freedom"] == 1
    assert not result["withheld_forms_predicted_from_fit"]
    assert result["separating_form"]["minimum"]["value"] == "1"
    assert result["separating_form"]["maximum"]["value"] == "3"
    assert result["status"] == "UNIQUE_MASSES_IN_DECLARED_ATOM_FAMILY"
    assert result["validation"]["masses"] == ["1/2","1/2"]
    assert result["candidate_snapshots"][0]["replay"]["checked_forms"] == 4
    assert not result["physical_source_identity_proved"]


def test_aggregate_cancellation_cannot_certify_a_source():
    # Every aggregate mark agrees. The two local mixed defects cancel.
    level = snapshot([atom("a",1),atom("b",3)],[atom("c",2,1),atom("d",2,1)],{"a":"c","b":"d"})
    result = search(level)
    assert result["fit"]["mass_freedom"] == 0
    assert result["status"] == "WITHHELD_EXCLUDED_IN_DECLARED_ATOM_FAMILY"
    dual = result["validation"]["certificate"]
    assert dual["target_value"] == "-1"
    assert "candidate_snapshots" not in result


def test_independent_ratio_and_joint_time_can_each_refute_fit():
    level = snapshot([atom("a",1),atom("b",3)],[atom("c",2,1)],{"a":"c","b":"c"})
    probe = {"id":"independent-Schur-ratio-control","coefficients":{"a":1,"b":-2},"expected":0}
    result = search(level, independent_probes=[probe])
    assert result["status"] == "WITHHELD_EXCLUDED_IN_DECLARED_ATOM_FAMILY"
    second = snapshot([atom("a",1),atom("b",3)],[atom("c","3/2",1)],{"a":"c","b":"c"},"t2")
    assert search(second)["status"] == "UNIQUE_MASSES_IN_DECLARED_ATOM_FAMILY"
    assert search(level,second)["status"] == "WITHHELD_EXCLUDED_IN_DECLARED_ATOM_FAMILY"


def test_support_functional_certifies_coefficient_uniqueness():
    level = snapshot([atom("a",1),atom("b",3)],[atom("c",1,1)],{"a":"c","b":"c"})
    result = search(level)
    assert result["validation"]["support"] == [0]
    assert result["validation"]["support_certificate"]["atom_values"] == ["0","1"]


def test_positive_cone_exclusion_and_indistinguishable_atoms():
    mapping = {"a":"c","b":"c"}
    outside = search(snapshot([atom("a",1),atom("b",3)],[atom("c",4,1)],mapping))
    assert outside["status"] == "WITHHELD_EXCLUDED_IN_DECLARED_ATOM_FAMILY"
    assert outside["validation"]["certificate"]["target_value"] == "-1"
    same = search(snapshot([atom("a",2),atom("b",2)],[atom("c",2,1)],mapping))
    assert same["status"] == "MASSES_UNDERDETERMINED_IN_DECLARED_ATOM_FAMILY"
    assert same["withheld_forms_predicted_from_fit"]
    assert [item["atom"] for item in same["proportional_atom_classes"][0]] == ["a","b"]


def test_exact_coarse_replay_can_leave_the_fine_source_undetermined():
    first,second = atom("a",2),atom("b",2)
    for record, hidden in ((first,1),(second,3)):
        for form in record["ordered_source_forms"].values():
            form[0].append([0,0])
            form.append([[0,0],[hidden,0]])
    level = snapshot([first,second],[atom("c",2,1)],{"a":"c","b":"c"})
    level["fine_size"] = 2
    level["coarse_graining"]["word_embedding"] = [[[1,0]],[[0,0]]]
    result = search(level)
    assert result["candidate_snapshots"][0]["replay"]["status"] == "EXACT_SUPPLIED_MARKED_PUSHFORWARD"
    assert not result["candidate_snapshots"][0]["fine_Z0_fixed_by_full_constraints"]
