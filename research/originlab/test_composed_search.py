from composed_search import run, square_control


def test_mixed_words_select_operational_class():
    case = square_control()
    selected = run(case)
    assert selected["fit_survivors"] == 20
    assert selected["withheld_survivors"] == 2
    assert selected["unique_operational_class"]
    assert selected["all_word_classes"][0]["size"] == 2
    assert not selected["physical_source_identity_proved"]

    open_case = run({**case,"withheld":[]})
    witness = open_case["next_distinguishing_word"]
    assert open_case["equivalence_complete"] and not open_case["unique_operational_class"]
    assert witness["word"] not in [observation["word"] for observation in case["fit"]]
    assert witness["responses"][0] != witness["responses"][1]
    unresolved = run({**case,"equivalence_budget":0})
    assert not unresolved["equivalence_complete"] and not unresolved["unique_operational_class"]
    assert unresolved["withheld_survivors"] == 2


def test_shared_bindings_prevent_cross_source_combinations():
    def slot(name):
        return {"letters":[name],"choices":[{"id":value, "bindings":{"chart":chart},
                                               "operators":{name:[[value]]}}
                                              for chart,value in (("a",1),("b",2))]}
    result = run({"gram":[[1]],"trace_weight":[[1]],"slots":[slot("U"),slot("Y")],
                  "withheld":[{"word":["U","Y"],"value":2}]})
    assert result["fit_survivors"] == 2
    assert result["pruned_prefixes"]["binding"] == 2
    assert result["status"] == "NO_MATCH_IN_DECLARED_COMPOSITION_CLASS"


def test_boundary_state_and_dynamics_are_selected_together():
    case = {"boundary_choices":[{"id":name,"gram":[[1]],"trace_weight":[[mass]],"bindings":{"seam":name}}
                                for name,mass in (("a",1),("b",2))],
            "slots":[{"letters":["T"],"choices":[{"id":name,"bindings":{"seam":name},
                                                  "operators":{"T":[[step]]}}
                                                 for name,step in (("a",1),("b","1/2"))]}],
            "positive_steps":["T"],"fit":[{"word":["T"],"value":1}],
            "withheld":[{"word":["T","T"],"value":"1/2"}]}
    result = run(case)
    assert result["fit_survivors"] == 2
    assert result["all_word_classes"] == [{"representative":[{"boundary":"b"},"b"],"size":1}]
