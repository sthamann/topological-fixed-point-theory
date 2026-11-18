from complex_moments import analyze, control
from marked_search import phase_and_clock, run


def test_marked_campaign_and_probe_feedback():
    found = run()
    marks, cells = found["marks"], found["cells"]
    assert marks["phase"]["remaining"] == [0,1]
    assert marks["clock"]["initial_clocks"] == 30
    assert marks["clock"]["initial_involutions"] == 15
    assert marks["clock"]["next_probe"]["worst_case_total_probes"] == 3
    assert abs(int(marks["meromorphic_clock"]["change_of_basis_determinant"])) == 4
    assert marks["integral_glue"]["root_counts"] == [40,50,50,50,50]
    assert set(cells["remaining"]) == {"I","-I","R","-R"}
    assert cells["next_probe"]["separated_pairs"] == 6
    replay = run({"cells":{"observations":{"mixed_1_1":"1/2"}}})
    assert replay["cells"]["remaining"] == ["R"]
    pruned = phase_and_clock({"observations":{"clock":{"image_of_point_1":2}}})
    assert pruned["clock"]["remaining_clocks"] == 6


def test_complex_phase_witness_and_fixed_observations():
    result = control()
    assert result["next_probe"]["step"] == 1
    assert result["next_probe"]["component"] == "imaginary"
    assert {w["mixed_value"] for w in result["witnesses"]} == {"I/4","-I/4"}
    assert all(w["transfer_replay"]["matched_held_out_moments"] == 2 for w in result["witnesses"])
    instance = {"spectrum":["1/2"],"diagonal_observations":[[1,1]],"mixed_observations":{"0":[0,1]}}
    assert analyze(instance)["status"] == "MIXED_MOMENTS_FIXED"
    instance["mixed_observations"]["0"] = [1,1]
    assert analyze(instance)["status"] == "EXCLUDED_COMPLEX_MOMENTS"
