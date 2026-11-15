import sympy as sp

from complex_moments import analyze, control
from mark_transfer import analyze as analyze_transfer
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
    transfer = found["transfer"]
    assert transfer["next_probe"] == {"step":3,"sector":"triplet","interval":["25/128","31/128"]}
    assert transfer["minimum_additional_directions_from_first_three"] == 3
    for witness, third in zip(transfer["sectors"]["triplet"]["spectral_witnesses"], (25, 31)):
        spectral_data = [(sp.Rational(row["transfer_eigenvalue"]),sp.Rational(row["weight"])) for row in witness]
        assert all(0 <= value <= 1 and weight > 0 for value, weight in spectral_data)
        assert [sum(value**k*weight for value, weight in spectral_data) for k in range(4)] == [
            1, sp.Rational(1,2), sp.Rational(5,16), sp.Rational(third,128)]
    replay = run({"cells":{"observations":{"mixed_1_1":"1/2"}}})
    assert replay["cells"]["remaining"] == ["R"]
    pruned = phase_and_clock({"observations":{"clock":{"image_of_point_1":2}}})
    assert pruned["clock"]["remaining_clocks"] == 6


def test_mark_transfer_observations_and_charge_selection():
    p = sp.ones(4)/4
    q = sp.eye(4)-p

    def forms(singlet, triplet):
        return [[list(map(str, row)) for row in (s*p+t*q).tolist()] for s, t in zip(singlet, triplet)]

    data = {"moments":forms([1]*4, [1,sp.Rational(1,2),sp.Rational(5,16),sp.Rational(7,32)])}
    assert analyze_transfer(data)["sectors"]["triplet"]["status"] == "OBSERVED_THIRD_COMPATIBLE"
    assert analyze_transfer({**data,"conserved_distinct_weights":True})["first_defect_step"] == 1
    data["moments"] = forms([1]*4, [1,sp.Rational(1,2),sp.Rational(5,16),sp.Rational(1,4)])
    assert run({"transfer":data})["status"] == "REFUTED_SUPPLIED_MARKS"
    data["moments"] = forms([1]*3, [1,sp.Rational(1,2),sp.Rational(1,8)])
    assert analyze_transfer(data)["sectors"]["triplet"]["status"] == "INFEASIBLE_MOMENTS"
    data["moments"][2][0][1] = 0
    assert analyze_transfer(data)["status"] == "REFUTED_TETRAHEDRAL_FORM"

    moments = [3,sp.Rational(3,2),sp.Rational(3,4)]
    closed = analyze_transfer({"moments":forms(moments,moments),"conserved_distinct_weights":True})
    assert closed["sectors"]["four_weight_lines"]["third_interval"] == ["3/8","3/8"]
    assert closed["minimum_additional_directions_from_first_three"] == 0
    assert closed["next_probe"] is None


def test_mark_transfer_handles_a_null_sector_and_boundary_spectrum():
    p = sp.ones(4)/4
    forms = [p, p/2, p/2]
    result = analyze_transfer({"moments":[[list(map(str,row)) for row in form.tolist()] for form in forms]})
    assert result["sectors"]["triplet"]["spectral_witnesses"] == [[],[]]
    assert result["sectors"]["singlet"]["third_interval"] == ["1/2","1/2"]
    assert result["minimum_additional_directions_from_first_three"] == 1


def test_strictly_positive_transfers_can_differ_at_the_third_step():
    rational = sp.Rational
    spectra = [((rational(1,4), rational(1,2)), (rational(3,4), rational(1,2))),
               ((rational(1,8), rational(4,13)), (rational(2,3), rational(9,13)))]
    moments = [[sum(x**k*w for x,w in spectrum) for k in range(4)] for spectrum in spectra]
    assert moments[0][:3] == moments[1][:3] == [1, rational(1,2), rational(5,16)]
    assert [row[3] for row in moments] == [rational(7,32), rational(79,384)]


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
