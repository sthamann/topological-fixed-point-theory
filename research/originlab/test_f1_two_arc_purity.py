import sympy as sp

from f1_two_arc_purity import controls, guards, ratio, search, x


def test_two_rules_survive_and_are_told_apart():
    result = controls()
    assert result["integer_spin_survivor"]["rule"] == (1, 1, 1)
    assert result["free_seam_survivor"]["rule"] == (1, 0, 0)
    assert result["integer_spin_survivor"]["square"] == "3/2"
    assert result["free_seam_survivor"]["square"] == "2"
    assert result["integer_spin_survivor"]["counterwitness"] == "21/20"
    assert result["free_seam_survivor"]["counterwitness"] == "5/4"


def test_family_between_the_two_is_discrete():
    rows = search(4)
    pairs = sorted(int(row["fermion_pairs_exchanged"]) for row in rows)
    assert pairs[0] == 0 and pairs[-1] == 16
    assert all(row["rule"][1] == row["rule"][2] for row in rows)


def test_unequal_side_pairings_fail_decoupling():
    assert guards((1, 1, 0)) == {"distant_arcs_decouple": False}
    assert guards((0, 1, 2)) == {"distant_arcs_decouple": False}


def test_e8_has_no_fermion_exchange():
    assert sp.series(ratio((1, 1, 1)), x, 0, 2).removeO() == 1


def test_original_not_proved():
    result = controls()
    assert result["original_r113_values"] is None
    assert result["physical_source_identity_proved"] is False


def test_mark_types_and_reflections():
    from f1_two_arc_purity import reflections
    data = reflections()
    assert data["weights"] == {"1": "0", "S": "1", "P": "1", "PS": "1/2"}
    assert data["U_is_S_and_Y_is_P_on_equal_copies"].startswith("EXCLUDED")


def test_recorded_canonical_model_is_the_two_copy_algebra():
    from coherent_origin import marked_operators
    swap, phase, parity = marked_operators()
    exchange, sign = sp.Matrix([[0, 1], [1, 0]]), sp.diag(1, -1)
    assert swap[:2, :2] == exchange and parity[:2, :2] == sign
    assert phase[:2, :2] == sp.I*exchange*sign


def test_one_charge_function_for_flat_sector_and_quarter_law():
    from f1_two_arc_purity import twist_rule
    data = twist_rule()
    assert data["marks"]["PS"]["ground_states_per_fermion"] == [1, 1]
    assert data["marks"]["P"]["ground_states_per_fermion"] == [2, 1]
    assert data["matter_half"] == 128


def test_charged_seam_fixes_sixteen_fermions_at_the_quarter_angle():
    from f1_two_arc_purity import charged_seam
    data = charged_seam()
    assert data["quarter_angle_solution"] == {"fermions": 16, "angle_over_2pi": "1/4"}
    assert {row["fermions"] for row in data["weight_half_solutions"]} == {4, 9, 16, 25, 36, 49, 64}
