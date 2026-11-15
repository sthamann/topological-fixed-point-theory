import sympy as sp

from mark_transfer import analyze
from square_search import projective_stabilizer, reciprocal_bridge_sweep, run


def test_harmonic_seam_selectors_and_bounded_certificate():
    assert len(projective_stabilizer(-1, 1)) == 8
    assert len(projective_stabilizer(-1, 4)) == 4
    report = run({"max_denominator": 5, "transfer_grid": 8, "odd_green_ratio": "2",
                  "adjacent": "1/8"})
    assert report["moduli"]["projective_stabilizers"] == {"V4": 9, "D4": 1}
    assert report["moduli"]["quarter_survivors_t"] == ["1"]
    assert report["moduli"]["odd_green_ratio_survivors_t"] == ["1/2"]
    assert report["d4_transfers"]["matches_supplied_orbit_observations"] == 9
    unique = run({"max_denominator": 1, "transfer_grid": 8, "adjacent": "1/8", "opposite": "0"})
    assert unique["d4_transfers"]["matching_example_eigenvalues"] == [["1/2", "3/4"]]


def forms_for_spectrum(spectrum):
    characters = sp.Matrix([[1, 1, 1, 1], [1, -1, 1, -1],
                            [1, 1, -1, -1], [1, -1, -1, 1]])
    forms = [characters.T*sp.diag(*(x**step for x in spectrum))*characters/4 for step in range(4)]
    return [[list(map(str, row)) for row in form.tolist()] for form in forms]


def test_square_symmetry_preserves_two_independent_rates():
    data = {"moments": forms_for_spectrum([1, sp.Rational(1, 2), sp.Rational(3, 4), sp.Rational(3, 4)])}
    result = analyze({**data, "symmetry": "D4"})
    assert result["status"] == "EXACT_FOUR_MARK_TRANSFER_BOUNDS"
    assert result["sectors"]["alternating"]["third_interval"] == ["1/8", "1/8"]
    assert result["sectors"]["doublet"]["third_interval"] == ["27/64", "27/64"]
    assert analyze(data)["status"] == "REFUTED_TETRAHEDRAL_FORM"


def test_off_square_odd_doublet_can_split():
    forms = forms_for_spectrum([1, sp.Rational(1, 2), sp.Rational(3, 4), sp.Rational(1, 4)])
    assert analyze({"moments": forms, "symmetry": "D4"})["status"] == "REFUTED_SQUARE_FORM"
    result = analyze({"moments": forms, "symmetry": "V4"})
    assert result["sectors"]["odd_1"]["third_interval"] == ["27/64", "27/64"]
    assert result["sectors"]["odd_2"]["third_interval"] == ["1/64", "1/64"]


def test_reciprocal_bridge_needs_a_withheld_probe():
    report = reciprocal_bridge_sweep({"factor_limit": 5, "max_factors": 2, "max_power": 1,
                                      "fit": {"1/2": "16/25"}, "withheld": {"1/3": "25/49"}})
    assert report["fit_survivors"] >= 2
    assert report["withheld_survivors"] == 1
    assert any(model["factors"] == ["34/7"] and not model["passes_withheld"]
               for model in report["fit_survivor_examples"])
