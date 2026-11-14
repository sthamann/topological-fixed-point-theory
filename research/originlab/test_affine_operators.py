from fractions import Fraction

import pytest

import affine_operators as A
from affine_time import charge_signature as original_signature, marked_charges
import sympy as sp


@pytest.fixture(scope="module")
def receipt():
    return A.controls(word_size=4, algebra_roots=8)


def test_charge_signature_matches_the_original_chart():
    for row in marked_charges():
        weight = sp.Matrix([sp.Rational(x) for x in row["root"]])
        doubled = tuple(int(2*x) for x in weight)
        glue, triality, hyper = A.charge_signature(doubled)
        assert (glue, triality, str(hyper)) == original_signature(weight)


def test_operator_source_reproduces_the_affine_word(receipt):
    assert receipt["status"] == "EXACT_LEVEL_ONE_OPERATOR_SOURCE"
    assert receipt["level"] == 1 and receipt["current_relation_checks"] > 1000
    assert receipt["neutral_return_spectra"]["2"] == {"spectrum": {"1": 43, "4": 1}, "dimension": 44}
    assert receipt["weighted_r113_values"] is None and receipt["physical_source_identity_proved"] is False


def test_joint_return_extends_the_single_root_line(receipt):
    joint = receipt["joint_neutral_return"]
    assert joint["orthogonal_frame_commutes"] and joint["orthogonal_frame_excess_rank"] == 8
    assert joint["marked_plane_excess_rank"] == 3
    assert {(row["pairing"], row["commutator_rank"]) for row in joint["marked_plane_commutators"]} == {(1, 2), (-1, 0)}


def test_marked_plane_sees_only_a2_and_its_commutant_is_e6(receipt):
    assert receipt["marked_plane_system"] == "A2" and receipt["commutant_system"] == "E6"
    assert receipt["grade_one_vectors_generated_by_marked_currents"] == 8


def test_completion_needs_six_currents_and_one_weak_doublet(receipt):
    searches = receipt["mark_completion"]
    assert searches["centre_neutral_d8_pairs"]["minimal_extra_currents"] == 6
    assert receipt["blocked_completion_classes"] == [
        "centre_and_glue_neutral", "centre_and_hypercharge_neutral", "hypercharge_neutral"]
    for name in receipt["blocked_completion_classes"]:
        assert searches[name]["largest_generated_system"] == ["A1", "E7"]
        assert searches[name]["conserved_parities"] == [list(A.WEAK_ROOT)]
    for row in searches.values():
        if row["reaches_e8"]:
            assert any(A.dot(tuple(r), A.WEAK_ROOT) % 2 for r in row["witness"])
    assert receipt["weak_slot_parity"]["colour_neutral_pair_doublets_glue_and_abs_hypercharge"] == [[2, "1/2"]]


def test_five_currents_never_complete_the_marked_plane():
    # Negative control for the rank bound: a seven-dimensional span cannot contain all roots.
    marked = (A.MARKED_GAMMA, A.MARKED_P2)
    witness = A.completion_search(marked, A.ROOTS)["witness"]
    system = A.closure({*marked, *(tuple(r) for r in witness[:5])})
    assert A.rank(sorted(system)) == 7 and len(system) < 240


def test_wrong_kernel_and_wrong_level_are_rejected():
    frame, _ = A.lattice_frame()
    fock = A.LatticeFock([[A.dot(a, b) for b in frame] for a in frame])
    word = A.marked_word(fock, frame, 3)
    assert [abs(word[m, n]) for m in (1, 2, 3) for n in (1, 2, 3)] == [1, 1, 1, 1, 2, 2, 1, 2, 3]
    assert any(abs(word[m, n]) != m*n for m, n in word)           # not a product of separate factors
    gamma = A.integral(frame, A.MARKED_GAMMA)
    vector = fock.vertex(gamma, -2, fock.vacuum())
    assert fock.inner(vector, vector) == Fraction(2)               # level one, not the free-square minimum 30
    doubled = A.LatticeFock([[2*x for x in row] for row in fock.gram])
    assert doubled.vertex(gamma, -1, doubled.vacuum()) != {}       # weight two field, no current algebra
    assert A.current_algebra(doubled, frame, sorted(A.closure((A.MARKED_P1, A.MARKED_P2))), (-1, 0, 1),
                             [doubled.vacuum()]) is None
