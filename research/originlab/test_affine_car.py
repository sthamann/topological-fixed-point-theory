from fractions import Fraction

import affine_car as C
from affine_operators import DOUBLET_ALPHA, MARKED_GAMMA, MARKED_P2, WEAK_ROOT, neg


def test_free_source_reproduces_the_doublet_word():
    result = C.controls(size=4)
    assert result["status"] == "EXACT_FREE_CAR_DOUBLET_SOURCE"
    assert result["agrees_with_lattice_operators"] and not result["marked_word_in_free_source"]
    assert result["triangles_through_marked_gamma"] == {
        "weak_even_pair_currents": 8, "weak_doublet_pair_currents": 4,
        "weak_even_spinor_ends": 8, "weak_doublet_spinor_ends": 8}
    assert result["weighted_r113_values"] is None


def test_canonical_anticommutator_and_spinor_rejection():
    one = C.fermion(False, 0, -1, C.VACUUM)
    assert C.fermion(True, 0, 1, one) == C.VACUUM                    # psi*_{1/2} psi_{-1/2} on the vacuum
    assert C.fermion(False, 0, -1, one) == {}                        # Pauli exclusion
    try:
        C.bilinear(MARKED_P2)
    except ArithmeticError:
        pass
    else:
        raise AssertionError("a spinor current must not be a fermion bilinear")


def test_wrong_weak_assignment_is_detected():
    # A doublet built on the wrong slot pair is a singlet of the weak root e4 - e5.
    wrong = (0, 0, 0, 0, 0, 2, -2, 0)
    assert C.weak_casimir(C.current(wrong, -1, C.VACUUM)) == {}
    doublet = C.current(DOUBLET_ALPHA, -1, C.VACUUM)
    assert C.weak_casimir(doublet) == C.scale(Fraction(3, 4), doublet)
    lowered = C.current(neg(WEAK_ROOT), 0, doublet)
    assert lowered and C.inner(lowered, lowered) == 1
    assert C.weak_casimir(C.current(MARKED_GAMMA, -1, C.VACUUM)) == {}
