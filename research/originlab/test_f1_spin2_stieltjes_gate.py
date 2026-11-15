"""Exact native TFPT4.5 Stieltjes TT pole counterexample tests."""
import unittest
import sympy as sp
from f1_spin2_stieltjes_gate import controls


class StieltjesSpin2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=controls()

    def test_extra_TT_pole_exists_with_positive_stieltjes(self):
        self.assertEqual(self.r['new_TT_propagator_pole'],'z=-(m+mu/a)')
        self.assertEqual(self.r['concrete_additional_spin2_mass_squared'],2)

    def test_valid_positive_residues_and_explicit_partial_fractions(self):
        z=sp.symbols('z')
        D=(z+1)/(z*(z+2))
        self.assertEqual(sp.simplify(
            sp.apart(D,z)-(1/(2*z)+1/(2*(z+2)))),0)
        self.assertTrue(self.r['positive_stieltjes_precludes_ghost_sign_for_this_example'])

    def test_zero_free_in_upper_half_not_on_negative_cut(self):
        self.assertFalse(self.r['zero_free_in_open_complex_plane_implies_no_poles_on_negative_axis'])

    def test_no_f1_spin2_physics_claimed(self):
        self.assertFalse(self.r['original_3plus1D_gravity_spin2_proved'])
        self.assertFalse(self.r['F1_solved'])


if __name__=='__main__':
    unittest.main()
