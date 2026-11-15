"""Regression for the original 5+3 geometric seam Dirac candidate."""
import unittest
import sympy as sp
from p1_seam_dirac import controls


class SeamDiracTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=controls()

    def test_full_D4_SU8_determinant_cancellation(self):
        self.assertEqual(self.result["geometric_D4_rotation_determinants"],
                         ["-1","-1","1"])
        self.assertEqual(self.result["geometric_D4_reflection_determinants"],
                         ["-1","-1","1"])
        self.assertTrue(self.result["complex_D4_action_lifts_to_SU8"])

    def test_exact_three_pairs_plus_two_endpoint_zeros(self):
        self.assertEqual(self.result["original_five_plus_three_mode_rank"],8)
        self.assertEqual(self.result["canonical_inclusion_rank"],3)
        self.assertEqual(self.result["unreduced_geometric_index"],2)
        self.assertEqual(self.result["Dirac_zero_modes"],
                         "two at geometric quarter fixed points 0 and infinity")

    def test_mass_spectral_flow_and_reflection(self):
        m=sp.symbols("m",real=True)
        self.assertTrue(self.result["rotation_invariant"])
        self.assertEqual(self.result["reflection_odd_mass_intertwiner"],
                         "S H(m) S = -H(-m)")
        self.assertEqual(self.result["Dirac_mass_determinant"],"-m^2(1+m^2)^3")
        self.assertEqual(self.result["one_selected_endpoint_index"],1)

    def test_one_sided_positive_time_has_unique_vacuum(self):
        self.assertEqual(self.result["one_sided_positive_time_H_squared_spectrum"],
                         ["0 (x1)","1 (x6)"])
        self.assertEqual(self.result["one_sided_vacuum_rank"],1)
        self.assertEqual(self.result["one_sided_transfer_formula"],
                         "T(t)=P_vac+exp(-t)*(I-P_vac)")
        self.assertEqual(self.result["one_sided_two_branch_intertwinement"],
                         "S P_0 S=P_infinity")
        self.assertFalse(self.result["one_sided_cap_selection_and_time_scale_derived"])

    def test_source_and_time_not_proven(self):
        self.assertFalse(self.result["single_endpoint_selection_from_original_P1_proved"])
        self.assertFalse(self.result["native_Dirac_selfadjoint_identified"])
        self.assertFalse(self.result["physically_independent_8_CAR_sources_proved"])
        self.assertFalse(self.result["physical_nonnegative_time_generator_derived"])
        self.assertFalse(self.result["F1_solved"])


if __name__=="__main__":
    unittest.main()
