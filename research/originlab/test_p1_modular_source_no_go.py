"""Regression of chiral E8 square saddle and synthetic positive modulus controls."""
import unittest
from p1_modular_source_no_go import controls

class P1ModularSourceNoGoTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=controls()

    def test_exact_E8_character_square_is_twelve_but_not_maximum(self):
        d=self.data["E8_chiral_modular"]
        self.assertEqual(d["chi_E8_at_tau_i"],"12")
        self.assertEqual(d["E6_at_tau_i"],"0")
        self.assertFalse(d["chiral_E8_square_selected_by_maximizing_character"])
        self.assertLess(float(d["abs_chi_i_plus_0p01"]),12)
        self.assertGreater(float(d["abs_chi_1p01i"]),12)

    def test_explicit_modular_hessian_is_saddle(self):
        d=self.data["E8_chiral_modular"]
        self.assertEqual(d["log_chi_holomorphic_second_derivative_at_i"],
                         "-(2*pi^2/3)*E4(i)")
        self.assertTrue(float(d["conditional_rest_curvature_threshold_numerical"])>19)
        self.assertFalse(d["actual_TFPT_factorization_into_E8_character_times_rest_proved"])

    def test_both_mirror_positive_OS_sources_disagree_on_square(self):
        d=self.data["positive_modulus_countermodels"]
        self.assertFalse(d["symmetry_reflection_positivity_and_sewing_select_square"])
        m=d["source_models"]
        self.assertEqual(m["selects_square"]["square_log_weight_curvature"],"-2/3")
        self.assertEqual(m["rejects_square"]["square_log_weight_curvature"],"2")
        self.assertEqual(m["selects_square"]["w_half_over_w_one"],"22/25")
        self.assertEqual(m["rejects_square"]["w_half_over_w_one"],"34/25")
        self.assertEqual(m["selects_square"]["reflection_positive_time_Gram_rank"],2)
        self.assertEqual(m["rejects_square"]["reflection_positive_time_Gram_rank"],2)

    def test_original_source_not_derived(self):
        self.assertFalse(self.data["original_blocked_fiber_packet_produced"])
        self.assertFalse(self.data["original_R113_Z0_ZU_ZY_ZUY_produced"])
        self.assertFalse(self.data["F1_solved"])

if __name__=="__main__":
    unittest.main()
