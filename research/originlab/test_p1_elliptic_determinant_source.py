"""Source-conditional elliptic double-cover determinant and E8 Hessian tests."""
import unittest
from p1_elliptic_determinant_source import controls

class P1EllipticDeterminantSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=controls()

    def test_original_v512_fourmark_j_uniformizes_rectangular_elliptic_cover(self):
        d=self.data
        self.assertTrue(d["native_v512_j_equals_elliptic_j_of_x"])
        self.assertEqual(d["numeric_samples"]["1"]["rectangular_modulus_y"],"1.0")
        self.assertEqual(d["numeric_samples"]["1/2"]["elliptic_lambda"],"0.2")

    def test_unit_area_scalar_determinant_square_maximizes_rectangular_family(self):
        d=self.data
        self.assertTrue(d["positive_det_power_maximizes_square_in_rectangular_class"])
        self.assertFalse(d["free_boson_negative_half_det_power_maximizes_square"])
        self.assertLess(float(d["numeric_samples"]["1/2"]["unit_area_det"]),
                        float(d["numeric_samples"]["1"]["unit_area_det"]))

    def test_E8_chiral_squared_competition_requires_at_least_15_positive_powers(self):
        d=self.data
        self.assertEqual(d["conditional_chiral_E8_times_D_power_min_integer_for_local_square_max"],15)
        self.assertTrue(14.75<float(d["conditional_real_power_threshold"])<14.76)
        self.assertLess(float(d["det_log_second_y_derivative_at_square"]),-1.29)
        self.assertGreater(float(d["chiral_E8_abs_character_squared_log_second_y_derivative_at_square"]),19)

    def test_not_a_native_physical_source_or_spin2_parent(self):
        d=self.data
        self.assertFalse(d["source_physical_elliptic_Dirac_Laplacian_identified"])
        self.assertFalse(d["source_net_ghost_fermion_boson_exponent_identified"])
        self.assertFalse(d["source_E8_character_factorization_proved"])
        self.assertFalse(d["native_Q0_and_R113_produced"])
        self.assertFalse(d["F1_solved"])

if __name__=="__main__":
    unittest.main()
