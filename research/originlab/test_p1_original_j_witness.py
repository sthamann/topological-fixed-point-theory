"""Exact original v512 j values and monotone-E8-character regression."""
import unittest
from p1_original_j_witness import controls

class P1OriginalJWitnessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=controls()
    def test_exact_nonsquare_original_v512_j(self):
        self.assertEqual(self.data["nonsquare_v512_t_half_j"],"148176/25")
        self.assertEqual(self.data["exact_heldout_j"],{
             "1/3":"48228544/2025","1/2":"148176/25",
             "2/3":"37642192/13689","1":"1728"})
    def test_chiral_character_minimum_on_entire_antipodal_source_family(self):
        self.assertTrue(self.data["j_strictly_decreases_for_t_in_0_to_1"])
        self.assertEqual(self.data["square_j"],"1728")
        self.assertEqual(self.data["square_abs_chi_E8"],"12")
        self.assertTrue(self.data["chiral_E8_magnitude_minimum_not_maximum_at_square"])
    def test_no_native_R113_promotion(self):
        self.assertFalse(self.data["original_weighted_source_factors_as_chiral_E8_proved"])
        self.assertFalse(self.data["original_R113_Z0_ZU_ZY_ZUY_computed"])
        self.assertFalse(self.data["F1_solved"])

if __name__=="__main__":
    unittest.main()
