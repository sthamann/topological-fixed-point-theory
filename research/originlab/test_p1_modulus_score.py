"""Original-source modulus stress-versus-fluctuation exact decision tests."""
import unittest
from p1_modulus_score import controls

class P1ModulusScoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=controls()

    def test_exact_source_second_variation(self):
        self.assertEqual(self.data["exact_second_derivative"],
                         "d2log(Z)/dx2=Var(Sprime)-E(Ssecond)")
        self.assertEqual(self.data["square_curvature"],"a^2-c")
        self.assertEqual(self.data["necessary_square_maximum"],
                         "E(Ssecond)>Var(Sprime)")

    def test_positive_mirror_sources_give_opposite_square_selection(self):
        self.assertEqual(self.data["positive_source_counterexamples"],{
          "square_maximum":"-1","square_minimum":"2",
          "square_flat_to_second_order":"0"})
        self.assertFalse(self.data["mirror_only_selects_square"])

    def test_no_missing_original_terms_got_promoted(self):
        self.assertFalse(self.data["original_modulus_stress_fluctuations_computed"])
        self.assertFalse(self.data["original_moving_domain_projector_and_eta_variations_included"])
        self.assertFalse(self.data["original_R113_four_values_produced"])
        self.assertFalse(self.data["F1_solved"])

if __name__=="__main__":
    unittest.main()
