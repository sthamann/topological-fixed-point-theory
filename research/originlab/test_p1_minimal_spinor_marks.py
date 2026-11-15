"""Exact n-mark Riemann-Roch/Spin(4n) local simple-current selection."""
import unittest
from p1_minimal_spinor_marks import puncture_class, e8_weight_one_root_census, controls


class MinimalSpinorMarksTests(unittest.TestCase):
    def test_all_n_two_to_thirty_two_are_classified_exactly(self):
        result=controls()
        self.assertEqual(result["integer_local_extension_mark_counts_to_bound"],
                         list(range(4,33,4)))
        self.assertEqual(result["unique_n_with_new_spinor_weight_one_currents"],4)

    def test_three_and_five_arise_without_positing_E8(self):
        c=puncture_class(4)
        self.assertEqual((c["meromorphic_sections"],c["logarithmic_differentials"]),
                         (5,3))
        self.assertEqual((c["candidate_real_majoranas"],c["chiral_central_charge"]),
                         (16,8))
        self.assertEqual(c["spinor_conformal_weight"],"1")

    def test_lower_mark_counts_cannot_have_bosonic_spinor_weight_one(self):
        self.assertEqual([puncture_class(n)["spinor_conformal_weight"]
                          for n in (2,3,4,5)],["1/2","3/4","1","5/4"])
        for n in (2,3):
            self.assertFalse(puncture_class(n)["integral_spin_simple_current_allowed"])

    def test_higher_integral_spin_does_not_produce_new_weight_one_currents(self):
        c=puncture_class(8)
        self.assertTrue(c["integral_spin_simple_current_allowed"])
        self.assertFalse(c["new_spinor_weight_one_currents"])
        self.assertEqual(c["SO_bilinear_current_dimension"],496)

    def test_native_E8_root_census_is_emergent_at_four_marks(self):
        self.assertEqual(e8_weight_one_root_census(),(112,128,240))
        self.assertEqual(controls()["at_n4"]["full_E8_currents_including_Cartan"],248)

    def test_physical_source_assumptions_are_not_disguised_as_theorems(self):
        r=controls()
        self.assertFalse(r["candidate_original_source_independence_of_2n_CAR_modes_proved"])
        self.assertFalse(r["physical_spinor_extension_locality_and_GSO_proved"])
        self.assertFalse(r["raw_P1_square_modulus_selected"])
        self.assertFalse(r["original_weighted_R113_produced"])
        self.assertFalse(r["F1_solved"])


if __name__=="__main__":
    unittest.main()
