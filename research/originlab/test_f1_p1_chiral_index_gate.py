"""Exact original four-mark log zero-mode index and source/flux obligations."""
import unittest
import sympy as sp
from f1_p1_chiral_index_gate import controls, line_bundle_cohomology


class P1ChiralIndexTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=controls()

    def test_original_four_marks_yield_log_bundle_O2(self):
        self.assertEqual(self.r["original_P1_canonical_line_bundle"],"K=O(-2)")
        self.assertEqual(self.r["logarithmic_one_form_line_bundle"],"K(D)=O(2)")
        self.assertEqual(self.r["logarithmic_one_form_cohomology"],
                         {"degree":2,"h0":3,"h1":0,"index":3})

    def test_spin_twist_O3_yields_three_one_chirality_Dirac_zero_modes(self):
        self.assertEqual(self.r["positive_spin_bundle"],"S+=O(-1)")
        self.assertEqual(self.r["twisting_line_bundle_for_same_log_spinor"],"E=O(3)")
        self.assertEqual(self.r["twisted_two_dimensional_spin_Dirac_index"],3)
        self.assertEqual(self.r["positive_flux_zero_modes"],3)
        self.assertEqual(self.r["negative_chirality_zero_modes_for_positive_flux"],0)

    def test_opposite_monopole_flux_changes_chiral_index(self):
        self.assertEqual(self.r["opposite_flux_O_minus_three_index"],-3)
        self.assertEqual(line_bundle_cohomology(-4),
                         {"degree":-4,"h0":0,"h1":3,"index":-3})

    def test_family_index_is_topologically_stable_but_does_not_select_square(self):
        self.assertTrue(self.r["index_stable_under_moving_four_distinct_marks"])
        self.assertEqual(self.r["n_mark_index_formula"],
                         "ind(dbar_K(D))=n-1 on P1")
        self.assertFalse(self.r["square_modulus_from_index_selected"])

    def test_requires_physical_parent_and_does_not_claim_three_SM_generations(self):
        self.assertFalse(self.r["P1_log_forms_source_to_actual_fermion_field_map_proved"])
        self.assertFalse(self.r["higher_dimensional_chiral_Dirac_parent_built"])
        self.assertFalse(self.r["original_E8_R_and_R_conjugate_chirality_index_verified"])
        self.assertFalse(self.r["three_physical_3plus1D_chiral_families_proved"])
        self.assertFalse(self.r["original_R113_Q0_os_source_identified"])
        self.assertFalse(self.r["F1_solved"])


if __name__=="__main__":
    unittest.main()
