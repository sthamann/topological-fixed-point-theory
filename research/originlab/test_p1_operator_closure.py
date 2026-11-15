"""Exact three-state matrix closure and native affine current-level exclusion."""
import unittest
from p1_operator_closure import controls


class P1OperatorClosureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report=controls()

    def test_full_matrix_closure_selects_exactly_four_marks(self):
        r=self.report
        self.assertEqual(r["unique_n_at_least_2"],[4])
        self.assertEqual(r["n4_operator_algebra_dim"],9)
        self.assertEqual(r["spin_one_operator_decomposition"],[1,3,5])
        self.assertEqual(r["matrix_associative_products_checked"],81)

    def test_compact_su3_lie_and_positive_trace_metric(self):
        r=self.report
        self.assertEqual(r["traceless_compact_Lie_algebra"],"su3")
        self.assertEqual(r["positive_Hilbert_Schmidt_metric_spectrum"],
                         {"3":1,"2":7,"6":1})

    def test_eight_complex_adjoint_fermions_are_su3_level_six_not_native_one(self):
        r=self.report
        self.assertEqual(r["adjoint_trace_index_over_fundamental"],6)
        self.assertEqual(r["affine_level_if_8_COMPLEX_chiral_fermions_in_su3_adjoint"],6)
        self.assertEqual(r["native_regular_E8_root_su3_level"],1)
        self.assertFalse(r["can_identify_both_as_the_same_level1_native_colour_current"])

    def test_source_not_conflated_with_mathematical_matrix_model(self):
        r=self.report
        self.assertFalse(r["native_positive_weighted_source_and_R113_identified"])
        self.assertFalse(r["physical_principle_all_operators_from_original_sections_derived"])
        self.assertFalse(r["F1_solved"])


if __name__=="__main__":
    unittest.main()
