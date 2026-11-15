"""Exact fixed-divisor no-go, moving PSL2 lift and primitive n=4 Lie gate."""
import unittest
from p1_moving_divisor_origin import controls, polynomial_sl2


class MovingDivisorOriginTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report=controls()

    def test_original_regional_mobius_requires_moving_marks(self):
        r=self.report
        self.assertTrue(r["original_four_regional_fields_span_sl2"])
        self.assertEqual(r["fixed_D_naive_Lie_derivative_double_pole"],
                         "-1/4 at z=1 for translation of f0")
        self.assertEqual(r["fixed_four_point_SL2_invariant_vector_field_dimension"],0)
        self.assertEqual(r["moving_PSL2_cross_ratio_derivatives"],[0,0,0])

    def test_correction_closes_five_and_three_mode_polynomial_transport(self):
        self.assertEqual(self.report["moving_divisor_numerator_corrections_for_1_z_z2"],
                         ["0","-4","-4*z"])
        for n in (2,4):
            E,F,H,G=polynomial_sl2(n)
            self.assertEqual(E*F-F*E,2*H)
            self.assertEqual(E.T*G,-G*F)

    def test_old_D4_exact_sequence_not_continuous_sl2_equivariant(self):
        r=self.report
        self.assertFalse(r["canonical_F3_to_C5_inclusion_full_sl2_equivariant"])
        self.assertFalse(r["f2_compensator_invariant_line_under_full_regional_sl2"])
        self.assertTrue(r["original_four_regional_fields_span_sl2"])

    def test_no_nontrivial_positive_finite_fixed_noncompact_sl2(self):
        self.assertFalse(self.report["positive_nontrivial_finite_fixed_SL2R_unitary"])
        self.assertEqual(self.report["C5_positive_compact_SU2_metric_diagonal"],
                         ["1","1/4","1/6","1/4","1"])

    def test_quartic_to_quadratic_third_transvectant_and_eight_point_jacobi(self):
        r=self.report
        self.assertEqual(r["third_transvectant_sl2_equivariance_checks"],75)
        self.assertEqual(r["nonzero_quartic_binary_third_transvectant"],"576*x*y")
        self.assertEqual(r["n8_unique_wedge_bracket_Jacobi_failure"],"5184000*x**6")
        self.assertEqual(r["nontrivial_sl2_equivariant_internal_Lie_on_log_Fn_forces_n"],4)

    def test_unital_operator_completion_unique_four_without_physical_promotion(self):
        self.assertEqual(self.report["unital_End_Fn_equals_scalar_plus_Fn_plus_Cn_for_n"],[4])
        self.assertFalse(self.report["physical_moving_divisor_connection_identified"])
        self.assertFalse(self.report["native_weighted_Q0_R113_produced"])
        self.assertFalse(self.report["F1_solved"])


if __name__=="__main__":
    unittest.main()
