"""Exact tests of the conditional S4 covariant oriented-seam CP source."""
import unittest

from f1_oriented_seam_source import controls


class OrientedSeamSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=controls()

    def test_full_M4_source_is_strict_GKLS_with_unique_scalar_fixed_point(self):
        self.assertEqual(self.r["full_M4_generator_spectrum"],
                         {"0":1,"-gamma/2":6,"-gamma":9})
        self.assertEqual(self.r["exact_CCPositivity_tangent_spectrum"],
                         {"1":3,"0":13})
        self.assertEqual(self.r["unique_stationary_density"],
                         "P0=J4/4, rank one")

    def test_single_S4_covariant_jump_scale_and_sewing(self):
        self.assertTrue(self.r["finite_time_CPTP_and_exact_sewing"])
        self.assertTrue(self.r["full_S4_covariance"])
        self.assertEqual(
            self.r["S4_orthogonal_sector_schur"]["commutant_linear_constraint_rank"],8)
        self.assertEqual(
            self.r["S4_orthogonal_sector_schur"]["covariant_positive_jump_rate_free_scales"],1)

    def test_minimal_Kraus_histories_saturate_and_do_not_create_long_time(self):
        self.assertEqual(
            self.r["minimal_history_register"]["one_step_minimal_environment_dimension"],6)
        self.assertEqual(
            self.r["minimal_history_register"]["two_step_minimal_environment_dimension"],10)
        self.assertEqual(
            self.r["minimal_history_register"]["three_step_minimal_environment_dimension"],10)
        self.assertFalse(self.r["minimal_history_register"][
            "post_second_step_unbounded_finite_matrix_history_growth"])

    def test_seam_jump_not_in_original_six_cycle_operator_span(self):
        self.assertTrue(self.r["history_permutations_cannot_mix_P0_and_Q"])
        self.assertFalse(self.r["faithful_tracial_KMS_vacuum"])
        self.assertFalse(self.r["original_R113_source_selected_jumps_or_arrow"])
        self.assertFalse(self.r["F1_solved"])


if __name__=="__main__":
    unittest.main()
