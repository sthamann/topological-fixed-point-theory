"""Independent exact original TFPT source theorem gate regression."""
import unittest
import sympy as sp
from f1_source_theorem_audit import (
    controls, uv_rg_contract, blocked_dirac_and_reflection,
    bernstein_moment_nonuniqueness, integrable_kernel_repair)

class F1OriginalSourceTheoremTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report=controls()

    def test_uniform_kernel_on_unbounded_t_is_not_finite_RG(self):
        r=self.report["original_UV_fixed_point"]
        self.assertFalse(r["unique_finite_fixed_measure_derived_from_as_stated_assumptions"])
        self.assertIn("infinity",r["total_mass_T_nu"])

    def test_original_Birkhoff_diameter_bound_is_too_small(self):
        r=self.report["original_UV_fixed_point"]
        self.assertTrue(r["published_contraction_bound_falsified"])
        self.assertEqual(r["actual_projective_contraction"],"1/3")
        self.assertEqual(r["actual_valid_uniform_diameter_bound"],"2*log(M/m)")
        self.assertTrue(sp.Rational(1,3)>3-2*sp.sqrt(2))

    def test_positive_normalized_K_does_not_preserve_branch_moments(self):
        r=self.report["original_UV_fixed_point"]
        self.assertFalse(r["positive_normalized_kernel_preserves_prescribed_low_moments"])
        self.assertEqual(r["branch_moment_counterexample"],
                         "mean 1/2 changes to 1/3")

    def test_exact_fine_Dirac_compression_but_no_measure_pushforward(self):
        r=self.report["spectral_measure_and_OS"]
        self.assertTrue(r["fine_Dirac_compresses_exactly_to_coarse_for_every_y"])
        self.assertFalse(r["spectral_weight_projective_consistency_from_block_Dirac"])
        self.assertEqual(r["fine_trace_actions"],
                         {"0,0":"2","0,1":"1","1,0":"3/2","1,1":"1"})
        self.assertEqual(r["coarse_trace_actions"],["1","1/2"])

    def test_pointwise_positive_spectral_weight_fails_OS(self):
        r=self.report["spectral_measure_and_OS"]
        self.assertFalse(r["OS_positivity_from_positive_pointwise_spectral_weight"])
        q=sp.exp(-sp.Rational(1,2))
        self.assertTrue(q*q-q<0)
        self.assertTrue((sp.Matrix([[q*q,q],[q,q*q]])).det()<0)

    def test_two_strictly_positive_measures_share_four_moments(self):
        r=self.report["finite_moment_nonidentifiability"]
        self.assertTrue(r["both_strict_positive_and_normalized"])
        self.assertTrue(r["moments_k0_k1_k2_k3_equal"])
        self.assertFalse(r["first_four_moments_derive_unique_UV_profile"])
        self.assertEqual(r["complete_Laplace_profiles_differ"],
                         "(1-exp(-t))^4/30")
        self.assertTrue(r["both_normalize_mass_and_have_unique_stationary_measures"])
        self.assertTrue(r["both_preserve_the_same_four_prescribed_moments"])
        self.assertTrue(r["both_contract_projectively_to_one_ray_in_one_step"])
        self.assertTrue(r["still_distinct_complete_Bernstein_profiles"])

    def test_integrable_RG_repair_has_valid_unique_fixed_profile(self):
        r=self.report["repaired_positive_RG_comparison"]
        self.assertEqual(r["moment_fixed_point"],"u*=(sqrt(13)-2)/3")
        self.assertEqual(r["moment_map_contraction"],
                         "|F'(u)|<=1/12 on u in [0,1]")
        self.assertTrue(r["different_valid_integrable_kernels_have_different_fixed_profiles"])
        self.assertFalse(r["native_physical_RG_kernel_derived"])

    def test_no_original_R113_or_F1_scientific_promotion(self):
        self.assertFalse(self.report["F1_solved"])
        self.assertFalse(self.report["original_source_Q0_Z0_ZU_ZY_ZUY_computed"])
        self.assertFalse(self.report["repaired_positive_RG_comparison"]
                         ["native_branch_moment_slice_compatibility_proved"])

if __name__=="__main__":
    unittest.main()
