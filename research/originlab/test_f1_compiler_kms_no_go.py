"""Exact no-go for a KMS/trace-symmetric full quantum time at D4 x=1/2."""
import unittest
import sympy as sp
from f1_compiler_kms_no_go import controls


class CompilerKMSSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=controls()

    def test_original_D4_two_step_full_quantum_mixture_is_CPTP_and_invertible(self):
        self.assertTrue(self.r["M4_two_step_transfer_invertible"])
        self.assertEqual(self.r["M4_two_step_eigenspace_counts"],
                         {"1":2,"1/4":4,"1/16":10})

    def test_exact_KMS_symmetric_logarithm_requires_an_invalid_Choi_tangent(self):
        self.assertEqual(self.r["CCP_Choi_tangent_witness_norm"],8)
        self.assertEqual(self.r["CCP_Choi_tangent_witness_quadratic_over_log2"],
                         "-16/3")
        self.assertEqual(self.r["CCP_Choi_tangent_negative_eigenvalue_over_log2"],
                         "-2/3")
        self.assertFalse(self.r["trace_symmetric_quantum_markov_embedding_exists"])

    def test_equal_six_source_has_even_stronger_zero_mode_obstruction(self):
        self.assertTrue(self.r["full_S4_democratic_quantum_M4_two_step_has_2_zero_modes"])

    def test_limit_is_strictly_KMS_and_original_source_not_claimed(self):
        self.assertFalse(self.r["non_trace_symmetric_other_logarithms_excluded"])
        self.assertFalse(self.r["physical_original_Padm_and_coherent_amplitudes_evaluated"])
        self.assertFalse(self.r["original_Q0_R113_all_fields_and_OS_time_produced"])
        self.assertFalse(self.r["F1_solved"])


if __name__=="__main__":
    unittest.main()
