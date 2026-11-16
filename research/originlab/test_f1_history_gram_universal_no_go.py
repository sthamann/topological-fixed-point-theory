"""Independent exact regressions for six-history source Gram and all-x CCP no-go."""
import unittest
from f1_history_gram_universal_no_go import controls


class HistorySourceKernelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report=controls()

    def test_all_D4_weights_excluded_for_tracial_KMS_full_M4_source(self):
        a=self.report["all_weights_no_go"]
        self.assertTrue(a["all_x_in_closed_interval_excluded_in_tracial_KMS_class"])
        self.assertEqual(a["singular_source_weights"],["0","1/3","1"])
        self.assertEqual(a["M4_two_step_generic_eigenvalues"],
                         {"1":2,"x^2":4,"(1-x)^2/4":8,"(3*x-1)^2/4":2})
        self.assertFalse(a["nontracial_or_nonselfadjoint_or_dilated_sources_excluded"])

    def test_Gram_kernel_dimension_and_CP_TP_open_family(self):
        g=self.report["history_kernel"]
        self.assertEqual(g["D4_Hermitian_Gram_real_dimension"],7)
        self.assertEqual(g["trace_preserving_constraint_rank"],3)
        self.assertEqual(g["remaining_affine_positive_TP_dimensions"],4)
        self.assertTrue(g["six_cycle_Kraus_vectors_independent"])

    def test_only_one_untwisted_D4_symmetric_coherent_clock(self):
        c=self.report["coherent_symmetry"]
        self.assertEqual(c["relative_quarter_phase_assignments"],1024)
        self.assertEqual(c["unitary_history_sums"],9)
        self.assertEqual(c["untwisted_D4_invariant_rank_one_Grams"],1)
        self.assertFalse(c["projective_native_spin_action_evaluated"])

    def test_supplied_16D_spectra_are_not_a_quantum_CP_certificate(self):
        s=self.report["spectral_typing"]
        self.assertTrue(s["supplied_spectrum_matched_by_T_tensor_involution"])
        self.assertFalse(s["that_spectrum_implies_CP_or_quantum_time"])

    def test_physical_F1_not_fabricated(self):
        self.assertFalse(self.report["original_R113_full_source_bundle_available"])
        self.assertFalse(self.report["F1_solved"])


if __name__=="__main__":
    unittest.main()
