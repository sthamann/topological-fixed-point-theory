"""Exact regression suite: pure P1 constraint state vs native projective consistency."""
import unittest

from f1_seam_projective_consistency import controls


class P1ProjectiveConsistencyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=controls()

    def test_pure_GHZ_seam_stages_cannot_be_ordinary_quantum_marginals(self):
        z=self.r["ordinary_coarse_graining"]
        self.assertEqual(z["d4_source_discrepancy_trace_distance"],"3/4")
        for d in (2,3,4):
            t=z["dimensions"][str(d)]
            self.assertEqual(t["seam_pure_rank"],1)
            self.assertEqual(t["trace_fine_to_coarse_rank"],d)
            self.assertFalse(t["ordinary_partial_trace_projective_consistent"])
            self.assertTrue(t["dephased_mixtures_partial_trace_consistent"])

    def test_original_cyclic_invariant_joint_mark_detects_marginal_failure(self):
        z=self.r["ordinary_coarse_graining"]
        self.assertTrue(z["observable_joint_shift_commutes_with_constraint"])
        self.assertEqual(z["d4_shift_invariant_joint_probe_pure"],"1")
        self.assertEqual(z["d4_shift_invariant_joint_probe_partial_trace"],"0")

    def test_all_four_quarter_characters_survive_constraints_without_physical_phase(self):
        z=self.r["time_character_selection"]
        self.assertEqual(z["Z4_character_sector_count"],4)
        self.assertEqual(z["joint_rank_each"],1)
        self.assertFalse(z["neutral_character_source_physically_selected"])
        self.assertFalse(z["nontrivial_quarter_spin_source_physically_selected"])

    def test_cap_norm_versus_valid_CPTP_nonlocal_decoder(self):
        z=self.r["physical_cap_vs_categorical_gluing"]
        self.assertEqual(z["physical_normalized_cap_single_step_success"],"1/4")
        self.assertTrue(z["nonlocal_CPTP_decoder_on_fine_algebra_constructed"])
        self.assertTrue(z["copy_map_delta_isometric_and_coassociative"])
        self.assertFalse(z["nonlocal_decoder_preserves_old_local_AB_marginal"])
        self.assertFalse(z["unnormalized_categorical_counit_is_physical_effect"])

    def test_partial_which_history_information_dephases_source(self):
        z=self.r["source_overlap_tradeoff"]
        self.assertEqual(z["coarse_2state_fidelity"],"(1+3lambda)/4")
        self.assertEqual(z["coarse_2state_trace_distance"],"3(1-lambda)/4")
        self.assertEqual(z["connected_label_correlation_at_any_distance"],1)

    def test_source_boundaries_remain_open(self):
        self.assertFalse(self.r["native_original_R113_full_Q0_available"])
        self.assertFalse(self.r["actual_Padm_and_absolute_Z0_evaluated"])
        self.assertFalse(self.r["physical_3plus1D_locality_or_E8_current_or_gravity_derived"])
        self.assertFalse(self.r["F1_solved"])


if __name__=="__main__":
    unittest.main()
