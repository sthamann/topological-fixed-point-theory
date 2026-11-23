"""Direct regression on the *actual* TFPT R113 cap and Higgs kinetic sources."""
import unittest
import sympy as sp
from f1_single_original_source_gate import controls


class F1SingleOriginalSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=controls()

    def test_original_R113_domain_is_not_physical_Q0_or_OS_time(self):
        r=self.r["original_native_status"]
        # Reflects the current domain file; a new bundle is flagged for review, not rejected.
        self.assertEqual(r["full_original_source_packet_present"],r["source_packet_review_required"])
        self.assertFalse(r["original_spatial_coordinate_is_OS_time"])
        self.assertEqual(r["original_pre_projector_zero_fermion_cap"],"1")
        self.assertFalse(r["conditional_center_trace_is_original_empty_cap"])

    def test_even_Fock_projector_can_make_original_empty_source_R_dependent(self):
        r=self.r["same_source_empty_cap_and_selector"]
        self.assertEqual(r["empty_original_fermion_cap"],"1 for any R")
        self.assertEqual(r["allowed_even_parity_projected_cap"],"(1+r*s)/4")
        self.assertTrue(r["parity_preserved_but_number_mixed"])
        self.assertFalse(r["physical_original_P_adm_number_action_computed"])

    def test_full_trace_phase_scan_not_equivalent_to_all_ordered_R113_marks(self):
        r=self.r["same_source_empty_cap_and_selector"]
        self.assertEqual(r["unmarked_empty_cap_phase_scan"],"1")
        self.assertEqual(r["same_Fock_source_ordered_even_pair"],"r*s*z")
        self.assertFalse(r["all_original_R113_marks_follow_from_unmarked_phasescan_alone"])

    def test_NINE_native_full_field_Higgs_gradient_terms_replay_original(self):
        r=self.r["native_full_Higgs_kinetic_direction"]
        self.assertEqual(r["original_U7_readout"],"13*a/100")
        self.assertEqual(r["nine_original_nonzero_gradient_directions"],9)
        self.assertEqual(r["source_gradient_coefficients_full"],
                         {"-1":"5117/5000","1":"5117/5000"})
        self.assertEqual(r["full_source_over_one_mode"],"40936/169")
        self.assertEqual(r["full_source_kinetic_local_coefficient"],
                         "(10234/625)*pi/R")

    def test_kinetic_form_alone_has_no_actual_physical_time_generator(self):
        r=self.r["native_full_Higgs_kinetic_direction"]
        self.assertFalse(r["physical_source_P_adm_reduced_energy_form_derived"])
        self.assertFalse(r["physical_OS_time_generator_derived_from_local_Z_H"])

    def test_common_energy_moments_expose_missing_source_state_directions(self):
        r=self.r["hidden_source_generator_moments"]
        self.assertEqual(r["finite_exact_comparison_hidden_rank"],1)
        self.assertEqual(r["native_spatial_full96_first_hidden_rank"],192)
        self.assertFalse(r["spatial_full96_return_proves_physical_OS_time"])
        self.assertFalse(r["original_common_Higgs_family_clock_stress_M0_M1_M2_available"])

    def test_tracial_M3_source_modular_flow_is_trivial(self):
        r=self.r["state_and_modular_clock"]
        self.assertFalse(r["tracial_M3_modular_flow_nontrivial"])
        self.assertEqual(r["nontrivial_modular_E01_frequency"],"log(3/2)")
        self.assertFalse(r["physical_source_nontracial_state_derived"])

    def test_original_full_theory_not_derived_by_a_comparison_model(self):
        self.assertFalse(self.r["original_full_R113_four_numbers_computed"])
        self.assertFalse(self.r["original_common_positive_generator_constructed"])
        self.assertFalse(self.r["original_physical_3plus1D_E8_SM_spin2_derived"])
        self.assertFalse(self.r["F1_solved"])


if __name__=="__main__":
    unittest.main()
