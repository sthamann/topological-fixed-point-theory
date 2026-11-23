"""Independent P1 origin physics constraint tests: do not promote to F1."""
import unittest
from f1_p1_physical_necessity_routes import controls


class P1PhysicalNecessityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=controls()

    def test_proper_finite_modular_inclusion_is_impossible(self):
        x=self.r["finite_modular_translation"]
        self.assertFalse(x["finite_common_cyclic_separating_algebras_can_be_properly_nested"])
        self.assertFalse(x["finite_positive_nonzero_translation_satisfies_exact_modular_dilation"])
        self.assertTrue(x["standard_modular_inclusion_requires_nontrivial_infinite_source"])

    def test_gauge_matching_leaves_coherent_history_phases_undetermined(self):
        x=self.r["gauge_interface"]
        self.assertEqual(x["edge_matching_projection_rank"],4)
        self.assertTrue(x["both_history_states_admissible"])
        self.assertTrue(x["both_one_side_reduced_states_identical"])
        self.assertEqual(x["cross_mark_values"],["+1","-1"])
        self.assertFalse(x["gauss_matching_fixes_original_coherent_amplitudes"])

    def test_stationary_four_clock_gives_relational_not_unique_time(self):
        x=self.r["relational_time"]
        self.assertTrue(x["invariant_global_state"])
        self.assertEqual(x["constraint_fixed_subspace_dimension"],4)
        self.assertEqual(x["clock_reduced_state"],"I4/4")
        self.assertTrue(x["conditionally_nontrivial_system_evolution"])
        self.assertFalse(x["unique_p1_or_arrow_from_stationarity"])

    def test_independent_matching_and_cyclic_constraints_have_unique_joint_state(self):
        x=self.r["constraint_intersection"]
        self.assertEqual(x["matching_projector_rank"],4)
        self.assertEqual(x["global_cycle_constraint_projector_rank"],4)
        self.assertEqual(x["commuting_constraint_intersection_rank"],1)
        self.assertTrue(x["relative_four_amplitudes_unique_under_assumed_constraints"])
        self.assertFalse(x["cyclic_orientation_selected"])
        self.assertFalse(x["original_six_D4_history_orbit_weight_selected"])

    def test_three_region_coherence_cannot_be_inferred_from_pair_marginals(self):
        x=self.r["quantum_three_region_memory"]
        self.assertTrue(x["same_AB_BC_AC_pair_marginals"])
        self.assertEqual(x["three_region_matching_rank"],4)
        self.assertEqual(x["three_region_cyclic_constraint_rank"],16)
        self.assertEqual(x["unique_three_region_joint_constraint_rank"],1)
        self.assertEqual(x["conditional_mutual_information_GHZ_bits"],2)
        self.assertEqual(x["conditional_mutual_information_dephased_bits"],0)
        self.assertFalse(x["global_coherence_recoverable_from_pair_marginals_only"])
        self.assertFalse(x["exact_quantum_Markov_Petz_recovery_for_GHZ"])
        self.assertTrue(x["exact_quantum_Markov_Petz_recovery_for_dephased"])

    def test_topology_leaves_square_vs_hexagonal_modulus_open(self):
        x=self.r["elliptic_modulus"]
        self.assertEqual(x["four_simple_branch_points_genus"],1)
        self.assertEqual(x["square_harmonic_cross_ratio_j"],1728)
        self.assertEqual(x["hexagonal_equi_anharmonic_cross_ratio_j"],0)
        self.assertFalse(x["genus_one_and_modular_symmetry_uniquely_force_square"])

    def test_bounce_and_Zeno_are_not_pregeometric_F1_primitives(self):
        x=self.r["bounce_and_zeno"]
        self.assertFalse(x["flat_GR_bounce_under_NEC"])
        self.assertFalse(x["effective_bounce_is_pregeometric_P1_derivation"])
        self.assertFalse(x["repeated_rank_one_projection_generates_nontrivial_clock"])

    def test_native_original_source_state_is_reported(self):
        x=self.r["original_source"]
        self.assertEqual(x["native_R113_source_bundle"] is not None,x["native_R113_source_review_required"])
        self.assertFalse(x["diagnostic_circle_is_physical_OS_time"])
        self.assertFalse(self.r["F1_solved"])


if __name__=="__main__":
    unittest.main()
