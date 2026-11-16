"""Original five-point odd-quadratic compiler: 30 clocks, six anchored, 
paired positive dynamics and exact quantum mixture obstruction.
"""
import unittest
import sympy as sp
from f1_pregeometric_history_compiler import controls, anchored_clocks


class PregeometricHistoryCompilerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=controls()

    def test_existing_original_q0_code_has_five_points_and_30_clocks(self):
        r=self.result
        self.assertEqual(len(r["original_five_q0_points"]),5)
        self.assertEqual(r["all_unanchored_quarter_clocks"],30)
        self.assertEqual(len(anchored_clocks()[0]),30)

    def test_unanchored_democratic_compiler_erases_all_nontrivial_states(self):
        self.assertEqual(self.result["uniform_unanchored_sum"],
                         "J5/5 (rank one, all nontrivial modes lost)")

    def test_distinguished_anchor_selects_six_equal_symmetry_related_clocks(self):
        self.assertEqual(self.result["one_distinguished_anchor_compatible_clocks"],6)
        self.assertEqual(self.result["one_step_spectrum"],{"1":1,"-1/3":3})
        self.assertEqual(len(anchored_clocks()[1]),6)

    def test_first_OS_nonnegative_return_has_one_vacuum_and_three_modes(self):
        self.assertEqual(self.result["minimum_even_positive_return_steps"],2)
        self.assertEqual(self.result["two_step_spectrum"],{"1":1,"1/9":3})
        self.assertTrue(self.result["paired_two_step_OS_semigroup_compatible"])
        self.assertTrue(self.result["natural_three_dimensional_nontrivial_standard_sector"])

    def test_one_time_generator_and_positive_energy_are_exact(self):
        self.assertEqual(self.result["finite_source_generator"],
                         "H=log(9)*(I4-J4/4)")
        self.assertEqual(self.result["positive_three_mode_energy_probe"],
                         "log(9)/2")

    def test_full_quantum_mixture_loses_modes_and_has_no_bounded_generator(self):
        r=self.result
        self.assertEqual(r["quantum_random_unitary_one_step_eigenspaces"],
                         {"1":2,"-1/3":9,"1/3":3,"0":2})
        self.assertEqual(r["quantum_random_unitary_two_step_eigenspaces"],
                         {"1":2,"1/9":12,"0":2})
        self.assertFalse(r["incoherent_quantum_average_has_bounded_continuous_generator"])

    def test_no_unauthorized_P1_clock_phase_or_three_SM_generation_claims(self):
        r=self.result
        self.assertFalse(r["physical_seam_anchor_selected_by_original_Q0"])
        self.assertFalse(r["three_dimensional_sector_is_three_physical_SM_families"])
        self.assertFalse(r["original_quadratic_refinement_spin_phase_selected"])
        self.assertFalse(r["coherent_average_vs_CP_history_mixture_physically_resolved"])

    def test_not_full_original_Q0_or_3plus1D_F1(self):
        r=self.result
        self.assertFalse(r["physical_R113_Q0_time_and_reference_derived"])
        self.assertFalse(r["actual_pregeometric_3plus1D_gravity_and_E8_derived"])
        self.assertFalse(r["F1_solved"])


if __name__=="__main__":
    unittest.main()
