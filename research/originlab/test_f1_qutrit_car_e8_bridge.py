"""Exact source-architecture bridge: M3, CAR8, E8 native color, continuum clock."""
import unittest
import sympy as sp
from f1_qutrit_car_e8_bridge import (
    controls, matrix_algebra_type_gate, colour_branch, continuum_time_gate)


class QutritCARE8BridgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report=controls()

    def test_no_unital_M3_embedded_inside_one_irreducible_256_fermion_Fock(self):
        r=self.report["unital_qutrit_CAR_source_embedding"]
        self.assertEqual(r["irreducible_CAR8_Fock_dimension"],256)
        self.assertFalse(r["unital_M3_inside_irreducible_CAR8_Fock_exists"])
        self.assertFalse(r["unital_M3_inside_each_CAR8_even_parity_128_block_exists"])
        self.assertEqual(r["minimum_common_finite_unital_representation_dimension"],768)
        self.assertEqual(r["exact_full_product_algebra"],"M3 tensor M256=M768")

    def test_exhaustive_CAR8_Clifford_relations_and_independent_M3_action(self):
        r=self.report["unital_qutrit_CAR_source_embedding"]
        self.assertEqual(r["exhaustive_CAR_anticommutator_basis_checks"],49152)
        self.assertEqual(r["exhaustive_commuting_qutrit_CAR8_basis_checks"],55296)
        self.assertFalse(r["code_rank_three_compression_preserves_full_CAR_products"])

    def test_native_colour_requires_one_pair_and_two_E8_spinors(self):
        r=self.report["E8_native_color_and_E6_backtransport"]
        self.assertEqual(r["geometric_E8_inverse_weyl_color_positive_roots"],
                         [[0,0,0,0,-2,-2,0,0],
                          [1,1,1,1,1,1,-1,-1],
                          [1,1,1,1,-1,-1,-1,-1]])
        self.assertEqual(r["geometric_E8_inverse_weyl_colour_positive_types"],
                         ["bilinear_D8","spinor","spinor"])
        self.assertEqual(r["all_six_colour_root_types"],
                         {"bilinear_D8":2,"spinor":4})

    def test_colour_E6_commutant_and_six_27_weight_groups(self):
        r=self.report["E8_native_color_and_E6_backtransport"]
        self.assertEqual(r["native_colour_E6_commutant_roots"],72)
        self.assertEqual(r["native_colour_E6_commutant_dimension_including_Cartan"],78)
        self.assertEqual(r["charged_E8_color_weights_each_multiplicity"],27)
        self.assertEqual(r["charged_E8_color_weight_count"],6)
        self.assertEqual(r["exact_E8_branching"],
                         "248=(8,1)+(1,78)+(3,27)+(3bar,27bar)")

    def test_geometric_three_mode_source_colour_not_quarter_neutral(self):
        r=self.report["E8_native_color_and_E6_backtransport"]
        self.assertEqual(r["source_F3_complex_mode_axes"],[3,5,7])
        self.assertEqual(r["source_F3_complex_mode_quarter_phases"],[1,2,3])
        self.assertEqual(r["natural_F3_fundamental_SU3_quarter_current_profile"],
                         [2,2,2,2])
        self.assertEqual(r["native_pinned_color_SU3_quarter_current_profile"],
                         [8,0,0,0])
        self.assertFalse(r["scalar_f2_orientation_character_removes_nontrivial_colour_root_grades"])

    def test_current_level_of_one_fundamental_triplet_vs_eight_adjoint_fermions(self):
        r=self.report["E8_native_color_and_E6_backtransport"]
        self.assertEqual(r["index_of_3_fundamental_plus_5_singlets_SU3_complex_fermions"],1)
        self.assertEqual(r["index_of_8_adjoint_SU3_complex_fermions"],6)

    def test_e8_cannot_be_full_level_one_on_any_finite_host(self):
        r=self.report["unital_qutrit_CAR_source_embedding"]
        self.assertFalse(r["finite_dimension_nonzero_affine_Kac_Moody_level_possible"])
        self.assertTrue(r["infinite_mode_extension_required_for_E8_level_one"])

    def test_fixed_q_depolarization_is_not_strongly_continuous_in_UV(self):
        r=self.report["continuum_microscopic_clock_gate"]
        self.assertFalse(r["strongly_continuous_nontrivial_time_from_fixed_q"])
        self.assertEqual(r["fixed_q_at_any_positive_time_continuum_limit"],"0 (instant depolarization)")
        self.assertEqual(r["microscopic_scaling_for_finite_continuum_rate"],"q(a)=exp(-rate*a)")
        self.assertFalse(r["physical_time_rate_and_absolute_clock_derived"])

    def test_no_F1_promotion_from_algebraic_bridge(self):
        self.assertFalse(self.report["physical_original_Q0_R113_produced"])
        self.assertFalse(self.report["physical_E8_level_one_VOA_produced"])
        self.assertFalse(self.report["physical_chiral_SM_family_count_produced"])
        self.assertFalse(self.report["physical_3plus1D_gravity_produced"])
        self.assertFalse(self.report["F1_solved"])


if __name__=="__main__":
    unittest.main()
