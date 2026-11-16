"""Exact Klein exchangeable infinite-source and spin refinement controls."""
import unittest
from f1_klein_twisted_source import controls


class P1KleinTwistedSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=controls()

    def test_hermitian_isotropy_rule_selects_Klein_not_alternating_pairing(self):
        x=self.r["bicharacter"]
        self.assertEqual(x["G"],"Z2 x Z2")
        self.assertEqual(x["klein_self_isotropic_subgroup"],["00","11"])
        self.assertTrue(x["klein_u_on_isotropic_subgroup_trivial"])
        self.assertFalse(x["alternating_symplectic_u_on_its_isotropic_subgroup_trivial"])
        self.assertFalse(x["klein_directly_isometric_to_original_alternating_TFPT_pairing"])
        self.assertFalse(x["cyclic_Z4_primitive_quarter_bicharacter_Hermitian"])

    def test_exact_four_quadratic_refinements_and_two_semionic_orientations(self):
        x=self.r["quadratic_spin"]
        self.assertEqual(x["q_twisted_quadratic_refinement_count"],4)
        self.assertEqual(x["double_semion_c0_mod8_refinements"],2)
        self.assertEqual(x["other_chiral_c_mod8_classes"],[2,6])
        self.assertTrue(x["reflection_exchanges_two_double_semion_orientations"])
        self.assertFalse(x["native_P1_orienting_spin_character_derived"])
        self.assertFalse(x["absolute_E8_stacking_number_fixed"])

    def test_2026_external_theorem_is_not_a_physical_native_F1_proof(self):
        self.assertFalse(self.r["infinite_twisted_Cstar_algebra_embedded_in_original_TFPT"])
        self.assertFalse(self.r["actual_original_R113_Padm_Q0_time_and_framing_selected"])
        self.assertFalse(self.r["physical_3plus1D_locality_or_spin2_derived"])
        self.assertFalse(self.r["F1_solved"])


if __name__=="__main__":
    unittest.main()
