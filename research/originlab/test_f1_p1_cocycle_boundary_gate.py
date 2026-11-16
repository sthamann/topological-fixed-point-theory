"""Tests for independent local-gauge P1 and twisted-bulk boundary source."""
import unittest
from f1_p1_cocycle_boundary_gate import controls


class P1CocycleBoundaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=controls()

    def test_four_edge_Gauss_plus_flux_remains_fourfold_not_unique(self):
        p=self.r["local_gauge_seam"]
        self.assertTrue(p["four_edge_star_Gauss_law_commutes_with_even_flux_overlap"])
        self.assertTrue(p["odd_overlap_anticommutes"])
        self.assertEqual(p["two_commuting_checks_joint_rank"],4)
        self.assertFalse(p["unique_P1_from_local_Gauss_plus_one_flux"])

    def test_twisted_Z2_cocycle_has_nontrivial_boundary_Frobenius_obstruction(self):
        p=self.r["Z2_nontrivial_associator"]
        self.assertEqual(p["Z2_cocycle_identity_checked"],16)
        self.assertEqual(p["omega_111"],-1)
        self.assertEqual(p["normalized_two_cochain_coboundary_111"],"+1")
        self.assertFalse(p["full_subgroup_K_Z2_cancels_omega_with_boundary_two_cochain"])
        self.assertTrue(p["trivial_subgroup_K_identity_restriction_unobstructed"])

    def test_semionic_quarter_spin_does_not_select_E8_chiral_absolute_source(self):
        p=self.r["double_semion_and_E8"]
        self.assertEqual(p["double_semion_topological_spins"],["1","i","-i","1"])
        self.assertEqual(p["double_semion_chiral_c"],0)
        self.assertEqual(p["double_semion_stacked_E8_chiral_c"],8)
        self.assertEqual(p["E8_determinant"],1)
        self.assertEqual(p["double_semion_stacked_E8_anyon_count"],4)
        self.assertFalse(p["anyon_spins_distinguish_E8_stack"])
        self.assertTrue(p["requires_original_framing_eta_or_absolute_reference"])

    def test_physical_R113_and_F1_remain_open(self):
        self.assertFalse(self.r["physical_P1_from_pregeometric_original_Q0"])
        self.assertFalse(self.r["original_R113_outer_Padm_and_all_fields_evaluated"])
        self.assertFalse(self.r["F1_solved"])


if __name__=="__main__":
    unittest.main()
