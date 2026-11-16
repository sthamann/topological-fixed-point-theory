"""Exact 4-sector Klein to 5-null-point odd code and missing P1 anchor tests."""
import unittest
from f1_klein_phase_space_double import controls


class KleinPhaseSpaceDoubleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=controls()

    def test_five_point_odd_code_recovered_from_four_sector_dual(self):
        self.assertEqual(self.r["nonzero_null_point_count"],5)
        self.assertEqual(self.r["odd_q0_count_including_zero"],6)
        self.assertEqual(self.r["sum_of_five_null_points"],"0000")
        self.assertTrue(self.r["pairwise_symplectic_one"])

    def test_120_symmetries_make_anchor_selection_impossible_without_S4_break(self):
        self.assertEqual(self.r["full_isometry_count"],120)
        self.assertEqual(self.r["one_named_point_stabilizer_order"],24)
        self.assertEqual(self.r["one_named_point_stabilizer"],"S4")
        self.assertFalse(self.r["fifth_P1_point_selected_by_full_symmetry"])

    def test_odd_spin_choice_and_physical_source_still_missing(self):
        self.assertFalse(self.r["extra_spin_polarization_derived"])
        self.assertFalse(self.r["original_R113_source_or_physical_symmetry_breaking_found"])
        self.assertFalse(self.r["F1_solved"])


if __name__=="__main__":
    unittest.main()
