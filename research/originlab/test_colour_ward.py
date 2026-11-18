import unittest

from colour_ward import controls, native_carrier_check


class ColourWardTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = controls()

    def test_neutral_spinor_connects_the_charged_ward_relations(self):
        self.assertEqual(self.report["grade_one_dimensions"],[86,81,81])
        self.assertEqual(self.report["D8_pair_roots_by_colour"],[46,33,33])
        self.assertEqual(self.report["spinor_roots_by_colour"],[32,48,48])
        self.assertEqual(self.report["neutral_root_counts_along_extension"],[46,78,240])
        for raw,complete in zip(self.report["D8_neutral_ward"],self.report["charged_ward"]):
            self.assertEqual(raw["neutral_root_orbit_sizes"],[3,30,48])
            self.assertEqual(raw["cyclic_ward_component_sizes"],[90,720])
            self.assertEqual(complete["relative_complex_parameters"],1)
            self.assertEqual(complete["relation_rank"],809)
            self.assertEqual(complete["nonzero_ward_rows"],21870)
            self.assertEqual(complete["oriented_cubic_entries"],1620)
            self.assertEqual(complete["unordered_cubic_triangles"],270)
        self.assertEqual(self.report["cubic_scale_modulus_squared"],1)
        with self.assertRaisesRegex(ValueError,"fixed physical colour weights"):
            native_carrier_check([(1,0,-1,0),(1,0,1,-1),(1,0,0,1)])

    def test_neutral_current_and_charged_primary_have_different_positive_modes(self):
        self.assertEqual(self.report["positive_modes"]["neutral_current_positive_mode_vacuum_coefficient"],1)
        self.assertTrue(self.report["positive_modes"]["charged_grade_one_annihilated_by_neutral_positive_modes"])
        strings = self.report["neutral_marked_lowering"]
        self.assertEqual([s["lowering_rank"] for s in strings],[22,18,18])
        self.assertEqual([s["lowering_square_rank"] for s in strings],[1,0,0])
        self.assertEqual([s["casimir_trace"] for s in strings],["36","27","27"])
        self.assertEqual((strings[2]["spinor_to_D8"],strings[2]["D8_to_spinor"]),(15,3))

    def test_cyclic_reordering_and_charge_sensitive_reverse_control(self):
        word = self.report["marked_ordering"]
        self.assertEqual(word["tested_mode_pairs"],16)
        self.assertEqual(word["early_to_terminal_Y"],["1/3","1/3","5/6"])
        self.assertEqual(word["reverse_two_current_word_norm_squared"],0)
        self.assertEqual(word["forward_two_current_word_norm_squared_at_mode_m"],"m")
        self.assertIsNone(self.report["original_values"])
