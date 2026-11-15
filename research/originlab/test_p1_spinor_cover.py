"""Independent verification of the meromorphic 16-spinor Z4 cover."""
import unittest
from collections import Counter
from itertools import combinations
from p1_spinor_cover import controls


class SpinorCoverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=controls()

    def test_five_mode_even_spinor_has_four_of_each_phase(self):
        degrees=[0,1,2,3,0]
        phases=Counter(sum(degrees[i] for i in subset)%4
                       for size in (0,2,4)
                       for subset in combinations(range(5),size))
        self.assertEqual([phases[k] for k in range(4)],[4,4,4,4])
        self.assertEqual(self.result["even_exterior_16_dimension"],16)

    def test_end_spinor_gives_minimal_balanced_256_cover(self):
        self.assertEqual(self.result["endomorphism_256_dimension"],256)
        self.assertEqual(self.result["endomorphism_C4_phase_counts"],[64]*4)
        self.assertEqual(self.result["native_affine_E8_grade_one_phase_counts"],
                         [60,64,60,64])
        self.assertEqual(self.result["minimal_even_grade_complement_dimensions"],
                         [4,0,4,0])
        self.assertEqual(self.result["minimal_total_null_complement"],8)

    def test_character_embedding_not_an_e8_16_dimensional_action(self):
        self.assertTrue(self.result["graded_vector_space_embedding_exists"])
        self.assertFalse(self.result["e8_lie_bracket_embedded_in_gl16"])
        self.assertFalse(self.result["canonical_current_intertwiner_identified"])
        self.assertEqual(self.result["native_current_parity_trace"],-8)
        self.assertEqual(self.result["cover_current_parity_trace"],0)

    def test_no_same_source_claim(self):
        self.assertFalse(self.result["spinor16_equals_16_Majorana_fields_proved"])
        self.assertFalse(self.result["physical_R113_source_proved"])
        self.assertFalse(self.result["F1_solved"])


if __name__=="__main__":
    unittest.main()
