"""Independent exact source Gram controls for original D4 five plus three."""
import unittest
import sympy as sp
from p1_geometric_gram import controls


class GeometricGramSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=controls()

    def test_D4_cross_block_has_one_complex_parameter(self):
        self.assertEqual(self.result["mixed_complex_intertwiner_dimension"],1)
        self.assertEqual(self.result["normalized_mixed_matrix_nonzero_entries"],
                         [[1,0,"z"],[3,2,"-z"]])
        self.assertTrue(self.result["reflection_forbidden_middle_grade_mixing"])

    def test_full_gram_rank_8_vs_rank_6(self):
        witnesses=self.result["regression_witnesses"]
        self.assertEqual(witnesses["0"]["rank"],8)
        self.assertEqual(witnesses["1/2"]["rank"],8)
        self.assertEqual(witnesses["1"]["rank"],6)
        self.assertFalse(witnesses["3/2"]["positive_semidefinite"])
        self.assertTrue(self.result
                         ["rank_5_double_counting_impossible_under_given_D4_and_positive_marginals"])

    def test_exact_determinant_from_two_character_pairs(self):
        z=sp.symbols("z")
        X=sp.zeros(5,3)
        X[1,0]=z
        X[3,2]=-z
        G=sp.eye(5).row_join(X).col_join(X.H.row_join(sp.eye(3)))
        self.assertEqual(sp.factor(G.det()-(1-z*sp.conjugate(z))**2),0)

    def test_source_outstanding(self):
        self.assertFalse(self.result["actual_original_weighted_mixed_gram_z_computed"])
        self.assertFalse(self.result["source_actual_fermion_CAR_assignments_proved"])
        self.assertFalse(self.result["original_R113_Q0_and_time_computed"])
        self.assertFalse(self.result["F1_solved"])


if __name__=="__main__":
    unittest.main()
