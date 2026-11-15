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

    def test_exact_D4_equivariant_gram_whitening_in_positive_phase(self):
        report=self.result
        self.assertTrue(report["source_D4_equivariant_whitening_exists_for_every_strictly_positive_z"])
        self.assertTrue(report["whitening_preserves_geometric_C4_D4_operators"])
        self.assertTrue(report["source_parameter_z_does_not_affect_abstract_D4_CAR_character_within_positive_rank_8"])
        self.assertTrue(report["source_parameter_z_still_needed_for_weighted_Q0_and_time"])
        self.assertEqual(report["full_8_CAR_rank_condition"],"|z|<1")

    def test_exact_Gamma_source_Ward_isometry_suffices_for_rank_eight(self):
        r=self.result
        self.assertTrue(r["form_parity_ward_orthogonality_forces_z_zero"])
        self.assertEqual(r["source_graded_eight_mode_positive_rank_when_Gamma_isometric"],8)
        self.assertTrue(r["candidate_D_physical_G_selfadjoint_only_for_real_z"])
        self.assertFalse(r["original_spectral_triple_Gamma_identified_with_form_degree"])
        z=sp.symbols("z")
        M=sp.Matrix(5,3,lambda i,j:z*(int(i==1 and j==0)-
                                      int(i==3 and j==2)))
        G=sp.eye(5).row_join(M).col_join(M.H.row_join(sp.eye(3)))
        gam=sp.diag(*([1]*5+[-1]*3))
        self.assertEqual((gam*G*gam-G).subs(z,0),sp.zeros(8))
        self.assertNotEqual((gam*G*gam-G).subs(z,sp.Rational(1,2)),sp.zeros(8))

    def test_source_outstanding(self):
        self.assertFalse(self.result["actual_original_weighted_mixed_gram_z_computed"])
        self.assertFalse(self.result["source_actual_fermion_CAR_assignments_proved"])
        self.assertFalse(self.result["original_R113_Q0_and_time_computed"])
        self.assertFalse(self.result["F1_solved"])


if __name__=="__main__":
    unittest.main()
