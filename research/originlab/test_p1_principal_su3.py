"""Exact original geometric 5+3 SU3 Lie closure and chirality-twist controls."""
import unittest
import sympy as sp

from p1_principal_su3 import controls, _basis, _E


class P1PrincipalSU3Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=controls()

    def test_compact_su3_eight_generators_and_all_lie_brackets(self):
        self.assertEqual(self.data["compact_eight_generator_lie_algebra"],
                         "su3 (real) / sl3(C) (complex)")
        self.assertEqual(self.data["full_commutator_pairs"],28)
        self.assertEqual(self.data["jacobi_generator_triples_checked"],56)
        self.assertEqual(self.data["derived_lie_rank"],8)
        self.assertEqual(self.data["positive_compact_metric_eigenvalues"],{"2":7,"6":1})

    def test_principal_SO3_decomposes_as_spin_2_plus_spin_1(self):
        s,a=_basis()
        self.assertEqual((len(s),len(a)),(5,3))
        self.assertEqual(self.data["principal_so3_spin1_dimension"],3)
        self.assertEqual(self.data["principal_so3_spin2_dimension"],5)
        # In su3, commutator of imaginary symmetric traceless generators is
        # a real antisymmetric one, absent in the abelian comparison class.
        x,y=sp.I*s[0],sp.I*s[3]
        self.assertNotEqual(x*y-y*x,sp.zeros(3))
        self.assertTrue((x*y-y*x).T==-(x*y-y*x))

    def test_exact_D4_sign_character_and_internal_geometric_f2_line(self):
        self.assertEqual(self.data["compensator_character"],"chi(r)=chi(s)=-1")
        self.assertTrue(self.data["twisted_geometric_D4_can_act_as_principal_SU3_automorphism"])
        self.assertFalse(self.data["untwisted_geometric_D4_is_an_SU3_bracket_automorphism"])
        self.assertEqual(self.data["geometric_five_intertwiner_determinant"],"32768*I")
        self.assertEqual(self.data["geometric_three_intertwiner_determinant"],"-512*I")
        z=sp.symbols("z")
        f2=z*z/(z**4-1)
        self.assertEqual(sp.simplify(f2.subs(z,sp.I*z)+f2),0)
        self.assertEqual(sp.simplify(f2.subs(z,1/z)+f2),0)
        self.assertEqual(sp.limit(f2/z**2,z,0),-1)
        self.assertEqual(sp.limit(f2*z**2,z,sp.oo),1)

    def test_rank_two_principal_selection_requires_explicit_source_premise(self):
        self.assertEqual(self.data["principal_rank_two_simple_n_mark_matches"],[("A2",4)])
        self.assertFalse(self.data["source_bracket_and_cubic_OPE_identified"])

    def test_SU3_embedding_in_E8_is_not_selected_by_dimensions(self):
        self.assertEqual(self.data["all_E8_A2_root_subalgebras"],1120)
        self.assertEqual(self.data["quarter_fixed_E8_A2_root_subalgebras"],44)
        self.assertFalse(self.data["physical_colour_SU3_embedding_selected"])

    def test_native_c4_profiles_classify_all_1120_regular_su3_embeddings(self):
        r=self.data
        self.assertEqual(r["E8_A2_quarter_grade_profiles"],
                         {"(2, 2, 2, 2)":480,"(4, 0, 4, 0)":180,
                          "(4, 2, 0, 2)":416,"(8, 0, 0, 0)":44})
        self.assertEqual(r["quarter_balanced_E8_A2_subalgebras"],480)

    def test_geometric_balanced_su3_cannot_commute_with_full_native_SM(self):
        r=self.data
        self.assertEqual(r["quarter_balanced_A2_commuting_with_native_colour"],48)
        self.assertEqual(r["A2_commuting_with_native_full_colour_weak_hypercharge"],10)
        self.assertEqual(r["quarter_balanced_A2_commuting_with_native_full_SM"],0)
        self.assertEqual(r["native_full_SM_A2_centralizer_quarter_profiles"],
                         {"(4, 2, 0, 2)":6,"(8, 0, 0, 0)":4})
        self.assertEqual(r["full_SM_no_go_scope"],
                         "regular E8 root A2 embedding in pinned physical D8 charge chart")

    def test_same_D4_character_does_not_choose_between_CAR_and_su3(self):
        self.assertFalse(self.data["8_complex_free_fermion_interpretation_excluded"])
        self.assertFalse(self.data["original_weighted_Q0_and_R113_produced"])
        self.assertFalse(self.data["F1_solved"])


if __name__=="__main__":
    unittest.main()
