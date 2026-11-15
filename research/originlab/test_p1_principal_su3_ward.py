"""Exact independent tests of the original regional sl2 -> compact SU3 Ward gate."""
import unittest
import sympy as sp
from p1_principal_su3_ward import controls, rep, wedge_rep
from p1_principal_su3 import _basis


class PrincipalSU3WardTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=controls()

    def test_original_regional_four_generators_span_sl2_not_3_commuters(self):
        from source_bridge import regional_transport
        result=regional_transport()
        self.assertEqual(result["span_rank"],3)
        self.assertEqual(result["max_independent_commuting_rank"],1)
        x=sp.symbols("x")
        fields=[sp.Integer(1),x,x*x]
        def bracket(f,g):return sp.expand(f*sp.diff(g,x)-g*sp.diff(f,x))
        self.assertEqual([bracket(fields[i],fields[j])
                          for i,j in ((0,1),(0,2),(1,2))],
                         [1,2*x,x*x])

    def test_unique_SO3_equivariant_two_spin2_to_one_spin1(self):
        self.assertEqual(self.data["principal_SO3_equivariant_wedge2_5_to_3_dimension"],1)
        self.assertEqual(self.data["full_basis_Ward_triples_checked"],512)

    def test_positivity_selects_compact_real_su3(self):
        a,b,l=sp.symbols("a b l",positive=True)
        self.assertEqual(sp.simplify(4*(a*l-b).subs(b,a*l)),0)
        self.assertTrue(self.data["positive_metric_forces_compact_SU3_lambda_positive"])
        self.assertFalse(self.data["lambda_zero_semidirect_compact_ad_invariant_metric"])
        self.assertFalse(self.data["lambda_negative_split_real_SL3_positive_ad_invariant_metric"])

    def test_geometric_cubic_requires_one_compensating_fourth_mark(self):
        self.assertEqual(self.data["geometric_D4_cubic_generator_checks"],112)
        self.assertEqual(self.data["nonzero_original_representation_cubic_sample"],"4")
        self.assertEqual(self.data["unmarked_geometric_cubic_character"],
                         "chi(r)=chi(s)=-1 (forbidden in invariant scalar cap)")
        self.assertTrue(self.data["single_f2_twist_cancels_cubic_character"])
        self.assertFalse(self.data["physical_original_f2_mark_and_cubic_source_bound"])

    def test_cubic_Ward_is_one_missing_original_source_number(self):
        self.assertFalse(self.data["native_original_cubic_source_lambda_computed"])
        self.assertFalse(self.data["native_original_Gram_ratio_b_over_a_computed"])
        self.assertFalse(self.data["native_physical_regional_sl2_action_on_C5_F3_proved"])
        self.assertFalse(self.data["native_R113_and_parent_proved"])
        self.assertFalse(self.data["F1_solved"])


if __name__=="__main__":
    unittest.main()
