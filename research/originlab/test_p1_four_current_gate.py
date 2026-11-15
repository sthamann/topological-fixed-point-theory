"""Independent exact chiral four-current and same-cap parity controls."""
import unittest
import sympy as sp
from affine_operators import MARKED_P2, MARKED_GAMMA, ROOT_SET, neg
from p1_four_current_gate import free_vertex, controls


class FourCurrentF1Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=controls()

    def test_all_four_e8_current_insertions_nonzero_and_neutral(self):
        a,b=MARKED_P2,MARKED_GAMMA
        for word in ((a,b,neg(a),neg(b)),(a,neg(b),neg(a),b)):
            self.assertTrue(all(v in ROOT_SET for v in word))
            self.assertTrue(all(sum(v[k] for v in word)==0 for k in range(8)))
        self.assertEqual(self.result["native_E8_nonzero_currents_in_each_word"],4)

    def test_root_fourpoint_ratio_is_crossratio_squared(self):
        a,b=MARKED_P2,MARKED_GAMMA
        u=sp.symbols("u",nonzero=True)
        z=(1,u,-1,-u)
        a4=free_vertex(z,(a,b,neg(a),neg(b)))
        b4=free_vertex(z,(a,neg(b),neg(a),b))
        x=sp.factor((z[0]-z[1])*(z[2]-z[3])/
                    ((z[0]-z[2])*(z[1]-z[3])))
        self.assertEqual(sp.factor(a4/b4-((1-x)/x)**2),0)

    def test_original_v512_neutral_odd_split_predicts_charged_ratio(self):
        t=sp.symbols("t",positive=True)
        self.assertEqual(sp.simplify(sp.exp(2*sp.pi*(2/sp.pi)*sp.log(1/t))-t**(-4)),0)
        self.assertEqual(self.result["non_square_original_v512_t_half_channel_ratio"],"16")
        self.assertEqual(self.result["square_t_one_channel_ratio"],"1")

    def test_crossing_does_not_mean_same_frame_amplitudes_equal(self):
        self.assertFalse(self.result["crossing_symmetry_alone_selects_square"])
        self.assertFalse(self.result["same_frame_channel_equality_from_original_source_proved"])
        self.assertFalse(self.result["three_point_carries_fourth_mark_invariant"])

    def test_even_parity_trace_constant_is_not_determinant_product(self):
        self.assertEqual(self.result["even_fermion_projection_trace"],"2 (all moduli)")
        self.assertEqual(self.result["ordinary_times_supertrace"],
                         "16*t*t/(1+t*t)^2")
        self.assertFalse(self.result["original_Padm_equals_literal_even_Fock_projector_proved"])

    def test_no_original_R113_or_F1_promoted(self):
        self.assertFalse(self.result["original_R113_UY_sources_computed"])
        self.assertFalse(self.result["physical_3plus1_parent_derived"])
        self.assertFalse(self.result["F1_solved"])


if __name__=="__main__":
    unittest.main()
