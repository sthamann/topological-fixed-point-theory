"""Exact source-adjacent E8 root-word modulus falsifiers."""
import unittest
import sympy as sp
from affine_operators import ROOTS, MARKED_P1, MARKED_P2, MARKED_GAMMA, dot, neg, add
from p1_charged_modulus import (root_triple_no_max,
                                paired_determinant_candidate,
                                two_face_fermion_sewing,
                                neutral_chord_exponents,
                                controls)


class ChargedP1ModulusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=controls()

    def test_complete_13440_native_root_triangles(self):
        r=self.result["free_charge_triangle"]
        self.assertEqual((r["e8_roots"],r["partners_per_root"],
                          r["ordered_root_pairs_with_dot_minus_one"]), (240,56,13440))
        self.assertEqual(r["original_marked_word"],
                         ["p2","gamma","minus_p1","vacuum"])

    def test_native_three_root_word_is_neutral_and_symmetric(self):
        self.assertEqual(add(MARKED_P2,MARKED_GAMMA),MARKED_P1)
        roots=[MARKED_P2,MARKED_GAMMA,neg(MARKED_P1)]
        self.assertEqual([dot(a,b) for i,a in enumerate(roots) for b in roots[i+1:]],
                         [-1,-1,-1])
        q=[sp.Matrix(v)/2 for v in roots]
        q.append(sp.zeros(8,1))
        self.assertEqual(neutral_chord_exponents(q),
                         (sp.Integer(-1),sp.Integer(-1),sp.Integer(-1)))

    def test_general_neutral_exponent_identity(self):
        q=[sp.Matrix(v) for v in
           ([1,2],[1,-1],[-2,1],[0,-2])]
        a,b,c=neutral_chord_exponents(q)
        self.assertEqual(a+b,-(q[0]+q[2]).dot(q[0]+q[2]))
        self.assertLessEqual(a+b,0)

    def test_free_root_word_square_is_minimum(self):
        delta=sp.symbols("delta",positive=True,real=True)
        f=1/(4*sp.sin(delta))
        self.assertEqual(sp.diff(f,delta).subs(delta,sp.pi/2),0)
        self.assertEqual(sp.diff(f,delta,2).subs(delta,sp.pi/2),sp.Rational(1,4))
        self.assertEqual(self.result["free_charge_triangle"]["square_stationary_type"],
                         "strict minimum")

    def test_paired_determinants_make_square_max_in_candidate(self):
        t=sp.symbols("t",positive=True)
        w=2*t/(1+t*t)
        self.assertEqual(sp.diff(w,t).subs(t,1),0)
        self.assertEqual(sp.diff(w,t,2).subs(t,1),-1)
        self.assertEqual(self.result["paired_relative_determinant"]
                         ["withheld_ratios_relative_to_square"],
                         {"1/2":"4/5","1/3":"3/5","2/3":"12/13","1":"1"})

    def test_required_net_determinant_exponent_threshold(self):
        x=sp.symbols("x",real=True,positive=True)
        for m, second in ((0,1),(sp.Rational(1,2),0),(1,-1)):
            sample=sp.sin(x)**(2*m-1)
            self.assertEqual(sp.simplify(sp.diff(sample,x,2).subs(x,sp.pi/2)),second)
        self.assertEqual(self.result["paired_relative_determinant"]
                         ["source_net_determinant_exponent_threshold"],
                         "m>1/2 in (det_even*det_odd)^m")

    def test_single_eight_majorana_pfaffian_realizes_paired_determinants(self):
        sewing=self.result["two_face_fermionic_sewing"]
        self.assertEqual(sewing["fermionic_variables_in_this_control"],8)
        self.assertEqual(sewing["same_forward_unit_coupling_pfaffian"],
                         "Pf([[0,K],[-K.T,0]])=det(K)=4*sin(delta)^2")
        self.assertEqual(sewing["charged_root_dressed_at_g_one"],
                         "sin(delta), unique global maximum square")
        self.assertFalse(sewing["actual_original_Dirac_eta_ghost_source_identified"])
        self.assertFalse(sewing["OS_reflection_positivity_proved"])

    def test_forward_and_reflection_reverse_sewing_are_not_equivalent(self):
        sewing=self.result["two_face_fermionic_sewing"]
        self.assertEqual(sewing["reverse_orientation_forward_det"],
                         "(1-g*g)^2 (delta-independent)")
        self.assertTrue(sewing["g_below_one_pole_uncanceled"])
        self.assertEqual(sewing["local_maximum_square_threshold"],
                         "g>sqrt(2)-1 (does not give global maximum if g<1)")
        self.assertFalse(sewing["physical_unit_coupling_and_CAR_orientation_proved"])

    def test_no_r113_or_calderon_source_promoted(self):
        self.assertFalse(self.result["F1_solved"])
        self.assertFalse(self.result["physical_source_identity_proved"])
        candidate=self.result["paired_relative_determinant"]
        self.assertFalse(candidate["original_R113_values_computed"])
        self.assertFalse(candidate["source_Pfaffian_and_ghost_exponent_derived"])
        self.assertFalse(candidate["original_Calderon_deck_operators_identified"])


if __name__=="__main__":
    unittest.main()
