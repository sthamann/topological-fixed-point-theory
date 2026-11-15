"""Original-source countermodel and determinant regression checks."""
import unittest
import sympy as sp
from p1_source_falsifiers import controls, four_mark_family, gaussian_arc_determinant, full_gaussian_gluing, relative_polarization_clock, calderon_hardy_obstruction


class OriginalP1SourceGateTests(unittest.TestCase):
    def test_universal_four_mark_three_plus_two_is_not_square_selector(self):
        result=four_mark_family()
        self.assertTrue(result["three_plus_two_equivariant_for_all_noncolliding_moduli"])
        self.assertTrue(result["unit_dlog_winding_for_all_moduli"])
        self.assertFalse(result["square_modulus_selected"])

    def test_exact_non_square_witness_and_quarter_failure(self):
        z=sp.symbols("z")
        u=(3+4*sp.I)/5
        p=sp.expand((z*z-1)*(z*z-u*u))
        self.assertEqual(sp.simplify(p.subs(z,-z)-p),0)
        self.assertEqual(sp.simplify(p.subs(z,sp.I*z)-p-2*z*z*(1+u*u)),0)
        self.assertNotEqual(sp.simplify(1+u*u),0)

    def test_gaussian_determinant_selects_square_conditionally(self):
        result=gaussian_arc_determinant()
        self.assertEqual(result["unique_maximum_in_0_delta_pi"],"delta=pi/2 (square)")
        self.assertFalse(result["physical_original_boson_action_identified"])
        self.assertFalse(result["fermion_ghost_eta_factors_included"])

    def test_original_arc_matrix_tree_exact_formula(self):
        a,b=sp.symbols("a b",positive=True)
        lap=sp.Matrix([[a+b,-a,0,-b],[-a,a+b,-b,0],
                       [0,-b,a+b,-a],[-b,0,-a,a+b]])
        self.assertEqual(sp.factor((lap+sp.ones(4,4)/4).det()),8*a*b*(a+b))
        self.assertEqual(sp.factor(lap[:3,:3].det()),2*a*b*(a+b))

    def test_complete_gaussian_gluing_cancels_vertex_square_bias(self):
        result=full_gaussian_gluing()
        self.assertTrue(result["physical_free_boson_partition_independent_of_subdivision"])
        self.assertTrue(result["vertex_only_square_saddle_canceled_by_interior_modes"])
        self.assertEqual(result["combined_determinant"],
                         "64*perimeter (128*pi when perimeter=2*pi)")
        self.assertFalse(result["physical_F1_modulus_selected"])

    def test_cross_ratio_as_relative_polarization_trace(self):
        result=relative_polarization_clock()
        self.assertEqual(result["cross_ratio_from_trace"],
                         "lambda=(trace(U)-2)/(trace(U)+2)=-t^2")
        self.assertEqual(result["square_iff"],
                         "trace(U)=0 iff t=1 (t>0); then U^2=-I")
        self.assertFalse(result["microlocal_source_intertwiner_derived"])
        self.assertFalse(result["physical_R113_values_produced"])

    def test_two_fourier_modes_refute_local_calderon_deck_equality(self):
        result=calderon_hardy_obstruction()
        self.assertTrue(result["no_frequency_independent_boundary_deck_matrix"])
        self.assertFalse(result["full_TFPT_original_boundary_conditional_relation_proved"])

    def test_preserves_f1_physical_gate(self):
        result=controls()
        self.assertFalse(result["F1_solved"])
        self.assertFalse(result["physical_source_identity_proved"])
        self.assertFalse(result["gaussian_arc_determinant"]["original_R113_computed"])


if __name__=="__main__":
    unittest.main()
