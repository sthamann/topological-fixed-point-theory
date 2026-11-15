"""Independent exact regression for the geometric 5+3 -> native E8 Weyl lift."""
import unittest
from collections import Counter
from itertools import combinations
from affine_operators import ROOTS, charge_signature, components
from p1_e8_geometric_lift import (
    controls, GEOM8, GEOM_TORUS, NATIVE_GLUE_TORUS,
    SIMPLE_DOUBLED, WEYL_WORD, character, reflect)


class P1E8GeometricLiftTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report=controls()

    def test_original_eight_complex_and_sixteen_real_modes(self):
        self.assertEqual(len(GEOM8),8)
        self.assertEqual(sum(GEOM8)%4,0)
        phases=Counter([k%4 for k in GEOM8]+[(-k)%4 for k in GEOM8])
        self.assertEqual([phases[g] for g in range(4)],[4]*4)

    def test_SO16_and_halfspinor_native_E8_counts(self):
        report=self.report
        self.assertEqual(report["so16_bilinear_current_phase_counts"],[28,32,28,32])
        self.assertEqual(report["positive_128_spinor_phase_counts"],[32]*4)
        self.assertEqual(report["full_E8_248_phase_counts"],[60,64,60,64])
        self.assertEqual(report["native_E8_248_phase_counts"],[60,64,60,64])

    def test_fixed_lie_root_system_is_D5_plus_A3(self):
        doubled=tuple(2*x for x in GEOM_TORUS)
        zeros=[root for root in ROOTS if character(root,doubled)==0]
        self.assertEqual(len(zeros),52)
        self.assertEqual(components(zeros),["A3","D5"])

    def test_explicit_31_E8_reflections_reproduce_all_240_native_grades(self):
        self.assertEqual(len(WEYL_WORD),31)
        doubled=tuple(2*x for x in GEOM_TORUS)
        for index in WEYL_WORD:
            doubled=reflect(doubled,SIMPLE_DOUBLED[index-1])
        native=tuple(2*x for x in NATIVE_GLUE_TORUS)
        r7=SIMPLE_DOUBLED[6]
        self.assertEqual(doubled,tuple(native[i]-4*r7[i] for i in range(8)))
        self.assertTrue(all(character(root,doubled)==charge_signature(root)[0]
                            for root in ROOTS))
        self.assertTrue(self.report["all_240_native_root_phase_equalities_proved"])

    def test_exact_weyl_word_is_not_plain_SO16_bogoliubov(self):
        r=self.report
        self.assertEqual(r["Weyl_root_type_conversion"],
                         {"D8_to_D8":48,"D8_to_spinor":64,
                          "spinor_to_D8":64,"spinor_to_spinor":64})
        self.assertEqual(r["bilinear_root_to_spinor_example"],
                         {"from":[2,0,2,0,0,0,0,0],
                          "to":[-1,-1,1,-1,1,1,1,-1]})
        self.assertFalse(r["this_Weyl_transformation_can_be_a_plain_O16_Bogoliubov_rotation"])
        self.assertTrue(r["full_E8_spinor_extension_required_for_this_Weyl_word"])

    def test_geometric_CAR_and_full_source_identification_remain_open(self):
        report=self.report
        self.assertFalse(report["five_plus_three_independent_physical_modes_proved"])
        self.assertFalse(report["geometry_to_fermion_CAR_source_identified"])
        self.assertFalse(report["original_Q0_and_R113_produced"])
        self.assertFalse(report["physical_local_3plus1D_parent_proved"])
        self.assertFalse(report["F1_solved"])


if __name__=="__main__":
    unittest.main()
