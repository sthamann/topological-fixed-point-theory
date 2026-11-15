"""Exact original Ramond/Spin16/E8 weight-one compatibility regression."""
import unittest
from f1_ramond_e8_grade_bridge import controls


class RamondE8Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=controls()

    def test_exact_E8_root_invariant(self):
        self.assertEqual(self.r['whole_E8_root_system_size'],240)
        self.assertEqual(self.r['simple_root_Gram_determinant'],1)
        self.assertEqual(self.r['root_pairs_checked'],57600)

    def test_all_128_Ramond_spinor_fields_weight_one(self):
        self.assertTrue(self.r['Ramond_even_spinor_weights_all_one'])
        self.assertEqual(self.r['Ramond_ground_spin_field_weight'],'1')
        self.assertEqual(self.r['Ramond_positive_chirality_spin_fields'],128)

    def test_true_E8_weight_one_count_and_quarter(self):
        self.assertEqual(self.r['SO16_NS_bilinear_weight_one_currents'],120)
        self.assertEqual(self.r['E8_weight_one_currents'],248)
        self.assertEqual(self.r['quarter_grading'],[60,64,60,64])

    def test_not_same_as_finite_M3_number_generator(self):
        self.assertFalse(self.r['finite_source_number_operator_automatically_equals_VOA_L0'])

    def test_original_physical_spin_structure_not_selected(self):
        self.assertFalse(self.r['Ramond_spin_structure_GSO_from_original_P1_derived'])
        self.assertFalse(self.r['original_local_chiral_fermion_OPE_and_OS_source_constructed'])
        self.assertFalse(self.r['F1_solved'])


if __name__=='__main__':
    unittest.main()
