"""49,152 CAR and exact M3->Spin16/E8 phase / affine-grade no-go tests."""
import unittest
from f1_single_source_clifford_bridge import controls


class SourceCliffordBridgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=controls()

    def test_positive_exact_oneparticle_metric(self):
        self.assertEqual(self.r['source_one_particle_hilbert'],
                         'traceless M3(C), Hilbert-Schmidt metric, dimension 8')
        self.assertEqual(self.r['positive_metric_diagonal'],[2]*7+[6])

    def test_direct_49152_CAR_verifications(self):
        self.assertEqual(self.r['all_CAR_pairs_on_256_Fock_masks_verified'],49152)
        self.assertEqual(self.r['real_Majorana_generators'],16)

    def test_E8_abstract_representation_counts(self):
        self.assertEqual(self.r['positive_even_spinor_dimension'],128)
        self.assertEqual(self.r['SO16_bilinear_dimension'],120)
        self.assertEqual(self.r['E8_abstract_quarter_phase_counts'],[60,64,60,64])

    def test_zeroFock_cap_is_not_trace_or_supertrace(self):
        self.assertEqual(self.r['original_zero_Fock_cap_of_lift'],'1 (no determinant)')
        self.assertEqual(self.r['lifted_transfer_trace_at_q_half'],'6561/256')
        self.assertEqual(self.r['lifted_transfer_supertrace_at_q_half'],'1/256')

    def test_affine_E8_weight_one_does_not_follow_from_number_operator(self):
        self.assertEqual(self.r['different_even_spinor_source_energies'],[0,2,4,6,8])
        self.assertFalse(self.r['nontrivial_common_weight_one_possible_by_affine_rescaling'])

    def test_original_OPE_and_F1_still_open(self):
        self.assertFalse(self.r['original_physical_CAR_or_E8_current_OPE_derived'])
        self.assertFalse(self.r['F1_solved'])


if __name__=='__main__':
    unittest.main()
