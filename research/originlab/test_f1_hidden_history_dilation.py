"""Regress exact D4 hidden Markov source, OS positivity, blocking and projection memory."""
import unittest
import sympy as sp
from f1_hidden_history_dilation import controls


class HiddenHistoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=controls()

    def test_four_state_D4(self):
        self.assertEqual(self.r['global_states'],4)
        self.assertEqual(self.r['global_spectrum'],{'1':1,'3/5':1,'7/10':2})

    def test_weighted_two_level_sewing(self):
        self.assertEqual(self.r['exact_weighted_two_level_boundary_pairs'],16)
        self.assertTrue(self.r['global_source_is_reversible_and_finite_block_consistent'])

    def test_full_history_OS(self):
        self.assertTrue(self.r['OS_positive_on_full_history'])

    def test_observed_nonmarkov(self):
        self.assertFalse(self.r['observed_autonomous_markov_transfer_exists_from_one_step'])
        self.assertEqual(self.r['observed_semigroup_defect'],
                         str([[sp.Rational(1,600),-sp.Rational(1,600)],
                              [-sp.Rational(1,1800),sp.Rational(1,1800)]]))

    def test_not_original_no_lindblad_theorem(self):
        self.assertFalse(self.r['native_no_Lindblad_log_proved_by_this_toy'])
        self.assertFalse(self.r['original_TFPT_20_state_or_R113_physical_embedding_proved'])

    def test_F1_not_solved(self):
        self.assertFalse(self.r['F1_solved'])


if __name__=='__main__':
    unittest.main()
