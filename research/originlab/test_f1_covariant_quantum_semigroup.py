"""Regression of exact 648-equation SU3 conjugation-covariant quantum semigroup."""
import unittest
from f1_covariant_quantum_semigroup import controls


class CovariantSemigroupTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=controls()

    def test_81_dim_linear_commutant_reduced_to_two(self):
        self.assertEqual(self.r['commutant_linear_rank'],79)
        self.assertEqual(self.r['commutant_dimension'],2)

    def test_su3_adjoint_and_scalar(self):
        self.assertEqual(self.r['scalar_subrepresentation_dim'],1)
        self.assertEqual(self.r['traceless_adjoint_subrepresentation_dim'],8)

    def test_CP_exact_interval(self):
        self.assertEqual(self.r['complete_positivity_exact_interval'],'-1/8<=q<=1')

    def test_positive_continuous_semigroup_unique_up_to_clock(self):
        self.assertEqual(self.r['continuous_CP_semigroup_form'],
                         'q(t)=exp(-kappa*t), kappa>=0')

    def test_source_clock_scale_not_selected(self):
        self.assertFalse(self.r['source_covariance_and_markov_requirements_choose_kappa'])

    def test_physical_e8_and_spacetime_not_claimed(self):
        self.assertFalse(self.r['physical_P1_mark_geometry_or_E8_reconstructed'])
        self.assertFalse(self.r['F1_solved'])


if __name__=='__main__':
    unittest.main()
