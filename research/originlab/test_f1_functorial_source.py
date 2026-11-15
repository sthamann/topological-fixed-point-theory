"""Exact pre-spatial, fully positive, reflected and sewn finite source controls."""
import unittest
import sympy as sp
from f1_functorial_source import controls, transition, quantum_transfer, rank_and_gram


class TestFunctorialSource(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=controls()

    def test_CP_positive(self):
        self.assertEqual(self.r['Kraus_operators'],10)
        self.assertEqual(self.r['OS_Gram_rank'],3)

    def test_exact_two_level_block(self):
        self.assertEqual(self.r['fine_to_coarse_exact_endpoint_pairs_verified'],9)
        self.assertEqual(transition(3,sp.Rational(1,2))**2,
                         transition(3,sp.Rational(1,4)))

    def test_OS_Gram(self):
        G=rank_and_gram(sp.Rational(1,2))
        self.assertTrue(G.is_positive_definite)

    def test_all_observables_quantum_sewing(self):
        q=sp.Rational(1,3)
        a=sp.Matrix([[3,2,5],[0,1,4],[3,0,-2]])
        self.assertEqual(quantum_transfer(quantum_transfer(a,q),q),
                         quantum_transfer(a,q*q))

    def test_native_readouts_not_faked(self):
        self.assertFalse(self.r['original_R113_computed'])
        self.assertFalse(self.r['original_TFPT_Q0_identified'])

    def test_observations_unidentifiable_without_rate(self):
        self.assertTrue(self.r['same_lower_three_source_marks_different_connected_UY'])
        self.assertEqual(self.r['control_q_half']['CUY'],'1/9')
        self.assertEqual(self.r['control_q_quarter']['CUY'],'1/18')

    def test_exact_fixed_scale_only_trivial_solutions(self):
        x=sp.Symbol('q')
        self.assertEqual(sp.solve(x*x-x,x),[0,1])

    def test_permutation_symmetry_and_full_classical_reflection(self):
        self.assertEqual(self.r['state_permutation_automorphism_group'],
                         'S3 (all 6 permutations)')
        self.assertEqual(self.r['full_reflected_classical_history_one_step_eigenvalues'],
                         {'1/3':1,'1/12':2})

    def test_zero_dimension_in_uv_of_finite_M3(self):
        self.assertEqual(self.r['ultraviolet_spectral_dimension_limit'],0)
        self.assertFalse(self.r['four_dimensional_continuum_derived'])

    def test_twenty_state_analogue(self):
        self.assertEqual(self.r['twenty_state_analog_stationary'],'I20/20')

    def test_no_spacetime_or_TOE(self):
        self.assertFalse(self.r['emergent_spatial_geometry_and_spin2_proved'])
        self.assertFalse(self.r['F1_solved'])


if __name__=='__main__':
    unittest.main()
