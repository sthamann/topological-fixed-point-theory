"""Original domain inference and trace-to-field inference must pass adversaries."""
import unittest
import sympy as s
from origin_operator_checks import original_perturbation_check, fixed_field_ward_check


class OriginOperatorTests(unittest.TestCase):
    def test_original_order_one_bound_keeps_quarter_and_parity_but_is_not_infinitesimal(self):
        d=original_perturbation_check()
        self.assertTrue(d['anchor_reflection_even'])
        self.assertTrue(d['perturbation_reflection_odd'])
        self.assertEqual(d['relative_bound_ratio_limit_per_unit_amplitude'],'1')
        self.assertEqual(d['at_amplitude_one_slow_eigenvalue_limit'],'0')
        self.assertFalse(d['actual_original_source_excluded'])

    def test_new_ward_data_do_not_determine_fixed_field_four_point_time(self):
        d=fixed_field_ward_check()
        self.assertTrue(d['unitary_exact'])
        self.assertTrue(d['vacuum_and_entire_one_particle_sector_fixed'])
        self.assertTrue(d['same_full_charge_operator'])
        self.assertTrue(d['same_charged_quarter_operator'])
        self.assertTrue(d['thermal_alpha_Ward_operator_residual_zero'])
        self.assertEqual(d['first_neutral_pair_return_at_t_log2'],
                         {'original':'1/2','fixed_fields_conjugated_time':'13/50'})
        self.assertTrue(d['fixed_field_linear_equation_fails_after_conjugation'])
        self.assertFalse(d['geometric_locality_of_the_many_body_conjugation_proved'])


if __name__=='__main__':
    unittest.main()
