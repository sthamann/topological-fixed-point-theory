"""Exact E8 lattice bulk-edge anomaly and SM gauge-I6 non-equivalence tests."""
import unittest
from fractions import Fraction
import sympy as sp
from f1_e8_gravitational_inflow import controls, standard_model_left_weyl_anomalies


class E8InflowSourceGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=controls()

    def test_e8_native_root_cartan_is_positive_even_unimodular(self):
        r=self.r
        self.assertEqual(r["bulk_K_determinant"],1)
        self.assertEqual(r["bulk_K_signature"],8)
        self.assertEqual(r["bulk_K_principal_minors"],[2,4,6,5,4,3,2,1])
        self.assertTrue(r["bulk_K_even_integral_unimodular"])
        self.assertEqual(r["bulk_anyonic_ground_state_deg_on_torus"],1)

    def test_single_e8_chiral_character_carries_nontrivial_modular_T(self):
        phase=sp.exp(-2*sp.pi*sp.I/3)
        self.assertNotEqual(sp.simplify(phase),1)
        self.assertEqual(sp.simplify(phase**3),1)
        self.assertEqual(self.r["E8_1_modular_T_phase"],"exp(-2*pi*i/3)")
        self.assertEqual(self.r["local_two_dimensional_gravitational_I4_coefficient_magnitude"],"p1/3")

    def test_two_oppositely_oriented_E8_sectors_cancel_local_and_global_chiral_anomaly(self):
        phase=sp.exp(-2*sp.pi*sp.I/3)
        self.assertEqual(sp.simplify(phase*sp.conjugate(phase)),1)
        self.assertEqual(self.r["opposite_chirality_E8_double_chiral_c_minus"],0)
        self.assertEqual(self.r["opposite_chirality_E8_double_partition"],"abs(chi_E8)^2")

    def test_three_same_chiral_E8_make_T_trivial_but_not_gravitationally_anomaly_free(self):
        self.assertEqual(self.r["three_same_chirality_E8_central_charge"],24)
        self.assertEqual(self.r["three_same_chirality_E8_modular_T_phase"],"1")
        self.assertFalse(self.r["three_same_chirality_E8_local_gravitational_anomaly_vanishes"])

    def test_standard_bosonic_E8_bulk_Chern_Simons_signature_matches_chiral_edge(self):
        self.assertEqual(self.r["bulk_boundary_chiral_c_minus"],8)
        self.assertEqual(self.r["source_required_inverse_T_holonomy"],"exp(+2*pi*i/3)")
        self.assertTrue(self.r["bulk_boundary_matching_gravitational_anomaly_inflow_known"])

    def test_anomaly_free_one_family_SM_I6_does_not_prove_e8_edge_I4(self):
        anomalies=standard_model_left_weyl_anomalies()
        self.assertTrue(anomalies["SM_I6_gauge_and_mixed_anomaly_zero"])
        self.assertFalse(anomalies["this_proves_chiral_E8_1_gravitational_I4_zero"])
        self.assertEqual(anomalies["SM_SU2_Witten_number_of_doublets"],4)

    def test_native_eta_bulk_topology_and_fourD_remain_unproved(self):
        self.assertFalse(self.r["source_actual_eta_Pfaffian_holonomy_measured"])
        self.assertFalse(self.r["original_2D_chiral_E8_worldsheet_identified"])
        self.assertFalse(self.r["original_2plus1D_E8_bulk_or_opposite_sector_derived"])
        self.assertFalse(self.r["four_dimensional_chiral_SM_index_from_this_double_proved"])
        self.assertFalse(self.r["actual_original_R113_absolute_all_field_Q0_computed"])
        self.assertFalse(self.r["F1_solved"])


if __name__=="__main__":
    unittest.main()
