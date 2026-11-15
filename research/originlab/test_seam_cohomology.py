"""Exact independent source formula and D4 orientation controls."""
import unittest
import sympy as sp
from seam_cohomology import controls


class SeamCohomologyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = controls()

    def test_original_rational_sections_transform_as_claimed(self):
        z = sp.symbols("z")
        f = [z**k/(z**4-1) for k in range(5)]
        for k in range(5):
            self.assertEqual(sp.simplify(f[k].subs(z,sp.I*z)-sp.I**k*f[k]),0)
            self.assertEqual(sp.simplify(f[k].subs(z,1/z)+f[4-k]),0)
        self.assertEqual(sp.simplify(f[4]-f[0]),1)

    def test_log_differential_pullback_has_opposite_reflection_sign(self):
        z = sp.symbols("z")
        omega = [z**j/(z**4-1) for j in range(3)]
        for j in range(3):
            # Pulling back a 1-form multiplies its coefficient by dz'/dz.
            self.assertEqual(sp.simplify(sp.I*omega[j].subs(z,sp.I*z)
                                         -sp.I**(j+1)*omega[j]),0)
            self.assertEqual(sp.simplify(-omega[j].subs(z,1/z)/z**2
                                         -omega[2-j]),0)

    def test_two_fixed_points_give_equivariant_exact_sequence(self):
        r=sp.diag(1,sp.I,-1,-sp.I,1)
        s=-sp.Matrix(5,5,lambda i,j:int(i+j==4))
        er=sp.eye(2)
        es=sp.Matrix([[0,1],[1,0]])
        ev=sp.Matrix([[-1,0,0,0,0],[0,0,0,0,1]])
        inc=sp.Matrix(5,3,lambda i,j:int(i==j+1))
        fr=sp.diag(sp.I,-1,-sp.I)
        fs=sp.Matrix(3,3,lambda i,j:int(i+j==2))
        self.assertEqual((ev.rank(),inc.rank(),(ev*inc).rank()),(2,3,0))
        self.assertEqual(ev*r,er*ev)
        self.assertEqual(ev*s,es*ev)
        self.assertEqual(r*inc,inc*fr)
        self.assertEqual(s*inc,-inc*fs)
        self.assertTrue(self.report["D4_orientation_twisted_three_plus_two_exact_sequence"])
        self.assertEqual(self.report["fixed_point_evaluation_rank"],2)

    def test_new_character_five_to_three_plus_two_decomposition(self):
        result = self.report
        self.assertEqual(result["five_mode_c4_multiplicities"],
                         {"1":2,"i":1,"-1":1,"-i":1})
        self.assertEqual(result["log_family_c4_multiplicities"],
                         {"1":0,"i":1,"-1":1,"-i":1})
        self.assertEqual(result["minimal_distinct_character_compensator_dimension"],3)

    def test_D4_geometry_vacuum_gap_and_trace_do_not_fix_time(self):
        r=self.report
        self.assertEqual(r["D4_irreducible_carrier_dimensions"],[1,1,1,2])
        self.assertEqual(r["D4_invariant_generator_excited_free_parameters"],3)
        self.assertEqual(r["matched_gap_and_trace_witness_spectra"],
                         [["0","1","2","3","3"],["0","1","2","2","4"]])
        self.assertEqual(r["same_gap_trace_distinguishing_return_at_log2"],
                         ["1/4","1/16"])

    def test_orientation_is_a_nontrivial_physical_gate(self):
        result=self.report
        self.assertEqual((result["cyclic_neutral_compensated_rank"],
                          result["full_dihedral_invariant_rank"],
                          result["orientation_twisted_full_invariant_rank"]),
                         (3,1,2))
        self.assertEqual(result["uncompensated_colour_commutant_complex_dimension"],3)
        self.assertEqual(result["uncompensated_colour_full_d4_commutant_complex_dimension"],2)
        self.assertEqual(result["carrier_quarter_commutant_complex_dimension"],7)

    def test_does_not_promote_origin_or_r113(self):
        result=self.report
        self.assertFalse(result["F1_solved"])
        self.assertFalse(result["source_orientation_and_eta_identified"])
        self.assertFalse(result["original_weighted_Q0_produced"])
        self.assertFalse(result["original_R113_values_produced"])
        self.assertFalse(result["physical_SU3_SU2_identification_proved"])


if __name__=="__main__":
    unittest.main()
