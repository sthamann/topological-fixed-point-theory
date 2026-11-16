"""Exact and negative controls for hypothetical spectral P1 deck reader."""
import math
import unittest

from p1_torus_deck_spectral_reader import (
    torus_heat, branched_even_odd, branched_double_branch_count,
    branched_scalar_uv_trace, heat_residue_estimate,
    klein_glide_trace, klein_glide_poisson, spin_pi_lift,
    spin_structure_inversion_trace, two_sheet_branch_points, controls, overlapped_elliptic_covers_genus,
    diagonal_affine_level, parity_twirl_qubit, dephasing_qubit,
    dephasing_plus_state_eigenvalues,
)


class TorusDeckReaderTests(unittest.TestCase):
    def test_four_points_not_initial_marks(self):
        self.assertEqual(set(two_sheet_branch_points()), {
            (0.,0.), (.5,0.), (0.,.5), (.5,.5)
        })
        self.assertEqual(branched_double_branch_count(1), 4)

    def test_exact_branched_even_minus_odd_is_one_at_every_cutoff(self):
        for aspect in (1., 1.7):
            for t in (.005, .02, .1):
                even, odd = branched_even_odd(t, aspect)
                self.assertAlmostEqual(even-odd, 1., places=12)
                self.assertGreaterEqual(odd, 0.)

    def test_half_heat_residue_but_not_square_selector(self):
        target = 1/(8*math.pi)
        for aspect in (1., 1.7):
            self.assertAlmostEqual(heat_residue_estimate(.005, aspect),
                                   target, delta=2e-11)

    def test_poisson_residue_survives_small_t_beyond_fixed_fourier_cutoff(self):
        target = 1/(8*math.pi)
        for t in (1e-8, 1e-10):
            for aspect in (1., 1.7):
                self.assertAlmostEqual(heat_residue_estimate(t, aspect),
                                       target, delta=2e-15)
                finite_residue = t*torus_heat(t, aspect, cutoff=40)/2
                self.assertLess(finite_residue, target/100)
                # Every finite-mode summand is <=1, so this upper bound
                # goes to zero when t->0 at a fixed Fourier cutoff.
                self.assertLessEqual(finite_residue, t*(2*40+1)**2/2)

    def test_poisson_residue_matches_resolved_direct_fourier_sums(self):
        for aspect in (.6, 1.7, 3.):
            for t in (.002, .02, .2):
                direct = t*torus_heat(t, aspect, cutoff=100)/2
                self.assertAlmostEqual(heat_residue_estimate(t, aspect),
                                       direct, delta=2e-13)

    def test_free_klein_cover_has_different_uv_signature(self):
        for t in (.005, .02, .07):
            direct = klein_glide_trace(t, 1.7)
            poisson = klein_glide_poisson(t, 1.7)
            self.assertAlmostEqual(direct, poisson, delta=2e-11)
        self.assertLess(klein_glide_poisson(.005, 1.7), .004)
        self.assertGreater(klein_glide_poisson(.02, 1.7), .48)

    def test_scale_independent_scalar_deck_trace_is_genus_one_neccessary(self):
        # For branched degree-two S² covers: UV fixed-point index=(g+1)/2,
        # while IR index=1 from the unique constant scalar zero mode.
        self.assertEqual(branched_scalar_uv_trace(1),1.)
        self.assertEqual(branched_scalar_uv_trace(0),.5)
        self.assertEqual(branched_scalar_uv_trace(2),1.5)
        self.assertEqual(branched_scalar_uv_trace(3),2.)

    def test_spin_lift_has_order_four_not_two(self):
        phase=spin_pi_lift()
        self.assertAlmostEqual(phase.real, 0., places=14)
        self.assertAlmostEqual(phase.imag, 1., places=14)
        self.assertAlmostEqual((phase**2).real, -1., places=14)
        self.assertAlmostEqual((phase**4).real, 1., places=14)

    def test_only_periodic_spin_structure_has_unpaired_fourier_zero_mode(self):
        traces=[spin_structure_inversion_trace(a,b)
                for a in (0.,.5) for b in (0.,.5)]
        self.assertCountEqual(traces,[1,0,0,0])

    def test_two_independent_elliptic_overlays_are_generically_higher_genus(self):
        self.assertEqual([overlapped_elliptic_covers_genus(k)
                          for k in range(5)], [5,4,3,2,None])

    def test_larger_E8_parent_exchange_diagonal_is_not_level_one(self):
        self.assertEqual(diagonal_affine_level((1,1)),2)
        self.assertEqual(diagonal_affine_level((1,1,1)),3)
        self.assertEqual(diagonal_affine_level((1,)),1)

    def test_holographic_twirl_is_noninjective_on_unrestricted_states(self):
        plus=((.5,.5),(.5,.5))
        minus=((.5,-.5),(-.5,.5))
        self.assertNotEqual(plus,minus)
        self.assertEqual(parity_twirl_qubit(plus),parity_twirl_qubit(minus))

    def test_deck_derived_dephasing_is_CPTP_semigroup_but_leaves_rate_free(self):
        rho=((.5,.5),(.5,.5))
        for t in (0.,.3,1.,4.):
            eig=dephasing_plus_state_eigenvalues(t,2.)
            self.assertGreaterEqual(min(eig),0.)
            self.assertAlmostEqual(sum(eig),1.)
            self.assertAlmostEqual(dephasing_qubit(rho,t,2.)[0][1],
                                   math.exp(-2*t)/2)
        r1=dephasing_qubit(dephasing_qubit(rho,.2,2.),.3,2.)
        r2=dephasing_qubit(rho,.5,2.)
        self.assertAlmostEqual(r1[0][1],r2[0][1])

    def test_no_original_source_or_spacetime_is_claimed(self):
        data=controls()
        self.assertFalse(data["original_all_field_Q0_found"])
        self.assertFalse(data["physical_OS_time_recovered"])
        self.assertFalse(data["F1_solved"])

if __name__ == "__main__":
    unittest.main()
