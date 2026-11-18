import unittest
import sympy as s
from source_sewing import source_lift,determinant_lift,covariance_residual,creator,checks


class SewingTests(unittest.TestCase):
    def test_noncommuting_chronological_composition(self):
        self.assertTrue(checks()['noncommuting_transfer_composition_exact'])

    def test_recurrence_agrees_with_all_three_mode_minors(self):
        r=s.Matrix([[s.Rational(1,2),s.I/7,0],[-s.I/7,s.Rational(1,3),s.Rational(1,11)],[0,s.Rational(1,11),s.Rational(1,4)]])
        self.assertEqual(source_lift(r),determinant_lift(r))

    def test_covariance_is_stronger_than_vacuum_intertwining(self):
        r=s.diag(s.Rational(1,4),s.Rational(1,4))
        wrong=s.diag(1,s.Rational(1,4),s.Rational(1,4),s.Rational(3,4))
        for j in range(2):
            residual=wrong*creator(2,j)-s.Rational(1,4)*creator(2,j)*wrong
            self.assertEqual(residual[:,0],s.zeros(4,1))
        self.assertEqual(covariance_residual(wrong,r),[0,1])

    def test_euclidean_covariance_is_not_naive_star_automorphism(self):
        a=s.Rational(1,4);q=source_lift(s.diag(a,a));c=creator(2,0).adjoint()
        self.assertEqual(c*q,q*a*c)
        self.assertNotEqual(q*c,a*c*q)

    def test_signed_collar_needs_explicit_sea_relative_normalization(self):
        raw=source_lift(s.diag(s.Rational(1,4),64))
        negative_creator=creator(2,1)
        particle_hole=negative_creator+negative_creator.adjoint()
        relative=particle_hole*raw*particle_hole.adjoint()/64
        self.assertEqual(relative,source_lift(s.diag(s.Rational(1,4),s.Rational(1,64))))
        self.assertEqual(max(raw.diagonal()),64)
        self.assertEqual(max(relative.diagonal()),1)


if __name__=='__main__':unittest.main()
