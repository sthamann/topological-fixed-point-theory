import unittest
import sympy as s
from boundary_amplitude import state_from_projection,oriented_amplitude,graph_projection,verify_graph_route
from source_sewing import source_lift


class BoundaryAmplitudeTests(unittest.TestCase):
    def test_two_mode_source_graph_recovers_product_response(self):
        q,receipt=verify_graph_route(s.eye(2)/4)
        self.assertEqual(q,s.diag(1,s.Rational(1,4),s.Rational(1,4),s.Rational(1,16)))
        self.assertTrue(receipt['state_ray_unique_from_CAR_projection'])
        self.assertFalse(receipt['operator_scalar_recovered'])

    def test_complex_nondiagonal_source_preserves_orientation_and_phase(self):
        r=s.Matrix([[s.I/2,s.Rational(1,3)],[s.Rational(1,5),-s.I/4]])
        self.assertEqual(oriented_amplitude(graph_projection(r))[0],source_lift(r))

    def test_mixed_covariance_does_not_get_pretended_pure(self):
        with self.assertRaises(ValueError):state_from_projection(s.eye(2)/2)

    def test_quarter_phase_survives_field_products(self):
        q,_=verify_graph_route(s.I*s.eye(2))
        self.assertEqual(q,s.diag(1,s.I,s.I,-1))

    def test_labeled_face_changes_transport_input_to_output(self):
        r=s.Matrix([[s.I/2,s.Rational(1,3)],[s.Rational(1,5),-s.I/4]])
        u=s.Matrix([[0,s.I],[1,0]]);v=s.diag(-s.I,1)
        g=graph_projection(r);w=s.diag(u,v)
        q=oriented_amplitude(g)[0]
        actual=oriented_amplitude(w*g*w.adjoint())[0]
        self.assertEqual(actual,source_lift(v)*q*source_lift(u.adjoint()))
        self.assertNotEqual(actual,source_lift(u)*q*source_lift(v.adjoint()))

    def test_cauchy_projection_does_not_imply_positive_amplitude(self):
        q=oriented_amplitude(graph_projection(s.Matrix([[-1]])))[0]
        self.assertEqual(q,s.diag(1,-1))


if __name__=='__main__':unittest.main()
