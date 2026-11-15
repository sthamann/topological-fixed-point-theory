"""Exact negative controls for P1 seam and source hypotheses.

These are finite countermodels, not a reconstruction of the original
weighted R113 source, an E8 representation, or a full F1 theorem.
"""
import unittest
import sympy as s


class SeamSourceObstructions(unittest.TestCase):
    def test_reflection_polarization_and_quarter_return(self):
        r = s.Matrix([[0, 1], [1, 0]])
        e = s.diag(1, -1)
        i = s.eye(2)
        self.assertEqual(r*e+e*r, s.zeros(2))
        p, q = (i+e)/2, (i+r)/2
        self.assertEqual(p*q*p, p/2)
        p2, q2 = s.kronecker_product(p, p), s.kronecker_product(q, q)
        self.assertEqual(p2*q2*p2, p2/4)
        c = r*e
        self.assertEqual(c*c, -i)
        j = s.Matrix([[1, 2], [3, 4]])
        self.assertEqual(c*c*j*c.inv()**2, j)

    def test_reflection_covariant_spectral_flow(self):
        r = s.Matrix([[1, 0, 0], [0, 0, 1], [0, 1, 0]])
        t = s.symbols("t", real=True)
        b = lambda v: s.diag(v, 2, -2)
        self.assertEqual(r*b(t)*r, -b(-t))
        self.assertEqual(len(b(0).nullspace()), 1)
        # One level crosses 0 upward: spectral flow +1, eta jump 2.
        eta = lambda v: sum(s.sign(x) for x in b(v).diagonal())
        self.assertEqual(eta(1)-eta(-1), 2)

    def test_noncentral_square_gives_real_quarter_adjoint(self):
        u, parity = s.diag(1, s.I), s.diag(1, -1)
        creation = s.Matrix([[0, 0], [1, 0]])
        self.assertEqual(u*u, parity)
        self.assertEqual(u*creation*u.H, s.I*creation)
        self.assertEqual(u*creation.T*u.H, -s.I*creation.T)

    def test_positive_weight_not_os_positive(self):
        w = s.Matrix([[s.Rational(1,2), s.Rational(3,4)],
                      [s.Rational(3,4), s.Rational(1,2)]])
        self.assertTrue(all(0 < x <= 1 for x in w))
        v = s.Matrix([1, -1])
        self.assertEqual((v.T*w*v)[0], -s.Rational(1,2))

    def test_fiber_integration_is_not_block_restriction(self):
        x, y = s.symbols("x y", real=True)
        integral = s.integrate(s.exp(-(3*x*x+2*x*y+y*y)/2),
                              (y, -s.oo, s.oo))
        self.assertEqual(s.simplify(integral/(s.sqrt(2*s.pi)*s.exp(-x*x))), 1)

    def test_uniform_positive_noncompact_kernel_has_infinite_mass(self):
        t = s.symbols("t", positive=True)
        self.assertEqual(s.integrate(1, (t, 0, s.oo)), s.oo)


    def test_joint_z12_marked_grades_do_not_fix_transfer(self):
        # Latest marked word: 1+2+9 = 0 mod 12. Grades and the
        # vacuum return may agree while the positive transfer differs.
        z = s.exp(2*s.pi*s.I/12)
        phase = s.diag(1, z, z**3)
        a, b, c = s.zeros(3), s.zeros(3), s.zeros(3)
        a[1, 0], b[2, 1], c[0, 2] = 1, 1, 1
        for j, grade in ((a, 1), (b, 2), (c, 9)):
            self.assertEqual((phase*j*phase.inv()-z**grade*j).applyfunc(s.simplify), s.zeros(3))
        self.assertEqual(c*b*a, s.diag(1, 0, 0))
        vacuum = s.Matrix([1, 0, 0])
        results = []
        for high_energy in (2, 3):
            transfer = s.diag(1, s.Rational(1,2), s.Rational(1,2)**high_energy)
            results.append((vacuum.T*c*transfer*b*transfer*a*vacuum)[0])
        self.assertEqual(results, [s.Rational(1,8), s.Rational(1,16)])
        # This is a graded 3-state control, not a full affine E8 module.


if __name__ == "__main__":
    unittest.main()
