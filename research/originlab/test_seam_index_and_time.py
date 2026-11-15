"""Finite seam inclusion and time no-go tests (not a physical R113 source)."""
import unittest
import sympy as sp


def units(n, i, j):
    return sp.Matrix(n, n, lambda r, c: int(r == i and c == j))


def twirl(a, generators):
    group = [sp.eye(a.rows)]
    for g in generators:
        group += [u*g for u in group[:]]
    return sp.simplify(sum((u*a*u.inv() for u in group),
                           sp.zeros(a.rows))/len(group))


class P1IndexAndDynamics(unittest.TestCase):
    def test_one_cell_index_four(self):
        x = sp.Matrix([[0, 1], [1, 0]])
        z = sp.diag(1, -1)
        a = sp.Matrix(2, 2, sp.symbols("a b c d"))
        self.assertEqual(twirl(a, [x, z]), sp.trace(a)*sp.eye(2)/2)
        reconstructed, idx = sp.zeros(2), sp.zeros(2)
        for i in range(2):
            for j in range(2):
                u = sp.sqrt(2)*units(2, i, j)
                v = sp.sqrt(2)*units(2, j, i)
                reconstructed += u*(sp.trace(v*a)/2)
                idx += u*v
        self.assertEqual((reconstructed-a).applyfunc(sp.simplify), sp.zeros(2))
        self.assertEqual(idx, 4*sp.eye(2))

    def test_single_cell_invariant_hamiltonian_is_scalar(self):
        x = sp.Matrix([[0, 1], [1, 0]])
        z = sp.diag(1, -1)
        i, sy = sp.eye(2), sp.Matrix([[0, -sp.I], [sp.I, 0]])
        a, b, c, d = sp.symbols("a b c d", real=True)
        h = a*i+b*x+c*sy+d*z
        self.assertEqual(sp.solve(list(h*x-x*h)+list(h*z-z*h),
                                  [b,c,d], dict=True),
                         [{b: 0, c: 0, d: 0}])
        # The sole density matrix commuting with both and of trace 1
        # is I/2; hence its modular operator is trivial.

    def test_two_cell_bell_ground_and_spectrum(self):
        x, z = sp.Matrix([[0, 1], [1, 0]]), sp.diag(1, -1)
        sx, sz = sp.kronecker_product(x, x), sp.kronecker_product(z, z)
        ident = sp.eye(4)
        self.assertEqual(sx*sz, sz*sx)
        vacuum = sp.Matrix([1, 0, 0, 1])/sp.sqrt(2)
        projector = (ident+sx)*(ident+sz)/4
        self.assertEqual(projector, vacuum*vacuum.T)
        self.assertEqual((ident-(sx+sz)/2).eigenvals(), {0: 1, 1: 2, 2: 1})
        a, b = sp.symbols("a b", positive=True)
        general = a*(ident-sx)/2+b*(ident-sz)/2
        self.assertEqual(set(general.eigenvals()), {0, a, b, a+b})

    def test_two_cell_bell_expectation_index_four(self):
        x, z = sp.Matrix([[0, 1], [1, 0]]), sp.diag(1, -1)
        sx, sz = sp.kronecker_product(x, x), sp.kronecker_product(z, z)
        bell = sp.Matrix.hstack(sp.Matrix([1, 0, 0, 1])/sp.sqrt(2),
                               sp.Matrix([1, 0, 0, -1])/sp.sqrt(2),
                               sp.Matrix([0, 1, 1, 0])/sp.sqrt(2),
                               sp.Matrix([0, 1, -1, 0])/sp.sqrt(2))
        a = sp.Matrix(4, 4, sp.symbols("a0:16"))
        projected = twirl(a, [sx, sz])
        in_basis = bell.T*projected*bell
        target = sp.diag(*[(bell.T*a*bell)[i,i] for i in range(4)])
        self.assertEqual((in_basis-target).applyfunc(sp.simplify), sp.zeros(4))
        reconstructed, idx = sp.zeros(4), sp.zeros(4)
        for i in range(4):
            for j in range(4):
                u = bell*units(4,i,j)*bell.T
                v = u.T
                reconstructed += u*twirl(v*a, [sx, sz])
                idx += u*v
        self.assertEqual((reconstructed-a).applyfunc(sp.simplify), sp.zeros(4))
        self.assertEqual(idx.applyfunc(sp.simplify), 4*sp.eye(4))

    def test_two_cell_twirl_is_not_finite_lindblad_step(self):
        x, z = sp.Matrix([[0, 1], [1, 0]]), sp.diag(1, -1)
        sx, sz = sp.kronecker_product(x, x), sp.kronecker_product(z, z)
        basis = [units(4,i,j) for i in range(4) for j in range(4)]
        columns = [sp.Matrix(twirl(b, [sx,sz])).reshape(16,1)
                   for b in basis]
        self.assertEqual(sp.Matrix.hstack(*columns).rank(), 4)
        # Every finite-dimensional exp(t*L) is invertible for finite t.
        # The twirl has rank 4 on the 16-dimensional operator space.

    def test_z12_graded_returns_depend_on_positive_generator(self):
        zeta = sp.exp(2*sp.pi*sp.I/12)
        u = sp.diag(1, zeta, zeta**3)
        a, b, c = units(3,1,0), units(3,2,1), units(3,0,2)
        for op, g in ((a,1),(b,2),(c,9)):
            self.assertEqual((u*op*u.inv()-zeta**g*op).applyfunc(sp.simplify),
                             sp.zeros(3))
        vacuum = sp.Matrix([1,0,0])
        returns = []
        for high in (2, 3):
            transfer = sp.diag(1, sp.Rational(1,2),
                               sp.Rational(1,2)**high)
            returns.append((vacuum.T*c*transfer*b*transfer*a*vacuum)[0])
        self.assertEqual(returns, [sp.Rational(1,8), sp.Rational(1,16)])


if __name__ == "__main__":
    unittest.main()
