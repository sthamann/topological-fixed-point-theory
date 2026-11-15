"""Exact finite checks linking original v512 geometry to candidate Clifford time.

These controls do not supply the original all-field R113 Q0 or a 3+1D parent.
"""
import itertools
import unittest
import sympy as sp


def cross_ratio(pts):
    a, b, c, d = pts
    return sp.simplify((a-c)*(b-d)/((a-d)*(b-c)))


class P1GeometryTimeBridge(unittest.TestCase):
    def test_square_has_d4_not_a4_projective_permutations(self):
        pts = [sp.Integer(1), -sp.Integer(1), sp.I, -sp.I]
        self.assertEqual(cross_ratio(pts), -1)
        stabilizer = [perm for perm in itertools.permutations(range(4))
                      if cross_ratio([pts[i] for i in perm]) == -1]
        self.assertEqual(len(stabilizer), 8)
        self.assertFalse(any(sum(perm[i] == i for i in range(4)) == 1
                             for perm in stabilizer))
        cw = (3+4*sp.I)/5
        pts = [sp.Integer(1), -sp.Integer(1), cw, -cw]
        value = cross_ratio(pts)
        self.assertEqual(sum(cross_ratio([pts[i] for i in perm]) == value
                             for perm in itertools.permutations(range(4))), 4)

    def test_fourier_volume_selects_square_in_v512_family(self):
        u = sp.symbols("u", nonzero=True)
        points = [1, -1, u, -u]
        v = sp.Matrix([[point**k for k in range(4)] for point in points])
        self.assertEqual(sp.simplify(v.det()-4*u*(u*u-1)**2), 0)
        self.assertEqual(sp.simplify((v.subs(u,sp.I)/2)*
                                     (v.subs(u,sp.I)/2).H), sp.eye(4))
        cw = (3+4*sp.I)/5
        determinant_squared = sp.simplify(sp.expand_complex(
            v.det().subs(u,cw)*sp.conjugate(v.det().subs(u,cw))))
        self.assertEqual(determinant_squared/256, sp.Rational(256,625))
        delta_odd = 2*sp.log(2)/sp.pi
        self.assertEqual(sp.simplify(sp.cosh(sp.pi*delta_odd/2)**(-4)),
                         sp.Rational(256,625))

    def test_weighted_isometry_does_not_select_square(self):
        # Explicit original v512 non-square witness, rational complex entries.
        u = (3+4*sp.I)/5
        pts = [1, -1, u, -u]
        f = sp.Matrix([[p**k for k in range(4)] for p in pts])/2
        self.assertEqual(sp.simplify((f*f.H).det()), sp.Rational(256,625))
        # Every invertible sampling map is an isometry for its pullback
        # metric. Its existence is therefore not a physical modulus selector.
        metric = (f*f.H).inv()
        self.assertEqual((metric-metric.H).applyfunc(sp.simplify), sp.zeros(4))
        self.assertEqual(metric.eigenvals(), {sp.Rational(5,2):2,
                                              sp.Rational(5,8):2})
        self.assertEqual((f.H*metric*f).applyfunc(sp.simplify), sp.eye(4))
        self.assertNotEqual((f*f.H).applyfunc(sp.simplify), sp.eye(4))

    def test_square_laplacian_intertwines_with_bell_clifford(self):
        x = sp.Matrix([[0,1],[1,0]])
        z = sp.diag(1,-1)
        ident = sp.eye(4)
        ta = sp.kronecker_product(x,sp.eye(2))
        tb = sp.kronecker_product(sp.eye(2),x)
        sx = sp.kronecker_product(x,x)
        sz = sp.kronecker_product(z,z)
        h = (sp.Matrix([[1,1],[1,-1]])/sp.sqrt(2))
        f = sp.kronecker_product(h,h)
        bell = sp.Matrix.hstack(sp.Matrix([1,0,0,1])/sp.sqrt(2),
                                sp.Matrix([0,1,1,0])/sp.sqrt(2),
                                sp.Matrix([1,0,0,-1])/sp.sqrt(2),
                                sp.Matrix([0,1,-1,0])/sp.sqrt(2))
        w = bell*f.T
        self.assertEqual(sp.simplify(w*ta*w.H),sx)
        self.assertEqual(sp.simplify(w*tb*w.H),sz)
        square = ident-(ta+tb)/2
        self.assertEqual(square.eigenvals(),{0:1,1:2,2:1})
        self.assertEqual(sp.simplify(w*square*w.H),ident-(sx+sz)/2)
        proj = {(a,b):(ident+a*ta)*(ident+b*tb)/4
                for a,b in itertools.product((-1,1),repeat=2)}
        transfer = proj[1,1]+(proj[1,-1]+proj[-1,1])/2+proj[-1,-1]/4
        self.assertEqual(transfer,sp.Matrix([[9,3,3,1],[3,9,1,3],
                                              [3,1,9,3],[1,3,3,9]])/16)

    def test_square_symmetry_without_locality_keeps_free_jump(self):
        x = sp.Matrix([[0,1],[1,0]])
        ta = sp.kronecker_product(x,sp.eye(2))
        tb = sp.kronecker_product(sp.eye(2),x)
        a,b = sp.symbols("a b",positive=True)
        h = a*(2*sp.eye(4)-ta-tb)+b*(sp.eye(4)-ta*tb)
        l = sp.Symbol("lambda")
        self.assertEqual(sp.factor(h.charpoly(l).as_expr()),
                         -l*(4*a-l)*(2*a+2*b-l)**2)

    def test_puncture_and_tetrahedron_cycles_have_different_character(self):
        edges = list(itertools.combinations(range(4),2))
        boundary = sp.Matrix(4,len(edges),lambda v,j:
                             -1 if v==edges[j][0] else
                             1 if v==edges[j][1] else 0)
        cycles = sp.Matrix.hstack(*boundary.nullspace())
        perm = (1,2,3,0)
        vp = sp.zeros(4)
        ep = sp.zeros(6)
        for i in range(4):
            vp[perm[i],i] = 1
        for j,(a,b) in enumerate(edges):
            aa,bb = perm[a],perm[b]
            ep[edges.index(tuple(sorted((aa,bb)))),j] = 1 if aa<bb else -1
        standard = sp.Matrix.hstack(*(sp.eye(4)[:,i]-sp.eye(4)[:,3]
                                      for i in range(3)))
        a = sp.Matrix.hstack(*(standard.gauss_jordan_solve(vp*standard[:,i])[0]
                               for i in range(3)))
        b = sp.Matrix.hstack(*(cycles.gauss_jordan_solve(ep*cycles[:,i])[0]
                               for i in range(3)))
        self.assertEqual((sp.trace(a),sp.trace(b)),(-1,1))

    def test_actual_free_collar_covariance_is_not_c4_heat(self):
        q = sp.Rational(1,2)
        c0 = -sp.log((1-q)**2)/(2*sp.pi)
        c1 = -sp.log(1+q*q)/(2*sp.pi)
        c2 = -sp.log((1+q)**2)/(2*sp.pi)
        g = sp.Matrix([[c0,c1,c2,c1],[c1,c0,c1,c2],
                       [c2,c1,c0,c1],[c1,c2,c1,c0]])
        samples = [(sp.ones(4,1),sp.log(sp.Rational(16,15))/sp.pi),
                   (sp.Matrix([1,sp.I,-1,-sp.I]),sp.log(3)/sp.pi),
                   (sp.Matrix([1,-1,1,-1]),sp.log(sp.Rational(5,3))/sp.pi)]
        for vec,eigenvalue in samples:
            self.assertEqual(sp.simplify(sp.expand_log((g*vec)[0],force=True)
                                         -eigenvalue),0)
        self.assertGreater(float(sp.log(3)/sp.log(sp.Rational(16,15))),1)


if __name__ == "__main__":
    unittest.main()
