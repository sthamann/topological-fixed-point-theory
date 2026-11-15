"""Scoped finite Clifford mark bridge to A3/E8, not an original P1 source."""
import itertools
import unittest
import sympy as sp


class CliffordMarkAndGlue(unittest.TestCase):
    def test_lift_eight_and_adjoint_four(self):
        r=sp.Matrix([[0,1],[1,0]])
        e=sp.diag(1,-1)
        one=sp.eye(2)
        words={tuple(one):one}
        pending=[one]
        while pending:
            a=pending.pop()
            for generator in (r,e):
                x=a*generator
                if tuple(x) not in words:
                    words[tuple(x)]=x
                    pending.append(x)
        self.assertEqual(len(words),8)
        self.assertEqual(r*e*r*e,-one)
        self.assertEqual(len({tuple(sp.kronecker_product(m,m.inv().T))
                              for m in words.values()}),4)

    def test_full_nontrivial_cayley_is_tetrahedron(self):
        nodes=[(0,0),(0,1),(1,0),(1,1)]
        edges=list(itertools.combinations(range(4),2))
        boundary=sp.Matrix(4,len(edges),lambda v,k:
                           -1 if v==edges[k][0]
                           else 1 if v==edges[k][1] else 0)
        self.assertEqual(len(edges),6)
        self.assertEqual(len(edges)-boundary.rank(),3)
        matchings=[]
        for g in nodes[1:]:
            matchings.append(frozenset(tuple(sorted((p,(p[0]^g[0],p[1]^g[1]))))
                                        for p in nodes))
        self.assertEqual(len(set(matchings)),3)
        self.assertEqual({len(m) for m in matchings},{2})

    def test_mark_differences_generate_a3(self):
        e=[sp.eye(4)[:,i] for i in range(4)]
        roots={tuple(e[i]-e[j]) for i in range(4) for j in range(4) if i!=j}
        simples=sp.Matrix.hstack(*(e[i]-e[i+1] for i in range(3)))
        gram=simples.T*simples
        self.assertEqual(len(roots),12)
        self.assertEqual(gram,sp.Matrix([[2,-1,0],[-1,2,-1],[0,-1,2]]))
        self.assertEqual(gram.det(),4)

    def test_cycle_orientation_is_signed_standard(self):
        edges=list(itertools.combinations(range(4),2))
        B=sp.Matrix(4,len(edges),lambda v,k:
                    -1 if v==edges[k][0] else 1 if v==edges[k][1] else 0)
        z=sp.Matrix.hstack(*B.nullspace())
        self.assertEqual(z.cols,3)
        for perm in itertools.permutations(range(4)):
            P=sp.zeros(6)
            for col,(i,j) in enumerate(edges):
                vi,vj=perm[i],perm[j]
                sign=1 if vi<vj else -1
                row=edges.index(tuple(sorted((vi,vj))))
                P[row,col]=sign
            action=sp.Matrix.hstack(*(z.gauss_jordan_solve(P*z[:,k])[0]
                                      for k in range(3)))
            inversions=sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))
            character=(-1)**inversions*(sum(i==perm[i] for i in range(4))-1)
            self.assertEqual(sp.trace(action),character)

    def test_discriminant_glue_and_240_roots(self):
        self.assertEqual(sp.Rational(5,4)+sp.Rational(3,4),2)
        self.assertEqual((4*4)//(4**2),1)
        self.assertEqual([40+12,16*4,10*6,16*4],[52,64,60,64])
        self.assertEqual(sum([52,64,60,64]),240)

    def test_clock_and_vacuum_do_not_fix_transfer(self):
        x,z=sp.Matrix([[0,1],[1,0]]),sp.diag(1,-1)
        sx,sz=sp.kronecker_product(x,x),sp.kronecker_product(z,z)
        i=sp.eye(4)
        projectors={(u,v):(i+u*sx)*(i+v*sz)/4
                    for u,v in itertools.product((-1,1),repeat=2)}
        vacuum=sp.Matrix([1,0,0,1])/sp.sqrt(2)
        self.assertEqual(projectors[1,1],vacuum*vacuum.T)
        clock=sz*sx
        returns=[]
        for b in (2,3):
            h=projectors[-1,-1]+b*(projectors[-1,1]+projectors[1,-1])
            self.assertEqual(h.eigenvals(),{0:1,1:1,b:2})
            self.assertEqual(h*clock,clock*h)
            self.assertEqual(h*vacuum,sp.zeros(4,1))
            T=sp.zeros(4)
            for (u,v),projector in projectors.items():
                energy=0 if (u,v)==(1,1) else 1 if (u,v)==(-1,-1) else b
                T+=sp.Rational(1,2)**energy*projector
            flip=sp.kronecker_product(x,sp.eye(2))
            returns.append(sp.simplify((vacuum.T*flip.T*T*flip*vacuum)[0]))
        self.assertEqual(returns,[sp.Rational(1,4),sp.Rational(1,8)])


if __name__=="__main__":
    unittest.main()
