"""Exact group check of marked v117; no fitted flavor or spatial dynamics."""
import sys
from pathlib import Path
from itertools import permutations
import sympy as s
sys.path.insert(0, str(Path('/Users/stefanhamann/Projekte/tfpt-theoryv4/verification')))
from v117_monodromy_weyl_a3 import M0_EXACT as M, U_EXACT as U, enumerate_group

I = s.eye(3)
simp = lambda x: x.applyfunc(s.simplify)
V = [tuple(x) for x in permutations(range(4))]

def aug(p):
    return s.Matrix([[int(i == p[j]) - int(i == p[3])
                      for j in range(3)] for i in range(3)])

# The marked order-4 action is a sign-twisted 4-cycle. Test an exact
# simultaneous intertwiner with the actual matrices, not a dimension match.
def rep(p):
    return int(s.combinatorics.Permutation(list(p)).signature()) * aug(p)

from sympy.combinatorics import Permutation
a = (1,2,3,0)
X = s.Matrix(3, 3, s.symbols('x:9'))
intertwiner = None
for b in V:
    if Permutation(list(b)).order() != 3:
        continue
    rels = list(U*X-X*rep(a))+list(M*X-X*rep(b))
    for z in s.linear_eq_to_matrix(rels, tuple(X))[0].nullspace():
        C = s.Matrix(3,3,z)
        if s.simplify(C.det()) != 0:
            intertwiner = C
            break
    if intertwiner is not None:
        break
assert intertwiner is not None
C = intertwiner
Cinv = C.inv()
G = {p: simp(C*rep(p)*Cinv) for p in V}
assert len({tuple(x) for x in G.values()}) == 24
assert {tuple(x) for x in G.values()} == {tuple(simp(x)) for x in enumerate_group([M,U])}
assert all(simp(g.H*g) == I for g in G.values())
assert all(g != g[0,0]*I for p,g in G.items() if p != tuple(range(4)))

for k, target in enumerate([(24,1),(6,2),(2,5),(1,9)]):
    subgroup = [g for p,g in G.items() if all(p[i] == i for i in range(k))]
    eq = [v for g in subgroup for v in list(g*X-X*g)]
    dim = 9 - s.linear_eq_to_matrix(eq, tuple(X))[0].rank()
    assert (len(subgroup),dim) == target
    print(f'{k} ordered marks: |stabilizer|={len(subgroup)}, commutant dimension={dim}')
    if k == 2:
        transposition = next(g for g in subgroup if g != I)
        sizes = sorted(transposition.eigenvals().values())
        assert sizes == [1,2]
        print('   common two-mark involution splits both quark sectors into 1+2')

for power, expected in [(1, (1,1,1)), (2, (1,2))]:
    assert tuple(sorted((U**power).eigenvals().values())) == expected
    print(f'clock power {power}: multiplicities {expected}')
oriented_pair_stabilizer = [p for p in V if p[0]==0 and p[1]==1
                            and Permutation(list(p)).signature()==1]
assert oriented_pair_stabilizer == [tuple(range(4))]
print('Two ordered marks plus an orientation have trivial stabilizer.')

def proj(v):
    return v*v.T/(v.T*v)[0]

r, q = s.Rational(2,3)**6, s.Rational(1,3)**6
T = proj(s.Matrix([1,1,1]))+r*proj(s.Matrix([1,-1,0]))+q*proj(s.Matrix([1,1,-2]))
assert sorted(T.eigenvals()) == [q,r,1]
assert simp(U*T-T*U) != s.zeros(3)
norm = s.simplify(s.trace((U*T-T*U).H*(U*T-T*U)))
assert norm > 0
# The exact one-step v330 covariance of the existing diagonal mark.
D = s.diag(1,-1,0)
v = s.Matrix([1,1,1])
value = s.simplify((v.T*D*T*D*v)[0]/3)
assert value == 2*r/3
print(f'v330: exact marked one-step covariance {value}; same-chart clock commutator norm^2 {norm}')
print('Distinct v330 eigenvalues forbid commutation with v117 in ANY equivalent irreducible three-dimensional chart.')
