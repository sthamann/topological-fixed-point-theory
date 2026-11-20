"""An exact operational nonidentifiability test. Run with Python 3 and SymPy."""

from itertools import combinations, product
from functools import lru_cache
from math import factorial

import sympy as sp


def axes(dimension):
    basis = list(sp.eye(dimension).columnspace())
    return basis + [left + right for left, right in combinations(basis, 2)]


def response_rank(dimension):
    # Lengths cancel, so the chosen probe vectors need no square roots.
    probes = axes(dimension)
    kernel = sp.Matrix([
        [left.dot(right)**2 / (left.dot(left)*right.dot(right))
         - sp.Rational(1, dimension) for right in probes]
        for left in probes
    ])
    return kernel.rank()


def expectations(alice, bob, state):
    return sp.Matrix([
        [(state.T * sp.kronecker_product(a, b) * state)[0] for b in bob]
        for a in alice
    ])


def tomography():
    n, r = sp.symbols("n r", positive=True, integer=True)

    def dim(size, scalar):
        return size + scalar*size*(size-1)/2

    gap = sp.factor(dim(n, r)**2 - dim(n**2, r))
    expected = r*(r-2)*(n*(n-1)/2)**2
    assert sp.simplify(gap - expected) == 0
    # The quaternion entry compares matrix ranks only, without positing a
    # canonical quaternion tensor product.
    counts = [(k, int(dim(2, k)**2), int(dim(4, k))) for k in (1, 2, 4)]
    assert counts == [(1, 9, 10), (2, 16, 16), (4, 36, 28)]
    return gap, counts


def hidden_real_coordinate():
    # A legitimate real bipartite observable invisible to all local
    # symmetric-real product measurements.
    one = sp.eye(2)
    x = sp.Matrix([[0, 1], [1, 0]])
    z = sp.diag(1, -1)
    j = x*z
    hidden = sp.kronecker_product(j, j)
    states = [(sp.eye(4) + s*hidden)/4 for s in (-1, 1)]
    for rho in states:
        assert rho.trace() == 1 and rho.eigenvals() == {0: 2, sp.Rational(1, 2): 2}
    for a, b in product((one, x, z), repeat=2):
        assert ((states[0]-states[1])*sp.kronecker_product(a, b)).trace() == 0
    return [(rho*hidden).trace() for rho in states]


def outcomes_a3():
    # In the vectorized Bell frame the Bell state is |I>/sqrt(2).
    # With the normal as z, axial conjugation fixes I and Z and assigns
    # opposite weights to the two transition operators.
    e12 = sp.Matrix([[0, 1], [0, 0]])
    operators = [sp.eye(2), sp.diag(1, -1), e12, e12.T]
    vectors = [a.reshape(4, 1) for a in operators]
    projectors = [v*v.H/(v.H*v)[0] for v in vectors]
    assert sum(projectors, sp.zeros(4)) == sp.eye(4)
    assert sp.Matrix([[(p*q).trace() for q in projectors] for p in projectors]) == sp.eye(4)
    roots = [projectors[i]-projectors[i+1] for i in range(3)]
    cartan = sp.Matrix([[(a*b).trace() for b in roots] for a in roots])
    assert cartan == sp.Matrix([[2, -1, 0], [-1, 2, -1], [0, -1, 2]])
    root_count = sum(sum(t)==0 and sum(x*x for x in t)==2
                     for t in product((-1, 0, 1), repeat=4))
    assert root_count == 12 and cartan.det() == 4
    # R represents a spatial quarter turn up to an irrelevant phase.
    quarter = sp.diag(1, sp.I)
    adjoint = sp.kronecker_product(quarter, quarter.conjugate())
    assert all(adjoint*p*adjoint.H == p for p in projectors)
    return cartan, root_count


@lru_cache(maxsize=None)
def labelled_graph_count(degrees):
    """Count all simple labelled graphs with this prescribed degree list."""
    if not degrees or degrees[0] == 0:
        return 1
    degree, others = degrees[0], degrees[1:]
    valid = [i for i, d in enumerate(others) if d > 0]
    count = 0
    for near in combinations(valid, degree):
        remaining = list(others)
        for i in near:
            remaining[i] -= 1
        # Sorting only memoizes label permutations; each choice of labelled
        # neighbors still contributes once.
        count += labelled_graph_count(tuple(sorted(remaining, reverse=True)))
    return count


def dodecahedral_control():
    # Added choice of finite probe design, independent of the Bell optimum.
    # The convention of the labelled Atlas comparison is + for meeting
    # K5 edges and - for disjoint edges.
    golden = (1+sp.sqrt(5))/2
    probes = [sp.Matrix([1, a, b]) for a, b in product((-1, 1), repeat=2)]
    for sign in (-1, 1):
        coords = [0, golden-1, sign*golden]
        for i in range(3):
            probes.append(sp.Matrix(coords[i:]+coords[:i]))
    assert len(probes) == 10 and all(sp.simplify(p.dot(p)-3) == 0 for p in probes)
    cos2 = sp.Matrix([[sp.expand(left.dot(right)**2)/9 for right in probes]
                      for left in probes])
    assert set(cos2) == {1, sp.Rational(1, 9), sp.Rational(5, 9)}
    gram = (3*cos2-sp.ones(10))/2
    edges = {(a, b) for a, b in combinations(range(10), 2)
             if cos2[a, b] == sp.Rational(5, 9)}
    # The five independent four-sets give the canonical Petersen K5 labels:
    # each probe belongs to two of the five sets. No graph fit is needed.
    tetrahedra = [q for q in combinations(range(10), 4)
                  if all(e not in edges for e in combinations(q, 2))]
    assert len(tetrahedra) == 5
    labels = [{i for i, t in enumerate(tetrahedra) if k in t} for k in range(10)]
    assert len(set(map(frozenset, labels))) == 10 and all(len(e) == 2 for e in labels)
    conference = sp.Matrix(10, 10, lambda a, b: 0 if a == b else
                           (1 if labels[a] & labels[b] else -1))
    atlas = sp.eye(10)+conference/3
    assert gram == sp.eye(10)-conference/3
    assert atlas*gram == sp.zeros(10) and atlas.rank() == gram.rank() == 5
    assert atlas+gram == 2*sp.eye(10)
    graph = sp.Matrix(10, 10, lambda a, b: int((min(a,b), max(a,b)) in edges))
    assert graph**2 == 2*sp.eye(10)-graph+sp.ones(10)
    assert gram*sp.ones(10, 1) == sp.zeros(10, 1) and gram**2 == 2*gram
    cubic = labelled_graph_count((3,)*10)
    # Each graph automorphism permutes the five distinguished four-sets
    # faithfully; conversely all five-set permutations act on their pairs.
    petersen_count = factorial(10)//factorial(5)
    assert (cubic, petersen_count) == (11180820, 30240)
    print("Fixed-degree labelled graph baseline:", petersen_count, "/", cubic)
    return len(edges), len(tetrahedra), (atlas*gram).rank()


def rotation_weights():
    x, y, z, t = sp.symbols("x y z t")
    basis = (x*x - y*y, 2*z*z - x*x - y*y, x*y, x*z, y*z)
    # J = x*d_y - y*d_x. The half turn sends x,y to -x,-y.
    coefficients = sp.symbols("c:5")
    template = sum(c * f for c, f in zip(coefficients, basis))

    def in_basis(f):
        equations = sp.Poly(template-f, x, y, z).coeffs()
        solution = sp.solve(equations, coefficients, dict=True)[0]
        return sp.Matrix([solution[c] for c in coefficients])

    generator = sp.Matrix.hstack(*[
        in_basis(x*sp.diff(f, y) - y*sp.diff(f, x)) for f in basis
    ])
    turn = sp.Matrix.hstack(*[
        in_basis(f.subs({x: -x, y: -y}, simultaneous=True)) for f in basis
    ])
    polynomial = sp.factor(generator.charpoly(t).as_expr())
    assert polynomial == t*(t*t + 1)*(t*t + 4)
    assert turn.eigenvals() == {1: 3, -1: 2}
    return polynomial, turn.eigenvals()


def main():
    gap, counts = tomography()
    print("Local-product dimension minus uniform composite dimension:", gap)
    print("Scalar dimension, local-product count, joint count:", counts)
    print("Real locally indistinguishable states; global witness:", hidden_real_coordinate())
    gram, root_count = outcomes_a3()
    print("Bell outcome difference lattice:", gram.tolist(),
          "; determinant =", gram.det(), "; roots =", root_count)
    print("Axial quarter turn fixes all four outcome projectors.")
    edges, tetrahedra, cross_rank = dodecahedral_control()
    print("Dodecahedral control: Petersen edges =", edges, "; independent four-sets =",
          tetrahedra, "; rank of product with labelled Atlas Gram =", cross_rank)
    x = sp.Matrix([[0, 1], [1, 0]])
    z = sp.diag(1, -1)
    psi = sp.Matrix([1, 0, 0, 1])/sp.sqrt(2)
    cor = expectations([z, x], [(z+x)/sp.sqrt(2), (z-x)/sp.sqrt(2)], psi)
    bell = cor[0, 0] + cor[0, 1] + cor[1, 0] - cor[1, 1]
    assert bell == 2*sp.sqrt(2)
    assert cor.applyfunc(lambda a: a*a) == sp.ones(2)/2
    print("Both fields: CHSH =", bell, "; squared table =", cor.applyfunc(lambda a: a*a))

    ranks = {d: response_rank(d) for d in (2, 3)}
    assert ranks == {2: 2, 3: 5}
    print("Ranks of sign-blind two-trial responses (real bit, complex bit):", ranks)
    polynomial, parity = rotation_weights()
    print("Axial generator characteristic polynomial:", polynomial)
    print("Half-turn eigenspaces:", parity)

    # Uniform sign matrices at identical syntax cost, before row/column relabelling.
    tables = list(product((-1, 1), repeat=4))
    frustrated = sum(sp.prod(c) == -1 for c in tables)
    assert frustrated == 8
    print("Nonclassical sign games:", frustrated, "/", len(tables))
    print("Game x two actually realized Bell fields:", frustrated, "/", 2*len(tables))
    print("After conditioning on nonclassical games (R,C):", frustrated, "/", 2*frustrated)
    print("Game x three candidate scalar algebras with tomography filter:",
          frustrated, "/", 3*len(tables))


if __name__ == "__main__":
    main()
