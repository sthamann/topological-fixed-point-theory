"""Exact controls for exchange records, meromorphic modes and the hidden torus history."""
from collections import Counter
from itertools import permutations, product
from math import factorial
import sympy as sp

A = sp.Matrix([[1, 1], [0, 1]])
B = sp.Matrix([[1, 0], [-1, 1]])
I = sp.eye(2)


def compose(p, q):
    return tuple(p[i] for i in q)


def inverse(p):
    return tuple(p.index(i) for i in range(len(p)))


def mod(M, m):
    return M.applyfunc(lambda x: x % m)


def order_mod(M, m, limit):
    power = I
    for order in range(1, limit + 1):
        power = mod(power * M, m)
        if power == I:
            return order
    raise AssertionError("Order exceeds the certified search bound")


def commutator(X, Y):
    return X * Y * X.inv() * Y.inv()


def exchange_and_modes():
    outcomes = list(permutations(range(3)))
    identity = outcomes[0]
    nontrivial = outcomes[1:]
    parity = {p: sum(p[i] > p[j] for i in range(3) for j in range(i + 1, 3)) % 2
              for p in outcomes}
    sizes = Counter(parity[p] for p in nontrivial)
    assert sizes == {1: 3, 0: 2}
    assert [factorial(n) - 1 for n in range(1, 6)] == [0, 1, 5, 23, 119]

    def conjugation_character(p):
        return sum(compose(compose(p, g), inverse(p)) == g for g in nontrivial)

    def meromorphic_character(p):
        return 2 + sum(i == j for i, j in enumerate(p))

    representatives = (identity, (1, 0, 2), (1, 2, 0))
    characters = [[f(p) for p in representatives]
                  for f in (conjugation_character, meromorphic_character)]
    assert characters == [[5, 1, 2], [5, 3, 2]]
    # Multiply the five rational functions by Q; a nonzero polynomial
    # determinant certifies independence on the entire distinct-root locus.
    z, a, b, c = sp.symbols("z a b c")
    roots = (a, b, c)
    Q = sp.prod(z - p for p in roots)
    numerators = [Q, (z - sum(roots) / 3) * Q]
    numerators += [sp.prod(z - p for p in roots if p != root) for root in roots]
    coeffs = sp.Matrix([[sp.expand(f).coeff(z, degree) for f in numerators]
                        for degree in range(5)])
    vandermonde = (a-b)*(a-c)*(b-c)
    assert sp.factor(coeffs.det() / vandermonde) in (1, -1)

    def j(t):
        return sp.cancel(256 * (1-t+t*t)**3 / (t*t*(1-t)**2))

    js = [j(sp.Integer(-1)), j(sp.Integer(3)), sp.simplify(j((1+sp.sqrt(-3))/2))]
    assert js == [1728, sp.Rational(21952, 9), 0]
    print("Nonidentity endpoints, odd/even:", sizes[1], sizes[0])
    print("Characters (identity, swap, 3-cycle), endpoints vs meromorphic:", characters)
    print("Five pole-mode determinant / Vandermonde:", sp.factor(coeffs.det()/vandermonde))
    print("Elliptic j at t=-1, 3, (1+i sqrt(3))/2:", js)


def torus_probe():
    assert A * B * A == B * A * B
    C = commutator(A**2, B**2)
    assert C == sp.Matrix([[13, 8], [8, 5]])
    u = sp.symbols("u")
    assert C.charpoly(u).as_expr() == u**2 - 18*u + 1
    assert mod(C, 2) == I and mod(C, 4) == I
    assert mod(C**2, 3) == mod(-I, 3)
    orders = {m: order_mod(C, m, m**4) for m in range(2, 8)}
    # Compute fixed-product local comparison with exactly the same word grammar.
    letters = [A**2, A**-2, B**2, B**-2]
    control = Counter(order_mod(commutator(x, y), 3, 24)
                      for x, y in product(letters, repeat=2))
    group = [sp.Matrix([[a, b], [c, d]]) for a, b, c, d in product(range(3), repeat=4)
             if (a*d - b*c) % 3 == 1]
    assert control == {1: 8, 4: 8}
    assert len(group) == 24 and sum(order_mod(g, 3, 24) == 4 for g in group) == 6
    # A and B have faithful independent transpositions on the three nonzero mod-2 vectors.
    points = [sp.Matrix(v) for v in [(1, 0), (0, 1), (1, 1)]]
    permutations_mod2 = [tuple(points.index(mod(g * v, 2)) for v in points) for g in (A, B)]
    s, t = permutations_mod2
    assert compose(s, s) == compose(t, t) == (0, 1, 2) and compose(s, t) != compose(t, s)
    print("Hidden commutator: endpoint=id, winding=0, torus=", C.tolist(), "trace=", C.trace())
    print("Orders of hidden commutator on torsion:", orders)
    print("Identical-cost pure exchange commutators on 3-torsion:", dict(control))
    print("Order four in SL(2,F3): 6/24")


def matrix_key(M, m):
    return tuple(int(x % m) for x in M)


def subgroup(generators, m):
    seen = {matrix_key(I, m)}
    todo = [I]
    while todo:
        g = todo.pop()
        for h in generators:
            gh = mod(g * h, m)
            key = matrix_key(gh, m)
            if key not in seen:
                seen.add(key)
                todo.append(gh)
    return seen


def native_lines_from_braid():
    # 0,1,w,w+1 represent F2[w]/(w²+w+1); only the two order-three
    # endpoint operations need the quadratic splitting field.
    from hermitian_origin import times, bar
    rays = [(1,0),(0,1),(1,1),(2,1),(3,1)]
    def line(x,y):
        return (times(x,bar(y)),1) if y else (1,0)
    def act(M, v):
        a,b,c,d = M
        x,y = v
        return line(times(a,x)^times(b,y),times(c,x)^times(d,y))
    def h(v, w):
        x,y = v
        z,t = w
        return times(bar(x),t)^times(bar(y),z)
    operators = subgroup([A,B],2)
    endpoints = sorted(operators - {matrix_key(I,2)})
    fixed = {M: {v for v in rays if act(M,v)==v} for M in endpoints}
    rational = {v for v in rays if v==tuple(bar(x) for x in v)}
    single = {next(iter(v)) for v in fixed.values() if len(v)==1}
    pair = {v for lines in fixed.values() if len(lines)==2 for v in lines}
    assert single == rational == {v for v in rays if not h(v,v)}
    assert pair == set(rays)-rational and len(pair)==2
    assert Counter(map(len,fixed.values()))=={1:3,2:2}

    # Symplectic (hyperbolic) form from E[2] transformed to the diagonal
    # Hermitian h used in the existing five-ray audit. Columns form an exact
    # orthonormal basis; both choices exchanged by Frobenius are valid.
    f=(2,1)
    g=(3,1)
    assert (h(f,f), h(g,g), h(f,g))==(1,1,0)
    def conjugate(M,N):
        # GL(2,F2)=SL(2,F2); adjugate reduces to inverse.
        m,n=sp.Matrix(2,2,M),sp.Matrix(2,2,N)
        return matrix_key(m*n*m.adjugate(),2)
    matched = sum(all(bijection[conjugate(M,N)]==act(M,bijection[N])
                      for M in operators for N in endpoints)
                  for images in permutations(rays)
                  for bijection in [dict(zip(endpoints,images))])
    assert matched==2
    print("Braid E[2] eigendirections over minimal splitting field F4:",len(single),"+",len(pair))
    print("Exact h=Weil pairing after Frobenius; endpoint-equivariant bijections:",matched,"/ 120")
    # The normalized winding fiber x*y*(x+y)=1 has another elliptic
    # completion. Check its cubic equation modulo the fiber relation.
    x,y=sp.symbols("x y")
    X,Y=1/(x+y),(x-y)/(x+y)
    assert sp.cancel((Y**2-1+4*X**3)*(x+y)**3 + 4*(x*y*(x+y)-1)) == 0
    print("Ordered Vandermonde fiber: Y²=1−4X³, j=0")


def shortest_invisible():
    letters = {1: A, -1: A.inv(), 2: B, -2: B.inv()}
    words = [((), I)]
    hits = []
    for length in range(1, 5):
        words = [(prefix + (i,), M * matrix) for prefix, M in words
                 for i, matrix in letters.items() if not prefix or i != -prefix[-1]]
        hits = [M for word, M in words if sum(1 if i>0 else -1 for i in word) == 0
                and mod(M, 2) == I and M != I]
        assert len(hits) == (8 if length == 4 else 0)
    assert len(words) == 108 and {M.trace() for M in hits} == {6}
    assert {order_mod(M, 3, 24) for M in hits} == {4}
    group3 = subgroup([A, B], 3)
    assert len(group3) == 24
    matrices = [sp.Matrix(2, 2, key) for key in group3]
    # All determinants are 1 mod 3; use the adjugate for inversion in F3.
    derived = subgroup([mod(x*y*x.adjugate()*y.adjugate(), 3)
                        for x in matrices for y in matrices], 3)
    detected = subgroup(hits, 3)
    assert detected == derived and len(derived) == 8
    orders = Counter(order_mod(sp.Matrix(2, 2, key), 3, 24) for key in detected)
    assert orders == {1: 1, 2: 1, 4: 6}
    assert len({matrix_key(M, 3) for M in hits}) == 6
    print("Shortest winding-zero, endpoint-neutral torus witnesses: 8/108, length 4, trace 6")
    print("All hidden words on 3-torsion: Q8; element orders", dict(sorted(orders.items())))


def finite_reflection_witness():
    # Normalized character of the actual permutation on all nine E[3] points.
    vectors = [sp.Matrix(v) for v in product(range(3), repeat=2)]
    def phi(M):
        return sp.Rational(sum(mod(M*v, 3) == v for v in vectors), len(vectors))
    G = sp.Matrix([[phi(A**(i+j)) for j in (1, 2)] for i in (1, 2)])
    assert G == sp.Matrix([[sp.Rational(1, 3), 1], [1, sp.Rational(1, 3)]])
    assert G.eigenvals() == {sp.Rational(-2,3): 1, sp.Rational(4,3): 1}
    print("Torsion permutation trace, reflected Gram on a,a^2:", G.tolist(),
          "eigenvalues=", sorted(G.eigenvals()))


def integral_lorentz_probe():
    # The smooth completion of the ordered winding fiber has j=0. Its
    # order-three automorphism gives the maximal Eisenstein order R.
    # Hermitian forms on R tensor the rank-two intersection lattice have
    # determinant of real signature (1,3). Which are invariant under BOTH
    # standard monodromies? Hermitian form matrices transform by congruence.
    a,d,x,y = sp.symbols("a d x y",real=True)
    H = sp.Matrix([[a,x+sp.I*y],[x-sp.I*y,d]])
    equations = list(A.T*H*A - H) + list(B.T*H*B - H)
    invariant = sp.linsolve(equations,[a,d,x,y])
    assert invariant==sp.FiniteSet((0,0,0,y))
    assert H.det().subs({a:0,d:0,x:0})==-y**2
    primitive = sp.Matrix([[0,sp.sqrt(-3)],[-sp.sqrt(-3),0]])
    assert primitive.det()==-3 and A.T*primitive*A==B.T*primitive*B==primitive
    print("Fixed real Hermitian forms (a,d,Re z,Im z):", invariant,
          "primitive Eisenstein determinant =", primitive.det())
    print("Invariant positive semidefinite forms: zero only")
    u=sp.symbols("u",real=True)
    v=sp.symbols("v",positive=True)
    tau=u+sp.I*v
    def hodge(t):
        t=sp.expand_complex(t)
        x,y=sp.re(t),sp.im(t)
        return sp.Matrix([[1,x],[x,x*x+y*y]])/y
    g=hodge(tau)
    assert sp.simplify(g.det())==1
    for M in (A,B):
        a,b,c,d=list(M)
        new=(b+d*tau)/(a+c*tau)
        assert (hodge(new)-M.T*g*M).applyfunc(sp.simplify)==sp.zeros(2)
    print("Moving positive Hodge metric: determinant 1, covariant under both braid generators")


def reflection_probe():
    # Normalized positive-definite full-group functions and their OS Gram
    # matrices on the positive braid monoid, checked exactly on words <=3.
    r_values = [sp.Rational(1, 2), sp.Rational(2, 3), sp.Integer(1)]
    e = (0, 1, 2)
    letters = ((1, 0, 2), (0, 2, 1))
    records = [(e, 0)]
    level = [(e, 0)]
    for _ in range(3):
        level = [(compose(g, s), n + 1) for (g, n), s in product(level, letters)]
        records.extend(level)
    for r in r_values:
        endpoints = list(permutations(range(3)))
        embedding = sp.Matrix([[int(g == h)*r**n for h in endpoints] for g, n in records])
        reflected = sp.Matrix([[int(g == h)*r**(n+k) for h,k in records] for g,n in records])
        assert reflected == embedding * embedding.T
        assert reflected.rank() == 6
    print("Reflection Gram on positive words <=3: rank 6 at r=1/2, 2/3, 1")
    print("Same double exchange response:", [r*r for r in r_values])


if __name__ == "__main__":
    exchange_and_modes()
    torus_probe()
    native_lines_from_braid()
    shortest_invisible()
    integral_lorentz_probe()
    finite_reflection_witness()
    reflection_probe()
