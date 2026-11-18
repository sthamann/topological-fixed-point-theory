"""Independent exact replays of the accompanying F1 research report."""

from itertools import combinations, permutations, product

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

from reconstruction import require


def self_dual_relations():
    require(all((1<<rank)-1 < 2*rank for rank in (1,2)),
            "two times the rank exceeds available distinct nonzero binary columns")
    nonzero = [sp.Matrix(column) for column in product((0,1), repeat=3) if any(column)]
    candidates = list(combinations(range(len(nonzero)),6))
    for indices in candidates:
        columns = [nonzero[i] for i in indices]
        require(sum(int(sp.Matrix.hstack(*basis).det()) % 2 for basis in combinations(columns,3)) == 16,
                "wrong rank-three basis census")

    edges = tuple(combinations(range(4),2))
    incidence = sp.zeros(4,6)
    for column,(i,j) in enumerate(edges):
        incidence[i,column], incidence[j,column] = 1,-1
    bases = {frozenset(indices) for indices in combinations(range(6),3)
             if int(incidence[:3, list(indices)].det()) % 2}
    opposite = [edges.index(tuple(sorted(set(range(4))-set(edge)))) for edge in edges]
    require({frozenset(opposite[i] for i in basis) for basis in bases}
            == {frozenset(set(range(6))-basis) for basis in bases}, "K4 dual map")
    laplacian = incidence*incidence.T
    covariance_projector = sp.eye(4)-sp.ones(4)/4
    require(laplacian == 4*covariance_projector, "K4 Laplacian")
    smith = smith_normal_form(laplacian[:3,:3],domain=sp.ZZ).diagonal()
    require([abs(x) for x in smith] == [1,4,4],"integral discriminant is Z4²")
    reflected_gram = (covariance_projector/4).extract([0,1],[2,3])
    require(reflected_gram.eigenvals() == {-sp.Rational(1,8):1,0:1}, "scalar reflection must fail")
    refinements = [sum((a*c+b*d+t*(a+c))%2 == 0 for a,b,c,d in product((0,1),repeat=4)) for t in (0,1)]
    require(refinements == [10,6], "even and odd refinements are distinct")
    return {"minimum_nonempty_simple_binary_self_dual_size":6,"rank":3,"rank_three_presentations":len(candidates),
            "bases":len(bases),"unique_up_to_binary_linear_relabelling":"GL(3,2) transitive on omitted nonzero column",
            "smith":[1,4,4],"critical_group":"Z4 × Z4", "even_and_odd_null_counts_including_zero":refinements,
            "scalar_reflection_positive":False,"odd_reflection_quotient_rank":1,
            "scope":"minimal simple binary self-dual relations; no native code, scalar normalizer or time is inferred"}


def orientation_flags():
    def conjugate(move, cycle):
        return tuple(move[cycle[move.index(i)]] for i in range(len(move)))

    def inverse(cycle):
        return tuple(cycle.index(i) for i in range(len(cycle)))

    def cyclic_steps(cycle,start):
        a = cycle[start]
        b = cycle[a]
        return a != start and b != start and cycle[cycle[b]] == start

    moves = tuple(permutations(range(5)))
    flags = {(next(i for i in range(5) if p[i] == i),p) for p in moves
             if sum(i==p[i] for i in range(5)) == 1 and all(cyclic_steps(p,i) for i in range(5) if p[i] != i)}
    first,cycle = min(flags)
    orbit = {(move[first],conjugate(move,cycle)) for move in moves}
    require(len(flags) == 30 and orbit == flags and sum((move[first],conjugate(move,cycle)) == (first,cycle)
                for move in moves) == 4,"flags orbit and stabilizer")
    histories = [p for p in permutations(range(4)) if all(cyclic_steps(p,i) for i in range(4))]
    pairs = [p for p in histories if p < inverse(p)]
    representatives = ((0,1,2,3),(1,0,2,3),(1,0,3,2),(1,2,0,3),(1,2,3,0))
    character = [sum(int(conjugate(move,p)==p)-int(conjugate(move,p)==inverse(p)) for p in pairs)
                 for move in representatives]
    require(character == [3,-1,-1,0,1], "oriented history character")
    return {"flags":len(flags),"orbits":1,"stabilizer_order":4,"groupoid_cardinality":"1/4",
            "orientation_free_cardinality":"1/8","anchor_only_cardinality":"1/24",
            "odd_history_character":character,"standard_A3_character":[3,1,-1,0,-1],
            "scope":"the two 3-dimensional S4 representations differ; one orbit of flags does not fix the state"}


def reciprocal_sources():
    t, weight = sp.symbols("t w",real=True)
    q = t**2-5*t/2+sp.Rational(15,16)

    def laplace_image(polynomial):
        return -t*sp.diff(polynomial,t,2)+(2*t-sp.Rational(3,2))*sp.diff(polynomial,t)+(sp.Rational(3,2)-t)*polynomial

    one = sp.S.One
    f_q = laplace_image(laplace_image(one))-5*laplace_image(one)/2+sp.Rational(15,16)
    require(sp.expand(f_q-q)==0 and q.subs(t,sp.Rational(5,4))==-sp.Rational(5,8),"positive Fourier self-reciprocal family")

    def radial_moment(k):
        return sp.rf(sp.Rational(3,2),k)

    def normalized_radial_moment(k):
        raw = radial_moment(k)+weight*(radial_moment(k+2)-5*radial_moment(k+1)/2+15*radial_moment(k)/16)
        return sp.simplify(raw/(1+15*weight/16))

    variance = normalized_radial_moment(1)/(12*sp.pi)
    cumulant = normalized_radial_moment(2)/(80*sp.pi**2)-3*variance**2
    sample = sp.Rational(16,25)
    require(variance.subs(weight,0)==1/(8*sp.pi) and variance.subs(weight,sample)==1/(4*sp.pi)
            and sp.simplify(cumulant.subs(weight,sample))==-3/(320*sp.pi**2),"distinct self-dual variances and fourth cumulant")
    d = tuple(product(range(4),repeat=2))
    f = sp.Matrix([[sp.I**((2*a*c+a*e+b*c+2*b*e)%4)/4 for c,e in d] for a,b in d])
    parity_space = sp.Matrix([int(a%2==b%2==0) for a,b in d])
    vacuum = sp.eye(16)[:,0]+sp.ones(16,1)/4

    def probability(v):
        return sum(v[i]**2 for i in range(len(v)) if parity_space[i])/v.dot(v)

    require(f*f.H==sp.eye(16) and f**4==sp.eye(16) and f**2!=sp.eye(16)
            and f*vacuum==vacuum and f*parity_space==parity_space,"discrete Fourier-fixed sources")
    require([probability(v) for v in (vacuum,vacuum+parity_space)] == [sp.Rational(7,10),sp.Rational(13,14)],
            "different positive Fourier-fixed cap probabilities")
    return {"continuous_normalizations":[str(variance.subs(weight,v)) for v in (0,sample)],
            "second_source_fourth_cumulant":str(sp.simplify(cumulant.subs(weight,sample))),
            "discrete_cap_probabilities":["7/10","13/14"],
            "independent_dyadic_composition":"forces Gaussian for finite covariance; verified alternative fails by its nonzero fourth cumulant",
            "scope":"positive self reciprocity alone leaves both covariance and discrete cap unfixed"}
