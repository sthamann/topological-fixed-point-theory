"""All relabeling-covariant exchanges preserving the product of two group events."""
from collections import Counter
from functools import reduce
from itertools import combinations, permutations, product

from sympy import Matrix, Rational, ZZ, diag, eye, simplify, sqrt, symbols, zeros
from sympy.matrices.normalforms import smith_normal_form

G = tuple(permutations(range(3)))
identity = G[0]
events = G[1:]
pairs = tuple(product(events, repeat=2))


def multiply(g, h):
    return tuple(g[i] for i in h)


def inverse(g):
    return tuple(g.index(i) for i in range(len(g)))


def conjugate(k, g):
    return multiply(multiply(k, g), inverse(k))


def conjugate_pair(k, pair):
    return tuple(conjugate(k, g) for g in pair)


def exchange_options(representative):
    # A product-preserving permutation is a tuple of fiber permutations.
    # Covariance propagates each representative's permutation to its conjugate fibers.
    # Propagation is well-defined precisely for stabilizer-equivariant permutations.
    fiber = [p for p in pairs if multiply(*p) == representative]
    stabilizer = [k for k in G if conjugate(k, representative) == representative]
    for permuted in permutations(fiber):
        local = dict(zip(fiber, permuted))
        if all(local[conjugate_pair(k, p)] == conjugate_pair(k, local[p])
               for k in stabilizer for p in fiber):
            yield {conjugate_pair(k, p): conjugate_pair(k, local[p]) for k in G for p in fiber}


def braid_relations(exchange):
    """For each basis triple, yield final words and successive pair arguments on both sides."""
    for a,b,c in product(events, repeat=3):
        d,e = exchange[a,b]
        f,g = exchange[e,c]
        h,i = exchange[d,f]
        j,k = exchange[b,c]
        l,m = exchange[a,j]
        n,o = exchange[m,k]
        yield (h,i,g), (l,n,o), ((a,b),(e,c),(d,f)), ((b,c),(a,j),(m,k))


def cycle_lengths(exchange):
    unseen = set(exchange)
    lengths = []
    while unseen:
        point = next(iter(unseen))
        size = 0
        while point in unseen:
            unseen.remove(point)
            point = exchange[point]
            size += 1
        lengths.append(size)
    return sorted(lengths)


def phase_lattice(exchange, relations, orbit_index):
    """Integer rows for the unitary characters satisfying involution and the braid relation."""
    def exponents(positive, negative=()):
        row = [0] * len(set(orbit_index.values()))
        for p in positive:
            row[orbit_index[p]] += 1
        for p in negative:
            row[orbit_index[p]] -= 1
        return tuple(row)

    constraints = {exponents((p, exchange[p])) for p in pairs}
    constraints.update(exponents(left, right) for _,_,left,right in relations)
    relations_matrix = Matrix(sorted(constraints))
    diagonal = smith_normal_form(relations_matrix, domain=ZZ).diagonal()
    return tuple(abs(int(x)) for x in diagonal)


def sign_solutions(exchange, relations, pair_orbits):
    """Full real covariant sign search: one sign per simultaneous conjugacy orbit."""
    for bits in product(range(2),repeat=len(pair_orbits)):
        phase = {p:bit for orbit,bit in zip(pair_orbits,bits) for p in orbit}
        if any(phase[p] != phase[exchange[p]] for p in pairs):
            continue
        if any((sum(phase[p] for p in left)-sum(phase[p] for p in right)) % 2
               for _,_,left,right in relations):
            continue
        yield phase


def alternating_dimensions(exchange, phase, limit):
    """Trace of a k-cycle on V^k gives the Newton power sum for alternating projectors."""
    cycle_traces = []
    dimensions = [1]
    for k in range(1, limit + 1):
        trace = 0
        for initial in product(events, repeat=k):
            word = list(initial)
            sign = 0
            for i in range(k - 1):
                p = (word[i],word[i+1])
                word[i],word[i+1] = exchange[p]
                sign += phase[p]
            if tuple(word) == initial:
                trace += (-1)**sign
        cycle_traces.append(trace)
        numerator = sum((-1)**(j-1)*cycle_traces[j-1]*dimensions[k-j] for j in range(1,k+1))
        assert numerator % k == 0
        dimensions.append(numerator // k)
    return tuple(dimensions)


def audit_recorded_exterior():
    """Terminal ordered-product memory and the oriented-set quotient.

    This is an algebraic quotient for creation operators. A coherent sum over
    orthogonal memory labels is not the physical channel that traces them out.
    """
    attainable = {}
    ordinary = []
    retained = []
    capacities = []
    for degree in range(len(events)+1):
        layer = list(combinations(events,degree))
        for subset in layer:
            attainable[subset] = {
                reduce(multiply, word, identity) for word in permutations(subset)
            }
        ordinary.append(len(layer))
        retained.append(sum(len(attainable[s]) for s in layer))
        capacities.append(max(len(attainable[s]) for s in layer))
    for subset, outputs in attainable.items():
        for g in events:
            if g in subset:
                continue
            extended = tuple(sorted((*subset,g)))
            # The orientation sign for both the lifted map and the quotient is
            # (-1)^#{s in subset | s > g}, independent of the stored product.
            assert all(multiply(p,g) in attainable[extended] for p in outputs)
    assert ordinary == [1,5,10,10,5,1] and retained == [1,5,19,30,15,3]
    assert all(len(outputs) == 3 for subset,outputs in attainable.items() if len(subset)>=3)
    print("basis dimensions of (occupied set, ordered product), and oriented quotient:",retained,ordinary)
    print("worst-case distinguishable memory states by event number:",capacities)
    # F|S,p> = |S> implies FF*|S> = number_of_products(S) |S>.
    # The largest uniform success probability on a basis history for K=cF is 1/max r(S).
    print("optimal squared scale of the coherent erasure Kraus operator:", "1/"+str(max(capacities)))


def native_channel_control():
    """Compare the source channel with the class-averaged group operation probabilities.

    Formulas are transcribed from the Atlas-registered source
    experiments/theory-contracts/landkarte-hypothesen-20261119/check_h09_scharf.py,
    lines 1 through 17. This is a downstream discrimination test, not a primitive.
    """
    q = symbols("q",real=True)
    X = Matrix([[0,1],[1,0]])
    Z = diag(1,-1)
    R = Matrix([[1,-sqrt(3)],[sqrt(3),1]])/2
    U = (eye(2),R,R.T,Z,(-Z+sqrt(3)*X)/2,(-Z-sqrt(3)*X)/2)
    a = Rational(1,18)-q/6
    b = Rational(2,9)-q/6
    weights = (Rational(1,2)+q/6, q/6, q/6, a,b,b)
    mean = simplify((a+2*b)/3)
    central = (*weights[:3],mean,mean,mean)

    def channel(w,op):
        return simplify(sum((p*u*op*u.T for p,u in zip(w,U)),zeros(2)))
    assert channel(weights,X)==Rational(2,3)*X and channel(weights,Z)==Z/3
    assert channel(central,X)==X/2 and channel(central,Z)==Z/2
    assert simplify(b-a)==Rational(1,6)
    variation = simplify((mean-a + 2*(b-mean))/2)
    assert variation==Rational(1,9)
    print("native and class-averaged X,Z factors:",("2/3","1/3"),("1/2","1/2"))
    print("native distance from class-averaged label distribution:",variation)


def run():
    involution = next(g for g in events if multiply(g,g)==identity)
    rotation = next(g for g in events if multiply(g,g)!=identity)
    options = [list(exchange_options(r)) for r in (identity,involution,rotation)]
    orbits = sorted({frozenset(conjugate_pair(k,p) for k in G) for p in pairs}, key=sorted)
    bare_ranks = Counter()
    braid_ranks = Counter()
    signed_braid_ranks = Counter()
    rank_ten_degrees = Counter()
    real_components_matching_square = 0
    component_shapes = Counter()
    orbit_index = {p:i for i,orbit in enumerate(orbits) for p in orbit}
    nondegenerate = 0
    for fiber_maps in product(*options):
        B = {p:q for f in fiber_maps for p,q in f.items()}
        assert len(B)==len(pairs)==len(set(B.values()))
        if any(B[B[p]]!=p for p in pairs):
            continue
        if all(len({B[a,g][0] for g in events})==len(events)
               and len({B[g,a][1] for g in events})==len(events) for a in events):
            nondegenerate += 1
        rank = len(pairs)-len(cycle_lengths(B))
        bare_ranks[rank]+=1
        relations=list(braid_relations(B))
        if all(left==right for left,right,_,_ in relations):
            braid_ranks[rank]+=1
            diag = phase_lattice(B, relations, orbit_index)
            component_shapes[diag] += 1
            # The torsion has exponent at most 2. Each torus component therefore
            # contains exactly 2^free_rank real sign points in unimodular coordinates.
            assert all(d in (0,1,2) for d in diag)
            signs_in_component = 2 ** diag.count(0)
            local_rank_ten = 0
            for phase in sign_solutions(B, relations, orbits):
                signed_rank = rank + sum(phase[p] for p in pairs if B[p] == p)
                signed_braid_ranks[signed_rank] += 1
                if signed_rank == 10:
                    local_rank_ten += 1
                    rank_ten_degrees[alternating_dimensions(B,phase,6)] += 1
            assert local_rank_ten % signs_in_component == 0
            real_components_matching_square += local_rank_ten // signs_in_component

    hurwitz={(g,h):(conjugate(g,h),g) for g,h in pairs}
    swap={(g,h):(h,g) for g,h in pairs}
    assert all(multiply(*p)==multiply(*q) for p,q in hurwitz.items())
    print("choices for one representative per product class:", [len(x) for x in options])
    print("involutive exchanges by rank(I-B):", dict(sorted(bare_ranks.items())))
    print("involutive Yang-Baxter exchanges by rank(I-B):", dict(sorted(braid_ranks.items())))
    print("nondegenerate involutive exchanges:",nondegenerate)
    print("signed involutive Yang-Baxter exchanges by rank(I-B):",dict(sorted(signed_braid_ranks.items())))
    print("Smith diagonals for all U(1) phase extensions:", dict(component_shapes))
    total_components = sum(count*2**diag.count(2) for diag,count in component_shapes.items())
    print("all U(1) components / rank-ten components:", total_components, real_components_matching_square)
    print("alternating dimensions 0..6 for real rank-ten representatives:",dict(rank_ten_degrees))
    assert (total_components,real_components_matching_square)==(26,2)
    assert rank_ten_degrees=={(1,5,10,20,40,80,160):4}
    print("Hurwitz cycle lengths:",cycle_lengths(hurwitz))
    print("ranks for Hurwitz and ordinary flip:",len(pairs)-len(cycle_lengths(hurwitz)),len(pairs)-len(cycle_lengths(swap)))
    print("ordinary flip, products changed:",sum(multiply(*p)!=multiply(*q) for p,q in swap.items()),"/",len(pairs))
    assert sum(bare_ranks.values())==12 and 10 not in bare_ranks
    assert sum(braid_ranks.values())==5 and nondegenerate==0
    audit_recorded_exterior()
    native_channel_control()


if __name__=="__main__":
    run()
