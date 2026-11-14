"""Enumerate finite labelled simplicial complexes and rational process seeds."""
from fractions import Fraction
from itertools import combinations, product


def subsets(items):
    for mask in range(1<<len(items)):
        yield tuple(value for i,value in enumerate(items) if mask>>i&1)


def simplicial_complexes(vertices, max_dimension):
    if not 1<=vertices<=8 or not 0<=max_dimension<vertices:
        raise ValueError("invalid finite complex bounds")
    levels = [tuple((v,) for v in range(vertices))]

    def extend(dimension):
        if dimension>max_dimension:
            yield tuple(levels)
            return
        boundary = set(levels[-1])
        faces = tuple(c for c in combinations(range(vertices),dimension+1)
                      if all(face in boundary for face in combinations(c,dimension)))
        for present in subsets(faces):
            levels.append(present)
            yield from extend(dimension+1)
            levels.pop()

    yield from extend(1)


def chain_report(levels):
    import sympy as sp
    from sympy.polys.matrices import DomainMatrix
    ranks = [0]
    for lower,upper in zip(levels,levels[1:]):
        indices = {face:i for i,face in enumerate(lower)}
        boundary = sp.zeros(len(lower),len(upper))
        for col,cell in enumerate(upper):
            for k in range(len(cell)):
                boundary[indices[cell[:k]+cell[k+1:]],col] = (-1)**k
        ranks.append(DomainMatrix.from_Matrix(boundary).convert_to(sp.QQ).rank())
    ranks.append(0)
    cells = list(map(len,levels))
    return {"cells":cells,"betti_Q":[n-ranks[i]-ranks[i+1] for i,n in enumerate(cells)],
            "euler":sum((-1)**i*n for i,n in enumerate(cells))}


def finite_space(vertices,dimension):
    counts = {}
    examples = {}
    tested = 0
    for levels in simplicial_complexes(vertices,dimension):
        result = chain_report(levels)
        key = str(result["betti_Q"])
        tested += 1
        counts[key] = counts.get(key,0)+1
        examples.setdefault(key,{"complex":levels,"invariants":result})
    return {"status":"EXHAUSTED_LABELLED_COMPLEX_CLASS","vertices":vertices,
            "max_dimension":dimension,"count":tested,"homology_counts":counts,
            "representatives_by_homology":examples,
            "scope":"all labelled complexes at these bounds; equal homology is not topological equivalence"}


def rational_seeds(dimension, height):
    """All square matrices over rationals of bounded numerator and denominator."""
    if dimension<1 or height<1:
        raise ValueError("positive bounds required")
    values = sorted({Fraction(p,q) for p in range(-height,height+1) for q in range(1,height+1)})
    for entries in product(values, repeat=dimension**2):
        yield tuple(tuple(str(x) for x in entries[j:j+dimension])
                    for j in range(0,len(entries),dimension))


def finite_process_space(dimension,height):
    import sympy as sp
    from language import positive
    from reconstruction import matrix, equivalent
    if dimension>2 or height>2:
        raise ValueError("configured exact enumeration supports dimension and height at most two")
    tested, accepted = 0, []
    for spec in rational_seeds(dimension,height):
        tested += 1
        t = matrix([list(row) for row in spec])
        if positive(t) and positive(sp.eye(dimension)-t):
            accepted.append((spec,t))
    state = sp.eye(dimension)[:,0]
    groups = []
    # Compare all words over a genuinely marked alphabet, including the
    # boundary reversal. The bare transfer is not enough to compare sources.
    reversal = sp.eye(dimension)[::-1,:]
    for spec,t in accepted:
        model = (state,state.H,{"transfer":t,"reversal":reversal})
        for group in groups:
            proof = equivalent(group["representative"],model)
            if proof["status"]=="ALL_WORDS_EQUAL":
                group["members"].append(spec)
                break
        else:
            groups.append({"representative":model,"members":[spec]})
    return {"status":"EXHAUSTED_BOUNDED_RATIONAL_CLASS","dimension":dimension,"height":height,
            "tested":tested,"positive_contractions":len(accepted),
            "equivalence_classes":[group["members"] for group in groups],
            "alphabet":["transfer","reversal"],"physical_solution":False}


def custom_complex(cells):
    if not isinstance(cells,list) or not cells:
        raise ValueError("nonempty list of simplices required")
    if any(not isinstance(cell,list) or not cell or any(type(x) is not int for x in cell)
           or len(set(cell))!=len(cell) for cell in cells):
        raise ValueError("each simplex requires distinct integer vertices")
    closed = set()
    for cell in cells:
        for size in range(1,len(cell)+1):
            closed.update(combinations(sorted(cell),size))
    levels = [tuple(sorted(c for c in closed if len(c)==size))
              for size in range(1,max(map(len,closed))+1)]
    return {"status":"EXACT_CUSTOM_COMPLEX","invariants":chain_report(levels),
            "closure_added_faces":len(closed)-len({tuple(sorted(c)) for c in cells})}


def coverage():
    return {"builtin_topologies":"configured families of labelled triangulations",
            "generated_topologies":"exhaustive labelled finite complexes within each declared bound",
            "rational_processes":"exhaustive matrices in the configured dimensions and rational height",
            "custom_inputs":"exact coherent Gram and marked operator forms over Gaussian rationals",
            "global_reconstruction":"complete equality decision for supplied finite word representations",
            "physical_reconstruction":"exact finite null quotient, marked maps, transfer and regional checks",
            "open_classes":["unbounded or noncomputable dynamics","general interacting continuum",
                            "arbitrary infinite local nets","unknown source and observables"],
            "universal_exhaustiveness_claimed":False}
