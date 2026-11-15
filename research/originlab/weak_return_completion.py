"""Joint neutral return spaces for the D8 and two half-spinor completions."""

from itertools import combinations, product

import sympy as sp
from sympy.polys.matrices import DomainMatrix

from affine_operators import LatticeFock, dot, integral, lattice_frame, neg, return_matrix
from language import positive
from reconstruction import require
from weak_source import triangle_components


BRANCHES = ("d8", "weak_spinor_minus", "weak_spinor_plus")


def late_orbits():
    tensors = triangle_components()
    return {name: tuple(neg(word[0]) for word in tensors[name]) for name in BRANCHES}


def _return_gram(roots):
    # Positive highest-to-neutral vectors use two affine lowering modes.
    # In degree two their Heisenberg polynomial is a(-1)^2 + a(-2),
    # whose normalized overlap at a·b=d is (d^2+d)/6.
    def overlap(a,b):
        d = dot(a,b)
        return sp.Rational(d*d+d,6)

    return sp.Matrix([[overlap(a,b) for b in roots] for a in roots])


def _row(operator):
    return list(operator.reshape(1,operator.rows*operator.cols))


def _words(gram):
    """Shortest representative for each ordered endpoint, up to length three."""
    size = gram.rows
    pivots = gram.rref()[1]
    metric = gram.extract(pivots,pivots)
    bras = gram.extract(range(size),pivots)
    kets = metric.inv()*bras.T
    adjacent = [{j for j in range(size) if gram[i,j]} for i in range(size)]
    words = {(i,i): (i,) for i in range(size)}
    layers = []
    for length in (1,2,3):
        if length > 1:
            for word in tuple(words.values()):
                for j in sorted(adjacent[word[-1]]):
                    words.setdefault((word[0],j),(*word,j))
        matrix = sp.Matrix([_row(kets[:,i]*bras[j,:]) for i,j in words])
        layers.append((len(words),DomainMatrix.from_Matrix(matrix).rank()))
    rows = DomainMatrix.from_Matrix(matrix.T).rref()[1]
    ordered = list(words.values())
    basis = [ordered[i] for i in rows]
    holdout = next((*word,j) for word in ordered for j in sorted(adjacent[word[-1]])
                   if (word[0],j) not in words)
    return metric, pivots, kets, bras, basis, layers, holdout


def _operator(word, gram, kets, bras):
    coefficient = sp.prod(gram[i,j] for i,j in zip(word,word[1:]))
    return coefficient*kets[:,word[0]]*bras[word[-1],:]


def _grade_excess(fock, frame, orbits, degree):
    # The affine root sl2 has weights -1, +1 on the return. In grades two
    # and three there is only one nontrivial return eigenvalue (4).
    # The entire exceptional part comes from its weight-three space,
    # the charge 2*a with oscillator degree d-2. Map it down twice.
    # The next possible high weight is at d=6.
    keys = [next(iter(vector)) for vector in fock.basis((0,)*8,degree)]
    vectors, columns = {}, {}
    for name, roots in orbits.items():
        vectors[name] = []
        for root in roots:
            a = integral(frame,root)
            high_weight = tuple(2*x for x in a)
            vectors[name].append([
                fock.vertex(neg(a),1,fock.vertex(neg(a),1,vector))
                for vector in fock.basis(high_weight,degree-2)])
        columns[name] = [[vector.get(key,0) for key in keys]
                         for line in vectors[name] for vector in line]

    def rank(names):
        return DomainMatrix.from_Matrix(sp.Matrix([v for name in names for v in columns[name]])).rank()

    branches = {name: rank([name]) for name in orbits}
    pairs = [{"branches":list(pair),"span_rank":rank(pair)} for pair in combinations(orbits,2)]
    return vectors, {"neutral_grade":degree,"neutral_dimension":len(keys),
                     "separate_span_ranks":branches,"pairwise_span_ranks":pairs,
                     "combined_span_rank":rank(orbits),
                     "single_root_excess_rank":len(next(iter(vectors.values()))[0])}


def controls():
    orbits = late_orbits()
    roots = sum(orbits.values(),())
    gram = _return_gram(roots)
    frame,_ = lattice_frame()
    fock = LatticeFock([[dot(a,b) for b in frame] for a in frame])
    vectors, grade_two = _grade_excess(fock,frame,orbits,2)
    higher_vectors, grade_three = _grade_excess(fock,frame,orbits,3)
    flat = [line[0] for orbit in vectors.values() for line in orbit]
    require(gram*12 == sp.Matrix([[fock.inner(a,b) for b in flat] for a in flat]),
            "affine modes disagree with the neutral grade-two Gram matrix")
    # One explicit matrix per orbit also verifies the positive return convention.
    keys = fock.basis((0,)*8,2)
    for name, orbit in orbits.items():
        vector = vectors[name][0][0]
        ket = sp.Matrix([vector.get(next(iter(key)),0) for key in keys])
        bra = sp.Matrix([[fock.inner(vector,key) for key in keys]])
        require(sp.Matrix(return_matrix(fock,frame,orbit[0],2))-sp.eye(len(keys))==ket*bra/4,
                "root return did not agree with the weight-three construction")
        root = integral(frame,orbit[0])
        for excited in higher_vectors[name][0]:
            composed = fock.vertex(neg(root),1,fock.vertex(root,-1,excited))
            require(not fock.combine((-1,composed),(-4,excited)),
                    "grade-three weight-three descendants must have return eigenvalue four")
    metric,_,_,_,basis,layers,holdout = _words(gram)
    require(grade_two["combined_span_rank"]==metric.rows==14 and grade_three["combined_span_rank"]==98 and
            list(grade_two["separate_span_ranks"].values())==[6]*3 and
            list(grade_three["separate_span_ranks"].values())==[47]*3 and
            [row["span_rank"] for row in grade_two["pairwise_span_ranks"]]==[10]*3 and
            [row["span_rank"] for row in grade_three["pairwise_span_ranks"]]==[74]*3 and
            layers==[(18,18),(108,108),(252,196)], "joint excess overlap or ordered span changed")
    return {"status":"EXACT_AFFINE_JOINT_WEAK_RETURN_SPANS",
            "roots_doubled":[list(root) for root in roots],"branch_order":list(orbits),
            "support_by_grade":[grade_two,grade_three],
            "ordered_operator_span_by_maximum_length":[rank for _,rank in layers],
            "continuous_unitary_normalizer_of_the_finite_marked_alphabet_is_scalar_on_grade_two_support":True,
            "word_basis":[list(word) for word in basis],
            "example_length_four_holdout":list(holdout),
            "scope":"resolved charged root legs in the fixed affine frame; physical dressing, time and regional transport remain open",
            "physical_source_identity_proved":False}


def checked_grade_two(data, contract):
    """Unique positive functional on the declared grade-two excess algebra."""
    roots = sum(late_orbits().values(),())
    require(data.get("roots_doubled")==[list(root) for root in roots],"use the pinned order of the 18 late roots")
    gram = _return_gram(roots)
    metric,pivots,kets,bras,basis,_,_ = _words(gram)
    supplied = data["basis_words"]
    require([row["word"] for row in supplied]==[list(word) for word in basis],
            "supply the full emitted grade-two word basis in order")
    operators = [_operator(word,gram,kets,bras) for word in basis]
    holdout = data["held_out_words"]
    words = [row["word"] for row in holdout]
    require(words and all(isinstance(word,list) and 4<=len(word)<=8 and
                          all(type(i) is int and 0<=i<gram.rows for i in word) for word in words) and
            len({tuple(word) for word in words})==len(words), "supply distinct held-out words of length four to eight")
    withheld = [_operator(word,gram,kets,bras) for word in words]
    require(any(operator and not any(gram[word[0],j] and gram[j,word[-1]] for j in range(gram.rows))
                for word,operator in zip(words,withheld)),
            "hold back a nonzero word whose endpoints require three or more overlap steps")
    # Tr(D A) = vec(D)^T vec(A^T). D represents the reconstructed
    # functional in the non-orthogonal frame and is positive iff G D is PSD.
    design = sp.Matrix([_row(operator.T) for operator in operators])
    values = sp.Matrix([contract(row["word_form"]) for row in supplied])
    density = design.inv(method="DM")*values
    density = density.reshape(metric.rows,metric.cols)
    observed = [contract(row["word_form"]) for row in holdout]
    predicted = [(density*operator).trace() for operator in withheld]
    residuals = [sp.simplify(a-b) for a,b in zip(observed,predicted)]
    form = (metric*density).applyfunc(sp.simplify)
    compatible = positive(form)
    nonzero = density.trace()!=0
    if not compatible or any(residuals):
        status = "REFUTED_SUPPLIED_GRADE_TWO_EXCESS_MODEL"
    elif nonzero:
        status = "RECONSTRUCTED_POSITIVE_GRADE_TWO_EXCESS_FUNCTIONAL"
    else:
        status = "INCONCLUSIVE_ZERO_EXCESS_FUNCTIONAL"
    return {"status":status,"operator_span_dimension":len(basis),
            "positive_functional":compatible,"absolute_excess_mass":str(density.trace()),
            "model_basis_root_indices":list(pivots),
            "model_gram":[list(map(str,row)) for row in metric.tolist()],
            "linear_candidate_form":[list(map(str,row)) for row in form.tolist()],
            "reference": "the supplied common absolute endpoint reference",
            "held_out":[{"word":word,"absolute":str(a),"predicted":str(b),"residual":str(r)}
                        for word,a,b,r in zip(words,observed,predicted,residuals)],
            "physical_source_identity_proved":False}
