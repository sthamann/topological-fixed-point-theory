"""Exact test of the six charged D8 return directions on neutral affine grade two."""
from itertools import combinations, product

import sympy as sp

from affine_operators import ROOT_SET, LatticeFock, charge_signature, dot, integral, lattice_frame, neg, return_matrix


def controls():
    # The late D8 weak leg is the triplet-doublet orbit e_i + e_j.
    roots = [tuple(2 if k in (i,j) else 0 for k in range(8))
             for i in range(3) for j in (3,4)]
    assert all(r in ROOT_SET and charge_signature(r)[1:]==(1,sp.Rational(1,6)) for r in roots)
    frame, _ = lattice_frame()
    fock = LatticeFock([[dot(a,b) for b in frame] for a in frame])
    vectors = []
    for root in roots:
        a = integral(frame,root)
        lower = neg(a)
        highest = fock.vacuum(tuple(2*x for x in a))
        vectors.append(fock.vertex(lower,1,fock.vertex(lower,1,highest)))
    raw = sp.Matrix([[fock.inner(left,right) for right in vectors] for left in vectors])
    # The overlapping roots form the triangular prism. Its adjacency
    # eigenvalues are 3, 1, 0, 0, -2, -2.
    expected = sp.Matrix(6,6,lambda i,j: 1 if i==j else sp.Rational(dot(roots[i],roots[j]),3))
    assert raw == 12*expected
    assert expected.det()==sp.Rational(8,27)
    assert expected.eigenvals()=={2:1,sp.Rational(4,3):1,1:2,sp.Rational(1,3):2}
    for grade in (0,1):
        assert all(sp.Matrix(return_matrix(fock,frame,r,grade))==sp.eye(1 if grade==0 else 8)
                   for r in roots)

    # The actual excess E*E-I acts on each highest-derived line as a
    # rank-one projector times three, in the original oscillator metric.
    basis = fock.basis((0,)*8,2)
    matrices = [sp.Matrix(return_matrix(fock,frame,r,2))-sp.eye(len(basis)) for r in roots]
    for vector, excess in zip(vectors,matrices):
        ket = sp.Matrix([vector.get(next(iter(element)),0) for element in basis])
        bra = sp.Matrix([[fock.inner(vector,element) for element in basis]])
        assert excess==ket*bra/4
    assert sum(matrices,sp.zeros(len(basis))).rank()==6
    noncommuting = [(i,j) for i,j in combinations(range(6),2) if matrices[i]*matrices[j]!=matrices[j]*matrices[i]]
    assert noncommuting==[(i,j) for i,j in combinations(range(6),2) if dot(roots[i],roots[j])==1]
    assert len(noncommuting)==9
    # Each rail is a triangle. An independent commuting frame can choose
    # at most one root per triangle, and there are nonadjacent cross pairs.
    assert any(dot(roots[i],roots[j])==0 for i,j in combinations(range(6),2))
    assert all(not all(dot(roots[i],roots[j])==0 for i,j in combinations(group,2))
               for group in combinations(range(6),3))

    # A word in the projectors is a multiple of the outer product
    # of its endpoint vectors. It is nonzero if each overlap is
    # nonzero. Invertibility of the Gram makes these outer products
    # independent for distinct ordered endpoints.
    endpoints = set()
    word_spans = []
    for length in (1,2,3):
        for word in product(range(6),repeat=length):
            if all(raw[i,j]!=0 for i,j in zip(word,word[1:])):
                endpoints.add((word[0],word[-1]))
        word_spans.append(len(endpoints))
    assert word_spans==[6,24,36]

    # The two positive densities differ in a coherence between a
    # nonadjacent pair. In the normalized vector basis the dual vectors
    # are the columns of G^-1. Their imaginary antisymmetric combination
    # has trace zero and vanishes on every one- or two-projector word.
    i,j = next((i,j) for i,j in combinations(range(6),2) if raw[i,j]==0)
    intermediate = next(k for k in range(6) if raw[i,k] and raw[k,j])
    units = sp.eye(6)
    dual = expected.inv()
    difference = sp.I*(dual[:,i]*units[j,:]-dual[:,j]*units[i,:])
    assert difference.trace()==0 and (expected*difference).is_hermitian
    step = 1/(6*(1+(difference*difference).trace()))
    densities = [units/6+sign*step*difference for sign in (-1,1)]
    for density in densities:
        form = expected*density
        assert form.is_hermitian and all(form[:k,:k].det()>0 for k in range(1,7))
    projectors = [units[:,k]*expected[k,:] for k in range(6)]
    assert all((difference*p).trace()==0 for p in projectors)
    assert all((difference*p*q).trace()==0 for p,q in product(projectors,repeat=2))
    witness = projectors[i]*projectors[intermediate]*projectors[j]
    contrast = ((densities[1]-densities[0])*witness).trace()
    assert contrast!=0 and sp.re(contrast)==0
    return {
        "status":"EXACT_AFFINE_SIX_WEAK_RETURNS",
        "normalized_gram":[[str(v) for v in row] for row in expected.tolist()],
        "roots_doubled":[list(r) for r in roots],
        "gram_determinant":"8/27", "excess_ranks_by_neutral_grade":[0,0,6],
        "neutral_excess_vector_norm_squared":str(raw[0,0]),
        "noncommuting_return_pairs":[list(pair) for pair in noncommuting],
        "maximum_pairwise_commuting_subset":2,
        "return_word_span_by_maximum_length":word_spans,
        "positive_density_holdout":{"word_indices":[i,intermediate,j],
                                    "contrast":str(contrast),
                                    "all_shorter_words_agree":True},
        "physical_source_identity_proved":False,
        "weighted_r113_values":None,
        "next_original_observation":"Read all six marked returns with one absolute reference; separate "
                                    "neutral grade two, charged grade three, and the physical support transport.",
    }
