"""Classify native Fourier histories under constant gauge and spatial chart changes."""
from collections import defaultdict
from itertools import combinations, combinations_with_replacement

import sympy as sp

from higgs_moment_search import field_coefficients


def gram_for_field(field, modes):
    zero = sp.zeros(4,1)
    doublets = {}
    for n in modes:
        a,b,c,d = field.get(n,zero)
        doublets[n] = sp.Matrix([-a-sp.I*b,c+sp.I*d])
    return {(n,m):sp.expand((v.H*w)[0]) for n,v in doublets.items() for m,w in doublets.items()}


def phase_compatible(left, right, reflect=False):
    """Check G'_{nm}=z^(m-n) G_{s n,s m}, including a common unit z."""
    phases = {}
    direction = -1 if reflect else 1
    for (n,m), actual in right.items():
        if n > m:
            continue
        expected = left[direction*n,direction*m]
        if n == m or not (actual and expected):
            if actual != expected:
                return False
            continue
        exponent = m-n
        phase = sp.cancel(actual/expected)
        if sp.expand(phase*sp.conjugate(phase))!=1 or phases.setdefault(exponent,phase)!=phase:
            return False
    if not phases:
        return True
    generator,phase = next(iter(phases.items()))
    for exponent, value in phases.items():
        a,b,gcd = sp.gcdex(generator,exponent)
        phase = sp.cancel(phase**a*value**b)
        generator = gcd
    return all(sp.cancel(phase**(n//generator)-value)==0 for n,value in phases.items())


def majorana_factor(couplings):
    def matrix(name):
        return sp.Matrix([[sp.Rational(str(a))+sp.I*sp.Rational(str(b)) for a,b in row]
                          for row in couplings[name]])
    family = matrix("Ynu").conjugate()*matrix("MR").H
    return sp.simplify(sp.trace(family.H*family))


def orbit_classes(grams, reflection):
    orbits = []
    orbit_by_gene = {}
    # Sorted diagonal weights are unchanged under rotation and reflection.
    buckets = defaultdict(list)
    for genes,g in grams.items():
        key = tuple(sorted(v for (n,m),v in g.items() if n==m))
        for index in buckets[key]:
            representative = grams[orbits[index][0]]
            if phase_compatible(representative,g) or (reflection and phase_compatible(representative,g,True)):
                orbits[index].append(genes)
                break
        else:
            index = len(orbits)
            orbits.append([genes])
            buckets[key].append(index)
        orbit_by_gene[genes] = index
    return orbits,orbit_by_gene


def even_word(gram, first, second):
    n,m = first
    i,j = second
    direct = gram[n,m]*sp.conjugate(gram[i,j])
    reflected = gram[-n,-m]*sp.conjugate(gram[-i,-j])
    return sp.expand((direct+reflected)/2)


def orbit_plan(proposal, norm_classes):
    genes = [tuple(g) for group in norm_classes for g in group]
    fields = {g:field_coefficients(proposal,g) for g in genes}
    modes = sorted(set().union(*(f.keys() for f in fields.values())))
    grams = {g:gram_for_field(field,modes) for g,field in fields.items()}
    oriented,_ = orbit_classes(grams,False)
    orbits,index = orbit_classes(grams,True)
    ambiguous = {}
    for group in norm_classes:
        reps = {}
        for g in map(tuple,group):
            reps.setdefault(index[g],g)
        for i,j in combinations(sorted(reps),2):
            ambiguous.setdefault((i,j),(reps[i],reps[j]))
    edges = defaultdict(list)
    for pair in combinations_with_replacement(modes,2):
        if any(gram[pair] for gram in grams.values()):
            edges[pair[1]-pair[0]].append(pair)
    best = None
    for bucket in edges.values():
        for first,second in combinations_with_replacement(bucket,2):
            values = {g:even_word(gram,first,second) for g,gram in grams.items()}
            pairs = [(a,b) for a,b in ambiguous.values() if values[a]!=values[b]]
            if pairs and (best is None or len(pairs)>best["separated_pairs"]):
                a,b = pairs[0]
                best = {"first":first,"second":second,"separated_pairs":len(pairs),
                        "example_genes":[a,b],"example_values":[str(values[a]),str(values[b])]}
    # This closed word directly tests the assumed quarter-translation
    # identification when the native sine and cosine amplitudes differ.
    sign_word = None
    if (len(proposal["slots"]) == 3 and (-1,1) in grams[genes[0]] and (0,2) in grams[genes[0]]):
        a,b = (2,2,0),(2,2,2)
        if a in grams and b in grams:
            values = [even_word(grams[g],(-1,1),(0,2)) for g in (a,b)]
            sign_word = {"genes":[a,b],"values":list(map(str,values)),
                         "same_gauge_spatial_orbit":index[a]==index[b]}
    return {"constant_gauge_and_rotation_orbits":len(oriented),
            "including_spatial_reflection_orbits":len(orbits),"norm_equal_distinct_orbit_pairs":len(ambiguous),
            "best_quartic_word":best,"second_harmonic_sign_witness":sign_word,
            "scope":"fixed flat transport chart on the spatial Higgs doublet; neither a physical OS weight nor a local gauge quotient"}
