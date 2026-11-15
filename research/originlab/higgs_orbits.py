"""Classify native Fourier histories under constant gauge and spatial chart changes."""
from collections import defaultdict
from itertools import combinations, combinations_with_replacement
import json

import sympy as sp

from higgs_moment_search import basis_rows, field_coefficients, probe_plan
from reconstruction import require
from source import MARKS, as_form


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


def support_polynomial(points, labels, basis, constant_coordinates):
    """Describe a nonnegative polynomial vanishing on the planned signatures."""
    used = sorted({j for word,_ in basis for j in word})
    if not used and not constant_coordinates:
        return None
    scales = [max(x[j] for x in points)-min(x[j] for x in points) for j in used]
    return {
        "coordinates":[labels[j] for j in used],
        "coordinate_scales":list(map(str,scales)),
        "zeros":[[str(x[j]) for j in used] for x in points],
        "constant_coordinates":constant_coordinates,
        "formula":("sum of squared deviations on constant coordinates plus the product over zeros s "
                   "of squared scaled distances; omit the product when there are no varying coordinates"),
        "degree_in_coordinates":max(2 if constant_coordinates else 0,2*len(points) if used else 0),
        "bounded_equivalent":"P/(1+P)",
        "required_measurement":"integrate the nonnegative test against the same positive scalar weight",
        "zero_scope":"support of the pushforward on the listed coordinates only",
    }


def orbit_measurements(grams, orbits, edges):
    """Build an exact moment inverse, allowing unresolved orbit classes."""
    generators = []
    signatures = set()
    constants = []

    def add(label, values):
        for half,part in enumerate(("real","imaginary")):
            values_by_orbit = []
            for group in orbits:
                orbit_values = {values[gene].as_real_imag()[half] for gene in group}
                if len(orbit_values)!=1:
                    raise ArithmeticError("a probe depends on the spatial or gauge coordinates")
                values_by_orbit.append(orbit_values.pop())
            column = tuple(values_by_orbit)
            if len(set(column))>1 and column not in signatures:
                signatures.add(column)
                generators.append(({**label,"component":part},column))
            elif len(set(column))==1 and label["kind"]=="reflected_diagonal" and part=="real":
                constants.append({**label,"component":part,"value":str(column[0])})

    for n,m in edges.get(0, []):
        add({"kind":"reflected_diagonal","mode":n},
            {genes:(g[n,m]+g[-n,-m])/2 for genes,g in grams.items()})
    for bucket in edges.values():
        for first,second in combinations_with_replacement(bucket,2):
            add({"kind":"reflected_gram_product","first":first,"second":second},
                {genes:even_word(g,first,second) for genes,g in grams.items()})

    classes = {}
    for i,orbit in enumerate(orbits):
        signature = tuple(column[i] for _,column in generators)
        classes.setdefault(signature,[]).append(orbit)
    labels = [label for label,_ in generators]
    basis = basis_rows(list(classes),labels)
    matrix = sp.Matrix([values for _,values in basis])
    return {"separates_all_orbits":len(classes)==len(orbits),
            "probe_classes":list(classes.values()),
            "independent_moments_including_empty":len(basis),
            "moments":[{"word":[labels[j] for j in word],"value_per_class":list(map(str,values))}
                       for word,values in basis],
            "class_totals_from_ordered_measured_moments":[list(map(str,row)) for row in matrix.inv().tolist()],
            "positive_support_test":support_polynomial(list(classes),labels,basis,constants),
            "scope":"each measured moment uses one insertion of the entire product on the same weighted fibre"}


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
            "weighted_orbit_measurement_plan":orbit_measurements(grams,orbits,edges),
            "scope":"fixed flat transport chart on the spatial Higgs doublet; neither a physical OS weight nor a local gauge quotient"}


def reconstruct_marked_forms(data):
    """Reconstruct the four forms per probe class from jointly supplied moments."""
    proposal = data["field_proposal"]
    plan = orbit_plan(proposal,probe_plan(proposal)["norm_signature_classes"])["weighted_orbit_measurement_plan"]
    measurements = data["measurements"]
    size = data["word_space_size"]
    require(type(size) is int and size>0,"positive word-space dimension required")
    expected = [row["word"] for row in plan["moments"]]
    require(json.dumps([row["word"] for row in measurements],sort_keys=True)==json.dumps(expected,sort_keys=True),
            "measurements must identify the exact words in plan order")
    require(all(set(row["ordered_source_forms"])==set(MARKS) for row in measurements),
            "each measurement must supply all four marked forms")
    inverse = sp.Matrix(plan["class_totals_from_ordered_measured_moments"]).applyfunc(sp.Rational)
    forward = sp.Matrix([row["value_per_class"] for row in plan["moments"]]).applyfunc(sp.Rational)
    recovered = [dict() for _ in plan["probe_classes"]]
    for mark in MARKS:
        observed = [as_form(row["ordered_source_forms"][mark],size) for row in measurements]
        forms = [sum((weight*item for weight,item in zip(row,observed)),sp.zeros(size))
                 for row in inverse.tolist()]
        for row,target in zip(forward.tolist(),observed):
            if sum((weight*item for weight,item in zip(row,forms)),sp.zeros(size))!=target:
                raise ArithmeticError("marked form inversion failed replay")
        for record,form in zip(recovered,forms):
            record[mark] = [[[str(x) for x in value.as_real_imag()] for value in row] for row in form.tolist()]
    return {"status":"EXACT_FORMS_ON_DECLARED_PROBE_CLASSES","probe_classes":plan["probe_classes"],
            "reconstructed_marked_forms":recovered,"physical_source_identity_proved":False}
