"""Exhaust the declared Full96 partner and order choices, without fitting weights."""

from fractions import Fraction
from itertools import combinations, product
from math import ceil, comb

import numpy as np
import sympy as sp

from f1_reflected_cap import load_assembler
from reconstruction import require
from r113_three_window_census import census as three_window_census
from r113_two_window_census import numeric_fixed_width, series_census


PRIME = 65537


def residue(value):
    value = Fraction(value)
    return value.numerator*pow(value.denominator, -1, PRIME) % PRIME


def gaussian_residues(matrix):
    converted = np.zeros(matrix.shape, dtype=np.int64)
    for i, j in zip(*np.nonzero(matrix)):
        value = complex(matrix[i, j])
        # 256² = -1 mod 65537; exact dyadics come from the original fixture.
        converted[i, j] = (residue(value.real)+256*residue(value.imag)) % PRIME
    return converted


def generators_mod_prime(assembler):
    domain = assembler.domain
    nodes = domain["modes"]
    identity = np.eye(len(nodes), dtype=np.int64)
    momenta = np.diag([residue(Fraction(n)+Fraction(domain["quarter"])) for n in nodes])
    gamma = gaussian_residues(assembler.gamma)
    variations = [gaussian_residues(-1j*assembler.gamma@v) for v in assembler.V]
    constant = gaussian_residues(-1j*assembler.gamma@assembler.DM)
    for coefficient, variation in zip(domain["field_proposal"]["base_constant"], variations):
        constant = (constant+residue(coefficient)*variation) % PRIME
    inputs = [(np.kron(momenta,gamma)+np.kron(identity,constant)) % PRIME]
    differences = np.subtract.outer(nodes, nodes)
    for slot in domain["field_proposal"]["slots"]:
        harmonic = slot["harmonic"]
        unit = Fraction(slot["unit"])
        if slot["kind"] == "constant":
            spatial = residue(unit)*identity
        else:
            positive, negative = differences == harmonic, differences == -harmonic
            if slot["kind"] == "cos":
                spatial = residue(unit/2)*(positive.astype(np.int64)+negative)
            else:
                spatial = 256*residue(unit/2)*(negative.astype(np.int64)-positive)
        inputs.append(np.kron(spatial, variations[slot["component"]]) % PRIME)
    return inputs


def commutator_certificate(generators):
    pairs = list(combinations(range(len(generators)), 2))
    columns = np.column_stack([
        ((generators[i]@generators[j]-generators[j]@generators[i]) % PRIME).ravel()
        for i,j in pairs
    ])
    basis = columns.copy()
    pivots = []
    for j in range(len(pairs)):
        nonzero = np.flatnonzero(basis[:, j])
        if not len(nonzero):
            break
        pivot = int(nonzero[0])
        pivots.append(pivot)
        divisor = pow(int(basis[pivot,j]), -1, PRIME)
        basis[:,j] = (divisor*basis[:,j]) % PRIME
        for k in range(j+1,len(pairs)):
            basis[:,k] = (basis[:,k]-basis[pivot,k]*basis[:,j]) % PRIME
    witness = columns[pivots,:].tolist()
    if len(pivots) == len(pairs):
        require(int(sp.Matrix(witness).det()) % PRIME != 0, "noninvertible commutator witness")
    return {"prime": PRIME, "independent_pairs_certified": len(pivots),
            "pair_labels": [list(pair) for pair in pairs], "matrix_entry_indices": pivots,
            "residue_minor": witness}


def certified_excess_bound(assembler, multipliers, width):
    """Conservative rational degree-one lower bound, with no spectral rounding."""
    domain = assembler.domain
    count = len(domain["modes"])
    require(sorted(domain["modes"]) == list(range(min(domain["modes"]), min(domain["modes"])+count)),
            "Frobenius weights require consecutive spatial Fourier modes")
    settings = domain["field_proposal"]
    frequencies = [(slot["kind"],slot["harmonic"]) for slot in settings["slots"]]
    require(len(set(frequencies)) == len(frequencies), "separate cosine and sine directions required for the gap certificate")
    require(all((kind == "constant" and harmonic == 0)
                or (kind in ("cos","sin") and 0 < harmonic < count) for kind,harmonic in frequencies),
            "gap certificate requires nonzero present Fourier diagonals")

    def l1(value):
        return abs(Fraction(float(value.real)))+abs(Fraction(float(value.imag)))

    def row_bound(matrix):
        return max(sum(map(l1, row)) for row in matrix)

    peak = max(abs(mode+Fraction(domain["quarter"])) for mode in domain["modes"])
    peak += row_bound(assembler.DM)
    norms = [row_bound(v) for v in assembler.V]
    peak += sum(abs(Fraction(c))*norm for c,norm in zip(settings["base_constant"],norms))
    peak += max(map(abs,multipliers))*sum(abs(Fraction(slot["unit"]))*norms[slot["component"]]
                                          for slot in settings["slots"])
    spacing = min(abs(a-b) for a,b in combinations(multipliers,2))
    squared_distances = []
    for slot in settings["slots"]:
        harmonic = slot["harmonic"]
        matrix = assembler.V[slot["component"]]
        frobenius_squared = sum(Fraction(float(v.real))**2+Fraction(float(v.imag))**2
                                for v in matrix.flat if v)
        repetitions = Fraction(count) if harmonic == 0 else Fraction(count-harmonic,2)
        squared_distances.append(repetitions*(spacing*Fraction(slot["unit"]))**2*frobenius_squared)
    # Different harmonics are Frobenius orthogonal. At one harmonic, cos
    # and sin cross terms cancel between the two ordered Fourier diagonals.
    distance_squared = min(squared_distances)
    bound = width**2*distance_squared/(2*count*len(assembler.gamma)*2**ceil(20*width*peak))
    require(bound > 0, "profile mismatch needs a positive certified gap")
    return {"rational_lower_bound":str(bound),"lower_bound_decimal":float(bound),
            "generator_norm_upper_bound":str(peak),"formula":"w^2 * min_distance_squared * 2^-ceil(20*w*norm_upper) / (2*d)",
            "mixture_inversion":"mismatch_mass <= min(1, (mean_degree_one_overlap - 1) / rational_lower_bound)"}


def reverse_order_negative_control(dimension):
    a = sp.diag(2,sp.Rational(1,2))
    b = sp.Matrix([[5,3],[3,5]])/4
    require(all(m.det() == 1 and m.is_positive_definite for m in (a,b)), "positive determinant-one control")
    reversed_product = a*b*a.inv()*b.inv()
    overlap = sp.trace(reversed_product)/2
    require(overlap == sp.Rational(47,128), "multiwindow order control")
    single_overlap = sp.trace(a)/2
    mass = (1-overlap)/(single_overlap-overlap)
    identity = sp.eye(2)
    single_defect, reversed_defect = (matrix-identity for matrix in (a,reversed_product))
    defect = (mass*sp.trace(single_defect.H*single_defect)
              +(1-mass)*sp.trace(reversed_defect.H*reversed_defect))/2
    require(0<mass<1 and mass*single_overlap+(1-mass)*overlap==1 and defect>0,
            "unit mean overlap must coexist with a strictly positive same-history defect")
    return {"two_dimensional_overlap":str(overlap),
            "embedded_overlap":str((dimension-2+2*overlap)/dimension),
            "unit_mean_overlap_wrong_history_masses":[str(mass),str(1-mass)],
            "unit_mean_overlap_squared_defect":str(defect)}


def pairing_defect_control(transfers, profile_index, mirrors, pairing, witness):
    def ordered_word(indices):
        early, middle, late = (transfers[i] for i in indices)
        return late@middle@early

    left, other = ([profile_index[tuple(profile)] for profile in witness[key]]
                   for key in ("left","other_left"))
    original = ordered_word(left)
    errors = {}
    for name, history in (("ordered_mirror",left),("quadratically_aliased_mirror",other)):
        right = ordered_word([mirrors[i] for i in history])
        delta = original@pairing@right.T-pairing
        errors[name] = float(np.vdot(delta,delta).real/len(pairing))
    require(errors["ordered_mirror"] < 1e-22 and errors["quadratically_aliased_mirror"] > 1e-12,
            "positive squared-defect control disagrees with the exact cubic witness")
    return {"normalized_squared_defects":errors,"comparison_width":"1/10",
            "exact_selection_rule":"zero expected squared defect under a positive common history law",
            "physical_source_measure_evaluated":False}


def mirrored_dyadic_control(slots):
    signs = [-1 if slot["kind"] == "sin" else 1 for slot in slots]
    mirror = sp.diag(*signs)
    identity = sp.eye(len(signs))
    proposed_merge = (identity+mirror)/sp.sqrt(2)
    odd = (identity-mirror)/2
    even = (identity+mirror)/2
    require(odd*proposed_merge == sp.zeros(len(signs)) and even*proposed_merge == sp.sqrt(2)*even,
            "a matched dyadic merge must erase odd coordinates and scale even coordinates")
    # The fixed law is even. Its second moment would then be twice itself.
    return {"mirror_signs":signs,"second_moment_factors":[(1+s)**2//2 for s in signs],
            "sole_finite_second_moment_fixed_law":"delta_at_zero_profile_fluctuation",
            "scope":"linear normalized sum of the two matched spatial coefficient vectors"}


def search(root):
    assembler = load_assembler(root)
    fixture = assembler.domain
    slots = fixture["field_proposal"]["slots"]
    coefficients = tuple(map(Fraction, fixture["field_proposal"]["allowed_multipliers"]))
    profiles = tuple(product(range(len(coefficients)), repeat=len(slots)))
    generators = generators_mod_prime(assembler)
    certificate = commutator_certificate(generators)
    require(certificate["independent_pairs_certified"] == comb(len(slots)+1, 2),
            "cannot classify all reversed histories without a full commutator certificate")
    index = {genes: i for i, genes in enumerate(profiles)}
    mirrors = [index[tuple(coefficients.index(-coefficients[g]) if slot["kind"] == "sin" else g
                           for g, slot in zip(genes, slots))] for genes in profiles]
    require(len(set(mirrors)) == len(profiles) and all(mirrors[j] == i for i,j in enumerate(mirrors)),
            "reflection must close the input menu")
    width = Fraction(1,10)
    transfer_width = float(width)
    transfers = []
    for genes in profiles:
        h = assembler.assemble({"schema":"tfpt.field-configuration.v1",
                                 "domain_sha256":fixture["_sha256"],"genes":list(genes)})["h_N"]
        require(np.trace(h) == 0 and np.array_equal(h,h.conj().T), "positive determinant-one comparison requires traceless Hermitian h")
        eigenvalues, basis = np.linalg.eigh(h)
        transfers.append((basis*np.exp(-transfer_width*eigenvalues))@basis.conj().T)
    pairing = np.kron(np.eye(len(fixture["modes"])), assembler.gamma@assembler.K)
    dimension = len(pairing)
    left = np.stack(transfers).reshape(len(profiles), -1)
    weighted_right = np.stack([(pairing@r.T@pairing.T).T for r in transfers]).reshape(len(profiles),-1)
    overlaps = (left@weighted_right.T)/dimension
    mask = np.ones(overlaps.shape,dtype=bool)
    mask[np.arange(len(profiles)),mirrors] = False
    gap = float(np.min(overlaps.real[mask])-1)
    sharpness = certified_excess_bound(assembler, coefficients, width)
    require(np.max(np.abs(overlaps.imag)) < 1e-12 and gap > 1e-8
            and np.max(np.abs(overlaps.real[~mask]-1)) < 1e-12 and gap >= sharpness["lower_bound_decimal"],
            "native degree-one overlap disagrees with exact matching")

    rows = [{"left":list(genes),"unique_partner":list(profiles[mirrors[i]])}
            for i,genes in enumerate(profiles)]
    n = len(rows)
    histories, joint_census = series_census(profiles, mirrors, coefficients, certificate["independent_pairs_certified"])
    joint_census["fixed_width_replay"] = numeric_fixed_width(transfers, histories, certificate["matrix_entry_indices"])
    vectors = np.array([(1,*(residue(coefficients[g]) for g in genes)) for genes in profiles],dtype=np.int64)
    triple_census = three_window_census(generators,vectors,profiles,mirrors,PRIME)
    # Each exact affine generator is over Q(i). Spectral projectors and
    # eigenvalues are algebraic. A separating jet makes at least one
    # grouped coefficient of a finite exponential sum nonzero in a
    # matrix entry. At nonzero algebraic width the distinct algebraic
    # exponents give independent exponentials by Lindemann–Weierstrass.
    algebraic_width_certificate = {
        "theorem":"Lindemann–Weierstrass",
        "scope":"exact affine lift with dyadic internal entries and declared rational slot weights",
        "width_class":"positive algebraic common width",
        "tested_width":str(width),
        "partner_counts_by_window_count":{"1":n,"2":len(histories),"3":triple_census["equal_length_partner_count"]},
        "independent_legs_allow_nontrivial_mixture_if_pairing_holds_per_history":False,
        "independence_reason":"a product of two supports lies on the graph of a bijection only for singleton supports",
    }
    triple_census["pairing_defect_control"] = pairing_defect_control(
        transfers,index,mirrors,pairing,triple_census["quadratic_negative_control"])
    return {"status":"EXACT_SCOPED_SEARCH_COMPLETED", "physical_solution":False,
            "native_domain_digest":fixture["_sha256"], "one_particle_dimension":dimension,
            "one_step":{"candidate_pairs":n*n,"preserving_pairs":n,"excluded_pairs":n*n-n,
                        "same_profile_preserving":sum(i==j for i,j in enumerate(mirrors)),"partners":rows},
            "two_step_individually_matched":{"left_histories":n*n,"same_order_preserving":n*n,
                                            "reverse_order_preserving":n,"reverse_order_excluded":n*n-n,
                                            "commutator_certificate":certificate},
            "two_step_independently_chosen_partners":joint_census,
            "three_step_independently_chosen_partners":triple_census,
            "algebraic_width_certificate":algebraic_width_certificate,
            "mirrored_dyadic_composition":mirrored_dyadic_control(slots),
            "multiwindow_positivity_negative_control":reverse_order_negative_control(dimension),
            "positive_single_step_ensembles":{"strict_excess_for_every_interior_pure_paired_degree_on_mismatch":True,
                                              "proof":"AM-GM of all k-subset eigenvalue products, equality only for identity pairing",
                                              "numeric_smallest_mismatch_degree_one_excess":gap,
                                              "degree_one_mismatch_mass_certificate":sharpness,
                                              "comparison_width":transfer_width},
            "scope":"fixed S, declared spatial fixture and common step width; reverse-order census has individually matched partners; no OS time or native outer-cap claim"}
