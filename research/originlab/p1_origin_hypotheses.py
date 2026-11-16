"""Replay origin hypotheses at the stage where their premises apply."""

from itertools import product

import sympy as sp

from affine_operators import require


def reconstruction():
    # Enumerate the two-element unital multiplication laws at the first
    # cardinality with any element distinct from the unit.
    monoids = []
    for square in (0, 1):
        table = [[0, 1], [1, square]]
        if all(table[table[a][b]][c] == table[a][table[b][c]]
               for a, b, c in product(range(2), repeat=3)):
            monoids.append(table)
    require(len(monoids) == 2 and monoids[0][1][1] != monoids[1][1][1],
            "minimal unital laws must retain both inequivalent squares")
    # The unit must be fixed, so no relabelling identifies the two tables.
    return {"minimum_nontrivial_cardinality_in_class": 2,
            "inequivalent_associative_laws": monoids,
            "repeated_rewrite_orbit_lengths_from_generator": [2, 1]}


def error_correction_readout():
    identity, x, z = sp.eye(2), sp.Matrix([[0, 1], [1, 0]]), sp.diag(1, -1)
    axes = [sp.kronecker_product(a, b) for a, b in
            ((x, identity), (z, identity), (identity, x), (x, x))]
    effects = [(sp.eye(4) + a/4)/5 for a in axes]
    effects.append(sp.eye(4) - sum(effects, sp.zeros(4)))
    hidden = sp.kronecker_product(identity, z)
    states = [(sp.eye(4) + sign*hidden)/4 for sign in (1, -1)]
    rank = sp.Matrix.hstack(*(effect.reshape(16, 1) for effect in effects)).rank()
    require(rank == 5 and all(effect.is_positive_semidefinite for effect in effects),
            "five-effect control must be a positive rank-five POVM")
    require(all(state.is_positive_semidefinite and sp.trace(state) == 1 for state in states),
            "valid logical states")
    outputs = [[sp.trace(effect*state) for effect in effects] for state in states]
    require(outputs[0] == outputs[1] and states[0] != states[1], "tomography witness")
    return {"readout_rank": rank, "full_logical_operator_dimension": 16,
            "indistinguishable_probabilities": [str(value) for value in outputs[0]],
            "separating_logical_expectations": [str(sp.trace(hidden*state)) for state in states]}


def relational_order():
    forward = {(0, 1), (0, 2)}
    reversed_order = {(b, a) for a, b in forward}
    support = lambda edges: {frozenset(pair) for pair in edges}
    minimal = lambda edges: sorted(set(range(3)) - {b for a, b in edges})
    require(support(forward) == support(reversed_order)
            and len(minimal(forward)) != len(minimal(reversed_order)), "orientation witness")
    return {"common_undirected_edges": [sorted(edge) for edge in sorted(forward)],
            "different_minimal_event_sets": [minimal(forward), minimal(reversed_order)]}


def modular_normalization():
    orbit = sp.diag(1, sp.I, -1, -sp.I)
    probe = sp.diag(1, 0, 0, 0)
    beta = 2*sp.pi
    responses = []
    for ratio in (sp.Rational(1, 2), sp.Rational(1, 3)):
        energy = -sp.log(ratio)/beta
        weights = sp.diag(1, *(sp.exp(-beta*energy) for _ in range(3)))
        rho = weights/sp.trace(weights)
        require(rho.is_positive_definite and orbit*rho == rho*orbit and orbit**4 == sp.eye(4)
                and orbit**2 != sp.eye(4), "faithful Gibbs state with Z4 symmetry")
        responses.append(str(sp.trace(probe*rho)))
    require(responses == ["2/5", "1/2"], "fixed beta does not fix the response")
    return {"beta": str(beta), "symmetry_order": 4,
            "same_mark_different_expectations": responses,
            "scope": "Z4 symmetry and fixed Gibbs beta alone, with independently varied gap"}


def conic_defects(max_order):
    # On C/C_n the trace contribution from rotation j is
    # 1/(n det(I-R_j)). Fourier transformation to character k gives
    # the displayed rational. Its centered second difference is
    # 1/n off k=0 and -(n-1)/n at k=0 by character orthogonality.
    census = []
    for order in range(2, max_order + 1):
        coeffs = [sp.Rational(order**2 - 1 - 6*k*(order-k), 12*order) for k in range(order)]
        require(sum(coeffs) == 0 and all(
            coeffs[(k+1) % order] - 2*coeffs[k] + coeffs[(k-1) % order]
            == sp.Rational(1-order if k == 0 else 1, order)
            for k in range(order)), "equivariant cone trace recurrence")
        spin = [sign**order * -1 for sign in (1, -1)]
        census.append({"order": order, "deficit_in_turns": str(1 - sp.Rational(1, order)),
                       "fixed_point_coefficients": list(map(str, coeffs)),
                       "plain_spin_lift_descends": 1 in spin})
    return {
        "tested_orders": [2, max_order],
        "twisted_sectors": sum(row["order"] for row in census),
        "census": census,
        "order_four": next(row for row in census if row["order"] == 4),
        "scope": "local flat cyclic cone, scalar character heat trace and bare Spin(2) lift",
    }


def orientation_cover():
    from scenarios import boundaries, signature

    mobius = signature(*boundaries("moebius", 3))
    cover = signature(*boundaries("cylinder", 6))
    require(mobius["betti_Q"] == cover["betti_Q"] == [1, 1, 0]
            and not mobius["orientable"] and cover["orientable"]
            and [mobius["boundary_components"], cover["boundary_components"]] == [1, 2],
            "Möbius orientation double cover invariants")
    return {"mobius": mobius, "cylinder_cover": cover}


def hopf_readout():
    def hopf(pair):
        left, right = pair
        return (2*left*right.conjugate(), left.norm()**2 - right.norm()**2)

    a = sp.Quaternion(1, 1, 0, 0)/2
    b = sp.Quaternion(0, 0, 1, 1)/2
    fiber_phases = [sp.Quaternion(-1, 0, 0, 0), sp.Quaternion(0, 1, 0, 0)]
    require(a.norm()**2 + b.norm()**2 == 1 and all(
        hopf((a*phase, b*phase)) == hopf((a, b)) for phase in fiber_phases
    ), "Hopf fiber witness")
    return {"fiber_states_with_equal_readout": 3,
            "fiber_phase_omitted": True,
            "reason": "For any unit quaternion h, (a h) conjugate(b h)=a conjugate(b)."}


def projection_fold():
    y = sp.symbols("y", positive=True)
    z = sp.symbols("z", positive=True)
    jacobian = (y*y).diff(y)
    preimage = sp.sqrt(z)
    branches = [sp.simplify(1/jacobian.subs(y, sign*preimage)) for sign in (1, -1)]
    require(sum(branches) == 0 and sum(map(abs, branches)) == 1/sp.sqrt(z),
            "signed and absolute fold responses")
    return {"map": "y -> y^2", "signed_sum": "0", "unit_density_response": "1/sqrt(z)",
            "free_source_input": "density on the two sheets"}


FAMILIES = {
    "self_interpretation": ("reconstruction", "unique minimal physical quotient"),
    "global_obstruction": ("prephysical", "phase-stable commutator plus a native boundary map"),
    "larger_parent": ("geometry", "visible marked quotient and parent discriminator"),
    "self_reproduction": ("reconstruction", "inequivalent minimal rewrite laws"),
    "quantum_error_correction": ("logical_readout", "recovery on the specified logical algebra"),
    "relational_time": ("relational", "orientation and positive physical time"),
    "protected_defect": ("cone", "primitive spin and the complete character heat response"),
    "symmetric_singularity": ("cone", "independent conversion from defect trace to native c3"),
    "torus_and_double_cover": ("geometry", "full deck trace, spin lift and genus discriminator"),
    "spherical_holographic_reader": ("geometry", "full covariance and radial energy at the same scale"),
    "p1_as_transformation": ("geometry", "ordered quotient-switch measurements and the native clock"),
    "modular_index_kms": ("modular", "absolute marked response at fixed symmetry and temperature"),
    "mobius_orientation_cover": ("orientation_cover", "lift of the oriented marked boundary pairing"),
    "hopf_holographic_reader": ("hopf", "coherent recovery of a fiber phase"),
    "projection_singularity": ("fold", "signed versus weighted pushforward across both sheets"),
    "minimal_self_dual_relations": ("self_dual_relations", "an independently justified binary self-dual composition and a positive reconstructed state"),
    "directed_anchor_groupoid": ("orientation_flags", "native action on the oriented flag and its absolute state"),
    "fourier_reciprocal_sources": ("reciprocal_sources", "a native fourth cumulant or independently justified composition law"),
    "independent_dyadic_composition": ("reciprocal_sources", "independence and exact composition on the same source, including its ordered markers"),
    "coherent_fourier_composition": ("coherent_gaussian", "native pairing, purity, two-copy mixer and the primitive relative marker"),
    "coherent_additive_generator": ("coherent_number", "native additive generator, physical time and absolute empty transfer"),
    "coherent_ordered_markers": ("coherent_words", "match original U and Y with the outer cap and even Majorana intermediates"),
    "coherent_four_code_decoder": ("coherent_decoder", "energy-covariant decoder and an observable map retaining reservoir coherence"),
}


def geometry():
    from f1_origin_scenarios import controls

    result = controls()
    require(result["F1_solved"] is False, "comparison replay must preserve its source scope")
    return result


def search(hypotheses=None, max_cone_order=32):
    if hypotheses is None:
        hypotheses = list(FAMILIES)
    if (not isinstance(hypotheses, list) or not hypotheses
            or any(not isinstance(key, str) or key not in FAMILIES for key in hypotheses)):
        raise ValueError("select one or more registered P1 hypotheses")
    if type(max_cone_order) is not int or not 4 <= max_cone_order <= 128:
        raise ValueError("maximum cone order must be between four and 128")
    from f1_prephysical_interpretation import controls as prephysical
    from coherent_origin import decoder_controls, gaussian_covariance, number_sectors, words_control
    from p1_origin_report_controls import orientation_flags, reciprocal_sources, self_dual_relations

    checks = {
        "reconstruction": reconstruction, "prephysical": prephysical,
        "geometry": geometry, "logical_readout": error_correction_readout,
        "relational": relational_order, "cone": lambda: conic_defects(max_cone_order),
        "modular": modular_normalization, "orientation_cover": orientation_cover,
        "hopf": hopf_readout, "fold": projection_fold,
        "self_dual_relations": self_dual_relations, "orientation_flags": orientation_flags,
        "reciprocal_sources": reciprocal_sources,
        "coherent_gaussian":gaussian_covariance,"coherent_number":lambda: number_sectors(6),
        "coherent_words":words_control,"coherent_decoder":decoder_controls,
    }
    reports, alternatives = {}, {}
    for hypothesis in dict.fromkeys(hypotheses):
        evidence, probe = FAMILIES[hypothesis]
        if evidence not in reports:
            reports[evidence] = checks[evidence]()
        alternatives[hypothesis] = {"comparison": evidence, "next_probe": probe}
    return {
        "status": "SCOPED_P1_HYPOTHESES_TESTED",
        "hypotheses": alternatives,
        "evidence": reports,
        "selection": "Compare finite prephysical laws first; require positivity and time at physical interpretation.",
        "original_source_binding": None,
        "physical_solution": False,
    }
