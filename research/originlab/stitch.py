"""Compose reviewed conditional results using conjunctive and source-bound premises."""
from collections import defaultdict

RULES = (
    ("tetrahedral_readout",("a3_standard_weight_frame",),"tetrahedral_sic_readout"),
    ("ordered_charge",("original_weighted_form","original_absolute_reference",
                      "admissible_charged_field","quarter_covariance_on_field",
                      "quarter_preserves_state_and_transfer"),"charged_quarter_response"),
    ("primitive_monodromy",("original_weighted_form","admissible_charged_field",
                            "full_order_family_character"),"primitive_family_insertion"),
    ("r102",("primitive_family_insertion","one_sided_positive_standard_integral_pairing",
             "local_weight_one_chiral_current","original_operator_product_closure"),"physical_e8_current_extension"),
    ("modular_regional",("original_weighted_form","original_absolute_reference",
                        "common_cyclic_separating_state","compatible_half_sided_modular_inclusions",
                        "regional_algebra_transport_and_joint_commutation"),"strongly_commuting_regional_translations"),
    ("rotation_spectrum",("strongly_commuting_regional_translations","positive_null_direction",
                          "continuous_transitive_spatial_rotations"),"joint_forward_spectrum"),
    ("joint_block_spectrum",("strongly_commuting_regional_translations",
                             "original_operator_valued_cone_positivity"),"joint_forward_spectrum"),
    ("common_structure",("physical_e8_current_extension","charged_quarter_response",
                         "joint_forward_spectrum"),"common_chiral_and_regional_origin"),
    ("complete_candidate",("common_chiral_and_regional_origin","complete_local_operator_generation",
                           "controlled_interacting_continuum","physical_chiral_matter_no_mirrors",
                           "massless_universally_coupled_spin2","nonlinear_ward_identities",
                           "independent_empirical_validation"),"full_tfpt_candidate"),
)

FOUR_MARK_RULES = (
    ("selected_pairing",("actual_face_operations","same_source_positive_cyclic_space",
                         "native_oriented_face_intertwiner","original_clifford_relations"),
     "four_mark_quantum_pairing"),
    ("selected_charge",("actual_charged_loop_operations","same_source_positive_cyclic_space",
                        "original_integral_cartan_action","full_order_family_character"),
     "primitive_family_insertion"),
    ("common_four_mark_bridge",("four_mark_quantum_pairing","primitive_family_insertion",
                               "nonzero_original_ordered_return","original_weighted_form",
                               "original_absolute_reference"),
     "physical_four_mark_bridge"),
)

ORDERED_RETURN_RULES = (
    ("affine_line_in_original_process",
     ("original_weighted_form","original_absolute_reference","native_primitive_charged_current",
      "intertwines_native_charge_adjoint_and_physical_time"),
     "original_primitive_current_lift"),
    ("measured_ordered_return",
     ("original_primitive_current_lift","original_complement_time_transfer",
      "nonzero_independently_matched_ordered_return"),
     "original_dynamical_return"),
    ("same_source_regional_charged_dynamics",
     ("original_dynamical_return","one_sided_positive_standard_integral_pairing",
      "original_reference_and_relational_time_ratio","compatible_regional_gluing"),
     "original_ordered_dynamical_sewing"),
)


def derive(facts, rules=RULES):
    """Facts are (source identity, proposition, evidence) tuples.

    The caller supplies reviewed facts. Finite diagnostic output must use
    its own source identity and must not be renamed to the original source.
    """
    grouped = defaultdict(dict)
    for source, proposition, evidence in facts:
        if not source or not proposition or not evidence:
            raise ValueError("source, proposition and evidence are required")
        grouped[source][proposition] = {"given":evidence}
    for propositions in grouped.values():
        changed = True
        while changed:
            changed = False
            for name, requires, output in rules:
                if output not in propositions and all(k in propositions for k in requires):
                    propositions[output] = {"rule":name,"requires":list(requires)}
                    changed = True
    return dict(grouped)


def route(goal, available, rules=RULES, visiting=frozenset()):
    """Generate inclusion-minimal sets of outstanding obligations for a goal."""
    if goal in available:
        return [frozenset()]
    if goal in visiting:
        return []
    producers = [(name,requires) for name,requires,output in rules if goal==output]
    if not producers:
        return [frozenset({goal})]
    next_visited = visiting | {goal}
    alternatives = set()
    for _,requires in producers:
        combined = {frozenset()}
        for premise in requires:
            sub = route(premise,available,rules,next_visited)
            combined = {a|b for a in combined for b in sub}
        alternatives.update(combined)
    return sorted((a for a in alternatives if not any(b<a for b in alternatives)),
                  key=lambda a:(len(a),sorted(a)))


def frontier():
    # The local shape result is available as mathematics. The physical
    # original source has no reviewed bridge facts in this campaign.
    established = derive([("standard_A3_frame","a3_standard_weight_frame",
                           "finding:session-a3-null-gram-20261113")])
    alternatives = [sorted(row) for row in route("full_tfpt_candidate",{})]
    blockers = [sorted(row) for row in route("common_chiral_and_regional_origin",{})]
    frequency = defaultdict(int)
    for row in blockers:
        for item in row:
            frequency[item] += 1
    return {"typed_derivations":established,"full_candidate_routes":alternatives,
            "common_origin_routes":blockers,
            "four_mark_route":[sorted(row) for row in route("physical_four_mark_bridge",{},FOUR_MARK_RULES)],
            "ordered_return_route":[sorted(row) for row in route(
                "original_ordered_dynamical_sewing",{},ORDERED_RETURN_RULES)],
            "shared_source_obligations":sorted(k for k,v in frequency.items() if v==len(blockers)),
            "physical_gate_passes":[]}
