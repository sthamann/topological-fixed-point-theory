"""Unpromoted exact quotient-switch model in the declared Hilbert space C[G48].

Four conjugate cubic decks give overlapping invariant subspaces, not four
orthogonal physical sectors.  Their A4 permutation action has no quarter cycle.
The conditional compression is accompanied by a separately declared complete
quantum channel: choose a deck uniformly and perform its unread binary
measurement.  Poisson events give a completely positive continuous semigroup.
The state space, measurement rule and event parameter are additional choices;
no native R113 kernel, physical clock or F1 source is derived here.
"""

from collections import Counter
import json

import sympy as sp

if __package__:
    from .origin_picard import ExactPicardGroup, require
else:
    from origin_picard import ExactPicardGroup, require


def _cycle_lengths(permutation):
    seen, lengths = set(), []
    for start in range(len(permutation)):
        if start in seen:
            continue
        current, length = start, 0
        while current not in seen:
            seen.add(current)
            length += 1
            current = permutation[current]
        lengths.append(length)
    return tuple(sorted(lengths))


def _complete_comparison_dynamics(group, projectors, subgroup_action):
    """Construct a complete channel and certify its Poisson semigroup.

    Phi(X) = (1/4) sum_a [Pa X Pa + (I-Pa) X (I-Pa)].
    Its eight Hermitian Kraus operators are Pa/2 and (I-Pa)/2.
    L = Phi-id is a GKSL generator because sum_j K_j^dagger K_j = I:
      L(X) = sum_j [K_j X K_j^dagger - {K_j^dagger K_j,X}/2].
    For s >= 0, T_s = exp(s L) = exp(-s) sum_n s^n Phi^n/n!.
    Every summand is CPTP and the weights are nonnegative and sum to one.
    Norm convergence is automatic in finite dimension.  This gives a complete
    semigroup construction without diagonalizing a 2304-by-2304 superoperator.
    s is a declared dimensionless event parameter, not native physical time.
    """
    dimension = len(group.elements)
    identity, zero = sp.eye(dimension), sp.zeros(dimension)
    require(len(projectors) == 4, "The uniform deck-choice rule requires four decks")
    kraus = []
    for projector in projectors:
        complement = identity - projector
        require(projector.T == projector and complement.T == complement and
                projector * projector == projector and
                complement * complement == complement and
                projector * complement == zero,
                "A binary pinching is not a self-adjoint orthogonal pinching")
        kraus.extend((projector / 2, complement / 2))
    require(len(kraus) == 8 and all(k.H == k for k in kraus),
            "The eight Hermitian Kraus operators changed")
    require(sum((k.H * k for k in kraus), zero) == identity,
            "The complete pinching channel is not trace preserving")
    require(sum((k * k.H for k in kraus), zero) == identity,
            "The complete pinching channel is not unital")

    def phi(operator):
        require(operator.shape == identity.shape, "Operator has the wrong state-space dimension")
        return sum((k * operator * k.H for k in kraus), sp.zeros(dimension))

    # Each group generator permutes the projectors.  Since the listed
    # generators generate G48, this establishes covariance for every g in G48.
    covariance_checks = 0
    for generator in group.generators:
        unitary = group.regular(generator)
        permutation = subgroup_action(generator)
        for index, projector in enumerate(projectors):
            require(unitary * projector * unitary.T == projectors[permutation[index]],
                    "The pinching family is not covariant under the full group")
            covariance_checks += 1
    require(covariance_checks == 12, "The covariance generator checks were incomplete")

    state = projectors[0] / 16
    image = phi(state)
    initial_purity = sp.trace(state * state)
    final_purity = sp.trace(image * image)
    distance_squared = sp.trace((image - state).H * (image - state))
    require(sp.trace(state) == sp.trace(image) == 1,
            "The test state lost trace under the complete channel")
    require(initial_purity == sp.Rational(1, 16) and
            final_purity == sp.Rational(55, 1296) and
            initial_purity - final_purity == sp.Rational(13, 648) and
            distance_squared == sp.Rational(5, 648),
            "The complete channel's exact nontrivial action changed")

    central = group.regular(group.rho)
    require(phi(identity) == identity and phi(central) == central,
            "The unital channel failed to preserve a central observable")
    central_involution = central * central
    require(central_involution.H == central_involution and
            central_involution * central_involution == identity,
            "The central involution is not Hermitian unitary")
    fixed_states = []
    for sign in (1, -1):
        support = (identity + sign * central_involution) / 2
        fixed_state = support / 24
        require(support * support == support and support.H == support and
                support.rank() == 24,
                "The stationary-state support is not a positive rank-24 projector")
        require(sp.trace(fixed_state) == 1 and phi(fixed_state) == fixed_state and
                sp.trace(fixed_state * fixed_state) == sp.Rational(1, 24),
                "The two explicit stationary density matrices changed")
        fixed_states.append(fixed_state)
    require(fixed_states[0] != fixed_states[1],
            "Non-primitivity requires distinct stationary states")

    event_parameter = sp.symbols("s", positive=True)
    number = sp.symbols("n", integer=True, nonnegative=True)
    poisson_total = sp.summation(
        sp.exp(-event_parameter) * event_parameter**number / sp.factorial(number),
        (number, 0, sp.oo),
    )
    require(sp.simplify(poisson_total - 1) == 0,
            "The Poisson weights do not normalize the complete channel mixture")

    return {
        "state_space": "End(C[G48])",
        "declared_rule": "Uniform deck choice followed by an unread binary measurement",
        "channel": "Phi(X)=(1/4)*sum_a[Pa*X*Pa+(I-Pa)*X*(I-Pa)]",
        "kraus_operators": "Pa/2 and (I-Pa)/2, for a=0,1,2,3",
        "kraus_count": len(kraus),
        "sum_Kdagger_K": "I48",
        "sum_K_Kdagger": "I48",
        "completely_positive": True,
        "trace_preserving": True,
        "unital": True,
        "each_pinching_Hilbert_Schmidt_self_adjoint_and_idempotent": True,
        "channel_Hilbert_Schmidt_self_adjoint": True,
        "G48_covariant": True,
        "projector_covariance_checks_on_generators": covariance_checks,
        "test_state": "P0/16",
        "test_state_output_trace": "1",
        "test_state_input_purity": str(initial_purity),
        "test_state_output_purity": str(final_purity),
        "test_state_purity_decrease": str(initial_purity - final_purity),
        "test_state_Hilbert_Schmidt_change_squared": str(distance_squared),
        "fixed_point_controls": {
            "central_Urho_is_fixed": True,
            "two_stationary_states": "(I+U_rho^2)/48 and (I-U_rho^2)/48",
            "stationary_support_ranks": [24, 24],
            "stationary_state_purities": ["1/24", "1/24"],
            "primitive": False,
            "complete_fixed_algebra_classified": False,
        },
        "continuous_semigroup": {
            "generator": "L=Phi-id",
            "GKSL_identity": "L(X)=sum_j[Kj*X*Kj^dagger-{Kj^dagger*Kj,X}/2]",
            "definition": "T_s=exp(s*(Phi-id))=exp(-s)*sum_n(s^n/n!)*Phi^n",
            "parameter_domain": "s>=0; declared dimensionless event parameter",
            "poisson_weight_sum": "1",
            "T_0_is_identity": True,
            "completely_positive_and_trace_preserving_for_all_s": True,
            "unital_for_all_s": True,
            "G48_covariant_for_all_s": True,
            "generator_is_log_of_discrete_channel": False,
            "closed_unitary_evolution_claimed": False,
            "native_physical_rate_derived": False,
            "native_KMS_clock_identified": False,
        },
        "scope": {
            "complete_open_comparison_dynamics_constructed": True,
            "additional_choices": [
                "Regular Hilbert space C[G48]",
                "Uniform distribution on four deck choices",
                "Unread binary pinching measurement",
                "Poisson event parameter with a chosen unit rate",
            ],
            "original_R113_or_native_M20_generator_identified": False,
            "physical_c3_normalization_derived": False,
            "F1_closed": False,
        },
    }


def controls():
    group = ExactPicardGroup()
    subgroups = group.cubic_deck_subgroups()
    subgroup_index = {subgroup: i for i, subgroup in enumerate(subgroups)}

    def subgroup_action(element):
        return tuple(subgroup_index[frozenset(group.conjugate(element, h)
                                             for h in subgroup)]
                     for subgroup in subgroups)

    permutations = {subgroup_action(g) for g in group.elements}
    cycle_census = Counter(_cycle_lengths(p) for p in permutations)
    require(len(permutations) == 12 and
            dict(cycle_census) == {(1, 1, 1, 1): 1, (2, 2): 3, (1, 3): 8},
            "The four-quotient permutation action changed")
    require(subgroup_action(group.rho) == (0, 1, 2, 3),
            "The central quarter action must fix every cubic deck")
    require(all(4 not in _cycle_lengths(p) for p in permutations),
            "A forbidden quarter cycle appeared in A4")

    identity = sp.eye(48)
    projectors = [sum((group.regular(h) for h in subgroup), sp.zeros(48)) / 3
                  for subgroup in subgroups]
    for projector in projectors:
        require(projector * projector == projector and projector.T == projector,
                "Subgroup averaging ceased to be an orthogonal projector")
        require(projector.rank() == 16, "A cubic-deck invariant subspace changed dimension")
    P = projectors[0]
    U = group.regular(group.switch)
    require(U.T * U == identity, "The group switch is not unitary in the regular representation")
    Q = U * P * U.T
    require(any(Q == projector for projector in projectors) and Q != P,
            "The additional automorphism did not change the cubic quotient")
    compressed = P * Q * P
    t = sp.symbols("t")
    compressed_expected = t**32 * (t - 1)**2 * (t - sp.Rational(1, 3))**8 * (
        t - sp.Rational(1, 9))**6
    require(sp.factor(compressed.charpoly(t).as_expr() - compressed_expected) == 0,
            "The exact quotient-switch spectrum changed")
    commutator_rank = (P * Q - Q * P).rank()
    require(commutator_rank == 28, "The two quotient readouts unexpectedly commute")
    survival = sp.trace(compressed) / 16
    require(survival == sp.Rational(1, 3), "Conditional mean survival changed")

    total = sum(projectors, sp.zeros(48))
    total_expected = t**12 * (t - 4)**2 * (t - 2)**16 * (t - sp.Rational(4, 3))**18
    require(sp.factor(total.charpoly(t).as_expr() - total_expected) == 0,
            "The sum of deck projections changed spectrum")
    require(total.rank() == 36 and sp.trace(total) == 64,
            "The overlapping readout supports were miscounted")
    require(total != identity and total * total != total,
            "The four projections were mistaken for a complete partition")
    pair_overlaps = [sp.trace(projectors[i] * projectors[j])
                     for i in range(4) for j in range(i + 1, 4)]
    require(all(value == sp.Rational(16, 3) for value in pair_overlaps),
            "The deck subspaces were mistaken for orthogonal sectors")

    return {
        "declared_state_space": "C[G48] with its regular unitary group action",
        "hilbert_dimension": 48,
        "native_TFPT_state_space_identification_derived": False,
        "quotient_permutations": {
            "image_order": len(permutations),
            "cycle_type_census": {str(k): v for k, v in sorted(cycle_census.items())},
            "central_rho_action": list(subgroup_action(group.rho)),
            "four_cycle_present": False,
            "native_quarter_clock_is_quotient_permutation": False,
        },
        "conditional_switch": {
            "projector_rank": 16,
            "switch_U_is_unitary": True,
            "commutator_rank": commutator_rank,
            "PQP_spectrum": {"1": 2, "1/3": 8, "1/9": 6, "0": 32},
            "maximally_mixed_P_sector_mean_survival": str(survival),
            "conditional_compression_is_complete_unitary_dynamics": False,
        },
        "four_readouts": {
            "sum_spectrum": {"4": 2, "2": 16, "4/3": 18, "0": 12},
            "sum_rank": 36,
            "sum_trace": 64,
            "pairwise_trace_overlap": [str(v) for v in pair_overlaps],
            "orthogonal_sectors": False,
            "sum_is_identity": False,
            "sum_is_projector": False,
            "four_projections_are_a_complete_POVM": False,
        },
        "complete_comparison_dynamics": _complete_comparison_dynamics(
            group, projectors, subgroup_action),
        "scope": {
            "status": "unpromoted exact declared-model control",
            "first_additional_choice": "Regular Hilbert space and a sequence of selected deck readouts",
            "complete_open_comparison_dynamics_constructed": True,
            "physical_history_weight_or_clock_selected": False,
            "native_R113_kernel_and_charge_intertwiner_derived": False,
            "physical_c3_normalization_derived": False,
            "F1_closed": False,
        },
    }


if __name__ == "__main__":
    print(json.dumps(controls(), indent=2))
