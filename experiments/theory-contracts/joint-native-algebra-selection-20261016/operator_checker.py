#!/usr/bin/env python3
"""Small exact audit of the native 4/5/20 operator chain.

This is deliberately not a new group census.  It imports the existing native
constructors, checks the corrected symmetric-power statement, reconstructs the
actual erasure isometry, and tests only the marked operators and six F20 lifts
needed by operator_audit.md.
"""

from __future__ import annotations

import hashlib
import itertools as it
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp


HERE = Path(__file__).resolve().parent
OUT = HERE / "operator_results.json"
REPO = Path("/Users/stefanhamann/Projekte/tfpt-theoryv4")
OLD = Path("/Users/stefanhamann/Documents/Codex/2026-10-07/che-2/work")
sys.path.insert(0, str(REPO))

from tfpt_explorer.cartan_source import _lattice_data, _rays  # noqa: E402
from tfpt_explorer.charge_selection import _parse_matrix, _process_frame  # noqa: E402
from tfpt_explorer.marker_selection import build_marker_data  # noqa: E402
from tfpt_explorer.process import _code_basis, _invariants  # noqa: E402
from tfpt_explorer.quartic_selection import _pauli_group, _symmetric_characters  # noqa: E402


SOURCE_FILES = [
    REPO / "tfpt_explorer/process.py",
    REPO / "tfpt_explorer/quartic_selection.py",
    REPO / "tfpt_explorer/marker_selection.py",
    REPO / "verification/v634_st31_structure.py",
    REPO / "verification/v689_gaussian_code_bridge.py",
    REPO / "verification/v419_seam_galois_carrier.py",
    OLD / "common_reflection_source_20261008/hyper_common_reflection_check.py",
    OLD / "common_reflection_source_20261008/hyper_mark_compatible_search.py",
    OLD / "common_reflection_source_20261008/audit_source_first.py",
    OLD / "common_reflection_source_20261008/hidden_sgn_check.py",
    OLD / "joint_synthesis/marked_rr_completion.py",
]


def clean(matrix: sp.Matrix) -> sp.Matrix:
    return matrix.applyfunc(sp.simplify)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def native_reflection(root: tuple[int, ...]) -> sp.Matrix:
    vector = sp.Matrix(
        [sp.Rational(root[2 * k], 2) + sp.I * sp.Rational(root[2 * k + 1], 2) for k in range(4)]
    )
    return clean(sp.eye(4) - 2 * vector * vector.conjugate().T / (vector.conjugate().T * vector)[0])


def quartic_compression(matrix: sp.Matrix, supports: list[list[tuple[int, ...]]]) -> sp.Matrix:
    out = sp.zeros(5)
    for row, targets in enumerate(supports):
        for column, sources in enumerate(supports):
            value = sum(
                sp.prod(matrix[a, b] for a, b in zip(target, source))
                for target in targets
                for source in sources
            )
            out[row, column] = sp.simplify(value / sp.sqrt(len(targets) * len(sources)))
    return clean(out)


def tensor_power(matrix: sp.Matrix, count: int) -> sp.Matrix:
    out = sp.ones(1, 1)
    for _ in range(count):
        out = sp.kronecker_product(out, matrix)
    return out


def erasure_isometry(supports: list[list[tuple[int, ...]]], erased: int = 0) -> sp.Matrix:
    out = sp.zeros(64, 20)
    for logical, support in enumerate(supports):
        amplitude = sp.simplify(2 / sp.sqrt(len(support)))
        for word in support:
            gauge = word[erased]
            remaining = tuple(word[k] for k in range(4) if k != erased)
            row = remaining[0] * 16 + remaining[1] * 4 + remaining[2]
            out[row, 4 * logical + gauge] += amplitude
    return out


def permutation_matrix(permutation: tuple[int, ...]) -> sp.Matrix:
    out = sp.zeros(len(permutation))
    for column, row in enumerate(permutation):
        out[row, column] = 1
    return out


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[right[i]] for i in range(len(left)))


def inverse(permutation: tuple[int, ...]) -> tuple[int, ...]:
    out = [0] * len(permutation)
    for i, value in enumerate(permutation):
        out[value] = i
    return tuple(out)


def matrix_record(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(sp.simplify(matrix[i, j])) for j in range(matrix.cols)] for i in range(matrix.rows)]


def main() -> None:
    paulis = _pauli_group()
    assert len(paulis) == len({tuple(matrix) for matrix in paulis}) == 64

    # Correct typing: full tensor invariants and symmetric quartic invariants are
    # different representations.  The latter is the polynomial/code carrier.
    full_tensor_fixed = sp.simplify(sum(sp.trace(g) ** 4 for g in paulis) / 64)
    sym_coefficients = [
        sp.simplify(sum(_symmetric_characters(g, 4)[degree] for g in paulis) / 64)
        for degree in range(5)
    ]
    assert full_tensor_fixed == 16
    assert sym_coefficients == [1, 0, 0, 0, 5]

    _, supports = _code_basis()
    # Each source support is an S4 orbit, so the five sparse columns lie in Sym^4.
    symmetric_supports = all(
        set(support) == {tuple(word[i] for i in p) for word in support for p in it.permutations(range(4))}
        for support in supports
    )
    assert symmetric_supports

    # Exact monomial action proves that the actual five sparse code columns are
    # pointwise P64-fixed, without forming 64 dense 256 by 256 matrices.
    pauli_code_fixed = []
    for g in paulis:
        target = []
        phase = []
        for source in range(4):
            nz = [(row, sp.simplify(g[row, source])) for row in range(4) if g[row, source] != 0]
            assert len(nz) == 1
            target.append(nz[0][0])
            phase.append(nz[0][1])
        fixed = True
        for support in supports:
            support_set = set(support)
            for word in support:
                image = tuple(target[index] for index in word)
                coefficient = sp.prod(phase[index] for index in word)
                fixed &= image in support_set and sp.simplify(coefficient - 1) == 0
        pauli_code_fixed.append(bool(fixed))
    assert all(pauli_code_fixed)

    invariants = _invariants()
    roots = _rays(_lattice_data()["roots"])
    chosen_ids = [1, 21, 32]
    chosen4 = [native_reflection(roots[index]) for index in chosen_ids]
    chosen5 = [quartic_compression(matrix, supports) for matrix in chosen4]
    assert [list(invariants["transposition_pairs"][invariants["logical_to_pair"][invariants["labels"][i]]]) for i in chosen_ids] == [[2, 3], [2, 4], [3, 4]]

    W = erasure_isometry(supports)
    assert clean(W.T * W - sp.eye(20)) == sp.zeros(20)
    erasure_intertwiners = []
    for event4, event5 in zip(chosen4, chosen5):
        residual = clean(
            tensor_power(event4, 3) * W
            - W * sp.kronecker_product(event5, event4.conjugate())
        )
        erasure_intertwiners.append(residual == sp.zeros(64, 20))
    assert erasure_intertwiners == [True, True, True]

    # The old S3 sign twin, now in the actual marked bases.
    s4 = chosen4[0]
    c4 = clean(chosen4[0] * chosen4[1])
    s5 = chosen5[0]
    c5 = clean(chosen5[0] * chosen5[1])
    rep4 = [sp.eye(4), c4, c4**2, s4, clean(s4 * c4), clean(s4 * c4**2)]
    rep5 = [sp.eye(5), c5, c5**2, s5, clean(s5 * c5), clean(s5 * c5**2)]
    signs = [1, 1, 1, -1, -1, -1]
    e_sgn = clean(
        sum(
            (
                sign * sp.kronecker_product(event5, event4.conjugate())
                for sign, event5, event4 in zip(signs, rep5, rep4)
            ),
            sp.zeros(20),
        )
        / 6
    )
    z20 = sp.eye(20) - 2 * e_sgn
    assert e_sgn.rank() == 1 and clean(e_sgn**2 - e_sgn) == sp.zeros(20)
    assert clean(z20**2 - sp.eye(20)) == sp.zeros(20)

    marker = build_marker_data()["data"]
    rotation = _parse_matrix(marker["basis_alignment"]["rotation_matrix"])
    native_y = _parse_matrix(marker["native_polar_transport"]["native_charge_matrix"])
    y = clean(rotation.T * native_y * rotation)
    frame = _process_frame()

    def standard(permutation: tuple[int, ...]) -> sp.Matrix:
        return clean(frame.T * permutation_matrix(permutation) * frame)

    sigma5 = standard((0, 1, 3, 4, 2, 5))
    tau5 = standard((1, 0, 2, 3, 4, 5))
    p2 = clean((6 * y + 2 * sp.eye(5)) / 5)
    p3 = sp.eye(5) - p2
    f_sigma = clean((sp.eye(5) + sigma5 + sigma5**2) / 3)
    marked_R = clean(p3 * (sp.eye(5) - 2 * f_sigma) + sp.I * p2 * tau5)
    marked_R20 = sp.kronecker_product(marked_R, sp.eye(4))
    z_commutes_marked_R = clean(z20 * marked_R20 - marked_R20 * z20) == sp.zeros(20)
    assert z_commutes_marked_R
    marked_stabilizer_permutations = [
        (1, 0, 2, 3, 4, 5),  # weak S2 generator (01)
        (0, 1, 3, 2, 4, 5),  # color S3 generators (23),(24)
        (0, 1, 4, 3, 2, 5),
    ]
    marked_stabilizer_matrices = [standard(p) for p in marked_stabilizer_permutations]
    r_stabilizer_commutators_exact = [
        clean(marked_R * matrix - matrix * marked_R) == sp.zeros(5)
        for matrix in marked_stabilizer_matrices
    ]
    assert r_stabilizer_commutators_exact == [True, True, True]

    # E(rho)=rho tensor I4/4.  A rank-one color-standard state gives an exact
    # branchwise failure for z although all six genuine S3 branches are covariant.
    pstd5 = clean((2 * sp.eye(5) - c5 - c5**2) / 3)
    rho5 = clean(pstd5 * (sp.eye(5) + s5) / 2)
    assert rho5.rank() == 1 and sp.trace(rho5) == 1 and clean(rho5**2 - rho5) == sp.zeros(5)
    erased = sp.kronecker_product(rho5, sp.eye(4) / 4)
    z_difference = clean(z20 * erased * z20 - erased)
    z_branch_fixes_erasure = z_difference == sp.zeros(20)
    z_difference_hs2 = sp.simplify(sp.trace(z_difference.conjugate().T * z_difference))
    assert not z_branch_fixes_erasure and z_difference.rank() == 2 and z_difference_hs2 == sp.Rational(1, 8)
    # Hermitian, traceless and rank two imply eigenvalues +/-1/4 here.
    assert clean(z_difference - z_difference.conjugate().T) == sp.zeros(20)
    assert sp.trace(z_difference) == 0
    z_trace_distance = sp.sqrt(z_difference_hs2 / 2)
    assert z_trace_distance == sp.Rational(1, 4)
    genuine_branch_covariance = []
    for event5, event4 in zip(rep5, rep4):
        joint = sp.kronecker_product(event5, event4.conjugate())
        left = clean(joint * erased * joint.conjugate().T)
        right = sp.kronecker_product(clean(event5 * rho5 * event5.conjugate().T), sp.eye(4) / 4)
        genuine_branch_covariance.append(left == right)
    assert genuine_branch_covariance == [True] * 6

    # R is a coherent marked completion.  It has the same spectrum as a
    # 4-cycle in the standard S6 module, but is not any of the 720 actual S6
    # matrices in the fixed process frame.
    permutations6 = list(it.permutations(range(6)))
    stabilizer_centralizer_permutations = [
        p for p in permutations6
        if all(compose(p, h) == compose(h, p) for h in marked_stabilizer_permutations)
    ]
    assert stabilizer_centralizer_permutations == [
        (0, 1, 2, 3, 4, 5),
        (1, 0, 2, 3, 4, 5),
    ]
    # The matrices themselves are exact algebraic numbers.  Convert only this
    # finite 720-way equality/commutator census to complex128 to avoid spending
    # minutes repeatedly simplifying the same radicals; the minimum nonzero
    # residual is recorded below and is macroscopically separated from zero.
    standard_s6 = [standard(p) for p in permutations6]
    marked_R_np = np.array(marked_R, dtype=complex)
    s6_np = [np.array(matrix, dtype=complex) for matrix in standard_s6]
    equality_residuals = [float(np.linalg.norm(matrix - marked_R_np)) for matrix in s6_np]
    commutator_residuals = [float(np.linalg.norm(matrix @ marked_R_np - marked_R_np @ matrix)) for matrix in s6_np]
    marked_R_matches = [p for p, residual in zip(permutations6, equality_residuals) if residual < 1e-12]
    assert not marked_R_matches
    r_centralizer = [
        p for p, residual in zip(permutations6, commutator_residuals)
        if residual < 1e-12
    ]
    assert len(r_centralizer) == 12
    assert clean(sigma5**3 - sp.eye(5)) == sp.zeros(5)
    assert clean(sigma5 - marked_R) != sp.zeros(5)

    # Canonical 3+2 orbit already present in the original process: every perfect
    # matching is a fixed-point-free involution tau in S6.  Its standard trace is
    # -1, hence multiplicities (+1,-1)=(2,3) and (I+tau)/2 has rank two.
    matchings = invariants["matchings"]
    response_ranks = []
    physical_p2_equalities = []
    for matching in matchings:
        permutation = list(range(6))
        for left, right in matching:
            permutation[left], permutation[right] = permutation[right], permutation[left]
        tau = standard(tuple(permutation))
        response = clean((sp.eye(5) + tau) / 2)
        response_ranks.append(response.rank())
        physical_p2_equalities.append(clean(response - p2) == sp.zeros(5))
    assert len(matchings) == 15 and response_ranks == [2] * 15
    assert physical_p2_equalities == [False] * 15
    color_generators = marked_stabilizer_permutations[1:]
    p2_color_commutators = [
        clean(p2 * standard(generator) - standard(generator) * p2) == sp.zeros(5)
        for generator in color_generators
    ]
    assert p2_color_commutators == [True, True]
    color_centralizer = [
        permutation for permutation in permutations6
        if all(compose(permutation, generator) == compose(generator, permutation)
               for generator in color_generators)
    ]
    color_centralizer_fixed_point_free_involutions = [
        permutation for permutation in color_centralizer
        if permutation != tuple(range(6))
        and compose(permutation, permutation) == tuple(range(6))
        and all(permutation[index] != index for index in range(6))
    ]
    assert len(color_centralizer) == 6
    assert color_centralizer_fixed_point_free_involutions == []
    # A double transposition also makes a 3+2 split, but with the complementary
    # orientation: (I+tau)/2 has rank three.  It is not the response constructor.
    double = standard((1, 0, 3, 2, 4, 5))
    double_positive_rank = clean((sp.eye(5) + double) / 2).rank()
    assert double_positive_rank == 3

    # Six honest native F20 lifts.  The fixed marked label is 5; S5 on the other
    # labels contains six Sylow-5 subgroups.  A BFS over the 15 already lifted
    # transpositions supplies one native lift of each quotient permutation; the
    # listed Pauli corrections make all three F20 relations exact.
    identity_perm = tuple(range(6))
    pairs = list(invariants["transposition_pairs"])
    pair_lifts: dict[tuple[int, int], sp.Matrix] = {}
    pair_ray_ids: dict[str, int] = {}
    for pair in pairs:
        pair_id = pairs.index(pair)
        ray_id = next(
            index
            for index, label in enumerate(invariants["labels"])
            if invariants["logical_to_pair"][label] == pair_id
        )
        pair_lifts[pair] = native_reflection(roots[ray_id])
        pair_ray_ids[str(pair)] = ray_id
    quotient_lifts = {identity_perm: sp.eye(4)}
    frontier = [identity_perm]
    while frontier:
        current = frontier.pop(0)
        for pair, lift in pair_lifts.items():
            transposition = list(identity_perm)
            left, right = pair
            transposition[left], transposition[right] = transposition[right], transposition[left]
            image = compose(tuple(transposition), current)
            if image not in quotient_lifts:
                quotient_lifts[image] = clean(lift * quotient_lifts[current])
                frontier.append(image)
    assert len(quotient_lifts) == 720

    sylow5: list[tuple[frozenset[tuple[int, ...]], tuple[int, ...]]] = []
    for p5 in it.permutations(range(5)):
        permutation = p5 + (5,)
        subgroup = []
        power = identity_perm
        for _ in range(5):
            subgroup.append(power)
            power = compose(permutation, power)
        if power != identity_perm or any(permutation[i] == i for i in range(5)):
            continue
        frozen = frozenset(subgroup)
        if frozen not in [old[0] for old in sylow5]:
            sylow5.append((frozen, min(subgroup[1:])))
    assert len(sylow5) == 6

    def permutation_subgroup(generators: list[tuple[int, ...]]) -> set[tuple[int, ...]]:
        subgroup = {identity_perm}
        frontier = [identity_perm]
        while frontier:
            current = frontier.pop()
            for generator in generators:
                image = compose(generator, current)
                if image not in subgroup:
                    subgroup.add(image)
                    frontier.append(image)
        return subgroup

    weak_swap = (1, 0, 2, 3, 4, 5)
    color_cycle = (0, 1, 3, 4, 2, 5)
    color_swap = (0, 1, 3, 2, 4, 5)
    marked_H0 = permutation_subgroup([weak_swap, color_cycle])
    marked_H = permutation_subgroup([weak_swap, color_cycle, color_swap])
    assert len(marked_H0) == 6 and len(marked_H) == 12
    sylow_index = {subgroup: index for index, (subgroup, _) in enumerate(sylow5)}

    def sylow_action(group: set[tuple[int, ...]]) -> dict[int, list[int]]:
        action = {}
        for h_index, h in enumerate(sorted(group)):
            hi = inverse(h)
            image = frozenset(compose(compose(h, element), hi) for element in sylow5[0][0])
            action[h_index] = [sylow_index[image]]
        return action

    H0_images = [row[0] for row in sylow_action(marked_H0).values()]
    H_images = [row[0] for row in sylow_action(marked_H).values()]
    H0_stabilizer_order = H0_images.count(0)
    H_stabilizer_order = H_images.count(0)
    assert sorted(H0_images) == list(range(6)) and H0_stabilizer_order == 1
    assert sorted(set(H_images)) == list(range(6)) and H_stabilizer_order == 2
    pauli_corrections = [(3, 7), (0, 3), (1, 14), (1, 14), (2, 4), (0, 5)]
    # Standalone v419 carrier, reconstructed verbatim: multiplication by zeta5
    # in the power basis and Frobenius zeta5 -> zeta5^2.
    x = sp.symbols("x")
    phi5 = x**4 + x**3 + x**2 + x + 1
    cyclotomic_C5 = sp.Matrix([[0, 0, 0, -1], [1, 0, 0, -1], [0, 1, 0, -1], [0, 0, 1, -1]])

    def cyclotomic_column(exponent: int) -> sp.Matrix:
        polynomial = sp.Poly(sp.rem(x**exponent, phi5, x), x)
        return sp.Matrix([polynomial.coeff_monomial(x**degree) for degree in range(4)])

    cyclotomic_G4 = sp.Matrix.hstack(*(cyclotomic_column(2 * degree) for degree in range(4)))
    assert cyclotomic_G4**4 == sp.eye(4)
    assert clean(cyclotomic_G4 * cyclotomic_C5 * cyclotomic_G4.inv() - cyclotomic_C5**2) == sp.zeros(4)
    f20_rows = []
    for index, ((_, c_perm), (c_pauli, g_pauli)) in enumerate(zip(sylow5, pauli_corrections)):
        g_candidates = []
        for g5 in it.permutations(range(5)):
            g_perm = g5 + (5,)
            if compose(compose(g_perm, c_perm), inverse(g_perm)) == compose(c_perm, c_perm):
                g_candidates.append(g_perm)
        g_perm = min(g_candidates)
        C = clean(paulis[c_pauli] * quotient_lifts[c_perm])
        G = clean(paulis[g_pauli] * quotient_lifts[g_perm])
        checks = {
            "C5": clean(C**5 - sp.eye(4)) == sp.zeros(4),
            "G4": clean(G**4 - sp.eye(4)) == sp.zeros(4),
            "frobenius": clean(G * C * G.inv() - C**2) == sp.zeros(4),
        }
        assert all(checks.values())
        # Solve the actual representation intertwiner, with no chosen embedding.
        variables = sp.symbols(f"a{index}_0:16")
        A_variable = sp.Matrix(4, 4, variables)
        equations = list(C * A_variable - A_variable * cyclotomic_C5)
        equations += list(G * A_variable - A_variable * cyclotomic_G4)
        coefficient_matrix, _ = sp.linear_eq_to_matrix(equations, variables)
        nullspace = coefficient_matrix.nullspace()
        assert len(nullspace) == 1
        A = sp.Matrix(4, 4, nullspace[0])
        first_nonzero = next(entry for entry in A if entry != 0)
        A = clean(A / first_nonzero)
        intertwiner_checks = {
            "nullity": len(nullspace),
            "rank": A.rank(),
            "C_intertwines": clean(C * A - A * cyclotomic_C5) == sp.zeros(4),
            "G_intertwines": clean(G * A - A * cyclotomic_G4) == sp.zeros(4),
        }
        assert intertwiner_checks == {
            "nullity": 1,
            "rank": 4,
            "C_intertwines": True,
            "G_intertwines": True,
        }
        f20_rows.append(
            {
                "sylow_index": index,
                "C5_quotient_permutation": list(c_perm),
                "G4_quotient_permutation": list(g_perm),
                "pauli_indices": [c_pauli, g_pauli],
                "checks": checks,
                "cyclotomic_intertwiner": matrix_record(A),
                "cyclotomic_intertwiner_checks": intertwiner_checks,
                "C": matrix_record(C),
                "G": matrix_record(G),
            }
        )

    result = {
        "status": "PASS_EXACT_NATIVE_OPERATOR_AUDIT",
        "source_sha256": {str(path): sha256(path) for path in SOURCE_FILES},
        "corrected_pauli_invariants": {
            "P64_order": len(paulis),
            "full_tensor_fourth_power_fixed_dimension": int(full_tensor_fixed),
            "symmetric_degrees_0_to_4": [int(value) for value in sym_coefficients],
            "Sym4_fixed_dimension": int(sym_coefficients[4]),
            "actual_sparse_code_supports_S4_symmetric": symmetric_supports,
            "actual_sparse_code_pointwise_P64_fixed": all(pauli_code_fixed),
        },
        "erasure": {
            "chosen_native_ray_ids": chosen_ids,
            "logical_pairs": [[2, 3], [2, 4], [3, 4]],
            "W_shape": list(W.shape),
            "W_isometry_exact": True,
            "eventwise_intertwiners": erasure_intertwiners,
            "E_formula": "E(rho5)=rho5 tensor I4/4; physical retained state is W E(rho5) W^*",
        },
        "old_S3_sign_twin": {
            "e_sgn_rank_on_H20": e_sgn.rank(),
            "z20_involution": True,
            "commutes_with_actual_marked_R_tensor_I4": z_commutes_marked_R,
            "genuine_S3_branches_erasure_covariant": genuine_branch_covariance,
            "z_branch_fixes_erasure_witness": z_branch_fixes_erasure,
            "witness_rho_rank": rho5.rank(),
            "z_erasure_difference_rank": z_difference.rank(),
            "z_erasure_difference_HS_squared": str(z_difference_hs2),
            "z_erasure_trace_distance": str(z_trace_distance),
        },
        "marked_R": {
            "formula": "P3(I-2F_sigma)+i P2 tau",
            "characteristic_polynomial": str(sp.factor(marked_R.charpoly().as_expr())),
            "trace": str(sp.simplify(sp.expand_complex(sp.trace(marked_R)))),
            "matches_among_720_actual_S6_matrices": len(marked_R_matches),
            "minimum_Frobenius_distance_to_actual_S6": min(equality_residuals),
            "centralizer_in_actual_S6_size": len(r_centralizer),
            "minimum_nonzero_commutator_Frobenius_norm": min(
                residual for residual in commutator_residuals if residual >= 1e-12
            ),
            "compiler_degree4_image": "sigma5",
            "compiler_degree4_image_order": 3,
            "marked_R_order": 4,
            "compiler_image_equals_marked_R": False,
            "commutes_exactly_with_marked_S2_times_S3_generators": r_stabilizer_commutators_exact,
            "centralizer_of_marked_S2_times_S3_in_S6_order": len(stabilizer_centralizer_permutations),
            "exact_nonmembership_proof": "If R were rho5(p), p would lie in Cent_S6(S2_(01) x S3_(234))={<identity>,(01)}, whose elements have order at most 2; R has order 4.",
        },
        "canonical_3_plus_2": {
            "perfect_matching_count": len(matchings),
            "response_projector_ranks": response_ranks,
            "perfect_matching_involution_spectrum": "+1^2,-1^3",
            "physical_P2_rank": p2.rank(),
            "physical_P2_equals_each_response_projector": physical_p2_equalities,
            "physical_P2_in_perfect_matching_orbit": any(physical_p2_equalities),
            "physical_P2_commutes_with_color_S3_generators": p2_color_commutators,
            "centralizer_of_color_S3_in_S6_order": len(color_centralizer),
            "fixed_point_free_involutions_in_color_S3_centralizer": len(
                color_centralizer_fixed_point_free_involutions
            ),
            "exact_nonidentification_proof": "If P2=(I+rho5(tau))/2 for a perfect-matching involution tau, then rho5(tau)=2P2-I would commute with color S3. Faithfulness forces tau into Cent_S6(S3_(234))=S3_(015), which has no fixed-point-free involution.",
            "double_transposition_positive_projector_rank": double_positive_rank,
            "selection_boundary": "the original response constructor supplies a canonical 15-orbit of 2+3 splits, but the actual physical P2 mark is not a member; equal dimensions are not an identification",
        },
        "native_F20": {
            "marked_fixed_label": 5,
            "sylow5_subgroups_on_remaining_five_labels": len(sylow5),
            "all_six_have_exact_honest_F20_lift": True,
            "all_six_intertwine_with_v419_carrier_uniquely_up_to_scalar": True,
            "v419_C5": matrix_record(cyclotomic_C5),
            "v419_G4": matrix_record(cyclotomic_G4),
            "rows": f20_rows,
            "marked_stabilizer_orbit": {
                "H0": "S2_(01) x C3_(234)",
                "H0_order": len(marked_H0),
                "H0_images_of_first_copy": H0_images,
                "H0_stabilizer_order": H0_stabilizer_order,
                "H0_action_regular_on_six": True,
                "H": "S2_(01) x S3_(234)",
                "H_order": len(marked_H),
                "H_images_of_first_copy": H_images,
                "H_stabilizer_order": H_stabilizer_order,
                "H_action_transitive_on_six": True,
            },
            "selection_boundary": "the six literal Sylow-5 representatives form one orbit under the existing marked stabilizer (already regular under H0); modulo that native symmetry there is one equivalence class, although no basis representative is singled out",
        },
        "scope": (
            "Finite exact checks of existing constructors only.  This does not derive the marked R, "
            "choose one F20 embedding, identify an internal module with spacetime, or supply dynamics."
        ),
    }
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": result["status"],
        "corrected_pauli_invariants": result["corrected_pauli_invariants"],
        "erasure": result["erasure"],
        "old_S3_sign_twin": result["old_S3_sign_twin"],
        "marked_R": result["marked_R"],
        "canonical_3_plus_2": result["canonical_3_plus_2"],
        "native_F20_count": len(f20_rows),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
