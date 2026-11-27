#!/usr/bin/env python3
"""Exact structure test for a four-terminal reciprocal response origin.

The starting object is not a TFPT seam, carrier, clock, or lattice.  It is the
Dirichlet-to-Neumann response of a passive reciprocal network on an unlabeled
finite terminal set.  Such a response is encoded by one effective conductance
for every unordered terminal pair.  We ask for the smallest terminal number
whose response *shape* (overall conductance scale removed) has more than one
nontrivial irreducible permutation sector.

The positive test is representation-theoretic, not a dimension lookup:

    R^{2-subsets of {1,2,3,4}} = 1 + 3 + 2.

We construct the three orthogonal equivariant projectors and check them against
all 24 permutations.  The negative controls are terminal counts 2, 3 and 5,
the natural integral response lattice, and the symmetry loss after a circular
order is fixed.

This is an exploration check.  It does not establish a physical TFPT origin.
"""

from __future__ import annotations

from itertools import combinations, permutations, product
from math import isqrt

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def edges(n: int) -> tuple[tuple[int, int], ...]:
    return tuple(combinations(range(n), 2))


def incidence(n: int) -> sp.Matrix:
    es = edges(n)
    return sp.Matrix(n, len(es), lambda i, j: int(i in es[j]))


def edge_permutation_matrix(p: tuple[int, ...]) -> sp.Matrix:
    es = edges(len(p))
    index = {e: i for i, e in enumerate(es)}
    out = sp.zeros(len(es))
    for source, (i, j) in enumerate(es):
        target_edge = tuple(sorted((p[i], p[j])))
        out[index[target_edge], source] = 1
    return out


def projectors_n4() -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Projectors to the S4 irreps 1, [31] and [22]."""
    b = incidence(4)
    q_vertex = b.T * (b * b.T).inv() * b
    one = sp.ones(6, 1)
    p1 = one * one.T / 6
    p3 = sp.simplify(q_vertex - p1)
    p2 = sp.eye(6) - q_vertex
    return p1, p3, p2


def check_projectors() -> dict[str, object]:
    p1, p3, p2 = projectors_n4()
    ps = (p1, p3, p2)
    ranks = tuple(p.rank() for p in ps)
    require(ranks == (1, 3, 2), f"wrong projector ranks: {ranks}")
    require(sum(ps, sp.zeros(6)) == sp.eye(6), "projectors do not resolve identity")
    for i, p in enumerate(ps):
        require(p * p == p, f"projector {i} is not idempotent")
        require(p.T == p, f"projector {i} is not orthogonal")
        for j, q in enumerate(ps):
            if i != j:
                require(p * q == sp.zeros(6), f"projectors {i},{j} overlap")

    commute_checks = 0
    for perm in permutations(range(4)):
        rho = edge_permutation_matrix(perm)
        for p in ps:
            require(rho * p == p * rho, f"projector not S4-equivariant for {perm}")
            commute_checks += 1

    # Character values on the five S4 conjugacy classes
    reps = {
        "1^4": (0, 1, 2, 3),
        "2,1,1": (1, 0, 2, 3),
        "2,2": (1, 0, 3, 2),
        "3,1": (1, 2, 0, 3),
        "4": (1, 2, 3, 0),
    }
    chars = tuple(int(sp.trace(edge_permutation_matrix(p))) for p in reps.values())
    require(chars == (6, 2, 2, 0, 0), f"unexpected two-subset character: {chars}")
    require(
        chars
        == tuple(
            a + b + c
            for a, b, c in zip(
                (1, 1, 1, 1, 1),
                (3, 1, -1, 0, -1),
                (2, 0, 2, -1, 0),
            )
        ),
        "character does not decompose as 1 + 3 + 2",
    )
    return {
        "ranks": ranks,
        "commute_checks": commute_checks,
        "characters": chars,
        "projectors": ps,
    }


def sector_dimensions(n: int) -> tuple[int, ...]:
    """Irrep dimensions of the permutation module on two-subsets.

    For n >= 4 this is [n] + [n-1,1] + [n-2,2].  The small cases are
    written separately because the last partition is absent.
    """
    if n == 2:
        return (1,)
    if n == 3:
        return (1, 2)
    return (1, n - 1, n * (n - 3) // 2)


def check_minimality_and_counterfamilies() -> dict[int, tuple[int, ...]]:
    table = {n: sector_dimensions(n) for n in range(2, 9)}
    for n, dims in table.items():
        require(sum(dims) == len(edges(n)), f"sector count fails for n={n}")
    first_multi = min(n for n, dims in table.items() if len(dims[1:]) >= 2)
    require(first_multi == 4, f"first multi-sector response is n={first_multi}, not 4")
    require(table[3][1:] == (2,), "three-terminal negative control failed")
    require(table[4][1:] == (3, 2), "four-terminal 3+2 failed")
    require(table[5][1:] == (4, 5), "five-terminal negative control failed")
    return table


def response_from_conductances(c: sp.Matrix) -> sp.Matrix:
    es = edges(4)
    require(c.shape == (6, 1), "need six effective conductances")
    lam = sp.zeros(4)
    for weight, (i, j) in zip(c, es):
        lam[i, i] += weight
        lam[j, j] += weight
        lam[i, j] -= weight
        lam[j, i] -= weight
    return lam


def check_passive_response() -> dict[str, object]:
    conductances = sp.Matrix([1, 2, 3, 5, 7, 11])
    lam = response_from_conductances(conductances)
    require(lam == lam.T, "reciprocity failed")
    require(lam * sp.ones(4, 1) == sp.zeros(4, 1), "gauge kernel failed")
    require(lam.rank() == 3, "connected response must have rank three")
    recovered = sp.Matrix([-lam[i, j] for i, j in edges(4)])
    require(recovered == conductances, "pair response is not recovered from Lambda")

    # Exact positivity on a finite spanning test set, plus a positive tree cofactor.
    energies = []
    for x_raw in product(range(-2, 3), repeat=4):
        if len(set(x_raw)) == 1:
            continue
        x = sp.Matrix(x_raw)
        energy = (x.T * lam * x)[0]
        require(energy > 0, f"passivity failed on {x_raw}: {energy}")
        energies.append(energy)
    cofactor = lam[:3, :3].det()
    require(cofactor > 0, "matrix-tree cofactor is not positive")

    p1, p3, p2 = projectors_n4()
    shape = conductances / sum(conductances)
    require(p1 * shape + p3 * shape + p2 * shape == shape, "shape split failed")
    require(p1 * shape == sp.ones(6, 1) / 6, "scale direction not isolated")
    require(p3 * shape != sp.zeros(6, 1), "sample missed the 3-sector")
    require(p2 * shape != sp.zeros(6, 1), "sample missed the 2-sector")
    return {
        "conductances": conductances,
        "response": lam,
        "tree_cofactor": cofactor,
        "min_test_energy": min(energies),
        "shape3": p3 * shape,
        "shape2": p2 * shape,
    }


def check_integral_lattice_countercontrol() -> dict[str, object]:
    """Show that the natural response lattice is A5, not D5."""
    # A5 = {z in Z^6 : sum z = 0}; columns are e_i-e_6.
    a5_basis = sp.Matrix.vstack(sp.eye(5), -sp.ones(1, 5))
    gram_a5 = a5_basis.T * a5_basis
    det_a5 = int(gram_a5.det())
    require(det_a5 == 6, f"natural response lattice is not A5: det={det_a5}")

    # Integral basis of the standard 3-sector: B^T applied to A3 roots.
    a3_roots = sp.Matrix(
        [
            [1, 0, 0],
            [-1, 1, 0],
            [0, -1, 1],
            [0, 0, -1],
        ]
    )
    v3_basis = incidence(4).T * a3_roots
    gram_v3 = v3_basis.T * v3_basis
    require(int(gram_v3.det()) == 32, "3-sector lattice determinant changed")

    # Integral kernel of incidence B: a scaled A2 lattice.
    v2_basis = sp.Matrix(
        [
            [1, 0],
            [0, 1],
            [-1, -1],
            [-1, -1],
            [0, 1],
            [1, 0],
        ]
    )
    require(incidence(4) * v2_basis == sp.zeros(4, 2), "2-sector is not ker B")
    gram_v2 = v2_basis.T * v2_basis
    require(int(gram_v2.det()) == 12, "2-sector lattice determinant changed")

    combined = v3_basis.row_join(v2_basis)
    gram_combined = combined.T * combined
    det_combined = int(gram_combined.det())
    require(det_combined == 384, f"combined sector lattice determinant={det_combined}")
    index_sq = det_combined // det_a5
    index = isqrt(index_sq)
    require(index * index == index_sq and index == 8, f"wrong gluing index: {index_sq}")

    d5_cartan = sp.Matrix(
        [
            [2, -1, 0, 0, 0],
            [-1, 2, -1, 0, 0],
            [0, -1, 2, -1, -1],
            [0, 0, -1, 2, 0],
            [0, 0, -1, 0, 2],
        ]
    )
    det_d5 = int(d5_cartan.det())
    require(det_d5 == 4, "D5 control determinant changed")
    require(det_a5 != det_d5, "natural response lattice accidentally equals D5")
    return {
        "det_A5": det_a5,
        "det_D5": det_d5,
        "det_3_sector": int(gram_v3.det()),
        "det_2_sector": int(gram_v2.det()),
        "combined_det": det_combined,
        "gluing_index_in_A5": index,
    }


def check_oriented_boundary_countercontrol() -> dict[str, object]:
    """Fixing a cyclic order reduces S4 to D4 and refines the 3+2 split.

    Class order: e, r^2, {r,r^3}, diagonal reflections (transpositions),
    edge reflections (double transpositions).
    """
    restricted_3 = (3, -1, -1, 1, -1)
    restricted_2 = (2, 2, 0, 0, 2)
    irrep_e = (2, -2, 0, 0, 0)
    irrep_b1 = (1, 1, -1, 1, -1)
    irrep_a1 = (1, 1, 1, 1, 1)
    irrep_b2 = (1, 1, -1, -1, 1)
    require(
        restricted_3 == tuple(a + b for a, b in zip(irrep_e, irrep_b1)),
        "3-sector restriction is not E+B1",
    )
    require(
        restricted_2 == tuple(a + b for a, b in zip(irrep_a1, irrep_b2)),
        "2-sector restriction is not A1+B2",
    )
    return {
        "S4_3_restricts_to_D4": "2+1",
        "S4_2_restricts_to_D4": "1+1",
        "oriented_shape_irreps": "2+1+1+1",
    }


def representation_on_basis(
    ambient_action: sp.Matrix, basis: sp.Matrix
) -> sp.Matrix:
    """Restrict an invariant ambient action to a full-column-rank basis."""
    left_inverse = (basis.T.conjugate() * basis).inv() * basis.T.conjugate()
    out = sp.simplify(left_inverse * ambient_action * basis)
    require(
        sp.simplify(ambient_action * basis - basis * out) == sp.zeros(*basis.shape),
        "basis is not invariant",
    )
    return out


def hom_space(left: list[sp.Matrix], right: list[sp.Matrix]) -> list[sp.Matrix]:
    """Exact basis of intertwiners T with left[g] T = T right[g]."""
    require(len(left) == len(right), "generator lists differ")
    n = left[0].rows
    variables = sp.symbols(f"t0:{n*n}")
    t = sp.Matrix(n, n, variables)
    equations = []
    for a, b in zip(left, right):
        equations.extend(list(a * t - t * b))
    coeff, rhs = sp.linear_eq_to_matrix(equations, variables)
    require(rhs == sp.zeros(rhs.rows, 1), "intertwiner equations became affine")
    return [sp.Matrix(n, n, list(v)) for v in coeff.nullspace()]


def maximum_rank_and_witness(basis: list[sp.Matrix]) -> tuple[int, sp.Matrix | None]:
    """Search small exact coefficients; enough for these multiplicity-free modules."""
    best_rank = -1
    best = None
    for coefficients in product((-2, -1, 1, 2), repeat=len(basis)):
        candidate = sum((c * b for c, b in zip(coefficients, basis)), sp.zeros(5))
        rank = candidate.rank()
        if rank > best_rank:
            best_rank = rank
            best = candidate
        if rank == 5:
            return rank, candidate
    return best_rank, best


def check_h6_carrier_representation() -> dict[str, object]:
    """Compare the response shape with H6 z34's actual V4/D4 carrier.

    H6 notation:
      F       = the three-dimensional family representation;
      Lambda2F = the actual colour representation;
      weak     = span(1-z^4, z^2).

    The network 3-sector is F, not Lambda2F.  They agree on V4 but differ on
    the odd D4 elements.  Tensoring the whole response shape with the real
    sign line det(F) repairs the full D4 representation exactly.
    """
    # D4 on four terminals/marks in cyclic order 0,1,2,3.
    d4_perms = {
        "e": (0, 1, 2, 3),
        "r": (1, 2, 3, 0),
        "r2": (2, 3, 0, 1),
        "r3": (3, 0, 1, 2),
        "s": (0, 3, 2, 1),
        "rs": (1, 0, 3, 2),
        "r2s": (2, 1, 0, 3),
        "r3s": (3, 2, 1, 0),
    }
    sign = {
        name: int(sp.det(sp.Matrix(4, 4, lambda i, j: int(i == p[j]))))
        for name, p in d4_perms.items()
    }
    require(sign == {"e": 1, "r": -1, "r2": 1, "r3": -1,
                     "s": -1, "rs": 1, "r2s": -1, "r3s": 1},
            f"wrong terminal sign character: {sign}")

    # Network shape basis = standard 3-sector plus incidence-kernel 2-sector.
    a3_roots = sp.Matrix(
        [[1, 0, 0], [-1, 1, 0], [0, -1, 1], [0, 0, -1]]
    )
    network_3_basis = incidence(4).T * a3_roots
    network_2_basis = sp.Matrix(
        [[1, 0], [0, 1], [-1, -1], [-1, -1], [0, 1], [1, 0]]
    )
    network_basis = network_3_basis.row_join(network_2_basis)
    network = {
        name: representation_on_basis(edge_permutation_matrix(p), network_basis)
        for name, p in d4_perms.items()
    }

    # Actual H6 carrier f=Q/(z^4-1), in the z34 V4 split.
    # Pullback by r:z->iz and s:z->1/z on monomials 1,z,...,z^4.
    carrier_r = sp.diag(1, sp.I, -1, -sp.I, 1)
    carrier_s = sp.zeros(5)
    for k in range(5):
        carrier_s[4 - k, k] = -1
    carrier_ambient = {
        "e": sp.eye(5),
        "r": carrier_r,
        "r2": carrier_r**2,
        "r3": carrier_r**3,
        "s": carrier_s,
        "rs": carrier_r * carrier_s,
        "r2s": carrier_r**2 * carrier_s,
        "r3s": carrier_r**3 * carrier_s,
    }
    e = sp.eye(5)
    actual_colour_basis = sp.Matrix.hstack(e[:, 1], e[:, 3], e[:, 0] + e[:, 4])
    actual_weak_basis = sp.Matrix.hstack(e[:, 0] - e[:, 4], e[:, 2])
    actual_basis = actual_colour_basis.row_join(actual_weak_basis)
    actual = {
        name: representation_on_basis(action, actual_basis)
        for name, action in carrier_ambient.items()
    }

    class_names = ("e", "r2", "r", "s", "rs")
    char_network = tuple(sp.trace(network[g]) for g in class_names)
    char_actual = tuple(sp.trace(actual[g]) for g in class_names)
    char_f = tuple(sp.trace(network[g][:3, :3]) for g in class_names)
    char_network_weak = tuple(sp.trace(network[g][3:, 3:]) for g in class_names)
    char_l2f = tuple(sp.trace(actual[g][:3, :3]) for g in class_names)
    char_actual_weak = tuple(sp.trace(actual[g][3:, 3:]) for g in class_names)
    char_sign = tuple(sign[g] for g in class_names)

    require(char_f == (3, -1, -1, 1, -1), f"network 3-sector is not F: {char_f}")
    require(char_l2f == (3, -1, 1, -1, -1), f"actual colour is not Lambda2F: {char_l2f}")
    require(char_network_weak == char_actual_weak == (2, 2, 0, 0, 2),
            "weak doublets differ")
    require(char_actual == tuple(a * d for a, d in zip(char_network, char_sign)),
            "actual carrier is not the sign twist of network shape")
    require(char_l2f == tuple(a * d for a, d in zip(char_f, char_sign)),
            "Lambda2F is not det(F) tensor F")

    # Untwisted full D4: no invertible intertwiner, maximum rank four.
    full_names = ("r", "s")
    hom_d4 = hom_space([actual[g] for g in full_names], [network[g] for g in full_names])
    rank_d4, _ = maximum_rank_and_witness(hom_d4)
    require(rank_d4 == 4, f"untwisted D4 maximum rank should be 4, got {rank_d4}")

    # On the even Klein four group the sign line is trivial: exact full-rank match.
    v4_names = ("r2", "rs")
    hom_v4 = hom_space([actual[g] for g in v4_names], [network[g] for g in v4_names])
    rank_v4, witness_v4 = maximum_rank_and_witness(hom_v4)
    require(rank_v4 == 5 and witness_v4 is not None, "V4 intertwiner is not invertible")

    # With the real sign/determinant line, the full D4 intertwiner is invertible.
    twisted_network = {g: sign[g] * network[g] for g in d4_perms}
    hom_twisted = hom_space(
        [actual[g] for g in full_names], [twisted_network[g] for g in full_names]
    )
    rank_twisted, witness_twisted = maximum_rank_and_witness(hom_twisted)
    require(rank_twisted == 5 and witness_twisted is not None,
            "sign-twisted D4 intertwiner is not invertible")
    for g in d4_perms:
        require(
            sp.simplify(actual[g] * witness_twisted - witness_twisted * twisted_network[g])
            == sp.zeros(5),
            f"twisted intertwiner fails on {g}",
        )

    # Sign/lift typing: a real 1D determinant line, not a mu4 phase or Spin lift.
    require(all(v in (-1, 1) for v in sign.values()), "sign line is not real")
    require(sign["r"] == sign["s"] == -1 and sign["r2"] == sign["rs"] == 1,
            "sign line is not det(F)")
    return {
        "class_order": class_names,
        "network_character": char_network,
        "actual_character": char_actual,
        "F_character": char_f,
        "Lambda2F_character": char_l2f,
        "weak_character": char_actual_weak,
        "untwisted_D4_hom_dimension": len(hom_d4),
        "untwisted_D4_max_rank": rank_d4,
        "V4_hom_dimension": len(hom_v4),
        "V4_max_rank": rank_v4,
        "sign_twisted_D4_hom_dimension": len(hom_twisted),
        "sign_twisted_D4_max_rank": rank_twisted,
        "sign_character": char_sign,
        "twisted_intertwiner_det": sp.simplify(witness_twisted.det()),
    }


def main() -> None:
    proj = check_projectors()
    families = check_minimality_and_counterfamilies()
    passive = check_passive_response()
    lattice = check_integral_lattice_countercontrol()
    oriented = check_oriented_boundary_countercontrol()
    h6 = check_h6_carrier_representation()

    print("STRUCTURE TEST: PASS")
    print(f"S4 projector ranks: {proj['ranks']} (scale + 3 + 2)")
    print(f"S4 equivariance checks: {proj['commute_checks']}")
    print(f"two-subset character: {proj['characters']} = 1 + 3 + 2")
    print("terminal counterfamilies:")
    for n in (2, 3, 4, 5, 6):
        print(f"  n={n}: {families[n]} -> shape {families[n][1:]}")
    print(
        "passive response: rank=3, tree cofactor=%s, min sampled energy=%s"
        % (passive["tree_cofactor"], passive["min_test_energy"])
    )
    print(
        "integral countercontrol: det(A5)=%s != det(D5)=%s; "
        "orthogonal sector sublattice det=%s, index=%s"
        % (
            lattice["det_A5"],
            lattice["det_D5"],
            lattice["combined_det"],
            lattice["gluing_index_in_A5"],
        )
    )
    print(
        "orientation countercontrol: S4 3 -> D4 %s; S4 2 -> D4 %s; shape -> %s"
        % (
            oriented["S4_3_restricts_to_D4"],
            oriented["S4_2_restricts_to_D4"],
            oriented["oriented_shape_irreps"],
        )
    )
    print(f"H6 class order: {h6['class_order']}")
    print(
        "H6 characters: network=%s, actual=%s; F=%s, Lambda2F=%s, weak=%s"
        % (
            h6["network_character"],
            h6["actual_character"],
            h6["F_character"],
            h6["Lambda2F_character"],
            h6["weak_character"],
        )
    )
    print(
        "H6 intertwiners: untwisted D4 max rank=%s; V4 max rank=%s; "
        "det(F)-twisted D4 max rank=%s, witness det=%s"
        % (
            h6["untwisted_D4_max_rank"],
            h6["V4_max_rank"],
            h6["sign_twisted_D4_max_rank"],
            h6["twisted_intertwiner_det"],
        )
    )
    print("VERDICT: structural stage-1 candidate; not P1/P2 or a full TFPT origin")


if __name__ == "__main__":
    main()
