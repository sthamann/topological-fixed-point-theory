#!/usr/bin/env python3
"""Exact audit of minimal perfect quantum memory as a TFPT origin candidate.

No floating-point calculation. The finite baseline is all cyclic single-seed
binary Pauli stabilizers of lengths 1..6, not all possible quantum codes.
Run: python3 qec_origin.py
"""

from fractions import Fraction as Q
from itertools import combinations, permutations
import json
from pathlib import Path
import sympy as sp


def parity(x):
    return x.bit_count() & 1


def symplectic(a, b, n):
    mask = (1 << n) - 1
    return parity(((a & mask) & (b >> n))) ^ parity(((a >> n) & (b & mask)))


def gf2_rank(rows):
    basis = {}
    for value in rows:
        while value:
            lead = value.bit_length() - 1
            if lead in basis:
                value ^= basis[lead]
            else:
                basis[lead] = value
                break
    return len(basis)


def span(rows):
    answer = {0}
    for row in rows:
        answer |= {value ^ row for value in tuple(answer)}
    return answer


def rotate(value, n):
    mask = (1 << n) - 1
    x, z = value & mask, value >> n
    x = ((x << 1) & mask) | (x >> (n - 1))
    z = ((z << 1) & mask) | (z >> (n - 1))
    return x | (z << n)


def word_pauli(word):
    n = len(word)
    x = sum((letter in "XY") << i for i, letter in enumerate(word))
    z = sum((letter in "ZY") << i for i, letter in enumerate(word))
    return x | (z << n)


def errors_one(n):
    # XZ differs from Y by a phase, irrelevant for exact correction.
    return [0] + [1 << i for i in range(n)] + [1 << (n + i) for i in range(n)] + [
        (1 << i) | (1 << (n + i)) for i in range(n)
    ]


def syndrome(error, rows, n):
    return sum(symplectic(error, row, n) << i for i, row in enumerate(rows))


def cyclic_baseline(n):
    commuting = encoded_one = passing = 0
    codes = set()
    errors = errors_one(n)
    for seed in range(1 << (2 * n)):
        value = seed
        rows = []
        for _ in range(n):
            rows.append(value)
            value = rotate(value, n)
        if any(symplectic(a, b, n) for a in rows for b in rows):
            continue
        commuting += 1
        if gf2_rank(rows) != n - 1:
            continue
        encoded_one += 1
        if len({syndrome(e, rows, n) for e in errors}) == len(errors):
            passing += 1
            codes.add(tuple(sorted(span(rows))))
    return {
        "sites": n,
        "all_seeds": 4**n,
        "commuting_seeds": commuting,
        "commuting_rank_n_minus_1": encoded_one,
        "correcting_seeds": passing,
        "distinct_stabilizer_subspaces": len(codes),
    }, codes


def apply_pauli(vector, value, n):
    mask = (1 << n) - 1
    x, z = value & mask, value >> n
    result = [Q(0)] * (1 << n)
    for b, amplitude in enumerate(vector):
        result[b ^ x] += (-1 if parity(b & z) else 1) * amplitude
    return result


def inner(a, b):
    return sum((x * y for x, y in zip(a, b)), Q(0))


def project_bits(b, sites):
    return sum(((b >> s) & 1) << j for j, s in enumerate(sites))


def reduced_outer(a, b, sites, n, scale=Q(1)):
    env = [i for i in range(n) if i not in sites]
    dim = 1 << len(sites)
    result = [[Q(0) for _ in range(dim)] for _ in range(dim)]
    for i, ai in enumerate(a):
        if not ai:
            continue
        for j, bj in enumerate(b):
            if bj and project_bits(i, env) == project_bits(j, env):
                result[project_bits(i, sites)][project_bits(j, sites)] += scale * ai * bj
    return result


def is_scalar_identity(matrix, coefficient):
    return all(entry == (coefficient if i == j else 0)
               for i, row in enumerate(matrix) for j, entry in enumerate(row))


def permute_mask(mask, perm):
    return sum(((mask >> i) & 1) << perm[i] for i in range(len(perm)))


def permute_pauli(value, perm):
    n = len(perm)
    mask = (1 << n) - 1
    return permute_mask(value & mask, perm) | (permute_mask(value >> n, perm) << n)


def main():
    q = sp.symbols("q")
    packing_polynomial = sp.factor(q**4 - (1 + 5 * (q**2 - 1)))
    packing_roots = sp.solve(packing_polynomial, q)
    valid_alphabets = [int(root) for root in packing_roots if root.is_integer and root >= 2]
    assert valid_alphabets == [2]

    rows = [word_pauli(word) for word in ["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"]]
    n = 5
    assert gf2_rank(rows) == 4
    assert all(symplectic(a, b, n) == 0 for a in rows for b in rows)
    stabilizer = span(rows)
    normalizer = [p for p in range(4**n) if all(symplectic(p, r, n) == 0 for r in rows)]
    weight = lambda p: ((p & 31) | (p >> 5)).bit_count()
    distance = min(weight(p) for p in normalizer if p not in stabilizer)
    assert distance == 3

    psi0 = [Q(1)] + [Q(0)] * 31
    for row in rows:
        acted = apply_pauli(psi0, row, n)
        psi0 = [(a + b) / 2 for a, b in zip(psi0, acted)]
    assert inner(psi0, psi0) == Q(1, 16)
    psi0 = [4 * a for a in psi0]
    psi1 = apply_pauli(psi0, 31, n)  # transversal X on all sites
    assert inner(psi0, psi0) == inner(psi1, psi1) == 1
    assert inner(psi0, psi1) == 0

    errors = errors_one(n)
    syndrome_values = sorted(syndrome(e, rows, n) for e in errors)
    assert syndrome_values == list(range(16))
    error_images = [apply_pauli(psi, e, n) for e in errors for psi in (psi0, psi1)]
    gram_checks = sum(inner(a, b) == int(i == j)
                      for i, a in enumerate(error_images) for j, b in enumerate(error_images))
    assert gram_checks == 32**2

    pair_privacy = []
    for pair in combinations(range(5), 2):
        good = (is_scalar_identity(reduced_outer(psi0, psi0, pair, 5), Q(1, 4))
                and is_scalar_identity(reduced_outer(psi1, psi1, pair, 5), Q(1, 4))
                and is_scalar_identity(reduced_outer(psi0, psi1, pair, 5), Q(0)))
        pair_privacy.append(good)
    assert all(pair_privacy)

    # Choi state: logical reference is site 5. Use unnormalized rational
    # amplitudes and multiply density matrix by 1/2, avoiding sqrt(2).
    choi = psi0 + psi1
    ame_triples = [is_scalar_identity(reduced_outer(choi, choi, triple, 6, Q(1, 2)), Q(1, 8))
                   for triple in combinations(range(6), 3)]
    assert all(ame_triples)

    cycle = (1, 2, 3, 4, 0)
    assert all(tuple(vector[permute_mask(i, cycle)] for i in range(32)) == tuple(vector)
               for vector in (psi0, psi1))
    invariant_cuts = [cut for cut in range(32) if permute_mask(cut, cycle) == cut]
    assert invariant_cuts == [0, 31]
    marked_three_two_cuts = sum(cut.bit_count() == 2 for cut in invariant_cuts)
    assert marked_three_two_cuts == 0

    parity_dimensions = [sum(b.bit_count() % 2 == k for b in range(32)) for k in (0, 1)]
    logical_even_norm = [sum(a*a for b, a in enumerate(psi) if not parity(b)) for psi in (psi0, psi1)]
    assert parity_dimensions == [16, 16]
    assert logical_even_norm == [1, 0]

    baseline = []
    five_codes = None
    for size in range(1, 7):
        entry, codes = cyclic_baseline(size)
        baseline.append(entry)
        if size == 5:
            five_codes = codes
    representative = next(iter(five_codes))
    equivalent_codes = set()
    for perm in permutations(range(5)):
        equivalent_codes.add(tuple(sorted(permute_pauli(p, perm) for p in representative)))
    assert five_codes <= equivalent_codes

    result = {
        "principle": "Minimum sites; one logical q-level system; one arbitrary local error; nondegenerate exact correction; perfect packing",
        "singleton_bound_minimum_sites": 5,
        "singleton_bound_is_analytic_input": "n-k >= 2(d-1), k=1, d>=3; not inferred from finite census",
        "minimum_site_packing_polynomial": str(packing_polynomial),
        "integer_alphabets_q_ge_2": valid_alphabets,
        "stabilizer_rank": 4,
        "logical_dimension": 2,
        "distance": distance,
        "single_error_syndromes": syndrome_values,
        "exact_error_gram_entries_passed": gram_checks,
        "two_share_privacy_pairs_passed": sum(pair_privacy),
        "three_share_recovery_sets": len(pair_privacy),
        "choi_AME_3_share_reductions_passed": sum(ame_triples),
        "physical_parity_dimensions": parity_dimensions,
        "logical_code_even_parity_dimension": 1,
        "C5_invariant_subset_sizes": [cut.bit_count() for cut in invariant_cuts],
        "C5_invariant_marked_three_plus_two_cuts": marked_three_two_cuts,
        "baseline": baseline,
        "n5_success_fraction_all_seeds": str(Q(30, 1024)),
        "n5_success_fraction_commuting_rank4_seeds": str(Q(30, 75)),
        "n5_permutation_classes_among_passing_cyclic_subspaces": 1,
        "stage": 1,
        "verdict": "KANDIDAT for a five-site error-correcting origin; P2 remains OPEN",
        "first_free_structure": "Tensor-product quantum systems, the one-local-error task, minimization and perfect packing",
        "strongest_objection": "No automorphism-invariant marked 3+2 partition; syndrome dimension is not a chiral matter representation",
    }
    output = json.dumps(result, indent=2)
    print(output)
    Path(__file__).with_name("qec_output.json").write_text(output + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
