"""Exact support and rank check for the marked duad word law."""

from itertools import combinations, permutations

import sympy as sp


def duad_image(duad, permutation):
    return tuple(sorted(permutation[point] for point in duad))


def permutation_matrix(duads, permutation):
    index = {duad: position for position, duad in enumerate(duads)}
    matrix = sp.zeros(len(duads))
    for position, duad in enumerate(duads):
        matrix[position, index[duad_image(duad, permutation)]] = 1
    return matrix


def verify_word_law():
    points = tuple(range(6))
    mark = points[-1]
    identity = points
    duads = list(combinations(points, 2))
    outer = [duad for duad in duads if mark in duad]
    inner = [duad for duad in duads if mark not in duad]
    index = {duad: position for position, duad in enumerate(duads)}

    doily_adjacency = sp.Matrix(
        15, 15, lambda row, col: int(not set(duads[row]).intersection(duads[col]))
    )
    petersen_adjacency = sp.Matrix(
        10, 10, lambda row, col: int(not set(inner[row]).intersection(inner[col]))
    )
    restrict = sp.Matrix(10, 15, lambda row, col: int(inner[row] == duads[col]))
    petersen_five = (petersen_adjacency + 2 * sp.eye(10)) / 3 - sp.ones(10) / 6
    corrected = (
        (doily_adjacency + 3 * sp.eye(15)) / 9
        - restrict.T * petersen_five * restrict / 3
    )

    # A positive mixture cannot use a permutation at a forbidden entry.
    candidates = []
    for permutation in permutations(points):
        if all(
            corrected[index[duad], index[duad_image(duad, permutation)]] > 0
            for duad in duads
        ):
            candidates.append(permutation)
    assert len(candidates) == 46
    assert identity in candidates
    assert all(
        sum(duad_image(duad, permutation) == duad for duad in outer) == 1
        for permutation in candidates
        if permutation != identity
    )

    # The outer diagonal sum fixes the identity weight. It already
    # accounts for the whole inner diagonal, excluding any more fixed duads.
    outer_diagonal_sum = sum(corrected[index[duad], index[duad]] for duad in outer)
    identity_weight = (outer_diagonal_sum - 1) / 4
    assert identity_weight == sp.Rational(1, 6)
    assert all(corrected[index[duad], index[duad]] == identity_weight for duad in inner)
    support = [
        permutation
        for permutation in candidates
        if permutation == identity
        or all(duad_image(duad, permutation) != duad for duad in inner)
    ]

    expected = {identity}
    for partner in points[:-1]:
        remaining = [point for point in points if point not in (mark, partner)]
        first, *tail = remaining
        for order in permutations(tail):
            cycle = (first, *order)
            word = list(identity)
            word[mark], word[partner] = partner, mark
            for point, successor in zip(cycle, cycle[1:] + cycle[:1]):
                word[point] = successor
            expected.add(tuple(word))
    assert set(support) == expected and len(support) == 31

    matrices = [permutation_matrix(duads, permutation) for permutation in support]
    columns = sp.Matrix.hstack(*(matrix.reshape(225, 1) for matrix in matrices))
    assert columns.rank() == 31
    weights = sp.Matrix(
        [identity_weight if permutation == identity else sp.Rational(1, 36)
         for permutation in support]
    )
    assert sum(weights) == 1 and columns * weights == corrected.reshape(225, 1)

    print("Exact 225-by-31 rank 31; unique law on the 31 supported permutations.")


if __name__ == "__main__":
    verify_word_law()
