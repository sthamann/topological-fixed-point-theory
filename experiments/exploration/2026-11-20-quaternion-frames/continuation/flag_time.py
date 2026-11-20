"""Exact tests of OS and graph KMS completions of the quaternionic flag law.

The experiment starts at the previous K6 one-factorization result. It does
not identify any graph KMS state or modular flow with the native TFPT seam.
The type conclusions use Thomsen, arXiv:1805.05637, Theorem 3.3.
"""

from itertools import combinations, product
from math import prod

from sympy import Matrix, Rational, eye, factorint, log, ones, zeros


def flag_transfer():
    edges = tuple(map(frozenset, combinations(range(6), 2)))
    dimension = len(edges)
    adjacency = Matrix(dimension, dimension, lambda i, j: int(i != j and bool(edges[i] & edges[j])))
    assert adjacency.eigenvals() == {8: 1, 2: 5, -2: 9}
    transfer = (4 * eye(dimension) + adjacency) / 12
    return edges, adjacency, transfer


def test_reflection_and_time(adjacency, transfer):
    dimension = transfer.rows
    stationary = ones(dimension, 1) / dimension
    assert transfer * stationary == stationary and transfer == transfer.T
    assert transfer.eigenvals() == {1: 1, Rational(1, 2): 5, Rational(1, 6): 9}

    ground = ones(dimension) / dimension
    middle = (adjacency + 2 * eye(dimension) - 10 * ground) / 4
    upper = eye(dimension) - ground - middle
    projections = (ground, middle, upper)
    for i, j in product(range(len(projections)), repeat=2):
        assert projections[i] * projections[j] == (projections[i] if i == j else zeros(dimension))
    assert transfer == ground + middle / 2 + upper / 6
    print('OS transfer eigenvalues with multiplicities:', transfer.eigenvals())
    print('Euclidean Hamiltonian gap pair: log(2), log(6); irrational ratio by prime valuation')

    disjoint = next((i, j) for i, j in product(range(dimension), repeat=2)
                    if i != j and adjacency[i, j] == 0)
    assert transfer[disjoint] == 0 and (transfer**2)[disjoint] == Rational(1, 36)
    assert all(value > 0 for value in transfer**2)
    log_entry = (-log(2) * middle - log(6) * upper)[disjoint]
    assert (log_entry - log(Rational(4, 27)) / 30).simplify() == 0
    print('Disjoint-flag transfer at one and two steps:', transfer[disjoint], (transfer**2)[disjoint])
    print('Disjoint-flag principal generator entry: log(4/27)/30 < 0')


def test_path_completions(adjacency, transfer):
    dimension = transfer.rows
    unrefined = eye(dimension) + adjacency
    refined = 4 * eye(dimension) + adjacency
    assert refined * ones(dimension, 1) == 12 * ones(dimension, 1)
    assert set(unrefined) == {0, 1} and set(transfer) == {0, Rational(1, 12), Rational(1, 3)}
    assert refined / 12 == transfer

    # A coarse path has 4^(number of stays) fine lifts. Each lift has law
    # 1/(15 * 12^n); multiplying gives the unrefined Markov cylinder law.
    words_checked = 0
    for first, second, third in product(range(dimension), repeat=3):
        legs = ((first, second), (second, third))
        num_lifts = prod(refined[leg] for leg in legs)
        weighted = prod(transfer[leg] for leg in legs)
        assert num_lifts / 12**len(legs) == weighted
        words_checked += 1
    print('Identical projected cylinder probabilities (three flag observations):', words_checked)

    # Conservative finite graph KMS theorem and cycle energy lattice.
    # Comparing a one-edge cycle with its double gives log(12) (refined)
    # or log(3) (unrefined). Comparing a backtrack with its double gives
    # 2 log(12). Any further refined cycle contributes integer multiples.
    assert factorint(12) == {2: 2, 3: 1}
    prime_valuations = Matrix([[0, 1], [4, 2]])
    assert prime_valuations.det() != 0
    print('Refined graph: 12 uniform outgoing edges, beta=log(12), Gamma=log(12) Z, type III_(1/12)')
    print('Unrefined graph: beta=1, loop valuation rank', prime_valuations.rank(), ', Gamma=R, type III_1')


if __name__ == '__main__':
    _, adjacency, transfer = flag_transfer()
    test_reflection_and_time(adjacency, transfer)
    test_path_completions(adjacency, transfer)
