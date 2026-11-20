"""Exact controls for reversible hidden-memory clocks and erased permutation histories."""
import argparse
from collections import Counter
from fractions import Fraction
from itertools import permutations, product
import json
from math import factorial, prod
from pathlib import Path


def compositions(total):
    for first, second, third in product(range(total + 1), repeat=3):
        fourth = total - first - second - third
        if fourth >= 0:
            yield first, second, third, fourth


def conjugation_sectors(multiplicities):
    return tuple(sum(multiplicities[index] * multiplicities[(index + residue) % 4]
                     for index in range(4)) for residue in range(4))


def clock_control():
    native = (23, 14, 17, 14)
    target_tr2 = sum((-1) ** index * count for index, count in enumerate(native))
    combined = []
    for multiplicities in compositions(8):
        adjoint = conjugation_sectors(multiplicities)
        candidate = (adjoint[0] - 1 + 5, *adjoint[1:])
        combined.append(candidate)
    half_turn_traces = sorted({sum((-1) ** k * count for k, count in enumerate(candidate))
                               for candidate in combined})
    assert len(combined) == 165
    assert half_turn_traces == [4, 8, 20, 40, 68]
    assert target_tr2 == 12 and target_tr2 not in half_turn_traces
    print('Clock traces J² on visible 5 plus traceless matrices of size 8:', half_turn_traces)
    print('R023 J² trace:', target_tr2, 'compatible phase allocations at size 8: 0 of', len(combined))

    regular_grade = (5 + 16 - 1, 16, 16, 16)
    assert regular_grade != native
    print('Visible 5 plus nonconstant functions on a regular 64-state C4 grade:', regular_grade)

    # Ordinary geometric quarter-turns in finite SO(3) modules obey an
    # additive mod-8 congruence. Spin l+4 adds 8 to dimension and keeps both
    # traces, so the check on residues l=0,1,2,3 proves the general identity.
    for l in range(4):
        weights = range(-l,l+1)
        quarter_trace = sum((1,0,-1,0)[m%4] for m in weights)
        half_trace = sum((1,-1)[m%2] for m in weights)
        assert (2*l+1 + 2*quarter_trace - 3*half_trace) % 8 == 0
    assert (68 + 2*6 - 3*12) % 8 == 4
    print('Native SO(3) quarter-turn congruence: residue 4, required 0 modulo 8')

    host_options = []
    for multiplicities in compositions(9):
        adjoint = conjugation_sectors(multiplicities)
        if all(found >= needed for found, needed in zip(adjoint, native)):
            host_options.append(multiplicities)
    assert host_options
    print('Size 9 host capacities, with extra invariant choices allowed:', len(host_options), 'of', len(list(compositions(9))))
    print('Size 9 example phase counts and full adjoint capacity:', host_options[0], conjugation_sectors(host_options[0]))


def native_response(evidence_directory):
    import sympy as sp

    source_path = evidence_directory / 'RECURSIVE-WORD-SOURCE.json'
    source = json.loads(source_path.read_text())
    operator = source['operators']
    transfer = sp.Matrix(operator['L']).applyfunc(sp.sympify)
    clock = sp.Matrix(operator['J']).applyfunc(sp.sympify)
    metric = sp.diag(*map(sp.sympify, source['metric_diagonal']))
    reference = sp.eye(68)[:, 4]
    assert (reference.conjugate().T * metric * reference)[0] == 1
    assert clock * reference == reference
    lifted = transfer * reference
    history = lifted
    moments = []
    for _ in range(4):
        moments.append(sp.simplify((lifted.conjugate().T * metric * history)[0]))
        history = clock * history
    assert history == lifted
    assert moments == [6144, 0, -2048, 0]
    first_return = sp.simplify((reference.conjugate().T * metric * lifted)[0])
    assert first_return == 32
    assert tuple(moment / first_return ** 2 for moment in moments) == (6, 0, -2, 0)
    phase_weights = [sp.simplify(sum(sp.I ** (-phase * step) * response
                                    for step, response in enumerate(moments)) / 4)
                     for phase in range(4)]
    assert phase_weights == [1024, 2048, 1024, 2048]
    plus_port = [sp.Rational(1, 2) * (1 + moment / moments[0]) for moment in moments]
    assert plus_port == [1, sp.Rational(1, 2), sp.Rational(1, 3), sp.Rational(1, 2)]
    print('Native f4 exact words L J^k L:', moments)
    print('Native f4 relative to its first return squared:', [moment / first_return ** 2 for moment in moments])
    print('Native f4 exact clock-sector weights:', phase_weights)
    print('Native f4 coherent X-basis plus-port probabilities:', plus_port)


def cross_block_response():
    import sympy as sp

    # Geometric quarter turn on the traceless diagonal and off-diagonal
    # parts of a selected real quadrupole frame. Derive both from rotations
    # on symmetric tensors, before comparing any native amplitudes.
    spatial_rotation = sp.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    diagonal = (sp.diag(1, -1, 0), sp.diag(1, 1, -2))
    off_diagonal = tuple(sp.eye(3)[:, i] * sp.eye(3)[j, :] + sp.eye(3)[:, j] * sp.eye(3)[i, :]
                         for i, j in ((0, 1), (0, 2), (1, 2)))

    def induced_action(basis):
        return sp.Matrix([[sp.trace(left * spatial_rotation * right * spatial_rotation.T)
                           / sp.trace(left ** 2) for right in basis] for left in basis])

    block_actions = (induced_action(diagonal), induced_action(off_diagonal))
    all_moments = {}
    unnormalized_cross = None
    for source, source_action in zip((2, 3), block_actions):
        for target, target_action in zip((2, 3), block_actions):
            action = sp.kronecker_product(target_action, source_action.inv().T)
            assert action ** 4 == sp.eye(source * target)
            traces = tuple(sp.trace(action ** order) for order in range(4))
            all_moments[source, target] = tuple(value / (source * target) for value in traces)
            if (source, target) == (2, 3):
                unnormalized_cross = traces
    assert unnormalized_cross == (6, 0, -2, 0)
    matching_blocks = [block for block, moments in all_moments.items() if moments == (1, 0, -sp.Rational(1, 3), 0)]
    assert matching_blocks == [(2, 3), (3, 2)]
    print('Quarter-turn normalized traces on the four ordered blocks:', all_moments)
    print('Native normalized first-return match:', matching_blocks, 'of 4 ordered blocks, or 1 of 3 types')


def distinct_cycles(mapping):
    unvisited = set(range(len(mapping)))
    cycles = 0
    while unvisited:
        next_vertex = min(unvisited)
        cycles += 1
        while next_vertex in unvisited:
            unvisited.remove(next_vertex)
            next_vertex = mapping[next_vertex]
    return cycles


def erase_vertex(mapping):
    removed = len(mapping) - 1
    return tuple(mapping[removed] if destination == removed else destination
                 for destination in mapping[:-1])


def weight(mapping, concentration):
    return concentration ** distinct_cycles(mapping) / prod(concentration + step for step in range(len(mapping)))


def erasure_control():
    cycle = tuple(range(1, 7)) + (0,)
    for deletions in range(1, len(cycle) + 1):
        targets = cycle[:deletions]
        completions = [chosen + cycle[deletions:] for chosen in permutations(targets)]
        cyclic = [rule for rule in completions if distinct_cycles(rule) == 1]
        assert len(completions) == factorial(deletions)
        assert len(cyclic) == factorial(deletions - 1)
        print('Erased outgoing transitions:', deletions, 'reversible completions:', len(completions),
              'one-cycle completions:', len(cyclic))
    for parameter in (Fraction(1, 2), Fraction(1), Fraction(2)):
        for size in range(2, 7):
            coarse = Counter()
            for mapping in permutations(range(size)):
                coarse[erase_vertex(mapping)] += weight(mapping, parameter)
            direct = {mapping: weight(mapping, parameter) for mapping in permutations(range(size - 1))}
            assert coarse == direct and sum(coarse.values()) == 1
        print('Deletion-consistent concentration:', parameter, 'two-vertex swap probability:', weight((1, 0), parameter))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native', type=Path, help='Directory containing RECURSIVE-WORD-SOURCE.json')
    options = parser.parse_args()
    clock_control()
    cross_block_response()
    if options.native is not None:
        native_response(options.native)
    erasure_control()
