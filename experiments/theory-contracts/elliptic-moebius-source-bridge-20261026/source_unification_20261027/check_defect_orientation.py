"""Exact three-mode CAR orientation test; not a physical source model.

All matrices act on the complete exterior algebra of the given three-point
defect. Phases are exponents of zeta_8, so no floating-point tolerance enters.
The infinite Hardy argument is separate from this finite action check.
"""
import json
from pathlib import Path

N = 3
DIM = 1 << N


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def create(j, mask):
    if mask & (1 << j):
        return None
    sign = (mask & ((1 << j) - 1)).bit_count() % 2
    return (mask | (1 << j), 4 * sign)


def annihilate(j, mask):
    if not mask & (1 << j):
        return None
    sign = (mask & ((1 << j) - 1)).bit_count() % 2
    return (mask ^ (1 << j), 4 * sign)


def compose_action(first, second, mask):
    intermediate = second(mask)
    if intermediate is None:
        return None
    result = first(intermediate[0])
    if result is None:
        return None
    return result[0], (result[1] + intermediate[1]) % 8


def exterior(permutation, phase):
    output = []
    for mask in range(DIM):
        occupied = [j for j in range(N) if mask & (1 << j)]
        image = [permutation[j] for j in occupied]
        inversions = sum(image[i] > image[j] for i in range(len(image)) for j in range(i + 1, len(image)))
        output.append((sum(1 << j for j in image), (phase * len(image) + 4 * inversions) % 8))
    return output


def multiply(left, right):
    return [(left[t][0], (e + left[t][1]) % 8) for t, e in right]


def inverse(matrix):
    output = [None] * DIM
    for m, (t, e) in enumerate(matrix):
        output[t] = (m, -e % 8)
    return output


identity = [(m, 0) for m in range(DIM)]
c = exterior([1, 2, 0], 0)
q = exterior([0, 2, 1], 1)
beta = exterior([0, 1, 2], 2)
parity = exterior([0, 1, 2], 4)
need(multiply(multiply(c, c), c) == identity, 'c^3')
need(multiply(q, q) == beta, 'Q^2=beta on full exterior algebra')
need(multiply(beta, beta) == parity, 'beta^2=fermion parity')
need(multiply(multiply(q, c), inverse(q)) == inverse(c), 'Q c Q^-1=c^-1')

car_tests = 0
for i in range(N):
    for j in range(N):
        for mask in range(DIM):
            t1 = compose_action(lambda m: annihilate(i, m), lambda m: create(j, m), mask)
            t2 = compose_action(lambda m: create(j, m), lambda m: annihilate(i, m), mask)
            if i == j:
                nonzero = [x for x in (t1, t2) if x is not None]
                need(nonzero == [(mask, 0)], 'CAR diagonal')
            else:
                need((t1 is None and t2 is None) or
                     (t1 is not None and t2 is not None and t1[0] == t2[0] and (t1[1] - t2[1]) % 8 == 4), 'CAR off diagonal')
            car_tests += 1

for matrix, perm, phase in [(q, [0, 2, 1], 1), (c, [1, 2, 0], 0)]:
    inv = inverse(matrix)
    for j in range(N):
        for mask in range(DIM):
            t, e = inv[mask]
            result = create(j, t)
            expected = create(perm[j], mask)
            if result is None:
                need(expected is None, 'zero covariance')
            else:
                out, finalphase = matrix[result[0]]
                actual = out, (e + result[1] + finalphase) % 8
                need(actual == (expected[0], (expected[1] + phase) % 8), 'creation covariance')

full = DIM - 1
need(beta[full] == (full, 6), 'filled determinant line beta=-i')
need(q[full] == (full, 7), 'filled determinant Q=zeta8^-1')
need(c[full] == (full, 0), 'filled determinant c=1')
need(parity[full] == (full, 4), 'filled determinant odd')
need((-beta[full][1]) % 8 == 2, 'dual determinant beta=+i')
need((-q[full][1]) % 8 == 1, 'dual determinant Q=zeta8')

# A state sees no ray phase; charged transition amplitudes can retain it.
need(all((matrix[full][1] - matrix[full][1]) % 8 == 0 for matrix in (c, q, beta)), 'density phase cancellation')
state = (0, 0)
for j in range(N):
    nxt = create(j, state[0])
    state = nxt[0], (state[1] + nxt[1]) % 8
need(state[0] == full, 'three-particle excitation')
for j in range(N):
    need(annihilate(j, 0) is None and create(j, full) is None, 'empty/filled covariance')

result = {
    'verdict': 'EXACT_THREE_MODE_KET_AND_FREDHOLM_DUAL_DISTINCT_WITH_SAME_ODD_PARITY',
    'source_scope': 'Given three-dimensional defect representation; no P1 physical source identification.',
    'dimension_defect': N, 'dimension_full_exterior_algebra': DIM,
    'car_basis_checks': car_tests, 'phases': 'exponents modulo 8 of zeta_8',
    'ket': {'c': c[full][1], 'Q': q[full][1], 'beta': beta[full][1], 'parity': 'odd'},
    'dual': {'c': -c[full][1] % 8, 'Q': -q[full][1] % 8, 'beta': -beta[full][1] % 8, 'parity': 'odd'},
    'beta_square_equals_full_fermion_parity': True,
    'all_creation_operators_covariant': True,
    'density_matrix_forgets_ray_phase': True,
}
print(json.dumps(result, sort_keys=True, indent=2))
