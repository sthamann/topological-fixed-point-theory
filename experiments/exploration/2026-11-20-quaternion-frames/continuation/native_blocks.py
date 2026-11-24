#!/usr/bin/env python3
"""Exact small-block derivation of the native R023 return rank.

The icosian frame supplies b,c. The native vertex rule, quartics and Fock metric
are the declared comparison interface.
The proof uses spectator conservation and active c-parity, without the saved
68-vector basis or the complete 726-dimensional neutral Fock enumeration.
"""

import argparse
import json
from collections import Counter, defaultdict
from functools import lru_cache
from math import factorial, prod
from pathlib import Path

from sympy import sympify
from sympy.polys.domains import QQ_I as G

from icosian_frame import derive_frame


UNIT = G.one
ZERO = G.zero
I = G(0, 1)


def add(*terms):
    out = defaultdict(lambda: ZERO)
    for factor, poly in terms:
        for monomial, coefficient in poly.items():
            out[monomial] += factor * coefficient
    return {monomial: coefficient for monomial, coefficient in out.items() if coefficient}


def multiply(a, b):
    out = defaultdict(lambda: ZERO)
    for first, x in a.items():
        for second, y in b.items():
            out[tuple(sorted(first + second))] += x * y
    return {monomial: coefficient for monomial, coefficient in out.items() if coefficient}


def conjugate(x):
    return G(x.x, -x.y)


def mode(index):
    return 1 + index // 4


@lru_cache(None)
def weight(monomial):
    return prod(factorial(n) for n in Counter(monomial).values()) * prod(mode(k) for k in monomial)


def active_basis(energy, c_parity, first=0):
    if energy == 0:
        if c_parity == 0:
            yield ()
        return
    for index in range(first, 4 * energy):
        for rest in active_basis(energy - mode(index), c_parity ^ (index % 2), index):
            yield (index,) + rest


def source_operators(b, c):
    def linear(direction, level):
        return {(4 * (level - 1) + j,): conjugate(value)
                for j, value in enumerate(direction) if value}

    exponent = [{(): UNIT}]
    for n in range(1, 4):
        exponent.append(add(*[(UNIT / n, multiply(linear(b, j), exponent[n - j]))
                               for j in range(1, n + 1)]))

    @lru_cache(None)
    def translate(monomial):
        result = {0: {(): UNIT}}
        for index in monomial:
            terms = defaultdict(list)
            for power, poly in result.items():
                terms[power].append((UNIT, multiply(poly, {(index,): UNIT})))
                if b[index % 4]:
                    terms[power - mode(index)].append((-b[index % 4], poly))
            result = {power: add(*polys) for power, polys in terms.items()}
        return result

    def coefficient(monomial, power):
        return add(*[(UNIT, multiply(poly, exponent[power - shift]))
                     for shift, poly in translate(monomial).items()
                     if 0 <= power - shift < len(exponent)])

    def source(monomial):
        terms = [(UNIT, multiply(linear(c, n), coefficient(monomial, -n - 1)))
                 for n in range(1, 4)]
        for i, index in enumerate(monomial):
            if c[index % 4]:
                shorter = monomial[:i] + monomial[i + 1:]
                terms.append((mode(index) * c[index % 4], coefficient(shorter, mode(index) - 1)))
        return add(*terms)

    return source


def return_operator(source, energy, c_parity):
    adjoint = defaultdict(list)
    images = {}
    for before in active_basis(energy, c_parity):
        images[before] = source(before)
        for after, coefficient in images[before].items():
            assert sum(map(mode, after)) == energy - 1
            assert sum(k % 2 for k in after) % 2 == 1 - c_parity
            adjoint[after].append((conjugate(coefficient) * G(weight(after)) / weight(before),
                                   {before: UNIT}))
    return_images = {before: add(*[(coefficient * scalar, basis)
                                   for after, coefficient in images[before].items()
                                   for scalar, basis in adjoint[after]])
                     for before in images}

    def apply(poly):
        return add(*[(coefficient, return_images[monomial])
                     for monomial, coefficient in poly.items() if monomial in return_images])

    return apply


def row_basis(polys):
    rows = {}
    for vector in polys:
        insert(rows, vector)
    return rows


def insert(rows, candidate):
    while candidate:
        pivot = min(candidate)
        leading = candidate[pivot]
        if pivot not in rows:
            rows[pivot] = add((UNIT / leading, candidate))
            return rows[pivot]
        candidate = add((UNIT, candidate), (-leading, rows[pivot]))
    return None


def audit(energy, parity, seeds, b, c):
    def eigenspace(monomial):
        if monomial[0] < 0:
            return 0  # Frozen all-positive spectator fourth-grade directions.
        return (4 - energy + sum(1 if k % 4 < 2 else -1 for k in monomial)) % 4

    def clock(poly):
        return {monomial: coefficient * I ** eigenspace(monomial)
                for monomial, coefficient in poly.items()}

    source = source_operators(b, c)
    repeated = return_operator(source, energy, parity)
    rows = row_basis(seeds)
    queue = list(rows.values())
    for vector in queue:
        for operator in (repeated, clock):
            new = insert(rows, operator(vector))
            if new is not None:
                queue.append(new)
    character = tuple(len(row_basis({m: v for m, v in vector.items() if eigenspace(m) == q}
                                    for vector in queue)) for q in range(4))
    active_projection_rank = len(row_basis({m: v for m, v in vector.items() if m[0] >= 0}
                                           for vector in queue))
    assert sum(character) == len(queue)
    return character, active_projection_rank, repeated, clock


def native_operator_control(native_json, b, c):
    data = json.loads(native_json.read_text())
    vector_polynomials = [{tuple(word): G.from_sympy(sympify(value))
                           for word, value in entry['polynomial']}
                          for entry in data['original_vector_basis']]
    source = source_operators(b, c)
    active_indices = (2, 3, 6, 7)

    @lru_cache(None)
    def source_return(energy, parity):
        return return_operator(source, energy, parity)

    @lru_cache(None)
    def return_monomial(original):
        spectators = tuple(k for k in original if k % 8 not in active_indices)
        active = tuple(4*(k//8) + active_indices.index(k%8) for k in original
                       if k % 8 in active_indices)
        if not active:
            return {}
        energy = sum(map(mode, active))
        parity = sum(k % 2 for k in active) % 2
        images = source_return(energy, parity)({active: UNIT})
        return {tuple(sorted(spectators + tuple(8*(k//4) + active_indices[k%4]
                                                 for k in monomial))): value
                for monomial, value in images.items()}

    def ret(poly):
        return add(*[(value, return_monomial(monomial)) for monomial, value in poly.items()])

    def clock(poly, direction):
        return {monomial: value * I ** (direction*sum(1 if k%8 < 4 else -1 for k in monomial) % 4)
                for monomial, value in poly.items()}

    count = 0
    actions = {'L': ret, 'J': lambda poly: clock(poly, 1), 'Jdagger': lambda poly: clock(poly,-1)}
    for label, action in actions.items():
        matrix = data['operators'][label]
        for i, poly in enumerate(vector_polynomials):
            expected = add(*[(G.from_sympy(sympify(row[i])), component)
                             for row, component in zip(matrix, vector_polynomials) if row[i] != '0'])
            assert action(poly) == expected, (label, i)
            count += 1
    assert count == 204
    print('Complete comparison: 204 native L/J/Jdagger images agree as exact Fock polynomials')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native', type=Path, help='Path to RECURSIVE-WORD-SOURCE.json')
    args = parser.parse_args()
    def active(data):
        assert all(pair == (0,0) for pair in data[:2])
        positive = [G(*data[i]) for i in (2,3)]
        return tuple(positive + [conjugate(x) for x in positive])

    b, c = map(active, derive_frame())
    # The quartet variables x0,x1 are spectators. Active variables 0,1 here
    # denote native x2,x3. Negative indices encode the frozen spectator terms.
    quartics = [{(0,) * 4: UNIT, (1,) * 4: UNIT, (-1,): UNIT},
                {(0, 0, 1, 1): UNIT, (-2,): UNIT}]
    plus = [{(0, 0): UNIT, (1, 1): UNIT}]
    minus = [{(0, 0): UNIT, (1, 1): -UNIT}]
    mixed = [{(0, 1): UNIT}]
    parts = ((4, 0, quartics), (2, 0, plus), (2, 0, minus), (2, 1, mixed))
    results = [audit(*part, b=b, c=c) for part in parts]
    characters = [result[0] for result in results]
    expected = [(15, 11, 13, 11), (3, 1, 1, 1), (3, 1, 1, 1), (2, 1, 2, 1)]
    assert characters == expected and results[0][1] == 50
    assert audit(4, 0, [{m: v for m, v in seed.items() if m[0] >= 0}
                        for seed in quartics], b=b, c=c)[0] == expected[0]
    print('Native ranks:', [sum(chars) for chars in characters], 'total:', sum(map(sum, characters)))
    print('Native clock dimensions by block:', characters)
    print('Native total clock dimensions:', tuple(map(sum, zip(*characters))))
    print('The frozen spectator terms add zero independent directions to the 50-block.')

    repeated, clock = results[-1][2:]
    f4 = (0, 1)
    after = repeated(mixed[0])
    first = after[f4]
    second = repeated(after)[f4]
    with_clock = repeated(clock(after)).get(f4, ZERO)
    assert (first, second, with_clock) == (G(32), G(6144), ZERO)
    print('Unit f4 amplitudes L, L², LJL:', tuple(map(str, (first, second, with_clock))))

    no_root = [audit(*part, b=(ZERO,) * 4, c=c)[0] for part in parts]
    print('Control with root translation b removed:', [sum(chars) for chars in no_root])
    assert sum(map(sum, no_root)) < 68
    if args.native:
        native_operator_control(args.native, b, c)


if __name__ == '__main__':
    main()
