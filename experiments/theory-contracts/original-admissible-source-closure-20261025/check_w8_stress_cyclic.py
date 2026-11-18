#!/usr/bin/env python3
"""Exact ORIGINAL c1 stress reconstruction in the existing W8 quarter source.

Reconstruct from all old quantum-completion coefficients, independently of
the old sparse PBW expansion. No source, state, current or coupling is added.
The folded c8 obstruction is retained ONLY in its smaller-generator scope.
The all-level original c1 cyclicity is proved in the accompanying report.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
import json
from math import factorial
from pathlib import Path
import sympy as sp


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), F(0))


def roots_e8():
    roots = []
    for i, j in combinations(range(8), 2):
        for a, b in product((-1, 1), repeat=2):
            root = [F(0)] * 8
            root[i], root[j] = F(a), F(b)
            roots.append(tuple(root))
    for signs in product((-1, 1), repeat=8):
        if sum(x < 0 for x in signs) % 2 == 0:
            roots.append(tuple(F(x, 2) for x in signs))
    need(len(set(roots)) == 240, 'E8 root count')
    need(all(dot(a, a) == 2 for a in roots), 'root norm')
    return roots


def tangential(vector, d):
    multiple = dot(d, vector) / dot(d, d)
    return tuple(x - multiple * y for x, y in zip(vector, d))


def contraction_constant_and_first(pairs, d):
    """Q_ab -> (-1)^(a+b-2)d²; only a=1 or b=1 creates level1."""
    d2 = dot(d, d)
    constants = [F((-1) ** (a + b - 2)) * d2 for a, b in pairs]
    first_multiples = [F(((-1) ** (b - 1) if a == 1 else 0)
                         + ((-1) ** (a - 1) if b == 1 else 0))
                       for a, b in pairs]
    constant = F(1)
    for value in constants:
        constant *= value
    first_multiple = F(0)
    for i, value in enumerate(first_multiples):
        term = value
        for j, other in enumerate(constants):
            if i != j:
                term *= other
        first_multiple += term
    first = tuple(first_multiple * value for value in d)
    need(tangential(first, d) == (F(0),) * 8, 'nonradial contraction')
    return constant, first_multiple


@lru_cache(None)
def partitions(n, largest=8):
    if n == 0:
        return ((),)
    return tuple((k,) + p for k in range(min(n, largest), 0, -1)
                 for p in partitions(n - k, k))


def fermion_move(mask, source, target):
    if not (mask >> source) & 1 or (mask >> target) & 1:
        return None
    sign = (-1) ** (mask & ((1 << source) - 1)).bit_count()
    intermediate = mask ^ (1 << source)
    sign *= (-1) ** (intermediate & ((1 << target) - 1)).bit_count()
    return intermediate | (1 << target), sign


def original_stress_action(vector, n, alpha):
    """Literal old CAR mode Lcal_-n=sum_r(r+alpha+n/2):a*_(r+n)a_r:."""
    out = {}
    for mask, coefficient in vector.items():
        for i in range(16 - n):
            moved = fermion_move(mask, i, i + n)
            if moved is not None:
                state, sign = moved
                factor = F(i - 8) + alpha + F(n, 2)
                out[state] = out.get(state, F(0)) + coefficient * sign * factor
    return {key: value for key, value in out.items() if value}


def boson_stress_action(vector, n, p):
    """Independent c1 Sugawara action in the one-current PBW basis."""
    out = {}
    def add(key, value):
        key = tuple(sorted(key, reverse=True))
        out[key] = out.get(key, F(0)) + value
    for key, coefficient in vector.items():
        add(key + (n,), coefficient * p)
        for k in range(1, n):
            add(key + (k, n - k), coefficient / 2)
        for k, count in Counter(key).items():
            remainder = list(key)
            remainder.remove(k)
            add(tuple(remainder) + (k + n,), coefficient * k * count)
    return {key: value for key, value in out.items() if value}


def rank_mod_prime(matrix, prime):
    rows = [[int(v.p) * pow(int(v.q), -1, prime) % prime for v in row]
            for row in matrix.tolist()]
    rank = 0
    for column in range(matrix.cols):
        pivot = next((i for i in range(rank, matrix.rows) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        rows[rank] = [v * inverse % prime for v in rows[rank]]
        for i in range(rank + 1, matrix.rows):
            factor = rows[i][column]
            rows[i] = [(a - factor * b) % prime for a, b in zip(rows[i], rows[rank])]
        rank += 1
    return rank


def reconstruct_in_original_stress(quarter, source_vector, fold, d_standard):
    alpha = F(quarter)
    sea = (1 << 8) - 1
    car_basis = []
    for occupied in combinations(range(16), 8):
        if sum(i - 8 for i in occupied) - sum(range(-8, 0)) == 8:
            car_basis.append(sum(1 << i for i in occupied))
    words = partitions(8)
    need(len(car_basis) == len(words) == 22, 'actual grade8 dimensions')
    car_columns, boson_columns = [], []
    for word in words:
        car, boson = {sea: F(1)}, {(): F(1)}
        for n in reversed(word):
            car = original_stress_action(car, n, alpha)
            boson = boson_stress_action(boson, n, alpha - F(1, 2))
        need(set(car) <= set(car_basis), 'energy escapes the exact grade8 sector')
        car_columns.append(sp.Matrix([sp.Rational(car.get(mask, F(0))) for mask in car_basis]))
        boson_columns.append(sp.Matrix([sp.Rational(boson.get(key, F(0))) for key in words]))
    car_matrix = sp.Matrix.hstack(*car_columns)
    boson_matrix = sp.Matrix.hstack(*boson_columns)
    boson_norms = []
    for key in words:
        norm = 1
        for k, multiplicity in Counter(key).items():
            norm *= k**multiplicity * factorial(multiplicity)
        boson_norms.append(norm)
    need(car_matrix.T * car_matrix ==
         boson_matrix.T * sp.diag(*boson_norms) * boson_matrix,
         'original CAR and independent c1 Sugawara Gram matrices differ')
    target_source = fold.T * sp.Matrix(source_vector)
    actual_delta = sp.Matrix([sp.Rational(F(k) + alpha, 8) - sp.Rational(1, 2)
                              for k in range(8)])
    need(fold.T * sp.Matrix(d_standard) == actual_delta, 'original source chart changed')
    target, folded_stress_target = {}, {}
    for k in range(8):
        mask, sign = fermion_move(sea, k, k + 8)
        target[mask] = sign * target_source[k]
        folded_stress_target[mask] = sign * actual_delta[k]
    target_vector = sp.Matrix([target.get(mask, 0) for mask in car_basis])
    need(car_columns[0] / 8 == sp.Matrix([folded_stress_target.get(mask, 0)
                                         for mask in car_basis]),
         'Lhat_-1=Lcal_-8/8 fails on the actual sea')
    solution = car_matrix.inv(method='DM') * target_vector
    need(car_matrix * solution == target_vector, 'exact original stress reconstruction')
    need((target_vector.T * target_vector)[0] == sum(x*x for x in source_vector),
         'CAR chart lost the W8 response norm')
    ranks = {str(prime): rank_mod_prime(car_matrix, prime)
             for prime in (65537, 65539, 1000003)}
    need(all(rank == 22 for rank in ranks.values()), 'original stress not full rank')
    determinant = car_matrix.det(method='domain-ge')
    need(determinant != 0, 'zero original stress determinant')
    return {'original_c': 1, 'original_momentum': str(alpha-F(1, 2)),
            'original_ground_weight': '1/32', 'source_grade': 8,
            'neutral_car_basis_masks': car_basis, 'dimension': 22,
            'literal_CAR_and_boson_Sugawara_Grams_equal': True,
            'rank_certificates_mod_prime': ranks,
            'CAR_coefficient_matrix_determinant': str(determinant),
            'W8_first_response_in_source_channel_chart': list(map(str, target_source)),
            'stress_words': [list(word) for word in words],
            'exact_W8_first_response_coefficients': list(map(str, solution)),
            'reconstruction_residual': '0', 'same_state_time_and_adjoint': True}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source-dir', type=Path, default=Path(
        '/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/'
        'original-admissible-source-closure-20261025'))
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    inputs = {name: (args.source_dir / name).read_bytes() for name in (
        'w8-source-twopoint-results.json', 'quantum-octic-completion-results.json')}
    fold_path = args.source_dir.parent / 'v16-operator-pullback-20261014/fold-results.json'
    fold_bytes = fold_path.read_bytes()
    fold = sp.Matrix(json.loads(fold_bytes)['explicit_Weyl_matrix'])
    need(fold.T * fold == sp.eye(8), 'old M is not orthogonal')
    old = json.loads(inputs['w8-source-twopoint-results.json'])
    completion = json.loads(inputs['quantum-octic-completion-results.json'])
    coefficients = list(map(F, completion['coefficients']))
    pairs = completion['basis_after_H8']
    need(len(coefficients) == 17 and len(pairs) == 16 and coefficients[0] == 1,
         'the declared old W8 ansatz changed')
    roots = roots_e8()
    cases = {}
    for quarter, case in old['cases'].items():
        d = tuple(map(F, case['d_standard']))
        d2 = dot(d, d)
        h8 = sum((dot(a, d)**8 for a in roots), F(0)) - 30 * d2**4
        gradient = tuple(sum((8 * dot(a, d)**7 * a[j] for a in roots), F(0))
                         - 240 * d2**3 * d[j] for j in range(8))
        need(dot(d, gradient) == 8 * h8, 'Euler identity for H8')
        onepoint, radial_correction = h8, F(0)
        for coefficient, monomial in zip(coefficients[1:], pairs):
            constant, first = contraction_constant_and_first(monomial, d)
            onepoint += coefficient * constant
            radial_correction += coefficient * first
        v1 = tuple(x + radial_correction * y for x, y in zip(gradient, d))
        need(onepoint == F(case['onepoint']), 'onepoint differs from sparse PBW')
        need(v1 == tuple(map(F, case['level_1_vector'])), 'v1 differs from sparse PBW')
        need(dot(d, v1) == 8 * onepoint, 'primary Ward identity')
        norm = dot(v1, v1)
        need(norm == F(case['level_1_norm_squared']), 'first spectral norm')
        perpendicular = tangential(v1, d)
        need(perpendicular == tangential(gradient, d), 'intrinsic tangent gradient')
        need(dot(perpendicular, d) == 0, 'orthogonality')
        residual = dot(perpendicular, perpendicular)
        parallel = dot(d, v1)**2 / d2
        need(norm == parallel + residual, 'orthogonal norm decomposition')
        need(residual == F(51311452775625, 9345848836096) > 0,
             'the expected nonzero transverse response')
        need(residual == F(case['level_1_excess_over_cauchy_bound']),
             'old Cauchy excess differs')
        # A purely radial comparison polynomial has no such component.
        radial_control = tuple(8 * d2**3 * x for x in d)
        need(tangential(radial_control, d) == (F(0),) * 8, 'radial control')
        cases[quarter] = {
            'd_squared': str(d2), 'h_d': str(d2 / 2),
            'H8_value': str(h8), 'H8_gradient': list(map(str, gradient)),
            'all_16_Q_corrections_radial': True,
            'quantum_radial_correction': str(radial_correction),
            'onepoint': str(onepoint), 'level_one_vector': list(map(str, v1)),
            'level_one_norm_squared': str(norm),
            'folded_c8_Virasoro_level_one_space': 'span{Lhat_-1^d Omega_d}=span{d.h_-1 Omega_d}',
            'folded_c8_projection_norm_squared': str(parallel),
            'transverse_vector': list(map(str, perpendicular)),
            'transverse_norm_squared': str(residual),
            'transverse_fraction_of_level_one_norm': str(residual / norm),
            'not_a_TFPT_completion_fraction': True,
            'full_original_c1_stress_reconstruction': reconstruct_in_original_stress(
                quarter, v1, fold, d),
        }
    need(cases['1/4']['transverse_norm_squared'] == cases['3/4']['transverse_norm_squared'],
         'quarter norm comparison')
    result = {
        'claim': 'MATH.ORIGINAL-SEAM.W8-ORIGINAL-STRESS-RECONSTRUCTION.01',
        'status': 'EXACT_ORIGINAL_C1_RECONSTRUCTION_AND_SCOPED_FOLDED_C8_EXCLUSION',
        'input_sha256': {name: sha256(data).hexdigest() for name, data in inputs.items()},
        'original_fold_sha256': sha256(fold_bytes).hexdigest(),
        'root_count': len(roots), 'cases': cases,
        'analytic_consequence': 'The full ORIGINAL c1 Virasoro cycle equals the entire fixed '
            'CAR charge sector because h=1/32 is not a c1 Kac weight m^2/4. All W8 insertion '
            'levels and their positive-semigroup-smoothed response are in its closure. Only '
            'the smaller folded c8 cycle misses the tangential first-level direction.',
        'scope': 'No local P1 operator identification or single implemented stress pulse. '
            'No new source or state, no all-P1 exclusion, no absolute physical time calibration.'}
    encoded = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(encoded)
    else:
        print(encoded, end='')


if __name__ == '__main__':
    main()
