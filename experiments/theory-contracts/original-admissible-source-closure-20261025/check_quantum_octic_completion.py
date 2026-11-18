#!/usr/bin/env python3
"""Complete the supplied octic in the existing c=8 current algebra.

Exact finite polynomial identities, not a new source model.  The search
space is explicitly H8 plus the sixteen quadratic-contraction monomials
of conformal weight eight.  Uniqueness is asserted in this space only.
"""

from collections import Counter
from itertools import combinations, combinations_with_replacement, product
from math import factorial, prod
from pathlib import Path
import importlib.util
import json
import sys

import sympy as s

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))


def need(condition, message):
    if not condition:
        raise ArithmeticError(message)


def main():
    h = {n: s.symbols(f'h{n}_0:8') for n in range(1, 8)}
    variables = tuple(v for n in range(1, 8) for v in h[n])
    x = h[1]
    # Sum over the 128 half roots: sign parity kills every exponent
    # pattern except all-even exponents and (1,1,1,1,1,1,1,1).
    p8 = s.S.Zero
    for indices in combinations_with_replacement(range(8), 4):
        powers = Counter(indices)
        coefficient = s.Rational(factorial(8), 2 * prod(
            factorial(2 * powers[j]) for j in range(8)))
        p8 += coefficient * prod(x[j] ** (2 * powers[j]) for j in range(8))
    p8 += 20160 * prod(x)
    # The remaining 112 integer roots are +-e_i +-e_j.
    for i, j in combinations(range(8), 2):
        for k in (0, 2, 4, 6, 8):
            p8 += 4 * s.binomial(8, k) * x[i] ** (8-k) * x[j] ** k
    q = sum(v*v for v in x)
    H8 = s.expand(p8 - 30*q**4)
    need(s.expand(sum(s.diff(H8, v, 2) for v in x)) == 0,
         'primitive root moment is not harmonic')

    Q = {(a, b): sum(h[a][j]*h[b][j] for j in range(8))
         for a in range(1, 8) for b in range(a, 8) if a+b <= 8}
    keys = sorted(Q, key=lambda k: (sum(k), k))
    monomials = [ks for n in range(1, 5)
                 for ks in combinations_with_replacement(keys, n)
                 if sum(a+b for a, b in ks) == 8]
    need(len(monomials) == 16, 'the declared ansatz changed')
    basis = [H8] + [s.expand(prod(Q[k] for k in ks)) for ks in monomials]
    w = s.Symbol('w')
    substitution = {h[n][j]: h[n][j] + (1 if j == 1 else -1)*w**n
                    for n in range(1, 8) for j in (0, 1)}
    exponential = [s.S.One]
    for n in range(1, 8):
        exponential.append(s.expand(sum(
            (h[m][0]-h[m][1])*exponential[n-m] for m in range(1, n+1))/n))

    def screening(polynomial):
        shifted = s.Poly(s.expand(polynomial.subs(substitution, simultaneous=True)), w)
        return s.expand(sum(shifted.nth(r)*exponential[r-1] for r in range(1, 9)))

    def vir(polynomial, k):
        linear = sum((n+k)*h[n][j]*s.diff(polynomial, h[n+k][j])
                     for n in range(1, 8-k) for j in range(8))
        quadratic = sum(s.Rational(m*(k-m), 2)*s.diff(
            polynomial, h[m][j], h[k-m][j])
            for m in range(1, k) for j in range(8) if m < 8 and k-m < 8)
        return s.expand(linear + quadratic)

    rows = {}
    for column, polynomial in enumerate(basis):
        for label, expression in (('screen', screening(polynomial)),
                                  ('L1', vir(polynomial, 1)),
                                  ('L2', vir(polynomial, 2))):
            for mon, coefficient in s.Poly(expression, *variables).terms():
                if coefficient:
                    rows.setdefault((label, mon), {})[column] = coefficient
    entries = {(i, j): value for i, row in enumerate(rows.values())
               for j, value in row.items()}
    matrix = s.SparseMatrix(len(rows), len(basis), entries)
    nullspace = matrix.to_DM().nullspace().to_Matrix()
    need(nullspace.rows == 1 and nullspace[0, 0] != 0,
         'completion is absent or not unique in the declared ansatz')
    coefficients = [s.cancel(nullspace[0, j]/nullspace[0, 0]) for j in range(17)]
    expected = [1, *[s.Rational(v, 1333) for v in (
        28800, -100800, 201600, -126000, -80640, 201600, -130368,
        0, -12096, 34944, -46368, 21924, 40656, -33516, 3024, -2178)]]
    need(coefficients == expected, 'exact completion coefficients changed')
    state = s.expand(sum(c*p for c, p in zip(coefficients, basis)))
    need(screening(state) == 0, 'actual root zero mode does not annihilate W8')
    need(all(vir(state, k) == 0 for k in range(1, 9)),
         'a positive Virasoro mode does not annihilate W8')
    terms = s.Poly(state, *variables).terms()
    norm = s.cancel(sum(c*c*prod((1+i//8)**a*factorial(a)
                      for i, a in enumerate(mon)) for mon, c in terms))
    need(norm == s.Rational(65028096000, 43), 'positive full Fock norm changed')

    # Check the actual old lattice chart, not just an abstract E8 name.
    old_path = ROOT/'experiments/theory-contracts/native-octic-seam-origin-20261025/native_folded_octic_origin.py'
    spec = importlib.util.spec_from_file_location('old_octic_source', old_path)
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    lattice, O, U, M = old.original_transport()
    from tfpt_explorer import cartan_source as cs
    T = s.Matrix.hstack(*map(s.Matrix, cs._TARGET_SIMPLE))
    original_roots2 = {tuple(T*s.Matrix(r)) for r in lattice['root_coordinates']}
    standard_roots2 = set()
    for i, j in combinations(range(8), 2):
        for a, b in product((-2, 2), repeat=2):
            root = [0]*8
            root[i], root[j] = a, b
            standard_roots2.add(tuple(root))
    standard_roots2.update(r for r in product((-1, 1), repeat=8)
                           if sum(v < 0 for v in r) % 2 == 0)
    need(original_roots2 == standard_roots2 and len(original_roots2) == 240,
         'original root isometry no longer reaches this exact root set')

    z = s.symbols('z0:4')
    i8 = sum(v**8 for v in z) + 14*sum(z[a]**4*z[b]**4
          for a, b in combinations(range(4), 2)) + 168*prod(v*v for v in z)
    native_moment = s.expand(sum(
        sum((U.T*s.Matrix(r)/s.sqrt(2))[j]*z[j] for j in range(4))**8
        for r in lattice['roots']))
    need(s.expand(native_moment-s.Rational(15, 2)*i8) == 0,
         'the original Gaussian octic restriction was lost')
    native_norm = sum(c*c*prod(factorial(a) for a in mon)
                      for mon, c in s.Poly(i8, *z).terms())
    need(native_norm == 96**2*140, 'original octic norm changed')
    overlap_squared = s.Rational(15, 2)**2 * native_norm
    need(overlap_squared/norm == s.Rational(43, 896), 'native overlap changed')

    # The old complete S8 root OPE coefficient is not replaced by W8.
    # Pairing a state with its full root coherent exponential evaluates
    # each h_j(-n) at beta_j; all 22 partitions are retained by that identity.
    beta = [s.Integer(1), s.Integer(-1), *[s.Integer(0)]*6]
    evaluation = {h[n][j]: beta[j] for n in range(1, 8) for j in range(8)}
    overlap_s8 = s.cancel(240*state.subs(evaluation, simultaneous=True))
    need(overlap_s8 == s.Rational(48384000, 1333), 'full root OPE overlap changed')
    need(overlap_s8/norm == s.Rational(1, 41664), 'S8 component coefficient changed')

    # Actual original quarter folds: no direction, profile or coupling fit.
    source_cases = {}
    for alpha in (s.Rational(1, 4), s.Rational(3, 4)):
        delta = s.Matrix([(j+alpha)/8-s.Rational(1, 2) for j in range(8)])
        d = s.simplify(O.T*M*delta)
        length2 = s.simplify((d.T*d)[0])
        moment = s.simplify(sum(((s.Matrix(r).T*d)[0]/s.sqrt(2))**8
                               for r in lattice['roots']))
        # Y(h(-n),z)=partial^(n-1)h(z)/(n-1)!; at h(z)=d/z
        # every contraction monomial has total sign + at weight eight.
        charge = s.factor(moment-30*length2**4 + sum(
            c*length2**len(ks) for c, ks in zip(coefficients[1:], monomials)))
        need(charge == s.Rational(25123488075, 357824462848),
             'the given source one-point value vanished or changed')
        source_cases[str(alpha)] = {'d_squared': str(length2),
            'root_moment_8': str(moment), 'W8_onepoint_times_z8': str(charge)}

    print(json.dumps({
        'claim': 'MATH.ORIGINAL-SEAM.QUANTUM-OCTIC-COMPLETION.01',
        'status': 'PASS_EXACT_IN_GIVEN_SOURCE_AND_DECLARED_ANSATZ',
        'leading_harmonic': 'H8=sum_roots(alpha.h1)^8-30(h1.h1)^4',
        'basis_after_H8': [[list(k) for k in ks] for ks in monomials],
        'coefficients': [str(c) for c in coefficients],
        'constraint_rows': len(rows), 'ansatz_dimension': 17,
        'solution_dimension': nullspace.rows, 'constraint_rank': 17-nullspace.rows,
        'positive_Virasoro_modes_checked': list(range(1, 9)),
        'root_zero_mode': 'e_(e0-e1),0 W8=0; all roots follow by Weyl covariance',
        'old_root_chart_isometry': 'all 240 original roots match exactly',
        'full_PBW_terms': len(terms), 'norm_squared': str(norm),
        'normalized_native_F8_overlap_squared': str(overlap_squared/norm),
        'full_old_S8_overlap': str(overlap_s8),
        'coefficient_of_W8_in_old_S8': str(overlap_s8/norm),
        'normalized_old_S8_W8_overlap_squared': str(s.cancel(overlap_s8**2/(norm*15600))),
        'original_quarter_source_cases': source_cases,
        'globalization': 'Fixed by the connected compact finite-dimensional E8 zero-mode group, hence by K transitions and native inner clocks; conformal weight remains eight',
        'not_claimed': ['annihilation by all affine positive modes',
            'a physical projection or new prepared state',
            'P1 collar/reflection/reference identification',
            'a complete physical TFPT solution']
    }, indent=2)+'\n', end='')


if __name__ == '__main__':
    main()
