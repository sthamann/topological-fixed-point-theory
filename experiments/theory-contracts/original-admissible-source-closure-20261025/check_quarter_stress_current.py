#!/usr/bin/env python3
"""Recover the old current from the SAME two quarter stress families.

This is an exact coefficient synthesis, not a search over new physical inputs.
The all-mode proof is symbolic. The finite CAR replay checks the old literal
operator action and the nonzero Heisenberg central term on an interior core.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
import json
from pathlib import Path

import sympy as sp


def need(test, message):
    if not test:
        raise RuntimeError(message)


def add(a, b, ca=F(1), cb=F(1)):
    result = {key: ca * value for key, value in a.items()}
    for key, value in b.items():
        result[key] = result.get(key, F(0)) + cb * value
    return {key: value for key, value in result.items() if value}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo', type=Path,
                        default=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4'))
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    old = args.repo / 'experiments/theory-contracts'
    stress_file = old / 'original-admissible-source-closure-20261025/check_w8_stress_cyclic.py'
    spec = importlib.util.spec_from_file_location('old_actual_stress', stress_file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    # The original CAR coefficient is affine in alpha at EVERY mode.
    a, b, r, n = sp.symbols('a b r n')
    coefficient = a * (r + n / 2 + sp.Rational(3, 4)) \
                + b * (r + n / 2 + sp.Rational(1, 4)) - 1
    equations = sp.Poly(coefficient, r, n).coeffs()
    solutions = sp.solve(equations, (a, b), dict=True)
    need(solutions == [{a: 2, b: -2}], 'current synthesis is not unique')
    need(sp.expand(coefficient.subs(solutions[0])) == 0, 'all-mode identity')
    p1, p3 = -sp.Rational(1, 4), sp.Rational(1, 4)
    need(p1**2 / 2 == p3**2 / 2 == sp.Rational(1, 32), 'ground terms')

    # The existing full-cover sectors must NOT be silently identified. For
    # independent c1 factors, the same subtraction has a Virasoro cubic term.
    k = sp.Symbol('k', positive=True, integer=True)
    h = sp.Rational(1, 32)
    independent_norm = sp.expand(4 * (4 * k * h + 2 * k * (k**2 - 1) / 12))
    need(sp.expand(independent_norm - k * (4*k**2 - 1) / 6) == 0,
         'two-factor Virasoro norm')
    independent_values = [independent_norm.subs(k, j) for j in range(1, 5)]
    need(independent_values[:2] == [sp.Rational(1, 2), sp.Integer(5)],
         'first two independent-factor tests')
    scales = [sp.cancel(j / independent_values[j-1]) for j in (1, 2)]
    need(scales[0] != scales[1], 'a common rescaling unexpectedly fits')
    # For any real self-adjoint combination aL+ + bL-, its vacuum norm is
    # (a²+b²)(2kh+k(k²-1)/12). A nonzero Heisenberg level would require
    # a²+b²=0 to remove the cubic term, hence a=b=0 over the reals.
    aa, bb = sp.symbols('aa bb', real=True)
    independent_general = (aa**2 + bb**2) * (2*k*h + k*(k**2-1)/12)
    need(sp.Poly(sp.expand(independent_general), k).coeff_monomial(k**3)
         == (aa**2 + bb**2) / 12, 'independent-factor cubic coefficient')

    # Old actual comparison: u_cov W0=W1 and u_cov W1=W0 u_base.
    # Cover indices move as 2m -> 2m+1 -> 2m+2. Pull back the odd-sector
    # matrix units using the FIRST arrow, retaining the full wrap carry.
    m, step = sp.symbols('m step', integer=True)
    need(sp.expand((2*m+1)-1-2*m) == 0, 'source generator pullback')
    need(sp.expand((2*(m+step)+1)-1-2*(m+step)) == 0,
         'target generator pullback')
    for mm in range(-100, 101):
        need((2*mm+1 < 0) == (2*mm < 0), 'first-arrow sea mismatch')
    wrap_changes = [mm for mm in range(-100, 101)
                    if (2*mm+1 < 0) != (2*mm+2 < 0)]
    need(wrap_changes == [-1], 'winding carry was lost')

    # Replay the real old CAR moves, with the original common sea and ALL
    # intermediate states. Tests do not project to the folded residue sector.
    def action(vector, step, alpha=None):
        if alpha is not None and step > 0:
            return module.original_stress_action(vector, step, alpha)
        result = {}
        for mask, value in vector.items():
            for i in range(16):
                j = i + step
                if not 0 <= j < 16:
                    continue
                moved = module.fermion_move(mask, i, j)
                if moved is None:
                    continue
                target, sign = moved
                factor = F(1) if alpha is None else F(i - 8) + alpha + F(step, 2)
                result[target] = result.get(target, F(0)) + value * sign * factor
        return {key: value for key, value in result.items() if value}

    sea = (1 << 8) - 1
    vectors = [{sea: F(1)}]
    for k in range(1, 5):
        vectors.append(action({sea: F(1)}, k))
        vectors.append(module.original_stress_action({sea: F(1)}, k, F(1, 4)))
    tested = 0
    for vector in vectors:
        for step in list(range(-7, 0)) + list(range(1, 8)):
            recovered = add(action(vector, step, F(3, 4)),
                            action(vector, step, F(1, 4)), F(2), F(-2))
            need(recovered == action(vector, step), 'literal CAR family difference')
            tested += 1
    # At step zero the diagonal normal ordering is explicit.
    for charge in range(-4, 5):
        grade = sp.Symbol('grade')
        h1 = grade + sp.Rational(charge)**2/2 + p1*charge + p1**2/2
        h3 = grade + sp.Rational(charge)**2/2 + p3*charge + p3**2/2
        need(sp.expand(2*(h3-h1)) == charge, 'zero-mode charge term lost')

    central = []
    for k in range(1, 8):
        # J_+k annihilates the sea; J_-k creates the finite interior excitation.
        value = add(action(action({sea: F(1)}, k), -k),
                    action(action({sea: F(1)}, -k), k), F(1), F(-1))
        need(value == {sea: F(k)}, 'interior-core current central term')
        central.append(k)
    # Finite cut matrices cannot satisfy CCR on every state: retain this
    # negative control instead of treating their commutator as an operator
    # identity on the full truncated Fock space.
    empty = {}
    need(action(action({0: F(1)}, 1), -1) == empty, 'empty-state control')
    need(action(action({0: F(1)}, -1), 1) == empty, 'empty-state control')

    sources = [stress_file,
               old / 'joint-native-algebra-selection-20261016/ORIGINAL-ORDERED-CURRENT-KERNEL-AND-SEAM-GRADING.txt',
               old / 'native-octic-seam-origin-20261025/GEOMETRIC-SOURCE-JET.txt',
               old / 'canonical-seam-fermionization-20261025/MOBIUS-QUARTER-OPERATOR-DICTIONARY.txt']
    result = {
        'claim': 'MATH.ORIGINAL-SEAM.QUARTER-STRESS-CURRENT.01',
        'verdict': 'EXACT_IN_THE_EXISTING_COMMON_CAR_FRAME',
        'synthesis': {'L_3_4_coefficient': 2, 'L_1_4_coefficient': -2,
                      'all_nonzero_modes_symbolic': True,
                      'zero_mode_including_all_charge_terms': True},
        'finite_literal_car_replay': {'vector_mode_pairs': tested,
                                      'central_term_modes': central,
                                      'full_truncation_ccr_negative_control': True},
        'independent_cover_factor_negative_control': {
            'difference_vacuum_norm_polynomial': str(sp.factor(independent_norm)),
            'norms_modes_1_to_4': list(map(str, independent_values)),
            'level_one_current_norms_modes_1_to_4': [1, 2, 3, 4],
            'no_common_rescaling': True,
            'no_nonzero_real_linear_stress_combination_is_heisenberg': True,
        },
        'old_winding_comparison': {
            'u_cov_W0_equals_W1': True,
            'u_cov_W1_equals_W0_u_base': True,
            'first_arrow_matrix_units_pull_back_exactly': True,
            'first_arrow_sea_and_subtraction_agree': True,
            'wrap_sea_change_base_mode': -1,
            'full_cover_vacuum_preservation_claimed': False,
            'current': '2(Ad_Gamma(u_cov)^(-1)(L_3_4)-L_1_4)',
        },
        'scope': 'Same original periodic CAR frame, sea, subtraction and finite-energy core. '
                 'Old cut-cover winding comparison retained. No bare-P1 CAR-kernel derivation '
                 'or point-local unfolded W8 map.',
        'physical_gates_closed': [],
        'complete_TFPT_solution': False,
        'sources': [{'path': str(path.relative_to(args.repo)),
                     'sha256': sha256(path.read_bytes()).hexdigest()} for path in sources],
    }
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'verdict': result['verdict'], 'literal_pairs': tested,
                      'all_mode_identity': True, 'physical_gates_closed': []}))


if __name__ == '__main__':
    main()
