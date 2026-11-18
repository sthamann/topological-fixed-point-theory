"""Executable inference controls on existing TFPT source functions.

These check explicitly named consequences, not whole origin hypotheses.
Analytic identities are recorded separately from finite precision replays.
"""
import argparse
import json
import mpmath as mp


def evaluate(check_id):
    mp.mp.dps = 65
    theta = lambda v, t: mp.jtheta(3, mp.pi*v, mp.exp(-mp.pi*t))
    theta0 = lambda t: theta(0, t)
    close = lambda a, b: abs(a-b) < mp.mpf('1e-55')
    if check_id == 'quillen-relative-normalization':
        # Same geometry and reference on numerator/denominator: N(0,t)=1.
        norm = lambda v, t: (theta(v, t)/theta0(t))**2
        dv2 = mp.diff(lambda v: mp.log(norm(v, 1)), 0, 2)
        dt = mp.diff(lambda t: mp.log(norm(0, t)), 1)
        if not close(dt, 0) or not close(dv2, -2*mp.pi):
            raise ArithmeticError('Original relative response replay failed')
        facts = {'shape_derivative_log_N_at_zero_source': str(dt),
                 'holonomy_second_derivative_log_N': str(dv2),
                 'analytic_reason': 'N(0,t)=1 identically for the co-moving zero-source reference'}
        claims = {'normalized_zero_source_has_shape_response': False,
                  'normalized_charge_response_alone_fixes_absolute_replica_entropy': False}
        scope = 'The given flat elliptic relative function; t variation is not an identified physical replica operation.'
    elif check_id == 'theta-heat-shape-response':
        dt = mp.diff(lambda t: mp.log(theta0(t)), 1)
        ratio = mp.diff(lambda v: theta(v, 1), 0, 2)/theta0(1)
        if not close(dt, -mp.mpf(1)/4) or not close(ratio/(4*mp.pi), dt):
            raise ArithmeticError('Theta heat identity failed')
        facts = {'theta_shape_log_derivative': str(dt),
                 'theta_holonomy_second_derivative_ratio': str(ratio),
                 'analytic_reason': 'Poisson Z(t)=t^-1/2 Z(1/t) and theta heat equation'}
        claims = {'charge_curvature_fixes_theta_shape_derivative': True,
                  'theta_shape_derivative_is_already_full_master_area_entropy': False}
        scope = 'Unnormalized theta factor only; eta, metric, reference and sewing terms remain part of the actual determinant.'
    elif check_id == 'doubled-source-counting':
        beta = 4*mp.pi
        def variance(a):
            weights = [(q, mp.exp(-beta*(mp.mpf(q*q)/2+(a-mp.mpf(1)/2)*q))) for q in range(-16, 17)]
            z = sum(w for q, w in weights)
            mean = sum(q*w for q, w in weights)/z
            return sum((q-mean)**2*w for q, w in weights)/z
        variances = [variance(mp.mpf(1)/4), variance(mp.mpf(3)/4)]
        etas = [mp.mpf(1)-2*a for a in (mp.mpf(1)/4, mp.mpf(3)/4)]
        if not all(close(v, 1/(8*mp.pi)) for v in variances) or sum(etas) != 0:
            raise ArithmeticError('Quarter source identities failed')
        facts = {'quarter_variances': [str(v) for v in variances],
                 'pair_variance': str(sum(variances)),
                 'beta_times_each_variance': [str(beta*v) for v in variances],
                 'eta_quarter_pair': [str(x) for x in etas], 'eta_sum': str(sum(etas)),
                 'analytic_reason': 'Original quarter pair reduces to exp(-pi p^2); eta_a(0)=1-2a for n+a.'}
        claims = {'two_quarter_blocks_are_two_faces': False,
                  'nonzero_total_eta_is_needed_for_original_half_susceptibility': False,
                  'quarter_pair_variance_equals_one_quarter_variance': False}
        scope = 'Original independent quarter blocks on one marginal, not the two faces; eta statement is for the unweighted direct sum, not the whole physical master.'
    elif check_id == 'variable-mark-completion':
        # Candidate class, not a deduction of these Gram forms from P1.
        from fractions import Fraction
        import sympy as sp
        from sympy.matrices.normalforms import hermite_normal_form
        finite_replays = []
        for r in (1, 2, 3):
            n, g = 4*r*r, 4*r*r+1
            d = sp.zeros(g)
            for j in range(g-1):
                d[j, j], d[j+1, j] = 1, -1
            d[g-2, g-1], d[g-1, g-1] = 1, 1
            a = sp.zeros(n, n-1)
            for j in range(n-1):
                a[j, j], a[j+1, j] = 1, -1
            base = sp.diag(d, a)
            v = sp.Matrix([sp.Rational(1, 2)]*g +
                          [sp.Rational(r*(n-1), n)]+[sp.Rational(-r, n)]*(n-1))
            coordinates = base.gauss_jordan_solve(v)[0]
            den = 4*r
            integer_generators = (den*sp.eye(2*n)).row_join(den*coordinates)
            if any(x.q != 1 for x in integer_generators):
                raise ArithmeticError('Nonintegral presentation')
            integer_generators = sp.Matrix(integer_generators.rows, integer_generators.cols,
                [sp.Integer(x) for x in integer_generators])
            h = hermite_normal_form(integer_generators)
            gram = h.T*(base.T*base)*h/den**2
            if (gram.det() != 1 or any(x.q != 1 for x in gram)
                    or any(gram[j, j] % 2 for j in range(2*n))):
                raise ArithmeticError('Lattice completion replay failed')
            finite_replays.append({'marks': n, 'rank': 2*n, 'base': f'D{g}+A{n-1}',
                'glue_order': 4*r, 'generator_squared_norm': str((v.T*v)[0]),
                'gram_determinant': str(gram.det()), 'integral_even': True})
        facts = {'conditional_class': 'standard positive D_(n+1)+A_(n-1), complete even unimodular extension',
                 'all_admissible_mark_counts': 'n=4*r^2, integer r>=1',
                 'analytic_necessity': 'det(base)=4n must be a square; even unimodular rank 2n is divisible by8',
                 'analytic_sufficiency': 'glue generator (spinor_D,r*lambda1_A) has order4r and norm2r^2, yielding index4r and determinant1 for every r',
                 'finite_independent_hnf_replays': finite_replays,
                 'minimal_marks': 4, 'minimal_rank': 8,
                 'source_pairing_or_minimality_derived_from_P1': False,
                 'primitive_family_projection_condition': 'A local glue element with primitive A_(n-1) discriminant projection has order divisible by n; n divides |H|=4r, while n=4r^2, so r=1.',
                 'primitive_projection_selected_marks': 4,
                 'primitive_local_charge_derived_from_P1': False}
        claims = {'higher_rank_completion_forces_residual_discriminant_sectors': False,
                  'local_complete_gluing_alone_selects_four_marks': False,
                  'minimal_member_of_this_declared_class_has_four_marks': True,
                  'primitive_family_projection_allows_higher_marks': False}
        scope = 'An explicit candidate family with assumed standard pairing and even selfdual completion. Minimality is additional. No derivation of physical P1, spacetime, or all TFPT outputs.'
    else:
        raise ValueError('Unknown curated check ID')
    return {'id': check_id, 'status': 'pass', 'facts': facts,
            'inference_claims': [{'id': k, 'holds': v} for k, v in claims.items()],
            'scope': scope, 'precision_digits': 65,
            'new_physical_proof': False, 'fundamental_solution': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', required=True)
    args = parser.parse_args()
    print(json.dumps(evaluate(args.check), ensure_ascii=False, indent=2))
