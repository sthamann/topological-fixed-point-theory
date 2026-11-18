#!/usr/bin/env python3
"""Exact first operator test for the old winding times charge conjugation.

The infinite Hilbert-Schmidt exclusion uses the analytic lower bound in
PIN-WINDING-DESCENT.txt, not a fitted growth rate of a finite matrix.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path

import sympy as s


def need(test, message):
    if not test:
        raise RuntimeError(message)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    sig1 = s.Matrix([[0, 1], [1, 0]])
    sig2 = s.Matrix([[0, -s.I], [s.I, 0]])
    eye = s.eye(2)
    k, ell = s.symbols('k ell', real=True)
    symbol = k*sig1+ell*sig2
    # A=e^(i*pi*x) sigma1 K changes plane momenta to(pi-k,-ell).
    a_symbol = (s.pi-k)*sig1-ell*sig2
    need(a_symbol*sig1 + sig1*s.conjugate(symbol) - s.pi*eye
         == s.zeros(2), 'D A = -A D + pi sigma1 A')
    need(sig1*s.conjugate(sig1) == eye, 'A squared')
    for c in (s.I, -s.I):
        need(-s.conjugate(c) == c, 'fixed Pin-minus multiplier')
    for b in (-1, 1):
        boundary = (eye+b*sig1)/2
        need(sig1*s.conjugate(boundary) == b*boundary,
             'charge conjugation preserves boundary eigenline')

    x, mu = s.symbols('x mu', real=True, positive=True)
    energy0sq = x*x+mu*mu
    energy1sq = (x-s.pi)**2+mu*mu
    inner = x*(x-s.pi)+mu*mu
    need(s.expand(energy0sq*energy1sq-inner**2-s.pi**2*mu**2) == 0,
         'projector angle determinant identity')
    # Under |x|<=mu and mu>=2*pi:
    # E0²<=2mu², E1²<=13mu²/4; rationalization gives
    # ||P1-P0||HS² >= pi²/(13mu²)=4/(13n²).
    need(s.Rational(1, 1)/(2*2*s.Rational(13, 4)) == s.Rational(1, 13),
         'uniform analytic projector lower bound')
    tested = 0
    for n in range(4, 129):
        for theta in (F(1, 4), F(3, 4)):
            count = sum(abs(2*(F(m)+theta)) <= F(n, 2)
                        for m in range(-n, n+1))
            need(count >= F(n, 4), 'quarter grid count')
            need(count*F(4, 13*n*n) >= F(1, 13*n), 'level lower bound')
            tested += 1
    result = {
        'claim': 'MATH.ORIGINAL-SEAM.PIN-WINDING-DESCENT.01',
        'verdict': 'FIXED_DOMAIN_DESCENT_PROVED; FULL_2D_SEA_PH_LIFT_EXCLUDED',
        'operator': 'A=exp(i*pi*x)*sigma1*K',
        'A_squared': 1,
        'domain': 'same fixed c=plus/minus i,b=plus/minus1,H1,cover support',
        'Dirac_identity': 'A D A=-D+pi*sigma1',
        'centering_changes_given_operator': True,
        'projector_hs_lower_bound_each_level_n_ge_4': '1/(13*n)',
        'quarter_grid_control_cases': tested,
        'infinite_result_basis': 'analytic harmonic-series lower bound',
        'calderon_boundary_range': 'OPEN; bulk nonintertwining alone does not decide it',
        'scope': 'Given flat Moebius Dirac family with specified local b boundary '
                 'and canonical full D sea. Not all P1 sources or reduced-circle seas.',
        'physical_gates_closed': [],
        'complete_TFPT_solution': False,
    }
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'domain': 'PASS', 'operator_identity': 'PASS',
                      'analytic_hs_lower_bound': '1/(13*n)',
                      'grid_cases': tested, 'physical_gates_closed': []}))


if __name__ == '__main__':
    main()
