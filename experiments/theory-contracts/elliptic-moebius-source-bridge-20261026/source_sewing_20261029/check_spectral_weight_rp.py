"""Decisive countercheck of a spectral-weight RP inference, not a TFPT model.

The analytic sign certificate uses rational inequalities. Matrix exponentials
only give an independent numerical reproduction of the explicit finite case.
No original verification, paper, ledger, or physical source is modified.
"""
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def calculate():
    q = sp.symbols('q', positive=True)
    a, g, z = sp.symbols('a g z', real=True)
    D_same = sp.Matrix([[a, g], [g, a]])
    D_opposite = sp.Matrix([[a, g], [g, -a]])
    require(sp.expand((z*sp.eye(2)-D_same).det()
                      - ((z-a)**2-g**2)) == 0, 'same spectrum')
    require(D_opposite**2 == (a*a+g*g)*sp.eye(2), 'opposite square')
    S_same = q + q**9
    S_opposite = 2*q**5
    delta = q*(1-q**4)**2
    require(sp.expand(S_same-S_opposite-delta) == 0, 'action difference')
    base = sp.Matrix([[10, 9], [9, 10]])
    require(base.eigenvals() == {19: 1, 1: 1}, 'strict initial RP')
    reflection = sp.Matrix([[0, 1], [1, 0]])
    gaps = []
    for x in [-2, 2]:
        for y in [-2, 2]:
            D = sp.Matrix([[x, 1], [1, y]])
            swapped = sp.Matrix([[y, 1], [1, x]])
            require(reflection*D*reflection == swapped, 'reflection covariance')
            squares = [(e*e).simplify() for e in D.eigenvals()]
            require(all(s >= 1 for s in squares), 'nonzero gap')
            gaps.append({'configuration': [x, y],
                         'eigenvalues': sorted(str(v) for v in D.eigenvals())})

    # 2 < e < 3 follows from its positive series and a geometric tail bound.
    # Hence 1/3 < q=e^-1 < 1/2 and delta > (1/3)(15/16)^2.
    e_upper = sp.Rational(5, 2) + sp.Rational(1, 6)/(1-sp.Rational(1, 4))
    lower_delta = sp.Rational(1, 3)*(1-sp.Rational(1, 2)**4)**2
    require(e_upper == sp.Rational(49, 18) and e_upper < 3, 'e upper bound')
    require(sp.Rational(5, 2) > 2, 'e lower bound')
    require(lower_delta == sp.Rational(75, 256), 'delta rational lower bound')
    require(lower_delta > sp.Rational(1, 9), 'delta exceeds log threshold bound')
    # log(1+x)<x at x=1/9; thus delta>log(10/9) analytically.

    mp.mp.dps = 80
    weights = mp.matrix(2, 2)
    actions = mp.matrix(2, 2)
    signs = [-1, 1]
    for i, s in enumerate(signs):
        for j, t in enumerate(signs):
            D = mp.matrix([[2*s, 1], [1, 2*t]])
            heat = mp.expm(-(D*D))
            action = heat[0, 0]+heat[1, 1]
            actions[i, j] = action
            weights[i, j] = mp.exp(-action)
    same = mp.exp(-1)+mp.exp(-9)
    opposite = 2*mp.exp(-5)
    require(abs(actions[0, 0]-same) < mp.mpf('1e-70'), 'independent same action')
    require(abs(actions[0, 1]-opposite) < mp.mpf('1e-70'), 'independent opposite action')
    kernel = mp.matrix([[10*weights[0, 0], 9*weights[0, 1]],
                        [9*weights[1, 0], 10*weights[1, 1]]])
    normalizer = sum(kernel[i, j] for i in range(2) for j in range(2))
    direct_response = sum(signs[i]*signs[j]*kernel[i, j]
                          for i in range(2) for j in range(2))/normalizer
    eigen_minus = 10*mp.exp(-same)-9*mp.exp(-opposite)
    formula_response = eigen_minus/(10*mp.exp(-same)+9*mp.exp(-opposite))
    require(normalizer > 0 and min(weights) > 0, 'ordinary positive weights')
    require(abs(direct_response-formula_response) < mp.mpf('1e-70'), 'direct observable')
    require(direct_response < 0, 'RP violation')

    # A source-independent reference denominator cancels in a normalized state.
    x, y, ref = sp.symbols('x y ref', positive=True)
    c0, c1, c2, c12 = sp.symbols('c0 c1 c2 c12', real=True)
    Z = c0+c1*x+c2*y+c12*x*y
    require(sp.simplify((Z/ref)/(c0/ref)-Z/c0) == 0, 'reference normalization')
    require(sp.simplify(sp.diff(sp.log(Z/ref), x, y)
                         - sp.diff(sp.log(Z), x, y)) == 0, 'connected cross response')

    return {
        'verdict': 'POINTWISE_POSITIVE_SPECTRAL_WEIGHT_DOES_NOT_PRESERVE_RP',
        'exact_profile': 'f(u)=exp(-u)',
        'operator': 'D(x,y)=[[x,1],[1,y]], x,y in {-2,2}',
        'base_reflection_kernel': [[10, 9], [9, 10]],
        'base_eigenvalues': [19, 1],
        'spectra': gaps,
        'action_difference': str(delta),
        'analytic_sign_certificate': 'Delta>75/256>1/9>log(10/9)',
        'action_same': mp.nstr(same, 40),
        'action_opposite': mp.nstr(opposite, 40),
        'negative_kernel_eigenvalue': mp.nstr(eigen_minus, 40),
        'normalized_reflected_response': mp.nstr(direct_response, 40),
        'base_normalized_reflected_response': '1/19',
        'relative_reference_cancels_in_normalized_state': True,
        'general_profile_qualification':
            'Every nonconstant finite positive Laplace profile has Delta>0; '
            'a general counterexample uses r in (exp(-Delta),1). The fixed '
            'ratio 9/10 needs Delta>log(10/9).',
        'scope': 'Generic inference in the old measure proof only. This is '
                 'not an original P1/P2 source, a physical TFPT prediction, '
                 'or a continuum geometric/BRST countermodel.',
        'complete_TFPT_solution': False,
        'ledger_changed': False,
    }


if __name__ == '__main__':
    result = calculate()
    destination = Path(__file__).with_name('spectral-weight-rp-results.json')
    destination.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(result['verdict'])
    print('Reflected response:', result['normalized_reflected_response'])
    print('Sign proof:', result['analytic_sign_certificate'])
    print('Physical P1/P2 source closure: not established')
