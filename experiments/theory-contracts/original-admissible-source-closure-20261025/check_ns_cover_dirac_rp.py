#!/usr/bin/env python3
"""Exact checks of the actual NS cover Dirac Green calculation.

The infinite dimensional positivity proof is in NS-COVER-FULL-DIRAC-RP.txt.
No finite calculation is used to infer an unrestricted spectral theorem.
"""

import argparse
import hashlib
import json
from pathlib import Path

import sympy as s


def require(value, message):
    if not value:
        raise ValueError(message)


def zero(expr):
    return s.simplify(expr) == 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    r, t, u = s.symbols('r t u', positive=True, real=True)
    # r=exp(E), beta=2.  All formulas retain both signs of E.
    g_positive = r**(-t) / (1 + r**-2)
    g_negative = -r**(-(t + 2)) / (1 + r**-2)
    require(zero(g_positive.subs(t, 0) - g_negative.subs(t, 0) - 1), 'jump')
    require(zero(g_negative + g_positive.subs(t, t + 2)), 'NS period')
    k = -g_negative.subs(t, -t-u)
    bt = r**t / s.sqrt(1+r**2)
    bu = r**u / s.sqrt(1+r**2)
    require(zero(k - bt*bu), 'all real energy Gram identity')
    btilde_t = r**(1-t) / s.sqrt(1+r**2)
    btilde_u = r**(1-u) / s.sqrt(1+r**2)
    require(zero(g_positive.subs(t, t+u) - btilde_t*btilde_u),
            'conjugate charge all real energy Gram identity')

    # Spectral finite Gram controls have rational entries by using integer
    # quarter-time powers of a=exp(E/4), rather than rounding exponentials.
    modes = [s.Rational(1, 3), s.Rational(1, 2), s.Integer(1), s.Integer(2), s.Integer(3)]
    times = [1, 2, 3]  # times t_i/4 in (0,1)
    gram = s.Matrix([[sum(a**(i+j)/(1+a**8) for a in modes)
                      for j in times] for i in times])
    minors = [s.factor(gram[:n, :n].det()) for n in range(1, 4)]
    require(all(m > 0 for m in minors), 'positive exact finite Gram')

    sx = s.Matrix([[0, 1], [1, 0]])
    sz = s.diag(1, -1)
    w = s.symbols('w', real=True)
    spatial = []
    for b in (-1, 1):
        v = s.Matrix([1, b]) / s.sqrt(2)
        require(sx*v == b*v, 'zero boundary vector')
        require((v.conjugate().T*sz*v)[0] == 0, 'isotropic endpoint')
        for n in range(-4, 5):
            e = n*s.pi/2
            phi = s.diag(s.exp(-s.I*e*(w+1)), s.exp(s.I*e*(w+1)))*v
            require(all(zero(x) for x in s.I*sz*phi.diff(w)-e*phi), 'spatial eigenvector')
            for endpoint in (-1, 1):
                p = phi.subs(w, endpoint)
                require(all(zero(x) for x in sx*p-b*p), 'endpoint line')
            reflected = sx*phi.subs(w, -w)
            j = b*(-1)**n if n >= 0 else b*(-1)**(-n)
            require(all(zero(x) for x in reflected-j*phi), 'transverse parity')
            spatial.append({'b': b, 'n': n, 'energy': str(e), 'J_eigenvalue': j})

    # At a pair of deck sectors A has eigenvalues +i,-i.  Linear Clifford
    # reflection exchanges them; antiunitary KRx preserves their complex
    # eigenspaces because K conjugates the eigenvalue at the same time.
    A = s.diag(s.I, -s.I)
    reflection = sx
    p = (s.eye(2)-s.I*A)/2
    q = (s.eye(2)+s.I*A)/2
    require(A*A == -s.eye(2), 'deck square')
    require(p*p == p and q*q == q and p*q == s.zeros(2), 'Pin projectors')
    require(reflection*A*reflection == -A, 'linear reflection deck action')
    require(reflection*p*reflection == q, 'linear reflection exchange')
    # Here KRx on the sector labels is K: its action on A is bar(A)=-A.
    require(p.conjugate() == p and q.conjugate() == q, 'antiunitary sector preservation')

    plus = lambda h: s.simplify(h/(h-1))
    minus = lambda h: s.simplify(1/(h-1))
    pair = {
        '0_to_1': s.simplify((plus(s.I)+plus(-s.I))/2),
        '1_to_2': s.simplify((s.I*plus(s.I)-s.I*plus(-s.I))/2),
        'minus_1_to_0': s.simplify((minus(s.I)+minus(-s.I))/2),
        'minus_2_to_minus_1': s.simplify((-s.I*minus(s.I)+s.I*minus(-s.I))/2),
    }
    require(pair['0_to_1'] == pair['1_to_2'] == s.Rational(1, 2), 'positive zero pairing')
    require(pair['minus_1_to_0'] == pair['minus_2_to_minus_1'] == -s.Rational(1, 2), 'negative zero pairing')

    data = {
        'scope': 'given free charged Green field on the actual NS orientation cover',
        'beta_in_core_length_units': 2,
        'all_energy_green_jump': True,
        'all_energy_NS_antiperiodicity': True,
        'all_energy_reflected_gram_factorization': True,
        'conjugate_charge_all_energy_reflected_gram_factorization': True,
        'exact_positive_gram_leading_minors': [str(m) for m in minors],
        'transverse_mode_controls': spatial,
        'linear_clifford_reflection_exchanges_pin_lifts': True,
        'antiunitary_KRx_preserves_pin_lifts': True,
        'zero_mode_kernel_pairing': {k: str(v) for k, v in pair.items()},
        'thermal_zero_occupation': '1/2',
        'thermal_covariance_is_not_a_projection': True,
        'physical_gates_closed': [],
        'complete_TFPT_solution': False,
    }
    body = json.dumps(data, indent=2, sort_keys=True)+'\n'
    args.output.write_text(body)
    print(hashlib.sha256(body.encode()).hexdigest())


if __name__ == '__main__':
    main()
