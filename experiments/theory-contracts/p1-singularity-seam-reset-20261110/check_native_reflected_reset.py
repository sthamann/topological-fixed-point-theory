"""Exact two-mode restriction of the existing thermal HS source, in even parity."""

import sympy as s


def main():
    # |00>, |11>: the ratio is exp[-beta*(E1+E2)/2] in the original source.
    t = s.Symbol('t', positive=True)
    d = s.diag(1, t)
    z = s.trace(d**2)
    rho = d**2 / z
    n = s.diag(0, 1)
    y = s.Matrix(2, 2, s.symbols('y:4'))
    e00, e11 = s.diag(1, 0), n
    e01 = s.Matrix([[0, 1], [0, 0]])
    e10 = e01.T
    probes = [e00, e11, e01 + e10, s.I * (e10 - e01)]
    rows = s.Matrix([s.trace(d * a * y) for a in probes])
    det = s.factor(rows.jacobian(list(y)).det())
    p = s.trace(rho * n)
    cov = s.factor(s.trace(d * n * d * n) / z - p**2)
    q = s.Symbol('q', positive=True)
    phase = s.diag(1, q)
    ht = phase * y * phase.inv()
    checks = {
        'four_reference_cross_probes_are_invertible': det == -2 * s.I * t**2,
        'both_reflected_occupation_marginals_match': s.trace(d * n * d) / z == p,
        'centered_cross_occupation_is_strictly_positive': cov == t**2 / (1 + t**2)**2 and cov.is_positive,
        'diagonal_time_cannot_send_an_occupation_to_the_other_face': phase * n == n * phase,
        'both_occupations_are_conserved_by_conjugation': phase * n * y * phase.inv() == n * ht
        and phase * y * n * phase.inv() == ht * n,
        'infinite_strip_weight_limit_collapses_the_even_block': d.subs(t, 0).rank() == 1,
    }
    for name, passed in checks.items():
        print(f'{name}: {"PASS" if passed else "FAIL"}')
    print(f'{sum(checks.values())}/{len(checks)}; determinant={det}; covariance={cov}')
    if not all(checks.values()):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
