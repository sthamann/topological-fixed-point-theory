"""Exact decision test on the existing D5+A3 discriminant modules.

This tests chiral cycle-exchange closure, not a physical P1 sewing law.
"""
import itertools
import json
import sympy as s


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


labels = list(itertools.product(range(4), repeat=2))
idx = {a: i for i, a in enumerate(labels)}
S = s.Matrix(16, 16, lambda i, j: s.I ** (-(5*labels[i][0]*labels[j][0]
                                           + 3*labels[i][1]*labels[j][1]) % 4) / 4)
T = s.diag(*[s.exp(s.pi*s.I*s.Rational(5*a*a+3*b*b, 4)).expand(complex=True)
             for a, b in labels])
C = s.zeros(16)
for a in labels:
    C[idx[a], idx[(-a[0] % 4, -a[1] % 4)]] = 1
require(S.conjugate().T*S == s.eye(16), 'unitary finite Fourier action')
require(S*S == C, 'S squared is charge conjugation')
require(s.simplify(T.conjugate().T*T) == s.eye(16), 'unitary twists')

supports = {
    'U_D5_A3': [(0, 0)],
    'D8': [(0, 0), (2, 2)],
    'E8_plus': [(a, a) for a in range(4)],
    'E8_minus': [(a, (-a) % 4) for a in range(4)],
}
out = {}
for name, H in supports.items():
    v = s.Matrix([int(a in H) for a in labels])
    perp = [a for a in labels if all((5*a[0]*h[0]+3*a[1]*h[1]) % 4 == 0 for h in H)]
    expected = s.Rational(len(H), 4)*s.Matrix([int(a in perp) for a in labels])
    require(S*v == expected, 'Fourier support identity '+name)
    require(T*v == v, 'local bosonic support '+name)
    fixed = S*v == v
    require(fixed == (len(H) == 4), 'closure classification '+name)
    out[name] = {
        'support': H,
        'support_size': len(H),
        'orthogonal_support_size': len(perp),
        'vacuum_coefficient_after_S': str((S*v)[0]),
        'single_chiral_vacuum_fixed': fixed,
        'contains_marked_h_1_1': (1, 1) in H,
    }

# A complete scalar paired amplitude can be modular invariant without a
# one-dimensional chiral block: M=I pairs every module with its conjugate.
M = s.eye(16)
require(S.conjugate().T*M*S == M, 'paired S covariance')
require(s.simplify(T.conjugate().T*M*T) == M, 'paired T covariance')
require(M.rank() == 16 and M[0, 0] == 1, 'paired amplitude is not chiral rank one')

print(json.dumps({
    'verdict': 'CHIRAL_VACUUM_CYCLE_EXCHANGE_SELECTS_TWO_E8_EXTENSIONS',
    'formula': 'S v_H = (|H|/4) v_H_perp',
    'extensions': out,
    'paired_scalar_control': {'rank': 16, 'vacuum_multiplicity': 1,
                             'S_and_T_invariant': True},
    'scope': 'existing unitary connected same-conformal extension class; no P1 origin of chiral cycle exchange, no physical countermodel',
}, indent=2, sort_keys=True))
