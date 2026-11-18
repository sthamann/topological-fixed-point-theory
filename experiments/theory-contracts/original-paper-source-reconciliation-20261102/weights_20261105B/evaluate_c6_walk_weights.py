"""Evaluate existing original C6 amplitudes; no physical-time/source promotion."""
from pathlib import Path
from fractions import Fraction as F
from math import ceil, isqrt, log
import hashlib
import json
import sys

import mpmath as mp

ROOT = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'verification'))
from v17_hexagonal_resolvent import U6_power

mp.mp.dps = 75
delta = ((mp.mpf(794) - 7 * mp.sqrt(9961)) / 2187) ** (mp.mpf(1) / 6)

# A rational enclosure of the original algebraic root, verified by exact
# integer/rational operations against an enclosure of sqrt(9961).
den = 10**50
delta_lo = F(int(mp.floor(delta * den)), den)
delta_hi = delta_lo + F(1, den)
sd = 10**80
slo = F(isqrt(9961 * sd * sd), sd)
shi = slo + F(1, sd)
zlo = (794 - 7 * shi) / 2187
zhi = (794 - 7 * slo) / 2187
assert delta_lo**6 < zlo <= zhi < delta_hi**6

U = mp.matrix(U6_power(1).tolist())
I = mp.eye(6)
A = U + U.T
assert U**6 == I
assert all(sum(A[i, j] for j in range(6)) == 2 for i in range(6))


def m(q):
    return mp.mpf(q.numerator) / q.denominator


def exact_tail_below(c, q, n, target):
    # All positive, and q is a certified rational upper bound.
    return c * q**(n + 1) < target * (1 - q)


def evaluate(y_exact):
    y = m(y_exact)
    a = y * y + delta * delta
    b = y * delta
    q = 2 * b / a
    assert 0 < q < 1
    a_min = y_exact**2 + delta_lo**2
    # q(delta) increases for delta<y and decreases for delta>y.
    assert delta_hi < y_exact or delta_lo > y_exact
    d_for_max = delta_hi if delta_hi < y_exact else delta_lo
    qmax = 2 * y_exact * d_for_max / (y_exact**2 + d_for_max**2)
    qden = 10**9
    q_upper = F(ceil(qmax * qden), qden)
    assert qmax < q_upper < 1
    pref_green = ceil(1 / a_min)
    sqrt_den = 10**50
    sqrt_amin_lower = F(isqrt(a_min.numerator * sqrt_den**2 // a_min.denominator), sqrt_den)
    assert sqrt_amin_lower**2 <= a_min
    pref_sqrt = ceil(1 / sqrt_amin_lower)
    pref = max(pref_green, pref_sqrt)
    target = F(1, 10**35)
    nmax = max(0, ceil(log(float(target * (1 - q_upper) / pref)) / log(float(q_upper))) - 1)
    while not exact_tail_below(pref, q_upper, nmax, target):
        nmax += 1

    # Circulant powers retain all oriented + and - walks on the original C6.
    # v_n is the first row of ((b/a) A)^n; no stochastic normalization occurs.
    v = [mp.mpf(1), *[mp.mpf(0)] * 5]
    sum_green = [mp.mpf(0)] * 6
    sum_sqrt = [mp.mpf(0)] * 6
    c = mp.mpf(1)
    for n in range(nmax + 1):
        for j in range(6):
            sum_green[j] += v[j] / a
            sum_sqrt[j] += c * v[j] / mp.sqrt(a)
        v = [(b / a) * (v[(j - 1) % 6] + v[(j + 1) % 6]) for j in range(6)]
        c *= mp.mpf(2 * n + 1) / (2 * n + 2)

    G = a * I - b * A
    direct_green = G**-1
    eig, Q = mp.eigsy(G)
    assert min(eig) > 0
    direct_sqrt = Q * mp.diag([1 / mp.sqrt(x) for x in eig]) * Q.T
    walk_green = mp.matrix(6, 6)
    walk_sqrt = mp.matrix(6, 6)
    for i in range(6):
        for j in range(6):
            walk_green[i, j] = sum_green[(j - i) % 6]
            walk_sqrt[i, j] = sum_sqrt[(j - i) % 6]
    # Numerical errors are evaluated separately from the exact rational tail.
    err_green = max(abs(walk_green[i,j] - direct_green[i,j]) for i in range(6) for j in range(6))
    err_sqrt = max(abs(walk_sqrt[i,j] - direct_sqrt[i,j]) for i in range(6) for j in range(6))
    identity_green = mp.norm(G * direct_green - I, p=mp.inf)
    identity_sqrt = mp.norm(direct_sqrt * G * direct_sqrt - I, p=mp.inf)
    assert err_green < m(target) and err_sqrt < m(target)
    pf_green = 1 / (y - delta)**2
    pf_sqrt = 1 / abs(y - delta)
    assert abs(sum(direct_green[0,j] for j in range(6)) - pf_green) < mp.mpf('1e-65')
    assert abs(sum(direct_sqrt[0,j] for j in range(6)) - pf_sqrt) < mp.mpf('1e-65')
    sample_weights = []
    for n in [0,1,2,3,6]:
        cn = mp.binomial(2*n,n) / 4**n
        sample_weights.append(dict(length=n, green_per_oriented_walk=mp.nstr((b/a)**n/a,45),
                                   inverse_sqrt_per_oriented_walk=mp.nstr(cn*(b/a)**n/mp.sqrt(a),45)))
    return dict(y=str(y_exact), epsilon=0, a=mp.nstr(a,55),b=mp.nstr(b,55),q=mp.nstr(q,55),
                certified_q_upper=str(q_upper), certified_prefactor_green_upper=pref_green,
                certified_prefactor_inverse_sqrt_upper=pref_sqrt, N=nmax,
                exact_rational_tail_verified_below='1e-35 for both operator norms',
                tail_formula='C*q_upper^(N+1)/(1-q_upper)',
                green_max_entry_error=mp.nstr(err_green,10),inverse_sqrt_max_entry_error=mp.nstr(err_sqrt,10),
                green_identity_residual=mp.nstr(identity_green,10),inverse_sqrt_identity_residual=mp.nstr(identity_sqrt,10),
                green_absolute_PF=mp.nstr(pf_green,55),inverse_sqrt_absolute_PF=mp.nstr(pf_sqrt,55),
                green_first_row=[mp.nstr(direct_green[0,j],45) for j in range(6)],
                inverse_sqrt_first_row=[mp.nstr(direct_sqrt[0,j],45) for j in range(6)],
                sample_weights=sample_weights)


sources = ['verification/v17_hexagonal_resolvent.py',
           'documents/archive--tfpt-45--source_extracts--03_em_flavor_source/sections/003-subsection-discrete-phase-operator-and-transport-grammar.tex',
           'documents/archive--tfpt-45--source_extracts--03_em_flavor_source/sections/005-source.tex']
data = dict(scope='Existing full local C6 Green and structural inverse-square-root walk amplitudes. No physical time, Gibbs law, primitive mixed field source or new parameter selected.',
            delta=mp.nstr(delta,65), delta_interval=[str(delta_lo),str(delta_hi)],
            exact_delta_interval_verified=True, matrix_producer='Original v17_hexagonal_resolvent.U6_power(1)',
            green_weight_rule='a^-1*(b/a)^n per oriented walk of length n',
            inverse_sqrt_weight_rule='binomial(2n,n)/(4^n*sqrt(a))*(b/a)^n per oriented walk of length n',
            source_sha256={s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources},
            rows=[evaluate(y) for y in [F(1),F(1,3),F(2,3)]])
(OUT/'c6-walk-weight-evaluation.json').write_text(json.dumps(data,indent=2)+'\n')
for row in data['rows']:
    print(json.dumps({k:row[k] for k in ['y','N','q','green_absolute_PF','inverse_sqrt_absolute_PF','green_max_entry_error','inverse_sqrt_max_entry_error','exact_rational_tail_verified_below']}))
print('PASS: all three original cusp matrices; Green and inverse square root kept distinct; exact rational tails below1e-35; absolute factors retained.')
