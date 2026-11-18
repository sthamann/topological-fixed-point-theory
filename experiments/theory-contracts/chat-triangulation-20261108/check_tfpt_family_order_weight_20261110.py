"""Evaluate the two orders with the actual old Hamiltonian and q reference."""
from contextlib import redirect_stdout
import io
import json
import runpy
from pathlib import Path
import sys

repo = Path(sys.argv[1]).resolve()
path = repo / 'experiments/theory-contracts/chat-triangulation-20261108/check_tfpt_family_middle_ordered_cut_20261110.py'
with redirect_stdout(io.StringIO()):
    s = runpy.run_path(str(path))

sp, q, r, d1, Q = (s[k] for k in ('sp','q','r','d1','Q'))
old, vertex, middle, gram, pairing, combine = (s[k] for k in ('old','vertex','middle','gram','pairing','combine'))
t, delay = s['t'], s['delta_t']
verdict = {}
for alpha in map(sp.Rational, ('1/4','3/4')):
    drift = sp.Matrix([j+alpha-4 for j in range(8)])
    assert drift.dot(q) == 0 and 2*q == q+q
    inputs = [old['series_response'](f,[Q.from_sympy(z) for z in old['cl'](s['V'].T*drift/8)],4)
              for f in s['polys']]
    states = [old['vertex_zero'](f) for f in inputs]
    g_in = gram(states)
    late_states = [vertex(vertex(v,d1,q),q,r) for v in states]
    early_states = [vertex(vertex(v,q,q),d1,2*q) for v in states]
    w = [middle(v) for v in states]
    assert [bool(v) for v in early_states] == [True,False,False,False,False]
    w_diff = [combine([*a.items(),*((m,-x) for m,x in b.items())])
              for a,b in zip(late_states,early_states,strict=True)]
    sign = 1 if w_diff == w else -1
    assert w_diff == [{m:sign*x for m,x in v.items()} for v in w]
    after = sign*pairing(w,late_states)
    before = sign*pairing(w,early_states)
    assert after-before == gram(w)
    assert s['rank'](after)==5 and s['rank'](before)==1
    e = 32+drift.dot(r)
    late_factor, early_factor = sp.exp(-e*(t+delay)), sp.exp(-e*t-32*delay)
    d_sub = (after-before)[1:,1:]
    g_sub = g_in[1:,1:]
    assert all(sp.simplify(g_sub[:n,:n].det()) > 0 for n in range(1,5))
    # Exact, shared bound in original q-packet norm. No fit or inserted mass.
    bound = sp.Rational(1,2)
    certificate = d_sub-bound*g_sub
    positive_minors = [sp.factor(certificate[:n,:n].det()) for n in range(1,5)]
    assert all(v>0 for v in positive_minors)
    spin = lambda p: sp.exp(2*sp.pi*sp.I*p[0])
    assert tuple(map(sp.simplify,map(spin,(q,r,2*q,q+r)))) == (-1,-1,1,1)
    assert q + d1 + q == q + r
    fold = s['M']
    assert all(all(component == 0 for component in (fold * p)[:3])
               for p in (q, r, 2*q, q+r))
    weak = s['weak']
    hypercharge = s['hypercharge']
    assert (fold*q).dot(weak) == (fold*q).dot(hypercharge) == 0
    assert (fold*d1).dot(weak) == (fold*r).dot(weak) == -1
    assert (fold*d1).dot(hypercharge) == (fold*r).dot(hypercharge) == sp.Rational(1,2)
    def as_complex(z):
        return [float(sp.re(z).evalf(20)),float(sp.im(z).evalf(20))]
    def diagonal(a,factor):
        return [as_complex(a[i,i]*factor/g_in[i,i]) for i in range(5)]
    verdict[str(alpha)] = {
        'common_reference': [str(g_in[i,i]) for i in range(5)],
        'initial_energy': 32,
        'charged_intermediate_energy': str(e),
        'spins_input_early_intermediate_late_intermediate_output': [-1,1,-1,1],
        'late_exact_diagonal_coefficients_in_reference': [str(sp.expand(after[i,i]/g_in[i,i])) for i in range(5)],
        'early_exact_diagonal_coefficients_in_reference': [str(sp.expand(before[i,i]/g_in[i,i])) for i in range(5)],
        'weighted_late_cross_return': diagonal(after,late_factor),
        'weighted_early_cross_return': diagonal(before,early_factor),
        'weighted_late_reflected_norm': diagonal(gram(late_states),late_factor**2),
        'weighted_early_reflected_norm': diagonal(gram(early_states),early_factor**2),
        'exact_mixed_late_rank': 5,
        'exact_mixed_early_rank': 1,
        'exact_early_null_packets': [1,2,3,4],
        'late_reduced_gram_lower_bound_in_q_norm_before_heat_weight': str(bound),
        'sylvester_certificate': list(map(str,positive_minors)),
    }

print(json.dumps({'source':'old exact lattice current/vertex and H_alpha from existing checker',
                  'window_times':{'return':str(t),'between_insertions':str(delay)},
                  'quarters':verdict,
                  'physical_returns':None},indent=2))
