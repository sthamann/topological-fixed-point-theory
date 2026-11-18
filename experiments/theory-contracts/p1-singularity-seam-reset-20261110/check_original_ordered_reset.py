"""Exact ordered reset contrast on the two original Ramond preparations."""

from contextlib import redirect_stdout
from pathlib import Path
import hashlib
import io
import json
import runpy
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[3]
CONTRACTS = ROOT / 'experiments/theory-contracts'


def folded_forms():
    path = CONTRACTS / 'chat-triangulation-20261108/check_tfpt_family_middle_ordered_cut_20261110.py'
    saved_args = sys.argv
    try:
        sys.argv = [str(path), str(ROOT)]
        with redirect_stdout(io.StringIO()):
            source = runpy.run_path(str(path))
    finally:
        sys.argv = saved_args
    q, r, b, axes = (source[k] for k in ('q', 'r', 'd1', 'V'))
    old, field, gram, pairing = (source[k] for k in ('old', 'vertex', 'gram', 'pairing'))
    for alpha in (sp.Rational(1, 4), sp.Rational(3, 4)):
        drift = sp.Matrix([j + alpha - 4 for j in range(8)])
        shift = [source['Q'].from_sympy(z) for z in old['cl'](axes.T * drift / 8)]
        packets = [old['vertex_zero'](old['series_response'](p, shift, 4))
                   for p in source['polys'][:2]]
        late = [field(field(v, b, q), q, r, 1) for v in packets]
        early = [field(field(v, q, q, 1), b, 2*q) for v in packets]
        assert q.dot(b) == -2 and drift.dot(r) == 2*alpha - 2 and drift.dot(q) == 0
        assert source['rank'](gram(late + early)) == 4
        yield alpha, [gram(packets), gram(late), gram(early), pairing(late, early)]


def compression_system(forms):
    base, late, early, cross = forms
    x, y, h = sp.symbols('x y h', real=True)
    marked = late - (x+sp.I*y)*cross - (x-sp.I*y)*cross.H + h*early
    residual = (base[0, 0]*marked - marked[0, 0]*base).applyfunc(sp.expand)
    equations = [residual[1, 1], sp.re(residual[0, 1]), sp.im(residual[0, 1])]
    matrix, target = sp.linear_eq_to_matrix(equations, [x, y, h])
    det = sp.factor(matrix.det())
    numerators = [sp.factor(sp.Matrix.hstack(matrix[:, :i], target, matrix[:, i+1:]).det())
                  for i in range(3)]
    # Positive covariance requires h>=x²+y² when the first diagonal is 1.
    delta_numerator = sp.factor(numerators[2]*det - numerators[0]**2 - numerators[1]**2)
    return det, delta_numerator


def probe_original_adapter():
    path = CONTRACTS / 'original-paper-source-reconciliation-20261102/contribution_20261106Q'
    domain_path = path / 'domain/native-diagnostic-v1.json'
    domain = json.loads(domain_path.read_text())
    domain['_sha256'] = hashlib.sha256(domain_path.read_bytes()).hexdigest()
    adapter = runpy.run_path(str(path / 'runtime/original_source_adapter.py'))
    plan = adapter['integration_plan'](domain)
    attempt = adapter['evaluate_contribution'](
        {'schema': 'tfpt.original-source-unit.v1', 'original_configuration_id': None},
        {'native_inputs_sha256': domain['_sha256'], 'domain_id': domain['domain_id']}, domain)
    assert not plan['ready'] and plan['items'] == []
    assert attempt['missing_operations'] == ['original_blocked_fiber_packet']
    assert attempt['original_measure_mass'] is None
    assert all(v is None for v in attempt['contributions'].values())
    return {'atoms': 0, 'missing': attempt['missing_operations'],
            'mass': attempt['original_measure_mass'], 'contributions': attempt['contributions']}


def main():
    quarters = list(folded_forms())
    t = sp.Symbol('t', real=True)
    for alpha, forms in quarters:
        base, late, early, cross = forms
        assert base.is_positive_definite and sp.re(base[0, 1]) == 0
        det, delta = compression_system(forms)
        assert det != 0 and delta < 0
        # A pure early word cannot be scalar either, so a zero late weight is covered.
        assert sp.re(early[0, 1]) != 0
        real_contrast = sp.Poly(sp.re((late-t*(cross+cross.H)+t*t*early)[0, 1]), t)
        assert real_contrast.nth(1) < 0 and real_contrast.nth(2) < 0
        minimum = real_contrast.eval(1)
        assert minimum > 0
        print(f'alpha={alpha} rank=4; min real contrast={minimum}; '
              f'h-|z|² at formal cancellation={sp.N(delta/det**2, 12)}')

    w = sp.Symbol('w', nonnegative=True)
    joint = [a + w*b for a, b in zip(quarters[0][1], quarters[1][1], strict=True)]
    det, delta = compression_system(joint)
    det_poly = sp.Poly(det, w)
    delta_poly = sp.Poly(delta, w)
    assert all(c > 0 for c in det_poly.all_coeffs())
    assert all(c < 0 for c in delta_poly.all_coeffs())
    assert all(sp.re(forms[2][0, 1]) < 0 for _, forms in quarters)
    print(f'joint quarter weights with common covariance: exact degree {det_poly.degree()} '
          f'and {delta_poly.degree()} coefficient certificates')
    print('installed original adapter:', json.dumps(probe_original_adapter()))


if __name__ == '__main__':
    main()
