"""Exact native-source response identifiability; no physical closure claim."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

import sympy as s

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from tfpt_explorer.neutral_source_response import _native_neutral_branch


def pairings(items):
    if not items:
        yield []
        return
    for index in range(1, len(items)):
        for rest in pairings(items[1:index] + items[index + 1:]):
            yield [(items[0], items[index])] + rest


def wick_cycles(n):
    counts = Counter()
    for matching in pairings(list(range(2 * n))):
        edges = [tuple(sorted((i // 2, j // 2))) for i, j in matching]
        if any(i == j for i, j in edges):
            continue
        seen = {0}
        for _ in range(n):
            seen |= {j for i, j in edges if i in seen}
            seen |= {i for i, j in edges if j in seen}
        if len(seen) == n:
            counts[tuple(sorted(edges))] += 1
    return counts


def run():
    manifest=json.loads(Path(__file__).with_name('source_manifest.json').read_text())
    for source in manifest['sources']:
        actual_hash=hashlib.sha256((ROOT/source['path']).read_bytes()).hexdigest()
        if actual_hash != source['sha256']:
            raise RuntimeError('Pinned source changed: '+source['path'])
    native = _native_neutral_branch()
    p5, p2, p1 = [native[k] for k in ('P5', 'P2', 'P1')]
    actual = native['R_matrix']
    a, b = s.symbols('a b')
    c = -5*a - 2*b
    candidate = a*p5 + b*p2 + c*p1
    old_witness = s.Rational(16, 7)*p5-s.Rational(12, 7)*p2-8*p1
    cubic_witness = 2*p5-6*p2+2*p1
    checks = []

    def check(name, condition):
        checks.append({'name': name, 'ok': bool(condition)})

    check('native projectors orthogonal ranks5,2,1',
          [p.rank() for p in (p5,p2,p1)] == [5,2,1]
          and p5+p2+p1 == s.eye(8)
          and all(x*x == x for x in (p5,p2,p1))
          and all(x*y == s.zeros(8) for x,y in ((p5,p2),(p5,p1),(p2,p1))))
    check('native R equals 4P2-8P1', actual == 4*p2-8*p1)
    check('native primary and norm48', actual.trace() == 0
          and s.trace(actual**2)/2 == native['norms']['R'] == 48)
    check('native cubic-384 quartic4608', s.trace(actual**3) == -384
          and s.trace(actual**4) == 4608)
    check('prior witness rejected by cubic', s.trace(old_witness**2) == 96
          and s.trace(old_witness**3) == -s.Rational(22656,49))
    check('new witness retains primary norm and cubic', cubic_witness.trace() == 0
          and s.trace(cubic_witness**2) == 96 and s.trace(cubic_witness**3) == -384
          and s.trace(p5*cubic_witness)/5 == 2)
    check('quartic separates cubic witness', s.trace(cubic_witness**4) == 2688)
    equations = [s.expand(s.trace(candidate**k)-s.trace(actual**k)) for k in (2,3,4)]
    gb = s.groebner(equations, a, b, domain=s.QQ)
    check('exact response ideal is a,b-4', [x.as_expr() for x in gb.polys] == [a,b-4])
    resultants = [s.factor(s.resultant(equations[0], f, b)) for f in equations[1:]]
    check('independent elimination forces a zero',
          s.Poly(s.gcd(*resultants), a).monic().as_expr() == a)
    check('source satisfies every response equation',
          all(f.subs({a:0,b:4}) == 0 for f in equations))
    cycles = {n:wick_cycles(n) for n in (2,3,4)}
    check('two point Wick prefactor one half', list(cycles[2].values()) == [2])
    check('three point Wick prefactor one', list(cycles[3].values()) == [8])
    check('four point three cycles each prefactor one', len(cycles[4]) == 3
          and set(cycles[4].values()) == {16})
    z = [s.Integer(x) for x in (0,1,3,7)]
    kernel = sum(s.prod((z[i]-z[j])**-2 for i,j in edges) for edges in cycles[4])
    check('separated point four point discriminator nonzero', kernel > 0
          and (4608-2688)*kernel != 0)
    # Isospectral source fields retain every self-response but change carrier alignment.
    # Reflection swaps a native D5 unit vector with a native P2 unit vector.
    e = s.eye(8)[:,0]
    f = next(p2[:,i] for i in range(8) if p2[:,i] != s.zeros(8,1))
    f = f/s.sqrt((f.T*f)[0])
    u = e-f
    reflection = s.eye(8)-2*u*u.T/(u.T*u)[0]
    rotated = s.simplify(reflection*actual*reflection.T)
    check('isospectral carrier witness', reflection.T*reflection == s.eye(8)
          and s.trace(p5*rotated) == 4
          and all(s.trace(rotated**k) == s.trace(actual**k) for k in range(1,9)))
    check('carrier mixed response detects isospectral witness',
          s.trace(p5*actual)/2 == 0 and s.trace(p5*rotated)/2 == 2)
    m5, m2 = s.symbols('m5 m2')
    check('two native mixed probes reconstruct block coefficients',
          s.trace(p5*candidate)/2 == 5*a/2
          and s.trace(p2*candidate)/2 == b
          and s.simplify(c.subs({a:2*m5/5,b:m2}) + 2*(m5+m2)) == 0)
    result = {
        'verdict':'PARTIAL', 'mathematical_verdict':'SCOPED_RESPONSE_IDENTIFIABILITY_PROVED',
        'passed':sum(x['ok'] for x in checks), 'total':len(checks), 'checks':checks,
        'source_sha256':hashlib.sha256((ROOT/'tfpt_explorer/neutral_source_response.py').read_bytes()).hexdigest(),
        'native_moments':{'m1':0,'m2':96,'m3':-384,'m4':4608},
        'cubic_counterexample':{'coefficients':[2,-6,2],'moments':[0,96,-384,2688]},
        'response_ideal':[str(x.as_expr()) for x in gb.polys],
        'elimination_resultants':list(map(str,resultants)),
        'wick_connected_counts':{str(n):sorted(x.values()) for n,x in cycles.items()},
        'four_point_kernel_at_0_1_3_7':str(kernel),
        'block_reconstruction':{'a':'2*m5/5','b':'m2','c':'-2*(m5+m2)'},
        'scope':'Fixed native 5+2+1 projectors, real block-scalar quadratic field and original vacuum Wick algebra. Physical access to probes is not derived.',
        'not_claimed':['5/4 physical susceptibility','source selection by P1/P2','full matrix reconstruction from self responses','T1-T8 closure']}
    Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if result['passed'] != result['total']:
        raise SystemExit(1)


if __name__ == '__main__':
    run()
