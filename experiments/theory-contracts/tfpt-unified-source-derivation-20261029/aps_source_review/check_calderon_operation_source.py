"""R074: original graph replay and exact whole-source identity controls.

Infinite class statements are proved analytically in the companion text.
Numerical windows below only check the tail-controlled whole-sea formulas.
"""
import json
from pathlib import Path
import sys
import mpmath as mp
import sympy as s

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / 'experiments/theory-contracts/seam-collar-multifield-20261016'))
from boundary_amplitude import graph_projection


def run():
    checks = []
    def check(ok, label):
        if not ok:
            raise AssertionError(label)
        checks.append(label)
    def equal(a, b, label):
        d = a-b
        entries = list(d) if isinstance(d, s.MatrixBase) else [d]
        check(all(s.factor(x) == 0 for x in entries), label)

    t, z = s.symbols('t z', positive=True)
    G = graph_projection(s.Matrix([[1/t]]))
    expected = s.Matrix([[t*t,t],[t,1]])/(1+t*t)
    equal(G, expected, 'unchanged old graph builder and original face labels')
    equal(G*G, G, 'whole two-face projector')
    Gz = G.subs(t,z*t)
    C = 1/(1+t*t)
    delta = (Gz-G)**2
    equal(s.trace(delta), 2*(z-1)/(z+1)*(C-C.subs(t,z*t)),
          'exact full shift norm identity at arbitrary positive spectral ratio')
    a,b,c,d = s.symbols('a b c d')
    A = s.Matrix([[a,b],[c,d]])
    P = s.diag(1,0)
    tail = P*A-A*P
    equal(tail, s.Matrix([[0,b],[-c,0]]), 'UV tail forces both off-face blocks zero')
    diagA = s.diag(a,d)
    equal((G*diagA-diagA*G)[0,1], (d-a)*t/(1+t*t),
          'same-momentum cross coherence forces same face weight')
    equal((Gz*diagA-diagA*G)[1,1], d*(C.subs(t,z*t)-C),
          'different momentum strictly monotone diagonal condition')
    equal(s.trace((P*s.Matrix([[0,1],[1,0]])-s.Matrix([[0,1],[1,0]])*P).H
                  *(P*s.Matrix([[0,1],[1,0]])-s.Matrix([[0,1],[1,0]])*P)), 2,
          'bare face swap has nondecaying squared UV defect2')
    for m in [-2,-1,0,1,2]:
        remaining = [(k,l) for k in range(8) for l in range(8) if 8*m+k-l == 0]
        check(remaining == ([(k,k) for k in range(8)] if m == 0 else []),
              f'independent finite residue shift classification m{m}')
    fold = json.loads((REPO / 'experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json').read_text())
    M = s.Matrix(fold['explicit_Weyl_matrix'])
    root = s.Matrix([0,1,-1,0,0,0,0,0])
    q = M.T*root
    equal((q.T*q)[0], 2, 'unchanged original M gives actual q level2')

    # beta=4pi and original quarters only. Whole infinite identities have
    # analytic proofs; these finite sums have exponentially tiny tails.
    mp.mp.dps = 65
    beta = 4*mp.pi
    def occ(p): return 1/(1+mp.exp(beta*p))
    def block(p):
        v = mp.exp(-beta*p/2)
        return mp.matrix([[1,v],[v,v*v]])/(1+v*v)
    rows = []
    for alpha in [mp.mpf(1)/4,mp.mpf(3)/4]:
        for m in [-3,-1,1,3,8]:
            total = mp.mpf(0)
            telescope = mp.mpf(0)
            for r in range(-36,37):
                diff = block(r+alpha+m)-block(r+alpha)
                total += mp.fsum(abs(diff[i,j])**2 for i in range(2) for j in range(2))
                telescope += occ(r+alpha)-occ(r+alpha+m)
            target = 2*m*mp.tanh(beta*m/4)
            check(abs(total-target) < mp.mpf('1e-55'),f'original whole-source shift norm alpha{alpha} m{m}')
            check(abs(telescope-m) < mp.mpf('1e-55'),f'original signed whole-sea telescoping alpha{alpha} m{m}')
            rows.append({'alpha':str(alpha),'m':m,'norm2':str(total),'target':str(target)})
        for m in [-1,1,2]:
            total = mp.mpf(0)
            for k in range(8):
                norm = mp.mpf(0)
                for n in range(-12,13):
                    diff=block(8*(n+m)+k+alpha)-block(8*n+k+alpha)
                    norm+=mp.fsum(abs(diff[i,j])**2 for i in range(2) for j in range(2))
                total += mp.mpf(str(q[k]))**2 * norm
            target = 4*m*mp.tanh(2*beta*m)
            check(abs(total-target) < mp.mpf('1e-55'),f'original weighted folded common-source norm alpha{alpha} m{m}')
    data = {'result':'R074','passed':len(checks),'checks':checks,'q':[str(x) for x in q],
            'original_graph_builder':'seam-collar-multifield-20261016.boundary_amplitude.graph_projection',
            'whole_class_proof':'CALDERON-OPERATION-SOURCE-PROOF.txt',
            'norm_controls':rows,'finite_sums_are_not_all_mode_proofs':True,
            'complete_TFPT_solution':False,'ledger_changed':False}
    Path(__file__).with_name('CALDERON-OPERATION-SOURCE.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({'result':'R074','passed':len(checks),'complete_TFPT_solution':False}))


if __name__ == '__main__':
    run()
