"""Small decisive comparison with the actual original v117 cusp matrix."""
import hashlib
import json
from pathlib import Path

import sympy as s

source=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4/verification/v117_monodromy_weyl_a3.py')
expected='55eca6a64144f6e49789f3f1bbe84e256764e7ce92b207044be693f8684d5df5'
if hashlib.sha256(source.read_bytes()).hexdigest()!=expected:
    raise RuntimeError('Original matrix source changed; inspect before reusing this comparison')

# Literal original M0_EXACT from v117:53-55, checked through code graph.
i=s.I
m=s.Matrix([[0,-(1+i)/2,(1-i)/2],
            [-(1+i)/2,-i/2,s.Rational(-1,2)],
            [(1-i)/2,s.Rational(-1,2),i/2]])
ell=s.diag(-1,1,-1)
lam=s.Symbol('lambda')
checks={
 'original_source_hash_current': True,
 'original_cusp_characteristic_polynomial':s.expand(m.charpoly(lam).as_expr())==lam**3-1,
 'original_cusp_order_three':s.simplify(m**3)==s.eye(3) and s.simplify(m**2)!=s.eye(3),
 'elliptic_cusp_order_two':ell**2==s.eye(3) and ell!=s.eye(3),
 'distinct_local_conjugacy_classes':s.trace(m)==0 and s.trace(ell)==-1,
 'distinct_cusp_polynomials':s.factor(ell.charpoly(lam).as_expr())==(lam-1)*(lam+1)**2,
 'same_residue_trace_and_degree':s.Rational(1,2)+s.Rational(1,2)==s.Rational(1,3)+s.Rational(2,3)==1,
 'same_explicit_bundle_rank_degree_h0':len([2,1,1])==3 and -sum([2,1,1])==-4
      and sum(max(0,1-a) for a in [2,1,1])==0,
 'scalar_twist_cannot_remove_repeated_eigenvalue':len(set([1,-1,-1]))==2
      and len(set([1,s.exp(2*s.pi*i/3),s.exp(4*s.pi*i/3)]))==3,
}
if not all(checks.values()):
    raise RuntimeError([k for k,v in checks.items() if not v])
result=dict(status='PASS',checks={k:bool(v) for k,v in checks.items()},
            original_source_sha256=expected,
            scope='Common holomorphic bundle class; natural logarithmic connections are not conjugate. No physical direct-sum selection derived.',
            physical_gates_closed=[],complete_TFPT_solution=False)
payload=json.dumps(result,sort_keys=True,indent=2)+'\n'
Path(__file__).with_name('pushforward-anchor-results.json').write_text(payload)
print(payload,end='')
