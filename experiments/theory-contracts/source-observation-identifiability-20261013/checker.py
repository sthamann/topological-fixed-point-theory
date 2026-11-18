import json
from pathlib import Path
import sympy as s
n,A,lam=s.symbols('n A lam',real=True)
Qrep=4*s.pi*(1-n)*A
shift=s.simplify(-s.diff(lam*Qrep,n).subs(n,1))
if shift!=4*s.pi*lam*A: raise RuntimeError('Replica identity')
g,J=s.symbols('g J')
W=s.Function('W'); Q=s.Function('Q')
normalized=s.simplify((W(g,J)+lam*Q(g))-(W(g,0)+lam*Q(g))-(W(g,J)-W(g,0)))
if normalized!=0: raise RuntimeError('Normalization identity')
if s.diff(Q(g),J,g)!=0: raise RuntimeError('Mixed derivatives kill metric-only term')
p=s.Rational(5,4);f=s.Rational(3,4)
if p+f!=2: raise RuntimeError('Glue exponent')
r={'normalized_shift':str(normalized),'replica_entropy_shift':str(shift),'carrier_plus_family_exponent':str(p+f),'scope':'Symbolic identities only; source and continuum premises and expectation theorem are not computationally proved here.'}
Path(__file__).with_name('results.json').write_text(json.dumps(r,indent=2))
print(json.dumps(r,indent=2))
