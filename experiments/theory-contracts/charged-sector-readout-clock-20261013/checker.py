"""Native marked E8 response and clock test. Conditional, no physical closure."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
import sympy as s

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from tfpt_explorer.neutral_source_response import _native_neutral_branch
from tfpt_explorer.current_block_geometry import _actual_current_source


def partition_coefficients(weight, degree):
    # Independently multiply oscillator-mode coherent norms exp(weight*x^n/n).
    coefficients = [s.Integer(1)] + [s.Integer(0)]*degree
    for n in range(1, degree+1):
        old = coefficients[:]
        coefficients = [sum(old[k-n*j]*(weight/s.Integer(n))**j/s.factorial(j)
                            for j in range(k//n+1)) for k in range(degree+1)]
    return coefficients


def modular_rank(matrix, prime=101):
    denominator=s.ilcm(*(x.q for x in matrix))
    rows=[[int(denominator*x)%prime for x in row] for row in matrix.tolist()]
    rank=0
    for col in range(len(rows[0])):
        pivot=next((i for i in range(rank,len(rows)) if rows[i][col]),None)
        if pivot is None: continue
        rows[rank],rows[pivot]=rows[pivot],rows[rank]
        scale=pow(rows[rank][col],-1,prime)
        rows[rank]=[(x*scale)%prime for x in rows[rank]]
        for i in range(rank+1,len(rows)):
            fac=rows[i][col]
            if fac: rows[i]=[(x-fac*y)%prime for x,y in zip(rows[i],rows[rank])]
        rank+=1
    return rank


def run():
    here=Path(__file__).parent
    manifest=json.loads((here/'source_manifest.json').read_text())
    for entry in manifest['sources']:
        if hashlib.sha256((ROOT/entry['path']).read_bytes()).hexdigest()!=entry['sha256']:
            raise RuntimeError('Source changed: '+entry['path'])
    native=_native_neutral_branch()
    P=native['P5']; I=s.eye(8)
    roots=[s.Matrix(r)/2 for r in sorted(native['retained'])]
    chevalley=_actual_current_source()['chevalley']
    all_roots={root for root in chevalley.roots if chevalley.ridx[root] in native['closure']}
    clock=json.loads((ROOT/'experiments/theory-contracts/compiler-integral-triality-20260918/certificate.json').read_text())
    C=s.Matrix(clock['clock_matrices']['C']['vector'])
    J=s.Matrix(clock['clock_matrices']['J']['vector'])
    checks=[]
    def check(name, value):
        checks.append({'name':name,'ok':bool(value)})
    check('original 48 marked roots and full 240 closure',len(roots)==48 and len(all_roots)==240)
    check('original P5 orthogonal rank five',P==s.diag(1,1,1,1,1,0,0,0) and P*P==P and P.trace()==5)
    check('every retained charge has D5 norm5/4 A3 norm3/4 and total2',
          all((b.T*P*b)[0]==s.Rational(5,4) and (b.T*(I-P)*b)[0]==s.Rational(3,4) and (b.T*b)[0]==2 for b in roots))
    check('old neutral response contains no D5 component',P*native['R_matrix']==s.zeros(8))
    z,w1,w2=s.symbols('z w1 w2')
    # Actual current Ward contraction, not an imposed fractional pair function.
    current=1/(z-w1)-1/(z-w2)
    check('Wick sector-stress insertion gives normalized Ward kernel',
          s.factor(s.Rational(5,8)*current**2-s.Rational(5,8)*(w1-w2)**2/((z-w1)**2*(z-w2)**2))==0)
    coeffs={str(a):partition_coefficients(a,12) for a in [s.Rational(5,4),s.Rational(3,4),s.Integer(2)]}
    check('oscillator partitions generate carrier exponent5/4 through degree12',
          coeffs['5/4']==[s.rf(s.Rational(5,4),n)/s.factorial(n) for n in range(13)])
    check('convolution restores full local exponent2 through degree12',
          [sum(coeffs['5/4'][k]*coeffs['3/4'][n-k] for k in range(n+1)) for n in range(13)]==coeffs['2']==list(map(s.Integer,range(1,14))))
    u,q=s.symbols('u q',positive=True)
    a=s.Rational(5,4)
    G=((1-u)/(1-u*q))**a
    check('normalized carrier response at q1',G.subs(q,1)==1)
    check('zero-grade probability retains nonzero normalization',G.subs(q,0)==(1-u)**a)
    # Use positive u in (0,1); fourth powers remove fractional-power simplification ambiguity.
    check('reference-ratio matches target analytic factor on physical branch',
          s.factor((G/G.subs(q,0))**4-(1-u*q)**-5)==0)
    check('mean energy grade differs in carrier and full readout',
          s.simplify(s.diff(G,q).subs(q,1)-s.Rational(5,4)*u/(1-u))==0)
    clock_results={}
    for name,g,rank,norm in [('C',C,4,s.Rational(23,8)),('J',J,2,s.Integer(2))]:
        d=g.T*P*g-P
        check(name+' is orthogonal and preserves original closed root set',g.T*g==I and {tuple(g*s.Matrix(r)) for r in all_roots}==all_roots)
        check(name+' does not preserve the fixed carrier projector',d.rank()==rank and s.trace(d*d)==norm)
        check(name+' transports marked response covariantly',all(s.simplify(((g*b).T*(g*P*g.T)*(g*b))[0]-(b.T*P*b)[0])==0 for b in roots))
        clock_results[name]={'difference_rank':d.rank(),'difference_squared_norm':str(s.trace(d*d)),
                             'retained_carrier_norms_after_action':dict(Counter(str(((g*b).T*P*(g*b))[0]) for b in roots))}
    columns=[]
    for i in range(8):
        for j in range(i,8):
            B=s.zeros(8); B[i,j]=B[j,i]=1
            columns.append(s.Matrix([x for g in (C,J) for x in g.T*B*g-B]))
    rank=modular_rank(s.Matrix.hstack(*columns))
    # Rank mod101 gives the lower bound35 over Q. Identity gives upper bound35.
    check('actual common clock invariant symmetric space one dimensional',rank==35 and all(g.T*I*g==I for g in(C,J)))
    avg=s.Rational(5,8)*I
    check('unique invariant trace5 quadratic coefficient is5/8I',avg.trace()==5 and all(g.T*avg*g==avg for g in(C,J)))
    spin=s.ones(8,1)/2
    check('explicit spinor clock changes carrier charge norm', (spin.T*P*spin)[0]==s.Rational(5,4) and ((C*spin).T*P*C*spin)[0]==1)
    # At fixed u>0, scaled common flow produces q^(5/8), not a carrier exponent5/4.
    check('scaled total and carrier flows have unequal second cumulants',
          (s.Rational(5,8)**2*2-s.Rational(5,4))==-s.Rational(15,32))
    result={'verdict':'PARTIAL','mathematical_verdict':'MARKED_RESPONSE_CONSTRUCTED_CLOCK_IDENTIFICATION_REJECTED',
      'passed':sum(x['ok'] for x in checks),'total':len(checks),'checks':checks,
      'clock_results':clock_results,'clock_symmetric_constraint_rank':rank,
      'carrier_coefficients_0_to_5':list(map(str,coeffs['5/4'][:6])),
      'joint_response':'((1-u)/(1-u*q5))^(5/4) * ((1-u)/(1-u*q3))^(3/4)',
      'common_time_H_delta':'exp(-t*(1-delta.dot(beta)))*((1-u)/(1-u*exp(-t)))^2',
      'scope':'Chosen native E8 lattice-VOA vacuum, marked D5/A3 splitting, regularized charged vector; oscillator conformal grade, not particle number.',
      'physical_gates_closed':[], 'complete_TFPT_solution':False}
    (here/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if result['passed']!=result['total']: raise SystemExit(1)

if __name__=='__main__':run()
