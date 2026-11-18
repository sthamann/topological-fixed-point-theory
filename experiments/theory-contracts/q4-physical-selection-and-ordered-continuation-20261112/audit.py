from pathlib import Path
import argparse, hashlib, json
import sympy as s
parser=argparse.ArgumentParser()
parser.add_argument('--repo', required=True, type=Path)
parser.add_argument('--output', type=Path, default=Path(__file__).with_name('results.json'))
args=parser.parse_args()
R=args.repo
p=R/'_newest2/TFPT_Gesamtdokumentation2_20260927.md'
code=p.read_text().split('Originalmitglied: `TFPT_Rekursive_Bindung_20260926/audit.py`',1)[1].split('```python\n',1)[1].split('```',1)[0]
expected='4a19f8145d5ccaef08bfdc5893f64484feae630cf61e49639df52f89f56528c4'
assert expected in [hashlib.sha256((code.rstrip()+'\n'*k).encode()).hexdigest() for k in range(4)]
ns={'__file__':str(p)}
exec(code.split('# Exact complete 35 x 35 two-module model.',1)[0],ns)
G,F,Kraw=[ns[x] for x in ['Gs','Fnum','Recovery']]
M=s.kronecker_product(G**s.Rational(1,2),G**s.Rational(1,2))
K=(M*Kraw*M.inv()).applyfunc(s.simplify)
I=s.eye(25);vI=s.eye(5).reshape(25,1);B=vI*vI.T/5
P9=s.Rational(9,4)*(K-B); E=B+P9
S=s.Matrix(25,25,lambda i,j:int(i==5*(j%5)+j//5))
Pm=(I-S)/2; P5=(I+S)/2-E
ch=lambda A:s.Matrix(25,25,lambda a,b:A[5*(a%5)+b%5,5*(a//5)+b//5])
assert K==ch(K)
assert all(P*P==P for P in [B,P9,P5,Pm])
assert K*Pm==s.zeros(25) and Pm*vI==s.zeros(25,1)
assert K.eigenvals()=={s.Integer(1):1,s.Rational(4,9):9,s.Integer(0):15}
assert ch(E).eigenvals()=={s.Integer(2):1,s.Rational(3,4):9,-s.Rational(1,4):15}
x=s.symbols('x',real=True)
assert ch(B+x*P9)==(1-9*x/4)*I/5+9*x/4*K
# Original normalized effect projectors, exact Kraus and coherent dilation.
us=[(s.sqrt(2)*G**(-s.Rational(1,2))*F[i,:].T/4).applyfunc(s.simplify) for i in range(10)]
Ps=[u*u.T for u in us]
As=[P/s.sqrt(2) for P in Ps]
assert sum(Ps,s.zeros(5))==2*s.eye(5)
assert sum((A.T*A for A in As),s.zeros(5))==s.eye(5)
assert sum((s.kronecker_product(A,A) for A in As),s.zeros(25))==K
V=s.Matrix.vstack(*As)
assert V.T*V==s.eye(5) and V.rank()==5
assert all(E*P.reshape(25,1)==P.reshape(25,1) for P in Ps)
assert s.Matrix.hstack(*(P.reshape(25,1) for P in Ps)).rank()==10
comm=Ps[0]*Ps[1]-Ps[1]*Ps[0]
assert s.trace(comm.T*comm)==s.Rational(16,81)
assert K*comm.reshape(25,1)==s.zeros(25,1)
# Full response range generates M5 as an algebra; pair products already span it.
assert s.Matrix.hstack(*((P*Q).reshape(25,1) for P in Ps for Q in Ps)).rank()==25
follow=comm*Ps[0]
recovered_follow=s.Matrix(5,5,K*follow.reshape(25,1))
assert s.trace(Ps[1]*recovered_follow)==-s.Rational(32,729)
assert K*(Ps[0]*Ps[1]).reshape(25,1)==K*(Ps[1]*Ps[0]).reshape(25,1)
# Primitive CPTP powers start only at exponent one in the support-power family.
assert s.simplify((1-s.Rational(9,4)*s.Rational(2,3))/5)==-s.Rational(1,10)
res={
 'original_code_sha256':expected,'original_checks':len(ns['checks']),
 'recovery_spectrum':{'1':1,'4/9':9,'0':15},
 'active_symmetry':'Sym^2(C5): 1+9 active, 5 kernel; wedge^2(C5): all 10 kernel',
 'support_projection_choi_spectrum':{'2':1,'3/4':9,'-1/4':15},
 'sqrt_recovery_choi_min':'-1/10',
 'kraus_count':len(As),'choi_rank':K.rank(),'stinespring_shape':list(V.shape),
 'effect_span_rank':10,'effect_product_span_rank':25,
 'commutator_norm_squared':'16/81','recovered_commutator':'zero','ordered_continuation_witness':'-32/729',
 'source_values_computed':False,'F1_complete':False}
args.output.write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
