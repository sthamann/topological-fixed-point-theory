"""Direct differential and original H0(O1) checks; no physical profile selected."""
from pathlib import Path
import json
import sympy as s
x,y=s.symbols('x y',real=True)
sig=s.Function('sigma')(x,y);f=s.exp(-sig/2);W=1/f;h=s.diff(sig,x)
I=s.I;zero=s.zeros(2,1);c=s.Matrix([[0,1],[1,0]])
def matrix(prefix):return s.Matrix(2,2,lambda i,j:s.Function(prefix+str(i)+str(j))(x,y))
V=matrix('v');u=s.Matrix([s.Function('u0')(x,y),s.Function('u1')(x,y)])
def B(q):return -I*c*q.diff(y)+V*q
def Bs(q):return f*B(f*q)
def simp(q):return q.applyfunc(s.simplify)
checks=[]
def check(name,value):checks.append(dict(name=name,passed=bool(value)))
check('first order scalar commutator including full potential',simp(B(h*u)-h*B(u)+I*c*s.diff(h,y)*u)==zero)
lhs=s.diff(Bs(u),x)-Bs(u.diff(x))
rhs=f*V.diff(x)*(f*u)-h*Bs(u)/2-f*B(h*f*u)/2
check('full congruence derivative retains both ordered factors',simp(lhs-rhs)==zero)
check('unscaled normal operator has exact multiplier residual',simp(B(W*u)-W*B(u)+I*c*s.diff(W,y)*u)==zero)
check('whole outer insertion congruence identity',simp((Bs(u)+B(u))-f*(B(f*u)+W*B(W*f*u)))==zero)
check('outer insertion defect retains first order commutator',simp(B(u)-f*B(f*u)-(1-f*f)*B(u)+f*(B(f*u)-f*B(u)))==zero)
# Exact polynomial identity for the WHOLE original span(1,z).
a,b,d,e,zeta=s.symbols('a b d e zeta')
poly=s.Poly((a+b*zeta)*zeta-(d+e*zeta),zeta)
check('all smooth scalar multipliers preserving original O1 span are constant',s.solve(poly.all_coeffs(),[b,d,e],dict=True)==[{b:0,d:0,e:a}])
# General product-chart differential identity,including varying symbol and
# noncommuting horizontal/vertical connections; no compatibility imposed.
C=matrix('c');Amu=matrix('m');An=matrix('n');Phi=matrix('p')
def N(q):return q.diff(y)+An*q
def H(q):return q.diff(x)+Amu*q
def Bgen(q):return C*N(q)+Phi*q
dc=C.diff(x)+Amu*C-C*Amu
Fmix=An.diff(x)-Amu.diff(y)+Amu*An-An*Amu
dPhi=Phi.diff(x)+Amu*Phi-Phi*Amu
check('full horizontal normal derivative includes symbol and mixed curvature',simp(H(Bgen(u))-Bgen(H(u))-dc*N(u)-C*Fmix*u-dPhi*u)==zero)
# Clifford algebra check in one local two-dimensional orthonormal frame.
p,q=s.symbols('p q',real=True);pauli2=s.Matrix([[0,-I],[I,0]])
cl=p*c+q*pauli2
check('real nonzero normal covector has invertible Clifford action',s.simplify(cl.det()+p*p+q*q)==0 and simp(cl*cl-(p*p+q*q)*s.eye(2))==s.zeros(2))
result=dict(date='2026-10-31',result='R052',scope='Original fixed H0(O1) multiplier class and common-core first-order producer/placement identities. No chosen normal geometry,horizontal connection,physical profile or exact dynamical reduction.',all_passed=all(t['passed'] for t in checks),checks=checks)
path=Path(__file__).with_name('NORMAL-SOURCE-PLACEMENT.json');path.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(all_passed=result['all_passed'],checks=len(checks))))
assert result['all_passed']
