"""Symbolic tangent identities for the original R035 normal-kernel map.
Finite matrices are coordinate checks of universal identities, not physical models.
"""
from pathlib import Path
import json
import sympy as s
checks=[]
def need(name,condition):
 good=bool(condition);checks.append(dict(name=name,passed=good))
 if not good:raise AssertionError(name)
def eq(a,b):return all(s.simplify(v)==0 for v in a-b)
def herm(prefix):
 a,b,c,d=s.symbols(prefix+'0:4',real=True)
 return s.Matrix([[a,c+s.I*d],[c-s.I*d,b]])
def mat(prefix,n,m):
 ar=s.symbols(prefix+'r0:'+str(n*m),real=True);ai=s.symbols(prefix+'i0:'+str(n*m),real=True)
 return s.Matrix(n,m,[x+s.I*y for x,y in zip(ar,ai)])
I=s.eye(2);T=s.Matrix.vstack(I,s.zeros(2));P=T*T.H;Q=s.eye(4)-P
A0=s.I*herm('a');A1=s.I*herm('b');C0=mat('c',2,2);C1=mat('d',2,2)
J0=T*A0+s.Matrix.vstack(s.zeros(2),C0);J1=T*A1+s.Matrix.vstack(s.zeros(2),C1)
need('unitary-frame tangent constraints',eq(T.H*J0+J0.H*T,s.zeros(2)))
need('whole moving-frame derivative splits orthogonally',eq(T*(T.H*J0)+Q*J0,J0))
need('connection antihermitian',eq(A0.H,-A0))
h=mat('h',2,1);dh=mat('v',2,1)
v=J0*h+T*dh
need('whole kinetic restriction retains positive normal Gram',eq(v.H*v,(dh+A0*h).H*(dh+A0*h)+h.H*C0.H*C0*h))
curv=J0.H*J1-J1.H*J0+A0*A1-A1*A0
need('moving-kernel curvature sign and full matrix order',eq(curv,C0.H*C1-C1.H*C0))
need('off-diagonal Gram blocks need not be Hermitian',s.I*s.diag(1,-1)!=(s.I*s.diag(1,-1)).H)

V=s.Matrix([[0,s.I],[s.I,0]])
need('passive frame covariance of connection',eq((T*V).H*(J0*V),V.H*A0*V))
need('passive frame covariance of normal Gram',eq((Q*J0*V).H*(Q*J0*V),V.H*C0.H*C0*V))
Gframe=s.I*herm('g')
need('varying unitary-frame derivative has affine connection term',eq(T.H*(J0+T*Gframe),A0+Gframe))
need('varying frame does not change orthogonal derivative',eq(Q*(J0+T*Gframe),Q*J0))

# Arbitrary Hermitian invertible paired-normal block and tangent off-diagonal.
M=herm('m');X0=mat('x',2,2);X1=mat('y',2,2);Minv=M.inv()
B=s.diag(s.zeros(2),M)
dB=s.Matrix.vstack(s.Matrix.hstack(s.zeros(2),X0.H),s.Matrix.hstack(X0,s.zeros(2)))
Coperator=-Minv*X0;Joperator=s.Matrix.vstack(A0,Coperator)
need('differentiated actual normal kernel equation',eq(dB*T+B*Joperator,s.zeros(4,2)))
need('reduced normal inverse gives whole positive kinetic Gram',eq(Coperator.H*Coperator,X0.H*Minv*Minv*X0))
need('reduced normal inverse gives mixed ordered Gram',eq(Coperator.H*(-Minv*X1),X0.H*Minv*Minv*X1))

# Noncommuting Gram square root, with arbitrary positive branch understood.
r0,r1,r2=s.symbols('r0 r1 r2',real=True);R=s.Matrix([[r0,r1],[r1,r2]])
dR=herm('t');G=R*R;dG=R*dR+dR*R;Ri=R.inv()
A=s.Rational(1,2)*Ri*dG*Ri-dR*Ri
need('noncommuting square-root derivative keeps its Sylvester order',eq(A,s.Rational(1,2)*(Ri*dR-dR*Ri)))
need('resulting scaled connection is antihermitian',eq(A.H,-A))

# Whole positive scalar multiplier class: a compressed multiplication algebra.
# S is a passive orthonormal coordinate frame; no S2 metric/profile is selected.
S=s.Matrix([[1/s.sqrt(2),1/s.sqrt(6)],[-1/s.sqrt(2),1/s.sqrt(6)],[0,-2/s.sqrt(6)]])
ww=s.symbols('w0:3',positive=True);v0=s.symbols('v0:3',real=True);v1=s.symbols('q0:3',real=True)
W=s.diag(*ww);dW0=s.diag(*v0);dW1=s.diag(*v1)
GG=S.H*W*W*S;Ginv=GG.inv();PP=W*S*Ginv*S.H*W;QQ=s.eye(3)-PP
need('scalar-multiplier compressed projector is exact',eq(PP*PP,PP))
need('whole normal covariance including mixed matrix order',eq(S.H*dW0*QQ*dW1*S,S.H*dW0*dW1*S-(S.H*dW0*W*S)*Ginv*(S.H*W*dW1*S)))
need('uniform normal multiplier derivative has zero covariance',eq(S.H*W*QQ*W*S,s.zeros(2)))

out=dict(date='2026-10-31',result='R051',scope='Original R035 moving normal kernel with supplied metric connection and smooth family. Universal tangent/positive-multiplier identities,not selection of physical geometry or exact band decoupling.',all_passed=all(x['passed'] for x in checks),checks=checks)
Path(__file__).with_name('MOVING-NORMAL-KERNEL.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(checks=len(checks),all_passed=out['all_passed'])))
