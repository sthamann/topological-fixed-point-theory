#!/usr/bin/env python3
"""Exact checks for the Hodge channel versus a continuous identity start.

Uses the original epsilon amplitudes. This is a finite operator result, not a
physical identification of time, a source instrument, or a continuum limit.
"""
import importlib.util
from pathlib import Path
import sympy as s

source = Path(__file__).resolve().parents[1] / 'hodge_complement_20261028' / 'hodge_complement_channel.py'
spec = importlib.util.spec_from_file_location('original_hodge', source)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
C = [s.Matrix(a.tolist()) for a in module.hodge_amplitudes()]
checks = []
def check(ok, description):
    if not ok:
        raise AssertionError(description)
    checks.append(description)

z = s.symbols('z0:5')
w = s.symbols('w0:5') # w is the conjugate of z on the actual complex slice
V = sum((z[i]*C[i] for i in range(5)), s.zeros(10))
VH = sum((w[i]*C[i].H for i in range(5)), s.zeros(10))
r = sum(z[i]*w[i] for i in range(5))
P = VH*V
check(all(s.expand(x)==0 for x in P*P-r*P), 'all complex linear combinations have normalized projection square')
check(s.expand(s.trace(P)-6*r)==0, 'trace is six times squared coefficient norm')
check(all((C[i].H*C[j]).trace()==6*int(i==j) for i in range(5) for j in range(5)), 'orthogonality of five actual epsilon amplitudes')

def vec(A):
    return s.Matrix([A[i,j] for i in range(10) for j in range(10)])
P5=sum((vec(A)*vec(A).H for A in C),s.zeros(100))/6
check(P5*P5==P5 and s.trace(P5)==5, 'P5 is rank five projector')
J=sum((vec(A)*vec(A).H for A in C),s.zeros(100))/3
check(J==2*P5 and s.trace(J)==10, 'actual Hodge Choi operator is 2 P5')
# C0 is a signed partial permutation. Complete it to an exact unitary.
in_missing=[j for j in range(10) if C[0][:,j]==s.zeros(10,1)]
out_missing=[i for i in range(10) if C[0][i,:]==s.zeros(1,10)]
U=C[0].copy()
for i,j in zip(out_missing,in_missing): U[i,j]=1
check(U.H*U==s.eye(10) and U*U.H==s.eye(10), 'unitary completion of the actual rank six partial isometry')
v=vec(U)
check((v.H*P5*v)[0]==6, 'overlap upper bound six is attained')
check(s.trace(P5*J/10)==1 and (v.H*P5*v)[0]/10==s.Rational(3,5), 'normalized Choi test distinguishes channel and unitary by 2/5')
q,t=s.symbols('q t',real=True)
block=s.Matrix([[s.Rational(1,5)-q, -t],[-t,q-1]])
discriminant=s.expand(s.trace(block)**2-4*block.det()).subs(t*t,q*(1-q)).expand()
check(s.expand(discriminant-(s.Rational(36,25)-s.Rational(4,5)*q))==0, 'two-dimensional Choi discriminant after t squared equals q(1-q)')
e=s.zeros(10,1); e[in_missing[0]]=1
check(all((e.H*U.H*A*e)[0]==0 for A in C), 'a concrete state perfectly distinguishes this unitary from the Hodge output')
# Unitarly identify H* with H using U*. The next checks prove trace
# preservation of the one- and two-jump terms in the Poisson embedding.
K=[U.H*A/s.sqrt(3) for A in C]
check(sum((a.H*a for a in K),s.zeros(10))==s.eye(10), 'one-jump channel is trace preserving')
K2=[a*b for a in K for b in K]
check(sum((a.H*a for a in K2),s.zeros(10))==s.eye(10), 'two-jump channel is trace preserving')
check(s.Matrix.hstack(*(a.reshape(100,1) for a in C)).rank()==5, 'five epsilon amplitudes span a five-dimensional space')
check(s.Matrix.hstack(*(a.reshape(100,1) for a in C),U.reshape(100,1)).rank()==6, 'unitary no-event amplitude is outside strict epsilon span')
print('PASS',len(checks),'exact checks')
for item in checks: print(' ',item)
