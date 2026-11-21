#!/usr/bin/env python3
"""Exact finite witness: one ancilla qubit permits nontrivial global modular
action with precisely tracial reduced M5 state. This does NOT select rho."""
from sympy import Matrix,eye,diag,Rational,log,kronecker_product,simplify,zeros
Z5=diag(1,-1,0,0,0); Z2=diag(1,-1)
eps=Rational(1,2)  # synthetic control only
rho=(eye(10)+eps*kronecker_product(Z5,Z2))/10
assert rho.trace()==1 and all(v>0 for v in rho.diagonal())
rho5=Matrix(5,5,lambda i,j:sum(rho[2*i+a,2*j+a] for a in range(2)))
assert rho5==eye(5)/5
K=diag(*[log(x) for x in rho.diagonal()])
E=Matrix.zeros(5);E[0,1]=1
A=kronecker_product(E,eye(2))
comm=K*A-A*K
assert comm!=zeros(10)
assert all(simplify(v)==0 for v in comm-log(3)*kronecker_product(E,Z2))
print("EXACT: rho_5=I5/5, [log(rho),E01 x I2]=log(3)(E01 x Z2) !=0")
print("Outside N=M5 x I2, modular flow changes N by creating a memory operator.")
print("With memory dimension1 and tracial N, modular flow is identity on N.")
print("eps is FREE, so this is a witness, NOT P1 or selected physical time.")
print("PASS")
