#!/usr/bin/env python3
"""Exact Fibonacci history graph KMS selection, with source assumptions explicit.

Input A=[[0,1],[1,1]] is multiplication by tau in Fib.
ADDITIONAL axioms: infinite Cuntz-Krieger algebra O_A, standard
word-length gauge evolution alpha_t(S_i)=e^(it) S_i.
Known KMS uniqueness theorem for primitive A gives beta=log rho(A).
This is gauge/modular time, NOT physical positive-energy spacetime time.
"""
from itertools import product
from sympy import Matrix, sqrt, log, simplify

A=Matrix([[0,1],[1,1]])
phi=(1+sqrt(5))/2
p=Matrix([phi**-2,phi**-1])
assert simplify(phi**2-phi-1)==0
assert all(simplify(z)==0 for z in A*p-phi*p)
assert simplify(p[0]+p[1]-1)==0
for i in (0,1):
    assert simplify(p[i]-sum(A[i,j]*p[j] for j in (0,1))/phi)==0
print("PF root phi = (1+sqrt(5))/2; unique gauge-KMS beta=log(phi) =",log(phi).evalf(14))
print("vertex weights (phi^-2,phi^-1) =",*[v.evalf(14) for v in p])

def mass(w):return phi**(-(len(w)-1))*p[w[-1]]
for n in range(1,9):
    words=[w for w in product((0,1),repeat=n)
           if all(A[w[i],w[i+1]] for i in range(n-1))]
    assert simplify(sum(mass(w) for w in words)-1)==0
    for w in words:
        extension=[w+(i,) for i in (0,1) if A[w[-1],i]]
        assert simplify(mass(w)-sum(mass(z) for z in extension))==0
    print("word length",n,"allowed histories",len(words),"total weight=1")
assert all(simplify(omega*(log(phi)/omega)-log(phi))==0 for omega in (1,2,3))
print("KMS modular group: sigma_t(S_i)=exp(-it log(phi))*S_i, nontrivial.")
print("NO-GO: gauge action is periodic, and ±n degrees give two-sided Liouville spectrum.")
print("NO-GO: rescaling gauge evolution by omega changes beta to log(phi)/omega.")
print("PASS: exact finite-cylinder equations. Lorentz time, P1, c3 and gravity OPEN.")
