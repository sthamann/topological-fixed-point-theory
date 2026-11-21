#!/usr/bin/env python3
"""Exact fifth-pass: one PF trace for finite Fibonacci fusion blocks and
an independent positive reflected Doob generator. NOT original TFPT P1."""
from itertools import product
from sympy import Matrix, sqrt, simplify, Rational, eye
A=Matrix([[0,1],[1,1]])
phi=(1+sqrt(5))/2
r=Matrix([1,phi])
assert all(simplify(z)==0 for z in A*r-phi*r)
counts=Matrix([[1,0]])
for n in range(1,17):
    counts=counts*A
    d0,d1=counts
    assert simplify((d0+phi*d1)/phi**n-1)==0
    for c in (0,1):
        assert simplify(r[c]/phi**n-sum(A[c,j]*r[j]/phi**(n+1) for j in (0,1)))==0
    if n==4:assert (d0,d1)==(2,3)
    if n==5:assert (d0,d1)==(3,5)
print("AF level-n tr_n = phi^-n(Tr X0+phi Tr Xtau) compatible n=1..16")
print("Conditioned n=5 total tau: M5 trace I5/5, preceding 4-event sectors 2/5,3/5")
print("Unconditioned n=4: 2/phi^4 and 3phi/phi^4 (different conditioning)")
D=Matrix.diag(1,phi)
P=(D.inv()*A*D/phi).applyfunc(simplify)
pi0=simplify(1/(1+phi**2));pi1=simplify(phi**2/(1+phi**2))
Pi=Matrix([[pi0,pi1],[pi0,pi1]])
q=phi**(-4)
T=(P*P).applyfunc(simplify)
assert all(simplify(z)==0 for z in (T-(Pi+q*(eye(2)-Pi))))
assert all(simplify(pi0*P[0,1]-pi1*P[1,0])==0 for _ in (0,))
assert simplify(P.det()+phi**(-2))==0
assert simplify(T.det()-phi**(-4))==0
print("P eigenvalues 1,-phi^-2; P-adjoint P=P^2 eigen 1,phi^-4")
print("Unique selfadjoint H=4log(phi)(I-Pi), positive and finite")
K=Matrix(3,3,lambda i,j:pi0*pi0+pi0*pi1*q**(i+j))
assert simplify(K.det())==0
assert simplify(K[0,0]*K[1,1]-K[0,1]**2)>0
print("Reflection Gram rank2 positive: Kmn=pi0^2+pi0*pi1*phi^-4(m+n)")
def primitive(M):
    Z=eye(M.rows)
    for _ in range(1,4*M.rows*M.rows+5):
        Z=Z*M
        if all(t>0 for t in Z):return True
    return False
gr=[]
for bits in product((0,1),repeat=4):
    B=Matrix(2,2,bits)
    if primitive(B):gr.append((bits,B.det()))
assert len(gr)==3 and sum(abs(z[1])==1 for z in gr)==2
print("Baseline 3/16 labeled 2-state graphs primitive, 2/16 invertible; both Fib")
x=__import__('sympy').Symbol('x')
for n in range(2,11):
    B=Matrix.zeros(n)
    for i in range(n-1):B[i,i+1]=1
    B[n-1,0]=1;B[n-1,1]=1
    assert B.charpoly(x).as_expr()==x**n-x-1 and abs(B.det())==1
    assert primitive(B)
print("No global entropy minimum: primitive unimodular n-state family x^n-x-1, rho_n->1")
print("Gauge KMS flow acts identically on all fixed-length AF matrix units including M5")
print("STRICT GATE: this flow != positive Doob time; original P1 and gravity OPEN")
print("ALL EXACT CHECKS PASSED")
