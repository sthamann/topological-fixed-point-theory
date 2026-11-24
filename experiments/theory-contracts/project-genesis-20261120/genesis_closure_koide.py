#!/usr/bin/env python3
"""Project Genesis: n=4 closure and Koide/C3 source constraints.

Analytical identities and one retrospective empirical consistency number.
Does not derive the origin of any symmetry or mass pattern.
"""
from sympy import Matrix, eye, diag, symbols, factor, sqrt, Rational, simplify, S, N

n=symbols('n',integer=True,positive=True)
fn=n-1
cn=n+1
closure=factor(fn*fn-1-(fn+cn))
assert closure==n*(n-4)
print('Self-description closure: dim End0(F_n)-(dim F_n + dim C_n) =',closure)
print('Unique nonzero n satisfying equality:', 4)
assert all(((k-1)**2-1 == (k-1)+(k+1)) == (k==4) for k in range(2,20))

x1,x2,x3=symbols('x1 x2 x3',real=True)
x=Matrix([x1,x2,x3]); ones=Matrix([1,1,1])
P0=ones*ones.T/3
P2=eye(3)-P0
sing=simplify((P0*x).dot(P0*x))
doub=simplify((P2*x).dot(P2*x))
K=simplify(x.dot(x)/(x.dot(ones)**2))
assert simplify(K-Rational(1,3)*(1+doub/sing))==0
assert simplify((K-Rational(2,3))*3*sing-(doub-sing))==0
print('Koide: K = 1/3 (1 + ||doublet||^2/||singlet||^2) EXACT')
print('Koide K=2/3 iff singlet norm = doublet norm EXACT')

C=Matrix([[0,0,1],[1,0,0],[0,1,0]])
assert C**3==eye(3) and C*P0==P0*C and C*P2==P2*C
assert (P0*C).trace()==1 and (P2*C).trace()==-1
print('C3: family permutation real decomposition = 1 singlet + 2 real doublet')

A=diag(-1,1,-1);B=diag(1,-1,-1)
M=Matrix(3,3,symbols('z0:9'))
system=[]
for G in (A,B,C):
    system += list(M*G-G*M)
from sympy import linear_eq_to_matrix
mat,_=linear_eq_to_matrix(system, list(M))
assert mat.rank()==8 # 9 unknowns - 8 eq rank -> 1 commuting direction
print('A4 exact residual symmetry: hermitian 3x3 masses are proportional to I; no hierarchy or mixing')
assert K.subs({x1:1,x2:1,x3:1})==Rational(1,3)

m=[S('0.51099895'),S('105.6583755'),S('1776.86')]
v=Matrix([sqrt(z) for z in m])
num=N(v.dot(v)/(sum(v)**2),16)
print('Retrospective example charged-lepton pole masses (MeV): K =',num)
print('Delta K from 2/3 =',N(num-Rational(2,3),12))
assert abs(float(num)-2/3)<1e-4
print('ALL CLOSURE AND KOIDE CHECKS PASSED')
