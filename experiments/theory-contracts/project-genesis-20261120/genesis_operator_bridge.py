#!/usr/bin/env python3
"""Exact representation dictionary: V4 theta 8 ↔ traceless M3(C).

The dictionary is an isomorphism of the stated finite V4 x deck modules.
No geometric/physical identification of operators, couplings or clock is proved.
"""
from sympy import Matrix, zeros, eye, diag, I, symbols, simplify, sqrt, Rational
from itertools import combinations

# Three-dimensional real tetrahedral representation of V4 ⋊ C3 = A4.
a=diag(-1,1,-1); b=diag(1,-1,-1)
c=Matrix([[0,0,1],[1,0,0],[0,1,0]])
E=lambda i,j: Matrix(3,3,lambda k,l: int(k==i and l==j))
D1=diag(1,-1,0); D2=diag(1,1,-2)
S=[E(i,j)+E(j,i) for i,j in [(0,1),(0,2),(1,2)]]
K=[E(i,j)-E(j,i) for i,j in [(0,1),(0,2),(1,2)]]
Q=[D1,D2]+S       # symmetric-traceless 5 = carrier
F=K                # antisymmetric 3 = family
basis=Q+F
assert all(M.trace()==0 for M in basis)
assert all(Qi==Qi.T for Qi in Q) and all(Fi==-Fi.T for Fi in F)
W=Matrix.hstack(*[Matrix(9,1,list(M)) for M in basis]);assert W.rank()==8

def coords(M):
    ans, params=W.gauss_jordan_solve(Matrix(9,1,list(M)))
    assert not params.rows
    return ans

def act(g):
    return Matrix.hstack(*[coords(g*M*g.T) for M in basis])

def parity(M): return M.T

A,B,C=act(a),act(b),act(c)
Deck=diag(*([1]*5+[-1]*3))
assert A**2==eye(8) and B**2==eye(8) and C**3==eye(8)
assert A*B==B*A and C*A*C.inv() in [A,B,A*B] and C*B*C.inv() in [A,B,A*B]
assert all(G*Deck==Deck*G for G in [A,B,C])
print('M3 traceless V4 traces: total', [G.trace() for G in [eye(8),A,B,A*B]])
print('Carrier Sym0 traces:', [G[:5,:5].trace() for G in [eye(8),A,B,A*B]])
print('Family so3 traces:', [G[5:,5:].trace() for G in [eye(8),A,B,A*B]])
print('C3 trace on carrier/family:',C[:5,:5].trace(),C[5:,5:].trace())
assert [G.trace() for G in [eye(8),A,B,A*B]]==[8,0,0,0]
assert [G[:5,:5].trace() for G in [eye(8),A,B,A*B]]==[5,1,1,1]
assert [G[5:,5:].trace() for G in [eye(8),A,B,A*B]]==[3,-1,-1,-1]
assert C[:5,:5].trace()==-1 and C[5:,5:].trace()==0

# Verify the compact real su3 Lie algebra symmetry pair 3 + 5:
# so(3) + i Sym0(3). A symmetry-respecting algebraic bridge, not physical gluons.
anti_basis=K+[I*M for M in Q]
for M in anti_basis:
    assert M.conjugate().T==-M and M.trace()==0
for ii,jj in combinations(range(8),2):
    x,y=anti_basis[ii],anti_basis[jj]
    z=x*y-y*x
    assert z.conjugate().T==-z
    cc,params=Matrix.hstack(*[Matrix(9,1,list(M)) for M in anti_basis]).gauss_jordan_solve(Matrix(9,1,list(z)))
    assert not params.rows
    if ii<3 and jj<3: assert all(cc[t]==0 for t in range(3,8))
    if ii>=3 and jj>=3: assert all(cc[t]==0 for t in range(3,8))
    if ii<3 and jj>=3: assert all(cc[t]==0 for t in range(3))
print('su3 symmetric-pair closure: all 28 independent commutator pairs passed')

# Build explicit rational intertwiner for the torus 8 mode representation in GESAMT2.
# T1(k)=(-1)^k, T2 maps k to k+4, deck maps k to -k modulo 8.
T1=diag(*[(-1)**k for k in range(8)])
T2=zeros(8); d=zeros(8)
for k in range(8): T2[(k+4)%8,k]=1;d[(-k)%8,k]=1
assert T1*T2==T2*T1 and d*T1==T1*d and d*T2==T2*d

# Construct common eigenvectors for (T1,T2,deck), in the same order as basis.
eigvals=[]
for M in basis:
    eigvals.append(((a*M*a.T-M)==zeros(3), (b*M*b.T-M)==zeros(3), M==M.T))
# Convert boolean to sign ±1 and project to common eigenspaces.
signs=[tuple(1 if v else -1 for v in t) for t in eigvals]
cols=[]
for alpha,beta,delta in signs:
    P=(eye(8)+alpha*T1)*(eye(8)+beta*T2)*(eye(8)+delta*d)/8
    v=next((P[:,j] for j in range(8)
            if P[:,j]!=zeros(8,1) and Matrix.hstack(*(cols+[P[:,j]])).rank()==len(cols)+1),None)
    assert v is not None
    cols.append(v)
V=Matrix.hstack(*cols)
assert V.det()!=0
J=V.inv()
assert J*T1*J.inv()==A and J*T2*J.inv()==B and J*d*J.inv()==Deck
print('EXPLICIT INTERTWINER J from theta 8 to End0(C3): determinant =',J.det())
print('Verified J T1 J^-1=A, J T2 J^-1=B, J deck J^-1=transpose parity')

# Four ordered punctures have cross ratio lambda; order-3 geometrical symmetry
# requires lambda^2-lambda+1=0 (equianharmonic) while square is harmonic.
lam=symbols('lam')
j=256*(1-lam+lam**2)**3/(lam**2*(1-lam)**2)
assert simplify(j.subs(lam,-1))==1728
assert simplify(j.subs(lam,2))==1728
assert simplify(j.subs(lam,Rational(1,2)))==1728
assert simplify(j.subs(lam,(1+I*sqrt(3))/2))==0
assert simplify(j.subs(lam,(1-I*sqrt(3))/2))==0
print('GEOMETRIC OBSTRUCTION: harmonic square j=1728; order-3 fixed cross ratio j=0')
print('ALL BRIDGE TESTS PASSED (scope: representations, not physical origin)')
