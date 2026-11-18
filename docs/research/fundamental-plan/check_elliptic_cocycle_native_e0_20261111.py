"""Exact old elliptic normal spin times the existing E0/C3 family lifts."""
import sympy as s
w=(-1+s.I*s.sqrt(3))/2
simple=lambda x:x.applyfunc(s.simplify)
D=s.Matrix([[0,-1],[1,0]])
Q=s.diag(1,-1)
A=s.Matrix([[0,0,-1],[0,1,0],[1,0,0]])
A4=s.diag(A,1)
F=s.diag(w,w,-w,-1)
assert s.simplify(F.det()-1)==0
K4=simple(A4*F*A4.inv()*F.inv())
assert K4==s.diag(-1,1,-1,1)
total=simple(s.kronecker_product(Q,A4)*s.kronecker_product(D,F)*
             s.kronecker_product(Q,A4).inv()*s.kronecker_product(D,F).inv())
assert total==s.kronecker_product(s.eye(2),s.diag(1,-1,1,-1))
P=s.diag(1,0,1)
Pc=s.Matrix([[0,0,1],[1,0,0],[0,1,0]])
assert Pc*P!=P*Pc
assert s.Matrix.hstack(*[Pc**i*P for i in range(3)]).rank()==3
assert P*A==A*P and P*F[:3,:3]==F[:3,:3]*P
print("Original C3 quarter and E0 compensator on q four-slot family:", list(K4.diagonal()))
print("With elliptic normal spin:", list((s.eye(4)-K4).diagonal()/2),
      "are the family slots with trivial total commutator")
print("Cancellation on L0,L2; no cancellation on L1 or the fourth line.")
print("L0,L2 fails downstairs: cubic orbit of the cancellation plane spans the triplet.")
