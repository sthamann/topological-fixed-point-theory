import sympy as s
I=s.eye(4); Z=s.zeros(4); ii=s.I
X=s.Matrix([[0,1],[1,0]]); Y=s.Matrix([[0,-ii],[ii,0]]); Z2=s.diag(1,-1); I2=s.eye(2)
k=s.kronecker_product
G=[k(X,I2),k(Y,I2),k(Z2,X),k(Z2,Y),k(Z2,Z2)]
# A,F choice any anticomm. modes with proper eigen assignments. use A=G1*G2 bivector times i? and F=G1 vector
A=ii*G[0]*G[1]; F=G[0]
assert A.H==A and F.H==F and A*F==-F*A
E=[(I+s.sqrt(3)*A/2+F/2)/3,(I-s.sqrt(3)*A/2+F/2)/3,(I-F)/3]
def phi(x): return s.simplify(s.Rational(7,12)*x+sum((g*x*g for g in G),Z)/12)
def branch(rho,j):
 # sqrt(Ej) = sqrt(3/2) Ej because Ej^2=(2/3)Ej
 un=s.Rational(3,2)*E[j]*rho*E[j]
 p=s.simplify(s.trace(un))
 return s.simplify(un/p),p
rho=I/4
# first outcome i, then Phi, then j, then Phi, then k
B=s.Matrix([[13,1,4],[1,13,4],[4,4,10]])/18
J=s.ones(3)/3
M=s.simplify((B+J)/2)
Gm=s.Matrix(3,3,lambda i,j:s.trace(E[i]*E[j])/4)
M_first=s.zeros(3)
third={}
for i in range(3):
 r1,p1=branch(rho,i)
 for j in range(3):
  pj=s.simplify(s.trace(E[j]*phi(r1)))
  M_first[i,j]=pj
  # normalized state after outcome j from phi(r1)
  r2,_=branch(phi(r1),j)
  third[(i,j)]=[s.simplify(s.trace(E[kk]*phi(r2))) for kk in range(3)]
print('G=',Gm)
print('3G=',3*Gm)
print('B=',B)
print('Mformula=',M)
print('Mfirst=',M_first)
print('equal=',M_first==M)
for key,val in third.items():
 if s.Matrix(val)-s.Matrix(list(M.row(key[1]))) != s.zeros(3,1):
  print('counterexample previous/current=',key,'third=',val,'autonomous row=',list(M.row(key[1])))
  break
print('all_third_autonomous=',all(s.Matrix(v)==s.Matrix(list(M.row(k[1]))) for k,v in third.items()))
print('third00=',third[(0,0)])
print('third10=',third[(1,0)])
# sqrt(Ej) Phi(Ek) sqrt(Ej) = M_jk Ej, sqrt(Ej)=sqrt(3/2)Ej
print('operator_markov_identity=',all(s.simplify(s.Rational(3,2)*E[j]*phi(E[k])*E[j]-M[j,k]*E[j])==Z for j in range(3) for k in range(3)))
print('J3=',J)
