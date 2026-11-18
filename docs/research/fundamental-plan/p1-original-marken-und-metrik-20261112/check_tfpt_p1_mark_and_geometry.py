"""Original E8 charge stabilizers and local collar metric identities."""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, product
import sympy as s

D8 = [tuple(u*(k==i)+v*(k==j) for k in range(8))
      for i,j in combinations(range(8),2)
      for u,v in product((-1,1),repeat=2)]
spin = [tuple(Q(x,2) for x in v)
        for v in product((-1,1),repeat=8) if sum(v)%4==0]
roots = set(D8+spin)
dot = lambda x,y:sum(a*b for a,b in zip(x,y))
add = lambda x,y:tuple(a+b for a,b in zip(x,y))
unit = lambda i:tuple(int(k==i-1) for k in range(8))
neg = lambda x:tuple(-t for t in x)
sub = lambda x,y:add(x,neg(y))

p1=(Q(1,2),)*8
p2=tuple(Q(-1 if i in (0,5) else 1,2) for i in range(8))
gamma=sub(p1,p2)
w=(0,2,4,-2,0,4,2,0)
Y=(-Q(1,3),)*3+(Q(1,2),)*2+(0,)*3
colour1=sub(unit(1),unit(2))
colour2=sub(unit(2),unit(3))
weak=sub(unit(4),unit(5))
assert len(roots)==240 and gamma==add(unit(1),unit(6))

def commutes_su2(root, simple):
    return dot(root,simple)==0 and add(root,simple) not in roots and sub(root,simple) not in roots

central_roots={r for r in roots if dot(Y,r)==0 and
               all(commutes_su2(r,x) for x in (colour1,colour2,weak,gamma))}
a=sub(unit(7),unit(8))
b=neg(add(unit(7),unit(8)))
h=add(a,b)
assert h==neg(add(unit(8),unit(8)))
assert {r for r in central_roots if dot(w,r)%2==0}=={a,b,neg(a),neg(b)}
assert not {r for r in central_roots if dot(w,r)%8==0}
assert all(commutes_su2(p,a) and dot(p,b)==-1 and add(p,b) in roots for p in (p1,p2))
assert all(dot(k,p)==dot(k,add(p,b)) for k in (Y,colour1,colour2,weak) for p in (p1,p2))
assert sub(add(p1,b),add(p2,b))==gamma
assert [dot(w,p)%8 for p in (p1,p2,add(p1,b),add(p2,b))]==[5,1,3,7]
assert s.Matrix([gamma]).rank()==1
print('Roots commuting with the SM gauge generators and middle-root SU(2):',len(central_roots))
print('Adding old half parity leaves exactly the four D2 roots:',sorted({r for r in central_roots if dot(w,r)%2==0}))
print('Adding the full old eighth-clock character leaves 0 roots and the connected 7D torus.')
print('Family rotation preserves all SM Cartan charges and the middle root; old eighth-clock endpoint residues (5,1) become (3,7).')

# The exact reflection-even collar metric retains both marks and the same
# reference at the boundary. Test curvature independently from Christoffels.
n,theta=s.symbols('n theta',real=True)
a_fn=s.Function('a')(n,theta)
coords=(n,theta)
metric=s.diag(1,a_fn**2)
inverse=metric.inv()
conn={(i,j,k):s.simplify(sum(inverse[i,l]*(s.diff(metric[l,k],coords[j])+s.diff(metric[l,j],coords[k])-s.diff(metric[j,k],coords[l])) for l in range(2))/2)
      for i,j,k in product(range(2),repeat=3)}
ric=s.Matrix(2,2,lambda i,j:s.simplify(sum(
     s.diff(conn[k,i,j],coords[k])-s.diff(conn[k,i,k],coords[j])
     +sum(conn[k,i,j]*conn[l,k,l]-conn[l,i,k]*conn[k,j,l] for l in range(2))
     for k in range(2))))
scalar=s.simplify(s.trace(inverse*ric))
assert s.simplify(scalar+2*s.diff(a_fn,n,2)/a_fn)==0
assert s.simplify(ric-metric*scalar/2)==s.zeros(2)
assert s.simplify(scalar*a_fn+2*s.diff(a_fn,n,2))==0
print('For the reflection-even 2D collar, R=-2 a_nn/a; Einstein tensor exactly zero; R dV is a normal boundary term.')
eps,r=s.symbols('epsilon r',real=True,positive=True)
u=(1-n**2)**4
second_coeff=s.simplify(4*r*s.integrate(s.diff(u,n,2)**2,(n,-1,1))*s.integrate(s.cos(4*theta)**2,(theta,0,2*s.pi)))
assert second_coeff>0
print('Illustrative curved reflection-even test profile: leading integral R^2 term =',second_coeff,'times epsilon^2.')
