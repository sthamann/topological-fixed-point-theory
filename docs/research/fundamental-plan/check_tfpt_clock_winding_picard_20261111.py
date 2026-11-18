"""Compare proposed winding and Picard operators with the original Q4 code."""
from pathlib import Path
import runpy
from itertools import combinations, product
from collections import Counter, deque
import sympy as s
from sympy.combinatorics import Permutation

P=Path(__file__).with_name('check_tfpt_q4_picard_selection_20261111.py')
state=runpy.run_path(str(P))
W=state['W']; K=state['K']; gram=state['G']; B=state['B']
Wslot=state['sv']; labels=state['native_labels']
I=s.I
simp=lambda m:m.applyfunc(s.simplify)
vec=lambda m:s.Matrix(25,1,list(m))
def perm_matrix(p):
    return s.Matrix(len(p),len(p),lambda i,j:int(i==p[j]))
def comp(p,q):
    return tuple(p[q[k]] for k in range(len(p)))
e=s.eye(5)
C=perm_matrix(tuple((i+1)%5 for i in range(5)))
G=perm_matrix(tuple(2*i%5 for i in range(5)))
Wind=[sum((e[:,i]*e[(i+d)%5,:]-e[:,(i+d)%5]*e[i,:] for i in range(5)),s.zeros(5))
      for d in (1,2)]
w1,w2=Wind
assert C*w1*C.T==w1 and C*w2*C.T==w2
assert G*w1*G.T==w2 and G*w2*G.T==-w1
a=s.Matrix(5,5,s.symbols('x:25'))
assert len(s.linear_eq_to_matrix(list(a+a.T)+list(C*a*C.T-a),tuple(a))[0].nullspace())==2
assert K*vec(gram.inv())==vec(gram.inv())
assert all((K*vec(simp(W*A*W.T))).applyfunc(s.simplify)==s.zeros(25,1) for A in Wind)
omega=(vec(w1)+I*vec(w2))/s.sqrt(20)
assert (s.kronecker_product(G,G)*omega).applyfunc(s.simplify)==-I*omega
Sw=state['Sw']
candidate=s.Rational(4,9)*(s.eye(25)-Sw)/2+s.Rational(5,9)*omega*omega.H
assert candidate*vec(e)==s.zeros(25,1)
assert candidate*omega==omega
assert ((s.eye(25)-Sw)*K)==s.zeros(25) and (s.eye(25)+Sw)*candidate==s.zeros(25)
print('Native windings: C5 fixes, G rotates, +/- i lines are exact; original Q4 kills both and fixes identity.')
print('The proposed positive matrix annihilates the order unit: as a nonzero End(C5) map it cannot even be positive.')
print('Original and proposed operator images are orthogonal in the actual same slot coordinates.')

# Recover both windings from genuine ordered Q4 effect products.
projs=[simp(Wslot[:,j]*Wslot[:,j].T) for j in range(10)]
for i,j in combinations(range(10),2):
    A=simp(projs[i]*projs[j]-projs[j]*projs[i])
    V=simp(sum((C**k*A*(C.T)**k for k in range(5)),s.zeros(5))/5)
    if V!=s.zeros(5):break
v1,v2=[s.simplify(s.trace(v.T*V)/10) for v in Wind]
assert simp(V-v1*w1-v2*w2)==s.zeros(5)
assert s.simplify(v1**2+v2**2)>0
assert simp(G*V*G.T+I*V)!=s.zeros(5)
assert s.Matrix.hstack(vec(V),vec(simp(G*V*G.T))).rank()==2
print('Native Bell indices',i,j,'Q4 effect indices',state['align'][0][i],state['align'][0][j],
      'give C5-averaged commutator coefficients',v1,v2)
print('This native ordered word and its G image generate both windings.')
# Fixed 3+2 gauge partition on the slot candidate. Independent of choice of
# pair: neither real nor complex Wind line commutes with the native Y diagonal.
spectra=Counter()
for weak in combinations(range(5),2):
    Y=s.diag(*[s.Rational(1,2) if j in weak else -s.Rational(1,3) for j in range(5)])
    L=[Y*A-A*Y for A in Wind]
    assert s.Matrix.hstack(*(vec(a) for a in L)).rank()==2
    spectra[tuple(s.trace(a.T*a) for a in L)]+=1
print('All ten fixed slot charges: no neutral wind combination; commutator norm-squared census',dict(spectra))
# Normalization of the actual clock algebra excludes inversion of the
# multiplier under any operator preserving C5 as subgroup, including anti-linear.
assert all(pow(2,1,5)*j%5!=pow(2,-1,5)*j%5 for j in range(1,5))
print('Reversing G while normalizing the same C5 is inconsistent with G C5 G^-1=C5^2.')


# The setwise stabilizer of any weak duad in AGL(1,F5) has two elements.
# In either complex character it gives exactly the same phase.
aff=[(b,2**j%5,j) for b,j in product(range(5),range(4))]
for weak in combinations(range(5),2):
    keep=[(b,a,j) for b,a,j in aff if {(a*k+b)%5 for k in weak}==set(weak)]
    assert len(keep)==2 and sorted(j for b,a,j in keep)==[0,2]
    assert all(I**j==(-I)**j for b,a,j in keep)
print('Any affine return to a fixed weak duad has equal phases on both winding lines, +1 or -1.')

# Picard group on the submitted twenty points. Exact projective transformations.
rt=s.sqrt(3); z=(-1+I*rt)/2
root=[I**k for k in range(4)]
pic=[(s.Integer(0),s.Integer(1),s.Integer(0))]+[(0,-z**k,1) for k in range(3)]
pic += [(x,0,1) for x in root]+[(rt*x,2*z**k,1) for k,x in product(range(3),root)]
def project(p):
    j=2 if s.simplify(p[2])!=0 else 1
    return tuple(s.simplify(s.expand(s.sympify(x)/s.sympify(p[j]))) for x in p)
pic=[project(p) for p in pic]
assert len(set(pic))==20 and all(s.simplify(y**3*h+h**4-x**4)==0 for x,y,h in pic)
def sigma(x,y,h): return x,z*y,h
def tau(x,y,h): return I*x,y,h
def rho(x,y,h): return rt*x,2*h-y,y+h
def real(x,y,h): return s.conjugate(x),s.conjugate(y),s.conjugate(h)
def pic_index(q):
    return next(i for i,p in enumerate(pic) if all(s.simplify(a-b)==0 for a,b in zip(p,q)))
gens=[tuple(pic_index(project(fn(*p))) for p in pic) for fn in [sigma,tau,rho,real]]
unit=tuple(range(20))
def group(gs):
    todo=[unit]; seen={unit}
    for p in todo:
        for q in gs:
            pq=comp(p,q)
            if pq not in seen:
                todo.append(pq);seen.add(pq)
    return seen
def order(p):
    return Permutation(list(p)).order()
def cycle(p):
    perm=Permutation(list(p))
    return tuple(sorted([1]*(20-sum(map(len,perm.cyclic_form)))+list(map(len,perm.cyclic_form))))
hol=group(gens[:3]);full=group(gens)
assert len(hol)==48 and len(full)==96
for j in (0,4):
    orbit={p[j] for p in hol}
    assert len(orbit)==(4 if j==0 else 16)
c3={frozenset(k for k,j in enumerate(p) if k==j) for p in hol if order(p)==3}
assert len(c3)==4 and set.union(*(set(k) for k in c3))==set(range(20))
assert all(len(d)==5 and len(d&set(range(4)))==1 for d in c3)
assert cycle(gens[1])==(1,1,1,1,4,4,4,4)
print('Picard 20 = 4+16 orbits; four disjoint C3 fixed sets of five, each with one special point.')
print('Holomorphic order-four cycle types',dict(Counter(cycle(p) for p in hol if order(p)==4)))
full_quarters=[p for p in full if order(p)==4 and cycle(p)==(4,)*5]
print('Fixed-point-free extended order-four permutations:',len(full_quarters),
      'holomorphic:',sum(p in hol for p in full_quarters))
assert not any(p in hol for p in full_quarters)
quarter_marks=[p for p in full if len({p[0],p[p[0]],p[p[p[0]]],0})==4]
print('Full Picard lifts of a four-cycle on the four special marks:',dict(Counter((order(p),cycle(p)) for p in quarter_marks)))
assert not any(order(p)==5 for p in full)
# Native quotient by shifts is four slopes. It can match an extended
# Picard quarter only as a permutation of records, never the joint affine action.
assert all(len(d&set(range(4)))==1 for d in c3)

# Exact principal parts in order (1,-1,i,-i,infinity), with y and x/y
# as local coordinates. Rows are evaluations of three regular differentials.
al=[1,-1,I,-I]
R=s.Matrix([[s.Rational(3,4)/a**2 for a in al]+[0],[0,0,0,0,-3],
            [s.Rational(3,4)/a**3 for a in al]+[0]])
v0=s.Matrix([1,1,1,1,0]);v1=s.Matrix([1,-1,I,-I,0])
V=s.Matrix.hstack(v0,v1)
assert R.rank()==3 and R*V==s.zeros(3,2) and V.H*V==4*s.eye(2)
T=s.zeros(5)
for j,a in enumerate(al):T[j,al.index(I*a)]=1
T[4,4]=-I
assert T*V==V*s.diag(1,I)
Basis=s.Matrix.hstack(v0,v1,*[s.Matrix([a**m for a in al]+[0]) for m in (2,3)],e[:,4])
assert Basis.inv()*T*Basis==s.diag(1,I,-1,-I,-I)
assert V.row_join(s.diag(*[-s.Rational(1,3)]*2,*[s.Rational(1,2)]*2,-s.Rational(1,3))*v0).rank()==3
weak=s.Matrix.hstack(e[:,2],e[:,3])
assert simp(weak.H*V*V.H*weak/4)==s.eye(2)/2
print('Mittag-Leffler principal parts: exact kernel rank two, tau phases 1,i; quotient -1,-i,-i.')
print('The weak coordinate pair has both squared overlaps 1/2 with that kernel; Y does not preserve the kernel.')
print('No scalar phase turns the tau kernel spectrum into the native winding spectrum +i,-i.')
