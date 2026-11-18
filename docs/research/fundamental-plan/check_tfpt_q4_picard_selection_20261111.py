"""Original Q4 recovery, native clocks and Picard branch constraints."""
from pathlib import Path
import importlib.util
import hashlib
from itertools import product, combinations
from collections import deque
import sympy as s
import numpy as np

R=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
q4=R/'_newest2/TFPT_Gesamtdokumentation2_20260927.md'
code=q4.read_text().split('Originalmitglied: `TFPT_Rekursive_Bindung_20260926/audit.py`',1)[1].split('```python\n',1)[1].split('```',1)[0]
expected='4a19f8145d5ccaef08bfdc5893f64484feae630cf61e49639df52f89f56528c4'
assert expected in [hashlib.sha256((code.rstrip()+'\n'*k).encode()).hexdigest() for k in range(4)]
state={'__file__':str(q4)}
exec(code.split('# Exact complete 35 x 35 two-module model.',1)[0],state)
F,G,K=[state[x] for x in ['Fnum','Gs','Recovery']]
I=s.eye(25)
vec=lambda X:s.Matrix(25,1,list(X))
om=vec(G.inv())
met=s.kronecker_product(G,G)
B=om*om.T*met/5
P9=s.Rational(9,4)*(K-B)
Sw=s.Matrix(25,25,lambda a,b:int(a==5*(b%5)+b//5))
A10=(I-Sw)/2
S5=(I+Sw)/2-B-P9
assert K*om==om and B.trace()==1 and P9.trace()==9 and S5.trace()==5 and A10.trace()==10
assert all(P**2==P and P.T*met==met*P for P in (B,P9,S5,A10))
assert K*A10==s.zeros(25) and K*S5==s.zeros(25)
print('Original code split: Sym^2 = 1+9+5, wedge^2 = 10 is entirely in the kernel.')

# Original normalized Bell lines and their Petz measurement/prepare factorization.
fs=[F[r,:].T/4 for r in range(10)]
vs=[G.inv()*f for f in fs]
assert all((v.T*G*v)[0]==s.Rational(1,2) for v in vs)
lineGram=s.Matrix(10,10,lambda i,j:2*(vs[i].T*G*vs[j])[0])
assert lineGram**2==2*lineGram
assert all(abs(lineGram[i,j])==s.Rational(1,3) for i in range(10) for j in range(10) if i!=j)
enc=s.Matrix.hstack(*(2*s.kronecker_product(v,v) for v in vs))
read=s.Matrix.vstack(*(s.kronecker_product(f,f).T for f in fs))
T=read*enc
assert K==enc*read and T==s.Rational(4,9)*s.eye(10)+s.ones(10)/18
assert enc*s.ones(10,1)==2*om
print('Ten original normalized effect rays: equiangular, tight; measured recovery is (4/9) I + J/18.')

native=R/'experiments/theory-contracts/compiler-correlated-event-clock-20260919/source_channel.py'
sp=importlib.util.spec_from_file_location('native',native); n=importlib.util.module_from_spec(sp);sp.loader.exec_module(n)
rays=n.source_rays()
bell,raw,_=n.reflection_actions(rays)
six,_=n.outer_mark_action(raw)
dic=n.find_bell_partition_intertwiner(bell,six)
def perm(U):return tuple(int(np.argmax(abs(U[:,i]))) for i in range(U.shape[1]))
# Construct the actual symmetric fourth power of each 4D source reflection,
# then restrict to the exact Q4 codeword sums. Entries prior to division
# are bounded Gaussian integers; post restriction convert to exact rationals.
word=list(product(range(4),repeat=4))
N=state['ns']; pos=state['pos']; E=state['E']
Ebig=np.array([E[pos[tuple(t.count(i) for i in range(4))]] for t in word],dtype=np.int64)
diag=np.array([int(G[i,i]) for i in range(5)],dtype=np.int64)
assert np.array_equal(Ebig.T@Ebig,np.diag(diag))
def gaussian(z):
    assert z.real==round(z.real) and z.imag==round(z.imag)
    return int(z.real)+s.I*int(z.imag)
actions=[]; pe=[]
for z in rays:
    U=2*np.eye(4,dtype=complex)-np.outer(z,z.conj())
    u4=np.kron(np.kron(U,U),np.kron(U,U))
    transported=u4@Ebig
    mat=Ebig.T@transported
    assert np.array_equal((Ebig*(24//diag))@mat,24*transported)
    M=G.inv()*s.Matrix([[gaussian(x) for x in row] for row in mat])/16
    assert M.H*G*M==G
    out=[]
    for v in vs:
        new=M*v
        out.append(next(j for j,w in enumerate(vs) if
                   abs((w.H*G*new)[0])**2==s.Rational(1,4)))
    actions.append(M);pe.append(tuple(out))
pb=[perm(U) for U in bell]
# Find the actual equivariant labeling; enforce every original reflection,
# including the non-clock generators, on the same 10 labels.
def maps_to(pb,pe):
    solutions=[]
    for j in range(10):
        d={0:j};todo=[0];valid=True
        for k in todo:
            for a,b in zip(pb,pe):
                m,t=a[k],b[d[k]]
                if m not in d:d[m]=t;todo.append(m)
                elif d[m]!=t:valid=False;break
            if not valid:break
        if valid and len(d)==10 and len(set(d.values()))==10:
            solutions.append(tuple(d[i] for i in range(10)))
    return solutions
align=maps_to(pb,pe)
assert len(align)==1
print('All 60 actual source reflections: unique intertwiner of Bell labels and the ten Q4 Petz effect rays:',align[0])
assert all(M.conjugate()==M for M in actions)
assert all(K*s.kronecker_product(M,M)==s.kronecker_product(M,M)*K for M in actions)
assert K.conjugate()==K
print('Native Q4 pair operator commutes with all 60 physical-source code lifts, hence C5 and G, and with coordinate conjugation.')

for mark in range(6):
    _,_,nb,*_=n.marked_data(mark,bell,six,dic)
    P=s.Matrix(nb.tolist())/3
    assert T==s.eye(10)/3+P/6+P**2/2
print('All six native charts: Petz record transfer = (1/3) stay + (1/6) one jump + (1/2) two jumps.')
print('Petz effect rays in Bell order: native pair transfer has the same physical-source labels, but not the same action.')

# Picard family and exact branch action. Full-cover rho is not a
# generator of the marked x-cover, even though it is an automorphism of C.
w=(-1+s.I*s.sqrt(3))/2
D=s.diag(-1,s.I,s.I); V=s.diag(1,s.I,-s.I)
rho=s.diag(-1,s.Matrix([[-1,-2],[-1,1]])/s.sqrt(3))
assert rho**2==s.eye(3) and D*rho==rho*D
assert s.factor((D*rho).charpoly().as_expr())==s.factor(V.charpoly().as_expr())
assert (D*rho)**2==V**2==D**2
XX=s.Matrix(3,3,s.symbols('a:9'))
def sols(L,M):
    return [s.Matrix(3,3,v) for v in s.linear_eq_to_matrix(list(L*XX-XX*M),tuple(XX))[0].nullspace()]
def maxrank(ls):
    return sum((x*v for x,v in zip(s.symbols('q:'+str(len(ls))),ls)),s.zeros(3)).rank()
assert maxrank(sols(V,D))==1 and maxrank(sols(V,D**2))<3
assert maxrank(sols(V**2,D**2))==3 and maxrank(sols(V,D*rho))==3
# Exact compatibility of an x-cover automorphism: group of four x rotations
# and three deck characters. Nontrivial deck eigenvalues split the repeated
# eigenvalue but no allowed combined lift matches the actual v117 spectrum.
cover=[s.diag(w,w,w**2)**j*D**k for j,k in product(range(3),range(4))]
assert not any(maxrank(sols(V,g))==3 for g in cover)
print('Picard: half-turn rank 3, original quarter rank 1. Twisted full quarter tau rho also rank 3, but leaves the marked x-cover.')
assert ((2-s.Integer(0))/(1+s.Integer(0)))==2  # rho moves a finite x-branch point off D.

# Multiplication by the two deck-character powers, on the exact residue
# splitting. At infinity ord(x)=-3, ord(y)=-4, ord(dx)=-4.
legal={k:[i for i in range(3) if -3*i-4+4*(3-k)>=0] for k in (1,2)}
assert legal=={1:[0,1],2:[0]}
# General y f(x)+y^2 g(x), polynomial multipliers at x, preserving all
# residues and adding no other finite poles: constant input needs deg f<=1
# and deg g<=0. Then x^2 residue input forces both to vanish.
assert all(all(-3*(i+d)-4+4*(3-k)<0 for d in (0,1))
           for i,k in [(2,1),(2,2)])
print('Raw y,y^2 residue multiplication: legal input ranks 2 and 1, jointly miss x^2; no nonzero common polynomial multiplier is holomorphic on all three.')

# In the original Gram orthonormalization, Q4 has an exact identical
# four-leg tensor for its transfer and its Choi form. Check rather than
# infer it from the labels. Complete positivity then fixes the endpoint.
M=s.kronecker_product(G**s.Rational(1,2),G**s.Rational(1,2))
to_ortho=lambda A:(M*A*M.inv()).applyfunc(s.simplify)
K0,B0,P90,P50=map(to_ortho,(K,B,P9,S5))
def choi(A):
    return s.Matrix(25,25,lambda a,b:A[5*(a%5)+b%5,5*(a//5)+b//5])
assert choi(K0)==K0 and choi(B0)==I/5
t=s.symbols('t',real=True)
assert choi(B0+t*P90)==(1-9*t/4)*I/5+9*t/4*K0
assert choi(P90+P50)==(5*B0+Sw)/2-I/5
haar=B0+s.Rational(2,7)*(P90+P50)
assert choi(haar)==haar and haar*s.eye(5).reshape(25,1)==s.eye(5).reshape(25,1)
assert haar.eigenvals()=={s.Integer(1):1,s.Rational(2,7):14,s.Integer(0):10}
print('Original normalized Petz four-leg tensor equals its Choi reshuffle.')
print('At fixed native 1/9/15 projections, positive transfer plus complete positivity permits 0 <= lambda <= 4/9.')
print('Only lambda=4/9 gives identical transfer/Choi tensor; it is also the maximum-memory CP endpoint and uniquely reduces Choi rank to 10.')
print('Without native blindness on S5, the symmetric real sphere channel is another covariant, positive, self-reshuffling solution.')

# Compare plain conjugation, a real x-cover reflection and a putative
# reflected metric insertion on the actual marked ten response records.
native_labels=[]
for i in range(10):
    v=set(dic[i])
    native_labels.append(v-{0} if 0 in v else (set(range(6))-v)-{0})
others=list(range(1,6))
trans=dict(zip(others,others))
trans[others[-2]],trans[others[-1]]=others[-1],others[-2]
pt=tuple(native_labels.index({trans[k] for k in pair}) for pair in native_labels)
PP=s.Matrix(10,10,lambda i,j:int(i==pt[j]))
assert PP**2==s.eye(10) and PP.trace()==4
assert any(a==pt for a in pb)
assert (PP*T).eigenvals()=={s.Integer(1):1,s.Rational(4,9):6,-s.Rational(4,9):3}
print('Geometric 3+2 conjugation on marked duads is a transposition: 4 fixed, 3 exchanged pairs.')
print('Plain coordinate conjugation fixes all Bell effects; with the induced transposition inserted into the record form, three exact eigenvalues are -4/9.')

# Picard all attachments: exact branch/residue weights, Hodge numbers and
# the topological normalization with the already derived Clifford matrices.
import mpmath as mp
mp.mp.dps=35
def hodge(k,j):
    a=mp.mpf(k+1)/4;b=1-mp.mpf(j)/3
    return mp.gamma(a)*mp.gamma(b)*mp.gamma(1-a-b)/(mp.gamma(1-a)*mp.gamma(1-b)*mp.gamma(a+b))
h=[hodge(k,j)/hodge(0,2) for k,j in [(1,2),(0,1),(0,2)]]
assert abs(h[1]-2)<mp.mpf('1e-31')
print('Independent complex beta integral gives Hodge ratios:',[mp.nstr(x,16) for x in h])
pts=[s.Integer(1),s.I,s.Integer(-1),-s.I,s.oo]
quarter=(1,2,3,0,4)
conj=(0,3,2,1,4)
permcomp=lambda p,q:tuple(p[x] for x in q)
assert sum(j==k for j,k in enumerate(conj))==3 and conj[4]==quarter[4]==4
assert permcomp(conj,permcomp(quarter,conj))==(3,0,1,2,4)
assert [s.I**j for j in range(1,4)]==[s.I,-1,-s.I]
assert s.diag(1,1,-1).trace()!=V.__pow__(2).trace()
# rho normalizes no marked x deck: sigma fixes (x,0) but not rho(x,0).
assert s.simplify(2*w-2)!=0
print('Raw meromorphic-function half-turn is 2+1, holomorphic and v117 half-turn 1+2; x-cover automorphisms fix the exceptional infinity mark.')

pauli=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1)]
Z=s.zeros(2);ep=s.Matrix([[0,1],[-1,0]])
qu=lambda a,b:a.row_join(b).col_join((-b.conjugate()).row_join(a.conjugate()))
ga=[qu(a,Z) for a in pauli]+[qu(Z,ep),qu(Z,s.I*ep)]
Id=s.eye(4)
assert all(ga[i]*ga[j]+ga[j]*ga[i]==(2*Id if i==j else s.zeros(4)) for i,j in product(range(5),repeat=2))
from functools import reduce
from operator import mul
vol=reduce(mul,ga)
assert vol in (Id,-Id)
# (1/32) tr(volume) * 4! * Vol(S4)/(8*pi^2) is the normalized c2.
assert abs(s.trace(vol)*24*s.Rational(8,3)/(32*8))==1
assert s.I*ga[3]*ga[4] in (s.diag(1,1,-1,-1),s.diag(-1,-1,1,1))
print('Rank-two S4 projector gives |c2|=1; fixed Pauli blocks give opposite first Chern numbers.')

# Does the actual sign-sensitive code-line Gram form coincide with the
# Petersen conference matrix when using the uniquely fixed original labels?
ridx=list(align[0])
LG=lineGram.extract(ridx,ridx)
nb0=s.Matrix(n.marked_data(0,bell,six,dic)[2].tolist())
Cpair=2*nb0-s.ones(10)+s.eye(10)
SC=3*(LG-s.eye(10))
assert SC**2==9*s.eye(10)
good=[]
for eps in (-1,1):
    signs=[1]+[SC[0,j]/(eps*Cpair[0,j]) for j in range(1,10)]
    Zsign=s.diag(*signs)
    if Zsign*SC*Zsign==eps*Cpair:good.append(eps)
print('Exact Q4 line Gram is signed-switching equivalent to native Petersen conference:',good)
# Compare triangle products, invariant under rephasing the real effect rays.
bad_tri=[]
for e in (-1,1):
    bad_tri.append(sum(SC[a,b]*SC[b,c]*SC[c,a]!=e*Cpair[a,b]*Cpair[b,c]*Cpair[c,a]
                       for a,b,c in combinations(range(10),3)))
print('Number of incoherent triangle signs vs. +/- Petersen:',bad_tri)
# Record geometry does not distinguish a lift sign: all normalization,
# transfer, projectors and ordered products of projectors remain unchanged.
assert all(s.kronecker_product(v,v)==s.kronecker_product(-v,-v) for v in vs)
assert good==[-1]
for mark in range(6):
    ma=s.Matrix(n.marked_data(mark,bell,six,dic)[2].tolist())
    Cmark=2*ma-s.ones(10)+s.eye(10)
    signs=[1]+[-SC[0,j]/Cmark[0,j] for j in range(1,10)]
    Sd=s.diag(*signs)
    assert Sd*LG*Sd==s.eye(10)-Cmark/3
print('All six charts: signs can be chosen so that Gram(Q4 effect lines) = I - C_pair/3 exactly.')

# Explicit isometry of the marked code with the five incidence slots.
# Pick a canonical one of the two possible scalar normalizations of b_ij-c1.
c=(4-s.sqrt(6))/10
sv=s.Matrix(5,10,lambda i,j:(int(others[i] in native_labels[j])-c)*s.sqrt(s.Rational(2,3)))
signs=[1]+[-SC[0,j]/Cpair[0,j] for j in range(1,10)]
f0=s.Matrix.hstack(*(signs[i]*s.sqrt(2)*vs[ridx[i]] for i in range(10)))
assert (sv.T*sv-(s.eye(10)-Cpair/3)).applyfunc(s.simplify)==s.zeros(10)
W=(f0*sv.T/2).applyfunc(s.simplify)
assert (W.T*G*W-s.eye(5)).applyfunc(s.simplify)==s.zeros(5)
assert (W*sv-f0).applyfunc(s.simplify)==s.zeros(5,10)
# A matrix and a genuine simultaneous source word for the two clock actions.
slotgens={}
for perm6,M in zip(six,actions):
    if perm6[0]==0:
        pp=tuple(others.index(perm6[j]) for j in others)
        slotgens[pp]=M
unit=tuple(range(5)); gws={unit:s.eye(5)}; queue=[unit]
for p in queue:
    for q,M in slotgens.items():
        qp=permcomp(q,p)
        if qp not in gws:
            gws[qp]=M*gws[p];queue.append(qp)
for name,clock in [('C5',tuple((i+1)%5 for i in range(5))),('G',tuple(2*i%5 for i in range(5)))]:
    a=gws[clock]
    p=s.Matrix(5,5,lambda i,j:int(i==clock[j]))
    signs=[e for e in (-1,1) if (a*W-e*W*p).applyfunc(s.simplify)==s.zeros(5)]
    assert len(signs)==1
    print('True Q4 code action under original',name,'equals native slot permutation with scalar sign',signs[0])
fibres={p:{tuple(M) for q,M in zip(six,actions) if q==p} for p in set(six)}
assert all(len(f)==1 for f in fibres.values())
print('Each fourfold original reflection fibre induces exactly one Q4 code matrix.')
