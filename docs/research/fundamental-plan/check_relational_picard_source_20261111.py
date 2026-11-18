"""Exact original event and pair tests, and Picard actions."""
from pathlib import Path
from collections import Counter, deque
from itertools import combinations, permutations, product
import importlib.util
import numpy as np
import sympy as s

R=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
path=R/'experiments/theory-contracts/compiler-correlated-event-clock-20260919/source_channel.py'
spec=importlib.util.spec_from_file_location('original_source',path)
original=importlib.util.module_from_spec(spec)
spec.loader.exec_module(original)
rays=original.source_rays()
bell,raw,_=original.reflection_actions(rays)
six,_=original.outer_mark_action(raw)
mapping=original.find_bell_partition_intertwiner(bell,six)
checks=0
def check(test,name):
    global checks
    assert test,name
    checks+=1
    print('PASS',name)

def compose(a,b):
    return tuple(a[i] for i in b)

def power(p,n):
    r=tuple(range(len(p)))
    for _ in range(n):r=compose(p,r)
    return r

def order(p):
    i=1
    while power(p,i)!=tuple(range(len(p))):i+=1
    return i

def key(M):
    return np.stack((np.rint(M.real),np.rint(M.imag)),axis=0).astype(np.int8).tobytes()

lookup={key(a):i for i,a in enumerate(bell)}
assert len(lookup)==60
char_results=[]
for mark in range(6):
    others=[i for i in range(6) if i!=mark]
    indices=[i for i,p in enumerate(six) if p[mark]!=mark]
    native_index={i:k for k,i in enumerate(indices)}
    backgrounds=[i for i,p in enumerate(six) if p[mark]==mark]
    generators={}
    for i in backgrounds:
        U=bell[i]
        pe=tuple(native_index[lookup[key(U@bell[j]@U.conj().T)]] for j in indices)
        ps=tuple(others.index(six[i][j]) for j in others)
        generators[(ps,pe)]=U
    base=(tuple(range(5)),tuple(range(20)))
    grp={base}; todo=deque([base]); representatives={base:np.eye(10,dtype=complex)}
    while todo:
        x=todo.popleft()
        for a in generators:
            y=(compose(a[0],x[0]),compose(a[1],x[1]))
            if y not in grp:
                grp.add(y)
                todo.append(y)
                if mark==0:
                    representatives[y]=generators[a]@representatives[x]
    shift=tuple((i+1)%5 for i in range(5))
    gear=tuple(2*i%5 for i in range(5))
    fixed_C=[a[1] for a in grp if a[0]==shift and order(a[1])==5]
    fixed_G=[a[1] for a in grp if a[0]==gear and order(a[1])==4]
    matched=[]
    for C,G in product(fixed_C,fixed_G):
        if compose(G,C)!=compose(power(C,2),G):continue
        elements={compose(power(C,i),power(G,j)) for i,j in product(range(5),range(4))}
        traces=tuple(sorted(sum(i==j for i,j in enumerate(p)) for p in elements))
        matched.append((traces,C,G))
    cases=Counter(tr for tr,_,_ in matched)
    char_results.append((len(grp),len(fixed_C),len(fixed_G),cases))
    print('native mark',mark,'group size, eligible C and G, solution characters',char_results[-1])
    if mark==0:
        regular=[(C,G) for tr,C,G in matched if tr==(0,)*19+(20,)]
        check(len(regular)==16 and len(matched)==16,'sixteen original native affine lifts are all regular')
        C,G=regular[0]
        inverse=lambda x:tuple(x.index(j) for j in range(len(x)))
        conjug=lambda x,y:compose(compose(x,y),inverse(x))
        orbit={(conjug(k[1],C),conjug(k[1],G)) for k in grp if k[0]==base[0]}
        check(orbit==set(regular),'native label kernel acts transitively on the sixteen valid lift pairs')
        cm=representatives[(shift,C)]
        gm=representatives[(gear,G)]
        def scalar(matrix):
            phase=matrix[0,0]
            return phase if np.array_equal(matrix,phase*np.eye(10)) else None
        cphase=scalar(np.linalg.matrix_power(cm,5))
        gphase=scalar(np.linalg.matrix_power(gm,4))
        twist=scalar(gm@cm@gm.conj().T@np.linalg.matrix_power(cm.conj().T,2))
        check(cphase is not None and gphase is not None and twist is not None
              and cphase*twist**5==1,
              'native Bell lifts realize the same affine relation up to removable central phases')
        print('Bell lift powers and relation phases',cphase,gphase,twist)

        # Label the native regular orbit by ordered pairs. Choose base event
        # of complement 3, the midpoint of (0,1), in the native chart.
        native_comp=[others.index(six[j][mark]) for j in indices]
        event0=native_comp.index(3)
        labels_event={}
        for i,j in product(range(5),range(4)):
            el=compose(power(C,i),power(G,j))
            labels_event[el[event0]]=(i,(i+2**j)%5)
        check(len(labels_event)==20 and all(native_comp[i]==3*(a+b)%5
                for i,(a,b) in labels_event.items()),
              'native event complement equals midpoint under the regular ordered-pair labeling')
        orientation={v:k for k,v in labels_event.items()}
        twins=[(i,orientation[(b,a)]) for i,(a,b) in labels_event.items()]
        jumps={}
        for k,j in enumerate(indices):
            U=bell[j]
            stay=np.diag([int(t==u) for t,u in enumerate(original.bell_permutation(U))])
            jumps[k]=U@(np.eye(10)-stay)
        distinguish=[]
        for i,j in twins:
            K,L=jumps[i],jumps[j]
            same_pop=np.array_equal(abs(K)**2,abs(L)**2)
            unequal_word=not np.array_equal(np.kron(K.conj(),K),np.kron(L.conj(),L))
            distinguish.append(same_pop and unequal_word)
        check(all(distinguish),'opposite native orientations agree on populations and differ as coherent jump maps')
        # Exact graph match under the originally pinned marked Bell dictionary.
        _,_,native_A,*_=original.marked_data(mark,bell,six,mapping)
        labels=[]
        for j in range(10):
            side=set(mapping[j])
            if mark not in side:side=set(range(6))-side
            labels.append({others.index(i) for i in side-{mark}})
        A=s.Matrix([[int(e.isdisjoint(f)) for f in labels] for e in labels])
        check(A==s.Matrix(native_A.tolist()),'actual Bell pair permutation intertwines every original adjacency entry')
        C_pair=2*A-s.ones(10)+s.eye(10)
        check(C_pair**2==9*s.eye(10),'original Bell adjacency gives the proposed exact Petersen conference form')
        inc=s.Matrix([[int(i in e) for i in range(5)] for e in labels])
        pm=(s.eye(10)-C_pair/3)/2
        check(inc*(inc.T*inc).inv()*inc.T==pm,'negative projector is exactly the native slot incidence subspace')
        def pair_perm(p):
            return tuple(labels.index({p[i] for i in e}) for e in labels)
        def mat_perm(p):
            return s.Matrix(len(p),len(p),lambda i,j:int(i==p[j]))
        gp=mat_perm(pair_perm(gear)); cp=mat_perm(pair_perm(shift)); P=s.eye(10)-pm
        check(s.trace(gp*P)==-1 and s.trace(gp*pm)==1,
              'native G acts differently on positive five-plane and actual carrier')
        vars=s.Matrix(10,5,s.symbols('t:50'))
        lin=s.linear_eq_to_matrix(list(cp*vars-vars*mat_perm(shift))+list(gp*vars+vars*mat_perm(gear))+list(P*vars-vars),tuple(vars))[0]
        ns=lin.nullspace()
        check(len(ns)==2 and any(s.Matrix(10,5,x).rank()==5 for x in ns),
              'positive plane is the sign twist of the slot representation under full affine group')
        # A pair choice made from native incidence gives an actual P2 mark.
        check(all(sum(inc.row(j))==2 for j in range(10)),'all ten pair labels recover a two-slot mark by incidence')

check(all(x==char_results[0] for x in char_results),'same native affine lift census on all six Bell charts')

# Universal graph selection under a fixed incidence sign rule.
n=s.symbols('n',integer=True)
lam0=(n*n-9*n+14)/2
lam1=7-2*n
lam2=s.Integer(3)
check(s.solve(lam1**2-lam2**2,n)==[2,5] and lam0.subs(n,5)==-3,
      'n>=4 is uniquely five if the same pair-sign operator has scalar square')
# Kernel of forgetting orientation fails the actual pair matrix-unit product.
E=lambda i,j:s.eye(5)[:,i]*s.eye(5)[:,j].T
check((E(0,1)-E(1,0))*E(1,2)==E(0,2),'orientation forgetting does not descend to sequential composition')

# Exceptional Picard curve. Generators from Borner--Bouw--Wewers (4.1.2).
x,y=s.symbols('x y', nonzero=True)
rt=s.sqrt(3)
w=(-1+s.I*rt)/2
simple=lambda M:M.applyfunc(s.simplify)
xp=rt*x/(y+1)
yp=(2-y)/(y+1)
fac0=s.simplify(s.diff(yp,y)*(x/xp)**2)
fac1=s.simplify(s.diff(yp,y)*(x/xp)**3)
check(s.simplify(xp**4-yp**3-1-9*(x**4-y**3-1)/(y+1)**4)==0 and
      fac0==-1 and s.simplify(fac1+(1+y)/rt)==0 and
      s.simplify(yp*fac1-(y-2)/rt)==0,
      'Picard rho preserves the curve and has the stated exact differential pullback')
sig=s.diag(w,s.conjugate(w))
r=s.Matrix([[-1,-2],[-1,1]])/rt
D=s.diag(1,2)
S=s.diag(1,s.sqrt(2))
rho=simple(S*r*S.inv())
tau=s.I*s.eye(2)
check(simple(r.H*D*r)==D and simple(sig.H*sig)==s.eye(2) and
      simple(rho.H*rho)==s.eye(2) and s.simplify(rho.det())==-1 and
      tau**2==-s.eye(2),
      'actual Prym is unitary in the invariant metric, and the quarter is scalar i')

identity4=tuple(range(4))
branch=[s.oo,-s.Integer(1),-w,-w**2]
branch_index=lambda z:next(i for i,u in enumerate(branch) if z==u or
                         (z!=s.oo and u!=s.oo and s.simplify(z-u)==0))
p_sig=tuple(branch_index(s.oo if z==s.oo else w*z) for z in branch)
p_rho=tuple(branch_index(-1 if z==s.oo else s.oo if z==-1 else
                         (2-z)/(z+1)) for z in branch)
p_k=tuple(branch_index(s.conjugate(z)) for z in branch)
def gen_group(gens, n):
    unit=tuple(range(n))
    group={unit}; queue=[unit]
    for a in queue:
        for b in gens:
            c=compose(a,b)
            if c not in group:
                group.add(c)
                queue.append(c)
    return group
check((p_sig,p_rho,p_k)==((0,2,3,1),(1,0,3,2),(0,1,3,2)) and
      len(gen_group([p_sig,p_rho],4))==12 and
      len(gen_group([p_sig,p_rho,p_k],4))==24,
      'actual degree-four branch markers give tetrahedral A4 then S4 with conjugation')
p_x_tau=(1,2,3,0)
p_x_k=(0,3,2,1)
cross=lambda a,b,c,d:s.simplify((a-c)*(b-d)/((a-d)*(b-c)))
lam_x=cross(1,s.I,-1,-s.I)
lam_y=s.simplify((branch[1]-branch[3])/(branch[1]-branch[2]))
check(lam_x==2 and s.simplify(lam_y**2-lam_y+1)==0 and
      len(gen_group([p_x_tau,p_x_k],4))==8,
      'harmonic x quartet and equianharmonic y branch set carry different actions')

pauli=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1)]
def adjoint(u):
    return s.Matrix(3,3,lambda i,j:s.simplify(s.trace(pauli[i]*u*pauli[j]*u.H)/2))
H={identity4:s.eye(3)}
q=[identity4]
K=s.diag(1,-1,1)
for a in q:
    for p,M in [(p_sig,adjoint(sig)),(p_rho,adjoint(rho)),(p_k,K)]:
        z=compose(p,a)
        Mz=simple(M*H[a])
        if z not in H:
            H[z]=Mz
            q.append(z)
        else:
            assert H[z]==Mz
from sympy.combinatorics import Permutation
sgn=lambda p:int(Permutation(list(p)).signature())
check(len(H)==24 and all(M.det()==sgn(p) and
      s.trace(M)==sum(j==k for j,k in enumerate(p))-1 for p,M in H.items()),
      'Prym traceless Hermitian adjoint is ordinary, not signed, S4 augmentation')
import sys
sys.path.insert(0,str(R/'verification'))
from v117_monodromy_weyl_a3 import M0_EXACT as vM, U_EXACT as vU
Hsigned={p:sgn(p)*M for p,M in H.items()}
b=compose(compose(p_sig,p_rho),p_k)
assert order(b)==4
X=s.Matrix(3,3,s.symbols('m:9'))
matches=[]
for a in H:
    if order(a)!=3 or len(gen_group([a,b],4))!=24:continue
    eq=list(vU*X-X*Hsigned[b])+list(vM*X-X*Hsigned[a])
    for v in s.linear_eq_to_matrix(eq,tuple(X))[0].nullspace():
        T=s.Matrix(3,3,v)
        if s.simplify(T.det())!=0:matches.append((a,b))
check(bool(matches) and s.trace(H[b])==-1 and s.trace(vU)==1,
      'exact rank-three intertwiner with native v117 requires the parity twist')
print('Picard signed v117 order-three and order-four marker choices',matches)

# Quaternionification of the actual Prym quarter. tau is central and
# scalar on S; it becomes a relative phase on S plus conjugate(S).
I2=s.eye(2); Z=s.zeros(2); I4=s.eye(4)
block=lambda a,b,c,d:a.row_join(b).col_join(c.row_join(d))
lift=lambda a:s.diag(a,a.conjugate())
anti=block(Z,-I2,I2,Z)
Q=lambda a,b:block(a,b,-b.conjugate(),a.conjugate())
eps=s.Matrix([[0,1],[-1,0]])
G=[Q(a,Z) for a in pauli]+[Q(Z,eps),Q(Z,s.I*eps)]
U=lift(tau)
check(anti*anti.conjugate()==-I4 and all(A.H==A and A*anti==anti*A.conjugate()
      for A in G) and all(G[i]*G[j]+G[j]*G[i]==(2*I4 if i==j else s.zeros(4))
      for i,j in product(range(5),repeat=2)),
      'doubled actual Prym yields five Hermitian quaternionic Clifford directions')
check(G[3]*G[4]==-U and [U*A*U.H for A in G]==G[:3]+[-A for A in G[3:]],
      'same Picard quarter fixes the three Pauli directions and reverses the two normals')
u=s.symbols('u:5',real=True)
V=sum((u[i]*G[i] for i in range(5)),s.zeros(4))
check(simple(V**2)==sum(v**2 for v in u)*I4 and V.trace()==0,
      'same cone has Lorentz quadratic form in five directions and in the fixed three')
projector=(I4+s.I*G[3]*G[4])/2
check(projector==s.diag(1,1,0,0) and
      (s.I*G[3]+G[4])*projector==s.zeros(4) and
      all((projector*A*projector)[:2,:2]==p for A,p in zip(G[:3],pauli)),
      'same quarter fixes chirality of sign wall and gives the three tangential Weyl matrices')
check(96%5!=0 and 48%5!=0,
      'no faithful native order-five element can act by holomorphic or extended Picard automorphisms')
print(f'{checks} passed checks')

# Put the incidence slot representation and the Picard Clifford system
# on the same five-plane by an isometry. The sign choice is substantive.
vol=G[0]*G[1]*G[2]*G[3]*G[4]
check(vol in (I4,-I4) and simple(vol.conjugate())==vol and sgn(gear)==-1,
      'ordinary slot quarter cannot lift by conjugation on one irreducible Clifford module')
normal1=(G[1]-G[4])/s.sqrt(2)
normal2=(G[2]-G[3])/s.sqrt(2)
cliff_transposition=lambda i,j:(G[i]-G[j])/s.sqrt(2)
gspin=simple(cliff_transposition(1,3)*cliff_transposition(1,4)*cliff_transposition(1,2))
cspin=simple(cliff_transposition(0,4)*cliff_transposition(0,3)*
             cliff_transposition(0,2)*cliff_transposition(0,1))
if simple(cspin**5)==-I4:cspin=-cspin
check(all(simple(gspin*G[i]*gspin.H)==-G[gear[i]] and
          simple(cspin*G[i]*cspin.H)==G[shift[i]] for i in range(5)),
      'the sign-twisted full native affine action has explicit Clifford lifts')
check(simple(cspin**5)==I4 and simple(gspin*cspin*gspin.H)==simple(cspin**2) and
      simple(gspin**2) in (simple(normal1*normal2),simple(-normal1*normal2)) and
      simple(gspin**4)==-I4,
      'spin lift has the Picard-type normal-plane square and a nontrivial central order-eight gear')
spin_set={tuple(simple(cspin**i*gspin**j)) for i,j in product(range(5),range(8))}
check(len(spin_set)==40 and
      all(simple(A*anti-anti*A.conjugate())==s.zeros(4) for A in [gspin,cspin]),
      'common Clifford realization has forty quaternionic operators over the twenty native permutations')
z8=(1+s.I)/s.sqrt(2)
rephased=simple(z8*gspin)
check(simple(rephased**4)==I4 and simple(rephased**2)==simple(s.I*gspin**2) and
      simple(rephased*anti-anti*rephased.conjugate())!=s.zeros(4),
      'scalar repair of the Bell fourth power destroys the same quaternionic lift and changes the square')
print(f'{checks} passed checks including the common spin lift')

# The central discrepancy can be invisible on even observables. Test the
# actual ten-dimensional Bell word, not only its unsigned permutation.
duads=list(combinations(range(5),2))
def wedge2_perm(p):
    return np.array([[((1 if (p[a],p[b])==(i,j) else -1) if
                       {p[a],p[b]}=={i,j} else 0)
                      for a,b in duads] for i,j in duads],dtype=complex)
cw,gw=wedge2_perm(shift),wedge2_perm(gear)
pows=lambda M,k:np.linalg.matrix_power(M,k)
sig_bell=[np.trace(pows(cm,i)@pows(gm,j)) for j,i in product(range(4),range(5))]
sig_wedge=[np.trace(pows(cw,i)@pows(gw,j)) for j,i in product(range(4),range(5))]
print('Bell vs exterior square character by j,i:',sig_bell,sig_wedge)
samechar=all(a==b for a,b in zip(sig_bell,sig_wedge))
print('Bell vs exterior square equivalent as complex linear affine representations:',samechar)
# Stronger test: exploit the pinned marked Bell dictionary to solve for a
# diagonal change of basis. A general linear equivalence need not pass it.
ids=[duads.index(tuple(sorted(e))) for e in labels]
cw,gw=cw[np.ix_(ids,ids)],gw[np.ix_(ids,ids)]
d=[None]*10
d[0]=1
queue=[0]; ok=True
for i in queue:
    for Uo,Ue in [(cm,cw),(gm,gw)]:
        j=int(np.argmax(abs(Ue[:,i])))
        if Uo[j,i]==0:
            ok=False
            continue
        v=Uo[j,i]*d[i]/Ue[j,i]
        if d[j] is None:
            d[j]=v
            queue.append(j)
        elif d[j]!=v:
            ok=False
print('Bell vs exterior square through pinned Bell labels and diagonal rephasing:',ok and len(queue)==10)
print('Original Bell lift character class averages:',[(j,sig_bell[5*j:5*j+5]) for j in range(4)])
def pinned_twist(t):
    d=[None]*10; d[0]=1
    queue=[0]
    for i in queue:
        for Uo,Ue in [(cm,cw),(gm,t*gw)]:
            j=int(np.argmax(abs(Ue[:,i])))
            if Uo[j,i]==0:return False
            v=Uo[j,i]*d[i]/Ue[j,i]
            if d[j] is None:
                d[j]=v
                queue.append(j)
            elif d[j]!=v:return False
    return len(queue)==10 and all(abs(z)==1 for z in d)
twists={str(t):pinned_twist(t) for t in (1,-1,1j,-1j)}
print('Actual Bell vs Lambda^2 pinned monomial rephasings by gear character:',twists)
check(not samechar and twists=={'1':False,'-1':False,'1j':True,'(-0-1j)':True},
      'actual Bell pair amplitudes equal the native exterior square with either fourth-root character, not bare')
print(f'{checks} passed checks including the original exterior-square amplitudes')

# Common normalization over x: y_E=t^3, y_P=t^2, t=y_E/y_P.
# At a finite branch point ord(t)=1 and ord(dx)=5; at each infinity
# ord(x)=-3, ord(t)=-2, ord(dx)=-4. The monomials below exhaust genus 7.
forms=[(i,j) for j in range(1,6) for i in range(3) if 2*j-3*i-4>=0]
genus=1+(-2*6+4*(6-1)+(6-2))//2
e_forms=[(i,j) for i,j in forms if j%3==0]
p_forms=[(i,j) for i,j in forms if j%2==0]
new_forms=[f for f in forms if f not in e_forms+p_forms]
check(genus==len(forms)==7 and (len(e_forms),len(p_forms),len(new_forms))==(1,3,3),
      'common normalized sixth cover has genus seven and separates elliptic, Picard and primitive forms')
def phase_char(i,j):
    # Elliptic deck fixing Picard, and shared x -> ix.
    return ((-1)**j,s.I**(i+1))
check([phase_char(*f) for f in p_forms]==[(1,s.I),(1,s.I),(1,-1)] and
      [phase_char(*f) for f in e_forms]==[(-1,s.I)],
      'shared x-quarter and elliptic deck act diagonally on actual pullback differentials')
De=s.Matrix([[0,-1],[1,0]])
Qe=s.diag(1,-1)
check(De**2==-I2 and Qe*De*Qe.inv()*De.inv()==-I2 and
      all(s.kronecker_product(Qe,s.diag(*[phase_char(*f)[1] for f in fs]))*
          s.kronecker_product(De,s.diag(*[phase_char(*f)[0] for f in fs])) ==
          -s.kronecker_product(De,s.diag(*[phase_char(*f)[0] for f in fs]))*
          s.kronecker_product(Qe,s.diag(*[phase_char(*f)[1] for f in fs]))
          for fs in [p_forms,new_forms]),
      'neither Picard nor primitive pullback differentials cancel the old elliptic normal spin commutator')
print('Common sixth-cover differential slots (i,j) and phases:',
      [(f,phase_char(*f)) for f in forms])
print(f'{checks} passed checks including the older elliptic source bridge')
