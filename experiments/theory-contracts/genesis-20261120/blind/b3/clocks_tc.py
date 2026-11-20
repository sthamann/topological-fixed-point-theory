"""H2 test: k cyclic clocks x_1..x_k of orders p_1..p_k, single law x_1 x_2 ... x_k = 1 (joint return).
Group = <x_1..x_{k-1} | x_i^{p_i}, (x_1...x_{k-1})^{p_k}>.  Todd-Coxeter (HLT) over the trivial subgroup
decides finiteness when it terminates; 'CAP' = >CAP cosets defined (enumeration cannot prove infinity;
the theorem 1/p+1/q+1/r<=1 => infinite is quoted, the TC result is consistent evidence)."""
import itertools, math, sys, time
import numpy as np
CAP=30000
class Overflow(Exception): pass
def tc(ngen,rels,cap=CAP):
    A=2*ngen; inv=lambda s:s^1
    tab=[[None]*A]; par=[0]
    def find(c):
        while par[c]!=c: par[c]=par[par[c]]; c=par[c]
        return c
    def define(c,s):
        if len(tab)>cap: raise Overflow
        d=len(tab); tab.append([None]*A); par.append(d); tab[c][s]=d; tab[d][inv(s)]=c
    def coincide(a,b):
        q=[(a,b)]
        while q:
            a,b=q.pop(); a,b=find(a),find(b)
            if a==b: continue
            if a>b: a,b=b,a
            par[b]=a
            for s in range(A):
                d=tab[b][s]
                if d is None: continue
                tab[b][s]=None
                d=find(d)
                if tab[d][inv(s)] is not None and find(tab[d][inv(s)])==b: tab[d][inv(s)]=None
                e=tab[a][s]
                if e is None:
                    tab[a][s]=d
                    if tab[d][inv(s)] is None: tab[d][inv(s)]=a
                    else: q.append((find(tab[d][inv(s)]),a))
                else: q.append((find(e),d))
    def scan(c,w):
        f=c;i=0;b=c;j=len(w)-1
        while True:
            while i<=j and tab[f][w[i]] is not None: f=find(tab[f][w[i]]); i+=1
            if i>j:
                if f!=b: coincide(f,b)
                return
            while j>=i and tab[b][inv(w[j])] is not None: b=find(tab[b][inv(w[j])]); j-=1
            if j<i: coincide(f,b); return
            if j==i:
                tab[f][w[i]]=b; tab[b][inv(w[i])]=f; return
            define(f,w[i])
    c=0
    try:
        while c<len(tab):
            if find(c)==c:
                for w in rels:
                    if find(c)!=c: break
                    scan(c,w)
                if find(c)==c:
                    for s in range(A):
                        if tab[c][s] is None: define(c,s)
            c+=1
    except Overflow: return None
    return sum(1 for x in range(len(tab)) if find(x)==x)
def clock_group(orders):
    k=len(orders); ng=k-1
    rels=[[2*i]*orders[i] for i in range(ng)]+[[2*i for i in range(ng)]*orders[-1]]
    return tc(ng,rels)
cop=lambda t: all(math.gcd(a,b)==1 for a,b in itertools.combinations(t,2))
t0=time.time()
print("== sanity: Z5 =",tc(1,[[0]*5]),"; S3=<a,b|a^2,b^3,(ab)^2> =",clock_group((2,3,2)),"; <a,b|a^2,b^3,(ab)^7> (infinite) =",clock_group((2,3,7)))
print("\n== A. k=3 clocks, all 2<=p<=q<=r<=9 (120 triples)")
fin=[];mism=0
for t in itertools.combinations_with_replacement(range(2,10),3):
    o=clock_group(t); s=sum(1/x for x in t)
    pred=round(2/(s-1)) if s>1+1e-12 else None
    if o!=pred: mism+=1; print("MISMATCH with spherical criterion",t,o,pred)
    if o is not None: fin.append((t,o,cop(t)))
print("finite groups found (triple, order, pairwise coprime):")
for f in fin: print("  ",f)
print("mismatches with 1/p+1/q+1/r>1 criterion:",mism)
print("pairwise-coprime triples tested:",sum(1 for t in itertools.combinations_with_replacement(range(2,10),3) if cop(t)),
      " finite & nontrivial among them:",[f for f in fin if f[2]])
print("\n== B. NEGATIVE CONTROLS on the number of clocks")
print("k=2 (x^p=1=x^q), coprime pairs p<q<=9: orders",sorted({tc(1,[[0]*p,[0]*q]) for p in range(2,10) for q in range(p+1,10) if math.gcd(p,q)==1}))
for t in [(2,3,5,7),(2,3,5,11),(2,3,7,11),(2,5,7,9),(3,4,5,7),(2,3,5,7,11)]:
    print(" k=%d coprime"%len(t),t,"->",clock_group(t) or "CAP exceeded (infinite by Fuchsian criterion)")
print(" k=4 non-coprime (2,2,2,2) ->",clock_group((2,2,2,2)) or "CAP exceeded (Euclidean, infinite)")
print("TC time %.1fs"%(time.time()-t0))

print("\n== C. the forced group: structure of <2,3,5> and its double cover inside the unit quaternions")
def qm(a,b):
    a0,a1,a2,a3=a;b0,b1,b2,b3=b
    return (a0*b0-a1*b1-a2*b2-a3*b3,a0*b1+a1*b0+a2*b3-a3*b2,a0*b2-a1*b3+a2*b0+a3*b1,a0*b3+a1*b2-a2*b1+a3*b0)
key=lambda q: tuple(round(x,6)+0.0 for x in q)
phi=(1+5**.5)/2
gens=[key((.5,.5,.5,.5)),key((phi/2,1/(2*phi),.5,0))]
def closure(gens,mul):
    S=set(gens);fr=list(gens)
    while fr:
        nf=[]
        for a in fr:
            for g in gens:
                c=mul(a,g)
                if c not in S:S.add(c);nf.append(c)
        fr=nf
    return sorted(S)
m2=lambda a,b:key(qm(a,b))
def canon(q):
    for x in q:
        if abs(x)>1e-9: return q if x>0 else key(tuple(-y for y in q))
m1=lambda a,b:canon(key(qm(a,b)))
def analyse(name,el,mul,chiV=None):
    n=len(el); idx={g:i for i,g in enumerate(el)}; e=[g for g in el if all(mul(g,h)==h for h in el[:5])][0]
    T=np.array([[idx[mul(a,b)] for b in el] for a in el]); invs=np.array([np.where(T[i]==idx[e])[0][0] for i in range(n)])
    cls=-np.ones(n,int);reps=[]
    for i in range(n):
        if cls[i]<0:
            for h in range(n): cls[T[T[h,i],invs[h]]]=len(reps)
            reps.append(i)
    r=len(reps); size=np.bincount(cls)
    comm={T[T[a,b],T[invs[a],invs[b]]] for a in range(n) for b in range(n)}
    der=closure([el[i] for i in comm],mul)
    # Burnside-Dixon: class constants a_ijk, common eigenvectors
    a=np.zeros((r,r,r))
    for x in range(n):
        for y in range(n): a[cls[x],cls[y],cls[T[x,y]]]+=1
    a/=size[None,None,:]
    rng=np.random.default_rng(1); M=sum(rng.normal()*a[i] for i in range(r))
    w,V=np.linalg.eig(M); V=V/V[cls[idx[e]]][None,:]
    chars=[]
    for c in range(r):
        om=V[:,c]; d=np.sqrt(n/np.sum(np.abs(om)**2/size)).real; chars.append(d*om/size)
    chars=sorted(chars,key=lambda c:(round(c[cls[idx[e]]].real),tuple(np.round(c.real,3))))
    deg=[int(round(c[cls[idx[e]]].real)) for c in chars]
    sq=np.array([cls[T[i,i]] for i in range(n)])
    fs=[int(round(np.sum(c[sq]).real/n)) for c in chars]
    print(f"{name}: order={n} classes={r} |[G,G]|={len(der)} (perfect={len(der)==n}) irrep degrees={deg} sum d^2={sum(d*d for d in deg)} Frobenius-Schur={fs}")
    if chiV:
        cv=np.array([chiV(el[reps[k]]) for k in range(r)])
        A=np.array([[np.sum(size*cv*ci*np.conj(cj)).real/n for cj in chars] for ci in chars]); A=np.round(A).astype(int)
        print(" McKay graph (tensor with the 2-dim quaternion rep): adjacency\n",A)
        print(" symmetric:",(A==A.T).all()," edges:",A.sum()//2," node degrees:",sorted(A.sum(1))," (2-A)@dims =",(2*np.eye(r,dtype=int)-A)@np.array(deg))
        b=[i for i in range(r) if A[i].sum()==3][0]; arms=[]
        for s in np.where(A[b])[0]:
            L=1;prev=b;cur=s
            while True:
                nx=[j for j in np.where(A[cur])[0] if j!=prev]
                if not nx: break
                prev,cur=cur,nx[0];L+=1
            arms.append(L)
        print(" tree with one branch node, arm lengths",sorted(arms),"-> affine E8 iff [1,2,5]")
    return el,T,idx,e
G2=closure(gens,m2); el,T,idx,e=analyse("binary group 2I = <a,b>",G2,m2,chiV=lambda q:2*q[0])
G1=closure([canon(g) for g in gens],m1); el1,T1,idx1,e1=analyse("quotient 2I/{+-1}",G1,m1)
def order(i,T,ie):
    k=1;x=i
    while x!=ie: x=T[x,i];k+=1
    return k
ie=idx1[e1]; ords=[order(i,T1,ie) for i in range(len(el1))]
print("element orders in quotient:",{o:ords.count(o) for o in sorted(set(ords))})
ok=[(a,b) for a in range(60) for b in range(60) if ords[a]==2 and ords[b]==3 and ords[T1[a,b]]==5]
a,b=ok[0]; gen=closure([el1[a],el1[b]],m1)
print("pairs (a,b) with orders (2,3) and ab of order 5:",len(ok),"; such a pair generates",len(gen),"elements -> quotient realises the (2,3,5) clock group of order 60")
