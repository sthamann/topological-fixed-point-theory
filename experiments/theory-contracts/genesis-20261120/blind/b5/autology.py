"""H1 Autology: structure (V,B) with  Lie(Aut(V,B)) ~= V  as Aut-modules.
Part A: exhaustive over single bilinear forms B (not assumed symmetric / nondegenerate).
Part B: named classical structures (controls).
Part C: 'smallest nontrivial irrep = adjoint' over all simple Lie algebras rank<=8."""
import numpy as np, itertools, collections, time
from fractions import Fraction
TOL=1e-8
def nullspace(M):
    if M.shape[0]==0: return np.eye(M.shape[1])
    u,s,vt=np.linalg.svd(M)
    r=int((s>TOL).sum()); return vt[r:].T
def comm(n):
    K=np.zeros((n*n,n*n))
    for i in range(n):
        for j in range(n): K[i*n+j,j*n+i]=1   # row-major vec: vec(X^T)=K vec(X)
    return K
def stab_ops(n):
    I=np.eye(n);K=comm(n)
    # row-major vec: vec(AXC)=(A kron C^T) vec X
    return I,K
def stab_matrix(B,n,I,K):
    # X^T B + B X = 0 ; vec(BX)=(B kron I)vecX ; vec(X^T B)=(I kron B^T) K vecX
    return np.kron(B,I)+np.kron(I,B.T)@K
def module_iso(Xs,n):
    """Xs: list of d matrices spanning g (d must equal n). Returns (dimHom, iso?)"""
    d=len(Xs)
    T=np.array([X.reshape(-1) for X in Xs])           # d x n^2
    br=np.array([[ (Xs[a]@Xs[i]-Xs[i]@Xs[a]).reshape(-1) for i in range(d)] for a in range(d)]) # a,i,n^2
    rows=[]
    # unknown P (d x n): phi(e_j)=sum_i P[i,j] X_i
    cols=[]
    for r in range(d):
        for s in range(n):
            P=np.zeros((d,n));P[r,s]=1
            res=[]
            for a in range(d):
                L=(P@Xs[a]).T@T          # n x n^2 : row j = phi(X_a e_j)
                R=P.T@br[a]              # n x n^2 : row j = [X_a,phi(e_j)]
                res.append((L-R).reshape(-1))
            cols.append(np.concatenate(res))
    A=np.array(cols).T
    N=nullspace(A)
    if N.shape[1]==0: return 0,False
    rng=np.random.default_rng(1)
    best=0
    for _ in range(6):
        P=(N@rng.standard_normal(N.shape[1])).reshape(d,n)
        best=max(best,abs(np.linalg.det(P)))
    return N.shape[1],best>1e-6
def killing_sig(Xs):
    d=len(Xs);T=np.array([X.reshape(-1) for X in Xs])
    pinv=np.linalg.pinv(T.T)
    ad=[]
    for a in range(d):
        M=np.array([pinv@((Xs[a]@Xs[i]-Xs[i]@Xs[a]).reshape(-1)) for i in range(d)]).T
        ad.append(M)
    Kf=np.array([[np.trace(ad[a]@ad[b]) for b in range(d)] for a in range(d)])
    ev=np.linalg.eigvalsh((Kf+Kf.T)/2)
    return (int((ev<-1e-7).sum()),int((abs(ev)<=1e-7).sum()),int((ev>1e-7).sum()))
def describe(B):
    n=B.shape[0];S=(B+B.T)/2;A=(B-B.T)/2
    rk=np.linalg.matrix_rank(B)
    sym=np.allclose(A,0);anti=np.allclose(S,0)
    ev=np.linalg.eigvalsh(S)
    sig=(int((ev>1e-9).sum()),int((ev<-1e-9).sum()))
    return rk,('sym' if sym else 'antisym' if anti else 'mixed'),sig
def scan(n,entries,sample=None,seed=0):
    I,K=stab_ops(n)
    stats=collections.Counter();hits=collections.Counter();dimmatch=0;total=0
    if sample is None: it=itertools.product(entries,repeat=n*n)
    else:
        rng=np.random.default_rng(seed)
        it=(tuple(rng.choice(entries,size=n*n)) for _ in range(sample))
    for t in it:
        B=np.array(t,float).reshape(n,n)
        if not B.any(): continue
        total+=1
        M=stab_matrix(B,n,I,K)
        s=np.linalg.svd(M,compute_uv=False)
        d=n*n-int((s>TOL).sum())
        if d!=n: continue
        dimmatch+=1
        N=nullspace(M);Xs=[N[:,i].reshape(n,n) for i in range(n)]
        h,iso=module_iso(Xs,n)
        rk,ty,sig=describe(B)
        key=(rk,ty,sig,killing_sig(Xs))
        stats[key]+=1
        if iso: hits[key]+=1
    return total,dimmatch,stats,hits
print("=== PART A: single bilinear form B on R^n, g={X: X^T B + B X=0}; test dim g=n and V~=g as g-module ===")
for n,entries,sample in [(1,(-1,0,1),None),(2,(-1,0,1),None),(3,(-1,0,1),None),(4,(0,1),None),(4,(-1,0,1),60000),(5,(-1,0,1),20000)]:
    t0=time.time()
    total,dm,stats,hits=scan(n,entries,sample)
    print(f"\nn={n} entries={entries} {'EXHAUSTIVE' if sample is None else 'RANDOM SAMPLE '+str(sample)}: forms={total}, dim g = n for {dm}, time {time.time()-t0:.1f}s")
    print("  class key = (rank B, symmetry type, signature of sym part (+,-), Killing signature (-,0,+)) : #dim-match / #module-iso")
    for k in sorted(stats): print("   ",k,":",stats[k],"/",hits.get(k,0))
# ---------- Part B
print("\n=== PART B: named structures (controls). cols: structure, dim V, dim g, Hom_g(V,g) dim, V~=g ?, Killing sig (-,0,+) ===")
def stab_general(n,constraints):
    """constraints: list of matrices C (m x n^2) such that C vecX=0"""
    M=np.vstack(constraints);N=nullspace(M)
    return [N[:,i].reshape(n,n) for i in range(N.shape[1])]
def report(name,n,Xs):
    d=len(Xs)
    if d==n:
        h,iso=module_iso(Xs,n)
    else: h,iso='-',False
    ks=killing_sig(Xs) if d>0 else None
    print(f"  {name:28s} dimV={n:2d} dimg={d:2d} Hom={h} iso={iso} Killing={ks}")
def formC(B):
    n=B.shape[0];I,K=stab_ops(n);return stab_matrix(B,n,I,K)
def commC(J):
    n=J.shape[0];I=np.eye(n);return np.kron(I,J.T)-np.kron(J,I)   # XJ-JX
for n in range(2,8):
    for q in range(0,n//2+1):
        B=np.diag([1.]*(n-q)+[-1.]*q)
        report(f"so({n-q},{q})",n,stab_general(n,[formC(B)]))
for m in (1,2,3):
    J=np.block([[np.zeros((m,m)),np.eye(m)],[-np.eye(m),np.zeros((m,m))]])
    report(f"sp({2*m},R)",2*m,stab_general(2*m,[formC(J)]))
for m in (1,2,3):
    J=np.block([[np.zeros((m,m)),-np.eye(m)],[np.eye(m),np.zeros((m,m))]])
    for q in range(0,m//2+1):
        g=np.diag(([1.]*(m-q)+[-1.]*q)*2)
        report(f"u({m-q},{q}) on R^{2*m}",2*m,stab_general(2*m,[formC(g),commC(J)]))
    g=np.diag([1.]*m+[-1.]*m)   # Re of sum z_i^2
    report(f"so({m},C) on R^{2*m}",2*m,stab_general(2*m,[formC(g),commC(J)]))
m=4;J=np.block([[np.zeros((m,m)),-np.eye(m)],[np.eye(m),np.zeros((m,m))]]);g=np.diag([1.]*m+[-1.]*m)
report("so(4,C) on R^8",8,stab_general(8,[formC(g),commC(J)]))
for n in (1,2,3):
    report(f"gl({n},R) (no structure)",n,[np.eye(n*n)[i].reshape(n,n) for i in range(n*n)])
    if n>1: report(f"sl({n},R) (volume)",n,stab_general(n,[np.eye(n).reshape(1,-1)]))
# G2: stabiliser of the 3-form on R^7
tri=[(0,1,2),(0,3,4),(0,5,6),(1,3,5),(1,4,6),(2,3,6),(2,4,5)];sg=[1,1,1,1,-1,-1,-1]
phi=np.zeros((7,7,7))
for (a,b,c),s in zip(tri,sg):
    for p,sn in [((a,b,c),1),((b,c,a),1),((c,a,b),1),((b,a,c),-1),((a,c,b),-1),((c,b,a),-1)]: phi[p]=s*sn
rows=[]
for k in range(49):
    X=np.zeros(49);X[k]=1;X=X.reshape(7,7)
    r=np.einsum('ai,abc->ibc',X,phi)+np.einsum('bi,abc->aic',X,phi)+np.einsum('ci,abc->abi',X,phi)
    rows.append(r.reshape(-1))
report("g2 (3-form on R^7)",7,stab_general(7,[np.array(rows).T]))
# ---------- Part C
print("\n=== PART C: all simple Lie algebras rank<=8: dim of smallest nontrivial irrep vs dim g (Weyl formula, exact) ===")
def cartan(t,n):
    A=[[0]*n for _ in range(n)]
    for i in range(n): A[i][i]=2
    def link(i,j,a=-1,b=-1): A[i][j]=a;A[j][i]=b
    if t in 'ABC':
        for i in range(n-1): link(i,i+1)
        if t=='B' and n>1: A[n-2][n-1]=-2
        if t=='C' and n>1: A[n-1][n-2]=-2
    elif t=='D':
        for i in range(n-2): link(i,i+1)
        link(n-3,n-1)
    elif t=='E':
        for i in range(n-2): link(i,i+1)
        link(2,n-1)
    elif t=='F':
        link(0,1);link(1,2);link(2,3);A[1][2]=-2
    elif t=='G':
        link(0,1);A[0][1]=-3
    return A
def analyse(t,n):
    A=cartan(t,n)
    # D_j with A[i][j]*D[j]=A[j][i]*D[i]
    D=[None]*n;D[0]=Fraction(1);stack=[0]
    while stack:
        i=stack.pop()
        for j in range(n):
            if j!=i and A[i][j]!=0 and D[j] is None:
                D[j]=D[i]*Fraction(A[j][i],A[i][j]);stack.append(j)
    roots=set(tuple(1 if k==i else 0 for k in range(n)) for i in range(n));frontier=list(roots)
    while frontier:
        b=frontier.pop()
        for i in range(n):
            c=sum(b[j]*A[j][i] for j in range(n))
            nb=list(b);nb[i]-=c;nb=tuple(nb)
            if nb not in roots: roots.add(nb);frontier.append(nb)
    pos=[r for r in roots if all(x>=0 for x in r)]
    assert len(roots)==2*len(pos)
    def dim(m):
        num=Fraction(1)
        for a in pos:
            num*=Fraction(sum(a[i]*(m[i]+1)*D[i] for i in range(n)),sum(a[i]*D[i] for i in range(n)))
        assert num.denominator==1;return int(num)
    fund=[dim([1 if k==i else 0 for k in range(n)]) for i in range(n)]
    return 2*len(pos)+n,min(fund)
types=[('A',n) for n in range(1,9)]+[('B',n) for n in range(2,9)]+[('C',n) for n in range(3,9)]+[('D',n) for n in range(4,9)]+[('E',6),('E',7),('E',8),('F',4),('G',2)]
for t,n in types:
    dg,dmin=analyse(t,n)
    print(f"  {t}{n}: dim g={dg:4d}  min nontrivial irrep={dmin:4d}  {'<== ADJOINT IS MINIMAL' if dmin==dg else ''}")
