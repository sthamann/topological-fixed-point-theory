"""C1: 'syndrome = generator' count over the class (K, p, N).
K in {R,C,H}, p-norm on K^N, N=2,3,4.
g    = dim of Lie algebra of K-linear isometries of the p-norm (numerical nullspace)
r    = dim of span of outcome statistics f_{U,i}(x)=|(Ux)_i|^p, U in isometry group
t    = dim of generators acting trivially on every statistic
Rule: g - t == r - 1  (effective continuous reversible moves == nontrivial statistics)
"""
import numpy as np, itertools
rng=np.random.default_rng(7)
Lq=[np.eye(4),
    np.array([[0,-1,0,0],[1,0,0,0],[0,0,0,-1],[0,0,1,0.]]),
    np.array([[0,0,-1,0],[0,0,0,1],[1,0,0,0],[0,-1,0,0.]]),
    np.array([[0,0,0,-1],[0,0,-1,0],[0,1,0,0],[1,0,0,0.]])]
Lc=[np.eye(2),np.array([[0,-1],[1,0.]])]
UN={'R':[np.eye(1)],'C':Lc,'H':Lq}
def basis(K,N):
    B=[]
    for a in range(N):
        for b in range(N):
            E=np.zeros((N,N));E[a,b]=1
            for L in UN[K]: B.append(np.kron(E,L))
    return B
def expm(A):
    n=np.linalg.norm(A,2); s=max(0,int(np.ceil(np.log2(max(n,1e-16))))+4)
    A=A/2**s; X=np.eye(len(A)); T=np.eye(len(A))
    for k in range(1,25): T=T@A/k; X=X+T
    for _ in range(s): X=X@X
    return X
def null(M,tol=1e-9):
    u,s,vt=np.linalg.svd(M,full_matrices=True)
    r=int((s>tol*s[0]).sum()) if s.size and s[0]>0 else 0
    return vt[r:].T,r
def blocks(x,d): return x.reshape(-1,d)
def run(K,N,p):
    d={'R':1,'C':2,'H':4}[K]; B=basis(K,N); nb=len(B)
    rows=[]
    for _ in range(6*nb+20):
        x=rng.normal(size=d*N); xb=blocks(x,d); w=np.linalg.norm(xb,axis=1)**(p-2)
        rows.append([float((w[:,None]*xb*blocks(Bk@x,d)).sum()) for Bk in B])
    V,_=null(np.array(rows)); g=V.shape[1]
    Lb=[sum(V[k,j]*B[k] for k in range(nb)) for j in range(g)]
    # check: closure sanity — exp of generator preserves norm
    Us=[np.eye(d*N)]
    for _ in range(60):
        if g==0: break
        c=rng.normal(size=g); Us.append(expm(sum(ci*Li for ci,Li in zip(c,Lb))))
    err=0.
    for U in Us[1:4]:
        x=rng.normal(size=d*N)
        nrm=lambda v:(np.linalg.norm(blocks(v,d),axis=1)**p).sum()
        err=max(err,abs(nrm(U@x)-nrm(x))/nrm(x))
    X=rng.normal(size=(300,d*N))
    F=[]
    for U in Us:
        Y=X@U.T
        Yb=np.linalg.norm(Y.reshape(300,N,d),axis=2)**p
        F.append(Yb)
    F=np.concatenate(F,axis=1)
    s=np.linalg.svd(F,compute_uv=False); r=int((s>1e-9*s[0]).sum())
    t=0
    if g>0:
        rows=[]
        for U in Us[:12]:
            for _ in range(8):
                x=rng.normal(size=d*N); yb=blocks(U@x,d); w=np.linalg.norm(yb,axis=1)**(p-2)
                Z=[blocks(U@(Lk@x),d) for Lk in Lb]
                for i in range(N):
                    rows.append([float(w[i]*(yb[i]*Zk[i]).sum()) for Zk in Z])
        R=np.array(rows)
        if np.abs(R).max()<1e-9: t=g
        else:
            s2=np.linalg.svd(R,compute_uv=False); t=g-int((s2>1e-9*s2[0]).sum())
    return g,t,r,err
print(f"{'K':>2} {'N':>2} {'p':>4} | {'g':>3} {'t':>3} {'eff=g-t':>7} | {'r':>3} {'r-1':>4} | rule(eff==r-1) raw(g==r)  isometry_err")
hits=[]
for K in 'RCH':
    for N in (2,3,4):
        for p in (1.0,1.5,2.0,3.0,4.0):
            g,t,r,err=run(K,N,p)
            ok=(g-t==r-1); raw=(g==r)
            if ok: hits.append((K,N,p))
            print(f"{K:>2} {N:>2} {p:>4} | {g:>3} {t:>3} {g-t:>7} | {r:>3} {r-1:>4} | {str(ok):>14} {str(raw):>9}  {err:.1e}")
print("\nRule satisfied for:",hits)
print("Closed forms expected at p=2: R: g=N(N-1)/2,r=N(N+1)/2; C: g=N^2,r=N^2,t=1; H: g=N(2N+1),r=N(2N-1)")
