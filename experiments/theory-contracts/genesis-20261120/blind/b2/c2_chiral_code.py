"""C2: smallest chiral, anomaly-free charge code.
A 'state' is a joint eigenvalue string q in A^n of n commuting charges (alphabet A).
A chiral multiset is an ODD integer function x on A^n (x(-q) = -x(q)); x(q)=+1 means a
state q with no partner at -q.  Condition ("invisible to all probes of order <= m"):
  sum_q x(q) * (monomial of degree d in q) = 0   for all monomials of odd degree d <= m.
m=3 is the cubic (triangle) condition incl. the linear one.
Exhaustive over ALL x with values in {-1,0,1} (not only linear codes).
"""
import numpy as np, itertools, sympy as sp, time
def reps(A,n):
    pts=[q for q in itertools.product(A,repeat=n)]
    R=[]
    for q in pts:
        nq=tuple(-a for a in q)
        if q==nq: continue
        if q>nq: R.append(q)
    return R
def monos(n,m):
    M=[]
    for d in range(1,m+1,2):
        M+=list(itertools.combinations_with_replacement(range(n),d))
    return M
def cmatrix(A,n,m):
    R=reps(A,n); M=monos(n,m)
    C=np.array([[int(np.prod([q[i] for i in mo])) for q in R] for mo in M],dtype=np.int64)
    return R,C
def enum(A,n,m):
    """count all nonzero x in {-1,0,1}^pairs with C x = 0, meet in the middle."""
    R,C=cmatrix(A,n,m); k=len(R)
    h=k//2
    def table(cols):
        t={}
        for v in itertools.product((-1,0,1),repeat=len(cols)):
            key=(C[:,cols]@np.array(v,dtype=np.int64)).tobytes() if cols else np.zeros(C.shape[0],dtype=np.int64).tobytes()
            t.setdefault(key,[]).append(v)
        return t
    return R,C,h
def enum_fast(A,n,m,keep=4):
    R,C=cmatrix(A,n,m); k=len(R); h=k//2
    def allv(c):
        g=np.array(np.meshgrid(*[(-1,0,1)]*c,indexing='ij'),dtype=np.int64).reshape(c,-1) if c else np.zeros((0,1),dtype=np.int64)
        return g
    VL=allv(h); VR=allv(k-h)
    SL=C[:,:h]@VL; SR=-(C[:,h:]@VR)
    # hash rows
    rng=np.random.default_rng(1); w=rng.integers(1,2**40,size=C.shape[0]).astype(np.int64)
    hl=w@SL; hr=w@SR
    order=np.argsort(hr); hrs=hr[order]
    lo=np.searchsorted(hrs,hl,'left'); hi=np.searchsorted(hrs,hl,'right')
    sols=[]; cnt=0
    for i in np.nonzero(hi>lo)[0]:
        for j in order[lo[i]:hi[i]]:
            if np.array_equal(SL[:,i],SR[:,j]):
                x=np.concatenate([VL[:,i],VR[:,j]])
                if np.any(x):
                    cnt+=1
                    if len(sols)<keep: sols.append(x)
    return R,C,cnt,sols
def qrank(C): return sp.Matrix(C.tolist()).rank() if C.size else 0

print("=== Part A: binary alphabet {+1,-1}^n, cubic condition m=3 (exhaustive, x in {-1,0,1}) ===")
for n in range(1,6):
    t=time.time(); R,C,cnt,sols=enum_fast((1,-1),n,3)
    print(f"n={n}: antipodal pairs={len(R)}, conditions={C.shape[0]}, rank={qrank(C)}, "
          f"Q-nullity={len(R)-qrank(C)}, nonzero solutions={cnt}  ({time.time()-t:.1f}s)")
    for x in sols:
        S=[tuple(int(s*a) for a in q) for q,s in zip(R,x) if s!=0]
        par=set(int(np.prod(q)) for q in S)
        print(f"   solution: #states={len(S)}, product-of-signs over states={par}  (single parity class => even/odd-weight code)")
print("n=6: nullity and exact count via x = chi_all * sum a_i s_i, a_i in {0,+-1/2,+-1}:")
R,C=cmatrix((1,-1),6,3); print("   Q-nullity =",len(R)-qrank(C))
pts=np.array(list(itertools.product((1,-1),repeat=6)))
cnt=0; sizes={}
for a in itertools.product((0,.5,-.5,1,-1),repeat=6):
    v=pts@np.array(a)
    if np.any(v) and np.all(np.isin(v,(-1,0,1))):
        cnt+=1; k=int(np.count_nonzero(v))//2; sizes[k]=sizes.get(k,0)+1
print("   nonzero {-1,0,1} solutions =",cnt," by #chiral states:",sizes)

print("\n=== Part B: order m of the invisibility condition vs minimal length (Q-nullity>0) ===")
for m in (1,3,5):
    row=[]
    for n in range(1,9 if m<5 else 9):
        R,C=cmatrix((1,-1),n,m)
        Cn=C.astype(float); r=np.linalg.matrix_rank(Cn)
        row.append((n,len(R)-r))
    nmin=[n for n,z in row if z>0][0]
    print(f"m={m}: (n, nullity) = {row}  -> minimal n = {nmin}, states = 2^{nmin-1} = {2**(nmin-1)}")
print("If m = D/2+1 (order of the chiral anomaly in D dims, D=0 mod 4): n_min = D/2+3.")
for Dm in (4,8,12,16):
    print(f"   D={Dm}: n_min={Dm//2+3}, D+1={Dm+1}, equal={Dm//2+3==Dm+1}")

print("\n=== Part C: negative controls (other alphabets, cubic condition m=3, exhaustive) ===")
for A,nmax in (((1,0,-1),3),((3,1,-1,-3),2),((2,1,-1,-2),2),((2,1,0,-1,-2),2)):
    for n in range(1,nmax+1):
        R,C=cmatrix(A,n,3)
        if len(R)>26: print(f"alphabet {A} n={n}: skipped ({len(R)} pairs)"); continue
        t=time.time(); R,C,cnt,sols=enum_fast(A,n,3,keep=1)
        msg=f"alphabet {A} n={n}: pairs={len(R)}, nonzero solutions={cnt}"
        if cnt:
            # smallest solution
            msg+=f"; example with {int(np.count_nonzero(sols[0]))} states: "+str([tuple(int(s*a) for a in q) for q,s in zip(R,sols[0]) if s!=0])
        print(msg,f"({time.time()-t:.1f}s)")
