"""H2 Parallax: X finite set, projections = partitions; rule: every pair of projections is
jointly complete and non-redundant (X -> P_i x P_j bijective); take a MAXIMAL family.
=> |X|=n^2; k partitions <=> k-2 mutually orthogonal Latin squares (MOLS). k<=n+1.
Exhaustive over ALL reduced Latin squares of order n=2..6 (every MOLS set is isotopic to one
whose first square is reduced). Negative controls: n=6, and the cyclic square Z4."""
import itertools, time, collections
def reduced_squares(n):
    L=[[None]*n for _ in range(n)]
    for j in range(n): L[0][j]=j
    for i in range(n): L[i][0]=i
    rows=[set(r)-{None} for r in L]; cols=[{L[i][j] for i in range(n) if L[i][j] is not None} for j in range(n)]
    cells=[(i,j) for i in range(1,n) for j in range(1,n)]
    def rec(k):
        if k==len(cells):
            yield tuple(tuple(r) for r in L);return
        i,j=cells[k]
        for s in range(n):
            if s not in rows[i] and s not in cols[j]:
                L[i][j]=s;rows[i].add(s);cols[j].add(s)
                yield from rec(k+1)
                rows[i].remove(s);cols[j].remove(s);L[i][j]=None
    yield from rec(0)
def transversals(L,n):
    out=[]
    def rec(i,usedc,useds,cur):
        if i==n: out.append(tuple(cur));return
        for j in range(n):
            if not usedc>>j&1 and not useds>>L[i][j]&1:
                cur.append(j);rec(i+1,usedc|1<<j,useds|1<<L[i][j],cur);cur.pop()
    rec(0,0,0,[])
    return out
def mates(L,n):
    """all orthogonal mates, normalised: symbol of mate at (0,j) is j"""
    T=transversals(L,n)
    by0=collections.defaultdict(list)
    for t in T: by0[t[0]].append(t)
    res=[]
    def rec(j,chosen,occ):
        if j==n:
            M=[[None]*n for _ in range(n)]
            for s,t in enumerate(chosen):
                for i in range(n): M[i][t[i]]=s
            res.append(tuple(tuple(r) for r in M));return
        for t in by0[j]:
            if all(not occ[i]>>t[i]&1 for i in range(n)):
                rec(j+1,chosen+[t],[occ[i]|1<<t[i] for i in range(n)])
    rec(0,[],[0]*n)
    return len(T),res
def orth(A,B,n): return len({(A[i][j],B[i][j]) for i in range(n) for j in range(n)})==n*n
def maxclique(ms,n):
    adj=[[orth(a,b,n) for b in ms] for a in ms]
    best=0
    def rec(cl,cand):
        nonlocal best
        best=max(best,len(cl))
        for idx,c in enumerate(cand):
            rec(cl+[c],[d for d in cand[idx+1:] if adj[c][d]])
    rec([],list(range(len(ms))))
    return best
def loop_info(L,n):
    assoc=all(L[L[a][b]][c]==L[a][L[b][c]] for a in range(n) for b in range(n) for c in range(n))
    if not assoc: return "non-associative loop"
    def order(a):
        k=1;x=a
        while x!=0: x=L[x][a];k+=1
        return k
    comm=all(L[a][b]==L[b][a] for a in range(n) for b in range(n))
    return "group, element orders "+str(sorted(order(a) for a in range(n)))+(" abelian" if comm else " nonabelian")
print("n | #reduced squares | #with >=1 transversal | #with orthogonal mate | max MOLS | max k (projections) | n+1")
summary=[]
for n in range(2,7):
    t0=time.time();cnt=0;wt=0;wm=0;best=1;detail=collections.Counter()
    for L in reduced_squares(n):
        cnt+=1
        nt,ms=mates(L,n)
        if nt: wt+=1
        mc=0
        if ms:
            wm+=1; mc=maxclique(ms,n) if n<=5 else 1
        best=max(best,1+mc)
        detail[(loop_info(L,n),nt,len(ms),1+mc)]+=1
    print(f"{n} | {cnt} | {wt} | {wm} | {best} | {best+2} | {n+1}   ({time.time()-t0:.1f}s)")
    summary.append((n,detail))
print("\nDetail: (algebraic type of the reduced square as a loop, #transversals, #normalised mates, size of largest MOLS set containing it) : #squares")
for n,detail in summary:
    print(f" n={n}")
    if n<=5:
        for k in sorted(detail,key=str): print("    ",k,":",detail[k])
    else:
        agg=collections.Counter()
        for (ty,nt,nm,mc),c in detail.items(): agg[(ty.split(',')[0],'transversals:'+str(nt),'mates:'+str(nm))]+=c
        for k in sorted(agg,key=str): print("    ",k,":",agg[k])
