"""C2b: negative control for C2. Alphabet {-1,0,1}, n=4, cubic condition.
(i) verify the 15-state set = weights of SU(5) 5bar+10 in Dynkin-label basis is a solution;
(ii) exhaustive search for ALL chiral solutions with <= 8 states (signed subsets of the 40 antipodal pairs)."""
import numpy as np, itertools, time
from c2_chiral_code import cmatrix
R,C=cmatrix((1,0,-1),4,3); idx={q:i for i,q in enumerate(R)}
print("pairs",len(R),"conditions",C.shape[0],"rank",np.linalg.matrix_rank(C.astype(float)))
f=[(1,0,0,0),(-1,1,0,0),(0,-1,1,0),(0,0,-1,1),(0,0,0,-1)]
ten=[tuple(a+b for a,b in zip(f[i],f[j])) for i in range(5) for j in range(i+1,5)]
S=ten+[tuple(-a for a in w) for w in f]
x=np.zeros(len(R),dtype=np.int64)
for q in S:
    if q in idx: x[idx[q]]+=1
    else: x[idx[tuple(-a for a in q)]]-=1
print("SU(5) 5bar+10: states",len(S),"alphabet ok",all(abs(a)<=1 for q in S for a in q),
      "max|x|",abs(x).max(),"anomaly vector zero:",not np.any(C@x))
# exhaustive supports <= 8 via meet in the middle on signed subsets of size <= 4
t=time.time(); k=len(R)
w=np.random.default_rng(3).integers(1,2**36,size=C.shape[0]).astype(np.int64)
hc=w@C  # hash per column
tab={}
for s in range(0,5):
    for sub in itertools.combinations(range(k),s):
        for sg in itertools.product((1,-1),repeat=s):
            h=int(sum(g*hc[i] for g,i in zip(sg,sub)))
            tab.setdefault(h,[]).append((sub,sg))
best={}
for h,lst in tab.items():
    if -h in tab:
        for (sa,ga) in lst:
            for (sb,gb) in tab[-h]:
                x=np.zeros(k,dtype=np.int64)
                for g,i in zip(ga,sa): x[i]+=g
                for g,i in zip(gb,sb): x[i]+=g
                if np.any(x) and abs(x).max()<=1 and not np.any(C@x):
                    best.setdefault(int(np.count_nonzero(x)),set()).add(tuple(x))
print("solutions with <=8 states, by size:",{s:len(v) for s,v in sorted(best.items())},f"({time.time()-t:.0f}s)")
if best:
    s=min(best); x=np.array(sorted(best[s])[0])
    print("smallest example:",[tuple(int(g*a) for a in q) for q,g in zip(R,x) if g])
