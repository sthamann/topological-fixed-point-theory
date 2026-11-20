import numpy as np, itertools, math
from collections import Counter
n=4; N=16
def key(v):
    i=np.flatnonzero(np.abs(v)>1e-9)[0]; w=v/v[i]; w=w/np.linalg.norm(w)
    return tuple(np.round(w.real,6))+tuple(np.round(w.imag,6))
H=np.array([[1,1],[1,-1]])/math.sqrt(2); S=np.diag([1,1j])
def one(g,q):
    m=np.eye(1)
    for i in range(n): m=np.kron(m,g if i==q else np.eye(2))
    return m
def cz(a,b):
    d=np.ones(N)
    for x in range(N):
        if (x>>(n-1-a))&1 and (x>>(n-1-b))&1: d[x]=-1
    return np.diag(d)
gens=[one(H,q) for q in range(n)]+[one(S,q) for q in range(n)]+[cz(a,b) for a,b in itertools.combinations(range(n),2)]
v0=np.zeros(N,complex); v0[0]=1
seen={key(v0):v0}; front=[v0]
while front:
    nxt=[]
    for v in front:
        for g in gens:
            w=g@v; k=key(w)
            if k not in seen: seen[k]=w; nxt.append(w)
    front=nxt
print("Stabilisatorzustände:",len(seen))
perms=list(itertools.permutations(range(n)))
def perm(v,p):
    w=np.zeros(N,complex)
    for x in range(N):
        bits=[(x>>(n-1-i))&1 for i in range(n)]; y=0
        for i in range(n): y=(y<<1)|bits[p[i]]
        w[y]=v[x]
    return w
res=Counter(); ex={}
for v in seen.values():
    v=v/np.linalg.norm(v)
    if not all(np.allclose(perm(v,p),v,atol=1e-7) for p in [(1,0,2,3),(1,2,3,0)]): continue
    c=[sum(v[x] for x in range(N) if bin(x).count("1")==k)/math.comb(4,k) for k in range(5)]  # a_k
    e,d,cc,b,a=c
    I=a*e-4*b*d+3*cc*cc; J=a*cc*e+2*b*cc*d-a*d*d-b*b*e-cc**3; disc=I**3-27*J*J
    if abs(disc)<1e-9: t="mehrfache Marken"
    elif abs(J)<1e-9: t="vier verschiedene, harmonisch (Quadrat, j=1728)"
    elif abs(I)<1e-9: t="vier verschiedene, Tetraeder (j=0)"
    else: t="vier verschiedene, sonst j=%.3f"%(1728*(I**3/disc)).real
    res[t]+=1; ex.setdefault(t,np.round(np.array(c)/max(abs(x) for x in c),3))
print("symmetrische Stabilisatorzustände:",sum(res.values()))
for t,k in res.items(): print(k,t,ex[t])
# Gegenprobe: Tetraederzustand ist kein Stabilisatorzustand (nicht in Liste, da I=0 nicht vorkommt)
