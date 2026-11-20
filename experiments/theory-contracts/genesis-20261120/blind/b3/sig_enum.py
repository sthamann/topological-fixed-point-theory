"""H1 test: signature of the fold-back form Q(x)=tr(x o x) over classes of algebras.
Question: does 'exactly one direction of odd sign' (min(p,q)=1) single out a dimension?"""
import numpy as np, itertools
def sig(G,tol=1e-9):
    w=np.linalg.eigvalsh((G+G.T)/2); return (int((w>tol).sum()),int((w<-tol).sum()),int((abs(w)<=tol).sum()))
I2=np.eye(2); J=np.array([[0,-1.],[1,0]])
def Lq():
    # left mult matrices of 1,i,j,k on H=R^4
    def mul(a,b):
        a0,a1,a2,a3=a;b0,b1,b2,b3=b
        return np.array([a0*b0-a1*b1-a2*b2-a3*b3,a0*b1+a1*b0+a2*b3-a3*b2,a0*b2-a1*b3+a2*b0+a3*b1,a0*b3+a1*b2-a2*b1+a3*b0])
    E=np.eye(4); return [np.array([mul(E[u],E[v]) for v in range(4)]).T for u in range(4)]
UN={'R':[np.eye(1)],'C':[I2,J],'H':Lq()}
def basis(n,K):
    B=[]
    for i in range(n):
        for j in range(n):
            E=np.zeros((n,n));E[i,j]=1
            for u in UN[K]: B.append(np.kron(E,u))
    return B
def trace_form(B): 
    F=np.array([b.flatten() for b in B]); Ft=np.array([b.T.flatten() for b in B]); return F@Ft.T   # tr(a b)
def regular_form(B):
    F=np.array([b.flatten() for b in B]).T; d=len(B)
    L=[np.linalg.lstsq(F,np.array([(a@b).flatten() for b in B]).T,rcond=None)[0] for a in B]
    return np.array([[np.trace(L[a]@L[b]) for b in range(d)] for a in range(d)])
print("== A. simple real associative algebras M_n(K): signature (pos,neg,null) of tr(x^2)")
hits=[];simple={}
for K in 'RCH':
    for n in range(1,9):
        B=basis(n,K); s=sig(trace_form(B)); simple[(n,K)]=(len(B),s)
        extra=''
        if len(B)<=64:
            s2=sig(regular_form(B)); extra=' regular-rep check '+('OK' if s2==s else 'MISMATCH %s'%(s2,))
        flag=' <== min(p,q)=1' if min(s[0],s[1])==1 else ''
        if flag: hits.append((n,K,len(B),s))
        print(f"M_{n}({K}) dim={len(B):3d} sig={s}{extra}{flag}")
print("index-1 simple algebras:",hits)
print("closed forms: M_n(R):(n(n+1)/2,n(n-1)/2)  M_n(C):(n^2,n^2)  M_n(H):(n(2n-1),n(2n+1))")
print("\n== B. NEGATIVE CONTROL 1: all semisimple (products) of real dim<=8 with min(p,q)=1")
blocks=[(f"M{n}({K})",d,s) for (n,K),(d,s) in simple.items() if d<=8]
res={}
def rec(start,dim,p,q,names):
    if names and min(p,q)==1: res.setdefault(dim,[]).append(('x'.join(names),(p,q)))
    for i in range(start,len(blocks)):
        nm,d,s=blocks[i]
        if dim+d<=8: rec(i,dim+d,p+s[0],q+s[1],names+[nm])
rec(0,0,0,0,[])
for d in sorted(res): print(f"dim {d}: {len(res[d])} algebras:",res[d])
print("\n== C. NEGATIVE CONTROL 2: drop associativity (Cayley-Dickson tower), Q=Re(x^2)")
def cd_mul(a,b):
    n=len(a)
    if n==1: return a*b
    h=n//2; a1,a2,b1,b2=a[:h],a[h:],b[:h],b[h:]
    conj=lambda x: np.concatenate([x[:1],-x[1:]])
    return np.concatenate([cd_mul(a1,b1)-cd_mul(conj(b2),a2), cd_mul(b2,a1)+cd_mul(a2,conj(b1))])
for k,name in enumerate(['R','C','H','O','S(sedenions)']):
    n=2**k;E=np.eye(n)
    G=np.array([[ (cd_mul(E[a],E[b])[0]+cd_mul(E[b],E[a])[0])/2 for b in range(n)] for a in range(n)])
    assoc=all(np.allclose(cd_mul(cd_mul(E[a],E[b]),E[c]),cd_mul(E[a],cd_mul(E[b],E[c]))) for a in range(n) for b in range(n) for c in range(n))
    print(f"{name:13s} dim={n:2d} sig={sig(G)} associative={assoc}")
print("\n== D. NEGATIVE CONTROL 3: a different natural quadratic form e2(x)=((tr x)^2-tr x^2)/2 on M_n(R)")
for n in range(2,7):
    B=basis(n,'R'); t=np.array([np.trace(b) for b in B]); print(f"M_{n}(R) e2 sig={sig((np.outer(t,t)-trace_form(B))/2)}")
print("\n== E. NEGATIVE CONTROL 4: group ontology. Q(x)=coeff of e in x*x on R[G], all groups of order<=8")
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
pm=lambda a,b: tuple(a[i] for i in b)
def cyc(n,off=0,N=None):
    N=N or n; p=list(range(N))
    for i in range(n): p[off+i]=off+(i+1)%n
    return tuple(p)
def qm(a,b):
    a0,a1,a2,a3=a;b0,b1,b2,b3=b
    return (a0*b0-a1*b1-a2*b2-a3*b3,a0*b1+a1*b0+a2*b3-a3*b2,a0*b2-a1*b3+a2*b0+a3*b1,a0*b3+a1*b2-a2*b1+a3*b0)
G={f"Z{n}":([cyc(n)],pm) for n in range(1,9)}
G["Z2xZ2"]=([cyc(2,0,4),cyc(2,2,4)],pm);G["Z2xZ4"]=([cyc(2,0,6),cyc(4,2,6)],pm)
G["Z2^3"]=([cyc(2,0,6),cyc(2,2,6),cyc(2,4,6)],pm);G["S3"]=([(1,0,2),(1,2,0)],pm)
G["D4"]=([(1,2,3,0),(0,3,2,1)],pm);G["Q8"]=([(0,1,0,0),(0,0,1,0)],qm)
for nm,(gens,mul) in G.items():
    el=closure(gens,mul); e=[x for x in el if all(mul(x,y)==y for y in el)][0]
    M=np.array([[1.0 if mul(a,b)==e else 0 for b in el] for a in el]); s=sig(M)
    print(f"{nm:6s} |G|={len(el)} sig={s}"+(' <== min(p,q)=1' if min(s[:2])==1 else ''))
