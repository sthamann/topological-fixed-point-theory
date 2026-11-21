#!/usr/bin/env python3
"""Exact Fibonacci B4 -> B5 enlargement and invariant-state obstruction.

Input: tau*tau = 1+tau, its unitary Fibonacci F/R data, a chosen final
tau charge. a=phi^-1 obeys a^2+a-1=0. A diagonal similarity removes
sqrt(a), leaving projectors over Z[a]/(a^2+a-1).
Modular primes certify LOWER BOUNDS on the characteristic-zero algebra;
the 2+3 block upper bound and 5x5 full dimension give exact equalities.
All F/TL algebra identities are checked symbolically over Z[a].
"""
from collections import deque
Z=(0,0); U=(1,0); a=(0,1); a2=(1,-1)
def plus(x,y):return (x[0]+y[0],x[1]+y[1])
def times(x,y):
    u,v=x;w,z=y
    return u*w+v*z,u*z+v*w-v*z
def ident(n):return [[U if i==j else Z for j in range(n)] for i in range(n)]
def prod(X,Y):
    n=len(X)
    return [[sum2(times(X[i][k],Y[k][j]) for k in range(n)) for j in range(n)] for i in range(n)]
def sum2(xs):
    t=Z
    for x in xs:t=plus(t,x)
    return t
def scaled(X,t):return [[times(x,t) for x in row] for row in X]
def paths(n=5,end=1):
    p=[(0,)]
    for _ in range(n):p=[r+(v,) for r in p for v in ((1,) if r[-1]==0 else (0,1))]
    return sorted(r for r in p if r[-1]==end)
def projector(i,n=5,end=1):
    pp=paths(n,end); idx={v:j for j,v in enumerate(pp)}; N=len(pp)
    M=[[Z for _ in range(N)] for _ in range(N)]
    for col,p in enumerate(pp):
        if i==1:M[col][col]=U if p[2]==0 else Z
        elif p[i-1]==0 and p[i+1]==0:M[col][col]=U
        elif p[i-1]==p[i+1]==1:
            for mid in (0,1):
                row=idx[p[:i]+(mid,)+p[i+1:]]
                M[row][col]=a2 if p[i]==0 else a
    return M
def finite_field_rank(gens,prime,root):
    assert (root*root+root-1)%prime==0
    n=len(gens[0]); basis={}; q=deque([(ident(n),0)]); depths=[]
    while q:
        M,depth=q.popleft()
        v=[(u+root*w)%prime for row in M for u,w in row]
        for pivot in sorted(basis):
            if v[pivot]:
                c=v[pivot]; v=[(x-c*y)%prime for x,y in zip(v,basis[pivot])]
        pivot=next((i for i,x in enumerate(v) if x),None)
        if pivot is None:continue
        inv=pow(v[pivot],-1,prime)
        basis[pivot]=[(x*inv)%prime for x in v]
        depths.append(depth)
        if len(depths)==n*n:break
        for G in gens:q.append((prod(M,G),depth+1))
    return len(depths),max(depths)
def main():
    pp=paths(); P=[projector(i) for i in range(1,5)]
    assert len(pp)==5 and [p[4] for p in pp].count(0)==2
    for i,Q in enumerate(P):
        assert prod(Q,Q)==Q
        for j in range(i+2,4):assert prod(Q,P[j])==prod(P[j],Q)
    for i in range(3):
        assert prod(prod(P[i],P[i+1]),P[i])==scaled(P[i],a2)
        assert prod(prod(P[i+1],P[i]),P[i+1])==scaled(P[i+1],a2)
    print("Exact golden-field idempotence/locality/Temperley-Lieb: PASS")
    for prime,root in [(11,3),(19,4),(29,5)]:
        r3,d3=finite_field_rank(P[:3],prime,root)
        r4,d4=finite_field_rank(P,prime,root)
        assert (r3,r4)==(13,25)
        print(f"GF({prime}), a={root}: B4 algebra rank={r3}, B5 rank={r4}; max word lengths={d3},{d4}")
    split=[p[4] for p in pp]
    cross=sum(1 for i in range(5) for j in range(5) if split[i]!=split[j] and P[3][i][j]!=Z)
    assert cross==4
    print("Braid B4 mixes charge sectors; nonzero off-block entries:",cross)
    print("Exactly: Alg(B1,B2,B3)=M2(C)+M3(C); Alg(B1..B4)=M5(C)")
    print("NO-GO: full braid-invariant rho must be I5/5, so sigma_t^rho=id.")
    print("NO-GO: B4 invariance alone permits a continuum of faithful block weights.")
    print("PASS — native P1, 3+1, physical time and gravity not derived.")
if __name__=="__main__":main()
