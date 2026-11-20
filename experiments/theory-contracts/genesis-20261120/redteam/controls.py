import sys, itertools as it
from fractions import Fraction as F
from collections import Counter
import numpy as np, sympy
sys.path.insert(0, "/Users/stefanhamann/Projekte/tfpt-theoryv4/research/originlab")
sys.dont_write_bytecode = True

print("== C1 kernel dimension of odd functions on {+-1}^n blind to given odd orders")
def kern(n, orders):
    words=[w for w in it.product((1,-1),repeat=n) if w[0]==1]
    rows=[[int(np.prod([w[i] for i in idx])) for w in words] for m in orders for idx in it.combinations(range(n),m)]
    return len(words)-(np.linalg.matrix_rank(np.array(rows)) if rows else 0)
for orders in [(1,),(1,3),(1,3,5)]:
    print(orders,{n:int(kern(n,orders)) for n in range(1,9)})

print("== C2 five signs, 16 even words: EVERY linear charge has sum Y = sum Y^3 = 0 (random integer charges)")
W=[w for w in it.product((1,-1),repeat=5) if np.prod(w)==1]
rng=np.random.default_rng(1)
for _ in range(3):
    a=[int(x) for x in rng.integers(-9,10,5)]
    ys=[sum(F(ai*wi,2) for ai,wi in zip(a,w)) for w in W]
    print(a,"sumY",sum(ys),"sumY3",sum(y**3 for y in ys))

print("== C3 multiplicities by split (generic a,b) and with the SM coefficients")
def mult(k,a,b):
    cnt=Counter(a*sum(F(w[i],2) for i in range(k))+b*sum(F(w[i],2) for i in range(k,5)) for w in W)
    return {str(y):cnt[y] for y in sorted(cnt)}
for k in (5,4,3):
    print(f"split {k}+{5-k} generic (a=7,b=11):",sorted(mult(k,F(7),F(11)).values()))
print("3+2 a=-1/3 b=1/2 :",mult(3,F(-1,3),F(1,2)))
print("3+2 a=-1/3 b=-1/2 (flipped SU(5)-like):",mult(3,F(-1,3),F(-1,2)))
print("3+2 a=1/3 b=0 (B-L/2... ):",mult(3,F(1,3),F(0)))
print("3+2 a=1 b=1  :",mult(3,F(1),F(1)))
print("4+1 a=-1/4 b=1:",mult(4,F(-1,4),F(1)))
print("5+0 a=1      :",mult(5,F(1),F(0)))
SM=sorted([F(1,6)]*6+[F(-2,3)]*3+[F(1,3)]*3+[F(-1,2)]*2+[F(1)]+[F(0)])
hits=[]
for k in (3,4,5):
  for an in range(-12,13):
    for bn in range(-12,13):
      a,b=F(an,6),F(bn,6)
      ys=sorted(a*sum(F(w[i],2) for i in range(k))+b*sum(F(w[i],2) for i in range(k,5)) for w in W)
      if ys==SM or sorted(-y for y in ys)==SM: hits.append((k,str(a),str(b)))
print("grid a,b in Z/6, |.|<=2: charges reproducing SM multiset (up to conjugation):",hits)

print("== C4 ternary alphabet {-1,0,1}^n: odd functions blind to all monomials of degree 1 and 3")
def tern(n):
    words=[w for w in it.product((1,0,-1),repeat=n) if any(w) and next(x for x in w if x)==1]
    mons=[e for e in it.product(range(4),repeat=n) if sum(e) in (1,3)]
    A=np.array([[int(np.prod([w[i]**e[i] for i in range(n)])) for w in words] for e in mons])
    return words,A
for n in (1,2,3):
    words,A=tern(n); r=np.linalg.matrix_rank(A)
    print("n",n,"mirror pairs",len(words),"kernel dim",len(words)-r)
words,A=tern(3)
ns=[v*sympy.ilcm(*[x.q for x in v]) for v in sympy.Matrix(A.tolist()).nullspace()]
best=None
for co in it.product(range(-3,4),repeat=len(ns)):
    if not any(co): continue
    v=sum((c*n for c,n in zip(co,ns)),sympy.zeros(len(words),1))
    k=sum(abs(int(x)) for x in v)
    if best is None or k<best[0]: best=(k,[int(x) for x in v])
print("smallest integer solution at n=3 (states with multiplicity):",best[0],[(w,ci) for w,ci in zip(words,best[1]) if ci])
print("== C5 one probe, integer charges: chiral set with sum q = sum q^3 = 0")
q=[1,5,-7,-8,9]; print(q,sum(q),sum(x**3 for x in q), "-> 5 states, not 16")

print("== C6 rule R1 applied to the output of R2: 128 glue vectors of E8 = (16,4)+(16bar,4bar)")
G=[w for w in it.product((1,-1),repeat=8) if np.prod(w)==1]
c=lambda w:int(np.prod(w[:5]))
print("glue vectors",len(G),"class +",sum(c(w)==1 for w in G),"class -",sum(c(w)==-1 for w in G))
print("unlabelled set mirror symmetric:",set(G)=={tuple(-x for x in w) for w in G})
bad={}
for m in (1,3):
    for idx in it.combinations(range(8),m):
        v=sum(c(w)*int(np.prod([w[i] for i in idx])) for w in G)
        if v: bad[idx]=v
print("nonzero probes of order 1,3 for chirality = D5 class:",bad)
fam3=[w for w in G if not (w[5]==w[6]==w[7])]  # drop the family singlet direction (4 -> 3+1), one choice
bad3={idx:v for m in (1,3) for idx in it.combinations(range(8),m) if (v:=sum(c(w)*int(np.prod([w[i] for i in idx])) for w in fam3))}
print("after 4 -> 3 (remove weights with s6=s7=s8):",len(fam3),"states; nonzero probes:",bad3)

print("== C7 trace form signatures of small simple real algebras (regular trace x,y -> tr L_x L_y)")
def sig(basis):
    d=len(basis); n=basis[0].shape[0]
    M=np.array([[np.trace(a@b).real for b in basis] for a in basis])
    ev=np.linalg.eigvalsh((M+M.T)/2); return int((ev>1e-9).sum()),int((ev<-1e-9).sum())
def E(n,i,j,v=1):
    m=np.zeros((n,n),dtype=complex); m[i,j]=v; return m
for n in (2,3): print(f"M_{n}(R)",sig([E(n,i,j) for i in range(n) for j in range(n)]))
for n in (1,2): print(f"M_{n}(C) as real",sig([E(n,i,j,v) for i in range(n) for j in range(n) for v in (1,1j)]))
I2=np.eye(2); qi=np.array([[1j,0],[0,-1j]]); qj=np.array([[0,1],[-1,0]],dtype=complex); qk=qi@qj
print("H",sig([I2.astype(complex),qi,qj,qk]))
H=[I2.astype(complex),qi,qj,qk]
print("M_2(H)",sig([np.kron(E(2,i,j),h) for i in range(2) for j in range(2) for h in H]))

print("== C8 cut universe of the lane")
import genesis_rules as g
for o in g.cuts(8): print({k:(sorted(v) if isinstance(v,(set,frozenset)) else v) for k,v in o.items()})
