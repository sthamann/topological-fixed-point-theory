"""Evaluate the original early and late contact via only even total charge."""
from collections import defaultdict
from contextlib import redirect_stdout
from functools import lru_cache
from itertools import product
from math import factorial
from pathlib import Path
import io,json,runpy,sys

repo=Path(sys.argv[1])
base=repo/"experiments/theory-contracts/chat-triangulation-20261108"
with redirect_stdout(io.StringIO()):
    src=runpy.run_path(str(base/"check_tfpt_family_middle_ordered_cut_20261110.py"))
sp,K,q,r,B,M,V,old=(src[x] for x in ("sp","Q","q","r","d1","M","V","old"))
axes=src["axes"]
axis_metric=src["metric"]
d=q+r
assert d==src["d0"] and sum(x!=0 for x in d)==4
assert all(x.is_Integer for x in d) and sum(d)%2==0
assert q.dot(B)==-2 and q.dot(q)==2

def add_term(dst,key,value):
    if value:dst[key]=dst.get(key,K.zero)+value

def sum_degree(w):
    return sum(n for n,_ in w)

def inner(v,w):
    return src["inner"](v,w)

def add(a,b):
    return src["combine"]([*a.items(),*b.items()])

def subtract(a,b):
    return add(a,{k:-v for k,v in b.items()})

def binomial_factor(s,k):
    if k<0:return 0
    if s==-2:return k+1
    assert s==2
    return (1,-2,1)[k] if k<=2 else 0

@lru_cache(None)
def creation_terms(late,degree):
    # Coefficients of exp(sum_(i,n) (a_i(-n)/n) z_i^n).
    charges=[(1,0),(-1,1),(1,0)] if late else [(-1,1),(1,0),(1,0)]
    cells=[(i,n,j,c) for i,a in enumerate(charges) for n in range(1,degree+1)
           for j,c in zip((4,5),a,strict=True) if c]
    ans={}
    def generate(start,left,word,exponent,value):
        if left==0:
            add_term(ans,(tuple(sorted(word)),exponent),value)
            return
        for k in range(start,len(cells)):
            i,n,j,c=cells[k]
            if n>left:continue
            for multiplicity in range(1,left//n+1):
                shifted=list(exponent)
                shifted[i]+=n*multiplicity
                ratio=K.from_sympy(sp.Rational(c**multiplicity,
                                   factorial(multiplicity)*n**multiplicity))
                generate(k+1,left-n*multiplicity,word+((n,j),)*multiplicity,
                         tuple(shifted),value*ratio)
    generate(0,degree,(),(0,0,0),K.one)
    return {k:v for k,v in ans.items() if v}

@lru_cache(None)
def scalar_coefficient(late,delta):
    # Pair expansions for radial order |z0|>|z1|>|z2|.
    if sum(delta):return 0
    s01,s02,s12=(-2,2,-2) if late else (-2,-2,2)
    total=0
    for k01 in range(max(0,delta[1]),-delta[0]+1):
        k02=-delta[0]-k01
        k12=k01-delta[1]
        total+=binomial_factor(s01,k01)*binomial_factor(s02,k02)*binomial_factor(s12,k12)
    return total

def fused(source,m,late):
    # Initial charge zero. Extract physical conformal modes (m,0,0)
    # or (0,m,0), with field weights (1,2,1) or (2,1,1).
    requested=(-1-m,-2,-1) if late else (-2,-1-m,-1)
    scalar_prefactor=(0,-2,0) if late else (-4,2,0)
    charges=(q,B,q) if late else (B,q,q)
    contractions=[[K.from_sympy(sp.expand(x)) for x in axes.T*a] for a in charges]
    final={}
    for word,c in source.items():
        needed=sum_degree(word)-m-2
        for assignment in product((-1,0,1,2),repeat=len(word)):
            uncontracted=[]
            powers=list(scalar_prefactor)
            coefficient=c
            for (n,j),i in zip(word,assignment,strict=True):
                if i==-1:
                    uncontracted.append((n,j))
                else:
                    coefficient*=-contractions[i][j]
                    powers[i]-=n
                if not coefficient:break
            remaining=needed-sum_degree(uncontracted)
            if not coefficient or remaining<0:continue
            for (created,exps),ratio in creation_terms(late,remaining).items():
                delta=tuple(a-b-c for a,b,c in zip(requested,powers,exps,strict=True))
                n=scalar_coefficient(late,delta)
                if n:add_term(final,tuple(sorted(tuple(uncontracted)+created)),n*ratio*coefficient)
    return {k:v for k,v in final.items() if v}

def compare(a,b):
    v=subtract(a,b)
    # Redundant oscillator axes are allowed, so compare exact Hilbert norm.
    return inner(v,v)==0

rows={}
for alpha in (sp.Rational(1,4),sp.Rational(3,4)):
    drift=sp.Matrix([j+alpha-4 for j in range(8)])
    co=[K.from_sympy(x) for x in old["cl"](V.T*drift/8)]
    preparations=[old["series_response"](f,co,4) for f in src["polys"]]
    inputs=[old["input_word_map"](v) for v in preparations]
    packets=[old["vertex_zero"](v) for v in preparations]
    loaded=[src["vertex"](p,B,q) for p in packets]
    table={}
    for m in range(-2,3):
        early=[fused(v,m,False) for v in inputs]
        late=[fused(v,m,True) for v in inputs]
        ref_l=[src["vertex"](u,q,r,-m) for u in loaded]
        ref_e=[src["vertex"](src["vertex"](v,q,q,-m),B,2*q) for v in packets]
        for l,e,rl,re in zip(late,early,ref_l,ref_e,strict=True):
            assert compare(l,rl) and compare(e,re),(alpha,m)
        differences=[subtract(l,e) for l,e in zip(late,early,strict=True)]
        G=src["gram"](differences)
        table[str(m)]={"P1_early_raw_norm":str(inner(early[1],early[1])),
                       "P1_late_raw_norm":str(inner(late[1],late[1])),
                       "P1_contrast_raw_norm":str(G[1,1]),
                       "contrast_gram_rank":src["rank"](G)}
    rows[str(alpha)]=table
print(json.dumps({"fused_total_old_momentum":list(map(str,d)),
                  "mode_words_tested":100,
                  "no_half_integer_intermediate_created_by_fused_algorithm":True,
                  "cases":rows},indent=2))
