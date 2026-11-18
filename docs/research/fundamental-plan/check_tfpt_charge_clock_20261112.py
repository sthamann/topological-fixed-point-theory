"""Exact Wick expansion on the original E8 root shell, no state fit."""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import factorial, prod
from pathlib import Path
import ast
import sys

repo = Path(sys.argv[1])
original = repo / 'docs/research/fundamental-plan/check_tfpt_five_clock_source_20261112.py'
tree = ast.parse(original.read_text())
# Reuse the unmodified original calculation up to and including both root and
# polynomial data. The orbit checks only need those data, not a chosen lift.
end = next(i for i,node in enumerate(tree.body) if isinstance(node,ast.FunctionDef) and node.name == 'apply')
data = {}
exec(compile(ast.Module(body=tree.body[:end],type_ignores=[]),str(original),'exec'),data)
roots = sorted(data['roots'],key=lambda r: tuple((z.real,z.imag) for z in r))
assert len(roots) == 240
a,b = data['a'],data['b']
q,v = 4,b[4]
W = data['W']

def monomials(j):
    terms = {tuple(4*(k==i) for k in range(4)): W[0][j] for i in range(4)}
    for n,pairs in enumerate((((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2))),1):
        for pair in pairs:
            terms[tuple(2*(k in pair) for k in range(4))] = 6*W[n][j]
    terms[(1,1,1,1)] = 24*W[4][j]
    return {e:c for e,c in terms.items() if c}

P,Q = monomials(q),monomials(v)
indices = [g for g in product(range(5),repeat=4) if sum(g)<=4]
fact = lambda ex: prod(factorial(n) for n in ex)

def evaluate(poly,r):
    return sum(c*prod(x**n for x,n in zip(r,e)) for e,c in poly.items())

def derivative(poly,g,r):
    return sum(c*prod(factorial(e0)//factorial(e0-g0)*x**(e0-g0)
                      for x,e0,g0 in zip(r,e,g))
               for e,c in poly.items() if all(g0<=e0 for e0,g0 in zip(e,g)))

def pair_add(target,k,z,denom):
    assert z.real == int(z.real) and z.imag == int(z.imag)
    target[k][0] += Fraction(int(z.real),denom)
    target[k][1] += Fraction(int(z.imag),denom)

def paired(poly1,poly2,r):
    out = [[Fraction(0),Fraction(0)] for _ in range(5)]
    for g in indices:
        k=sum(g)
        value=derivative(poly1,g,r).conjugate()*derivative(poly2,g,r)
        pair_add(out,k,value,fact(g)*2**(2*(4-k)))
    return tuple(tuple(p) for p in out)

def separate_count(poly1,poly2,r):
    # Explicit partial bijections between Wick legs, grouped only by the
    # chosen subsets and equal-colour multiplicities, not by derivatives.
    out = [[Fraction(0),Fraction(0)] for _ in range(5)]
    for e,c in poly1.items():
        L=[i for i,n in enumerate(e) for _ in range(n)]
        for f,d in poly2.items():
            R=[i for i,n in enumerate(f) for _ in range(n)]
            for im in range(16):
                I=[n for j,n in enumerate(L) if (im>>j)&1]
                K=Counter(I); k=len(I)
                for jm in range(16):
                    J=[n for j,n in enumerate(R) if (jm>>j)&1]
                    if len(J)!=k or Counter(J)!=K: continue
                    z=c*d*prod(factorial(n) for n in K.values())
                    z*=prod(r[n].conjugate() for j,n in enumerate(L) if not ((im>>j)&1))
                    z*=prod(r[n] for j,n in enumerate(R) if not ((jm>>j)&1))
                    pair_add(out,k,z,2**(2*(4-k)))
    return tuple(tuple(p) for p in out)

def contrast(A,B):
    return tuple(tuple(b+a/5 for a,b in zip(z,w)) for z,w in zip(A,B))

def value(poly,x,connected=False):
    return tuple(sum((poly[k][c]*x**(4-k) for k in range(int(connected),5)),Fraction(0)) for c in range(2))

S=Counter()
tot=[(Fraction(0),Fraction(0)) for _ in range(5)]
first={}
for r in roots:
    A,B = paired(P,P,r),paired(P,Q,r)
    C=contrast(A,B)
    assert A[4]==(960,0) and B[4]==(-192,0)
    assert evaluate(P,r)==data['polynomial_values'](r)[q]
    assert evaluate(Q,r)==data['polynomial_values'](r)[v]
    S[C]+=1
    first.setdefault(C,(r,A,B))
    for k in range(5):
        tot[k]=tuple(a0+c0 for a0,c0 in zip(tot[k],C[k]))
assert not any(any(p) for p in tot)
assert all(paired(P,P,r)==separate_count(P,P,r) and
           paired(P,Q,r)==separate_count(P,Q,r)
           for r,A,B in first.values())
print('Fixed labels', q, v, 'on original five-clock data; root shell size',len(roots))
print('Distinct contrast polynomials (C0,...,C4):',len(S))
for C,n in sorted(S.items()):
    r,A,B=first[C]
    print('count',n,'root',r,'contrast',C,'A',A,'B',B)
    print('connected F_A,F_B,contrast at X=1/4 and 1/16:',
          [(x,value(A,x,True),value(B,x,True),value(C,x,True)) for x in (Fraction(1,4),Fraction(1,16))])
print('Unresolved full shell cancels the contrast coefficient by coefficient.')

from collections import defaultdict

def same_gaussian_class(r,s):
    # Construction A: (r-s)/(1+i) lies in L iff integral and its odd support
    # belongs to the original binary RM(1,3) code.
    diff=[(u-v)/(1+1j) for u,v in zip(r,s)]
    as_real=[x for y in diff for x in (y.real,y.imag)]
    if not all(x==int(x) for x in as_real): return False
    odd=tuple(i for i,x in enumerate(as_real) if int(x)%2)
    return not odd or odd==tuple(range(8)) or odd in data['supports']

by_pair=defaultdict(list)
for r in roots:
    mags=[abs(z)**2 for z in data['polynomial_values'](r)]
    assert sorted(mags)==[256]*4+[1024]*2
    pair=tuple(i for i in range(6) if mags[i]==1024)
    by_pair[pair].append(r)

assert len(by_pair)==15
for pair, rs in by_pair.items():
    assert len(rs)==16
    assert all(same_gaussian_class(rs[0],r)==(r in rs) for r in roots)
    values={contrast(paired(P,P,r),paired(P,Q,r)) for r in rs}
    assert len(values)==1
    C=next(iter(values))
    print('Native quotient pair',pair, 'C2',C[2][0],'C1',C[1][0])

import sympy as sy
i=sy.I
z,az,bz=1,-1,-i
X=lambda w: sy.simplify((z-w)**2*(az-bz)**2/((z-az)*(z-bz)*(w-az)*(w-bz)))
positions=(i,(3+4*i)/5)
coords=tuple(X(w) for w in positions)
print('With the other two original quarter marks as charge endpoints, moved w and X:',list(zip(positions,coords)))
for pair in ((0,4),(0,1),(0,2),(4,5)):
    r=by_pair[pair][0]
    A,B=paired(P,P,r),paired(P,Q,r)
    C=contrast(A,B)
    print('Moved point with charges on pair',pair,[(x,value(B,Fraction(x),True)[0]/value(A,Fraction(x),True)[0],
          value(C,Fraction(x),True)[0]/960) for x in coords])

# Derive the full unconditioned *reflected* selector algebra using the actual
# root class partition. For each of 15 nonzero vectors its 16 roots agree
# on all response coefficients, before any choice of state weights.
assert all(contrast(paired(P,P,rs[0]),paired(P,Q,rs[0]))[3]==(0,0)
           for rs in by_pair.values())

# No assumption about which compatible quarter lift is the physical one.
# Read the full original enumerator without changing its lattice data.
import contextlib, io
full = {'__name__':'__tfpt_original_check__'}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(tree,str(original),'exec'),full)
apply,mp,mc=full['apply'],full['mat_power'],full['mat_conj']
compatible=[]
for k in full['kernel']:
    g=full['matmul'](k,full['qa'])
    if any(mc(g,u)==mp(u,2) for u in full['as_five']) and any(mc(g,u)==mp(u,2) for u in full['bs_five']):
        compatible.append(g)
assert len(compatible)==32
G=full['quarter']
fixed_pair=[pair for pair in by_pair if tuple(sorted(G[j] for j in pair))==pair]
assert fixed_pair==[(4,5)]
for g in compatible:
    for pair,rs in by_pair.items():
        out=tuple(sorted(G[j] for j in pair))
        assert {apply(g,r) for r in rs}==set(by_pair[out])
fix_by_g = Counter(sum(apply(g,r)==r for r in by_pair[(4,5)]) for g in compatible)
print('Only G-fixed nonzero class:', fixed_pair[0])
print('Compatible G lifts grouped by number of pointwise fixed roots in that class:',dict(sorted(fix_by_g.items())))
assert not any(apply(g,r)==r for g in compatible for r in roots if r not in by_pair[(4,5)])
# If the action on charges is a finite permutation with trivial phases,
# invariant orbit sums include each nontrivial orbit; this is different
# from requiring the charge label itself to be fixed.

# Conformal shell coefficient check in the unscaled integral coordinates.
n=4
root_list=[]
choices=[(0,)]*8
for support in [(), *data['supports'], tuple(range(8))]:
    choices=[tuple(x for x in range(-4,5) if x%2==(i in support)) for i in range(8)]
    for r8 in product(*choices):
        sq=sum(x*x for x in r8)
        if not 0<sq<=4*n: continue
        assert sq%4==0
        root_list.append((sq//4,tuple(complex(x,y) for x,y in zip(r8[::2],r8[1::2]))))
theta=Counter()
degree=Counter()
for h,r in root_list:
    if same_gaussian_class(r,by_pair[(4,5)][0]):
        theta[h]+=1
        for k,c in enumerate(contrast(paired(P,P,r),paired(P,Q,r))):
            degree[(h,k)]+=c[0]
assert dict(sorted(theta.items()))=={1:16,2:128,3:448,4:1024}
print('Canonical fixed-class theta coefficients h=1..4:',dict(sorted(theta.items())))
print('Fixed-class summed contrast coefficients by h:',[[degree[h,k] for k in range(5)] for h in range(1,5)])

u=sy.symbols('u')
fixed_dims=Counter()
fixed_cycles=Counter()
for g in compatible:
    m=sy.Matrix(g).applyfunc(lambda x:sy.nsimplify(x))
    fixed_dims[4-(m-sy.eye(4)).rank()]+=1
    seen=set()
    cycles=[]
    for root in by_pair[(4,5)]:
        if root in seen: continue
        orbit=set()
        x=root
        while x not in orbit:
            orbit.add(x)
            x=apply(g,x)
        seen|=orbit
        cycles.append(len(orbit))
    fixed_cycles[tuple(sorted(Counter(cycles).items()))]+=1
print('Fixed complex dimensions of G lifts:',dict(sorted(fixed_dims.items())))
print('G cycle types on 16 roots of class (4,5):',dict(sorted(fixed_cycles.items())))

from decimal import Decimal as D, getcontext
getcontext().prec=45
q0=D(2)**-16
dz=lambda x:D(x.numerator)/D(x.denominator)
osc=prod((1-q0**j)**-8 for j in range(1,13))
approx_z=osc*sum(D(theta[h])*q0**h for h in range(1,5))
d_ratio= D('0.005') # Ratio bound for all j>=5 of the polynomial q-tail.
ratio = (D(6)/5)**4*(D(7)/6)**4*q0
assert ratio < d_ratio
for eu in (D(1)/16,D(1)/256):
    K=sum(D(j)*(eu**j+(q0/eu)**j)/(1-q0**j) for j in range(1,33))
    # Closed geometric-tail bound, K < 1 also including the omitted terms.
    m=max(eu,q0/eu)
    kt=2*m**33*(D(33)-D(32)*m)/((1-m)**2*(1-q0))
    assert K+kt<1
    mark=sum(q0**h * sum(dz(degree[h,k])*K**k for k in range(5)) for h in range(1,5))*osc
    tail = D(1152*240*5**4*6**4)*q0**5/(1-ratio)
    tail += D(4608*240)*kt*sum(q0**h*D(h**4*(h+1)**4) for h in range(1,5))
    tail *= osc
    # Control oscillator omissions by bounding their positive log-tail.
    eosclog=8*q0**13/((1-q0)*(1-q0**13))
    tail += 2*eosclog*(abs(mark)+tail)
    assert tail<D('2e-13')
    assert mark+tail < 0
    print('Absolute conformal trace q=2^-16, exp(-t)=',eu,'Z45=',approx_z,
          'marked contrast=',mark, 'error <',tail)

# The Gaussian quotient has its canonical alternating form.
# Coordinates are unscaled integral r, so the real lattice inner product
# is half the real part of the complex Hermitian product.
def omega(r,s):
    prod2=sum(((1+1j)*x).conjugate()*y for x,y in zip(r,s)).real
    assert prod2%2==0
    return int(prod2/2)%2

rep=by_pair[(4,5)][0]
assert all(len({omega(r,s) for r in R for s in S})==1
           for R,S in product(by_pair.values(),repeat=2))
assert all(omega(R[0],R[0])==0 and any(omega(R[0],s) for s in roots)
           for R in by_pair.values())
parity={p:omega(rep,rs[0]) for p,rs in by_pair.items()}
assert Counter(parity.values())=={0:7,1:8}
for g in compatible:
    assert all(omega(rep,r)==omega(rep,apply(g,r)) for r in roots)
a_lifts=full['as_five']
assert all(any(omega(rep,r)!=omega(rep,apply(ga,r)) for r in roots) for ga in a_lifts)
for par in (0,1):
    rs=[r for r in roots if omega(rep,r)==par]
    result=tuple(sum(contrast(paired(P,P,r),paired(P,Q,r))[k][0] for r in rs)
                 for k in range(5))
    print('Quarter-fixed parity',par,'native root count',len(rs),
          'pairs',[p for p,c in sorted(parity.items()) if c==par],
          'sum contrast coefficients',result)
    if par==0:
        # Connected even crystallographic root system of rank eight, with
        # 112 roots: D8. Irreducibility is checked by nonorthogonality.
        vecs=[tuple(int(v) for x in r for v in (x.real,x.imag)) for r in rs]
        assert sy.Matrix(vecs).rank()==8
        seen={vecs[0]}
        queue=list(seen)
        for v0 in queue:
            for w in vecs:
                if w not in seen and sum(a0*b0 for a0,b0 in zip(v0,w)):
                    seen.add(w)
                    queue.append(w)
        assert len(seen)==112

# Compare the actual older order-eight character. Every isometric D8
# dictionary is a signed permutation of its eight vector weight directions.
from itertools import combinations
even_vecs=[tuple(int(y) for z0 in r for y in (z0.real,z0.imag)) for r in roots if omega(rep,r)==0]
dot=lambda a0,b0:sum(i*j for i,j in zip(a0,b0))
vweights=set()
for r,s in combinations(even_vecs,2):
    if dot(r,s): continue
    p=tuple(x+y for x,y in zip(r,s))
    if p in vweights: continue
    if all(dot(p,k)%4==0 for k in even_vecs):
        vweights.add(p)
assert len(vweights)==16
axes=sorted(p for p in vweights if p>tuple(-x for x in p))
assert len(axes)==8 and all(dot(u0,v0)==8*(j==k) for j,u0 in enumerate(axes) for k,v0 in enumerate(axes))
def cart8(r):
    # Twice the standard D8 coordinates in the signed vector weight frame.
    vec=tuple(int(y) for z0 in r for y in (z0.real,z0.imag))
    v=[dot(p,vec) for p in axes]
    assert all(x%2==0 for x in v)
    return tuple(x//2 for x in v)
if sum(cart8(next(r for r in roots if omega(rep,r)==1)))%4:
    axes[0]=tuple(-x for x in axes[0])
standard_d8={tuple(2*(s*(k==i)+t*(k==j)) for k in range(8)) for i,j in combinations(range(8),2)
             for s,t in product((-1,1),repeat=2)}
assert {cart8(r) for r in roots if omega(rep,r)==0}==standard_d8
assert {cart8(r) for r in roots if omega(rep,r)==1}=={v for v in product((-1,1),repeat=8) if sum(v)%4==0}
g_perm=set()
for g in compatible:
    cols=[]
    for axis in axes:
        # g applied to a twice-raw vector, so no fractions enter.
        cv=apply(g,tuple(complex(x,y) for x,y in zip(axis[::2],axis[1::2])))
        y=tuple(int(i) for z0 in cv for i in (z0.real,z0.imag))
        col=tuple(dot(aj,y)//8 for aj in axes)
        assert sum(abs(c) for c in col)==1
        cols.append(col)
    g_perm.add(tuple(next((j,v0) for j,v0 in enumerate(row) if v0) for row in zip(*cols)))
ws=[]
for hpos in combinations(range(8),2):
    for lpos in combinations([j for j in range(8) if j not in hpos],3):
        for sig in product((-1,1),repeat=5):
            wi=[0]*8
            for j,s in zip(hpos+lpos,sig):
                wi[j]=s*(4 if j in hpos else 2)
            ws.append(tuple(wi))
assert len(set(ws))==len(ws)==17920
def residue_check(gv,w,k):
    differences=[s*w[j]-k*x for (j,s),x in zip(gv,w)]
    remainder=differences[0]%8
    return remainder in (0,4) and all(x%8==remainder for x in differences) and sum(differences)%16==0
assert all(sum(x*x for x in w)==44 and sum(w)%2==0 for w in ws)
norm_counts=Counter()
w_examples={}
for gp in sorted(g_perm):
    result=[]
    for k in (1,3,5,7):
        ok=[w for w in ws if residue_check(gp,w,k)]
        result.append(len(ok))
        if ok:
            w_examples.setdefault(k,(gp,ok[0]))
    norm_counts[tuple(result)]+=1
print('Older order-eight character on the unique G-invariant D8, counts for G w = k w modulo 8L, k=1,3,5,7:',dict(norm_counts))
print('Representative exact old-character match examples:',w_examples)

import numpy as np
weights=np.asarray(ws,dtype=np.int16)
charge_roots=np.asarray([cart8(r) for r in roots],dtype=np.int16)
nominal=(weights@charge_roots.T//2)%8
odd_idx=np.array([omega(rep,r)==1 for r in roots])
same_counts=Counter()
orbit_hist=Counter()
sew_max=set()
word_examples={}
for gp in sorted(g_perm):
    p=np.asarray([a0 for a0,_ in gp])
    signs=np.asarray([b0 for _,b0 in gp])
    wp=weights.copy()
    residues=[nominal]
    for j in range(1,5):
        wp=wp[:,p]*signs
        if j<4:
            residues.append((wp@charge_roots.T//2)%8)
    assert np.array_equal(wp,weights)
    # Equality on all roots is exactly equality of characters on L.
    squared=np.asarray([np.any(np.all(residues[2]==(k*nominal)%8,axis=1)) for k in (1,3,5,7)])
    same_counts[tuple(bool(x) for x in squared)]+=1
    orbit_two=np.any([np.all(residues[2]==(k*nominal)%8,axis=1) for k in (1,3,5,7)],axis=0)
    orbit_hist[int(orbit_two.sum())]+=1
    common=np.logical_and.reduce([rr==nominal for rr in residues[1:]])
    co=common[:,odd_idx]
    noco=co.sum(axis=1)
    sew_max.add(int(noco.max()))
    for count in sorted(set(noco.tolist())):
        if count in word_examples: continue
        n0=int(np.flatnonzero(noco==count)[0])
        word_examples[count]=(tuple(int(x) for x in weights[n0]),tuple(sorted(set(int(i) for i in nominal[n0,odd_idx][co[n0]]))))
assert same_counts and orbit_hist and sew_max
print('Half-turn normalizing possibilities k=1,3,5,7, by compatible quarter:',dict(same_counts))
print('Number of old frames with two-step cyclic character orbit, by quarter:',dict(orbit_hist))
print('Odd norm-two roots retaining exactly their eighth-tick residue through all four steps:',sorted(word_examples),
      '; best count per quarter:',sorted(sew_max))
print('Example old frame and preserved odd residues for each survivor count:',word_examples)

for p,m in full['transpositions'].items():
    for pair,rs in by_pair.items():
        out=tuple(sorted(p[j] for j in pair))
        assert {apply(m,r) for r in rs}==set(by_pair[out])

par_deg=Counter()
par_theta=Counter()
for h,r in root_list:
    par=omega(rep,r)
    par_theta[par,h]+=1
    for k,c in enumerate(contrast(paired(P,P,r),paired(P,Q,r))):
        par_deg[par,h,k]+=c[0]
assert [par_theta[1,h] for h in range(1,5)]==[128,1024,3584,8192]
assert all(par_deg[0,h,k]==-par_deg[1,h,k] for h,k in product(range(1,5),range(5)))
odd_z=osc*sum(D(par_theta[1,h])*q0**h for h in range(1,5))
for eu in (D(1)/16,D(1)/256):
    K=sum(D(j)*(eu**j+(q0/eu)**j)/(1-q0**j) for j in range(1,33))
    mark=sum(q0**h * sum(dz(par_deg[1,h,k])*K**k for k in range(5)) for h in range(1,5))*osc
    assert mark>tail
    print('Absolute odd conformal trace q=2^-16, exp(-t)=',eu,'Z_odd=',odd_z,
          'marked contrast=',mark, 'error <',tail)

# Exhaust all binary characters on L, including those that do not factor
# through L/(1+i)L, to test uniqueness of the D8 parity at the full level.
support_sets=[(),*data['supports'],tuple(range(8))]
bits=[tuple(int(i in s) for i in range(8)) for s in support_sets]
classes2=set()
for b0 in bits:
    for c0 in product((0,2),repeat=8):
        v=tuple(x+y for x,y in zip(b0,c0))
        classes2.add(min(tuple((x+2*t)%4 for x,t in zip(v,d0)) for d0 in bits))
assert len(classes2)==256
R=np.asarray([tuple(int(y) for x in r for y in (x.real,x.imag)) for r in roots],dtype=np.int16)
characters=(np.asarray(sorted(classes2),dtype=np.int16)@R.T//2)%2
assert len(set(tuple(x) for x in characters))==256
types=Counter()
indices_by_r={r:i for i,r in enumerate(roots)}
for g in compatible:
    permutation=[indices_by_r[apply(g,r)] for r in roots]
    eq=characters[np.all(characters==characters[:,permutation],axis=1)]
    types[tuple(sorted(Counter(np.sum(row==0) for row in eq).items()))]+=1
    even112=eq[np.sum(eq==0,axis=1)==112]
    assert len(even112)>=1 and any(np.array_equal(row,odd_idx) for row in even112)
print('All 256 binary E8 characters, fixed types (#zero roots, multiplicity):',dict(types))

# Standard charged lattice lift, not just its action on Cartan currents.
holonomy=Counter()
for g in compatible:
    square=mp(g,2)
    par4=[]
    for r in roots:
        gr=apply(square,r)
        pair=sum(x.conjugate()*y for x,y in zip(r,gr)).real
        assert pair%2==0
        par4.append(int(pair/2)%2)
    assert any(np.array_equal(par4,c) for c in characters)
    holonomy[(sum(par4),bool(np.array_equal(par4,odd_idx)))]+=1
print('Standard lattice lift fourth-power character (#negative roots, equal to canonical parity):',dict(holonomy))

def chart_for_binary(par0):
    even=R[~par0]
    mini=set()
    for r,s in combinations(map(tuple,even.tolist()),2):
        if dot(r,s): continue
        vec=tuple(x+y for x,y in zip(r,s))
        if vec not in mini and np.all((even@np.asarray(vec))%4==0):
            mini.add(vec)
    assert len(mini)==16
    ax=np.asarray(sorted(v for v in mini if v>tuple(-x for x in v)),dtype=np.int16)
    assert np.array_equal(ax@ax.T,8*np.eye(8,dtype=np.int16))
    if (R[par0][0]@ax.T//2).sum()%4:
        ax[0]*=-1
    assert set(map(tuple,(even@ax.T//2).tolist()))==standard_d8
    assert all(np.all(np.abs(c)==1) and c.sum()%4==0 for c in R[par0]@ax.T//2)
    return ax

def get_perm(g, ax):
    rows=[]
    for x in ax:
        cv=apply(g,tuple(complex(y,z) for y,z in zip(x[::2],x[1::2])))
        yy=[int(y) for z0 in cv for y in (z0.real,z0.imag)]
        rows.append(ax@yy//8)
    rows=np.asarray(rows).T
    assert np.all(np.abs(rows).sum(axis=0)==1) and np.all(np.abs(rows).sum(axis=1)==1)
    p=np.argmax(np.abs(rows),axis=1)
    return p,rows[np.arange(8),p]

whole_compare=Counter()
alt_part=Counter()
charts={}
for g in compatible:
    perm=[indices_by_r[apply(g,r)] for r in roots]
    eq=characters[np.all(characters==characters[:,perm],axis=1)]
    for ch in eq[np.sum(eq==0,axis=1)==112]:
        key=tuple(ch)
        if key not in charts:
            charts[key]=chart_for_binary(ch.astype(bool))
        p,s=get_perm(g,charts[key])
        k0=[]
        for wg in (weights[:,p]*s,weights[:,p][:,p]*s[p]*s):
            counts=[]
            for k in (1,3,5,7):
                diff=wg-k*weights
                rem=diff%8
                cond=np.all(rem==rem[:,[0]],axis=1) & ((rem[:,0]%4)==0) & ((diff.sum(axis=1)%16)==0)
                counts.append(int(cond.sum()))
            k0.append(tuple(counts))
        whole_compare[tuple(k0)]+=1
        if key!=tuple(odd_idx):
            alt_part[tuple(k0)]+=1
print('All invariant D8 charts, quarter and half-turn normalization counts:',dict(whole_compare))
print('Same comparison excluding the canonical Gaussian parity:',dict(alt_part))

assert holonomy=={(128,True):32}
assert whole_compare=={((0,0,0,0),(0,0,0,0)):32,((0,0,0,0),(128,128,128,128)):64}
assert alt_part=={((0,0,0,0),(128,128,128,128)):64}
shell=Counter(h for h,r in root_list)
Ztw=osc*(1+sum(D(shell[h]-2*par_theta[1,h])*q0**h for h in range(1,5)))
assert D('0.9998')<Ztw<D(1)
for eu in (D(1)/16,D(1)/256):
    K=sum(D(j)*(eu**j+(q0/eu)**j)/(1-q0**j) for j in range(1,33))
    Dtw=-2*sum(q0**h * sum(dz(par_deg[1,h,k])*K**k for k in range(5)) for h in range(1,5))*osc
    assert Dtw<-2*tail
    print('Full standard-lift closed twist q=2^-16, exp(-t)=',eu,
          'Z_tw=',Ztw,'contrast_tw=',Dtw,'error(contrast) <',2*tail)
print('G transport with covariant residue projectors preserves all odd root labels; fourth power contributes their exact minus sign.')
