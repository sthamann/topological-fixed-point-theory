"""Actual 60-ray context Fourier transform and the first native/OS test.

Only the original source constructor is used. Dyadic entries and integer
Pauli matrices make all asserted equalities exact; no fitted tolerance.
"""
from pathlib import Path
import importlib.util, itertools as it, json, hashlib
import numpy as np
import sympy as sp

ROOT=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
HERE=Path(__file__).resolve().parent
source=ROOT/'experiments/theory-contracts/compiler-correlated-event-clock-20260919/source_channel.py'
spec=importlib.util.spec_from_file_location('original_source',source)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
zs=m.source_rays(); Z=np.column_stack(zs)
Es=[np.outer(z,z.conj())/4 for z in zs]
paulis=m.PAULIS
zero=m.V4.index((0,0,0,0))
chars=np.array([[np.trace(P@E).real for E in Es] for P in paulis])
assert set(chars.ravel())=={-1,0,1}
ctxs=sorted({tuple(v for v in range(16) if v!=zero and chars[v,r]!=0) for r in range(60)})
assert len(ctxs)==15 and all(len(c)==3 for c in ctxs)
groups=[[r for r in range(60) if tuple(v for v in range(16) if v!=zero and chars[v,r]!=0)==c] for c in ctxs]
assert all(len(g)==4 for g in groups)
for c,g in zip(ctxs,groups):
    assert np.array_equal(sum(Es[r] for r in g),np.eye(4))
    for v in c:
        assert np.array_equal(sum(chars[v,r]*Es[r] for r in g),paulis[v])
    for r in g:
        a,b=c[:2]
        assert np.array_equal((np.eye(4)+chars[a,r]*paulis[a])@(np.eye(4)+chars[b,r]*paulis[b])/4,Es[r])

# First15 rows: normalized context identity modes. Last45: signed incidences.
F=np.zeros((60,60),dtype=int)
for i,g in enumerate(groups): F[i,g]=1
inc=[]
for i,(c,g) in enumerate(zip(ctxs,groups)):
    for v in c:
        j=15+len(inc);F[j,g]=chars[v,g].astype(int);inc.append((i,v))
assert np.array_equal(F@F.T,4*np.eye(60))
Q=F/2
S=np.column_stack([E.ravel() for E in Es])
G=(S.conj().T@S).real
assert np.array_equal(G@(G-3*np.eye(60))@(G-15*np.eye(60)),np.zeros((60,60)))
assert np.trace(G)==60 and np.trace(G@G)==360
def exact_rank(A): return sp.Matrix(A).applyfunc(sp.nsimplify).rank()
assert exact_rank(G)==16
assert exact_rank(S@Q[:15].T)==1
assert exact_rank(S@Q[15:].T)==15

def key(A): return tuple(A.ravel())
lookup={key(E):i for i,E in enumerate(Es)}
reflections=[np.eye(4)-2*E for E in Es]
fixed=[]
for R in reflections:
    perm=[lookup[key(R@E@R.conj().T)] for E in Es]
    P=np.zeros((60,60),int);P[perm,np.arange(60)]=1
    transformed=Q@P@Q.T
    assert np.array_equal(transformed[:15,15:],np.zeros((15,45)))
    assert np.array_equal(transformed[15:,:15],np.zeros((45,15)))
    assert set(transformed[:15,:15].ravel()) <= {0,1}
    assert set(transformed[15:,15:].ravel()) <= {-1,0,1}
    assert np.array_equal(abs(transformed[15:,15:]).sum(0),np.ones(45))
    assert np.array_equal(abs(transformed[15:,15:]).sum(1),np.ones(45))
    fixed.append(int(np.trace(transformed[:15,:15])))
assert set(fixed)=={7}

avg=np.zeros((60,60))
signatures=[]
for W in paulis:
    perm=[lookup[key(W@E@W.conj().T)] for E in Es]
    P=np.zeros((60,60),int);P[perm,np.arange(60)]=1
    A=Q@P@Q.T
    assert np.array_equal(A[:15,:15],np.eye(15))
    assert np.array_equal(A[15:,15:],np.diag(np.diag(A[15:,15:])))
    signatures.append(tuple(np.diag(A[15:,15:])))
    avg+=A/16
assert len(set(signatures))==16
assert np.array_equal(avg,np.diag([1]*15+[0]*45))

# Actual native quartic5 already constructed in prior contract: recover its
# six coefficient vectors, then test the NEW context Fourier location.
_,raw,_=m.reflection_actions(zs)
outer,_=m.outer_mark_action(raw)
v=np.array([[3*int(p[q]!=q)-1 for q in range(6)] for p in outer])
assert np.array_equal(v.T@v,144*np.eye(6)-24*np.ones((6,6)))
assert np.array_equal(Q[15:]@v,np.zeros((45,6)))
assert exact_rank(Q[:15]@v)==5
assert np.array_equal(S@v,np.zeros((16,6)))
# Its actual fourth-power image survives, using the inherited native phase lift.
coherent4=(Z.conj().T@Z/4)**4
assert np.array_equal(coherent4@v,3*v)

# Exact input-attachment controls, not new proposed models.
X,Zp,I=m.X2,m.Z2,m.I2;Y=m.Y2
square=[[np.kron(Zp,I),np.kron(I,Zp),np.kron(Zp,Zp)],
        [np.kron(I,X),np.kron(X,I),np.kron(X,X)],
        [np.kron(Zp,X),np.kron(X,Zp),np.kron(Y,Y)]]
def prod(xs):
    out=np.eye(4,dtype=complex)
    for x in xs:out=out@x
    return out
signs=[int(np.trace(prod(row)).real/4) for row in square]
signs+=[int(np.trace(prod([square[i][j] for i in range(3)])).real/4) for j in range(3)]
assert signs==[1,1,1,1,1,-1]

duads=list(it.combinations(range(6),2)); synthemes=set()
for p in it.permutations(range(6)):
    synthemes.add(tuple(sorted(tuple(sorted(p[j:j+2])) for j in (0,2,4))))
synthemes=sorted(synthemes)
a,b=sp.symbols('a b',positive=True)
L=sp.zeros(30)
edges=[]
for j,syn in enumerate(synthemes):
    for d in syn:
        i=duads.index(d);w=a if 0 in d else b;edges.append((i,15+j,w))
        L[i,i]+=w;L[15+j,15+j]+=w;L[i,15+j]-=w;L[15+j,i]-=w
assert len(edges)==45
assert sp.trace(L)==30*a+60*b
assert sp.expand(sp.trace(L*L))==90*a*a+60*a*b+210*b*b
x=sp.zeros(30,1);x[duads.index((0,1))]=1;x-=sp.ones(30,1)/30
assert sp.simplify((x.T*L*x)[0]/(x.T*x)[0])==90*a/29

# Strictly positive, reflection invariant joint measure. Both separate
# classical marginals are reflection positive; the joint product is not.
# This tests the logical commuting-factor inference, not full TFPT CAR.
vv=sp.Matrix([1,-1,-1,1]);W=sp.ones(4)-vv*vv.T/2
assert all(z>0 for z in W)
assert (vv.T*W*vv)[0]/sum(W)==-sp.Rational(1,2)
for which in (0,1):
    marg=sp.zeros(2)
    for i,j in it.product(range(4),repeat=2):marg[(i//2 if which==0 else i%2),(j//2 if which==0 else j%2)]+=W[i,j]
    assert marg==sp.ones(2)*4

paths=[str(source.relative_to(ROOT)),
 'verification/v783_two_qubit_clifford.py',
 'experiments/theory-contracts/atlas-defect-synthesis-20261024/DOILY-QUARTIC.txt',
 'experiments/theory-contracts/compiler-quartic-carry-20260919/PROOF.txt',
 'experiments/theory-contracts/joint-source-clock-field-closure-20261102/PROOF.txt']
result=dict(original_rays=60,contexts=15,incidences=45,
    Fourier_exact=True,one_copy_Gram_spectrum={'15':1,'3':15,'0':44},
    identity_sector_synthesis_rank=1,signed_sector_synthesis_rank=15,
    all_60_reflections_context_fixed_count=7,Pauli_signed_signatures=16,
    Pauli_fixed_formal_dimensions={'context':15,'signed':0},
    native_quartic_five={'context_rank':5,'signed_rank':0,'one_copy_synthesis_rank':0,'coherent_fourth_Gram_eigenvalue':3},
    Mermin_products=signs,Doily_trace=str(sp.trace(L)),Doily_trace_square=str(sp.expand(sp.trace(L*L))),
    Doily_Rayleigh_bound='90*a/29',positive_joint_counterexample={'separate_reflection_Grams':'4*J2','joint_product_reflected_expectation':'-1/2','scope':'Commuting marginal positivity only; not a full TFPT CAR countermodel'},
    source_pins={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
(HERE/'results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(result,ensure_ascii=False,indent=2))
