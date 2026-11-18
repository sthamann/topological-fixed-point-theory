from pathlib import Path
import itertools as it, importlib.util, json
import numpy as np
import sympy as s

ROOT=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
HERE=Path(__file__).resolve().parent
path=ROOT/'experiments/theory-contracts/compiler-correlated-event-clock-20260919/source_channel.py'
spec=importlib.util.spec_from_file_location('native',path)
n=importlib.util.module_from_spec(spec);spec.loader.exec_module(n)
duads=list(it.combinations(range(5),2))
A=np.array([[int(set(a).isdisjoint(b)) for b in duads] for a in duads],int)
edges=[(i,j) for i in range(10) for j in range(i+1,10) if A[i,j]]
def matching(left):
    if not left:
        yield ();return
    i=left[0]
    for j in left[1:]:
        if A[i,j]:
            for tail in matching([x for x in left if x not in (i,j)]):
                yield tuple(sorted(((i,j),)+tail))
ms=list(matching(list(range(10))))
assert len(ms)==6
perms=list(it.permutations(range(5)))
def wedge(g):
    U=np.zeros((10,10),int)
    for i,(a,b) in enumerate(duads):
        c,d=g[a],g[b];U[duads.index(tuple(sorted((c,d)))),i]=1 if c<d else -1
    return U
Us=[wedge(g) for g in perms]
def act_matching(m,U):
    p=np.argmax(abs(U),axis=0)
    return tuple(sorted(tuple(sorted((int(p[a]),int(p[b])))) for a,b in m))
Rs=[];solutions=[]
for m in ms:
    stab=[U for U in Us if act_matching(m,U)==m]
    good=[]
    for signs in it.product([-1,1],repeat=5):
        R=np.zeros((10,10),int)
        for (a,b),z in zip(m,signs):R[a,b]=R[b,a]=z
        if all(np.array_equal(U@R@U.T,R) or np.array_equal(U@R@U.T,-R) for U in stab):good.append(R)
    solutions.append(len(good));assert len(good)==2
    Rs.append(good[0])
for U in Us:
    for m,R in zip(ms,Rs):
        target=Rs[ms.index(act_matching(m,U))]
        assert np.array_equal(U@R@U.T,target) or np.array_equal(U@R@U.T,-target)
def cycles(p):
    remaining=set(range(len(p)));cs=[]
    while remaining:
        a=min(remaining);b=a;k=0
        while b in remaining:remaining.remove(b);k+=1;b=p[b]
        cs.append(k)
    return sorted(cs)
cycle_results={}
for g in [(1,2,0,3,4),(1,0,2,3,4)]:
    U=wedge(g);cycle_results[str(g)]=cycles([ms.index(act_matching(m,U)) for m in ms])
S6=sum(np.kron(R,R) for R in Rs)
G6=s.Matrix([[int(np.trace(R.T@T)) for T in Rs] for R in Rs])
assert int(np.trace(S6@S6))==720

# Same marked duad chart of the ORIGINAL sixty matrices.
bell,raw,_=n.reflection_actions(n.source_rays())
outer,_=n.outer_mark_action(raw)
mapping=n.find_bell_partition_intertwiner(bell,outer)
q=0;others=list(range(1,6));labels=[]
for b in range(10):
    side=set(mapping[b]);side=side if q in side else set(range(6))-side
    labels.append(tuple(others.index(x) for x in sorted(side-{q})))
order=[labels.index(d) for d in duads]
bell=[U[np.ix_(order,order)] for U in bell]
selected=[k for k,p in enumerate(outer) if p[q]!=q]
S=np.zeros((100,100),complex);J=S.copy()
for k in selected:
    U=bell[k];D=np.diag(np.diag(U)!=0).astype(int)
    Ks=U@D;Kj=U@(np.eye(10)-D)
    S+=np.kron(Ks.conj(),Ks);J+=np.kron(Kj.conj(),Kj)
assert np.array_equal(S,np.diag(np.diag(S)))
T=J@np.diag(1/(20-np.diag(S)))
assert np.max(abs(T.imag))==0
T=s.Matrix(T.real).applyfunc(s.nsimplify)
assert s.trace(T*T)==s.Rational(68,15)
def vec(K):return s.Matrix(np.array(K).reshape(-1,order='F'))
def choi(S):
    return s.Matrix(100,100,lambda a,b:S[(a%10)+10*(b%10),(a//10)+10*(b//10)])
C=choi(T)
# The actual G240/A5 action must be intertwined, not identified by labels.
D=np.diag([1,-1,1j,1j,1,-1j,-1j,-1j,1j,1j])
def comp(p,q):return tuple(p[q[a]] for a in range(6))
pg=comp(outer[0],outer[46]);ph=comp(outer[1],outer[21])
g5=tuple(pg[a]-1 for a in range(1,6));h5=tuple(ph[a]-1 for a in range(1,6))
assert np.array_equal(D@(-bell[0]@bell[46])@D.conj().T,wedge(g5))
assert np.array_equal(D@(bell[1]@bell[21])@D.conj().T,wedge(h5))
W=s.Matrix(np.kron(D.conj(),D)).applyfunc(s.nsimplify)
C=W*C*W.conjugate().T
active=[a+10*b for a in range(10) for b in range(10) if A[a,b]]
Cedge=C.extract(active,active)
assert Cedge.eigenvals()=={s.Rational(2,15):15,s.Rational(8,15):15}
C6=sum((vec(R)*vec(R).T for R in Rs),s.zeros(100))/6

# Existing five-complement source in the same declared ordered wedge chart.
def signperm(p):return (-1)**sum(p[i]>p[j] for i in range(5) for j in range(i+1,5))
Cs=[]
for m in range(5):
    H=np.zeros((10,10),int)
    for a,(i,j) in enumerate(duads):
        for b,(k,l) in enumerate(duads):
            if len({m,i,j,k,l})==5:H[b,a]=signperm((m,i,j,k,l))
    Cs.append(H)
C5=sum((vec(H)*vec(H).T for H in Cs),s.zeros(100))/3
assert sum((H.T@H for H in Cs),np.zeros((10,10),int)).tolist()==(3*np.eye(10,dtype=int)).tolist()

# Exact successful-branch bounds inside the native STOPPED channel.
# Generalized eigenvalues computed on the small Kraus space.
inv=Cedge.inv()
embeddings={}
for label,kraus,den in [('six',Rs,6),('five',Cs,3)]:
    V=s.Matrix.hstack(*(vec(K).extract(active,[0]) for K in kraus))
    E=V.T*inv*V/den
    ev=E.eigenvals()
    maxeval=max(ev,key=lambda x:float(x))
    t=s.simplify(1/maxeval)
    # Congruence reduces positivity and rank to these small Gram eigenvalues.
    assert all(1-t*lam>=0 for lam in ev)
    embeddings[label]={'generalized_eigenvalues':{str(k):v for k,v in ev.items()},'max_branch_probability':str(t),'residual_rank':30-ev[maxeval]}
assert embeddings['six']['max_branch_probability']=='1/25'
assert embeddings['five']['max_branch_probability']=='1/15'
# Exact common-sector decomposition proves the sharp joint region.
wv=sum((vec(H) for H in Cs),s.zeros(100,1))
P1=wv*wv.T/30
P4=C5/2-P1
P5=3*C6/4-5*P1/2
for P,r in [(P1,1),(P4,4),(P5,5)]:assert P*P==P and s.trace(P)==r
assert P1*P4==P1*P5==P4*P5==s.zeros(100)
Psym=s.zeros(100);Panti=s.zeros(100)
for a,b in edges:
    Sedge=s.zeros(10);Sedge[a,b]=Sedge[b,a]=1
    Aedge=s.zeros(10);Aedge[a,b]=s.I;Aedge[b,a]=-s.I
    sv=s.Matrix(list(Sedge.T));av=s.Matrix(list(Aedge.T))
    Psym+=sv*sv.conjugate().T/2;Panti+=av*av.conjugate().T/2
assert C==s.Rational(2,15)*Psym+s.Rational(8,15)*Panti
Psi=s.Matrix(J.real).applyfunc(lambda v:s.Integer(int(v)))/12
assert W*choi(Psi)*W.conjugate().T==s.Rational(2,3)*Panti
assert C-s.Rational(3,5)*(s.Rational(2,3)*Panti)==s.Rational(2,15)*(Psym+Panti)
assert Psym*(P1+P4+P5)==P1+P4+P5
assert C6==s.Rational(10,3)*P1+s.Rational(4,3)*P5
assert C5==2*(P1+P4)
# Both maximal residuals are PSD by the small Gram calculation above.
combined=C-s.Rational(1,50)*C6-s.Rational(1,30)*C5
assert combined==(C-C6/25)/2+(C-C5/15)/2
partial_trace=s.Matrix(10,10,lambda a,b:sum(combined[j+10*a,j+10*b] for j in range(10)))
assert partial_trace==s.Rational(71,75)*s.eye(10)
# Exact fixed-space check of the new six-matching channel.
assert n.rank_mod_p(S6-6*np.eye(100,dtype=int),3,1)==99
result={'six_matching_count':len(ms),'sign_solutions_per_matching':solutions,'matching_cycle_types':cycle_results,'six_kraus_gram_spectrum':{str(k):v for k,v in G6.eigenvals().items()},'six_channel_fixed_space_dimension':1,'six_channel_trace_square':'20','native_stopped_trace_square':'68/15','native_stopped_choi_spectrum':{'2/15':15,'8/15':15},'actual_A5_phase_intertwiner':[str(x) for x in np.diag(D)],'conditional_instruments':embeddings,'joint_region':'a>=0,b>=0,25*a+15*b<=1','joint_example':{'six_probability':'1/50','five_probability':'1/30','remainder_probability':'71/75','choice':'explicit boundary example, not selected by P1'},'common_symmetric_Choi_sectors':[1,4,5,5],'scope':'Actual A5-equivariant finite chart and chosen coherent environment measurement; not a derived physical P1 source or selection of instrument. A measured classical record alone does not authorize arbitrary coherent POVMs.'}
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
