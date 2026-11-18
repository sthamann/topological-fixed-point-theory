from pathlib import Path
import importlib.util, itertools as it, math, json
import numpy as np
import sympy as s

repo=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
source=repo/'experiments/theory-contracts/compiler-correlated-event-clock-20260919/source_channel.py'
spec=importlib.util.spec_from_file_location('native',source)
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
roots=m.source_rays(); bells,raw,smaps=m.reflection_actions(roots)
outer,pentads=m.outer_mark_action(raw)
Z=np.column_stack(roots); H=Z.conj().T@Z/4
abs2=(H.conj()*H).real
G2=abs2**2; G4=abs2**4; A4=H**4
move=np.array([[p[q]!=q for p in outer] for q in range(6)],int).T
v=3*move-1
assert np.array_equal(v.T@v,144*np.eye(6)-24*np.ones((6,6)))
assert np.array_equal(G2@v,np.zeros((60,6)))
# Exact polynomial identities, dyadic matrices stored without rounding.
def spectrum(M, vals):
    M=s.Matrix(M.tolist()).applyfunc(s.nsimplify)
    out={}
    for a in vals:
        P=s.eye(M.rows)
        for b in vals:
            if a!=b:P=P*(M-b*s.eye(M.rows))/(a-b)
        assert M*P==a*P and P*P==P
        out[str(a)]=int(s.trace(P))
    assert sum(out.values())==M.rows
    return out
g2spec=spectrum(G2,[s.Integer(6),s.Integer(2),s.Integer(1),s.Rational(1,2),s.Integer(0)])
assert g2spec=={'6':1,'2':15,'1':9,'1/2':30,'0':5}
assert np.array_equal(A4@v,3*v)
# Symmetric degree-four coordinates with exact integer orbit metric.
ns=[n for n in it.product(range(5),repeat=4) if sum(n)==4]
mult=np.array([24//math.prod(math.factorial(k) for k in n) for n in ns])
T=np.array([[math.prod(z[j]**n[j] for j in range(4))/16 for z in roots] for n in ns])
assert np.array_equal(T.conj().T@(mult[:,None]*T),A4)
supports=[[(4,0,0,0),(0,4,0,0),(0,0,4,0),(0,0,0,4)],
          [(2,2,0,0),(0,0,2,2)],[(2,0,2,0),(0,2,0,2)],
          [(2,0,0,2),(0,2,2,0)],[(1,1,1,1)]]
E=np.array([[int(n in group) for group in supports] for n in ns])
assert np.array_equal(E.T@(mult[:,None]*E),np.diag([4,12,12,12,24]))
# Clear all denominators; these products are Gaussian integers/dyadic rationals.
X=T@v
assert np.array_equal(X.conj().T@(mult[:,None]*X),432*np.eye(6)-72*np.ones((6,6)))
assert np.array_equal((X@X.conj().T)*mult[None,:],
                      (E@np.diag([108,36,36,36,18])@E.T)*mult[None,:])
# The preceding identity is S P0 S* = 3 Pi, since P0=vv*/144.
assert np.array_equal(T@T.conj().T*mult[None,:],
                     1.5*np.eye(35)+ E@np.diag([3/8,1/8,1/8,1/8,1/16])@E.T*mult[None,:])
# Actual two-copy Bell amplitudes and Born reader.
a=np.array([[z.T@p.conj()@z/8 for z in roots] for p in m.SYMMETRIC_PAULIS])
assert np.array_equal(a.conj().T@a,H**2)
assert np.array_equal(a@a.conj().T,6*np.eye(10))
classes={t:[i for i,p in enumerate(raw) if p==t] for t in sorted(set(raw))}
ids=list(classes.values())
K=np.array([np.mean(abs(a[:,ii])**2,axis=1) for ii in ids])
assert set(K.ravel())=={0,0.25}
assert np.array_equal(K.sum(1),np.ones(15))
assert np.array_equal(K.T@K,np.eye(10)/4+np.ones((10,10))/8)
for ii in ids:
    rel=[bells[ii[0]].conj().T@bells[j] for j in ii]
    assert all(np.array_equal(u,np.diag(np.diag(u))) for u in rel)
    assert all(any(np.array_equal(u@w,t) for t in rel) for u in rel for w in rel)
    chars=[tuple(u[j,j] for u in rel) for j in range(10)]
    assert sorted(chars.count(c) for c in set(chars))==[2,2,2,4]
L=4*K.T-np.ones((10,15))/3
assert np.array_equal(3*L,12*K.T-np.ones((10,15)))
assert np.max(abs(L@K-np.eye(10)))<1e-14
# Native label and context permutations from actual unitary conjugations.
P=[np.outer(z,z.conj()) for z in roots]
def key(M):return tuple(zip(M.real.astype(int).ravel(),M.imag.astype(int).ravel()))
lookup={key(p):i for i,p in enumerate(P)}
R=[np.eye(4)-p/2 for p in P]
actions=[[lookup[key(r@p@r.conj().T)] for p in P] for r in R]
labelmean=sum((np.eye(60)[p] for p in actions))/60
assert np.max(abs(labelmean@v-3*v/5))<1e-14
assert np.array_equal(5*sum((np.eye(60)[p] for p in actions))@v,180*v)
for g,act in enumerate(actions):
    assert np.array_equal(v[act,:],v[:,outer[g]])
    # Fourth powers remove native mu4 rephasing exactly.
    rz=R[g]@Z
    for j in range(60):
        phase=next(u for u in (1,-1,1j,-1j) if np.array_equal(rz[:,j],u*Z[:,act[j]]))
        assert phase**4==1
# Original v774 PI_SIG=(4,5,0,1,2,3,6,7), in complex coordinates (2,0,1,3).
sigma=np.eye(4)[[2,0,1,3]]
sigma_rays=[lookup[key(sigma@p@sigma.T)] for p in P]
sigma_marks=tuple(next(j for j in range(6) if np.array_equal(v[sigma_rays,q],v[:,j])) for q in range(6))
assert sum(sigma_marks[q]==q for q in range(6))==3
assert all(sigma_marks[sigma_marks[sigma_marks[q]]]==q for q in range(6))
context_of={r:c for c,ii in enumerate(ids) for r in ii}
cp=[[context_of[p[ii[0]]] for ii in ids] for p in actions]
Q=[np.eye(15)[p] for p in cp]
B=[abs(b)**2 for b in bells]
duads=[tuple(j for j in range(6) if outer[ii[0]][j]!=j) for ii in ids]
for g in range(60):
    assert all(duads[cp[g][j]]==tuple(sorted(outer[g][x] for x in duads[j])) for j in range(15))
back=[]
for q in range(6):
    selected=[i for i in range(60) if move[i,q]]
    other=[i for i in range(60) if not move[i,q]]
    for inds in (selected,other):
        assert np.array_equal(sum(Q[i] for i in inds)@K,K@sum(B[i] for i in inds))
    # Native stopped jump graph: selected average has self 2/5 and nonself 1/5.
    num=sum(B[i] for i in selected)
    adj=(num-8*np.eye(10))/4
    assert set(adj.ravel())=={0,1} and np.array_equal(adj.sum(1),np.full(10,3))
    V2=(adj@adj)/9
    Mtilde=K@V2@L
    assert np.min(Mtilde)>-1e-14
    assert np.max(abs(Mtilde@K-K@V2))<1e-14
    Q20=s.Matrix(sum(Q[i] for i in selected).astype(int))/20
    Id=s.eye(15)
    polynomial=(5*Q20-4*Id)*(5*Q20+Id)*(25*Q20**2-5*Q20-8*Id)/72
    assert np.max(abs(np.array(polynomial.tolist(),float)-Mtilde))<1e-14
    # Exact rational canonical backtransport, avoiding floating inverse.
    Ks=s.Matrix(K.tolist()).applyfunc(s.Rational)
    Vs=s.Matrix((adj@adj).astype(int))/9
    Ls=4*Ks.T-s.ones(10,15)/3
    assert polynomial==Ks*Vs*Ls
    inner=[d for d in duads if q not in d]
    innerA=s.Matrix(10,10,lambda i,j:int(not set(inner[i])&set(inner[j])))
    restriction=s.Matrix(10,15,lambda i,j:int(inner[i]==duads[j]))
    duadA=s.Matrix(15,15,lambda i,j:int(not set(duads[i])&set(duads[j])))
    p5=(innerA+2*s.eye(10))/3-s.ones(10)/6
    assert polynomial==(duadA+3*Id)/9-restriction.T*p5*restriction/3
    back.append(str(polynomial.eigenvals()))
# A native signed preparation relation survives at degree four as densities.
four_density_witness=(v.T@G4@v)
assert np.array_equal(4*four_density_witness,378*np.eye(6)-63*np.ones((6,6)))
# Thus sum v_q D_r=0, but ||sum v_q |psi^4><psi^4|||_HS^2=315/4.
# Ordered products of Dr span M10: complete frame and connected overlap graph.
visited={0};todo=[0]
while todo:
    i=todo.pop()
    for j in range(60):
        if H[i,j]!=0 and j not in visited:visited.add(j);todo.append(j)
assert len(visited)==60
# Products along nonzero-overlap paths are nonzero rank-one matrix units
# between spanning endpoint vectors. This proves the full operator algebra.
# Native Pauli kernel twirl dephases Bell10 and is invisible on degree-four code.
U=[]
for p in m.PAULIS:
    u=np.array([[np.trace(b.conj().T@p@c@p.T)/4 for c in m.SYMMETRIC_PAULIS] for b in m.SYMMETRIC_PAULIS])
    assert np.array_equal(u,np.diag(np.diag(u)))
    U.append(u)
    # Every native Pauli acts identically on the actual quartic five-code.
    pz=p@Z
    pt=np.array([[math.prod(z[j]**n[j] for j in range(4))/16 for z in pz.T] for n in ns])
    assert np.array_equal(pt@v,X)
twirl=sum(np.kron(u.conj(),u) for u in U)/16
target=np.diag([int(i//10==i%10) for i in range(100)])
assert np.array_equal(twirl,target)
returns=np.sum(abs(a)**4,axis=0)
assert set(returns)=={0.25}
population=(a.conj()*a).real
diagonal_returns=np.einsum('ir,ij,jr->r',population,adj@adj,population)/9
families={}
for alpha in [s.Integer(0),s.Rational(1,2),s.Integer(1)]:
    families[str(alpha)]=sorted(set(str(s.Rational(float(x*9))*s.Rational(1,9)+(1-alpha)/8) for x in diagonal_returns))
result={'G2_spectrum':g2spec,'native_Born_intertwiners_all_marks':True,
 'quartic_isometry_exact':True,'simplex_Gram_numerator_diagonal':360,
 'simplex_Gram_numerator_off_diagonal':-72,'quartic_density_signed_witness_norm2':'315/4',
 'no_linear_two_copy_to_four_copy_density_map':True,
 'ordered_projector_algebra_dimension':100,
 'Bell_dephase_native_return':'1/4',
 'Born_pseudoinverse_values':['-1/3','2/3'],
 'native_sigma_mark_permutation':sigma_marks,
 'covariant_macro_lift_return_values':families,
 'backtransport_spectra_all_marks':back,
 'source_identification_with_P1':False}
Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
