from pathlib import Path
import runpy,contextlib,io,json,itertools as it
import sympy as s
with contextlib.redirect_stdout(io.StringIO()):b=runpy.run_path(str(Path(__file__).with_name('common_construct.py')))
I=s.I
g=I*s.Matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,-1]])
h=s.Matrix([[-I,1,-I,-1],[1,I,-1,I],[1,-I,-1,-I],[I,1,I,-1]])/2
Es=[]
for a,c in it.combinations(range(4),2):
 E=s.zeros(4);E[a,c]=1;E[c,a]=-1;Es.append(E)
L=s.Matrix.hstack(*[s.Matrix(list(g.T*E*g-E)+list(h.T*E*h-E)) for E in Es])
ker=L.nullspace();assert len(ker)==1
O=sum((ker[0][j]*Es[j] for j in range(6)),s.zeros(4))
assert O.conjugate().T*O==3*s.eye(4)
paulis=[s.Matrix(P).applyfunc(s.nsimplify) for P in b['n'].PAULIS[1:]]
constraints=[]
for P in paulis:
 R=P.T*O-O*P
 constraints.append(s.Matrix([s.re(z) for z in R]+[s.im(z) for z in R]))
vs=s.Matrix.hstack(*constraints).nullspace();assert len(vs)==5
Xs=[sum((v[j]*paulis[j] for j in range(15)),s.zeros(4)) for v in vs]
Gram=s.Matrix([[s.trace(X.conjugate().T*Y)/4 for Y in Xs] for X in Xs])
Ps=[s.Matrix(P).applyfunc(s.nsimplify) for P in b['n'].SYMMETRIC_PAULIS]
D=s.Matrix(b['D']).applyfunc(s.nsimplify);order=b['order'];Ls=[]
for X in Xs:
 K=s.Matrix(10,10,lambda a,c:s.trace(Ps[a].conjugate().T*(X*Ps[c]+Ps[c]*X.T))/8)
 K=D*K.extract(order,order)*D.conjugate().T
 Ls.append(K)
V=s.Matrix.hstack(*(b['vec'](K) for K in Ls))
C=V*Gram.inv()*V.conjugate().T/3
P=C/2
assert P*P==P and s.trace(P)==5
partial=s.Matrix(10,10,lambda a,c:sum(C[j+10*a,j+10*c] for j in range(10)))
assert partial==s.eye(10)
Ugs=[s.Matrix(b['wedge'](p)) for p in [b['g5'],b['h5']]]
chars=[]
for U in Ugs:
 A=s.kronecker_product(U.conjugate(),U)
 assert A*P*A.conjugate().T==P
 chars.append(s.trace(P*A))
PH=b['C5']/2
Ps_edge=b['Psym'];Pa_edge=b['Panti']
overlapH=s.simplify(s.trace(P*PH))
overlapedge=s.simplify(s.trace(P*(Ps_edge+Pa_edge)))
Vchar=s.Matrix.hstack(*(s.Matrix(list(g*X*g.conjugate().T)) for X in Xs))
# Inner commutant dimension: invariant endomorphisms of the current five.
coeff=s.Matrix.hstack(*(s.Matrix(list(X)) for X in Xs))
acts=[]
for U in [g,h]:
 M=s.Matrix.hstack(*(s.Matrix(list(U*X*U.conjugate().T)) for X in Xs))
 acts.append((coeff.conjugate().T*coeff).inv()*coeff.conjugate().T*M)
constraints5=s.Matrix.vstack(*(s.kronecker_product(s.eye(5),A)-s.kronecker_product(A.T,s.eye(5)) for A in acts))
commdim=25-constraints5.rank();assert commdim==1
# Compare with the original UNREFINED sixty-event source, before Stay/Jump.
C60=s.zeros(100)
for U0 in b['bell']:
    U=D*s.Matrix(U0).applyfunc(s.nsimplify)*D.conjugate().T
    v=b['vec'](U);C60+=v*v.conjugate().T/60
lam=s.simplify(s.trace(P*C60)/5)
eigenspace=(C60*P==lam*P)
if eigenspace:
    assert lam>0
out={'unrefined_C60_eigenvalue_on_current_five':str(lam),'unrefined_current_five_is_eigenspace':eigenspace,'unrefined_max_current_branch':str(lam/2) if eigenspace else 'not yet determined','unique_alternating_line':1,'Omega':str(O),'Omega_dagger_Omega':'3I4','current_five_dimension':5,'current_five_commutant_dimension':commdim,'current_channel_TP':True,'current_choi_rank':5,'current_vs_slot_hodge_projector_overlap':str(overlapH),'current_choi_support_overlap_with_30_edge_support':str(overlapedge),'current_symmetric_edge_overlap':str(s.trace(P*Ps_edge)),'current_antisymmetric_edge_overlap':str(s.trace(P*Pa_edge)),'current_matching_five_overlap':str(s.trace(P*b['P5'])),'A5_current_characters_g_h':list(map(str,chars)),'scope':'Fixed native G240 selects symplectic line and irreducible current five; not physical P1 selection; current five differs from permutation five'}
Path(__file__).with_name('symplectic_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
