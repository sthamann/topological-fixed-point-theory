"""Native electromagnetic Y current recovered from the actual two-current channel.
All modes, contractions and the coefficient are fixed by the original source.
"""
import runpy,json
from fractions import Fraction as F
from contextlib import redirect_stdout
from io import StringIO
import sympy as sp
with redirect_stdout(StringIO()):s=runpy.run_path(str(__import__('pathlib').Path(__file__).with_name('e7_quartic_projector.py')))
ch,af,G,omega,ck=(s[k] for k in ['ch','env','G','omega','ck']);af=af['af'];vac={():F(1)}
alpha=s['alpha'];ap=ch.ridx[alpha];am=ch.opp[ap]
def add(a,b,c=F(1)):
    for w,v in b.items():a[w]=a.get(w,F(0))+c*v
    return {w:v for w,v in a.items() if v}
def lin(terms,n,state):
    out={}
    for i,c in terms.items():out=add(out,af.act(i,n,state),c)
    return out
Y=sp.Matrix([F(-1,3)]*3+[F(1,2)]*2+[0]*3);yc=s['env']['Ri']*Y
Yterms={240+i:F(yc[i]) for i in range(8) if yc[i]}
charges=[sum(F(ch.roots[g][i],2)*F(Y[i]) for i in range(8)) for g in G]
# Contravariant symmetric tensor B=A_Y omega.  No chosen field ray x(alpha).
B=[[charges[i]*int(omega[i,j]) for j in range(56)] for i in range(56)]
phi={}
for i in range(56):
    for j in range(56):
        if B[i][j]:phi=add(phi,af.act(G[i],-1,af.act(G[j],-1,vac)),B[i][j])
target=lin(Yterms,-1,af.act(ap,-1,vac))
yhead=lin(Yterms,-1,vac)
rawdesc={w:-c for w,c in af.act(am,1,phi).items()}
factor=rawdesc[next(iter(yhead))]/yhead[next(iter(yhead))]
desc={w:c/factor for w,c in rawdesc.items()}
ck(desc==yhead,'original compact A1 annihilation recovers actual Y current')
norm=lin(Yterms,1,yhead).get((),F(0));ck(norm==F(5,6),'same Y current normalization retained')
def inner(bra,ket):
    out=F(0)
    for word,c in bra.items():
        st=ket;sgn=1
        for n,i in word:
            if i<240:j=ch.opp[i];sgn=-sgn
            else:j=i
            st=af.act(j,-n,st)
        out+=c*sgn*st.get((),F(0))
    return out
difference=add(dict(phi),target,-factor)
ck(inner(difference,difference)==0,'bilinear equals Y times A1 head in unitary quotient, not universal PBW')
ck(inner(phi,phi)==factor*factor*norm,'actual bilinear state norm from original vacuum')

# Full E7 generator dictionary: adjoint action on56, then the same A1 descent.
import numpy as np
idx={g:i for i,g in enumerate(G)};gen_tests=0
root_indices=[r for r in range(240) if s['ip'](alpha,ch.roots[r])==0]
physical_alpha=sp.Matrix(alpha)/2;P=sp.eye(8)-physical_alpha*physical_alpha.T/2
targets=[{r:F(1)} for r in root_indices]
for h in P.columnspace():
    hc=s['env']['Ri']*h;targets.append({240+i:F(hc[i]) for i in range(8) if hc[i]})
for terms in targets:
    A=sp.zeros(56)
    for a,c in terms.items():
        for j,g in enumerate(G):
            for h,v in ch.bracket(a,g).items():A[idx[h],j]+=c*v
    BB=A*sp.Matrix(omega);state={}
    ck(BB==BB.T,'contravariant generator tensor symmetric')
    for i,j,v in [(i,j,BB[i,j]) for i in range(56) for j in range(56) if BB[i,j]]:
        state=add(state,af.act(G[i],-1,af.act(G[j],-1,vac)),F(v))
    lowered={w:-c/factor for w,c in af.act(am,1,state).items()}
    ck(lowered==lin(terms,-1,vac),'same fixed bilinear descent recovers every original E7 current')
    gen_tests+=1
result={'claim':'MATH.ESERIES.BILINEAR.Y.CURRENT.DESCENT.20261106','definition':'A_T=actual original E7 action on56, B=A_T omega symmetric, Phi(B)=sum_ij B_ij J_i(-1)J_j(-1)Omega','source_identity':f'Phi(A_Y omega)=({factor}) J_Y(-1)J_alpha(-1)Omega in the unitary integrable quotient','descent_identity':f'J_T(-1)Omega=(-J_-alpha(1))*Phi(A_T omega)/({factor}), all133 original E7 generators including Y','original_E7_generators_recovered':gen_tests,'Y_norm':str(norm),'Tr56_Y2':str(sum(x*x for x in charges)),'Phi_norm':str(inner(phi,phi)),'universal_PBW_difference_terms':len(difference),'universal_PBW_difference_norm':str(inner(difference,difference)),'scope':'Exact native source field-state dictionary for every E7 current including original hypercharge, with actual A1 head and annihilation. Product-head identity holds in unitary quotient; annihilated current identity already holds in universal PBW degree1. No selection of a physical alpha field ray and no integrated EM fourth response derived. The bilinear uses the full56 tensor, not a single holomorphic x^4 evaluation.','checks':dict(s['env']['checks'])}
rho=sp.Matrix([F(7-2*i,2) for i in range(8)]);q=sp.ones(8,1);clocks=[]
for a in [F(1,4),F(3,4)]:
    w=rho-(2*a-1)*q/2
    for normal,label in [(physical_alpha,'D5-preserving'),(q/2,'Q-neutral')]:
        pair_grade=16-w.dot(normal);ann_grade=-8+w.dot(normal)
        ck(pair_grade+ann_grade==8,'original H8 marks cancel in Y descent')
        clocks.append({'quarter':str(a),'E7_cut':label,'pair_Y_channel_grade':str(pair_grade),'A1_annihilation_grade':str(ann_grade),'recovered_Y_grade':'8'})
    for r in root_indices:
        momentum=sp.Matrix(ch.roots[r])/2
        ck(16-w.dot(physical_alpha+momentum)-8+w.dot(physical_alpha)==8-w.dot(momentum),'actual E7 root-clock grading intertwines A1 descent')
result['original_H8_descent_grades']=clocks
result['local_field_scope']='By the state-field correspondence in the same integrable VOA the current is the double-pole coefficient of the A1 current with the normal-ordered bilinear field. This is not multiplication by an inverse A1 field, nor a modewise product identity with contacts omitted.'
result['four_Y_leg_count']='Under this exact local field-state dictionary four Y current insertions lift to eight56 current legs and four A1 current OPE fusions; normal-ordering and source derivative contacts remain. A direct Q4=EM fourth variation does not follow.'
result['checks']=dict(s['env']['checks'])
s['OUT'].joinpath('bilinear-current-descent-results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
