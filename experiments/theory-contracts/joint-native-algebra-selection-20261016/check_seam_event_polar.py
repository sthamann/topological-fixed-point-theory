"""Actual native events, full local quartic tensors and polar normalization.

No field coupling, physical record/instrument or intermediate projection.
Previous intrinsic ten-ray result is reused, then literal original labels.
"""
from pathlib import Path
from collections import Counter
import itertools as it
import json
import math
import sympy as s
import check_mixed_return as m
from check_conformal_erasure_carry import pd, L

Q=m.Q
def cl(A): return A.applyfunc(s.simplify)
def hs(A): return s.simplify(s.trace(A.H*A))
def spectrum(A): return {str(k):v for k,v in A.eigenvals().items()}
def md(p,mon):
    for a in mon:p=pd(p,a)
    return p
def neg(p): return {tuple(a+4 for a in mon):m.conj(c) for mon,c in p.items()}
def wick(p,q,r):
    return m.add(*[(Q.one/math.prod(math.factorial(n) for n in Counter(mon).values()),
                   m.mul(md(p,mon),neg(md(q,mon))))
                  for mon in it.combinations_with_replacement(range(4),r)])

def main():
    ps,ns=m.code()
    old=json.loads((m.BASE/'full_quartic_ope_results.json').read_text())
    H=s.Matrix(old['derived_pair_eigenvalue_rows'])
    six=cl(6*H*H.T);solutions=[]
    ints=[[int(v) for v in six.row(i)] for i in range(10)]
    for signs in it.product((1,-1),repeat=10):
        if signs[0]!=1:continue
        if all(sum(ints[i][j]*signs[j] for j in range(10))==6*signs[i] for i in range(10)):
            solutions.append(cl(H.T*s.Matrix(signs)/s.sqrt(10)))
    m.need(len(solutions)==6,'Reuse intrinsic source ten-ray reconstruction, no target fit')
    S=[solutions[0]]+[w if (solutions[0].T*w)[0]<0 else -w for w in solutions[1:]]
    m.need(cl(s.Matrix.hstack(*S)*s.ones(6,1))==s.zeros(5,1),'Coherent simplex relative signs fixed by zero sum')
    raw=s.Matrix([[2,2,-1,-1,-1,-1],[0,0,-1,-1,1,1],
                  [0,0,-1,1,-1,1],[0,0,1,-1,-1,1],[1,-1,0,0,0,0]])
    N=cl(s.diag(2,s.sqrt(12),s.sqrt(12),s.sqrt(12),s.sqrt(24))*raw/s.sqrt(48))
    m.need(all(any(cl(w*w.T-s.Rational(6,5)*N[:,j]*N[:,j].T)==s.zeros(5) for w in S) for j in range(6)),
           'Original N is only literal existing label bridge to reconstructed axes')
    D=cl(N[:,:5]-N[:,5]*s.ones(1,5));G=s.eye(5)+s.ones(5)
    m.need(cl(D.T*D)==G,'Actual coherent native event response Gram')
    # A_j=sqrt(2) F(n_j) has rational Fock coefficients in original P_i.
    coeff=[s.Rational(1,24),s.Rational(1,4),s.Rational(1,4),s.Rational(1,4),s.Integer(1)]
    A=[m.add(*[(Q.from_sympy(raw[i,j]*coeff[i]),ps[i]) for i in range(5)]) for j in range(6)]
    B=[m.add((Q.one,A[j]),(-Q.one,A[5])) for j in range(5)]
    m.need(s.Matrix(5,5,lambda j,k:Q.to_sympy(m.inner(B[j],B[k]))/2)==G,'Literal original quartic event fields, not abstract rays')
    # Every transposition on reconstructed simplex, four native events per image.
    all_native=s.zeros(5);all_det=s.zeros(5);incidents=0
    for j,k in it.combinations(range(6),2):
        P=s.eye(6);P.row_swap(j,k);T=cl(N*P*N.T)
        native=cl(T*N[:,5]-N[:,5]);det=cl(-T*N[:,5]-N[:,5])
        if 5 in (j,k):
            target=k if j==5 else j
            m.need(native==D[:,target],'Actual incident transposition difference')
            incidents+=1
        else:m.need(native==s.zeros(5,1),'Nonincident native event fixes marked field')
        all_native+=4*native*native.T;all_det+=4*det*det.T
    QC=s.diag(0,0,1,1,1);QW=s.eye(5)-QC
    KC=cl(D*QC*D.T);KW=cl(D*QW*D.T);R=cl(KC+KW)
    e0=cl(s.sqrt(s.Rational(6,5))*N[:,5]);mm=N[:,2]+N[:,3]+N[:,4]
    e1=cl(s.sqrt(s.Rational(5,6))*(mm+s.Rational(3,5)*N[:,5]))
    a=cl((N[:,0]-N[:,1])/s.sqrt(2));E=s.Matrix.hstack(e0,e1)
    m.need(cl(E.T*KC*E)==s.Matrix([[s.Rational(18,5),-s.Rational(6,5)],[-s.Rational(6,5),s.Rational(2,5)]]),'Full raw C response on unresolved multiplicity plane')
    m.need(cl(E.T*KW*E)==s.Matrix([[s.Rational(12,5),s.Rational(6,5)],[s.Rational(6,5),s.Rational(3,5)]]),'Full raw W response on unresolved multiplicity plane')
    m.need(spectrum(KC)=={'4':1,'1':2,'0':2} and spectrum(KW)=={'3':1,'1':1,'0':3},'Raw responses are not orthogonal projectors')
    m.need(hs(KC*KW-KW*KC)==72,'Noncommuting full source coefficient responses')
    # Actual leading components of FULL local OPE; all Taylor terms prescribed.
    tensors={};tensor_grams={}
    for r in range(5):
        tc=m.add(*[(Q.one/2,wick(B[j],B[j],r)) for j in (2,3,4)])
        tw=m.add(*[(Q.one/2,wick(B[j],B[j],r)) for j in (0,1)])
        tensors[r]=(tc,tw)
        tensor_grams[str(r)]=[[str(Q.to_sympy(m.inner(x,y))) for y in (tc,tw)] for x in (tc,tw)]
    m.need(tensor_grams['0']==[['18','6'],['6','10']],'Uncontracted weight8 source tensors retain distinct full response')
    T=m.add(*[(Q.one,{(j,j+4):Q.one}) for j in range(4)])
    m.need(tensors[3][0]==m.add((Q(6),T)) and tensors[3][1]==m.add((Q(4),T)),'r3 tensor gives actual canonical T with derived traces')
    m.need(tensors[4]==({():Q(6)},{():Q(4)}),'r4 vacuum contractions match trace6 and4')
    # All vacuum translation descendants: no unphysical finite-locality claim.
    current=ps
    desc=[]
    for k in range(4):
        gram=s.Matrix(5,5,lambda i,j:Q.to_sympy(m.inner(current[i],current[j]))/s.sqrt(ns[i]*ns[j]))
        expected=math.factorial(k)*s.rf(8,k)
        m.need(gram==expected*s.eye(5),'Actual finite checks of analytic all-translation descendant norm')
        desc.append(str(expected));current=[L(p,-1) for p in current]
    # Strong state control: full regular local product acting on normalized mark.
    marked=[m.mul(A[5],b) for b in B]
    Gmark=s.Matrix(5,5,lambda j,k:s.Rational(3,10)*Q.to_sympy(m.inner(marked[j],marked[k])))
    m.need(Gmark==s.Rational(43,5)*s.eye(5)+s.Rational(73,5)*s.ones(5),'Full weight8 marked-state response includes all monomials')
    Gi=cl(s.eye(5)+(1/s.sqrt(6)-1)*s.ones(5)/5)
    U=cl(D*Gi);Ppol=cl(U*QW*U.T)
    cp=cl(-s.sqrt(s.Rational(3,5))*e0+s.sqrt(s.Rational(2,5))*e1)
    F=cl(a*a.T+e0*e0.T+e1*e1.T)
    m.need(cl(U.T*U)==s.eye(5) and cl(Ppol-F+cp*cp.T)==s.zeros(5),'Positive Gram normalization reproduces literal original polar marking')
    Ri=cl(s.eye(5)+(1/s.sqrt(6)-1)*e0*e0.T)
    m.need(cl(Ri*D-U)==s.zeros(5) and cl(Ri*KW*Ri-Ppol)==s.zeros(5),'Original mark is intrinsic relative response kernel')
    Gmarkwhite=cl(Gi*Gmark*Gi)
    m.need(spectrum(Gmarkwhite)=={'68/5':1,'43/5':4},'Vacuum orthonormality is not state-independent field-operation normalization')
    oblique=cl(D*QW*D.inv())
    m.need(cl(oblique*oblique)==oblique and hs(oblique-oblique.H)==2,'Original event-label cut is not selfadjoint in source metric')
    m.need(hs(G*QW-QW*G)==12,'GNS positivity does not make label cut a physical orthogonal observable')
    support=cl(D[:,:2]*(D[:,:2].T*D[:,:2]).inv()*D[:,:2].T)
    m.need(hs(support-Ppol)==(28-8*s.sqrt(6))/25,'Raw weak span not original mark')
    comms=[hs(KC*Ppol-Ppol*KC),hs(KW*Ppol-Ppol*KW)]
    m.need(comms==[(3144-384*s.sqrt(6))/625,(1044+216*s.sqrt(6))/625],'Raw spectral responses do not preserve original mark')
    # Generated response algebra is C+ C+ M2, no distinguished 3+2 central split.
    words=[s.eye(5),KC,KW];basis=[]
    for B0 in words:
        if s.Matrix.hstack(*[s.Matrix(v).reshape(25,1) for v in basis+[B0]]).rank()>len(basis):basis.append(B0)
    while True:
        added=False
        for V in list(basis):
            for B0 in (KC,KW):
                W=cl(V*B0)
                if s.Matrix.hstack(*[s.Matrix(v).reshape(25,1) for v in basis+[W]]).rank()>len(basis):basis.append(W);added=True
        if not added:break
    m.need(len(basis)==6,'Actual response algebra dimension6, with 2+1+2 decomposition')
    # Natural signed contrast also has 3+2 inertia, but a different projector.
    contrast=cl(KC-KW);negproj=cl(a*a.T+E*(3*s.eye(2)-E.T*contrast*E)*E.T/5)
    m.need(cl(negproj**2)==negproj and s.trace(negproj)==2,'Contrast negative subspace rank2')
    m.need(hs(negproj-Ppol)==(118-48*s.sqrt(6))/125,'Same 3+2 counts do not fix operator orientation')
    # Distinguish native reflection from already supplied determinant cover.
    Ddet=cl(-(N[:,:5]+N[:,5]*s.ones(1,5)))
    m.need(cl(Ddet.T*Ddet)==s.eye(5)+s.ones(5)/3,'Determinant phase changes coherent response, not a removable ray convention')
    Udet=cl(Ddet*(s.eye(5)+(s.sqrt(s.Rational(3,8))-1)*s.ones(5)/5))
    m.need(hs(Udet*QW*Udet.T-Ppol)==s.Rational(48,25),'Even determinant-cover positive polar mark is different')
    # All60 native r preserve roots; lambda*r fails original lattice already.
    from tfpt_explorer.marked_source_process import _source_operators
    _,_,rays=_source_operators()
    LD=m._lattice_data();rawroots={tuple(v) for v in LD['roots']}
    root=s.Matrix([2,0,0,0,0,0,0,0]);witness=None
    for pp in rays:
        r4=cl(s.eye(4)-2*pp);real=s.zeros(8)
        for i in range(4):
            for j in range(4):
                re,im=s.re(r4[i,j]),s.im(r4[i,j]);real[2*i,2*j]=re;real[2*i,2*j+1]=-im;real[2*i+1,2*j]=im;real[2*i+1,2*j+1]=re
        m.need(all(tuple(real*s.Matrix(v)) in rawroots for v in rawroots),'Original native event preserves all240 lattice roots')
        lam=s.zeros(8)
        for j in range(4):lam[2*j,2*j]=1/s.sqrt(2);lam[2*j,2*j+1]=1/s.sqrt(2);lam[2*j+1,2*j]=-1/s.sqrt(2);lam[2*j+1,2*j+1]=1/s.sqrt(2)
        image=cl(lam*real*root)
        m.need(tuple(image) not in rawroots,'Every actual determinant lift fails this raw root, hence Cartan-preserving VOA extension')
        if witness is None:witness=[str(x) for x in image]
    out={'status':'PASS_EXACT_ORIGINAL_EVENT_TENSORS_AND_POLAR_SCOPE','date':'2026-10-17','contract_id':'UR.NATIVE.JOINT_OPERATOR.01',
         'native_incident_transposition_images':incidents,'native_fibre_size':4,'vacuum_event_gram':m.mat(G),
         'same_profile_local_field_gram':'c_f (I+J), analytically all translation descendants, all smooth supported profiles with nonzero response',
         'translation_descendant_squared_norms_k0to3':desc,
         'raw_three_response_spectrum':spectrum(KC),'raw_two_response_spectrum':spectrum(KW),
         'raw_response_commutator_HS_squared':'72','full_wick_leading_component_grams':tensor_grams,
         'all_Taylor_terms_retained':True,'leading_r2_not_entire_pole4':True,
         'marked_state_full_regular_weight8_gram':m.mat(Gmark),'vacuum_whitened_marked_state_gram_spectrum':spectrum(Gmarkwhite),
         'positive_relative_kernel':'R^(-1/2) K_W R^(-1/2) = U Q_W U† = literal original P2; R=K_C+K_W',
         'polar_rule_physically_derived':False,'positive_normalization_is_additional_execution_rule':True,
         'raw_label_oblique_projection_antihermiticity_HS_squared':'2','label_cut_G_commutator_HS_squared':'12',
         'raw_weak_support_distance_squared':str(hs(support-Ppol)),
         'raw_C_W_commutator_squared_with_original_P2':[str(x) for x in comms],
         'actual_response_algebra':'C on color-difference2 + C on weak-difference1 + M2 on two trivial multiplicities; dimension6',
         'contrast_negative_rank':2,'contrast_mark_distance_squared':str(hs(negproj-Ppol)),
         'native_60_raw_response_sum_spectrum':spectrum(cl(all_native)),
         'determinant_60_raw_response_sum_spectrum':spectrum(cl(all_det)),
         'determinant_event_gram':'I+J/3','determinant_polar_mark_distance_squared':'48/25',
         'native_lattice_event_checks':len(rays)*len(rawroots),'determinant_lift_lattice_failures':len(rays),
         'determinant_lift_raw_root_witness':{'input':[2,0,0,0,0,0,0,0],'output':witness},
         'scope':'Given actual E8 lattice/current realization and reconstructed quartic simplex. Full local tensor identity, no intermediate projections. Event automorphy and difference fields are mathematical response, not derived coherent physical intervention, record or instrument. Vacuum and marked-source-state controls are distinct. Native seam J fixes quartic fields and all conformal descendants; boundary reflection positivity supplies a metric, not positive inverse-root processing. No all-source/P1 no-go, physical3+2, fermion/Higgs, time/rate or empirical promotion.'}
    Path(__file__).with_name('seam_event_polar_results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(out,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
