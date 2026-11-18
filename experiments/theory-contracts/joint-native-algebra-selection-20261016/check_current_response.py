"""Actual E8 brackets, clock lifts and the native tensor-dictionary test."""
from pathlib import Path
import sys,json,importlib.util,itertools as it
import sympy as s
from fractions import Fraction
REPO=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
BASE=REPO/'experiments/theory-contracts/joint-native-algebra-selection-20261016'
sys.path.insert(0,str(REPO))
from tfpt_explorer.current_block_geometry import _actual_current_source,_current,_bracket,_kappa,_ward
from tfpt_explorer.process import _code_basis,_invariants
from tfpt_explorer.charge_selection import _process_frame

def clean(m):return m.applyfunc(s.simplify)
def zero(m):return clean(m)==s.zeros(*m.shape)
def need(ok,why):
    if not ok:raise ArithmeticError(why)
def canon(v):return {k:s.simplify(x) for k,x in v.items() if s.simplify(x)!=0}
def add(*terms):
    out={}
    for coeff,v in terms:
        for k,x in v.items():out[k]=out.get(k,0)+coeff*x
    return canon(out)
def eq(a,b):return not add((1,a),(-1,b))
def module(name,path):
    sp=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m

def main():
    src=_actual_current_source();chev=src['chevalley']
    X={(i,a):_current(src,i,a) for i,a in it.product(range(5),range(4))}
    F={(i,a):_current(src,i,a,True) for i,a in X}
    Cij={(i,j):add(*[(s.Rational(1,4),_bracket(src,X[i,a],F[j,a])) for a in range(4)]) for i,j in it.product(range(5),repeat=2)}
    T=add(*[(1,Cij[i,i]) for i in range(5)])
    S={(i,j):add((1,v),(-s.Rational(i==j,5),T)) for (i,j),v in Cij.items()}
    Z=add((-s.Rational(4,9),T))
    def Q(H):
        tr=s.trace(H);h0=H-tr*s.eye(5)/5
        return add(*[(-h0[i,j],S[j,i]) for i,j in it.product(range(5),repeat=2)],(tr/5,Z))
    # Original C4 phase convention; no independent register appended.
    cw=[tuple([-1]*5+list(src['w20'][0,a][5:])) for a in range(4)]
    phases=[Fraction(1)]
    for a in range(3):
        simple=tuple(cw[a+1][j]-cw[a][j] for j in range(8))
        phases.append(phases[-1]*chev.bracket(chev.ridx[simple],chev.ridx[cw[a]])[chev.ridx[cw[a+1]]])
    C=[{chev.ridx[root]:phase} for root,phase in zip(cw,phases)]
    need(all(eq(_bracket(src,T,c),add((-s.Rational(15,4),c))) for c in C),'Actual T on original C4')
    need(all(eq(_bracket(src,Z,c),add((s.Rational(5,3),c))) for c in C),'Actual Z on original C4')
    need(all(eq(_bracket(src,Z,f),f) for f in F.values()),'Actual Z on opposite charged currents')
    # Construct current generator for every native matrix unit, then test all
    # actual charged fields. This independently establishes full index action.
    qs={};intertwining=0
    for i,j in it.product(range(5),repeat=2):
        E=s.zeros(5);E[i,j]=1;q=Q(E);qs[i,j]=q
        for k,a in F:
            expected=F[i,a] if k==j else {}
            need(eq(_bracket(src,q,F[k,a]),expected),'Actual full gl5 action on F, with compact dagger')
            intertwining+=1
        need(all(eq(_bracket(src,q,c),add((s.Rational(i==j,3),c))) for c in C),'Carrier action on actual original register')
    # Metric: from the same kappa, not a chosen coupling normalization.
    for i,j,k,l in it.product(range(5),repeat=4):
        actual=s.simplify(_kappa(src,qs[i,j],qs[k,l]));expected=s.Rational(int(j==k and i==l))-s.Rational(int(i==j and k==l),9)
        need(actual==expected,'Induced u5 current metric tr(HK)-trHtrK/9')
    data=json.loads((BASE/'polar_fibre_results.json').read_text())
    e0=s.Matrix(data['fixed_basis']['e0']);e1=s.Matrix(data['fixed_basis']['e1']);aa=s.Matrix(data['fixed_basis']['a'])
    cp=clean(-s.sqrt(s.Rational(3,5))*e0+s.sqrt(s.Rational(2,5))*e1)
    fix=clean(aa*aa.H+e0*e0.H+e1*e1.H)
    P=clean(fix-cp*cp.H);Palt=clean(fix-e0*e0.H)
    tau=clean(s.eye(5)-2*aa*aa.H)
    h=clean(cp*cp.H+P*tau/2)
    R=clean(s.eye(5)-fix-cp*cp.H+s.I*P*tau)
    need(s.trace(h)==1 and zero(R.H*R-s.eye(5)) and zero(R**4-s.eye(5)),'Actual original quarter and normalized logarithm')
    # Spectral construction of exp(i*pi*h), no numerical matrix exponential.
    need(zero(R-((s.eye(5)-fix)-cp*cp.H+s.I*(P-aa*aa.H)-s.I*aa*aa.H)),'Exact spectral clock exponential')
    qh=Q(h)
    need(s.simplify(_kappa(src,qh,qh))==s.Rational(25,18),'Exact normalized clock-current norm')
    # Noncommuting actual E+ response directions, not stipulated matrices alone.
    HX=clean(e0*e1.H+e1*e0.H);HZ=clean(e0*e0.H-e1*e1.H)
    comm=_bracket(src,Q(HX),Q(HZ))
    need(eq(comm,Q(HX*HZ-HZ*HX)) and bool(comm),'Noncommuting actual source currents')
    # Actual three-current Ward coefficients in the original root basis.
    def field(v,a,dagger):return add(*[(v[i],F[i,a] if dagger else X[i,a]) for i in range(5)])
    ypol=clean((5*P-2*s.eye(5))/6);yalt=clean((5*Palt-2*s.eye(5))/6)
    f0=field(e0,0,True);x0=field(e0,0,False);x1=field(e1,0,False)
    points=(Fraction(-1),Fraction(0),Fraction(1))
    response={}
    for label,y in [('polar',ypol),('counter',yalt)]:
        qp=Q(y)
        r00=s.simplify(_ward(src,(qp,f0,x0),points));r01=s.simplify(_ward(src,(qp,f0,x1),points))
        # Ward F3=kappa([Q,F],X)/(z12*z13*z23); this order has
        # kinematic factor -1/2 at (-1,0,1).
        kin=Fraction(1,(points[0]-points[1])*(points[0]-points[2])*(points[1]-points[2]))
        expected00=s.simplify(kin*(e0.H*y*e0)[0]);expected01=s.simplify(kin*(e1.H*y*e0)[0])
        need(r00==expected00 and r01==expected01,'Actual separated-point Ward response')
        response[label]={'q_q':str(r00),'q_even':str(r01),'matrix_q_q':str(r00/kin),'matrix_q_even':str(r01/kin)}
    need(response['polar']['q_q']=='0' and response['counter']['q_q']=='1/6','Exact discriminating source response')
    # Complete source lift branches: C phase cubed =det R = -1. On all F,
    # each branch is the same R. Full24 generates E8, so the field action
    # fixes the full automorphism; no unverified group-order inference.
    branches=[]
    for k in range(3):
        c=s.exp(s.I*s.pi*(2*k+1)/3).expand(complex=True)
        need(s.simplify(c**3)==-1,'Exact cubic source-lift constraint')
        period=next(n for n in range(1,13) if n%4==0 and s.simplify(c**n-1)==0)
        branches.append({'trace_h':str(2*k+1),'C_phase':str(c),'full_source_order':period,'triple_C_phase':str(s.simplify(c**3))})
    need([v['full_source_order'] for v in branches]==[12,4,12],'Unique full-source order4 branch')
    # Actual native field dictionary: W is exactly Sym3, but independent
    # source SU5 acts trivially on C and nontrivially on F.
    audit=module('native_audit_current',BASE/'operator_checker.py');_,supports=_code_basis();W=audit.erasure_isometry(supports)
    need(zero(W.H*W-s.eye(20)),'Actual W rank20')
    need(all(set(it.permutations(word))<=set(support) for support in supports for word in support),'Actual symmetric code supports')
    residual=clean(-W-W*s.kronecker_product(R,s.eye(4)))
    need(residual.rank()>0,'Every source clock branch fails the literal native tensor dictionary')
    mismatch2=s.simplify(s.trace(residual.H*residual))
    # All sixty native source-word factors: distinguish finite determinant
    # path from the direct original order2 word. Verify the actual existing
    # event chart too, not an assumed same-dimensional representation.
    rankcheck=module('native_integer_current',BASE/'check_hidden_gauge.py');inv=_invariants();frame=_process_frame()
    norms=s.diag(*[s.sqrt(len(v)) for v in supports]);invnorm=norms.inv()
    matches=0;events=0
    phase=(1-s.I)/s.sqrt(2)
    for ray,label in zip(inv['hamming_rays'],inv['labels']):
        zc=[(int(v.real),int(v.imag)) for v in ray]
        twice=[]
        for i in range(4):
            row=[]
            for j in range(4):
                real,imag=rankcheck.gaussian_multiply(zc[i],(zc[j][0],-zc[j][1]));row.append((2*int(i==j)-real,-imag))
            twice.append(row)
        num,ok=rankcheck.exact_quartic_code_action(twice,supports)
        need(ok,'Actual native quartic source-word action')
        logical=clean(norms*num*invnorm/16)
        pair=inv['transposition_pairs'][inv['logical_to_pair'][label]]
        pp=s.eye(6);pp.row_swap(*pair);existing=clean(frame.T*pp*frame)
        if zero(logical-existing):matches+=1
        # Native h^3 phase and opposite-current determinant-path phase.
        need(s.simplify(phase**3+ s.conjugate(phase))==0,'Actual lifted opposite-current phase equals native cubic phase')
        events+=1
    need(matches==60,'All actual source current event charts match native quartic action exactly')
    out={'status':'PASS_EXACT_ORIGINAL_CURRENT_RESPONSE_AND_TENSOR_SCOPE','actual_gl5_field_brackets':intertwining,'matrix_unit_metric_checks':625,'induced_current_metric':'kappa(Q(H),Q(K))=tr(HK)-tr(H)tr(K)/9','source_clock_current_norm_squared':'25/18','noncommuting_Eplus_source_currents':True,'ward_positions':['-1','0','1'],'actual_three_current_responses':response,'source_lift_branches':branches,'actual_native_event_chart_matches':matches,'finite_link_scope':'Same oriented determinant-path lift h=exp(-i*pi/4)r on the coupled finite cover. Opposite E8 currents F transform as A5 tensor bar(h); direct native order2 word is not the same charged lift.','native_tensor_dictionary':{'W_image':'Sym3 C4, dimension20','quarter_mismatch_HS_squared':str(mismatch2),'all_source_clock_branches_fail':True,'infinitesimal_counterexample':'Every nonzero traceless H: Q(H) acts0 on original C4 and all its products, but H tensor I4 on the opposite charged F20.','scope':'Rules out identifying source C4 with the unchanged native register AND preserving the actual quartic/erasure tensor map under the whole source current operation. Not a general no-go for nonlinear/composite field maps or other collective source dynamics.'},'scope':'Existing stipulated E8_1 source, original Chevalley brackets, compact dagger and Ward vacuum. No raw-P1 identification, physical source time, detector instrument or coupling inserted/derived.'}
    Path(__file__).with_name('current_response_results.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
