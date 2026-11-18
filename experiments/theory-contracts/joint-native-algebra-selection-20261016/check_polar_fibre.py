"""Actual native simplex: complete clock-compatible polar-mark fibre.
No physical intervention is stipulated. Exact rational/radical arithmetic.
"""
from pathlib import Path
import importlib.util,sys,json,itertools as it
import sympy as s
REPO=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
BASE=REPO/'experiments/theory-contracts/joint-native-algebra-selection-20261016'
sys.path.insert(0,str(REPO))
from tfpt_explorer.process import _code_basis

def clean(x): return x.applyfunc(lambda z:s.simplify(s.expand_complex(z)))
def need(x,label):
    if not x:raise ArithmeticError(label)
def zero(x):return clean(x)==s.zeros(*x.shape)
def rec(x):return [[str(v) for v in row] for row in x.tolist()]
def main():
    _,supports=_code_basis();sizes=[len(v) for v in supports]
    # Original W6, original actual code norms, no fitted orientation.
    raw=s.Matrix([[2,2,-1,-1,-1,-1],[0,0,-1,-1,1,1],[0,0,-1,1,-1,1],[0,0,1,-1,-1,1],[1,-1,0,0,0,0]])
    N=clean(s.diag(*[s.sqrt(v) for v in sizes])*raw/s.sqrt(48))
    need(zero(N.H*N-s.eye(6)+s.ones(6)/6),'Original simplex Gram')
    def action(p):
        pm=s.zeros(6)
        for j,k in enumerate(p):pm[k,j]=1
        return clean(N*pm*N.H)
    sigma=action((0,1,3,4,2,5));tau=action((1,0,2,3,4,5))
    F=clean((s.eye(5)+sigma+sigma**2)/3)
    m=N[:,2]+N[:,3]+N[:,4]
    e0=clean(s.sqrt(s.Rational(6,5))*N[:,5])
    e1=clean(s.sqrt(s.Rational(5,6))*(m+s.Rational(3,5)*N[:,5]))
    a=clean((N[:,0]-N[:,1])/s.sqrt(2));basis=a.row_join(e0).row_join(e1)
    need(zero(basis.H*basis-s.eye(3)),'Actual ON clock-fixed basis')
    need(zero(basis*basis.H-F),'Complete Fix sigma basis')
    need(zero(basis.H*tau*basis-s.diag(-1,1,1)),'Actual tau character decomposition')
    x,y=s.symbols('x y',real=True);z=x+s.I*y;v=e0+z*e1
    cc=clean(v*v.H/(1+x*x+y*y));P=clean(F-cc)
    need(zero(P*P-P) and s.simplify(s.trace(P))==2,'All finite CP1 chart projectors')
    need(zero(F*P-P) and zero(P*tau-tau*P),'Clock and weak transposition compatibility')
    R=clean((s.eye(5)-P)*(s.eye(5)-2*F)+s.I*P*tau)
    need(zero(R.H*R-s.eye(5)),'Whole rational CP1 chart R is unitary')
    need(zero(R*R-s.eye(5)+2*P),'Whole rational CP1 chart quarter-square')
    need(zero(R**4-s.eye(5)),'Whole rational CP1 chart quarter-order')
    need(s.simplify(s.trace(R))==1,'Whole rational CP1 chart trace')
    need(zero(P.subs({x:0,y:0})-(F-e0*e0.H)),'Countermark c=e0')
    # Point at infinity, absent from finite complex chart.
    Pinfty=clean(F-e1*e1.H);Rinf=clean((s.eye(5)-Pinfty)*(s.eye(5)-2*F)+s.I*Pinfty*tau)
    need(zero(Rinf.H*Rinf-s.eye(5)) and zero(Rinf**2-s.eye(5)+2*Pinfty),'CP1 infinity')
    cpol=clean(-s.sqrt(s.Rational(3,5))*e0+s.sqrt(s.Rational(2,5))*e1)
    Ppol=clean(F-cpol*cpol.H);Palt=clean(F-e0*e0.H)
    prior=json.loads((BASE/'marked_lift_results.json').read_text())
    cprior=s.Matrix(prior['original_mark']['polar_color_lines']['5'])
    need(zero(cpol-cprior),'Exact original polar mark, previous native source certificate')
    need(zero(P.subs({x:-s.sqrt(s.Rational(2,3)),y:0})-Ppol),'Original polar point in CP1 chart')
    delta=clean(Ppol-Palt);hs=s.simplify(s.trace(delta.H*delta))
    polar_response=s.simplify((e0.H*Ppol*e0)[0]);alt_response=s.simplify((e0.H*Palt*e0)[0])
    need(hs==s.Rational(4,5) and polar_response==s.Rational(2,5) and alt_response==0,'Exact distinguishing original marked-ray response')
    others=[j for j in range(6) if j!=5]
    D=N[:,others]-N[:,5]*s.ones(1,5)
    S=s.eye(5)+(s.sqrt(6)-1)*s.ones(5)/5
    Upol=clean(D*(s.eye(5)+(1/s.sqrt(6)-1)*s.ones(5)/5))
    need(zero(D.H*D-s.eye(5)-s.ones(5)) and zero(S*S-D.H*D),'Actual polar Gram and positive square root')
    need(zero(Upol.H*Upol-s.eye(5)) and zero(Upol*S-D),'Actual positive polar decomposition')
    mincost=s.simplify(s.trace((D-Upol).H*(D-Upol)))
    need(mincost==7-2*s.sqrt(6),'Actual unweighted Frobenius minimum value')
    # The original native parent is Hnat=I-(Pi5+Psym)/2=I-20M4.
    # Pi5 V=V and Psym V=V: exact native code, not a new Hamiltonian.
    V=s.zeros(256,5)
    for j,support in enumerate(supports):
        for word in support:V[64*word[0]+16*word[1]+4*word[2]+word[3],j]=1/s.sqrt(sizes[j])
    need(zero(V.H*V-s.eye(5)),'Original normalized code')
    need(all(set(it.permutations(w))<=set(support) for support in supports for w in support),'Original code is symmetric')
    PiV=V*(V.H*V)
    HnatV=V-(PiV+V)/2
    need(zero(HnatV),'Original moment parent annihilates the entire actual code')
    # On E+, the original code-level sigma and tau are both identity.
    Eplus=e0.row_join(e1)
    need(zero(Eplus.H*sigma*Eplus-s.eye(2)) and zero(Eplus.H*tau*Eplus-s.eye(2)),'Original clock algebra is scalar on the unresolved plane')
    # All q families covariant under native clock-centralizer quotient S3fix.
    marks={}
    for q in (0,1,5):
        b0=clean(s.sqrt(s.Rational(6,5))*N[:,q]);b1=clean(s.sqrt(s.Rational(5,6))*(m+s.Rational(3,5)*N[:,q]))
        vq=b0+z*b1;marks[q]=clean(F-vq*vq.H/(1+x*x+y*y))
    for i,j in ((0,1),(0,5),(1,5)):
        p=list(range(6));p[i],p[j]=p[j],p[i];g=action(p)
        for q in marks:need(zero(g*marks[q]*g.H-marks[p[q]]),'Whole CP1 q-family equivariance')
    # Original stabilizer of q and sigma normalizer fixes both e0 and e1.
    stabilizer_count=0
    for color in it.permutations((2,3,4)):
        for weak in ((0,1),(1,0)):
            p=list(range(6));p[0],p[1]=weak
            for j,k in zip((2,3,4),color):p[j]=k
            g=action(p)
            need(zero(g*e0-e0) and zero(g*e1-e1),'Retained mark/clock-normalizer action on even plane')
            stabilizer_count+=1
    # Exact single-site compression for every matrix unit, every site.
    onebody_checked=0
    for site in range(4):
        for left,right in it.product(range(4),repeat=2):
            cols=[]
            for support in supports:
                image={}
                for word in support:
                    if word[site]==right:
                        w=list(word);w[site]=left;w=tuple(w)
                        image[w]=image.get(w,0)+1
                cols.append(image)
            compressed=s.Matrix(5,5,lambda i,j:sum(cols[j].get(word,0) for word in supports[i])/s.sqrt(sizes[i]*sizes[j]))
            need(zero(compressed-(s.eye(5)/4 if left==right else s.zeros(5))),'Every original one-site matrix-unit compression')
            onebody_checked+=1
    result={'status':'PASS_EXACT_NATIVE_POLAR_FIBRE_AND_PRODUCT_PATH_CONNECTION','source':'actual W6, original _code_basis supports; no SM charges used','fixed_basis':{'a':rec(a),'e0':rec(e0),'e1':rec(e1)},'classification':'All orthogonal rank2 P inside Fix(sigma) with [P,tau]=0 and one+ / one- tau eigenvalue are F-cc*, [c] in CP1(span(e0,e1)).','finite_chart':'c line proportional e0+(x+i*y)e1; real x,y; infinity e1','clock_checks':'Full symbolic x,y: projector, tau compatibility, unitary R, R^2=I-2P, R^4=I, trace1. Infinity checked too. Spectrum follows from orthogonal block decomposition.','original_polar_chart':{'x':'-sqrt(2/3)','y':'0'},'countermark':{'c':'e0','P_HS_difference_squared':str(hs),'original_marked_ray_response_polar':str(polar_response),'original_marked_ray_response_counter':str(alt_response)},'original_polar_minimum':{'D_Gram':'I5+J5','positive_square_root':'I5+(sqrt(6)-1)J5/5','unique_unweighted_Frobenius_cost':str(mincost),'scope':'Canonical closest-isometry characterization of the original rule, not a P1-derived physical variational law.'},'native_parent':{'formula':'I-(Pi5+Psym4)/2=I-20M4','Hnat_V_exactly_zero':True,'all_logical_marks_same_ground_energy':True,'sigma_tau_on_Eplus':'I2, I2','scope':'Native moment parent and original clock algebra, not every conceivable source-dependent interaction.'},'full_family_fixedpoint_covariance':True,'q_and_clock_normalizer_elements_fix_even_plane':stabilizer_count,'single_site_matrix_units_checked':onebody_checked,'single_site_compression':'V* X_j V=(tr X)/4 I5','product_path_connection':'A=(1/4) sum_j tr(U_j* dU_j) I5; projectively flat; closed code-loop holonomy lies in projective S6 by the existing full product stabilizer theorem','scope':'Compatibility-class classification and native geometry. No physical permission for every CP1 mark, no complete P1 countermodels, no general collective/ancilla/code-deformation no-go.'}
    Path(__file__).with_name('polar_fibre_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','countermark','single_site_matrix_units_checked','q_and_clock_normalizer_elements_fix_even_plane']},indent=2))
if __name__=='__main__':main()
