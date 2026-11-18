"""Actual native code versus all neutral local weight-two source zero modes.
Exact Fock calculation; no chosen drive coefficients or physical coupling.
"""
from pathlib import Path
import sys,json,math
from collections import Counter
from fractions import Fraction as F
import sympy as s
REPO=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4');sys.path.insert(0,str(REPO))
from tfpt_explorer.process import _code_basis
from tfpt_explorer.neutral_source_response import _native_neutral_branch
from tfpt_explorer.current_block_geometry import _actual_current_source

def need(ok,msg):
    if not ok:raise ArithmeticError(msg)

def gcdict(poly,ys):
    out={}
    for m,c in s.Poly(s.expand(poly),*ys).terms():
        cr,ci=s.expand(c).as_real_imag()
        need(cr.is_Integer and ci.is_Integer,'Gaussian integer Fock coefficients')
        if c:out[m]=(int(cr),int(ci))
    return out

def inner(a,b):
    re=im=0
    for m,(ar,ai) in a.items():
        br,bi=b.get(m,(0,0));w=math.prod(math.factorial(v) for v in m)
        re+=(ar*br+ai*bi)*w;im+=(ar*bi-ai*br)*w
    return re,im

def main():
    native=_native_neutral_branch();AR=native['R_matrix']
    need(s.trace(AR)==0 and s.trace(AR*AR)==96,'Actual original neutral R, not a fitted quadratic field')
    need(native['A64']==16*s.eye(8) and native['A48']==12*s.eye(8)+AR,'Original B64/B48 current-product identities')
    moments=[]
    for a in range(4):
        subset=[r for r in native['spinor'] if tuple(r[5:])==tuple(native['family'][a])]
        A16=sum((s.Matrix(r)*s.Matrix(r).T/4 for r in subset),s.zeros(8))
        Ra=4*s.eye(8)-A16
        need(len(subset)==16 and s.trace(Ra)==0 and s.trace(Ra*Ra)==96,'Each original family omission, no selected new coefficients')
        moments.append(Ra)
    need(sum(moments,s.zeros(8))==s.zeros(8),'Four original neutral R responses sum to zero')
    _,supports=_code_basis();xs=s.symbols('x0:4');ys=s.symbols('y0:8')
    # Drop only each original nonzero normalization factor, and the common
    # (1/sqrt2)^4=1/4 from the original J-positive creation coordinates.
    patterns=[sum(v**4 for v in xs),xs[0]**2*xs[1]**2+xs[2]**2*xs[3]**2,
              xs[0]**2*xs[2]**2+xs[1]**2*xs[3]**2,
              xs[0]**2*xs[3]**2+xs[1]**2*xs[2]**2,s.prod(xs)]
    norm2=[F(1,96),F(1,8),F(1,8),F(1,8),F(1)]
    for j,support in enumerate(supports):
        raw=sum(s.prod(xs[k] for k in word) for word in support)
        ratio=s.cancel(raw/patterns[j]);need(not ratio.free_symbols and ratio!=0,'Original code supports')
        need(s.simplify(ratio**2/(24*len(support))-s.Rational(norm2[j].numerator,norm2[j].denominator))==0,'Actual code normalization')
    polynomials=[s.expand(f.subs({xs[k]:ys[2*k]-s.I*ys[2*k+1] for k in range(4)})) for f in patterns]
    code=[gcdict(f,ys) for f in polynomials];weights=[v/16 for v in norm2]
    for j in range(5):
        need(F(inner(code[j],code[j])[0])*weights[j]==1,'Original normalized Fock Gram')
        need(all(inner(code[j],code[k])==(0,0) for k in range(5) if j!=k),'Code orthogonality')
    basis=[];labels=[]
    for i in range(8):
        for j in range(i,8):
            A=s.zeros(8);A[i,j]=1;A[j,i]=1;basis.append(A);labels.append([i,j])
    effects=[]
    for A in basis:
        v=[]
        for f in polynomials:
            df=s.expand(sum(A[i,j]*ys[i]*s.diff(f,ys[j]) for i in range(8) for j in range(8) if A[i,j]))
            v.append(gcdict(df,ys))
        effects.append(v)
    # Logical compression: the actual polynomial calculation must agree with
    # P dGamma(A) P=(tr A)/2 I5; no projection rule is made physical here.
    for a,A in enumerate(basis):
        tr=F(int(s.trace(A)),2)
        for j in range(5):
            for k in range(5):
                re,im=inner(code[j],effects[a][k]);need(im==0,'Real logical matrix element')
                need(re==0 if j!=k else F(re)*weights[j]==tr,'Actual logical compression')
    gram=s.zeros(36);imaggram=s.zeros(36);metric=s.diag(*[int(s.trace(A*A)) for A in basis])
    J=s.zeros(8)
    for i in range(4):J[2*i,2*i+1]=-1;J[2*i+1,2*i]=1
    for a,A in enumerate(basis):
        Ac=(A-J*A*J)/2;Aa=(A+J*A*J)/2
        for b in range(a,36):
            B=basis[b];total=F(0);imag=F(0)
            for k in range(5):
                re,im=inner(effects[a][k],effects[b][k]);total+=weights[k]*re;imag+=weights[k]*im
            # Real drive coefficients use the real quadratic form. The
            # imaginary, antisymmetric part can be nonzero in the J-polarized
            # code and must not be discarded as a claimed complex identity.
            imaggram[a,b]=s.Rational(imag.numerator,imag.denominator)
            imaggram[b,a]=-imaggram[a,b]
            total-=F(5*int(s.trace(A))*int(s.trace(B)),4)
            gram[a,b]=gram[b,a]=s.Rational(total.numerator,total.denominator)
            Bc=(B-J*B*J)/2;Ba=(B+J*B*J)/2
            expected=4*s.trace(Ac*Bc)-s.trace(A)*s.trace(B)/2+s.Rational(5,2)*s.trace(Aa*Ba)
            need(gram[a,b]==expected,'Full analytic J-split leakage formula')
    eigen=(metric.inv()*gram).eigenvals()
    need(eigen=={s.Integer(0):1,s.Integer(4):15,s.Rational(5,2):20},'Complete positive leakage classification')
    kernel=gram.nullspace();need(len(kernel)==1,'Only scalar real A preserves code infinitesimally')
    A0=sum((kernel[0][a]*A for a,A in enumerate(basis)),s.zeros(8));need(A0==A0[0,0]*s.eye(8),'Scalar kernel, not rank-only evidence')
    # The genuinely combined source is not covered by the neutral Gaussian
    # theorem. Derive its first actual local-field response instead.
    response_rows=[];response_norms=Counter()
    for raw in _actual_current_source()['chevalley'].roots:
        beta=s.Matrix(raw)/2
        need((beta.T*beta)[0]==2,'Original doubled Chevalley root convention')
        lam=(beta.T*AR*beta)[0]/2
        v=-AR*beta+lam*beta
        need((beta.T*v)[0]==0,'Actual charged derivative is transverse')
        v2=s.simplify((v.T*v)[0]);response_norms[str(v2)]+=1
        response_rows.append({'root_doubled':list(raw),'v':[str(c) for c in v],'norm_squared':str(v2)})
    need(dict(response_norms)=={'0':46,'14':96,'30':32,'32':60,'64':6},'All original mixed source-response channels')
    out={'status':'PASS_EXACT_ALL_NEUTRAL_WEIGHT_TWO_ZERO_MODES','original_code_support_sizes':[len(v) for v in supports],
         'real_symmetric_coefficient_dimension':36,'logical_compression':'P o(T_A) P=(tr A)/2 I5',
         'leakage_formula':'sum_j ||(I-P)o(T_A) f_j||^2 = 4||Acomm-(tr A)/8 I8||HS^2 + (5/2)||Aanti||HS^2',
         'J_split':'Acomm=(A-JAJ)/2; Aanti=(A+JAJ)/2',
         'leakage_spectrum_relative_to_matrix_HS_metric':{str(k):v for k,v in eigen.items()},
         'kernel':'A=a I8','R_source_trace':0,'R_source_HS_norm_squared':96,
         'R_code_compression':'0','R_total_leakage_bound':'240 <= sum_j ||(I-P)o(R)f_j||^2 <= 384, for every real isometric Cartan chart',
         'B48_code_compression':'48 I5','B64_code_compression':'64 I5','four_omitted_family_R_sum':'0; the unmarked full64 field is16 T_E8',
         'finite_operation_scope':'Any unitary evolution/sequence generated solely by neutral local weight-two zero modes preserves mode-one occupation; endpoint code preservation implies restriction in U(1)G31, hence no marked R5. Intermediate exit from Hplus is included.',
         'neutral_grade_two_dimension':44,'derivative_descendants':8,'derivative_zero_mode':'o(h(-2)Omega)=-h(0), trivial on neutral oscillator sector',
         'combined_source_first_response':{'formula':'E_beta,0 T_A=(v_beta.h(-1)) e_beta; v_beta=-A beta+(beta^T A beta)beta/2',
         'original_R_response_norm_distribution':dict(sorted(response_norms.items())),
         'nonzero_root_directions':194,'code_compression':'0 by nonzero lattice charge',
         'code_leakage_partition3':'f(b) v_beta(-3)e_beta - (D_c f)(b) beta(-3)e_beta/3; b=U^T beta, c=U^T v_beta',
         'scope':'Actual first source derivative fixed by original R and root. Does not exclude combined source evolutions or supply a physically selected interaction history.'},
         'not_excluded':['nonneutral weight-two fields','combinations with root-current zero modes','higher-weight or nonlinear mode-product generators','single-mode roundtrip K48 (not o(B48))','non-automorphism operations'],
         'basis':labels,'real_leakage_gram':[[str(v) for v in row] for row in gram.tolist()],
         'imaginary_leakage_gram':[[str(v) for v in row] for row in imaggram.tolist()]}
    out['combined_response_witnesses']=response_rows
    Path(__file__).with_name('neutral_quadratic_results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('basis','real_leakage_gram','imaginary_leakage_gram','combined_response_witnesses')},indent=2))
if __name__=='__main__':main()
