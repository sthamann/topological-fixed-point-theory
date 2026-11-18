"""Exact positive trace-gap theorem for all original adjoint-return profiles."""
from pathlib import Path
import json
import sympy as s
import check_mixed_return as m
from tfpt_explorer.cartan_source import _rays,_native_reflection

def main():
    Q=m.Q;polys,norms=m.code();derivs=[m.deriv(f,[Q(int(j==a)) for j in range(4)]) for f in polys for a in range(4)]
    roots=m._lattice_data()['roots'];rays=_rays(roots);start=(2,0,0,0,0,0,0,0);orbit={start};todo=[start]
    while todo:
        r=todo.pop()
        for j in (0,13,2,3,1):
            t=_native_reflection(r,rays[j])
            if t not in orbit:orbit.add(t);todo.append(t)
    m.need(orbit==set(roots),'Actual native code group is transitive on all240 roots')
    b=[Q.one,Q.zero,Q.zero,Q.zero];A=[];B=[]
    for k in range(1,8):
        c=[Q.zero]*4;c[k//2]=Q(0,-1) if k%2 else Q.one
        A.append([m.response(f,b,c) for f in polys]);B.append([m.response(f,b,c) for f in derivs])
    C5=s.zeros(7);C20=s.zeros(7);I5=s.zeros(7);I20=s.zeros(7)
    for i in range(7):
        for j in range(7):
            a=sum((m.inner(A[i][k],A[j][k])/norms[k] for k in range(5)),Q.zero)
            b20=sum((m.inner(B[i][k],B[j][k])/norms[k//4] for k in range(20)),Q.zero)
            C5[i,j]=s.Rational(a.x);I5[i,j]=s.Rational(a.y);C20[i,j]=s.Rational(b20.x);I20[i,j]=s.Rational(b20.y)
    m.need(C5==s.diag(s.Rational(45,8),*[s.Rational(203,24)]*6),'Complete exact real code trace form')
    m.need(C20==s.diag(s.Rational(46,3),*[s.Rational(64,3)]*6),'Complete exact real trimer trace form')
    gap=4*C5-C20;m.need(gap==s.diag(s.Rational(43,6),*[s.Rational(25,2)]*6),'Strictly positive trace gap for every real transverse vector')
    # Transfer by actual native covariance. Independently compare the implied
    # trace with EVERY previously computed original code-return matrix.
    old=json.loads(Path(__file__).with_name('mixed_return_results.json').read_text());O=s.Matrix(old['original_isometry']);AN=s.Matrix(old['original_R_native_chart']);J=s.zeros(8)
    for i in range(4):J[2*i,2*i+1]=-1;J[2*i+1,2*i]=1
    sum5=0;sum20=0;gaps=[]
    for row in old['individual_return_grams']:
        beta=m.clean(O.T*s.Matrix(row['root_doubled'])/2);v=m.clean(-AN*beta+((beta.T*AN*beta)[0]/2)*beta)
        nv=s.simplify((v.T*v)[0]);par=s.simplify((v.T*J*beta)[0]**2/2);trans=nv-par
        tr5=s.Rational(45,16)*par+s.Rational(203,48)*trans
        tr20=s.Rational(23,3)*par+s.Rational(32,3)*trans
        G=s.Matrix(row['raw_gram']);actual=sum(G[j,j]/norms[j] for j in range(5))
        m.need(s.simplify(actual-tr5)==0,'Universal trace form matches every actual root-return trace')
        delta=s.simplify(4*tr5-tr20);m.need(delta>0 if nv else delta==0,'All nonzero original channels have positive trace gap')
        sum5+=tr5;sum20+=tr20;gaps.append(str(delta))
    ext=json.loads(Path(__file__).with_name('mixed_return_extensions_results.json').read_text());w=s.Matrix(ext['witness_trimer_return'])
    m.need(s.trace(w)==s.Rational(2048,3),'Actual20matrix independently matches trace theorem witness')
    out={'status':'PASS_EXACT_ALL_POSITIVE_ORIGINAL_ADJOINT_RETURN_PROFILES','actual_native_root_orbit_size':len(orbit),
         'real_transverse_test_basis':'beta=sqrt2 e0; v_j=sqrt2 e_j, j1..7, each squared norm2',
         'real_code_trace_form':m.mat(C5),'real_trimer_trace_form':m.mat(C20),'imaginary_code_trace_form':m.mat(I5),'imaginary_trimer_trace_form':m.mat(I20),
         'universal_real_v_formulas':{'trK5':'(45/16)||v_parallel_Jbeta||^2+(203/48)||v_orthogonal_beta_Jbeta||^2','trK20':'(23/3)||v_parallel_Jbeta||^2+(32/3)||v_orthogonal_beta_Jbeta||^2','gap':'4 trK5-trK20=(43/12)||v_parallel_Jbeta||^2+(25/4)||v_orthogonal_beta_Jbeta||^2'},
         'all_original_code_traces_verified':240,'aggregate_code_trace':str(sum5),'aggregate_trimer_trace_from_theorem':str(sum20),'aggregate_gap':str(4*sum5-sum20),
         'positive_profile_conclusion':'For any nonzero complex coefficients a_beta of the original real transverse W_beta, charge orthogonality gives Kd(a)=sum |a_beta|^2 Kd,beta. Strict positive gap excludes K20(a)=K5(a) tensor I4 for every such profile, without fitting root weights.',
         'scope':'Positive actual adjoint two-insertion response with SAME source normalization and original real v_beta. Not arbitrary nonadjoint two-field kernels, complex transverse fields from different source terms, higher histories or every physical channel. Responses are not selected instruments.'}
    Path(__file__).with_name('mixed_return_trace_gap_results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if 'form' not in k or k=='universal_real_v_formulas'},indent=2))
if __name__=='__main__':main()
