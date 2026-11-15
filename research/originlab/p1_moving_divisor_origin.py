"""Exact original regional sl2 moving-divisor and four-mark closure tests.

Fixed four-mark sections do not carry the full infinitesimal Möbius action:
raw Lie derivatives create double poles. The correct source-adjacent
action transports the divisor. All results are finite mathematical
controls; original weighted Q0/R113, physical compact form and F1 open.
"""
from itertools import combinations
from math import comb
import sympy as sp


def require(ok, message):
    if not ok:
        raise ArithmeticError(message)


def polynomial_sl2(n):
    E,F=sp.zeros(n+1),sp.zeros(n+1)
    for k in range(1,n+1):
        E[k-1,k]=k
    for k in range(n):
        F[k+1,k]=k-n
    H=sp.diag(*[sp.Rational(k)-sp.Rational(n,2) for k in range(n+1)])
    G=sp.diag(*[sp.Rational(1,comb(n,k)) for k in range(n+1)])
    require(E*F-F*E==2*H and H*E-E*H==-E and H*F-F*H==F,
            "transport sl2 polynomial representation failed")
    require(E.T*G==-G*F and H.T*G==G*H,
            "compact SU2 binomial metric failed")
    return E,F,H,G


def transvectant(f,g,order,x,y):
    return sp.expand(sum((-1)**k*sp.binomial(order,k)*
        sp.diff(f,x,order-k,y,k)*sp.diff(g,x,k,y,order-k)
        for k in range(order+1)))


def controls():
    z=sp.symbols("z")
    marks=(sp.S.One,-sp.S.One,sp.I,-sp.I)
    P=sp.prod(z-p for p in marks)
    require(sp.expand(P)==z**4-1,"original mu4 divisor changed")
    source_fields=(z*z,(z-1)**2/2,sp.S.One,(z+1)**2/2)
    source_rank=sp.Matrix([[sp.Poly(v,z).nth(k) for v in source_fields]
                           for k in range(3)]).rank()
    require(source_rank==3,"original four regional fields do not span sl2")
    require(sp.limit((z-1)**2*sp.diff(1/P,z),z,1)==-sp.Rational(1,4),
            "fixed-divisor translation has no double pole")
    A,B,C=sp.symbols("A B C")
    require(sp.linsolve([A+B*p+C*p*p for p in marks],(A,B,C))==
            sp.FiniteSet((0,0,0)),
            "four-mark divisor has infinitesimal PSL2 stabilizer")

    a,b,c,d=sp.symbols("a b c d")
    pts=(a,b,c,d)
    cross=(a-c)*(b-d)/((a-d)*(b-c))
    for v in (sp.S.One,z,z*z):
        require(sp.factor(sum(sp.diff(cross,p)*v.subs(z,p) for p in pts))==0,
                "regional Möbius flow changes projective cross ratio")

    # a_i'=v(a_i) in the SAME physical moving-divisor family.
    # The transport of f_k=z^k/P and omega_j=z^j dz/P has
    # correction sum_i [v(a_i)-v(z)]/[z-a_i].
    moving=[]
    for v in (sp.S.One,z,z*z):
        sig=sp.cancel(sum((v.subs(z,p)-v)/(z-p) for p in marks))
        expected={sp.S.One:sp.S.Zero,z:-4,z*z:-4*z}[v]
        require(sig==expected,"moving-divisor pole cancellation changed")
        for k in range(5):
            numerator=sp.cancel(v*sp.diff(z**k,z)+sig*z**k)
            require(sp.degree(numerator,z)<=4,
                    "transported meromorphic mode leaves C5")
        for j in range(3):
            numerator=sp.cancel(v*sp.diff(z**j,z)+sp.diff(v,z)*z**j+sig*z**j)
            require(sp.degree(numerator,z)<=2,
                    "transported logarithmic form leaves F3")
        moving.append(str(sig))

    E4,F4,H4,G4=polynomial_sl2(4)
    E2,F2,H2,G2=polynomial_sl2(2)
    I=sp.Matrix(5,3,lambda i,j:int(i==j+1))
    require(H4*I==I*H2 and
            E4*I-I*E2==sp.Matrix(5,3,lambda i,j:int(i==j)) and
            F4*I-I*F2==-sp.Matrix(5,3,lambda i,j:int(i==j+2)),
            "D4 evaluation sequence incorrectly became SL2 equivariant")
    require(E4[:,2]!=sp.zeros(5,1) and F4[:,2]!=sp.zeros(5,1),
            "f2 compensator incorrectly treated as SL2 invariant")
    r5=sp.diag(1,sp.I,-1,-sp.I,1)
    s5=-sp.Matrix(5,5,lambda i,j:int(i+j==4))
    r3=sp.diag(sp.I,-1,-sp.I)
    s3=sp.Matrix(3,3,lambda i,j:int(i+j==2))
    require(all(r.H*G*r==G and s.T*G*s==G
                for r,s,G in ((r5,s5,G4),(r3,s3,G2))),
            "positive compact binomial metrics violate D4")
    require(H4[0,0]==-2 and H2[0,0]==-1 and
            H4.T*G4+G4*H4!=sp.zeros(5),
            "noncompact SL2 dilation incorrectly declared unitary")

    x,y=sp.symbols("x y")
    quartics=[x**(4-k)*y**k for k in range(5)]
    ops=(lambda f:sp.expand(x*sp.diff(f,y)),
         lambda f:sp.expand(y*sp.diff(f,x)),
         lambda f:sp.expand(x*sp.diff(f,x)-y*sp.diff(f,y)))
    checks=0
    for f in quartics:
        for g in quartics:
            t=transvectant(f,g,3,x,y)
            require(t==-transvectant(g,f,3,x,y),
                    "third transvectant lost antisymmetry")
            for op in ops:
                require(sp.expand(op(t)-transvectant(op(f),g,3,x,y)-
                                  transvectant(f,op(g),3,x,y))==0,
                        "third transvectant is not SL2 equivariant")
                checks+=1
    require(checks==75 and transvectant(x**4,y**4,3,x,y)==576*x*y,
            "quartic-to-quadratic equivariant bracket changed")

    sext=(x**6,x**5*y,x*y**5)
    jac=sp.expand(sum(transvectant(sext[i],
                 transvectant(sext[(i+1)%3],sext[(i+2)%3],3,x,y),3,x,y)
                 for i in range(3)))
    require(jac==5184000*x**6,"n=8 higher bracket Jacobi witness changed")
    exact_marks=[n for n in range(2,33)
                 if (n-1)**2==1+(n-1)+(n+1)]
    require(exact_marks==[4],"unital End(F_n) closure not unique at n=4")
    return {
      "status":"EXACT_MOVING_DIVISOR_AND_PRIMITIVE_F1_SELECTION",
      "original_four_regional_fields_span_sl2":True,
      "fixed_D_naive_Lie_derivative_double_pole":"-1/4 at z=1 for translation of f0",
      "fixed_four_point_SL2_invariant_vector_field_dimension":0,
      "moving_PSL2_cross_ratio_derivatives":[0,0,0],
      "moving_divisor_numerator_corrections_for_1_z_z2":moving,
      "C5_generators":"E=partial_z, F=z^2partial_z-4z, H=zpartial_z-2",
      "F3_generators":"E=partial_z, F=z^2partial_z-2z, H=zpartial_z-1",
      "C5_positive_compact_SU2_metric_diagonal":[str(v) for v in G4.diagonal()],
      "F3_positive_compact_SU2_metric_diagonal":[str(v) for v in G2.diagonal()],
      "positive_nontrivial_finite_fixed_SL2R_unitary":False,
      "canonical_F3_to_C5_inclusion_full_sl2_equivariant":False,
      "f2_compensator_invariant_line_under_full_regional_sl2":False,
      "third_transvectant_sl2_equivariance_checks":checks,
      "nonzero_quartic_binary_third_transvectant":"576*x*y",
      "n8_unique_wedge_bracket_Jacobi_failure":"5184000*x**6",
      "nontrivial_sl2_equivariant_internal_Lie_on_log_Fn_forces_n":4,
      "unital_End_Fn_equals_scalar_plus_Fn_plus_Cn_for_n":[4],
      "physical_moving_divisor_connection_identified":False,
      "physical_compact_real_form_and_source_Gram_derived":False,
      "native_weighted_Q0_R113_produced":False,
      "F1_solved":False,
    }


if __name__=="__main__":
    import json
    print(json.dumps(controls(),indent=2,sort_keys=True))
