"""Exact geometric admission and response test of the existing elliptic source.

This changes complex structure, unlike the earlier Weyl probes. It retains
the old spin partition and tracks the original real cuts. Fixed-fiber Q is
NOT preserved away from the square point. No primitive measure is selected.
"""
from pathlib import Path
import json
import mpmath as mp
import sympy as sp

HERE=Path(__file__).resolve().parent


def require(condition,message):
    if not condition:raise RuntimeError(message)


def geometry():
    a,x,t,u,rho=sp.symbols('a x t u rho',real=True)
    f=lambda b,z:z**4+b*z*z-1
    require(sp.expand(f(-a,sp.I*x)-f(a,x))==0,'Q acts between a and -a')
    obstruction=sp.expand(f(a,sp.I*x)-f(a,x))
    require(obstruction==-2*a*x*x,'fixed quarter symmetry requires a=0')
    require(sp.simplify(-f(a,x)/x**4-f(-a,1/x))==0,'original inversion acts between a and -a')
    e=(1+sp.I)/sp.sqrt(2)
    require(sp.simplify(f(sp.I*rho,e*u)+u**4+rho*u*u+1)==0,'first real oval')
    require(sp.simplify(f(sp.I*rho,sp.conjugate(e)*u)+u**4-rho*u*u+1)==0,'second real oval')
    # Independent cross-ratio derivation of the j invariant.
    cross=((t-1)/(t+1))**2
    jlambda=256*(1-cross+cross**2)**3/(cross**2*(1-cross)**2)
    ja=16*(12-a*a)**3/(a*a+4)**2
    require(sp.factor(jlambda-ja.subs(a,-sp.I*(t+1/t)))==0,'j invariant from marked cross-ratio')
    marked=(a+2*sp.I)/(a-2*sp.I)
    require(sp.simplify(marked.subs(a,-a)-1/marked)==0,'quarter acts by inverse marked cross-ratio')
    lam=sp.symbols('lambda')
    require(sp.factor(lam*lam-1)==(lam-1)*(lam+1),'fixed quarter has no local marked modulus at lambda=-1')
    legendre=1/(1-marked)
    require(sp.simplify(legendre-(sp.Rational(1,2)+sp.I*a/4))==0,'actual Legendre parameter')
    require(sp.factor(256*(1-marked+marked**2)**3/(marked**2*(1-marked)**2)-ja)==0,
            'independent cross-ratio invariant')
    # At a0: u_branch=1,v_branch=i and their derivatives -1/4,+i/4.
    cross0=((1-sp.I)/(1+sp.I))**2
    dcross=sp.simplify(cross0*2*((-sp.Rational(1,4)-sp.I/4)/(1-sp.I)
                                              -(-sp.Rational(1,4)+sp.I/4)/(1+sp.I)))
    require(sp.simplify(cross0+1)==0 and dcross==sp.I,'nontrivial marked modulus')
    return dict(curve='E_rho: y^2=x^4+i rho x^2-1, rho real, |rho|<2',
                fixed_Q_obstruction=str(obstruction),family_Q_action='rho -> -rho',
                marked_cross_ratio='lambda(a)=(a+2i)/(a-2i); lambda(0)=-1; partial_a lambda(0)=i; partial_rho lambda(0)=-1',
                j_invariant='16(12+rho^2)^3/(4-rho^2)^2',
                real_cuts='hQkappa,hQ^3kappa remain dividing real forms; fixed oval equations v^2=u^4 +/- rho u^2+1',
                free_glide='hQ^2:(x,y)->(-x,-y) remains free',
                limitation='Q itself is a map between fibers, not the original fixed-fiber quarter clock')


def period(rho):
    return 2*mp.quad(lambda u:1/mp.sqrt(u**4+rho*u*u+1),[0,1,mp.inf])


def ratio(rho):return period(-rho)/period(rho)


def occupation(rho):return 1/(1+mp.exp(mp.pi*ratio(rho)))


def spin_kernel(a,x,z):
    d=mp.sqrt(a*a+4);rp=(-a+d)/2;rm=(-a-d)/2
    def g(q):return ((q*q-rp)/(q*q-rm))**mp.mpf('.25')
    r=g(x)/g(z)
    return (r+1/r)/(2*(x-z))


def main():
    mp.mp.dps=70
    result=geometry()
    i0=mp.gamma(mp.mpf(1)/4)**2/(2*mp.sqrt(mp.pi))
    slope=2*mp.gamma(mp.mpf(3)/4)**2/mp.gamma(mp.mpf(1)/4)**2
    require(abs(period(0)-i0)<mp.mpf('1e-60'),'primitive period beta integral')
    require(abs(mp.diff(ratio,0)-slope)<mp.mpf('1e-60'),'derived modulus variation')
    c0=occupation(0);cp=-mp.pi*slope*c0*(1-c0)
    require(abs(mp.diff(occupation,0)-cp)<mp.mpf('1e-60'),'actual boundary covariance derivative')
    # Direct derivative of the retained global Szego kernel in the fixed x
    # half-density frame. The real CAR comparison above uses normalized periods.
    xx,zz=mp.mpf(2),mp.mpf(4)
    g=lambda x:((x*x-1)/(x*x+1))**mp.mpf('.25')
    rr=g(xx)/g(zz)
    expected=(1/(4*(xx**4-1))-1/(4*(zz**4-1)))*(rr-1/rr)/(2*(xx-zz))
    derivative=mp.diff(lambda aa:spin_kernel(aa,xx,zz),0)
    require(abs(derivative-expected)<mp.mpf('1e-60'),'retained holomorphic spin kernel variation')
    # Exact two-mode density comparison, not a substitute Gaussian source.
    p,q=sp.symbols('p q',real=True);mean=(p+q)/2
    density=lambda c:sp.Matrix([c*(1-c),(1-c)**2,c*c,c*(1-c)])
    mixed=(density(p)+density(q))/2
    defect=sp.factor(mixed[3]-mean*(1-mean))
    require(defect==-(p-q)**2/4,'exact connected neutral four-field defect')
    rows=[]
    for eps in ['.1','.05','.025','.0125']:
        er=mp.mpf(eps);rp=ratio(er);rm=ratio(-er)
        legendre=(2-er)/4
        require(abs(period(er)-2*mp.ellipk(legendre))<mp.mpf('1e-60'),'independent complete elliptic period')
        require(abs(rp-mp.ellipk(1-legendre)/mp.ellipk(legendre))<mp.mpf('1e-60'),'actual marked torus modulus')
        require(abs(rp*rm-1)<mp.mpf('1e-60'),'geometric quarter exchanges periods')
        p=occupation(er);q=occupation(-er);v=(p-q)**2/4
        diag=[(p*(1-p)+q*(1-q))/2,((1-p)**2+(1-q)**2)/2,
              (p*p+q*q)/2,(p*(1-p)+q*(1-q))/2]
        require(min(diag)>0 and abs(sum(diag)-1)<mp.mpf('1e-60'),'positive normalized mixed CAR state')
        rows.append(dict(diagnostic_rho=eps,period_ratio=float(rp),c_plus=float(p),c_minus=float(q),
                         connected_four_field=float(-v),coefficient_divided_by_rho_squared=float(-v/er**2),
                         density_eigenvalues=[float(x) for x in diag]))
    result.update(primitive_period_at_square=float(i0),period_ratio_derivative=float(slope),
                  boundary_modulus='r(rho)=I(-rho)/I(rho), W=r/2; I(rho)=2 integral_0^infinity du/sqrt(u^4+rho u^2+1)',
                  independent_period='I(rho)=2K((2-rho)/4); tau=i K(1-m)/K(m); Q gives tau -> -1/tau, not modular T',
                  covariance='c_n,alpha(rho)=[1+exp(4 pi r(rho)(n+alpha))]^-1, alpha=1/4,3/4',
                  derivative_covariance='c_prime=-4 pi r_prime(0)(n+alpha)c(1-c)',
                  lowest_positive_mode_at_square=float(c0),lowest_mode_derivative=float(cp),
                  spin_kernel_derivative_at_x2_z4=float(derivative),
                  exact_mixed_four_field='kappa4=-Var(c); equal diagnostic weights at +/-rho give -(c_plus-c_minus)^2/4',
                  leading_four_field_coefficient=float(-cp**2),diagnostic_mixtures=rows,
                  new_premise='Complex structure may vary along the stated real marked-spin family; fixed-fiber Q is relaxed to its explicit family action',
                  not_derived=['primitive P1 measure support on this parameter','original weight or variance of rho',
                               'full ramified marked-equator trace','physical time or EM coupling identification',
                               'fixed-fiber Q clock away from rho=0','full interacting native or spacetime theory'],
                  source_scope='Existing free elliptic Cauchy quantization; equal-time CAR family, not a replacement of the pure Dirac-sea machine',
                  retained_compiler_scope='The supplied marked base mu4 and its fixed-fiber Q exclude rho != 0. This reuses the existing square-rigidity result; primitive P1 identification is still conditional.',
                  local_fixed_Q_consequence='Marked genus-one parity-cover deformations have one complex cross-ratio; Q sends lambda to 1/lambda. At lambda=-1 its derivative is -1, so no invariant tangent. With fixed spin/connection all remaining compatible metric variations are Weyl and preserve transported free Cauchy/Wick data.',
                  physical_gates_closed=[],complete_TFPT_solution=False)
    out=HERE/'runs/elliptic_moduli_response.json';out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(result,indent=2,ensure_ascii=False))


if __name__=='__main__':main()
