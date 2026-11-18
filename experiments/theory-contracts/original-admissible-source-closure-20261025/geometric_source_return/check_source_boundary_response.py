"""Actual NS cover: half-collar Cauchy projection and its first response.

Uses the already stipulated free Green quantization. It does not identify
the x-time cut with the primitive P1 seam or with the old quarter Hardy time.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def zero(expression):
    if isinstance(expression, sp.MatrixBase):
        return all(sp.simplify(x) == 0 for x in expression)
    return sp.simplify(expression) == 0


def annihilator(index, count):
    result = sp.zeros(2**count)
    for bits in range(2**count):
        if (bits >> index) & 1:
            sign = (-1)**((bits & ((1 << index)-1)).bit_count())
            result[bits ^ (1 << index), bits] = sign
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parent/'runs/source_boundary_response.json')
    args = parser.parse_args()
    r = sp.symbols('r', positive=True)
    energy = sp.symbols('E', real=True)
    identity = sp.eye(2)
    P = sp.Matrix([[1, r], [r, r*r]])/(1+r*r)
    v = sp.Matrix([1, r])/sp.sqrt(1+r*r)
    other = sp.Matrix([-r, 1])/sp.sqrt(1+r*r)
    C, cross = r*r/(1+r*r), r/(1+r*r)
    need(zero(P-v*v.T), 'Cauchy graph projection')
    need(zero(P*P-P) and zero(P.T-P), 'Orthogonal projection')
    need(zero(other*other.T-(identity-P)), 'Other NS half is complement')
    need(zero(P*other), 'NS halves are orthogonal')
    need(zero(P[1,1]-C) and zero(P[0,1]-cross), 'Restriction and cross response')
    need(zero(cross*cross-C*(1-C)), 'First exit and return response')
    need(zero(P[0,0]**2+cross**2-P[0,0]), 'Full first return identity')
    R = 2*P-identity
    need(zero(R*R-identity), 'Cauchy polarization involution')
    need(zero(R-sp.Matrix([[(1-r*r)/(1+r*r),2*r/(1+r*r)],
                          [2*r/(1+r*r),-(1-r*r)/(1+r*r)]])), 'Tanh/sech form')
    # E_s=e^-s E, hence dr/ds=E r at the original geometry.
    derivative = sp.diff(P, r)*energy*r
    dC = sp.diff(C,r)*energy*r
    need(zero(dC-2*energy*C*(1-C)), 'Source-derived covariance response')
    need(zero(P*derivative+derivative*P-derivative), 'Projection tangent identity')
    need(zero(derivative.trace()), 'Rank trace response')
    need(zero(P.subs(r,1)-sp.ones(2)/2), 'Actual zero mode')
    # Negative controls: covariance is not itself a pure projection; two common
    # coordinate shortcuts must fail away from zero energy.
    need(not zero(C*C-C), 'Mixed covariance mistakenly called pure')
    swap = sp.Matrix([[0,1],[1,0]])
    need(not zero(swap-(2*P-identity)), 'Raw endpoint exchange mistaken for polarization')
    periodic_wrong = sp.Matrix([[1,r],[r,r*r]])/(1-r*r)
    need(not zero(periodic_wrong**2-periodic_wrong), 'NS sign removed undetected')
    frequency=sp.symbols('k', real=True)
    G=1/(sp.I*frequency+energy)
    need(zero(-energy*sp.diff(G,energy)-G*energy*G), 'Full Green first variation')
    need(zero(sp.diag(energy,energy)*P-P*sp.diag(energy,energy)), 'Boundary time covariance')

    # Independent CAR realization on the actual two boundary traces for two
    # spectral modes. Rational r values verify identities without tolerance.
    a = [annihilator(i,4) for i in range(4)]
    for i in range(4):
        for j in range(4):
            need(zero(a[i]*a[j].T+a[j].T*a[i]-(sp.eye(16) if i==j else sp.zeros(16))), 'CAR')
    r1, r2 = sp.Rational(1,2), sp.Rational(2,3)
    vacuum = sp.zeros(16,1); vacuum[0]=1
    state = (a[0].T+r1*a[1].T)*(a[2].T+r2*a[3].T)*vacuum/sp.sqrt((1+r1*r1)*(1+r2*r2))
    need(zero((state.T*state)[0]-1), 'Pure state norm')
    expected = sp.diag(P.subs(r,r1),P.subs(r,r2))
    for i in range(4):
        for j in range(4):
            need(zero((state.T*a[i].T*a[j]*state)[0]-expected[i,j]), 'Boundary two point')
    for i in range(4):
        for j in range(4):
            for k in range(4):
                for l in range(4):
                    actual=(state.T*a[i].T*a[j].T*a[l]*a[k]*state)[0]
                    wick=expected[i,k]*expected[j,l]-expected[i,l]*expected[j,k]
                    need(zero(actual-wick), 'Boundary four point Wick identity')

    # Original transverse eigenfunctions, same b boundary domain, all-mode
    # formula checked on both b and several signed modes.
    w = sp.symbols('w', real=True)
    sigma1=sp.Matrix([[0,1],[1,0]]); sigma3=sp.diag(1,-1)
    for b in [-1,1]:
        for n in range(-3,4):
            lam=sp.pi*n/2
            phi=sp.Matrix([sp.exp(-sp.I*lam*w),b*(-1)**n*sp.exp(sp.I*lam*w)])/2
            need(zero(sp.I*sigma3*sp.diff(phi,w)-lam*phi),'Original spatial operator')
            for endpoint in [-1,1]:
                need(zero((sigma1-b*identity)*phi.subs(w,endpoint)), 'Fixed boundary domain')

    sigma2=sp.Matrix([[0,-sp.I],[sp.I,0]])
    need(zero(sigma1*sigma2*sigma1+sigma2), 'Normal metric variation preserves Pin covariance')
    # The existing transverse zero branch of the PURE D-sea realizes the old
    # quarter circle. It is a different source quantization from thermal C.
    x=sp.symbols('x', real=True)
    scale=sp.symbols('a', positive=True)
    pure_bridge_checks=0
    for b in [-1,1]:
        phi0=sp.Matrix([1,b])/2
        need(zero(sp.integrate((phi0.T*phi0)[0],(w,-1,1))-1),'Zero branch isometry')
        need(zero(sp.integrate(w*(phi0.T*phi0)[0],(w,-1,1))),'Zero mean normal coordinate')
        need(zero(sp.integrate(w*w*(phi0.T*phi0)[0],(w,-1,1))-sp.Rational(1,3)), 'Normal nonlocality witness')
        for pin in [sp.I,-sp.I]:
            quarter=pin*b
            alpha=sp.Rational(1,4) if quarter==sp.I else sp.Rational(3,4)
            for nlong in [-1,0,1]:
                f=sp.exp(2*sp.pi*sp.I*(nlong+alpha)*x)
                psi=phi0*f
                applied=-sp.I*(sigma1*sp.diff(psi,x)+scale*sigma2*sp.diff(psi,w))
                need(zero(applied-phi0*(-sp.I*b*sp.diff(f,x))), 'Pure source same-time intertwiner')
                need(zero(psi.subs(x,x+1)-pin*sigma1*psi.subs(w,-w)), 'Actual quarter Pin domain')
                need(zero(sp.conjugate(psi.subs(x,-x))-phi0*sp.conjugate(f.subs(x,-x))), 'Same longitudinal reflection')
                need(zero(sp.diff(applied,scale)), 'No exterior response to uniform normal scale')
                pure_bridge_checks+=1
    # First response of an actual, deck-even conformal metric variation.
    # D_eps=e^(-eps*w^2/2) D_0 e^(-eps*w^2/2), in the original L2 density.
    kk=sp.symbols('k0', real=True)
    eps=sp.symbols('epsilon', real=True)
    for b in [-1,1]:
        phi0=sp.Matrix([1,b])/2
        orth=sp.Matrix([1,-b])/2
        psi=phi0*sp.exp(sp.I*kk*x)
        u=sp.exp(-eps*w*w/2)
        direct=u*(-sp.I)*(sigma1*sp.diff(u*psi,x)+sigma2*sp.diff(u*psi,w))
        predicted=b*sp.exp(-eps*w*w)*(kk*phi0+eps*w*orth)*sp.exp(sp.I*kk*x)
        need(zero(direct-predicted),'Actual conformal Dirac first response')
        variation=sp.diff(direct,eps).subs(eps,0)*sp.exp(-sp.I*kk*x)
        compression=sp.integrate((phi0.T*variation)[0],(w,-1,1))
        need(zero(compression+b*kk/3),'Geometric input compression')
        exterior=sp.simplify(variation-phi0*compression)
        exterior_norm=sp.integrate((sp.conjugate(exterior).T*exterior)[0],(w,-1,1))
        need(zero(exterior_norm-(4*kk*kk/45+sp.Rational(1,3))), 'All-mode geometric return form')
    spectral_sum=4/sp.pi**4*(2*kk**2*sp.zeta(4)+sp.pi**2*sp.zeta(2)/2)
    need(zero(spectral_sum-(4*kk**2/45+sp.Rational(1,3))),'Independent full Fourier return sum')
    # The coupling law holds for any smooth real even conformal profile.
    # Check distinct profiles directly; the accompanying proof uses no ansatz.
    phi0=sp.Matrix([1,1])/2;orth=sp.Matrix([1,-1])/2
    for fprofile in [sp.Integer(1),w*w,w**4,sp.cos(sp.pi*w)]:
        average=sp.integrate(fprofile,(w,-1,1))/2
        exterior=-kk*(fprofile-average)*phi0+sp.diff(fprofile,w)*orth/2
        actual=sp.integrate((sp.conjugate(exterior).T*exterior)[0],(w,-1,1))
        variance=sp.integrate((fprofile-average)**2,(w,-1,1))/2
        gradient=sp.integrate(sp.diff(fprofile,w)**2,(w,-1,1))/2
        need(zero(actual-kk*kk*variance-gradient/4),'General conformal profile return law')
    # Independent Gibbs computation of the neutral negative-frequency witness.
    d=[annihilator(i,2) for i in range(2)]
    n=[x.T*x for x in d]
    c1=sp.symbols('c1', positive=True)
    rho=sp.diag((1-c1)/2,(1-c1)/2,c1/2,c1/2)
    neutral=d[0].T*d[1]
    ham=energy*n[1]
    need(zero(ham*neutral-neutral*ham+energy*neutral),'Neutral energy frequency')
    need(zero((rho*neutral.T*neutral).trace()-c1/2),'Neutral witness from Gibbs density')
    # Actual spectral values: finite numerical validation of the analytic
    # infinite-mode formulas. No claim of an infinite cutoff proof here.
    finite_modes=[]
    for n in range(-8,9):
        lam=float(sp.pi*n/2)
        rr=np.exp(-lam)
        p=np.array([[1,rr],[rr,rr*rr]])/(1+rr*rr)
        need(np.linalg.norm(p@p-p)<1e-14,'Numeric original projection')
        c=rr*rr/(1+rr*rr)
        eps=1e-5
        cp=1/(1+np.exp(2*np.exp(-eps)*lam))
        cm=1/(1+np.exp(2*np.exp(eps)*lam))
        need(abs((cp-cm)/(2*eps)-2*lam*c*(1-c))<1e-9,'Original geometric response')
        finite_modes.append(dict(n=n,energy=lam,covariance=c,cross_response=float(p[0,1]),
                                 geometric_derivative=float(2*lam*c*(1-c))))

    # A neutral energy-lowering operator already detects that a thermal state
    # is not a positive-energy ground state for the SAME H time.
    # B=a_0^* a_m: <B^*B>=C_m(1-C_0)>0, frequency=-E_m.
    witnesses=[]
    for m in [1,2]:
        lam=sp.pi*m/2
        norm=1/(2*(1+sp.exp(2*lam)))
        witnesses.append(dict(m=m,energy_frequency=str(-lam),norm_squared=str(norm),
                              norm_squared_numeric=float(norm),
                              spatial_J_invariant=(m%2==0)))
    result=dict(
        analytic_identities='Graph, purity, restriction, NS complement, full return and first geometric response',
        exact_four_point_checks=4**4,
        actual_transverse_domain_checks=14,
        pure_D_zero_branch={'exact_mode_checks':pure_bridge_checks,
                           'circle_twist':'q=c*b','circle_generator':'b*(-i partial_x)',
                           'same_state_time_reflection':True,
                           'uniform_normal_scale_return':'B=Q D_a P0=0 for every a>0',
                           'normal_locality_witness':'||Q M_w phi0||^2=1/3',
                           'scope':'Reducing zero sector of pure D sea; not identified with thermal C or P1 admissibility'},
        conformal_geometric_response={'metric':'g_eps=exp(2*eps*w^2)*g0; diagnostic profile explicitly chosen',
            'actual_D':'exp(-eps*w^2/2) D0 exp(-eps*w^2/2)',
            'first_exterior_amplitude':'b*(-k*(w^2-1/3)*phi0+w*phi_orth)',
            'first_return_form':'||Q dD/deps phi0 e^(ikx)||^2=4*k^2/45+1/3',
            'full_Fourier_amplitudes':'-2*b*(-1)^m*(k+i*pi*m/2)/(pi^2*m^2), m!=0',
            'general_profile_return':'k^2 Var(f)+<f_prime^2>/4; zero exactly for constant f',
            'physical_geometry_selection':False},
        zero_mode={'C':'1/2','cross':'1/2','first_return':'1/4'},
        negative_controls=['Mixed C is not a projection','Raw endpoint exchange is not Cauchy polarization',
                           'Periodic replacement fails the stipulated NS formula'],
        actual_modes=finite_modes, same_time_ground_state_witnesses=witnesses,
        scope='Actual comparison source: pure D sea and free thermal Green quantization kept separate; no P1 seam or primitive history selection claimed',
        physical_gates_closed=[],complete_TFPT_solution=False)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='actual_modes'},indent=2))


if __name__=='__main__':
    main()
