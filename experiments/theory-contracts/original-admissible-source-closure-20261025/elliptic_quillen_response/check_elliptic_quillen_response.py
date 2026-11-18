"""Relative flat-connection determinant on the retained square NS/NS curve.

No new geometry, heat profile, EM identification or P1 measure is selected.
The relative heat integral independently checks the theta determinant.
"""
from pathlib import Path
import json
import mpmath as mp

HERE=Path(__file__).resolve().parent


def require(condition,message):
    if not condition:raise RuntimeError(message)


def theta(z):
    # Our Jacobi argument has period one; mpmath's argument has period pi.
    return mp.jtheta(3,mp.pi*z,mp.exp(-mp.pi))


def log_norm(z):
    return -2*mp.pi*mp.im(z)**2+2*mp.log(abs(theta(z)/theta(0)))


def heat_circle(t,a):
    a=a-mp.nint(a)
    if t<mp.mpf('.03'):
        return (1+2*sum(mp.exp(-k*k/(4*t))*mp.cos(2*mp.pi*k*a)
                        for k in range(1,21)))/mp.sqrt(4*mp.pi*t)
    return sum(mp.exp(-4*mp.pi**2*t*(n+a)**2) for n in range(-20,21))


def relative_heat_integral(z):
    u,v=mp.im(z),mp.re(z)
    def difference(t):
        if t==0:return mp.mpf(0)
        return (heat_circle(t,mp.mpf('.5')+v)*heat_circle(t,mp.mpf('.5')+u)
                -heat_circle(t,mp.mpf('.5'))**2)/t
    return -mp.quad(difference,[0,mp.mpf('.03'),mp.mpf('.2'),1,mp.inf])


def product_theta_ratio(z):
    # Independent Jacobi product for theta3(z)/theta3(0).
    q=mp.exp(-mp.pi);w=mp.exp(2*mp.pi*1j*z)
    return mp.fprod((1+w*q**(2*n-1))*(1+q**(2*n-1)/w)
                    /(1+q**(2*n-1))**2 for n in range(1,81))


def main():
    mp.mp.dps=60
    probes=[mp.mpc('.17','.23'),mp.mpc('.5','0'),mp.mpc('.11','-.19')]
    rows=[]
    for z in probes:
        target=log_norm(z);heat=relative_heat_integral(z)
        require(abs(heat-target)<mp.mpf('1e-45'),'independent relative spectral heat determinant')
        require(abs(product_theta_ratio(z)-theta(z)/theta(0))<mp.mpf('1e-55'),'independent Jacobi product')
        err=max(abs(log_norm(z+1)-target),abs(log_norm(z+1j)-target),abs(log_norm(1j*z)-target))
        require(err<mp.mpf('1e-55'),'both large gauges and square action preserve Quillen norm')
        rows.append(dict(z=[float(mp.re(z)),float(mp.im(z))],relative_log_determinant=float(target),
                         heat_integral_error=float(abs(heat-target)),gauge_and_quarter_error=float(err)))
    hx=mp.diff(lambda x:log_norm(mp.mpc(x,0)),0,2)
    hy=mp.diff(lambda y:log_norm(mp.mpc(0,y)),0,2)
    require(abs(hx+2*mp.pi)<mp.mpf('1e-50') and abs(hy+2*mp.pi)<mp.mpf('1e-50'),'isotropic exact square susceptibility')
    zero=mp.mpc('.5','.5')
    require(abs(theta(zero))<mp.mpf('1e-55'),'explicit null-mode divisor')
    circle=[theta(zero+mp.mpf('.07')*mp.exp(2*mp.pi*1j*j/128)) for j in range(129)]
    divisor_winding=sum(mp.arg(circle[j+1]/circle[j]) for j in range(128))/(2*mp.pi)
    require(abs(divisor_winding-1)<mp.mpf('1e-55'),'independent determinant-section divisor winding')
    z=mp.mpc('.17','.23')
    coefficient=lambda w:mp.exp(-mp.pi*mp.im(w)**2)*theta(w)
    phase=coefficient(z+1j)/coefficient(z)
    require(abs(phase-mp.exp(-2*mp.pi*1j*mp.re(z)))<mp.mpf('1e-55'),'unitary theta coefficient transition')
    # Same old full-NS winding, finite open-window telescoping only.
    # q is a diagnostic counting parameter, not alpha or a gauge holonomy.
    q=mp.mpf(1)/3
    occ=lambda n:1/(1+mp.exp(2*mp.pi*(mp.mpf(n)+mp.mpf('.5'))))
    relative=(1-q*occ(-33))/(1-q*occ(32))
    require(abs(relative-(1-q))<mp.mpf('1e-55'),'retained thermal winding factor')
    c3=1/(8*mp.pi);q0=48*c3**4
    target_slope=-mp.mpf('2.5')*q0/(1-q0)
    require(target_slope<0,'old minimal EM transfer has nonzero slope at alpha0')
    result=dict(curve='E=C/(Z+iZ), original untwisted nontrivial NS/NS spin',
                operator='K_z=(-i partial_x+2 pi v)+i(-i partial_y+2 pi u), z=v+i u; fixed H1_NSNS domain',
                positive_spectrum='K_z* K_z: 4 pi^2[(m+1/2+v)^2+(n+1/2+u)^2]',
                relative_log_determinant='log[det_zeta(K_z* K_z)/det_zeta(K_0* K_0)]=-2 pi u^2+2 log|theta3(z|i)/theta3(0|i)|',
                heat_contract='-integral_0^infinity [Tr exp(-t L_z)-Tr exp(-t L_0)]dt/t; flat UV coefficients cancel; invertible probes only',
                probes=rows,quillen_log_norm_hessian=[[float(hx),0.0],[0.0,float(hy)]],
                vector_dirac_effective_action='Gamma_vec=-log det_zeta L_z; Hessian at0=2 pi I2; chiral modulus has half this coefficient',
                null_divisor='z=(1+i)/2 modulo Z+iZ; one complex zero mode of K, simple theta zero, quadratic norm zero',
                independently_counted_divisor_winding=float(divisor_winding),
                line_convention='positive theta/determinant-section line; explicitly dualize if comparing with occupied Berry convention',
                line_unitary_frame_connection='A=2 pi i u dv; F=-2 pi i dv wedge du; (i/2pi) integral F=+1',
                frame_transitions='e_unit(z+1)=e_unit(z), e_unit(z+i)=exp(2pi i v)e_unit(z)',
                line_holonomies='v-cycle exp(-2 pi i u); u-cycle exp(2 pi i v), including endpoint frame identification',
                old_counting_q=float(q),old_relative_winding_determinant=float(relative),
                full_spectral_large_gauge_ratio=1.0,
                smooth_zero_twist_EM_test=dict(condition='Only if a proposed smooth z(alpha) has z(0)=0 and uses this norm alone',
                                               quillen_first_derivative=0.0,old_minimal_transfer_first_derivative=float(target_slope),
                                               verdict='They cannot agree up to a constant reference or fixed multiplicity under that condition; no exclusion for a derived nonzero background or additional original insertions'),
                conclusion='Full spectral norm and the old occupation-pair counting determinant have distinct insertions despite sharing the old circle and winding.',
                not_derived=['primitive P1 charged source/measure and cross-face embedding',
                             'finite v472 Berry-to-continuum metric/connection limit',
                             'identification z=z(alpha_EM) or determinant power5/4',
                             'positive CAR/KMS state for arbitrary temporal gauge twist',
                             'physical time, interacting native source, or complete TFPT solution'],
                physical_gates_closed=[],complete_TFPT_solution=False)
    out=HERE/'runs/elliptic_quillen_response.json';out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(result,indent=2,ensure_ascii=False))


if __name__=='__main__':main()
