"""Original R032/R063 transport tests using the unchanged v252 inputs."""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import expit

path = Path(sys.argv[1])
spec = importlib.util.spec_from_file_location("r032", path)
r032 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r032)
beta = float(r032.time)
M = r032.mass_matrix
Y = r032.higgs_matrix
gamma = r032.grading
m2 = np.array([float(m)**2 for m in r032.masses])


def scaled_majorana_energy():
    def scores(coupling, modes):
        p = (2*np.arange(modes)+1)[:,None]/4
        e = np.sqrt(p*p+coupling**2*m2)
        de = coupling*m2/e
        raw = beta*np.sum(de*np.tanh(beta*e/2))
        sea = -beta*np.sum(de)
        whole = -2*beta*np.sum(de*expit(-beta*e))
        assert abs(raw+sea-whole) < 2e-10
        return float(raw), float(sea), float(whole)

    coupling = 0.5
    rows = [[n,*scores(coupling,n)] for n in (32,128,512,2048,8192,32768)]
    uv = 2*beta*coupling*m2.sum()
    slope = (rows[-1][1]-rows[-2][1])/np.log(4)
    assert abs(slope/uv-1) < 1e-6
    integral, error = quad(lambda lam: scores(lam,96)[2],0,1,epsabs=1e-12,epsrel=1e-12)
    target = float(r032.log_ratio(r032.mp.mpf(1),32))
    assert abs(integral-target) < 5e-13
    p = 0.25
    mark = -1j*gamma@M
    eigenvalues, basis = np.linalg.eigh(p*gamma+coupling*mark)
    diag = np.diag(basis.conj().T@mark@basis).real
    raw = -beta*np.dot(expit(-beta*eigenvalues),diag)
    sea = beta*diag[eigenvalues<0].sum()
    raw_one, sea_one, whole_one = scores(coupling,1)
    assert max(abs(raw-raw_one),abs(sea-sea_one),abs(raw+sea-whole_one)) < 1e-12
    return dict(coupling=coupling, rows_cutoff_raw_sea_combined=rows, predicted_uv_log_slope=uv,
                measured_uv_log_slope=slope, integrated_combined_score=integral,
                previous_direct_log_ratio=target, quad_error=error,
                native_first_mode=[float(raw),float(sea),float(raw+sea)])


def opposite_twist_phase():
    n = np.arange(-128,128)
    left = n+0.25
    right = (-n-1)+0.75
    assert np.array_equal(left,-right)
    # eta(0)=zeta_H(0,a)-zeta_H(0,1-a), with zeta_H(0,x)=1/2-x.
    eta = [1-2*a for a in (0.25,0.75)]
    assert eta == [0.5,-0.5]
    relative = np.exp(-1j*np.pi*(eta[0]-eta[1])/2)
    assert abs(relative+1j)<1e-15
    return dict(modulus_spectra_equal=True, individual_eta=eta,
                oriented_weight_phase=[float(relative.real),float(relative.imag)])


def coordinate_empty_cap():
    radius = 1.0
    expected = radius*m2.sum()
    rows = []
    for modes in (32,128,512,2048,8192,32768,131072):
        p = (2*np.arange(modes)+1)[:,None]/(4*radius)
        energy_shift = -np.sum(m2/(np.sqrt(p*p+m2)+p))
        t = 1/np.log(modes)
        rows.append([modes,float(energy_shift),float(beta*energy_shift),
                     float(np.exp(t*energy_shift)),float(t)])
    slope = (rows[-1][1]-rows[-2][1])/np.log(4)
    assert abs(slope+expected)<1e-6
    assert abs(np.log(rows[-1][3])+expected)<abs(np.log(rows[0][3])+expected)
    # Direct original matrix and empty exterior coefficient Γ(exp(-t h))_00=1.
    p = 0.25
    ev = np.linalg.eigvalsh(p*gamma-1j*gamma@M)
    native_one = ev[ev<0].sum()+48*p
    expected_one = -np.sum(np.sqrt(p*p+m2)-p)
    assert abs(native_one-expected_one)<1e-13
    return dict(rows_modes_sea_delta_log_cap_at_beta_cap_at_t_times_logN_1_t=rows,
                predicted_logN_slope=-float(expected), measured_logN_slope=float(slope),
                native_first_mode_delta=float(native_one),
                joint_limit_at_t_logN_1=float(np.exp(-expected)))


def reference_slater_cap():
    p = 0.25
    native_time = 0.3
    energies, basis = np.linalg.eigh(p*gamma-1j*gamma@M)
    _, reference = np.linalg.eigh(p*gamma)
    overlap = reference[:,:48].conj().T@basis
    compressed = (overlap*np.exp(-native_time*energies))@overlap.conj().T
    sign, logdet = np.linalg.slogdet(compressed)
    assert abs(sign-1)<1e-12
    native = logdet+native_time*energies[:48].sum()

    def factors(momentum, time):
        energy = np.sqrt(momentum*momentum+m2)
        missing_overlap = m2/(2*energy*(energy+momentum))
        damping = 1 if np.isinf(time) else -np.expm1(-2*time*energy)
        return np.log1p(-missing_overlap*damping)

    one = factors(p,native_time).sum()
    assert abs(native-one)<1e-12
    log_semigroup_defect = factors(p,2*native_time).sum()-2*one
    assert log_semigroup_defect>0.1
    rows = []
    times = (0,1e-4,0.1,beta,np.inf)
    for modes in (128,2048,32768,131072):
        p = (2*np.arange(modes)+1)[:,None]/4
        # Uniform in t and independent of the choice of truncated time window.
        tail_bound = 2*m2.sum()/(modes-0.5)
        rows.append([modes,float(tail_bound),*[float(factors(p,t).sum()) for t in times]])
    for shorter,longer in zip(rows,rows[1:]):
        assert np.all(np.abs(np.array(shorter[2:])-longer[2:])<=shorter[1])
    assert rows[-1][2]==0 and rows[-1][3]>-0.1
    assert np.all(np.diff(rows[-1][2:])<=0)
    return dict(time_grid=[0,1e-4,0.1,beta,"infinity"],
                rows_modes_tail_bound_and_log_caps=rows,
                first_mode_log_semigroup_defect=float(log_semigroup_defect),
                native_first_mode_log_cap=float(native),analytic_first_mode_log_cap=float(one))


def cap_infrared_response():
    diagonal, offdiagonal = r032.grouped_higgs_coefficients()
    diagonal = [(float(lam),float(c)) for lam,c in diagonal]
    offdiagonal = [(float(lam),float(mu),float(c)) for lam,mu,c in offdiagonal if c>1e-18]
    kernel_weight = diagonal[0][1]-sum(c/mu for lam,mu,c in offdiagonal if lam==0 and mu>0)
    slope = -beta*kernel_weight/2

    def derivatives(p, x, cap):
        e = np.sqrt(x)
        if not cap:
            f = expit(-beta*e)
            return -beta*f/(2*e), beta*f/(4*e**3)+beta**2*f*(1-f)/(4*e**2)
        a = (x-p*p)/(2*e*(e+p))
        da = p/(4*e**3)
        dda = -3*p/(8*e**5)
        g = -np.expm1(-2*beta*e)
        dg = beta*np.exp(-2*beta*e)/e
        ddg = -dg*(beta/e+0.5/e**2)
        w = 1-a*g
        dw = -da*g-a*dg
        ddw = -dda*g-2*da*dg-a*ddg
        return dw/(2*w), (ddw/w-(dw/w)**2)/2

    def response(p, cap):
        ds = {lam: derivatives(p,p*p+lam,cap) for lam,_ in diagonal}
        result = sum(c*ds[lam][0] for lam,c in diagonal)
        for lam,mu,c in offdiagonal:
            term = ds[lam][1] if lam==mu else (ds[lam][0]-ds[mu][0])/(lam-mu)
            result += c*term/2
        return result

    p = 0.25
    vR = np.linalg.eigh(p*gamma)[1][:,:48]

    def cap_score(x):
        eigenvalues, basis = np.linalg.eigh(p*gamma-1j*gamma@(M+x*Y))
        overlap = vR.conj().T@basis
        # Evaluate at the original time in the reference-paired invariant planes,
        # since exponentiation and determinant of the dense unscaled matrix are ill conditioned.
        e2 = np.linalg.eigvalsh(p*p*np.eye(96)+(M+x*Y)@(M+x*Y))
        e = np.sqrt(e2)
        assert np.max(abs(e-np.sort(abs(eigenvalues))))<1e-12
        assert np.linalg.norm(overlap@overlap.conj().T-np.eye(48))<1e-12
        a = (e2-p*p)/(2*e*(e+p))
        return np.log1p(a*np.expm1(-2*beta*e)).sum()/2

    baseline = cap_score(0)
    def symmetric_difference(h):
        return (cap_score(h)+cap_score(-h)-2*baseline)/(2*h*h)

    finite_difference = (4*symmetric_difference(5e-5)-symmetric_difference(1e-4))/3
    assert abs(finite_difference-response(p,True))<1e-5
    integrand = lambda p: 2*(response(p,False)-response(p,True))
    integral, error = quad(integrand,0,np.inf,epsabs=2e-9,epsrel=2e-11,limit=250)
    tail, tail_error = quad(integrand,64,np.inf,epsabs=2e-9,epsrel=2e-11,limit=250)
    rows = []
    for r in (32,64,128,256,512,1024):
        p = (2*np.arange(128*r)+1)/(4*r)
        trace = response(p,False).sum()/r
        cap = response(p,True).sum()/r
        rows.append([r,float(trace),float(cap),float(trace-cap+tail)])
    trace_slope, cap_slope = (np.array(rows[-1][1:3])-rows[-2][1:3])/np.log(2)
    assert max(abs(trace_slope-slope),abs(cap_slope-slope))<0.01
    extrapolated = (4*rows[-1][-1]-rows[-2][-1])/3
    assert abs(extrapolated-integral)<1e-6
    return dict(predicted_log_radius_slope=float(slope),measured_trace_slope=float(trace_slope),
                measured_slater_slope=float(cap_slope),
                finite_trace_minus_cap_density=float(integral),quadrature_error=error,
                extrapolated_trace_minus_cap_density=float(extrapolated),
                tail_correction_at_p_64=float(tail),tail_error=tail_error,
                first_mode_response=float(response(0.25,True)),
                first_mode_finite_difference=float(finite_difference),rows=rows)


def ordered_occupation_word():
    h0 = 0.5
    p = 0.25
    _, _, vs = np.linalg.svd(r032.yukawa[0])
    vector = np.zeros(96,dtype=complex)
    for i,value in enumerate(vs.conj().T[:,0]):
        vector[r032.native.IDX[("u","R",0,i)]]=value
    n = np.outer(vector,vector.conj())
    y = -1j*gamma@(h0*Y)
    h_mass = p*gamma-1j*gamma@M
    energies, basis = np.linalg.eigh(h_mass+y)
    f = expit(-beta*energies)
    a = basis.conj().T@n@basis
    b = basis.conj().T@y@basis
    products = a*b.T
    diff = energies[:,None]-energies[None,:]
    weights = f[:,None]*expit(beta*energies)[None,:]

    def word(tau):
        value = np.sum(weights*np.exp(tau*diff)*products)
        assert abs(value.imag)<5e-12
        return float(value.real)

    close = abs(diff)<1e-10
    integral_kernel = np.empty_like(diff)
    integral_kernel[close] = beta*weights[close]
    difference = f[None,:]-f[:,None]
    integral_kernel[~close] = difference[~close]/diff[~close]
    analytic_integral = np.sum(integral_kernel*products).real
    integrated, error = quad(word,0,beta,epsabs=1e-13,epsrel=1e-13)
    assert abs(integrated-analytic_integral)<5e-13

    def expectation(k):
        e, v = np.linalg.eigh(h_mass+k*y)
        return np.dot(np.abs(v.conj().T@vector)**2,expit(-beta*e)).real

    eps = 1e-5
    derivative = (expectation(1+eps)-expectation(1-eps))/(2*eps)
    assert abs(derivative+analytic_integral)<2e-10
    rows = [[float(tau/beta),word(tau)] for tau in np.linspace(0,beta,5)]
    assert abs(rows[0][1]-rows[2][1])>0.05
    return dict(p=p,higgs_amplitude=h0, rows_relative_time_connected_word=rows,
                integrated_word=float(analytic_integral), direct_occupation_derivative=float(derivative),
                quadrature_error=error, external_higgs_fixed=True)


def moving_projector():
    p = 0.25
    eigM,vM = np.linalg.eigh(M)
    active = abs(eigM)>1e-10
    Q = vM[:,active]@vM[:,active].conj().T
    P0 = np.eye(96)-Q
    pinv = (vM[:,active]/eigM[active])@vM[:,active].conj().T
    theta = p*p+min(eigM[active]**2)/2
    x = P0@Y@pinv
    analytic = x+x.conj().T
    W = np.diag([-1 if field[1]=="L" else 1 for field in r032.native.P]*2)
    assert np.linalg.norm(W@M@W-M)<1e-14 and np.linalg.norm(W@Y@W+Y)<1e-14
    assert np.linalg.norm(analytic)>1
    L = M@Y+Y@M
    insertion = np.trace(analytic@L).real
    expected_insertion = 2*np.linalg.norm(P0@Y@Q)**2
    assert abs(insertion-expected_insertion)<2e-12

    def projector(a):
        phi=M+a*Y
        eigenvalues,basis=np.linalg.eigh(p*p*np.eye(96)+phi@phi)
        keep=eigenvalues>=theta
        assert keep.sum()==6
        return basis[:,keep]@basis[:,keep].conj().T

    assert np.linalg.norm(projector(0)-Q)<1e-12
    rows = []
    for eps in (0.002,0.0002,0.00002):
        plus, minus = projector(eps), projector(-eps)
        derivative = (plus-minus)/(2*eps)
        assert np.linalg.norm(W@plus@W-minus)<1e-12
        rows.append([eps,np.linalg.norm(derivative-analytic),
                     np.linalg.norm(plus@W-W@plus)/eps,
                     np.trace(derivative@L).real])
    assert rows[-1][1]<1e-6 and abs(rows[-1][2]-2*np.linalg.norm(analytic))<1e-6
    return dict(p=p,diagnostic_squared_threshold=float(theta),rank=6,
                projector_derivative_hs_norm=float(np.linalg.norm(analytic)),
                predicted_weak_center_commutator_slope=float(2*np.linalg.norm(analytic)),
                neutral_squared_generator_insertion=float(insertion),rows=rows)


if __name__ == "__main__":
    print(json.dumps(dict(energy=scaled_majorana_energy(),twist_phase=opposite_twist_phase(),
                          coordinate_empty_cap=coordinate_empty_cap(),
                          reference_slater_cap=reference_slater_cap(),
                          cap_infrared_response=cap_infrared_response(),
                          ordered_word=ordered_occupation_word(),projector=moving_projector()),indent=2))
