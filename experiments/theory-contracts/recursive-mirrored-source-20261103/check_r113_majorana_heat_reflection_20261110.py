"""Check the mixed reflected spectral action on the unchanged v252 matrices."""
import json
from collections import Counter
from itertools import product
from pathlib import Path
import sys

from mpmath import iv
import numpy as np
from scipy.integrate import quad
from scipy.linalg import eigh, expm
from scipy.special import ive
import sympy as sp

repo = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(repo / 'verification'))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'chat-triangulation-20261108'))
from check_quarter_reference import grading, higgs_matrix, mass_matrix
from v252_full_finite_triple import IDX, NP

heavy = [IDX[('nu','R',0,j)] + offset for offset in (0,NP) for j in range(3)]
light = np.array([j for j in range(96) if j not in heavy])
assert np.count_nonzero(higgs_matrix[np.ix_(heavy,heavy)]) == 0
assert np.count_nonzero(mass_matrix[light,:]) == 0
signed_m, rotation = eigh(mass_matrix[np.ix_(heavy,heavy)])
assert min(abs(signed_m)) > 0.4
coupling = rotation.conj().T @ higgs_matrix[np.ix_(heavy,light)]
massive_weights = 2 * np.sum(abs(coupling)**2,axis=1)  # Both insertion orders.
massless_weight = np.linalg.norm(higgs_matrix[np.ix_(light,light)])**2
assert abs(sum(massive_weights)+massless_weight-np.linalg.norm(higgs_matrix)**2) < 1e-11


def exact_matrix(matrix):
    return sp.Matrix([[sp.Rational(float(z.real))+sp.I*sp.Rational(float(z.imag))
                       for z in row] for row in matrix])


majorana = exact_matrix(mass_matrix[np.ix_(heavy,heavy)])
cross = exact_matrix(higgs_matrix[np.ix_(heavy,light)])
Y = exact_matrix(higgs_matrix)
assert sp.trace(Y.H*Y).expand() < 367
B = majorana**2
C = cross*cross.H
roots = [bounds for bounds,_ in B.charpoly().as_poly().intervals(eps=sp.Rational(1,10**8))]
assert len(roots)==3
(a0,a1),(b0,b1),(c0,c1)=roots
assert b0>sp.Rational(6,5) and b1<sp.Rational(3,2)
tc,tbc,tbbc=(sp.expand(sp.trace(B**k*C)) for k in (0,1,2))
# Both insertions are counted in the weight. Numerator and denominator
# are both negative because the selected eigenvalue is the middle one.
neg_numerator_lower=(a0+c0)*tbc-a1*c1*tc-tbbc
neg_denominator_upper=(b1-a0)*(c1-b0)
weight_lower=2*neg_numerator_lower/neg_denominator_upper
assert weight_lower>40


def f_prime(t,m):
    x = m*m*t/2
    return m*m*(ive(1,x)-ive(0,x))/2


def mixed_spatial_kernel_direct(t,d,m):
    """Direct Duhamel convolution differentiated twice in the separation d."""
    def integrand(theta):
        u = (1-np.cos(theta))/2
        var = t*np.sin(theta)**2/2
        return np.exp(-d*d/(2*var)-m*m*t*u) * (m*m + 1/var - d*d/var**2)
    knots = [d/np.sqrt(t)*2.0**j for j in range(-2,14)]
    points = [x for x in knots if x < np.pi/2]
    points += [np.pi-x for x in points]
    return quad(integrand,0,np.pi,points=points,epsabs=2e-8,epsrel=2e-13,limit=260)[0]/(4*np.pi)


def kernel_zero(t):
    return -sum(w*f_prime(t,m)/2 for w,m in zip(massive_weights,signed_m))


def kernel(t,d):
    massless = -massless_weight*d*np.exp(-d*d/t)/np.sqrt(np.pi*t**3)
    return massless + sum(w*mixed_spatial_kernel_direct(t,d,m)
                          for w,m in zip(massive_weights,signed_m))


values={}
for t in (0.1,1.0,5.0):
    val=kernel_zero(t)
    at_d=kernel(t,1e-5)
    assert at_d > 0 and abs(at_d-val)<0.075
    values[str(t)]={'exact_limit_numerically_evaluated':float(val),
                    'direct_separated_d_1e-5':float(at_d),
                    'direct_separated_d_point01':float(kernel(t,.01))}


def cap_rates(momentum):
    generator = momentum*grading-1j*grading@mass_matrix
    variation = -1j*grading@higgs_matrix
    energy, frame = eigh(generator)
    neg, pos = energy<0, energy>0
    matrix = frame[:,pos].conj().T@variation@frame[:,neg]
    rates = energy[pos,None]-energy[None,neg]
    distinct_masses = sorted(set(np.round(abs(signed_m),12)))
    required = np.array([2*abs(momentum)]+[abs(momentum)+np.hypot(momentum,m) for m in distinct_masses])
    active = rates[abs(matrix)>1e-10]
    assert max(min(abs(required-r)) for r in active)<1e-11
    assert all(min(abs(active-r))<1e-11 for r in required)
    return required, float(2*np.linalg.norm(matrix)**2)


def conjugated_kernel(d, rates):
    """Return exp(d*d) times the geometric mixed kernel after annihilating the cap rates."""
    coeff = np.polynomial.polynomial.polyfromroots(np.repeat(-rates,2))
    def one_mass(m):
        def integrand(theta):
            u = (1-np.cos(theta))/2
            beta = 1/np.sin(theta)**2
            deriv = [1., -2*beta*d]
            for k in range(1,len(coeff)+1):
                deriv.append(-2*beta*(d*deriv[-1]+k*deriv[-2]))
            response = m*m*np.array(deriv[:-2])-np.array(deriv[2:])
            return np.exp(-d*d*(beta-1)-m*m*u)*(coeff@response)
        return quad(integrand,0,np.pi,epsabs=1e-6,epsrel=1e-10,limit=120)[0]/(4*np.pi)
    return massless_weight*one_mass(0.)+sum(w*one_mass(m) for w,m in zip(massive_weights,signed_m))


rates, cap_mixed_seam = cap_rates(.25)
far_d, spacing = 10., .75
modified = np.array([conjugated_kernel(far_d+j*spacing,rates) for j in range(3)])
assert max(modified)<0
minor_ratio = float(np.exp(-2*spacing**2)*modified[0]*modified[2]/modified[1]**2)
assert 0 < minor_ratio < .5


def poisson_direct(m, p, a, d):
    """Mixed kernel with reversed sign after the original even exp(-a|D|) profile."""
    def theta_integrand(theta):
        u = (1-np.cos(theta))/2
        var = (np.sin(theta)/2)**2
        q = p*p + m*m*u
        k = np.sqrt(q/var)
        c2 = a*a*var
        R = np.sqrt(d*d+c2)
        return np.exp(-k*R)*(k*k*d*d/R**2 - k*c2/R**3-m*m)/np.sqrt(q)
    return a*quad(theta_integrand,0,np.pi,epsabs=5e-11,epsrel=5e-11)[0]/(2*np.pi)


def poisson_laplace(m, p, a, multiplier):
    plus = np.hypot(p,m)+abs(p)
    minus = np.hypot(p,m)-abs(p)
    def integrand(z):
        lam = plus+z*z
        smooth = np.sqrt((lam+plus)*(lam*lam-minus*minus))
        phase = a*z*smooth/(2*lam)
        return multiplier(lam)*(lam*lam-m*m)*np.cos(phase)/smooth
    return 4*a*quad(integrand,0,np.inf,epsabs=5e-10,epsrel=5e-10,limit=200)[0]/np.pi


poisson_errors=[]
for d in (.4,1.,2.):
    for m in (0., *abs(signed_m[3:])):
        direct_value=poisson_direct(m,.25,1.,d)
        inverse_value=poisson_laplace(m,.25,1.,lambda lam:np.exp(-d*lam))
        poisson_errors.append(abs(direct_value-inverse_value))
assert max(poisson_errors) < 1e-10


def poisson_density(lam,m,p,a):
    plus=np.hypot(p,m)+abs(p)
    minus=np.hypot(p,m)-abs(p)
    if lam<=plus:
        return 0.
    delta=np.sqrt((lam*lam-plus*plus)*(lam*lam-minus*minus))
    return 2*a*(lam*lam-m*m)*np.cos(a*delta/(2*lam))/(np.pi*delta)


lambda0=2*np.pi
poisson_density_at_lambda0=(massless_weight*poisson_density(lambda0,0,.25,1.)
        +sum(w*poisson_density(lambda0,m,.25,1.) for w,m in zip(massive_weights,signed_m)))
assert poisson_density_at_lambda0 < -180

order=32
distance=2*order/lambda0
def probe_multiplier(lam):
    gamma_density_scaled = np.exp(2*order*np.log(lam/lambda0)-distance*(lam-lambda0))
    return gamma_density_scaled*np.prod((lam-rates)/(lambda0-rates))**2
probe = massless_weight*poisson_laplace(0,.25,1.,probe_multiplier)
probe += sum(w*poisson_laplace(m,.25,1.,probe_multiplier) for w,m in zip(massive_weights,signed_m))
assert probe < 0


def rational_interval(x):
    return iv.mpf(str(x.p))/iv.mpf(str(x.q))


iv.dps=36
iv_roots=[iv.mpf([rational_interval(lo).a,rational_interval(hi).b]) for lo,hi in roots]
trace_C, trace_BC, trace_BBC=(rational_interval(x) for x in (tc,tbc,tbbc))
total_W=rational_interval(sp.trace(Y.H*Y).expand())
groups=[(iv.mpf(0),total_W-2*trace_C)]
for i, mass2 in enumerate(iv_roots):
    other_1,other_2=[x for j,x in enumerate(iv_roots) if j!=i]
    weight=2*(trace_BBC-(other_1+other_2)*trace_BC+other_1*other_2*trace_C)/((mass2-other_1)*(mass2-other_2))
    assert weight.a>0
    groups.append((mass2,weight))

lambda_all=8
max_p2_numerator=(2*lambda_all)**2  # p^2 in sixteenths, |p| < lambda/2.
counts=Counter((4*n+1)**2+16*(j*j+k*k) for n,j,k in product(range(-4,5),repeat=3))
all_modes_density=iv.mpf(0)
active_momentum_count=0
for p2_numerator,count in counts.items():
    if p2_numerator>=max_p2_numerator:
        continue
    p_squared=iv.mpf(p2_numerator)/16
    p=iv.sqrt(p_squared)
    active_momentum_count+=count
    for mass2,weight in groups:
        shifted=lambda_all**2-mass2
        discriminant_side=shifted-2*lambda_all*p
        if discriminant_side.b<0:
            continue
        assert discriminant_side.a>0  # No actual threshold at the certified point.
        delta=iv.sqrt(discriminant_side*(shifted+2*lambda_all*p))
        all_modes_density += 2*count*weight*shifted*iv.cos(delta/(2*lambda_all))/(iv.pi*delta)
assert all_modes_density.b < -24000


def whole_hessian_mode(t, wave, m2):
    """Before the common positive 4t and normal/tangential heat factors."""
    if wave==0:
        return -2*np.exp(-t*m2)
    value=quad(lambda u:np.exp(-t*(m2*u+wave*wave*u*(1-u))),
               0.,1.,epsabs=1e-11,epsrel=1e-10)[0]
    return t*(m2+wave*wave)*value-1-np.exp(-t*m2)


def poisson_hessian_mode(wave,m2):
    def integrand(t):
        scalar=np.exp(-1/(4*t)-t/16)/(np.pi*t)
        if scalar < 1e-23:
            return 0.
        return scalar*whole_hessian_mode(t,wave,m2)
    return quad(integrand,0.,np.inf,epsabs=2e-9,epsrel=2e-9,limit=120)[0]


q0_low=poisson_hessian_mode(0,0)
q0_high=poisson_hessian_mode(16,0)
qM_low=massless_weight*q0_low+sum(w*poisson_hessian_mode(0,m*m)
                                for w,m in zip(massive_weights,signed_m))
assert q0_low<0<q0_high and qM_low<0


def paired_aps_green_jet(distance):
    total=0.
    for p in (-.25,.25):
        H=p*grading-1j*grading@mass_matrix
        V=-1j*grading@higgs_matrix
        energy, U=eigh(H)
        gap=abs(energy)
        G_forward=(U*np.where(energy>0,np.exp(-distance*gap),0.))@U.conj().T
        G_backward=-(U*np.where(energy<0,np.exp(-distance*gap),0.))@U.conj().T
        total-=np.trace(G_forward@V@G_backward@V).real
    return total


def direct_cap_two_instant_insertions(x,y,distance):
    total=0.
    for p in (-.25,.25):
        H=p*grading-1j*grading@mass_matrix
        V=-1j*grading@higgs_matrix
        energy,U=eigh(H)
        W=U[:,energy<0]
        half0=(U*np.exp(-distance*energy/2))@U.conj().T
        # Divide each input column by its unperturbed cap growth.
        A=half0@expm(-x*V)@W@np.diag(np.exp(distance*energy[energy<0]/2))
        B=half0@expm(-y*V)@W@np.diag(np.exp(distance*energy[energy<0]/2))
        total+=np.linalg.slogdet(A.conj().T@B)[1]
    return total


green_results={}
for distance in (.3,.9):
    def centered(dx):
        return (direct_cap_two_instant_insertions(dx,dx,distance)
                -direct_cap_two_instant_insertions(dx,-dx,distance)
                -direct_cap_two_instant_insertions(-dx,dx,distance)
                +direct_cap_two_instant_insertions(-dx,-dx,distance))/(4*dx*dx)
    independent=(4*centered(1e-4)-centered(2e-4))/3
    green=paired_aps_green_jet(distance)
    assert abs(independent-green)<7e-6
    green_results[str(distance)]=[green,independent]

print(json.dumps({'zero_to_zero_weight':float(massless_weight),
                  'middle_mass_squared_verified_interval':[float(b0),float(b1)],
                  'middle_ordered_massive_weight_certified_lower':float(weight_lower),
                  'signed_mass_eigenvalues':signed_m.tolist(),
                  'ordered_massive_weights':massive_weights.tolist(),
                  'positive_mixed_heat_kernel_excluding_common_4t_tangential_factor':values,
                  'v252_cap_rates_at_p_quarter':rates.tolist(),
                  'v252_paired_cap_mixed_log_contrast_at_seam_per_unit_squared_area':cap_mixed_seam,
                  'remaining_spectral_hessian_minor_ratio_after_cap_annihilation':minor_ratio,
                  'v259_even_poisson_inverse_laplace_checks_max_error':max(poisson_errors),
                  'v259_even_poisson_signed_density_at_2pi':float(poisson_density_at_lambda0),
                  'v259_even_poisson_normalized_cap_null_probe':float(probe),
                  'poisson_cap_null_probe_extra_derivative_order':order,
                  'poisson_cap_null_probe_center_pair_distance':float(distance),
                  'v259_even_poisson_all_torus_modes_density_at_8_interval':
                        [float(all_modes_density.a),float(all_modes_density.b)],
                  'v259_even_poisson_all_torus_modes_momentum_count_at_8':active_momentum_count,
                  'moving_R032_reference_hessian_checks_at_p_quarter':
                       {'massless_per_weight_at_zero_normal_frequency':float(q0_low),
                        'massless_per_weight_at_normal_frequency_16':float(q0_high),
                        'all_96_at_zero_normal_frequency':float(qM_low)},
                  'chiral_aps_green_vs_independent_capped_transfer_jets':green_results},indent=2))
