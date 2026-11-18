"""Probe three proposed uses of the already given quarter and R063 charts."""
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import expit

source = Path(sys.argv[1])
spec = importlib.util.spec_from_file_location("r032", source)
r032 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r032)
beta = float(r032.time)


def quarter_time(t, radius=1.0):
    # The truncation is on t*p, independently for every time.
    modes = int(80 * radius/t) + 1
    p = (2*np.arange(modes)+1)/(4*radius)
    f = expit(-t*p)
    value = 0.0
    derivative = 0.0
    for mass in r032.masses:
        mass = float(mass)
        delta = mass**2/(np.hypot(p, mass)+p)
        d = np.expm1(-t*delta)
        normalization = 1+f*d
        value += 2*np.log1p(f*d).sum()
        derivative += 2*np.sum(-p*f*(1-f)*d/normalization-delta*f*(1+d)/normalization)
    return value, derivative


def short_time_anchor():
    r = 1.0
    m2 = sum(float(m)**2 for m in r032.masses)
    rows = []
    for t in (0.01, 0.002, 0.0005, 0.0001):
        value, derivative = quarter_time(t, r)
        bound = 4*r*t*m2*np.arctanh(np.exp(-t/(4*r)))
        assert 0 < -value < bound
        rows.append([t, value, derivative, bound])
    small = rows[-1][0]
    leading = (quarter_time(small)[0]/small-quarter_time(2*small)[0]/(2*small))/np.log(2)
    assert abs(leading+m2)/m2 < 0.0002
    lower = 0.01
    start = quarter_time(lower)[0]
    integral, err = quad(lambda x: np.exp(x)*quarter_time(np.exp(x))[1],
                         np.log(lower), np.log(beta), epsabs=2e-11, epsrel=2e-11)
    target = quarter_time(beta)[0]
    assert abs(start+integral-target) < 5e-11
    previous = float(r032.log_ratio(r032.mp.mpf(1),32))
    assert abs(previous-target) < 2e-13
    return dict(rows_t_log_ratio_score_bound=rows, predicted_short_t_log_coefficient=-m2,
                observed_short_t_log_coefficient=leading, direct=target,
                energy_integral_from_0p01=integral, reconstructed_from_0p01=start+integral,
                quadrature_error=err, infinite_tower_anchor_at_0=1)


def right_family(name, yukawa):
    # The right-family singular direction and all original left/right slots.
    _, singular, vstar = np.linalg.svd(yukawa)
    vector = np.zeros(96, dtype=complex)
    for family, coeff in enumerate(vstar.conj().T[:,0]):
        vector[r032.native.IDX[(name,"R",0,family)]] = coeff
    return singular[0]**2, vector


def direct_occupation(momentum, amplitude, vector):
    h = momentum*r032.grading - 1j*r032.grading @ (r032.mass_matrix+amplitude*r032.higgs_matrix)
    energies, basis = np.linalg.eigh(h)
    return np.dot(np.abs(basis.conj().T@vector)**2, expit(-beta*energies)).real


def occupation_response():
    u = 0.25
    h = np.sqrt(u)
    step = 2e-5
    rows = []
    for name,yukawa in (("u",r032.yukawa[0]),("e",r032.yukawa[2])):
        y2, vector = right_family(name,yukawa)
        for p in (0.25,0.75,1.25,2.25,4.25):
            e = np.sqrt(p*p+y2*u)
            x = beta*e/2
            tanh = np.tanh(x)
            base = 0.5+p*tanh/(2*e)
            slope = p*y2*u*(x*(1-tanh*tanh)-tanh)/(2*e**3)
            direct = direct_occupation(p,h,vector)
            diff = (direct_occupation(p,(1+step)*h,vector)
                    -direct_occupation(p,(1-step)*h,vector))/(2*step)
            scaled_u = (1+step)**2*u
            scaled_e = np.sqrt(p*p+y2*scaled_u)
            pushed = 0.5+p*np.tanh(beta*scaled_e/2)/(2*scaled_e)
            assert abs(direct-base) < 1e-13
            assert abs(diff-slope) < 2e-9
            assert abs(direct_occupation(p,(1+step)*h,vector)-pushed) < 1e-13
            rows.append([name,p,direct,diff,slope])
    return dict(u=u, rows_channel_p_occupation_direct_derivative_analytic_derivative=rows,
                signed_marginal_change_in_fixed_background=0)


def threshold_contact():
    # The tested step is deliberately on a single exact original charged-lepton
    # pair. Its application to the whole physical selector is a separate claim.
    y2, _ = right_family("e",r032.yukawa[2])
    delta = 0.25
    g = lambda p: np.logaddexp(0,-beta*p)
    dg = lambda p: -beta*expit(-beta*p)
    bulk = quad(lambda p: 2*y2*dg(p)/p,delta,20,epsabs=1e-12)[0]
    edge = 2*y2*g(delta)/delta
    moving = bulk+edge
    moved_variable = quad(lambda p: 2*y2*g(p)/p**2,delta,20,epsabs=1e-12)[0]
    assert abs(moving-moved_variable) < 1e-11 and moving > 0 > bulk
    continuum = []
    for u in (1e-5,5e-6,2.5e-6):
        lo = np.sqrt(delta*delta-y2*u)
        score = (4/u)*(
            quad(lambda p: g(np.hypot(p,np.sqrt(y2*u)))-g(p),delta,20,epsabs=1e-14)[0]
            +quad(lambda p: g(np.hypot(p,np.sqrt(y2*u))),lo,delta,epsabs=1e-14)[0])
        continuum.append([u,score])
    extrapolated = 2*continuum[-1][1]-continuum[-2][1]
    assert abs(extrapolated-moving) < 1e-4
    lattice = []
    for r in (32,128,512,2048,8192,32768,131072):
        p = (2*np.arange(8*r)+1)/(4*r)
        kept = p>=delta
        frozen = np.sum(y2*dg(p[kept])/p[kept])/r
        crossing_count = math.isqrt(r)
        lower = delta-crossing_count/(2*r)
        u = (delta**2-lower**2)/y2
        crossed = (p>=lower) & ~kept
        assert int(crossed.sum()) == crossing_count
        before = p[kept]
        new_energy = np.hypot(before,np.sqrt(y2*u))
        diff = (2/(u*r)) * (
            np.sum(g(new_energy)-g(before))
            + np.sum(g(np.hypot(p[crossed],np.sqrt(y2*u))))
        )
        lattice.append([r,frozen,diff,u,crossing_count])
    assert abs(lattice[-1][1]-bulk) < 2e-5
    assert abs(lattice[-1][2]-moving) < 0.012
    critical = delta**2/y2
    cusp = []
    for ratio in (1e-2,1e-4,1e-6,1e-8,1e-10,1e-12):
        u = critical*(1-ratio)
        lower = delta*np.sqrt(ratio)
        bulk_moving = quad(lambda p: 2*y2*dg(np.hypot(p,np.sqrt(y2*u)))/np.hypot(p,np.sqrt(y2*u)),
                           lower,20,epsabs=1e-11)[0]
        cusp.append([ratio, lower*bulk_moving + 2*y2*g(delta)])
    expected_cusp = 2*y2*g(delta)
    assert abs(cusp[-1][1]/expected_cusp - 1) < 0.0002
    return dict(diagnostic_threshold=delta,y2=y2,frozen_bulk=bulk,boundary=edge,
                moving_continuum=moving,continuum_differences=continuum,
                lattice_radius_frozen_hard_difference_u_crossings=lattice,
                field_threshold=critical, predicted_cusp=expected_cusp,
                cusp_rows=cusp)


if __name__ == "__main__":
    print(json.dumps(dict(absolute_reference=short_time_anchor(),
                          occupation_split=occupation_response(),
                          gap_contact=threshold_contact()),indent=2))
